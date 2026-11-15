//! Streaming ONT deduplication, mirroring the Illumina `DedupContext`.
//!
//! Reads arrive one at a time. Each read is compared only against cluster
//! leaders whose seeds match it at both ends. It joins the first leader
//! it matches; otherwise it starts a new cluster and becomes that cluster's
//! leader. Leaders never change, so cluster membership is decided against a
//! fixed sequence, and clusters never merge. The best read of each cluster
//! (by quality, then length) is tracked and resolved in `finalize()`.

use super::compare::{
    OntDedupParams, end_windows, ends_match, lengths_compatible, whole_read_matches,
};
use super::leaders::LeaderStore;
use super::reverse_complement_into;
use super::seeds::{EndSeeds, OntMinimizerParams, end_seeds};
use ahash::AHashSet;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

/// Most seeds belong to a single leader, so a bucket holds one posting
/// inline and allocates only when it grows.
type Bucket = SmallVec<[Posting; 1]>;

/// A leader's entry in a seed's bucket: which end of the leader the seed
/// came from, and its offset from that end (see `Seed`).
#[derive(Debug, Clone, Copy)]
struct Posting {
    cluster: u32,
    offset: u16,
    at_end: bool,
}

// Evidence bits for a candidate leader: which end of the read matched a seed
// from which end of the leader. Same orientation needs START_START and
// END_END; reverse complement needs START_END and END_START.
const START_START: u8 = 1;
const END_END: u8 = 2;
const START_END: u8 = 4;
const END_START: u8 = 8;
const SAME_STRAND: u8 = START_START | END_END;
const OPPOSITE_STRAND: u8 = START_END | END_START;
// Set once a leader has been aligned in that orientation.
const TRIED_SAME: u8 = 16;
const TRIED_OPPOSITE: u8 = 32;

#[derive(Debug, Clone)]
struct Cluster {
    best_read_idx: u32,
    best_score: Score,
}

/// Whether a seed at `a` bases from its read end in one read and `b` bases
/// in another could belong to duplicates. Between duplicates a seed's offset
/// differs by the endpoint displacement plus the net indels between the seed
/// and the read end. Within `end_check_len` of the end, and assuming the
/// shared seed is the same (homologous) k-mer in both reads, this bound is
/// exact: the endpoint check allows at most `max_end_offset` +
/// `max_end_errors`. In repetitive sequence a shared seed may be a different
/// copy of the k-mer, and the filter can then skip a duplicate.
/// Further in, indels beyond the end window are only bounded by the
/// whole-read error rate, so the allowance grows with distance; this part is
/// a heuristic, not a guarantee, and a rare dense indel cluster could hide a
/// duplicate.
fn offsets_compatible(a: u16, b: u16, params: &OntDedupParams) -> bool {
    let drift = (params.max_error_frac * f64::from(a.max(b))).ceil() as usize;
    usize::from(a.abs_diff(b)) <= params.max_end_offset + params.max_end_errors + drift
}

/// Per-leader evidence bits for the current read: a dense array indexed by
/// cluster id. In dense groups a read can touch thousands of postings, and
/// this is measurably faster than a hash map. Entries from earlier reads are
/// recognized by their stamp, so nothing needs clearing between reads.
#[derive(Default)]
struct Evidence {
    stamp: u32,
    entries: Vec<(u32, u8)>,
}

impl Evidence {
    fn start_read(&mut self, clusters: usize) {
        self.entries.resize(clusters, (0, 0));
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            // Wrapped: stale entries could now look current, so reset them.
            self.entries.fill((0, 0));
            self.stamp = 1;
        }
    }

    /// Evidence bits for `cluster` (zero if untouched for this read).
    fn bits(&mut self, cluster: u32) -> &mut u8 {
        let entry = &mut self.entries[cluster as usize];
        if entry.0 != self.stamp {
            *entry = (self.stamp, 0);
        }
        &mut entry.1
    }
}

/// Buffers for unpacking a leader during comparison.
#[derive(Default)]
struct LeaderScratch {
    start: Vec<u8>,
    end: Vec<u8>,
    full: Vec<u8>,
}

/// `is_duplicate(read, leader)` for a packed leader, unpacking only what
/// each stage needs: the end windows first, the whole leader only if the
/// ends match.
fn matches_leader(
    read: &[u8],
    leaders: &LeaderStore,
    id: u32,
    scratch: &mut LeaderScratch,
    params: &OntDedupParams,
) -> bool {
    let len = leaders.len(id);
    if !lengths_compatible(read.len(), len, params) {
        return false;
    }
    let (start, end) = end_windows(len, params);
    leaders.unpack_into(id, start, &mut scratch.start);
    leaders.unpack_into(id, end, &mut scratch.end);
    if !ends_match(read, &scratch.start, &scratch.end, params) {
        return false;
    }
    leaders.unpack_into(id, 0..len, &mut scratch.full);
    whole_read_matches(read, &scratch.full, params)
}

/// The offsets `b` compatible with `a` under `offsets_compatible`, which
/// form one contiguous range: below `a` the allowance is fixed, and above
/// `a` it grows by at most one per base when `max_error_frac` < 1.
fn compatible_offset_range(a: u16, params: &OntDedupParams) -> (u16, u16) {
    let base = params.max_end_offset + params.max_end_errors;
    let drift_at_a = (params.max_error_frac * f64::from(a)).ceil() as usize;
    let lo = usize::from(a).saturating_sub(base + drift_at_a) as u16;
    let mut hi = if params.max_error_frac < 1.0 {
        // Just above the largest b with b - a <= base + frac * b + 1.
        let estimate = (f64::from(a) + base as f64 + 1.0) / (1.0 - params.max_error_frac);
        estimate.min(f64::from(u16::MAX)) as u16
    } else {
        u16::MAX
    };
    while hi > a && !offsets_compatible(a, hi, params) {
        hi -= 1;
    }
    (lo, hi)
}

/// The postings in an offset-sorted bucket whose offsets are compatible with
/// `offset`.
fn compatible_postings<'a>(
    bucket: &'a [Posting],
    offset: u16,
    params: &OntDedupParams,
) -> &'a [Posting] {
    let (lo, hi) = compatible_offset_range(offset, params);
    let first = bucket.partition_point(|p| p.offset < lo);
    let last = bucket.partition_point(|p| p.offset <= hi);
    &bucket[first..last]
}

/// Exemplar ranking: mean Phred quality first, then length, compared
/// lexicographically. The arithmetic mean of Phred scores deliberately
/// matches the Illumina path (issue #6 discusses better quality summaries).
type Score = (f64, usize);

fn score(seq: &[u8], qual: &[u8]) -> Score {
    let mean_q = if qual.is_empty() {
        0.0
    } else {
        qual.iter()
            .map(|&q| u32::from(q.saturating_sub(33)))
            .sum::<u32>() as f64
            / qual.len() as f64
    };
    (mean_q, seq.len())
}

/// Counters describing the work done, for tuning and diagnostics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OntDedupStats {
    /// Reads processed.
    pub reads: usize,
    /// Clusters created (one per leader).
    pub clusters: usize,
    /// Reads without seeds at both ends (too short or no valid k-mers);
    /// each becomes a singleton cluster.
    pub unseeded_reads: usize,
    /// Bucket entries examined, summed over reads.
    pub postings_visited: usize,
    /// Alignments against leaders with seed evidence at both ends (one per
    /// orientation tried), summed over reads.
    pub candidates_compared: usize,
}

/// Streaming deduplicator for single-end ONT reads; see the module docs and
/// `docs/ont-design.md`. Feed reads with `process_read_by_index`, call
/// `finalize`, then query exemplars.
pub struct OntDedupContext {
    dedup_params: OntDedupParams,
    minimizer_params: OntMinimizerParams,

    /// Seed hash -> leaders with that seed at either end, sorted by offset
    /// (then insertion order) so a lookup can visit just the compatible
    /// offset range. Dense groups make some buckets very large.
    buckets: FxHashMap<u64, Bucket>,
    /// Leader sequences, indexed by cluster id.
    leaders: LeaderStore,

    /// Indexed by cluster id (the order clusters were created).
    clusters: Vec<Cluster>,
    /// Read index -> cluster id while streaming; read index -> exemplar read
    /// index after `finalize()`.
    results: Vec<u32>,

    stats: OntDedupStats,
    finalized: bool,

    // Per-read scratch space, kept to avoid reallocating.
    hash_scratch: Vec<u64>,
    read_scratch: Vec<u8>,
    rc_scratch: Vec<u8>,
    leader_scratch: LeaderScratch,
    /// Evidence bits per leader for the current read, and leaders ready to
    /// align.
    evidence: Evidence,
    candidates: Vec<u32>,
}

impl OntDedupContext {
    /// Panics if `minimizer_params` would fail `OntMinimizerParams::new`, or
    /// if `dedup_params.max_error_frac` is not between 0 and 1.
    pub fn new(dedup_params: OntDedupParams, minimizer_params: OntMinimizerParams) -> Self {
        let p = &minimizer_params;
        if let Err(e) = OntMinimizerParams::new(p.kmer_len, p.window_len, p.num_windows) {
            panic!("invalid OntMinimizerParams: {e}");
        }
        assert!(
            (0.0..=1.0).contains(&dedup_params.max_error_frac),
            "max_error_frac must be between 0 and 1, got {}",
            dedup_params.max_error_frac
        );
        Self {
            dedup_params,
            minimizer_params,
            buckets: FxHashMap::default(),
            leaders: LeaderStore::default(),
            clusters: Vec::new(),
            results: Vec::new(),
            stats: OntDedupStats::default(),
            finalized: false,
            hash_scratch: Vec::new(),
            read_scratch: Vec::new(),
            rc_scratch: Vec::new(),
            leader_scratch: LeaderScratch::default(),
            evidence: Evidence::default(),
            candidates: Vec::new(),
        }
    }

    /// Process one read. `idx` must be the read's position in the input:
    /// indices start at 0, each is used once, and they should be contiguous.
    /// (As in the Illumina path, an index that is skipped is still counted by
    /// `stats` and reports cluster 0's exemplar.) Quality is Phred+33.
    ///
    /// Returns the cluster id the read was assigned to.
    pub fn process_read_by_index(&mut self, idx: usize, seq: &[u8], qual: &[u8]) -> u32 {
        assert!(
            !self.finalized,
            "process_read_by_index called after finalize()"
        );
        let read_idx = u32::try_from(idx).expect("read index must fit in u32");
        let score = score(seq, qual);
        self.stats.reads += 1;

        // Comparisons are case-sensitive, so normalize case once.
        let mut read = std::mem::take(&mut self.read_scratch);
        read.clear();
        read.extend(seq.iter().map(u8::to_ascii_uppercase));

        let seeds = end_seeds(&read, &self.minimizer_params, &mut self.hash_scratch);
        if seeds.start.is_empty() || seeds.end.is_empty() {
            self.stats.unseeded_reads += 1;
        }
        let matched = self.find_matching_cluster(&read, &seeds);

        let cluster_id = match matched {
            Some(cluster_id) => {
                let cluster = &mut self.clusters[cluster_id as usize];
                if score > cluster.best_score {
                    cluster.best_read_idx = read_idx;
                    cluster.best_score = score;
                }
                cluster_id
            }
            None => {
                let cluster_id = u32::try_from(self.clusters.len()).expect("too many clusters");
                self.clusters.push(Cluster {
                    best_read_idx: read_idx,
                    best_score: score,
                });
                self.leaders.push(&read);
                // A leader needs seeds at both ends to ever be found.
                if !seeds.start.is_empty() && !seeds.end.is_empty() {
                    for (at_end, end) in [(false, &seeds.start), (true, &seeds.end)] {
                        for seed in end {
                            let bucket = self.buckets.entry(seed.hash).or_default();
                            let pos = bucket.partition_point(|p| p.offset <= seed.offset);
                            let posting = Posting {
                                cluster: cluster_id,
                                offset: seed.offset,
                                at_end,
                            };
                            bucket.insert(pos, posting);
                        }
                    }
                }
                cluster_id
            }
        };
        self.read_scratch = read;

        if self.results.len() <= idx {
            self.results.resize(idx + 1, 0);
        }
        self.results[idx] = cluster_id;
        cluster_id
    }

    /// First leader that `read` duplicates in an allowed orientation.
    ///
    /// Seeds are consulted from the read ends inward, one window per end at
    /// a time. Each leader posting with a compatible offset adds evidence
    /// for the (read end, leader end) pair it links. As soon as a leader has
    /// evidence at both ends in a consistent orientation it is aligned, and
    /// the search stops at the first match. Seeds near the ends are the most
    /// precise, and in dense groups most reads match early, which avoids
    /// scanning the large buckets of inner seeds.
    fn find_matching_cluster(&mut self, read: &[u8], seeds: &EndSeeds) -> Option<u32> {
        self.evidence.start_read(self.clusters.len());
        let mut rc_ready = false;
        for window in 0..seeds.start.len().max(seeds.end.len()) {
            self.collect_evidence(seeds, window);
            if let Some(cluster) = self.align_ready(read, &mut rc_ready) {
                return Some(cluster);
            }
        }
        None
    }

    /// Add evidence from both ends' seeds in `window`, queueing leaders that
    /// just gained two-end evidence in some orientation.
    fn collect_evidence(&mut self, seeds: &EndSeeds, window: usize) {
        let params = &self.dedup_params;
        let evidence = &mut self.evidence;
        let ready = &mut self.candidates;
        ready.clear();
        for (read_at_end, end) in [(false, &seeds.start), (true, &seeds.end)] {
            let Some(seed) = end.get(window) else {
                continue;
            };
            let Some(bucket) = self.buckets.get(&seed.hash) else {
                continue;
            };
            let in_range = compatible_postings(bucket, seed.offset, params);
            self.stats.postings_visited += in_range.len();
            for posting in in_range {
                let bit = match (read_at_end, posting.at_end) {
                    (false, false) => START_START,
                    (true, true) => END_END,
                    (false, true) => START_END,
                    (true, false) => END_START,
                };
                let bits = evidence.bits(posting.cluster);
                let before = *bits;
                *bits |= bit;
                let completed = |b: u8, mask: u8| b & mask != mask && *bits & mask == mask;
                // Either pattern queues the leader, even without RC tolerance:
                // misleading seeds can produce crossed evidence for a forward
                // duplicate, and `align_ready` applies the strand restriction.
                if completed(before, SAME_STRAND) || completed(before, OPPOSITE_STRAND) {
                    ready.push(posting.cluster);
                }
            }
        }
    }

    /// Align `read` against the queued leaders, in queue order; returns the
    /// first match. Seed evidence picks the orientation to try first, but
    /// the other allowed orientation is tried too: in terminal inverted
    /// repeats, seeds from one end can match the other end's sequence and
    /// point the wrong way. Each orientation is tried at most once per
    /// leader.
    fn align_ready(&mut self, read: &[u8], rc_ready: &mut bool) -> Option<u32> {
        let params = &self.dedup_params;
        for &cluster in self.candidates.iter() {
            let bits = self.evidence.bits(cluster);
            let order = if *bits & SAME_STRAND == SAME_STRAND {
                [(false, TRIED_SAME), (true, TRIED_OPPOSITE)]
            } else {
                [(true, TRIED_OPPOSITE), (false, TRIED_SAME)]
            };
            for (reverse, tried) in order {
                let bits = self.evidence.bits(cluster);
                if *bits & tried != 0 || (reverse && !params.rc_tolerant) {
                    continue;
                }
                *bits |= tried;
                self.stats.candidates_compared += 1;
                let oriented = if reverse {
                    // Computed once per read, and only if needed.
                    if !*rc_ready {
                        reverse_complement_into(read, &mut self.rc_scratch);
                        *rc_ready = true;
                    }
                    &self.rc_scratch[..]
                } else {
                    read
                };
                let scratch = &mut self.leader_scratch;
                if matches_leader(oriented, &self.leaders, cluster, scratch, params) {
                    return Some(cluster);
                }
            }
        }
        None
    }

    /// Resolve every read to its cluster's best read. Must be called before
    /// querying exemplars; frees the seed index and leader sequences.
    pub fn finalize(&mut self) {
        if self.finalized {
            return;
        }
        for result in &mut self.results {
            *result = self.clusters[*result as usize].best_read_idx;
        }
        self.finalized = true;
        self.buckets = FxHashMap::default();
        self.leaders = LeaderStore::default();
    }

    /// (total reads, clusters), as in the Illumina `DedupContext::stats`.
    pub fn stats(&self) -> (usize, usize) {
        (self.results.len(), self.clusters.len())
    }

    /// Work counters for tuning and diagnostics.
    pub fn detailed_stats(&self) -> OntDedupStats {
        OntDedupStats {
            clusters: self.clusters.len(),
            ..self.stats.clone()
        }
    }

    /// Indices of all exemplar reads (the best read in each cluster).
    /// Must be called after `finalize()`.
    pub fn get_exemplar_indices(&self) -> AHashSet<usize> {
        assert!(
            self.finalized,
            "get_exemplar_indices must be called after finalize()"
        );
        self.clusters
            .iter()
            .map(|c| c.best_read_idx as usize)
            .collect()
    }

    /// Exemplar read index for a read. Must be called after `finalize()`.
    pub fn get_exemplar_index(&self, read_idx: usize) -> usize {
        assert!(
            self.finalized,
            "get_exemplar_index must be called after finalize()"
        );
        self.results[read_idx] as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn offset_allowance_grows_with_distance_from_the_end() {
        let params = OntDedupParams {
            max_end_offset: 2,
            max_end_errors: 5,
            max_error_frac: 0.04,
            ..Default::default()
        };
        // Near the end: endpoint allowance only (2 + 5 + ceil(0.04 * 10)).
        assert!(offsets_compatible(10, 2, &params));
        assert!(!offsets_compatible(10, 1, &params));
        // 250 bases in: plus ceil(0.04 * 250) = 10 bases of indel drift.
        assert!(offsets_compatible(250, 233, &params));
        assert!(!offsets_compatible(250, 232, &params));
        // Symmetric.
        assert!(offsets_compatible(233, 250, &params));
    }

    #[test]
    fn offset_range_is_exactly_the_compatible_offsets() {
        for (frac, base_errors) in [(0.0, 0), (0.04, 5), (0.2, 10), (1.0, 3)] {
            let params = OntDedupParams {
                max_end_offset: 2,
                max_end_errors: base_errors,
                max_error_frac: frac,
                ..Default::default()
            };
            for a in (0..2000).step_by(7) {
                let (lo, hi) = compatible_offset_range(a, &params);
                for b in 0..2600 {
                    let in_range = (lo..=hi).contains(&b);
                    assert_eq!(
                        in_range,
                        offsets_compatible(a, b, &params),
                        "a={a} b={b} frac={frac}"
                    );
                }
            }
        }
    }

    #[test]
    fn misleading_seed_orientation_falls_back_to_the_other_strand() {
        // Leader = P + M + rc(Q) and read = Q + M + rc(P), where Q is P with
        // one base changed inside each of P's start minimizers. The two are
        // forward duplicates (a few substitutions at each end), but no start
        // seed is shared in the forward orientation, while the read's start
        // (Q) matches the leader's end (rc(Q)) exactly. So all seed
        // evidence points to the reverse complement.
        use crate::ont::reverse_complement;
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(1);
        let mut random =
            |n: usize| -> Vec<u8> { (0..n).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect() };
        let minimizer_params = OntMinimizerParams::default();
        let p = random(256);
        let m = random(1500);
        let mut q = p.clone();
        for seed in &end_seeds(&p, &minimizer_params, &mut Vec::new()).start {
            let pos = usize::from(seed.offset) + 7;
            q[pos] = if q[pos] == b'A' { b'C' } else { b'A' };
        }
        let leader = [&p[..], &m, &reverse_complement(&q)].concat();
        let read = [&q[..], &m, &reverse_complement(&p)].concat();
        for rc_tolerant in [true, false] {
            let params = OntDedupParams {
                rc_tolerant,
                ..Default::default()
            };
            let mut ctx = OntDedupContext::new(params, minimizer_params.clone());
            let qual = vec![b'5'; leader.len()];
            ctx.process_read_by_index(0, &leader, &qual);
            ctx.process_read_by_index(1, &read, &qual);
            ctx.finalize();
            // Found via reverse-complement evidence, matched forward, with
            // or without RC tolerance.
            assert_eq!(ctx.stats(), (2, 1), "rc_tolerant={rc_tolerant}");
        }
    }

    #[test]
    fn bucket_lookup_returns_exactly_the_compatible_postings() {
        let params = OntDedupParams::default();
        let posting = |offset| Posting {
            cluster: u32::from(offset),
            offset,
            at_end: false,
        };
        let bucket: Vec<Posting> = (0..400).map(posting).collect();
        for a in [0, 1, 50, 200, 399] {
            let found: Vec<u16> = compatible_postings(&bucket, a, &params)
                .iter()
                .map(|p| p.offset)
                .collect();
            let expected: Vec<u16> = (0..400)
                .filter(|&b| offsets_compatible(a, b, &params))
                .collect();
            assert_eq!(found, expected, "a={a}");
        }
    }

    #[test]
    fn evidence_resets_when_its_stamp_wraps() {
        let mut evidence = Evidence::default();
        evidence.start_read(2);
        *evidence.bits(0) = START_START;
        // Force the stamp to wrap on the next read, then again after a cycle.
        evidence.stamp = u32::MAX;
        evidence.entries[1] = (1, END_END); // would look current after a naive wrap
        evidence.start_read(2);
        assert_eq!(evidence.stamp, 1);
        assert_eq!(*evidence.bits(0), 0);
        assert_eq!(*evidence.bits(1), 0);
    }

    #[test]
    #[should_panic(expected = "max_error_frac must be between 0 and 1")]
    fn invalid_error_fraction_is_rejected() {
        let params = OntDedupParams {
            max_error_frac: -0.1,
            ..Default::default()
        };
        OntDedupContext::new(params, OntMinimizerParams::default());
    }
}
