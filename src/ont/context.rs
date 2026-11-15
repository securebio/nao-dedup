//! Streaming ONT deduplication, mirroring the Illumina `DedupContext`.
//!
//! Reads arrive one at a time. Each read is compared only against cluster
//! leaders that share a seed-pair bucket with it. It joins the first leader
//! it matches; otherwise it starts a new cluster and becomes that cluster's
//! leader. Leaders never change, so cluster membership is decided against a
//! fixed sequence, and clusters never merge. The best read of each cluster
//! (by quality, then length) is tracked and resolved in `finalize()`.

use super::compare::{OntDedupParams, is_duplicate};
use super::reverse_complement_into;
use super::seeds::{EndSeeds, OntMinimizerParams, end_seeds};
use ahash::AHashSet;
use rustc_hash::FxHashMap;

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

/// Bounds `(lo, hi)` containing every offset compatible with `a` under
/// `offsets_compatible` (which is still checked exactly, since the drift
/// allowance depends on the larger offset).
fn compatible_offset_range(a: u16, params: &OntDedupParams) -> (u16, u16) {
    let base = (params.max_end_offset + params.max_end_errors) as f64;
    let a = f64::from(a);
    let frac = params.max_error_frac;
    // b <= a: a - b <= base + ceil(frac * a).
    let lo = a - base - (frac * a).ceil();
    // b >= a: b - a <= base + ceil(frac * b) < base + frac * b + 1, which
    // bounds b only when frac < 1.
    let hi = if frac < 1.0 {
        ((a + base + 1.0) / (1.0 - frac)).ceil()
    } else {
        f64::from(u16::MAX)
    };
    (lo.max(0.0) as u16, hi.min(f64::from(u16::MAX)) as u16)
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
    /// Reads that yielded no seed-pair keys (too short or no valid k-mers);
    /// each becomes a singleton cluster.
    pub unseeded_reads: usize,
    /// Bucket entries examined, summed over reads.
    pub postings_visited: usize,
    /// Leaders with seed evidence at both ends that were aligned against,
    /// summed over reads.
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
    buckets: FxHashMap<u64, Vec<Posting>>,
    /// Leader sequences, concatenated; leader `c` is
    /// `leader_bases[leader_starts[c]..leader_starts[c + 1]]`.
    leader_bases: Vec<u8>,
    leader_starts: Vec<usize>,

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
    /// Candidate leader -> evidence bits, and leaders ready to align.
    evidence: FxHashMap<u32, u8>,
    candidates: Vec<u32>,
}

impl OntDedupContext {
    /// Panics if `minimizer_params` would fail `OntMinimizerParams::new`.
    pub fn new(dedup_params: OntDedupParams, minimizer_params: OntMinimizerParams) -> Self {
        let p = &minimizer_params;
        if let Err(e) = OntMinimizerParams::new(p.kmer_len, p.window_len, p.num_windows) {
            panic!("invalid OntMinimizerParams: {e}");
        }
        Self {
            dedup_params,
            minimizer_params,
            buckets: FxHashMap::default(),
            leader_bases: Vec::new(),
            leader_starts: vec![0],
            clusters: Vec::new(),
            results: Vec::new(),
            stats: OntDedupStats::default(),
            finalized: false,
            hash_scratch: Vec::new(),
            read_scratch: Vec::new(),
            rc_scratch: Vec::new(),
            evidence: FxHashMap::default(),
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
                self.leader_bases.extend_from_slice(&read);
                self.leader_starts.push(self.leader_bases.len());
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
    /// evidence at both ends in a consistent orientation it is aligned, in
    /// that orientation only, and the search stops at the first match. Seeds
    /// near the ends are the most precise, and in dense groups most reads
    /// match early, which avoids scanning the large buckets of inner seeds.
    fn find_matching_cluster(&mut self, read: &[u8], seeds: &EndSeeds) -> Option<u32> {
        let params = &self.dedup_params;
        let evidence = &mut self.evidence;
        let ready = &mut self.candidates;
        evidence.clear();
        let mut rc_ready = false;
        for window in 0..seeds.start.len().max(seeds.end.len()) {
            ready.clear();
            for (read_at_end, end) in [(false, &seeds.start), (true, &seeds.end)] {
                let Some(seed) = end.get(window) else {
                    continue;
                };
                let Some(bucket) = self.buckets.get(&seed.hash) else {
                    continue;
                };
                let (lo, hi) = compatible_offset_range(seed.offset, params);
                let first = bucket.partition_point(|p| p.offset < lo);
                let last = bucket.partition_point(|p| p.offset <= hi);
                let in_range = &bucket[first..last];
                self.stats.postings_visited += in_range.len();
                for posting in in_range {
                    if !offsets_compatible(posting.offset, seed.offset, params) {
                        continue;
                    }
                    let bit = match (read_at_end, posting.at_end) {
                        (false, false) => START_START,
                        (true, true) => END_END,
                        (false, true) => START_END,
                        (true, false) => END_START,
                    };
                    let bits = evidence.entry(posting.cluster).or_insert(0);
                    let before = *bits;
                    *bits |= bit;
                    let completed = |b: u8, mask: u8| b & mask != mask && *bits & mask == mask;
                    if completed(before, SAME_STRAND)
                        || (params.rc_tolerant && completed(before, OPPOSITE_STRAND))
                    {
                        ready.push(posting.cluster);
                    }
                }
            }

            for &cluster in ready.iter() {
                // A leader completing both orientations in one window is
                // listed twice; the second visit finds nothing new to try.
                let bits = evidence.get_mut(&cluster).unwrap();
                let same = *bits & (SAME_STRAND | TRIED_SAME) == SAME_STRAND;
                let opposite = params.rc_tolerant
                    && *bits & (OPPOSITE_STRAND | TRIED_OPPOSITE) == OPPOSITE_STRAND;
                if same {
                    *bits |= TRIED_SAME;
                }
                if opposite {
                    *bits |= TRIED_OPPOSITE;
                }
                if !same && !opposite {
                    continue;
                }
                self.stats.candidates_compared += 1;
                let c = cluster as usize;
                let leader = &self.leader_bases[self.leader_starts[c]..self.leader_starts[c + 1]];
                if same && is_duplicate(read, leader, params) {
                    return Some(cluster);
                }
                if opposite {
                    // Computed once per read, and only if needed.
                    if !rc_ready {
                        reverse_complement_into(read, &mut self.rc_scratch);
                        rc_ready = true;
                    }
                    if is_duplicate(&self.rc_scratch, leader, params) {
                        return Some(cluster);
                    }
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
        self.leader_bases = Vec::new();
        self.leader_starts = Vec::new();
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
    fn offset_range_contains_every_compatible_offset() {
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
                    if offsets_compatible(a, b, &params) {
                        assert!((lo..=hi).contains(&b), "a={a} b={b} frac={frac}");
                    }
                }
            }
        }
    }
}
