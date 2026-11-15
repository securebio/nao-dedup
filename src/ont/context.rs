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
use super::seeds::{OntMinimizerParams, PairKey, end_seeds};
use ahash::AHashSet;
use rustc_hash::{FxHashMap, FxHashSet};

/// A leader's entry in a bucket, with the offsets of the bucket key's seeds
/// in the leader (see `PairKey`).
#[derive(Debug, Clone, Copy)]
struct Posting {
    cluster: u32,
    offsets: (u16, u16),
}

#[derive(Debug, Clone)]
struct Cluster {
    best_read_idx: u32,
    best_score: Score,
}

/// Exemplar ranking, as in the Illumina path: mean Phred quality first,
/// then length. Compared lexicographically.
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
    /// Distinct leaders compared against, summed over reads.
    pub candidates_compared: usize,
}

pub struct OntDedupContext {
    dedup_params: OntDedupParams,
    minimizer_params: OntMinimizerParams,

    /// Seed-pair key -> leaders with that key.
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
    checked: FxHashSet<u32>,
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
            checked: FxHashSet::default(),
        }
    }

    /// Process one read. `idx` must be the read's position in the input:
    /// indices start at 0, each is used once, and they should be contiguous
    /// (sparse indices work but waste memory). Quality is Phred+33.
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

        let keys = end_seeds(&read, &self.minimizer_params, &mut self.hash_scratch).pair_keys();
        if keys.is_empty() {
            self.stats.unseeded_reads += 1;
        }
        let matched = self.find_matching_cluster(&read, &keys);

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
                for key in &keys {
                    self.buckets.entry(key.key).or_default().push(Posting {
                        cluster: cluster_id,
                        offsets: key.offsets,
                    });
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

    /// First leader, in bucket order, that `read` duplicates in an allowed
    /// orientation.
    fn find_matching_cluster(&mut self, read: &[u8], keys: &[PairKey]) -> Option<u32> {
        let mut checked = std::mem::take(&mut self.checked);
        checked.clear();
        // A seed's offset differs between duplicates by the endpoint
        // displacement plus net indels between the seed and the read end.
        // Leaders further off than the endpoint allowance plus the whole-read
        // error rate over that distance are skipped without aligning.
        let params = &self.dedup_params;
        let base_shift = params.max_end_offset + params.max_end_errors;
        let near = |a: u16, b: u16| {
            let drift = (params.max_error_frac * f64::from(a.max(b))).ceil() as usize;
            usize::from(a.abs_diff(b)) <= base_shift + drift
        };
        let mut rc_ready = false;
        let mut matched = None;
        'search: for key in keys {
            let Some(bucket) = self.buckets.get(&key.key) else {
                continue;
            };
            for posting in bucket {
                if !near(posting.offsets.0, key.offsets.0)
                    || !near(posting.offsets.1, key.offsets.1)
                    || !checked.insert(posting.cluster)
                {
                    continue;
                }
                self.stats.candidates_compared += 1;
                let c = posting.cluster as usize;
                let leader = &self.leader_bases[self.leader_starts[c]..self.leader_starts[c + 1]];
                let mut rc_matches = || {
                    params.rc_tolerant && {
                        // Computed once per read, and only if needed.
                        if !rc_ready {
                            reverse_complement_into(read, &mut self.rc_scratch);
                            rc_ready = true;
                        }
                        is_duplicate(&self.rc_scratch, leader, params)
                    }
                };
                if is_duplicate(read, leader, params) || rc_matches() {
                    matched = Some(posting.cluster);
                    break 'search;
                }
            }
        }
        self.checked = checked;
        matched
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
        self.stats.clusters = self.clusters.len();
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
