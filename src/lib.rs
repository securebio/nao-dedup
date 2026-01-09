use ahash::{AHashMap, AHashSet};
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

// ============================================================================
// Configuration Parameters
// ============================================================================

/// Parameters for deduplication matching.
#[derive(Debug, Clone)]
pub struct DedupParams {
    pub max_offset: usize,
    pub max_error_frac: f64,
}

impl Default for DedupParams {
    fn default() -> Self {
        Self {
            max_offset: 1,
            max_error_frac: 0.01,
        }
    }
}

/// Parameters for minimizer extraction.
#[derive(Debug, Clone)]
pub struct MinimizerParams {
    pub kmer_len: usize,
    pub window_len: usize,
    pub num_windows: usize,
}

impl MinimizerParams {
    /// Create new MinimizerParams with validation.
    pub fn new(kmer_len: usize, window_len: usize, num_windows: usize) -> Result<Self, String> {
        if kmer_len > 32 {
            return Err(format!(
                "k-mer length must be <= 32 for 2-bit encoding (fits in 64 bits), got {}",
                kmer_len
            ));
        }
        if kmer_len > window_len {
            return Err(format!(
                "kmer_len ({}) must be <= window_len ({})",
                kmer_len, window_len
            ));
        }
        Ok(Self {
            kmer_len,
            window_len,
            num_windows,
        })
    }
}

impl Default for MinimizerParams {
    fn default() -> Self {
        Self {
            kmer_len: 15,
            window_len: 25,
            num_windows: 4,
        }
    }
}

// ============================================================================
// Read Pair
// ============================================================================

#[derive(Debug, Clone)]
pub struct ReadPair {
    pub read_id: String,
    pub fwd_seq: String,
    pub rev_seq: String,
    pub fwd_qual: String,
    pub rev_qual: String,
}

/// Calculate mean quality score from forward and reverse quality strings.
fn mean_quality(fwd_qual: &str, rev_qual: &str) -> f64 {
    let total: u32 = fwd_qual.bytes().chain(rev_qual.bytes())
        .map(|b| (b - 33) as u32)
        .sum();
    let count = (fwd_qual.len() + rev_qual.len()) as f64;
    if count == 0.0 { 0.0 } else { total as f64 / count }
}

impl ReadPair {
    pub fn mean_quality(&self) -> f64 {
        mean_quality(&self.fwd_qual, &self.rev_qual)
    }
}

/// Lightweight representation of an exemplar for similarity checking.
/// Only stores sequences (not quality strings) to reduce memory footprint.
#[derive(Clone)]
struct StoredExemplar {
    fwd_seq: String,
    rev_seq: String,
}

/// ID registry for interning read IDs to compact u32 indices.
/// Dramatically reduces memory usage and improves hash/comparison performance.
struct IDRegistry {
    id_to_index: AHashMap<String, u32>,
    index_to_id: Vec<String>,
}

impl IDRegistry {
    fn new() -> Self {
        Self {
            id_to_index: AHashMap::new(),
            index_to_id: Vec::new(),
        }
    }

    /// Get or create an index for a read ID
    fn get_or_create(&mut self, id: &str) -> u32 {
        if let Some(&idx) = self.id_to_index.get(id) {
            return idx;
        }
        let idx = self.index_to_id.len() as u32;
        self.index_to_id.push(id.to_string());
        self.id_to_index.insert(id.to_string(), idx);
        idx
    }

    /// Convert index back to read ID
    #[inline]
    fn get_id(&self, idx: u32) -> &str {
        &self.index_to_id[idx as usize]
    }
}

// ============================================================================
// Minimizer Extraction
//
// Strategy: Use 2-bit encoding (A=0,C=1,G=2,T=3) to pack k-mers into u64.
// This allows fast rolling hash computation and comparison.
// ============================================================================

// Lookup table for base encoding (faster than match statement)
// Maps ASCII byte values to 2-bit encodings: A/a=0, C/c=1, G/g=2, T/t=3
// Invalid bases (including N) are marked with u64::MAX
const ENCODE_LOOKUP: [u64; 256] = {
    let mut table = [u64::MAX; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

#[inline(always)]
fn encode_base(b: u8) -> Option<u64> {
    let encoded = ENCODE_LOOKUP[b as usize];
    if encoded == u64::MAX {
        None
    } else {
        Some(encoded)
    }
}

/// Extract one minimizer per window from a sequence.
///
/// Uses rolling hash to efficiently compute k-mer hashes.
/// Windows are adjacent starting from position 0, focusing on the most reliable
/// portion of the read (quality typically degrades toward the end).
fn extract_minimizers(seq: &str, params: &MinimizerParams) -> SmallVec<[u64; 8]> {
    let seq_bytes = seq.as_bytes();
    let seq_len = seq_bytes.len();

    if seq_len < params.kmer_len {
        return SmallVec::new();
    }

    // Mask to keep only the rightmost k*2 bits (each base uses 2 bits)
    let mask: u64 = if params.kmer_len == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * params.kmer_len)) - 1
    };

    let mut minimizers = SmallVec::with_capacity(params.num_windows);


    // Use adjacent windows starting from the beginning of the read
    for i in 0..params.num_windows {
        let window_start = i * params.window_len;

        if window_start >= seq_len {
            break;
        }

        let window_end = (window_start + params.window_len).min(seq_len);

        // Stop if this window would start beyond the sequence
        if window_end - window_start < params.kmer_len {
            break;
        }

        let mut min_hash = u64::MAX;
        let mut hash: u64 = 0;
        let mut valid_len: usize = 0;  // number of consecutive valid ACGT bases

        for pos in window_start..window_end {
            if let Some(encoded) = encode_base(seq_bytes[pos]) {
                // Update forward hash: shift left, add new base
                hash = ((hash << 2) | encoded) & mask;
                valid_len += 1;

                if valid_len >= params.kmer_len {
                    min_hash = min_hash.min(hash);
                }
            } else {
                // Non-ACGT base: reset
                hash = 0;
                valid_len = 0;
            }
        }

        if min_hash != u64::MAX {
            minimizers.push(min_hash);
        }
    }

    minimizers
}

// ============================================================================
// Similarity Checking
//
// Allow sequences to match with small alignment shifts (indels) and mismatches.
// The offset counts as error: e.g., 1bp offset + 1 mismatch = 2 errors total.
// ============================================================================

fn check_similarity(
    seq1: &str,
    seq2: &str,
    max_offset: usize,
    max_error_frac: f64,
) -> bool {
    let s1 = seq1.as_bytes();
    let s2 = seq2.as_bytes();

    // Optimized helper function with early exit for hot path performance
    #[inline]
    fn check_one_way(seqa: &[u8], seqb: &[u8], off: usize, max_error_frac: f64) -> bool {
        if off >= seqa.len() {
            return false;
        }
        let overlap_len = (seqa.len() - off).min(seqb.len());
        if overlap_len == 0 {
            return false;
        }

        // Pre-calculate error budget to avoid floating-point division in the loop
        let max_errors = (max_error_frac * overlap_len as f64).floor() as usize;
        if off > max_errors {
            return false; // Offset alone exceeds budget
        }

        let allowed_mismatches = max_errors - off;
        let mut mismatches = 0;

        let a_part = &seqa[off..off + overlap_len];
        let b_part = &seqb[..overlap_len];

        // Manual loop for early exit when error budget is exceeded
        for i in 0..overlap_len {
            if a_part[i] != b_part[i] {
                mismatches += 1;
                if mismatches > allowed_mismatches {
                    return false; // Short-circuit early!
                }
            }
        }
        true
    }

    for offset in 0..=max_offset {
        // Check with s1 shifted left relative to s2
        if check_one_way(s1, s2, offset, max_error_frac) {
            return true;
        }
        // Check with s2 shifted left relative to s1 (equivalent to s1 shifted
        // right)
        if offset > 0 && check_one_way(s2, s1, offset, max_error_frac) {
            return true;
        }
    }

    false
}

/// Check if two read pairs are similar enough to be duplicates.
///
/// Checks two orientations:
/// 1. Standard: (Fwd, Rev) vs (Fwd, Rev)
/// 2. Swapped: (Fwd, Rev) vs (Rev, Fwd)
///
/// The swapped check handles the case where adapters attached in the opposite
/// orientation, causing the same DNA fragment to be sequenced with forward/reverse
/// swapped.
fn reads_are_similar(
    fwd_seq: &str,
    rev_seq: &str,
    exemplar: &StoredExemplar,
    dedup_params: &DedupParams,
) -> bool {
    if check_similarity(fwd_seq, &exemplar.fwd_seq, dedup_params.max_offset, dedup_params.max_error_frac)
        && check_similarity(rev_seq, &exemplar.rev_seq, dedup_params.max_offset, dedup_params.max_error_frac)
    {
        return true;
    }

    if check_similarity(fwd_seq, &exemplar.rev_seq, dedup_params.max_offset, dedup_params.max_error_frac)
        && check_similarity(rev_seq, &exemplar.fwd_seq, dedup_params.max_offset, dedup_params.max_error_frac)
    {
        return true;
    }

    false
}

// ============================================================================
// Deduplication Context
//
// Streaming algorithm: processes reads one at a time, storing only unique
// sequences (exemplars) rather than all reads. Memory usage scales with
// unique sequences, not total input size.
//
// Key invariant: Each cluster is identified by the read_id of its FIRST member
// (the initial exemplar). As we see better reads, we update best_read_id in
// ClusterStats, but the cluster key remains unchanged. This is crucial for
// lookups to work correctly.
// ============================================================================

#[derive(Debug, Clone)]
struct ClusterStats {
    best_read_idx: u32,  // Index of best read (can change as we see better reads)
    best_score: f64,
    count: usize,
}

pub struct DedupContext {
    dedup_params: DedupParams,
    minimizer_params: MinimizerParams,

    // ID interning: read IDs -> compact u32 indices for faster hashing/comparison
    id_registry: IDRegistry,

    // HashMap choices:
    // - FxHashMap for integer keys: u64/u32 keys are well-distributed
    //   integers, so we can use the ultra-fast FxHash (just a multiply + XOR)

    // minimizer -> list of read indices (instead of read IDs)
    buckets: FxHashMap<u64, Vec<u32>>,

    // read_idx -> read sequences (only for exemplars, quality strings omitted)
    exemplar_store: Vec<Option<StoredExemplar>>,

    // read_idx -> cluster_leader_idx (grows linearly with total reads)
    results: Vec<u32>,

    // cluster_leader_idx -> ClusterStats
    clusters: FxHashMap<u32, ClusterStats>,

    finalized: bool,
}

impl DedupContext {
    pub fn new(dedup_params: DedupParams, minimizer_params: MinimizerParams) -> Self {
        Self {
            dedup_params,
            minimizer_params,
            id_registry: IDRegistry::new(),
            buckets: FxHashMap::default(),
            exemplar_store: Vec::new(),
            results: Vec::new(),
            clusters: FxHashMap::default(),
            finalized: false,
        }
    }

    /// Process one read pair by index (more efficient for indexed reads).
    ///
    /// This method is optimized for cases where reads are naturally indexed
    /// (e.g., from enumerate), avoiding string allocation and hash lookups.
    ///
    /// Algorithm:
    /// 1. Extract minimizers and look up matching exemplars in buckets
    /// 2. Compare this read against candidates until we find a match
    /// 3a. If match found: add to existing cluster, potentially updating best_read_idx
    /// 3b. If no match: create new cluster with this read as initial exemplar
    ///
    /// Returns the cluster leader index.
    pub fn process_read_by_index(
        &mut self,
        idx: usize,
        fwd_seq: String,
        rev_seq: String,
        fwd_qual: String,
        rev_qual: String,
    ) -> u32 {
        let read_idx = idx as u32;

        // Calculate mean quality
        let mean_q = mean_quality(&fwd_qual, &rev_qual);

        // Calculate score: quality is primary (scaled by 1000), length is secondary
        let length = (fwd_seq.len() + rev_seq.len()) as f64;
        let score = mean_q * 1000.0 + length;

        let fwd_mins = extract_minimizers(&fwd_seq, &self.minimizer_params);
        let rev_mins = extract_minimizers(&rev_seq, &self.minimizer_params);

        let mut all_mins = fwd_mins;
        all_mins.extend(rev_mins);

        // Track which candidates we've already checked
        let mut checked_indices = AHashSet::new();
        let mut matching_cluster_idx: Option<u32> = None;

        'outer: for &min_hash in &all_mins {
            if let Some(bucket_reads) = self.buckets.get(&min_hash) {
                for &candidate_idx in bucket_reads {
                    if !checked_indices.insert(candidate_idx) {
                        continue;  // Already checked this candidate
                    }

                    if let Some(candidate) = self.exemplar_store.get(candidate_idx as usize).and_then(|opt| opt.as_ref()) {
                        if reads_are_similar(&fwd_seq, &rev_seq, candidate, &self.dedup_params) {
                            // candidate_idx from buckets is always a cluster leader
                            matching_cluster_idx = Some(candidate_idx);
                            break 'outer;
                        }
                    }
                }
            }
        }

        let cluster_leader_idx = if let Some(cluster_idx) = matching_cluster_idx {
            // Found a match - add to existing cluster
            if let Some(cluster) = self.clusters.get_mut(&cluster_idx) {
                cluster.count += 1;
                // Update best read if this one is better
                if score > cluster.best_score {
                    cluster.best_read_idx = read_idx;
                    cluster.best_score = score;
                }
            }
            cluster_idx
        } else {
            // New unique sequence - create new cluster with this read as leader
            self.clusters.insert(
                read_idx,
                ClusterStats {
                    best_read_idx: read_idx,
                    best_score: score,
                    count: 1,
                },
            );

            // Ensure exemplar_store has space for this index
            if self.exemplar_store.len() <= read_idx as usize {
                self.exemplar_store.resize(read_idx as usize + 1, None);
            }

            // Store only sequences (not quality strings) to reduce memory footprint
            self.exemplar_store[read_idx as usize] = Some(StoredExemplar {
                fwd_seq,
                rev_seq,
            });

            // Add read to minimizer buckets (only for new exemplars)
            for &min_hash in &all_mins {
                self.buckets.entry(min_hash).or_insert_with(Vec::new).push(read_idx);
            }

            read_idx
        };

        // Track this read's cluster assignment
        if self.results.len() <= read_idx as usize {
            self.results.resize(read_idx as usize + 1, 0);
        }
        self.results[read_idx as usize] = cluster_leader_idx;

        cluster_leader_idx
    }

    /// Process one read pair. Returns the cluster ID it was assigned to.
    ///
    /// This method interns the read ID to an index and delegates to
    /// process_read_by_index. See that method for algorithm details.
    pub fn process_read(&mut self, read_pair: ReadPair) -> String {
        // Intern the read ID to a compact u32 index
        let read_idx = self.id_registry.get_or_create(&read_pair.read_id);

        // Delegate to the index-based implementation
        let cluster_leader_idx = self.process_read_by_index(
            read_idx as usize,
            read_pair.fwd_seq,
            read_pair.rev_seq,
            read_pair.fwd_qual,
            read_pair.rev_qual,
        );

        // Return the cluster leader's ID (as a String)
        self.id_registry.get_id(cluster_leader_idx).to_string()
    }

    /// Finalize results: resolve all reads to their cluster's best_read_idx.
    ///
    /// During streaming, reads point to the cluster leader index (first exemplar).
    /// After finalization, they point to the best exemplar found for that cluster.
    pub fn finalize(&mut self) {
        // Update each read's cluster assignment to point to the best exemplar
        for read_idx in 0..self.results.len() {
            let cluster_leader_idx = self.results[read_idx];
            if let Some(cluster) = self.clusters.get(&cluster_leader_idx) {
                self.results[read_idx] = cluster.best_read_idx;
            }
        }

        self.finalized = true;

        // Free memory for intermediate data structures
        self.buckets.clear();
        self.exemplar_store.clear();
    }

    pub fn get_cluster_id(&self, read_id: &str) -> String {
        if let Some(&read_idx) = self.id_registry.id_to_index.get(read_id) {
            if let Some(&cluster_idx) = self.results.get(read_idx as usize) {
                return self.id_registry.get_id(cluster_idx).to_string();
            }
        }
        read_id.to_string()
    }

    pub fn stats(&self) -> (usize, usize) {
        let total_reads = self.results.len();
        let unique_clusters = self.clusters.len();
        (total_reads, unique_clusters)
    }

    /// Get indices of all exemplar reads (the best read in each cluster).
    ///
    /// Must be called after finalize(). Returns a HashSet of read indices where
    /// each index represents an exemplar (a read that is the best in its cluster).
    ///
    /// This is much more efficient than calling get_cluster_id in a loop, as it
    /// directly collects the best_read_idx from each cluster.
    pub fn get_exemplar_indices(&self) -> AHashSet<usize> {
        assert!(
            self.finalized,
            "get_exemplar_indices must be called after finalize()"
        );
        self.clusters
            .values()
            .map(|stats| stats.best_read_idx as usize)
            .collect()
    }
}

// ============================================================================
// Public API
// ============================================================================

pub fn deduplicate_read_pairs(
    read_pairs: Vec<ReadPair>,
    dedup_params: Option<DedupParams>,
    minimizer_params: Option<MinimizerParams>,
) -> AHashMap<String, String> {
    let dedup_params = dedup_params.unwrap_or_default();
    let minimizer_params = minimizer_params.unwrap_or_default();

    let mut ctx = DedupContext::new(dedup_params, minimizer_params);

    for rp in read_pairs {
        ctx.process_read(rp);
    }

    ctx.finalize();

    // Convert internal index-based results to HashMap<String, String> for public API
    let mut result_map = AHashMap::with_capacity(ctx.results.len());
    for read_idx in 0..ctx.results.len() {
        let cluster_idx = ctx.results[read_idx];
        let read_id = ctx.id_registry.get_id(read_idx as u32);
        let cluster_id = ctx.id_registry.get_id(cluster_idx);
        result_map.insert(read_id.to_string(), cluster_id.to_string());
    }

    result_map
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    // Sentinel value for empty k-mer windows
    const EMPTY_KMER_SENTINEL_HASH: u64 = u64::MAX;

    /// Generate a random DNA sequence of specified length
    fn random_seq(length: usize, rng: &mut StdRng) -> String {
        const BASES: &[u8] = b"ACGT";
        (0..length)
            .map(|_| BASES[rng.gen_range(0..4)] as char)
            .collect()
    }

    /// Reverse complement a DNA sequence
    fn reverse_complement(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                'A' | 'a' => 'T',
                'T' | 't' => 'A',
                'C' | 'c' => 'G',
                'G' | 'g' => 'C',
                'N' | 'n' => 'N',
                _ => c,
            })
            .collect()
    }

    /// Count mismatches between two sequences (compares only up to shorter length)
    fn mismatch_count(seq1: &str, seq2: &str) -> usize {
        seq1.chars()
            .zip(seq2.chars())
            .filter(|(a, b)| a != b)
            .count()
    }

    // ========================================================================
    // Helper Functions Tests
    // ========================================================================

    #[test]
    fn test_reverse_complement_standard_bases() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("TTTT"), "AAAA");
        assert_eq!(reverse_complement("GCGC"), "GCGC");
    }

    #[test]
    fn test_reverse_complement_with_n() {
        assert_eq!(reverse_complement("ACGTN"), "NACGT");
        assert_eq!(reverse_complement("NNNNN"), "NNNNN");
    }

    #[test]
    fn test_mismatch_count_equal_length() {
        assert_eq!(mismatch_count("AAAA", "AAAA"), 0);
        assert_eq!(mismatch_count("AAAA", "TTTT"), 4);
        assert_eq!(mismatch_count("AAAA", "AAAT"), 1);
        assert_eq!(mismatch_count("ACGT", "TGCA"), 4);
    }

    #[test]
    fn test_mismatch_count_unequal_length() {
        // Truncates to shorter length
        assert_eq!(mismatch_count("AAAA", "AA"), 0); // Only compares first 2
        assert_eq!(mismatch_count("AA", "AAAA"), 0); // Only compares first 2
        assert_eq!(mismatch_count("AAAA", "TT"), 2); // Compares first 2, both differ
        assert_eq!(mismatch_count("AAAT", "TT"), 2); // AA vs TT
    }

    // ========================================================================
    // Minimizer Extraction Tests
    // ========================================================================

    #[test]
    fn test_extract_minimizer_normal_window() {
        let params = MinimizerParams::new(7, 20, 2).unwrap();
        let seq = format!("{}{}", "A".repeat(20), "C".repeat(20)); // 40bp sequence

        // Test first window (all A's - should give consistent result)
        let minimizers1 = extract_minimizers(&seq, &params);
        let minimizers2 = extract_minimizers(&seq, &params);
        assert_eq!(minimizers1.len(), 2);
        assert_eq!(minimizers1[0], minimizers2[0]); // Should be deterministic

        // Different windows should give different results
        assert_ne!(minimizers1[0], minimizers1[1]);
    }

    #[test]
    fn test_extract_minimizer_with_n_bases() {
        let params = MinimizerParams::new(3, 10, 1).unwrap();
        let seq_with_n = "AANAAAANAA";
        let seq_without_n = "AAGAAAAGAA";

        let mins_with_n = extract_minimizers(seq_with_n, &params);
        let mins_without_n = extract_minimizers(seq_without_n, &params);

        // Should skip N-containing kmers and find valid ones
        assert!(!mins_with_n.is_empty());
        assert!(!mins_without_n.is_empty());
    }

    #[test]
    fn test_extract_minimizer_window_too_short() {
        let params = MinimizerParams::new(7, 10, 1).unwrap();
        let seq = "AAAAA"; // 5bp sequence, need 7bp kmer

        let mins = extract_minimizers(seq, &params);
        assert!(mins.is_empty());
    }

    #[test]
    fn test_extract_minimizer_sequence_too_short() {
        // Collected windows longer than sequence, should not add sentinels (Rust behavior)
        let params = MinimizerParams::new(7, 10, 2).unwrap();
        let seq = "AAAAACCGGTT"; // 11bp sequence, second window is too short

        let mins = extract_minimizers(seq, &params);
        // Rust returns only the first window's minimizer, skips the second
        assert_eq!(mins.len(), 1);
    }

    #[test]
    fn test_extract_minimizer_sequence_matches_window_matches_kmer() {
        let params = MinimizerParams::new(11, 11, 1).unwrap();
        let seq = "AAAAACCGGTT"; // 11bp sequence

        let mins = extract_minimizers(seq, &params);
        assert_eq!(mins.len(), 1);
        assert_ne!(mins[0], EMPTY_KMER_SENTINEL_HASH);
    }

    #[test]
    fn test_extract_minimizer_all_n_window() {
        let params = MinimizerParams::new(3, 10, 1).unwrap();
        let seq = "NNNNNNNNNN";

        let mins = extract_minimizers(seq, &params);
        assert!(mins.is_empty()); // Rust skips N-only windows
    }

    // ========================================================================
    // Sequence Matching Tests
    // ========================================================================

    #[test]
    fn test_sequences_match_exact() {
        let params = DedupParams {
            max_offset: 1,
            max_error_frac: 0.01,
        };

        assert!(check_similarity("AAAA", "AAAA", params.max_offset, params.max_error_frac));
        assert!(check_similarity("ACGT", "ACGT", params.max_offset, params.max_error_frac));
    }

    #[test]
    fn test_sequences_match_with_offset_1() {
        let params = DedupParams {
            max_offset: 1,
            max_error_frac: 0.25,
        };
        // need a large max_error_frac; for short test seqs, an offset of 1 is a large relative error

        // Left shift: XAAAA vs AAAA (X removed)
        assert!(check_similarity("GAAAA", "AAAA", params.max_offset, params.max_error_frac));
        // Right shift: AAAA vs XAAAA (X added at start)
        assert!(check_similarity("AAAA", "GAAAA", params.max_offset, params.max_error_frac));
    }

    #[test]
    fn test_sequences_match_no_match_large_offset() {
        let params = DedupParams {
            max_offset: 1,
            max_error_frac: 0.01,
        };

        // Should not match with offset > 1
        assert!(!check_similarity("GGAAAA", "AAAA", params.max_offset, params.max_error_frac));
        assert!(!check_similarity("AAAA", "GGAAAA", params.max_offset, params.max_error_frac));
    }

    #[test]
    fn test_sequences_error_threshold() {
        let params = DedupParams {
            max_offset: 0,
            max_error_frac: 0.1,
        }; // 10% error allowed

        // 1 error in 10bp = 10% error rate
        assert!(check_similarity("AAAAAAAAAA", "AAAAAAAAAG", params.max_offset, params.max_error_frac));
        // 2 errors in 10bp = 20% error rate (should fail)
        assert!(!check_similarity("AAAAAAAAAA", "AAAAAAAGGG", params.max_offset, params.max_error_frac));
    }

    #[test]
    fn test_read_pairs_equivalent_standard_orientation() {
        let params = DedupParams {
            max_offset: 1,
            max_error_frac: 0.01,
        };

        let exemplar1 = StoredExemplar {
            fwd_seq: "AAAA".to_string(),
            rev_seq: "TTTT".to_string(),
        };
        let exemplar2 = StoredExemplar {
            fwd_seq: "AAAA".to_string(),
            rev_seq: "CCCC".to_string(),
        };

        assert!(reads_are_similar("AAAA", "TTTT", &exemplar1, &params));
        assert!(!reads_are_similar("AAAA", "TTTT", &exemplar2, &params));
    }

    #[test]
    fn test_read_pairs_equivalent_no_match() {
        let params = DedupParams {
            max_offset: 1,
            max_error_frac: 0.01,
        };

        let exemplar = StoredExemplar {
            fwd_seq: "GGGG".to_string(),
            rev_seq: "CCCC".to_string(),
        };

        assert!(!reads_are_similar("AAAA", "TTTT", &exemplar, &params));
    }

    // ========================================================================
    // ReadPair Validation Tests
    // ========================================================================

    #[test]
    fn test_mean_qual_calculation() {
        // Phred 33: '!' = 0, 'I' = 40
        let rp = ReadPair {
            read_id: "test".to_string(),
            fwd_seq: "AAAA".to_string(),
            rev_seq: "TTTT".to_string(),
            fwd_qual: "!!!!".to_string(),
            rev_qual: "IIII".to_string(),
        };
        let expected_mean = (0.0 + 0.0 + 0.0 + 0.0 + 40.0 + 40.0 + 40.0 + 40.0) / 8.0;
        assert_eq!(rp.mean_quality(), expected_mean);
    }

    #[test]
    fn test_mean_qual_empty_qualities() {
        let rp = ReadPair {
            read_id: "test".to_string(),
            fwd_seq: "AAAA".to_string(),
            rev_seq: "TTTT".to_string(),
            fwd_qual: "".to_string(),
            rev_qual: "".to_string(),
        };
        assert_eq!(rp.mean_quality(), 0.0);
    }

    // ========================================================================
    // Parameter Validation Tests
    // ========================================================================

    #[test]
    fn test_minimizer_params_kmer_too_large() {
        let result = MinimizerParams::new(7, 5, 1);
        assert!(result.is_err());
        let err_msg = result.unwrap_err();
        assert!(err_msg.contains("kmer_len"));
        assert!(err_msg.contains("window_len"));
    }

    #[test]
    fn test_minimizer_params_valid() {
        let params = MinimizerParams::new(5, 10, 1).unwrap();
        assert_eq!(params.window_len, 10);
        assert_eq!(params.kmer_len, 5);

        // it's legal if surprising to have window and kmer length match
        // Note: Rust has max kmer_len of 32, so we use 32 instead of 47
        let params2 = MinimizerParams::new(32, 32, 1).unwrap();
        assert_eq!(params2.window_len, 32);
        assert_eq!(params2.kmer_len, 32);
    }

    // ========================================================================
    // Bucketing Test (1 test)
    // ========================================================================

    #[test]
    fn test_dups_share_buckets() {
        // Test with realistic-ish data that two read pairs which are duplicates
        // with a single base error in each of fwd/rev read share a bucket.
        // This test is probabilistic, though with fixed random seed should be consistent.
        let mut rng = StdRng::seed_from_u64(1234);
        let read_len = 150;

        fn seq_with_error(seq: &str, rng: &mut StdRng, read_len: usize) -> String {
            let error_loc = rng.gen_range(0..read_len);
            let start = if rng.gen_bool(0.5) { 1 } else { 0 };
            let bases = ['A', 'C', 'G', 'T'];
            let new_base = bases[rng.gen_range(0..4)];

            // Handle the offset case: trim from start if needed
            if start == 1 && seq.len() > 1 && error_loc > 0 {
                // Remove first base, introduce error at error_loc-1 in the remaining sequence
                let adjusted_loc = error_loc.saturating_sub(1).min(seq.len() - 2);
                format!("{}{}{}", &seq[1..=adjusted_loc], new_base, &seq[adjusted_loc + 2..])
            } else if error_loc < seq.len() {
                // No offset, just introduce error at error_loc
                format!("{}{}{}", &seq[0..error_loc], new_base, &seq[error_loc + 1..].trim_end())
            } else {
                seq.to_string()
            }
        }

        for _ in 0..100 {
            // create a pair of reads that are dups but not perfect dups
            let fwd1 = random_seq(read_len, &mut rng);
            let rev1 = random_seq(read_len, &mut rng);

            let fwd2 = seq_with_error(&fwd1, &mut rng, read_len);
            let rev2 = seq_with_error(&rev1, &mut rng, read_len);

            let fwd2_len = fwd2.len();
            let rev2_len = rev2.len();

            let rp1 = ReadPair {
                read_id: "pair1".to_string(),
                fwd_seq: fwd1,
                rev_seq: rev1,
                fwd_qual: "I".repeat(read_len),
                rev_qual: "I".repeat(read_len),
            };
            let rp2 = ReadPair {
                read_id: "pair2".to_string(),
                fwd_seq: fwd2,
                rev_seq: rev2,
                fwd_qual: "I".repeat(fwd2_len),
                rev_qual: "I".repeat(rev2_len),
            };

            let m_params = MinimizerParams::default();

            // Extract minimizers from both reads
            let mins1_fwd = extract_minimizers(&rp1.fwd_seq, &m_params);
            let mins1_rev = extract_minimizers(&rp1.rev_seq, &m_params);
            let mins2_fwd = extract_minimizers(&rp2.fwd_seq, &m_params);
            let mins2_rev = extract_minimizers(&rp2.rev_seq, &m_params);

            // Check if they share at least one minimizer
            let mut all_mins1 = mins1_fwd;
            all_mins1.extend(mins1_rev);
            let mut all_mins2 = mins2_fwd;
            all_mins2.extend(mins2_rev);

            let set1: AHashSet<u64> = all_mins1.iter().copied().collect();
            let set2: AHashSet<u64> = all_mins2.iter().copied().collect();

            let shared = set1.intersection(&set2).count();
            assert!(shared > 0, "Duplicates with errors should share at least one bucket");
        }
    }

    // ========================================================================
    // Base Encoding Tests
    // ========================================================================

    #[test]
    fn test_encode_base_valid_uppercase() {
        assert_eq!(encode_base(b'A'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'G'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
    }

    #[test]
    fn test_encode_base_valid_lowercase() {
        assert_eq!(encode_base(b'a'), Some(0));
        assert_eq!(encode_base(b'c'), Some(1));
        assert_eq!(encode_base(b'g'), Some(2));
        assert_eq!(encode_base(b't'), Some(3));
    }

    #[test]
    fn test_encode_base_invalid_n() {
        assert_eq!(encode_base(b'N'), None);
        assert_eq!(encode_base(b'n'), None);
    }

    #[test]
    fn test_encode_base_invalid_characters() {
        assert_eq!(encode_base(b'X'), None);
        assert_eq!(encode_base(b'Y'), None);
        assert_eq!(encode_base(b'-'), None);
        assert_eq!(encode_base(b' '), None);
        assert_eq!(encode_base(b'0'), None);
    }

    // ========================================================================
    // IDRegistry Tests
    // ========================================================================

    #[test]
    fn test_id_registry_get_or_create_new_ids() {
        let mut registry = IDRegistry::new();

        let idx1 = registry.get_or_create("read1");
        let idx2 = registry.get_or_create("read2");
        let idx3 = registry.get_or_create("read3");

        assert_eq!(idx1, 0);
        assert_eq!(idx2, 1);
        assert_eq!(idx3, 2);
    }

    #[test]
    fn test_id_registry_duplicate_ids_same_index() {
        let mut registry = IDRegistry::new();

        let idx1 = registry.get_or_create("read1");
        let idx2 = registry.get_or_create("read1");
        let idx3 = registry.get_or_create("read1");

        assert_eq!(idx1, idx2);
        assert_eq!(idx2, idx3);
    }

    #[test]
    fn test_id_registry_get_id_retrieval() {
        let mut registry = IDRegistry::new();

        registry.get_or_create("read1");
        registry.get_or_create("read2");
        registry.get_or_create("read3");

        assert_eq!(registry.get_id(0), "read1");
        assert_eq!(registry.get_id(1), "read2");
        assert_eq!(registry.get_id(2), "read3");
    }

    #[test]
    fn test_id_registry_roundtrip() {
        let mut registry = IDRegistry::new();

        let original_ids = vec!["read_A", "read_B", "read_C", "read_D"];

        // Store all IDs and get their indices
        let mut indices = Vec::new();
        for id in &original_ids {
            indices.push(registry.get_or_create(id));
        }

        // Retrieve IDs and verify they match
        for (idx, original_id) in indices.iter().zip(&original_ids) {
            assert_eq!(registry.get_id(*idx), *original_id);
        }
    }

    // ========================================================================
    // DedupContext Public API Tests
    // ========================================================================

    #[test]
    fn test_dedup_context_stats() {
        let dedup_params = DedupParams::default();
        let minimizer_params = MinimizerParams::default();

        let read_pairs = vec![
            ReadPair {
                read_id: "read1".to_string(),
                fwd_seq: "ACGTACGTACGTACGT".to_string(),
                rev_seq: "TGCATGCATGCATGCA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
            ReadPair {
                read_id: "read2".to_string(),
                fwd_seq: "ACGTACGTACGTACGT".to_string(), // Duplicate
                rev_seq: "TGCATGCATGCATGCA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
            ReadPair {
                read_id: "read3".to_string(),
                fwd_seq: "GGGGAAAACCCCTTTT".to_string(), // Different
                rev_seq: "TTTTGGGGCCCCAAAA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
        ];

        let mut ctx = DedupContext::new(dedup_params, minimizer_params);

        for rp in read_pairs {
            ctx.process_read(rp);
        }

        ctx.finalize();

        let (total_reads, unique_clusters) = ctx.stats();
        assert_eq!(total_reads, 3);
        assert_eq!(unique_clusters, 2); // read1 and read2 cluster together, read3 is separate
    }

    #[test]
    fn test_dedup_context_get_cluster_id() {
        let dedup_params = DedupParams::default();
        let minimizer_params = MinimizerParams::default();

        let read_pairs = vec![
            ReadPair {
                read_id: "read1".to_string(),
                fwd_seq: "ACGTACGTACGTACGT".to_string(),
                rev_seq: "TGCATGCATGCATGCA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
            ReadPair {
                read_id: "read2".to_string(),
                fwd_seq: "ACGTACGTACGTACGT".to_string(), // Duplicate
                rev_seq: "TGCATGCATGCATGCA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
        ];

        let mut ctx = DedupContext::new(dedup_params, minimizer_params);

        for rp in read_pairs {
            ctx.process_read(rp);
        }

        ctx.finalize();

        // Both reads should have the same cluster ID
        let cluster1 = ctx.get_cluster_id("read1");
        let cluster2 = ctx.get_cluster_id("read2");
        assert_eq!(cluster1, cluster2);
    }

    #[test]
    fn test_dedup_context_get_cluster_id_unknown_read() {
        let dedup_params = DedupParams::default();
        let minimizer_params = MinimizerParams::default();

        let mut ctx = DedupContext::new(dedup_params, minimizer_params);
        ctx.finalize();

        // Unknown read ID should return the ID itself
        assert_eq!(ctx.get_cluster_id("unknown_read"), "unknown_read");
    }

    #[test]
    fn test_dedup_context_get_exemplar_indices() {
        let dedup_params = DedupParams::default();
        let minimizer_params = MinimizerParams::default();

        let read_pairs = vec![
            ReadPair {
                read_id: "read1".to_string(),
                fwd_seq: "ACGTACGTACGTACGT".to_string(),
                rev_seq: "TGCATGCATGCATGCA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
            ReadPair {
                read_id: "read2".to_string(),
                fwd_seq: "ACGTACGTACGTACGT".to_string(), // Duplicate
                rev_seq: "TGCATGCATGCATGCA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
            ReadPair {
                read_id: "read3".to_string(),
                fwd_seq: "GGGGAAAACCCCTTTT".to_string(), // Different
                rev_seq: "TTTTGGGGCCCCAAAA".to_string(),
                fwd_qual: "IIIIIIIIIIIIIIII".to_string(),
                rev_qual: "IIIIIIIIIIIIIIII".to_string(),
            },
        ];

        let mut ctx = DedupContext::new(dedup_params, minimizer_params);

        for rp in read_pairs {
            ctx.process_read(rp);
        }

        ctx.finalize();

        let exemplar_indices = ctx.get_exemplar_indices();

        // Should have 2 exemplars (one for cluster of read1/read2, one for read3)
        assert_eq!(exemplar_indices.len(), 2);

        // Exemplar indices should be 0 (for read1) and 2 (for read3)
        assert!(exemplar_indices.contains(&0) || exemplar_indices.contains(&1));
        assert!(exemplar_indices.contains(&2));
    }

    #[test]
    #[should_panic(expected = "get_exemplar_indices must be called after finalize()")]
    fn test_dedup_context_get_exemplar_indices_before_finalize() {
        let dedup_params = DedupParams::default();
        let minimizer_params = MinimizerParams::default();

        let ctx = DedupContext::new(dedup_params, minimizer_params);

        // Should panic because finalize() hasn't been called
        ctx.get_exemplar_indices();
    }

    // ========================================================================
    // Default Parameter Tests
    // ========================================================================

    #[test]
    fn test_dedup_params_default() {
        let params = DedupParams::default();

        assert_eq!(params.max_offset, 1);
        assert_eq!(params.max_error_frac, 0.01);
    }

    #[test]
    fn test_minimizer_params_default() {
        let params = MinimizerParams::default();

        assert_eq!(params.kmer_len, 15);
        assert_eq!(params.window_len, 25);
        assert_eq!(params.num_windows, 4);
    }
}
