//! Seeds for bucketing ONT reads.
//!
//! The Illumina path takes one minimizer from each of several adjacent windows
//! at the start of each mate. ONT reads are single-end, but a duplicate is
//! defined by both fragment endpoints, so the two ends of the read play the
//! role of the two mates. Windows are laid out inward from the start and
//! inward from the end.
//!
//! Two differences from the Illumina minimizers:
//! - K-mers are canonical (the smaller of the k-mer and its reverse
//!   complement), and windows at the two ends mirror each other. A read and
//!   its reverse complement therefore yield the same seeds, with start and
//!   end swapped.
//! - Minimizers are chosen by a mixing hash, not by the raw 2-bit value, so
//!   A-rich k-mers (poly-A tails) aren't systematically preferred.
//!   Single-base homopolymer k-mers are never seeds.

use smallvec::SmallVec;

/// Parameters for ONT seed extraction.
#[derive(Debug, Clone)]
pub struct OntMinimizerParams {
    pub kmer_len: usize,
    pub window_len: usize,
    /// Number of windows at each end of the read.
    pub num_windows: usize,
}

impl OntMinimizerParams {
    /// Create new OntMinimizerParams with validation.
    pub fn new(kmer_len: usize, window_len: usize, num_windows: usize) -> Result<Self, String> {
        // k = 1 would leave no seeds: every 1-mer is a homopolymer.
        if !(2..=32).contains(&kmer_len) {
            return Err(format!(
                "k-mer length must be between 2 and 32 for 2-bit encoding, got {kmer_len}"
            ));
        }
        if kmer_len > window_len {
            return Err(format!(
                "kmer_len ({kmer_len}) must be <= window_len ({window_len})"
            ));
        }
        if num_windows == 0 {
            return Err("num_windows must be at least 1".to_string());
        }
        if num_windows.saturating_mul(window_len) > usize::from(u16::MAX) {
            return Err(format!(
                "num_windows * window_len must be at most {}, got {num_windows} * {window_len}",
                u16::MAX
            ));
        }
        Ok(Self {
            kmer_len,
            window_len,
            num_windows,
        })
    }
}

/// Defaults chosen by sweeps on simulated reads of roughly 1% per-read
/// error. Smaller windows give more minimizers per end, which raised
/// retrieval recall more than larger k-mers or fewer, larger windows.
impl Default for OntMinimizerParams {
    fn default() -> Self {
        Self {
            kmer_len: 15,
            window_len: 32,
            num_windows: 8,
        }
    }
}

/// A minimizer and its distance in bases from the read end it was taken
/// from: from the read start to the k-mer start for start windows, from the
/// k-mer end to the read end for end windows. Under reverse complementation
/// both the hash and the offset carry over to the other end.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Seed {
    pub hash: u64,
    pub offset: u16,
}

pub(crate) type Minimizers = SmallVec<[Seed; 8]>;

/// Minimizers from the windows at each end of a read, ordered from the read
/// boundary inward. Each seed is a bucket key; a read and its reverse
/// complement have the same seeds with the ends swapped.
#[derive(Debug, Default, PartialEq)]
pub(crate) struct EndSeeds {
    pub start: Minimizers,
    pub end: Minimizers,
}

const INVALID: u64 = u64::MAX;

// Maps ASCII to 2-bit codes (A/a=0, C/c=1, G/g=2, T/t=3); anything else is
// invalid and breaks k-mers.
const ENCODE: [u8; 256] = {
    let mut table = [4u8; 256];
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

/// SplitMix64 finalizer: a fast bijective mix so that minimizer order is
/// effectively random rather than lexicographic.
#[inline]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// Hash of the canonical k-mer starting at each position of `region`, or
/// INVALID for k-mers containing a non-ACGT base or consisting of one
/// repeated base. `out[i]` describes `region[i..i + k]`.
fn kmer_hashes(region: &[u8], k: usize, out: &mut Vec<u64>) {
    out.clear();
    if region.len() < k {
        return;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    // 0b0101...01 over k bases: the encoding of C...C.
    let ones = mask / 3;
    let rc_shift = 2 * (k as u32 - 1);
    let (mut fwd, mut rc, mut valid) = (0u64, 0u64, 0usize);
    for (pos, &byte) in region.iter().enumerate() {
        let code = ENCODE[byte as usize];
        if code > 3 {
            valid = 0;
        } else {
            fwd = ((fwd << 2) | code as u64) & mask;
            rc = (rc >> 2) | ((3 - code as u64) << rc_shift);
            valid += 1;
        }
        if pos + 1 >= k {
            let canonical = fwd.min(rc);
            // Canonical homopolymers are A...A (0) and C...C (`ones`).
            let hash = if valid >= k && canonical != 0 && canonical != ones {
                // Keep INVALID reserved; the remap merges two hash values.
                mix(canonical).min(INVALID - 1)
            } else {
                INVALID
            };
            out.push(hash);
        }
    }
}

/// Smallest valid hash among k-mers lying entirely within `[start, end)` of
/// the region whose k-mer hashes are `hashes`, with its k-mer's position.
/// Ties go to the k-mer nearest the window's outer edge: the start for start
/// windows, the end (`from_end`) for end windows, so that reverse
/// complementing a read selects the same k-mer.
fn window_min(
    hashes: &[u64],
    start: usize,
    end: usize,
    k: usize,
    from_end: bool,
) -> Option<(u64, usize)> {
    if end < start + k {
        return None;
    }
    let mut best: Option<(u64, usize)> = None;
    let mut consider = |pos: usize| {
        let hash = hashes[pos];
        if hash != INVALID && best.is_none_or(|(b, _)| hash < b) {
            best = Some((hash, pos));
        }
    };
    if from_end {
        (start..=end - k).rev().for_each(&mut consider);
    } else {
        (start..=end - k).for_each(&mut consider);
    }
    best
}

/// Extract the minimizers of `num_windows` adjacent windows inward from each
/// end of `seq`. Windows are clipped at the read boundary, and a window
/// shorter than k ends the sequence of windows at that end. On reads shorter
/// than `num_windows * window_len` the two ends' windows overlap.
pub(crate) fn end_seeds(
    seq: &[u8],
    params: &OntMinimizerParams,
    scratch: &mut Vec<u64>,
) -> EndSeeds {
    let (k, w) = (params.kmer_len, params.window_len);
    let span = (params.num_windows * w).min(seq.len());
    let mut seeds = EndSeeds::default();
    let seed = |hash, offset: usize| Seed {
        hash,
        offset: u16::try_from(offset)
            .expect("seeded span exceeds u16; construct params with OntMinimizerParams::new"),
    };

    kmer_hashes(&seq[..span], k, scratch);
    for i in 0..params.num_windows {
        let (start, end) = (i * w, ((i + 1) * w).min(span));
        match window_min(scratch, start, end, k, false) {
            Some((hash, pos)) => seeds.start.push(seed(hash, pos)),
            None if end < start + k => break,
            None => {}
        }
    }

    kmer_hashes(&seq[seq.len() - span..], k, scratch);
    for i in 0..params.num_windows {
        let end = span.saturating_sub(i * w);
        let start = span.saturating_sub((i + 1) * w);
        match window_min(scratch, start, end, k, true) {
            Some((hash, pos)) => seeds.end.push(seed(hash, span - (pos + k))),
            None if end < start + k => break,
            None => {}
        }
    }
    seeds
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ont::{reverse_complement, test_util::random_seq};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn seeds(seq: &[u8], params: &OntMinimizerParams) -> EndSeeds {
        end_seeds(seq, params, &mut Vec::new())
    }

    #[test]
    fn params_validation() {
        assert!(OntMinimizerParams::new(17, 64, 4).is_ok());
        assert!(OntMinimizerParams::new(33, 64, 4).is_err());
        assert!(OntMinimizerParams::new(0, 64, 4).is_err());
        assert!(OntMinimizerParams::new(1, 64, 4).is_err());
        assert!(OntMinimizerParams::new(65, 64, 4).is_err());
        assert!(OntMinimizerParams::new(17, 64, 0).is_err());
    }

    #[test]
    fn canonical_hash_is_strand_independent() {
        let mut rng = StdRng::seed_from_u64(1);
        let kmer = random_seq(17, &mut rng);
        let (mut a, mut b) = (Vec::new(), Vec::new());
        kmer_hashes(&kmer, 17, &mut a);
        kmer_hashes(&reverse_complement(&kmer), 17, &mut b);
        assert_eq!(a.len(), 1);
        assert_eq!(a, b);
    }

    #[test]
    fn reverse_complement_swaps_ends() {
        let mut rng = StdRng::seed_from_u64(2);
        let params = OntMinimizerParams::default();
        // Includes lengths where the two ends' windows overlap or are clipped.
        for len in [20, 63, 64, 100, 300, 511, 512, 513, 2000] {
            let seq = random_seq(len, &mut rng);
            let fwd = seeds(&seq, &params);
            let rc = seeds(&reverse_complement(&seq), &params);
            assert_eq!(fwd.start, rc.end, "len {len}");
            assert_eq!(fwd.end, rc.start, "len {len}");
        }
    }

    #[test]
    fn equal_minima_resolve_toward_the_read_end() {
        // Each 6-base window of this repeat contains the same 2-mers several
        // times, so the minimum is always tied.
        let params = OntMinimizerParams::new(2, 6, 2).unwrap();
        let seq = b"ACACACACACAC";
        let fwd = seeds(seq, &params);
        let offsets = |m: &Minimizers| m.iter().map(|s| s.offset).collect::<Vec<_>>();
        // Ties go to the occurrence nearest each window's outer edge.
        assert_eq!(offsets(&fwd.start), [0, 6]);
        assert_eq!(offsets(&fwd.end), [0, 6]);
        let rc = seeds(&reverse_complement(seq), &params);
        assert_eq!((fwd.start.clone(), fwd.end.clone()), (rc.end, rc.start));
    }

    #[test]
    fn span_validation_bounds() {
        assert!(OntMinimizerParams::new(15, 65535, 1).is_ok());
        assert!(OntMinimizerParams::new(15, 32768, 2).is_err());
        assert!(OntMinimizerParams::new(15, usize::MAX, 2).is_err());
    }

    #[test]
    fn window_count_and_short_reads() {
        let mut rng = StdRng::seed_from_u64(3);
        let params = OntMinimizerParams::new(17, 64, 4).unwrap();
        let s = seeds(&random_seq(1000, &mut rng), &params);
        assert_eq!((s.start.len(), s.end.len()), (4, 4));
        // 100 bp: one full window and a 36 bp window at each end.
        let s = seeds(&random_seq(100, &mut rng), &params);
        assert_eq!((s.start.len(), s.end.len()), (2, 2));
        // Shorter than k: no seeds.
        let s = seeds(&random_seq(16, &mut rng), &params);
        assert!(s.start.is_empty() && s.end.is_empty());
    }

    #[test]
    fn invalid_and_homopolymer_kmers_are_not_seeds() {
        let params = OntMinimizerParams::new(5, 10, 1).unwrap();
        assert_eq!(seeds(b"AAAAAAAAAA", &params), EndSeeds::default());
        assert_eq!(seeds(b"NNNNNNNNNN", &params), EndSeeds::default());
        // The only valid non-homopolymer 5-mers lie entirely after the N.
        let with_n = seeds(b"AAAANACGTA", &params);
        let tail = seeds(b"ACGTA", &OntMinimizerParams::new(5, 5, 1).unwrap());
        assert_eq!(with_n.start[0].hash, tail.start[0].hash);
        assert_eq!(with_n.start[0].offset, 5);
        assert_eq!(with_n.end[0].offset, 0);
    }

    #[test]
    fn duplicates_with_errors_share_seeds_at_both_ends() {
        // Substitutions and indels near both ends, in either orientation:
        // each end still shares a seed with the original, at an offset
        // within the indel size. For the reverse complement, the read's
        // start corresponds to the original's end.
        let mut rng = StdRng::seed_from_u64(4);
        let params = OntMinimizerParams::default();
        let seq = random_seq(1500, &mut rng);
        let mut read = seq.clone();
        read[10] = if seq[10] == b'A' { b'C' } else { b'A' };
        read.insert(70, b'G'); // shifts later start seeds by one
        read.remove(1400); // shifts end seeds inward of it by one
        let original = seeds(&seq, &params);
        let shares = |a: &Minimizers, b: &Minimizers| {
            a.iter().any(|x| {
                b.iter()
                    .any(|y| x.hash == y.hash && x.offset.abs_diff(y.offset) <= 1)
            })
        };
        let fwd = seeds(&read, &params);
        assert!(shares(&original.start, &fwd.start) && shares(&original.end, &fwd.end));
        let rc = seeds(&reverse_complement(&read), &params);
        assert!(shares(&original.start, &rc.end) && shares(&original.end, &rc.start));
    }

    /// Independent per-k-mer reference for `kmer_hashes`.
    fn kmer_hash_oracle(kmer: &[u8]) -> u64 {
        let code = |b: u8| match b.to_ascii_uppercase() {
            b'A' => Some(0u64),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };
        let Some(fwd) = kmer
            .iter()
            .try_fold(0u64, |acc, &b| Some((acc << 2) | code(b)?))
        else {
            return INVALID;
        };
        let rc = kmer
            .iter()
            .rev()
            .fold(0u64, |acc, &b| (acc << 2) | (3 - code(b).unwrap()));
        let homopolymer = kmer.iter().all(|b| b.eq_ignore_ascii_case(&kmer[0]));
        if homopolymer {
            INVALID
        } else {
            mix(fwd.min(rc)).min(INVALID - 1)
        }
    }

    #[test]
    fn rolling_hashes_agree_with_oracle() {
        let mut rng = StdRng::seed_from_u64(5);
        let alphabet = b"ACGTACGTacgtN";
        for k in [2, 3, 17, 31, 32] {
            let mut seq: Vec<u8> = (0..300)
                .map(|_| alphabet[rng.gen_range(0..alphabet.len())])
                .collect();
            for (i, base) in b"ACGT".iter().enumerate() {
                // Include each homopolymer.
                let at = 40 + 60 * i;
                seq[at..at + 32].fill(*base);
            }
            let mut hashes = Vec::new();
            kmer_hashes(&seq, k, &mut hashes);
            assert_eq!(hashes.len(), seq.len() - k + 1);
            for (i, &hash) in hashes.iter().enumerate() {
                assert_eq!(hash, kmer_hash_oracle(&seq[i..i + k]), "k={k} pos={i}");
            }
        }
    }

    #[test]
    fn valid_kmer_never_hashes_to_sentinel() {
        // This 32-mer's mixed hash is u64::MAX before remapping.
        let kmer = b"TATTGCGGAACAGGTTTTGGCGGTGGTCTAAA";
        let mut hashes = Vec::new();
        kmer_hashes(kmer, 32, &mut hashes);
        assert_eq!(hashes, [INVALID - 1]);
        let params = OntMinimizerParams::new(32, 32, 1).unwrap();
        assert_eq!(seeds(kmer, &params).start.len(), 1);
    }

    #[test]
    fn seed_offsets_track_endpoint_displacement() {
        // Trimming the start leaves the end seeds unchanged, and seeds that
        // survive at the start move 3 bases closer to the read start.
        // (Window boundaries move too, so not every start seed survives.)
        let mut rng = StdRng::seed_from_u64(6);
        let params = OntMinimizerParams::default();
        let seq = random_seq(2000, &mut rng);
        let (full, trimmed) = (seeds(&seq, &params), seeds(&seq[3..], &params));
        assert_eq!(full.end, trimmed.end);
        let mut survivors = 0;
        for b in &trimmed.start {
            if let Some(a) = full.start.iter().find(|a| a.hash == b.hash) {
                assert_eq!(a.offset, b.offset + 3);
                survivors += 1;
            }
        }
        assert!(survivors > 0);
    }
}
