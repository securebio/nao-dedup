//! Deciding whether two ONT reads (already in the same orientation) are
//! duplicates.
//!
//! The Illumina comparator tries small offsets, counts the offset as error,
//! and compares by Hamming distance. ONT errors are indel-heavy, and
//! duplicates are defined by shared fragment endpoints, so the check has two
//! parts:
//!
//! 1. **Endpoints.** At each end, the first `end_check_len` bases are aligned
//!    with the outer boundary anchored. Endpoint displacement up to
//!    `max_end_offset` is free. The inner end of the window is free, so an
//!    indel inside the window doesn't push the alignment off its end. At most
//!    `max_end_errors` edits are allowed, and a strictly cheaper alignment at
//!    a larger displacement rejects the pair (see `end_window_matches`). An
//!    indel near an end is tolerated, but a different fragment endpoint is
//!    not.
//! 2. **Whole read.** The edit distance between the reads must be at most
//!    `max_error_frac` of the longer read.
//!
//! The endpoint check runs first because it is cheap and rejects most
//! non-duplicates.

use smallvec::SmallVec;
use std::ops::Range;

/// Parameters for ONT duplicate matching.
///
/// The defaults suit reads of roughly 1% per-read error (typical of current
/// ONT chemistry and high-accuracy basecalling): two such reads differ by
/// about 2%, comfortably inside `max_error_frac`. They were chosen by sweeps
/// on simulated reads. Noisier reads need a larger `max_error_frac` and
/// `max_end_errors`; otherwise true duplicates are silently kept as unique.
#[derive(Debug, Clone)]
pub struct OntDedupParams {
    /// Endpoint displacement (bases) allowed at each end at no cost.
    pub max_end_offset: usize,
    /// Number of bases inward from each end used for the endpoint check.
    pub end_check_len: usize,
    /// Edits allowed within each endpoint check window.
    pub max_end_errors: usize,
    /// Whole-read edit distance allowed, as a fraction of the longer read.
    pub max_error_frac: f64,
    /// Whether a read may match the reverse complement of another.
    pub rc_tolerant: bool,
}

impl Default for OntDedupParams {
    fn default() -> Self {
        Self {
            max_end_offset: 2,
            end_check_len: 64,
            max_end_errors: 5,
            max_error_frac: 0.04,
            rc_tolerant: true,
        }
    }
}

/// Whether `a` and `b`, in the same orientation, are duplicates. Comparison
/// is case-sensitive; callers pass uppercase sequences.
#[cfg(test)]
pub(crate) fn is_duplicate(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
    let (b_start, b_end) = end_windows(b.len(), params);
    lengths_compatible(a.len(), b.len(), params)
        && ends_match(a, &b[b_start], &b[b_end], params)
        && whole_read_matches(a, b, params)
}

// The stages of `is_duplicate`, in order, for callers that store `b` packed
// and want to unpack only what each stage needs.

/// Edit distance is at least the length difference, so this cheap check can
/// reject before any alignment.
pub(crate) fn lengths_compatible(a_len: usize, b_len: usize, params: &OntDedupParams) -> bool {
    a_len.abs_diff(b_len) <= whole_read_budget(a_len, b_len, params)
}

/// Ranges of the start and end windows of a sequence of length `len`.
pub(crate) fn end_windows(len: usize, params: &OntDedupParams) -> (Range<usize>, Range<usize>) {
    let n = params.end_check_len.min(len);
    (0..n, len - n..len)
}

/// Whether `a`'s ends match `b`'s, given `b`'s `end_windows`.
pub(crate) fn ends_match(a: &[u8], b_start: &[u8], b_end: &[u8], params: &OntDedupParams) -> bool {
    let (a_start, a_end) = end_windows(a.len(), params);
    if !end_window_matches(&a[a_start], b_start, params) {
        return false;
    }
    // Align the 3' ends from the boundary inward too.
    let reversed = |s: &[u8]| -> SmallVec<[u8; 256]> { s.iter().rev().copied().collect() };
    end_window_matches(&reversed(&a[a_end]), &reversed(b_end), params)
}

pub(crate) fn whole_read_matches(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
    editdistancek::edit_distance_bounded(a, b, whole_read_budget(a.len(), b.len(), params))
        .is_some()
}

fn whole_read_budget(a_len: usize, b_len: usize, params: &OntDedupParams) -> usize {
    let longer = a_len.max(b_len);
    // Edit distance never exceeds the longer length, so cap there (this also
    // keeps out-of-range fractions from overflowing).
    ((params.max_error_frac * longer as f64).floor() as usize).min(longer)
}

/// Whether two end windows, read from the boundary inward, agree.
///
/// Edit distance alone can't tell a displaced endpoint from indels at the
/// boundary: d extra leading bases cost the same as d insertions. Allowing
/// errors would then quietly widen the endpoint tolerance by up to
/// `max_end_errors`. So two alignments are compared: the cheapest with free
/// displacement up to `max_end_offset` must fit the error allowance, and must
/// cost no more than the cheapest alignment allowing displacement up to
/// `max_end_offset + max_end_errors`. A strictly cheaper alignment at larger
/// displacement means the endpoints differ.
///
/// Low-complexity ends (homopolymers, short tandem repeats) are inherently
/// ambiguous: shifting along a repeat can change which alignment is
/// cheapest, so such ends may be rejected despite matching, or accepted
/// despite a displacement within the repeat.
fn end_window_matches(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
    let (m, max_errors) = (params.max_end_offset, params.max_end_errors);
    let Some(within) = anchored_end_cost(a, b, m, max_errors as u32) else {
        return false;
    };
    if within == 0 {
        return true;
    }
    // Look for a cheaper alignment with free displacement up to
    // m + max_errors: the range a within-budget alignment could have paid
    // for with errors. (Paid skips beyond that radius are still possible; in
    // tandem repeats they can find a cheaper alignment, which is part of the
    // documented ambiguity of low-complexity ends.)
    let wider = anchored_end_cost(a, b, m + max_errors, within).unwrap_or(within);
    within <= wider
}

/// Banded edit distance between two end windows, read from the anchored
/// (outer) boundary inward, or None if it exceeds `max_cost`. Up to
/// `free_skip` leading bases of either window may be skipped for free;
/// further skipping costs one per base. The alignment may end anywhere once
/// it has consumed all of either window.
fn anchored_end_cost(a: &[u8], b: &[u8], free_skip: usize, max_cost: u32) -> Option<u32> {
    const INF: u32 = u32::MAX / 2;
    // Any path costing <= max_cost stays within this many diagonals of the
    // main one: at most `free_skip` free skips plus `max_cost` paid edits.
    let band = free_skip + max_cost as usize;
    let nb = b.len();

    // Two DP rows: cost of aligning the first i bases of `a` with the first
    // j bases of `b`. Only cells in the band [i - band, i + band] are
    // computed; bands only move right, so cells right of the band are still
    // INF from initialization.
    let mut prev: SmallVec<[u32; 256]> = (0..=nb)
        .map(|j| {
            if j <= band {
                j.saturating_sub(free_skip) as u32
            } else {
                INF
            }
        })
        .collect();
    let mut cur: SmallVec<[u32; 256]> = SmallVec::from_elem(INF, nb + 1);
    let mut best = prev[nb]; // alignments that consumed all of `b`
    let (mut lo, mut hi) = (0, band.min(nb));
    let result = |cost: u32| (cost <= max_cost).then_some(cost);

    for (i, &ai) in (1usize..).zip(a) {
        lo = i.saturating_sub(band);
        hi = (i + band).min(nb);
        if lo > hi {
            // `a` has run past the band's reach into `b`: only alignments
            // that already consumed all of `b` (counted in `best`) remain.
            return result(best);
        }
        if lo == 0 {
            cur[0] = i.saturating_sub(free_skip) as u32;
        } else {
            cur[lo - 1] = INF; // may hold a stale value from two rows ago
        }
        let mut row_min = if lo == 0 { cur[0] } else { INF };
        for j in lo.max(1)..=hi {
            let diag = prev[j - 1] + u32::from(ai != b[j - 1]);
            let cost = diag.min(prev[j] + 1).min(cur[j - 1] + 1);
            cur[j] = cost;
            row_min = row_min.min(cost);
        }
        if row_min > max_cost {
            // Costs never decrease down the rows.
            return result(best);
        }
        if hi == nb {
            best = best.min(cur[nb]);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    // All of `a` consumed: the alignment may end anywhere in `b`.
    let last_row = prev[lo..=hi].iter().copied().min().unwrap_or(INF);
    result(best.min(last_row))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn random_seq(len: usize, rng: &mut StdRng) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
    }

    fn params(max_end_offset: usize, max_end_errors: usize) -> OntDedupParams {
        OntDedupParams {
            max_end_offset,
            end_check_len: 64,
            max_end_errors,
            max_error_frac: 0.05,
            rc_tolerant: true,
        }
    }

    fn ends_match_full(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
        let (start, end) = end_windows(b.len(), params);
        ends_match(a, &b[start], &b[end], params)
    }

    fn substitute(seq: &mut [u8], pos: usize) {
        seq[pos] = if seq[pos] == b'A' { b'C' } else { b'A' };
    }

    #[test]
    fn identical_reads_match() {
        let mut rng = StdRng::seed_from_u64(1);
        let seq = random_seq(500, &mut rng);
        assert!(is_duplicate(&seq, &seq, &params(0, 0)));
    }

    #[test]
    fn endpoint_displacement_within_offset_is_free() {
        let mut rng = StdRng::seed_from_u64(2);
        let seq = random_seq(500, &mut rng);
        // Shift the start by 3 and the end by 2 in each direction.
        for (a, b) in [(&seq[3..498], &seq[..]), (&seq[..], &seq[3..498])] {
            assert!(ends_match_full(a, b, &params(3, 0)));
            assert!(!ends_match_full(a, b, &params(2, 0)));
            // The error allowance doesn't widen the endpoint tolerance.
            assert!(!ends_match_full(a, b, &params(2, 5)));
        }
    }

    #[test]
    fn different_endpoints_are_rejected_despite_whole_read_budget() {
        let mut rng = StdRng::seed_from_u64(3);
        let seq = random_seq(2000, &mut rng);
        // 20 bp shorter at the start: within a 5% whole-read budget (100
        // edits), but a different fragment endpoint.
        let shifted = &seq[20..];
        assert!(whole_read_matches(shifted, &seq, &params(3, 6)));
        assert!(!is_duplicate(shifted, &seq, &params(3, 6)));
        // Same at the 3' end.
        let shifted = &seq[..1980];
        assert!(!is_duplicate(shifted, &seq, &params(3, 6)));
    }

    #[test]
    fn displacement_beyond_offset_is_rejected_at_either_end() {
        let mut rng = StdRng::seed_from_u64(10);
        let seq = random_seq(1000, &mut rng);
        for d in [4, 9] {
            let trimmed_start = &seq[d..];
            let trimmed_end = &seq[..seq.len() - d];
            for (a, b) in [
                (trimmed_start, &seq[..]),
                (&seq[..], trimmed_start),
                (trimmed_end, &seq[..]),
                (&seq[..], trimmed_end),
            ] {
                assert!(!ends_match_full(a, b, &params(3, 6)), "d={d}");
                assert!(ends_match_full(a, b, &params(d, 0)), "d={d}");
            }
        }
    }

    #[test]
    fn repeat_period_displacement_is_rejected() {
        // A period-8 repeat at the start, truncated by one period, with one
        // substitution: cheaper to explain by displacing the endpoint 8 bases.
        let mut rng = StdRng::seed_from_u64(11);
        let mut seq = b"ACGTGCTA".repeat(16);
        seq.extend(random_seq(1000, &mut rng));
        let mut truncated = seq[8..].to_vec();
        truncated[0] = b'C';
        assert!(!ends_match_full(&truncated, &seq, &params(3, 5)));
    }

    #[test]
    fn indels_near_ends_are_tolerated() {
        let mut rng = StdRng::seed_from_u64(4);
        let seq = random_seq(1000, &mut rng);
        let mut read = seq.clone();
        read.remove(10); // deletion near the start
        read.insert(40, b'G'); // insertion inside the start window
        read.remove(980); // deletion near the end
        substitute(&mut read, 500);
        assert!(is_duplicate(&read, &seq, &params(0, 3)));
        assert!(is_duplicate(&seq, &read, &params(0, 3)));
    }

    #[test]
    fn errors_at_the_boundary_are_tolerated_with_defaults() {
        // Errors within the first max_end_offset + max_end_errors bases are
        // where the within/wider comparison could misfire.
        let mut rng = StdRng::seed_from_u64(12);
        let seq = random_seq(1500, &mut rng);
        let defaults = OntDedupParams::default();
        let edits: [&dyn Fn(&mut Vec<u8>); 4] = [
            &|r| substitute(r, 0),
            &|r| substitute(r, 2),
            &|r| r.insert(1, b'T'),
            &|r| {
                r.remove(3);
                substitute(r, 6);
            },
        ];
        for edit in edits {
            let mut read = seq.clone();
            edit(&mut read);
            assert!(is_duplicate(&read, &seq, &defaults));
            assert!(is_duplicate(&seq, &read, &defaults));
        }
    }

    #[test]
    fn end_errors_over_allowance_are_rejected() {
        let mut rng = StdRng::seed_from_u64(5);
        let seq = random_seq(1000, &mut rng);
        let mut read = seq.clone();
        for pos in [5, 15, 25, 35] {
            substitute(&mut read, pos);
        }
        assert!(!ends_match_full(&read, &seq, &params(3, 3)));
        assert!(ends_match_full(&read, &seq, &params(3, 4)));
    }

    #[test]
    fn whole_read_budget_is_enforced() {
        let mut rng = StdRng::seed_from_u64(6);
        let seq = random_seq(1000, &mut rng);
        let mut read = seq.clone();
        // 60 substitutions in the middle: ends are clean, budget is 50.
        for pos in (200..800).step_by(10) {
            substitute(&mut read, pos);
        }
        assert!(ends_match_full(&read, &seq, &params(0, 0)));
        assert!(!is_duplicate(&read, &seq, &params(0, 0)));
        let loose = OntDedupParams {
            max_error_frac: 0.06,
            ..params(0, 0)
        };
        assert!(is_duplicate(&read, &seq, &loose));
    }

    #[test]
    fn short_reads_use_whole_read_as_end_window() {
        let mut rng = StdRng::seed_from_u64(7);
        let seq = random_seq(40, &mut rng);
        assert!(ends_match_full(&seq, &seq, &params(0, 0)));
        assert!(!ends_match_full(&seq[5..], &seq, &params(2, 0)));
        assert!(!is_duplicate(&seq, b"", &params(3, 6)));
    }

    #[test]
    fn unrelated_reads_are_rejected() {
        let mut rng = StdRng::seed_from_u64(8);
        for _ in 0..50 {
            let (a, b) = (random_seq(800, &mut rng), random_seq(800, &mut rng));
            assert!(!is_duplicate(&a, &b, &params(3, 6)));
        }
    }

    /// Unbanded reference implementation of `anchored_end_cost`.
    fn anchored_end_oracle(a: &[u8], b: &[u8], free_skip: usize) -> usize {
        let (na, nb) = (a.len(), b.len());
        let mut d = vec![vec![0usize; nb + 1]; na + 1];
        for (j, cell) in d[0].iter_mut().enumerate() {
            *cell = j.saturating_sub(free_skip);
        }
        for i in 1..=na {
            d[i][0] = i.saturating_sub(free_skip);
            for j in 1..=nb {
                let diag = d[i - 1][j - 1] + usize::from(a[i - 1] != b[j - 1]);
                d[i][j] = diag.min(d[i - 1][j] + 1).min(d[i][j - 1] + 1);
            }
        }
        let last_row = d[na].iter().min().copied().unwrap();
        let last_col = d.iter().map(|row| row[nb]).min().unwrap();
        last_row.min(last_col)
    }

    #[test]
    fn banded_end_cost_agrees_with_full_dp() {
        let mut rng = StdRng::seed_from_u64(9);
        for trial in 0..3000 {
            let len = rng.gen_range(0..80);
            let base = random_seq(len, &mut rng);
            let mut other = base.clone();
            // A mix of near-duplicates and heavily edited windows.
            for _ in 0..rng.gen_range(0..12) {
                if other.is_empty() {
                    break;
                }
                let pos = rng.gen_range(0..other.len());
                match rng.gen_range(0..4) {
                    0 => {
                        other.remove(pos);
                    }
                    1 => other.insert(pos, b"ACGT"[rng.gen_range(0..4)]),
                    2 => substitute(&mut other, pos),
                    _ => {
                        other.drain(..rng.gen_range(0..=pos.min(8)));
                    }
                }
            }
            let (skip, max_cost) = (rng.gen_range(0..12), rng.gen_range(0..8));
            let expected = anchored_end_oracle(&base, &other, skip) as u32;
            assert_eq!(
                anchored_end_cost(&base, &other, skip, max_cost),
                (expected <= max_cost).then_some(expected),
                "trial {trial}: skip={skip} max_cost={max_cost}"
            );
        }
    }
}
