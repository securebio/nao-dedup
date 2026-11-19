//! Deciding whether two ONT reads (already in the same orientation) are
//! duplicates.
//!
//! The Illumina comparator tries small offsets, counts the offset as error,
//! and compares by Hamming distance. ONT errors are indel-heavy, and
//! duplicates are defined by shared fragment endpoints, so the check has two
//! parts:
//!
//! 1. **Endpoints.** At each end, the first `END_CHECK_LEN` bases are aligned
//!    with the outer boundary anchored. Endpoint displacement up to
//!    `max_end_offset` is free. The inner end of the window is free, so an
//!    indel inside the window doesn't push the alignment off its end. At most
//!    `MAX_END_ERRORS` edits are allowed, and a strictly cheaper alignment at
//!    a larger displacement rejects the pair (see `end_window_matches`). An
//!    indel near an end is tolerated, but a different fragment endpoint is
//!    not.
//! 2. **Whole read.** The edit distance between the reads must be at most
//!    a quality-derived lower-tail cutoff for a protected strain divergence,
//!    capped at `MAX_ERROR_FRAC` of the longer read.
//!
//! After a cheap length check against the cap, the endpoint check rejects
//! most non-duplicates before calculating the pair budget.

use smallvec::SmallVec;
use std::ops::Range;

// Fixed engineering defaults for ONT matching.
/// Bases inward from each boundary used for the endpoint check.
const END_CHECK_LEN: usize = 32;
/// Edits allowed within each endpoint check window.
const MAX_END_ERRORS: usize = 3;
/// Upper cap on edit distance as a fraction of the longer read.
const MAX_ERROR_FRAC: f64 = 0.12;
/// Multiplier on usable-quality expected errors.
const QUALITY_CALIBRATION: f64 = 1.2;
/// Additive errors per base with usable qualities.
const QUALITY_ERROR_FLOOR: f64 = 0.002;
/// Extra relative uncertainty in pair distance from sequencing noise.
const PAIR_DISTANCE_CV: f64 = 0.25;
/// Per-comparison lower-tail probability for the pair budget.
const PAIR_TAIL_PROBABILITY: f64 = 0.05;
/// Fractional indel-drift allowance for seed offsets.
const SEED_DRIFT_FRAC: f64 = 0.04;

/// Parameters for ONT duplicate matching.
///
/// Quality scores estimate sequencing noise. Reads merge only when their
/// distance is below a fixed 5% lower-tail cutoff for strains
/// `protect_divergence` apart, subject to a whole-read cap.
#[derive(Debug, Clone)]
pub struct OntDedupParams {
    /// Endpoint displacement (bases) allowed at each end at no cost.
    pub max_end_offset: usize,
    /// Strain divergence against which to protect (0 < value <= 1).
    pub protect_divergence: f64,
    /// Expected error probability for empty or exactly constant qualities.
    /// Default 0.01, unscaled by quality calibration.
    pub assumed_error_rate: f64,
    /// Whether a read may match the reverse complement of another.
    pub rc_tolerant: bool,
}

impl Default for OntDedupParams {
    fn default() -> Self {
        Self {
            max_end_offset: 2,
            protect_divergence: 0.015,
            assumed_error_rate: 0.01,
            rc_tolerant: true,
        }
    }
}

impl OntDedupParams {
    /// Validate parameters before processing reads.
    pub fn validate(&self) -> Result<(), String> {
        if self.max_end_offset > END_CHECK_LEN {
            return Err(format!(
                "max_end_offset must be at most {END_CHECK_LEN}, got {}",
                self.max_end_offset
            ));
        }
        for (name, value, positive) in [
            ("protect_divergence", self.protect_divergence, true),
            ("assumed_error_rate", self.assumed_error_rate, false),
        ] {
            if !(0.0..=1.0).contains(&value) || (positive && value == 0.0) {
                let lower = if positive { "greater than" } else { "at least" };
                return Err(format!(
                    "{name} must be {lower} 0 and at most 1, got {value}"
                ));
            }
        }
        Ok(())
    }
}

/// Whether a seed at `a` bases from its read end in one read and `b` bases
/// in another could belong to duplicates. Between duplicates a seed's offset
/// differs by the endpoint displacement plus the net indels between the seed
/// and the read end. Within `END_CHECK_LEN` of the end, and assuming the
/// shared seed is the same (homologous) k-mer in both reads, this bound is
/// exact: the endpoint check allows at most `max_end_offset` +
/// `MAX_END_ERRORS`. In repetitive sequence a shared seed may be a different
/// copy of the k-mer, and the filter can then skip a duplicate.
/// Further in, indels beyond the end window are only bounded by the
/// independent `SEED_DRIFT_FRAC`, so the allowance grows with distance; this part is
/// a heuristic, not a guarantee, and a rare dense indel cluster could hide a
/// duplicate.
pub(super) fn offsets_compatible(a: u16, b: u16, params: &OntDedupParams) -> bool {
    let drift = (SEED_DRIFT_FRAC * f64::from(a.max(b))).ceil() as usize;
    usize::from(a.abs_diff(b)) <= params.max_end_offset + MAX_END_ERRORS + drift
}

/// The offsets `b` compatible with `a` under `offsets_compatible`, which
/// form one contiguous range: below `a` the allowance is fixed, and above
/// `a` it grows by at most one per base.
pub(super) fn compatible_offset_range(a: u16, params: &OntDedupParams) -> (u16, u16) {
    let base = params.max_end_offset + MAX_END_ERRORS;
    let drift_at_a = (SEED_DRIFT_FRAC * f64::from(a)).ceil() as usize;
    let lo = usize::from(a).saturating_sub(base + drift_at_a) as u16;
    // Just above the largest b with b - a <= base + frac * b + 1.
    let estimate = (f64::from(a) + base as f64 + 1.0) / (1.0 - SEED_DRIFT_FRAC);
    let mut hi = estimate.min(f64::from(u16::MAX)) as u16;
    while hi > a && !offsets_compatible(a, hi, params) {
        hi -= 1;
    }
    (lo, hi)
}

/// Whether every byte is a valid Phred+33 quality (33..=126).
/// Empty qualities are valid and use the assumed error rate.
pub fn validate_quality(qual: &[u8]) -> bool {
    // The non-short-circuit fold lets the compiler vectorize byte validation.
    !qual
        .iter()
        .fold(false, |invalid, &q| invalid | !(33..=126).contains(&q))
}

pub(crate) struct QualitySummary {
    pub errors: f32,
    pub phred_sum: u64,
}

/// Summarize qualities in one pass over the read. Histogram reduction makes
/// expected errors independent of base order, including reverse complements.
pub(crate) fn quality_summary(len: usize, qual: &[u8]) -> QualitySummary {
    use std::sync::OnceLock;
    static PHRED: OnceLock<[f64; 94]> = OnceLock::new();
    assert!(
        qual.is_empty() || qual.len() == len,
        "nonempty quality length must match sequence length"
    );
    let mut counts = [0u32; 256];
    for &q in qual {
        counts[usize::from(q)] += 1;
    }
    assert!(
        counts[..33].iter().chain(&counts[127..]).all(|&n| n == 0),
        "quality bytes must be valid Phred+33 (33..=126)"
    );
    let table = PHRED.get_or_init(|| std::array::from_fn(|q| 10_f64.powf(-(q as f64) / 10.0)));
    let (mut errors, mut phred_sum, mut distinct) = (0.0, 0u64, 0);
    for (q, &count) in counts[33..=126].iter().enumerate() {
        errors += table[q] * f64::from(count);
        phred_sum += q as u64 * u64::from(count);
        distinct += usize::from(count != 0);
    }
    QualitySummary {
        errors: if distinct > 1 { errors as f32 } else { -1.0 },
        phred_sum,
    }
}

#[cfg(test)]
pub(crate) fn expected_errors(len: usize, qual: &[u8]) -> f32 {
    quality_summary(len, qual).errors
}

/// Usable-quality noise; unavailable summaries have no quality-derived estimate.
fn read_noise(len: usize, errors: f32) -> Option<f64> {
    (errors >= 0.0)
        .then(|| QUALITY_CALIBRATION * f64::from(errors) + QUALITY_ERROR_FLOOR * len as f64)
}

/// Extra spread models noise misprediction, while divergence adds Poisson variance.
fn noise_variance(mu: f64, noise: f64, cv: f64) -> f64 {
    mu + (cv * noise).powi(2)
}

/// Upper cap; safe for a cheap precheck before calculating the pair budget.
pub(crate) fn whole_read_cap(a_len: usize, b_len: usize) -> usize {
    let longer = a_len.max(b_len);
    (MAX_ERROR_FRAC * longer as f64).floor() as usize
}

/// None represents T=-1. Unavailable summaries get an assumed expected-error
/// count, not a cap-sized acceptance budget. Arithmetic after summary rounding
/// is f64 and symmetric, including when only one read has usable qualities.
/// Usable qualities contribute c * expected_errors + f * length.
/// Unavailable qualities contribute assumed_error_rate * length, without
/// quality terms.
pub(crate) fn pair_budget(
    a_len: usize,
    b_len: usize,
    a_errors: f32,
    b_errors: f32,
    params: &OntDedupParams,
) -> Option<usize> {
    if a_len == 0 || b_len == 0 {
        return None;
    }
    let errors = |len, e| read_noise(len, e).unwrap_or(params.assumed_error_rate * len as f64);
    let noise = errors(a_len, a_errors) + errors(b_len, b_errors);
    let mu = noise + params.protect_divergence * a_len.min(b_len) as f64;
    // The extra spread models misprediction of the noise, so it scales with
    // the noise only: variance = mu + (PAIR_DISTANCE_CV * noise)^2. Applied to
    // all of mu, it would keep the cutoff a fixed fraction below the mean at
    // every length, and longer reads would gain no power.
    lower_tail_budget(
        mu,
        noise_variance(mu, noise, PAIR_DISTANCE_CV),
        whole_read_cap(a_len, b_len),
        PAIR_TAIL_PROBABILITY,
    )
}

/// Largest k <= cap with F(k) < tail_probability for NB(mean, variance), or Poisson
/// when variance equals the mean.
/// Sum relative masses starting at one, keeping log(P(0)) separately. Rescale
/// before masses grow too large. Never initialize a recurrence with exp(logP0),
/// which can underflow and silently turn an entire CDF into zero.
fn lower_tail_budget(mu: f64, variance: f64, cap: usize, tail_probability: f64) -> Option<usize> {
    if !mu.is_finite() || mu < 0.0 || !variance.is_finite() || variance < mu {
        return None;
    }
    if mu == 0.0 {
        return None;
    }
    let h2 = (variance - mu) / mu / mu;
    let scaled_mean = mu * h2;
    if !h2.is_finite() || !scaled_mean.is_finite() {
        return None;
    }
    let mut log_scale = if h2 == 0.0 || scaled_mean == 0.0 {
        -mu
    } else {
        -scaled_mean.ln_1p() / h2
    };
    if !log_scale.is_finite() {
        return None;
    }
    let mut threshold = (tail_probability.ln() - log_scale).exp();
    let (mut mass, mut sum) = (1.0, 1.0);
    if sum >= threshold {
        return None;
    }
    let step = mu / (1.0 + scaled_mean);
    for k in 1..=cap {
        let ratio = step * (1.0 + (k - 1) as f64 * h2) / k as f64;
        mass *= ratio;
        sum += mass;
        if !sum.is_finite() {
            return None;
        }
        if sum >= threshold {
            return Some(k - 1);
        }
        if sum > 1e100 {
            log_scale += sum.ln();
            mass /= sum;
            sum = 1.0;
            threshold = (tail_probability.ln() - log_scale).exp();
        }
    }
    Some(cap)
}

/// Edit distance is at least the length difference, so this cheap check can
/// reject before any alignment.
pub(crate) fn lengths_compatible(a_len: usize, b_len: usize, budget: usize) -> bool {
    a_len.abs_diff(b_len) <= budget
}

/// Ranges of the start and end windows of a sequence of length `len`.
pub(crate) fn end_windows(len: usize) -> (Range<usize>, Range<usize>) {
    let n = END_CHECK_LEN.min(len);
    (0..n, len - n..len)
}

/// Whether `a`'s ends match `b`'s, given `b`'s `end_windows`.
pub(crate) fn ends_match(a: &[u8], b_start: &[u8], b_end: &[u8], params: &OntDedupParams) -> bool {
    let (a_start, a_end) = end_windows(a.len());
    if !end_window_matches(&a[a_start], b_start, params) {
        return false;
    }
    // Align the 3' ends from the boundary inward too.
    let reversed = |s: &[u8]| -> SmallVec<[u8; 256]> { s.iter().rev().copied().collect() };
    end_window_matches(&reversed(&a[a_end]), &reversed(b_end), params)
}

#[cfg(test)]
pub(crate) fn whole_read_matches(a: &[u8], b: &[u8], budget: usize) -> bool {
    editdistancek::edit_distance_bounded(a, b, budget).is_some()
}

/// Whether two end windows, read from the boundary inward, agree.
///
/// Edit distance alone can't tell a displaced endpoint from indels at the
/// boundary: d extra leading bases cost the same as d insertions. Allowing
/// errors would then quietly widen the endpoint tolerance by up to
/// `MAX_END_ERRORS`. So two alignments are compared: the cheapest with free
/// displacement up to `max_end_offset` must fit the error allowance, and must
/// cost no more than the cheapest alignment allowing displacement up to
/// `max_end_offset + MAX_END_ERRORS`. A strictly cheaper alignment at larger
/// displacement means the endpoints differ.
///
/// Low-complexity ends (homopolymers, short tandem repeats) are inherently
/// ambiguous: shifting along a repeat can change which alignment is
/// cheapest, so such ends may be rejected despite matching, or accepted
/// despite a displacement within the repeat.
fn end_window_matches(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
    let (m, max_errors) = (params.max_end_offset, MAX_END_ERRORS);
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
    use crate::ont::test_util::{random_seq, substitute};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn params(max_end_offset: usize) -> OntDedupParams {
        OntDedupParams {
            max_end_offset,
            ..Default::default()
        }
    }

    fn ends_match_full(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
        let (start, end) = end_windows(b.len());
        ends_match(a, &b[start], &b[end], params)
    }

    #[test]
    fn endpoint_displacement_within_offset_is_free() {
        let mut rng = StdRng::seed_from_u64(2);
        let seq = random_seq(500, &mut rng);
        // Shift the start by 3 and the end by 2 in each direction.
        for (a, b) in [(&seq[3..498], &seq[..]), (&seq[..], &seq[3..498])] {
            assert!(ends_match_full(a, b, &params(3)));
            assert!(!ends_match_full(a, b, &params(2)));
        }
    }

    #[test]
    fn different_endpoints_are_rejected_despite_whole_read_budget() {
        let mut rng = StdRng::seed_from_u64(3);
        let seq = random_seq(2000, &mut rng);
        // 20 bp shorter at the start: within a 5% whole-read budget (100
        // edits), but a different fragment endpoint.
        let shifted = &seq[20..];
        assert!(whole_read_matches(shifted, &seq, 100));
        assert!(!ends_match_full(shifted, &seq, &params(3)));
        // Same at the 3' end.
        let shifted = &seq[..1980];
        assert!(!ends_match_full(shifted, &seq, &params(3)));
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
                assert!(!ends_match_full(a, b, &params(3)), "d={d}");
                assert!(ends_match_full(a, b, &params(d)), "d={d}");
            }
        }
    }

    #[test]
    fn repeat_period_displacement_is_rejected() {
        // A period-6 repeat at the start, truncated by one period, with one
        // substitution: cheaper to explain by displacing the endpoint 6 bases.
        let mut rng = StdRng::seed_from_u64(11);
        let mut seq = b"ACGTGC".repeat(16);
        seq.extend(random_seq(1000, &mut rng));
        let mut truncated = seq[6..].to_vec();
        truncated[0] = b'C';
        assert!(!ends_match_full(&truncated, &seq, &params(3)));
    }

    #[test]
    fn indels_near_ends_are_tolerated() {
        let mut rng = StdRng::seed_from_u64(4);
        let seq = random_seq(1000, &mut rng);
        let mut read = seq.clone();
        read.remove(10); // deletion near the start
        read.insert(20, b'G'); // insertion inside the start window
        read.remove(980); // deletion near the end
        substitute(&mut read, 500);
        assert!(ends_match_full(&read, &seq, &params(0)));
        assert!(ends_match_full(&seq, &read, &params(0)));
    }

    #[test]
    fn default_endpoint_rule() {
        let d = OntDedupParams::default();
        assert_eq!(
            (d.max_end_offset, END_CHECK_LEN, MAX_END_ERRORS),
            (2, 32, 3)
        );
        assert_eq!((MAX_ERROR_FRAC, d.rc_tolerant), (0.12, true));
        // With the defaults: 2 bases of displacement pass, 3 don't (even
        // though 3 edits are allowed), and 3 substitutions in the 32-base
        // window pass while 4 don't.
        let mut rng = StdRng::seed_from_u64(13);
        let seq = random_seq(1000, &mut rng);
        assert!(ends_match_full(&seq[2..], &seq, &d));
        assert!(!ends_match_full(&seq[3..], &seq, &d));
        let mut read = seq.clone();
        for pos in [4, 12, 20] {
            substitute(&mut read, pos);
        }
        assert!(ends_match_full(&read, &seq, &d));
        substitute(&mut read, 28);
        assert!(!ends_match_full(&read, &seq, &d));
    }

    #[test]
    fn errors_at_the_boundary_are_tolerated_with_defaults() {
        // Errors within the first max_end_offset + max_end_errors bases are
        // where the within/wider comparison could misfire.
        let mut rng = StdRng::seed_from_u64(12);
        let seq = random_seq(1500, &mut rng);
        let defaults = OntDedupParams::default();
        let edits: [fn(&mut Vec<u8>); 4] = [
            |r| substitute(r, 0),
            |r| substitute(r, 2),
            |r| r.insert(1, b'T'),
            |r| {
                r.remove(3);
                substitute(r, 6);
            },
        ];
        for edit in edits {
            let mut read = seq.clone();
            edit(&mut read);
            assert!(ends_match_full(&read, &seq, &defaults));
            assert!(ends_match_full(&seq, &read, &defaults));
        }
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
        assert!(ends_match_full(&read, &seq, &params(0)));
        assert!(!whole_read_matches(&read, &seq, 50));
        assert!(whole_read_matches(&read, &seq, 60));
    }

    #[test]
    fn short_reads_use_whole_read_as_end_window() {
        let mut rng = StdRng::seed_from_u64(7);
        let seq = random_seq(24, &mut rng);
        assert!(ends_match_full(&seq, &seq, &params(0)));
        assert!(!ends_match_full(&seq[5..], &seq, &params(2)));
    }

    #[test]
    fn unrelated_reads_are_rejected() {
        let mut rng = StdRng::seed_from_u64(8);
        for _ in 0..50 {
            let (a, b) = (random_seq(800, &mut rng), random_seq(800, &mut rng));
            assert!(!ends_match_full(&a, &b, &params(3)));
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

#[cfg(test)]
mod quality_tests {
    use super::*;
    fn cutoff(mu: f64, cv: f64, cap: usize) -> Option<usize> {
        lower_tail_budget(mu, noise_variance(mu, mu, cv), cap, 0.05)
    }

    /// The budget for a pair with this noise and divergence term, cv 0.25.
    fn expected(noise: f64, divergence: f64, cap: usize) -> Option<usize> {
        let mu = noise + divergence;
        cutoff(mu, 0.25 * noise / mu, cap)
    }

    #[test]
    fn calibration_and_floor_apply_only_to_usable_qualities() {
        let p = OntDedupParams {
            protect_divergence: 0.01,
            ..Default::default()
        };
        assert_eq!(read_noise(1000, -1.0), None);
        assert_eq!(read_noise(1000, 10.0), Some(14.0));
        let budget = |ea, eb| pair_budget(1000, 1200, ea, eb, &p);
        assert_eq!(budget(-1.0, -1.0), expected(22.0, 10.0, 144));
        assert_eq!(budget(10.0, -1.0), expected(14.0 + 12.0, 10.0, 144));
        assert_eq!(
            budget(-1.0, 12.0),
            expected(10.0 + 1.2 * 12.0 + 2.4, 10.0, 144)
        );
        assert_eq!(budget(10.0, 12.0), expected(1.2 * 22.0 + 4.4, 10.0, 144));
    }

    #[test]
    fn tail_probability_controls_strictness() {
        let budget =
            |tail| lower_tail_budget(100.0, noise_variance(100.0, 80.0, 0.25), 384, tail).unwrap();
        assert!(budget(0.01) < budget(0.05));
        assert!(budget(0.1) > budget(0.05));
    }

    #[test]
    fn validates_statistical_parameters() {
        type Setter = fn(&mut OntDedupParams, f64);
        // (setter, must be positive, upper bound)
        let fields: [(Setter, bool, f64); 2] = [
            (|p, v| p.protect_divergence = v, true, 1.0),
            (|p, v| p.assumed_error_rate = v, false, 1.0),
        ];
        for (set, positive, max) in fields {
            for value in [
                f64::NAN,
                f64::INFINITY,
                f64::NEG_INFINITY,
                -0.1,
                0.0,
                0.04,
                1.0,
                1.1,
                100.0,
                100.1,
                f64::MAX,
            ] {
                let mut p = OntDedupParams::default();
                set(&mut p, value);
                let valid =
                    value.is_finite() && value >= 0.0 && (!positive || value > 0.0) && value <= max;
                assert_eq!(p.validate().is_ok(), valid, "{p:?}");
            }
        }
        for (offset, valid) in [(0, true), (END_CHECK_LEN, true), (END_CHECK_LEN + 1, false)] {
            let p = OntDedupParams {
                max_end_offset: offset,
                ..Default::default()
            };
            assert_eq!(p.validate().is_ok(), valid, "{p:?}");
        }
    }

    #[test]
    fn phred_sums_availability_and_no_floor() {
        assert!((f64::from(expected_errors(3, b"!+5")) - 1.11).abs() < 1e-7);
        assert_eq!(expected_errors(10, b""), -1.0);
        assert_eq!(expected_errors(3, b"555"), -1.0);
        assert_eq!(expected_errors(0, b""), -1.0);
        assert!(expected_errors(2, b"}~") < 2e-9);
        let q = b"+5?I";
        assert_eq!(
            expected_errors(4, q),
            expected_errors(4, &q.iter().rev().copied().collect::<Vec<_>>())
        );
        for (len, q) in [(1, &b"55"[..]), (2, &b"  "[..]), (2, &b"5\x7f"[..])] {
            assert!(std::panic::catch_unwind(|| expected_errors(len, q)).is_err());
        }
    }

    #[test]
    fn discrete_goldens_and_strict_boundary() {
        assert_eq!(cutoff(3.0, 0.25, 1000), None);
        assert_eq!(cutoff(4.0, 0.25, 1000), Some(0));
        assert_eq!(cutoff(30.0, 0.25, 1000), Some(15));
        assert_eq!(cutoff(30.0, 0.25, 7), Some(7));
        assert_eq!(cutoff(-0.05_f64.ln(), 0.0, 100), None);
        assert_eq!(cutoff(-0.05_f64.ln() + 1e-10, 0.0, 100), Some(0));
        assert_eq!(cutoff(0.0, 0.0, 100), None);
        // Unlike a floor at mu0, the lower tail deliberately sacrifices recall.
        let p = OntDedupParams {
            protect_divergence: 0.01,
            ..Default::default()
        };
        // Expected noise is 1.2 * (10 + 10) + 2 * 0.002 * 1000 = 28 edits.
        assert!(pair_budget(1000, 1000, 10.0, 10.0, &p).unwrap() < 28);
    }

    #[test]
    fn large_means_and_small_dispersion_do_not_underflow_to_cap() {
        // Poisson(1000) lower 5% tail ends at 947; P(0) underflows.
        for cv in [0.0, 1e-12, 1e-6] {
            assert_eq!(cutoff(1000.0, cv, 2000), Some(947));
        }
        assert_eq!(cutoff(10000.0, 0.0, 20000), Some(9835));
        assert_eq!(cutoff(f64::INFINITY, 0.25, 1000), None);
        assert_eq!(cutoff(f64::NAN, 0.25, 1000), None);
    }

    #[test]
    fn budgets_are_symmetric_and_monotone() {
        let p = OntDedupParams {
            protect_divergence: 0.01,
            ..Default::default()
        };
        let (a, b) = (expected_errors(3, b"!+5"), expected_errors(4, b"+5?I"));
        assert_eq!(
            pair_budget(1000, 1200, a, b, &p),
            pair_budget(1200, 1000, b, a, &p)
        );
        assert_eq!(
            pair_budget(1000, 1200, -1.0, b, &p),
            pair_budget(1200, 1000, b, -1.0, &p)
        );
        // Unavailable qualities use the assumed rate independently of quality terms.
        assert_eq!(
            pair_budget(1000, 1200, -1.0, -1.0, &p),
            expected(22.0, 10.0, 144)
        );
        for cv in [0.0, 0.01, 0.25, 1.0] {
            let mut last = None;
            for i in 0..500 {
                let next = cutoff(i as f64 * 0.5, cv, 300);
                assert!(next >= last, "cv={cv} i={i}");
                last = next;
            }
        }
        type Set = fn(&mut OntDedupParams, &mut f32, &mut f32, usize);
        let setters: [Set; 3] = [
            |_, a, _, i| *a = i as f32,
            |_, _, b, i| *b = i as f32,
            |p, _, _, i| p.protect_divergence = i as f64 / 1000.0,
        ];
        for (field, set) in setters.into_iter().enumerate() {
            let mut last = None;
            for i in 1..=100 {
                let mut q = p.clone();
                let (mut ea, mut eb) = (10.0, 12.0);
                set(&mut q, &mut ea, &mut eb, i);
                let next = pair_budget(1000, 1200, ea, eb, &q);
                assert!(next >= last, "field={field} i={i}");
                last = next;
            }
        }
        // Divergence evidence uses the shorter read; don't add length difference.
        assert_eq!(
            pair_budget(1000, 1200, 10.0, 12.0, &p),
            expected(1.2 * 22.0 + 4.4, 0.01 * 1000.0, 144)
        );
    }

    #[test]
    fn lengths_use_exact_integer_budget() {
        assert_eq!(whole_read_cap(1000, 900), 120);
        assert_eq!(whole_read_cap(900, 999), 119);
        assert!(lengths_compatible(1000, 880, whole_read_cap(1000, 880)));
        assert!(!lengths_compatible(1000, 879, whole_read_cap(1000, 879)));
        assert!(lengths_compatible(1000, 1015, 15));
        assert!(!lengths_compatible(1000, 1016, 15));
        assert!(lengths_compatible(1015, 1000, 15));
        assert_eq!(
            pair_budget(0, 0, -1.0, -1.0, &OntDedupParams::default()),
            None
        );
    }
}
