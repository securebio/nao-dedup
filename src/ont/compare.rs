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
//!    a quality-derived lower-tail cutoff for a protected strain divergence,
//!    capped at `max_error_frac` of the longer read.
//!
//! After a cheap length check against the cap, the endpoint check rejects
//! most non-duplicates before calculating the pair budget.

use smallvec::SmallVec;
use std::ops::Range;

/// Parameters for ONT duplicate matching.
///
/// The whole-read edit budget is set per pair (see `pair_budget`): quality
/// scores estimate the pair's sequencing noise, and reads merge only if their
/// distance is unlikely (< 5%) for reads of strains `protect_divergence`
/// apart. `max_error_frac` only caps the budget. `pair_distance_cv` is a
/// provisional allowance for extra variability.
#[derive(Debug, Clone)]
pub struct OntDedupParams {
    /// Endpoint displacement (bases) allowed at each end at no cost.
    pub max_end_offset: usize,
    /// Number of bases inward from each end used for the endpoint check.
    pub end_check_len: usize,
    /// Edits allowed within each endpoint check window.
    pub max_end_errors: usize,
    /// Upper cap on whole-read edit distance, as a fraction of the longer read.
    pub max_error_frac: f64,
    /// Strain divergence against which to protect (0 < value <= 1).
    pub protect_divergence: f64,
    /// Positive multiplier on usable-quality expected errors (default 1.2).
    /// Does not scale the assumed-rate contribution from unavailable qualities.
    pub quality_calibration: f64,
    /// Positive multiplier on usable-quality noise (default 1.0).
    pub noise_scale: f64,
    /// Additive errors per base with usable qualities (default 0.002).
    pub quality_error_floor: f64,
    /// Additional relative uncertainty of pair distance (finite, nonnegative).
    pub pair_distance_cv: f64,
    /// Seed-offset drift allowance, independent of the whole-read cutoff.
    pub seed_drift_frac: f64,
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
            end_check_len: 32,
            max_end_errors: 3,
            max_error_frac: 0.12,
            protect_divergence: 0.015,
            quality_calibration: 1.2,
            noise_scale: 1.0,
            quality_error_floor: 0.002,
            pair_distance_cv: 0.25,
            seed_drift_frac: 0.04,
            assumed_error_rate: 0.01,
            rc_tolerant: true,
        }
    }
}

impl OntDedupParams {
    /// Validate all statistical parameters before processing reads.
    pub fn validate(&self) -> Result<(), String> {
        // Multipliers are capped generously; past these, noise estimates
        // overflow and merging would silently stop.
        for (name, value, positive, max) in [
            ("max_error_frac", self.max_error_frac, false, 1.0),
            ("protect_divergence", self.protect_divergence, true, 1.0),
            ("quality_calibration", self.quality_calibration, true, 100.0),
            ("noise_scale", self.noise_scale, true, 100.0),
            ("quality_error_floor", self.quality_error_floor, false, 1.0),
            ("pair_distance_cv", self.pair_distance_cv, false, 100.0),
            ("seed_drift_frac", self.seed_drift_frac, false, 1.0),
            ("assumed_error_rate", self.assumed_error_rate, false, 1.0),
        ] {
            if !value.is_finite() || value < 0.0 || value > max || (positive && value == 0.0) {
                let lower = if positive { "greater than" } else { "at least" };
                return Err(format!(
                    "{name} must be {lower} 0 and at most {max}, got {value}"
                ));
            }
        }
        Ok(())
    }
}

/// Whether every byte is a valid Phred+33 quality (33..=126).
/// Empty qualities are valid and use the assumed error rate.
pub fn validate_quality(qual: &[u8]) -> bool {
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
pub(crate) fn read_noise(len: usize, errors: f32, params: &OntDedupParams) -> Option<f64> {
    (errors >= 0.0).then(|| {
        params.noise_scale
            * (params.quality_calibration * f64::from(errors)
                + params.quality_error_floor * len as f64)
    })
}

/// Extra spread models noise misprediction, while divergence adds Poisson variance.
pub(crate) fn noise_variance(mu: f64, noise: f64, cv: f64) -> f64 {
    mu + (cv * noise).powi(2)
}

/// Upper cap; safe for a cheap precheck before calculating the pair budget.
pub(crate) fn whole_read_cap(a_len: usize, b_len: usize, params: &OntDedupParams) -> usize {
    let longer = a_len.max(b_len);
    (params.max_error_frac * longer as f64).floor() as usize
}

/// None represents T=-1. Unavailable summaries get an assumed expected-error
/// count, not a cap-sized acceptance budget. Arithmetic after summary rounding
/// is f64 and symmetric, including when only one read has usable qualities.
/// Usable qualities contribute noise_scale * (c * expected_errors + f * length).
/// Unavailable qualities contribute assumed_error_rate * length, without
/// quality terms or noise_scale.
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
    let errors =
        |len, e| read_noise(len, e, params).unwrap_or(params.assumed_error_rate * len as f64);
    let noise = errors(a_len, a_errors) + errors(b_len, b_errors);
    let mu = noise + params.protect_divergence * a_len.min(b_len) as f64;
    // The extra spread models misprediction of the noise, so it scales with
    // the noise only: variance = mu + (pair_distance_cv * noise)^2. Applied to
    // all of mu, it would keep the cutoff a fixed fraction below the mean at
    // every length, and longer reads would gain no power.
    lower_tail_budget(
        mu,
        noise_variance(mu, noise, params.pair_distance_cv),
        whole_read_cap(a_len, b_len, params),
    )
}

/// Largest k <= cap with F(k) < .05 for NB(mean, variance), or Poisson
/// when variance equals the mean.
/// Sum relative masses starting at one, keeping log(P(0)) separately. Rescale
/// before masses grow too large. Never initialize a recurrence with exp(logP0),
/// which can underflow and silently turn an entire CDF into zero.
fn lower_tail_budget(mu: f64, variance: f64, cap: usize) -> Option<usize> {
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
    let mut threshold = (0.05_f64.ln() - log_scale).exp();
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
            threshold = (0.05_f64.ln() - log_scale).exp();
        }
    }
    Some(cap)
}

/// Assumed-rate comparator, with the production order: cap, ends, budget,
/// length compatibility, then whole-read alignment.
#[cfg(test)]
pub(crate) fn is_duplicate(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
    let cap = whole_read_cap(a.len(), b.len(), params);
    if a.len().abs_diff(b.len()) > cap {
        return false;
    }
    let (b_start, b_end) = end_windows(b.len(), params);
    if !ends_match(a, &b[b_start], &b[b_end], params) {
        return false;
    }
    let budget = pair_budget(a.len(), b.len(), -1.0, -1.0, params);
    budget.is_some_and(|t| lengths_compatible(a.len(), b.len(), t) && whole_read_matches(a, b, t))
}

/// Edit distance is at least the length difference, so this cheap check can
/// reject before any alignment.
pub(crate) fn lengths_compatible(a_len: usize, b_len: usize, budget: usize) -> bool {
    a_len.abs_diff(b_len) <= budget
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

#[cfg(test)]
pub(crate) fn whole_read_matches(a: &[u8], b: &[u8], budget: usize) -> bool {
    editdistancek::edit_distance_bounded(a, b, budget).is_some()
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
    use crate::ont::test_util::{cap_only, random_seq, substitute};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn params(max_end_offset: usize, max_end_errors: usize) -> OntDedupParams {
        OntDedupParams {
            max_end_offset,
            end_check_len: 64,
            max_end_errors,
            ..cap_only(0.05)
        }
    }

    fn ends_match_full(a: &[u8], b: &[u8], params: &OntDedupParams) -> bool {
        let (start, end) = end_windows(b.len(), params);
        ends_match(a, &b[start], &b[end], params)
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
        assert!(whole_read_matches(shifted, &seq, 100));
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
    fn default_endpoint_rule() {
        let d = OntDedupParams::default();
        assert_eq!(
            (d.max_end_offset, d.end_check_len, d.max_end_errors),
            (2, 32, 3)
        );
        assert_eq!((d.max_error_frac, d.rc_tolerant), (0.12, true));
        // With the defaults: 2 bases of displacement pass, 3 don't (even
        // though 3 edits are allowed), and 3 substitutions in the 32-base
        // window pass while 4 don't.
        let mut rng = StdRng::seed_from_u64(13);
        let seq = random_seq(1000, &mut rng);
        assert!(is_duplicate(&seq[2..], &seq, &d));
        assert!(!is_duplicate(&seq[3..], &seq, &d));
        let mut read = seq.clone();
        for pos in [4, 12, 20] {
            substitute(&mut read, pos);
        }
        assert!(is_duplicate(&read, &seq, &d));
        substitute(&mut read, 28);
        assert!(!is_duplicate(&read, &seq, &d));
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
        assert!(!whole_read_matches(&read, &seq, 50));
        assert!(whole_read_matches(&read, &seq, 60));
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

#[cfg(test)]
mod quality_tests {
    use super::*;
    fn cutoff(mu: f64, cv: f64, cap: usize) -> Option<usize> {
        lower_tail_budget(mu, noise_variance(mu, mu, cv), cap)
    }

    /// The budget for a pair with this noise and divergence term, cv 0.25.
    fn expected(noise: f64, divergence: f64, cap: usize) -> Option<usize> {
        let mu = noise + divergence;
        cutoff(mu, 0.25 * noise / mu, cap)
    }

    #[test]
    fn calibration_and_floor_apply_only_to_usable_qualities() {
        let mut p = OntDedupParams {
            protect_divergence: 0.01,
            ..Default::default()
        };
        assert_eq!(p.quality_calibration, 1.2);
        for c in [0.5, 1.0, 2.0, 5.0] {
            p.quality_calibration = c;
            let budget = |ea, eb| pair_budget(1000, 1200, ea, eb, &p);
            assert_eq!(budget(-1.0, -1.0), expected(22.0, 10.0, 144));
            assert_eq!(
                budget(10.0, -1.0),
                expected(c * 10.0 + 2.0 + 12.0, 10.0, 144)
            );
            assert_eq!(
                budget(-1.0, 12.0),
                expected(10.0 + c * 12.0 + 2.4, 10.0, 144)
            );
            assert_eq!(budget(10.0, 12.0), expected(c * 22.0 + 4.4, 10.0, 144));
        }
    }

    #[test]
    fn noise_scale_multiplies_both_quality_terms_only() {
        let mut p = OntDedupParams::default();
        for scale in [0.5, 1.0, 4.0] {
            p.noise_scale = scale;
            assert_eq!(read_noise(1000, -1.0, &p), None);
            assert_eq!(read_noise(1000, 10.0, &p), Some(scale * 14.0));
            assert_eq!(
                pair_budget(1000, 1000, -1.0, -1.0, &p),
                pair_budget(1000, 1000, -1.0, -1.0, &OntDedupParams::default())
            );
        }
    }

    #[test]
    fn validates_statistical_parameters() {
        type Setter = fn(&mut OntDedupParams, f64);
        // (setter, must be positive, upper bound)
        let fields: [(Setter, bool, f64); 8] = [
            (|p, v| p.max_error_frac = v, false, 1.0),
            (|p, v| p.protect_divergence = v, true, 1.0),
            (|p, v| p.quality_calibration = v, true, 100.0),
            (|p, v| p.noise_scale = v, true, 100.0),
            (|p, v| p.quality_error_floor = v, false, 1.0),
            (|p, v| p.pair_distance_cv = v, false, 100.0),
            (|p, v| p.seed_drift_frac = v, false, 1.0),
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
            quality_calibration: 1.0,
            protect_divergence: 0.01,
            ..Default::default()
        };
        // Expected noise is 10 + 10 + 2 * 0.002 * 1000 = 24 edits.
        assert!(pair_budget(1000, 1000, 10.0, 10.0, &p).unwrap() < 24);
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
            quality_error_floor: 0.0,
            quality_calibration: 1.0,
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
        // Unavailable qualities get E=.01L and the same protected cutoff.
        assert_eq!(
            pair_budget(1000, 1200, -1.0, -1.0, &p),
            pair_budget(1000, 1200, 10.0, 12.0, &p)
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
        let setters: [Set; 7] = [
            |_, a, _, i| *a = i as f32,
            |_, _, b, i| *b = i as f32,
            |p, _, _, i| p.quality_calibration = i as f64 / 20.0,
            |p, _, _, i| p.quality_error_floor = i as f64 / 5000.0,
            |p, _, _, i| p.protect_divergence = i as f64 / 1000.0,
            |p, _, _, i| p.max_error_frac = i as f64 / 1000.0,
            |p, _, _, i| p.noise_scale = i as f64 / 20.0,
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
            expected(22.0, 0.01 * 1000.0, 144)
        );
    }

    #[test]
    fn lengths_use_exact_integer_budget() {
        assert!(lengths_compatible(1000, 1015, 15));
        assert!(!lengths_compatible(1000, 1016, 15));
        assert!(lengths_compatible(1015, 1000, 15));
        assert_eq!(
            pair_budget(0, 0, -1.0, -1.0, &OntDedupParams::default()),
            None
        );
    }
}
