//! Replayable calibration of one multiplier on the configured noise shape.

use super::compare::{noise_variance, read_noise};
use super::{OntDedupContext, OntDedupParams, OntMinimizerParams};

/// Whole-read distance cap in the pilot, independent of production parameters.
const PILOT_CAP_FRAC: f64 = 0.12;
const QUALITY_PROBE_READS: usize = 20_000;
const MIN_BASELINE_NOISE: f64 = 10.0;

/// Collection limits for a calibration pilot.
#[derive(Debug, Clone)]
pub struct CalibrationOptions {
    /// Retained pairs required to stop (default 500).
    pub target_pairs: usize,
    /// Retained pairs required for an estimate (default 200).
    pub min_pairs: usize,
    /// Maximum examined reads (default 200,000).
    pub max_reads: usize,
    /// Sequence-plus-quality bytes before stopping (default 256 MiB).
    /// The last read can cross this limit.
    pub max_bytes: usize,
}

impl Default for CalibrationOptions {
    fn default() -> Self {
        Self {
            target_pairs: 500,
            min_pairs: 200,
            max_reads: 200_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

impl CalibrationOptions {
    pub fn validate(&self) -> Result<(), String> {
        if self.min_pairs == 0
            || self.target_pairs < self.min_pairs
            || self.max_reads == 0
            || self.max_bytes == 0
        {
            return Err("invalid calibration limits or pair thresholds".into());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationStatus {
    Estimated,
    NoUsableQualities,
    InsufficientPairs,
    NumericalFailure,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    NoUsableQualities,
    TargetPairs,
    ReadLimit,
    ByteLimit,
}

/// Noise multiplier and pilot support. Apply the multiplier as `noise_scale`.
#[derive(Debug, Clone, PartialEq)]
pub struct CalibrationReport {
    pub multiplier: f64,
    pub pairs_used: usize,
    pub reads_examined: usize,
    pub status: CalibrationStatus,
    /// None means input ended before a stopping rule fired.
    pub stop_reason: Option<StopReason>,
}

impl std::fmt::Display for CalibrationReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = match self.status {
            CalibrationStatus::Estimated => "estimated",
            CalibrationStatus::NoUsableQualities => "fallback(no usable qualities)",
            CalibrationStatus::InsufficientPairs => "fallback(insufficient eligible pairs)",
            CalibrationStatus::NumericalFailure => "fallback(numerical failure)",
        };
        let stop = match self.stop_reason {
            Some(StopReason::NoUsableQualities) => "no-usable-qualities",
            Some(StopReason::TargetPairs) => "target-pairs",
            Some(StopReason::ReadLimit) => "read-limit",
            Some(StopReason::ByteLimit) => "byte-limit",
            None => "none",
        };
        write!(
            f,
            "Noise scale: status={status}; a={:.4}; pairs={}; reads={}; stop={stop}",
            self.multiplier, self.pairs_used, self.reads_examined
        )
    }
}

/// Incremental pilot. After finishing, replay sampled reads into a fresh
/// context using the report's multiplier as `noise_scale`.
pub struct NoiseScalePilot {
    context: OntDedupContext,
    opts: CalibrationOptions,
    /// (Distance / baseline noise, baseline noise / standard deviation).
    pairs: Vec<(f64, f64)>,
    has_pair: Vec<bool>,
    reads: usize,
    bytes: usize,
    has_usable_qualities: bool,
    baseline: OntDedupParams,
}

impl NoiseScalePilot {
    /// Panics on invalid parameters or options.
    pub fn new(
        dedup: &OntDedupParams,
        minimizer: &OntMinimizerParams,
        opts: &CalibrationOptions,
    ) -> Self {
        opts.validate().expect("invalid CalibrationOptions");
        dedup.validate().expect("invalid OntDedupParams");
        let baseline = OntDedupParams {
            noise_scale: 1.0,
            ..dedup.clone()
        };
        Self {
            context: OntDedupContext::new(
                OntDedupParams {
                    // Divergence and zero CV keep the distance cap binding.
                    // Exact endpoints make copies of one molecule much more
                    // likely; allowing offsets admits more chance coincidences.
                    protect_divergence: 1.0,
                    pair_distance_cv: 0.0,
                    max_end_offset: 0,
                    max_error_frac: PILOT_CAP_FRAC,
                    ..baseline.clone()
                },
                minimizer.clone(),
            ),
            opts: opts.clone(),
            pairs: Vec::new(),
            has_pair: Vec::new(),
            reads: 0,
            bytes: 0,
            has_usable_qualities: false,
            baseline,
        }
    }

    fn stop_reason(&self) -> Option<StopReason> {
        Some(if self.pairs.len() >= self.opts.target_pairs {
            StopReason::TargetPairs
        } else if self.reads >= self.opts.max_reads {
            StopReason::ReadLimit
        } else if self.bytes >= self.opts.max_bytes {
            StopReason::ByteLimit
        } else if self.reads >= QUALITY_PROBE_READS && !self.has_usable_qualities {
            StopReason::NoUsableQualities
        } else {
            return None;
        })
    }

    pub fn is_finished(&self) -> bool {
        self.stop_reason().is_some()
    }

    /// Process one read; panics after a stopping rule fires or on malformed qualities.
    pub fn process_read(&mut self, sequence: &[u8], quality: &[u8]) {
        assert!(!self.is_finished(), "calibration pilot already finished");
        let (_, matched, errors) = self.context.process_read(self.reads, sequence, quality);
        self.reads += 1;
        self.bytes = self
            .bytes
            .saturating_add(sequence.len())
            .saturating_add(quality.len());
        self.has_usable_qualities |= errors >= 0.0;
        let Some(matched) = matched else { return };
        let (Some(read_noise), Some(leader_noise)) = (
            read_noise(sequence.len(), errors, &self.baseline),
            read_noise(matched.len, matched.errors, &self.baseline),
        ) else {
            return;
        };
        let noise = read_noise + leader_noise;
        if noise < MIN_BASELINE_NOISE || !noise.is_finite() {
            return;
        }
        self.retain(matched.cluster, matched.distance as f64, noise);
    }

    fn retain(&mut self, leader: u32, distance: f64, noise: f64) {
        let leader = leader as usize;
        if self.has_pair.len() <= leader {
            self.has_pair.resize(leader + 1, false);
        }
        if !self.has_pair[leader] {
            self.has_pair[leader] = true;
            let s = noise_variance(noise, noise, self.baseline.pair_distance_cv)
                .sqrt()
                .max(1.0);
            self.pairs.push((distance / noise, noise / s));
        }
    }

    /// Resolve the estimate after input ends or a stopping rule fires.
    pub fn finish(mut self) -> CalibrationReport {
        let mut report = CalibrationReport {
            multiplier: 1.0,
            pairs_used: self.pairs.len(),
            reads_examined: self.reads,
            status: if !self.has_usable_qualities {
                CalibrationStatus::NoUsableQualities
            } else {
                CalibrationStatus::InsufficientPairs
            },
            stop_reason: self.stop_reason(),
        };
        if self.pairs.len() < self.opts.min_pairs {
            return report;
        }
        let total_weight: f64 = self.pairs.iter().map(|p| p.1).sum();
        if total_weight <= 0.0
            || !total_weight.is_finite()
            || self
                .pairs
                .iter()
                .any(|p| !p.0.is_finite() || !p.1.is_finite() || p.1 <= 0.0)
        {
            report.status = CalibrationStatus::NumericalFailure;
            return report;
        }
        self.pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        let mut cumulative = 0.0;
        for &(ratio, weight) in &self.pairs {
            cumulative += weight;
            if cumulative >= total_weight / 2.0 {
                let n = self.pairs.len() as f64;
                let shrunk = (n * ratio.ln() / (n + 50.0)).exp();
                // Favor precision: very noisy inputs may retain duplicates
                // rather than loosen thresholds enough to merge strains.
                report.multiplier = shrunk.clamp(0.5, 2.0);
                report.status = CalibrationStatus::Estimated;
                break;
            }
        }
        report
    }
}

/// Estimate from sequence/quality slices, stopping at the pilot's limits or EOF.
/// Invalid inputs panic. Callers must buffer sampled reads themselves to replay
/// them into a fresh production context using the report's multiplier.
pub fn estimate_noise_scale<'a>(
    reads: impl IntoIterator<Item = (&'a [u8], &'a [u8])>,
    dedup: &OntDedupParams,
    minimizer: &OntMinimizerParams,
    opts: &CalibrationOptions,
) -> CalibrationReport {
    let mut pilot = NoiseScalePilot::new(dedup, minimizer, opts);
    for (sequence, quality) in reads {
        pilot.process_read(sequence, quality);
        if pilot.is_finished() {
            break;
        }
    }
    pilot.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ont::compare::expected_errors;
    use crate::ont::test_util::{random_seq, substitute, varied_quality};
    use rand::{SeedableRng, rngs::StdRng};

    fn pilot(opts: &CalibrationOptions) -> NoiseScalePilot {
        NoiseScalePilot::new(
            &OntDedupParams::default(),
            &OntMinimizerParams::default(),
            opts,
        )
    }

    #[test]
    fn eligibility_uses_configured_coefficients() {
        let mut rng = StdRng::seed_from_u64(19);
        let seq = random_seq(1200, &mut rng);
        let quality = varied_quality(1200, b'5');
        let params = OntDedupParams {
            quality_calibration: 0.01,
            quality_error_floor: 0.0,
            ..Default::default()
        };
        let mut p = NoiseScalePilot::new(
            &params,
            &OntMinimizerParams::default(),
            &CalibrationOptions::default(),
        );
        p.process_read(&seq, &quality);
        p.process_read(&seq, &quality);
        assert_eq!(p.context.stats().1, 1);
        assert!(p.pairs.is_empty());
    }

    #[test]
    fn pilot_uses_independent_cap_and_unscaled_baseline() {
        let mut rng = StdRng::seed_from_u64(91);
        let seq = random_seq(1200, &mut rng);
        let quality = varied_quality(1200, b'5');
        let params = OntDedupParams {
            max_error_frac: 0.0,
            noise_scale: 4.0,
            quality_calibration: 0.01,
            quality_error_floor: 0.006,
            ..Default::default()
        };
        let mut p = NoiseScalePilot::new(
            &params,
            &OntMinimizerParams::default(),
            &CalibrationOptions::default(),
        );
        p.process_read(&seq, &quality);
        let mut altered = seq.clone();
        for i in (400..800).step_by(10) {
            substitute(&mut altered, i);
        }
        p.process_read(&altered, &quality);
        let observation = p.pairs[0];
        let noise = params.quality_calibration * 2.0 * f64::from(expected_errors(1200, &quality))
            + params.quality_error_floor * 2400.0;
        assert!((observation.0 - 40.0 / noise).abs() < 1e-12);
        assert_eq!(p.pairs.len(), 1);
        assert_eq!(p.context.stats().1, 1);
        altered[300..900].fill(b'N');
        p.process_read(&altered, &quality);
        assert_eq!(p.context.stats().1, 2);
        assert_eq!(p.pairs.len(), 1);
    }

    #[test]
    fn pilot_assigns_only_zero_offset_pairs() {
        let mut rng = StdRng::seed_from_u64(78);
        let seq = random_seq(1200, &mut rng);
        let quality = varied_quality(1200, b'5');
        let params = OntDedupParams {
            max_end_offset: 10,
            ..Default::default()
        };
        let mut p = NoiseScalePilot::new(
            &params,
            &OntMinimizerParams::default(),
            &CalibrationOptions::default(),
        );
        p.process_read(&seq, &quality);
        p.process_read(&seq[5..], &quality[5..]);
        assert_eq!(p.context.stats().1, 2);
        assert!(p.pairs.is_empty());
        let rc = crate::ont::reverse_complement(&seq);
        p.process_read(&rc, &quality);
        assert_eq!(p.pairs[0].0, 0.0);
        assert_eq!(p.context.stats().1, 2);
        // A genuine indel inside the read contributes to D without subtraction.
        let mut deletion = seq.clone();
        deletion.remove(600);
        let mut p = pilot(&CalibrationOptions::default());
        p.process_read(&seq, &quality);
        p.process_read(&deletion, &quality[1..]);
        assert_eq!(p.context.stats().1, 1);
        let observation = p.pairs[0];
        let params = OntDedupParams::default();
        let noise = params.quality_calibration
            * (f64::from(expected_errors(seq.len(), &quality))
                + f64::from(expected_errors(deletion.len(), &quality[1..])))
            + params.quality_error_floor * (seq.len() + deletion.len()) as f64;
        assert!((observation.0 - 1.0 / noise).abs() < 1e-12);
    }

    #[test]
    fn pilot_does_not_search_past_an_accepted_but_ineligible_leader() {
        let mut rng = StdRng::seed_from_u64(401);
        let seq = random_seq(1200, &mut rng);
        let quality = varied_quality(1200, b'5');
        let mut other = seq.clone();
        other[400..560].fill(b'N');
        let mut query = seq.clone();
        query[400..480].fill(b'N');
        let mut p = pilot(&CalibrationOptions::default());
        p.process_read(&seq, b"");
        p.process_read(&other, &quality);
        assert_eq!(p.context.stats().1, 2);
        // The first acceptable leader wins even when its qualities are unavailable.
        p.process_read(&query, &quality);
        assert_eq!(p.context.stats().1, 2);
        assert!(p.pairs.is_empty());
    }

    fn weighted_report(cv: f64, ratio: Option<f64>, count: usize) -> CalibrationReport {
        let params = OntDedupParams {
            quality_calibration: 2.3,
            quality_error_floor: 0.0007,
            pair_distance_cv: cv,
            ..Default::default()
        };
        let mut p = NoiseScalePilot::new(
            &params,
            &OntMinimizerParams::default(),
            &CalibrationOptions::default(),
        );
        p.has_usable_qualities = true;
        for i in 0..count {
            let noise = if i < 120 { 10.0 } else { 100.0 };
            let ratio = ratio.unwrap_or(if i < 120 { 0.8 } else { 2.0 });
            p.retain(i as u32, ratio * noise, noise);
            p.retain(i as u32, 1000.0, noise);
        }
        p.finish()
    }

    #[test]
    fn multiplier_uses_weighted_median_configured_cv_and_log_shrinkage() {
        for (cv, raw) in [(0.0, 2.0_f64), (0.7, 0.8_f64)] {
            let report = weighted_report(cv, None, 200);
            let expected = raw.powf(200.0 / 250.0);
            assert_eq!(report.status, CalibrationStatus::Estimated);
            assert_eq!(report.pairs_used, 200);
            assert!((report.multiplier - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn insufficient_support_returns_one_and_bounds_apply_after_shrinkage() {
        let fallback = weighted_report(0.25, Some(2.0), 199);
        assert_eq!(fallback.status, CalibrationStatus::InsufficientPairs);
        assert_eq!(fallback.multiplier, 1.0);
        for (ratio, expected) in [
            (0.0, 0.5),
            (0.25, 0.5),
            (2.0, 2.0_f64.powf(0.8)),
            (4.0, 2.0),
            (8.0, 2.0),
        ] {
            let report = weighted_report(0.25, Some(ratio), 200);
            assert!((report.multiplier - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn stopping_rules_do_not_overconsume_and_input_end_has_no_stop_reason() {
        use std::cell::Cell;
        let mut rng = StdRng::seed_from_u64(8);
        let seq = random_seq(3000, &mut rng);
        let q = varied_quality(3000, b'I');
        for (target_pairs, max_reads, max_bytes, expected, stop) in [
            (1, 10, usize::MAX, 2, StopReason::TargetPairs),
            (2, 3, usize::MAX, 3, StopReason::ReadLimit),
            (2, 10, 6001, 2, StopReason::ByteLimit),
        ] {
            let opts = CalibrationOptions {
                target_pairs,
                max_reads,
                max_bytes,
                min_pairs: 1,
            };
            let consumed = Cell::new(0);
            let reads = std::iter::repeat((seq.as_slice(), q.as_slice()))
                .inspect(|_| consumed.set(consumed.get() + 1));
            let report = estimate_noise_scale(
                reads,
                &OntDedupParams::default(),
                &OntMinimizerParams::default(),
                &opts,
            );
            assert_eq!(consumed.get(), expected);
            assert_eq!(report.stop_reason, Some(stop));
            assert_eq!(report.pairs_used, 1);
        }
        let mut p = pilot(&CalibrationOptions::default());
        p.process_read(&seq, &q);
        assert_eq!(p.finish().stop_reason, None);
        let report = estimate_noise_scale(
            [(seq.as_slice(), q.as_slice()); 3],
            &OntDedupParams::default(),
            &OntMinimizerParams::default(),
            &CalibrationOptions::default(),
        );
        assert_eq!(report.stop_reason, None);
        assert_eq!(report.reads_examined, 3);
    }

    #[test]
    fn quality_probe_stops_without_usable_qualities_but_accepts_a_mixed_prefix() {
        for usable in [false, true] {
            let mut p = pilot(&CalibrationOptions::default());
            p.process_read(b"ACGT", if usable { b"IJIJ" } else { b"IIII" });
            for _ in 1..QUALITY_PROBE_READS {
                p.process_read(b"ACGT", b"IIII");
            }
            assert_eq!(p.is_finished(), !usable);
            let report = p.finish();
            assert_eq!(
                report.stop_reason,
                if usable {
                    None
                } else {
                    Some(StopReason::NoUsableQualities)
                }
            );
        }
    }

    #[test]
    fn invalid_options_and_malformed_qualities_are_errors() {
        for opts in [
            CalibrationOptions {
                max_reads: 0,
                ..Default::default()
            },
            CalibrationOptions {
                max_bytes: 0,
                ..Default::default()
            },
            CalibrationOptions {
                min_pairs: 0,
                ..Default::default()
            },
            CalibrationOptions {
                target_pairs: 199,
                ..Default::default()
            },
        ] {
            assert!(std::panic::catch_unwind(|| pilot(&opts)).is_err());
        }
        assert!(
            std::panic::catch_unwind(|| {
                NoiseScalePilot::new(
                    &OntDedupParams {
                        max_error_frac: f64::NAN,
                        ..Default::default()
                    },
                    &OntMinimizerParams::default(),
                    &CalibrationOptions::default(),
                );
            })
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| {
                pilot(&CalibrationOptions::default()).process_read(b"ACGT", b"  II");
            })
            .is_err()
        );
    }
}
