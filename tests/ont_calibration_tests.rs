mod common;
use common::random_seq;
use nao_dedup::ont::{
    CalibrationOptions, CalibrationStatus, OntDedupContext, OntDedupParams, OntMinimizerParams,
    estimate_noise_scale,
};
use ont_sim::{QualityConfig, add_quality_errors};
use rand::{SeedableRng, rngs::StdRng, seq::SliceRandom};

type Read = (Vec<u8>, Vec<u8>);

fn read(source: &[u8], median: f64, factor: f64, rng: &mut StdRng) -> Read {
    let model = QualityConfig::enabled(median, 0.1, 1.0, factor, 0.0);
    let (sequence, quality, _) = add_quality_errors(source, &model, [0.4, 0.3, 0.3], 1.0, rng);
    (sequence, quality)
}

#[test]
fn indel_noise_level_uses_emitted_quality_predictors() {
    let minimizer = OntMinimizerParams::default();
    let opts = CalibrationOptions::default();
    for factor in [0.7, 2.0] {
        let mut rng = StdRng::seed_from_u64(209);
        let mut fitting = Vec::new();
        for family in 0..240 {
            let source = random_seq(if family % 2 == 0 { 1200 } else { 2400 }, &mut rng);
            fitting.push(read(&source, 0.004, factor, &mut rng));
            fitting.push(read(&source, 0.008, factor, &mut rng));
        }
        let emitted_rate = fitting
            .iter()
            .flat_map(|(_, q)| q)
            .map(|q| 10_f64.powf(-f64::from(q - 33) / 10.0))
            .sum::<f64>()
            / fitting.iter().map(|(s, _)| s.len()).sum::<usize>() as f64;
        fitting.shuffle(&mut rng);
        let params = OntDedupParams::default();
        let report = estimate_noise_scale(
            fitting.iter().map(|(s, q)| (s.as_slice(), q.as_slice())),
            &params,
            &minimizer,
            &opts,
        );
        assert_eq!(
            report.status,
            CalibrationStatus::Estimated,
            "{factor}: {report:?}"
        );
        // Read-to-read edit distance, cap/endpoint selection and shrinkage
        // make this a robust heuristic rather than an unbiased mean fit.
        assert!(
            (report.multiplier
                * (params.quality_calibration * emitted_rate + params.quality_error_floor)
                - factor * emitted_rate)
                .abs()
                < 0.0035,
            "{factor}: {report:?}"
        );
        assert_eq!(params.pair_distance_cv, 0.25);

        // Independent reads assess the complete calibrated pipeline. Endpoint
        // coincidence is deliberate for both duplicates and strain near misses.
        let resolved = OntDedupParams {
            noise_scale: report.multiplier,
            ..params
        };
        let (mut recalled, mut false_merges) = (0, 0);
        for _ in 0..40 {
            let source = random_seq(6000, &mut rng);
            let mut variant = source.clone();
            // A strain exactly `protect_divergence` away.
            let step = (1.0 / resolved.protect_divergence).round() as usize;
            for i in (step..6000 - step).step_by(step) {
                variant[i] = if variant[i] == b'A' { b'C' } else { b'A' };
            }
            let a = read(&source, 0.001, factor, &mut rng);
            let b = read(&source, 0.002, factor, &mut rng);
            let near = read(&variant, 0.002, factor, &mut rng);
            let merges = |other: &Read| {
                let mut ctx = OntDedupContext::new(resolved.clone(), minimizer.clone());
                ctx.process_read_by_index(0, &a.0, &a.1);
                ctx.process_read_by_index(1, &other.0, &other.1) == 0
            };
            recalled += usize::from(merges(&b));
            false_merges += usize::from(merges(&near));
        }
        assert!(recalled >= 12, "factor={factor}: recall {recalled}/40");
        // Exactly at the protected divergence the nominal merge rate is 5%;
        // allow sampling noise around it.
        assert!(
            false_merges <= 6,
            "factor={factor}: false merges {false_merges}/40"
        );
    }
}

#[test]
fn empty_flat_and_mixed_qualities_keep_only_usable_pairs() {
    let mut rng = StdRng::seed_from_u64(77);
    let source = random_seq(1200, &mut rng);
    let varied: Vec<u8> = (0..1200).map(|i| b'5' + (i % 2) as u8).collect();
    let flat = vec![b'5'; 1200];
    for qualities in [
        vec![&b""[..], &b""[..]],
        vec![flat.as_slice(), varied.as_slice()],
        vec![varied.as_slice(), flat.as_slice(), varied.as_slice()],
    ] {
        let expected = usize::from(qualities.len() == 3);
        let report = estimate_noise_scale(
            qualities.into_iter().map(|q| (source.as_slice(), q)),
            &OntDedupParams::default(),
            &OntMinimizerParams::default(),
            &CalibrationOptions::default(),
        );
        assert_eq!(report.pairs_used, expected);
        assert!(report.status != CalibrationStatus::Estimated);
    }
    let empty = estimate_noise_scale(
        [],
        &OntDedupParams::default(),
        &OntMinimizerParams::default(),
        &CalibrationOptions::default(),
    );
    assert_eq!(empty.reads_examined, 0);
    assert_eq!(empty.multiplier, 1.0);
}

#[test]
fn prefix_order_changes_support_and_later_pairs_can_reach_the_minimum() {
    let mut rng = StdRng::seed_from_u64(720);
    let reads: Vec<_> = (0..400).map(|_| random_seq(1200, &mut rng)).collect();
    let quality: Vec<u8> = (0..1200).map(|i| b'5' + (i % 2) as u8).collect();
    let adjacent: Vec<_> = reads.iter().flat_map(|s| [s, s]).collect();
    let separated: Vec<_> = reads.iter().chain(reads.iter()).collect();
    let opts = CalibrationOptions {
        max_reads: 400,
        target_pairs: 200,
        ..Default::default()
    };
    let estimate = |order: &[&Vec<u8>], opts: &CalibrationOptions| {
        estimate_noise_scale(
            order.iter().map(|s| (s.as_slice(), quality.as_slice())),
            &OntDedupParams::default(),
            &OntMinimizerParams::default(),
            opts,
        )
    };
    let supported = estimate(&adjacent, &opts);
    assert_eq!(supported.status, CalibrationStatus::Estimated);
    assert_eq!(supported.reads_examined, 400);
    assert_eq!(supported.pairs_used, 200);
    let fallback = estimate(&separated, &opts);
    assert!(fallback.status != CalibrationStatus::Estimated);
    assert_eq!(fallback.pairs_used, 0);
    let longer = estimate(
        &separated,
        &CalibrationOptions {
            max_reads: 600,
            ..opts
        },
    );
    assert_eq!(longer.reads_examined, 600);
    assert_eq!(longer.pairs_used, 200);
    assert_eq!(longer.status, CalibrationStatus::Estimated);
}
