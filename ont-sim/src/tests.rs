use super::*;
use clap::Parser;

#[derive(Parser)]
struct Args {
    #[command(flatten)]
    config: SimConfig,
}

fn config() -> SimConfig {
    Args::parse_from([
        "sim",
        "--sub-rate",
        "0",
        "--ins-rate",
        "0",
        "--del-rate",
        "0",
    ])
    .config
}

fn truth(i: usize, molecule: usize) -> Truth {
    Truth {
        read_idx: i,
        molecule_id: molecule,
        genome_id: 0,
        mol_start: molecule * 10,
        mol_end: molecule * 10 + 100,
        is_rc: false,
        n_sub: 0,
        n_ins: 0,
        n_del: 0,
    }
}

#[test]
fn evidence_metrics_count_molecules_and_genomes_once() {
    let t: Vec<_> = (0..4).map(|i| truth(i, i / 2)).collect();
    let mut cross = t.clone();
    cross[2].genome_id = 1;
    cross[3].genome_id = 1;
    for (truth, assignments, counts, rates) in [
        (
            &t,
            vec![(0, 0), (1, 0), (2, 2), (3, 2)],
            [0, 0, 0, 0],
            ["0", "0"],
        ),
        (
            &t,
            vec![(0, 0), (1, 0), (2, 0), (3, 0)],
            [1, 0, 0, 0],
            ["0.5", "0"],
        ),
        (
            &cross,
            vec![(0, 0), (1, 0), (2, 0), (3, 0)],
            [1, 1, 1, 0],
            ["0.5", "0"],
        ),
        (
            &t,
            vec![(0, 0), (1, 1), (2, 2), (3, 2)],
            [0, 0, 0, 1],
            ["0", "0.5"],
        ),
        // Molecule 0 is split, but still represented by exemplar 0.
        (
            &cross,
            vec![(0, 0), (1, 2), (2, 2), (3, 2)],
            [0, 0, 0, 0],
            ["0", "0"],
        ),
        (&vec![], vec![], [0, 0, 0, 0], ["null", "null"]),
    ] {
        let e = evaluate(truth, &assignments).unwrap();
        assert_eq!(
            [
                e.lost_molecules,
                e.lost_molecules_cross_genome,
                e.lost_genomes,
                e.extra_exemplars
            ],
            counts
        );
        assert_eq!(e.lost_molecule_rate(), rates[0].parse::<f64>().ok());
        assert_eq!(e.missed_duplicate_rate(), rates[1].parse::<f64>().ok());
        let json = e.to_json();
        for (key, value) in [
            "lost_molecules",
            "lost_molecules_cross_genome",
            "lost_genomes",
            "extra_exemplars",
        ]
        .into_iter()
        .zip(counts)
        {
            assert!(json.contains(&format!("\"{key}\":{value}")), "{json}");
        }
        for (key, value) in ["lost_molecule_rate", "missed_duplicate_rate"]
            .into_iter()
            .zip(rates)
        {
            assert!(json.contains(&format!("\"{key}\":{value}")), "{json}");
        }
        assert!(e.to_string().starts_with("lost_molecules="));
    }
    let e = evaluate(&[truth(10, 0), truth(20, 1)], &[(20, 10), (10, 10)]).unwrap();
    assert!(e.to_json().contains("\"missed_duplicate_rate\":null"));
}

#[test]
fn lost_molecule_can_span_clusters_while_its_genome_survives() {
    let mut t: Vec<_> = (0..6).map(|i| truth(i, i / 2)).collect();
    for row in &mut t[2..] {
        row.genome_id = 1;
    }
    // Molecule 1 loses both reads, one across genomes and one within genome 1.
    // Genome 1 still has exemplar 4 from molecule 2.
    let e = evaluate(&t, &[(0, 0), (1, 0), (2, 0), (3, 4), (4, 4), (5, 4)]).unwrap();
    assert_eq!(
        (
            e.lost_molecules,
            e.lost_molecules_cross_genome,
            e.lost_genomes,
            e.extra_exemplars
        ),
        (1, 1, 0, 0)
    );
}

#[test]
fn related_genomes_copy_reference_with_substitutions_only() {
    let parsed = Args::try_parse_from([
        "sim",
        "--sub-rate",
        "0",
        "--ins-rate",
        "0",
        "--del-rate",
        "0",
        "--genome-divergence",
        "0.25",
    ]);
    assert!(parsed.is_ok(), "--genome-divergence must be accepted");
    let mut c = parsed.unwrap().config;
    c.n_genomes = 3;
    c.genome_len = 10_000;
    c.n_molecules = 0;
    let s = simulate(&c).unwrap();
    assert_eq!(s, simulate(&c).unwrap());
    for genome in &s.genomes[1..] {
        assert_eq!(genome.len(), c.genome_len);
        let changed = genome
            .iter()
            .zip(&s.genomes[0])
            .filter(|(a, b)| a != b)
            .count();
        assert!((changed as f64 / c.genome_len as f64 - 0.25).abs() < 0.02);
    }
    let changed = s.genomes[1]
        .iter()
        .zip(&s.genomes[2])
        .filter(|(a, b)| a != b)
        .count();
    // Independent substitutions can coincide: divergence is 2p - 4p²/3.
    assert!((changed as f64 / c.genome_len as f64 - (0.5 - 0.25 / 3.0)).abs() < 0.02);
    c.genome_divergence = 1.0;
    let all_changed = simulate(&c).unwrap();
    for genome in &all_changed.genomes[1..] {
        assert!(
            genome
                .iter()
                .zip(&all_changed.genomes[0])
                .all(|(a, b)| a != b)
        );
    }
    c.genome_divergence = 0.0;
    assert_eq!(s.genomes[0], simulate(&c).unwrap().genomes[0]);
    for p in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        c.genome_divergence = p;
        assert!(simulate(&c).is_err());
    }
}

#[test]
fn deterministic_and_reverse_complement_truth() {
    assert_eq!(simulate::reverse_complement(b"ACGTT"), b"AACGT");
    let mut c = config();
    c.n_molecules = 40;
    c.family_mean = 3.0;
    let a = simulate(&c).unwrap();
    assert_eq!(a, simulate(&c).unwrap());
    assert!(a.reads.iter().any(|r| r.truth.is_rc));
    assert!(a.reads.iter().any(|r| !r.truth.is_rc));
    for (i, r) in a.reads.iter().enumerate() {
        let t = &r.truth;
        let mut expected = a.genomes[t.genome_id][t.mol_start..t.mol_end].to_vec();
        if t.is_rc {
            expected = simulate::reverse_complement(&expected);
        }
        assert_eq!(r.sequence, expected);
        assert_eq!(r.quality.len(), r.sequence.len());
        assert_eq!(t.read_idx, i);
    }
    c.seed += 1;
    assert_ne!(a, simulate(&c).unwrap());
}

#[test]
fn errors_follow_requested_probabilities() {
    let mut c = config();
    c.n_molecules = 200;
    c.genome_len = 20_000;
    c.len_min = 10_000;
    c.len_max = 10_000;
    c.family_mean = 1.0;
    c.sub_rate = 0.07;
    c.ins_rate = 0.04;
    c.del_rate = 0.03;
    let s = simulate(&c).unwrap();
    let (mut sub, mut ins, mut del, mut bases) = (0, 0, 0, 0);
    for r in s.reads {
        let t = r.truth;
        assert_eq!(
            r.sequence.len(),
            t.mol_end - t.mol_start + t.n_ins - t.n_del
        );
        sub += t.n_sub;
        ins += t.n_ins;
        del += t.n_del;
        bases += t.mol_end - t.mol_start;
    }
    for (observed, expected) in [
        (sub as f64 / (bases - del) as f64, c.sub_rate),
        (ins as f64 / bases as f64, c.ins_rate),
        (del as f64 / bases as f64, c.del_rate),
    ] {
        assert!(
            (observed - expected).abs() < 0.002,
            "{observed} vs {expected}"
        );
    }
}

#[test]
fn hotspot_collisions() {
    let mut c = config();
    c.n_molecules = 30;
    c.endpoints = Endpoints::Hotspot;
    c.n_hotspots = 1;
    c.hotspot_jitter = 0;
    let s = simulate(&c).unwrap();
    let truth: Vec<_> = s.reads.iter().map(|r| r.truth.clone()).collect();
    let assignments: Vec<_> = truth.iter().map(|t| (t.read_idx, 0)).collect();
    assert_eq!(
        evaluate(&truth, &assignments)
            .unwrap()
            .identifiability_classes,
        1
    );
}

#[test]
fn pair_metrics_and_identifiability() {
    let t: Vec<_> = (0..4).map(|i| truth(i, i / 2)).collect();
    let perfect = evaluate(&t, &[(0, 0), (1, 0), (2, 2), (3, 2)]).unwrap();
    assert_eq!(perfect.molecule.precision(), Some(1.0));
    assert_eq!(perfect.molecule.recall(), Some(1.0));
    let merged = evaluate(&t, &[(0, 0), (1, 0), (2, 0), (3, 0)]).unwrap();
    assert_eq!(merged.molecule.precision(), Some(1.0 / 3.0));
    assert_eq!(merged.molecule.recall(), Some(1.0));
    let split = evaluate(&t, &[(0, 0), (1, 1), (2, 2), (3, 3)]).unwrap();
    assert_eq!(split.molecule.precision(), None);
    assert_eq!(split.molecule.recall(), Some(0.0));
    let mut identical = t.clone();
    for row in &mut identical {
        row.mol_start = 0;
        row.mol_end = 100;
    }
    let e = evaluate(&identical, &[(0, 0), (1, 0), (2, 0), (3, 0)]).unwrap();
    assert_eq!(e.identifiability.precision(), Some(1.0));
    assert_eq!(e.identifiability.recall(), Some(1.0));
    assert_eq!(e.molecule.precision(), Some(1.0 / 3.0));
    let e = evaluate(&identical, &[(0, 0), (1, 1), (2, 2), (3, 3)]).unwrap();
    assert_eq!(e.identifiability.precision(), None);
    assert_eq!(e.identifiability.recall(), Some(0.0));
    let empty = evaluate(&[], &[]).unwrap();
    assert_eq!(empty.molecule.precision(), None);
    assert_eq!(empty.molecule.recall(), None);
    assert_eq!(empty.identifiability.precision(), None);
    assert_eq!(empty.identifiability.recall(), None);
}

#[test]
fn assignment_order_and_validation() {
    let mut t: Vec<_> = (0..4).map(|i| truth(i, 0)).collect();
    let a = [(3, 2), (0, 0), (2, 2), (1, 0)]; // input order must not matter
    let e = evaluate(&t, &a).unwrap();
    assert_eq!(e.molecule.true_pairs, 6);
    assert_eq!(e.molecule.predicted_pairs, 2);
    assert_eq!(e.molecule.recovered_pairs, 2);
    for bad in [
        vec![(0, 0)],
        vec![(0, 0), (1, 0), (2, 2), (2, 2)],
        vec![(0, 1), (1, 2), (2, 2), (3, 2)],
        vec![(0, 0), (1, 0), (2, 2), (4, 2)],
    ] {
        assert!(evaluate(&t, &bad).is_err());
    }
    t[1].mol_end += 1;
    assert!(evaluate(&t, &a).is_err()); // one molecule cannot have inconsistent coordinates
}

#[test]
fn end_jitter_preserves_the_interior() {
    let mut c = config();
    c.n_molecules = 20;
    c.rc_prob = 0.0;
    c.end_jitter = 7;
    let s = simulate(&c).unwrap();
    for r in s.reads {
        let t = &r.truth;
        let template = &s.genomes[t.genome_id][t.mol_start..t.mol_end];
        assert!(r.sequence.len().abs_diff(template.len()) <= 14);
        assert!(
            r.sequence
                .windows(template.len() - 14)
                .any(|w| w == &template[7..template.len() - 7])
        );
    }
}

#[test]
fn rejects_invalid_parameters_and_truth() {
    let mut c = config();
    for rate in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        c.ins_rate = rate;
        assert!(simulate(&c).is_err());
    }
    c.ins_rate = 0.0;
    c.family_mean = 0.5;
    assert!(simulate(&c).is_err());
    let mut t = vec![truth(10, 0), truth(20, 0)];
    assert_eq!(
        evaluate(&t, &[(20, 10), (10, 10)])
            .unwrap()
            .molecule
            .recall(),
        Some(1.0)
    );
    t.push(t[0].clone());
    assert!(evaluate(&t, &[(10, 10), (20, 10)]).is_err());
}

#[test]
fn errors_change_bases_and_homopolymers_modify_indels() {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(7);
    let mut c = config();
    c.sub_rate = 1.0;
    let (s, counts) = add_errors(
        b"AAACGTTT",
        c.sub_rate,
        c.ins_rate,
        c.del_rate,
        c.homopolymer_indel_mult,
        &mut rng,
    );
    assert_eq!(counts, [8, 0, 0]);
    assert!(s.iter().zip(b"AAACGTTT").all(|(a, b)| a != b));
    c.sub_rate = 0.0;
    c.ins_rate = 1.0;
    c.homopolymer_indel_mult = 0.0;
    let (s, counts) = add_errors(
        b"AAACGTTT",
        c.sub_rate,
        c.ins_rate,
        c.del_rate,
        c.homopolymer_indel_mult,
        &mut rng,
    );
    assert_eq!(counts, [0, 2, 0]);
    assert_eq!(
        (&s[..3], s[4], s[6], &s[7..]),
        (&b"AAA"[..], b'C', b'G', &b"TTT"[..])
    );
    c.ins_rate = 0.0;
    c.del_rate = 1.0;
    let (s, counts) = add_errors(
        b"AAACGTTT",
        c.sub_rate,
        c.ins_rate,
        c.del_rate,
        c.homopolymer_indel_mult,
        &mut rng,
    );
    assert_eq!((s, counts), (b"AAATTT".to_vec(), [0, 0, 2]));
}

#[test]
fn quality_mode_requires_explicit_model_and_emits_varied_qualities() {
    let base = [
        "sim",
        "--sub-rate",
        "1",
        "--ins-rate",
        "0",
        "--del-rate",
        "0",
    ];
    let mut argv = base.to_vec();
    argv.extend([
        "--quality-conditioned",
        "--error-median",
        "0.01",
        "--error-sigma",
        "0.2",
        "--local-q-sigma",
        "1.5",
        "--miscalibration",
        "1",
        "--latent-sigma",
        "0.1",
    ]);
    let parsed = Args::try_parse_from(argv);
    assert!(
        parsed.is_ok(),
        "quality-conditioned model flags must be accepted"
    );
    let mut c = parsed.unwrap().config;
    c.n_molecules = 20;
    c.family_mean = 1.0;
    let a = simulate(&c).unwrap();
    assert_eq!(a, simulate(&c).unwrap());
    for r in &a.reads {
        assert_eq!(r.sequence.len(), r.quality.len());
        assert!(r.quality.windows(2).any(|w| w[0] != w[1]));
    }
    let mut missing = base.to_vec();
    missing.push("--quality-conditioned");
    assert!(Args::try_parse_from(missing).is_err());
}

fn quality_config() -> QualityConfig {
    QualityConfig::enabled(0.01, 0.2, 1.5, 1.0, 0.15)
}

#[test]
fn quality_calibration_uses_emitted_bases_and_batch_scaling() {
    use rand::{SeedableRng, rngs::StdRng};
    let template = b"ACGT".repeat(500);
    for (weights, hp, source) in [
        ([1.0, 0.0, 0.0], 1.0, template.clone()),
        ([0.3, 0.35, 0.35], 1.0, template),
        ([0.3, 0.35, 0.35], 2.0, b"AAAACCCCGGGGTTTT".repeat(125)),
    ] {
        for factor in [1.0, 2.0] {
            let mut c = quality_config();
            c.miscalibration = Some(factor);
            let mut rng = StdRng::seed_from_u64(92);
            let (mut expected, mut actual) = (0.0, 0);
            for _ in 0..500 {
                let (seq, qual, stats) = add_quality_errors(&source, &c, weights, hp, &mut rng);
                assert_eq!(seq.len(), qual.len());
                assert_eq!(seq.len(), source.len() + stats.events[1] - stats.events[2]);
                assert!(qual.iter().all(|q| (33..=126).contains(q)));
                expected += qual
                    .iter()
                    .map(|q| 10.0_f64.powf(-f64::from(q - 33) / 10.0))
                    .sum::<f64>();
                actual += stats.events.iter().sum::<usize>();
            }
            let ratio = actual as f64 / expected;
            assert!(
                (ratio / factor - 1.0).abs() < 0.06,
                "weights={weights:?} hp={hp} factor={factor}: {ratio}"
            );
        }
    }
}

#[test]
fn quality_model_validation_and_disabled_rng_identity() {
    let good = quality_config();
    let setters: [fn(&mut QualityConfig, f64); 6] = [
        |c, v| c.error_median = Some(v),
        |c, v| c.error_sigma = Some(v),
        |c, v| c.local_q_sigma = Some(v),
        |c, v| c.miscalibration = Some(v),
        |c, v| c.latent_sigma = Some(v),
        |c, v| c.unreported_error_floor = v,
    ];
    for set in setters {
        for bad in [f64::NAN, f64::INFINITY, -0.1] {
            let mut c = good.clone();
            set(&mut c, bad);
            assert!(c.validate([1.0, 0.0, 0.0], 1.0).is_err());
        }
    }
    assert!(good.validate([0.0; 3], 1.0).is_err());
    let mut c = config();
    c.n_molecules = 10;
    c.sub_rate = 0.02;
    c.ins_rate = 0.01;
    c.del_rate = 0.03;
    let original = simulate(&c).unwrap();
    c.quality = good;
    c.quality.quality_conditioned = false;
    assert_eq!(original, simulate(&c).unwrap());
}

#[test]
fn additive_floor_events_follow_quality_plus_floor_and_latent_clipping() {
    use rand::{SeedableRng, rngs::StdRng};
    for (factor, floor, latent) in [
        (0.0, 0.02, 0.0),
        (1.7, 0.003, 0.0),
        (1.0, 0.006, 0.4),
        (4.0, 2.0, 0.0),
    ] {
        let mut c = quality_config();
        c.miscalibration = Some(factor);
        c.unreported_error_floor = floor;
        c.latent_sigma = Some(latent);
        c.error_sigma = Some(0.0);
        for nominal in [0.001, 0.02] {
            c.error_median = Some(nominal);
            let mut rng = StdRng::seed_from_u64(791);
            let (mut expected, mut events, mut positions, mut saturated) = (0.0, 0, 0, 0);
            let mut q_sum = 0.0;
            for _ in 0..400 {
                let (s, q, stats) =
                    add_quality_errors(&b"ACGT".repeat(500), &c, [1.0, 0.0, 0.0], 1.0, &mut rng);
                positions += s.len();
                saturated += stats.saturated_positions;
                events += stats.events[0];
                for &q in &q {
                    expected += (factor * 10_f64.powf(-f64::from(q - 33) / 10.0) + floor).min(1.0);
                    q_sum += f64::from(q - 33);
                }
            }
            assert!(
                (events as f64 / expected - 1.0).abs() < if latent > 0.0 { 0.10 } else { 0.04 },
                "{factor} {floor} {latent} {nominal}: events={events}, expected={expected}"
            );
            let nominal_q = -10.0 * nominal.log10();
            assert!((q_sum / positions as f64 - nominal_q).abs() < 0.1);
            if floor > 1.0 {
                assert_eq!(saturated, positions);
                assert_eq!(events, positions);
            } else {
                assert_eq!(saturated, 0);
            }
        }
    }
}

#[test]
fn latent_multiplier_applies_to_floor_before_probability_clipping() {
    use rand::{SeedableRng, rngs::StdRng};
    let mut c = quality_config();
    c.miscalibration = Some(0.0);
    c.unreported_error_floor = 0.5;
    c.latent_sigma = Some(1.0);
    let mut rng = StdRng::seed_from_u64(991);
    let source = b"ACGT".repeat(100);
    let (mut saturated, mut events) = (0, 0);
    for _ in 0..1000 {
        let (_, _, stats) = add_quality_errors(&source, &c, [1.0, 0.0, 0.0], 1.0, &mut rng);
        assert!(stats.saturated_positions == 0 || stats.saturated_positions == source.len());
        saturated += stats.saturated_positions;
        events += stats.events[0];
    }
    assert!(saturated > 0 && saturated < source.len() * 1000 / 2);
    assert!((0.30..0.48).contains(&(events as f64 / (source.len() * 1000) as f64)));
}
