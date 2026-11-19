mod common;
use clap::Parser;
use common::random_seq;
use nao_dedup::ont::{OntDedupContext, OntDedupParams, OntMinimizerParams};
use ont_sim::{SimConfig, evaluate, reverse_complement, simulate};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn qual(len: usize, q: u8) -> Vec<u8> {
    vec![q + 33; len]
}

/// Deduplicate (sequence, quality) reads; returns each read's exemplar.
fn run(reads: &[(Vec<u8>, Vec<u8>)], params: OntDedupParams) -> Vec<usize> {
    let mut ctx = OntDedupContext::new(params, OntMinimizerParams::default());
    for (idx, (seq, q)) in reads.iter().enumerate() {
        ctx.process_read_by_index(idx, seq, q);
    }
    ctx.finalize();
    (0..reads.len())
        .map(|i| ctx.get_exemplar_index(i))
        .collect()
}

#[test]
fn duplicates_cluster_and_best_quality_is_exemplar() {
    let mut rng = StdRng::seed_from_u64(1);
    let seq = random_seq(1500, &mut rng);
    let other = random_seq(1500, &mut rng);
    let reads = vec![
        (seq.clone(), qual(1500, 10)),
        (other.clone(), qual(1500, 20)),
        (seq.clone(), qual(1500, 30)),
        (seq, qual(1500, 15)),
    ];
    assert_eq!(run(&reads, OntDedupParams::default()), [2, 1, 2, 2]);
}

#[test]
fn reverse_complement_matches_only_when_tolerated() {
    let mut rng = StdRng::seed_from_u64(2);
    let seq = random_seq(1200, &mut rng);
    let reads = vec![
        (seq.clone(), qual(1200, 20)),
        (reverse_complement(&seq), qual(1200, 20)),
    ];
    assert_eq!(run(&reads, OntDedupParams::default()), [0, 0]);
    let strict = OntDedupParams {
        rc_tolerant: false,
        ..Default::default()
    };
    assert_eq!(run(&reads, strict), [0, 1]);
}

#[test]
fn different_endpoints_do_not_cluster() {
    let mut rng = StdRng::seed_from_u64(3);
    let genome = random_seq(3000, &mut rng);
    let reads = vec![
        (genome[100..2100].to_vec(), qual(2000, 20)),
        (genome[120..2100].to_vec(), qual(1980, 20)), // start differs by 20
        (genome[100..2080].to_vec(), qual(1980, 20)), // end differs by 20
        (genome[100..2100].to_vec(), qual(2000, 20)), // duplicate of read 0
    ];
    assert_eq!(run(&reads, OntDedupParams::default()), [0, 1, 2, 0]);
}

#[test]
fn case_is_ignored() {
    let mut rng = StdRng::seed_from_u64(4);
    let seq = random_seq(800, &mut rng);
    let reads = vec![
        (seq.clone(), qual(800, 20)),
        (seq.to_ascii_lowercase(), qual(800, 20)),
    ];
    assert_eq!(run(&reads, OntDedupParams::default()), [0, 0]);
}

/// `seq` with substitutions at every 10th position in `range`.
fn with_substitutions(seq: &[u8], range: std::ops::Range<usize>) -> Vec<u8> {
    let mut out = seq.to_vec();
    for pos in range.step_by(10) {
        out[pos] = if out[pos] == b'A' { b'C' } else { b'A' };
    }
    out
}

#[test]
fn members_are_judged_against_the_fixed_leader() {
    // At 1000 bp with constant qualities (assumed 1% error per read) and 1%
    // protected divergence, the cutoff lies between 10 and 20 edits. A and C
    // each differ from B by 10 edits, but A and C differ by 20.
    let params = OntDedupParams {
        protect_divergence: 0.01,
        ..Default::default()
    };
    let mut rng = StdRng::seed_from_u64(5);
    let b = random_seq(1000, &mut rng);
    let a = with_substitutions(&b, 200..300);
    let c = with_substitutions(&b, 550..650);
    let reads = vec![
        (a.clone(), qual(1000, 10)), // leader of the first cluster
        (b.clone(), qual(1000, 30)), // joins A and becomes its exemplar
        (c.clone(), qual(1000, 20)), // fails against leader A: new cluster
    ];
    assert_eq!(run(&reads, params.clone()), [1, 1, 2]);

    // A read matching two existing leaders joins the earlier one; the
    // clusters stay separate.
    let reads = vec![
        (a, qual(1000, 20)),
        (c, qual(1000, 20)),
        (b, qual(1000, 10)),
    ];
    assert_eq!(run(&reads, params.clone()), [0, 1, 0]);
}

#[test]
fn indel_inside_seeded_region_does_not_hide_duplicate() {
    // A 9-base insertion just past the endpoint window shifts later seed
    // offsets; the candidate filter must still let the pair through.
    let mut rng = StdRng::seed_from_u64(8);
    let seq = random_seq(1000, &mut rng);
    let mut read = seq.clone();
    read.splice(100..100, b"GATTACAGA".iter().copied());
    for other in [read.clone(), reverse_complement(&read)] {
        let reads = vec![(seq.clone(), qual(1000, 20)), (other, qual(1009, 20))];
        let exemplars = run(&reads, OntDedupParams::default());
        assert_eq!(exemplars[0], exemplars[1]);
    }
}

#[test]
fn displaced_leaders_are_skipped_without_aligning() {
    // Thirty distinct molecules, each shifted 12 bases along the genome
    // from the previous one. They share seeds, but at offsets too far apart
    // for duplicates, so the offset filter should skip nearly all of them.
    let mut rng = StdRng::seed_from_u64(9);
    let genome = random_seq(3000, &mut rng);
    let mut ctx = OntDedupContext::new(OntDedupParams::default(), OntMinimizerParams::default());
    for j in 0..30 {
        let start = 100 + 12 * j;
        ctx.process_read_by_index(j, &genome[start..start + 2000], &qual(2000, 20));
    }
    ctx.finalize();
    assert_eq!(ctx.stats(), (30, 30));
    let compared = ctx.detailed_stats().candidates_compared;
    // Without the filter this is over 200.
    assert!(compared < 30, "{compared} candidates compared");
}

#[test]
fn quality_outranks_length() {
    // The longer read has slightly lower mean quality, so it loses even
    // though it is 10 bases longer.
    let mut rng = StdRng::seed_from_u64(6);
    let long = random_seq(3010, &mut rng);
    let short = long[..3000].to_vec();
    let params = OntDedupParams {
        max_end_offset: 10,
        ..cap_only()
    };
    let mut long_qual = qual(3010, 30);
    long_qual[0] = b'!'; // mean quality just under 30
    let reads = vec![(short, qual(3000, 30)), (long, long_qual)];
    let exemplars = run(&reads, params);
    assert_eq!(exemplars[0], exemplars[1]);
    assert_eq!(exemplars[0], 0);
}

#[test]
fn finalize_is_idempotent() {
    let mut rng = StdRng::seed_from_u64(7);
    let seq = random_seq(800, &mut rng);
    let mut ctx = OntDedupContext::new(OntDedupParams::default(), OntMinimizerParams::default());
    ctx.process_read_by_index(0, &seq, &qual(800, 10));
    ctx.process_read_by_index(1, &seq, &qual(800, 30));
    ctx.finalize();
    ctx.finalize();
    assert_eq!(
        (ctx.get_exemplar_index(0), ctx.get_exemplar_index(1)),
        (1, 1)
    );
}

#[test]
fn unseedable_reads_are_singletons() {
    let reads = [
        (b"ACGT".to_vec(), qual(4, 20)),
        (b"ACGT".to_vec(), qual(4, 20)),
        (Vec::new(), Vec::new()),
        (vec![b'N'; 500], qual(500, 20)),
    ];
    let mut ctx = OntDedupContext::new(OntDedupParams::default(), OntMinimizerParams::default());
    for (idx, (seq, q)) in reads.iter().enumerate() {
        ctx.process_read_by_index(idx, seq, q);
    }
    assert_eq!(ctx.detailed_stats().unseeded_reads, 4);
    ctx.finalize();
    assert_eq!(ctx.stats(), (4, 4));
    assert_eq!(ctx.get_exemplar_indices().len(), 4);
}

#[test]
#[should_panic(expected = "get_exemplar_index must be called after finalize()")]
fn exemplar_query_before_finalize_panics() {
    let ctx = OntDedupContext::new(OntDedupParams::default(), OntMinimizerParams::default());
    ctx.get_exemplar_index(0);
}

/// Simulate reads with substitutions and indels.
fn simulated(extra: &[&str]) -> ont_sim::Simulation {
    #[derive(Parser)]
    struct Args {
        #[command(flatten)]
        config: SimConfig,
    }
    let mut args = vec![
        "sim",
        "--seed",
        "7",
        "--n-molecules",
        "1500",
        "--sub-rate",
        "0.004",
        "--ins-rate",
        "0.004",
        "--del-rate",
        "0.004",
    ];
    args.extend_from_slice(extra);
    simulate(&Args::parse_from(args).config).unwrap()
}

fn evaluate_dedup(sim: &ont_sim::Simulation, params: OntDedupParams) -> ont_sim::Evaluation {
    let reads: Vec<_> = sim
        .reads
        .iter()
        .map(|r| (r.sequence.clone(), r.quality.clone()))
        .collect();
    let exemplars = run(&reads, params);
    let truth: Vec<_> = sim.reads.iter().map(|r| r.truth.clone()).collect();
    let assignments: Vec<_> = exemplars.into_iter().enumerate().collect();
    evaluate(&truth, &assignments).unwrap()
}

// Isolate seeds and endpoint matching with a cap-only whole-read budget.
fn cap_only() -> OntDedupParams {
    OntDedupParams {
        protect_divergence: 1.0,
        ..Default::default()
    }
}

#[test]
fn simulated_uniform_endpoints() {
    let result = evaluate_dedup(&simulated(&[]), cap_only());
    let (precision, recall) = (
        result.molecule.precision().unwrap(),
        result.molecule.recall().unwrap(),
    );
    assert!(precision > 0.99, "precision {precision}");
    assert!(recall > 0.93, "recall {recall}");
}

#[test]
fn simulated_nearby_endpoints() {
    // Many distinct molecules whose endpoints lie within 10 bp of each other.
    let sim = simulated(&[
        "--endpoints",
        "hotspot",
        "--n-hotspots",
        "60",
        "--hotspot-jitter",
        "10",
    ]);
    let strict = OntDedupParams {
        max_end_offset: 0,
        ..cap_only()
    };
    let result = evaluate_dedup(&sim, strict);
    let precision = result.identifiability.precision().unwrap();
    let recall = result.identifiability.recall().unwrap();
    assert!(precision > 0.98, "precision {precision}");
    assert!(recall > 0.9, "recall {recall}");
}

#[test]
fn default_params_recall_on_quality_conditioned_reads() {
    let sim = simulated(&[
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
        "--len-min",
        "2000",
        "--len-max",
        "3000",
        "--genome-len",
        "1000000",
    ]);
    let result = evaluate_dedup(&sim, OntDedupParams::default());
    let recall = result.molecule.recall().unwrap();
    assert!(recall > 0.95, "recall {recall}");
}

#[test]
fn deterministic() {
    let sim = simulated(&[]);
    let reads: Vec<_> = sim
        .reads
        .iter()
        .map(|r| (r.sequence.clone(), r.quality.clone()))
        .collect();
    assert_eq!(
        run(&reads, OntDedupParams::default()),
        run(&reads, OntDedupParams::default())
    );
}

#[test]
fn unavailable_qualities_are_counted_and_use_protected_assumed_noise() {
    let mut rng = StdRng::seed_from_u64(43);
    let seq = random_seq(1000, &mut rng);
    let mut ctx = OntDedupContext::new(OntDedupParams::default(), OntMinimizerParams::default());
    ctx.process_read_by_index(0, &seq, &[]);
    ctx.process_read_by_index(1, &seq, &qual(1000, 20));
    let mut varied = qual(1000, 20);
    varied[0] += 1;
    ctx.process_read_by_index(2, &seq, &varied);
    assert_eq!(ctx.detailed_stats().unavailable_quality_reads, 2);
    assert_eq!(ctx.detailed_stats().assumed_rate_comparisons, 2);
    ctx.finalize();
    assert_eq!(ctx.stats(), (3, 1));
    // At 50 bp even an identical pair lacks enough divergence evidence.
    let mut ctx = OntDedupContext::new(OntDedupParams::default(), OntMinimizerParams::default());
    ctx.process_read_by_index(0, &seq[..50], &[]);
    ctx.process_read_by_index(1, &seq[..50], &[]);
    ctx.finalize();
    assert_eq!(ctx.stats(), (2, 2));
}

#[test]
fn leader_quality_summary_survives_exemplar_replacement() {
    let mut rng = StdRng::seed_from_u64(44);
    let seq = random_seq(2000, &mut rng);
    let altered = with_substitutions(&seq, 500..700);
    let mut moderate = qual(2000, 20);
    moderate[0] += 1;
    let mut clean = qual(2000, 40);
    clean[0] += 1;
    // Clean+clean would reject twenty edits; leader+clean allows them.
    let params = OntDedupParams {
        protect_divergence: 0.01,
        ..Default::default()
    };
    let mut ctx = OntDedupContext::new(params.clone(), OntMinimizerParams::default());
    ctx.process_read_by_index(0, &seq, &moderate);
    ctx.process_read_by_index(1, &seq, &clean);
    ctx.process_read_by_index(2, &altered, &clean);
    assert_eq!(ctx.detailed_stats().unavailable_quality_reads, 0);
    assert_eq!(ctx.detailed_stats().assumed_rate_comparisons, 0);
    ctx.finalize();
    assert_eq!(ctx.stats(), (3, 1));
    assert_eq!(ctx.get_exemplar_index(2), 1);
    assert_ne!(
        run(&[(seq, clean.clone()), (altered, clean)], params),
        [0, 0]
    );
}

#[test]
fn quality_conditioned_duplicate_recall_and_protected_near_misses() {
    use ont_sim::{QualityConfig, add_quality_errors};
    let model = |median| QualityConfig::enabled(median, 0.15, 1.5, 1.0, 0.1);
    let mut rng = StdRng::seed_from_u64(9281);
    let params = OntDedupParams::default();
    let trials = 150;
    let mut duplicates = 0;
    // Assess both fixed variant counts and independent random substitutions,
    // unequal accuracies, arrival orders, and opposite-strand comparisons.
    for fixed in [false, true] {
        let mut false_merges = 0;
        for trial in 0..trials {
            let source = random_seq(3000, &mut rng);
            let mut variant = source.clone();
            // Strains 1.5x the protected divergence apart: the nominal 5% applies
            // exactly at the protected divergence, and merges fall off quickly
            // beyond it.
            let positions: Vec<usize> = if fixed {
                use rand::seq::index::sample;
                sample(
                    &mut rng,
                    source.len(),
                    (1.5 * params.protect_divergence * source.len() as f64).ceil() as usize,
                )
                .into_vec()
            } else {
                (0..source.len())
                    .filter(|_| rng.gen_bool(1.5 * params.protect_divergence))
                    .collect()
            };
            for pos in positions {
                variant[pos] = if variant[pos] == b'A' { b'C' } else { b'A' };
            }
            let (a, qa, _) =
                add_quality_errors(&source, &model(0.001), [0.4, 0.3, 0.3], 1.0, &mut rng);
            let (mut b, mut qb, _) =
                add_quality_errors(&source, &model(0.002), [0.4, 0.3, 0.3], 1.0, &mut rng);
            let (mut near, mut qn, _) =
                add_quality_errors(&variant, &model(0.002), [0.4, 0.3, 0.3], 1.0, &mut rng);
            if trial % 2 == 0 {
                b = reverse_complement(&b);
                qb.reverse();
                near = reverse_complement(&near);
                qn.reverse();
            }
            let mut dup_reads = vec![(a.clone(), qa.clone()), (b, qb)];
            let mut near_reads = vec![(a, qa), (near, qn)];
            if trial % 3 == 0 {
                dup_reads.reverse();
                near_reads.reverse();
            }
            let assignments = run(&dup_reads, params.clone());
            duplicates += usize::from(assignments[0] == assignments[1]);
            let assignments = run(&near_reads, params.clone());
            false_merges += usize::from(assignments[0] == assignments[1]);
        }
        // With n=150, <=2 merges has one-sided 95% Wilson bound below 5%.
        assert!(
            false_merges <= 2,
            "fixed={fixed}: {false_merges}/{trials} false merges"
        );
    }
    assert!(
        duplicates as f64 / (2 * trials) as f64 > 0.9,
        "{duplicates}/{} duplicates merged",
        2 * trials
    );
}

#[test]
fn default_params_protect_related_strains_with_good_duplicate_recall() {
    use ont_sim::{QualityConfig, add_quality_errors};
    use rand::seq::index::sample;

    let params = OntDedupParams::default();
    let model = QualityConfig::enabled(0.01, 0.2, 1.5, 1.0, 0.1);
    let mut rng = StdRng::seed_from_u64(5801);
    let trials = 150;
    let (mut duplicates, mut false_merges) = (0, 0);
    for trial in 0..trials {
        let source = random_seq(3000, &mut rng);
        let mut strain = source.clone();
        let changes = (1.5 * params.protect_divergence * source.len() as f64).ceil() as usize;
        for pos in sample(&mut rng, source.len(), changes) {
            strain[pos] = if strain[pos] == b'A' { b'C' } else { b'A' };
        }
        let mut reads = Vec::new();
        for template in [&source, &source, &strain] {
            let (mut seq, mut q, _) =
                add_quality_errors(template, &model, [0.4, 0.3, 0.3], 1.0, &mut rng);
            if rng.gen_bool(0.5) {
                seq = reverse_complement(&seq);
                q.reverse();
            }
            reads.push((seq, q));
        }
        // Each pair starts a fresh deduplicator, including seed retrieval.
        for (other, count) in [(1, &mut duplicates), (2, &mut false_merges)] {
            let mut pair = vec![reads[0].clone(), reads[other].clone()];
            if trial % 2 == 0 {
                pair.reverse();
            }
            let assignments = run(&pair, params.clone());
            *count += usize::from(assignments[0] == assignments[1]);
        }
    }
    assert!(duplicates >= 135, "{duplicates}/{trials} duplicates merged");
    // Allow sampling margin without claiming the model's 5% tail is a
    // pipeline-wide bound, or testing only specially calibrated parameters.
    assert!(
        false_merges <= 7,
        "{false_merges}/{trials} strain pairs merged"
    );
}

#[test]
fn default_params_preserve_rare_strains_in_dense_mixed_clusters() {
    use ont_sim::{QualityConfig, add_quality_errors};
    use rand::seq::index::sample;

    struct Read {
        seq: Vec<u8>,
        qual: Vec<u8>,
        molecule: usize,
        strain_b: bool,
    }

    let params = OntDedupParams::default();
    let model = QualityConfig::enabled(0.01, 0.2, 1.5, 1.0, 0.1);
    let mut rng = StdRng::seed_from_u64(5802);
    let molecules = 12;
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for molecule in 0..molecules {
        let source = random_seq(3000, &mut rng);
        let mut strain = source.clone();
        let changes = (1.5 * params.protect_divergence * source.len() as f64).ceil() as usize;
        for pos in sample(&mut rng, source.len(), changes) {
            strain[pos] = if strain[pos] == b'A' { b'C' } else { b'A' };
        }
        // One A and one B molecule have coincident endpoints at each site;
        // PCR abundance differs sixteenfold. Every copy gets fresh errors.
        for (strain_b, template, copies) in [(false, &source, 64), (true, &strain, 4)] {
            for _ in 0..copies {
                let (mut seq, mut qual, _) =
                    add_quality_errors(template, &model, [0.4, 0.3, 0.3], 1.0, &mut rng);
                if rng.gen_bool(0.5) {
                    seq = reverse_complement(&seq);
                    qual.reverse();
                }
                let read = Read {
                    seq,
                    qual,
                    molecule,
                    strain_b,
                };
                if strain_b {
                    b.push(read);
                } else {
                    a.push(read);
                }
            }
        }
    }
    let interleaved: Vec<_> = a
        .chunks(16)
        .zip(&b)
        .flat_map(|(chunk, read)| chunk.iter().chain(std::iter::once(read)))
        .collect();
    for (name, order) in [
        ("B first", b.iter().chain(&a).collect::<Vec<_>>()),
        ("interleaved", interleaved),
        ("B last", a.iter().chain(&b).collect()),
    ] {
        assert_eq!(order.len(), a.len() + b.len());
        let mut ctx = OntDedupContext::new(params.clone(), OntMinimizerParams::default());
        for (idx, read) in order.iter().enumerate() {
            ctx.process_read_by_index(idx, &read.seq, &read.qual);
        }
        ctx.finalize();
        let a_clusters: std::collections::HashSet<_> = order
            .iter()
            .enumerate()
            .filter(|(_, read)| !read.strain_b)
            .map(|(idx, _)| ctx.get_exemplar_index(idx))
            .collect();
        let mut represented = vec![false; molecules];
        let mut absorbed = 0;
        let mut b_duplicates = 0;
        for (idx, read) in order.iter().enumerate().filter(|(_, r)| r.strain_b) {
            let exemplar = ctx.get_exemplar_index(idx);
            let own_cluster = order[exemplar].strain_b && !a_clusters.contains(&exemplar);
            absorbed += usize::from(!own_cluster);
            if own_cluster {
                assert_eq!(order[exemplar].molecule, read.molecule, "{name}");
                represented[read.molecule] = true;
                b_duplicates += usize::from(exemplar != idx);
            }
        }
        assert!(
            absorbed <= 2,
            "{name}: {absorbed}/{} B reads absorbed",
            b.len()
        );
        assert!(represented.iter().all(|&r| r), "{name}: B molecule lost");
        // Protection must retain real clusters rather than split every read.
        assert!(
            b_duplicates >= 30,
            "{name}: only {b_duplicates} B duplicates merged"
        );
    }
}
