use clap::Parser;
use nao_dedup::ont::{OntDedupContext, OntDedupParams, OntMinimizerParams};
use ont_sim::{SimConfig, evaluate, simulate};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn random_seq(len: usize, rng: &mut StdRng) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        })
        .collect()
}

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
    // 1000 bp with a 4% budget allows 40 edits. A and C each differ from B
    // by 25 substitutions in different regions, so A~B and B~C but not A~C.
    let mut rng = StdRng::seed_from_u64(5);
    let b = random_seq(1000, &mut rng);
    let a = with_substitutions(&b, 200..450);
    let c = with_substitutions(&b, 550..800);
    let reads = vec![
        (a.clone(), qual(1000, 10)), // leader of the first cluster
        (b.clone(), qual(1000, 30)), // joins A and becomes its exemplar
        (c.clone(), qual(1000, 20)), // fails against leader A: new cluster
    ];
    assert_eq!(run(&reads, OntDedupParams::default()), [1, 1, 2]);

    // A read matching two existing leaders joins the earlier one; the
    // clusters stay separate.
    let reads = vec![
        (a, qual(1000, 20)),
        (c, qual(1000, 20)),
        (b, qual(1000, 10)),
    ];
    assert_eq!(run(&reads, OntDedupParams::default()), [0, 1, 0]);
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
    // though it is 100 bases longer.
    let mut rng = StdRng::seed_from_u64(6);
    let long = random_seq(3100, &mut rng);
    let short = long[..3000].to_vec();
    let params = OntDedupParams {
        max_end_offset: 100,
        max_end_errors: 0,
        ..Default::default()
    };
    let mut long_qual = qual(3100, 30);
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

/// Simulate with illustrative error rates (not a model of any instrument).
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

#[test]
fn simulated_uniform_endpoints() {
    let result = evaluate_dedup(&simulated(&[]), OntDedupParams::default());
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
        ..Default::default()
    };
    let result = evaluate_dedup(&sim, strict);
    let precision = result.identifiability.precision().unwrap();
    let recall = result.identifiability.recall().unwrap();
    assert!(precision > 0.98, "precision {precision}");
    assert!(recall > 0.9, "recall {recall}");
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
