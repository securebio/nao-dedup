use super::*;
use clap::Parser;

fn config() -> SimConfig {
    #[derive(Parser)]
    struct Args {
        #[command(flatten)]
        config: SimConfig,
    }
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
        read_name: format!("r{i}"),
        molecule_id: molecule,
        genome_id: 0,
        mol_start: molecule * 10,
        mol_end: molecule * 10 + 100,
        is_rc: false,
        foldback_len: 0,
        polya_len: 0,
        n_sub: 0,
        n_ins: 0,
        n_del: 0,
        mean_q: 20.0,
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
        assert_eq!(t.read_name, format!("r{i}"));
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
fn foldbacks_tails_and_hotspot_collisions() {
    let mut c = config();
    c.n_molecules = 30;
    c.endpoints = Endpoints::Hotspot;
    c.n_hotspots = 1;
    c.hotspot_jitter = 0;
    c.foldback_prob = 1.0;
    c.foldback_min_frac = 0.5;
    c.foldback_max_frac = 0.5;
    c.polya_prob = 1.0;
    c.polya_min = 7;
    c.polya_max = 7;
    let s = simulate(&c).unwrap();
    for r in &s.reads {
        let t = &r.truth;
        let n = t.mol_end - t.mol_start;
        assert_eq!(t.foldback_len, n / 2);
        assert_eq!(
            &r.sequence[n..n + n / 2],
            simulate::reverse_complement(&r.sequence[n - n / 2..n])
        );
        assert_eq!(&r.sequence[n + n / 2..], b"AAAAAAA");
        assert_eq!(t.polya_len, 7);
    }
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
    assert_eq!((merged.merges, merged.splits), (1, 0));
    assert_eq!(merged.dedup_count_bias, Some(-0.5));
    let split = evaluate(&t, &[(0, 0), (1, 1), (2, 2), (3, 3)]).unwrap();
    assert_eq!(split.molecule.precision(), None);
    assert_eq!(split.molecule.recall(), Some(0.0));
    assert_eq!((split.merges, split.splits), (0, 2));
    let mut identical = t.clone();
    for row in &mut identical {
        row.mol_start = 0;
        row.mol_end = 100;
    }
    let e = evaluate(&identical, &[(0, 0), (1, 0), (2, 0), (3, 0)]).unwrap();
    assert_eq!(e.identifiability.precision(), Some(1.0));
    assert_eq!(e.molecule.precision(), Some(1.0 / 3.0));
    assert_eq!((e.merges, e.splits), (1, 0));
    assert_eq!((e.identifiability_merges, e.identifiability_splits), (0, 0));
    let e = evaluate(&identical, &[(0, 0), (1, 1), (2, 2), (3, 3)]).unwrap();
    assert_eq!((e.merges, e.splits), (0, 2));
    assert_eq!((e.identifiability_merges, e.identifiability_splits), (0, 1));
    // Two reads from each of two distinct intervals in one cluster: one merge
    // (counted per cluster, not per read) under both truths.
    let e = evaluate(&t, &[(0, 0), (1, 0), (2, 0), (3, 0)]).unwrap();
    assert_eq!((e.identifiability_merges, e.identifiability_splits), (1, 0));
    assert_eq!(evaluate(&[], &[]).unwrap().dedup_count_bias, None);
}

#[test]
fn stratified_recall_and_assignment_validation() {
    let mut t: Vec<_> = (0..4).map(|i| truth(i, 0)).collect();
    t[2].foldback_len = 10;
    t[3].foldback_len = 10;
    t[1].is_rc = true;
    t[3].is_rc = true;
    let a = [(3, 2), (0, 0), (2, 2), (1, 0)]; // input order must not matter
    let e = evaluate(&t, &a).unwrap();
    assert_eq!(
        e.foldback.iter().map(|s| s.true_pairs).collect::<Vec<_>>(),
        [1, 4, 1]
    );
    assert_eq!(
        e.foldback
            .iter()
            .map(|s| s.recovered_pairs)
            .collect::<Vec<_>>(),
        [1, 0, 1]
    );
    assert_eq!(
        e.strand.iter().map(|s| s.true_pairs).collect::<Vec<_>>(),
        [2, 4]
    );
    assert_eq!(
        e.strand
            .iter()
            .map(|s| s.recovered_pairs)
            .collect::<Vec<_>>(),
        [0, 2]
    );
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
fn pcr_sublineage_and_end_jitter() {
    let mut c = config();
    c.n_molecules = 20;
    c.pcr_error_prob = 1.0;
    c.rc_prob = 0.0;
    let s = simulate(&c).unwrap();
    for molecule in 0..c.n_molecules {
        let family: Vec<_> = s
            .reads
            .iter()
            .filter(|r| r.truth.molecule_id == molecule)
            .collect();
        let t = &family[0].truth;
        let template = &s.genomes[t.genome_id][t.mol_start..t.mol_end];
        let changes: Vec<_> = family
            .iter()
            .filter_map(|r| {
                let differences: Vec<_> = r
                    .sequence
                    .iter()
                    .zip(template)
                    .enumerate()
                    .filter(|(_, (a, b))| a != b)
                    .map(|(i, (a, _))| (i, *a))
                    .collect();
                assert!(differences.len() <= 1);
                differences.first().copied()
            })
            .collect();
        assert_eq!(changes.len(), family.len() / 2);
        assert!(changes.windows(2).all(|w| w[0] == w[1]));
    }
    c.pcr_error_prob = 0.0;
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
fn rejects_invalid_parameters_and_tsv() {
    let mut c = config();
    for rate in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        c.ins_rate = rate;
        assert!(simulate(&c).is_err());
    }
    c.ins_rate = 0.0;
    c.family_mean = 0.5;
    assert!(simulate(&c).is_err());
    for input in [
        "",
        "read_idx,exemplar_idx\n0,0\n",
        "read_idx\texemplar_idx\n0\tbogus\n",
        "read_idx\texemplar_idx\n0\t0\textra\n",
    ] {
        assert!(read_assignments(input.as_bytes()).is_err());
    }
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
    let (s, counts) = simulate::errors(b"AAACGTTT", &c, &mut rng);
    assert_eq!(counts, [8, 0, 0]);
    assert!(s.iter().zip(b"AAACGTTT").all(|(a, b)| a != b));
    c.sub_rate = 0.0;
    c.ins_rate = 1.0;
    c.homopolymer_indel_mult = 0.0;
    let (s, counts) = simulate::errors(b"AAACGTTT", &c, &mut rng);
    assert_eq!(counts, [0, 2, 0]);
    assert_eq!(
        (&s[..3], s[4], s[6], &s[7..]),
        (&b"AAA"[..], b'C', b'G', &b"TTT"[..])
    );
    c.ins_rate = 0.0;
    c.del_rate = 1.0;
    let (s, counts) = simulate::errors(b"AAACGTTT", &c, &mut rng);
    assert_eq!((s, counts), (b"AAATTT".to_vec(), [0, 0, 2]));
}
