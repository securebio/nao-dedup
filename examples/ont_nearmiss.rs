//! Measure pair merge probability over sequence divergence and endpoint offsets.
//!
//!   cargo run --release --example ont_nearmiss -- \
//!     --sub-rate <P> --ins-rate <P> --del-rate <P>
//!
//! Error rates are required runtime inputs, using ont-sim's error semantics.
//! Divergence=0, offset=0 is a duplicate control measuring true-pair recall.
//! Other cells measure false merges between different molecules. Each trial
//! starts with a fresh context, so merge probability is for an isolated pair.

#[path = "common/dedup.rs"]
mod dedup;
#[path = "common/quality_noise.rs"]
mod quality_noise;

use clap::{Parser, ValueEnum};
use nao_dedup::ont::{OntDedupContext, OntMinimizerParams};
use ont_sim::{QualityConfig, add_errors, add_quality_errors, reverse_complement};
use rand::{Rng, SeedableRng, rngs::StdRng};

#[derive(Parser)]
struct Args {
    #[command(flatten)]
    quality: QualityConfig,
    #[arg(long)]
    sub_rate: f64,
    #[arg(long)]
    ins_rate: f64,
    #[arg(long)]
    del_rate: f64,
    #[arg(long, default_value_t = 1.0)]
    homopolymer_indel_mult: f64,
    #[arg(long, default_value_t = 0.1)]
    variant_indel_frac: f64,
    #[arg(long, default_value_t = 500)]
    len_min: usize,
    #[arg(long, default_value_t = 5000)]
    len_max: usize,
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "0,0.0025,0.005,0.01,0.02,0.03,0.05"
    )]
    divergences: Vec<f64>,
    #[arg(long, value_delimiter = ',', default_value = "0,1,2,3,4,6")]
    offsets: Vec<usize>,
    #[arg(long, value_delimiter = ',', value_enum, default_value = "one,both")]
    modes: Vec<Mode>,
    #[arg(long, default_value_t = 1000)]
    trials: usize,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    #[command(flatten)]
    dedup: dedup::DedupArgs,
    #[command(flatten)]
    noise: quality_noise::QualityNoiseArgs,
    #[arg(long, default_value_t = OntMinimizerParams::default().kmer_len)]
    kmer_len: usize,
    #[arg(long, default_value_t = OntMinimizerParams::default().window_len)]
    window_len: usize,
    #[arg(long, default_value_t = OntMinimizerParams::default().num_windows)]
    num_windows: usize,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    One,
    Both,
}

fn mutate_source(
    source: &[u8],
    divergence: f64,
    indel_frac: f64,
    rng: &mut StdRng,
) -> (Vec<u8>, Vec<usize>) {
    let mut out = Vec::with_capacity(source.len());
    let mut map = Vec::with_capacity(source.len() + 1);
    for &base in source {
        // Boundary i is before any insertion anchored to source base i.
        map.push(out.len());
        if !rng.gen_bool(divergence) {
            out.push(base);
        } else if !rng.gen_bool(indel_frac) {
            let code = b"ACGT".iter().position(|&b| b == base).unwrap();
            out.push(b"ACGT"[(code + rng.gen_range(1..4)) % 4]);
        } else if rng.gen_bool(0.5) {
            out.push(b"ACGT"[rng.gen_range(0..4)]);
            out.push(base);
        } // Otherwise delete this source base.
    }
    map.push(out.len());
    (out, map)
}

fn molecule_interval(
    len: usize,
    flank: usize,
    offset: usize,
    mode: Mode,
    rng: &mut StdRng,
) -> (usize, usize) {
    match mode {
        Mode::Both => (
            shift_endpoint(flank, offset, rng),
            shift_endpoint(flank + len, offset, rng),
        ),
        Mode::One => {
            if rng.gen_bool(0.5) {
                (shift_endpoint(flank, offset, rng), flank + len)
            } else {
                (flank, shift_endpoint(flank + len, offset, rng))
            }
        }
    }
}

fn shift_endpoint(pos: usize, offset: usize, rng: &mut StdRng) -> usize {
    if rng.gen_bool(0.5) {
        pos + offset
    } else {
        pos - offset
    }
}

fn wilson_upper(merges: usize, trials: usize) -> f64 {
    // Upper endpoint of the two-sided 95% Wilson score interval (z = Φ⁻¹(.975)).
    let z = 1.959963984540054;
    let n = trials as f64;
    let p = merges as f64 / n;
    let z2 = z * z;
    ((p + z2 / (2.0 * n) + z * (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt()) / (1.0 + z2 / n))
        .min(1.0)
}

fn validate(args: &Args) -> Result<(), String> {
    let rates = [args.sub_rate, args.ins_rate, args.del_rate];
    args.quality.validate(rates, args.homopolymer_indel_mult)?;
    if !args.quality.quality_conditioned
        && (rates.iter().any(|r| !(0.0..1.0).contains(r)) || rates.iter().sum::<f64>() >= 1.0)
    {
        return Err(
            "sequencing error rates must be finite, in [0,1), and sum to less than 1".into(),
        );
    }
    if !(0.0..=1.0).contains(&args.variant_indel_frac)
        || !args.homopolymer_indel_mult.is_finite()
        || args.homopolymer_indel_mult < 0.0
    {
        return Err("--variant-indel-frac must be in [0,1]; --homopolymer-indel-mult must be finite and nonnegative".into());
    }
    if args.divergences.is_empty() || args.offsets.is_empty() || args.modes.is_empty() {
        return Err("--divergences, --offsets, and --modes must be nonempty".into());
    }
    if args.divergences.iter().any(|d| !(0.0..1.0).contains(d)) {
        return Err("--divergences must be finite and in [0,1)".into());
    }
    if args.trials == 0 || args.len_min == 0 || args.len_min > args.len_max {
        return Err("--trials must be positive, and 0 < --len-min <= --len-max".into());
    }
    let flank = *args.offsets.iter().max().unwrap();
    let flanks = flank.checked_mul(2).ok_or("offsets are too large")?;
    // Guarantee non-crossing endpoints for every possible sign choice.
    let inward = if args.modes.iter().any(|m| matches!(m, Mode::Both)) {
        flanks
    } else {
        flank
    };
    if args.len_min <= inward {
        return Err("--len-min must exceed the maximum inward shift (twice the largest offset for mode both)".into());
    }
    if args
        .len_max
        .checked_add(flanks)
        .is_none_or(|n| n > isize::MAX as usize / std::mem::size_of::<usize>() - 1)
    {
        return Err("lengths and offsets are too large for the coordinate map".into());
    }
    args.noise.params(args.dedup.params()).validate()?;
    OntMinimizerParams::new(args.kmer_len, args.window_len, args.num_windows)?;
    Ok(())
}

fn trial(args: &Args, mode: Mode, divergence: f64, offset: usize, rng: &mut StdRng) -> bool {
    let len = rng.gen_range(args.len_min..=args.len_max);
    let flank = *args.offsets.iter().max().unwrap();
    let source: Vec<_> = (0..len + 2 * flank)
        .map(|_| b"ACGT"[rng.gen_range(0..4)])
        .collect();
    let (variant, map) = mutate_source(&source, divergence, args.variant_indel_frac, rng);
    let (start, end) = molecule_interval(len, flank, offset, mode, rng);
    let templates = [&source[flank..flank + len], &variant[map[start]..map[end]]];
    let mut ctx = OntDedupContext::new(
        args.noise.params(args.dedup.params()),
        OntMinimizerParams::new(args.kmer_len, args.window_len, args.num_windows).unwrap(),
    );
    // Under --no-rc both reads share one strand, so the duplicate control
    // still measures recall; otherwise each read's strand is independent.
    let shared_strand = rng.gen_bool(0.5);
    for (idx, template) in templates.into_iter().enumerate() {
        let (mut read, mut quality) = if args.quality.quality_conditioned {
            let (read, quality, _) = add_quality_errors(
                template,
                &args.quality,
                [args.sub_rate, args.ins_rate, args.del_rate],
                args.homopolymer_indel_mult,
                rng,
            );
            (read, quality)
        } else {
            let (read, _) = add_errors(
                template,
                args.sub_rate,
                args.ins_rate,
                args.del_rate,
                args.homopolymer_indel_mult,
                rng,
            );
            let quality = vec![b'5'; read.len()];
            (read, quality)
        };
        let rc = if args.dedup.no_rc {
            shared_strand
        } else {
            rng.gen_bool(0.5)
        };
        if rc {
            read = reverse_complement(&read);
            quality.reverse();
        }
        ctx.process_read_by_index(idx, &read, &quality);
    }
    ctx.finalize();
    ctx.get_exemplar_index(0) == ctx.get_exemplar_index(1)
}

fn main() {
    let args = Args::parse();
    if let Err(e) = validate(&args) {
        eprintln!("Error: {e}");
        std::process::exit(2);
    }
    let mut rng = StdRng::seed_from_u64(args.seed);
    println!("mode\tdivergence\toffset\ttrials\tmerges\tmerge_rate\tmerge_rate_upper95");
    for &mode in &args.modes {
        let name = match mode {
            Mode::One => "one",
            Mode::Both => "both",
        };
        for &divergence in &args.divergences {
            for &offset in &args.offsets {
                let merges = (0..args.trials)
                    .filter(|_| trial(&args, mode, divergence, offset, &mut rng))
                    .count();
                println!(
                    "{name}\t{divergence}\t{offset}\t{}\t{merges}\t{}\t{}",
                    args.trials,
                    merges as f64 / args.trials as f64,
                    wilson_upper(merges, args.trials)
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn variants_map_half_open_source_boundaries() {
        let source = b"ACGT".repeat(32);
        let mut rng = StdRng::seed_from_u64(7);
        let (same, map) = mutate_source(&source, 0.0, 1.0, &mut rng);
        assert_eq!(same, source);
        assert_eq!(map, (0..=source.len()).collect::<Vec<_>>());
        let (substituted, map) = mutate_source(&source, 1.0, 0.0, &mut rng);
        assert!(substituted.iter().zip(&source).all(|(a, b)| a != b));
        assert_eq!(map, (0..=source.len()).collect::<Vec<_>>());
        let (indels, map) = mutate_source(&source, 1.0, 1.0, &mut rng);
        assert_eq!(map.len(), source.len() + 1);
        assert_eq!((map[0], map[source.len()]), (0, indels.len()));
        let mut sizes = HashSet::new();
        for (i, &base) in source.iter().enumerate() {
            let emitted = &indels[map[i]..map[i + 1]];
            sizes.insert(emitted.len());
            match emitted {
                [] => (), // deleted source base
                [inserted, retained] => {
                    assert!(b"ACGT".contains(inserted));
                    assert_eq!(*retained, base); // insertion before this source base
                }
                _ => panic!("expected a single insertion or deletion"),
            }
        }
        assert_eq!(sizes, HashSet::from([0, 2]));
    }

    #[test]
    fn offsets_have_requested_magnitude_and_shifted_end_count() {
        let mut rng = StdRng::seed_from_u64(1);
        for mode in [Mode::One, Mode::Both] {
            let mut seen = HashSet::new();
            for _ in 0..100 {
                let (start, end) = molecule_interval(100, 6, 6, mode, &mut rng);
                let shifts = (start as isize - 6, end as isize - 106);
                match mode {
                    Mode::One => assert!(matches!(shifts, (-6 | 6, 0) | (0, -6 | 6))),
                    Mode::Both => assert!(matches!(shifts, (-6 | 6, -6 | 6))),
                }
                seen.insert(shifts);
                assert_eq!(molecule_interval(100, 6, 0, mode, &mut rng), (6, 106));
            }
            assert_eq!(seen.len(), 4);
        }
    }

    #[test]
    fn wilson_upper_includes_uncertainty_at_zero_merges() {
        assert!((wilson_upper(0, 100) - 0.0369934982).abs() < 1e-9);
        assert!((wilson_upper(50, 100) - 0.5961684696).abs() < 1e-9);
        assert!((wilson_upper(100, 100) - 1.0).abs() < 1e-12);
    }

    fn args() -> Args {
        Args::parse_from([
            "nearmiss",
            "--sub-rate",
            "0",
            "--ins-rate",
            "0",
            "--del-rate",
            "0",
        ])
    }

    fn with_flag(flag: &str, value: &str) -> Args {
        let mut argv = vec!["nearmiss".to_string()];
        for rate in ["--sub-rate", "--ins-rate", "--del-rate"] {
            if rate != flag {
                argv.push(format!("{rate}=0"));
            }
        }
        argv.push(format!("{flag}={value}"));
        Args::try_parse_from(argv).unwrap()
    }

    #[test]
    fn validates_rates_grids_and_endpoint_ranges() {
        assert!(validate(&args()).is_ok());
        for flag in [
            "--sub-rate",
            "--ins-rate",
            "--del-rate",
            "--max-error-frac",
            "--homopolymer-indel-mult",
            "--variant-indel-frac",
            "--divergences",
        ] {
            let a = with_flag(flag, "NaN");
            assert!(validate(&a).is_err(), "{flag} accepted NaN");
        }
        for (flag, value) in [
            ("--sub-rate", "-0.1"),
            ("--ins-rate", "1"),
            ("--del-rate", "inf"),
            ("--variant-indel-frac", "1.1"),
            ("--homopolymer-indel-mult", "-1"),
            ("--max-error-frac", "1.1"),
            ("--divergences", "1"),
            ("--len-min", "0"),
            ("--len-min", "6000"),
            ("--len-min", "12"),
            ("--trials", "0"),
            ("--kmer-len", "1"),
            ("--num-windows", "0"),
            ("--len-max", &usize::MAX.to_string()),
        ] {
            let a = with_flag(flag, value);
            assert!(validate(&a).is_err(), "{flag}={value} accepted");
        }
        let mut a = args();
        a.sub_rate = 0.5;
        a.ins_rate = 0.5;
        assert!(validate(&a).is_err());
        a.quality = QualityConfig::enabled(0.01, 0.2, 1.5, 1.0, 0.1);
        a.sub_rate = 4.0;
        a.ins_rate = 3.0;
        a.del_rate = 3.0;
        assert!(validate(&a).is_ok());
        let mut a = args();
        a.divergences.clear();
        assert!(validate(&a).is_err());
        let mut a = args();
        a.offsets.clear();
        assert!(validate(&a).is_err());
        let mut a = args();
        a.modes.clear();
        assert!(validate(&a).is_err());
    }

    #[test]
    fn identical_template_control_merges_with_and_without_rc() {
        let mut a = args();
        let mut rng = StdRng::seed_from_u64(a.seed);
        for _ in 0..20 {
            assert!(trial(&a, Mode::Both, 0.0, 0, &mut rng));
        }
        a.dedup.no_rc = true;
        for _ in 0..20 {
            assert!(trial(&a, Mode::One, 0.0, 0, &mut rng));
        }
    }

    #[test]
    fn strict_dedup_distinguishes_sequence_and_endpoint_near_misses() {
        let mut a = args();
        a.dedup.max_end_offset = 0;
        a.dedup.max_end_errors = 0;
        a.dedup.max_error_frac = 0.0;
        let mut rng = StdRng::seed_from_u64(a.seed);
        for _ in 0..20 {
            assert!(trial(&a, Mode::One, 0.0, 0, &mut rng));
            assert!(!trial(&a, Mode::One, 0.0, 2, &mut rng));
            assert!(!trial(&a, Mode::Both, 0.0, 2, &mut rng));
            a.variant_indel_frac = 0.0;
            assert!(!trial(&a, Mode::Both, 0.5, 0, &mut rng));
            a.variant_indel_frac = 1.0;
            assert!(!trial(&a, Mode::Both, 0.5, 0, &mut rng));
        }
    }
}
