//! Scale benchmark: stream many simulated ONT reads through OntDedupContext
//! and report deduplication time, generation time, and peak memory. Reads use
//! ont-sim's quality-conditioned errors and are never stored.
//!
//!   cargo run --release --example ont_scale -- --reads 10000000 \
//!     --sub-weight <W> --ins-weight <W> --del-weight <W> \
//!     --error-median <P> --error-sigma <S> --local-q-sigma <Q> \
//!     --miscalibration <M> --latent-sigma <S>
//!
//! Error and quality inputs are required runtime arguments for ont-sim's
//! `add_quality_errors` (see ont-sim/README.md), here without homopolymer
//! enrichment or an unreported error floor. Molecules are drawn uniformly over
//! the genomes, each family gets a geometric number of reads, and reads arrive
//! in random order. Use `--genomes 1 --genome-len 7000` to exercise saturated
//! seed groups. Generating quality-conditioned reads can take longer than
//! deduplicating them; the two are timed separately.

#[path = "common/dedup.rs"]
mod dedup;

use clap::Parser;
use nao_dedup::ont::{OntDedupContext, OntMinimizerParams};
use ont_sim::{QualityConfig, add_quality_errors};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use std::time::{Duration, Instant};

#[derive(Parser)]
struct Args {
    #[command(flatten)]
    dedup: dedup::DedupArgs,
    /// Total reads to generate.
    #[arg(long, default_value_t = 1_000_000)]
    reads: usize,
    /// Mean reads per molecule (geometric family sizes; 1 = no duplicates).
    #[arg(long, default_value_t = 3.0)]
    family_mean: f64,
    #[arg(long, default_value_t = 100)]
    genomes: usize,
    #[arg(long, default_value_t = 30_000)]
    genome_len: usize,
    /// Read lengths are uniform in [len_min, len_max].
    #[arg(long, default_value_t = 500)]
    len_min: usize,
    #[arg(long, default_value_t = 4500)]
    len_max: usize,
    /// Relative substitution event weight (normalized with indel weights).
    #[arg(long)]
    sub_weight: f64,
    /// Relative insertion event weight.
    #[arg(long)]
    ins_weight: f64,
    /// Relative deletion event weight.
    #[arg(long)]
    del_weight: f64,
    /// Median nominal per-read error probability (lognormal, clipped to 1).
    #[arg(long)]
    error_median: f64,
    /// Standard deviation of log nominal error probability.
    #[arg(long)]
    error_sigma: f64,
    /// Local standard deviation in Phred units, before rounding to Q0–Q93.
    #[arg(long)]
    local_q_sigma: f64,
    /// Batch multiplier on quality-derived error probabilities.
    #[arg(long)]
    miscalibration: f64,
    /// Log standard deviation of the mean-one per-read latent multiplier.
    #[arg(long)]
    latent_sigma: f64,
    #[arg(long, default_value_t = 1)]
    seed: u64,
}

fn peak_rss_mb() -> f64 {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    status
        .lines()
        .find_map(|l| l.strip_prefix("VmHWM:"))
        .and_then(|v| v.trim().trim_end_matches("kB").trim().parse::<f64>().ok())
        .map_or(f64::NAN, |kb| kb / 1024.0)
}

fn reverse_complement_in_place(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match *b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        };
    }
}

fn validate(args: &Args) -> Result<(), String> {
    if !(args.family_mean >= 1.0 && args.family_mean.is_finite()) {
        return Err("--family-mean must be at least 1".into());
    }
    if args.reads == 0 || args.genomes == 0 || args.len_min == 0 || args.len_min > args.len_max {
        return Err(
            "--reads and --genomes must be positive, and 0 < --len-min <= --len-max".into(),
        );
    }
    if args.len_min > args.genome_len {
        return Err("--len-min must not exceed --genome-len".into());
    }
    args.dedup.params().validate()?;
    Ok(())
}

fn main() {
    let args = Args::parse();
    let quality = QualityConfig::enabled(
        args.error_median,
        args.error_sigma,
        args.local_q_sigma,
        args.miscalibration,
        args.latent_sigma,
    );
    let weights = [args.sub_weight, args.ins_weight, args.del_weight];
    if let Err(e) = quality.validate(weights, 1.0).and_then(|_| validate(&args)) {
        eprintln!("Error: {e}");
        std::process::exit(2);
    }
    let mut rng = StdRng::seed_from_u64(args.seed);
    let genomes: Vec<Vec<u8>> = (0..args.genomes)
        .map(|_| {
            (0..args.genome_len)
                .map(|_| b"ACGT"[rng.gen_range(0..4)])
                .collect()
        })
        .collect();

    // Molecules (genome, start, len) and a shuffled list of which molecule
    // each read comes from.
    let mut molecules = Vec::new();
    let mut read_molecule: Vec<u32> = Vec::with_capacity(args.reads);
    while read_molecule.len() < args.reads {
        let len = rng
            .gen_range(args.len_min..=args.len_max)
            .min(args.genome_len);
        let genome = rng.gen_range(0..args.genomes) as u32;
        let start = rng.gen_range(0..=args.genome_len - len) as u32;
        let id = molecules.len() as u32;
        molecules.push((genome, start, len as u32));
        read_molecule.push(id);
        while read_molecule.len() < args.reads && !rng.gen_bool(1.0 / args.family_mean) {
            read_molecule.push(id);
        }
    }
    read_molecule.shuffle(&mut rng);
    let setup_rss = peak_rss_mb();

    let mut ctx = OntDedupContext::new(args.dedup.params(), OntMinimizerParams::default());
    let mut dedup_time = Duration::ZERO;
    let mut generation_time = Duration::ZERO;
    let started = Instant::now();
    for (idx, &m) in read_molecule.iter().enumerate() {
        let t = Instant::now();
        let (genome, start, len) = molecules[m as usize];
        let template = &genomes[genome as usize][start as usize..(start + len) as usize];
        let (mut read, mut qual, _) =
            add_quality_errors(template, &quality, weights, 1.0, &mut rng);
        if rng.gen_bool(0.5) {
            reverse_complement_in_place(&mut read);
            qual.reverse();
        }
        generation_time += t.elapsed();
        let t = Instant::now();
        ctx.process_read_by_index(idx, &read, &qual);
        dedup_time += t.elapsed();
        if (idx + 1) % 1_000_000 == 0 {
            eprintln!(
                "{} M reads, {:.1} us/read dedup, {:.1} us/read generation, peak RSS {:.0} MB",
                (idx + 1) / 1_000_000,
                dedup_time.as_secs_f64() * 1e6 / (idx + 1) as f64,
                generation_time.as_secs_f64() * 1e6 / (idx + 1) as f64,
                peak_rss_mb()
            );
        }
    }
    let t = Instant::now();
    ctx.finalize();
    dedup_time += t.elapsed();
    let stats = ctx.detailed_stats();
    println!(
        "reads={} molecules={} clusters={} dedup_s={:.1} us_per_read={:.1} generation_s={:.1} total_s={:.1} \
         setup_rss_mb={setup_rss:.0} peak_rss_mb={:.0} postings_per_read={:.0} candidates_per_read={:.2} unavailable_quality_reads={} assumed_rate_comparisons={}",
        args.reads,
        molecules.len(),
        stats.clusters,
        dedup_time.as_secs_f64(),
        dedup_time.as_secs_f64() * 1e6 / args.reads as f64,
        generation_time.as_secs_f64(),
        started.elapsed().as_secs_f64(),
        peak_rss_mb(),
        stats.postings_visited as f64 / args.reads as f64,
        stats.candidates_compared as f64 / args.reads as f64,
        stats.unavailable_quality_reads,
        stats.assumed_rate_comparisons,
    );
}
