//! Stream generated reads through the production deduplicator; retain only
//! genomes, molecule coordinates and shuffled molecule IDs, never all reads.
//! Error/quality parameters are required runtime inputs, not instrument defaults.
//! Use --genomes 1 --genome-len 10000 for saturated seed groups. Generation
//! can dominate runtime; the two are timed separately.
//! Peak RSS comes from /proc (Linux) and prints NaN elsewhere.
//!
//! cargo run --release --example ont_scale -- --reads 50000 \
//!   --sub-weight 0.4 --ins-weight 0.3 --del-weight 0.3 \
//!   --error-median 0.01 --error-sigma 0.2 --local-q-sigma 1.5 \
//!   --miscalibration 1 --latent-sigma 0.1

#[allow(dead_code)] // The benchmark only uses the streaming error generator.
#[path = "../tests/common/sim.rs"]
mod sim;

use clap::Parser;
use nao_dedup::ont::{OntDedupContext, OntDedupParams, OntMinimizerParams};
use rand::{Rng, SeedableRng, rngs::StdRng, seq::SliceRandom};
use sim::{QualityConfig, add_quality_errors, reverse_complement};
use std::time::{Duration, Instant};

#[derive(Parser)]
struct Args {
    #[arg(long, default_value_t = 1_000_000)]
    reads: usize,
    /// Geometric family mean (1 = no duplicates).
    #[arg(long, default_value_t = 3.0)]
    family_mean: f64,
    #[arg(long, default_value_t = 100)]
    genomes: usize,
    #[arg(long, default_value_t = 30_000)]
    genome_len: usize,
    /// Uniform molecule length bounds, clipped to genome length.
    #[arg(long, default_value_t = 500)]
    len_min: usize,
    #[arg(long, default_value_t = 4500)]
    len_max: usize,
    /// Relative substitution/insertion/deletion event weights.
    #[arg(long)]
    sub_weight: f64,
    #[arg(long)]
    ins_weight: f64,
    #[arg(long)]
    del_weight: f64,
    /// Median nominal per-read error probability (lognormal, clipped to one).
    #[arg(long)]
    error_median: f64,
    /// Log standard deviation of nominal error.
    #[arg(long)]
    error_sigma: f64,
    /// Local Phred noise standard deviation, before rounding/clamping.
    #[arg(long)]
    local_q_sigma: f64,
    /// Batch multiplier on quality-derived error probabilities.
    #[arg(long)]
    miscalibration: f64,
    /// Log standard deviation of the mean-one latent read multiplier.
    #[arg(long)]
    latent_sigma: f64,
    #[arg(long, default_value_t = 1)]
    seed: u64,
}

fn peak_rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .unwrap_or_default()
        .lines()
        .find_map(|l| l.strip_prefix("VmHWM:"))
        .and_then(|v| v.split_whitespace().next()?.parse::<f64>().ok())
        .map_or(f64::NAN, |kb| kb / 1024.0)
}

fn main() {
    let args = Args::parse();
    let quality = QualityConfig {
        error_median: args.error_median,
        error_sigma: args.error_sigma,
        local_q_sigma: args.local_q_sigma,
        miscalibration: args.miscalibration,
        latent_sigma: args.latent_sigma,
    };
    let weights = [args.sub_weight, args.ins_weight, args.del_weight];
    let valid_sizes = args.reads > 0
        && args.genomes > 0
        && args.len_min > 0
        && args.len_min <= args.len_max
        && args.len_min <= args.genome_len
        && args.family_mean.is_finite()
        && args.family_mean >= 1.0;
    if let Err(e) = quality.validate(weights) {
        eprintln!("Error: {e}");
        std::process::exit(2);
    }
    if !valid_sizes {
        eprintln!(
            "Error: require positive sizes, ordered length bounds within the genome, \
             and family mean >= 1"
        );
        std::process::exit(2);
    }
    let started = Instant::now();
    let mut rng = StdRng::seed_from_u64(args.seed);
    let genomes: Vec<Vec<u8>> = (0..args.genomes)
        .map(|_| {
            (0..args.genome_len)
                .map(|_| b"ACGT"[rng.gen_range(0..4)])
                .collect()
        })
        .collect();
    let mut molecules = Vec::new();
    let mut read_molecule = Vec::with_capacity(args.reads);
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
    let mut generation_time = started.elapsed();
    let setup_rss = peak_rss_mb();
    let mut ctx = OntDedupContext::new(OntDedupParams::default(), OntMinimizerParams::default());
    let mut dedup_time = Duration::ZERO;
    for (idx, &m) in read_molecule.iter().enumerate() {
        let t = Instant::now();
        let (genome, start, len) = molecules[m as usize];
        let template = &genomes[genome as usize][start as usize..(start + len) as usize];
        let (mut read, mut qual, _) = add_quality_errors(template, &quality, weights, &mut rng);
        if rng.gen_bool(0.5) {
            read = reverse_complement(&read);
            qual.reverse();
        }
        generation_time += t.elapsed();
        let t = Instant::now();
        ctx.process_read_by_index(idx, &read, &qual);
        dedup_time += t.elapsed();
        if (idx + 1) % 1_000_000 == 0 {
            eprintln!(
                "{} M reads, {:.1} us/read dedup, peak RSS {:.0} MB",
                (idx + 1) / 1_000_000,
                dedup_time.as_secs_f64() * 1e6 / (idx + 1) as f64,
                peak_rss_mb()
            );
        }
    }
    let t = Instant::now();
    ctx.finalize();
    dedup_time += t.elapsed();
    let stats = ctx.detailed_stats();
    println!(
        "reads={} molecules={} clusters={} dedup_s={:.2} us_per_read={:.1} generation_s={:.2} total_s={:.2} \
         setup_rss_mb={setup_rss:.1} peak_rss_mb={:.1} postings_per_read={:.1} candidates_per_read={:.2}",
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
    );
}
