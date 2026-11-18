//! Simulate ONT reads, deduplicate them, and score the result against truth.
//!
//! For local tuning sweeps, e.g.:
//!   cargo run --release --example ont_eval -- \
//!     --sub-rate <P> --ins-rate <P> --del-rate <P> --n-molecules 20000
//! (Error rates are required; choose them to suit the reads being modeled.)

#[path = "common/dedup.rs"]
mod dedup;
#[path = "common/quality_noise.rs"]
mod quality_noise;

use clap::Parser;
use nao_dedup::ont::{
    CalibrationOptions, OntDedupContext, OntMinimizerParams, estimate_noise_scale,
};
use ont_sim::{SimConfig, evaluate, simulate};
use std::time::Instant;

#[derive(Parser)]
struct Args {
    #[command(flatten)]
    sim: SimConfig,
    #[command(flatten)]
    dedup: dedup::DedupArgs,
    #[command(flatten)]
    noise: quality_noise::QualityNoiseArgs,
    /// Estimate a noise multiplier before deduplication
    #[arg(long)]
    auto_calibrate: bool,

    #[arg(long, default_value_t = OntMinimizerParams::default().kmer_len)]
    kmer_len: usize,
    #[arg(long, default_value_t = OntMinimizerParams::default().window_len)]
    window_len: usize,
    #[arg(long, default_value_t = OntMinimizerParams::default().num_windows)]
    num_windows: usize,
    /// Print metrics as one JSON object instead of text.
    #[arg(long)]
    json: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut dedup_params = args.noise.params(args.dedup.params());
    dedup_params.validate()?;
    let minimizer_params =
        OntMinimizerParams::new(args.kmer_len, args.window_len, args.num_windows)?;

    let sim = simulate(&args.sim)?;
    if args.auto_calibrate {
        let report = estimate_noise_scale(
            sim.reads
                .iter()
                .map(|r| (r.sequence.as_slice(), r.quality.as_slice())),
            &dedup_params,
            &minimizer_params,
            &CalibrationOptions::default(),
        );
        eprintln!("{report}");
        dedup_params.noise_scale = report.multiplier;
    }
    let start = Instant::now();
    let mut ctx = OntDedupContext::new(dedup_params, minimizer_params);
    for (idx, read) in sim.reads.iter().enumerate() {
        ctx.process_read_by_index(idx, &read.sequence, &read.quality);
    }
    ctx.finalize();
    let seconds = start.elapsed().as_secs_f64();

    let assignments: Vec<(usize, usize)> = (0..sim.reads.len())
        .map(|idx| (idx, ctx.get_exemplar_index(idx)))
        .collect();
    let truth: Vec<_> = sim.reads.iter().map(|r| r.truth.clone()).collect();
    let result = evaluate(&truth, &assignments)?;
    let stats = ctx.detailed_stats();
    if args.json {
        println!(
            "{{\"metrics\":{},\"seconds\":{seconds},\"unseeded_reads\":{},\"postings_visited\":{},\"candidates_compared\":{},\"unavailable_quality_reads\":{},\"assumed_rate_comparisons\":{}}}",
            result.to_json(),
            stats.unseeded_reads,
            stats.postings_visited,
            stats.candidates_compared,
            stats.unavailable_quality_reads,
            stats.assumed_rate_comparisons
        );
    } else {
        print!("{result}");
        println!(
            "seconds={seconds:.3} us_per_read={:.2} unseeded_reads={} postings_per_read={:.1} candidates_per_read={:.3} unavailable_quality_reads={} assumed_rate_comparisons={}",
            seconds * 1e6 / sim.reads.len().max(1) as f64,
            stats.unseeded_reads,
            stats.postings_visited as f64 / sim.reads.len().max(1) as f64,
            stats.candidates_compared as f64 / sim.reads.len().max(1) as f64,
            stats.unavailable_quality_reads,
            stats.assumed_rate_comparisons
        );
    }
    Ok(())
}
