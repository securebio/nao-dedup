//! Simulate ONT reads, deduplicate them, and score the result against truth.
//!
//! For local tuning sweeps, e.g.:
//!   cargo run --release --example ont_eval -- \
//!     --sub-rate 0.004 --ins-rate 0.004 --del-rate 0.004 --n-molecules 20000
//! (Error rates here are illustrative, not a model of any instrument.)

use clap::Parser;
use nao_dedup::ont::{OntDedupContext, OntDedupParams, OntMinimizerParams};
use ont_sim::{SimConfig, evaluate, simulate};
use std::time::Instant;

#[derive(Parser)]
struct Args {
    #[command(flatten)]
    sim: SimConfig,
    #[arg(long, default_value_t = OntDedupParams::default().max_end_offset)]
    max_end_offset: usize,
    #[arg(long, default_value_t = OntDedupParams::default().end_check_len)]
    end_check_len: usize,
    #[arg(long, default_value_t = OntDedupParams::default().max_end_errors)]
    max_end_errors: usize,
    #[arg(long, default_value_t = OntDedupParams::default().max_error_frac)]
    max_error_frac: f64,
    #[arg(long, default_value_t = OntMinimizerParams::default().kmer_len)]
    kmer_len: usize,
    #[arg(long, default_value_t = OntMinimizerParams::default().window_len)]
    window_len: usize,
    #[arg(long, default_value_t = OntMinimizerParams::default().num_windows)]
    num_windows: usize,
    /// Only match reads in the same orientation.
    #[arg(long)]
    no_rc: bool,
    /// Print metrics as one JSON object instead of text.
    #[arg(long)]
    json: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let sim = simulate(&args.sim)?;
    let dedup_params = OntDedupParams {
        max_end_offset: args.max_end_offset,
        end_check_len: args.end_check_len,
        max_end_errors: args.max_end_errors,
        max_error_frac: args.max_error_frac,
        rc_tolerant: !args.no_rc,
    };
    let minimizer_params =
        OntMinimizerParams::new(args.kmer_len, args.window_len, args.num_windows)?;

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
            "{{\"metrics\":{},\"seconds\":{seconds},\"unseeded_reads\":{},\"postings_visited\":{},\"candidates_compared\":{}}}",
            result.to_json(),
            stats.unseeded_reads,
            stats.postings_visited,
            stats.candidates_compared
        );
    } else {
        print!("{result}");
        println!(
            "seconds={seconds:.3} us_per_read={:.2} unseeded_reads={} candidates_per_read={:.3}",
            seconds * 1e6 / sim.reads.len().max(1) as f64,
            stats.unseeded_reads,
            stats.candidates_compared as f64 / sim.reads.len().max(1) as f64
        );
    }
    Ok(())
}
