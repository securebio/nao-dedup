//! Scale benchmark: stream many simulated ONT reads through OntDedupContext
//! and report throughput and peak memory. Unlike `ont-sim`, reads are
//! generated on the fly and never stored, so tens of millions are feasible.
//!
//!   cargo run --release --example ont_scale -- --reads 10000000 \
//!     --sub-rate <P> --ins-rate <P> --del-rate <P>
//!
//! Constant qualities exercise the assumed-error-rate path.
//!
//! Error rates are required and illustrative. Molecules are drawn uniformly
//! over the genomes, each family gets a geometric number of reads, and reads
//! arrive in random order. Errors are placed by geometric skipping, so
//! generation is cheap relative to deduplication. The rates are per-event
//! probabilities at each step through the template; an insertion doesn't
//! advance the template, so insertions per template base come out slightly
//! higher than `--ins-rate` (rate / (1 - rate)).

#[path = "common/dedup.rs"]
mod dedup;

use clap::Parser;
use nao_dedup::ont::{OntDedupContext, OntMinimizerParams};
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
    #[arg(long)]
    sub_rate: f64,
    #[arg(long)]
    ins_rate: f64,
    #[arg(long)]
    del_rate: f64,
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

/// Copy `template` into `out` with sequencing errors, placing each kind of
/// error by geometric skipping over template positions.
fn sequence(template: &[u8], args: &Args, rng: &mut StdRng, out: &mut Vec<u8>) {
    let total = args.sub_rate + args.ins_rate + args.del_rate;
    out.clear();
    let mut pos = 0;
    while pos < template.len() {
        // Bases until the next error of any kind.
        let skip = if total > 0.0 {
            let u: f64 = rng.r#gen();
            ((1.0 - u).ln() / (1.0 - total).ln()) as usize
        } else {
            usize::MAX
        };
        let run_end = pos.saturating_add(skip).min(template.len());
        out.extend_from_slice(&template[pos..run_end]);
        pos = run_end;
        if pos >= template.len() {
            break;
        }
        let kind = rng.r#gen::<f64>() * total;
        if kind < args.sub_rate {
            let code = b"ACGT".iter().position(|&b| b == template[pos]).unwrap();
            out.push(b"ACGT"[(code + rng.gen_range(1..4)) % 4]);
            pos += 1;
        } else if kind < args.sub_rate + args.ins_rate {
            out.push(b"ACGT"[rng.gen_range(0..4)]);
        } else {
            pos += 1; // deletion
        }
    }
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
    let rates = [args.sub_rate, args.ins_rate, args.del_rate];
    if rates.iter().any(|r| !r.is_finite() || *r < 0.0) || rates.iter().sum::<f64>() >= 1.0 {
        return Err("error rates must be finite, nonnegative, and sum to less than 1".into());
    }
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
    if let Err(e) = validate(&args) {
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
    let (mut read, mut qual) = (Vec::new(), Vec::new());
    let mut dedup_time = Duration::ZERO;
    let started = Instant::now();
    for (idx, &m) in read_molecule.iter().enumerate() {
        let (genome, start, len) = molecules[m as usize];
        let template = &genomes[genome as usize][start as usize..(start + len) as usize];
        sequence(template, &args, &mut rng, &mut read);
        if rng.gen_bool(0.5) {
            reverse_complement_in_place(&mut read);
        }
        qual.resize(read.len(), b'5');
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
        "reads={} molecules={} clusters={} dedup_s={:.1} us_per_read={:.1} total_s={:.1} \
         setup_rss_mb={setup_rss:.0} peak_rss_mb={:.0} postings_per_read={:.0} candidates_per_read={:.2} unavailable_quality_reads={} assumed_rate_comparisons={}",
        args.reads,
        molecules.len(),
        stats.clusters,
        dedup_time.as_secs_f64(),
        dedup_time.as_secs_f64() * 1e6 / args.reads as f64,
        started.elapsed().as_secs_f64(),
        peak_rss_mb(),
        stats.postings_visited as f64 / args.reads as f64,
        stats.candidates_compared as f64 / args.reads as f64,
        stats.unavailable_quality_reads,
        stats.assumed_rate_comparisons,
    );
}
