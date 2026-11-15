use clap::Parser;
use ont_sim::{SimConfig, simulate, write_outputs};
use std::{error::Error, path::PathBuf};

#[derive(Parser)]
#[command(about = "Generate synthetic ONT reads and PCR-molecule truth")]
struct Args {
    #[command(flatten)]
    config: SimConfig,
    #[arg(long)]
    fastq: PathBuf,
    #[arg(long)]
    truth: PathBuf,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let sim = simulate(&args.config)?;
    write_outputs(&sim, &args.fastq, &args.truth)?;
    eprintln!(
        "generated {} reads from {} molecules",
        sim.reads.len(),
        args.config.n_molecules
    );
    Ok(())
}
