use clap::Parser;
use ont_sim::{evaluate, read_assignments, read_truth};
use std::{
    error::Error,
    fs::{self, File},
    io::BufReader,
    path::PathBuf,
};

#[derive(Parser)]
#[command(about = "Score deduplication assignments against simulated molecule truth")]
struct Args {
    #[arg(long)]
    truth: PathBuf,
    #[arg(long)]
    assignments: PathBuf,
    #[arg(long)]
    json: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let truth = read_truth(BufReader::new(File::open(args.truth)?))?;
    let assignments = read_assignments(BufReader::new(File::open(args.assignments)?))?;
    let result = evaluate(&truth, &assignments)?;
    print!("{result}");
    if let Some(path) = args.json {
        fs::write(path, format!("{}\n", result.to_json()))?;
    }
    Ok(())
}
