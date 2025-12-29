use clap::Parser;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use nao_dedup::{DedupContext, DedupParams, MinimizerParams, ReadPair};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "dedup_interleaved_fastq")]
#[command(about = "Deduplicate interleaved paired-end FASTQ files", long_about = None)]
struct Cli {
    /// Input FASTQ.gz file (interleaved R1/R2)
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Output FASTQ.gz file (exemplars only)
    #[arg(value_name = "OUTPUT")]
    output: PathBuf,

    /// Maximum alignment offset (default: 1)
    #[arg(long, default_value_t = 1)]
    max_offset: usize,

    /// Maximum error fraction (default: 0.01)
    #[arg(long, default_value_t = 0.01)]
    max_error_frac: f64,

    /// K-mer length for minimizers (default: 15)
    #[arg(long, default_value_t = 15)]
    kmer_len: usize,

    /// Window length for minimizers (default: 25)
    #[arg(long, default_value_t = 25)]
    window_len: usize,

    /// Number of windows for minimizers (default: 4)
    #[arg(long, default_value_t = 4)]
    num_windows: usize,
}

#[derive(Debug)]
struct FastqRecord {
    header: String,
    sequence: String,
    plus: String,
    quality: String,
}

fn read_fastq_record<R: BufRead>(reader: &mut R) -> std::io::Result<Option<FastqRecord>> {
    let mut header = String::new();
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();

    // Read header
    if reader.read_line(&mut header)? == 0 {
        return Ok(None);
    }

    // Read sequence
    if reader.read_line(&mut sequence)? == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "Incomplete FASTQ record: missing sequence",
        ));
    }

    // Read + line
    if reader.read_line(&mut plus)? == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "Incomplete FASTQ record: missing plus line",
        ));
    }

    // Read quality
    if reader.read_line(&mut quality)? == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "Incomplete FASTQ record: missing quality",
        ));
    }

    Ok(Some(FastqRecord {
        header,
        sequence,
        plus,
        quality,
    }))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    // Set up parameters
    let dedup_params = DedupParams {
        max_offset: cli.max_offset,
        max_error_frac: cli.max_error_frac,
    };

    let minimizer_params = MinimizerParams {
        kmer_len: cli.kmer_len,
        window_len: cli.window_len,
        num_windows: cli.num_windows,
    };

    eprintln!("Pass 1: Building deduplication index...");

    // Pass 1: Read all pairs and build deduplication index
    let input_file = File::open(&cli.input)?;
    let gz_decoder = GzDecoder::new(input_file);
    let mut reader = BufReader::new(gz_decoder);

    let mut ctx = DedupContext::new(dedup_params, minimizer_params);
    let mut pair_index = 0;

    loop {
        // Read R1
        let r1 = match read_fastq_record(&mut reader)? {
            Some(record) => record,
            None => break,
        };

        // Read R2
        let r2 = match read_fastq_record(&mut reader)? {
            Some(record) => record,
            None => {
                eprintln!("Warning: Odd number of reads in file. Last read ignored.");
                break;
            }
        };

        // Create ReadPair with index as ID
        let read_pair = ReadPair {
            read_id: pair_index.to_string(),
            fwd_seq: r1.sequence.trim().to_string(),
            rev_seq: r2.sequence.trim().to_string(),
            fwd_qual: r1.quality.trim().to_string(),
            rev_qual: r2.quality.trim().to_string(),
        };

        ctx.process_read(read_pair);
        pair_index += 1;

        if pair_index % 100_000 == 0 {
            eprintln!("  Processed {} read pairs...", pair_index);
        }
    }

    eprintln!("  Total read pairs: {}", pair_index);

    // Finalize deduplication
    eprintln!("Finalizing deduplication...");
    ctx.finalize();

    let (total_reads, unique_clusters) = ctx.stats();
    eprintln!("  Total reads: {}", total_reads);
    eprintln!("  Unique clusters: {}", unique_clusters);
    eprintln!(
        "  Deduplication rate: {:.2}%",
        (1.0 - unique_clusters as f64 / total_reads as f64) * 100.0
    );

    // Build set of exemplar indices
    let mut exemplar_indices = HashSet::new();
    for i in 0..pair_index {
        let exemplar = ctx.get_cluster_id(&i.to_string());
        if exemplar == i.to_string() {
            exemplar_indices.insert(i);
        }
    }

    eprintln!("Pass 2: Writing exemplars to output...");

    // Pass 2: Write exemplars
    let input_file = File::open(&cli.input)?;
    let gz_decoder = GzDecoder::new(input_file);
    let mut reader = BufReader::new(gz_decoder);

    let output_file = File::create(&cli.output)?;
    let gz_encoder = GzEncoder::new(output_file, Compression::default());
    let mut writer = BufWriter::new(gz_encoder);

    let mut current_index = 0;
    let mut written = 0;

    loop {
        // Read R1
        let r1 = match read_fastq_record(&mut reader)? {
            Some(record) => record,
            None => break,
        };

        // Read R2
        let r2 = match read_fastq_record(&mut reader)? {
            Some(record) => record,
            None => break,
        };

        // Write if this is an exemplar
        if exemplar_indices.contains(&current_index) {
            write!(writer, "{}", r1.header)?;
            write!(writer, "{}", r1.sequence)?;
            write!(writer, "{}", r1.plus)?;
            write!(writer, "{}", r1.quality)?;
            write!(writer, "{}", r2.header)?;
            write!(writer, "{}", r2.sequence)?;
            write!(writer, "{}", r2.plus)?;
            write!(writer, "{}", r2.quality)?;
            written += 1;
        }

        current_index += 1;

        if current_index % 100_000 == 0 {
            eprintln!("  Processed {} read pairs...", current_index);
        }
    }

    eprintln!("  Wrote {} exemplar pairs", written);
    eprintln!("Done!");

    Ok(())
}
