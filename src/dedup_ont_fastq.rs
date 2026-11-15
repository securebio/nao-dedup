//! Deduplicate a single-end ONT FASTQ file.
//!
//! Two passes, like `dedup_interleaved_fastq`: the first builds clusters, the
//! second writes each cluster's exemplar record unchanged.

use clap::Parser;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use nao_dedup::ont::{OntDedupContext, OntDedupParams, OntMinimizerParams};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name = "dedup_ont_fastq")]
#[command(about = "Deduplicate single-end ONT FASTQ files", long_about = None)]
struct Cli {
    /// Input FASTQ file (gzip-compressed or plain; must be a regular file
    /// since it is read twice)
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Output FASTQ file (exemplars only); gzip-compressed if it ends in .gz
    #[arg(value_name = "OUTPUT")]
    output: PathBuf,

    /// Also write each read's exemplar as TSV (read_idx, exemplar_idx;
    /// 0-based positions in the input)
    #[arg(long, value_name = "PATH")]
    assignments: Option<PathBuf>,

    /// Endpoint displacement (bases) allowed at each end
    #[arg(long, default_value_t = OntDedupParams::default().max_end_offset)]
    max_end_offset: usize,

    /// Bases inward from each end used for the endpoint check
    #[arg(long, default_value_t = OntDedupParams::default().end_check_len)]
    end_check_len: usize,

    /// Edits allowed within each endpoint check window
    #[arg(long, default_value_t = OntDedupParams::default().max_end_errors)]
    max_end_errors: usize,

    /// Whole-read edit distance allowed, as a fraction of the longer read
    #[arg(long, default_value_t = OntDedupParams::default().max_error_frac, value_parser = parse_fraction)]
    max_error_frac: f64,

    /// Only match reads in the same orientation (for directional libraries)
    #[arg(long)]
    no_rc: bool,

    /// K-mer length for seeds
    #[arg(long, default_value_t = OntMinimizerParams::default().kmer_len)]
    kmer_len: usize,

    /// Seed window length
    #[arg(long, default_value_t = OntMinimizerParams::default().window_len)]
    window_len: usize,

    /// Seed windows at each end of the read
    #[arg(long, default_value_t = OntMinimizerParams::default().num_windows)]
    num_windows: usize,
}

fn parse_fraction(s: &str) -> Result<f64, String> {
    let value: f64 = s.parse().map_err(|e| format!("{e}"))?;
    if (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(format!("must be between 0 and 1, got {s}"))
    }
}

/// One FASTQ record. Lines are stored without their line endings and written
/// back with `\n`.
#[derive(Default)]
struct FastqRecord {
    header: Vec<u8>,
    sequence: Vec<u8>,
    plus: Vec<u8>,
    quality: Vec<u8>,
}

/// Reads FASTQ records, reusing one record's buffers.
struct FastqReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> FastqReader<R> {
    fn new(reader: R) -> Self {
        Self { reader }
    }

    /// Read the next record into `record`; returns false at end of input.
    fn next_into(&mut self, record: &mut FastqRecord) -> io::Result<bool> {
        if !read_line(&mut self.reader, &mut record.header)? {
            return Ok(false);
        }
        if !record.header.starts_with(b"@") {
            return Err(invalid("FASTQ header must start with '@'"));
        }
        let complete = read_line(&mut self.reader, &mut record.sequence)?
            && read_line(&mut self.reader, &mut record.plus)?
            && read_line(&mut self.reader, &mut record.quality)?;
        if !complete {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "incomplete FASTQ record",
            ));
        }
        if !record.plus.starts_with(b"+") {
            return Err(invalid("FASTQ separator line must start with '+'"));
        }
        if record.sequence.len() != record.quality.len() {
            return Err(invalid("FASTQ sequence and quality lengths differ"));
        }
        Ok(true)
    }
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

/// Read one line without its line ending; returns false at end of input.
fn read_line(reader: &mut impl BufRead, line: &mut Vec<u8>) -> io::Result<bool> {
    line.clear();
    if reader.read_until(b'\n', line)? == 0 {
        return Ok(false);
    }
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
    Ok(true)
}

/// Open a FASTQ file, decompressing it if it starts with the gzip magic bytes.
fn open_fastq(path: &Path) -> io::Result<FastqReader<Box<dyn BufRead>>> {
    let mut file = BufReader::new(File::open(path)?);
    let gzipped = file.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    let reader: Box<dyn BufRead> = if gzipped {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(file)
    };
    Ok(FastqReader::new(reader))
}

enum Output {
    Plain(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<File>>),
}

impl Output {
    /// Gzip-compressed if the path ends in .gz.
    fn create(path: &Path) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(if path.extension().is_some_and(|e| e == "gz") {
            Self::Gzip(BufWriter::new(GzEncoder::new(file, Compression::default())))
        } else {
            Self::Plain(BufWriter::new(file))
        })
    }

    fn writer(&mut self) -> &mut dyn Write {
        match self {
            Self::Plain(w) => w,
            Self::Gzip(w) => w,
        }
    }

    /// Flush and, for gzip, write the trailer, surfacing any I/O error.
    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Gzip(w) => w
                .into_inner()
                .map_err(|e| e.into_error())?
                .finish()?
                .flush(),
        }
    }
}

/// Refuse to run if an output path would overwrite the input or another
/// output, which would corrupt the second pass or lose results.
fn check_distinct_paths(cli: &Cli) -> io::Result<()> {
    // Resolve symlinks and relative paths; outputs may not exist yet, so
    // resolve their parent directory instead.
    let resolve = |path: &Path| -> io::Result<PathBuf> {
        if path.exists() {
            return path.canonicalize();
        }
        let parent = match path.parent() {
            Some(p) if !p.as_os_str().is_empty() => p.canonicalize()?,
            _ => std::env::current_dir()?,
        };
        Ok(parent.join(path.file_name().unwrap_or_default()))
    };
    let mut paths = vec![resolve(&cli.input)?, resolve(&cli.output)?];
    if let Some(assignments) = &cli.assignments {
        paths.push(resolve(assignments)?);
    }
    for (i, a) in paths.iter().enumerate() {
        if paths[i + 1..].contains(a) {
            return Err(invalid(&format!(
                "input, output, and assignments must be different files ({})",
                a.display()
            )));
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let dedup_params = OntDedupParams {
        max_end_offset: cli.max_end_offset,
        end_check_len: cli.end_check_len,
        max_end_errors: cli.max_end_errors,
        max_error_frac: cli.max_error_frac,
        rc_tolerant: !cli.no_rc,
    };
    let minimizer_params = OntMinimizerParams::new(cli.kmer_len, cli.window_len, cli.num_windows)?;
    check_distinct_paths(&cli)?;
    let mut record = FastqRecord::default();

    eprintln!("Pass 1: Building deduplication index...");
    let mut ctx = OntDedupContext::new(dedup_params, minimizer_params);
    let mut reader = open_fastq(&cli.input)?;
    let mut n_reads = 0;
    while reader.next_into(&mut record)? {
        ctx.process_read_by_index(n_reads, &record.sequence, &record.quality);
        n_reads += 1;
        if n_reads % 100_000 == 0 {
            eprintln!("  Processed {n_reads} reads...");
        }
    }
    ctx.finalize();
    let (total, clusters) = ctx.stats();
    let stats = ctx.detailed_stats();
    eprintln!("  Total reads: {total}");
    eprintln!("  Unique clusters: {clusters}");
    eprintln!("  Reads without seeds: {}", stats.unseeded_reads);
    eprintln!("  Candidate comparisons: {}", stats.candidates_compared);
    if total > 0 {
        eprintln!(
            "  Deduplication rate: {:.2}%",
            (1.0 - clusters as f64 / total as f64) * 100.0
        );
    } else {
        eprintln!("Warning: No reads found in input file");
    }

    if let Some(path) = &cli.assignments {
        let mut out = BufWriter::new(File::create(path)?);
        writeln!(out, "read_idx\texemplar_idx")?;
        for idx in 0..total {
            writeln!(out, "{idx}\t{}", ctx.get_exemplar_index(idx))?;
        }
        out.flush()?;
    }

    eprintln!("Pass 2: Writing exemplars to output...");
    let exemplars = ctx.get_exemplar_indices();
    let mut reader = open_fastq(&cli.input)?;
    let mut output = Output::create(&cli.output)?;
    let writer = output.writer();
    let (mut idx, mut written) = (0, 0);
    while reader.next_into(&mut record)? {
        if exemplars.contains(&idx) {
            for line in [
                &record.header,
                &record.sequence,
                &record.plus,
                &record.quality,
            ] {
                writer.write_all(line)?;
                writer.write_all(b"\n")?;
            }
            written += 1;
        }
        idx += 1;
    }
    if idx != n_reads {
        return Err(format!("input changed between passes ({n_reads} then {idx} reads)").into());
    }
    output.finish()?;
    eprintln!("  Wrote {written} exemplar reads");
    eprintln!("Done!");
    Ok(())
}
