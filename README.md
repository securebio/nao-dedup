# nao-dedup

Sequencing read deduplication with error-tolerant matching.

This repository provides two things:

1. **Command-line tool**: Deduplicates interleaved FASTQ.gz files directly
2. **Rust library**: For integrating deduplication into other Rust programs

Both use minimizer-based bucketing for efficiency and tolerate small alignment
shifts and sequencing errors.

## Features

- **Error-tolerant matching**: Identifies duplicates even with sequencing errors
- **Alignment shift tolerance**: Handles reads that are slightly offset
- **Efficient bucketing**: Uses minimizers to avoid O(n²) comparisons
- **Quality-based selection**: Keeps the highest quality read from each cluster
- **Mate-pair swap handling**: Detects duplicates regardless of orientation

---

## Command-Line Tool

The quickest way to deduplicate your FASTQ files.

### Building

```bash
cargo build --release --bin dedup_interleaved_fastq
```

The binary will be at `target/release/dedup_interleaved_fastq`.

### Usage

Basic usage:

```bash
./target/release/dedup_interleaved_fastq input.fastq.gz output.fastq.gz
```

With custom parameters:

```bash
./target/release/dedup_interleaved_fastq \
  --max-offset 2 \
  --max-error-frac 0.02 \
  --kmer-len 15 \
  --window-len 25 \
  --num-windows 4 \
  input.fastq.gz output.fastq.gz
```

View all options:

```bash
./target/release/dedup_interleaved_fastq --help
```

### Input Format

Expects interleaved paired-end FASTQ files where R1 and R2 reads alternate
(R1, R2, R1, R2, ...). The input file must be gzip-compressed.

The input file must not be streamed, because we process it in two passes. So
you can't do:

```bash
./target/release/dedup_interleaved_fastq \
    <(aws s3 cp s3://.../input.fastq.gz -) output.fastq.gz
```

### Output

The output file contains only the exemplar read pairs (one representative per
duplicate cluster), maintaining the interleaved format.

### Parameters

- `--max-offset`: Maximum alignment offset in bases (default: 1)
- `--max-error-frac`: Maximum error fraction allowed (default: 0.01)
- `--kmer-len`: K-mer length for minimizers (default: 15)
- `--window-len`: Window length for minimizers (default: 25)
- `--num-windows`: Number of windows per read (default: 4)

### Performance

**Large file (685K reads, 456K alignment-unique)**:
- Time: 127 seconds
- Memory: 1.37 GB peak

Single-threaded; we parallelize by running multiple files concurrently.

---

## ONT (Nanopore) Deduplication

`dedup_ont_fastq` deduplicates single-end ONT reads. Copies of a fragment must
share both endpoints and match within an indel-aware edit budget. Opposite
strands match by default. See [docs/ont-design.md](docs/ont-design.md).

```bash
cargo build --release --bin dedup_ont_fastq
./target/release/dedup_ont_fastq input.fastq.gz output.fastq.gz \
  --assignments assignments.tsv
```

- **INPUT** accepts plain or gzip FASTQ, detected by content. The CLI buffers
  sequence and quality for a calibration prefix, replays it into clustering,
  then continues from the same reader. The default buffer limit is 256 MiB
  of sequence plus quality, with at most one read beyond the limit.
- **OUTPUT** contains one exemplar per group, preserving record text except
  for normalizing line endings to `\n`. A `.gz` suffix enables gzip compression;
  plain output can go to `/dev/stdout` for a parallel compressor.
  Writing FASTQ requires a second pass, so INPUT must be a regular file.
  Streams (including `-` for stdin and named pipes) require `--assignments`
  without OUTPUT; invalid combinations are rejected before input is read.
- **`--assignments PATH`** writes `read_idx`/`exemplar_idx` for every read
  using 0-based input positions. OUTPUT is optional with this flag.
- **`--protect-divergence`** sets the strain divergence against which to
  protect (default 0.015, or 1.5%). Each pair gets a cutoff based on expected
  sequencing noise and this divergence. Short or noisy reads can remain
  unmerged when the evidence is insufficient. Lower values protect closer
  strains and make merging stricter, keeping more duplicates; higher values
  merge more.
- **`--noise-scale auto|POSITIVE`** defaults to `auto`. The expected noise per
  read is `noise_scale × (scale × quality-predicted errors + per-base floor × length)`,
  with scale and floor defaulting to 1.2 and 0.002. The CLI measures one
  multiplier from the start of the input, shrinks it toward 1 and bounds it
  to [0.5, 2]. The 2x cap prioritizes protecting strains; noisier inputs may
  retain some duplicates. Calibration is best effort: even a very large
  input with little duplication can yield too few pilot pairs and fall back
  to 1.0. A number sets the multiplier directly and skips the pilot. A
  one-line report on stderr gives the automatic factor, support and stop
  reason. Fixed mode prints `Noise scale: fixed a=N`.
- **`--no-rc`** restricts matching to the same orientation.
- **Expert flags:** `--max-end-offset` sets endpoint tolerance (default 2
  bases). `--assumed-error-rate` sets per-base noise for empty or exactly
  constant qualities (default 0.01); the quality scale and floor do not apply
  to these reads. The CLI reports affected reads and comparisons.

For a single pass from plain or gzip stdin, write assignments only:

```bash
cat input.fastq.gz | ./target/release/dedup_ont_fastq - --assignments assignments.tsv
```

From Rust, use `nao_dedup::ont::OntDedupContext`: feed sequence/quality slices
with `process_read_by_index`, then `finalize` and query `get_exemplar_index`.
Library parameters include `quality_calibration` (default 1.2),
`quality_error_floor` (default 0.002), `noise_scale` (default 1.0), and endpoint,
distance and seed controls. `pair_tail_probability` overrides the per-comparison
tail cutoff (default 0.05). `noise_scale` multiplies both usable-quality noise
terms; assumed-rate noise is independent of it.
The library panics on invalid parameters, malformed nonempty Phred+33
qualities, or nonempty quality/sequence length mismatches. The CLI reports
these as errors.
`ont::validate_quality` checks the byte range without panicking.

`estimate_noise_scale` accepts an iterator of sequence/quality slices;
`NoiseScalePilot` supports incremental collection with `CalibrationOptions`
for pair, read and byte limits. Both return a `CalibrationReport`. Callers must
buffer sampled reads themselves to replay them through a fresh context with
`noise_scale = report.multiplier`.
Call `finish()` when input ends or a collection limit is reached. A report with
`stop_reason = None` means input ended before a stopping rule fired.

---

## Rust Library

For integrating deduplication into your own Rust programs.

### Building

```bash
cargo build --release
```

### API Overview

#### Requirements

- **Contiguous indices**: Pass reads with indices 0, 1, 2, ... for memory
  efficiency. Sparse indices work but waste memory.
- **Unique indices**: Each index must be used exactly once.
- **Finalize before querying**: You must call `finalize()` before calling
  `get_exemplar_index()` or `get_exemplar_indices()`.

Sequences can be upper or lowercase (both are handled).

#### Usage

```rust
use nao_dedup::{DedupContext, DedupParams, MinimizerParams};

// Create context with default or custom parameters
let dedup_params = DedupParams::default();
let minimizer_params = MinimizerParams::default();
let mut ctx = DedupContext::new(dedup_params, minimizer_params);

// Process each read pair by index
ctx.process_read_by_index(
    0,                        // read index
    "ACGTACGT".to_string(),   // forward sequence
    "TGCATGCA".to_string(),   // reverse sequence
    "IIIIIIII".to_string(),   // forward quality (Phred+33)
    "IIIIIIII".to_string(),   // reverse quality (Phred+33)
);
ctx.process_read_by_index(
    1,
    "ACGTACGT".to_string(),   // duplicate of read 0
    "TGCATGCA".to_string(),
    "IIIIIIII".to_string(),
    "IIIIIIII".to_string(),
);

// Finalize before querying
ctx.finalize();

// Query results
let exemplar_idx = ctx.get_exemplar_index(0);  // Returns 0 (itself)
let exemplar_idx = ctx.get_exemplar_index(1);  // Returns 0 (clustered with read 0)

// Get set of all exemplar indices (for filtering output)
let exemplar_indices = ctx.get_exemplar_indices();

// Get statistics
let (total_processed, unique_clusters) = ctx.stats();
```

### Integration Example

See `src/dedup_interleaved_fastq.rs` for a complete example of how to integrate
this library with file I/O.

---

## How It Works

### Two-Pass Architecture

The library processes input in two passes:

**Pass 1 (Clustering)**: Reads stream in one at a time. Each read is compared
against existing cluster exemplars using minimizer-based bucketing. If it
matches an existing cluster, it joins that cluster; otherwise it starts a new
cluster. The library tracks the best read (by quality) seen so far in each
cluster, but during this pass each read's result points to its *cluster leader*
(the first read that started the cluster), not the current best.

**Pass 2 (Finalization)**: After all reads are processed, the library resolves
each read's assignment from the cluster leader to the cluster's *best exemplar*.
This second pass is necessary because the best read in a cluster can change as
new reads arrive, and we can't update previously-processed reads during
streaming.

Reads are tracked by their linear position (index) rather than by sequence ID,
avoiding string allocation and hash lookups.

### 1. Minimizer Extraction

Each read is divided into windows, and the lexicographically smallest k-mer
(minimizer) is extracted from each window. This creates a signature for each
read pair.

#### K-mer Encoding

K-mers are encoded using a 2-bit DNA encoding (A=0, C=1, G=2, T=3) that
represents each base with exactly 2 bits. This provides several advantages when
used as a hash key:

- **Fast**: Just bit shifts and ORs, much faster than CRC32 or polynomial hashing
- **No collisions**: Bijective mapping for k-mers up to length 32 (fits in 64 bits)
- **Consistent**: Deterministic encoding produces identical hash values
- **DNA-aware**: Leverages the 4-base structure of DNA sequences

K-mers containing non-ACGT bases (primarily N) return a sentinel value,
ensuring they won't be selected as minimizers.

#### Window Placement

Windows are placed adjacently starting from the beginning of each read
(positions 0, window_len, 2× window_len, etc.). This strategy:
- Focuses on the most stable region of reads (the beginning)
- Avoids the tail region, which is most likely to be trimmed or contain
  sequencing errors

For example, with 3 windows of 25bp on a 150bp read:
- Window 0: positions [0, 25)
- Window 1: positions [25, 50)
- Window 2: positions [50, 75)
- Positions [75, 150) are not examined for minimizers

### 2. Bucketing

Read pairs with matching minimizers are assigned to the same buckets. This
dramatically reduces the number of pairwise comparisons needed.

- Reads are processed one at a time
- Uses individual minimizer hashes as keys
- Generates 2n keys per read (the minimizers from forward and reverse reads)
- Only stores unique exemplars, not all reads

### 3. Pairwise Comparison

Within each bucket, read pairs are compared to determine if they're
duplicates. Comparison allows for:
- Small alignment offsets (configurable via `max_offset`)
- Sequencing errors (configurable via `max_error_frac`)
- Mate-pair orientation swaps (always enabled)

### 4. Clustering

- No graph construction - purely streaming approach
- First read in a cluster identifies that cluster (becomes the cluster leader)
- As reads are processed, the best read (by quality and length) becomes the
  exemplar
- Subsequent reads matching the cluster are compared against the current
  exemplar
- If a better read is found, it replaces the previous exemplar

### 5. Exemplar Selection

Exemplars are selected based on:
1. Mean quality score (higher is better)
2. Total read length (longer is better)
3. First read in cluster serves as tie-breaker

## Development Setup

### Rust Toolchain

On macOS with Homebrew:
```bash
brew install rust
```

On Linux:
```bash
# Ubuntu/Debian
sudo apt install rustc cargo

# Redhat derived distributions (Fedora, Amazon Linux 2023, etc)
sudo dnf install rust cargo
```

Verify installation:
```bash
cargo --version
```

## Testing

Run all tests:

```bash
cargo test
```

The test suite includes:

- **Unit tests**: Located in `src/lib.rs` (inline tests), testing helper
  functions, minimizer extraction, sequence matching, and parameter validation
- **Integration tests**: Located in `tests/dedup_tests.rs`, testing the public
  API with various deduplication scenarios including edge cases with N's,
  offsets, and errors
- **CLI tests**: Located in `tests/cli_tests.rs`, testing the command-line
  binary with gzipped FASTQ files

To run specific test suites:
```bash
cargo test --lib           # Unit tests only
cargo test --test dedup_tests  # Integration tests
cargo test --test cli_tests    # CLI tests
cargo test --test ont_tests --test ont_cli_tests  # ONT library and CLI
```
