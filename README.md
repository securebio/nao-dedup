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
```
