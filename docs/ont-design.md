# ONT deduplication: design

The `ont` library and `dedup_ont_fastq` CLI identify duplicates among single-end
Oxford Nanopore reads. Two reads must share both fragment endpoints within
`max_end_offset` (default 2 bases) and match within an indel-aware edit budget.
Reverse-complement matching is enabled by default; directional libraries can
turn it off. Truncated copies sharing only one endpoint are outside this definition.

## Clustering and comparison

Reads stream into greedy clusters. Each cluster's first read is its immutable
comparison leader. A new read joins the first matching leader; clusters never
merge. The exemplar is the member with highest mean Phred quality, then longest
length. `finalize()` resolves assignments to exemplars and releases the seed
index, leader storage, evidence and comparison buffers. Members match their
leader, and need not match the final exemplar directly.

Canonical minimizers from adjacent windows at both ends retrieve candidates.
A candidate needs seed evidence at both ends in a consistent orientation.
Seeds are visited from the ends inward; evidence chooses which orientation to
try first, followed by the other allowed orientation if necessary. Single-base
homopolymer k-mers are excluded. Postings are sorted by offset, so only offsets
within the endpoint allowance plus `seed_drift_frac` drift are visited.
This drift filter is a heuristic in repeats and beyond the end-check window.

Each candidate passes these stages:

1. The length difference must fit `max_error_frac` of the longer read
   (default 0.12), before unpacking the leader's ends.
2. At each end, align `end_check_len` bases inward from the boundary
   (default 32). Up to `max_end_offset` bases of displacement are free;
   at most `max_end_errors` edits (default 3) are allowed. The inner boundary
   is free. Reject if allowing a larger displacement makes the alignment
   strictly cheaper, so the error allowance does not widen endpoint tolerance.
3. Calculate a quality-dependent pair budget `T`, capped at `max_error_frac`.
   The length difference must fit `T`.
4. Unpack the full leader and compute bounded whole-read edit distance with
   `editdistancek`; accept only distances at most `T`.

## Pair budget

The expected noise per read is
`noise_scale × (quality_calibration × quality-predicted errors + quality_error_floor × length)`.
Quality-predicted errors sum `10^(-Q/10)` over bases; library parameters
`quality_calibration` and `quality_error_floor` default to 1.2 and 0.002;
`noise_scale` defaults to 1.0 and can be overridden. Empty or exactly constant
qualities instead contribute `assumed_error_rate × length` (default rate 0.01),
without quality terms or `noise_scale`. Malformed nonempty qualities are errors.

Pair noise sums the two reads' contributions. A pair from strains
`protect_divergence` apart has mean distance
`mu = noise + protect_divergence × shorter_length`. `T` is the largest integer
distance with lower-tail probability below 5% under a negative binomial with
that mean and variance `mu + (pair_distance_cv × noise)^2` (CV default 0.25).
A negative budget rejects even identical reads. Longer reads give more evidence
about strain differences; short or noisy reads can remain unmerged.

## CLI and storage

The CLI exposes INPUT, optional OUTPUT, `--assignments`, `--protect-divergence`,
`--noise-scale` and `--no-rc`, plus expert `--max-end-offset` and
`--assumed-error-rate`. Other matching and seed controls remain library and
example parameters.

Input can be plain or gzip FASTQ, detected from the first two bytes even when
reads from a stream are short. Assignments are written after finalization in
input order. Writing exemplar FASTQ rereads INPUT, which must be a regular
file; the CLI checks this before consuming input. Streams, including `-` for
stdin and named pipes, require `--assignments` without OUTPUT. Pass 2 selects
records whose assignment points to themselves, preserving record text except
for normalized line endings. A `.gz` output suffix enables gzip compression.

Leaders use 2 bits per ACGT base with verbatim fallback for other bytes, plus
a 4-byte expected-error summary. Assignments use one `u32` per read. The seed
index maps hashes to offset-sorted 8-byte postings (leader, offset, end), up to
16 postings per leader by default. A bucket holds one posting inline.

## Limitations and evaluation tools

Independent molecules can share endpoints by chance, especially at deep
coverage of short genomes. Close strains can still merge, and the pair budget
does not guarantee a false-merge rate outside the noise model. Raw read ends,
including foldbacks, are compared without trimming. Low-complexity endpoints
can be ambiguous.

`ont-sim` generates reads and evaluates molecule/endpoint assignments;
`ont_eval` adds timing and work counters. `ont_nearmiss` measures pair merge
probability across divergence and endpoint offsets. `ont_scale` streams reads
to measure throughput and peak memory. Error rates are supplied at runtime.
