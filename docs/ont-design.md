# ONT deduplication: design

The `ont` library and `dedup_ont_fastq` CLI identify duplicates among single-end
Oxford Nanopore reads. Two reads must share both fragment endpoints within
`max_end_offset` (default 2 bases) and match within an indel-aware edit budget.
Reverse-complement matching is enabled by default; directional libraries can
turn it off. Truncated copies sharing only one endpoint are outside this definition.

Endpoints are the main evidence. Two independent molecules from the same
genome region can be nearly identical in sequence. A tight endpoint
tolerance is what keeps them apart. At deep coverage of a short genome,
independent molecules can share endpoints by chance, so groups are operational
("same endpoints and sequence"), not guaranteed molecular counts.

## Read ends

Raw read ends are compared as they are. Some ONT reads end in sequence that
isn't part of the template's forward strand, most notably foldback reads,
which continue into the reverse complement of the molecule. Reads are not
trimmed for this. If a foldback forms before PCR, all copies of the
molecule carry it, and they match without trimming; trimming each read
independently would only add jitter at the cut. Libraries where foldbacks
arise independently in each read would need a trimming step.

Single-base homopolymer k-mers are never used as seeds, so poly-A/T tails
don't create giant buckets.

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
   (fixed 0.12), before unpacking the leader's ends.
2. At each end, align `end_check_len` bases inward from the boundary
   (fixed 32). Up to `max_end_offset` bases of displacement are free;
   at most `max_end_errors` edits (fixed 3) are allowed. The inner boundary
   is free. Reject if allowing a larger displacement makes the alignment
   strictly cheaper, so the error allowance does not widen endpoint tolerance.
3. Calculate a quality-dependent pair budget `T`, capped at `max_error_frac`.
   The length difference must fit `T`.
4. Unpack the full leader and compute bounded whole-read edit distance with
   `editdistancek`; accept only distances at most `T`.

## Pair budget

The expected noise per read is
`quality_calibration × quality-predicted errors + quality_error_floor × length`.
Quality-predicted errors sum `10^(-Q/10)` over bases. Empty or exactly constant
qualities instead contribute `assumed_error_rate × length` (default rate 0.01),
without quality terms. Malformed nonempty qualities are errors.

`quality_calibration = 1.2`, `quality_error_floor = 0.002`,
`pair_distance_cv = 0.25` and the 5% tail cutoff
(`pair_tail_probability = 0.05`) are fixed engineering defaults, as are the
endpoint window, endpoint error allowance, whole-read cap and seed drift
(`seed_drift_frac = 0.04`).

Pair noise sums the two reads' contributions. A pair from strains
`protect_divergence` apart has mean distance
`mu = noise + protect_divergence × shorter_length`. `T` is the largest integer
distance with lower-tail probability below `pair_tail_probability` under a
negative binomial with that mean and variance
`mu + (pair_distance_cv × noise)^2` (fixed CV 0.25).
A negative budget rejects even identical reads. Longer reads give more evidence
about strain differences; short or noisy reads can remain unmerged.

The default 5% is per comparison under the working model, not a bound on final
cluster contamination. A read can encounter multiple leaders, and greedy
clustering and model mismatch can change the aggregate false-merge rate.

## CLI and storage

The CLI exposes INPUT, optional OUTPUT, `--assignments`, `--protect-divergence`,
and `--no-rc`, plus expert `--max-end-offset` and
`--assumed-error-rate`. Matching constants are fixed; seed controls remain
library and example parameters.

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

In synthetic tests with reads of 0.5–4.5 kb (mean ~2.5 kb), peak memory was
about 0.9 KB per leader when leaders shared seeds and about 1.4 KB when every
seed was distinct. Packed sequence costs a quarter byte per base, so longer
reads cost proportionally more. At that read length with shared seeds, even
50M unique reads fit in 64 GB. `ont_scale` streams synthetic reads to measure
throughput and peak memory at that scale.

## Validation

- **Evidence loss and near misses.** The evaluator reports lost molecules,
  lost genomes and extra exemplars. `examples/ont_nearmiss.rs` measures pair
  merge probability across sequence divergence and endpoint offsets through
  the full deduplicator, with an identical-template duplicate control. Each
  trial uses a fresh context; mixed-strain clustering tests also check whether
  a rare strain survives many comparisons in different arrival orders.
- **Simulator.** `ont-sim` generates genomes, optionally related by a set
  divergence, and molecules with defined endpoints, including dense coverage
  and nearby distinct endpoints. PCR families get independent substitution,
  insertion and deletion errors, with homopolymer indel enrichment and optional
  reverse complements. Quality-conditioned simulation draws
  qualities before errors, with configurable calibration, unreported errors
  and read-level variability. Simulation error-model parameters are supplied
  at runtime.
- **Metrics.** The evaluator reports pair-level precision and recall against
  molecule truth and endpoint-identifiability truth. `ont_eval` adds timing
  and work counters.
- **Scale test.** `examples/ont_scale.rs` streams quality-conditioned reads
  without storing them, reporting deduplication time, generation time and peak
  memory. It requires `--sub-weight`, `--ins-weight`, `--del-weight`,
  `--error-median`, `--error-sigma`, `--local-q-sigma`, `--miscalibration` and
  `--latent-sigma`; add `--genomes 1 --genome-len 7000` for saturated seed groups.

## Limitations

Independent molecules can share endpoints by chance, especially at deep
coverage of short genomes. Close strains can still merge, and the pair budget
does not guarantee a false-merge rate outside the noise model. Raw read ends,
including foldbacks, are compared without trimming (see "Read ends").
Low-complexity endpoints can be ambiguous.

Saturated groups also cost more per read: each seed matches many leaders,
and the cost grows with the number of leaders. In the original constant-quality
scale test with all reads from one 7 kb genome (`ont_scale --genomes 1
--genome-len 7000`, default parameters), per-read cost went from about 43 to 74 µs
as reads went from 1M to 4M; other synthetic inputs ran at 19–33 µs per read.
These figures illustrate bucketing and storage costs; runtime depends on the
machine and comparison parameters.
