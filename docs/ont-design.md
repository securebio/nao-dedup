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
`noise_scale` defaults to 1.0. Empty or exactly constant qualities instead
contribute `assumed_error_rate × length` (default rate 0.01), without quality
terms or `noise_scale`. Malformed nonempty qualities are errors.

`quality_calibration = 1.2`, `quality_error_floor = 0.002`,
`pair_distance_cv = 0.25` and the 5% tail cutoff
(`pair_tail_probability = 0.05`) are tuned engineering defaults. Library
callers can override all four; the CLI also exposes `noise_scale` for adjusting
the quality-derived noise contribution. The CLI estimates `noise_scale` before
clustering by default; it stays fixed throughout production clustering.

Pair noise sums the two reads' contributions. A pair from strains
`protect_divergence` apart has mean distance
`mu = noise + protect_divergence × shorter_length`. `T` is the largest integer
distance with lower-tail probability below `pair_tail_probability` under a
negative binomial with that mean and variance
`mu + (pair_distance_cv × noise)^2` (CV default 0.25).
A negative budget rejects even identical reads. Longer reads give more evidence
about strain differences; short or noisy reads can remain unmerged.

The default 5% is per comparison under the working model, not a bound on final
cluster contamination. A read can encounter multiple leaders, and greedy
clustering and model mismatch can change the aggregate false-merge rate.

## Automatic quality calibration

The CLI defaults to `--noise-scale auto`: a pilot at the start of the input
estimates one multiplier on the configured quality scale and floor.
`--noise-scale POSITIVE` sets it directly and skips the pilot. Calibration
does not change assumed-rate noise for reads without usable qualities.

The pilot uses greedy clustering with a 12% distance cap and only exact-endpoint
(zero-offset) pairs. Such pairs are almost surely copies of one molecule, so
their distances measure sequencing noise; tolerating offsets would admit more
chance coincidences. Pilot comparisons use `protect_divergence = 1` and
`pair_distance_cv = 0` so the cap binds. The fit's weights use the configured
CV. Only the first eligible pair per leader contributes, with usable qualities
in both reads and at least 10 expected baseline errors per pair. A weighted
median of distance-to-baseline-noise ratios is shrunk toward 1 and bounded to
[0.5, 2]. The 2x ceiling prioritizes precision: noisier inputs keep some
duplicates rather than risk merging strains.

Calibration is best effort. Fewer than 200 eligible pairs, unavailable qualities
or numerical failure fall back to the baseline multiplier 1.0. Even very large
inputs with little duplication may yield too few pairs in the sampled prefix.
Coincident endpoints from related strains, inputs concatenated from several
runs, and input order can shift the bounded, shrunk estimate. The stderr report
gives the multiplier, support and fallback or stop reason.

Collection stops at 500 retained pairs, 200,000 reads, 256 MiB of sequence plus
quality, or EOF; a prefix without usable qualities stops at 20,000 reads.
The CLI buffers sequence and quality, drops the pilot, then replays the prefix
in input order into a fresh production context and continues on the same
reader. Buffered payload is bounded to 256 MiB plus the final read; vector
capacity, record metadata, reader buffers and the pilot's seed/leader index
add overhead. The read limit also bounds per-record overhead. Library callers
using `NoiseScalePilot` or `estimate_noise_scale` must buffer sampled reads
themselves for replay; `CalibrationOptions` configures collection limits.

## CLI and storage

The CLI exposes INPUT, optional OUTPUT, `--assignments`, `--protect-divergence`,
`--noise-scale` and `--no-rc`, plus expert `--max-end-offset` and
`--assumed-error-rate`. Other matching and seed controls remain library and
example parameters. `ont_eval --auto-calibrate` estimates the same multiplier
on its configured baseline shape.

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
  insertion and deletion errors, with optional foldback tails, homopolymer
  tails and reverse complements. Quality-conditioned simulation draws
  qualities before errors, with configurable calibration, unreported errors
  and read-level variability. Simulation error-model parameters are supplied
  at runtime.
- **Metrics.** The evaluator reports pair-level precision and recall against
  molecule truth and endpoint-identifiability truth, merges, splits,
  duplicate-count bias, and recall by strand and foldback status. `ont_eval`
  adds timing and work counters.
- **Scale test.** `examples/ont_scale.rs` streams synthetic reads without
  storing them to measure throughput and peak memory at tens of millions of
  reads.

## Limitations

Independent molecules can share endpoints by chance, especially at deep
coverage of short genomes. Close strains can still merge, and the pair budget
does not guarantee a false-merge rate outside the noise model. Raw read ends,
including foldbacks, are compared without trimming (see "Read ends").
Low-complexity endpoints can be ambiguous.

Saturated groups also cost more per read: each seed matches many leaders,
and the cost grows with the number of leaders. In the synthetic scale test
with all reads from one 7 kb genome (`ont_scale --genomes 1 --genome-len 7000`,
default parameters), per-read cost went from about 43 to 74 µs as reads went
from 1M to 4M; other synthetic inputs ran at 19–33 µs per read. These figures
illustrate bucketing and storage costs; runtime depends on the machine and
comparison parameters.
