# ONT deduplication: design

The `ont` library and `dedup_ont_fastq` CLI identify duplicates among single-end
Oxford Nanopore reads. Two reads are duplicates when they share both fragment
endpoints within `max_end_offset` (default 2 bases) and their sequences match
within a quality-dependent, indel-aware edit budget. Reverse-complement
matching is on by default; directional libraries can turn it off. Truncated
copies sharing only one endpoint are outside this definition.

Endpoints are the main evidence. Independent molecules from the same genome
region can be nearly identical in sequence, and a tight endpoint tolerance is
what keeps them apart. At deep coverage of a short genome, independent
molecules can share endpoints by chance, so groups are operational ("same
endpoints and sequence"), not guaranteed molecular counts.

## Read ends

Read ends are taken to be fragment ends, so adapters, barcodes and primers
should be trimmed first. Sequence shared by many reads at their ends (an
untrimmed adapter, say) defeats the endpoint check: if it covers the 32-base
endpoint window, distinct fragments whose true ends are a few bases apart can
merge on the whole-read budget alone, and every read lands in the same seed
buckets. Ends are otherwise compared raw.

Some reads continue into the reverse complement of their own sequence. This
can arise in two ways, and reads are not trimmed for either:

- **At the pore (foldback reads).** The complementary strand follows the first
  strand through the pore. This happens independently in each read, so a
  foldback read has a different end and length from other copies of its
  molecule and stays unmerged: a missed duplicate, not a false merge.
- **During library preparation.** A hairpin or inverted chimera formed before
  amplification is part of the molecule, so every PCR copy carries it and the
  copies match as they are. Trimming each read independently would only add
  jitter at the cut.

Libraries with many pore foldbacks would need a trimming step to recover those
duplicates.

Single-base homopolymer k-mers are never used as seeds, so poly-A/T tails don't
create giant buckets.

## Clustering and comparison

Reads stream into greedy clusters. Each cluster's first read is its immutable
comparison leader; a new read joins the first matching leader in the order
candidates are found (below), and clusters never merge. The exemplar is the
member with the highest mean Phred quality, then the longest. Members match
their leader, not necessarily the exemplar. `finalize()` resolves assignments to
exemplars and frees the seed index, leader sequences and comparison buffers,
keeping only assignments and each cluster's exemplar.

Reads are clustered in input order, without sorting first. Sorting longest
first helps when shorter reads are absorbed into longer ones; here duplicates
share both endpoints, so cluster members are already nearly the same length.
Sorting would also give up single-pass streaming.

New reads are compared with the fixed leader rather than the current exemplar.
If the comparison target moved to each better member, a cluster could drift:
every replacement is a member within tolerance of the previous target, so
endpoints and sequence could creep further from the first read with each step,
widening the effective endpoint tolerance. A fixed leader keeps every member
within one comparison of the same read, and its seeds are indexed once. A noisy
leader mainly costs recall: the pair budget accounts for the leader's predicted
noise, so the protected-divergence test still applies under the noise model.

### Finding candidates

Each read yields seeds from both ends: the canonical minimizer (15-mer) of each
of up to 8 adjacent 32-base windows (by default) laid out inward from the start
and inward from the end. Short reads have fewer windows, and a window without a
valid non-homopolymer ACGT k-mer gives no seed. Each seed keeps its offset from
its end: k-mer start minus read start, or read end minus k-mer end. Canonical
k-mers make a read and its reverse complement produce the same seeds with the
ends swapped.

Leaders' seeds go into one index, a hash map from seed to a bucket of
*postings*. A posting records one occurrence of the seed in a leader: the
leader, which of its ends the seed came from, and the seed's offset from that
end. Each bucket is sorted by offset.

The Illumina path buckets reads by minimizers from fixed adjacent windows at
each mate's start. Here endpoints may shift by a few bases and indels move
seeds, so a shared seed is only evidence, and offsets decide which postings
count. For a new read, seeds are looked up window by window, from the ends
inward. For each seed only the postings whose offset is close to the read's are
visited: within `max_end_offset + 3` plus 4% of the larger offset, rounded up
(indel drift). Each visited posting records which read end matched which leader
end. A leader becomes a candidate once it has evidence at both ends in a
consistent orientation: read start with leader start and read end with leader
end (same strand), or the crossed pairing (opposite strand).

Candidates are not ranked. Seeds are processed window by window from the ends
inward; in each window the start seed's postings are visited before the end
seed's, each bucket in offset order, then leader creation order. Leaders queue
as their evidence completes, the queue is aligned after both lookups, and the
first match ends the search. Seeds near the ends are the most position-specific,
so true duplicates usually qualify after a few lookups, without scanning the
large buckets of inner seeds. Evidence picks the orientation to try first; the
other allowed orientation is tried if needed. Beyond the endpoint window and in
repeats the offset filter is a heuristic and can miss a duplicate.

### Comparing a candidate

Each candidate passes these stages, cheapest first:

1. The length difference must fit the whole-read cap, 12% of the longer read.
   Sequencing indels change length, so duplicates' lengths can differ by more
   than the `2 × max_end_offset` that endpoint displacement allows. Every
   base of length difference costs an edit, though, so the edit budget bounds
   it.
2. At each end, align up to 32 bases inward from the boundary, with the outer
   boundary anchored and the inner one free. Up to `max_end_offset` bases of
   displacement are free and at most 3 edits are allowed. Reject if allowing
   displacement up to `max_end_offset + 3` makes the alignment strictly
   cheaper, so the error allowance cannot widen the endpoint tolerance.
3. Compute the pair budget `T` (below), capped at 12%. The length difference
   must fit `T`.
4. Compute the bounded whole-read edit distance; accept if it is at most `T`.

## Pair budget

Each read's expected sequencing noise is `1.2 × Σ 10^(-Q/10) + 0.002 ×
length`: quality-predicted errors with a calibration factor and a per-base
floor for errors that qualities miss. Reads with empty or exactly constant
qualities instead contribute `assumed_error_rate × length` (default 0.01).
Malformed nonempty qualities are errors.

The budget asks whether the pair could come from two strains
`protect_divergence` apart (default 1.5%). Such a pair has mean distance
`mu = noise_a + noise_b + protect_divergence × shorter_length`. `T` is the
largest integer up to the 12% cap with `P(D ≤ T) < 5%` when `D` follows a
negative binomial with mean `mu` and variance
`mu + (0.25 × (noise_a + noise_b))^2` (Poisson when the noise is zero); the
extra term allows for noise misprediction. If no such `T ≥ 0` exists, the pair
cannot merge. Longer reads carry more evidence about strain differences and
get more room relative to their noise; short or noisy reads may stay unmerged,
even when identical.

The 5% is per comparison under this model, not a bound on cluster
contamination: a read can meet several leaders, and model mismatch changes the
real rate. The quality coefficients, dispersion, tail probability, endpoint
window and its error allowance, whole-read cap and seed drift are fixed
engineering defaults. They are documented constants at the top of
`src/ont/compare.rs`. The tests in `tests/ont_tests.rs` named
`default_params_*` pin the behavior they produce. The simulator, parameter
sweeps (`examples/ont_eval.rs`) and near-miss merge-probability grid
(`examples/ont_nearmiss.rs`) used to choose them are on the
`coding-agent/ont-tuning-toolkit` branch.

## CLI and storage

The CLI takes INPUT, optional OUTPUT, `--assignments`, `--protect-divergence`
and `--no-rc`, plus expert `--max-end-offset` and `--assumed-error-rate` (see
the README). Input is plain or gzip FASTQ, detected from content. Assignments
are written after finalization in input order. Writing exemplar FASTQ rereads
INPUT, so it must be a regular file; streams such as stdin require
`--assignments` without OUTPUT.

Leaders are packed at 2 bits per ACGT base; a leader containing any other byte
is stored verbatim. Each also keeps a 4-byte expected-error summary. Assignments
use one `u32` per read. The seed index maps hashes to offset-sorted 8-byte
postings, up to 16 per leader by default. In synthetic tests with 0.5–4.5 kb
reads, peak memory was about 0.9 KB per leader when leaders shared seeds and
1.4 KB when every seed was distinct, so tens of millions of unique reads of that
length fit in 64 GB. Packed sequence costs a quarter byte per base, so longer
reads cost proportionally more.

## Validation

- Unit tests check endpoint DP against a full-DP oracle, rolling hashes against
  an oracle, and negative-binomial tail goldens.
- Integration tests use generated reads with explicit illustrative error inputs
  to check retrieval and endpoint precision/recall, quality-conditioned
  duplicate recall at defaults, related-strain protection at defaults, and rare
  strains in dense mixed clusters. The shared generator is in
  `tests/common/sim.rs`.
- CLI tests cover FASTQ, gzip, streams, assignments, and input validation.

The simulator, parameter sweeps and near-miss experiments used to tune the
defaults live on the separate `coding-agent/ont-tuning-toolkit` branch, along
with the `ont_scale` throughput and memory benchmark.

## Limitations

Independent molecules can share endpoints by chance, especially at deep
coverage of short genomes. Close strains can still merge; the pair budget does
not guarantee a false-merge rate outside its noise model. Low-complexity
endpoints can be ambiguous. Untrimmed shared end sequence weakens the
endpoint check (see "Read ends"). Some duplicates are always missed: reads
without seeds at both ends become singletons, and candidates are found only
through seeds in the outer 256 bases of each end (8 windows of 32 by default).

Runtime is single-threaded. In dense groups (many reads from one short genome)
each read meets more candidate leaders as the group grows, so the cost per
read rises with group size; reads that start new clusters must check every
candidate.
