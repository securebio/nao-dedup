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

Raw read ends are compared as they are. Some ONT reads end in sequence that
isn't part of the template's forward strand, most notably foldback reads,
which continue into the reverse complement of the molecule. If a foldback
forms before PCR, all copies of the molecule carry it and match without
trimming; trimming each read independently would only add jitter at the cut.
Libraries where foldbacks arise independently in each read would need a
trimming step. Single-base homopolymer k-mers are never used as seeds, so
poly-A/T tails don't create giant buckets.

## Clustering and comparison

Reads stream into greedy clusters. Each cluster's first read is its immutable
comparison leader; a new read joins the first leader it matches, and clusters
never merge. The exemplar is the member with the highest mean Phred quality,
then the longest. Members match their leader, not necessarily the exemplar.
`finalize()` resolves assignments to exemplars and frees the seed index,
leader sequences and comparison buffers, keeping only assignments and each
cluster's exemplar.

Canonical minimizers from adjacent windows at both ends retrieve candidates. A
candidate needs seed evidence at both ends in a consistent orientation. Seeds
are visited from the ends inward, and postings are sorted by offset, so only
postings whose offset from the end differs by at most `max_end_offset + 3`
plus 4% of the larger offset, rounded up (indel drift) are visited. Beyond the endpoint window
and in repeats this filter is a heuristic and can miss a duplicate. Evidence
picks the orientation to try first; the other allowed orientation is tried if
needed.

Each candidate passes these stages, cheapest first:

1. The length difference must fit the whole-read cap, 12% of the longer read.
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
engineering defaults.

## CLI and storage

The CLI takes INPUT, optional OUTPUT, `--assignments`, `--protect-divergence`
and `--no-rc`, plus expert `--max-end-offset` and `--assumed-error-rate` (see
the README). Input is plain or gzip FASTQ, detected from content. Assignments
are written after finalization in input order. Writing exemplar FASTQ rereads
INPUT, so it must be a regular file; streams such as stdin require
`--assignments` without OUTPUT.

Leaders are packed at 2 bits per ACGT base, with verbatim fallback for other
bytes, plus a 4-byte expected-error summary. Assignments use one `u32` per
read. The seed index maps hashes to offset-sorted 8-byte postings, up to 16
per leader by default. In synthetic tests with 0.5–4.5 kb reads, peak memory
was about 0.9 KB per leader when leaders shared seeds and 1.4 KB when every
seed was distinct, so tens of millions of unique reads of that length fit in
64 GB. Packed sequence costs a quarter byte per base, so longer reads cost
proportionally more.

## Validation

- **Near misses.** `examples/ont_nearmiss.rs` measures pair merge probability
  across sequence divergence and endpoint offset through the full
  deduplicator, with an identical-template duplicate control. Mixed-strain
  tests check that a rare strain survives many comparisons in different
  arrival orders.
- **Simulation.** `ont-sim` generates related genomes, molecules with dense or
  nearby endpoints, PCR families, and reads with substitutions, indels,
  reverse complements and quality-conditioned errors. Its evaluator reports
  pair precision and recall and evidence loss (lost molecules and genomes,
  extra exemplars). Error-model parameters are supplied at runtime.
- **Scale.** `examples/ont_scale.rs` streams quality-conditioned reads without
  storing them and reports deduplication time and peak memory.

## Limitations

Independent molecules can share endpoints by chance, especially at deep
coverage of short genomes. Close strains can still merge; the pair budget does
not guarantee a false-merge rate outside its noise model. Low-complexity
endpoints can be ambiguous.

Runtime is single-threaded. In dense groups (many reads from one short genome)
each read meets more candidate leaders as the group grows, so the cost per
read rises with group size; reads that start new clusters must check every
candidate.
