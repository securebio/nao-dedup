# ONT deduplication: design

Status: proposed design for the `ont` module. It will be updated as the
implementation and measurements progress.

## Goal

Identify PCR duplicates among Oxford Nanopore (ONT) single-end reads and
choose one exemplar per duplicate group, in the style of the Illumina
`DedupContext`. Inputs are read sets of up to tens of millions of reads
(lengths up to ~15 kb). The target is to fit in 64 GB of RAM and finish in
under an hour at that scale.

## What counts as a duplicate

PCR copies of one template molecule share that molecule's two endpoints. Two
reads are duplicates when:

1. **Both endpoints agree** within a small tolerance (`max_end_offset`, a few
   bases).
2. **The sequences are compatible** with independent sequencing errors: the
   whole-read edit distance is at most `max_error_frac` of the longer read's
   length.
3. **The orientation is allowed**: the same strand, or either strand when
   `rc_tolerant` is set (the default). Non-directional library preps sequence
   PCR copies from both strands. Directional preps can disable
   `rc_tolerant`.

Endpoints are the main evidence. Two independent molecules from the same
genome region can be nearly identical in sequence. A tight endpoint
tolerance, not a tight identity threshold, is what keeps them apart.

One limitation is inherent to any reference-free or coordinate-based
method. When coverage of a short genome is very deep, independent molecules
often share endpoints by chance. Groups are therefore operational ("same
endpoints and sequence"), not guaranteed molecular counts.

Separately, this version deliberately does not merge truncated copies, which
share only one endpoint.

## Read normalization

Raw ONT read ends can include sequence that is not part of the template
fragment. Examples are foldback reads, where the read continues into the
reverse complement of the molecule, and homopolymer tails.
Requiring raw ends to match would miss these duplicates. Loosening
`max_end_offset` enough to absorb them would give up endpoint specificity.

Instead, each read gets a *virtual* biological interval. Conservatively
detected foldback tails are excluded from it. Seeding, comparison, and
quality scoring all use this interval, but output keeps the original record
unchanged. If detection is uncertain, the read is left untrimmed. Homopolymer
tails are not trimmed; low-complexity k-mers are simply not used as seeds.

## Algorithm

The structure is the same as the Illumina path: bucket reads by minimizers,
compare each new read only against cluster leaders that share a bucket,
cluster greedily, and resolve the best exemplar at finalize.

1. **Seeds.** Illumina takes minimizers from adjacent windows at the start of
   each mate. ONT treats the read's two ends like two mates. It takes one
   minimizer from each of `num_windows` adjacent windows of `window_len`
   bases, counted inward from both the start and the end of the interval.
   - K-mers are canonical, so each read and its reverse complement give the
     same seeds.
   - Minimizers are chosen by a mixing hash, not the lexicographic 2-bit
     value, so poly-A k-mers aren't systematically preferred.
2. **Buckets.** The primary keys are pairs of seeds, one from each end.
   Every start-window seed is combined with every end-window seed. A pair is
   stored unordered, because reverse complementing a read swaps its ends.
   This differs from Illumina, where each single minimizer is a bucket.
   Candidates need evidence at both ends, which keeps buckets small at
   endpoint hotspots. Whether a rescue path using single-end seeds is
   needed will be decided by measured recall.
3. **Comparison.** Candidates go through a staged check:
   - **Orientation.** Same strand, or reverse complement when allowed.
   - **Endpoint check.** At each end, the first bases inward from the
     boundary are aligned with the outer boundary anchored. The boundaries
     may be displaced by at most `max_end_offset`, and the inner end of the
     aligned window is free. That window gets its own small error allowance,
     so an indel near an end is tolerated but a different fragment endpoint
     is not.
   - **Whole-read check.** A bounded, banded edit distance over the whole
     interval. A pair whose alignment cannot fit in the band is rejected.
     The band is sized from the error budget.
4. **Clustering.** As in the Illumina path:
   - The first read of a cluster becomes its immutable comparison leader.
   - Later reads join the first leader they match, and clusters never merge.
   - The best-scoring read (quality on the normalized interval) becomes the
     exemplar at finalize.

   Members are verified against the leader, not necessarily against the
   final exemplar.

## Data structures

- **Leader sequences.** Packed at 2 bits per base in a chunked arena, with a
  fallback for non-ACGT bases, and indexed by leader id rather than read
  index. Memory scales with the bases of unique reads. The 64 GB target
  assumes inputs of at most ~100 Gbp; at that size, even fully unique input
  packs to ~25 GB, leaving room for the index.
- **Assignments.** One `u32` assignment per read.
- **Seed index.** A compact map from pair key to leader ids, deduplicated per
  leader.

## Validation

- **Simulator.** Genomes yield molecules with defined endpoints, including
  dense coverage and nearby-but-distinct endpoints. Molecules form PCR
  families of various sizes, and each read gets independent substitution,
  insertion and deletion errors, with optional foldback tails, homopolymer
  tails and reverse complements. All error parameters are supplied at runtime,
  with defaults from published literature.
- **Metrics.** Seed-retrieval recall, comparator recall and precision, and
  duplicate-count bias are measured separately, alongside run time and peak
  memory.
- **Scale test.** The simulator generates a synthetic input at the target
  scale.
