# ONT deduplication: design

Status: describes the implemented `ont` module (`src/ont/`) and the
`dedup_ont_fastq` CLI.

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

## Algorithm

The structure is the same as the Illumina path: bucket reads by minimizers,
compare each new read only against cluster leaders that share a bucket,
cluster greedily, and resolve the best exemplar at finalize.

1. **Seeds.** Illumina takes minimizers from adjacent windows at the start of
   each mate. ONT treats the read's two ends like two mates. It takes one
   minimizer from each of `num_windows` adjacent windows of `window_len`
   bases, counted inward from both the start and the end of the read.
   - K-mers are canonical, so each read and its reverse complement give the
     same seeds.
   - Minimizers are chosen by a mixing hash, not the lexicographic 2-bit
     value, so poly-A k-mers aren't systematically preferred.
2. **Buckets.** As in the Illumina path, each seed is a bucket key: a
   leader is posted under each of its end seeds, with the end each came
   from and its distance from that end. There are up to 16 seeds with the
   defaults; fewer for short reads or reads with non-ACGT or homopolymer
   runs. A leader with no seed at one end gets no postings and stays a
   singleton. Unlike Illumina, a candidate needs
   evidence at *both* ends. A read's start seeds must hit one end of the
   leader and its end seeds the other, in a consistent orientation (start
   to start and end to end, or swapped for the reverse complement). Only
   then is it aligned. The evidence picks which orientation to try first;
   if that fails, the other allowed orientation is tried too, because in
   terminal inverted repeats seeds can point the wrong way.
   - **Offset filter.** Postings whose offsets differ from the read's by
     more than the endpoint allowance (plus indel drift proportional to the
     distance from the end) are ignored, since such a leader can't pass the
     endpoint check. For non-repetitive sequence this is exact within the
     end-check window; further in, and in repeats, it is a heuristic.
     Buckets are sorted by offset, so a lookup visits only the compatible
     range.
   - **Search order.** Seeds are consulted from the read ends inward, and
     the search stops at the first matching leader. This matters in dense
     groups, where inner seeds can have very large compatible ranges.
3. **Comparison.** Candidates go through a staged check:
   - **Orientation.** Same strand, or reverse complement when allowed.
   - **Length check.** The length difference can't exceed the whole-read
     edit budget.
   - **Endpoint check.** At each end, the first `end_check_len` bases
     inward from the boundary are aligned with the outer boundary anchored
     and the inner end free. Up to `max_end_offset` bases of displacement
     are free, and at most `max_end_errors` edits are allowed. Edit distance
     alone can't tell a displaced endpoint from indels at the boundary, so
     an end may be aligned twice. The first alignment gets `max_end_offset`
     bases of free displacement; displacement beyond that costs one edit
     per base. The second gets `max_end_offset + max_end_errors` bases of
     free displacement, and is only needed when the first costs something
     but fits the allowance. If the second is strictly cheaper, the
     endpoints differ and the pair is rejected. This
     keeps the error allowance from silently widening the endpoint
     tolerance. Low-complexity ends are inherently ambiguous under this
     rule.
   - **Whole-read check.** Edit distance over the whole read, bounded by
     `max_error_frac` of the longer read, computed with a diagonal-transition
     algorithm (the `editdistancek` crate) whose cost grows with the actual
     distance.
4. **Clustering.** As in the Illumina path:
   - The first read of a cluster becomes its immutable comparison leader.
   - Later reads join the first leader they match, and clusters never merge.
   - The best-scoring read (mean quality, then length) becomes the
     exemplar at finalize.

   Members are verified against the leader, not necessarily against the
   final exemplar.

## Data structures

- **Leader sequences.** Packed at 2 bits per base, with a verbatim fallback
  for the rare leader containing non-ACGT bases, and indexed by leader id
  rather than read index. A candidate's end windows are unpacked first, and
  the whole leader only if its ends match.
- **Assignments.** One `u32` assignment per read.
- **Seed index.** A hash map from seed to an offset-sorted list of 8-byte
  postings (leader id, offset, end), up to 16 postings per leader. Most seeds
  belong to one leader, so a bucket holds one posting inline. Insertion
  keeps each bucket sorted, so it is linear in the bucket's size, which
  only matters for seeds shared by very many leaders.
- **Measured footprint.** In synthetic tests with reads of 0.5-4.5 kb
  (mean ~2.5 kb), peak memory was about 0.9 KB per leader when leaders
  share seeds and about 1.4 KB when every seed is distinct. Packed sequence
  accounts for a quarter byte per base, so longer reads cost
  proportionally more. With tens of millions of reads, the distinct genome
  sequence present is usually far smaller than the bases sequenced, so
  seeds are shared. At that read length, even 50M unique reads then fit
  in 64 GB.
  `examples/ont_scale.rs` streams synthetic reads to measure throughput and
  peak memory at that scale.

## Validation

- **Simulator.** Genomes yield molecules with defined endpoints, including
  dense coverage and nearby-but-distinct endpoints. Molecules form PCR
  families of various sizes, and each read gets independent substitution,
  insertion and deletion errors, with optional foldback tails, homopolymer
  tails and reverse complements. All error rates are supplied at runtime;
  nothing in the code encodes a particular instrument's error profile.
- **Metrics.** The evaluator reports pair-level precision and recall
  (against both molecule truth and endpoint-identifiability truth), merges,
  splits, duplicate-count bias, and recall by strand and foldback status.
  `ont_eval` adds run time and work counters. Seed retrieval and the
  comparator were also measured separately during development, with a
  local diagnostic.
- **Scale test.** `examples/ont_scale.rs` streams synthetic reads, never
  storing them, to measure throughput and peak memory at tens of millions
  of reads.
- **Real data.** Parameters were also checked locally against real data,
  using reads with identical aligned endpoints as a stand-in for truth.
  That data is not part of this repository.

## Limitations

- **Saturated groups.** When very many reads cover a short genome,
  independent molecules often share endpoints within the tolerance, so
  counts in such groups are approximate in both directions. These groups
  are also slower per read: each seed matches many leaders, and the cost
  grows with the number of leaders. For example, with all reads from one
  7 kb genome (`ont_scale --genomes 1 --genome-len 7000`, default
  parameters), per-read cost on our test machine went from about 43 to 74
  µs as reads went from 1M to 4M. Typical inputs ran at 19-33 µs per read.
- **Accuracy regime.** The defaults are a design assumption: reads of
  roughly 1% per-read error, typical of current high-accuracy ONT
  basecalling. Noisier reads need larger `max_error_frac` and `max_end_errors`.
- **Truncated copies.** Copies that share only one endpoint are not merged.
- **Foldbacks** are not trimmed (see "Read ends").
