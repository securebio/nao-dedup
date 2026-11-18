# ont-sim

Dev-only, in-memory simulator and evaluator for ONT deduplication experiments.
This crate is independent of `nao_dedup`, unpublished, and available to its
integration tests as a dev dependency. Use `simulate(&SimConfig)` and
`evaluate(&[Truth], &[(read_idx, exemplar_idx)])` from Rust.

```sh
cargo run --release -p ont-sim --bin ont_simulate -- \
  --seed 42 --n-genomes 2 --genome-len 100000 --n-molecules 1000 \
  --len-median 1000 --len-sigma 0.5 --len-min 100 --len-max 15000 \
  --family-dist geometric --family-mean 3 \
  --endpoints hotspot --n-hotspots 20 --hotspot-jitter 3 \
  --sub-rate 0.01 --ins-rate 0.01 --del-rate 0.01 \
  --fastq reads.fastq.gz --truth truth.tsv

# Supply assignments from a deduplicator, including singleton exemplars.
cargo run --release -p ont-sim --bin ont_evaluate -- \
  --truth truth.tsv --assignments assignments.tsv --json metrics.json
```

The 1%/1%/1% rates above are **illustrative**, not an instrument error profile.
All three sequencing error rates are required. `--help` lists the remaining
parameters. The same seed and parameters reproduce output within the same
build/dependency versions, including shuffle order. All genomes and reads are
held in memory; this is intended for tuning datasets, not streaming at scale.

Truth and model semantics:

- Genomes are independent random sequences by default. With
  `--genome-divergence D > 0`, genome 0 is random and each other genome copies
  it with independent per-base substitutions of probability D, always to a
  different base. D must be finite and in [0,1]. Coordinates remain comparable
  across genomes. For small D, pairwise divergence between two non-reference
  genomes is about twice D (exactly `2D - 4D²/3` in expectation). D=0 draws
  independent genomes.
- Each `molecule_id` is a separately drawn molecule, even when another molecule
  has identical endpoints. Coordinates are zero-based, half-open `[start,end)`
  on the forward genome, before PCR mutations, strand choice, and artifacts.
- Uniform placement draws a rounded log-normal length, clips it to the bounds
  and genome length, then draws a valid start uniformly. Hotspot placement
  first draws `n_hotspots` **paired start/end anchors per genome** using that
  model; each molecule chooses an anchor uniformly and jitters both endpoints
  independently by a uniform integer in `[-J,J]`. Bounds are enforced again.
  Zero jitter intentionally yields identical-endpoint hard negatives.
- Family size is geometric on `1,2,...`, with the requested expectation. With
  `pcr_error_prob`, one position and different base are chosen per family;
  exactly `floor(k/2)` randomly selected reads inherit this shared substitution.
  These remain one true molecule. Singletons inherit no PCR substitution.
- Reverse complementation precedes foldback, poly-A, end jitter, and errors.
  Foldback appends the reverse complement of the terminal `floor(length*fraction)`
  bases; poly-A follows it. Each end is then independently trimmed or extended
  with equal probability by `0..=end_jitter` random bases (trim stops at empty).
- Insertions (one random base before a source base) and deletions are independent
  Bernoulli events per source base. Substitution probability is conditional on
  retention and always changes the base. Runs of at least three identical bases
  multiply both indel probabilities, capped at 1. Runs include added artifacts.
  The error counts describe sequencing errors only; a deleted base cannot also
  count as a substitution. All-deleted reads may be empty.
  The public `add_errors(seq, sub, ins, del, homopolymer_indel_mult, rng)`
  function shares this model and returns the sequence and `[sub, ins, del]`
  event counts; callers supply the probabilities and RNG.
- `is_rc` records strand choice; `foldback_len` and `polya_len` record appended
  lengths **before** jitter and sequencing errors, not surviving tail lengths.
  Mean quality is drawn as an integer, with independent uniform ±3 noise per
  emitted base, clamped to Q2–Q50. `mean_q` is the actual emitted mean (0 if empty).
- Reads are shuffled before assigning `read_idx` and `read_name` (`r{idx}`). These
  indices, also used for exemplar indices, refer to FASTQ output order.

Truth TSV has the exact header exported as `TRUTH_HEADER`; `is_rc` is `true` or
`false`. Assignments require `read_idx<TAB>exemplar_idx` and exactly one row for
every truth index. Rows may be unordered; indices may be sparse. Every exemplar
must occur in truth and map to itself. Headers and column counts are strict.

Evidence-level metrics are the headline in text and are also available in JSON:

- `lost_molecules`: true molecules none of whose reads is an exemplar (a read
  assigned to itself). This counts evidence lost to false merges.
- `lost_molecules_cross_genome`: lost molecules with at least one read assigned
  to an exemplar from a different `genome_id`.
- `lost_genomes`: observed genome IDs with no exemplar read.
- `extra_exemplars`: predicted clusters minus represented true molecules,
  where represented molecules = `true_molecules - lost_molecules`. This
  nonnegative count measures missed duplicates beyond one exemplar per molecule.
- `lost_molecule_rate` divides lost molecules by true molecules;
  `missed_duplicate_rate` divides extra exemplars by `reads - true_molecules`.
  A molecule split across clusters is still represented if any of its reads is
  an exemplar. Counts give each molecule equal weight regardless of family size.

The existing pair precision/recall use unordered, distinct read pairs and contingency counts.
Both molecule truth and identifiability truth (equal `genome_id,start,end`) are
reported. Merges and splits are reported under both truths; count bias and
foldback/strand recall use molecule truth. Artifact strata describe pairs with neither/one/both nonzero foldbacks,
and same/opposite strand. Missing denominators yield `None` in Rust and `null`
in text/JSON; pair counts accompany metrics. No family-size weighting is applied
beyond the natural number of pairs per family. Header-only inputs are valid.

The root example `examples/ont_nearmiss.rs` measures isolated-pair merge
probability using the full deduplicator, including seed retrieval:

```sh
cargo run --release --example ont_nearmiss -- \
  --sub-rate <P> --ins-rate <P> --del-rate <P>
```

It sweeps `--divergences`, `--offsets`, and `--modes one,both`, with `--trials`
independent pairs per cell and a reproducible `--seed`. Template lengths are
uniform in `--len-min..=--len-max`. Variants are substitutions or single-base
indels (`--variant-indel-frac`), with endpoint offsets measured in source
coordinates via a boundary map. Insertions belong to the following source base;
deletions collapse adjacent boundaries. Mode `one` shifts one random end, and
mode `both` shifts each end independently, with random signs. Lengths must
exceed the largest inward shift (twice the maximum offset for mode `both`), so
endpoints cannot cross. Sequencing errors use `add_errors`; each read then has
an independent 50% reverse-complement probability; with `--no-rc` both reads of
a pair share one random strand, so the control still measures recall.

TSV stdout reports merge counts, rates, and the upper endpoint of the two-sided
95% Wilson score interval. Divergence=0, offset=0 is a duplicate control: its
merge rate measures true-pair recall. Other cells measure false merges between
distinct molecules. Deduplication/minimizer flags match `ont_eval`.

### Quality-conditioned errors (opt in)

Add `--quality-conditioned --error-median <P> --error-sigma <S>
--local-q-sigma <Q> --miscalibration <M> --latent-sigma <S>` to `ont_simulate`
or `ont_eval`. All quality-model distribution inputs are required. Without this
flag, error rates set event probabilities directly.

In this mode `--sub-rate`, `--ins-rate`, and `--del-rate` are relative event
weights, normalized to sum to one. A read's nominal probability is lognormal
with the requested median and log standard deviation (clipped to one). Local
normal variation in Phred units is rounded and clipped to Q0–Q93 **before**
drawing errors. Each quality-derived probability is multiplied by the batch
miscalibration factor, then `--unreported-error-floor` (default 0) is added.
The sum is multiplied by an independent mean-one lognormal read multiplier
(log standard deviation `--latent-sigma`) and capped at one. Zero latent sigma
removes residual read dispersion. These multipliers affect errors, not qualities.
Homopolymer enrichment changes the local event weights and nominal qualities.
Local variation can occasionally round to constant qualities on very short reads.

Events are categorical: a substitution changes the base, an insertion emits a
random base before the retained source base, and a deletion emits nothing.
Insertions share the pre-sampled source quality; deletions remove its quality.
`add_quality_errors(seq, &QualityConfig, weights, homopolymer_mult, rng)` exposes
this model publicly and returns `(sequence, quality, QualityErrorStats)`.
The stats contain `[sub, ins, del]` event counts and the number of source positions
whose event probability exceeds one; `Simulation` sums the saturation counts
and `ont_simulate` reports them on stderr.
`QualityConfig::enabled(median, sigma, local_q_sigma, miscalibration, latent_sigma)`
constructs an enabled model with a zero unreported error floor.
Call `QualityConfig::validate` first. Qualities always match emitted sequence
lengths. Emitted-quality expected errors approximate actual event counts at low
error rates; indel selection and changed lengths introduce second-order bias.
Calibration of read-to-read edit distance needs separate validation because
coincident errors and alternative alignments can mask events.

Tests check determinism, seeded output of the independent-error model, substitution-only
and indel calibration against **emitted** qualities, homopolymer enrichment, and
batch scaling. This simple simulator does not model spatially correlated error
bursts or systematic sequence-specific substitutions; those remain validation
work, along with held-out protected strata, saturated groups, and broader length
and accuracy sweeps. The near-miss report is a complete-pipeline check, not a
comparator-only estimate or a guarantee for every stratum.
