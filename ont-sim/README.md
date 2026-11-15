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

Pair precision/recall use unordered, distinct read pairs and contingency counts.
Both molecule truth and identifiability truth (equal `genome_id,start,end`) are
reported. Merges, splits, count bias, and foldback/strand recall use molecule
truth. Artifact strata describe pairs with neither/one/both nonzero foldbacks,
and same/opposite strand. Missing denominators yield `None` in Rust and `null`
in text/JSON; pair counts accompany metrics. No family-size weighting is applied
beyond the natural number of pairs per family. Header-only inputs are valid.
