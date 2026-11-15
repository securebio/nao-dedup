# Agent notes: ONT deduplication feature branch

> Branch-specific working notes for the ONT dedup effort. Trim or remove this
> file before the feature branch is proposed for merge into `main`.

## Goal

Add Oxford Nanopore (ONT) read deduplication to `nao-dedup`, alongside the
existing Illumina paired-end support, as a library API plus a CLI. The end
state is a feature branch that is ready for human review.

The existing Illumina path (`DedupContext`, `dedup_interleaved_fastq`) must
keep its current behavior and performance.

### Why ONT differs from the Illumina path

- Single-end reads of highly variable length (hundreds of bp to tens of kb),
  not fixed-length mate pairs.
- The error profile is indel-heavy, so Hamming distance with a small offset
  (the Illumina comparator) is the wrong similarity model. Comparison needs
  to be indel-aware, for example banded edit distance.
- Reads can come from either strand, so matching must be
  reverse-complement aware.
- Read ends are less reliable: variable adapter/barcode trimming and
  truncated reads.

## Priorities

Good performance > readability > great performance > clever code.

## Decisions so far

- **Scope**: the pipeline's virus-hit read subsets, not whole-run FASTQs.
  Design target: tens of millions of reads per input (reads up to ~15 kb),
  within 64 GB RAM and under 1 hour on a 32-core machine. Single-threaded
  first; add parallelism only if measurements call for it.
- **Duplicate definition (v1)**: same fragment, meaning both endpoints agree
  within a tolerance, with sequence similarity consistent
  with independent ONT errors. Truncated copies that share only one endpoint
  are out of scope for v1. Measure how common they are before considering
  an opt-in containment mode.
- **Reverse-complement tolerance is optional** (an `rc_tolerant` parameter,
  on by default). With non-directional preps, PCR copies of one molecule can
  be sequenced from either strand. Directional preps can turn it off so that
  only same-orientation matches count.
- **Foldbacks are not trimmed** (see the design doc's "Read ends"). A
  virtual-trimming detector was built and evaluated, and is archived outside
  the repo. Revisit only for libraries where foldbacks arise independently
  in each read.
- **Stay close to the Illumina design** where it doesn't compromise the
  goals above: bucket reads (with an ONT-appropriate bucketing strategy),
  compare within buckets, greedy streaming clustering, and resolve exemplars
  at finalize. This keeps behavior predictable and allows code reuse.
- **Validation**: public tests use simulated data with literature or runtime
  parameters. Tuning and sanity checks on real data stay local and never go
  into the repo.
- **Branching**: stacked PRs (`gh stack`), with base `ejf/ont-dedup`
  (pushed by the lead). The lead reads PR descriptions and tracks the overall
  shape of the solution, so each description must convey the design clearly.

## Working rules

- **Roles**: Claude drives the implementation. Codex is a thought partner,
  debugger, adversarial reviewer, and domain-expert sounding board.
  - Nuanced questions and design brainstorming: `gpt-6-astra` at
    medium/high effort (high for design; ask it for diverse options).
  - Well-scoped dispatched tasks: `gpt-6.1-sol` at medium/high effort.
- **Cross-family review**: all code is reviewed by the model family that did
  not write it. Claude-written code goes to Codex (CLI), and vice versa.
- Use `/simplify` and `/review-pr` liberally.
- **Branches and PRs**: PRs target the ONT feature base branch, never `main`.
  Never merge, commit, or push to `main`. The bot can only push
  `coding-agent/*` branches. Keep PRs under ~400 changed lines and stack
  them when the work spans several concerns.
- **Questions**: take project goals and constraints to the human lead. Settle
  technical questions by experiment and research with Codex. High-level
  design questions can go to the lead too.

## Public repo hygiene

This repository is **public**. Never put the following in commits, PR
descriptions, issues, comments, code defaults, or test fixtures:

- secrets or credentials;
- data or metadata read from AWS/S3 (sample names, paths, read content,
  statistics derived from internal data);
- internal details such as our measured sequencing error model.

Simulation and test parameters in public code should come from published
literature, or be supplied at runtime.

## Codebase map

- `src/lib.rs`: Illumina library (`DedupContext`). Minimizer bucketing over
  fixed adjacent windows at the read start, greedy streaming clustering
  against cluster leaders, best-quality exemplar resolved in `finalize()`.
- `src/dedup_interleaved_fastq.rs`: two-pass CLI for gzip interleaved
  FASTQ.
- `src/ont/`: ONT library (`seeds`, `compare`, `leaders` for packed leader
  storage, `context`), and
  `src/dedup_ont_fastq.rs`: the ONT CLI. Tests are in `tests/ont_tests.rs`
  and `tests/ont_cli_tests.rs`; `examples/ont_eval.rs` runs simulation
  sweeps, `examples/ont_nearmiss.rs` measures near-miss pair merge probability,
  and `examples/ont_scale.rs` measures throughput and memory at
  tens of millions of reads.
- `tests/`: integration tests (`dedup_tests.rs`) and CLI tests
  (`cli_tests.rs`).
- `ont-sim/`: dev-only ONT read simulator and dedup evaluator (see
  `ont-sim/README.md`). Use it for all accuracy measurements.
- Build/test: `cargo build --release && cargo test --release`. Format new
  code with `cargo fmt -p <crate>` or rustfmt on new files only; the existing
  Illumina code is not rustfmt-clean, and reformatting it is out of scope.
