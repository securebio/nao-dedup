//! Deduplication of Oxford Nanopore (ONT) single-end reads.
//!
//! See `docs/ont-design.md` for the design. This module is independent of the
//! Illumina path in the crate root.

mod compare;
mod seeds;

pub use compare::OntDedupParams;
pub use seeds::OntMinimizerParams;

/// Reverse complement of an ACGT sequence; other bytes (e.g. N) are kept.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        })
        .collect()
}
