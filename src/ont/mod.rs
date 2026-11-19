//! Deduplication of Oxford Nanopore (ONT) single-end reads.
//!
//! See `docs/ont-design.md` for the design. This module is independent of the
//! Illumina path in the crate root.
//!
//! Read ends are taken to be fragment ends: trim adapters, barcodes and primers
//! first. Sequence shared at the ends of many reads defeats the endpoint check.

mod compare;
mod context;
mod leaders;
mod seeds;

pub use compare::{OntDedupParams, validate_quality};
pub use context::{OntDedupContext, OntDedupStats};
pub use seeds::OntMinimizerParams;

/// Reverse complement of an ACGT sequence; other bytes (e.g. N) are kept.
#[cfg(test)]
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    reverse_complement_into(seq, &mut out);
    out
}

/// Write the reverse complement of `seq` into `out`, replacing its contents.
pub(crate) fn reverse_complement_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }));
}

#[cfg(test)]
pub(crate) mod test_util {
    use rand::{Rng, rngs::StdRng};

    pub fn random_seq(len: usize, rng: &mut StdRng) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
    }

    pub fn substitute(seq: &mut [u8], pos: usize) {
        seq[pos] = if seq[pos] == b'A' { b'C' } else { b'A' };
    }
}
