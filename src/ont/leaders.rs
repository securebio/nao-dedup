//! Compact storage for cluster leader sequences.
//!
//! Leaders are the only sequences kept in memory, and at the largest inputs
//! they dominate memory use, so they are packed at 2 bits per base. A leader
//! containing any other byte (e.g. N, which ONT basecallers rarely emit) is
//! kept verbatim instead.

use rustc_hash::FxHashMap;
use std::ops::Range;

const CODE: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table
};
const BASE: [u8; 4] = *b"ACGT";

#[derive(Default)]
pub(crate) struct LeaderStore {
    /// All packed leaders back to back, 4 bases per byte, base `i` of the
    /// stream in bits `2 * (i % 4)` of byte `i / 4`.
    packed: Vec<u8>,
    /// Leader `id` occupies stream bases `starts[id]..starts[id + 1]`.
    starts: Vec<u64>,
    /// Leaders stored verbatim (zero-length in the packed stream).
    verbatim: FxHashMap<u32, Box<[u8]>>,
}

impl LeaderStore {
    /// Store an uppercase sequence; returns its id (0, 1, 2, ...).
    pub(crate) fn push(&mut self, seq: &[u8]) -> u32 {
        if self.starts.is_empty() {
            self.starts.push(0);
        }
        let id = u32::try_from(self.starts.len() - 1).expect("too many leaders");
        let start = *self.starts.last().unwrap();
        if seq.iter().all(|&b| CODE[b as usize] < 4) {
            let end = start + seq.len() as u64;
            self.packed.resize(end.div_ceil(4) as usize, 0);
            for (i, &b) in (start..).zip(seq) {
                self.packed[(i / 4) as usize] |= CODE[b as usize] << (2 * (i % 4));
            }
            self.starts.push(end);
        } else {
            self.verbatim.insert(id, seq.into());
            self.starts.push(start);
        }
        id
    }

    pub(crate) fn len(&self, id: u32) -> usize {
        match self.verbatim.get(&id) {
            Some(seq) => seq.len(),
            None => (self.starts[id as usize + 1] - self.starts[id as usize]) as usize,
        }
    }

    /// Write bases `range` of leader `id` into `out`, replacing its contents.
    pub(crate) fn unpack_into(&self, id: u32, range: Range<usize>, out: &mut Vec<u8>) {
        out.clear();
        if let Some(seq) = self.verbatim.get(&id) {
            out.extend_from_slice(&seq[range]);
            return;
        }
        let start = self.starts[id as usize];
        out.extend((range.start as u64..range.end as u64).map(|i| {
            let i = start + i;
            BASE[((self.packed[(i / 4) as usize] >> (2 * (i % 4))) & 3) as usize]
        }));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ont::test_util::random_seq;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn round_trips_packed_and_verbatim_leaders() {
        let mut rng = StdRng::seed_from_u64(1);
        let mut store = LeaderStore::default();
        let mut seqs = Vec::new();
        for len in [0, 1, 3, 4, 5, 17, 1000] {
            seqs.push(random_seq(len, &mut rng));
        }
        seqs.push(b"ACGTNACGT".to_vec());
        seqs.push(b"GATTACA".to_vec());
        for (i, seq) in seqs.iter().enumerate() {
            assert_eq!(store.push(seq), i as u32);
        }
        let mut out = Vec::new();
        for (i, seq) in seqs.iter().enumerate() {
            let id = i as u32;
            assert_eq!(store.len(id), seq.len());
            store.unpack_into(id, 0..seq.len(), &mut out);
            assert_eq!(&out, seq);
            if seq.len() > 4 {
                store.unpack_into(id, 2..seq.len() - 1, &mut out);
                assert_eq!(out, seq[2..seq.len() - 1]);
            }
        }
        assert!(store.verbatim.contains_key(&7));
    }
}
