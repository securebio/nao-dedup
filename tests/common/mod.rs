use rand::{Rng, rngs::StdRng};

pub fn random_seq(len: usize, rng: &mut StdRng) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
}
