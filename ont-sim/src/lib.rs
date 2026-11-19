//! Deterministic synthetic ONT reads and evidence/pair evaluation; see ../README.md.
mod evaluate;
mod simulate;
pub use evaluate::*;
pub use simulate::*;

#[derive(Clone, Debug, PartialEq)]
pub struct Truth {
    pub read_idx: usize,
    pub molecule_id: usize,
    pub genome_id: usize,
    pub mol_start: usize,
    pub mol_end: usize,
    pub is_rc: bool,
    pub n_sub: usize,
    pub n_ins: usize,
    pub n_del: usize,
}

#[cfg(test)]
mod tests;
