//! Deterministic synthetic ONT reads and pair-count evaluation; see ../README.md.
mod evaluate;
mod io;
mod simulate;
pub use evaluate::*;
pub use io::*;
pub use simulate::*;

#[derive(Clone, Debug, PartialEq)]
pub struct Truth {
    pub read_idx: usize,
    pub read_name: String,
    pub molecule_id: usize,
    pub genome_id: usize,
    pub mol_start: usize,
    pub mol_end: usize,
    pub is_rc: bool,
    pub foldback_len: usize,
    pub polya_len: usize,
    pub n_sub: usize,
    pub n_ins: usize,
    pub n_del: usize,
    pub mean_q: f64,
}

#[cfg(test)]
mod tests;
