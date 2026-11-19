use crate::Truth;
use std::{
    collections::{HashMap, HashSet},
    fmt,
};

fn pairs(n: u64) -> u64 {
    n * n.saturating_sub(1) / 2
}
fn ratio(n: u64, d: u64) -> Option<f64> {
    (d > 0).then(|| n as f64 / d as f64)
}
fn number(n: Option<f64>) -> String {
    n.map_or("null".into(), |n| n.to_string())
}

#[derive(Debug)]
pub struct PairMetrics {
    pub true_pairs: u64,
    pub predicted_pairs: u64,
    pub recovered_pairs: u64,
}
impl PairMetrics {
    pub fn precision(&self) -> Option<f64> {
        ratio(self.recovered_pairs, self.predicted_pairs)
    }
    pub fn recall(&self) -> Option<f64> {
        ratio(self.recovered_pairs, self.true_pairs)
    }
    fn json(&self) -> String {
        format!(
            "{{\"precision\":{},\"recall\":{},\"true_pairs\":{},\"predicted_pairs\":{},\"recovered_pairs\":{}}}",
            number(self.precision()),
            number(self.recall()),
            self.true_pairs,
            self.predicted_pairs,
            self.recovered_pairs
        )
    }
}
#[derive(Debug)]
pub struct Evaluation {
    /// True molecules with no read selected as an exemplar.
    pub lost_molecules: usize,
    /// Lost molecules with at least one read assigned across genomes.
    pub lost_molecules_cross_genome: usize,
    /// Observed genomes with no exemplar read.
    pub lost_genomes: usize,
    /// Exemplars beyond one per represented true molecule.
    pub extra_exemplars: usize,
    pub reads: usize,
    pub true_molecules: usize,
    pub identifiability_classes: usize,
    pub predicted_clusters: usize,
    pub molecule: PairMetrics,
    pub identifiability: PairMetrics,
}

/// Score evidence loss and unordered pairs in O(reads) expected time.
/// Indices may be sparse; row order is irrelevant. Undefined ratios are None.
pub fn evaluate(truth: &[Truth], assignments: &[(usize, usize)]) -> Result<Evaluation, String> {
    let rows: HashMap<_, _> = truth.iter().map(|t| (t.read_idx, t)).collect();
    let assigned: HashMap<_, _> = assignments.iter().copied().collect();
    if rows.len() != truth.len()
        || assigned.len() != assignments.len()
        || assigned.len() != rows.len()
        || assigned
            .iter()
            .any(|(i, e)| !rows.contains_key(i) || assigned.get(e) != Some(e))
    {
        return Err("truth indices must be unique; assign every truth read exactly once to a self-mapped exemplar".into());
    }

    type Coord = (usize, usize, usize); // genome, start, end
    // molecule_id -> (interval, reads in the molecule)
    let mut molecules: HashMap<usize, (Coord, u64)> = HashMap::new();
    // (molecule_id, exemplar) -> reads of that molecule in that cluster
    let mut cells: HashMap<(usize, usize), u64> = HashMap::new();
    // interval -> reads with that interval (identifiability classes)
    let mut identities: HashMap<Coord, u64> = HashMap::new();
    // (interval, exemplar) -> reads of that class in that cluster
    let mut identity_cells: HashMap<(Coord, usize), u64> = HashMap::new();
    // exemplar -> cluster size
    let mut predicted: HashMap<usize, u64> = HashMap::new();
    let mut represented_molecules = HashSet::new();
    let mut cross_genome_molecules = HashSet::new();
    let mut genomes: HashMap<usize, bool> = HashMap::new();
    for t in truth {
        let coord = (t.genome_id, t.mol_start, t.mol_end);
        let e = assigned[&t.read_idx];
        let represented_genome = genomes.entry(t.genome_id).or_default();
        if e == t.read_idx {
            represented_molecules.insert(t.molecule_id);
            *represented_genome = true;
        }
        if rows[&e].genome_id != t.genome_id {
            cross_genome_molecules.insert(t.molecule_id);
        }
        let molecule = molecules.entry(t.molecule_id).or_insert((coord, 0));
        if coord != molecule.0 || t.mol_start >= t.mol_end {
            return Err("each molecule must have one consistent, nonempty genomic interval".into());
        }
        molecule.1 += 1;
        *cells.entry((t.molecule_id, e)).or_default() += 1;
        *predicted.entry(e).or_default() += 1;
        *identities.entry(coord).or_default() += 1;
        *identity_cells.entry((coord, e)).or_default() += 1;
    }

    let predicted_pairs = predicted.values().map(|n| pairs(*n)).sum();
    Ok(Evaluation {
        lost_molecules: molecules.len() - represented_molecules.len(),
        lost_molecules_cross_genome: cross_genome_molecules
            .difference(&represented_molecules)
            .count(),
        lost_genomes: genomes
            .values()
            .filter(|represented| !**represented)
            .count(),
        extra_exemplars: predicted.len() - represented_molecules.len(),
        reads: truth.len(),
        true_molecules: molecules.len(),
        identifiability_classes: identities.len(),
        predicted_clusters: predicted.len(),
        molecule: PairMetrics {
            true_pairs: molecules.values().map(|(_, n)| pairs(*n)).sum(),
            predicted_pairs,
            recovered_pairs: cells.values().map(|n| pairs(*n)).sum(),
        },
        identifiability: PairMetrics {
            true_pairs: identities.values().map(|n| pairs(*n)).sum(),
            predicted_pairs,
            recovered_pairs: identity_cells.values().map(|n| pairs(*n)).sum(),
        },
    })
}

impl Evaluation {
    pub fn lost_molecule_rate(&self) -> Option<f64> {
        ratio(self.lost_molecules as u64, self.true_molecules as u64)
    }

    pub fn missed_duplicate_rate(&self) -> Option<f64> {
        ratio(
            self.extra_exemplars as u64,
            (self.reads - self.true_molecules) as u64,
        )
    }

    /// Numeric JSON only; no external serialization dependency is needed.
    pub fn to_json(&self) -> String {
        format!(
            concat!(
                "{{\"lost_molecules\":{},\"lost_molecules_cross_genome\":{},\"lost_genomes\":{},",
                "\"extra_exemplars\":{},\"lost_molecule_rate\":{},\"missed_duplicate_rate\":{},",
                "\"reads\":{},\"true_molecules\":{},\"identifiability_classes\":{},",
                "\"predicted_clusters\":{},\"molecule\":{},\"identifiability\":{}}}"
            ),
            self.lost_molecules,
            self.lost_molecules_cross_genome,
            self.lost_genomes,
            self.extra_exemplars,
            number(self.lost_molecule_rate()),
            number(self.missed_duplicate_rate()),
            self.reads,
            self.true_molecules,
            self.identifiability_classes,
            self.predicted_clusters,
            self.molecule.json(),
            self.identifiability.json()
        )
    }
}

impl fmt::Display for Evaluation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "lost_molecules={} lost_molecule_rate={} lost_molecules_cross_genome={} lost_genomes={}",
            self.lost_molecules,
            number(self.lost_molecule_rate()),
            self.lost_molecules_cross_genome,
            self.lost_genomes
        )?;
        writeln!(
            f,
            "extra_exemplars={} missed_duplicate_rate={}",
            self.extra_exemplars,
            number(self.missed_duplicate_rate())
        )?;
        writeln!(
            f,
            "reads={} true_molecules={} identifiability_classes={} predicted_clusters={}",
            self.reads, self.true_molecules, self.identifiability_classes, self.predicted_clusters
        )?;
        for (name, m) in [
            ("molecule", &self.molecule),
            ("identifiability", &self.identifiability),
        ] {
            writeln!(
                f,
                "{name}: precision={} recall={} (recovered={} predicted={} true={} pairs)",
                number(m.precision()),
                number(m.recall()),
                m.recovered_pairs,
                m.predicted_pairs,
                m.true_pairs
            )?;
        }
        Ok(())
    }
}
