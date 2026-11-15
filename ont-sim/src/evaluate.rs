use crate::Truth;
use std::{
    collections::{HashMap, HashSet},
    fmt,
    hash::Hash,
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
pub struct Recall {
    pub true_pairs: u64,
    pub recovered_pairs: u64,
}
impl Recall {
    pub fn recall(&self) -> Option<f64> {
        ratio(self.recovered_pairs, self.true_pairs)
    }
    fn json(&self) -> String {
        format!(
            "{{\"recall\":{},\"true_pairs\":{},\"recovered_pairs\":{}}}",
            number(self.recall()),
            self.true_pairs,
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
    pub dedup_count_bias: Option<f64>,
    pub molecule: PairMetrics,
    pub identifiability: PairMetrics,
    /// Predicted clusters containing more than one true molecule.
    pub merges: usize,
    /// True molecules spread over more than one predicted cluster.
    pub splits: usize,
    /// As `merges`, but molecules with identical endpoints count as one class,
    /// so merges no endpoint-based method could avoid are excluded.
    pub identifiability_merges: usize,
    /// As `splits`, under identifiability truth.
    pub identifiability_splits: usize,
    pub foldback: [Recall; 3], // neither, one, both; molecule truth
    pub strand: [Recall; 2],   // same, opposite; molecule truth
}

/// Unordered read-pair counts in each artifact stratum.
#[derive(Clone, Copy, Default)]
struct Strata {
    foldback_neither: u64,
    foldback_one: u64,
    foldback_both: u64,
    strand_same: u64,
    strand_opposite: u64,
}

impl std::ops::AddAssign for Strata {
    fn add_assign(&mut self, o: Self) {
        self.foldback_neither += o.foldback_neither;
        self.foldback_one += o.foldback_one;
        self.foldback_both += o.foldback_both;
        self.strand_same += o.strand_same;
        self.strand_opposite += o.strand_opposite;
    }
}

/// Read counts for one group of reads (a molecule, or a molecule/cluster cell).
#[derive(Default)]
struct Histogram {
    n: u64,
    with_foldback: u64,
    rc: u64,
}

impl Histogram {
    fn add(&mut self, t: &Truth) {
        self.n += 1;
        self.with_foldback += u64::from(t.foldback_len > 0);
        self.rc += u64::from(t.is_rc);
    }

    /// Pair counts within this group, by stratum.
    fn strata(&self) -> Strata {
        let (fold, plain) = (self.with_foldback, self.n - self.with_foldback);
        let (rc, fwd) = (self.rc, self.n - self.rc);
        Strata {
            foldback_neither: pairs(plain),
            foldback_one: plain * fold,
            foldback_both: pairs(fold),
            strand_same: pairs(fwd) + pairs(rc),
            strand_opposite: fwd * rc,
        }
    }
}

/// Count (merges, splits) from the distinct (truth class, predicted cluster)
/// cells: clusters holding more than one class, and classes spread over more
/// than one cluster.
fn merges_and_splits<T: Hash + Eq, C: Hash + Eq>(
    cells: impl Iterator<Item = (T, C)>,
) -> (usize, usize) {
    let mut classes_per_cluster: HashMap<C, usize> = HashMap::new();
    let mut clusters_per_class: HashMap<T, usize> = HashMap::new();
    for (class, cluster) in cells {
        *classes_per_cluster.entry(cluster).or_default() += 1;
        *clusters_per_class.entry(class).or_default() += 1;
    }
    let count = |n: &usize| *n > 1;
    (
        classes_per_cluster.values().filter(|n| count(n)).count(),
        clusters_per_class.values().filter(|n| count(n)).count(),
    )
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
    let mut molecules: HashMap<usize, (Coord, Histogram)> = HashMap::new();
    // (molecule_id, exemplar) -> reads of that molecule in that cluster
    let mut cells: HashMap<(usize, usize), Histogram> = HashMap::new();
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
        let molecule = molecules
            .entry(t.molecule_id)
            .or_insert_with(|| (coord, Histogram::default()));
        if coord != molecule.0 || t.mol_start >= t.mol_end {
            return Err("each molecule must have one consistent, nonempty genomic interval".into());
        }
        molecule.1.add(t);
        cells.entry((t.molecule_id, e)).or_default().add(t);
        *predicted.entry(e).or_default() += 1;
        *identities.entry(coord).or_default() += 1;
        *identity_cells.entry((coord, e)).or_default() += 1;
    }

    let mut strata = Strata::default();
    for (_, h) in molecules.values() {
        strata += h.strata();
    }
    let mut recovered = Strata::default();
    for h in cells.values() {
        recovered += h.strata();
    }
    let recall = |t: u64, r: u64| Recall {
        true_pairs: t,
        recovered_pairs: r,
    };
    let (merges, splits) = merges_and_splits(cells.keys().copied());
    let (identifiability_merges, identifiability_splits) =
        merges_and_splits(identity_cells.keys().copied());
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
        dedup_count_bias: (!molecules.is_empty())
            .then(|| predicted.len() as f64 / molecules.len() as f64 - 1.0),
        molecule: PairMetrics {
            true_pairs: molecules.values().map(|(_, h)| pairs(h.n)).sum(),
            predicted_pairs,
            recovered_pairs: cells.values().map(|h| pairs(h.n)).sum(),
        },
        identifiability: PairMetrics {
            true_pairs: identities.values().map(|n| pairs(*n)).sum(),
            predicted_pairs,
            recovered_pairs: identity_cells.values().map(|n| pairs(*n)).sum(),
        },
        merges,
        splits,
        identifiability_merges,
        identifiability_splits,
        foldback: [
            recall(strata.foldback_neither, recovered.foldback_neither),
            recall(strata.foldback_one, recovered.foldback_one),
            recall(strata.foldback_both, recovered.foldback_both),
        ],
        strand: [
            recall(strata.strand_same, recovered.strand_same),
            recall(strata.strand_opposite, recovered.strand_opposite),
        ],
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
                "\"predicted_clusters\":{},\"dedup_count_bias\":{},\"merges\":{},\"splits\":{},",
                "\"identifiability_merges\":{},\"identifiability_splits\":{},",
                "\"molecule\":{},\"identifiability\":{},\"foldback\":{{\"neither\":{},\"one\":{},\"both\":{}}},",
                "\"strand\":{{\"same\":{},\"opposite\":{}}}}}"
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
            number(self.dedup_count_bias),
            self.merges,
            self.splits,
            self.identifiability_merges,
            self.identifiability_splits,
            self.molecule.json(),
            self.identifiability.json(),
            self.foldback[0].json(),
            self.foldback[1].json(),
            self.foldback[2].json(),
            self.strand[0].json(),
            self.strand[1].json()
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
        writeln!(
            f,
            "dedup_count_bias={} merges={} splits={} identifiability_merges={} identifiability_splits={}",
            number(self.dedup_count_bias),
            self.merges,
            self.splits,
            self.identifiability_merges,
            self.identifiability_splits
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
        for (name, s) in [
            "foldback_neither",
            "foldback_one",
            "foldback_both",
            "strand_same",
            "strand_opposite",
        ]
        .into_iter()
        .zip(self.foldback.iter().chain(&self.strand))
        {
            writeln!(
                f,
                "{name}: recall={} ({}/{} pairs)",
                number(s.recall()),
                s.recovered_pairs,
                s.true_pairs
            )?;
        }
        Ok(())
    }
}
