use crate::Truth;
use std::{collections::HashMap, fmt};

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
    pub reads: usize,
    pub true_molecules: usize,
    pub identifiability_classes: usize,
    pub predicted_clusters: usize,
    pub dedup_count_bias: Option<f64>,
    pub molecule: PairMetrics,
    pub identifiability: PairMetrics,
    pub merges: usize,
    pub splits: usize,
    pub foldback: [Recall; 3], // neither, one, both; molecule truth
    pub strand: [Recall; 2],   // same, opposite; molecule truth
}

#[derive(Default)]
struct Histogram {
    n: u64,
    foldback: [u64; 2],
    strand: [u64; 2],
}
impl Histogram {
    fn add(&mut self, t: &Truth) {
        self.n += 1;
        self.foldback[usize::from(t.foldback_len > 0)] += 1;
        self.strand[usize::from(t.is_rc)] += 1;
    }
    fn strata(&self) -> [u64; 5] {
        let [a, b] = self.foldback;
        let [f, r] = self.strand;
        [pairs(a), a * b, pairs(b), pairs(f) + pairs(r), f * r]
    }
}

/// Score unordered pairs in O(reads) expected time using contingency tables.
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
    let (mut molecules, mut cells, mut identities, mut identity_cells, mut predicted) = (
        HashMap::new(),
        HashMap::new(),
        HashMap::new(),
        HashMap::new(),
        HashMap::new(),
    );
    for t in truth {
        let coord = (t.genome_id, t.mol_start, t.mol_end);
        let e = assigned[&t.read_idx];
        let group = molecules
            .entry(t.molecule_id)
            .or_insert((coord, Histogram::default(), 0));
        if coord != group.0 || t.mol_start >= t.mol_end {
            return Err("each molecule must have one consistent, nonempty genomic interval".into());
        }
        group.1.add(t);
        let cell = cells
            .entry((t.molecule_id, e))
            .or_insert_with(Histogram::default);
        let cluster = predicted.entry(e).or_insert((0u64, 0));
        if cell.n == 0 {
            group.2 += 1;
            cluster.1 += 1;
        }
        cell.add(t);
        cluster.0 += 1;
        *identities.entry(coord).or_insert(0u64) += 1;
        *identity_cells.entry((coord, e)).or_insert(0u64) += 1;
    }
    let mut strata = [0; 5];
    let mut recovered = [0; 5];
    for (_, h, _) in molecules.values() {
        for (sum, n) in strata.iter_mut().zip(h.strata()) {
            *sum += n;
        }
    }
    for h in cells.values() {
        for (sum, n) in recovered.iter_mut().zip(h.strata()) {
            *sum += n;
        }
    }
    let predicted_pairs = predicted.values().map(|(n, _)| pairs(*n)).sum();
    Ok(Evaluation {
        reads: truth.len(),
        true_molecules: molecules.len(),
        identifiability_classes: identities.len(),
        predicted_clusters: predicted.len(),
        dedup_count_bias: (!molecules.is_empty())
            .then(|| predicted.len() as f64 / molecules.len() as f64 - 1.0),
        molecule: PairMetrics {
            true_pairs: molecules.values().map(|(_, h, _)| pairs(h.n)).sum(),
            predicted_pairs,
            recovered_pairs: cells.values().map(|h| pairs(h.n)).sum(),
        },
        identifiability: PairMetrics {
            true_pairs: identities.values().map(|n| pairs(*n)).sum(),
            predicted_pairs,
            recovered_pairs: identity_cells.values().map(|n| pairs(*n)).sum(),
        },
        merges: predicted.values().filter(|(_, n)| *n > 1).count(),
        splits: molecules.values().filter(|(_, _, n)| *n > 1).count(),
        foldback: std::array::from_fn(|i| Recall {
            true_pairs: strata[i],
            recovered_pairs: recovered[i],
        }),
        strand: std::array::from_fn(|i| Recall {
            true_pairs: strata[i + 3],
            recovered_pairs: recovered[i + 3],
        }),
    })
}

impl Evaluation {
    /// Numeric JSON only; no external serialization dependency is needed.
    pub fn to_json(&self) -> String {
        format!(
            concat!(
                "{{\"reads\":{},\"true_molecules\":{},\"identifiability_classes\":{},",
                "\"predicted_clusters\":{},\"dedup_count_bias\":{},\"merges\":{},\"splits\":{},",
                "\"molecule\":{},\"identifiability\":{},\"foldback\":{{\"neither\":{},\"one\":{},\"both\":{}}},",
                "\"strand\":{{\"same\":{},\"opposite\":{}}}}}"
            ),
            self.reads,
            self.true_molecules,
            self.identifiability_classes,
            self.predicted_clusters,
            number(self.dedup_count_bias),
            self.merges,
            self.splits,
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
            "reads={} true_molecules={} identifiability_classes={} predicted_clusters={}",
            self.reads, self.true_molecules, self.identifiability_classes, self.predicted_clusters
        )?;
        writeln!(
            f,
            "dedup_count_bias={} merges={} splits={}",
            number(self.dedup_count_bias),
            self.merges,
            self.splits
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
