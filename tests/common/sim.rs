//! Small deterministic generator for integration tests and the scale example.
//! All error parameters are explicit illustrative inputs, not instrument profiles.
//! Dataset generation stores reads; the benchmark calls `add_quality_errors`
//! one read at a time instead. All randomness comes from the caller's `StdRng`
//! (or `SimConfig::seed`); hash iteration never influences generated reads.

use rand::{Rng, SeedableRng, distributions::Open01, rngs::StdRng, seq::SliceRandom};
use std::{collections::HashMap, hash::Hash};

/// Quality-first model: no homopolymer enrichment or unreported-error floor.
#[derive(Clone, Copy, Debug)]
pub struct QualityConfig {
    /// Median nominal per-read event probability, clipped to one.
    pub error_median: f64,
    /// Log standard deviation of nominal error.
    pub error_sigma: f64,
    /// Standard deviation of per-base Phred noise before rounding/clamping.
    pub local_q_sigma: f64,
    /// Batch multiplier on quality-derived error probabilities.
    pub miscalibration: f64,
    /// Log standard deviation of the mean-one latent read multiplier.
    pub latent_sigma: f64,
}

impl QualityConfig {
    pub fn validate(&self, weights: [f64; 3]) -> Result<(), String> {
        if !(self.error_median > 0.0
            && self.error_median <= 1.0
            && self.local_q_sigma.is_finite()
            && self.local_q_sigma > 0.0)
            || [self.error_sigma, self.miscalibration, self.latent_sigma]
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
            || weights.iter().any(|x| !x.is_finite() || *x < 0.0)
            || !weights.iter().sum::<f64>().is_finite()
            || weights.iter().sum::<f64>() <= 0.0
        {
            return Err("invalid quality model or event weights".into());
        }
        Ok(())
    }
}

fn normal(rng: &mut StdRng) -> f64 {
    (-2.0 * rng.sample::<f64, _>(Open01).ln()).sqrt()
        * (std::f64::consts::TAU * rng.gen_range(0.0..1.0)).cos()
}

fn base(rng: &mut StdRng) -> u8 {
    b"ACGT"[rng.gen_range(0..4)]
}

fn substitution(b: u8, rng: &mut StdRng) -> u8 {
    loop {
        let other = base(rng);
        if other != b {
            return other;
        }
    }
}

/// Reverse complement of uppercase ACGT; other symbols are unsupported.
pub fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("expected uppercase ACGT"),
        })
        .collect()
}

/// Draw a lognormal nominal error per read (clipped to one), then local Phred
/// noise, rounded/clamped to Q0–Q93. Actual categorical event probabilities are
/// quality-derived, times batch miscalibration and a mean-one lognormal latent
/// read multiplier, capped at one. Weights specify sub/ins/del proportions.
/// Insertions precede their retained source base and share its quality;
/// deletions drop both base and quality. Returns emitted bases, Phred+33, counts.
pub fn add_quality_errors(
    s: &[u8],
    c: &QualityConfig,
    weights: [f64; 3],
    rng: &mut StdRng,
) -> (Vec<u8>, Vec<u8>, [usize; 3]) {
    c.validate(weights).expect("invalid quality model");
    let nominal_log = (c.error_median.ln() + c.error_sigma * normal(rng)).min(0.0);
    let latent_log = c.latent_sigma * (normal(rng) - 0.5 * c.latent_sigma);
    let by_q: [f64; 94] = std::array::from_fn(|q| {
        (-(q as f64) * std::f64::consts::LN_10 / 10.0 + c.miscalibration.ln() + latent_log)
            .exp()
            .min(1.0)
    });
    let total = weights.iter().sum::<f64>();
    let weights = weights.map(|w| w / total);
    let (mut out, mut qual, mut counts) = (
        Vec::with_capacity(s.len()),
        Vec::with_capacity(s.len()),
        [0; 3],
    );
    for &b in s {
        let q = (-10.0 / std::f64::consts::LN_10 * nominal_log + c.local_q_sigma * normal(rng))
            .round()
            .clamp(0.0, 93.0) as u8;
        if !rng.gen_bool(by_q[usize::from(q)]) {
            out.push(b);
            qual.push(q + 33);
            continue;
        }
        let event = rng.gen_range(0.0..1.0);
        if event < weights[0] {
            out.push(substitution(b, rng));
            qual.push(q + 33);
            counts[0] += 1;
        } else if event < weights[0] + weights[1] {
            out.extend([base(rng), b]);
            qual.extend([q + 33, q + 33]);
            counts[1] += 1;
        } else {
            counts[2] += 1;
        }
    }
    (out, qual, counts)
}

#[derive(Clone, Copy, Debug)]
pub enum Errors {
    /// Per-base sub/ins/del rates: independent insertion/deletion trials,
    /// substitution only if the source base is retained. Qualities are
    /// independent of errors: one integer mean per read in Q10–Q30, plus
    /// uniform [-3,3] noise per emitted base, clamped to Q2–Q50.
    Simple { rates: [f64; 3] },
    /// Quality-first errors; see `add_quality_errors`.
    Quality {
        model: QualityConfig,
        weights: [f64; 3],
    },
}

fn sequence_read(s: &[u8], errors: Errors, rng: &mut StdRng) -> (Vec<u8>, Vec<u8>) {
    let [sub, ins, del] = match errors {
        Errors::Quality { model, weights } => {
            let (seq, qual, _) = add_quality_errors(s, &model, weights, rng);
            return (seq, qual);
        }
        Errors::Simple { rates } => rates,
    };
    let mut seq = Vec::with_capacity(s.len());
    for &b in s {
        if rng.gen_bool(ins) {
            seq.push(base(rng));
        }
        if !rng.gen_bool(del) {
            seq.push(if rng.gen_bool(sub) {
                substitution(b, rng)
            } else {
                b
            });
        }
    }
    let q = rng.gen_range(10_u8..=30) as i16;
    let qual = (0..seq.len())
        .map(|_| (q + rng.gen_range(-3..=3)).clamp(2, 50) as u8 + 33)
        .collect();
    (seq, qual)
}

#[derive(Clone, Copy, Debug)]
pub enum Endpoints {
    Uniform,
    /// Paired anchors per genome; each molecule jitters both ends independently.
    Hotspot {
        count: usize,
        jitter: usize,
    },
}

/// Dataset inputs. Lengths are rounded lognormal draws (median 1000, log sd
/// 0.5), clipped to the supplied bounds and genome length. Families are
/// geometric on 1,2,... with the given mean, and each read is reverse
/// complemented with probability 0.5. Genomes are independent random ACGT;
/// related strains in tests are made explicitly. Use positive sizes, ordered
/// length bounds and family_mean >= 1.
pub struct SimConfig {
    pub seed: u64,
    pub n_genomes: usize,
    pub genome_len: usize,
    pub n_molecules: usize,
    pub len_min: usize,
    pub len_max: usize,
    pub endpoints: Endpoints,
    pub family_mean: f64,
    pub errors: Errors,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Truth {
    pub molecule_id: usize,
    pub genome_id: usize,
    pub start: usize,
    pub end: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Read {
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub truth: Truth,
}

fn interval(c: &SimConfig, rng: &mut StdRng) -> (usize, usize) {
    let len = (1000_f64.ln() + 0.5 * normal(rng)).exp().round() as usize;
    let len = len.clamp(c.len_min.min(c.genome_len), c.len_max.min(c.genome_len));
    let start = rng.gen_range(0..=c.genome_len - len);
    (start, start + len)
}

fn jitter(pos: usize, amount: usize, limit: usize, rng: &mut StdRng) -> usize {
    let offset = rng.gen_range(-(amount as i128)..=amount as i128);
    (pos as i128 + offset).clamp(0, limit as i128) as usize
}

/// Generate PCR families and shuffle reads. Truth coordinates stay in the
/// forward genome orientation, regardless of each read's strand or errors.
/// Hotspot jitter is uniform integer [-J,J], then clipped to valid intervals.
pub fn simulate(c: &SimConfig) -> Vec<Read> {
    assert!(c.n_genomes > 0 && c.genome_len > 0 && c.len_min > 0);
    assert!(c.len_min <= c.len_max && c.family_mean.is_finite() && c.family_mean >= 1.0);
    let mut rng = StdRng::seed_from_u64(c.seed);
    let genomes: Vec<Vec<u8>> = (0..c.n_genomes)
        .map(|_| (0..c.genome_len).map(|_| base(&mut rng)).collect())
        .collect();
    let count = match c.endpoints {
        Endpoints::Uniform => 0,
        Endpoints::Hotspot { count, .. } => {
            assert!(count > 0);
            count
        }
    };
    let hotspots: Vec<Vec<_>> = (0..c.n_genomes)
        .map(|_| (0..count).map(|_| interval(c, &mut rng)).collect())
        .collect();
    let mut reads = Vec::new();
    for molecule_id in 0..c.n_molecules {
        let genome_id = rng.gen_range(0..c.n_genomes);
        let (start, end) = match c.endpoints {
            Endpoints::Uniform => interval(c, &mut rng),
            Endpoints::Hotspot { jitter: j, .. } => {
                let &(s, e) = hotspots[genome_id].choose(&mut rng).unwrap();
                let min = c.len_min.min(c.genome_len);
                let s = jitter(s, j, c.genome_len - min, &mut rng);
                let e = jitter(e, j, c.genome_len, &mut rng)
                    .clamp(s + min, s.saturating_add(c.len_max).min(c.genome_len));
                (s, e)
            }
        };
        let template = &genomes[genome_id][start..end];
        let mut copies = 1;
        while !rng.gen_bool(1.0 / c.family_mean) {
            copies += 1;
        }
        for _ in 0..copies {
            let seq = if rng.gen_bool(0.5) {
                reverse_complement(template)
            } else {
                template.to_vec()
            };
            let (sequence, quality) = sequence_read(&seq, c.errors, &mut rng);
            reads.push(Read {
                sequence,
                quality,
                truth: Truth {
                    molecule_id,
                    genome_id,
                    start,
                    end,
                },
            });
        }
    }
    reads.shuffle(&mut rng);
    reads
}

pub struct PairMetrics {
    pub true_pairs: u64,
    pub predicted_pairs: u64,
    pub recovered_pairs: u64,
}

impl PairMetrics {
    /// Undefined precision/recall (zero denominator) is None.
    pub fn precision(&self) -> Option<f64> {
        (self.predicted_pairs > 0)
            .then(|| self.recovered_pairs as f64 / self.predicted_pairs as f64)
    }
    pub fn recall(&self) -> Option<f64> {
        (self.true_pairs > 0).then(|| self.recovered_pairs as f64 / self.true_pairs as f64)
    }
}

pub struct Evaluation {
    pub molecule: PairMetrics,
    pub identifiability: PairMetrics,
}

fn pair_counts<L: Copy + Eq + Hash>(labels: impl Iterator<Item = L>, ids: &[usize]) -> PairMetrics {
    let (mut counts, mut cells, mut predicted) = (HashMap::new(), HashMap::new(), HashMap::new());
    for (label, &cluster) in labels.zip(ids) {
        *counts.entry(label).or_insert(0_u64) += 1;
        *cells.entry((label, cluster)).or_insert(0_u64) += 1;
        *predicted.entry(cluster).or_insert(0_u64) += 1;
    }
    let pairs = |n: &u64| n * n.saturating_sub(1) / 2;
    PairMetrics {
        true_pairs: counts.values().map(pairs).sum(),
        predicted_pairs: predicted.values().map(pairs).sum(),
        recovered_pairs: cells.values().map(pairs).sum(),
    }
}

/// Score unordered pairs in expected O(reads) time by hashing label and
/// (label, cluster) counts. Slices share input order; cluster IDs need not be
/// dense. Endpoint identifiability means exactly the same (genome,start,end).
pub fn evaluate(truth: &[Truth], clusters: &[usize]) -> Evaluation {
    assert_eq!(truth.len(), clusters.len());
    Evaluation {
        molecule: pair_counts(truth.iter().map(|t| t.molecule_id), clusters),
        identifiability: pair_counts(
            truth.iter().map(|t| (t.genome_id, t.start, t.end)),
            clusters,
        ),
    }
}
