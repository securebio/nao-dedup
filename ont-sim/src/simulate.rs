use crate::Truth;
use clap::{Args, ValueEnum};
use rand::{Rng, SeedableRng, distributions::Open01, rngs::StdRng, seq::SliceRandom};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Endpoints {
    Uniform,
    Hotspot,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum FamilyDist {
    Geometric,
}

/// Opt-in quality-conditioned model; all distribution inputs are runtime supplied.
#[derive(Clone, Debug, Default, Args)]
pub struct QualityConfig {
    #[arg(long)]
    pub quality_conditioned: bool,
    /// Median nominal per-read error probability (lognormal, clipped to 1).
    #[arg(long, required_if_eq("quality_conditioned", "true"))]
    pub error_median: Option<f64>,
    /// Standard deviation of log nominal error probability.
    #[arg(long, required_if_eq("quality_conditioned", "true"))]
    pub error_sigma: Option<f64>,
    /// Local standard deviation in Phred units, before rounding to Q0–Q93.
    #[arg(long, required_if_eq("quality_conditioned", "true"))]
    pub local_q_sigma: Option<f64>,
    /// Batch multiplier on quality-derived error probabilities.
    #[arg(long, required_if_eq("quality_conditioned", "true"))]
    pub miscalibration: Option<f64>,
    /// Log standard deviation of the mean-one per-read latent multiplier.
    #[arg(long, required_if_eq("quality_conditioned", "true"))]
    pub latent_sigma: Option<f64>,
    /// Additive event probability per source base, without changing reported Q.
    #[arg(long, default_value_t = 0.0)]
    pub unreported_error_floor: f64,
}

impl QualityConfig {
    /// Construct an enabled quality model with no unreported error floor.
    pub fn enabled(
        error_median: f64,
        error_sigma: f64,
        local_q_sigma: f64,
        miscalibration: f64,
        latent_sigma: f64,
    ) -> Self {
        Self {
            quality_conditioned: true,
            error_median: Some(error_median),
            error_sigma: Some(error_sigma),
            local_q_sigma: Some(local_q_sigma),
            miscalibration: Some(miscalibration),
            latent_sigma: Some(latent_sigma),
            unreported_error_floor: 0.0,
        }
    }

    pub fn validate(&self, proportions: [f64; 3], homopolymer_mult: f64) -> Result<(), String> {
        if !self.quality_conditioned {
            return Ok(());
        }
        let positive = |v: Option<f64>| v.is_some_and(|x| x.is_finite() && x > 0.0);
        let nonnegative = |v: Option<f64>| v.is_some_and(|x| x.is_finite() && x >= 0.0);
        if !positive(self.error_median)
            || self.error_median.is_some_and(|x| x > 1.0)
            || !nonnegative(self.error_sigma)
            || !positive(self.local_q_sigma)
            || !nonnegative(self.miscalibration)
            || !nonnegative(self.latent_sigma)
            || !self.unreported_error_floor.is_finite()
            || self.unreported_error_floor < 0.0
            || proportions.iter().any(|x| !x.is_finite() || *x < 0.0)
            || !proportions.iter().sum::<f64>().is_finite()
            || proportions.iter().sum::<f64>() <= 0.0
            || !homopolymer_mult.is_finite()
            || homopolymer_mult < 0.0
        {
            return Err("quality mode requires 0 < error-median <= 1, finite error-sigma/latent-sigma/miscalibration/unreported-error-floor >= 0, local-q-sigma > 0, and finite nonnegative event weights with positive sum and homopolymer multiplier".into());
        }
        Ok(())
    }
}

fn normal(rng: &mut StdRng) -> f64 {
    (-2.0 * rng.sample::<f64, _>(Open01).ln()).sqrt()
        * (std::f64::consts::TAU * rng.gen_range(0.0..1.0)).cos()
}

/// Sequencing event counts and source positions whose probability exceeded one.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QualityErrorStats {
    pub events: [usize; 3],
    pub saturated_positions: usize,
}

/// Generate qualities first, then categorical sub/ins/del events from their
/// probabilities times the batch multiplier, plus the per-base floor, all
/// multiplied by the mean-one lognormal read multiplier and capped at one.
/// Event weights are normalized. Insertions precede the retained source base and
/// share its pre-sampled quality; deletions emit neither base nor quality.
/// Homopolymer indel enrichment is reflected in the nominal qualities too.
/// Panics unless quality mode is enabled and `config.validate(weights, mult)` passes.
pub fn add_quality_errors(
    s: &[u8],
    config: &QualityConfig,
    weights: [f64; 3],
    mult: f64,
    rng: &mut StdRng,
) -> (Vec<u8>, Vec<u8>, QualityErrorStats) {
    config
        .validate(weights, mult)
        .expect("invalid quality model");
    assert!(config.quality_conditioned, "quality mode must be enabled");
    let nominal_log =
        (config.error_median.unwrap().ln() + config.error_sigma.unwrap() * normal(rng)).min(0.0);
    let sigma = config.latent_sigma.unwrap();
    let latent_log = sigma * (normal(rng) - 0.5 * sigma);
    let batch_log = config.miscalibration.unwrap().ln();
    let total = weights.iter().sum::<f64>();
    let weights = weights.map(|w| w / total);
    let (mut out, mut qual, mut counts) = (
        Vec::with_capacity(s.len()),
        Vec::with_capacity(s.len()),
        [0; 3],
    );
    let mut saturated_positions = 0;
    let mut start = 0;
    while start < s.len() {
        let mut end = start + 1;
        while end < s.len() && s[end] == s[start] {
            end += 1;
        }
        let hp = if end - start >= 3 { mult } else { 1.0 };
        let local_weights = [weights[0], weights[1] * hp, weights[2] * hp];
        let local_total = local_weights.iter().sum::<f64>();
        for &b in &s[start..end] {
            let q = (-10.0 / std::f64::consts::LN_10 * (nominal_log + local_total.ln())
                + config.local_q_sigma.unwrap() * normal(rng))
            .round()
            .clamp(0.0, 93.0) as u8;
            let quality_probability =
                (-f64::from(q) * std::f64::consts::LN_10 / 10.0 + batch_log + latent_log).exp();
            let unbounded =
                quality_probability + (config.unreported_error_floor.ln() + latent_log).exp();
            saturated_positions += usize::from(unbounded > 1.0);
            let probability = unbounded.min(1.0);
            if local_total == 0.0 || !rng.gen_bool(probability) {
                out.push(b);
                qual.push(q + 33);
                continue;
            }
            let event = rng.gen_range(0.0..local_total);
            if event < local_weights[0] {
                out.push(substitution(b, rng));
                qual.push(q + 33);
                counts[0] += 1;
            } else if event < local_weights[0] + local_weights[1] {
                out.push(base(rng));
                out.push(b);
                qual.extend([q + 33, q + 33]);
                counts[1] += 1;
            } else {
                counts[2] += 1;
            }
        }
        start = end;
    }
    (
        out,
        qual,
        QualityErrorStats {
            events: counts,
            saturated_positions,
        },
    )
}

/// Runtime model parameters. Sequencing error probabilities are deliberately required.
#[derive(Clone, Debug, Args)]
pub struct SimConfig {
    #[command(flatten)]
    pub quality: QualityConfig,
    #[arg(long, default_value_t = 1)]
    pub seed: u64,
    #[arg(long, default_value_t = 1)]
    pub n_genomes: usize,
    #[arg(long, default_value_t = 100_000)]
    pub genome_len: usize,
    /// Per-base substitution probability from genome 0; zero draws independent genomes.
    #[arg(long, default_value_t = 0.0)]
    pub genome_divergence: f64,
    #[arg(long, default_value_t = 1_000)]
    pub n_molecules: usize,
    #[arg(long, default_value_t = 1_000.0)]
    pub len_median: f64,
    #[arg(long, default_value_t = 0.5)]
    pub len_sigma: f64,
    #[arg(long, default_value_t = 100)]
    pub len_min: usize,
    #[arg(long, default_value_t = 15_000)]
    pub len_max: usize,
    #[arg(long, value_enum, default_value = "uniform")]
    pub endpoints: Endpoints,
    #[arg(long, default_value_t = 10)]
    pub n_hotspots: usize,
    #[arg(long, default_value_t = 3)]
    pub hotspot_jitter: usize,
    #[arg(long, value_enum, default_value = "geometric")]
    pub family_dist: FamilyDist,
    #[arg(long, default_value_t = 3.0)]
    pub family_mean: f64,
    #[arg(long, default_value_t = 0.0)]
    pub pcr_error_prob: f64,
    #[arg(long, default_value_t = 0.5)]
    pub rc_prob: f64,
    #[arg(long)]
    pub sub_rate: f64,
    #[arg(long)]
    pub ins_rate: f64,
    #[arg(long)]
    pub del_rate: f64,
    #[arg(long, default_value_t = 1.0)]
    pub homopolymer_indel_mult: f64,
    #[arg(long, default_value_t = 0.0)]
    pub foldback_prob: f64,
    #[arg(long, default_value_t = 0.1)]
    pub foldback_min_frac: f64,
    #[arg(long, default_value_t = 0.5)]
    pub foldback_max_frac: f64,
    #[arg(long, default_value_t = 0.0)]
    pub polya_prob: f64,
    #[arg(long, default_value_t = 10)]
    pub polya_min: usize,
    #[arg(long, default_value_t = 100)]
    pub polya_max: usize,
    #[arg(long, default_value_t = 0)]
    pub end_jitter: usize,
    #[arg(long, default_value_t = 10)]
    pub qual_min: u8,
    #[arg(long, default_value_t = 30)]
    pub qual_max: u8,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Read {
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub truth: Truth,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Simulation {
    pub genomes: Vec<Vec<u8>>,
    pub reads: Vec<Read>,
    /// Source positions with a quality-conditioned event probability above one.
    pub saturated_positions: usize,
}

impl SimConfig {
    pub fn validate(&self) -> Result<(), String> {
        self.quality.validate(
            [self.sub_rate, self.ins_rate, self.del_rate],
            self.homopolymer_indel_mult,
        )?;
        let probabilities = [
            self.genome_divergence,
            self.sub_rate,
            self.ins_rate,
            self.del_rate,
            self.rc_prob,
            self.pcr_error_prob,
            self.foldback_prob,
            self.polya_prob,
            self.foldback_min_frac,
            self.foldback_max_frac,
        ];
        if probabilities.iter().any(|p| !(0.0..=1.0).contains(p)) {
            return Err("probabilities and foldback fractions must be finite and in [0,1]".into());
        }
        if self.n_genomes == 0
            || self.genome_len == 0
            || self.len_min == 0
            || self.len_min > self.len_max
            || (matches!(self.endpoints, Endpoints::Hotspot) && self.n_hotspots == 0)
            || !self.len_median.is_finite()
            || self.len_median <= 0.0
            || !self.len_sigma.is_finite()
            || self.len_sigma < 0.0
            || !self.family_mean.is_finite()
            || self.family_mean < 1.0
            || !self.homopolymer_indel_mult.is_finite()
            || self.homopolymer_indel_mult < 0.0
            || self.foldback_min_frac > self.foldback_max_frac
            || self.polya_min > self.polya_max
            || self.qual_min < 2
            || self.qual_max > 50
            || self.qual_min > self.qual_max
        {
            return Err("invalid sizes, ranges, or distribution parameters".into());
        }
        Ok(())
    }
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
/// Reverse complement of an uppercase ACGT sequence.
pub fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect()
}

fn interval(c: &SimConfig, rng: &mut StdRng) -> (usize, usize) {
    let z = normal(rng);
    let len = (c.len_median.ln() + c.len_sigma * z).exp().round() as usize;
    let len = len.clamp(c.len_min.min(c.genome_len), c.len_max.min(c.genome_len));
    let start = rng.gen_range(0..=c.genome_len - len);
    (start, start + len)
}

fn jitter(pos: usize, amount: usize, limit: usize, rng: &mut StdRng) -> usize {
    let offset = rng.gen_range(-(amount as i128)..=amount as i128);
    (pos as i128 + offset).clamp(0, limit as i128) as usize
}

fn jitter_ends(s: Vec<u8>, amount: usize, rng: &mut StdRng) -> Vec<u8> {
    if amount == 0 {
        return s;
    }
    let left = rng.gen_range(0..=amount);
    let right = rng.gen_range(0..=amount);
    let trim_left = rng.gen_bool(0.5);
    let trim_right = rng.gen_bool(0.5);
    let start = if trim_left { left.min(s.len()) } else { 0 };
    let end = if trim_right {
        s.len().saturating_sub(right).max(start)
    } else {
        s.len()
    };
    let mut out = Vec::new();
    if !trim_left {
        out.extend((0..left).map(|_| base(rng)));
    }
    out.extend_from_slice(&s[start..end]);
    if !trim_right {
        out.extend((0..right).map(|_| base(rng)));
    }
    out
}

/// Add independent per-source-base ONT errors to an uppercase ACGT sequence.
/// Returns the emitted sequence and [substitution, insertion, deletion] counts.
/// Insertions precede their source base; substitutions are conditional on retention.
/// Runs of at least three identical bases multiply indel probabilities (capped at 1).
/// Callers must supply rates in [0,1] and a finite, nonnegative multiplier.
pub fn add_errors(
    s: &[u8],
    sub_rate: f64,
    ins_rate: f64,
    del_rate: f64,
    homopolymer_indel_mult: f64,
    rng: &mut StdRng,
) -> (Vec<u8>, [usize; 3]) {
    let (mut out, mut counts, mut start) = (Vec::with_capacity(s.len()), [0; 3], 0);
    while start < s.len() {
        let mut end = start + 1;
        while end < s.len() && s[end] == s[start] {
            end += 1;
        }
        let mult = if end - start >= 3 {
            homopolymer_indel_mult
        } else {
            1.0
        };
        for &b in &s[start..end] {
            if rng.gen_bool((ins_rate * mult).min(1.0)) {
                out.push(base(rng));
                counts[1] += 1;
            }
            if rng.gen_bool((del_rate * mult).min(1.0)) {
                counts[2] += 1;
            } else if rng.gen_bool(sub_rate) {
                out.push(substitution(b, rng));
                counts[0] += 1;
            } else {
                out.push(b);
            }
        }
        start = end;
    }
    (out, counts)
}

/// Generate and shuffle an in-memory dataset. Indices and names refer to output order.
pub fn simulate(c: &SimConfig) -> Result<Simulation, String> {
    c.validate()?;
    let mut rng = StdRng::seed_from_u64(c.seed);
    // Keep the zero-divergence path's draws exactly as before.
    let mut genomes: Vec<Vec<u8>> = Vec::with_capacity(c.n_genomes);
    for i in 0..c.n_genomes {
        let genome = if i == 0 || c.genome_divergence == 0.0 {
            (0..c.genome_len).map(|_| base(&mut rng)).collect()
        } else {
            genomes[0]
                .iter()
                .map(|&b| {
                    if rng.gen_bool(c.genome_divergence) {
                        substitution(b, &mut rng)
                    } else {
                        b
                    }
                })
                .collect()
        };
        genomes.push(genome);
    }
    let hotspots: Vec<Vec<_>> = (0..c.n_genomes)
        .map(|_| {
            (0..if matches!(c.endpoints, Endpoints::Hotspot) {
                c.n_hotspots
            } else {
                0
            })
                .map(|_| interval(c, &mut rng))
                .collect()
        })
        .collect();
    let mut reads = Vec::new();
    let mut saturated_positions = 0;
    for molecule_id in 0..c.n_molecules {
        let genome_id = rng.gen_range(0..c.n_genomes);
        let (start, end) = match c.endpoints {
            Endpoints::Uniform => interval(c, &mut rng),
            Endpoints::Hotspot => {
                let &(s, e) = hotspots[genome_id].choose(&mut rng).unwrap();
                let min = c.len_min.min(c.genome_len);
                let s = jitter(s, c.hotspot_jitter, c.genome_len - min, &mut rng);
                let e = jitter(e, c.hotspot_jitter, c.genome_len, &mut rng)
                    .clamp(s + min, s.saturating_add(c.len_max).min(c.genome_len));
                (s, e)
            }
        };
        let template = &genomes[genome_id][start..end];
        let mut k = 1;
        while !rng.gen_bool(1.0 / c.family_mean) {
            k += 1;
        }
        let mut lineage = vec![false; k];
        let mutation = if rng.gen_bool(c.pcr_error_prob) {
            lineage[..k / 2].fill(true);
            lineage.shuffle(&mut rng);
            let pos = rng.gen_range(0..template.len());
            Some((pos, substitution(template[pos], &mut rng)))
        } else {
            None
        };
        for mutated in lineage {
            let mut s = template.to_vec();
            if let Some((pos, b)) = mutation.filter(|_| mutated) {
                s[pos] = b;
            }
            let is_rc = rng.gen_bool(c.rc_prob);
            if is_rc {
                s = reverse_complement(&s);
            }
            let foldback_len = if rng.gen_bool(c.foldback_prob) {
                (s.len() as f64 * rng.gen_range(c.foldback_min_frac..=c.foldback_max_frac)).floor()
                    as usize
            } else {
                0
            };
            s.extend(reverse_complement(&s[s.len() - foldback_len..]));
            let polya_len = if rng.gen_bool(c.polya_prob) {
                rng.gen_range(c.polya_min..=c.polya_max)
            } else {
                0
            };
            s.extend(std::iter::repeat_n(b'A', polya_len));
            s = jitter_ends(s, c.end_jitter, &mut rng);
            let (sequence, quality, [n_sub, n_ins, n_del]) = if c.quality.quality_conditioned {
                let (sequence, quality, stats) = add_quality_errors(
                    &s,
                    &c.quality,
                    [c.sub_rate, c.ins_rate, c.del_rate],
                    c.homopolymer_indel_mult,
                    &mut rng,
                );
                saturated_positions += stats.saturated_positions;
                (sequence, quality, stats.events)
            } else {
                let (sequence, counts) = add_errors(
                    &s,
                    c.sub_rate,
                    c.ins_rate,
                    c.del_rate,
                    c.homopolymer_indel_mult,
                    &mut rng,
                );
                let q = rng.gen_range(c.qual_min..=c.qual_max) as i16;
                let quality: Vec<u8> = (0..sequence.len())
                    .map(|_| (q + rng.gen_range(-3..=3)).clamp(2, 50) as u8 + 33)
                    .collect();
                (sequence, quality, counts)
            };
            let mean_q = if quality.is_empty() {
                0.0
            } else {
                quality.iter().map(|b| (b - 33) as f64).sum::<f64>() / quality.len() as f64
            };
            reads.push(Read {
                sequence,
                quality,
                truth: Truth {
                    read_idx: 0,
                    read_name: String::new(),
                    molecule_id,
                    genome_id,
                    mol_start: start,
                    mol_end: end,
                    is_rc,
                    foldback_len,
                    polya_len,
                    n_sub,
                    n_ins,
                    n_del,
                    mean_q,
                },
            });
        }
    }
    reads.shuffle(&mut rng);
    for (i, r) in reads.iter_mut().enumerate() {
        r.truth.read_idx = i;
        r.truth.read_name = format!("r{i}");
    }
    Ok(Simulation {
        genomes,
        reads,
        saturated_positions,
    })
}
