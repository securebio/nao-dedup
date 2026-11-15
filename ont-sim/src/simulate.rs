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

/// Runtime model parameters. Sequencing error probabilities are deliberately required.
#[derive(Clone, Debug, Args)]
pub struct SimConfig {
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
}

impl SimConfig {
    pub fn validate(&self) -> Result<(), String> {
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
    let z = (-2.0 * rng.sample::<f64, _>(Open01).ln()).sqrt()
        * (std::f64::consts::TAU * rng.gen_range(0.0..1.0)).cos();
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
            let (sequence, [n_sub, n_ins, n_del]) = add_errors(
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
    Ok(Simulation { genomes, reads })
}
