use nao_dedup::ont::OntDedupParams;

#[derive(clap::Args)]
pub struct DedupArgs {
    #[arg(long, default_value_t = OntDedupParams::default().max_end_offset)]
    pub max_end_offset: usize,
    #[arg(long, default_value_t = OntDedupParams::default().end_check_len)]
    pub end_check_len: usize,
    #[arg(long, default_value_t = OntDedupParams::default().max_end_errors)]
    pub max_end_errors: usize,
    /// Upper cap on whole-read edit distance as a fraction of the longer read
    #[arg(long, default_value_t = OntDedupParams::default().max_error_frac)]
    pub max_error_frac: f64,
    /// Strain divergence against which to protect
    #[arg(long, default_value_t = OntDedupParams::default().protect_divergence)]
    pub protect_divergence: f64,
    /// Provisional extra relative uncertainty in pair distance
    #[arg(long, default_value_t = OntDedupParams::default().pair_distance_cv)]
    pub pair_distance_cv: f64,
    /// Independent fractional indel-drift allowance for seed offsets
    #[arg(long, default_value_t = OntDedupParams::default().seed_drift_frac)]
    pub seed_drift_frac: f64,
    /// Expected error rate for empty or exactly constant qualities
    #[arg(long, default_value_t = OntDedupParams::default().assumed_error_rate)]
    pub assumed_error_rate: f64,

    /// Only match reads in the same orientation.
    #[arg(long)]
    pub no_rc: bool,
}

impl DedupArgs {
    pub fn params(&self) -> OntDedupParams {
        OntDedupParams {
            max_end_offset: self.max_end_offset,
            end_check_len: self.end_check_len,
            max_end_errors: self.max_end_errors,
            max_error_frac: self.max_error_frac,
            protect_divergence: self.protect_divergence,
            pair_distance_cv: self.pair_distance_cv,
            seed_drift_frac: self.seed_drift_frac,
            assumed_error_rate: self.assumed_error_rate,
            rc_tolerant: !self.no_rc,
            ..Default::default()
        }
    }
}
