use nao_dedup::ont::OntDedupParams;

#[derive(clap::Args)]
pub struct DedupArgs {
    #[arg(long, default_value_t = OntDedupParams::default().max_end_offset)]
    pub max_end_offset: usize,
    /// Strain divergence against which to protect
    #[arg(long, default_value_t = OntDedupParams::default().protect_divergence)]
    pub protect_divergence: f64,
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
            protect_divergence: self.protect_divergence,
            assumed_error_rate: self.assumed_error_rate,
            rc_tolerant: !self.no_rc,
        }
    }
}
