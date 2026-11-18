use nao_dedup::ont::OntDedupParams;

/// Options used by examples with varied qualities.
#[derive(clap::Args)]
pub struct QualityNoiseArgs {
    /// Multiplier on usable-quality expected errors
    #[arg(long, default_value_t = OntDedupParams::default().quality_calibration)]
    quality_calibration: f64,
    /// Additive errors per base with usable qualities
    #[arg(long, default_value_t = OntDedupParams::default().quality_error_floor)]
    quality_error_floor: f64,
    /// Multiplier on usable-quality noise
    #[arg(long, default_value_t = OntDedupParams::default().noise_scale)]
    noise_scale: f64,
}

impl QualityNoiseArgs {
    pub fn params(&self, base: OntDedupParams) -> OntDedupParams {
        OntDedupParams {
            quality_calibration: self.quality_calibration,
            quality_error_floor: self.quality_error_floor,
            noise_scale: self.noise_scale,
            ..base
        }
    }
}
