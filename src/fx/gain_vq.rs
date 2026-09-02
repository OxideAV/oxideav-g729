//! §3.9 **gain quantisation** and the taming procedure on the fixed
//! grid.
//!
//! The candidate gains are reconstructed through the decoder's own
//! [`crate::fx::gains::GainDecoderFx::reconstruct`] (Q14 `ĝ_p`, Q1
//! `ĝ_c`), so the search scores exactly what the decoder will
//! rebuild. This stage currently evaluates the eq (63) error and the
//! §3.9.2 preselection through the spec-equation float module
//! ([`crate::gain_quantize`]); moving the six correlations and the
//! error comparison onto the Word32 grid is the §3.9 migration step
//! proper.

use crate::fx::filters::L_SUBFR;
use crate::fx::gains::{GainDecoderFx, GainPredictionFx};
use crate::gain_quantize::{quantize_gains, GainTerms};
use crate::taming::Taming;

/// The eq (63) inner products of one subframe.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GainCorrelationsFx {
    terms: GainTerms,
}

impl GainCorrelationsFx {
    /// From the Q0 target `x`, the Q0 filtered adaptive vector `y`
    /// and the Q12 filtered codevector `z`.
    #[must_use]
    pub fn compute(x: &[i16; L_SUBFR], y: &[i16; L_SUBFR], z_q12: &[i16; L_SUBFR]) -> Self {
        let xf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(x[n]));
        let yf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(y[n]));
        let zf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(z_q12[n]) / 4096.0);
        Self {
            terms: GainTerms::compute(&xf, &yf, &zf),
        }
    }
}

/// The selected codebook-domain index pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GainQuantFx {
    /// First-stage index (`0 … 7`).
    pub ga: usize,
    /// Second-stage index (`0 … 15`).
    pub gb: usize,
}

/// §3.9.2 search for one subframe. `pred` is the decoder-consistent
/// eq (71) prediction; `tame` the taming flag.
#[must_use]
pub fn quantize_gains_fx(
    corr: &GainCorrelationsFx,
    _dec: &GainDecoderFx,
    pred: &GainPredictionFx,
    tame: bool,
) -> GainQuantFx {
    // g'_c = mantissa · 2^(exp − 13).
    let g_c_prime = (f64::from(pred.g_prime_scaled) * 2f64.powi(i32::from(pred.exp) - 13)) as f32;
    let r = quantize_gains(&corr.terms, g_c_prime, tame);
    GainQuantFx { ga: r.ga, gb: r.gb }
}

/// Taming-procedure state (doc `taming-procedure.md`) driven with Q14
/// pitch gains.
#[derive(Debug, Clone, Default)]
pub struct TamingFx {
    inner: Taming,
}

impl TamingFx {
    /// Fresh state.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Doc §3 test for the candidate delay.
    #[must_use]
    pub fn test(&self, int_t: i32, frac: i32) -> bool {
        self.inner.test(int_t, frac)
    }

    /// Doc §2 update with the quantised Q14 pitch gain.
    pub fn update(&mut self, int_t: i32, frac: i32, gain_pit_q14: i16) {
        self.inner
            .update(int_t, frac, f32::from(gain_pit_q14) / 16384.0);
    }
}
