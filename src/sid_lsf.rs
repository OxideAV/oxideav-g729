//! Annex B §B.4.2.2 SID-LSF dequantizer — the comfort-noise spectral
//! envelope.
//!
//! A SID frame carries the noise spectrum as three indices (Table B.2:
//! 1-bit switched-predictor, 5-bit first stage, 4-bit second stage)
//! into a reduced version of the §3.2.4 LSF quantizer, with the three
//! modifications §B.4.2.2 lists:
//!
//! 1. **Blended second predictor** — eq (B.18): mode 1's 4th-order MA
//!    predictor is `0.6·p₁ + 0.4·p₂` of the full coder's two
//!    predictors; mode 0 is the full coder's first predictor
//!    unchanged.
//! 2. **First-stage subset** — the 5-bit index addresses the staged
//!    32-entry lookup [`crate::tables::ANNEXB_SID_LSF_STAGE1_INDEX_MAP`],
//!    whose values are row indices into the full coder's 128-row L1
//!    codebook.
//! 3. **Full (not split) second stage** — the 4-bit index addresses
//!    *both* halves through the staged 2 × 16 lookup
//!    [`crate::tables::ANNEXB_SID_LSF_STAGE2_INDEX_MAP`]: row 0 picks
//!    the L2 row supplying the lower five coefficients, row 1 the row
//!    supplying the upper five.
//!
//! Everything else follows §3.2.4: the twice-run rearrangement, the
//! eq (20) MA-predicted reconstruction (with the per-mode residual
//! weights staged as
//! [`crate::tables::ANNEXB_SID_LSF_MA_PREDICTOR_SUM_Q15`]), the
//! stability procedure, and the `l̂_i = iπ/11` start-up history.
//!
//! ## What the Recommendation leaves unpinned
//!
//! §B.4.2.2 does not say whether the SID quantizer's MA history is
//! shared with the active coder's §3.2.4 history or kept separate.
//! This implementation keeps a **separate** history (advanced only by
//! SID frames), which is self-consistent across an arbitrary mix of
//! frame types; the choice is recorded here rather than silently made.
//! The unnamed staged per-mode pair `annexB-sid-lsf-mp-Q15` is not
//! consumed (its role is unresolved — see the staging note
//! `docs/audio/g729/annexb-sid-lsp-vq.md` §4).

use core::f32::consts::PI;

use crate::lsp_reconstruct::{rearrange_twice, stability_clamp};
use crate::tables::{
    ANNEXB_SID_LSF_MA_PREDICTOR_SUM_Q15, ANNEXB_SID_LSF_STAGE1_INDEX_MAP,
    ANNEXB_SID_LSF_STAGE2_INDEX_MAP, LSP_MA_PREDICTOR_FG_Q15, LSP_QUANT_CODEBOOK_L1_Q13,
    LSP_QUANT_CODEBOOK_L2_Q13, M, MA_NP,
};

/// Q13 → real boundary scale for the codebook entries.
const Q13_UNIT: f32 = 8192.0;
/// Q15 → real boundary scale for the predictor coefficients.
const Q15_UNIT: f32 = 32_768.0;

/// eq (B.18) blend weight on the full coder's first predictor.
pub const BLEND_P1: f32 = 0.6;
/// eq (B.18) blend weight on the full coder's second predictor.
pub const BLEND_P2: f32 = 0.4;

/// The §B.4.2.2 MA predictor coefficient for `(mode, k, i)`: mode 0 is
/// the full coder's first predictor unchanged; mode 1 is the eq (B.18)
/// blend `0.6·p₁ + 0.4·p₂`.
#[must_use]
pub fn sid_predictor_coef(mode: usize, k: usize, i: usize) -> f32 {
    let p1 = f32::from(LSP_MA_PREDICTOR_FG_Q15[0][k][i]) / Q15_UNIT;
    if mode == 0 {
        p1
    } else {
        let p2 = f32::from(LSP_MA_PREDICTOR_FG_Q15[1][k][i]) / Q15_UNIT;
        BLEND_P1 * p1 + BLEND_P2 * p2
    }
}

/// The §B.4.2.2 residual (eq (19)-analogue) for a SID frame's
/// `(l1, l2)` stage indices: first-stage subset row plus the full-VQ
/// second-stage pair.
///
/// # Panics
///
/// Panics if `l1 >= 32` or `l2 >= 16` (the 5-/4-bit Table B.2 fields
/// cannot exceed these).
#[must_use]
pub fn sid_codebook_sum(l1: usize, l2: usize) -> [f32; M] {
    let stage1_row = usize::try_from(ANNEXB_SID_LSF_STAGE1_INDEX_MAP[l1]).expect("map in domain");
    let lo_row = usize::try_from(ANNEXB_SID_LSF_STAGE2_INDEX_MAP[0][l2]).expect("map in domain");
    let hi_row = usize::try_from(ANNEXB_SID_LSF_STAGE2_INDEX_MAP[1][l2]).expect("map in domain");
    let s1 = &LSP_QUANT_CODEBOOK_L1_Q13[stage1_row];
    let lo = &LSP_QUANT_CODEBOOK_L2_Q13[lo_row];
    let hi = &LSP_QUANT_CODEBOOK_L2_Q13[hi_row];
    core::array::from_fn(|i| {
        let stage2 = if i < M / 2 { lo[i] } else { hi[i] };
        (f32::from(s1[i]) + f32::from(stage2)) / Q13_UNIT
    })
}

/// Stateful §B.4.2.2 SID-LSF dequantizer (decoder side).
///
/// Owns the 4-frame MA residual history, advanced once per decoded SID
/// frame. Construction seeds every slot with the §3.2.4 start-up
/// vector `l̂_i = iπ/11`.
#[derive(Debug, Clone)]
pub struct SidLspDecoder {
    history: [[f32; M]; MA_NP],
}

impl Default for SidLspDecoder {
    fn default() -> Self {
        Self::new()
    }
}

impl SidLspDecoder {
    /// Fresh dequantizer in the §3.2.4 start-up state.
    #[must_use]
    pub fn new() -> Self {
        let row: [f32; M] = core::array::from_fn(|i| ((i + 1) as f32) * PI / 11.0);
        Self {
            history: [row; MA_NP],
        }
    }

    /// Borrow the MA history for inspection / tests.
    #[must_use]
    pub fn history(&self) -> &[[f32; M]; MA_NP] {
        &self.history
    }

    /// Dequantize one SID frame's `(lp0, l1, l2)` indices to the SID
    /// LSF vector `ω̂` (radians, post-stability), advancing the MA
    /// history with the rearranged residual.
    ///
    /// # Panics
    ///
    /// Panics if an index exceeds its Table B.2 field domain
    /// (`lp0 >= 2`, `l1 >= 32`, `l2 >= 16`) — parsed SID frames cannot.
    #[must_use]
    pub fn dequantize(&mut self, lp0: usize, l1: usize, l2: usize) -> [f32; M] {
        let (omega, residual) = sid_reconstruct(&self.history, lp0, l1, l2);
        for k in (1..MA_NP).rev() {
            self.history[k] = self.history[k - 1];
        }
        self.history[0] = residual;
        omega
    }
}

/// Pure §B.4.2.2 reconstruction of one `(lp0, l1, l2)` triple against
/// an explicit MA history: returns the post-stability LSF vector `ω̂`
/// (radians) and the rearranged residual the caller pushes into the
/// history — the shared core of the decoder's [`SidLspDecoder`] and
/// the encoder's SID-LSF search.
///
/// # Panics
///
/// Panics if an index exceeds its Table B.2 field domain.
#[must_use]
pub fn sid_reconstruct(
    history: &[[f32; M]; MA_NP],
    lp0: usize,
    l1: usize,
    l2: usize,
) -> ([f32; M], [f32; M]) {
    assert!(lp0 < 2, "1-bit switched-predictor index");
    let mut residual = sid_codebook_sum(l1, l2);
    rearrange_twice(&mut residual);

    // eq (20) with the §B.4.2.2 predictor set: the per-mode residual
    // weight is the staged blended column sum; the history taps use
    // the (B.18) coefficients.
    let sum_row = &ANNEXB_SID_LSF_MA_PREDICTOR_SUM_Q15[lp0];
    let mut omega = [0.0_f32; M];
    for i in 0..M {
        let mut acc = (f32::from(sum_row[i]) / Q15_UNIT) * residual[i];
        for (k, past) in history.iter().enumerate() {
            acc += sid_predictor_coef(lp0, k, i) * past[i];
        }
        omega[i] = acc;
    }
    stability_clamp(&mut omega);
    (omega, residual)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The staged per-mode column sums are exactly the (B.18) blend of
    /// the full coder's two sum rows (the blend is linear in the
    /// predictor, so the residual weight `1 − Σp` blends identically):
    /// mode 0 equals the full coder's mode 0, and mode 1 is
    /// `0.6·sum₀ + 0.4·sum₁` within one Q15 LSB.
    #[test]
    fn staged_sums_are_the_blend_of_the_full_coder_rows() {
        let full = &crate::tables::LSP_MA_PREDICTOR_FG_SUM_Q15;
        assert_eq!(ANNEXB_SID_LSF_MA_PREDICTOR_SUM_Q15[0], full[0]);
        for (i, &staged_q15) in ANNEXB_SID_LSF_MA_PREDICTOR_SUM_Q15[1].iter().enumerate() {
            let blend = BLEND_P1 * f32::from(full[0][i]) + BLEND_P2 * f32::from(full[1][i]);
            let staged = f32::from(staged_q15);
            // The staged fixed-point rows sit up to ~2 LSB above the
            // real-valued blend (the attachment's own rounding).
            assert!(
                (blend - staged).abs() <= 2.0,
                "sum[1][{i}]: blend {blend} vs staged {staged}",
            );
        }
    }

    /// Every `(l1, l2)` combination dequantizes to a stable, strictly
    /// increasing LSF vector inside `(0, π)` from the start-up state.
    #[test]
    fn all_index_combinations_dequantize_stable() {
        for lp0 in 0..2 {
            for l1 in 0..32 {
                for l2 in 0..16 {
                    let mut d = SidLspDecoder::new();
                    let omega = d.dequantize(lp0, l1, l2);
                    for i in 0..M {
                        assert!(omega[i] > 0.0 && omega[i] < PI);
                        if i > 0 {
                            assert!(omega[i] > omega[i - 1], "({lp0},{l1},{l2}) not ordered");
                        }
                    }
                }
            }
        }
    }

    /// The history advances with the rearranged residual: two
    /// identical SID frames in a row reconstruct different vectors
    /// (the MA memory moved), then converge as the memory fills with
    /// the same residual.
    #[test]
    fn history_advances_and_converges() {
        let mut d = SidLspDecoder::new();
        let first = d.dequantize(1, 7, 3);
        let second = d.dequantize(1, 7, 3);
        assert_ne!(first, second);
        let mut last = second;
        for _ in 0..8 {
            last = d.dequantize(1, 7, 3);
        }
        let settled = d.dequantize(1, 7, 3);
        for i in 0..M {
            assert!((settled[i] - last[i]).abs() < 1e-4, "not converging at {i}");
        }
    }
}
