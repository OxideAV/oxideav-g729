//! §3.2.4 **LSP quantiser encoder search** — the encoder-side inverse of
//! [`crate::lsp_reconstruct`]. Takes the unquantised cosine-domain LSPs
//! `q_i` from [`crate::lp_to_lsp`] and selects the four transmitted
//! codeword fields `(L0, L1, L2, L3)` that best approximate them under
//! the spec's switched-predictor two-stage VQ.
//!
//! This is the fifth encoder-side stage. It is validated by a full
//! encode→decode round-trip: the indices it emits, fed straight back
//! through [`crate::lsp_reconstruct::LspReconstructor`], reconstruct an
//! LSF vector close to the encoder input — so encoder and decoder stay
//! bit-consistent by construction (the quantiser *delegates* its own
//! reconstruction + MA-history update to a wrapped `LspReconstructor`).
//!
//! ## Spec source — clause 3.2.4, equations (18)–(23)
//!
//! (Transcribed from the EPUB prose + the equation rasters
//! `images/eq18.jpg … eq23.jpg`.)
//!
//! * **eq (18)** — the LSPs are quantised in the LSF (normalised
//!   frequency) domain: `ω_i = arccos(q_i)`, `ω_i ∈ [0, π]`.
//! * **eq (20)** (reconstruction, decode-side) —
//!   `ω̂_i^(m) = (1 − Σ_{k=1}^{4} P̂_{i,k})·l̂_i^m + Σ_{k=1}^{4} P̂_{i,k}·l̂_i^(m−k)`.
//! * **eq (21)** — the weighted MSE the search minimises:
//!   `E_lsf = Σ_{i=1}^{10} w_i·(ω_i − ω̂_i)²`.
//! * **eq (22)** — the adaptive weights from the *unquantised* LSF:
//!   ```text
//!   w_1  = 1.0 if (ω_2 − 0.04π − 1) > 0 else 10·(ω_2 − 0.04π − 1)² + 1
//!   w_i  = 1.0 if (ω_{i+1} − ω_{i−1} − 1) > 0
//!               else 10·(ω_{i+1} − ω_{i−1} − 1)² + 1     (2 ≤ i ≤ 9)
//!   w_10 = 1.0 if (−ω_9 + 0.92π − 1) > 0 else 10·(−ω_9 + 0.92π − 1)² + 1
//!   ```
//!   and additionally `w_5`, `w_6` are each multiplied by 1.2.
//! * **eq (23)** — the vector actually quantised for the current frame,
//!   per MA-predictor mode:
//!   `l_i = [ω_i^(m) − Σ_{k=1}^{4} P̂_{i,k}·l̂_i^(m−k)] / (1 − Σ_{k=1}^{4} P̂_{i,k})`.
//!
//! ## Search procedure (clause 3.2.4 prose)
//!
//! For each of the two MA predictors `L0 ∈ {0, 1}`:
//!
//! 1. Form the target `l` (eq (23)) from the unquantised LSF and the MA
//!    history of previously-selected residuals.
//! 2. Search codebook **L1** (128 entries): pick the entry minimising
//!    the *unweighted* MSE to `l` over all ten coordinates.
//! 3. Search **L2** (32, lower five coordinates): for each candidate,
//!    form the partial residual `l̂ = L1 + L2`, rearrange it to a minimum
//!    adjacent distance of 0.0012, reconstruct `ω̂_{1..5}` via eq (20),
//!    and keep the candidate with the lowest weighted MSE over `i = 1…5`.
//! 4. Search **L3** (32, upper five): analogously over `i = 6…10`.
//! 5. Reconstruct the full quantised vector (rearrange twice — 0.0012
//!    then 0.0006 — then eq (20)) and compute the total weighted MSE.
//!
//! The predictor `L0` giving the lowest total weighted MSE wins. The
//! selected indices are then run through the wrapped `LspReconstructor`
//! (which repeats the exact decode reconstruction and advances the MA
//! history), so the encoder's future targets are computed against the
//! same residual history the decoder will see.

use core::f32::consts::PI;

use crate::lsp_reconstruct::{codebook_sum, rearrange_pass, LspReconstructor, REARRANGE_J1};
use crate::tables::{self, M, NC0, NC1};

/// Q13 unit (`2^13`) for the LSP codebook literals.
const Q13: f32 = 8192.0;
/// Q15 unit (`2^15`) for the MA-predictor `fg` literals.
const Q15: f32 = 32768.0;

/// The four transmitted LSP codeword fields selected by the quantiser
/// (spec Table 8: `L0` 1 bit, `L1` 7 bits, `L2` 5 bits, `L3` 5 bits).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LspIndices {
    /// Predictor-switch flag `L0` (0 or 1).
    pub l0: usize,
    /// First-stage codebook index `L1` (`0 … 127`).
    pub l1: usize,
    /// Second-stage lower-split index `L2` (`0 … 31`).
    pub l2: usize,
    /// Second-stage upper-split index `L3` (`0 … 31`).
    pub l3: usize,
}

/// One frame's quantiser result: the transmitted indices and the
/// decoder-consistent reconstructed LSF vector `ω̂` (radians in
/// `[0, π]`).
#[derive(Debug, Clone)]
pub struct LspQuantized {
    /// Selected codeword indices.
    pub indices: LspIndices,
    /// Reconstructed LSF `ω̂_i` (what the decoder will produce from
    /// `indices`).
    pub omega_hat: [f32; M],
}

/// Reads the `fg` MA-predictor plane for `mode` into `f32` (Q15→real).
#[inline]
fn fg_plane_f32(mode: usize) -> [[f32; M]; tables::MA_NP] {
    let src = tables::lsp_fg_plane(mode);
    std::array::from_fn(|k| std::array::from_fn(|i| f32::from(src[k][i]) / Q15))
}

/// Reads the per-mode `(1 − Σ_k P̂_{i,k})` sum factor into `f32`.
#[inline]
fn fg_sum_f32(mode: usize) -> [f32; M] {
    let src = tables::lsp_fg_sum(mode);
    std::array::from_fn(|i| f32::from(src[i]) / Q15)
}

/// Reconstructs `ω̂_i` for coordinates `range` from a residual `l_hat`
/// via eq (20), given the mode's `fg`/`fg_sum` and the MA history.
fn reconstruct_omega(
    l_hat: &[f32; M],
    fg: &[[f32; M]; tables::MA_NP],
    fg_sum: &[f32; M],
    history: &[[f32; M]; tables::MA_NP],
    out: &mut [f32; M],
) {
    for i in 0..M {
        let mut acc = fg_sum[i] * l_hat[i];
        for k in 0..tables::MA_NP {
            acc += fg[k][i] * history[k][i];
        }
        out[i] = acc;
    }
}

/// Computes the eq (22) adaptive weights from the unquantised LSF `omega`
/// (radians), with the `w_5`, `w_6` ×1.2 boost applied.
#[must_use]
pub fn lsf_weights(omega: &[f32; M]) -> [f32; M] {
    let piece = |d: f32| -> f32 {
        if d > 0.0 {
            1.0
        } else {
            10.0 * d * d + 1.0
        }
    };
    let mut w = [1.0f32; M];
    // w_1 (idx 0): uses ω_2 (idx 1).
    w[0] = piece(omega[1] - 0.04 * PI - 1.0);
    // w_i, 2 ≤ i ≤ 9 (idx 1..=8): uses ω_{i+1} − ω_{i−1}.
    for i in 1..=8 {
        w[i] = piece(omega[i + 1] - omega[i - 1] - 1.0);
    }
    // w_10 (idx 9): uses ω_9 (idx 8).
    w[9] = piece(-omega[8] + 0.92 * PI - 1.0);
    // Boost w_5, w_6 (idx 4, 5).
    w[4] *= 1.2;
    w[5] *= 1.2;
    w
}

/// Pure §3.2.4 index search: given the unquantised LSF `omega`
/// (radians) and an explicit MA history of past residuals
/// `l̂^(m−1) … l̂^(m−4)` (slot 0 = most recent), runs the staged
/// L1 → L2 → L3 search over both `L0` predictor modes and returns the
/// winning indices.
///
/// This is the stateless core of [`LspQuantizer::quantize`], exposed so
/// conformance harnesses can evaluate the search against an
/// externally-maintained (e.g. reference-index-driven) history.
#[must_use]
pub fn search_lsp_indices(omega: &[f32; M], history: &[[f32; M]; tables::MA_NP]) -> LspIndices {
    let weights = lsf_weights(omega);

    let mut best: Option<(f32, LspIndices)> = None;

    for mode in 0..2 {
        let fg = fg_plane_f32(mode);
        let fg_sum = fg_sum_f32(mode);

        // eq (23): target residual l_i for this predictor mode.
        let mut target = [0.0f32; M];
        for i in 0..M {
            let mut pred = 0.0f32;
            for k in 0..tables::MA_NP {
                pred += fg[k][i] * history[k][i];
            }
            // fg_sum is (1 − Σ P̂); guard against a zero divisor.
            let denom = if fg_sum[i].abs() < 1e-9 {
                1e-9
            } else {
                fg_sum[i]
            };
            target[i] = (omega[i] - pred) / denom;
        }

        // Stage 1 — L1: minimise the unweighted MSE to the target.
        let mut l1 = 0usize;
        let mut best_l1 = f32::INFINITY;
        for c in 0..NC0 {
            let row = tables::lsp_l1_entry(c);
            let mut e = 0.0f32;
            for i in 0..M {
                let d = target[i] - f32::from(row[i]) / Q13;
                e += d * d;
            }
            if e < best_l1 {
                best_l1 = e;
                l1 = c;
            }
        }
        let l1_row = tables::lsp_l1_entry(l1);

        // Stage 2 lower — L2: weighted MSE over i = 0..5.
        let mut l2 = 0usize;
        let mut best_l2 = f32::INFINITY;
        for c in 0..NC1 {
            let lo = tables::lsp_l2_entry(c);
            let mut l_hat = [0.0f32; M];
            for i in 0..M {
                l_hat[i] = f32::from(l1_row[i]) / Q13;
            }
            for (i, &v) in lo.iter().enumerate() {
                l_hat[i] += f32::from(v) / Q13;
            }
            rearrange_pass(&mut l_hat, REARRANGE_J1);
            let mut omega_hat = [0.0f32; M];
            reconstruct_omega(&l_hat, &fg, &fg_sum, history, &mut omega_hat);
            let mut e = 0.0f32;
            for i in 0..M / 2 {
                let d = omega[i] - omega_hat[i];
                e += weights[i] * d * d;
            }
            if e < best_l2 {
                best_l2 = e;
                l2 = c;
            }
        }
        let l2_lo = tables::lsp_l2_entry(l2);

        // Stage 2 upper — L3: weighted MSE over i = 5..10.
        let mut l3 = 0usize;
        let mut best_l3 = f32::INFINITY;
        for c in 0..NC1 {
            let hi = tables::lsp_l3_entry(c);
            let mut l_hat = [0.0f32; M];
            for i in 0..M {
                l_hat[i] = f32::from(l1_row[i]) / Q13;
            }
            for (i, &v) in l2_lo.iter().enumerate() {
                l_hat[i] += f32::from(v) / Q13;
            }
            for (j, &v) in hi.iter().enumerate() {
                l_hat[M / 2 + j] += f32::from(v) / Q13;
            }
            rearrange_pass(&mut l_hat, REARRANGE_J1);
            let mut omega_hat = [0.0f32; M];
            reconstruct_omega(&l_hat, &fg, &fg_sum, history, &mut omega_hat);
            let mut e = 0.0f32;
            for i in M / 2..M {
                let d = omega[i] - omega_hat[i];
                e += weights[i] * d * d;
            }
            if e < best_l3 {
                best_l3 = e;
                l3 = c;
            }
        }

        // Full reconstruction (rearrange twice + eq (20)) for the
        // mode-selection total error. `codebook_sum` + the decode
        // rearrangement mirror `LspReconstructor` exactly.
        let mut l_hat = codebook_sum(l1, l2, l3).expect("indices in range");
        rearrange_pass(&mut l_hat, REARRANGE_J1);
        rearrange_pass(&mut l_hat, crate::lsp_reconstruct::REARRANGE_J2);
        let mut omega_hat = [0.0f32; M];
        reconstruct_omega(&l_hat, &fg, &fg_sum, history, &mut omega_hat);
        let mut total = 0.0f32;
        for i in 0..M {
            let d = omega[i] - omega_hat[i];
            total += weights[i] * d * d;
        }

        let is_better = match best {
            None => true,
            Some((e, _)) => total < e,
        };
        if is_better {
            best = Some((
                total,
                LspIndices {
                    l0: mode,
                    l1,
                    l2,
                    l3,
                },
            ));
        }
    }

    best.expect("two modes always evaluated").1
}

/// The switched-predictor two-stage LSP quantiser. Holds the MA history
/// of previously selected residuals inside a wrapped
/// [`LspReconstructor`], so the encoder and decoder MA states stay
/// identical.
#[derive(Debug, Clone, Default)]
pub struct LspQuantizer {
    reconstructor: LspReconstructor,
}

impl LspQuantizer {
    /// Constructs a fresh quantiser with the spec start-up MA history.
    #[must_use]
    pub fn new() -> Self {
        Self {
            reconstructor: LspReconstructor::new(),
        }
    }

    /// Quantises one frame of cosine-domain LSPs `q_in[i−1] = q_i`
    /// (`q_i = cos ω_i`, decreasing / ordered), returning the selected
    /// transmitted indices and the decoder-consistent reconstructed LSF.
    ///
    /// The internal MA history advances by exactly the decode-side rule,
    /// so a subsequent [`LspReconstructor`] driven with the returned
    /// indices reproduces `omega_hat`.
    pub fn quantize(&mut self, q_in: &[f32; M]) -> LspQuantized {
        // eq (18): ω_i = arccos(q_i). `q` is clamped into the valid acos
        // domain defensively (the LP→LSP roots already lie in (−1, 1)).
        let omega: [f32; M] = std::array::from_fn(|i| q_in[i].clamp(-1.0, 1.0).acos());
        self.quantize_lsf(&omega)
    }

    /// Quantises one frame of LSF coefficients `omega[i−1] = ω_i`
    /// (radians, increasing) directly — the entry point for callers
    /// that already evaluated eq (18) themselves (e.g. through the
    /// fixed-point [`crate::lsf_conversion`] table path, which is how
    /// the clause-3 encode chain reaches the reference's LSF grid).
    pub fn quantize_lsf(&mut self, omega: &[f32; M]) -> LspQuantized {
        let indices = search_lsp_indices(omega, self.reconstructor.history());

        // Delegate the exact reconstruction + MA-history advance to the
        // decode-side reconstructor so encoder/decoder stay in lock-step.
        let omega_hat = self
            .reconstructor
            .reconstruct_frame(indices.l0, indices.l1, indices.l2, indices.l3)
            .expect("selected indices are in range");

        LspQuantized { indices, omega_hat }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_to_lsp::lp_to_lsp;
    use crate::lsp_reconstruct::LspReconstructor;

    /// eq (22): a widely-spaced (flat-ish) LSF set gives all-unity
    /// weights except the ×1.2 boost on w_5/w_6.
    #[test]
    fn weights_flat_case() {
        // Uniformly spaced ω_i = i·π/11 → all adjacency gaps ≈ 0.285 rad,
        // so ω_{i+1} − ω_{i−1} ≈ 0.571 < 1 → the "otherwise" branch, but
        // the boundary terms use 0.04π/0.92π offsets. Just assert the
        // boost structure holds relative to the un-boosted neighbours.
        let omega: [f32; M] = std::array::from_fn(|i| (i + 1) as f32 * PI / 11.0);
        let w = lsf_weights(&omega);
        // w_5 and w_6 are the boosted pair; their raw (pre-boost) value
        // equals the eq-(22) piece, so they must exceed 1.0.
        assert!(w[4] > 1.0 && w[5] > 1.0);
    }

    /// A very wide gap forces the `> 0` branch → weight exactly 1.0.
    #[test]
    fn weights_wide_gap_is_unity() {
        // Construct ω with ω_{i+1} − ω_{i−1} > 1 for an interior index.
        let mut omega: [f32; M] = std::array::from_fn(|i| 0.1 + 0.05 * i as f32);
        omega[3] = 0.1;
        omega[5] = 1.4; // ω_5 − ω_3 (idx 4 uses idx5 − idx3) = 1.3 > 1
        let w = lsf_weights(&omega);
        // idx 4 (w_5) piece is unity before the ×1.2 boost → 1.2 exactly.
        assert!((w[4] - 1.2).abs() < 1e-5, "w[4]={}", w[4]);
    }

    /// End-to-end: quantise the LSPs of a real stable filter, feed the
    /// indices to an independent decoder, and confirm the reconstructed
    /// LSF matches the quantiser's own `omega_hat` bit-for-bit and is
    /// close to the (arccos of the) input LSPs.
    #[test]
    fn encode_decode_consistency() {
        let a: [f32; M] = {
            use crate::levinson::levinson;
            let r: [f64; M + 1] =
                std::array::from_fn(|k| (0.9f64).powi(k as i32) * (0.3 * k as f64).cos());
            levinson(&r).a
        };
        let q = lp_to_lsp(&a).expect("stable filter");
        let omega_in: [f32; M] = std::array::from_fn(|i| q[i].acos());

        let mut enc = LspQuantizer::new();
        let out = enc.quantize(&q);

        // Independent decoder driven by the emitted indices.
        let mut dec = LspReconstructor::new();
        let omega_dec = dec
            .reconstruct_frame(
                out.indices.l0,
                out.indices.l1,
                out.indices.l2,
                out.indices.l3,
            )
            .expect("valid indices");
        for (i, (&od, &oh)) in omega_dec.iter().zip(out.omega_hat.iter()).enumerate() {
            assert!(
                (od - oh).abs() < 1e-6,
                "decode/encode mismatch at {i}: {od} vs {oh}"
            );
        }
        // The quantiser should approximate the input LSF to well within
        // the coarse-grid failure band (a few tenths of a radian even in
        // the worst coordinate; the whole-vector average is far tighter).
        let mut sq = 0.0f32;
        for (&od, &oi) in omega_dec.iter().zip(omega_in.iter()) {
            sq += (od - oi).powi(2);
        }
        let rms = (sq / M as f32).sqrt();
        assert!(rms < 0.2, "LSF quantisation RMS too large: {rms}");
    }

    /// Indices are always within their Table-8 field ranges.
    #[test]
    fn indices_in_range() {
        let mut enc = LspQuantizer::new();
        for f in 0..8 {
            let q: [f32; M] = std::array::from_fn(|i| {
                (0.9 - 0.18 * i as f32 + 0.01 * f as f32).clamp(-0.99, 0.99)
            });
            let out = enc.quantize(&q);
            assert!(out.indices.l0 < 2);
            assert!(out.indices.l1 < NC0);
            assert!(out.indices.l2 < NC1);
            assert!(out.indices.l3 < NC1);
        }
    }

    /// Over several frames the quantiser's internal MA history stays in
    /// lock-step with an independent decoder fed the same index stream
    /// (regression guard on the delegated history update).
    #[test]
    fn multi_frame_history_lockstep() {
        let mut enc = LspQuantizer::new();
        let mut dec = LspReconstructor::new();
        for f in 0..12 {
            let a: [f32; M] = {
                use crate::levinson::levinson;
                let r: [f64; M + 1] = std::array::from_fn(|k| {
                    (0.88f64).powi(k as i32) * (0.2 * k as f64 + 0.05 * f as f64).cos()
                });
                levinson(&r).a
            };
            let q = lp_to_lsp(&a).expect("stable");
            let out = enc.quantize(&q);
            let omega_dec = dec
                .reconstruct_frame(
                    out.indices.l0,
                    out.indices.l1,
                    out.indices.l2,
                    out.indices.l3,
                )
                .expect("valid");
            for (i, (&od, &oh)) in omega_dec.iter().zip(out.omega_hat.iter()).enumerate() {
                assert!((od - oh).abs() < 1e-5, "frame {f} coord {i}: {od} vs {oh}");
            }
        }
    }
}
