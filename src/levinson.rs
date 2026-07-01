//! §3.2.2 **Levinson-Durbin recursion** — converts the modified
//! autocorrelations `r'(0 … M)` from [`crate::lp_analysis`] into the
//! 10th-order LP coefficients `a_i` (`i = 1 … M`).
//!
//! This is the third encoder-side stage. Its output feeds §3.2.3
//! LP→LSP conversion (a follow-up module); its second-order by-product
//! (the reflection coefficients) feeds the §3.3 adaptive perceptual
//! weighting.
//!
//! ## Spec source — clause 3.2.2, equations (8), (9)
//!
//! Per clause 3.2.2 (transcribed from the EPUB prose) the LP filter
//! coefficients are solved from the modified autocorrelations `r'(k)`
//! using the standard Levinson-Durbin recursion. Denoting the order-`i`
//! solution `a_j^{(i)}` and the prediction-error energy `E^{(i)}`:
//!
//! ```text
//! E^{(0)} = r'(0)
//! for i = 1 … M:
//!     k_i = −( r'(i) + Σ_{j=1}^{i−1} a_j^{(i−1)}·r'(i−j) ) / E^{(i−1)}   (8)
//!     a_i^{(i)} = k_i
//!     a_j^{(i)} = a_j^{(i−1)} + k_i·a_{i−j}^{(i−1)}   j = 1 … i−1        (9)
//!     E^{(i)} = (1 − k_i²)·E^{(i−1)}
//! ```
//!
//! The final `a_j = a_j^{(M)}` are the LP coefficients of the analysis
//! filter `A(z) = 1 + Σ_{i=1}^{M} a_i·z^{−i}`, with `a_0 = 1` implicit.
//! `k_i` are the reflection (PARCOR) coefficients; the spec (clause 3.3)
//! uses the **second-order** reflection pair to derive the adaptive
//! bandwidth-expansion factors of the perceptual weighting filter.
//!
//! ## Sign convention
//!
//! This module returns the `a_i` in the same convention the decoder's
//! [`crate::lsp_to_lp`] produces and consumes: `A(z) = 1 + Σ a_i z^{−i}`
//! (the coefficients that multiply the *synthesis* filter feedback with
//! a leading `+`). The recursion above yields exactly that convention
//! (a minimum-phase `A(z)` for a positive-definite `r'`).
//!
//! ## Numeric domain
//!
//! `f64` throughout, matching the [`crate::lp_analysis`] autocorrelation
//! output; the resulting `a_i` are returned as `f32` to match the
//! decoder LP-coefficient surface.

use crate::tables::M;

/// Result of the §3.2.2 recursion. `a[i − 1]` holds `a_i` for
/// `i = 1 … M` (with `a_0 = 1` implicit); `reflection[i − 1]` holds the
/// order-`i` reflection coefficient `k_i`; `residual_energy` is the
/// final prediction-error energy `E^{(M)}`.
#[derive(Debug, Clone, PartialEq)]
pub struct LevinsonResult {
    /// LP coefficients `a_1 … a_M` of `A(z) = 1 + Σ a_i z^{−i}`.
    pub a: [f32; M],
    /// Reflection (PARCOR) coefficients `k_1 … k_M`.
    pub reflection: [f32; M],
    /// Final prediction-error energy `E^{(M)}`.
    pub residual_energy: f64,
}

/// Runs the §3.2.2 Levinson-Durbin recursion on the modified
/// autocorrelations `r'(0 … M)` (as produced by
/// [`crate::lp_analysis::analyze_window`]) and returns the LP
/// coefficients, reflection coefficients, and residual energy.
///
/// If `r'(0)` is non-positive (which the clause-3.2.1 `r(0) ≥ 1.0`
/// guard prevents in the normal pipeline) or the recursion becomes
/// ill-conditioned (`E ≤ 0`), the recursion stops early and the
/// remaining coefficients stay at their last stable values — a
/// defensive fallback that keeps the output a valid (if lower-order)
/// filter rather than producing NaNs.
#[must_use]
pub fn levinson(r: &[f64; M + 1]) -> LevinsonResult {
    let mut a = [0.0f64; M + 1]; // a[0] = 1 implicit; a[1..=M] filled.
    let mut reflection = [0.0f64; M];
    let mut energy = r[0];

    if energy <= 0.0 {
        // Degenerate: return the trivial all-zero predictor.
        return LevinsonResult {
            a: [0.0; M],
            reflection: [0.0; M],
            residual_energy: energy,
        };
    }

    for i in 1..=M {
        // Eq (8): k_i = −(r'(i) + Σ_{j=1}^{i−1} a_j·r'(i−j)) / E.
        let mut acc = r[i];
        for j in 1..i {
            acc += a[j] * r[i - j];
        }
        let k = -acc / energy;
        reflection[i - 1] = k;

        // Eq (9): a_j^{(i)} = a_j^{(i−1)} + k·a_{i−j}^{(i−1)} for
        // j = 1 … i−1, computed from a snapshot of the previous order so
        // the in-place update never reads an already-updated coefficient
        // (correct for both even and odd i, including the middle term).
        let prev = a;
        for j in 1..i {
            a[j] = prev[j] + k * prev[i - j];
        }
        a[i] = k;

        // Eq (9) energy update. A reflection magnitude ≥ 1 means the
        // autocorrelation was not positive-definite; stop to avoid a
        // negative / zero energy.
        energy *= 1.0 - k * k;
        if energy <= 0.0 {
            break;
        }
    }

    let a_out: [f32; M] = std::array::from_fn(|i| a[i + 1] as f32);
    let refl_out: [f32; M] = std::array::from_fn(|i| reflection[i] as f32);
    LevinsonResult {
        a: a_out,
        reflection: refl_out,
        residual_energy: energy,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lsp_to_lp::lsp_to_lp;

    /// White noise: `r'(k) = 0` for `k ≥ 1` yields an all-zero predictor
    /// (a flat spectrum) with residual energy equal to `r'(0)`.
    #[test]
    fn white_noise_gives_flat_predictor() {
        let mut r = [0.0f64; M + 1];
        r[0] = 100.0;
        let res = levinson(&r);
        for &ai in &res.a {
            assert!(ai.abs() < 1e-6, "expected flat predictor, got {ai}");
        }
        assert!((res.residual_energy - 100.0).abs() < 1e-9);
    }

    /// A first-order AR process `x(n) = ρ·x(n−1) + e(n)` has
    /// autocorrelation `r(k) = ρ^{|k|}·r(0)`; Levinson-Durbin must
    /// recover `a_1 ≈ −ρ` (analysis filter `1 − ρ·z⁻¹`) with all higher
    /// coefficients ≈ 0.
    #[test]
    fn first_order_ar_recovers_rho() {
        let rho = 0.8f64;
        let r0 = 1.0f64;
        let r: [f64; M + 1] = std::array::from_fn(|k| r0 * rho.powi(k as i32));
        let res = levinson(&r);
        assert!(
            (f64::from(res.a[0]) - (-rho)).abs() < 1e-6,
            "a_1={} expected {}",
            res.a[0],
            -rho
        );
        for &ai in &res.a[1..] {
            assert!(ai.abs() < 1e-6, "higher coefficient should vanish: {ai}");
        }
        // The first reflection coefficient equals −ρ for an AR(1).
        assert!((f64::from(res.reflection[0]) - (-rho)).abs() < 1e-6);
        assert!((res.residual_energy - (1.0 - rho * rho)).abs() < 1e-6);
    }

    /// A positive-definite autocorrelation yields a stable (all
    /// reflection magnitudes < 1) minimum-phase filter, and the residual
    /// energy is positive and no larger than `r'(0)`.
    #[test]
    fn stable_minimum_phase() {
        // Autocorrelation of a low-pass process (decaying + oscillating).
        let r: [f64; M + 1] =
            std::array::from_fn(|k| (0.9f64).powi(k as i32) * (0.3 * k as f64).cos());
        let res = levinson(&r);
        for (i, &k) in res.reflection.iter().enumerate() {
            assert!(
                k.abs() < 1.0,
                "reflection {i} magnitude ≥ 1 (not min-phase): {k}"
            );
        }
        assert!(res.residual_energy > 0.0);
        assert!(res.residual_energy <= r[0] + 1e-9);
    }

    /// Round-trip against the decoder's §3.2.6 LSP→LP conversion is
    /// deferred to the §3.2.3 LP→LSP module; here we cross-check that the
    /// LP coefficients Levinson produces, when re-expressed as a
    /// synthesis filter, reproduce the input autocorrelation via the
    /// normal-equations residual (a direct algebraic identity of the
    /// recursion): `r'(0) + Σ a_i·r'(i) = E^{(M)}`.
    #[test]
    fn normal_equations_residual_identity() {
        let r: [f64; M + 1] =
            std::array::from_fn(|k| (0.85f64).powi(k as i32) * (0.2 * k as f64 + 1.0).cos());
        let res = levinson(&r);
        let mut e = r[0];
        for (i, &ai) in res.a.iter().enumerate() {
            e += f64::from(ai) * r[i + 1];
        }
        assert!(
            (e - res.residual_energy).abs() < 1e-6,
            "normal-equations residual {e} vs E^(M) {}",
            res.residual_energy
        );
    }

    /// Sanity that the produced `a_i` are consumable by the existing
    /// decoder path: feeding a synthesized LSP set through `lsp_to_lp`
    /// and (separately) Levinson on the same order both yield
    /// finite 10-tap vectors of the same length/type.
    #[test]
    fn output_shape_matches_decoder_surface() {
        let q: [f32; M] = std::array::from_fn(|i| (0.95 - 0.19 * i as f32).clamp(-0.99, 0.99));
        let a_dec = lsp_to_lp(&q);
        let r: [f64; M + 1] = std::array::from_fn(|k| (0.8f64).powi(k as i32));
        let a_enc = levinson(&r).a;
        assert_eq!(a_dec.len(), a_enc.len());
        assert!(a_enc.iter().all(|v| v.is_finite()));
    }
}
