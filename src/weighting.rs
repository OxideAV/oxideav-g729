//! ITU-T G.729 §3.3 perceptual weighting + §3.5 weighted-synthesis-filter
//! impulse response.
//!
//! These are encoder-analysis building blocks. The analysis-by-synthesis
//! search of the adaptive and fixed codebooks (§3.7, §3.8) is driven by
//! the impulse response `h(n)` of the weighted synthesis filter
//! `W(z)/Â(z)`. This module provides that `h(n)` (§3.5) together with the
//! perceptual weighting filter `W(z)` it is derived from (§3.3).
//!
//! ## §3.3 Perceptual weighting filter (eq 27)
//!
//! ```text
//!            A(z/γ1)     1 + Σ_{i=1..10} γ1^i · a_i · z^-i
//!   W(z) = ---------- = -----------------------------------
//!            A(z/γ2)     1 + Σ_{i=1..10} γ2^i · a_i · z^-i
//! ```
//!
//! `a_i` are the **unquantised** LP coefficients. The pair (γ1, γ2) is
//! adapted once per 10 ms frame to the spectral shape of the input:
//!
//! - The spectrum is classified flat / tilted via a hysteresis decision
//!   (eq 30) on log-area-ratio coefficients `o_1, o_2` (eq 28), derived
//!   from the first two reflection coefficients (a by-product of the
//!   Levinson-Durbin recursion, §3.2.2).
//! - **Flat** (`flat = 1`): `γ1 = 0.94`, `γ2 = 0.6`.
//! - **Tilted** (`flat = 0`): `γ1 = 0.98`; `γ2 = −6.0·d_min + 1.0`
//!   bounded to `[0.4, 0.7]` (eq 31, 32), where `d_min` is the smallest
//!   gap between two successive LSP frequencies of the current subframe.
//!
//! ## §3.5 Impulse response of `W(z)/Â(z)`
//!
//! Per the §3.5 prose, `h(n)` is computed for each subframe by filtering
//! the coefficients of `A(z/γ1)` — i.e. the sequence
//! `[1, γ1·a_1, γ1²·a_2, …, γ1^10·a_10]` — extended by zeros, through the
//! two all-pole filters `1/Â(z)` and `1/A(z/γ2)`, where `Â(z)` is the
//! **quantised** LP inverse filter used for synthesis.
//!
//! All-pole filtering of an input `x(n)` through `1/B(z)`
//! (`B(z) = 1 + Σ b_k z^-k`) is the standard direct-form recursion
//! `y(n) = x(n) − Σ_{k≥1} b_k · y(n-k)`, with zero initial state.

use crate::{LPC_ORDER, SUBFRAME_SAMPLES};

/// Weight factor γ1 for a spectrum classified flat (§3.3).
pub const GAMMA1_FLAT: f32 = 0.94;
/// Weight factor γ2 for a spectrum classified flat (§3.3).
pub const GAMMA2_FLAT: f32 = 0.6;
/// Weight factor γ1 for a spectrum classified tilted (§3.3).
pub const GAMMA1_TILTED: f32 = 0.98;
/// Lower bound on the adapted γ2 for a tilted spectrum (eq 32).
pub const GAMMA2_MIN: f32 = 0.4;
/// Upper bound on the adapted γ2 for a tilted spectrum (eq 32).
pub const GAMMA2_MAX: f32 = 0.7;

/// Length of the impulse response retained for the codebook search:
/// one subframe (40 samples), per §3.5 / §3.7 (the search correlations
/// run over the subframe).
pub const H_LEN: usize = SUBFRAME_SAMPLES;

/// Bandwidth-expand a direct-form LP filter `A(z)` into `A(z/γ)`
/// (ITU-T G.729 §3.3, denominator/numerator of eq 27).
///
/// `a[0]` is the leading `1.0`; `a[k]` (k = 1..=10) are the LP
/// coefficients in `A(z) = 1 + Σ a_k z^-k`. The expanded filter has
/// coefficients `a[k]·γ^k`.
pub fn bandwidth_expand(a: &[f32; LPC_ORDER + 1], gamma: f32) -> [f32; LPC_ORDER + 1] {
    let mut out = [0.0f32; LPC_ORDER + 1];
    out[0] = a[0];
    let mut g = 1.0f32;
    for k in 1..=LPC_ORDER {
        g *= gamma;
        out[k] = a[k] * g;
    }
    out
}

/// Log-area-ratio coefficient from a reflection coefficient `k`
/// (ITU-T G.729 §3.3 eq 28): `o = log((1 + k) / (1 - k))`.
///
/// `k` is clamped just inside `(-1, 1)` so the ratio stays finite even
/// for a marginally-unstable reflection coefficient.
pub fn log_area_ratio(k: f32) -> f32 {
    let k = k.clamp(-0.999_999, 0.999_999);
    ((1.0 + k) / (1.0 - k)).ln()
}

/// Smallest gap between two successive LSP frequencies (ITU-T G.729 §3.3
/// eq 31): `d_min = min_{i=1..9} (ω_{i+1} − ω_i)`.
///
/// `lsp` holds the line spectral pairs in the **cosine domain**
/// (`lsp[k] = cos(ω_k)`), strictly decreasing as produced by the
/// encoder's `lpc_to_lsp`. The successive distance is a frequency-domain
/// (radian) quantity, so each entry is mapped back through `acos` before
/// differencing. `acos` is monotonically decreasing, so the
/// strictly-decreasing cosine sequence maps to a strictly-increasing
/// frequency sequence and every gap is non-negative.
pub fn lsp_min_distance(lsp: &[f32; LPC_ORDER]) -> f32 {
    let mut prev = lsp[0].clamp(-1.0, 1.0).acos();
    let mut d_min = f32::INFINITY;
    for k in 1..LPC_ORDER {
        let w = lsp[k].clamp(-1.0, 1.0).acos();
        let gap = w - prev;
        if gap < d_min {
            d_min = gap;
        }
        prev = w;
    }
    d_min
}

/// Adapt γ2 to the strength of the resonances in the LP synthesis filter
/// for a *tilted* spectrum (ITU-T G.729 §3.3 eq 32):
/// `γ2 = −6.0·d_min + 1.0`, bounded to `[0.4, 0.7]`.
///
/// A strong resonance produces a small `d_min` (two LSPs close together),
/// which drives γ2 toward the upper bound.
pub fn adapt_gamma2(d_min: f32) -> f32 {
    (-6.0 * d_min + 1.0).clamp(GAMMA2_MIN, GAMMA2_MAX)
}

/// Spectral-flatness classifier with hysteresis (ITU-T G.729 §3.3 eq 30).
///
/// Given the two interpolated LAR coefficients `o1`, `o2` of the current
/// subframe and the previous subframe's `flat` flag, returns the new
/// `flat` flag:
///
/// ```text
///   flat = 0   if o1 < -1.74 and o2 > 0.65 and prev_flat == 1
///   flat = 1   if (o1 > -1.52 or o2 < 0.43) and prev_flat == 0
///   flat = prev_flat   otherwise
/// ```
pub fn classify_flat(o1: f32, o2: f32, prev_flat: bool) -> bool {
    if prev_flat && o1 < -1.74 && o2 > 0.65 {
        false
    } else if !prev_flat && (o1 > -1.52 || o2 < 0.43) {
        true
    } else {
        prev_flat
    }
}

/// Choose the perceptual-weighting (γ1, γ2) pair for a subframe
/// (ITU-T G.729 §3.3): flat → (0.94, 0.6); tilted → (0.98, adapted).
pub fn weighting_gammas(flat: bool, lsp: &[f32; LPC_ORDER]) -> (f32, f32) {
    if flat {
        (GAMMA1_FLAT, GAMMA2_FLAT)
    } else {
        (GAMMA1_TILTED, adapt_gamma2(lsp_min_distance(lsp)))
    }
}

/// Filter an input sequence through an all-pole filter `1/B(z)`,
/// `B(z) = 1 + Σ_{k=1..10} b_k z^-k`, with zero initial state.
///
/// Direct-form recursion: `y(n) = x(n) − Σ_{k≥1} b_k · y(n-k)`.
fn all_pole_filter(x: &[f32; H_LEN], b: &[f32; LPC_ORDER + 1]) -> [f32; H_LEN] {
    let mut y = [0.0f32; H_LEN];
    for n in 0..H_LEN {
        let mut acc = x[n];
        for k in 1..=LPC_ORDER {
            if n >= k {
                acc -= b[k] * y[n - k];
            }
        }
        y[n] = acc;
    }
    y
}

/// Impulse response `h(n)` of the weighted synthesis filter
/// `W(z)/Â(z) = A(z/γ1) / [Â(z)·A(z/γ2)]` (ITU-T G.729 §3.5).
///
/// `a_q` is the **quantised** per-subframe LP filter `Â(z)` (used for
/// synthesis); `a_unq` is the **unquantised** LP filter the weighting
/// filter `W(z)` is built from (§3.3). The (γ1, γ2) pair is the adapted
/// perceptual-weighting pair for this subframe.
///
/// Computed per §3.5 by filtering the coefficients of `A(z/γ1)` extended
/// by zeros through `1/Â(z)` then `1/A(z/γ2)`. The result is truncated to
/// [`H_LEN`] (one subframe).
pub fn weighted_synthesis_impulse_response(
    a_q: &[f32; LPC_ORDER + 1],
    a_unq: &[f32; LPC_ORDER + 1],
    gamma1: f32,
    gamma2: f32,
) -> [f32; H_LEN] {
    // Numerator A(z/γ1): its coefficients form the input signal, zero-padded.
    let aw1 = bandwidth_expand(a_unq, gamma1);
    let mut input = [0.0f32; H_LEN];
    input[..=LPC_ORDER].copy_from_slice(&aw1[..=LPC_ORDER]);
    // 1/Â(z): all-pole synthesis filter (quantised).
    let stage1 = all_pole_filter(&input, a_q);
    // 1/A(z/γ2): all-pole denominator of the weighting filter (unquantised,
    // bandwidth-expanded).
    let aw2 = bandwidth_expand(a_unq, gamma2);
    all_pole_filter(&stage1, &aw2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_a() -> [f32; LPC_ORDER + 1] {
        // A(z) = 1 (all higher coefficients zero) — a unit-gain all-pass.
        let mut a = [0.0f32; LPC_ORDER + 1];
        a[0] = 1.0;
        a
    }

    #[test]
    fn bandwidth_expand_scales_by_gamma_powers() {
        let mut a = [0.0f32; LPC_ORDER + 1];
        a[0] = 1.0;
        for k in 1..=LPC_ORDER {
            a[k] = 1.0; // every coefficient = 1 so γ^k is isolated
        }
        let g = 0.9f32;
        let out = bandwidth_expand(&a, g);
        assert_eq!(out[0], 1.0);
        let mut expect = 1.0f32;
        for k in 1..=LPC_ORDER {
            expect *= g;
            assert!(
                (out[k] - expect).abs() < 1e-6,
                "coeff {k}: got {} want {expect}",
                out[k]
            );
        }
    }

    #[test]
    fn bandwidth_expand_identity_for_gamma_one() {
        let mut a = [0.0f32; LPC_ORDER + 1];
        a[0] = 1.0;
        for k in 1..=LPC_ORDER {
            a[k] = -0.1 * k as f32;
        }
        let out = bandwidth_expand(&a, 1.0);
        for k in 0..=LPC_ORDER {
            assert!((out[k] - a[k]).abs() < 1e-6, "coeff {k} changed at γ=1");
        }
    }

    #[test]
    fn log_area_ratio_sign_and_zero() {
        // k = 0 → o = 0; positive k → positive o; symmetric in sign.
        assert!(log_area_ratio(0.0).abs() < 1e-6);
        assert!(log_area_ratio(0.5) > 0.0);
        assert!(log_area_ratio(-0.5) < 0.0);
        assert!((log_area_ratio(0.5) + log_area_ratio(-0.5)).abs() < 1e-6);
    }

    #[test]
    fn log_area_ratio_clamps_near_unit() {
        // Must stay finite even at k = ±1.
        assert!(log_area_ratio(1.0).is_finite());
        assert!(log_area_ratio(-1.0).is_finite());
        assert!(log_area_ratio(1.0) > 0.0);
    }

    #[test]
    fn lsp_min_distance_uniform_spread() {
        // Uniformly spaced LSP frequencies → equal gaps; d_min equals the
        // common step within float tolerance.
        let mut lsp = [0.0f32; LPC_ORDER];
        let step = core::f32::consts::PI / (LPC_ORDER as f32 + 1.0);
        for k in 0..LPC_ORDER {
            lsp[k] = (step * (k as f32 + 1.0)).cos();
        }
        let d = lsp_min_distance(&lsp);
        assert!((d - step).abs() < 1e-4, "d_min {d} should ≈ step {step}");
    }

    #[test]
    fn lsp_min_distance_detects_close_pair() {
        // One pair of LSPs placed very close → small d_min.
        let mut lsp = [0.0f32; LPC_ORDER];
        let step = core::f32::consts::PI / (LPC_ORDER as f32 + 1.0);
        for k in 0..LPC_ORDER {
            lsp[k] = (step * (k as f32 + 1.0)).cos();
        }
        // Squeeze lsp[5] toward lsp[4] in the frequency domain.
        let w4 = lsp[4].acos();
        lsp[5] = (w4 + 0.01).cos();
        let d = lsp_min_distance(&lsp);
        assert!(d < 0.02, "close pair should drive d_min small: {d}");
        assert!(d >= 0.0, "d_min must be non-negative: {d}");
    }

    #[test]
    fn adapt_gamma2_bounds_and_resonance() {
        // d_min = 0 (strong resonance) → upper bound 0.7.
        assert!((adapt_gamma2(0.0) - 0.7).abs() < 1e-6);
        // Large d_min (very flat) → lower bound 0.4.
        assert!((adapt_gamma2(1.0) - 0.4).abs() < 1e-6);
        // Mid value follows the linear law before clamping:
        // -6*0.1 + 1.0 = 0.4 (right at the lower bound).
        assert!((adapt_gamma2(0.1) - 0.4).abs() < 1e-6);
        // -6*0.08 + 1.0 = 0.52, inside the band.
        assert!((adapt_gamma2(0.08) - 0.52).abs() < 1e-5);
    }

    #[test]
    fn classify_flat_hysteresis() {
        // From flat (prev=1): only flips to tilted when BOTH strong
        // conditions hold (o1 < -1.74 AND o2 > 0.65).
        assert!(!classify_flat(-2.0, 0.7, true)); // flips to tilted (0)
        assert!(classify_flat(-2.0, 0.5, true)); // o2 not > 0.65 → stays flat
        assert!(classify_flat(-1.0, 0.7, true)); // o1 not < -1.74 → stays flat
                                                 // From tilted (prev=0): flips to flat when (o1 > -1.52 OR o2 < 0.43).
        assert!(classify_flat(-1.0, 0.9, false)); // o1 > -1.52 → flat
        assert!(classify_flat(-2.0, 0.2, false)); // o2 < 0.43 → flat
        assert!(!classify_flat(-2.0, 0.9, false)); // neither → stays tilted
    }

    #[test]
    fn weighting_gammas_flat_branch() {
        let lsp = {
            let mut l = [0.0f32; LPC_ORDER];
            let step = core::f32::consts::PI / (LPC_ORDER as f32 + 1.0);
            for k in 0..LPC_ORDER {
                l[k] = (step * (k as f32 + 1.0)).cos();
            }
            l
        };
        let (g1, g2) = weighting_gammas(true, &lsp);
        assert_eq!(g1, GAMMA1_FLAT);
        assert_eq!(g2, GAMMA2_FLAT);
        let (g1t, g2t) = weighting_gammas(false, &lsp);
        assert_eq!(g1t, GAMMA1_TILTED);
        assert!((GAMMA2_MIN..=GAMMA2_MAX).contains(&g2t));
    }

    #[test]
    fn impulse_response_starts_at_unity() {
        // For Â(z) = A(z) = 1 (trivial filters), h(n) = A(z/γ1) coeffs
        // through two unit filters = the bandwidth-expanded numerator,
        // so h[0] = 1.0 exactly and h[k>10] = 0.
        let a = flat_a();
        let h = weighted_synthesis_impulse_response(&a, &a, 0.94, 0.6);
        assert!(
            (h[0] - 1.0).abs() < 1e-6,
            "h[0] should be 1.0, got {}",
            h[0]
        );
        for n in (LPC_ORDER + 1)..H_LEN {
            assert!(h[n].abs() < 1e-6, "h[{n}] should be ~0 for trivial filters");
        }
    }

    #[test]
    fn impulse_response_h0_is_one_for_general_filters() {
        // h[0] = aw1[0] / a_q[0] / aw2[0] = 1/1/1 = 1, independent of the
        // tap values, because all three leading coefficients are 1.0 and
        // the all-pole recursion has zero history at n = 0.
        let mut a_q = [0.0f32; LPC_ORDER + 1];
        let mut a_unq = [0.0f32; LPC_ORDER + 1];
        a_q[0] = 1.0;
        a_unq[0] = 1.0;
        for k in 1..=LPC_ORDER {
            a_q[k] = -0.05 * k as f32;
            a_unq[k] = -0.04 * k as f32;
        }
        let h = weighted_synthesis_impulse_response(&a_q, &a_unq, 0.94, 0.6);
        assert!((h[0] - 1.0).abs() < 1e-6, "h[0] = {}", h[0]);
        for v in h.iter() {
            assert!(v.is_finite(), "impulse response must stay finite");
        }
    }

    #[test]
    fn impulse_response_matches_two_stage_reference_path() {
        // Cross-check the fused routine against an explicit two-pass
        // reference: build A(z/γ1) coefficients, filter through 1/Â(z),
        // then through 1/A(z/γ2), all with zero state.
        let mut a_q = [0.0f32; LPC_ORDER + 1];
        let mut a_unq = [0.0f32; LPC_ORDER + 1];
        a_q[0] = 1.0;
        a_unq[0] = 1.0;
        for k in 1..=LPC_ORDER {
            a_q[k] = -0.03 * k as f32 + 0.01;
            a_unq[k] = -0.02 * k as f32;
        }
        let (g1, g2) = (0.98f32, 0.55f32);

        // Reference: manual direct-form recursions.
        let aw1 = bandwidth_expand(&a_unq, g1);
        let aw2 = bandwidth_expand(&a_unq, g2);
        let mut x = [0.0f32; H_LEN];
        x[..=LPC_ORDER].copy_from_slice(&aw1[..=LPC_ORDER]);
        let mut y1 = [0.0f32; H_LEN];
        for n in 0..H_LEN {
            let mut acc = x[n];
            for k in 1..=LPC_ORDER {
                if n >= k {
                    acc -= a_q[k] * y1[n - k];
                }
            }
            y1[n] = acc;
        }
        let mut y2 = [0.0f32; H_LEN];
        for n in 0..H_LEN {
            let mut acc = y1[n];
            for k in 1..=LPC_ORDER {
                if n >= k {
                    acc -= aw2[k] * y2[n - k];
                }
            }
            y2[n] = acc;
        }

        let h = weighted_synthesis_impulse_response(&a_q, &a_unq, g1, g2);
        for n in 0..H_LEN {
            assert!(
                (h[n] - y2[n]).abs() < 1e-5,
                "h[{n}] {} != reference {}",
                h[n],
                y2[n]
            );
        }
    }
}
