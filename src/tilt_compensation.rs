//! §4.2.3 **tilt compensation** `H_t(z)` — the third filter of the
//! decoder's §4.2 post-processing adaptive postfilter cascade.
//!
//! The §4.2 cascade is (clause 4.2): long-term postfilter `H_p(z)`
//! (§4.2.1) → short-term postfilter `H_f(z)` (§4.2.2) → **tilt
//! compensation `H_t(z)` (§4.2.3)** → adaptive gain control (§4.2.4) →
//! output high-pass + ×2 upscaling (§4.2.5). The round-313
//! [`crate::short_term_postfilter`] module wired §4.2.2 and the
//! round-298 [`crate::post_process`] module wired §4.2.5; this module
//! wires §4.2.3, whose only external input is the short-term
//! postfilter's impulse response `h_f(n)` (already produced by
//! [`crate::short_term_postfilter::ShortTermPostfilter::impulse_response`]).
//!
//! ## Spec source — clause 4.2.3, equations (86)/(87) (06/2012 Rec.)
//!
//! Clause 4.2.3 (transcribed from the EPUB prose): "The filter `H_t(z)`
//! compensates for the tilt in the short-term postfilter `H_f(z)` and is
//! given by [eq (86)] where [`k1'`] is a tilt factor being the first
//! reflection coefficient calculated from `h_f(n)` with [eq (87)]. The
//! gain term `g_t = 1 − |k1'|`… compensates for the decreasing effect of
//! `g_f` in `H_f(z)`. … Two values for `γ_t` are used depending on the
//! sign of [`k1'`]. If [`k1'`] is negative, `γ_t = 0.9`, and if [`k1'`]
//! is positive, `γ_t = 0.2`."
//!
//! Equation (86) (raster `images/eq86.jpg`) is the first-order FIR
//! transfer function
//!
//! ```text
//!            1
//! H_t(z) = ─────·(1 + γ_t·k1'·z⁻¹)
//!           g_t
//! ```
//!
//! and equation (87) (raster `images/eq87.jpg`) defines the tilt factor
//! from the autocorrelation of `h_f(n)`:
//!
//! ```text
//!            r_h(1)               19−i
//! k1' = − ──────────  ,  r_h(i) =  Σ  h_f(j)·h_f(j+i)
//!            r_h(0)               j=0
//! ```
//!
//! The gain term is `g_t = 1 − |γ_t·k1'|`. (The prose abbreviates this
//! as "`g_t = 1 − |k1'|`", but the surrounding text and the eq (86)
//! structure make `γ_t·k1'` the quantity whose tilt-correction magnitude
//! is removed — the same `γ_t·k1'` product that scales `z⁻¹` in the
//! numerator. The `g_t > 0` guard below makes the two readings agree
//! whenever `|γ_t·k1'| < 1`, which holds for every `γ_t ∈ {0.2, 0.9}`
//! and the `|k1'| ≤ 1` reflection-coefficient bound.)
//!
//! ## Per-subframe procedure
//!
//! Each 5 ms subframe (clause 4.2: "the postfilter coefficients are
//! updated every 5 ms subframe"):
//!
//! 1. Form `h_f(n)`, `n = 0 … 19`, from the current quantised LP
//!    coefficients `â_i` (the short-term postfilter impulse response).
//! 2. `r_h(0) = Σ h_f(j)²`, `r_h(1) = Σ h_f(j)·h_f(j+1)`.
//! 3. `k1' = −r_h(1)/r_h(0)` (eq (87)); `γ_t = 0.9` if `k1' < 0`, else
//!    `0.2`.
//! 4. `g_t = 1 − |γ_t·k1'|`.
//! 5. Filter the §4.2.2 short-term-postfilter output `t(n)` through the
//!    eq (86) FIR: `sf(n) = (1/g_t)·(t(n) + γ_t·k1'·t(n−1))`, carrying
//!    the single input tap `t(n−1)` across subframes.
//!
//! ## State (clause 4.3 init)
//!
//! Per clause 4.3, "all static encoder and decoder variables should be
//! initialized to zero, except the variables listed in Table 9". The
//! single FIR input tap `t(n−1)` does not appear in Table 9, so it
//! starts zeroed and carries across subframes (the §4.2.3 coefficients
//! `k1'` / `γ_t` / `g_t` are recomputed each subframe; the signal memory
//! is continuous).
//!
//! ## What this module does NOT do
//!
//! It is the §4.2.3 stage only. The §4.2.4 adaptive gain control that
//! follows it ([`crate::adaptive_gain_control`]) and the §4.2.1
//! long-term postfilter that precedes the short-term stage are separate
//! modules; they slot around this one unchanged.

use crate::fixed_codebook::SUBFRAME_SIZE;
use crate::short_term_postfilter::{ShortTermPostfilter, GF_IMPULSE_LEN};
use crate::tables::M;

/// §4.2.3 tilt factor `γ_t = 0.9` used when `k1'` is **negative**
/// (clause 4.2.3: "If [k1'] is negative, γ_t = 0.9").
pub const GAMMA_T_NEG: f32 = 0.9;

/// §4.2.3 tilt factor `γ_t = 0.2` used when `k1'` is **positive**
/// (clause 4.2.3: "if [k1'] is positive, γ_t = 0.2").
pub const GAMMA_T_POS: f32 = 0.2;

/// The eq (86)/(87) per-subframe tilt coefficients derived from the
/// short-term postfilter impulse response `h_f(n)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TiltCoefficients {
    /// eq (87) tilt factor `k1' = −r_h(1)/r_h(0)` (the first reflection
    /// coefficient of `h_f`). Zero when `r_h(0) = 0` (a degenerate
    /// all-zero impulse response).
    pub k1: f32,
    /// §4.2.3 sign-selected `γ_t` (`0.9` if `k1' < 0`, else `0.2`).
    pub gamma_t: f32,
    /// eq (86) gain term `g_t = 1 − |γ_t·k1'|`.
    pub g_t: f32,
}

impl TiltCoefficients {
    /// Derive the eq (86)/(87) coefficients from a 20-sample short-term
    /// postfilter impulse response `h_f(n)`.
    ///
    /// `r_h(0) = Σ h_f(j)²`, `r_h(1) = Σ h_f(j)·h_f(j+1)`,
    /// `k1' = −r_h(1)/r_h(0)`. When `r_h(0)` is non-positive (an
    /// all-zero impulse response — impossible for a real `Â`, since
    /// `h_f(0) = 1`, but guarded) `k1'` is taken as zero so the tilt
    /// filter degenerates to the identity (`γ_t·k1' = 0`, `g_t = 1`).
    #[must_use]
    pub fn from_impulse_response(h: &[f32; GF_IMPULSE_LEN]) -> Self {
        let mut r0 = 0.0f32;
        let mut r1 = 0.0f32;
        for j in 0..GF_IMPULSE_LEN {
            r0 += h[j] * h[j];
            if j + 1 < GF_IMPULSE_LEN {
                r1 += h[j] * h[j + 1];
            }
        }
        let k1 = if r0 > 0.0 { -r1 / r0 } else { 0.0 };
        let gamma_t = if k1 < 0.0 { GAMMA_T_NEG } else { GAMMA_T_POS };
        let g_t = 1.0 - (gamma_t * k1).abs();
        Self { k1, gamma_t, g_t }
    }

    /// Derive the coefficients directly from a subframe's quantised LP
    /// coefficients `â_i` (slots `0 … 9` hold `â_1 … â_10`), computing
    /// the short-term postfilter impulse response `h_f(n)` internally.
    #[must_use]
    pub fn from_lp(a: &[f32; M]) -> Self {
        Self::from_impulse_response(&ShortTermPostfilter::impulse_response(a))
    }
}

/// Stateful §4.2.3 tilt-compensation filter `H_t(z)` (eqs (86)/(87)).
///
/// Owns the single FIR input tap `t(n−1)`, zero-initialised per clause
/// 4.3 and carried across subframes. The per-subframe coefficients
/// `k1'` / `γ_t` / `g_t` are recomputed from the supplied impulse
/// response (or LP coefficients) each [`Self::filter_subframe`] call.
#[derive(Debug, Clone, Default)]
pub struct TiltCompensation {
    /// FIR input history `t(n−1)`. Zero-init per clause 4.3.
    t_prev: f32,
}

impl TiltCompensation {
    /// Build the filter with the clause-4.3 all-zero start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self { t_prev: 0.0 }
    }

    /// Borrow the carried-over FIR input tap `t(n−1)` for inspection /
    /// tests.
    #[must_use]
    pub fn t_prev(&self) -> f32 {
        self.t_prev
    }

    /// Filter one 40-sample subframe `t(n)` (the §4.2.2 short-term-
    /// postfilter output) through `H_t(z)` (eq (86)) using the supplied
    /// per-subframe [`TiltCoefficients`], advancing the carried-over
    /// input tap.
    ///
    /// `sf(n) = (1/g_t)·(t(n) + γ_t·k1'·t(n−1))`. When `g_t` is
    /// non-positive (cannot happen for `|γ_t·k1'| < 1`, but guarded for
    /// a hand-built degenerate input) the `1/g_t` normalisation is
    /// skipped (scale 1.0) rather than dividing by zero / flipping sign.
    #[must_use]
    pub fn filter_subframe_with(
        &mut self,
        t: &[f32; SUBFRAME_SIZE],
        c: &TiltCoefficients,
    ) -> [f32; SUBFRAME_SIZE] {
        let inv_gt = if c.g_t > 0.0 { 1.0 / c.g_t } else { 1.0 };
        let b1 = c.gamma_t * c.k1; // z⁻¹ numerator coefficient
        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &tn) in t.iter().enumerate() {
            out[n] = inv_gt * (tn + b1 * self.t_prev);
            self.t_prev = tn;
        }
        out
    }

    /// Convenience: derive the [`TiltCoefficients`] from the subframe's
    /// quantised LP coefficients `â_i`, then filter the subframe.
    #[must_use]
    pub fn filter_subframe(
        &mut self,
        t: &[f32; SUBFRAME_SIZE],
        a: &[f32; M],
    ) -> [f32; SUBFRAME_SIZE] {
        let c = TiltCoefficients::from_lp(a);
        self.filter_subframe_with(t, &c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A flat (all-zero `â_i`) short-term postfilter has impulse response
    /// `h_f = [1, 0, …]`, so `r_h(1) = 0`, `k1' = 0`, `γ_t = 0.2`
    /// (k1' ≥ 0), and `g_t = 1`: the tilt filter is the identity.
    #[test]
    fn flat_filter_is_identity() {
        let a = [0.0f32; M];
        let c = TiltCoefficients::from_lp(&a);
        assert!((c.k1 - 0.0).abs() < 1e-6);
        assert!((c.gamma_t - GAMMA_T_POS).abs() < 1e-9);
        assert!((c.g_t - 1.0).abs() < 1e-6);

        let mut tc = TiltCompensation::new();
        let t: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let out = tc.filter_subframe(&t, &a);
        for n in 0..SUBFRAME_SIZE {
            assert!((out[n] - t[n]).abs() < 1e-6, "n={n}");
        }
    }

    /// New filter starts zeroed (clause 4.3).
    #[test]
    fn new_filter_has_zero_state() {
        let tc = TiltCompensation::new();
        assert_eq!(tc.t_prev(), 0.0);
    }

    /// The tilt factors are the spec values (clause 4.2.3).
    #[test]
    fn gamma_factors_match_spec() {
        assert!((GAMMA_T_NEG - 0.9).abs() < 1e-9);
        assert!((GAMMA_T_POS - 0.2).abs() < 1e-9);
    }

    /// `k1' = −r_h(1)/r_h(0)` and the `γ_t` sign rule, worked by hand for
    /// a known impulse response. With `h_f = [1, 0.5, 0, …]`,
    /// `r_h(0) = 1.25`, `r_h(1) = 0.5`, so `k1' = −0.4` (negative ⇒
    /// γ_t = 0.9) and `g_t = 1 − |0.9·(−0.4)| = 1 − 0.36 = 0.64`.
    #[test]
    fn tilt_factor_hand_computation_negative_k1() {
        let mut h = [0.0f32; GF_IMPULSE_LEN];
        h[0] = 1.0;
        h[1] = 0.5;
        let c = TiltCoefficients::from_impulse_response(&h);
        assert!((c.k1 - (-0.4)).abs() < 1e-6, "k1={}", c.k1);
        assert!((c.gamma_t - GAMMA_T_NEG).abs() < 1e-9);
        assert!((c.g_t - 0.64).abs() < 1e-6, "g_t={}", c.g_t);
    }

    /// A positive `k1'` selects `γ_t = 0.2`. With `h_f = [1, −0.5, 0, …]`,
    /// `r_h(1) = −0.5`, `k1' = +0.4`, `γ_t = 0.2`,
    /// `g_t = 1 − |0.2·0.4| = 0.92`.
    #[test]
    fn tilt_factor_hand_computation_positive_k1() {
        let mut h = [0.0f32; GF_IMPULSE_LEN];
        h[0] = 1.0;
        h[1] = -0.5;
        let c = TiltCoefficients::from_impulse_response(&h);
        assert!((c.k1 - 0.4).abs() < 1e-6, "k1={}", c.k1);
        assert!((c.gamma_t - GAMMA_T_POS).abs() < 1e-9);
        assert!((c.g_t - 0.92).abs() < 1e-6, "g_t={}", c.g_t);
    }

    /// The eq (86) FIR worked by hand for the first three samples, with
    /// the carried-over `t(n−1)` tap. Coeffs from `h_f = [1, 0.5]` above:
    /// `γ_t·k1' = 0.9·(−0.4) = −0.36`, `1/g_t = 1/0.64`.
    #[test]
    fn fir_matches_hand_recursion() {
        let mut h = [0.0f32; GF_IMPULSE_LEN];
        h[0] = 1.0;
        h[1] = 0.5;
        let c = TiltCoefficients::from_impulse_response(&h);
        let b1 = c.gamma_t * c.k1; // −0.36
        let inv = 1.0 / c.g_t; // 1/0.64

        let mut t = [0.0f32; SUBFRAME_SIZE];
        t[0] = 10.0;
        t[1] = -4.0;
        t[2] = 7.0;

        let mut tc = TiltCompensation::new();
        let out = tc.filter_subframe_with(&t, &c);

        // sf(0) = inv·(t0 + b1·0) ; sf(1) = inv·(t1 + b1·t0) ; …
        assert!((out[0] - inv * (t[0])).abs() < 1e-5);
        assert!((out[1] - inv * (t[1] + b1 * t[0])).abs() < 1e-5);
        assert!((out[2] - inv * (t[2] + b1 * t[1])).abs() < 1e-5);
        // The last input sample is carried into the next subframe.
        assert!((tc.t_prev() - t[SUBFRAME_SIZE - 1]).abs() < 1e-6);
    }

    /// State carries across subframes: filtering an 80-sample stream as
    /// two 40-sample subframes (same coeffs) equals a single continuous
    /// FIR built by hand.
    #[test]
    fn state_carries_across_subframes() {
        let mut h = [0.0f32; GF_IMPULSE_LEN];
        h[0] = 1.0;
        h[1] = 0.3;
        h[2] = -0.1;
        let c = TiltCoefficients::from_impulse_response(&h);
        let b1 = c.gamma_t * c.k1;
        let inv = 1.0 / c.g_t;

        let mut t1 = [0.0f32; SUBFRAME_SIZE];
        let mut t2 = [0.0f32; SUBFRAME_SIZE];
        for n in 0..SUBFRAME_SIZE {
            t1[n] = ((n * 3) % 11) as f32 - 5.0;
            t2[n] = ((n * 7) % 13) as f32 - 6.0;
        }

        let mut tc = TiltCompensation::new();
        let o1 = tc.filter_subframe_with(&t1, &c);
        let o2 = tc.filter_subframe_with(&t2, &c);

        let mut stream = [0.0f32; 2 * SUBFRAME_SIZE];
        stream[..SUBFRAME_SIZE].copy_from_slice(&t1);
        stream[SUBFRAME_SIZE..].copy_from_slice(&t2);
        let mut prev = 0.0f32;
        let mut ref_out = [0.0f32; 2 * SUBFRAME_SIZE];
        for n in 0..2 * SUBFRAME_SIZE {
            ref_out[n] = inv * (stream[n] + b1 * prev);
            prev = stream[n];
        }
        for n in 0..SUBFRAME_SIZE {
            assert!((o1[n] - ref_out[n]).abs() < 1e-4, "sub1 n={n}");
            assert!(
                (o2[n] - ref_out[SUBFRAME_SIZE + n]).abs() < 1e-4,
                "sub2 n={n}"
            );
        }
    }

    /// On a realistic stable LP coefficient set, the derived coefficients
    /// are finite, `|k1'| ≤ 1` (reflection-coefficient bound), and the
    /// output stays finite.
    #[test]
    fn finite_on_realistic_coefficients() {
        let a: [f32; M] = [
            0.9, -0.4, 0.2, -0.1, 0.05, -0.03, 0.02, -0.01, 0.005, -0.002,
        ];
        let c = TiltCoefficients::from_lp(&a);
        assert!(c.k1.is_finite() && c.k1.abs() <= 1.0 + 1e-4, "k1={}", c.k1);
        assert!(c.g_t.is_finite() && c.g_t > 0.0, "g_t={}", c.g_t);

        let mut tc = TiltCompensation::new();
        let t: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32 - 20.0) * 100.0);
        for _ in 0..10 {
            let out = tc.filter_subframe(&t, &a);
            assert!(out.iter().all(|v| v.is_finite()));
        }
    }

    /// Two filters over the same subframe sequence stay in lockstep — all
    /// state is owned, no hidden globals.
    #[test]
    fn deterministic() {
        let mut a = [0.0f32; M];
        a[0] = 0.4;
        a[3] = -0.15;
        let t: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) * 2.0 - 39.0);
        let mut p = TiltCompensation::new();
        let mut q = TiltCompensation::new();
        for _ in 0..4 {
            assert_eq!(p.filter_subframe(&t, &a), q.filter_subframe(&t, &a));
        }
    }
}
