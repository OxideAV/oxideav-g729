//! §4.2.4 **adaptive gain control** — the fourth function of the
//! decoder's §4.2 post-processing, applied after the tilt-compensation
//! filter to match the postfiltered signal's energy back to the
//! reconstructed speech.
//!
//! The §4.2 cascade is (clause 4.2): long-term postfilter `H_p(z)`
//! (§4.2.1) → short-term postfilter `H_f(z)` (§4.2.2) → tilt
//! compensation `H_t(z)` (§4.2.3) → **adaptive gain control (§4.2.4)** →
//! output high-pass + ×2 upscaling (§4.2.5). The round-313
//! [`crate::short_term_postfilter`], round-319
//! [`crate::tilt_compensation`], and round-298 [`crate::post_process`]
//! modules wired §4.2.2 / §4.2.3 / §4.2.5; this module wires §4.2.4,
//! whose inputs are the reconstructed speech `ŝ(n)` (from
//! [`crate::lp_synthesis`]) and the tilt-compensated postfilter output
//! `sf(n)` (from [`crate::tilt_compensation`]).
//!
//! ## Spec source — clause 4.2.4, equations (88)/(89)/(90) (06/2012 Rec.)
//!
//! Clause 4.2.4 (transcribed from the EPUB prose): "Adaptive gain
//! control is used to compensate for gain differences between the
//! reconstructed speech signal `ŝ(n)` and the postfiltered signal
//! `sf(n)`. The gain scaling factor `G` for the present subframe is
//! computed by [eq (88)]. The gain-scaled postfiltered signal `sf′(n)`
//! is given by [eq (89)] where `g(n)` is updated on a sample-by-sample
//! basis and given by [eq (90)]. The initial value of `g(−1) = 1.0` is
//! used. Then for each new subframe, `g(−1)` is set equal to `g(39)` of
//! the previous subframe."
//!
//! Equation (88) (raster `images/eq88.jpg`) — the per-subframe energy
//! ratio over the 40 samples:
//!
//! ```text
//!        Σ_{n=0}^{39} |ŝ(n)|
//! G  =  ──────────────────────
//!        Σ_{n=0}^{39} |sf(n)|
//! ```
//!
//! Equation (89) (raster `images/eq89.jpg`) — the sample-by-sample
//! gain-scaling, `n = 0 … 39`:
//!
//! ```text
//! sf′(n) = g(n)·sf(n)
//! ```
//!
//! Equation (90) (raster `images/eq90.jpg`) — the per-sample smoothing
//! of the gain toward `G`, `n = 0 … 39`:
//!
//! ```text
//! g(n) = 0.85·g(n−1) + 0.15·G
//! ```
//!
//! ## State (clause 4.3 / Table 9 init)
//!
//! Per Table 9, `g(−1)` (referenced from clause 4.2.4) has the non-zero
//! initial value `1.0`. The smoothed gain `g` is the single piece of
//! cross-subframe state: after each subframe, `g(−1)` for the next
//! subframe is `g(39)` of this one (clause 4.2.4), which this module
//! realises by simply carrying the running `g` across calls.
//!
//! ## Numerical guard
//!
//! Eq (88) divides by `Σ |sf(n)|`. On a silent subframe (all-zero
//! `sf(n)`) that denominator is zero; the spec's fixed-point path keeps
//! the previous `G` in that case (no new ratio is formed). This module
//! takes `G = g(−1)` (the carried-over gain) when the denominator is
//! non-positive, so `g(n)` holds steady through silence rather than
//! exploding — the energy match is undefined when there is no
//! postfiltered energy to match against.

use crate::fixed_codebook::SUBFRAME_SIZE;

/// eq (90) gain-smoothing pole weight on the previous sample's gain
/// (`0.85`).
pub const AGC_PREV_WEIGHT: f32 = 0.85;

/// eq (90) gain-smoothing weight on the new target ratio `G` (`0.15`).
pub const AGC_TARGET_WEIGHT: f32 = 0.15;

/// Table 9 initial value of the §4.2.4 smoothed gain `g(−1) = 1.0`.
pub const G_INIT: f32 = 1.0;

/// Stateful §4.2.4 adaptive gain control (eqs (88)/(89)/(90)).
///
/// Owns the running smoothed gain `g`, initialised to [`G_INIT`] (Table
/// 9 `g(−1) = 1.0`) and carried across subframes so that `g(−1)` for the
/// next subframe is `g(39)` of this one (clause 4.2.4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveGainControl {
    /// The running smoothed gain. Holds `g(n)` after the most recent
    /// sample; at a subframe boundary it is `g(39)` of the previous
    /// subframe, i.e. `g(−1)` for the next (clause 4.2.4). Table 9
    /// init `1.0`.
    g: f32,
}

impl Default for AdaptiveGainControl {
    fn default() -> Self {
        Self::new()
    }
}

impl AdaptiveGainControl {
    /// Build the AGC with the Table 9 start-up gain `g(−1) = 1.0`.
    #[must_use]
    pub fn new() -> Self {
        Self { g: G_INIT }
    }

    /// Borrow the current smoothed gain (`g(−1)` for the next subframe).
    #[must_use]
    pub fn gain(&self) -> f32 {
        self.g
    }

    /// eq (88) per-subframe target ratio `G = Σ|ŝ(n)| / Σ|sf(n)|` over
    /// the 40 samples. Returns `None` when the denominator `Σ|sf(n)|` is
    /// non-positive (a silent postfiltered subframe), in which case the
    /// caller holds the previous gain.
    #[must_use]
    pub fn target_ratio(s_hat: &[f32; SUBFRAME_SIZE], sf: &[f32; SUBFRAME_SIZE]) -> Option<f32> {
        let num: f32 = s_hat.iter().map(|v| v.abs()).sum();
        let den: f32 = sf.iter().map(|v| v.abs()).sum();
        if den > 0.0 {
            Some(num / den)
        } else {
            None
        }
    }

    /// Apply §4.2.4 adaptive gain control to one 40-sample subframe.
    ///
    /// `s_hat` is the reconstructed speech `ŝ(n)` (§4.1.6 output, the
    /// energy reference); `sf` is the tilt-compensated postfilter output
    /// `sf(n)` (§4.2.3 output). Returns the gain-scaled `sf′(n)`
    /// (eq (89)). The running gain advances per eq (90) for every
    /// sample and carries into the next subframe.
    #[must_use]
    pub fn process_subframe(
        &mut self,
        s_hat: &[f32; SUBFRAME_SIZE],
        sf: &[f32; SUBFRAME_SIZE],
    ) -> [f32; SUBFRAME_SIZE] {
        // eq (88): G for the present subframe (held at the previous gain
        // when sf is silent — see module "Numerical guard").
        let g_target = Self::target_ratio(s_hat, sf).unwrap_or(self.g);

        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &sfn) in sf.iter().enumerate() {
            // eq (90): g(n) = 0.85·g(n−1) + 0.15·G.
            self.g = AGC_PREV_WEIGHT * self.g + AGC_TARGET_WEIGHT * g_target;
            // eq (89): sf′(n) = g(n)·sf(n).
            out[n] = self.g * sfn;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Table 9 start-up gain.
    #[test]
    fn new_has_table9_unit_gain() {
        let agc = AdaptiveGainControl::new();
        assert_eq!(agc.gain(), 1.0);
        assert_eq!(G_INIT, 1.0);
    }

    /// The smoothing weights are the spec values and sum to 1 (eq (90)).
    #[test]
    fn weights_match_spec() {
        assert!((AGC_PREV_WEIGHT - 0.85).abs() < 1e-9);
        assert!((AGC_TARGET_WEIGHT - 0.15).abs() < 1e-9);
        assert!((AGC_PREV_WEIGHT + AGC_TARGET_WEIGHT - 1.0).abs() < 1e-9);
    }

    /// eq (88) ratio computed by hand: `Σ|ŝ| = 40·2 = 80`,
    /// `Σ|sf| = 40·1 = 40` ⇒ `G = 2`.
    #[test]
    fn target_ratio_hand_computation() {
        let s_hat = [2.0f32; SUBFRAME_SIZE];
        let sf = [1.0f32; SUBFRAME_SIZE];
        let g = AdaptiveGainControl::target_ratio(&s_hat, &sf).unwrap();
        assert!((g - 2.0).abs() < 1e-6);
    }

    /// A silent postfiltered subframe (`Σ|sf| = 0`) yields no ratio.
    #[test]
    fn target_ratio_none_on_silence() {
        let s_hat = [5.0f32; SUBFRAME_SIZE];
        let sf = [0.0f32; SUBFRAME_SIZE];
        assert!(AdaptiveGainControl::target_ratio(&s_hat, &sf).is_none());
    }

    /// When `ŝ = sf` the target ratio is exactly 1, and with the Table 9
    /// init `g(−1) = 1` the gain stays at 1 (eq (90):
    /// `g = 0.85·1 + 0.15·1 = 1`), so the output equals the input.
    #[test]
    fn unity_gain_when_energies_match() {
        let mut agc = AdaptiveGainControl::new();
        let x: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let out = agc.process_subframe(&x, &x);
        for n in 0..SUBFRAME_SIZE {
            assert!((out[n] - x[n]).abs() < 1e-5, "n={n}");
        }
        assert!((agc.gain() - 1.0).abs() < 1e-5);
    }

    /// eq (90) sample-by-sample smoothing worked by hand for the first
    /// three samples. With `G = 2` and `g(−1) = 1`:
    /// `g(0) = 0.85·1 + 0.15·2 = 1.15`,
    /// `g(1) = 0.85·1.15 + 0.3 = 1.2775`,
    /// `g(2) = 0.85·1.2775 + 0.3 = 1.385875`; `sf′(n) = g(n)·sf(n)`.
    #[test]
    fn gain_smoothing_matches_hand_recursion() {
        let s_hat = [2.0f32; SUBFRAME_SIZE]; // Σ = 80
        let sf = [1.0f32; SUBFRAME_SIZE]; // Σ = 40 ⇒ G = 2
        let mut agc = AdaptiveGainControl::new();
        let out = agc.process_subframe(&s_hat, &sf);

        let g0 = 0.85 * 1.0 + 0.15 * 2.0;
        let g1 = 0.85 * g0 + 0.15 * 2.0;
        let g2 = 0.85 * g1 + 0.15 * 2.0;
        assert!((out[0] - g0 * 1.0).abs() < 1e-5, "g0");
        assert!((out[1] - g1 * 1.0).abs() < 1e-5, "g1");
        assert!((out[2] - g2 * 1.0).abs() < 1e-5, "g2");
    }

    /// Over many subframes with a constant target `G`, the smoothed gain
    /// converges to `G` (the eq (90) fixed point: `g = 0.85g + 0.15G ⇒
    /// g = G`).
    #[test]
    fn gain_converges_to_constant_target() {
        let s_hat = [3.0f32; SUBFRAME_SIZE]; // Σ = 120
        let sf = [1.0f32; SUBFRAME_SIZE]; // Σ = 40 ⇒ G = 3
        let mut agc = AdaptiveGainControl::new();
        for _ in 0..50 {
            let _ = agc.process_subframe(&s_hat, &sf);
        }
        assert!((agc.gain() - 3.0).abs() < 1e-3, "g={}", agc.gain());
    }

    /// On a silent postfiltered subframe the gain holds steady (no
    /// division by zero) — `G` falls back to the carried gain.
    #[test]
    fn silent_subframe_holds_gain() {
        let mut agc = AdaptiveGainControl::new();
        // Drive the gain to ~2 first.
        let s_hat = [2.0f32; SUBFRAME_SIZE];
        let sf = [1.0f32; SUBFRAME_SIZE];
        for _ in 0..50 {
            let _ = agc.process_subframe(&s_hat, &sf);
        }
        let g_before = agc.gain();
        // Now a silent postfiltered subframe.
        let silent = [0.0f32; SUBFRAME_SIZE];
        let out = agc.process_subframe(&s_hat, &silent);
        assert!(out.iter().all(|v| v.abs() < 1e-6), "silent sf → silent out");
        // G == g_before ⇒ g stays at g_before (fixed point).
        assert!((agc.gain() - g_before).abs() < 1e-4, "g drifted");
    }

    /// State carries across subframes: the gain at the end of subframe 1
    /// is `g(−1)` for subframe 2 (clause 4.2.4 "g(−1) := g(39)").
    #[test]
    fn gain_carries_across_subframes() {
        let s_hat = [4.0f32; SUBFRAME_SIZE];
        let sf = [1.0f32; SUBFRAME_SIZE];
        let mut agc = AdaptiveGainControl::new();
        let _ = agc.process_subframe(&s_hat, &sf);
        let g_after_sub1 = agc.gain();

        // A fresh AGC whose init gain is g_after_sub1 must produce the
        // identical second subframe.
        let mut cont = AdaptiveGainControl::new();
        let _ = cont.process_subframe(&s_hat, &sf);
        let a = agc.process_subframe(&s_hat, &sf);
        let b = cont.process_subframe(&s_hat, &sf);
        assert_eq!(a, b);
        assert!(g_after_sub1.is_finite());
    }

    /// Two AGCs over the same sequence stay in lockstep — all state is
    /// owned, no hidden globals.
    #[test]
    fn deterministic() {
        let s_hat: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32).sin() * 100.0);
        let sf: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32).cos() * 80.0 + 1.0);
        let mut p = AdaptiveGainControl::new();
        let mut q = AdaptiveGainControl::new();
        for _ in 0..5 {
            assert_eq!(
                p.process_subframe(&s_hat, &sf),
                q.process_subframe(&s_hat, &sf)
            );
        }
    }
}
