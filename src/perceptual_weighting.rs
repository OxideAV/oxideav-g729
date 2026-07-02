//! §3.3 **perceptual weighting** — the adaptive weighting filter
//! `W(z) = A(z/γ₁)/A(z/γ₂)` that shapes the encoder's error criterion
//! and produces the weighted speech `sw(n)` used by the §3.4 open-loop
//! pitch estimation.
//!
//! Sixth encoder-side stage. Consumes the *unquantised* LP coefficients
//! from [`crate::levinson`] plus its 2nd-order reflection-coefficient
//! by-product, and the unquantised subframe LSF from
//! [`crate::lp_to_lsp`] (for the resonance criterion).
//!
//! ## Spec source — clause 3.3, equations (27)–(33)
//!
//! (All transcribed from the EPUB prose + equation rasters
//! `images/eq27.jpg … eq33.jpg`.)
//!
//! * **eq (27)** — `W(z) = A(z/γ₁)/A(z/γ₂)
//!   = (1 + Σ γ₁ⁱ·a_i·z⁻ⁱ) / (1 + Σ γ₂ⁱ·a_i·z⁻ⁱ)`, built from the
//!   **unquantised** `a_i`.
//! * **eq (28)** — the first two reflection coefficients `k₁, k₂` (the
//!   2nd-order Levinson by-product) are converted to log-area-ratio
//!   coefficients `o_i = log((1.0 + k_i)/(1.0 − k_i))`, `i = 1, 2`.
//! * **eq (29)** — per-subframe LAR interpolation: subframe 1 uses
//!   `0.5·o_i^(previous) + 0.5·o_i^(current)`, subframe 2 uses the
//!   current frame's LARs directly.
//! * **eq (30)** — flat/tilted classification with hysteresis:
//!   ```text
//!   flat^(m) = 0            if o₁ < −1.74 and o₂ > 0.65 and flat^(m−1) = 1
//!   flat^(m) = 1            if (o₁ > −1.52 or o₂ < 0.43) and flat^(m−1) = 0
//!   flat^(m) = flat^(m−1)   otherwise
//!   ```
//! * flat ⇒ `γ₁ = 0.94, γ₂ = 0.6`; tilted ⇒ `γ₁ = 0.98` and `γ₂` from:
//! * **eq (31)** — `d_min = min[ω_{i+1} − ω_i]`, `i = 1 … 9` (the
//!   current subframe's LSF adjacency minimum), and
//! * **eq (32)** — `γ₂ = −6.0·d_min + 1.0`, bounded `0.4 ≤ γ₂ ≤ 0.7`.
//! * **eq (33)** — the weighted speech for the subframe:
//!   `sw(n) = s(n) + Σ a_i·γ₁ⁱ·s(n−i) − Σ a_i·γ₂ⁱ·sw(n−i)`,
//!   `n = 0 … 39` — the FIR (numerator) part runs on the pre-processed
//!   speech `s(n)`, the IIR (denominator) part on `sw(n)` itself, with
//!   both 10-tap memories carried across subframes.
//!
//! ## State (clause 4.3)
//!
//! The previous frame's LAR pair (for eq (29)), the previous subframe's
//! `flat` classification (for the eq (30) hysteresis), and the two
//! 10-tap eq (33) filter memories all persist across calls and are
//! zero-initialised (clause 4.3 lists no Table-9 exception for them;
//! `flat` starts at 1, the neutral "flat" state that selects the fixed
//! `γ` pair — the hysteresis then adapts from the first real frame).

use crate::tables::M;

/// Samples per subframe (clause 2.1).
pub const L_SUBFR: usize = 40;

/// γ pair for a **flat** spectral envelope (clause 3.3).
pub const GAMMA_FLAT: (f32, f32) = (0.94, 0.6);
/// γ₁ for a **tilted** spectral envelope (clause 3.3).
pub const GAMMA1_TILTED: f32 = 0.98;
/// Bounds on the tilted-case γ₂ (eq (32)).
pub const GAMMA2_BOUNDS: (f32, f32) = (0.4, 0.7);

/// The per-subframe weighting decision: the chosen `(γ₁, γ₂)` pair and
/// the flat/tilted classification that produced it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WeightingDecision {
    /// eq (27) numerator bandwidth-expansion factor `γ₁`.
    pub gamma1: f32,
    /// eq (27) denominator bandwidth-expansion factor `γ₂`.
    pub gamma2: f32,
    /// eq (30) classification (`true` = flat).
    pub flat: bool,
}

/// eq (28): converts a reflection coefficient to a log-area ratio.
/// `k` is clamped away from ±1 to keep the log finite on degenerate
/// input.
#[must_use]
pub fn lar(k: f32) -> f32 {
    let k = k.clamp(-0.9999, 0.9999);
    ((1.0 + k) / (1.0 - k)).ln()
}

/// eq (31): the minimum adjacent LSF distance `d_min` of one subframe's
/// LSF vector `ω` (radians, ascending).
#[must_use]
pub fn min_lsf_distance(omega: &[f32; M]) -> f32 {
    omega
        .windows(2)
        .map(|w| w[1] - w[0])
        .fold(f32::INFINITY, f32::min)
}

/// §3.3 adaptive perceptual weighting state machine + eq (33) filter.
#[derive(Debug, Clone)]
pub struct PerceptualWeighting {
    /// Previous frame's LAR pair `o₁, o₂` (eq (29) interpolation input).
    prev_lar: [f32; 2],
    /// Previous subframe's eq (30) classification.
    prev_flat: bool,
    /// eq (33) numerator memory — the last `M` pre-processed samples
    /// `s(n−1) … s(n−M)` (most recent first).
    s_mem: [f32; M],
    /// eq (33) denominator memory — the last `M` weighted samples
    /// `sw(n−1) … sw(n−M)` (most recent first).
    sw_mem: [f32; M],
}

impl Default for PerceptualWeighting {
    fn default() -> Self {
        Self::new()
    }
}

impl PerceptualWeighting {
    /// Fresh state: zero LARs / filter memories, `flat = true` start.
    #[must_use]
    pub fn new() -> Self {
        Self {
            prev_lar: [0.0; 2],
            prev_flat: true,
            s_mem: [0.0; M],
            sw_mem: [0.0; M],
        }
    }

    /// Runs the eq (29)/(30)/(31)/(32) γ-adaptation for both subframes
    /// of one frame. `reflection2` is the `(k₁, k₂)` pair from the
    /// Levinson by-product; `omega_sub` holds each subframe's
    /// (unquantised, interpolated) LSF vector for the eq (31) resonance
    /// criterion. Advances the LAR / flat state.
    pub fn adapt_frame(
        &mut self,
        reflection2: (f32, f32),
        omega_sub: &[[f32; M]; 2],
    ) -> [WeightingDecision; 2] {
        let cur_lar = [lar(reflection2.0), lar(reflection2.1)];
        // eq (29): subframe 1 = mean of previous and current LARs;
        // subframe 2 = current LARs.
        let lar_sub1 = [
            0.5 * self.prev_lar[0] + 0.5 * cur_lar[0],
            0.5 * self.prev_lar[1] + 0.5 * cur_lar[1],
        ];
        let lar_sub2 = cur_lar;
        self.prev_lar = cur_lar;

        let mut out = [WeightingDecision {
            gamma1: GAMMA_FLAT.0,
            gamma2: GAMMA_FLAT.1,
            flat: true,
        }; 2];

        for (sub, (lars, omega)) in [lar_sub1, lar_sub2]
            .iter()
            .zip(omega_sub.iter())
            .enumerate()
        {
            // eq (30) hysteresis.
            let flat = if self.prev_flat && lars[0] < -1.74 && lars[1] > 0.65 {
                false
            } else if !self.prev_flat && (lars[0] > -1.52 || lars[1] < 0.43) {
                true
            } else {
                self.prev_flat
            };
            self.prev_flat = flat;

            let (gamma1, gamma2) = if flat {
                GAMMA_FLAT
            } else {
                // eq (31)/(32).
                let d_min = min_lsf_distance(omega);
                let g2 = (-6.0 * d_min + 1.0).clamp(GAMMA2_BOUNDS.0, GAMMA2_BOUNDS.1);
                (GAMMA1_TILTED, g2)
            };
            out[sub] = WeightingDecision {
                gamma1,
                gamma2,
                flat,
            };
        }
        out
    }

    /// eq (33): filters one 40-sample subframe of pre-processed speech
    /// `s` through `W(z) = A(z/γ₁)/A(z/γ₂)` built from the unquantised
    /// LP coefficients `a` (`a[i−1] = a_i`), returning the weighted
    /// speech `sw(n)`. Both 10-tap memories persist across subframes.
    #[must_use]
    pub fn weight_subframe(
        &mut self,
        s: &[f32; L_SUBFR],
        a: &[f32; M],
        decision: WeightingDecision,
    ) -> [f32; L_SUBFR] {
        // Pre-scale the coefficient rows: aγ₁ⁱ and aγ₂ⁱ.
        let mut num = [0.0f32; M];
        let mut den = [0.0f32; M];
        let mut g1 = 1.0f32;
        let mut g2 = 1.0f32;
        for i in 0..M {
            g1 *= decision.gamma1;
            g2 *= decision.gamma2;
            num[i] = a[i] * g1;
            den[i] = a[i] * g2;
        }

        let mut sw = [0.0f32; L_SUBFR];
        for n in 0..L_SUBFR {
            // FIR part: s(n) + Σ a_i γ₁ⁱ s(n−i).
            let mut acc = s[n];
            for i in 0..M {
                let s_prev = if n > i {
                    s[n - 1 - i]
                } else {
                    self.s_mem[i - n]
                };
                acc += num[i] * s_prev;
            }
            // IIR part: − Σ a_i γ₂ⁱ sw(n−i).
            for i in 0..M {
                let sw_prev = if n > i {
                    sw[n - 1 - i]
                } else {
                    self.sw_mem[i - n]
                };
                acc -= den[i] * sw_prev;
            }
            sw[n] = acc;
        }

        // Advance the two memories: most recent first.
        for i in 0..M {
            self.s_mem[i] = s[L_SUBFR - 1 - i];
            self.sw_mem[i] = sw[L_SUBFR - 1 - i];
        }
        sw
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// eq (28): k = 0 → LAR 0; the function is odd; a positive k gives a
    /// positive LAR.
    #[test]
    fn lar_basics() {
        assert_eq!(lar(0.0), 0.0);
        assert!((lar(0.5) + lar(-0.5)).abs() < 1e-6);
        assert!(lar(0.5) > 0.0);
        // Degenerate |k| ≥ 1 stays finite.
        assert!(lar(1.0).is_finite());
        assert!(lar(-1.0).is_finite());
    }

    /// eq (31): minimum adjacency of a uniform grid is the grid step.
    #[test]
    fn min_distance_uniform_grid() {
        let omega: [f32; M] = std::array::from_fn(|i| (i + 1) as f32 * std::f32::consts::PI / 11.0);
        let d = min_lsf_distance(&omega);
        assert!((d - std::f32::consts::PI / 11.0).abs() < 1e-6);
    }

    /// eq (30): from the flat start, strongly tilted LARs (o₁ < −1.74,
    /// o₂ > 0.65) flip to tilted; mildly flat LARs flip back.
    #[test]
    fn flat_hysteresis_transitions() {
        let mut w = PerceptualWeighting::new();
        let omega: [[f32; M]; 2] = {
            let row: [f32; M] =
                std::array::from_fn(|i| (i + 1) as f32 * std::f32::consts::PI / 11.0);
            [row, row]
        };
        // LAR(k) = ln((1+k)/(1−k)): k₁ = −0.75 → o₁ ≈ −1.945 < −1.74;
        // k₂ = 0.4 → o₂ ≈ 0.847 > 0.65 ⇒ tilted. On the FIRST frame the
        // eq (29) subframe-1 interpolation averages with the zero
        // start-up LARs (o₁ ≈ −0.97 > −1.74 → stays flat), so the flip
        // to tilted lands on subframe 2, which uses the full current
        // LARs.
        let d = w.adapt_frame((-0.75, 0.4), &omega);
        assert!(d[0].flat, "subframe 1 still flat (interpolated LARs)");
        assert!(!d[1].flat, "subframe 2 should classify tilted");
        assert!((d[1].gamma1 - GAMMA1_TILTED).abs() < 1e-6);
        assert!(d[1].gamma2 >= GAMMA2_BOUNDS.0 && d[1].gamma2 <= GAMMA2_BOUNDS.1);

        // Now feed clearly flat LARs (o₁ > −1.52 via k₁ = 0): the
        // subframe-1 interpolated o₁ = −0.97 > −1.52 already releases
        // the hysteresis back to flat.
        let d2 = w.adapt_frame((0.0, 0.4), &omega);
        assert!(d2[0].flat, "hysteresis should release to flat");
        assert!((d2[0].gamma1 - GAMMA_FLAT.0).abs() < 1e-6);
        assert!((d2[0].gamma2 - GAMMA_FLAT.1).abs() < 1e-6);
    }

    /// eq (32): a tightly clustered LSF pair (strong resonance) pushes
    /// γ₂ to the 0.7 upper bound; a wide-open spectrum pins the 0.4
    /// lower bound.
    #[test]
    fn gamma2_resonance_mapping() {
        let mut w = PerceptualWeighting::new();
        // Force tilted classification.
        let tight: [f32; M] = {
            let mut v: [f32; M] = std::array::from_fn(|i| 0.3 + 0.28 * i as f32);
            v[5] = v[4] + 0.01; // d_min = 0.01 → γ₂ = 0.94 → clamp 0.7
            v
        };
        // On the first frame the tilt flip lands on subframe 2 (the
        // subframe-1 LARs are interpolated with the zero start-up pair).
        let d = w.adapt_frame((-0.75, 0.4), &[tight, tight]);
        assert!(!d[1].flat);
        assert!((d[1].gamma2 - 0.7).abs() < 1e-6, "γ₂={}", d[1].gamma2);

        let mut w2 = PerceptualWeighting::new();
        let wide: [f32; M] = std::array::from_fn(|i| 0.2 + 0.30 * i as f32);
        // d_min = 0.30 → γ₂ = 1 − 1.8 = −0.8 → clamp 0.4.
        let d2 = w2.adapt_frame((-0.75, 0.4), &[wide, wide]);
        assert!(!d2[1].flat);
        assert!((d2[1].gamma2 - 0.4).abs() < 1e-6, "γ₂={}", d2[1].gamma2);
    }

    /// eq (33) with γ₁ = γ₂ reduces W(z) to the identity — the weighted
    /// speech equals the input (numerator and denominator cancel).
    #[test]
    fn equal_gammas_identity() {
        let mut w = PerceptualWeighting::new();
        let a: [f32; M] = [-0.5, 0.2, -0.1, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let dec = WeightingDecision {
            gamma1: 0.8,
            gamma2: 0.8,
            flat: true,
        };
        let s: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 13 % 100) as f32) - 50.0);
        let sw = w.weight_subframe(&s, &a, dec);
        for (n, (&x, &y)) in s.iter().zip(sw.iter()).enumerate() {
            assert!(
                (x - y).abs() < 1e-3,
                "identity violated at n={n}: {x} vs {y}"
            );
        }
    }

    /// eq (33) memory continuity: filtering 80 samples as two subframes
    /// equals filtering them with a single fresh filter run of the same
    /// recursion (state threading check).
    #[test]
    fn subframe_state_continuity() {
        let a: [f32; M] = [-0.6, 0.3, -0.15, 0.07, -0.03, 0.01, 0.0, 0.0, 0.0, 0.0];
        let dec = WeightingDecision {
            gamma1: 0.94,
            gamma2: 0.6,
            flat: true,
        };
        let frame: [f32; 80] = std::array::from_fn(|n| ((n * 29 % 211) as f32) - 105.0);

        // Two-subframe path.
        let mut w = PerceptualWeighting::new();
        let s1: [f32; L_SUBFR] = std::array::from_fn(|n| frame[n]);
        let s2: [f32; L_SUBFR] = std::array::from_fn(|n| frame[L_SUBFR + n]);
        let sw1 = w.weight_subframe(&s1, &a, dec);
        let sw2 = w.weight_subframe(&s2, &a, dec);

        // Whole-frame direct recursion.
        let mut num = [0.0f32; M];
        let mut den = [0.0f32; M];
        let (mut g1, mut g2) = (1.0f32, 1.0f32);
        for i in 0..M {
            g1 *= dec.gamma1;
            g2 *= dec.gamma2;
            num[i] = a[i] * g1;
            den[i] = a[i] * g2;
        }
        let mut sw_ref = [0.0f32; 80];
        for n in 0..80 {
            let mut acc = frame[n];
            for i in 0..M {
                if n > i {
                    acc += num[i] * frame[n - 1 - i];
                    acc -= den[i] * sw_ref[n - 1 - i];
                }
            }
            sw_ref[n] = acc;
        }
        for n in 0..L_SUBFR {
            assert!(
                (sw1[n] - sw_ref[n]).abs() < 1e-3,
                "subframe 1 mismatch at {n}"
            );
            assert!(
                (sw2[n] - sw_ref[L_SUBFR + n]).abs() < 1e-3,
                "subframe 2 mismatch at {n}"
            );
        }
    }
}
