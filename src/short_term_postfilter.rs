//! §4.2.2 **short-term postfilter** `H_f(z)` — the spectral-shaping
//! stage of the decoder's §4.2 post-processing cascade.
//!
//! The §4.2 cascade is (clause 4.2): long-term postfilter `H_p(z)`
//! (§4.2.1) → **short-term postfilter `H_f(z)` (§4.2.2)** → tilt
//! compensation `H_t(z)` (§4.2.3) → adaptive gain control (§4.2.4) →
//! output high-pass + ×2 upscaling (§4.2.5). The round-298
//! [`crate::post_process`] module wired the §4.2.5 tail; this module
//! wires the §4.2.2 short-term postfilter, which is the stage that is
//! both fully self-contained (it depends only on the per-subframe
//! quantised LP coefficients `â_i`, the two constant weight factors,
//! and its own filter memories) and the structural backbone of the
//! cascade.
//!
//! ## Spec source — clause 4.2.2, equations (84)/(85) (06/2012 Rec.)
//!
//! Clause 4.2.2 (transcribed from the EPUB prose): "The short-term
//! postfilter is given by [eq (84)] where `Â(z)` is the received
//! quantized LP inverse filter (LP analysis is not done at the decoder)
//! and the factors `γ_n` and `γ_d` control the amount of short-term
//! postfiltering, and are set to `γ_n = 0.55`, and `γ_d = 0.7`. The
//! gain term `g_f` is calculated on the truncated impulse response
//! `h_f(n)` of the filter `Â(z/γ_n)/Â(z/γ_d)` and is given by
//! [eq (85)]."
//!
//! Equation (84) (transcribed from the EPUB's equation raster
//! `images/eq84.jpg`) is the transfer function
//!
//! ```text
//!            1   Â(z/γ_n)    1   1 + Σ_{i=1}^{10} γ_n^i·â_i·z^−i
//! H_f(z) = ─── · ──────── = ─── · ────────────────────────────────
//!           g_f  Â(z/γ_d)   g_f  1 + Σ_{i=1}^{10} γ_d^i·â_i·z^−i
//! ```
//!
//! and equation (85) (raster `images/eq85.jpg`) is the gain-
//! normalisation term
//!
//! ```text
//!         19
//! g_f =   Σ  | h_f(n) |
//!        n=0
//! ```
//!
//! where `h_f(n)` is the impulse response of the *un-normalised* filter
//! `Â(z/γ_n)/Â(z/γ_d)` (i.e. eq (84) without the leading `1/g_f`),
//! truncated to its first 20 samples (`n = 0 … 19`).
//!
//! The filter is therefore applied per subframe as three steps:
//!
//! 1. **Numerator (all-zero) `Â(z/γ_n)`** — the residual signal
//!    `r(n) = x(n) + Σ_{i=1}^{10} γ_n^i·â_i·x(n−i)`. (Clause 4.2.1
//!    names this same residual: it "is obtained by filtering the speech
//!    `ŝ(n)` through `Â(z/γ_n)`, which is the numerator of the
//!    short-term postfilter".)
//! 2. **Denominator (all-pole) `1/Â(z/γ_d)`** —
//!    `y(n) = r(n) − Σ_{i=1}^{10} γ_d^i·â_i·y(n−i)`.
//! 3. **Gain normalisation** — scale by `1/g_f`, with `g_f` recomputed
//!    each subframe from the current `â_i` per eq (85).
//!
//! ## Constants (clause 4.2.2 / Table 5)
//!
//! Table 5 ("Glossary of most relevant constants") lists
//! `γ_n = 0.55` and `γ_d = 0.70` as the postfilter weight factors;
//! clause 4.2.2 repeats both values inline.
//!
//! ## State (clause 4.3 init)
//!
//! Per clause 4.3, "all static encoder and decoder variables should be
//! initialized to zero, except the variables listed in Table 9". Neither
//! the numerator input history `x(n−i)` nor the denominator output
//! history `y(n−i)` appears in Table 9, so both 10-tap memories start
//! zeroed and carry across subframes (the postfilter coefficients are
//! "updated every 5 ms subframe", clause 4.2, but the signal memory is
//! continuous).
//!
//! ## What this module does NOT do
//!
//! It is the §4.2.2 stage only. The §4.2.1 long-term (harmonic)
//! postfilter that precedes it (its residual-domain two-pass pitch
//! search) is a separate follow-up round. The §4.2.3 tilt compensation
//! that follows it consumes this stage's [`Self::impulse_response`]
//! (its first reflection coefficient `k1'`, eq (87)); see
//! [`crate::tilt_compensation`]. Until §4.2.1 lands, the input to this
//! stage is the raw §4.1.6 reconstructed speech `ŝ(n)` from
//! [`crate::lp_synthesis`] rather than the long-term-postfiltered
//! residual the full cascade would feed in.

use crate::fixed_codebook::SUBFRAME_SIZE;
use crate::tables::M;

/// §4.2.2 / Table 5 numerator weight factor `γ_n = 0.55`.
pub const GAMMA_N: f32 = 0.55;

/// §4.2.2 / Table 5 denominator weight factor `γ_d = 0.70`.
pub const GAMMA_D: f32 = 0.70;

/// eq (85) impulse-response truncation length — the gain term sums the
/// magnitudes of the first 20 samples `h_f(0 … 19)`.
pub const GF_IMPULSE_LEN: usize = 20;

/// Stateful §4.2.2 short-term postfilter `H_f(z)` (eqs (84)/(85)).
///
/// Owns the two 10-tap filter memories (numerator input history and
/// denominator output history), both zero-initialised per clause 4.3 and
/// carried across subframes. The per-subframe quantised LP coefficients
/// `â_i` are supplied to [`Self::filter_subframe`] each call (the
/// coefficients change every subframe; the memory is continuous).
#[derive(Debug, Clone)]
pub struct ShortTermPostfilter {
    /// Numerator `Â(z/γ_n)` input history `[x(n−1) … x(n−10)]`
    /// (`x_hist[0]` = `x(n−1)`). Zero-init per clause 4.3.
    x_hist: [f32; M],
    /// Denominator `1/Â(z/γ_d)` output history `[y(n−1) … y(n−10)]`
    /// (`y_hist[0]` = `y(n−1)`). Zero-init per clause 4.3.
    y_hist: [f32; M],
}

impl Default for ShortTermPostfilter {
    fn default() -> Self {
        Self::new()
    }
}

impl ShortTermPostfilter {
    /// Build the postfilter with the clause-4.3 all-zero start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            x_hist: [0.0; M],
            y_hist: [0.0; M],
        }
    }

    /// Borrow the numerator input history `[x(n−1) … x(n−10)]`.
    #[must_use]
    pub fn x_hist(&self) -> &[f32; M] {
        &self.x_hist
    }

    /// Borrow the denominator output history `[y(n−1) … y(n−10)]`.
    #[must_use]
    pub fn y_hist(&self) -> &[f32; M] {
        &self.y_hist
    }

    /// The weighted numerator coefficients `γ_n^i·â_i` for `i = 1 … 10`
    /// (slot `j` holds the `i = j+1` coefficient), from the subframe
    /// `â_i` (slots `0 … 9` hold `â_1 … â_10`).
    ///
    /// Exposed crate-internally because the §4.2.1 long-term postfilter's
    /// eq (79) residual `r̂(n) = ŝ(n) + Σ γ_n^i·â_i·ŝ(n−i)` is, per
    /// clause 4.2.1, "the numerator of the short-term postfilter" applied
    /// to `ŝ(n)` — the very `γ_n^i·â_i` weights computed here. Sharing the
    /// helper keeps the two stages on one definition (see
    /// [`crate::long_term_postfilter`]).
    pub(crate) fn weighted_num(a: &[f32; M]) -> [f32; M] {
        let mut w = [0.0f32; M];
        let mut g = GAMMA_N;
        for i in 0..M {
            w[i] = g * a[i];
            g *= GAMMA_N;
        }
        w
    }

    /// The weighted denominator coefficients `γ_d^i·â_i` for
    /// `i = 1 … 10`.
    fn weighted_den(a: &[f32; M]) -> [f32; M] {
        let mut w = [0.0f32; M];
        let mut g = GAMMA_D;
        for i in 0..M {
            w[i] = g * a[i];
            g *= GAMMA_D;
        }
        w
    }

    /// The first [`GF_IMPULSE_LEN`] samples of the impulse response
    /// `h_f(n)` of the un-normalised filter `Â(z/γ_n)/Â(z/γ_d)`
    /// (numerator over denominator, no `1/g_f`).
    ///
    /// Computed by driving a *local, zero-state* copy of the
    /// numerator-then-denominator filter with a unit impulse for
    /// [`GF_IMPULSE_LEN`] samples. This never touches the carried-over
    /// [`Self::x_hist`] / [`Self::y_hist`] state — `h_f(n)` is a property
    /// of the coefficients alone. It is the same response whose absolute
    /// values eq (85) sums for `g_f`, and whose autocorrelation eq (87)
    /// uses for the §4.2.3 tilt factor `k1'`.
    #[must_use]
    pub fn impulse_response(a: &[f32; M]) -> [f32; GF_IMPULSE_LEN] {
        let wn = Self::weighted_num(a);
        let wd = Self::weighted_den(a);
        // Local zero-state impulse response of Â(z/γ_n)/Â(z/γ_d).
        let mut xh = [0.0f32; M]; // numerator input history
        let mut yh = [0.0f32; M]; // denominator output history
        let mut h = [0.0f32; GF_IMPULSE_LEN];
        for (n, hn) in h.iter_mut().enumerate() {
            let x = if n == 0 { 1.0 } else { 0.0 };
            // Numerator: r = x + Σ wn[i]·x(n−1−i).
            let mut r = x;
            for (&w, &xv) in wn.iter().zip(xh.iter()) {
                r += w * xv;
            }
            // Denominator: y = r − Σ wd[i]·y(n−1−i).
            let mut y = r;
            for (&w, &yv) in wd.iter().zip(yh.iter()) {
                y -= w * yv;
            }
            *hn = y;
            // Advance local histories (most-recent at index 0).
            for i in (1..M).rev() {
                xh[i] = xh[i - 1];
                yh[i] = yh[i - 1];
            }
            xh[0] = x;
            yh[0] = y;
        }
        h
    }

    /// eq (85) gain term `g_f = Σ_{n=0}^{19} |h_f(n)|`, where `h_f(n)`
    /// is the impulse response of the un-normalised filter
    /// `Â(z/γ_n)/Â(z/γ_d)` (numerator over denominator, no `1/g_f`).
    ///
    /// `g_f` is a property of the coefficients alone — it never touches
    /// the carried-over [`Self::x_hist`] / [`Self::y_hist`] state.
    #[must_use]
    pub fn gain_term(a: &[f32; M]) -> f32 {
        Self::impulse_response(a).iter().map(|v| v.abs()).sum()
    }

    /// Filter one 40-sample subframe `x(n)` through `H_f(z)` (eq (84))
    /// with the per-subframe quantised LP coefficients `a` (slots
    /// `0 … 9` hold `â_1 … â_10`), advancing the carried-over state.
    ///
    /// Returns the gain-normalised postfiltered subframe `H_f(z)·x`.
    /// `g_f` (eq (85)) is recomputed from `a` each call; if it is
    /// non-positive (cannot happen for a stable `Â`, but guarded for a
    /// hand-built degenerate `a`) the gain normalisation is skipped
    /// (scale 1.0) rather than dividing by zero.
    #[must_use]
    pub fn filter_subframe(
        &mut self,
        x: &[f32; SUBFRAME_SIZE],
        a: &[f32; M],
    ) -> [f32; SUBFRAME_SIZE] {
        let wn = Self::weighted_num(a);
        let wd = Self::weighted_den(a);
        let g_f = Self::gain_term(a);
        let inv_gf = if g_f > 0.0 { 1.0 / g_f } else { 1.0 };

        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &xn) in x.iter().enumerate() {
            // Numerator Â(z/γ_n): r(n) = x(n) + Σ wn[i]·x(n−1−i).
            let mut r = xn;
            for (&w, &xh) in wn.iter().zip(self.x_hist.iter()) {
                r += w * xh;
            }
            // Denominator 1/Â(z/γ_d): y(n) = r(n) − Σ wd[i]·y(n−1−i).
            let mut y = r;
            for (&w, &yh) in wd.iter().zip(self.y_hist.iter()) {
                y -= w * yh;
            }
            // Advance state (most-recent at index 0).
            for i in (1..M).rev() {
                self.x_hist[i] = self.x_hist[i - 1];
                self.y_hist[i] = self.y_hist[i - 1];
            }
            self.x_hist[0] = xn;
            self.y_hist[0] = y;
            // Gain normalisation 1/g_f (eq (84) leading factor).
            out[n] = inv_gf * y;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A flat (all-zero `â_i`) filter is the identity numerator and
    /// denominator, so its impulse response is the unit impulse and
    /// `g_f = 1`; the postfilter then passes the input through unchanged.
    #[test]
    fn flat_filter_is_identity() {
        let a = [0.0f32; M];
        assert!((ShortTermPostfilter::gain_term(&a) - 1.0).abs() < 1e-6);

        let mut pf = ShortTermPostfilter::new();
        let mut x = [0.0f32; SUBFRAME_SIZE];
        for (n, v) in x.iter_mut().enumerate() {
            *v = (n as f32) - 20.0;
        }
        let out = pf.filter_subframe(&x, &a);
        for n in 0..SUBFRAME_SIZE {
            assert!((out[n] - x[n]).abs() < 1e-6, "n={n}");
        }
    }

    /// New filter starts zeroed (clause 4.3).
    #[test]
    fn new_filter_has_zero_state() {
        let pf = ShortTermPostfilter::new();
        assert_eq!(pf.x_hist(), &[0.0; M]);
        assert_eq!(pf.y_hist(), &[0.0; M]);
    }

    /// The weight factors are the spec values (Table 5 / clause 4.2.2).
    #[test]
    fn weight_factors_match_spec() {
        assert!((GAMMA_N - 0.55).abs() < 1e-9);
        assert!((GAMMA_D - 0.70).abs() < 1e-9);
        assert_eq!(GF_IMPULSE_LEN, 20);
    }

    /// Weighted coefficients are `γ^i·â_i` with `i = 1 … 10`. With a
    /// single non-zero `â_1 = 1`, the numerator weight is `γ_n^1` and
    /// the denominator weight is `γ_d^1`; with `â_2 = 1` they are
    /// `γ_n^2` / `γ_d^2`.
    #[test]
    fn weighting_is_gamma_power_i() {
        let mut a = [0.0f32; M];
        a[0] = 1.0; // â_1
        let wn = ShortTermPostfilter::weighted_num(&a);
        let wd = ShortTermPostfilter::weighted_den(&a);
        assert!((wn[0] - GAMMA_N).abs() < 1e-9);
        assert!((wd[0] - GAMMA_D).abs() < 1e-9);

        let mut a = [0.0f32; M];
        a[1] = 1.0; // â_2
        let wn = ShortTermPostfilter::weighted_num(&a);
        let wd = ShortTermPostfilter::weighted_den(&a);
        assert!((wn[1] - GAMMA_N * GAMMA_N).abs() < 1e-9);
        assert!((wd[1] - GAMMA_D * GAMMA_D).abs() < 1e-9);
    }

    /// `g_f` (eq (85)) reproduced by an independent direct impulse
    /// recursion for a single-tap `â_1`. The un-normalised filter is
    /// `(1 + γ_n·z^−1)/(1 + γ_d·z^−1)`; its impulse response is
    /// `h(0) = 1`, `h(n) = (γ_n − γ_d)·(−γ_d)^{n−1}` for `n ≥ 1`.
    #[test]
    fn gain_term_matches_hand_recursion_single_tap() {
        let mut a = [0.0f32; M];
        a[0] = 1.0;
        let gf = ShortTermPostfilter::gain_term(&a);

        let mut expected = 1.0f32; // |h(0)|
        for n in 1..GF_IMPULSE_LEN {
            let h = (GAMMA_N - GAMMA_D) * (-GAMMA_D).powi((n - 1) as i32);
            expected += h.abs();
        }
        assert!((gf - expected).abs() < 1e-5, "got {gf}, want {expected}");
    }

    /// `gain_term` is the sum of `|h_f(n)|` over the impulse response, so
    /// the two helpers stay consistent: `g_f = Σ |impulse_response|`.
    #[test]
    fn gain_term_is_sum_of_impulse_magnitudes() {
        let mut a = [0.0f32; M];
        a[0] = 0.7;
        a[1] = -0.25;
        a[5] = 0.1;
        let h = ShortTermPostfilter::impulse_response(&a);
        let sum: f32 = h.iter().map(|v| v.abs()).sum();
        assert!((sum - ShortTermPostfilter::gain_term(&a)).abs() < 1e-6);
        // h_f(0) of Â(z/γ_n)/Â(z/γ_d) is the unit impulse passed through
        // both an FIR with leading 1 and an all-pole with leading 1 → 1.
        assert!((h[0] - 1.0).abs() < 1e-6);
    }

    /// `h_f(1)` of the single-tap filter `(1 + γ_n·z⁻¹)/(1 + γ_d·z⁻¹)` is
    /// `γ_n − γ_d` (the order-1 step of the recursion in the gain-term
    /// hand check above), pinning the impulse-response sign convention.
    #[test]
    fn impulse_response_single_tap_first_sample() {
        let mut a = [0.0f32; M];
        a[0] = 1.0;
        let h = ShortTermPostfilter::impulse_response(&a);
        assert!((h[0] - 1.0).abs() < 1e-6);
        assert!((h[1] - (GAMMA_N - GAMMA_D)).abs() < 1e-6);
    }

    /// `g_f` depends only on the coefficients, not on the carried state:
    /// the static [`ShortTermPostfilter::gain_term`] equals the value
    /// the filter would use after the state has been driven by prior
    /// data.
    #[test]
    fn gain_term_is_state_independent() {
        let mut a = [0.0f32; M];
        a[0] = 0.6;
        a[2] = -0.2;
        let g_clean = ShortTermPostfilter::gain_term(&a);

        let mut pf = ShortTermPostfilter::new();
        // Drive arbitrary data through to populate the state.
        let x: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32).sin() * 500.0);
        let _ = pf.filter_subframe(&x, &a);
        // The gain term recomputed for the same coefficients is identical.
        assert_eq!(g_clean, ShortTermPostfilter::gain_term(&a));
    }

    /// State carries across subframes: filtering an 80-sample stream as
    /// two 40-sample subframes (same `â_i`) equals filtering the
    /// concatenation through a single continuous recursion built by
    /// hand.
    #[test]
    fn state_carries_across_subframes() {
        let mut a = [0.0f32; M];
        a[0] = 0.5;
        a[1] = -0.3;
        a[4] = 0.1;

        let mut x1 = [0.0f32; SUBFRAME_SIZE];
        let mut x2 = [0.0f32; SUBFRAME_SIZE];
        for n in 0..SUBFRAME_SIZE {
            x1[n] = ((n * 3) % 11) as f32 - 5.0;
            x2[n] = ((n * 7) % 13) as f32 - 6.0;
        }

        let mut pf = ShortTermPostfilter::new();
        let o1 = pf.filter_subframe(&x1, &a);
        let o2 = pf.filter_subframe(&x2, &a);

        // Independent continuous reference over the 80-sample stream.
        let wn = ShortTermPostfilter::weighted_num(&a);
        let wd = ShortTermPostfilter::weighted_den(&a);
        let inv_gf = 1.0 / ShortTermPostfilter::gain_term(&a);
        let mut stream = [0.0f32; 2 * SUBFRAME_SIZE];
        stream[..SUBFRAME_SIZE].copy_from_slice(&x1);
        stream[SUBFRAME_SIZE..].copy_from_slice(&x2);
        let mut xh = [0.0f32; M];
        let mut yh = [0.0f32; M];
        let mut ref_out = [0.0f32; 2 * SUBFRAME_SIZE];
        for n in 0..2 * SUBFRAME_SIZE {
            let mut r = stream[n];
            for i in 0..M {
                r += wn[i] * xh[i];
            }
            let mut y = r;
            for i in 0..M {
                y -= wd[i] * yh[i];
            }
            for i in (1..M).rev() {
                xh[i] = xh[i - 1];
                yh[i] = yh[i - 1];
            }
            xh[0] = stream[n];
            yh[0] = y;
            ref_out[n] = inv_gf * y;
        }

        for n in 0..SUBFRAME_SIZE {
            assert!((o1[n] - ref_out[n]).abs() < 1e-4, "sub1 n={n}");
            assert!(
                (o2[n] - ref_out[SUBFRAME_SIZE + n]).abs() < 1e-4,
                "sub2 n={n}"
            );
        }
    }

    /// Two filters over the same subframe sequence stay in lockstep —
    /// all state is owned, no hidden globals.
    #[test]
    fn deterministic() {
        let mut a = [0.0f32; M];
        a[0] = 0.4;
        a[3] = -0.15;
        let x: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) * 2.0 - 39.0);

        let mut p = ShortTermPostfilter::new();
        let mut q = ShortTermPostfilter::new();
        for _ in 0..4 {
            assert_eq!(p.filter_subframe(&x, &a), q.filter_subframe(&x, &a));
        }
    }

    /// On a realistic stable LP coefficient set the postfilter output
    /// stays finite and the gain term is positive and finite.
    #[test]
    fn finite_on_realistic_coefficients() {
        // A mildly-resonant 10th-order set (values in a plausible
        // decoder range; magnitudes < 1).
        let a: [f32; M] = [
            0.9, -0.4, 0.2, -0.1, 0.05, -0.03, 0.02, -0.01, 0.005, -0.002,
        ];
        let g_f = ShortTermPostfilter::gain_term(&a);
        assert!(g_f.is_finite() && g_f > 0.0, "g_f = {g_f}");

        let mut pf = ShortTermPostfilter::new();
        let x: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32 - 20.0) * 100.0);
        for _ in 0..10 {
            let out = pf.filter_subframe(&x, &a);
            assert!(out.iter().all(|v| v.is_finite()));
        }
    }
}
