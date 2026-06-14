//! §4.2.5 **output high-pass filtering and upscaling** — the final
//! stage of the decoder's §4.2 post-processing cascade.
//!
//! The §4.2 cascade is (clause 4.2): long-term postfilter `H_p(z)`
//! (4.2.1) → short-term postfilter `H_f(z)` (4.2.2) → tilt compensation
//! `H_t(z)` (4.2.3) → adaptive gain control (4.2.4) → **output
//! high-pass + ×2 upscaling (4.2.5)**. The four front stages are
//! follow-up rounds; this module wires the cascade *tail*, which is the
//! one stage that is both fully self-contained (a fixed 2nd-order IIR
//! plus a constant scale) and entirely table-backed by an already-
//! compiled coefficient set.
//!
//! ## Spec source — clause 4.2.5, equation (91) (06/2012 Recommendation)
//!
//! Clause 4.2.5 (transcribed from the EPUB prose): "A high-pass filter
//! with a cut-off frequency of 100 Hz is applied to the reconstructed
//! postfiltered speech `sf′(n)`. The filter is given by [eq (91)]. The
//! filtered signal is multiplied by a factor 2 to restore the input
//! signal level."
//!
//! Equation (91) (transcribed from the EPUB's equation raster
//! `images/eq91.jpg`) is the transfer function
//!
//! ```text
//!            0.93980581 − 1.8795834·z⁻¹ + 0.93980581·z⁻²
//! H_h2(z) = ─────────────────────────────────────────────
//!            1 − 1.9330735·z⁻¹ + 0.93589199·z⁻²
//! ```
//!
//! i.e. the difference equation
//!
//! ```text
//! y(n) = b0·x(n) + b1·x(n−1) + b2·x(n−2)
//!              + 1.9330735·y(n−1) − 0.93589199·y(n−2)
//! ```
//!
//! The 16-bit fixed-point coefficients (clause 2.5: the prose decimals
//! are "rounded versions of the values used in the 16-bit fixed-point
//! implementation") are the already-compiled Q13 tables
//! [`crate::tables::HPF_PREPROC_100HZ_B_Q13`] (`{7699, −15398, 7699}`)
//! and [`crate::tables::HPF_PREPROC_100HZ_A_Q13`] (`{8192, 15836,
//! −7667}`). Converting to real values:
//!
//! * `b = {7699, −15398, 7699} / 2^13 = {0.939819, −1.879639,
//!   0.939819}` — matches the eq (91) numerator.
//! * `a = {8192, 15836, −7667} / 2^13 = {1.0, 1.933105, −0.935791}`.
//!   The first entry is the (unit) `a0`; the remaining two are the
//!   **feedback gains already arranged for an additive recursion**:
//!   `a[1] = +1.933105` is the `y(n−1)` coefficient and `a[2] =
//!   −0.935791` is the `y(n−2)` coefficient, so the recursion *adds*
//!   `a[1]·y(n−1) + a[2]·y(n−2)`. This reproduces the eq (91)
//!   denominator `1 − 1.9330735·z⁻¹ + 0.93589199·z⁻²` (the
//!   denominator's leading `−` on `z⁻¹` flips to a `+` feedback gain,
//!   and the `+` on `z⁻²` flips to a `−` feedback gain), to within the
//!   Q13 rounding step.
//!
//! ## State (clause 4.3 init)
//!
//! Per clause 4.3, "all static encoder and decoder variables should be
//! initialized to zero, except the variables listed in Table 9". This
//! filter's memory does not appear in Table 9, so both the two input
//! taps `x(n−1)`, `x(n−2)` and the two output taps `y(n−1)`, `y(n−2)`
//! start zeroed.
//!
//! ## What this module does NOT do
//!
//! The four front stages of the §4.2 cascade (long-/short-term
//! postfilter, tilt compensation, adaptive gain control) are not yet
//! wired, so the input to this stage is, for now, the raw §4.1.6
//! reconstructed speech `ŝ(n)` from [`crate::lp_synthesis`] rather than
//! the postfiltered `sf′(n)` the spec names. The eq (91) filter + ×2
//! scale is identical either way; once the front stages land they slot
//! in front of this module unchanged.

use crate::tables::{HPF_PREPROC_100HZ_A_Q13, HPF_PREPROC_100HZ_B_Q13};

/// Q13 fixed-point scale of the eq (91) coefficient tables (`2^13`).
const Q13_SCALE: f32 = 8_192.0;

/// The ×2 output upscaling factor (clause 4.2.5: "multiplied by a
/// factor 2 to restore the input signal level"). The §3.1 pre-processor
/// scales the input down by 2 as an overflow precaution; this restores
/// it.
pub const OUTPUT_UPSCALE: f32 = 2.0;

/// Stateful §4.2.5 output high-pass filter (eq (91)) with the ×2
/// upscaling folded into each output sample.
///
/// The two input taps (`x(n−1)`, `x(n−2)`) and two output taps
/// (`y(n−1)`, `y(n−2)`) are carried across calls so a stream can be fed
/// in arbitrary chunks (per sample, per subframe, or per frame) with no
/// boundary discontinuity. The stored output history is the **pre-×2**
/// filtered value `y(n)` so the recursion stays numerically identical to
/// the spec filter regardless of how the ×2 is applied at the boundary;
/// the ×2 is applied only to the returned samples.
#[derive(Debug, Clone)]
pub struct OutputHighPass {
    /// b-coefficients `{b0, b1, b2}` (eq (91) numerator), real-valued.
    b: [f32; 3],
    /// Feedback gains `{_, a1, a2}` (eq (91) denominator, sign-arranged
    /// for an additive recursion), real-valued. Slot 0 is the unit `a0`
    /// and is unused in the recursion.
    a: [f32; 3],
    /// Input history `[x(n−1), x(n−2)]`. Zero-init per clause 4.3.
    x_hist: [f32; 2],
    /// Pre-×2 output history `[y(n−1), y(n−2)]`. Zero-init per clause 4.3.
    y_hist: [f32; 2],
}

impl Default for OutputHighPass {
    fn default() -> Self {
        Self::new()
    }
}

impl OutputHighPass {
    /// Build the filter with the clause-4.3 all-zero start-up state and
    /// the eq (91) coefficients taken from the compiled Q13 tables.
    #[must_use]
    pub fn new() -> Self {
        let b = [
            f32::from(HPF_PREPROC_100HZ_B_Q13[0]) / Q13_SCALE,
            f32::from(HPF_PREPROC_100HZ_B_Q13[1]) / Q13_SCALE,
            f32::from(HPF_PREPROC_100HZ_B_Q13[2]) / Q13_SCALE,
        ];
        let a = [
            f32::from(HPF_PREPROC_100HZ_A_Q13[0]) / Q13_SCALE,
            f32::from(HPF_PREPROC_100HZ_A_Q13[1]) / Q13_SCALE,
            f32::from(HPF_PREPROC_100HZ_A_Q13[2]) / Q13_SCALE,
        ];
        Self {
            b,
            a,
            x_hist: [0.0; 2],
            y_hist: [0.0; 2],
        }
    }

    /// Filter one input sample `x(n)` through eq (91), advance the
    /// state, and return the ×2-upscaled output `2·y(n)`.
    ///
    /// The recursion is
    /// `y(n) = b0·x(n) + b1·x(n−1) + b2·x(n−2) + a1·y(n−1) + a2·y(n−2)`
    /// with the sign-arranged feedback gains from the Q13 table (see the
    /// module docs).
    #[must_use]
    pub fn filter_sample(&mut self, x: f32) -> f32 {
        let y = self.b[0] * x
            + self.b[1] * self.x_hist[0]
            + self.b[2] * self.x_hist[1]
            + self.a[1] * self.y_hist[0]
            + self.a[2] * self.y_hist[1];
        // Advance taps: x(n−2) ← x(n−1) ← x(n); same for y.
        self.x_hist[1] = self.x_hist[0];
        self.x_hist[0] = x;
        self.y_hist[1] = self.y_hist[0];
        self.y_hist[0] = y;
        OUTPUT_UPSCALE * y
    }

    /// Filter a slice of input samples in place, replacing each with its
    /// ×2-upscaled high-pass output. State carries across the whole
    /// slice and into subsequent calls.
    pub fn filter_in_place(&mut self, samples: &mut [f32]) {
        for s in samples.iter_mut() {
            *s = self.filter_sample(*s);
        }
    }

    /// Filter a slice into a freshly-allocated output vector, leaving the
    /// input untouched. State carries into subsequent calls.
    #[must_use]
    pub fn filter(&mut self, samples: &[f32]) -> Vec<f32> {
        samples.iter().map(|&x| self.filter_sample(x)).collect()
    }

    /// Borrow the real-valued b-coefficients `{b0, b1, b2}` (eq (91)
    /// numerator) for inspection / tests.
    #[must_use]
    pub fn b_coeffs(&self) -> &[f32; 3] {
        &self.b
    }

    /// Borrow the real-valued feedback gains `{a0, a1, a2}` (eq (91)
    /// denominator, sign-arranged for an additive recursion) for
    /// inspection / tests. Slot 0 is the unit `a0`.
    #[must_use]
    pub fn a_coeffs(&self) -> &[f32; 3] {
        &self.a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// New filter starts zeroed (clause 4.3) and reads the eq (91)
    /// coefficients from the compiled Q13 tables.
    #[test]
    fn new_filter_has_zero_state_and_table_coeffs() {
        let f = OutputHighPass::new();
        assert_eq!(f.x_hist, [0.0; 2]);
        assert_eq!(f.y_hist, [0.0; 2]);
        // b = {7699, -15398, 7699} / 8192.
        assert!((f.b_coeffs()[0] - 7699.0 / 8192.0).abs() < 1e-9);
        assert!((f.b_coeffs()[1] - (-15398.0) / 8192.0).abs() < 1e-9);
        assert!((f.b_coeffs()[2] - 7699.0 / 8192.0).abs() < 1e-9);
        // a = {8192, 15836, -7667} / 8192.
        assert!((f.a_coeffs()[0] - 1.0).abs() < 1e-9);
        assert!((f.a_coeffs()[1] - 15836.0 / 8192.0).abs() < 1e-9);
        assert!((f.a_coeffs()[2] - (-7667.0) / 8192.0).abs() < 1e-9);
    }

    /// The real coefficients match the eq (91) decimals (clause 2.5: the
    /// prose decimals are rounded versions of the fixed-point values, so
    /// the Q13 round-trip must agree to within one Q13 step ≈ 1.2e-4).
    #[test]
    fn coeffs_match_eq91_decimals() {
        let f = OutputHighPass::new();
        let q13 = 1.0 / 8192.0;
        // eq (91) numerator (decimals truncated to f32 precision; the
        // q13 tolerance ≈ 1.2e-4 dwarfs the dropped digits). eq (91):
        // 0.93980581, −1.8795834, 0.93980581.
        assert!((f.b_coeffs()[0] - 0.939_805_8).abs() < q13);
        assert!((f.b_coeffs()[1] - (-1.879_583_4)).abs() < q13);
        assert!((f.b_coeffs()[2] - 0.939_805_8).abs() < q13);
        // eq (91) denominator, expressed as additive feedback gains:
        // y(n−1) gain = +1.9330735, y(n−2) gain = -0.93589199.
        assert!((f.a_coeffs()[1] - 1.933_073_5).abs() < q13);
        assert!((f.a_coeffs()[2] - (-0.935_892)).abs() < q13);
    }

    /// The first output sample is `2·b0·x(0)` — with zero state every
    /// history tap is zero, so `y(0) = b0·x(0)` and the returned value
    /// is the ×2-upscaled `2·b0·x(0)`.
    #[test]
    fn first_sample_is_upscaled_b0_gain() {
        let mut f = OutputHighPass::new();
        let x0 = 100.0;
        let out = f.filter_sample(x0);
        let expected = OUTPUT_UPSCALE * (7699.0 / 8192.0) * x0;
        assert!((out - expected).abs() < 1e-3, "got {out}, want {expected}");
    }

    /// A constant DC input decays to (near) zero output — a high-pass
    /// filter rejects DC. Equivalently, the steady-state gain at `z = 1`
    /// is `(b0+b1+b2)·2 / (1 − a1 − a2)`; with the eq (91) coefficients
    /// the numerator `b0+b1+b2 ≈ 0` (a true zero at DC), so the output
    /// settles to ~0 regardless of the denominator.
    #[test]
    fn rejects_dc() {
        let mut f = OutputHighPass::new();
        let mut last = 0.0;
        for _ in 0..2_000 {
            last = f.filter_sample(1000.0);
        }
        // b0 + b1 + b2 = (7699 - 15398 + 7699)/8192 = 0 exactly, so the
        // DC zero is exact and the steady-state output is ~0.
        assert!(last.abs() < 1.0, "DC not rejected: settled at {last}");
    }

    /// The b-coefficients sum to exactly zero (the DC zero of eq (91)):
    /// `7699 − 15398 + 7699 = 0`. This is the defining property of the
    /// output high-pass and a sign-error tripwire on `b1`.
    #[test]
    fn numerator_has_exact_dc_zero() {
        let f = OutputHighPass::new();
        let sum: f32 = f.b_coeffs().iter().sum();
        assert!(sum.abs() < 1e-6, "b coeffs should sum to 0, got {sum}");
    }

    /// `filter_in_place` and `filter` agree with a per-sample drive on
    /// the same input, and both advance state identically.
    #[test]
    fn batch_apis_match_per_sample() {
        let input: Vec<f32> = (0..40).map(|n| (n as f32 - 20.0) * 3.0).collect();

        let mut f_sample = OutputHighPass::new();
        let per_sample: Vec<f32> = input.iter().map(|&x| f_sample.filter_sample(x)).collect();

        let mut f_batch = OutputHighPass::new();
        let batch = f_batch.filter(&input);

        let mut f_inplace = OutputHighPass::new();
        let mut inplace = input.clone();
        f_inplace.filter_in_place(&mut inplace);

        assert_eq!(per_sample, batch);
        assert_eq!(per_sample, inplace);
        // State matches too.
        assert_eq!(f_sample.x_hist, f_batch.x_hist);
        assert_eq!(f_sample.y_hist, f_batch.y_hist);
        assert_eq!(f_sample.x_hist, f_inplace.x_hist);
        assert_eq!(f_sample.y_hist, f_inplace.y_hist);
    }

    /// State carries across calls: feeding a stream in two halves
    /// produces the same output as feeding it whole.
    #[test]
    fn state_carries_across_calls() {
        let input: Vec<f32> = (0..80).map(|n| ((n * 7) % 13) as f32 - 6.0).collect();

        let mut whole = OutputHighPass::new();
        let out_whole = whole.filter(&input);

        let mut split = OutputHighPass::new();
        let mut out_split = split.filter(&input[..40]);
        out_split.extend(split.filter(&input[40..]));

        assert_eq!(out_whole, out_split);
    }

    /// eq (91) recursion worked by hand for the first three samples of a
    /// unit impulse, pinning the additive-feedback sign convention.
    #[test]
    fn impulse_response_matches_hand_recursion() {
        let mut f = OutputHighPass::new();
        let b0 = 7699.0 / 8192.0;
        let b1 = -15398.0 / 8192.0;
        let b2 = 7699.0 / 8192.0;
        let a1 = 15836.0 / 8192.0;
        let a2 = -7667.0 / 8192.0;

        // x = [1, 0, 0, ...]
        // y(0) = b0
        // y(1) = b1 + a1·y(0)
        // y(2) = b2 + a1·y(1) + a2·y(0)
        let y0 = b0;
        let y1 = b1 + a1 * y0;
        let y2 = b2 + a1 * y1 + a2 * y0;

        let o0 = f.filter_sample(1.0);
        let o1 = f.filter_sample(0.0);
        let o2 = f.filter_sample(0.0);

        assert!((o0 - OUTPUT_UPSCALE * y0).abs() < 1e-6);
        assert!((o1 - OUTPUT_UPSCALE * y1).abs() < 1e-6);
        assert!((o2 - OUTPUT_UPSCALE * y2).abs() < 1e-6);
    }

    /// The filter is BIBO-stable (its poles are inside the unit circle):
    /// a bounded ±8000 input never produces a non-finite or runaway
    /// output across a long stream.
    #[test]
    fn stable_on_bounded_input() {
        let mut f = OutputHighPass::new();
        let mut peak = 0.0f32;
        for n in 0..10_000 {
            let x = if (n / 50) % 2 == 0 { 8000.0 } else { -8000.0 };
            let y = f.filter_sample(x);
            assert!(y.is_finite(), "non-finite output at n={n}");
            peak = peak.max(y.abs());
        }
        // A 100 Hz high-pass on a square wave passes the edges; the
        // ×2-upscaled transient peak stays well bounded (< ~5× the input
        // amplitude), confirming no instability.
        assert!(peak < 5.0 * 8000.0 * OUTPUT_UPSCALE, "runaway peak {peak}");
    }

    /// Two filters over the same input stay in lockstep — all state is
    /// owned, no hidden globals.
    #[test]
    fn deterministic() {
        let input: Vec<f32> = (0..200).map(|n| (n as f32).sin() * 1000.0).collect();
        let mut a = OutputHighPass::new();
        let mut b = OutputHighPass::new();
        assert_eq!(a.filter(&input), b.filter(&input));
    }
}
