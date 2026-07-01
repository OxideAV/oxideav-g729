//! §3.1 **encoder pre-processing** — the two operations applied to the
//! input speech before any analysis stage runs.
//!
//! This is the first *encoder*-side stage of the crate (every earlier
//! round grew the decoder). Per clause 3.1 (transcribed from the EPUB
//! prose): "Two pre-processing functions are applied prior to the
//! encoding process: 1) signal scaling, and 2) high-pass filtering. The
//! scaling consists of dividing the input by a factor 2 to reduce the
//! possibility of overflows in the fixed-point implementation. The
//! high-pass filter serves as a precaution against undesired low
//! frequency components. A second-order pole/zero filter with a cut-off
//! frequency of 140 Hz is used."
//!
//! ## Spec source — clause 3.1, equation (1)
//!
//! Equation (1) (transcribed from the EPUB's equation raster
//! `images/eq1.jpg`) is the transfer function
//!
//! ```text
//!            0.46363718 − 0.92724705·z⁻¹ + 0.46363718·z⁻²
//! H_h1(z) = ────────────────────────────────────────────────
//!            1 − 1.9059465·z⁻¹ + 0.9114024·z⁻²
//! ```
//!
//! The clause-3.1 ÷2 scaling is **folded into the numerator**: the eq (1)
//! numerator coefficients are already the halved values (a non-scaled
//! 140 Hz high-pass of the same shape would have a numerator of
//! `{0.92727…, −1.8544…, 0.92727…}`). So the difference equation
//!
//! ```text
//! s(n) = b0·x(n) + b1·x(n−1) + b2·x(n−2)
//!              + 1.9059465·s(n−1) − 0.9114024·s(n−2)
//! ```
//!
//! consumes the raw 16-bit PCM input `x(n)` directly (no separate ÷2
//! step) and produces the pre-processed signal `s(n)` that clause 3.2
//! (LP analysis) and every later encoder stage operates on.
//!
//! ## Coefficient tables (clause 2.5 Q-format)
//!
//! The 16-bit fixed-point coefficients (clause 2.5: the prose decimals
//! are "rounded versions of the values used in the 16-bit fixed-point
//! implementation") are the already-compiled Q12 tables
//! [`crate::tables::HPF_PREPROC_140HZ_B_Q12`] (`{1899, −3798, 1899}`)
//! and [`crate::tables::HPF_PREPROC_140HZ_A_Q12`] (`{4096, 7807,
//! −3733}`). Converting to real values:
//!
//! * `b = {1899, −3798, 1899} / 2^12 = {0.463623, −0.927246, 0.463623}`
//!   — matches the eq (1) numerator (÷2 folded in).
//! * `a = {4096, 7807, −3733} / 2^12 = {1.0, 1.905762, −0.911377}`. The
//!   first entry is the (unit) `a0`; the remaining two are the feedback
//!   gains **already arranged for an additive recursion** — the
//!   recursion *adds* `a[1]·s(n−1) + a[2]·s(n−2)`, reproducing the eq (1)
//!   denominator `1 − 1.9059465·z⁻¹ + 0.9114024·z⁻²` (the leading `−`
//!   on `z⁻¹` flips to a `+` feedback gain, the `+` on `z⁻²` flips to a
//!   `−` feedback gain), to within the Q12 rounding step. This mirrors
//!   the additive-recursion convention documented for the decoder's
//!   §4.2.5 output high-pass ([`crate::post_process`]).
//!
//! ## State (clause 4.3 init)
//!
//! Per clause 4.3, "all static encoder and decoder variables should be
//! initialized to zero, except the variables listed in Table 9". This
//! filter's memory does not appear in Table 9, so both the two input
//! taps `x(n−1)`, `x(n−2)` and the two output taps `s(n−1)`, `s(n−2)`
//! start zeroed.

use crate::tables::{HPF_PREPROC_140HZ_A_Q12, HPF_PREPROC_140HZ_B_Q12};

/// Q12 fixed-point scale of the eq (1) coefficient tables (`2^12`).
const Q12_SCALE: f32 = 4_096.0;

/// Number of samples in one 10 ms G.729 frame (clause 2).
pub const FRAME_LEN: usize = 80;

/// §3.1 pre-processing high-pass filter `H_h1(z)` — a stateful
/// second-order pole/zero IIR (eq (1)) whose four memory taps carry
/// across frames per clause 4.3 (zero-initialised).
///
/// Feed one frame of raw 16-bit input PCM (as `f32` sample values) to
/// [`Preprocessor::process_frame`]; the returned frame is the
/// pre-processed signal `s(n)` consumed by §3.2 LP analysis.
#[derive(Debug, Clone)]
pub struct Preprocessor {
    /// Numerator coefficients `b0, b1, b2` (real-valued, ÷2 folded in).
    b: [f32; 3],
    /// Additive feedback gains `a1, a2` (real-valued); `a0 = 1` implicit.
    a: [f32; 2],
    /// Input delay line `x(n−1)`, `x(n−2)`.
    x_hist: [f32; 2],
    /// Output delay line `s(n−1)`, `s(n−2)`.
    y_hist: [f32; 2],
}

impl Default for Preprocessor {
    fn default() -> Self {
        Self::new()
    }
}

impl Preprocessor {
    /// Constructs a fresh pre-processor with zeroed state (clause 4.3)
    /// and the eq (1) coefficients dequantised from the compiled Q12
    /// tables.
    #[must_use]
    pub fn new() -> Self {
        let b = [
            f32::from(HPF_PREPROC_140HZ_B_Q12[0]) / Q12_SCALE,
            f32::from(HPF_PREPROC_140HZ_B_Q12[1]) / Q12_SCALE,
            f32::from(HPF_PREPROC_140HZ_B_Q12[2]) / Q12_SCALE,
        ];
        // a[0] is the unit a0 (== 4096 / 4096 == 1.0); the additive
        // feedback gains are entries [1] and [2].
        let a = [
            f32::from(HPF_PREPROC_140HZ_A_Q12[1]) / Q12_SCALE,
            f32::from(HPF_PREPROC_140HZ_A_Q12[2]) / Q12_SCALE,
        ];
        Self {
            b,
            a,
            x_hist: [0.0; 2],
            y_hist: [0.0; 2],
        }
    }

    /// Processes one input sample through the eq (1) high-pass filter,
    /// advancing the four delay taps. Returns the pre-processed sample
    /// `s(n)`.
    #[must_use]
    pub fn process_sample(&mut self, x: f32) -> f32 {
        // s(n) = b0·x(n) + b1·x(n−1) + b2·x(n−2)
        //             + a1·s(n−1) + a2·s(n−2)
        let s = self.b[0] * x
            + self.b[1] * self.x_hist[0]
            + self.b[2] * self.x_hist[1]
            + self.a[0] * self.y_hist[0]
            + self.a[1] * self.y_hist[1];
        // Advance delay lines.
        self.x_hist[1] = self.x_hist[0];
        self.x_hist[0] = x;
        self.y_hist[1] = self.y_hist[0];
        self.y_hist[0] = s;
        s
    }

    /// Processes one 80-sample frame in place, returning the
    /// pre-processed `s(n)` frame. Convenience wrapper over
    /// [`Preprocessor::process_sample`] preserving inter-frame state.
    #[must_use]
    pub fn process_frame(&mut self, frame: &[f32; FRAME_LEN]) -> [f32; FRAME_LEN] {
        let mut out = [0.0f32; FRAME_LEN];
        for (o, &x) in out.iter_mut().zip(frame.iter()) {
            *o = self.process_sample(x);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The dequantised coefficients must match the eq (1) decimals to
    /// within the Q12 rounding step (`1/4096 ≈ 2.4e-4`).
    #[test]
    fn coefficients_match_eq1_decimals() {
        let p = Preprocessor::new();
        assert!((p.b[0] - 0.46363718).abs() < 3e-4);
        assert!((p.b[1] - (-0.92724705)).abs() < 3e-4);
        assert!((p.b[2] - 0.46363718).abs() < 3e-4);
        // Additive feedback gains reproduce the eq (1) denominator
        // 1 − 1.9059465 z⁻¹ + 0.9114024 z⁻²: a1 = +1.9059465,
        // a2 = −0.9114024.
        assert!((p.a[0] - 1.9059465).abs() < 3e-4);
        assert!((p.a[1] - (-0.9114024)).abs() < 3e-4);
    }

    /// A constant DC input must be rejected by the high-pass: after the
    /// transient decays the output settles toward zero.
    #[test]
    fn dc_is_rejected() {
        let mut p = Preprocessor::new();
        let dc = 4000.0f32;
        let mut last = 0.0;
        for _ in 0..4000 {
            last = p.process_sample(dc);
        }
        assert!(
            last.abs() < 1.0,
            "steady-state DC response should decay to ~0, got {last}"
        );
    }

    /// The filter is stable: a bounded input yields a bounded output
    /// (BIBO). A full-scale square wave must not blow up.
    #[test]
    fn bounded_input_bounded_output() {
        let mut p = Preprocessor::new();
        let mut peak = 0.0f32;
        for n in 0..8000 {
            let x = if (n / 40) % 2 == 0 { 30000.0 } else { -30000.0 };
            let s = p.process_sample(x);
            peak = peak.max(s.abs());
        }
        // A 140 Hz high-pass has near-unity passband gain; the square
        // wave's fundamental (100 Hz here) is only mildly attenuated, so
        // the peak stays within a couple of input amplitudes.
        assert!(peak < 90_000.0, "unstable / excessive gain, peak={peak}");
        assert!(peak > 1_000.0, "filter should pass the AC content");
    }

    /// The passband gain near the mid-band is close to **0.5**, not
    /// unity: the clause-3.1 ÷2 input scaling is folded into the eq (1)
    /// numerator, so a 1 kHz sinusoid (well above the 140 Hz corner)
    /// passes at roughly half its input amplitude.
    #[test]
    fn passband_gain_near_half() {
        let mut p = Preprocessor::new();
        let amp = 10000.0f32;
        let mut peak = 0.0f32;
        // 1 kHz at 8 kHz Fs → period 8 samples.
        for n in 0..2000 {
            let x = amp * (std::f32::consts::TAU * (n as f32) * 1000.0 / 8000.0).sin();
            let s = p.process_sample(x);
            if n > 200 {
                peak = peak.max(s.abs());
            }
        }
        let gain = peak / amp;
        assert!(
            (gain - 0.5).abs() < 0.1,
            "1 kHz passband gain should be ~0.5 (÷2 folded in), got {gain}"
        );
    }

    /// Processing a frame sample-by-sample and via `process_frame` must
    /// give identical results (state threading equivalence).
    #[test]
    fn frame_matches_sample_loop() {
        let frame: [f32; FRAME_LEN] = std::array::from_fn(|i| ((i * 37 % 500) as f32) - 250.0);

        let mut p1 = Preprocessor::new();
        let out_frame = p1.process_frame(&frame);

        let mut p2 = Preprocessor::new();
        let mut out_loop = [0.0f32; FRAME_LEN];
        for (o, &x) in out_loop.iter_mut().zip(frame.iter()) {
            *o = p2.process_sample(x);
        }
        assert_eq!(out_frame, out_loop);
    }
}
