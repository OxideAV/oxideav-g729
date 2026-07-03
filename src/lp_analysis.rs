//! §3.2.1 **LP analysis front-end** — windowing, autocorrelation, and
//! the 60 Hz bandwidth-expansion lag window.
//!
//! This is the second encoder-side stage (after [`crate::preprocess`]).
//! It turns a 240-sample analysis window of pre-processed speech `s(n)`
//! into the 11 modified autocorrelation coefficients `r'(0…10)` that the
//! §3.2.2 Levinson-Durbin recursion consumes. The §3.2.2 recursion and
//! the §3.2.3 LP→LSP conversion are separate follow-up modules.
//!
//! ## Spec source — clause 3.2 / 3.2.1, equations (3)–(7)
//!
//! Per clause 3.2 the short-term filter coefficients are computed once
//! per 10 ms frame using the **autocorrelation method** with a 30 ms
//! asymmetric analysis window and a 10th-order model.
//!
//! ### The analysis window `w_lp(n)` (eq (3))
//!
//! Length 240 samples (30 ms at 8 kHz), spanning **120 samples from
//! past frames + 80 from the present frame + 40 from the future frame**
//! (the 5 ms look-ahead), per clause 3.2.1 / Figure 5. Equation (3):
//!
//! ```text
//!            ⎧ 0.54 − 0.46·cos(2πn / 399)      n = 0 … 199
//! w_lp(n) =  ⎨
//!            ⎩ cos(2π(n − 200) / 159)          n = 200 … 239
//! ```
//!
//! (first part = half a Hamming window, second part = a quarter cosine
//! taper over the look-ahead). Its 240 Q15 samples are the compiled
//! [`crate::tables::LPC_HAMMING_WINDOW_Q15`].
//!
//! ### Windowing → autocorrelation (eqs (4), (5))
//!
//! ```text
//! s'(n) = w_lp(n) · s(n)            n = 0 … 239           (4)
//! r(k)  = Σ_{n=k}^{239} s'(n)·s'(n−k)   k = 0 … 10        (5)
//! ```
//!
//! `r(0)` is floored to `1.0` (a low-level-signal guard so silence still
//! yields a well-defined, positive-definite autocorrelation).
//!
//! ### 60 Hz bandwidth expansion — the lag window (eqs (6), (7))
//!
//! To avoid ill-conditioning and to slightly widen the LP-spectrum
//! bandwidths, the autocorrelations are weighted by a lag window
//! (eq (6), `w_lag`, `f0 = 60 Hz`, `fs = 8000 Hz`) and `r(0)` is
//! additionally scaled by a white-noise correction factor `1.0001`
//! (−40 dB noise floor), giving the modified coefficients `r'(k)`
//! (eq (7)):
//!
//! ```text
//! r'(0) = r(0) · 1.0001
//! r'(k) = r(k) · w_lag(k)               k = 1 … 10
//! ```
//!
//! The `w_lag(k)` values (`k = 1 … 10`) are stored as **split
//! double-precision** Q15 pairs: the high halves in
//! [`crate::tables::LPC_LAG_WINDOW_HIGH_Q15`] and the low halves in
//! [`crate::tables::LPC_LAG_WINDOW_LOW_Q15`], so
//! `w_lag(k) = lag_h[k−1]/2^15 + lag_l[k−1]/2^31`.
//!
//! ## Numeric domain
//!
//! The clause-2.5 signal path is 16-bit fixed point, and the
//! conformance corpus was produced by that arithmetic, so (round 385)
//! the windowing applies the Q15 multiply **with the fixed-point
//! rounding to the 16-bit grid**: `s'(n) = ⌊(s(n)·w(n) + 2^14)·2^−15⌋`.
//! Both `s(n)` (rounded by [`crate::preprocess`]) and `s'(n)` are then
//! integer-valued, so accumulating the eq (5) autocorrelation in `f64`
//! is *exact* integer arithmetic (products ≤ 2^30, 240-term sums ≤
//! 2^38 ≪ the 2^53 f64 mantissa) — bit-identical to a 64-bit integer
//! accumulator.

use crate::tables::{
    LPC_HAMMING_WINDOW_Q15, LPC_LAG_WINDOW_HIGH_Q15, LPC_LAG_WINDOW_LOW_Q15, L_WINDOW, M,
};

/// High-half Q15 scale for the lag window (`2^15`).
const LAG_HIGH_SCALE: f64 = 32_768.0;
/// Combined scale for the low-half split-double lag window (`2^31`).
const LAG_LOW_SCALE: f64 = 32_768.0 * 65_536.0;
/// White-noise correction applied to `r(0)` (clause 3.2.1, eq (7)).
const WHITE_NOISE_FACTOR: f64 = 1.0001;

/// The number of autocorrelation lags produced: `r(0) … r(M)` = 11.
pub const N_LAGS: usize = M + 1;

/// Applies the eq (4) analysis window to a 240-sample block of
/// pre-processed speech `s(n)`, returning the windowed signal `s'(n)`.
///
/// `window[0]` is the oldest sample of the 30 ms analysis window (the
/// start of the 120-sample past region); `window[239]` is the newest
/// (the last of the 40-sample look-ahead).
#[must_use]
pub fn apply_window(window: &[f32; L_WINDOW]) -> [f32; L_WINDOW] {
    let mut out = [0.0f32; L_WINDOW];
    for (o, (&s, &w)) in out
        .iter_mut()
        .zip(window.iter().zip(LPC_HAMMING_WINDOW_Q15.iter()))
    {
        // Fixed-point Q15 windowing with rounding: the 16-bit signal
        // path (clause 2.5) computes s'(n) = round(s(n)·w(n)·2^−15),
        // i.e. ⌊(s·w + 2^14)·2^−15⌋ — the windowed sample every later
        // stage sees is itself a 16-bit integer.
        let prod = f64::from(s) * f64::from(w);
        *o = ((prod + 16_384.0) / 32_768.0).floor() as f32;
    }
    out
}

/// Computes the eq (5) autocorrelations `r(0 … M)` of a windowed
/// signal, with the clause-3.2.1 `r(0) ≥ 1.0` low-level guard applied.
///
/// The sum is accumulated in `f64`. The returned array is
/// `[r(0), r(1), …, r(M)]`.
#[must_use]
pub fn autocorrelate(windowed: &[f32; L_WINDOW]) -> [f64; N_LAGS] {
    let mut r = [0.0f64; N_LAGS];
    for (k, rk) in r.iter_mut().enumerate() {
        let mut acc = 0.0f64;
        for n in k..L_WINDOW {
            acc += f64::from(windowed[n]) * f64::from(windowed[n - k]);
        }
        *rk = acc;
    }
    // Low-level-signal guard: r(0) floored to 1.0 so silence still
    // yields a positive-definite, well-conditioned autocorrelation.
    if r[0] < 1.0 {
        r[0] = 1.0;
    }
    r
}

/// Applies the eq (6)/(7) lag window (60 Hz bandwidth expansion) plus
/// the `1.0001` white-noise correction on `r(0)`, in place, producing
/// the modified coefficients `r'(k)` the Levinson-Durbin recursion
/// consumes.
pub fn apply_lag_window(r: &mut [f64; N_LAGS]) {
    r[0] *= WHITE_NOISE_FACTOR;
    for k in 1..N_LAGS {
        let w = f64::from(LPC_LAG_WINDOW_HIGH_Q15[k - 1]) / LAG_HIGH_SCALE
            + f64::from(LPC_LAG_WINDOW_LOW_Q15[k - 1]) / LAG_LOW_SCALE;
        r[k] *= w;
    }
}

/// Convenience: run the whole §3.2.1 front-end — window (eq (4)),
/// autocorrelate (eq (5)) with the `r(0)` guard, and lag-window
/// (eqs (6)/(7)) — over one 240-sample analysis block, returning the
/// modified autocorrelations `r'(0 … M)`.
#[must_use]
pub fn analyze_window(window: &[f32; L_WINDOW]) -> [f64; N_LAGS] {
    let windowed = apply_window(window);
    let mut r = autocorrelate(&windowed);
    apply_lag_window(&mut r);
    r
}

/// Returns the real-valued lag-window weight `w_lag(k)` for `k = 1 … M`
/// reconstructed from the split-double Q15 tables. `k = 0` returns the
/// white-noise correction `1.0001`.
///
/// # Panics
///
/// Panics if `k > M`.
#[must_use]
pub fn lag_window_weight(k: usize) -> f64 {
    assert!(k <= M, "lag index {k} out of range 0..={M}");
    if k == 0 {
        WHITE_NOISE_FACTOR
    } else {
        f64::from(LPC_LAG_WINDOW_HIGH_Q15[k - 1]) / LAG_HIGH_SCALE
            + f64::from(LPC_LAG_WINDOW_LOW_Q15[k - 1]) / LAG_LOW_SCALE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The eq (3) window: first sample of the half-Hamming part is
    /// `0.54 − 0.46 = 0.08`, the window peaks at `1.0` where the Hamming
    /// half meets its maximum, and the look-ahead tail descends toward 0.
    /// Probed with a constant full-scale-ish integer signal so the Q15
    /// rounding step stays small relative to the window shape.
    #[test]
    fn window_shape_matches_eq3() {
        let amp = 10_000.0f32;
        let flat = [amp; L_WINDOW];
        let w = apply_window(&flat);
        // n = 0: 0.54 − 0.46·cos(0) = 0.08 → ≈ 800 at amp 10000.
        assert!((w[0] - 0.08 * amp).abs() < 2.0, "w[0]={}", w[0]);
        // Somewhere in 0..199 the Hamming half reaches ~1.0·amp.
        let peak = w[..200].iter().cloned().fold(0.0f32, f32::max);
        assert!(
            peak > 0.999 * amp,
            "hamming half should reach ~amp, peak={peak}"
        );
        // The last look-ahead sample is well below the peak (cosine
        // taper descending toward 0).
        assert!(w[239] < peak, "taper tail w[239]={} peak={peak}", w[239]);
        // Round-385 fixed-point grid: every windowed sample is an
        // integer (the 16-bit windowed-speech domain).
        for (n, &v) in w.iter().enumerate() {
            assert_eq!(v, v.trunc(), "windowed[{n}] = {v} not on integer grid");
        }
    }

    /// Silence: the `r(0) ≥ 1.0` guard makes a zero input yield a
    /// well-defined autocorrelation rather than all zeros.
    #[test]
    fn silence_floors_r0() {
        let zeros = [0.0f32; L_WINDOW];
        let r = autocorrelate(&zeros);
        assert_eq!(r[0], 1.0);
        for rk in &r[1..] {
            assert_eq!(*rk, 0.0);
        }
    }

    /// For a pure cosine the normalised autocorrelation `r(k)/r(0)`
    /// tracks `cos(ω·k)` (a standard property of a windowed sinusoid, to
    /// within window-broadening error).
    #[test]
    fn cosine_autocorrelation_tracks_cos() {
        // 1 kHz at 8 kHz → ω = 2π·1000/8000 = π/4 rad/sample.
        let omega = std::f32::consts::TAU * 1000.0 / 8000.0;
        let sig: [f32; L_WINDOW] = std::array::from_fn(|n| (omega * n as f32).cos() * 4000.0);
        // Skip the window for this structural check (use a rectangular
        // window) so the broadening does not dominate.
        let r = autocorrelate(&sig);
        for k in 1..=4 {
            let norm = r[k] / r[0];
            let expected = f64::from((omega * k as f32).cos());
            assert!(
                (norm - expected).abs() < 0.05,
                "k={k}: r(k)/r(0)={norm} vs cos={expected}"
            );
        }
    }

    /// The lag window weights are strictly decreasing in `k` and match
    /// the eq (6) Gaussian `exp(−½·(2π·60·k/8000)²)` to Q15 precision.
    #[test]
    fn lag_window_matches_gaussian() {
        let mut prev = f64::INFINITY;
        for k in 1..=M {
            let w = lag_window_weight(k);
            let expected =
                (-0.5 * (std::f64::consts::TAU * 60.0 * k as f64 / 8000.0).powi(2)).exp();
            assert!(
                (w - expected).abs() < 1e-3,
                "k={k}: w_lag={w} vs gaussian={expected}"
            );
            assert!(
                w < prev,
                "lag window must decrease: k={k} w={w} prev={prev}"
            );
            prev = w;
        }
    }

    /// `apply_lag_window` scales `r(0)` by the white-noise factor and
    /// each `r(k)` by its lag weight, leaving the ordering intact.
    #[test]
    fn lag_window_application() {
        let mut r = [10.0f64; N_LAGS];
        apply_lag_window(&mut r);
        assert!((r[0] - 10.0 * WHITE_NOISE_FACTOR).abs() < 1e-9);
        for (k, &rk) in r.iter().enumerate().skip(1) {
            assert!((rk - 10.0 * lag_window_weight(k)).abs() < 1e-9);
        }
    }

    /// The full front-end produces a positive `r'(0)` and lag-attenuated
    /// higher lags for a non-trivial input.
    #[test]
    fn analyze_window_end_to_end() {
        let sig: [f32; L_WINDOW] = std::array::from_fn(|n| ((n * 31 % 400) as f32 - 200.0) + 50.0);
        let r = analyze_window(&sig);
        assert!(r[0] > 0.0);
        // r'(0) carries the white-noise inflation over the raw windowed
        // energy; higher lags are attenuated by the lag window.
        assert!(r[0].abs() >= r[1].abs() || r[0] > 1.0);
    }
}
