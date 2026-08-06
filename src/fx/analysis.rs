//! §3.2.1 **fixed-point LP-analysis front end** — windowing,
//! autocorrelation, and the 60 Hz lag window on the clause-5
//! Word16/Word32 operator grid.
//!
//! The float module [`crate::lp_analysis`] emulates the 16-bit grid by
//! rounding the eq (4) windowed samples to integers and accumulating
//! the eq (5) sums exactly in `f64`, then normalising the finished
//! sums onto the Word32 grid. A genuine 16/32-bit pipeline cannot do
//! that: the eq (5) accumulator *is* a saturating Word32, so a loud
//! analysis window overflows mid-sum and the signal itself must be
//! scaled down before the sum is re-run. The surviving low bits of
//! each product therefore differ from the float emulation on exactly
//! the loud frames where the encoder's §3.2.4 index agreement drops —
//! this module implements the genuine protocol so that latitude can
//! be pinned black-box against the conformance corpus.
//!
//! ## Spec source
//!
//! - eq (4) windowing `s'(n) = s(n)·w_lp(n)` on the Q15 grid with
//!   rounding (`mult_r`), per the clause-2.5 16-bit signal path.
//! - eq (5) autocorrelation `r(k) = Σ s'(n)·s'(n−k)` in a Word32
//!   accumulator (`l_mac`), with the low-level `r(0) ≥ 1` guard of
//!   clause 3.2.1.
//! - eqs (6)/(7) lag window applied as a DPF 32×32 product
//!   (`mpy_32`) against the staged split-double Q15 tables; the
//!   white-noise correction on `r(0)` stays at its corpus-pinned
//!   effective unity (see `WHITE_NOISE_FACTOR` in
//!   [`crate::lp_analysis`]).
//!
//! ## The overflow-rescale protocol (unstated latitude)
//!
//! The Recommendation prints the equations, not the 16-bit scaling
//! schedule. The protocol implemented here — detect Word32 saturation
//! of the `r(0)` accumulation, right-shift the windowed signal, retry
//! — is the only scheme expressible in the Table 10/11 operator set
//! without widening the accumulator; the *step* of that shift is
//! unstated and exposed as [`FrontEndLatitude::overflow_shift`],
//! pinned black-box against the corpus (see the conformance tests).

use crate::fx::dsp::{l_comp, l_extract, mpy_32, mpy_32_16};
use crate::fx::ops::{l_add, l_mac, l_mult, l_shl, mult_r, norm_l, round, shr};
use crate::lp_analysis::N_LAGS;
use crate::tables::{
    HPF_PREPROC_140HZ_A_Q12, HPF_PREPROC_140HZ_B_Q12, LPC_HAMMING_WINDOW_Q15,
    LPC_LAG_WINDOW_HIGH_Q15, LPC_LAG_WINDOW_LOW_Q15, L_WINDOW,
};

/// A Word32 value carried in the DPF `(hi, lo)` split
/// (`L = hi·2^16 + lo·2`).
pub type Dpf = (i16, i16);

/// Unstated fixed-point latitude of the front end, pinned black-box.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FrontEndLatitude {
    /// Right-shift applied to the windowed signal on each Word32
    /// overflow of the eq (5) `r(0)` accumulation before the retry.
    pub overflow_shift: i16,
}

impl Default for FrontEndLatitude {
    fn default() -> Self {
        Self { overflow_shift: 2 }
    }
}

/// §3.1 fixed-point pre-processing — the eq (1) 140 Hz high-pass
/// (÷2 folded into the Q12 numerator) on the Word16/Word32 grid,
/// mirroring the structure the decoder's §4.2.5 output high-pass
/// pinned black-box: a Q13 accumulator (`l_mult`/`l_mac` over the
/// staged Q12 coefficients), the feedback memory kept as the
/// **unrounded Word32** recursion value on the Q16 grid
/// (`mpy_32_16` against the Q12 feedback gains), and the output
/// rounded to the saturated 16-bit PCM grid.
#[derive(Debug, Clone, Default)]
pub struct PreprocessorFx {
    /// Input taps `x(n−1)`, `x(n−2)`.
    x_hist: [i16; 2],
    /// Unrounded feedback `y(n−1)`, `y(n−2)` on Word32 Q16.
    y_hist: [i32; 2],
}

impl PreprocessorFx {
    /// Fresh filter with the clause-4.3 zero state.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Processes one raw 16-bit PCM sample, returning `s(n)` on the
    /// 16-bit grid.
    #[must_use]
    pub fn process_sample(&mut self, x: i16) -> i16 {
        let b = &HPF_PREPROC_140HZ_B_Q12;
        let a = &HPF_PREPROC_140HZ_A_Q12;
        // Numerator on the Q13 accumulator: b (Q12) × x (Q0) doubled.
        let mut acc = l_mult(x, b[0]);
        acc = l_mac(acc, self.x_hist[0], b[1]);
        acc = l_mac(acc, self.x_hist[1], b[2]);
        // Feedback: a (Q12) × y (Q16 Word32) lands on Q13 directly
        // (corpus-pinned: a Q28 DPF coefficient variant from the
        // printed eq (1) decimals costs ~8/10 points of LSP/SPEECH,
        // and a Q13 y-storage costs ~12/16 — the staged Q12 grid with
        // three extra feedback fraction bits is the reference
        // arrangement).
        let (h1, l1) = l_extract(self.y_hist[0]);
        acc = l_add(acc, mpy_32_16(h1, l1, a[1]));
        let (h2, l2) = l_extract(self.y_hist[1]);
        acc = l_add(acc, mpy_32_16(h2, l2, a[2]));
        self.y_hist[1] = self.y_hist[0];
        self.y_hist[0] = l_shl(acc, 3);
        self.x_hist[1] = self.x_hist[0];
        self.x_hist[0] = x;
        // Q13 → Q16 high word: rounded, saturated 16-bit output.
        round(l_shl(acc, 3))
    }

    /// Processes one 80-sample frame, preserving inter-frame state.
    #[must_use]
    pub fn process_frame<const N: usize>(&mut self, frame: &[i16; N]) -> [i16; N] {
        std::array::from_fn(|n| self.process_sample(frame[n]))
    }
}

/// eq (4) windowing on the 16-bit grid: `s'(n) = mult_r(s(n), w(n))`
/// (round-to-nearest Q15 product, exactly the float emulation's
/// `⌊(s·w + 2^14)·2^−15⌋` for in-range inputs).
#[must_use]
pub fn window_fx(speech: &[i16; L_WINDOW]) -> [i16; L_WINDOW] {
    std::array::from_fn(|n| mult_r(speech[n], LPC_HAMMING_WINDOW_Q15[n]))
}

/// Result of the eq (5) autocorrelation on the Word32 grid.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AutocorrFx {
    /// `r(0 … M)` as normalised DPF pairs (common scale; `r(0)`'s
    /// Word32 recomposition sits in `[2^30, 2^31)`).
    pub r: [Dpf; N_LAGS],
    /// Total right-shift applied to the windowed signal by the
    /// overflow-rescale protocol (0 when no overflow occurred).
    pub signal_shift: i16,
    /// The `norm_l` normalisation applied to the finished sums.
    pub norm: i16,
}

/// eq (5) autocorrelation with the overflow-rescale protocol and the
/// clause-3.2.1 `r(0) ≥ 1` low-level guard, returning normalised DPF
/// coefficients.
///
/// The `r(0)` sum runs in a Word32 accumulator; if it saturates, the
/// windowed signal is right-shifted by `lat.overflow_shift` and the
/// sum re-run (repeating as needed). The higher lags are then bounded
/// by the finished `r(0)` (Cauchy-Schwarz on identical inputs), so
/// they accumulate without further rescue, and all eleven sums are
/// left-shifted by `norm_l(r(0))` onto the normalised grid.
#[must_use]
pub fn autocorr_fx(windowed: &[i16; L_WINDOW], lat: &FrontEndLatitude) -> AutocorrFx {
    let mut y = *windowed;
    let mut signal_shift: i16 = 0;

    // r(0) with overflow detection: the exact integer sum is compared
    // against the Word32 ceiling — behaviourally identical to running
    // the saturating `l_mac` chain and testing its overflow flag.
    let r0: i32 = loop {
        let exact: i64 = y.iter().map(|&v| i64::from(v) * i64::from(v) * 2).sum();
        if exact <= i64::from(i32::MAX) {
            break exact as i32;
        }
        for v in &mut y {
            *v = shr(*v, lat.overflow_shift);
        }
        signal_shift += lat.overflow_shift;
    };
    // Low-level guard: clause 3.2.1 floors r(0) so silence still
    // yields a positive-definite autocorrelation.
    let r0 = r0.max(1);

    let norm = norm_l(r0);
    let mut r = [(0i16, 0i16); N_LAGS];
    r[0] = l_extract(l_shl(r0, norm));
    for k in 1..N_LAGS {
        let mut acc: i32 = 0;
        for n in k..L_WINDOW {
            // l_mac without rescue: |partial| ≤ r(0) by Cauchy-Schwarz.
            acc += i32::from(y[n]) * i32::from(y[n - k]) * 2;
        }
        r[k] = l_extract(l_shl(acc, norm));
    }
    AutocorrFx {
        r,
        signal_shift,
        norm,
    }
}

/// eqs (6)/(7) lag window on the DPF grid: each `r(k)` (`k ≥ 1`) is
/// multiplied by the staged split-double weight via `mpy_32`; `r(0)`
/// keeps the corpus-pinned effective-unity white-noise correction.
pub fn lag_window_fx(r: &mut [Dpf; N_LAGS]) {
    for k in 1..N_LAGS {
        let (hi, lo) = r[k];
        let prod = mpy_32(
            hi,
            lo,
            LPC_LAG_WINDOW_HIGH_Q15[k - 1],
            LPC_LAG_WINDOW_LOW_Q15[k - 1],
        );
        r[k] = l_extract(prod);
    }
}

/// Convenience: the whole §3.2.1 fixed-point front end — eq (4)
/// windowing, eq (5) autocorrelation with the overflow protocol,
/// eqs (6)/(7) lag window — over one 240-sample analysis block.
#[must_use]
pub fn analyze_window_fx(speech: &[i16; L_WINDOW], lat: &FrontEndLatitude) -> AutocorrFx {
    let windowed = window_fx(speech);
    let mut out = autocorr_fx(&windowed, lat);
    lag_window_fx(&mut out.r);
    out
}

/// Recomposes a DPF autocorrelation vector into `f64` (test/AB use).
#[must_use]
pub fn dpf_to_f64(r: &[Dpf; N_LAGS]) -> [f64; N_LAGS] {
    std::array::from_fn(|k| f64::from(l_comp(r[k].0, r[k].1)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_analysis::{apply_window, autocorrelate, lag_window_weight};

    fn test_signal(amp: f32) -> [f32; L_WINDOW] {
        std::array::from_fn(|n| {
            let n = n as f32;
            (amp * (0.07 * n).sin() + 0.25 * amp * (0.31 * n).cos()).round()
        })
    }

    /// §3.1: the fixed pre-processor tracks the float filter within
    /// one PCM LSB (the Word32-vs-f32 feedback truncation) on a
    /// speech-band test signal, and both reject DC.
    #[test]
    fn preprocessor_tracks_float_and_rejects_dc() {
        use crate::preprocess::Preprocessor;
        let mut fx = PreprocessorFx::new();
        let mut fl = Preprocessor::new();
        let mut max_d = 0i32;
        for n in 0..800 {
            let x = (8_000.0 * (0.4 * n as f32).sin() + 3_000.0 * (0.9 * n as f32).cos()).round();
            let a = i32::from(fx.process_sample(x as i16));
            let b = fl.process_sample(x) as i32;
            max_d = max_d.max((a - b).abs());
        }
        assert!(max_d <= 1, "fx vs float preprocess max delta {max_d}");
        // DC rejection: a long constant input decays toward zero.
        let mut fx = PreprocessorFx::new();
        let mut last = 0i16;
        for _ in 0..2_000 {
            last = fx.process_sample(8_192);
        }
        assert!(last.abs() <= 1, "DC leak {last}");
    }

    /// eq (4): the fixed windowing matches the float emulation's
    /// integer grid sample for sample.
    #[test]
    fn window_matches_float_emulation() {
        let sig = test_signal(11_000.0);
        let sig_i: [i16; L_WINDOW] = std::array::from_fn(|n| sig[n] as i16);
        let float = apply_window(&sig);
        let fixed = window_fx(&sig_i);
        for n in 0..L_WINDOW {
            assert_eq!(
                f32::from(fixed[n]),
                float[n],
                "windowed sample {n} off the float-emulated grid"
            );
        }
    }

    /// No-overflow case: the DPF pipeline reproduces the float
    /// emulation's normalised sums exactly (both are exact integer
    /// arithmetic under the same normalisation).
    #[test]
    fn quiet_signal_is_exact() {
        let sig = test_signal(400.0);
        let sig_i: [i16; L_WINDOW] = std::array::from_fn(|n| sig[n] as i16);
        let windowed = apply_window(&sig);
        let float_r = autocorrelate(&windowed);
        let fx = autocorr_fx(&window_fx(&sig_i), &FrontEndLatitude::default());
        assert_eq!(fx.signal_shift, 0, "quiet signal must not trigger rescue");
        let got = dpf_to_f64(&fx.r);
        let scale = (2f64).powi(i32::from(fx.norm)) * 2.0;
        for k in 0..N_LAGS {
            assert_eq!(
                got[k],
                float_r[k] * scale,
                "lag {k}: fixed {} vs float {}",
                got[k],
                float_r[k] * scale
            );
        }
    }

    /// Loud case: the overflow protocol engages, `r(0)` stays
    /// normalised, and the result tracks the exact sums within the
    /// truncation loss of the rescue shift.
    #[test]
    fn loud_signal_engages_overflow_protocol() {
        let sig: [f32; L_WINDOW] =
            std::array::from_fn(|n| test_signal(30_000.0)[n].clamp(-32_768.0, 32_767.0));
        let sig_i: [i16; L_WINDOW] = std::array::from_fn(|n| sig[n] as i16);
        let fx = autocorr_fx(&window_fx(&sig_i), &FrontEndLatitude::default());
        assert!(fx.signal_shift > 0, "loud signal must trigger rescue");
        let r0 = f64::from(l_comp(fx.r[0].0, fx.r[0].1));
        assert!(
            ((2f64).powi(30)..(2f64).powi(31)).contains(&r0),
            "r(0) = {r0} not normalised"
        );
        // Against the exact (float-emulated) sums the relative error is
        // bounded by the rescue truncation (well under 2%).
        let windowed = apply_window(&sig);
        let exact = autocorrelate(&windowed);
        let got = dpf_to_f64(&fx.r);
        for k in 0..N_LAGS {
            let want = exact[k] / exact[0];
            let have = got[k] / got[0];
            assert!(
                (want - have).abs() < 2e-2,
                "lag {k}: normalised {have} vs exact {want}"
            );
        }
    }

    /// Silence: the r(0) guard yields a well-defined normalised value.
    #[test]
    fn silence_guard() {
        let fx = autocorr_fx(&[0i16; L_WINDOW], &FrontEndLatitude::default());
        let r0 = l_comp(fx.r[0].0, fx.r[0].1);
        assert!(r0 >= 1 << 30, "guarded r(0) must normalise, got {r0}");
        for k in 1..N_LAGS {
            assert_eq!(l_comp(fx.r[k].0, fx.r[k].1), 0);
        }
    }

    /// The DPF lag window tracks the real-valued weights within the
    /// `mpy_32` truncation (the dropped `lo·lo` cross term).
    #[test]
    fn lag_window_tracks_weights() {
        let sig = test_signal(6_000.0);
        let sig_i: [i16; L_WINDOW] = std::array::from_fn(|n| sig[n] as i16);
        let mut fx = autocorr_fx(&window_fx(&sig_i), &FrontEndLatitude::default());
        let before = dpf_to_f64(&fx.r);
        lag_window_fx(&mut fx.r);
        let after = dpf_to_f64(&fx.r);
        assert_eq!(after[0], before[0], "r(0) must keep effective unity");
        for k in 1..N_LAGS {
            let want = before[k] * lag_window_weight(k);
            // mpy_32 truncates each partial product: bounded absolute
            // error of a few Word32 LSB·2^16.
            assert!(
                (after[k] - want).abs() <= 3.0 * 65_536.0,
                "lag {k}: {} vs {}",
                after[k],
                want
            );
        }
    }
}
