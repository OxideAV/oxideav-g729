//! §4.1.6 **LP synthesis** — the first decoder stage to emit
//! reconstructed-speech PCM.
//!
//! This module sits at the tail of the clause-4.1 decode order. The
//! round-282 [`crate::decode_chain::FrameDecoder`] turns one ITU serial
//! frame into a [`crate::decode_chain::DecodedFrame`] carrying, per
//! subframe, the §4.1.3 fractional pitch delay, the §4.1.4 (post-eq
//! (48)) fixed-codebook vector `c(n)`, the §4.1.5 gains `ĝ_p` / `ĝ_c`,
//! and the §4.1.1 interpolated LP filter coefficients `â_i`. This
//! module consumes those products and runs the spec §4.1.3 → §3.10 →
//! §4.1.6 chain to produce the reconstructed speech `ŝ(n)`.
//!
//! ## Spec source — clauses 3.7.1, 3.10, 4.1.6 (06/2012 Recommendation)
//!
//! The three load-bearing equations (transcribed from the EPUB's
//! equation-JPGs, which the prose renders as raster images):
//!
//! * **eq (40)** (clause 3.7.1) — the adaptive-codebook vector `v(n)`,
//!   built by interpolating the past excitation `u(n)` at integer
//!   delay `k` and fraction `t ∈ {0, 1, 2}` through the 31-tap `b_30`
//!   synthesis filter (Hamming-windowed sinc truncated at ±29, padded
//!   with a zero at ±30, −3 dB at 3600 Hz):
//!
//!   ```text
//!   v(n) = Σ_{i=0}^{9} u(n − k − i)·b30(t + 3i)
//!        + Σ_{i=0}^{9} u(n − k + 1 + i)·b30(3 − t + 3i)
//!        for n = 0 … 39,  t = 0, 1, 2
//!   ```
//!
//! * **eq (75)** (clause 3.10) — the subframe excitation, the sum of
//!   the scaled adaptive- and fixed-codebook contributions:
//!
//!   ```text
//!   u(n) = ĝ_p·v(n) + ĝ_c·c(n)   for n = 0 … 39
//!   ```
//!
//!   where `c(n)` is "the fixed-codebook vector including harmonic
//!   enhancement" (the eq (48)-modified vector the decode chain
//!   already produced).
//!
//! * **eq (77)** (clause 4.1.6) — the reconstructed speech, the
//!   excitation filtered through the 10th-order LP synthesis filter
//!   `1/Â(z)`:
//!
//!   ```text
//!   ŝ(n) = u(n) − Σ_{i=1}^{10} â_i·ŝ(n − i)   for n = 0 … 39
//!   ```
//!
//!   where `â_i` are "the interpolated LP filter coefficients for the
//!   current subframe".
//!
//! ## Cross-subframe state (clause 4.3 init)
//!
//! Per clause 4.3, "all static encoder and decoder variables should be
//! initialized to zero, except the variables listed in Table 9".
//! Neither of this module's two state pieces appears in Table 9, so
//! both start zeroed:
//!
//! | state | role | init |
//! |-------|------|------|
//! | [`Synthesizer::exc_history`] | eq (40) past-excitation buffer `u(n)`, `n < 0` | all zero |
//! | [`Synthesizer::syn_mem`] | eq (77) 10th-order `1/Â(z)` filter memory `ŝ(n)`, `n < 0` | all zero |
//!
//! The past-excitation buffer holds [`EXC_HISTORY`] = 153 samples —
//! the eq (40) sum reaches `u(n − k − i)` with `n = 0`, the maximum
//! integer delay `k = 143` (spec §3.7 / §4.1.3 cap), and `i = 9`,
//! i.e. index `−152`; 153 history samples (indices `−152 … −1`) cover
//! every access. The spec excitation buffer is `u(n)`,
//! `n = −143 … 39` (clause 3.7); the extra 9 samples below `−143`
//! are the interpolation tail of eq (40).
//!
//! ## What this module does NOT do
//!
//! The §4.2 post-processing cascade (long-term + short-term postfilter,
//! tilt compensation, adaptive gain control, output high-pass, ×2
//! upscaling) is a follow-up round; this module emits the raw §4.1.6
//! `ŝ(n)` that feeds the §4.2 postfilter. The §3.10 weighting-filter
//! memory update (eq (76) `ew(n)`) is an *encoder*-side memory and is
//! not part of the decoder synthesis path.

use crate::decode_chain::DecodedFrame;
use crate::fixed_codebook::SUBFRAME_SIZE;
use crate::tables::{M, PITCH_INTERP_FILTER_SYNTHESIS_Q15};

/// Number of past-excitation samples retained across subframes for the
/// eq (40) `b_30` interpolation. The deepest access is `u(n − k − i)`
/// with `n = 0`, `k = 143` (the §3.7 / §4.1.3 maximum integer pitch
/// delay), `i = 9`, i.e. index `−152`; 153 history samples (indices
/// `−152 … −1`) cover it.
pub const EXC_HISTORY: usize = 153;

/// eq (40) interpolation half-length: each of the two sums runs
/// `i = 0 … 9` (10 taps), matching the `b_30` filter's ±29-truncated
/// sinc decimated by 3.
const INTERP_TAPS: usize = 10;

/// `b_30` Q15 fixed-point scale (the synthesis interpolation filter is
/// stored as Q15 `i16`; eq (40) divides the integer convolution by
/// `2^15`).
const B30_Q15_SCALE: f32 = 32_768.0;

/// Everything §4.1.3 → §3.10 → §4.1.6 produces for **one subframe**.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SynthesizedSubframe {
    /// eq (40) adaptive-codebook vector `v(n)` (interpolated past
    /// excitation), `n = 0 … 39`.
    pub adaptive: [f32; SUBFRAME_SIZE],
    /// eq (75) excitation `u(n) = ĝ_p·v(n) + ĝ_c·c(n)`, `n = 0 … 39`.
    pub excitation: [f32; SUBFRAME_SIZE],
    /// eq (77) reconstructed speech `ŝ(n)`, `n = 0 … 39`.
    pub speech: [f32; SUBFRAME_SIZE],
}

/// One 10 ms frame's worth of reconstructed speech — two 40-sample
/// subframes, in spec order (index 0 = subframe 1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SynthesizedFrame {
    /// Per-subframe synthesis products in spec §4.1 order.
    pub subframes: [SynthesizedSubframe; 2],
}

impl SynthesizedFrame {
    /// The 80 reconstructed-speech samples of the frame in time order
    /// (subframe 1 then subframe 2).
    #[must_use]
    pub fn speech(&self) -> [f32; 2 * SUBFRAME_SIZE] {
        let mut out = [0.0f32; 2 * SUBFRAME_SIZE];
        out[..SUBFRAME_SIZE].copy_from_slice(&self.subframes[0].speech);
        out[SUBFRAME_SIZE..].copy_from_slice(&self.subframes[1].speech);
        out
    }
}

/// Stateful §4.1.6 LP synthesizer.
///
/// Owns the two cross-subframe state pieces (the eq (40) past-
/// excitation buffer and the eq (77) 10th-order filter memory), both
/// zero-initialised per clause 4.3. Call [`Self::synthesize_frame`]
/// once per decoded frame, in stream order.
#[derive(Debug, Clone)]
pub struct Synthesizer {
    /// eq (40) past excitation `u(n)`, `n < 0`. `exc_history[EXC_HISTORY
    /// − 1]` is the most recent past sample `u(−1)`; index
    /// `EXC_HISTORY − 1 − j` is `u(−1 − j)`. Zero-init per clause 4.3.
    exc_history: [f32; EXC_HISTORY],
    /// eq (77) synthesis-filter memory `ŝ(n)`, `n < 0`. `syn_mem[M − 1]`
    /// is `ŝ(−1)`; index `M − 1 − j` is `ŝ(−1 − j)`. Zero-init per
    /// clause 4.3.
    syn_mem: [f32; M],
}

impl Default for Synthesizer {
    fn default() -> Self {
        Self::new()
    }
}

impl Synthesizer {
    /// Build a synthesizer with the clause-4.3 all-zero start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            exc_history: [0.0; EXC_HISTORY],
            syn_mem: [0.0; M],
        }
    }

    /// Borrow the eq (40) past-excitation buffer for inspection / tests
    /// (oldest sample first; `[EXC_HISTORY − 1]` is `u(−1)`).
    #[must_use]
    pub fn exc_history(&self) -> &[f32; EXC_HISTORY] {
        &self.exc_history
    }

    /// Borrow the eq (77) 10th-order filter memory for inspection /
    /// tests (`[M − 1]` is `ŝ(−1)`).
    #[must_use]
    pub fn syn_mem(&self) -> &[f32; M] {
        &self.syn_mem
    }

    /// Map a decoder fractional pitch delay `(int_t, frac)` —
    /// `frac ∈ {−1, 0, 1}`, `T = int_t + frac/3` — onto the eq (40)
    /// `(k, t)` form the [`Self::adaptive_sample`] read geometry
    /// consumes.
    ///
    /// The two-branch b30 read in `adaptive_sample` evaluates the
    /// band-limited past excitation at position `n − k + t/3` (write
    /// `b30(j) ≈ sinc(j/3)` and expand either branch: the coefficient
    /// of `u(m)` is `sinc(m − (n − k + t/3))` in both), so the phase
    /// axis runs *against* the delay axis — the effective delay is
    /// `k − t/3`, and `T = int_t + frac/3` folds with a negated
    /// fraction and an upward borrow:
    ///
    /// * `frac = 0`  → `k = int_t`,     `t = 0`
    /// * `frac = −1` → `k = int_t`,     `t = 1`  (`T = int_t − 1/3`)
    /// * `frac = +1` → `k = int_t + 1`, `t = 2`  (`T = (int_t+1) − 2/3`)
    ///
    /// (Round 410: the previous `T = k + t/3` reading was measured
    /// against the conformance corpus to interpolate every fractional
    /// subframe 2/3 of a sample off, decorrelating voiced material.)
    #[inline]
    fn delay_to_kt(int_t: i32, frac: i32) -> (i32, i32) {
        match frac {
            0 => (int_t, 0),
            -1 => (int_t, 1),
            // frac == +1
            _ => (int_t + 1, 2),
        }
    }

    /// Read past excitation `u(idx)` for an index `idx` that the eq (40)
    /// sum may take negative (history) or, within the current subframe,
    /// non-negative (the `cur` slice already built for `n' < n`).
    ///
    /// For `idx < 0` the value comes from [`Self::exc_history`]; for
    /// `idx ≥ 0` it comes from the partially-built current-subframe
    /// excitation `cur`. Indices below `−EXC_HISTORY` cannot occur for
    /// a spec-domain delay (see [`EXC_HISTORY`]); a defensive clamp
    /// returns zero rather than panicking on a hand-built out-of-domain
    /// delay.
    #[inline]
    fn u_at(&self, cur: &[f32; SUBFRAME_SIZE], idx: i32) -> f32 {
        if idx >= 0 {
            // Within the current subframe (only reachable for very small
            // delays, T < 40, where the adaptive vector folds onto the
            // already-computed part of u(n)).
            let i = idx as usize;
            if i < SUBFRAME_SIZE {
                cur[i]
            } else {
                0.0
            }
        } else {
            // idx in {-1, -2, ...}: exc_history[EXC_HISTORY - 1] is u(-1).
            let back = (-idx) as usize; // 1, 2, ...
            if back <= EXC_HISTORY {
                self.exc_history[EXC_HISTORY - back]
            } else {
                0.0
            }
        }
    }

    /// eq (40): interpolate the past excitation `u` at delay `(k, t)`
    /// for one output sample `n`, reading history (`idx < 0`) and the
    /// partially-built current excitation (`0 ≤ idx < n`) through
    /// [`Self::u_at`]. Returns `v(n)`.
    #[inline]
    fn adaptive_sample(&self, cur: &[f32; SUBFRAME_SIZE], n: i32, k: i32, t: i32) -> f32 {
        let b30 = &PITCH_INTERP_FILTER_SYNTHESIS_Q15;
        let mut acc = 0.0f32;
        for i in 0..INTERP_TAPS as i32 {
            // Σ u(n − k − i)·b30(t + 3i)
            let c0 = f32::from(b30[(t + 3 * i) as usize]);
            acc += self.u_at(cur, n - k - i) * c0;
            // Σ u(n − k + 1 + i)·b30(3 − t + 3i)
            let c1 = f32::from(b30[(3 - t + 3 * i) as usize]);
            acc += self.u_at(cur, n - k + 1 + i) * c1;
        }
        acc / B30_Q15_SCALE
    }

    /// eq (77): filter one subframe's excitation `u(n)` through the
    /// 10th-order LP synthesis filter `1/Â(z)`, advancing the filter
    /// memory. `a` holds `â_i` for `i = 1 … 10` in slots `0 … 9`.
    fn synthesis_filter(
        &mut self,
        excitation: &[f32; SUBFRAME_SIZE],
        a: &[f32; M],
    ) -> [f32; SUBFRAME_SIZE] {
        let mut speech = [0.0f32; SUBFRAME_SIZE];
        for n in 0..SUBFRAME_SIZE {
            // ŝ(n) = u(n) − Σ_{i=1}^{10} â_i·ŝ(n − i)
            let mut acc = excitation[n];
            for i in 1..=M {
                let prev = if n >= i {
                    speech[n - i]
                } else {
                    // ŝ(n − i) with n − i < 0 → filter memory.
                    // syn_mem[M - 1] is ŝ(-1); ŝ(-j) is syn_mem[M - j].
                    let j = i - n; // 1, 2, ...
                    self.syn_mem[M - j]
                };
                acc -= a[i - 1] * prev;
            }
            speech[n] = acc;
        }
        speech
    }

    /// Advance both state buffers after a subframe: append the subframe
    /// excitation to the past-excitation buffer and the reconstructed
    /// speech to the synthesis-filter memory.
    fn advance_state(&mut self, excitation: &[f32; SUBFRAME_SIZE], speech: &[f32; SUBFRAME_SIZE]) {
        // Shift the past-excitation buffer left by one subframe and
        // append the 40 new excitation samples as the most-recent past.
        self.exc_history.copy_within(SUBFRAME_SIZE.., 0);
        let tail = EXC_HISTORY - SUBFRAME_SIZE;
        self.exc_history[tail..].copy_from_slice(excitation);

        // The new synthesis memory is the last M reconstructed samples.
        self.syn_mem.copy_from_slice(&speech[SUBFRAME_SIZE - M..]);
    }

    /// Synthesize one subframe: eq (40) `v(n)`, eq (75) `u(n)`,
    /// eq (77) `ŝ(n)`, then advance both state buffers.
    fn synthesize_subframe(
        &mut self,
        int_t: i32,
        frac: i32,
        codevector: &[f32; SUBFRAME_SIZE],
        g_p: f32,
        g_c: f32,
        a: &[f32; M],
    ) -> SynthesizedSubframe {
        let (k, t) = Self::delay_to_kt(int_t, frac);

        // eq (75) excitation built sample-by-sample. For delays T < 40
        // the eq (40) interpolation of v(n) reads u(n') at non-negative
        // n' < n, so u(n) must already be available; we therefore
        // interleave: compute v(n), then u(n), advancing the current
        // excitation slice as we go.
        let mut excitation = [0.0f32; SUBFRAME_SIZE];
        let mut adaptive = [0.0f32; SUBFRAME_SIZE];
        for n in 0..SUBFRAME_SIZE {
            let vn = self.adaptive_sample(&excitation, n as i32, k, t);
            adaptive[n] = vn;
            // eq (75): u(n) = ĝ_p·v(n) + ĝ_c·c(n).
            excitation[n] = g_p * vn + g_c * codevector[n];
        }

        // eq (77): filter the excitation through 1/Â(z).
        let speech = self.synthesis_filter(&excitation, a);

        // Advance cross-subframe state.
        self.advance_state(&excitation, &speech);

        SynthesizedSubframe {
            adaptive,
            excitation,
            speech,
        }
    }

    /// Synthesize one decoded frame into reconstructed speech, running
    /// the §4.1.3 → §3.10 → §4.1.6 chain for both subframes and
    /// advancing the cross-subframe state.
    ///
    /// Consumes a round-282 [`DecodedFrame`]: for each subframe it reads
    /// the §4.1.3 pitch delay, the §4.1.4 (post-eq (48)) codevector, the
    /// §4.1.5 gains `ĝ_p` / `ĝ_c`, and the §4.1.1 LP coefficients.
    #[must_use]
    pub fn synthesize_frame(&mut self, frame: &DecodedFrame) -> SynthesizedFrame {
        let mut subs: [Option<SynthesizedSubframe>; 2] = [None, None];
        for (i, sub) in frame.subframes.iter().enumerate() {
            subs[i] = Some(self.synthesize_subframe(
                sub.pitch.int_t,
                sub.pitch.frac,
                &sub.codevector,
                sub.gains.g_p_hat,
                sub.g_c_hat,
                &sub.lp,
            ));
        }
        let [s0, s1] = subs;
        SynthesizedFrame {
            subframes: [
                s0.expect("subframe 1 synthesized above"),
                s1.expect("subframe 2 synthesized above"),
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decode_chain::FrameDecoder;
    use crate::parameters::Parameters;

    /// A frame whose pitch parity holds (mirrors the decode-chain test
    /// fixture).
    fn parity_ok_params() -> Parameters {
        Parameters {
            l0: 1,
            l1: 17,
            l2: 5,
            l3: 9,
            p1: 0b1100_0000,
            p0: 1,
            c1: 0x0AB3,
            s1: 0b1010,
            ga1: 3,
            gb1: 7,
            p2: 11,
            c2: 0x15C2,
            s2: 0b0101,
            ga2: 6,
            gb2: 12,
        }
    }

    #[test]
    fn new_synthesizer_has_zero_state() {
        let s = Synthesizer::new();
        assert!(s.exc_history().iter().all(|&x| x == 0.0));
        assert!(s.syn_mem().iter().all(|&x| x == 0.0));
        assert_eq!(EXC_HISTORY, 153);
    }

    /// eq (40) with an all-zero history and a delay deep enough that
    /// every interpolation tap reads the (zero) history gives a zero
    /// adaptive vector. With zero gains the synthesized excitation and
    /// speech are zero too.
    #[test]
    fn zero_history_gives_zero_adaptive_vector() {
        let mut s = Synthesizer::new();
        let cv = [0.0f32; SUBFRAME_SIZE];
        let a = [0.0f32; M];
        // int_t = 80 ≥ 40 + 10 → every eq (40) access stays in history.
        let sub = s.synthesize_subframe(80, 0, &cv, 0.7, 0.0, &a);
        assert!(sub.adaptive.iter().all(|&x| x == 0.0));
        assert!(sub.excitation.iter().all(|&x| x == 0.0));
    }

    /// `delay_to_kt` realises the eq (40) read-geometry fold: the
    /// effective delay of the `(k, t)` pair is `k − t/3` (see the
    /// helper's doc comment), so the transmitted fraction negates and
    /// the `frac = +1` case borrows upward.
    #[test]
    fn delay_to_kt_fraction_convention() {
        assert_eq!(Synthesizer::delay_to_kt(50, 0), (50, 0));
        assert_eq!(Synthesizer::delay_to_kt(50, -1), (50, 1));
        assert_eq!(Synthesizer::delay_to_kt(50, 1), (51, 2));
    }

    /// eq (77) with all-zero LP coefficients is the identity
    /// `ŝ(n) = u(n)`.
    #[test]
    fn synthesis_filter_identity_with_zero_lp() {
        let mut s = Synthesizer::new();
        let mut exc = [0.0f32; SUBFRAME_SIZE];
        for (n, e) in exc.iter_mut().enumerate() {
            *e = n as f32;
        }
        let a = [0.0f32; M];
        let sp = s.synthesis_filter(&exc, &a);
        assert_eq!(sp, exc);
    }

    /// eq (77) single-tap recurrence pins the spec sign and memory
    /// threading: with `â_1 = 0.5` and the rest zero,
    /// `ŝ(n) = u(n) − 0.5·ŝ(n − 1)`, `ŝ(−1) = 0`.
    #[test]
    fn synthesis_filter_single_tap_recurrence() {
        let mut s = Synthesizer::new();
        let mut exc = [0.0f32; SUBFRAME_SIZE];
        exc[0] = 1.0;
        exc[1] = 1.0;
        let mut a = [0.0f32; M];
        a[0] = 0.5;
        let sp = s.synthesis_filter(&exc, &a);
        // ŝ(0) = 1 − 0.5·0     = 1
        // ŝ(1) = 1 − 0.5·1     = 0.5
        // ŝ(2) = 0 − 0.5·0.5   = -0.25
        assert!((sp[0] - 1.0).abs() < 1e-6);
        assert!((sp[1] - 0.5).abs() < 1e-6);
        assert!((sp[2] - (-0.25)).abs() < 1e-6);
    }

    /// eq (75): excitation is exactly `ĝ_p·v(n) + ĝ_c·c(n)`. With zero
    /// history `v(n) = 0`, so `u(n) = ĝ_c·c(n)`.
    #[test]
    fn excitation_is_eq75_combination() {
        let mut s = Synthesizer::new();
        let mut cv = [0.0f32; SUBFRAME_SIZE];
        cv[0] = 1.0;
        cv[7] = -1.0;
        let g_c = 2.5;
        let a = [0.0f32; M];
        // A delay ≥ 40 + 10 keeps every eq (40) access in the
        // (zero) history, so v(n) = 0 across the whole subframe.
        let sub = s.synthesize_subframe(80, 0, &cv, 0.7, g_c, &a);
        // Zero history → v(n) = 0 → u(n) = g_c·c(n).
        for (n, &c) in cv.iter().enumerate() {
            assert!((sub.excitation[n] - g_c * c).abs() < 1e-5);
        }
        // Zero LP → ŝ(n) = u(n).
        assert_eq!(sub.speech, sub.excitation);
    }

    /// State advance threads the excitation into the past buffer: after
    /// one subframe the most-recent 40 history samples are that
    /// subframe's excitation (so the next subframe's adaptive vector
    /// can reference it).
    #[test]
    fn state_advances_past_excitation_buffer() {
        let mut s = Synthesizer::new();
        let mut cv = [0.0f32; SUBFRAME_SIZE];
        cv[0] = 1.0;
        let a = [0.0f32; M];
        let sub = s.synthesize_subframe(40, 0, &cv, 0.0, 1.0, &a);
        let hist = s.exc_history();
        // The last 40 history samples == this subframe's excitation.
        for n in 0..SUBFRAME_SIZE {
            assert!((hist[EXC_HISTORY - SUBFRAME_SIZE + n] - sub.excitation[n]).abs() < 1e-6);
        }
        // syn_mem holds the last M reconstructed samples.
        for j in 0..M {
            assert!((s.syn_mem()[j] - sub.speech[SUBFRAME_SIZE - M + j]).abs() < 1e-6);
        }
    }

    /// A short-delay (`T < 40`) subframe makes eq (40) fold onto the
    /// already-built current excitation: `v(n)` for `n ≥ k` reads
    /// `u(n − k)` which lives in the current subframe. The synthesis
    /// stays finite and the adaptive vector is non-zero once the
    /// excitation builds up.
    #[test]
    fn short_delay_folds_onto_current_excitation() {
        let mut s = Synthesizer::new();
        let mut cv = [0.0f32; SUBFRAME_SIZE];
        cv[0] = 1.0;
        let a = [0.0f32; M];
        // int_t = 20 (< 40), frac = 0 → k = 20, t = 0.
        let sub = s.synthesize_subframe(20, 0, &cv, 0.8, 1.0, &a);
        assert!(sub.adaptive.iter().all(|x| x.is_finite()));
        assert!(sub.excitation.iter().all(|x| x.is_finite()));
        assert!(sub.speech.iter().all(|x| x.is_finite()));
        // Some adaptive-vector sample past n = 20 must be non-zero once
        // u(0) = g_c·c(0) > 0 feeds back through the interpolation.
        assert!(sub.adaptive[20..].iter().any(|&x| x != 0.0));
    }

    /// End-to-end: the decode chain → synthesizer produces finite
    /// reconstructed speech for both subframes of a real-shaped frame.
    #[test]
    fn decode_chain_to_synthesis_is_finite() {
        let mut chain = FrameDecoder::new();
        let mut synth = Synthesizer::new();
        let frame = chain
            .decode_parameters(&parity_ok_params())
            .expect("in-domain");
        let out = synth.synthesize_frame(&frame);
        for sub in &out.subframes {
            assert!(sub.adaptive.iter().all(|x| x.is_finite()));
            assert!(sub.excitation.iter().all(|x| x.is_finite()));
            assert!(sub.speech.iter().all(|x| x.is_finite()));
        }
        assert_eq!(out.speech().len(), 80);
    }

    /// Two synthesizers over the same decoded-frame sequence stay in
    /// lockstep — all state is owned, no hidden globals.
    #[test]
    fn synthesizer_is_deterministic() {
        let mut chain_a = FrameDecoder::new();
        let mut chain_b = FrameDecoder::new();
        let mut a = Synthesizer::new();
        let mut b = Synthesizer::new();
        for _ in 0..4 {
            let fa = chain_a.decode_parameters(&parity_ok_params()).unwrap();
            let fb = chain_b.decode_parameters(&parity_ok_params()).unwrap();
            assert_eq!(a.synthesize_frame(&fa), b.synthesize_frame(&fb));
        }
    }
}
