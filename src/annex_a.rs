//! **Annex A** — the reduced-complexity 8 kbit/s CS-ACELP codec's
//! decoder-side deltas (§A.4), bit-stream interoperable with the main
//! body.
//!
//! Annex A (Recommendation §A.1–A.5, transcribed from the in-repo EPUB
//! `back01.xhtml` prose + equation rasters `bm1eq11 … bm1eq15`) leaves
//! the whole §4.1 parameter decode and §4.1.6 synthesis identical to
//! the main body — "the only change in the decoder is in the
//! postfilter" (§A.4). The §A.4.2 adaptive postfilter differs from
//! clause 4.2 in five pinned ways:
//!
//! * **§A.4.2.1 long-term postfilter** — eq (A.11) is the eq (78)
//!   harmonic filter with **integer delays only**: `T` is searched in
//!   `[T_cl − 3, T_cl + 3]` where `T_cl` is the integer part of the
//!   *current subframe's* transmitted pitch delay, bounded
//!   `T_cl ≤ 140` (the main body anchors on `int(T_1)` and refines to
//!   1/8 resolution).
//! * **§A.4.2.2 short-term postfilter** — eq (A.12) is the eq (84)
//!   pole/zero pair `Â(z/γ_n)/Â(z/γ_d)` (γ_n = 0.55, γ_d = 0.7)
//!   **without** the `1/g_f` gain normalisation.
//! * **§A.4.2.3 tilt compensation** — eq (A.13) `H_t(z) = 1 + γ_t·k′_1·z⁻¹`
//!   without the `1/g_t` normalisation; eq (A.14) computes `k′_1 =
//!   −r_h(1)/r_h(0)` from the **length-22** truncated impulse response
//!   of `Â(z/γ_n)/Â(z/γ_d)` (`r_h(i) = Σ_{j=0}^{21−i} h_f(j)·h_f(j+i)`);
//!   `γ_t = 0.8` if `k′_1 < 0`, else `γ_t = 0` (the main body uses
//!   0.9 / 0.2).
//! * **Cascade order** — "the compensation filtering is performed
//!   before synthesis filtering through `1/Â(z/γ_d)`" (§A.4.2): the
//!   per-sample chain is numerator `Â(z/γ_n)` → tilt `H_t` → synthesis
//!   `1/Â(z/γ_d)` (the main body runs the full `H_f` then `H_t`).
//! * **§A.4.2.4 adaptive gain control** — eq (A.15) is the **energy**
//!   ratio `G = √(Σ ŝ²(n) / Σ sf²(n))` (the main body eq (88) uses
//!   absolute sums), smoothed as `g(n) = 0.9·g(n−1) + 0.1·G`
//!   (main body: 0.85 / 0.15).
//!
//! §A.4.2.5 high-pass + ×2 upscaling is unchanged
//! ([`crate::post_process::OutputHighPass`] is reused).
//!
//! [`AnnexADecoder`] wraps the main-body [`FrameDecoder`]'s §4.1
//! parameter chain + §4.1.6 synthesis and swaps in the
//! [`AnnexAPostfilterCascade`]. It is validated black-box against the
//! staged `g729a` conformance corpus
//! (`docs/audio/g729/conformance/g729a/`, `tests/annex_a_conformance.rs`):
//! the same first-subframe agreement band and bounded-RMS-ratio guard
//! as the base-codec `.PST` harness (the §3.9.1 gain-predictor
//! fixed-point saturation gap applies identically here).
//!
//! Encoder-side §A.3 deltas (γ = 0.75 quantized-LP weighting, fast
//! open-loop / closed-loop pitch, depth-first ACELP search) are not
//! implemented in this module; the §A.3.8.1 depth-first tree search's
//! exact pulse-combination schedule is prose-unpinned (Annex A defers
//! it to the barred reference C).

use crate::decode_chain::{DecodedFrame, FrameDecodeError, FrameDecoder};
use crate::fixed_codebook::SUBFRAME_SIZE;
use crate::long_term_postfilter::{ENABLE_THRESHOLD, GAMMA_P, MAX_DELAY};
use crate::lp_synthesis::SynthesizedFrame;
use crate::parameters::Parameters;
use crate::post_process::OutputHighPass;
use crate::serial::parse_frame;
use crate::short_term_postfilter::{GAMMA_D, GAMMA_N};
use crate::tables::M;

/// §A.4.2.1 delay-search half range: `T ∈ [T_cl − 3, T_cl + 3]`.
pub const DELAY_HALF_RANGE: usize = 3;

/// §A.4.2.1 bound on the search anchor: `T_cl ≤ 140` (so the search
/// never reaches past delay 143, the eq (40) maximum).
pub const TCL_MAX: usize = 140;

/// §A.4.2.3 / eq (A.14) truncated impulse-response length of
/// `Â(z/γ_n)/Â(z/γ_d)` (`j = 0 … 21`).
pub const TILT_IMPULSE_LEN: usize = 22;

/// §A.4.2.3 tilt weight when `k′_1 < 0` (`γ_t = 0.8`; zero otherwise).
pub const TILT_GAMMA_T_NEG: f32 = 0.8;

/// §A.4.2.4 gain-smoothing weight on the previous gain (`0.9`).
pub const AGC_PREV_WEIGHT: f32 = 0.9;

/// §A.4.2.4 gain-smoothing weight on the eq (A.15) target (`0.1`).
pub const AGC_TARGET_WEIGHT: f32 = 0.1;

/// History window long enough to reach the maximum searched delay
/// (`T_cl + 3 ≤ 143 < MAX_DELAY`).
const HIST_LEN: usize = MAX_DELAY;

/// The §A.4.2.1 integer-delay decision for one subframe.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnnexALongTermDecision {
    /// The winning integer delay `T`.
    pub delay: usize,
    /// eq (83) long-term gain `g_l`, bounded `[0, 1]`; exactly `0.0`
    /// when the eq (82) 3 dB disable test fires.
    pub gain: f32,
}

/// §A.4.2.1 integer-delay long-term (harmonic) postfilter.
#[derive(Debug, Clone)]
struct AnnexALongTerm {
    /// Reconstructed-speech history (`s_hist[0]` = most recent
    /// `ŝ(n−1)`), zero-init per clause 4.3.
    s_hist: [f32; HIST_LEN],
    /// eq (79) residual history, same layout.
    r_hist: [f32; HIST_LEN],
}

impl AnnexALongTerm {
    fn new() -> Self {
        Self {
            s_hist: [0.0; HIST_LEN],
            r_hist: [0.0; HIST_LEN],
        }
    }

    /// eq (79): the residual of the current subframe through the
    /// short-term numerator `Â(z/γ_n)`, using the carried speech
    /// history for the first ten samples.
    fn residual(&self, s: &[f32; SUBFRAME_SIZE], a: &[f32; M]) -> [f32; SUBFRAME_SIZE] {
        let mut wn = [0.0f32; M];
        let mut g = 1.0f32;
        for (w, &ai) in wn.iter_mut().zip(a.iter()) {
            g *= GAMMA_N;
            *w = g * ai;
        }
        let mut r = [0.0f32; SUBFRAME_SIZE];
        for (n, out) in r.iter_mut().enumerate() {
            let mut acc = s[n];
            for (i, &w) in wn.iter().enumerate() {
                let lag = i + 1;
                let sample = if n >= lag {
                    s[n - lag]
                } else {
                    self.s_hist[lag - n - 1]
                };
                acc += w * sample;
            }
            *out = acc;
        }
        r
    }

    /// Residual sample at `n − k` (reaching into the carried history
    /// when `n < k`).
    fn residual_at(&self, r: &[f32; SUBFRAME_SIZE], n: usize, k: usize) -> f32 {
        if n >= k {
            r[n - k]
        } else {
            self.r_hist[k - n - 1]
        }
    }

    /// Speech sample at `n − k`.
    fn speech_at(&self, s: &[f32; SUBFRAME_SIZE], n: usize, k: usize) -> f32 {
        if n >= k {
            s[n - k]
        } else {
            self.s_hist[k - n - 1]
        }
    }

    /// Filter one subframe: eq (80) integer search over
    /// `[T_cl − 3, T_cl + 3]`, the eq (82) disable test and eq (83)
    /// bounded gain, then eq (A.11) applied to `ŝ(n)` at the integer
    /// delay. Advances the carried histories.
    fn filter_subframe(
        &mut self,
        s: &[f32; SUBFRAME_SIZE],
        a: &[f32; M],
        t_cl: usize,
    ) -> ([f32; SUBFRAME_SIZE], AnnexALongTermDecision) {
        let r = self.residual(s, a);
        let t_cl = t_cl.min(TCL_MAX);

        // eq (80) first pass: the correlation-maximising integer delay.
        let lo = t_cl.saturating_sub(DELAY_HALF_RANGE).max(1);
        let hi = (t_cl + DELAY_HALF_RANGE).min(HIST_LEN - 1);
        let mut t0 = lo;
        let mut best_corr = f32::NEG_INFINITY;
        for k in lo..=hi {
            let mut corr = 0.0f32;
            for (n, &rn) in r.iter().enumerate() {
                corr += rn * self.residual_at(&r, n, k);
            }
            if corr > best_corr {
                best_corr = corr;
                t0 = k;
            }
        }

        // eq (82)/(83): disable test + bounded gain at the winner.
        let energy: f32 = r.iter().map(|v| v * v).sum();
        let mut den = 0.0f32;
        for n in 0..SUBFRAME_SIZE {
            let rt = self.residual_at(&r, n, t0);
            den += rt * rt;
        }
        let enabled = best_corr > 0.0
            && den > 0.0
            && energy > 0.0
            && best_corr * best_corr >= ENABLE_THRESHOLD * energy * den;
        let gain = if enabled {
            (best_corr / den).clamp(0.0, 1.0)
        } else {
            0.0
        };

        // eq (A.11): out(n) = (ŝ(n) + γ_p·g_l·ŝ(n−T)) / (1 + γ_p·g_l).
        let scale = GAMMA_P * gain;
        let inv = 1.0 / (1.0 + scale);
        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &sn) in s.iter().enumerate() {
            let delayed = if gain == 0.0 {
                0.0
            } else {
                self.speech_at(s, n, t0)
            };
            out[n] = inv * (sn + scale * delayed);
        }

        // Advance histories (most-recent at index 0).
        self.s_hist
            .copy_within(0..HIST_LEN - SUBFRAME_SIZE, SUBFRAME_SIZE);
        self.r_hist
            .copy_within(0..HIST_LEN - SUBFRAME_SIZE, SUBFRAME_SIZE);
        for j in 0..SUBFRAME_SIZE {
            self.s_hist[j] = s[SUBFRAME_SIZE - 1 - j];
            self.r_hist[j] = r[SUBFRAME_SIZE - 1 - j];
        }

        (out, AnnexALongTermDecision { delay: t0, gain })
    }
}

/// §A.4.2.2 + §A.4.2.3 short-term postfilter and tilt compensation in
/// the Annex A order (numerator → tilt → synthesis), with no `1/g_f` /
/// `1/g_t` normalisation.
#[derive(Debug, Clone)]
struct AnnexAShortTermTilt {
    /// Input history for the numerator `Â(z/γ_n)` (most recent first).
    x_hist: [f32; M],
    /// The tilt filter's `u(n−1)` memory.
    u_prev: f32,
    /// Output history for the synthesis `1/Â(z/γ_d)`.
    y_hist: [f32; M],
}

impl AnnexAShortTermTilt {
    fn new() -> Self {
        Self {
            x_hist: [0.0; M],
            u_prev: 0.0,
            y_hist: [0.0; M],
        }
    }

    /// eq (A.14): `k′_1` from the length-22 truncated impulse response
    /// of `Â(z/γ_n)/Â(z/γ_d)`.
    fn tilt_k1(wn: &[f32; M], wd: &[f32; M]) -> f32 {
        let mut h = [0.0f32; TILT_IMPULSE_LEN];
        for n in 0..TILT_IMPULSE_LEN {
            let mut acc = if n == 0 { 1.0 } else { 0.0 };
            if (1..=M).contains(&n) {
                acc += wn[n - 1];
            }
            for (i, &w) in wd.iter().enumerate() {
                let lag = i + 1;
                if n >= lag {
                    acc -= w * h[n - lag];
                }
            }
            h[n] = acc;
        }
        let mut rh0 = 0.0f32;
        let mut rh1 = 0.0f32;
        for j in 0..TILT_IMPULSE_LEN {
            rh0 += h[j] * h[j];
            if j + 1 < TILT_IMPULSE_LEN {
                rh1 += h[j] * h[j + 1];
            }
        }
        if rh0 > 0.0 {
            -rh1 / rh0
        } else {
            0.0
        }
    }

    /// Run one subframe through numerator → tilt → synthesis.
    fn filter_subframe(&mut self, x: &[f32; SUBFRAME_SIZE], a: &[f32; M]) -> [f32; SUBFRAME_SIZE] {
        let mut wn = [0.0f32; M];
        let mut wd = [0.0f32; M];
        let (mut gn, mut gd) = (1.0f32, 1.0f32);
        for i in 0..M {
            gn *= GAMMA_N;
            gd *= GAMMA_D;
            wn[i] = gn * a[i];
            wd[i] = gd * a[i];
        }
        // §A.4.2.3: γ_t = 0.8 iff k′_1 < 0, else the tilt is disabled.
        let k1 = Self::tilt_k1(&wn, &wd);
        let b1 = if k1 < 0.0 { TILT_GAMMA_T_NEG * k1 } else { 0.0 };

        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &xn) in x.iter().enumerate() {
            // Numerator Â(z/γ_n): u(n) = x(n) + Σ wn_i·x(n−i).
            let mut u = xn;
            for (&w, &xh) in wn.iter().zip(self.x_hist.iter()) {
                u += w * xh;
            }
            // eq (A.13) tilt: t(n) = u(n) + γ_t·k′_1·u(n−1).
            let t = u + b1 * self.u_prev;
            self.u_prev = u;
            // Synthesis 1/Â(z/γ_d): y(n) = t(n) − Σ wd_i·y(n−i).
            let mut y = t;
            for (&w, &yh) in wd.iter().zip(self.y_hist.iter()) {
                y -= w * yh;
            }
            for i in (1..M).rev() {
                self.x_hist[i] = self.x_hist[i - 1];
                self.y_hist[i] = self.y_hist[i - 1];
            }
            self.x_hist[0] = xn;
            self.y_hist[0] = y;
            out[n] = y;
        }
        out
    }
}

/// §A.4.2.4 adaptive gain control on the eq (A.15) energy ratio.
#[derive(Debug, Clone)]
struct AnnexAAgc {
    g: f32,
}

impl AnnexAAgc {
    fn new() -> Self {
        // Table 9 start-up gain g(−1) = 1.0 (clause A.4.3 defers to 4.3).
        Self { g: 1.0 }
    }

    fn process_subframe(
        &mut self,
        s_hat: &[f32; SUBFRAME_SIZE],
        sf: &[f32; SUBFRAME_SIZE],
    ) -> [f32; SUBFRAME_SIZE] {
        let num: f32 = s_hat.iter().map(|v| v * v).sum();
        let den: f32 = sf.iter().map(|v| v * v).sum();
        // eq (A.15): G = √(Σ ŝ² / Σ sf²); hold the previous gain on a
        // silent postfiltered subframe.
        let g_target = if den > 0.0 {
            (num / den).sqrt()
        } else {
            self.g
        };
        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &sfn) in sf.iter().enumerate() {
            self.g = AGC_PREV_WEIGHT * self.g + AGC_TARGET_WEIGHT * g_target;
            out[n] = self.g * sfn;
        }
        out
    }
}

/// The §A.4.2 post-processing cascade (long-term integer-delay
/// postfilter → short-term numerator → tilt → synthesis → adaptive
/// gain control → output high-pass + ×2 upscaling).
#[derive(Debug, Clone)]
pub struct AnnexAPostfilterCascade {
    long_term: AnnexALongTerm,
    short_term: AnnexAShortTermTilt,
    agc: AnnexAAgc,
    high_pass: OutputHighPass,
}

impl Default for AnnexAPostfilterCascade {
    fn default() -> Self {
        Self::new()
    }
}

impl AnnexAPostfilterCascade {
    /// Build the cascade with the clause-4.3 start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            long_term: AnnexALongTerm::new(),
            short_term: AnnexAShortTermTilt::new(),
            agc: AnnexAAgc::new(),
            high_pass: OutputHighPass::new(),
        }
    }

    /// Run the §A.4.2 cascade on one decoded frame, advancing every
    /// stage's cross-subframe state. Returns the 80 post-processed,
    /// ×2-upscaled output samples and the two per-subframe long-term
    /// decisions.
    #[must_use]
    pub fn process_frame(
        &mut self,
        synth: &SynthesizedFrame,
        frame: &DecodedFrame,
    ) -> ([f32; 2 * SUBFRAME_SIZE], [AnnexALongTermDecision; 2]) {
        let mut out = [0.0f32; 2 * SUBFRAME_SIZE];
        let mut decisions = [AnnexALongTermDecision {
            delay: 0,
            gain: 0.0,
        }; 2];
        for i in 0..2 {
            let s_hat = &synth.subframes[i].speech;
            let a: &[f32; M] = &frame.subframes[i].lp;
            // §A.4.2.1: T_cl = the *current subframe's* transmitted
            // integer delay, bounded ≤ 140.
            let t_cl = frame.subframes[i].pitch.int_t.max(0) as usize;

            let (hp, decision) = self.long_term.filter_subframe(s_hat, a, t_cl);
            let sf = self.short_term.filter_subframe(&hp, a);
            let mut scaled = self.agc.process_subframe(s_hat, &sf);
            self.high_pass.filter_in_place(&mut scaled);

            out[i * SUBFRAME_SIZE..(i + 1) * SUBFRAME_SIZE].copy_from_slice(&scaled);
            decisions[i] = decision;
        }
        (out, decisions)
    }
}

/// The Annex A decoder: the main-body §4.1 parameter chain + §4.1.6
/// synthesis (unchanged per §A.4.1), post-processed by the §A.4.2
/// reduced-complexity cascade.
#[derive(Debug, Clone, Default)]
pub struct AnnexADecoder {
    inner: FrameDecoder,
    cascade: AnnexAPostfilterCascade,
}

impl AnnexADecoder {
    /// Build a decoder in the clause-4.3 start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            inner: FrameDecoder::new(),
            cascade: AnnexAPostfilterCascade::new(),
        }
    }

    /// Borrow the wrapped main-body decoder (parameter chain state).
    #[must_use]
    pub fn inner(&self) -> &FrameDecoder {
        &self.inner
    }

    /// Decode one 164-byte ITU serial frame to the §A.4.2
    /// post-processed output (80 ×2-upscaled samples).
    ///
    /// # Errors
    ///
    /// As [`FrameDecoder::decode_serial_frame`]; an erasure sentinel
    /// returns [`FrameDecodeError::Erased`] (the §A.4.4 concealment —
    /// clause 4.4 without voicing detection — is not wired here).
    pub fn decode_serial_frame_to_postfiltered(
        &mut self,
        frame_bytes: &[u8],
    ) -> Result<[f32; 2 * SUBFRAME_SIZE], FrameDecodeError> {
        let kind = parse_frame(frame_bytes)?;
        let params = crate::parameters::unpack_parameters(&kind)?;
        self.decode_parameters_to_postfiltered(&params)
    }

    /// Decode one frame's unpacked Table-8 codewords to the §A.4.2
    /// post-processed output.
    ///
    /// # Errors
    ///
    /// As [`FrameDecoder::decode_parameters`].
    pub fn decode_parameters_to_postfiltered(
        &mut self,
        params: &Parameters,
    ) -> Result<[f32; 2 * SUBFRAME_SIZE], FrameDecodeError> {
        let (frame, speech) = self.inner.decode_parameters_with_speech(params)?;
        let (out, _) = self.cascade.process_frame(&speech, &frame);
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A flat (all-zero `â_i`) short-term/tilt stage is the identity:
    /// `h_f = δ`, `k′_1 = 0` → tilt disabled, numerator and synthesis
    /// both degenerate to 1.
    #[test]
    fn flat_short_term_tilt_is_identity() {
        let mut st = AnnexAShortTermTilt::new();
        let a = [0.0f32; M];
        let x: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let out = st.filter_subframe(&x, &a);
        for n in 0..SUBFRAME_SIZE {
            assert!((out[n] - x[n]).abs() < 1e-6, "n={n}");
        }
    }

    /// eq (A.14) on a flat filter: `k′_1 = 0` (impulse response is the
    /// unit impulse, so `r_h(1) = 0`).
    #[test]
    fn tilt_k1_flat_is_zero() {
        let wn = [0.0f32; M];
        let wd = [0.0f32; M];
        assert!(AnnexAShortTermTilt::tilt_k1(&wn, &wd).abs() < 1e-9);
    }

    /// The §A.4.2.1 integer search locks onto the true period of a
    /// synthetic periodic residual when `T_cl` points near it, and the
    /// eq (83) gain for a perfectly periodic signal is close to 1.
    #[test]
    fn long_term_locks_on_period() {
        let mut lt = AnnexALongTerm::new();
        let a = [0.0f32; M]; // residual = speech
        let period = 40usize;
        // Prime the history with two subframes of a period-40 pulse
        // train, then measure on a third.
        let s: [f32; SUBFRAME_SIZE] =
            core::array::from_fn(|n| if n % period == 0 { 100.0 } else { 1.0 });
        let _ = lt.filter_subframe(&s, &a, period);
        let _ = lt.filter_subframe(&s, &a, period);
        let (_, decision) = lt.filter_subframe(&s, &a, period);
        assert_eq!(decision.delay, period);
        assert!(
            (decision.gain - 1.0).abs() < 1e-3,
            "gain = {}",
            decision.gain
        );
    }

    /// A `T_cl` above the §A.4.2.1 bound is clamped to 140 (search
    /// window `[137, 143]`).
    #[test]
    fn long_term_clamps_tcl() {
        let mut lt = AnnexALongTerm::new();
        let a = [0.0f32; M];
        let s = [1.0f32; SUBFRAME_SIZE];
        let (_, decision) = lt.filter_subframe(&s, &a, 500);
        assert!((137..=143).contains(&decision.delay));
    }

    /// eq (A.15) AGC drives a scaled copy back toward the reference
    /// energy: filtering `2·ŝ` against reference `ŝ` converges to gain
    /// ≈ 0.5.
    #[test]
    fn agc_energy_ratio_converges() {
        let mut agc = AnnexAAgc::new();
        let s_hat: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| ((n * 7) % 13) as f32 - 6.0);
        let sf: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| 2.0 * s_hat[n]);
        let mut last = [0.0f32; SUBFRAME_SIZE];
        for _ in 0..8 {
            last = agc.process_subframe(&s_hat, &sf);
        }
        // After several subframes g(n) ≈ G = 0.5, so out ≈ ŝ.
        for n in 0..SUBFRAME_SIZE {
            assert!(
                (last[n] - s_hat[n]).abs() < 0.05 * (s_hat[n].abs() + 1.0),
                "n={n}: {} vs {}",
                last[n],
                s_hat[n]
            );
        }
    }

    /// End-to-end: an [`AnnexADecoder`] decodes a hand-built valid
    /// frame to finite post-processed output from the start-up state.
    #[test]
    fn annex_a_decoder_decodes_finite() {
        let params = Parameters {
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
        };
        let mut dec = AnnexADecoder::new();
        let out = dec
            .decode_parameters_to_postfiltered(&params)
            .expect("valid frame");
        assert!(out.iter().all(|v| v.is_finite()));
        assert!(out.iter().any(|&v| v != 0.0));
    }
}
