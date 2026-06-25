//! §4.2.1 **long-term (harmonic) postfilter** `H_p(z)` — the first
//! filter of the decoder's §4.2 post-processing adaptive postfilter
//! cascade.
//!
//! The §4.2 cascade is (clause 4.2): **long-term postfilter `H_p(z)`
//! (§4.2.1)** → short-term postfilter `H_f(z)` (§4.2.2,
//! [`crate::short_term_postfilter`]) → tilt compensation `H_t(z)`
//! (§4.2.3, [`crate::tilt_compensation`]) → adaptive gain control
//! (§4.2.4, [`crate::adaptive_gain_control`]) → output high-pass + ×2
//! upscaling (§4.2.5, [`crate::post_process`]). The four trailing stages
//! were wired in earlier rounds; this module wires the **head** of the
//! cascade, the harmonic postfilter that re-emphasises the pitch
//! periodicity of the reconstructed speech.
//!
//! ## Spec source — clause 4.2.1, equations (78)–(83) (06/2012 Rec.)
//!
//! Clause 4.2.1 (transcribed from the EPUB prose): "The long-term
//! postfilter is given by [eq (78)] where `T` is the pitch delay, and
//! `g_l` is the gain coefficient. Note that `g_l` is bounded by 1, and
//! it is set to zero if the long-term prediction gain is less than 3 dB.
//! The factor `γ_p` controls the amount of long-term postfiltering and
//! has the value of `γ_p = 0.5`. The long-term delay and gain are
//! computed from the residual signal `r̂(n)` obtained by filtering the
//! speech `ŝ(n)` through `Â(z/γ_n)`, which is the numerator of the
//! short-term postfilter (see clause 4.2.2)."
//!
//! Equation (78) (raster `images/eq78.jpg`) is the transfer function
//!
//! ```text
//!              1
//! H_p(z) = ──────────·(1 + γ_p·g_l·z⁻ᵀ)        γ_p = 0.5
//!           1 + γ_p·g_l
//! ```
//!
//! Equation (79) (raster `images/eq79.jpg`) is the residual
//!
//! ```text
//!           10
//! r̂(n) = ŝ(n) + Σ γ_n^i·â_i·ŝ(n−i)            γ_n = 0.55
//!           i=1
//! ```
//!
//! — i.e. the numerator `Â(z/γ_n)` of the §4.2.2 short-term postfilter
//! applied to the reconstructed speech (`γ_n` = [`GAMMA_N`], the
//! short-term numerator weight, clause 4.2.2). This is exactly
//! [`ShortTermPostfilter::weighted_num`] driven without the §4.2.2
//! denominator / gain normalisation, and this module reuses that helper
//! so the two stages share one definition of `γ_n^i·â_i`.
//!
//! The long-term delay is a **two-pass** procedure (clause 4.2.1):
//!
//! 1. **Integer pass** — pick the best integer `T_0` in the range
//!    `[int(T_1) − 1, int(T_1) + 1]`, where `int(T_1)` is the integer
//!    part of the transmitted first-subframe pitch delay `T_1`. "The
//!    best integer delay is the one that maximizes the correlation"
//!    (eq (80), raster `images/eq80.jpg`):
//!
//!    ```text
//!             39
//!    R(k) =   Σ  r̂(n)·r̂(n−k)
//!            n=0
//!    ```
//!
//! 2. **Fractional pass** — refine `T_0` to a 1/8-resolution fractional
//!    delay `T` maximising the pseudo-normalised correlation `R′(k)`
//!    (eq (81)), using an interpolation filter of length 33 then 129
//!    (the staged `tab_hup_s` / `tab_hup_l` tables).
//!
//! Once the delay is fixed, the long-term prediction gain test (eq (82),
//! raster `images/eq82.jpg`) disables the filter (`g_l = 0`) if
//!
//! ```text
//!     R′(T)²
//!   ────────── < 0.5
//!    Σ r̂(n)²
//! ```
//!
//! and otherwise the gain is (eq (83), raster `images/eq83.jpg`)
//!
//! ```text
//!         Σ r̂(n)·r̂_k(n)
//! g_l = ─────────────────         bounded 0 ≤ g_l ≤ 1.0
//!        Σ r̂_k(n)·r̂_k(n)
//! ```
//!
//! where `r̂_k(n)` is the residual delayed by the chosen delay (`r̂(n−T)`
//! for an integer delay).
//!
//! ## Scope of this module — the **two-pass** postfilter (eqs (78)–(83))
//!
//! This module realises the full clause-4.2.1 two-pass long-term delay
//! search: the residual (eq (79)), the integer delay search (eq (80)
//! over the three candidates `int(T_1) − 1 … int(T_1) + 1`), the
//! **1/8-resolution fractional second pass** (eq (81)), the long-term
//! prediction-gain disable test (eq (82)), the bounded gain (eq (83)),
//! and the eq (78) harmonic filter at the chosen fractional delay.
//!
//! The fractional second pass (clause 4.2.1) refines the integer anchor
//! `T_0` to `T = T_0 + frac/8` (`frac ∈ {0 … 7}`) by maximising the
//! pseudo-normalised correlation `R(T)² / E_T`. The interpolation uses
//! the staged `tab_hup_s` (length-33, 28 stored coefficients = 7 phases
//! × 4 taps) filter first; the chosen non-integer fraction is then
//! re-evaluated with the longer `tab_hup_l` (length-129, 112 = 7 phases
//! × 16 taps) filter, and the long-filter fraction replaces the short
//! one only if it increases `R(T)²/E_T` ("the new signal replaces the
//! previous one only if the longer filter increases the value of
//! `R′(T)`"). The per-phase tap layout of the two tables — previously a
//! recorded docs-gap — is the 7-phase polyphase decomposition of the
//! staged CSVs (see [`crate::tables::postfilter_interp_short`] /
//! [`crate::tables::postfilter_interp_long`]); phase `p` and phase
//! `8 − p` are mirror images.
//!
//! ## State (clause 4.3 init)
//!
//! Per clause 4.3, "all static encoder and decoder variables should be
//! initialized to zero, except the variables listed in Table 9". Neither
//! the reconstructed-speech history `ŝ(n−i)` (needed by eq (79) and by
//! the eq (78) `z⁻ᵀ` delay) nor the residual history `r̂(n−k)` (needed by
//! the eq (80) correlation) appears in Table 9, so both start zeroed and
//! carry across subframes.

use crate::fixed_codebook::SUBFRAME_SIZE;
use crate::short_term_postfilter::ShortTermPostfilter;
use crate::tables::{
    postfilter_interp_long, postfilter_interp_short, M, POSTFILTER_FRAC_RES,
    POSTFILTER_INTERP_LONG_TAPS, POSTFILTER_INTERP_SHORT_TAPS,
};

/// §4.2.1 long-term-postfilter weight factor `γ_p = 0.5` (clause 4.2.1:
/// "The factor `γ_p` controls the amount of long-term postfiltering and
/// has the value of `γ_p = 0.5`").
pub const GAMMA_P: f32 = 0.5;

/// eq (82) long-term-prediction-gain disable threshold. The filter is
/// disabled (`g_l = 0`) when the squared normalised correlation
/// `R′(T)² / Σ r̂(n)²` falls below this value (the spec's "less than
/// 3 dB" long-term prediction gain, expressed as the `0.5` ratio in
/// eq (82)).
pub const ENABLE_THRESHOLD: f32 = 0.5;

/// Largest integer pitch delay the search can reach. The transmitted
/// first-subframe delay `T_1` has integer part in `[19, 143]` (§3.7);
/// the integer pass extends one sample beyond the top of that range
/// (`int(T_1) + 1`), so the longest delay the postfilter can apply is
/// `143 + 1 = 144`. The carried history buffers reach back this far.
pub const MAX_DELAY: usize = 144;

/// The maximum length of the carried reconstructed-speech / residual
/// history, sized to address `ŝ(n − T)` / `r̂(n − k)` for the longest
/// reachable delay [`MAX_DELAY`] across a whole 40-sample subframe.
const HIST_LEN: usize = MAX_DELAY;

/// The eq (82)/(83) outcome of the delay search for one subframe: the
/// chosen (possibly fractional) delay and the (possibly zero) long-term
/// gain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LongTermDecision {
    /// The best **integer** delay `T_0` (eq (80)) in
    /// `[int(T_1) − 1, int(T_1) + 1]`. Always reported even when the
    /// filter is disabled (so callers can inspect the search outcome).
    /// The full fractional delay the postfilter applies is
    /// `T_0 + frac/8`.
    pub delay: usize,
    /// The §4.2.1 second-pass **fractional** part of the delay, in
    /// eighths of a sample (`0 … 7`, resolution 1/8 around `T_0`). The
    /// applied long-term delay is `delay + frac/8`. `0` for the
    /// integer-only outcome (and whenever the filter is disabled).
    pub frac: usize,
    /// eq (83) long-term gain `g_l`, bounded `0 ≤ g_l ≤ 1`. Exactly
    /// `0.0` when the eq (82) disable test fires (long-term prediction
    /// gain below 3 dB) — in which case eq (78) degenerates to the
    /// identity `H_p(z) = 1`.
    pub gain: f32,
}

/// Stateful §4.2.1 long-term (harmonic) postfilter `H_p(z)`
/// (eqs (78)–(83), integer-delay form).
///
/// Owns the reconstructed-speech history (for the eq (79) residual and
/// the eq (78) `z⁻ᵀ` delay) and the residual history (for the eq (80)
/// correlation), both zero-initialised per clause 4.3 and carried across
/// subframes.
#[derive(Debug, Clone)]
pub struct LongTermPostfilter {
    /// Reconstructed-speech history `[ŝ(n−1) … ŝ(n−HIST_LEN)]`
    /// (`s_hist[0]` = most recent `ŝ(n−1)`). Zero-init per clause 4.3.
    s_hist: [f32; HIST_LEN],
    /// Residual history `[r̂(n−1) … r̂(n−HIST_LEN)]` (`r_hist[0]` = most
    /// recent `r̂(n−1)`). Zero-init per clause 4.3.
    r_hist: [f32; HIST_LEN],
}

impl Default for LongTermPostfilter {
    fn default() -> Self {
        Self::new()
    }
}

impl LongTermPostfilter {
    /// Build the postfilter with the clause-4.3 all-zero start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            s_hist: [0.0; HIST_LEN],
            r_hist: [0.0; HIST_LEN],
        }
    }

    /// Borrow the reconstructed-speech history `[ŝ(n−1) … ŝ(n−HIST_LEN)]`.
    #[must_use]
    pub fn s_hist(&self) -> &[f32; HIST_LEN] {
        &self.s_hist
    }

    /// Borrow the residual history `[r̂(n−1) … r̂(n−HIST_LEN)]`.
    #[must_use]
    pub fn r_hist(&self) -> &[f32; HIST_LEN] {
        &self.r_hist
    }

    /// eq (79) residual for one subframe: `r̂(n) = ŝ(n) + Σ γ_n^i·â_i·ŝ(n−i)`.
    ///
    /// `s` is the 40-sample reconstructed-speech subframe `ŝ(n)`; `a`
    /// holds the per-subframe quantised LP coefficients (slots `0 … 9` =
    /// `â_1 … â_10`). The `ŝ(n−i)` history for the first few samples is
    /// drawn from the carried [`Self::s_hist`]. The weighting `γ_n^i·â_i`
    /// is the short-term postfilter numerator
    /// [`ShortTermPostfilter::weighted_num`] (`γ_n` = [`GAMMA_N`]), shared
    /// so the two stages cannot drift.
    fn residual(&self, s: &[f32; SUBFRAME_SIZE], a: &[f32; M]) -> [f32; SUBFRAME_SIZE] {
        let wn = ShortTermPostfilter::weighted_num(a);
        let mut r = [0.0f32; SUBFRAME_SIZE];
        for (n, &sn) in s.iter().enumerate() {
            let mut acc = sn;
            for (i, &w) in wn.iter().enumerate() {
                // ŝ(n − (i+1)): in-subframe for n > i, else from history.
                let lag = i + 1;
                let sv = if n >= lag {
                    s[n - lag]
                } else {
                    self.s_hist[lag - 1 - n]
                };
                acc += w * sv;
            }
            r[n] = acc;
        }
        r
    }

    /// Fetch the residual at absolute subframe index `idx − k`, drawing
    /// in-subframe samples from `r` and earlier ones from the carried
    /// [`Self::r_hist`]. `idx` is `0 … 39`, `k ≥ 1`.
    fn residual_at(&self, r: &[f32; SUBFRAME_SIZE], idx: usize, k: usize) -> f32 {
        if idx >= k {
            r[idx - k]
        } else {
            // r̂(n−k) with n−k < 0 → history index (k − idx − 1).
            let h = k - idx - 1;
            if h < HIST_LEN {
                self.r_hist[h]
            } else {
                0.0
            }
        }
    }

    /// Fetch the reconstructed speech at absolute subframe index `idx − k`,
    /// drawing in-subframe samples from `s` and earlier ones from the
    /// carried [`Self::s_hist`]. `idx` is `0 … 39`, `k ≥ 1`. Mirrors
    /// [`Self::residual_at`] for the speech history (needed by the eq (78)
    /// `z^{-T}` delay of the postfilter input).
    fn speech_at(&self, s: &[f32; SUBFRAME_SIZE], idx: usize, k: usize) -> f32 {
        if idx >= k {
            s[idx - k]
        } else {
            let h = k - idx - 1;
            if h < HIST_LEN {
                self.s_hist[h]
            } else {
                0.0
            }
        }
    }

    /// Interpolate the reconstructed **speech** at the fractional delay
    /// `T_0 + frac/8` for output sample `idx`, using the same kernel /
    /// tap-offset convention as [`Self::residual_interp`]. This realises
    /// the eq (78) `z^{-T}` fractional delay of the postfilter input. The
    /// long (length-129) filter is used for a non-integer `frac` so the
    /// applied delay matches the fraction the search settled on.
    fn speech_interp(&self, s: &[f32; SUBFRAME_SIZE], idx: usize, t0: usize, frac: usize) -> f32 {
        if frac == 0 {
            return self.speech_at(s, idx, t0);
        }
        let taps = postfilter_interp_long(frac);
        let half = POSTFILTER_INTERP_LONG_TAPS / 2;
        let mut acc = 0.0f32;
        for (j, &h) in taps.iter().enumerate() {
            let k = (t0 + half).wrapping_sub(j);
            acc += (h as f32 / 32768.0) * self.speech_at(s, idx, k);
        }
        acc
    }

    /// Interpolate the **residual** at the fractional delay `T_0 + frac/8`
    /// for output sample `idx` (`0 … 39`), drawing residual samples from
    /// the current subframe `r` and the carried history.
    ///
    /// `frac == 0` is the integer case, returning `r̂(idx − T_0)` directly.
    /// For `frac = 1 … 7` the interpolation kernel ([`postfilter_interp_short`]
    /// when `long == false`, length-33 first pass; [`postfilter_interp_long`]
    /// when `long == true`, length-129 refinement) is convolved with the
    /// residual samples centred on the integer delay: short taps span the
    /// offsets `−2 … +1`, long taps span `−8 … +7` (tap `j` ↔ residual
    /// sample `idx − T_0 − (HALF − j)`, where `HALF` is half the tap
    /// count). A larger `frac` pulls the dominant tap toward the
    /// next-deeper integer delay, realising the 1/8-resolution fractional
    /// delay of clause 4.2.1.
    fn residual_interp(
        &self,
        r: &[f32; SUBFRAME_SIZE],
        idx: usize,
        t0: usize,
        frac: usize,
        long: bool,
    ) -> f32 {
        if frac == 0 {
            return self.residual_at(r, idx, t0);
        }
        let (taps, half): (&[i16], usize) = if long {
            (
                postfilter_interp_long(frac),
                POSTFILTER_INTERP_LONG_TAPS / 2,
            )
        } else {
            (
                postfilter_interp_short(frac),
                POSTFILTER_INTERP_SHORT_TAPS / 2,
            )
        };
        let mut acc = 0.0f32;
        for (j, &h) in taps.iter().enumerate() {
            // tap j applies to r̂(idx − (T_0 + (half − j))) = the residual
            // delayed by k = T_0 + half − j. `half − j` ranges over the
            // tap offsets (e.g. short: +2 … −1).
            let k = (t0 + half).wrapping_sub(j);
            let rk = self.residual_at(r, idx, k);
            acc += (h as f32 / 32768.0) * rk;
        }
        acc
    }

    /// Compute the eq (80)/(81) correlation `R(T) = Σ r̂(n)·r̂_T(n)` and the
    /// delayed energy `E_T = Σ r̂_T(n)²` for a candidate fractional delay
    /// `T_0 + frac/8`, using the short (`long == false`) or long
    /// (`long == true`) interpolation kernel.
    fn fractional_corr_energy(
        &self,
        r: &[f32; SUBFRAME_SIZE],
        t0: usize,
        frac: usize,
        long: bool,
    ) -> (f32, f32) {
        let mut corr = 0.0f32;
        let mut energy = 0.0f32;
        for (idx, &rn) in r.iter().enumerate() {
            let rk = self.residual_interp(r, idx, t0, frac, long);
            corr += rn * rk;
            energy += rk * rk;
        }
        (corr, energy)
    }

    /// Run the clause-4.2.1 **two-pass** delay search (integer + 1/8
    /// fractional refinement) + gain decision for one subframe.
    ///
    /// `r` is the eq (79) residual for the subframe; `int_t1` is the
    /// integer part `int(T_1)` of the transmitted first-subframe pitch
    /// delay.
    ///
    /// **First pass** (eq (80)): the best integer `T_0` over
    /// `int(T_1) − 1 … int(T_1) + 1` (clamped to `≥ 1` and to
    /// [`MAX_DELAY`]) maximising the correlation `R(k)`.
    ///
    /// **Second pass** (clause 4.2.1, eq (81)): the best fractional offset
    /// `frac ∈ {0 … 7}` (resolution 1/8 around `T_0`) maximising the
    /// pseudo-normalised correlation `R(T)² / E_T`, where `E_T` is the
    /// energy of the interpolated delayed residual. The search uses the
    /// length-33 short filter first; the chosen non-integer `frac` is then
    /// re-evaluated with the length-129 long filter, and the long-filter
    /// fraction replaces the short-filter one **only if it increases**
    /// `R(T)² / E_T` (the spec's "new signal replaces the previous one
    /// only if the longer filter increases the value of `R′(T)`").
    ///
    /// The eq (82) disable test then sets `g_l = 0` if the squared
    /// normalised correlation is below [`ENABLE_THRESHOLD`], else the
    /// eq (83) gain (bounded `[0, 1]`) is returned. Returns the
    /// [`LongTermDecision`] (integer delay + fractional part + gain)
    /// without mutating state.
    #[must_use]
    pub fn decide(&self, r: &[f32; SUBFRAME_SIZE], int_t1: usize) -> LongTermDecision {
        // Energy of the current residual, Σ r̂(n)² (eq (82) denominator).
        let energy: f32 = r.iter().map(|v| v * v).sum();

        // --- First pass: best integer T_0 (eq (80)). ---
        let lo = int_t1.saturating_sub(1).max(1);
        let hi = (int_t1 + 1).min(MAX_DELAY);

        let mut t0 = lo.max(1);
        let mut best_int_corr = f32::NEG_INFINITY;
        for k in lo..=hi {
            let mut corr = 0.0f32;
            for (idx, &rn) in r.iter().enumerate() {
                corr += rn * self.residual_at(r, idx, k);
            }
            if corr > best_int_corr {
                best_int_corr = corr;
                t0 = k;
            }
        }

        // --- Second pass: best fractional frac/8 around T_0 (eq (81)). ---
        // Maximise the pseudo-normalised correlation R(T)²/E_T over the
        // eight 1/8 offsets, short (length-33) filter first.
        let mut best_frac = 0usize;
        let mut best_num = 0.0f32; // Σ r̂·r̂_T  (eq (83) numerator)
        let mut best_den = 0.0f32; // Σ r̂_T²    (eq (83) denominator)
        let mut best_score = f32::NEG_INFINITY; // R(T)²/E_T (only when R>0)
        for frac in 0..POSTFILTER_FRAC_RES {
            let (corr, e_t) = self.fractional_corr_energy(r, t0, frac, false);
            // A negative correlation cannot beat a disabled filter; only
            // positive-correlation candidates compete for the maximum.
            let score = if corr > 0.0 && e_t > 0.0 {
                (corr * corr) / e_t
            } else {
                f32::NEG_INFINITY
            };
            if score > best_score {
                best_score = score;
                best_frac = frac;
                best_num = corr;
                best_den = e_t;
            }
        }

        // Long-filter refinement: re-evaluate the chosen non-integer frac
        // with the length-129 filter; keep it only if R(T)²/E_T rises.
        if best_frac != 0 {
            let (corr_l, e_l) = self.fractional_corr_energy(r, t0, best_frac, true);
            if corr_l > 0.0 && e_l > 0.0 {
                let score_l = (corr_l * corr_l) / e_l;
                if score_l > best_score {
                    best_score = score_l;
                    best_num = corr_l;
                    best_den = e_l;
                }
            }
        }
        let _ = best_score;

        // eq (82) long-term-prediction-gain disable: g_l = 0 when
        // R′(T)² / Σ r̂(n)² < 0.5, where R′(T)² = (Σ r̂·r̂_T)² / (Σ r̂_T²).
        let enabled = best_num > 0.0
            && best_den > 0.0
            && energy > 0.0
            && (best_num * best_num) >= ENABLE_THRESHOLD * energy * best_den;

        let (frac, gain) = if enabled {
            // eq (83) g_l = Σ r̂·r̂_T / Σ r̂_T², bounded [0, 1].
            (best_frac, (best_num / best_den).clamp(0.0, 1.0))
        } else {
            // Disabled → eq (78) is the identity; report the integer
            // anchor with no fractional offset.
            (0, 0.0)
        };

        LongTermDecision {
            delay: t0,
            frac,
            gain,
        }
    }

    /// Filter one 40-sample reconstructed-speech subframe `ŝ(n)` through
    /// the §4.2.1 long-term postfilter (eqs (78)–(83), integer-delay
    /// form), advancing the carried history.
    ///
    /// `s` is the §4.1.6 reconstructed speech for the subframe; `a` holds
    /// the per-subframe quantised LP coefficients (`â_1 … â_10`);
    /// `int_t1` is the integer part `int(T_1)` of the transmitted
    /// first-subframe pitch delay (clause 4.2.1 anchors the integer
    /// search on `int(T_1)`).
    ///
    /// Returns `(out, decision)` where `out` is the long-term-postfiltered
    /// subframe `H_p(z)·ŝ` (eq (78)) and `decision` is the
    /// [`LongTermDecision`] (chosen delay + gain) for inspection / tests.
    #[must_use]
    pub fn filter_subframe(
        &mut self,
        s: &[f32; SUBFRAME_SIZE],
        a: &[f32; M],
        int_t1: usize,
    ) -> ([f32; SUBFRAME_SIZE], LongTermDecision) {
        let r = self.residual(s, a);
        let decision = self.decide(&r, int_t1);

        // eq (78): out(n) = (ŝ(n) + γ_p·g_l·ŝ(n−T)) / (1 + γ_p·g_l), where
        // T = T_0 + frac/8 is the chosen fractional delay; ŝ(n−T) is the
        // long-filter-interpolated reconstructed speech at that delay.
        let gl = decision.gain;
        let scale = GAMMA_P * gl;
        let inv = 1.0 / (1.0 + scale);
        let t0 = decision.delay;
        let frac = decision.frac;

        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &sn) in s.iter().enumerate() {
            // ŝ(n − T): the (possibly fractional) delayed speech.
            let s_delayed = if gl == 0.0 {
                0.0 // unused when the filter is disabled (scale = 0)
            } else {
                self.speech_interp(s, n, t0, frac)
            };
            out[n] = inv * (sn + scale * s_delayed);
        }

        // Advance both histories by the 40 new samples (most-recent at
        // index 0). Shift the old window down by SUBFRAME_SIZE, then write
        // the new subframe in reverse so index 0 holds the last sample.
        self.push_history(s, &r);

        (out, decision)
    }

    /// Push a 40-sample reconstructed-speech subframe and its residual
    /// into the carried histories (most-recent at index 0).
    fn push_history(&mut self, s: &[f32; SUBFRAME_SIZE], r: &[f32; SUBFRAME_SIZE]) {
        // Shift existing history down by one subframe.
        self.s_hist
            .copy_within(0..HIST_LEN - SUBFRAME_SIZE, SUBFRAME_SIZE);
        self.r_hist
            .copy_within(0..HIST_LEN - SUBFRAME_SIZE, SUBFRAME_SIZE);
        // Newest sample ŝ(last) at index 0, ŝ(last−1) at index 1, …
        for j in 0..SUBFRAME_SIZE {
            self.s_hist[j] = s[SUBFRAME_SIZE - 1 - j];
            self.r_hist[j] = r[SUBFRAME_SIZE - 1 - j];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::short_term_postfilter::GAMMA_N;

    /// New filter starts zeroed (clause 4.3 — neither history is in
    /// Table 9).
    #[test]
    fn new_filter_has_zero_state() {
        let pf = LongTermPostfilter::new();
        assert_eq!(pf.s_hist(), &[0.0; HIST_LEN]);
        assert_eq!(pf.r_hist(), &[0.0; HIST_LEN]);
    }

    /// The weight factors are the spec values (clause 4.2.1 / 4.2.2).
    #[test]
    fn weight_factors_match_spec() {
        assert!((GAMMA_P - 0.5).abs() < 1e-9);
        assert!((GAMMA_N - 0.55).abs() < 1e-9);
        assert!((ENABLE_THRESHOLD - 0.5).abs() < 1e-9);
    }

    /// eq (79) residual with all-zero `â_i` is the identity: `r̂(n) = ŝ(n)`
    /// (the numerator `Â(z/γ_n)` degenerates to `1`).
    #[test]
    fn residual_flat_filter_is_identity() {
        let pf = LongTermPostfilter::new();
        let a = [0.0f32; M];
        let s: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let r = pf.residual(&s, &a);
        for n in 0..SUBFRAME_SIZE {
            assert!((r[n] - s[n]).abs() < 1e-6, "n={n}");
        }
    }

    /// eq (79) residual matches an independent direct convolution of the
    /// weighted numerator over `ŝ`, with zeroed history for the leading
    /// samples.
    #[test]
    fn residual_matches_direct_convolution() {
        let pf = LongTermPostfilter::new();
        let mut a = [0.0f32; M];
        a[0] = 0.7;
        a[1] = -0.3;
        a[4] = 0.1;
        let s: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| ((n * 5) % 17) as f32 - 8.0);
        let r = pf.residual(&s, &a);

        let wn = ShortTermPostfilter::weighted_num(&a);
        for n in 0..SUBFRAME_SIZE {
            let mut expect = s[n];
            for (i, &w) in wn.iter().enumerate() {
                let lag = i + 1;
                let sv = if n >= lag { s[n - lag] } else { 0.0 };
                expect += w * sv;
            }
            assert!((r[n] - expect).abs() < 1e-4, "n={n}");
        }
    }

    /// A purely periodic residual at exactly the search delay maximises
    /// eq (80) and yields a long-term gain near 1 (eq (83) numerator ≈
    /// denominator) — the filter is enabled.
    #[test]
    fn periodic_residual_enables_filter_with_unit_gain() {
        let mut pf = LongTermPostfilter::new();
        let period = 40usize;
        // Pre-load the residual history with the same period so r̂(n−T)
        // is well-defined for the whole subframe.
        for j in 0..HIST_LEN {
            // history index j = lag (j+1); value at ŝ(n−(j+1)).
            let phase = (HIST_LEN - 1 - j) % period;
            pf.r_hist[j] = (phase as f32 * core::f32::consts::TAU / period as f32).sin();
        }
        // Residual for this subframe continues the same sine.
        let r: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| {
            ((HIST_LEN + n) as f32 * core::f32::consts::TAU / period as f32).sin()
        });
        let d = pf.decide(&r, period);
        assert_eq!(d.delay, period);
        assert!(d.gain > 0.9, "gain {} should be near 1", d.gain);
        assert!(d.gain <= 1.0);
    }

    /// A residual uncorrelated with its delayed copy fails the eq (82)
    /// disable test → `g_l = 0` and eq (78) is the identity.
    #[test]
    fn uncorrelated_residual_disables_filter() {
        let mut pf = LongTermPostfilter::new();
        // White-ish residual, history left zero.
        let r: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| if n % 2 == 0 { 1.0 } else { -1.0 });
        let d = pf.decide(&r, 40);
        assert_eq!(d.gain, 0.0, "uncorrelated residual must disable filter");

        // eq (78) with g_l = 0 passes ŝ through unchanged.
        let a = [0.0f32; M];
        let s: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let (out, dec) = pf.filter_subframe(&s, &a, 40);
        assert_eq!(dec.gain, 0.0);
        for n in 0..SUBFRAME_SIZE {
            assert!((out[n] - s[n]).abs() < 1e-6, "n={n}");
        }
    }

    /// The eq (82) threshold gates enabling. A residual that is an exact
    /// `T`-periodic copy of its history (`r̂(n) = r̂(n−T)`) has unit
    /// normalised correlation and enables; an alternating-sign residual
    /// with zero history is anti-correlated at the search lags and
    /// disables.
    #[test]
    fn enable_threshold_is_one_half() {
        let t = 40usize;
        // Exact-period residual: history one full period of a clean wave,
        // current subframe its exact continuation → R′(T)²/Σr̂² ≈ 1 > 0.5.
        let mut pf = LongTermPostfilter::new();
        let wave = |n: usize| (n as f32 * core::f32::consts::TAU / t as f32).sin() + 0.3;
        for j in 0..HIST_LEN {
            pf.r_hist[j] = wave((HIST_LEN - 1 - j) % t);
        }
        let r: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| wave((HIST_LEN + n) % t));
        assert!(
            pf.decide(&r, t).gain > 0.0,
            "exact-period residual must enable (normalised correlation ≈ 1 ≥ 0.5)"
        );

        // No long-term structure (history zero, white residual) → the
        // delayed energy is ~0 / uncorrelated → below 0.5 → disabled.
        let pf2 = LongTermPostfilter::new();
        let white: [f32; SUBFRAME_SIZE] =
            core::array::from_fn(|n| if n % 2 == 0 { 1.0 } else { -1.0 });
        assert_eq!(
            pf2.decide(&white, t).gain,
            0.0,
            "no long-term structure must disable"
        );
    }

    /// eq (78) with a known gain applies the harmonic emphasis:
    /// out(n) = (ŝ(n) + γ_p·g_l·ŝ(n−T)) / (1 + γ_p·g_l).
    #[test]
    fn filter_applies_eq78_with_history() {
        let mut pf = LongTermPostfilter::new();
        // Force a periodic ŝ so the integer search locks a known delay
        // and a positive gain; verify the eq (78) algebra at one sample.
        let period = 40usize;
        for j in 0..HIST_LEN {
            let phase = (HIST_LEN - 1 - j) % period;
            let v = (phase as f32 * core::f32::consts::TAU / period as f32).sin();
            pf.s_hist[j] = v;
            pf.r_hist[j] = v;
        }
        let a = [0.0f32; M]; // identity residual → r̂ = ŝ
        let s: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| {
            ((HIST_LEN + n) as f32 * core::f32::consts::TAU / period as f32).sin()
        });
        let s_snapshot = s;
        let s_hist_before = pf.s_hist;
        let (out, dec) = pf.filter_subframe(&s, &a, period);
        // An exact-integer-period signal scores highest at the integer
        // delay (the fractional kernels slightly attenuate the energy), so
        // the second pass settles on frac = 0 and the eq (78) algebra is
        // the plain integer-delay form.
        assert_eq!(dec.frac, 0, "exact period must keep the integer delay");
        let scale = GAMMA_P * dec.gain;
        let inv = 1.0 / (1.0 + scale);
        for n in 0..SUBFRAME_SIZE {
            let s_delayed = if n >= dec.delay {
                s_snapshot[n - dec.delay]
            } else {
                s_hist_before[dec.delay - n - 1]
            };
            let expect = inv * (s_snapshot[n] + scale * s_delayed);
            assert!((out[n] - expect).abs() < 1e-5, "n={n}");
        }
    }

    /// State carries across subframes: history after a call holds the
    /// just-filtered subframe with the last sample at index 0.
    #[test]
    fn history_advances_most_recent_first() {
        let mut pf = LongTermPostfilter::new();
        let a = [0.0f32; M];
        let s: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n + 1) as f32);
        let _ = pf.filter_subframe(&s, &a, 40);
        // ŝ(last) = 40 at index 0; ŝ(last−1) = 39 at index 1, …
        assert!((pf.s_hist[0] - 40.0).abs() < 1e-6);
        assert!((pf.s_hist[1] - 39.0).abs() < 1e-6);
        // residual = ŝ for the flat filter.
        assert!((pf.r_hist[0] - 40.0).abs() < 1e-6);
        // A second subframe pushes the first one down by 40.
        let s2: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n + 101) as f32);
        let _ = pf.filter_subframe(&s2, &a, 40);
        assert!((pf.s_hist[0] - 140.0).abs() < 1e-6);
        // ŝ(last) of the *first* subframe (40) is now SUBFRAME_SIZE down.
        assert!((pf.s_hist[SUBFRAME_SIZE] - 40.0).abs() < 1e-6);
    }

    /// Gain is always bounded to `[0, 1]` (eq (83) "bounded by 1") and
    /// delay stays inside the clamped `[int(T1)−1, int(T1)+1]` window for
    /// arbitrary inputs.
    #[test]
    fn gain_bounded_and_delay_in_window() {
        let mut pf = LongTermPostfilter::new();
        let a: [f32; M] = [
            0.9, -0.4, 0.2, -0.1, 0.05, -0.03, 0.02, -0.01, 0.005, -0.002,
        ];
        for k in 0..20 {
            let s: [f32; SUBFRAME_SIZE] =
                core::array::from_fn(|n| (((n + k) as f32 * 0.3).sin()) * 1000.0);
            let int_t1 = 30 + k;
            let (out, d) = pf.filter_subframe(&s, &a, int_t1);
            assert!((0.0..=1.0).contains(&d.gain), "gain {}", d.gain);
            let lo = int_t1.saturating_sub(1).max(1);
            let hi = (int_t1 + 1).min(MAX_DELAY);
            assert!(
                (lo..=hi).contains(&d.delay),
                "delay {} window [{lo},{hi}]",
                d.delay
            );
            assert!(out.iter().all(|v| v.is_finite()));
        }
    }

    /// Two filters fed identical subframes stay in lockstep — all state
    /// is owned, no hidden globals.
    #[test]
    fn deterministic() {
        let a: [f32; M] = [0.5, -0.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) * 2.0 - 39.0);
        let mut p = LongTermPostfilter::new();
        let mut q = LongTermPostfilter::new();
        for _ in 0..6 {
            assert_eq!(p.filter_subframe(&s, &a, 50), q.filter_subframe(&s, &a, 50));
        }
    }

    /// A disabled filter (g_l = 0) is exactly the identity over the whole
    /// subframe regardless of the delay — eq (78) collapses to `H_p = 1`.
    #[test]
    fn disabled_filter_is_identity() {
        let mut pf = LongTermPostfilter::new();
        let a = [0.0f32; M];
        // Alternating sign residual → never correlated → disabled.
        let s: [f32; SUBFRAME_SIZE] =
            core::array::from_fn(|n| if n % 2 == 0 { 500.0 } else { -500.0 });
        let (out, d) = pf.filter_subframe(&s, &a, 40);
        assert_eq!(d.gain, 0.0);
        assert_eq!(d.frac, 0, "disabled filter reports no fractional offset");
        for n in 0..SUBFRAME_SIZE {
            assert!((out[n] - s[n]).abs() < 1e-6, "n={n}");
        }
    }

    /// The reported fractional offset is always a valid 1/8 phase
    /// (`0 … 7`) for arbitrary inputs.
    #[test]
    fn fractional_offset_in_range() {
        let mut pf = LongTermPostfilter::new();
        let a: [f32; M] = [0.6, -0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        for k in 0..30 {
            let s: [f32; SUBFRAME_SIZE] =
                core::array::from_fn(|n| (((n + k) as f32 * 0.21).sin()) * 800.0);
            let (_out, d) = pf.filter_subframe(&s, &a, 35 + k);
            assert!(d.frac < POSTFILTER_FRAC_RES, "frac {} out of range", d.frac);
        }
    }

    /// A signal whose period is a non-integer number of samples drives
    /// the second pass to a **non-zero** fractional offset: when the best
    /// integer anchor is the *floor* of the period, the 1/8 refinement
    /// adds a positive fraction `T_0 + frac/8` to reach the true period.
    /// (Exercises the fractional path end-to-end.)
    #[test]
    fn non_integer_period_selects_fractional_delay() {
        let mut pf = LongTermPostfilter::new();
        // Period 40.25 samples: the integer correlation favours 40 (the
        // floor) and the fractional pass adds ≈ 2/8 to reach 40.25.
        let period = 40.25f32;
        for j in 0..HIST_LEN {
            let phase = (HIST_LEN - 1 - j) as f32;
            let v = (phase * core::f32::consts::TAU / period).sin();
            pf.s_hist[j] = v;
            pf.r_hist[j] = v;
        }
        let r: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| {
            ((HIST_LEN + n) as f32 * core::f32::consts::TAU / period).sin()
        });
        // Anchor the integer window on 40 so T_0 is the period floor.
        let d = pf.decide(&r, 40);
        assert!(d.gain > 0.0, "fractional match enables the filter");
        // The applied delay T_0 + frac/8 should approach the true 40.25.
        let applied = d.delay as f32 + d.frac as f32 / 8.0;
        assert!(
            (applied - period).abs() <= 0.5,
            "applied delay {applied} (T_0={}, frac={}) should be near {period}",
            d.delay,
            d.frac
        );
        assert_ne!(
            d.frac, 0,
            "a non-integer period must pick a non-zero 1/8 fraction (got frac=0)"
        );
    }

    /// Fractional interpolation at `frac = 0` is exactly the integer
    /// residual fetch (no kernel applied), so the two paths agree.
    #[test]
    fn frac_zero_interp_is_integer_fetch() {
        let mut pf = LongTermPostfilter::new();
        for j in 0..HIST_LEN {
            pf.r_hist[j] = ((j as f32) * 0.37).cos() * 100.0;
        }
        let r: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| ((n as f32) * 0.5).sin() * 100.0);
        for idx in 0..SUBFRAME_SIZE {
            let a = pf.residual_interp(&r, idx, 30, 0, false);
            let b = pf.residual_at(&r, idx, 30);
            assert!((a - b).abs() < 1e-6, "idx={idx}");
        }
    }

    /// The half-sample fractional phase (`frac = 4`) interpolates a slowly
    /// varying residual to roughly the midpoint of its two integer
    /// neighbours (the kernel is symmetric and near-unity-gain there).
    #[test]
    fn half_sample_interp_is_near_midpoint() {
        let mut pf = LongTermPostfilter::new();
        // A smooth ramp in the history so the midpoint is well-defined.
        for j in 0..HIST_LEN {
            // r_hist[j] holds r̂(n−(j+1)); make it a smooth linear ramp.
            pf.r_hist[j] = 1000.0 - (j as f32) * 3.0;
        }
        let r = [0.0f32; SUBFRAME_SIZE];
        // At idx = 0, integer delay T0 = 30 fetches r̂(−30) = r_hist[29];
        // its neighbour r̂(−31) = r_hist[30]. The half-sample (frac=4)
        // interpolation should land near the average of the two.
        let t0 = 30usize;
        let neighbour_lo = pf.residual_at(&r, 0, t0); // r̂(−30)
        let neighbour_hi = pf.residual_at(&r, 0, t0 + 1); // r̂(−31)
        let mid = 0.5 * (neighbour_lo + neighbour_hi);
        let interp = pf.residual_interp(&r, 0, t0, 4, false);
        assert!(
            (interp - mid).abs() < 30.0,
            "frac=4 interp {interp} should be near midpoint {mid}"
        );
    }

    /// Known-answer sinc-interpolation check: the §4.2.1 interpolation
    /// kernels reproduce a band-limited sinusoid sampled at the exact
    /// fractional delay `T_0 + frac/8`, and the **long** (length-129)
    /// filter is more accurate than the **short** (length-33) one. This
    /// validates the polyphase tap-offset layout end-to-end from the
    /// staged data — independent of any reference — and confirms the
    /// spec's rationale for the long-filter refinement.
    #[test]
    fn fractional_kernels_reproduce_fractional_delay() {
        // Sampled sinusoid as the residual history; r̂(n−k) = sin(2π f (−k))
        // relative to the output index. Load the history with the sine so
        // residual_at / residual_interp draw the right samples.
        let f = 0.07f32;
        let mut pf = LongTermPostfilter::new();
        for j in 0..HIST_LEN {
            // r_hist[j] holds r̂(n−(j+1)); for output index 0 that is the
            // sample at continuous offset −(j+1).
            let n = -((j + 1) as f32);
            pf.r_hist[j] = (core::f32::consts::TAU * f * n).sin();
        }
        let r = [0.0f32; SUBFRAME_SIZE];
        let t0 = 20usize;
        for frac in 1..POSTFILTER_FRAC_RES {
            let truth = (core::f32::consts::TAU * f * (-(t0 as f32) - frac as f32 / 8.0)).sin();
            let short = pf.residual_interp(&r, 0, t0, frac, false);
            let long = pf.residual_interp(&r, 0, t0, frac, true);
            let err_s = (short - truth).abs();
            let err_l = (long - truth).abs();
            assert!(err_s < 0.01, "frac={frac} short err {err_s}");
            assert!(err_l < 0.01, "frac={frac} long err {err_l}");
            // The long filter is at least as accurate as the short one
            // (within a tight margin) — the spec keeps the long result
            // only when it improves the correlation.
            assert!(
                err_l <= err_s + 1e-3,
                "frac={frac}: long err {err_l} should not exceed short err {err_s}"
            );
        }
    }
}
