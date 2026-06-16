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
//! ## Scope of this module — the **integer-delay** postfilter (eqs
//! (78)–(83) with `T = T_0`)
//!
//! This module realises the integer-delay form of clause 4.2.1: the
//! full residual (eq (79)), the integer delay search (eq (80) over the
//! three candidates `int(T_1) − 1 … int(T_1) + 1`), the long-term
//! prediction-gain disable test (eq (82)), the bounded gain (eq (83)),
//! and the eq (78) harmonic filter. For an **integer** delay `r̂_k(n)`
//! reduces to `r̂(n − T_0)`, so eqs (81)–(83) collapse to sums over the
//! integer-delayed residual, all of which are spelled out by the spec
//! prose with no table dependency.
//!
//! The 1/8-resolution **fractional** refinement (the second pass) and
//! its `tab_hup_s` (length 33, 28 stored coefficients) / `tab_hup_l`
//! (length 129, 112 stored coefficients) interpolation filters are NOT
//! wired here: the spec prose names the filter lengths but does not give
//! the per-phase tap layout or the convolution indexing of those tables
//! (it defers them to the electronic-attachment reference), so the
//! fractional pass is a documented docs-gap and a separate follow-up.
//! The integer pass is the spec-complete, prose-sourced core of the
//! stage; an integer-only long-term postfilter is itself a valid
//! harmonic postfilter (the fractional pass only refines the delay /
//! correlation by interpolation, it cannot newly enable a filter the
//! integer test disabled, since the integer delay is one of the search
//! candidates).
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
use crate::tables::M;

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
/// chosen integer delay and the (possibly zero) long-term gain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LongTermDecision {
    /// The best integer delay `T_0` (eq (80)) in
    /// `[int(T_1) − 1, int(T_1) + 1]`. Always reported even when the
    /// filter is disabled (so callers can inspect the search outcome).
    pub delay: usize,
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

    /// Run the clause-4.2.1 integer delay search + gain decision for one
    /// subframe.
    ///
    /// `r` is the eq (79) residual for the subframe; `int_t1` is the
    /// integer part `int(T_1)` of the transmitted first-subframe pitch
    /// delay. The search maximises eq (80) `R(k)` over the three integer
    /// candidates `int(T_1) − 1 … int(T_1) + 1` (clamped to `≥ 1` and to
    /// [`MAX_DELAY`]), retaining the candidate with the largest positive
    /// `R(k)`. The eq (82) disable test then sets `g_l = 0` if the
    /// squared normalised correlation is below [`ENABLE_THRESHOLD`], else
    /// the eq (83) gain (bounded `[0, 1]`) is returned.
    ///
    /// Returns the [`LongTermDecision`] (delay + gain) without mutating
    /// state — the caller applies eq (78) and then advances the history.
    #[must_use]
    pub fn decide(&self, r: &[f32; SUBFRAME_SIZE], int_t1: usize) -> LongTermDecision {
        // Energy of the current residual, Σ r̂(n)² (eq (82) denominator).
        let energy: f32 = r.iter().map(|v| v * v).sum();

        // Integer candidate window [int(T1)−1, int(T1)+1], clamped so the
        // delay stays a usable positive lag within the history reach.
        let lo = int_t1.saturating_sub(1).max(1);
        let hi = (int_t1 + 1).min(MAX_DELAY);

        let mut best_delay = lo.max(1);
        let mut best_corr = f32::NEG_INFINITY; // R(k) (eq (80))
        let mut best_num = 0.0f32; //  Σ r̂(n)·r̂(n−T)  (eq (83) numerator)
        let mut best_den = 0.0f32; //  Σ r̂(n−T)²       (eq (83) denominator)

        for k in lo..=hi {
            // R(k) = Σ_{n=0}^{39} r̂(n)·r̂(n−k); also accumulate the
            // delayed-energy Σ r̂(n−k)² for the eq (83) gain.
            let mut corr = 0.0f32;
            let mut den = 0.0f32;
            for (idx, &rn) in r.iter().enumerate() {
                let rk = self.residual_at(r, idx, k);
                corr += rn * rk;
                den += rk * rk;
            }
            if corr > best_corr {
                best_corr = corr;
                best_delay = k;
                best_num = corr;
                best_den = den;
            }
        }

        // eq (82) long-term-prediction-gain disable: g_l = 0 when
        // R′(T)² / Σ r̂(n)² < 0.5. For an integer delay,
        // R′(T)² = (Σ r̂·r̂_k)² / (Σ r̂_k²) = best_num² / best_den.
        // The test compares (best_num² / best_den) / energy < 0.5.
        let enabled = best_num > 0.0
            && best_den > 0.0
            && energy > 0.0
            && (best_num * best_num) >= ENABLE_THRESHOLD * energy * best_den;

        let gain = if enabled {
            // eq (83) g_l = Σ r̂·r̂_k / Σ r̂_k², bounded [0, 1].
            (best_num / best_den).clamp(0.0, 1.0)
        } else {
            0.0
        };

        LongTermDecision {
            delay: best_delay,
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

        // eq (78): out(n) = (ŝ(n) + γ_p·g_l·ŝ(n−T)) / (1 + γ_p·g_l).
        let gl = decision.gain;
        let scale = GAMMA_P * gl;
        let inv = 1.0 / (1.0 + scale);
        let t = decision.delay;

        let mut out = [0.0f32; SUBFRAME_SIZE];
        for (n, &sn) in s.iter().enumerate() {
            // ŝ(n − T): in-subframe for n ≥ T, else from history.
            let s_delayed = if gl == 0.0 {
                0.0 // unused when the filter is disabled (scale = 0)
            } else if n >= t {
                s[n - t]
            } else {
                let h = t - n - 1;
                if h < HIST_LEN {
                    self.s_hist[h]
                } else {
                    0.0
                }
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
        for n in 0..SUBFRAME_SIZE {
            assert!((out[n] - s[n]).abs() < 1e-6, "n={n}");
        }
    }
}
