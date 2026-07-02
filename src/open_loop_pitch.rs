//! §3.4 **open-loop pitch analysis** — estimates the once-per-frame
//! candidate pitch delay `T_op` that bounds the §3.7 closed-loop
//! adaptive-codebook search.
//!
//! Seventh encoder-side stage. Runs on the weighted speech `sw(n)` from
//! [`crate::perceptual_weighting`] (eq (33)).
//!
//! ## Spec source — clause 3.4, equations (34), (35) + selection logic
//!
//! (Transcribed from the EPUB prose + rasters `images/eq34.jpg`,
//! `eq35.jpg`, `f0017-01.jpg`.)
//!
//! * **eq (34)** — the correlation over the 80-sample frame:
//!   `R(k) = Σ_{n=0}^{79} sw(n)·sw(n−k)` (negative indices reach into
//!   the weighted-speech history).
//! * Three maxima of `R(k)` are found in the ranges `k = 80 … 143`
//!   (`t₁`), `k = 40 … 79` (`t₂`), and `k = 20 … 39` (`t₃`).
//! * **eq (35)** — each retained maximum is normalised:
//!   `R'(t_i) = R(t_i) / √(Σ_n sw²(n − t_i))` (the denominator sums the
//!   lagged energy over the same `n = 0 … 79` window).
//! * **Selection** (spec procedure `f0017-01`, favouring lower delays to
//!   avoid pitch multiples):
//!   ```text
//!   T_op = t₁; R'(T_op) = R'(t₁)
//!   if R'(t₂) ≥ 0.85·R'(T_op):  R'(T_op) = R'(t₂); T_op = t₂
//!   if R'(t₃) ≥ 0.85·R'(T_op):  R'(T_op) = R'(t₃); T_op = t₃
//!   ```
//!
//! Note the ranges are searched from the **longest** delays (`t₁` in
//! `80…143`) down; the ≥ 0.85 tests progressively prefer the shorter
//! candidates `t₂` then `t₃` when their normalised correlation is
//! competitive.
//!
//! ## History requirement
//!
//! `R(k)` with `k` up to 143 needs 143 samples of weighted speech
//! *before* the current frame. The caller supplies a
//! [`PIT_BUFFER`]-sample buffer laid out as
//! `[history (143) | current frame (80)]`; the analysis window
//! `n = 0 … 79` addresses the current-frame region.

/// Samples per frame (clause 2).
pub const L_FRAME: usize = 80;
/// Maximum open-loop delay searched (clause 3.4 range `80 … 143`).
pub const PIT_MAX: usize = 143;
/// Minimum open-loop delay searched (clause 3.4 range `20 … 39`).
pub const PIT_MIN: usize = 20;
/// Total buffer the analysis consumes: `PIT_MAX` history samples then
/// the current 80-sample frame.
pub const PIT_BUFFER: usize = PIT_MAX + L_FRAME;

/// The three (range, result) delay sections of clause 3.4, longest
/// first: `t₁ ∈ [80, 143]`, `t₂ ∈ [40, 79]`, `t₃ ∈ [20, 39]`.
const SECTIONS: [(usize, usize); 3] = [(80, 143), (40, 79), (20, 39)];

/// The favour-lower-delays comparison factor (spec selection procedure).
const FAVOUR_FACTOR: f32 = 0.85;

/// Result of the §3.4 analysis: the selected open-loop delay and its
/// normalised correlation (useful for voicing decisions downstream).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpenLoopPitch {
    /// The best open-loop delay `T_op` (`20 … 143`).
    pub t_op: usize,
    /// The eq (35) normalised correlation `R'(T_op)`.
    pub r_norm: f32,
}

/// eq (34): `R(k) = Σ_{n=0}^{79} sw(n)·sw(n−k)` where `sw` is addressed
/// through `buf` with the current frame starting at index [`PIT_MAX`].
fn correlation(buf: &[f32; PIT_BUFFER], k: usize) -> f32 {
    let mut acc = 0.0f32;
    for n in 0..L_FRAME {
        acc += buf[PIT_MAX + n] * buf[PIT_MAX + n - k];
    }
    acc
}

/// eq (35) denominator: the lagged energy `Σ_{n=0}^{79} sw²(n − k)`.
fn lagged_energy(buf: &[f32; PIT_BUFFER], k: usize) -> f32 {
    let mut acc = 0.0f32;
    for n in 0..L_FRAME {
        let v = buf[PIT_MAX + n - k];
        acc += v * v;
    }
    acc
}

/// §3.4 open-loop pitch estimation over one frame.
///
/// `buf` is `[143 samples of weighted-speech history | 80 samples of the
/// current frame's weighted speech]`. Returns the selected `T_op` and
/// its normalised correlation.
#[must_use]
pub fn open_loop_pitch(buf: &[f32; PIT_BUFFER]) -> OpenLoopPitch {
    // Step 1: the raw-correlation maximum in each of the three sections.
    let mut cand = [(0usize, 0.0f32); 3];
    for (s, &(lo, hi)) in SECTIONS.iter().enumerate() {
        let mut best_k = lo;
        let mut best_r = f32::NEG_INFINITY;
        for k in lo..=hi {
            let r = correlation(buf, k);
            if r > best_r {
                best_r = r;
                best_k = k;
            }
        }
        cand[s] = (best_k, best_r);
    }

    // Step 2: eq (35) normalisation of the three retained maxima.
    let mut norm = [0.0f32; 3];
    for (s, &(k, r)) in cand.iter().enumerate() {
        let e = lagged_energy(buf, k);
        norm[s] = if e > 0.0 { r / e.sqrt() } else { 0.0 };
    }

    // Step 3: favour-lower-delays selection (t₁ → t₂ → t₃).
    let mut t_op = cand[0].0;
    let mut r_best = norm[0];
    if norm[1] >= FAVOUR_FACTOR * r_best {
        r_best = norm[1];
        t_op = cand[1].0;
    }
    if norm[2] >= FAVOUR_FACTOR * r_best {
        r_best = norm[2];
        t_op = cand[2].0;
    }

    OpenLoopPitch {
        t_op,
        r_norm: r_best,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a periodic buffer with the given integer period.
    fn periodic_buf(period: usize) -> [f32; PIT_BUFFER] {
        std::array::from_fn(|n| {
            let phase = (n % period) as f32 / period as f32;
            (std::f32::consts::TAU * phase).sin()
                + 0.3 * (2.0 * std::f32::consts::TAU * phase).sin()
        })
    }

    /// A signal with period 30 (inside the t₃ range) is picked at 30 —
    /// the favour-lower-delays rule must reject the 60/90/120 multiples
    /// that score equally raw.
    #[test]
    fn period_30_beats_its_multiples() {
        let buf = periodic_buf(30);
        let res = open_loop_pitch(&buf);
        assert_eq!(res.t_op, 30, "expected fundamental, got {}", res.t_op);
        assert!(res.r_norm > 0.8, "strongly periodic: R'={}", res.r_norm);
    }

    /// Period 55 lives in the t₂ range; its double 110 is in t₁. The
    /// selection must land on 55.
    #[test]
    fn period_55_beats_double() {
        let buf = periodic_buf(55);
        let res = open_loop_pitch(&buf);
        assert_eq!(res.t_op, 55, "expected fundamental, got {}", res.t_op);
    }

    /// A long-period signal (period 100, only representable in t₁) is
    /// found near there. Because the eq (34) raw-correlation maximum is
    /// taken **before** the eq (35) normalisation and the 80-sample
    /// analysis window covers less than one full period, the raw peak
    /// can land a couple of samples off the exact period (the lagged
    /// window's energy varies with position) — that is the spec's own
    /// two-step behaviour, so the assertion allows a small neighbourhood
    /// of the fundamental (or of its half-alias 50, which the
    /// favour-lower-delays rule may legitimately prefer).
    #[test]
    fn long_period_stays_near_fundamental() {
        let buf = periodic_buf(100);
        let res = open_loop_pitch(&buf);
        let near = |c: usize, w: usize| res.t_op >= c - w && res.t_op <= c + w;
        assert!(
            near(100, 3) || near(50, 3),
            "expected ~100 or its half-alias ~50, got {}",
            res.t_op
        );
    }

    /// Silence yields a well-defined result (no NaN), r_norm = 0.
    #[test]
    fn silence_is_safe() {
        let buf = [0.0f32; PIT_BUFFER];
        let res = open_loop_pitch(&buf);
        assert!(res.r_norm == 0.0);
        assert!((PIT_MIN..=PIT_MAX).contains(&res.t_op));
    }

    /// The returned delay is always inside the legal open-loop range.
    #[test]
    fn delay_always_in_range() {
        for period in [21, 37, 45, 70, 90, 130] {
            let buf = periodic_buf(period);
            let res = open_loop_pitch(&buf);
            assert!(
                (PIT_MIN..=PIT_MAX).contains(&res.t_op),
                "period {period}: T_op {} out of range",
                res.t_op
            );
        }
    }
}
