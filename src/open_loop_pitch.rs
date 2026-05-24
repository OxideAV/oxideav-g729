//! ITU-T G.729 §3.4 open-loop pitch analysis.
//!
//! The open-loop pitch search is performed once per 10 ms frame on the
//! **perceptually-weighted speech signal** `sw(n)` of equation (33). Its
//! output `T_op` is the candidate delay that anchors the per-subframe
//! closed-loop adaptive-codebook search of §3.7: subframe 1 searches a
//! six-sample window around `T_op`, subframe 2 searches around the
//! subframe-1 result.
//!
//! ## Equation (33) — weighted speech `sw(n)`
//!
//! For each subframe (`n = 0..39`):
//!
//! ```text
//!   sw(n) = s(n) + Σ_{i=1..10} a_i·γ1^i · s(n-i)
//!                − Σ_{i=1..10} a_i·γ2^i · sw(n-i)
//! ```
//!
//! This is the perceptual weighting filter `W(z) = A(z/γ1)/A(z/γ2)`
//! applied to the (pre-emphasised) input speech, with `(a_i)` the
//! **unquantised** LP coefficients of the current frame's analysis and
//! `(γ1, γ2)` the adapted weighting gammas from §3.3 (see
//! [`crate::weighting::weighting_gammas`]). The filter is run with
//! persistent state across subframes and frames: the encoder retains
//! enough history of past `s` and `sw` to evaluate the recursion at the
//! start of every new subframe.
//!
//! ## Equations (34) + (35) — pitch correlation
//!
//! For each candidate integer delay `k`:
//!
//! ```text
//!   R(k)  = Σ_{n=0..79} sw(n) · sw(n-k)             (eq 34)
//!   R'(k) = R(k) / sqrt( Σ_{n=0..79} sw²(n-k) )     (eq 35)
//! ```
//!
//! `R(k)` is a **frame-level** (80-sample) correlation despite eq (33)
//! defining `sw` on a subframe basis: by the time the open-loop search
//! runs, both subframes of `sw(n)` for the current frame have been
//! produced and are concatenated with the trailing portion of the
//! previous frame's `sw` to form the lag-up-to-143 window.
//!
//! ## Three-range search + sub-multiple bias
//!
//! Three local maxima are found over the three disjoint delay ranges:
//!
//! ```text
//!   i = 1:  80 … 143
//!   i = 2:  40 …  79
//!   i = 3:  20 …  39
//! ```
//!
//! (The spec text mislabels the third row "i = 1"; the surrounding
//! prose makes the intended `i = 3` unambiguous.)
//!
//! The winner is selected with a sub-multiple bias (spec §3.4 prose,
//! "favouring the delays with the values in the lower range"):
//!
//! ```text
//!   T_op = t_1; R'(T_op) = R'(t_1)
//!   if R'(t_2) ≥ 0.85 · R'(T_op):
//!       R'(T_op) = R'(t_2);  T_op = t_2
//!   if R'(t_3) ≥ 0.85 · R'(T_op):
//!       R'(T_op) = R'(t_3);  T_op = t_3
//! ```
//!
//! Lower-range candidates only win when their correlation is within 15 %
//! of the current best, which prevents the search from latching onto a
//! pitch *multiple* (an integer doubling of the true period whose
//! correlation will look just as strong as the fundamental's).

use crate::{LPC_ORDER, SUBFRAME_SAMPLES};

/// Total samples per frame on which `R(k)` (eq 34) sums.
pub const PITCH_FRAME_SAMPLES: usize = 2 * SUBFRAME_SAMPLES; // 80

/// Maximum integer pitch delay searched by §3.4 (and used by the
/// adaptive-codebook in §3.7). Spec range is `[20, 143]`.
pub const MAX_PITCH_LAG: usize = 143;

/// Minimum integer pitch delay searched by §3.4. Spec range is
/// `[20, 143]`.
pub const MIN_PITCH_LAG: usize = 20;

/// Sub-multiple bias threshold from the §3.4 selection tree. A
/// shorter-range candidate `t_i` (i > 1) only displaces the current
/// winner `T_op` when `R'(t_i) ≥ BIAS · R'(T_op)`. The spec uses
/// `BIAS = 0.85`.
pub const SUBMULTIPLE_BIAS: f32 = 0.85;

/// Persistent state for the §3.3/§3.4 weighted-speech path. Holds
/// enough history of `s(n)` and `sw(n)` to evaluate equation (33)
/// continuously across subframe and frame boundaries, plus enough
/// history of `sw(n)` to evaluate equation (34) at lags up to
/// [`MAX_PITCH_LAG`].
///
/// All buffers store newest sample at the *high* end. Each call to
/// [`Self::run_subframe`] shifts the buffer left by
/// [`SUBFRAME_SAMPLES`] and appends the new subframe at the tail.
pub struct WeightedSpeechState {
    /// History of input speech `s(n)`. Need the latest [`LPC_ORDER`]
    /// samples to evaluate eq (33) at the start of a new subframe;
    /// the size is rounded up to one subframe for clean roll-in/out.
    s_hist: [f32; LPC_ORDER],
    /// History of weighted speech `sw(n)`. Need at most
    /// `LPC_ORDER` samples for the eq (33) recursion **and** up to
    /// [`MAX_PITCH_LAG`] samples for the eq (34) correlation window.
    /// Sized to `MAX_PITCH_LAG` because that's the larger requirement.
    sw_hist: [f32; MAX_PITCH_LAG],
}

impl Default for WeightedSpeechState {
    fn default() -> Self {
        Self::new()
    }
}

impl WeightedSpeechState {
    /// Construct a state with all-zero history.
    pub fn new() -> Self {
        Self {
            s_hist: [0.0; LPC_ORDER],
            sw_hist: [0.0; MAX_PITCH_LAG],
        }
    }

    /// Apply equation (33) to a single subframe of pre-emphasised input
    /// speech `s_sub` (length [`SUBFRAME_SAMPLES`]) using the supplied
    /// unquantised per-subframe LP coefficients `a` (length
    /// `LPC_ORDER + 1`) and the perceptual-weighting gammas
    /// `(gamma1, gamma2)`. Writes the resulting weighted speech samples
    /// into `sw_out` and advances the internal history.
    ///
    /// Implementation detail: the spec writes the recursion in terms of
    /// `a_i · γ1^i` and `a_i · γ2^i`. These per-frame coefficients are
    /// the bandwidth-expanded LP filters
    /// `A(z/γ1) = 1 + Σ a_i·γ1^i·z^-i` and
    /// `A(z/γ2) = 1 + Σ a_i·γ2^i·z^-i`. To avoid recomputing the
    /// γ^i powers per sample, the caller passes pre-expanded arrays —
    /// see [`Self::run_subframe`].
    pub fn run_subframe_expanded(
        &mut self,
        s_sub: &[f32; SUBFRAME_SAMPLES],
        aw1: &[f32; LPC_ORDER + 1],
        aw2: &[f32; LPC_ORDER + 1],
        sw_out: &mut [f32; SUBFRAME_SAMPLES],
    ) {
        // Evaluate the recursion sample-by-sample. For each `n` we need
        // s(n-i) and sw(n-i) for i = 1..LPC_ORDER. Indices with n-i < 0
        // come from `s_hist` / `sw_hist`; indices with n-i >= 0 come
        // from the current subframe (`s_sub` / `sw_out`).
        for n in 0..SUBFRAME_SAMPLES {
            let s_n = s_sub[n];
            // Σ_{i=1..10} aw1[i] · s(n-i)
            let mut acc_num = s_n;
            for i in 1..=LPC_ORDER {
                let s_im = if (n as isize) - (i as isize) >= 0 {
                    s_sub[n - i]
                } else {
                    // s_hist holds the latest LPC_ORDER input samples;
                    // index `LPC_ORDER - 1` is the most recent.
                    let off = (i as isize) - (n as isize) - 1;
                    if (off as usize) < LPC_ORDER {
                        self.s_hist[LPC_ORDER - 1 - off as usize]
                    } else {
                        0.0
                    }
                };
                acc_num += aw1[i] * s_im;
            }
            // Σ_{i=1..10} aw2[i] · sw(n-i)
            let mut acc_den = 0.0f32;
            for i in 1..=LPC_ORDER {
                let sw_im = if (n as isize) - (i as isize) >= 0 {
                    sw_out[n - i]
                } else {
                    let off = (i as isize) - (n as isize) - 1;
                    // sw_hist's most recent sample sits at index
                    // MAX_PITCH_LAG - 1.
                    if (off as usize) < MAX_PITCH_LAG {
                        self.sw_hist[MAX_PITCH_LAG - 1 - off as usize]
                    } else {
                        0.0
                    }
                };
                acc_den += aw2[i] * sw_im;
            }
            sw_out[n] = acc_num - acc_den;
        }

        // Roll history: append the new subframe at the tail. `s_hist`
        // only retains the last LPC_ORDER samples, so the slide is
        // effectively "keep the last LPC_ORDER samples of s_sub".
        // `sw_hist` retains MAX_PITCH_LAG samples.
        // -- s_hist: keep the most recent LPC_ORDER input samples --
        self.s_hist
            .copy_from_slice(&s_sub[SUBFRAME_SAMPLES - LPC_ORDER..]);
        // -- sw_hist: shift left by SUBFRAME_SAMPLES, then append sw_out --
        let keep = MAX_PITCH_LAG - SUBFRAME_SAMPLES;
        self.sw_hist.copy_within(SUBFRAME_SAMPLES.., 0);
        self.sw_hist[keep..].copy_from_slice(sw_out);
    }

    /// Convenience wrapper around [`Self::run_subframe_expanded`] that
    /// derives `A(z/γ1)` and `A(z/γ2)` from the raw LP coefficients
    /// `a` and the gamma pair, then runs the subframe.
    pub fn run_subframe(
        &mut self,
        s_sub: &[f32; SUBFRAME_SAMPLES],
        a: &[f32; LPC_ORDER + 1],
        gamma1: f32,
        gamma2: f32,
        sw_out: &mut [f32; SUBFRAME_SAMPLES],
    ) {
        let aw1 = crate::weighting::bandwidth_expand(a, gamma1);
        let aw2 = crate::weighting::bandwidth_expand(a, gamma2);
        self.run_subframe_expanded(s_sub, &aw1, &aw2, sw_out);
    }

    /// Borrow the up-to-[`MAX_PITCH_LAG`] history of past weighted-speech
    /// samples retained by the state. Newest sample at the high end.
    pub fn sw_history(&self) -> &[f32; MAX_PITCH_LAG] {
        &self.sw_hist
    }
}

/// Result of [`open_loop_pitch_search`]: the chosen integer pitch
/// delay `T_op` and the normalised correlation `R'(T_op)` at that
/// delay. The latter is informational (useful for downstream
/// confidence-of-pitch heuristics); the search itself uses it only
/// internally.
#[derive(Clone, Copy, Debug, Default)]
pub struct OpenLoopPitch {
    /// Chosen integer pitch delay, in samples.
    pub t_op: usize,
    /// Normalised correlation `R'(T_op)` (eq 35).
    pub r_norm: f32,
}

/// Search ranges for the three-section open-loop pitch search of §3.4.
/// Order matters: ranges are processed `i = 1 → 2 → 3`, so the highest
/// (`80..=143`) is evaluated first.
const PITCH_RANGES: [(usize, usize); 3] = [(80, 143), (40, 79), (20, 39)];

/// Open-loop pitch search per ITU-T G.729 §3.4 (equations 34, 35 + the
/// sub-multiple bias decision tree).
///
/// * `sw_frame` is the current frame's weighted speech, `sw(0..=79)`.
/// * `sw_hist` holds the trailing portion of the previous frame's
///   weighted speech, with the newest sample at the high end (index
///   `MAX_PITCH_LAG - 1`). This is what
///   [`WeightedSpeechState::sw_history`] returns.
///
/// Conceptually `sw_full(n)` is the concatenation `[ sw_hist | sw_frame ]`,
/// indexed so that the latest sample is `sw_full(79)` (== `sw_frame[79]`).
/// Equations (34) and (35) index this stitched buffer at lags 1..=143
/// from the current frame.
pub fn open_loop_pitch_search(
    sw_frame: &[f32; PITCH_FRAME_SAMPLES],
    sw_hist: &[f32; MAX_PITCH_LAG],
) -> OpenLoopPitch {
    // Stitch sw_hist + sw_frame so we can index uniformly.
    let mut stitched = [0.0f32; MAX_PITCH_LAG + PITCH_FRAME_SAMPLES];
    stitched[..MAX_PITCH_LAG].copy_from_slice(sw_hist);
    stitched[MAX_PITCH_LAG..].copy_from_slice(sw_frame);
    // The current-frame indices 0..=79 map onto stitched indices
    // [MAX_PITCH_LAG..MAX_PITCH_LAG + 80].
    let base = MAX_PITCH_LAG;

    // For each range, find the integer lag `t` maximising the *signed*
    // normalised correlation R'(t) = R(t) / sqrt(Σ sw²(n-t)). We track
    // the maximum signed value (per spec: the open-loop search retains
    // the algebraic sign — a negative R(t) means the candidate is
    // anti-correlated and gets rejected against the i=1 baseline).
    let mut t_per_range = [0usize; 3];
    let mut r_per_range = [f32::NEG_INFINITY; 3];

    for (range_idx, &(lo, hi)) in PITCH_RANGES.iter().enumerate() {
        let mut best_t = lo;
        let mut best_r = f32::NEG_INFINITY;
        for t in lo..=hi {
            // R(t) and Σ sw²(n-t) over n = 0..=79.
            let mut num = 0.0f32;
            let mut den = 0.0f32;
            for n in 0..PITCH_FRAME_SAMPLES {
                let sw_n = stitched[base + n];
                let sw_nt = stitched[base + n - t]; // base + n >= base; n - t may be negative within stitched but stays >= 0 by construction (t <= MAX_PITCH_LAG, base = MAX_PITCH_LAG, so base + n - t >= 0).
                num += sw_n * sw_nt;
                den += sw_nt * sw_nt;
            }
            let denom = den.sqrt();
            let r_norm = if denom > 1e-6 { num / denom } else { 0.0 };
            if r_norm > best_r {
                best_r = r_norm;
                best_t = t;
            }
        }
        t_per_range[range_idx] = best_t;
        r_per_range[range_idx] = best_r;
    }

    // Sub-multiple bias decision tree (§3.4 prose).
    let mut t_op = t_per_range[0];
    let mut r_top = r_per_range[0];
    if r_per_range[1] >= SUBMULTIPLE_BIAS * r_top {
        t_op = t_per_range[1];
        r_top = r_per_range[1];
    }
    if r_per_range[2] >= SUBMULTIPLE_BIAS * r_top {
        t_op = t_per_range[2];
        r_top = r_per_range[2];
    }
    OpenLoopPitch {
        t_op,
        r_norm: r_top,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::weighting::bandwidth_expand;

    /// A trivial LP filter A(z) = 1 (no past memory contribution). With
    /// this filter the weighting recursion reduces to `sw(n) = s(n)` —
    /// the perceptual filter is a unit-gain all-pass.
    fn trivial_a() -> [f32; LPC_ORDER + 1] {
        let mut a = [0.0f32; LPC_ORDER + 1];
        a[0] = 1.0;
        a
    }

    #[test]
    fn run_subframe_trivial_filter_is_identity() {
        // With A(z) = 1, sw(n) = s(n) regardless of γ1, γ2.
        let mut state = WeightedSpeechState::new();
        let a = trivial_a();
        let mut s = [0.0f32; SUBFRAME_SAMPLES];
        for n in 0..SUBFRAME_SAMPLES {
            s[n] = (n as f32 * 0.1).sin();
        }
        let mut sw = [0.0f32; SUBFRAME_SAMPLES];
        state.run_subframe(&s, &a, 0.94, 0.6, &mut sw);
        for n in 0..SUBFRAME_SAMPLES {
            assert!(
                (sw[n] - s[n]).abs() < 1e-6,
                "sw[{n}] {} != s[{n}] {}",
                sw[n],
                s[n]
            );
        }
    }

    #[test]
    fn run_subframe_continuous_state_across_subframes() {
        // Running two contiguous subframes through `run_subframe` must
        // produce the same output as running an 80-sample block through
        // an independent reference implementation — verifies that the
        // state hand-off (s_hist + sw_hist) is correct.
        let mut a = [0.0f32; LPC_ORDER + 1];
        a[0] = 1.0;
        for k in 1..=LPC_ORDER {
            a[k] = -0.05 * k as f32;
        }
        let g1 = 0.98f32;
        let g2 = 0.55f32;

        // Block input.
        let mut s_full = [0.0f32; PITCH_FRAME_SAMPLES];
        for n in 0..PITCH_FRAME_SAMPLES {
            s_full[n] = ((n as f32) * 0.07).sin() + 0.3 * ((n as f32) * 0.21).cos();
        }

        // ---- Path A: two subframes through the stateful API.
        let mut state = WeightedSpeechState::new();
        let mut sub0 = [0.0f32; SUBFRAME_SAMPLES];
        sub0.copy_from_slice(&s_full[..SUBFRAME_SAMPLES]);
        let mut sub1 = [0.0f32; SUBFRAME_SAMPLES];
        sub1.copy_from_slice(&s_full[SUBFRAME_SAMPLES..]);
        let mut sw0 = [0.0f32; SUBFRAME_SAMPLES];
        let mut sw1 = [0.0f32; SUBFRAME_SAMPLES];
        state.run_subframe(&sub0, &a, g1, g2, &mut sw0);
        state.run_subframe(&sub1, &a, g1, g2, &mut sw1);
        let mut sw_a = [0.0f32; PITCH_FRAME_SAMPLES];
        sw_a[..SUBFRAME_SAMPLES].copy_from_slice(&sw0);
        sw_a[SUBFRAME_SAMPLES..].copy_from_slice(&sw1);

        // ---- Path B: independent reference recursion over the full 80
        // samples with zero initial state.
        let aw1 = bandwidth_expand(&a, g1);
        let aw2 = bandwidth_expand(&a, g2);
        let mut sw_b = [0.0f32; PITCH_FRAME_SAMPLES];
        for n in 0..PITCH_FRAME_SAMPLES {
            let mut num = s_full[n];
            for i in 1..=LPC_ORDER {
                if n >= i {
                    num += aw1[i] * s_full[n - i];
                }
            }
            let mut den = 0.0f32;
            for i in 1..=LPC_ORDER {
                if n >= i {
                    den += aw2[i] * sw_b[n - i];
                }
            }
            sw_b[n] = num - den;
        }

        for n in 0..PITCH_FRAME_SAMPLES {
            assert!(
                (sw_a[n] - sw_b[n]).abs() < 1e-4,
                "sw mismatch at n={n}: stateful={} reference={}",
                sw_a[n],
                sw_b[n]
            );
        }
    }

    /// Build a synthetic weighted-speech buffer with a known periodic
    /// structure, then verify the open-loop search picks the period.
    fn build_periodic(period: usize, n: usize) -> Vec<f32> {
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            // Periodic pulse train with period `period` plus a small
            // smooth envelope so the correlation is finite everywhere.
            let phase = (i % period) as f32 / period as f32;
            let pulse = (-((phase - 0.0) * (phase - 0.0) * 50.0)).exp();
            out.push(pulse + 0.01 * ((i as f32) * 0.3).sin());
        }
        out
    }

    #[test]
    fn open_loop_picks_period_in_high_range() {
        // Period 100 (in the 80..=143 range): all three normalised
        // correlations should peak near multiples of 100, but only
        // t=100 is in range 1; the other ranges peak elsewhere with a
        // lower correlation. Sub-multiple bias should NOT flip the
        // answer since the lower ranges aren't within 15 %.
        let period = 100;
        let total = MAX_PITCH_LAG + PITCH_FRAME_SAMPLES;
        let buf = build_periodic(period, total);
        let mut hist = [0.0f32; MAX_PITCH_LAG];
        hist.copy_from_slice(&buf[..MAX_PITCH_LAG]);
        let mut frame = [0.0f32; PITCH_FRAME_SAMPLES];
        frame.copy_from_slice(&buf[MAX_PITCH_LAG..]);
        let pitch = open_loop_pitch_search(&frame, &hist);
        assert!(
            (pitch.t_op as isize - period as isize).abs() <= 2,
            "expected T_op ~{period}, got {} (R'={:.3})",
            pitch.t_op,
            pitch.r_norm
        );
    }

    #[test]
    fn open_loop_picks_period_in_mid_range() {
        // Period 60 (in 40..=79). Range-1 (80..143) will see a multiple
        // (120) with correlation near the period-60 correlation; the
        // 0.85 sub-multiple bias must steer the result to range 2.
        let period = 60;
        let total = MAX_PITCH_LAG + PITCH_FRAME_SAMPLES;
        let buf = build_periodic(period, total);
        let mut hist = [0.0f32; MAX_PITCH_LAG];
        hist.copy_from_slice(&buf[..MAX_PITCH_LAG]);
        let mut frame = [0.0f32; PITCH_FRAME_SAMPLES];
        frame.copy_from_slice(&buf[MAX_PITCH_LAG..]);
        let pitch = open_loop_pitch_search(&frame, &hist);
        assert!(
            (40..=79).contains(&pitch.t_op),
            "expected T_op in mid range, got {} (R'={:.3})",
            pitch.t_op,
            pitch.r_norm
        );
        assert!(
            (pitch.t_op as isize - period as isize).abs() <= 2,
            "expected T_op ~{period}, got {}",
            pitch.t_op
        );
    }

    #[test]
    fn open_loop_picks_period_in_low_range() {
        // Period 25 (in 20..=39). Higher ranges will see multiples at
        // 50, 75, 100, 125 — sub-multiple bias must steer to range 3.
        let period = 25;
        let total = MAX_PITCH_LAG + PITCH_FRAME_SAMPLES;
        let buf = build_periodic(period, total);
        let mut hist = [0.0f32; MAX_PITCH_LAG];
        hist.copy_from_slice(&buf[..MAX_PITCH_LAG]);
        let mut frame = [0.0f32; PITCH_FRAME_SAMPLES];
        frame.copy_from_slice(&buf[MAX_PITCH_LAG..]);
        let pitch = open_loop_pitch_search(&frame, &hist);
        assert!(
            (20..=39).contains(&pitch.t_op),
            "expected T_op in low range, got {}",
            pitch.t_op
        );
        assert!(
            (pitch.t_op as isize - period as isize).abs() <= 2,
            "expected T_op ~{period}, got {}",
            pitch.t_op
        );
    }

    #[test]
    fn open_loop_zero_input_is_finite() {
        let frame = [0.0f32; PITCH_FRAME_SAMPLES];
        let hist = [0.0f32; MAX_PITCH_LAG];
        let pitch = open_loop_pitch_search(&frame, &hist);
        // With all-zero input every R'(k) is 0, so the i=1 baseline
        // wins by default at the bottom of range 1 (t=80). The exact
        // tie-break value doesn't matter; we only require the search
        // to terminate with a valid in-range delay and a finite
        // correlation.
        assert!((MIN_PITCH_LAG..=MAX_PITCH_LAG).contains(&pitch.t_op));
        assert!(pitch.r_norm.is_finite());
    }

    #[test]
    fn submultiple_bias_constant_value() {
        // Sanity-check the constant pulled from §3.4 prose: 0.85.
        // Catches a silent value drift in future edits.
        assert!((SUBMULTIPLE_BIAS - 0.85).abs() < 1e-6);
    }

    #[test]
    fn pitch_ranges_cover_spec_window() {
        // The three ranges, concatenated, must cover [20, 143] exactly
        // once with no gaps and no overlap — eq-34 search window per
        // §3.4.
        let mut covered = [false; MAX_PITCH_LAG + 1];
        for &(lo, hi) in PITCH_RANGES.iter() {
            for t in lo..=hi {
                assert!(!covered[t], "range overlap at t={t}");
                covered[t] = true;
            }
        }
        for t in MIN_PITCH_LAG..=MAX_PITCH_LAG {
            assert!(covered[t], "lag t={t} not covered by any range");
        }
        for t in 0..MIN_PITCH_LAG {
            assert!(!covered[t], "lag t={t} should not be covered");
        }
    }
}
