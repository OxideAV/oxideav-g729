//! §3.7 **closed-loop pitch (adaptive-codebook) search** — finds the
//! per-subframe fractional pitch delay `T1` / `T2` that maximises the
//! weighted-correlation criterion, searching a small window around the
//! §3.4 open-loop delay (subframe 1) or around `int(T1)` (subframe 2).
//!
//! Ninth encoder-side stage. Consumes the §3.6 target `x(n)`, the §3.5
//! impulse response `h(n)`, and the excitation buffer `u(n)` (past
//! excitation extended by the LP residual, per §3.7: "in the search
//! stage, the samples `u(n), n = 0 … 39` are not known … the LP
//! residual is copied to `u(n)` to make the relation in equation (38)
//! valid for all delays").
//!
//! ## Spec source — clause 3.7, equations (37)–(39) + window procedures
//!
//! (Transcribed from the EPUB prose + rasters `f0018-01/02.jpg`,
//! `eq37/38/39.jpg`.)
//!
//! * **Subframe-1 window** (`f0018-01`, around the open-loop `T_op`):
//!   ```text
//!   t_min = T_op − 3;  if t_min < 20 { t_min = 20 }
//!   t_max = t_min + 6
//!   if t_max > 143 { t_max = 143; t_min = t_max − 6 }
//!   ```
//! * **Subframe-2 window** (`f0018-02`, around `int(T1)`):
//!   ```text
//!   t_min = int(T1) − 5;  if t_min < 20 { t_min = 20 }
//!   t_max = t_min + 9
//!   if t_max > 143 { t_max = 143; t_min = t_max − 9 }
//!   ```
//!   (searched at 1/3 resolution between `t_min − 2/3` and
//!   `t_max + 2/3`).
//! * **eq (37)** — the criterion:
//!   `R(k) = Σ_{n=0}^{39} x(n)·y_k(n) / √(Σ_{n=0}^{39} y_k(n)·y_k(n))`
//!   where `y_k(n)` is the past excitation at delay `k` convolved with
//!   `h(n)`.
//! * **eq (38)** — the delay recursion:
//!   `y_k(n) = y_{k−1}(n−1) + u(−k)·h(n)`, `n = 39 … 0`, with
//!   `y_{k−1}(−1) = 0` — each next (longer) delay shifts the previous
//!   filtered excitation and adds one impulse-response contribution.
//! * **eq (39)** — fractional refinement by interpolating the
//!   normalised correlation with the 1/3-resolution FIR `b12` (Hamming
//!   windowed sinc, ±11, `b12(12) = 0`, −3 dB at 3600 Hz):
//!   `R(k)_t = Σ_{i=0}^{3} R(k−i)·b12(t+3i) + Σ_{i=0}^{3} R(k+1+i)·b12(3−t+3i)`,
//!   `t = 0, 1, 2` ↔ fractions `0, 1/3, 2/3`. The correlations are
//!   computed over `t_min − 4 … t_max + 4` to feed the interpolator.
//! * **Fraction policy** (clause 3.7): for `T1`, fractions are tested
//!   only when the optimum integer delay is `< 85` (delays ≥ 85 are
//!   integer-only); for `T2` fractions are always tested.
//!
//! The transmitted-index encode (`P1`/`P2`, eqs (41)/(42a-b)) already
//! exists as [`crate::pitch_decode::encode_p1`] /
//! [`crate::pitch_decode::encode_p2`]; this module returns
//! [`crate::pitch_decode::PitchDelay`] values compatible with them.

use crate::pitch_decode::PitchDelay;
use crate::tables::PITCH_INTERP_FILTER_ANALYSIS_Q15;

/// Samples per subframe (clause 2.1).
pub const L_SUBFR: usize = 40;
/// Minimum pitch delay (clause 3.7).
pub const PIT_MIN: usize = 20;
/// Maximum pitch delay (clause 3.7).
pub const PIT_MAX: usize = 143;
/// The excitation buffer layout: `u(n)` for `n = −143 … 39` = 183
/// samples; index `PIT_MAX + n` addresses `u(n)`.
pub const EXC_BUFFER: usize = PIT_MAX + L_SUBFR;

/// Q15 scale of the `b12` interpolation filter.
const Q15: f32 = 32768.0;

/// Integer-only threshold for `T1` (clause 3.7: fractional resolution
/// applies below 85; `[85, 143]` is integer-only).
pub const T1_FRACTIONAL_LIMIT: i32 = 85;

/// Result of one subframe's closed-loop search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosedLoopResult {
    /// The selected fractional pitch delay.
    pub delay: PitchDelay,
    /// The eq (37) normalised correlation at the selected delay (the
    /// interpolated value when a non-zero fraction won).
    pub score: f32,
}

/// The §3.7 subframe-1 search window `[t_min, t_max]` around the
/// open-loop delay `T_op` (spec procedure `f0018-01`).
#[must_use]
pub fn subframe1_window(t_op: i32) -> (i32, i32) {
    let mut t_min = t_op - 3;
    if t_min < PIT_MIN as i32 {
        t_min = PIT_MIN as i32;
    }
    let mut t_max = t_min + 6;
    if t_max > PIT_MAX as i32 {
        t_max = PIT_MAX as i32;
        t_min = t_max - 6;
    }
    (t_min, t_max)
}

/// The §3.7 subframe-2 search window `[t_min, t_max]` around `int(T1)`
/// (spec procedure `f0018-02`).
#[must_use]
pub fn subframe2_window(int_t1: i32) -> (i32, i32) {
    let mut t_min = int_t1 - 5;
    if t_min < PIT_MIN as i32 {
        t_min = PIT_MIN as i32;
    }
    let mut t_max = t_min + 9;
    if t_max > PIT_MAX as i32 {
        t_max = PIT_MAX as i32;
        t_min = t_max - 9;
    }
    (t_min, t_max)
}

/// eq (38) driver: computes the filtered past excitation `y_k(n)` for
/// every integer delay `k ∈ [k_lo, k_hi]`, returning them as rows of a
/// `Vec` (row 0 ↔ `k_lo`). `exc` is the `u(n)` buffer (`u(−143) …
/// u(39)`, index `PIT_MAX + n`), `h` the §3.5 impulse response.
///
/// The first delay is computed as a direct convolution; each subsequent
/// delay uses the eq (38) shift-and-add recursion.
fn filtered_excitation_rows(
    exc: &[f32; EXC_BUFFER],
    h: &[f32; L_SUBFR],
    k_lo: usize,
    k_hi: usize,
) -> Vec<[f32; L_SUBFR]> {
    let mut rows = Vec::with_capacity(k_hi - k_lo + 1);

    // Direct convolution for k_lo:
    // y(n) = Σ_{j=0}^{n} u(j − k_lo)·h(n − j).
    let mut y = [0.0f32; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = 0.0f32;
        for j in 0..=n {
            acc += exc[PIT_MAX + j - k_lo] * h[n - j];
        }
        y[n] = acc;
    }
    rows.push(y);

    // eq (38): y_k(n) = y_{k−1}(n−1) + u(−k)·h(n), n = 39 … 0.
    for k in (k_lo + 1)..=k_hi {
        let mut next = [0.0f32; L_SUBFR];
        let u_mk = exc[PIT_MAX - k];
        for n in (0..L_SUBFR).rev() {
            let shifted = if n > 0 { y[n - 1] } else { 0.0 };
            next[n] = shifted + u_mk * h[n];
        }
        y = next;
        rows.push(y);
    }
    rows
}

/// eq (37) at one integer delay row.
fn normalised_correlation(x: &[f32; L_SUBFR], y: &[f32; L_SUBFR]) -> f32 {
    let mut num = 0.0f32;
    let mut den = 0.0f32;
    for n in 0..L_SUBFR {
        num += x[n] * y[n];
        den += y[n] * y[n];
    }
    if den > 0.0 {
        num / den.sqrt()
    } else {
        0.0
    }
}

/// `b12(i)` — the eq (39) interpolation filter read symmetrically from
/// the compiled 13-entry one-sided table (`|i| ≤ 12`; `b12(12) = 0`).
fn b12(i: i32) -> f32 {
    let idx = i.unsigned_abs() as usize;
    if idx >= PITCH_INTERP_FILTER_ANALYSIS_Q15.len() {
        0.0
    } else {
        f32::from(PITCH_INTERP_FILTER_ANALYSIS_Q15[idx]) / Q15
    }
}

/// eq (39): interpolates the normalised-correlation sequence `r` (whose
/// index `i` holds `R(k_base + i)`) at delay `k` and fraction `t/3`
/// (`t = 0, 1, 2`).
fn interpolate_correlation(r: &[f32], k_base: i32, k: i32, t: i32) -> f32 {
    let at = |kk: i32| -> f32 {
        let idx = kk - k_base;
        if idx < 0 || idx as usize >= r.len() {
            0.0
        } else {
            r[idx as usize]
        }
    };
    let mut acc = 0.0f32;
    for i in 0..4 {
        acc += at(k - i) * b12(t + 3 * i) + at(k + 1 + i) * b12(3 - t + 3 * i);
    }
    acc
}

/// §3.7 closed-loop search over one subframe.
///
/// * `x` — the §3.6 target signal.
/// * `h` — the §3.5 impulse response.
/// * `exc` — the excitation buffer `u(−143) … u(39)` with the current
///   subframe region filled with the LP residual (per clause 3.7).
/// * `(t_min, t_max)` — the integer search window
///   ([`subframe1_window`] / [`subframe2_window`]).
/// * `allow_fractions` — pass `false` only when the caller knows the
///   whole window is integer-only; when `true` the fraction test is
///   still skipped for `T1` delays `≥ 85` per the clause-3.7 rule (the
///   caller communicates that by the window; this function applies the
///   rule via `frac_limit`).
/// * `frac_limit` — delays `≥ frac_limit` are searched integer-only
///   (pass [`T1_FRACTIONAL_LIMIT`] for subframe 1, `i32::MAX` — no
///   limit — for subframe 2).
#[must_use]
pub fn closed_loop_search(
    x: &[f32; L_SUBFR],
    h: &[f32; L_SUBFR],
    exc: &[f32; EXC_BUFFER],
    t_min: i32,
    t_max: i32,
    frac_limit: i32,
) -> ClosedLoopResult {
    // Correlations over the widened range t_min−4 … t_max+4 feed the
    // eq (39) interpolator; delays outside [PIT_MIN, PIT_MAX] are
    // clamped out of the correlation window (eq (39) sees 0 there —
    // the interpolation only ever evaluates fractions of in-window
    // delays, whose ±4 neighbourhood stays within the buffer's reach
    // because u() extends back to −143 and t_max ≤ 143 − … guarded by
    // the window procedures).
    let k_lo = (t_min - 4).max(PIT_MIN as i32 - 4).max(1) as usize;
    let k_hi = (t_max + 4).min(PIT_MAX as i32) as usize;

    let rows = filtered_excitation_rows(exc, h, k_lo, k_hi);
    let r: Vec<f32> = rows.iter().map(|y| normalised_correlation(x, y)).collect();
    let k_base = k_lo as i32;
    let r_at = |k: i32| -> f32 {
        let idx = k - k_base;
        if idx < 0 || idx as usize >= r.len() {
            f32::NEG_INFINITY
        } else {
            r[idx as usize]
        }
    };

    // Integer pass over the actual window.
    let mut best_k = t_min;
    let mut best_r = f32::NEG_INFINITY;
    for k in t_min..=t_max {
        let v = r_at(k);
        if v > best_r {
            best_r = v;
            best_k = k;
        }
    }

    // Fractional pass (clause 3.7): only when the optimum integer delay
    // is below the integer-only threshold.
    let mut best = ClosedLoopResult {
        delay: PitchDelay {
            int_t: best_k,
            frac: 0,
        },
        score: best_r,
    };
    if best_k < frac_limit {
        // Candidate fractions around best_k: −1/3 (= fraction 2/3 of
        // k−1) and +1/3, +2/3 of k itself; the +2/3 case re-expresses
        // as −1/3 of k+1 in the {−1, 0, +1} transmission convention.
        // Test t = 1, 2 at k−1 (fractions −2/3, −1/3) and t = 1, 2 at
        // k (fractions +1/3, +2/3), mapping each winner back into the
        // PitchDelay {int_t, frac ∈ {−1, 0, 1}} convention.
        let candidates: [(i32, i32, i32, i32); 4] = [
            // (eval_k, eval_t, out_int, out_frac)
            (best_k - 1, 1, best_k - 1, 1), // k − 2/3 = (k−1) + 1/3
            (best_k - 1, 2, best_k, -1),    // k − 1/3 = (k−1) + 2/3
            (best_k, 1, best_k, 1),         // k + 1/3
            (best_k, 2, best_k + 1, -1),    // k + 2/3 = (k+1) − 1/3
        ];
        for &(ek, et, oi, of) in &candidates {
            if ek < PIT_MIN as i32 || oi > PIT_MAX as i32 || oi < PIT_MIN as i32 {
                continue;
            }
            let v = interpolate_correlation(&r, k_base, ek, et);
            if v > best.score {
                best = ClosedLoopResult {
                    delay: PitchDelay {
                        int_t: oi,
                        frac: of,
                    },
                    score: v,
                };
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Windows: interior, low clamp, and high clamp of both procedures.
    #[test]
    fn search_windows() {
        // Subframe 1: interior.
        assert_eq!(subframe1_window(60), (57, 63));
        // Low clamp.
        assert_eq!(subframe1_window(20), (20, 26));
        // High clamp.
        assert_eq!(subframe1_window(143), (137, 143));

        // Subframe 2: interior.
        assert_eq!(subframe2_window(60), (55, 64));
        // Low clamp.
        assert_eq!(subframe2_window(21), (20, 29));
        // High clamp.
        assert_eq!(subframe2_window(143), (134, 143));
    }

    /// eq (38) recursion equals the direct convolution at every delay
    /// in a window.
    #[test]
    fn eq38_recursion_matches_direct() {
        let exc: [f32; EXC_BUFFER] = std::array::from_fn(|n| ((n * 13 % 97) as f32) - 48.0);
        let h: [f32; L_SUBFR] =
            std::array::from_fn(|n| if n == 0 { 1.0 } else { 0.8f32.powi(n as i32) });
        let (k_lo, k_hi) = (30usize, 45usize);
        let rows = filtered_excitation_rows(&exc, &h, k_lo, k_hi);
        for (row, k) in rows.iter().zip(k_lo..=k_hi) {
            for n in 0..L_SUBFR {
                let mut want = 0.0f32;
                for j in 0..=n {
                    want += exc[PIT_MAX + j - k] * h[n - j];
                }
                assert!(
                    (row[n] - want).abs() < 1e-2 * (1.0 + want.abs()),
                    "k={k} n={n}: {} vs {}",
                    row[n],
                    want
                );
            }
        }
    }

    /// A synthetic periodic excitation: the target is the filtered
    /// excitation at a known integer delay; the search recovers exactly
    /// that delay with fraction 0.
    #[test]
    fn recovers_known_integer_delay() {
        let period = 58usize;
        let exc: [f32; EXC_BUFFER] = std::array::from_fn(|n| {
            let ph = (n % period) as f32 / period as f32;
            (std::f32::consts::TAU * ph).sin()
        });
        let h: [f32; L_SUBFR] =
            std::array::from_fn(|n| if n == 0 { 1.0 } else { 0.7f32.powi(n as i32) });
        // Target = filtered excitation at delay `period`.
        let rows = filtered_excitation_rows(&exc, &h, period, period);
        let x = rows[0];

        let (t_min, t_max) = subframe1_window(period as i32);
        let res = closed_loop_search(&x, &h, &exc, t_min, t_max, T1_FRACTIONAL_LIMIT);
        assert_eq!(res.delay.int_t, period as i32, "got {:?}", res.delay);
        assert_eq!(res.delay.frac, 0);
    }

    /// Delays ≥ 85 stay integer-only for T1 even when fractions would
    /// score higher; for the same setup a subframe-2 style search (no
    /// limit) may return fractions.
    #[test]
    fn t1_integer_only_above_limit() {
        let period = 100usize;
        let exc: [f32; EXC_BUFFER] = std::array::from_fn(|n| {
            let ph = (n % period) as f32 / period as f32;
            (std::f32::consts::TAU * ph).sin() + 0.2 * (2.0 * std::f32::consts::TAU * ph).cos()
        });
        let h: [f32; L_SUBFR] =
            std::array::from_fn(|n| if n == 0 { 1.0 } else { 0.6f32.powi(n as i32) });
        let rows = filtered_excitation_rows(&exc, &h, period, period);
        let x = rows[0];
        let (t_min, t_max) = subframe1_window(period as i32);
        let res = closed_loop_search(&x, &h, &exc, t_min, t_max, T1_FRACTIONAL_LIMIT);
        assert_eq!(
            res.delay.frac, 0,
            "integer-only region must not emit fractions"
        );
    }

    /// The emitted delay is always encodable: `encode_p1` accepts every
    /// subframe-1 result across a sweep of open-loop anchors.
    #[test]
    fn results_are_p1_encodable() {
        use crate::pitch_decode::encode_p1;
        let exc: [f32; EXC_BUFFER] = std::array::from_fn(|n| {
            let ph = (n % 45) as f32 / 45.0;
            (std::f32::consts::TAU * ph).sin()
        });
        let h: [f32; L_SUBFR] =
            std::array::from_fn(|n| if n == 0 { 1.0 } else { 0.5f32.powi(n as i32) });
        for t_op in [20i32, 45, 84, 90, 120, 143] {
            let (t_min, t_max) = subframe1_window(t_op);
            let rows = filtered_excitation_rows(&exc, &h, t_min as usize, t_min as usize);
            let x = rows[0];
            let res = closed_loop_search(&x, &h, &exc, t_min, t_max, T1_FRACTIONAL_LIMIT);
            assert!(
                encode_p1(res.delay).is_some(),
                "unencodable T1 {:?} for t_op {t_op}",
                res.delay
            );
        }
    }
}
