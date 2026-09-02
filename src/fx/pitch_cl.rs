//! §3.7 **closed-loop pitch search** and the §3.7.3 adaptive-codebook
//! gain on the fixed grid.
//!
//! ## Number grids
//!
//! - Target `x(n)` — Q0; impulse response `h(n)` — Q12; excitation
//!   `u(n)` — Q0 with the current subframe holding the LP residual
//!   (clause 3.7).
//! - `y_k(n)` — the eq (38) recursion runs on the Word32 Q13
//!   accumulator grid (`l_mult`/`l_mac` of Q0 × Q12), which keeps the
//!   shift-and-add exact across the delay range; each row lands on Q0
//!   (`<< 3`, round) for the eq (37) sums.
//! - eq (37) — numerator `Σ x·y_k` and energy `Σ y_k²` on an exact
//!   wide accumulator, then the normalised correlation as a Word16
//!   mantissa product (`norm_l` of the numerator, the Q30 `inv_sqrt`
//!   mantissa of the energy) with tracked exponents, re-aligned to a
//!   common Word32 scale before the comparisons and the eq (39)
//!   interpolation (Q15 `b12` taps through `mpy_32_16`).
//! - eq (43) — `g_p = xᵗy / yᵗy` through `div_s` on Word16 mantissas,
//!   landed on Q14 and bounded to `[0, 1.2]` (`19661`).
//!
//! The prose pins the equations, the windows and the fraction policy;
//! the accumulator scaling, the normalisation precision and the exact
//! candidate set of the fractional pass are the fixed chain's own
//! composition, exposed as latitude and pinned black-box.

use crate::fx::dsp::{inv_sqrt, l_extract, mpy_32_16};
use crate::fx::excitation::EXC_BUF;
use crate::fx::filters::L_SUBFR;
use crate::fx::ops::{div_s, extract_h, l_add, l_mac, l_mult, l_shl, l_shr, norm_l, round, shr};
use crate::pitch_decode::PitchDelay;
use crate::tables::PITCH_INTERP_FILTER_ANALYSIS_Q15;

/// Minimum pitch delay.
pub const PIT_MIN: i32 = 20;
/// Maximum pitch delay.
pub const PIT_MAX: i32 = 143;
/// Integer-only threshold for `T1` (clause 3.7).
pub const T1_FRACTIONAL_LIMIT: i32 = 85;
/// `1.2` on Q14 — the eq (43) upper bound.
pub const GP_MAX_Q14: i16 = 19661;

/// Which fractional candidates the eq (39) pass evaluates around the
/// best integer delay `k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FracMode {
    /// `k − 2/3, k − 1/3, k + 1/3, k + 2/3` interpolated against the
    /// raw `R(k)` (the float module's reading).
    RawCentre,
    /// `k − 2/3 … k + 2/3` in 1/3 steps, all five interpolated
    /// (including the fraction-0 value).
    #[default]
    Five,
    /// `k − 1/3, k, k + 1/3`, all interpolated.
    Three,
}

/// Unstated latitude of the closed-loop search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClosedLoopLatitude {
    /// Fractional candidate set.
    pub frac: FracMode,
    /// Exact double-precision `R(k)` (the spec-equation oracle) instead
    /// of the Word16 mantissa form.
    pub wide: bool,
    /// eq (43) on an exact quotient instead of `div_s` mantissas.
    pub gain_exact: bool,
    /// A later integer delay must beat the running maximum strictly.
    pub strict_max: bool,
    /// eq (39) phase geometry: `R(k)_t` is the correlation at delay
    /// `k − t/3` (`true`, the same negated fold the decoder's eq (40)
    /// interpolator was corpus-pinned to) or `k + t/3` (`false`, the
    /// literal reading).
    pub negated_phase: bool,
}

impl Default for ClosedLoopLatitude {
    fn default() -> Self {
        Self {
            frac: FracMode::Five,
            wide: false,
            gain_exact: false,
            strict_max: true,
            negated_phase: false,
        }
    }
}

/// One subframe's search result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClosedLoopResultFx {
    /// The selected fractional delay.
    pub delay: PitchDelay,
}

/// The §3.7 subframe-1 window around `T_op`.
#[must_use]
pub fn subframe1_window(t_op: i32) -> (i32, i32) {
    let mut t_min = t_op - 3;
    if t_min < PIT_MIN {
        t_min = PIT_MIN;
    }
    let mut t_max = t_min + 6;
    if t_max > PIT_MAX {
        t_max = PIT_MAX;
        t_min = t_max - 6;
    }
    (t_min, t_max)
}

/// The §3.7 subframe-2 window around `int(T1)`.
#[must_use]
pub fn subframe2_window(int_t1: i32) -> (i32, i32) {
    let mut t_min = int_t1 - 5;
    if t_min < PIT_MIN {
        t_min = PIT_MIN;
    }
    let mut t_max = t_min + 9;
    if t_max > PIT_MAX {
        t_max = PIT_MAX;
        t_min = t_max - 9;
    }
    (t_min, t_max)
}

/// eq (37) numerator and energy at one delay row (Q0 `y`).
fn sums(x: &[i16; L_SUBFR], y: &[i16; L_SUBFR]) -> (i64, i64) {
    let mut num = 0i64;
    let mut den = 0i64;
    for n in 0..L_SUBFR {
        num += i64::from(x[n]) * i64::from(y[n]);
        den += i64::from(y[n]) * i64::from(y[n]);
    }
    (num, den)
}

/// A signed wide value as a Word16 mantissa with a base-2 exponent:
/// `v ≈ mant · 2^exp`.
fn mantissa16(v: i64) -> (i16, i32) {
    if v == 0 {
        return (0, 0);
    }
    let mut shift = 0i32;
    let mut w = v;
    while w > i64::from(i32::MAX) || w < i64::from(i32::MIN) {
        w >>= 1;
        shift += 1;
    }
    let w32 = w as i32;
    let n = norm_l(w32);
    (extract_h(l_shl(w32, n)), shift + 16 - i32::from(n))
}

/// The eq (37) normalised correlation `num / √den` as `(mant, exp)`
/// on the Word16 mantissa grid (`mant` a Word32 product of two Word16
/// mantissas).
fn normalised(num: i64, den: i64) -> (i32, i32) {
    if den <= 0 {
        return (0, 0);
    }
    // Energy as a Word32 with an even shift so the square root of the
    // scale is exact.
    let mut shift = 0i32;
    let mut d = den;
    while d > i64::from(i32::MAX) {
        d >>= 2;
        shift += 2;
    }
    let inv = inv_sqrt(d as i32); // (1/√d) · 2^30
    let ni = norm_l(inv);
    let inv16 = extract_h(l_shl(inv, ni));
    // 1/√den = inv · 2^(−30) · 2^(−shift/2); inv = inv16 · 2^(16 − ni).
    let e_inv = -30 - shift / 2 + 16 - i32::from(ni);
    let (m_num, e_num) = mantissa16(num);
    (l_mult(m_num, inv16), e_num + e_inv)
}

/// `b12(i)` from the compiled one-sided table (`|i| ≤ 12`).
fn b12(i: i32) -> i16 {
    let idx = i.unsigned_abs() as usize;
    if idx >= PITCH_INTERP_FILTER_ANALYSIS_Q15.len() {
        0
    } else {
        PITCH_INTERP_FILTER_ANALYSIS_Q15[idx]
    }
}

/// The eq (39) evaluation point `(k, t)` for the delay `k + d/3`
/// (`d ∈ −2 … 2`) under either phase geometry.
fn eval_point(k: i32, d: i32, negated: bool) -> (i32, i32) {
    if negated {
        // R(k)_t ↔ delay k − t/3.
        match d {
            -2 => (k, 2),
            -1 => (k, 1),
            0 => (k, 0),
            1 => (k + 1, 2),
            _ => (k + 1, 1),
        }
    } else {
        // R(k)_t ↔ delay k + t/3.
        match d {
            -2 => (k - 1, 1),
            -1 => (k - 1, 2),
            0 => (k, 0),
            1 => (k, 1),
            _ => (k, 2),
        }
    }
}

/// eq (39) on the aligned Word32 correlation sequence.
fn interpolate(r: &[i32], k_base: i32, k: i32, t: i32) -> i32 {
    let at = |kk: i32| -> i32 {
        let idx = kk - k_base;
        if idx < 0 || idx as usize >= r.len() {
            0
        } else {
            r[idx as usize]
        }
    };
    let mut acc = 0i32;
    for i in 0..4 {
        let (a_hi, a_lo) = l_extract(at(k - i));
        acc = l_add(acc, mpy_32_16(a_hi, a_lo, b12(t + 3 * i)));
        let (b_hi, b_lo) = l_extract(at(k + 1 + i));
        acc = l_add(acc, mpy_32_16(b_hi, b_lo, b12(3 - t + 3 * i)));
    }
    acc
}

/// The eq (37) normalised correlations for every integer delay in
/// `k_lo ..= k_hi`, on a common Word32 scale (row 0 ↔ `k_lo`).
fn correlation_rows(
    x: &[i16; L_SUBFR],
    h_q12: &[i16; L_SUBFR],
    exc: &[i16; EXC_BUF],
    off: usize,
    k_lo: i32,
    k_hi: i32,
    lat: &ClosedLoopLatitude,
) -> Vec<i32> {
    let n_rows = (k_hi - k_lo + 1) as usize;

    // eq (38) rows on the Q13 accumulator grid, landed on Q0.
    let mut acc = [0i32; L_SUBFR];
    let u = |n: i32| -> i16 { exc[(off as i32 + n) as usize] };
    for n in 0..L_SUBFR as i32 {
        let mut a = 0i32;
        for j in 0..=n {
            a = l_mac(a, u(j - k_lo), h_q12[(n - j) as usize]);
        }
        acc[n as usize] = a;
    }
    let mut corr: Vec<(i64, i64)> = Vec::with_capacity(n_rows);
    let mut y = [0i16; L_SUBFR];
    for k in k_lo..=k_hi {
        if k > k_lo {
            let u_mk = u(-k);
            for n in (1..L_SUBFR).rev() {
                acc[n] = l_mac(acc[n - 1], u_mk, h_q12[n]);
            }
            acc[0] = l_mult(u_mk, h_q12[0]);
        }
        for n in 0..L_SUBFR {
            y[n] = round(l_shl(acc[n], 3));
        }
        corr.push(sums(x, &y));
    }

    // Normalised correlations on a common Word32 scale.
    if lat.wide {
        // Exact values scaled into the Word32 range (common scale).
        let vals: Vec<f64> = corr
            .iter()
            .map(|&(num, den)| {
                if den > 0 {
                    num as f64 / (den as f64).sqrt()
                } else {
                    0.0
                }
            })
            .collect();
        let peak = vals.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-9);
        let scale = (2f64.powi(29)) / peak;
        vals.iter().map(|v| (v * scale).round() as i32).collect()
    } else {
        let pairs: Vec<(i32, i32)> = corr.iter().map(|&(n, d)| normalised(n, d)).collect();
        let e_max = pairs.iter().map(|p| p.1).max().unwrap_or(0);
        pairs
            .iter()
            .map(|&(m, e)| {
                let s = e_max - e;
                if s >= 31 {
                    0
                } else {
                    l_shr(m, s as i16)
                }
            })
            .collect()
    }
}

/// Diagnostic: the eq (37)/(39) scores of two fractional delays on one
/// common scale (the fraction-0 score is the interpolated value, as in
/// [`FracMode::Five`]).
#[doc(hidden)]
#[must_use]
pub fn compare_delays_fx(
    x: &[i16; L_SUBFR],
    h_q12: &[i16; L_SUBFR],
    exc: &[i16; EXC_BUF],
    off: usize,
    a: PitchDelay,
    b: PitchDelay,
    lat: &ClosedLoopLatitude,
) -> (i32, i32) {
    let k_lo = (a.int_t.min(b.int_t) - 6).max(1);
    let k_hi = (a.int_t.max(b.int_t) + 6).min(PIT_MAX);
    let r = correlation_rows(x, h_q12, exc, off, k_lo, k_hi, lat);
    let score = |d: PitchDelay| -> i32 {
        let (ek, et) = eval_point(d.int_t, d.frac, lat.negated_phase);
        interpolate(&r, k_lo, ek, et)
    };
    (score(a), score(b))
}

/// §3.7 over one subframe. `exc` is the decoder-layout excitation
/// buffer with the current subframe at `off` already holding the LP
/// residual (clause 3.7); `h` the Q12 impulse response.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn closed_loop_search_fx(
    x: &[i16; L_SUBFR],
    h_q12: &[i16; L_SUBFR],
    exc: &[i16; EXC_BUF],
    off: usize,
    t_min: i32,
    t_max: i32,
    frac_limit: i32,
    lat: &ClosedLoopLatitude,
) -> ClosedLoopResultFx {
    let k_lo = (t_min - 4).max(PIT_MIN - 4).max(1);
    let k_hi = (t_max + 4).min(PIT_MAX);
    let r = correlation_rows(x, h_q12, exc, off, k_lo, k_hi, lat);
    let r_at = |k: i32| -> i32 {
        let idx = k - k_lo;
        if idx < 0 || idx as usize >= r.len() {
            i32::MIN
        } else {
            r[idx as usize]
        }
    };

    // Integer pass over the window.
    let mut best_k = t_min;
    let mut best_r = i32::MIN;
    for k in t_min..=t_max {
        let v = r_at(k);
        let better = if lat.strict_max {
            v > best_r
        } else {
            v >= best_r
        };
        if better {
            best_r = v;
            best_k = k;
        }
    }

    let mut best = PitchDelay {
        int_t: best_k,
        frac: 0,
    };
    if best_k < frac_limit {
        // Candidate delays as thirds around k: −2 … +2.
        let thirds: &[i32] = match lat.frac {
            FracMode::RawCentre => &[-2, -1, 1, 2],
            FracMode::Five => &[-2, -1, 0, 1, 2],
            FracMode::Three => &[-1, 0, 1],
        };
        let mut best_score = match lat.frac {
            FracMode::RawCentre => best_r,
            _ => i32::MIN,
        };
        for &d in thirds {
            // Transmitted form of k + d/3.
            let (oi, of) = match d {
                -2 => (best_k - 1, 1),
                -1 => (best_k, -1),
                0 => (best_k, 0),
                1 => (best_k, 1),
                _ => (best_k + 1, -1),
            };
            // The transmitted-domain reach: T1 ≥ 19⅓, both ≤ 143.
            if 3 * oi + of < 3 * PIT_MIN - 2 || oi > PIT_MAX {
                continue;
            }
            let (ek, et) = eval_point(best_k, d, lat.negated_phase);
            let v = interpolate(&r, k_lo, ek, et);
            if v > best_score {
                best_score = v;
                best = PitchDelay {
                    int_t: oi,
                    frac: of,
                };
            }
        }
    }
    ClosedLoopResultFx { delay: best }
}

/// eq (43): the adaptive-codebook gain `g_p = xᵗy / yᵗy` bounded to
/// `[0, 1.2]`, on Q14, under the default latitude.
#[must_use]
pub fn pitch_gain_q14(x: &[i16; L_SUBFR], y: &[i16; L_SUBFR]) -> i16 {
    pitch_gain_q14_lat(x, y, &ClosedLoopLatitude::default())
}

/// eq (43) under an explicit latitude.
#[must_use]
pub fn pitch_gain_q14_lat(x: &[i16; L_SUBFR], y: &[i16; L_SUBFR], lat: &ClosedLoopLatitude) -> i16 {
    let (xy, yy) = sums(x, y);
    if xy <= 0 || yy <= 0 {
        return 0;
    }
    if lat.gain_exact {
        let g = xy as f64 / yy as f64;
        return (g * 16384.0).round().clamp(0.0, f64::from(GP_MAX_Q14)) as i16;
    }
    // Word16 mantissas: xy = mx · 2^ex, yy = my · 2^ey; div_s needs
    // 0 ≤ mx ≤ my, so halve mx (exponent + 1) when it is larger.
    let (mut mx, mut ex) = mantissa16(xy);
    let (my, ey) = mantissa16(yy);
    if mx > my {
        mx = shr(mx, 1);
        ex += 1;
    }
    let q = div_s(mx, my); // (mx/my) on Q15
                           // g_p = q · 2^(ex − ey), wanted on Q14: shift by (ex − ey − 1).
    let shift = ex - ey - 1;
    let g = if shift >= 0 {
        l_shl(i32::from(q), shift as i16)
    } else {
        l_shr(i32::from(q), (-shift) as i16)
    };
    g.clamp(0, i32::from(GP_MAX_Q14)) as i16
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fx::excitation::PIT_MAX as PIT_MAX_U;

    const OFF: usize = PIT_MAX_U + 10;

    /// Windows: interior, low clamp, high clamp.
    #[test]
    fn search_windows() {
        assert_eq!(subframe1_window(60), (57, 63));
        assert_eq!(subframe1_window(20), (20, 26));
        assert_eq!(subframe1_window(143), (137, 143));
        assert_eq!(subframe2_window(60), (55, 64));
        assert_eq!(subframe2_window(21), (20, 29));
        assert_eq!(subframe2_window(143), (134, 143));
    }

    fn decaying_h() -> [i16; L_SUBFR] {
        std::array::from_fn(|n| (4096.0 * 0.7f64.powi(n as i32)) as i16)
    }

    /// A periodic excitation whose target is its own filtered image at
    /// a known integer delay is recovered at that delay with fraction 0.
    #[test]
    fn recovers_known_integer_delay() {
        let period = 58i32;
        let mut exc = [0i16; EXC_BUF];
        for (n, slot) in exc.iter_mut().enumerate() {
            let ph = (n % period as usize) as f64 / period as f64;
            *slot = (3000.0 * (std::f64::consts::TAU * ph).sin()) as i16;
        }
        let h = decaying_h();
        // Target = filtered excitation at delay `period`.
        let v: [i16; L_SUBFR] = std::array::from_fn(|n| exc[OFF + n - period as usize]);
        let x = crate::fx::filters::convolve_h_q0(&v, &h);
        let (t_min, t_max) = subframe1_window(period);
        let r = closed_loop_search_fx(
            &x,
            &h,
            &exc,
            OFF,
            t_min,
            t_max,
            T1_FRACTIONAL_LIMIT,
            &ClosedLoopLatitude::default(),
        );
        assert_eq!(r.delay.int_t, period, "{:?}", r.delay);
        assert_eq!(r.delay.frac, 0);
    }

    /// Delays ≥ 85 stay integer-only for `T1`.
    #[test]
    fn t1_integer_only_above_limit() {
        let period = 100i32;
        let mut exc = [0i16; EXC_BUF];
        for (n, slot) in exc.iter_mut().enumerate() {
            let ph = (n % period as usize) as f64 / period as f64;
            *slot = (3000.0 * (std::f64::consts::TAU * ph).sin()
                + 600.0 * (2.0 * std::f64::consts::TAU * ph).cos()) as i16;
        }
        let h = decaying_h();
        let v: [i16; L_SUBFR] = std::array::from_fn(|n| exc[OFF + n - period as usize + 1]);
        let x = crate::fx::filters::convolve_h_q0(&v, &h);
        let (t_min, t_max) = subframe1_window(period);
        let r = closed_loop_search_fx(
            &x,
            &h,
            &exc,
            OFF,
            t_min,
            t_max,
            T1_FRACTIONAL_LIMIT,
            &ClosedLoopLatitude::default(),
        );
        assert_eq!(r.delay.frac, 0);
    }

    /// eq (43): `x = y` gives 1.0 (Q14 16384); a doubled target clamps
    /// at 1.2; an anti-correlated target gives 0; the `div_s` form
    /// tracks the exact quotient within a Q14 LSB or two.
    #[test]
    fn pitch_gain_basics() {
        let y: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 7 % 13) as i16) * 300 - 1800);
        let g1 = pitch_gain_q14(&y, &y);
        assert!((i32::from(g1) - 16384).abs() <= 2, "{g1}");
        let x2: [i16; L_SUBFR] = std::array::from_fn(|n| y[n] * 2);
        assert_eq!(pitch_gain_q14(&x2, &y), GP_MAX_Q14);
        let xn: [i16; L_SUBFR] = std::array::from_fn(|n| -y[n]);
        assert_eq!(pitch_gain_q14(&xn, &y), 0);
        let x3: [i16; L_SUBFR] = std::array::from_fn(|n| (i32::from(y[n]) * 3 / 7) as i16);
        let exact = pitch_gain_q14_lat(
            &x3,
            &y,
            &ClosedLoopLatitude {
                gain_exact: true,
                ..Default::default()
            },
        );
        let fx = pitch_gain_q14(&x3, &y);
        assert!(
            (i32::from(exact) - i32::from(fx)).abs() <= 2,
            "{exact} vs {fx}"
        );
    }

    /// The Word16 mantissa normalisation orders candidates like the
    /// exact quotient on a spread of (num, den) pairs.
    #[test]
    fn normalised_orders_like_exact() {
        let pairs: [(i64, i64); 5] = [
            (1_000_000, 4_000_000_000),
            (-500, 30_000),
            (123_456_789_012, 987_654_321_000),
            (7, 9),
            (250_000, 62_500_000),
        ];
        let fx: Vec<(i32, i32)> = pairs.iter().map(|&(n, d)| normalised(n, d)).collect();
        let e_max = fx.iter().map(|p| p.1).max().unwrap();
        let aligned: Vec<f64> = fx
            .iter()
            .map(|&(m, e)| f64::from(m) * 2f64.powi(e - e_max))
            .collect();
        let exact: Vec<f64> = pairs
            .iter()
            .map(|&(n, d)| n as f64 / (d as f64).sqrt())
            .collect();
        let scale = aligned[0] / exact[0];
        for i in 0..5 {
            assert!(
                (aligned[i] - exact[i] * scale).abs() <= 1e-3 * exact[i].abs() * scale + 1.0,
                "{i}: {} vs {}",
                aligned[i],
                exact[i] * scale
            );
        }
    }
}
