//! §3.2.3 **fixed-point LP→LSP conversion** — the Chebyshev-series
//! grid search of eqs (11)–(17) on the clause-5 Word16/Word32
//! operator grid.
//!
//! The float module [`crate::lp_to_lsp`] evaluates the sum/difference
//! half-polynomials `G1`/`G2` in `f64`; here the whole conversion
//! runs on the fixed grid: Q11 half-polynomial coefficients built
//! from the Q12 LP set, a Clenshaw recurrence carried in Word32 on
//! the Q23 grid (`mpy_32_16` against the Q15 abscissa), the
//! 60-interval staged cosine grid with four interval bisections, and
//! a final `div_s` linear interpolation inside the refined bracket.
//!
//! ## Grids
//!
//! - LP input `a_i` on Q12.
//! - Half-polynomial coefficients `f1/f2` on **Q11** (the divide-out
//!   recursion doubles values up to ≈ ±16, one integer bit beyond the
//!   Q12 ceiling).
//! - Chebyshev evaluation on **Q23** Word32 (Q11 coefficients ×
//!   `2^12`), abscissae on the staged Q15 cosine grid
//!   [`crate::tables::LSF_SEARCH_GRID_COS_Q15`].
//! - Output LSPs on Q15 (cosine domain).

use crate::fx::analysis::Dpf;
use crate::fx::dsp::{l_extract, mpy_32_16};
use crate::fx::ops::{add, div_s, l_abs, l_shl, l_sub, mult, norm_l, shr, shr_r, sub};
use crate::tables::{LSF_SEARCH_GRID_COS_Q15, M};

/// Half the LP order — the root count of each of `F1` / `F2`.
const NC: usize = M / 2;
/// Bisection passes inside a sign-change bracket (clause 3.2.3).
const N_BISECT: usize = 4;

/// Builds the eq (15)-family `f1[0…5]` / `f2[0…5]` half-polynomial
/// coefficients on the Q11 grid from Q12 LP coefficients
/// (`a_in[i−1] = a_i`, `a_0 = 1` implicit):
///
/// `f1[i] = (a_i + a_{M+1−i}) − f1[i−1]`,
/// `f2[i] = (a_i − a_{M+1−i}) + f2[i−1]` — the Q12→Q11 landing uses
/// the rounding shift (`shr_r`).
fn build_polys_fx(a: &[i16; M]) -> ([i16; NC + 1], [i16; NC + 1]) {
    let mut f1 = [0i16; NC + 1];
    let mut f2 = [0i16; NC + 1];
    f1[0] = 2048; // 1.0 on Q11
    f2[0] = 2048;
    for i in 1..=NC {
        // a_i ± a_{M+1−i} on Q12, halved onto Q11 with rounding.
        let sum = shr_r(add(a[i - 1], a[M - i]), 1);
        let diff = shr_r(sub(a[i - 1], a[M - i]), 1);
        f1[i] = sub(sum, f1[i - 1]);
        f2[i] = add(diff, f2[i - 1]);
    }
    (f1, f2)
}

/// Clenshaw evaluation of the Chebyshev series `G(x)` for a Q11
/// coefficient set at a Q15 abscissa, returning the value on the Q23
/// Word32 grid:
///
/// `b_k = 2x·b_{k+1} − b_{k+2} + f[5−k]`, `G = x·b_1 − b_2 + ½·f[5]`.
fn cheb_fx(x_q15: i16, f: &[i16; NC + 1]) -> i32 {
    let mut b1: i32 = 0;
    let mut b2: i32 = 0;
    for k in (1..=NC).rev() {
        let (h, l) = l_extract(b1);
        // 2·x·b1 on Q23: mpy_32_16 lands on Q23, one left shift doubles.
        let t = l_shl(mpy_32_16(h, l, x_q15), 1);
        // + f[5−k] on Q23 (l_mac's doubling folded into the 2^11 factor).
        let b0 = crate::fx::ops::l_mac(l_sub(t, b2), f[NC - k], 2048);
        b2 = b1;
        b1 = b0;
    }
    let (h, l) = l_extract(b1);
    let t = mpy_32_16(h, l, x_q15);
    // ½·f[5] on Q23.
    crate::fx::ops::l_mac(l_sub(t, b2), f[NC], 1024)
}

/// Linear interpolation for the root inside a refined bracket
/// `[x_hi, x_lo]` (grid order: `x_lo > x_hi`) with the polynomial
/// values `y_lo`/`y_hi` of opposite sign: `x = x_lo + (x_hi −
/// x_lo)·|y_lo|/|y_hi − y_lo|`, the ratio through `div_s` on aligned
/// high words.
fn interpolate_root(x_lo: i16, x_hi: i16, y_lo: i32, y_hi: i32) -> i16 {
    let dy = l_sub(y_hi, y_lo);
    if dy == 0 {
        return x_lo;
    }
    let n = norm_l(dy);
    let den = crate::fx::ops::extract_h(l_shl(l_abs(dy), n));
    let num = crate::fx::ops::extract_h(l_shl(l_abs(y_lo), n));
    let ratio = if num >= den { 32767 } else { div_s(num, den) };
    add(x_lo, mult(sub(x_hi, x_lo), ratio))
}

/// §3.2.3 fixed-point LP→LSP conversion: Q12 LP coefficients to ten
/// Q15 cosine-domain LSPs (decreasing in value / increasing in `ω`),
/// walking the staged 60-interval cosine grid with four bisections
/// and a `div_s` interpolation per root, alternating `F1`/`F2`.
///
/// Returns `None` when fewer than ten sign changes are found (the
/// clause-3.2.3 fallback: keep the previous frame's LSPs). Accepts
/// the Q27 DPF LP set from [`crate::fx::levinson`] via
/// [`lp_to_lsp_fx_q27`] or a plain Q12 set here.
#[must_use]
pub fn lp_to_lsp_fx(a_q12: &[i16; M]) -> Option<[i16; M]> {
    let (f1, f2) = build_polys_fx(a_q12);

    let mut q = [0i16; M];
    let mut found = 0usize;
    let mut on_f1 = true;
    let eval = |x: i16, on_f1: bool| -> i32 {
        if on_f1 {
            cheb_fx(x, &f1)
        } else {
            cheb_fx(x, &f2)
        }
    };

    let mut x_lo = LSF_SEARCH_GRID_COS_Q15[0];
    let mut y_lo = eval(x_lo, on_f1);
    for &x_grid in &LSF_SEARCH_GRID_COS_Q15[1..] {
        let x_hi = x_grid;
        let mut y_hi = eval(x_hi, on_f1);
        if (i64::from(y_lo) * i64::from(y_hi)) <= 0 {
            // Refine by four interval bisections (mean via halves —
            // the Q15 endpoints straddle ±2^14 so a direct add could
            // saturate).
            let (mut xl, mut yl) = (x_lo, y_lo);
            let (mut xh, mut yh) = (x_hi, y_hi);
            for _ in 0..N_BISECT {
                let xm = add(shr(xl, 1), shr(xh, 1));
                let ym = eval(xm, on_f1);
                if (i64::from(yl) * i64::from(ym)) <= 0 {
                    xh = xm;
                    yh = ym;
                } else {
                    xl = xm;
                    yl = ym;
                }
            }
            q[found] = interpolate_root(xl, xh, yl, yh);
            found += 1;
            if found == M {
                break;
            }
            on_f1 = !on_f1;
            y_hi = eval(x_hi, on_f1);
        }
        x_lo = x_hi;
        y_lo = y_hi;
    }

    if found == M {
        Some(q)
    } else {
        None
    }
}

/// Convenience for the Q27 DPF LP set of
/// [`crate::fx::levinson::LevinsonFxResult`]: rounds each tap onto
/// Q12 (the §3.2.3 search precision) and converts.
#[must_use]
pub fn lp_to_lsp_fx_q27(a_q27: &[Dpf; M]) -> Option<[i16; M]> {
    let a_q12: [i16; M] = std::array::from_fn(|j| {
        crate::fx::ops::round(l_shl(crate::fx::dsp::l_comp(a_q27[j].0, a_q27[j].1), 1))
    });
    lp_to_lsp_fx(&a_q12)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fx::lsp::STARTUP_LSP_Q15;
    use crate::lp_to_lsp::lp_to_lsp;

    /// Valid LSP set: strictly decreasing, inside (−1, 1).
    fn assert_valid(q: &[i16; M]) {
        for w in q.windows(2) {
            assert!(w[0] > w[1], "LSPs must decrease: {q:?}");
        }
    }

    /// The flat predictor's LSPs are the uniform cosine grid — the
    /// same grid Table 9 initialises (within the grid-search
    /// resolution of the 60-interval walk).
    #[test]
    fn flat_predictor_hits_uniform_grid() {
        let q = lp_to_lsp_fx(&[0i16; M]).expect("flat spectrum");
        assert_valid(&q);
        for (i, (&got, &want)) in q.iter().zip(STARTUP_LSP_Q15.iter()).enumerate() {
            assert!(
                (i32::from(got) - i32::from(want)).abs() <= 96,
                "flat LSP {i}: {got} vs Table 9 {want}"
            );
        }
    }

    /// Against the float search on identical Q12 inputs the Q15 roots
    /// agree within the bisection-step tolerance across a spread of
    /// stable spectra.
    #[test]
    fn tracks_float_search() {
        use crate::fx::levinson::levinson_fx;
        use crate::lp_analysis::N_LAGS;

        for (decay, freq) in [(0.9f64, 0.3f64), (0.85, 0.2), (0.92, 0.35), (0.8, 0.0)] {
            let rf: [f64; N_LAGS] =
                std::array::from_fn(|k| decay.powi(k as i32) * (freq * k as f64).cos());
            // Onto the normalised Word32 grid.
            let mut scale = 0i32;
            while rf[0] * (2f64).powi(scale) < (2f64).powi(30) {
                scale += 1;
            }
            let dpf: [Dpf; N_LAGS] =
                std::array::from_fn(|k| l_extract((rf[k] * (2f64).powi(scale)).floor() as i32));
            let lev = levinson_fx(&dpf).expect("stable");

            let q_fx = lp_to_lsp_fx(&lev.a_q12).expect("roots");
            assert_valid(&q_fx);
            let a_f32: [f32; M] = std::array::from_fn(|j| f32::from(lev.a_q12[j]) / 4096.0);
            let q_float = lp_to_lsp(&a_f32).expect("float roots");
            for j in 0..M {
                let want = (f64::from(q_float[j]) * 32_768.0).round() as i32;
                let got = i32::from(q_fx[j]);
                assert!(
                    (got - want).abs() <= 24,
                    "({decay},{freq}) root {j}: fx {got} vs float {want}"
                );
            }
            // The Q27 entry point agrees with the Q12 one.
            assert_eq!(lp_to_lsp_fx_q27(&lev.a_q27).expect("q27"), q_fx);
        }
    }

    /// The interpolation helper lands inside the bracket and on the
    /// exact crossing for a symmetric bracket.
    #[test]
    fn interpolation_midpoint() {
        let x = interpolate_root(1000, 900, -(1 << 20), 1 << 20);
        assert!((949..=951).contains(&x), "symmetric crossing: {x}");
        let x = interpolate_root(1000, 900, 0, 1 << 20);
        assert_eq!(x, 1000, "root at the low edge");
    }

    /// An LP set whose polynomial never changes sign ten times inside
    /// the grid reports failure (caller keeps previous LSPs).
    #[test]
    fn degenerate_reports_none() {
        // A wildly unstable coefficient set (far outside minimum
        // phase) — the grid walk cannot find ten interlaced roots.
        let a: [i16; M] = [
            32767, -32768, 32767, -32768, 32767, -32768, 32767, -32768, 32767, -32768,
        ];
        assert!(lp_to_lsp_fx(&a).is_none());
    }
}
