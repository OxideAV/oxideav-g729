//! §3.2.3 **LP-to-LSP conversion** (forward) — the encoder-side inverse
//! of the decoder's [`crate::lsp_to_lp`]. Turns the 10 LP coefficients
//! `a_i` produced by [`crate::levinson`] into the cosine-domain line
//! spectral pairs `q_i = cos(ω_i)` (`0 < ω_1 < … < ω_10 < π`) that the
//! §3.2.4 LSP quantiser encodes.
//!
//! ## Spec source — clause 3.2.3, equations (11)–(17)
//!
//! Per clause 3.2.3 the LSPs are the roots on the unit circle of the sum
//! and difference polynomials
//!
//! ```text
//! F1'(z) = A(z) + z^-(M+1)·A(z^-1)        (11)  — symmetric, root at z = −1
//! F2'(z) = A(z) − z^-(M+1)·A(z^-1)        (12)  — antisymmetric, root at z = +1
//! ```
//!
//! The known roots at `z = ±1` are divided out (`F1 = F1'/(1+z^-1)`,
//! `F2 = F2'/(1−z^-1)`), leaving two symmetric degree-10 polynomials
//! whose five conjugate root pairs each give five LSPs. Because both are
//! symmetric in `z`, evaluating them on the unit circle `z = e^{jω}`
//! collapses to a real Chebyshev series in `x = cos ω` (eqs (16)/(17)):
//! writing the polynomial coefficients `f(0 … 5)`,
//!
//! ```text
//! F1(e^{jω}) / z^-5 = 2·G1(x),
//!   G1(x) = f1[0]·T5(x) + f1[1]·T4(x) + … + f1[4]·T1(x) + ½·f1[5]·T0(x)
//! ```
//!
//! with `T_k(x) = cos(kω)` the Chebyshev polynomials. The coefficients
//! `f1(i)` / `f2(i)` are obtained from `a_i` by the divide-out recursion
//!
//! ```text
//! f1[0] = 1,  f1[i] = (a_i + a_{M+1−i}) − f1[i−1]     i = 1 … 5
//! f2[0] = 1,  f2[i] = (a_i − a_{M+1−i}) + f2[i−1]     i = 1 … 5
//! ```
//!
//! (the `− f1[i−1]` / `+ f2[i−1]` terms are exactly the polynomial
//! division by `(1 + z^-1)` / `(1 − z^-1)`). `G1`/`G2` are evaluated with
//! the Clenshaw recurrence for a Chebyshev sum.
//!
//! ## Root search (clause 3.2.3)
//!
//! The ten LSPs interlace between `F1` (giving `q_1, q_3, q_5, q_7, q_9`)
//! and `F2` (giving `q_2, q_4, q_6, q_8, q_10`), the smallest-`ω`
//! (largest-`x`) root belonging to `F1`. The search walks the
//! **60-interval cosine grid** [`crate::tables::LSF_SEARCH_GRID_COS_Q15`]
//! (61 samples from `x ≈ +1` down to `x ≈ −1`), evaluating the *current*
//! polynomial and detecting a sign change; on each sign change the root
//! is refined by **four interval bisections** followed by one linear
//! interpolation (clause 3.2.3), the found `q` is recorded, and the
//! active polynomial toggles `F1 ↔ F2`. This yields `q_1 … q_10` in
//! decreasing-`x` (increasing-`ω`) order — exactly the layout
//! [`crate::lsp_to_lp`] consumes (`q_in[i−1] = q_i`).
//!
//! ## Numeric domain
//!
//! `f64` internally (the grid search needs the precision), input `a_i`
//! and output `q_i` are `f32` to match the decoder LP / LSP surfaces.

use crate::tables::{LSF_SEARCH_GRID_COS_Q15, M};

/// Q15 scale for the cosine search grid (`2^15`).
const Q15: f64 = 32_768.0;
/// Number of interval bisections used to refine each root (clause
/// 3.2.3, "4× interval bisection").
const N_BISECT: usize = 4;
/// Half the LP order — the number of roots each of `F1` / `F2` carries.
const NC: usize = M / 2;

/// Evaluates the Chebyshev-series image `G(x)` of a symmetric
/// half-polynomial with coefficients `f[0 … 5]` at `x = cos ω`, via the
/// Clenshaw recurrence:
///
/// `G(x) = ½·f[5] + Σ_{k=1}^{5} f[5−k]·T_k(x)` evaluated as
/// `G = g0 + x·b1 − b2` with `g0 = ½·f[5]`, `g_k = f[5−k]`.
fn cheb(x: f64, f: &[f64; NC + 1]) -> f64 {
    let mut b1 = 0.0f64; // b_{k+1}
    let mut b2 = 0.0f64; // b_{k+2}
                         // k = 5 down to 1, with g_k = f[5 − k].
    for k in (1..=NC).rev() {
        let g_k = f[NC - k];
        let b0 = 2.0 * x * b1 - b2 + g_k;
        b2 = b1;
        b1 = b0;
    }
    0.5 * f[NC] + x * b1 - b2
}

/// Builds the `f1[0 … 5]` (sum) and `f2[0 … 5]` (difference) Chebyshev
/// half-polynomial coefficients from the LP coefficients `a_i`
/// (`a_in[i−1] = a_i`, `a_0 = 1` implicit) via the divide-out recursion.
fn build_polys(a_in: &[f32; M]) -> ([f64; NC + 1], [f64; NC + 1]) {
    // a[0] = 1, a[1..=M] from a_in.
    let mut a = [0.0f64; M + 1];
    a[0] = 1.0;
    for i in 1..=M {
        a[i] = f64::from(a_in[i - 1]);
    }
    let mut f1 = [0.0f64; NC + 1];
    let mut f2 = [0.0f64; NC + 1];
    f1[0] = 1.0;
    f2[0] = 1.0;
    for i in 1..=NC {
        f1[i] = (a[i] + a[M + 1 - i]) - f1[i - 1];
        f2[i] = (a[i] - a[M + 1 - i]) + f2[i - 1];
    }
    (f1, f2)
}

/// §3.2.3 LP→LSP conversion. Converts the 10 LP coefficients `a_i` of
/// `A(z) = 1 + Σ a_i z^{−i}` (layout `a_in[i−1] = a_i`) into the 10
/// cosine-domain LSPs `q_i = cos(ω_i)` (layout `q_out[i−1] = q_i`,
/// decreasing in value / increasing in `ω`).
///
/// Returns `None` if the grid search fails to locate all ten roots
/// (which does not happen for a stable minimum-phase `A(z)`; the caller
/// then keeps the previous frame's LSPs, per clause 3.2.3).
#[must_use]
pub fn lp_to_lsp(a_in: &[f32; M]) -> Option<[f32; M]> {
    let (f1, f2) = build_polys(a_in);

    let mut q = [0.0f64; M];
    let mut found = 0usize;
    // `true` → currently searching F1 (odd LSPs), else F2.
    let mut on_f1 = true;

    let eval = |x: f64, on_f1: bool| -> f64 {
        if on_f1 {
            cheb(x, &f1)
        } else {
            cheb(x, &f2)
        }
    };

    let mut x_lo = f64::from(LSF_SEARCH_GRID_COS_Q15[0]) / Q15;
    let mut y_lo = eval(x_lo, on_f1);

    for &grid_q15 in &LSF_SEARCH_GRID_COS_Q15[1..] {
        let x_hi = f64::from(grid_q15) / Q15;
        let mut y_hi = eval(x_hi, on_f1);

        if y_lo * y_hi <= 0.0 {
            // Sign change in [x_hi, x_lo]: refine by bisection.
            let (mut xl, mut yl) = (x_lo, y_lo);
            let (mut xh, mut yh) = (x_hi, y_hi);
            for _ in 0..N_BISECT {
                let xm = 0.5 * (xl + xh);
                let ym = eval(xm, on_f1);
                if yl * ym <= 0.0 {
                    xh = xm;
                    yh = ym;
                } else {
                    xl = xm;
                    yl = ym;
                }
            }
            // Linear interpolation for the root inside the refined bracket.
            let denom = yh - yl;
            let root = if denom.abs() > f64::EPSILON {
                xl - yl * (xh - xl) / denom
            } else {
                0.5 * (xl + xh)
            };
            q[found] = root.clamp(-1.0, 1.0);
            found += 1;
            if found == M {
                break;
            }
            // Toggle polynomial and restart the running low endpoint at
            // the current grid point re-evaluated with the new
            // polynomial (`x_hi` already holds this grid point).
            on_f1 = !on_f1;
            y_hi = eval(x_hi, on_f1);
        }
        x_lo = x_hi;
        y_lo = y_hi;
    }

    if found != M {
        return None;
    }
    Some(std::array::from_fn(|i| q[i] as f32))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::levinson::levinson;
    use crate::lsp_to_lp::lsp_to_lp;

    /// Ten roots found, strictly decreasing in value (increasing ω), all
    /// inside `(−1, 1)`, for a stable filter.
    fn assert_valid_lsp(q: &[f32; M]) {
        for w in q.windows(2) {
            assert!(
                w[0] > w[1],
                "LSPs must be strictly decreasing (ordered ω): {:?}",
                q
            );
        }
        for &qi in q {
            assert!(qi > -1.0 && qi < 1.0, "LSP {qi} out of (−1,1)");
        }
    }

    /// Round-trip: for a stable filter derived from a real
    /// autocorrelation, `lp_to_lsp` followed by the decoder's
    /// `lsp_to_lp` reproduces the original LP coefficients to grid /
    /// bisection precision.
    #[test]
    fn roundtrip_lp_lsp_lp() {
        // A voiced-like autocorrelation (decaying oscillation).
        let r: [f64; M + 1] =
            std::array::from_fn(|k| (0.92f64).powi(k as i32) * (0.35 * k as f64).cos());
        let a = levinson(&r).a;
        let q = lp_to_lsp(&a).expect("stable filter must yield 10 LSPs");
        assert_valid_lsp(&q);
        let a_rt = lsp_to_lp(&q);
        for (i, (&orig, &rt)) in a.iter().zip(a_rt.iter()).enumerate() {
            assert!(
                (orig - rt).abs() < 2e-3,
                "a[{i}] round-trip: {orig} vs {rt}"
            );
        }
    }

    /// A second, spectrally different filter (low-pass, all-positive
    /// autocorrelation) also round-trips.
    #[test]
    fn roundtrip_lowpass() {
        let r: [f64; M + 1] = std::array::from_fn(|k| (0.8f64).powi(k as i32));
        let a = levinson(&r).a;
        let q = lp_to_lsp(&a).expect("10 LSPs");
        assert_valid_lsp(&q);
        let a_rt = lsp_to_lp(&q);
        for (&orig, &rt) in a.iter().zip(a_rt.iter()) {
            assert!((orig - rt).abs() < 2e-3);
        }
    }

    /// The flat spectrum (all-zero predictor `A(z) = 1`) has its LSPs at
    /// the evenly spaced grid `ω_i = iπ/11` → `q_i = cos(iπ/11)`
    /// (clause 4.3 Table 9 default LSP set). The forward transform of the
    /// zero predictor reproduces that set.
    #[test]
    fn flat_predictor_uniform_lsps() {
        let a = [0.0f32; M];
        let q = lp_to_lsp(&a).expect("10 LSPs for flat spectrum");
        assert_valid_lsp(&q);
        for (i, &qi) in q.iter().enumerate() {
            let expected = ((i + 1) as f64 * std::f64::consts::PI / (M as f64 + 1.0)).cos();
            assert!(
                (f64::from(qi) - expected).abs() < 5e-3,
                "flat LSP {i}: {qi} vs cos({}π/11)={expected}",
                i + 1
            );
        }
    }

    /// Interlacing property: the odd-indexed LSPs (`q_1, q_3, …`) are the
    /// roots of `F1`, the even-indexed of `F2`. Evaluating `F1` at an
    /// odd-index LSP must be ~0.
    #[test]
    fn roots_satisfy_polynomials() {
        let r: [f64; M + 1] =
            std::array::from_fn(|k| (0.9f64).powi(k as i32) * (0.25 * k as f64 + 0.5).cos());
        let a = levinson(&r).a;
        let (f1, f2) = build_polys(&a);
        let q = lp_to_lsp(&a).expect("10 LSPs");
        for (i, &qi) in q.iter().enumerate() {
            let on_f1 = i % 2 == 0; // q_1 (i=0) belongs to F1.
            let val = if on_f1 {
                cheb(f64::from(qi), &f1)
            } else {
                cheb(f64::from(qi), &f2)
            };
            assert!(val.abs() < 5e-3, "residual at LSP {i}: {val}");
        }
    }
}
