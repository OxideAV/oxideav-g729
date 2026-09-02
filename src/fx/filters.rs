//! Encoder-side fixed-point filter primitives shared by the clause-3
//! stages: the eq (36) LP inverse filter (`residu`), the γ-weighting
//! of an LP coefficient set (`weight_az`, the `A(z/γ)` construction
//! of eq (27)), the zero-state convolution with an impulse response
//! (eqs (44)/(64)), and a zero-memory synthesis run.
//!
//! ## Number grids
//!
//! - LP coefficients `a_i` — Q12 with `a_0 = 4096`, the grid the
//!   §3.2.2 recursion ([`crate::fx::levinson`]) and the §3.2.6
//!   conversion ([`crate::fx::lsp`]) both land on.
//! - Signals (speech, residual, weighted speech, target) — Q0, the
//!   halved 16-bit PCM grid the §3.1 pre-processor produces.
//! - Impulse response `h(n)` — Q12 (`h(0)` of the identity filter is
//!   `4096`).
//!
//! The FIR/IIR recursions accumulate `Q12 × Q0` products on the
//! doubled Q13 grid of `l_mult`/`l_mac`, shift up 3 and round onto Q0
//! — the same landing [`crate::fx::excitation::syn_filt`] uses on the
//! decode side, so an encoder-side residual followed by the decoder's
//! synthesis is an exact inverse pair on the integer grid whenever the
//! accumulator stays in range.

use crate::fx::ops::{extract_h, l_mac, l_msu, l_mult, l_shl, mult, mult_r, round};
use crate::tables::M;

/// Samples per subframe.
pub const L_SUBFR: usize = 40;

/// Unstated Q0 landing latitude of the encoder-side filter runs
/// (round-to-nearest vs truncation of the `<< 3` accumulator).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FilterLatitude {
    /// Truncate the eq (36) residual landing.
    pub residu_trunc: bool,
    /// Truncate the all-pole (`1/Â(z)`, `1/A(z/γ₂)`) landings.
    pub syn_trunc: bool,
    /// Truncate the eq (44) convolution landing.
    pub conv_trunc: bool,
    /// Truncating (`mult`) instead of rounding (`mult_r`) products in
    /// the γ-weighting of the LP coefficients.
    pub weight_trunc: bool,
}

/// The `<< 3` Q0 landing of a Q13 accumulator.
#[inline]
fn land(acc: i32, trunc: bool) -> i16 {
    let v = l_shl(acc, 3);
    if trunc {
        extract_h(v)
    } else {
        round(v)
    }
}

/// eq (27) / §3.3: bandwidth-expands an LP coefficient set,
/// `a'_i = γ^i · a_i`, with `γ` on Q15. The running power `γ^i` is
/// kept on Q15 and refreshed by a rounding multiply each tap, so the
/// tenth-order term carries the rounded — not truncated — product.
#[must_use]
pub fn weight_az(a: &[i16; M + 1], gamma_q15: i16) -> [i16; M + 1] {
    weight_az_lat(a, gamma_q15, false)
}

/// [`weight_az`] with a selectable product rounding.
#[must_use]
pub fn weight_az_lat(a: &[i16; M + 1], gamma_q15: i16, trunc: bool) -> [i16; M + 1] {
    let mut out = [0i16; M + 1];
    out[0] = a[0];
    let mut fac = gamma_q15;
    for i in 1..=M {
        if trunc {
            out[i] = mult(a[i], fac);
            fac = mult(fac, gamma_q15);
        } else {
            out[i] = mult_r(a[i], fac);
            fac = mult_r(fac, gamma_q15);
        }
    }
    out
}

/// eq (36): the LP residual `r(n) = Σ_{i=0}^{M} a_i · x(n − i)` of one
/// subframe. `x_hist` holds the samples preceding the subframe
/// (`x_hist[M − 1]` = `x(−1)`, most recent last) and `x` the subframe
/// itself.
#[must_use]
pub fn residu(a: &[i16; M + 1], x_hist: &[i16; M], x: &[i16; L_SUBFR]) -> [i16; L_SUBFR] {
    residu_lat(a, x_hist, x, false)
}

/// [`residu`] with a selectable landing.
#[must_use]
pub fn residu_lat(
    a: &[i16; M + 1],
    x_hist: &[i16; M],
    x: &[i16; L_SUBFR],
    trunc: bool,
) -> [i16; L_SUBFR] {
    let mut out = [0i16; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = l_mult(x[n], a[0]);
        for i in 1..=M {
            let s = if n >= i {
                x[n - i]
            } else {
                x_hist[M - (i - n)]
            };
            acc = l_mac(acc, a[i], s);
        }
        out[n] = land(acc, trunc);
    }
    out
}

/// eq (77)-shaped all-pole run `y(n) = x(n) − Σ_{i=1}^{M} a_i · y(n−i)`
/// with an explicit output history (`y_hist[M − 1]` = `y(−1)`),
/// returning the subframe output without touching the caller's
/// memory. The same Q13 accumulator / `<< 3` / round landing as
/// [`crate::fx::excitation::syn_filt`].
#[must_use]
pub fn syn_filt_mem(a: &[i16; M + 1], y_hist: &[i16; M], x: &[i16; L_SUBFR]) -> [i16; L_SUBFR] {
    syn_filt_mem_lat(a, y_hist, x, false)
}

/// [`syn_filt_mem`] with a selectable landing.
#[must_use]
pub fn syn_filt_mem_lat(
    a: &[i16; M + 1],
    y_hist: &[i16; M],
    x: &[i16; L_SUBFR],
    trunc: bool,
) -> [i16; L_SUBFR] {
    let mut out = [0i16; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = l_mult(x[n], a[0]);
        for i in 1..=M {
            let s = if n >= i {
                out[n - i]
            } else {
                y_hist[M - (i - n)]
            };
            acc = l_msu(acc, a[i], s);
        }
        out[n] = land(acc, trunc);
    }
    out
}

/// Zero-state convolution with an impulse response on Q12:
/// `y(n) = Σ_{i=0}^{n} x(i) · h(n − i)`, the eq (44) / eq (64)
/// filtering of a Q0 vector through `h`. The Q13 (doubled Q12 × Q0)
/// accumulator shifts up 3 and rounds onto Q0.
#[must_use]
pub fn convolve_h_q0(x: &[i16; L_SUBFR], h_q12: &[i16; L_SUBFR]) -> [i16; L_SUBFR] {
    convolve_h_q0_lat(x, h_q12, false)
}

/// [`convolve_h_q0`] with a selectable landing.
#[must_use]
pub fn convolve_h_q0_lat(
    x: &[i16; L_SUBFR],
    h_q12: &[i16; L_SUBFR],
    trunc: bool,
) -> [i16; L_SUBFR] {
    let mut out = [0i16; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = 0i32;
        for i in 0..=n {
            acc = l_mac(acc, x[i], h_q12[n - i]);
        }
        out[n] = land(acc, trunc);
    }
    out
}

/// Zero-state convolution of a Q13 codevector with the Q12 impulse
/// response, landing on Q12: `z(n) = Σ c(i) · h(n − i)` (eq (64)).
/// The doubled Q26 accumulator shifts up 2 and rounds to the high
/// word.
#[must_use]
pub fn convolve_code_q12(c_q13: &[i16; L_SUBFR], h_q12: &[i16; L_SUBFR]) -> [i16; L_SUBFR] {
    let mut out = [0i16; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = 0i32;
        for i in 0..=n {
            if c_q13[i] != 0 {
                acc = l_mac(acc, c_q13[i], h_q12[n - i]);
            }
        }
        out[n] = round(l_shl(acc, 2));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `weight_az` with γ = 1 is the identity; with γ = 0.5 the taps
    /// halve geometrically.
    #[test]
    fn weight_az_scales_geometrically() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2048;
        a[2] = 1024;
        a[3] = -512;
        let same = weight_az(&a, 32767);
        assert_eq!(same[0], 4096);
        assert!((i32::from(same[1]) + 2048).abs() <= 1);
        let half = weight_az(&a, 16384);
        assert_eq!(half[1], -1024);
        assert_eq!(half[2], 256);
        assert_eq!(half[3], -64);
    }

    /// The residual of a synthesised signal recovers the excitation
    /// (eq (36) inverts eq (77) on the integer grid).
    #[test]
    fn residu_inverts_synthesis() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2867; // −0.7
        a[2] = 1638; // 0.4
        a[3] = -819; // −0.2
        let exc: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 17 % 41) as i16) * 50 - 1000);
        let hist = [0i16; M];
        let s = syn_filt_mem(&a, &hist, &exc);
        let r = residu(&a, &hist, &s);
        for n in 0..L_SUBFR {
            assert!(
                (i32::from(r[n]) - i32::from(exc[n])).abs() <= 1,
                "n={n}: {} vs {}",
                r[n],
                exc[n]
            );
        }
    }

    /// Convolving a unit impulse with `h` reproduces `h` on the Q0
    /// grid (scaled by the Q12 unit).
    #[test]
    fn convolution_identity() {
        let h: [i16; L_SUBFR] = std::array::from_fn(|n| (4096 >> n.min(12)) as i16);
        let mut x = [0i16; L_SUBFR];
        x[0] = 1000;
        let y = convolve_h_q0(&x, &h);
        assert_eq!(y[0], 1000);
        assert_eq!(y[1], 500);
        assert_eq!(y[2], 250);
        let mut c = [0i16; L_SUBFR];
        c[3] = 8192;
        let z = convolve_code_q12(&c, &h);
        assert_eq!(z[3], 4096);
        assert_eq!(z[4], 2048);
        assert_eq!(z[0], 0);
    }
}
