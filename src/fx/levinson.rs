//! §3.2.2 **fixed-point Levinson-Durbin recursion** — eqs (8)/(9) on
//! the clause-5 Word16/Word32 operator grid in the double-precision
//! (hi, lo) format.
//!
//! The float module [`crate::levinson`] emulates the fixed-point grid
//! by rounding each order's reflection coefficient to Q31 and each LP
//! coefficient to Q27; a genuine DPF pipeline instead *truncates*
//! inside every 32×32 product (`mpy_32` drops the `lo·lo` cross term)
//! and divides through the Newton-refined [`div_32`], so its rounding
//! pattern differs from the emulation by a few LSB per order — the
//! same family as the encoder's residual front-end ω divergence.
//!
//! ## Numeric layout
//!
//! - `r'(k)` arrive as normalised DPF pairs from
//!   [`crate::fx::analysis`] (common scale `Q(n)`, `r(0)` in
//!   `[2^30, 2^31)`).
//! - LP coefficients are carried as Word32 on the **Q27** grid
//!   (`a_j·2^27`, i.e. hi half on Q12 — 3 integer bits of headroom
//!   for the transient over-unity coefficients of the recursion).
//! - Reflection coefficients are carried on the **Q31** grid
//!   (`k_i·2^31`, hi half on Q15).
//! - The prediction-error energy `E^{(i)}` is kept **renormalised**
//!   after every order (`norm_l` shift, cumulative exponent
//!   `alp_exp`), so the eq (8) division always sees a normalised
//!   denominator; the quotient is re-referred to the `Q(n)` grid by
//!   shifting the accumulated exponent back in.
//!
//! ## Order-`i` schedule
//!
//! ```text
//! t0   = Σ_{j=1}^{i−1} mpy_32(a_j, r'(i−j))      (Q(n−4))
//! t0   = l_shl(t0, 4) + r'(i)                    (Q(n))
//! k_i  = −t0 / E^{(i−1)}     via div_32 + l_shl(·, alp_exp)  (Q31)
//! a_j += mpy_32(k_i, a_{i−j})     j = 1 … i−1    (Q27, snapshot)
//! a_i  = l_shr(k_i, 4)                           (Q27)
//! E    = E − mpy_32(k_i², E) ; renormalise, alp_exp += shift
//! ```
//!
//! The recursion reports failure (caller keeps the previous frame's
//! LP set, the clause-3.2.3 fallback) when the autocorrelation is not
//! positive definite: `|t0| > E` (a reflection magnitude ≥ 1) or a
//! vanishing energy.

use crate::fx::analysis::Dpf;
use crate::fx::dsp::{div_32, l_comp, l_extract, mpy_32};
use crate::fx::ops::{l_abs, l_negate, l_shl, l_shr, l_sub, norm_l, round};
use crate::lp_analysis::N_LAGS;
use crate::tables::M;

/// `1.0` on the Q27 LP-coefficient grid.
const ONE_Q27: i32 = 0x0800_0000;

/// Result of the fixed-point §3.2.2 recursion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LevinsonFxResult {
    /// LP coefficients `a_1 … a_M` rounded to Q12
    /// (`A(z) = 1 + Σ a_i z^{−i}`).
    pub a_q12: [i16; M],
    /// The same coefficients on the full Q27 DPF grid (hi, lo) —
    /// the precision later fixed-point stages consume.
    pub a_q27: [Dpf; M],
    /// Reflection coefficients `k_1 … k_M` rounded to Q15 (the §3.3
    /// perceptual weighting consumes the first two).
    pub rc_q15: [i16; M],
}

/// Runs eqs (8)/(9) on normalised DPF autocorrelations (as produced
/// by [`crate::fx::analysis::analyze_window_fx`]). Returns `None`
/// when the recursion is ill-conditioned (non-positive-definite
/// input); the caller keeps the previous frame's LP set per the
/// clause-3.2.3 fallback.
#[must_use]
pub fn levinson_fx(r: &[Dpf; N_LAGS]) -> Option<LevinsonFxResult> {
    // a on Q27; a[0] = 1 implicit (kept for the update indexing).
    let mut a = [0i32; M + 1];
    a[0] = ONE_Q27;
    let mut rc = [0i16; M];

    // E^{(0)} = r'(0), already normalised; alp_exp tracks how far the
    // stored E has been up-shifted relative to the Q(n) grid.
    let mut e = l_comp(r[0].0, r[0].1);
    let mut alp_exp: i16 = 0;
    if e <= 0 {
        return None;
    }

    for i in 1..=M {
        // eq (8) numerator: t0 = r'(i) + Σ a_j·r'(i−j) on Q(n).
        let mut t0: i32 = 0;
        for j in 1..i {
            let (ah, al) = l_extract(a[j]);
            t0 = t0.wrapping_add(mpy_32(ah, al, r[i - j].0, r[i - j].1));
        }
        let t0 = l_shl(t0, 4).wrapping_add(l_comp(r[i].0, r[i].1));

        // eq (8) division: k = −t0/E. On the Q(n) grid the true energy
        // is the stored E down-shifted by alp_exp, so |t0| reaching it
        // means a reflection magnitude ≥ 1 (not positive definite).
        // The check also guarantees div_32's `num < den` contract
        // (l_shr(e, alp_exp) ≤ e).
        let t1 = l_abs(t0);
        if t1 >= l_shr(e, alp_exp) {
            return None;
        }
        let (e_hi, e_lo) = l_extract(e);
        let quot = div_32(t1, e_hi, e_lo);
        // Re-refer the quotient from the stored-E grid to Q(n).
        let mut k = l_shl(quot, alp_exp);
        if t0 > 0 {
            k = l_negate(k);
        }
        let (k_hi, k_lo) = l_extract(k);
        if k_hi.unsigned_abs() >= 32750 {
            // Reflection magnitude at the saturation rim — treat as
            // ill-conditioned rather than committing a clipped k.
            return None;
        }
        rc[i - 1] = round(k);

        // eq (9) coefficient update from an order-(i−1) snapshot.
        let prev = a;
        for j in 1..i {
            let (ph, pl) = l_extract(prev[i - j]);
            a[j] = prev[j].wrapping_add(mpy_32(k_hi, k_lo, ph, pl));
        }
        a[i] = l_shr(k, 4);

        // Energy update E ← E·(1 − k²), renormalised.
        let k2 = mpy_32(k_hi, k_lo, k_hi, k_lo);
        let (k2_hi, k2_lo) = l_extract(l_abs(k2));
        let e_drop = mpy_32(k2_hi, k2_lo, e_hi, e_lo);
        e = l_sub(e, e_drop);
        if e <= 0 {
            return None;
        }
        let sh = norm_l(e);
        e = l_shl(e, sh);
        alp_exp += sh;
    }

    // Q27 → Q12 rounding: one left shift lands the high word on Q12.
    let a_q12: [i16; M] = std::array::from_fn(|j| round(l_shl(a[j + 1], 1)));
    let a_q27: [Dpf; M] = std::array::from_fn(|j| l_extract(a[j + 1]));
    Some(LevinsonFxResult {
        a_q12,
        a_q27,
        rc_q15: rc,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fx::analysis::dpf_to_f64;
    use crate::levinson::levinson;

    /// Builds a normalised DPF autocorrelation from real-valued lags.
    fn to_dpf(r: &[f64; N_LAGS]) -> [Dpf; N_LAGS] {
        let mut scale = 0i32;
        while r[0] * (2f64).powi(scale) >= (2f64).powi(31) {
            scale -= 1;
        }
        while r[0] * (2f64).powi(scale) < (2f64).powi(30) {
            scale += 1;
        }
        let f = (2f64).powi(scale);
        std::array::from_fn(|k| {
            let w = (r[k] * f).floor() as i32;
            l_extract(w)
        })
    }

    /// AR(1): `r(k) = ρ^k` must recover `a_1 ≈ −ρ`, higher taps ≈ 0.
    #[test]
    fn ar1_recovery() {
        let rho = 0.8f64;
        let r: [f64; N_LAGS] = std::array::from_fn(|k| rho.powi(k as i32));
        let res = levinson_fx(&to_dpf(&r)).expect("stable");
        assert!(
            (f64::from(res.a_q12[0]) / 4096.0 + rho).abs() < 2e-3,
            "a_1 = {} expected ≈ {}",
            res.a_q12[0],
            -rho * 4096.0
        );
        for &ai in &res.a_q12[1..] {
            assert!(ai.abs() <= 4, "higher tap should vanish: {ai}");
        }
        assert!((f64::from(res.rc_q15[0]) / 32768.0 + rho).abs() < 2e-3);
    }

    /// Flat spectrum: `r(k) = 0` for `k ≥ 1` gives the zero predictor.
    #[test]
    fn flat_predictor() {
        let mut r = [0.0f64; N_LAGS];
        r[0] = 5e8;
        let res = levinson_fx(&to_dpf(&r)).expect("stable");
        assert_eq!(res.a_q12, [0i16; M]);
        assert_eq!(res.rc_q15, [0i16; M]);
    }

    /// Against the float-emulated recursion the Q12 coefficients agree
    /// within a few LSB across a spread of stable spectra (the DPF
    /// truncation pattern differs from the emulation's
    /// round-to-grid — that residual is what the corpus pin measures).
    #[test]
    fn tracks_float_emulation() {
        for (decay, freq, phase) in [
            (0.9f64, 0.3f64, 0.0f64),
            (0.85, 0.2, 1.0),
            (0.92, 0.35, 0.0),
            (0.8, 0.0, 0.0),
            (0.95, 0.7, 0.0),
        ] {
            let r: [f64; N_LAGS] =
                std::array::from_fn(|k| decay.powi(k as i32) * (freq * k as f64 + phase).cos());
            let dpf = to_dpf(&r);
            let fx = levinson_fx(&dpf).expect("stable");
            let float = levinson(&dpf_to_f64(&dpf));
            for j in 0..M {
                let f_q12 = (float.a[j] * 4096.0).round() as i32;
                let d = (i32::from(fx.a_q12[j]) - f_q12).abs();
                assert!(
                    d <= 1,
                    "decay {decay} freq {freq}: a[{j}] fx {} vs float {f_q12}",
                    fx.a_q12[j]
                );
            }
        }
    }

    /// The Q27 DPF output recomposes to the Q12 rounding it reports.
    #[test]
    fn q27_consistent_with_q12() {
        let r: [f64; N_LAGS] =
            std::array::from_fn(|k| (0.9f64).powi(k as i32) * (0.25 * k as f64).cos());
        let res = levinson_fx(&to_dpf(&r)).expect("stable");
        for j in 0..M {
            let w = l_comp(res.a_q27[j].0, res.a_q27[j].1);
            assert_eq!(round(l_shl(w, 1)), res.a_q12[j], "tap {j}");
        }
    }

    /// A phase-shifted damped cosine whose spectrum goes negative is
    /// not a valid autocorrelation; the fixed recursion reports the
    /// ill-conditioning honestly (where the float emulation silently
    /// early-stops with the offending reflection already committed).
    #[test]
    fn non_pd_phase_shifted_cosine_bails() {
        let r: [f64; N_LAGS] =
            std::array::from_fn(|k| (0.95f64).powi(k as i32) * (0.7 * k as f64 + 0.4).cos());
        assert!(levinson_fx(&to_dpf(&r)).is_none());
    }

    /// A non-positive-definite input reports failure instead of a
    /// clipped reflection coefficient.
    #[test]
    fn non_positive_definite_bails() {
        let mut r = [0.0f64; N_LAGS];
        r[0] = 2e9 / 2.0;
        r[1] = 1.6e9; // |r(1)| > r(0) → |k_1| > 1.
        assert!(levinson_fx(&to_dpf(&r)).is_none());
    }

    /// Degenerate zero energy reports failure.
    #[test]
    fn zero_energy_bails() {
        let r = [(0i16, 0i16); N_LAGS];
        assert!(levinson_fx(&r).is_none());
    }
}
