//! Clause-5 / Table 15 mathematical functions and the extended
//! double-precision operator family.
//!
//! ## The lookup-table functions
//!
//! Table 12 stages three data tables for the mathematical functions —
//! `tablog` (33 entries, "lookup table in base 2 logarithm
//! computation"), `tabpow` (33, "lookup table in 2^x computation") and
//! `tabsqr` (49, "lookup table in inverse square root computation") —
//! compiled here from `tables/basic-op-{log2,pow2,invsqrt}-table.csv`.
//! Their shapes pin the algorithm: 33 = 2^5 + 1 segment endpoints over
//! one octave (5-bit segment index + linear interpolation on the
//! residual Q15 fraction), 49 = 2^5 + 2^4 + 1 endpoints over the
//! two-octave `[0.25, 1)` mantissa range of the inverse square root.
//! The interpolation arithmetic below runs entirely on the Table 11
//! operator set (`l_deposit_h` of the endpoint, `l_msu` of the
//! endpoint difference against the residual fraction), so every
//! rounding matches the 16-bit grid.
//!
//! - [`log2`] — `log2(L_x)` split into a Word16 integer exponent and a
//!   Q15 fraction; the mantissa is normalised with `norm_l` so the
//!   segment index is the 5 bits below the leading one.
//! - [`pow2`] — `2^(exponent + fraction)` as a Word32 integer, the
//!   inverse split: 5 index bits off the top of the Q15 fraction,
//!   interpolation over `tabpow` (Q14 over `[1, 2)`), then a rounding
//!   shift `l_shr_r` by `30 - exponent`.
//! - [`inv_sqrt`] — `1/√L_x` as a Word32 Q30, with the exponent
//!   forced even before halving so the mantissa lands in `[0.25, 1)`.
//!
//! ## The double-precision (DPF) helpers
//!
//! Table 15 names an extended operator file; the conventional
//! double-precision format it denotes represents a Word32 as
//! `hi·2^16 + lo·2` (both Word16), which composes through `l_mult` /
//! `l_mac` without extra shifts:
//!
//! - [`l_extract`] / [`l_comp`] — split / recompose.
//! - [`mpy_32_16`] — `(hi, lo) × n` keeping 32 bits.
//! - [`mpy_32`] — full 32 × 32 keeping the top 32 bits.
//! - [`div_32`] — normalised 32/32 division by reciprocal
//!   approximation plus one Newton refinement step.

use super::ops::{
    div_s, extract_h, extract_l, l_add, l_deposit_h, l_mac, l_msu, l_mult, l_shl, l_shr, l_shr_r,
    l_sub, mult, norm_l, sub,
};
use crate::tables::{INV_SQRT_TABLE_Q15, LOG2_TABLE_Q15, POW2_TABLE_Q15};

/// Base-2 logarithm of a positive Word32.
///
/// Returns `(exponent, fraction)` with `fraction` on the Q15 grid,
/// such that `L_x ≈ 2^(exponent + fraction/32768)`. Non-positive
/// input returns `(0, 0)` (the callers guard the domain; this matches
/// the defensive zero of the table-driven convention).
#[must_use]
pub fn log2(l_x: i32) -> (i16, i16) {
    if l_x <= 0 {
        return (0, 0);
    }
    let norm = norm_l(l_x);
    let normalized = l_shl(l_x, norm);
    let exponent = sub(30, norm);

    // After normalisation bit 30 is set; drop 9 bits so the high word
    // holds `32 + segment` and the low word the residual fraction.
    let shifted = l_shr(normalized, 9);
    let i = usize::try_from(extract_h(shifted) - 32).expect("normalised segment in 0..32");
    let a = extract_l(l_shr(shifted, 1)) & 0x7fff;

    let mut l_y = l_deposit_h(LOG2_TABLE_Q15[i]);
    let step = sub(LOG2_TABLE_Q15[i], LOG2_TABLE_Q15[i + 1]);
    l_y = l_msu(l_y, step, a);
    (exponent, extract_h(l_y))
}

/// `2^(exponent + fraction/32768)` as a Word32 integer.
///
/// `fraction` is a non-negative Q15 value; `exponent` the integer
/// part. The result saturates on the Word32 grid for `exponent > 30`.
#[must_use]
pub fn pow2(exponent: i16, fraction: i16) -> i32 {
    // Top 5 bits of the Q15 fraction select the tabpow segment, the
    // remaining bits interpolate.
    let l_x = l_mult(fraction, 32); // fraction << 6
    let i = usize::try_from(extract_h(l_x)).expect("pow2 segment in 0..32");
    let a = extract_l(l_shr(l_x, 1)) & 0x7fff;

    let mut l_y = l_deposit_h(POW2_TABLE_Q15[i]);
    let step = sub(POW2_TABLE_Q15[i], POW2_TABLE_Q15[i + 1]);
    l_y = l_msu(l_y, step, a);

    l_shr_r(l_y, sub(30, exponent))
}

/// `1/√L_x` as a Word32 on the Q30 grid.
///
/// Non-positive input returns the saturation ceiling `0x3fffffff`
/// (≈ 1.0 in Q30), matching the defensive convention of the callers
/// (energy terms are non-negative by construction).
#[must_use]
pub fn inv_sqrt(l_x: i32) -> i32 {
    if l_x <= 0 {
        return 0x3fff_ffff;
    }
    let mut exp = norm_l(l_x);
    let mut m = l_shl(l_x, exp);
    exp = sub(30, exp);
    // Force the exponent even so halving it is exact; the mantissa
    // then lands in [0.25, 1) covered by the 49-entry table.
    if (exp & 1) == 0 {
        m = l_shr(m, 1);
    }
    let exp_half = l_add(i32::from(exp >> 1), 1) as i16;

    let shifted = l_shr(m, 9);
    let i = usize::try_from(extract_h(shifted) - 16).expect("inv_sqrt segment in 0..48");
    let a = extract_l(l_shr(shifted, 1)) & 0x7fff;

    let mut l_y = l_deposit_h(INV_SQRT_TABLE_Q15[i]);
    let step = sub(INV_SQRT_TABLE_Q15[i], INV_SQRT_TABLE_Q15[i + 1]);
    l_y = l_msu(l_y, step, a);

    l_shr(l_y, exp_half)
}

/// Split a Word32 into the DPF `(hi, lo)` pair with
/// `L = hi·2^16 + lo·2`.
#[must_use]
pub fn l_extract(l: i32) -> (i16, i16) {
    let hi = extract_h(l);
    let lo = extract_l(l_msu(l_shr(l, 1), hi, 16384));
    (hi, lo)
}

/// Recompose a DPF `(hi, lo)` pair into a Word32.
#[must_use]
pub fn l_comp(hi: i16, lo: i16) -> i32 {
    l_mac(l_deposit_h(hi), lo, 1)
}

/// DPF 32×16 multiplication: `(hi, lo) × n`, keeping 32 bits.
#[must_use]
pub fn mpy_32_16(hi: i16, lo: i16, n: i16) -> i32 {
    l_mac(l_mult(hi, n), mult(lo, n), 1)
}

/// DPF 32×32 multiplication keeping the top 32 bits.
#[must_use]
pub fn mpy_32(hi1: i16, lo1: i16, hi2: i16, lo2: i16) -> i32 {
    let mut l = l_mult(hi1, hi2);
    l = l_mac(l, mult(hi1, lo2), 1);
    l_mac(l, mult(lo1, hi2), 1)
}

/// Normalised 32/32 division `L_num / L_denom` on the Q31 fraction
/// grid, for `0 < L_num < L_denom` with `L_denom` normalised (bit 30
/// set). Reciprocal seed from [`div_s`] plus one refinement pass —
/// the conventional composition of the Table 11 operator set.
#[must_use]
pub fn div_32(l_num: i32, denom_hi: i16, denom_lo: i16) -> i32 {
    // First approximation of 1/denom on the Q14-ish grid.
    let approx = div_s(0x3fff, denom_hi);
    // denom × approx (≈ 0.5 when exact).
    let l_32 = mpy_32_16(denom_hi, denom_lo, approx);
    // Residual 2 − denom·approx (Q30-domain Newton step).
    let residual = l_sub(0x7fff_ffff, l_32);
    let (r_hi, r_lo) = l_extract(residual);
    // Refined reciprocal: approx × residual.
    let recip = mpy_32_16(r_hi, r_lo, approx);
    // num × recip, renormalised.
    let (n_hi, n_lo) = l_extract(l_num);
    let (c_hi, c_lo) = l_extract(recip);
    let q = mpy_32(n_hi, n_lo, c_hi, c_lo);
    l_shl(q, 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Exact powers of two: `log2(2^k)` must produce exponent `k`,
    /// fraction 0, for every representable power.
    #[test]
    fn log2_exact_powers() {
        for k in 0..31 {
            let (e, f) = log2(1i32 << k);
            assert_eq!((e, f), (k as i16, 0), "log2(2^{k})");
        }
    }

    /// `pow2` inverts `log2` on the integer grid: round-tripping a
    /// spread of values reproduces them within the linear-interpolation
    /// tolerance of one table segment (relative error < 2^-11).
    #[test]
    fn pow2_inverts_log2() {
        for &v in &[
            3i32,
            7,
            40,
            100,
            1000,
            12345,
            65536,
            1 << 20,
            0x7fff_0000,
            0x1234_5678,
        ] {
            let (e, f) = log2(v);
            let back = pow2(e, f);
            let err = (f64::from(back) - f64::from(v)).abs() / f64::from(v);
            assert!(err < 5e-4, "pow2(log2({v})) = {back} (rel err {err})");
        }
    }

    /// `pow2` reference points: integer exponents with zero fraction
    /// are exact; the half-step fraction lands on √2 within one LSB of
    /// the interpolated grid.
    #[test]
    fn pow2_reference_points() {
        assert_eq!(pow2(0, 0), 1);
        assert_eq!(pow2(10, 0), 1024);
        assert_eq!(pow2(30, 0), 1 << 30);
        let sqrt2 = pow2(15, 16384); // 2^15.5 = 46340.95
        assert!(
            (f64::from(sqrt2) - 46340.95).abs() < 4.0,
            "2^15.5 = {sqrt2}"
        );
    }

    /// `inv_sqrt` reference points on the Q30 grid: exact powers of
    /// four have exact reciprocal square roots.
    #[test]
    fn inv_sqrt_reference_points() {
        // 1/sqrt(4) = 0.5 → Q30 0x2000_0000, within interpolation LSBs.
        let v = inv_sqrt(4);
        assert!(
            (f64::from(v) / (1i64 << 30) as f64 - 0.5).abs() < 1e-3,
            "inv_sqrt(4) = {v:#x}"
        );
        for &(x, expect) in &[(16i32, 0.25f64), (64, 0.125), (100, 0.1), (10000, 0.01)] {
            let got = f64::from(inv_sqrt(x)) / (1i64 << 30) as f64;
            assert!(
                (got - expect).abs() / expect < 1e-3,
                "inv_sqrt({x}) = {got} expected {expect}"
            );
        }
        assert_eq!(inv_sqrt(0), 0x3fff_ffff);
        assert_eq!(inv_sqrt(-5), 0x3fff_ffff);
    }

    /// DPF split/recompose is lossless for even low bits (the format
    /// carries 31 significant bits: `hi·2^16 + lo·2`).
    #[test]
    fn dpf_extract_comp_roundtrip() {
        for &v in &[0i32, 2, -2, 0x1234_5678 & !1, -0x1234_5678 & !1, 1 << 30] {
            let (hi, lo) = l_extract(v);
            assert_eq!(l_comp(hi, lo), v & !1, "roundtrip {v:#x}");
        }
    }

    /// `mpy_32_16` multiplies a DPF value by a Q15 factor: scaling a
    /// known Word32 by 0.5 (16384) halves it within one LSB.
    #[test]
    fn mpy_32_16_halves_with_q15_half() {
        for &v in &[0x0004_0000i32, 0x1234_5678, 0x7000_0000] {
            let (hi, lo) = l_extract(v);
            let half = mpy_32_16(hi, lo, 16384);
            assert!(
                (i64::from(half) - i64::from(v) / 2).abs() <= 2,
                "half of {v:#x} = {half:#x}"
            );
        }
    }

    /// `mpy_32` squares a DPF value: (0.5)² = 0.25 on the Q31 grid.
    #[test]
    fn mpy_32_squares_half() {
        let half = 0x4000_0000i32; // 0.5 in Q31
        let (hi, lo) = l_extract(half);
        let sq = mpy_32(hi, lo, hi, lo);
        assert!((i64::from(sq) - 0x2000_0000).abs() <= 2, "0.5^2 = {sq:#x}");
    }

    /// `div_32` divides normalised fractions: 0.25 / 0.5 = 0.5 within
    /// the refinement tolerance.
    #[test]
    fn div_32_quarter_over_half() {
        let denom = 0x4000_0000i32; // 0.5 in Q31, normalised
        let (dh, dl) = l_extract(denom);
        let q = div_32(0x2000_0000, dh, dl);
        let got = f64::from(q) / (1u64 << 31) as f64;
        assert!((got - 0.5).abs() < 1e-4, "0.25/0.5 = {got}");
    }

    /// The three staged tables have the endpoint values their roles
    /// require (already shape-asserted at build time; the numeric ends
    /// pin the Q formats).
    #[test]
    fn table_endpoints() {
        assert_eq!(LOG2_TABLE_Q15[0], 0);
        assert_eq!(LOG2_TABLE_Q15[32], 32767);
        assert_eq!(POW2_TABLE_Q15[0], 16384);
        assert_eq!(POW2_TABLE_Q15[32], 32767);
        assert_eq!(INV_SQRT_TABLE_Q15[0], 32767);
        assert_eq!(INV_SQRT_TABLE_Q15[48], 16384);
    }
}
