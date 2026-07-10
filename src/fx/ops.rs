//! The clause-5 / Table 11 basic operator set.
//!
//! Table 10 pins the two data types (`Word16` = signed 16-bit,
//! `Word32` = signed 32-bit two's complement); Table 11 enumerates the
//! operators and their one-line roles ("Short addition", "Multiply and
//! accumulate", …). The semantics implemented here are the
//! conventional saturating two's-complement fixed-point semantics that
//! operator set denotes:
//!
//! - 16-bit results saturate to `[-32768, 32767]`, 32-bit results to
//!   `[-2147483648, 2147483647]`.
//! - `mult` is the Q15 fractional product `(a·b) >> 15` with the
//!   `(-1)·(-1)` corner saturated to `0x7fff`.
//! - `l_mult` is the doubled product `(a·b) << 1` (Q15×Q15 → Q31),
//!   saturating only on the `(-32768)²` corner.
//! - `round` adds `2^15` and takes the high word (round-half-up on
//!   the Q16 boundary).
//! - shifts saturate on the way left and shift arithmetically (sign
//!   propagating) on the way right; negative shift counts reverse
//!   direction.
//! - `norm_s` / `norm_l` return the left-shift count that normalises
//!   the operand (0 for zero input).
//! - `div_s` is the non-negative fractional division `a/b` in Q15,
//!   defined for `0 <= a <= b`, `b > 0`.
//!
//! Where a G.729 clause depends on one of these corners (e.g. the
//! §4.1.6 synthesis-filter saturation on the OVERFLOW vector), the
//! conformance harness exercises it end-to-end; the unit tests below
//! pin every corner directly.

/// Saturate a wide intermediate onto the Word16 grid.
#[inline]
#[must_use]
pub fn sature(x: i32) -> i16 {
    x.clamp(i32::from(i16::MIN), i32::from(i16::MAX)) as i16
}

/// Saturate a 64-bit intermediate onto the Word32 grid.
#[inline]
#[must_use]
pub fn sature32(x: i64) -> i32 {
    x.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
}

/// Short addition with saturation.
#[inline]
#[must_use]
pub fn add(a: i16, b: i16) -> i16 {
    sature(i32::from(a) + i32::from(b))
}

/// Short subtraction with saturation.
#[inline]
#[must_use]
pub fn sub(a: i16, b: i16) -> i16 {
    sature(i32::from(a) - i32::from(b))
}

/// Short absolute value; `abs_s(-32768) = 32767`.
#[inline]
#[must_use]
pub fn abs_s(a: i16) -> i16 {
    if a == i16::MIN {
        i16::MAX
    } else {
        a.abs()
    }
}

/// Short shift left with saturation; negative counts shift right.
#[inline]
#[must_use]
pub fn shl(a: i16, b: i16) -> i16 {
    if b < 0 {
        return shr(a, sature(-i32::from(b)));
    }
    if b >= 15 {
        return match a.cmp(&0) {
            core::cmp::Ordering::Greater => i16::MAX,
            core::cmp::Ordering::Less => i16::MIN,
            core::cmp::Ordering::Equal => 0,
        };
    }
    sature(i32::from(a) << b)
}

/// Short arithmetic shift right; negative counts shift left.
#[inline]
#[must_use]
pub fn shr(a: i16, b: i16) -> i16 {
    if b < 0 {
        return shl(a, sature(-i32::from(b)));
    }
    if b >= 15 {
        return if a < 0 { -1 } else { 0 };
    }
    a >> b
}

/// Q15 fractional multiplication `(a·b) >> 15`, saturating the
/// `(-32768)·(-32768)` corner.
#[inline]
#[must_use]
pub fn mult(a: i16, b: i16) -> i16 {
    sature((i32::from(a) * i32::from(b)) >> 15)
}

/// Long (doubled) multiplication `(a·b) << 1` — Q15 × Q15 → Q31.
/// Saturates only on the `(-32768)²` corner.
#[inline]
#[must_use]
pub fn l_mult(a: i16, b: i16) -> i32 {
    let p = i32::from(a) * i32::from(b);
    if p == 0x4000_0000 {
        i32::MAX
    } else {
        p << 1
    }
}

/// Short negate; `negate(-32768) = 32767`.
#[inline]
#[must_use]
pub fn negate(a: i16) -> i16 {
    if a == i16::MIN {
        i16::MAX
    } else {
        -a
    }
}

/// Extract the high word of a Word32.
#[inline]
#[must_use]
pub fn extract_h(l: i32) -> i16 {
    (l >> 16) as i16
}

/// Extract the low word of a Word32 (no saturation — raw bits).
#[inline]
#[must_use]
pub fn extract_l(l: i32) -> i16 {
    (l & 0xffff) as i16
}

/// Round: add `2^15` (saturating) and extract the high word.
#[inline]
#[must_use]
pub fn round(l: i32) -> i16 {
    extract_h(l_add(l, 0x0000_8000))
}

/// Multiply-accumulate: `l + (a·b) << 1` with both steps saturating.
#[inline]
#[must_use]
pub fn l_mac(l: i32, a: i16, b: i16) -> i32 {
    l_add(l, l_mult(a, b))
}

/// Multiply-subtract: `l - (a·b) << 1` with both steps saturating.
#[inline]
#[must_use]
pub fn l_msu(l: i32, a: i16, b: i16) -> i32 {
    l_sub(l, l_mult(a, b))
}

/// Long addition with saturation.
#[inline]
#[must_use]
pub fn l_add(a: i32, b: i32) -> i32 {
    sature32(i64::from(a) + i64::from(b))
}

/// Long subtraction with saturation.
#[inline]
#[must_use]
pub fn l_sub(a: i32, b: i32) -> i32 {
    sature32(i64::from(a) - i64::from(b))
}

/// Long negate; `l_negate(0x80000000) = 0x7fffffff`.
#[inline]
#[must_use]
pub fn l_negate(a: i32) -> i32 {
    if a == i32::MIN {
        i32::MAX
    } else {
        -a
    }
}

/// Q15 fractional multiplication with rounding:
/// `(a·b + 2^14) >> 15`, saturated.
#[inline]
#[must_use]
pub fn mult_r(a: i16, b: i16) -> i16 {
    sature((i32::from(a) * i32::from(b) + 0x4000) >> 15)
}

/// Long shift left with saturation; negative counts shift right.
#[inline]
#[must_use]
pub fn l_shl(a: i32, b: i16) -> i32 {
    if b <= 0 {
        return l_shr(a, sature(-i32::from(b)));
    }
    let mut out = a;
    for _ in 0..b {
        if out > 0x3fff_ffff {
            return i32::MAX;
        }
        if out < -0x4000_0000 {
            return i32::MIN;
        }
        out <<= 1;
    }
    out
}

/// Long arithmetic shift right; negative counts shift left.
#[inline]
#[must_use]
pub fn l_shr(a: i32, b: i16) -> i32 {
    if b < 0 {
        return l_shl(a, sature(-i32::from(b)));
    }
    if b >= 31 {
        return if a < 0 { -1 } else { 0 };
    }
    a >> b
}

/// Short shift right with rounding: shifts by `b-1`, adds one, then
/// shifts the final bit out.
#[inline]
#[must_use]
pub fn shr_r(a: i16, b: i16) -> i16 {
    if b > 15 {
        return 0;
    }
    if b <= 0 {
        return shr(a, b);
    }
    let pre = shr(a, sub(b, 1));
    shr(add(pre, 1), 1)
}

/// Mac with rounding: `round(l_mac(l, a, b))`.
#[inline]
#[must_use]
pub fn mac_r(l: i32, a: i16, b: i16) -> i16 {
    round(l_mac(l, a, b))
}

/// Msu with rounding: `round(l_msu(l, a, b))`.
#[inline]
#[must_use]
pub fn msu_r(l: i32, a: i16, b: i16) -> i16 {
    round(l_msu(l, a, b))
}

/// Deposit a Word16 into the high word of a Word32.
#[inline]
#[must_use]
pub fn l_deposit_h(a: i16) -> i32 {
    i32::from(a) << 16
}

/// Deposit a Word16 into the low word of a Word32 (sign-extended).
#[inline]
#[must_use]
pub fn l_deposit_l(a: i16) -> i32 {
    i32::from(a)
}

/// Long shift right with rounding.
#[inline]
#[must_use]
pub fn l_shr_r(a: i32, b: i16) -> i32 {
    if b > 31 {
        return 0;
    }
    if b <= 0 {
        return l_shr(a, b);
    }
    let pre = l_shr(a, sub(b, 1));
    l_shr(l_add(pre, 1), 1)
}

/// Long absolute value; `l_abs(0x80000000) = 0x7fffffff`.
#[inline]
#[must_use]
pub fn l_abs(a: i32) -> i32 {
    if a == i32::MIN {
        i32::MAX
    } else {
        a.abs()
    }
}

/// Short norm: the left-shift count that normalises `a` into
/// `[0x4000, 0x7fff]` (positive) or `[0x8000, 0xc000)` (negative);
/// `norm_s(0) = 0`, `norm_s(-1) = 15`.
#[inline]
#[must_use]
pub fn norm_s(a: i16) -> i16 {
    if a == 0 {
        return 0;
    }
    if a == -1 {
        return 15;
    }
    let mut v = if a < 0 { !a } else { a };
    let mut n = 0i16;
    while v < 0x4000 {
        v <<= 1;
        n += 1;
    }
    n
}

/// Long norm: the left-shift count that normalises `a` into
/// `[0x40000000, 0x7fffffff]` (positive) or its negative mirror;
/// `norm_l(0) = 0`, `norm_l(-1) = 31`.
#[inline]
#[must_use]
pub fn norm_l(a: i32) -> i16 {
    if a == 0 {
        return 0;
    }
    if a == -1 {
        return 31;
    }
    let mut v = if a < 0 { !a } else { a };
    let mut n = 0i16;
    while v < 0x4000_0000 {
        v <<= 1;
        n += 1;
    }
    n
}

/// Short fractional division `a/b` in Q15 for `0 <= a <= b`, `b > 0`.
/// Saturates to `0x7fff` when `a == b`.
///
/// # Panics
///
/// Panics if the domain contract (`b > 0`, `0 <= a <= b`) is violated
/// — every G.729 call site establishes it beforehand (the callers
/// clamp or normalise), so a violation is an internal logic error.
#[inline]
#[must_use]
pub fn div_s(a: i16, b: i16) -> i16 {
    assert!(b > 0, "div_s: divisor must be positive, got {b}");
    assert!((0..=b).contains(&a), "div_s: dividend {a} outside [0, {b}]");
    if a == b {
        return i16::MAX;
    }
    // Long division on the Q15 grid: (a << 15) / b, truncating —
    // identical to the conventional 15-step conditional-subtract loop.
    ((i32::from(a) << 15) / i32::from(b)) as i16
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_sub_saturate() {
        assert_eq!(add(30000, 10000), 32767);
        assert_eq!(add(-30000, -10000), -32768);
        assert_eq!(add(1, 2), 3);
        assert_eq!(sub(-30000, 10000), -32768);
        assert_eq!(sub(30000, -10000), 32767);
        assert_eq!(sub(5, 7), -2);
    }

    #[test]
    fn abs_negate_corners() {
        assert_eq!(abs_s(-32768), 32767);
        assert_eq!(abs_s(-5), 5);
        assert_eq!(negate(-32768), 32767);
        assert_eq!(negate(7), -7);
    }

    #[test]
    fn mult_q15_semantics() {
        // 0.5 * 0.5 = 0.25 on the Q15 grid.
        assert_eq!(mult(16384, 16384), 8192);
        // (-1) * (-1) saturates to 0x7fff.
        assert_eq!(mult(-32768, -32768), 32767);
        // Truncation toward -inf for the >> operator on negatives:
        // (-1 * 3) >> 15 == -1 on arithmetic shift.
        assert_eq!(mult(-1, 3), -1);
    }

    #[test]
    fn l_mult_is_doubled_product() {
        assert_eq!(l_mult(16384, 16384), 0x2000_0000);
        assert_eq!(l_mult(-32768, -32768), i32::MAX);
        assert_eq!(l_mult(1, 1), 2);
        assert_eq!(l_mult(-1, 1), -2);
    }

    #[test]
    fn round_is_high_word_after_bias() {
        assert_eq!(round(0x0000_8000), 1);
        assert_eq!(round(0x0000_7fff), 0);
        assert_eq!(round(-0x0000_8000), 0);
        assert_eq!(round(-0x0000_8001), -1);
        // Saturating bias: rounding the max Word32 cannot wrap.
        assert_eq!(round(i32::MAX), 32767);
    }

    #[test]
    fn mac_msu_chain_saturates() {
        assert_eq!(l_mac(0, 16384, 16384), 0x2000_0000);
        assert_eq!(l_mac(i32::MAX, 1, 1), i32::MAX);
        assert_eq!(l_msu(i32::MIN, 1, 1), i32::MIN);
        assert_eq!(l_msu(10, 1, 2), 6);
    }

    #[test]
    fn shifts_left_saturate_right_propagate_sign() {
        assert_eq!(shl(16384, 1), 32767);
        assert_eq!(shl(-16385, 1), -32768);
        assert_eq!(shl(1, 3), 8);
        assert_eq!(shl(1, -1), 0);
        assert_eq!(shr(-2, 1), -1);
        assert_eq!(shr(-1, 5), -1);
        assert_eq!(shr(3, 1), 1);
        assert_eq!(shr(3, -1), 6);
        assert_eq!(l_shl(0x4000_0000, 1), i32::MAX);
        assert_eq!(l_shl(-0x4000_0001, 1), i32::MIN);
        assert_eq!(l_shl(3, 2), 12);
        assert_eq!(l_shr(-4, 1), -2);
        assert_eq!(l_shr(-1, 40), -1);
        assert_eq!(l_shr(5, -1), 10);
    }

    #[test]
    fn rounding_shifts() {
        assert_eq!(shr_r(3, 1), 2);
        assert_eq!(shr_r(1, 1), 1);
        assert_eq!(shr_r(-3, 1), -1);
        assert_eq!(shr_r(5, 0), 5);
        assert_eq!(l_shr_r(3, 1), 2);
        assert_eq!(l_shr_r(-3, 1), -1);
        assert_eq!(l_shr_r(7, 2), 2);
    }

    #[test]
    fn mult_r_rounds_half_up() {
        assert_eq!(mult_r(16384, 16384), 8192);
        // 1·1 = 1 → (1 + 16384) >> 15 = 0 (below half).
        assert_eq!(mult_r(1, 1), 0);
        // 1·16384 = 16384 → (16384 + 16384) >> 15 = 1 (at half).
        assert_eq!(mult_r(1, 16384), 1);
    }

    #[test]
    fn deposit_and_extract_are_inverses() {
        assert_eq!(extract_h(l_deposit_h(-1234)), -1234);
        assert_eq!(extract_l(l_deposit_l(-1234)), -1234);
        assert_eq!(extract_l(0x7fff_8000u32 as i32), -32768);
        assert_eq!(extract_h(0x7fff_8000u32 as i32), 32767);
    }

    #[test]
    fn norms_count_left_shifts() {
        assert_eq!(norm_s(0), 0);
        assert_eq!(norm_s(-1), 15);
        assert_eq!(norm_s(0x4000), 0);
        assert_eq!(norm_s(1), 14);
        assert_eq!(norm_s(-32768), 0);
        assert_eq!(norm_l(0), 0);
        assert_eq!(norm_l(-1), 31);
        assert_eq!(norm_l(0x4000_0000), 0);
        assert_eq!(norm_l(1), 30);
        assert_eq!(norm_l(i32::MIN), 0);
        // The defining property: shifting by the norm count lands the
        // value in the normalised band.
        for &v in &[1i32, 2, 3, 100, 12345, 0x1234_5678, -7, -100, -65536] {
            let n = norm_l(v);
            let s = l_shl(v, n);
            assert!(
                (0x4000_0000..=0x7fff_ffff).contains(&s) || (i32::MIN..-0x4000_0000).contains(&s),
                "norm_l({v}) = {n} lands at {s:#x}",
            );
        }
    }

    #[test]
    fn div_s_fractional_division() {
        assert_eq!(div_s(1, 2), 16384);
        assert_eq!(div_s(1, 4), 8192);
        assert_eq!(div_s(3, 4), 24576);
        assert_eq!(div_s(5, 5), 32767);
        assert_eq!(div_s(0, 9), 0);
        // Truncation: 1/3 on the Q15 grid.
        assert_eq!(div_s(1, 3), 10922);
    }

    #[test]
    fn l_abs_l_negate_corners() {
        assert_eq!(l_abs(i32::MIN), i32::MAX);
        assert_eq!(l_abs(-5), 5);
        assert_eq!(l_negate(i32::MIN), i32::MAX);
        assert_eq!(l_negate(9), -9);
    }
}
