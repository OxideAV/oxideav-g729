//! §3.8 / §4.1.4 **pitch sharpening** of the fixed-codebook codevector
//! `c(n)` (the spec eq (48) modification applied when the integer
//! pitch delay `int(T) < 40`).
//!
//! This module sits between the round-266 [`crate::fixed_codebook`]
//! codevector builder and the §3.10 / §4.1.6 per-subframe excitation
//! mixer. It applies the spec §3.8 adaptive pre-filter
//! `P(z) = 1/(1 − β·z^−T)` to the unmodified algebraic codevector,
//! which — for the decode path of clause 4.1.4 — reduces to the eq (48)
//! recursive in-place modification.
//!
//! ## Spec source — clause 3.8 / 4.1.4 (06/2012 Recommendation)
//!
//! Per clause 3.8 verbatim:
//!
//! > "A special feature incorporated in the codebook is that the
//! > selected codebook vector is filtered through an adaptive
//! > pre-filter *P(z)* that enhances harmonic components to improve the
//! > quality of the reconstructed speech. Here the filter \[eq (46)\]
//! > is used, where *T* is the integer component of the pitch delay of
//! > the current subframe, and *β* is a pitch gain. The value of *β* is
//! > made adaptive by using the quantized adaptive-codebook gain from
//! > the previous subframe \[eq (47)\] ... For delays less than 40, the
//! > codebook *c(n)* of equation (45) is modified according to \[eq
//! > (48)\]."
//!
//! Spec eq (47) (EPUB image `eq47.jpg`) defines `β`:
//!
//! ```text
//! β = ĝ_p^(m−1)     bounded by   0.2 ≤ β ≤ 0.8                  (eq 47)
//! ```
//!
//! i.e. `β` is the quantized adaptive-codebook gain `ĝ_p` of the
//! **previous** subframe, clamped to the closed interval `[0.2, 0.8]`.
//!
//! Spec eq (48) (EPUB image `eq48.jpg`) gives the modification:
//!
//! ```text
//!          ⎧ c(n)                  n = 0, …, T − 1
//! c(n) =   ⎨                                                     (eq 48)
//!          ⎩ c(n) + β·c(n − T)     n = T, …, 39
//! ```
//!
//! where `T` is the integer part of the current subframe's pitch delay
//! (spec clause 3.8: "*T* is the integer component of the pitch delay
//! of the current subframe").
//!
//! Per clause 4.1.4 verbatim, this is the decoder's final fixed-vector
//! step:
//!
//! > "Once the pulse positions and signs are decoded, the
//! > fixed-codebook vector *c(n)* is constructed using equation (45).
//! > If the integer part of the pitch delay *T* is less than the
//! > subframe size 40, *c(n)* is modified according to equation (48)."
//!
//! ## Recurrence note
//!
//! Equation (48) is the time-domain expansion of the all-pole filter
//! `1/(1 − β·z^−T)`: the `c(n − T)` term on the right-hand side refers
//! to the **already-modified** sample (the filter is recursive), so for
//! a small delay `T` the contribution propagates forward in steps of
//! `T`. Iterating `n` from `T` to `39` in increasing order — reading
//! `c(n − T)` *after* it has itself been updated — realises the
//! recursion exactly. With the round-266 codevector being four ±1
//! impulses, the first pulse at position `m_0` seeds the recursion and
//! the sharpened taps appear at `m_0 + T`, `m_0 + 2T`, … with
//! geometrically decaying weights `β`, `β²`, ….

use crate::fixed_codebook::SUBFRAME_SIZE;

/// Spec eq (47) lower clamp on the pitch gain `β`.
pub const BETA_MIN: f32 = 0.2;

/// Spec eq (47) upper clamp on the pitch gain `β`.
pub const BETA_MAX: f32 = 0.8;

/// Clamp the previous-subframe quantized adaptive-codebook gain
/// `ĝ_p^(m−1)` into the spec eq (47) pitch-gain `β` range `[0.2, 0.8]`.
///
/// Per spec eq (47): `β = ĝ_p^(m−1)` bounded by `0.2 ≤ β ≤ 0.8`. The
/// input is the quantized adaptive-codebook gain reconstructed for the
/// previous subframe (see [`crate::gain_reconstruct`]); the output is
/// the pitch gain used by the eq (48) sharpening of the current
/// subframe.
///
/// A non-finite input (`NaN`) is clamped to [`BETA_MIN`] so the
/// downstream recurrence stays finite; this never arises from a
/// well-formed gain reconstruction but keeps the sharpening total even
/// when driven from raw values.
#[inline]
pub fn clamp_beta(g_p_prev: f32) -> f32 {
    if g_p_prev.is_nan() {
        return BETA_MIN;
    }
    g_p_prev.clamp(BETA_MIN, BETA_MAX)
}

/// Apply the spec eq (48) pitch sharpening to a fixed-codebook
/// codevector, returning the modified `[f32; 40]`.
///
/// `c` is the unmodified codevector from
/// [`crate::fixed_codebook::build_codevector`] (four ±1 impulses).
/// `int_t` is the integer part of the current subframe's pitch delay
/// `T` (see [`crate::pitch_decode::PitchDelay::int_t`]). `g_p_prev` is
/// the quantized adaptive-codebook gain of the **previous** subframe;
/// it is clamped to `[0.2, 0.8]` per eq (47) internally.
///
/// Per clause 4.1.4 the modification is applied **only** when
/// `int(T) < 40` (= [`SUBFRAME_SIZE`]). For `int(T) ≥ 40` the
/// codevector is returned unchanged (promoted to `f32`), matching the
/// decoder branch "if the integer part of the pitch delay *T* is less
/// than the subframe size 40, *c(n)* is modified".
///
/// A non-positive `int_t` (which cannot occur for a well-formed
/// [`crate::pitch_decode`] output — `int(T) ≥ 19`) leaves the vector
/// unmodified: the recurrence `c(n) + β·c(n − T)` is only defined for a
/// strictly positive delay.
#[inline]
pub fn sharpen(c: &[i8; SUBFRAME_SIZE], int_t: i32, g_p_prev: f32) -> [f32; SUBFRAME_SIZE] {
    let mut out = [0.0f32; SUBFRAME_SIZE];
    for (o, &ci) in out.iter_mut().zip(c.iter()) {
        *o = f32::from(ci);
    }
    // Per clause 4.1.4: only modify when int(T) < subframe size 40.
    // A non-positive delay has no defined recurrence; leave unchanged.
    if int_t <= 0 || int_t >= SUBFRAME_SIZE as i32 {
        return out;
    }
    let beta = clamp_beta(g_p_prev);
    let t = int_t as usize;
    // Per eq (48): for n = T … 39, c(n) += β·c(n − T), reading the
    // already-modified c(n − T) so the recursion of P(z) = 1/(1 − β·z^−T)
    // is realised by the forward in-place sweep.
    for n in t..SUBFRAME_SIZE {
        out[n] += beta * out[n - t];
    }
    out
}

/// Energy of a sharpened codevector — `Σ c(n)²` over the subframe.
///
/// Exposed so the §3.9.1 gain-prediction energy term (spec eq (66))
/// and the unit tests can read the post-sharpening codevector energy
/// without re-deriving the sum. The unmodified four-pulse codevector
/// has energy exactly `4.0`; sharpening with `int(T) < 40` adds the
/// forward-propagated taps, so the energy grows monotonically with
/// `β` (for a non-trivial overlap).
#[inline]
pub fn codevector_energy(c: &[f32; SUBFRAME_SIZE]) -> f32 {
    c.iter().map(|x| x * x).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_codebook::{build_codevector, decode_pulses};

    /// Spec eq (47) clamp: below `0.2` → `0.2`, above `0.8` → `0.8`,
    /// in-range passes through.
    #[test]
    fn clamp_beta_eq_47() {
        assert_eq!(clamp_beta(0.0), BETA_MIN);
        assert_eq!(clamp_beta(0.1), BETA_MIN);
        assert_eq!(clamp_beta(0.2), 0.2);
        assert_eq!(clamp_beta(0.5), 0.5);
        assert_eq!(clamp_beta(0.8), 0.8);
        assert_eq!(clamp_beta(0.9), BETA_MAX);
        assert_eq!(clamp_beta(1.0), BETA_MAX);
        assert_eq!(clamp_beta(2.5), BETA_MAX);
        assert_eq!(clamp_beta(-1.0), BETA_MIN);
        assert_eq!(clamp_beta(f32::NAN), BETA_MIN);
    }

    /// `int(T) ≥ 40` leaves the codevector unmodified (only promoted
    /// to `f32`), per the clause-4.1.4 "less than the subframe size 40"
    /// guard.
    #[test]
    fn no_sharpen_when_delay_ge_subframe() {
        let pulses = decode_pulses(0b0001, 0b1111).unwrap();
        let c = build_codevector(&pulses);
        let out = sharpen(&c, 40, 0.5);
        for (o, &ci) in out.iter().zip(c.iter()) {
            assert_eq!(*o, f32::from(ci));
        }
        // Energy unchanged: still exactly four ±1 pulses.
        assert_eq!(codevector_energy(&out), 4.0);
    }

    /// A non-positive delay is a no-op (the recurrence is undefined).
    #[test]
    fn no_sharpen_for_non_positive_delay() {
        let pulses = decode_pulses(0, 0b1111).unwrap();
        let c = build_codevector(&pulses);
        for bad in [0, -1, -40] {
            let out = sharpen(&c, bad, 0.5);
            for (o, &ci) in out.iter().zip(c.iter()) {
                assert_eq!(*o, f32::from(ci));
            }
        }
    }

    /// Worked eq (48) example with a hand-built impulse train.
    ///
    /// Place a single +1 impulse at position 0 (C = 0 → m_0 = 0 with
    /// s_0 = +1) and delay `T = 10`. The sharpened taps should appear
    /// at 0, 10, 20, 30 with weights `1, β, β², β³`.
    #[test]
    fn sharpen_eq_48_geometric_train() {
        // C = 0 → positions [0, 1, 2, 3]; isolate the pulse at 0 by
        // building a vector with only that pulse (signs chosen so the
        // others land elsewhere, then we reason on position 0's tail).
        // Simpler: hand-build the codevector directly.
        let mut c = [0i8; SUBFRAME_SIZE];
        c[0] = 1;
        let beta = 0.5f32;
        let out = sharpen(&c, 10, beta);
        // n = 0 untouched (n < T).
        assert!((out[0] - 1.0).abs() < 1e-6);
        // n = 10: c(10) += β·c(0) = 0 + 0.5·1 = 0.5.
        assert!((out[10] - 0.5).abs() < 1e-6);
        // n = 20: c(20) += β·c(10) = 0 + 0.5·0.5 = 0.25.
        assert!((out[20] - 0.25).abs() < 1e-6);
        // n = 30: c(30) += β·c(20) = 0 + 0.5·0.25 = 0.125.
        assert!((out[30] - 0.125).abs() < 1e-6);
        // Everything else zero.
        for (n, v) in out.iter().enumerate() {
            if ![0, 10, 20, 30].contains(&n) {
                assert!(v.abs() < 1e-6, "non-zero at n={n}: {v}");
            }
        }
    }

    /// The recurrence reads the already-modified `c(n − T)` so a second
    /// pulse inside the tail interacts. With a +1 pulse at 0 and a +1
    /// pulse at 5 and T = 5, the taps accumulate:
    /// c(5) = 1 + β·1, c(10) = β·(1+β), …
    #[test]
    fn sharpen_eq_48_overlapping_pulses_recursive() {
        let mut c = [0i8; SUBFRAME_SIZE];
        c[0] = 1;
        c[5] = 1;
        let beta = 0.5f32;
        let out = sharpen(&c, 5, beta);
        // n = 0 untouched.
        assert!((out[0] - 1.0).abs() < 1e-6);
        // n = 5: original 1 + β·c(0) = 1 + 0.5 = 1.5.
        assert!((out[5] - 1.5).abs() < 1e-6);
        // n = 10: 0 + β·c(5) = 0.5·1.5 = 0.75 (reads modified c(5)).
        assert!((out[10] - 0.75).abs() < 1e-6);
        // n = 15: 0 + β·c(10) = 0.5·0.75 = 0.375.
        assert!((out[15] - 0.375).abs() < 1e-6);
    }

    /// Sharpening increases codevector energy for an overlapping delay
    /// (the forward taps are added in phase with the seed pulses here),
    /// confirming the modification is wired in the right direction.
    #[test]
    fn sharpen_increases_energy_for_short_delay() {
        let mut c = [0i8; SUBFRAME_SIZE];
        c[0] = 1;
        let unmod = {
            let mut f = [0.0f32; SUBFRAME_SIZE];
            for (o, &ci) in f.iter_mut().zip(c.iter()) {
                *o = f32::from(ci);
            }
            codevector_energy(&f)
        };
        let out = sharpen(&c, 10, 0.8);
        let sharpened = codevector_energy(&out);
        assert!(
            sharpened > unmod,
            "sharpened energy {sharpened} should exceed unmodified {unmod}",
        );
    }

    /// β = 0 (after clamping a tiny gain to the 0.2 floor still adds a
    /// tap, but a literal β = 0.2 yields the expected first tap) — and
    /// the eq (47) clamp is applied inside `sharpen`, so a previous
    /// gain of 5.0 behaves identically to 0.8.
    #[test]
    fn sharpen_applies_eq_47_clamp_internally() {
        let mut c = [0i8; SUBFRAME_SIZE];
        c[0] = 1;
        let high = sharpen(&c, 10, 5.0);
        let at_max = sharpen(&c, 10, 0.8);
        for (a, b) in high.iter().zip(at_max.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
        let low = sharpen(&c, 10, -3.0);
        let at_min = sharpen(&c, 10, 0.2);
        for (a, b) in low.iter().zip(at_min.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }

    /// Realistic codevector from the decode path: decode pulses from
    /// `(C, S)`, build the codevector, sharpen with a representative
    /// short delay. Confirms the four ±1 seeds are preserved at
    /// `n < T` and the tail is non-trivial.
    #[test]
    fn sharpen_real_codevector_preserves_head() {
        // C = 0b0001 → positions [5, 1, 2, 3]; S = 0b1111 → all +1.
        let pulses = decode_pulses(0b0001, 0b1111).unwrap();
        let c = build_codevector(&pulses);
        let t = 19; // representative minimal int(T1).
        let out = sharpen(&c, t, 0.6);
        // Head (n < 19) equals the promoted original pulses.
        for n in 0..(t as usize) {
            assert!((out[n] - f32::from(c[n])).abs() < 1e-6, "head drift n={n}");
        }
        // Tail picks up β·c(n − T): positions 1,2,3,5 are < 19 so their
        // sharpened echoes land at 1+19=20, 2+19=21, 3+19=22, 5+19=24.
        assert!((out[20] - 0.6 * f32::from(c[1])).abs() < 1e-6);
        assert!((out[21] - 0.6 * f32::from(c[2])).abs() < 1e-6);
        assert!((out[22] - 0.6 * f32::from(c[3])).abs() < 1e-6);
        assert!((out[24] - 0.6 * f32::from(c[5])).abs() < 1e-6);
    }

    /// The unmodified four-pulse codevector has energy exactly 4.0
    /// (sanity tie to the round-266 invariant).
    #[test]
    fn unmodified_codevector_energy_is_four() {
        let pulses = decode_pulses(512, 0b1010).unwrap();
        let c = build_codevector(&pulses);
        let promoted = sharpen(&c, 40, 0.5); // ≥40 → unmodified
        assert_eq!(codevector_energy(&promoted), 4.0);
    }
}
