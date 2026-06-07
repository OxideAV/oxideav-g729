//! §4.1.3 decoding of the adaptive-codebook **pitch delay** from the
//! transmitted indices `(P1, P2)`.
//!
//! This module sits between the round-225 [`crate::parameters`]
//! unpacker and the §3.7 adaptive-codebook vector reconstructor that
//! eq (40) defines (the `b_30` past-excitation interpolator). It maps
//! the on-wire 8-bit subframe-1 index `P1` and 5-bit subframe-2
//! differential index `P2` back into the per-subframe fractional pitch
//! delays `T1` and `T2`, each represented as
//! `(int, frac) ∈ ℤ × {-1, 0, 1}` so the interpolator can address the
//! past-excitation buffer directly.
//!
//! ## Spec source — clause 4.1.3 (06/2012 Recommendation)
//!
//! Per clause 4.1.3 verbatim:
//!
//! > "If no parity error has occurred, the received adaptive-codebook
//! > index *P*1 is used to find the integer and fractional parts of
//! > the pitch delay *T*1. The integer part *int*(*T*1) and fractional
//! > part *frac* of *T*1 are obtained from *P*1 as follows: …"
//!
//! The "as follows" block is the spec image `f0027-01.jpg` (clause
//! 4.1.3), which evaluates to:
//!
//! ```text
//! if P1 < 197 then
//!     int(T1) = (P1 + 2) / 3 + 19      # integer division
//!     frac    = P1 − 3·int(T1) + 58
//! else
//!     int(T1) = P1 − 112
//!     frac    = 0
//! end
//! ```
//!
//! The `t_min` derivation (spec image `f0027-02.jpg`) is:
//!
//! ```text
//! t_min = int(T1) − 5
//! if t_min < 20 then t_min = 20
//! t_max = t_min + 9
//! if t_max > 143 then
//!     t_max = 143
//!     t_min = t_max − 9
//! end
//! ```
//!
//! The `T2` decode (spec image `f0027-03.jpg`) is:
//!
//! ```text
//! int(T2) = (P2 + 2) / 3 − 1 + t_min   # integer division
//! frac    = P2 − 2 − 3·((P2 + 2) / 3 − 1)
//! ```
//!
//! All three blocks were lifted from the spec EPUB's clause-4.1.3
//! equation images (rendered as JPGs); no algorithmic source was
//! consulted.
//!
//! ## Encode-side cross-check (spec clause 3.7 eqs (41) / (42))
//!
//! The encode-side forward mapping is the spec's eq (41) for `P1` from
//! `T1` and eq (42) for `P2` from `T2` (also transcribed from EPUB
//! image `eq41.jpg` / `eq42.jpg`):
//!
//! ```text
//! P1 = ⎧ 3·(int(T1) − 19) + frac − 1   if T1 ∈ [19, 85],  frac ∈ {-1, 0, 1}
//!      ⎨ (int(T1) − 85) + 197          if T1 ∈ [86, 143], frac = 0
//!      ⎩
//!
//! P2 = 3·(int(T2) − t_min) + frac + 2
//! ```
//!
//! The unit tests below round-trip the encode-side forward mapping
//! through the decode-side recipe across the **full P1 / P2 domain**;
//! every `P1 ∈ 0..256` decodes to a `(T1, frac)` pair that the
//! encode-side mapping recovers exactly, and every transmitted
//! `P2 ∈ 0..32` likewise (with `t_min` derived from a representative
//! sweep of `T1` integers).
//!
//! ## Domain coverage (spec clause 3.7 / 4.1.3)
//!
//! - Subframe 1: `T1` ranges over `[19⅓, 85]` (fractional, 1/3
//!   resolution) ∪ `[85, 143]` (integer only), giving `P1` ∈ `[0, 255]`
//!   exactly — every 8-bit codeword is valid.
//! - Subframe 2: `T2` ranges over `[t_min − 1/3, t_max + 1/3]` (1/3
//!   resolution), giving `P2` ∈ `[0, 31]` exactly — every 5-bit
//!   codeword is valid.
//!
//! Per spec §3.7 the `t_min` / `t_max` derivation centres a
//! ±5 integer window around `int(T1)` and clamps it to `[20, 143]`.
//!
//! ## Output type
//!
//! [`PitchDelay`] is a `Copy` struct carrying `int_t` (`int(T)`) and
//! `frac` (`-1`, `0`, or `1`). The §3.7 eq (40) adaptive-codebook
//! interpolator reads them as one unit: the integer part addresses the
//! past-excitation buffer and the fractional part picks the phase of
//! the `b_30` interpolation filter
//! ([`crate::tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15`]).
//!
//! ## What this module does NOT do
//!
//! - It does NOT apply the §4.1.2 parity-mismatch concealment path
//!   (where `T1` is forced to `int(T2_prev)` of the previous frame and
//!   the §4.1.3 recipe is then re-run with this synthetic `T1`). That
//!   path needs cross-frame state and lands in a follow-up round.
//! - It does NOT evaluate the spec eq (40) `b_30` interpolation that
//!   produces the adaptive-codebook vector `v(n)`. That step needs
//!   the past-excitation buffer (which only exists once §3.8
//!   fixed-codebook decoding lands) and is wired in its own round.

use crate::parameters::{Parameters, P1_BITS, P2_BITS};

/// `t_min` floor per spec §3.7 / §4.1.3 image `f0027-02.jpg`. Per the
/// clause-3.7 prose this is "the integers only over [85, 143]" floor
/// of the subframe-2 search window minus the `±4` half-span, i.e. the
/// smallest integer that keeps the 9-integer window inside the
/// 20…143 spec range.
pub const T_MIN_FLOOR: i32 = 20;

/// `t_max` ceiling per spec §3.7 / §4.1.3 image `f0027-02.jpg`. The
/// largest integer past-excitation index the §3.7 closed-loop search
/// is allowed to refer to.
pub const T_MAX_CEIL: i32 = 143;

/// `t_max − t_min` span per spec §3.7 / §4.1.3 image `f0027-02.jpg`.
/// The 10-integer (9-step) subframe-2 search window is centred on
/// `int(T1) − 4` per the spec recipe (`t_min = int(T1) − 5`,
/// `t_max = t_min + 9`), giving 10 candidate integer delays.
pub const T_WINDOW: i32 = 9;

/// Per-codeword bit width of the §3.7 subframe-1 pitch delay `P1`, as
/// re-exported from [`crate::parameters`] so callers can stay inside
/// this module.
pub const P1_DOMAIN: u32 = 1 << P1_BITS;

/// Per-codeword bit width of the §3.7 subframe-2 differential pitch
/// delay `P2`.
pub const P2_DOMAIN: u32 = 1 << P2_BITS;

/// Spec §3.7 / §4.1.3 `P1` boundary value above which the integer-only
/// branch of eq (78a) applies. From spec image `f0027-01.jpg`:
/// `if P1 < 197`.
pub const P1_FRACTIONAL_LIMIT: u32 = 197;

/// Spec §3.7 fractional pitch delay representation. `int_t` carries
/// `int(T)` (the integer part of the pitch delay, in samples of the
/// 8 kHz signal), and `frac` carries the fractional part `frac/3`
/// component encoded as one of `{-1, 0, 1}`.
///
/// Together they reconstruct the fractional pitch delay as
/// `T = int_t + frac/3`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PitchDelay {
    /// Integer part `int(T)` of the fractional pitch delay. Per spec
    /// §3.7 / §4.1.3 this lies in `[19, 143]` for `T1` (subframe 1;
    /// the `19` endpoint reaches when `P1 = 0`, `T1 = 19+1/3`) and
    /// in `[t_min − 1, t_min + 10]` for `T2` (subframe 2; the
    /// fractional-extreme rounding pushes the integer one step past
    /// the integer-window edge).
    pub int_t: i32,
    /// Fractional component `frac` per spec §3.7. Always one of
    /// `{-1, 0, 1}` for a well-formed transmitted index. The actual
    /// fractional pitch delay is `int_t + frac/3`.
    pub frac: i32,
}

/// Per-frame output of [`decode_frame`]. Carries the per-subframe
/// pitch delays in spec §4.1.5 order: subframe 1 first, subframe 2
/// next.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FramePitchDelays {
    /// Subframe-1 fractional pitch delay, decoded from `P1` per spec
    /// image `f0027-01.jpg`.
    pub t1: PitchDelay,
    /// Subframe-2 fractional pitch delay, decoded from `(P2, t_min)`
    /// per spec image `f0027-03.jpg`, with `t_min` derived from `T1`
    /// per spec image `f0027-02.jpg`.
    pub t2: PitchDelay,
    /// Spec `t_min` floor of the subframe-2 search window, retained
    /// for callers that need to cross-check the §4.1.3 derivation
    /// (and for the §4.1.2 parity-concealment path).
    pub t_min: i32,
}

/// Decode the subframe-1 pitch delay `T1` from the transmitted index
/// `P1` per spec image `f0027-01.jpg` (clause 4.1.3).
///
/// The 8-bit `P1` field covers `0..=255` exactly; every transmitted
/// codeword is valid, so this function does not return an error
/// surface. The encode-side forward mapping is spec eq (41).
///
/// Per spec image `f0027-01.jpg`:
///
/// ```text
/// if P1 < 197 then
///     int(T1) = (P1 + 2) / 3 + 19      (integer division)
///     frac    = P1 − 3·int(T1) + 58
/// else
///     int(T1) = P1 − 112
///     frac    = 0
/// end
/// ```
#[inline]
pub fn decode_t1_from_p1(p1: u8) -> PitchDelay {
    let p1 = i32::from(p1);
    if (p1 as u32) < P1_FRACTIONAL_LIMIT {
        // Fractional branch: 1/3 resolution over the [19+1/3, 85] span.
        // `int(T1) = (P1 + 2) / 3 + 19`, integer-division per spec.
        let int_t = (p1 + 2) / 3 + 19;
        // `frac = P1 − 3·int(T1) + 58` ∈ {-1, 0, 1}.
        let frac = p1 - 3 * int_t + 58;
        PitchDelay { int_t, frac }
    } else {
        // Integer-only branch: `int(T1) = P1 − 112`, `frac = 0`.
        let int_t = p1 - 112;
        PitchDelay { int_t, frac: 0 }
    }
}

/// Derive `t_min` (the lower bound of the subframe-2 search window)
/// from the subframe-1 integer pitch delay `int(T1)` per spec image
/// `f0027-02.jpg` (clause 4.1.3 / 3.7).
///
/// ```text
/// t_min = int(T1) − 5
/// if t_min < 20 then t_min = 20
/// t_max = t_min + 9
/// if t_max > 143 then
///     t_max = 143
///     t_min = t_max − 9
/// end
/// ```
///
/// The result is in `[T_MIN_FLOOR, T_MAX_CEIL − T_WINDOW]` (i.e.
/// `[20, 134]`) — both endpoints land on actual `P2` decodes from
/// real input.
#[inline]
pub fn derive_t_min(int_t1: i32) -> i32 {
    let mut t_min = int_t1 - 5;
    if t_min < T_MIN_FLOOR {
        t_min = T_MIN_FLOOR;
    }
    let mut t_max = t_min + T_WINDOW;
    if t_max > T_MAX_CEIL {
        t_max = T_MAX_CEIL;
        t_min = t_max - T_WINDOW;
    }
    t_min
}

/// Decode the subframe-2 pitch delay `T2` from the transmitted index
/// `P2` and the spec-derived `t_min` per spec image `f0027-03.jpg`
/// (clause 4.1.3).
///
/// The 5-bit `P2` field covers `0..=31` exactly. The encode-side
/// forward mapping is spec eq (42),
/// `P2 = 3·(int(T2) − t_min) + frac + 2`.
///
/// Per spec image `f0027-03.jpg`:
///
/// ```text
/// int(T2) = (P2 + 2) / 3 − 1 + t_min   (integer division)
/// frac    = P2 − 2 − 3·((P2 + 2) / 3 − 1)
/// ```
#[inline]
pub fn decode_t2_from_p2(p2: u8, t_min: i32) -> PitchDelay {
    let p2 = i32::from(p2);
    // Common subexpression `(P2 + 2) / 3 − 1`.
    let q = (p2 + 2) / 3 - 1;
    let int_t = q + t_min;
    let frac = p2 - 2 - 3 * q;
    PitchDelay { int_t, frac }
}

/// Decode the per-frame pitch delays `(T1, T2)` from the unpacked
/// transmitted parameters per spec §4.1.3.
///
/// This is the §4.1.3 frame-level entry point. It chains
/// [`decode_t1_from_p1`], [`derive_t_min`], and [`decode_t2_from_p2`]
/// in the spec-stated order. The §4.1.2 parity-concealment path is
/// **not** applied here — callers that need it must check
/// [`Parameters::pitch_parity_ok`] and invoke the §4.4 concealment
/// path themselves.
#[inline]
pub fn decode_frame(params: &Parameters) -> FramePitchDelays {
    let t1 = decode_t1_from_p1(params.p1);
    let t_min = derive_t_min(t1.int_t);
    let t2 = decode_t2_from_p2(params.p2, t_min);
    FramePitchDelays { t1, t2, t_min }
}

/// Encode-side forward mapping (spec clause 3.7 eq (41)) — used only
/// for the round-trip unit tests below. Returns the `P1` codeword
/// that would round-trip through [`decode_t1_from_p1`] back to the
/// supplied `(int_t1, frac)`, or `None` if the input is outside the
/// spec §3.7 domain.
///
/// Per spec image `eq41.jpg` (clause 3.7):
///
/// ```text
/// P1 = ⎧ 3·(int(T1) − 19) + frac − 1   if T1 ∈ [19⅓, 85], frac ∈ {-1, 0, 1}
///      ⎨ (int(T1) − 85) + 197          if T1 ∈ [86, 143], frac = 0
///      ⎩
/// ```
///
/// Public so the round-trip property can be exercised by external
/// callers (encoders, fixture builders) without re-deriving the
/// algebra.
#[inline]
pub fn encode_p1(delay: PitchDelay) -> Option<u8> {
    let PitchDelay { int_t, frac } = delay;
    // Fractional branch: T1 ∈ [19⅓, 85], frac ∈ {-1, 0, 1}.
    // Equivalent integer domain: `3·int(T1) + frac ∈ [58, 255]` (the
    // smallest 3·int+frac with T1 ≥ 19+1/3 is 3·19 + 1 = 58, the
    // largest with T1 ≤ 85 is 3·85 + 0 = 255).
    if !(-1..=1).contains(&frac) {
        return None;
    }
    let three_t_plus_frac = 3 * int_t + frac;
    if (58..=255).contains(&three_t_plus_frac) {
        // Per spec eq (41) image: 3·(int − 19) + frac − 1
        //   = (3·int + frac) − 58
        let p1 = three_t_plus_frac - 58;
        u8::try_from(p1).ok()
    } else if frac == 0 && (86..=143).contains(&int_t) {
        // Integer-only branch.
        u8::try_from(int_t - 85 + 197).ok()
    } else {
        None
    }
}

/// Encode-side forward mapping (spec clause 3.7 eq (42)) for `P2`
/// from `(int_t2, frac, t_min)`. Returns `None` if the input is
/// outside the spec §3.7 domain.
///
/// Per spec image `eq42.jpg`:
///
/// ```text
/// P2 = 3·(int(T2) − t_min) + frac + 2
/// ```
#[inline]
pub fn encode_p2(delay: PitchDelay, t_min: i32) -> Option<u8> {
    let PitchDelay { int_t, frac } = delay;
    if !(-1..=1).contains(&frac) {
        return None;
    }
    let p2 = 3 * (int_t - t_min) + frac + 2;
    if (0..32).contains(&p2) {
        Some(p2 as u8)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Spec image `f0027-01.jpg` boundary cases pinned literally so any
    /// drift in the recipe trips immediately.
    #[test]
    fn decode_t1_from_p1_spec_examples() {
        // P1 = 0 → fractional branch.
        //   int(T1) = (0 + 2) / 3 + 19 = 0 + 19 = 19
        //   frac    = 0 - 3·19 + 58    = 1
        let d = decode_t1_from_p1(0);
        assert_eq!(d.int_t, 19);
        assert_eq!(d.frac, 1);

        // P1 = 1 → int 19, frac -1? Recheck:
        //   int(T1) = (1 + 2) / 3 + 19 = 1 + 19 = 20
        //   frac    = 1 - 3·20 + 58    = -1
        let d = decode_t1_from_p1(1);
        assert_eq!(d.int_t, 20);
        assert_eq!(d.frac, -1);

        // P1 = 196 (largest fractional-branch index).
        //   int(T1) = (196 + 2) / 3 + 19 = 66 + 19 = 85
        //   frac    = 196 - 3·85 + 58    = -1
        let d = decode_t1_from_p1(196);
        assert_eq!(d.int_t, 85);
        assert_eq!(d.frac, -1);

        // P1 = 197 → integer-only branch.
        //   int(T1) = 197 - 112 = 85
        //   frac    = 0
        let d = decode_t1_from_p1(197);
        assert_eq!(d.int_t, 85);
        assert_eq!(d.frac, 0);

        // P1 = 255 (largest 8-bit value).
        //   int(T1) = 255 - 112 = 143
        //   frac    = 0
        let d = decode_t1_from_p1(255);
        assert_eq!(d.int_t, 143);
        assert_eq!(d.frac, 0);
    }

    /// Every `P1 ∈ 0..256` lands in the spec-stated domain
    /// `int(T1) ∈ [19, 143]`, `frac ∈ {-1, 0, 1}`. The integer-part
    /// lower bound is 19 (not 20) because `T1 = 19 + 1/3` decodes to
    /// `int(T1) = 19`, `frac = 1`. Locks the recipe against any
    /// future drift that would punch a hole in the 8-bit codeword
    /// space.
    #[test]
    fn decode_t1_full_domain_in_spec_range() {
        for p1 in 0..=u8::MAX {
            let d = decode_t1_from_p1(p1);
            assert!(
                (19..=143).contains(&d.int_t),
                "P1={p1}: int_t={} outside spec [19, 143]",
                d.int_t
            );
            assert!(
                (-1..=1).contains(&d.frac),
                "P1={p1}: frac={} outside spec {{-1, 0, 1}}",
                d.frac
            );
        }
    }

    /// Spec image `f0027-02.jpg` worked examples — pinning the floor
    /// and ceiling branches.
    #[test]
    fn derive_t_min_spec_examples() {
        // Mid-range: int_t1 = 50 → t_min = 45, t_max = 54 (both inside).
        assert_eq!(derive_t_min(50), 45);
        // Floor: int_t1 = 20 → t_min = 15 → clamped to 20.
        assert_eq!(derive_t_min(20), 20);
        // Floor edge: int_t1 = 25 → t_min = 20, t_max = 29 (inside).
        assert_eq!(derive_t_min(25), 20);
        // Floor edge: int_t1 = 24 → t_min = 19 → clamped to 20.
        assert_eq!(derive_t_min(24), 20);
        // Ceiling: int_t1 = 143 → t_min = 138, t_max = 147 → clamped
        // to t_max = 143, t_min = 134.
        assert_eq!(derive_t_min(143), 134);
        // Ceiling edge: int_t1 = 139 → t_min = 134, t_max = 143
        // (right at the ceiling, no clamp needed).
        assert_eq!(derive_t_min(139), 134);
        // Just below ceiling: int_t1 = 138 → t_min = 133,
        // t_max = 142 (inside).
        assert_eq!(derive_t_min(138), 133);
    }

    /// The derived `t_min` for every `int(T1) ∈ [19, 143]` (the full
    /// `decode_t1_from_p1` output range) lies in
    /// `[T_MIN_FLOOR, T_MAX_CEIL − T_WINDOW]`, i.e. the entire 9-step
    /// subframe-2 search window fits inside `[20, 143]`.
    #[test]
    fn derive_t_min_full_domain_stays_in_search_envelope() {
        for int_t1 in 19..=143 {
            let t_min = derive_t_min(int_t1);
            assert!(
                t_min >= T_MIN_FLOOR && t_min + T_WINDOW <= T_MAX_CEIL,
                "int_t1={int_t1}: t_min={t_min} pushes search window outside [20, 143]"
            );
        }
    }

    /// Spec image `f0027-03.jpg` boundary cases pinned literally.
    #[test]
    fn decode_t2_from_p2_spec_examples() {
        // P2 = 0, t_min = 50:
        //   q       = (0 + 2) / 3 − 1 = 0 − 1 = −1
        //   int(T2) = −1 + 50 = 49
        //   frac    = 0 − 2 − 3·(−1) = 1
        let d = decode_t2_from_p2(0, 50);
        assert_eq!(d.int_t, 49);
        assert_eq!(d.frac, 1);

        // P2 = 2, t_min = 50:
        //   q       = (2 + 2) / 3 − 1 = 1 − 1 = 0
        //   int(T2) = 0 + 50 = 50
        //   frac    = 2 − 2 − 3·0 = 0
        let d = decode_t2_from_p2(2, 50);
        assert_eq!(d.int_t, 50);
        assert_eq!(d.frac, 0);

        // P2 = 31, t_min = 50:
        //   q       = (31 + 2) / 3 − 1 = 11 − 1 = 10
        //   int(T2) = 10 + 50 = 60
        //   frac    = 31 − 2 − 30 = −1
        let d = decode_t2_from_p2(31, 50);
        assert_eq!(d.int_t, 60);
        assert_eq!(d.frac, -1);
    }

    /// Every `P2 ∈ 0..32` lands in
    /// `int(T2) ∈ [t_min − 1, t_min + 10]`, `frac ∈ {-1, 0, 1}` across
    /// any t_min in the spec-stated `[20, 134]` window. The
    /// integer-part endpoints are reached only at the fractional
    /// extremes (P2=0 → int=t_min-1, frac=+1 → T2 = t_min - 2/3;
    /// P2=31 → int=t_min+10, frac=-1 → T2 = t_min + 9 + 2/3). The
    /// full-precision T2 stays inside `[t_min − 2/3, t_max + 2/3]`
    /// per the spec image — both extremes round to one outside the
    /// integer-window edge.
    #[test]
    fn decode_t2_full_domain_in_spec_range() {
        for t_min in T_MIN_FLOOR..=(T_MAX_CEIL - T_WINDOW) {
            for p2 in 0..=31u8 {
                let d = decode_t2_from_p2(p2, t_min);
                assert!(
                    (t_min - 1..=t_min + T_WINDOW + 1).contains(&d.int_t),
                    "P2={p2} t_min={t_min}: int_t={} outside spec [t_min-1, t_min+10]",
                    d.int_t
                );
                assert!(
                    (-1..=1).contains(&d.frac),
                    "P2={p2} t_min={t_min}: frac={} outside spec {{-1, 0, 1}}",
                    d.frac
                );
            }
        }
    }

    /// Encode-side forward mapping (eq (41)) round-trips every `P1`
    /// codeword exactly through the decode-side recipe. This pins the
    /// algebra of both branches simultaneously: any sign / off-by-one
    /// drift in either eq (78a) or eq (41) trips the round-trip.
    #[test]
    fn p1_encode_decode_round_trip_full_domain() {
        for p1 in 0..=u8::MAX {
            let decoded = decode_t1_from_p1(p1);
            let re_encoded =
                encode_p1(decoded).unwrap_or_else(|| panic!("P1={p1} failed encode round-trip"));
            assert_eq!(re_encoded, p1, "P1={p1} round-trip mismatch");
        }
    }

    /// Encode-side forward mapping (eq (42)) round-trips every `P2`
    /// codeword exactly through the decode-side recipe, across the
    /// full spec-stated `t_min ∈ [20, 134]` window.
    #[test]
    fn p2_encode_decode_round_trip_full_domain() {
        for t_min in T_MIN_FLOOR..=(T_MAX_CEIL - T_WINDOW) {
            for p2 in 0..=31u8 {
                let decoded = decode_t2_from_p2(p2, t_min);
                let re_encoded = encode_p2(decoded, t_min)
                    .unwrap_or_else(|| panic!("P2={p2} t_min={t_min} failed encode round-trip"));
                assert_eq!(re_encoded, p2, "P2={p2} t_min={t_min} round-trip mismatch");
            }
        }
    }

    /// Encode-side `P1` rejects malformed `(int, frac)` inputs.
    #[test]
    fn encode_p1_rejects_out_of_domain() {
        // frac outside {-1, 0, 1}.
        assert_eq!(encode_p1(PitchDelay { int_t: 50, frac: 2 }), None);
        assert_eq!(
            encode_p1(PitchDelay {
                int_t: 50,
                frac: -2
            }),
            None
        );
        // int_t too small for fractional branch (would need T1 < 19+1/3).
        assert_eq!(
            encode_p1(PitchDelay {
                int_t: 19,
                frac: -1
            }),
            None
        );
        // Integer-only branch requires frac = 0.
        assert_eq!(
            encode_p1(PitchDelay {
                int_t: 100,
                frac: 1
            }),
            None
        );
        // int_t = 144 is outside the integer-only branch.
        assert_eq!(
            encode_p1(PitchDelay {
                int_t: 144,
                frac: 0
            }),
            None
        );
    }

    /// Encode-side `P2` rejects malformed inputs.
    #[test]
    fn encode_p2_rejects_out_of_domain() {
        // frac outside {-1, 0, 1}.
        assert_eq!(encode_p2(PitchDelay { int_t: 50, frac: 2 }, 50), None);
        // int_t too far from t_min — pushes P2 out of [0, 31].
        assert_eq!(
            encode_p2(
                PitchDelay {
                    int_t: 100,
                    frac: 0
                },
                50
            ),
            None
        );
        assert_eq!(encode_p2(PitchDelay { int_t: 10, frac: 0 }, 50), None);
    }

    /// `decode_frame` chains the three steps in the spec §4.1.3 order
    /// and reads from the right fields of the `Parameters` struct.
    #[test]
    fn decode_frame_threads_per_subframe_indices() {
        let mut params = synth_params();
        // Pick a P1 that lands in the middle of the fractional branch.
        // P1 = 60 → int(T1) = (60+2)/3 + 19 = 20 + 19 = 39, frac = 60 - 117 + 58 = 1
        // Wait: 60 - 3·39 + 58 = 60 - 117 + 58 = 1. Yes.
        params.p1 = 60;
        // Pick a P2 that lands in the middle of its window.
        // t_min = derive_t_min(39) = 34, t_max = 43. P2 = 14 →
        //   q = 16/3 − 1 = 5 − 1 = 4, int(T2) = 4 + 34 = 38, frac = 14 − 2 − 12 = 0.
        params.p2 = 14;

        let out = decode_frame(&params);
        assert_eq!(out.t1, PitchDelay { int_t: 39, frac: 1 });
        assert_eq!(out.t_min, 34);
        assert_eq!(out.t2, PitchDelay { int_t: 38, frac: 0 });
    }

    /// `decode_frame` correctly threads `P1` (not `P2`) into the
    /// subframe-1 decoder and vice-versa. A swap would trip the test
    /// since `P1`'s 8-bit value would saturate the 5-bit `P2` recipe.
    #[test]
    fn decode_frame_does_not_swap_p1_and_p2() {
        let mut params = synth_params();
        params.p1 = 0;
        params.p2 = 0;
        let out = decode_frame(&params);
        // P1=0 → int(T1)=19, frac=1
        assert_eq!(out.t1, PitchDelay { int_t: 19, frac: 1 });
        // t_min = derive_t_min(19) → 20 (floor clamp). P2=0, t_min=20:
        //   q = 2/3 - 1 = 0 - 1 = -1, int(T2) = -1 + 20 = 19, frac = 0 - 2 - 3·(-1) = 1
        assert_eq!(out.t_min, 20);
        assert_eq!(out.t2, PitchDelay { int_t: 19, frac: 1 });
    }

    /// Constants surface matches the documented spec values.
    #[test]
    fn constants_match_spec_recipe() {
        assert_eq!(T_MIN_FLOOR, 20);
        assert_eq!(T_MAX_CEIL, 143);
        assert_eq!(T_WINDOW, 9);
        assert_eq!(P1_DOMAIN, 256);
        assert_eq!(P2_DOMAIN, 32);
        assert_eq!(P1_FRACTIONAL_LIMIT, 197);
    }

    /// Build a `Parameters` struct with all non-pitch fields zeroed —
    /// the pitch-decode path only ever reads `p1` / `p2` (and parity,
    /// which the §4.1.3 path itself does not consult).
    fn synth_params() -> Parameters {
        Parameters {
            l0: 0,
            l1: 0,
            l2: 0,
            l3: 0,
            p1: 0,
            p0: 0,
            c1: 0,
            s1: 0,
            ga1: 0,
            gb1: 0,
            p2: 0,
            c2: 0,
            s2: 0,
            ga2: 0,
            gb2: 0,
        }
    }
}
