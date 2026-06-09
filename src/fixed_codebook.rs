//! §3.8 / §4.1.4 decoding of the **fixed (algebraic) codebook** from
//! the transmitted indices `(C, S)`.
//!
//! This module sits between the round-225 [`crate::parameters`]
//! unpacker and the per-subframe excitation builder. It maps the
//! on-wire 13-bit `C` codeword into the 4 pulse positions
//! `(m_0, m_1, m_2, m_3)` per spec Table 7 / eq (62), the on-wire
//! 4-bit `S` codeword into the 4 pulse signs `(s_0, s_1, s_2, s_3)`
//! per spec eq (61), and constructs the 40-sample codevector `c(n)`
//! per spec eq (45).
//!
//! ## Spec source — clause 3.8 / 3.8.2 / 4.1.4 (06/2012 Recommendation)
//!
//! Per clause 3.8 verbatim:
//!
//! > "The fixed codebook is based on an algebraic codebook structure
//! > using an interleaved single-pulse permutation (ISPP) design. In
//! > this codebook, each codebook vector contains four non-zero
//! > pulses. Each pulse can have either the amplitudes +1 or –1, and
//! > can assume the positions given in Table 7."
//!
//! Spec Table 7 gives the per-pulse allowed positions:
//!
//! ```text
//! Pulse  Sign   Allowed positions
//! i_0    ±1     0, 5, 10, 15, 20, 25, 30, 35
//! i_1    ±1     1, 6, 11, 16, 21, 26, 31, 36
//! i_2    ±1     2, 7, 12, 17, 22, 27, 32, 37
//! i_3    ±1     3, 8, 13, 18, 23, 28, 33, 38  (jx = 0)
//!         OR    4, 9, 14, 19, 24, 29, 34, 39  (jx = 1)
//! ```
//!
//! Per clause 3.8.2 (spec EPUB images `eq61.jpg` / `eq62.jpg`):
//!
//! ```text
//! S = s_0 + 2·s_1 + 4·s_2 + 8·s_3                          (eq 61)
//! C = (m_0 / 5)
//!   + 8·(m_1 / 5)
//!   + 64·(m_2 / 5)
//!   + 512·(2·(m_3 / 5) + jx)                               (eq 62)
//! ```
//!
//! where `s_k = 1` if the sign of pulse `k` is positive, `s_k = 0`
//! otherwise (clause 3.8.2 prose), and `jx = 0` if `m_3 ∈ {3, 8, …,
//! 38}`, `jx = 1` if `m_3 ∈ {4, 9, …, 39}` (clause 3.8.2 prose
//! immediately following eq (62)).
//!
//! Per clause 3.8 (spec EPUB image `eq45.jpg`):
//!
//! ```text
//! c(n) = s_0·δ(n − m_0) + s_1·δ(n − m_1)
//!      + s_2·δ(n − m_2) + s_3·δ(n − m_3)         n = 0, …, 39    (eq 45)
//! ```
//!
//! where `s_k ∈ {-1, +1}` is the signed unit amplitude (clause 3.8
//! prose: "Each pulse can have either the amplitudes +1 or –1").
//!
//! Per clause 4.1.4 verbatim:
//!
//! > "The received fixed-codebook index *C* is used to extract the
//! > positions of the excitation pulses. The pulse signs are obtained
//! > from *S*. This is done by reversing the process described in
//! > clause 3.8.2."
//!
//! That reversal is exactly the inverse of eq (61) / eq (62) above,
//! and is what this module wires.
//!
//! ## What this module does NOT do
//!
//! - It does NOT apply the spec §3.8 eq (48) / (48a) **pitch
//!   sharpening** that the decoder applies when the integer pitch
//!   delay `int(T) < 40`. Per clause 4.1.4: "If the integer part of
//!   the pitch delay *T* is less than the subframe size 40, *c*(*n*)
//!   is modified according to equation (48)." That step needs the
//!   round-255 [`crate::pitch_decode`] output AND the quantised
//!   adaptive-codebook gain `β` of the previous subframe; it is left
//!   for a follow-up round.
//! - It does NOT apply the §3.8 eq (46) harmonic pre-filter
//!   `P(z) = 1/(1 − β·z^−T)`. That filter is applied by the encoder
//!   during the algebraic-codebook search; the decode path described
//!   by clause 4.1.4 only references eqs (45) and (48).

use crate::parameters::{Parameters, C_BITS, S_BITS};

/// G.729 subframe size in samples — the length of the fixed-codebook
/// codevector `c(n)` per spec clause 2.1 / eq (45).
pub const SUBFRAME_SIZE: usize = 40;

/// Number of pulses in one fixed-codebook codevector per spec
/// clause 3.8 / Table 7 (`i_0, i_1, i_2, i_3`).
pub const NUM_PULSES: usize = 4;

/// Pulses-per-track interleave factor per spec clause 3.8 / Table 7.
/// All four tracks have 8 allowed positions on the 5-step lattice.
pub const TRACK_STRIDE: usize = 5;

/// Per-track field width in the `C` codeword per spec eq (62).
/// Tracks 0 / 1 / 2 carry the 3-bit position index `m_k / 5`; track 3
/// carries an extra `jx` track-selector bit on top of its 3-bit
/// position index, giving 4 bits.
pub const C_FIELD_BITS_TRACK_0: usize = 3;
/// Per-track field width for track 1.
pub const C_FIELD_BITS_TRACK_1: usize = 3;
/// Per-track field width for track 2.
pub const C_FIELD_BITS_TRACK_2: usize = 3;
/// Per-track field width for track 3 (3 position bits + 1 `jx` track
/// selector bit per spec eq (62)).
pub const C_FIELD_BITS_TRACK_3: usize = 4;

/// Decoded pulse — a `(position, sign)` pair per spec eq (45). The
/// codevector `c(n)` is the sum of [`NUM_PULSES`] such impulses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pulse {
    /// Pulse position `m_k` per spec Table 7. Lies in
    /// `0..SUBFRAME_SIZE` and is congruent to `k mod TRACK_STRIDE`
    /// for tracks 0..=2; track 3 is congruent to `3` or `4`
    /// (`mod TRACK_STRIDE`) depending on the `jx` bit.
    pub position: u8,
    /// Pulse sign per spec clause 3.8: either `+1` or `-1`.
    pub sign: i8,
}

/// Decoded per-subframe fixed-codebook pulses per spec eq (45) and
/// Table 7. Indexed in spec order `(i_0, i_1, i_2, i_3)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FixedCodebookPulses {
    /// The four signed pulses that make up the codevector. Indexed
    /// `[0..NUM_PULSES]` in spec `(i_0, i_1, i_2, i_3)` order.
    pub pulses: [Pulse; NUM_PULSES],
    /// `jx` per spec eq (62) — `0` if track 3 selected positions
    /// `{3, 8, …, 38}`, `1` if track 3 selected `{4, 9, …, 39}`.
    /// Preserved so a downstream encoder can rebuild the `C`
    /// codeword from the decoded `(positions, jx)`.
    pub jx: u8,
}

/// Decoded fixed-codebook excitation for one transmitted G.729 frame
/// per spec §4.1.4. Carries the per-subframe pulse positions and
/// signs in spec §4.1.5 order: subframe 1 first, subframe 2 next.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FrameFixedCodebook {
    /// Subframe-1 decoded pulses from `(C1, S1)`.
    pub subframe_1: FixedCodebookPulses,
    /// Subframe-2 decoded pulses from `(C2, S2)`.
    pub subframe_2: FixedCodebookPulses,
}

/// Errors returned by the §4.1.4 fixed-codebook decode entry points.
///
/// The transmitted `C` and `S` codewords carry 13 + 4 = 17 bits per
/// subframe; a well-formed unpack from
/// [`crate::parameters::unpack_parameters`] never produces an
/// out-of-domain input. The bounds-checked entry points are provided
/// so callers driving the decoder from raw integers (fuzzing,
/// fixture builders) cannot panic on invalid input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixedCodebookError {
    /// `C` carried more than [`C_BITS`] bits.
    CTooWide { value: u16 },
    /// `S` carried more than [`S_BITS`] bits.
    STooWide { value: u8 },
}

impl core::fmt::Display for FixedCodebookError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::CTooWide { value } => write!(
                f,
                "g729 §4.1.4 fixed-codebook decode: C = {value} exceeds the {C_BITS}-bit field",
            ),
            Self::STooWide { value } => write!(
                f,
                "g729 §4.1.4 fixed-codebook decode: S = {value} exceeds the {S_BITS}-bit field",
            ),
        }
    }
}

impl std::error::Error for FixedCodebookError {}

/// Decode the four pulse signs from the 4-bit `S` codeword per spec
/// eq (61) (clause 3.8.2).
///
/// Per spec eq (61): `S = s_0 + 2·s_1 + 4·s_2 + 8·s_3`, where each
/// `s_k = 1` if the sign of pulse `k` is positive, `s_k = 0`
/// otherwise (clause 3.8.2 prose). Bit `k` of `S` therefore carries
/// the sign indicator for pulse `k`; the returned `i8` is `+1` for
/// a positive sign and `-1` for a negative sign per spec eq (45).
///
/// Returns [`FixedCodebookError::STooWide`] if `s` carries more
/// than [`S_BITS`] bits — `unpack_parameters` never produces such
/// input, but the bounds check exists so a caller driving the
/// decode from raw integers cannot bypass the field width.
#[inline]
pub fn decode_signs(s: u8) -> Result<[i8; NUM_PULSES], FixedCodebookError> {
    if (s as u32) >= (1u32 << S_BITS) {
        return Err(FixedCodebookError::STooWide { value: s });
    }
    let mut signs = [0i8; NUM_PULSES];
    for (k, slot) in signs.iter_mut().enumerate() {
        // Per eq (61) prose: s_k = 1 ⇒ positive (+1), s_k = 0 ⇒ negative (-1).
        let bit = (s >> k) & 1;
        *slot = if bit == 1 { 1 } else { -1 };
    }
    Ok(signs)
}

/// Decode the four pulse positions and the `jx` track-selector bit
/// from the 13-bit `C` codeword per spec eq (62) (clause 3.8.2).
///
/// Per spec eq (62):
/// `C = (m_0 / 5) + 8·(m_1 / 5) + 64·(m_2 / 5) + 512·(2·(m_3 / 5) + jx)`,
/// so the bit-field layout (LSB → MSB) is:
///
/// ```text
/// bits  0..3   →   k_0 = m_0 / 5         (3 bits, 0..=7)
/// bits  3..6   →   k_1 = m_1 / 5         (3 bits, 0..=7)
/// bits  6..9   →   k_2 = m_2 / 5         (3 bits, 0..=7)
/// bits  9..13  →   2·(m_3 / 5) + jx      (4 bits, 0..=15)
///                          where jx is the LSB of the field
/// ```
///
/// Per spec Table 7 / clause 3.8.2 prose, the `m_3` track is
/// `{3, 8, …, 38}` if `jx = 0` and `{4, 9, …, 39}` if `jx = 1`.
/// Pulses on tracks 0..=2 occupy fixed `(0, 1, 2)` modular residues.
///
/// Returns the `[m_0, m_1, m_2, m_3]` positions in spec
/// `(i_0, i_1, i_2, i_3)` order alongside `jx`.
///
/// Returns [`FixedCodebookError::CTooWide`] if `c` carries more
/// than [`C_BITS`] bits.
#[inline]
pub fn decode_positions(c: u16) -> Result<([u8; NUM_PULSES], u8), FixedCodebookError> {
    if (c as u32) >= (1u32 << C_BITS) {
        return Err(FixedCodebookError::CTooWide { value: c });
    }
    let c = u32::from(c);
    // Field extractions per eq (62).
    let k0 = (c & 0b111) as u8;
    let k1 = ((c >> C_FIELD_BITS_TRACK_0) & 0b111) as u8;
    let k2 = ((c >> (C_FIELD_BITS_TRACK_0 + C_FIELD_BITS_TRACK_1)) & 0b111) as u8;
    let track3_field = ((c >> (C_FIELD_BITS_TRACK_0 + C_FIELD_BITS_TRACK_1 + C_FIELD_BITS_TRACK_2))
        & 0b1111) as u8;
    // Per eq (62): the 4-bit track-3 field is `2·(m_3 / 5) + jx`.
    let jx = track3_field & 1;
    let k3 = track3_field >> 1;
    // Per Table 7 the per-track positions are `k * TRACK_STRIDE + offset`:
    // track 0 → +0, track 1 → +1, track 2 → +2, track 3 → +(3 + jx).
    let m0 = k0 * TRACK_STRIDE as u8;
    let m1 = k1 * TRACK_STRIDE as u8 + 1;
    let m2 = k2 * TRACK_STRIDE as u8 + 2;
    let m3 = k3 * TRACK_STRIDE as u8 + 3 + jx;
    Ok(([m0, m1, m2, m3], jx))
}

/// Decode the four `(position, sign)` pulses from the 13-bit `C`
/// codeword and the 4-bit `S` codeword per spec eqs (61) and (62).
///
/// This is the §4.1.4 decode-side reversal of clause 3.8.2: it
/// reconstructs the algebraic-codebook pulse layout that the
/// encoder packed into the on-wire `(C, S)` codewords.
///
/// Returns [`FixedCodebookError::CTooWide`] or
/// [`FixedCodebookError::STooWide`] if either input exceeds its
/// spec-stated field width.
#[inline]
pub fn decode_pulses(c: u16, s: u8) -> Result<FixedCodebookPulses, FixedCodebookError> {
    let (positions, jx) = decode_positions(c)?;
    let signs = decode_signs(s)?;
    let mut pulses = [Pulse {
        position: 0,
        sign: 1,
    }; NUM_PULSES];
    for (k, slot) in pulses.iter_mut().enumerate() {
        *slot = Pulse {
            position: positions[k],
            sign: signs[k],
        };
    }
    Ok(FixedCodebookPulses { pulses, jx })
}

/// Build the 40-sample fixed-codebook codevector `c(n)` per spec
/// eq (45) from the decoded pulses.
///
/// Per spec eq (45):
/// `c(n) = s_0·δ(n − m_0) + s_1·δ(n − m_1) + s_2·δ(n − m_2)
///         + s_3·δ(n − m_3)` for `n = 0, …, 39`.
///
/// The output is an `[i8; 40]` carrying exactly four non-zero
/// entries; `c[m_k] = s_k ∈ {-1, +1}`. All other samples are `0`.
/// This is the **unmodified** codevector: the spec eq (48) / (48a)
/// pitch sharpening (applied when `int(T) < 40`) is NOT performed
/// here — see the module-level documentation.
#[inline]
pub fn build_codevector(pulses: &FixedCodebookPulses) -> [i8; SUBFRAME_SIZE] {
    let mut c = [0i8; SUBFRAME_SIZE];
    for p in &pulses.pulses {
        // The pulse position is always in [0, SUBFRAME_SIZE) by
        // construction (decode_positions clamps via k ∈ 0..=7 and
        // 5·k + offset ≤ 5·7 + 4 = 39). The defensive bound below
        // keeps the indexing safe even if a caller hand-builds a
        // `FixedCodebookPulses` with an out-of-range position.
        if (p.position as usize) < SUBFRAME_SIZE {
            c[p.position as usize] = p.sign;
        }
    }
    c
}

/// Decode the per-frame fixed-codebook excitation per spec §4.1.4
/// from the unpacked transmitted parameters.
///
/// Threads `(C1, S1)` into subframe 1 and `(C2, S2)` into
/// subframe 2 per spec §4.1.5. The codewords produced by
/// [`crate::parameters::unpack_parameters`] are always in-domain,
/// so this entry point cannot fail on well-formed input — but the
/// `Result` is returned for consistency with the per-subframe
/// decode entry points (a hand-built `Parameters` struct could
/// still carry out-of-domain integers).
#[inline]
pub fn decode_frame(params: &Parameters) -> Result<FrameFixedCodebook, FixedCodebookError> {
    let subframe_1 = decode_pulses(params.c1, params.s1)?;
    let subframe_2 = decode_pulses(params.c2, params.s2)?;
    Ok(FrameFixedCodebook {
        subframe_1,
        subframe_2,
    })
}

/// Encode the four pulse signs into the 4-bit `S` codeword per spec
/// eq (61). Returns `None` if any sign is not `+1` or `-1`.
///
/// Exposed so the unit tests can round-trip the encode-side mapping
/// without re-deriving the algebra; encoder wire-up will use it.
///
/// Per spec eq (61): `S = s_0 + 2·s_1 + 4·s_2 + 8·s_3` with
/// `s_k = 1` ⇔ positive sign.
#[inline]
pub fn encode_signs(signs: &[i8; NUM_PULSES]) -> Option<u8> {
    let mut s: u8 = 0;
    for (k, sign) in signs.iter().enumerate() {
        let bit = match sign {
            1 => 1u8,
            -1 => 0u8,
            _ => return None,
        };
        s |= bit << k;
    }
    Some(s)
}

/// Encode the four pulse positions into the 13-bit `C` codeword per
/// spec eq (62). Returns `None` if any position is not on its
/// spec-stated Table-7 track.
///
/// The `jx` bit is derived from `m_3` per the spec prose: `jx = 0`
/// if `m_3 ∈ {3, 8, …, 38}`, `jx = 1` if `m_3 ∈ {4, 9, …, 39}`.
///
/// Per spec eq (62):
/// `C = (m_0/5) + 8·(m_1/5) + 64·(m_2/5) + 512·(2·(m_3/5) + jx)`.
#[inline]
pub fn encode_positions(positions: &[u8; NUM_PULSES]) -> Option<u16> {
    // Tracks 0..=2: position must equal `5·k + offset` with offset
    // = track index.
    for (track, expected_offset) in (0..=2).zip([0u8, 1, 2]) {
        let m = positions[track];
        if (m as usize) >= SUBFRAME_SIZE {
            return None;
        }
        if m % TRACK_STRIDE as u8 != expected_offset {
            return None;
        }
    }
    // Track 3: position must equal `5·k + 3` or `5·k + 4`.
    let m3 = positions[3];
    if (m3 as usize) >= SUBFRAME_SIZE {
        return None;
    }
    let m3_resid = m3 % TRACK_STRIDE as u8;
    let jx = match m3_resid {
        3 => 0u8,
        4 => 1u8,
        _ => return None,
    };
    let k0 = u16::from(positions[0] / TRACK_STRIDE as u8);
    let k1 = u16::from(positions[1] / TRACK_STRIDE as u8);
    let k2 = u16::from(positions[2] / TRACK_STRIDE as u8);
    let k3 = u16::from(positions[3] / TRACK_STRIDE as u8);
    let c = k0
        + (k1 << C_FIELD_BITS_TRACK_0)
        + (k2 << (C_FIELD_BITS_TRACK_0 + C_FIELD_BITS_TRACK_1))
        + ((k3 * 2 + u16::from(jx))
            << (C_FIELD_BITS_TRACK_0 + C_FIELD_BITS_TRACK_1 + C_FIELD_BITS_TRACK_2));
    Some(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Spec eq (61) worked examples — the LSB-first sign packing
    /// pinned for boundary values so any drift in the recipe trips
    /// immediately.
    #[test]
    fn decode_signs_spec_examples() {
        // S = 0b0000 → all negative.
        assert_eq!(decode_signs(0b0000).unwrap(), [-1, -1, -1, -1]);
        // S = 0b1111 → all positive.
        assert_eq!(decode_signs(0b1111).unwrap(), [1, 1, 1, 1]);
        // S = 0b0001 → s_0 positive, s_1..3 negative.
        assert_eq!(decode_signs(0b0001).unwrap(), [1, -1, -1, -1]);
        // S = 0b1000 → s_3 positive, others negative.
        assert_eq!(decode_signs(0b1000).unwrap(), [-1, -1, -1, 1]);
        // S = 0b0101 → s_0 + s_2 positive, s_1 + s_3 negative.
        assert_eq!(decode_signs(0b0101).unwrap(), [1, -1, 1, -1]);
    }

    /// Rejects an out-of-domain `S`.
    #[test]
    fn decode_signs_rejects_too_wide() {
        // Any value with the high nibble set must be rejected since
        // `S_BITS = 4`.
        assert_eq!(
            decode_signs(0b10000),
            Err(FixedCodebookError::STooWide { value: 0b10000 }),
        );
        assert_eq!(
            decode_signs(0xFF),
            Err(FixedCodebookError::STooWide { value: 0xFF }),
        );
    }

    /// Full-domain encode↔decode round-trip on `S`. Pins both
    /// eq (61) (decode) and the symmetric encode-side mapping
    /// simultaneously.
    #[test]
    fn signs_round_trip_full_domain() {
        for s in 0u8..(1u8 << S_BITS) {
            let signs = decode_signs(s).unwrap();
            let s_back = encode_signs(&signs).unwrap();
            assert_eq!(s, s_back, "round-trip drift on S = {s}");
        }
    }

    /// Spec eq (62) worked examples — boundary values pinned
    /// literally.
    #[test]
    fn decode_positions_spec_examples() {
        // C = 0 → k_0 = k_1 = k_2 = 0, track-3 field = 0
        //       → m_0 = 0, m_1 = 1, m_2 = 2, m_3 = 3 (jx = 0).
        let (m, jx) = decode_positions(0).unwrap();
        assert_eq!(m, [0, 1, 2, 3]);
        assert_eq!(jx, 0);

        // C = 1 → k_0 = 1, rest zero
        //       → m_0 = 5, m_1 = 1, m_2 = 2, m_3 = 3.
        let (m, jx) = decode_positions(1).unwrap();
        assert_eq!(m, [5, 1, 2, 3]);
        assert_eq!(jx, 0);

        // C with track-3 field LSB set → jx = 1, k_3 = 0
        //   field bits 9..13 = 0b0001 → C = 1 << 9 = 512.
        let (m, jx) = decode_positions(512).unwrap();
        assert_eq!(m, [0, 1, 2, 4]);
        assert_eq!(jx, 1);

        // Largest valid C = 0x1FFF = 8191 → all fields saturated:
        //   k_0 = 7, k_1 = 7, k_2 = 7, track-3 = 15 → k_3 = 7, jx = 1
        //   → m_0 = 35, m_1 = 36, m_2 = 37, m_3 = 39.
        let (m, jx) = decode_positions(0x1FFF).unwrap();
        assert_eq!(m, [35, 36, 37, 39]);
        assert_eq!(jx, 1);
    }

    /// Rejects an out-of-domain `C`.
    #[test]
    fn decode_positions_rejects_too_wide() {
        // `C_BITS = 13`; bit 13 set is out-of-domain.
        let c = 1u16 << C_BITS;
        assert_eq!(
            decode_positions(c),
            Err(FixedCodebookError::CTooWide { value: c }),
        );
        assert_eq!(
            decode_positions(0xFFFF),
            Err(FixedCodebookError::CTooWide { value: 0xFFFF }),
        );
    }

    /// Full-domain envelope on `decode_positions`: every C in the
    /// `0..8192` domain produces 4 positions inside `[0, 40)` and
    /// on the spec-stated Table-7 tracks.
    #[test]
    fn decode_positions_full_domain_envelope() {
        for c in 0..(1u16 << C_BITS) {
            let (m, jx) = decode_positions(c).unwrap();
            // Track 0: positions ∈ {0, 5, 10, …, 35}.
            assert!(m[0] < SUBFRAME_SIZE as u8);
            assert_eq!(m[0] % TRACK_STRIDE as u8, 0);
            // Track 1: positions ∈ {1, 6, 11, …, 36}.
            assert!(m[1] < SUBFRAME_SIZE as u8);
            assert_eq!(m[1] % TRACK_STRIDE as u8, 1);
            // Track 2: positions ∈ {2, 7, 12, …, 37}.
            assert!(m[2] < SUBFRAME_SIZE as u8);
            assert_eq!(m[2] % TRACK_STRIDE as u8, 2);
            // Track 3: positions ∈ {3, 8, …, 38} if jx = 0 else
            // {4, 9, …, 39}.
            assert!(m[3] < SUBFRAME_SIZE as u8);
            let resid = m[3] % TRACK_STRIDE as u8;
            assert!(resid == 3 || resid == 4);
            assert_eq!(resid, 3 + jx);
            assert!(jx <= 1);
        }
    }

    /// Full-domain encode↔decode round-trip on `C`. Pins eq (62)
    /// (decode) and the symmetric encode-side mapping simultaneously.
    #[test]
    fn positions_round_trip_full_domain() {
        for c in 0..(1u16 << C_BITS) {
            let (m, _jx) = decode_positions(c).unwrap();
            let c_back = encode_positions(&m).unwrap();
            assert_eq!(c, c_back, "round-trip drift on C = {c}");
        }
    }

    /// `encode_positions` rejects positions not on the Table-7
    /// tracks.
    #[test]
    fn encode_positions_rejects_off_track() {
        // Position on track 0 placed at track-1 residue → reject.
        assert_eq!(encode_positions(&[1, 1, 2, 3]), None);
        // Track-3 position with residue {0, 1, 2} → reject.
        assert_eq!(encode_positions(&[0, 1, 2, 0]), None);
        assert_eq!(encode_positions(&[0, 1, 2, 1]), None);
        assert_eq!(encode_positions(&[0, 1, 2, 2]), None);
        // Out-of-range position → reject.
        assert_eq!(encode_positions(&[40, 1, 2, 3]), None);
    }

    /// `encode_signs` rejects non-`±1` inputs.
    #[test]
    fn encode_signs_rejects_invalid() {
        assert_eq!(encode_signs(&[0, 1, -1, 1]), None);
        assert_eq!(encode_signs(&[2, 1, -1, 1]), None);
    }

    /// `decode_pulses` ties decode_positions + decode_signs together
    /// in spec order.
    #[test]
    fn decode_pulses_threads_positions_and_signs() {
        // C = 0b0001 (k_0 = 1, m_0 = 5, others base), S = 0b1010
        //   → s_0 = -1, s_1 = +1, s_2 = -1, s_3 = +1.
        let p = decode_pulses(0b0001, 0b1010).unwrap();
        assert_eq!(p.jx, 0);
        assert_eq!(
            p.pulses,
            [
                Pulse {
                    position: 5,
                    sign: -1,
                },
                Pulse {
                    position: 1,
                    sign: 1,
                },
                Pulse {
                    position: 2,
                    sign: -1,
                },
                Pulse {
                    position: 3,
                    sign: 1,
                },
            ],
        );
    }

    /// `build_codevector` constructs eq (45) — 4 signed impulses,
    /// zero elsewhere.
    #[test]
    fn build_codevector_eq_45() {
        let p = decode_pulses(0b0001, 0b1010).unwrap();
        let c = build_codevector(&p);
        assert_eq!(c.len(), SUBFRAME_SIZE);
        // Positions 5, 1, 2, 3 with signs -1, +1, -1, +1.
        assert_eq!(c[5], -1);
        assert_eq!(c[1], 1);
        assert_eq!(c[2], -1);
        assert_eq!(c[3], 1);
        // Exactly 4 non-zero samples; all others zero.
        let nonzero = c.iter().filter(|x| **x != 0).count();
        assert_eq!(nonzero, NUM_PULSES);
        // Every non-zero sample is ±1 (spec eq (45) prose).
        for x in &c {
            assert!(*x == 0 || x.unsigned_abs() == 1);
        }
    }

    /// Energy of `c(n)` is exactly `NUM_PULSES` per spec eq (45)
    /// (each ±1 impulse contributes 1 to the sum of squares).
    #[test]
    fn build_codevector_energy_pins_pulse_count() {
        // Sweep a representative slice of the C-domain.
        for c in [0u16, 1, 7, 64, 512, 1023, 4095, 8191] {
            for s in 0u8..16 {
                let pulses = decode_pulses(c, s).unwrap();
                let cv = build_codevector(&pulses);
                let energy: i32 = cv.iter().map(|x| i32::from(*x) * i32::from(*x)).sum();
                assert_eq!(energy, NUM_PULSES as i32, "C = {c}, S = {s}");
            }
        }
    }

    /// All four pulses sit on distinct positions per spec Table 7
    /// (the four tracks are disjoint modulo 5, so even with the
    /// same `k` index no two pulses can collide).
    #[test]
    fn pulses_occupy_distinct_positions_full_domain() {
        // Spot-check across the full C domain (8192 codewords).
        for c in 0..(1u16 << C_BITS) {
            let (m, _jx) = decode_positions(c).unwrap();
            assert_ne!(m[0], m[1]);
            assert_ne!(m[0], m[2]);
            assert_ne!(m[0], m[3]);
            assert_ne!(m[1], m[2]);
            assert_ne!(m[1], m[3]);
            assert_ne!(m[2], m[3]);
        }
    }

    /// `decode_frame` correctly threads `(C1, S1)` into subframe 1
    /// and `(C2, S2)` into subframe 2 (not swapped).
    #[test]
    fn decode_frame_threads_per_subframe_indices() {
        // Build a Parameters struct with distinct C1/C2 + S1/S2
        // so a swap is detectable. Use raw struct construction
        // since the test is verifying the frame-level wiring.
        let params = Parameters {
            l0: 0,
            l1: 0,
            l2: 0,
            l3: 0,
            p1: 0,
            p0: 0,
            c1: 0b0001, // m_0 = 5, others base
            s1: 0b0000, // all negative
            ga1: 0,
            gb1: 0,
            p2: 0,
            c2: 512,    // jx = 1, others base
            s2: 0b1111, // all positive
            ga2: 0,
            gb2: 0,
        };
        let frame = decode_frame(&params).unwrap();
        // Subframe 1: m = [5, 1, 2, 3], all negative.
        assert_eq!(frame.subframe_1.jx, 0);
        assert_eq!(frame.subframe_1.pulses[0].position, 5);
        assert_eq!(frame.subframe_1.pulses[0].sign, -1);
        // Subframe 2: m_3 = 4 (jx = 1), all positive.
        assert_eq!(frame.subframe_2.jx, 1);
        assert_eq!(frame.subframe_2.pulses[3].position, 4);
        assert_eq!(frame.subframe_2.pulses[3].sign, 1);
    }

    /// Static-asserts on the per-track field widths sum to
    /// `C_BITS`. A drift in any of the per-track constants would
    /// trip the build.
    #[test]
    fn c_field_widths_sum_to_c_bits() {
        assert_eq!(
            C_FIELD_BITS_TRACK_0
                + C_FIELD_BITS_TRACK_1
                + C_FIELD_BITS_TRACK_2
                + C_FIELD_BITS_TRACK_3,
            C_BITS,
        );
    }

    /// `NUM_PULSES` matches the spec-stated `S_BITS` field width
    /// (one sign bit per pulse).
    #[test]
    fn one_sign_bit_per_pulse() {
        assert_eq!(NUM_PULSES, S_BITS);
    }
}
