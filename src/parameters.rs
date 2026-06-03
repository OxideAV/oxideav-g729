//! §4.1 transmitted-parameter extraction from the 80-bit per-frame
//! payload that [`crate::serial::parse_frame`] returns as a
//! [`crate::serial::FrameKind::Active`] bit array.
//!
//! Round 191 wired the framing layer (sync / bits-header / per-word
//! `0x007F`/`0x0081` packing for the staged conformance corpus).
//! This module wires the next link: turning that bit array into the
//! 15 typed parameter indices that the §4.1 decode-side procedure
//! consumes.
//!
//! ## Spec source — clause 4.1, Table 8 (06/2012 Recommendation)
//!
//! The 80-bit payload is the concatenation of the 15 codewords listed
//! in spec Table 8, in the order shown there (top-to-bottom). The
//! Table-8 NOTE pins the per-codeword bit-ordering: "the bit stream
//! ordering is reflected by the order in the table; for each
//! parameter, the most significant bit (MSB) is transmitted first."
//!
//! The 15-codeword breakdown, per Table 8, is:
//!
//! ```text
//!  #  codeword  bits  role
//!  1  L0        1     §3.2.4 switched MA-predictor selector
//!  2  L1        7     §3.2.4 LSP stage-1 codebook index
//!  3  L2        5     §3.2.4 LSP stage-2 lower-half index
//!  4  L3        5     §3.2.4 LSP stage-2 upper-half index
//!  5  P1        8     §3.7   subframe-1 pitch delay
//!  6  P0        1     §3.7.2 pitch-delay parity bit
//!  7  C1        13    §3.8   subframe-1 fixed-codebook positions
//!  8  S1        4     §3.8   subframe-1 fixed-codebook signs
//!  9  GA1       3     §3.9.2 subframe-1 GA (stage-1) gain index
//! 10  GB1       4     §3.9.2 subframe-1 GB (stage-2) gain index
//! 11  P2        5     §3.7   subframe-2 differential pitch delay
//! 12  C2        13    §3.8   subframe-2 fixed-codebook positions
//! 13  S2        4     §3.8   subframe-2 fixed-codebook signs
//! 14  GA2       3     §3.9.2 subframe-2 GA (stage-1) gain index
//! 15  GB2       4     §3.9.2 subframe-2 GB (stage-2) gain index
//!                    ----
//!                    80 bits  ←  matches the [`crate::tables::BITS_PER_FRAME`] budget
//! ```
//!
//! The frame-level grouping (LSP first, then subframe-1 block, then
//! subframe-2 block) lines up with the decoder's §4.1.1 → §4.1.3 →
//! §4.1.4 → §4.1.5 control flow.
//!
//! ## Wire-bit ordering inside one codeword
//!
//! The §4.1 unpacker walks the bit array left-to-right (slot 0 first)
//! per Table-8 NOTE; for a `B`-bit codeword starting at bit slot
//! `start`, the integer value is
//! `Σ_{k = 0..B} bits[start + k] · 2^(B - 1 - k)`,
//! i.e. slot `start` carries the MSB and slot `start + B - 1` carries
//! the LSB. The 80-bit array layout in
//! [`crate::serial::FrameKind::Active`] places word `i + 2` of the
//! input at array index `i`, matching the wire order used by the
//! Table-8 NOTE statement.
//!
//! ## What this module does NOT do
//!
//! The output [`Parameters`] struct carries the **raw codeword
//! integers** the spec defines. Downstream §4.1.1 / §4.1.3 / §4.1.4 /
//! §4.1.5 procedures interpret these integers further (e.g. mapping
//! `P1` → fractional pitch delay `(int, frac)` via eq (78a), mapping
//! `C1` → pulse positions). Those interpretation steps are wired in
//! their respective modules as they land.

use crate::serial::FrameKind;
use crate::tables::{self, BITS_PER_FRAME, L0_BITS, L1_BITS, L2_BITS, L3_BITS, LSP_TOTAL_BITS};

/// Per-codeword bit width of the §3.7 subframe-1 pitch delay `P1` per
/// spec Table 8.
pub const P1_BITS: usize = 8;

/// Per-codeword bit width of the §3.7.2 pitch-delay parity bit `P0`
/// per spec Table 8.
pub const P0_BITS: usize = 1;

/// Per-codeword bit width of the §3.8 fixed-codebook position field
/// `C1` / `C2` per spec Table 8 (per subframe). The 13 bits split
/// per §3.8 into 4 pulse-position fields (3 + 3 + 3 + 4 with the last
/// pulse carrying the extra `jx` track-selector bit).
pub const C_BITS: usize = 13;

/// Per-codeword bit width of the §3.8 fixed-codebook sign field
/// `S1` / `S2` per spec Table 8 (one bit per pulse, four pulses per
/// subframe).
pub const S_BITS: usize = 4;

/// Per-codeword bit width of the §3.9.2 stage-1 gain-VQ index
/// `GA1` / `GA2` per spec Table 8.
pub const GA_BITS: usize = 3;

/// Per-codeword bit width of the §3.9.2 stage-2 gain-VQ index
/// `GB1` / `GB2` per spec Table 8.
pub const GB_BITS: usize = 4;

/// Per-codeword bit width of the §3.7 subframe-2 differential pitch
/// delay `P2` per spec Table 8.
pub const P2_BITS: usize = 5;

/// Per-frame fixed-codebook total: `C` + `S` per subframe,
/// `2 × (C_BITS + S_BITS) = 2 × (13 + 4) = 34`.
pub const FIXED_CODEBOOK_BITS_PER_FRAME: usize = 2 * (C_BITS + S_BITS);

/// Per-frame gain-VQ total: `GA + GB` per subframe,
/// `2 × (GA_BITS + GB_BITS) = 2 × (3 + 4) = 14`.
pub const GAIN_QUANT_BITS_PER_FRAME: usize = 2 * (GA_BITS + GB_BITS);

/// Per-frame pitch total: `P1 + P0 + P2 = 8 + 1 + 5 = 14`.
pub const PITCH_BITS_PER_FRAME: usize = P1_BITS + P0_BITS + P2_BITS;

/// Decoded indices for the 15 codewords carried by one transmitted
/// G.729 frame, per spec Table 8.
///
/// Each field carries the raw integer value of the codeword as read
/// MSB-first from the wire. The struct is `Copy` since every field is
/// a small integer; pass it by value through the §4.1 procedures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Parameters {
    /// `L0` — 1-bit predictor-mode selector for §3.2.4 LSP
    /// reconstruction. Domain: `0..=1` (1-bit field).
    pub l0: u8,
    /// `L1` — 7-bit first-stage VQ index into
    /// [`crate::tables::LSP_QUANT_CODEBOOK_L1_Q13`]
    /// (`NC0 = 128` rows). Domain: `0..crate::tables::NC0`.
    pub l1: u8,
    /// `L2` — 5-bit second-stage lower-half VQ index. Selects the
    /// lower-5 contribution of [`crate::tables::lsp_l2_entry`].
    /// Domain: `0..crate::tables::NC1`.
    pub l2: u8,
    /// `L3` — 5-bit second-stage upper-half VQ index. Selects the
    /// upper-5 contribution via [`crate::tables::lsp_l3_entry`].
    /// Domain: `0..crate::tables::NC1`.
    pub l3: u8,
    /// `P1` — 8-bit subframe-1 pitch-delay codeword. Maps to a
    /// fractional pitch delay `T1` via spec eq (78a). Domain:
    /// `0..=255`.
    pub p1: u8,
    /// `P0` — 1-bit pitch-delay parity. Per spec §3.7.2 this is
    /// the XOR of the 6 MSBs of `P1`. The §4.1.2 procedure
    /// recomputes it and signals frame erasure on mismatch.
    /// Domain: `0..=1`.
    pub p0: u8,
    /// `C1` — 13-bit subframe-1 fixed-codebook position field.
    /// Per §3.8.2 the field splits as 3+3+3+4 bits across four
    /// pulse-position slots (the trailing pulse carries an extra
    /// `jx` track-selector bit). Domain: `0..=8191`.
    pub c1: u16,
    /// `S1` — 4-bit subframe-1 fixed-codebook sign field, one bit
    /// per pulse. Domain: `0..=15`.
    pub s1: u8,
    /// `GA1` — 3-bit stage-1 gain-VQ index into the GA codebook
    /// (8 entries). Domain: `0..=7`.
    pub ga1: u8,
    /// `GB1` — 4-bit stage-2 gain-VQ index into the GB codebook
    /// (16 entries). Domain: `0..=15`.
    pub gb1: u8,
    /// `P2` — 5-bit subframe-2 differential pitch-delay codeword.
    /// Maps to `T2` relative to `int(T1)` per spec eqs (79)/(80).
    /// Domain: `0..=31`.
    pub p2: u8,
    /// `C2` — 13-bit subframe-2 fixed-codebook position field.
    /// Same 3+3+3+4 split as `C1`. Domain: `0..=8191`.
    pub c2: u16,
    /// `S2` — 4-bit subframe-2 fixed-codebook sign field.
    /// Domain: `0..=15`.
    pub s2: u8,
    /// `GA2` — 3-bit subframe-2 stage-1 gain-VQ index.
    /// Domain: `0..=7`.
    pub ga2: u8,
    /// `GB2` — 4-bit subframe-2 stage-2 gain-VQ index.
    /// Domain: `0..=15`.
    pub gb2: u8,
}

/// Errors returned by [`unpack_parameters`].
///
/// `unpack_parameters` itself only ever consumes the spec-stated
/// 80-bit payload of an active frame; the only failure mode is the
/// caller passing the wrong sentinel ([`FrameKind::Erased`], for
/// which a separate §4.4 concealment path applies).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParameterError {
    /// The frame is a §4.4 erasure sentinel. The §4.1 decode-side
    /// procedure does not consume bits for an erasure; callers must
    /// invoke the §4.4 concealment path instead. Returned from
    /// [`unpack_parameters`] when the input is
    /// [`FrameKind::Erased`].
    Erased,
}

impl core::fmt::Display for ParameterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Erased => write!(
                f,
                "g729 §4.1 parameter unpack: frame is an erasure sentinel; §4.4 concealment path applies"
            ),
        }
    }
}

impl std::error::Error for ParameterError {}

/// Per-codeword §4.1 Table-8 layout: `(bit_width, start_offset)`
/// where `start_offset` is the slot index inside the 80-bit array at
/// which the codeword begins (MSB first; spec Table 8 NOTE).
///
/// The layout is built once as a `const` and is verified at compile
/// time to sum to [`BITS_PER_FRAME`] (a static assert at the module
/// surface; see [`Layout::WIDTHS`] and the build-time `const`
/// computation).
struct Layout;

impl Layout {
    /// Bit widths in Table-8 top-to-bottom order
    /// (`L0, L1, L2, L3, P1, P0, C1, S1, GA1, GB1, P2, C2, S2, GA2, GB2`).
    const WIDTHS: [usize; 15] = [
        L0_BITS, L1_BITS, L2_BITS, L3_BITS, // LSP block
        P1_BITS, P0_BITS, C_BITS, S_BITS, GA_BITS, GB_BITS, // subframe-1 block
        P2_BITS, C_BITS, S_BITS, GA_BITS, GB_BITS, // subframe-2 block
    ];

    /// Compile-time start offsets of each codeword in the 80-bit
    /// array (built from [`Self::WIDTHS`] by a running prefix sum).
    const OFFSETS: [usize; 15] = {
        let mut out = [0; 15];
        let mut acc = 0;
        let mut i = 0;
        while i < 15 {
            out[i] = acc;
            acc += Self::WIDTHS[i];
            i += 1;
        }
        out
    };

    /// Sum of [`Self::WIDTHS`] — must equal [`BITS_PER_FRAME`].
    /// Verified at compile time by the static assertion below.
    const TOTAL: usize = {
        let mut acc = 0;
        let mut i = 0;
        while i < 15 {
            acc += Self::WIDTHS[i];
            i += 1;
        }
        acc
    };
}

// Static spec invariants — verified at compile time so the layout
// can never silently drift.
const _: () = assert!(
    Layout::TOTAL == BITS_PER_FRAME,
    "Table-8 codeword widths must sum to BITS_PER_FRAME = 80",
);
const _: () = assert!(
    LSP_TOTAL_BITS == L0_BITS + L1_BITS + L2_BITS + L3_BITS,
    "LSP_TOTAL_BITS must equal L0+L1+L2+L3",
);
const _: () = assert!(
    PITCH_BITS_PER_FRAME == P1_BITS + P0_BITS + P2_BITS,
    "pitch total = P1 + P0 + P2",
);
const _: () = assert!(
    FIXED_CODEBOOK_BITS_PER_FRAME == 2 * (C_BITS + S_BITS),
    "fixed-codebook total = 2 · (C + S)",
);
const _: () = assert!(
    GAIN_QUANT_BITS_PER_FRAME == 2 * (GA_BITS + GB_BITS),
    "gain-VQ total = 2 · (GA + GB)",
);

/// Reads `width` bits MSB-first from `bits` starting at `start`. The
/// caller is responsible for `start + width <= bits.len()`; for the
/// spec Table-8 layout this is guaranteed by the static assertion
/// `Layout::TOTAL == BITS_PER_FRAME` above and by the fact that the
/// `Active` array length equals `BITS_PER_FRAME` (a struct invariant
/// of [`FrameKind::Active`]).
#[inline]
fn read_msb_first(bits: &[bool; BITS_PER_FRAME], start: usize, width: usize) -> u16 {
    let mut acc: u16 = 0;
    for k in 0..width {
        acc = (acc << 1) | u16::from(bits[start + k]);
    }
    acc
}

/// Unpacks one 80-bit transmitted G.729 frame into the 15 codeword
/// indices defined by spec Table 8.
///
/// The Table-8 NOTE pins the per-codeword bit ordering: for each
/// codeword the MSB comes first. This function walks the 80-bit
/// array left-to-right, slicing off the spec-stated width of each
/// codeword in turn.
///
/// # Errors
///
/// Returns [`ParameterError::Erased`] if the input is a §4.4 erasure
/// sentinel ([`FrameKind::Erased`]). For an erasure the §4.4
/// concealment path must be invoked instead — no bits are consumed.
///
/// # Examples
///
/// ```no_run
/// use oxideav_g729::parameters::unpack_parameters;
/// use oxideav_g729::serial::parse_frame;
///
/// # fn doctest_ignore(buf: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
/// let frame = parse_frame(buf)?;
/// let params = unpack_parameters(&frame)?;
/// // params.l0, params.l1, params.l2, params.l3 feed the §3.2.4 LSP
/// // reconstructor; params.p1, params.p0, params.p2 feed the §3.7
/// // pitch-delay decoder; etc.
/// # let _ = params;
/// # Ok(())
/// # }
/// ```
pub fn unpack_parameters(frame: &FrameKind) -> Result<Parameters, ParameterError> {
    let bits: &[bool; BITS_PER_FRAME] = match frame {
        FrameKind::Active(bits) => bits.as_ref(),
        FrameKind::Erased => return Err(ParameterError::Erased),
    };
    Ok(unpack_bit_array(bits))
}

/// Lower-level variant of [`unpack_parameters`] that takes the
/// 80-bit array directly. Useful for unit-testing the unpacker
/// without spinning the serial-framing layer.
#[must_use]
pub fn unpack_bit_array(bits: &[bool; BITS_PER_FRAME]) -> Parameters {
    let o = Layout::OFFSETS;
    let w = Layout::WIDTHS;

    // The cast back to `u8` / `u16` per field is bounded by the
    // spec-stated codeword width: `L0 = 1` fits `u8`, `L1 = 7` fits
    // `u8`, `C1 = 13` fits `u16`, etc. The widest field is 13 bits
    // (`C1` / `C2`) so `u16` is sufficient throughout; the narrower
    // fields are narrowed to `u8` at the boundary.
    Parameters {
        l0: read_msb_first(bits, o[0], w[0]) as u8,
        l1: read_msb_first(bits, o[1], w[1]) as u8,
        l2: read_msb_first(bits, o[2], w[2]) as u8,
        l3: read_msb_first(bits, o[3], w[3]) as u8,
        p1: read_msb_first(bits, o[4], w[4]) as u8,
        p0: read_msb_first(bits, o[5], w[5]) as u8,
        c1: read_msb_first(bits, o[6], w[6]),
        s1: read_msb_first(bits, o[7], w[7]) as u8,
        ga1: read_msb_first(bits, o[8], w[8]) as u8,
        gb1: read_msb_first(bits, o[9], w[9]) as u8,
        p2: read_msb_first(bits, o[10], w[10]) as u8,
        c2: read_msb_first(bits, o[11], w[11]),
        s2: read_msb_first(bits, o[12], w[12]) as u8,
        ga2: read_msb_first(bits, o[13], w[13]) as u8,
        gb2: read_msb_first(bits, o[14], w[14]) as u8,
    }
}

impl Parameters {
    /// §3.7.2 pitch-delay parity check.
    ///
    /// The spec text (§3.7.2) reads: "the parity bit is generated
    /// through an XOR operation on the six most-significant bits of
    /// `P1`". The bit-exact reading consistent with every active
    /// frame of the staged conformance corpus is **odd parity over
    /// the seven bits `{P0, six MSBs of P1}`** — equivalently,
    /// `P0 = 1 XOR (six-MSB XOR-reduction of P1)`. Empirically this
    /// match-rate is 100% on every active frame of the
    /// `g729-core/SPEECH.BIT` and `g729a/SPEECH.BIT` clean encoder-
    /// output vectors (3 750 frames each), and the dedicated
    /// `PARITY.BIT` mismatch-exerciser carries a non-zero number of
    /// mismatches under exactly this rule (the §4.1.2 concealment
    /// path is its intended target).
    ///
    /// Returns `true` when the parity matches (no transmission
    /// error detected on the pitch-delay field) and `false` when
    /// the §4.1.2 concealment path must be invoked.
    ///
    /// The "six MSBs of `P1`" are bits 7..=2 of the byte
    /// (i.e. `P1 >> 2`).
    ///
    /// # Spec gap
    ///
    /// The spec prose is ambiguous about the parity initialiser
    /// value (is the running XOR initialised to 0 or 1?). The
    /// corpus-validated reading above is odd parity (initialiser
    /// 1). This convention is pinned against the staged corpus,
    /// not invented from a third-party reference.
    #[must_use]
    pub fn pitch_parity_ok(&self) -> bool {
        let six_msbs = self.p1 >> 2;
        // XOR-reduce the low 6 bits, then invert (odd-parity convention).
        let parity = ((six_msbs & 0x3F).count_ones() ^ 1) & 1;
        u32::from(self.p0) == parity
    }
}

/// Re-export the spec frame size at this module's surface so callers
/// staying in `parameters` do not have to dip into [`crate::tables`]
/// just for the 80-bit budget.
pub use tables::BITS_PER_FRAME as FRAME_BITS;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tables::{NC0, NC1};

    /// All-zeros input must produce all-zero codewords.
    #[test]
    fn all_zero_bits_yield_all_zero_parameters() {
        let bits = Box::new([false; BITS_PER_FRAME]);
        let p = unpack_bit_array(&bits);
        assert_eq!(p.l0, 0);
        assert_eq!(p.l1, 0);
        assert_eq!(p.l2, 0);
        assert_eq!(p.l3, 0);
        assert_eq!(p.p1, 0);
        assert_eq!(p.p0, 0);
        assert_eq!(p.c1, 0);
        assert_eq!(p.s1, 0);
        assert_eq!(p.ga1, 0);
        assert_eq!(p.gb1, 0);
        assert_eq!(p.p2, 0);
        assert_eq!(p.c2, 0);
        assert_eq!(p.s2, 0);
        assert_eq!(p.ga2, 0);
        assert_eq!(p.gb2, 0);
    }

    /// All-ones input must saturate every codeword to its maximum
    /// (`2^width - 1`).
    #[test]
    fn all_ones_bits_saturate_every_codeword() {
        let bits = Box::new([true; BITS_PER_FRAME]);
        let p = unpack_bit_array(&bits);
        assert_eq!(u32::from(p.l0), (1u32 << L0_BITS) - 1);
        assert_eq!(u32::from(p.l1), (1u32 << L1_BITS) - 1);
        assert_eq!(u32::from(p.l2), (1u32 << L2_BITS) - 1);
        assert_eq!(u32::from(p.l3), (1u32 << L3_BITS) - 1);
        assert_eq!(u32::from(p.p1), (1u32 << P1_BITS) - 1);
        assert_eq!(u32::from(p.p0), (1u32 << P0_BITS) - 1);
        assert_eq!(u32::from(p.c1), (1u32 << C_BITS) - 1);
        assert_eq!(u32::from(p.s1), (1u32 << S_BITS) - 1);
        assert_eq!(u32::from(p.ga1), (1u32 << GA_BITS) - 1);
        assert_eq!(u32::from(p.gb1), (1u32 << GB_BITS) - 1);
        assert_eq!(u32::from(p.p2), (1u32 << P2_BITS) - 1);
        assert_eq!(u32::from(p.c2), (1u32 << C_BITS) - 1);
        assert_eq!(u32::from(p.s2), (1u32 << S_BITS) - 1);
        assert_eq!(u32::from(p.ga2), (1u32 << GA_BITS) - 1);
        assert_eq!(u32::from(p.gb2), (1u32 << GB_BITS) - 1);
    }

    /// L0/L1/L2/L3 land within their spec-stated codebook domains.
    /// `NC0 = 2^L1_BITS = 128`, `NC1 = 2^L2_BITS = 2^L3_BITS = 32`.
    #[test]
    fn lsp_codeword_domains_match_codebook_sizes() {
        assert_eq!(NC0, 1usize << L1_BITS);
        assert_eq!(NC1, 1usize << L2_BITS);
        assert_eq!(NC1, 1usize << L3_BITS);
        // Saturated value is in-bounds for the codebooks.
        let bits = Box::new([true; BITS_PER_FRAME]);
        let p = unpack_bit_array(&bits);
        assert!((p.l1 as usize) < NC0);
        assert!((p.l2 as usize) < NC1);
        assert!((p.l3 as usize) < NC1);
    }

    /// One-bit-at-a-time test — flipping bit `k` of the 80-bit array
    /// changes exactly one codeword's value, and that codeword is
    /// the one whose `[start, start + width)` window contains `k`.
    /// This locks in the Table-8 NOTE bit-ordering rule (MSB first;
    /// slot 0 of the array is the L0 field's only bit; slot 79 is
    /// the GB2 field's LSB).
    #[test]
    fn single_bit_flip_changes_exactly_one_codeword() {
        let baseline = Box::new([false; BITS_PER_FRAME]);
        let zero = unpack_bit_array(&baseline);
        let widths = Layout::WIDTHS;
        let offsets = Layout::OFFSETS;
        for k in 0..BITS_PER_FRAME {
            let mut bits = baseline.clone();
            bits[k] = true;
            let p = unpack_bit_array(&bits);
            // Identify which codeword `k` belongs to.
            let mut owner = None;
            for (i, &off) in offsets.iter().enumerate() {
                if k >= off && k < off + widths[i] {
                    owner = Some(i);
                    break;
                }
            }
            let owner = owner.expect("Layout::OFFSETS covers every k in 0..BITS_PER_FRAME");
            let expected_owner_field = 1u16 << (widths[owner] - 1 - (k - offsets[owner]));
            // Each field of `p` matches `zero` except `owner`.
            for i in 0..15 {
                let actual = field_value(&p, i);
                let baseline_v = field_value(&zero, i);
                if i == owner {
                    assert_eq!(
                        actual, expected_owner_field,
                        "bit {k} should set {i}-th codeword to MSB-shifted value"
                    );
                } else {
                    assert_eq!(
                        actual, baseline_v,
                        "bit {k} should not affect codeword #{i}"
                    );
                }
            }
        }
    }

    fn field_value(p: &Parameters, i: usize) -> u16 {
        match i {
            0 => u16::from(p.l0),
            1 => u16::from(p.l1),
            2 => u16::from(p.l2),
            3 => u16::from(p.l3),
            4 => u16::from(p.p1),
            5 => u16::from(p.p0),
            6 => p.c1,
            7 => u16::from(p.s1),
            8 => u16::from(p.ga1),
            9 => u16::from(p.gb1),
            10 => u16::from(p.p2),
            11 => p.c2,
            12 => u16::from(p.s2),
            13 => u16::from(p.ga2),
            14 => u16::from(p.gb2),
            _ => unreachable!(),
        }
    }

    /// Packed-then-unpacked round trip for a hand-chosen parameter
    /// vector. Locks the MSB-first slot mapping bidirectionally.
    #[test]
    fn round_trip_pack_unpack() {
        let target = Parameters {
            l0: 1,
            l1: 0b101_0101,  // 7 bits = 85
            l2: 0b1_1010,    // 5 bits = 26
            l3: 0b0_1011,    // 5 bits = 11
            p1: 0b1100_0011, // 8 bits = 195
            p0: 0,
            c1: 0b1_0101_1100_1011, // 13 bits = 0x15CB
            s1: 0b1001,             // 4 bits
            ga1: 0b110,             // 3 bits
            gb1: 0b1010,            // 4 bits
            p2: 0b1_0110,           // 5 bits = 22
            c2: 0b0_1100_0110_1010, // 13 bits
            s2: 0b0011,             // 4 bits
            ga2: 0b001,             // 3 bits
            gb2: 0b1111,            // 4 bits = 15
        };
        // Pack target into a bit array.
        let bits = pack(&target);
        // Unpack and compare.
        let round = unpack_bit_array(&bits);
        assert_eq!(round, target);
    }

    /// Same as above but verifies the explicit bit positions
    /// claimed by Layout::OFFSETS: L0 at slot 0, L1 at slots 1..=7,
    /// L2 at 8..=12, L3 at 13..=17, P1 at 18..=25, P0 at 26, C1 at
    /// 27..=39, S1 at 40..=43, GA1 at 44..=46, GB1 at 47..=50, P2
    /// at 51..=55, C2 at 56..=68, S2 at 69..=72, GA2 at 73..=75,
    /// GB2 at 76..=79.
    #[test]
    fn layout_offsets_match_documented_table8_positions() {
        assert_eq!(Layout::OFFSETS[0], 0); // L0
        assert_eq!(Layout::OFFSETS[1], 1); // L1
        assert_eq!(Layout::OFFSETS[2], 8); // L2
        assert_eq!(Layout::OFFSETS[3], 13); // L3
        assert_eq!(Layout::OFFSETS[4], 18); // P1
        assert_eq!(Layout::OFFSETS[5], 26); // P0
        assert_eq!(Layout::OFFSETS[6], 27); // C1
        assert_eq!(Layout::OFFSETS[7], 40); // S1
        assert_eq!(Layout::OFFSETS[8], 44); // GA1
        assert_eq!(Layout::OFFSETS[9], 47); // GB1
        assert_eq!(Layout::OFFSETS[10], 51); // P2
        assert_eq!(Layout::OFFSETS[11], 56); // C2
        assert_eq!(Layout::OFFSETS[12], 69); // S2
        assert_eq!(Layout::OFFSETS[13], 73); // GA2
        assert_eq!(Layout::OFFSETS[14], 76); // GB2
                                             // Final slot (after GB2) lands exactly on the 80-bit budget.
        assert_eq!(Layout::OFFSETS[14] + Layout::WIDTHS[14], BITS_PER_FRAME);
    }

    /// Helper for the round-trip test: encode a `Parameters` value
    /// back into the 80-bit array using the same MSB-first rule
    /// that `unpack_bit_array` reads with.
    fn pack(p: &Parameters) -> Box<[bool; BITS_PER_FRAME]> {
        let mut bits = Box::new([false; BITS_PER_FRAME]);
        let widths = Layout::WIDTHS;
        let offsets = Layout::OFFSETS;
        let vals: [u16; 15] = [
            u16::from(p.l0),
            u16::from(p.l1),
            u16::from(p.l2),
            u16::from(p.l3),
            u16::from(p.p1),
            u16::from(p.p0),
            p.c1,
            u16::from(p.s1),
            u16::from(p.ga1),
            u16::from(p.gb1),
            u16::from(p.p2),
            p.c2,
            u16::from(p.s2),
            u16::from(p.ga2),
            u16::from(p.gb2),
        ];
        for i in 0..15 {
            let off = offsets[i];
            let w = widths[i];
            let v = vals[i];
            for k in 0..w {
                let bit = ((v >> (w - 1 - k)) & 1) == 1;
                bits[off + k] = bit;
            }
        }
        bits
    }

    /// §3.7.2 parity rule under the odd-parity reading: P0 is the
    /// complement of the XOR of the 6 MSBs of P1, equivalently
    /// `(P0 XOR XOR_reduce(six_MSBs)) == 1`. Worked checks on
    /// three crafted P1 values.
    #[test]
    fn pitch_parity_xor_of_p1_six_msbs_with_odd_parity_initialiser() {
        // P1 = 0b1100_0000 → six MSBs = 0b110000 → XOR-reduce = 0
        // → complement = 1, so P0 = 1 holds, P0 = 0 fails.
        let mut p = baseline_parameters();
        p.p1 = 0b1100_0000;
        p.p0 = 1;
        assert!(p.pitch_parity_ok());
        p.p0 = 0;
        assert!(!p.pitch_parity_ok());

        // P1 = 0b1110_0011 → six MSBs = 0b111000 → XOR-reduce = 1
        // → complement = 0, so P0 = 0 holds.
        p.p1 = 0b1110_0011;
        p.p0 = 0;
        assert!(p.pitch_parity_ok());
        p.p0 = 1;
        assert!(!p.pitch_parity_ok());

        // P1 = 0b1111_1111 → six MSBs = 0b111111 → XOR-reduce = 0
        // → complement = 1, so P0 = 1 holds.
        p.p1 = 0b1111_1111;
        p.p0 = 1;
        assert!(p.pitch_parity_ok());
    }

    fn baseline_parameters() -> Parameters {
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

    /// An erasure-sentinel frame should be rejected by
    /// `unpack_parameters`. The §4.4 concealment path applies
    /// instead.
    #[test]
    fn unpack_parameters_rejects_erasure() {
        let err = unpack_parameters(&FrameKind::Erased).unwrap_err();
        assert_eq!(err, ParameterError::Erased);
    }

    /// `unpack_parameters` on an active frame matches the
    /// lower-level `unpack_bit_array` on the same bits.
    #[test]
    fn unpack_parameters_matches_unpack_bit_array() {
        let bits = Box::new([true; BITS_PER_FRAME]);
        let from_array = unpack_bit_array(&bits);
        let from_frame =
            unpack_parameters(&FrameKind::Active(bits.clone())).expect("active frame unpacks");
        assert_eq!(from_array, from_frame);
    }
}
