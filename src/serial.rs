//! ITU-T G.729 serial bitstream framing — the on-wire format used by
//! the ITU electronic-attachment reference encoder/decoder for the
//! `.bit` conformance test sequences staged under
//! `docs/audio/g729/conformance/`.
//!
//! ## Frame layout
//!
//! Each frame is a fixed `164`-byte little-endian `Word16` packet:
//!
//! ```text
//!   word index 0   : sync word         = 0x6B21         (2 bytes)
//!   word index 1   : per-frame bitcount = 0x0050 (= 80) (2 bytes)
//!   word index 2..82 : 80 bit-value words, each       (160 bytes)
//!                      either 0x007F ("0") or 0x0081 ("1")
//!                                                       --------
//!                                                       164 bytes
//! ```
//!
//! Total payload: `2 + 80 = 82` Word16 entries = `164` bytes per
//! 10 ms / 80-PCM-sample G.729 frame, matching the per-frame bit budget
//! of clause 2 / Table 1 of ITU-T G.729 (06/2012). The `80` bit-count
//! header word equals the spec frame size [`crate::tables::BITS_PER_FRAME`].
//!
//! ## Provenance of the framing constants
//!
//! [`SYNC_WORD`], [`BIT_ZERO`], [`BIT_ONE`] are **empirically observed
//! constants** of the staged `docs/audio/g729/conformance/` byte stream
//! itself — the literals occupy fixed word positions in every `.bit`
//! file (`ALGTHM.BIT`, `FIXED.BIT`, `LSP.BIT`, `PITCH.BIT`, `SPEECH.BIT`,
//! `TAME.BIT` of the base codec; the Annex-A and Annex-B sets share the
//! same wrapper). They are not algorithmic in nature — they are the
//! fixed labels the reference encoder writes alongside its 80
//! transmitted bits so the reference decoder can sanity-check a frame
//! boundary before consuming the bits. The 164-byte cadence is also
//! documented in `docs/audio/g729/conformance/README.md` (the docs
//! collaborator's self-consistency cross-check). No algorithmic source
//! from the ITU electronic attachment was read to determine these
//! constants — only the data files themselves and the docs
//! collaborator's framing-byte description.
//!
//! ## What this module does NOT do
//!
//! The 80 transmitted bits per frame are an **opaque payload** at this
//! layer: this module unpacks the per-frame bit array but does not
//! interpret the bits per the §4.1 Table-8 parameter packing. Wiring
//! the bit-array to LSP / pitch / fixed-codebook indices is a future
//! round, gated on the same docs-collaborator specifier pass that
//! gates the remaining bit-exact codebook tables.

/// Sync word at word-0 of every G.729 ITU serial bitstream frame.
///
/// Empirically the first two bytes of every `.bit` file in
/// `docs/audio/g729/conformance/` are `0x21 0x6B` (little-endian),
/// decoding to the `Word16` value `0x6B21`. The same literal appears
/// at the start of every subsequent 164-byte frame.
pub const SYNC_WORD: u16 = 0x6B21;

/// Per-frame bit-count header at word-1 of every frame. Equals the
/// G.729 per-frame bit budget [`crate::tables::BITS_PER_FRAME`] = 80,
/// encoded as the `Word16` value `0x0050`.
pub const BITS_HEADER: u16 = 80;

/// Word value the reference encoder writes for a transmitted **0**
/// bit, observed at every non-sync, non-length word position in the
/// staged `.bit` files where the spec-defined bit is `0`.
pub const BIT_ZERO: u16 = 0x007F;

/// Word value the reference encoder writes for a transmitted **1**
/// bit, observed at every non-sync, non-length word position in the
/// staged `.bit` files where the spec-defined bit is `1`.
pub const BIT_ONE: u16 = 0x0081;

/// Per-word marker the reference encoder writes for a **frame
/// erasure** sentinel: when every word in the 80-word payload of an
/// otherwise-valid frame (sync + bits-header still 0x6B21 / 0x0050) is
/// this value, the frame is signalling "drop and conceal" rather than
/// carrying real bits. Empirically observed in the staged
/// `docs/audio/g729/conformance/{g729-core,g729a}/ERASURE.BIT`
/// sequences (60 of 300 frames marked this way), exercising the §4.4
/// concealment path of [`crate::tables`] §4.4.
///
/// The docs collaborator's `docs/audio/g729/conformance/README.md`
/// classifies `ERASURE` as a decoder-only test but does not document
/// the on-wire marker shape; that gap is reported in the round-191
/// follow-up notes.
pub const BIT_ERASED: u16 = 0x0000;

/// Fixed Word16 count of one ITU serial frame: 2 header words +
/// [`crate::tables::BITS_PER_FRAME`] bit-value words.
pub const FRAME_WORDS: usize = 2 + crate::tables::BITS_PER_FRAME;

/// Fixed byte count of one ITU serial frame:
/// [`FRAME_WORDS`] × 2 = 164 bytes.
pub const FRAME_BYTES: usize = FRAME_WORDS * 2;

/// Outcome of [`parse_frame`] for a well-formed 164-byte ITU serial
/// frame. The variant distinguishes a normal bit-carrying frame from
/// a frame-erasure sentinel ([`BIT_ERASED`] payload).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FrameKind {
    /// A normal transmitted frame; the inner array carries the 80
    /// bits in wire order (word index 2 = bit 0, etc.).
    Active(Box<[bool; crate::tables::BITS_PER_FRAME]>),
    /// A frame-erasure sentinel: sync + header still valid, but every
    /// bit-payload word equals [`BIT_ERASED`]. Indicates the §4.4
    /// concealment path should produce this frame's output samples
    /// without consuming any actual transmitted bits.
    Erased,
}

/// Errors returned by [`parse_frame`] when a 164-byte input does not
/// conform to the ITU serial framing of `docs/audio/g729/conformance/`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SerialError {
    /// The input slice was not exactly [`FRAME_BYTES`] = 164 bytes.
    /// Carries the actual length the caller passed.
    WrongLength(usize),
    /// Word-0 was not [`SYNC_WORD`]. Carries the observed value.
    BadSync(u16),
    /// Word-1 was not [`BITS_HEADER`] (= 80). Carries the observed value.
    BadBitsHeader(u16),
    /// A bit-payload word (word indices 2..82) was neither
    /// [`BIT_ZERO`], [`BIT_ONE`], nor [`BIT_ERASED`] (in the
    /// all-erased-payload pattern of an erasure-sentinel frame).
    /// Carries the offending word index
    /// (`0..crate::tables::BITS_PER_FRAME`) and the observed value.
    BadBitWord { index: usize, value: u16 },
    /// The frame had a mix of erasure-marker words ([`BIT_ERASED`])
    /// and normal `0x007F`/`0x0081` bit words. The ITU encoder either
    /// writes all-erased payload or no erased words at all; a
    /// mid-frame mix is a corruption signal.
    MixedErasure { erased: usize, normal: usize },
}

impl core::fmt::Display for SerialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::WrongLength(n) => write!(
                f,
                "g729 serial frame must be {FRAME_BYTES} bytes; got {n}"
            ),
            Self::BadSync(v) => write!(
                f,
                "g729 serial sync word mismatch: expected 0x{SYNC_WORD:04X}, got 0x{v:04X}"
            ),
            Self::BadBitsHeader(v) => write!(
                f,
                "g729 serial bits-per-frame header mismatch: expected {BITS_HEADER}, got 0x{v:04X}"
            ),
            Self::BadBitWord { index, value } => write!(
                f,
                "g729 serial bit-word #{index} not 0x{BIT_ZERO:04X}/0x{BIT_ONE:04X}: got 0x{value:04X}"
            ),
            Self::MixedErasure { erased, normal } => write!(
                f,
                "g729 serial frame mixes erasure-marker words ({erased}) with normal bits ({normal}); ITU encoder writes either all-erased or none"
            ),
        }
    }
}

impl std::error::Error for SerialError {}

/// Parses one 164-byte ITU serial G.729 frame, distinguishing a
/// normal bit-carrying frame ([`FrameKind::Active`]) from a frame-
/// erasure sentinel ([`FrameKind::Erased`]).
///
/// For an active frame, the bit at index `i` of the inner array comes
/// from word `i + 2` of the input (word indices 0 and 1 are the
/// [`SYNC_WORD`] / [`BITS_HEADER`] header). The ordering matches the
/// wire ordering; per the spec Table-8 NOTE this is the §4.1-Table-8
/// top-to-bottom MSB-first packing for the 11 transmitted parameters,
/// but **this function does not unpack that further** — it is a pure
/// framing parser. Callers that want parameter indices apply the
/// Table-8 split themselves (a future-round wiring task).
///
/// # Errors
///
/// Returns [`SerialError`] if any framing constant is violated, or
/// [`SerialError::MixedErasure`] if the 80-word payload mixes
/// [`BIT_ERASED`] markers with normal `0x007F`/`0x0081` words.
pub fn parse_frame(frame_bytes: &[u8]) -> Result<FrameKind, SerialError> {
    if frame_bytes.len() != FRAME_BYTES {
        return Err(SerialError::WrongLength(frame_bytes.len()));
    }

    let word =
        |i: usize| -> u16 { u16::from_le_bytes([frame_bytes[i * 2], frame_bytes[i * 2 + 1]]) };

    let sync = word(0);
    if sync != SYNC_WORD {
        return Err(SerialError::BadSync(sync));
    }

    let header = word(1);
    if header != BITS_HEADER {
        return Err(SerialError::BadBitsHeader(header));
    }

    // First pass: classify each payload word. The ITU encoder writes
    // an entire frame either as bits or as an erasure sentinel — a
    // mid-frame mix would be a corruption.
    let mut erased = 0usize;
    let mut normal = 0usize;
    for i in 0..crate::tables::BITS_PER_FRAME {
        let w = word(i + 2);
        match w {
            BIT_ZERO | BIT_ONE => normal += 1,
            BIT_ERASED => erased += 1,
            other => {
                return Err(SerialError::BadBitWord {
                    index: i,
                    value: other,
                })
            }
        }
    }
    if erased > 0 && normal > 0 {
        return Err(SerialError::MixedErasure { erased, normal });
    }
    if erased == crate::tables::BITS_PER_FRAME {
        return Ok(FrameKind::Erased);
    }

    // All-normal payload — decode the bits.
    let mut bits = Box::new([false; crate::tables::BITS_PER_FRAME]);
    for i in 0..crate::tables::BITS_PER_FRAME {
        let w = word(i + 2);
        bits[i] = match w {
            BIT_ZERO => false,
            BIT_ONE => true,
            // Unreachable: classified as Normal above.
            other => {
                return Err(SerialError::BadBitWord {
                    index: i,
                    value: other,
                })
            }
        };
    }
    Ok(FrameKind::Active(bits))
}

/// Convenience: number of complete frames in a byte buffer whose
/// length is an exact multiple of [`FRAME_BYTES`]. Returns
/// [`SerialError::WrongLength`] if the buffer length is not a multiple
/// of 164 bytes (a partial trailing frame is a corruption signal under
/// the ITU encoder's fixed-cadence output).
pub fn frame_count(bitstream_bytes: &[u8]) -> Result<usize, SerialError> {
    if bitstream_bytes.len() % FRAME_BYTES != 0 {
        return Err(SerialError::WrongLength(bitstream_bytes.len()));
    }
    Ok(bitstream_bytes.len() / FRAME_BYTES)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A synthetic single-frame buffer matching the framing rules:
    /// sync, header(=80), 80 bit-zero words.
    fn synth_all_zero_frame() -> Vec<u8> {
        let mut buf = Vec::with_capacity(FRAME_BYTES);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for _ in 0..crate::tables::BITS_PER_FRAME {
            buf.extend_from_slice(&BIT_ZERO.to_le_bytes());
        }
        buf
    }

    /// Same but with all-ones payload.
    fn synth_all_one_frame() -> Vec<u8> {
        let mut buf = Vec::with_capacity(FRAME_BYTES);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for _ in 0..crate::tables::BITS_PER_FRAME {
            buf.extend_from_slice(&BIT_ONE.to_le_bytes());
        }
        buf
    }

    #[test]
    fn frame_bytes_is_164() {
        assert_eq!(FRAME_BYTES, 164);
        assert_eq!(FRAME_WORDS, 82);
    }

    #[test]
    fn header_equals_bits_per_frame() {
        // Tie the header constant to the spec constant.
        assert_eq!(BITS_HEADER as usize, crate::tables::BITS_PER_FRAME);
    }

    #[test]
    fn parse_all_zero_frame_returns_all_false() {
        let buf = synth_all_zero_frame();
        match parse_frame(&buf).expect("synth frame parses") {
            FrameKind::Active(bits) => assert!(bits.iter().all(|b| !*b)),
            FrameKind::Erased => panic!("expected Active, got Erased"),
        }
    }

    #[test]
    fn parse_all_one_frame_returns_all_true() {
        let buf = synth_all_one_frame();
        match parse_frame(&buf).expect("synth frame parses") {
            FrameKind::Active(bits) => assert!(bits.iter().all(|b| *b)),
            FrameKind::Erased => panic!("expected Active, got Erased"),
        }
    }

    /// A synthetic erasure-sentinel frame: sync + header valid,
    /// every payload word is [`BIT_ERASED`].
    fn synth_erasure_frame() -> Vec<u8> {
        let mut buf = Vec::with_capacity(FRAME_BYTES);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for _ in 0..crate::tables::BITS_PER_FRAME {
            buf.extend_from_slice(&BIT_ERASED.to_le_bytes());
        }
        buf
    }

    #[test]
    fn parse_erasure_frame_returns_erased_variant() {
        let buf = synth_erasure_frame();
        match parse_frame(&buf).expect("erasure frame parses") {
            FrameKind::Erased => {}
            FrameKind::Active(_) => panic!("expected Erased, got Active"),
        }
    }

    #[test]
    fn parse_rejects_mixed_erasure() {
        let mut buf = synth_erasure_frame();
        // Flip word index 2's first byte to 0x7F so word 2 = 0x007F
        // (a normal BIT_ZERO) — payload now has 79 erased + 1 normal.
        let word2_byte0 = 2 * 2;
        buf[word2_byte0] = 0x7F;
        buf[word2_byte0 + 1] = 0x00;
        assert_eq!(
            parse_frame(&buf).unwrap_err(),
            SerialError::MixedErasure {
                erased: 79,
                normal: 1,
            },
        );
    }

    #[test]
    fn parse_rejects_wrong_length() {
        assert_eq!(
            parse_frame(&[0u8; 163]).unwrap_err(),
            SerialError::WrongLength(163),
        );
        assert_eq!(
            parse_frame(&[0u8; 165]).unwrap_err(),
            SerialError::WrongLength(165),
        );
        assert_eq!(parse_frame(&[]).unwrap_err(), SerialError::WrongLength(0),);
    }

    #[test]
    fn parse_rejects_bad_sync() {
        let mut buf = synth_all_zero_frame();
        // Corrupt word 0 to 0x0000.
        buf[0] = 0;
        buf[1] = 0;
        assert_eq!(parse_frame(&buf).unwrap_err(), SerialError::BadSync(0));
    }

    #[test]
    fn parse_rejects_bad_header() {
        let mut buf = synth_all_zero_frame();
        // Corrupt word 1 to 0x0040 (= 64, not 80).
        buf[2] = 0x40;
        buf[3] = 0x00;
        assert_eq!(
            parse_frame(&buf).unwrap_err(),
            SerialError::BadBitsHeader(0x0040),
        );
    }

    #[test]
    fn parse_rejects_invalid_bit_word() {
        let mut buf = synth_all_zero_frame();
        // Corrupt the 5th bit-payload word (i.e. word-index 7) to 0x00AA.
        let offset = 7 * 2;
        buf[offset] = 0xAA;
        buf[offset + 1] = 0x00;
        assert_eq!(
            parse_frame(&buf).unwrap_err(),
            SerialError::BadBitWord {
                index: 5,
                value: 0x00AA,
            },
        );
    }

    #[test]
    fn frame_count_multiple_frames() {
        let mut buf = synth_all_zero_frame();
        buf.extend(synth_all_one_frame());
        buf.extend(synth_all_zero_frame());
        assert_eq!(frame_count(&buf).unwrap(), 3);
    }

    #[test]
    fn frame_count_rejects_partial() {
        let mut buf = synth_all_zero_frame();
        buf.push(0); // 165 bytes
        assert_eq!(
            frame_count(&buf).unwrap_err(),
            SerialError::WrongLength(165)
        );
    }

    #[test]
    fn frame_count_empty_is_zero_frames() {
        assert_eq!(frame_count(&[]).unwrap(), 0);
    }
}
