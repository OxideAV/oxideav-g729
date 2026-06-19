//! ITU-T G.729 **Annex B** (silence compression: VAD / DTX / CNG)
//! decoder-side framing and SID parameter decode.
//!
//! Annex B reduces the average bit rate during silence by transmitting
//! the full 80-bit G.729 frame only for active speech. Non-active frames
//! are either a 15-bit **SID** (silence insertion descriptor) frame that
//! refreshes the comfort-noise spectrum/energy, or an **untransmitted**
//! frame carrying no payload at all (the decoder reuses the last SID
//! parameters). This module implements the parts of the §B decoder path
//! that are fully determined by the staged conformance corpus and the
//! Annex B prose:
//!
//! * **Variable-length serial framing** (§B.4.3 / the staged
//!   `docs/audio/g729/conformance/g729b/*.bit` corpus) — the per-frame
//!   bit-count header word selects the frame type.
//! * **SID bitstream unpack** (Table B.2) — the four SID indices
//!   (switched-predictor / first-stage LSF / second-stage LSF / gain).
//! * **§B.4.2.1 energy dequantization** — the 5-bit non-uniform
//!   logarithmic gain quantizer, which the spec states explicitly "does
//!   not need the storage of a quantizer table" (so it is fully
//!   prose-sourced).
//! * **§B.4 / §B.4.5 decoder frame-type classification** including the
//!   erased-frame inheritance rule.
//!
//! ## Annex B serial framing
//!
//! The Annex B `.bit` files use the same `Word16` serial convention as
//! the base codec ([`crate::serial`]) but with a **variable** per-frame
//! payload. Each frame is:
//!
//! ```text
//!   word 0          : sync word      (0x6B21 normal, 0x6B20 erased)
//!   word 1          : payload bit count  n_bits ∈ {0, 16, 80}
//!   word 2 .. 2+n   : n_bits bit-value words (0x007F = 0, 0x0081 = 1)
//! ```
//!
//! The frame type is read directly from the header:
//!
//! | sync     | n_bits | Frame type      | Spec |
//! |----------|-------:|-----------------|------|
//! | `0x6B21` |     80 | active speech   | Table 8 |
//! | `0x6B21` |     16 | SID             | Table B.2 (15 used + 1 pad) |
//! | `0x6B21` |      0 | untransmitted   | §B.4.1 |
//! | `0x6B20` |    any | erased frame    | §B.4.5 |
//!
//! A frame can also be signalled erased a **second** way: the normal
//! sync `0x6B21` with a non-empty payload whose every word is the
//! [`crate::serial::BIT_ERASED`] marker (`0x0000`). Both shapes appear in
//! the staged `g729b/tstseq6.bit` (one erased untransmitted frame with
//! the `0x6B20` sync, and one erased SID-sized frame whose 16 payload
//! words are all `0x0000`); §B.4.5 resolves both to a concealed frame
//! whose type is inherited from the preceding frame.
//!
//! All four header shapes are **empirically observed** across every
//! staged `docs/audio/g729/conformance/g729b/*.bit` sequence
//! (`tstseq1..6` and their Annex-A+B `…a` variants): the sync/erased-sync
//! word values, the `{0, 16, 80}` bit-count set, and the
//! `0x007F`/`0x0081` bit-value markers all occupy fixed word positions in
//! the byte stream. They are not algorithmic in nature — they are the
//! fixed labels the reference encoder writes so the decoder can route a
//! frame to the right path. No algorithmic source from the ITU
//! electronic attachment was read to determine them: only the staged data
//! files themselves and the `conformance/README.md` framing description.
//!
//! ## SID bit-stream (Table B.2)
//!
//! The SID frame carries **15** meaningful bits (transmitted in a 16-word
//! payload; the 16th word is padding), in the Table B.2 order, MSB-first
//! per field:
//!
//! | Field                                  | Bits |
//! |----------------------------------------|-----:|
//! | Switched predictor index of LSF VQ     |    1 |
//! | First stage vector of LSF VQ           |    5 |
//! | Second stage vector of LSF VQ          |    4 |
//! | Gain (energy)                          |    5 |
//!
//! ## What this module does NOT do
//!
//! The §B.4.2.2 SID-LSP vector dequantization (the modified MA predictor,
//! the 32-address first-stage subset, the two 16-address second-stage
//! subset) and the §B.4.4 CNG excitation generation (the Gaussian mixture
//! excitation, the random pitch/fixed-codebook selection, the
//! energy-smoothing constants) require numeric tables that are **not**
//! present under `docs/audio/g729/tables/` — only the CNG spectrum
//! factor/shift and VAD margin tables are staged; the
//! `annexB-cng-lsp-sid-reset-Q15.csv` (`lspSid_reset`) listed in the
//! tables README is **absent**, and the SID-LSP VQ subset codebooks are
//! not extracted at all. Those stages are reported as a precise docs-gap
//! and are out of scope here. This module decodes the SID *energy* (which
//! needs no table) and the SID *indices* (raw codewords), and classifies
//! every frame type so a CNG synthesis stage can be slotted in once the
//! missing tables are staged.

use crate::serial::{BIT_ERASED, BIT_ONE, BIT_ZERO, SYNC_WORD};

/// Sync word the Annex B reference encoder writes for a **frame
/// erasure**: the normal [`SYNC_WORD`] (`0x6B21`) with its low bit
/// cleared (`0x6B20`).
///
/// Empirically observed at the head of the erased frames in the staged
/// `docs/audio/g729/conformance/g729b/tstseq6.bit` sequence (one frame
/// of 100 carries this marker, immediately preceding an untransmitted
/// frame). Per §B.4.5 a decoder treats the erased frame's type as
/// inherited from the preceding frame.
pub const ERASED_SYNC_WORD: u16 = 0x6B20;

/// Payload bit count of an Annex B **active speech** frame (the full
/// G.729 frame; equals [`crate::tables::BITS_PER_FRAME`]).
pub const ACTIVE_BITS: usize = crate::tables::BITS_PER_FRAME;

/// Payload word count of an Annex B **SID** frame. Table B.2 defines 15
/// meaningful bits; the reference encoder transmits them in a 16-word
/// payload (the 16th word is a padding bit), as observed in every staged
/// `g729b/*.bit` SID frame.
pub const SID_WORDS: usize = 16;

/// Number of **meaningful** SID payload bits (Table B.2:
/// `1 + 5 + 4 + 5`); the 16th transmitted word is padding.
pub const SID_BITS: usize = SID_LP0_BITS + SID_L1_BITS + SID_L2_BITS + SID_GAIN_BITS;

/// Table B.2 field width — switched-predictor index of the SID LSF VQ.
pub const SID_LP0_BITS: usize = 1;
/// Table B.2 field width — first-stage vector of the SID LSF VQ.
pub const SID_L1_BITS: usize = 5;
/// Table B.2 field width — second-stage vector of the SID LSF VQ.
pub const SID_L2_BITS: usize = 4;
/// Table B.2 field width — SID gain (energy) index.
pub const SID_GAIN_BITS: usize = 5;

/// §B.4.1 frame-type classification of one Annex B serial frame.
///
/// The DTX module (§B.4.1) outputs one of three transmission types per
/// frame — active speech, SID, or untransmitted — encoded in the serial
/// stream by the per-frame bit count. A fourth case, [`Self::Erased`],
/// is the §B.4.5 frame-erasure sentinel (the [`ERASED_SYNC_WORD`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnnexBFrame {
    /// Active speech frame (`Ftyp = 1`): the full 80-bit G.729 payload.
    /// The inner array carries the 80 bits in wire order, matching
    /// [`crate::serial::FrameKind::Active`]'s payload.
    Active(Box<[bool; ACTIVE_BITS]>),
    /// SID frame (`Ftyp = 2`): a 15-bit silence-insertion descriptor.
    Sid(SidFrame),
    /// Untransmitted frame (`Ftyp = 0`): no payload; the decoder reuses
    /// the last received SID parameters to continue comfort-noise.
    Untransmitted,
    /// §B.4.5 frame-erasure sentinel — either the [`ERASED_SYNC_WORD`]
    /// (`0x6B20`) sync, or a normal-sync frame whose entire payload is the
    /// [`crate::serial::BIT_ERASED`] (`0x0000`) marker. The decoder
    /// derives the concealed frame's type from the preceding frame:
    /// active→active, else→untransmitted.
    Erased,
}

/// The four Table B.2 SID indices, decoded from a SID frame's 15
/// meaningful payload bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SidFrame {
    /// Switched-predictor index of the SID LSF VQ (1 bit, §B.4.2.2
    /// modification 1 — selects which MA predictor combination is used).
    pub lp0: u8,
    /// First-stage SID-LSF VQ index (5 bits; addresses the 32-entry
    /// first-stage subset of §B.4.2.2 modification 2).
    pub l1: u8,
    /// Second-stage SID-LSF VQ index (4 bits; addresses the 16-entry
    /// second-stage subset of §B.4.2.2 modification 3).
    pub l2: u8,
    /// SID gain (energy) index (5 bits; dequantized by
    /// [`dequant_sid_energy_db`]).
    pub gain: u8,
}

impl SidFrame {
    /// Dequantize this SID frame's gain index to the §B.4.2.1 energy
    /// level in dB. See [`dequant_sid_energy_db`].
    #[must_use]
    pub fn energy_db(&self) -> f32 {
        dequant_sid_energy_db(self.gain)
    }
}

/// Errors from [`parse_annex_b_frame`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnexBError {
    /// The input slice was too short to hold even the 2-word header.
    TooShort(usize),
    /// The header advertised a payload that runs past the slice end.
    /// Carries `(needed_bytes, available_bytes)`.
    Truncated { needed: usize, available: usize },
    /// Word-0 was neither [`SYNC_WORD`] nor [`ERASED_SYNC_WORD`].
    /// Carries the observed value.
    BadSync(u16),
    /// The header bit count was not one of the Annex B frame sizes
    /// (`0`, [`SID_WORDS`], or [`ACTIVE_BITS`]). Carries the observed
    /// value.
    BadBitCount(usize),
    /// A bit-payload word was neither [`BIT_ZERO`], [`BIT_ONE`], nor
    /// [`BIT_ERASED`]. Carries `(word_index, observed_value)`.
    BadBitWord { index: usize, value: u16 },
    /// The payload mixed [`BIT_ERASED`] (`0x0000`) markers with normal
    /// `0x007F`/`0x0081` bit words. The reference encoder writes either
    /// an all-erased payload or no erased words at all; a mid-frame mix
    /// is a corruption signal. Carries the erased/normal word counts.
    MixedErasure { erased: usize, normal: usize },
}

impl core::fmt::Display for AnnexBError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooShort(n) => {
                write!(f, "g729 Annex B frame: need ≥4 header bytes, got {n}")
            }
            Self::Truncated { needed, available } => write!(
                f,
                "g729 Annex B frame: header advertises {needed} bytes, only {available} available"
            ),
            Self::BadSync(v) => write!(
                f,
                "g729 Annex B sync word mismatch: expected 0x{SYNC_WORD:04X} or 0x{ERASED_SYNC_WORD:04X}, got 0x{v:04X}"
            ),
            Self::BadBitCount(n) => write!(
                f,
                "g729 Annex B bit-count header {n} not one of 0 / {SID_WORDS} / {ACTIVE_BITS}"
            ),
            Self::BadBitWord { index, value } => write!(
                f,
                "g729 Annex B bit-word #{index} not 0x{BIT_ZERO:04X}/0x{BIT_ONE:04X}/0x{BIT_ERASED:04X}: got 0x{value:04X}"
            ),
            Self::MixedErasure { erased, normal } => write!(
                f,
                "g729 Annex B frame mixes erasure-marker words ({erased}) with normal bits ({normal})"
            ),
        }
    }
}

impl std::error::Error for AnnexBError {}

/// Read a little-endian `Word16` at word index `i` of `bytes`.
#[inline]
fn word(bytes: &[u8], i: usize) -> u16 {
    u16::from_le_bytes([bytes[i * 2], bytes[i * 2 + 1]])
}

/// Decode one Annex B serial bit word (`0x007F`/`0x0081`) at word index
/// `i`, mapping it to a boolean. Returns the offending value on a bad
/// word so the caller can report which word index failed.
#[inline]
fn bit_at(bytes: &[u8], i: usize) -> Result<bool, u16> {
    match word(bytes, i) {
        BIT_ZERO => Ok(false),
        BIT_ONE => Ok(true),
        other => Err(other),
    }
}

/// Parse one Annex B serial frame from the head of `bytes`, returning the
/// classified [`AnnexBFrame`] and the number of bytes the frame consumed.
///
/// Unlike the base-codec [`crate::serial::parse_frame`] (which assumes a
/// fixed 164-byte frame), Annex B frames are **variable length**: the
/// caller walks the stream by feeding the remaining slice and advancing
/// by the returned byte count.
///
/// # Errors
///
/// Returns [`AnnexBError`] for a short/truncated slice, an unknown sync
/// word, an unexpected bit-count header, or an invalid bit-payload word.
pub fn parse_annex_b_frame(bytes: &[u8]) -> Result<(AnnexBFrame, usize), AnnexBError> {
    if bytes.len() < 4 {
        return Err(AnnexBError::TooShort(bytes.len()));
    }

    let sync = word(bytes, 0);
    let erased = match sync {
        SYNC_WORD => false,
        ERASED_SYNC_WORD => true,
        other => return Err(AnnexBError::BadSync(other)),
    };

    let n_bits = usize::from(word(bytes, 1));
    // Validate the bit count against the Annex B frame-size set.
    let payload_words = match n_bits {
        0 => 0,
        SID_WORDS => SID_WORDS,
        ACTIVE_BITS => ACTIVE_BITS,
        other => return Err(AnnexBError::BadBitCount(other)),
    };

    let needed = (2 + payload_words) * 2;
    if bytes.len() < needed {
        return Err(AnnexBError::Truncated {
            needed,
            available: bytes.len(),
        });
    }

    // An erased-sync frame's payload is ignored; its concealed type is
    // resolved against the preceding frame by the caller (§B.4.5).
    if erased {
        return Ok((AnnexBFrame::Erased, needed));
    }

    // A normal-sync frame whose entire (non-empty) payload is the
    // BIT_ERASED marker is the second §B.4.5 erasure shape. Classify the
    // payload words first: all-erased ⇒ Erased; a mix ⇒ corruption.
    if payload_words > 0 {
        let mut erased_words = 0usize;
        let mut normal_words = 0usize;
        for i in 0..payload_words {
            match word(bytes, i + 2) {
                BIT_ERASED => erased_words += 1,
                BIT_ZERO | BIT_ONE => normal_words += 1,
                other => {
                    return Err(AnnexBError::BadBitWord {
                        index: i,
                        value: other,
                    })
                }
            }
        }
        if erased_words > 0 && normal_words > 0 {
            return Err(AnnexBError::MixedErasure {
                erased: erased_words,
                normal: normal_words,
            });
        }
        if erased_words == payload_words {
            return Ok((AnnexBFrame::Erased, needed));
        }
    }

    let frame = match n_bits {
        0 => AnnexBFrame::Untransmitted,
        SID_WORDS => {
            let sid = unpack_sid(bytes)?;
            AnnexBFrame::Sid(sid)
        }
        ACTIVE_BITS => {
            let mut bits = Box::new([false; ACTIVE_BITS]);
            for (i, b) in bits.iter_mut().enumerate() {
                *b = bit_at(bytes, i + 2)
                    .map_err(|value| AnnexBError::BadBitWord { index: i, value })?;
            }
            AnnexBFrame::Active(bits)
        }
        // Unreachable: payload_words was validated above.
        _ => return Err(AnnexBError::BadBitCount(n_bits)),
    };

    Ok((frame, needed))
}

/// Unpack a SID frame's 15 Table-B.2 bits (word indices 2..17 of the
/// frame) into the four typed indices. The 16th transmitted word
/// (word index 17) is padding and is not consumed.
fn unpack_sid(bytes: &[u8]) -> Result<SidFrame, AnnexBError> {
    // Read the 15 meaningful bits, MSB-first per field, in Table B.2
    // order.
    let mut idx = 2usize; // first payload word
    let mut read = |n: usize| -> Result<u8, AnnexBError> {
        let mut v = 0u8;
        for _ in 0..n {
            let bit = bit_at(bytes, idx).map_err(|value| AnnexBError::BadBitWord {
                index: idx - 2,
                value,
            })?;
            v = (v << 1) | u8::from(bit);
            idx += 1;
        }
        Ok(v)
    };

    let lp0 = read(SID_LP0_BITS)?;
    let l1 = read(SID_L1_BITS)?;
    let l2 = read(SID_L2_BITS)?;
    let gain = read(SID_GAIN_BITS)?;

    Ok(SidFrame { lp0, l1, l2, gain })
}

/// Walk an entire Annex B serial bitstream, returning every frame in
/// stream order paired with the byte offset it started at.
///
/// # Errors
///
/// Surfaces the first [`AnnexBError`] from [`parse_annex_b_frame`]; on
/// success the whole buffer is consumed (a trailing partial frame is a
/// [`AnnexBError::Truncated`] / [`AnnexBError::TooShort`]).
pub fn parse_annex_b_stream(mut bytes: &[u8]) -> Result<Vec<AnnexBFrame>, AnnexBError> {
    let mut frames = Vec::new();
    while !bytes.is_empty() {
        let (frame, consumed) = parse_annex_b_frame(bytes)?;
        frames.push(frame);
        bytes = &bytes[consumed..];
    }
    Ok(frames)
}

/// §B.4.2.1 SID **energy** dequantizer: map a 5-bit gain index to its
/// reconstructed log-energy level in dB.
///
/// Per §B.4.2.1 the energy is quantized "with a 5-bit non-uniform
/// quantizer in the logarithmic domain in the range of –12 dB to 66 dB":
///
/// * a **2 dB** step between 16 dB and 66 dB,
/// * a **4 dB** step between –4 dB and 16 dB,
/// * and below –4 dB a single **8 dB** step giving a level of –12 dB.
///
/// The spec states this quantizer "is straightforward and does not need
/// the storage of a quantizer table", so the level set is fully
/// determined by the prose. The 32 levels (5 bits) partition as:
///
/// | Index    | Level (dB)            | Step |
/// |----------|-----------------------|------|
/// | 0        | −12                   |  —   |
/// | 1 … 5    | −4, 0, 4, 8, 12       | 4 dB |
/// | 6 … 31   | 16, 18, …, 66         | 2 dB |
///
/// (index 5 = 12 dB, index 6 = 16 dB; 26 levels of 2 dB from 16 to 66
/// inclusive = `(66−16)/2 + 1`, 5 levels of 4 dB from −4 to 12, plus the
/// single −12 dB floor = 32.)
///
/// The `index` is taken modulo 32 (its 5-bit field width) so an
/// out-of-range value cannot panic.
#[must_use]
pub fn dequant_sid_energy_db(index: u8) -> f32 {
    let q = i32::from(index & 0x1F);
    if q == 0 {
        // Below −4 dB: single 8 dB step to the −12 dB floor.
        -12.0
    } else if q <= 5 {
        // −4 dB to 12 dB in 4 dB steps (q = 1 → −4, q = 5 → 12).
        -4.0 + 4.0 * (q - 1) as f32
    } else {
        // 16 dB to 66 dB in 2 dB steps (q = 6 → 16, q = 31 → 66).
        16.0 + 2.0 * (q - 6) as f32
    }
}

/// The resolved (post-§B.4.5) transmission type of a decoded Annex B
/// frame, after the erasure-inheritance rule has been applied.
///
/// This differs from [`AnnexBFrame`] in that [`AnnexBFrame::Erased`] has
/// been resolved against the preceding frame: an erasure after an active
/// frame becomes [`Self::Active`]; an erasure after silence becomes
/// [`Self::Untransmitted`]. The variant carries the data the comfort-
/// noise / speech synthesis stage needs.
///
/// Carries `f32` energies, so it is [`PartialEq`] but not [`Eq`].
#[derive(Debug, Clone, PartialEq)]
pub enum ResolvedFrame {
    /// Active speech: decode through the §4.1 base chain. Carries the
    /// 80 transmitted bits.
    Active(Box<[bool; ACTIVE_BITS]>),
    /// A SID frame refreshing the comfort-noise parameters. Carries the
    /// freshly-decoded [`SidFrame`] and the §B.4.2.1 energy in dB.
    Sid { sid: SidFrame, energy_db: f32 },
    /// A non-active frame with no fresh parameters — comfort noise
    /// continues from the last [`Self::Sid`]. Carries that last SID's
    /// energy (`None` before any SID has been seen, e.g. the leading
    /// frames of a stream that opens on silence).
    Untransmitted { last_energy_db: Option<f32> },
    /// An active frame was erased: the §4.4 base-codec concealment path
    /// runs. (Distinct from [`Self::Untransmitted`], which is the
    /// silence-side concealment.)
    ErasedActive,
}

/// Stateful Annex B decoder driver implementing the §B.4.1 frame-type
/// continuity and the §B.4.5 erasure-inheritance rule.
///
/// Feed each parsed [`AnnexBFrame`] in stream order to [`Self::resolve`];
/// the driver tracks the previous resolved type and the last SID
/// parameters so it can:
///
/// * resolve [`AnnexBFrame::Erased`] per §B.4.5 — "if the preceding frame
///   was active, then the current frame is considered as active; else if
///   the preceding frame was either a SID frame or an untransmitted
///   frame, the current erased frame is considered as untransmitted";
/// * carry the last SID energy forward across [`AnnexBFrame::Untransmitted`]
///   frames (§B.4.4 — "the non-active voice signal is generated … according
///   to the last received energy and spectral shape information").
///
/// This is the routing layer; the actual §4.1 speech synthesis and the
/// §B.4.4 CNG synthesis (blocked on absent tables) plug in downstream by
/// matching on the [`ResolvedFrame`] this returns.
#[derive(Debug, Clone)]
pub struct AnnexBDecoder {
    /// The previous frame's resolved transmission class (whether it was
    /// active or silence), needed for the §B.4.5 erasure rule. `None`
    /// before the first frame.
    prev_active: Option<bool>,
    /// The last received SID frame, reused across untransmitted frames
    /// (§B.4.4). `None` before any SID has been seen.
    last_sid: Option<SidFrame>,
}

impl Default for AnnexBDecoder {
    fn default() -> Self {
        Self::new()
    }
}

impl AnnexBDecoder {
    /// A fresh driver with no prior frame and no SID parameters.
    #[must_use]
    pub fn new() -> Self {
        Self {
            prev_active: None,
            last_sid: None,
        }
    }

    /// The last received SID frame, or `None` if none has been seen yet.
    #[must_use]
    pub fn last_sid(&self) -> Option<SidFrame> {
        self.last_sid
    }

    /// Resolve one parsed [`AnnexBFrame`] into a [`ResolvedFrame`],
    /// advancing the §B.4.1 / §B.4.5 state.
    ///
    /// The returned variant tells the synthesis stage which path to run;
    /// see [`ResolvedFrame`].
    pub fn resolve(&mut self, frame: &AnnexBFrame) -> ResolvedFrame {
        match frame {
            AnnexBFrame::Active(bits) => {
                self.prev_active = Some(true);
                ResolvedFrame::Active(bits.clone())
            }
            AnnexBFrame::Sid(sid) => {
                self.last_sid = Some(*sid);
                self.prev_active = Some(false);
                ResolvedFrame::Sid {
                    sid: *sid,
                    energy_db: sid.energy_db(),
                }
            }
            AnnexBFrame::Untransmitted => {
                self.prev_active = Some(false);
                ResolvedFrame::Untransmitted {
                    last_energy_db: self.last_sid.map(|s| s.energy_db()),
                }
            }
            AnnexBFrame::Erased => {
                // §B.4.5: an erased frame inherits the preceding frame's
                // class. Active→active concealment; else→untransmitted.
                // Before any frame (prev_active None), treat as silence.
                if self.prev_active == Some(true) {
                    // Stays "active" for the inheritance chain.
                    self.prev_active = Some(true);
                    ResolvedFrame::ErasedActive
                } else {
                    self.prev_active = Some(false);
                    ResolvedFrame::Untransmitted {
                        last_energy_db: self.last_sid.map(|s| s.energy_db()),
                    }
                }
            }
        }
    }
}

/// Number of PCM samples in one G.729 frame (10 ms at 8 kHz).
pub const FRAME_SAMPLES: usize = 2 * crate::fixed_codebook::SUBFRAME_SIZE;

/// What the §B stream decoder produced for one frame.
#[derive(Debug, Clone, PartialEq)]
pub enum AnnexBOutput {
    /// Active speech, reconstructed bit-exactly through the §4.1 →
    /// §4.1.6 → §4.2 base decode chain. Carries 80 post-filtered samples.
    Speech(Box<[f32; FRAME_SAMPLES]>),
    /// A non-active (SID / untransmitted / silence-erasure) frame. The
    /// §B.4.4 comfort-noise synthesis is **blocked on absent tables**
    /// (the SID-LSP VQ subset codebooks), so this is a documented
    /// placeholder: 80 samples of silence, tagged with the comfort-noise
    /// energy that *would* drive CNG once the tables are staged.
    ComfortNoisePlaceholder { energy_db: Option<f32> },
    /// An active frame was erased; the §4.4 base-codec concealment path
    /// should run. Carries no synthesized samples here because the base
    /// chain surfaces erasure separately — the caller drives concealment.
    ErasedActivePlaceholder,
}

/// End-to-end Annex B stream decoder: routes each parsed frame to the
/// right path and drives the base §4.1 → §4.2 chain for active speech.
///
/// This threads two pieces of cross-frame state together — the §B.4.1 /
/// §B.4.5 routing ([`AnnexBDecoder`]) and the base-codec synthesis
/// ([`crate::decode_chain::FrameDecoder`]) — so an Annex B `.bit` stream
/// decodes to per-frame PCM blocks in one walk. Active frames are
/// reconstructed bit-exactly; non-active frames return a documented
/// [`AnnexBOutput::ComfortNoisePlaceholder`] because the §B.4.4 CNG
/// synthesis is blocked on absent numeric tables (see the module-level
/// docs-gap note).
#[derive(Debug, Clone)]
pub struct AnnexBStreamDecoder {
    router: AnnexBDecoder,
    base: crate::decode_chain::FrameDecoder,
}

impl Default for AnnexBStreamDecoder {
    fn default() -> Self {
        Self::new()
    }
}

impl AnnexBStreamDecoder {
    /// A fresh decoder with clause-4.3 / §B start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            router: AnnexBDecoder::new(),
            base: crate::decode_chain::FrameDecoder::new(),
        }
    }

    /// Borrow the §B.4.1/§B.4.5 routing state for inspection / tests.
    #[must_use]
    pub fn router(&self) -> &AnnexBDecoder {
        &self.router
    }

    /// Decode one parsed [`AnnexBFrame`], advancing all cross-frame state.
    ///
    /// Active frames run the base §4.1 → §4.1.6 → §4.2 chain and return
    /// [`AnnexBOutput::Speech`]; non-active frames return a comfort-noise
    /// placeholder; an erased-active frame returns
    /// [`AnnexBOutput::ErasedActivePlaceholder`].
    ///
    /// # Errors
    ///
    /// Surfaces a [`crate::decode_chain::FrameDecodeError`] only on an
    /// active frame whose 80-bit payload fails the base decode (which
    /// cannot happen on a well-formed stream).
    pub fn decode_frame(
        &mut self,
        frame: &AnnexBFrame,
    ) -> Result<AnnexBOutput, crate::decode_chain::FrameDecodeError> {
        match self.router.resolve(frame) {
            ResolvedFrame::Active(bits) => {
                let kind = crate::serial::FrameKind::Active(bits);
                let pf = self.base.decode_frame_kind_to_postfiltered(&kind)?;
                Ok(AnnexBOutput::Speech(Box::new(pf.output())))
            }
            ResolvedFrame::Sid { energy_db, .. } => Ok(AnnexBOutput::ComfortNoisePlaceholder {
                energy_db: Some(energy_db),
            }),
            ResolvedFrame::Untransmitted { last_energy_db } => {
                Ok(AnnexBOutput::ComfortNoisePlaceholder {
                    energy_db: last_energy_db,
                })
            }
            ResolvedFrame::ErasedActive => Ok(AnnexBOutput::ErasedActivePlaceholder),
        }
    }

    /// Walk an entire Annex B serial bitstream, decoding every frame in
    /// stream order.
    ///
    /// # Errors
    ///
    /// Surfaces the first [`AnnexBError`] from the framing parse (wrapped
    /// in [`crate::decode_chain::FrameDecodeError::Serial`] is *not* used
    /// — framing errors are returned via [`StreamDecodeError::Framing`]),
    /// or a base-chain decode error via [`StreamDecodeError::Decode`].
    pub fn decode_stream(
        &mut self,
        mut bytes: &[u8],
    ) -> Result<Vec<AnnexBOutput>, StreamDecodeError> {
        let mut out = Vec::new();
        while !bytes.is_empty() {
            let (frame, consumed) =
                parse_annex_b_frame(bytes).map_err(StreamDecodeError::Framing)?;
            let decoded = self
                .decode_frame(&frame)
                .map_err(StreamDecodeError::Decode)?;
            out.push(decoded);
            bytes = &bytes[consumed..];
        }
        Ok(out)
    }
}

/// Error from [`AnnexBStreamDecoder::decode_stream`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StreamDecodeError {
    /// The serial framing parse failed.
    Framing(AnnexBError),
    /// An active frame's base-codec decode failed.
    Decode(crate::decode_chain::FrameDecodeError),
}

impl core::fmt::Display for StreamDecodeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Framing(e) => write!(f, "g729 Annex B framing: {e}"),
            Self::Decode(e) => write!(f, "g729 Annex B base decode: {e}"),
        }
    }
}

impl std::error::Error for StreamDecodeError {}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic Annex B serial frame: `sync` header word,
    /// `n_bits` payload-word count header, then the supplied bit words.
    fn synth_frame(sync: u16, n_bits: u16, bits: &[bool]) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&sync.to_le_bytes());
        buf.extend_from_slice(&n_bits.to_le_bytes());
        for &b in bits {
            let w = if b { BIT_ONE } else { BIT_ZERO };
            buf.extend_from_slice(&w.to_le_bytes());
        }
        buf
    }

    #[test]
    fn sid_bits_total_15() {
        assert_eq!(SID_BITS, 15);
        assert_eq!(SID_LP0_BITS + SID_L1_BITS + SID_L2_BITS + SID_GAIN_BITS, 15);
    }

    #[test]
    fn parse_untransmitted_frame() {
        let buf = synth_frame(SYNC_WORD, 0, &[]);
        let (frame, consumed) = parse_annex_b_frame(&buf).expect("parses");
        assert_eq!(frame, AnnexBFrame::Untransmitted);
        assert_eq!(consumed, 4);
    }

    #[test]
    fn parse_active_frame() {
        let bits = [true; ACTIVE_BITS];
        let buf = synth_frame(SYNC_WORD, ACTIVE_BITS as u16, &bits);
        let (frame, consumed) = parse_annex_b_frame(&buf).expect("parses");
        match frame {
            AnnexBFrame::Active(b) => assert!(b.iter().all(|x| *x)),
            other => panic!("expected Active, got {other:?}"),
        }
        assert_eq!(consumed, (2 + ACTIVE_BITS) * 2);
    }

    #[test]
    fn parse_erased_frame() {
        // Erased frame: erased sync, 0 payload.
        let buf = synth_frame(ERASED_SYNC_WORD, 0, &[]);
        let (frame, consumed) = parse_annex_b_frame(&buf).expect("parses");
        assert_eq!(frame, AnnexBFrame::Erased);
        assert_eq!(consumed, 4);
    }

    #[test]
    fn parse_sid_field_layout() {
        // Table B.2: lp0(1) l1(5) l2(4) gain(5), MSB-first per field.
        // Choose distinctive values: lp0=1, l1=0b10101=21,
        // l2=0b1100=12, gain=0b01011=11.
        let mut bits = Vec::new();
        bits.push(true); // lp0 = 1
        bits.extend_from_slice(&[true, false, true, false, true]); // l1 = 21
        bits.extend_from_slice(&[true, true, false, false]); // l2 = 12
        bits.extend_from_slice(&[false, true, false, true, true]); // gain = 11
        bits.push(false); // 16th padding word
        assert_eq!(bits.len(), SID_WORDS);

        let buf = synth_frame(SYNC_WORD, SID_WORDS as u16, &bits);
        let (frame, consumed) = parse_annex_b_frame(&buf).expect("parses");
        match frame {
            AnnexBFrame::Sid(sid) => {
                assert_eq!(sid.lp0, 1);
                assert_eq!(sid.l1, 21);
                assert_eq!(sid.l2, 12);
                assert_eq!(sid.gain, 11);
            }
            other => panic!("expected Sid, got {other:?}"),
        }
        assert_eq!(consumed, (2 + SID_WORDS) * 2);
    }

    #[test]
    fn parse_all_erased_payload_frame() {
        // Normal sync, 16-word payload all 0x0000 — the second §B.4.5
        // erasure shape (as in tstseq6's erased SID frame).
        let mut buf = Vec::new();
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&(SID_WORDS as u16).to_le_bytes());
        for _ in 0..SID_WORDS {
            buf.extend_from_slice(&crate::serial::BIT_ERASED.to_le_bytes());
        }
        let (frame, consumed) = parse_annex_b_frame(&buf).expect("parses");
        assert_eq!(frame, AnnexBFrame::Erased);
        assert_eq!(consumed, (2 + SID_WORDS) * 2);
    }

    #[test]
    fn reject_mixed_erasure_payload() {
        // 16-word payload: 15 erased markers + 1 normal bit.
        let mut buf = Vec::new();
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&(SID_WORDS as u16).to_le_bytes());
        buf.extend_from_slice(&BIT_ONE.to_le_bytes());
        for _ in 1..SID_WORDS {
            buf.extend_from_slice(&crate::serial::BIT_ERASED.to_le_bytes());
        }
        assert_eq!(
            parse_annex_b_frame(&buf).unwrap_err(),
            AnnexBError::MixedErasure {
                erased: 15,
                normal: 1,
            },
        );
    }

    #[test]
    fn reject_bad_sync() {
        let buf = synth_frame(0x1234, 0, &[]);
        assert_eq!(
            parse_annex_b_frame(&buf).unwrap_err(),
            AnnexBError::BadSync(0x1234)
        );
    }

    #[test]
    fn reject_bad_bit_count() {
        let buf = synth_frame(SYNC_WORD, 40, &[false; 40]);
        assert_eq!(
            parse_annex_b_frame(&buf).unwrap_err(),
            AnnexBError::BadBitCount(40),
        );
    }

    #[test]
    fn reject_truncated_payload() {
        // Header claims 80 payload words but supply only 4.
        let mut buf = Vec::new();
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&(ACTIVE_BITS as u16).to_le_bytes());
        buf.extend_from_slice(&BIT_ZERO.to_le_bytes());
        buf.extend_from_slice(&BIT_ZERO.to_le_bytes());
        match parse_annex_b_frame(&buf).unwrap_err() {
            AnnexBError::Truncated { needed, available } => {
                assert_eq!(needed, (2 + ACTIVE_BITS) * 2);
                assert_eq!(available, buf.len());
            }
            other => panic!("expected Truncated, got {other:?}"),
        }
    }

    #[test]
    fn reject_too_short() {
        assert_eq!(
            parse_annex_b_frame(&[0u8; 3]).unwrap_err(),
            AnnexBError::TooShort(3)
        );
    }

    #[test]
    fn energy_dequant_breakpoints() {
        // §B.4.2.1 partition.
        assert_eq!(dequant_sid_energy_db(0), -12.0); // floor
        assert_eq!(dequant_sid_energy_db(1), -4.0); // start of 4 dB band
        assert_eq!(dequant_sid_energy_db(5), 12.0); // end of 4 dB band
        assert_eq!(dequant_sid_energy_db(6), 16.0); // start of 2 dB band
        assert_eq!(dequant_sid_energy_db(31), 66.0); // top of range
    }

    #[test]
    fn energy_dequant_monotonic() {
        // The reconstructed energy levels must strictly increase with the
        // index across the whole 5-bit range.
        let mut prev = f32::NEG_INFINITY;
        for q in 0..32u8 {
            let e = dequant_sid_energy_db(q);
            assert!(e > prev, "level for q={q} ({e}) not > prev ({prev})");
            prev = e;
        }
    }

    #[test]
    fn energy_dequant_index_masked_to_5_bits() {
        // Index is taken mod 32; q=32 aliases to q=0.
        assert_eq!(dequant_sid_energy_db(32), dequant_sid_energy_db(0));
        assert_eq!(dequant_sid_energy_db(37), dequant_sid_energy_db(5));
    }

    fn sid(gain: u8) -> SidFrame {
        SidFrame {
            lp0: 0,
            l1: 0,
            l2: 0,
            gain,
        }
    }

    #[test]
    fn resolve_erased_after_active_is_erased_active() {
        let mut dec = AnnexBDecoder::new();
        let active = AnnexBFrame::Active(Box::new([false; ACTIVE_BITS]));
        let _ = dec.resolve(&active);
        let r = dec.resolve(&AnnexBFrame::Erased);
        assert_eq!(r, ResolvedFrame::ErasedActive);
    }

    #[test]
    fn resolve_erased_after_silence_is_untransmitted() {
        let mut dec = AnnexBDecoder::new();
        // SID then untransmitted then erased — the erased frame inherits
        // the silence class and becomes untransmitted.
        let _ = dec.resolve(&AnnexBFrame::Sid(sid(10)));
        let _ = dec.resolve(&AnnexBFrame::Untransmitted);
        let r = dec.resolve(&AnnexBFrame::Erased);
        match r {
            ResolvedFrame::Untransmitted { last_energy_db } => {
                assert_eq!(last_energy_db, Some(dequant_sid_energy_db(10)));
            }
            other => panic!("expected Untransmitted, got {other:?}"),
        }
    }

    #[test]
    fn resolve_leading_erased_is_untransmitted_no_energy() {
        // An erased frame at the very start (no prior frame) defaults to
        // silence with no carried energy.
        let mut dec = AnnexBDecoder::new();
        let r = dec.resolve(&AnnexBFrame::Erased);
        assert_eq!(
            r,
            ResolvedFrame::Untransmitted {
                last_energy_db: None
            }
        );
    }

    #[test]
    fn resolve_untransmitted_carries_last_sid_energy() {
        let mut dec = AnnexBDecoder::new();
        let _ = dec.resolve(&AnnexBFrame::Sid(sid(20)));
        // Two untransmitted frames keep the SID-20 energy.
        for _ in 0..2 {
            match dec.resolve(&AnnexBFrame::Untransmitted) {
                ResolvedFrame::Untransmitted { last_energy_db } => {
                    assert_eq!(last_energy_db, Some(dequant_sid_energy_db(20)));
                }
                other => panic!("expected Untransmitted, got {other:?}"),
            }
        }
        // A new SID updates the carried energy.
        let _ = dec.resolve(&AnnexBFrame::Sid(sid(6)));
        match dec.resolve(&AnnexBFrame::Untransmitted) {
            ResolvedFrame::Untransmitted { last_energy_db } => {
                assert_eq!(last_energy_db, Some(dequant_sid_energy_db(6)));
            }
            other => panic!("expected Untransmitted, got {other:?}"),
        }
        assert_eq!(dec.last_sid(), Some(sid(6)));
    }

    #[test]
    fn resolve_active_then_sid_then_active_chain() {
        let mut dec = AnnexBDecoder::new();
        let active = AnnexBFrame::Active(Box::new([true; ACTIVE_BITS]));
        assert!(matches!(dec.resolve(&active), ResolvedFrame::Active(_)));
        assert!(matches!(
            dec.resolve(&AnnexBFrame::Sid(sid(8))),
            ResolvedFrame::Sid { .. }
        ));
        // Active after silence is still active; a following erasure is
        // ErasedActive again.
        assert!(matches!(dec.resolve(&active), ResolvedFrame::Active(_)));
        assert_eq!(
            dec.resolve(&AnnexBFrame::Erased),
            ResolvedFrame::ErasedActive
        );
    }

    #[test]
    fn stream_decoder_routes_silence_to_placeholder() {
        let mut dec = AnnexBStreamDecoder::new();
        // SID then untransmitted → comfort-noise placeholders carrying
        // the SID energy.
        let out = dec
            .decode_frame(&AnnexBFrame::Sid(sid(12)))
            .expect("sid decodes");
        assert_eq!(
            out,
            AnnexBOutput::ComfortNoisePlaceholder {
                energy_db: Some(dequant_sid_energy_db(12)),
            }
        );
        let out = dec
            .decode_frame(&AnnexBFrame::Untransmitted)
            .expect("untx decodes");
        assert_eq!(
            out,
            AnnexBOutput::ComfortNoisePlaceholder {
                energy_db: Some(dequant_sid_energy_db(12)),
            }
        );
    }

    #[test]
    fn stream_decoder_erased_active_placeholder() {
        let mut dec = AnnexBStreamDecoder::new();
        let active = AnnexBFrame::Active(Box::new([false; ACTIVE_BITS]));
        let out = dec.decode_frame(&active).expect("active decodes");
        assert!(matches!(out, AnnexBOutput::Speech(_)));
        let out = dec
            .decode_frame(&AnnexBFrame::Erased)
            .expect("erased decodes");
        assert_eq!(out, AnnexBOutput::ErasedActivePlaceholder);
    }

    #[test]
    fn stream_decoder_active_produces_80_samples() {
        let mut dec = AnnexBStreamDecoder::new();
        let active = AnnexBFrame::Active(Box::new([false; ACTIVE_BITS]));
        match dec.decode_frame(&active).expect("decodes") {
            AnnexBOutput::Speech(s) => {
                assert_eq!(s.len(), FRAME_SAMPLES);
                assert!(s.iter().all(|x| x.is_finite()));
            }
            other => panic!("expected Speech, got {other:?}"),
        }
    }

    #[test]
    fn stream_walk_mixed_types() {
        let mut buf = Vec::new();
        buf.extend(synth_frame(
            SYNC_WORD,
            ACTIVE_BITS as u16,
            &[true; ACTIVE_BITS],
        ));
        buf.extend(synth_frame(
            SYNC_WORD,
            SID_WORDS as u16,
            &[false; SID_WORDS],
        ));
        buf.extend(synth_frame(SYNC_WORD, 0, &[]));
        buf.extend(synth_frame(ERASED_SYNC_WORD, 0, &[]));
        let frames = parse_annex_b_stream(&buf).expect("walks");
        assert_eq!(frames.len(), 4);
        assert!(matches!(frames[0], AnnexBFrame::Active(_)));
        assert!(matches!(frames[1], AnnexBFrame::Sid(_)));
        assert_eq!(frames[2], AnnexBFrame::Untransmitted);
        assert_eq!(frames[3], AnnexBFrame::Erased);
    }
}
