//! **Registry codec surface** — the `oxideav-core` trait wiring that
//! turns the crate's clause-3 encoder ([`crate::encode_chain`]) and
//! clause-4 decoder ([`crate::decode_chain`]) into registrable
//! [`Decoder`] / [`Encoder`] implementations, plus the crate's
//! `register` hook and the direct `make_decoder` / `make_encoder`
//! factory endpoints (the workspace dual-API convention).
//!
//! ## On-wire framing
//!
//! The payload format is the raw G.729 frame: **10 octets per 10 ms
//! frame** (80 bits — the spec Table 8 bit budget), transmitted
//! parameters packed in Table-8 order, MSB-first within each octet —
//! the natural dense octet packing of the Table-8 bit sequence.
//! A packet may carry any whole number of frames back-to-back; each
//! 10-octet slice decodes through the stateful chain in order.
//!
//! The 164-byte ITU *serial* format of the conformance corpus
//! ([`crate::serial`]) is a test-harness encoding of the same 80 bits
//! (one 16-bit word per bit); [`crate::parameters::pack_bit_array`] /
//! [`crate::parameters::unpack_bit_array`] are the shared bit-order
//! layer both framings sit on, so serial `.BIT` material converts
//! losslessly to the wire format (validated in
//! `tests/codec_registry.rs`).
//!
//! ## Sample format
//!
//! 8000 Hz, mono, S16 interleaved (the clause-2 16-bit PCM interface).
//! Decoded floats are rounded-and-saturated onto the 16-bit grid; the
//! encoder consumes S16 bytes directly.

use std::collections::VecDeque;

use oxideav_core::{
    AudioFrame, CodecCapabilities, CodecId, CodecInfo, CodecParameters, Error, Frame, Packet,
    Result, RuntimeContext, SampleFormat, TimeBase,
};

use crate::decode_chain::FrameDecoder;
use crate::encode_chain::{FrameEncoder, L_FRAME};
use crate::parameters::{pack_bit_array, unpack_bit_array};
use crate::tables::BITS_PER_FRAME;

/// Registry codec id.
pub const CODEC_ID: &str = "g729";
/// The only sample rate G.729 defines (clause 1.2).
pub const SAMPLE_RATE: u32 = 8000;
/// G.729 is mono.
pub const CHANNELS: u16 = 1;
/// Octets per encoded 10 ms frame (80 bits, Table 8).
pub const PACKED_FRAME_BYTES: usize = BITS_PER_FRAME / 8;
/// PCM samples per frame (10 ms at 8 kHz).
pub const SAMPLES_PER_FRAME: usize = L_FRAME;
/// Nominal bit rate (8 kbit/s).
pub const BIT_RATE: u64 = 8000;

/// Expands one 10-octet payload frame into the Table-8 bit order
/// (MSB-first within each octet).
fn payload_to_bits(payload: &[u8]) -> [bool; BITS_PER_FRAME] {
    std::array::from_fn(|i| (payload[i / 8] >> (7 - (i % 8))) & 1 == 1)
}

/// Packs a Table-8 bit array into the 10-octet wire frame (MSB-first
/// within each octet).
fn bits_to_payload(bits: &[bool; BITS_PER_FRAME]) -> [u8; PACKED_FRAME_BYTES] {
    let mut out = [0u8; PACKED_FRAME_BYTES];
    for (i, &b) in bits.iter().enumerate() {
        if b {
            out[i / 8] |= 1 << (7 - (i % 8));
        }
    }
    out
}

/// Rounds and saturates a decoded float sample onto the 16-bit PCM
/// grid.
fn to_i16(s: f32) -> i16 {
    s.round().clamp(f32::from(i16::MIN), f32::from(i16::MAX)) as i16
}

/// Validates the audio-side parameters shared by both factories.
fn check_params(params: &CodecParameters) -> Result<()> {
    if let Some(sr) = params.sample_rate {
        if sr != SAMPLE_RATE {
            return Err(Error::unsupported(format!(
                "g729: sample rate {sr} unsupported (G.729 is {SAMPLE_RATE} Hz only)"
            )));
        }
    }
    if let Some(ch) = params.channels {
        if ch != CHANNELS {
            return Err(Error::unsupported(format!(
                "g729: {ch} channels unsupported (G.729 is mono)"
            )));
        }
    }
    Ok(())
}

// ───────────────────────── decoder ─────────────────────────

/// Registry [`Decoder`]: 10-octet G.729 frames in, S16 mono 8 kHz
/// PCM frames out (one [`AudioFrame`] per input packet, all of the
/// packet's frames concatenated).
pub struct G729Decoder {
    id: CodecId,
    chain: FrameDecoder,
    pending: VecDeque<Frame>,
    flushed: bool,
}

impl std::fmt::Debug for G729Decoder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("G729Decoder")
            .field("pending", &self.pending.len())
            .field("flushed", &self.flushed)
            .finish_non_exhaustive()
    }
}

impl G729Decoder {
    /// Fresh decoder in the clause-4.3 start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            id: CodecId::new(CODEC_ID),
            chain: FrameDecoder::new(),
            pending: VecDeque::new(),
            flushed: false,
        }
    }
}

impl Default for G729Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl oxideav_core::Decoder for G729Decoder {
    fn codec_id(&self) -> &CodecId {
        &self.id
    }

    fn send_packet(&mut self, packet: &Packet) -> Result<()> {
        if packet.data.is_empty() || packet.data.len() % PACKED_FRAME_BYTES != 0 {
            return Err(Error::invalid(format!(
                "g729: packet length {} is not a positive multiple of the {PACKED_FRAME_BYTES}-octet frame size",
                packet.data.len()
            )));
        }
        let n_frames = packet.data.len() / PACKED_FRAME_BYTES;
        let mut pcm = Vec::with_capacity(n_frames * SAMPLES_PER_FRAME * 2);
        for chunk in packet.data.chunks_exact(PACKED_FRAME_BYTES) {
            let bits = payload_to_bits(chunk);
            let params = unpack_bit_array(&bits);
            let out = self
                .chain
                .decode_parameters_to_postfiltered(&params)
                .map_err(|e| Error::invalid(format!("g729: frame decode failed: {e:?}")))?;
            for s in out.output() {
                pcm.extend_from_slice(&to_i16(s).to_le_bytes());
            }
        }
        self.pending.push_back(Frame::Audio(AudioFrame {
            samples: (n_frames * SAMPLES_PER_FRAME) as u32,
            pts: packet.pts,
            data: vec![pcm],
        }));
        Ok(())
    }

    fn receive_frame(&mut self) -> Result<Frame> {
        match self.pending.pop_front() {
            Some(f) => Ok(f),
            None if self.flushed => Err(Error::Eof),
            None => Err(Error::NeedMore),
        }
    }

    fn flush(&mut self) -> Result<()> {
        self.flushed = true;
        Ok(())
    }

    fn reset(&mut self) -> Result<()> {
        self.chain = FrameDecoder::new();
        self.pending.clear();
        self.flushed = false;
        Ok(())
    }
}

// ───────────────────────── encoder ─────────────────────────

/// Registry [`Encoder`]: S16 mono 8 kHz PCM frames in, one 10-octet
/// packet per encoded 10 ms frame out.
pub struct G729Encoder {
    params: CodecParameters,
    chain: FrameEncoder,
    /// Carry-over samples (< one frame) between `send_frame` calls.
    buf: Vec<f32>,
    pending: VecDeque<Packet>,
    /// Running sample clock for output pts (1/8000 time base).
    next_pts: i64,
    flushed: bool,
}

impl std::fmt::Debug for G729Encoder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("G729Encoder")
            .field("buffered_samples", &self.buf.len())
            .field("pending", &self.pending.len())
            .field("flushed", &self.flushed)
            .finish_non_exhaustive()
    }
}

impl G729Encoder {
    /// Fresh encoder with the clause-4.3 start-up state and canonical
    /// output parameters (8 kHz / mono / S16 / 8 kbit/s).
    #[must_use]
    pub fn new() -> Self {
        let mut params = CodecParameters::audio(CodecId::new(CODEC_ID));
        params.sample_rate = Some(SAMPLE_RATE);
        params.channels = Some(CHANNELS);
        params.sample_format = Some(SampleFormat::S16);
        params.bit_rate = Some(BIT_RATE);
        Self {
            params,
            chain: FrameEncoder::new(),
            buf: Vec::new(),
            pending: VecDeque::new(),
            next_pts: 0,
            flushed: false,
        }
    }

    /// Encodes every complete 80-sample block in `self.buf` into
    /// pending packets.
    fn drain_full_frames(&mut self) {
        let mut off = 0;
        while self.buf.len() - off >= SAMPLES_PER_FRAME {
            let mut frame = [0.0f32; SAMPLES_PER_FRAME];
            frame.copy_from_slice(&self.buf[off..off + SAMPLES_PER_FRAME]);
            off += SAMPLES_PER_FRAME;
            let encoded = self.chain.encode_frame(&frame);
            let payload = bits_to_payload(&pack_bit_array(&encoded.params));
            let pkt = Packet::new(
                0,
                TimeBase::new(1, i64::from(SAMPLE_RATE)),
                payload.to_vec(),
            )
            .with_pts(self.next_pts)
            .with_duration(SAMPLES_PER_FRAME as i64)
            .with_keyframe(true);
            self.next_pts += SAMPLES_PER_FRAME as i64;
            self.pending.push_back(pkt);
        }
        self.buf.drain(..off);
    }
}

impl Default for G729Encoder {
    fn default() -> Self {
        Self::new()
    }
}

impl oxideav_core::Encoder for G729Encoder {
    fn codec_id(&self) -> &CodecId {
        &self.params.codec_id
    }

    fn output_params(&self) -> &CodecParameters {
        &self.params
    }

    fn send_frame(&mut self, frame: &Frame) -> Result<()> {
        let Frame::Audio(audio) = frame else {
            return Err(Error::invalid("g729: encoder accepts audio frames only"));
        };
        let [plane] = audio.data.as_slice() else {
            return Err(Error::invalid(
                "g729: expected one interleaved S16 mono plane",
            ));
        };
        if plane.len() % 2 != 0 {
            return Err(Error::invalid("g729: odd S16 plane length"));
        }
        self.buf.extend(
            plane
                .chunks_exact(2)
                .map(|c| f32::from(i16::from_le_bytes([c[0], c[1]]))),
        );
        self.drain_full_frames();
        Ok(())
    }

    fn receive_packet(&mut self) -> Result<Packet> {
        match self.pending.pop_front() {
            Some(p) => Ok(p),
            None if self.flushed => Err(Error::Eof),
            None => Err(Error::NeedMore),
        }
    }

    fn flush(&mut self) -> Result<()> {
        // Zero-pad a trailing partial block to a whole 10 ms frame —
        // G.729 has no partial-frame syntax.
        if !self.buf.is_empty() {
            self.buf.resize(SAMPLES_PER_FRAME, 0.0);
            self.drain_full_frames();
        }
        self.flushed = true;
        Ok(())
    }
}

// ───────────────────────── factories + registration ─────────────────────────

/// Direct decoder factory (dual-API endpoint, also the registry
/// [`CodecInfo`] decoder factory).
///
/// # Errors
///
/// Rejects parameters naming a sample rate other than 8000 Hz or a
/// channel count other than 1.
pub fn make_decoder(params: &CodecParameters) -> Result<Box<dyn oxideav_core::Decoder>> {
    check_params(params)?;
    Ok(Box::new(G729Decoder::new()))
}

/// Direct encoder factory (dual-API endpoint, also the registry
/// [`CodecInfo`] encoder factory).
///
/// # Errors
///
/// As [`make_decoder`].
pub fn make_encoder(params: &CodecParameters) -> Result<Box<dyn oxideav_core::Encoder>> {
    check_params(params)?;
    Ok(Box::new(G729Encoder::new()))
}

/// Installs the G.729 codec into the runtime context's codec registry
/// under the id `"g729"` (decode + encode, lossy, 8 kHz mono).
pub fn register(ctx: &mut RuntimeContext) {
    let mut caps = CodecCapabilities::audio("oxideav-g729");
    caps.decode = true;
    caps.encode = true;
    caps.lossy = true;
    caps.max_sample_rate = Some(SAMPLE_RATE);
    caps.max_channels = Some(CHANNELS);
    ctx.codecs.register(
        CodecInfo::new(CodecId::new(CODEC_ID))
            .capabilities(caps)
            .decoder(make_decoder)
            .encoder(make_encoder),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxideav_core::{Decoder as _, Encoder as _};

    /// Bit/octet packing round-trips for every bit position.
    #[test]
    fn payload_bits_roundtrip() {
        for i in 0..BITS_PER_FRAME {
            let mut bits = [false; BITS_PER_FRAME];
            bits[i] = true;
            let payload = bits_to_payload(&bits);
            assert_eq!(payload_to_bits(&payload), bits, "bit {i}");
            assert_eq!(payload.iter().map(|b| b.count_ones()).sum::<u32>(), 1);
        }
    }

    /// Encode → packet → decode through the registry surface: a tone
    /// encodes into 10-octet packets that the decoder accepts and
    /// reconstructs with real energy.
    #[test]
    fn encoder_decoder_loop_via_traits() {
        let mut enc = G729Encoder::new();
        let mut dec = G729Decoder::new();

        // 0.5 s of a voiced-ish tone as one big S16 frame.
        let n = 4000usize;
        let pcm: Vec<u8> = (0..n)
            .flat_map(|t| {
                let s = (5000.0 * (std::f32::consts::TAU * t as f32 / 57.0).sin()) as i16;
                s.to_le_bytes()
            })
            .collect();
        enc.send_frame(&Frame::Audio(AudioFrame {
            samples: n as u32,
            pts: Some(0),
            data: vec![pcm],
        }))
        .unwrap();
        enc.flush().unwrap();

        let mut n_packets = 0usize;
        let mut energy = 0.0f64;
        loop {
            match enc.receive_packet() {
                Ok(p) => {
                    assert_eq!(p.data.len(), PACKED_FRAME_BYTES);
                    assert_eq!(p.pts, Some((n_packets * SAMPLES_PER_FRAME) as i64));
                    n_packets += 1;
                    dec.send_packet(&p).unwrap();
                    let Frame::Audio(a) = dec.receive_frame().unwrap() else {
                        panic!("expected audio frame");
                    };
                    assert_eq!(a.samples as usize, SAMPLES_PER_FRAME);
                    for c in a.data[0].chunks_exact(2) {
                        let s = f64::from(i16::from_le_bytes([c[0], c[1]]));
                        energy += s * s;
                    }
                }
                Err(Error::Eof) => break,
                Err(e) => panic!("unexpected encoder error: {e}"),
            }
        }
        assert_eq!(n_packets, n / SAMPLES_PER_FRAME);
        assert!(energy > 1.0e6, "decoded energy {energy}");
        // Drained encoder and decoder both signal correctly.
        assert!(matches!(dec.receive_frame(), Err(Error::NeedMore)));
        dec.flush().unwrap();
        assert!(matches!(dec.receive_frame(), Err(Error::Eof)));
    }

    /// A packet carrying several frames decodes into one AudioFrame
    /// with the concatenated samples.
    #[test]
    fn multi_frame_packet() {
        let mut enc = G729Encoder::new();
        let pcm: Vec<u8> = (0..240usize)
            .flat_map(|t| ((800.0 * (t as f32 / 13.0).sin()) as i16).to_le_bytes())
            .collect();
        enc.send_frame(&Frame::Audio(AudioFrame {
            samples: 240,
            pts: None,
            data: vec![pcm],
        }))
        .unwrap();
        let mut wire = Vec::new();
        for _ in 0..3 {
            wire.extend_from_slice(&enc.receive_packet().unwrap().data);
        }
        let mut dec = G729Decoder::new();
        let pkt = Packet::new(0, TimeBase::new(1, 8000), wire);
        dec.send_packet(&pkt).unwrap();
        let Frame::Audio(a) = dec.receive_frame().unwrap() else {
            panic!("expected audio");
        };
        assert_eq!(a.samples, 240);
        assert_eq!(a.data[0].len(), 480);
    }

    /// Registration exposes decode+encode under "g729" and the
    /// registry factories build working instances.
    #[test]
    fn registry_roundtrip() {
        let mut ctx = RuntimeContext::new();
        register(&mut ctx);
        let id = CodecId::new(CODEC_ID);
        assert!(ctx.codecs.has_decoder(&id));
        assert!(ctx.codecs.has_encoder(&id));
        let params = CodecParameters::audio(id);
        let dec = ctx.codecs.first_decoder(&params);
        assert!(dec.is_ok(), "registry should build a decoder");
        let enc = ctx.codecs.first_encoder(&params);
        assert!(enc.is_ok(), "registry should build an encoder");
    }

    /// Parameter validation: wrong rate / channel count is rejected,
    /// unspecified is accepted.
    #[test]
    fn factory_param_validation() {
        let mut params = CodecParameters::audio(CodecId::new(CODEC_ID));
        assert!(make_decoder(&params).is_ok());
        params.sample_rate = Some(16_000);
        assert!(make_decoder(&params).is_err());
        params.sample_rate = Some(SAMPLE_RATE);
        params.channels = Some(2);
        assert!(make_encoder(&params).is_err());
        params.channels = Some(1);
        assert!(make_encoder(&params).is_ok());
    }

    /// Malformed packet lengths are rejected up front.
    #[test]
    fn bad_packet_length_rejected() {
        let mut dec = G729Decoder::new();
        for len in [1usize, 9, 11, 15] {
            let pkt = Packet::new(0, TimeBase::new(1, 8000), vec![0u8; len]);
            assert!(dec.send_packet(&pkt).is_err(), "len {len} should fail");
        }
    }
}
