//! Robustness: hostile / random input must never panic — the registry
//! decoder over arbitrary 10-octet payloads (every Table-8 bit pattern
//! is a decodable frame, so garbage drives the whole §4.1 + §4.2 chain
//! through unstable filters, saturation and the overflow-rescale
//! path), arbitrary packet lengths, the registry encoder over
//! full-scale noise, and the serial / Annex B stream parsers and
//! decoder over random bytes.

use oxideav_core::{AudioFrame, Decoder as _, Encoder as _, Frame, Packet, TimeBase};
use oxideav_g729::annex_b::{AnnexBStreamDecoder, parse_annex_b_stream};
use oxideav_g729::codec::{G729Decoder, G729Encoder, PACKED_FRAME_BYTES};
use oxideav_g729::serial;

/// Deterministic xorshift generator (no external dependency).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn bytes(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| (self.next() >> 24) as u8).collect()
    }
}

fn drain(dec: &mut G729Decoder) -> usize {
    let mut n = 0;
    while let Ok(Frame::Audio(a)) = dec.receive_frame() {
        n += a.samples as usize;
    }
    n
}

#[test]
fn registry_decoder_survives_random_frames() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut dec = G729Decoder::new();
    let tb = TimeBase::new(1, 8000);
    let mut samples = 0usize;
    // Random payloads, biased patterns, and constant fills.
    for round in 0..6000usize {
        let payload = match round % 6 {
            0 => vec![0xFF; PACKED_FRAME_BYTES],
            1 => vec![0x00; PACKED_FRAME_BYTES],
            2 => {
                // Maximum gains / shortest delays held for many frames:
                // drives the synthesis into the overflow-rescale path.
                let mut p = rng.bytes(PACKED_FRAME_BYTES);
                p[4] |= 0x0F;
                p[5] = 0xFF;
                p[8] |= 0x7F;
                p
            }
            _ => rng.bytes(PACKED_FRAME_BYTES),
        };
        let _ = dec.send_packet(&Packet::new(0, tb, payload));
        samples += drain(&mut dec);
    }
    assert!(samples >= 80 * 5000, "decoded only {samples} samples");

    // Arbitrary packet lengths (empty, short, multi-frame, ragged).
    for len in [0usize, 1, 2, 9, 11, 19, 20, 21, 30, 100, 1001] {
        let _ = dec.send_packet(&Packet::new(0, tb, rng.bytes(len)));
        let _ = drain(&mut dec);
    }
    let _ = dec.flush();
    dec.reset().expect("reset");
}

#[test]
fn registry_encoder_survives_hostile_pcm() {
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let mut enc = G729Encoder::new();
    let mut dec = G729Decoder::new();
    let tb = TimeBase::new(1, 8000);
    for round in 0..300usize {
        let pcm: Vec<i16> = (0..80)
            .map(|n| match round % 5 {
                0 => i16::MAX,
                1 => i16::MIN,
                2 => {
                    if n % 2 == 0 {
                        i16::MAX
                    } else {
                        i16::MIN
                    }
                }
                3 => 0,
                _ => (rng.next() >> 16) as i16,
            })
            .collect();
        let bytes: Vec<u8> = pcm.iter().flat_map(|s| s.to_le_bytes()).collect();
        let frame = Frame::Audio(AudioFrame {
            samples: 80,
            pts: None,
            data: vec![bytes],
        });
        if enc.send_frame(&frame).is_ok() {
            while let Ok(pkt) = enc.receive_packet() {
                let _ = dec.send_packet(&Packet::new(0, tb, pkt.data.clone()));
                let _ = drain(&mut dec);
            }
        }
    }
    let _ = enc.flush();
}

#[test]
fn serial_and_annex_b_parsers_survive_random_bytes() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    for len in [0usize, 1, 3, 4, 163, 164, 165, 328, 1000, 5000] {
        for _ in 0..20 {
            let bytes = rng.bytes(len);
            let _ = serial::parse_frame(&bytes);
            let _ = parse_annex_b_stream(&bytes);
            let mut dec = AnnexBStreamDecoder::new();
            let _ = dec.decode_stream(&bytes);
        }
    }
    // Well-formed Annex B framing with random payload words: sync word
    // + length ∈ {0, 16, 80} + that many random bit words.
    let mut stream = Vec::new();
    for k in 0..400usize {
        let n_bits: u16 = [80, 16, 0, 80, 80][k % 5];
        stream.extend_from_slice(&0x6B21u16.to_le_bytes());
        stream.extend_from_slice(&n_bits.to_le_bytes());
        for _ in 0..n_bits {
            let word: u16 = if rng.next() & 1 == 0 { 0x007F } else { 0x0081 };
            stream.extend_from_slice(&word.to_le_bytes());
        }
    }
    let mut dec = AnnexBStreamDecoder::new();
    let _ = dec.decode_stream(&stream);
}
