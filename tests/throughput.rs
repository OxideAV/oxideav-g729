//! Throughput report: registry decode and encode speed over the
//! clause-5 conformance vectors, in multiples of real time (one frame
//! = 10 ms). Prints the measurement; asserts only that both
//! directions run faster than real time. Skips when the corpus is
//! absent.
//!
//! ```text
//! cargo test --release --test throughput -- --nocapture
//! ```

use std::path::{Path, PathBuf};
use std::time::Instant;

use oxideav_core::{AudioFrame, Decoder as _, Encoder as _, Frame, Packet, TimeBase};
use oxideav_g729::codec::{G729Decoder, G729Encoder, PACKED_FRAME_BYTES};
use oxideav_g729::parameters::{pack_bit_array, unpack_parameters};
use oxideav_g729::serial::{FRAME_BYTES, FrameKind, parse_frame};

fn conformance_root() -> Option<PathBuf> {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut cursor: &Path = &manifest;
    loop {
        let cand = cursor.join("docs/audio/g729/conformance");
        if cand.join("README.md").is_file() {
            return Some(cand);
        }
        cursor = cursor.parent()?;
    }
}

fn wire_frames(bit: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for chunk in bit.as_chunks::<FRAME_BYTES>().0 {
        let Ok(frame) = parse_frame(chunk) else { break };
        let FrameKind::Active(_) = frame else { break };
        let bits = pack_bit_array(&unpack_parameters(&frame).unwrap());
        let mut payload = [0u8; PACKED_FRAME_BYTES];
        for (i, &b) in bits.iter().enumerate() {
            if b {
                payload[i / 8] |= 1 << (7 - (i % 8));
            }
        }
        out.extend_from_slice(&payload);
    }
    out
}

#[test]
fn throughput_report() {
    let Some(root) = conformance_root() else {
        eprintln!("g729 conformance corpus absent — skipping throughput report");
        return;
    };
    // Best-of-N timing; `G729_THROUGHPUT_REPS` overrides N (profiling).
    let reps = std::env::var("G729_THROUGHPUT_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(if cfg!(debug_assertions) { 1 } else { 5 });
    let tb = TimeBase::new(1, 8000);
    for name in ["SPEECH", "PITCH", "LSP"] {
        // Decode.
        let wire = wire_frames(&std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap());
        let n_frames = wire.len() / PACKED_FRAME_BYTES;
        let mut best_dec = f64::MAX;
        // FNV-1a over every decoded PCM byte and every encoded octet:
        // the bit-identity witness for performance work.
        let mut checksum = 0xcbf2_9ce4_8422_2325u64;
        let mut absorb = |bytes: &[u8]| {
            for &b in bytes {
                checksum = (checksum ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        for rep in 0..reps {
            let mut dec = G729Decoder::new();
            let t0 = Instant::now();
            for payload in wire.as_chunks::<PACKED_FRAME_BYTES>().0 {
                dec.send_packet(&Packet::new(0, tb, payload.to_vec()))
                    .unwrap();
                if let Ok(Frame::Audio(a)) = dec.receive_frame()
                    && rep == 0
                {
                    absorb(&a.data[0]);
                }
            }
            best_dec = best_dec.min(t0.elapsed().as_secs_f64());
        }
        let audio_s = n_frames as f64 * 0.010;

        // Encode.
        let pcm = std::fs::read(root.join(format!("g729-core/{name}.IN"))).unwrap();
        let mut best_enc = f64::MAX;
        for rep in 0..reps {
            let mut enc = G729Encoder::new();
            let t0 = Instant::now();
            for frame in pcm.as_chunks::<160>().0 {
                let f = Frame::Audio(AudioFrame {
                    samples: 80,
                    pts: None,
                    data: vec![frame.to_vec()],
                });
                enc.send_frame(&f).unwrap();
                while let Ok(pkt) = enc.receive_packet() {
                    if rep == 0 {
                        absorb(&pkt.data);
                    }
                }
            }
            best_enc = best_enc.min(t0.elapsed().as_secs_f64());
        }
        let enc_audio_s = (pcm.len() / 160) as f64 * 0.010;
        let (dx, ex) = (audio_s / best_dec, enc_audio_s / best_enc);
        eprintln!(
            "{name}: decode {dx:.0}x realtime ({n_frames} frames in {:.1} ms), encode {ex:.0}x realtime ({} frames in {:.1} ms)  [fnv {checksum:016x}]",
            best_dec * 1e3,
            pcm.len() / 160,
            best_enc * 1e3
        );
        assert!(dx > 1.0 && ex > 1.0, "{name}: slower than real time");
    }
}
