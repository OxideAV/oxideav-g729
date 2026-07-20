//! Registry-surface conformance: the 10-octet wire framing of the
//! [`oxideav_g729::codec`] module carries exactly the same 80
//! Table-8 bits as the ITU serial format of the staged conformance
//! corpus, and the registry decoder reproduces the crate's
//! fixed-point decode path (`fx` §4.1 chain + §4.2 cascade)
//! sample-for-sample on real reference bitstreams.

use std::path::{Path, PathBuf};

use oxideav_core::{CodecId, CodecParameters, Decoder as _, Frame, Packet, TimeBase};
use oxideav_g729::codec::{G729Decoder, PACKED_FRAME_BYTES, SAMPLES_PER_FRAME};
use oxideav_g729::fx::decoder::FrameDecoderFx;
use oxideav_g729::fx::postfilter::PostfilterFx;
use oxideav_g729::parameters::{pack_bit_array, unpack_parameters};
use oxideav_g729::serial::{parse_frame, FrameKind, FRAME_BYTES};

fn conformance_root() -> Option<PathBuf> {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut cursor: &Path = &manifest;
    loop {
        let cand = cursor.join("docs/audio/g729/conformance");
        if cand.is_dir() {
            return Some(cand);
        }
        cursor = cursor.parent()?;
    }
}

/// Converts one serial frame's bit array into the 10-octet wire
/// payload (MSB-first octet packing of the Table-8 bit order).
fn serial_to_wire(bits: &[bool; 80]) -> [u8; PACKED_FRAME_BYTES] {
    let mut out = [0u8; PACKED_FRAME_BYTES];
    for (i, &b) in bits.iter().enumerate() {
        if b {
            out[i / 8] |= 1 << (7 - (i % 8));
        }
    }
    out
}

/// Every active frame of every clean corpus vector, converted from the
/// ITU serial format to the wire format and decoded through the
/// registry `Decoder`, reproduces the fixed-point decode path
/// (`fx` §4.1 chain + §4.2 cascade) exactly (same parameters, same
/// state trajectory, same rounding).
#[test]
fn registry_decoder_matches_fx_chain_on_corpus() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    for name in ["ALGTHM", "FIXED", "LSP", "PITCH", "SPEECH", "TAME"] {
        let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
        let mut wire = G729Decoder::new();
        let mut fx = FrameDecoderFx::new();
        let mut pf = PostfilterFx::new();
        let mut compared = 0usize;
        for chunk in bit_bytes.chunks_exact(FRAME_BYTES) {
            let frame = parse_frame(chunk).unwrap();
            let FrameKind::Active(_) = frame else {
                // The wire format has no erasure marker; stop at the
                // first non-active frame (state stays comparable up to
                // that point).
                break;
            };
            let params = unpack_parameters(&frame).unwrap();
            let payload = serial_to_wire(&pack_bit_array(&params));

            let pkt = Packet::new(0, TimeBase::new(1, 8000), payload.to_vec());
            wire.send_packet(&pkt).unwrap();
            let Frame::Audio(a) = wire.receive_frame().unwrap() else {
                panic!("expected audio frame");
            };
            assert_eq!(a.samples as usize, SAMPLES_PER_FRAME);

            let dec = fx.decode_frame(&params);
            let int_t1 = usize::try_from(dec.sub[0].t_int.max(1)).unwrap();
            let mut expect = Vec::with_capacity(SAMPLES_PER_FRAME);
            for s in 0..2 {
                let speech: [i16; 40] = std::array::from_fn(|n| dec.speech[s * 40 + n]);
                let (out, _) = pf.process_subframe(&speech, &dec.sub[s].a_q12, int_t1);
                expect.extend_from_slice(&out);
            }
            for (n, (got, want)) in a.data[0].chunks_exact(2).zip(expect.iter()).enumerate() {
                let got = i16::from_le_bytes([got[0], got[1]]);
                assert_eq!(got, *want, "{name}: sample {n} of frame {compared}");
            }
            compared += 1;
        }
        assert!(compared > 30, "{name}: compared only {compared} frames");
        println!("{name}: registry decode matches the fx chain on {compared} frames");
    }
}

/// The registry decoder's `reset` returns it to the clause-4.3
/// start-up state: decoding the same stream before and after a reset
/// yields identical PCM.
#[test]
fn registry_decoder_reset_restores_startup_state() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let bit_bytes = std::fs::read(root.join("g729-core/ALGTHM.BIT")).unwrap();
    let mut payloads = Vec::new();
    for chunk in bit_bytes.chunks_exact(FRAME_BYTES) {
        let frame = parse_frame(chunk).unwrap();
        let FrameKind::Active(_) = frame else { break };
        let params = unpack_parameters(&frame).unwrap();
        payloads.push(serial_to_wire(&pack_bit_array(&params)));
    }
    assert!(payloads.len() > 10);

    let mut dec = G729Decoder::new();
    let run = |dec: &mut G729Decoder, payloads: &[[u8; PACKED_FRAME_BYTES]]| -> Vec<Vec<u8>> {
        payloads
            .iter()
            .map(|p| {
                let pkt = Packet::new(0, TimeBase::new(1, 8000), p.to_vec());
                dec.send_packet(&pkt).unwrap();
                let Frame::Audio(a) = dec.receive_frame().unwrap() else {
                    panic!("expected audio");
                };
                a.data.into_iter().next().unwrap()
            })
            .collect()
    };
    let first = run(&mut dec, &payloads);
    dec.reset().unwrap();
    let second = run(&mut dec, &payloads);
    assert_eq!(first, second, "reset must restore the start-up state");
}

/// Dual-API endpoints exist with the historical module paths and
/// build working instances.
#[test]
fn dual_api_endpoints() {
    let params = CodecParameters::audio(CodecId::new("g729"));
    assert!(oxideav_g729::decoder::make_decoder(&params).is_ok());
    assert!(oxideav_g729::encoder::make_encoder(&params).is_ok());
}
