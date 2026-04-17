//! Integration tests for the Annex B VAD/DTX/CNG path.
//!
//! These exercise the encoder + decoder coupling at the packet boundary:
//! with Annex B enabled on the encoder, feed silence and tone through the
//! encode/decode pipeline and assert the expected emission pattern.

use oxideav_codec::Encoder;
use oxideav_core::{AudioFrame, CodecId, CodecParameters, Frame, SampleFormat};
use oxideav_g729::{ANNEX_B_ENABLE_EXTRADATA, CODEC_ID_STR, FRAME_SAMPLES, SAMPLE_RATE};

fn make_params(annex_b: bool) -> CodecParameters {
    let mut p = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
    p.sample_rate = Some(SAMPLE_RATE);
    p.channels = Some(1);
    p.sample_format = Some(SampleFormat::S16);
    if annex_b {
        p.extradata = vec![ANNEX_B_ENABLE_EXTRADATA];
    }
    p
}

fn pack_audio_frame(samples: &[i16]) -> AudioFrame {
    let mut bytes = Vec::with_capacity(samples.len() * 2);
    for &s in samples {
        bytes.extend_from_slice(&s.to_le_bytes());
    }
    AudioFrame {
        format: SampleFormat::S16,
        channels: 1,
        sample_rate: SAMPLE_RATE,
        samples: samples.len() as u32,
        pts: None,
        time_base: oxideav_core::TimeBase::new(1, SAMPLE_RATE as i64),
        data: vec![bytes],
    }
}

/// Drain all pending packets from the encoder.
fn collect_packets(enc: &mut Box<dyn Encoder>) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    while let Ok(pkt) = enc.receive_packet() {
        out.push(pkt.data.clone());
    }
    out
}

/// Test: feed 100 ms of silence followed by a 1-kHz tone; the encoder
/// must emit at least one non-10-byte frame during the silence section
/// and switch to 10-byte speech frames once the tone starts.
#[test]
fn vad_classifies_silence_vs_speech() {
    let params = make_params(true);
    let mut enc = oxideav_g729::encoder::make_encoder(&params).expect("encoder");

    let frames_silence = 10; // 100 ms
    let frames_tone = 20; // 200 ms

    // --- Silence ---
    let silence = vec![0i16; FRAME_SAMPLES * frames_silence];
    enc.send_frame(&Frame::Audio(pack_audio_frame(&silence)))
        .expect("send silence");
    let silence_pkts = collect_packets(&mut enc);
    assert_eq!(silence_pkts.len(), frames_silence);

    // Must see at least one non-speech packet (SID or NODATA).
    let non_speech = silence_pkts.iter().filter(|p| p.len() != 10).count();
    assert!(
        non_speech >= 1,
        "expected at least one SID/NODATA during silence, got packet sizes {:?}",
        silence_pkts.iter().map(|p| p.len()).collect::<Vec<_>>()
    );

    // --- 1 kHz tone ---
    let sr = SAMPLE_RATE as f32;
    let mut tone = Vec::with_capacity(FRAME_SAMPLES * frames_tone);
    for i in 0..(FRAME_SAMPLES * frames_tone) {
        let t = i as f32 / sr;
        let v = 8_000.0 * (2.0 * core::f32::consts::PI * 1000.0 * t).sin();
        tone.push(v.round() as i16);
    }
    enc.send_frame(&Frame::Audio(pack_audio_frame(&tone)))
        .expect("send tone");
    enc.flush().expect("flush");
    let tone_pkts = collect_packets(&mut enc);

    // Expect all post-hang-over frames to be full 10-byte speech frames.
    // Allow up to 3 frames of hang-over bleed-over at the head of the
    // tone section where the VAD may still be deciding.
    let tail_start = 5;
    assert!(
        tone_pkts.len() >= tail_start + 5,
        "not enough tone packets: {}",
        tone_pkts.len()
    );
    for (i, pkt) in tone_pkts.iter().enumerate().skip(tail_start) {
        assert_eq!(
            pkt.len(),
            10,
            "tone packet {i} should be full 10-byte speech, got {}",
            pkt.len()
        );
    }
}

/// Test: encode silence with Annex B on; verify the output contains SID
/// (2-byte) or NODATA (0-byte) frames; feed those to the decoder and
/// confirm it produces non-zero (comfort-noise) output.
#[test]
fn sid_frame_roundtrip() {
    let params = make_params(true);
    let mut enc = oxideav_g729::encoder::make_encoder(&params).expect("encoder");
    let mut dec = oxideav_g729::decoder::make_decoder(&params).expect("decoder");

    // 200 ms of pure silence.
    let silence = vec![0i16; FRAME_SAMPLES * 20];
    enc.send_frame(&Frame::Audio(pack_audio_frame(&silence)))
        .expect("send silence");
    enc.flush().expect("flush");
    let packets = collect_packets(&mut enc);
    assert!(!packets.is_empty());

    // Every packet must be 0, 2, or 10 bytes.
    for p in &packets {
        assert!(
            p.is_empty() || p.len() == 2 || p.len() == 10,
            "unexpected packet size {}",
            p.len()
        );
    }
    // At least one SID (2-byte) frame must appear.
    let sid_count = packets.iter().filter(|p| p.len() == 2).count();
    assert!(
        sid_count >= 1,
        "expected ≥1 SID frame, packet sizes: {:?}",
        packets.iter().map(|p| p.len()).collect::<Vec<_>>()
    );

    // Feed packets to the decoder and collect output.
    let mut total_samples = 0usize;
    let mut had_nonzero_cng = false;
    let tb = oxideav_core::TimeBase::new(1, SAMPLE_RATE as i64);
    for p in &packets {
        let pkt = oxideav_core::Packet::new(0, tb, p.clone());
        dec.send_packet(&pkt).expect("send_packet");
        let frame = dec.receive_frame().expect("receive_frame");
        let Frame::Audio(a) = frame else {
            panic!("expected audio");
        };
        total_samples += a.samples as usize;
        // Look at output samples.
        for chunk in a.data[0].chunks_exact(2) {
            let s = i16::from_le_bytes([chunk[0], chunk[1]]);
            if s.unsigned_abs() > 0 {
                had_nonzero_cng = true;
            }
        }
    }

    assert_eq!(total_samples, FRAME_SAMPLES * packets.len());
    assert!(
        had_nonzero_cng,
        "decoder CNG produced all-zero output — comfort noise did not run"
    );
}

/// Test: NODATA frame before any SID must be rejected (defensive check;
/// a NODATA without priming state is malformed per Annex B §5).
#[test]
fn nodata_without_sid_is_rejected() {
    let params = make_params(false);
    let mut dec = oxideav_g729::decoder::make_decoder(&params).expect("decoder");
    let tb = oxideav_core::TimeBase::new(1, SAMPLE_RATE as i64);
    let pkt = oxideav_core::Packet::new(0, tb, Vec::new());
    dec.send_packet(&pkt).expect("send_packet");
    assert!(dec.receive_frame().is_err());
}

/// Test: Annex B disabled -> encoder emits only 10-byte packets, never
/// SID or NODATA, even on pure silence.
#[test]
fn annex_b_off_never_emits_sid() {
    let params = make_params(false);
    let mut enc = oxideav_g729::encoder::make_encoder(&params).expect("encoder");
    let silence = vec![0i16; FRAME_SAMPLES * 10];
    enc.send_frame(&Frame::Audio(pack_audio_frame(&silence)))
        .expect("send silence");
    enc.flush().expect("flush");
    let packets = collect_packets(&mut enc);
    assert!(!packets.is_empty());
    for p in &packets {
        assert_eq!(p.len(), 10, "Annex B off: all packets must be 10 bytes");
    }
}
