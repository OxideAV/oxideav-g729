//! Encoder conformance harness: runs the round-382 clause-3
//! `FrameEncoder` over the staged ITU `.IN` input vectors (black-box
//! validators under `docs/audio/g729/conformance/g729-core/`) and
//! measures per-parameter agreement against the reference `.BIT`
//! encoder-output streams.
//!
//! The float encoder is not yet bit-exact against the fixed-point
//! reference (each analysis stage makes float-precision decisions on
//! razor-thin comparisons), so this harness pins **regression floors**
//! on the measured agreement rather than equality:
//!
//! * frame alignment is exactly 1:1 (a ±1-frame shift collapses the
//!   agreement — measured, and asserted here via the L1 rate);
//! * the §3.2.4 first-stage LSP index `L1` matches exactly on a large
//!   fraction of frames (measured 33–80% per vector);
//! * the §3.7 subframe-1 pitch delay `int(T1)` lands within ±2 samples
//!   of the reference on most frames (measured 56–91%);
//! * the §3.2.4 predictor switch `L0` matches on most frames
//!   (measured 81–95%);
//! * every one of our own emitted frames decodes cleanly through the
//!   crate's decoder (`decode_parameters_to_speech`) with finite
//!   output.
//!
//! Ratchet these floors upward as the encoder closes on the reference.
//! When the corpus directory is absent (published-crate build) every
//! test logs a skip and exits clean, mirroring the sibling harnesses.

use std::path::{Path, PathBuf};

use oxideav_g729::decode_chain::FrameDecoder;
use oxideav_g729::encode_chain::FrameEncoder;
use oxideav_g729::parameters::{unpack_parameters, Parameters};
use oxideav_g729::pitch_decode::decode_t1_from_p1;
use oxideav_g729::serial::{parse_frame, FrameKind, FRAME_BYTES};

const SAMPLES_PER_FRAME: usize = 80;

/// The clean encoder-input vectors of the base-codec corpus.
const VECTORS: [&str; 6] = ["ALGTHM", "FIXED", "LSP", "PITCH", "SPEECH", "TAME"];

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

/// Reads a 16-bit little-endian PCM `.IN` file as f32 sample values.
fn read_pcm(path: &Path) -> Vec<f32> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    bytes
        .chunks_exact(2)
        .map(|c| f32::from(i16::from_le_bytes([c[0], c[1]])))
        .collect()
}

/// Encodes a whole `.IN` vector and collects `(ours, reference)`
/// parameter pairs for every active reference frame.
fn encode_vector(root: &Path, name: &str) -> Vec<(Parameters, Parameters)> {
    let samples = read_pcm(&root.join(format!("g729-core/{name}.IN")));
    let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
    let n_frames = (samples.len() / SAMPLES_PER_FRAME).min(bit_bytes.len() / FRAME_BYTES);

    let mut enc = FrameEncoder::new();
    let mut pairs = Vec::new();
    for f in 0..n_frames {
        let mut frame = [0.0f32; SAMPLES_PER_FRAME];
        frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
        let out = enc.encode_frame(&frame);
        let rf = parse_frame(&bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES]).unwrap();
        if let FrameKind::Active(_) = rf {
            let rp = unpack_parameters(&rf).unwrap();
            pairs.push((out.params, rp));
        }
    }
    pairs
}

/// Percentage helper.
fn pct(hits: usize, total: usize) -> f64 {
    100.0 * hits as f64 / total.max(1) as f64
}

/// Whole-corpus parameter-agreement floors. Measured rates (round 382):
/// L0 81–95%, L1 33–80%, T1±2 56–91% — floors set with headroom so a
/// regression (not float jitter) trips them.
#[test]
fn parameter_agreement_floors() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    for name in VECTORS {
        let pairs = encode_vector(&root, name);
        let total = pairs.len();
        assert!(total > 0, "{name}: no active frames");

        let l0 = pairs.iter().filter(|(o, r)| o.l0 == r.l0).count();
        let l1 = pairs.iter().filter(|(o, r)| o.l1 == r.l1).count();
        let t1 = pairs
            .iter()
            .filter(|(o, r)| {
                let ot = decode_t1_from_p1(o.p1);
                let rt = decode_t1_from_p1(r.p1);
                (ot.int_t - rt.int_t).abs() <= 2
            })
            .count();

        println!(
            "{name}: frames={total} L0={:.1}% L1={:.1}% T1±2={:.1}%",
            pct(l0, total),
            pct(l1, total),
            pct(t1, total)
        );
        assert!(
            pct(l0, total) >= 75.0,
            "{name}: L0 agreement regressed to {:.1}%",
            pct(l0, total)
        );
        assert!(
            pct(l1, total) >= 25.0,
            "{name}: L1 agreement regressed to {:.1}%",
            pct(l1, total)
        );
        assert!(
            pct(t1, total) >= 45.0,
            "{name}: T1 agreement regressed to {:.1}%",
            pct(t1, total)
        );
    }
}

/// Frame alignment is 1:1: shifting our parameter stream by ±1 frame
/// against the reference must strictly reduce the L1 agreement on the
/// long SPEECH vector.
#[test]
fn frame_alignment_is_one_to_one() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let pairs = encode_vector(&root, "SPEECH");
    let n = pairs.len();
    let rate_at = |shift: i64| -> f64 {
        let mut hits = 0usize;
        let mut total = 0usize;
        for i in 0..n {
            let j = i as i64 + shift;
            if j < 0 || j as usize >= n {
                continue;
            }
            total += 1;
            if pairs[i].0.l1 == pairs[j as usize].1.l1 {
                hits += 1;
            }
        }
        pct(hits, total)
    };
    let aligned = rate_at(0);
    let minus = rate_at(-1);
    let plus = rate_at(1);
    println!("L1 agreement: shift -1 {minus:.1}%, 0 {aligned:.1}%, +1 {plus:.1}%");
    assert!(
        aligned > minus + 10.0 && aligned > plus + 10.0,
        "aligned agreement should dominate: {minus:.1} / {aligned:.1} / {plus:.1}"
    );
}

/// Every frame the encoder emits over the whole SPEECH vector decodes
/// through the crate's own decoder with finite output.
#[test]
fn own_stream_decodes_cleanly() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let samples = read_pcm(&root.join("g729-core/SPEECH.IN"));
    let n_frames = samples.len() / SAMPLES_PER_FRAME;

    let mut enc = FrameEncoder::new();
    let mut dec = FrameDecoder::new();
    let mut peak = 0.0f32;
    for f in 0..n_frames {
        let mut frame = [0.0f32; SAMPLES_PER_FRAME];
        frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
        let out = enc.encode_frame(&frame);
        let d = dec
            .decode_parameters_to_speech(&out.params)
            .unwrap_or_else(|e| panic!("frame {f}: decoder rejected our own stream: {e:?}"));
        for s in d.speech() {
            assert!(s.is_finite(), "frame {f}: non-finite decode output");
            peak = peak.max(s.abs());
        }
    }
    println!("decoded {n_frames} own-encoded frames, peak |s| = {peak}");
    assert!(peak > 0.0, "reconstruction should be non-silent");
}
