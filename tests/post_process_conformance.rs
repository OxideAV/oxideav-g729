//! Conformance harness for the §4.2.5 output high-pass + ×2 upscaling
//! stage (`oxideav_g729::post_process`) against the staged ITU-T G.729
//! test vectors under `docs/audio/g729/conformance/`.
//!
//! The four front stages of the §4.2 cascade (long-/short-term
//! postfilter, tilt compensation, adaptive gain control) are not yet
//! wired, so this harness cannot do a PCM bit-exact comparison against
//! the reference output `.PST` / `.OUT` sequences (those are the output
//! of the *full* cascade). What it validates, over every active frame of
//! the base + Annex-A `.BIT` corpus, is that running the §4.1.6
//! reconstructed speech through the eq (91) filter + ×2 upscaling keeps
//! every output sample **finite** — i.e. the 2nd-order IIR is stable on
//! real decoded-speech excursions and never diverges across a whole
//! sequence with state carried frame-to-frame.
//!
//! When the corpus directory is absent (published-crate build) every
//! test logs a skip and exits clean, mirroring the sibling harnesses.

use std::path::{Path, PathBuf};

use oxideav_g729::decode_chain::FrameDecoder;
use oxideav_g729::lp_synthesis::Synthesizer;
use oxideav_g729::post_process::OutputHighPass;
use oxideav_g729::serial::{self, FrameKind, FRAME_BYTES};

/// Walks parent directories from `CARGO_MANIFEST_DIR` looking for
/// `docs/audio/g729/conformance/`. Returns `None` if not found.
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

/// Every `.BIT` file in one corpus directory, sorted for determinism.
fn bit_files(dir: &Path) -> Vec<PathBuf> {
    let mut out: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("read_dir {}: {e}", dir.display()))
        .filter_map(|entry| {
            let p = entry.ok()?.path();
            let ext = p.extension()?.to_str()?;
            ext.eq_ignore_ascii_case("bit").then_some(p)
        })
        .collect();
    out.sort();
    out
}

/// Runs decode → synthesis → eq (91) output filter over one `.BIT`
/// stream, asserting every post-filtered sample is finite. Returns the
/// active-frame count.
fn run_post(label: &str, bit_bytes: &[u8]) -> usize {
    let n = serial::frame_count(bit_bytes)
        .unwrap_or_else(|e| panic!("{label}: frame_count failed: {e}"));
    let mut chain = FrameDecoder::new();
    let mut synth = Synthesizer::new();
    let mut hpf = OutputHighPass::new();
    let mut active = 0usize;
    for f in 0..n {
        let bytes = &bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        match serial::parse_frame(bytes) {
            Ok(FrameKind::Erased) => continue,
            Ok(FrameKind::Active(_)) => {}
            Err(e) => panic!("{label}: frame #{f} framing parse failed: {e}"),
        }
        let frame = chain
            .decode_serial_frame(bytes)
            .unwrap_or_else(|e| panic!("{label}: frame #{f} decode failed: {e}"));
        let synthesized = synth.synthesize_frame(&frame);
        let mut speech = synthesized.speech();
        hpf.filter_in_place(&mut speech);
        for (i, &s) in speech.iter().enumerate() {
            assert!(
                s.is_finite(),
                "{label}: frame #{f} output sample {i} not finite",
            );
        }
        active += 1;
    }
    active
}

#[test]
fn post_process_full_corpus_finite() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let mut walked = 0usize;
    let mut total_active = 0usize;
    for dir in ["g729-core", "g729a"] {
        for bit_path in bit_files(&root.join(dir)) {
            let label = format!(
                "{dir}/{}",
                bit_path.file_name().unwrap_or_default().to_string_lossy(),
            );
            let bytes =
                std::fs::read(&bit_path).unwrap_or_else(|e| panic!("{label}: read failed: {e}"));
            let active = run_post(&label, &bytes);
            eprintln!("{label}: {active} active frames post-filtered (all finite)");
            total_active += active;
            walked += 1;
        }
    }
    assert!(walked >= 18, "expected ≥ 18 .BIT vectors, walked {walked}");
    assert!(
        total_active > 2 * 3_750,
        "expected > 7 500 active frames post-filtered, got {total_active}",
    );
    eprintln!("post-process: {total_active} active frames across {walked} vectors");
}

/// Two output filters fed the same synthesized stream stay in lockstep —
/// all eq (91) state is owned, no hidden globals.
#[test]
fn post_process_is_deterministic() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let bytes = std::fs::read(root.join("g729-core/ALGTHM.BIT")).expect("ALGTHM.BIT staged");
    let n = serial::frame_count(&bytes).expect("well-formed stream");
    let mut chain = FrameDecoder::new();
    let mut synth = Synthesizer::new();
    let mut a = OutputHighPass::new();
    let mut b = OutputHighPass::new();
    for f in 0..n {
        let frame = &bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        if matches!(serial::parse_frame(frame), Ok(FrameKind::Erased)) {
            continue;
        }
        let decoded = chain.decode_serial_frame(frame).expect("active frame");
        let speech = synth.synthesize_frame(&decoded).speech();
        assert_eq!(
            a.filter(&speech),
            b.filter(&speech),
            "frame #{f} diverged between identical output filters",
        );
    }
}
