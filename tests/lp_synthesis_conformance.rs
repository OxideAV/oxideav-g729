//! Conformance harness for the §4.1.6 LP synthesis stage
//! (`oxideav_g729::lp_synthesis`) against the staged ITU-T G.729 test
//! vectors under `docs/audio/g729/conformance/`.
//!
//! The §4.2 post-processing cascade is not yet wired, so this harness
//! cannot do a PCM bit-exact comparison against the reference output
//! `.PST` / `.OUT` sequences. What it does validate, over every active
//! frame of the full base-codec + Annex-A `.BIT` corpus, is that the
//! §4.1.3 → §3.10 → §4.1.6 chain keeps every intermediate **finite**:
//!
//! * the eq (40) adaptive-codebook vector `v(n)`;
//! * the eq (75) excitation `u(n) = ĝ_p·v(n) + ĝ_c·c(n)`;
//! * the eq (77) reconstructed speech `ŝ(n)`.
//!
//! A non-finite value anywhere would mean the stateful past-excitation
//! buffer or the 10th-order synthesis-filter memory diverged — exactly
//! the failure mode this stage's recursion can exhibit if the b30
//! interpolation indices or the eq (77) sign were wrong.
//!
//! When the corpus directory is absent (published-crate build) every
//! test logs a skip and exits clean, mirroring the sibling harnesses.

use std::path::{Path, PathBuf};

use oxideav_g729::decode_chain::FrameDecoder;
use oxideav_g729::lp_synthesis::Synthesizer;
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

/// Runs decode → synthesis over one `.BIT` byte stream, asserting every
/// §4.1.6 intermediate is finite. Returns the active-frame count.
fn run_synthesis(label: &str, bit_bytes: &[u8]) -> usize {
    let n = serial::frame_count(bit_bytes)
        .unwrap_or_else(|e| panic!("{label}: frame_count failed: {e}"));
    let mut chain = FrameDecoder::new();
    let mut synth = Synthesizer::new();
    let mut active = 0usize;
    for f in 0..n {
        let bytes = &bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        // Skip §4.4 erasure sentinels (concealment not wired); the
        // synthesizer never sees them, mirroring the decode chain.
        match serial::parse_frame(bytes) {
            Ok(FrameKind::Erased) => continue,
            Ok(FrameKind::Active(_)) => {}
            Err(e) => panic!("{label}: frame #{f} framing parse failed: {e}"),
        }
        let frame = chain
            .decode_serial_frame(bytes)
            .unwrap_or_else(|e| panic!("{label}: frame #{f} decode failed: {e}"));
        let out = synth.synthesize_frame(&frame);
        for (s, sub) in out.subframes.iter().enumerate() {
            for (i, &v) in sub.adaptive.iter().enumerate() {
                assert!(
                    v.is_finite(),
                    "{label}: frame #{f} sub {s} v({i}) not finite",
                );
            }
            for (i, &u) in sub.excitation.iter().enumerate() {
                assert!(
                    u.is_finite(),
                    "{label}: frame #{f} sub {s} u({i}) not finite",
                );
            }
            for (i, &sp) in sub.speech.iter().enumerate() {
                assert!(
                    sp.is_finite(),
                    "{label}: frame #{f} sub {s} ŝ({i}) not finite",
                );
            }
        }
        active += 1;
    }
    active
}

#[test]
fn lp_synthesis_full_corpus_finite() {
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
            let active = run_synthesis(&label, &bytes);
            eprintln!("{label}: {active} active frames synthesized (all finite)");
            total_active += active;
            walked += 1;
        }
    }
    assert!(walked >= 18, "expected ≥ 18 .BIT vectors, walked {walked}");
    assert!(
        total_active > 2 * 3_750,
        "expected > 7 500 active frames synthesized, got {total_active}",
    );
    eprintln!("lp synthesis: {total_active} active frames across {walked} vectors");
}

/// Two synthesizers over the same stream produce identical
/// reconstructed speech frame-by-frame — all §4.1.6 state is owned.
#[test]
fn lp_synthesis_is_deterministic() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let bytes = std::fs::read(root.join("g729-core/ALGTHM.BIT")).expect("ALGTHM.BIT staged");
    let n = serial::frame_count(&bytes).expect("well-formed stream");
    let mut chain_a = FrameDecoder::new();
    let mut chain_b = FrameDecoder::new();
    let mut a = Synthesizer::new();
    let mut b = Synthesizer::new();
    for f in 0..n {
        let frame = &bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        if matches!(serial::parse_frame(frame), Ok(FrameKind::Erased)) {
            continue;
        }
        let fa = chain_a.decode_serial_frame(frame).expect("active frame");
        let fb = chain_b.decode_serial_frame(frame).expect("active frame");
        assert_eq!(
            a.synthesize_frame(&fa),
            b.synthesize_frame(&fb),
            "frame #{f} diverged between identical synthesizers",
        );
    }
}
