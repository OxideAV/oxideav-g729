//! Conformance harness for the §4.2.2 short-term postfilter
//! (`oxideav_g729::short_term_postfilter`) against the staged ITU-T
//! G.729 test vectors under `docs/audio/g729/conformance/`.
//!
//! The §4.2.1 long-term postfilter that precedes this stage is not yet
//! wired, so this harness cannot do a PCM bit-exact comparison against
//! the reference output `.PST` / `.OUT` sequences (those are the output
//! of the *full* §4.2 cascade). What it validates, over every active
//! frame of the base + Annex-A `.BIT` corpus, is that running the
//! §4.1.6 reconstructed speech through the eq (84)/(85) short-term
//! postfilter (per-subframe `â_i`, continuous filter memory) keeps every
//! output sample **finite** — i.e. the weighted pole/zero filter pair is
//! stable on real decoded-speech excursions and never diverges across a
//! whole sequence — and that the eq (85) gain term `g_f` is positive and
//! finite for every real subframe's LP coefficients.
//!
//! When the corpus directory is absent (published-crate build) every
//! test logs a skip and exits clean, mirroring the sibling harnesses.

use std::path::{Path, PathBuf};

use oxideav_g729::decode_chain::FrameDecoder;
use oxideav_g729::lp_synthesis::Synthesizer;
use oxideav_g729::serial::{self, FrameKind, FRAME_BYTES};
use oxideav_g729::short_term_postfilter::ShortTermPostfilter;

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

/// Runs decode → synthesis → §4.2.2 short-term postfilter over one
/// `.BIT` stream, asserting every postfiltered sample and every
/// per-subframe `g_f` is finite/positive. Returns the active-frame
/// count.
fn run_postfilter(label: &str, bit_bytes: &[u8]) -> usize {
    let n = serial::frame_count(bit_bytes)
        .unwrap_or_else(|e| panic!("{label}: frame_count failed: {e}"));
    let mut chain = FrameDecoder::new();
    let mut synth = Synthesizer::new();
    let mut pf = ShortTermPostfilter::new();
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
        for (s, sub) in synthesized.subframes.iter().enumerate() {
            let a = &frame.subframes[s].lp;
            // eq (85) gain term must be finite and positive for a real
            // (stable) LP coefficient set.
            let g_f = ShortTermPostfilter::gain_term(a);
            assert!(
                g_f.is_finite() && g_f > 0.0,
                "{label}: frame #{f} sub {s} g_f = {g_f} not finite/positive",
            );
            let out = pf.filter_subframe(&sub.speech, a);
            for (i, &v) in out.iter().enumerate() {
                assert!(
                    v.is_finite(),
                    "{label}: frame #{f} sub {s} output sample {i} not finite",
                );
            }
        }
        active += 1;
    }
    active
}

#[test]
fn short_term_postfilter_full_corpus_finite() {
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
            let active = run_postfilter(&label, &bytes);
            eprintln!("{label}: {active} active frames short-term-postfiltered (all finite)");
            total_active += active;
            walked += 1;
        }
    }
    assert!(walked >= 18, "expected ≥ 18 .BIT vectors, walked {walked}");
    assert!(
        total_active > 2 * 3_750,
        "expected > 7 500 active frames postfiltered, got {total_active}",
    );
    eprintln!("short-term postfilter: {total_active} active frames across {walked} vectors");
}

/// Two postfilters fed the same synthesized stream stay in lockstep —
/// all eq (84) filter state is owned, no hidden globals.
#[test]
fn short_term_postfilter_is_deterministic() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let bytes = std::fs::read(root.join("g729-core/ALGTHM.BIT")).expect("ALGTHM.BIT staged");
    let n = serial::frame_count(&bytes).expect("well-formed stream");
    let mut chain = FrameDecoder::new();
    let mut synth = Synthesizer::new();
    let mut a = ShortTermPostfilter::new();
    let mut b = ShortTermPostfilter::new();
    for f in 0..n {
        let frame = &bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        if matches!(serial::parse_frame(frame), Ok(FrameKind::Erased)) {
            continue;
        }
        let decoded = chain.decode_serial_frame(frame).expect("active frame");
        let synthesized = synth.synthesize_frame(&decoded);
        for (s, sub) in synthesized.subframes.iter().enumerate() {
            let lp = &decoded.subframes[s].lp;
            assert_eq!(
                a.filter_subframe(&sub.speech, lp),
                b.filter_subframe(&sub.speech, lp),
                "frame #{f} sub {s} diverged between identical postfilters",
            );
        }
    }
}
