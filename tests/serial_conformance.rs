//! Structural conformance harness for the staged ITU-T G.729 test
//! vectors under `docs/audio/g729/conformance/`.
//!
//! Each test reads a `.bit` / `.in` / `.pst` triple and validates the
//! ITU serial bitstream geometry the docs collaborator records in
//! `docs/audio/g729/conformance/README.md` (164 bytes per frame, equal
//! frame count between `.in`, `.bit`, `.pst`, modulo encoder
//! look-ahead). For each frame: `serial::parse_frame` must succeed
//! and yield 80 well-formed bits.
//!
//! ## Why a workspace-relative path
//!
//! The conformance corpus is **not** shipped with the published crate
//! (`Cargo.toml` `include = […]` omits it; the vectors are ITU
//! electronic-attachment data that the docs collaborator stages
//! separately under `docs/audio/g729/conformance/`). When this test is
//! built from a workspace checkout (or any checkout that has the
//! `docs/audio/g729/conformance/` directory beside the crate), the
//! harness walks the staged corpus. When the path does not exist (e.g.
//! a build from `cargo install oxideav-g729`), the harness logs a
//! skip and exits clean — no panic, no false negative. This keeps the
//! crate's public-test surface healthy in both shipping modes.
//!
//! No reference decoder is run here: this is a pure-structure check
//! against the framing layout described in the conformance corpus's
//! own `README.md`.

use std::path::{Path, PathBuf};

use oxideav_g729::serial::{self, FrameKind, FRAME_BYTES, FRAME_WORDS};
use oxideav_g729::tables::BITS_PER_FRAME;

/// Walks parent directories from `CARGO_MANIFEST_DIR` looking for
/// `docs/audio/g729/conformance/`. Returns `None` if not found
/// (published-crate mode).
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

/// Returns `(.in, .bit, .pst)` byte buffers for a base/Annex-A triple
/// whose `.in` file exists. (`ERASURE`, `OVERFLOW`, `PARITY` are
/// decoder-only — no `.in` — and are handled separately.)
fn read_triple(dir: &Path, stem: &str) -> Option<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let in_path = dir.join(format!("{stem}.IN"));
    let bit_path = dir.join(format!("{stem}.BIT"));
    let pst_path_lower = dir.join(format!("{stem}.pst"));
    let pst_path_upper = dir.join(format!("{stem}.PST"));
    let pst_path = if pst_path_lower.is_file() {
        pst_path_lower
    } else {
        pst_path_upper
    };
    let in_bytes = std::fs::read(&in_path).ok()?;
    let bit_bytes = std::fs::read(&bit_path).ok()?;
    let pst_bytes = std::fs::read(&pst_path).ok()?;
    Some((in_bytes, bit_bytes, pst_bytes))
}

/// Reads a decoder-only `(.bit, .pst)` pair for `ERASURE` / `OVERFLOW`
/// / `PARITY`.
fn read_decoder_only(dir: &Path, stem: &str) -> Option<(Vec<u8>, Vec<u8>)> {
    let bit_bytes = std::fs::read(dir.join(format!("{stem}.BIT"))).ok()?;
    let pst_path_lower = dir.join(format!("{stem}.pst"));
    let pst_path_upper = dir.join(format!("{stem}.PST"));
    let pst_path = if pst_path_lower.is_file() {
        pst_path_lower
    } else {
        pst_path_upper
    };
    let pst_bytes = std::fs::read(&pst_path).ok()?;
    Some((bit_bytes, pst_bytes))
}

/// Validates that `bit_bytes` is a sequence of well-formed ITU serial
/// frames. Returns `(total_frames, erased_frames)`.
fn validate_bit_stream(label: &str, bit_bytes: &[u8]) -> (usize, usize) {
    let n = serial::frame_count(bit_bytes)
        .unwrap_or_else(|e| panic!("{label}: frame_count failed: {e}"));
    assert!(
        n > 0,
        "{label}: expected at least one ITU serial frame in {} bytes",
        bit_bytes.len(),
    );
    let mut erased = 0usize;
    for f in 0..n {
        let start = f * FRAME_BYTES;
        let end = start + FRAME_BYTES;
        let frame = &bit_bytes[start..end];
        match serial::parse_frame(frame) {
            Ok(FrameKind::Active(bits)) => assert_eq!(bits.len(), BITS_PER_FRAME),
            Ok(FrameKind::Erased) => erased += 1,
            Err(e) => panic!("{label}: frame #{f} parse failed: {e}"),
        }
    }
    (n, erased)
}

/// Cross-checks an `.in` / `.bit` / `.pst` triple per the
/// `conformance/README.md` self-consistency table:
///
/// * `.in` and `.pst` are 160-byte-per-frame PCM (Word16 LE at 8 kHz,
///   80 samples × 2 bytes); `.in` may exceed `.pst` by up to one
///   subframe's worth of look-ahead (40 samples = 80 bytes).
/// * `.bit` is exactly `164 * frames` bytes.
fn cross_check_triple(label: &str, in_bytes: &[u8], bit_bytes: &[u8], pst_bytes: &[u8]) {
    let (bit_frames, erased) = validate_bit_stream(label, bit_bytes);
    // Active/encoder-output sequences should not carry erasure
    // sentinels — those only appear in the ERASURE.BIT decoder-only
    // sequence. Catch a future drift loudly.
    assert_eq!(
        erased, 0,
        "{label}: unexpected erasure-sentinel frames in encoder-output bit stream ({erased})",
    );

    assert_eq!(
        pst_bytes.len() % 160,
        0,
        "{label}: .pst length {} not a multiple of 160 bytes",
        pst_bytes.len(),
    );
    let pst_frames = pst_bytes.len() / 160;
    assert_eq!(
        pst_frames, bit_frames,
        "{label}: .pst frame count ({pst_frames}) != .bit frame count ({bit_frames})",
    );

    // `.in` PCM has 80 samples per frame; the encoder's 5 ms look-ahead
    // (clause 2.3) can produce up to 40 extra trailing samples in `.in`
    // relative to `.pst`. The conformance README's table shows up to
    // 64-byte excess for SPEECH, which is 32 samples — within one
    // subframe.
    let in_min = pst_bytes.len();
    let in_max = pst_bytes.len() + 80; // 40 samples * 2 bytes/sample
    assert!(
        in_bytes.len() >= in_min && in_bytes.len() <= in_max,
        "{label}: .in length {} outside expected window [{in_min}, {in_max}] given .pst {}",
        in_bytes.len(),
        pst_bytes.len(),
    );
}

#[test]
fn corpus_path_is_present_or_skipped() {
    // Always passes — this test exists so the harness's skip-vs-run
    // mode is visible in the test report.
    match conformance_root() {
        Some(root) => eprintln!("conformance corpus found at {}", root.display()),
        None => eprintln!(
            "conformance corpus not present at any ancestor of CARGO_MANIFEST_DIR; \
             structural harness will be skipped (published-crate build mode)"
        ),
    }
}

#[test]
fn base_codec_full_triples_pass_structural_checks() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let dir = root.join("g729-core");
    if !dir.is_dir() {
        eprintln!("skip: {} not present", dir.display());
        return;
    }
    for stem in ["ALGTHM", "FIXED", "LSP", "PITCH", "SPEECH", "TAME"] {
        let label = format!("g729-core/{stem}");
        let (in_b, bit_b, pst_b) =
            read_triple(&dir, stem).unwrap_or_else(|| panic!("{label}: triple missing"));
        cross_check_triple(&label, &in_b, &bit_b, &pst_b);
    }
}

#[test]
fn base_codec_decoder_only_sequences_pass_structural_checks() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let dir = root.join("g729-core");
    if !dir.is_dir() {
        eprintln!("skip: {} not present", dir.display());
        return;
    }
    for stem in ["ERASURE", "OVERFLOW", "PARITY"] {
        let label = format!("g729-core/{stem}");
        let (bit_b, pst_b) = read_decoder_only(&dir, stem)
            .unwrap_or_else(|| panic!("{label}: decoder-only pair missing"));
        let (bit_frames, erased) = validate_bit_stream(&label, &bit_b);
        // ERASURE.BIT is dominated by erasure-sentinel frames (it's
        // the dedicated concealment-coverage sequence). OVERFLOW
        // contains a single spliced erasure to exercise that branch
        // alongside its main overflow-handling vectors. PARITY has
        // no erasures. The exact counts (60/1/0 of 300/384/300) are
        // pinned so a future drift in the staged corpus surfaces.
        let expected_erased = match stem {
            "ERASURE" => 60,
            "OVERFLOW" => 1,
            "PARITY" => 0,
            _ => unreachable!(),
        };
        assert_eq!(
            erased, expected_erased,
            "{label}: erasure-sentinel frame count drifted ({erased} vs expected {expected_erased})",
        );
        assert_eq!(
            pst_b.len() % 160,
            0,
            "{label}: .pst length {} not a multiple of 160 bytes",
            pst_b.len(),
        );
        assert_eq!(
            pst_b.len() / 160,
            bit_frames,
            "{label}: .pst frame count != .bit frame count",
        );
    }
}

#[test]
fn annex_a_full_triples_pass_structural_checks() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let dir = root.join("g729a");
    if !dir.is_dir() {
        eprintln!("skip: {} not present", dir.display());
        return;
    }
    // Annex A shares the base codec's eight named sequences plus an
    // extra `TEST` triple per the conformance/README.md.
    for stem in ["ALGTHM", "FIXED", "LSP", "PITCH", "SPEECH", "TAME", "TEST"] {
        let label = format!("g729a/{stem}");
        let (in_b, bit_b, pst_b) =
            read_triple(&dir, stem).unwrap_or_else(|| panic!("{label}: triple missing"));
        cross_check_triple(&label, &in_b, &bit_b, &pst_b);
    }
}

#[test]
fn annex_a_decoder_only_sequences_pass_structural_checks() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let dir = root.join("g729a");
    if !dir.is_dir() {
        eprintln!("skip: {} not present", dir.display());
        return;
    }
    for stem in ["ERASURE", "OVERFLOW", "PARITY"] {
        let label = format!("g729a/{stem}");
        let (bit_b, pst_b) = read_decoder_only(&dir, stem)
            .unwrap_or_else(|| panic!("{label}: decoder-only pair missing"));
        let (bit_frames, erased) = validate_bit_stream(&label, &bit_b);
        // Annex A decoder-only erasure counts match the base codec
        // set exactly (same encoder fixture inputs, same per-frame
        // splice points).
        let expected_erased = match stem {
            "ERASURE" => 60,
            "OVERFLOW" => 1,
            "PARITY" => 0,
            _ => unreachable!(),
        };
        assert_eq!(
            erased, expected_erased,
            "{label}: erasure-sentinel frame count drifted ({erased} vs expected {expected_erased})",
        );
        assert_eq!(
            pst_b.len() % 160,
            0,
            "{label}: .pst length {} not a multiple of 160 bytes",
            pst_b.len(),
        );
        assert_eq!(
            pst_b.len() / 160,
            bit_frames,
            "{label}: .pst frame count != .bit frame count",
        );
    }
}

#[test]
fn frame_geometry_constants_match_readme_table() {
    // Tie the parser constants to the conformance/README.md
    // self-consistency table (10 ms frame = 80 samples = 160 bytes
    // PCM; 2-word header + 80 bit-words = 164 bytes per .bit frame).
    assert_eq!(FRAME_WORDS, 82);
    assert_eq!(FRAME_BYTES, 164);
    assert_eq!(BITS_PER_FRAME, 80);
    assert_eq!(BITS_PER_FRAME * 2, 160); // PCM frame bytes
}
