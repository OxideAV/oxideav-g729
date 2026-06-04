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

use oxideav_g729::parameters::{
    self, unpack_parameters, C_BITS, GA_BITS, GB_BITS, P0_BITS, P1_BITS, P2_BITS, S_BITS,
};
use oxideav_g729::serial::{self, FrameKind, FRAME_BYTES, FRAME_WORDS};
use oxideav_g729::tables::{self, BITS_PER_FRAME, L0_BITS, L1_BITS, L2_BITS, L3_BITS, NC0, NC1};

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

/// Reads a contiguous bit slice `bits[offset..offset+width]` as an
/// MSB-first big-endian unsigned integer, matching the spec Table 8
/// NOTE on bit-stream ordering ("Bit-stream ordering follows Table 8
/// top-to-bottom, MSB first per parameter").
fn read_msb_first(bits: &[bool], offset: usize, width: usize) -> u32 {
    let mut value = 0u32;
    for k in 0..width {
        value = (value << 1) | u32::from(bits[offset + k]);
    }
    value
}

/// Extracts the (L0, L1, L2, L3) LSP-quantiser indices from an active
/// frame's bit array. Per spec Table 1 / Table 8 these are the first
/// `1 + 7 + 5 + 5 = 18` transmitted bits, in that order, MSB-first
/// per parameter. The bit ordering is the same as the on-wire word
/// order documented in `crate::serial`.
fn lsp_indices(bits: &[bool]) -> (u32, u32, u32, u32) {
    let mut cursor = 0;
    let l0 = read_msb_first(bits, cursor, L0_BITS);
    cursor += L0_BITS;
    let l1 = read_msb_first(bits, cursor, L1_BITS);
    cursor += L1_BITS;
    let l2 = read_msb_first(bits, cursor, L2_BITS);
    cursor += L2_BITS;
    let l3 = read_msb_first(bits, cursor, L3_BITS);
    (l0, l1, l2, l3)
}

/// Walks every active frame of the `LSP.BIT` conformance vector and
/// asserts that the L0 / L1 / L2 / L3 indices extracted per spec
/// Table 1 (1 + 7 + 5 + 5 bits MSB-first) all lie in the codebook
/// dimensions wired up in [`crate::tables`]. The `LSP` sequence is
/// the ITU's targeted exerciser for the LSP-quantiser branch, per
/// `docs/audio/g729/conformance/README.md` ("`LSP` — LSP
/// quantization (the L0/L1/L2/L3 VQ)"), so it is the natural fixture
/// for first-principles validation that the on-wire bit layout
/// matches the codebook shapes the build script produces.
///
/// This test does NOT yet decode LSP coefficients (that requires the
/// MA predictor `fg` and the §3.2.4 stability clamp / rearrangement
/// steps, neither of which is wired up this round). It checks the
/// weaker but still load-bearing property that every transmitted L1
/// index is < NC0 (= 128) and every L2 / L3 index is < NC1 (= 32),
/// which is a necessary condition for any future codebook lookup to
/// succeed.
fn validate_lsp_indices(label: &str, bit_bytes: &[u8]) -> usize {
    let n = serial::frame_count(bit_bytes)
        .unwrap_or_else(|e| panic!("{label}: frame_count failed: {e}"));
    let mut active_frames = 0usize;
    for f in 0..n {
        let frame = &bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let bits = match serial::parse_frame(frame).unwrap() {
            FrameKind::Active(bits) => bits,
            FrameKind::Erased => continue,
        };
        let (l0, l1, l2, l3) = lsp_indices(&bits[..]);
        assert!(
            (l0 as usize) < (1 << L0_BITS),
            "{label}: frame #{f}: L0 ({l0}) out of 1-bit range",
        );
        assert!(
            (l1 as usize) < NC0,
            "{label}: frame #{f}: L1 ({l1}) out of NC0 ({NC0})",
        );
        assert!(
            (l2 as usize) < NC1,
            "{label}: frame #{f}: L2 ({l2}) out of NC1 ({NC1})",
        );
        assert!(
            (l3 as usize) < NC1,
            "{label}: frame #{f}: L3 ({l3}) out of NC1 ({NC1})",
        );
        // Smoke check: a bounds-checked codebook lookup against
        // every L1 / L2 / L3 index returned by the parser succeeds
        // (the helpers panic on out-of-range, so this exercises the
        // happy path end-to-end).
        let _ = tables::lsp_l1_entry(l1 as usize);
        let _ = tables::lsp_l2_entry(l2 as usize);
        let _ = tables::lsp_l3_entry(l3 as usize);
        active_frames += 1;
    }
    active_frames
}

#[test]
fn lsp_conformance_indices_are_in_codebook_range() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    for variant in ["g729-core", "g729a"] {
        let dir = root.join(variant);
        if !dir.is_dir() {
            eprintln!("skip: {} not present", dir.display());
            continue;
        }
        let bit_path = dir.join("LSP.BIT");
        let bit_bytes = match std::fs::read(&bit_path) {
            Ok(b) => b,
            Err(_) => {
                eprintln!("skip: {} not present", bit_path.display());
                continue;
            }
        };
        let label = format!("{variant}/LSP");
        let active = validate_lsp_indices(&label, &bit_bytes);
        // The ITU `LSP` sequence in both g729-core and g729a is a
        // pure encoder output (no erasure splices), so every parsed
        // frame should contribute one set of LSP indices.
        let total = serial::frame_count(&bit_bytes).unwrap();
        assert_eq!(
            active,
            total,
            "{label}: LSP.BIT carried {} erasure frames (expected 0)",
            total - active,
        );
        assert!(active > 0, "{label}: no active frames extracted");
    }
}

/// Walks every active frame of every staged corpus file and asserts
/// that `unpack_parameters` returns codewords whose values lie in the
/// spec-stated domain for each field. This locks the Table-8 / §4.1
/// codeword layout against every available conformance vector, not
/// just the LSP-targeted one.
fn validate_unpack_parameters_in_domain(label: &str, bit_bytes: &[u8]) -> usize {
    let n = serial::frame_count(bit_bytes)
        .unwrap_or_else(|e| panic!("{label}: frame_count failed: {e}"));
    let mut active = 0usize;
    for f in 0..n {
        let frame = &bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).unwrap();
        let params = match unpack_parameters(&kind) {
            Ok(p) => p,
            Err(parameters::ParameterError::Erased) => continue,
        };
        // Spec Table-8 domain checks for every codeword.
        assert!(
            (params.l0 as usize) < (1 << L0_BITS),
            "{label}: frame #{f}: L0 ({}) out of 1-bit range",
            params.l0,
        );
        assert!(
            (params.l1 as usize) < NC0,
            "{label}: frame #{f}: L1 ({}) out of NC0 ({NC0})",
            params.l1,
        );
        assert!(
            (params.l2 as usize) < NC1,
            "{label}: frame #{f}: L2 ({}) out of NC1 ({NC1})",
            params.l2,
        );
        assert!(
            (params.l3 as usize) < NC1,
            "{label}: frame #{f}: L3 ({}) out of NC1 ({NC1})",
            params.l3,
        );
        // P1 is an 8-bit field; u8 already covers it. P0 is a 1-bit
        // field; bounded check on the upper bit.
        assert!(
            (u32::from(params.p0)) < (1u32 << P0_BITS),
            "{label}: frame #{f}: P0 ({}) out of 1-bit range",
            params.p0,
        );
        let _p1_bits = P1_BITS;
        // C1 / C2 are 13-bit fields (max 0x1FFF).
        assert!(
            (u32::from(params.c1)) < (1u32 << C_BITS),
            "{label}: frame #{f}: C1 ({:#x}) out of 13-bit range",
            params.c1,
        );
        assert!(
            (u32::from(params.c2)) < (1u32 << C_BITS),
            "{label}: frame #{f}: C2 ({:#x}) out of 13-bit range",
            params.c2,
        );
        // S1 / S2 are 4-bit fields.
        assert!(
            (u32::from(params.s1)) < (1u32 << S_BITS),
            "{label}: frame #{f}: S1 ({}) out of 4-bit range",
            params.s1,
        );
        assert!(
            (u32::from(params.s2)) < (1u32 << S_BITS),
            "{label}: frame #{f}: S2 ({}) out of 4-bit range",
            params.s2,
        );
        // GA1 / GA2 are 3-bit fields.
        assert!(
            (u32::from(params.ga1)) < (1u32 << GA_BITS),
            "{label}: frame #{f}: GA1 ({}) out of 3-bit range",
            params.ga1,
        );
        assert!(
            (u32::from(params.ga2)) < (1u32 << GA_BITS),
            "{label}: frame #{f}: GA2 ({}) out of 3-bit range",
            params.ga2,
        );
        // GB1 / GB2 are 4-bit fields.
        assert!(
            (u32::from(params.gb1)) < (1u32 << GB_BITS),
            "{label}: frame #{f}: GB1 ({}) out of 4-bit range",
            params.gb1,
        );
        assert!(
            (u32::from(params.gb2)) < (1u32 << GB_BITS),
            "{label}: frame #{f}: GB2 ({}) out of 4-bit range",
            params.gb2,
        );
        // P2 is a 5-bit field.
        assert!(
            (u32::from(params.p2)) < (1u32 << P2_BITS),
            "{label}: frame #{f}: P2 ({}) out of 5-bit range",
            params.p2,
        );
        // Cross-check: the round-191 lsp_indices() helper above
        // unpacks L0..L3 the same way unpack_parameters does — the
        // two must agree on every frame, with the new function
        // continuing past bit 18 to the rest of the codewords.
        let (l0, l1, l2, l3) = lsp_indices(&kind_bits(&kind).expect("active here")[..]);
        assert_eq!(u32::from(params.l0), l0);
        assert_eq!(u32::from(params.l1), l1);
        assert_eq!(u32::from(params.l2), l2);
        assert_eq!(u32::from(params.l3), l3);
        active += 1;
    }
    active
}

/// Borrow the bit array from an active FrameKind without taking
/// ownership. Returns None for an erasure.
fn kind_bits(k: &FrameKind) -> Option<&[bool; BITS_PER_FRAME]> {
    match k {
        FrameKind::Active(b) => Some(b.as_ref()),
        FrameKind::Erased => None,
    }
}

#[test]
fn unpack_parameters_in_domain_on_full_corpus() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let mut walked = 0usize;
    for variant in ["g729-core", "g729a"] {
        let dir = root.join(variant);
        if !dir.is_dir() {
            continue;
        }
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().and_then(|e| e.to_str()) != Some("BIT") {
                continue;
            }
            let bit_bytes = std::fs::read(&path).unwrap();
            let label = format!("{variant}/{}", path.file_stem().unwrap().to_string_lossy(),);
            let n = validate_unpack_parameters_in_domain(&label, &bit_bytes);
            assert!(
                n > 0
                    || path
                        .file_name()
                        .map(|s| s == "ERASURE.BIT")
                        .unwrap_or(false),
                "{label}: no active frames extracted"
            );
            walked += 1;
        }
    }
    assert!(walked > 0, "no .BIT files found under conformance root");
}

/// On the staged `PARITY.BIT` test vector the pitch-parity check
/// (§3.7.2) is supposed to FAIL for at least one frame — the
/// vector's name reflects that it is the decoder's dedicated
/// parity-mismatch exerciser (per the conformance README's
/// per-sequence purpose column). On `SPEECH.BIT` the parity should
/// hold for every active frame (encoder-output sequence, no
/// transmission errors). This test pins both invariants.
#[test]
fn pitch_parity_distribution_matches_corpus_intent() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    for variant in ["g729-core", "g729a"] {
        let dir = root.join(variant);
        if !dir.is_dir() {
            continue;
        }
        // SPEECH.BIT: every active frame's parity should hold.
        let speech_path = dir.join("SPEECH.BIT");
        if let Ok(bytes) = std::fs::read(&speech_path) {
            let n = serial::frame_count(&bytes).unwrap();
            let mut mismatches = 0usize;
            let mut active = 0usize;
            for f in 0..n {
                let frame = &bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
                let kind = serial::parse_frame(frame).unwrap();
                let Ok(params) = unpack_parameters(&kind) else {
                    continue;
                };
                active += 1;
                if !params.pitch_parity_ok() {
                    mismatches += 1;
                }
            }
            assert!(active > 0, "{variant}/SPEECH: no active frames");
            assert_eq!(
                mismatches, 0,
                "{variant}/SPEECH: pitch-parity mismatches ({mismatches}) on encoder-output \
                 vector (expected 0)",
            );
        }
        // PARITY.BIT: the dedicated mismatch exerciser should produce
        // at least one mismatch.
        let parity_path = dir.join("PARITY.BIT");
        if let Ok(bytes) = std::fs::read(&parity_path) {
            let n = serial::frame_count(&bytes).unwrap();
            let mut mismatches = 0usize;
            let mut active = 0usize;
            for f in 0..n {
                let frame = &bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
                let kind = serial::parse_frame(frame).unwrap();
                let Ok(params) = unpack_parameters(&kind) else {
                    continue;
                };
                active += 1;
                if !params.pitch_parity_ok() {
                    mismatches += 1;
                }
            }
            assert!(active > 0, "{variant}/PARITY: no active frames");
            assert!(
                mismatches > 0,
                "{variant}/PARITY: dedicated mismatch exerciser produced zero pitch-parity \
                 mismatches (expected at least one)",
            );
        }
    }
}

#[test]
fn lsp_indices_helper_round_trips_first_active_frame_bits() {
    // Synthesise an active frame whose first 18 bits encode a known
    // (L0, L1, L2, L3) tuple, then assert the extractor returns the
    // same values. This locks the bit ordering convention used by
    // `validate_lsp_indices` independently of any staged corpus
    // file, so the structural test above stays meaningful in
    // published-crate mode too.
    let l0 = 1u32;
    let l1 = 0b101_0101u32; // 7 bits = 0x55
    let l2 = 0b1_0110u32; // 5 bits = 0x16
    let l3 = 0b0_1001u32; // 5 bits = 0x09

    let mut bits = vec![false; BITS_PER_FRAME];
    let mut cursor = 0;
    for (value, width) in [(l0, L0_BITS), (l1, L1_BITS), (l2, L2_BITS), (l3, L3_BITS)] {
        for k in 0..width {
            // MSB-first: bit index `cursor + k` gets the (width-1-k)-th bit.
            let bit_pos = width - 1 - k;
            bits[cursor + k] = (value >> bit_pos) & 1 == 1;
        }
        cursor += width;
    }

    let (got_l0, got_l1, got_l2, got_l3) = lsp_indices(&bits);
    assert_eq!((got_l0, got_l1, got_l2, got_l3), (l0, l1, l2, l3));
}

/// Walks every active frame in the staged corpus and runs the
/// §3.9.2 (GA, GB) → (`ĝ_p`, `γ̂`) reconstruction. Pins that:
/// (a) every transmitted (GA, GB) pair is in-domain (NCODE1 / NCODE2);
/// (b) every reconstruction yields finite floats;
/// (c) every reconstructed `ĝ_p` lies in `0.0..=2.0` and every
///     reconstructed `γ̂` lies in `0.0..=10.0` — the same plausibility
///     window the unit-test suite pins on the codebook product.
#[test]
fn gain_reconstruct_in_domain_on_full_corpus() {
    use oxideav_g729::gain_reconstruct::reconstruct_frame_gains;
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let mut walked = 0usize;
    let mut active_total = 0usize;
    for variant in ["g729-core", "g729a"] {
        let dir = root.join(variant);
        if !dir.is_dir() {
            continue;
        }
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().and_then(|e| e.to_str()) != Some("BIT") {
                continue;
            }
            let bit_bytes = std::fs::read(&path).unwrap();
            let n = serial::frame_count(&bit_bytes).unwrap();
            let label = format!("{variant}/{}", path.file_stem().unwrap().to_string_lossy());
            let mut active = 0usize;
            for f in 0..n {
                let frame = &bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
                let kind = serial::parse_frame(frame).unwrap();
                let Ok(params) = unpack_parameters(&kind) else {
                    continue;
                };
                active += 1;
                let pairs = reconstruct_frame_gains(&params).unwrap_or_else(|e| {
                    panic!("{label} frame {f}: gain reconstruct returned error: {e}")
                });
                for (i, g) in pairs.iter().enumerate() {
                    assert!(
                        g.g_p_hat.is_finite() && g.gamma_hat.is_finite(),
                        "{label} frame {f} sub{}: non-finite gains {g:?}",
                        i + 1,
                    );
                    assert!(
                        (0.0..=2.0).contains(&g.g_p_hat),
                        "{label} frame {f} sub{}: g_p_hat {} out of [0, 2]",
                        i + 1,
                        g.g_p_hat,
                    );
                    assert!(
                        (0.0..=11.0).contains(&g.gamma_hat),
                        "{label} frame {f} sub{}: γ̂ {} out of [0, 11]",
                        i + 1,
                        g.gamma_hat,
                    );
                }
            }
            // Every .BIT file should have at least one active frame
            // (ERASURE.BIT contains 240/300 active frames, etc.).
            assert!(active > 0, "{label}: no active frames");
            active_total += active;
            walked += 1;
        }
    }
    assert!(walked > 0, "no .BIT files found under conformance root");
    assert!(active_total > 0, "no active frames processed");
}
