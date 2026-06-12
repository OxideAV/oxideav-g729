//! Conformance harness for the §4.1 per-frame decode parameter chain
//! (`oxideav_g729::decode_chain`) against the staged ITU-T G.729 test
//! vectors under `docs/audio/g729/conformance/`.
//!
//! The chain does not yet produce PCM (§4.1.6 synthesis and §4.2
//! post-processing are future rounds), so this harness validates the
//! **parameter-domain** outputs over every active frame of the full
//! base-codec + Annex-A `.BIT` corpus:
//!
//! * every frame decodes without error (the codeword-domain error
//!   variants must be unreachable on real streams);
//! * the §4.1.1 LSF output respects the §3.2.4 stability clamp
//!   (floor 0.005, ceiling 3.135, strictly increasing);
//! * the §4.1.3 delays land in the spec windows (`int(T1)` ∈
//!   [19, 143] on parity-clean frames, `int(T2)` ∈
//!   [t_min − 1, t_min + 10]);
//! * the §4.1.4 codevector carries the eq (45) four-pulse structure
//!   (post-eq (48) energy ≥ the 4-pulse floor when sharpening adds
//!   energy, always > 0) with finite samples;
//! * the §4.1.5 gains and the §3.9.1 prediction path stay finite;
//! * erasure sentinels (`ERASURE.BIT`) surface
//!   `FrameDecodeError::Erased` and active frames around them keep
//!   decoding;
//! * `PARITY.BIT` exercises the §4.1.2 concealment substitution on a
//!   non-zero number of frames.
//!
//! When the corpus directory is absent (published-crate build) every
//! test logs a skip and exits clean, mirroring
//! `tests/serial_conformance.rs`.

use std::path::{Path, PathBuf};

use oxideav_g729::decode_chain::{DecodedFrame, FrameDecodeError, FrameDecoder};
use oxideav_g729::lsp_reconstruct::{CLAMP_CEIL, CLAMP_FLOOR};
use oxideav_g729::serial::{self, FrameKind, FRAME_BYTES};
use oxideav_g729::tables::M;

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

/// Asserts every parameter-domain invariant on one decoded frame.
fn check_frame(label: &str, f: usize, frame: &DecodedFrame) {
    // §4.1.1 / §3.2.4 stability clamp: ω̂ strictly increasing within
    // [CLAMP_FLOOR, CLAMP_CEIL].
    for i in 0..M {
        assert!(
            frame.omega[i].is_finite(),
            "{label}: frame #{f} ω̂[{i}] not finite",
        );
        assert!(
            frame.omega[i] >= CLAMP_FLOOR - 1e-6 && frame.omega[i] <= CLAMP_CEIL + 1e-6,
            "{label}: frame #{f} ω̂[{i}] = {} outside [{CLAMP_FLOOR}, {CLAMP_CEIL}]",
            frame.omega[i],
        );
        if i > 0 {
            assert!(
                frame.omega[i] > frame.omega[i - 1],
                "{label}: frame #{f} ω̂ not strictly increasing at {i}",
            );
        }
    }

    // §4.1.3 delay windows.
    let t1 = frame.subframes[0].pitch;
    let t2 = frame.subframes[1].pitch;
    if frame.parity_ok {
        assert!(
            (19..=143).contains(&t1.int_t),
            "{label}: frame #{f} int(T1) = {} outside [19, 143]",
            t1.int_t,
        );
    }
    assert!(
        (-1..=1).contains(&t1.frac) && (-1..=1).contains(&t2.frac),
        "{label}: frame #{f} fractional parts outside {{-1, 0, 1}}",
    );
    assert!(
        (20..=134).contains(&frame.t_min),
        "{label}: frame #{f} t_min = {} outside [20, 134]",
        frame.t_min,
    );
    assert!(
        t2.int_t >= frame.t_min - 1 && t2.int_t <= frame.t_min + 10,
        "{label}: frame #{f} int(T2) = {} outside [t_min − 1, t_min + 10] (t_min = {})",
        t2.int_t,
        frame.t_min,
    );

    for (s, sub) in frame.subframes.iter().enumerate() {
        // §4.1.4 codevector: finite, non-zero energy (eq (45) places
        // four ±1 pulses; eq (48) only ever adds copies scaled by
        // β ∈ [0.2, 0.8], so a few of the four taps may overlap but
        // the energy can never reach zero).
        let mut energy = 0.0f32;
        for c in &sub.codevector {
            assert!(
                c.is_finite(),
                "{label}: frame #{f} sub {s} codevector sample not finite",
            );
            energy += c * c;
        }
        assert!(
            energy > 0.0,
            "{label}: frame #{f} sub {s} codevector energy is zero",
        );

        // Pulse positions respect the Table-7 tracks.
        for (k, p) in sub.fixed.pulses.iter().enumerate() {
            assert!(
                (p.position as usize) < 40,
                "{label}: frame #{f} sub {s} pulse {k} position out of subframe",
            );
            let residue = p.position % 5;
            let expected = match k {
                0..=2 => k as u8,
                _ => 3 + sub.fixed.jx,
            };
            assert_eq!(
                residue, expected,
                "{label}: frame #{f} sub {s} pulse {k} off its Table-7 track",
            );
            assert!(p.sign == 1 || p.sign == -1);
        }

        // §4.1.5 gains + §3.9.1 path: finite throughout, β in clamp.
        assert!(
            sub.gains.g_p_hat.is_finite() && sub.gains.gamma_hat.is_finite(),
            "{label}: frame #{f} sub {s} (ĝ_p, γ̂) not finite",
        );
        assert!(
            (0.2..=0.8).contains(&sub.beta),
            "{label}: frame #{f} sub {s} β = {} outside eq (47) clamp",
            sub.beta,
        );
        assert!(
            sub.g_c_hat.is_finite(),
            "{label}: frame #{f} sub {s} ĝ_c not finite",
        );
        assert!(
            sub.predicted.e_db.is_finite()
                && sub.predicted.e_tilde_db.is_finite()
                && sub.predicted.g_c_prime.is_finite()
                && sub.predicted.g_c_prime > 0.0,
            "{label}: frame #{f} sub {s} §3.9.1 prediction path not finite/positive",
        );

        // §4.1.1 LP coefficients finite.
        for (i, a) in sub.lp.iter().enumerate() {
            assert!(
                a.is_finite(),
                "{label}: frame #{f} sub {s} a_{} not finite",
                i + 1,
            );
        }
    }
}

/// Runs the chain over one `.BIT` byte stream. Returns
/// `(active_frames, erased_frames, parity_mismatch_frames)`.
fn run_chain(label: &str, bit_bytes: &[u8]) -> (usize, usize, usize) {
    let n = serial::frame_count(bit_bytes)
        .unwrap_or_else(|e| panic!("{label}: frame_count failed: {e}"));
    let mut chain = FrameDecoder::new();
    let mut active = 0usize;
    let mut erased = 0usize;
    let mut parity_bad = 0usize;
    for f in 0..n {
        let bytes = &bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        match serial::parse_frame(bytes) {
            Ok(FrameKind::Erased) => {
                // §4.4 concealment is not wired: the chain must reject
                // the sentinel with the dedicated variant and keep
                // decoding subsequent active frames.
                let err = chain
                    .decode_serial_frame(bytes)
                    .expect_err("erasure sentinel must not decode");
                assert_eq!(err, FrameDecodeError::Erased, "{label}: frame #{f}");
                erased += 1;
            }
            Ok(FrameKind::Active(_)) => {
                let frame = chain
                    .decode_serial_frame(bytes)
                    .unwrap_or_else(|e| panic!("{label}: frame #{f} chain failed: {e}"));
                check_frame(label, f, &frame);
                if !frame.parity_ok {
                    parity_bad += 1;
                }
                active += 1;
            }
            Err(e) => panic!("{label}: frame #{f} framing parse failed: {e}"),
        }
    }
    (active, erased, parity_bad)
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

#[test]
fn decode_chain_full_corpus_parameter_domain() {
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
            let (active, erased, parity_bad) = run_chain(&label, &bytes);
            eprintln!(
                "{label}: {active} active frames decoded ({erased} erased, {parity_bad} parity mismatches)",
            );
            let stem = bit_path
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .to_uppercase();
            match stem.as_str() {
                // ERASURE carries the §4.4 sentinels by design.
                "ERASURE" => assert!(erased > 0, "{label}: expected erasure sentinels"),
                // OVERFLOW / PARITY are the other decoder-only
                // streams; OVERFLOW empirically carries one sentinel
                // frame, so no constraint beyond "decodes cleanly".
                "OVERFLOW" | "PARITY" => {}
                _ => assert_eq!(erased, 0, "{label}: unexpected erasure sentinels"),
            }
            match stem.as_str() {
                // PARITY is the dedicated §4.1.2 exerciser; every
                // other sequence is clean encoder output.
                "PARITY" => assert!(
                    parity_bad > 0,
                    "{label}: expected §4.1.2 concealment activations",
                ),
                _ => assert_eq!(parity_bad, 0, "{label}: unexpected parity mismatches"),
            }
            total_active += active;
            walked += 1;
        }
    }
    assert!(walked >= 18, "expected ≥ 18 .BIT vectors, walked {walked}");
    // SPEECH alone is 3 750 frames per variant; the full two-variant
    // corpus must clear well past that.
    assert!(
        total_active > 2 * 3_750,
        "expected > 7 500 active frames across the corpus, got {total_active}",
    );
    eprintln!("decode chain: {total_active} active frames across {walked} vectors");
}

/// Two independent chains over the same stream produce identical
/// frame-by-frame output — the chain is deterministic and all state
/// is owned (no hidden globals).
#[test]
fn decode_chain_is_deterministic() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let bytes = std::fs::read(root.join("g729-core/ALGTHM.BIT")).expect("ALGTHM.BIT staged");
    let n = serial::frame_count(&bytes).expect("well-formed stream");
    let mut a = FrameDecoder::new();
    let mut b = FrameDecoder::new();
    for f in 0..n {
        let frame = &bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let fa = a.decode_serial_frame(frame).expect("active frame");
        let fb = b.decode_serial_frame(frame).expect("active frame");
        assert_eq!(fa, fb, "frame #{f} diverged between identical chains");
    }
}

/// The §4.1.2 concealment path on `PARITY.BIT`: on every
/// parity-mismatch frame the substituted subframe-1 delay is integer
/// (`frac = 0`) and equals the previous frame's `int(T2)`.
#[test]
fn parity_bit_concealment_substitution_is_exact() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let bytes = std::fs::read(root.join("g729-core/PARITY.BIT")).expect("PARITY.BIT staged");
    let n = serial::frame_count(&bytes).expect("well-formed stream");
    let mut chain = FrameDecoder::new();
    let mut prev_int_t2: i32 = 0; // clause 4.3 zero default
    let mut concealed = 0usize;
    for f in 0..n {
        let frame = chain
            .decode_serial_frame(&bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES])
            .unwrap_or_else(|e| panic!("PARITY frame #{f}: {e}"));
        if !frame.parity_ok {
            assert_eq!(
                frame.subframes[0].pitch.int_t, prev_int_t2,
                "PARITY frame #{f}: concealment delay must be previous int(T2)",
            );
            assert_eq!(frame.subframes[0].pitch.frac, 0);
            concealed += 1;
        }
        prev_int_t2 = frame.subframes[1].pitch.int_t;
    }
    assert!(
        concealed > 0,
        "PARITY.BIT should activate the §4.1.2 concealment at least once",
    );
    eprintln!("PARITY.BIT: {concealed} concealment activations over {n} frames");
}
