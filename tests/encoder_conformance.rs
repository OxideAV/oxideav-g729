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

/// Whole-corpus parameter-agreement floors, pinned per vector.
///
/// Measured rates after the round-390 fixed-point Q13 quantiser
/// search (on top of the round-388 residual-domain split-stage
/// metric): L0 62.5–97.1%, L1 32.9–100%, T1±2 75.0–90.5%. The largest
/// end-to-end moves vs round 388 are FIXED L0 88.3 → 90.8% and
/// L1 64.2 → 70.0%; LSP L1 gives back 0.9 (48.1 → 47.2%,
/// error-propagation reshuffling — its locked-history number rises).
/// The floors sit ~3–6 points under the measured rates so a
/// regression (not float jitter) trips them; ratchet them upward as
/// the encoder closes on the reference.
///
/// TAME is the outlier (L0 62.5%): its stationary extreme-spectra
/// input locks the reference into a 5-frame index cycle whose L0
/// choice is systematically ~2% away from the eq (21) ω-domain
/// argmin — measured to be outside the staged fixed-point latitude
/// (round-390 CHANGELOG) and left as the documented gap of the LSF
/// chain.
#[test]
fn parameter_agreement_floors() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    // (vector, L0 floor, L1 floor, T1±2 floor).
    let floors: [(&str, f64, f64, f64); 6] = [
        ("ALGTHM", 91.0, 74.0, 77.0),
        ("FIXED", 87.0, 65.0, 74.0),
        ("LSP", 84.0, 44.0, 87.0),
        ("PITCH", 91.0, 29.0, 79.0),
        ("SPEECH", 82.0, 47.0, 75.0),
        ("TAME", 59.0, 96.0, 71.0),
    ];
    for (name, f_l0, f_l1, f_t1) in floors {
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
        // §3.9.2 gain-codeword agreement (both subframes pooled) — the
        // taming procedure acts directly on these decisions.
        let ga = pairs
            .iter()
            .map(|(o, r)| usize::from(o.ga1 == r.ga1) + usize::from(o.ga2 == r.ga2))
            .sum::<usize>();
        let gb = pairs
            .iter()
            .map(|(o, r)| usize::from(o.gb1 == r.gb1) + usize::from(o.gb2 == r.gb2))
            .sum::<usize>();

        println!(
            "{name}: frames={total} L0={:.1}% L1={:.1}% T1±2={:.1}% GA={:.1}% GB={:.1}%",
            pct(l0, total),
            pct(l1, total),
            pct(t1, total),
            pct(ga, 2 * total),
            pct(gb, 2 * total)
        );
        assert!(
            pct(l0, total) >= f_l0,
            "{name}: L0 agreement regressed to {:.1}% (floor {f_l0})",
            pct(l0, total)
        );
        assert!(
            pct(l1, total) >= f_l1,
            "{name}: L1 agreement regressed to {:.1}% (floor {f_l1})",
            pct(l1, total)
        );
        assert!(
            pct(t1, total) >= f_t1,
            "{name}: T1 agreement regressed to {:.1}% (floor {f_t1})",
            pct(t1, total)
        );
    }
}

/// Reference-locked §3.2.4 agreement: re-runs the LSP front end +
/// quantiser search with the MA history driven by the **reference's
/// own transmitted indices** every frame (the exact quantiser state
/// the reference encoder had, since its MA feedback is built from its
/// chosen indices). This removes error propagation and measures pure
/// per-frame front-end + search fidelity — the number the fixed-point
/// work ratchets. Round 390 moves the harness onto the fixed-point
/// Q13 search (`search_lsp_indices_q13` under the corpus-pinned
/// `FxLatitude::default`, integer MA history). Measured:
/// ALL-four-indices agreement ALGTHM 82.9%, FIXED 97.5%, LSP 77.7%,
/// PITCH 92.8%, SPEECH 81.2%, TAME 80.5% (round-388 float baselines
/// 82.9 / 97.5 / 77.5 / 92.8 / 80.9 / 80.5). TAME's residual misses
/// stay 24/25 pure `L0` mode flips — measured to be *systematic*
/// ~2.1% ω-domain margins inside a self-sustaining 5-frame reference
/// index cycle, not fixed-point near-ties (see the round-390
/// CHANGELOG entry for the full negative-result matrix).
#[test]
fn locked_history_lsp_agreement_floors() {
    use oxideav_g729::levinson::levinson;
    use oxideav_g729::lp_analysis::analyze_window;
    use oxideav_g729::lp_to_lsp::lp_to_lsp;
    use oxideav_g729::lsf_conversion::lsp_to_lsf_q13;
    use oxideav_g729::lsp_quantize::{
        advance_history_q13, search_lsp_indices_q13, startup_history_q13, FxLatitude,
    };
    use oxideav_g729::preprocess::Preprocessor;
    use oxideav_g729::tables::{L_WINDOW, M};

    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    // (vector, ALL-four-indices agreement floor %).
    let floors: [(&str, f64); 6] = [
        ("ALGTHM", 80.0),
        ("FIXED", 96.0),
        ("LSP", 76.5),
        ("PITCH", 91.5),
        ("SPEECH", 80.0),
        ("TAME", 78.0),
    ];
    let lat = FxLatitude::default();
    for (name, floor) in floors {
        let samples = read_pcm(&root.join(format!("g729-core/{name}.IN")));
        let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
        let n_frames = (samples.len() / SAMPLES_PER_FRAME).min(bit_bytes.len() / FRAME_BYTES);

        // The clause-3.2 front end, stage for stage as FrameEncoder
        // runs it (pre-processing, Figure-5 buffer, eqs (3)-(9), Q12
        // rounding, root search, table eq (18)).
        let mut preproc = Preprocessor::new();
        let mut speech = [0.0f32; L_WINDOW];
        let mut prev_q: [f32; M] =
            std::array::from_fn(|i| ((i + 1) as f32 * std::f32::consts::PI / 11.0).cos());
        let mut history = startup_history_q13();

        let (mut tot, mut hit) = (0usize, 0usize);
        for f in 0..n_frames {
            let mut frame = [0.0f32; SAMPLES_PER_FRAME];
            frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
            let s_new = preproc.process_frame(&frame);
            speech.copy_within(SAMPLES_PER_FRAME.., 0);
            speech[L_WINDOW - SAMPLES_PER_FRAME..].copy_from_slice(&s_new);
            let r = analyze_window(&speech);
            let lev = levinson(&r);
            let a_q12: [f32; M] = std::array::from_fn(|i| (lev.a[i] * 4096.0).round() / 4096.0);
            let q = match lp_to_lsp(&a_q12) {
                Some(q) => {
                    prev_q = q;
                    q
                }
                None => prev_q,
            };
            let omega_q13: [i32; M] = std::array::from_fn(|i| lsp_to_lsf_q13(q[i]));

            let rf = parse_frame(&bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES]).unwrap();
            let FrameKind::Active(_) = rf else { continue };
            let rp = unpack_parameters(&rf).unwrap();

            let ours = search_lsp_indices_q13(&omega_q13, &history, &lat);
            tot += 1;
            if ours.l0 == rp.l0 as usize
                && ours.l1 == rp.l1 as usize
                && ours.l2 == rp.l2 as usize
                && ours.l3 == rp.l3 as usize
            {
                hit += 1;
            }

            // Lock the MA history to the reference's transmitted indices.
            advance_history_q13(
                &mut history,
                rp.l1 as usize,
                rp.l2 as usize,
                rp.l3 as usize,
                &lat,
            );
        }
        let rate = pct(hit, tot);
        println!("{name}: locked-history LSP ALL-indices agreement {rate:.1}% ({hit}/{tot})");
        assert!(
            rate >= floor,
            "{name}: locked LSP agreement regressed to {rate:.1}% (floor {floor})"
        );
    }
}

/// Taming-procedure fingerprint of the reference encoder, read off its
/// own `.BIT` streams (docs/audio/g729/taming-procedure.md leaves the
/// trigger threshold unpinned; this test pins the corpus-consistent
/// region).
///
/// Replaying the documented accumulator recursion
/// (`E ← 1 + ĝ_p²·max(E_spanned)`, per-zone via `tab_zone`) over the
/// reference's own decoded `(pitch delay, ĝ_p)` sequence:
///
/// * the reference keeps choosing `ĝ_p > GPCLIP` while the simulated
///   accumulator climbs past 15 000 (peak ≈ 18 200 on TAME) — so the
///   reference's threshold in these units is **above** 18 200 and any
///   lower-threshold variant is contradicted by the corpus;
/// * the simulated accumulator never reaches 60 000 anywhere in the
///   corpus — so with the commonly-attributed 60 000 threshold the
///   reference **never tames** on these vectors, and neither do we
///   (behaviourally identical, which is what the flat agreement deltas
///   with taming wired confirm).
#[test]
fn reference_taming_fingerprint() {
    use oxideav_g729::gain_index_map::{demap_ga, demap_gb};
    use oxideav_g729::gain_reconstruct::reconstruct_gains;
    use oxideav_g729::pitch_decode::{decode_t2_from_p2, derive_t_min};
    use oxideav_g729::taming::{Taming, GPCLIP, THRESH_ERR};

    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let mut corpus_peak = 0.0f64;
    for name in VECTORS {
        let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
        let mut tam = Taming::new();
        let mut peak_acc_with_high_gp = 0.0f64;
        for chunk in bit_bytes.chunks_exact(FRAME_BYTES) {
            let frame = parse_frame(chunk).unwrap();
            let FrameKind::Active(_) = frame else {
                continue;
            };
            let p = unpack_parameters(&frame).unwrap();
            let t1 = decode_t1_from_p1(p.p1);
            let t2 = decode_t2_from_p2(p.p2, derive_t_min(t1.int_t));
            for (d, ga_tx, gb_tx) in [(t1, p.ga1, p.gb1), (t2, p.ga2, p.gb2)] {
                let ga = demap_ga(ga_tx as usize).unwrap();
                let gb = demap_gb(gb_tx as usize).unwrap();
                let gp = reconstruct_gains(ga, gb).unwrap().g_p_hat;
                let acc = tam.accumulators().iter().cloned().fold(0.0f64, f64::max);
                if gp > GPCLIP {
                    peak_acc_with_high_gp = peak_acc_with_high_gp.max(acc);
                }
                corpus_peak = corpus_peak.max(acc);
                tam.update(d.int_t, d.frac, gp);
            }
        }
        println!(
            "{name}: peak simulated accumulator with ĝ_p > GPCLIP = {peak_acc_with_high_gp:.1}"
        );
        if name == "TAME" {
            assert!(
                peak_acc_with_high_gp > 15_000.0,
                "TAME should pin the threshold lower bound above 15000, got {peak_acc_with_high_gp:.1}"
            );
        }
    }
    println!("corpus peak simulated accumulator = {corpus_peak:.1}");
    assert!(
        corpus_peak < THRESH_ERR,
        "reference never tames on the staged corpus; peak {corpus_peak:.1} crossed {THRESH_ERR}"
    );
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

/// Bitstream writer validation: parsing every active frame of every
/// reference `.BIT` vector, unpacking to parameters, repacking, and
/// re-serialising reproduces the original 164 bytes **byte-exactly** —
/// the pack/write path is the true inverse of the parse/unpack path on
/// the whole corpus.
#[test]
fn bitstream_writer_is_byte_exact_inverse() {
    use oxideav_g729::parameters::pack_bit_array;
    use oxideav_g729::serial::write_frame;

    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let mut checked = 0usize;
    for name in VECTORS {
        let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
        for (f, chunk) in bit_bytes.chunks_exact(FRAME_BYTES).enumerate() {
            let frame = parse_frame(chunk).unwrap();
            let FrameKind::Active(_) = frame else {
                continue;
            };
            let params = unpack_parameters(&frame).unwrap();
            let rewritten = write_frame(&pack_bit_array(&params));
            assert_eq!(
                &rewritten[..],
                chunk,
                "{name} frame {f}: re-serialisation not byte-exact"
            );
            checked += 1;
        }
    }
    println!("re-serialised {checked} frames byte-exactly");
    assert!(checked > 8000, "corpus should cover thousands of frames");
}

/// The `.IN` → `.BIT` convenience path emits parseable 164-byte frames
/// that survive the parse→unpack round-trip.
#[test]
fn encode_to_serial_roundtrips() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let samples = read_pcm(&root.join("g729-core/ALGTHM.IN"));
    let mut enc = FrameEncoder::new();
    let mut enc2 = FrameEncoder::new();
    for f in 0..samples.len() / SAMPLES_PER_FRAME {
        let mut frame = [0.0f32; SAMPLES_PER_FRAME];
        frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
        let wire = enc.encode_frame_to_serial(&frame);
        let expect = enc2.encode_frame(&frame).params;
        let parsed = parse_frame(&wire).unwrap();
        let got = unpack_parameters(&parsed).unwrap();
        assert_eq!(got, expect, "frame {f}: serial round-trip mismatch");
    }
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
