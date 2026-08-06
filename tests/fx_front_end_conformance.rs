//! Fixed-point §3.2.1–§3.2.3 front-end conformance A/B: runs the
//! reference-locked LSP-index agreement harness (MA history driven by
//! the reference's own transmitted indices — the exact quantiser
//! state the reference encoder had) once through the float-emulated
//! front end and once through the genuine Word16/Word32 chain
//! (`fx::analysis` → `fx::levinson` → `fx::lp_to_lsp`), and sweeps
//! the unstated overflow-rescale latitude black-box against the
//! corpus.
//!
//! Skips clean when the conformance corpus is absent (published-crate
//! build), mirroring the sibling harnesses.

use std::path::{Path, PathBuf};

use oxideav_g729::fx::analysis::{analyze_window_fx, FrontEndLatitude};
use oxideav_g729::fx::levinson::levinson_fx;
use oxideav_g729::fx::lp_to_lsp::lp_to_lsp_fx;
use oxideav_g729::levinson::levinson;
use oxideav_g729::lp_analysis::analyze_window;
use oxideav_g729::lp_to_lsp::lp_to_lsp;
use oxideav_g729::lsf_conversion::{acos_q15_to_lsf_q13, lsp_to_lsf_q13};
use oxideav_g729::lsp_quantize::{
    advance_history_q13, search_lsp_indices_q13, startup_history_q13, FxLatitude,
};
use oxideav_g729::parameters::unpack_parameters;
use oxideav_g729::preprocess::Preprocessor;
use oxideav_g729::serial::{parse_frame, FrameKind, FRAME_BYTES};
use oxideav_g729::tables::{L_WINDOW, M};

const SAMPLES_PER_FRAME: usize = 80;
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

fn read_pcm(path: &Path) -> Vec<f32> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    bytes
        .chunks_exact(2)
        .map(|c| f32::from(i16::from_le_bytes([c[0], c[1]])))
        .collect()
}

/// Which front end produces the unquantised ω vector for the frame.
#[derive(Clone, Copy)]
enum FrontEnd {
    Float,
    Fixed(FrontEndLatitude),
    /// Fixed front end fed by the §3.1 fixed-point pre-processor
    /// (Word32 feedback) instead of the float filter.
    FixedPre(FrontEndLatitude),
}

/// Reference-locked ALL-four-indices agreement for one vector under
/// the chosen front end. The MA history advances with the reference's
/// transmitted `(L1, L2, L3)` every frame; only the per-frame front
/// end + search fidelity is measured.
fn locked_agreement(root: &Path, name: &str, front: FrontEnd) -> (usize, usize) {
    let samples = read_pcm(&root.join(format!("g729-core/{name}.IN")));
    let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
    let n_frames = (samples.len() / SAMPLES_PER_FRAME).min(bit_bytes.len() / FRAME_BYTES);

    let mut preproc = Preprocessor::new();
    let mut preproc_fx = oxideav_g729::fx::analysis::PreprocessorFx::new();
    let mut speech = [0.0f32; L_WINDOW];
    let mut prev_q_float: [f32; M] =
        std::array::from_fn(|i| ((i + 1) as f32 * std::f32::consts::PI / 11.0).cos());
    let mut prev_q_fx: [i16; M] = oxideav_g729::fx::lsp::STARTUP_LSP_Q15;
    let mut history = startup_history_q13();
    let lat = FxLatitude::default();

    let (mut tot, mut hit) = (0usize, 0usize);
    for f in 0..n_frames {
        let mut frame = [0.0f32; SAMPLES_PER_FRAME];
        frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
        speech.copy_within(SAMPLES_PER_FRAME.., 0);
        if matches!(front, FrontEnd::FixedPre(_)) {
            let frame_i: [i16; SAMPLES_PER_FRAME] = std::array::from_fn(|n| frame[n] as i16);
            let s_new = preproc_fx.process_frame(&frame_i);
            for (o, &v) in speech[L_WINDOW - SAMPLES_PER_FRAME..]
                .iter_mut()
                .zip(&s_new)
            {
                *o = f32::from(v);
            }
        } else {
            let s_new = preproc.process_frame(&frame);
            speech[L_WINDOW - SAMPLES_PER_FRAME..].copy_from_slice(&s_new);
        }

        let omega_q13: [i32; M] = match front {
            FrontEnd::Float => {
                let r = analyze_window(&speech);
                let lev = levinson(&r);
                let a_q12: [f32; M] = std::array::from_fn(|i| (lev.a[i] * 4096.0).round() / 4096.0);
                let q = match lp_to_lsp(&a_q12) {
                    Some(q) => {
                        prev_q_float = q;
                        q
                    }
                    None => prev_q_float,
                };
                std::array::from_fn(|i| lsp_to_lsf_q13(q[i]))
            }
            FrontEnd::Fixed(fe_lat) | FrontEnd::FixedPre(fe_lat) => {
                // The preprocessed samples sit on the 16-bit integer
                // grid (round-385 pin), so the cast is lossless.
                let speech_i: [i16; L_WINDOW] = std::array::from_fn(|n| speech[n] as i16);
                let ac = analyze_window_fx(&speech_i, &fe_lat);
                let q = levinson_fx(&ac.r)
                    .and_then(|lev| lp_to_lsp_fx(&lev.a_q12))
                    .inspect(|q| prev_q_fx = *q)
                    .unwrap_or(prev_q_fx);
                std::array::from_fn(|i| acos_q15_to_lsf_q13(i32::from(q[i])))
            }
        };

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
        advance_history_q13(
            &mut history,
            rp.l1 as usize,
            rp.l2 as usize,
            rp.l3 as usize,
            &lat,
        );
    }
    (hit, tot)
}

/// The A/B: fixed-point front end vs the float emulation, plus the
/// overflow-shift latitude sweep and the §3.1 fixed pre-processor.
/// Measured under the pinned configuration (overflow shift 2,
/// truncating Q11 poly landing, rounding interpolation step,
/// Levinson numerator pre-shift; `fx+pre` adds the Q16-feedback
/// Word32 §3.1 filter — the production `FrameEncoder` chain):
///
/// | vector | float | fx    | fx+pre |
/// |--------|-------|-------|--------|
/// | ALGTHM | 82.9% | 97.1% | 94.3%  |
/// | FIXED  | 97.5% | 97.5% | 97.5%  |
/// | LSP    | 77.7% | 81.3% | 85.8%  |
/// | PITCH  | 92.8% | 93.0% | 93.6%  |
/// | SPEECH | 81.2% | 83.4% | 91.7%  |
/// | TAME   | 80.5% | 99.2% | 99.2%  |
///
/// TAME's jump retires the round-390 "unstaged structural element of
/// the reference's final mode compare" hypothesis: the 24
/// locked-history `L0` flips were front-end ω divergence — the
/// genuine Word32 overflow-rescale truncation pattern moves ω onto
/// the reference's grid and the printed eq (21) mode compare then
/// picks the reference's predictor. The `overflow_shift = 4` variant
/// collapses TAME to 60.9% (over-truncation), pinning the default.
#[test]
fn fx_front_end_locked_agreement() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };

    // Float-emulation baselines (the round-390 measured rates).
    let mut float_rate = [0.0f64; VECTORS.len()];
    for (i, name) in VECTORS.iter().enumerate() {
        let (hit, tot) = locked_agreement(&root, name, FrontEnd::Float);
        float_rate[i] = 100.0 * hit as f64 / tot as f64;
        println!("{name}: float {:.1}% ({hit}/{tot})", float_rate[i]);
    }

    // Regression floors ~2 points under the measured fx rates.
    let fx_floors: [f64; VECTORS.len()] = [94.0, 95.5, 79.0, 91.0, 81.4, 97.0];

    // §3.1 A/B: the fixed pre-processor feeding the fixed front end
    // (the production FrameEncoder chain) with its own floors.
    let pre_floors: [f64; VECTORS.len()] = [91.0, 95.5, 83.0, 91.5, 89.0, 97.0];
    for (i, name) in VECTORS.iter().enumerate() {
        let (hit, tot) =
            locked_agreement(&root, name, FrontEnd::FixedPre(FrontEndLatitude::default()));
        let rate = 100.0 * hit as f64 / tot as f64;
        println!(
            "{name}: fx+pre {rate:.1}% ({hit}/{tot})  [float {:.1}%]",
            float_rate[i]
        );
        assert!(
            rate >= pre_floors[i],
            "{name}: fx+pre {rate:.1}% under floor {:.1}%",
            pre_floors[i]
        );
    }

    for shift in [1i16, 2, 4] {
        let fe = FrontEndLatitude {
            overflow_shift: shift,
        };
        for (i, name) in VECTORS.iter().enumerate() {
            let (hit, tot) = locked_agreement(&root, name, FrontEnd::Fixed(fe));
            let rate = 100.0 * hit as f64 / tot as f64;
            println!(
                "{name}: fx(shift {shift}) {rate:.1}% ({hit}/{tot})  [float {:.1}%]",
                float_rate[i]
            );
            if shift == FrontEndLatitude::default().overflow_shift {
                assert!(
                    rate >= fx_floors[i],
                    "{name}: fx front end {rate:.1}% under floor {:.1}%",
                    fx_floors[i]
                );
                assert!(
                    rate >= float_rate[i] - 1.0,
                    "{name}: fx front end {rate:.1}% regresses the float baseline {:.1}%",
                    float_rate[i]
                );
            }
        }
    }
}
