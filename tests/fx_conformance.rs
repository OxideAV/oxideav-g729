//! Conformance tracking for the fixed-point (`fx`) decode path against
//! the staged ITU `.PST` references.
//!
//! Stage 1 (this harness): the §4.1 chain runs on the clause-5
//! Word16/Word32 grid (`fx::decoder::FrameDecoderFx`) and its
//! pre-postfilter speech `ŝ(n)` is fed through the *float* §4.2
//! cascade — isolating the excitation/synthesis fidelity from the
//! not-yet-fixed-point postfilter. Metrics per clean vector:
//! correlation, RMS ratio, bit-exact share.
//!
//! When the corpus is absent (published-crate build) the tests skip.

use std::path::{Path, PathBuf};

use oxideav_g729::decode_chain::{DecodedFrame, SubframeDecode};
use oxideav_g729::fixed_codebook::decode_pulses;
use oxideav_g729::fx::decoder::FrameDecoderFx;
use oxideav_g729::lp_synthesis::{SynthesizedFrame, SynthesizedSubframe};
use oxideav_g729::parameters::unpack_parameters;
use oxideav_g729::pitch_decode::PitchDelay;
use oxideav_g729::postfilter_cascade::PostfilterCascade;
use oxideav_g729::serial::{self, FrameKind, FRAME_BYTES};

const SAMPLES_PER_FRAME: usize = 80;
const SUBFRAME: usize = 40;
const M: usize = 10;

const CLEAN_VECTORS: [&str; 6] = ["ALGTHM", "FIXED", "LSP", "PITCH", "SPEECH", "TAME"];
const CORPORA: [&str; 2] = ["g729-core", "g729a"];

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

fn read_pst(path: &Path) -> Vec<i16> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    bytes
        .chunks_exact(2)
        .map(|c| i16::from_le_bytes([c[0], c[1]]))
        .collect()
}

/// Decode a `.BIT` stream through the fx §4.1 chain + the float §4.2
/// cascade. Erasure sentinels are skipped (clean vectors carry none).
fn decode_fx_hybrid(label: &str, bit: &[u8]) -> Vec<f32> {
    let mut fx = FrameDecoderFx::new();
    let mut pf = PostfilterCascade::new();
    let n_frames = bit.len() / FRAME_BYTES;
    let mut out = Vec::with_capacity(n_frames * SAMPLES_PER_FRAME);

    for f in 0..n_frames {
        let frame = &bit[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).unwrap_or_else(|e| panic!("{label}: {e:?}"));
        if matches!(kind, FrameKind::Erased) {
            continue;
        }
        let params = unpack_parameters(&kind).unwrap_or_else(|e| panic!("{label}: {e:?}"));
        let dec = fx.decode_frame(&params);

        // Feed the float cascade: it consumes ŝ, the per-subframe LP
        // sets, and int(T1) of subframe 1.
        let synth = SynthesizedFrame {
            subframes: [0, 1].map(|s| SynthesizedSubframe {
                adaptive: [0.0; SUBFRAME],
                excitation: [0.0; SUBFRAME],
                speech: std::array::from_fn(|n| f32::from(dec.speech[s * SUBFRAME + n])),
            }),
        };
        let pulses = decode_pulses(0, 0).unwrap();
        let frame_info = DecodedFrame {
            params,
            parity_ok: true,
            omega: [0.0; M],
            t_min: 0,
            subframes: [0, 1].map(|s| SubframeDecode {
                lsp_q: [0.0; M],
                lp: std::array::from_fn(|i| f32::from(dec.sub[s].a_q12[i + 1]) / 4096.0),
                pitch: PitchDelay {
                    int_t: dec.sub[s].t_int,
                    frac: 0,
                },
                fixed: pulses,
                beta: 0.0,
                codevector: [0.0; SUBFRAME],
                gains: oxideav_g729::gain_reconstruct::QuantisedGains {
                    g_p_hat: f32::from(dec.sub[s].gains.gain_pit_q14) / 16384.0,
                    gamma_hat: 0.0,
                },
                predicted: oxideav_g729::gain_predict::PredictedGain {
                    e_db: 0.0,
                    e_tilde_db: 0.0,
                    g_c_prime: 0.0,
                },
                g_c_hat: f32::from(dec.sub[s].gains.gain_code_q1) / 2.0,
            }),
        };
        let post = pf.process_frame(&synth, &frame_info);
        out.extend_from_slice(&post.output());
    }
    out
}

struct Metrics {
    corr: f64,
    rms_ratio: f64,
    exact_pct: f64,
    max_delta: f64,
}

fn metrics(out: &[f32], reference: &[i16]) -> Metrics {
    let n = out.len().min(reference.len());
    let (mut dot, mut oe, mut re) = (0.0f64, 0.0f64, 0.0f64);
    let mut exact = 0usize;
    let mut max_delta = 0.0f64;
    for i in 0..n {
        let o = f64::from(out[i]);
        let r = f64::from(reference[i]);
        dot += o * r;
        oe += o * o;
        re += r * r;
        let d = (o - r).abs();
        max_delta = max_delta.max(d);
        if out[i].round().clamp(-32768.0, 32767.0) as i16 == reference[i] {
            exact += 1;
        }
    }
    Metrics {
        corr: if oe > 0.0 && re > 0.0 {
            dot / (oe.sqrt() * re.sqrt())
        } else {
            0.0
        },
        rms_ratio: (oe / re.max(1.0)).sqrt(),
        exact_pct: 100.0 * exact as f64 / n as f64,
        max_delta,
    }
}

/// Pinned correlation floors for the fx-§4.1 + float-§4.2 hybrid
/// (round-410 measurements 0.995–0.9998, floors a safety margin
/// below; ratchet upward as the fixed-point §4.2 cascade lands).
fn corr_floor(_corpus: &str, name: &str) -> f64 {
    match name {
        "TAME" => 0.995,
        _ => 0.99,
    }
}

#[test]
fn fx_hybrid_exactness_metrics() {
    let Some(root) = conformance_root() else {
        eprintln!("g729 conformance corpus absent — skipping fx hybrid metrics");
        return;
    };

    let mut checked = 0usize;
    for corpus in CORPORA {
        for name in CLEAN_VECTORS {
            let label = format!("{corpus}/{name}");
            let bit_path = root.join(format!("{corpus}/{name}.BIT"));
            let pst_path = root.join(format!("{corpus}/{name}.PST"));
            if !bit_path.is_file() || !pst_path.is_file() {
                continue;
            }
            let bit = std::fs::read(&bit_path).unwrap();
            let reference = read_pst(&pst_path);
            let out = decode_fx_hybrid(&label, &bit);
            let m = metrics(&out, &reference);
            eprintln!(
                "{label}: fx-hybrid corr {:.4}  rms {:.2}x  exact {:.2}%  max|d| {:.0}",
                m.corr, m.rms_ratio, m.exact_pct, m.max_delta
            );
            assert!(
                m.corr >= corr_floor(corpus, name),
                "{label}: fx hybrid correlation {:.4} under floor",
                m.corr
            );
            checked += 1;
        }
    }
    assert!(checked >= 12, "checked only {checked} vectors");
}

/// Decode a `.BIT` stream that may contain §4.4 erasure sentinels:
/// erased frames route through the fx concealment primitives, active
/// frames through the normal chain; the float §4.2 cascade runs on
/// every frame's output.
fn decode_fx_hybrid_concealed(label: &str, bit: &[u8]) -> Vec<f32> {
    let mut fx = FrameDecoderFx::new();
    let mut pf = PostfilterCascade::new();
    let n_frames = bit.len() / FRAME_BYTES;
    let mut out = Vec::with_capacity(n_frames * SAMPLES_PER_FRAME);

    for f in 0..n_frames {
        let frame = &bit[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).unwrap_or_else(|e| panic!("{label}: {e:?}"));
        let (dec, params) = if matches!(kind, FrameKind::Erased) {
            (fx.decode_erased_frame(), zero_params())
        } else {
            let params = unpack_parameters(&kind).unwrap_or_else(|e| panic!("{label}: {e:?}"));
            (fx.decode_frame(&params), params)
        };
        let synth = SynthesizedFrame {
            subframes: [0, 1].map(|s| SynthesizedSubframe {
                adaptive: [0.0; SUBFRAME],
                excitation: [0.0; SUBFRAME],
                speech: std::array::from_fn(|n| f32::from(dec.speech[s * SUBFRAME + n])),
            }),
        };
        let frame_info = frame_info_for(&dec, params);
        let post = pf.process_frame(&synth, &frame_info);
        out.extend_from_slice(&post.output());
    }
    out
}

fn zero_params() -> oxideav_g729::parameters::Parameters {
    oxideav_g729::parameters::Parameters {
        l0: 0,
        l1: 0,
        l2: 0,
        l3: 0,
        p1: 0,
        p0: 0,
        c1: 0,
        s1: 0,
        ga1: 0,
        gb1: 0,
        p2: 0,
        c2: 0,
        s2: 0,
        ga2: 0,
        gb2: 0,
    }
}

fn frame_info_for(
    dec: &oxideav_g729::fx::decoder::DecodedFrameFx,
    params: oxideav_g729::parameters::Parameters,
) -> DecodedFrame {
    let pulses = decode_pulses(0, 0).unwrap();
    DecodedFrame {
        params,
        parity_ok: true,
        omega: [0.0; M],
        t_min: 0,
        subframes: [0, 1].map(|s| SubframeDecode {
            lsp_q: [0.0; M],
            lp: std::array::from_fn(|i| f32::from(dec.sub[s].a_q12[i + 1]) / 4096.0),
            pitch: PitchDelay {
                int_t: dec.sub[s].t_int,
                frac: 0,
            },
            fixed: pulses,
            beta: 0.0,
            codevector: [0.0; SUBFRAME],
            gains: oxideav_g729::gain_reconstruct::QuantisedGains {
                g_p_hat: f32::from(dec.sub[s].gains.gain_pit_q14) / 16384.0,
                gamma_hat: 0.0,
            },
            predicted: oxideav_g729::gain_predict::PredictedGain {
                e_db: 0.0,
                e_tilde_db: 0.0,
                g_c_prime: 0.0,
            },
            g_c_hat: f32::from(dec.sub[s].gains.gain_code_q1) / 2.0,
        }),
    }
}

/// The decoder-only stress vectors through the fx chain: PARITY (the
/// §4.1.2 T1-substitution path), OVERFLOW (the 16-bit saturation
/// behaviour), and ERASURE (the §4.4 concealment primitives). Floors
/// pinned from the round-410 measurements.
#[test]
fn fx_hybrid_stress_vectors() {
    let Some(root) = conformance_root() else {
        eprintln!("g729 conformance corpus absent — skipping fx stress vectors");
        return;
    };

    // (vector, per-corpus correlation floors [g729-core, g729a])
    let cases = [
        // PARITY exercises the §4.1.2 T1 substitution — measured 0.9987.
        ("PARITY", [0.99, 0.99]),
        // OVERFLOW exercises the 16-bit overflow-rescale protocol —
        // measured 0.78 (base) / 0.25 (g729a decodes with the Annex A
        // reduced decoder whose overflow behaviour this base chain
        // does not model).
        ("OVERFLOW", [0.70, 0.15]),
        // ERASURE exercises §4.4 concealment (voicing class defaults
        // to non-periodic in this hybrid) — measured 0.91 / 0.94.
        ("ERASURE", [0.88, 0.90]),
    ];
    for (name, floors) in cases {
        for (ci, corpus) in CORPORA.iter().enumerate() {
            let label = format!("{corpus}/{name}");
            let bit_path = root.join(format!("{corpus}/{name}.BIT"));
            let pst_path = root.join(format!("{corpus}/{name}.PST"));
            if !bit_path.is_file() || !pst_path.is_file() {
                continue;
            }
            let bit = std::fs::read(&bit_path).unwrap();
            let reference = read_pst(&pst_path);
            let out = decode_fx_hybrid_concealed(&label, &bit);
            let m = metrics(&out, &reference);
            eprintln!(
                "{label}: fx-hybrid corr {:.4}  rms {:.2}x  exact {:.2}%  max|d| {:.0}",
                m.corr, m.rms_ratio, m.exact_pct, m.max_delta
            );
            for (i, &s) in out.iter().enumerate() {
                assert!(s.is_finite(), "{label}: sample {i} not finite");
            }
            assert!(
                m.corr >= floors[ci],
                "{label}: corr {:.4} under pinned floor {:.4}",
                m.corr,
                floors[ci]
            );
        }
    }
}
