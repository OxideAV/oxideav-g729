//! Full fixed-point decode conformance: the §4.1 chain
//! (`fx::decoder::FrameDecoderFx`) **plus** the §4.2 cascade
//! (`fx::postfilter::PostfilterFx`) against the staged ITU `.PST`
//! references — the byte-exactness ratchet.
//!
//! Beyond the aggregate metrics this harness reports the **first
//! diverging sample** of every vector (absolute index, frame,
//! subframe, in-subframe offset, ours vs reference) — the bisection
//! anchor for the clause-5 residual-divergence hunt.
//!
//! When the corpus is absent (published-crate build) the tests skip.

use std::path::{Path, PathBuf};

use oxideav_g729::fx::decoder::FrameDecoderFx;
use oxideav_g729::fx::postfilter::PostfilterFx;
use oxideav_g729::parameters::unpack_parameters;
use oxideav_g729::serial::{self, FrameKind, FRAME_BYTES};

const SAMPLES_PER_FRAME: usize = 80;
const SUBFRAME: usize = 40;

const CLEAN_VECTORS: [&str; 6] = ["ALGTHM", "FIXED", "LSP", "PITCH", "SPEECH", "TAME"];

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

/// Decode a `.BIT` stream through the full fixed-point chain
/// (§4.1 + §4.2), routing erasure sentinels through the §4.4
/// concealment primitives and latching the §4.4 voicing class from
/// the §4.2.1 long-term decisions of good frames.
fn decode_fx_full(label: &str, bit: &[u8]) -> Vec<i16> {
    let mut fx = FrameDecoderFx::new();
    let mut pf = PostfilterFx::new();
    let n_frames = bit.len() / FRAME_BYTES;
    let mut out = Vec::with_capacity(n_frames * SAMPLES_PER_FRAME);

    for f in 0..n_frames {
        let frame = &bit[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).unwrap_or_else(|e| panic!("{label}: {e:?}"));
        let dec = if matches!(kind, FrameKind::Erased) {
            fx.decode_erased_frame()
        } else {
            let params = unpack_parameters(&kind).unwrap_or_else(|e| panic!("{label}: {e:?}"));
            fx.decode_frame(&params)
        };

        // Clause 4.2.1: both subframes anchor on int(T_1) of subframe 1.
        let int_t1 = usize::try_from(dec.sub[0].t_int.max(1)).unwrap();
        let mut periodic = false;
        for s in 0..2 {
            let speech: [i16; SUBFRAME] = std::array::from_fn(|n| dec.speech[s * SUBFRAME + n]);
            let (pcm, decision) = pf.process_subframe(&speech, &dec.sub[s].a_q12, int_t1);
            periodic |= decision.gain_q15 > 0;
            out.extend_from_slice(&pcm);
        }
        // §4.4 voicing classifier: periodic iff any subframe's
        // long-term prediction gain cleared the eq (82) threshold.
        if !matches!(kind, FrameKind::Erased) {
            fx.erasure_periodic = periodic;
        }
    }
    out
}

struct Metrics {
    corr: f64,
    rms_ratio: f64,
    exact_pct: f64,
    max_delta: i32,
    first_div: Option<usize>,
    /// (start, length) of the longest byte-exact run.
    longest_run: (usize, usize),
    /// Frames whose full 80 samples are byte-exact.
    clean_frames: usize,
    total_frames: usize,
}

fn metrics(out: &[i16], reference: &[i16]) -> Metrics {
    let n = out.len().min(reference.len());
    let (mut dot, mut oe, mut re) = (0.0f64, 0.0f64, 0.0f64);
    let mut exact = 0usize;
    let mut max_delta = 0i32;
    let mut first_div = None;
    let mut longest_run = (0usize, 0usize);
    let mut run_start = 0usize;
    let mut in_run = false;
    for i in 0..n {
        let o = f64::from(out[i]);
        let r = f64::from(reference[i]);
        dot += o * r;
        oe += o * o;
        re += r * r;
        let d = (i32::from(out[i]) - i32::from(reference[i])).abs();
        max_delta = max_delta.max(d);
        if d == 0 {
            exact += 1;
            if !in_run {
                run_start = i;
                in_run = true;
            }
            if i + 1 - run_start > longest_run.1 {
                longest_run = (run_start, i + 1 - run_start);
            }
        } else {
            in_run = false;
            if first_div.is_none() {
                first_div = Some(i);
            }
        }
    }
    let total_frames = n / SAMPLES_PER_FRAME;
    let mut clean_frames = 0usize;
    for f in 0..total_frames {
        let base = f * SAMPLES_PER_FRAME;
        if (0..SAMPLES_PER_FRAME).all(|k| out[base + k] == reference[base + k]) {
            clean_frames += 1;
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
        first_div,
        longest_run,
        clean_frames,
        total_frames,
    }
}

fn report(label: &str, out: &[i16], reference: &[i16]) -> Metrics {
    let m = metrics(out, reference);
    match m.first_div {
        Some(i) => {
            let frame = i / SAMPLES_PER_FRAME;
            let sub = (i % SAMPLES_PER_FRAME) / SUBFRAME;
            let off = i % SUBFRAME;
            eprintln!(
                "{label}: corr {:.5}  rms {:.3}x  exact {:.2}%  max|d| {}  first-div @{} (frame {frame} sub {sub} n {off}: {} vs {})  run {}@{}  clean {}/{}",
                m.corr,
                m.rms_ratio,
                m.exact_pct,
                m.max_delta,
                i,
                out[i],
                reference[i],
                m.longest_run.1,
                m.longest_run.0,
                m.clean_frames,
                m.total_frames
            );
        }
        None => eprintln!(
            "{label}: corr {:.5}  rms {:.3}x  exact {:.2}%  BYTE-EXACT",
            m.corr, m.rms_ratio, m.exact_pct
        ),
    }
    m
}

/// Per-vector pinned floors `(corr, exact%)` for the full fixed-point
/// chain on the clean vectors, both corpora (measured r419:
/// ALGTHM 0.9927/0.9944, FIXED 0.9502/0.9756, LSP 0.996+,
/// PITCH 0.9961+, SPEECH 0.9973+, TAME 0.9994+; exact 0.4–27%).
fn floors(name: &str) -> (f64, f64) {
    match name {
        // FIXED is the §4.2.1 stress case: the long-term postfilter's
        // enable/gain decisions on onset material still diverge from
        // the reference (LT-off measures 0.986/0.992) — the top open
        // divergence of the hunt.
        "FIXED" => (0.94, 10.0),
        "TAME" => (0.995, 0.3),
        _ => (0.99, 1.0),
    }
}

#[test]
fn fx_full_clean_vectors() {
    let Some(root) = conformance_root() else {
        eprintln!("g729 conformance corpus absent — skipping fx full metrics");
        return;
    };

    let mut checked = 0usize;
    for corpus in ["g729-core", "g729a"] {
        for name in CLEAN_VECTORS {
            let label = format!("{corpus}/{name}");
            let bit_path = root.join(format!("{corpus}/{name}.BIT"));
            let pst_path = root.join(format!("{corpus}/{name}.PST"));
            if !bit_path.is_file() || !pst_path.is_file() {
                continue;
            }
            let bit = std::fs::read(&bit_path).unwrap();
            let reference = read_pst(&pst_path);
            let out = decode_fx_full(&label, &bit);
            let m = report(&label, &out, &reference);
            let (corr_floor, exact_floor) = floors(name);
            assert!(
                m.corr >= corr_floor,
                "{label}: corr {:.4} under floor {corr_floor}",
                m.corr
            );
            assert!(
                m.exact_pct >= exact_floor,
                "{label}: exact {:.2}% under floor {exact_floor}%",
                m.exact_pct
            );
            checked += 1;
        }
    }
    assert!(checked >= 12, "checked only {checked} vectors");
}

/// Stage-by-stage dump of the first frames of one vector — the manual
/// bisection instrument. Gated on `G729_FX_TRACE=<corpus>/<name>:<frames>`
/// (e.g. `G729_FX_TRACE=g729-core/FIXED:2`); a no-op otherwise.
#[test]
fn fx_full_trace_dump() {
    let Ok(spec) = std::env::var("G729_FX_TRACE") else {
        return;
    };
    let Some(root) = conformance_root() else {
        return;
    };
    let (vector, n_frames) = spec.split_once(':').unwrap_or((spec.as_str(), "2"));
    let n_frames: usize = n_frames.parse().unwrap();
    let bit = std::fs::read(root.join(format!("{vector}.BIT"))).unwrap();
    let reference = read_pst(&root.join(format!("{vector}.PST")));

    let mut fx = FrameDecoderFx::new();
    let mut pf = PostfilterFx::new();
    for f in 0..n_frames.min(bit.len() / FRAME_BYTES) {
        let frame = &bit[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).unwrap();
        let dec = if matches!(kind, FrameKind::Erased) {
            fx.decode_erased_frame()
        } else {
            let params = unpack_parameters(&kind).unwrap();
            eprintln!("--- {vector} frame {f} params {params:?}");
            fx.decode_frame(&params)
        };
        let exc = fx.last_frame_excitation();
        let int_t1 = usize::try_from(dec.sub[0].t_int.max(1)).unwrap();
        for s in 0..2 {
            let speech: [i16; SUBFRAME] = std::array::from_fn(|n| dec.speech[s * SUBFRAME + n]);
            let t = pf.process_subframe_traced(&speech, &dec.sub[s].a_q12, int_t1);
            let base = f * SAMPLES_PER_FRAME + s * SUBFRAME;
            let refsl = &reference[base..(base + SUBFRAME).min(reference.len())];
            eprintln!(
                "=== {vector} frame {f} sub {s} (t_int {} int_t1 {int_t1})",
                dec.sub[s].t_int
            );
            eprintln!("a_q12   {:?}", dec.sub[s].a_q12);
            eprintln!("gains   {:?}", dec.sub[s].gains);
            eprintln!("lt-dec  {:?}  gf_q12 {}", t.decision, t.gf_q12);
            eprintln!("agc g   {} -> {}", t.agc_gain_in_q12, t.agc_gain_out_q12);
            eprintln!("exc     {:?}", &exc[s * SUBFRAME..s * SUBFRAME + 12]);
            eprintln!("speech  {:?}", &speech[..12]);
            eprintln!("lt      {:?}", &t.long_term[..12]);
            eprintln!("st      {:?}", &t.short_term[..12]);
            eprintln!("tilt    {:?}", &t.tilt[..12]);
            eprintln!("agc     {:?}", &t.agc[..12]);
            eprintln!("out     {:?}", &t.output[..12]);
            eprintln!("ref     {:?}", &refsl[..12.min(refsl.len())]);
        }
    }
}

/// The decoder-only stress vectors through the full fx chain: PARITY
/// (§4.1.2 T1 substitution), OVERFLOW (16-bit saturation behaviour),
/// ERASURE (§4.4 concealment with the §4.2.1-latched voicing class).
#[test]
fn fx_full_stress_vectors() {
    let Some(root) = conformance_root() else {
        eprintln!("g729 conformance corpus absent — skipping fx full stress");
        return;
    };

    // (vector, per-corpus correlation floors [g729-core, g729a]).
    // g729a's OVERFLOW reference decodes with the Annex A reduced
    // decoder whose §A.4 overflow behaviour this base chain does not
    // model — its correlation is near zero by construction (measured
    // −0.05; the float-§4.2 hybrid measures 0.25).
    let cases = [
        ("PARITY", [0.99, 0.99]),
        ("OVERFLOW", [0.70, -0.30]),
        // ERASURE measured 0.922 / 0.887 (the concealed stretches
        // re-sync a little differently through the fx cascade than
        // the float hybrid's 0.91/0.94).
        ("ERASURE", [0.88, 0.86]),
    ];
    for (name, floors) in cases {
        for (ci, corpus) in ["g729-core", "g729a"].iter().enumerate() {
            let label = format!("{corpus}/{name}");
            let bit_path = root.join(format!("{corpus}/{name}.BIT"));
            let pst_path = root.join(format!("{corpus}/{name}.PST"));
            if !bit_path.is_file() || !pst_path.is_file() {
                continue;
            }
            let bit = std::fs::read(&bit_path).unwrap();
            let reference = read_pst(&pst_path);
            let out = decode_fx_full(&label, &bit);
            let m = report(&label, &out, &reference);
            assert!(
                m.corr >= floors[ci],
                "{label}: corr {:.4} under pinned floor {:.4}",
                m.corr,
                floors[ci]
            );
        }
    }
}
