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
    let annex_a = label.starts_with("g729a");
    let mut fx = FrameDecoderFx::new();
    let mut pf = if annex_a {
        PostfilterFx::new_annex_a()
    } else {
        PostfilterFx::new()
    };
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

        // Clause 4.2.1: both subframes anchor on int(T_1) of subframe 1;
        // §A.4.2.1 anchors each subframe on its own int(T).
        let int_t1 = usize::try_from(dec.sub[0].t_int.max(1)).unwrap();
        let mut periodic = false;
        for s in 0..2 {
            let anchor = if annex_a {
                usize::try_from(dec.sub[s].t_int.max(1)).unwrap()
            } else {
                int_t1
            };
            let speech: [i16; SUBFRAME] = std::array::from_fn(|n| dec.speech[s * SUBFRAME + n]);
            let (pcm, decision) = pf.process_subframe(&speech, &dec.sub[s].a_q12, anchor);
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
    diff_count: usize,
    total: usize,
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
        diff_count: n - exact,
        total: n,
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
                "{label}: corr {:.5}  rms {:.3}x  exact {:.2}% (≠ {} of {})  max|d| {}  first-div @{} (frame {frame} sub {sub} n {off}: {} vs {})  run {}@{}  clean {}/{}",
                m.corr,
                m.rms_ratio,
                m.exact_pct,
                m.diff_count,
                m.total,
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
/// chain on the clean vectors, per corpus (measured r461 after the
/// eq (78) fix, the eq (90) re-pin, the eq (83) clamp and the Annex A
/// cascade — base / g729a: ALGTHM 7.82 / 6.43, FIXED 57.07 / 27.79,
/// LSP 5.21 / 5.32, PITCH 7.13 / 5.94, SPEECH 24.51 / 18.04, TAME
/// 1.10 / 0.69 exact %; corr 0.9988–0.99999). The r452 numbers were
/// ALGTHM 4.04 / 3.61, FIXED 34.45 / 20.14, LSP 4.01 / 3.72, PITCH
/// 1.98 / 1.90, SPEECH 21.80 / 15.76, TAME 0.81 / 0.48.
fn floors(corpus: &str, name: &str) -> (f64, f64) {
    match (corpus, name) {
        ("g729-core", "ALGTHM") => (0.9999, 7.0),
        ("g729-core", "FIXED") => (0.9999, 50.0),
        ("g729-core", "LSP") => (0.998, 4.8),
        ("g729-core", "PITCH") => (0.9999, 6.5),
        ("g729-core", "SPEECH") => (0.9999, 23.0),
        ("g729-core", "TAME") => (0.9999, 1.0),
        (_, "ALGTHM") => (0.9998, 6.0),
        (_, "FIXED") => (0.9998, 25.0),
        (_, "LSP") => (0.999, 5.0),
        (_, "PITCH") => (0.9999, 5.5),
        (_, "SPEECH") => (0.9999, 17.0),
        (_, "TAME") => (0.9999, 0.6),
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
            let (corr_floor, exact_floor) = floors(corpus, name);
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

/// Exact inversion of the fixed-point §4.2.5 output stage: recovers
/// the reference decoder's AGC output `2·sf′(n)` (the Q1 grid the
/// high-pass consumes) from a `.PST` sample sequence, running the
/// crate's own high-pass model forward over a beam of candidate input
/// paths and keeping every path that reproduces the reference output
/// exactly. Because the eq (91) filter has a double zero at DC, the
/// input is recoverable only up to a slowly drifting offset; the beam
/// is ranked by distance from `guide` (our own AGC output) so the
/// physically plausible path wins the ties.
///
/// Returns `(recovered, first_inconsistent)`: `recovered[n]` is the
/// best path's input at `n`; `first_inconsistent` is the first index
/// at which no candidate reproduced the reference (the high-pass model
/// itself is wrong there) or `None` when the whole vector inverts.
fn recover_agc_output(reference: &[i16], guide: &[i16]) -> (Vec<i32>, Option<usize>) {
    const B: [i32; 3] = [7699, -15398, 7699];
    const A: [i32; 3] = [8192, 15836, -7667];
    const BEAM: usize = 32;
    #[derive(Clone)]
    struct Path {
        hx: [i32; 2],
        hy: [i32; 2],
        xs: Vec<i32>,
        cost: i64,
    }
    fn sat32(v: i64) -> i32 {
        v.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
    }
    fn round16(l: i32) -> i32 {
        i32::from((sat32(i64::from(l) + 0x8000) >> 16) as i16)
    }
    fn step(p: &Path, x: i32) -> ([i32; 2], [i32; 2], i32) {
        // acc = L_mult(x, b0) + L_mac(hx0, b1) + L_mac(hx1, b2) + wide feedback.
        let mut acc = sat32(2 * i64::from(x) * i64::from(B[0]));
        acc = sat32(i64::from(acc) + 2 * i64::from(p.hx[0]) * i64::from(B[1]));
        acc = sat32(i64::from(acc) + 2 * i64::from(p.hx[1]) * i64::from(B[2]));
        let fb1 = ((i64::from(p.hy[0]) * i64::from(A[1])) >> 15) << 2;
        let fb2 = ((i64::from(p.hy[1]) * i64::from(A[2])) >> 15) << 2;
        acc = sat32(i64::from(acc) + i64::from(sat32(fb1)));
        acc = sat32(i64::from(acc) + i64::from(sat32(fb2)));
        let y = round16(sat32(i64::from(acc) << 2));
        ([x, p.hx[0]], [acc, p.hy[0]], y)
    }
    let mut paths = vec![Path {
        hx: [0; 2],
        hy: [0; 2],
        xs: Vec::with_capacity(reference.len()),
        cost: 0,
    }];
    // Running gain ratio recovered/guide (Q16), refreshed every
    // subframe from the best path, so the beam is ranked against the
    // guide's SHAPE rather than its level (an offset path is invisible
    // to the double-zero-at-DC filter and must be priced explicitly).
    let mut gain_q16: i64 = 65536;
    for (n, &r) in reference.iter().enumerate() {
        if n > 0 && n % SUBFRAME == 0 {
            let best = paths.iter().min_by_key(|p| p.cost).unwrap();
            let b = n - SUBFRAME;
            let num: i64 = (0..SUBFRAME)
                .map(|i| i64::from(best.xs[b + i]) * i64::from(guide[b + i]))
                .sum();
            let den: i64 = (0..SUBFRAME)
                .map(|i| i64::from(guide[b + i]) * i64::from(guide[b + i]))
                .sum();
            if den > 40 * 100 * 100 {
                gain_q16 = (num * 65536 / den).clamp(32768, 131072);
            }
        }
        let mut next: Vec<Path> = Vec::new();
        for p in &paths {
            // Solve round((2·b0·x + c) << 2) == r for x near the linear estimate.
            let (_, _, y0) = step(p, 0);
            let est = (i64::from(r) - i64::from(y0)) * 65536 / (8 * i64::from(B[0]));
            for x in (est - 3)..=(est + 3) {
                let x = x as i32;
                let (hx, hy, y) = step(p, x);
                if y == i32::from(r) {
                    let mut xs = p.xs.clone();
                    xs.push(x);
                    let target = (gain_q16 * i64::from(guide[n.min(guide.len() - 1)])) >> 16;
                    let cost = p.cost + (i64::from(x) - target).abs();
                    next.push(Path { hx, hy, xs, cost });
                }
            }
        }
        if next.is_empty() {
            let best = paths.into_iter().min_by_key(|p| p.cost).unwrap();
            return (best.xs, Some(n));
        }
        next.sort_by_key(|p| p.cost);
        next.truncate(BEAM);
        paths = next;
    }
    let best = paths.into_iter().min_by_key(|p| p.cost).unwrap();
    (best.xs, None)
}

/// Stage-isolated §4.2 scoring against the reference: inverts the
/// §4.2.5 output stage on the `.PST` (validating the high-pass model
/// on the way — an inconsistent inversion means the eq (91) schedule
/// is wrong at that sample) and scores our AGC output `2·sf′(n)`
/// against the recovered reference sequence — exact share, first
/// divergence, and the per-subframe least-squares gain ratio between
/// the reference's postfiltered signal and ours (the AGC-trajectory
/// instrument: the ratio's per-sample smoothness measures the shape
/// agreement of the pre-AGC cascade independently of the gain).
#[test]
fn fx_full_agc_output_oracle() {
    let Some(root) = conformance_root() else {
        eprintln!("g729 conformance corpus absent — skipping agc-output oracle");
        return;
    };
    let only = std::env::var("G729_FX_ORACLE").ok();
    let mut checked = 0usize;
    for (corpus, name) in ["g729-core", "g729a"]
        .iter()
        .flat_map(|c| CLEAN_VECTORS.iter().map(move |n| (*c, *n)))
    {
        if let Some(o) = &only {
            if o != name && o.as_str() != format!("{corpus}/{name}") {
                continue;
            }
        }
        let label = format!("{corpus}/{name}");
        let bit = std::fs::read(root.join(format!("{corpus}/{name}.BIT"))).unwrap();
        let reference = read_pst(&root.join(format!("{corpus}/{name}.PST")));
        let (agc, tilt) = stage_signals(&bit, corpus == "g729a");
        let n = reference.len().min(agc.len());
        // A saturated reference sample frees the ramp ambiguity of the
        // inversion (any input above the clip point reproduces it), so
        // the recovered sequence is trusted only before the first one.
        let first_sat = reference[..n]
            .iter()
            .position(|&v| v == i16::MAX || v == i16::MIN);
        let n = first_sat.unwrap_or(n);
        // The beam search is O(samples × beam); CI runs the head of
        // every vector, `G729_FX_ORACLE_FULL=1` the whole corpus.
        let n = if std::env::var("G729_FX_ORACLE_FULL").is_ok() {
            n
        } else {
            let cap = std::env::var("G729_FX_ORACLE_CAP")
                .ok()
                .and_then(|c| c.parse().ok())
                .unwrap_or(24_000);
            n.min(cap)
        };
        let (recovered, inconsistent) = recover_agc_output(&reference[..n], &agc[..n]);
        if let Ok(dir) = std::env::var("G729_FX_ORACLE_DUMP") {
            let mut w = String::new();
            for i in 0..recovered.len() {
                use std::fmt::Write as _;
                let _ = writeln!(w, "{},{},{}", recovered[i], agc[i], tilt[i]);
            }
            std::fs::write(
                format!("{dir}/oracle_{}_{name}.csv", corpus.replace('-', "_")),
                w,
            )
            .unwrap();
        }
        let exact = (0..recovered.len())
            .filter(|&i| recovered[i] == i32::from(agc[i]))
            .count();
        let first_div = (0..recovered.len()).find(|&i| recovered[i] != i32::from(agc[i]));
        // Per-subframe least-squares gain ratio reference/ours on the
        // pre-AGC signal (ours: tilt output; theirs: recovered/2/g).
        let mut worst_shape = 0.0f64;
        let mut shape_sum = 0.0f64;
        let mut shape_cnt = 0usize;
        let mut bad_shape = 0usize;
        for k in 0..recovered.len() / SUBFRAME {
            let b = k * SUBFRAME;
            let num: f64 = (0..SUBFRAME)
                .map(|i| f64::from(recovered[b + i]) * f64::from(tilt[b + i]))
                .sum();
            let den: f64 = (0..SUBFRAME)
                .map(|i| 4.0 * f64::from(tilt[b + i]) * f64::from(tilt[b + i]))
                .sum();
            if den < 4.0 * 40.0 * 100.0 * 100.0 {
                continue;
            }
            let gain = 2.0 * num / den;
            // Residual after removing the LS gain, relative.
            let res: f64 = (0..SUBFRAME)
                .map(|i| {
                    let e = f64::from(recovered[b + i]) - gain * 2.0 * f64::from(tilt[b + i]);
                    e * e
                })
                .sum();
            let rel = (res / (num * num / den).max(1.0)).sqrt();
            shape_sum += rel;
            shape_cnt += 1;
            worst_shape = worst_shape.max(rel);
            if rel > 0.1 {
                bad_shape += 1;
            }
        }
        eprintln!(
            "{label}: hp-inversion {}  (trusted range {} samples{})  agc-input exact {}/{} ({:.2}%)  first-div {:?}  shape-residual mean {:.4} worst {:.4} bad(>0.1) {} over {} loud subframes",
            match inconsistent {
                None => "consistent".to_string(),
                Some(i) => format!("INCONSISTENT at {i}"),
            },
            n,
            first_sat.map_or(String::new(), |i| format!(", reference clips at {i}")),
            exact,
            recovered.len(),
            100.0 * exact as f64 / recovered.len().max(1) as f64,
            first_div,
            if shape_cnt > 0 { shape_sum / shape_cnt as f64 } else { 0.0 },
            worst_shape,
            bad_shape,
            shape_cnt
        );
        checked += 1;
        // The eq (91) model is pinned: the inversion must stay
        // consistent over every clean base vector.
        assert!(
            inconsistent.is_none(),
            "{label}: §4.2.5 model inconsistent at sample {inconsistent:?}"
        );
    }
    assert!(checked >= 1);
}

/// Our own (AGC output, tilt output) stage signals for a `.BIT` stream.
fn stage_signals(bit: &[u8], annex_a: bool) -> (Vec<i16>, Vec<i16>) {
    let mut fx = FrameDecoderFx::new();
    let mut pf = if annex_a {
        PostfilterFx::new_annex_a()
    } else {
        PostfilterFx::new()
    };
    if let Ok(spec) = std::env::var("G729_FX_DEC") {
        // Decoder-side (§4.1) latitude: `exc_mode=N,energy_plain=B,
        // code_trunc=B,code_q0=B,recon_ga=N,recon_gb=N,push_ga=N,push_gb=N`.
        let mut grid = oxideav_g729::fx::gains::GainGridFx::default();
        for item in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            let (k, v) = item.split_once('=').expect("field=value");
            let b = v == "1";
            let i: i16 = v.parse().unwrap_or(0);
            match k {
                "exc_mode" => fx.exc_mode = i as u8,
                "energy_plain" => fx.energy_plain = b,
                "code_trunc" => grid.code_trunc = b,
                "code_q0" => grid.code_q0 = b,
                "recon_ga" => grid.recon_ga = i,
                "recon_gb" => grid.recon_gb = i,
                "push_ga" => grid.push_ga = i,
                "push_gb" => grid.push_gb = i,
                other => panic!("unknown decoder latitude field {other}"),
            }
        }
        fx.set_gain_grid(grid);
    }
    if let Ok(spec) = std::env::var("G729_FX_LAT") {
        pf.set_latitude(
            oxideav_g729::fx::postfilter::PfLatitudeFx::default().with_overrides(&spec),
        );
    }
    let mut agc = Vec::new();
    let mut tilt = Vec::new();
    for f in 0..bit.len() / FRAME_BYTES {
        let frame = &bit[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).unwrap();
        let dec = if matches!(kind, FrameKind::Erased) {
            fx.decode_erased_frame()
        } else {
            let params = unpack_parameters(&kind).unwrap();
            fx.decode_frame(&params)
        };
        let int_t1 = usize::try_from(dec.sub[0].t_int.max(1)).unwrap();
        let mut periodic = false;
        for s in 0..2 {
            let anchor = if annex_a {
                usize::try_from(dec.sub[s].t_int.max(1)).unwrap()
            } else {
                int_t1
            };
            let speech: [i16; SUBFRAME] = std::array::from_fn(|n| dec.speech[s * SUBFRAME + n]);
            let t = pf.process_subframe_traced(&speech, &dec.sub[s].a_q12, anchor);
            periodic |= t.decision.gain_q15 > 0;
            agc.extend_from_slice(&t.agc);
            tilt.extend_from_slice(&t.tilt);
        }
        if !matches!(kind, FrameKind::Erased) {
            fx.erasure_periodic = periodic;
        }
    }
    (agc, tilt)
}

/// Whole-vector stage dump for the offline schedule-fitting rig.
/// Gated on `G729_FX_DUMP=<corpus>/<name>:<out.csv>`; a no-op otherwise.
/// One CSV row per sample: `speech,lt,st,tilt,agc,out,ref,gf,gain_in,
/// gain_out,lt_gain,lt_delay,lt_frac,lt_long,t_int,a1..a10`.
#[test]
fn fx_full_stage_dump() {
    let Ok(spec) = std::env::var("G729_FX_DUMP") else {
        return;
    };
    let Some(root) = conformance_root() else {
        return;
    };
    let (vector, out_path) = spec.split_once(':').expect("<vector>:<path>");
    let annex_a = vector.starts_with("g729a");
    let bit = std::fs::read(root.join(format!("{vector}.BIT"))).unwrap();
    let reference = read_pst(&root.join(format!("{vector}.PST")));
    let mut w = String::new();
    let mut fx = FrameDecoderFx::new();
    let mut pf = if annex_a {
        PostfilterFx::new_annex_a()
    } else {
        PostfilterFx::new()
    };
    for f in 0..bit.len() / FRAME_BYTES {
        let frame = &bit[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).unwrap();
        let dec = if matches!(kind, FrameKind::Erased) {
            fx.decode_erased_frame()
        } else {
            let params = unpack_parameters(&kind).unwrap();
            fx.decode_frame(&params)
        };
        let int_t1 = usize::try_from(dec.sub[0].t_int.max(1)).unwrap();
        let mut periodic = false;
        for s in 0..2 {
            let anchor = if annex_a {
                usize::try_from(dec.sub[s].t_int.max(1)).unwrap()
            } else {
                int_t1
            };
            let speech: [i16; SUBFRAME] = std::array::from_fn(|n| dec.speech[s * SUBFRAME + n]);
            let t = pf.process_subframe_traced(&speech, &dec.sub[s].a_q12, anchor);
            periodic |= t.decision.gain_q15 > 0;
            let base = f * SAMPLES_PER_FRAME + s * SUBFRAME;
            for (n, &sp) in speech.iter().enumerate() {
                let r = reference.get(base + n).copied().unwrap_or(0);
                use std::fmt::Write as _;
                let _ = write!(
                    w,
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    sp,
                    t.long_term[n],
                    t.short_term[n],
                    t.tilt[n],
                    t.agc[n],
                    t.output[n],
                    r,
                    t.gf_q12,
                    t.agc_gain_in_q12,
                    t.agc_gain_out_q12,
                    t.decision.gain_q15,
                    t.decision.delay,
                    t.decision.frac,
                    u8::from(t.decision.use_long),
                    dec.sub[s].t_int
                );
                for a in &dec.sub[s].a_q12[1..] {
                    let _ = write!(w, ",{a}");
                }
                w.push('\n');
            }
        }
        if !matches!(kind, FrameKind::Erased) {
            fx.erasure_periodic = periodic;
        }
    }
    std::fs::write(out_path, w).unwrap();
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
        // r461: 0.99996 / 0.99996.
        ("PARITY", [0.9999, 0.9999]),
        ("OVERFLOW", [0.70, -0.30]),
        // ERASURE measured 0.923 / 0.887 (the concealed stretches
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
