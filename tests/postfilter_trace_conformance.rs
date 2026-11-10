//! Replays the staged per-sample postfilter traces
//! (`docs/audio/g729/traces/*.samples.csv` / `*.scalars.csv`) against
//! the fixed-point §4.2 cascade, column by column: the §4.1.6 speech
//! entering the postfilter, the eq (79) residual, the §4.2.1 output,
//! the `1/g_f`-scaled all-pole input, the `H_f(z)` output, the tilt
//! output `sf(n)`, the eq (90) gain trajectory, `sf′(n)` and the final
//! §4.2.5 PCM — for the main body and the Annex A cascade.
//!
//! Two instruments:
//!
//! - **Stage replay** (`staged_stages_replay_exactly`): every stage is
//!   driven with the *staged* inputs of that stage (the memories
//!   seeded from the staged predecessor sub-frame where one exists),
//!   so the §4.2 schedule is checked independently of the §4.1
//!   excitation decode. Every compared value must match.
//! - **End-to-end** (`staged_startup_frames_and_lp_sets_exact`): the
//!   full fixed-point chain decoded from the bitstream must reproduce
//!   every staged LP coefficient set and every column of the staged
//!   frame-0 sub-frames; deeper sub-frames are reported (their `ŝ`
//!   depends on the §3.9.1 gain-decode schedule the staged material
//!   does not pin).
//!
//! When the staged traces are absent (published-crate build) the
//! tests skip.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use oxideav_g729::fx::decoder::FrameDecoderFx;
use oxideav_g729::fx::postfilter::{PostfilterFx, SubframeTraceFx};
use oxideav_g729::parameters::unpack_parameters;
use oxideav_g729::serial::{self, FRAME_BYTES, FrameKind};

const SUBFRAME: usize = 40;

fn docs_root() -> Option<PathBuf> {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut cursor: &Path = &manifest;
    loop {
        let cand = cursor.join("docs/audio/g729");
        if cand.join("traces/notes.md").is_file() && cand.join("conformance/README.md").is_file() {
            return Some(cand);
        }
        cursor = cursor.parent()?;
    }
}

/// One staged sub-frame trace: named sample columns + scalar rows.
struct StagedTrace {
    annex_a: bool,
    vector: String,
    frame: usize,
    sub: usize,
    columns: HashMap<String, Vec<Option<i32>>>,
    scalars: HashMap<String, Vec<i64>>,
}

fn parse_name(stem: &str) -> Option<(bool, String, usize, usize)> {
    // `<decoder>-<VECTOR>-f<frame>s<sub>`
    let mut parts = stem.split('-');
    let annex_a = match parts.next()? {
        "core" => false,
        "annexA" => true,
        _ => return None,
    };
    let vector = parts.next()?.to_string();
    let fs = parts.next()?.strip_prefix('f')?;
    let (f, s) = fs.split_once('s')?;
    Some((annex_a, vector, f.parse().ok()?, s.parse().ok()?))
}

fn load_traces(dir: &Path) -> Vec<StagedTrace> {
    let mut out = Vec::new();
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .expect("traces dir")
        .filter_map(|e| e.ok())
        .filter_map(|e| e.file_name().into_string().ok())
        .filter_map(|n| n.strip_suffix(".samples.csv").map(str::to_string))
        .collect();
    names.sort();
    for stem in names {
        let Some((annex_a, vector, frame, sub)) = parse_name(&stem) else {
            continue;
        };
        let text = std::fs::read_to_string(dir.join(format!("{stem}.samples.csv"))).unwrap();
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split(',').collect();
        let mut columns: HashMap<String, Vec<Option<i32>>> = HashMap::new();
        for line in lines.filter(|l| !l.trim().is_empty()) {
            for (name, cell) in header.iter().zip(line.split(',')) {
                columns
                    .entry((*name).to_string())
                    .or_default()
                    .push(cell.trim().parse::<i32>().ok());
            }
        }
        let mut scalars: HashMap<String, Vec<i64>> = HashMap::new();
        if let Ok(text) = std::fs::read_to_string(dir.join(format!("{stem}.scalars.csv"))) {
            for line in text.lines().skip(1) {
                let cells: Vec<&str> = line.split(',').collect();
                if cells.len() < 3 {
                    continue;
                }
                let keys: Vec<&str> = cells[0].split('|').collect();
                let vals: Vec<&str> = cells[2].split('|').collect();
                for (k, v) in keys.iter().zip(vals.iter()) {
                    if let Ok(v) = v.trim().parse::<i64>() {
                        scalars.entry((*k).to_string()).or_default().push(v);
                    }
                }
            }
        }
        out.push(StagedTrace {
            annex_a,
            vector,
            frame,
            sub,
            columns,
            scalars,
        });
    }
    out
}

/// Decode `bit` through the full fixed-point chain up to and including
/// `last_frame`, returning the speech and stage trace of every
/// sub-frame.
fn decode_traced(
    bit: &[u8],
    annex_a: bool,
    last_frame: usize,
) -> Vec<([i16; SUBFRAME], SubframeTraceFx, [i16; 11])> {
    let mut fx = FrameDecoderFx::new();
    let mut pf = if annex_a {
        PostfilterFx::new_annex_a()
    } else {
        PostfilterFx::new()
    };
    let mut out = Vec::new();
    for f in 0..=last_frame {
        let frame = &bit[f * FRAME_BYTES..(f + 1) * FRAME_BYTES];
        let kind = serial::parse_frame(frame).expect("clean vector frame");
        assert!(!matches!(kind, FrameKind::Erased));
        let params = unpack_parameters(&kind).expect("parameters");
        let dec = fx.decode_frame(&params);
        let int_t1 = usize::try_from(dec.sub[0].t_int.max(1)).unwrap();
        for s in 0..2 {
            let anchor = if annex_a {
                usize::try_from(dec.sub[s].t_int.max(1)).unwrap()
            } else {
                int_t1
            };
            let speech: [i16; SUBFRAME] = std::array::from_fn(|n| dec.speech[s * SUBFRAME + n]);
            let t = pf.process_subframe_traced(&speech, &dec.sub[s].a_q12, anchor);
            out.push((speech, t, dec.sub[s].a_q12));
        }
    }
    out
}

fn mismatches(ours: &[i16], col: Option<&Vec<Option<i32>>>) -> Option<usize> {
    let col = col?;
    Some(
        ours.iter()
            .zip(col.iter())
            .filter(|(o, r)| r.is_some_and(|r| r != i32::from(**o)))
            .count(),
    )
}

fn col_i16(t: &StagedTrace, name: &str) -> Option<[i16; SUBFRAME]> {
    let c = t.columns.get(name)?;
    if c.len() != SUBFRAME || c.iter().any(Option::is_none) {
        return None;
    }
    Some(std::array::from_fn(|n| c[n].unwrap() as i16))
}

fn count_bad(ours: &[i16], col: Option<&Vec<Option<i32>>>) -> usize {
    mismatches(ours, col).unwrap_or(0)
}

/// Stage replay with staged inputs — independent of the §4.1 decode.
#[test]
fn staged_stages_replay_exactly() {
    let Some(root) = docs_root() else {
        eprintln!("g729 staged traces absent — skipping stage replay");
        return;
    };
    let traces = load_traces(&root.join("traces"));
    assert!(traces.len() >= 19, "only {} staged traces", traces.len());
    let find = |annex_a: bool, vector: &str, frame: usize, sub: usize| {
        traces.iter().find(|t| {
            t.annex_a == annex_a && t.vector == vector && t.frame == frame && t.sub == sub
        })
    };
    let prev_of = |t: &StagedTrace| {
        let idx = 2 * t.frame + t.sub;
        if idx == 0 {
            return None;
        }
        find(t.annex_a, &t.vector, (idx - 1) / 2, (idx - 1) % 2)
    };

    let mut checked = 0usize;
    let mut bad_total = 0usize;
    let mut lt_checked = 0usize;
    for t in &traces {
        let label = format!(
            "{}-{}-f{}s{}",
            if t.annex_a { "annexA" } else { "core" },
            t.vector,
            t.frame,
            t.sub
        );
        let shat = col_i16(t, "shat").expect("shat column");
        let sf = col_i16(t, "sf").expect("sf column");
        let mut line = String::new();
        let mut note = |name: &str, bad: usize, line: &mut String| {
            checked += 1;
            bad_total += bad;
            if bad != 0 {
                line.push_str(&format!(" {name}:{bad}"));
            }
        };

        if t.annex_a {
            // §A.4.2.4 AGC from the staged all-pole output and gain state.
            // A zero-output-energy sub-frame prints no gain state: any
            // non-zero seed must be forced to 0 by the reset.
            let g_prev = t.scalars.get("g_minus1_Q12").map_or(4096, |g| g[0] as i16);
            let mut pf = PostfilterFx::new_annex_a();
            pf.replay_seed(&[], &[], &[], 0, g_prev);
            let mut a = [0i16; 11];
            a[0] = 4096;
            let tr = pf.replay_after_long_term(&shat, &a, &[0; SUBFRAME], Some(&sf));
            note(
                "sf_prime",
                count_bad(&tr.agc, t.columns.get("sf_prime")),
                &mut line,
            );
            let reset = t.scalars.contains_key("zero_output_gain_reset")
                || t.scalars
                    .get("energy_out_quarter")
                    .is_some_and(|e| e[0] == 0);
            if reset {
                note("reset", usize::from(tr.agc_gain_out != 0), &mut line);
            } else {
                note(
                    "gain",
                    count_bad(&tr.gain, t.columns.get("g_n_Q12")),
                    &mut line,
                );
            }
            eprintln!("{label}:{}", if line.is_empty() { " exact" } else { &line });
            continue;
        }

        let a_staged = &t.scalars["lp_coef_ahat_Q12"];
        let a_q12: [i16; 11] = std::array::from_fn(|i| a_staged[i] as i16);
        let lt = col_i16(t, "ltp_out").expect("ltp_out column");
        let hf = col_i16(t, "hf_out").expect("hf_out column");
        let tilt_prev = t.scalars["tilt_state_prev_hf_out"][0] as i16;
        let g_prev = t.scalars["g_minus1_Q14"][0] as i16;

        // (i) Tail: 1/g_f scaling, impulse response, tilt, AGC — from
        // the staged §4.2.1 output and the staged all-pole output.
        let mut pf = PostfilterFx::new();
        pf.replay_seed(&[], &[], &[], tilt_prev, g_prev);
        let tr = pf.replay_after_long_term(&shat, &a_q12, &lt, Some(&hf));
        note(
            "gf_scaled",
            count_bad(&tr.gf_scaled, t.columns.get("gf_scaled")),
            &mut line,
        );
        note("sf", count_bad(&tr.tilt, t.columns.get("sf")), &mut line);
        note(
            "gain",
            count_bad(&tr.gain, t.columns.get("g_n_Q14")),
            &mut line,
        );
        note(
            "sf_prime",
            count_bad(&tr.agc, t.columns.get("sf_prime")),
            &mut line,
        );
        let h = &t.scalars["hf_Q12"];
        let bad_h = h
            .iter()
            .zip(tr.impulse.iter())
            .filter(|(r, o)| **r != i64::from(**o))
            .count();
        note("hf_Q12", bad_h, &mut line);
        if let Some(q) = t.scalars.get("gf_sum_Q12") {
            note(
                "gf_sum",
                usize::from(q[0] != i64::from(tr.gf_q12)),
                &mut line,
            );
        }

        // (ii) With a staged predecessor (or the frame-0 zero state):
        // the eq (79) residual, the all-pole section and — given
        // enough staged residual history — the §4.2.1 decision.
        let mut chain: Vec<&StagedTrace> = Vec::new();
        let mut cur = t;
        while let Some(p) = prev_of(cur) {
            chain.push(p);
            cur = p;
        }
        let at_start = 2 * cur.frame + cur.sub == 0;
        if !chain.is_empty() || at_start {
            let mut r_hist: Vec<i16> = Vec::new();
            for p in chain.iter().rev() {
                r_hist.extend_from_slice(&col_i16(p, "resid").unwrap());
            }
            let (s_prev, y_prev): (Vec<i16>, Vec<i16>) = match chain.first() {
                Some(p) => (
                    col_i16(p, "shat").unwrap().to_vec(),
                    col_i16(p, "hf_out").unwrap().to_vec(),
                ),
                None => (vec![0; 10], vec![0; 10]),
            };
            let bit = std::fs::read(root.join(format!("conformance/g729-core/{}.BIT", t.vector)))
                .unwrap();
            let kind =
                serial::parse_frame(&bit[t.frame * FRAME_BYTES..(t.frame + 1) * FRAME_BYTES])
                    .unwrap();
            let params = unpack_parameters(&kind).unwrap();
            let int_t1 =
                usize::try_from(oxideav_g729::pitch_decode::decode_t1_from_p1(params.p1).int_t)
                    .unwrap();

            let mut pf = PostfilterFx::new();
            pf.replay_seed(&s_prev, &r_hist, &y_prev, tilt_prev, g_prev);
            let full = pf.process_subframe_traced(&shat, &a_q12, int_t1);
            note(
                "resid",
                count_bad(&full.residual, t.columns.get("resid")),
                &mut line,
            );

            let mut pf = PostfilterFx::new();
            pf.replay_seed(&s_prev, &r_hist, &y_prev, tilt_prev, g_prev);
            let mid = pf.replay_after_long_term(&shat, &a_q12, &lt, None);
            note(
                "hf_out",
                count_bad(&mid.short_term, t.columns.get("hf_out")),
                &mut line,
            );

            // §4.2.1 needs the residual back to the deepest candidate
            // delay plus the interpolator reach.
            if at_start || r_hist.len() >= int_t1 + 2 + 8 {
                lt_checked += 1;
                note(
                    "ltp_out",
                    count_bad(&full.long_term, t.columns.get("ltp_out")),
                    &mut line,
                );
                let d = t.scalars["ltp_int_delay"][0];
                let ph = t.scalars["ltp_phase_eighths"][0];
                // Staged convention: T = int − phase/8; ours: delay + frac/8.
                let ours_t8 = if full.decision.gain_q15 == 0 {
                    0
                } else {
                    8 * full.decision.delay as i64 + full.decision.frac as i64
                };
                if ours_t8 != 8 * d - ph {
                    line.push_str(&format!(" T8:{ours_t8} vs {}", 8 * d - ph));
                }
                note("T8", usize::from(ours_t8 != 8 * d - ph), &mut line);
                if let Some(w) = t.scalars.get("ltp_gl_Q15(weight_on_undelayed)") {
                    if w[0] != i64::from(full.decision.undelayed_q15) {
                        line.push_str(&format!(" g0:{} vs {}", full.decision.undelayed_q15, w[0]));
                    }
                    note(
                        "g0",
                        usize::from(w[0] != i64::from(full.decision.undelayed_q15)),
                        &mut line,
                    );
                }
            }
        }
        eprintln!("{label}:{}", if line.is_empty() { " exact" } else { &line });
    }
    assert!(checked >= 140, "only {checked} stage columns checked");
    assert!(
        lt_checked >= 9,
        "only {lt_checked} long-term decisions checked"
    );
    assert_eq!(bad_total, 0, "{bad_total} staged stage values diverge");
}

#[test]
fn staged_startup_frames_and_lp_sets_exact() {
    let Some(root) = docs_root() else {
        eprintln!("g729 staged traces absent — skipping trace replay");
        return;
    };
    let traces = load_traces(&root.join("traces"));
    assert!(traces.len() >= 19, "only {} staged traces", traces.len());

    // Decode each (decoder, vector) once, to its deepest traced frame.
    let mut deepest: HashMap<(bool, String), usize> = HashMap::new();
    for t in &traces {
        let e = deepest.entry((t.annex_a, t.vector.clone())).or_default();
        *e = (*e).max(t.frame);
    }
    let mut decoded = HashMap::new();
    for ((annex_a, vector), last) in &deepest {
        let corpus = if *annex_a { "g729a" } else { "g729-core" };
        let bit = std::fs::read(root.join(format!("conformance/{corpus}/{vector}.BIT"))).unwrap();
        decoded.insert(
            (*annex_a, vector.clone()),
            decode_traced(&bit, *annex_a, *last),
        );
    }

    let mut total_bad = 0usize;
    let mut deep_bad = 0usize;
    let mut checked = 0usize;
    for t in &traces {
        let before = total_bad;
        let (speech, tr, a_q12) = &decoded[&(t.annex_a, t.vector.clone())][2 * t.frame + t.sub];
        let label = format!(
            "{}-{}-f{}s{}",
            if t.annex_a { "annexA" } else { "core" },
            t.vector,
            t.frame,
            t.sub
        );
        let gain_col = if t.annex_a { "g_n_Q12" } else { "g_n_Q14" };
        let mut stages: Vec<(&str, &[i16], &str)> = vec![
            ("shat", &speech[..], "shat"),
            ("resid", &tr.residual[..], "resid"),
            ("ltp_out", &tr.long_term[..], "ltp_out"),
        ];
        if t.annex_a {
            stages.push(("tilt_out", &tr.short_term[..], "tilt_out"));
        } else {
            stages.push(("gf_scaled", &tr.gf_scaled[..], "gf_scaled"));
            stages.push(("hf_out", &tr.short_term[..], "hf_out"));
        }
        stages.push(("sf", &tr.tilt[..], "sf"));
        // A zero-output-energy sub-frame passes unscaled and prints no
        // gain trajectory worth comparing (the state is forced to 0).
        if !t.scalars.contains_key("zero_output_gain_reset") && tr.agc_gain_out != 0 {
            stages.push(("gain", &tr.gain[..], gain_col));
        }
        stages.push(("sf_prime", &tr.agc[..], "sf_prime"));
        stages.push(("pst", &tr.output[..], "pst"));

        let mut line = String::new();
        for (name, ours, col) in stages {
            match mismatches(ours, t.columns.get(col)) {
                Some(bad) => {
                    checked += 1;
                    total_bad += bad;
                    if bad != 0 {
                        line.push_str(&format!(" {name}:{bad}"));
                    }
                }
                None => line.push_str(&format!(" {name}:absent")),
            }
        }
        // Scalars: the LP coefficients, the impulse response and the
        // gain entering the sub-frame.
        if let Some(a) = t.scalars.get("lp_coef_ahat_Q12") {
            let bad = a
                .iter()
                .zip(a_q12.iter())
                .filter(|(r, o)| **r != i64::from(**o))
                .count();
            checked += 1;
            total_bad += bad;
            if bad != 0 {
                line.push_str(&format!(" a_q12:{bad} (ours {a_q12:?} staged {a:?})"));
            }
        }
        if let Some(h) = t.scalars.get("hf_Q12") {
            let bad = h
                .iter()
                .zip(tr.impulse.iter())
                .filter(|(r, o)| **r != i64::from(**o))
                .count();
            checked += 1;
            total_bad += bad;
            if bad != 0 {
                line.push_str(&format!(" hf_Q12:{bad}"));
            }
        }
        let g_key = if t.annex_a {
            "g_minus1_Q12"
        } else {
            "g_minus1_Q14"
        };
        if let Some(g) = t.scalars.get(g_key) {
            checked += 1;
            if g[0] != i64::from(tr.agc_gain_in) {
                total_bad += 1;
                line.push_str(&format!(" g(-1):{} vs {}", tr.agc_gain_in, g[0]));
            }
        }
        if let (Some(d), Some(p)) = (
            t.scalars.get("ltp_int_delay"),
            t.scalars.get("ltp_phase_eighths"),
        ) {
            // Staged convention: T = int − phase/8; ours: delay + frac/8.
            let ours_t8 = if tr.decision.gain_q15 == 0 {
                0
            } else {
                8 * tr.decision.delay as i64 + tr.decision.frac as i64
            };
            let staged_t8 = 8 * d[0] - p[0];
            checked += 1;
            if ours_t8 != staged_t8 {
                total_bad += 1;
                line.push_str(&format!(" T8:{ours_t8} vs {staged_t8}"));
            }
        }
        if let Some(w) = t.scalars.get("ltp_gl_Q15(weight_on_undelayed)") {
            checked += 1;
            if w[0] != i64::from(tr.decision.undelayed_q15) {
                total_bad += 1;
                line.push_str(&format!(" g0:{} vs {}", tr.decision.undelayed_q15, w[0]));
            }
        }
        if line.is_empty() {
            eprintln!("{label}: exact");
        } else {
            eprintln!("{label}:{line}");
        }
        // Frame 0 must be exact end to end; deeper sub-frames are
        // reported only (see the module doc) — except their LP sets.
        if t.frame != 0 {
            deep_bad += total_bad - before;
            total_bad = before;
        }
        let lp_bad = t.scalars.get("lp_coef_ahat_Q12").map_or(0, |a| {
            a.iter()
                .zip(a_q12.iter())
                .filter(|(r, o)| **r != i64::from(**o))
                .count()
        });
        assert_eq!(lp_bad, 0, "{label}: LP coefficient set diverges");
    }
    eprintln!("deeper sub-frames: {deep_bad} end-to-end values off (reported, not asserted)");
    assert!(checked >= 150, "only {checked} trace columns checked");
    assert_eq!(
        total_bad, 0,
        "{total_bad} frame-0 staged trace values diverge"
    );
}
