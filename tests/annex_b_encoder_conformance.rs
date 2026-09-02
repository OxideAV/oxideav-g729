//! Annex B **encoder-side** conformance: the §B.4.1 DTX decision and
//! the §B.4.2 SID quantisation (`oxideav_g729::annex_b_encoder`)
//! against the staged `g729b/tstseq1–4` `.bin` → `.bit` pairs under a
//! **locked VAD** — the voice-activity flag of every frame is read off
//! the reference bitstream's own frame type (active vs SID /
//! untransmitted), since the §B.3 VAD is not implementable from the
//! prose (see the module docs). What is measured:
//!
//! * on inactive frames, whether our DTX emits a SID frame exactly
//!   where the reference does (`Ftyp` agreement);
//! * on frames where both emit SID, the agreement of the four Table
//!   B.2 indices and the gain-index distance;
//! * the eq (B.15) energy-scale pin (the SID energy is on the
//!   un-halved input scale).
//!
//! Skips clean when the corpus is absent.

use std::path::{Path, PathBuf};

use oxideav_g729::annex_b::{parse_annex_b_stream, AnnexBFrame};
use oxideav_g729::annex_b_encoder::{
    absolute_autocorrelation, DtxDecision, DtxEncoder, DtxLatitude,
};
use oxideav_g729::fx::analysis::{analyze_window_fx, FrontEndLatitude, PreprocessorFx};
use oxideav_g729::tables::L_WINDOW;

const SAMPLES_PER_FRAME: usize = 80;

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

fn read_pcm(path: &Path) -> Vec<i16> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    bytes
        .chunks_exact(2)
        .map(|c| i16::from_le_bytes([c[0], c[1]]))
        .collect()
}

#[derive(Default, Debug)]
struct Tally {
    inactive: usize,
    ref_sid: usize,
    our_sid: usize,
    ftyp_agree: usize,
    both_sid: usize,
    lp0: usize,
    l1: usize,
    l2: usize,
    gain: usize,
    gain_within1: usize,
    all4: usize,
    /// Our `E'` (dB) minus the reference SID's decoded energy, per
    /// reference SID frame.
    e_offsets: Vec<f64>,
}

fn run(root: &Path, name: &str, raw_input: bool, lat: DtxLatitude) -> Tally {
    let samples = read_pcm(&root.join(format!("g729b/{name}.bin")));
    let bit = std::fs::read(root.join(format!("g729b/{name}.bit"))).unwrap();
    let frames = parse_annex_b_stream(&bit).expect("reference stream parses");
    let n_frames = (samples.len() / SAMPLES_PER_FRAME).min(frames.len());

    let mut preproc = PreprocessorFx::new();
    let mut speech = [0i16; L_WINDOW];
    let mut dtx = DtxEncoder::with_latitude(lat);
    let mut t = Tally::default();
    let mut prev_active = true;
    for f in 0..n_frames {
        let mut frame = [0i16; SAMPLES_PER_FRAME];
        frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
        let s_new = if raw_input {
            frame
        } else {
            preproc.process_frame(&frame)
        };
        speech.copy_within(SAMPLES_PER_FRAME.., 0);
        speech[L_WINDOW - SAMPLES_PER_FRAME..].copy_from_slice(&s_new);
        let ac = analyze_window_fx(&speech, &FrontEndLatitude::default());
        let r = absolute_autocorrelation(&ac);

        let vad = match &frames[f] {
            AnnexBFrame::Active(_) => true,
            AnnexBFrame::Sid(_) | AnnexBFrame::Untransmitted => false,
            AnnexBFrame::Erased => prev_active,
        };
        prev_active = vad;
        let reference = match &frames[f] {
            AnnexBFrame::Sid(sid) => Some(sid),
            _ => None,
        };
        let ours = if vad {
            dtx.process_frame(vad, &r)
        } else {
            dtx.process_frame_locked(vad, &r, reference)
        };
        if vad {
            continue;
        }
        if let AnnexBFrame::Sid(rs) = &frames[f] {
            t.e_offsets
                .push(dtx.last_e_prime_db - f64::from(rs.energy_db()));
        }
        t.inactive += 1;
        let ref_sid = matches!(frames[f], AnnexBFrame::Sid(_));
        let our_sid = matches!(ours, DtxDecision::Sid(_));
        t.ref_sid += usize::from(ref_sid);
        t.our_sid += usize::from(our_sid);
        t.ftyp_agree += usize::from(ref_sid == our_sid);
        if let (AnnexBFrame::Sid(r), DtxDecision::Sid(o)) = (&frames[f], ours) {
            t.both_sid += 1;
            t.lp0 += usize::from(o.lp0 == r.lp0);
            t.l1 += usize::from(o.l1 == r.l1);
            t.l2 += usize::from(o.l2 == r.l2);
            t.gain += usize::from(o.gain == r.gain);
            t.gain_within1 += usize::from((i32::from(o.gain) - i32::from(r.gain)).abs() <= 1);
            t.all4 += usize::from(o == *r);
        }
    }
    t
}

/// Locked-VAD DTX / SID agreement with regression floors.
///
/// Measured (round 455, default latitude): Ftyp agreement 86.5 / 69.9
/// / 85.7 / 77.1 %, SID gain index exact 89.5 / 100 / 100 / 83.7 %
/// (±1: 100 % everywhere), lp0 70 / 60 / 53 / 49 %, l1 63 / 20 / 47 /
/// 26 %, l2 62 / 20 / 34 / 35 %, all four 42 / 0 / 17 / 16 %
/// (tstseq1–4; tstseq2 has only 13 reference SID frames).
#[test]
fn locked_vad_dtx_sid_agreement() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let pct = |h: usize, n: usize| 100.0 * h as f64 / n.max(1) as f64;
    // (sequence, Ftyp floor %, gain-exact floor %, gain-within-1 floor %)
    let floors: [(&str, f64, f64, f64); 4] = [
        ("tstseq1", 83.0, 85.0, 97.0),
        ("tstseq2", 66.0, 80.0, 97.0),
        ("tstseq3", 82.0, 95.0, 97.0),
        ("tstseq4", 73.0, 79.0, 97.0),
    ];
    let mut total = Tally::default();
    for (name, f_ftyp, f_gain, f_gain1) in floors {
        let t = run(&root, name, false, DtxLatitude::default());
        println!(
            "{name}: inactive={} refSID={} ourSID={} Ftyp={:.1}% | bothSID={} lp0={:.1} l1={:.1} l2={:.1} gain={:.1} gain±1={:.1} all4={:.1}",
            t.inactive,
            t.ref_sid,
            t.our_sid,
            pct(t.ftyp_agree, t.inactive),
            t.both_sid,
            pct(t.lp0, t.both_sid),
            pct(t.l1, t.both_sid),
            pct(t.l2, t.both_sid),
            pct(t.gain, t.both_sid),
            pct(t.gain_within1, t.both_sid),
            pct(t.all4, t.both_sid),
        );
        assert!(
            t.inactive > 0,
            "{name}: no inactive frames under the reference VAD"
        );
        let mut off = t.e_offsets.clone();
        off.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if !off.is_empty() {
            println!(
                "{name}: E' offset vs reference SID energy: median {:.2} dB, p10 {:.2}, p90 {:.2}",
                off[off.len() / 2],
                off[off.len() / 10],
                off[off.len() * 9 / 10]
            );
        }
        assert!(
            pct(t.ftyp_agree, t.inactive) >= f_ftyp,
            "{name}: Ftyp under floor"
        );
        assert!(
            pct(t.gain, t.both_sid) >= f_gain,
            "{name}: SID gain under floor"
        );
        assert!(
            pct(t.gain_within1, t.both_sid) >= f_gain1,
            "{name}: SID gain ±1 under floor"
        );
        total.inactive += t.inactive;
        total.ref_sid += t.ref_sid;
        total.our_sid += t.our_sid;
        total.ftyp_agree += t.ftyp_agree;
        total.both_sid += t.both_sid;
        total.gain_within1 += t.gain_within1;
        total.gain += t.gain;
    }
    println!(
        "ALL: inactive={} refSID={} ourSID={} Ftyp={:.1}% bothSID={} gain={:.1} gain±1={:.1}",
        total.inactive,
        total.ref_sid,
        total.our_sid,
        pct(total.ftyp_agree, total.inactive),
        total.both_sid,
        pct(total.gain, total.both_sid),
        pct(total.gain_within1, total.both_sid),
    );
}

/// The eq (B.15) energy scale is the un-halved input scale: computed
/// on the §3.1 pre-processed signal without the pinned ×4 the SID
/// energies sit ≈ 6 dB under the reference on every sequence with
/// reference SID frames below the ladder ceiling.
#[test]
fn energy_scale_pin() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let literal = DtxLatitude {
        energy_scale_x100: 100,
        ..Default::default()
    };
    for name in ["tstseq1", "tstseq3", "tstseq4"] {
        let mut lit = run(&root, name, false, literal).e_offsets;
        let mut pin = run(&root, name, false, DtxLatitude::default()).e_offsets;
        lit.sort_by(|a, b| a.partial_cmp(b).unwrap());
        pin.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let m_lit = lit[lit.len() / 2];
        let m_pin = pin[pin.len() / 2];
        println!("{name}: median E' offset literal {m_lit:.2} dB, pinned {m_pin:.2} dB");
        assert!(
            m_lit < -4.5,
            "{name}: literal reading should sit ≈ 6 dB low"
        );
        assert!(
            m_pin.abs() < 1.5,
            "{name}: pinned scale should centre on the reference"
        );
    }
}
