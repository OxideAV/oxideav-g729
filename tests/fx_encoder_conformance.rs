//! Fixed-point encoder conformance: the clause-3 `FrameEncoderFx`
//! against the staged ITU `.IN` → `.BIT` vectors, measured two ways.
//!
//! * **Locked** — every frame's searches run from the state the
//!   reference encoder had (the local decoder state is committed from
//!   the reference's own transmitted parameters, see
//!   `fx::encoder`), so each stage's per-parameter agreement is
//!   measured in isolation from upstream error propagation.
//! * **Free** — the ordinary end-to-end encode, whose agreement is the
//!   wire-exactness ratchet.
//!
//! When the corpus is absent (published-crate build) the tests skip.

use std::path::{Path, PathBuf};

use oxideav_g729::fx::encoder::FrameEncoderFx;
use oxideav_g729::parameters::{unpack_parameters, Parameters};
use oxideav_g729::pitch_decode::decode_t1_from_p1;
use oxideav_g729::serial::{parse_frame, FrameKind, FRAME_BYTES};

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

fn read_pcm(path: &Path) -> Vec<i16> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    bytes
        .chunks_exact(2)
        .map(|c| i16::from_le_bytes([c[0], c[1]]))
        .collect()
}

/// Per-parameter hit counts over the active frames of one vector.
#[derive(Debug, Default, Clone, Copy)]
pub struct Agreement {
    pub frames: usize,
    pub lsp_all: usize,
    pub l1: usize,
    pub t1: usize,
    pub t1_int: usize,
    pub t2: usize,
    pub c1: usize,
    pub s1: usize,
    pub c2: usize,
    pub s2: usize,
    pub ga: usize,
    pub gb: usize,
    /// All 15 codewords of the frame agree.
    pub frame_exact: usize,
    /// Reference `int(T1)` outside our §3.7 subframe-1 window (an
    /// open-loop miss rather than a closed-loop one).
    pub t1_window_miss: usize,
    /// `C1` agreement conditional on `T1` agreement (hits, total).
    pub c1_given_t1: (usize, usize),
    /// `GA1`/`GB1` agreement conditional on `T1` + `C1` + `S1`.
    pub g1_given_t1c1: (usize, usize),
}

impl Agreement {
    fn tally(&mut self, o: &Parameters, r: &Parameters, t_op: i32) {
        self.frames += 1;
        let (lo, hi) = oxideav_g729::fx::pitch_cl::subframe1_window(t_op);
        let rt1 = decode_t1_from_p1(r.p1).int_t;
        self.t1_window_miss += usize::from(rt1 < lo || rt1 > hi);
        if o.p1 == r.p1 {
            self.c1_given_t1.1 += 1;
            self.c1_given_t1.0 += usize::from(o.c1 == r.c1);
            if o.c1 == r.c1 && o.s1 == r.s1 {
                self.g1_given_t1c1.1 += 1;
                self.g1_given_t1c1.0 += usize::from(o.ga1 == r.ga1 && o.gb1 == r.gb1);
            }
        }
        self.lsp_all += usize::from(o.l0 == r.l0 && o.l1 == r.l1 && o.l2 == r.l2 && o.l3 == r.l3);
        self.l1 += usize::from(o.l1 == r.l1);
        self.t1 += usize::from(o.p1 == r.p1);
        self.t1_int += usize::from(decode_t1_from_p1(o.p1).int_t == decode_t1_from_p1(r.p1).int_t);
        self.t2 += usize::from(o.p2 == r.p2);
        self.c1 += usize::from(o.c1 == r.c1);
        self.s1 += usize::from(o.s1 == r.s1);
        self.c2 += usize::from(o.c2 == r.c2);
        self.s2 += usize::from(o.s2 == r.s2);
        self.ga += usize::from(o.ga1 == r.ga1) + usize::from(o.ga2 == r.ga2);
        self.gb += usize::from(o.gb1 == r.gb1) + usize::from(o.gb2 == r.gb2);
        self.frame_exact += usize::from(o == r);
    }

    fn pct(hits: usize, total: usize) -> f64 {
        100.0 * hits as f64 / total.max(1) as f64
    }

    fn line(&self, name: &str, mode: &str) -> String {
        let n = self.frames;
        format!(
            "{name:6} {mode:6} n={n:4} LSP={:5.1} T1={:5.1} (int {:5.1}) T2={:5.1} C1={:5.1} S1={:5.1} C2={:5.1} S2={:5.1} GA={:5.1} GB={:5.1} frame={:5.1}",
            Self::pct(self.lsp_all, n),
            Self::pct(self.t1, n),
            Self::pct(self.t1_int, n),
            Self::pct(self.t2, n),
            Self::pct(self.c1, n),
            Self::pct(self.s1, n),
            Self::pct(self.c2, n),
            Self::pct(self.s2, n),
            Self::pct(self.ga, 2 * n),
            Self::pct(self.gb, 2 * n),
            Self::pct(self.frame_exact, n),
        ) + &format!(
            " | OLmiss={:4.1} C1|T1={:5.1} G1|T1C1={:5.1}",
            Self::pct(self.t1_window_miss, n),
            Self::pct(self.c1_given_t1.0, self.c1_given_t1.1),
            Self::pct(self.g1_given_t1c1.0, self.g1_given_t1c1.1),
        )
    }
}

/// Runs one vector in the requested mode.
fn run(root: &Path, name: &str, locked: bool) -> Agreement {
    let samples = read_pcm(&root.join(format!("g729-core/{name}.IN")));
    let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
    let n_frames = (samples.len() / SAMPLES_PER_FRAME).min(bit_bytes.len() / FRAME_BYTES);
    let mut enc = FrameEncoderFx::new();
    let mut agg = Agreement::default();
    for f in 0..n_frames {
        let mut frame = [0i16; SAMPLES_PER_FRAME];
        frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
        let rf = parse_frame(&bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES]).unwrap();
        let rp = match rf {
            FrameKind::Active(_) => Some(unpack_parameters(&rf).unwrap()),
            _ => None,
        };
        let out = match (locked, &rp) {
            (true, Some(r)) => enc.encode_frame_locked(&frame, r),
            _ => enc.encode_frame(&frame),
        };
        if let Some(r) = rp {
            agg.tally(&out.params, &r, out.t_op);
        }
    }
    agg
}

/// Whole-corpus per-parameter agreement in both modes, with
/// regression floors on the locked and free rates.
#[test]
fn fx_encoder_agreement() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    // (vector, locked T1 %, locked C1 %, locked GA %, free frame-exact %)
    //
    // Measured after the §3.3 migration (table-log2 LAR decision with
    // the eq (28) logarithm read as base 10, γ₂ from the Q13 d_min):
    // locked T1 71.4/34.2/81.5/81.9/69.9/84.4, C1 68.6/70.8/39.3/90.2/
    // 68.8/55.5, GA 81.4/97.5/69.5/90.1/83.1/60.2 (ALL: T1 75.5, T2 76.3,
    // C1 65.4, GA 80.8, frame-exact 24.0); free frame-exact ≈ 0.
    let floors: [(&str, f64, f64, f64, f64); 6] = [
        ("ALGTHM", 68.0, 65.0, 78.0, 0.0),
        ("FIXED", 31.0, 67.0, 94.0, 0.0),
        ("LSP", 78.0, 36.0, 66.0, 0.0),
        ("PITCH", 78.0, 87.0, 87.0, 0.0),
        ("SPEECH", 66.0, 65.0, 80.0, 0.0),
        ("TAME", 81.0, 52.0, 57.0, 0.0),
    ];
    let mut total_locked = Agreement::default();
    let mut total_free = Agreement::default();
    for (name, f_t1, f_c1, f_ga, f_frame) in floors {
        let locked = run(&root, name, true);
        let free = run(&root, name, false);
        println!("{}", locked.line(name, "locked"));
        println!("{}", free.line(name, "free"));
        let n = locked.frames;
        assert!(
            Agreement::pct(locked.t1, n) >= f_t1,
            "{name}: locked T1 under floor"
        );
        assert!(
            Agreement::pct(locked.c1, n) >= f_c1,
            "{name}: locked C1 under floor"
        );
        assert!(
            Agreement::pct(locked.ga, 2 * n) >= f_ga,
            "{name}: locked GA under floor"
        );
        assert!(
            Agreement::pct(free.frame_exact, n) >= f_frame,
            "{name}: free frame-exact under floor"
        );
        for (t, s) in [(&mut total_locked, locked), (&mut total_free, free)] {
            t.frames += s.frames;
            t.lsp_all += s.lsp_all;
            t.l1 += s.l1;
            t.t1 += s.t1;
            t.t1_int += s.t1_int;
            t.t2 += s.t2;
            t.c1 += s.c1;
            t.s1 += s.s1;
            t.c2 += s.c2;
            t.s2 += s.s2;
            t.ga += s.ga;
            t.gb += s.gb;
            t.frame_exact += s.frame_exact;
            t.t1_window_miss += s.t1_window_miss;
            t.c1_given_t1.0 += s.c1_given_t1.0;
            t.c1_given_t1.1 += s.c1_given_t1.1;
            t.g1_given_t1c1.0 += s.g1_given_t1c1.0;
            t.g1_given_t1c1.1 += s.g1_given_t1c1.1;
        }
    }
    println!("{}", total_locked.line("ALL", "locked"));
    println!("{}", total_free.line("ALL", "free"));
}

/// The fixed encoder's stream decodes cleanly through the fixed
/// decoder and its `.IN` → serial path round-trips.
#[test]
fn fx_encoder_stream_is_well_formed() {
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    let samples = read_pcm(&root.join("g729-core/ALGTHM.IN"));
    let mut enc = FrameEncoderFx::new();
    let mut enc2 = FrameEncoderFx::new();
    let mut dec = oxideav_g729::fx::decoder::FrameDecoderFx::new();
    for f in 0..samples.len() / SAMPLES_PER_FRAME {
        let mut frame = [0i16; SAMPLES_PER_FRAME];
        frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
        let wire = enc.encode_frame_to_serial(&frame);
        let expect = enc2.encode_frame(&frame).params;
        let parsed = parse_frame(&wire).unwrap();
        let got = unpack_parameters(&parsed).unwrap();
        assert_eq!(got, expect, "frame {f}: serial round-trip mismatch");
        assert!(got.pitch_parity_ok(), "frame {f}: parity");
        let d = dec.decode_frame(&got);
        assert!(d.speech.iter().all(|&s| s != i16::MIN));
    }
}

/// Stage-isolation diagnostic for the §3.8.1 search: on every locked
/// subframe where our pulse set differs from the reference's, score
/// both sets under OUR `d(n)` / `φ` (eq (53) `C²/E`). A reference set
/// that scores higher under our own criterion is a search-structure
/// miss (threshold / budget / loop order); one that scores lower means
/// the inputs (`x′`, `h`) differ from the reference's.
#[test]
fn acelp_stage_isolation_diagnostic() {
    use oxideav_g729::fixed_codebook::decode_pulses;
    use oxideav_g729::fixed_codebook_search::{correlation_d, phi_matrix};
    let Some(root) = conformance_root() else {
        eprintln!("skip: conformance corpus not present");
        return;
    };
    for name in VECTORS {
        let samples = read_pcm(&root.join(format!("g729-core/{name}.IN")));
        let bit_bytes = std::fs::read(root.join(format!("g729-core/{name}.BIT"))).unwrap();
        let n_frames = (samples.len() / SAMPLES_PER_FRAME).min(bit_bytes.len() / FRAME_BYTES);
        let mut enc = FrameEncoderFx::new();
        let (mut total, mut miss, mut ref_better, mut ours_better) =
            (0usize, 0usize, 0usize, 0usize);
        for f in 0..n_frames {
            let mut frame = [0i16; SAMPLES_PER_FRAME];
            frame.copy_from_slice(&samples[f * SAMPLES_PER_FRAME..(f + 1) * SAMPLES_PER_FRAME]);
            let rf = parse_frame(&bit_bytes[f * FRAME_BYTES..(f + 1) * FRAME_BYTES]).unwrap();
            let FrameKind::Active(_) = rf else {
                enc.encode_frame(&frame);
                continue;
            };
            let r = unpack_parameters(&rf).unwrap();
            let out = enc.encode_frame_locked(&frame, &r);
            for sub in 0..2 {
                total += 1;
                let (oc, os, rc, rs) = if sub == 0 {
                    (out.params.c1, out.params.s1, r.c1, r.s1)
                } else {
                    (out.params.c2, out.params.s2, r.c2, r.s2)
                };
                if oc == rc && os == rs {
                    continue;
                }
                miss += 1;
                let p = &enc.probe[sub];
                let xf: [f32; 40] = std::array::from_fn(|n| f32::from(p.x_prime[n]));
                let hf: [f32; 40] = std::array::from_fn(|n| f32::from(p.h_pre[n]) / 4096.0);
                let d = correlation_d(&xf, &hf);
                let phi = phi_matrix(&hf);
                let score = |c: u16, s: u8| -> f64 {
                    let pulses = decode_pulses(c, s).unwrap();
                    let mut num = 0.0f64;
                    let mut en = 0.0f64;
                    for i in 0..4 {
                        let mi = usize::from(pulses.pulses[i].position);
                        let si = f64::from(pulses.pulses[i].sign);
                        num += si * f64::from(d[mi]);
                        for j in 0..4 {
                            let mj = usize::from(pulses.pulses[j].position);
                            let sj = f64::from(pulses.pulses[j].sign);
                            en += si * sj * f64::from(phi[mi][mj]);
                        }
                    }
                    if en <= 0.0 {
                        0.0
                    } else {
                        num * num / en
                    }
                };
                let so = score(oc, os);
                let sr = score(rc, rs);
                if sr > so * (1.0 + 1e-6) {
                    ref_better += 1;
                } else {
                    ours_better += 1;
                }
            }
        }
        println!(
            "{name}: subframes={total} C/S misses={miss} ref-scores-higher={ref_better} ours-scores-higher={ours_better}"
        );
    }
}
