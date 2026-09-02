//! Annex B **encoder side** — the §B.4.1 DTX decision and the §B.4.2
//! SID parameter evaluation / quantisation (energy + SID-LSF VQ).
//!
//! The VAD (§B.3) is **not** implemented here: its multi-boundary
//! constants are printed (Table B.1) but the parameter grids they
//! apply to, the AR coefficient sets of §B.3.7 ("different sets …
//! according to the value of `C_n`" — no values are printed) and the
//! §B.3.3 minimum-buffer segmentation are not, so the voice-activity
//! input `Vad_t` is a caller-supplied flag. Against the staged Annex B
//! corpus the flag is read off the reference bitstream's own frame
//! types (a locked-VAD drive), which measures the DTX decision and the
//! SID quantisation in isolation.
//!
//! ## §B.4.1 DTX (per inactive frame)
//!
//! - eq (B.9): `R_t(j) = Σ_{i=t−N_cur+1}^{t} r'_i(j)` over the stored
//!   per-frame autocorrelations (`N_cur = 2`), then Levinson-Durbin →
//!   `A_t(z)`, residual energy `E_t`.
//! - eq (B.10): the first inactive frame of a zone is a SID frame,
//!   `Ē = E_t`, `k_E = 1`.
//! - Otherwise `flag_chang` is raised when the LPC filters differ
//!   (§B.4.1.3, eqs (B.12)/(B.13): `Σ_j R_a(j)·R_t(j) ≥ E_t·thr1`,
//!   `thr1 = 1.20226`, `R_a` the autocorrelation of the previous
//!   SID filter's coefficients) or the energies differ (§B.4.1.4:
//!   `k_E ← min(k_E + 1, N_g = 2)`, eq (B.14) `Ē = Σ E_i` over the last
//!   `k_E` frames, eq (B.15) `E' = α_w·Ē/(k_E·N_cur·80)`, `α_w = 0.125`,
//!   quantised through the §B.4.2.1 ladder and compared to the previous
//!   SID's decoded energy against `thr2 = 2 dB`).
//! - eq (B.11): `count_fr ≥ N_min = 2` and `flag_chang` ⇒ SID; a SID
//!   resets both.
//!
//! ## §B.4.2 SID parameters
//!
//! - Energy: eq (B.15) `E'` quantised to the nearest level of the
//!   −12 … 66 dB ladder (`dequant_sid_energy_db` is the decoder's
//!   reading of the same ladder).
//! - Filter: eq (B.16) past-average autocorrelation over `N_p = 6`
//!   frames ending at `t'` → `Ā_p(z)`; eq (B.17) selects `A_t(z)` when
//!   `distance(A_t, Ā_p) ≥ thr3 = 1.12202` (the (B.12) form with `R_a`
//!   from `Ā_p`), else `Ā_p`. The selected filter's LSFs go through the
//!   §B.4.2.2 two-stage switched-predictive VQ (staged subset maps,
//!   eq (B.18) blended mode-1 predictor, eq (22) weights, exhaustive
//!   joint search — a superset of the unpinned "few candidates"
//!   delayed decision).
//!
//! ## Unpinned points (recorded, exposed as [`DtxLatitude`])
//!
//! - §B.4.2.2: "`t'` varies in `[t − 1, t − N_cur]` depending on the
//!   rest of the Euclidean division of `t` by `N_cur`" — which
//!   remainder maps to which end is not stated.
//! - §B.4.1.3: whether `A_sid(z)` for the comparison is the quantised
//!   (decoder-side) SID filter or the pre-quantisation estimate.
//! - §B.4.2.1: whether the energy index is the nearest ladder level or
//!   the level at or below `E'`.

use crate::annex_b::{dequant_sid_energy_db, SidFrame};
use crate::levinson::levinson;
use crate::lp_to_lsp::lp_to_lsp;
use crate::lsp_quantize::lsf_weights;
use crate::sid_lsf::sid_reconstruct;
use crate::tables::{M, MA_NP};

/// Autocorrelation lags stored per frame (`r'(0 … 10)`).
pub const N_LAGS: usize = M + 1;
/// Frames summed for the current inactive-frame LPC (`N_cur`).
pub const N_CUR: usize = 2;
/// Frames summed for the past-average SID filter (`N_p`).
pub const N_P: usize = 6;
/// Maximum frames in the energy sum (`N_g`).
pub const N_G: usize = 2;
/// Minimum frame interval between SID frames (`N_min`).
pub const N_MIN: usize = 2;
/// §B.4.1.3 LPC-change threshold.
pub const THR1: f64 = 1.202_26;
/// §B.4.1.4 energy-change threshold (dB).
pub const THR2_DB: f64 = 2.0;
/// §B.4.2.2 stationarity threshold.
pub const THR3: f64 = 1.122_02;
/// eq (B.15) windowing / bandwidth-expansion scaling `α_w`.
pub const ALPHA_W: f64 = 0.125;
/// Corpus-pinned scale of the eq (B.15) energy relative to the §3.1
/// pre-processed (halved) signal's autocorrelation: the reference's
/// SID energies sit 6 dB above `α_w·Ē/(k_E·N_cur·80)` computed on the
/// halved signal on every sequence (medians −5.9, −5.9, −6.4 dB), i.e.
/// the energy is on the un-halved input scale.
pub const ENERGY_INPUT_SCALE: f64 = 4.0;

/// Unstated points of the DTX / SID clauses, exposed for black-box
/// sweeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DtxLatitude {
    /// `t' = t − 1 − (t mod N_cur)` (`true`) or `t' = t − N_cur + (t mod
    /// N_cur)` (`false`).
    pub t_prime_high_on_even: bool,
    /// Compare against the quantised SID filter (`true`) or the
    /// pre-quantisation estimate (`false`).
    pub compare_quantised: bool,
    /// Nearest ladder level (`true`) or the level at or below `E'`.
    pub energy_nearest: bool,
    /// Scale applied to the eq (B.15) energy ([`ENERGY_INPUT_SCALE`]
    /// by default; `1.0` is the literal reading on the halved signal).
    pub energy_scale_x100: u32,
}

impl Default for DtxLatitude {
    fn default() -> Self {
        Self {
            t_prime_high_on_even: true,
            compare_quantised: true,
            energy_nearest: true,
            energy_scale_x100: (ENERGY_INPUT_SCALE * 100.0) as u32,
        }
    }
}

/// One frame's DTX decision (`Ftyp_t`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DtxDecision {
    /// `Ftyp = 1`: active speech, the normal encoder runs.
    Active,
    /// `Ftyp = 2`: a SID frame with these Table B.2 indices.
    Sid(SidFrame),
    /// `Ftyp = 0`: untransmitted.
    Untransmitted,
}

/// The absolute (un-normalised) §3.2.1 autocorrelation `r'(0 … 10)` of
/// a fixed-point front-end result: the DPF values are on the
/// `2^norm`-normalised doubled-product grid of a signal right-shifted
/// by `signal_shift`, so `r' = value · 2^(2·shift − norm) / 2`.
#[must_use]
pub fn absolute_autocorrelation(ac: &crate::fx::analysis::AutocorrFx) -> [f64; N_LAGS] {
    let v = crate::fx::analysis::dpf_to_f64(&ac.r);
    let scale = 2f64.powi(2 * i32::from(ac.signal_shift) - i32::from(ac.norm) - 1);
    std::array::from_fn(|k| v[k] * scale)
}

/// §B.4.2.1: the 5-bit logarithmic energy quantiser (index of the
/// ladder level for `e_prime_db`).
#[must_use]
pub fn quantize_sid_energy_db(e_prime_db: f64, nearest: bool) -> u8 {
    let mut best = 0u8;
    if nearest {
        let mut best_d = f64::INFINITY;
        for q in 0..32u8 {
            let d = (f64::from(dequant_sid_energy_db(q)) - e_prime_db).abs();
            if d < best_d {
                best_d = d;
                best = q;
            }
        }
    } else {
        for q in 0..32u8 {
            if f64::from(dequant_sid_energy_db(q)) <= e_prime_db {
                best = q;
            }
        }
    }
    best
}

/// eq (B.13): the autocorrelation `R_a(j)` of an LP coefficient set
/// `a(0 … 10)` (`a(0) = 1`), with the `j ≠ 0` terms doubled.
#[must_use]
pub fn coefficient_autocorrelation(a: &[f64; N_LAGS]) -> [f64; N_LAGS] {
    std::array::from_fn(|j| {
        let mut acc = 0.0;
        for k in 0..N_LAGS - j {
            acc += a[k] * a[k + j];
        }
        if j == 0 {
            acc
        } else {
            2.0 * acc
        }
    })
}

/// eq (B.12) left-hand side: `Σ_j R_a(j)·R(j)` — the Itakura-style
/// distance numerator of a filter (via its coefficient autocorrelation
/// `ra`) against an autocorrelation `r`.
#[must_use]
pub fn itakura_numerator(ra: &[f64; N_LAGS], r: &[f64; N_LAGS]) -> f64 {
    ra.iter().zip(r.iter()).map(|(a, b)| a * b).sum()
}

/// One evaluated SID-LSF candidate: weighted error, indices,
/// reconstruction and the residual to push.
type SidCandidate = (f32, (usize, usize, usize), [f32; M], [f32; M]);

/// The §B.4.2.2 SID-LSF quantiser: exhaustive joint search over the
/// two predictor modes, the 32-entry first-stage subset and the
/// 16-entry full-VQ second stage under the eq (22) weighted MSE, with
/// its own MA history (advanced only by SID frames, mirroring the
/// decoder's [`crate::sid_lsf::SidLspDecoder`]).
#[derive(Debug, Clone)]
pub struct SidLsfQuantizer {
    history: [[f32; M]; MA_NP],
}

impl Default for SidLsfQuantizer {
    fn default() -> Self {
        Self::new()
    }
}

impl SidLsfQuantizer {
    /// Fresh quantiser in the §3.2.4 start-up state.
    #[must_use]
    pub fn new() -> Self {
        let row: [f32; M] =
            std::array::from_fn(|i| ((i + 1) as f32) * core::f32::consts::PI / 11.0);
        Self {
            history: [row; MA_NP],
        }
    }

    /// Quantises the LSF vector `omega` (radians, ascending), returning
    /// `(lp0, l1, l2)` and the reconstructed LSF, and advancing the
    /// history with the winner's residual.
    pub fn quantize(&mut self, omega: &[f32; M]) -> ((u8, u8, u8), [f32; M]) {
        let w = lsf_weights(omega);
        let mut best: Option<SidCandidate> = None;
        for lp0 in 0..2 {
            for l1 in 0..32 {
                for l2 in 0..16 {
                    let (rec, residual) = sid_reconstruct(&self.history, lp0, l1, l2);
                    let mut e = 0.0f32;
                    for i in 0..M {
                        let d = omega[i] - rec[i];
                        e += w[i] * d * d;
                    }
                    let better = match &best {
                        None => true,
                        Some((be, ..)) => e < *be,
                    };
                    if better {
                        best = Some((e, (lp0, l1, l2), rec, residual));
                    }
                }
            }
        }
        let (_, (lp0, l1, l2), rec, residual) = best.expect("non-empty search");
        for k in (1..MA_NP).rev() {
            self.history[k] = self.history[k - 1];
        }
        self.history[0] = residual;
        ((lp0 as u8, l1 as u8, l2 as u8), rec)
    }
}

/// Stateful §B.4.1 / §B.4.2 encoder-side DTX + SID module.
#[derive(Debug, Clone)]
pub struct DtxEncoder {
    /// Frame counter `t` (the next frame's number).
    t: usize,
    /// Stored per-frame autocorrelations, most recent last;
    /// `r_hist[0]` is frame `base`.
    r_hist: Vec<[f64; N_LAGS]>,
    base: usize,
    /// Stored per-frame residual energies `E_i` (inactive frames),
    /// most recent last.
    e_hist: Vec<f64>,
    prev_vad: bool,
    k_e: usize,
    count_fr: usize,
    flag_chang: bool,
    /// Previous SID filter `a_sid(0 … 10)` (for the §B.4.1.3 compare).
    sid_a: Option<[f64; N_LAGS]>,
    /// Previous SID's decoded energy (dB).
    sid_eq_db: f64,
    lsf_quant: SidLsfQuantizer,
    lat: DtxLatitude,
    /// The last inactive frame's eq (B.15) `E'` in dB (instrument).
    #[doc(hidden)]
    pub last_e_prime_db: f64,
}

impl Default for DtxEncoder {
    fn default() -> Self {
        Self::new()
    }
}

impl DtxEncoder {
    /// Fresh module (no SID sent yet).
    #[must_use]
    pub fn new() -> Self {
        Self::with_latitude(DtxLatitude::default())
    }

    /// Fresh module under an explicit latitude.
    #[must_use]
    pub fn with_latitude(lat: DtxLatitude) -> Self {
        Self {
            t: 0,
            r_hist: Vec::new(),
            base: 0,
            e_hist: Vec::new(),
            prev_vad: true,
            k_e: 0,
            count_fr: 0,
            flag_chang: false,
            sid_a: None,
            sid_eq_db: 0.0,
            lsf_quant: SidLsfQuantizer::new(),
            lat,
            last_e_prime_db: 0.0,
        }
    }

    /// eq (B.15) with the corpus-pinned input scale.
    fn e_prime(&self, e_bar: f64) -> f64 {
        f64::from(self.lat.energy_scale_x100) / 100.0 * ALPHA_W * e_bar
            / ((self.k_e.max(1) * N_CUR * 80) as f64)
    }

    /// Sum of the stored autocorrelations over frames `lo ..= hi`
    /// (absolute frame numbers; frames before 0 or already dropped
    /// contribute nothing).
    fn sum_abs(&self, lo: isize, hi: isize) -> [f64; N_LAGS] {
        let mut acc = [0.0f64; N_LAGS];
        for f in lo.max(0)..=hi {
            let idx = f as usize;
            if idx < self.base {
                continue;
            }
            if let Some(r) = self.r_hist.get(idx - self.base) {
                for j in 0..N_LAGS {
                    acc[j] += r[j];
                }
            }
        }
        acc
    }

    /// Processes frame `t` with the caller's voice-activity flag and
    /// the frame's §3.2.1 autocorrelation `r'(0 … 10)` (bandwidth
    /// expansion and noise correction applied).
    pub fn process_frame(&mut self, vad: bool, r_t: &[f64; N_LAGS]) -> DtxDecision {
        self.process_frame_inner(vad, r_t, None)
    }

    /// Locked drive: returns the module's own decision for frame `t`
    /// but commits the **reference's** frame to the state afterwards —
    /// a reference SID frame pushes its residual into the SID-LSF
    /// history, becomes the comparison filter / energy and resets the
    /// counters; a reference untransmitted frame advances `count_fr`.
    /// `reference` is `Some(sid)` for a reference SID frame, `None`
    /// for a reference untransmitted frame (call with `vad = false`).
    #[doc(hidden)]
    pub fn process_frame_locked(
        &mut self,
        vad: bool,
        r_t: &[f64; N_LAGS],
        reference: Option<&SidFrame>,
    ) -> DtxDecision {
        self.process_frame_inner(vad, r_t, Some(reference))
    }

    fn process_frame_inner(
        &mut self,
        vad: bool,
        r_t: &[f64; N_LAGS],
        lock: Option<Option<&SidFrame>>,
    ) -> DtxDecision {
        let t = self.t;
        self.t += 1;
        self.r_hist.push(*r_t);
        // Keep enough history for the eq (B.16) window.
        let keep = N_P + N_CUR + 2;
        if self.r_hist.len() > keep + 8 {
            let drop = self.r_hist.len() - keep;
            self.r_hist.drain(..drop);
            self.base += drop;
        }
        if self.e_hist.len() > 8 {
            let drop = self.e_hist.len() - 4;
            self.e_hist.drain(..drop);
        }

        if vad {
            self.prev_vad = true;
            self.e_hist.push(0.0);
            return DtxDecision::Active;
        }

        // eq (B.9) current LPC over N_cur frames.
        let r_cur = self.sum_abs(t as isize - N_CUR as isize + 1, t as isize);
        let lev = levinson(&r_cur);
        let e_t = lev.residual_energy;
        let mut a_t = [0.0f64; N_LAGS];
        a_t[0] = 1.0;
        for i in 0..M {
            a_t[i + 1] = f64::from(lev.a[i]);
        }
        self.e_hist.push(e_t);

        let first_inactive = self.prev_vad;
        self.prev_vad = false;
        // Under lock the counters / comparison state are restored and
        // re-committed from the reference after the own decision.
        let snapshot = (
            self.k_e,
            self.count_fr,
            self.flag_chang,
            self.sid_a,
            self.sid_eq_db,
            self.lsf_quant.clone(),
        );

        let send_sid = if first_inactive {
            // eq (B.10).
            self.k_e = 1;
            true
        } else {
            // §B.4.1.3 LPC comparison against the previous SID filter.
            if let Some(sid_a) = &self.sid_a {
                let ra = coefficient_autocorrelation(sid_a);
                if itakura_numerator(&ra, &r_cur) >= e_t * THR1 {
                    self.flag_chang = true;
                }
            }
            // §B.4.1.4 energy comparison.
            self.k_e = (self.k_e + 1).min(N_G);
            let e_bar: f64 = self.e_hist.iter().rev().take(self.k_e).sum();
            let e_prime = self.e_prime(e_bar);
            let eq_db = f64::from(dequant_sid_energy_db(quantize_sid_energy_db(
                10.0 * e_prime.max(1e-30).log10(),
                self.lat.energy_nearest,
            )));
            if (eq_db - self.sid_eq_db).abs() > THR2_DB {
                self.flag_chang = true;
            }
            self.count_fr += 1;
            self.count_fr >= N_MIN && self.flag_chang
        };

        let e_bar_now: f64 = self.e_hist.iter().rev().take(self.k_e).sum();
        self.last_e_prime_db = 10.0 * self.e_prime(e_bar_now).max(1e-30).log10();
        if !send_sid {
            if let Some(reference) = lock {
                self.commit_reference(t, snapshot, reference);
            }
            return DtxDecision::Untransmitted;
        }
        self.count_fr = 0;
        self.flag_chang = false;

        // §B.4.2.1 energy.
        let e_bar: f64 = self.e_hist.iter().rev().take(self.k_e).sum();
        let e_prime = self.e_prime(e_bar);
        let gain =
            quantize_sid_energy_db(10.0 * e_prime.max(1e-30).log10(), self.lat.energy_nearest);
        self.sid_eq_db = f64::from(dequant_sid_energy_db(gain));

        // §B.4.2.2 filter: past average over N_p frames ending at t'.
        let rem = t % N_CUR;
        let t_prime = if self.lat.t_prime_high_on_even {
            t as isize - 1 - rem as isize
        } else {
            t as isize - N_CUR as isize + rem as isize
        };
        let r_past = self.sum_abs(t_prime - N_P as isize, t_prime);
        let lev_p = levinson(&r_past);
        let mut a_p = [0.0f64; N_LAGS];
        a_p[0] = 1.0;
        for i in 0..M {
            a_p[i + 1] = f64::from(lev_p.a[i]);
        }
        let ra_p = coefficient_autocorrelation(&a_p);
        let use_current = itakura_numerator(&ra_p, &r_cur) >= e_t * THR3;
        let a_sid = if use_current { a_t } else { a_p };

        // LSF of the SID filter through the §3.2.3 root search.
        let a_f32: [f32; M] = std::array::from_fn(|i| a_sid[i + 1] as f32);
        let omega: [f32; M] = match lp_to_lsp(&a_f32) {
            Some(q) => std::array::from_fn(|i| q[i].clamp(-1.0, 1.0).acos()),
            None => std::array::from_fn(|i| ((i + 1) as f32) * core::f32::consts::PI / 11.0),
        };
        let ((lp0, l1, l2), rec) = self.lsf_quant.quantize(&omega);

        // The filter the next comparisons use.
        self.sid_a = Some(if self.lat.compare_quantised {
            let q: [f32; M] = std::array::from_fn(|i| rec[i].cos());
            let a_q = crate::lsp_to_lp::lsp_to_lp(&q);
            let mut a = [0.0f64; N_LAGS];
            a[0] = 1.0;
            for i in 0..M {
                a[i + 1] = f64::from(a_q[i]);
            }
            a
        } else {
            a_sid
        });

        let own = DtxDecision::Sid(SidFrame { lp0, l1, l2, gain });
        if let Some(reference) = lock {
            self.commit_reference(t, snapshot, reference);
        }
        own
    }

    /// Locked drive: restore the pre-decision state and commit the
    /// reference frame's effect instead of our own.
    #[allow(clippy::type_complexity)]
    fn commit_reference(
        &mut self,
        t: usize,
        snapshot: (
            usize,
            usize,
            bool,
            Option<[f64; N_LAGS]>,
            f64,
            SidLsfQuantizer,
        ),
        reference: Option<&SidFrame>,
    ) {
        let (k_e, count_fr, flag_chang, sid_a, sid_eq_db, lsf_quant) = snapshot;
        self.k_e = k_e;
        self.count_fr = count_fr;
        self.flag_chang = flag_chang;
        self.sid_a = sid_a;
        self.sid_eq_db = sid_eq_db;
        self.lsf_quant = lsf_quant;
        let _ = t;
        match reference {
            Some(sid) => {
                self.count_fr = 0;
                self.flag_chang = false;
                self.sid_eq_db = f64::from(dequant_sid_energy_db(sid.gain));
                let (rec, residual) = sid_reconstruct(
                    &self.lsf_quant.history,
                    usize::from(sid.lp0),
                    usize::from(sid.l1),
                    usize::from(sid.l2),
                );
                for k in (1..MA_NP).rev() {
                    self.lsf_quant.history[k] = self.lsf_quant.history[k - 1];
                }
                self.lsf_quant.history[0] = residual;
                let q: [f32; M] = std::array::from_fn(|i| rec[i].cos());
                let a_q = crate::lsp_to_lp::lsp_to_lp(&q);
                let mut a = [0.0f64; N_LAGS];
                a[0] = 1.0;
                for i in 0..M {
                    a[i + 1] = f64::from(a_q[i]);
                }
                self.sid_a = Some(a);
            }
            None => {
                self.count_fr += 1;
            }
        }
    }
}
