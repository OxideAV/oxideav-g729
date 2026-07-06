//! **Annex A encoder** — the reduced-complexity §A.3 analysis chain,
//! bit-stream interoperable with the main body (Table 8 / Table 1 bit
//! allocation unchanged).
//!
//! Transcribed from the in-repo Recommendation Annex A prose +
//! equation rasters (eqs (A.1)–(A.10)). Everything Annex A marks
//! "same as clause X" reuses the main-body module directly; the five
//! prose-pinned simplifications implemented here are:
//!
//! * **§A.3.3 perceptual weighting** — `W(z) = Â(z)/Â(z/γ)` with a
//!   *fixed* `γ = 0.75` on the **quantized** LP coefficients (no §3.3
//!   LAR/γ adaptation, no unquantised track). The combination
//!   `W(z)/Â(z) = 1/Â(z/γ)` collapses the impulse-response and target
//!   computations to a single all-pole filter.
//! * **§A.3.3/§A.3.4 low-pass weighted speech** — the open-loop signal
//!   is the LP residual (eq (A.3), quantized `â`) filtered through
//!   `1/[Â(z/γ)(1 − 0.7·z⁻¹)]` (eq (A.2); the printed sum runs
//!   `i = 1…10`, so the degree-11 product polynomial's last tap is
//!   dropped exactly as printed).
//! * **§A.3.4 fast open-loop pitch** — eq (A.4) computes the three
//!   range maxima over the **even samples only**
//!   (`R(k) = Σ_{n=0}^{39} s_w(2n)·s_w(2n−k)`), eq (A.5) normalises by
//!   the lagged even-sample energy, and in the third range
//!   `[80, 143]` only **even delays** are correlated in the first
//!   pass, then the winner's `±1` neighbours are tested. The winner
//!   among the three ranges favours lower delays (the exact
//!   augmentation rule for submultiples is prose-unpinned; the
//!   clause-3.4 `0.85` cascade of the main-body module is used —
//!   black-box latitude).
//! * **§A.3.7 fast adaptive-codebook search** — eq (A.7) maximises the
//!   correlation `Σ x_b(n)·u_k(n)` between the **backward-filtered
//!   target** `x_b(n) = Σ_{m=n}^{39} x(m)·h(m−n)` and the *unfiltered*
//!   past excitation (the eq (A.6) energy denominator is dropped).
//!   Fractions `−⅓, 0, +⅓` around the integer optimum are tested by
//!   interpolating the past excitation with the clause-3.7 `b30`
//!   filter (eq (A.8)), for `T_1 < 85` and always for `T_2`.
//! * **§A.3.10 memory update** — the `1/Â(z/γ)` target-filter states
//!   are updated without filtering, via eq (A.10)
//!   `e_w(n) = x(n) − ĝ_p·y(n) − ĝ_c·z(n)` evaluated at `n = 30 … 39`.
//!
//! **Not implemented (documented gap):** the §A.3.8.1 depth-first
//! ACELP tree search — Annex A pins only its *existence* ("a smaller
//! number of pulse position combinations is tested and it has fixed
//! complexity") and defers the pulse-combination schedule to the
//! barred reference C. The main-body §3.8.1 focused search is used
//! instead (bit-stream legal; the emitted `C/S` codewords are a valid
//! encoder choice but will not match the G.729A reference's).

use crate::closed_loop_pitch::{
    subframe1_window, subframe2_window, EXC_BUFFER, PIT_MAX, T1_FRACTIONAL_LIMIT,
};
use crate::encode_target::lp_residual;
use crate::fixed_codebook::{encode_positions, encode_signs};
use crate::fixed_codebook_search::{
    adaptive_gain, correlation_d, filter_through_h, phi_matrix, prefilter_impulse_response,
    search_fixed_codebook, update_target, MAX_LOOP4_PER_FRAME,
};
use crate::gain_index_map::{map_ga, map_gb};
use crate::gain_predict::GainPredictor;
use crate::gain_quantize::{quantize_gains, GainTerms};
use crate::levinson::levinson;
use crate::lp_analysis::analyze_window;
use crate::lp_to_lsp::lp_to_lsp;
use crate::lsp_interpolate::{omega_to_q, LspInterpolator};
use crate::lsp_quantize::LspQuantizer;
use crate::lsp_to_lp::lsp_to_lp;
use crate::open_loop_pitch::PIT_BUFFER;
use crate::parameters::Parameters;
use crate::pitch_decode::{derive_t_min, encode_p1, encode_p2, PitchDelay};
use crate::pitch_sharpen::{clamp_beta, sharpen};
use crate::preprocess::Preprocessor;
use crate::serial::write_frame;
use crate::tables::{L_WINDOW, M, PITCH_INTERP_FILTER_SYNTHESIS_Q15};
use crate::taming::Taming;

/// Samples per frame / subframe (clause 2 framing).
pub const L_FRAME: usize = 80;
/// Samples per subframe.
pub const L_SUBFR: usize = 40;

/// §A.3.3 fixed weighting factor `γ = 0.75`.
pub const GAMMA_A: f32 = 0.75;

/// §A.3.3 low-pass pole of the open-loop weighting (`1 − 0.7·z⁻¹`).
pub const LOWPASS_POLE: f32 = 0.7;

/// The three §A.3.4 delay ranges (inclusive), longest first — the
/// same clause-3.4 section order as the main-body module, so the
/// favour cascade walks toward the **lower** delay ranges.
const RANGES: [(usize, usize); 3] = [(80, 143), (40, 79), (20, 39)];

/// Favour-lower-delays cascade factor (clause-3.4 value; the Annex A
/// submultiple-augmentation rule is prose-unpinned).
const FAVOUR_FACTOR: f32 = 0.85;

/// Weighted coefficients `γ^i·â_i` of `Â(z/γ)`.
fn weighted_a(a_hat: &[f32; M]) -> [f32; M] {
    let mut w = [0.0f32; M];
    let mut g = 1.0f32;
    for (wi, &ai) in w.iter_mut().zip(a_hat.iter()) {
        g *= GAMMA_A;
        *wi = g * ai;
    }
    w
}

/// §A.3.5: the impulse response `h(n)` of `W(z)/Â(z) = 1/Â(z/γ)`.
#[must_use]
pub fn impulse_response_a(a_hat: &[f32; M]) -> [f32; L_SUBFR] {
    let wa = weighted_a(a_hat);
    let mut h = [0.0f32; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = if n == 0 { 1.0 } else { 0.0 };
        for (i, &w) in wa.iter().enumerate() {
            if n > i {
                acc -= w * h[n - 1 - i];
            }
        }
        h[n] = acc;
    }
    h
}

/// The §A.3.6 target filter `1/Â(z/γ)`: ten-sample output memory,
/// updated per eq (A.10) without any filtering.
#[derive(Debug, Clone, Default)]
struct TargetFilterA {
    /// `x(n−1) … x(n−10)` (most recent first).
    mem: [f32; M],
}

impl TargetFilterA {
    /// `x(n) = r(n) − Σ γ^i·â_i·x(n−i)` over one subframe.
    fn target(&self, residual: &[f32; L_SUBFR], a_hat: &[f32; M]) -> [f32; L_SUBFR] {
        let wa = weighted_a(a_hat);
        let mut x = [0.0f32; L_SUBFR];
        for n in 0..L_SUBFR {
            let mut acc = residual[n];
            for (i, &w) in wa.iter().enumerate() {
                let past = if n > i { x[n - 1 - i] } else { self.mem[i - n] };
                acc -= w * past;
            }
            x[n] = acc;
        }
        x
    }

    /// eq (A.10): the filter states for the next subframe are the
    /// weighted error `e_w(n) = x(n) − ĝ_p·y(n) − ĝ_c·z(n)` at
    /// `n = 30 … 39` (no filtering operations).
    fn update(
        &mut self,
        x: &[f32; L_SUBFR],
        y: &[f32; L_SUBFR],
        z: &[f32; L_SUBFR],
        gp: f32,
        gc: f32,
    ) {
        for i in 0..M {
            let n = L_SUBFR - 1 - i;
            self.mem[i] = x[n] - gp * y[n] - gc * z[n];
        }
    }
}

/// §A.3.4 fast open-loop pitch over one frame of low-pass weighted
/// speech (`buf` layout: 143 history samples then the 80-sample
/// frame, as [`crate::open_loop_pitch::open_loop_pitch`]).
#[must_use]
pub fn open_loop_pitch_fast(buf: &[f32; PIT_BUFFER]) -> usize {
    // eq (A.4): decimated correlation over the even frame samples.
    let corr = |k: usize| -> f32 {
        let mut acc = 0.0f32;
        for n in 0..L_FRAME / 2 {
            acc += buf[PIT_MAX + 2 * n] * buf[PIT_MAX + 2 * n - k];
        }
        acc
    };
    // eq (A.5) denominator: lagged even-sample energy.
    let energy = |k: usize| -> f32 {
        let mut acc = 0.0f32;
        for n in 0..L_FRAME / 2 {
            let v = buf[PIT_MAX + 2 * n - k];
            acc += v * v;
        }
        acc
    };

    let mut cand = [(0usize, 0.0f32); 3];
    for (s, &(lo, hi)) in RANGES.iter().enumerate() {
        let (mut best_k, mut best_r) = (lo, f32::NEG_INFINITY);
        if s == 0 {
            // Longest range `[80, 143]`: even delays first …
            let mut k = lo;
            while k <= hi {
                let r = corr(k);
                if r > best_r {
                    best_r = r;
                    best_k = k;
                }
                k += 2;
            }
            // … then the winner's odd neighbours.
            for k in [best_k.saturating_sub(1), best_k + 1] {
                if (lo..=hi).contains(&k) {
                    let r = corr(k);
                    if r > best_r {
                        best_r = r;
                        best_k = k;
                    }
                }
            }
        } else {
            for k in lo..=hi {
                let r = corr(k);
                if r > best_r {
                    best_r = r;
                    best_k = k;
                }
            }
        }
        cand[s] = (best_k, best_r);
    }

    // eq (A.5) normalisation of the three retained maxima.
    let mut norm = [0.0f32; 3];
    for (s, &(k, r)) in cand.iter().enumerate() {
        let e = energy(k);
        norm[s] = if e > 0.0 { r / e.sqrt() } else { 0.0 };
    }

    // Favour lower delays (clause-3.4 cascade; see module docs).
    let mut t_op = cand[0].0;
    let mut r_best = norm[0];
    if norm[1] >= FAVOUR_FACTOR * r_best {
        r_best = norm[1];
        t_op = cand[1].0;
    }
    if norm[2] >= FAVOUR_FACTOR * r_best {
        t_op = cand[2].0;
    }
    t_op
}

/// eq (A.8) / eq (40): the past excitation interpolated at delay
/// `(int_t, frac ∈ {−1, 0, 1})` with the clause-3.7 `b30` filter.
/// `exc` uses the [`EXC_BUFFER`] layout (`u(n)` at `PIT_MAX + n`).
fn interpolated_excitation(exc: &[f32; EXC_BUFFER], delay: PitchDelay) -> [f32; L_SUBFR] {
    const Q15: f32 = 32_768.0;
    let (k, t) = match delay.frac {
        0 => (delay.int_t, 0),
        1 => (delay.int_t, 1),
        _ => (delay.int_t - 1, 2),
    };
    let u = |idx: i32| -> f32 {
        let pos = PIT_MAX as i32 + idx;
        if pos >= 0 && (pos as usize) < EXC_BUFFER {
            exc[pos as usize]
        } else {
            0.0
        }
    };
    std::array::from_fn(|n| {
        let n = n as i32;
        let mut acc = 0.0f32;
        for i in 0..10i32 {
            acc +=
                u(n - k - i) * f32::from(PITCH_INTERP_FILTER_SYNTHESIS_Q15[(t + 3 * i) as usize]);
            acc += u(n - k + 1 + i)
                * f32::from(PITCH_INTERP_FILTER_SYNTHESIS_Q15[(3 - t + 3 * i) as usize]);
        }
        acc / Q15
    })
}

/// §A.3.7 fast adaptive-codebook search: maximise the eq (A.7)
/// correlation of the backward-filtered target with the past
/// excitation, integers first, then the `±⅓` fractions around the
/// optimum (when the fractional resolution applies).
fn fast_closed_loop_search(
    x: &[f32; L_SUBFR],
    h: &[f32; L_SUBFR],
    exc: &[f32; EXC_BUFFER],
    t_min: i32,
    t_max: i32,
    frac_limit: i32,
) -> PitchDelay {
    // x_b(n) = Σ_{m=n}^{39} x(m)·h(m−n) — the backward-filtered target.
    let mut xb = [0.0f32; L_SUBFR];
    for (n, out) in xb.iter_mut().enumerate() {
        let mut acc = 0.0f32;
        for m in n..L_SUBFR {
            acc += x[m] * h[m - n];
        }
        *out = acc;
    }

    // Integer pass: R_N(k) = Σ x_b(n)·u(n−k).
    let mut best_k = t_min;
    let mut best_r = f32::NEG_INFINITY;
    for k in t_min..=t_max {
        let mut acc = 0.0f32;
        for (n, &xbn) in xb.iter().enumerate() {
            acc += xbn * exc[(PIT_MAX as i32 + n as i32 - k) as usize];
        }
        if acc > best_r {
            best_r = acc;
            best_k = k;
        }
    }

    // Fractional pass (eq (A.8) interpolation), when applicable.
    let mut best = PitchDelay {
        int_t: best_k,
        frac: 0,
    };
    if best_k < frac_limit {
        for frac in [-1i32, 1] {
            let cand = PitchDelay {
                int_t: best_k,
                frac,
            };
            // A −⅓ fraction reaches into delay `best_k − 1 + ⅔`; keep
            // the candidate only when it stays inside the window.
            if frac < 0 && best_k - 1 < t_min.min(20) {
                continue;
            }
            let ukt = interpolated_excitation(exc, cand);
            let mut acc = 0.0f32;
            for (n, &xbn) in xb.iter().enumerate() {
                acc += xbn * ukt[n];
            }
            if acc > best_r {
                best_r = acc;
                best = cand;
            }
        }
    }
    best
}

/// One encoded Annex A frame: the Table-8 codewords plus the local
/// decoded gains and delays for inspection.
#[derive(Debug, Clone)]
pub struct AnnexAEncodedFrame {
    /// The 15 transmitted codewords.
    pub params: Parameters,
    /// Per-subframe quantised `(ĝ_p, ĝ_c)`.
    pub gains: [(f32, f32); 2],
    /// Per-subframe selected pitch delays.
    pub delays: [PitchDelay; 2],
}

/// The stateful Annex A frame encoder (§A.3).
#[derive(Debug, Clone)]
pub struct AnnexAEncoder {
    preproc: Preprocessor,
    speech: [f32; L_WINDOW],
    lsp_quant: LspQuantizer,
    /// §A.3.2.5: only the *quantized* LSP track is interpolated.
    interp_q: LspInterpolator,
    prev_q_unq: [f32; M],
    /// Low-pass weighted-speech history for the §A.3.4 correlation.
    sw_hist: [f32; PIT_MAX],
    /// eq (A.2) filter memory (`s_w(n−1) … s_w(n−10)`, most recent
    /// first) — the low-pass weighting filter is recursive across
    /// subframes.
    sw_mem: [f32; M],
    /// Past-speech tail `s(−1) … s(−10)` for the eq (A.3) residual.
    s_past: [f32; M],
    exc_hist: [f32; PIT_MAX],
    target_filter: TargetFilterA,
    gain_pred: GainPredictor,
    prev_gp_hat: f32,
    taming: Taming,
}

impl Default for AnnexAEncoder {
    fn default() -> Self {
        Self::new()
    }
}

impl AnnexAEncoder {
    /// Build an encoder in the clause-4.3 start-up state (§A.4.3
    /// defers to clause 4.3).
    #[must_use]
    pub fn new() -> Self {
        Self {
            preproc: Preprocessor::new(),
            speech: [0.0; L_WINDOW],
            lsp_quant: LspQuantizer::new(),
            interp_q: LspInterpolator::new(),
            prev_q_unq: std::array::from_fn(|i| {
                ((i + 1) as f32 * std::f32::consts::PI / 11.0).cos()
            }),
            sw_hist: [0.0; PIT_MAX],
            sw_mem: [0.0; M],
            s_past: [0.0; M],
            exc_hist: [0.0; PIT_MAX],
            target_filter: TargetFilterA::default(),
            gain_pred: GainPredictor::new(),
            prev_gp_hat: crate::decode_chain::BETA_INIT,
            taming: Taming::new(),
        }
    }

    /// eq (A.2): one subframe of low-pass weighted speech from the
    /// eq (A.3) residual, through `1/[Â(z/γ)(1 − 0.7·z⁻¹)]` with the
    /// printed ten-tap denominator sum.
    fn lowpass_weighted(&mut self, residual: &[f32; L_SUBFR], a_hat: &[f32; M]) -> [f32; L_SUBFR] {
        let wa = weighted_a(a_hat);
        // A′(z) = Â(z/γ)(1 − 0.7z⁻¹), truncated to taps 1 … 10 as
        // printed in eq (A.2).
        let mut ap = [0.0f32; M];
        for i in 0..M {
            let prev = if i == 0 { 1.0 } else { wa[i - 1] };
            ap[i] = wa[i] - LOWPASS_POLE * prev;
        }
        let mut sw = [0.0f32; L_SUBFR];
        for n in 0..L_SUBFR {
            let mut acc = residual[n];
            for (i, &c) in ap.iter().enumerate() {
                let past = if n > i {
                    sw[n - 1 - i]
                } else {
                    self.sw_mem[i - n]
                };
                acc -= c * past;
            }
            sw[n] = acc;
        }
        for i in 0..M {
            self.sw_mem[i] = sw[L_SUBFR - 1 - i];
        }
        sw
    }

    /// Encode one 80-sample frame of 16-bit PCM (as `f32` sample
    /// values) into the Table-8 codewords.
    #[allow(clippy::too_many_lines)]
    pub fn encode_frame(&mut self, pcm: &[f32; L_FRAME]) -> AnnexAEncodedFrame {
        // §A.3.1 pre-processing + Figure-5 buffer slide.
        let s_new = self.preproc.process_frame(pcm);
        self.speech.copy_within(L_FRAME.., 0);
        self.speech[L_WINDOW - L_FRAME..].copy_from_slice(&s_new);

        // §A.3.2.1/2/3 — unchanged LP analysis (Q12 grid as the main
        // chain; Annex A's 50-point / 2-bisection root search is a
        // complexity reduction, not a different result grid — the
        // main-body root search is kept).
        let r = analyze_window(&self.speech);
        let lev = levinson(&r);
        let a_q12: [f32; M] = std::array::from_fn(|i| (lev.a[i] * 4096.0).round() / 4096.0);
        let q_unq = match lp_to_lsp(&a_q12) {
            Some(q) => {
                self.prev_q_unq = q;
                q
            }
            None => self.prev_q_unq,
        };

        // §A.3.2.4 — unchanged LSP quantisation.
        let omega_unq: [f32; M] = std::array::from_fn(|i| {
            crate::lsf_conversion::lsp_to_lsf_q13(q_unq[i]) as f32 / 8192.0
        });
        let lsp_out = self.lsp_quant.quantize_lsf(&omega_unq);
        let q_quant = omega_to_q(&lsp_out.omega_hat);

        // §A.3.2.5 — only the quantized track is interpolated.
        let q_sub = self.interp_q.interpolate(&q_quant);
        let a_hat: [[f32; M]; 2] = [lsp_to_lp(&q_sub[0]), lsp_to_lp(&q_sub[1])];

        let s_frame: [f32; L_FRAME] = std::array::from_fn(|n| self.speech[120 + n]);

        // §A.3.3/§A.3.4: per-subframe residual → low-pass weighted
        // speech → frame-level fast open-loop pitch.
        let mut sw_buf = [0.0f32; PIT_BUFFER];
        sw_buf[..PIT_MAX].copy_from_slice(&self.sw_hist);
        let mut residuals = [[0.0f32; L_SUBFR]; 2];
        let mut s_past = self.s_past;
        for sub in 0..2 {
            let s_sub: [f32; L_SUBFR] = std::array::from_fn(|n| s_frame[sub * L_SUBFR + n]);
            residuals[sub] = lp_residual(&s_sub, &s_past, &a_hat[sub]);
            let sw = self.lowpass_weighted(&residuals[sub], &a_hat[sub]);
            sw_buf[PIT_MAX + sub * L_SUBFR..PIT_MAX + (sub + 1) * L_SUBFR].copy_from_slice(&sw);
            s_past = std::array::from_fn(|i| s_sub[L_SUBFR - 1 - i]);
        }
        let t_op = open_loop_pitch_fast(&sw_buf);
        self.sw_hist.copy_within(L_FRAME.., 0);
        self.sw_hist[PIT_MAX - L_FRAME..].copy_from_slice(&sw_buf[PIT_MAX..]);

        // Per-subframe closed-loop analysis.
        let mut loop4_budget = MAX_LOOP4_PER_FRAME;
        let mut params = Parameters {
            l0: lsp_out.indices.l0 as u8,
            l1: lsp_out.indices.l1 as u8,
            l2: lsp_out.indices.l2 as u8,
            l3: lsp_out.indices.l3 as u8,
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
        };
        let mut gains = [(0.0f32, 0.0f32); 2];
        let mut delays = [PitchDelay { int_t: 0, frac: 0 }; 2];
        let mut t1_delay: Option<PitchDelay> = None;
        let mut exc_frame = [0.0f32; L_FRAME];

        for sub in 0..2 {
            let s_sub: [f32; L_SUBFR] = std::array::from_fn(|n| s_frame[sub * L_SUBFR + n]);
            let residual = residuals[sub];

            // §A.3.5 impulse response + §A.3.6 target (single filter).
            let h = impulse_response_a(&a_hat[sub]);
            let x = self.target_filter.target(&residual, &a_hat[sub]);

            // Excitation buffer (clause 3.7: residual copied into the
            // current-subframe region so delays < 40 are searchable).
            let mut exc = [0.0f32; EXC_BUFFER];
            exc[..PIT_MAX].copy_from_slice(&self.exc_hist);
            if sub == 1 {
                exc[..PIT_MAX - L_SUBFR].copy_from_slice(&self.exc_hist[L_SUBFR..]);
                exc[PIT_MAX - L_SUBFR..PIT_MAX].copy_from_slice(&exc_frame[..L_SUBFR]);
            }
            exc[PIT_MAX..PIT_MAX + L_SUBFR].copy_from_slice(&residual);

            // §A.3.7 fast adaptive-codebook search.
            let (t_min, t_max) = if sub == 0 {
                subframe1_window(t_op as i32)
            } else {
                subframe2_window(t1_delay.expect("subframe order").int_t)
            };
            let frac_limit = if sub == 0 {
                T1_FRACTIONAL_LIMIT
            } else {
                i32::MAX
            };
            let delay = fast_closed_loop_search(&x, &h, &exc, t_min, t_max, frac_limit);
            delays[sub] = delay;

            // §A.3.7.1–3 (same as main): adaptive vector, filtered
            // vector, adaptive gain, fixed-codebook target.
            let v = interpolated_excitation(&exc, delay);
            let y = filter_through_h(&v, &h);
            let gp_ol = adaptive_gain(&x, &y);
            let x_prime = update_target(&x, &y, gp_ol);

            // §A.3.8: main-body focused search (see module docs — the
            // Annex A depth-first schedule is prose-unpinned).
            let beta = clamp_beta(self.prev_gp_hat);
            let mut h_pre = h;
            prefilter_impulse_response(&mut h_pre, delay.int_t.max(0) as usize, beta);
            let d = correlation_d(&x_prime, &h_pre);
            let phi = phi_matrix(&h_pre);
            let choice = search_fixed_codebook(&d, &phi, &mut loop4_budget);
            let mut c_i8 = [0i8; L_SUBFR];
            for (pos, sign) in choice.positions.iter().zip(choice.signs.iter()) {
                c_i8[*pos as usize] += *sign;
            }
            let c = sharpen(&c_i8, delay.int_t, self.prev_gp_hat);
            let z = filter_through_h(&c, &h);

            // §A.3.9 — unchanged gain quantisation (incl. taming).
            let pred = self.gain_pred.predict_only(&c);
            let terms = GainTerms::compute(&x, &y, &z);
            let tame = self.taming.test(delay.int_t, delay.frac);
            let gq = quantize_gains(&terms, pred.g_c_prime, tame);
            self.gain_pred.push_quantised_error(gq.gamma_hat);
            self.prev_gp_hat = gq.gp_hat;
            gains[sub] = (gq.gp_hat, gq.gc_hat);
            self.taming.update(delay.int_t, delay.frac, gq.gp_hat);

            // §A.3.10 memory updates: excitation + eq (A.10) states.
            let u: [f32; L_SUBFR] = std::array::from_fn(|n| gq.gp_hat * v[n] + gq.gc_hat * c[n]);
            exc_frame[sub * L_SUBFR..(sub + 1) * L_SUBFR].copy_from_slice(&u);
            self.target_filter.update(&x, &y, &z, gq.gp_hat, gq.gc_hat);
            self.s_past = std::array::from_fn(|i| s_sub[L_SUBFR - 1 - i]);

            // Transmitted codewords.
            let c_code = encode_positions(&choice.positions).expect("track-valid positions");
            let s_code = encode_signs(&choice.signs).expect("±1 signs");
            let ga_tx = map_ga(gq.ga).expect("in-range") as u8;
            let gb_tx = map_gb(gq.gb).expect("in-range") as u8;
            if sub == 0 {
                let p1 = encode_p1(delay).expect("window-bounded T1");
                params.p1 = p1;
                params.p0 = (((p1 >> 2).count_ones() ^ 1) & 1) as u8;
                params.c1 = c_code;
                params.s1 = s_code;
                params.ga1 = ga_tx;
                params.gb1 = gb_tx;
                t1_delay = Some(delay);
            } else {
                let t_min2 = derive_t_min(t1_delay.expect("subframe order").int_t);
                params.p2 = encode_p2(delay, t_min2).expect("window-bounded T2");
                params.c2 = c_code;
                params.s2 = s_code;
                params.ga2 = ga_tx;
                params.gb2 = gb_tx;
            }
        }

        self.exc_hist.copy_within(L_FRAME.., 0);
        self.exc_hist[PIT_MAX - L_FRAME..].copy_from_slice(&exc_frame);

        AnnexAEncodedFrame {
            params,
            gains,
            delays,
        }
    }

    /// Encode one frame straight to the 164-byte ITU serial format.
    pub fn encode_frame_to_serial(&mut self, pcm: &[f32; L_FRAME]) -> [u8; 164] {
        let out = self.encode_frame(pcm);
        write_frame(&crate::parameters::pack_bit_array(&out.params))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `1/Â(z/γ)` impulse response of a flat filter is the unit
    /// impulse.
    #[test]
    fn impulse_response_flat_is_delta() {
        let a = [0.0f32; M];
        let h = impulse_response_a(&a);
        assert!((h[0] - 1.0).abs() < 1e-6);
        assert!(h[1..].iter().all(|&v| v.abs() < 1e-6));
    }

    /// The fast open-loop estimator finds the fundamental of a
    /// periodic signal and rejects its double (favour-lower rule).
    #[test]
    fn fast_open_loop_finds_period() {
        for period in [30usize, 55] {
            let buf: [f32; PIT_BUFFER] = std::array::from_fn(|n| {
                let phase = (n % period) as f32 / period as f32;
                (std::f32::consts::TAU * phase).sin()
                    + 0.3 * (2.0 * std::f32::consts::TAU * phase).sin()
            });
            let t = open_loop_pitch_fast(&buf);
            assert!(
                t.abs_diff(period) <= 1,
                "period {period}: got {t} (decimated estimate)"
            );
        }
    }

    /// A long odd period is located to within the decimated
    /// estimator's resolution: the even-first pass lands on a
    /// neighbouring even delay and the ±1 refinement pulls toward the
    /// true period. (The §A.3.4 estimate is only the anchor of the
    /// ±-samples closed-loop window, so ±2 is the meaningful band —
    /// eq (A.4) sees half the samples.)
    #[test]
    fn fast_open_loop_third_range_odd_delay() {
        let period = 101usize;
        let buf: [f32; PIT_BUFFER] = std::array::from_fn(|n| {
            let phase = (n % period) as f32 / period as f32;
            (std::f32::consts::TAU * phase).sin()
                + 0.4 * (2.0 * std::f32::consts::TAU * phase).sin()
        });
        let t = open_loop_pitch_fast(&buf);
        assert!(t.abs_diff(period) <= 2, "got {t}");
    }

    /// eq (A.10): the target-filter memory after an update equals the
    /// weighted error tail, with no filtering.
    #[test]
    fn target_memory_is_weighted_error_tail() {
        let mut tf = TargetFilterA::default();
        let x: [f32; L_SUBFR] = std::array::from_fn(|n| n as f32);
        let y: [f32; L_SUBFR] = std::array::from_fn(|n| 0.5 * n as f32);
        let z: [f32; L_SUBFR] = std::array::from_fn(|n| 0.25 * n as f32);
        tf.update(&x, &y, &z, 0.8, 2.0);
        for i in 0..M {
            let n = L_SUBFR - 1 - i;
            let expect = x[n] - 0.8 * y[n] - 2.0 * z[n];
            assert!((tf.mem[i] - expect).abs() < 1e-6, "i={i}");
        }
    }

    /// Every emitted codeword stays inside its Table-8 field domain
    /// over a multi-frame synthetic input, and the serial path emits
    /// parseable frames.
    #[test]
    fn codewords_stay_in_field_domains() {
        let mut enc = AnnexAEncoder::new();
        for f in 0..12 {
            let pcm: [f32; L_FRAME] = std::array::from_fn(|n| {
                let t = (f * L_FRAME + n) as f32;
                4000.0 * (0.07 * t).sin() + 1500.0 * (0.19 * t).sin()
            });
            let out = enc.encode_frame(&pcm);
            let p = &out.params;
            assert!(p.l0 < 2 && p.l1 < 128 && p.l2 < 32 && p.l3 < 32);
            assert!(p.p1 < 244, "p1 = {}", p.p1);
            assert!(p.p2 < 32 && p.p0 < 2);
            assert!(p.c1 < 8192 && p.c2 < 8192 && p.s1 < 16 && p.s2 < 16);
            assert!(p.ga1 < 8 && p.gb1 < 16 && p.ga2 < 8 && p.gb2 < 16);
            assert!(p.pitch_parity_ok(), "frame {f}: parity");
        }
    }

    /// The Annex A encoder's stream decodes cleanly through the
    /// main-body decoder (bit-stream interoperability, §A.1) and the
    /// Annex A decoder.
    #[test]
    fn own_stream_decodes_cleanly() {
        let mut enc = AnnexAEncoder::new();
        let mut dec = crate::decode_chain::FrameDecoder::new();
        let mut dec_a = crate::annex_a::AnnexADecoder::new();
        for f in 0..10 {
            let pcm: [f32; L_FRAME] = std::array::from_fn(|n| {
                let t = (f * L_FRAME + n) as f32;
                3000.0 * (0.05 * t).sin()
            });
            let wire = enc.encode_frame_to_serial(&pcm);
            let synth = dec
                .decode_serial_frame_to_speech(&wire)
                .expect("main decoder accepts");
            assert!(synth.speech().iter().all(|v| v.is_finite()));
            let out = dec_a
                .decode_serial_frame_to_postfiltered(&wire)
                .expect("Annex A decoder accepts");
            assert!(out.iter().all(|v| v.is_finite()));
        }
    }
}
