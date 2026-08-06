//! §3 **per-frame encode chain** — the stateful glue that drives every
//! encoder stage in clause-3 order, turning 80 samples of 16-bit PCM
//! into the 15 transmitted Table-8 codewords per 10 ms frame.
//!
//! The capstone of the round-382 encoder build-out. Per frame:
//!
//! 1. §3.1 pre-processing ([`crate::preprocess`]).
//! 2. §3.2.1 windowing/autocorrelation/lag window
//!    ([`crate::lp_analysis`]) on the 240-sample buffer (120 past + 80
//!    present + 40 look-ahead, Figure 5) → §3.2.2 Levinson-Durbin
//!    ([`crate::levinson`]) → §3.2.3 LP→LSP ([`crate::lp_to_lsp`]).
//! 3. §3.2.4 LSP quantisation ([`crate::lsp_quantize`]) → `L0 … L3`.
//! 4. §3.2.5 per-subframe LSP interpolation (both the quantised and the
//!    unquantised sets, [`crate::lsp_interpolate`]) → §3.2.6 LSP→LP
//!    ([`crate::lsp_to_lp`]) → per-subframe `â` / `a`.
//! 5. §3.3 perceptual weighting ([`crate::perceptual_weighting`]) →
//!    weighted speech `sw(n)`.
//! 6. §3.4 open-loop pitch ([`crate::open_loop_pitch`]) once per frame.
//! 7. Per subframe: §3.5 impulse response + §3.6 residual/target
//!    ([`crate::encode_target`]); §3.7 closed-loop pitch
//!    ([`crate::closed_loop_pitch`]) → `P1`/`P2` (+ §3.7.2 parity
//!    `P0`); eq (40) adaptive vector; eq (43) gain; eq (50) updated
//!    target; eq (49) pre-filtered `h`; §3.8.1 fixed-codebook search
//!    ([`crate::fixed_codebook_search`]) → `C`/`S`; §3.9 gain
//!    quantisation ([`crate::gain_predict`] +
//!    [`crate::gain_quantize`]) → `GA`/`GB`; §3.10 memory update
//!    (eq (75) excitation, filter-state update).
//!
//! ## Frame timing (clause 3.2.1 / Figure 5)
//!
//! The LP window spans 120 past + 80 present + 40 future samples, so
//! the encoder consumes one 80-sample input frame per call but encodes
//! with 5 ms algorithmic look-ahead: the frame being encoded occupies
//! buffer positions `120 … 199`, and the caller's newest 80 samples
//! fill `160 … 239` — i.e. the encoded frame is the last 40 samples of
//! the previous input plus the first 40 of the current input.
//!
//! ## Decoder consistency
//!
//! Every quantised quantity is produced by the *decode-side* modules
//! (LSP reconstruction inside [`crate::lsp_quantize`], gain
//! reconstruction inside [`crate::gain_quantize`], the §3.9.1
//! predictor via [`crate::gain_predict::GainPredictor`]), so the
//! encoder's local decoded state tracks a real decoder bit-for-bit
//! given the same parameter stream.

use crate::closed_loop_pitch::{
    closed_loop_search, subframe1_window, subframe2_window, EXC_BUFFER, PIT_MAX,
    T1_FRACTIONAL_LIMIT,
};
use crate::encode_target::{impulse_response, lp_residual, TargetFilter, WeightedSynthesisCoefs};
use crate::fixed_codebook::{encode_positions, encode_signs};
use crate::fixed_codebook_search::{
    adaptive_gain, correlation_d, filter_through_h, phi_matrix, prefilter_impulse_response,
    search_fixed_codebook, update_target, MAX_LOOP4_PER_FRAME,
};
use crate::fx::analysis::PreprocessorFx;
use crate::gain_index_map::{map_ga, map_gb};
use crate::gain_predict::GainPredictor;
use crate::gain_quantize::{quantize_gains, GainTerms};
use crate::lp_analysis::analyze_window;
use crate::lp_to_lsp::lp_to_lsp;
use crate::lsp_interpolate::{omega_to_q, LspInterpolator};
use crate::lsp_quantize::LspQuantizer;
use crate::lsp_to_lp::lsp_to_lp;
use crate::open_loop_pitch::{open_loop_pitch, PIT_BUFFER};
use crate::parameters::Parameters;
use crate::perceptual_weighting::PerceptualWeighting;
use crate::pitch_decode::{derive_t_min, encode_p1, encode_p2, PitchDelay};
use crate::pitch_sharpen::{clamp_beta, sharpen};
use crate::tables::{L_WINDOW, M};
use crate::taming::Taming;

use crate::levinson::levinson;

/// Samples per frame.
pub const L_FRAME: usize = 80;
/// Samples per subframe.
pub const L_SUBFR: usize = 40;

/// The stateful G.729 frame encoder.
#[derive(Debug, Clone)]
pub struct FrameEncoder {
    /// §3.1 pre-processing filter (fixed-point, round 438).
    preproc: PreprocessorFx,
    /// The 240-sample analysis buffer of pre-processed speech
    /// (Figure 5 layout: `0 … 119` past, `120 … 199` present frame,
    /// `200 … 239` look-ahead).
    speech: [f32; L_WINDOW],
    /// §3.2.4 LSP quantiser (owns the decoder-consistent MA history).
    lsp_quant: LspQuantizer,
    /// §3.2.5 interpolator for the *quantised* LSP track.
    interp_q: LspInterpolator,
    /// §3.2.5 interpolator for the *unquantised* LSP track.
    interp_u: LspInterpolator,
    /// Last frame's unquantised cosine-domain LSPs (root-search
    /// fallback per clause 3.2.3).
    prev_q_unq: [f32; M],
    /// §3.3 weighting state.
    weighting: PerceptualWeighting,
    /// Weighted-speech history for the §3.4 open-loop correlation
    /// (`sw(−143) … sw(−1)`, oldest first).
    sw_hist: [f32; PIT_MAX],
    /// Excitation history `u(−143) … u(−1)` (oldest first).
    exc_hist: [f32; PIT_MAX],
    /// §3.6 target-filter memories.
    target_filter: TargetFilter,
    /// §3.9.1 gain predictor (decoder-consistent 4-tap history).
    gain_pred: GainPredictor,
    /// Previous subframe's quantised adaptive gain (eq (47) `β`
    /// source; Table 9 init 0.8).
    prev_gp_hat: f32,
    /// Past-speech tail `s(−1) … s(−10)` for the eq (36) residual.
    s_past: [f32; M],
    /// Taming-procedure per-zone worst-case excitation-error state
    /// (`docs/audio/g729/taming-procedure.md`).
    taming: Taming,
}

/// One encoded frame: the 15 Table-8 codewords (as the shared
/// [`Parameters`] struct) plus the encoder's local decoded gains for
/// inspection.
#[derive(Debug, Clone)]
pub struct EncodedFrame {
    /// The transmitted parameter set.
    pub params: Parameters,
    /// Per-subframe `(ĝ_p, ĝ_c)` the encoder's local decoder used.
    pub gains: [(f32, f32); 2],
    /// Per-subframe selected pitch delays.
    pub delays: [PitchDelay; 2],
}

impl Default for FrameEncoder {
    fn default() -> Self {
        Self::new()
    }
}

impl FrameEncoder {
    /// Fresh encoder with clause-4.3 start-up state.
    #[must_use]
    pub fn new() -> Self {
        // Clause 3.2.3 / 4.3 default LSP set: q_i = cos(iπ/11).
        let prev_q: [f32; M] =
            std::array::from_fn(|i| ((i + 1) as f32 * std::f32::consts::PI / 11.0).cos());
        Self {
            preproc: PreprocessorFx::new(),
            speech: [0.0; L_WINDOW],
            lsp_quant: LspQuantizer::new(),
            interp_q: LspInterpolator::new(),
            interp_u: LspInterpolator::new(),
            prev_q_unq: prev_q,
            weighting: PerceptualWeighting::new(),
            sw_hist: [0.0; PIT_MAX],
            exc_hist: [0.0; PIT_MAX],
            target_filter: TargetFilter::new(),
            gain_pred: GainPredictor::new(),
            prev_gp_hat: 0.8,
            s_past: [0.0; M],
            taming: Taming::new(),
        }
    }

    /// eq (40): builds the adaptive-codebook vector `v(n)` by
    /// interpolating the excitation buffer at delay `(int_t, frac)`.
    /// `exc` uses the [`crate::closed_loop_pitch::EXC_BUFFER`] layout
    /// (`u(n)` at index `PIT_MAX + n`), with the current-subframe
    /// region already filled with the LP residual (clause 3.7).
    fn adaptive_vector(exc: &[f32; EXC_BUFFER], delay: PitchDelay) -> [f32; L_SUBFR] {
        use crate::tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15 as B30;
        const Q15: f32 = 32_768.0;
        // Map {int_t, frac ∈ {−1,0,1}} onto (k, t ∈ {0,1,2}).
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
                acc += u(n - k - i) * f32::from(B30[(t + 3 * i) as usize]);
                acc += u(n - k + 1 + i) * f32::from(B30[(3 - t + 3 * i) as usize]);
            }
            acc / Q15
        })
    }

    /// Encodes one frame and serialises it straight to the 164-byte
    /// ITU serial wire format (Table-8 MSB-first packing +
    /// sync/bit-count framing) — the full `.IN` → `.BIT` path.
    pub fn encode_frame_to_serial(
        &mut self,
        pcm: &[f32; L_FRAME],
    ) -> [u8; crate::serial::FRAME_BYTES] {
        let out = self.encode_frame(pcm);
        let bits = crate::parameters::pack_bit_array(&out.params);
        crate::serial::write_frame(&bits)
    }

    /// Encodes one 80-sample frame of 16-bit PCM (as `f32` sample
    /// values). Returns the transmitted parameter set. The first call
    /// primes the 5 ms look-ahead, so output frame `k` corresponds to
    /// input samples `80k − 40 … 80k + 39` (clause 3.2.1 timing).
    pub fn encode_frame(&mut self, pcm: &[f32; L_FRAME]) -> EncodedFrame {
        // §3.1 pre-process the new input and slide it into the
        // analysis buffer (Figure 5: new samples fill 160 … 239).
        // Raw input PCM is integral; the §3.1 filter runs on the
        // Word16/Word32 grid and its output is stored back as f32 for
        // the float mid-chain stages.
        let pcm_i: [i16; L_FRAME] = std::array::from_fn(|n| pcm[n] as i16);
        let s_new = self.preproc.process_frame(&pcm_i);
        self.speech.copy_within(L_FRAME.., 0);
        for (o, &v) in self.speech[L_WINDOW - L_FRAME..].iter_mut().zip(&s_new) {
            *o = f32::from(v);
        }

        // §3.2.1/2/3: LP analysis on the clause-5 fixed-point grid
        // (round 438): genuine Word32 autocorrelation with the
        // overflow-rescale protocol, DPF Levinson-Durbin, and the
        // fixed Chebyshev root search. The corpus-pinned configuration
        // (see `tests/fx_front_end_conformance.rs`) beats the float
        // grid emulation on every locked-history vector — most
        // dramatically TAME 80.5 → 99.2%. The float chain remains the
        // fallback for an ill-conditioned frame (never triggered on
        // the corpus) and the spec-equation oracle.
        let speech_i: [i16; L_WINDOW] = std::array::from_fn(|n| self.speech[n] as i16);
        let ac = crate::fx::analysis::analyze_window_fx(
            &speech_i,
            &crate::fx::analysis::FrontEndLatitude::default(),
        );
        let fx_front = crate::fx::levinson::levinson_fx(&ac.r).and_then(|lev_fx| {
            crate::fx::lp_to_lsp::lp_to_lsp_fx(&lev_fx.a_q12).map(|q_fx| (lev_fx, q_fx))
        });
        let (reflection, q_unq) = match fx_front {
            Some((lev_fx, q_fx)) => {
                let q: [f32; M] = std::array::from_fn(|i| f32::from(q_fx[i]) / 32_768.0);
                self.prev_q_unq = q;
                (
                    (
                        f32::from(lev_fx.rc_q15[0]) / 32_768.0,
                        f32::from(lev_fx.rc_q15[1]) / 32_768.0,
                    ),
                    q,
                )
            }
            None => {
                // Defensive float fallback: Q12-rounded Levinson +
                // float root search; on a failed search keep the
                // previous frame's LSPs (clause 3.2.3).
                let r = analyze_window(&self.speech);
                let lev = levinson(&r);
                let a_q12: [f32; M] = std::array::from_fn(|i| (lev.a[i] * 4096.0).round() / 4096.0);
                if let Some(q) = lp_to_lsp(&a_q12) {
                    self.prev_q_unq = q;
                }
                ((lev.reflection[0], lev.reflection[1]), self.prev_q_unq)
            }
        };

        // §3.2.4: quantise (advances the decoder-consistent MA state).
        // eq (18) runs through the fixed-point arccos lookup so the
        // quantiser sees LSFs on the reference's Q13 grid.
        let omega_unq: [f32; M] = std::array::from_fn(|i| {
            crate::lsf_conversion::lsp_to_lsf_q13(q_unq[i]) as f32 / 8192.0
        });
        let lsp_out = self.lsp_quant.quantize_lsf(&omega_unq);
        let q_quant = omega_to_q(&lsp_out.omega_hat);

        // §3.2.5: per-subframe interpolation of both tracks.
        let q_sub_quant = self.interp_q.interpolate(&q_quant);
        let q_sub_unq = self.interp_u.interpolate(&q_unq);

        // Per-subframe LP coefficient sets.
        let a_hat: [[f32; M]; 2] = [lsp_to_lp(&q_sub_quant[0]), lsp_to_lp(&q_sub_quant[1])];
        let a_unq: [[f32; M]; 2] = [lsp_to_lp(&q_sub_unq[0]), lsp_to_lp(&q_sub_unq[1])];

        // §3.3: γ adaptation needs the per-subframe unquantised LSF (in
        // radians).
        let omega_sub: [[f32; M]; 2] = [
            crate::lsp_interpolate::q_to_omega(&q_sub_unq[0]),
            crate::lsp_interpolate::q_to_omega(&q_sub_unq[1]),
        ];
        let decisions = self.weighting.adapt_frame(reflection, &omega_sub);

        // Weighted speech for both subframes (the present frame is
        // buffer 120 … 199).
        let s_frame: [f32; L_FRAME] = std::array::from_fn(|n| self.speech[120 + n]);
        let sub1: [f32; L_SUBFR] = std::array::from_fn(|n| s_frame[n]);
        let sub2: [f32; L_SUBFR] = std::array::from_fn(|n| s_frame[L_SUBFR + n]);
        let sw1 = self
            .weighting
            .weight_subframe(&sub1, &a_unq[0], decisions[0]);
        let sw2 = self
            .weighting
            .weight_subframe(&sub2, &a_unq[1], decisions[1]);

        // §3.4: open-loop pitch over the whole frame's weighted speech.
        let mut sw_buf = [0.0f32; PIT_BUFFER];
        sw_buf[..PIT_MAX].copy_from_slice(&self.sw_hist);
        sw_buf[PIT_MAX..PIT_MAX + L_SUBFR].copy_from_slice(&sw1);
        sw_buf[PIT_MAX + L_SUBFR..].copy_from_slice(&sw2);
        let olp = open_loop_pitch(&sw_buf);
        // Advance the weighted-speech history.
        self.sw_hist.copy_within(L_FRAME.., 0);
        self.sw_hist[PIT_MAX - L_FRAME..PIT_MAX - L_SUBFR].copy_from_slice(&sw1);
        self.sw_hist[PIT_MAX - L_SUBFR..].copy_from_slice(&sw2);

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
        // Fresh excitation this frame (filled subframe by subframe).
        let mut exc_frame = [0.0f32; L_FRAME];

        for sub in 0..2 {
            let s_sub: [f32; L_SUBFR] = std::array::from_fn(|n| s_frame[sub * L_SUBFR + n]);
            let coefs = WeightedSynthesisCoefs {
                a_hat: a_hat[sub],
                a: a_unq[sub],
                gamma1: decisions[sub].gamma1,
                gamma2: decisions[sub].gamma2,
            };

            // §3.6: residual + target.
            let residual = lp_residual(&s_sub, &self.s_past, &a_hat[sub]);
            let x = self.target_filter.target(&residual, &coefs);
            // §3.5: impulse response.
            let h = impulse_response(&coefs);

            // Excitation buffer for the search: history + this frame's
            // already-built excitation + the residual in the current
            // subframe region (clause 3.7).
            let mut exc = [0.0f32; EXC_BUFFER];
            exc[..PIT_MAX].copy_from_slice(&self.exc_hist);
            if sub == 1 {
                // Subframe 1's final excitation shifts into history
                // range: u(−40) … u(−1) of subframe 2 live at the tail
                // of the frame's history. Rebuild: history minus the
                // oldest 40 + subframe-1 excitation.
                exc[..PIT_MAX - L_SUBFR].copy_from_slice(&self.exc_hist[L_SUBFR..]);
                exc[PIT_MAX - L_SUBFR..PIT_MAX].copy_from_slice(&exc_frame[..L_SUBFR]);
            }
            exc[PIT_MAX..PIT_MAX + L_SUBFR].copy_from_slice(&residual);

            // §3.7 closed-loop pitch.
            let (t_min, t_max) = if sub == 0 {
                subframe1_window(olp.t_op as i32)
            } else {
                let t1 = t1_delay.expect("subframe order");
                subframe2_window(t1.int_t)
            };
            let frac_limit = if sub == 0 {
                T1_FRACTIONAL_LIMIT
            } else {
                i32::MAX
            };
            let cl = closed_loop_search(&x, &h, &exc, t_min, t_max, frac_limit);
            delays[sub] = cl.delay;

            // eq (40) adaptive vector + eq (43) gain.
            let v = Self::adaptive_vector(&exc, cl.delay);
            let y = filter_through_h(&v, &h);
            let gp_ol = adaptive_gain(&x, &y);
            // eq (50) fixed-codebook target.
            let x_prime = update_target(&x, &y, gp_ol);
            // eq (49): fold the harmonic pre-filter into h.
            let beta = clamp_beta(self.prev_gp_hat);
            let mut h_pre = h;
            prefilter_impulse_response(&mut h_pre, cl.delay.int_t.max(0) as usize, beta);
            // §3.8.1 search.
            let d = correlation_d(&x_prime, &h_pre);
            let phi = phi_matrix(&h_pre);
            let choice = search_fixed_codebook(&d, &phi, &mut loop4_budget);
            // eq (45) codevector + eq (48) sharpening (the eq (47)
            // clamp is applied inside `sharpen`).
            let mut c_i8 = [0i8; L_SUBFR];
            for (pos, sign) in choice.positions.iter().zip(choice.signs.iter()) {
                c_i8[*pos as usize] += *sign;
            }
            let c = sharpen(&c_i8, cl.delay.int_t, self.prev_gp_hat);
            // eq (64): filtered fixed vector (through the un-prefiltered
            // h — the sharpening is already inside c).
            let z = filter_through_h(&c, &h);

            // §3.9: gains. The taming test bounds the quantised
            // adaptive gain when the long-term loop for this delay has
            // accumulated too much worst-case decoder error
            // (docs/audio/g729/taming-procedure.md §3/§4).
            let pred = self.gain_pred.predict_only(&c);
            let terms = GainTerms::compute(&x, &y, &z);
            let tame = self.taming.test(cl.delay.int_t, cl.delay.frac);
            let gq = quantize_gains(&terms, pred.g_c_prime, tame);
            self.gain_pred.push_quantised_error(gq.gamma_hat);
            self.prev_gp_hat = gq.gp_hat;
            gains[sub] = (gq.gp_hat, gq.gc_hat);
            // Taming update: propagate the worst-case error one pitch
            // period through the loop at the chosen gain (doc §2).
            self.taming.update(cl.delay.int_t, cl.delay.frac, gq.gp_hat);

            // §3.10: excitation + memory updates.
            let u: [f32; L_SUBFR] = std::array::from_fn(|n| gq.gp_hat * v[n] + gq.gc_hat * c[n]);
            exc_frame[sub * L_SUBFR..(sub + 1) * L_SUBFR].copy_from_slice(&u);
            self.target_filter.update(&residual, &u, &coefs);
            self.s_past = std::array::from_fn(|i| s_sub[L_SUBFR - 1 - i]);

            // Transmitted codewords for this subframe.
            let c_code = encode_positions(&choice.positions).expect("track-valid positions");
            let s_code = encode_signs(&choice.signs).expect("±1 signs");
            let ga_tx = map_ga(gq.ga).expect("in-range") as u8;
            let gb_tx = map_gb(gq.gb).expect("in-range") as u8;
            if sub == 0 {
                let p1 = encode_p1(cl.delay).expect("window-bounded T1");
                params.p1 = p1;
                // §3.7.2 parity: odd parity over the six MSBs (the
                // corpus-pinned convention in `Parameters::pitch_parity_ok`).
                params.p0 = (((p1 >> 2).count_ones() ^ 1) & 1) as u8;
                params.c1 = c_code;
                params.s1 = s_code;
                params.ga1 = ga_tx;
                params.gb1 = gb_tx;
                t1_delay = Some(cl.delay);
            } else {
                let t1 = t1_delay.expect("subframe order");
                let t_min2 = derive_t_min(t1.int_t);
                params.p2 = encode_p2(cl.delay, t_min2).expect("window-bounded T2");
                params.c2 = c_code;
                params.s2 = s_code;
                params.ga2 = ga_tx;
                params.gb2 = gb_tx;
            }
        }

        // Advance the excitation history by one frame.
        self.exc_hist.copy_within(L_FRAME.., 0);
        self.exc_hist[PIT_MAX - L_FRAME..].copy_from_slice(&exc_frame);

        EncodedFrame {
            params,
            gains,
            delays,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Feeds a synthetic voiced-like tone through the encoder and
    /// checks every emitted codeword is inside its Table-8 field
    /// domain, the parity predicate holds, and the state machine
    /// survives many frames.
    #[test]
    fn codewords_stay_in_field_domains() {
        let mut enc = FrameEncoder::new();
        for f in 0..20 {
            let pcm: [f32; L_FRAME] = std::array::from_fn(|n| {
                let t = (f * L_FRAME + n) as f32;
                6000.0 * (std::f32::consts::TAU * t / 57.0).sin()
                    + 1500.0 * (std::f32::consts::TAU * t / 19.0).sin()
            });
            let out = enc.encode_frame(&pcm);
            let p = &out.params;
            assert!(p.l0 < 2 && p.l1 < 128 && p.l2 < 32 && p.l3 < 32);
            assert!(p.p0 < 2 && p.p2 < 32);
            assert!(p.c1 < 8192 && p.c2 < 8192);
            assert!(p.s1 < 16 && p.s2 < 16);
            assert!(p.ga1 < 8 && p.ga2 < 8 && p.gb1 < 16 && p.gb2 < 16);
            assert!(
                p.pitch_parity_ok(),
                "frame {f}: parity predicate must accept our own P0"
            );
            for (gp, gc) in out.gains {
                assert!(gp.is_finite() && gc.is_finite());
            }
        }
    }

    /// Silence encodes without NaN/panic and produces valid codewords.
    #[test]
    fn silence_is_stable() {
        let mut enc = FrameEncoder::new();
        let pcm = [0.0f32; L_FRAME];
        for _ in 0..5 {
            let out = enc.encode_frame(&pcm);
            assert!(out.params.pitch_parity_ok());
            for (gp, gc) in out.gains {
                assert!(gp.is_finite() && gc.is_finite());
            }
        }
    }

    /// End-to-end encode→decode: the decoder must accept every frame
    /// the encoder emits (parameter-domain validity) and reconstruct
    /// speech whose energy tracks the input's onset (a structural
    /// closed-loop sanity check, not yet a conformance bound).
    #[test]
    fn encode_decode_loop_is_consistent() {
        use crate::decode_chain::FrameDecoder;

        let mut enc = FrameEncoder::new();
        let mut dec = FrameDecoder::new();
        let mut out_energy = 0.0f64;
        for f in 0..25 {
            let pcm: [f32; L_FRAME] = std::array::from_fn(|n| {
                let t = (f * L_FRAME + n) as f32;
                5000.0 * (std::f32::consts::TAU * t / 63.0).sin()
            });
            let e = enc.encode_frame(&pcm);
            let d = dec
                .decode_parameters_to_speech(&e.params)
                .expect("decoder accepts encoder output");
            for s in d.speech() {
                out_energy += f64::from(s) * f64::from(s);
            }
        }
        assert!(
            out_energy > 1.0e6,
            "decoded speech should carry real energy, got {out_energy}"
        );
    }

    /// The encoder's local decoded gains match what a decoder
    /// reconstructs from the emitted (GA, GB) codewords — pinning the
    /// map/demap + reconstruction lock-step.
    #[test]
    fn gains_are_decoder_consistent() {
        use crate::gain_index_map::{demap_ga, demap_gb};
        use crate::gain_reconstruct::reconstruct_gains;

        let mut enc = FrameEncoder::new();
        let pcm: [f32; L_FRAME] =
            std::array::from_fn(|n| 4000.0 * (std::f32::consts::TAU * (n as f32) / 41.0).sin());
        let out = enc.encode_frame(&pcm);
        let ga = demap_ga(out.params.ga1 as usize).unwrap();
        let gb = demap_gb(out.params.gb1 as usize).unwrap();
        let q = reconstruct_gains(ga, gb).unwrap();
        assert!((q.g_p_hat - out.gains[0].0).abs() < 1e-6);
    }
}
