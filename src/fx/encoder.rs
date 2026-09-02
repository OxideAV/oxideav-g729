//! Clause-3 **fixed-point frame encoder** — the encoder analysis
//! chain driven on the clause-5 Word16/Word32 grid, sharing every
//! decoder-side primitive with [`crate::fx::decoder::FrameDecoderFx`]
//! so the encoder's local reconstruction is the decoder's, bit for
//! bit.
//!
//! ## Stage map
//!
//! | clause | stage | module |
//! |---|---|---|
//! | §3.1 | pre-processing | [`crate::fx::analysis::PreprocessorFx`] |
//! | §3.2.1–3.2.3 | LP analysis, Levinson, LP → LSP | `fx::analysis` / `fx::levinson` / `fx::lp_to_lsp` |
//! | §3.2.4 | LSP quantisation | [`crate::lsp_quantize::search_lsp_indices_q13`] |
//! | §3.2.5–3.2.6 | interpolation, LSP → LP | [`crate::fx::lsp`] (local decoder) |
//! | §3.3 | perceptual weighting | [`crate::fx::weighting`] |
//! | §3.4 | open-loop pitch | [`crate::fx::pitch_ol`] |
//! | §3.5–3.6 | impulse response, residual, target | [`crate::fx::target`] |
//! | §3.7 | closed-loop pitch + gain | [`crate::fx::pitch_cl`] |
//! | §3.8 | algebraic codebook search | [`crate::fx::acelp`] |
//! | §3.9 | gain quantisation + taming | [`crate::fx::gain_vq`] |
//! | §3.10 | memory update | this module (decoder primitives) |
//!
//! ## Locked-history drive
//!
//! [`FrameEncoderFx::encode_frame_locked`] runs every stage's search
//! exactly as the free encoder does but commits the **reference's**
//! transmitted parameters to the local decoder state afterwards
//! (LSP MA history, excitation, gain predictor, synthesis and
//! weighting memories). Against a conformance `.BIT` stream that puts
//! the encoder in the state the reference encoder had at every frame,
//! so each stage's decision can be measured in isolation from the
//! error propagation of earlier misses.

use crate::fixed_codebook::{decode_pulses, encode_positions, encode_signs};
use crate::fx::acelp::search_acelp_fx;
use crate::fx::analysis::{analyze_window_fx, FrontEndLatitude, PreprocessorFx};
use crate::fx::excitation::{
    build_codevector_q13, build_excitation_mode, pred_lt3, syn_filt_check, EXC_BUF, PIT_MAX,
    SHARP_INIT_Q14, SHARP_MAX_Q14, SHARP_MIN_Q14,
};
use crate::fx::filters::{convolve_code_q12, convolve_h_q0, residu, weight_az};
use crate::fx::gain_vq::{quantize_gains_fx, GainCorrelationsFx, TamingFx};
use crate::fx::gains::{DecodedGainsFx, GainDecoderFx};
use crate::fx::levinson::levinson_fx;
use crate::fx::lp_to_lsp::lp_to_lsp_fx;
use crate::fx::lsp::{interpolate_lsp_q15, lsp_to_lp_q12, LspDecoderFx, STARTUP_LSP_Q15};
use crate::fx::ops::{l_mac, l_mult, l_shl, round};
use crate::fx::pitch_cl::{
    closed_loop_search_fx, pitch_gain_q14, subframe1_window, subframe2_window, ClosedLoopLatitude,
    T1_FRACTIONAL_LIMIT,
};
use crate::fx::pitch_ol::{open_loop_pitch_fx, PIT_BUFFER};
use crate::fx::target::{impulse_response_fx, TargetChainFx};
use crate::fx::weighting::{WeightingFx, WeightingLatitude};
use crate::gain_index_map::{demap_ga, demap_gb, map_ga, map_gb};
use crate::lsf_conversion::acos_q15_to_lsf_q13;
use crate::lsp_quantize::{
    advance_history_q13, search_lsp_indices_q13, startup_history_q13, FxLatitude,
};
use crate::parameters::Parameters;
use crate::pitch_decode::{
    decode_t1_from_p1, decode_t2_from_p2, derive_t_min, encode_p1, encode_p2, PitchDelay,
};
use crate::tables::{L_WINDOW, M, MA_NP};

/// Samples per frame.
pub const L_FRAME: usize = 80;
/// Samples per subframe.
pub const L_SUBFR: usize = 40;
/// History region kept ahead of the working frame in the excitation
/// buffer (maximum delay + interpolator margin) — the decoder's layout.
const EXC_HIST: usize = PIT_MAX + 10;
/// Offset of the present frame inside the Figure-5 analysis buffer.
const FRAME_OFFSET: usize = 120;

/// One encoded frame: the 15 Table-8 codewords plus the encoder's
/// local decoded side data.
#[derive(Debug, Clone, Copy)]
pub struct EncodedFrameFx {
    /// The transmitted parameter set.
    pub params: Parameters,
    /// Per-subframe selected pitch delays.
    pub delays: [PitchDelay; 2],
    /// Per-subframe quantised gains the local decoder used.
    pub gains: [DecodedGainsFx; 2],
    /// The §3.4 open-loop delay.
    pub t_op: i32,
}

/// Per-subframe intermediate signals of the last encoded frame — the
/// bisection instrument for the stage-isolation harness.
#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct SubframeProbe {
    /// §3.6 target `x(n)` (Q0).
    pub x: [i16; L_SUBFR],
    /// eq (50) fixed-codebook target `x′(n)` (Q0).
    pub x_prime: [i16; L_SUBFR],
    /// §3.5 impulse response (Q12).
    pub h: [i16; L_SUBFR],
    /// eq (49) pre-filtered impulse response (Q12).
    pub h_pre: [i16; L_SUBFR],
    /// eq (44) filtered adaptive vector at the committed delay (Q0).
    pub y: [i16; L_SUBFR],
    /// eq (64) filtered committed codevector (Q12).
    pub z_q12: [i16; L_SUBFR],
    /// eq (43) gain (Q14).
    pub gp_q14: i16,
    /// eq (71) prediction the gain search used.
    pub pred: crate::fx::gains::GainPredictionFx,
    /// Taming flag the gain search used.
    pub tame: bool,
    /// The §3.3 `(γ₁, γ₂)` pair (Q15) the subframe ran with.
    pub gamma_q15: (i16, i16),
}

impl Default for SubframeProbe {
    fn default() -> Self {
        Self {
            x: [0; L_SUBFR],
            x_prime: [0; L_SUBFR],
            h: [0; L_SUBFR],
            h_pre: [0; L_SUBFR],
            y: [0; L_SUBFR],
            z_q12: [0; L_SUBFR],
            gp_q14: 0,
            pred: crate::fx::gains::GainPredictionFx {
                g_prime_scaled: 0,
                exp: 0,
            },
            tame: false,
            gamma_q15: (0, 0),
        }
    }
}

/// Unstated fixed-point latitude of the encoder mid-chain, pinned
/// black-box against the conformance corpus.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EncoderLatitude {
    /// §3.7 closed-loop search latitude.
    pub cl: ClosedLoopLatitude,
    /// Force a `(γ₁, γ₂)` pair (Q15) on every subframe instead of the
    /// §3.3 adaptation (sensitivity probe).
    pub gamma_override: Option<(i16, i16)>,
    /// §3.3 decision latitude.
    pub weighting: WeightingLatitude,
}

/// The stateful fixed-point G.729 frame encoder.
#[derive(Debug, Clone)]
pub struct FrameEncoderFx {
    preproc: PreprocessorFx,
    /// Figure-5 analysis buffer (`0 … 119` past, `120 … 199` present,
    /// `200 … 239` look-ahead), pre-processed speech on Q0.
    speech: [i16; L_WINDOW],
    /// Previous frame's unquantised LSPs (Q15) — the §3.2.5
    /// interpolation memory of the unquantised track and the §3.2.3
    /// root-search fallback.
    prev_lsp_unq: [i16; M],
    /// §3.2.4 MA history driving the index search.
    history_q13: [[i16; M]; MA_NP],
    /// The local §4.1.1 LP-parameter decoder (quantised track).
    lsp_dec: LspDecoderFx,
    /// §3.3 weighting state.
    weighting: WeightingFx,
    /// Weighted-speech history + present frame for the §3.4 search.
    sw_buf: [i16; PIT_BUFFER],
    /// Excitation buffer in the decoder's layout.
    exc: [i16; EXC_BUF],
    /// Local synthesis memory (`ŝ(n−1) … ŝ(n−M)`, most recent last).
    syn_mem: [i16; M],
    /// §3.6 / §3.10 target-chain memories (`e(n)`, `ew(n)`).
    target: TargetChainFx,
    /// §3.9.1 gain predictor (decoder-consistent).
    gains: GainDecoderFx,
    /// eq (47)/(49) sharpening gain β (Q14) from the previous
    /// subframe's quantised pitch gain.
    sharp_q14: i16,
    /// Taming-procedure state.
    taming: TamingFx,
    /// Fixed-point latitude.
    lat: EncoderLatitude,
    /// §3.8.1 fourth-loop budget for the current frame.
    loop4_budget: u32,
    /// Last frame's per-subframe intermediates (instrument).
    #[doc(hidden)]
    pub probe: [SubframeProbe; 2],
}

impl Default for FrameEncoderFx {
    fn default() -> Self {
        Self::new()
    }
}

impl FrameEncoderFx {
    /// Fresh encoder in the clause-4.3 start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self::with_latitude(EncoderLatitude::default())
    }

    /// Fresh encoder under an explicit latitude (black-box sweeps).
    #[doc(hidden)]
    #[must_use]
    pub fn with_latitude(lat: EncoderLatitude) -> Self {
        Self {
            preproc: PreprocessorFx::new(),
            speech: [0; L_WINDOW],
            prev_lsp_unq: STARTUP_LSP_Q15,
            history_q13: startup_history_q13(),
            lsp_dec: LspDecoderFx::new(),
            weighting: WeightingFx::with_latitude(lat.weighting),
            sw_buf: [0; PIT_BUFFER],
            exc: [0; EXC_BUF],
            syn_mem: [0; M],
            target: TargetChainFx::new(),
            gains: GainDecoderFx::new(),
            sharp_q14: SHARP_INIT_Q14,
            taming: TamingFx::new(),
            lat,
            loop4_budget: 0,
            probe: [SubframeProbe::default(); 2],
        }
    }

    /// Encodes one 80-sample frame of 16-bit PCM. The first call
    /// primes the 5 ms look-ahead, so output frame `k` corresponds to
    /// input samples `80k − 40 … 80k + 39` (clause 3.2.1 timing).
    pub fn encode_frame(&mut self, pcm: &[i16; L_FRAME]) -> EncodedFrameFx {
        self.encode_frame_inner(pcm, None)
    }

    /// Encodes one frame and serialises it straight to the 164-byte
    /// ITU serial wire format.
    pub fn encode_frame_to_serial(
        &mut self,
        pcm: &[i16; L_FRAME],
    ) -> [u8; crate::serial::FRAME_BYTES] {
        let out = self.encode_frame(pcm);
        let bits = crate::parameters::pack_bit_array(&out.params);
        crate::serial::write_frame(&bits)
    }

    /// Locked-history drive: runs every search, returns the encoder's
    /// own decisions, but commits `reference`'s parameters to the local
    /// decoder state (see the module docs).
    #[doc(hidden)]
    pub fn encode_frame_locked(
        &mut self,
        pcm: &[i16; L_FRAME],
        reference: &Parameters,
    ) -> EncodedFrameFx {
        self.encode_frame_inner(pcm, Some(reference))
    }

    fn encode_frame_inner(
        &mut self,
        pcm: &[i16; L_FRAME],
        lock: Option<&Parameters>,
    ) -> EncodedFrameFx {
        // §3.1: pre-process and slide into the analysis buffer.
        let s_new = self.preproc.process_frame(pcm);
        self.speech.copy_within(L_FRAME.., 0);
        self.speech[L_WINDOW - L_FRAME..].copy_from_slice(&s_new);

        // §3.2.1–§3.2.3 on the fixed grid; a failed root search keeps
        // the previous frame's LSPs (clause 3.2.3).
        let ac = analyze_window_fx(&self.speech, &FrontEndLatitude::default());
        let (a_unq_q12, rc_q15, lsp_unq) = match levinson_fx(&ac.r) {
            Some(lev) => {
                let lsp = lp_to_lsp_fx(&lev.a_q12).unwrap_or(self.prev_lsp_unq);
                let mut a = [0i16; M + 1];
                a[0] = 4096;
                a[1..].copy_from_slice(&lev.a_q12);
                (a, [lev.rc_q15[0], lev.rc_q15[1]], lsp)
            }
            None => {
                // Ill-conditioned frame: reuse the previous LSPs and a
                // flat reflection pair.
                (lsp_to_lp_q12(&self.prev_lsp_unq), [0, 0], self.prev_lsp_unq)
            }
        };

        // §3.2.4: eq (18) through the fixed arccos lookup, Q13 search.
        let lat_q = FxLatitude::default();
        let lsf_unq_q13: [i32; M] =
            std::array::from_fn(|i| acos_q15_to_lsf_q13(i32::from(lsp_unq[i])));
        let idx = search_lsp_indices_q13(&lsf_unq_q13, &self.history_q13, &lat_q);
        let (l0, l1, l2, l3) = match lock {
            Some(r) => (
                usize::from(r.l0),
                usize::from(r.l1),
                usize::from(r.l2),
                usize::from(r.l3),
            ),
            None => (idx.l0, idx.l1, idx.l2, idx.l3),
        };
        advance_history_q13(&mut self.history_q13, l1, l2, l3, &lat_q);
        // §3.2.5–§3.2.6 quantised track through the local decoder.
        let lp = self.lsp_dec.decode_frame(l0, l1, l2, l3);
        let a_hat: [[i16; M + 1]; 2] = [lp.a_sub1, lp.a_sub2];

        // Unquantised track: cosine-domain interpolation for subframe
        // 1, the current LSPs for subframe 2.
        let lsp_sub1 = interpolate_lsp_q15(&self.prev_lsp_unq, &lsp_unq);
        self.prev_lsp_unq = lsp_unq;
        let a_unq: [[i16; M + 1]; 2] = [lsp_to_lp_q12(&lsp_sub1), a_unq_q12];
        let lsf_sub: [[i16; M]; 2] = [
            std::array::from_fn(|i| acos_q15_to_lsf_q13(i32::from(lsp_sub1[i])) as i16),
            std::array::from_fn(|i| lsf_unq_q13[i] as i16),
        ];

        // §3.3: γ adaptation + weighted speech for both subframes.
        let mut decisions = self.weighting.adapt_frame(rc_q15, &lsf_sub);
        if let Some((g1, g2)) = self.lat.gamma_override {
            for d in &mut decisions {
                d.gamma1_q15 = g1;
                d.gamma2_q15 = g2;
            }
        }
        let mut ap1 = [[0i16; M + 1]; 2];
        let mut ap2 = [[0i16; M + 1]; 2];
        let mut sw_frame = [0i16; L_FRAME];
        for sub in 0..2 {
            ap1[sub] = weight_az(&a_unq[sub], decisions[sub].gamma1_q15);
            ap2[sub] = weight_az(&a_unq[sub], decisions[sub].gamma2_q15);
            let base = FRAME_OFFSET + sub * L_SUBFR;
            let s_hist: [i16; M] = std::array::from_fn(|i| self.speech[base - M + i]);
            let s_sub: [i16; L_SUBFR] = std::array::from_fn(|n| self.speech[base + n]);
            let sw = self
                .weighting
                .weight_subframe(&s_hist, &s_sub, &ap1[sub], &ap2[sub]);
            sw_frame[sub * L_SUBFR..(sub + 1) * L_SUBFR].copy_from_slice(&sw);
        }

        // §3.4: open-loop pitch over the frame's weighted speech.
        self.sw_buf[PIT_MAX..].copy_from_slice(&sw_frame);
        let olp = open_loop_pitch_fx(&self.sw_buf);
        self.sw_buf.copy_within(L_FRAME.., 0);

        let mut params = Parameters {
            l0: idx.l0 as u8,
            l1: idx.l1 as u8,
            l2: idx.l2 as u8,
            l3: idx.l3 as u8,
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
        let mut delays = [PitchDelay { int_t: 0, frac: 0 }; 2];
        let mut gains = [DecodedGainsFx {
            gain_pit_q14: 0,
            gain_code_q1: 0,
        }; 2];
        self.loop4_budget = crate::fx::acelp::MAX_LOOP4_PER_FRAME;
        // The committed subframe-1 delay (own or reference) anchors
        // the subframe-2 window.
        let mut t1_committed: Option<PitchDelay> = None;

        for sub in 0..2 {
            let off = EXC_HIST + sub * L_SUBFR;
            let base = FRAME_OFFSET + sub * L_SUBFR;
            let s_hist: [i16; M] = std::array::from_fn(|i| self.speech[base - M + i]);
            let s_sub: [i16; L_SUBFR] = std::array::from_fn(|n| self.speech[base + n]);

            // §3.5: impulse response of W(z)/Â(z) on Q12.
            let h = impulse_response_fx(&a_hat[sub], &ap1[sub], &ap2[sub]);
            // §3.6: residual + target.
            let res = residu(&a_hat[sub], &s_hist, &s_sub);
            let x = self.target.target(&res, &a_hat[sub], &ap1[sub], &ap2[sub]);

            // Extend the excitation with the residual (clause 3.7).
            self.exc[off..off + L_SUBFR].copy_from_slice(&res);

            // §3.7: closed-loop pitch.
            let (t_min, t_max, frac_limit) = if sub == 0 {
                let (lo, hi) = subframe1_window(olp.t_op);
                (lo, hi, T1_FRACTIONAL_LIMIT)
            } else {
                let t1 = t1_committed.expect("subframe order");
                let (lo, hi) = subframe2_window(t1.int_t);
                (lo, hi, i32::MAX)
            };
            let cl = closed_loop_search_fx(
                &x,
                &h,
                &self.exc,
                off,
                t_min,
                t_max,
                frac_limit,
                &self.lat.cl,
            );
            let own_delay = cl.delay;

            // The committed delay for this subframe.
            let delay = match lock {
                Some(r) => {
                    if sub == 0 {
                        decode_t1_from_p1(r.p1)
                    } else {
                        let t1 = t1_committed.expect("subframe order");
                        decode_t2_from_p2(r.p2, derive_t_min(t1.int_t))
                    }
                }
                None => own_delay,
            };
            if sub == 0 {
                t1_committed = Some(delay);
            }
            delays[sub] = own_delay;

            // eq (40) adaptive vector (in place, decoder geometry) at
            // the committed delay, eq (44) filtered vector, eq (43)
            // gain. Under lock every downstream stage runs from the
            // reference's upstream decisions, so each stage's own
            // decision is measured in isolation.
            pred_lt3(&mut self.exc, off, delay.int_t, delay.frac);
            let v: [i16; L_SUBFR] = std::array::from_fn(|n| self.exc[off + n]);
            let y = convolve_h_q0(&v, &h);
            let gp_q14 = pitch_gain_q14(&x, &y);

            // eq (50) target update + eq (49) pre-filter on h.
            let x_prime: [i16; L_SUBFR] = std::array::from_fn(|n| {
                let acc = l_mult(y[n], gp_q14); // Q15
                crate::fx::ops::sub(x[n], round(l_shl(acc, 1)))
            });
            let mut h_pre = h;
            prefilter_h(&mut h_pre, delay.int_t, self.sharp_q14);

            // §3.8.1 search.
            let choice = search_acelp_fx(&x_prime, &h_pre, &mut self.loop4_budget);

            // Committed codevector.
            let (positions, signs) = match lock {
                Some(r) => {
                    let (c, s) = if sub == 0 { (r.c1, r.s1) } else { (r.c2, r.s2) };
                    let pulses = decode_pulses(c, s).expect("reference fields in domain");
                    let positions: [u8; 4] = std::array::from_fn(|k| pulses.pulses[k].position);
                    let signs: [i8; 4] = std::array::from_fn(|k| pulses.pulses[k].sign);
                    (positions, signs)
                }
                None => (choice.positions, choice.signs),
            };

            // §3.9: gains on the committed codevector (sharpened at the
            // committed delay).
            let code = build_codevector_q13(&positions, &signs, self.sharp_q14, delay.int_t);
            let z = convolve_code_q12(&code, &h);
            let pred = self.gains.predict(&code);
            let corr = GainCorrelationsFx::compute(&x, &y, &z);
            let tame = self.taming.test(delay.int_t, delay.frac);
            let gq = quantize_gains_fx(&corr, &self.gains, &pred, tame);

            // Committed gains + state.
            let (ga, gb) = match lock {
                Some(r) => {
                    let (tga, tgb) = if sub == 0 {
                        (r.ga1, r.gb1)
                    } else {
                        (r.ga2, r.gb2)
                    };
                    (
                        demap_ga(usize::from(tga)).expect("3-bit field"),
                        demap_gb(usize::from(tgb)).expect("4-bit field"),
                    )
                }
                None => (gq.ga, gq.gb),
            };
            let g = self.gains.decode(ga, gb, &code);
            gains[sub] = g;
            self.probe[sub] = SubframeProbe {
                x,
                x_prime,
                h,
                h_pre,
                y,
                z_q12: z,
                gp_q14,
                pred,
                tame,
                gamma_q15: (decisions[sub].gamma1_q15, decisions[sub].gamma2_q15),
            };

            // eq (48)/(47): next subframe's sharpening gain.
            self.sharp_q14 = g.gain_pit_q14.clamp(SHARP_MIN_Q14, SHARP_MAX_Q14);
            self.taming.update(delay.int_t, delay.frac, g.gain_pit_q14);

            // §3.10: excitation, local synthesis, filter memories.
            build_excitation_mode(&mut self.exc, off, &code, g.gain_pit_q14, g.gain_code_q1, 3);
            let u: [i16; L_SUBFR] = std::array::from_fn(|n| self.exc[off + n]);
            let mut s_hat = [0i16; L_SUBFR];
            let ovf = syn_filt_check(&u, &a_hat[sub], &mut self.syn_mem, &mut s_hat, false);
            if ovf {
                for slot in self.exc.iter_mut() {
                    *slot = crate::fx::ops::shr(*slot, 2);
                }
                for slot in self.syn_mem.iter_mut() {
                    *slot = crate::fx::ops::shr(*slot, 2);
                }
                let u2: [i16; L_SUBFR] = std::array::from_fn(|n| self.exc[off + n]);
                let _ = syn_filt_check(&u2, &a_hat[sub], &mut self.syn_mem, &mut s_hat, true);
            } else {
                let _ = syn_filt_check(&u, &a_hat[sub], &mut self.syn_mem, &mut s_hat, true);
            }
            self.target
                .update(&s_sub, &s_hat, &x, &y, &z, g.gain_pit_q14, g.gain_code_q1);

            // Transmitted codewords (the encoder's own decisions).
            let ga_tx = map_ga(gq.ga).expect("in-range") as u8;
            let gb_tx = map_gb(gq.gb).expect("in-range") as u8;
            if sub == 0 {
                let p1 = encode_p1(own_delay).expect("window-bounded T1");
                params.p1 = p1;
                params.p0 = (((p1 >> 2).count_ones() ^ 1) & 1) as u8;
                params.c1 = encode_positions(&choice.positions).expect("track-valid positions");
                params.s1 = encode_signs(&choice.signs).expect("±1 signs");
                params.ga1 = ga_tx;
                params.gb1 = gb_tx;
            } else {
                let t1 = t1_committed.expect("subframe order");
                // The own T2 is encoded relative to the committed T1's
                // window (the free encoder commits its own T1).
                params.p2 = encode_p2(own_delay, derive_t_min(t1.int_t)).unwrap_or(0);
                params.c2 = encode_positions(&choice.positions).expect("track-valid positions");
                params.s2 = encode_signs(&choice.signs).expect("±1 signs");
                params.ga2 = ga_tx;
                params.gb2 = gb_tx;
            }
        }

        // Advance the excitation history by one frame.
        self.exc.copy_within(2 * L_SUBFR.., 0);
        for slot in &mut self.exc[EXC_HIST..] {
            *slot = 0;
        }

        EncodedFrameFx {
            params,
            delays,
            gains,
            t_op: olp.t_op,
        }
    }

    /// eq (75) excitation of the last encoded frame (bisection
    /// instrument).
    #[doc(hidden)]
    #[must_use]
    pub fn last_frame_excitation(&self) -> [i16; L_FRAME] {
        std::array::from_fn(|n| self.exc[EXC_HIST - L_FRAME + n])
    }
}

/// eq (49): folds the harmonic pre-filter into the Q12 impulse
/// response, `h(n) += β·h(n−T)` for `n ≥ T` (only when `T < 40`), with
/// `β` on Q14 — the encoder-side mirror of the decoder's eq (48)
/// codevector sharpening.
pub fn prefilter_h(h: &mut [i16; L_SUBFR], t: i32, beta_q14: i16) {
    if t <= 0 || t as usize >= L_SUBFR {
        return;
    }
    let t = t as usize;
    for n in t..L_SUBFR {
        // Q12 × Q14 → Q26 doubled = Q27; back to Q12 with one left
        // shift and a rounding high-word extraction.
        let contrib = round(l_shl(l_mac(0, h[n - t], beta_q14), 1));
        h[n] = crate::fx::ops::add(h[n], contrib);
    }
}
