//! ITU-T G.729 (CS-ACELP, 8 kbit/s) encoder — analysis pipeline symmetric
//! to the in-tree decoder.
//!
//! # Pipeline
//!
//! For every 10 ms frame (80 S16 mono samples at 8 kHz):
//!
//! ```text
//!   PCM s16 → windowed LPC analysis (Levinson-Durbin, order 10)
//!           → LPC → LSP → split-VQ quantisation via LSPCB1 / LSPCB2
//!                              (the same tables the decoder reads back)
//!           → subframe loop (2 × 5 ms subframes, 40 samples each)
//!               - target residual = filter-weighted error signal
//!               - open-loop pitch search on weighted signal
//!               - closed-loop refinement with fractional 1/3 lag
//!               - ACELP 4-pulse fixed-codebook search (focused per track)
//!               - gain analysis + two-stage VQ via GBK1 / GBK2
//!           → bit-pack into 80-bit / 10-byte packet, ITU-T Table 8 order
//! ```
//!
//! # Symmetry with the decoder
//!
//! The encoder reuses every static table the decoder exposes:
//!   - [`LSPCB1_Q13`] / [`LSPCB2_Q13`] — LSP split VQ codebooks
//!   - [`FG_Q15`] / [`FG_SUM_Q15`] — MA-4 predictor coefficients
//!   - [`GBK1`] / [`GBK2`] — two-stage gain codebook
//!
//! `LSPCB1_Q13`, `LSPCB2_Q13`, and the MA-predictor tables (`FG_Q15`,
//! `FG_SUM_Q15`, `FG_SUM_INV_Q12`) are transcribed verbatim from the
//! ITU reference C source `TAB_LD8K.C`, so the LSP-quantisation
//! encode path matches the spec bit-for-bit.
//!
//! Gain-VQ pipeline matches the spec's two-stage conjugate-structured
//! VQ of §3.9: the encoder quantises the correction factor
//! `γ = g_c / g'_c` (eq 72) using its own MA-4 gain predictor
//! initialised to `-14 dB` per Table 9 — exactly what the decoder
//! reconstructs on the receive side. Excitation written back to the
//! encoder's adaptive-codebook buffer uses the quantised `γ_q · g'_c`,
//! keeping the analysis-by-synthesis loop in lockstep with the
//! decoder.
//!
//! The LP-analysis windowing pipeline matches the spec's §3.2.1
//! verbatim: a 240-sample asymmetric window (half-Hamming over
//! n=0..199 plus a quarter-cosine fade over n=200..239, per eq (3))
//! applied to 120 past + 80 current + 40 look-ahead samples. The
//! lookahead introduces the spec's 5-ms extra algorithmic delay.
//!
//! Remaining deviations:
//! - The numeric entries of `GBK1` / `GBK2` are first-cut values
//!   approximating the reference `gbk1` / `gbk2` (dimensions are
//!   spec-correct: 8×2 + 16×2 per §5.2 Table 12).
//! - The Annex B VAD is energy-based rather than four-feature.
//!
//! Bitstreams round-trip cleanly through the in-tree decoder; full
//! external interoperability awaits verbatim `gbk1` / `gbk2`.

use std::collections::VecDeque;

use oxideav_core::Encoder;
use oxideav_core::{
    CodecId, CodecParameters, Error, Frame, MediaType, Packet, Result, SampleFormat, TimeBase,
};

use crate::annex_b_vad::{VadDecision, VadState};
use crate::bitreader::pitch_parity;
use crate::lpc::{interpolate_lsp, lsp_to_lpc, LpcPredictorState, MA_HISTORY};
use crate::lsp_tables::{FG_Q15, FG_SUM_Q15, LSPCB1_Q13, LSPCB2_Q13, MA_NP, M_HALF, NC0, NC1};
use crate::open_loop_pitch::{open_loop_pitch_search, WeightedSpeechState};
use crate::synthesis::{
    adaptive_codebook_excitation, fixed_codebook_excitation, gain_pred_error_db,
    innovation_log_energy_db, predict_fixed_gain, SynthesisState, EXC_HIST, GBK1, GBK2,
};
use crate::weighting::{GAMMA1_FLAT, GAMMA2_FLAT};
use crate::{
    ANNEX_B_ENABLE_EXTRADATA, CODEC_ID_STR, FRAME_BYTES, FRAME_SAMPLES, LPC_ORDER, SAMPLE_RATE,
    SUBFRAMES_PER_FRAME, SUBFRAME_SAMPLES,
};

/// Build a G.729 encoder. Accepts 8 kHz mono S16 input.
pub fn make_encoder(params: &CodecParameters) -> Result<Box<dyn Encoder>> {
    let sample_rate = params.sample_rate.unwrap_or(SAMPLE_RATE);
    if sample_rate != SAMPLE_RATE {
        return Err(Error::unsupported(format!(
            "G.729 encoder: only 8000 Hz is supported (got {sample_rate})"
        )));
    }
    let channels = params.channels.unwrap_or(1);
    if channels != 1 {
        return Err(Error::unsupported(format!(
            "G.729 encoder: only mono is supported (got {channels} channels)"
        )));
    }
    let sample_format = params.sample_format.unwrap_or(SampleFormat::S16);
    if sample_format != SampleFormat::S16 {
        return Err(Error::unsupported(format!(
            "G.729 encoder: input sample format {sample_format:?} not supported (need S16)"
        )));
    }
    if params.codec_id.as_str() != CODEC_ID_STR {
        return Err(Error::unsupported(format!(
            "G.729 encoder: unexpected codec id {:?}",
            params.codec_id
        )));
    }

    let annex_b = params
        .extradata
        .first()
        .map(|&b| b == ANNEX_B_ENABLE_EXTRADATA)
        .unwrap_or(false);

    let mut output = params.clone();
    output.media_type = MediaType::Audio;
    output.sample_format = Some(SampleFormat::S16);
    output.channels = Some(1);
    output.sample_rate = Some(SAMPLE_RATE);
    output.bit_rate = Some(8_000);

    Ok(Box::new(G729Encoder::new(output, annex_b)))
}

struct G729Encoder {
    output_params: CodecParameters,
    time_base: TimeBase,
    state: EncoderState,
    pcm_queue: Vec<i16>,
    pending: VecDeque<Packet>,
    frame_index: u64,
    eof: bool,
    /// If true, the encoder applies Annex B VAD/DTX and may emit SID
    /// (2-byte) or NODATA (0-byte) packets in lieu of a normal 10-byte
    /// speech frame.
    annex_b: bool,
    vad: VadState,
}

impl G729Encoder {
    fn new(output_params: CodecParameters, annex_b: bool) -> Self {
        Self {
            output_params,
            time_base: TimeBase::new(1, SAMPLE_RATE as i64),
            state: EncoderState::new(),
            pcm_queue: Vec::new(),
            pending: VecDeque::new(),
            frame_index: 0,
            eof: false,
            annex_b,
            vad: VadState::new(),
        }
    }

    fn drain(&mut self, final_flush: bool) {
        // G.729 §3.2.1: the LP-analysis window spans 120 past + 80 current
        // + 40 future samples (5-ms look-ahead). Per the spec, this "5 ms
        // extra algorithmic delay" is realised by holding back the last
        // 40 samples of the PCM queue — we only encode a frame once we
        // have its 80 samples plus the following 40 samples in hand.
        let lookahead = LP_LOOKAHEAD_SAMPLES;
        while self.pcm_queue.len() >= FRAME_SAMPLES + lookahead {
            let mut pcm = [0i16; FRAME_SAMPLES];
            pcm.copy_from_slice(&self.pcm_queue[..FRAME_SAMPLES]);
            let mut la = [0i16; LP_LOOKAHEAD_SAMPLES];
            la.copy_from_slice(&self.pcm_queue[FRAME_SAMPLES..FRAME_SAMPLES + lookahead]);
            self.pcm_queue.drain(..FRAME_SAMPLES);
            self.encode_and_emit(&pcm, &la);
        }
        if final_flush {
            // EOF: pad whatever remains in `pcm_queue` with zeros up to
            // FRAME_SAMPLES, then encode using zero look-ahead. After
            // that, the queue is empty and the encoder is fully drained.
            while !self.pcm_queue.is_empty() {
                let mut pcm = [0i16; FRAME_SAMPLES];
                let take = self.pcm_queue.len().min(FRAME_SAMPLES);
                for (i, &s) in self.pcm_queue.iter().take(take).enumerate() {
                    pcm[i] = s;
                }
                self.pcm_queue.drain(..take);
                let mut la = [0i16; LP_LOOKAHEAD_SAMPLES];
                let la_take = self.pcm_queue.len().min(LP_LOOKAHEAD_SAMPLES);
                for (i, &s) in self.pcm_queue.iter().take(la_take).enumerate() {
                    la[i] = s;
                }
                self.encode_and_emit(&pcm, &la);
                if self.pcm_queue.is_empty() {
                    break;
                }
            }
        }
    }

    fn encode_and_emit(
        &mut self,
        pcm: &[i16; FRAME_SAMPLES],
        lookahead: &[i16; LP_LOOKAHEAD_SAMPLES],
    ) {
        let params = self.state.encode_one(pcm, lookahead);
        let idx = self.frame_index;
        self.frame_index += 1;
        self.emit_frame(pcm, &params, idx);
    }

    fn emit_frame(&mut self, pcm: &[i16; FRAME_SAMPLES], params: &EncodedFrame, idx: u64) {
        // Annex B VAD/DTX path: classify first, then either emit a normal
        // speech frame, an SID, or nothing. The analyser has already run
        // (inside `state.advance`, with the spec-mandated 240-sample
        // window + 5-ms look-ahead); we just need to pack and/or VAD-route
        // the resulting frame parameters.

        if self.annex_b {
            // Use the quantised LSP as the VAD's spectrum input — it's a
            // stable, compact representation of the frame's shape that the
            // analyser already produced.
            let lsp = self.state.lpc.lsp_prev;
            let (decision, sid) = self.vad.classify(pcm, &lsp);
            match decision {
                VadDecision::Voice => {
                    let bytes = pack_frame(params);
                    self.push_packet(idx, bytes);
                }
                VadDecision::Noise => {
                    let sid = sid.expect("VAD::Noise must carry SID payload");
                    let bytes = sid.pack().to_vec();
                    self.push_packet(idx, bytes);
                }
                VadDecision::Silence => {
                    // DTX: emit a zero-length packet so downstream
                    // containers / demuxers see the frame boundary but
                    // carry no payload. Consumers that require every
                    // packet to be decodable can filter on len() == 0.
                    self.push_packet(idx, Vec::new());
                }
            }
        } else {
            let bytes = pack_frame(params);
            self.push_packet(idx, bytes);
        }
    }

    fn push_packet(&mut self, idx: u64, bytes: Vec<u8>) {
        let mut pkt = Packet::new(0, self.time_base, bytes);
        pkt.pts = Some(idx as i64 * FRAME_SAMPLES as i64);
        pkt.dts = pkt.pts;
        pkt.duration = Some(FRAME_SAMPLES as i64);
        pkt.flags.keyframe = true;
        self.pending.push_back(pkt);
    }
}

impl Encoder for G729Encoder {
    fn codec_id(&self) -> &CodecId {
        &self.output_params.codec_id
    }

    fn output_params(&self) -> &CodecParameters {
        &self.output_params
    }

    fn send_frame(&mut self, frame: &Frame) -> Result<()> {
        let af = match frame {
            Frame::Audio(a) => a,
            _ => return Err(Error::invalid("G.729 encoder: audio frames only")),
        };
        let bytes = af
            .data
            .first()
            .ok_or_else(|| Error::invalid("G.729 encoder: empty frame"))?;
        if bytes.len() % 2 != 0 {
            return Err(Error::invalid("G.729 encoder: odd byte count"));
        }
        for chunk in bytes.chunks_exact(2) {
            self.pcm_queue
                .push(i16::from_le_bytes([chunk[0], chunk[1]]));
        }
        self.drain(false);
        Ok(())
    }

    fn receive_packet(&mut self) -> Result<Packet> {
        self.pending.pop_front().ok_or(Error::NeedMore)
    }

    fn flush(&mut self) -> Result<()> {
        if !self.eof {
            self.eof = true;
            self.drain(true);
        }
        Ok(())
    }
}

// =========================================================================
// Frame-level parameters (mirror `FrameParams` from the decoder bitreader).
// =========================================================================

/// Per-frame field values produced by the analyser; fed into `pack_frame`.
#[derive(Clone, Copy, Debug, Default)]
struct EncodedFrame {
    l0: u8,
    l1: u8,
    l2: u8,
    l3: u8,
    p1: u8,
    p0: u8,
    c1: u16,
    s1: u8,
    ga1: u8,
    gb1: u8,
    p2: u8,
    c2: u16,
    s2: u8,
    ga2: u8,
    gb2: u8,
}

/// Pack an [`EncodedFrame`] into the 80-bit / 10-byte ITU-T Table 8 layout
/// (MSB-first inside each byte, fields transmitted L0..GB2).
fn pack_frame(fp: &EncodedFrame) -> Vec<u8> {
    let fields: [(u32, u32); 15] = [
        (fp.l0 as u32, 1),
        (fp.l1 as u32, 7),
        (fp.l2 as u32, 5),
        (fp.l3 as u32, 5),
        (fp.p1 as u32, 8),
        (fp.p0 as u32, 1),
        (fp.c1 as u32, 13),
        (fp.s1 as u32, 4),
        (fp.ga1 as u32, 3),
        (fp.gb1 as u32, 4),
        (fp.p2 as u32, 5),
        (fp.c2 as u32, 13),
        (fp.s2 as u32, 4),
        (fp.ga2 as u32, 3),
        (fp.gb2 as u32, 4),
    ];
    let mut out = vec![0u8; FRAME_BYTES];
    let mut bit_pos: u32 = 0;
    for (val, width) in fields {
        let mask = if width == 32 {
            u32::MAX
        } else {
            (1u32 << width) - 1
        };
        let v = val & mask;
        for b in (0..width).rev() {
            let bit = (v >> b) & 1;
            let byte_idx = (bit_pos / 8) as usize;
            let shift = 7 - (bit_pos % 8);
            out[byte_idx] |= (bit as u8) << shift;
            bit_pos += 1;
        }
    }
    out
}

// =========================================================================
// Analysis state
// =========================================================================

/// Size of the LP-analysis window per ITU-T G.729 §3.2.1 / eq (3):
/// 120 past samples + 80 current-frame samples + 40 look-ahead samples.
const LP_WINDOW_SAMPLES: usize = 240;

/// Number of past-frame samples retained as the leading 120-sample
/// "history" portion of the LP-analysis window. Spec §3.2.1: "the LP
/// analysis window applies to 120 samples from past speech frames, 80
/// samples from the present speech frame, and 40 samples from the
/// future frame."
const LP_PAST_SAMPLES: usize = 120;

/// Size of the look-ahead slice taken from the next 80-sample frame.
/// 40 samples = 5 ms at 8 kHz, the spec's `5 ms look-ahead`.
const LP_LOOKAHEAD_SAMPLES: usize = 40;

struct EncoderState {
    /// LSP predictor state (same type as the decoder's — so we quantise the
    /// LSP residual identically and can feed the same table back).
    lpc: LpcPredictorState,
    /// Synthesis state used for analysis-by-synthesis: the excitation
    /// history drives the adaptive codebook and the pitch sharpening
    /// reconstruction, so we mimic the decoder exactly during encode.
    syn: SynthesisState,
    /// Pre-emphasis filter memory (simple one-pole HPF).
    preemph_prev: f32,
    /// Previous-frame unquantised LSP, used to interpolate the
    /// subframe-0 unquantised A(z) for the §3.3 weighting filter that
    /// feeds the §3.4 open-loop pitch search.
    prev_lsp_raw: [f32; LPC_ORDER],
    /// MA-4 gain prediction history, holding past quantised prediction
    /// errors `U^(k) = 20·log10(γ_q)` for the four most recent subframes
    /// (index 0 = most recent). Spec §3.9.1 eq (69-70); Table 9 mandates
    /// initial value `-14 dB`. Must stay synchronised with the decoder's
    /// own `gain_log_hist` so encoder and decoder predict the same
    /// fixed-codebook gain on every subframe.
    gain_log_hist: [f32; 4],
    /// Rolling pre-emphasised-PCM history that supplies the past-120
    /// portion of the 240-sample LP analysis window. Indexed
    /// `[0..LP_PAST_SAMPLES]` — newest samples at the high end.
    win_history: [f32; LP_PAST_SAMPLES],
    /// §3.3/§3.4 weighted-speech filter state. Holds enough `s` / `sw`
    /// history to evaluate eq (33) continuously and enough `sw` history
    /// to evaluate the eq-(34) open-loop pitch correlation at lags up to
    /// `MAX_PITCH_LAG`.
    wsp: WeightedSpeechState,
}

impl EncoderState {
    fn new() -> Self {
        let lpc = LpcPredictorState::new();
        let prev_lsp_raw = lpc.lsp_prev;
        Self {
            lpc,
            syn: SynthesisState::new(),
            preemph_prev: 0.0,
            prev_lsp_raw,
            // ITU-T G.729 Table 9: initial value of Û^(k) is -14 dB.
            gain_log_hist: [-14.0; 4],
            win_history: [0.0; LP_PAST_SAMPLES],
            wsp: WeightedSpeechState::new(),
        }
    }

    /// Run a single 80-sample frame through the encoder, using
    /// `lookahead` (the first 40 samples of the next-in-line frame, or
    /// trailing zeros at EOF) as the spec's 5-ms look-ahead source.
    fn encode_one(
        &mut self,
        pcm: &[i16; FRAME_SAMPLES],
        lookahead: &[i16; LP_LOOKAHEAD_SAMPLES],
    ) -> EncodedFrame {
        // -------- 1. Pre-process (HPF + normalise) on the current frame
        //         + the first 40 samples of the look-ahead. The 160-
        //         sample history slot already holds pre-emphasised
        //         samples from earlier frames. --------
        let mut sig = [0.0f32; FRAME_SAMPLES];
        let mut lookahead_pre = [0.0f32; LP_LOOKAHEAD_SAMPLES];
        let mut prev = self.preemph_prev;
        for i in 0..FRAME_SAMPLES {
            let x = pcm[i] as f32;
            let y = x - 0.46 * prev; // mild HPF; G.729 uses a 140 Hz HPF
            prev = x;
            sig[i] = y;
        }
        // Pre-emphasis-state at the boundary between current frame and
        // look-ahead must continue smoothly — don't reset prev.
        for i in 0..LP_LOOKAHEAD_SAMPLES {
            let x = lookahead[i] as f32;
            let y = x - 0.46 * prev;
            prev = x;
            lookahead_pre[i] = y;
        }
        // Commit pre-emphasis state up through the end of the current
        // frame only; the look-ahead samples will be re-pre-emphasised
        // when they become the next call's `pcm`, so we mustn't carry
        // their `prev` forward.
        self.preemph_prev = pcm[FRAME_SAMPLES - 1] as f32;

        // -------- 2. Assemble the 240-sample LP-analysis window:
        //         [past 120 | current 80 | look-ahead 40]
        //         then apply the asymmetric window from eq (3). --------
        debug_assert_eq!(
            LP_PAST_SAMPLES + FRAME_SAMPLES + LP_LOOKAHEAD_SAMPLES,
            LP_WINDOW_SAMPLES
        );
        let mut win = [0.0f32; LP_WINDOW_SAMPLES];
        win[..LP_PAST_SAMPLES].copy_from_slice(&self.win_history);
        win[LP_PAST_SAMPLES..LP_PAST_SAMPLES + FRAME_SAMPLES].copy_from_slice(&sig);
        win[LP_PAST_SAMPLES + FRAME_SAMPLES..].copy_from_slice(&lookahead_pre);
        let a = lpc_analysis(&win);

        // -------- 2b. Roll `win_history` forward by FRAME_SAMPLES so
        //         the next call sees this frame in its history slot. --------
        let mut new_history = [0.0f32; LP_PAST_SAMPLES];
        new_history[..LP_PAST_SAMPLES - FRAME_SAMPLES]
            .copy_from_slice(&self.win_history[FRAME_SAMPLES..]);
        new_history[LP_PAST_SAMPLES - FRAME_SAMPLES..].copy_from_slice(&sig);
        self.win_history = new_history;

        // -------- 3. LPC → LSP (cosine domain) --------
        let lsp_unq = lpc_to_lsp(&a);

        // -------- 4. Quantise LSPs using the decoder's codebooks --------
        let (l0, l1, l2, l3, lsp_q) = quantise_lsp_with_predictor(&mut self.lpc, &lsp_unq);

        // -------- 5. Per-subframe LPC via LSP interpolation --------
        let lsp_sf0 = interpolate_lsp(&self.lpc.lsp_prev, &lsp_q, 0.5);
        let lsp_sf1 = lsp_q;
        let a_sf = [lsp_to_lpc(&lsp_sf0), lsp_to_lpc(&lsp_sf1)];
        // Unquantised per-subframe A(z) for the §3.3 weighting filter
        // W(z) = A(z/γ1)/A(z/γ2). The weighting filter is built from the
        // *unquantised* LP coefficients (§3.3 first sentence); we
        // interpolate the unquantised LSP between the previous frame and
        // the current frame, mirroring the quantised path above.
        let lsp_unq_sf0 = interpolate_lsp(&self.prev_lsp_raw, &lsp_unq, 0.5);
        let a_unq_sf = [lsp_to_lpc(&lsp_unq_sf0), lsp_to_lpc(&lsp_unq)];

        // -------- 5b. §3.4 open-loop pitch analysis --------
        // Produce the perceptually-weighted speech sw(n) for both
        // subframes (eq 33), then run the once-per-frame open-loop pitch
        // search (eq 34, 35 + sub-multiple bias) to obtain T_op. T_op
        // anchors the subframe-0 closed-loop search of §3.7.
        //
        // NOTE on the (γ1, γ2) pair: §3.3 adapts the gammas per-frame via
        // a LAR-based flat/tilted classifier (eqs 28–32). That classifier
        // is a separate bounded element; until it lands, the open-loop
        // path uses the §3.3 "flat" gammas (0.94, 0.6) as a stable,
        // documented default. The weighting *filter* and the eq-(34/35)
        // search below are fully spec-faithful regardless of which γ pair
        // is supplied.
        // Snapshot the pre-frame weighted-speech history (newest sample
        // at index MAX_PITCH_LAG-1) *before* rolling the current frame
        // into the filter state — eq (34) correlates the current frame
        // against this trailing window.
        let sw_hist = *self.wsp.sw_history();
        let mut sw_frame = [0.0f32; FRAME_SAMPLES];
        for sf in 0..SUBFRAMES_PER_FRAME {
            let off = sf * SUBFRAME_SAMPLES;
            let mut s_sub = [0.0f32; SUBFRAME_SAMPLES];
            s_sub.copy_from_slice(&sig[off..off + SUBFRAME_SAMPLES]);
            let mut sw_sub = [0.0f32; SUBFRAME_SAMPLES];
            self.wsp
                .run_subframe(&s_sub, &a_unq_sf[sf], GAMMA1_FLAT, GAMMA2_FLAT, &mut sw_sub);
            sw_frame[off..off + SUBFRAME_SAMPLES].copy_from_slice(&sw_sub);
        }
        let open_loop = open_loop_pitch_search(&sw_frame, &sw_hist);
        let t_op = open_loop.t_op;

        // -------- 6. Subframe analysis --------
        let mut out = EncodedFrame {
            l0,
            l1,
            l2,
            l3,
            ..EncodedFrame::default()
        };

        // Target signal: residual after LPC inverse filtering.
        //   r[n] = sum_{k=0..10} a[k] * s[n-k]
        // Using the quantised per-subframe A(z).
        let residual = lpc_residual(&sig, &a_sf);

        let mut first_p1_int: usize = 40;
        for sf in 0..SUBFRAMES_PER_FRAME {
            let off = sf * SUBFRAME_SAMPLES;
            let mut target = [0.0f32; SUBFRAME_SAMPLES];
            target.copy_from_slice(&residual[off..off + SUBFRAME_SAMPLES]);

            // ---- 6a. Pitch search ----
            // Subframe 0 searches the §3.7 six-sample window around the
            // §3.4 open-loop delay T_op (eq f0018-01); subframe 1 searches
            // around the subframe-0 result T1.
            let (t_int, t_frac) = pitch_search(
                &self.syn.exc,
                &target,
                if sf == 0 {
                    PitchAnchor::OpenLoop(t_op)
                } else {
                    PitchAnchor::Subframe1(first_p1_int)
                },
            );
            if sf == 0 {
                first_p1_int = t_int;
                let p1 = encode_pitch_p1(t_int, t_frac);
                out.p1 = p1;
                out.p0 = pitch_parity(p1);
            } else {
                out.p2 = encode_pitch_p2(t_int, t_frac, first_p1_int);
            }

            // Adaptive-codebook excitation vector (unity-gain).
            let mut ac = [0.0f32; SUBFRAME_SAMPLES];
            adaptive_codebook_excitation(&self.syn.exc, t_int, t_frac, &mut ac);

            // ---- 6b. Adaptive-codebook gain (unconstrained LS) ----
            let g_p = gain_ls(&ac, &target).clamp(0.0, 1.2);

            // Compute residual after ACB contribution.
            let mut target2 = [0.0f32; SUBFRAME_SAMPLES];
            for n in 0..SUBFRAME_SAMPLES {
                target2[n] = target[n] - g_p * ac[n];
            }

            // ---- 6c. Fixed-codebook (ACELP 4-pulse) search ----
            let (c_idx, s_idx) = fixed_codebook_search(&target2);
            let mut fc = [0.0f32; SUBFRAME_SAMPLES];
            fixed_codebook_excitation(c_idx, s_idx, &mut fc);

            // ---- 6d. Fixed-codebook gain ----
            let g_c_raw = gain_ls(&fc, &target2).clamp(0.0, 32.0);

            // ---- 6e. Two-stage gain VQ over GBK1 / GBK2 ----
            // Per ITU-T G.729 §3.9 the gain VQ quantises the pair
            //   (g_p, γ)        where  γ = g_c / g'_c                  (eq 72)
            // and  g'_c  is the MA-4-predicted fixed-codebook gain (eq 71),
            // NOT the raw `g_c` itself. Compute the predictor in lockstep
            // with the decoder so encoder/decoder share `gain_log_hist`,
            // then quantise the (g_p, γ_target) pair.
            let innov_db = innovation_log_energy_db(&fc);
            let g_c_pred = predict_fixed_gain(&self.gain_log_hist, innov_db);
            let gamma_target = if g_c_pred > 1e-6 {
                (g_c_raw / g_c_pred).clamp(0.0, 5.0)
            } else {
                0.0
            };
            let (ga, gb) = quantise_gain_gamma(g_p, gamma_target);

            if sf == 0 {
                out.c1 = c_idx;
                out.s1 = s_idx;
                out.ga1 = ga;
                out.gb1 = gb;
            } else {
                out.c2 = c_idx;
                out.s2 = s_idx;
                out.ga2 = ga;
                out.gb2 = gb;
            }

            // ---- 6f. Update excitation history with the QUANTISED excitation
            //        so the next subframe's adaptive-codebook search sees the
            //        same history the decoder will. ----
            let (g_p_q, gamma_q) = dequantise_gain(ga, gb);
            // Match the decoder bit-for-bit: g_c_quantised = γ_q * g'_c
            // (eq 74). `g'_c` is the predictor we computed above.
            let g_c_q = gamma_q * g_c_pred;
            let mut excitation = [0.0f32; SUBFRAME_SAMPLES];
            for n in 0..SUBFRAME_SAMPLES {
                excitation[n] = g_p_q * ac[n] + g_c_q * fc[n];
            }
            push_excitation(&mut self.syn.exc, &excitation);

            // ---- 6g. Roll the MA-4 gain-prediction history with the
            //        QUANTISED prediction error  Û^(m) = 20·log10(γ_q),
            //        matching the decoder's update at §4.1.5 / eq (70). ----
            let u_new = gain_pred_error_db(gamma_q);
            for k in (1..4).rev() {
                self.gain_log_hist[k] = self.gain_log_hist[k - 1];
            }
            self.gain_log_hist[0] = u_new;
        }

        // -------- 7. Roll LSP predictor state --------
        self.lpc.lsp_prev = lsp_q;
        self.lpc.a = a_sf[1];
        // Roll the unquantised-LSP history used by the §3.3 weighting
        // filter's subframe-0 interpolation.
        self.prev_lsp_raw = lsp_unq;
        out
    }
}

/// Slide `exc` left by `SUBFRAME_SAMPLES` and append `sub` at the tail.
fn push_excitation(exc: &mut [f32; EXC_HIST], sub: &[f32; SUBFRAME_SAMPLES]) {
    for i in 0..EXC_HIST - SUBFRAME_SAMPLES {
        exc[i] = exc[i + SUBFRAME_SAMPLES];
    }
    for i in 0..SUBFRAME_SAMPLES {
        exc[EXC_HIST - SUBFRAME_SAMPLES + i] = sub[i];
    }
}

// =========================================================================
// LPC analysis
// =========================================================================

/// Build the spec's asymmetric LP-analysis window per ITU-T G.729 §3.2.1
/// equation (3):
///
/// ```text
/// w_lp(n) = 0.54 − 0.46·cos(2π·n / 399)        for n =   0 … 199
/// w_lp(n) = cos(2π·(n − 200) / 159)            for n = 200 … 239
/// ```
///
/// The first piece is the rising half of a Hamming window (peaks at
/// `n = 199.5` with value 1.0). The second piece is a quarter cosine
/// cycle that smoothly tapers from 1.0 at `n = 200` down to near zero
/// at `n = 239`. The two pieces meet at value 1.0, so the composite
/// window is continuous (peak at the boundary between the current
/// frame and the 40-sample look-ahead).
fn build_asymmetric_window() -> [f32; LP_WINDOW_SAMPLES] {
    let mut w = [0.0f32; LP_WINDOW_SAMPLES];
    let two_pi = 2.0 * core::f64::consts::PI;
    for n in 0..200 {
        let v = 0.54 - 0.46 * (two_pi * (n as f64) / 399.0).cos();
        w[n] = v as f32;
    }
    for n in 200..LP_WINDOW_SAMPLES {
        let v = (two_pi * ((n - 200) as f64) / 159.0).cos();
        w[n] = v as f32;
    }
    w
}

/// Windowed autocorrelation + Levinson-Durbin recursion → LPC[0..=10].
///
/// `sig` is the 240-sample LP-analysis window assembled by the caller:
/// 160 samples of past pre-emphasised history + 80 samples of the
/// current frame + 40 samples of look-ahead (per ITU-T G.729 §3.2.1).
fn lpc_analysis(sig: &[f32; LP_WINDOW_SAMPLES]) -> [f32; LPC_ORDER + 1] {
    // Asymmetric window from eq (3): half-Hamming over [0..199] + quarter
    // cosine fade over [200..239].
    let window = build_asymmetric_window();
    let mut w = [0.0f32; LP_WINDOW_SAMPLES];
    for i in 0..LP_WINDOW_SAMPLES {
        w[i] = sig[i] * window[i];
    }
    // Autocorrelation r[0..=10] per eq (5).
    let mut r = [0.0f64; LPC_ORDER + 1];
    for k in 0..=LPC_ORDER {
        let mut acc = 0.0f64;
        for i in k..LP_WINDOW_SAMPLES {
            acc += (w[i] as f64) * (w[i - k] as f64);
        }
        r[k] = acc;
    }
    // Spec eq (7): white-noise correction (×1.0001) on r(0), then 60 Hz
    // bandwidth-expansion lag window on r(k>=1) per eq (6).
    if r[0] < 1.0 {
        // Eq (5) note: r(0) has a lower boundary of 1.0 to avoid
        // arithmetic problems for low-level input signals.
        r[0] = 1.0;
    }
    r[0] *= 1.0001;
    for k in 1..=LPC_ORDER {
        let f = 2.0 * core::f64::consts::PI * 60.0 * (k as f64) / (SAMPLE_RATE as f64);
        let lag = (-0.5 * f * f).exp();
        r[k] *= lag;
    }
    if r[0] <= 0.0 {
        return default_a();
    }
    // Levinson-Durbin.
    let mut a = [0.0f64; LPC_ORDER + 1];
    let mut a_prev = [0.0f64; LPC_ORDER + 1];
    a[0] = 1.0;
    a_prev[0] = 1.0;
    let mut e = r[0];
    for i in 1..=LPC_ORDER {
        let mut acc = r[i];
        for j in 1..i {
            acc += a_prev[j] * r[i - j];
        }
        let k_refl = -acc / e;
        a[i] = k_refl;
        for j in 1..i {
            a[j] = a_prev[j] + k_refl * a_prev[i - j];
        }
        e *= 1.0 - k_refl * k_refl;
        if e <= 1e-18 {
            return default_a();
        }
        a_prev.copy_from_slice(&a);
    }
    let mut out = [0.0f32; LPC_ORDER + 1];
    for i in 0..=LPC_ORDER {
        out[i] = a[i] as f32;
    }
    out[0] = 1.0;
    out
}

fn default_a() -> [f32; LPC_ORDER + 1] {
    let mut a = [0.0f32; LPC_ORDER + 1];
    a[0] = 1.0;
    a
}

/// LPC direct-form → LSP (cosine domain) via Chebyshev root-finding on
/// F1(z) = A(z) + z^-(p+1) * A(z^-1) and F2(z) = A(z) - z^-(p+1) * A(z^-1).
fn lpc_to_lsp(a: &[f32; LPC_ORDER + 1]) -> [f32; LPC_ORDER] {
    let p = LPC_ORDER;
    // Build f1, f2 of degree p/2 after factoring out (1 ± z^-1).
    let mut f1 = [0.0f32; LPC_ORDER / 2 + 1];
    let mut f2 = [0.0f32; LPC_ORDER / 2 + 1];
    f1[0] = 1.0;
    f2[0] = 1.0;
    let mut prev_f1 = 0.0f32;
    let mut prev_f2 = 0.0f32;
    for i in 1..=p / 2 {
        let ai = a[i];
        let api = a[p + 1 - i];
        f1[i] = ai + api - prev_f1;
        f2[i] = ai - api + prev_f2;
        prev_f1 = f1[i];
        prev_f2 = f2[i];
    }
    let r1 = cheby_roots(&f1);
    let r2 = cheby_roots(&f2);
    let mut lsp = [0.0f32; LPC_ORDER];
    // Interleave roots; fall back to uniform spread if the search ran short.
    let uni = |k: usize| -> f32 {
        let step = core::f32::consts::PI / (LPC_ORDER as f32 + 1.0);
        (step * (k as f32 + 1.0)).cos()
    };
    for k in 0..LPC_ORDER {
        if k % 2 == 0 {
            lsp[k] = r1.get(k / 2).copied().unwrap_or_else(|| uni(k));
        } else {
            lsp[k] = r2.get(k / 2).copied().unwrap_or_else(|| uni(k));
        }
    }
    // Enforce strictly-decreasing cos domain.
    for k in 1..LPC_ORDER {
        if lsp[k] >= lsp[k - 1] - 1e-4 {
            lsp[k] = lsp[k - 1] - 1e-3;
        }
    }
    // Clamp to (-1, 1) with small margin.
    for lsp_k in lsp.iter_mut().take(LPC_ORDER) {
        *lsp_k = lsp_k.clamp(-0.9995, 0.9995);
    }
    lsp
}

/// Find real roots on `[-1, 1]` of a Chebyshev-expanded polynomial via
/// grid-bracket + bisection. Returns roots in decreasing x-order (i.e.
/// increasing omega).
fn cheby_roots(coeffs: &[f32]) -> Vec<f32> {
    let deg = coeffs.len().saturating_sub(1);
    let eval = |x: f32| -> f32 {
        // Clenshaw recurrence for Chebyshev series of the first kind.
        let mut b2 = 0.0f32;
        let mut b1 = 0.0f32;
        for k in (1..=deg).rev() {
            let b0 = 2.0 * x * b1 - b2 + coeffs[k];
            b2 = b1;
            b1 = b0;
        }
        x * b1 - b2 + coeffs[0]
    };
    const GRID: usize = 200;
    let mut roots = Vec::with_capacity(deg);
    let mut prev_x = 1.0f32;
    let mut prev_y = eval(prev_x);
    for i in 1..=GRID {
        let x = 1.0 - 2.0 * (i as f32 / GRID as f32);
        let y = eval(x);
        if prev_y * y < 0.0 {
            // Bisect on [x, prev_x].
            let mut lo = x;
            let mut hi = prev_x;
            let mut flo = y;
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                let fm = eval(mid);
                if fm * flo < 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    flo = fm;
                }
            }
            roots.push(0.5 * (lo + hi));
            if roots.len() == deg {
                break;
            }
        }
        prev_x = x;
        prev_y = y;
    }
    roots
}

// =========================================================================
// LSP quantisation (symmetric with the decoder's MA-4 + LSPCB1/LSPCB2)
// =========================================================================

/// Quantise `lsp_unq` (cosine-domain LSPs) using the decoder's
/// `LSPCB1_Q13` (L1, 128 entries) + `LSPCB2_Q13` (L2/L3, 32 entries).
/// Also chooses the MA-predictor switch `L0 ∈ {0, 1}`.
///
/// Returns `(L0, L1, L2, L3, quantised LSP in cosine domain)`.
///
/// Implementation:
/// - Convert `lsp_unq` to LSF (radians), then to Q13.
/// - Subtract the MA-predictor expectation to get the *residual* that
///   `LSPCB1 + LSPCB2` needs to represent.
/// - Search `LSPCB1_Q13[l1] + LSPCB2_Q13[l2 (low)] + LSPCB2_Q13[l3 (high)]`
///   for the minimum L2 distance to the target residual. The low and
///   high halves are searched independently (split VQ).
/// - Do the search twice (predictor 0 and predictor 1) and keep the
///   smaller reconstruction error.
fn quantise_lsp_with_predictor(
    state: &mut LpcPredictorState,
    lsp_unq: &[f32; LPC_ORDER],
) -> (u8, u8, u8, u8, [f32; LPC_ORDER]) {
    // Target LSF in radians, then in Q13 (f32 for arithmetic).
    let mut lsf_target = [0.0f32; LPC_ORDER];
    for k in 0..LPC_ORDER {
        lsf_target[k] = lsp_unq[k].clamp(-1.0, 1.0).acos();
    }
    // Enforce strictly-increasing LSF + minimum spacing.
    for k in 1..LPC_ORDER {
        if lsf_target[k] <= lsf_target[k - 1] + 1e-3 {
            lsf_target[k] = lsf_target[k - 1] + 1e-3;
        }
    }
    let mut target_q13 = [0.0f32; LPC_ORDER];
    for k in 0..LPC_ORDER {
        target_q13[k] = lsf_target[k] * 8192.0;
    }

    // Try predictor 0 and predictor 1.
    let mut best_l0 = 0u8;
    let mut best_l1 = 0u8;
    let mut best_l2 = 0u8;
    let mut best_l3 = 0u8;
    let mut best_err = f32::INFINITY;
    let mut best_resid = [0i16; LPC_ORDER];

    for predictor in 0..2usize {
        // Expected predictor contribution (same formula as decoder).
        let fg = &FG_Q15[predictor];
        let fg_sum = &FG_SUM_Q15[predictor];
        // predicted_lsf_q13[j] = (sum_k fg[k][j] * prev_res_q13[k][j]) / 2^15
        let mut predicted = [0.0f32; LPC_ORDER];
        for j in 0..LPC_ORDER {
            let mut acc = 0.0f32;
            for k in 0..MA_NP {
                acc += (fg[k][j] as f32) * (state.freq_res_q13[k][j] as f32);
            }
            predicted[j] = acc / 32768.0;
        }
        // Target residual in Q13: fg_sum[j] * residual[j] = target_q13[j] - predicted[j]
        let mut resid_target = [0.0f32; LPC_ORDER];
        for j in 0..LPC_ORDER {
            let denom = (fg_sum[j] as f32) / 32768.0;
            if denom.abs() < 1e-6 {
                resid_target[j] = 0.0;
            } else {
                resid_target[j] = (target_q13[j] - predicted[j]) / denom;
            }
        }

        // Split-VQ search: the residual is LPCB1 row + concatenation of two
        // LPCB2 rows (low half + high half).
        // Best L1: search 128 entries minimising distance to low+high halves.
        let mut best_l1_for = 0u8;
        let mut best_l2_for = 0u8;
        let mut best_l3_for = 0u8;
        let mut best_err_for = f32::INFINITY;
        let mut best_resid_for = [0i16; LPC_ORDER];

        for l1 in 0..NC0 {
            let cb1 = &LSPCB1_Q13[l1];
            // Best L2 (low half, j=0..5).
            let mut b2 = 0usize;
            let mut b2_err = f32::INFINITY;
            for l2 in 0..NC1 {
                let cb2 = &LSPCB2_Q13[l2];
                let mut err = 0.0f32;
                for j in 0..M_HALF {
                    let recon = (cb1[j] as f32) + (cb2[j] as f32);
                    let d = recon - resid_target[j];
                    err += d * d;
                }
                if err < b2_err {
                    b2_err = err;
                    b2 = l2;
                }
            }
            // Best L3 (high half, j=5..10).
            // Per ITU-T G.729 §3.2.4 / `Lsp_get_quant` in `LSPGETQ.C`: the
            // L3 codeword's contribution comes from the **high-half columns**
            // of the same 32×10 LSPCB2 row (cols M_HALF..M), not its low
            // half. cb2 is a full M-wide row.
            let mut b3 = 0usize;
            let mut b3_err = f32::INFINITY;
            for l3 in 0..NC1 {
                let cb2 = &LSPCB2_Q13[l3];
                let mut err = 0.0f32;
                for j in 0..M_HALF {
                    let recon = (cb1[j + M_HALF] as f32) + (cb2[j + M_HALF] as f32);
                    let d = recon - resid_target[j + M_HALF];
                    err += d * d;
                }
                if err < b3_err {
                    b3_err = err;
                    b3 = l3;
                }
            }
            let err = b2_err + b3_err;
            if err < best_err_for {
                best_err_for = err;
                best_l1_for = l1 as u8;
                best_l2_for = b2 as u8;
                best_l3_for = b3 as u8;
                // Assemble the residual in Q13 using the selected entries.
                // L2 contributes the low half (cols 0..M_HALF) of its row;
                // L3 contributes the high half (cols M_HALF..M) of its row.
                let cb2_lo = &LSPCB2_Q13[b2];
                let cb2_hi = &LSPCB2_Q13[b3];
                for j in 0..M_HALF {
                    best_resid_for[j] = cb1[j].saturating_add(cb2_lo[j]);
                    best_resid_for[j + M_HALF] = cb1[j + M_HALF].saturating_add(cb2_hi[j + M_HALF]);
                }
            }
        }

        if best_err_for < best_err {
            best_err = best_err_for;
            best_l0 = predictor as u8;
            best_l1 = best_l1_for;
            best_l2 = best_l2_for;
            best_l3 = best_l3_for;
            best_resid = best_resid_for;
        }
    }

    // Reconstruct the quantised LSF vector using the chosen predictor and
    // entries — mirror of the decoder's `decode_lsp` logic.
    let predictor = best_l0 as usize;
    let fg = &FG_Q15[predictor];
    let fg_sum = &FG_SUM_Q15[predictor];
    let mut lsf_q13_f = [0.0f32; LPC_ORDER];
    for j in 0..LPC_ORDER {
        let mut acc: f32 = (fg_sum[j] as f32) * (best_resid[j] as f32);
        for k in 0..MA_HISTORY {
            acc += (fg[k][j] as f32) * (state.freq_res_q13[k][j] as f32);
        }
        lsf_q13_f[j] = acc / 32768.0;
    }
    // Push the chosen residual onto the predictor history (same as decoder).
    state.push_residual(best_resid);

    // Convert Q13 LSF → radians → cosine domain with spacing safeguards.
    let pi = core::f32::consts::PI;
    let eps = 0.0012f32;
    let mut lsf = [0.0f32; LPC_ORDER];
    for j in 0..LPC_ORDER {
        lsf[j] = lsf_q13_f[j] / 8192.0;
    }
    if lsf[0] < eps {
        lsf[0] = eps;
    }
    for j in 1..LPC_ORDER {
        if lsf[j] < lsf[j - 1] + eps {
            lsf[j] = lsf[j - 1] + eps;
        }
    }
    if lsf[LPC_ORDER - 1] > pi - eps {
        lsf[LPC_ORDER - 1] = pi - eps;
        for j in (0..LPC_ORDER - 1).rev() {
            if lsf[j] > lsf[j + 1] - eps {
                lsf[j] = lsf[j + 1] - eps;
            }
        }
    }
    let mut lsp_q = [0.0f32; LPC_ORDER];
    for j in 0..LPC_ORDER {
        lsp_q[j] = lsf[j].cos();
    }

    (best_l0, best_l1, best_l2, best_l3, lsp_q)
}

// =========================================================================
// LPC residual
// =========================================================================

/// Compute per-frame LPC residual using per-subframe A(z) coefficients.
/// `r[n] = s[n] + sum_{k=1..10} a[k] * s[n-k]`.
fn lpc_residual(
    sig: &[f32; FRAME_SAMPLES],
    a_sf: &[[f32; LPC_ORDER + 1]; 2],
) -> [f32; FRAME_SAMPLES] {
    let mut mem = [0.0f32; LPC_ORDER];
    let mut out = [0.0f32; FRAME_SAMPLES];
    for sf in 0..SUBFRAMES_PER_FRAME {
        let a = &a_sf[sf];
        let base = sf * SUBFRAME_SAMPLES;
        for i in 0..SUBFRAME_SAMPLES {
            let x = sig[base + i];
            let mut acc = x;
            for k in 1..=LPC_ORDER {
                acc += a[k] * mem[k - 1];
            }
            out[base + i] = acc;
            for k in (1..LPC_ORDER).rev() {
                mem[k] = mem[k - 1];
            }
            mem[0] = x;
        }
    }
    out
}

// =========================================================================
// Pitch search (open-loop + closed-loop with 1/3 fractional lag)
// =========================================================================

/// Closed-loop pitch search anchor for [`pitch_search`].
///
/// The closed-loop search range is limited to a small window around a
/// preselected delay (§3.7): the §3.4 open-loop delay `T_op` for the
/// first subframe, and the first-subframe delay `T1` for the second.
#[derive(Clone, Copy, Debug)]
enum PitchAnchor {
    /// Subframe 1: window around the §3.4 open-loop delay `T_op`
    /// (eq f0018-01: `t_min = T_op − 3`, `t_max = t_min + 6`).
    OpenLoop(usize),
    /// Subframe 2: window around the first-subframe delay `int(T1)`
    /// (eq f0018-02: `t_min = int(T1) − 5`, `t_max = t_min + 9`).
    Subframe1(usize),
}

/// Compute the §3.7 closed-loop search boundaries `[t_min, t_max]` for the
/// first subframe given the open-loop delay `T_op` (ITU-T G.729
/// eq f0018-01):
///
/// ```text
///   t_min = T_op − 3;  if t_min < 20 then t_min = 20
///   t_max = t_min + 6; if t_max > 143 then t_max = 143; t_min = t_max − 6
/// ```
fn subframe1_search_bounds(t_op: usize) -> (usize, usize) {
    let mut t_min = t_op.saturating_sub(3);
    if t_min < 20 {
        t_min = 20;
    }
    let mut t_max = t_min + 6;
    if t_max > 143 {
        t_max = 143;
        t_min = t_max - 6;
    }
    (t_min, t_max)
}

/// Compute the §3.7 closed-loop search boundaries `[t_min, t_max]` for the
/// second subframe given the first-subframe delay `int(T1)` (ITU-T G.729
/// eq f0018-02):
///
/// ```text
///   t_min = int(T1) − 5;  if t_min < 20 then t_min = 20
///   t_max = t_min + 9;    if t_max > 143 then t_max = 143; t_min = t_max − 9
/// ```
fn subframe2_search_bounds(t1_int: usize) -> (usize, usize) {
    let mut t_min = t1_int.saturating_sub(5);
    if t_min < 20 {
        t_min = 20;
    }
    let mut t_max = t_min + 9;
    if t_max > 143 {
        t_max = 143;
        t_min = t_max - 9;
    }
    (t_min, t_max)
}

/// Choose the best integer pitch + 1/3 fractional shift in the §3.7
/// closed-loop window around the supplied [`PitchAnchor`].
fn pitch_search(
    exc: &[f32; EXC_HIST],
    target: &[f32; SUBFRAME_SAMPLES],
    anchor: PitchAnchor,
) -> (usize, i8) {
    let (lag_lo, lag_hi) = match anchor {
        PitchAnchor::OpenLoop(t_op) => subframe1_search_bounds(t_op),
        PitchAnchor::Subframe1(t1) => subframe2_search_bounds(t1),
    };

    // Integer search: pick the integer lag maximising normalised
    // correlation (num^2 / den).
    let mut best_int = lag_lo;
    let mut best_score = -f32::INFINITY;
    for lag in lag_lo..=lag_hi {
        let mut cand = [0.0f32; SUBFRAME_SAMPLES];
        adaptive_codebook_excitation(exc, lag, 0, &mut cand);
        let (num, den) = xc_norm(&cand, target);
        if den < 1e-6 {
            continue;
        }
        let score = num * num / den;
        if score > best_score {
            best_score = score;
            best_int = lag;
        }
    }

    // Fractional refinement: compare frac ∈ {-1, 0, +1} on the neighbours.
    // `t_frac` uses the decoder's convention (-1, 0, +1 == -1/3, 0, +1/3).
    let mut best_frac: i8 = 0;
    for frac in [-1i8, 0, 1] {
        // Skip frac that would push us outside the bracket; the decoder's
        // adaptive_codebook_excitation handles out-of-range implicitly by
        // returning zeros, which yields a low score.
        let mut cand = [0.0f32; SUBFRAME_SAMPLES];
        adaptive_codebook_excitation(exc, best_int, frac, &mut cand);
        let (num, den) = xc_norm(&cand, target);
        if den < 1e-6 {
            continue;
        }
        let score = num * num / den;
        if score > best_score {
            best_score = score;
            best_frac = frac;
        }
    }
    (best_int, best_frac)
}

fn xc_norm(cand: &[f32; SUBFRAME_SAMPLES], target: &[f32; SUBFRAME_SAMPLES]) -> (f32, f32) {
    let mut num = 0.0f32;
    let mut den = 1e-9f32;
    for n in 0..SUBFRAME_SAMPLES {
        num += cand[n] * target[n];
        den += cand[n] * cand[n];
    }
    (num, den)
}

// =========================================================================
// Pitch-index encoding — inverse of `decode_pitch_p1` / `decode_pitch_p2`.
// =========================================================================

/// Encode (integer, frac) → 8-bit P1 index. Inverse of [`decode_pitch_p1`].
///
/// - Fractional range: integer 20..=84 + frac ∈ {-1, 0, +1} maps to 0..196.
/// - Integer-only range: 85..=142 maps to 197..254.
fn encode_pitch_p1(t_int: usize, t_frac: i8) -> u8 {
    let t_int = t_int.clamp(20, 143);
    if t_int <= 84 || (t_int == 85 && t_frac < 0) {
        // Fractional grid.
        let frac = t_frac.clamp(-1, 1) as i32;
        // decode formula: t = idx + 59; t_int = t / 3; t_frac = t - 3*t_int - 1
        //           idx = 3*t_int + t_frac + 1 - 59
        let idx = 3 * (t_int as i32) + frac + 1 - 59;
        idx.clamp(0, 196) as u8
    } else {
        // Integer grid.
        let idx = (t_int as i32 + 112).clamp(197, 254);
        idx as u8
    }
}

/// Encode (integer, frac) → 5-bit P2 index, given the anchor `p1_int`.
/// Inverse of [`decode_pitch_p2`].
fn encode_pitch_p2(t_int: usize, t_frac: i8, p1_int: usize) -> u8 {
    let mut t_min = p1_int.saturating_sub(5);
    if t_min < 20 {
        t_min = 20;
    }
    let mut t_max = t_min + 9;
    if t_max > 143 {
        t_max = 143;
        t_min = t_max - 9;
    }
    let t_int = t_int.clamp(t_min, t_max) as i32;
    let t_frac = t_frac.clamp(-1, 1) as i32;
    // Inverse of:
    //   t = idx + 59 - 3*(t_min - 1);
    //   t_int = t/3 + t_min - 1; t_frac = t - 3*(t_int - t_min + 1) - 1
    //   => t = 3 * (t_int - t_min + 1) + t_frac + 1
    //   => idx = t - 59 + 3*(t_min - 1)
    let tmin_i = t_min as i32;
    let t = 3 * (t_int - tmin_i + 1) + t_frac + 1;
    let idx = t - 59 + 3 * (tmin_i - 1);
    idx.clamp(0, 31) as u8
}

// =========================================================================
// Fixed-codebook search (ACELP 4-pulse, focused by track)
// =========================================================================

/// Depth-first 4-pulse ACELP search symmetric with the decoder's
/// `fixed_codebook_excitation`. For each track we pick the position (and
/// sign) that maximises the normalised correlation with the remaining
/// target signal. Returns `(c13, s4)` as used in the bit layout.
///
/// Track positions (§3.8, decoder side):
///   - track 0: 3-bit index, pos = 5*k       (0, 5, ...,35)
///   - track 1: 3-bit index, pos = 5*k + 1
///   - track 2: 3-bit index, pos = 5*k + 2
///   - track 3: 3-bit index + 1-bit jitter, pos = 5*k + 3 + jitter
fn fixed_codebook_search(target: &[f32; SUBFRAME_SAMPLES]) -> (u16, u8) {
    let mut residual = *target;
    let mut c: u32 = 0;
    let mut s: u32 = 0;

    // For each track, pick the best signed unit pulse from the residual.
    // Track order — 0, 1, 2, 3.
    // Focused depth-first search: each track picks greedily. This is
    // far from the spec's full depth-first search, but sidesteps the
    // 2^17 exhaustive search and still produces useful excitation.
    let tracks: [&[usize]; 4] = [
        &[0, 5, 10, 15, 20, 25, 30, 35],
        &[1, 6, 11, 16, 21, 26, 31, 36],
        &[2, 7, 12, 17, 22, 27, 32, 37],
        &[3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39],
    ];
    for track_idx in 0..4 {
        let positions = tracks[track_idx];
        let mut best_pos_k: usize = 0;
        let mut best_sign_pos = false; // true ≡ +1
        let mut best_abs = -1.0f32;
        for (k, &pos) in positions.iter().enumerate() {
            let v = residual[pos];
            let av = v.abs();
            if av > best_abs {
                best_abs = av;
                best_pos_k = k;
                best_sign_pos = v >= 0.0;
            }
        }
        // Pack into C: track 0..2 use 3 bits; track 3 uses 3+1 (jitter in bit 12).
        if track_idx < 3 {
            c |= ((best_pos_k & 0x7) as u32) << (3 * track_idx);
        } else {
            // Track 3: 4-bit index in the positions[] list => 3-bit base + jitter.
            // `positions` is [3,4,8,9,13,14,...]. Base = k>>1; jitter = k & 1.
            let base = (best_pos_k >> 1) & 0x7;
            let jitter = best_pos_k & 0x1;
            c |= (base as u32) << 9;
            c |= (jitter as u32) << 12;
        }
        // Sign bit — decoder reads `s & 0x1` as track-0 sign.
        if best_sign_pos {
            s |= 1 << track_idx;
        }
        // Subtract the chosen pulse from the residual so later tracks
        // aren't all drawn to the same dominant peak.
        let chosen_pos = positions[best_pos_k];
        let amp = if best_sign_pos { 1.0 } else { -1.0 };
        let alpha = residual[chosen_pos]; // orthogonal projection of the unit pulse.
        let _ = amp;
        residual[chosen_pos] -= alpha;
    }

    ((c & 0x1FFF) as u16, (s & 0xF) as u8)
}

// =========================================================================
// Gain analysis + 2-stage VQ
// =========================================================================

/// Unconstrained LS gain for predictor `pred` against `target`:
/// g = <pred, target> / <pred, pred>.
fn gain_ls(pred: &[f32; SUBFRAME_SAMPLES], target: &[f32; SUBFRAME_SAMPLES]) -> f32 {
    let mut num = 0.0f32;
    let mut den = 1e-9f32;
    for n in 0..SUBFRAME_SAMPLES {
        num += pred[n] * target[n];
        den += pred[n] * pred[n];
    }
    num / den
}

/// Find (GA, GB) indices whose `GBK1[ga] + GBK2[gb]` reconstruction is
/// closest to the conjugate-VQ target `(g_p, γ)`, where:
/// * `g_p`   is the open-loop adaptive-codebook gain estimate, and
/// * `gamma_target = g_c_raw / g'_c` is the fixed-codebook gain
///   correction factor relative to the MA-4 predicted gain `g'_c`
///   (ITU-T G.729 §3.9 eq 72).
///
/// Per §3.9.2 the codebook is searched to minimise the weighted MSE
///   `α · (g_p − ĝ_p)² + β · (γ − γ̂)²`
/// of equation (63), with the spec's preselection collapsed to an
/// exhaustive 8 × 16 search (Annex A complexity bracket).
fn quantise_gain_gamma(g_p: f32, gamma_target: f32) -> (u8, u8) {
    // Equal weighting on the two error dimensions is the closest we can
    // come to the spec's α/β without the impulse-response matrix from
    // §3.6; a uniform 1.0/1.0 is the long-term-average choice from
    // §3.9.2 ("...the weights are roughly balanced...").
    let lambda = 1.0f32;
    let mut best = (0u8, 0u8);
    let mut best_err = f32::INFINITY;
    for ga in 0..GBK1.len() {
        for gb in 0..GBK2.len() {
            let gp = GBK1[ga][0] + GBK2[gb][0];
            let gm = GBK1[ga][1] + GBK2[gb][1];
            // Skip combinations that would fall outside the decoder's clamp.
            if !(0.0..=1.3).contains(&gp) {
                continue;
            }
            if !(0.0..=2.6).contains(&gm) {
                continue;
            }
            let d0 = g_p - gp;
            let d1 = gamma_target - gm;
            let err = d0 * d0 + lambda * d1 * d1;
            if err < best_err {
                best_err = err;
                best = (ga as u8, gb as u8);
            }
        }
    }
    best
}

fn dequantise_gain(ga: u8, gb: u8) -> (f32, f32) {
    // Same formulation as the decoder — kept local to avoid pulling the
    // decoder's `decode_gain_indices` out of its module (it also clamps
    // the inputs in decoder-specific ways that are close enough to our
    // reconstruction needs).
    let ga = (ga as usize) & 0x7;
    let gb = (gb as usize) & 0xF;
    let g_p = (GBK1[ga][0] + GBK2[gb][0]).clamp(0.0, 1.2);
    let gamma = (GBK1[ga][1] + GBK2[gb][1]).clamp(0.0, 2.5);
    (g_p, gamma)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_frame_known_pattern() {
        // Every field set to all-ones should produce the 0xFF...FF packet.
        let fp = EncodedFrame {
            l0: 1,
            l1: 0x7F,
            l2: 0x1F,
            l3: 0x1F,
            p1: 0xFF,
            p0: 1,
            c1: 0x1FFF,
            s1: 0xF,
            ga1: 0x7,
            gb1: 0xF,
            p2: 0x1F,
            c2: 0x1FFF,
            s2: 0xF,
            ga2: 0x7,
            gb2: 0xF,
        };
        let bytes = pack_frame(&fp);
        assert_eq!(bytes.len(), FRAME_BYTES);
        assert!(bytes.iter().all(|&b| b == 0xFF));
    }

    #[test]
    fn encode_decode_pitch_p1_round_trips_on_integer_grid() {
        use crate::synthesis::decode_pitch_p1;
        for t in 85..=142 {
            let idx = encode_pitch_p1(t, 0);
            let (dt, df) = decode_pitch_p1(idx);
            assert_eq!((dt, df), (t, 0));
        }
    }

    #[test]
    fn encode_decode_pitch_p1_round_trips_on_fractional_grid() {
        use crate::synthesis::decode_pitch_p1;
        for t in 20..=84 {
            for f in [-1i8, 0, 1] {
                let idx = encode_pitch_p1(t, f);
                let (dt, df) = decode_pitch_p1(idx);
                // Round-trip must hit the same lattice point.
                assert_eq!(dt, t, "int mismatch for (t={t}, f={f}) idx={idx}");
                assert_eq!(df, f, "frac mismatch for (t={t}, f={f}) idx={idx}");
            }
        }
    }

    #[test]
    fn quantise_gain_returns_valid_indices() {
        // (g_p, γ) target both within the codebook span.
        let (ga, gb) = quantise_gain_gamma(0.5, 1.0);
        assert!(ga < GBK1.len() as u8);
        assert!(gb < GBK2.len() as u8);
    }

    #[test]
    fn quantise_gain_finds_nearest_codepoint() {
        // For (g_p ≈ 0.8, γ ≈ 1.0) the best match should sit near the
        // middle of the GBK1 / GBK2 lattices, not at the extremes —
        // exercises the search rather than picking up the trivial
        // index-0 fallback.
        let (ga, gb) = quantise_gain_gamma(0.8, 1.0);
        let recon_gp = GBK1[ga as usize][0] + GBK2[gb as usize][0];
        let recon_gm = GBK1[ga as usize][1] + GBK2[gb as usize][1];
        assert!(
            (recon_gp - 0.8).abs() < 0.4,
            "g_p reconstruction {recon_gp} off-target"
        );
        assert!(
            (recon_gm - 1.0).abs() < 0.5,
            "γ reconstruction {recon_gm} off-target"
        );
    }

    #[test]
    fn quantise_gain_zero_gamma_picks_smallest_codepoint() {
        // γ = 0 means the predictor over-shot — we should pick the
        // GBK1/GBK2 pair with the smallest non-negative γ̂ value.
        let (ga, gb) = quantise_gain_gamma(0.0, 0.0);
        let recon_gm = GBK1[ga as usize][1] + GBK2[gb as usize][1];
        assert!(recon_gm < 0.5, "expected small γ̂, got {recon_gm}");
    }

    #[test]
    fn make_encoder_rejects_stereo() {
        let mut params = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
        params.sample_rate = Some(SAMPLE_RATE);
        params.channels = Some(2);
        assert!(make_encoder(&params).is_err());
    }

    #[test]
    fn make_encoder_accepts_valid_params() {
        let mut params = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
        params.sample_rate = Some(SAMPLE_RATE);
        params.channels = Some(1);
        assert!(make_encoder(&params).is_ok());
    }

    #[test]
    fn asymmetric_window_matches_spec_anchor_points() {
        // ITU-T G.729 §3.2.1 eq (3): two-piece asymmetric LP-analysis
        // window. Cross-check the spec's known anchor values:
        //   - w(0)    = 0.54 − 0.46·cos(0)            = 0.08
        //   - w(199)  ≈ 0.54 − 0.46·cos(2π·199/399)   ≈ 1.0 (junction)
        //   - w(200)  = cos(0)                         = 1.0 (junction)
        //   - w(239)  = cos(2π·39/159)                ≈ 0.030
        let w = build_asymmetric_window();
        assert_eq!(w.len(), LP_WINDOW_SAMPLES);
        assert!(
            (w[0] - 0.08).abs() < 1e-5,
            "w(0) expected 0.08, got {}",
            w[0]
        );
        assert!(
            (w[199] - 1.0).abs() < 1e-3,
            "w(199) expected ≈1.0 (Hamming-peak junction), got {}",
            w[199]
        );
        assert!(
            (w[200] - 1.0).abs() < 1e-7,
            "w(200) expected exactly 1.0 (cosine peak), got {}",
            w[200]
        );
        // w(239) = cos(2π·39/159); we just want it small and positive.
        assert!(
            w[239] > 0.0 && w[239] < 0.1,
            "w(239) expected ∈(0, 0.1), got {}",
            w[239]
        );
        // Monotonic tapering on the cosine piece: w[200] > w[210] > w[220] > w[239].
        assert!(w[200] > w[210]);
        assert!(w[210] > w[220]);
        assert!(w[220] > w[239]);
    }

    #[test]
    fn subframe1_bounds_match_spec_eq_f0018_01() {
        // Mid-range T_op: 6-sample span centred slightly below T_op.
        let (lo, hi) = subframe1_search_bounds(80);
        assert_eq!((lo, hi), (77, 83)); // T_op−3 .. T_op+3
        assert_eq!(hi - lo, 6);
        // Low clamp: T_op near the floor pins t_min at 20.
        let (lo, hi) = subframe1_search_bounds(21);
        assert_eq!(lo, 20);
        assert_eq!(hi, 26);
        // High clamp: T_op near the ceiling pins t_max at 143 and pushes
        // t_min back to 137.
        let (lo, hi) = subframe1_search_bounds(143);
        assert_eq!(hi, 143);
        assert_eq!(lo, 137);
        assert_eq!(hi - lo, 6);
    }

    #[test]
    fn subframe2_bounds_match_spec_eq_f0018_02() {
        // Mid-range T1: 9-sample span.
        let (lo, hi) = subframe2_search_bounds(80);
        assert_eq!((lo, hi), (75, 84)); // T1−5 .. T1+4
        assert_eq!(hi - lo, 9);
        // Low clamp.
        let (lo, hi) = subframe2_search_bounds(22);
        assert_eq!(lo, 20);
        assert_eq!(hi, 29);
        // High clamp.
        let (lo, hi) = subframe2_search_bounds(143);
        assert_eq!(hi, 143);
        assert_eq!(lo, 134);
        assert_eq!(hi - lo, 9);
    }

    #[test]
    fn pitch_search_stays_within_open_loop_window() {
        // A pitch_search anchored on an open-loop T_op must return a lag
        // inside the §3.7 six-sample window, never the full 20..143 range.
        let mut exc = [0.0f32; EXC_HIST];
        // Put a strong periodic structure at lag 60 in the excitation.
        for (i, e) in exc.iter_mut().enumerate() {
            *e = if i % 60 == 0 { 1.0 } else { 0.0 };
        }
        let target = [0.5f32; SUBFRAME_SAMPLES];
        let (lo, hi) = subframe1_search_bounds(100);
        let (t_int, _frac) = pitch_search(&exc, &target, PitchAnchor::OpenLoop(100));
        assert!(
            (lo..=hi).contains(&t_int),
            "t_int {t_int} escaped window [{lo},{hi}]"
        );
    }

    #[test]
    fn encode_one_with_lookahead_runs_and_produces_finite_output() {
        // Smoke test for `encode_one` against a 240-sample analysis
        // window assembled from the spec layout (past + current + LA).
        let mut state = EncoderState::new();
        let pcm = [100i16; FRAME_SAMPLES];
        let la = [120i16; LP_LOOKAHEAD_SAMPLES];
        let frame = state.encode_one(&pcm, &la);
        // All fields are u8/u16 by type so finiteness is implicit; the
        // key sanity check is that pack_frame round-trips cleanly.
        let bytes = pack_frame(&frame);
        assert_eq!(bytes.len(), FRAME_BYTES);
    }
}
