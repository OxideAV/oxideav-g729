//! §4.1 per-frame **decode parameter chain** — the glue that turns one
//! ITU serial frame into the fully-populated parameter structs every
//! already-landed decode piece consumes.
//!
//! This module wires the round-191 framing layer ([`crate::serial`]),
//! the round-225 Table-8 unpacker ([`crate::parameters`]), and the
//! per-clause decode pieces (§3.2.4 LSP reconstruction, §3.2.5
//! interpolation, §3.2.6 LSP→LP, §4.1.3 pitch-delay decode, §4.1.4
//! fixed-codebook decode + eq (48) sharpening, §3.9 / §4.1.5 gain
//! decode) into one stateful per-frame call.
//!
//! ## Spec source — clause 4.1 (06/2012 Recommendation)
//!
//! Clause 4.1 opens: "The decoding process is done in the following
//! order", and the per-clause breakdown is:
//!
//! * **§4.1.1 — LP filter parameters.** "The received indices L0, L1,
//!   L2 and L3 of the LSP quantizer are used to reconstruct the
//!   quantized LSP coefficients using the procedure described in
//!   clause 3.2.4. The interpolation procedure described in clause
//!   3.2.5 is used to obtain two sets of interpolated LSP coefficients
//!   (corresponding to two subframes). For each subframe, the
//!   interpolated LSP coefficients are converted to LP filter
//!   coefficients a_i." Then "the following steps are repeated for
//!   each subframe: 1) decoding of the adaptive-codebook vector;
//!   2) decoding of the fixed-codebook vector; 3) decoding of the
//!   adaptive and fixed-codebook gains; and 4) computation of the
//!   reconstructed speech."
//! * **§4.1.2 — parity.** "The parity bit is recomputed from the
//!   adaptive-codebook delay index P1 … If a parity error occurs on
//!   P1, the delay value T1 is set to the integer part of the delay
//!   value T2 of the previous frame. The value T2 is derived with the
//!   procedure outlined in clause 4.1.3, using this new value of T1."
//! * **§4.1.3 — adaptive-codebook delays** (`P1 → T1`,
//!   `(P2, t_min) → T2`), wired in [`crate::pitch_decode`].
//! * **§4.1.4 — fixed-codebook vector** (`(C, S) →` pulses, eq (45)
//!   codevector, eq (48) modification when `int(T) < 40`), wired in
//!   [`crate::fixed_codebook`] + [`crate::pitch_sharpen`].
//! * **§4.1.5 — gains.** "The received gain-codebook index gives the
//!   adaptive-codebook gain ĝ_p and the fixed-codebook gain
//!   correction factor [γ̂] … The estimated fixed-codebook gain g′_c
//!   is found using equation (71). The fixed-codebook [gain] is
//!   obtained from the product of the quantized gain correction
//!   factor with this predicted gain [equation (74)]. The
//!   adaptive-codebook gain is reconstructed using equation (73)."
//!   Wired in [`crate::gain_index_map`] + [`crate::gain_reconstruct`]
//!   + [`crate::gain_predict`].
//!
//! §4.1.6 (LP synthesis to PCM) and §4.2 (post-processing) are **not**
//! part of this chain — they are the remaining wire-up rounds. The
//! §4.4 frame-erasure concealment is likewise not wired: an erasure
//! sentinel frame returns [`FrameDecodeError::Erased`].
//!
//! ## Chain state (clause 4.3 / Table 9 initialization)
//!
//! Per clause 4.3, "all static encoder and decoder variables should be
//! initialized to zero, except the variables listed in Table 9". The
//! chain owns four pieces of cross-frame state:
//!
//! | state | role | init |
//! |-------|------|------|
//! | [`crate::lsp_reconstruct::LspReconstructor`] | §3.2.4 4-frame MA residual history | `l̂_i = iπ/11` (Table 9 `î_i`) |
//! | [`crate::lsp_interpolate::LspInterpolator`] | §3.2.5 previous-frame LSPs | `q_i = cos(iπ/11)` (Table 9 `q_i`) |
//! | [`crate::gain_predict::GainPredictor`] | §3.9.1 4-tap `Û` history | `−14 dB` (Table 9 `Û(k)`) |
//! | `g_p_prev` (this module) | eq (47) `β = ĝ_p^(m−1)` source | `0.8` (Table 9 `β`) |
//!
//! plus the previous frame's `int(T2)` for the §4.1.2 parity
//! concealment, which Table 9 does not list and therefore starts at
//! zero per the clause-4.3 default.
//!
//! ## Energy term of the gain prediction
//!
//! The eq (66) codevector energy that feeds the §3.9.1 predictor is
//! computed on the **post-sharpening** codevector: per §3.10, `c(n)`
//! in the excitation equation (75) is "the fixed-codebook vector
//! including harmonic enhancement", i.e. the eq (48)-modified vector,
//! and §3.9.1 defines `E` over that same `c(n)`.

use crate::concealment::Concealer;
use crate::fixed_codebook::{self, FixedCodebookError, FixedCodebookPulses, SUBFRAME_SIZE};
use crate::gain_predict::{GainPredictor, PredictedGain};
use crate::gain_reconstruct::{
    reconstruct_gains_from_transmitted, GainReconstructError, QuantisedGains,
};
use crate::lp_synthesis::{SynthesizedFrame, Synthesizer};
use crate::lsp_interpolate::{omega_to_q, LspInterpolator, SUBFRAMES_PER_FRAME};
use crate::lsp_reconstruct::{LspReconstructError, LspReconstructor};
use crate::lsp_to_lp::{lsp_to_lp, LpCoefficients};
use crate::parameters::{unpack_parameters, ParameterError, Parameters};
use crate::pitch_decode::{decode_t1_from_p1, decode_t2_from_p2, derive_t_min, PitchDelay};
use crate::pitch_sharpen::{clamp_beta, sharpen};
use crate::postfilter_cascade::{PostfilterCascade, PostfilteredFrame};
use crate::serial::{parse_frame, FrameKind, SerialError};
use crate::tables::M;

/// Clause 4.3 / Table 9 initial value of the eq (47) pitch-sharpening
/// gain `β` (spec row: "β | 3.8 | 0.8").
pub const BETA_INIT: f32 = 0.8;

/// Everything the §4.1 procedure decodes for **one subframe**, fully
/// typed. Field order follows the clause-4.1 per-subframe step list.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SubframeDecode {
    /// §3.2.5-interpolated cosine-domain LSP vector `q_i` for this
    /// subframe (§4.1.1).
    pub lsp_q: [f32; M],
    /// §3.2.6 LP filter coefficients `a_i` converted from
    /// [`Self::lsp_q`] (§4.1.1: "for each subframe, the interpolated
    /// LSP coefficients are converted to LP filter coefficients").
    pub lp: LpCoefficients,
    /// §4.1.3 fractional pitch delay `T` for this subframe
    /// (`T1` for subframe 1, `T2` for subframe 2). When a §4.1.2
    /// parity error is detected, subframe 1 carries the concealment
    /// delay (previous frame's `int(T2)`, fraction zero) instead of
    /// the `P1` decode.
    pub pitch: PitchDelay,
    /// §4.1.4 decoded fixed-codebook pulses from `(C, S)`.
    pub fixed: FixedCodebookPulses,
    /// The eq (47) pitch gain `β = ĝ_p^(m−1)` (clamped to
    /// `[0.2, 0.8]`) that the eq (48) modification used. Recorded
    /// even when `int(T) ≥ 40` left the codevector unmodified.
    pub beta: f32,
    /// The fixed-codebook vector `c(n)` after the §4.1.4 eq (48)
    /// modification ("including harmonic enhancement", §3.10) —
    /// the vector that `ĝ_c` scales in eq (75).
    pub codevector: [f32; SUBFRAME_SIZE],
    /// §4.1.5 quantised gain pair: the adaptive-codebook gain `ĝ_p`
    /// (eq (73)) and the fixed-codebook gain correction factor `γ̂`
    /// (eq (74)), reconstructed from the transmitted `(GA, GB)`
    /// through the §3.9.3 inverse permutation and the §3.9.2
    /// conjugate-structure codebooks.
    pub gains: QuantisedGains,
    /// §3.9.1 prediction-path intermediates (`E` from eq (66),
    /// `Ẽ^(m)` from eq (69), `g′_c` from eq (71)) for this subframe.
    pub predicted: PredictedGain,
    /// §4.1.5 quantised fixed-codebook gain `ĝ_c = γ̂ · g′_c`
    /// (eq (74)) — the gain that scales [`Self::codevector`] in the
    /// eq (75) excitation.
    pub g_c_hat: f32,
}

/// Everything the §4.1 procedure decodes for **one 10 ms frame**.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DecodedFrame {
    /// The 15 raw Table-8 codewords the frame transmitted.
    pub params: Parameters,
    /// §4.1.2 parity outcome: `true` when the recomputed parity of
    /// `P1` matches the transmitted `P0`. When `false`, subframe 1's
    /// pitch delay in [`Self::subframes`] is the §4.1.2 concealment
    /// value (previous frame's `int(T2)`), not the `P1` decode.
    pub parity_ok: bool,
    /// §4.1.1 reconstructed LSF vector `ω̂^(m)` (via §3.2.4), before
    /// the §3.2.5 per-subframe interpolation.
    pub omega: [f32; M],
    /// §4.1.3 `t_min` — the lower edge of the subframe-2 delay
    /// window derived from `int(T1)`.
    pub t_min: i32,
    /// Per-subframe decode products in spec order (index 0 =
    /// subframe 1, index 1 = subframe 2).
    pub subframes: [SubframeDecode; SUBFRAMES_PER_FRAME],
}

/// Errors surfaced by the §4.1 decode chain.
///
/// For bits that came off the wire through
/// [`crate::parameters::unpack_parameters`] the codeword-domain
/// variants (`Lsp` / `FixedCodebook` / `Gain`) are unreachable — every
/// Table-8 codeword width matches its codebook domain. They surface
/// only when a caller hands [`FrameDecoder::decode_parameters`] a
/// hand-built [`Parameters`] struct with out-of-domain integers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameDecodeError {
    /// The 164-byte serial framing was malformed
    /// (only from [`FrameDecoder::decode_serial_frame`]).
    Serial(SerialError),
    /// The frame is a §4.4 erasure sentinel. The §4.4 concealment
    /// path is not wired yet; callers must skip / conceal upstream.
    Erased,
    /// §3.2.4 LSP reconstruction rejected a codeword index.
    Lsp(LspReconstructError),
    /// §4.1.4 fixed-codebook decode rejected a codeword.
    FixedCodebook(FixedCodebookError),
    /// §4.1.5 gain decode rejected a codeword index.
    Gain(GainReconstructError),
}

impl core::fmt::Display for FrameDecodeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Serial(e) => write!(f, "g729 §4.1 decode chain: {e}"),
            Self::Erased => write!(
                f,
                "g729 §4.1 decode chain: frame is a §4.4 erasure sentinel; concealment not wired"
            ),
            Self::Lsp(e) => write!(f, "g729 §4.1 decode chain: {e}"),
            Self::FixedCodebook(e) => write!(f, "g729 §4.1 decode chain: {e}"),
            Self::Gain(e) => write!(f, "g729 §4.1 decode chain: {e}"),
        }
    }
}

impl std::error::Error for FrameDecodeError {}

impl From<SerialError> for FrameDecodeError {
    fn from(e: SerialError) -> Self {
        Self::Serial(e)
    }
}

impl From<ParameterError> for FrameDecodeError {
    fn from(e: ParameterError) -> Self {
        match e {
            ParameterError::Erased => Self::Erased,
        }
    }
}

impl From<LspReconstructError> for FrameDecodeError {
    fn from(e: LspReconstructError) -> Self {
        Self::Lsp(e)
    }
}

impl From<FixedCodebookError> for FrameDecodeError {
    fn from(e: FixedCodebookError) -> Self {
        Self::FixedCodebook(e)
    }
}

impl From<GainReconstructError> for FrameDecodeError {
    fn from(e: GainReconstructError) -> Self {
        Self::Gain(e)
    }
}

/// Stateful §4.1 per-frame decode chain.
///
/// Owns the four cross-frame state pieces listed in the module
/// documentation (LSP MA history, LSP interpolation memory, gain-
/// predictor history, previous-subframe `ĝ_p`) plus the previous
/// frame's `int(T2)` for the §4.1.2 parity concealment. Call
/// [`Self::decode_serial_frame`] once per 164-byte ITU serial frame,
/// in stream order.
#[derive(Debug, Clone)]
pub struct FrameDecoder {
    /// §3.2.4 reconstruction state (4-frame MA residual history).
    lsp: LspReconstructor,
    /// §3.2.5 interpolation state (previous frame's cosine-domain LSPs).
    interpolator: LspInterpolator,
    /// §3.9.1 gain-prediction state (4-tap `Û` history).
    gain_predictor: GainPredictor,
    /// `ĝ_p^(m−1)` — the previous subframe's quantised adaptive-
    /// codebook gain, the eq (47) `β` source. Initialised to
    /// [`BETA_INIT`] per clause 4.3 / Table 9.
    g_p_prev: f32,
    /// Previous frame's `int(T2)`, used by the §4.1.2 parity
    /// concealment ("the delay value T1 is set to the integer part of
    /// the delay value T2 of the previous frame"). Not listed in
    /// Table 9, so it starts at zero per the clause-4.3 default.
    prev_int_t2: i32,
    /// §3.10 / §4.1.6 LP synthesizer. Owns the eq (40) past-excitation
    /// buffer and the eq (77) 10th-order filter memory, both zero-init
    /// per clause 4.3. Driven only by the `*_to_speech` entry points;
    /// the parameter-only `decode_*` calls leave it untouched.
    synthesizer: Synthesizer,
    /// §4.2 post-processing cascade (long-term + short-term postfilters,
    /// tilt compensation, adaptive gain control, output high-pass + ×2).
    /// Owns the five stages' clause-4.3 state; driven only by the
    /// `*_to_postfiltered` entry points (the `*_to_speech` calls leave it
    /// untouched).
    postfilter: PostfilterCascade,
    /// §4.4 frame-erasure concealer (voicing class, eq (96) random seed,
    /// carried periodic pitch delay). Advanced only by the `*_concealed`
    /// entry points; the strict `decode_*` / `*_to_speech` /
    /// `*_to_postfiltered` calls (which error on an erasure) leave it at
    /// its clause-4.3 start-up state.
    concealer: Concealer,
    /// The last successfully-decoded **good** frame, retained for the
    /// §4.4.1 LP-parameter repeat and the §4.4.2 gain attenuation that an
    /// erased frame is built from. `None` until the first good frame is
    /// decoded through a `*_concealed` entry point (a leading erasure with
    /// no prior good frame synthesizes silence).
    last_good: Option<DecodedFrame>,
}

impl Default for FrameDecoder {
    fn default() -> Self {
        Self::new()
    }
}

impl FrameDecoder {
    /// Build a chain with the clause 4.3 / Table 9 start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            lsp: LspReconstructor::new(),
            interpolator: LspInterpolator::new(),
            gain_predictor: GainPredictor::new(),
            g_p_prev: BETA_INIT,
            prev_int_t2: 0,
            synthesizer: Synthesizer::new(),
            postfilter: PostfilterCascade::new(),
            concealer: Concealer::new(),
            last_good: None,
        }
    }

    /// Borrow the §4.4 [`Concealer`] for inspection / tests. Only the
    /// `*_concealed` entry points advance it.
    #[must_use]
    pub fn concealer(&self) -> &Concealer {
        &self.concealer
    }

    /// Borrow the current eq (47) `β` source (`ĝ_p^(m−1)`) for
    /// inspection / tests.
    #[must_use]
    pub fn g_p_prev(&self) -> f32 {
        self.g_p_prev
    }

    /// Borrow the previous frame's `int(T2)` (the §4.1.2 parity-
    /// concealment delay) for inspection / tests.
    #[must_use]
    pub fn prev_int_t2(&self) -> i32 {
        self.prev_int_t2
    }

    /// Borrow the §3.10 / §4.1.6 [`Synthesizer`] for inspection / tests.
    ///
    /// The synthesizer only advances when a `*_to_speech` entry point is
    /// used; the parameter-only `decode_*` calls leave it at its
    /// clause-4.3 start-up state.
    #[must_use]
    pub fn synthesizer(&self) -> &Synthesizer {
        &self.synthesizer
    }

    /// Borrow the §4.2 [`PostfilterCascade`] for inspection / tests.
    ///
    /// The cascade only advances when a `*_to_postfiltered` entry point is
    /// used; the `decode_*` and `*_to_speech` calls leave it at its
    /// clause-4.3 start-up state.
    #[must_use]
    pub fn postfilter(&self) -> &PostfilterCascade {
        &self.postfilter
    }

    /// Decode one 164-byte ITU serial frame: framing parse → Table-8
    /// unpack → full §4.1 parameter chain.
    ///
    /// # Errors
    ///
    /// [`FrameDecodeError::Serial`] on malformed framing,
    /// [`FrameDecodeError::Erased`] on a §4.4 erasure sentinel. The
    /// codeword-domain variants cannot occur on a well-formed frame.
    pub fn decode_serial_frame(
        &mut self,
        frame_bytes: &[u8],
    ) -> Result<DecodedFrame, FrameDecodeError> {
        let kind = parse_frame(frame_bytes)?;
        self.decode_frame_kind(&kind)
    }

    /// Decode one already-parsed [`FrameKind`].
    ///
    /// # Errors
    ///
    /// [`FrameDecodeError::Erased`] on [`FrameKind::Erased`]; the
    /// codeword-domain variants cannot occur on a well-formed frame.
    pub fn decode_frame_kind(
        &mut self,
        frame: &FrameKind,
    ) -> Result<DecodedFrame, FrameDecodeError> {
        let params = unpack_parameters(frame)?;
        self.decode_parameters(&params)
    }

    /// Run the §4.1 parameter chain on one frame's unpacked Table-8
    /// codewords, advancing all cross-frame state.
    ///
    /// # Errors
    ///
    /// Surfaces the codeword-domain variants only for a hand-built
    /// [`Parameters`] struct whose integers exceed their Table-8
    /// widths; output from
    /// [`crate::parameters::unpack_parameters`] never fails. On such
    /// an error the cross-frame state may already be partially
    /// advanced (the chain runs in spec §4.1 order and does not
    /// roll back) — callers feeding hand-built structs should treat
    /// the chain as poisoned after a codeword-domain error.
    pub fn decode_parameters(
        &mut self,
        params: &Parameters,
    ) -> Result<DecodedFrame, FrameDecodeError> {
        // --- §4.1.1: LP filter parameters. -------------------------
        // L0..L3 → ω̂^(m) (§3.2.4), interpolate (§3.2.5), convert to
        // LP coefficients per subframe (§3.2.6). The LSP MA history
        // must NOT advance on a codeword-domain error, so reconstruct
        // first and only then touch the interpolator.
        let omega = self.lsp.reconstruct_frame(
            usize::from(params.l0),
            usize::from(params.l1),
            usize::from(params.l2),
            usize::from(params.l3),
        )?;
        let q_current = omega_to_q(&omega);
        let lsp_q = self.interpolator.interpolate(&q_current);
        let lp = [lsp_to_lp(&lsp_q[0]), lsp_to_lp(&lsp_q[1])];

        // --- §4.1.2: parity check + concealment delay. -------------
        // "If a parity error occurs on P1, the delay value T1 is set
        // to the integer part of the delay value T2 of the previous
        // frame. The value T2 is derived with the procedure outlined
        // in clause 4.1.3, using this new value of T1."
        let parity_ok = params.pitch_parity_ok();
        let t1 = if parity_ok {
            decode_t1_from_p1(params.p1)
        } else {
            PitchDelay {
                int_t: self.prev_int_t2,
                frac: 0,
            }
        };

        // --- §4.1.3: adaptive-codebook delays. ---------------------
        let t_min = derive_t_min(t1.int_t);
        let t2 = decode_t2_from_p2(params.p2, t_min);
        let pitch = [t1, t2];

        // --- per-subframe §4.1.4 + §4.1.5. -------------------------
        let mut subframes = [None, None];
        let cs: [(u16, u8, u8, u8); SUBFRAMES_PER_FRAME] = [
            (params.c1, params.s1, params.ga1, params.gb1),
            (params.c2, params.s2, params.ga2, params.gb2),
        ];
        for (i, &(c, s, ga, gb)) in cs.iter().enumerate() {
            // §4.1.4: (C, S) → pulses → eq (45) codevector → eq (48)
            // modification when int(T) < 40, with β = ĝ_p^(m−1)
            // clamped per eq (47).
            let fixed = fixed_codebook::decode_pulses(c, s)?;
            let raw = fixed_codebook::build_codevector(&fixed);
            let beta = clamp_beta(self.g_p_prev);
            let codevector = sharpen(&raw, pitch[i].int_t, self.g_p_prev);

            // §4.1.5: (GA, GB) → §3.9.3 demap → eqs (73)/(74)
            // (ĝ_p, γ̂), then ĝ_c = γ̂ · g′_c with g′_c from the
            // §3.9.1 predictor (eq (71)) over the eq (66) energy of
            // the harmonic-enhanced c(n) (§3.10). The predictor's
            // 4-tap history advances per eq (72).
            let gains = reconstruct_gains_from_transmitted(usize::from(ga), usize::from(gb))?;
            let (g_c_hat, predicted) = self
                .gain_predictor
                .predict_and_update(&codevector, gains.gamma_hat);

            // eq (47) state advance: this subframe's ĝ_p is the next
            // subframe's β source.
            self.g_p_prev = gains.g_p_hat;

            subframes[i] = Some(SubframeDecode {
                lsp_q: lsp_q[i],
                lp: lp[i],
                pitch: pitch[i],
                fixed,
                beta,
                codevector,
                gains,
                predicted,
                g_c_hat,
            });
        }

        // §4.1.2 state advance: this frame's int(T2) is next frame's
        // parity-concealment delay.
        self.prev_int_t2 = t2.int_t;

        let [sub1, sub2] = subframes;
        Ok(DecodedFrame {
            params: *params,
            parity_ok,
            omega,
            t_min,
            subframes: [
                sub1.expect("subframe 1 populated above"),
                sub2.expect("subframe 2 populated above"),
            ],
        })
    }

    /// Decode one 164-byte ITU serial frame all the way to reconstructed
    /// speech `ŝ(n)`.
    ///
    /// Runs the §4.1 parameter chain ([`Self::decode_serial_frame`]) and
    /// then the §3.10 → §4.1.6 LP synthesis (eq (40) adaptive-codebook
    /// interpolation, eq (75) excitation `u(n) = ĝ_p·v(n) + ĝ_c·c(n)`,
    /// eq (77) `1/Â(z)` filter), advancing both the §4.1 cross-frame
    /// state and the synthesizer's eq (40)/(77) memories in stream
    /// order. The returned [`SynthesizedFrame`] exposes the 80
    /// reconstructed-speech samples via
    /// [`SynthesizedFrame::speech`].
    ///
    /// The §4.2 post-processing cascade (postfilters, tilt, AGC, output
    /// high-pass + ×2) is **not** applied here — those modules exist
    /// standalone and are wired separately. This entry point yields the
    /// pre-postfilter reconstructed speech.
    ///
    /// # Errors
    ///
    /// As [`Self::decode_serial_frame`].
    pub fn decode_serial_frame_to_speech(
        &mut self,
        frame_bytes: &[u8],
    ) -> Result<SynthesizedFrame, FrameDecodeError> {
        let kind = parse_frame(frame_bytes)?;
        self.decode_frame_kind_to_speech(&kind)
    }

    /// Decode one already-parsed [`FrameKind`] to reconstructed speech.
    ///
    /// See [`Self::decode_serial_frame_to_speech`] for the chain.
    ///
    /// # Errors
    ///
    /// As [`Self::decode_frame_kind`].
    pub fn decode_frame_kind_to_speech(
        &mut self,
        frame: &FrameKind,
    ) -> Result<SynthesizedFrame, FrameDecodeError> {
        let params = unpack_parameters(frame)?;
        self.decode_parameters_to_speech(&params)
    }

    /// Run the §4.1 parameter chain followed by the §3.10 → §4.1.6 LP
    /// synthesis on one frame's unpacked Table-8 codewords, advancing
    /// all cross-frame state plus the synthesizer's memories.
    ///
    /// See [`Self::decode_serial_frame_to_speech`] for the chain.
    ///
    /// # Errors
    ///
    /// As [`Self::decode_parameters`]. On a codeword-domain error the
    /// synthesizer is left untouched (synthesis runs only after a
    /// successful parameter decode).
    pub fn decode_parameters_to_speech(
        &mut self,
        params: &Parameters,
    ) -> Result<SynthesizedFrame, FrameDecodeError> {
        let frame = self.decode_parameters(params)?;
        Ok(self.synthesizer.synthesize_frame(&frame))
    }

    /// As [`Self::decode_parameters_to_speech`], additionally returning
    /// the §4.1 [`DecodedFrame`] alongside the synthesized speech — the
    /// pair an external post-processing cascade (e.g. the Annex A
    /// [`crate::annex_a::AnnexAPostfilterCascade`]) needs to run in
    /// place of the built-in §4.2 cascade.
    ///
    /// # Errors
    ///
    /// As [`Self::decode_parameters`].
    pub fn decode_parameters_with_speech(
        &mut self,
        params: &Parameters,
    ) -> Result<(DecodedFrame, SynthesizedFrame), FrameDecodeError> {
        let frame = self.decode_parameters(params)?;
        let speech = self.synthesizer.synthesize_frame(&frame);
        Ok((frame, speech))
    }

    /// Decode one 164-byte ITU serial frame all the way to the §4.2
    /// post-processed decoder output `sf′(n)`.
    ///
    /// Runs the §4.1 parameter chain, the §3.10 → §4.1.6 LP synthesis, and
    /// then the full §4.2 post-processing cascade (long-term postfilter,
    /// short-term postfilter, tilt compensation, adaptive gain control,
    /// output high-pass + ×2 up-scaling — see [`PostfilterCascade`]),
    /// advancing the §4.1 cross-frame state, the synthesizer memories, and
    /// the five cascade stages in stream order. The returned
    /// [`PostfilteredFrame`] exposes the 80 ×2-upscaled output samples via
    /// [`PostfilteredFrame::output`].
    ///
    /// The §4.2.1 long-term postfilter uses its full **two-pass** form
    /// (integer anchor + 1/8-resolution fractional second pass; see
    /// [`crate::long_term_postfilter`]).
    ///
    /// # Errors
    ///
    /// As [`Self::decode_serial_frame`].
    pub fn decode_serial_frame_to_postfiltered(
        &mut self,
        frame_bytes: &[u8],
    ) -> Result<PostfilteredFrame, FrameDecodeError> {
        let kind = parse_frame(frame_bytes)?;
        self.decode_frame_kind_to_postfiltered(&kind)
    }

    /// Decode one already-parsed [`FrameKind`] to the §4.2 post-processed
    /// decoder output.
    ///
    /// See [`Self::decode_serial_frame_to_postfiltered`] for the chain.
    ///
    /// # Errors
    ///
    /// As [`Self::decode_frame_kind`].
    pub fn decode_frame_kind_to_postfiltered(
        &mut self,
        frame: &FrameKind,
    ) -> Result<PostfilteredFrame, FrameDecodeError> {
        let params = unpack_parameters(frame)?;
        self.decode_parameters_to_postfiltered(&params)
    }

    /// Run the §4.1 parameter chain, the §3.10 → §4.1.6 LP synthesis, and
    /// the §4.2 post-processing cascade on one frame's unpacked Table-8
    /// codewords, advancing all cross-frame state.
    ///
    /// See [`Self::decode_serial_frame_to_postfiltered`] for the chain.
    ///
    /// # Errors
    ///
    /// As [`Self::decode_parameters`]. On a codeword-domain error neither
    /// the synthesizer nor the postfilter cascade is touched (both run
    /// only after a successful parameter decode).
    pub fn decode_parameters_to_postfiltered(
        &mut self,
        params: &Parameters,
    ) -> Result<PostfilteredFrame, FrameDecodeError> {
        let frame = self.decode_parameters(params)?;
        let speech = self.synthesizer.synthesize_frame(&frame);
        Ok(self.postfilter.process_frame(&speech, &frame))
    }

    // ---------------------------------------------------------------
    // §4.4 frame-erasure concealment — the `*_concealed` entry points.
    //
    // Unlike the strict `*_to_postfiltered` / `*_to_speech` calls (which
    // return `FrameDecodeError::Erased` on a §4.4 erasure sentinel so the
    // caller can decide its own policy), these entry points implement the
    // clause-4.4 concealment: an erased frame is reconstructed from the
    // last good frame's parameters with the §4.4.2 gain attenuation and
    // §4.4.4 replacement excitation, then run through the same §4.1.6
    // synthesis + §4.2 cascade as a good frame. A good frame decodes
    // normally and additionally latches the concealer's voicing class
    // (§4.4) from the §4.2.1 long-term-postfilter decisions and retains
    // itself as the new last-good frame.
    // ---------------------------------------------------------------

    /// Decode one 164-byte ITU serial frame to §4.2 post-processed output,
    /// applying §4.4 **frame-erasure concealment** on an erasure sentinel
    /// instead of erroring.
    ///
    /// A good frame decodes exactly as
    /// [`Self::decode_serial_frame_to_postfiltered`] and additionally
    /// updates the §4.4 concealer (voicing class + carried pitch delay)
    /// and the retained last-good frame. An erased frame
    /// ([`FrameKind::Erased`]) is concealed per clause 4.4: the §4.4.1 LP
    /// repeat, the §4.4.2 gain attenuation (eqs (93)/(94)), the §4.4.3
    /// gain-predictor-memory attenuation (eq (95)), and the §4.4.4
    /// replacement excitation (periodic repeat or eq (96) random fixed
    /// codebook by the inherited voicing class).
    ///
    /// A leading erasure before any good frame has been seen synthesizes
    /// from the clause-4.3 zero state (silence), matching the reference
    /// decoder's start-up behaviour.
    ///
    /// # Errors
    ///
    /// [`FrameDecodeError::Serial`] on malformed framing. A well-formed
    /// active or erased frame never errors here (the erasure that
    /// `*_to_postfiltered` would reject is concealed).
    pub fn decode_serial_frame_to_postfiltered_concealed(
        &mut self,
        frame_bytes: &[u8],
    ) -> Result<PostfilteredFrame, FrameDecodeError> {
        let kind = parse_frame(frame_bytes)?;
        match kind {
            FrameKind::Erased => Ok(self.conceal_to_postfiltered()),
            FrameKind::Active(_) => {
                let params = unpack_parameters(&kind)?;
                let frame = self.decode_parameters(&params)?;
                let speech = self.synthesizer.synthesize_frame(&frame);
                let pf = self.postfilter.process_frame(&speech, &frame);
                self.observe_good_frame(&frame, &pf);
                Ok(pf)
            }
        }
    }

    /// Decode one 164-byte ITU serial frame to the §4.1.6 reconstructed
    /// speech `ŝ(n)`, applying §4.4 concealment on an erasure sentinel.
    ///
    /// The reconstructed-speech analogue of
    /// [`Self::decode_serial_frame_to_postfiltered_concealed`] — runs the
    /// §4.1.6 synthesis but not the §4.2 cascade. Because the §4.4 voicing
    /// class is derived from the §4.2.1 long-term-postfilter decisions, a
    /// good frame here runs the cascade internally **only** to obtain that
    /// classification (its postfiltered output is discarded); the returned
    /// value is the pre-postfilter `ŝ(n)`.
    ///
    /// # Errors
    ///
    /// As [`Self::decode_serial_frame_to_postfiltered_concealed`].
    pub fn decode_serial_frame_to_speech_concealed(
        &mut self,
        frame_bytes: &[u8],
    ) -> Result<SynthesizedFrame, FrameDecodeError> {
        let kind = parse_frame(frame_bytes)?;
        match kind {
            FrameKind::Erased => Ok(self.conceal_to_speech()),
            FrameKind::Active(_) => {
                let params = unpack_parameters(&kind)?;
                let frame = self.decode_parameters(&params)?;
                let speech = self.synthesizer.synthesize_frame(&frame);
                // Classify via the cascade (postfiltered output discarded)
                // so the concealer's §4.4 voicing class stays in sync.
                let pf = self.postfilter.process_frame(&speech, &frame);
                self.observe_good_frame(&frame, &pf);
                Ok(speech)
            }
        }
    }

    /// Latch the §4.4 concealer's voicing class + carried pitch delay from
    /// a good frame's §4.2.1 long-term decisions, and retain the frame as
    /// the new last-good frame for any subsequent erasure.
    fn observe_good_frame(&mut self, frame: &DecodedFrame, pf: &PostfilteredFrame) {
        let decisions = [pf.subframes[0].long_term, pf.subframes[1].long_term];
        // §4.4.4 periodic repeat starts from int(T1) of the good frame.
        let int_t1 = frame.subframes[0].pitch.int_t.max(0) as usize;
        self.concealer.observe_good_frame(&decisions, int_t1);
        self.last_good = Some(*frame);
    }

    /// Build a concealed [`DecodedFrame`] for an erased frame from the
    /// retained last-good frame (clause 4.4). Returns `None` before any
    /// good frame has been seen (leading-erasure silence).
    fn build_concealed_frame(&mut self) -> Option<DecodedFrame> {
        let mut frame = self.last_good?;
        let periodic = self.concealer.is_periodic();

        for sub in &mut frame.subframes {
            // §4.4.2 gain attenuation (eqs (93)/(94)): decay the previous
            // subframe's gains.
            let g_p = Concealer::attenuate_adaptive_gain(sub.gains.g_p_hat);
            let g_c = Concealer::attenuate_fixed_gain(sub.g_c_hat);

            if periodic {
                // §4.4.4 periodic case: repeat the adaptive codebook at the
                // carried pitch delay (+1 per subframe, bounded 143) and
                // drop the fixed-codebook contribution.
                let delay = self.concealer.periodic_pitch_delay();
                sub.pitch = PitchDelay {
                    int_t: delay as i32,
                    frac: 0,
                };
                sub.codevector = [0.0; SUBFRAME_SIZE];
                sub.g_c_hat = 0.0;
            } else {
                // §4.4.4 non-periodic case: eq (96) random fixed codebook,
                // no adaptive contribution.
                let rnd = self.concealer.next_random_excitation();
                if let Ok(pulses) = fixed_codebook::decode_pulses(rnd.index, rnd.sign as u8) {
                    let raw = fixed_codebook::build_codevector(&pulses);
                    sub.fixed = pulses;
                    sub.codevector = sharpen(&raw, sub.pitch.int_t, self.g_p_prev);
                }
                sub.g_c_hat = g_c;
            }
            sub.gains.g_p_hat = if periodic { g_p } else { 0.0 };

            // §4.4.3 gain-predictor-memory attenuation (eq (95)): push the
            // averaged-minus-4 dB value onto the predictor so a sustained
            // erasure run keeps decaying the predicted gain.
            let u_db = Concealer::attenuated_predictor_memory_db(&self.gain_predictor);
            self.gain_predictor
                .push_quantised_error(10.0_f32.powf(u_db / 20.0));

            // Advance the eq (47) β source with the attenuated adaptive
            // gain, matching the good-frame state thread.
            self.g_p_prev = sub.gains.g_p_hat;
        }

        // Retain the concealed frame as the basis for a further erasure
        // (the §4.4.2 attenuation compounds across an erasure run).
        self.last_good = Some(frame);
        Some(frame)
    }

    /// Synthesize + postfilter a concealed erased frame (or silence before
    /// the first good frame).
    fn conceal_to_postfiltered(&mut self) -> PostfilteredFrame {
        match self.build_concealed_frame() {
            Some(frame) => {
                let speech = self.synthesizer.synthesize_frame(&frame);
                self.postfilter.process_frame(&speech, &frame)
            }
            None => {
                // Leading erasure: silence frame from the zero state.
                let silent = self.silence_frame();
                let speech = self.synthesizer.synthesize_frame(&silent);
                self.postfilter.process_frame(&speech, &silent)
            }
        }
    }

    /// Synthesize a concealed erased frame to reconstructed speech.
    fn conceal_to_speech(&mut self) -> SynthesizedFrame {
        let frame = self
            .build_concealed_frame()
            .unwrap_or_else(|| self.silence_frame());
        self.synthesizer.synthesize_frame(&frame)
    }

    /// A zero-excitation frame used for a leading erasure (no last-good
    /// frame yet). The LP filter is the identity (`â_i = 0`), so the
    /// synthesizer emits silence — matching the clause-4.3 zero start-up.
    fn silence_frame(&self) -> DecodedFrame {
        let zero_sub = SubframeDecode {
            lsp_q: [0.0; M],
            lp: [0.0; M],
            pitch: PitchDelay { int_t: 40, frac: 0 },
            fixed: FixedCodebookPulses::default(),
            beta: BETA_INIT,
            codevector: [0.0; SUBFRAME_SIZE],
            gains: QuantisedGains {
                g_p_hat: 0.0,
                gamma_hat: 0.0,
            },
            predicted: PredictedGain {
                e_db: 0.0,
                e_tilde_db: 0.0,
                g_c_prime: 0.0,
            },
            g_c_hat: 0.0,
        };
        DecodedFrame {
            params: Parameters::default(),
            parity_ok: true,
            omega: [0.0; M],
            t_min: 20,
            subframes: [zero_sub, zero_sub],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gain_predict::GainPredictor;
    use crate::parameters::FRAME_BITS;
    use crate::pitch_sharpen::{BETA_MAX, BETA_MIN};
    use crate::serial::{BITS_HEADER, BIT_ONE, BIT_ZERO, SYNC_WORD};

    /// A frame whose pitch parity holds under the corpus-validated
    /// odd-parity rule: P1 = 0b1100_0000 → six MSBs XOR-reduce to 0
    /// → P0 = 1.
    fn parity_ok_params() -> Parameters {
        Parameters {
            l0: 1,
            l1: 17,
            l2: 5,
            l3: 9,
            p1: 0b1100_0000,
            p0: 1,
            c1: 0x0AB3,
            s1: 0b1010,
            ga1: 3,
            gb1: 7,
            p2: 11,
            c2: 0x15C2,
            s2: 0b0101,
            ga2: 6,
            gb2: 12,
        }
    }

    /// Same frame with P0 flipped → parity mismatch.
    fn parity_bad_params() -> Parameters {
        let mut p = parity_ok_params();
        p.p0 = 0;
        p
    }

    /// Pack a [`Parameters`] struct into a 164-byte ITU serial frame
    /// (sync word + bit-count header + 80 bit-words), mirroring the
    /// inline packers used by the serial-entry-point tests.
    fn serial_bytes(params: &Parameters) -> Vec<u8> {
        let vals: [(u16, usize); 15] = [
            (u16::from(params.l0), 1),
            (u16::from(params.l1), 7),
            (u16::from(params.l2), 5),
            (u16::from(params.l3), 5),
            (u16::from(params.p1), 8),
            (u16::from(params.p0), 1),
            (params.c1, 13),
            (u16::from(params.s1), 4),
            (u16::from(params.ga1), 3),
            (u16::from(params.gb1), 4),
            (u16::from(params.p2), 5),
            (params.c2, 13),
            (u16::from(params.s2), 4),
            (u16::from(params.ga2), 3),
            (u16::from(params.gb2), 4),
        ];
        let mut buf = Vec::with_capacity(164);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for &(v, w) in &vals {
            for k in 0..w {
                let bit = ((v >> (w - 1 - k)) & 1) == 1;
                let word = if bit { BIT_ONE } else { BIT_ZERO };
                buf.extend_from_slice(&word.to_le_bytes());
            }
        }
        buf
    }

    /// A 164-byte erasure-sentinel serial frame (the all-`0x0000`-payload
    /// §4.4.5 shape) the concealment path must conceal rather than error
    /// on.
    fn erased_serial_bytes() -> Vec<u8> {
        let mut buf = Vec::with_capacity(164);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for _ in 0..crate::parameters::FRAME_BITS {
            buf.extend_from_slice(&crate::serial::BIT_ERASED.to_le_bytes());
        }
        buf
    }

    /// Clause 4.3 / Table 9 start-up state of the chain-owned slots.
    #[test]
    fn new_chain_has_table9_startup_state() {
        let d = FrameDecoder::new();
        assert_eq!(d.g_p_prev(), BETA_INIT);
        assert_eq!(d.prev_int_t2(), 0);
    }

    /// The chain's per-piece outputs equal calling the pieces by hand
    /// in the §4.1 order with the same start-up state — the chain
    /// adds sequencing + state threading, never different arithmetic.
    #[test]
    fn chain_matches_hand_sequenced_pieces() {
        let params = parity_ok_params();
        let mut chain = FrameDecoder::new();
        let frame = chain.decode_parameters(&params).expect("in-domain");

        // Hand-sequenced reference.
        let mut lsp = LspReconstructor::new();
        let omega = lsp
            .reconstruct_frame(1, 17, 5, 9)
            .expect("in-domain LSP codewords");
        assert_eq!(frame.omega, omega);

        let mut interp = LspInterpolator::new();
        let lsp_q = interp.interpolate(&omega_to_q(&omega));
        assert_eq!(frame.subframes[0].lsp_q, lsp_q[0]);
        assert_eq!(frame.subframes[1].lsp_q, lsp_q[1]);
        assert_eq!(frame.subframes[0].lp, lsp_to_lp(&lsp_q[0]));
        assert_eq!(frame.subframes[1].lp, lsp_to_lp(&lsp_q[1]));

        // §4.1.3 with parity OK.
        assert!(frame.parity_ok);
        let t1 = decode_t1_from_p1(params.p1);
        let t_min = derive_t_min(t1.int_t);
        let t2 = decode_t2_from_p2(params.p2, t_min);
        assert_eq!(frame.subframes[0].pitch, t1);
        assert_eq!(frame.subframes[1].pitch, t2);
        assert_eq!(frame.t_min, t_min);

        // §4.1.4 subframe 1: β source is the Table 9 init 0.8.
        let fixed1 = fixed_codebook::decode_pulses(params.c1, params.s1).expect("in-domain");
        assert_eq!(frame.subframes[0].fixed, fixed1);
        assert_eq!(frame.subframes[0].beta, BETA_INIT);
        let cv1 = sharpen(
            &fixed_codebook::build_codevector(&fixed1),
            t1.int_t,
            BETA_INIT,
        );
        assert_eq!(frame.subframes[0].codevector, cv1);

        // §4.1.5 subframe 1.
        let gains1 =
            reconstruct_gains_from_transmitted(usize::from(params.ga1), usize::from(params.gb1))
                .expect("in-domain");
        assert_eq!(frame.subframes[0].gains, gains1);
        let mut predictor = GainPredictor::new();
        let (g_c1, path1) = predictor.predict_and_update(&cv1, gains1.gamma_hat);
        assert_eq!(frame.subframes[0].g_c_hat, g_c1);
        assert_eq!(frame.subframes[0].predicted, path1);

        // §4.1.4 / §4.1.5 subframe 2: β source is subframe 1's ĝ_p.
        let fixed2 = fixed_codebook::decode_pulses(params.c2, params.s2).expect("in-domain");
        assert_eq!(
            frame.subframes[1].beta,
            clamp_beta(gains1.g_p_hat),
            "subframe-2 β must come from subframe-1 ĝ_p per eq (47)",
        );
        let cv2 = sharpen(
            &fixed_codebook::build_codevector(&fixed2),
            t2.int_t,
            gains1.g_p_hat,
        );
        assert_eq!(frame.subframes[1].codevector, cv2);
        let gains2 =
            reconstruct_gains_from_transmitted(usize::from(params.ga2), usize::from(params.gb2))
                .expect("in-domain");
        let (g_c2, path2) = predictor.predict_and_update(&cv2, gains2.gamma_hat);
        assert_eq!(frame.subframes[1].g_c_hat, g_c2);
        assert_eq!(frame.subframes[1].predicted, path2);

        // Cross-frame state advance.
        assert_eq!(chain.g_p_prev(), gains2.g_p_hat);
        assert_eq!(chain.prev_int_t2(), t2.int_t);
    }

    /// §4.1.2 concealment: on a parity mismatch, subframe 1's delay
    /// is the previous frame's int(T2) with zero fraction, and T2 is
    /// re-derived from that delay.
    #[test]
    fn parity_mismatch_substitutes_previous_frame_t2() {
        let mut chain = FrameDecoder::new();
        // Frame 1 (parity OK) establishes prev_int_t2.
        let f1 = chain
            .decode_parameters(&parity_ok_params())
            .expect("in-domain");
        let expected_conceal = f1.subframes[1].pitch.int_t;
        assert_eq!(chain.prev_int_t2(), expected_conceal);

        // Frame 2 (parity mismatch) must use it.
        let f2 = chain
            .decode_parameters(&parity_bad_params())
            .expect("in-domain");
        assert!(!f2.parity_ok);
        assert_eq!(
            f2.subframes[0].pitch,
            PitchDelay {
                int_t: expected_conceal,
                frac: 0,
            },
        );
        // T2 re-derived per §4.1.3 from the concealment T1.
        let t_min = derive_t_min(expected_conceal);
        assert_eq!(f2.t_min, t_min);
        assert_eq!(
            f2.subframes[1].pitch,
            decode_t2_from_p2(parity_bad_params().p2, t_min),
        );
    }

    /// First-frame parity mismatch falls back to the zero-initialised
    /// previous int(T2) (clause 4.3: unlisted variables start at
    /// zero); the t_min derivation clamps it into [20, 134] and the
    /// non-positive delay leaves the codevector unsharpened.
    #[test]
    fn first_frame_parity_mismatch_uses_zero_init_delay() {
        let mut chain = FrameDecoder::new();
        let f = chain
            .decode_parameters(&parity_bad_params())
            .expect("in-domain");
        assert!(!f.parity_ok);
        assert_eq!(f.subframes[0].pitch, PitchDelay { int_t: 0, frac: 0 });
        assert_eq!(f.t_min, derive_t_min(0));
        assert_eq!(f.t_min, 20);
    }

    /// An erasure sentinel is rejected with the dedicated variant —
    /// no state is advanced.
    #[test]
    fn erased_frame_is_rejected_without_state_advance() {
        let mut chain = FrameDecoder::new();
        let before_gp = chain.g_p_prev();
        let before_t2 = chain.prev_int_t2();
        let err = chain.decode_frame_kind(&FrameKind::Erased).unwrap_err();
        assert_eq!(err, FrameDecodeError::Erased);
        assert_eq!(chain.g_p_prev(), before_gp);
        assert_eq!(chain.prev_int_t2(), before_t2);
    }

    /// Hand-built out-of-domain codewords surface the typed variants
    /// (unreachable from `unpack_parameters` output).
    #[test]
    fn out_of_domain_hand_built_parameters_surface_typed_errors() {
        let mut chain = FrameDecoder::new();
        let mut p = parity_ok_params();
        p.l1 = 200; // >= NC0 = 128
        assert!(matches!(
            chain.decode_parameters(&p).unwrap_err(),
            FrameDecodeError::Lsp(_),
        ));

        let mut p = parity_ok_params();
        p.c1 = 0x3FFF; // > 13 bits
        assert!(matches!(
            chain.decode_parameters(&p).unwrap_err(),
            FrameDecodeError::FixedCodebook(_),
        ));

        let mut p = parity_ok_params();
        p.ga1 = 9; // > 3 bits
        assert!(matches!(
            chain.decode_parameters(&p).unwrap_err(),
            FrameDecodeError::Gain(_),
        ));
    }

    /// `decode_serial_frame` on a synthetic 164-byte frame matches
    /// `decode_parameters` on the equivalent unpacked codewords.
    #[test]
    fn serial_entry_point_matches_parameter_entry_point() {
        // Pack `parity_ok_params()` into a Table-8 bit array, then
        // into the 164-byte serial framing.
        let params = parity_ok_params();
        let vals: [(u16, usize); 15] = [
            (u16::from(params.l0), 1),
            (u16::from(params.l1), 7),
            (u16::from(params.l2), 5),
            (u16::from(params.l3), 5),
            (u16::from(params.p1), 8),
            (u16::from(params.p0), 1),
            (params.c1, 13),
            (u16::from(params.s1), 4),
            (u16::from(params.ga1), 3),
            (u16::from(params.gb1), 4),
            (u16::from(params.p2), 5),
            (params.c2, 13),
            (u16::from(params.s2), 4),
            (u16::from(params.ga2), 3),
            (u16::from(params.gb2), 4),
        ];
        let mut bits = Vec::with_capacity(FRAME_BITS);
        for &(v, w) in &vals {
            for k in 0..w {
                bits.push(((v >> (w - 1 - k)) & 1) == 1);
            }
        }
        assert_eq!(bits.len(), FRAME_BITS);
        let mut buf = Vec::with_capacity(164);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for b in &bits {
            let w = if *b { BIT_ONE } else { BIT_ZERO };
            buf.extend_from_slice(&w.to_le_bytes());
        }

        let mut chain_a = FrameDecoder::new();
        let from_serial = chain_a.decode_serial_frame(&buf).expect("well-formed");
        let mut chain_b = FrameDecoder::new();
        let from_params = chain_b.decode_parameters(&params).expect("in-domain");
        assert_eq!(from_serial, from_params);
        assert_eq!(from_serial.params, params);
    }

    /// The recorded per-subframe β always lies in the eq (47) clamp
    /// range, whatever the gain history does.
    #[test]
    fn beta_always_within_eq47_clamp() {
        let mut chain = FrameDecoder::new();
        for ga in 0..8u8 {
            let mut p = parity_ok_params();
            p.ga1 = ga;
            p.ga2 = 7 - ga;
            let f = chain.decode_parameters(&p).expect("in-domain");
            for sub in &f.subframes {
                assert!(
                    (BETA_MIN..=BETA_MAX).contains(&sub.beta),
                    "β = {} outside eq (47) clamp",
                    sub.beta,
                );
            }
        }
    }

    /// `decode_parameters_to_speech` must equal running the §4.1 chain
    /// then handing the [`DecodedFrame`] to a standalone [`Synthesizer`]:
    /// the speech entry point is exactly that composition, with the
    /// synthesizer state living inside the chain.
    #[test]
    fn speech_entry_point_matches_decode_then_synthesize() {
        use crate::lp_synthesis::Synthesizer;

        let params = parity_ok_params();

        let mut chain = FrameDecoder::new();
        let speech = chain
            .decode_parameters_to_speech(&params)
            .expect("in-domain");

        let mut ref_chain = FrameDecoder::new();
        let frame = ref_chain.decode_parameters(&params).expect("in-domain");
        let mut synth = Synthesizer::new();
        let ref_speech = synth.synthesize_frame(&frame);

        assert_eq!(speech, ref_speech);
        // The parameter-level cross-frame state advanced identically.
        assert_eq!(chain.g_p_prev(), ref_chain.g_p_prev());
        assert_eq!(chain.prev_int_t2(), ref_chain.prev_int_t2());
    }

    /// The 80 reconstructed-speech samples are finite for an in-domain
    /// frame, and the speech entry point advances the chain's synthesizer
    /// memory off its clause-4.3 all-zero start-up state.
    #[test]
    fn speech_entry_point_produces_finite_samples_and_advances_state() {
        let mut chain = FrameDecoder::new();
        assert!(chain.synthesizer().exc_history().iter().all(|&x| x == 0.0));

        let speech = chain
            .decode_parameters_to_speech(&parity_ok_params())
            .expect("in-domain");
        let pcm = speech.speech();
        assert_eq!(pcm.len(), 2 * SUBFRAME_SIZE);
        assert!(pcm.iter().all(|s| s.is_finite()));

        // A non-silent frame leaves a non-zero excitation tail in the
        // synthesizer's eq (40) history.
        assert!(chain.synthesizer().exc_history().iter().any(|&x| x != 0.0));
    }

    /// Two consecutive frames through the speech path are not identical
    /// in general because the synthesizer carries state — and the second
    /// frame still matches the standalone decode-then-synthesize
    /// composition that shares the same advancing state.
    #[test]
    fn speech_path_carries_synthesizer_state_across_frames() {
        use crate::lp_synthesis::Synthesizer;

        let params = parity_ok_params();

        let mut chain = FrameDecoder::new();
        let _ = chain.decode_parameters_to_speech(&params).expect("frame 1");
        let speech2 = chain.decode_parameters_to_speech(&params).expect("frame 2");

        let mut ref_chain = FrameDecoder::new();
        let mut synth = Synthesizer::new();
        let f1 = ref_chain.decode_parameters(&params).expect("frame 1");
        let _ = synth.synthesize_frame(&f1);
        let f2 = ref_chain.decode_parameters(&params).expect("frame 2");
        let ref_speech2 = synth.synthesize_frame(&f2);

        assert_eq!(speech2, ref_speech2);
    }

    /// The serial / frame-kind / parameter speech entry points all agree
    /// for the same well-formed frame.
    #[test]
    fn speech_serial_entry_point_matches_parameter_entry_point() {
        let params = parity_ok_params();
        let vals: [(u16, usize); 15] = [
            (u16::from(params.l0), 1),
            (u16::from(params.l1), 7),
            (u16::from(params.l2), 5),
            (u16::from(params.l3), 5),
            (u16::from(params.p1), 8),
            (u16::from(params.p0), 1),
            (params.c1, 13),
            (u16::from(params.s1), 4),
            (u16::from(params.ga1), 3),
            (u16::from(params.gb1), 4),
            (u16::from(params.p2), 5),
            (params.c2, 13),
            (u16::from(params.s2), 4),
            (u16::from(params.ga2), 3),
            (u16::from(params.gb2), 4),
        ];
        let mut bits = Vec::with_capacity(FRAME_BITS);
        for &(v, w) in &vals {
            for k in 0..w {
                bits.push(((v >> (w - 1 - k)) & 1) == 1);
            }
        }
        let mut buf = Vec::with_capacity(164);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for b in &bits {
            let w = if *b { BIT_ONE } else { BIT_ZERO };
            buf.extend_from_slice(&w.to_le_bytes());
        }

        let mut chain_a = FrameDecoder::new();
        let from_serial = chain_a
            .decode_serial_frame_to_speech(&buf)
            .expect("well-formed");
        let mut chain_b = FrameDecoder::new();
        let from_params = chain_b
            .decode_parameters_to_speech(&params)
            .expect("in-domain");
        assert_eq!(from_serial, from_params);
    }

    /// An erased frame is rejected by the speech path without advancing
    /// the synthesizer state (synthesis runs only after a successful
    /// parameter decode).
    #[test]
    fn speech_path_rejects_erased_frame_without_advancing_synthesizer() {
        let mut chain = FrameDecoder::new();
        let err = chain
            .decode_frame_kind_to_speech(&FrameKind::Erased)
            .unwrap_err();
        assert_eq!(err, FrameDecodeError::Erased);
        assert!(chain.synthesizer().exc_history().iter().all(|&x| x == 0.0));
        assert!(chain.synthesizer().syn_mem().iter().all(|&x| x == 0.0));
    }

    /// `decode_parameters_to_postfiltered` equals running the §4.1 chain,
    /// the synthesizer, then the §4.2 cascade by hand — the cascade entry
    /// point is exactly that composition, with the cascade state living
    /// inside the chain.
    #[test]
    fn postfiltered_entry_point_matches_decode_synthesize_postfilter() {
        use crate::lp_synthesis::Synthesizer;
        use crate::postfilter_cascade::PostfilterCascade;

        let params = parity_ok_params();

        let mut chain = FrameDecoder::new();
        let pf = chain
            .decode_parameters_to_postfiltered(&params)
            .expect("in-domain");

        let mut ref_chain = FrameDecoder::new();
        let frame = ref_chain.decode_parameters(&params).expect("in-domain");
        let mut synth = Synthesizer::new();
        let speech = synth.synthesize_frame(&frame);
        let mut cascade = PostfilterCascade::new();
        let ref_pf = cascade.process_frame(&speech, &frame);

        assert_eq!(pf, ref_pf);
    }

    /// The 80 post-filtered output samples are finite for an in-domain
    /// frame, and the cascade entry point advances the chain's cascade
    /// state off its clause-4.3 all-zero start-up state.
    #[test]
    fn postfiltered_entry_point_produces_finite_output_and_advances_state() {
        let mut chain = FrameDecoder::new();
        assert!(chain
            .postfilter()
            .long_term()
            .s_hist()
            .iter()
            .all(|&x| x == 0.0));

        let pf = chain
            .decode_parameters_to_postfiltered(&parity_ok_params())
            .expect("in-domain");
        let out = pf.output();
        assert_eq!(out.len(), 2 * SUBFRAME_SIZE);
        assert!(out.iter().all(|s| s.is_finite()));

        // A non-silent frame leaves a non-zero history in the long-term
        // postfilter stage.
        assert!(chain
            .postfilter()
            .long_term()
            .s_hist()
            .iter()
            .any(|&x| x != 0.0));
    }

    /// The serial / frame-kind / parameter post-filtered entry points all
    /// agree for the same well-formed frame.
    #[test]
    fn postfiltered_serial_entry_point_matches_parameter_entry_point() {
        let params = parity_ok_params();
        let vals: [(u16, usize); 15] = [
            (u16::from(params.l0), 1),
            (u16::from(params.l1), 7),
            (u16::from(params.l2), 5),
            (u16::from(params.l3), 5),
            (u16::from(params.p1), 8),
            (u16::from(params.p0), 1),
            (params.c1, 13),
            (u16::from(params.s1), 4),
            (u16::from(params.ga1), 3),
            (u16::from(params.gb1), 4),
            (u16::from(params.p2), 5),
            (params.c2, 13),
            (u16::from(params.s2), 4),
            (u16::from(params.ga2), 3),
            (u16::from(params.gb2), 4),
        ];
        let mut bits = Vec::with_capacity(FRAME_BITS);
        for &(v, w) in &vals {
            for k in 0..w {
                bits.push(((v >> (w - 1 - k)) & 1) == 1);
            }
        }
        let mut buf = Vec::with_capacity(164);
        buf.extend_from_slice(&SYNC_WORD.to_le_bytes());
        buf.extend_from_slice(&BITS_HEADER.to_le_bytes());
        for b in &bits {
            let w = if *b { BIT_ONE } else { BIT_ZERO };
            buf.extend_from_slice(&w.to_le_bytes());
        }

        let mut chain_a = FrameDecoder::new();
        let from_serial = chain_a
            .decode_serial_frame_to_postfiltered(&buf)
            .expect("well-formed");
        let mut chain_b = FrameDecoder::new();
        let from_params = chain_b
            .decode_parameters_to_postfiltered(&params)
            .expect("in-domain");
        assert_eq!(from_serial, from_params);
    }

    /// An erased frame is rejected by the post-filter path without
    /// advancing the cascade state (the cascade runs only after a
    /// successful parameter decode).
    #[test]
    fn postfiltered_path_rejects_erased_frame_without_advancing_cascade() {
        let mut chain = FrameDecoder::new();
        let err = chain
            .decode_frame_kind_to_postfiltered(&FrameKind::Erased)
            .unwrap_err();
        assert_eq!(err, FrameDecodeError::Erased);
        assert!(chain
            .postfilter()
            .long_term()
            .s_hist()
            .iter()
            .all(|&x| x == 0.0));
        assert_eq!(
            chain.postfilter().agc().gain(),
            crate::adaptive_gain_control::G_INIT
        );
    }

    /// §4.4: the concealed post-filter path decodes an erasure sentinel
    /// instead of erroring. After a good frame, an erased frame produces a
    /// finite 80-sample concealed output block (where the strict path
    /// returns `FrameDecodeError::Erased`).
    #[test]
    fn concealed_path_synthesizes_erased_frame() {
        let good = serial_bytes(&parity_ok_params());
        let erased = erased_serial_bytes();

        // Strict path errors on the sentinel.
        let mut strict = FrameDecoder::new();
        strict
            .decode_serial_frame_to_postfiltered(&good)
            .expect("good frame");
        assert_eq!(
            strict
                .decode_serial_frame_to_postfiltered(&erased)
                .unwrap_err(),
            FrameDecodeError::Erased,
        );

        // Concealed path synthesizes both.
        let mut dec = FrameDecoder::new();
        let g = dec
            .decode_serial_frame_to_postfiltered_concealed(&good)
            .expect("good frame conceals identically");
        assert!(g.output().iter().all(|s| s.is_finite()));
        let e = dec
            .decode_serial_frame_to_postfiltered_concealed(&erased)
            .expect("erased frame is concealed, not errored");
        let out = e.output();
        assert_eq!(out.len(), 2 * SUBFRAME_SIZE);
        assert!(out.iter().all(|s| s.is_finite()));
    }

    /// §4.4 leading erasure: an erasure before any good frame synthesizes
    /// silence (the clause-4.3 zero state), not an error.
    #[test]
    fn concealed_leading_erasure_is_silence() {
        let mut dec = FrameDecoder::new();
        let out = dec
            .decode_serial_frame_to_postfiltered_concealed(&erased_serial_bytes())
            .expect("leading erasure concealed")
            .output();
        // Zero excitation + identity LP → silence (within high-pass
        // transient bounds, which start at zero state → exactly zero).
        assert!(
            out.iter().all(|&s| s == 0.0),
            "leading erasure must be silent"
        );
    }

    /// §4.4 / §4.4.4: a sustained erasure run stays finite and bounded —
    /// the §4.4.2 gain attenuation (eqs (93)/(94)) compounds across the
    /// run so the energy decays rather than diverging.
    #[test]
    fn concealed_sustained_erasure_run_is_bounded() {
        let good = serial_bytes(&parity_ok_params());
        let erased = erased_serial_bytes();
        let mut dec = FrameDecoder::new();
        // Prime with a few good frames so a voicing class + last-good
        // frame are latched.
        for _ in 0..3 {
            dec.decode_serial_frame_to_postfiltered_concealed(&good)
                .expect("good frame");
        }
        // Long erasure run.
        let mut max_abs = 0.0f32;
        for _ in 0..20 {
            let out = dec
                .decode_serial_frame_to_postfiltered_concealed(&erased)
                .expect("concealed")
                .output();
            for &s in &out {
                assert!(s.is_finite());
                max_abs = max_abs.max(s.abs());
            }
        }
        // The attenuation keeps the run bounded well below any runaway.
        assert!(
            max_abs < 1.0e6,
            "erasure run diverged: max |sample| = {max_abs}"
        );
    }

    /// The concealed path is deterministic: two chains over the same
    /// good/erased sequence stay byte-identical (the eq (96) random
    /// generator is owned chain state).
    #[test]
    fn concealed_path_is_deterministic() {
        let good = serial_bytes(&parity_ok_params());
        let erased = erased_serial_bytes();
        let mut a = FrameDecoder::new();
        let mut b = FrameDecoder::new();
        for step in 0..8 {
            let f = if step % 3 == 2 { &erased } else { &good };
            let oa = a
                .decode_serial_frame_to_postfiltered_concealed(f)
                .unwrap()
                .output();
            let ob = b
                .decode_serial_frame_to_postfiltered_concealed(f)
                .unwrap()
                .output();
            assert_eq!(oa, ob, "step {step} diverged");
        }
    }

    /// A good frame through the concealed entry point produces identical
    /// output to the strict entry point (concealment is a superset — it
    /// only changes behaviour on erasures).
    #[test]
    fn concealed_good_frame_matches_strict_path() {
        let good = serial_bytes(&parity_ok_params());
        let mut strict = FrameDecoder::new();
        let mut conceal = FrameDecoder::new();
        for _ in 0..4 {
            let s = strict
                .decode_serial_frame_to_postfiltered(&good)
                .unwrap()
                .output();
            let c = conceal
                .decode_serial_frame_to_postfiltered_concealed(&good)
                .unwrap()
                .output();
            assert_eq!(s, c, "concealed path must match strict path on good frames");
        }
    }
}
