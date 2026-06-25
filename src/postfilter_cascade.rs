//! §4.2 decoder **post-processing cascade** — the five-stage adaptive
//! postfilter that turns the §4.1.6 reconstructed speech `ŝ(n)` into the
//! final decoder output `sf′(n)` (×2-upscaled, output-high-passed PCM).
//!
//! ## Spec source — clause 4.2 (06/2012 Recommendation)
//!
//! Clause 4.2 opens: "Postprocessing consists of three functions:
//! adaptive postfiltering, high-pass filtering and signal up-scaling. The
//! adaptive postfilter is the cascade of three filters: a long-term
//! postfilter `H_p(z)`, a short-term postfilter `H_f(z)` and a tilt
//! compensation filter `H_t(z)`, followed by an adaptive gain control
//! procedure." The cascade order (clause 4.2, fig. / eq references):
//!
//! 1. **§4.2.1 long-term postfilter `H_p(z)`** — harmonic postfilter on
//!    `ŝ(n)`, anchored on the **first subframe's** transmitted integer
//!    pitch delay `int(T_1)` (clause 4.2.1: "int(T_1) is the integer part
//!    of the (transmitted) pitch delay T_1 in the first subframe"); the
//!    same `int(T_1)` anchors the integer search in **both** subframes.
//!    Two-pass form (integer anchor + 1/8-resolution fractional
//!    refinement) — see [`crate::long_term_postfilter`].
//! 2. **§4.2.2 short-term postfilter `H_f(z)`** — `Â(z/γ_n)/Â(z/γ_d)`
//!    with gain normalisation `1/g_f` ([`crate::short_term_postfilter`]).
//! 3. **§4.2.3 tilt compensation `H_t(z)`** — first-order tilt corrector
//!    using the short-term postfilter impulse response
//!    ([`crate::tilt_compensation`]).
//! 4. **§4.2.4 adaptive gain control** — eq (88)–(90), matching the
//!    cascade output energy back to `ŝ(n)`
//!    ([`crate::adaptive_gain_control`]).
//! 5. **§4.2.5 output high-pass `H_h2(z)` + ×2 up-scaling** — eq (91)
//!    ([`crate::post_process`]).
//!
//! The §4.2.4 gain control compares the energy of the §4.2.3 tilt output
//! against the energy of the **reconstructed speech `ŝ(n)`** (clause
//! 4.2.4: "the gain control unit … compares the energies of the signal
//! `ŝ(n)` and the postfiltered signal `sf(n)`"), so this cascade keeps
//! the original `ŝ(n)` subframe as the AGC reference, not any
//! intermediate stage.
//!
//! ## State (clause 4.3 init)
//!
//! Every stage owns its own clause-4.3 start-up state (all-zero except
//! the AGC's `g(−1) = 1.0`, Table 9). This cascade simply holds the five
//! stage objects and threads them in spec order; it adds no arithmetic of
//! its own beyond the per-subframe sequencing.

use crate::adaptive_gain_control::AdaptiveGainControl;
use crate::decode_chain::DecodedFrame;
use crate::fixed_codebook::SUBFRAME_SIZE;
use crate::long_term_postfilter::{LongTermDecision, LongTermPostfilter};
use crate::lp_synthesis::SynthesizedFrame;
use crate::post_process::OutputHighPass;
use crate::short_term_postfilter::ShortTermPostfilter;
use crate::tables::M;
use crate::tilt_compensation::TiltCompensation;

/// Number of 5 ms subframes in a 10 ms frame.
const SUBFRAMES_PER_FRAME: usize = 2;

/// The §4.2 post-processing products for **one subframe**, kept for
/// inspection / tests alongside the final PCM.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PostfilteredSubframe {
    /// §4.2.1 long-term-postfilter decision (chosen integer anchor +
    /// 1/8 fractional offset + gain `g_l`) for this subframe.
    pub long_term: LongTermDecision,
    /// §4.2.5 final output samples `sf′(n)` (×2-upscaled, output
    /// high-passed), `n = 0 … 39`.
    pub output: [f32; SUBFRAME_SIZE],
}

/// The §4.2 post-processing products for **one 10 ms frame** — the two
/// 40-sample output subframes plus their per-subframe decisions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PostfilteredFrame {
    /// Per-subframe post-processing products in spec §4.1 order.
    pub subframes: [PostfilteredSubframe; SUBFRAMES_PER_FRAME],
}

impl PostfilteredFrame {
    /// The 80 final output samples of the frame in time order (subframe 1
    /// then subframe 2). These are the §4.2.5 `sf′(n)`: the fully
    /// post-filtered, ×2-upscaled decoder output.
    #[must_use]
    pub fn output(&self) -> [f32; SUBFRAMES_PER_FRAME * SUBFRAME_SIZE] {
        let mut out = [0.0f32; SUBFRAMES_PER_FRAME * SUBFRAME_SIZE];
        out[..SUBFRAME_SIZE].copy_from_slice(&self.subframes[0].output);
        out[SUBFRAME_SIZE..].copy_from_slice(&self.subframes[1].output);
        out
    }
}

/// Stateful §4.2 post-processing cascade.
///
/// Owns the five cascade stages, each carrying its own clause-4.3
/// cross-subframe state. Call [`Self::process_frame`] once per decoded
/// frame, in stream order, handing it the §4.1.6 [`SynthesizedFrame`] and
/// the originating §4.1 [`DecodedFrame`] (used for the per-subframe LP
/// coefficients and the first-subframe integer pitch delay).
#[derive(Debug, Clone, Default)]
pub struct PostfilterCascade {
    /// §4.2.1 long-term (harmonic) postfilter `H_p(z)`.
    long_term: LongTermPostfilter,
    /// §4.2.2 short-term postfilter `H_f(z)`.
    short_term: ShortTermPostfilter,
    /// §4.2.3 tilt-compensation filter `H_t(z)`.
    tilt: TiltCompensation,
    /// §4.2.4 adaptive gain control.
    agc: AdaptiveGainControl,
    /// §4.2.5 output high-pass `H_h2(z)` + ×2 up-scaling.
    high_pass: OutputHighPass,
}

impl PostfilterCascade {
    /// Build the cascade with every stage at its clause-4.3 start-up
    /// state.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Borrow the §4.2.1 long-term postfilter stage for inspection /
    /// tests.
    #[must_use]
    pub fn long_term(&self) -> &LongTermPostfilter {
        &self.long_term
    }

    /// Borrow the §4.2.4 adaptive-gain-control stage for inspection /
    /// tests.
    #[must_use]
    pub fn agc(&self) -> &AdaptiveGainControl {
        &self.agc
    }

    /// Run the §4.2 cascade on one decoded frame, advancing every stage's
    /// cross-subframe state.
    ///
    /// `synth` is the §4.1.6 reconstructed speech for the frame; `frame`
    /// is the originating §4.1 parameter decode (its per-subframe `lp`
    /// drives `H_f` / `H_t`, and the **first** subframe's `int(T)` anchors
    /// the §4.2.1 integer-delay search for both subframes per clause
    /// 4.2.1).
    #[must_use]
    pub fn process_frame(
        &mut self,
        synth: &SynthesizedFrame,
        frame: &DecodedFrame,
    ) -> PostfilteredFrame {
        // Clause 4.2.1: the integer-delay search is anchored on the
        // integer part of the *first-subframe* transmitted delay T_1,
        // used for both subframes of the frame.
        let int_t1 = frame.subframes[0].pitch.int_t.max(0) as usize;

        let mut subs = [None, None];
        for (i, slot) in subs.iter_mut().enumerate() {
            let s_hat: &[f32; SUBFRAME_SIZE] = &synth.subframes[i].speech;
            let a: &[f32; M] = &frame.subframes[i].lp;

            // §4.2.1 H_p(z): long-term (harmonic) postfilter on ŝ(n).
            let (hp, decision) = self.long_term.filter_subframe(s_hat, a, int_t1);
            // §4.2.2 H_f(z): short-term postfilter.
            let hf = self.short_term.filter_subframe(&hp, a);
            // §4.2.3 H_t(z): tilt compensation (coefficients from a).
            let ht = self.tilt.filter_subframe(&hf, a);
            // §4.2.4 adaptive gain control: match ht's energy to ŝ(n).
            let agc = self.agc.process_subframe(s_hat, &ht);
            // §4.2.5 output high-pass H_h2(z) + ×2 up-scaling.
            let mut output = agc;
            self.high_pass.filter_in_place(&mut output);

            *slot = Some(PostfilteredSubframe {
                long_term: decision,
                output,
            });
        }

        let [s0, s1] = subs;
        PostfilteredFrame {
            subframes: [
                s0.expect("subframe 0 populated above"),
                s1.expect("subframe 1 populated above"),
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adaptive_gain_control::G_INIT;
    use crate::decode_chain::FrameDecoder;
    use crate::parameters::Parameters;

    /// A frame whose pitch parity holds under the corpus-validated
    /// odd-parity rule (mirrors `decode_chain`'s fixture).
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

    /// A fresh cascade is all-zero state except the AGC's Table-9 unit
    /// gain.
    #[test]
    fn new_cascade_has_clause43_startup_state() {
        let c = PostfilterCascade::new();
        assert!(c.long_term().s_hist().iter().all(|&x| x == 0.0));
        assert!(c.long_term().r_hist().iter().all(|&x| x == 0.0));
        assert_eq!(c.agc().gain(), G_INIT);
    }

    /// Decode a real frame to speech, then run the cascade: the 80 output
    /// samples are finite and the cascade advances off its start-up state.
    #[test]
    fn cascade_produces_finite_output_and_advances_state() {
        let params = parity_ok_params();
        let mut chain = FrameDecoder::new();
        let frame = chain.decode_parameters(&params).expect("in-domain");
        let mut synth = crate::lp_synthesis::Synthesizer::new();
        let speech = synth.synthesize_frame(&frame);

        let mut cascade = PostfilterCascade::new();
        let out = cascade.process_frame(&speech, &frame);
        let pcm = out.output();
        assert_eq!(pcm.len(), 80);
        assert!(pcm.iter().all(|s| s.is_finite()));
        // The long-term postfilter history is no longer all-zero.
        assert!(cascade.long_term().s_hist().iter().any(|&x| x != 0.0));
    }

    /// With an all-zero (flat) short-term filter and an all-zero input
    /// `ŝ(n)`, every cascade stage is a linear filter of zero, so the
    /// output is identically zero — a sanity check that no stage injects
    /// a constant.
    #[test]
    fn silent_input_yields_silent_output() {
        // Build a DecodedFrame-shaped fixture with flat LP and zero speech.
        let params = parity_ok_params();
        let mut chain = FrameDecoder::new();
        let mut frame = chain.decode_parameters(&params).expect("in-domain");
        for sub in &mut frame.subframes {
            sub.lp = [0.0; M];
        }
        let speech = SynthesizedFrame {
            subframes: [
                crate::lp_synthesis::SynthesizedSubframe {
                    adaptive: [0.0; SUBFRAME_SIZE],
                    excitation: [0.0; SUBFRAME_SIZE],
                    speech: [0.0; SUBFRAME_SIZE],
                },
                crate::lp_synthesis::SynthesizedSubframe {
                    adaptive: [0.0; SUBFRAME_SIZE],
                    excitation: [0.0; SUBFRAME_SIZE],
                    speech: [0.0; SUBFRAME_SIZE],
                },
            ],
        };
        let mut cascade = PostfilterCascade::new();
        let out = cascade.process_frame(&speech, &frame);
        assert!(out.output().iter().all(|&s| s == 0.0));
    }

    /// The per-subframe long-term decision is reported and its gain is in
    /// the bounded `[0, 1]` range of eq (83).
    #[test]
    fn long_term_decision_gain_is_bounded() {
        let params = parity_ok_params();
        let mut chain = FrameDecoder::new();
        let frame = chain.decode_parameters(&params).expect("in-domain");
        let mut synth = crate::lp_synthesis::Synthesizer::new();
        let speech = synth.synthesize_frame(&frame);
        let mut cascade = PostfilterCascade::new();
        let out = cascade.process_frame(&speech, &frame);
        for sub in &out.subframes {
            assert!((0.0..=1.0).contains(&sub.long_term.gain));
        }
    }
}
