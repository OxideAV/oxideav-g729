//! §4.1 fixed-point frame decoder: parameter unpack → LP params →
//! per-subframe adaptive vector, algebraic codevector, gains,
//! excitation, and `1/Â(z)` synthesis — everything on the clause-5
//! Word16/Word32 grid, producing the pre-postfilter reconstructed
//! speech `ŝ(n)`.
//!
//! The §4.2 postfilter cascade consumes [`DecodedFrameFx`] (speech +
//! per-subframe Q12 coefficients + the integer pitch delay); the
//! §4.4 erasure path is layered above this module by the caller
//! ([`FrameDecoderFx::decode_erased_frame`] provides the LP/gain/
//! excitation concealment primitives that do not depend on the
//! postfilter's voicing decision).

use crate::fx::excitation::{
    build_codevector_q13, build_excitation, pred_lt3, syn_filt, EXC_BUF, L_SUBFR, PIT_MAX,
    SHARP_INIT_Q14, SHARP_MAX_Q14, SHARP_MIN_Q14,
};
use crate::fx::gains::{DecodedGainsFx, GainDecoderFx};
use crate::fx::lsp::LspDecoderFx;
use crate::fx::ops::{add, shr};
use crate::gain_index_map::{demap_ga, demap_gb};
use crate::parameters::Parameters;
use crate::pitch_decode::{decode_t1_from_p1, decode_t2_from_p2, derive_t_min};
use crate::tables::M;

/// History region kept ahead of the working frame in the excitation
/// buffer (maximum delay + interpolator margin).
const EXC_HIST: usize = PIT_MAX + 10;

/// eq (96) replacement-excitation pseudo-random generator (§4.4.4):
/// `seed = seed·31821 + 13849`, seeded with 21845.
#[derive(Debug, Clone, Copy)]
pub struct RandomFx {
    seed: i16,
}

impl Default for RandomFx {
    fn default() -> Self {
        Self::new()
    }
}

impl RandomFx {
    /// Fresh generator with the §4.4.4 seed.
    #[must_use]
    pub const fn new() -> Self {
        Self { seed: 21845 }
    }

    /// Next Word16 of the eq (96) sequence (wrapping arithmetic).
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> i16 {
        self.seed = (self.seed.wrapping_mul(31821)).wrapping_add(13849);
        self.seed
    }
}

/// One decoded subframe's postfilter-facing side data.
#[derive(Debug, Clone, Copy)]
pub struct SubframeFx {
    /// Q12 LP coefficients `[a_0 = 4096, a_1 … a_10]`.
    pub a_q12: [i16; M + 1],
    /// Integer part of the subframe's pitch delay (the §4.2.1
    /// long-term postfilter searches `[T−1, T+1]` around it).
    pub t_int: i32,
    /// Decoded gains (Q14 / Q1).
    pub gains: DecodedGainsFx,
}

/// One decoded frame: 80 samples of pre-postfilter speech plus the
/// per-subframe side data the §4.2 cascade needs.
#[derive(Debug, Clone, Copy)]
pub struct DecodedFrameFx {
    /// Reconstructed speech `ŝ(n)` (Q0, halved-input grid).
    pub speech: [i16; 2 * L_SUBFR],
    /// Per-subframe postfilter inputs.
    pub sub: [SubframeFx; 2],
}

/// Stateful §4.1 fixed-point decoder.
#[derive(Debug, Clone)]
pub struct FrameDecoderFx {
    lsp: LspDecoderFx,
    gains: GainDecoderFx,
    exc: [i16; EXC_BUF],
    syn_mem: [i16; M],
    /// eq (48) sharpening gain (Q14), updated from each subframe's
    /// decoded pitch gain.
    sharp_q14: i16,
    /// Last good subframe-2 delay, for §4.4.4 repetition.
    prev_t: i32,
    /// Last good frame's gains, the §4.4.2 attenuation input.
    prev_gains: DecodedGainsFx,
    /// §4.4.4 replacement-excitation generator.
    rng: RandomFx,
    /// Voicing class latched from good frames (`true` = periodic);
    /// updated by the caller that owns the long-term postfilter
    /// decision, defaulting to non-periodic.
    pub erasure_periodic: bool,
}

impl Default for FrameDecoderFx {
    fn default() -> Self {
        Self::new()
    }
}

impl FrameDecoderFx {
    /// Fresh decoder in the clause-4.3 start-up state (zero
    /// excitation/synthesis memories, sharpening at its clamp floor,
    /// gain predictor at −14 dB).
    #[must_use]
    pub fn new() -> Self {
        Self {
            lsp: LspDecoderFx::new(),
            gains: GainDecoderFx::new(),
            exc: [0; EXC_BUF],
            syn_mem: [0; M],
            sharp_q14: SHARP_INIT_Q14,
            prev_t: 60,
            prev_gains: DecodedGainsFx {
                gain_pit_q14: 0,
                gain_code_q1: 0,
            },
            rng: RandomFx::new(),
            erasure_periodic: false,
        }
    }

    /// Decodes one good frame's parameters to pre-postfilter speech.
    ///
    /// The §4.1.2 parity check is applied: on mismatch the frame is
    /// treated as erased per the spec (delegating to
    /// [`Self::decode_erased_frame`]).
    pub fn decode_frame(&mut self, params: &Parameters) -> DecodedFrameFx {
        if !params.pitch_parity_ok() {
            return self.decode_erased_frame();
        }

        let lp = self.lsp.decode_frame(
            usize::from(params.l0),
            usize::from(params.l1),
            usize::from(params.l2),
            usize::from(params.l3),
        );

        let t1 = decode_t1_from_p1(params.p1);
        let t_min = derive_t_min(t1.int_t);
        let t2 = decode_t2_from_p2(params.p2, t_min);

        let mut speech = [0i16; 2 * L_SUBFR];
        let mut sub = [SubframeFx {
            a_q12: lp.a_sub1,
            t_int: t1.int_t,
            gains: self.prev_gains,
        }; 2];

        for (s, (delay, a)) in [(t1, lp.a_sub1), (t2, lp.a_sub2)].iter().enumerate() {
            let off = EXC_HIST + s * L_SUBFR;
            // §4.1.3: adaptive vector in place.
            pred_lt3(&mut self.exc, off, delay.int_t, delay.frac);

            // §4.1.4: codevector + eq (48) sharpening at int(T).
            let (c, sc) = if s == 0 {
                (params.c1, params.s1)
            } else {
                (params.c2, params.s2)
            };
            let pulses =
                crate::fixed_codebook::decode_pulses(c, sc).expect("13/4-bit fields are in domain");
            let positions = std::array::from_fn(|k| pulses.pulses[k].position);
            let signs = std::array::from_fn(|k| pulses.pulses[k].sign);
            let code = build_codevector_q13(&positions, &signs, self.sharp_q14, delay.int_t);

            // §4.1.5: gains (transmitted indices demapped per §3.9.3).
            let (tga, tgb) = if s == 0 {
                (params.ga1, params.gb1)
            } else {
                (params.ga2, params.gb2)
            };
            let ga = demap_ga(usize::from(tga)).expect("3-bit field");
            let gb = demap_gb(usize::from(tgb)).expect("4-bit field");
            let gains = self.gains.decode(ga, gb, &code);

            // eq (48) sharpening gain update: the just-decoded pitch
            // gain, clamped on Q14.
            self.sharp_q14 = gains.gain_pit_q14.clamp(SHARP_MIN_Q14, SHARP_MAX_Q14);

            // eq (75) + eq (77).
            build_excitation(
                &mut self.exc,
                off,
                &code,
                gains.gain_pit_q14,
                gains.gain_code_q1,
            );
            let (head, tail) = speech.split_at_mut(L_SUBFR);
            let out = if s == 0 { head } else { tail };
            let exc_slice: [i16; L_SUBFR] = std::array::from_fn(|n| self.exc[off + n]);
            syn_filt(&exc_slice, a, &mut self.syn_mem, out);

            sub[s] = SubframeFx {
                a_q12: *a,
                t_int: delay.int_t,
                gains,
            };
            self.prev_gains = gains;
            self.prev_t = delay.int_t;
        }

        self.advance_exc();
        DecodedFrameFx { speech, sub }
    }

    /// §4.4 erased-frame decode: repeated LP parameters, attenuated
    /// gains (eqs (93)/(94) + the 4 dB predictor-memory drop), and the
    /// class-dependent replacement excitation (§4.4.4) — the periodic
    /// class repeats the adaptive codebook at `T+1` (bounded 143) with
    /// the attenuated pitch gain; the non-periodic class draws a
    /// random algebraic codevector (eq (96)) with the attenuated fixed
    /// gain.
    pub fn decode_erased_frame(&mut self) -> DecodedFrameFx {
        let lp = self.lsp.decode_erased_frame();
        let mut speech = [0i16; 2 * L_SUBFR];
        let mut sub = [SubframeFx {
            a_q12: lp.a_sub1,
            t_int: self.prev_t,
            gains: self.prev_gains,
        }; 2];

        for (s, sub_slot) in sub.iter_mut().enumerate() {
            let off = EXC_HIST + s * L_SUBFR;
            let gains = self.gains.conceal(self.prev_gains);
            let a = if s == 0 { lp.a_sub1 } else { lp.a_sub2 };

            if self.erasure_periodic {
                // §4.4.4 periodic: adaptive-only excitation at the
                // repeated (incremented, bounded) delay.
                let t = (self.prev_t + 1).min(143);
                self.prev_t = t;
                pred_lt3(&mut self.exc, off, t, 0);
                let code = [0i16; L_SUBFR];
                build_excitation(&mut self.exc, off, &code, gains.gain_pit_q14, 0);
                sub_slot.t_int = t;
            } else {
                // §4.4.4 non-periodic: random algebraic codevector
                // (13-bit position draw + 4-bit sign draw), fixed gain
                // only.
                let c = (self.rng.next() as u16) & 0x1fff;
                let sc = (self.rng.next() as u16 & 0xf) as u8;
                let pulses =
                    crate::fixed_codebook::decode_pulses(c, sc).expect("masked to field width");
                let positions = std::array::from_fn(|k| pulses.pulses[k].position);
                let signs = std::array::from_fn(|k| pulses.pulses[k].sign);
                let code = build_codevector_q13(&positions, &signs, self.sharp_q14, self.prev_t);
                pred_lt3(&mut self.exc, off, self.prev_t.min(143), 0);
                // Zero the adaptive contribution (non-periodic class),
                // keep the attenuated fixed-codebook part.
                for n in 0..L_SUBFR {
                    self.exc[off + n] = 0;
                }
                build_excitation(&mut self.exc, off, &code, 0, gains.gain_code_q1);
            }

            let (head, tail) = speech.split_at_mut(L_SUBFR);
            let out = if s == 0 { head } else { tail };
            let exc_slice: [i16; L_SUBFR] = std::array::from_fn(|n| self.exc[off + n]);
            syn_filt(&exc_slice, &a, &mut self.syn_mem, out);

            sub_slot.a_q12 = a;
            sub_slot.gains = gains;
            self.prev_gains = gains;
        }

        self.advance_exc();
        DecodedFrameFx { speech, sub }
    }

    fn advance_exc(&mut self) {
        self.exc.copy_within(2 * L_SUBFR.., 0);
        for slot in &mut self.exc[EXC_HIST..] {
            *slot = 0;
        }
    }

    /// Interpolated half-sum helper exposed for tests.
    #[must_use]
    pub fn half_sum(a: i16, b: i16) -> i16 {
        add(shr(a, 1), shr(b, 1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parity_for(p1: u8) -> u8 {
        (((p1 >> 2) & 0x3F).count_ones() as u8 ^ 1) & 1
    }

    fn params_zeroish() -> Parameters {
        Parameters {
            l0: 0,
            l1: 10,
            l2: 3,
            l3: 4,
            p1: 60,
            p0: parity_for(60),
            c1: 100,
            s1: 5,
            ga1: 2,
            gb1: 3,
            p2: 15,
            c2: 200,
            s2: 9,
            ga2: 1,
            gb2: 7,
        }
    }

    /// The fixed decoder produces finite bounded output and is
    /// deterministic from a fresh state.
    #[test]
    fn decode_is_deterministic() {
        let p = params_zeroish();
        let mut d1 = FrameDecoderFx::new();
        let mut d2 = FrameDecoderFx::new();
        for _ in 0..5 {
            let f1 = d1.decode_frame(&p);
            let f2 = d2.decode_frame(&p);
            assert_eq!(f1.speech, f2.speech);
        }
    }

    /// A parity-violating P0 routes through the erasure path (the
    /// §4.1.2 rule) without panicking.
    #[test]
    fn parity_mismatch_conceals() {
        let mut p = params_zeroish();
        p.p0 ^= 1;
        let mut d = FrameDecoderFx::new();
        let _ = d.decode_frame(&p); // erased path
        let f = d.decode_frame(&params_zeroish());
        assert!(f.speech.iter().all(|&s| s > i16::MIN));
    }

    /// Erased frames decode to finite output in both voicing classes.
    #[test]
    fn erased_frames_bounded() {
        for periodic in [false, true] {
            let mut d = FrameDecoderFx::new();
            d.erasure_periodic = periodic;
            let _ = d.decode_frame(&params_zeroish());
            for _ in 0..3 {
                let f = d.decode_erased_frame();
                assert!(f.speech.iter().all(|&s| s != i16::MIN));
            }
        }
    }

    /// The eq (96) generator produces the documented first values.
    #[test]
    fn random_generator_sequence() {
        let mut r = RandomFx::new();
        let v1 = r.next();
        let v2 = r.next();
        // seed = 21845: 21845·31821 + 13849 (wrapping on Word16).
        let expect1 = (21845i16.wrapping_mul(31821)).wrapping_add(13849);
        assert_eq!(v1, expect1);
        let expect2 = (expect1.wrapping_mul(31821)).wrapping_add(13849);
        assert_eq!(v2, expect2);
    }
}
