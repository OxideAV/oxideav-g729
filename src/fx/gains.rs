//! §3.9.1 / §4.1.5 fixed-point gain decode: the 4th-order MA gain
//! predictor on the Q10 dB grid, the eq (71) predicted gain through
//! the clause-5 [`crate::fx::dsp::log2`] / [`crate::fx::dsp::pow2`]
//! functions, and the eq (73)/(74) codebook reconstruction on the
//! integer grid (`ĝ_p` Q14, `ĝ_c` Q1).
//!
//! ## Schedule
//!
//! The prose pins the *values* (eqs (66)–(74), Ē = 30 dB, Table 9
//! init) but not the operator schedule. The schedule below is fitted
//! against the **reference excitation recovered from the staged
//! traces** (the staged `ŝ(n)` columns of consecutive sub-frames,
//! inverted through the eq (77) synthesis with the staged `â_i`,
//! give the reference `u(n)` exactly; at the pulse positions
//! `u(n) − ĝ_p·v(n)` brackets the reference `ĝ_c` to a fraction of a
//! Q1 step):
//!
//! 1. eq (66): `Σc²` over the Q13 codevector in a Word32 `L_mac`
//!    accumulator (Q27), then `log2`.
//! 2. `−10·log10(Σc²)` on Q14 as `mpy_32_16(e, f, −24660)`
//!    (`10·log10(2)` on Q13), plus the folded constant
//!    `Ē + 10·log10(40) + 27·10·log10(2) ≈ 127.3 dB` on Q14, shifted
//!    to Q24.
//! 3. eq (69): `L_mac` of the staged Q13 MA coefficients against the
//!    Q10 memory lands on the same Q24 grid; the predicted energy is
//!    **truncated to Q8 dB**.
//! 4. eq (71): the dB → base-2 exponent conversion multiplies by
//!    **`0.166` on Q15 (5439)** — the three-digit decimal, not the
//!    exact `log2(10)/20 = 0.16610` (a 0.07 % exponent difference
//!    that is level dependent: at a 55 dB prediction it moves `g′_c`
//!    by 0.42 %, which is the whole decode-side excitation mismatch
//!    the exact constant produced on the staged sub-frames). The Q16
//!    exponent splits into the integer / Q15-fraction pair `pow2`
//!    consumes; the mantissa is `pow2(14, frac)`.
//! 5. eq (74): the `γ̂` column sum is halved onto a Word16 (Q12),
//!    multiplied against the mantissa and **truncated** onto Q1.
//! 6. eq (72): `Û = 20·log10(γ̂)` truncated onto Q10.
//!
//! With this schedule the excitation of the nine staged consecutive
//! PITCH sub-frames reproduces sample-exactly except where `ĝ_c`
//! lands one Q1 step high (two of nine sub-frames) — the residual is
//! a sub-0.1 % truncation inside steps 3–5 that the staged material
//! does not pin (recorded docs ask: a trace of the gain-decode
//! intermediates).

use crate::fx::dsp::{l_extract, log2, mpy_32_16, pow2};
use crate::fx::ops::{
    add, extract_h, extract_l, l_add, l_deposit_l, l_mac, l_mult, l_shl, l_shr, mult, sature, sub,
};
use crate::tables::{GAIN_QUANT_MA_PREDICTOR_Q13, MA_NP, gain_ga_entry, gain_gb_entry};
use crate::tables::{GAIN_VQ_COL_GC, GAIN_VQ_COL_GP};

/// Subframe length the eq (66) energy averages over.
pub const L_SUBFR: usize = 40;

/// Clause-4.3 / Table 9 predictor start-up value: −14 dB on Q10.
pub const PAST_QUA_EN_INIT_Q10: i16 = -14336;

/// The eq (71) dB → base-2 exponent factor `0.166` on Q15 (the
/// three-digit decimal of `log2(10)/20`; see the module doc).
pub const DB_TO_LOG2_Q15: i16 = 5439;

/// `10·log10(2)` on the Q13 grid — converts a base-2 log to dB
/// (energy form).
pub const DB_PER_LOG2_ENERGY_Q13: i16 = 24_660;

/// `Ē + 10·log10(40) + 27·10·log10(2)` on Q14 — the folded eq (71)
/// constant for a Q27 energy accumulator (`= L_mult(32588, 32)`).
pub const ENERGY_CONST_Q14: i32 = 2_085_632;

/// `20·log10(2)` on the Q12 grid — converts a base-2 log to dB
/// (amplitude form).
const DB_PER_LOG2_Q12: i16 = 24_660;

/// Sweep hooks of the §3.9.1 schedule (black-box probes; the defaults
/// are the fitted schedule of the module doc).
#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct GainGridFx {
    /// Folded energy constant on Q14 dB.
    pub energy_const_q14: i32,
    /// dB → base-2 exponent factor on Q15.
    pub db_to_log2_q15: i16,
}

impl Default for GainGridFx {
    fn default() -> Self {
        Self {
            energy_const_q14: ENERGY_CONST_Q14,
            db_to_log2_q15: DB_TO_LOG2_Q15,
        }
    }
}

/// Fixed-point §3.9.1 gain predictor + §4.1.5 gain decoder state.
#[derive(Debug, Clone)]
pub struct GainDecoderFx {
    /// `[Û^(m−1) … Û^(m−4)]` on the Q10 dB grid.
    past_qua_en: [i16; MA_NP],
    /// γ̂ grid latitude (defaults to the uniform 2^13 base).
    grid: GainGridFx,
}

/// One subframe's decoded gains.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecodedGainsFx {
    /// Quantised adaptive-codebook gain `ĝ_p`, Q14.
    pub gain_pit_q14: i16,
    /// Quantised fixed-codebook gain `ĝ_c`, Q1.
    pub gain_code_q1: i16,
}

/// The eq (71) predicted gain `g'_c` in the form the eq (74)
/// reconstruction consumes: a `pow2` mantissa (`g'_c · 2^(13 − exp)`,
/// a Q14 value in `[16384, 32767]`) plus its integer exponent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GainPredictionFx {
    /// `g'_c · 2^(13 − exp)` as a Word32 (full-scale mantissa).
    pub g_prime_scaled: i32,
    /// The integer exponent of the eq (71) base-2 form.
    pub exp: i16,
}

impl Default for GainDecoderFx {
    fn default() -> Self {
        Self::new()
    }
}

impl GainDecoderFx {
    /// Fresh predictor in the clause-4.3 start-up state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            past_qua_en: [PAST_QUA_EN_INIT_Q10; MA_NP],
            grid: GainGridFx::default(),
        }
    }

    /// Override the γ̂ grid latitude (black-box sweep hook).
    #[doc(hidden)]
    pub fn set_grid(&mut self, grid: GainGridFx) {
        self.grid = grid;
    }

    /// Borrow the Q10 predictor memory (most recent first).
    #[must_use]
    pub fn past_qua_en(&self) -> &[i16; MA_NP] {
        &self.past_qua_en
    }

    /// eq (66)–(71): the predicted fixed-codebook gain for this
    /// subframe's Q13 codevector, as the `pow2` mantissa/exponent pair
    /// the eq (74) reconstruction consumes. Pure — the predictor
    /// memory is only advanced by [`Self::push`].
    #[must_use]
    pub fn predict(&self, code_q13: &[i16; L_SUBFR]) -> GainPredictionFx {
        // eq (66) energy: Σ code² on Q27 (Q13 × Q13 doubled).
        let mut l_ener = 0i32;
        for &c in code_q13.iter() {
            l_ener = l_mac(l_ener, c, c);
        }
        let (e, f) = log2(l_ener);
        // −10·log10(Σc²) on Q14, plus the folded constant, to Q24.
        let mut acc = mpy_32_16(e, f, -DB_PER_LOG2_ENERGY_Q13);
        acc = l_add(acc, self.grid.energy_const_q14);
        acc = l_shl(acc, 10);
        // eq (69): Σ b_i·Û^(m−i), Q13 × Q10 doubled = Q24.
        for (b, u) in GAIN_QUANT_MA_PREDICTOR_Q13
            .iter()
            .take(MA_NP)
            .zip(self.past_qua_en.iter())
        {
            acc = l_mac(acc, *b, *u);
        }
        // Predicted gain in dB, truncated to Q8.
        let db_q8 = extract_h(acc);
        // eq (71): base-2 exponent on Q16 → integer + Q15 fraction.
        let x_q16 = l_shr(l_mult(db_q8, self.grid.db_to_log2_q15), 8);
        let (x_int, x_frac) = l_extract(x_q16);
        GainPredictionFx {
            g_prime_scaled: pow2(14, x_frac),
            exp: sub(x_int, 1),
        }
    }

    /// eqs (73)/(74): the quantised gains for the **codebook-domain**
    /// `(ga, gb)` pair given the subframe's prediction — `ĝ_p` on Q14,
    /// `ĝ_c = γ̂ · g'_c` landed on Q1. Pure (no memory advance), so the
    /// encoder's §3.9.2 search can score every candidate pair on the
    /// decoder's exact grid.
    ///
    /// # Panics
    ///
    /// Panics if `ga`/`gb` exceed the codebook domains (3-/4-bit
    /// fields cannot).
    #[must_use]
    pub fn reconstruct(&self, pred: &GainPredictionFx, ga: usize, gb: usize) -> DecodedGainsFx {
        let ga_row = gain_ga_entry(ga);
        let gb_row = gain_gb_entry(gb);

        // eq (73): ĝ_p on Q14 (the column extremes sum well inside
        // Word16).
        let gain_pit_q14 = add(ga_row[GAIN_VQ_COL_GP], gb_row[GAIN_VQ_COL_GP]);
        // eq (74): γ̂ on the 2^13 grid — accumulated on Word32: the
        // worst-case column sum (27162 + 14276 = 41438, γ̂ ≈ 5.06)
        // overflows Word16, so the reconstruction keeps the two stage
        // contributions apart until they are combined against g'_c.
        // eq (74): γ̂ column sum (2^13 grid, up to 41438) halved onto
        // a Word16, times the Q14 mantissa: with
        // g′_c = mantissa·2^(exp − 13) the doubled product lands on Q1
        // after a right shift of `25 − exp`, truncating.
        let gamma_half = extract_l(l_shr(
            l_add(
                l_deposit_l(ga_row[GAIN_VQ_COL_GC]),
                l_deposit_l(gb_row[GAIN_VQ_COL_GC]),
            ),
            1,
        ));
        let mantissa = sature(pred.g_prime_scaled);
        let prod = l_mult(gamma_half, mantissa);
        let gain_code_q1 = sature(l_shr(prod, sub(25, pred.exp))).max(0);

        DecodedGainsFx {
            gain_pit_q14,
            gain_code_q1,
        }
    }

    /// eq (72) memory advance for the chosen `(ga, gb)` pair:
    /// `Û^(m) = 20·log10(γ̂)` pushed onto the Q10 history.
    pub fn push(&mut self, ga: usize, gb: usize) {
        let ga_row = gain_ga_entry(ga);
        let gb_row = gain_gb_entry(gb);
        let gamma_q13_32 = l_add(
            l_deposit_l(ga_row[GAIN_VQ_COL_GC]),
            l_deposit_l(gb_row[GAIN_VQ_COL_GC]),
        );

        // eq (72) decode form: Û^(m) = 20·log10(γ̂) pushed on Q10,
        // truncated. log2(γ̂) = log2(gamma_q13_32) − 13; dB = 6.0206·log2.
        let (ge, gf) = log2(gamma_q13_32);
        let lg_q15 = l_add(l_shl(l_deposit_l(sub(ge, 13)), 15), l_deposit_l(gf));
        let (lg_hi, lg_lo) = l_extract(l_shl(lg_q15, 9)); // Q24
        let u_q21 = mpy_32_16(lg_hi, lg_lo, DB_PER_LOG2_Q12); // Q24·Q12 → Q21
        let u_q10 = extract_h(l_shl(u_q21, 5)); // Q21 → Q26 high word = Q10

        for k in (1..MA_NP).rev() {
            self.past_qua_en[k] = self.past_qua_en[k - 1];
        }
        self.past_qua_en[0] = u_q10;
    }

    /// Decodes one subframe's gains from the **codebook-domain**
    /// `(ga, gb)` indices and the Q13 codevector, advancing the
    /// predictor memory (eq (72) decode form) — [`Self::predict`] →
    /// [`Self::reconstruct`] → [`Self::push`].
    ///
    /// # Panics
    ///
    /// Panics if `ga`/`gb` exceed the codebook domains (3-/4-bit
    /// fields cannot).
    pub fn decode(&mut self, ga: usize, gb: usize, code_q13: &[i16; L_SUBFR]) -> DecodedGainsFx {
        let pred = self.predict(code_q13);
        let gains = self.reconstruct(&pred, ga, gb);
        self.push(ga, gb);
        gains
    }

    /// §4.4.2 erased-frame gain attenuation (eqs (93)/(94)) plus the
    /// §4.4.3 4 dB predictor-memory attenuation with the −14 dB floor:
    /// `ĝ_p ← 0.9·ĝ_p` (bounded by the caller's previous value),
    /// `ĝ_c ← 0.98·ĝ_c`, and every memory slot moves toward the floor.
    pub fn conceal(&mut self, prev: DecodedGainsFx) -> DecodedGainsFx {
        // eq (93): fixed-codebook gain × 0.98 (Q15 32113).
        let gain_code_q1 = mult(prev.gain_code_q1, 32113);
        // eq (94): adaptive gain × 0.9 (Q15 29491).
        let gain_pit_q14 = mult(prev.gain_pit_q14, 29491);

        // §4.4.3: mean of the memory minus 4 dB, floored at −14 dB —
        // the spec prose says the predictor memory is "updated with"
        // an attenuated average; the float module shifts every slot,
        // and this path mirrors it: each slot loses 4 dB (Q10 4096)
        // with the −14 dB floor.
        for slot in &mut self.past_qua_en {
            let dropped = sub(*slot, 4096);
            *slot = if dropped < PAST_QUA_EN_INIT_Q10 {
                PAST_QUA_EN_INIT_Q10
            } else {
                dropped
            };
        }

        DecodedGainsFx {
            gain_pit_q14,
            gain_code_q1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A four-pulse ±8192 codevector (the minimum-energy §3.8 shape)
    /// from the start-up state must reproduce the float predictor's
    /// ĝ_c within a fraction of a percent (both sides implement the
    /// same eqs (66)–(74); only sub-LSB grids differ).
    #[test]
    fn matches_float_predictor_from_startup() {
        let mut code = [0i16; L_SUBFR];
        code[3] = 8192;
        code[11] = -8192;
        code[22] = 8192;
        code[38] = -8192;

        for (ga, gb) in [(0usize, 0usize), (3, 7), (5, 11), (7, 15), (2, 4)] {
            let mut fx = GainDecoderFx::new();
            let got = fx.decode(ga, gb, &code);

            // Float oracle.
            let q = crate::gain_reconstruct::reconstruct_gains(ga, gb).unwrap();
            let mut fl = crate::gain_predict::GainPredictor::new();
            let cf: [f32; L_SUBFR] = std::array::from_fn(|i| f32::from(code[i]) / 8192.0);
            let (gc_hat, _) = fl.predict_and_update(&cf, q.gamma_hat);

            let got_gc = f32::from(got.gain_code_q1) / 2.0;
            // Tolerance: 1% relative plus the Q1 grid step (0.5) —
            // ĝ_c lands on the half-integer grid by truncation.
            let tol = 0.01 * gc_hat.abs() + 0.5;
            assert!(
                (got_gc - gc_hat).abs() <= tol,
                "(ga={ga}, gb={gb}) fx ĝ_c {got_gc} vs float {gc_hat}"
            );
            let got_gp = f32::from(got.gain_pit_q14) / 16384.0;
            assert!((got_gp - q.g_p_hat).abs() < 1e-4);
        }
    }

    /// The eq (72) memory push tracks the float predictor's dB history
    /// within a Q10 LSB across a run of subframes.
    #[test]
    fn memory_tracks_float_history() {
        let mut code = [0i16; L_SUBFR];
        code[0] = 8192;
        code[10] = 8192;
        code[20] = -8192;
        code[30] = -8192;

        let mut fx = GainDecoderFx::new();
        let mut fl = crate::gain_predict::GainPredictor::new();
        let cf: [f32; L_SUBFR] = std::array::from_fn(|i| f32::from(code[i]) / 8192.0);

        for (ga, gb) in [(1usize, 2usize), (6, 13), (4, 9), (0, 15), (7, 0), (3, 3)] {
            let q = crate::gain_reconstruct::reconstruct_gains(ga, gb).unwrap();
            let _ = fx.decode(ga, gb, &code);
            let _ = fl.predict_and_update(&cf, q.gamma_hat);
            let got_db = f32::from(fx.past_qua_en()[0]) / 1024.0;
            let want_db = fl.history_db()[0];
            assert!(
                (got_db - want_db).abs() < 0.02,
                "(ga={ga}, gb={gb}) Û fx {got_db} vs float {want_db}"
            );
        }
    }

    /// Concealment attenuates both gains and drains the memory toward
    /// the −14 dB floor.
    #[test]
    fn conceal_attenuates() {
        let mut fx = GainDecoderFx::new();
        let prev = DecodedGainsFx {
            gain_pit_q14: 16000,
            gain_code_q1: 2000,
        };
        let g1 = fx.conceal(prev);
        assert!(g1.gain_pit_q14 < prev.gain_pit_q14);
        assert!(g1.gain_code_q1 < prev.gain_code_q1);
        for _ in 0..10 {
            let _ = fx.conceal(g1);
        }
        for &s in fx.past_qua_en() {
            assert_eq!(s, PAST_QUA_EN_INIT_Q10);
        }
    }
}
