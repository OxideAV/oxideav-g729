//! §3.9.1 / §4.1.5 fixed-point gain decode: the 4th-order MA gain
//! predictor on the Q10 log grid, the eq (71) predicted gain through
//! the clause-5 [`crate::fx::dsp::log2`] / [`crate::fx::dsp::pow2`]
//! functions, and the eq (73)/(74) codebook reconstruction on the
//! integer grid (`ĝ_p` Q14, `ĝ_c` Q1).
//!
//! ## Number grids
//!
//! - **Predictor memory** `Û^(m−k)` — Q10 dB values; the clause-4.3
//!   start-up value −14 dB is exactly `−14336`.
//! - **MA coefficients** — the staged Q13 `[b1 b2 b3 b4]`; the eq (69)
//!   dot product lands on Q24.
//! - **`ĝ_p`** — Q14, the direct eq (73) column sum (the GA/GB `g_p`
//!   columns are Q14).
//! - **`γ̂`** — the eq (74) column sum on the black-box-pinned 2^13
//!   grid (see [`crate::gain_reconstruct`]).
//! - **`ĝ_c`** — Q1 (a Word16 up to ≈ 16383.5), the grid the eq (75)
//!   excitation build multiplies against the Q13 codevector.
//!
//! ## The eq (71) exponent identity
//!
//! `g'_c = 10^((Ẽ + Ē − E)/20)` with `E = 10·log10(Σc²/40)` reduces on
//! the base-2 grid to
//!
//! ```text
//!   g'_c = 2^( K·Ẽ + C − log2(Σc²)/2 )
//!   K = log2(10)/20 (Q15 literal 5443)
//!   C = (27 + log2(40))/2 + K·30 − 27/2 … folded with the Q27 energy
//! ```
//!
//! because `K · 10·log10(2) = 1/2` exactly — the fixed-codebook energy
//! enters the exponent with weight −1/2 (a square root). The Q27
//! codevector energy accumulator (Q13 codevector, doubled products)
//! folds its `2^27` scale into the constant. The constant is kept in
//! Q24; the composed exponent is split into the integer/Q15-fraction
//! pair [`crate::fx::dsp::pow2`] consumes.
//!
//! The prose pins the *values* of every quantity here (eqs (66)–(74),
//! Ē = 30 dB, Table 9 init); the reference's exact operator schedule
//! inside this computation is not published, so this module's schedule
//! is the crate's own composition of the Table 11 operators and is
//! validated numerically against the float spec transcription and
//! end-to-end against the conformance corpus.

use crate::fx::dsp::{l_extract, log2, mpy_32_16, pow2};
use crate::fx::ops::{
    add, extract_h, extract_l, l_add, l_deposit_l, l_mac, l_shl, l_shr, l_shr_r, l_sub, mult,
    round, sub,
};
use crate::tables::{gain_ga_entry, gain_gb_entry, GAIN_QUANT_MA_PREDICTOR_Q13, MA_NP};
use crate::tables::{GAIN_VQ_COL_GC, GAIN_VQ_COL_GP};

/// Subframe length the eq (66) energy averages over.
pub const L_SUBFR: usize = 40;

/// Clause-4.3 / Table 9 predictor start-up value: −14 dB on Q10.
pub const PAST_QUA_EN_INIT_Q10: i16 = -14336;

/// `log2(10)/20` on the Q15 grid.
const K_LOG2_PER_DB_Q15: i16 = 5443;

/// `(27 + log2(40))/2 + 30·log2(10)/20` on the Q24 grid — the folded
/// eq (71) constant for a Q27 energy accumulator (see module doc).
const EXP_CONST_Q24: i32 = 354_735_042;

/// `20·log10(2)` on the Q12 grid — converts a base-2 log to dB.
const DB_PER_LOG2_Q12: i16 = 24_660;

/// Black-box latitude of the eq (74)/(72) γ̂ grids: per-codebook
/// binary-scale offsets relative to the 2^13 base, for the ĝ_c
/// reconstruction and the Û memory push separately. All-zero is the
/// round-410 uniform-2^13 pin; the residual-divergence hunt sweeps
/// this space against the conformance corpus.
#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GainGridFx {
    /// Extra right-shift of the GA γ column in the ĝ_c reconstruction.
    pub recon_ga: i16,
    /// Extra right-shift of the GB γ column in the ĝ_c reconstruction.
    pub recon_gb: i16,
    /// Extra right-shift of the GA γ column in the eq (72) Û push.
    pub push_ga: i16,
    /// Extra right-shift of the GB γ column in the eq (72) Û push.
    pub push_gb: i16,
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

    /// eq (69): predicted energy `Ẽ^(m)` on the Q24 grid.
    fn predicted_energy_q24(&self) -> i32 {
        let mut acc = 0i32;
        for (b, u) in GAIN_QUANT_MA_PREDICTOR_Q13
            .iter()
            .take(MA_NP)
            .zip(self.past_qua_en.iter())
        {
            acc = l_mac(acc, *b, *u);
        }
        acc
    }

    /// Decodes one subframe's gains from the **codebook-domain**
    /// `(ga, gb)` indices and the Q13 codevector, advancing the
    /// predictor memory (eq (72) decode form).
    ///
    /// # Panics
    ///
    /// Panics if `ga`/`gb` exceed the codebook domains (3-/4-bit
    /// fields cannot).
    pub fn decode(&mut self, ga: usize, gb: usize, code_q13: &[i16; L_SUBFR]) -> DecodedGainsFx {
        let ga_row = gain_ga_entry(ga);
        let gb_row = gain_gb_entry(gb);

        // eq (73): ĝ_p on Q14 (the column extremes sum well inside
        // Word16).
        let gain_pit_q14 = add(ga_row[GAIN_VQ_COL_GP], gb_row[GAIN_VQ_COL_GP]);
        // eq (74): γ̂ on the 2^13 grid — accumulated on Word32: the
        // worst-case column sum (27162 + 14276 = 41438, γ̂ ≈ 5.06)
        // overflows Word16, so the reconstruction keeps the two stage
        // contributions apart until they are combined against g'_c
        // (and takes the eq (72) logarithm of the 32-bit sum).
        let gamma_ga = ga_row[GAIN_VQ_COL_GC];
        let gamma_gb = gb_row[GAIN_VQ_COL_GC];
        // The eq (72) push γ̂ on its (possibly split) grid, Q13 base.
        let gamma_q13_32 = l_add(
            l_shr(l_deposit_l(gamma_ga), self.grid.push_ga),
            l_shr(l_deposit_l(gamma_gb), self.grid.push_gb),
        );

        // eq (66) energy: Σ code² on Q27 (Q13 × Q13 doubled).
        let mut l_ener = 0i32;
        for &c in code_q13.iter() {
            l_ener = l_mac(l_ener, c, c);
        }

        // eq (71) exponent x = K·Ẽ + C − log2(Σc²·2^27)/2, on Q24.
        let e_tilde_q24 = self.predicted_energy_q24();
        let (et_hi, et_lo) = l_extract(e_tilde_q24);
        let k_e_tilde_q24 = mpy_32_16(et_hi, et_lo, K_LOG2_PER_DB_Q15);

        let (log_e, log_f) = log2(l_ener);
        // (log2 value)/2 on Q24: ((e<<15) + f) << 9 is Q24, halved with
        // a rounding shift.
        let log_q24 = l_shl(l_add(l_shl(l_deposit_l(log_e), 15), l_deposit_l(log_f)), 9);
        let half_log_q24 = l_shr_r(log_q24, 1);

        let x_q24 = l_sub(l_add(k_e_tilde_q24, EXP_CONST_Q24), half_log_q24);

        // Split for pow2: integer exponent + Q15 fraction.
        let exp = extract_h(l_shr(x_q24, 8)); // x >> 24
        let frac = extract_l(l_shr(l_sub(x_q24, l_shl(l_deposit_l(exp), 24)), 9));

        // g'_c scaled so the γ̂ product lands on Q1:
        //   pow2(13, frac) = g'_c · 2^(13−exp)   (mantissa full-scale)
        //   mpy_32_16 with γ̂ (2^13 grid)        → γ̂·g'_c·2^(11−exp)
        //   reposition to value·2^17 (Q1 in the high word) and round —
        //   the saturating left shift makes an over-range gain clip to
        //   the Word16 ceiling instead of wrapping.
        let g_prime_scaled = pow2(13, frac); // g'_c · 2^(13−exp), Word32
        let (gp_hi, gp_lo) = l_extract(g_prime_scaled);
        // γ̂·g'_c as the sum of the two per-stage products (full Q13
        // precision, no Word16 γ̂ saturation), each column on its
        // (possibly split) grid.
        let prod = l_add(
            l_shr(mpy_32_16(gp_hi, gp_lo, gamma_ga), self.grid.recon_ga),
            l_shr(mpy_32_16(gp_hi, gp_lo, gamma_gb), self.grid.recon_gb),
        );
        let gain_code_q1 = round(l_shl(prod, add(exp, 6)));

        // eq (72) decode form: Û^(m) = 20·log10(γ̂) pushed on Q10.
        //   log2(γ̂) = log2(gamma_q13_32) − 13; dB = 6.0206·log2.
        let (ge, gf) = log2(gamma_q13_32);
        let lg_q15 = l_add(l_shl(l_deposit_l(sub(ge, 13)), 15), l_deposit_l(gf));
        let (lg_hi, lg_lo) = l_extract(l_shl(lg_q15, 9)); // Q24
        let u_q21 = mpy_32_16(lg_hi, lg_lo, DB_PER_LOG2_Q12); // Q24·Q12 → Q21
        let u_q10 = round(l_shl(u_q21, 5)); // Q21 → Q26 high word = Q10

        for k in (1..MA_NP).rev() {
            self.past_qua_en[k] = self.past_qua_en[k - 1];
        }
        self.past_qua_en[0] = u_q10;

        DecodedGainsFx {
            gain_pit_q14,
            gain_code_q1,
        }
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
            // Tolerance: 1% relative or the Q1 grid step (0.5), the
            // coarser of the two — small gains land on the half-integer
            // Q1 grid by construction.
            let tol = (0.01 * gc_hat.abs()).max(0.3);
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
