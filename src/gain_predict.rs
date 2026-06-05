//! §3.9.1 / §4.1.5 4th-order MA gain prediction (decode side).
//!
//! This module ties the round-231 [`crate::gain_reconstruct`] correction
//! factor `γ̂` into the **actual** quantised fixed-codebook gain
//! `ĝ_c = γ̂ · g'_c`, where `g'_c` is the gain predicted from the
//! recent fixed-codebook log-energy history via a stateful 4-tap MA
//! predictor (spec eqs (65) / (66) / (67) / (69) / (71) / (72) at
//! clause 3.9.1, with the `(-14, -14, -14, -14)` history start-up
//! pinned by spec Table 9, clause 4.3).
//!
//! ## Spec source — clauses 3.9.1 + 4.1.5 + 4.3 (06/2012 Recommendation)
//!
//! Clause 3.9.1 defines the relationship between the fixed-codebook
//! gain `g_c`, the predicted gain `g'_c`, and the correction factor
//! `γ`:
//!
//! - eq (65): `g_c = γ · g'_c`. The encoder transmits the factor `γ`
//!   through the §3.9.2 conjugate-structure VQ; on the decode side the
//!   reconstructed factor is `γ̂` (already produced by
//!   [`crate::gain_reconstruct::reconstruct_gains`]). The actual
//!   quantised fixed-codebook gain is then `ĝ_c = γ̂ · g'_c`.
//! - eq (66): the mean energy of the fixed-codebook contribution is
//!   `E = 10·log10((1/40)·Σ_{n=0..39} c(n)^2)`, computed over the
//!   40-sample subframe codevector.
//! - eq (67): the mean-removed energy `E^(m) = 20·log g_c + E − Ē`
//!   with `Ē = 30 dB` per clause 3.9.1.
//! - eq (69): the predicted `Ẽ^(m) = Σ_{i=1..=4} b_i · Û^(m−i)`,
//!   with `[b1 b2 b3 b4] = [0.68 0.58 0.34 0.19]` per clause 3.9.1
//!   (the Q13 form is staged in
//!   [`crate::tables::GAIN_QUANT_MA_PREDICTOR_Q13`]).
//! - eq (71): the predicted gain `g'_c = 10^((Ẽ^(m) + Ē − E)/20)`,
//!   obtained by substituting `Ẽ^(m)` for `E^(m)` in eq (68).
//! - eq (72): the prediction-error relationship
//!   `U^(m) = E^(m) − Ẽ^(m) = 20·log γ`, so the **quantised** version
//!   `Û^(m) = 20·log10 γ̂` is what the decoder pushes onto the
//!   predictor history each subframe.
//!
//! Clause 4.3 / Table 9 pins the start-up state for the gain
//! predictor: `Û^(k) = −14` for all four history slots.
//!
//! Clause 4.1.5 wires this into the decode order: the §3.9.2
//! conjugate-structure VQ first produces `(ĝ_p, γ̂)`; the §3.9.1
//! predictor then turns `γ̂` into `ĝ_c = γ̂ · g'_c` using the current
//! 40-sample codevector `c(n)` (energy via eq (66)) and the predictor
//! memory `[Û^(m−1), Û^(m−2), Û^(m−3), Û^(m−4)]`; afterwards the
//! predictor advances its history with `Û^(m) = 20·log10 γ̂` (eq (72)
//! decode form).
//!
//! ## Q-format conventions
//!
//! The §3.9.1 spec equations are expressed in floating-point /
//! decibel domains. The staged MA-prediction coefficients
//! `GAIN_QUANT_MA_PREDICTOR_Q13` carry Q13 scaling (each entry is
//! `b_i · 2^13`); this module converts to `f32` at the boundary to
//! match the floating-point boundary type the rest of the decode chain
//! already uses ([`crate::lsp_reconstruct`], [`crate::lsp_to_lp`],
//! [`crate::gain_reconstruct`]). Internally the predictor evaluates
//! eq (69) in `f32` so a future fixed-point variant can swap in by
//! changing the predictor-state type — the algorithm is independent.
//!
//! ## Stateful surface
//!
//! [`GainPredictor`] owns the 4-tap history `[Û^(m−1), Û^(m−2),
//! Û^(m−3), Û^(m−4)]` initialised per spec Table 9 to `(-14, -14,
//! -14, -14)`. [`GainPredictor::predict_and_update`] runs one
//! subframe's spec eq (66) energy, eq (69) predicted energy, eq (71)
//! predicted gain, eq (65) `ĝ_c = γ̂ · g'_c`, and eq (72) decode-form
//! history advance in sequence. The function returns the final
//! `ĝ_c` and an intermediate [`PredictedGain`] carrying `g'_c` and
//! `Ẽ^(m)` for callers that want to inspect the prediction path.
//!
//! Callers that need the predicted gain **without** advancing the
//! history (e.g. encoder-side dry-runs over candidate codevectors)
//! can use [`GainPredictor::predict_only`]; that path leaves the
//! 4-tap memory untouched.

use crate::tables::{GAIN_QUANT_MA_PREDICTOR_Q13, MA_NP};

/// G.729 §3.9.1 mean-energy reference for the fixed-codebook
/// contribution — `Ē = 30 dB`, per clause 3.9.1 directly after
/// eq (67). The same constant appears in eq (68) and eq (71) as the
/// re-mean step that converts mean-removed energy back to absolute
/// gain.
pub const FIXED_CODEBOOK_MEAN_ENERGY_DB: f32 = 30.0;

/// G.729 §4.3 / Table 9 start-up value for the gain predictor
/// history slots `Û^(k)`. Pinned to `−14` for all four slots; the
/// first frame's prediction therefore evaluates to
/// `Ẽ^(0) = −14 · (b1 + b2 + b3 + b4) ≈ −25.9 dB`.
pub const GAIN_PREDICTOR_INIT_DB: f32 = -14.0;

/// G.729 §3.9.1 codevector length over which eq (66) averages — 40
/// samples, matching the subframe size used by every other §3 stage.
pub const CODEVECTOR_LEN: usize = 40;

/// Q13 divisor used when converting the staged
/// [`GAIN_QUANT_MA_PREDICTOR_Q13`] entries to `f32` for the eq (69)
/// summation — `2^13 = 8192`.
const Q13_DIVISOR: f32 = 8192.0;

/// Numerical floor for the eq (66) `log10` argument — guards against
/// the all-zero `c(n)` corner case (a defensible decoder rather than a
/// `NaN`-producing one). The floor `1e-30` is well below the smallest
/// sum-of-squared `i16` codevector magnitude that the §3.8 ISPP
/// 4-pulse construction can produce (each pulse contributes `1` so
/// the minimum non-zero energy is `4 / 40 = 0.1`).
const ENERGY_LOG10_FLOOR: f32 = 1e-30;

/// Result of the §3.9.1 prediction step for one subframe.
///
/// Carries the intermediate quantities the spec equations name so
/// callers can inspect the prediction path (useful for trace tests
/// and conformance harnesses). The boundary types are all `f32` to
/// match the rest of the decode chain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PredictedGain {
    /// `E` from eq (66): mean energy of the current 40-sample
    /// fixed-codebook contribution `c(n)`, in dB.
    pub e_db: f32,
    /// `Ẽ^(m)` from eq (69): predicted mean-removed energy for the
    /// current subframe, in dB. The 4-tap MA sum
    /// `Σ_{i=1..=4} b_i · Û^(m−i)` over the history slots.
    pub e_tilde_db: f32,
    /// `g'_c` from eq (71): predicted fixed-codebook gain
    /// (dimensionless). The exponential conversion of the
    /// mean-removed predicted energy back to an absolute gain.
    pub g_c_prime: f32,
}

/// G.729 §3.9.1 stateful 4-tap MA gain predictor.
///
/// Owns the four-slot history `[Û^(m−1), Û^(m−2), Û^(m−3), Û^(m−4)]`
/// of quantised prediction errors (in dB; per eq (72) `Û` is the
/// quantised `20·log10 γ`). The history is initialised per spec
/// Table 9 to `[-14, -14, -14, -14]` and advanced per subframe by
/// [`GainPredictor::predict_and_update`] (or by the lower-level
/// [`GainPredictor::push_quantised_error`] helper).
///
/// The MA-prediction coefficients `[b1 b2 b3 b4] = [0.68 0.58 0.34
/// 0.19]` (clause 3.9.1) are read from
/// [`crate::tables::GAIN_QUANT_MA_PREDICTOR_Q13`] (Q13-staged) and
/// converted to `f32` at the boundary; the eq (69) dot product runs
/// in `f32`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GainPredictor {
    /// `[Û^(m−1), Û^(m−2), Û^(m−3), Û^(m−4)]` in dB. Index 0 is the
    /// most-recent slot the eq (69) sum weights with `b_1`; index 3
    /// is the oldest slot weighted with `b_4`.
    history_db: [f32; MA_NP],
}

impl Default for GainPredictor {
    fn default() -> Self {
        Self::new()
    }
}

impl GainPredictor {
    /// Builds a fresh predictor with the spec §4.3 / Table 9 start-up
    /// state — all four `Û^(k)` slots at `-14 dB`.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            history_db: [GAIN_PREDICTOR_INIT_DB; MA_NP],
        }
    }

    /// Borrows the four-slot prediction-error history in
    /// `[Û^(m−1), Û^(m−2), Û^(m−3), Û^(m−4)]` order (dB units).
    #[must_use]
    pub fn history_db(&self) -> &[f32; MA_NP] {
        &self.history_db
    }

    /// Evaluates spec eq (69): the 4-tap dot product
    /// `Ẽ^(m) = Σ_{i=1..=4} b_i · Û^(m−i)` over the staged Q13
    /// coefficients and the current history. Returned in dB.
    fn predict_energy_db(&self) -> f32 {
        let mut acc = 0.0_f32;
        for (i, &b_q13) in GAIN_QUANT_MA_PREDICTOR_Q13.iter().take(MA_NP).enumerate() {
            let b_i = f32::from(b_q13) / Q13_DIVISOR;
            acc += b_i * self.history_db[i];
        }
        acc
    }

    /// Computes spec eq (66): the mean energy (in dB) of one 40-sample
    /// fixed-codebook contribution `c(n)`.
    ///
    /// `c(n)` is the §3.8 algebraic codevector; the spec defines it
    /// over `n = 0..40` (one subframe). Callers that have only the
    /// energy itself (e.g. from a fixed-point-staged trace) can bypass
    /// this helper and feed the result directly into
    /// [`GainPredictor::predict_and_update_from_energy`].
    ///
    /// The function floors the `log10` argument at
    /// [`ENERGY_LOG10_FLOOR`] so an all-zero codevector produces a
    /// large-negative but finite `E` rather than `−∞`. This matters
    /// only for the defensive corner case — the §3.8 construction
    /// always has at least four ±1 pulses, so the codevector cannot
    /// be all-zero in normal operation.
    #[must_use]
    pub fn codevector_energy_db(c: &[f32; CODEVECTOR_LEN]) -> f32 {
        let mut sum_sq = 0.0_f32;
        for s in c.iter() {
            sum_sq += s * s;
        }
        let mean = sum_sq / (CODEVECTOR_LEN as f32);
        10.0 * mean.max(ENERGY_LOG10_FLOOR).log10()
    }

    /// Runs one subframe's §3.9.1 prediction path *without* advancing
    /// the predictor's 4-tap history.
    ///
    /// Returns the [`PredictedGain`] carrying `E`, `Ẽ^(m)`, and
    /// `g'_c` per eqs (66), (69), (71). Callers that just need the
    /// final quantised fixed-codebook gain `ĝ_c = γ̂ · g'_c` (eq (65))
    /// can multiply `γ̂` into [`PredictedGain::g_c_prime`].
    #[must_use]
    pub fn predict_only(&self, codevector: &[f32; CODEVECTOR_LEN]) -> PredictedGain {
        let e_db = Self::codevector_energy_db(codevector);
        let e_tilde_db = self.predict_energy_db();
        // eq (71): g'_c = 10^((Ẽ^(m) + Ē − E)/20).
        let exponent = (e_tilde_db + FIXED_CODEBOOK_MEAN_ENERGY_DB - e_db) / 20.0;
        let g_c_prime = 10.0_f32.powf(exponent);
        PredictedGain {
            e_db,
            e_tilde_db,
            g_c_prime,
        }
    }

    /// Same as [`Self::predict_only`] but skips the eq (66) energy
    /// computation by accepting `E` directly. Useful in tests +
    /// downstream code that has the energy from a parallel path.
    #[must_use]
    pub fn predict_only_from_energy(&self, e_db: f32) -> PredictedGain {
        let e_tilde_db = self.predict_energy_db();
        let exponent = (e_tilde_db + FIXED_CODEBOOK_MEAN_ENERGY_DB - e_db) / 20.0;
        let g_c_prime = 10.0_f32.powf(exponent);
        PredictedGain {
            e_db,
            e_tilde_db,
            g_c_prime,
        }
    }

    /// Per spec eq (72) decode form: pushes
    /// `Û^(m) = 20·log10(γ̂)` onto the history, shifting older slots
    /// one step deeper. The oldest slot `Û^(m−4)` is dropped.
    ///
    /// `gamma_hat` should be the §3.9.2 reconstructed correction
    /// factor (from [`crate::gain_reconstruct::reconstruct_gains`]).
    /// `gamma_hat` is floored to [`ENERGY_LOG10_FLOOR`] before the
    /// `log10` to defend against a zero / negative input that would
    /// otherwise produce a `NaN` history slot.
    pub fn push_quantised_error(&mut self, gamma_hat: f32) {
        let u_m_db = 20.0 * gamma_hat.max(ENERGY_LOG10_FLOOR).log10();
        // Shift: history[3] ← history[2], history[2] ← history[1],
        // history[1] ← history[0], history[0] ← Û^(m).
        for k in (1..MA_NP).rev() {
            self.history_db[k] = self.history_db[k - 1];
        }
        self.history_db[0] = u_m_db;
    }

    /// Runs one subframe's §3.9.1 prediction path *and* advances the
    /// 4-tap history (eq (72) decode form).
    ///
    /// The returned tuple `(ĝ_c, gain_path)` carries:
    /// - `ĝ_c = γ̂ · g'_c` (spec eq (65) / (74)) — the actual
    ///   quantised fixed-codebook gain that scales `c(n)` in the
    ///   §3.10 excitation equation (`u(n) = ĝ_p · v(n) + ĝ_c · c(n)`,
    ///   eq (75));
    /// - `gain_path` — the [`PredictedGain`] structure carrying the
    ///   intermediate `E`, `Ẽ^(m)`, and `g'_c` values for inspection
    ///   / trace.
    ///
    /// After the call the predictor's 4-tap history is advanced per
    /// eq (72) decode form: `Û^(m) = 20·log10(γ̂)` is pushed onto
    /// slot 0 and the oldest slot is dropped.
    pub fn predict_and_update(
        &mut self,
        codevector: &[f32; CODEVECTOR_LEN],
        gamma_hat: f32,
    ) -> (f32, PredictedGain) {
        let path = self.predict_only(codevector);
        let g_c_hat = gamma_hat * path.g_c_prime;
        self.push_quantised_error(gamma_hat);
        (g_c_hat, path)
    }

    /// Same as [`Self::predict_and_update`] but takes the pre-computed
    /// codevector energy `E` (in dB) instead of the 40-sample
    /// codevector. Convenient when the caller has already computed the
    /// energy on a parallel path.
    pub fn predict_and_update_from_energy(
        &mut self,
        e_db: f32,
        gamma_hat: f32,
    ) -> (f32, PredictedGain) {
        let path = self.predict_only_from_energy(e_db);
        let g_c_hat = gamma_hat * path.g_c_prime;
        self.push_quantised_error(gamma_hat);
        (g_c_hat, path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tables::GAIN_QUANT_MA_PREDICTOR_Q13;

    /// Per spec Table 9 / clause 4.3 the four-slot gain-predictor
    /// history starts at `[-14, -14, -14, -14]` dB. Verifies the
    /// const constructor + the [`GainPredictor::history_db`] borrow.
    #[test]
    fn new_predictor_matches_table9_init() {
        let p = GainPredictor::new();
        for slot in p.history_db() {
            assert!((slot - GAIN_PREDICTOR_INIT_DB).abs() < 1e-9);
        }
        assert_eq!(GAIN_PREDICTOR_INIT_DB, -14.0);
    }

    /// First subframe after init: eq (69) sum is
    /// `Ẽ^(0) = -14 · (b1 + b2 + b3 + b4)`. Computed from the staged
    /// Q13 coefficients directly so a future drift in the table's
    /// numeric content surfaces in this test.
    #[test]
    fn first_subframe_predicted_energy_matches_sum_of_taps() {
        let p = GainPredictor::new();
        let sum_b: f32 = (0..MA_NP)
            .map(|i| f32::from(GAIN_QUANT_MA_PREDICTOR_Q13[i]) / Q13_DIVISOR)
            .sum();
        let expected = -14.0 * sum_b;
        // Use a non-zero codevector so the energy path is also
        // exercised; we only assert on `e_tilde_db` here.
        let mut c = [0.0_f32; CODEVECTOR_LEN];
        c[0] = 1.0;
        c[10] = -1.0;
        c[20] = 1.0;
        c[35] = -1.0;
        let path = p.predict_only(&c);
        assert!(
            (path.e_tilde_db - expected).abs() < 1e-4,
            "Ẽ^(0) = {} expected {} (sum_b = {})",
            path.e_tilde_db,
            expected,
            sum_b,
        );
    }

    /// Eq (66): a 40-sample c(n) with exactly four ±1 pulses (the
    /// minimum-energy spec §3.8 codevector shape) has
    /// `E = 10·log10(4/40) = 10·log10(0.1) = -10 dB`. The helper
    /// must compute it inside 1e-5 dB.
    #[test]
    fn energy_eq66_matches_four_unit_pulses() {
        let mut c = [0.0_f32; CODEVECTOR_LEN];
        c[0] = 1.0;
        c[6] = -1.0;
        c[17] = 1.0;
        c[34] = -1.0;
        let e = GainPredictor::codevector_energy_db(&c);
        assert!((e - (-10.0_f32)).abs() < 1e-5, "E = {e} expected -10 dB");
    }

    /// Eq (66) with all-zero codevector falls back on the
    /// `ENERGY_LOG10_FLOOR` so the helper is finite (no `−∞`/`NaN`).
    /// `10·log10(1e-30) = -300 dB`.
    #[test]
    fn energy_eq66_zero_input_is_finite() {
        let c = [0.0_f32; CODEVECTOR_LEN];
        let e = GainPredictor::codevector_energy_db(&c);
        assert!(e.is_finite());
        // 1e-30 floor → 10 · -30 = -300 dB
        assert!((e - (-300.0_f32)).abs() < 1e-3, "E = {e}");
    }

    /// Eq (71): when `Ẽ^(m) = E - Ē`, the predicted gain is exactly
    /// `g'_c = 1.0`. This pins the sign convention of the
    /// `(Ẽ + Ē − E) / 20` exponent.
    #[test]
    fn predicted_gain_unity_when_predicted_energy_cancels() {
        // Build a custom predictor with a history that produces a
        // specific Ẽ^(m). We accept the predictor-internal Ẽ^(m) and
        // pass a matching codevector energy via the helper.
        let p = GainPredictor::new();
        let e_tilde = p.predict_energy_db();
        // E = Ẽ + Ē → exponent = 0 → g'_c = 1.0
        let e_db = e_tilde + FIXED_CODEBOOK_MEAN_ENERGY_DB;
        let path = p.predict_only_from_energy(e_db);
        assert!((path.g_c_prime - 1.0).abs() < 1e-5);
    }

    /// Eq (72) decode form: pushing `γ̂ = 1.0` produces `Û^(m) = 0
    /// dB`. The most-recent history slot must be 0 after the push,
    /// older slots must shift one step deeper, and the oldest slot
    /// must drop. `γ̂ = 10` produces `Û^(m) = 20 dB`.
    #[test]
    fn push_quantised_error_eq72_decode_form() {
        let mut p = GainPredictor::new();
        // Push γ̂ = 1.0 → Û^(m) = 0
        p.push_quantised_error(1.0);
        assert!(p.history_db()[0].abs() < 1e-6);
        // All other slots still at -14 (the original init state
        // shifted one step).
        for slot in &p.history_db()[1..] {
            assert!((slot - (-14.0)).abs() < 1e-6);
        }
        // Push γ̂ = 10 → Û^(m) = 20
        p.push_quantised_error(10.0);
        assert!((p.history_db()[0] - 20.0).abs() < 1e-5);
        // Slot 1 should now carry the previous Û^(m) = 0.
        assert!(p.history_db()[1].abs() < 1e-6);
        // Slot 2 should carry one of the original -14 slots.
        assert!((p.history_db()[2] - (-14.0)).abs() < 1e-6);
    }

    /// Per spec eq (72) decode form `Û^(m) = 20·log10(γ̂)`. A `γ̂ =
    /// 0.1` push must produce `Û^(m) = -20 dB`. Defends against a
    /// missing factor of 20 or a sign flip in the predictor update.
    #[test]
    fn push_quantised_error_decibel_scaling_is_20log10() {
        let mut p = GainPredictor::new();
        p.push_quantised_error(0.1);
        assert!((p.history_db()[0] - (-20.0)).abs() < 1e-5);
    }

    /// Zero / negative γ̂ floors at `ENERGY_LOG10_FLOOR` so the
    /// history slot is finite (very negative dB), not NaN. Defensive
    /// — well-formed §3.9.2 output is always positive.
    #[test]
    fn push_quantised_error_floor_for_nonpositive_input() {
        let mut p = GainPredictor::new();
        p.push_quantised_error(0.0);
        assert!(p.history_db()[0].is_finite());
        // 1e-30 floor → 20 · -30 = -600 dB
        assert!((p.history_db()[0] - (-600.0)).abs() < 1e-3);
    }

    /// End-to-end predict+update: γ̂ = 1.0 with a 4-pulse codevector
    /// (`E = -10 dB`) produces `ĝ_c = g'_c` (since γ̂ = 1). The
    /// history must advance with `Û^(m) = 0`.
    #[test]
    fn predict_and_update_consistent_with_pieces() {
        let mut c = [0.0_f32; CODEVECTOR_LEN];
        c[0] = 1.0;
        c[6] = -1.0;
        c[17] = 1.0;
        c[34] = -1.0;
        let mut p = GainPredictor::new();
        let path_before = p.predict_only(&c);
        let (g_c_hat, path) = p.predict_and_update(&c, 1.0);
        assert!((path.e_db - (-10.0)).abs() < 1e-5);
        assert!((path.e_tilde_db - path_before.e_tilde_db).abs() < 1e-7);
        assert!((path.g_c_prime - path_before.g_c_prime).abs() < 1e-7);
        // γ̂ = 1 → ĝ_c = g'_c
        assert!((g_c_hat - path.g_c_prime).abs() < 1e-7);
        // Eq (72) decode form: Û^(m) = 20·log10(1) = 0
        assert!(p.history_db()[0].abs() < 1e-6);
    }

    /// `predict_only` is side-effect-free: the history is unchanged
    /// after the call. Pinned for the encoder-side dry-run use case.
    #[test]
    fn predict_only_leaves_history_untouched() {
        let p = GainPredictor::new();
        let before = *p.history_db();
        let mut c = [0.0_f32; CODEVECTOR_LEN];
        c[0] = 1.0;
        let _ = p.predict_only(&c);
        // Borrow-check: `p` is immutable so the test is structural; the
        // real signal is that the call compiles against `&self`.
        let after = *p.history_db();
        assert_eq!(before, after);
    }

    /// Multi-subframe drift check: drive the predictor with a
    /// constant γ̂ for 6 subframes. After 4 pushes the history is
    /// fully populated with the steady-state Û = 20·log10(γ̂) and
    /// the next `Ẽ^(m)` is therefore exactly
    /// `Û · (b1 + b2 + b3 + b4)`. Verifies eq (69) drives the
    /// expected steady state once the start-up `-14` slots have
    /// cycled out.
    #[test]
    fn steady_state_after_four_pushes_matches_sum_of_taps() {
        let gamma = 0.5_f32;
        let u_steady = 20.0 * gamma.log10(); // = -6.0206...
        let sum_b: f32 = (0..MA_NP)
            .map(|i| f32::from(GAIN_QUANT_MA_PREDICTOR_Q13[i]) / Q13_DIVISOR)
            .sum();
        let mut p = GainPredictor::new();
        for _ in 0..MA_NP {
            p.push_quantised_error(gamma);
        }
        // Now every history slot is u_steady.
        for slot in p.history_db() {
            assert!((slot - u_steady).abs() < 1e-5);
        }
        // Eq (69): Ẽ^(m) = u_steady · sum_b
        let e_tilde = p.predict_only_from_energy(0.0).e_tilde_db;
        let expected = u_steady * sum_b;
        assert!((e_tilde - expected).abs() < 1e-4);
    }

    /// History order pin: index 0 carries the most-recent Û (the slot
    /// weighted by `b_1 ≈ 0.68`), index 3 carries the oldest (weighted
    /// by `b_4 ≈ 0.19`). Build a hand-set history where slot k carries
    /// the magic value `100 + k`; the expected `Ẽ^(m)` is then
    /// `Σ_k b_{k+1} · (100 + k)`, with `b_{k+1}` read from index `k`
    /// of the staged table. Defensive against a reversed history
    /// order or a 1-based vs 0-based index mistake.
    #[test]
    fn history_index_zero_is_most_recent_weighted_by_b1() {
        // Construct a predictor with a hand-set history. We can't set
        // it directly via a public field, so push four distinct γ̂
        // values whose Û^(m) = 20·log10(γ̂) lands on 103, 102, 101,
        // 100 (in chronological push order: 100 → oldest slot, 103 →
        // most-recent after 4 pushes).
        let mut p = GainPredictor::new();
        // Push 100 first (oldest slot after 4 pushes)
        p.push_quantised_error(10.0_f32.powf(100.0 / 20.0));
        p.push_quantised_error(10.0_f32.powf(101.0 / 20.0));
        p.push_quantised_error(10.0_f32.powf(102.0 / 20.0));
        p.push_quantised_error(10.0_f32.powf(103.0 / 20.0));
        // After 4 pushes the original -14 slots are gone. Slot 0 holds
        // the most-recent push (=103), slot 3 holds the oldest (=100).
        assert!((p.history_db()[0] - 103.0).abs() < 1e-3);
        assert!((p.history_db()[1] - 102.0).abs() < 1e-3);
        assert!((p.history_db()[2] - 101.0).abs() < 1e-3);
        assert!((p.history_db()[3] - 100.0).abs() < 1e-3);
        // Eq (69) weighting check: Σ_k b_{k+1} · history[k]
        let sum: f32 = (0..MA_NP)
            .map(|k| {
                let b = f32::from(GAIN_QUANT_MA_PREDICTOR_Q13[k]) / Q13_DIVISOR;
                b * p.history_db()[k]
            })
            .sum();
        let e_tilde = p.predict_only_from_energy(0.0).e_tilde_db;
        assert!((e_tilde - sum).abs() < 1e-3);
    }

    /// Eq (65): g_c = γ · g'_c. With a hand-controlled g'_c (forced
    /// to 1.0 via the energy-cancel trick) the returned ĝ_c equals
    /// γ̂ exactly across a sweep of γ̂ values.
    #[test]
    fn eq65_g_c_equals_gamma_times_g_c_prime() {
        for &gamma_hat in &[0.25_f32, 0.5, 1.0, 1.5, 2.5, 4.0] {
            let mut p = GainPredictor::new();
            let e_tilde = p.predict_energy_db();
            let e_db = e_tilde + FIXED_CODEBOOK_MEAN_ENERGY_DB; // g'_c = 1
            let (g_c_hat, path) = p.predict_and_update_from_energy(e_db, gamma_hat);
            assert!((path.g_c_prime - 1.0).abs() < 1e-5);
            assert!(
                (g_c_hat - gamma_hat).abs() < 1e-5,
                "ĝ_c = {g_c_hat} expected γ̂ = {gamma_hat}",
            );
        }
    }
}
