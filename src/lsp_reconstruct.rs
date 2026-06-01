//! §3.2.4 LSP-coefficient reconstruction.
//!
//! Implements the decoder-side reconstruction of the quantised LSF
//! coefficients `ω̂_i^(m)` from the four transmitted LSP codeword
//! fields `(L0, L1, L2, L3)`. This is the **inverse** of the encoder
//! quantisation described in spec §3.2.4 and is also re-used by the
//! decoder per spec §4.1.1 ("LP params: L0/L1/L2/L3 → LSP (reuse
//! clause 3.2.4)").
//!
//! Round 201 already wired the bit-exact numeric inputs:
//!
//! * the first-stage codebook [`crate::tables::LSP_QUANT_CODEBOOK_L1_Q13`]
//!   indexed by the 7-bit `L1` parameter;
//! * the packed second-stage codebook
//!   [`crate::tables::LSP_QUANT_CODEBOOK_L2_Q13`] split by
//!   [`crate::tables::lsp_l2_entry`] / [`crate::tables::lsp_l3_entry`]
//!   into the L2 (lower-5) and L3 (upper-5) 5-bit halves;
//! * the MA-predictor cube [`crate::tables::LSP_MA_PREDICTOR_FG_Q15`]
//!   with the per-mode [`crate::tables::LSP_MA_PREDICTOR_FG_SUM_Q15`]
//!   `(1 − Σ_k fg[mode][k][i])` factor.
//!
//! Round 207 (this module) wires the spec's §3.2.4 algorithmic steps
//! around those tables:
//!
//! 1. **Codebook sum** (spec eq (19)) — `l̂_i = L1_i(L1) + L2_i(L2)`
//!    for `i ∈ 1..=5`, `L1_i(L1) + L3_{i-5}(L3)` for `i ∈ 6..=10`.
//! 2. **Rearrangement** — for each adjacent pair `(l̂_{i-1}, l̂_i)`
//!    with `l̂_{i-1} > l̂_i − J`, replace them with the midpoint
//!    `((l̂_i + l̂_{i-1}) ± J) / 2`. Per spec §3.2.4 this is run
//!    **twice**: first with `J = 0.0012`, then `J = 0.0006`.
//! 3. **MA-prediction reconstruction** (spec eq (20)) —
//!    `ω̂_i^(m) = fg_sum[mode][i] · l̂_i^(m)
//!               + Σ_{k=1..=4} fg[mode][k-1][i] · l̂_i^(m-k)`
//!    where `fg_sum[mode] = 1 − Σ_k fg[mode][k]` is the round-201
//!    pre-tabulated factor.
//! 4. **Stability clamp** (spec §3.2.4 steps 1–4) — sort ascending,
//!    floor entry 0 at 0.005, enforce minimum adjacent gap 0.0391,
//!    ceil at 3.135.
//!
//! The inputs and outputs of this module are **real-valued** `f32`
//! LSF coefficients in the normalised frequency domain `[0, π]`. The
//! Q-format conversion happens at the codebook-lookup boundary: each
//! Q13 codebook entry is treated as `v / 8192.0`, each Q15 `fg`
//! coefficient as `v / 32768.0`. The clamp constants (`0.005`,
//! `0.0391`, `3.135`) are spec literals.
//!
//! At start-up the MA history `l̂_i^(m-k)` is initialised to
//! `i · π / 11` for all `k ≥ 1`, per spec §3.2.4 ("at start up the
//! initial values of l̂_i^(m-k) are given by `l̂_i = iπ/11` for all
//! `k < 0`"). [`LspReconstructor::new`] applies that initialisation;
//! [`LspReconstructor::reconstruct_frame`] consumes a fresh
//! `(L0, L1, L2, L3)` quadruple, updates the history, and returns
//! the freshly reconstructed `ω̂^(m)` vector.

use core::f32::consts::PI;

use crate::tables::{
    self, LSP_MA_PREDICTOR_FG_Q15, LSP_MA_PREDICTOR_FG_SUM_Q15, LSP_QUANT_CODEBOOK_L1_Q13, M,
    MA_NP, NC0, NC1,
};

/// Spec §3.2.4 step 2 floor on `ω̂_1` — `ω̂_1 < 0.005 ⇒ ω̂_1 := 0.005`.
pub const CLAMP_FLOOR: f32 = 0.005;

/// Spec §3.2.4 step 3 minimum adjacent gap (`ω̂_i − ω̂_{i-1} < 0.0391
/// ⇒ ω̂_i := ω̂_{i-1} + 0.0391`).
pub const CLAMP_MIN_GAP: f32 = 0.0391;

/// Spec §3.2.4 step 4 ceiling on `ω̂_{10}` — `ω̂_{10} > 3.135 ⇒
/// ω̂_{10} := 3.135`.
pub const CLAMP_CEIL: f32 = 3.135;

/// First rearrangement-pass minimum adjacent distance — spec §3.2.4.
pub const REARRANGE_J1: f32 = 0.0012;

/// Second rearrangement-pass minimum adjacent distance — spec §3.2.4.
pub const REARRANGE_J2: f32 = 0.0006;

/// Q15 unit (`2^15`) used to convert Q15 `fg` literals to `f32`.
const Q15_UNIT_F32: f32 = 32768.0;

/// Q13 unit (`2^13`) used to convert Q13 codebook literals to `f32`.
const Q13_UNIT_F32: f32 = 8192.0;

/// Errors that the §3.2.4 reconstruction can surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LspReconstructError {
    /// `L0` parameter outside `0..(1 << L0_BITS) == 0..2`.
    L0OutOfRange { value: usize },
    /// `L1` parameter outside `0..NC0` (NC0 = 128).
    L1OutOfRange { value: usize },
    /// `L2` parameter outside `0..NC1` (NC1 = 32).
    L2OutOfRange { value: usize },
    /// `L3` parameter outside `0..NC1` (NC1 = 32).
    L3OutOfRange { value: usize },
}

impl core::fmt::Display for LspReconstructError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::L0OutOfRange { value } => {
                write!(f, "G.729 §3.2.4: L0 = {value} outside 0..2")
            }
            Self::L1OutOfRange { value } => {
                write!(f, "G.729 §3.2.4: L1 = {value} outside 0..{NC0}")
            }
            Self::L2OutOfRange { value } => {
                write!(f, "G.729 §3.2.4: L2 = {value} outside 0..{NC1}")
            }
            Self::L3OutOfRange { value } => {
                write!(f, "G.729 §3.2.4: L3 = {value} outside 0..{NC1}")
            }
        }
    }
}

impl std::error::Error for LspReconstructError {}

/// Combine the L1 / L2 / L3 codebook entries into the current frame's
/// residual `l̂_i^(m)` per spec §3.2.4 eq (19).
///
/// Indices are validated against the codebook dimensions and an
/// in-range result is returned as a 10-coordinate `f32` vector in the
/// normalised frequency domain.
///
/// # Errors
///
/// Returns [`LspReconstructError::L1OutOfRange`] / `L2OutOfRange` /
/// `L3OutOfRange` when any index lies outside its codebook.
pub fn codebook_sum(l1: usize, l2: usize, l3: usize) -> Result<[f32; M], LspReconstructError> {
    if l1 >= NC0 {
        return Err(LspReconstructError::L1OutOfRange { value: l1 });
    }
    if l2 >= NC1 {
        return Err(LspReconstructError::L2OutOfRange { value: l2 });
    }
    if l3 >= NC1 {
        return Err(LspReconstructError::L3OutOfRange { value: l3 });
    }

    let l1_row = &LSP_QUANT_CODEBOOK_L1_Q13[l1];
    let l2_lo = tables::lsp_l2_entry(l2);
    let l3_hi = tables::lsp_l3_entry(l3);

    let mut residual = [0.0_f32; M];
    for i in 0..M / 2 {
        let stage1 = f32::from(l1_row[i]) / Q13_UNIT_F32;
        let stage2 = f32::from(l2_lo[i]) / Q13_UNIT_F32;
        residual[i] = stage1 + stage2;
    }
    for i in M / 2..M {
        let stage1 = f32::from(l1_row[i]) / Q13_UNIT_F32;
        let stage2 = f32::from(l3_hi[i - M / 2]) / Q13_UNIT_F32;
        residual[i] = stage1 + stage2;
    }
    Ok(residual)
}

/// In-place rearrangement pass — spec §3.2.4 figure F0013-01.
///
/// For each adjacent pair `(l̂_{i-1}, l̂_i)` where
/// `l̂_{i-1} > l̂_i − J`, replace both with the centred pair
/// `(l̂_i + l̂_{i-1} − J)/2` and `(l̂_i + l̂_{i-1} + J)/2`. The pair
/// is walked from `i = 1` to `i = M − 1` so each updated coefficient
/// participates in the next comparison; the post-pass invariant is
/// `l̂_i − l̂_{i-1} ≥ J` for every adjacent pair.
///
/// Per spec §3.2.4 this pass is run **twice**: first with
/// [`REARRANGE_J1`] = 0.0012, then with [`REARRANGE_J2`] = 0.0006.
/// The helper does one pass; [`rearrange_twice`] runs both.
pub fn rearrange_pass(coefs: &mut [f32; M], j: f32) {
    for i in 1..M {
        if coefs[i - 1] > coefs[i] - j {
            // Sequential per the spec figure: l̂_{i-1} is updated to
            // (old_l̂_i + old_l̂_{i-1} − J)/2 first, then l̂_i is
            // updated to (old_l̂_i + old_l̂_{i-1} + J)/2 = new_l̂_{i-1} + J.
            let sum = coefs[i] + coefs[i - 1];
            coefs[i - 1] = (sum - j) / 2.0;
            coefs[i] = (sum + j) / 2.0;
        }
    }
}

/// Run [`rearrange_pass`] twice per spec §3.2.4 — first with
/// `J = 0.0012`, then with `J = 0.0006`.
pub fn rearrange_twice(coefs: &mut [f32; M]) {
    rearrange_pass(coefs, REARRANGE_J1);
    rearrange_pass(coefs, REARRANGE_J2);
}

/// Apply the §3.2.4 4-step stability check to a freshly reconstructed
/// LSF vector in place.
///
/// 1. Sort `ω̂` in increasing value (the LSP-stability prerequisite
///    that `ω̂_1 < ω̂_2 < … < ω̂_10` per spec §3.2.4 step 1).
/// 2. If `ω̂_1 < 0.005`, set `ω̂_1 := 0.005`.
/// 3. For `i = 2..=10`, if `ω̂_i − ω̂_{i-1} < 0.0391`, set
///    `ω̂_i := ω̂_{i-1} + 0.0391`.
/// 4. If `ω̂_10 > 3.135`, set `ω̂_10 := 3.135`.
///
/// The clamp constants are spec literals and exposed as
/// [`CLAMP_FLOOR`] / [`CLAMP_MIN_GAP`] / [`CLAMP_CEIL`].
pub fn stability_clamp(coefs: &mut [f32; M]) {
    // Step 1 — sort ascending. f32 doesn't impl Ord so use partial_cmp;
    // any NaN sentinel is treated as "greater than" (placed at the
    // end), but real LSFs are finite-valued inputs.
    coefs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));

    // Step 2 — floor the first coefficient.
    if coefs[0] < CLAMP_FLOOR {
        coefs[0] = CLAMP_FLOOR;
    }

    // Step 3 — minimum adjacent gap.
    for i in 1..M {
        let min = coefs[i - 1] + CLAMP_MIN_GAP;
        if coefs[i] < min {
            coefs[i] = min;
        }
    }

    // Step 4 — ceiling on the last coefficient.
    if coefs[M - 1] > CLAMP_CEIL {
        coefs[M - 1] = CLAMP_CEIL;
    }
}

/// Apply spec eq (20) to a single residual + MA history pair.
///
/// `ω̂_i^(m) = fg_sum[mode][i] · l̂_i^(m)
///            + Σ_{k=1..=MA_NP} fg[mode][k-1][i] · l̂_i^(m-k)`
///
/// All inputs are interpreted as real-valued LSF coordinates; the
/// Q15 `fg` / `fg_sum` literals are converted at the boundary.
///
/// `mode` is the 1-bit `L0` field; callers must pass `0` or `1`.
///
/// # Panics
///
/// Panics if `mode >= 2`. Callers should validate `mode` against
/// the 1-bit `L0` field before invoking this helper.
fn ma_predict_one(mode: usize, residual: &[f32; M], history: &[[f32; M]; MA_NP]) -> [f32; M] {
    let fg_plane = &LSP_MA_PREDICTOR_FG_Q15[mode];
    let fg_sum_row = &LSP_MA_PREDICTOR_FG_SUM_Q15[mode];

    let mut out = [0.0_f32; M];
    for i in 0..M {
        let sum_factor = f32::from(fg_sum_row[i]) / Q15_UNIT_F32;
        let mut acc = sum_factor * residual[i];
        for k in 0..MA_NP {
            let coef = f32::from(fg_plane[k][i]) / Q15_UNIT_F32;
            acc += coef * history[k][i];
        }
        out[i] = acc;
    }
    out
}

/// §3.2.4 LSP-frame reconstruction state.
///
/// Holds the four-frame MA history of previous residuals
/// `l̂^(m-1)…l̂^(m-4)` per spec eq (20). At construction the history
/// is initialised to the spec-defined start-up vector
/// `l̂_i = i · π / 11` for `i ∈ 1..=10` and every history slot
/// `k ∈ 1..=MA_NP` (spec §3.2.4: "at start up the initial values of
/// `l̂_i^(m-k)` are given by `l̂_i = iπ/11` for all k < 0").
///
/// Call [`Self::reconstruct_frame`] once per 10 ms frame; the returned
/// `[f32; 10]` is the freshly reconstructed `ω̂^(m)` vector after
/// rearrangement and stability clamp. The internal history advances
/// by pushing the *unclamped* eq-(19) residual (not the post-clamp
/// LSF) into slot 0, per spec eq (20) which references the residuals
/// `l̂_i^(m-k)`, not the post-clamp output `ω̂_i^(m-k)`.
#[derive(Debug, Clone)]
pub struct LspReconstructor {
    /// `history[k][i] == l̂_i^(m-1-k)` — the residual from `1+k` frames
    /// ago, oldest at index `MA_NP - 1`. After every
    /// [`Self::reconstruct_frame`] call the array is shifted so the
    /// just-consumed current frame becomes `history[0]`.
    history: [[f32; M]; MA_NP],
}

impl Default for LspReconstructor {
    fn default() -> Self {
        Self::new()
    }
}

impl LspReconstructor {
    /// Build a reconstructor with the spec start-up history
    /// `l̂_i^(m-k) = i · π / 11` for every history slot.
    #[must_use]
    pub fn new() -> Self {
        let mut row = [0.0_f32; M];
        // Spec uses 1-based `i ∈ 1..=10`; our slice index is 0-based.
        for (i, entry) in row.iter_mut().enumerate() {
            *entry = ((i + 1) as f32) * PI / 11.0;
        }
        Self {
            history: [row; MA_NP],
        }
    }

    /// Borrow the current MA history for testing / inspection.
    #[must_use]
    pub fn history(&self) -> &[[f32; M]; MA_NP] {
        &self.history
    }

    /// Run one frame of §3.2.4 reconstruction.
    ///
    /// Performs (1) codebook sum, (2) rearrangement twice
    /// (J = 0.0012, then J = 0.0006), (3) MA-prediction reconstruction
    /// per spec eq (20) using the `L0`-selected predictor mode and
    /// the internal four-frame history, and (4) the 4-step stability
    /// clamp. Returns the reconstructed `ω̂^(m)` vector.
    ///
    /// The internal MA history is updated by pushing the post-rearrange
    /// residual `l̂^(m)` into slot 0 and discarding the oldest slot,
    /// per spec eq (20).
    ///
    /// # Errors
    ///
    /// Returns the appropriate [`LspReconstructError`] when any of
    /// `l0`, `l1`, `l2`, `l3` lies outside its respective codeword
    /// domain.
    pub fn reconstruct_frame(
        &mut self,
        l0: usize,
        l1: usize,
        l2: usize,
        l3: usize,
    ) -> Result<[f32; M], LspReconstructError> {
        if l0 >= (1 << tables::L0_BITS) {
            return Err(LspReconstructError::L0OutOfRange { value: l0 });
        }

        // Step 1 — codebook sum (eq (19)).
        let mut residual = codebook_sum(l1, l2, l3)?;

        // Step 2 — rearrangement twice.
        rearrange_twice(&mut residual);

        // Step 3 — MA-prediction reconstruction (eq (20)).
        let mut omega = ma_predict_one(l0, &residual, &self.history);

        // Step 4 — stability clamp.
        stability_clamp(&mut omega);

        // Advance MA history: the slot that was `l̂^(m-1)` becomes
        // `l̂^(m-2)`, etc., and the just-consumed residual `l̂^(m)`
        // becomes the new slot 0 (`l̂^(m-1)` for the next frame).
        for k in (1..MA_NP).rev() {
            self.history[k] = self.history[k - 1];
        }
        self.history[0] = residual;

        Ok(omega)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// At construction the MA history is the spec start-up vector
    /// `l̂_i = i · π / 11` for every slot — locks the initialisation
    /// step the spec calls out explicitly.
    #[test]
    fn new_history_matches_spec_startup_vector() {
        let r = LspReconstructor::new();
        for (k, plane) in r.history().iter().enumerate() {
            for (i, got) in plane.iter().enumerate() {
                let expected = ((i + 1) as f32) * PI / 11.0;
                assert!(
                    (got - expected).abs() < 1e-6,
                    "history[{k}][{i}] = {got}; expected {expected}",
                );
            }
        }
    }

    /// `codebook_sum(l1, l2, l3)` is the real-domain value of the
    /// Q13 entries summed: stage1[i] + (L2 lower-5 if i<5 else L3
    /// upper-5). Because the staged second-stage table packs both
    /// halves of the row side-by-side, with `(L2, L3) = (0, 0)` the
    /// per-coordinate stage-2 contribution is just the staged Q13
    /// row entry at index `i`. This pins the boundary conversion
    /// independent of the rearrangement / MA / clamp steps.
    #[test]
    fn codebook_sum_l1_l2_l3_zero_indices() {
        let residual = codebook_sum(0, 0, 0).expect("indices in range");
        for (i, got) in residual.iter().enumerate() {
            let stage1 = f32::from(LSP_QUANT_CODEBOOK_L1_Q13[0][i]) / 8192.0;
            // With L2 = L3 = 0, both halves come from row 0 of the
            // packed second-stage codebook; the i'th element is its
            // own Q13 literal in both halves of the row.
            let stage2 = f32::from(tables::LSP_QUANT_CODEBOOK_L2_Q13[0][i]) / 8192.0;
            let expected = stage1 + stage2;
            assert!(
                (got - expected).abs() < 1e-6,
                "residual[{i}] = {got}; expected {expected}",
            );
        }
    }

    /// Out-of-range indices are reported, not silently clamped.
    #[test]
    fn codebook_sum_rejects_out_of_range_indices() {
        assert_eq!(
            codebook_sum(NC0, 0, 0),
            Err(LspReconstructError::L1OutOfRange { value: NC0 }),
        );
        assert_eq!(
            codebook_sum(0, NC1, 0),
            Err(LspReconstructError::L2OutOfRange { value: NC1 }),
        );
        assert_eq!(
            codebook_sum(0, 0, NC1),
            Err(LspReconstructError::L3OutOfRange { value: NC1 }),
        );
    }

    /// An L0 outside the 1-bit domain is rejected — the start-up
    /// state guards against caller-side bit-masking drift.
    #[test]
    fn reconstruct_frame_rejects_l0_out_of_range() {
        let mut r = LspReconstructor::new();
        assert_eq!(
            r.reconstruct_frame(2, 0, 0, 0),
            Err(LspReconstructError::L0OutOfRange { value: 2 }),
        );
    }

    /// Rearrangement enforces the post-pass invariant
    /// `l̂_i − l̂_{i-1} ≥ J − ε` for every adjacent pair: the
    /// rearrangement is *equality-only* at the disturbed positions,
    /// so the resulting `l̂` is monotone increasing modulo numerical
    /// noise.
    #[test]
    fn rearrange_enforces_minimum_distance() {
        let j = 0.01_f32;
        // Construct a manifestly violating sequence.
        let mut v: [f32; M] = [0.10, 0.05, 0.30, 0.25, 0.40, 0.39, 0.60, 0.59, 0.80, 0.79];
        rearrange_pass(&mut v, j);
        for i in 1..M {
            // Allow a tiny float epsilon for the (sum ± j)/2 arithmetic.
            assert!(
                v[i] - v[i - 1] >= j - 1e-6,
                "pair ({}, {}) at i={i} violates min distance {j}: gap {}",
                v[i - 1],
                v[i],
                v[i] - v[i - 1],
            );
        }
    }

    /// A monotone-increasing input with adjacent gaps larger than `J`
    /// is left unchanged by the rearrangement.
    #[test]
    fn rearrange_is_a_noop_when_already_separated() {
        let j = 0.01_f32;
        let baseline: [f32; M] = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00];
        let mut v = baseline;
        rearrange_pass(&mut v, j);
        for i in 0..M {
            assert!(
                (v[i] - baseline[i]).abs() < 1e-6,
                "no-op rearrangement disturbed index {i}: {} vs {}",
                v[i],
                baseline[i],
            );
        }
    }

    /// `rearrange_twice` finishes the second pass with the tighter
    /// `J2 < J1` margin — every adjacent gap is at least `J2`.
    #[test]
    fn rearrange_twice_finishes_at_j2_margin() {
        let mut v: [f32; M] = [0.10, 0.05, 0.30, 0.25, 0.40, 0.39, 0.60, 0.59, 0.80, 0.79];
        rearrange_twice(&mut v);
        for i in 1..M {
            assert!(
                v[i] - v[i - 1] >= REARRANGE_J2 - 1e-6,
                "after rearrange_twice, gap at i={i} = {}; expected ≥ {}",
                v[i] - v[i - 1],
                REARRANGE_J2,
            );
        }
    }

    /// The stability clamp's 4 steps independently bring the head,
    /// the inner gaps, and the tail into range. Construct an input
    /// that violates every one of them and verify the post-clamp
    /// invariants.
    #[test]
    fn stability_clamp_enforces_floor_gap_and_ceil() {
        // Deliberate violations: unordered, head below floor, two
        // adjacent gaps below the min, tail above ceil.
        let mut v: [f32; M] = [
            0.0,   // < CLAMP_FLOOR — head violation
            0.002, // also < floor and < gap
            0.5,   // ok
            0.51,  // gap < 0.0391 from previous
            1.0,   // ok
            1.02,  // gap < 0.0391
            1.5,   // ok
            2.0,   // ok
            2.5,   // ok
            4.0,   // > CLAMP_CEIL — tail violation
        ];
        stability_clamp(&mut v);

        assert!(v[0] >= CLAMP_FLOOR, "head not floored: {}", v[0]);
        for i in 1..M {
            assert!(
                v[i] - v[i - 1] + 1e-6 >= CLAMP_MIN_GAP,
                "min-gap violated at i={i}: gap {} < {}",
                v[i] - v[i - 1],
                CLAMP_MIN_GAP,
            );
        }
        assert!(
            v[M - 1] <= CLAMP_CEIL + 1e-6,
            "tail not ceiled: {}",
            v[M - 1]
        );
    }

    /// An already-clean LSF vector is left intact by the clamp.
    #[test]
    fn stability_clamp_is_a_noop_on_clean_input() {
        let clean: [f32; M] = [0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00];
        let mut v = clean;
        stability_clamp(&mut v);
        for i in 0..M {
            assert!(
                (v[i] - clean[i]).abs() < 1e-6,
                "no-op clamp disturbed index {i}: {} vs {}",
                v[i],
                clean[i],
            );
        }
    }

    /// End-to-end reconstruction of one frame returns an LSF vector
    /// that satisfies all post-clamp invariants (sorted, in
    /// `[CLAMP_FLOOR, CLAMP_CEIL]`, min adjacent gap honoured).
    #[test]
    fn reconstruct_frame_yields_clamp_compliant_output() {
        let mut r = LspReconstructor::new();
        let omega = r.reconstruct_frame(0, 0, 0, 0).expect("valid indices");
        assert!(omega[0] >= CLAMP_FLOOR, "head violated: {}", omega[0]);
        for i in 1..M {
            assert!(
                omega[i] - omega[i - 1] + 1e-6 >= CLAMP_MIN_GAP,
                "gap at i={i}: {}",
                omega[i] - omega[i - 1],
            );
        }
        assert!(
            omega[M - 1] <= CLAMP_CEIL + 1e-6,
            "tail violated: {}",
            omega[M - 1]
        );
    }

    /// `reconstruct_frame` advances the MA history: after one call
    /// `history[0]` is the just-consumed *post-rearrange* residual,
    /// and the previous `history[0]` has migrated to `history[1]`.
    #[test]
    fn reconstruct_frame_advances_ma_history() {
        let mut r = LspReconstructor::new();
        let pre_h0 = r.history()[0];
        let pre_h1 = r.history()[1];
        let _ = r.reconstruct_frame(0, 0, 0, 0).expect("valid indices");
        // history[1] should equal the pre-call history[0].
        for (i, expected) in pre_h0.iter().enumerate() {
            let got = r.history()[1][i];
            assert!(
                (got - expected).abs() < 1e-6,
                "history[1][{i}] = {got}; expected pre-call history[0] = {expected}",
            );
        }
        // history[2] should equal the pre-call history[1] (same
        // start-up vector — verifies shift, not specific value).
        for (i, expected) in pre_h1.iter().enumerate() {
            let got = r.history()[2][i];
            assert!(
                (got - expected).abs() < 1e-6,
                "history[2][{i}] = {got}; expected pre-call history[1] = {expected}",
            );
        }
    }

    /// Both predictor modes (`L0 ∈ {0, 1}`) accept the same residual
    /// without panicking. The outputs need not coincide — the modes
    /// are the *whole point* of the switched predictor — but both
    /// must produce clamp-compliant LSF vectors.
    #[test]
    fn both_l0_modes_produce_clamp_compliant_output() {
        for mode in 0..2 {
            let mut r = LspReconstructor::new();
            let omega = r.reconstruct_frame(mode, 0, 0, 0).expect("valid indices");
            assert!(
                omega[0] >= CLAMP_FLOOR,
                "mode {mode} head violated: {}",
                omega[0],
            );
            for i in 1..M {
                assert!(
                    omega[i] - omega[i - 1] + 1e-6 >= CLAMP_MIN_GAP,
                    "mode {mode} gap at i={i}: {}",
                    omega[i] - omega[i - 1],
                );
            }
            assert!(
                omega[M - 1] <= CLAMP_CEIL + 1e-6,
                "mode {mode} tail violated: {}",
                omega[M - 1],
            );
        }
    }
}
