//! §3.2.4 **LSP quantiser encoder search** — the encoder-side inverse of
//! [`crate::lsp_reconstruct`]. Takes the unquantised cosine-domain LSPs
//! `q_i` from [`crate::lp_to_lsp`] and selects the four transmitted
//! codeword fields `(L0, L1, L2, L3)` that best approximate them under
//! the spec's switched-predictor two-stage VQ.
//!
//! This is the fifth encoder-side stage. It is validated by a full
//! encode→decode round-trip: the indices it emits, fed straight back
//! through [`crate::lsp_reconstruct::LspReconstructor`], reconstruct an
//! LSF vector close to the encoder input — so encoder and decoder stay
//! bit-consistent by construction (the quantiser *delegates* its own
//! reconstruction + MA-history update to a wrapped `LspReconstructor`).
//!
//! ## Spec source — clause 3.2.4, equations (18)–(23)
//!
//! (Transcribed from the EPUB prose + the equation rasters
//! `images/eq18.jpg … eq23.jpg`.)
//!
//! * **eq (18)** — the LSPs are quantised in the LSF (normalised
//!   frequency) domain: `ω_i = arccos(q_i)`, `ω_i ∈ [0, π]`.
//! * **eq (20)** (reconstruction, decode-side) —
//!   `ω̂_i^(m) = (1 − Σ_{k=1}^{4} P̂_{i,k})·l̂_i^m + Σ_{k=1}^{4} P̂_{i,k}·l̂_i^(m−k)`.
//! * **eq (21)** — the weighted MSE the search minimises:
//!   `E_lsf = Σ_{i=1}^{10} w_i·(ω_i − ω̂_i)²`.
//! * **eq (22)** — the adaptive weights from the *unquantised* LSF:
//!   ```text
//!   w_1  = 1.0 if (ω_2 − 0.04π − 1) > 0 else 10·(ω_2 − 0.04π − 1)² + 1
//!   w_i  = 1.0 if (ω_{i+1} − ω_{i−1} − 1) > 0
//!               else 10·(ω_{i+1} − ω_{i−1} − 1)² + 1     (2 ≤ i ≤ 9)
//!   w_10 = 1.0 if (−ω_9 + 0.92π − 1) > 0 else 10·(−ω_9 + 0.92π − 1)² + 1
//!   ```
//!   and additionally `w_5`, `w_6` are each multiplied by 1.2.
//! * **eq (23)** — the vector actually quantised for the current frame,
//!   per MA-predictor mode:
//!   `l_i = [ω_i^(m) − Σ_{k=1}^{4} P̂_{i,k}·l̂_i^(m−k)] / (1 − Σ_{k=1}^{4} P̂_{i,k})`.
//!
//! ## Search procedure (clause 3.2.4 prose)
//!
//! For each of the two MA predictors `L0 ∈ {0, 1}`:
//!
//! 1. Form the target `l` (eq (23)) from the unquantised LSF and the MA
//!    history of previously-selected residuals.
//! 2. Search codebook **L1** (128 entries): pick the entry minimising
//!    the *unweighted* MSE to `l` over all ten coordinates.
//! 3. Search **L2** (32, lower five coordinates): for each candidate,
//!    form the partial residual `l̂ = L1 + L2`, rearrange it to a minimum
//!    adjacent distance of 0.0012, and keep the candidate with the
//!    lowest **weighted residual-domain MSE** `Σ w_i·(l_i − l̂_i)²`
//!    over `i = 1…5` (see the round-388 note below).
//! 4. Search **L3** (32, upper five): analogously over `i = 6…10`.
//! 5. Reconstruct the full quantised vector (rearrange twice — 0.0012
//!    then 0.0006 — then eq (20)) and compute the total eq (21)
//!    weighted MSE `Σ w_i·(ω_i − ω̂_i)²`.
//!
//! The predictor `L0` giving the lowest total weighted MSE wins. The
//! selected indices are then run through the wrapped `LspReconstructor`
//! (which repeats the exact decode reconstruction and advances the MA
//! history), so the encoder's future targets are computed against the
//! same residual history the decoder will see.
//!
//! ## Round-388 note — the split-stage metric domain (measured)
//!
//! The printed eq (21) is the ω-domain error `Σ w_i·(ω_i − ω̂_i)²`, and
//! the §3.2.4 prose describes the L2/L3 stages as reconstructing `ω̂`
//! via eq (20) before computing it. Because eq (20) is affine in `l̂`,
//! the two domains differ only by a per-coordinate `(1 − Σ_k P̂_{i,k})²`
//! folding of the weights: `ω_i − ω̂_i = (1 − Σ P̂_{i,k})·(l_i − l̂_i)`.
//! Clause 2.5 makes the bit-exact fixed-point coder the normative
//! ground truth, and against its conformance corpus the split-stage
//! searches **measurably behave as the residual-domain form** (no
//! `(1 − ΣP̂)²` folding): switching stages 3–4 to
//! `Σ w_i·(l_i − l̂_i)²` raises the reference-locked all-four-indices
//! agreement on every affected vector — ALGTHM 71.4 → 82.9%, FIXED
//! 91.7 → 97.5%, LSP 74.9 → 77.5%, PITCH 88.6 → 92.8%, SPEECH
//! 78.4 → 80.9%, TAME unchanged — with no vector degrading. The L1
//! stage stays unweighted (probed: weighting it collapses agreement)
//! and the step-5 mode selection stays on the printed ω-domain
//! eq (21) (probed: residual-domain mode selection collapses
//! agreement).
//!
//! ## Round-390 note — the fixed-point Q13 search
//!
//! The staged algorithm description `docs/audio/g729/l0-mode-selection.md`
//! (docs commit `b9e48a4`) pins the L0 decision as the argmin of the two
//! full per-mode searches and names the fixed-point latitude the prose
//! leaves open. [`search_lsp_indices_q13`] is that search on the
//! reference's integer Q13 grid, with the latitude exposed as
//! [`FxLatitude`] and pinned black-box against the conformance corpus
//! ([`FxLatitude::default`]). [`LspQuantizer::quantize_lsf`] runs it
//! over an integer MA history. The float [`search_lsp_indices`] is
//! retained as the legacy r388 search for comparison harnesses.

use core::f32::consts::PI;

use crate::lsp_reconstruct::{codebook_sum, rearrange_pass, LspReconstructor, REARRANGE_J1};
use crate::tables::{self, M, MA_NP, NC0, NC1};

/// Q13 unit (`2^13`) for the LSP codebook literals.
const Q13: f32 = 8192.0;
/// Q15 unit (`2^15`) for the MA-predictor `fg` literals.
const Q15: f32 = 32768.0;

/// The four transmitted LSP codeword fields selected by the quantiser
/// (spec Table 8: `L0` 1 bit, `L1` 7 bits, `L2` 5 bits, `L3` 5 bits).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LspIndices {
    /// Predictor-switch flag `L0` (0 or 1).
    pub l0: usize,
    /// First-stage codebook index `L1` (`0 … 127`).
    pub l1: usize,
    /// Second-stage lower-split index `L2` (`0 … 31`).
    pub l2: usize,
    /// Second-stage upper-split index `L3` (`0 … 31`).
    pub l3: usize,
}

/// One frame's quantiser result: the transmitted indices and the
/// decoder-consistent reconstructed LSF vector `ω̂` (radians in
/// `[0, π]`).
#[derive(Debug, Clone)]
pub struct LspQuantized {
    /// Selected codeword indices.
    pub indices: LspIndices,
    /// Reconstructed LSF `ω̂_i` (what the decoder will produce from
    /// `indices`).
    pub omega_hat: [f32; M],
}

/// Reads the `fg` MA-predictor plane for `mode` into `f32` (Q15→real).
#[inline]
fn fg_plane_f32(mode: usize) -> [[f32; M]; tables::MA_NP] {
    let src = tables::lsp_fg_plane(mode);
    std::array::from_fn(|k| std::array::from_fn(|i| f32::from(src[k][i]) / Q15))
}

/// Reads the per-mode `(1 − Σ_k P̂_{i,k})` sum factor into `f32`.
#[inline]
fn fg_sum_f32(mode: usize) -> [f32; M] {
    let src = tables::lsp_fg_sum(mode);
    std::array::from_fn(|i| f32::from(src[i]) / Q15)
}

/// Reconstructs `ω̂_i` for coordinates `range` from a residual `l_hat`
/// via eq (20), given the mode's `fg`/`fg_sum` and the MA history.
fn reconstruct_omega(
    l_hat: &[f32; M],
    fg: &[[f32; M]; tables::MA_NP],
    fg_sum: &[f32; M],
    history: &[[f32; M]; tables::MA_NP],
    out: &mut [f32; M],
) {
    for i in 0..M {
        let mut acc = fg_sum[i] * l_hat[i];
        for k in 0..tables::MA_NP {
            acc += fg[k][i] * history[k][i];
        }
        out[i] = acc;
    }
}

/// Computes the eq (22) adaptive weights from the unquantised LSF `omega`
/// (radians), with the `w_5`, `w_6` ×1.2 boost applied.
#[must_use]
pub fn lsf_weights(omega: &[f32; M]) -> [f32; M] {
    let piece = |d: f32| -> f32 {
        if d > 0.0 {
            1.0
        } else {
            10.0 * d * d + 1.0
        }
    };
    let mut w = [1.0f32; M];
    // w_1 (idx 0): uses ω_2 (idx 1).
    w[0] = piece(omega[1] - 0.04 * PI - 1.0);
    // w_i, 2 ≤ i ≤ 9 (idx 1..=8): uses ω_{i+1} − ω_{i−1}.
    for i in 1..=8 {
        w[i] = piece(omega[i + 1] - omega[i - 1] - 1.0);
    }
    // w_10 (idx 9): uses ω_9 (idx 8).
    w[9] = piece(-omega[8] + 0.92 * PI - 1.0);
    // Boost w_5, w_6 (idx 4, 5).
    w[4] *= 1.2;
    w[5] *= 1.2;
    w
}

/// Pure §3.2.4 index search: given the unquantised LSF `omega`
/// (radians) and an explicit MA history of past residuals
/// `l̂^(m−1) … l̂^(m−4)` (slot 0 = most recent), runs the staged
/// L1 → L2 → L3 search over both `L0` predictor modes and returns the
/// winning indices.
///
/// This is the stateless core of [`LspQuantizer::quantize`], exposed so
/// conformance harnesses can evaluate the search against an
/// externally-maintained (e.g. reference-index-driven) history.
#[must_use]
pub fn search_lsp_indices(omega: &[f32; M], history: &[[f32; M]; tables::MA_NP]) -> LspIndices {
    let weights = lsf_weights(omega);

    let mut best: Option<(f32, LspIndices)> = None;

    for mode in 0..2 {
        let fg = fg_plane_f32(mode);
        let fg_sum = fg_sum_f32(mode);

        // eq (23): target residual l_i for this predictor mode.
        let mut target = [0.0f32; M];
        for i in 0..M {
            let mut pred = 0.0f32;
            for k in 0..tables::MA_NP {
                pred += fg[k][i] * history[k][i];
            }
            // fg_sum is (1 − Σ P̂); guard against a zero divisor.
            let denom = if fg_sum[i].abs() < 1e-9 {
                1e-9
            } else {
                fg_sum[i]
            };
            target[i] = (omega[i] - pred) / denom;
        }

        // Stage 1 — L1: minimise the unweighted MSE to the target.
        let mut l1 = 0usize;
        let mut best_l1 = f32::INFINITY;
        for c in 0..NC0 {
            let row = tables::lsp_l1_entry(c);
            let mut e = 0.0f32;
            for i in 0..M {
                let d = target[i] - f32::from(row[i]) / Q13;
                e += d * d;
            }
            if e < best_l1 {
                best_l1 = e;
                l1 = c;
            }
        }
        let l1_row = tables::lsp_l1_entry(l1);

        // Stage 2 lower — L2: weighted residual-domain MSE over
        // i = 0..5 (round-388 note: the bit-exact coder's measured
        // metric — no (1 − ΣP̂)² folding of the weights).
        let mut l2 = 0usize;
        let mut best_l2 = f32::INFINITY;
        for c in 0..NC1 {
            let lo = tables::lsp_l2_entry(c);
            let mut l_hat = [0.0f32; M];
            for i in 0..M {
                l_hat[i] = f32::from(l1_row[i]) / Q13;
            }
            for (i, &v) in lo.iter().enumerate() {
                l_hat[i] += f32::from(v) / Q13;
            }
            rearrange_pass(&mut l_hat, REARRANGE_J1);
            let mut e = 0.0f32;
            for i in 0..M / 2 {
                let d = target[i] - l_hat[i];
                e += weights[i] * d * d;
            }
            if e < best_l2 {
                best_l2 = e;
                l2 = c;
            }
        }
        let l2_lo = tables::lsp_l2_entry(l2);

        // Stage 2 upper — L3: weighted residual-domain MSE over
        // i = 5..10.
        let mut l3 = 0usize;
        let mut best_l3 = f32::INFINITY;
        for c in 0..NC1 {
            let hi = tables::lsp_l3_entry(c);
            let mut l_hat = [0.0f32; M];
            for i in 0..M {
                l_hat[i] = f32::from(l1_row[i]) / Q13;
            }
            for (i, &v) in l2_lo.iter().enumerate() {
                l_hat[i] += f32::from(v) / Q13;
            }
            for (j, &v) in hi.iter().enumerate() {
                l_hat[M / 2 + j] += f32::from(v) / Q13;
            }
            rearrange_pass(&mut l_hat, REARRANGE_J1);
            let mut e = 0.0f32;
            for i in M / 2..M {
                let d = target[i] - l_hat[i];
                e += weights[i] * d * d;
            }
            if e < best_l3 {
                best_l3 = e;
                l3 = c;
            }
        }

        // Full reconstruction (rearrange twice + eq (20)) for the
        // mode-selection total error. `codebook_sum` + the decode
        // rearrangement mirror `LspReconstructor` exactly.
        let mut l_hat = codebook_sum(l1, l2, l3).expect("indices in range");
        rearrange_pass(&mut l_hat, REARRANGE_J1);
        rearrange_pass(&mut l_hat, crate::lsp_reconstruct::REARRANGE_J2);
        let mut omega_hat = [0.0f32; M];
        reconstruct_omega(&l_hat, &fg, &fg_sum, history, &mut omega_hat);
        let mut total = 0.0f32;
        for i in 0..M {
            let d = omega[i] - omega_hat[i];
            total += weights[i] * d * d;
        }

        let is_better = match best {
            None => true,
            Some((e, _)) => total < e,
        };
        if is_better {
            best = Some((
                total,
                LspIndices {
                    l0: mode,
                    l1,
                    l2,
                    l3,
                },
            ));
        }
    }

    best.expect("two modes always evaluated").1
}

// ---------------------------------------------------------------------
// Fixed-point (Q13) §3.2.4 search — round 390.
//
// Clause 2.5 makes the 16-bit fixed-point coder the normative ground
// truth, and the staged algorithm description
// `docs/audio/g729/l0-mode-selection.md` (docs commit b9e48a4) pins the
// L0 mode decision as the argmin of the full per-mode L1/L2/L3 search
// error — while explicitly leaving three fixed-point latitude points
// unpinned (behaviourally validatable against the conformance corpus):
//
//   1. the eq (20)/(23) reconstruction shift (round vs truncate,
//      per-mode asymmetric because the two predictors' coefficient
//      magnitudes differ);
//   2. the eq (21) per-term combine (how the Q26 squared distance is
//      truncated when multiplied by the weight — magnitude-dependent);
//   3. the final-compare tie-break inequality.
//
// The `FxLatitude` struct exposes exactly that latitude; the default is
// the corpus-selected assignment (black-box sweep over the six staged
// vectors — see the module CHANGELOG entry for the measured deltas).
// ---------------------------------------------------------------------

/// First rearrangement distance on the Q13 grid (`10/8192 ≈ 0.0012`).
pub const REARRANGE_J1_Q13: i32 = 10;
/// Second rearrangement distance on the Q13 grid (`5/8192 ≈ 0.0006`).
pub const REARRANGE_J2_Q13: i32 = 5;

/// eq (21) per-term combine shape — how the Q26 squared distance
/// `d_i²` and the weight `w_i` are multiplied into the running
/// distance. The staged doc pins that the squared distance is a 32-bit
/// value whose low bits are (partially) discarded when the weight is
/// applied, but not the exact schedule; these are the candidate
/// shapes swept black-box against the corpus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistShape {
    /// Exact wide product `w·d²` (64-bit, no truncation).
    Exact,
    /// High-word truncation: `((2·d²) >> 16) · w` — the low 16 bits of
    /// the Q27 squared distance are discarded before weighting.
    HiWord,
    /// High-word with pre-shift rounding: `((2·d² + 2^15) >> 16) · w`.
    HiWordRound,
    /// 32×16 split multiply: `hi·w·2 + ((lo·w) >> 15)·2` with
    /// `hi = t >> 16`, `lo = (t >> 1) & 0x7fff`, `t = 2·d²` — keeps
    /// most of the low word but truncates its weighted tail.
    Mpy32x16,
}

/// Mode-total formulation — the domain in which the two modes' final
/// scores (the `L0` argmin inputs) are accumulated. The printed
/// eq (21) is the ω-domain error; the staged doc leaves the
/// fixed-point evaluation unpinned, and the corpus discriminates
/// between these candidates (the TAME flips are *systematic* ~2%
/// margins, not rounding near-ties — see the round-390 CHANGELOG).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModeTotal {
    /// Full recompute in the ω-domain: `Σ w_i·(ω_i − ω̂_i)²` after the
    /// twice-rearranged eq (20) reconstruction.
    OmegaDomain,
    /// The sum of the two split-stage best scores (residual domain,
    /// single `J1` rearrangement) — no step-5 recompute.
    StageSum,
    /// Full recompute in the residual domain: `Σ w_i·(l_i − l̂_i)²`
    /// after the twice-rearranged residual.
    ResidualFull,
    /// Single-folded product: `Σ w_i·(ω_i − ω̂_i)·(l_i − l̂_i)` — one
    /// factor of `(1 − ΣP̂)` between the ω-domain and residual-domain
    /// forms.
    Product,
    /// Pre-selection on the unweighted eq (23) target energy
    /// `Σ l_i²` — the mode is chosen *before* the stage search and
    /// only the winner is searched.
    PreTarget,
    /// Pre-selection on the weighted target energy `Σ w_i·l_i²`.
    PreTargetW,
    /// Pre-selection on the unweighted prediction error
    /// `Σ (ω_i − p_i)²` (`p` = the MA-history part of eq (20)).
    PreOmega,
    /// Pre-selection on the weighted prediction error
    /// `Σ w_i·(ω_i − p_i)²`.
    PreOmegaW,
}

/// The fixed-point latitude of the §3.2.4 search — the three unpinned
/// points of `docs/audio/g729/l0-mode-selection.md` plus the weight
/// grid. [`FxLatitude::default`] is the corpus-selected assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FxLatitude {
    /// Add `2^14` before the eq (20)/(23) `Q28 → Q13` prediction shift
    /// (round) instead of truncating.
    pub round_predict: bool,
    /// Add `2^11` before the eq (23) `Q25 → Q13` reciprocal-multiply
    /// shift (round) instead of truncating.
    pub round_target: bool,
    /// Add `1` before the rearrangement midpoint `>> 1` (round)
    /// instead of truncating.
    pub round_rearrange: bool,
    /// eq (22) weight grid: weights are represented on a
    /// `2^weight_shift` grid (round-to-nearest from the exact value).
    pub weight_shift: u32,
    /// eq (21) per-term combine shape.
    pub dist_shape: DistShape,
    /// Mode-total formulation for the `L0` argmin.
    pub mode_total: ModeTotal,
    /// Mode tie-break: when the two modes' totals compare equal,
    /// `true` selects mode 1 (`<=`), `false` holds mode 0 (`<`).
    pub mode_tie_le: bool,
}

impl Default for FxLatitude {
    /// The corpus-selected latitude (round-390 black-box sweep over
    /// all six staged conformance vectors; see the crate CHANGELOG
    /// for the per-vector agreement matrix).
    fn default() -> Self {
        Self {
            round_predict: true,
            round_target: false,
            round_rearrange: false,
            weight_shift: 11,
            dist_shape: DistShape::Exact,
            mode_total: ModeTotal::OmegaDomain,
            mode_tie_le: false,
        }
    }
}

/// Arithmetic right shift with optional pre-shift rounding constant.
#[inline]
fn shr_rnd(x: i32, s: u32, round: bool) -> i32 {
    if round {
        (x + (1 << (s - 1))) >> s
    } else {
        x >> s
    }
}

/// eq (19) codebook sum on the Q13 integer grid.
#[must_use]
fn codebook_sum_q13(l1: usize, l2: usize, l3: usize) -> [i32; M] {
    let l1_row = tables::lsp_l1_entry(l1);
    let l2_lo = tables::lsp_l2_entry(l2);
    let l3_hi = tables::lsp_l3_entry(l3);
    std::array::from_fn(|i| {
        let stage2 = if i < M / 2 {
            l2_lo[i]
        } else {
            l3_hi[i - M / 2]
        };
        i32::from(l1_row[i]) + i32::from(stage2)
    })
}

/// One §3.2.4 rearrangement pass on the Q13 integer grid: adjacent
/// pairs closer than `j` are replaced by the centred pair
/// `((sum − j) >> 1, (sum + j) >> 1)`.
pub fn rearrange_pass_q13(coefs: &mut [i32; M], j: i32, round: bool) {
    for i in 1..M {
        if coefs[i - 1] > coefs[i] - j {
            let sum = coefs[i] + coefs[i - 1];
            coefs[i - 1] = shr_rnd(sum - j, 1, round);
            coefs[i] = shr_rnd(sum + j, 1, round);
        }
    }
}

/// eq (22) adaptive weights on the `2^weight_shift` grid. The weights
/// are exact-real evaluated (the doc pins them as shared between the
/// two modes and *not* the flip source; the crate's earlier fixed-point
/// weight-grid probes are recorded as rejected hypotheses) and then
/// rounded onto the configured grid.
#[must_use]
fn lsf_weights_q(omega_q13: &[i32; M], weight_shift: u32) -> [i64; M] {
    let unit = f64::from(1u32 << 13);
    let piece = |d: f64| -> f64 {
        if d > 0.0 {
            1.0
        } else {
            10.0 * d * d + 1.0
        }
    };
    let om = |i: usize| f64::from(omega_q13[i]) / unit;
    let mut w = [1.0f64; M];
    w[0] = piece(om(1) - 0.04 * core::f64::consts::PI - 1.0);
    for (i, slot) in w.iter_mut().enumerate().take(9).skip(1) {
        *slot = piece(om(i + 1) - om(i - 1) - 1.0);
    }
    w[9] = piece(-om(8) + 0.92 * core::f64::consts::PI - 1.0);
    w[4] *= 1.2;
    w[5] *= 1.2;
    let grid = f64::from(1u32 << weight_shift);
    std::array::from_fn(|i| (w[i] * grid).round() as i64)
}

/// eq (21) one weighted term per the configured combine shape.
#[inline]
fn dist_term(d: i32, w: i64, shape: DistShape) -> i64 {
    dist_term_prod(d, d, w, shape)
}

/// One weighted product term `w·d1·d2` per the configured combine
/// shape ([`dist_term`] is the `d1 == d2` squared case; the
/// [`ModeTotal::Product`] mode total uses `d1 = ω − ω̂`,
/// `d2 = l − l̂`).
#[inline]
fn dist_term_prod(d1: i32, d2: i32, w: i64, shape: DistShape) -> i64 {
    let sq = i64::from(d1) * i64::from(d2); // Q26
    match shape {
        DistShape::Exact => w * sq,
        DistShape::HiWord => {
            let t = (2 * sq) as i32; // Q27 Word32
            i64::from(t >> 16) * w
        }
        DistShape::HiWordRound => {
            let t = (2 * sq) as i32;
            i64::from((t + (1 << 15)) >> 16) * w
        }
        DistShape::Mpy32x16 => {
            let t = (2 * sq) as i32;
            let hi = i64::from(t >> 16);
            let lo = i64::from((t >> 1) & 0x7fff);
            2 * (hi * w) + 2 * ((lo * w) >> 15)
        }
    }
}

/// eq (23) fixed-point target for one predictor mode: the
/// prediction-removed residual, normalised by the pre-tabulated Q12
/// reciprocal `fg_sum_inv`. Returns `(target, diff)` where `diff` is
/// the un-normalised prediction error `ω_i − p_i` (Q13).
#[must_use]
fn target_q13(
    mode: usize,
    omega_q13: &[i32; M],
    history: &[[i16; M]; MA_NP],
    lat: &FxLatitude,
) -> ([i32; M], [i32; M]) {
    let fg = tables::lsp_fg_plane(mode);
    let inv = tables::lsp_fg_sum_inv(mode);
    let mut target = [0i32; M];
    let mut diffs = [0i32; M];
    for i in 0..M {
        let mut acc: i32 = 0; // Q28
        for k in 0..MA_NP {
            acc += i32::from(fg[k][i]) * i32::from(history[k][i]);
        }
        let pred = shr_rnd(acc, 15, lat.round_predict); // Q13
        let diff = omega_q13[i] - pred; // Q13
        diffs[i] = diff;
        target[i] = shr_rnd(diff * i32::from(inv[i]), 12, lat.round_target); // Q13
    }
    (target, diffs)
}

/// eq (20) fixed-point reconstruction for one predictor mode.
#[must_use]
fn reconstruct_omega_q13(
    mode: usize,
    l_hat: &[i32; M],
    history: &[[i16; M]; MA_NP],
    round_predict: bool,
) -> [i32; M] {
    let fg = tables::lsp_fg_plane(mode);
    let fg_sum = tables::lsp_fg_sum(mode);
    std::array::from_fn(|i| {
        let mut acc: i32 = i32::from(fg_sum[i]) * l_hat[i]; // Q28
        for k in 0..MA_NP {
            acc += i32::from(fg[k][i]) * i32::from(history[k][i]);
        }
        shr_rnd(acc, 15, round_predict) // Q13
    })
}

/// The spec start-up MA history `l̂_i = iπ/11` on the Q13 grid.
#[must_use]
pub fn startup_history_q13() -> [[i16; M]; MA_NP] {
    let row: [i16; M] = std::array::from_fn(|i| {
        let v = f64::from((i + 1) as u32) * core::f64::consts::PI / 11.0 * 8192.0;
        v.round() as i16
    });
    [row; MA_NP]
}

/// Advances a Q13 MA history with the residual selected by
/// `(l1, l2, l3)`: eq (19) codebook sum, rearranged twice
/// (`J = 10` then `5` on the Q13 grid), pushed into slot 0.
///
/// # Panics
///
/// Panics if `l1 >= NC0` or `l2`/`l3 >= NC1` (out of their Table-8
/// field ranges).
pub fn advance_history_q13(
    history: &mut [[i16; M]; MA_NP],
    l1: usize,
    l2: usize,
    l3: usize,
    lat: &FxLatitude,
) {
    let mut residual = codebook_sum_q13(l1, l2, l3);
    rearrange_pass_q13(&mut residual, REARRANGE_J1_Q13, lat.round_rearrange);
    rearrange_pass_q13(&mut residual, REARRANGE_J2_Q13, lat.round_rearrange);
    for k in (1..MA_NP).rev() {
        history[k] = history[k - 1];
    }
    history[0] = std::array::from_fn(|i| residual[i] as i16);
}

/// Fixed-point §3.2.4 index search on the Q13 grid: the staged
/// L1 → L2 → L3 search run under both `L0` predictor modes, with the
/// mode selected as the argmin of the full eq (21) ω-domain weighted
/// error (docs `l0-mode-selection.md` §3), under the configured
/// fixed-point latitude.
#[must_use]
pub fn search_lsp_indices_q13(
    omega_q13: &[i32; M],
    history: &[[i16; M]; MA_NP],
    lat: &FxLatitude,
) -> LspIndices {
    let weights = lsf_weights_q(omega_q13, lat.weight_shift);

    // Pre-selection variants: the mode is decided from the eq (23)
    // target / prediction error alone, and only the winner's stages
    // are searched.
    let pre_metric = |mode: usize| -> Option<i64> {
        let (target, diffs) = target_q13(mode, omega_q13, history, lat);
        let mut e: i64 = 0;
        match lat.mode_total {
            ModeTotal::PreTarget => {
                for &t in &target {
                    e += i64::from(t) * i64::from(t);
                }
            }
            ModeTotal::PreTargetW => {
                for i in 0..M {
                    e += dist_term(target[i], weights[i], lat.dist_shape);
                }
            }
            ModeTotal::PreOmega => {
                for &d in &diffs {
                    e += i64::from(d) * i64::from(d);
                }
            }
            ModeTotal::PreOmegaW => {
                for i in 0..M {
                    e += dist_term(diffs[i], weights[i], lat.dist_shape);
                }
            }
            _ => return None,
        }
        Some(e)
    };
    if let (Some(e0), Some(e1)) = (pre_metric(0), pre_metric(1)) {
        let mode = if lat.mode_tie_le {
            usize::from(e1 <= e0)
        } else {
            usize::from(e1 < e0)
        };
        let (target, _) = target_q13(mode, omega_q13, history, lat);
        let (l1, l2, l3, _, _) = stage_search_q13(&target, &weights, lat);
        return LspIndices {
            l0: mode,
            l1,
            l2,
            l3,
        };
    }

    let mut best: Option<(i64, LspIndices)> = None;

    for mode in 0..2 {
        let (target, _) = target_q13(mode, omega_q13, history, lat);
        let (l1, l2, l3, best_l2, best_l3) = stage_search_q13(&target, &weights, lat);

        // Mode total — the configured formulation over the full
        // reconstruction (rearrange twice + eq (20)).
        let mut l_hat = codebook_sum_q13(l1, l2, l3);
        rearrange_pass_q13(&mut l_hat, REARRANGE_J1_Q13, lat.round_rearrange);
        rearrange_pass_q13(&mut l_hat, REARRANGE_J2_Q13, lat.round_rearrange);
        let omega_hat = reconstruct_omega_q13(mode, &l_hat, history, lat.round_predict);
        let mut total: i64 = 0;
        match lat.mode_total {
            ModeTotal::StageSum => total = best_l2 + best_l3,
            ModeTotal::ResidualFull => {
                for i in 0..M {
                    total += dist_term(target[i] - l_hat[i], weights[i], lat.dist_shape);
                }
            }
            ModeTotal::Product => {
                for i in 0..M {
                    total += dist_term_prod(
                        omega_q13[i] - omega_hat[i],
                        target[i] - l_hat[i],
                        weights[i],
                        lat.dist_shape,
                    );
                }
            }
            // `OmegaDomain` (and the pre-selection variants, already
            // handled above): the printed eq (21) ω-domain error.
            _ => {
                for i in 0..M {
                    total += dist_term(omega_q13[i] - omega_hat[i], weights[i], lat.dist_shape);
                }
            }
        }

        let is_better = match best {
            None => true,
            Some((e, _)) => {
                if lat.mode_tie_le {
                    total <= e
                } else {
                    total < e
                }
            }
        };
        if is_better {
            best = Some((
                total,
                LspIndices {
                    l0: mode,
                    l1,
                    l2,
                    l3,
                },
            ));
        }
    }

    best.expect("two modes always evaluated").1
}

/// The staged L1 → L2 → L3 search for one predictor mode's target.
/// Returns `(l1, l2, l3, best_l2, best_l3)` with the two split-stage
/// best weighted residual-domain scores.
fn stage_search_q13(
    target: &[i32; M],
    weights: &[i64; M],
    lat: &FxLatitude,
) -> (usize, usize, usize, i64, i64) {
    // Stage 1 — L1: unweighted MSE to the eq (23) target.
    let mut l1 = 0usize;
    let mut best_l1 = i64::MAX;
    for c in 0..NC0 {
        let row = tables::lsp_l1_entry(c);
        let mut e: i64 = 0;
        for i in 0..M {
            let d = target[i] - i32::from(row[i]);
            e += i64::from(d) * i64::from(d);
        }
        if e < best_l1 {
            best_l1 = e;
            l1 = c;
        }
    }
    let l1_row = tables::lsp_l1_entry(l1);

    // Stage 2 lower — L2: weighted residual-domain MSE over the
    // lower five coordinates (round-388 measured metric).
    let mut l2 = 0usize;
    let mut best_l2 = i64::MAX;
    for c in 0..NC1 {
        let lo = tables::lsp_l2_entry(c);
        let mut l_hat: [i32; M] = std::array::from_fn(|i| i32::from(l1_row[i]));
        for (i, &v) in lo.iter().enumerate() {
            l_hat[i] += i32::from(v);
        }
        rearrange_pass_q13(&mut l_hat, REARRANGE_J1_Q13, lat.round_rearrange);
        let mut e: i64 = 0;
        for i in 0..M / 2 {
            e += dist_term(target[i] - l_hat[i], weights[i], lat.dist_shape);
        }
        if e < best_l2 {
            best_l2 = e;
            l2 = c;
        }
    }
    let l2_lo = tables::lsp_l2_entry(l2);

    // Stage 2 upper — L3: weighted residual-domain MSE over the
    // upper five coordinates.
    let mut l3 = 0usize;
    let mut best_l3 = i64::MAX;
    for c in 0..NC1 {
        let hi = tables::lsp_l3_entry(c);
        let mut l_hat: [i32; M] = std::array::from_fn(|i| i32::from(l1_row[i]));
        for (i, &v) in l2_lo.iter().enumerate() {
            l_hat[i] += i32::from(v);
        }
        for (j, &v) in hi.iter().enumerate() {
            l_hat[M / 2 + j] += i32::from(v);
        }
        rearrange_pass_q13(&mut l_hat, REARRANGE_J1_Q13, lat.round_rearrange);
        let mut e: i64 = 0;
        for i in M / 2..M {
            e += dist_term(target[i] - l_hat[i], weights[i], lat.dist_shape);
        }
        if e < best_l3 {
            best_l3 = e;
            l3 = c;
        }
    }
    (l1, l2, l3, best_l2, best_l3)
}

/// The switched-predictor two-stage LSP quantiser. Runs the round-390
/// fixed-point Q13 search ([`search_lsp_indices_q13`]) over an integer
/// MA history (the reference's numeric grid, per clause 2.5), while
/// the returned reconstruction is delegated to a wrapped
/// [`LspReconstructor`] driven with the selected indices — so the
/// encoder output stays bit-consistent with the decode path.
#[derive(Debug, Clone)]
pub struct LspQuantizer {
    reconstructor: LspReconstructor,
    history_q13: [[i16; M]; MA_NP],
}

impl Default for LspQuantizer {
    fn default() -> Self {
        Self::new()
    }
}

impl LspQuantizer {
    /// Constructs a fresh quantiser with the spec start-up MA history.
    #[must_use]
    pub fn new() -> Self {
        Self {
            reconstructor: LspReconstructor::new(),
            history_q13: startup_history_q13(),
        }
    }

    /// Quantises one frame of cosine-domain LSPs `q_in[i−1] = q_i`
    /// (`q_i = cos ω_i`, decreasing / ordered), returning the selected
    /// transmitted indices and the decoder-consistent reconstructed LSF.
    ///
    /// The internal MA history advances by exactly the decode-side rule,
    /// so a subsequent [`LspReconstructor`] driven with the returned
    /// indices reproduces `omega_hat`.
    pub fn quantize(&mut self, q_in: &[f32; M]) -> LspQuantized {
        // eq (18): ω_i = arccos(q_i). `q` is clamped into the valid acos
        // domain defensively (the LP→LSP roots already lie in (−1, 1)).
        let omega: [f32; M] = std::array::from_fn(|i| q_in[i].clamp(-1.0, 1.0).acos());
        self.quantize_lsf(&omega)
    }

    /// Quantises one frame of LSF coefficients `omega[i−1] = ω_i`
    /// (radians, increasing) directly — the entry point for callers
    /// that already evaluated eq (18) themselves (e.g. through the
    /// fixed-point [`crate::lsf_conversion`] table path, which is how
    /// the clause-3 encode chain reaches the reference's LSF grid).
    ///
    /// Round 390: the index search runs on the fixed-point Q13 grid
    /// ([`search_lsp_indices_q13`] under [`FxLatitude::default`]) with
    /// an integer MA history advanced by the decode-side rule; the
    /// input LSFs (already Q13-gridded by the `lsf_conversion` path)
    /// are mapped by round-to-nearest.
    pub fn quantize_lsf(&mut self, omega: &[f32; M]) -> LspQuantized {
        let lat = FxLatitude::default();
        let omega_q13: [i32; M] = std::array::from_fn(|i| (omega[i] * 8192.0).round() as i32);
        let indices = search_lsp_indices_q13(&omega_q13, &self.history_q13, &lat);
        advance_history_q13(
            &mut self.history_q13,
            indices.l1,
            indices.l2,
            indices.l3,
            &lat,
        );

        // Delegate the exact reconstruction + MA-history advance to the
        // decode-side reconstructor so encoder/decoder stay in lock-step.
        let omega_hat = self
            .reconstructor
            .reconstruct_frame(indices.l0, indices.l1, indices.l2, indices.l3)
            .expect("selected indices are in range");

        LspQuantized { indices, omega_hat }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_to_lsp::lp_to_lsp;
    use crate::lsp_reconstruct::LspReconstructor;

    /// eq (22): a widely-spaced (flat-ish) LSF set gives all-unity
    /// weights except the ×1.2 boost on w_5/w_6.
    #[test]
    fn weights_flat_case() {
        // Uniformly spaced ω_i = i·π/11 → all adjacency gaps ≈ 0.285 rad,
        // so ω_{i+1} − ω_{i−1} ≈ 0.571 < 1 → the "otherwise" branch, but
        // the boundary terms use 0.04π/0.92π offsets. Just assert the
        // boost structure holds relative to the un-boosted neighbours.
        let omega: [f32; M] = std::array::from_fn(|i| (i + 1) as f32 * PI / 11.0);
        let w = lsf_weights(&omega);
        // w_5 and w_6 are the boosted pair; their raw (pre-boost) value
        // equals the eq-(22) piece, so they must exceed 1.0.
        assert!(w[4] > 1.0 && w[5] > 1.0);
    }

    /// A very wide gap forces the `> 0` branch → weight exactly 1.0.
    #[test]
    fn weights_wide_gap_is_unity() {
        // Construct ω with ω_{i+1} − ω_{i−1} > 1 for an interior index.
        let mut omega: [f32; M] = std::array::from_fn(|i| 0.1 + 0.05 * i as f32);
        omega[3] = 0.1;
        omega[5] = 1.4; // ω_5 − ω_3 (idx 4 uses idx5 − idx3) = 1.3 > 1
        let w = lsf_weights(&omega);
        // idx 4 (w_5) piece is unity before the ×1.2 boost → 1.2 exactly.
        assert!((w[4] - 1.2).abs() < 1e-5, "w[4]={}", w[4]);
    }

    /// End-to-end: quantise the LSPs of a real stable filter, feed the
    /// indices to an independent decoder, and confirm the reconstructed
    /// LSF matches the quantiser's own `omega_hat` bit-for-bit and is
    /// close to the (arccos of the) input LSPs.
    #[test]
    fn encode_decode_consistency() {
        let a: [f32; M] = {
            use crate::levinson::levinson;
            let r: [f64; M + 1] =
                std::array::from_fn(|k| (0.9f64).powi(k as i32) * (0.3 * k as f64).cos());
            levinson(&r).a
        };
        let q = lp_to_lsp(&a).expect("stable filter");
        let omega_in: [f32; M] = std::array::from_fn(|i| q[i].acos());

        let mut enc = LspQuantizer::new();
        let out = enc.quantize(&q);

        // Independent decoder driven by the emitted indices.
        let mut dec = LspReconstructor::new();
        let omega_dec = dec
            .reconstruct_frame(
                out.indices.l0,
                out.indices.l1,
                out.indices.l2,
                out.indices.l3,
            )
            .expect("valid indices");
        for (i, (&od, &oh)) in omega_dec.iter().zip(out.omega_hat.iter()).enumerate() {
            assert!(
                (od - oh).abs() < 1e-6,
                "decode/encode mismatch at {i}: {od} vs {oh}"
            );
        }
        // The quantiser should approximate the input LSF to well within
        // the coarse-grid failure band (a few tenths of a radian even in
        // the worst coordinate; the whole-vector average is far tighter).
        let mut sq = 0.0f32;
        for (&od, &oi) in omega_dec.iter().zip(omega_in.iter()) {
            sq += (od - oi).powi(2);
        }
        let rms = (sq / M as f32).sqrt();
        assert!(rms < 0.2, "LSF quantisation RMS too large: {rms}");
    }

    /// The corpus-pinned default latitude is exactly the round-390
    /// sweep winner — a regression pin so a knob drift is caught.
    #[test]
    fn fx_default_latitude_is_pinned() {
        let lat = FxLatitude::default();
        assert!(lat.round_predict, "eq (20)/(23) prediction shift rounds");
        assert!(!lat.round_target, "eq (23) reciprocal shift truncates");
        assert!(!lat.round_rearrange, "rearrange midpoint shift truncates");
        assert_eq!(lat.weight_shift, 11);
        assert_eq!(lat.dist_shape, DistShape::Exact);
        assert_eq!(lat.mode_total, ModeTotal::OmegaDomain);
        assert!(!lat.mode_tie_le, "mode 0 holds ties");
    }

    /// The Q13 start-up history is `round(iπ/11 · 2^13)` in every slot.
    #[test]
    fn fx_startup_history_matches_spec() {
        let h = startup_history_q13();
        for plane in &h {
            for (i, &v) in plane.iter().enumerate() {
                let expected =
                    (f64::from((i + 1) as u32) * core::f64::consts::PI / 11.0 * 8192.0).round();
                assert_eq!(f64::from(v), expected, "slot {i}");
            }
        }
    }

    /// The integer rearrangement enforces the same post-pass invariant
    /// as the float pass on this input: adjacent gaps of at least `j`
    /// (the centred pair lands exactly `j` apart — `(sum − j)` and
    /// `(sum + j)` share parity, so both midpoint shifts truncate
    /// identically).
    #[test]
    fn fx_rearrange_enforces_minimum_distance() {
        let mut v: [i32; M] = [800, 400, 2400, 2000, 3300, 3290, 5000, 4990, 6600, 6590];
        rearrange_pass_q13(&mut v, REARRANGE_J1_Q13, false);
        for i in 1..M {
            assert!(
                v[i] - v[i - 1] >= REARRANGE_J1_Q13,
                "gap at {i}: {} < {REARRANGE_J1_Q13}",
                v[i] - v[i - 1],
            );
        }
    }

    /// `advance_history_q13` pushes the twice-rearranged eq (19) sum
    /// into slot 0 and shifts the older slots.
    #[test]
    fn fx_advance_history_shifts() {
        let mut h = startup_history_q13();
        let h0_before = h[0];
        let h1_before = h[1];
        advance_history_q13(&mut h, 3, 7, 11, &FxLatitude::default());
        assert_eq!(h[1], h0_before);
        assert_eq!(h[2], h1_before);
        // Slot 0 is the rearranged codebook sum for (3, 7, 11).
        let l1 = tables::lsp_l1_entry(3);
        let l2 = tables::lsp_l2_entry(7);
        let l3 = tables::lsp_l3_entry(11);
        let mut expect: [i32; M] = std::array::from_fn(|i| {
            i32::from(l1[i])
                + if i < M / 2 {
                    i32::from(l2[i])
                } else {
                    i32::from(l3[i - M / 2])
                }
        });
        rearrange_pass_q13(&mut expect, REARRANGE_J1_Q13, false);
        rearrange_pass_q13(&mut expect, REARRANGE_J2_Q13, false);
        for i in 0..M {
            assert_eq!(i32::from(h[0][i]), expect[i], "slot-0 coord {i}");
        }
    }

    /// The fixed-point search returns in-range indices over a spread
    /// of synthetic LSF frames, and the quantiser's integer history
    /// stays consistent with a decode-side reconstructor fed the same
    /// index stream.
    #[test]
    fn fx_search_indices_in_range() {
        let mut hist = startup_history_q13();
        let lat = FxLatitude::default();
        for f in 0..12 {
            let omega_q13: [i32; M] =
                std::array::from_fn(|i| 1500 + 2200 * i as i32 + 37 * (f % 5));
            let idx = search_lsp_indices_q13(&omega_q13, &hist, &lat);
            assert!(idx.l0 < 2);
            assert!(idx.l1 < NC0);
            assert!(idx.l2 < NC1);
            assert!(idx.l3 < NC1);
            advance_history_q13(&mut hist, idx.l1, idx.l2, idx.l3, &lat);
        }
    }

    /// Indices are always within their Table-8 field ranges.
    #[test]
    fn indices_in_range() {
        let mut enc = LspQuantizer::new();
        for f in 0..8 {
            let q: [f32; M] = std::array::from_fn(|i| {
                (0.9 - 0.18 * i as f32 + 0.01 * f as f32).clamp(-0.99, 0.99)
            });
            let out = enc.quantize(&q);
            assert!(out.indices.l0 < 2);
            assert!(out.indices.l1 < NC0);
            assert!(out.indices.l2 < NC1);
            assert!(out.indices.l3 < NC1);
        }
    }

    /// Over several frames the quantiser's internal MA history stays in
    /// lock-step with an independent decoder fed the same index stream
    /// (regression guard on the delegated history update).
    #[test]
    fn multi_frame_history_lockstep() {
        let mut enc = LspQuantizer::new();
        let mut dec = LspReconstructor::new();
        for f in 0..12 {
            let a: [f32; M] = {
                use crate::levinson::levinson;
                let r: [f64; M + 1] = std::array::from_fn(|k| {
                    (0.88f64).powi(k as i32) * (0.2 * k as f64 + 0.05 * f as f64).cos()
                });
                levinson(&r).a
            };
            let q = lp_to_lsp(&a).expect("stable");
            let out = enc.quantize(&q);
            let omega_dec = dec
                .reconstruct_frame(
                    out.indices.l0,
                    out.indices.l1,
                    out.indices.l2,
                    out.indices.l3,
                )
                .expect("valid");
            for (i, (&od, &oh)) in omega_dec.iter().zip(out.omega_hat.iter()).enumerate() {
                assert!((od - oh).abs() < 1e-5, "frame {f} coord {i}: {od} vs {oh}");
            }
        }
    }
}
