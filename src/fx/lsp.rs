//! §4.1.1 fixed-point LP-parameter decode: LSF reconstruction on the
//! Q13 grid, the §3.2.4 stability procedure, the LSF → LSP cosine
//! conversion, the §3.2.5 subframe interpolation, and the §3.2.6
//! LSP → LP conversion to Q12 coefficients.
//!
//! ## Number grids
//!
//! - **LSF** `ω̂_i` — radians on the Q13 grid (`π ≈ 25736`); the
//!   quantiser codebooks are staged in this grid (L1 max entry 23533)
//!   and the round-390 encoder search already runs on it
//!   ([`crate::lsp_quantize::reconstruct_lsf_q13`] is the shared
//!   eq (19)/(20) integer arithmetic).
//! - **LSP** `q̂_i = cos(ω̂_i)` — Q15. The conversion uses the staged
//!   64-entry cosine table + Q19 slopes
//!   ([`crate::tables::LSF_LSP_COS_TABLE2_Q15`] /
//!   [`crate::tables::LSF_LSP_COS_SLOPE_Q19`]): the Q13 radian is
//!   first scaled by `2/π` (Q15 literal 20861) onto a `[0, 16384)`
//!   grid whose top 6 bits index the table and whose low 8 bits
//!   linearly interpolate with the slope.
//! - **LP coefficients** `â_i` — Q12 (`a_0 = 4096`), from the
//!   eq (25)/(26) sum/difference polynomials evaluated in the
//!   double-precision (hi, lo) format of [`crate::fx::dsp`], Q24
//!   accumulators.
//!
//! ## Stability constants (§3.2.4 steps 1–4)
//!
//! The prose gives the stability procedure in rounded decimals
//! (floor 0.005, minimum gap 0.0391, ceiling 3.135 — clause 2.5 marks
//! all prose decimals as rounded versions of the 16-bit values). On
//! the Q13 grid these are pinned as `41` (0.005005), `320` (0.039063,
//! the only neighbour that reprints as 0.0391) and `25681` (3.13489).
//!
//! ## Erasure behaviour (§4.4.1)
//!
//! An erased frame reuses the previous frame's LSP vector. The MA
//! predictor memory advance for an erased frame is not pinned by the
//! prose; this implementation re-pushes the last good residual (so a
//! following good frame predicts from a self-consistent history).

use crate::fx::dsp::{l_extract, mpy_32_16};
use crate::fx::ops::{
    add, extract_l, l_add, l_msu, l_mult, l_shl, l_shr, l_shr_r, l_sub, mult, shr,
};
use crate::lsp_quantize::{reconstruct_lsf_q13, startup_history_q13, FxLatitude};
use crate::tables::{LSF_LSP_COS_SLOPE_Q19, LSF_LSP_COS_TABLE2_Q15, M, MA_NP};

/// §3.2.4 stability floor on the Q13 radian grid (`41/8192 ≈ 0.005`).
pub const STABILITY_FLOOR_Q13: i16 = 41;
/// §3.2.4 stability minimum adjacent gap (`320/8192 ≈ 0.0391`).
pub const STABILITY_GAP_Q13: i16 = 320;
/// §3.2.4 stability ceiling (`25681/8192 ≈ 3.135`).
pub const STABILITY_CEIL_Q13: i16 = 25681;

/// `2/π` on the Q15 grid — scales a Q13 radian onto the 64-segment
/// cosine-table domain (`π ↦ 16384`).
const TWO_OVER_PI_Q15: i16 = 20861;

/// Table 9 start-up value of the §3.2.5 / §4.1.6 previous-frame LSP
/// memory (Q15 cosine domain), **corpus-pinned**.
///
/// The printed Table 9 row (`q_i = arccos(iπ/11)`) is domain-invalid,
/// and its textual resolution `cos(iπ/11)` (exact Q15 grid
/// `31441, 27566, 21458, 13612, 4663, …`, see
/// [`STARTUP_LSP_COS_GRID_Q15`]) does **not** reproduce the conformance
/// corpus: inverting the §4.2.5 output high-pass on the frame-0 `.PST`
/// samples of all six clean vectors (both corpora, whose two different
/// postfilters agree sample-for-sample on that frame) recovers the
/// reference's frame-0 `ŝ(n)` — `[1, 2, 1, 1, 1, 0, …]` on FIXED,
/// `[1, 2, 2, 2, 1, 0, …]` on ALGTHM/PITCH, `[1, 1, 0, …]` on TAME,
/// `[1, 1, 1, 1, 0, …]` on LSP, `[1, 2, 1, 1, 0, −1, …]` on SPEECH —
/// and the cosine grid lands the wrong sample on four of them, while
/// this coarse decreasing vector (a spectrum that is *not* flat: its
/// interpolated frame-0 LP has `a_1 ≈ −0.6`) reproduces all six
/// exactly. Round 452 (`tests/fx_full_conformance.rs`).
pub const STARTUP_LSP_Q15: [i16; M] = [
    30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000,
];

/// The flat-spectrum `q_i = cos(iπ/11)` reading of the Table 9 row on
/// the exact Q15 grid `round(cos(iπ/11)·2^15)` — the textual
/// resolution of the printed `arccos` erratum (§3.2.3 defines
/// `q_i = cos(ω_i)` and the row above initialises `ω̂_i = iπ/11`).
/// Kept as the spec-equation reference; the decoder starts from
/// [`STARTUP_LSP_Q15`] instead (see there).
pub const STARTUP_LSP_COS_GRID_Q15: [i16; M] = [
    31441, 27566, 21458, 13612, 4663, -4663, -13612, -21458, -27566, -31441,
];

/// §3.2.4 steps 1–4 stability procedure on the Q13 grid, in place.
pub fn stability_q13(lsf: &mut [i16; M]) {
    // Step 1 — ascending order.
    lsf.sort_unstable();
    // Step 2 — floor.
    if lsf[0] < STABILITY_FLOOR_Q13 {
        lsf[0] = STABILITY_FLOOR_Q13;
    }
    // Step 3 — minimum adjacent gap.
    for i in 1..M {
        let min = add(lsf[i - 1], STABILITY_GAP_Q13);
        if lsf[i] < min {
            lsf[i] = min;
        }
    }
    // Step 4 — ceiling.
    if lsf[M - 1] > STABILITY_CEIL_Q13 {
        lsf[M - 1] = STABILITY_CEIL_Q13;
    }
}

/// LSF (Q13 radians) → LSP (`cos ω̂`, Q15) via the staged 64-segment
/// cosine table with Q19 slope interpolation.
#[must_use]
pub fn lsf_to_lsp_q15(lsf: &[i16; M]) -> [i16; M] {
    std::array::from_fn(|i| {
        // Scale onto the [0, 16384) table domain.
        let f = mult(lsf[i], TWO_OVER_PI_Q15);
        let ind = usize::from(shr(f, 8).clamp(0, 63) as u16);
        let offset = i32::from(f) & 0xff;
        // table2[ind] + (slope[ind] · offset) — slope Q19 × offset Q8
        // is Q27 (doubled to Q28 by l_mult), landing on Q15 after the
        // 13-bit shift.
        let l_tmp = l_shr(l_mult(LSF_LSP_COS_SLOPE_Q19[ind], offset as i16), 13);
        add(LSF_LSP_COS_TABLE2_Q15[ind], extract_l(l_tmp))
    })
}

/// §3.2.5 subframe-1 interpolation in the LSP (cosine) domain:
/// `q̂_sub1 = (q̂_prev + q̂_curr) / 2` with truncating halves.
#[must_use]
pub fn interpolate_lsp_q15(prev: &[i16; M], curr: &[i16; M]) -> [i16; M] {
    std::array::from_fn(|i| add(shr(prev[i], 1), shr(curr[i], 1)))
}

/// One eq (25) polynomial recursion (spec §3.2.6 figure f0015-01)
/// over five Q15 LSP values, Q24 half-coefficients `f(0..=5)`:
///
/// ```text
///     f(i)  := −2·q·f(i−1) + 2·f(i−2)              (f(−1) = 0)
///     f(j)  := f(j) − 2·q·f(j−1) + f(j−2)          j = i−1 … 1
/// ```
///
/// This is the same recursion the float module
/// ([`crate::lsp_to_lp::lsp_to_lp`]) transcribes from the spec
/// diagram, expressed on the Table 11 operator grid: the `2·q·f(·)`
/// products run in the double-precision (hi, lo) format.
fn lsp_polynomial_q24(lsp: &[i16; 5]) -> [i32; 6] {
    // 2·q·x on the Q24 grid for a Q24 x and Q15 q.
    let mul_2q = |x: i32, q: i16| -> i32 {
        let (hi, lo) = l_extract(x);
        l_shl(mpy_32_16(hi, lo, q), 1)
    };
    let mut f = [0i32; 6];
    f[0] = 0x0100_0000; // 1.0 on the Q24 grid
    f[1] = l_msu(0, lsp[0], 512); // −2·q_1 (q·2^10 on Q24)
    for i in 2..=5 {
        let q = lsp[i - 1];
        // Leading coefficient: −2·q·f(i−1) + 2·f(i−2).
        let lead = l_sub(l_shl(f[i - 2], 1), mul_2q(f[i - 1], q));
        // Inner update, descending so old f(j−1)/f(j−2) are read.
        for j in (1..i).rev() {
            let jm2 = if j >= 2 { f[j - 2] } else { 0 };
            f[j] = l_add(l_sub(f[j], mul_2q(f[j - 1], q)), jm2);
        }
        f[i] = lead;
    }
    f
}

/// §3.2.6 LSP → LP conversion: Q15 LSP vector to Q12 LP coefficients
/// `[a_0 = 4096, a_1 … a_10]` via the eq (25) sum/difference
/// polynomials and the eq (26) combination.
#[must_use]
pub fn lsp_to_lp_q12(lsp: &[i16; M]) -> [i16; M + 1] {
    let even: [i16; 5] = std::array::from_fn(|k| lsp[2 * k]);
    let odd: [i16; 5] = std::array::from_fn(|k| lsp[2 * k + 1]);
    let mut f1 = lsp_polynomial_q24(&even);
    let mut f2 = lsp_polynomial_q24(&odd);

    // F1'(z)·(1 + z^-1) and F2'(z)·(1 − z^-1), in place, descending so
    // each step reads the not-yet-updated lower coefficient.
    for i in (1..=5).rev() {
        f1[i] = l_add(f1[i], f1[i - 1]);
        f2[i] = l_sub(f2[i], f2[i - 1]);
    }

    let mut a = [0i16; M + 1];
    a[0] = 4096;
    for i in 1..=5 {
        // eq (26): a_i = 0.5·(f'_1(i) + f'_2(i)),
        // a_{11−i} = 0.5·(f'_1(i) − f'_2(i)) — Q24 → Q12 is a 13-bit
        // rounding shift (12 bits plus the half).
        let sum = l_add(f1[i], f2[i]);
        a[i] = extract_l(l_shr_r(sum, 13));
        let diff = l_sub(f1[i], f2[i]);
        a[11 - i] = extract_l(l_shr_r(diff, 13));
    }
    a
}

/// Per-frame output of the fixed-point LP-parameter decode: the Q12
/// LP coefficient sets for both subframes plus the reconstructed
/// (post-stability) LSF vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LpParamsFx {
    /// Subframe-1 LP coefficients (interpolated LSPs), Q12.
    pub a_sub1: [i16; M + 1],
    /// Subframe-2 LP coefficients (current-frame LSPs), Q12.
    pub a_sub2: [i16; M + 1],
    /// The frame's reconstructed, stabilised LSF vector (Q13 radians).
    pub lsf_q13: [i16; M],
    /// The frame's LSP vector (Q15 cosines).
    pub lsp_q15: [i16; M],
}

/// Stateful §4.1.1 fixed-point LP-parameter decoder.
#[derive(Debug, Clone)]
pub struct LspDecoderFx {
    /// eq (20) MA residual history (Q13), most recent first.
    history: [[i16; M]; MA_NP],
    /// Previous frame's LSP vector (Q15) for §3.2.5 interpolation.
    prev_lsp_q15: [i16; M],
    /// Last good frame's stabilised LSF (Q13) for §4.4.1 erasure reuse.
    prev_lsf_q13: [i16; M],
    /// Last pushed residual, re-pushed on erased frames.
    last_residual: [i16; M],
    /// Fixed-point latitude shared with the encoder-side search.
    lat: FxLatitude,
}

impl Default for LspDecoderFx {
    fn default() -> Self {
        Self::new()
    }
}

impl LspDecoderFx {
    /// Fresh decoder in the clause-4.3 start-up state: MA history at
    /// the `l̂_i = iπ/11` vector, previous LSPs at `cos(iπ/11)`.
    ///
    /// The previous-LSP memory is the Table 9 `q_i` row on the exact
    /// Q15 grid `round(cos(iπ/11)·2^15)` (the printed `arccos` is a
    /// documented erratum for `cos` — §3.2.3 defines `q_i = cos(ω_i)`
    /// and the row above initialises `ω̂_i = iπ/11`), NOT the staged
    /// 64-segment cosine-table image of the startup LSF: the table
    /// interpolation lands up to 8 Q15 LSB off the true cosines, and
    /// that perturbation is measurable at frame 0 subframe 0 of every
    /// conformance vector through the §3.2.5 interpolated LP.
    #[must_use]
    pub fn new() -> Self {
        let history = startup_history_q13();
        let startup_lsf: [i16; M] = history[0];
        Self {
            history,
            prev_lsp_q15: STARTUP_LSP_Q15,
            prev_lsf_q13: startup_lsf,
            last_residual: startup_lsf,
            lat: FxLatitude::default(),
        }
    }

    /// Decodes one good frame's `(L0, L1, L2, L3)` into per-subframe
    /// Q12 LP coefficients, advancing the MA history and the
    /// interpolation state.
    pub fn decode_frame(&mut self, l0: usize, l1: usize, l2: usize, l3: usize) -> LpParamsFx {
        let (lsf_wide, residual) = reconstruct_lsf_q13(l0, l1, l2, l3, &self.history, &self.lat);
        // History advance with the rearranged residual.
        for k in (1..MA_NP).rev() {
            self.history[k] = self.history[k - 1];
        }
        self.history[0] = residual;
        self.last_residual = residual;

        let mut lsf: [i16; M] = std::array::from_fn(|i| lsf_wide[i] as i16);
        stability_q13(&mut lsf);
        self.finish_frame(lsf)
    }

    /// §4.4.1 erased frame: reuse the previous LSF vector and re-push
    /// the last residual so the MA history stays populated.
    pub fn decode_erased_frame(&mut self) -> LpParamsFx {
        for k in (1..MA_NP).rev() {
            self.history[k] = self.history[k - 1];
        }
        self.history[0] = self.last_residual;
        let lsf = self.prev_lsf_q13;
        self.finish_frame(lsf)
    }

    fn finish_frame(&mut self, lsf: [i16; M]) -> LpParamsFx {
        let lsp = lsf_to_lsp_q15(&lsf);
        let lsp_sub1 = interpolate_lsp_q15(&self.prev_lsp_q15, &lsp);
        let out = LpParamsFx {
            a_sub1: lsp_to_lp_q12(&lsp_sub1),
            a_sub2: lsp_to_lp_q12(&lsp),
            lsf_q13: lsf,
            lsp_q15: lsp,
        };
        self.prev_lsp_q15 = lsp;
        self.prev_lsf_q13 = lsf;
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cosine conversion tracks the real cosine within the
    /// table-interpolation tolerance across the whole stable LSF range.
    #[test]
    fn lsf_to_lsp_tracks_cosine() {
        let mut worst = 0.0f64;
        for w in (STABILITY_FLOOR_Q13..=STABILITY_CEIL_Q13).step_by(37) {
            let lsf = [w; M];
            let lsp = lsf_to_lsp_q15(&lsf);
            let expect = (f64::from(w) / 8192.0).cos() * 32768.0;
            let err = (f64::from(lsp[0]) - expect).abs();
            worst = worst.max(err);
        }
        assert!(worst <= 16.0, "worst cosine deviation {worst} Q15 LSBs");
    }

    /// The cosine-grid reading of Table 9 is `round(cos(iπ/11)·2^15)`
    /// exactly, antisymmetric and strictly decreasing as cosine on
    /// `(0, π)` demands; the corpus-pinned start-up memory is a valid
    /// (strictly decreasing, in-range) LSP vector.
    #[test]
    fn startup_lsp_grids_are_well_formed() {
        for (i, &q) in STARTUP_LSP_COS_GRID_Q15.iter().enumerate() {
            let w = f64::from((i + 1) as u32) * core::f64::consts::PI / 11.0;
            let expect = (w.cos() * 32768.0).round() as i16;
            assert_eq!(q, expect, "q_{} off the exact Q15 cosine grid", i + 1);
            assert_eq!(
                q,
                -STARTUP_LSP_COS_GRID_Q15[M - 1 - i],
                "antisymmetry at {i}"
            );
        }
        for i in 0..M {
            assert!(STARTUP_LSP_Q15[i].abs() < 32767);
            if i > 0 {
                assert!(
                    STARTUP_LSP_Q15[i] < STARTUP_LSP_Q15[i - 1],
                    "not strictly decreasing"
                );
                assert!(STARTUP_LSP_COS_GRID_Q15[i] < STARTUP_LSP_COS_GRID_Q15[i - 1]);
            }
        }
    }

    /// Monotonicity: increasing LSF maps to non-increasing LSP.
    #[test]
    fn lsf_to_lsp_monotone_decreasing() {
        let mut prev = i16::MAX;
        for w in (41..=25681).step_by(13) {
            let lsp = lsf_to_lsp_q15(&[w; M])[0];
            assert!(lsp <= prev, "cos not monotone at ω = {w}");
            prev = lsp;
        }
    }

    /// Stability procedure enforces its published invariants for
    /// arbitrary disordered inputs.
    #[test]
    fn stability_invariants() {
        let mut lsf: [i16; M] = [25000, 3, -50, 26000, 400, 380, 12000, 11990, 700, 700];
        stability_q13(&mut lsf);
        assert!(lsf[0] >= STABILITY_FLOOR_Q13);
        for i in 1..M {
            assert!(
                lsf[i] - lsf[i - 1] >= i32::from(STABILITY_GAP_Q13) as i16,
                "gap violated at {i}: {:?}",
                lsf
            );
        }
        assert!(lsf[M - 1] <= STABILITY_CEIL_Q13);
    }

    /// LSP → LP on a well-spread LSP vector reproduces the float
    /// conversion within a couple of Q12 LSBs (the float module is the
    /// spec-equation oracle; the fixed path must sit on its grid).
    #[test]
    fn lsp_to_lp_matches_float_reference() {
        // A stable, spread LSF vector (roughly iπ/11).
        let lsf: [i16; M] = std::array::from_fn(|i| (2340 * (i + 1)) as i16);
        let lsp = lsf_to_lsp_q15(&lsf);
        let a_fx = lsp_to_lp_q12(&lsp);

        let lsp_f: [f32; M] = std::array::from_fn(|i| f32::from(lsp[i]) / 32768.0);
        // Float module returns [a_1 … a_10] (0-based slots, no a_0).
        let a_float = crate::lsp_to_lp::lsp_to_lp(&lsp_f);
        assert_eq!(a_fx[0], 4096);
        for i in 1..=M {
            let expect = a_float[i - 1] * 4096.0;
            let got = f32::from(a_fx[i]);
            assert!(
                (got - expect).abs() <= 2.5,
                "a[{i}] fx {got} vs float {expect}"
            );
        }
    }

    /// Startup state: the first decoded frame from the clause-4.3
    /// state is deterministic and stable (no panics, ordered LSF).
    #[test]
    fn startup_frame_is_stable() {
        let mut dec = LspDecoderFx::new();
        let out = dec.decode_frame(0, 0, 0, 0);
        for i in 1..M {
            assert!(out.lsf_q13[i] > out.lsf_q13[i - 1]);
        }
        assert_eq!(out.a_sub1[0], 4096);
        assert_eq!(out.a_sub2[0], 4096);
    }

    /// Erased frames reuse the previous LSF and keep the decoder
    /// well-formed.
    #[test]
    fn erased_frame_repeats_lsf() {
        let mut dec = LspDecoderFx::new();
        let good = dec.decode_frame(1, 40, 7, 3);
        let erased = dec.decode_erased_frame();
        assert_eq!(erased.lsf_q13, good.lsf_q13);
        assert_eq!(erased.lsp_q15, good.lsp_q15);
        // Interpolation: prev == curr, so sub1 differs from sub2 only
        // by the truncating-halves LSB of `(x>>1) + (x>>1)`.
        for i in 0..=M {
            assert!((i32::from(erased.a_sub1[i]) - i32::from(erased.a_sub2[i])).abs() <= 2);
        }
    }
}
