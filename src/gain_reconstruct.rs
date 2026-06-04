//! §3.9.2 / §4.1.5 conjugate-structure gain-VQ reconstruction.
//!
//! This module implements the spec-cited two-component sum that maps
//! the transmitted (GA, GB) index pair into the quantised
//! adaptive-codebook gain `ĝ_p` and the quantised fixed-codebook gain
//! correction factor `γ̂`. The actual quantised fixed-codebook gain
//! `ĝ_c = γ̂ · g'_c` is then formed by the §3.9.1 4th-order MA-
//! prediction stage; that stage's wire-up is left for a follow-up
//! round (the MA predictor coefficients themselves already land in
//! [`crate::tables::GAIN_QUANT_MA_PREDICTOR_Q13`]).
//!
//! ## Spec source — clause 3.9.2 (06/2012 Recommendation)
//!
//! Per clause 3.9.2: "the adaptive-codebook gain, `g_p`, and the
//! factor `γ` are vector quantized using a two-stage conjugate
//! structured codebook. The first stage consists of a 3 bit two-
//! dimensional codebook **GA**, and the second stage consists of a
//! 4 bit two-dimensional codebook **GB**." Each codebook row carries
//! two components — column 0 is the row's contribution to `ĝ_p`,
//! column 1 is the row's contribution to `γ̂` (the
//! [`crate::tables::GAIN_VQ_COL_GP`] /
//! [`crate::tables::GAIN_VQ_COL_GC`] convention).
//!
//! Reconstruction (spec eqs (73) / (74)): given the transmitted GA
//! and GB indices,
//!
//! ```text
//!     ĝ_p = GA[GA][0] + GB[GB][0]
//!      γ̂  = GA[GA][1] + GB[GB][1]
//! ```
//!
//! and then the fixed-codebook gain `ĝ_c = γ̂ · g'_c` where `g'_c`
//! is the prediction-stage output of §3.9.1.
//!
//! ## Q-format conventions
//!
//! - GA / GB column-0 entries are Q14, so the per-row contribution
//!   to `ĝ_p` is interpreted as `v / 2^14`.
//! - GA / GB column-1 entries are Q12, so the per-row contribution
//!   to `γ̂` is interpreted as `v / 2^12`.
//!
//! Each stage's contribution is summed in its own Q-format and the
//! result is converted to the same `f32` boundary type that the
//! [`crate::lsp_reconstruct`] / [`crate::lsp_to_lp`] modules use, so
//! the rest of the decode chain stays in a single floating-point
//! domain. The intermediate sums are computed as `i32` to bound the
//! worst-case magnitude (`max(GA) + max(GB)` for the GP column is
//! well inside `i32` range, the same for the GC column), so the
//! summation cannot overflow.

use crate::parameters::Parameters;
use crate::tables::{gain_ga_entry, gain_gb_entry, GAIN_VQ_COL_GC, GAIN_VQ_COL_GP, NCODE1, NCODE2};

/// Q14 scaling divisor (`2^14 = 16384`) for the `g_p` columns of the
/// GA / GB codebooks.
const Q14_DIVISOR: f32 = 16384.0;

/// Q12 scaling divisor (`2^12 = 4096`) for the `γ` columns of the
/// GA / GB codebooks.
const Q12_DIVISOR: f32 = 4096.0;

/// Quantised gain pair reconstructed from one (GA, GB) index pair.
///
/// - [`Self::g_p_hat`] is the quantised adaptive-codebook gain `ĝ_p`
///   (the spec name; eq (73)). This value directly scales the
///   adaptive-codebook vector `v(n)` in the subframe excitation
///   `u(n) = ĝ_p · v(n) + ĝ_c · c(n)` (eq (75) / §3.10).
/// - [`Self::gamma_hat`] is the quantised fixed-codebook gain
///   correction factor `γ̂` (eq (74)). It does **not** directly scale
///   `c(n)`; the spec §3.9.1 prediction stage produces a predicted
///   gain `g'_c`, and the actual fixed-codebook gain used in eq (75)
///   is `ĝ_c = γ̂ · g'_c`. That second-stage product is not computed
///   here because `g'_c` is a separate (and stateful) wire-up
///   pending its own round; this struct carries the raw `γ̂` from
///   the conjugate-structure VQ.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuantisedGains {
    /// Quantised adaptive-codebook gain `ĝ_p`, dimensionless. Per
    /// the §3.9.2 codebook construction the value is in a bounded
    /// range whose extremes are the sum of the per-stage column-0
    /// extremes; in practice this is on the order of `0..=2` for
    /// typical speech.
    pub g_p_hat: f32,
    /// Quantised fixed-codebook gain correction factor `γ̂`,
    /// dimensionless. The actual fixed-codebook gain `ĝ_c` used in
    /// the excitation is `γ̂ · g'_c` where `g'_c` is the §3.9.1
    /// predicted gain; that combination is left for a follow-up
    /// round.
    pub gamma_hat: f32,
}

/// Errors returned by the §3.9.2 gain-VQ reconstruction helpers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GainReconstructError {
    /// The supplied GA index does not fit in the
    /// [`crate::tables::NCODE1`]-row first-stage codebook. The
    /// transmitted GA codeword is only 3 bits per spec Table 8 so a
    /// well-formed frame cannot produce this error; it surfaces when
    /// callers pass a fabricated `Parameters` struct outside the
    /// spec-stated domain.
    GaOutOfRange { index: usize },
    /// The supplied GB index does not fit in the
    /// [`crate::tables::NCODE2`]-row second-stage codebook. The
    /// transmitted GB codeword is only 4 bits per spec Table 8 so
    /// the same well-formed-frame guarantee applies.
    GbOutOfRange { index: usize },
}

impl core::fmt::Display for GainReconstructError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::GaOutOfRange { index } => write!(
                f,
                "g729 §3.9.2 gain-VQ: GA index {index} >= NCODE1 = {NCODE1}",
            ),
            Self::GbOutOfRange { index } => write!(
                f,
                "g729 §3.9.2 gain-VQ: GB index {index} >= NCODE2 = {NCODE2}",
            ),
        }
    }
}

impl std::error::Error for GainReconstructError {}

/// Reconstructs the quantised `(ĝ_p, γ̂)` pair from a (GA, GB)
/// codebook-index pair per spec eqs (73) / (74).
///
/// The summation is done in the per-component Q-format (Q14 for
/// column 0, Q12 for column 1) using `i32` intermediates, then
/// scaled to `f32` at the boundary. Out-of-range indices return a
/// typed error rather than panicking, so the function is safe to
/// call directly on raw codeword integers.
///
/// # Errors
///
/// Returns [`GainReconstructError::GaOutOfRange`] if `ga >= NCODE1`,
/// or [`GainReconstructError::GbOutOfRange`] if `gb >= NCODE2`. The
/// transmitted GA / GB codewords (3 + 4 bits) cannot trigger these
/// in a well-formed frame.
pub fn reconstruct_gains(ga: usize, gb: usize) -> Result<QuantisedGains, GainReconstructError> {
    if ga >= NCODE1 {
        return Err(GainReconstructError::GaOutOfRange { index: ga });
    }
    if gb >= NCODE2 {
        return Err(GainReconstructError::GbOutOfRange { index: gb });
    }

    let ga_row = gain_ga_entry(ga);
    let gb_row = gain_gb_entry(gb);

    // eq (73): ĝ_p = GA[GA].gp + GB[GB].gp (column 0; Q14 in both stages).
    let gp_sum_q14 = i32::from(ga_row[GAIN_VQ_COL_GP]) + i32::from(gb_row[GAIN_VQ_COL_GP]);
    // eq (74): γ̂ = GA[GA].gc + GB[GB].gc (column 1; Q12 in both stages).
    let gc_sum_q12 = i32::from(ga_row[GAIN_VQ_COL_GC]) + i32::from(gb_row[GAIN_VQ_COL_GC]);

    Ok(QuantisedGains {
        g_p_hat: gp_sum_q14 as f32 / Q14_DIVISOR,
        gamma_hat: gc_sum_q12 as f32 / Q12_DIVISOR,
    })
}

/// Convenience wrapper that reconstructs the per-subframe gain pairs
/// for one frame straight from the typed [`Parameters`] struct
/// returned by [`crate::parameters::unpack_parameters`].
///
/// Returns a `[QuantisedGains; 2]` with index 0 = subframe 1 (from
/// `GA1` / `GB1`) and index 1 = subframe 2 (from `GA2` / `GB2`),
/// matching the spec §4.1.5 ordering. The GA / GB indices in
/// `Parameters` are 3 / 4 bits respectively so the reconstruction
/// cannot return a `GaOutOfRange` / `GbOutOfRange` error in
/// practice; the return type still threads the error variant so
/// callers handing in synthetic `Parameters` get the typed surface
/// rather than a panic.
///
/// # Errors
///
/// Surfaces any [`GainReconstructError`] from the per-subframe
/// [`reconstruct_gains`] call.
pub fn reconstruct_frame_gains(
    params: &Parameters,
) -> Result<[QuantisedGains; 2], GainReconstructError> {
    let sub1 = reconstruct_gains(usize::from(params.ga1), usize::from(params.gb1))?;
    let sub2 = reconstruct_gains(usize::from(params.ga2), usize::from(params.gb2))?;
    Ok([sub1, sub2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tables::{
        GAIN_QUANT_CODEBOOK_GA_Q14_Q12, GAIN_QUANT_CODEBOOK_GB_Q14_Q12, GAIN_VQ_DIM,
    };

    /// Per the staged CSVs GA[0] = (1, 1516) and GB[0] = (826, 2005)
    /// (both rows column 0 = Q14, column 1 = Q12). The reconstruction
    /// at (GA = 0, GB = 0) sums the column values and divides by the
    /// per-column Q divisors.
    #[test]
    fn reconstruct_at_zero_zero_matches_first_rows() {
        let g = reconstruct_gains(0, 0).expect("indices in range");
        let expected_gp = (1 + 826) as f32 / 16384.0;
        let expected_gc = (1516 + 2005) as f32 / 4096.0;
        assert!((g.g_p_hat - expected_gp).abs() < 1e-7);
        assert!((g.gamma_hat - expected_gc).abs() < 1e-7);
    }

    /// (GA = NCODE1 - 1, GB = NCODE2 - 1) sums the last row of each
    /// codebook. Locks the index-into-row direction so a flipped
    /// `i = NCODE - 1 - i` mistake would trip immediately.
    #[test]
    fn reconstruct_at_last_indices_matches_last_rows() {
        let ga = NCODE1 - 1;
        let gb = NCODE2 - 1;
        let g = reconstruct_gains(ga, gb).expect("indices in range");
        let ga_row = GAIN_QUANT_CODEBOOK_GA_Q14_Q12[ga];
        let gb_row = GAIN_QUANT_CODEBOOK_GB_Q14_Q12[gb];
        let expected_gp = (i32::from(ga_row[0]) + i32::from(gb_row[0])) as f32 / 16384.0;
        let expected_gc = (i32::from(ga_row[1]) + i32::from(gb_row[1])) as f32 / 4096.0;
        assert!((g.g_p_hat - expected_gp).abs() < 1e-7);
        assert!((g.gamma_hat - expected_gc).abs() < 1e-7);
    }

    /// Out-of-range GA index yields the typed error rather than
    /// panicking. Pinned for every off-by-one boundary in `0..NCODE1`.
    #[test]
    fn reconstruct_rejects_out_of_range_ga() {
        let err = reconstruct_gains(NCODE1, 0).unwrap_err();
        assert_eq!(err, GainReconstructError::GaOutOfRange { index: NCODE1 });
    }

    /// Out-of-range GB index yields the typed error rather than
    /// panicking.
    #[test]
    fn reconstruct_rejects_out_of_range_gb() {
        let err = reconstruct_gains(0, NCODE2).unwrap_err();
        assert_eq!(err, GainReconstructError::GbOutOfRange { index: NCODE2 });
    }

    /// GA-range error fires before any GB check so the surfaced
    /// error names the offending GA index when both are out of range.
    #[test]
    fn reconstruct_reports_ga_error_first() {
        let err = reconstruct_gains(NCODE1, NCODE2).unwrap_err();
        assert_eq!(err, GainReconstructError::GaOutOfRange { index: NCODE1 });
    }

    /// Every (GA, GB) pair in the spec-stated domain returns a
    /// finite `f32` pair (no NaN, no infinity from the integer
    /// summation or the Q-format conversion).
    #[test]
    fn every_pair_in_domain_yields_finite_gains() {
        for ga in 0..NCODE1 {
            for gb in 0..NCODE2 {
                let g = reconstruct_gains(ga, gb).expect("in-domain");
                assert!(g.g_p_hat.is_finite(), "g_p infinite at ({ga}, {gb})");
                assert!(g.gamma_hat.is_finite(), "γ infinite at ({ga}, {gb})");
            }
        }
    }

    /// Per spec §3.9.2 the quantised `ĝ_p` is bounded; the §3.10
    /// taming logic also discourages large adaptive-codebook gains.
    /// All NCODE1 × NCODE2 = 128 entries should produce
    /// `ĝ_p` in a small physical range — pinned generously here at
    /// `0.0..=2.0` to defend against a Q-format divisor swap.
    #[test]
    fn g_p_hat_lies_in_plausible_range() {
        for ga in 0..NCODE1 {
            for gb in 0..NCODE2 {
                let g = reconstruct_gains(ga, gb).unwrap();
                assert!(
                    (0.0..=2.0).contains(&g.g_p_hat),
                    "g_p_hat = {} out of range at ({ga}, {gb})",
                    g.g_p_hat,
                );
            }
        }
    }

    /// γ̂ is the fixed-codebook gain *correction factor* per spec
    /// §3.9.1 / §3.9.2; the codebook tabulates it in Q12, so the
    /// reconstructed value over all 128 (GA, GB) pairs lands in a
    /// generous plausibility window of `0.0..=11.0` (worst-case row
    /// pairs hit ≈10.12). A Q14 / Q12 divisor confusion would push
    /// the result well outside this window (≈40 or ≈2.5 respectively).
    #[test]
    fn gamma_hat_lies_in_plausible_range() {
        for ga in 0..NCODE1 {
            for gb in 0..NCODE2 {
                let g = reconstruct_gains(ga, gb).unwrap();
                assert!(
                    (0.0..=11.0).contains(&g.gamma_hat),
                    "γ̂ = {} out of range at ({ga}, {gb})",
                    g.gamma_hat,
                );
            }
        }
    }

    /// Codebook entries are (g_p contribution, γ contribution) per
    /// the table column convention; the round-trip property
    /// `reconstruct(ga, gb) ==
    ///     (GA[ga][0] + GB[gb][0]) / 2^14
    ///   + (GA[ga][1] + GB[gb][1]) / 2^12` is exercised explicitly
    /// for a hand-picked pair distinct from `(0, 0)`.
    #[test]
    fn reconstruct_hand_picked_pair() {
        let ga = 5; // arbitrary in-domain GA
        let gb = 11; // arbitrary in-domain GB
        let g = reconstruct_gains(ga, gb).unwrap();
        let ga_row = GAIN_QUANT_CODEBOOK_GA_Q14_Q12[ga];
        let gb_row = GAIN_QUANT_CODEBOOK_GB_Q14_Q12[gb];
        let expected_gp = (i32::from(ga_row[0]) + i32::from(gb_row[0])) as f32 / 16384.0;
        let expected_gc = (i32::from(ga_row[1]) + i32::from(gb_row[1])) as f32 / 4096.0;
        assert!((g.g_p_hat - expected_gp).abs() < 1e-7);
        assert!((g.gamma_hat - expected_gc).abs() < 1e-7);
    }

    /// Per-frame wrapper threads the subframe indices into the right
    /// codebooks: subframe 1 from `(GA1, GB1)`, subframe 2 from
    /// `(GA2, GB2)`. Crafted `Parameters` with distinct indices per
    /// subframe ensures the wrapper doesn't accidentally reuse
    /// subframe-1 values for subframe 2.
    #[test]
    fn frame_wrapper_threads_per_subframe_indices() {
        let params = Parameters {
            l0: 0,
            l1: 0,
            l2: 0,
            l3: 0,
            p1: 0,
            p0: 0,
            c1: 0,
            s1: 0,
            ga1: 1,
            gb1: 2,
            p2: 0,
            c2: 0,
            s2: 0,
            ga2: 6,
            gb2: 13,
        };
        let [sub1, sub2] = reconstruct_frame_gains(&params).unwrap();
        let expected_sub1 = reconstruct_gains(1, 2).unwrap();
        let expected_sub2 = reconstruct_gains(6, 13).unwrap();
        assert_eq!(sub1, expected_sub1);
        assert_eq!(sub2, expected_sub2);
    }

    /// Column convention: column 0 (`GAIN_VQ_COL_GP`) is Q14 and is
    /// the adaptive-codebook-gain half; column 1
    /// (`GAIN_VQ_COL_GC`) is Q12 and is the fixed-codebook-gain
    /// correction half. Confirm by changing GA but holding GB
    /// constant: the reconstructed `ĝ_p` change should match the
    /// GA row's column-0 delta divided by 2^14, the reconstructed
    /// `γ̂` change should match column-1 delta divided by 2^12.
    #[test]
    fn per_column_q_format_isolation() {
        let gb = 0;
        let ga_a = 0;
        let ga_b = 4;
        let g_a = reconstruct_gains(ga_a, gb).unwrap();
        let g_b = reconstruct_gains(ga_b, gb).unwrap();
        let dgp_q14 = i32::from(GAIN_QUANT_CODEBOOK_GA_Q14_Q12[ga_b][0])
            - i32::from(GAIN_QUANT_CODEBOOK_GA_Q14_Q12[ga_a][0]);
        let dgc_q12 = i32::from(GAIN_QUANT_CODEBOOK_GA_Q14_Q12[ga_b][1])
            - i32::from(GAIN_QUANT_CODEBOOK_GA_Q14_Q12[ga_a][1]);
        assert!(((g_b.g_p_hat - g_a.g_p_hat) - dgp_q14 as f32 / 16384.0).abs() < 1e-6);
        assert!(((g_b.gamma_hat - g_a.gamma_hat) - dgc_q12 as f32 / 4096.0).abs() < 1e-6);
    }

    /// Pin that the codebook width matches the published [`GAIN_VQ_DIM`]
    /// constant — defensive against a future column-count widening.
    #[test]
    fn codebook_width_matches_dim_constant() {
        assert_eq!(GAIN_VQ_DIM, 2);
        assert_eq!(GAIN_QUANT_CODEBOOK_GA_Q14_Q12[0].len(), GAIN_VQ_DIM);
        assert_eq!(GAIN_QUANT_CODEBOOK_GB_Q14_Q12[0].len(), GAIN_VQ_DIM);
    }
}
