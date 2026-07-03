//! §3.2.3 / §3.2.4 **fixed-point LSP↔LSF conversion** — the eq (18)
//! `ω_i = arccos(q_i)` evaluated the way the 16-bit fixed-point
//! implementation evaluates it: a 64-segment table lookup with a
//! per-segment linear slope refinement, instead of a transcendental
//! call.
//!
//! ## Why this exists
//!
//! Clause 2.5 states the Recommendation's decimal values are "rounded
//! versions of the values used in the 16-bit fixed-point ANSI C
//! implementation" — the conformance bitstreams were produced by an
//! encoder whose LSFs live on a fixed-point grid. The §3.2.4 quantiser
//! makes razor-thin nearest-neighbour decisions on those LSFs, so the
//! *precision profile* of the arccos matters: evaluating eq (18) in
//! floating point yields LSFs that differ from the fixed-point grid by
//! a few Q13 steps, which measurably flips codebook decisions against
//! the reference corpus. Routing the encoder's eq (18) through this
//! table path (black-box-validated against the conformance corpus)
//! moves every downstream decision onto the reference grid.
//!
//! ## Table layout (compiled from the staged CSVs)
//!
//! * [`crate::tables::LSF_LSP_COS_TABLE_Q15`] — 65 samples of
//!   `cos(ω)` in Q15 on the uniform grid `ω = i·π/64`, `i = 0 … 64`
//!   (`+32767 … −32768`, strictly decreasing).
//! * [`crate::tables::LSF_LSP_ACOS_SLOPE_Q12`] — 64 per-segment
//!   arccos slopes in Q12: within segment `i` the normalised-frequency
//!   offset from the segment start is `((x − cos_table[i]) ·
//!   slope[i]) >> 12` in the Q15 normalised-frequency domain (`π` ↔
//!   `16384`, i.e. 256 units per segment).
//!
//! The returned LSF is converted to the crate-wide **Q13 radian**
//! domain (`π` ↔ `25736`, the domain of the §3.2.4 LSP codebooks):
//! `ω_q13 = (f_q15 · 25736) >> 14`.

use crate::tables::{LSF_LSP_ACOS_SLOPE_Q12, LSF_LSP_COS_TABLE_Q15};

/// π in the Q13 radian LSF domain (`⌊π·8192⌋`, the LSP-codebook unit).
pub const PI_Q13: i32 = 25_736;

/// π in the Q15 normalised-frequency domain of the lookup grid.
pub const PI_FREQ_Q15: i32 = 16_384;

/// Number of uniform arccos segments (`ω` step `π/64`).
pub const N_SEGMENTS: usize = 64;

/// Normalised-frequency units per segment (`16384 / 64`).
pub const SEGMENT_STEP: i32 = 256;

/// Fixed-point eq (18): arccos of a Q15 cosine-domain LSP `x_q15 ∈
/// [−32768, 32767]`, returned in the Q15 normalised-frequency domain
/// (`0 … 16384` for `ω ∈ [0, π]`).
///
/// Segment `i` is located by walking the strictly-decreasing cos table
/// (`cos_table[i] ≥ x > cos_table[i+1]`); the in-segment offset is the
/// Q12 slope refinement `((x − cos_table[i]) · slope[i]) >> 12`.
#[must_use]
pub fn acos_q15_to_freq_q15(x_q15: i32) -> i32 {
    let x = x_q15.clamp(-32_768, 32_767);
    let mut i = 0usize;
    while i < N_SEGMENTS - 1 && i32::from(LSF_LSP_COS_TABLE_Q15[i + 1]) >= x {
        i += 1;
    }
    let dx = x - i32::from(LSF_LSP_COS_TABLE_Q15[i]);
    let frac = ((i64::from(dx) * i64::from(LSF_LSP_ACOS_SLOPE_Q12[i])) >> 12) as i32;
    (i as i32) * SEGMENT_STEP + frac
}

/// Fixed-point eq (18) into the **Q13 radian** LSF domain (`π` ↔
/// [`PI_Q13`]): the [`acos_q15_to_freq_q15`] lookup followed by the
/// domain conversion `(f_q15 · 25736) >> 14`.
#[must_use]
pub fn acos_q15_to_lsf_q13(x_q15: i32) -> i32 {
    let f = acos_q15_to_freq_q15(x_q15);
    ((i64::from(f) * i64::from(PI_Q13)) >> 14) as i32
}

/// Convenience for the float-surfaced encoder chain: quantises a
/// cosine-domain LSP `q ∈ [−1, 1]` to Q15 and evaluates the
/// fixed-point eq (18), returning the LSF in Q13 radians.
#[must_use]
pub fn lsp_to_lsf_q13(q: f32) -> i32 {
    let x_q15 = (f64::from(q) * 32_768.0).round() as i32;
    acos_q15_to_lsf_q13(x_q15)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The compiled cos table is strictly decreasing over the full grid
    /// (the precondition of the segment walk).
    #[test]
    fn cos_table_strictly_decreasing() {
        for i in 1..=N_SEGMENTS {
            assert!(
                LSF_LSP_COS_TABLE_Q15[i] < LSF_LSP_COS_TABLE_Q15[i - 1],
                "cos table not decreasing at {i}"
            );
        }
    }

    /// The cos table matches `cos(i·π/64)` to within one Q15 step, and
    /// the endpoints are the exact Q15 saturation values.
    #[test]
    fn cos_table_matches_cosine() {
        assert_eq!(LSF_LSP_COS_TABLE_Q15[0], 32_767);
        assert_eq!(LSF_LSP_COS_TABLE_Q15[N_SEGMENTS], -32_768);
        for (i, &entry) in LSF_LSP_COS_TABLE_Q15.iter().enumerate() {
            let ideal = (f64::from(i as u32) * std::f64::consts::PI / 64.0).cos() * 32_768.0;
            let got = f64::from(entry);
            assert!(
                (got - ideal).abs() <= 1.5,
                "cos_table[{i}] = {got} vs ideal {ideal}"
            );
        }
    }

    /// At every segment boundary the lookup returns the exact grid
    /// frequency `i·256` (the slope term vanishes at `dx = 0`).
    #[test]
    fn exact_at_segment_boundaries() {
        for (i, &entry) in LSF_LSP_COS_TABLE_Q15.iter().take(N_SEGMENTS).enumerate() {
            let x = i32::from(entry);
            assert_eq!(
                acos_q15_to_freq_q15(x),
                (i as i32) * SEGMENT_STEP,
                "boundary {i}"
            );
        }
    }

    /// Whole-domain accuracy sweep. The 64-segment linear
    /// interpolation has an intrinsic error profile (this *is* the
    /// grid the fixed-point reference's LSFs live on): measured worst
    /// deviation from the ideal `round(ω·8192)` is 22 Q13 steps inside
    /// the eq (22) LSF operating band `[0.04π, 0.96π]` (where arccos
    /// curvature per segment is modest) and 115 steps at the extreme
    /// domain ends (where `d(arccos)/dx` blows up and Q15 can no
    /// longer resolve ω). The output must also be monotone
    /// non-decreasing in ω.
    #[test]
    fn round_trip_accuracy_and_monotonicity() {
        let lo = 0.04 * std::f64::consts::PI;
        let hi = 0.96 * std::f64::consts::PI;
        let mut prev = -1i32;
        for k in 0..=4096 {
            let omega = f64::from(k) * std::f64::consts::PI / 4096.0;
            let x = (omega.cos() * 32_768.0).round() as i32;
            let lsf = acos_q15_to_lsf_q13(x.clamp(-32_768, 32_767));
            let ideal = (omega * 8192.0).round() as i32;
            let bound = if (lo..=hi).contains(&omega) { 24 } else { 120 };
            assert!(
                (lsf - ideal).abs() <= bound,
                "omega {omega}: lsf {lsf} vs ideal {ideal} (bound {bound})"
            );
            assert!(lsf >= prev, "non-monotone at omega {omega}");
            prev = lsf;
        }
    }

    /// Domain endpoints: `x = +1` maps to ω = 0 and `x = −1` maps to
    /// ω = π in both output domains.
    #[test]
    fn endpoints() {
        assert_eq!(acos_q15_to_freq_q15(32_767), 0);
        assert_eq!(acos_q15_to_freq_q15(-32_768), PI_FREQ_Q15);
        assert_eq!(acos_q15_to_lsf_q13(-32_768), PI_Q13);
        assert_eq!(lsp_to_lsf_q13(1.0), 0);
        assert_eq!(lsp_to_lsf_q13(-1.0), PI_Q13);
    }

    /// Out-of-domain inputs are clamped, not wrapped.
    #[test]
    fn clamps_out_of_domain() {
        assert_eq!(acos_q15_to_freq_q15(40_000), acos_q15_to_freq_q15(32_767));
        assert_eq!(acos_q15_to_freq_q15(-40_000), acos_q15_to_freq_q15(-32_768));
    }
}
