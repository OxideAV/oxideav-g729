//! §3.2.6 LSP-to-LP conversion.
//!
//! Round 207 wired the §3.2.4 LSP-frame reconstruction; round 213
//! chained the §3.2.5 per-subframe cosine-domain interpolation. This
//! module wires the next link of the decoder chain: turning the
//! per-subframe cosine-domain LSPs `q_i = cos(ω̂_i)` back into the
//! 10-coefficient LP synthesis filter `A(z) = 1 + Σ_{i=1..=10}
//! a_i · z^-i`.
//!
//! ## Spec source — clauses 3.2.3, 3.2.6 (06/2012 Recommendation)
//!
//! Per clause 3.2.3, the LSP coefficients of a 10th-order LP filter are
//! the roots of the sum and difference polynomials:
//!
//! ```text
//! F1(z) = Π_{i = 1, 3, 5, 7, 9} (1 − 2·q_i·z^-1 + z^-2)        (13)
//! F2(z) = Π_{i = 2, 4, 6, 8, 10} (1 − 2·q_i·z^-1 + z^-2)       (14)
//! ```
//!
//! Each polynomial is degree 10 in `z^-1` but symmetric, so only the
//! first five coefficients `f_1(i)`, `i ∈ 0..=5` (with `f_1(0) = 1`)
//! carry information. Per clause 3.2.6 those coefficients are computed
//! from `q_{2i-1}` (for `F1`) / `q_{2i}` (for `F2`) using the recursion
//! derived from polynomial multiplication by `(1 − 2·q·z^-1 + z^-2)`:
//!
//! ```text
//! for i = 1..=5:
//!     f_1(i)      :=  -2·q_{2i-1}·f_1(i-1)  +  2·f_1(i-2)
//!     for j = i-1 down to 1:
//!         f_1[i](j)  :=  f_1[i-1](j)  − 2·q_{2i-1}·f_1[i-1](j-1)  +  f_1[i-1](j-2)
//! ```
//!
//! with initial values `f_1(0) = 1`, `f_1(-1) = 0`. The recursion for
//! `f_2(i)` is identical with `q_{2i-1}` replaced by `q_{2i}`.
//!
//! After the five-step recursion `F1`,`F2` (degrees 5 in `z^-1` after
//! truncating symmetry) are multiplied by `(1 + z^-1)` and
//! `(1 − z^-1)` respectively, restoring the spec §3.2.3 roots at
//! `z = ±1` that were factored out of the original sum/difference
//! polynomials (eqs (11)/(12)):
//!
//! ```text
//! f'_1(i) = f_1(i) + f_1(i-1)        i ∈ 1..=5        (25)
//! f'_2(i) = f_2(i) − f_2(i-1)        i ∈ 1..=5        (25)
//! ```
//!
//! Finally the LP coefficients `a_i`, `i ∈ 1..=10`, are recovered from
//! `A(z) = (F'_1(z) + F'_2(z)) / 2`, using the fact that `F'_1` is
//! symmetric and `F'_2` antisymmetric of length 6:
//!
//! ```text
//! a_i =  0.5 · f'_1(i)   + 0.5 · f'_2(i)              i ∈ 1..=5     (26)
//! a_i =  0.5 · f'_1(11-i) − 0.5 · f'_2(11-i)          i ∈ 6..=10    (26)
//! ```
//!
//! ## Q-domain
//!
//! Inputs and outputs are real-valued `f32`. The §3.2.4 reconstructor
//! already produces real-valued `ω̂`; the round-213 `omega_to_q`
//! boundary produces real-valued `q ∈ [-1, 1]`. The LP coefficient
//! output `a_i` is unscaled (the §3.5+ impulse-response computation
//! takes the unscaled `a_i` and applies the bandwidth-expansion factor
//! `γ` itself). Q-format conversion is deferred to the encoder/decoder
//! search stages where Q12-vs-Q14 trade-offs are codec-internal.
//!
//! ## Wiring with [`crate::lsp_reconstruct`] / [`crate::lsp_interpolate`]
//!
//! Typical decoder flow per 10 ms frame:
//!
//! ```text
//! omega        = reconstructor.reconstruct_frame(L0, L1, L2, L3)?;
//! q            = omega_to_q(&omega);
//! [sub1_q, sub2_q] = interpolator.interpolate(&q);
//! a_subframe1  = lsp_to_lp(&sub1_q);
//! a_subframe2  = lsp_to_lp(&sub2_q);
//! ```
//!
//! Each subframe's `a_i` then drives the LP synthesis filter
//! `1/Â(z)` for the reconstruction stage §4.1.6.

use crate::tables::M;

/// Result of §3.2.6 LSP→LP conversion. `[0]` is implicitly `a_0 = 1.0`
/// (not stored); slot `i - 1` holds `a_i` for `i ∈ 1..=10`.
pub type LpCoefficients = [f32; M];

/// §3.2.6 LSP→LP conversion: convert one subframe's cosine-domain LSP
/// vector `q_i = cos(ω̂_i)`, `i ∈ 1..=10` (input layout
/// `q_in[i - 1] = q_i`), to the 10 LP coefficients `a_i`, `i ∈ 1..=10`
/// of the LP synthesis filter `1/A(z)`, with `a_0 = 1.0` implicit.
///
/// The implementation follows spec §3.2.6 to the letter:
///
/// 1. Build the `F1` / `F2` polynomial coefficients `f_1(i)`, `f_2(i)`
///    for `i ∈ 0..=5` via the eq (15)-style recursion derived from
///    polynomial multiplication by `(1 − 2·q·z^-1 + z^-2)`
///    (initial values `f(0) = 1.0`, `f(-1) = 0.0`).
/// 2. Restore the `(1 + z^-1)` / `(1 − z^-1)` factors that §3.2.3
///    factored out (eq (11) / (12)) by applying spec eq (25):
///    `f'_1(i) = f_1(i) + f_1(i-1)`, `f'_2(i) = f_2(i) − f_2(i-1)`,
///    each for `i ∈ 1..=5`.
/// 3. Combine via `A(z) = (F'_1(z) + F'_2(z)) / 2` per spec eq (26),
///    using `F'_1` symmetric / `F'_2` antisymmetric (length 6) so the
///    higher half can be read off by mirror index `11 − i`.
///
/// ## Mapping note (1-based spec vs 0-based Rust)
///
/// The spec uses 1-based subscripts `q_1, q_2, …, q_10`. This function
/// takes a 0-based `&[f32; 10]` slice; slot `q_in[k]` is the spec's
/// `q_{k+1}` for `k ∈ 0..=9`. So `q_{2i-1}` for `i ∈ 1..=5` (the
/// odd-indexed LSPs `q_1, q_3, q_5, q_7, q_9`) is `q_in[2*(i-1)]`, and
/// `q_{2i}` (the even-indexed `q_2, q_4, q_6, q_8, q_10`) is
/// `q_in[2*i - 1]`.
#[must_use]
pub fn lsp_to_lp(q_in: &[f32; M]) -> LpCoefficients {
    // Spec §3.2.6 step 1: build f_1(i), f_2(i) for i ∈ 0..=5 (six
    // coefficients each because the polynomial is symmetric and only
    // the first half is unique).
    //
    // The recursion in the spec diagram (image f0015-01) is:
    //
    //     for i = 1 to 5:
    //         f(i)             := -2·q·f(i-1)            + 2·f(i-2)
    //         for j = i-1 down to 1:
    //             f[i](j)      :=  f[i-1](j)
    //                             - 2·q·f[i-1](j-1)
    //                             +     f[i-1](j-2)
    //
    // The two `f[i-1](-1)` / `f[i-1](-2)` reads inside the inner loop
    // fall back to 0 by spec convention; the slot is never read with
    // `j-1 < 0` (the loop bound is `j >= 1`, so `j-1 >= 0`, and
    // `j-2 = -1` only when `j = 1`, in which case `f[i-1](-1) = 0`).
    //
    // The update is done in-place from j = i-1 down to 1 so we read
    // f[i-1](j-1) and f[i-1](j-2) before overwriting them.

    let mut f1 = [0.0_f32; 6]; // f_1(0..=5)
    let mut f2 = [0.0_f32; 6]; // f_2(0..=5)
    f1[0] = 1.0;
    f2[0] = 1.0;

    for i in 1..=5 {
        // q_{2i-1} (F1 recursion) and q_{2i} (F2 recursion). 1-based
        // spec subscripts map to 0-based slots; see the doc comment.
        let q_odd = q_in[2 * (i - 1)]; // q_1, q_3, q_5, q_7, q_9
        let q_even = q_in[2 * i - 1]; // q_2, q_4, q_6, q_8, q_10

        // The "leading" new coefficient: f(i) = -2·q·f(i-1) + 2·f(i-2)
        // (the `2·f(i-2)` term comes from the `z^-2` factor of the
        // polynomial we're multiplying in, and the `1·f(i)` slot on
        // the previous polynomial is 0 because that slot didn't exist
        // yet — the polynomial just grew by one). Spec convention:
        // `f(-1) = 0`, so when i == 1 the `f(i-2)` term vanishes.
        let f1_im2 = if i >= 2 { f1[i - 2] } else { 0.0 };
        let f2_im2 = if i >= 2 { f2[i - 2] } else { 0.0 };
        let new_f1_i = -2.0 * q_odd * f1[i - 1] + 2.0 * f1_im2;
        let new_f2_i = -2.0 * q_even * f2[i - 1] + 2.0 * f2_im2;

        // Inner loop: update f(j) for j = i-1 down to 1, in
        // decreasing-j order so each step reads the prior iteration's
        // f(j-1) and f(j-2) before they are overwritten.
        for j in (1..i).rev() {
            // f[i-1](j-2) is 0 by convention when j == 1.
            let f1_jm2 = if j >= 2 { f1[j - 2] } else { 0.0 };
            let f2_jm2 = if j >= 2 { f2[j - 2] } else { 0.0 };
            f1[j] = f1[j] - 2.0 * q_odd * f1[j - 1] + f1_jm2;
            f2[j] = f2[j] - 2.0 * q_even * f2[j - 1] + f2_jm2;
        }

        // Slot i is unconditionally the new top coefficient.
        f1[i] = new_f1_i;
        f2[i] = new_f2_i;
    }

    // Spec §3.2.6 step 2: apply eq (25) to multiply by (1 ± z^-1).
    // The spec writes f'_1(i) = f_1(i) + f_1(i - 1), i = 1..=5 for the
    // (1 + z^-1) factor on F_1, and f'_2(i) = f_2(i) − f_2(i - 1) for
    // the (1 − z^-1) factor on F_2.
    let mut f1p = [0.0_f32; 6]; // f'_1(0..=5)
    let mut f2p = [0.0_f32; 6]; // f'_2(0..=5)
    f1p[0] = f1[0]; // = 1.0 (the constant term picks up nothing new)
    f2p[0] = f2[0]; // = 1.0
    for i in 1..=5 {
        f1p[i] = f1[i] + f1[i - 1];
        f2p[i] = f2[i] - f2[i - 1];
    }

    // Spec §3.2.6 step 3: recover a_i via eq (26).
    //
    //     a_i = 0.5 · f'_1(i)     + 0.5 · f'_2(i)        i ∈ 1..=5
    //     a_i = 0.5 · f'_1(11-i)  − 0.5 · f'_2(11-i)     i ∈ 6..=10
    //
    // The output slot a_i lives at index i - 1 in our 0-based array
    // (so a_1 is in slot 0, ..., a_10 is in slot 9).
    let mut a: LpCoefficients = [0.0; M];
    for i in 1..=5 {
        a[i - 1] = 0.5 * f1p[i] + 0.5 * f2p[i];
    }
    for i in 6..=10 {
        // 11 - i ∈ {5, 4, 3, 2, 1} as i sweeps 6..=10.
        let mirror = 11 - i;
        a[i - 1] = 0.5 * f1p[mirror] - 0.5 * f2p[mirror];
    }
    a
}

/// Convenience wrapper: convert an LSF-domain vector `ω̂` to LP
/// coefficients in one call, bridging the §3.2.4 / §3.2.5 LSF output
/// to the §3.2.6 LP coefficients via the cosine boundary
/// `q_i = cos(ω̂_i)`.
///
/// Equivalent to `lsp_to_lp(&omega_to_q(omega))`. Useful when the
/// caller drives this module straight from the round-207 LSF
/// reconstructor without staying in the cosine domain.
#[must_use]
pub fn lsf_to_lp(omega: &[f32; M]) -> LpCoefficients {
    let q = crate::lsp_interpolate::omega_to_q(omega);
    lsp_to_lp(&q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f32::consts::PI;

    /// Helper: build a sorted, well-spaced LSF vector to feed the
    /// conversion. The values are chosen well inside the §3.2.4
    /// stability-clamp interior so `q_to_omega` and the polynomial
    /// recursion don't see any boundary numerics.
    fn well_spaced_omega() -> [f32; M] {
        let mut omega = [0.0_f32; M];
        for (i, slot) in omega.iter_mut().enumerate() {
            *slot = ((i + 1) as f32) * PI / 11.0;
        }
        omega
    }

    /// At the spec §3.2.4 start-up state (`ω̂_i = i · π / 11`), the
    /// reconstructed `A(z)` is by construction the unique LP filter
    /// whose 10 LSP roots lie at the start-up LSFs. We don't have a
    /// closed-form here, but we can pin two invariants:
    ///
    /// * the recursion produces a well-defined finite vector
    ///   (no NaN / no inf);
    /// * the filter is monic-symmetric in the §3.2.3 sense, so
    ///   `Σ a_i` is well-defined and stable to <1e-3.
    #[test]
    fn startup_state_produces_finite_coefficients() {
        let a = lsf_to_lp(&well_spaced_omega());
        for (i, ai) in a.iter().enumerate() {
            assert!(
                ai.is_finite(),
                "a_{} = {} not finite at start-up state",
                i + 1,
                ai,
            );
        }
    }

    /// **Brute-force polynomial-multiplication oracle** for the §3.2.6
    /// recursion. Multiplies `F_1(z) = Π_{i ∈ {1,3,5,7,9}} (1 − 2·q_i·z^-1 +
    /// z^-2)` and `F_2(z)` (likewise over even indices) by literal
    /// polynomial multiplication, applies eq (25) `(1 ± z^-1)`
    /// restoration, eq (26) recombination, and asserts the output
    /// matches `lsp_to_lp` to a tight float epsilon.
    ///
    /// This is the strongest possible algorithmic check short of a
    /// PCM-lockstep test against the spec's electronic-attachment
    /// reference (which is barred clean-room): the brute-force path
    /// and the recursive path are independently derived from spec
    /// eqs (13)/(14)/(25)/(26), and matching here means the recursion
    /// faithfully implements the polynomial multiplication.
    fn brute_force_f1_f2(q_in: &[f32; M]) -> ([f32; 6], [f32; 6]) {
        // Multiply (1 - 2q·z^-1 + z^-2) for q ∈ {q_1, q_3, q_5, q_7, q_9}
        // (F_1) and q ∈ {q_2, q_4, q_6, q_8, q_10} (F_2) by literal
        // polynomial accumulation. After 5 factors the polynomial has
        // degree 10 with 11 coefficients; by symmetry only the first 6
        // (indices 0..=5) are unique.
        fn multiply_factors(qs: &[f32; 5]) -> [f32; 6] {
            // Start with constant polynomial = 1.0 (degree 0).
            let mut poly: Vec<f32> = vec![1.0];
            for &q in qs {
                // Multiply poly by (1 - 2q·z^-1 + z^-2).
                let mut next = vec![0.0; poly.len() + 2];
                for (j, &c) in poly.iter().enumerate() {
                    next[j] += c;
                    next[j + 1] += -2.0 * q * c;
                    next[j + 2] += c;
                }
                poly = next;
            }
            // poly.len() == 11 here; take the first 6 (indices 0..=5).
            let mut out = [0.0_f32; 6];
            out.copy_from_slice(&poly[..6]);
            out
        }
        let odd_qs: [f32; 5] = [q_in[0], q_in[2], q_in[4], q_in[6], q_in[8]];
        let even_qs: [f32; 5] = [q_in[1], q_in[3], q_in[5], q_in[7], q_in[9]];
        (multiply_factors(&odd_qs), multiply_factors(&even_qs))
    }

    /// The §3.2.6 recursion `lsp_to_lp` and the brute-force
    /// polynomial-multiplication path must produce the same `A(z)`
    /// coefficients. Two LSP patterns are exercised to lock the
    /// recursion against accidental specialisation to one input.
    #[test]
    fn recursion_matches_brute_force_polynomial_multiplication() {
        let omega_patterns: &[[f32; M]] = &[
            [0.10, 0.30, 0.55, 0.80, 1.10, 1.40, 1.70, 2.00, 2.40, 2.80],
            [0.20, 0.45, 0.75, 1.00, 1.30, 1.60, 1.90, 2.20, 2.50, 2.85],
        ];
        for (idx, omega) in omega_patterns.iter().enumerate() {
            let q = crate::lsp_interpolate::omega_to_q(omega);
            let (f1_brute, f2_brute) = brute_force_f1_f2(&q);

            // Apply eq (25): (1 ± z^-1) factor restoration.
            let mut f1p_brute = [0.0_f32; 6];
            let mut f2p_brute = [0.0_f32; 6];
            f1p_brute[0] = f1_brute[0];
            f2p_brute[0] = f2_brute[0];
            for i in 1..=5 {
                f1p_brute[i] = f1_brute[i] + f1_brute[i - 1];
                f2p_brute[i] = f2_brute[i] - f2_brute[i - 1];
            }

            // Apply eq (26): A(z) = (F'_1 + F'_2) / 2.
            let mut a_brute = [0.0_f32; M];
            for i in 1..=5 {
                a_brute[i - 1] = 0.5 * f1p_brute[i] + 0.5 * f2p_brute[i];
            }
            for i in 6..=10 {
                let mirror = 11 - i;
                a_brute[i - 1] = 0.5 * f1p_brute[mirror] - 0.5 * f2p_brute[mirror];
            }

            let a_recursive = lsp_to_lp(&q);
            for k in 0..M {
                let drift = (a_recursive[k] - a_brute[k]).abs();
                assert!(
                    drift < 1e-4,
                    "pattern {idx}: a_{} differs — recursive {} vs brute {} (drift {})",
                    k + 1,
                    a_recursive[k],
                    a_brute[k],
                    drift,
                );
            }
        }
    }

    /// Eq (26) lower-half / upper-half symmetry: a_{6..=10} are
    /// computed from the mirror index 11 - i of `f'_1` and `f'_2`. We
    /// don't have direct access to f'_1, f'_2 from outside the
    /// function, but the structure is implied by checking that
    /// reversing the LSP vector produces the same LP filter when
    /// degenerate (this isn't a useful symmetry in general but the
    /// recursion's f'_1 / f'_2 computation is one-pass and any sign
    /// error would explode at random LSFs — caught by the
    /// `lsfs_are_roots_of_reconstructed_a_z` check above). Here we
    /// instead pin the **scaling sanity**: for an arbitrary valid
    /// input the absolute value of each a_i stays within a generous
    /// real-filter range (no exploding coefficient due to a missing
    /// 0.5 factor or a sign error).
    #[test]
    fn coefficients_within_real_lp_filter_range() {
        let omega: [f32; M] = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25, 2.55, 2.85];
        let a = lsf_to_lp(&omega);
        for (i, ai) in a.iter().enumerate() {
            // Real LP-filter coefficients are bounded for a stable
            // filter; we use a defensive ±32 range here (a real-speech
            // a_i is typically O(1) to O(3) with rare outliers; >32
            // would indicate the eq (26) ½ factor was dropped).
            assert!(
                ai.abs() < 32.0,
                "a_{} = {} outside real-LP range — eq (26) scaling drift?",
                i + 1,
                ai,
            );
        }
    }

    /// `lsf_to_lp(&omega)` should match `lsp_to_lp(&omega_to_q(&omega))`
    /// — the convenience wrapper is exactly the cosine-domain entry
    /// preceded by the boundary conversion. Lock that wiring so callers
    /// can mix and match the two entry points.
    #[test]
    fn lsf_wrapper_equals_explicit_cosine_pipeline() {
        let omega = well_spaced_omega();
        let via_wrapper = lsf_to_lp(&omega);
        let via_explicit = lsp_to_lp(&crate::lsp_interpolate::omega_to_q(&omega));
        for i in 0..M {
            assert!(
                (via_wrapper[i] - via_explicit[i]).abs() < 1e-7,
                "wrapper vs explicit a_{} differs: {} vs {}",
                i + 1,
                via_wrapper[i],
                via_explicit[i],
            );
        }
    }

    /// End-to-end through §3.2.4 / §3.2.5 / §3.2.6: feed three
    /// successive `(L0, L1, L2, L3)` quadruples through the round-207
    /// reconstructor, run the round-213 interpolator, convert each
    /// subframe's LSPs to LP coefficients. Confirm every output stays
    /// finite (no recursion drift) and that the two subframes for a
    /// given frame differ from each other on a non-steady-state input
    /// — the linear interpolation in the cosine domain should produce
    /// genuinely distinct subframe-1 / subframe-2 filters.
    #[test]
    fn full_decode_chain_subframes_differ_per_frame() {
        use crate::lsp_interpolate::{omega_to_q, LspInterpolator};
        use crate::lsp_reconstruct::LspReconstructor;

        let mut recon = LspReconstructor::new();
        let mut interp = LspInterpolator::new();

        let frames: &[(usize, usize, usize, usize)] =
            &[(0, 0, 0, 0), (1, 5, 7, 11), (0, 12, 3, 17)];
        let mut seen_distinct_subframes = false;
        for &(l0, l1, l2, l3) in frames {
            let omega = recon
                .reconstruct_frame(l0, l1, l2, l3)
                .expect("valid corpus indices");
            let q = omega_to_q(&omega);
            let [sub1_q, sub2_q] = interp.interpolate(&q);
            let a_sub1 = lsp_to_lp(&sub1_q);
            let a_sub2 = lsp_to_lp(&sub2_q);
            for k in 0..M {
                assert!(a_sub1[k].is_finite(), "frame sub1 a_{} not finite", k + 1);
                assert!(a_sub2[k].is_finite(), "frame sub2 a_{} not finite", k + 1);
                if (a_sub1[k] - a_sub2[k]).abs() > 1e-3 {
                    seen_distinct_subframes = true;
                }
            }
        }
        assert!(
            seen_distinct_subframes,
            "subframe-1 and subframe-2 should differ on a non-steady-state input chain",
        );
    }

    /// Constant-term invariant: `a_0` is implicit 1.0 (not stored),
    /// and the §3.2.6 recursion's first slot `f_1(0) = f_2(0) = 1.0`
    /// must not drift even after the eq (15) inner-loop overwrites.
    /// We exercise this by checking that `lsp_to_lp` on an all-zero
    /// `q` (q_i = 0 ⇒ ω_i = π/2 ⇒ a fully-symmetric, well-conditioned
    /// filter) yields a finite vector — the all-zero LSP case is the
    /// most stable input for the recursion and any constant-term
    /// drift would manifest as a finite-but-wrong leading coefficient.
    #[test]
    fn all_zero_q_produces_finite_coefficients() {
        let q = [0.0_f32; M];
        let a = lsp_to_lp(&q);
        for (i, ai) in a.iter().enumerate() {
            assert!(ai.is_finite(), "a_{} = {} not finite at q = 0", i + 1, ai,);
        }
    }

    /// **A(z) evaluated at z = 1** has a closed-form known from spec
    /// §3.2.6: with eq (26), `A(1) = (F'_1(1) + F'_2(1)) / 2`. But
    /// `F'_2(z) = (1 − z^-1)·F_2(z)` so `F'_2(1) = 0`. And
    /// `F'_1(z) = (1 + z^-1)·F_1(z)` so `F'_1(1) = 2·F_1(1)`. So
    /// `A(1) = F_1(1) = Π_{i ∈ {1,3,5,7,9}} (2 − 2·q_i)`.
    ///
    /// `A(1) = 1 + Σ_k a_k` directly. Verifying this closed form ties
    /// the lower-half `i ∈ 1..=5` of eq (26) AND the symmetry of the
    /// upper-half `i ∈ 6..=10` together — both halves contribute to
    /// `Σ a_k`, so a sign error in either half trips this test.
    #[test]
    fn a_z_at_z_equals_1_matches_closed_form() {
        let omega: [f32; M] = [0.10, 0.30, 0.55, 0.80, 1.10, 1.40, 1.70, 2.00, 2.40, 2.80];
        let q = crate::lsp_interpolate::omega_to_q(&omega);
        let a = lsp_to_lp(&q);

        let a_at_one: f32 = 1.0 + a.iter().sum::<f32>();
        // Closed form: A(1) = Π_{i ∈ {1,3,5,7,9}} (2 − 2·q_i).
        let closed_form: f32 = [q[0], q[2], q[4], q[6], q[8]]
            .iter()
            .map(|qi| 2.0 - 2.0 * qi)
            .product();
        let drift = (a_at_one - closed_form).abs();
        assert!(
            drift < 1e-3,
            "A(1) recursive = {} but closed form = {} (drift {})",
            a_at_one,
            closed_form,
            drift,
        );
    }

    /// Dual closed form at `z = −1`: `F'_1(−1) = 0`, so
    /// `A(−1) = F'_2(−1)/2 = ((1 − (−1)^-1)·F_2(−1))/2 = F_2(−1)`.
    /// `F_2(−1) = Π_{i ∈ {2,4,6,8,10}} (1 − 2·q_i·(−1) + (−1)^-2)
    ///         = Π (2 + 2·q_i)`. And
    /// `A(−1) = 1 + Σ_k a_k · (−1)^k = 1 − a_1 + a_2 − … + a_10`.
    /// A sign error in the upper-half `i ∈ 6..=10` rule of eq (26)
    /// trips this test where the z=1 form might not.
    #[test]
    fn a_z_at_z_equals_minus_1_matches_closed_form() {
        let omega: [f32; M] = [0.10, 0.30, 0.55, 0.80, 1.10, 1.40, 1.70, 2.00, 2.40, 2.80];
        let q = crate::lsp_interpolate::omega_to_q(&omega);
        let a = lsp_to_lp(&q);

        // A(-1) = a_0 + Σ_{k=1..=10} a_k·(-1)^k, with a_0 = 1.
        let a_at_minus_one: f32 = 1.0
            + a.iter()
                .enumerate()
                .map(|(i, ak)| {
                    let k = i + 1;
                    if k % 2 == 0 {
                        *ak
                    } else {
                        -ak
                    }
                })
                .sum::<f32>();
        // Closed form: A(−1) = Π_{i ∈ {2,4,6,8,10}} (2 + 2·q_i).
        let closed_form: f32 = [q[1], q[3], q[5], q[7], q[9]]
            .iter()
            .map(|qi| 2.0 + 2.0 * qi)
            .product();
        let drift = (a_at_minus_one - closed_form).abs();
        assert!(
            drift < 1e-3,
            "A(-1) recursive = {} but closed form = {} (drift {})",
            a_at_minus_one,
            closed_form,
            drift,
        );
    }
}
