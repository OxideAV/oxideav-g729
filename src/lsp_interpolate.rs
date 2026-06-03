//! §3.2.5 LSP-coefficient interpolation.
//!
//! Round 207 wired the §3.2.4 LSP-frame reconstruction algorithm,
//! producing one `ω̂^(m)` LSF vector per 10 ms frame. The G.729
//! decoder, however, runs in 5 ms **sub**-frames, so the two
//! subframes within a frame need their own LP filter. The spec
//! §3.2.5 ("Interpolation of the LSP coefficients") describes how
//! to obtain the per-subframe coefficients:
//!
//! * the **second subframe** of frame `m` uses the frame's own
//!   reconstructed coefficients `q_i^(current)`;
//! * the **first subframe** of frame `m` uses the linear midpoint
//!   between this frame and the previous frame's coefficients:
//!
//!   ```text
//!   q_i^(1) = 0.5 · q_i^(previous) + 0.5 · q_i^(current)        i = 1..=10
//!   q_i^(2) =                          q_i^(current)            i = 1..=10
//!   ```
//!
//!   (spec eq (24))
//!
//! Per spec §3.2.5 the interpolation is done **in the cosine
//! domain**: the operand `q_i` is `cos(ω̂_i)`, not `ω̂_i` itself.
//! Linear interpolation of `cos(·)` is *not* the same as linear
//! interpolation of `ω̂` followed by `cos(·)`, so the conversion
//! ordering matters and is locked here.
//!
//! ## API
//!
//! [`LspInterpolator`] holds the previous frame's `q_i^(previous)`
//! vector and produces the per-subframe interpolated cosine-domain
//! LSPs for each new current frame. At start-up the `previous`
//! vector is the cosine of the spec §3.2.4 start-up LSFs `ω̂_i =
//! i · π / 11`, matching the [`crate::lsp_reconstruct::LspReconstructor`]
//! start-up state so the two front-end stages line up bit-cleanly
//! on frame 1.
//!
//! [`omega_to_q`] / [`q_to_omega`] expose the boundary conversion
//! (`q_i = cos(ω̂_i)`, `ω̂_i = acos(q_i)`) for callers that need to
//! drive the interpolator from an existing LSF vector (e.g. from the
//! round-207 reconstructor) or convert the per-subframe result back
//! into the LSF domain.
//!
//! ## Wiring with [`crate::lsp_reconstruct`]
//!
//! Typical decoder flow per 10 ms frame:
//!
//! ```text
//! omega = reconstructor.reconstruct_frame(L0, L1, L2, L3)?;  // §3.2.4 → ω̂^(m)
//! q     = omega_to_q(&omega);                                // boundary
//! [sub1_q, sub2_q] = interpolator.interpolate(&q);           // §3.2.5 → q^(1), q^(2)
//! ```
//!
//! The two `q^(s)` cosine-domain vectors then drive the §3.2.6
//! LSP→LP conversion (one LP filter per subframe).

use crate::tables::M;

/// Number of subframes per 10 ms frame (§2.1) — a G.729 frame splits
/// into two 5 ms / 40-sample subframes.
pub const SUBFRAMES_PER_FRAME: usize = 2;

/// Convert an LSF vector `ω̂` in the normalised frequency domain
/// `[0, π]` into the cosine-domain `q_i = cos(ω̂_i)`.
///
/// This is the boundary conversion §3.2.5 / §3.2.6 operate on: the
/// per-subframe interpolation and the LSP→LP polynomial expansion
/// both take `q_i` as input.
#[must_use]
pub fn omega_to_q(omega: &[f32; M]) -> [f32; M] {
    let mut out = [0.0_f32; M];
    for (i, w) in omega.iter().enumerate() {
        out[i] = w.cos();
    }
    out
}

/// Convert a cosine-domain vector `q_i = cos(ω̂_i)` back to the LSF
/// domain `ω̂_i = acos(q_i)`.
///
/// `q_i` outside `[-1, 1]` is clamped before the `acos` call so that
/// numerical noise (rounding past the boundary by ~1e-7) does not
/// produce NaN. Real-valued G.729 LSPs satisfy `|q_i| ≤ cos(0.005)
/// < 1` after the §3.2.4 stability clamp, so the clamp is a
/// defence-in-depth check only.
#[must_use]
pub fn q_to_omega(q: &[f32; M]) -> [f32; M] {
    let mut out = [0.0_f32; M];
    for (i, qi) in q.iter().enumerate() {
        let clamped = qi.clamp(-1.0, 1.0);
        out[i] = clamped.acos();
    }
    out
}

/// §3.2.5 per-subframe LSP interpolator.
///
/// Carries the cosine-domain LSPs of the previous 10 ms frame
/// (`q_i^(previous)`). One [`Self::interpolate`] call per frame
/// returns the two subframes' interpolated cosine-domain LSPs and
/// advances `previous := current` so the next call sees the right
/// history.
///
/// ## Start-up
///
/// [`Self::new`] initialises the previous frame to
/// `cos(ω̂_i)` where `ω̂_i = i · π / 11` (i ∈ 1..=10). This matches
/// the spec §3.2.4 start-up vector used by
/// [`crate::lsp_reconstruct::LspReconstructor::new`], so the first
/// frame's subframe-1 interpolation degenerates cleanly when the
/// decoder begins on a steady-state signal.
#[derive(Debug, Clone)]
pub struct LspInterpolator {
    /// `q_i^(previous)` — cosine-domain LSPs of the frame before
    /// the most recent [`Self::interpolate`] call. Updated to
    /// `q_i^(current)` at the end of every call.
    previous_q: [f32; M],
}

impl Default for LspInterpolator {
    fn default() -> Self {
        Self::new()
    }
}

impl LspInterpolator {
    /// Build an interpolator whose previous-frame state is the
    /// cosine of the spec §3.2.4 start-up LSF vector
    /// (`ω̂_i = i · π / 11`).
    #[must_use]
    pub fn new() -> Self {
        let mut omega = [0.0_f32; M];
        for (i, slot) in omega.iter_mut().enumerate() {
            *slot = ((i + 1) as f32) * core::f32::consts::PI / 11.0;
        }
        Self {
            previous_q: omega_to_q(&omega),
        }
    }

    /// Borrow the stored `q_i^(previous)` for testing / inspection.
    #[must_use]
    pub fn previous_q(&self) -> &[f32; M] {
        &self.previous_q
    }

    /// Apply spec eq (24): given the current frame's cosine-domain
    /// LSPs `q_i^(current)`, return the per-subframe interpolated
    /// vectors `[q^(1), q^(2)]` where:
    ///
    /// * `q^(1)[i] = 0.5 · q^(previous)[i] + 0.5 · q^(current)[i]`
    /// * `q^(2)[i] =                         q^(current)[i]`
    ///
    /// Internal `previous_q` advances to `current_q`.
    pub fn interpolate(&mut self, current_q: &[f32; M]) -> [[f32; M]; SUBFRAMES_PER_FRAME] {
        let mut sub1 = [0.0_f32; M];
        let mut sub2 = [0.0_f32; M];
        for (i, cur) in current_q.iter().enumerate() {
            sub1[i] = 0.5 * self.previous_q[i] + 0.5 * cur;
            sub2[i] = *cur;
        }
        self.previous_q = *current_q;
        [sub1, sub2]
    }

    /// Convenience: drive [`Self::interpolate`] from an LSF-domain
    /// (`ω̂`) input — converts to the cosine domain via
    /// [`omega_to_q`] internally and returns LSF-domain subframe
    /// outputs via [`q_to_omega`].
    ///
    /// Useful for plumbing the round-207 reconstructor's `ω̂^(m)`
    /// straight through §3.2.5 when the caller stays in the LSF
    /// domain for inspection / debugging. The actual interpolation
    /// is still done in the cosine domain per spec §3.2.5.
    pub fn interpolate_from_omega(
        &mut self,
        current_omega: &[f32; M],
    ) -> [[f32; M]; SUBFRAMES_PER_FRAME] {
        let q_current = omega_to_q(current_omega);
        let [sub1_q, sub2_q] = self.interpolate(&q_current);
        [q_to_omega(&sub1_q), q_to_omega(&sub2_q)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f32::consts::PI;

    /// `omega_to_q` followed by `q_to_omega` round-trips a real LSF
    /// vector to within a tight float epsilon — the boundary helpers
    /// are honest inverses of one another (modulo the `acos` clamp).
    #[test]
    fn omega_q_round_trip_is_identity_on_real_lsfs() {
        let omega: [f32; M] = [0.10, 0.30, 0.60, 0.90, 1.20, 1.50, 1.80, 2.10, 2.40, 2.80];
        let recovered = q_to_omega(&omega_to_q(&omega));
        for i in 0..M {
            assert!(
                (recovered[i] - omega[i]).abs() < 1e-5,
                "round-trip drift at i={i}: {} vs {}",
                recovered[i],
                omega[i],
            );
        }
    }

    /// `q_to_omega` clamps inputs that are *just* outside `[-1, 1]`
    /// rather than emitting NaN. This is a numerical-noise defence:
    /// after stability-clamp the real `q_i` satisfies `|q_i| < 1`,
    /// but float rounding can land it at 1.0000001 — we want a
    /// well-defined `acos` boundary value (`0` or `π`) instead.
    #[test]
    fn q_to_omega_clamps_just_past_boundary() {
        let q: [f32; M] = [
            1.000_001, 0.5, 0.0, -0.5, -1.000_001, 0.99, -0.99, 0.1, -0.1, 0.0,
        ];
        let omega = q_to_omega(&q);
        assert!(omega[0].is_finite(), "q=1+eps → acos must not NaN");
        assert!(omega[4].is_finite(), "q=-1-eps → acos must not NaN");
        // acos(1) = 0, acos(-1) = π exactly.
        assert!((omega[0] - 0.0).abs() < 1e-5, "acos(1+eps) ≈ 0");
        assert!((omega[4] - PI).abs() < 1e-5, "acos(-1-eps) ≈ π");
    }

    /// At construction the interpolator's `previous_q` equals the
    /// cosine of the spec §3.2.4 start-up LSFs. This locks the
    /// stage-1-on-frame-1 boundary value.
    #[test]
    fn new_previous_q_matches_spec_startup_cosines() {
        let interp = LspInterpolator::new();
        for i in 0..M {
            let expected = (((i + 1) as f32) * PI / 11.0).cos();
            let got = interp.previous_q()[i];
            assert!(
                (got - expected).abs() < 1e-6,
                "previous_q[{i}] = {got}; expected cos((i+1)·π/11) = {expected}",
            );
        }
    }

    /// Subframe 2 of every frame is exactly the current frame's
    /// cosine-domain LSPs — eq (24) lower row.
    #[test]
    fn subframe_two_is_current_frame_q() {
        let mut interp = LspInterpolator::new();
        let current_q: [f32; M] = [
            0.95, 0.80, 0.60, 0.30, 0.10, -0.10, -0.30, -0.60, -0.80, -0.95,
        ];
        let [_, sub2] = interp.interpolate(&current_q);
        for i in 0..M {
            assert!(
                (sub2[i] - current_q[i]).abs() < 1e-7,
                "subframe-2[{i}] = {} not current_q[{i}] = {}",
                sub2[i],
                current_q[i],
            );
        }
    }

    /// Subframe 1 of every frame is the per-coordinate midpoint of
    /// the previous and current cosine-domain LSPs — eq (24) upper
    /// row.
    #[test]
    fn subframe_one_is_midpoint_in_cosine_domain() {
        let mut interp = LspInterpolator::new();
        let pre = *interp.previous_q();
        let current_q: [f32; M] = [
            0.95, 0.80, 0.60, 0.30, 0.10, -0.10, -0.30, -0.60, -0.80, -0.95,
        ];
        let [sub1, _] = interp.interpolate(&current_q);
        for i in 0..M {
            let expected = 0.5 * pre[i] + 0.5 * current_q[i];
            assert!(
                (sub1[i] - expected).abs() < 1e-7,
                "subframe-1[{i}] = {} not midpoint {} ((pre={} + cur={}) / 2)",
                sub1[i],
                expected,
                pre[i],
                current_q[i],
            );
        }
    }

    /// After one `interpolate` call, `previous_q` equals the just-
    /// consumed `current_q`. The next call's subframe-1 midpoint
    /// uses the new previous — so back-to-back identical frames
    /// produce a subframe-1 == subframe-2 result on the second call.
    #[test]
    fn previous_q_advances_to_current_q() {
        let mut interp = LspInterpolator::new();
        let current_q: [f32; M] = [
            0.95, 0.80, 0.60, 0.30, 0.10, -0.10, -0.30, -0.60, -0.80, -0.95,
        ];
        let _ = interp.interpolate(&current_q);
        for (i, cur) in current_q.iter().enumerate() {
            let prev_i = interp.previous_q()[i];
            assert!(
                (prev_i - cur).abs() < 1e-7,
                "previous_q[{i}] = {prev_i} not advanced to current_q[{i}] = {cur}",
            );
        }
        // Second call with the same current_q: subframe-1 midpoint
        // collapses onto current_q because previous == current.
        let [sub1, sub2] = interp.interpolate(&current_q);
        for i in 0..M {
            assert!(
                (sub1[i] - sub2[i]).abs() < 1e-7,
                "steady-state frame {i}: sub1 {} ≠ sub2 {}",
                sub1[i],
                sub2[i],
            );
        }
    }

    /// LSF-domain convenience entry point produces results equivalent
    /// to the cosine-domain interpolation followed by `q_to_omega`.
    /// This pins the wiring boundary between the round-207 LSF
    /// reconstructor and the §3.2.5 interpolator.
    #[test]
    fn interpolate_from_omega_matches_cosine_domain_pipeline() {
        let mut interp_omega = LspInterpolator::new();
        let mut interp_q = LspInterpolator::new();
        let current_omega: [f32; M] = [0.10, 0.30, 0.60, 0.90, 1.20, 1.50, 1.80, 2.10, 2.40, 2.80];
        let from_omega = interp_omega.interpolate_from_omega(&current_omega);

        let current_q = omega_to_q(&current_omega);
        let [sub1_q, sub2_q] = interp_q.interpolate(&current_q);
        let sub1_w = q_to_omega(&sub1_q);
        let sub2_w = q_to_omega(&sub2_q);

        for i in 0..M {
            assert!(
                (from_omega[0][i] - sub1_w[i]).abs() < 1e-5,
                "sub1 omega-vs-q[{i}]: {} vs {}",
                from_omega[0][i],
                sub1_w[i],
            );
            assert!(
                (from_omega[1][i] - sub2_w[i]).abs() < 1e-5,
                "sub2 omega-vs-q[{i}]: {} vs {}",
                from_omega[1][i],
                sub2_w[i],
            );
        }
    }

    /// End-to-end with the round-207 LSP reconstructor: feed three
    /// successive `(L0, L1, L2, L3)` quadruples through the
    /// reconstructor, plumb each `ω̂^(m)` through the interpolator,
    /// and confirm structural invariants of the §3.2.5 output:
    ///
    /// * each subframe-2 LSF equals the frame's reconstructed LSF
    ///   (eq (24) lower row);
    /// * each subframe-1 LSF, after cosine inversion, lies between
    ///   the previous and current frame's LSFs at every coordinate
    ///   (eq (24) upper row — midpoint in q means *between* in
    ///   omega because cos is monotone on `[0, π]`).
    #[test]
    fn lsp_reconstructor_drives_interpolator_consistently() {
        use crate::lsp_reconstruct::LspReconstructor;

        let mut recon = LspReconstructor::new();
        let mut interp = LspInterpolator::new();

        let mut prev_omega: Option<[f32; M]> = None;
        let frames: &[(usize, usize, usize, usize)] =
            &[(0, 0, 0, 0), (1, 5, 7, 11), (0, 12, 3, 17)];
        for &(l0, l1, l2, l3) in frames {
            let omega = recon
                .reconstruct_frame(l0, l1, l2, l3)
                .expect("valid corpus indices");
            let [sub1_w, sub2_w] = interp.interpolate_from_omega(&omega);

            for i in 0..M {
                // Subframe 2 == current frame.
                assert!(
                    (sub2_w[i] - omega[i]).abs() < 1e-5,
                    "sub2 not current: {} vs {}",
                    sub2_w[i],
                    omega[i],
                );
            }
            if let Some(prev) = prev_omega {
                for i in 0..M {
                    // cos is monotone-decreasing on [0, π]: cosine-
                    // domain midpoint translates to *between* in omega.
                    let lo = prev[i].min(omega[i]);
                    let hi = prev[i].max(omega[i]);
                    // Allow a small float epsilon either side.
                    assert!(
                        sub1_w[i] + 1e-4 >= lo && sub1_w[i] - 1e-4 <= hi,
                        "sub1[{i}] = {} not between prev {} and curr {}",
                        sub1_w[i],
                        prev[i],
                        omega[i],
                    );
                }
            }
            prev_omega = Some(omega);
        }
    }
}
