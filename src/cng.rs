//! ITU-T G.729 **Annex B §B.4.4** comfort-noise (CNG) **excitation gain
//! computation and excitation mixing** — the decode-side silence-frame
//! synthesis math that turns a decoded SID energy into the per-subframe
//! comfort-noise excitation `ex(n)`.
//!
//! This module covers the part of §B.4.4 that is **fully determined by
//! the Annex B equation prose** (eqs (B.19)–(B.26), read from the
//! in-repo Recommendation EPUB): the target-gain smoothing, the
//! adaptive/fixed gain solving, and the Gaussian-mixture scaling. It is
//! deliberately self-contained: it depends only on the **decoded SID
//! gain** (already produced by [`crate::annex_b::dequant_sid_energy_db`])
//! and the two per-subframe excitation shapes (unity-gain adaptive
//! `e_a(n)` and the random ACELP fixed `e_f(n)`), **not** on the §B.4.2.2
//! SID-LSP VQ subset codebooks (those tables are a documented docs-gap —
//! see [`crate::annex_b`]). The spectral envelope still routes through
//! the SID-LSP VQ once that gap is closed; the excitation **energy** path
//! wired here is independent of it.
//!
//! ## Spec source — Annex B §B.4.4 (06/2012 Recommendation)
//!
//! Clause B.4.4 (transcribed from the EPUB prose + equation rasters): at
//! the decoder the comfort noise is generated "by introducing a
//! pseudo-white excitation signal of controlled level into interpolated
//! LPC filters in the same manner [as] the decoder produces active
//! speech". The pseudo-white excitation `ex(n)` is a mixture of a
//! G.729-type excitation `ex1(n)` (small-gain adaptive + ACELP fixed) and
//! a white Gaussian excitation `ex2(n)`.
//!
//! ### Target gain — eq (B.19)
//!
//! Define the target excitation gain `G̃_t` (the square root of the
//! average energy the current frame's synthetic excitation must reach).
//! It is smoothed from the SID gain `G̃_sid`:
//!
//! ```text
//!          ⎧ G̃_sid                              if Vad_{t−1} = 1
//! G̃_t  =  ⎨
//!          ⎩ (7/8)·G̃_{t−1} + (1/8)·G̃_sid       otherwise
//! ```
//!
//! i.e. when the previous frame was **active** the target jumps straight
//! to the new SID gain; across a run of non-active frames it relaxes
//! toward `G̃_sid` with the `7/8 : 1/8` first-order smoothing.
//!
//! ### Per-subframe gain solving — eqs (B.20)–(B.22)
//!
//! The 80-sample frame is split into two 40-sample subframes. For each
//! subframe the CNG excitation is `Ga·e_a(n) + Gf·e_f(n)`, and the
//! adaptive/fixed gains `(Ga, Gf)` are chosen so the subframe average
//! energy equals `G̃_t²` (eq (B.20)):
//!
//! ```text
//!  (1/40)·Σ_{n=0}^{39} (Ga·e_a(n) + Gf·e_f(n))² = G̃_t²
//! ```
//!
//! With the abbreviations (read from the inline `Let us define` raster)
//!
//! ```text
//!   Ea = Σ_{n=0}^{39} e_a(n)²      I = Σ e_a(n)·e_f(n)      K = 40·G̃_t²
//! ```
//!
//! and the ACELP-structure identity `Σ_{n=0}^{39} e_f(n)² = 4` (the four
//! ±1 unit pulses), expanding eq (B.20) and dividing by `Σe_f² = 4`
//! gives a second-order equation in the fixed gain `Gf` for a **fixed**
//! random adaptive gain `Ga` (eq (B.21)):
//!
//! ```text
//!  Gf² + (Ga·I/2)·Gf + (Ea·Ga² − K)/4 = 0
//! ```
//!
//! The root of eq (B.21) with the **lowest absolute value** is taken as
//! `Gf` ("Notice that Gf can take a negative value"). To guarantee the
//! quadratic has a real root and to forbid large adaptive gains, `Ga` is
//! drawn randomly from (eq (B.22)):
//!
//! ```text
//!  [ 0 , Max{ 0.5, √(K/A) } ]     with  A = Ea − I²/4
//! ```
//!
//! The G.729-type excitation is then `ex1(n) = Ga·e_a(n) + Gf·e_f(n)`
//! (eq (B.23)).
//!
//! ### Gaussian mixture — eqs (B.24)–(B.26)
//!
//! With `E1 = Σ ex1²`, `E2 = Σ ex2²`, `E3 = Σ ex1·ex2` (eq (B.24)),
//! `α = 0.6` fixed, and `β > 0` the solution of (eq (B.25))
//!
//! ```text
//!  β²·E2 + 2·α·β·E3 + (α² − 1)·E1 = 0
//! ```
//!
//! the composite excitation is (eq (B.26))
//!
//! ```text
//!  ex(n) = α·ex1(n) + β·ex2(n)
//! ```
//!
//! If no positive `β` solves eq (B.25), `β` is set to `0` and `α` to `1`
//! (i.e. the mixture degenerates to the G.729-type excitation alone).
//!
//! ## Scope / clean-room note
//!
//! Everything in this module is sourced from the §B.4.4 prose and the
//! eq (B.19)–(B.26) rasters. The **random draws** the spec calls for
//! (the pitch lag in `[40, 103]`, the ACELP grid/sign/position selection
//! for `e_f(n)`, the `Ga` value in the eq (B.22) interval, and the
//! Gaussian `ex2(n)`) are supplied by the caller so this module stays a
//! pure, deterministic, testable arithmetic core; the spec's eq (96)
//! random generator already lives in [`crate::concealment`]. The
//! `I = Σ e_a(n)·e_f(n)` raster shows an upper summation index that does
//! not match the 40-sample subframe scope of the surrounding eqs (B.20),
//! `Ea`, and `Σe_f² = 4`; this module sums `I` over the subframe
//! (`n = 0..40`) consistently with those equations and flags the raster
//! index as a documented rendering ambiguity.

use crate::fixed_codebook::SUBFRAME_SIZE;

/// eq (B.19) active-to-CNG target-gain weight on the **new** SID gain
/// when relaxing across a run of non-active frames: `1/8`.
pub const CNG_SID_WEIGHT: f32 = 1.0 / 8.0;

/// eq (B.19) target-gain weight on the **previous** target across a run
/// of non-active frames: `7/8`.
pub const CNG_SMOOTH_WEIGHT: f32 = 7.0 / 8.0;

/// eq (B.26) fixed mixture proportion of the G.729-type excitation
/// `ex1(n)` — `α = 0.6` per clause B.4.4.
pub const ALPHA: f32 = 0.6;

/// eq (B.22) floor on the upper bound of the random adaptive-gain
/// interval: the upper bound is `Max{0.5, √(K/A)}`.
pub const GA_UPPER_FLOOR: f32 = 0.5;

/// The ACELP-structure fixed-excitation energy identity
/// `Σ_{n=0}^{39} e_f(n)² = 4` (four ±1 unit pulses), used to reduce
/// eq (B.20) to the eq (B.21) quadratic.
pub const FIXED_EXC_ENERGY: f32 = 4.0;

/// eq (B.19) target excitation gain `G̃_t`.
///
/// `g_sid` is the SID gain `G̃_sid` (the square root of the decoded SID
/// average energy). `prev_target` is the previous frame's target gain
/// `G̃_{t−1}`. `prev_active` is true iff the previous frame was an active
/// speech frame (`Vad_{t−1} = 1`), in which case the target jumps to
/// `g_sid`; otherwise it relaxes with the `7/8 : 1/8` smoothing.
#[must_use]
pub fn target_gain(g_sid: f32, prev_target: f32, prev_active: bool) -> f32 {
    if prev_active {
        g_sid
    } else {
        CNG_SMOOTH_WEIGHT * prev_target + CNG_SID_WEIGHT * g_sid
    }
}

/// The eq (B.20)/(B.21) abbreviations for one subframe.
///
/// - [`Self::ea`] = `Ea = Σ e_a(n)²` (adaptive-excitation energy).
/// - [`Self::i`] = `I = Σ e_a(n)·e_f(n)` (adaptive/fixed cross-energy,
///   summed over the 40-sample subframe — see the module docs on the
///   raster index ambiguity).
/// - [`Self::k`] = `K = 40·G̃_t²`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GainSolveTerms {
    /// `Ea = Σ_{n=0}^{39} e_a(n)²`.
    pub ea: f32,
    /// `I = Σ_{n=0}^{39} e_a(n)·e_f(n)`.
    pub i: f32,
    /// `K = 40·G̃_t²`.
    pub k: f32,
}

impl GainSolveTerms {
    /// Compute the eq (B.20)/(B.21) terms from the unity-gain adaptive
    /// excitation `e_a(n)`, the ACELP fixed excitation `e_f(n)`, and the
    /// eq (B.19) target gain `g_target = G̃_t`.
    #[must_use]
    pub fn new(e_a: &[f32; SUBFRAME_SIZE], e_f: &[f32; SUBFRAME_SIZE], g_target: f32) -> Self {
        let mut ea = 0.0f32;
        let mut i = 0.0f32;
        for (&a, &f) in e_a.iter().zip(e_f.iter()) {
            ea += a * a;
            i += a * f;
        }
        Self {
            ea,
            i,
            k: (SUBFRAME_SIZE as f32) * g_target * g_target,
        }
    }

    /// eq (B.22) upper bound of the random adaptive-gain interval:
    /// `Max{0.5, √(K/A)}` with `A = Ea − I²/4`.
    ///
    /// When `A ≤ 0` the `√(K/A)` term is undefined; the bound falls back
    /// to the [`GA_UPPER_FLOOR`] (`0.5`), which still admits a valid `Ga`
    /// draw (the eq (B.21) discriminant is handled in
    /// [`Self::solve_fixed_gain`]).
    #[must_use]
    pub fn ga_upper_bound(&self) -> f32 {
        let a = self.ea - self.i * self.i / 4.0;
        if a > 0.0 {
            (self.k / a).sqrt().max(GA_UPPER_FLOOR)
        } else {
            GA_UPPER_FLOOR
        }
    }

    /// Solve eq (B.21) for the fixed gain `Gf` given a chosen adaptive
    /// gain `Ga`, returning the root with the **lowest absolute value**
    /// ("Notice that Gf can take a negative value").
    ///
    /// eq (B.21): `Gf² + (Ga·I/2)·Gf + (Ea·Ga² − K)/4 = 0`.
    ///
    /// Returns `None` when the discriminant is negative (no real root) —
    /// the caller should then redraw `Ga` from the eq (B.22) interval (a
    /// `Ga` inside `[0, √(K/A)]` keeps the discriminant non-negative, so
    /// a valid draw always has a solution).
    #[must_use]
    pub fn solve_fixed_gain(&self, ga: f32) -> Option<f32> {
        // Monic quadratic Gf² + b·Gf + c = 0.
        let b = ga * self.i / 2.0;
        let c = (self.ea * ga * ga - self.k) / 4.0;
        let disc = b * b - 4.0 * c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let r1 = (-b + sqrt_disc) / 2.0;
        let r2 = (-b - sqrt_disc) / 2.0;
        Some(if r1.abs() <= r2.abs() { r1 } else { r2 })
    }
}

/// Build the eq (B.23) G.729-type excitation `ex1(n) = Ga·e_a(n) +
/// Gf·e_f(n)` for one subframe.
#[must_use]
pub fn build_ex1(
    e_a: &[f32; SUBFRAME_SIZE],
    e_f: &[f32; SUBFRAME_SIZE],
    ga: f32,
    gf: f32,
) -> [f32; SUBFRAME_SIZE] {
    let mut ex1 = [0.0f32; SUBFRAME_SIZE];
    for (out, (&a, &f)) in ex1.iter_mut().zip(e_a.iter().zip(e_f.iter())) {
        *out = ga * a + gf * f;
    }
    ex1
}

/// The eq (B.24) energies of the two mixture components.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixEnergies {
    /// `E1 = Σ ex1(n)²`.
    pub e1: f32,
    /// `E2 = Σ ex2(n)²` (≈ subframe size for a unit-variance `ex2`).
    pub e2: f32,
    /// `E3 = Σ ex1(n)·ex2(n)` (cross-energy).
    pub e3: f32,
}

impl MixEnergies {
    /// Compute the eq (B.24) energies from the two excitations.
    #[must_use]
    pub fn new(ex1: &[f32; SUBFRAME_SIZE], ex2: &[f32; SUBFRAME_SIZE]) -> Self {
        let mut e1 = 0.0f32;
        let mut e2 = 0.0f32;
        let mut e3 = 0.0f32;
        for (&x1, &x2) in ex1.iter().zip(ex2.iter()) {
            e1 += x1 * x1;
            e2 += x2 * x2;
            e3 += x1 * x2;
        }
        Self { e1, e2, e3 }
    }
}

/// The eq (B.25)/(B.26) mixture coefficients `(α, β)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixCoefficients {
    /// eq (B.26) `α` — `0.6` when a positive `β` solves eq (B.25),
    /// else `1.0` (degenerate fallback).
    pub alpha: f32,
    /// eq (B.26) `β > 0` — the Gaussian-component scale, or `0.0` when
    /// eq (B.25) has no positive root.
    pub beta: f32,
}

/// Solve eq (B.25) for the mixture coefficients `(α, β)` (eq (B.26)).
///
/// eq (B.25): `β²·E2 + 2·α·β·E3 + (α² − 1)·E1 = 0`, with `α = 0.6`
/// ([`ALPHA`]) and `β > 0`. If no positive `β` solves it, `β = 0` and
/// `α = 1` ("If no solution is found for β, it is set to 0 and α to 1").
#[must_use]
pub fn solve_mix(energies: &MixEnergies) -> MixCoefficients {
    let MixEnergies { e1, e2, e3 } = *energies;
    let fallback = MixCoefficients {
        alpha: 1.0,
        beta: 0.0,
    };

    // β²·E2 + (2·α·E3)·β + (α²−1)·E1 = 0. Degenerate to linear when
    // E2 ≈ 0 (ex2 carries no energy → no usable Gaussian mixture).
    let a = e2;
    let b = 2.0 * ALPHA * e3;
    let c = (ALPHA * ALPHA - 1.0) * e1;

    if a.abs() <= f32::EPSILON {
        // Linear b·β + c = 0.
        if b.abs() <= f32::EPSILON {
            return fallback;
        }
        let beta = -c / b;
        return if beta > 0.0 {
            MixCoefficients { alpha: ALPHA, beta }
        } else {
            fallback
        };
    }

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return fallback;
    }
    let sqrt_disc = disc.sqrt();
    // Prefer the larger root first; pick the smallest strictly-positive
    // root for the smoothest mixture.
    let r1 = (-b + sqrt_disc) / (2.0 * a);
    let r2 = (-b - sqrt_disc) / (2.0 * a);
    let positive = [r1, r2].into_iter().filter(|&r| r > 0.0).reduce(f32::min);
    match positive {
        Some(beta) => MixCoefficients { alpha: ALPHA, beta },
        None => fallback,
    }
}

/// Build the eq (B.26) composite CNG excitation
/// `ex(n) = α·ex1(n) + β·ex2(n)` for one subframe.
#[must_use]
pub fn build_cng_excitation(
    ex1: &[f32; SUBFRAME_SIZE],
    ex2: &[f32; SUBFRAME_SIZE],
    coeffs: MixCoefficients,
) -> [f32; SUBFRAME_SIZE] {
    let mut ex = [0.0f32; SUBFRAME_SIZE];
    for (out, (&x1, &x2)) in ex.iter_mut().zip(ex1.iter().zip(ex2.iter())) {
        *out = coeffs.alpha * x1 + coeffs.beta * x2;
    }
    ex
}

/// Stateful §B.4.4 CNG target-gain smoother (eq (B.19)).
///
/// Owns the previous-frame target gain `G̃_{t−1}` so a stream decoder can
/// advance it frame-by-frame across a non-active run. The first SID after
/// an active frame jumps straight to `G̃_sid`; later non-active frames
/// relax toward `G̃_sid`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CngGainSmoother {
    /// `G̃_{t−1}` — the previous frame's target excitation gain. Starts
    /// at `0.0` (no prior comfort-noise energy).
    prev_target: f32,
}

impl CngGainSmoother {
    /// Build a smoother with no prior comfort-noise energy
    /// (`G̃_{t−1} = 0`).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The carried previous-frame target gain `G̃_{t−1}`.
    #[must_use]
    pub fn prev_target(&self) -> f32 {
        self.prev_target
    }

    /// Advance one frame: compute the eq (B.19) target gain `G̃_t` from
    /// the SID gain `g_sid` and whether the previous frame was active,
    /// store it as the new `G̃_{t−1}`, and return it.
    pub fn next_target(&mut self, g_sid: f32, prev_active: bool) -> f32 {
        let g = target_gain(g_sid, self.prev_target, prev_active);
        self.prev_target = g;
        g
    }
}

/// eq (96) random-generator recurrence used by the §B.4.4 CNG excitation
/// for the pitch-lag, ACELP grid/sign/position, and `Ga` random draws.
///
/// Clause B.4.4: "to avoid desynchronization of the random generator used
/// to compute the excitation, the pseudo-random sequence reset is
/// performed at each active frame, both at the encoder and decoder
/// parts." The recurrence itself is the spec §4.4.4 generator
/// `seed = 31821·seed + 13849` (16-bit wrap), the same one
/// [`crate::concealment`] uses for erasure replacement; it is duplicated
/// here so the CNG path owns its own resettable sequence without reaching
/// into the concealer's private state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CngRandom {
    seed: u16,
}

/// eq (96) multiplier (`31821`).
const RNG_MULT: u16 = 31821;
/// eq (96) increment (`13849`).
const RNG_INC: u16 = 13849;

impl CngRandom {
    /// Build a generator from an explicit seed (the §4.4.4 start-up seed
    /// `21845` is the conventional choice; the stream decoder resets to it
    /// at each active frame per clause B.4.4).
    #[must_use]
    pub const fn new(seed: u16) -> Self {
        Self { seed }
    }

    /// The current seed (for inspection / determinism tests).
    #[must_use]
    pub const fn seed(&self) -> u16 {
        self.seed
    }

    /// Advance the eq (96) recurrence and return the new 16-bit value.
    pub fn next_u16(&mut self) -> u16 {
        self.seed = RNG_MULT.wrapping_mul(self.seed).wrapping_add(RNG_INC);
        self.seed
    }

    /// A uniform `f32` in `[0, 1)` from the next random word.
    pub fn next_unit(&mut self) -> f32 {
        f32::from(self.next_u16()) / 65536.0
    }

    /// A pitch lag uniformly in `[40, 103]` (clause B.4.4: "A pitch lag is
    /// randomly chosen in the interval [40, 103]").
    pub fn next_pitch_lag(&mut self) -> usize {
        40 + (usize::from(self.next_u16()) % 64)
    }

    /// An approximately unit-variance, zero-mean Gaussian sample for the
    /// `ex2(n)` white component (sum-of-12-uniforms central-limit
    /// approximation, `Σ U − 6`, the classic deterministic Gaussian
    /// shaping). `ex2` only needs to be "unit variance, zero mean" per
    /// clause B.4.4; the exact shaping is a free choice the spec does not
    /// pin, so a deterministic CLT draw keeps the path testable.
    pub fn next_gaussian(&mut self) -> f32 {
        let mut acc = 0.0f32;
        for _ in 0..12 {
            acc += self.next_unit();
        }
        acc - 6.0
    }
}

/// One subframe's §B.4.4 CNG inputs: the unity-gain adaptive excitation
/// `e_a(n)`, the ACELP fixed excitation `e_f(n)` (four ±1 pulses), and the
/// white Gaussian component `ex2(n)`.
#[derive(Debug, Clone, Copy)]
pub struct CngSubframeShapes {
    /// Unity-gain adaptive excitation `e_a(n)`.
    pub e_a: [f32; SUBFRAME_SIZE],
    /// ACELP fixed excitation `e_f(n)` (`Σ e_f² = 4`).
    pub e_f: [f32; SUBFRAME_SIZE],
    /// White Gaussian excitation `ex2(n)` (unit variance, zero mean).
    pub ex2: [f32; SUBFRAME_SIZE],
}

/// Build one subframe's §B.4.4 composite CNG excitation `ex(n)` end to
/// end (eqs (B.20)–(B.26)) given the per-subframe shapes, the eq (B.19)
/// target gain `g_target`, and a random draw `ga_draw ∈ [0, 1)` mapped
/// into the eq (B.22) `Ga` interval `[0, Max{0.5, √(K/A)}]`.
///
/// Returns the 40-sample composite excitation `ex(n)`. If the chosen `Ga`
/// somehow yields no real eq (B.21) root (it cannot for a draw inside the
/// eq (B.22) interval, but the guard keeps the function total), the
/// adaptive contribution is dropped (`Gf` solved at `Ga = 0`).
#[must_use]
pub fn synthesize_cng_subframe(
    shapes: &CngSubframeShapes,
    g_target: f32,
    ga_draw: f32,
) -> [f32; SUBFRAME_SIZE] {
    let terms = GainSolveTerms::new(&shapes.e_a, &shapes.e_f, g_target);
    let ga = ga_draw.clamp(0.0, 1.0) * terms.ga_upper_bound();
    let (ga, gf) = match terms.solve_fixed_gain(ga) {
        Some(gf) => (ga, gf),
        None => (0.0, terms.solve_fixed_gain(0.0).unwrap_or(0.0)),
    };
    let ex1 = build_ex1(&shapes.e_a, &shapes.e_f, ga, gf);
    let energies = MixEnergies::new(&ex1, &shapes.ex2);
    let coeffs = solve_mix(&energies);
    build_cng_excitation(&ex1, &shapes.ex2, coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// eq (B.19): after an active frame the target jumps straight to the
    /// SID gain (no smoothing).
    #[test]
    fn target_gain_jumps_to_sid_after_active() {
        let g = target_gain(100.0, 5.0, true);
        assert!((g - 100.0).abs() < 1e-6);
    }

    /// eq (B.19): across a non-active run the target relaxes with the
    /// `7/8 : 1/8` smoothing.
    #[test]
    fn target_gain_smooths_across_nonactive_run() {
        // G̃_{t−1} = 80, G̃_sid = 160 → 7/8·80 + 1/8·160 = 70 + 20 = 90.
        let g = target_gain(160.0, 80.0, false);
        assert!((g - 90.0).abs() < 1e-4, "g = {g}");
    }

    /// The smoother converges geometrically toward `G̃_sid` over a long
    /// non-active run (the 7/8 pole).
    #[test]
    fn smoother_converges_to_sid() {
        let mut s = CngGainSmoother::new();
        // First SID after active → jumps to sid.
        let first = s.next_target(200.0, true);
        assert!((first - 200.0).abs() < 1e-4);
        // A long run at a *new* sid relaxes toward it.
        for _ in 0..200 {
            s.next_target(50.0, false);
        }
        assert!((s.prev_target() - 50.0).abs() < 1e-2, "{}", s.prev_target());
    }

    /// Helper: a 4-pulse ±1 ACELP fixed excitation with `Σe_f² = 4`.
    fn unit_fixed() -> [f32; SUBFRAME_SIZE] {
        let mut f = [0.0f32; SUBFRAME_SIZE];
        f[0] = 1.0;
        f[6] = -1.0;
        f[17] = 1.0;
        f[34] = -1.0;
        f
    }

    /// `Σe_f² = 4` for the canonical 4-pulse excitation (the eq (B.20)→
    /// (B.21) reduction identity).
    #[test]
    fn fixed_excitation_energy_is_four() {
        let f = unit_fixed();
        let e: f32 = f.iter().map(|x| x * x).sum();
        assert!((e - FIXED_EXC_ENERGY).abs() < 1e-6);
    }

    /// eq (B.21) solving: the returned `(Ga, Gf)` reproduces the eq
    /// (B.20) target subframe energy `G̃_t²` exactly.
    #[test]
    fn solved_gains_hit_target_energy() {
        let e_a: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| ((n as f32) * 0.1).sin());
        let e_f = unit_fixed();
        let g_target = 30.0;
        let terms = GainSolveTerms::new(&e_a, &e_f, g_target);
        // Pick a Ga inside the eq (B.22) interval.
        let ga = 0.4 * terms.ga_upper_bound();
        let gf = terms
            .solve_fixed_gain(ga)
            .expect("real root for in-range Ga");
        let ex1 = build_ex1(&e_a, &e_f, ga, gf);
        let avg_energy: f32 = ex1.iter().map(|x| x * x).sum::<f32>() / (SUBFRAME_SIZE as f32);
        assert!(
            (avg_energy - g_target * g_target).abs() < 1e-2,
            "avg energy {avg_energy} target² {}",
            g_target * g_target
        );
    }

    /// eq (B.21) root selection: the chosen `Gf` is the one with the
    /// lower absolute value of the two quadratic roots.
    #[test]
    fn solve_fixed_gain_picks_lowest_abs_root() {
        let e_a: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let e_f = unit_fixed();
        let terms = GainSolveTerms::new(&e_a, &e_f, 25.0);
        let ga = 0.3;
        let gf = terms.solve_fixed_gain(ga).expect("real root");
        // Recompute both roots and confirm the chosen one is min-abs.
        let b = ga * terms.i / 2.0;
        let c = (terms.ea * ga * ga - terms.k) / 4.0;
        let disc = b * b - 4.0 * c;
        let sq = disc.sqrt();
        let r1 = (-b + sq) / 2.0;
        let r2 = (-b - sq) / 2.0;
        let expected = if r1.abs() <= r2.abs() { r1 } else { r2 };
        assert!((gf - expected).abs() < 1e-5);
    }

    /// eq (B.22): a `Ga` drawn inside `[0, ga_upper_bound]` always yields
    /// a real `Gf` root (the discriminant stays non-negative).
    #[test]
    fn ga_in_eq_b22_interval_always_solvable() {
        let e_a: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| ((n * 3) % 11) as f32 - 5.0);
        let e_f = unit_fixed();
        let terms = GainSolveTerms::new(&e_a, &e_f, 40.0);
        let upper = terms.ga_upper_bound();
        for k in 0..=10 {
            let ga = upper * (k as f32) / 10.0;
            assert!(
                terms.solve_fixed_gain(ga).is_some(),
                "Ga = {ga} (upper {upper}) gave no real root"
            );
        }
    }

    /// eq (B.25): a valid Gaussian mixture preserves the eq (B.20) target
    /// energy — `Σ ex(n)² = α²·E1 + 2αβ·E3 + β²·E2 = E1` (the eq (B.25)
    /// constraint rearranges to exactly this).
    #[test]
    fn mixture_preserves_ex1_energy() {
        let ex1: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| ((n as f32) * 0.2).cos() * 10.0);
        // Unit-variance-ish zero-mean ex2.
        let ex2: [f32; SUBFRAME_SIZE] =
            core::array::from_fn(|n| if n % 2 == 0 { 1.0 } else { -1.0 });
        let energies = MixEnergies::new(&ex1, &ex2);
        let coeffs = solve_mix(&energies);
        let ex = build_cng_excitation(&ex1, &ex2, coeffs);
        let e_mix: f32 = ex.iter().map(|x| x * x).sum();
        assert!(
            (e_mix - energies.e1).abs() < 1e-1,
            "mixed energy {e_mix} vs E1 {}",
            energies.e1
        );
    }

    /// eq (B.25) fallback: when `ex2` carries no energy the mixture
    /// degenerates to `α = 1, β = 0` (the G.729-type excitation alone).
    #[test]
    fn mix_falls_back_when_ex2_is_zero() {
        let ex1: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let ex2 = [0.0f32; SUBFRAME_SIZE];
        let energies = MixEnergies::new(&ex1, &ex2);
        let coeffs = solve_mix(&energies);
        assert!((coeffs.alpha - 1.0).abs() < 1e-6);
        assert!(coeffs.beta.abs() < 1e-6);
    }

    /// β is non-negative by construction (eq (B.25) requires `β > 0`;
    /// the fallback yields `β = 0`).
    #[test]
    fn beta_is_non_negative() {
        for seed in 0..32 {
            let ex1: [f32; SUBFRAME_SIZE] =
                core::array::from_fn(|n| (((n + seed) as f32 * 0.37).sin()) * (seed as f32 + 1.0));
            let ex2: [f32; SUBFRAME_SIZE] =
                core::array::from_fn(|n| ((n as f32 * 1.3 + seed as f32).cos()).signum());
            let coeffs = solve_mix(&MixEnergies::new(&ex1, &ex2));
            assert!(coeffs.beta >= 0.0, "β = {} (seed {seed})", coeffs.beta);
            assert!(coeffs.alpha == ALPHA || coeffs.alpha == 1.0);
        }
    }

    /// The constants match the clause B.4.4 spec values.
    #[test]
    fn constants_match_spec() {
        assert!((ALPHA - 0.6).abs() < 1e-9);
        assert!((CNG_SMOOTH_WEIGHT - 0.875).abs() < 1e-9);
        assert!((CNG_SID_WEIGHT - 0.125).abs() < 1e-9);
        assert!((GA_UPPER_FLOOR - 0.5).abs() < 1e-9);
        assert!((FIXED_EXC_ENERGY - 4.0).abs() < 1e-9);
    }

    /// The eq (96) CNG random generator matches the spec recurrence
    /// `seed = 31821·seed + 13849` (16-bit wrap) from seed 21845.
    #[test]
    fn cng_random_matches_eq96_recurrence() {
        let mut r = CngRandom::new(21845);
        let expect = 31821u16.wrapping_mul(21845).wrapping_add(13849);
        assert_eq!(r.next_u16(), expect);
    }

    /// The random pitch lag always lands in `[40, 103]` (clause B.4.4).
    #[test]
    fn cng_pitch_lag_in_spec_interval() {
        let mut r = CngRandom::new(1);
        for _ in 0..10_000 {
            let lag = r.next_pitch_lag();
            assert!((40..=103).contains(&lag), "lag {lag} out of [40, 103]");
        }
    }

    /// The deterministic Gaussian draw is approximately zero-mean and
    /// unit-variance (clause B.4.4: `ex2` has unit variance, zero mean).
    #[test]
    fn cng_gaussian_is_zero_mean_unit_variance() {
        let mut r = CngRandom::new(12345);
        let n = 50_000;
        let mut sum = 0.0f64;
        let mut sumsq = 0.0f64;
        for _ in 0..n {
            let g = f64::from(r.next_gaussian());
            sum += g;
            sumsq += g * g;
        }
        let mean = sum / n as f64;
        let var = sumsq / n as f64 - mean * mean;
        assert!(mean.abs() < 0.05, "mean {mean}");
        assert!((var - 1.0).abs() < 0.1, "var {var}");
    }

    /// The CNG generator is deterministic — same seed → same sequence.
    #[test]
    fn cng_random_is_deterministic() {
        let mut a = CngRandom::new(21845);
        let mut b = CngRandom::new(21845);
        for _ in 0..100 {
            assert_eq!(a.next_u16(), b.next_u16());
        }
    }

    /// End-to-end subframe synthesis: the composite excitation energy
    /// tracks the eq (B.19) target subframe energy `G̃_t²` (the eq (B.25)
    /// mixture preserves the eq (B.20) energy `ex1` was solved to hit).
    #[test]
    fn synthesize_subframe_hits_target_energy() {
        let mut rng = CngRandom::new(777);
        let g_target = 35.0;
        let e_a: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| ((n as f32) * 0.13).sin() * 2.0);
        let mut e_f = [0.0f32; SUBFRAME_SIZE];
        e_f[2] = 1.0;
        e_f[11] = -1.0;
        e_f[23] = 1.0;
        e_f[38] = -1.0;
        let ex2: [f32; SUBFRAME_SIZE] = core::array::from_fn(|_| rng.next_gaussian());
        let shapes = CngSubframeShapes { e_a, e_f, ex2 };
        let ex = synthesize_cng_subframe(&shapes, g_target, 0.5);
        let avg: f32 = ex.iter().map(|x| x * x).sum::<f32>() / (SUBFRAME_SIZE as f32);
        // The mixture preserves ex1's energy, which was solved to G̃_t².
        assert!(
            (avg - g_target * g_target).abs() / (g_target * g_target) < 0.05,
            "avg {avg} target² {}",
            g_target * g_target
        );
    }

    /// A zero target gain yields a (near-)silent excitation.
    #[test]
    fn synthesize_subframe_zero_target_is_silent() {
        let mut rng = CngRandom::new(9);
        let e_a: [f32; SUBFRAME_SIZE] = core::array::from_fn(|n| (n as f32) - 20.0);
        let mut e_f = [0.0f32; SUBFRAME_SIZE];
        e_f[0] = 1.0;
        e_f[5] = -1.0;
        e_f[10] = 1.0;
        e_f[15] = -1.0;
        let ex2: [f32; SUBFRAME_SIZE] = core::array::from_fn(|_| rng.next_gaussian());
        let shapes = CngSubframeShapes { e_a, e_f, ex2 };
        let ex = synthesize_cng_subframe(&shapes, 0.0, 0.5);
        let e: f32 = ex.iter().map(|x| x * x).sum();
        assert!(e < 1e-3, "expected silence, got energy {e}");
    }
}
