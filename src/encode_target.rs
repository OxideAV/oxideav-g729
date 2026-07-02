//! §3.5 **impulse response** + §3.6 **target signal** computation — the
//! two per-subframe pre-computations the closed-loop adaptive- and
//! fixed-codebook searches operate on.
//!
//! Eighth encoder-side stage. Consumes the quantised LP coefficients
//! `â_i` (from the [`crate::lsp_quantize`] →
//! [`crate::lsp_interpolate`] → [`crate::lsp_to_lp`] chain) and the
//! unquantised `a_i` + `(γ₁, γ₂)` pair (from [`crate::levinson`] /
//! [`crate::perceptual_weighting`]).
//!
//! ## Spec source — clauses 3.5, 3.6, equation (36)
//!
//! * **§3.5** (impulse response): "The impulse response `h(n)` of the
//!   weighted synthesis filter `W(z)/Â(z)` is needed for the search of
//!   adaptive and fixed codebooks. The impulse response `h(n)` is
//!   computed for each subframe by filtering a signal consisting of the
//!   coefficients of the filter `A(z/γ₁)` extended by zeros through the
//!   two filters `1/Â(z)` and `1/A(z/γ₂)`."
//! * **§3.6** (target signal): the target `x(n)` for the
//!   adaptive-codebook search is the weighted speech minus the
//!   zero-input response of `W(z)/Â(z)`; "an equivalent procedure …
//!   used in this Recommendation, is the filtering of the LP residual
//!   signal `r(n)` through the combination of synthesis filter `1/Â(z)`
//!   and the weighting filter `A(z/γ₁)/A(z/γ₂)`. After determining the
//!   excitation for the subframe, the initial states of these filters
//!   are updated by filtering the difference between the residual and
//!   excitation signals."
//! * **eq (36)** (LP residual):
//!   `r(n) = s(n) + Σ_{i=1}^{10} â_i·s(n−i)`, `n = 0 … 39`.
//!
//! ## Filter-chain state
//!
//! The composite `W(z)/Â(z) = A(z/γ₁) / [Â(z)·A(z/γ₂)]` is realised as
//! three cascaded stages, each with its own 10-tap memory:
//!
//! 1. all-pole `1/Â(z)` (the synthesis filter, quantised coefficients),
//! 2. FIR `A(z/γ₁)` (unquantised, γ₁-scaled),
//! 3. all-pole `1/A(z/γ₂)` (unquantised, γ₂-scaled).
//!
//! [`TargetFilter::target`] runs the residual through the chain
//! *without* committing the memories (the target computation must not
//! disturb the states); [`TargetFilter::update`] then filters the
//! difference `r(n) − u(n)` (residual minus chosen excitation) through
//! the same chain and commits the resulting memories, per §3.6/§3.10.
//! When the excitation matches the residual perfectly the difference is
//! zero and the memories stay zero — the property the tests pin.

use crate::tables::M;

/// Samples per subframe (clause 2.1).
pub const L_SUBFR: usize = 40;

/// Per-subframe filter coefficient bundle for the weighted synthesis
/// chain: quantised `â`, unquantised `a`, and the adaptive `(γ₁, γ₂)`.
#[derive(Debug, Clone, Copy)]
pub struct WeightedSynthesisCoefs {
    /// Quantised LP coefficients `â_i` (`a_hat[i−1] = â_i`).
    pub a_hat: [f32; M],
    /// Unquantised LP coefficients `a_i`.
    pub a: [f32; M],
    /// Perceptual-weighting numerator factor `γ₁`.
    pub gamma1: f32,
    /// Perceptual-weighting denominator factor `γ₂`.
    pub gamma2: f32,
}

impl WeightedSynthesisCoefs {
    /// The γ-scaled FIR rows `a_i·γ₁ⁱ` and `a_i·γ₂ⁱ`.
    fn scaled(&self) -> ([f32; M], [f32; M]) {
        let mut num = [0.0f32; M];
        let mut den = [0.0f32; M];
        let (mut g1, mut g2) = (1.0f32, 1.0f32);
        for i in 0..M {
            g1 *= self.gamma1;
            g2 *= self.gamma2;
            num[i] = self.a[i] * g1;
            den[i] = self.a[i] * g2;
        }
        (num, den)
    }
}

/// §3.5: the 40-sample impulse response `h(n)` of the weighted
/// synthesis filter `W(z)/Â(z) = A(z/γ₁)/[Â(z)·A(z/γ₂)]`, computed with
/// all-zero filter states (the spec's "extended by zeros" input is the
/// `A(z/γ₁)` coefficient sequence, i.e. the FIR numerator applied to a
/// unit impulse).
#[must_use]
pub fn impulse_response(coefs: &WeightedSynthesisCoefs) -> [f32; L_SUBFR] {
    let (num, den) = coefs.scaled();

    // Input signal: the coefficients of A(z/γ₁) = [1, num…] extended by
    // zeros to the subframe length.
    let mut h = [0.0f32; L_SUBFR];
    h[0] = 1.0;
    h[1..=M].copy_from_slice(&num);

    // Through 1/Â(z): h(n) −= Σ â_i·h(n−i).
    for n in 0..L_SUBFR {
        for i in 0..M.min(n) {
            h[n] -= coefs.a_hat[i] * h[n - 1 - i];
        }
    }
    // Through 1/A(z/γ₂): h(n) −= Σ a_i γ₂ⁱ·h(n−i).
    for n in 0..L_SUBFR {
        for i in 0..M.min(n) {
            h[n] -= den[i] * h[n - 1 - i];
        }
    }
    h
}

/// eq (36): the LP residual of one subframe.
/// `r(n) = s(n) + Σ â_i·s(n−i)`; `s_past` carries `s(−1) … s(−10)`
/// (most recent first).
#[must_use]
pub fn lp_residual(s: &[f32; L_SUBFR], s_past: &[f32; M], a_hat: &[f32; M]) -> [f32; L_SUBFR] {
    let mut r = [0.0f32; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = s[n];
        for i in 0..M {
            let sv = if n > i { s[n - 1 - i] } else { s_past[i - n] };
            acc += a_hat[i] * sv;
        }
        r[n] = acc;
    }
    r
}

/// The three cascaded 10-tap memories of the target-computation chain
/// (§3.6): `1/Â(z)` → `A(z/γ₁)` → `1/A(z/γ₂)`.
#[derive(Debug, Clone, Default)]
pub struct TargetFilter {
    /// `1/Â(z)` output memory (most recent first).
    syn_mem: [f32; M],
    /// `A(z/γ₁)` input memory (most recent first).
    fir_mem: [f32; M],
    /// `1/A(z/γ₂)` output memory (most recent first).
    den_mem: [f32; M],
}

impl TargetFilter {
    /// Fresh chain with all-zero memories (clause 4.3).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Runs `input` through the three-stage chain starting from the
    /// supplied memories, returning the output and the advanced
    /// memories.
    fn run(
        &self,
        input: &[f32; L_SUBFR],
        coefs: &WeightedSynthesisCoefs,
    ) -> ([f32; L_SUBFR], Self) {
        let (num, den) = coefs.scaled();

        // Stage 1: all-pole 1/Â(z).
        let mut s1 = [0.0f32; L_SUBFR];
        for n in 0..L_SUBFR {
            let mut acc = input[n];
            for i in 0..M {
                let y = if n > i {
                    s1[n - 1 - i]
                } else {
                    self.syn_mem[i - n]
                };
                acc -= coefs.a_hat[i] * y;
            }
            s1[n] = acc;
        }
        // Stage 2: FIR A(z/γ₁).
        let mut s2 = [0.0f32; L_SUBFR];
        for n in 0..L_SUBFR {
            let mut acc = s1[n];
            for i in 0..M {
                let x = if n > i {
                    s1[n - 1 - i]
                } else {
                    self.fir_mem[i - n]
                };
                acc += num[i] * x;
            }
            s2[n] = acc;
        }
        // Stage 3: all-pole 1/A(z/γ₂).
        let mut s3 = [0.0f32; L_SUBFR];
        for n in 0..L_SUBFR {
            let mut acc = s2[n];
            for i in 0..M {
                let y = if n > i {
                    s3[n - 1 - i]
                } else {
                    self.den_mem[i - n]
                };
                acc -= den[i] * y;
            }
            s3[n] = acc;
        }

        let next = Self {
            syn_mem: std::array::from_fn(|i| s1[L_SUBFR - 1 - i]),
            fir_mem: std::array::from_fn(|i| s1[L_SUBFR - 1 - i]),
            den_mem: std::array::from_fn(|i| s3[L_SUBFR - 1 - i]),
        };
        (s3, next)
    }

    /// §3.6: computes the target signal `x(n)` by filtering the LP
    /// residual through the chain. The internal memories are **not**
    /// modified — the states advance only via [`Self::update`].
    #[must_use]
    pub fn target(
        &self,
        residual: &[f32; L_SUBFR],
        coefs: &WeightedSynthesisCoefs,
    ) -> [f32; L_SUBFR] {
        self.run(residual, coefs).0
    }

    /// §3.6/§3.10 memory update: filters the difference between the
    /// residual and the chosen subframe excitation `u(n)` through the
    /// chain and commits the resulting memories.
    pub fn update(
        &mut self,
        residual: &[f32; L_SUBFR],
        excitation: &[f32; L_SUBFR],
        coefs: &WeightedSynthesisCoefs,
    ) {
        let diff: [f32; L_SUBFR] = std::array::from_fn(|n| residual[n] - excitation[n]);
        let (_, next) = self.run(&diff, coefs);
        *self = next;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_coefs() -> WeightedSynthesisCoefs {
        WeightedSynthesisCoefs {
            a_hat: [-0.6, 0.3, -0.1, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            a: [-0.55, 0.28, -0.09, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            gamma1: 0.94,
            gamma2: 0.6,
        }
    }

    /// Trivial filters (`a = â = 0`): `W/Â = 1`, so `h` is the unit
    /// impulse.
    #[test]
    fn identity_impulse_response() {
        let coefs = WeightedSynthesisCoefs {
            a_hat: [0.0; M],
            a: [0.0; M],
            gamma1: 0.94,
            gamma2: 0.6,
        };
        let h = impulse_response(&coefs);
        assert!((h[0] - 1.0).abs() < 1e-6);
        for &v in &h[1..] {
            assert!(v.abs() < 1e-6);
        }
    }

    /// With `γ₁ = γ₂` and `a = â`, `W(z)/Â(z)` collapses to `1/Â(z)`;
    /// `h` must match the direct all-pole recursion.
    #[test]
    fn collapses_to_synthesis_filter() {
        let a: [f32; M] = [-0.5, 0.2, -0.08, 0.03, -0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
        let coefs = WeightedSynthesisCoefs {
            a_hat: a,
            a,
            gamma1: 0.7,
            gamma2: 0.7,
        };
        let h = impulse_response(&coefs);
        // Direct 1/Â(z) impulse response.
        let mut want = [0.0f32; L_SUBFR];
        want[0] = 1.0;
        for n in 0..L_SUBFR {
            for i in 0..M.min(n) {
                want[n] -= a[i] * want[n - 1 - i];
            }
        }
        for n in 0..L_SUBFR {
            assert!(
                (h[n] - want[n]).abs() < 1e-4,
                "h[{n}]={} want {}",
                h[n],
                want[n]
            );
        }
    }

    /// eq (36) inverts the synthesis filter: synthesising a signal from
    /// a known excitation through `1/Â(z)` and then taking the LP
    /// residual recovers the excitation.
    #[test]
    fn residual_inverts_synthesis() {
        let a_hat: [f32; M] = [-0.7, 0.4, -0.2, 0.1, -0.05, 0.02, 0.0, 0.0, 0.0, 0.0];
        // Build 50 samples of synthetic speech from a pseudo-random
        // excitation (10 warm-up + 40 subframe).
        let u: [f32; 50] = std::array::from_fn(|n| ((n * 17 % 41) as f32) - 20.0);
        let mut s = [0.0f32; 50];
        for n in 0..50 {
            let mut acc = u[n];
            for i in 0..M.min(n) {
                acc -= a_hat[i] * s[n - 1 - i];
            }
            s[n] = acc;
        }
        let sub: [f32; L_SUBFR] = std::array::from_fn(|n| s[10 + n]);
        let past: [f32; M] = std::array::from_fn(|i| s[9 - i]);
        let r = lp_residual(&sub, &past, &a_hat);
        for n in 0..L_SUBFR {
            assert!(
                (r[n] - u[10 + n]).abs() < 1e-3,
                "residual[{n}]={} want {}",
                r[n],
                u[10 + n]
            );
        }
    }

    /// With zero memories the target equals the convolution of the
    /// residual with the §3.5 impulse response (the chain is LTI).
    #[test]
    fn target_is_convolution_from_zero_state() {
        let coefs = test_coefs();
        let f = TargetFilter::new();
        let r: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 7 % 23) as f32) - 11.0);
        let x = f.target(&r, &coefs);
        let h = impulse_response(&coefs);
        for n in 0..L_SUBFR {
            let mut conv = 0.0f32;
            for k in 0..=n {
                conv += r[k] * h[n - k];
            }
            assert!(
                (x[n] - conv).abs() < 1e-3 * (1.0 + conv.abs()),
                "x[{n}]={} conv={conv}",
                x[n]
            );
        }
    }

    /// Perfect excitation (`u = r`) leaves the memories zero: the next
    /// subframe's target is again the pure convolution.
    #[test]
    fn perfect_excitation_keeps_zero_state() {
        let coefs = test_coefs();
        let mut f = TargetFilter::new();
        let r: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 11 % 31) as f32) - 15.0);
        f.update(&r, &r, &coefs);
        // Now the target of a fresh residual must equal the zero-state
        // convolution.
        let r2: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 5 % 17) as f32) - 8.0);
        let x = f.target(&r2, &coefs);
        let x_zero = TargetFilter::new().target(&r2, &coefs);
        for n in 0..L_SUBFR {
            assert!(
                (x[n] - x_zero[n]).abs() < 1e-4,
                "state leaked at n={n}: {} vs {}",
                x[n],
                x_zero[n]
            );
        }
    }

    /// `target` must not mutate the filter state: two consecutive calls
    /// give identical output.
    #[test]
    fn target_is_pure() {
        let coefs = test_coefs();
        let mut f = TargetFilter::new();
        let r: [f32; L_SUBFR] = std::array::from_fn(|n| (n as f32) - 20.0);
        // Disturb the state first so the check is non-trivial.
        let u = [0.0f32; L_SUBFR];
        f.update(&r, &u, &coefs);
        let x1 = f.target(&r, &coefs);
        let x2 = f.target(&r, &coefs);
        assert_eq!(x1, x2);
    }
}
