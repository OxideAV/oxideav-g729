//! §3.5 **impulse response** and §3.6 / §3.10 **target signal** chain
//! on the fixed grid.
//!
//! ## §3.5
//!
//! `h(n)` is the response of `W(z)/Â(z) = A(z/γ₁)/[Â(z)·A(z/γ₂)]`:
//! the Q12 coefficients of `A(z/γ₁)` extended by zeros, run through
//! `1/Â(z)` then `1/A(z/γ₂)` from zero state. Both all-pole runs use
//! the shared `<< 3` / round landing, so `h` sits on Q12 (`h(0) =
//! 4096` for the identity).
//!
//! ## §3.6 / §3.10
//!
//! The target `x(n)` is the LP residual filtered through `1/Â(z)`,
//! `A(z/γ₁)` and `1/A(z/γ₂)` with the memories the clause-3.10 update
//! leaves behind: the synthesis filter's state is the reconstruction
//! error `e(n) = s(n) − ŝ(n)` (it also feeds the FIR `A(z/γ₁)`), and
//! the weighting denominator's state is the weighted error `ew(n) =
//! x(n) − ĝ_p·y(n) − ĝ_c·z(n)` of eq (76), both over `n = 30 … 39`.
//!
//! ## Number grids
//!
//! Signals on Q0; `y(n)` (filtered adaptive vector) on Q0; `z(n)`
//! (filtered codevector) on Q12; `ĝ_p` Q14, `ĝ_c` Q1 — the eq (76)
//! products land on the doubled Q15 grid exactly like the decoder's
//! eq (75) excitation build, then round to Q0.

use crate::fx::filters::{residu, syn_filt_mem, L_SUBFR};
use crate::fx::ops::{l_add, l_mult, l_shl, round, sub};
use crate::tables::M;

/// §3.5: the Q12 impulse response of `W(z)/Â(z)` for one subframe.
#[must_use]
pub fn impulse_response_fx(
    a_hat: &[i16; M + 1],
    ap1: &[i16; M + 1],
    ap2: &[i16; M + 1],
) -> [i16; L_SUBFR] {
    let mut input = [0i16; L_SUBFR];
    input[..=M].copy_from_slice(ap1);
    let zero = [0i16; M];
    let stage1 = syn_filt_mem(a_hat, &zero, &input);
    syn_filt_mem(ap2, &zero, &stage1)
}

/// The §3.6 / §3.10 filter memories.
#[derive(Debug, Clone, Default)]
pub struct TargetChainFx {
    /// `e(n−M) … e(n−1)` — the `1/Â(z)` state (most recent last).
    mem_err: [i16; M],
    /// `ew(n−M) … ew(n−1)` — the `1/A(z/γ₂)` state.
    mem_w: [i16; M],
}

impl TargetChainFx {
    /// Fresh chain with zero memories (clause 4.3).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// §3.6: the target `x(n)` for the subframe's residual `res`.
    /// Pure — the memories advance only through [`Self::update`].
    #[must_use]
    pub fn target(
        &self,
        res: &[i16; L_SUBFR],
        a_hat: &[i16; M + 1],
        ap1: &[i16; M + 1],
        ap2: &[i16; M + 1],
    ) -> [i16; L_SUBFR] {
        let e = syn_filt_mem(a_hat, &self.mem_err, res);
        let r = residu(ap1, &self.mem_err, &e);
        syn_filt_mem(ap2, &self.mem_w, &r)
    }

    /// §3.10: commits the memories after the subframe's excitation is
    /// fixed — `e(n) = s(n) − ŝ(n)` and eq (76) `ew(n)`.
    #[allow(clippy::too_many_arguments)]
    pub fn update(
        &mut self,
        s: &[i16; L_SUBFR],
        s_hat: &[i16; L_SUBFR],
        x: &[i16; L_SUBFR],
        y: &[i16; L_SUBFR],
        z_q12: &[i16; L_SUBFR],
        gain_pit_q14: i16,
        gain_code_q1: i16,
    ) {
        for i in 0..M {
            let n = L_SUBFR - M + i;
            self.mem_err[i] = sub(s[n], s_hat[n]);
            // ĝ_p·y on the doubled Q15 grid; ĝ_c·z: Q1 × Q12 doubled
            // = Q14, shifted once onto the same Q15 grid.
            let pit = l_mult(y[n], gain_pit_q14);
            let code = l_shl(l_mult(z_q12[n], gain_code_q1), 1);
            let contribution = round(l_shl(l_add(pit, code), 1)); // Q16 → Q0
            self.mem_w[i] = sub(x[n], contribution);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Trivial filters give the unit impulse on Q12.
    #[test]
    fn identity_impulse_response() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        let h = impulse_response_fx(&a, &a, &a);
        assert_eq!(h[0], 4096);
        assert!(h[1..].iter().all(|&v| v == 0));
    }

    /// From zero state the target is the residual convolved with `h`.
    #[test]
    fn target_is_convolution_from_zero_state() {
        let mut a_hat = [0i16; M + 1];
        a_hat[0] = 4096;
        a_hat[1] = -2458;
        a_hat[2] = 1229;
        let mut ap1 = a_hat;
        ap1[1] = -2300;
        ap1[2] = 1100;
        let mut ap2 = a_hat;
        ap2[1] = -1400;
        ap2[2] = 400;
        let h = impulse_response_fx(&a_hat, &ap1, &ap2);
        let chain = TargetChainFx::new();
        let res: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 7 % 23) as i16) * 30 - 330);
        let x = chain.target(&res, &a_hat, &ap1, &ap2);
        let conv = crate::fx::filters::convolve_h_q0(&res, &h);
        for n in 0..L_SUBFR {
            assert!(
                (i32::from(x[n]) - i32::from(conv[n])).abs() <= 2 + n as i32 / 8,
                "n={n}: {} vs {}",
                x[n],
                conv[n]
            );
        }
    }

    /// A perfect reconstruction leaves both memories at zero.
    #[test]
    fn perfect_reconstruction_keeps_zero_state() {
        let mut chain = TargetChainFx::new();
        let s: [i16; L_SUBFR] = std::array::from_fn(|n| (n as i16) * 10);
        let x = [0i16; L_SUBFR];
        let y = [0i16; L_SUBFR];
        let z = [0i16; L_SUBFR];
        chain.update(&s, &s, &x, &y, &z, 8192, 100);
        assert!(chain.mem_err.iter().all(|&v| v == 0));
        assert!(chain.mem_w.iter().all(|&v| v == 0));
    }
}
