//! §3.8.1 **algebraic codebook search** on the fixed grid.
//!
//! This stage currently evaluates the eq (51)–(60) focused search
//! through the spec-equation float module
//! ([`crate::fixed_codebook_search`]) on the Q0 target and the Q12
//! pre-filtered impulse response; moving `d(n)`, `φ(i, j)` and the
//! `C²/E` comparison onto the Word16/Word32 grid is the §3.8
//! migration step proper.

use crate::fixed_codebook_search::{correlation_d, phi_matrix, search_fixed_codebook};
use crate::fx::filters::L_SUBFR;

/// Maximum number of fourth-loop entries per frame (clause 3.8.1).
pub const MAX_LOOP4_PER_FRAME: u32 = 180;

/// One selected codevector: Table-7 pulse positions and signs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcelpChoice {
    /// Pulse positions `m₀ … m₃`.
    pub positions: [u8; 4],
    /// Pulse signs `s₀ … s₃` (`±1`).
    pub signs: [i8; 4],
}

/// §3.8.1 focused search over one subframe. `budget` is the frame's
/// remaining fourth-loop budget (consumed in place).
#[must_use]
pub fn search_acelp_fx(
    x_prime: &[i16; L_SUBFR],
    h_pre_q12: &[i16; L_SUBFR],
    budget: &mut u32,
) -> AcelpChoice {
    let xf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(x_prime[n]));
    let hf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(h_pre_q12[n]) / 4096.0);
    let d = correlation_d(&xf, &hf);
    let phi = phi_matrix(&hf);
    let c = search_fixed_codebook(&d, &phi, budget);
    AcelpChoice {
        positions: c.positions,
        signs: c.signs,
    }
}
