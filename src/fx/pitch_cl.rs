//! §3.7 **closed-loop pitch search** and the §3.7.3 adaptive-codebook
//! gain on the fixed grid.
//!
//! This stage currently evaluates the eq (37)–(39) search through the
//! spec-equation float module ([`crate::closed_loop_pitch`]) on the
//! Q0 target / Q12 impulse response / Q0 excitation, and lands the
//! eq (43) gain on Q14; moving the correlations and the `b12`
//! interpolation onto the Word32 grid is the §3.7 migration step
//! proper.

use crate::closed_loop_pitch::{closed_loop_search, EXC_BUFFER};
use crate::fx::excitation::EXC_BUF;
use crate::fx::filters::L_SUBFR;
use crate::pitch_decode::PitchDelay;

/// Minimum pitch delay.
pub const PIT_MIN: i32 = 20;
/// Maximum pitch delay.
pub const PIT_MAX: i32 = 143;
/// Integer-only threshold for `T1` (clause 3.7).
pub const T1_FRACTIONAL_LIMIT: i32 = 85;

/// Unstated latitude of the closed-loop search.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ClosedLoopLatitude {}

/// One subframe's search result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClosedLoopResultFx {
    /// The selected fractional delay.
    pub delay: PitchDelay,
}

/// The §3.7 subframe-1 window around `T_op`.
#[must_use]
pub fn subframe1_window(t_op: i32) -> (i32, i32) {
    let mut t_min = t_op - 3;
    if t_min < PIT_MIN {
        t_min = PIT_MIN;
    }
    let mut t_max = t_min + 6;
    if t_max > PIT_MAX {
        t_max = PIT_MAX;
        t_min = t_max - 6;
    }
    (t_min, t_max)
}

/// The §3.7 subframe-2 window around `int(T1)`.
#[must_use]
pub fn subframe2_window(int_t1: i32) -> (i32, i32) {
    let mut t_min = int_t1 - 5;
    if t_min < PIT_MIN {
        t_min = PIT_MIN;
    }
    let mut t_max = t_min + 9;
    if t_max > PIT_MAX {
        t_max = PIT_MAX;
        t_min = t_max - 9;
    }
    (t_min, t_max)
}

/// §3.7 over one subframe. `exc` is the decoder-layout excitation
/// buffer with the current subframe at `off` already holding the LP
/// residual (clause 3.7); `h` the Q12 impulse response.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn closed_loop_search_fx(
    x: &[i16; L_SUBFR],
    h_q12: &[i16; L_SUBFR],
    exc: &[i16; EXC_BUF],
    off: usize,
    t_min: i32,
    t_max: i32,
    frac_limit: i32,
    _lat: &ClosedLoopLatitude,
) -> ClosedLoopResultFx {
    let xf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(x[n]));
    let hf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(h_q12[n]) / 4096.0);
    // u(−143) … u(39) → index PIT_MAX + n.
    let ef: [f32; EXC_BUFFER] = std::array::from_fn(|i| f32::from(exc[off - 143 + i]));
    let r = closed_loop_search(&xf, &hf, &ef, t_min, t_max, frac_limit);
    ClosedLoopResultFx { delay: r.delay }
}

/// eq (43): the adaptive-codebook gain `g_p = xᵗy / yᵗy` bounded to
/// `[0, 1.2]`, on Q14.
#[must_use]
pub fn pitch_gain_q14(x: &[i16; L_SUBFR], y: &[i16; L_SUBFR]) -> i16 {
    let xf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(x[n]));
    let yf: [f32; L_SUBFR] = std::array::from_fn(|n| f32::from(y[n]));
    let g = crate::fixed_codebook_search::adaptive_gain(&xf, &yf);
    (f64::from(g) * 16384.0).round().clamp(0.0, 19661.0) as i16
}
