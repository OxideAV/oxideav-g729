//! §3.4 **open-loop pitch analysis** on the fixed grid — the eq (34)
//! correlation maxima over the three delay sections, the eq (35)
//! normalisation, and the favour-lower-delays selection.
//!
//! This stage currently evaluates the search through the
//! spec-equation float module ([`crate::open_loop_pitch`]) on the Q0
//! weighted speech; moving the correlations onto the Word32 grid with
//! the overflow-rescale protocol is the §3.4 migration step proper.

use crate::open_loop_pitch::{open_loop_pitch, PIT_BUFFER as PIT_BUFFER_F};

/// Maximum open-loop delay.
pub const PIT_MAX: usize = 143;
/// Frame length.
pub const L_FRAME: usize = 80;
/// `[143 samples of history | 80 samples of the present frame]`.
pub const PIT_BUFFER: usize = PIT_MAX + L_FRAME;

const _: () = assert!(PIT_BUFFER == PIT_BUFFER_F);

/// The selected open-loop delay.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpenLoopPitchFx {
    /// `T_op ∈ [20, 143]`.
    pub t_op: i32,
}

/// §3.4 over one frame of weighted speech (Q0).
#[must_use]
pub fn open_loop_pitch_fx(sw: &[i16; PIT_BUFFER]) -> OpenLoopPitchFx {
    let buf: [f32; PIT_BUFFER] = std::array::from_fn(|n| f32::from(sw[n]));
    let r = open_loop_pitch(&buf);
    OpenLoopPitchFx {
        t_op: r.t_op as i32,
    }
}
