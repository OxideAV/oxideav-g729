//! # oxideav-g729
//!
//! **Status:** orphan-rebuild scaffold (reset 2026-05-24).
//!
//! The prior implementation was retired under the workspace clean-room
//! policy: several source modules transcribed numeric tables verbatim
//! from an external reference-software distribution and described
//! matching its behaviour by citing specific source files of that
//! distribution. The clean-room policy forbids consulting any external
//! implementation's source for any reason, regardless of licensing, so
//! the provenance of those tables could not be defended. The crate will
//! be re-implemented from scratch against the staged ITU-T G.729
//! Recommendation text in a future clean-room round.
//!
//! Every public API currently returns [`Error::NotImplemented`].

#![warn(missing_debug_implementations)]

use oxideav_core::RuntimeContext;

/// Crate-local error type. Until the clean-room rebuild lands every
/// public API path returns [`Error::NotImplemented`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// The crate has been reset to a scaffold pending clean-room
    /// rebuild; no decoder or encoder functionality is wired up yet.
    NotImplemented,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "oxideav-g729: orphan-rebuild scaffold — no codec wired up"
        )
    }
}

impl std::error::Error for Error {}

/// No-op codec registration — the orphan-rebuild scaffold registers
/// nothing into the runtime context.
pub fn register(_ctx: &mut RuntimeContext) {}

oxideav_core::register!("g729", register);
