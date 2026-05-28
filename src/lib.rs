//! # oxideav-g729
//!
//! **Status:** orphan-rebuild scaffold (reset 2026-05-24); early
//! data-only foundation landed 2026-05-29 (round 173).
//!
//! The prior implementation was retired under the workspace
//! clean-room policy: several source modules transcribed numeric
//! tables verbatim from an external reference-software distribution
//! and described matching its behaviour by citing specific source
//! files of that distribution. The clean-room policy forbids
//! consulting any external implementation's source for any reason,
//! regardless of licensing.
//!
//! ## What lives here today
//!
//! * A scaffold `register()` function that wires no codec into the
//!   runtime context yet.
//! * A [`tables`] module exposing a small subset of bit-exact ITU-T
//!   G.729 numeric tables as `pub const [i16; N]` arrays. The tables
//!   are emitted at build time from the spec-role-named CSV
//!   workspace at `docs/audio/g729/tables/`, carried in this crate
//!   under `tables/` for hermetic publishing. No algorithmic source
//!   is read by the extractor or by `build.rs`; the CSVs themselves
//!   are pure factual numeric content (see
//!   `docs/audio/g729/tables/README.md` for the provenance chain).
//!
//! The remaining decoder and encoder behaviour is still pending the
//! clean-room trace doc the docs collaborator is preparing.
//!
//! Every public codec API path currently returns
//! [`Error::NotImplemented`].

#![warn(missing_debug_implementations)]

use oxideav_core::RuntimeContext;

pub mod tables;

/// Crate-local error type. Until the clean-room rebuild lands every
/// public codec API path returns [`Error::NotImplemented`].
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
