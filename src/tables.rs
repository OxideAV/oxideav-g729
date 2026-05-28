//! Bit-exact numeric tables from the ITU-T G.729 specification.
//!
//! The constants in this module are generated at build time from CSV
//! files under the crate's `tables/` directory, which are verbatim
//! copies of the spec-role-named CSV workspace at
//! `docs/audio/g729/tables/`. See that workspace's `README.md` for
//! the full provenance chain.
//!
//! The constants are unconditionally typed `[i16; N]` because every
//! exposed table in the current scope is a `Word16` array in the
//! G.729 fixed-point convention. Q-format is recorded in each
//! constant's doc-string as a suffix on the identifier name (for
//! example `_Q13` means each value is a 16-bit signed integer
//! representing a value in `Q13` fixed-point — i.e. the implicit
//! divisor is `2^13 = 8192`).
//!
//! ## Current scope
//!
//! Only a small foundation subset is exposed today:
//!
//! - Pre-/post-processing IIR high-pass filter coefficients
//!   (§3.1 / §4.2) for both fc = 100 Hz and fc = 140 Hz variants.
//! - The §4.1 Table 8 bit-allocation raw extract.
//! - The three basic-op fixed-point math LUTs used by `Pow2`,
//!   `Log2`, and `Inv_sqrt`.
//!
//! Additional tables (LSP codebooks, gain quantizer codebooks,
//! interpolation filters, autocorrelation lag window, search grid,
//! taming zone table, Annex B DTX/CNG/VAD tables) remain available
//! in `docs/audio/g729/tables/` and will be wired into this module
//! in a future round once they are needed by a behavioural module.

include!(concat!(env!("OUT_DIR"), "/tables_data.rs"));
