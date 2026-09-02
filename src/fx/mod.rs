//! 16-bit fixed-point arithmetic layer for the bit-exact decode drive.
//!
//! Clause 2 of the Recommendation states that "the description of the
//! speech coding algorithm … is made in terms of bit-exact fixed-point
//! mathematical operations", and clause 5 defines the numeric model
//! that description runs on:
//!
//! - **Table 10** — exactly two fixed-point data types: `Word16`
//!   (signed two's-complement 16-bit, `0x8000..=0x7fff`) and `Word32`
//!   (signed two's-complement 32-bit, `0x80000000..=0x7fffffffL`).
//! - **Table 11** — the closed set of *basic operators* every
//!   computation is performed with (`add`, `sub`, `mult`, `L_mac`,
//!   `round`, `norm_s`, `div_s`, …).
//! - **Table 15** — the general routine groups, including the extended
//!   32-bit operator family and the mathematical functions (base-2
//!   logarithm, `2^x`, inverse square root) whose lookup tables are
//!   staged as data (`tables/basic-op-{log2,pow2,invsqrt}-table.csv`).
//!
//! [`ops`] implements the Table 11 operator set with the conventional
//! saturating two's-complement semantics those operator names denote;
//! [`dsp`] builds the Table 15 mathematical functions on top of the
//! staged lookup tables. Everything downstream of these two modules is
//! ordinary spec-clause arithmetic expressed through the operator set,
//! which is what makes bit-exactness against the ITU conformance
//! vectors reachable at all — the float decode chain tracks the same
//! equations but its rounding drift compounds through the
//! adaptive-codebook feedback until long vectors decorrelate.

pub mod acelp;
pub mod analysis;
pub mod decoder;
pub mod dsp;
pub mod encoder;
pub mod excitation;
pub mod filters;
pub mod gain_vq;
pub mod gains;
pub mod levinson;
pub mod lp_to_lsp;
pub mod lsp;
pub mod ops;
pub mod pitch_cl;
pub mod pitch_ol;
pub mod postfilter;
pub mod target;
pub mod weighting;
