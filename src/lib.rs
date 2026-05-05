//! ITU-T G.729 (CS-ACELP, 8 kbit/s) decoder — first real implementation.
//!
//! Pipeline (all pure-Rust):
//! - `bitreader`: 80-bit frame → 15 bit fields (L0..L3, P1/P0, C1/S1,
//!   GA1/GB1, P2, C2/S2, GA2/GB2), per §3.6 / Table 8.
//! - `lpc`: LSP decode from index quadruple via MA-4 predictor + safety
//!   monotonicity (§3.2.4), LSP ↔ LPC (§3.2.6), LSP interpolation
//!   between the two subframes (§3.2.5).
//! - `lsp_tables`: static codebook tables transcribed verbatim from the
//!   ITU-T G.729 reference C source `TAB_LD8K.C` (`LSPCB1_Q13`,
//!   `LSPCB2_Q13`, `FG_Q15`, `FG_SUM_Q15`, `FG_SUM_INV_Q12`). LSP
//!   quantisation is bit-exact against the reference. The gain-VQ
//!   tables (`GBK1` / `GBK2`) remain reduced — see the per-crate
//!   README for the full list of remaining gaps.
//! - `synthesis`: adaptive codebook (fractional-pitch, 1/3 resolution),
//!   algebraic fixed codebook (4-track × 4-pulse), two-stage gain VQ,
//!   10th-order all-pole synthesis filter, short-term + long-term
//!   postfilter with tilt compensation and AGC.
//! - `decoder`: orchestrates the pipeline per packet and exposes it via
//!   the `oxideav_core::Decoder` trait.
//!
//! Reference: ITU-T Recommendation G.729 (January 2007 edition) +
//! Annex A (`G.729a` simplified-complexity variant).

#![allow(
    clippy::needless_range_loop,
    clippy::unnecessary_cast,
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::doc_lazy_continuation,
    clippy::doc_overindented_list_items
)]

pub mod annex_b_cng;
pub mod annex_b_vad;
pub mod bitreader;
pub mod codec;
pub mod decoder;
pub mod encoder;
pub mod lpc;
pub mod lsp_tables;
pub mod synthesis;

use oxideav_core::CodecRegistry;

pub const CODEC_ID_STR: &str = "g729";

/// Extradata flag byte that enables Annex B (VAD/DTX/CNG) on the
/// encoder. Place this byte in [`oxideav_core::CodecParameters::extradata`]
/// to opt in to silence compression. Default is off.
///
/// Example:
/// ```ignore
/// use oxideav_core::{CodecId, CodecParameters};
/// use oxideav_g729::{CODEC_ID_STR, ANNEX_B_ENABLE_EXTRADATA};
/// let mut p = CodecParameters::audio(CodecId::new(CODEC_ID_STR));
/// p.extradata = vec![ANNEX_B_ENABLE_EXTRADATA];
/// ```
pub const ANNEX_B_ENABLE_EXTRADATA: u8 = 0xB1;

/// Number of samples per G.729 frame (10 ms @ 8 kHz).
pub const FRAME_SAMPLES: usize = 80;

/// Number of samples per 5-ms subframe.
pub const SUBFRAME_SAMPLES: usize = 40;

/// Number of subframes per frame.
pub const SUBFRAMES_PER_FRAME: usize = 2;

/// LPC order (10th-order short-term predictor).
pub const LPC_ORDER: usize = 10;

/// Encoded frame size in bytes (80 bits = 10 bytes).
pub const FRAME_BYTES: usize = 10;

/// Sample rate (Hz).
pub const SAMPLE_RATE: u32 = 8_000;

/// Register G.729 with the supplied [`CodecRegistry`]. Prefer the
/// unified [`register`] entry point when you have a
/// [`oxideav_core::RuntimeContext`] in hand.
pub fn register_codecs(reg: &mut CodecRegistry) {
    codec::register_codecs(reg);
}

/// Unified registration entry point — installs G.729 into the codec
/// sub-registry of the supplied [`oxideav_core::RuntimeContext`].
pub fn register(ctx: &mut oxideav_core::RuntimeContext) {
    register_codecs(&mut ctx.codecs);
}

oxideav_core::register!("g729", register);

#[cfg(test)]
mod register_tests {
    use super::*;
    use oxideav_core::CodecId;

    #[test]
    fn register_via_runtime_context_installs_codec_factory() {
        let mut ctx = oxideav_core::RuntimeContext::new();
        register(&mut ctx);
        let id = CodecId::new(CODEC_ID_STR);
        assert!(
            ctx.codecs.has_decoder(&id),
            "decoder factory not installed via RuntimeContext"
        );
        assert!(
            ctx.codecs.has_encoder(&id),
            "encoder factory not installed via RuntimeContext"
        );
    }
}
