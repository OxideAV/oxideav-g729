//! Bit-exact numeric tables required by the ITU-T G.729 (06/2012) main
//! 8 kbit/s coder. Each `pub const` here is generated at build time by
//! `build.rs` from `tables/<spec-name>.csv`, which is itself a staged
//! copy of `docs/audio/g729/tables/<spec-name>.csv` (extractor output,
//! see that directory's `README.md` for the provenance chain).
//!
//! Only data is reproduced. No algorithmic source distributed
//! alongside the ITU recommendation is read here or anywhere else in
//! this crate (`docs/IMPLEMENTOR_ROUND.md`, clean-room rule).
//!
//! ## Coverage
//!
//! Round 173 (foundation):
//!
//! * §3.1 / §4.2 pre/post high-pass filters (b100/a100 Q13;
//!   b140/a140 Q12).
//! * §4.1 spec Table 8 bit allocation per parameter (`bitsno`).
//! * basic_op() math LUTs (`Pow2`, `Log2`, `Inv_sqrt`).
//!
//! Round 189 (LP analysis / pitch / gain prediction):
//!
//! * §3.2.1 LP analysis windowing — `hamwindow` 240-sample Q15
//!   Hamming window, `lag_h` / `lag_l` 10-entry Q15 60 Hz
//!   bandwidth-expansion lag-window pair.
//! * §3.2.5 `az_lsf()` cosine grid (61 Q15 entries spanning the
//!   half-circle in 60 equal steps).
//! * §3.7 pitch interpolation filters — `inter_3` 13-tap analysis,
//!   `inter_3l` 31-tap synthesis (both Q15).
//! * §3.9 MA gain-prediction coefficients `pred` — 4 Q13 entries
//!   {0.68, 0.58, 0.34, 0.19}.
//!
//! Round 195 (LSP quantiser two-stage VQ codebooks):
//!
//! * §3.2.4 first-stage codebook (`lspcb1`) — [`LSP_QUANT_CODEBOOK_L1_Q13`],
//!   shape `[[i16; M]; NC0]` = `[[i16; 10]; 128]`, Q13 (7-bit index).
//! * §3.2.4 second-stage codebook (`lspcb2`) — [`LSP_QUANT_CODEBOOK_L2_Q13`],
//!   shape `[[i16; M]; NC1]` = `[[i16; 10]; 32]`, Q13. Per the staged
//!   trace doc §3.5 the spec describes two 5-D codebooks L2 + L3, but
//!   the staged 32 × 10 single-array packing exposes the lower 5
//!   coefficients (indices `0..M/2`) as the L2 contribution and the
//!   upper 5 (`M/2..M`) as L3, both 5-bit indexed.
//! * [`lsp_l1_entry`] / [`lsp_l2_entry`] / [`lsp_l3_entry`] — bounds-
//!   checked lookup helpers returning a borrowed 5- or 10-slice into
//!   the compiled codebook.
//!
//! Round 201 (LSP MA-predictor `fg` family — completes §3.2.4
//! reconstruction inputs per spec eqs (20) / (20a)):
//!
//! * [`LSP_MA_PREDICTOR_FG_Q15`] — shape `[[[i16; M]; MA_NP]; 2]` =
//!   `[[[i16; 10]; 4]; 2]`, Q15. Outer dim selects the L0 predictor
//!   mode (1 bit); the inner `[MA_NP][M]` plane carries the 4th-order
//!   MA coefficients across each LSP coordinate.
//! * [`LSP_MA_PREDICTOR_FG_SUM_Q15`] — shape `[[i16; M]; 2]`, Q15.
//!   Per-mode column sums of the `fg` plane (`Σ_k fg[mode][k][i]` for
//!   each `i ∈ 0..M`); the reconstruction equation factor used to
//!   normalise the MA prediction.
//! * [`LSP_MA_PREDICTOR_FG_SUM_INV_Q12`] — shape `[[i16; M]; 2]`, Q12.
//!   Pre-tabulated reciprocal of the per-mode column sums.
//! * [`lsp_fg_plane`] / [`lsp_fg_sum`] / [`lsp_fg_sum_inv`] —
//!   bounds-checked accessor helpers returning a borrowed
//!   per-mode plane / row.
//! * New spec-dimension constant `MA_NP = 4` (LSP MA prediction
//!   order); `L0_BITS == 1` matches the outer-dim count exactly,
//!   `1 << L0_BITS == 2`.
//!
//! Still NOT compiled (gated on the docs collaborator specifier
//! pass): gain GA/GB codebooks, postfilter interpolation
//! (`tab_hup_*`), taming (`tab_zone`), Annex B DTX/CNG,
//! LSF↔LSP cos/slope tables.
//!
//! ## Q-format convention reminder (G.729 §1.4)
//!
//! G.729's fixed-point arithmetic uses Q-format scaling where
//! `Q<n>` denotes that a Word16 value `v` represents the rational
//! number `v / 2^n`. The Q exponent is embedded in the constant
//! identifier (e.g. `Q13`, `Q15`) so callers do not have to look up
//! each table's scaling at the use site.

include!(concat!(env!("OUT_DIR"), "/preproc-highpass-100Hz-b-Q13.rs"));
include!(concat!(env!("OUT_DIR"), "/preproc-highpass-100Hz-a-Q13.rs"));
include!(concat!(env!("OUT_DIR"), "/preproc-highpass-140Hz-b-Q12.rs"));
include!(concat!(env!("OUT_DIR"), "/preproc-highpass-140Hz-a-Q12.rs"));
include!(concat!(
    env!("OUT_DIR"),
    "/bit-allocation-per-parameter-table8.rs"
));
include!(concat!(env!("OUT_DIR"), "/basic-op-pow2-table.rs"));
include!(concat!(env!("OUT_DIR"), "/basic-op-log2-table.rs"));
include!(concat!(env!("OUT_DIR"), "/basic-op-invsqrt-table.rs"));
include!(concat!(env!("OUT_DIR"), "/lpc-hamming-window-Q15.rs"));
include!(concat!(
    env!("OUT_DIR"),
    "/lpc-autocorr-lag-window-high-Q15.rs"
));
include!(concat!(
    env!("OUT_DIR"),
    "/lpc-autocorr-lag-window-low-Q15.rs"
));
include!(concat!(env!("OUT_DIR"), "/lsf-search-grid-cos-Q15.rs"));
include!(concat!(
    env!("OUT_DIR"),
    "/pitch-interpolation-filter-analysis-Q15.rs"
));
include!(concat!(
    env!("OUT_DIR"),
    "/pitch-interpolation-filter-synthesis-Q15.rs"
));
include!(concat!(
    env!("OUT_DIR"),
    "/gain-quantizer-ma-predictor-Q13.rs"
));
include!(concat!(
    env!("OUT_DIR"),
    "/lsp-quantizer-codebook-L1-Q13.rs"
));
include!(concat!(
    env!("OUT_DIR"),
    "/lsp-quantizer-codebook-L2-Q13.rs"
));
include!(concat!(env!("OUT_DIR"), "/lsp-ma-predictor-fg-Q15.rs"));
include!(concat!(env!("OUT_DIR"), "/lsp-ma-predictor-fg-sum-Q15.rs"));
include!(concat!(
    env!("OUT_DIR"),
    "/lsp-ma-predictor-fg-sum-inv-Q12.rs"
));

/// G.729 §4.1 transmitted parameter count per frame (spec `PRM_SIZE`).
/// The bit-allocation table [`BIT_ALLOCATION_TABLE8`] carries 13
/// entries (the literal count of the staged-CSV `bitsno` array); only
/// the first `PRM_SIZE` indices correspond to the 11 parameters
/// defined in spec Table 8. The remaining trailing entries are
/// preserved as-extracted to match the source array byte-for-byte and
/// to avoid silently dropping data — consumers indexing the table for
/// spec-defined parameters bound their loop by `PRM_SIZE`.
pub const PRM_SIZE: usize = 11;

/// G.729 8 kbit/s coder headline rate: 80 bits per 10 ms frame, per
/// spec §1.2 "Brief description of CS-ACELP". This is the
/// specification-stated frame size; the relationship between this
/// value and the 13-entry [`BIT_ALLOCATION_TABLE8`] extraction is
/// documented in the table-shape test (the literal CSV sum is 66,
/// not 80 — spec Table 8 itself lists 15 per-parameter rows, whereas
/// the staged `bitsno` array packs 13 values; reconciling the two is
/// deferred to the docs collaborator).
pub const BITS_PER_FRAME: usize = 80;

/// G.729 LP-analysis predictor order — 10, per spec §3.2.1. The
/// 10th-order LP filter is fit to each `L_WINDOW = 240`-sample
/// analysis frame using the autocorrelation method.
pub const M: usize = 10;

/// G.729 LP-analysis window length — 240 samples (`L_WINDOW`), per
/// spec §3.2.1. The Hamming window in [`LPC_HAMMING_WINDOW_Q15`]
/// spans this number of samples.
pub const L_WINDOW: usize = 240;

/// G.729 LSF root-search grid resolution — 60 evenly-spaced points
/// on the upper half-circle, per spec §3.2.5 (`GRID_POINTS`). The
/// staged cosine grid [`LSF_SEARCH_GRID_COS_Q15`] holds
/// `GRID_POINTS + 1 = 61` samples.
pub const GRID_POINTS: usize = 60;

/// G.729 §3.2.4 first-stage LSP VQ codebook size — `NC0 = 128`
/// entries indexed by the 7-bit `L1` parameter. The full codebook
/// [`LSP_QUANT_CODEBOOK_L1_Q13`] has shape `[[i16; M]; NC0]`.
pub const NC0: usize = 128;

/// G.729 §3.2.4 second-stage LSP VQ codebook size — `NC1 = 32`
/// entries per split (5-bit indices `L2` and `L3`). The packed
/// second-stage table [`LSP_QUANT_CODEBOOK_L2_Q13`] has shape
/// `[[i16; M]; NC1]`; the lower five coefficients of each row are
/// the L2 contribution, the upper five are L3.
pub const NC1: usize = 32;

/// G.729 §3.2.4 LSP MA-prediction order — 4 taps, per spec eq (20).
/// The `fg` matrix carries `MA_NP` per-tap coefficient rows for each
/// of the two predictor modes selected by `L0`.
pub const MA_NP: usize = 4;

/// Bit-width of the first-stage LSP-VQ index `L1` per spec Table 1
/// / §3.2.4 — 7 bits, so `1 << L1_BITS == NC0 == 128`.
pub const L1_BITS: usize = 7;

/// Bit-width of the second-stage L2 / L3 LSP-VQ indices per spec
/// Table 1 / §3.2.4 — 5 bits each, so `1 << L2_BITS == NC1 == 32`.
pub const L2_BITS: usize = 5;

/// Bit-width of `L3` (same split-VQ resolution as L2).
pub const L3_BITS: usize = 5;

/// Bit-width of the predictor-switch flag `L0` per spec Table 1 /
/// §3.2.4 — 1 bit selecting one of two MA-predictor histories. The
/// per-mode MA coefficient plane is [`LSP_MA_PREDICTOR_FG_Q15`] (the
/// outer dimension is indexed by `L0`); the helpers
/// [`lsp_fg_plane`] / [`lsp_fg_sum`] / [`lsp_fg_sum_inv`] borrow the
/// per-mode rows.
pub const L0_BITS: usize = 1;

/// Total transmitted LSP-quantiser bit count per frame:
/// `L0 + L1 + L2 + L3 = 1 + 7 + 5 + 5 = 18`, matching the spec
/// Table 1 entry for "Line spectrum pairs".
pub const LSP_TOTAL_BITS: usize = L0_BITS + L1_BITS + L2_BITS + L3_BITS;

/// Returns a borrowed slice of the §3.2.4 first-stage `lspcb1`
/// codebook row at index `l1`, which holds all `M = 10` Q13
/// coefficients of the first-stage contribution for that 7-bit
/// codeword.
///
/// # Panics
///
/// Panics if `l1 >= NC0` (i.e. the caller passed a value not
/// expressible in the 7-bit `L1` field). Callers that decode `L1`
/// from a transmitted frame can rely on the fact that masking with
/// `(1 << L1_BITS) - 1` always yields an in-range value.
#[must_use]
pub fn lsp_l1_entry(l1: usize) -> &'static [i16; M] {
    &LSP_QUANT_CODEBOOK_L1_Q13[l1]
}

/// Returns the lower-5 (L2) half of the packed `lspcb2` row at
/// index `l2`. Each value is the Q13 second-stage contribution to
/// LSP coefficients `0..M/2`.
///
/// # Panics
///
/// Panics if `l2 >= NC1`.
#[must_use]
pub fn lsp_l2_entry(l2: usize) -> &'static [i16] {
    &LSP_QUANT_CODEBOOK_L2_Q13[l2][..M / 2]
}

/// Returns the upper-5 (L3) half of the packed `lspcb2` row at
/// index `l3`. Each value is the Q13 second-stage contribution to
/// LSP coefficients `M/2..M`.
///
/// # Panics
///
/// Panics if `l3 >= NC1`.
#[must_use]
pub fn lsp_l3_entry(l3: usize) -> &'static [i16] {
    &LSP_QUANT_CODEBOOK_L2_Q13[l3][M / 2..]
}

/// Returns the §3.2.4 `fg` MA-predictor coefficient plane for the
/// `L0`-selected predictor mode — a `[MA_NP][M]` Q15 slab whose
/// `k`-th row holds the order-`k` MA tap across the 10 LSP
/// coordinates.
///
/// The 1-bit `L0` parameter selects between the two predictor modes
/// (eq (20)); the encoder picks whichever mode gives the lower
/// reconstruction error.
///
/// # Panics
///
/// Panics if `mode >= 2` (i.e. the caller passed a value outside the
/// 1-bit `L0` domain). Callers decoding `L0` from a transmitted frame
/// can rely on the fact that masking with `(1 << L0_BITS) - 1` always
/// yields an in-range value.
#[must_use]
pub fn lsp_fg_plane(mode: usize) -> &'static [[i16; M]; MA_NP] {
    &LSP_MA_PREDICTOR_FG_Q15[mode]
}

/// Returns the §3.2.4 `fg_sum` row for the `L0`-selected predictor
/// mode — the Q15 per-coordinate sum
/// `Σ_{k=0..MA_NP} fg[mode][k][i]` for `i ∈ 0..M`.
///
/// This per-mode sum participates in the reconstruction equation
/// (spec eq (20a)) that maps the L1/L2/L3 codebook contributions and
/// the past quantised LSP residuals onto the final reconstructed LSP
/// vector.
///
/// # Panics
///
/// Panics if `mode >= 2`.
#[must_use]
pub fn lsp_fg_sum(mode: usize) -> &'static [i16; M] {
    &LSP_MA_PREDICTOR_FG_SUM_Q15[mode]
}

/// Returns the §3.2.4 `fg_sum_inv` row for the `L0`-selected
/// predictor mode — the Q12 per-coordinate reciprocal of
/// [`lsp_fg_sum`], pre-tabulated so the reconstruction can avoid a
/// per-sample division.
///
/// # Panics
///
/// Panics if `mode >= 2`.
#[must_use]
pub fn lsp_fg_sum_inv(mode: usize) -> &'static [i16; M] {
    &LSP_MA_PREDICTOR_FG_SUM_INV_Q12[mode]
}
