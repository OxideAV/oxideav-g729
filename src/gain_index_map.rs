//! §3.9.3 gain-quantiser codeword index mapping (robustness layer).
//!
//! The §3.9.2 conjugate-structure gain VQ chooses two **codebook
//! indices** — `GA` in `0..NCODE1` (= `0..8`) and `GB` in `0..NCODE2`
//! (= `0..16`). Per spec §3.9.3 (clause heading "Codeword computation
//! for gain quantizer"), these codebook indices are **mapped** through
//! per-stage permutations before being placed on the wire:
//!
//! > "The codewords GA and GB for the gain quantizer are obtained from
//! > the indices corresponding to the best choice. To reduce the
//! > impact of single bit errors the codebook indices are mapped."
//! > — G.729 (06/2012), clause 3.9.3
//!
//! The spec prose stops there; the actual permutation is given by the
//! `map1` / `map2` (encoder-side, codebook→transmitted) and `imap1` /
//! `imap2` (decoder-side, transmitted→codebook) tables staged under
//! [`crate::tables`] as [`GAIN_QUANT_GA_PERMUTATION`] /
//! [`GAIN_QUANT_GB_PERMUTATION`] and
//! [`GAIN_QUANT_GA_INVERSE_PERMUTATION`] /
//! [`GAIN_QUANT_GB_INVERSE_PERMUTATION`]. Their mutual-inverse
//! property is pinned by the `tables_shape` integration test.
//!
//! ## What this module does
//!
//! - [`demap_ga`] and [`demap_gb`] are the decoder-side primitive:
//!   they take a **transmitted** index from the unpacked bit stream
//!   and return the corresponding **codebook** index that the
//!   §3.9.2 reconstruction (see [`crate::gain_reconstruct`])
//!   indexes the GA / GB codebook with. This is the step the spec
//!   §4.1.5 procedure implicitly invokes between
//!   [`crate::parameters::unpack_parameters`] and
//!   [`crate::gain_reconstruct::reconstruct_gains`].
//! - [`map_ga`] and [`map_gb`] are the symmetric encoder-side
//!   primitive: given the chosen codebook index, return the
//!   transmitted-domain value that goes on the wire. Useful for
//!   pack tests, future encoder wire-up, and the round-trip
//!   property exercised in the unit tests.
//! - [`demap_frame`] is the convenience wrapper that takes the
//!   typed [`Parameters`] from
//!   [`crate::parameters::unpack_parameters`] and yields a
//!   typed [`DemappedGainIndices`] holding the four codebook
//!   indices (subframe-1 GA / GB, subframe-2 GA / GB) ready to
//!   feed [`crate::gain_reconstruct::reconstruct_gains`].
//!
//! ## Why this matters for the decode path
//!
//! Without this module the round-225 unpacker's `ga1` / `gb1` / `ga2`
//! / `gb2` fields — which sit in the **transmitted** domain — were
//! being fed directly to the round-231 codebook lookup. Under the
//! §3.9.3 mapping that is bit-incorrect: the on-wire index `1` selects
//! codebook row `map1[1] = 1` (it happens to coincide by accident),
//! but on-wire index `0` actually selects codebook row `map1[0] = 5`,
//! not row `0`. The mismatch is invisible to the round-231 envelope
//! tests (every (GA, GB) pair stays inside the `ĝ_p ∈ [0, 2]` /
//! `γ̂ ∈ [0, 11]` window regardless of which permutation is applied),
//! so the corruption surfaces only when the reconstructed gains are
//! cross-checked against a known-good decoder output. Wiring the
//! §3.9.3 inverse permutation now is the smallest correct step that
//! puts the decode chain back in spec-conformance for the gain stage.
//!
//! ## Error surface
//!
//! All four primitives return [`GainIndexMapError`] on out-of-range
//! input. The transmitted GA / GB codewords are only 3 / 4 bits per
//! spec Table 8 so the decode-path entry points cannot produce these
//! errors in a well-formed frame; the typed error variant is
//! preserved so callers handing in fabricated `Parameters` structs
//! get the typed surface rather than a panic.

use crate::parameters::Parameters;
use crate::tables::{
    GAIN_QUANT_GA_INVERSE_PERMUTATION, GAIN_QUANT_GA_PERMUTATION,
    GAIN_QUANT_GB_INVERSE_PERMUTATION, GAIN_QUANT_GB_PERMUTATION, NCODE1, NCODE2,
};

/// Result of demapping the four transmitted gain-VQ indices in one
/// frame back into the §3.9.2 codebook-index domain.
///
/// Field order matches spec §4.1.5: subframe 1 first
/// ([`Self::ga1`], [`Self::gb1`]), then subframe 2
/// ([`Self::ga2`], [`Self::gb2`]). Every field lies in the
/// per-stage codebook domain (`0..NCODE1` for GA and `0..NCODE2`
/// for GB), so the values are safe to hand directly to
/// [`crate::gain_reconstruct::reconstruct_gains`] without any
/// further bounds check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DemappedGainIndices {
    /// Subframe-1 first-stage gain-VQ codebook index, after §3.9.3
    /// inverse mapping. Lies in `0..NCODE1`.
    pub ga1: usize,
    /// Subframe-1 second-stage gain-VQ codebook index, after §3.9.3
    /// inverse mapping. Lies in `0..NCODE2`.
    pub gb1: usize,
    /// Subframe-2 first-stage gain-VQ codebook index, after §3.9.3
    /// inverse mapping. Lies in `0..NCODE1`.
    pub ga2: usize,
    /// Subframe-2 second-stage gain-VQ codebook index, after §3.9.3
    /// inverse mapping. Lies in `0..NCODE2`.
    pub gb2: usize,
}

/// Errors returned by the §3.9.3 index-mapping helpers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GainIndexMapError {
    /// A GA-stage value (transmitted or codebook) did not fit in the
    /// [`NCODE1`]-entry permutation. The transmitted GA codeword is
    /// only 3 bits per spec Table 8 so a well-formed frame cannot
    /// produce this error; it surfaces when callers pass an
    /// out-of-domain index to one of [`map_ga`] / [`demap_ga`]
    /// directly.
    GaOutOfRange {
        /// The offending input index.
        index: usize,
    },
    /// A GB-stage value (transmitted or codebook) did not fit in the
    /// [`NCODE2`]-entry permutation. As with the GA variant the
    /// 4-bit transmitted GB codeword cannot trigger this from a
    /// well-formed frame.
    GbOutOfRange {
        /// The offending input index.
        index: usize,
    },
}

impl core::fmt::Display for GainIndexMapError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::GaOutOfRange { index } => write!(
                f,
                "g729 §3.9.3 gain-index map: GA index {index} >= NCODE1 = {NCODE1}",
            ),
            Self::GbOutOfRange { index } => write!(
                f,
                "g729 §3.9.3 gain-index map: GB index {index} >= NCODE2 = {NCODE2}",
            ),
        }
    }
}

impl std::error::Error for GainIndexMapError {}

/// Encoder-side §3.9.3 forward permutation for the first-stage gain
/// codebook: given the **codebook** index in `0..NCODE1`, return the
/// **transmitted** index that goes on the wire.
///
/// Implemented by indexing the [`GAIN_QUANT_GA_PERMUTATION`] table at
/// `codebook_index`. The companion [`demap_ga`] is the decode-side
/// inverse.
///
/// # Errors
///
/// Returns [`GainIndexMapError::GaOutOfRange`] if `codebook_index`
/// does not fit in `0..NCODE1`.
pub fn map_ga(codebook_index: usize) -> Result<usize, GainIndexMapError> {
    if codebook_index >= NCODE1 {
        return Err(GainIndexMapError::GaOutOfRange {
            index: codebook_index,
        });
    }
    Ok(GAIN_QUANT_GA_PERMUTATION[codebook_index] as usize)
}

/// Encoder-side §3.9.3 forward permutation for the second-stage gain
/// codebook: given the **codebook** index in `0..NCODE2`, return the
/// **transmitted** index that goes on the wire.
///
/// # Errors
///
/// Returns [`GainIndexMapError::GbOutOfRange`] if `codebook_index`
/// does not fit in `0..NCODE2`.
pub fn map_gb(codebook_index: usize) -> Result<usize, GainIndexMapError> {
    if codebook_index >= NCODE2 {
        return Err(GainIndexMapError::GbOutOfRange {
            index: codebook_index,
        });
    }
    Ok(GAIN_QUANT_GB_PERMUTATION[codebook_index] as usize)
}

/// Decoder-side §3.9.3 inverse permutation for the first-stage gain
/// codebook: given the **transmitted** index parsed off the wire,
/// return the **codebook** index that the §3.9.2 reconstruction
/// indexes the GA codebook with.
///
/// Implemented by indexing the [`GAIN_QUANT_GA_INVERSE_PERMUTATION`]
/// table at `transmitted_index`. The companion [`map_ga`] is the
/// encode-side inverse; round-tripping
/// `demap_ga(map_ga(i)?) == Ok(i)` for every `i` in `0..NCODE1` is
/// pinned by the unit tests below.
///
/// # Errors
///
/// Returns [`GainIndexMapError::GaOutOfRange`] if `transmitted_index`
/// does not fit in `0..NCODE1`. The 3-bit `GA` field of spec Table 8
/// guarantees the input is in range for a well-formed frame.
pub fn demap_ga(transmitted_index: usize) -> Result<usize, GainIndexMapError> {
    if transmitted_index >= NCODE1 {
        return Err(GainIndexMapError::GaOutOfRange {
            index: transmitted_index,
        });
    }
    Ok(GAIN_QUANT_GA_INVERSE_PERMUTATION[transmitted_index] as usize)
}

/// Decoder-side §3.9.3 inverse permutation for the second-stage gain
/// codebook: given the **transmitted** index parsed off the wire,
/// return the **codebook** index for the §3.9.2 reconstruction.
///
/// # Errors
///
/// Returns [`GainIndexMapError::GbOutOfRange`] if `transmitted_index`
/// does not fit in `0..NCODE2`. The 4-bit `GB` field of spec Table 8
/// guarantees the input is in range for a well-formed frame.
pub fn demap_gb(transmitted_index: usize) -> Result<usize, GainIndexMapError> {
    if transmitted_index >= NCODE2 {
        return Err(GainIndexMapError::GbOutOfRange {
            index: transmitted_index,
        });
    }
    Ok(GAIN_QUANT_GB_INVERSE_PERMUTATION[transmitted_index] as usize)
}

/// Per-frame decoder-side wrapper: take the four transmitted gain-VQ
/// indices from [`Parameters`] and return a typed
/// [`DemappedGainIndices`] with each index in the codebook domain.
///
/// The two subframes are handled independently — there is no
/// cross-subframe coupling in the §3.9.3 mapping. The returned
/// struct's field ordering matches spec §4.1.5 (subframe 1 first).
///
/// # Errors
///
/// Surfaces the first out-of-range index in
/// `ga1 → gb1 → ga2 → gb2` evaluation order. The 3 + 4 + 3 + 4 = 14
/// transmitted gain bits per spec Table 8 cannot drive any of these
/// out of range in a well-formed frame.
pub fn demap_frame(params: &Parameters) -> Result<DemappedGainIndices, GainIndexMapError> {
    let ga1 = demap_ga(usize::from(params.ga1))?;
    let gb1 = demap_gb(usize::from(params.gb1))?;
    let ga2 = demap_ga(usize::from(params.ga2))?;
    let gb2 = demap_gb(usize::from(params.gb2))?;
    Ok(DemappedGainIndices { ga1, gb1, ga2, gb2 })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every transmitted GA value in `0..NCODE1` round-trips through
    /// `map_ga ∘ demap_ga` to itself (encoder applies `map_ga`,
    /// decoder applies `demap_ga`).
    #[test]
    fn map_then_demap_round_trips_ga() {
        for i in 0..NCODE1 {
            let transmitted = map_ga(i).expect("in range");
            let recovered = demap_ga(transmitted).expect("in range");
            assert_eq!(recovered, i);
        }
    }

    /// Every transmitted GB value in `0..NCODE2` round-trips through
    /// `map_gb ∘ demap_gb` to itself.
    #[test]
    fn map_then_demap_round_trips_gb() {
        for i in 0..NCODE2 {
            let transmitted = map_gb(i).expect("in range");
            let recovered = demap_gb(transmitted).expect("in range");
            assert_eq!(recovered, i);
        }
    }

    /// And the other direction: every received GA value round-trips
    /// through `demap_ga ∘ map_ga` to itself. Pinned independently of
    /// the encoder-side test because the two compositions are only
    /// equal when the two permutations are true mutual inverses on
    /// the same domain.
    #[test]
    fn demap_then_map_round_trips_ga() {
        for i in 0..NCODE1 {
            let codebook = demap_ga(i).expect("in range");
            let transmitted = map_ga(codebook).expect("in range");
            assert_eq!(transmitted, i);
        }
    }

    /// Same property for the second-stage codebook.
    #[test]
    fn demap_then_map_round_trips_gb() {
        for i in 0..NCODE2 {
            let codebook = demap_gb(i).expect("in range");
            let transmitted = map_gb(codebook).expect("in range");
            assert_eq!(transmitted, i);
        }
    }

    /// Every output of `demap_ga` lies in the codebook domain
    /// `0..NCODE1` — guards against a CSV reader regression that
    /// silently widens the staged table values past the codebook
    /// dimension.
    #[test]
    fn demap_ga_output_stays_in_codebook_domain() {
        for i in 0..NCODE1 {
            let codebook = demap_ga(i).expect("in range");
            assert!(
                codebook < NCODE1,
                "demap_ga({i}) = {codebook} >= NCODE1 = {NCODE1}",
            );
        }
    }

    /// Same boundary check on the second-stage demap.
    #[test]
    fn demap_gb_output_stays_in_codebook_domain() {
        for i in 0..NCODE2 {
            let codebook = demap_gb(i).expect("in range");
            assert!(
                codebook < NCODE2,
                "demap_gb({i}) = {codebook} >= NCODE2 = {NCODE2}",
            );
        }
    }

    /// `demap_ga` is a bijection on `0..NCODE1` — every codebook
    /// index appears exactly once in the image. Equivalent to a
    /// "cover" check on the staged `imap1` table.
    #[test]
    fn demap_ga_covers_codebook_domain() {
        let mut seen = [false; NCODE1];
        for i in 0..NCODE1 {
            let codebook = demap_ga(i).expect("in range");
            assert!(
                !seen[codebook],
                "demap_ga maps two transmitted indices to codebook index {codebook}",
            );
            seen[codebook] = true;
        }
        assert!(seen.iter().all(|&s| s), "demap_ga does not cover NCODE1");
    }

    /// `demap_gb` is a bijection on `0..NCODE2`.
    #[test]
    fn demap_gb_covers_codebook_domain() {
        let mut seen = [false; NCODE2];
        for i in 0..NCODE2 {
            let codebook = demap_gb(i).expect("in range");
            assert!(
                !seen[codebook],
                "demap_gb maps two transmitted indices to codebook index {codebook}",
            );
            seen[codebook] = true;
        }
        assert!(seen.iter().all(|&s| s), "demap_gb does not cover NCODE2");
    }

    /// Out-of-range GA index yields the typed error rather than
    /// panicking. The boundary value `NCODE1` is just past the last
    /// in-domain index `NCODE1 - 1`.
    #[test]
    fn map_ga_rejects_out_of_range() {
        let err = map_ga(NCODE1).unwrap_err();
        assert_eq!(err, GainIndexMapError::GaOutOfRange { index: NCODE1 });
    }

    /// Same boundary check on the decode-side primitive.
    #[test]
    fn demap_ga_rejects_out_of_range() {
        let err = demap_ga(NCODE1).unwrap_err();
        assert_eq!(err, GainIndexMapError::GaOutOfRange { index: NCODE1 });
    }

    /// Out-of-range GB index yields the typed error rather than
    /// panicking.
    #[test]
    fn map_gb_rejects_out_of_range() {
        let err = map_gb(NCODE2).unwrap_err();
        assert_eq!(err, GainIndexMapError::GbOutOfRange { index: NCODE2 });
    }

    /// Same boundary check on the decode-side primitive.
    #[test]
    fn demap_gb_rejects_out_of_range() {
        let err = demap_gb(NCODE2).unwrap_err();
        assert_eq!(err, GainIndexMapError::GbOutOfRange { index: NCODE2 });
    }

    /// The §3.9.3 mapping is NOT the identity — at least one input is
    /// reshuffled. Locks the mapping against a silent CSV regression
    /// that emits `[0, 1, 2, ...]` for either of the two `imap`
    /// tables.
    #[test]
    fn demap_ga_is_not_identity() {
        let any_moved = (0..NCODE1).any(|i| demap_ga(i).unwrap() != i);
        assert!(any_moved, "demap_ga must be a non-identity permutation");
    }

    /// Same non-identity property for the GB stage.
    #[test]
    fn demap_gb_is_not_identity() {
        let any_moved = (0..NCODE2).any(|i| demap_gb(i).unwrap() != i);
        assert!(any_moved, "demap_gb must be a non-identity permutation");
    }

    /// Per-frame wrapper threads each subframe's pair through the
    /// corresponding scalar primitive. Crafted [`Parameters`] with
    /// four distinct gain indices ensures the wrapper doesn't
    /// accidentally reuse subframe-1 values for subframe 2 or swap
    /// the GA / GB columns.
    #[test]
    fn frame_wrapper_threads_per_subframe_indices() {
        let params = Parameters {
            l0: 0,
            l1: 0,
            l2: 0,
            l3: 0,
            p1: 0,
            p0: 0,
            c1: 0,
            s1: 0,
            ga1: 1,
            gb1: 2,
            p2: 0,
            c2: 0,
            s2: 0,
            ga2: 6,
            gb2: 13,
        };
        let out = demap_frame(&params).expect("all indices in range");
        assert_eq!(out.ga1, demap_ga(1).unwrap());
        assert_eq!(out.gb1, demap_gb(2).unwrap());
        assert_eq!(out.ga2, demap_ga(6).unwrap());
        assert_eq!(out.gb2, demap_gb(13).unwrap());
    }

    /// Sanity-pin a representative literal: the wire-index `0`
    /// demaps to `imap1[0]` (whatever the staged CSV says). The
    /// concrete value is read from the same table the production
    /// path uses so the test does not double-bake the literal;
    /// instead it ties the function's contract directly to the
    /// staged table.
    #[test]
    fn demap_ga_at_zero_matches_imap1_zero() {
        let expected = GAIN_QUANT_GA_INVERSE_PERMUTATION[0] as usize;
        assert_eq!(demap_ga(0).unwrap(), expected);
    }

    /// Same shape for the second-stage codebook.
    #[test]
    fn demap_gb_at_zero_matches_imap2_zero() {
        let expected = GAIN_QUANT_GB_INVERSE_PERMUTATION[0] as usize;
        assert_eq!(demap_gb(0).unwrap(), expected);
    }
}
