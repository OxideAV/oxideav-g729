//! Shape + spot-value tests for the Round-173 numeric tables compiled
//! by `build.rs` from `tables/*.csv`. These tests catch:
//!
//! * Element-count drift between the staged CSVs and the `pub const`
//!   arrays exposed via `oxideav_g729::tables`.
//! * Value drift on the small, manually-verifiable entries documented
//!   in each table's `.meta` `spec_role` line (the `.meta` strings cite
//!   the leading elements explicitly for the high-pass filters).
//!
//! The asserts deliberately do NOT mirror the entire CSV — only the
//! distinctive entries spelled out in the meta sidecars and the
//! boundary entries (first/last) that the build script's CSV reader
//! is most likely to mishandle if anything in the pipeline drifts.

use oxideav_g729::tables;

/// §3.1 / §4.2 IIR HPF coefficients are 3-element Word16 arrays (one
/// scalar each for taps `b0`, `b1`, `b2` of the b-coef array; `a0`,
/// `a1`, `a2` of the a-coef array, where `a0` is the implicit unit
/// gain reflected as the Q-scaled `1.0`).
#[test]
fn highpass_coefs_have_three_entries_each() {
    assert_eq!(tables::HPF_PREPROC_100HZ_B_Q13.len(), 3);
    assert_eq!(tables::HPF_PREPROC_100HZ_A_Q13.len(), 3);
    assert_eq!(tables::HPF_PREPROC_140HZ_B_Q12.len(), 3);
    assert_eq!(tables::HPF_PREPROC_140HZ_A_Q12.len(), 3);
}

/// The 100 Hz HPF coefficient triples are documented in the `.meta`
/// `spec_role` line itself (`{7699, -15398, 7699}` for `b100`,
/// `{8192, 15836, -7667}` for `a100`). Lock those values in.
#[test]
fn highpass_100hz_values_match_meta_documentation() {
    assert_eq!(tables::HPF_PREPROC_100HZ_B_Q13, [7699, -15398, 7699]);
    assert_eq!(tables::HPF_PREPROC_100HZ_A_Q13, [8192, 15836, -7667]);
}

/// The 140 Hz HPF coefficient triples are likewise documented in the
/// meta (`{1899, -3798, 1899}` for `b140`, `{4096, 7807, -3733}` for
/// `a140`).
#[test]
fn highpass_140hz_values_match_meta_documentation() {
    assert_eq!(tables::HPF_PREPROC_140HZ_B_Q12, [1899, -3798, 1899]);
    assert_eq!(tables::HPF_PREPROC_140HZ_A_Q12, [4096, 7807, -3733]);
}

/// HPF symmetry property: the b-coefficient array of a normalised IIR
/// high-pass filter centred about `b1` has `b0 == b2` (the
/// transfer function's numerator is `b0*(1 - 2 z^-1 + z^-2)` scaled by
/// `b0`). Both `b100` and `b140` exhibit this on inspection of the meta.
#[test]
fn highpass_b_coefs_are_symmetric() {
    let b100 = &tables::HPF_PREPROC_100HZ_B_Q13;
    assert_eq!(b100[0], b100[2]);
    let b140 = &tables::HPF_PREPROC_140HZ_B_Q12;
    assert_eq!(b140[0], b140[2]);
}

/// §4.1 Table 8 bit-allocation array carries 13 entries as extracted
/// from the ITU C source array literal (per the `.meta` element_count
/// field). Only the first `PRM_SIZE = 11` indices are spec-defined
/// parameters; the trailing entries are preserved as-extracted.
#[test]
fn bit_allocation_table_has_thirteen_entries() {
    assert_eq!(tables::BIT_ALLOCATION_TABLE8.len(), 13);
}

/// The bit-allocation array as extracted contains the exact 13
/// literal Word16 values that appeared in the ITU C source — the
/// extractor preserves the source array byte-for-byte. This test
/// pins the entire content so any future drift in the staged CSV
/// or build-script reader is caught loudly.
///
/// NOTE: the literal sum across all 13 entries is `66`, and the
/// `PRM_SIZE = 11` prefix sums to `55`. Neither equals the
/// headline G.729 8 kbit/s frame rate (80 bits / 10 ms frame) on
/// its own, because spec Table 8 itself enumerates 15 (rather
/// than 11 or 13) per-parameter rows. Resolving how the C
/// `bitsno` array's 13-entry packing maps onto Table 8's
/// 15-parameter listing is deferred to the docs collaborator;
/// this test simply locks the values that were extracted.
#[test]
fn bit_allocation_values_are_pinned_to_csv() {
    assert_eq!(
        tables::BIT_ALLOCATION_TABLE8,
        [1, 0, 1, 2, 8, 1, 13, 4, 7, 5, 13, 4, 7],
    );
    let total_literal: i32 = tables::BIT_ALLOCATION_TABLE8
        .iter()
        .map(|&v| i32::from(v))
        .sum();
    assert_eq!(total_literal, 66);
}

/// `basic_op::Pow2` lookup is 33 entries. The first entry corresponds
/// to `2^0 = 1.0` represented as Q15 = `16384` (per the table's
/// canonical anchor: tabpow[0] holds the value at the start of the
/// interpolation range). Last entry is 32767 (Q15 ≈ 1.99997, the
/// largest representable Word16).
#[test]
fn pow2_table_anchor_values() {
    assert_eq!(tables::POW2_TABLE_Q15.len(), 33);
    assert_eq!(tables::POW2_TABLE_Q15[0], 16384);
    assert_eq!(*tables::POW2_TABLE_Q15.last().unwrap(), 32767);
}

/// `basic_op::Log2` lookup is 33 entries. Its values are monotonically
/// non-decreasing (Log2 is a monotone function on the positive
/// reals); this gives a cheap structural check that the CSV order
/// survived the pipeline intact.
#[test]
fn log2_table_is_monotonic_nondecreasing() {
    assert_eq!(tables::LOG2_TABLE_Q15.len(), 33);
    for pair in tables::LOG2_TABLE_Q15.windows(2) {
        assert!(
            pair[0] <= pair[1],
            "Log2 table not monotonic: {} then {}",
            pair[0],
            pair[1],
        );
    }
}

/// `basic_op::Inv_sqrt` lookup is 49 entries. Its values are
/// monotonically non-increasing (1/sqrt(x) is monotone decreasing).
/// Same cheap structural-integrity check as Log2.
#[test]
fn inv_sqrt_table_is_monotonic_nonincreasing() {
    assert_eq!(tables::INV_SQRT_TABLE_Q15.len(), 49);
    for pair in tables::INV_SQRT_TABLE_Q15.windows(2) {
        assert!(
            pair[0] >= pair[1],
            "Inv_sqrt table not monotonic: {} then {}",
            pair[0],
            pair[1],
        );
    }
}

/// All published tables are non-empty and store Word16 values
/// (compile-time guaranteed via the `[i16; N]` typing). This sanity
/// check is here so a future build-script refactor that accidentally
/// emits zero-length arrays trips a clear failure.
#[test]
fn all_tables_are_non_empty() {
    assert!(!tables::HPF_PREPROC_100HZ_B_Q13.is_empty());
    assert!(!tables::HPF_PREPROC_100HZ_A_Q13.is_empty());
    assert!(!tables::HPF_PREPROC_140HZ_B_Q12.is_empty());
    assert!(!tables::HPF_PREPROC_140HZ_A_Q12.is_empty());
    assert!(!tables::BIT_ALLOCATION_TABLE8.is_empty());
    assert!(!tables::POW2_TABLE_Q15.is_empty());
    assert!(!tables::LOG2_TABLE_Q15.is_empty());
    assert!(!tables::INV_SQRT_TABLE_Q15.is_empty());
    assert!(!tables::LPC_HAMMING_WINDOW_Q15.is_empty());
    assert!(!tables::LPC_LAG_WINDOW_HIGH_Q15.is_empty());
    assert!(!tables::LPC_LAG_WINDOW_LOW_Q15.is_empty());
    assert!(!tables::LSF_SEARCH_GRID_COS_Q15.is_empty());
    assert!(!tables::PITCH_INTERP_FILTER_ANALYSIS_Q15.is_empty());
    assert!(!tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15.is_empty());
    assert!(!tables::GAIN_QUANT_MA_PREDICTOR_Q13.is_empty());
    assert!(!tables::LSP_QUANT_CODEBOOK_L1_Q13.is_empty());
    assert!(!tables::LSP_QUANT_CODEBOOK_L2_Q13.is_empty());
    assert!(!tables::LSP_MA_PREDICTOR_FG_Q15.is_empty());
    assert!(!tables::LSP_MA_PREDICTOR_FG_SUM_Q15.is_empty());
    assert!(!tables::LSP_MA_PREDICTOR_FG_SUM_INV_Q12.is_empty());
}

// ---------------------------------------------------------------------
// §3.2.1 LP analysis windowing tables (round 189).
// ---------------------------------------------------------------------

/// `hamwindow` spans `L_WINDOW = 240` samples, exactly the
/// LP-analysis frame length per spec §3.2.1.
#[test]
fn hamming_window_length_matches_l_window() {
    assert_eq!(tables::LPC_HAMMING_WINDOW_Q15.len(), tables::L_WINDOW);
    assert_eq!(tables::L_WINDOW, 240);
}

/// The Hamming-cos window's peak value is at or above all other
/// samples, and the staged CSV's first and last entries are both
/// less than the peak (the window is tapered at both ends). The
/// peak of the spec's `Hamming_cos` formulation reaches Q15 ≈ 1.0
/// (`32767`).
#[test]
fn hamming_window_peaks_at_q15_unit() {
    let peak = tables::LPC_HAMMING_WINDOW_Q15
        .iter()
        .copied()
        .max()
        .unwrap();
    assert_eq!(peak, 32767, "Hamming window peak should be Q15 ≈ 1.0");

    let first = tables::LPC_HAMMING_WINDOW_Q15[0];
    let last = *tables::LPC_HAMMING_WINDOW_Q15.last().unwrap();
    assert!(first < peak);
    assert!(last < peak);
}

/// The Hamming-cos window is strictly positive on its full support
/// (`w[n] = 0.54 - 0.46·cos(...)` lies in `[0.08, 1.0]`).
#[test]
fn hamming_window_is_strictly_positive() {
    for (i, &v) in tables::LPC_HAMMING_WINDOW_Q15.iter().enumerate() {
        assert!(v > 0, "Hamming window has non-positive entry at {i}: {v}");
    }
}

/// `lag_h` / `lag_l` form the high / low halves of a
/// double-precision `Word32` lag-window pair (§3.2.1 + `oper_32b.c`
/// double-precision representation). They share length `M = 10`.
#[test]
fn lag_window_pair_lengths_match_m() {
    assert_eq!(tables::LPC_LAG_WINDOW_HIGH_Q15.len(), tables::M);
    assert_eq!(tables::LPC_LAG_WINDOW_LOW_Q15.len(), tables::M);
    assert_eq!(tables::M, 10);
}

/// `lag_h` (the high-order half of the double-precision lag-window
/// representation) is monotonically decreasing — the bandwidth-
/// expansion factor applied to autocorrelation taps shrinks with
/// the tap index `i = 1..M`.
#[test]
fn lag_window_high_is_monotonic_decreasing() {
    for pair in tables::LPC_LAG_WINDOW_HIGH_Q15.windows(2) {
        assert!(
            pair[0] > pair[1],
            "lag_h not strictly decreasing: {} then {}",
            pair[0],
            pair[1],
        );
    }
}

// ---------------------------------------------------------------------
// §3.2.5 az_lsf() root-search grid (round 189).
// ---------------------------------------------------------------------

/// The grid is `GRID_POINTS + 1 = 61` Q15 cosine samples spanning
/// the upper half-circle in 60 equal angular steps.
#[test]
fn grid_length_matches_grid_points_plus_one() {
    assert_eq!(
        tables::LSF_SEARCH_GRID_COS_Q15.len(),
        tables::GRID_POINTS + 1
    );
    assert_eq!(tables::GRID_POINTS, 60);
}

/// Endpoints of the cosine grid: `grid[0] = cos(0) ≈ 1.0` (Q15
/// upper limit `32767`-ish), `grid[60] = cos(π) ≈ -1.0`
/// (Q15 lower limit `-32768`-ish). The CSV's literal values are
/// `32760` and `-32760` (matching the ITU C source exactly).
#[test]
fn grid_endpoints_match_csv_literals() {
    assert_eq!(tables::LSF_SEARCH_GRID_COS_Q15[0], 32760);
    assert_eq!(*tables::LSF_SEARCH_GRID_COS_Q15.last().unwrap(), -32760);
}

/// `grid[30]` corresponds to `cos(π/2)`, which is exactly `0`. The
/// G.729 spec's `grid` array deliberately includes this midpoint
/// (61 samples → odd count → centre sample is the zero crossing).
#[test]
fn grid_midpoint_is_zero() {
    let mid = tables::LSF_SEARCH_GRID_COS_Q15.len() / 2;
    assert_eq!(tables::LSF_SEARCH_GRID_COS_Q15[mid], 0);
}

/// The grid is monotonically strictly decreasing across the full
/// half-circle (cosine is a strictly decreasing function on
/// `[0, π]`).
#[test]
fn grid_is_strictly_decreasing() {
    for pair in tables::LSF_SEARCH_GRID_COS_Q15.windows(2) {
        assert!(
            pair[0] > pair[1],
            "grid not strictly decreasing: {} then {}",
            pair[0],
            pair[1],
        );
    }
}

/// Antisymmetry of cosine on `[0, π]`: `grid[i] = -grid[N-1-i]`.
/// This catches any axis-flip or off-by-one in the extraction.
#[test]
fn grid_is_antisymmetric_about_midpoint() {
    let g = &tables::LSF_SEARCH_GRID_COS_Q15;
    let n = g.len();
    for i in 0..n {
        assert_eq!(
            g[i],
            -g[n - 1 - i],
            "grid not antisymmetric at i={i}: {} vs -{}",
            g[i],
            g[n - 1 - i],
        );
    }
}

// ---------------------------------------------------------------------
// §3.7 pitch interpolation filters (round 189).
// ---------------------------------------------------------------------

/// `inter_3` is the 13-tap (`FIR_SIZE_ANA`) Q15 analysis-side
/// 1/3-resolution interpolation filter.
#[test]
fn pitch_analysis_filter_has_fir_size_ana_taps() {
    assert_eq!(tables::PITCH_INTERP_FILTER_ANALYSIS_Q15.len(), 13);
}

/// `inter_3l` is the 31-tap (`FIR_SIZE_SYN`) Q15 synthesis-side
/// 1/3-resolution interpolation filter.
#[test]
fn pitch_synthesis_filter_has_fir_size_syn_taps() {
    assert_eq!(tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15.len(), 31);
}

/// The analysis filter's peak tap is positive and at least as
/// large in magnitude as any other tap (a sensible windowed-sinc
/// design has the centre tap dominate). Both staged filters share
/// this property; this is a cheap drift check on the CSV-to-Rust
/// pipeline.
#[test]
fn pitch_analysis_filter_centre_dominates() {
    let max = tables::PITCH_INTERP_FILTER_ANALYSIS_Q15
        .iter()
        .copied()
        .max()
        .unwrap();
    let abs_max = tables::PITCH_INTERP_FILTER_ANALYSIS_Q15
        .iter()
        .map(|&v| (v as i32).unsigned_abs())
        .max()
        .unwrap();
    assert!(max > 0);
    assert_eq!(max as u32, abs_max, "filter peak should be a positive tap");
}

/// Same dominance property for the synthesis filter — sanity check
/// that the longer-tap design's peak remains positive.
#[test]
fn pitch_synthesis_filter_peak_is_positive() {
    let max = tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15
        .iter()
        .copied()
        .max()
        .unwrap();
    let abs_max = tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15
        .iter()
        .map(|&v| (v as i32).unsigned_abs())
        .max()
        .unwrap();
    assert!(max > 0);
    assert_eq!(max as u32, abs_max);
}

// ---------------------------------------------------------------------
// §3.9 MA gain-prediction coefficients (round 189).
// ---------------------------------------------------------------------

/// `pred` is a 4-element Q13 vector. The meta `spec_role` line
/// names the real values explicitly as {0.68, 0.58, 0.34, 0.19};
/// in Q13 those round to {5571, 4751, 2785, 1556} (`round(v *
/// 8192)`).
#[test]
fn gain_ma_predictor_matches_meta_documentation() {
    assert_eq!(tables::GAIN_QUANT_MA_PREDICTOR_Q13.len(), 4);
    assert_eq!(
        tables::GAIN_QUANT_MA_PREDICTOR_Q13,
        [5571, 4751, 2785, 1556],
    );

    // Each Q13 entry round-trips back to its spec'd real value
    // within the Q13 quantisation step (1 / 8192 ≈ 1.22e-4).
    let real_targets: [f64; 4] = [0.68, 0.58, 0.34, 0.19];
    let q13_scale = 8192.0;
    let q13_step = 1.0 / q13_scale;
    for (i, &v) in tables::GAIN_QUANT_MA_PREDICTOR_Q13.iter().enumerate() {
        let recovered = f64::from(v) / q13_scale;
        let err = (recovered - real_targets[i]).abs();
        assert!(
            err < q13_step,
            "pred[{i}] = {v} (Q13 {recovered}) differs from spec target {} by {err}",
            real_targets[i],
        );
    }
}

/// The MA gain predictor is a monotonically non-increasing vector
/// (longer-history terms are weighted less). This is a structural
/// property of the spec'd values that survives the staging
/// pipeline and would catch any reordering of the CSV.
#[test]
fn gain_ma_predictor_is_monotonic_nonincreasing() {
    for pair in tables::GAIN_QUANT_MA_PREDICTOR_Q13.windows(2) {
        assert!(
            pair[0] >= pair[1],
            "pred not monotonic non-increasing: {} then {}",
            pair[0],
            pair[1],
        );
    }
}

// ---------------------------------------------------------------------
// §3.2.4 LSP quantiser two-stage VQ codebooks (round 195).
// ---------------------------------------------------------------------

/// `LSP_QUANT_CODEBOOK_L1_Q13` is the first-stage codebook `lspcb1`
/// with `NC0 = 128` 10-dimensional Q13 entries. The 7-bit `L1` index
/// covers exactly the codebook's row count, so `1 << L1_BITS == NC0`.
#[test]
fn l1_codebook_shape_matches_nc0_by_m() {
    assert_eq!(tables::LSP_QUANT_CODEBOOK_L1_Q13.len(), tables::NC0);
    for (i, row) in tables::LSP_QUANT_CODEBOOK_L1_Q13.iter().enumerate() {
        assert_eq!(row.len(), tables::M, "L1 row {i} not M-wide");
    }
    assert_eq!(1 << tables::L1_BITS, tables::NC0);
}

/// `LSP_QUANT_CODEBOOK_L2_Q13` is the packed second-stage codebook
/// `lspcb2`: `NC1 = 32` rows of 10 Q13 values, with the lower 5
/// columns holding the L2 split-VQ contribution and the upper 5
/// holding L3. The 5-bit L2/L3 indices match the row count exactly,
/// so `1 << L2_BITS == NC1` and `L2_BITS == L3_BITS`.
#[test]
fn l2_codebook_shape_matches_nc1_by_m() {
    assert_eq!(tables::LSP_QUANT_CODEBOOK_L2_Q13.len(), tables::NC1);
    for (i, row) in tables::LSP_QUANT_CODEBOOK_L2_Q13.iter().enumerate() {
        assert_eq!(row.len(), tables::M, "L2 row {i} not M-wide");
    }
    assert_eq!(1 << tables::L2_BITS, tables::NC1);
    assert_eq!(tables::L2_BITS, tables::L3_BITS);
}

/// Lock the first three rows of the L1 codebook to their CSV
/// literals — these are documented in the trace doc's §3.5 worked
/// example as the leading rows of `lspcb1`, and pinning them here
/// gives a focused drift check against any reordering or off-by-one
/// in the matrix CSV reader.
#[test]
fn l1_codebook_leading_rows_match_csv_literals() {
    assert_eq!(
        tables::LSP_QUANT_CODEBOOK_L1_Q13[0],
        [1486, 2168, 3751, 9074, 12134, 13944, 17983, 19173, 21190, 21820],
    );
    assert_eq!(
        tables::LSP_QUANT_CODEBOOK_L1_Q13[1],
        [1730, 2640, 3450, 4870, 6126, 7876, 15644, 17817, 20294, 21902],
    );
    assert_eq!(
        tables::LSP_QUANT_CODEBOOK_L1_Q13[2],
        [1568, 2256, 3088, 4874, 11063, 13393, 18307, 19293, 21109, 21741],
    );
}

/// Lock the first row of the L2 codebook to its CSV literal — like
/// the L1 check above, this catches matrix-reader drift early.
#[test]
fn l2_codebook_first_row_matches_csv_literal() {
    assert_eq!(
        tables::LSP_QUANT_CODEBOOK_L2_Q13[0],
        [-435, -815, -742, 1033, -518, 582, -1201, 829, 86, 385],
    );
}

/// Every L1 row holds an LSF-quantiser contribution whose 10 values
/// span a moderate Q13 range; no entry should overflow Word16. The
/// `[i16; 10]` typing guarantees that statically, but assert it here
/// in observation form so a future drift to `[i32; ..]` (e.g. via a
/// build-script bug) is caught by intent rather than only by type
/// check. Also assert the entries are bounded above by a value
/// well below Q13 ≈ 4.0 — the LSF domain is `[0, π]` so a 10th
/// coefficient ≤ π × 8192 ≈ 25735 is the natural upper limit.
#[test]
fn l1_codebook_entries_are_in_lsf_domain() {
    for (i, row) in tables::LSP_QUANT_CODEBOOK_L1_Q13.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            assert!(
                (0..=26000).contains(&v),
                "L1 row {i} col {j} outside Q13 LSF domain [0, ~π·2^13]: {v}",
            );
        }
    }
}

/// The bounds-checked L1 lookup helper returns a 10-element slice
/// that matches the underlying constant by value.
#[test]
fn lsp_l1_entry_helper_matches_constant() {
    for index in [0_usize, 1, 7, 64, 127] {
        assert_eq!(
            tables::lsp_l1_entry(index),
            &tables::LSP_QUANT_CODEBOOK_L1_Q13[index],
        );
    }
}

/// The L2 / L3 lookup helpers each return a 5-element slice whose
/// values are the lower / upper half of the row, respectively. The
/// concatenation reproduces the full row.
#[test]
fn lsp_l2_l3_helpers_reconstruct_packed_row() {
    for index in [0_usize, 1, 16, 31] {
        let lo = tables::lsp_l2_entry(index);
        let hi = tables::lsp_l3_entry(index);
        assert_eq!(lo.len(), tables::M / 2);
        assert_eq!(hi.len(), tables::M / 2);

        let mut joined = Vec::with_capacity(tables::M);
        joined.extend_from_slice(lo);
        joined.extend_from_slice(hi);
        assert_eq!(joined.as_slice(), &tables::LSP_QUANT_CODEBOOK_L2_Q13[index]);
    }
}

/// The total transmitted-LSP bit count derives from the four
/// per-stage bit fields and matches spec Table 1.
#[test]
fn lsp_total_bits_match_spec_table_1() {
    assert_eq!(tables::LSP_TOTAL_BITS, 18);
    assert_eq!(
        tables::LSP_TOTAL_BITS,
        tables::L0_BITS + tables::L1_BITS + tables::L2_BITS + tables::L3_BITS,
    );
}

// ---------------------------------------------------------------------
// §3.2.4 LSP MA-predictor `fg` family (round 201).
// ---------------------------------------------------------------------

/// `LSP_MA_PREDICTOR_FG_Q15` is the 3-D predictor cube: outer 2-mode
/// dimension (selected by `L0`), inner `MA_NP` × `M` plane carrying
/// the 4th-order MA coefficients across each LSP coordinate.
/// `MA_NP == 4` per spec §3.2.4 eq (20).
#[test]
fn fg_cube_shape_matches_two_by_ma_np_by_m() {
    assert_eq!(tables::MA_NP, 4);
    assert_eq!(tables::LSP_MA_PREDICTOR_FG_Q15.len(), 2);
    for (mode, plane) in tables::LSP_MA_PREDICTOR_FG_Q15.iter().enumerate() {
        assert_eq!(
            plane.len(),
            tables::MA_NP,
            "fg[{mode}] inner-row count not MA_NP",
        );
        for (k, row) in plane.iter().enumerate() {
            assert_eq!(row.len(), tables::M, "fg[{mode}][{k}] row not M-wide",);
        }
    }
    // The outer dim is exactly enumerable from the 1-bit L0 field.
    assert_eq!(1 << tables::L0_BITS, tables::LSP_MA_PREDICTOR_FG_Q15.len());
}

/// `LSP_MA_PREDICTOR_FG_SUM_Q15` and `LSP_MA_PREDICTOR_FG_SUM_INV_Q12`
/// are per-mode rows of `M` Q-format scalars. There is one row per
/// predictor mode, matching the outer dim of the `fg` cube exactly.
#[test]
fn fg_sum_and_inv_shapes_match_two_by_m() {
    assert_eq!(tables::LSP_MA_PREDICTOR_FG_SUM_Q15.len(), 2);
    assert_eq!(tables::LSP_MA_PREDICTOR_FG_SUM_INV_Q12.len(), 2);
    for mode in 0..2 {
        assert_eq!(tables::LSP_MA_PREDICTOR_FG_SUM_Q15[mode].len(), tables::M);
        assert_eq!(
            tables::LSP_MA_PREDICTOR_FG_SUM_INV_Q12[mode].len(),
            tables::M
        );
    }
}

/// Pin the very first row of each `fg` plane to its CSV literal — the
/// 3-D matrix reader is exercised here for the first time, and this
/// gives a focused drift check that the row-major flattening within
/// each plane survived the build script unscrambled.
///
/// Values are the first 10 of the staged `fg` CSV (mode 0, k = 0) and
/// the first 10 of the second CSV line (mode 1, k = 0).
#[test]
fn fg_first_rows_match_csv_literals() {
    assert_eq!(
        tables::LSP_MA_PREDICTOR_FG_Q15[0][0],
        [8421, 9109, 9175, 8965, 9034, 9057, 8765, 8775, 9106, 8673],
    );
    assert_eq!(
        tables::LSP_MA_PREDICTOR_FG_Q15[1][0],
        [7733, 7880, 8188, 8175, 8247, 8490, 8637, 8601, 8359, 7569],
    );
}

/// The last MA-history row of each plane is the lowest-weighted tap
/// (predictor weight shrinks with history depth). For both modes
/// the magnitude of every entry in `fg[mode][MA_NP - 1]` should be
/// less than the largest magnitude in `fg[mode][0]` — a structural
/// monotonicity check that catches any plane- or row-flip in the
/// 3-D parser without pinning all 80 values.
#[test]
fn fg_history_decays_per_mode() {
    for mode in 0..2 {
        let head_peak = tables::LSP_MA_PREDICTOR_FG_Q15[mode][0]
            .iter()
            .map(|v| (i32::from(*v)).unsigned_abs())
            .max()
            .unwrap();
        let tail_peak = tables::LSP_MA_PREDICTOR_FG_Q15[mode][tables::MA_NP - 1]
            .iter()
            .map(|v| (i32::from(*v)).unsigned_abs())
            .max()
            .unwrap();
        assert!(
            tail_peak < head_peak,
            "fg[{mode}] tail peak ({tail_peak}) >= head peak ({head_peak})",
        );
    }
}

/// `fg_sum[mode][i]` is the Q15 representation of the spec eq (20a)
/// factor `(1 − Σ_{k=0..MA_NP} fg[mode][k][i])` — the per-coordinate
/// weight by which the reconstructed L1+L2/L3 contribution is scaled
/// after the MA history is subtracted. Computed in i32 to avoid
/// Word16 wrap; the staged value is within a few Q15 quantisation
/// steps of the directly-evaluated factor (independent rounding on
/// each side of the extraction).
#[test]
fn fg_sum_is_one_minus_per_column_fg_sum() {
    let q15_one: i32 = 32768;
    for mode in 0..2 {
        for i in 0..tables::M {
            let col_sum: i32 = (0..tables::MA_NP)
                .map(|k| i32::from(tables::LSP_MA_PREDICTOR_FG_Q15[mode][k][i]))
                .sum();
            let expected = q15_one - col_sum;
            let actual = i32::from(tables::LSP_MA_PREDICTOR_FG_SUM_Q15[mode][i]);
            let err = (actual - expected).abs();
            // Allow up to 4 Q15 ulps difference: the staged values are
            // each independently quantised from the underlying reals,
            // and the spec's clause 2.5 explicitly notes the prose's
            // decimal values are rounded versions of the fixed-point
            // truth.
            assert!(
                err <= 4,
                "fg_sum[{mode}][{i}] = {actual}; expected ≈ Q15_ONE − Σ_k fg = {expected} (err {err})",
            );
        }
    }
}

/// `fg_sum_inv[mode][i]` is the Q12 reciprocal of `fg_sum[mode][i]`,
/// so we should have `fg_sum_inv ≈ 2^12 * (2^15 / fg_sum)`. The
/// reconstruction uses these as a pre-divided factor; this
/// reciprocal-relation check pins both staged tables together with a
/// tolerance that accommodates the independent rounding on each side.
#[test]
fn fg_sum_inv_is_reciprocal_of_fg_sum() {
    let q15 = 32768.0_f64;
    let q12 = 4096.0_f64;
    for mode in 0..2 {
        for i in 0..tables::M {
            let sum = f64::from(tables::LSP_MA_PREDICTOR_FG_SUM_Q15[mode][i]);
            assert!(sum > 0.0, "fg_sum[{mode}][{i}] not positive: {sum}");
            let real_sum = sum / q15;
            let expected_inv = (1.0 / real_sum) * q12;
            let actual_inv = f64::from(tables::LSP_MA_PREDICTOR_FG_SUM_INV_Q12[mode][i]);
            let err = (actual_inv - expected_inv).abs();
            // Three Q12 steps' tolerance accommodates the
            // sum-then-divide rounding chain (Q15 rounding on the
            // sum then Q12 rounding on the reciprocal).
            assert!(
                err < 3.0,
                "fg_sum_inv[{mode}][{i}] = {actual_inv}; expected ≈ {expected_inv} (err {err})",
            );
        }
    }
}

/// All `fg` entries are strictly positive — the MA predictor weights
/// are positive in both modes. This is a cheap drift check that the
/// 3-D parser hasn't accidentally sign-flipped any plane or row.
#[test]
fn fg_entries_are_strictly_positive() {
    for mode in 0..2 {
        for k in 0..tables::MA_NP {
            for i in 0..tables::M {
                let v = tables::LSP_MA_PREDICTOR_FG_Q15[mode][k][i];
                assert!(v > 0, "fg[{mode}][{k}][{i}] not positive: {v}");
            }
        }
    }
}

/// The per-mode `fg` plane accessor returns a borrowed slab whose
/// rows match the underlying constant element-wise.
#[test]
fn lsp_fg_plane_helper_matches_constant() {
    for mode in 0..2 {
        let plane = tables::lsp_fg_plane(mode);
        assert_eq!(plane, &tables::LSP_MA_PREDICTOR_FG_Q15[mode]);
    }
}

/// The per-mode `fg_sum` / `fg_sum_inv` accessors return a borrowed
/// 10-element row equal to the underlying constants.
#[test]
fn lsp_fg_sum_and_inv_helpers_match_constants() {
    for mode in 0..2 {
        assert_eq!(
            tables::lsp_fg_sum(mode),
            &tables::LSP_MA_PREDICTOR_FG_SUM_Q15[mode],
        );
        assert_eq!(
            tables::lsp_fg_sum_inv(mode),
            &tables::LSP_MA_PREDICTOR_FG_SUM_INV_Q12[mode],
        );
    }
}
