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
}
