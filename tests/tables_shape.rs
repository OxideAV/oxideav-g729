//! Shape checks on the bit-exact tables module.
//!
//! These tests confirm the constants emitted by `build.rs` match the
//! sizes documented in the G.729 specification + the spec-role
//! `.meta` sidecars staged at `docs/audio/g729/tables/`. They do
//! NOT compare individual numeric values against any third-party
//! source — only against the spec-published structural invariants
//! (array lengths and a small handful of trivially-checkable
//! properties such as "high-pass filter must have a zero DC gain"
//! that follow directly from the spec's filter description).

use oxideav_g729::tables;

#[test]
fn highpass_100hz_has_three_taps_each() {
    assert_eq!(tables::PREPROC_HP_100HZ_B_Q13.len(), 3);
    assert_eq!(tables::PREPROC_HP_100HZ_A_Q13.len(), 3);
}

#[test]
fn highpass_140hz_has_three_taps_each() {
    assert_eq!(tables::PREPROC_HP_140HZ_B_Q12.len(), 3);
    assert_eq!(tables::PREPROC_HP_140HZ_A_Q12.len(), 3);
}

/// Spec §3.1 / §4.2 — both pre-processing high-pass filters are
/// 2nd-order IIR filters whose feedback path is normalised so that
/// `a[0]` is the implicit unity denominator anchor. The Q13 form
/// stores `1.0` as `8192` and the Q12 form stores `1.0` as `4096`.
#[test]
fn highpass_denominator_normalisation() {
    assert_eq!(tables::PREPROC_HP_100HZ_A_Q13[0], 8192);
    assert_eq!(tables::PREPROC_HP_140HZ_A_Q12[0], 4096);
}

/// Spec §3.1 — the numerator coefficients of a high-pass IIR
/// `b[z] = b0 + b1·z^-1 + b2·z^-2` cancel at DC (z = 1) so the
/// transfer function is exactly zero at 0 Hz. For both filters
/// this means `b0 + b1 + b2 = 0`.
#[test]
fn highpass_numerator_dc_zero() {
    let b100: i32 = tables::PREPROC_HP_100HZ_B_Q13
        .iter()
        .map(|&v| v as i32)
        .sum();
    assert_eq!(b100, 0, "100 Hz HP numerator b0+b1+b2 must equal 0 at DC");
    let b140: i32 = tables::PREPROC_HP_140HZ_B_Q12
        .iter()
        .map(|&v| v as i32)
        .sum();
    assert_eq!(b140, 0, "140 Hz HP numerator b0+b1+b2 must equal 0 at DC");
}

/// Spec §3.1 — the canonical Hp(z) form
/// `Hp(z) = b0 (1 − z^-1)(1 − z^-1) / a(z)` produces equal first
/// and last numerator taps with the middle tap at `−2·b0`. The
/// `.meta` sidecar records the Q13 100 Hz row as `{7699, -15398, 7699}`
/// — a sanity check that the build script did not silently swap the
/// order while transcribing the CSV.
#[test]
fn highpass_100hz_symmetry_and_middle_tap() {
    let b = &tables::PREPROC_HP_100HZ_B_Q13;
    assert_eq!(b[0], b[2], "first and last b-taps must match (HP symmetry)");
    assert_eq!(b[1], -2 * b[0], "middle b-tap must equal −2·b0");
}

/// Same canonical-form check for the 140 Hz variant — also produced
/// by `b0·(1 − z^-1)^2`. The `.meta` sidecar records `{1899, -3798, 1899}`.
#[test]
fn highpass_140hz_symmetry_and_middle_tap() {
    let b = &tables::PREPROC_HP_140HZ_B_Q12;
    assert_eq!(b[0], b[2], "first and last b-taps must match (HP symmetry)");
    assert_eq!(b[1], -2 * b[0], "middle b-tap must equal −2·b0");
}

/// Spec §4.1 Table 8 — bit allocation table CSV carries 13 entries.
/// (The first 11 correspond to the analysis parameters of the main
/// 8 kbit/s coder per the table's `.meta` sidecar; entries 11..13
/// remain raw extracted data for future interpretation by a
/// behavioural round once the clean-room trace doc lands.)
#[test]
fn bit_allocation_table_size() {
    assert_eq!(tables::BIT_ALLOCATION_TABLE8.len(), 13);
}

/// Each documented bit-count entry must be non-negative and fit in
/// a reasonable parameter range (≤ 32 bits per parameter, since the
/// largest analysis parameter in any G.729 profile is well under
/// that). This sanity-checks against a typo or sign flip during
/// extraction.
#[test]
fn bit_allocation_entries_in_reasonable_range() {
    for (i, &v) in tables::BIT_ALLOCATION_TABLE8.iter().enumerate() {
        assert!(
            (0..=32).contains(&v),
            "entry {i} = {v} outside [0, 32] bit-count range"
        );
    }
}

#[test]
fn pow2_table_length_is_33() {
    assert_eq!(tables::BASIC_OP_POW2_TABLE.len(), 33);
}

#[test]
fn log2_table_length_is_33() {
    assert_eq!(tables::BASIC_OP_LOG2_TABLE.len(), 33);
}

#[test]
fn invsqrt_table_length_is_49() {
    assert_eq!(tables::BASIC_OP_INVSQRT_TABLE.len(), 49);
}

/// `Pow2` is monotonically increasing on [0, 1] because `2^x` is.
/// The 33-entry table tabulates that interval, so successive entries
/// must satisfy `tab[i] ≤ tab[i+1]`.
#[test]
fn pow2_table_is_monotonically_non_decreasing() {
    let t = &tables::BASIC_OP_POW2_TABLE;
    for w in t.windows(2) {
        assert!(
            w[0] <= w[1],
            "Pow2 LUT must be non-decreasing; found {} > {}",
            w[0],
            w[1]
        );
    }
}

/// `Log2(1 + x)` is monotonically increasing on x ∈ [0, 1] for the
/// same reason. The first table entry is 0 (log2(1) = 0).
#[test]
fn log2_table_starts_at_zero_and_is_non_decreasing() {
    let t = &tables::BASIC_OP_LOG2_TABLE;
    assert_eq!(t[0], 0, "log2(1) must be 0");
    for w in t.windows(2) {
        assert!(
            w[0] <= w[1],
            "Log2 LUT must be non-decreasing; found {} > {}",
            w[0],
            w[1]
        );
    }
}

/// `1 / sqrt(x)` is monotonically *decreasing* in x, so the
/// `Inv_sqrt` LUT must be non-increasing across the tabulated x
/// range.
#[test]
fn invsqrt_table_is_non_increasing() {
    let t = &tables::BASIC_OP_INVSQRT_TABLE;
    for w in t.windows(2) {
        assert!(
            w[0] >= w[1],
            "Inv_sqrt LUT must be non-increasing; found {} < {}",
            w[0],
            w[1]
        );
    }
}

/// All exposed tables stay within Q15-style fixed-point bounds
/// (i.e. fit in an `i16`, which is enforced at build time too).
#[test]
fn all_tables_fit_in_i16() {
    // Compile-time constants are typed `[i16; N]` — this test just
    // confirms we can iterate them without surprises and that no
    // entry is `i16::MIN` (a sentinel that would indicate parsing
    // collapsed an unexpected token to `-32768`).
    let check = |name: &str, t: &[i16]| {
        for (i, &v) in t.iter().enumerate() {
            assert_ne!(v, i16::MIN, "{name}[{i}] is i16::MIN sentinel");
        }
    };
    check("PREPROC_HP_100HZ_B_Q13", &tables::PREPROC_HP_100HZ_B_Q13);
    check("PREPROC_HP_100HZ_A_Q13", &tables::PREPROC_HP_100HZ_A_Q13);
    check("PREPROC_HP_140HZ_B_Q12", &tables::PREPROC_HP_140HZ_B_Q12);
    check("PREPROC_HP_140HZ_A_Q12", &tables::PREPROC_HP_140HZ_A_Q12);
    check("BIT_ALLOCATION_TABLE8", &tables::BIT_ALLOCATION_TABLE8);
    check("BASIC_OP_POW2_TABLE", &tables::BASIC_OP_POW2_TABLE);
    check("BASIC_OP_LOG2_TABLE", &tables::BASIC_OP_LOG2_TABLE);
    check("BASIC_OP_INVSQRT_TABLE", &tables::BASIC_OP_INVSQRT_TABLE);
}
