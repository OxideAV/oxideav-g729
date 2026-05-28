//! Build script — parses the spec-role-named CSV tables under
//! `tables/` and emits them as compile-time `pub const` arrays in
//! `$OUT_DIR/tables_data.rs`. The Rust source includes the result
//! once via `include!`; nothing in `src/` retypes any numeric value.
//!
//! ## Provenance chain
//!
//! Each `tables/*.csv` is a verbatim copy of the canonical clean-room
//! workspace at `docs/audio/g729/tables/` (see that directory's
//! `README.md` for the provenance documentation). The paired `.meta`
//! sidecar records spec role + clause for the table.
//!
//! ## Clean-room guarantee
//!
//! No ITU reference C source is consulted, mirrored, or paraphrased
//! by this script. The CSVs were produced upstream by
//! `docs/audio/g729/tables/extract.py`, which reads only the data-
//! bearing `tab_*.c` files (numeric arrays + their inline comments)
//! from the ITU electronic-attachment package and emits one CSV per
//! top-level array, *named for the table's role in the G.729
//! specification* rather than its internal C identifier. The script
//! does not read any algorithmic source (`coder.c`, `decoder.c`,
//! `lpc.c`, …). This build script reads only the CSVs.
//!
//! ## Tables emitted in this scope
//!
//! All under `oxideav_g729::tables`, types and constants are
//! generated from the CSV headers and inferred Q-formats:
//!
//! * `PREPROC_HP_100HZ_B_Q13` / `PREPROC_HP_100HZ_A_Q13`   — §3.1 / §4.2
//! * `PREPROC_HP_140HZ_B_Q12` / `PREPROC_HP_140HZ_A_Q12`   — §3.1 alt.
//! * `BIT_ALLOCATION_TABLE8`                               — §4.1 Table 8
//! * `BASIC_OP_POW2_TABLE` / `BASIC_OP_LOG2_TABLE`         — basic_op LUTs
//! * `BASIC_OP_INVSQRT_TABLE`                              — basic_op LUT
//!
//! All values fit in `i16` (the C type is `Word16`, i.e. a 16-bit
//! signed integer).

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

/// One CSV-to-const binding. Each entry pins the on-disk file name,
/// the generated Rust identifier, and a one-line doc string drawn
/// from the spec-role description in the table's `.meta` sidecar.
struct CsvBinding {
    csv_name: &'static str,
    rust_ident: &'static str,
    doc: &'static str,
}

const BINDINGS: &[CsvBinding] = &[
    CsvBinding {
        csv_name: "preproc-highpass-100Hz-b-Q13.csv",
        rust_ident: "PREPROC_HP_100HZ_B_Q13",
        doc: "G.729 §3.1 / §4.2 pre_proc / post_pro — IIR high-pass filter \
              numerator (b) coefficients, fc = 100 Hz, Q13.",
    },
    CsvBinding {
        csv_name: "preproc-highpass-100Hz-a-Q13.csv",
        rust_ident: "PREPROC_HP_100HZ_A_Q13",
        doc: "G.729 §3.1 / §4.2 pre_proc / post_pro — IIR high-pass filter \
              denominator (a) coefficients, fc = 100 Hz, Q13.",
    },
    CsvBinding {
        csv_name: "preproc-highpass-140Hz-b-Q12.csv",
        rust_ident: "PREPROC_HP_140HZ_B_Q12",
        doc: "G.729 §3.1 alternate IIR high-pass filter numerator (b) \
              coefficients, fc = 140 Hz, Q12.",
    },
    CsvBinding {
        csv_name: "preproc-highpass-140Hz-a-Q12.csv",
        rust_ident: "PREPROC_HP_140HZ_A_Q12",
        doc: "G.729 §3.1 alternate IIR high-pass filter denominator (a) \
              coefficients, fc = 140 Hz, Q12.",
    },
    CsvBinding {
        csv_name: "bit-allocation-per-parameter-table8.csv",
        rust_ident: "BIT_ALLOCATION_TABLE8",
        doc: "G.729 §4.1 bit allocation per analysis parameter, raw \
              extract corresponding to spec Table 8 \
              (Bit allocation of 8 kbit/s ITU-T G.729). The CSV \
              carries 13 entries; consumers select the leading slice \
              that matches the parameter count required by the \
              chosen profile.",
    },
    CsvBinding {
        csv_name: "basic-op-pow2-table.csv",
        rust_ident: "BASIC_OP_POW2_TABLE",
        doc: "G.729 basic-op Pow2 lookup table (33 entries, Q15-style \
              fractional fixed-point). Used to evaluate `2^x` for \
              fractional `x ∈ [0, 1]` via piece-wise linear \
              interpolation between adjacent table entries.",
    },
    CsvBinding {
        csv_name: "basic-op-log2-table.csv",
        rust_ident: "BASIC_OP_LOG2_TABLE",
        doc: "G.729 basic-op Log2 lookup table (33 entries). Used to \
              evaluate `log2(1 + x)` for normalized inputs via \
              piece-wise linear interpolation between adjacent table \
              entries.",
    },
    CsvBinding {
        csv_name: "basic-op-invsqrt-table.csv",
        rust_ident: "BASIC_OP_INVSQRT_TABLE",
        doc: "G.729 basic-op Inv_sqrt lookup table (49 entries). Used \
              to evaluate `1 / sqrt(x)` for normalized inputs via \
              piece-wise linear interpolation between adjacent table \
              entries.",
    },
];

fn main() {
    let manifest_dir =
        PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR set by cargo"));
    let tables_dir = manifest_dir.join("tables");
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR set by cargo"));
    let out_path = out_dir.join("tables_data.rs");

    let mut emitted = String::new();
    emitted.push_str(
        "// Generated by build.rs from `tables/*.csv`. Do not edit by hand.\n\
         // See `build.rs` for the binding list and provenance chain.\n\n",
    );

    for b in BINDINGS {
        let csv_path = tables_dir.join(b.csv_name);
        println!("cargo:rerun-if-changed=tables/{}", b.csv_name);
        let raw = fs::read_to_string(&csv_path)
            .unwrap_or_else(|e| panic!("read {}: {e}", csv_path.display()));
        let values = parse_csv_i16(&raw, &csv_path);

        emitted.push_str("/// ");
        emitted.push_str(b.doc);
        emitted.push('\n');
        emitted.push_str(&format!(
            "///\n/// Source CSV: `tables/{}` (verbatim from \
             `docs/audio/g729/tables/`).\n",
            b.csv_name
        ));
        emitted.push_str(&format!(
            "pub const {ident}: [i16; {n}] = [\n",
            ident = b.rust_ident,
            n = values.len()
        ));
        for v in &values {
            emitted.push_str(&format!("    {v},\n"));
        }
        emitted.push_str("];\n\n");
    }

    fs::write(&out_path, emitted).expect("write generated tables_data.rs");
    println!("cargo:rerun-if-changed=build.rs");
}

/// Parses a CSV that holds one row per line, each row containing one
/// or more decimal integers separated by commas. The flattened
/// sequence of integers is returned. Empty lines and trailing
/// whitespace are tolerated; everything else must parse as an `i16`.
fn parse_csv_i16(raw: &str, path: &Path) -> Vec<i16> {
    let mut out = Vec::new();
    for (lineno, line) in raw.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        for tok in trimmed.split(',') {
            let tok = tok.trim();
            if tok.is_empty() {
                continue;
            }
            match tok.parse::<i32>() {
                Ok(v) if (i16::MIN as i32..=i16::MAX as i32).contains(&v) => {
                    out.push(v as i16);
                }
                Ok(v) => panic!(
                    "{}:{}: value {v} out of i16 range",
                    path.display(),
                    lineno + 1
                ),
                Err(e) => panic!(
                    "{}:{}: failed to parse `{tok}`: {e}",
                    path.display(),
                    lineno + 1
                ),
            }
        }
    }
    out
}
