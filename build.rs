//! Build script — compiles G.729 bit-exact numeric tables from staged
//! CSVs into `pub const [i16; N]` arrays for inclusion via
//! `src/tables.rs`.
//!
//! ## Provenance
//!
//! The CSVs under `tables/` are byte-for-byte copies of the
//! spec-role-named outputs under
//! `docs/audio/g729/tables/` produced by that directory's
//! `extract.py` against the ITU-T G.729 (06/2012) electronic
//! attachment ZIP (SHA-256
//! `979680ff3b52b13a5701453178efd32b53e340114638fd330fdb6669eec86620`,
//! verified per `docs/audio/g729/tables/README.md`). The Implementer
//! reads only the CSV + `.meta` sidecar pairs; no algorithmic source
//! from the ITU electronic attachment is read here or anywhere in
//! `src/`.
//!
//! ## Why a build script
//!
//! The numeric arrays are pure factual content that the spec defers to
//! the electronic attachment. Compiling them from CSV (rather than
//! retyping into `.rs`) guarantees the constants in the binary are
//! the exact values the spec normatively requires, with no chance of
//! transcription drift between `tables/*.csv` and `src/`. Each emitted
//! `pub const` carries the spec clause + Q-format in its leading doc
//! comment, derived from the paired `.meta` file's `spec_role` field.
//!
//! ## What is emitted
//!
//! For each `(csv, meta)` pair listed in `TABLES`, the script emits a
//! file `$OUT_DIR/<stem>.rs` containing a `pub const` Rust array typed
//! as the meta's `c_type` (Word16 -> i16; only i16 tables are wired up
//! this round). The wrapping `src/tables.rs` module then `include!`s
//! each file behind a public re-export.
//!
//! ## Round-173 scope
//!
//! This round wires up:
//! * §3.1 / §4.2 pre/post high-pass filter coefficients (b100, a100,
//!   b140, a140) — 4 × 3-entry Word16 arrays.
//! * §4.1 spec Table 8 bit-allocation per parameter (bitsno) — 13-entry
//!   Word16 array (the CSV literal count; spec PRM_SIZE constant is 11
//!   but the source array carries two trailing entries — preserved as
//!   extracted).
//! * basic_op() Pow2, Log2, Inv_sqrt lookup tables — 33, 33, 49 entries.
//!
//! Larger codebook tables (LSP L1/L2, gain GA/GB, MA predictor fg,
//! interpolation filters, etc.) are NOT compiled this round; their
//! addition is gated on the docs collaborator handoff (#859 per
//! workspace memory).

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

/// Static description of a table to compile: (file stem, declared
/// element count). The element count is asserted against the CSV's
/// observed literal count; a mismatch fails the build so a CSV that
/// drifts from the documented shape is caught immediately.
struct Table {
    stem: &'static str,
    elements: usize,
}

const TABLES: &[Table] = &[
    // §3.1 / §4.2 pre_proc()+post_pro() — IIR HPF coefs.
    Table {
        stem: "preproc-highpass-100Hz-b-Q13",
        elements: 3,
    },
    Table {
        stem: "preproc-highpass-100Hz-a-Q13",
        elements: 3,
    },
    Table {
        stem: "preproc-highpass-140Hz-b-Q12",
        elements: 3,
    },
    Table {
        stem: "preproc-highpass-140Hz-a-Q12",
        elements: 3,
    },
    // §4.1 spec Table 8 bit allocation per analysis parameter.
    Table {
        stem: "bit-allocation-per-parameter-table8",
        elements: 13,
    },
    // basic_op math lookup tables.
    Table {
        stem: "basic-op-pow2-table",
        elements: 33,
    },
    Table {
        stem: "basic-op-log2-table",
        elements: 33,
    },
    Table {
        stem: "basic-op-invsqrt-table",
        elements: 49,
    },
];

/// Maps a CSV stem to its Rust-side `pub const` identifier. Kept here
/// rather than parsed from the meta so the identifier is reviewable in
/// one place and tied to a specific spec-role naming choice.
fn const_ident(stem: &str) -> &'static str {
    match stem {
        "preproc-highpass-100Hz-b-Q13" => "HPF_PREPROC_100HZ_B_Q13",
        "preproc-highpass-100Hz-a-Q13" => "HPF_PREPROC_100HZ_A_Q13",
        "preproc-highpass-140Hz-b-Q12" => "HPF_PREPROC_140HZ_B_Q12",
        "preproc-highpass-140Hz-a-Q12" => "HPF_PREPROC_140HZ_A_Q12",
        "bit-allocation-per-parameter-table8" => "BIT_ALLOCATION_TABLE8",
        "basic-op-pow2-table" => "POW2_TABLE_Q15",
        "basic-op-log2-table" => "LOG2_TABLE_Q15",
        "basic-op-invsqrt-table" => "INV_SQRT_TABLE_Q15",
        _ => panic!("unknown table stem: {stem}"),
    }
}

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let tables_dir = manifest_dir.join("tables");
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    for table in TABLES {
        let csv = tables_dir.join(format!("{}.csv", table.stem));
        let meta = tables_dir.join(format!("{}.meta", table.stem));
        let dest = out_dir.join(format!("{}.rs", sanitize_stem(table.stem)));

        println!("cargo:rerun-if-changed={}", csv.display());
        println!("cargo:rerun-if-changed={}", meta.display());

        emit_table(table, &csv, &meta, &dest);
    }
}

/// Reads a CSV (one Word16 literal per line) plus its meta sidecar and
/// writes a `pub const NAME: [i16; N] = [...];` declaration to `dest`,
/// carrying the meta's spec_role + source provenance in a doc comment.
fn emit_table(table: &Table, csv: &Path, meta: &Path, dest: &Path) {
    let csv_text =
        fs::read_to_string(csv).unwrap_or_else(|e| panic!("read {} failed: {e}", csv.display()));
    let values: Vec<i16> = csv_text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.trim()
                .parse::<i32>()
                .unwrap_or_else(|e| panic!("parse `{}` from {} failed: {e}", l, csv.display()))
        })
        .map(|v| {
            i16::try_from(v).unwrap_or_else(|_| {
                panic!(
                    "value {v} in {} does not fit Word16 (i16) per meta's c_type",
                    csv.display()
                )
            })
        })
        .collect();

    assert_eq!(
        values.len(),
        table.elements,
        "{} CSV element count ({}) does not match TABLES declaration ({})",
        csv.display(),
        values.len(),
        table.elements,
    );

    let meta_text =
        fs::read_to_string(meta).unwrap_or_else(|e| panic!("read {} failed: {e}", meta.display()));
    let spec_role = meta_field(&meta_text, "spec_role").unwrap_or("(spec_role missing)");
    let source_file = meta_field(&meta_text, "source_file").unwrap_or("(source_file missing)");
    let source_sha256 = meta_field(&meta_text, "source_sha256").unwrap_or("(sha missing)");
    let zip_sha256 = meta_field(&meta_text, "source_zip_sha256").unwrap_or("(zip sha missing)");
    let c_identifier = meta_field(&meta_text, "c_identifier").unwrap_or("(c_identifier missing)");

    let ident = const_ident(table.stem);
    let mut body = String::new();
    body.push_str(&format!("/// {spec_role}\n"));
    body.push_str("///\n");
    body.push_str(&format!(
        "/// Compiled from `tables/{}.csv` (staged copy of `docs/audio/g729/tables/{}.csv`).\n",
        table.stem, table.stem
    ));
    body.push_str(&format!(
        "/// Original ITU C identifier: `{c_identifier}`.\n"
    ));
    body.push_str(&format!(
        "/// Source file inside ITU electronic attachment: `{source_file}`.\n"
    ));
    body.push_str(&format!("/// Source file SHA-256: `{source_sha256}`.\n"));
    body.push_str(&format!(
        "/// Electronic-attachment ZIP SHA-256: `{zip_sha256}`.\n"
    ));
    body.push_str(&format!(
        "pub const {ident}: [i16; {len}] = [\n",
        ident = ident,
        len = values.len(),
    ));
    for v in &values {
        body.push_str(&format!("    {v},\n"));
    }
    body.push_str("];\n");

    fs::write(dest, body).unwrap_or_else(|e| panic!("write {} failed: {e}", dest.display()));
}

/// Extracts a `key: value` line from a `.meta` text body. Returns
/// `None` if the key is not present. The match is whole-line and
/// trims the surrounding whitespace.
fn meta_field<'a>(meta_text: &'a str, key: &str) -> Option<&'a str> {
    for line in meta_text.lines() {
        if let Some(rest) = line.strip_prefix(&format!("{key}:")) {
            return Some(rest.trim());
        }
    }
    None
}

/// Converts a CSV stem (with dashes) into the `<stem>.rs` filename used
/// in `$OUT_DIR`. Dashes are kept; the `include!` site references the
/// same stem so this is just a pass-through with a single source of
/// truth.
fn sanitize_stem(stem: &str) -> String {
    stem.to_string()
}
