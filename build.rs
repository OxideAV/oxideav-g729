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
//! ## Cumulative scope
//!
//! Round 173 wired up the foundational HPF coefficients, the §4.1
//! Table-8 bit-allocation array, and the three `basic_op` math
//! lookups. Round 189 extends that with the §3.2.1 LP-analysis tables
//! (Hamming window + autocorrelation lag-window pair), the §3.2.5
//! LSF root-search cosine grid, the §3.7 pitch interpolation filters
//! (`inter_3` / `inter_3l`), and the §3.9 MA gain-prediction
//! coefficient vector (`pred`).
//!
//! Tables wired this build:
//!
//! * §3.1 / §4.2 pre/post high-pass filter coefficients (b100, a100,
//!   b140, a140) — 4 × 3-entry Word16 arrays. (r173)
//! * §4.1 spec Table 8 bit-allocation per parameter (bitsno) — 13-entry
//!   Word16 array (the CSV literal count; spec PRM_SIZE constant is 11
//!   but the source array carries two trailing entries — preserved as
//!   extracted). (r173)
//! * basic_op() Pow2, Log2, Inv_sqrt lookup tables — 33, 33, 49
//!   entries. (r173)
//! * §3.2.1 LP-analysis windowing: `hamwindow` (240 Q15 samples),
//!   `lag_h` / `lag_l` (10 + 10 Q15 entries forming the
//!   double-precision 60 Hz bandwidth-expansion lag window pair).
//!   (r189)
//! * §3.2.5 `az_lsf()` root-search grid: 61-entry Q15 cosine table
//!   over the half-circle. (r189)
//! * §3.7 pitch interpolation filters: `inter_3` (13-tap analysis)
//!   and `inter_3l` (31-tap synthesis), both Q15. (r189)
//! * §3.9 MA gain-prediction coefficients (`pred`) — 4 Q13 entries
//!   {0.68, 0.58, 0.34, 0.19}. (r189)
//!
//! Round 195 adds the §3.2.4 LSP quantiser two-stage VQ codebooks:
//!
//! * `LSP_QUANT_CODEBOOK_L1_Q13` — first-stage 10-D codebook
//!   (`lspcb1`), 128 × 10 Word16 (Q13); 7 transmitted bits per frame.
//! * `LSP_QUANT_CODEBOOK_L2_Q13` — second-stage 10-wide table
//!   (`lspcb2`) holding the two 5-D split codebooks side-by-side,
//!   32 × 10 Word16 (Q13); 5 + 5 transmitted bits per frame.
//!
//! Each codebook entry is one **row** of the source CSV — the build
//! script's [`Shape::Matrix`] path emits a 2-D `[[i16; cols]; rows]`
//! array so callers can index by `(stage_index, coefficient_index)`
//! directly.
//!
//! Round 201 adds the §3.2.4 LSP MA-predictor `fg` matrix that
//! completes the L0-switched reconstruction inputs (spec eq (20),
//! eq (20a)):
//!
//! * `LSP_MA_PREDICTOR_FG_Q15` — the [2][MA_NP=4][M=10] Word16 (Q15)
//!   coefficient cube. `L0` selects the outer slice; the inner
//!   4-tap row is the spec's 4th-order MA over past quantised LSP
//!   residuals.
//! * `LSP_MA_PREDICTOR_FG_SUM_Q15` — the [2][M] Q15 per-mode column
//!   sums (`Σ fg[mode][k][i]` over k) used by the reconstruction
//!   equation's denominator-style term.
//! * `LSP_MA_PREDICTOR_FG_SUM_INV_Q12` — the [2][M] Q12 reciprocals
//!   of the per-mode column sums, pre-tabulated for the
//!   reconstruction.
//!
//! `Shape::Cube { planes, rows, cols }` is added to support the 3-D
//! `fg` layout: the CSV stores rows in row-major order (`planes ×
//! rows` lines of `cols` comma-separated literals); the emitted
//! Rust array shape is `[[[i16; cols]; rows]; planes]`.
//!
//! Remaining codebook tables (gain GA/GB, postfilter interpolation
//! (`tab_hup_*`), taming (`tab_zone`), Annex B DTX/CNG,
//! LSF↔LSP cos/slope tables) are NOT compiled yet; their addition
//! is gated on the docs collaborator handoff (#859 per workspace
//! memory).

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

/// Static description of a table to compile.
///
/// `shape` selects the on-disk CSV layout AND the emitted Rust array
/// shape:
///
/// * [`Shape::Flat`] — one Word16 literal per line. Emits
///   `pub const NAME: [i16; N] = [...]`.
/// * [`Shape::Matrix { rows, cols }`] — `rows` lines, each carrying
///   `cols` comma-separated Word16 literals. Emits
///   `pub const NAME: [[i16; cols]; rows] = [[..], ..]`.
/// * [`Shape::Cube { planes, rows, cols }`] — `planes × rows` lines
///   each carrying `cols` comma-separated Word16 literals. The CSV
///   stores planes in row-major order (plane 0's rows first, then
///   plane 1's rows, etc.). Emits
///   `pub const NAME: [[[i16; cols]; rows]; planes] = [...]`.
///
/// In every case the build script asserts the observed literal counts
/// against the declared dimensions, so any CSV-vs-declaration drift
/// trips the build immediately.
struct Table {
    stem: &'static str,
    shape: Shape,
}

#[derive(Clone, Copy)]
enum Shape {
    /// 1-D table: `elements` literals on consecutive lines.
    Flat { elements: usize },
    /// 2-D codebook: `rows` rows × `cols` columns of literals (one row
    /// per line, `,`-separated columns).
    Matrix { rows: usize, cols: usize },
    /// 3-D coefficient cube: `planes × rows` lines × `cols` columns.
    /// Plane `p` occupies CSV lines `p * rows .. (p + 1) * rows`.
    Cube {
        planes: usize,
        rows: usize,
        cols: usize,
    },
}

const TABLES: &[Table] = &[
    // §3.1 / §4.2 pre_proc()+post_pro() — IIR HPF coefs.
    Table {
        stem: "preproc-highpass-100Hz-b-Q13",
        shape: Shape::Flat { elements: 3 },
    },
    Table {
        stem: "preproc-highpass-100Hz-a-Q13",
        shape: Shape::Flat { elements: 3 },
    },
    Table {
        stem: "preproc-highpass-140Hz-b-Q12",
        shape: Shape::Flat { elements: 3 },
    },
    Table {
        stem: "preproc-highpass-140Hz-a-Q12",
        shape: Shape::Flat { elements: 3 },
    },
    // §4.1 spec Table 8 bit allocation per analysis parameter.
    Table {
        stem: "bit-allocation-per-parameter-table8",
        shape: Shape::Flat { elements: 13 },
    },
    // basic_op math lookup tables.
    Table {
        stem: "basic-op-pow2-table",
        shape: Shape::Flat { elements: 33 },
    },
    Table {
        stem: "basic-op-log2-table",
        shape: Shape::Flat { elements: 33 },
    },
    Table {
        stem: "basic-op-invsqrt-table",
        shape: Shape::Flat { elements: 49 },
    },
    // §3.2.1 LP-analysis windowing.
    Table {
        stem: "lpc-hamming-window-Q15",
        shape: Shape::Flat { elements: 240 },
    },
    Table {
        stem: "lpc-autocorr-lag-window-high-Q15",
        shape: Shape::Flat { elements: 10 },
    },
    Table {
        stem: "lpc-autocorr-lag-window-low-Q15",
        shape: Shape::Flat { elements: 10 },
    },
    // §3.2.5 az_lsf() root-search grid.
    Table {
        stem: "lsf-search-grid-cos-Q15",
        shape: Shape::Flat { elements: 61 },
    },
    // §3.7 pitch interpolation filters.
    Table {
        stem: "pitch-interpolation-filter-analysis-Q15",
        shape: Shape::Flat { elements: 13 },
    },
    Table {
        stem: "pitch-interpolation-filter-synthesis-Q15",
        shape: Shape::Flat { elements: 31 },
    },
    // §3.9 MA gain-prediction coefficients.
    Table {
        stem: "gain-quantizer-ma-predictor-Q13",
        shape: Shape::Flat { elements: 4 },
    },
    // §3.2.4 LSP quantiser two-stage VQ codebooks. L1 is the
    // first-stage 10-D codebook (NC0 = 128 entries, 7 bits); L2 is the
    // second-stage 10-wide table that holds the two 5-D split codebooks
    // (NC1 = 32 entries, 5 + 5 bits) packed side-by-side per the
    // staged trace doc §3.5 (single `lspcb2` array in the ITU source).
    Table {
        stem: "lsp-quantizer-codebook-L1-Q13",
        shape: Shape::Matrix {
            rows: 128,
            cols: 10,
        },
    },
    Table {
        stem: "lsp-quantizer-codebook-L2-Q13",
        shape: Shape::Matrix { rows: 32, cols: 10 },
    },
    // §3.2.4 LSP MA-predictor `fg` family (round 201). `fg` is the
    // 3-D switched-predictor coefficient cube: outer dim selects the
    // L0 predictor mode (2 modes), middle dim is the MA history depth
    // (MA_NP = 4), inner dim is the LP order (M = 10). `fg_sum` and
    // `fg_sum_inv` are the per-mode column sum and its reciprocal,
    // pre-tabulated as the reconstruction inputs.
    Table {
        stem: "lsp-ma-predictor-fg-Q15",
        shape: Shape::Cube {
            planes: 2,
            rows: 4,
            cols: 10,
        },
    },
    Table {
        stem: "lsp-ma-predictor-fg-sum-Q15",
        shape: Shape::Matrix { rows: 2, cols: 10 },
    },
    Table {
        stem: "lsp-ma-predictor-fg-sum-inv-Q12",
        shape: Shape::Matrix { rows: 2, cols: 10 },
    },
    // §3.9.2 gain-quantizer conjugate-structure codebooks (round 231).
    // GA is the 3-bit first-stage 2-D codebook (8 entries × 2 components:
    // adaptive-codebook-gain contribution in Q14, fixed-codebook-gain
    // correction contribution in Q12). GB is the 4-bit second-stage 2-D
    // codebook (16 entries × 2 components, same per-column Q-formats).
    // Reconstruction per spec eqs (73)/(74) sums the per-row components
    // of the two codebooks indexed by GA / GB respectively.
    Table {
        stem: "gain-quantizer-codebook-GA-Q14-Q12",
        shape: Shape::Matrix { rows: 8, cols: 2 },
    },
    Table {
        stem: "gain-quantizer-codebook-GB-Q14-Q12",
        shape: Shape::Matrix { rows: 16, cols: 2 },
    },
    // §3.9.3 transmitted-index permutation + inverse permutation +
    // partial-search thresholds (the spec's robustness mapping that
    // reorders GA / GB indices before transmission so a single-bit
    // channel error lands on a perceptually-close codebook entry).
    Table {
        stem: "gain-quantizer-codebook-GA-permutation",
        shape: Shape::Flat { elements: 8 },
    },
    Table {
        stem: "gain-quantizer-codebook-GA-inverse-permutation",
        shape: Shape::Flat { elements: 8 },
    },
    Table {
        stem: "gain-quantizer-codebook-GA-thresholds-Q14",
        shape: Shape::Flat { elements: 4 },
    },
    Table {
        stem: "gain-quantizer-codebook-GB-permutation",
        shape: Shape::Flat { elements: 16 },
    },
    Table {
        stem: "gain-quantizer-codebook-GB-inverse-permutation",
        shape: Shape::Flat { elements: 16 },
    },
    Table {
        stem: "gain-quantizer-codebook-GB-thresholds-Q15",
        shape: Shape::Flat { elements: 8 },
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
        "lpc-hamming-window-Q15" => "LPC_HAMMING_WINDOW_Q15",
        "lpc-autocorr-lag-window-high-Q15" => "LPC_LAG_WINDOW_HIGH_Q15",
        "lpc-autocorr-lag-window-low-Q15" => "LPC_LAG_WINDOW_LOW_Q15",
        "lsf-search-grid-cos-Q15" => "LSF_SEARCH_GRID_COS_Q15",
        "pitch-interpolation-filter-analysis-Q15" => "PITCH_INTERP_FILTER_ANALYSIS_Q15",
        "pitch-interpolation-filter-synthesis-Q15" => "PITCH_INTERP_FILTER_SYNTHESIS_Q15",
        "gain-quantizer-ma-predictor-Q13" => "GAIN_QUANT_MA_PREDICTOR_Q13",
        "lsp-quantizer-codebook-L1-Q13" => "LSP_QUANT_CODEBOOK_L1_Q13",
        "lsp-quantizer-codebook-L2-Q13" => "LSP_QUANT_CODEBOOK_L2_Q13",
        "lsp-ma-predictor-fg-Q15" => "LSP_MA_PREDICTOR_FG_Q15",
        "lsp-ma-predictor-fg-sum-Q15" => "LSP_MA_PREDICTOR_FG_SUM_Q15",
        "lsp-ma-predictor-fg-sum-inv-Q12" => "LSP_MA_PREDICTOR_FG_SUM_INV_Q12",
        "gain-quantizer-codebook-GA-Q14-Q12" => "GAIN_QUANT_CODEBOOK_GA_Q14_Q12",
        "gain-quantizer-codebook-GB-Q14-Q12" => "GAIN_QUANT_CODEBOOK_GB_Q14_Q12",
        "gain-quantizer-codebook-GA-permutation" => "GAIN_QUANT_GA_PERMUTATION",
        "gain-quantizer-codebook-GA-inverse-permutation" => "GAIN_QUANT_GA_INVERSE_PERMUTATION",
        "gain-quantizer-codebook-GA-thresholds-Q14" => "GAIN_QUANT_GA_THRESHOLDS_Q14",
        "gain-quantizer-codebook-GB-permutation" => "GAIN_QUANT_GB_PERMUTATION",
        "gain-quantizer-codebook-GB-inverse-permutation" => "GAIN_QUANT_GB_INVERSE_PERMUTATION",
        "gain-quantizer-codebook-GB-thresholds-Q15" => "GAIN_QUANT_GB_THRESHOLDS_Q15",
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

/// Reads the CSV + meta sidecar pair for `table` and writes the
/// corresponding `pub const` declaration to `dest`. The emitted Rust
/// type follows [`Table::shape`] — 1-D arrays for [`Shape::Flat`],
/// 2-D arrays for [`Shape::Matrix`] (one inner array per row of the
/// CSV).
fn emit_table(table: &Table, csv: &Path, meta: &Path, dest: &Path) {
    let csv_text =
        fs::read_to_string(csv).unwrap_or_else(|e| panic!("read {} failed: {e}", csv.display()));
    let meta_text =
        fs::read_to_string(meta).unwrap_or_else(|e| panic!("read {} failed: {e}", meta.display()));
    let spec_role = meta_field(&meta_text, "spec_role").unwrap_or("(spec_role missing)");
    let zip_sha256 = meta_field(&meta_text, "source_zip_sha256").unwrap_or("(zip sha missing)");

    let ident = const_ident(table.stem);

    let mut body = String::new();
    body.push_str(&format!("/// {spec_role}\n"));
    body.push_str("///\n");
    body.push_str(&format!(
        "/// Compiled from `tables/{}.csv` (staged copy of `docs/audio/g729/tables/{}.csv`).\n",
        table.stem, table.stem
    ));
    body.push_str(
        "/// Per-file provenance (origin filename, byte ranges, file-level SHA-256) is\n",
    );
    body.push_str(
        "/// recorded in the paired `.meta` sidecar under `docs/audio/g729/tables/`; the\n",
    );
    body.push_str(
        "/// crate's `src/` deliberately does not name any algorithmic-source file of the\n",
    );
    body.push_str("/// staged electronic attachment.\n");
    body.push_str(&format!(
        "/// Electronic-attachment ZIP SHA-256: `{zip_sha256}`.\n"
    ));

    match table.shape {
        Shape::Flat { elements } => {
            let values = parse_flat_csv(&csv_text, csv);
            assert_eq!(
                values.len(),
                elements,
                "{} CSV element count ({}) does not match TABLES declaration ({})",
                csv.display(),
                values.len(),
                elements,
            );
            body.push_str(&format!(
                "pub const {ident}: [i16; {len}] = [\n",
                ident = ident,
                len = values.len(),
            ));
            for v in &values {
                body.push_str(&format!("    {v},\n"));
            }
            body.push_str("];\n");
        }
        Shape::Matrix { rows, cols } => {
            let matrix = parse_matrix_csv(&csv_text, csv, rows, cols);
            body.push_str(&format!(
                "pub const {ident}: [[i16; {cols}]; {rows}] = [\n",
                ident = ident,
                rows = rows,
                cols = cols,
            ));
            for row in &matrix {
                body.push_str("    [");
                for (i, v) in row.iter().enumerate() {
                    if i > 0 {
                        body.push_str(", ");
                    }
                    body.push_str(&format!("{v}"));
                }
                body.push_str("],\n");
            }
            body.push_str("];\n");
        }
        Shape::Cube { planes, rows, cols } => {
            let cube = parse_cube_csv(&csv_text, csv, planes, rows, cols);
            body.push_str(&format!(
                "pub const {ident}: [[[i16; {cols}]; {rows}]; {planes}] = [\n",
                ident = ident,
                planes = planes,
                rows = rows,
                cols = cols,
            ));
            for plane in &cube {
                body.push_str("    [\n");
                for row in plane {
                    body.push_str("        [");
                    for (i, v) in row.iter().enumerate() {
                        if i > 0 {
                            body.push_str(", ");
                        }
                        body.push_str(&format!("{v}"));
                    }
                    body.push_str("],\n");
                }
                body.push_str("    ],\n");
            }
            body.push_str("];\n");
        }
    }

    fs::write(dest, body).unwrap_or_else(|e| panic!("write {} failed: {e}", dest.display()));
}

/// Parses a 1-D CSV (one literal per non-empty line) into a `Vec<i16>`,
/// panicking with the CSV path on parse / overflow.
fn parse_flat_csv(csv_text: &str, csv: &Path) -> Vec<i16> {
    csv_text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| parse_word16(l.trim(), csv))
        .collect()
}

/// Parses a 2-D CSV: one row per non-empty line, `,`-separated literals
/// per row. Asserts the row count and per-row column count exactly
/// match the declared shape.
fn parse_matrix_csv(csv_text: &str, csv: &Path, rows: usize, cols: usize) -> Vec<Vec<i16>> {
    let parsed: Vec<Vec<i16>> = csv_text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let row: Vec<i16> = l
                .split(',')
                .map(|tok| parse_word16(tok.trim(), csv))
                .collect();
            assert_eq!(
                row.len(),
                cols,
                "{} row literal count ({}) does not match declared cols ({})",
                csv.display(),
                row.len(),
                cols,
            );
            row
        })
        .collect();
    assert_eq!(
        parsed.len(),
        rows,
        "{} row count ({}) does not match TABLES declaration ({})",
        csv.display(),
        parsed.len(),
        rows,
    );
    parsed
}

/// Parses a 3-D CSV: `planes` lines, each carrying
/// `rows * cols` comma-separated literals stored in row-major order
/// (row 0's cols, then row 1's cols, …). Line count and per-line
/// literal count are both asserted against the declared shape; any
/// drift trips the build with the offending stem in the error.
fn parse_cube_csv(
    csv_text: &str,
    csv: &Path,
    planes: usize,
    rows: usize,
    cols: usize,
) -> Vec<Vec<Vec<i16>>> {
    let flat_rows: Vec<Vec<i16>> = csv_text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let row: Vec<i16> = l
                .split(',')
                .map(|tok| parse_word16(tok.trim(), csv))
                .collect();
            assert_eq!(
                row.len(),
                rows * cols,
                "{} per-plane literal count ({}) does not match rows * cols ({} * {} = {})",
                csv.display(),
                row.len(),
                rows,
                cols,
                rows * cols,
            );
            row
        })
        .collect();
    assert_eq!(
        flat_rows.len(),
        planes,
        "{} plane count ({}) does not match declared planes ({})",
        csv.display(),
        flat_rows.len(),
        planes,
    );
    flat_rows
        .into_iter()
        .map(|plane_flat| {
            (0..rows)
                .map(|r| plane_flat[r * cols..(r + 1) * cols].to_vec())
                .collect()
        })
        .collect()
}

/// Parses one Word16 literal token (the input must already be
/// whitespace-trimmed). Panics with the CSV path on overflow.
fn parse_word16(tok: &str, csv: &Path) -> i16 {
    let v = tok
        .parse::<i32>()
        .unwrap_or_else(|e| panic!("parse `{tok}` from {} failed: {e}", csv.display()));
    i16::try_from(v).unwrap_or_else(|_| {
        panic!(
            "value {v} in {} does not fit Word16 (i16) per meta's c_type",
            csv.display()
        )
    })
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
