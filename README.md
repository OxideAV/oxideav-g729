# oxideav-g729

Pure-Rust ITU-T **G.729** (CS-ACELP, 8 kbit/s) narrowband speech codec.

> **Status: clean-room rebuild in progress.**
>
> The previous implementation was retired under the OxideAV workspace
> clean-room policy: several modules transcribed numeric tables verbatim
> from an external reference-software distribution and documented
> matching its behaviour by citing specific source files of it. The
> clean-room policy forbids consulting any external implementation's
> source for any reason, so the provenance of those tables could not
> be defended and master was force-erased to a register-only scaffold
> at `676fc46`. Re-growth is one spec-cited unit at a time.

## What is wired up

Round 173 landed the bit-exact numeric-tables foundation; round 189
extends it with the LP-analysis windowing, LSF cosine grid, pitch
interpolation filters, and the MA gain-prediction coefficients;
round 191 adds the ITU serial bitstream parser and a structural
harness validating it against the staged conformance corpus; round
195 wires in the §3.2.4 LSP-quantiser two-stage VQ codebooks (L1
128×10, L2/L3 packed 32×10), with bounds-checked lookup helpers
and per-frame conformance validation that every L1/L2/L3 index in
the staged `LSP.BIT` vector lies in the codebook dimensions; round
201 completes the §3.2.4 LSP-reconstruction inputs by wiring the
MA-predictor `fg` family (the 2×4×10 Q15 coefficient cube plus
its per-mode `fg_sum` Q15 / `fg_sum_inv` Q12 helpers) and adds a
`Shape::Cube` table type to `build.rs` so 3-D coefficient slabs
emit as `[[[i16; cols]; rows]; planes]` arrays with full
dimensional drift-checking.

All numeric values are compiled at build time by `build.rs` from CSVs
under `tables/`, themselves byte-for-byte copies of the spec-role-named
outputs under `docs/audio/g729/tables/`. Each `pub const` carries the
spec clause, the original ITU C identifier, and the source-file
SHA-256 in its generated doc comment.

### Round 191 — ITU serial bitstream parser

The on-wire `.bit` format used by the ITU conformance test sequences
under `docs/audio/g729/conformance/` is parsed by a new `serial`
module:

- `serial::SYNC_WORD` (`0x6B21`), `serial::BITS_HEADER` (= 80),
  `serial::BIT_ZERO` (`0x007F`), `serial::BIT_ONE` (`0x0081`), and
  `serial::BIT_ERASED` (`0x0000`) — the per-word framing literals
  observed in the staged corpus.
- `serial::FRAME_WORDS = 82`, `serial::FRAME_BYTES = 164` — the
  fixed Word16 / byte cadence per 10 ms / 80-bit frame.
- `serial::parse_frame(&[u8]) -> Result<FrameKind, SerialError>` —
  validates sync + header, distinguishes normal frames
  (`FrameKind::Active([bool; 80])`) from frame-erasure sentinels
  (`FrameKind::Erased`), and rejects mid-frame mixes of normal /
  erased payload words.
- `serial::frame_count(&[u8])` — byte-length sanity check that the
  buffer is a whole number of 164-byte frames.

A new `tests/serial_conformance.rs` integration test walks the
staged `docs/audio/g729/conformance/{g729-core,g729a}/` directories
when present, validating that every `.bit` file is a sequence of
well-formed ITU serial frames, that the frame count matches the
companion `.pst` PCM stream, and that the `.in` PCM stream's
length sits in the expected encoder-look-ahead window. Erasure-
sentinel counts are pinned exactly: `ERASURE.BIT` has 60 of 300
erased frames; `OVERFLOW.BIT` has 1 of 384; `PARITY.BIT` has 0;
all other sequences have 0. The harness logs a clean skip when the
corpus path is absent (published-crate build mode), so the test
surface remains green either way.

| spec clause | constant | shape | role |
|---|---|---|---|
| §3.1 / §4.2 | `HPF_PREPROC_100HZ_B_Q13` | `[i16; 3]` | pre/post-processing 100 Hz HPF b-coefficients (Q13) |
| §3.1 / §4.2 | `HPF_PREPROC_100HZ_A_Q13` | `[i16; 3]` | pre/post-processing 100 Hz HPF a-coefficients (Q13) |
| §3.1 | `HPF_PREPROC_140HZ_B_Q12` | `[i16; 3]` | alternate 140 Hz HPF b-coefficients (Q12) |
| §3.1 | `HPF_PREPROC_140HZ_A_Q12` | `[i16; 3]` | alternate 140 Hz HPF a-coefficients (Q12) |
| §4.1 (Table 8) | `BIT_ALLOCATION_TABLE8` | `[i16; 13]` | per-parameter bit allocation (`bitsno`) |
| `basic_op` | `POW2_TABLE_Q15` | `[i16; 33]` | `Pow2()` lookup |
| `basic_op` | `LOG2_TABLE_Q15` | `[i16; 33]` | `Log2()` lookup |
| `basic_op` | `INV_SQRT_TABLE_Q15` | `[i16; 49]` | `Inv_sqrt()` lookup |
| §3.2.1 | `LPC_HAMMING_WINDOW_Q15` | `[i16; 240]` | LP-analysis Hamming-cos window (`hamwindow`) |
| §3.2.1 | `LPC_LAG_WINDOW_HIGH_Q15` | `[i16; 10]` | 60 Hz bandwidth-expansion lag window — high half of double-precision pair (`lag_h`) |
| §3.2.1 | `LPC_LAG_WINDOW_LOW_Q15` | `[i16; 10]` | 60 Hz bandwidth-expansion lag window — low half (`lag_l`) |
| §3.2.5 | `LSF_SEARCH_GRID_COS_Q15` | `[i16; 61]` | `az_lsf()` root-search cosine grid (`grid`) |
| §3.7 | `PITCH_INTERP_FILTER_ANALYSIS_Q15` | `[i16; 13]` | adaptive-codebook 1/3-resolution analysis filter (`inter_3`) |
| §3.7 | `PITCH_INTERP_FILTER_SYNTHESIS_Q15` | `[i16; 31]` | adaptive-codebook 1/3-resolution synthesis filter (`inter_3l`) |
| §3.9 | `GAIN_QUANT_MA_PREDICTOR_Q13` | `[i16; 4]` | MA gain-prediction coefficients `pred` (≈ {0.68, 0.58, 0.34, 0.19}) |
| §3.2.4 | `LSP_QUANT_CODEBOOK_L1_Q13` | `[[i16; 10]; 128]` | first-stage LSP VQ codebook `lspcb1` (7-bit `L1` index, Q13) |
| §3.2.4 | `LSP_QUANT_CODEBOOK_L2_Q13` | `[[i16; 10]; 32]` | second-stage packed split-VQ codebook `lspcb2` (5-bit `L2` lower / `L3` upper, Q13) |
| §3.2.4 | `LSP_MA_PREDICTOR_FG_Q15` | `[[[i16; 10]; 4]; 2]` | LSP MA-predictor cube `fg` — outer dim is `L0` predictor mode, middle dim is MA history (`MA_NP = 4`), inner dim is LP order `M = 10`, Q15 |
| §3.2.4 | `LSP_MA_PREDICTOR_FG_SUM_Q15` | `[[i16; 10]; 2]` | per-mode `(Q15_ONE − Σ_k fg[mode][k][i])` factor, Q15 |
| §3.2.4 | `LSP_MA_PREDICTOR_FG_SUM_INV_Q12` | `[[i16; 10]; 2]` | per-mode reciprocal of `fg_sum`, Q12 — pre-tabulated reconstruction inputs (spec eq (20a)) |

Helper spec-dimension constants are also exposed:
`PRM_SIZE = 11`, `BITS_PER_FRAME = 80`, `M = 10` (LP order),
`L_WINDOW = 240`, `GRID_POINTS = 60`, `NC0 = 128`, `NC1 = 32`,
`MA_NP = 4`, `L0_BITS = 1`, `L1_BITS = 7`,
`L2_BITS = L3_BITS = 5`, `LSP_TOTAL_BITS = 18`.

Round-195 lookup helpers (bounds-checked):

- `lsp_l1_entry(l1: usize) -> &'static [i16; M]` — borrows the
  full 10-coefficient row of the first-stage codebook.
- `lsp_l2_entry(l2: usize) -> &'static [i16]` — borrows the lower
  5 coefficients of the packed second-stage codebook (the L2 split).
- `lsp_l3_entry(l3: usize) -> &'static [i16]` — borrows the upper
  5 coefficients of the same row (the L3 split).

Round-201 lookup helpers (bounds-checked):

- `lsp_fg_plane(mode: usize) -> &'static [[i16; M]; MA_NP]` —
  borrows the per-mode `[MA_NP][M]` Q15 predictor plane. `mode` is
  the 1-bit `L0` field.
- `lsp_fg_sum(mode: usize) -> &'static [i16; M]` — borrows the
  per-mode `(Q15_ONE − Σ_k fg[mode][k][i])` Q15 factor row.
- `lsp_fg_sum_inv(mode: usize) -> &'static [i16; M]` — borrows the
  per-mode Q12 reciprocal of the `fg_sum` row.

### Round 195 — LSP-quantiser two-stage VQ codebooks

The §3.2.4 first-stage (`lspcb1`, NC0 = 128 × M = 10) and packed
second-stage (`lspcb2`, NC1 = 32 × M = 10) Q13 codebooks compile to
2-D `[[i16; 10]; N]` arrays so callers index by `(stage_index,
coefficient_index)` directly. The build script's CSV reader gains a
matrix path that asserts both row count and per-row column count
against the declared shape, so any CSV drift trips the build with
the offending stem in the error.

The companion `tests/serial_conformance.rs` walks the staged
`LSP.BIT` vector (the ITU sequence whose `READMETV.txt` self-
documents as the "LSP quantization (the L0/L1/L2/L3 VQ)" exerciser)
in both `g729-core/` and `g729a/`, extracts (L0, L1, L2, L3) from
the first 18 bits of every active frame MSB-first per spec Table 8
NOTE, and asserts each index lies in the codebook dimensions. A
companion synthetic-frame test pins the MSB-first packing
convention so the bit ordering is locked even when the corpus is
absent (published-crate mode).

The companion `tests/tables_shape.rs` pins each round-173 table as
before, and additionally verifies the round-189 tables' structural
properties:

- `LPC_HAMMING_WINDOW_Q15.len() == L_WINDOW`; peak value is Q15 ≈ 1.0
  (`32767`); every entry is strictly positive.
- `LPC_LAG_WINDOW_HIGH_Q15.len() == LPC_LAG_WINDOW_LOW_Q15.len() == M`;
  `lag_h` is strictly monotonically decreasing.
- `LSF_SEARCH_GRID_COS_Q15.len() == GRID_POINTS + 1`; endpoints match
  the CSV literals (`32760`, `-32760`); midpoint is exactly `0`;
  the full grid is strictly decreasing and antisymmetric about the
  midpoint (`grid[i] == -grid[N-1-i]`).
- Both pitch interpolation filters: the maximum tap is positive and
  equals the maximum-magnitude tap (centre dominates).
- `GAIN_QUANT_MA_PREDICTOR_Q13` matches `[5571, 4751, 2785, 1556]`
  exactly and round-trips to the spec-documented real values
  {0.68, 0.58, 0.34, 0.19} within one Q13 quantisation step.

### Round 201 — LSP MA-predictor `fg` family

The §3.2.4 MA-predictor reconstruction inputs land this round. The
3-D `fg` cube is wired alongside its per-mode `fg_sum` Q15 factor
and the `fg_sum_inv` Q12 reciprocal — these are the staged tables
that, with the round-195 L1+L2/L3 codebooks, let the reconstruction
equation (spec eq (20) / (20a)) be evaluated end-to-end.

A new `Shape::Cube { planes, rows, cols }` `build.rs` table type
parses the 3-D CSV layout (`planes` lines, each carrying
`rows * cols` comma-separated literals in row-major order within
each plane) and asserts both the line count and per-line literal
count against the declared shape. The emitted Rust type is
`pub const NAME: [[[i16; cols]; rows]; planes]`, so callers index
by `(mode, history_step, coordinate)` directly.

The companion `tests/tables_shape.rs` grows with 8 new round-201
tests:

- shape (`2 × MA_NP × M` for `fg`; `2 × M` for `fg_sum` and
  `fg_sum_inv`) with the outer dim cross-checked against
  `1 << L0_BITS == 2`;
- first row of each `fg` plane pinned to CSV literals (a focused
  drift check on the 3-D cube reader's row-major flattening);
- history-depth peak decay (`fg[mode][MA_NP - 1]` peak magnitude
  is strictly less than the `fg[mode][0]` peak magnitude in both
  modes);
- strict positivity across all 80 `fg` entries (sign-flip drift
  check);
- `fg_sum` matches the spec-stated `(Q15_ONE − Σ_k fg[mode][k][i])`
  factor within 4 Q15 ulps;
- `fg_sum_inv` is the Q12 reciprocal of `fg_sum` within 3 Q12 ulps;
- `lsp_fg_plane` / `lsp_fg_sum` / `lsp_fg_sum_inv` helpers each
  return slices equal to the underlying constants.

`build.rs` also drops the per-table doc comment that previously
named the source file inside the staged electronic attachment. The
provenance chain itself is unchanged (it still lives in the
`.meta` sidecars under `docs/audio/g729/tables/`); only the per-
constant rustdoc emission is scrubbed to keep the in-`src/` doc
surface free of algorithmic-source filenames.

## What is NOT wired up

Every decode/encode entry point still returns `Error::NotImplemented`.
The remaining codebook tables (gain GA/GB, postfilter interpolation
`tab_hup_*`, taming `tab_zone`, Annex B DTX/CNG, LSF↔LSP cos/slope
tables) are staged under `docs/audio/g729/tables/` but are not yet
compiled in; the Implementer leaves them out until the docs
collaborator's specifier pass clarifies the per-clause wire-up
direction.

With round 201 the §3.2.4 reconstruction inputs (`fg`, `fg_sum`,
`fg_sum_inv`) are present, but the full LSP-from-bits
reconstruction is not yet implemented in code: that still requires
the §3.2.4 rearrangement steps that enforce a minimum adjacent
distance (twice, `J = 0.0012` then `J = 0.0006`) and the 4-step
stability clamp. The reconstruction function itself follows next
round.

The harness still validates the on-wire bit layout: the staged
`LSP.BIT` test vector (the ITU's dedicated L0/L1/L2/L3 exerciser)
is walked frame-by-frame, each frame's L1 7-bit index is checked
against `NC0`, and each L2/L3 5-bit index against `NC1` — a
necessary condition for any future LSP reconstruction to remain
in-bounds.

The 80 transmitted bits per frame that `serial::parse_frame` returns
are otherwise still an **opaque payload** at this layer; mapping the
non-LSP bits onto the §4.1 Table-8 parameters (P0/P1/P2 pitch,
C1/S1/GA1/GB1 fixed-codebook indices, etc.) is a future round.

## Clean-room provenance

| Step | Artefact | Verification |
|---|---|---|
| 1 | ITU recommendation page | <https://www.itu.int/rec/T-REC-G.729-201206-I/en> (manual fetch by docs collaborator; not consulted at runtime) |
| 2 | ITU electronic-attachment ZIP | SHA-256 `979680ff3b52b13a5701453178efd32b53e340114638fd330fdb6669eec86620` |
| 3 | `tab_*.c` files inside ZIP | each `.meta` records the source file's individual SHA-256 + the line range |
| 4 | `docs/audio/g729/tables/<spec-name>.csv` | extractor (`extract.py`) is deterministic; re-running against the same ZIP reproduces byte-identical CSVs |
| 5 | `crates/oxideav-g729/tables/<spec-name>.csv` (this crate) | byte-for-byte copy of step-4 output, included with the published crate so a from-crates.io build does not depend on the workspace `docs/` |
| 6 | `pub const` arrays in `oxideav_g729::tables` | `build.rs` parses the CSV literal-per-line, asserts shape vs the meta sidecar, and emits a `pub const [i16; N] = [...]` declaration into `$OUT_DIR` |

No algorithmic source from the ITU electronic attachment is read by
`build.rs`, `src/`, or `tests/`. The extractor under
`docs/audio/g729/tables/` reads only `tab_*.c` data files and `LD8K.H`
for symbolic dimensions; it does not read `coder.c`, `decoder.c`,
`lpc.c`, `acelp_co.c`, or any other algorithmic source.

## License

MIT — see [`LICENSE`](LICENSE).
