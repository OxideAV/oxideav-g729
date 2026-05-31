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
the staged `LSP.BIT` vector lies in the codebook dimensions.

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

Helper spec-dimension constants are also exposed:
`PRM_SIZE = 11`, `BITS_PER_FRAME = 80`, `M = 10` (LP order),
`L_WINDOW = 240`, `GRID_POINTS = 60`, `NC0 = 128`, `NC1 = 32`,
`L0_BITS = 1`, `L1_BITS = 7`, `L2_BITS = L3_BITS = 5`,
`LSP_TOTAL_BITS = 18`.

Round-195 lookup helpers (bounds-checked):

- `lsp_l1_entry(l1: usize) -> &'static [i16; M]` — borrows the
  full 10-coefficient row of the first-stage codebook.
- `lsp_l2_entry(l2: usize) -> &'static [i16]` — borrows the lower
  5 coefficients of the packed second-stage codebook (the L2 split).
- `lsp_l3_entry(l3: usize) -> &'static [i16]` — borrows the upper
  5 coefficients of the same row (the L3 split).

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

## What is NOT wired up

Every decode/encode entry point still returns `Error::NotImplemented`.
The remaining codebook tables (gain GA/GB, MA predictor `fg`,
postfilter interpolation `tab_hup_*`, taming `tab_zone`, Annex B
DTX/CNG, LSF↔LSP cos/slope tables) are staged under
`docs/audio/g729/tables/` but are not yet compiled in; the Implementer
leaves them out until the docs collaborator's specifier pass clarifies
the per-clause wire-up direction.

Round 195 wires the L1/L2 codebooks themselves but not the full
LSP-from-bits reconstruction: that requires the MA predictor `fg`
(switched by `L0`), the §3.2.4 rearrangement steps that enforce a
minimum adjacent distance, and the 4-step stability clamp. Those
follow on once `fg` is compiled. The harness still validates the
on-wire bit layout: the staged `LSP.BIT` test vector (the ITU's
dedicated L0/L1/L2/L3 exerciser) is walked frame-by-frame, each
frame's L1 7-bit index is checked against `NC0`, and each L2/L3
5-bit index against `NC1` — a necessary condition for any future
LSP reconstruction to remain in-bounds.

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
