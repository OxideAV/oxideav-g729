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
interpolation filters, and the MA gain-prediction coefficients. All
values are compiled at build time by `build.rs` from CSVs under
`tables/`, themselves byte-for-byte copies of the spec-role-named
outputs under `docs/audio/g729/tables/`. Each `pub const` carries
the spec clause, the original ITU C identifier, and the source-file
SHA-256 in its generated doc comment.

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

Helper spec-dimension constants are also exposed:
`PRM_SIZE = 11`, `BITS_PER_FRAME = 80`, `M = 10` (LP order),
`L_WINDOW = 240`, `GRID_POINTS = 60`.

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
The remaining codebook tables (LSP L1/L2, gain GA/GB, MA predictor
`fg`, postfilter interpolation `tab_hup_*`, taming `tab_zone`,
Annex B DTX/CNG, LSF↔LSP cos/slope tables) are staged under
`docs/audio/g729/tables/` but are not yet compiled in; the Implementer
leaves them out until the docs collaborator's specifier pass clarifies
the per-clause wire-up direction.

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
