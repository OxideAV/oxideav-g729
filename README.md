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
dimensional drift-checking; round 207 builds the §3.2.4
reconstruction algorithm itself on top of those tables — a new
`lsp_reconstruct` module that maps `(L0, L1, L2, L3)` to a
reconstructed `ω̂^(m)` LSF vector via codebook sum (eq (19)),
twice-applied rearrangement (`J = 0.0012` then `J = 0.0006`),
MA-prediction (eq (20)), and the 4-step stability clamp (floor
0.005, min-gap 0.0391, ceil 3.135); round 213 chains the §3.2.5
per-subframe LSP interpolation (spec eq (24), `q_i^(1) = ½(q^prev
+ q^curr)`, `q_i^(2) = q^curr` — *linear interpolation in the
cosine domain*) onto the §3.2.4 output, with `omega_to_q` /
`q_to_omega` boundary helpers and a stateful `LspInterpolator`
that advances `previous_q ← current_q` per frame; round 218
wires the §3.2.6 LSP→LP conversion (spec eqs (13) / (14)
sum/difference polynomial recursion `F_1 = Π(1 − 2q_{odd}·z^-1
+ z^-2)` / `F_2 = Π(1 − 2q_{even}·z^-1 + z^-2)`, eq (25) `(1 ±
z^-1)` factor restoration, eq (26) `A(z) = (F'_1 + F'_2) / 2`
recombination with `F'_1` symmetric / `F'_2` antisymmetric of
length 6) with a stateless `lsp_to_lp(&[f32; 10]) -> [f32; 10]`
per-subframe entry point and an `lsf_to_lp` boundary wrapper;
round 225 wires the §4.1 / Table-8 parameter unpacker, splitting
the round-191 serial 80-bit payload into the 15 typed codeword
indices (`L0` / `L1` / `L2` / `L3`, `P1` / `P0`, `C1` / `S1` /
`GA1` / `GB1`, `P2`, `C2` / `S2` / `GA2` / `GB2`) the §4.1
decode procedure consumes, and ties the §3.7.2 pitch-parity
predicate against the corpus — 0 mismatches on every active
frame of `SPEECH.BIT` (3 750 × 2 variants), non-zero mismatches
on `PARITY.BIT` (the dedicated §4.1.2 concealment exerciser),
every codeword in-domain on every active frame of the full
staged corpus; round 231 wires the §3.9.2 / §4.1.5
conjugate-structure gain-VQ decode-side reconstruction in a new
`gain_reconstruct` module that maps each transmitted `(GA, GB)`
index pair into the quantised `(ĝ_p, γ̂)` pair via spec
eqs (73) / (74) (`ĝ_p = GA[GA][0] + GB[GB][0]` in Q14
column 0, `γ̂ = GA[GA][1] + GB[GB][1]` in Q12 column 1),
with the GA (8 × 2 Q14/Q12), GB (16 × 2 Q14/Q12), per-stage
permutation, inverse permutation, and partial-search-threshold
tables all landing alongside in `tables`, plus an integration
test that runs the per-frame reconstruction across every
active frame of the staged conformance corpus and pins
`ĝ_p ∈ [0, 2]` and `γ̂ ∈ [0, 11]` end-to-end; round 239
chains the §3.9.1 4th-order MA gain prediction stage onto
the round-231 `γ̂` output via a new stateful `gain_predict`
module that owns the 4-tap history `[Û^(m-1), Û^(m-2),
Û^(m-3), Û^(m-4)]` initialised per spec Table 9 / §4.3 to
`(-14, -14, -14, -14)` dB, computes the codevector energy
`E` (eq (66), `10·log10((1/40)·Σ_{n=0..39} c(n)^2)`), the
predicted energy `Ẽ^(m)` (eq (69), the 4-tap MA dot product
`Σ_{i=1..=4} b_i · Û^(m-i)` with `[b1 b2 b3 b4] = [0.68
0.58 0.34 0.19]` reading from the staged
`GAIN_QUANT_MA_PREDICTOR_Q13` Q13 table), the predicted
gain `g'_c` (eq (71), `10^((Ẽ + Ē - E)/20)` with
`Ē = 30 dB`), and finally the quantised fixed-codebook
gain `ĝ_c = γ̂ · g'_c` (eq (65)), then advances the history
per eq (72) decode form `Û^(m) = 20·log10(γ̂)`; round 249
ties the round-225 transmitted-index unpacker to the round-231
codebook lookup by way of the §3.9.3 codeword-mapping layer —
a new `gain_index_map` module exposes `demap_ga` / `demap_gb`
(the `imap1` / `imap2` inverse permutations that turn on-wire
indices into codebook-domain indices), `map_ga` / `map_gb`
(the symmetric encoder-side forward permutations), and
`demap_frame` (the per-frame wrapper that demaps all four
gain indices in one call); the existing
`gain_reconstruct::reconstruct_frame_gains` is wired through
this layer so the decode chain is spec-conformant from the
on-wire bits, and a new
`gain_reconstruct::reconstruct_gains_from_transmitted` helper
exposes the same demap-then-reconstruct path as a per-pair
primitive; round 255 wires the §4.1.3 pitch-delay decode in a
new `pitch_decode` module that maps the transmitted `(P1, P2)`
indices into per-subframe fractional pitch delays via
`decode_t1_from_p1` (spec image `f0027-01.jpg` — `if P1 < 197:
int(T1) = (P1+2)/3 + 19`, `frac = P1 − 3·int(T1) + 58`, else
`int(T1) = P1 − 112`, `frac = 0`), `derive_t_min` (spec image
`f0027-02.jpg` — the `±5 / [20, 143]` subframe-2 search-window
derivation), `decode_t2_from_p2` (spec image `f0027-03.jpg` —
`int(T2) = (P2+2)/3 − 1 + t_min`, `frac = P2 − 2 − 3·((P2+2)/3
− 1)`), and the per-frame `decode_frame` wrapper that chains
them in spec §4.1.3 order, with the symmetric `encode_p1` /
`encode_p2` forward mappings (spec eqs (41) / (42)) exposed for
encoder wire-up and the corpus-walking round-trip test;
round 266 wires the §3.8 / §4.1.4 fixed (algebraic) codebook
decode in a new `fixed_codebook` module that maps the
transmitted `(C, S)` codewords into the four `(position,
sign)` pulses of spec eq (45) — `decode_positions` (spec
eq (62) — the 13-bit `C` field splits as 3+3+3+4 with track 3
carrying `2·(m_3/5) + jx` in 4 bits; track residues per spec
Table 7 are `(0, 1, 2, 3+jx) mod 5`), `decode_signs` (spec
eq (61) — `S = s_0 + 2·s_1 + 4·s_2 + 8·s_3` with `s_k = 1`
⇔ positive), `build_codevector` (spec eq (45) — four signed
unit impulses on the 40-sample subframe), `decode_pulses`
(per-subframe wrapper), and `decode_frame` (per-frame wrapper
threading `(C1, S1)` and `(C2, S2)` into their subframes),
with the symmetric `encode_positions` / `encode_signs`
forward mappings exposed for encoder wire-up and the
corpus-walking round-trip test against `FIXED.BIT`.

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
| §3.9.2 | `GAIN_QUANT_CODEBOOK_GA_Q14_Q12` | `[[i16; 2]; 8]` | first-stage gain VQ codebook `gbk1` (3-bit GA index; col 0 = `g_p` contribution Q14, col 1 = `γ` contribution Q12) |
| §3.9.2 | `GAIN_QUANT_CODEBOOK_GB_Q14_Q12` | `[[i16; 2]; 16]` | second-stage gain VQ codebook `gbk2` (4-bit GB index; same column convention) |
| §3.9.3 | `GAIN_QUANT_GA_PERMUTATION` / `..._INVERSE_PERMUTATION` | `[i16; 8]` each | transmission-side robustness mapping for GA |
| §3.9.3 | `GAIN_QUANT_GB_PERMUTATION` / `..._INVERSE_PERMUTATION` | `[i16; 16]` each | transmission-side robustness mapping for GB |
| §3.9.2 | `GAIN_QUANT_GA_THRESHOLDS_Q14` | `[i16; 4]` | encoder-side partial-search thresholds for GA (`NCODE1 − NCAN1`) |
| §3.9.2 | `GAIN_QUANT_GB_THRESHOLDS_Q15` | `[i16; 8]` | encoder-side partial-search thresholds for GB (`NCODE2 − NCAN2`) |

Helper spec-dimension constants are also exposed:
`PRM_SIZE = 11`, `BITS_PER_FRAME = 80`, `M = 10` (LP order),
`L_WINDOW = 240`, `GRID_POINTS = 60`, `NC0 = 128`, `NC1 = 32`,
`MA_NP = 4`, `L0_BITS = 1`, `L1_BITS = 7`,
`L2_BITS = L3_BITS = 5`, `LSP_TOTAL_BITS = 18`,
`NCODE1 = 8`, `NCODE2 = 16`, `GAIN_VQ_DIM = 2`,
`GAIN_VQ_COL_GP = 0`, `GAIN_VQ_COL_GC = 1`.

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

Round-231 lookup helpers (bounds-checked):

- `gain_ga_entry(ga: usize) -> &'static [i16; GAIN_VQ_DIM]` —
  borrows the 2-element `(g_p Q14, γ Q12)` row of the first-stage
  codebook at the supplied GA index.
- `gain_gb_entry(gb: usize) -> &'static [i16; GAIN_VQ_DIM]` —
  borrows the same shape row from the second-stage codebook.

Round-239 spec-constant surface:

- `gain_predict::FIXED_CODEBOOK_MEAN_ENERGY_DB = 30.0` — `Ē`
  per spec §3.9.1 (the mean energy of the fixed-codebook
  excitation used in eqs (67) / (68) / (71)).
- `gain_predict::GAIN_PREDICTOR_INIT_DB = -14.0` — `Û^(k)`
  start-up value per spec Table 9 / §4.3.
- `gain_predict::CODEVECTOR_LEN = 40` — spec §3.8 / eq (66)
  codevector length over which the energy is averaged.

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

### Round 207 — §3.2.4 LSP-frame reconstruction

A new `oxideav_g729::lsp_reconstruct` module ties the round-195
codebooks and the round-201 MA-predictor cube into the spec
§3.2.4 decode pipeline:

- `codebook_sum(l1, l2, l3) -> Result<[f32; 10], _>` evaluates spec
  eq (19) (`l̂_i = L1_i(L1) + L2_i(L2)` for `i ∈ 1..=5`,
  `L1_i(L1) + L3_{i-5}(L3)` for `i ∈ 6..=10`). Q13 codebook entries
  are converted to `f32` at the boundary; out-of-range indices
  surface as typed `LspReconstructError` variants (no panic).
- `rearrange_pass(coefs, j)` implements the spec §3.2.4 figure
  `F0013-01` fix-up; `rearrange_twice(coefs)` runs it with
  `J = REARRANGE_J1 = 0.0012` and then `J = REARRANGE_J2 = 0.0006`.
- `stability_clamp(coefs)` applies the spec §3.2.4 4-step
  invariants (sort ascending, `CLAMP_FLOOR = 0.005`,
  `CLAMP_MIN_GAP = 0.0391`, `CLAMP_CEIL = 3.135`).
- `LspReconstructor::new()` initialises the 4-frame MA history to
  the spec start-up vector `l̂_i = i · π / 11`.
  `reconstruct_frame(l0, l1, l2, l3)` runs eq (19) → rearrange
  twice → eq (20) MA prediction → stability clamp, advances the
  internal history, and returns the reconstructed `ω̂^(m)` LSF
  vector.

12 unit tests pin the algorithmic invariants (spec start-up
vector, codebook-sum boundary, every typed error variant,
rearrangement min-distance + no-op, two-pass `J2 = 0.0006`
finish, stability-clamp floor + gap + ceil + no-op, end-to-end
clamp compliance, MA-history shift, both `L0` modes).

### Round 218 — §3.2.6 LSP-to-LP conversion

A new `oxideav_g729::lsp_to_lp` module ties the round-213
per-subframe cosine-domain LSP output into the 10-coefficient LP
synthesis filter `A(z) = 1 + Σ_{i=1..=10} a_i · z^-i` that the
spec §3.5+ impulse-response / synthesis stages consume:

- `lsp_to_lp(q_in: &[f32; 10]) -> [f32; 10]` runs the three
  spec-cited steps of §3.2.6 in sequence: (1) build the
  `F_1` / `F_2` sum/difference polynomial coefficients `f_1(i)`,
  `f_2(i)` for `i ∈ 0..=5` via the spec's recursion derived from
  polynomial multiplication by `(1 − 2·q·z^-1 + z^-2)` (with the
  convention `f(-1) = 0`), (2) restore the `(1 ± z^-1)` factors
  via eq (25) (`f'_1(i) = f_1(i) + f_1(i-1)`,
  `f'_2(i) = f_2(i) − f_2(i-1)`), (3) recombine via eq (26)
  (`a_i = ½·f'_1(i) + ½·f'_2(i)` for `i ∈ 1..=5` and
  `a_i = ½·f'_1(11-i) − ½·f'_2(11-i)` for `i ∈ 6..=10`, using
  `F'_1` symmetric / `F'_2` antisymmetric of length 6).
- `lsf_to_lp(omega: &[f32; 10]) -> [f32; 10]` — convenience
  wrapper that converts LSF-domain `ω̂` to LP via the cosine
  boundary `q_i = cos(ω̂_i)` (so callers staying in the LSF domain
  don't have to drive `omega_to_q` explicitly).
- `LpCoefficients = [f32; 10]` — public alias for the output type;
  slot `i - 1` holds `a_i` for `i ∈ 1..=10` with `a_0 = 1.0`
  implicit (not stored).

8 new unit tests pin the algorithmic invariants:

- start-up state (`ω̂_i = i · π / 11`) produces finite
  coefficients (no NaN / no inf from the recursion);
- the recursion matches a brute-force polynomial-multiplication
  oracle (literally multiplying out `Π(1 − 2·q·z^-1 + z^-2)` by
  vector accumulation, then applying eqs (25) / (26) in the test)
  on two distinct LSP patterns to a ≤ 1e-4 drift — locks the
  in-place inner-loop ordering against any read-after-write bug;
- closed-form spot checks at `z = 1` and `z = −1`: from eq (26)
  `A(1) = F_1(1) = Π_{odd}(2 − 2·q_i)` (because `F'_2(1) = 0`)
  and `A(−1) = F_2(−1) = Π_{even}(2 + 2·q_i)` (because
  `F'_1(−1) = 0`), each pinned to a ≤ 1e-3 drift — a sign error
  in either half of eq (26) trips one or the other;
- coefficient range stays inside ±32 (defence against a missing
  ½ factor in eq (26) — real-speech `a_i` is O(1) to O(3));
- the `lsf_to_lp` convenience wrapper matches the explicit
  `lsp_to_lp(&omega_to_q(&omega))` pipeline to ≤ 1e-7;
- end-to-end through the §3.2.4 reconstructor + §3.2.5
  interpolator + §3.2.6 conversion on a 3-frame non-steady-state
  `(L0, L1, L2, L3)` chain produces finite per-subframe `a_i`
  vectors AND the two subframes within a frame differ on a
  changing input (the cosine-domain midpoint of `(previous,
  current)` LSPs in subframe 1 vs the bare `current` LSPs in
  subframe 2 yields distinct LP filters);
- the all-zero `q` corner case (all LSPs at `ω̂_i = π/2`)
  produces a finite vector — the well-conditioned input is a
  drift check on the eq (15) recursion's `f(-1) = 0` boundary.

### Round 225 — §4.1 / Table-8 parameter unpacker

The 80-bit payload that `serial::parse_frame` returns as a
`FrameKind::Active` bit array is split into its 15 spec-Table-8
codeword indices by a new `oxideav_g729::parameters` module:

- `Parameters` — `Copy` struct carrying the per-frame indices
  (`l0` / `l1` / `l2` / `l3` for the §3.2.4 LSP quantiser;
  `p1` / `p0` for §3.7 / §3.7.2 subframe-1 pitch + parity;
  `c1` / `s1` for §3.8 subframe-1 fixed codebook;
  `ga1` / `gb1` for §3.9.2 subframe-1 conjugate-structure gain
  VQ; same set with the `2` suffix for subframe 2). Field
  widths: 1+7+5+5 + 8+1+13+4+3+4 + 5+13+4+3+4 = 18 + 33 + 29 =
  **80** bits, matching the spec frame budget.
- `unpack_parameters(&FrameKind) -> Result<Parameters,
  ParameterError>` — the §4.1 frame-level entry point;
  rejects `FrameKind::Erased` with `ParameterError::Erased`
  (the §4.4 concealment path applies for an erasure-sentinel
  frame and consumes no transmitted bits).
- `unpack_bit_array(&[bool; 80]) -> Parameters` — lower-level
  variant for unit-testing the unpacker without spinning the
  framing layer.
- `Parameters::pitch_parity_ok(&self) -> bool` — §3.7.2 /
  §4.1.2 parity check, with **the parity-init value pinned to
  1 (odd-parity convention)** based on the corpus: every active
  frame of `g729-core/SPEECH.BIT` + `g729a/SPEECH.BIT` (3 750 ×
  2 = 7 500 frames of clean encoder output) has
  `P0 = 1 XOR XOR_reduce(six_MSBs(P1))`; under the alternative
  even-parity reading every frame fails. See the *Spec gap*
  note on the method docstring.

The per-codeword bit layout follows the spec **Table-8 NOTE**
("the bit stream ordering is reflected by the order in the
table; for each parameter, the most significant bit (MSB) is
transmitted first"). The 15 start-offsets are derived at
compile time from the codeword-width array and statically
asserted to sum to `BITS_PER_FRAME` so the layout can never
silently drift.

The published `parameters::C_BITS`, `S_BITS`, `GA_BITS`,
`GB_BITS`, `P0_BITS`, `P1_BITS`, `P2_BITS` constants pin the
per-codeword widths at the crate's public surface; their
aggregate constants `FIXED_CODEBOOK_BITS_PER_FRAME = 34`,
`GAIN_QUANT_BITS_PER_FRAME = 14`, `PITCH_BITS_PER_FRAME = 14`
plus the existing `LSP_TOTAL_BITS = 18` re-express the
frame-level grouping decoders rely on (`18 + 34 + 14 + 14 =
80`).

9 new unit tests pin the algorithmic invariants:

- all-zero bits map every codeword to 0; all-ones bits saturate
  every codeword to `(1 << width) - 1`;
- L1 saturation lies in `0..NC0` and L2 / L3 saturation lies in
  `0..NC1` (the codebook dimensions wire up cleanly);
- single-bit flip of each of the 80 array slots changes
  exactly one codeword, and that codeword is the one whose
  `[start, start + width)` window contains the slot — locks
  the MSB-first / Table-8-top-to-bottom convention against
  any off-by-one;
- the documented start offsets `(0, 1, 8, 13, 18, 26, 27, 40,
  44, 47, 51, 56, 69, 73, 76)` hold exactly, and slot 80 is
  the end-of-frame boundary (defensive against a future drift
  in the bit-width table);
- round-trip pack-then-unpack on a hand-chosen `Parameters`
  vector recovers every field bit-exactly;
- §3.7.2 parity-rule worked checks on three crafted P1 values
  pin the odd-parity-init convention;
- `unpack_parameters` rejects an erasure sentinel with
  `ParameterError::Erased`;
- the high-level and low-level entry points agree on the same
  bit array.

2 new integration tests against the staged corpus:

- `unpack_parameters_in_domain_on_full_corpus` walks every
  `.BIT` file in `g729-core/` + `g729a/` and asserts that
  every active frame's `Parameters` has every field in its
  spec-stated domain (L1 < NC0, L2/L3 < NC1, C1/C2 < 2^13,
  signs < 2^4, GA < 2^3, GB < 2^4, P2 < 2^5);
- `pitch_parity_distribution_matches_corpus_intent` asserts
  that `SPEECH.BIT` produces **zero** parity mismatches
  (clean encoder output) and that `PARITY.BIT` produces a
  **non-zero** number of mismatches (the dedicated §4.1.2
  concealment-path exerciser).

With round 225 the small piece of glue between the §3.2.4
reconstructor + §3.2.5 interpolator + §3.2.6 LSP→LP chain
and the on-wire bit stream is in place: `unpack_parameters`
yields the `(L0, L1, L2, L3)` tuple that
`LspReconstructor::reconstruct_frame` consumes directly.

### Round 231 — §3.9.2 conjugate-structure gain-VQ reconstruction

A new `oxideav_g729::gain_reconstruct` module ties the round-225
parameter unpacker output into the spec §3.9.2 / §4.1.5 decode-side
gain reconstruction:

- `reconstruct_gains(ga, gb) -> Result<QuantisedGains, _>` evaluates
  spec eqs (73) / (74): `ĝ_p = GA[GA][0] + GB[GB][0]` (column 0,
  Q14 in both stages) and `γ̂ = GA[GA][1] + GB[GB][1]` (column 1,
  Q12 in both stages); the summation runs in `i32` per Q-format and
  converts to `f32` at the boundary. Out-of-range indices surface
  as typed `GainReconstructError` variants — `GaOutOfRange { index }`
  / `GbOutOfRange { index }` — rather than panicking. The
  transmitted GA / GB codewords (3 + 4 bits) cannot trigger these
  in a well-formed frame.
- `reconstruct_frame_gains(&Parameters)` — per-frame wrapper that
  threads `(GA1, GB1)` into subframe 1 and `(GA2, GB2)` into
  subframe 2, matching the spec §4.1.5 ordering.
- `QuantisedGains` — `Copy` struct carrying `g_p_hat` (the
  quantised adaptive-codebook gain `ĝ_p`, which directly scales
  the adaptive-codebook vector in eq (75)) and `gamma_hat` (the
  quantised fixed-codebook gain correction factor `γ̂`). The
  actual quantised fixed-codebook gain `ĝ_c = γ̂ · g'_c` is left
  for a follow-up round — `g'_c` is the §3.9.1 4th-order MA
  prediction output and is stateful, so wiring it cleanly needs
  its own round.

10 new unit tests pin the algorithmic invariants: per-row CSV
literal at (0, 0) (`GA[0] = (1, 1516)`, `GB[0] = (826, 2005)` →
`ĝ_p = (1 + 826) / 2^14`, `γ̂ = (1516 + 2005) / 2^12`); per-row
literal match at `(NCODE1 - 1, NCODE2 - 1)`; both out-of-range
variants on each codebook; out-of-range-first-wins error
precedence (GA error reported when both indices are out of range);
every (GA, GB) pair in the 8 × 16 domain yields finite gains;
`ĝ_p` lies in `[0, 2]` and `γ̂` lies in `[0, 11]` across every
pair (worst-case row pairs reach ≈10.12) — a Q14 / Q12 divisor
swap would push results well outside this window; per-column
delta isolation (varying GA at fixed GB moves `ĝ_p` by the
Q14-scaled column-0 delta and `γ̂` by the Q12-scaled column-1
delta independently); hand-picked pair `(5, 11)` matches the
algebra; `reconstruct_frame_gains` correctly threads the per-
subframe indices into the right codebooks; codebook row width
matches the published `GAIN_VQ_DIM = 2` constant.

7 new structural tests in `tests/tables_shape.rs` for the staged
tables (`NCODE1 × GAIN_VQ_DIM` for GA, `NCODE2 × GAIN_VQ_DIM` for
GB, row counts cross-checked against `1 << GA_BITS == NCODE1` and
`1 << GB_BITS == NCODE2`; first-row CSV-literal pins; column-
constant convention; both permutations are complete covers of
`0..NCODE`; `imap ∘ map == id` inverse-permutation property;
threshold tables are strictly ascending — the §3.9.2 partial-
search band-ordering invariant; threshold lengths match
`NCODE - NCAN`; helper-vs-constant equivalence).

1 new integration test against the staged conformance corpus:
`gain_reconstruct_in_domain_on_full_corpus` walks every `.BIT`
file in `g729-core/` + `g729a/`, unpacks every active frame's
parameters, runs `reconstruct_frame_gains`, and pins that every
reconstructed pair is finite and lies in the documented
plausibility envelope (`ĝ_p ∈ [0, 2]`, `γ̂ ∈ [0, 11]`). With the
round-225 corpus walker this confirms that no transmitted (GA,
GB) index pair in the ITU conformance corpus ever drives the
reconstruction off the envelope.

### Round 239 — §3.9.1 4th-order MA gain prediction

A new `oxideav_g729::gain_predict` module ties the round-231
correction-factor output `γ̂` into the actual quantised
fixed-codebook gain `ĝ_c = γ̂ · g'_c` via the spec §3.9.1 4-tap MA
predictor:

- `GainPredictor::new()` owns the 4-slot history
  `[Û^(m-1), Û^(m-2), Û^(m-3), Û^(m-4)]` initialised per spec
  Table 9 / §4.3 to `[-14, -14, -14, -14]` dB. Slot 0 is the
  most-recent (eq (69) `b_1`-weighted) slot; slot 3 is the
  oldest (`b_4`-weighted).
- `GainPredictor::codevector_energy_db(c)` computes spec eq (66)
  `E = 10·log10((1/40) · Σ_{n=0..39} c(n)^2)` for one 40-sample
  fixed-codebook contribution. Defends against the all-zero
  corner case with a `1e-30` `log10` floor so output stays
  finite.
- `GainPredictor::predict_only(c)` runs eq (66) + eq (69) +
  eq (71) without advancing the history; returns a
  `PredictedGain { e_db, e_tilde_db, g_c_prime }` for inspection.
  `predict_only_from_energy(e_db)` is the same path with a
  pre-computed `E` so callers don't have to materialise `c(n)`.
- `GainPredictor::predict_and_update(c, γ̂)` runs the full
  per-subframe path: eq (66) energy, eq (69)
  `Ẽ^(m) = Σ_{i=1..=4} b_i · Û^(m-i)` (Q13 coefficients from
  the staged `GAIN_QUANT_MA_PREDICTOR_Q13` table converted to
  `f32` at the boundary), eq (71) `g'_c = 10^((Ẽ + Ē - E)/20)`
  with `Ē = 30 dB`, eq (65) `ĝ_c = γ̂ · g'_c`, and finally
  eq (72) decode-form history advance `Û^(m) = 20·log10(γ̂)`.
  Returns `(ĝ_c, gain_path)`.
- `GainPredictor::push_quantised_error(γ̂)` is the low-level
  eq (72) history advance, exposed so callers running custom
  loops (concealment paths) can drive the predictor explicitly.

13 new unit tests pin the algorithmic invariants:

- Table 9 init: every history slot starts at `-14 dB`;
- first-subframe `Ẽ^(0) = -14 · (b1 + b2 + b3 + b4)` matches
  the staged-Q13 dot product to 1e-4 dB;
- eq (66) matches `-10 dB` for the minimum-energy 4-pulse
  codevector (4/40 mean square) and stays finite on an
  all-zero codevector (the `1e-30` floor lands at `-300 dB`);
- eq (71) unity check: `Ẽ^(m) = E - Ē` exactly cancels the
  exponent to `g'_c = 1.0` (pins the sign convention of the
  `(Ẽ + Ē - E) / 20` exponent);
- eq (72) decode form: `push_quantised_error(γ̂ = 1.0)`
  produces `Û^(m) = 0 dB` with subsequent slots shifted one
  step deeper; `γ̂ = 10` produces `Û^(m) = 20 dB`;
- 20·log10 scaling: `γ̂ = 0.1` produces `Û^(m) = -20 dB`
  (defends against a missing factor of 20 or a sign flip);
- non-positive `γ̂` floors to a finite very-negative `Û^(m)`
  rather than NaN (defensive — well-formed §3.9.2 output is
  always positive);
- end-to-end predict-and-update consistency: the
  return-value `gain_path` matches a side-by-side `predict_only`
  call on the same predictor state, and `predict_and_update`
  with `γ̂ = 1.0` yields `ĝ_c = g'_c`;
- `predict_only` is side-effect-free: history is unchanged
  after the call (pinned for the encoder-side dry-run use case);
- steady-state convergence: after `MA_NP` pushes with constant
  `γ̂`, every history slot equals `20·log10(γ̂)` and `Ẽ^(m)`
  equals that times the sum of taps;
- history index 0 is the most-recent slot (the `b_1`-weighted
  one) — built by pushing four distinct `Û` values and verifying
  the eq (69) sum weights each slot with the staged
  `GAIN_QUANT_MA_PREDICTOR_Q13[k]` for `k ∈ 0..4`;
- eq (65) sweep: with `g'_c` forced to `1.0` via energy cancel,
  `ĝ_c` equals `γ̂` exactly across the sweep `{0.25, 0.5, 1.0,
  1.5, 2.5, 4.0}`.

1 new integration test against the staged conformance corpus:
`gain_predict_finite_on_full_corpus` walks every `.BIT` file in
`g729-core/` + `g729a/`, runs the round-231
`reconstruct_frame_gains` per active frame, then drives a fresh
`GainPredictor` through the resulting `γ̂` sequence with a
representative §3.8 minimum-energy 4-pulse codevector each
subframe. Asserts that every `(g'_c, ĝ_c)` pair stays finite,
that `ĝ_c` lies in the defensive envelope `[0, 1e6]`, and that
every history slot stays finite after the eq (72) update. The
envelope is generous because the synthetic codevector understates
real `E^(m)` (which would dampen `g'_c`); the actual §3.8 decode
path wiring is left to a follow-up round.

With round 239 the §3.9.1 / §4.1.5 gain-prediction stage chains
the round-231 `γ̂` output into the actual decoder-side fixed-
codebook gain. The full decoder gain pipeline is now wired:
`(GA, GB) → reconstruct_gains → (ĝ_p, γ̂) → gain_predict ·
predict_and_update(c, γ̂) → (ĝ_c, gain_path) → §3.10 excitation
`u(n) = ĝ_p · v(n) + ĝ_c · c(n)`. The §3.8 algebraic-codebook
decode path that produces `c(n)` is the remaining piece needed
to drive the predictor with real fixed-codebook contributions.

### Round 249 — §3.9.3 gain-quantiser codeword mapping

A new `oxideav_g729::gain_index_map` module ties the round-225
transmitted-index unpacker output to the round-231 codebook lookup
via the spec §3.9.3 robustness mapping. Per clause 3.9.3:

> "The codewords GA and GB for the gain quantizer are obtained from
> the indices corresponding to the best choice. To reduce the impact
> of single bit errors the codebook indices are mapped."

The spec prose stops there; the actual permutation is given by the
staged `map1` / `map2` (encoder-side, codebook→transmitted) and
`imap1` / `imap2` (decoder-side, transmitted→codebook) tables, whose
mutual-inverse property is pinned by the `tables_shape` integration
test.

- `demap_ga(transmitted: usize) -> Result<usize, _>` and
  `demap_gb(transmitted: usize) -> Result<usize, _>` — decoder-side
  primitives that index `GAIN_QUANT_GA_INVERSE_PERMUTATION` /
  `GAIN_QUANT_GB_INVERSE_PERMUTATION` with bounds-checked input;
  out-of-range surfaces as `GainIndexMapError::{GaOutOfRange,
  GbOutOfRange}` rather than panicking.
- `map_ga(codebook: usize) -> Result<usize, _>` and
  `map_gb(codebook: usize) -> Result<usize, _>` — the symmetric
  encoder-side forward permutations. Round-tripping
  `demap_ga(map_ga(i)?) == Ok(i)` (and the reverse composition) is
  pinned by the unit tests across the full per-stage domain.
- `demap_frame(&Parameters) -> Result<DemappedGainIndices, _>` —
  per-frame wrapper that demaps the four gain indices in one call.
  `DemappedGainIndices` is a `Copy` struct carrying
  `(ga1, gb1, ga2, gb2)` in the codebook-index domain, ready to
  feed `gain_reconstruct::reconstruct_gains`.

The existing `gain_reconstruct::reconstruct_frame_gains` is now
wired through `demap_frame` before the §3.9.2 codebook lookup, so
the on-wire `(GA, GB) → (ĝ_p, γ̂)` pipeline is spec-conformant
end-to-end. A new
`gain_reconstruct::reconstruct_gains_from_transmitted(t_ga, t_gb)`
helper exposes the same demap-then-reconstruct path as a per-pair
primitive for callers working off bare integers. The
`GainReconstructError` enum gains an `IndexMap(GainIndexMapError)`
variant + `From<GainIndexMapError>` conversion so demap failures
surface through the existing error type without bypassing it.

17 new unit tests in the `gain_index_map` module pin the algorithmic
invariants:

- forward / inverse round-trip on both stages and both compositions
  (`map ∘ demap == id` and `demap ∘ map == id` across the full
  per-stage domain);
- codebook-domain containment of every demap output;
- bijection check on each demap (cover of the codebook domain);
- out-of-range boundary on every entry point;
- non-identity assertion on both demaps (locks against a CSV
  regression that emits `[0, 1, ...]` for either `imap` table);
- per-frame wrapper threads the right transmitted indices into the
  right stages without swapping `(ga, gb)` columns or reusing
  subframe-1 values in subframe 2;
- per-stage zero-index match against the staged `imap1[0]` /
  `imap2[0]` literals (ties the function contract directly to the
  staged table rather than baking the literal in the test).

3 new unit tests in `gain_reconstruct` pin the round-249 wire-up:

- `frame_wrapper_demaps_before_codebook_lookup` finds a
  transmitted-domain pair `(t_ga, t_gb)` where the §3.9.3 inverse
  permutation maps to a different codebook-domain pair, then asserts
  the frame wrapper's output matches `reconstruct_gains(demap_ga(t_ga),
  demap_gb(t_gb))` and explicitly does NOT match the bare
  `reconstruct_gains(t_ga, t_gb)`;
- `reconstruct_from_transmitted_matches_manual_pipeline` agrees with
  the hand-composed `reconstruct_gains(demap_ga(t)?, demap_gb(t)?)`
  pipeline for every `(GA, GB)` pair in the transmitted domain;
- `reconstruct_from_transmitted_rejects_out_of_range_{ga,gb}` —
  out-of-range inputs surface through the `IndexMap` error variant.

The pre-existing `frame_wrapper_threads_per_subframe_indices` test
is updated to compose `demap_ga` / `demap_gb` in its expected-value
computation so it stays semantically correct under the new wiring.

The staged `gain_predict_finite_on_full_corpus` /
`gain_reconstruct_in_domain_on_full_corpus` integration tests
continue to pass with no envelope change: every (GA, GB) pair stays
inside the `ĝ_p ∈ [0, 2]` / `γ̂ ∈ [0, 11]` window regardless of
which permutation is applied to its row labels, so the corpus
walker was insensitive to the missing demap step. The bit-exact
match against the ITU conformance corpus's expected reconstructed-
gain trace is the cross-check that surfaces the demap correctness
itself; that integration is left for a follow-up round once the
§3.7 + §3.8 decode paths land and a full subframe-level decode
trace can be compared.

With round 249 the §4.1.5 gain decode is spec-conformant from the
on-wire bits forward: `serial::parse_frame` →
`parameters::unpack_parameters` → `(transmitted GA, GB)` →
`gain_index_map::demap_frame` → `(codebook GA, GB)` →
`gain_reconstruct::reconstruct_gains` → `(ĝ_p, γ̂)` →
`gain_predict::predict_and_update` → `(ĝ_p, ĝ_c)`.

### Round 255 — §4.1.3 pitch-delay decode

A new `oxideav_g729::pitch_decode` module ties the round-225
unpacker's transmitted `P1` / `P2` codewords into the per-subframe
fractional pitch delays `(T1, T2)` that the §3.7 eq (40) `b_30`
adaptive-codebook interpolator consumes:

- `decode_t1_from_p1(p1: u8) -> PitchDelay` evaluates spec image
  `f0027-01.jpg` (clause 4.1.3): the union of a fractional branch
  (`P1 < 197`, 1/3-resolution `T1 ∈ [19⅓, 85]`) and an
  integer-only branch (`P1 ≥ 197`, `T1 ∈ [86, 143]`). Every
  `P1 ∈ 0..256` codeword maps to a unique `(int_t, frac)` pair.
- `derive_t_min(int_t1: i32) -> i32` evaluates spec image
  `f0027-02.jpg`: clamps the subframe-2 search window's lower
  bound into `[T_MIN_FLOOR = 20, T_MAX_CEIL − T_WINDOW = 134]`
  per the spec recipe (`t_min = int(T1) − 5`, floor at 20,
  pulled down by `t_max = 143` − 9 if `int(T1)` is near the
  ceiling).
- `decode_t2_from_p2(p2: u8, t_min: i32) -> PitchDelay` evaluates
  spec image `f0027-03.jpg`: every `P2 ∈ 0..32` codeword maps to
  a unique `(int_t, frac)` pair inside `int(T2) ∈ [t_min − 1,
  t_min + 10]`, `frac ∈ {-1, 0, 1}`.
- `decode_frame(&Parameters) -> FramePitchDelays` — per-frame
  wrapper that chains the three steps in spec §4.1.3 order.
  Returns the per-subframe pitch delays AND the spec-derived
  `t_min` (preserved for callers driving the §4.1.2
  parity-concealment path).
- `encode_p1(delay) -> Option<u8>` and `encode_p2(delay, t_min)
  -> Option<u8>` — the spec §3.7 encode-side forward mappings
  (eqs (41) / (42) per spec images `eq41.jpg` / `eq42.jpg`)
  exposed publicly. The unit-test suite uses them for full-domain
  round-trip; encoder wire-up will use them in a future round.

13 unit tests pin the algorithmic invariants:

- Spec-image worked examples on `decode_t1_from_p1` (every
  branch boundary `P1 ∈ {0, 1, 196, 197, 255}` pinned literally);
- full-domain envelope on `decode_t1_from_p1` (every `P1 ∈
  0..256` lands in `int(T1) ∈ [19, 143]`, `frac ∈ {-1, 0, 1}`);
- spec-image worked examples on `derive_t_min` (mid-range, floor
  edge, ceiling edge);
- full-domain envelope on `derive_t_min` (the entire 9-step
  subframe-2 search window fits inside `[20, 143]` for every
  `int(T1) ∈ [19, 143]`);
- spec-image worked examples on `decode_t2_from_p2` (`P2 ∈ {0,
  2, 31}` at `t_min = 50`);
- full-domain envelope on `decode_t2_from_p2` (`P2 ∈ 0..32` ×
  `t_min ∈ [20, 134]` always lands inside the spec-stated
  `(int_t, frac)` window);
- **encode↔decode round-trip across the full `P1 ∈ 0..256`
  domain** — pins both the eq (78a) decode and the eq (41)
  encode recipes simultaneously; any sign / off-by-one drift
  in either direction trips the round-trip;
- **encode↔decode round-trip across the full `P2 ∈ 0..32` ×
  `t_min ∈ [20, 134]` domain** — same property for the
  subframe-2 differential mapping;
- encode-side out-of-domain rejection on both `encode_p1` and
  `encode_p2`;
- `decode_frame` threads the right field into the right
  subframe (P1↔P2 swap detection);
- constants surface (`T_MIN_FLOOR = 20`, `T_MAX_CEIL = 143`,
  `T_WINDOW = 9`, `P1_DOMAIN = 256`, `P2_DOMAIN = 32`,
  `P1_FRACTIONAL_LIMIT = 197`) matches the documented spec
  values.

2 integration tests against the staged conformance corpus:

- `pitch_decode_in_domain_on_full_corpus` walks every `.BIT`
  file in `g729-core/` + `g729a/`, unpacks each active frame's
  `Parameters`, runs `decode_frame`, and pins that every decoded
  `(T1, T2, t_min)` lies in the spec-stated domain.
- `pitch_decode_round_trips_pitch_corpus` walks the staged
  `PITCH.BIT` sequence (the ITU `READMETV.txt` self-documents
  this as the pitch-delay exerciser) for both `g729-core/` and
  `g729a/`, runs `decode_frame` on every active frame, and pins
  that `encode_p1(t1) == params.p1` and `encode_p2(t2, t_min) ==
  params.p2` exactly. This is the strongest in-corpus guarantee
  available without a known-good reference output: a single
  frame failing means the decode and encode recipes disagree on
  an ITU-encoded frame.

With round 255 the §4.1 decode chain extends another link from the
on-wire bits: `serial::parse_frame` →
`parameters::unpack_parameters` → `(P1, P2)` →
`pitch_decode::decode_frame` → `(T1, T2, t_min)`. The remaining
§4.1.4 path (`C` / `S` → pulse positions + signs → `c(n)`) and the
§3.7 eq (40) `b_30` interpolator (`(int(T), frac, past
excitation) → v(n)`) are the next links the decoder chain needs
before the per-subframe excitation `u(n) = ĝ_p · v(n) + ĝ_c ·
c(n)` (spec eq (75)) can be evaluated.

### Round 266 — §3.8 / §4.1.4 fixed (algebraic) codebook decode

A new `oxideav_g729::fixed_codebook` module ties the round-225
parameter unpacker's `(C1, S1, C2, S2)` codewords into the
per-subframe pulse layout that spec eq (45) builds into the
40-sample codevector `c(n)`:

- `decode_positions(c: u16) -> Result<([u8; 4], u8), _>` evaluates
  the inverse of spec eq (62): `C = (m_0/5) + 8·(m_1/5) +
  64·(m_2/5) + 512·(2·(m_3/5) + jx)`. The 13-bit field splits as
  3 + 3 + 3 + 4 bits across the four Table-7 tracks (LSB-first
  per track); the 4-bit track-3 field carries the 3-bit `m_3/5`
  position index in its upper 3 bits and the `jx` track-selector
  bit in its LSB. Returns the `[m_0, m_1, m_2, m_3]` positions
  in spec `(i_0, i_1, i_2, i_3)` order together with the
  decoded `jx ∈ {0, 1}`.
- `decode_signs(s: u8) -> Result<[i8; 4], _>` evaluates the
  inverse of spec eq (61): `S = s_0 + 2·s_1 + 4·s_2 + 8·s_3`.
  Bit `k` of `S` carries `s_k = 1` ⇔ positive sign per the
  clause-3.8.2 prose; the returned `i8` is `+1` for `s_k = 1`
  and `-1` for `s_k = 0`, matching spec eq (45)'s `s_k ∈ {-1, +1}`
  amplitude convention.
- `decode_pulses(c, s) -> Result<FixedCodebookPulses, _>` ties
  the position and sign primitives together — the per-subframe
  §4.1.4 entry point.
- `build_codevector(&FixedCodebookPulses) -> [i8; 40]` constructs
  spec eq (45) — four signed unit impulses placed on the 40-sample
  codevector, zeros elsewhere. The energy of the result is exactly
  `4` per the unit invariant.
- `decode_frame(&Parameters) -> Result<FrameFixedCodebook, _>` —
  per-frame wrapper that threads `(C1, S1)` into subframe 1 and
  `(C2, S2)` into subframe 2 per spec §4.1.5.
- `encode_positions(&[u8; 4]) -> Option<u16>` and
  `encode_signs(&[i8; 4]) -> Option<u8>` — the symmetric encode-
  side forward mappings (spec eqs (62) / (61)) exposed publicly
  for round-trip property tests and encoder wire-up.

12 new unit tests pin the algorithmic invariants:

- spec eq (61) worked examples on `decode_signs` (`S = 0b0000` /
  `0b1111` / `0b0001` / `0b1000` / `0b0101`);
- spec eq (62) worked examples on `decode_positions` (`C = 0` /
  `1` / `512` (`jx = 1, m_3 = 4`) / `0x1FFF` (full saturation,
  `m = [35, 36, 37, 39]`, `jx = 1`));
- out-of-domain rejection on both decode primitives (`S` with
  bit 4 set, `C` with bit 13 set);
- full-domain envelope on `decode_positions` (every `C ∈
  0..8192` produces 4 positions on the spec Table-7 tracks
  with the correct per-track residue);
- full-domain encode↔decode round-trip on `S` (every `S ∈
  0..16`) AND on `C` (every `C ∈ 0..8192`) — pins both the
  decode and encode recipes simultaneously;
- `encode_positions` rejects off-track inputs (track-0 position
  at track-1 residue, track-3 at residue {0, 1, 2}, position ≥ 40);
- `encode_signs` rejects non-`±1` inputs (`0`, `2`);
- `decode_pulses` threads positions + signs in spec order;
- `build_codevector` constructs eq (45) — exactly 4 non-zero
  ±1 samples;
- codevector energy is exactly `NUM_PULSES = 4` across a sweep
  of representative `(C, S)` pairs;
- the four pulses occupy distinct positions across the full `C`
  domain (the Table-7 tracks are disjoint modulo 5);
- `decode_frame` threads `(C1, S1)` and `(C2, S2)` into the
  right subframes (no swap detection).

2 new integration tests against the staged conformance corpus
(in `tests/serial_conformance.rs`):

- `fixed_codebook_in_domain_on_full_corpus` walks every `.BIT`
  file in `g729-core/` + `g729a/`, runs `decode_frame` on every
  active frame, and pins per-subframe invariants:
  every pulse position is in `[0, 40)` and on its spec Table-7
  track residue; every pulse sign is `±1`; the four pulses sit
  on distinct positions; the codevector energy `Σ_n c(n)²`
  equals exactly `4`; `jx ∈ {0, 1}`;
- `fixed_codebook_round_trips_fixed_corpus` walks the staged
  `FIXED.BIT` sequence (the ITU `READMETV.txt` self-documents
  this as the fixed-codebook search exerciser) for both
  `g729-core/` and `g729a/`, runs `decode_frame` on every active
  frame, and pins that `encode_positions(positions) == params.c1`
  / `c2` and `encode_signs(signs) == params.s1` / `s2` exactly.
  This is the strongest in-corpus guarantee available without a
  known-good reference output: a single frame failing means the
  eq (61) / (62) decode and encode recipes disagree on an
  ITU-encoded frame.

With round 266 the §4.1.4 decode-side fixed-codebook stage
chains the round-225 transmitted `(C, S)` codewords into the
spec eq (45) codevector `c(n)` that the round-239 gain
predictor's eq (66) energy step ultimately consumes. The
spec §3.8 eq (48) / (48a) pitch sharpening (applied when the
integer pitch delay `int(T) < 40`) needs the round-255
pitch delay AND the previous subframe's quantised
adaptive-codebook gain `β`, and is deferred to a follow-up
round; it modifies the codevector in place after `build_codevector`
but does not change the per-pulse layout this round wires.

### Round 282 — §4.1 per-frame decode parameter chain

Round 274 had landed the deferred §3.8 / §4.1.4 pitch sharpening
(`pitch_sharpen::clamp_beta` for eq (47), `pitch_sharpen::sharpen`
for the eq (48) recursive in-place modification when
`int(T) < 40`). Round 282 glues **everything above into one
stateful per-frame call** — a new `decode_chain` module whose
`FrameDecoder` owns all four cross-frame state pieces with their
clause 4.3 / Table 9 start-up values:

| state | role | init (Table 9) |
|---|---|---|
| `LspReconstructor` | §3.2.4 4-frame MA residual history | `î_i = iπ/11` |
| `LspInterpolator` | §3.2.5 previous-frame LSPs | `q_i = cos(iπ/11)` |
| `GainPredictor` | §3.9.1 4-tap `Û` history | `−14 dB` |
| `g_p_prev` | eq (47) `β = ĝ_p^(m−1)` source | `0.8` |

plus the previous frame's `int(T2)` (clause-4.3 zero default) for
the §4.1.2 parity concealment ("if a parity error occurs on P1,
the delay value T1 is set to the integer part of the delay value
T2 of the previous frame; the value T2 is derived with the
procedure outlined in clause 4.1.3, using this new value of T1" —
now wired). Three entry points: `decode_serial_frame` (164-byte
ITU serial frame), `decode_frame_kind` (parsed `FrameKind`), and
`decode_parameters` (unpacked Table-8 codewords). Each call runs
the clause-4.1 order — §4.1.1 LSP→LP per subframe, §4.1.2 parity
check + concealment substitution, §4.1.3 pitch delays, §4.1.4
fixed codebook + eq (48) sharpening (β threaded from the previous
subframe's `ĝ_p`), §4.1.5 gains (`ĝ_p`, `γ̂`, `ĝ_c = γ̂·g′_c`
with the eq (66) energy taken over the harmonic-enhanced `c(n)`
per §3.10) — and returns fully-typed `DecodedFrame` /
`SubframeDecode` structs. Erasure sentinels return
`FrameDecodeError::Erased` (§4.4 concealment is a future round).

8 new unit tests (Table-9 start-up state, hand-sequenced
piece-by-piece equivalence, §4.1.2 substitution incl. the
first-frame zero-default, erasure rejection without state
advance, typed out-of-domain errors, serial↔parameter entry-point
equality, eq (47) β clamp under gain sweeps) and a new
`tests/decode_chain_conformance.rs` integration harness that runs
the chain over **every** active frame of the staged base + Annex-A
corpus — 18 222 frames across 19 `.BIT` vectors — pinning the
§3.2.4 stability-clamp envelope on `ω̂`, the §4.1.3 delay windows,
the Table-7 track residues, finiteness of every gain/LP output,
exactly 60/300 §4.1.2 concealment activations on `PARITY.BIT`
(each substituting exactly the previous frame's `int(T2)`), and
frame-by-frame determinism of two independent chains.

Round 290 lands the **first decoder stage to emit reconstructed
speech** — the §4.1.6 LP synthesis path, in a new `lp_synthesis`
module. A stateful `Synthesizer` owns the two cross-subframe state
pieces clause 4.3 lists as zero-initialised: the eq (40)
past-excitation buffer (`EXC_HISTORY` = 153 samples — the deepest
eq (40) tap `u(n − k − i)` with `n = 0`, the §3.7 maximum integer
delay `k = 143`, `i = 9` lands at index `−152`) and the eq (77)
10th-order `1/Â(z)` synthesis-filter memory.
`synthesize_frame(&DecodedFrame)` consumes one round-282 decoded
frame and runs the spec §4.1.3 → §3.10 → §4.1.6 order per subframe:

- **eq (40)** interpolates the past excitation through the 31-tap
  `b_30` (`PITCH_INTERP_FILTER_SYNTHESIS_Q15`) at the decoded
  fractional pitch delay, mapping the decoder's `(int(T), frac)`
  (`frac ∈ {−1,0,1}`) onto the eq (39)/(40) `(k, t)` form
  (`t ∈ {0,1,2}` for fractions 0, 1/3, 2/3) to build `v(n)`;
- **eq (75)** forms the excitation `u(n) = ĝ_p·v(n) + ĝ_c·c(n)`
  (with `c(n)` the post-eq (48) harmonic-enhanced codevector);
- **eq (77)** filters it through the synthesis filter
  `ŝ(n) = u(n) − Σ_{i=1}^{10} â_i·ŝ(n − i)`.

The `v(n)`/`u(n)` build interleaves per sample so short delays
(`T < 40`) fold correctly onto the already-built current
excitation; both state buffers advance after every subframe.
Typed outputs `SynthesizedFrame` (two subframes + `speech()` →
the 80 time-ordered samples) and `SynthesizedSubframe` (the
`adaptive` `v(n)`, `excitation` `u(n)`, and `speech` `ŝ(n)`
40-sample vectors). 10 new unit tests pin the eq (40)/(75)/(77)
algebra, the fraction-convention mapping, the state advance, the
short-delay fold, and determinism; a new
`tests/lp_synthesis_conformance.rs` harness asserts every
`v(n)`/`u(n)`/`ŝ(n)` stays finite across all active frames of the
base + Annex-A corpus.

## What is NOT wired up

Every decode/encode entry point still returns `Error::NotImplemented`.
The remaining numeric tables (gain-quantizer coefficient matrix
`coef` / `L_coef`, postfilter interpolation `tab_hup_*`, taming
`tab_zone`, Annex B DTX/CNG, LSF↔LSP cos/slope tables) are staged
under `docs/audio/g729/tables/` but are not yet compiled in; the
Implementer leaves them out until the docs collaborator's
specifier pass clarifies the per-clause wire-up direction.

With round 266 the §4.1 / Table-8 parameter unpacker chains the
round-191 framing layer to the §3.2.4 / §3.2.5 / §3.2.6 LSP
decode chain AND the §3.9.2 / §4.1.5 gain VQ decode chain AND
the §3.9.1 / §4.1.5 gain-prediction stage AND the §3.9.3
codeword-mapping robustness layer AND the §4.1.3 pitch-delay
decode AND the §3.8 / §4.1.4 fixed-codebook decode end-to-end:
`serial::parse_frame` → `parameters::unpack_parameters` →
`(L0, L1, L2, L3)` → `LspReconstructor::reconstruct_frame` →
`LspInterpolator::interpolate` → `lsp_to_lp` per subframe,
with `(transmitted GA1, GB1, GA2, GB2)` →
`gain_index_map::demap_frame` → `(codebook GA, GB)` →
`gain_reconstruct::reconstruct_frame_gains` →
`gain_predict::GainPredictor::predict_and_update` yielding the
per-subframe `(ĝ_p, ĝ_c)` pairs, `(P1, P2)` →
`pitch_decode::decode_frame` → per-subframe `(T1, T2, t_min)`
fractional pitch delays, AND `(C1, S1, C2, S2)` →
`fixed_codebook::decode_frame` → per-subframe pulse layout +
`build_codevector` → `c(n)` codevector → `pitch_sharpen::sharpen`
(round 274) — and round 282's `decode_chain::FrameDecoder`
sequences all of it as one stateful per-frame call, whose output
round 290's `lp_synthesis::Synthesizer` turns into reconstructed
speech `ŝ(n)` via the §3.7 eq (40) `b_30` past-excitation
interpolator (`v(n)`), the §3.10 / §4.1.6 per-subframe excitation
`u(n) = ĝ_p · v(n) + ĝ_c · c(n)` build, and the eq (77) `1/Â(z)`
synthesis. The remaining decode-side work is the §4.2
post-processing cascade (long-/short-term postfilter, tilt
compensation, adaptive gain control, output high-pass + ×2
upscaling) and the §4.4 frame-erasure concealment. The
remaining numeric tables (gain-quantizer coefficient matrix
`coef` / `L_coef`, postfilter interpolation `tab_hup_*`, taming
`tab_zone`, Annex B DTX/CNG, LSF↔LSP cos/slope tables) are
staged under `docs/audio/g729/tables/` but not yet compiled in;
the Implementer leaves them out until the docs collaborator's
specifier pass clarifies the per-clause wire-up direction.

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
