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
per eq (72) decode form `Û^(m) = 20·log10(γ̂)`.

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

## What is NOT wired up

Every decode/encode entry point still returns `Error::NotImplemented`.
The remaining numeric tables (gain-quantizer coefficient matrix
`coef` / `L_coef`, postfilter interpolation `tab_hup_*`, taming
`tab_zone`, Annex B DTX/CNG, LSF↔LSP cos/slope tables) are staged
under `docs/audio/g729/tables/` but are not yet compiled in; the
Implementer leaves them out until the docs collaborator's
specifier pass clarifies the per-clause wire-up direction.

With round 239 the §4.1 / Table-8 parameter unpacker chains the
round-191 framing layer to the §3.2.4 / §3.2.5 / §3.2.6 LSP
decode chain AND the §3.9.2 / §4.1.5 gain VQ decode chain AND
the §3.9.1 / §4.1.5 gain-prediction stage end-to-end:
`serial::parse_frame` → `parameters::unpack_parameters` →
`(L0, L1, L2, L3)` → `LspReconstructor::reconstruct_frame` →
`LspInterpolator::interpolate` → `lsp_to_lp` per subframe,
with `(GA1, GB1, GA2, GB2)` →
`gain_reconstruct::reconstruct_frame_gains` →
`gain_predict::GainPredictor::predict_and_update` yielding the
per-subframe `(ĝ_p, ĝ_c)` pairs. The remaining transmitted
indices `(P1, P0, P2, C1, S1, C2, S2)` are *available* via the
same `Parameters` struct but the §3.7 / §3.8 decode-side
algorithms that consume them are not yet wired (§3.7 maps `P1`
→ fractional pitch delay; §3.8 maps `C1` → pulse positions —
producing the `c(n)` codevector that the round-239 gain
predictor's energy step consumes). The remaining numeric
tables (gain-quantizer coefficient matrix `coef` / `L_coef`,
postfilter interpolation `tab_hup_*`, taming `tab_zone`,
Annex B DTX/CNG, LSF↔LSP cos/slope tables) are staged under
`docs/audio/g729/tables/` but not yet compiled in; the
Implementer leaves them out until the docs collaborator's
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
