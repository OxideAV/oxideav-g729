# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the crate adheres
to [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Round 218 wires the §3.2.6 LSP-to-LP conversion on top of the
  round-213 per-subframe cosine-domain LSP output, in a new
  `oxideav_g729::lsp_to_lp` module:
  - `lsp_to_lp(q_in: &[f32; 10]) -> [f32; 10]` runs the spec
    §3.2.6 three-step recipe: (1) build the F1/F2 sum/difference
    polynomial coefficients `f_1(i)`, `f_2(i)` for `i ∈ 0..=5`
    via the spec recursion derived from polynomial multiplication
    by `(1 − 2·q·z^-1 + z^-2)` (spec convention `f(-1) = 0`);
    (2) restore the `(1 ± z^-1)` factors per eq (25)
    (`f'_1(i) = f_1(i) + f_1(i-1)`,
    `f'_2(i) = f_2(i) − f_2(i-1)`); (3) recombine via eq (26)
    (`a_i = ½·f'_1(i) + ½·f'_2(i)` for `i ∈ 1..=5`,
    `a_i = ½·f'_1(11-i) − ½·f'_2(11-i)` for `i ∈ 6..=10`, using
    `F'_1` symmetric / `F'_2` antisymmetric of length 6).
  - `lsf_to_lp(omega: &[f32; 10]) -> [f32; 10]` — boundary
    wrapper from the LSF domain `ω̂` via `q_i = cos(ω̂_i)`.
  - `pub type LpCoefficients = [f32; 10]` — output type alias;
    slot `i - 1` holds `a_i` for `i ∈ 1..=10`, with `a_0 = 1.0`
    implicit (not stored).
  - 8 new unit tests pin the algorithmic invariants: start-up
    state produces finite coefficients; recursion matches a
    brute-force polynomial-multiplication oracle on two LSP
    patterns to ≤ 1e-4 drift (independent re-derivation of
    eqs (13) / (14) by literal accumulation, then eqs (25) /
    (26) applied in the test, locks the in-place inner-loop
    ordering against read-after-write bugs); closed-form spot
    checks at `z = 1` (`A(1) = Π_{odd}(2 − 2·q_i)`, since
    `F'_2(1) = 0`) and at `z = -1`
    (`A(-1) = Π_{even}(2 + 2·q_i)`, since `F'_1(-1) = 0`) to
    ≤ 1e-3 drift; coefficient range stays in `±32` (defence
    against a missing ½ factor); `lsf_to_lp` wrapper matches
    the explicit `lsp_to_lp(&omega_to_q(&omega))` to ≤ 1e-7;
    full §3.2.4 → §3.2.5 → §3.2.6 chain on a 3-frame
    non-steady-state `(L0, L1, L2, L3)` input produces finite
    per-subframe `a_i` AND distinct subframe-1 / subframe-2
    filters; all-zero `q` corner case (`ω̂_i = π/2` for all i)
    produces finite coefficients (drift check on the spec
    `f(-1) = 0` boundary).

- Round 213 wires the §3.2.5 per-subframe LSP interpolation
  (spec eq (24)) on top of the round-207 §3.2.4 reconstructor, in
  a new `oxideav_g729::lsp_interpolate` module:
  - `omega_to_q(omega: &[f32; 10]) -> [f32; 10]` and
    `q_to_omega(q: &[f32; 10]) -> [f32; 10]` — the boundary
    helpers between the LSF domain `ω̂ ∈ [0, π]` and the cosine
    domain `q_i = cos(ω̂_i)`. `q_to_omega` clamps inputs to
    `[-1, 1]` before `acos` so float-rounding past the boundary
    by ~1e-7 cannot produce NaN.
  - `pub const SUBFRAMES_PER_FRAME: usize = 2` — spec §2.1
    sub-frame count exposed at the module surface.
  - `LspInterpolator` — stateful interpolator that carries the
    previous frame's cosine-domain LSPs `q_i^(previous)`.
    `new()` initialises `previous_q` to the cosine of the spec
    §3.2.4 start-up LSFs `ω̂_i = i · π / 11`, matching the
    round-207 reconstructor's start-up state.
    `interpolate(&[f32; 10]) -> [[f32; 10]; 2]` applies eq (24):
    `q^(1)[i] = 0.5 · q^(previous)[i] + 0.5 · q^(current)[i]` and
    `q^(2)[i] = q^(current)[i]`, then advances
    `previous_q := current_q`. `interpolate_from_omega(&[f32;
    10])` is a convenience entry point that wraps the boundary
    conversions for callers staying in the LSF domain (the actual
    interpolation is still done in the cosine domain per spec
    §3.2.5).
  - 8 new unit tests pin the algorithmic invariants: `omega↔q`
    round-trip identity on real LSFs (1e-5 tolerance); `acos`
    boundary clamp on `q = ±(1 + ε)` (`acos(1+ε) ≈ 0`,
    `acos(-1-ε) ≈ π`); start-up `previous_q` matches
    `cos((i+1)π/11)` to 1e-6; subframe 2 == current frame's `q`
    exactly; subframe 1 == per-coordinate midpoint of
    `(previous, current)` exactly; `previous_q` advances to
    `current` after each call (steady-state two-frame pass
    produces sub1 == sub2); `interpolate_from_omega` matches the
    explicit cosine-domain pipeline; and an end-to-end test that
    drives the §3.2.5 interpolator from three successive
    round-207 reconstructed frames `(L0, L1, L2, L3)` ∈ {(0, 0,
    0, 0), (1, 5, 7, 11), (0, 12, 3, 17)} and verifies subframe-2
    equals the frame's reconstructed LSF and that subframe-1
    (after `q_to_omega` inversion) lies between the previous and
    current frame's LSFs at every coordinate (cos is monotone on
    `[0, π]`, so the cosine-domain midpoint translates to a
    *between*-in-omega result).

- Round 207 wires the §3.2.4 LSP-frame reconstruction algorithm
  around the round-195 / round-201 tables, in a new
  `oxideav_g729::lsp_reconstruct` module:
  - `codebook_sum(l1, l2, l3) -> [f32; 10]` evaluates spec eq (19)
    (`l̂_i = L1_i(L1) + L2_i(L2)` for `i ∈ 1..=5`, `L1_i(L1) +
    L3_{i-5}(L3)` for `i ∈ 6..=10`). The Q13 codebook literals are
    converted to `f32` at the boundary (`v / 8192.0`); out-of-range
    indices surface as the typed `L1OutOfRange` / `L2OutOfRange` /
    `L3OutOfRange` variants of `LspReconstructError` rather than
    panicking.
  - `rearrange_pass(coefs, j)` performs the spec §3.2.4 figure
    `F0013-01` adjacent-pair fix-up — for `i = 1..10`, if
    `l̂_{i-1} > l̂_i − J` the pair is replaced with
    `((l̂_i + l̂_{i-1}) − J)/2` and `((l̂_i + l̂_{i-1}) + J)/2`.
    `rearrange_twice(coefs)` runs the two passes the spec calls
    for (`J = REARRANGE_J1 = 0.0012`, then
    `J = REARRANGE_J2 = 0.0006`).
  - `stability_clamp(coefs)` applies the spec §3.2.4 4-step
    stability check: (1) sort ascending; (2) floor `ω̂_1` at
    `CLAMP_FLOOR = 0.005`; (3) enforce a minimum adjacent gap
    `CLAMP_MIN_GAP = 0.0391` across `i = 2..=10`; (4) ceil
    `ω̂_10` at `CLAMP_CEIL = 3.135`. All four constants are the
    `pub const` spec literals.
  - `LspReconstructor` carries the 4-frame MA history (`[[f32; 10];
    MA_NP]`). `new()` initialises every history slot to the spec
    start-up vector `l̂_i = i · π / 11` for `i ∈ 1..=10` per spec
    §3.2.4 ("at start up the initial values of `l̂_i^(m-k)` are
    given by `l̂_i = iπ/11` for all `k < 0`").
    `reconstruct_frame(l0, l1, l2, l3)` runs codebook-sum →
    `rearrange_twice` → eq (20) MA prediction (using
    `LSP_MA_PREDICTOR_FG_Q15` and the round-201 `fg_sum` factor) →
    `stability_clamp`, advances the MA history (pushing the
    post-rearrange residual into slot 0), and returns the
    reconstructed `ω̂^(m)` LSF vector.
  - 12 new unit tests pin the algorithmic invariants: spec start-up
    vector, codebook-sum boundary conversion at `(0, 0, 0)`, every
    error variant on out-of-range indices, rearrange-pass minimum
    distance and no-op stability, `rearrange_twice` finishing at
    the `J2 = 0.0006` margin, stability clamp's floor + gap + ceil
    on a deliberately violating input and no-op on a clean LSF
    vector, end-to-end clamp compliance, MA-history shift across
    one frame, and clamp compliance under both `L0` predictor
    modes.
- Round 201 completes the §3.2.4 LSP-reconstruction inputs by
  wiring the MA-predictor `fg` family (spec eq (20) / (20a)):
  - `tables::LSP_MA_PREDICTOR_FG_Q15` — 3-D Q15 coefficient cube,
    shape `[[[i16; 10]; 4]; 2]`. Outer dim is the `L0` predictor
    mode; middle dim is the MA history depth (`MA_NP = 4`); inner
    dim is the LP order `M = 10`.
  - `tables::LSP_MA_PREDICTOR_FG_SUM_Q15` — per-mode `(Q15_ONE −
    Σ_k fg[mode][k][i])` Q15 factor, shape `[[i16; 10]; 2]`.
  - `tables::LSP_MA_PREDICTOR_FG_SUM_INV_Q12` — per-mode Q12
    reciprocal of `fg_sum`, shape `[[i16; 10]; 2]`. Pre-tabulated
    to avoid a per-sample division during reconstruction.
  - New spec-dimension constant `MA_NP = 4` (LSP MA prediction
    order); doc on `L0_BITS` updated to point at the new tables.
  - Bounds-checked lookup helpers `lsp_fg_plane(mode)`,
    `lsp_fg_sum(mode)`, `lsp_fg_sum_inv(mode)` returning borrowed
    per-mode plane / row references.
- `build.rs` gains a `Shape::Cube { planes, rows, cols }` table
  type and a `parse_cube_csv` helper, so 3-D coefficient slabs
  emit as `[[[i16; cols]; rows]; planes]` arrays. The CSV layout
  is `planes` lines, each carrying `rows * cols` comma-separated
  literals in row-major order within each plane; both the line
  count and the per-line literal count are asserted against the
  declared shape — any CSV drift trips the build with the
  offending stem in the error.
- `tests/tables_shape.rs` grows with 8 new round-201 tests:
  shape (`2 × MA_NP × M` for `fg`, `2 × M` for the two sum rows)
  cross-checked against `1 << L0_BITS == 2`; first row of each
  `fg` plane pinned to CSV literals (drift check on the 3-D
  cube reader's row-major flattening); history-depth peak decay
  (`fg[mode][MA_NP - 1]` peak magnitude is strictly less than the
  `fg[mode][0]` peak in both modes); strict positivity across all
  80 `fg` entries (sign-flip drift check); `fg_sum` matches the
  spec-stated `(Q15_ONE − Σ_k fg)` factor within 4 Q15 ulps;
  `fg_sum_inv` is the Q12 reciprocal of `fg_sum` within 3 Q12
  ulps; and the three new helpers each return slices equal to
  the underlying constants. The `all_tables_are_non_empty`
  smoke also gains the three new round-201 constants and the
  round-195 LSP codebook constants that it previously missed.

### Changed

- `build.rs` no longer emits per-`pub const` rustdoc lines naming
  the staged electronic attachment's source filename, original C
  identifier, or per-file SHA-256. The provenance chain itself
  is unchanged (it still lives in the `.meta` sidecars under
  `docs/audio/g729/tables/` and in the CHANGELOG's round notes);
  only the in-`src/` rustdoc emission is scrubbed to keep
  algorithmic-source filenames out of the generated documentation
  surface. Existing round-189 / round-195 doc comments lose the
  `Source file inside …` and `Original C identifier` lines as a
  side-effect of the same rebuild.

- Round 195 wires up the §3.2.4 LSP-quantiser two-stage VQ
  codebooks (the `lspcb1` / `lspcb2` tables of the staged trace
  doc) and lockable lookup helpers:
  - `tables::LSP_QUANT_CODEBOOK_L1_Q13` — first-stage codebook,
    shape `[[i16; 10]; 128]` Q13 (the 7-bit `L1` index).
  - `tables::LSP_QUANT_CODEBOOK_L2_Q13` — packed second-stage
    codebook, shape `[[i16; 10]; 32]` Q13 (5-bit `L2` lower-5
    coefficients / `L3` upper-5).
  - New spec-dimension constants: `NC0 = 128`, `NC1 = 32`,
    `L0_BITS = 1`, `L1_BITS = 7`, `L2_BITS = L3_BITS = 5`,
    `LSP_TOTAL_BITS = 18`.
  - Bounds-checked lookup helpers `lsp_l1_entry(l1)`,
    `lsp_l2_entry(l2)`, `lsp_l3_entry(l3)` returning borrowed
    slices into the compiled codebooks.
- `build.rs` gains a `Shape::Matrix { rows, cols }` table type and
  a comma-separated row parser, so 2-D codebooks emit as
  `[[i16; cols]; rows]` arrays directly. Row count and per-row
  column count are both asserted against the declared shape — any
  CSV drift trips the build with the offending stem in the error.
- `tests/tables_shape.rs` grows with 8 new round-195 tests:
  shape (NC0 × M, NC1 × M), bit-width derivations
  (`1 << L1_BITS == NC0`, same for L2 / L3), pinned literals for
  the first three L1 rows + first L2 row (matrix-reader drift
  check), Q13 LSF-domain range check across all 1280 L1 entries,
  helper-vs-constant equivalence, and L2 + L3 helper concatenation
  recovering the packed row.
- `tests/serial_conformance.rs` grows with 2 new round-195 tests:
  `lsp_conformance_indices_are_in_codebook_range` walks the
  staged `LSP.BIT` vector for both `g729-core/` and `g729a/`,
  extracts (L0, L1, L2, L3) per spec Table 8 NOTE (MSB-first per
  parameter), and asserts each frame's L1 / L2 / L3 lies in
  codebook range — also smoke-testing the bounds-checked lookup
  helpers across every transmitted index in the ITU corpus.
  `lsp_indices_helper_round_trips_first_active_frame_bits`
  synthesises a frame with known (L0, L1, L2, L3) and locks the
  MSB-first packing convention independently of the corpus, so
  the bit ordering is checked in published-crate mode too.

## [0.0.6](https://github.com/OxideAV/oxideav-g729/releases/tag/v0.0.6) - 2026-05-30

### Other

- r191 — ITU serial bitstream parser + conformance-corpus harness
- r189 — wire up LP-analysis / LSF grid / pitch interp / MA gain pred tables
- bit-exact tables foundation — HPF coefs, Table 8, basic_op LUTs
- orphan rebuild — clean-room reset to register-only scaffold

### Added

- Round 191 wires up the ITU serial bitstream format used by the
  staged conformance corpus at `docs/audio/g729/conformance/`:
  - New `serial` module with public `SYNC_WORD = 0x6B21`,
    `BITS_HEADER = 80`, `BIT_ZERO = 0x007F`, `BIT_ONE = 0x0081`,
    `BIT_ERASED = 0x0000`, `FRAME_WORDS = 82`, `FRAME_BYTES = 164`
    constants. The framing-literal values are empirically observed
    in the staged `.bit` files; the 164-byte cadence is documented
    in `docs/audio/g729/conformance/README.md`.
  - `serial::parse_frame(&[u8]) -> Result<FrameKind, SerialError>`
    distinguishes normal frames (`FrameKind::Active([bool; 80])`)
    from frame-erasure sentinels (`FrameKind::Erased`), and rejects
    wrong-length, wrong-sync, wrong-header, invalid-bit-word, and
    mid-frame mixed-erasure inputs.
  - `serial::frame_count(&[u8])` byte-length cross-check.
  - 13 new unit tests in `src/serial.rs` cover the happy paths and
    all five error variants.
  - 6 new integration tests in `tests/serial_conformance.rs` walk
    the staged `docs/audio/g729/conformance/{g729-core,g729a}/`
    corpus when present, validating sync + header per frame,
    matching `.bit` frame count against `.pst` PCM frame count, and
    pinning the erasure-sentinel frame count exactly for each
    decoder-only sequence (`ERASURE` 60/300, `OVERFLOW` 1/384,
    `PARITY` 0/300 — same counts on the Annex-A set). The harness
    cleanly skips when the corpus path is absent (published-crate
    build mode), so `cargo test` stays green either way.

- Round 189 extends the bit-exact numeric-tables wire-up with the
  §3.2.1 LP-analysis windowing tables, the §3.2.5 LSF cosine grid,
  the §3.7 pitch interpolation filters, and the §3.9 MA gain-
  prediction coefficient vector. Same CSV/meta provenance chain as
  the round-173 entries (`docs/audio/g729/tables/` → crate's
  `tables/` → `build.rs` → `OUT_DIR/<stem>.rs` →
  `src/tables.rs::include!`):
  - `tables::LPC_HAMMING_WINDOW_Q15` — §3.2.1 LP-analysis Hamming
    window (`hamwindow`, `[i16; 240]`).
  - `tables::LPC_LAG_WINDOW_HIGH_Q15` — §3.2.1 60 Hz lag-window
    high half (`lag_h`, `[i16; 10]`).
  - `tables::LPC_LAG_WINDOW_LOW_Q15` — §3.2.1 60 Hz lag-window
    low half (`lag_l`, `[i16; 10]`).
  - `tables::LSF_SEARCH_GRID_COS_Q15` — §3.2.5 `az_lsf()`
    root-search cosine grid (`grid`, `[i16; 61]`).
  - `tables::PITCH_INTERP_FILTER_ANALYSIS_Q15` — §3.7
    1/3-resolution pitch analysis filter (`inter_3`, `[i16; 13]`).
  - `tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15` — §3.7
    1/3-resolution pitch synthesis filter (`inter_3l`, `[i16; 31]`).
  - `tables::GAIN_QUANT_MA_PREDICTOR_Q13` — §3.9 MA gain
    predictor `pred` (`[i16; 4]` ≈ {0.68, 0.58, 0.34, 0.19}).
- New spec-dimension helper constants exposed alongside the table
  module: `M = 10` (LP order), `L_WINDOW = 240` (LP-analysis frame
  length), `GRID_POINTS = 60` (LSF root-search grid resolution).
- `tests/tables_shape.rs` grows from 10 to 26 tests, structurally
  verifying every newly-staged table:
  - Hamming window length matches `L_WINDOW`; peak equals Q15 ≈ 1.0
    (`32767`); every sample is strictly positive.
  - Lag-window pair lengths match `M`; `lag_h` is strictly
    monotonically decreasing.
  - Cosine grid length matches `GRID_POINTS + 1`; endpoints match
    the CSV literals (`32760` / `-32760`); midpoint is exactly `0`;
    grid is strictly decreasing and antisymmetric about the
    midpoint.
  - Pitch analysis / synthesis filters: peak tap is positive and
    equals the maximum-magnitude tap.
  - MA gain predictor matches `[5571, 4751, 2785, 1556]` literally
    and round-trips to {0.68, 0.58, 0.34, 0.19} within one Q13
    quantisation step; vector is monotonically non-increasing.

- Round 173 lands the bit-exact numeric-tables foundation. The
  following constants are compiled at build time by `build.rs` from
  CSVs under `tables/` (byte-for-byte copies of the spec-role-named
  outputs under `docs/audio/g729/tables/`):
  - `tables::HPF_PREPROC_100HZ_B_Q13` — §3.1 / §4.2 100 Hz pre/post
    HPF b-coefficients (`[i16; 3]`).
  - `tables::HPF_PREPROC_100HZ_A_Q13` — §3.1 / §4.2 100 Hz pre/post
    HPF a-coefficients (`[i16; 3]`).
  - `tables::HPF_PREPROC_140HZ_B_Q12` — §3.1 alternate 140 Hz HPF
    b-coefficients (`[i16; 3]`).
  - `tables::HPF_PREPROC_140HZ_A_Q12` — §3.1 alternate 140 Hz HPF
    a-coefficients (`[i16; 3]`).
  - `tables::BIT_ALLOCATION_TABLE8` — §4.1 per-parameter bit
    allocation (`[i16; 13]`, source `bitsno` array).
  - `tables::POW2_TABLE_Q15` — `basic_op::Pow2()` lookup
    (`[i16; 33]`).
  - `tables::LOG2_TABLE_Q15` — `basic_op::Log2()` lookup
    (`[i16; 33]`).
  - `tables::INV_SQRT_TABLE_Q15` — `basic_op::Inv_sqrt()` lookup
    (`[i16; 49]`).
- Companion test `tests/tables_shape.rs` pins each table's shape and
  spot-checks documented values (HPF triples, Pow2 anchor and
  boundary, Log2 monotonic non-decreasing, Inv_sqrt monotonic
  non-increasing, BIT_ALLOCATION_TABLE8 full literal sequence).
- `build.rs` carries the spec clause + Q-format + source-file
  SHA-256 + electronic-attachment ZIP SHA-256 into each generated
  `pub const`'s doc comment, so the in-crate provenance is
  reviewable without leaving `cargo doc`.

### Erased

- Prior master history was force-erased on **2026-05-24** under
  Hat-3 cold enforcement of the workspace clean-room policy
  (`docs/IMPLEMENTOR_ROUND.md`). The retired implementation
  transcribed numeric tables verbatim from an external
  reference-software distribution and described matching its
  behaviour by citing specific source files of that distribution.
  The clean-room policy forbids consulting any external
  implementation's source for any reason, regardless of licensing
  or technical merit.

### Reset

- Crate reduced to a minimal `oxideav_core::register!` stub. Every
  public API returns `Error::NotImplemented`. The crates.io version
  (`0.0.6`) is preserved on the new master to avoid breaking any
  downstream version pins; the published versions on crates.io will
  be yanked by the maintainer.

### Next

- Remaining codebook tables (LSP L1/L2 codebooks, gain GA/GB
  codebooks, MA predictor `fg`, postfilter interpolation
  `tab_hup_s` / `tab_hup_l`, taming `tab_zone`, LSF↔LSP cos/slope
  tables, Annex B DTX/CNG) — staged under
  `docs/audio/g729/tables/`, awaiting per-clause specifier pass.
- Decoder / encoder wire-up against the staged ITU-T G.729
  Recommendation text (spec PDF at `docs/audio/g729/`).
