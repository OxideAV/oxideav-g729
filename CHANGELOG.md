# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the crate adheres
to [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
