# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the crate adheres
to [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `tables` module exposing a small bit-exact subset of the G.729
  numeric tables as `pub const [i16; N]` arrays, emitted at build
  time by `build.rs` from the spec-role-named CSV workspace at
  `docs/audio/g729/tables/` (carried verbatim under this crate's
  own `tables/` directory for hermetic publishing). Initial scope:
  the §3.1 / §4.2 pre-/post-processing high-pass IIR filter
  coefficients (100 Hz Q13 and 140 Hz Q12 variants), the §4.1
  Table 8 bit-allocation raw extract, and the three basic-op
  fixed-point math LUTs (`Pow2`, `Log2`, `Inv_sqrt`).
- Shape tests confirming each constant's length matches its
  `.meta` sidecar plus spec-implied invariants (high-pass
  numerator symmetry, DC cancellation, LUT monotonicity).

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

- Clean-room re-implementation against the staged ITU-T G.729
  Recommendation text (numeric tables and decode behaviour read
  only from the standard) in a future round.
