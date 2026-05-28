# oxideav-g729

Pure-Rust ITU-T **G.729** (CS-ACELP, 8 kbit/s) narrowband speech codec.

> **Status: orphan-rebuild scaffold (reset 2026-05-24)** + early
> data-only foundation landed 2026-05-29 (round 173).
>
> The previous implementation was retired under the OxideAV workspace
> clean-room policy. Several of its source modules transcribed numeric
> tables verbatim from an external reference-software distribution and
> documented matching that distribution's behaviour by citing specific
> source files of it. The clean-room policy forbids consulting any
> external implementation's source for any reason — regardless of the
> distribution's licensing or the technical merit of the values — so
> the provenance of those tables could not be defended and the history
> was force-erased.
>
> The crate currently registers under the `g729` codec id with a
> no-op `register()` and every higher-level codec API returns
> `Error::NotImplemented`. The clean-room re-implementation against
> the staged ITU-T G.729 Recommendation text + the in-progress
> clean-room trace doc is pending.

## What landed in round 173

A `pub mod tables` foundation exposing a small subset of the
G.729 bit-exact numeric tables as `pub const [i16; N]` arrays:

- §3.1 / §4.2 pre-/post-processing high-pass IIR filter
  coefficients (100 Hz Q13 and 140 Hz Q12 variants).
- §4.1 Table 8 bit-allocation raw extract.
- The three basic-op fixed-point math LUTs used by `Pow2`, `Log2`,
  and `Inv_sqrt`.

All constants are emitted at build time from the spec-role-named CSV
workspace at `docs/audio/g729/tables/` (carried verbatim in this
crate's own `tables/` directory for hermetic publishing). The CSVs
themselves are produced by `docs/audio/g729/tables/extract.py`,
which reads only data-bearing `tab_*.c` files (numeric arrays + their
inline comments) from the ITU electronic-attachment package; it does
not read any algorithmic source. See
`docs/audio/g729/tables/README.md` for the full provenance chain.

Structural shape tests (15 of them) confirm each constant's length
matches its `.meta` sidecar, plus a handful of spec-implied
invariants — high-pass numerator symmetry, DC-cancellation, LUT
monotonicity — that catch any silent regression during table
copying. No test compares numeric values against any third-party
source.

## License

MIT — see [`LICENSE`](LICENSE).
