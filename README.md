# oxideav-g729

Pure-Rust ITU-T **G.729** (CS-ACELP, 8 kbit/s) narrowband speech codec.

> **Status: orphan-rebuild scaffold (reset 2026-05-24).**
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
> The crate currently registers under the `g729` codec id and every
> public API returns `Error::NotImplemented`. It will be
> re-implemented from scratch against the staged ITU-T G.729
> Recommendation text in a future clean-room round.

## License

MIT — see [`LICENSE`](LICENSE).
