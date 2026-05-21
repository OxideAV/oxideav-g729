# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Other

- Encoder LP-analysis windowing is now spec-exact per ITU-T G.729
  §3.2.1 / eq (3): a 240-sample asymmetric window (half-Hamming
  `0.54 − 0.46·cos(2π·n/399)` for n = 0…199, then quarter-cosine
  `cos(2π·(n−200)/159)` for n = 200…239) is applied to a buffer of
  120 past samples + the 80-sample current frame + 40 samples of
  look-ahead from the following frame. The autocorrelation `r(0)`
  lower bound (1.0), white-noise correction factor (×1.0001), and
  60 Hz bandwidth-expansion lag window from eqs (6, 7) are all
  applied verbatim. Replaces the previous 80-sample symmetric
  Hamming approximation.
- Encoder now exposes the spec's 5-ms extra algorithmic delay: 40
  samples of look-ahead PCM must be queued before a frame can be
  released. `flush()` drains any held samples by encoding against
  trailing zeros.
- `tests/annex_b.rs::vad_classifies_silence_vs_speech` rewritten to
  send silence + tone as a single contiguous stream (then partition
  by frame index) so the section boundary respects the encoder's
  new look-ahead latency.
- Encoder gain-VQ now quantises the spec's correction factor
  `γ = g_c / g'_c` (eq 72) instead of raw `g_c`, mirroring the
  decoder's `g_c = γ̂ · g'_c` reconstruction (eq 74).
- Encoder maintains a 4-deep MA-4 gain-prediction history `Û^(k)`
  initialised to `-14 dB` per ITU-T G.729 Table 9; history is rolled
  on each subframe with `Û^(m) = 20·log10(γ_q)` per eq (70) in
  lockstep with the decoder.
- Excitation history written back to the encoder's adaptive-codebook
  buffer now uses the quantised `γ_q · g'_c`, matching what the
  decoder reconstructs on the receive side.
- New integration test `steady_tone_gain_predictor_converges`
  validates that the gain-predictor state stays stable across
  5 seconds of input without runaway feedback or collapse.
- README: correct the gain-VQ table dimensions (spec says `gbk1` is
  8×2 and `gbk2` is 16×2 per §5.2 Table 12, not 8×8 + 16×8); call out
  the MA-4 gain predictor as spec-exact going forward.

## [0.0.6](https://github.com/OxideAV/oxideav-g729/compare/v0.0.5...v0.0.6) - 2026-05-06

### Other

- drop dead `linkme` dep
- auto-register via oxideav_core::register! macro (linkme distributed slice)

## [0.0.5](https://github.com/OxideAV/oxideav-g729/compare/v0.0.4...v0.0.5) - 2026-05-03

### Other

- replace never-match regex with semver_check = false
- cargo fmt: fix rustfmt --check CI gate
- migrate to centralized OxideAV/.github reusable workflows
- transcribe LSP codebooks verbatim from ITU TAB_LD8K.C, fix L3 high-half indexing
- adopt slim VideoFrame/AudioFrame shape
- pin release-plz to patch-only bumps

### Other

- LSP codebooks now transcribed verbatim from ITU reference
  `TAB_LD8K.C`: `LSPCB1_Q13` (128×10), `LSPCB2_Q13` (32×10), `FG_Q15`,
  `FG_SUM_Q15`, `FG_SUM_INV_Q12`. `synth_lspcb1_row` / `synth_lspcb2_row`
  procedural row synthesis dropped; `FREQ_PREV_RESET_Q13` reset vector
  added.
- Fix L3 column indexing in `lpc::decode_lsp` and the encoder split-VQ
  search: `cb2_hi` is a full 10-wide row, so the high-half contribution
  reads cols `M_HALF..M`, not `0..M_HALF`. Matches `Lsp_get_quant` in
  `LSPGETQ.C` (`lspcb2[code2][j]` for `j ∈ [NC, M)`).

## [0.0.4](https://github.com/OxideAV/oxideav-g729/compare/v0.0.3...v0.0.4) - 2026-04-25

### Other

- drop oxideav-codec/oxideav-container shims, import from oxideav-core
- frame-erasure concealment per §4.4
- per-sample AGC smoothing per §4.2.4
- adaptive spectral-tilt compensation per §4.2.3
- correct short-term post-filter gamma assignment per §4.2.2
- implement long-term (pitch) post-filter per §4.2.1
- proper MA-4 gain prediction per §3.9.2
- drop Cargo.lock — this crate is a library
- bump oxideav-core / oxideav-codec dep examples to "0.1"
- bump to oxideav-core 0.1.1 + codec 0.1.1
- migrate register() to CodecInfo builder
- bump oxideav-core + oxideav-codec deps to "0.1"
- document current G.729 capabilities and annex support
