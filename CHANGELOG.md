# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
