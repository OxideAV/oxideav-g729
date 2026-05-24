# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.7](https://github.com/OxideAV/oxideav-g729/compare/v0.0.6...v0.0.7) - 2026-05-24

### Other

- encoder §3.4 open-loop pitch analysis + §3.7 closed-loop windows
- encoder §3.3 perceptual weighting + §3.5 h(n) impulse response
- complete §4.2.5 decoder output stage + AGC/tilt spec fixes
- spec-exact 240-sample asymmetric LP-analysis window per §3.2.1
- align encoder gain-VQ with spec — γ correction factor + MA-4 predictor

### Other

- New `open_loop_pitch` module implements the encoder's §3.4 open-loop
  pitch analysis, spec-faithful (12 unit tests across the module + the
  encoder search-bounds tests):
  - `WeightedSpeechState` runs the perceptual-weighting filter
    `W(z) = A(z/γ1)/A(z/γ2)` to produce the weighted speech `sw(n)` of
    eq (33), with persistent `s`/`sw` history so the recursion is
    continuous across subframe and frame boundaries. Cross-checked
    against an independent block-mode reference recursion.
  - `open_loop_pitch_search` computes the three-range correlation
    `R(k) = Σ_{n=0..79} sw(n)·sw(n-k)` (eq 34), normalises it by
    `sqrt(Σ sw²(n-k))` (eq 35) over the three disjoint delay bands
    (80…143, 40…79, 20…39), and applies the sub-multiple bias decision
    tree (`R'(t_i) ≥ 0.85·R'(T_op)` favours shorter delays) to pick
    `T_op` without latching onto a pitch multiple.
- Encoder wires `T_op` into the §3.7 closed-loop pitch search: subframe 1
  now searches the spec's six-sample window `[T_op−3, T_op+3]`
  (eq f0018-01) and subframe 2 the ten-sample window
  `[int(T1)−5, int(T1)+4]` (eq f0018-02), both with the spec boundary
  clamps, via the new `PitchAnchor` enum and `subframe1_search_bounds` /
  `subframe2_search_bounds` helpers (unit-tested). Replaces the prior
  ad-hoc ±5 anchor windows.
- The §3.4 weighted-speech path builds the weighting filter from the
  *unquantised* per-subframe LP coefficients (interpolated from the
  unquantised LSP, mirroring the quantised interpolation), and uses the
  §3.3 "flat" gammas (0.94, 0.6) as a documented default until the
  LAR-based flat/tilted γ-adaptation (eqs 28–30) is wired in.
- New `weighting` module implements the encoder's §3.3 perceptual
  weighting filter and §3.5 weighted-synthesis-filter impulse response,
  all spec-faithful (12 unit tests):
  - `bandwidth_expand` — `A(z) → A(z/γ)` (coefficient `a_k → γ^k·a_k`,
    the numerator/denominator form of eq 27).
  - `log_area_ratio` (eq 28), `lsp_min_distance` `d_min` (eq 31),
    `adapt_gamma2` `γ2 = −6·d_min + 1` bounded to [0.4, 0.7] (eq 32),
    `classify_flat` flat/tilted hysteresis decision (eq 30), and
    `weighting_gammas` selecting (γ1, γ2) = (0.94, 0.6) flat or
    (0.98, adapted) tilted.
  - `weighted_synthesis_impulse_response` — `h(n)` of
    `W(z)/Â(z) = A(z/γ1)/[Â(z)·A(z/γ2)]`, computed per §3.5 by filtering
    the `A(z/γ1)` coefficients (zero-extended) through `1/Â(z)` then
    `1/A(z/γ2)`, truncated to one subframe (40 samples). Cross-checked
    against an explicit two-pass reference recursion.
  These are the building blocks for a true analysis-by-synthesis search;
  the §3.7/§3.8 search loops do not yet convolve candidates with `h(n)`.
- Decoder now completes the §4.2.5 output post-processing stage: the
  100 Hz high-pass filter `H_h2(z)` (eq 91, coefficients verbatim from
  the spec — b0 = b2 = 0.93980581, b1 = −1.8795834, a1 = −1.9330735,
  a2 = 0.93589199) followed by the ×2 level restore that undoes the
  encoder's §3.1 ×0.5 preprocessing scale. Applied to both the normal
  CS-ACELP path and the frame-erasure concealment path with a shared,
  continuous biquad memory. This was previously missing entirely — the
  decoder emitted the AGC output directly.
- Decoder AGC (§4.2.4) corrected to the spec's L1-norm gain ratio
  `G = Σ|ŝ(n)| / Σ|sf(n)|` (eq 88), replacing the prior RMS/energy
  ratio. Per-sample smoothing `g(n) = 0.85·g(n-1) + 0.15·G` (eq 90)
  unchanged.
- Decoder tilt-compensation filter (§4.2.3) now applies the eq-86
  leading scale `1/g_t` with `g_t = 1 − |k1'|`, so the product filter
  `H_f(z)·H_t(z)` has the spec's near-unit gain; the impulse-response
  truncation length for `k1'` corrected to L = 20 per eq 87 (was 22).
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
