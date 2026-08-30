# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the crate adheres
to [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Round 452 — frame-0 startup identified; every clean vector's first
  subframe byte-exact.** Inverting the §4.2.5 output high-pass on the
  frame-0 `.PST` samples of all six clean vectors (both corpora agree
  there sample-for-sample) recovers the reference's frame-0 `ŝ`; the
  §3.2.5/§4.1.6 previous-LSP start-up memory is corpus-pinned to
  `30000 26000 21000 15000 8000 0 −8000 −15000 −21000 −26000` (Q15) —
  *not* the flat `cos(iπ/11)` grid (kept as
  `STARTUP_LSP_COS_GRID_Q15`; the printed Table 9 row is
  domain-invalid either way). With the correct startup the r419
  §4.2/§3.9.1 pins were re-swept: synthesis lands by **rounding**;
  `1/g_f` scales the synthesis input; the §4.2.5 ×2 folds into the
  eq (89) AGC product (Q1 high-pass grid, exact 48-bit feedback
  products landed on Q15); `ĝ_c` lands on Q1 by **truncation**; the
  eq (75) excitation landing rounds ties toward zero. Measured
  (base-corpus exact%, r419 → r452): ALGTHM 1.89 → 4.04, FIXED
  19.09 → 34.45, LSP 3.49 → 4.01, PITCH 1.84 → 1.98, SPEECH
  14.02 → 21.80 (clean frames 162 → 285), TAME 0.78 → 0.81, PARITY
  14.38 → 26.59, ERASURE 13.94 → 20.43. Floors raised; six new
  doc(hidden) `PfLatitudeFx`/`GainGridFx` sweep hooks carry the
  measured rulings (`agc_lag` reproduces the reference's onset AGC
  trajectory exactly and is the recorded lead; `gf_exact_div` is
  refuted).
- **Round 452 — Annex B comfort noise gains its spectral envelope.**
  New `sid_lsf` module implements the §B.4.2.2 SID-LSF dequantizer
  from the newly staged tables: the eq (B.18) blended second predictor
  (staged per-mode column sums asserted to be the blend of the full
  coder's rows), the 32-address first-stage subset into L1, and the
  2 × 16 full-VQ second stage into L2, then the unchanged §3.2.4
  machinery. `decode_chain::decode_cng_frame_to_postfiltered` drives
  the §B.4.4 excitation through the interpolated SID-LSP synthesis
  filter and the §4.2 cascade **sharing the active chain's state**
  (interpolator, eq (40) excitation history, filter memories);
  `AnnexBOutput::ComfortNoise` now carries post-processed PCM and
  `AnnexBOutput::ErasedActive` (replacing the placeholder variant)
  carries §4.4-concealed PCM per §B.4.5. Measured vs the `g729b`
  reference `.out` (new `annex_b_stream_tracks_reference_output`,
  floors pinned): active-frame corr 0.9914/0.9886/0.9831/0.9716/
  0.6589/0.9965 on tstseq1–6, CN RMS envelope ratio 0.72–1.06.
- **Round 452 — the r452 docs stagings compiled**: the §3.2.4/Table 9
  MA-predictor reset vector (`trunc(iπ/11·2¹³)`, truncated not
  rounded) now seeds `startup_history_q13()`, and the seven Annex B
  SID-LSF/gain table pairs enter `build.rs` with shape/domain/
  monotonicity tests.

- Round 438 moves the **encoder's §3.1–§3.2.3 front end onto the
  clause-5 Word16/Word32 operator grid** (`fx::analysis`,
  `fx::levinson`, `fx::lp_to_lsp`) and wires it into the production
  `FrameEncoder`:
  - `fx::analysis` — eq (4) windowing via `mult_r`, the eq (5)
    autocorrelation in a genuine Word32 accumulator with the
    **overflow-rescale protocol** (saturation → right-shift the
    windowed signal, retry; step black-box-pinned at 2), the
    clause-3.2.1 `r(0) ≥ 1` guard, `norm_l` normalisation to DPF
    pairs, and the eqs (6)/(7) lag window as `mpy_32` split-double
    products over the staged tables. Where the float emulation sums
    exactly and normalises afterwards, the Word32 pipeline must
    down-scale the *input* on overflow — a different truncation
    pattern on exactly the loud frames where index agreement dropped.
  - `fx::analysis::PreprocessorFx` — the §3.1 eq (1) 140 Hz
    high-pass (÷2 folded into the staged Q12 numerator) with the
    feedback kept as the unrounded Word32 recursion value on the Q16
    grid (corpus-pinned: Q13 feedback storage costs ~12/16 locked
    points of LSP/SPEECH; a Q28 DPF coefficient variant quantised
    from the printed eq (1) decimals costs ~8/10).
  - `fx::levinson` — eqs (8)/(9) entirely in the double-precision
    (hi, lo) format: Q27 LP coefficients, Q31 reflections, `div_32`
    against a per-order renormalised energy with the numerator
    re-referred onto the stored-E grid *before* the division
    (corpus-pinned; the post-division quotient shift costs ~6 points
    of ALGTHM). Ill-conditioned input reports failure for the
    clause-3.2.3 keep-previous fallback.
  - `fx::lp_to_lsp` — the §3.2.3 Chebyshev search on the fixed grid:
    Q11 half-polynomial coefficients (**truncating** Q12→Q11 landing,
    corpus-pinned), Clenshaw in Word32 Q23, the staged 60-interval
    cosine grid with four bisections, and a `div_s` linear
    interpolation (rounding step).
  - **Measured** (`tests/fx_front_end_conformance.rs`, reference-
    locked all-four-LSP-indices agreement, float → fx chain):
    ALGTHM 82.9 → 94.3%, FIXED 97.5 → 97.5%, LSP 77.7 → 85.8%,
    PITCH 92.8 → 93.6%, SPEECH 81.2 → 91.7%, **TAME 80.5 → 99.2%**
    (corpus 83.1 → 90.7%). End-to-end own-history agreement:
    **TAME L0 62.5 → 96.1%** (T1±2 75.0 → 82.0%), SPEECH L0/L1
    85.5/50.7 → 87.0/53.9%, PITCH L1 32.9 → 35.6%, LSP L1
    47.2 → 48.5%.
  - **TAME's jump retires the round-390 negative result**: the 24
    locked-history `L0` flips attributed to an "unstaged structural
    element of the reference's final mode compare" (a documented
    docs-gap) were front-end ω divergence all along — the genuine
    fixed-point §3.1/§3.2.1 truncation pattern lands ω on the
    reference's grid and the printed eq (21) mode compare then picks
    the reference's predictor. No mode-compare docs-gap remains.
- The decoder's Table 9 start-up LSP memory is pinned to the **exact
  Q15 cosine grid** `round(cos(iπ/11)·2^15)` per the staged erratum
  resolution (`docs/audio/g729/table9-initialisation.md`: the printed
  `arccos(iπ/11)` row is domain-invalid for `i ≥ 4`; §3.2.3's
  `q_i = cos(ω_i)` applied to the adjacent `ω̂_i` row fixes the
  intended value), replacing the 64-segment cosine-table image (up
  to 8 Q15 LSB off). Measured: the frame-0 sub-1 interpolated Q12 LP
  moves 1–4 LSB but every conformance vector's PCM is byte-identical
  — the frame-0 startup divergence is *not* the startup LSP grid;
  the remaining schedule question is confirmed out of the
  Recommendation's scope by the same doc.

- Round 419 lands the **fixed-point §4.2 post-processing cascade**
  (`fx::postfilter`), completing the whole decode path on the clause-5
  Word16/Word32 operator grid. The cascade follows the clause-4.2
  signal order as printed — `ŝ` → `Â(z/γ_n)` → residual `r̂` →
  long-term postfilter `H_p(z)` **applied to the residual** →
  synthesis `1/[g_f·Â(z/γ_d)]` → tilt `H_t(z)` → adaptive gain
  control against `ŝ` → 100 Hz output high-pass + ×2 — which differs
  structurally from the float cascade's historical `H_p`-on-speech
  arrangement. The §4.2.1 two-pass search runs its correlations on a
  norm-scaled Word16 residual window (Word32-safe by construction),
  covers both sides of the integer anchor at 1/8 resolution, applies
  the clause's longer-filter replacement rule, and lands `g_l` on Q15
  by normalised division; `1/g_f` / `1/g_t` go through a shared
  normalised-reciprocal helper; the §4.2.4 gain runs on Q12 with the
  Q15 weight pair `27853/4915`; the §4.2.5 recursion keeps its
  feedback on the wide Word32 grid with the ×2 folded into the final
  rounding. Unpinned operator-schedule choices are exposed as a
  `#[doc(hidden)]` latitude struct (plus a γ̂-grid hook on the gain
  decoder) and were swept black-box against the corpus — the shipped
  defaults are the sweep optimum; the round-410 uniform-2^13 γ̂ grid
  measurably beats every per-codebook split (best challenger
  2–4× worse per-subframe energy delta on all 12 clean vectors).
- `tests/fx_full_conformance.rs` — the byte-exactness ratchet for the
  full fixed-point chain: per-vector correlation / RMS / exact-share /
  max|Δ| plus **first-diverging-sample** (frame/subframe/offset,
  ours-vs-reference), longest byte-exact run, and clean-frame counts,
  with pinned floors; an env-gated stage-by-stage trace dump
  (`G729_FX_TRACE`) and `#[doc(hidden)]` trace instrumentation on the
  cascade (`process_subframe_traced`) anchor the bisection workflow.
  Measured full-fx: corr 0.9927–0.9999 on all 12 clean vectors
  (FIXED 0.9855/0.9918), PARITY 0.998, ERASURE 0.92/0.89, OVERFLOW
  0.81 base.
- **§4.2.1 over-unity gain guard, pinned black-box** — when the raw
  eq (83) ratio exceeds 2 (`num > 2·den`, the ceiling of a Q14 gain
  grid) the long-term postfilter behaves as *disabled*, not clamped
  to 1. The corpus separates the readings cleanly: clamping keeps
  voiced material but wrecks the onset-heavy FIXED vectors
  (0.9502/0.9756 vs 0.9855/0.9918 with the guard); disabling every
  over-unity case instead costs SPEECH/PITCH (0.9953/0.9926). A
  per-subframe oracle probe confirms passthrough beats the clamped
  filter on 42 of 45 FIXED enables while SPEECH splits evenly.
- **Two more black-box §4.2 pins, ratcheted by whole-frame
  byte-exactness** — (1) the `1/[g_f·Â(z/γ_d)]` synthesis output
  lands on Q0 by *truncation* (`extract_h`), not rounding: the
  change lifts SPEECH's longest byte-exact run from 77 to **9155
  samples** (114 consecutive frames) and its fully-exact frame count
  from 0 to 159/3750 (g729a: 166), and the reference startup tails
  read back the same way (sub-unity decays land as zeros); rounding
  and toward-zero variants both collapse the clean-frame count.
  (2) a residual with mean square at most one Q0 LSB (`Σr̂² ≤ 40`)
  never enables the long-term postfilter — the eq (82) statistic
  over bare quantisation noise otherwise flickers the filter on
  during silence; the threshold is insensitive from 40 to 2560 and
  retires the silence-cluster divergence events (SPEECH's event
  census drops from 11 root events to 9, none LT-related).

### Changed

- **The registry decoder (`codec::G729Decoder`, the dual-API
  `decoder::make_decoder`) now decodes through the fixed-point chain**
  (`fx` §4.1 + §4.2) instead of the float `decode_chain` path. On the
  ITU conformance corpus this moves whole-vector correlation against
  the reference `.PST` from the float chain's long-vector drift
  (SPEECH ≈ 0.11) to 0.9855–0.9999, with byte-exact stretches over a
  hundred frames long. `reset()` still restores the clause-4.3
  start-up state byte-identically; the wire framing is unchanged.
- The in-flux internal modules (the per-clause §3/§4 processing-chain
  stages, the `fx` fixed-point operator grid, the numeric-table
  accessors, and the Annex B stream-decoder state machines) are now
  `#[doc(hidden)]`, shrinking the documented/semver-checked API to the
  stable surface: the registry codec (`codec`, root `register`), the
  dual-API `decoder`/`encoder` endpoints, the ITU serial-format
  helpers (`serial`), and the Annex B serial-framing parse surface.
  No signatures changed; the hidden items remain `pub` and callable.

### Added

- Round 410 lands the **fixed-point §4.1 decode chain** (`fx` module
  tree): the clause-5 Table 10/11 Word16/Word32 basic-operator set,
  the Table 15 `log2`/`pow2`/`inv_sqrt` over the staged 33/33/49-entry
  data tables plus the double-precision (hi, lo) helpers, the Q13
  LSF → Q15 LSP → Q12 LP parameter decode (sharing the round-390
  encoder-search integer arithmetic via
  `lsp_quantize::reconstruct_lsf_q13`; two newly compiled conversion
  tables `lsf-lsp-cos-table2-Q15` / `lsf-lsp-cos-slope-Q19`), the Q10
  gain predictor with the eq (71) exponent identity
  `(log2(10)/20)·(10·log10 2) = 1/2`, the eq (40) adaptive vector on
  the corrected fraction fold, the Q13/Q14 codevector + sharpening,
  the Q15-accumulator excitation and Q12 synthesis with the specified
  saturation, §4.1.2 parity T1-substitution, §4.4 concealment
  primitives, and the black-box-pinned §4.1.6 overflow-rescale
  protocol (excitation + synthesis memory ÷4 on Word32 overflow).
  Measured through the fx-§4.1 + float-§4.2 hybrid against the ITU
  `.PST` references (`tests/fx_conformance.rs`): corr 0.995–0.9998 on
  all 12 clean vectors, PARITY 0.9987, ERASURE 0.91/0.94, OVERFLOW
  0.78 (base corpus), RMS 0.96–1.07.

### Fixed

- Round 410 **fixes the eq (40) fraction fold** — the two-branch
  `b30` interpolation evaluates the past excitation at
  `n − k + t/3`, so the effective delay of the `(k, t)` pair is
  `k − t/3` and the transmitted fraction folds negated with an upward
  borrow (`frac = +1 → (t0 + 1, 2)`). The previous `T = k + t/3`
  fold read every fractional subframe 2/3 of a sample off; float
  decode-chain correlation against the `.PST` references moved
  SPEECH 0.11 → 0.93, PITCH 0.54 → 0.93, LSP 0.65 → 0.87/0.97
  (floors ratcheted).

- Round 410 **fixes the eq (74) γ̂ reconstruction grid** — the GA/GB
  codebook `γ` columns sum to the fixed-codebook gain correction factor
  on a 2^13 grid (not 2^12): the Recommendation never states the column
  scaling (Table 12 lists only dimensions), and the black-box divisor
  sweep over the ITU conformance corpus lands at 6.9–10× / ≈ 1.0× /
  0.15–0.2× whole-vector RMS ratio for 2^12 / 2^13 / 2^14 (the scale
  error compounds through the eq (69)/(72) MA feedback to exactly
  `2^(1+Σb_i)` ≈ 6.9×). This retires the historical "§3.9 gain-
  saturation gap": decode RMS ratios collapse from ≈ 7–10× (TAME 18.5×)
  to 0.97–1.36 (TAME 2.85). The encoder GA preselection ranks on the
  same grid. `tests/pcm_conformance.rs` now pins the RMS ratio inside a
  `[0.7, 3.2)` window and tracks per-vector correlation / max |Δ| /
  bit-exact share with pinned correlation floors.

### Added

- Round 390 lands the **Annex A (reduced-complexity) encoder**
  (`annex_a_encoder`): the §A.3 analysis chain assembled from the
  unchanged main-body stages plus the five prose-pinned
  simplifications (eqs (A.1)–(A.10)): the fixed `γ = 0.75`
  quantized-LP weighting (`W(z)/Â(z) = 1/Â(z/γ)` — single all-pole
  impulse response + target filter), the eq (A.2)/(A.3) low-pass
  weighted speech `1/[Â(z/γ)(1 − 0.7z⁻¹)]` for the open-loop
  correlation, the §A.3.4 **fast open-loop pitch** (eq (A.4)
  even-sample decimated correlation, eq (A.5) normalisation,
  even-delays-first + ±1 refinement in `[80, 143]`, favour-lower
  cascade), the §A.3.7 **fast adaptive-codebook search** (eq (A.7)
  backward-filtered target × unfiltered past excitation, ±⅓ fractions
  via the eq (A.8) `b30` interpolation), and the §A.3.10
  filtering-free eq (A.10) memory update. §A.3.2.5 interpolates only
  the quantized LSP track. Measured against the G.729A reference
  `.BIT` corpus (`annex_a_encoder_parameter_floors`): L0 62.5–94.6%,
  L1 33.8–100%, T1±2 57.8–86.4% (floors pinned per vector); the whole
  3750-frame g729a SPEECH vector encoded by us decodes finitely
  through both the main-body and Annex A decoders (§A.1 bit-stream
  interoperability). Documented gap: the §A.3.8.1 depth-first ACELP
  pulse schedule is prose-unpinned (Annex A defers it to the barred
  reference C); the main-body focused search is used instead.
- Round 390 lands the **Annex A (reduced-complexity) decoder**
  (`annex_a`): the §A.4 decoder-side deltas over the unchanged §4.1
  parameter chain + §4.1.6 synthesis, transcribed from the in-repo
  Recommendation Annex A prose + equation rasters (eqs (A.11)–(A.15)):
  the §A.4.2.1 **integer-delay** long-term postfilter (search
  `[T_cl − 3, T_cl + 3]` on the *current subframe's* transmitted
  delay, `T_cl ≤ 140`), the §A.4.2.2 short-term postfilter without
  `1/g_f`, the §A.4.2.3 tilt without `1/g_t` (eq (A.14) `k′_1` from
  the length-22 truncated impulse response; `γ_t = 0.8` iff
  `k′_1 < 0`, else disabled), the Annex A cascade order (numerator →
  tilt → synthesis), and the eq (A.15) **energy-ratio** AGC
  (`√(Σŝ²/Σsf²)`, smoothing 0.9/0.1). `AnnexADecoder` +
  `AnnexAPostfilterCascade` (with `FrameDecoder::
  decode_parameters_with_speech` exposing the frame + synthesis pair).
  Validated black-box against the staged `g729a` corpus
  (`tests/annex_a_conformance.rs`): first-subframe deviations 2.4–4.5
  PCM units across all six clean vectors (same 8-unit band as the base
  harness), whole-vector RMS ratios 7.1–10.6 under the same bounded
  §3.9-gain-gap ceiling, and per-frame gain-normalised shape distance
  at parity with the base cascade (3.577 vs 3.552 summed; the §3.9
  gain gap dominates both at current fidelity, so the harness pins
  parity + a regression ceiling). Encoder-side §A.3 deltas are not
  implemented: the §A.3.8.1 depth-first ACELP pulse schedule is
  prose-unpinned (deferred to the barred reference C — docs gap).
- Round 390 lands the **fixed-point Q13 §3.2.4 quantiser search**
  (`lsp_quantize::search_lsp_indices_q13`), implementing the L0
  MA-predictor mode selection exactly as described by the newly-staged
  clean-room algorithm doc `docs/audio/g729/l0-mode-selection.md`
  (docs commit `b9e48a4`): eq (23) targets via the Q12 `fg_sum_inv`
  reciprocal, integer split-stage searches, eq (20) reconstruction on
  the Q13 grid, and the ω-domain eq (21) argmin over the two
  predictors — with the doc's three unpinned fixed-point latitude
  points exposed as `FxLatitude` and pinned black-box against the
  conformance corpus (320+ configurations swept):
  - **Pinned latitude** (`FxLatitude::default`): the eq (20)/(23)
    `Q28 → Q13` prediction shift *rounds* (`+2^14`), the eq (23)
    reciprocal shift and the rearrangement midpoint *truncate*, the
    eq (21) per-term combine is the exact wide product, weights on a
    `2^11` grid, mode ties hold mode 0. `LspQuantizer::quantize_lsf`
    now runs this search over an integer MA history
    (`startup_history_q13` / `advance_history_q13`); the returned
    reconstruction still delegates to `LspReconstructor`, keeping the
    encode → decode index contract intact.
  - **Measured** (reference-locked all-four-indices agreement):
    LSP 77.5 → 77.7%, SPEECH 80.9 → 81.2%, ALGTHM 82.9%, FIXED 97.5%,
    PITCH 92.8%, TAME 80.5% (corpus 82.9 → 83.1%); end-to-end FIXED
    L0 88.3 → 90.8% and L1 64.2 → 70.0%, PITCH L1 31.3 → 32.9%, LSP
    L1 48.1 → 47.2% (error-propagation reshuffle; its locked number
    rises). Floors ratcheted in both conformance harnesses; the
    locked harness now drives the fixed-point search directly.
  - **Negative results, measured precisely** (the round's main
    finding): TAME's 24 residual L0 flips are *not* fixed-point
    rounding near-ties. Every swept rounding/truncation/tie-break
    combination leaves them intact; the per-term truncating combine
    shapes (`HiWord`, `Mpy32x16`) never beat the exact product; the
    alternative mode totals (`StageSum`, `ResidualFull`, the
    single-folded `Product`) and the `Lsp_pre_select`-style
    pre-selection metrics (`PreTarget[W]`, `PreOmega[W]`) each
    collapse corpus agreement to 34–65%. Dissection shows the flips
    sit in a self-sustaining 5-frame reference index cycle on the
    stationary TAME input — (0,72,12,5) → (0,80,8,5) → (0,80,29,5) ×2
    → (1,80,29,5) — where the reference emits mode 1 against a
    *systematic* ~2.1% ω-domain margin with identical stage indices
    under both modes; flipping it would need ~100 Q13-LSB input
    perturbations, three orders beyond any rounding latitude. The
    staged doc's round-vs-truncate/tie-break hypothesis is therefore
    measured-insufficient for these frames — an unstaged structural
    element of the reference's final mode compare remains (docs gap,
    reported upstream).
- Round 388 lands the **taming procedure**, the **residual-domain
  §3.2.4 split-stage metric**, and the **registry codec surface**:
  - **Taming procedure** (`taming`, from the newly-staged clean-room
    algorithm doc `docs/audio/g729/taming-procedure.md`) — the
    per-zone worst-case excitation-error state over the compiled
    `tab_zone` partition (`TAMING_ZONE_TABLE`, 153 entries, four
    zones split at 40/80/120), the per-subframe test (`tameflag`) and
    update (`E ← 1 + ĝ_p²·max(E_spanned)`), and the `GPCLIP = 15564/2^14`
    ceiling enforced inside the §3.9.2 gain search
    (`quantize_gains(…, tame)`), wired into `FrameEncoder`. The new
    `reference_taming_fingerprint` conformance test pins the doc's
    unpinned threshold black-box: replaying the recursion over the
    reference's own decoded `(delay, ĝ_p)` streams, the reference
    keeps choosing `ĝ_p > GPCLIP` up to a simulated accumulator of
    ≈ 18 186 (TAME) and never crosses 60 000 anywhere — the reference
    **never tames on the staged corpus**, and neither do we
    (behaviourally identical; a forced-lower-threshold sweep degrades
    TAME GB agreement). Retires the r385 hypothesis that missing
    taming caused TAME's end-to-end L0 gap.
  - **§3.2.4 split-stage searches on the residual-domain weighted
    MSE** — black-box metric-domain discovery: the L2/L3 stages of
    the bit-exact coder behave as `Σ w_i·(l_i − l̂_i)²` (no
    `(1 − ΣP̂)²` folding of the printed ω-domain eq (21)).
    Reference-locked all-four-indices agreement: ALGTHM 71.4 → 82.9%,
    FIXED 91.7 → 97.5%, LSP 74.9 → 77.5%, PITCH 88.6 → 92.8%, SPEECH
    78.4 → 80.9%, TAME flat; no vector degrades. L1 stays unweighted
    and mode selection stays on the printed eq (21) (both probed —
    each alternative collapses agreement). TAME's residual misses are
    now 24/25 pure `L0` mode flips. Also lands the Word32-normalised
    autocorrelation grid + truncated lag-window products and Q31/Q27
    Levinson intermediate rounding (flat-to-marginal, kept for
    clause-2.5 grid consistency). New pooled GA/GB gain-codeword
    agreement columns in the conformance harness; locked floors
    ratcheted to 78/93/73/89/77/74.
  - **Registry codec surface** (`codec`) — `G729Decoder` /
    `G729Encoder` implementing the `oxideav_core` traits over the raw
    10-octet wire framing (80 Table-8 bits, MSB-first octets;
    multi-frame packets), a real `register` hook (id `"g729"`,
    decode + encode, lossy, 8 kHz mono S16), and the dual-API
    endpoints `decoder::make_decoder` / `encoder::make_encoder`. The
    scaffold `Error::NotImplemented` type is retired. New
    `codec_registry` conformance test: all 8 100 active corpus frames
    convert losslessly from the ITU serial format to the wire format
    and decode through the registry `Decoder` sample-identical to the
    `decode_chain` path; `reset()` verified to restore the clause-4.3
    start-up state byte-identically.

- Round 385 moves the **encoder LSF chain onto the reference's
  fixed-point grid** — clause 2.5 makes the 16-bit arithmetic the
  conformance ground truth, and the §3.2.4 quantiser's razor-thin
  nearest-neighbour decisions measurably flip off it:
  - **§3.1 pre-processing on the 16-bit grid** — `process_sample`
    rounds its eq (1) output to the saturated 16-bit PCM grid
    (round-half-up) while the IIR feedback keeps the unrounded
    recursion value (the double-precision memory of a fixed-point
    pole/zero filter).
  - **§3.2.1 fixed-point windowing** — eq (4) evaluated as
    `⌊(s·w + 2^14)·2^−15⌋`, putting the windowed speech on the 16-bit
    grid and making the f64 eq (5) autocorrelation an exact 64-bit
    integer accumulation.
  - **Q12 LP coefficients** — the Levinson output is rounded to Q12
    before the §3.2.3 root search (the reference finds its LSP roots
    on the Q12-rounded `A(z)`; on TAME this single step takes the
    reference-locked all-indices agreement from 0.8% to 39.8%).
  - **eq (18) via the fixed-point arccos lookup** — new
    `lsf_conversion` module compiling the staged
    `lsf-lsp-cos-table-Q15` (65-sample `cos(ω)` grid) +
    `lsf-lsp-acos-slope-Q12` (per-segment slopes) CSVs;
    `acos_q15_to_freq_q15` / `acos_q15_to_lsf_q13` / `lsp_to_lsf_q13`,
    with unit tests pinning monotonicity, exact segment boundaries,
    endpoints, clamping, and the table's measured intrinsic error
    profile. `LspQuantizer` gains the LSF-domain `quantize_lsf` entry
    point; the §3.2.4 search core is extracted as the pure
    `search_lsp_indices` (explicit MA history).
  - **eq (7) white-noise correction at measured-effective unity** —
    black-box sweep of the `r(0)` factor over the whole corpus (the
    prose 1.0001, 1+2^−13 … 1+2^−19, 1.0) shows the literal 1.0001 at
    float precision over-inflates `r(0)` vs the 16-bit reference and
    was the single largest remaining LSF divergence; unity maximises
    corpus-wide agreement (locked ALGTHM 51.4→71.4%, LSP 47.4→74.9%,
    SPEECH 61.8→78.4%, TAME 39.8→80.5%; end-to-end TAME L1
    62.5→100%).
  - **Reference-locked conformance floors** — new
    `encoder_conformance` test driving the quantiser MA history with
    the reference's own transmitted indices (the exact state the
    reference encoder had), separating per-frame front-end fidelity
    from propagation: locked all-four-LSP-indices agreement ALGTHM
    71.4% / FIXED 91.7% / LSP 74.9% / PITCH 88.6% / SPEECH 78.4% /
    TAME 80.5%. End-to-end floors re-pinned per vector (L1 exact up
    to 100% on TAME; T1±2 TAME 56→79%).
  - Documented residuals (black-box-probed and excluded: search
    structure, window timing, root precision, F1/F2 quantisation,
    autocorrelation down-scaling, weight thresholds/caps): the
    fixed-point eq (21)/(22) mode-selection behaviour on extreme
    spectra (TAME prefers predictor 0 against a 3–10% float margin)
    and the DPF Levinson/normalisation precision chain (locked L1
    mismatch margins median ratio 1.21).
- Round 382 builds the **entire clause-3 encoder analysis chain** —
  thirteen encoder stages, spec-cited from the Recommendation prose and
  equation rasters, closing with a working `.IN` → `.BIT` path:
  - **§3.1 pre-processing** (`preprocess`) — the eq (1) 140 Hz
    pole/zero high-pass with the ÷2 scaling folded into the numerator,
    from the compiled Q12 tables.
  - **§3.2.1 LP-analysis front-end** (`lp_analysis`) — eq (4)
    windowing by the compiled eq (3) 240-sample asymmetric window,
    eq (5) autocorrelation (f64, `r(0) ≥ 1` guard), and the eq (6)/(7)
    60 Hz lag window rebuilt from the split-double Q15 pairs + the
    1.0001 white-noise correction.
  - **§3.2.2 Levinson-Durbin** (`levinson`) — eqs (8)/(9) with
    reflection-coefficient by-product and residual energy; validated by
    AR-recovery and the normal-equations residual identity.
  - **§3.2.3 LP→LSP** (`lp_to_lsp`) — sum/difference half-polynomials
    via the (1 ± z⁻¹) divide-out recursion, Chebyshev/Clenshaw
    evaluation, the 60-interval cosine-grid sign-change walk with 4×
    bisection + linear interpolation; full round-trip against the
    decode-side `lsp_to_lp` and the flat-predictor `cos(iπ/11)` set.
  - **§3.2.4 LSP quantiser search** (`lsp_quantize`) — eq (18) arccos,
    eq (22) adaptive weights (+×1.2 on w₅/w₆), eq (23) per-mode target,
    the L1 (unweighted) → L2 → L3 (weighted, eq (21), 0.0012
    rearrangement) staged search over both L0 predictors; wraps an
    `LspReconstructor` so encoder/decoder MA states stay in lock-step.
  - **§3.3 perceptual weighting** (`perceptual_weighting`) — eq (28)
    LARs from k₁/k₂, eq (29) subframe interpolation, the eq (30)
    flat/tilted hysteresis, eq (31)/(32) γ₂ resonance adaptation, and
    the eq (33) weighted-speech recursion with cross-subframe memories.
  - **§3.4 open-loop pitch** (`open_loop_pitch`) — eq (34) frame
    correlation, three-section maxima (80–143 / 40–79 / 20–39), eq (35)
    normalisation, and the 0.85-weighted favour-lower-delays selection.
  - **§3.5/§3.6 impulse response + target** (`encode_target`) — `h(n)`
    of `W(z)/Â(z)` exactly as the spec words it; eq (36) residual; the
    three-stage target filter with §3.10-style difference-update.
  - **§3.7 closed-loop pitch** (`closed_loop_pitch`) — both spec
    window procedures, the eq (38) shift-and-add recursion, eq (37)
    normalised correlation, and the eq (39) 1/3-resolution fractional
    refinement through the compiled 13-tap b12 filter (the T1 ≥ 85
    integer-only rule).
  - **§3.8.1 fixed-codebook focused search** (`fixed_codebook_search`)
    — eq (43) bounded adaptive gain, eq (50) target update, eq (49)
    pre-filter fold, eq (51)/(52) φ/d, eqs (56)–(59) sign-folded
    search, and the eq (60) `thr₃` gate (K₃ = 0.4, 180-entry frame
    budget); recovers planted pulse sets exactly.
  - **§3.9.2 gain VQ search** (`gain_quantize`) — the six eq (63)
    inner products, the 2×2 normal-equations optimum, the 4-of-8 GA /
    8-of-16 GB conjugate preselection, and the exhaustive 32-pair
    search scored through the decode-side eq (73)/(74) reconstruction.
  - **the frame encoder** (`encode_chain::FrameEncoder`) — glues all
    stages in clause-3 order with the Figure-5 240-sample buffer
    timing, dual quantised/unquantised LSP interpolation tracks, the
    §3.7.2 odd-parity P0, the §3.9.3 index mapping, and the §3.10
    memory updates; every quantised quantity is produced by the
    decode-side modules so the encoder's local decoded state tracks a
    real decoder bit-for-bit.
  - **Table-8 packer + serial writer** (`parameters::pack_bit_array`,
    `serial::write_frame`, `FrameEncoder::encode_frame_to_serial`) —
    validated byte-exact against 8100+ corpus frames
    (parse→unpack→repack→rewrite reproduces every reference `.BIT`
    frame bit-for-bit).
  - **Encoder conformance harness** (`tests/encoder_conformance.rs`)
    over the ITU `.IN` corpus: 1:1 frame alignment proven, per-vector
    agreement floors pinned (L0 81–95%, L1 exact 33–80%, int(T1) ±2 on
    56–91% of active frames — floors at 75/25/45% with ratchet
    headroom), and the whole 3750-frame SPEECH vector encoded then
    decoded by our own decoder with finite non-silent output.
- Round 371 implements the **§4.2.1 long-term-postfilter 1/8-resolution
  fractional second pass** (eq (81)), completing the two-pass delay
  search. The integer anchor `T_0` (eq (80)) is refined to a fractional
  delay `T_0 + frac/8` (`frac ∈ {0 … 7}`) by maximising the
  pseudo-normalised correlation `R(T)²/E_T`: the length-33 short
  interpolation filter is tried first, then the chosen non-integer
  fraction is re-evaluated with the length-129 long filter and kept only
  if it raises `R′(T)` (the spec's longer-filter-replacement rule). The
  eq (78) harmonic filter now delays the reconstructed speech by the
  chosen fractional delay (long-filter interpolated). `LongTermDecision`
  gains a `frac` field. The fractional pass engages on ~26 % of enabled
  subframes across the SPEECH corpus vector (verified end-to-end through
  `decode_serial_frame_to_postfiltered`). A known-answer test confirms
  the polyphase kernels reproduce a band-limited sinusoid at the exact
  fractional delay, with the long filter more accurate than the short —
  validating the tap layout independently of any reference.
- Round 371 compiles the **§4.2.1 long-term-postfilter fractional-delay
  interpolation filters** (`tab_hup_s` / `tab_hup_l`) into the
  `tables` module. The staged CSVs (28 / 112 Q15 entries) — previously
  carried under `docs/audio/g729/tables/` but not wired in — are now
  emitted as `POSTFILTER_PITCH_INTERP_SHORT_Q15` (length-33 filter, 7
  phases × 4 taps) and `POSTFILTER_PITCH_INTERP_LONG_Q15` (length-129
  filter, 7 phases × 16 taps), with bounds-checked per-phase accessors
  `postfilter_interp_short` / `postfilter_interp_long`. The seven stored
  phases are the non-integer fractional offsets `frac = 1 … 7` of the
  1/8-resolution second pass (clause 4.2.1); phase `p` and phase `8 − p`
  are mirror images, and each phase has unity Q15 DC gain. This unblocks
  the long-term-postfilter fractional refinement (the README's standing
  "docs-gap" was stale — both tables were staged with full provenance).
- Round 357 wires the **§4.4 frame-erasure concealment** into the
  streaming decode chain (`decode_chain`). New `*_concealed` entry points
  (`decode_serial_frame_to_postfiltered_concealed`,
  `decode_serial_frame_to_speech_concealed`) reconstruct an erasure
  sentinel into concealed PCM instead of returning
  `FrameDecodeError::Erased`: the §4.4.1 LP-parameter repeat (last good
  frame), the §4.4.2 gain attenuation (eqs (93)/(94), compounding across
  an erasure run), the §4.4.3 gain-predictor-memory attenuation (eq (95)),
  and the §4.4.4 replacement excitation — periodic (adaptive-codebook
  repeat at the carried pitch delay, +1/subframe, bounded 143) or
  non-periodic (eq (96) random fixed codebook) by the inherited voicing
  class, which is latched on every good frame from the §4.2.1 long-term-
  postfilter decisions. A leading erasure before any good frame
  synthesizes silence. The previously-standalone `concealment` module's
  primitives are now driven end-to-end; the whole ERASURE corpus (300
  frames / variant, 60 erased) decodes to finite bounded PCM. The good-
  frame path is unchanged (the concealed entry points match the strict
  ones bit-for-bit on active frames).
- Round 357 adds a **whole-corpus PCM conformance harness**
  (`tests/pcm_conformance.rs`) over both the base-codec and Annex-A `.PST`
  references: first-subframe bit-accuracy across all clean vectors of both
  corpora, and a regression guard asserting the documented §3.9 fixed-
  point gain-saturation gap is a **bounded** multiplicative over-gain
  (whole-vector RMS ratio < 40×, measured 7–18.5×), not a divergence —
  the ceiling can be ratcheted toward 1.0 as the gap closes.
- Round 346 adds the **§B.4.4 CNG (comfort-noise) excitation synthesis**
  (`cng`): the silence-frame *excitation* path, fully spec-prose-sourced
  from the Annex B eqs (B.19)–(B.26). The eq (B.19) target-gain smoothing
  (`G̃_t` jump-to-SID after active, `7/8 : 1/8` relax across a non-active
  run, stateful `CngGainSmoother`); the eqs (B.20)–(B.22) per-subframe
  adaptive/fixed gain solving (`GainSolveTerms` reducing eq (B.20) to the
  eq (B.21) monic quadratic via the ACELP `Σe_f²=4` identity,
  lowest-abs-root `solve_fixed_gain`, eq (B.22) `Max{0.5, √(K/A)}`
  `ga_upper_bound`); the eqs (B.24)–(B.26) Gaussian mixture (`MixEnergies`,
  `solve_mix` for `(α=0.6, β>0)` with the `α=1/β=0` fallback,
  `build_cng_excitation`); and `CngRandom` (the §4.4.4 eq (96) recurrence,
  reset per active frame per §B.4.4) driving the pitch-lag / `Ga` /
  Gaussian draws. The `AnnexBStreamDecoder` now synthesizes the
  energy-controlled comfort-noise excitation for every non-active frame
  (`AnnexBOutput::ComfortNoise`, 80 samples + target gain) instead of the
  prior silence placeholder — 3 102 comfort-noise frames synthesized
  across the staged `g729b` corpus. The §B.4.2.2 SID-LSP VQ subset
  codebooks (needed for the comfort-noise *spectral envelope* / LP-filtered
  PCM) remain a documented docs-gap.
- Round 343 adds the **Annex B (DTX / CNG) decoder framing + routing**
  surface (`annex_b`): variable-length Annex B serial framing (the
  per-frame bit-count header selecting untransmitted / SID / active,
  plus both §B.4.5 erasure shapes — the `0x6B20` erased-sync word and an
  all-`0x0000` payload); SID bitstream unpack (Table B.2:
  predictor 1 / first-stage LSF 5 / second-stage LSF 4 / gain 5);
  the §B.4.2.1 5-bit non-uniform log energy dequantizer
  (`dequant_sid_energy_db`, fully prose-sourced); the §B.4.1/§B.4.5
  frame-type state machine `AnnexBDecoder` (erasure inheritance + SID
  persistence); and the end-to-end `AnnexBStreamDecoder` that decodes a
  whole Annex B `.bit` stream to per-frame PCM (active speech bit-exact
  through the §4.1 → §4.2 base chain, non-active frames a documented
  comfort-noise placeholder). Validated over all 10 staged `g729b`
  conformance sequences in `tests/annex_b_conformance.rs`. The §B.4.2.2
  SID-LSP VQ dequant + §B.4.4 CNG excitation *synthesis* are blocked on
  absent numeric tables (the SID-LSP VQ subset codebooks +
  `annexB-cng-lsp-sid-reset-Q15.csv`) and reported as a docs-gap.
- Round 340 adds the **§4.2 decoder post-processing cascade** — the
  five-stage adaptive postfilter (`postfilter_cascade`) that turns the
  §4.1.6 reconstructed speech `ŝ(n)` into final decoder output: §4.2.1
  long-term harmonic postfilter `H_p(z)` (anchored on the first
  subframe's transmitted integer pitch delay `int(T_1)`), §4.2.2
  short-term postfilter `H_f(z)`, §4.2.3 tilt compensation `H_t(z)`,
  §4.2.4 adaptive gain control (eq (88)–(90)), and §4.2.5 output
  high-pass `H_h2(z)` + ×2 up-scaling (eq (91)). Wired into
  `decode_chain::FrameDecoder` as the post-synthesis stage; covered by
  `tests/postfilter_conformance.rs`.
- Round 337 wires the **§3.10 / §4.1.6 LP synthesis** stage into the
  `decode_chain::FrameDecoder`, turning the §4.1 parameter decoder into a
  decoder that emits reconstructed-speech PCM. The `Synthesizer` (eq (40)
  past-excitation interpolator, eq (75) excitation build
  `u(n) = ĝ_p·v(n) + ĝ_c·c(n)`, eq (77) `1/Â(z)` filter) is now owned by
  the chain as the fifth piece of clause-4.3 cross-frame state (its
  eq (40)/(77) memories, zero-init):
  - New `FrameDecoder::decode_serial_frame_to_speech`,
    `decode_frame_kind_to_speech`, and `decode_parameters_to_speech`
    entry points — the three existing parameter entry points each gain a
    `*_to_speech` sibling that runs the §4.1 chain then the §3.10 →
    §4.1.6 synthesis and returns a `SynthesizedFrame` (80 reconstructed
    samples via `SynthesizedFrame::speech`). The synthesizer advances
    only on these calls; the parameter-only `decode_*` calls leave it at
    its start-up state. A codeword-domain error (or a §4.4 erasure
    sentinel) is surfaced before any synthesis runs, so the synthesizer
    memory is not advanced on a rejected frame.
  - New `FrameDecoder::synthesizer` accessor for inspection / tests.
  - This is the pre-postfilter reconstructed speech `ŝ(n)`: the §4.2
    post-processing cascade (postfilters, tilt, AGC, output high-pass +
    ×2) stays standalone and is wired separately — the §4.2.1 long-term
    postfilter's 1/8-resolution fractional pass remains the documented
    gap that blocks an end-to-end conformance match.
  - 5 new tests: the speech path equals the standalone
    decode-then-synthesize composition, produces finite samples and
    advances the eq (40) history off zero, carries synthesizer state
    across consecutive frames, agrees across the serial / frame-kind /
    parameter entry points, and rejects an erased frame without advancing
    the synthesizer. Crate unit tests 219 → 224.

- Round 331 wires the **§4.4 frame-erasure concealment** machinery — the
  decoder-side error recovery that reconstructs a frame whose 80-bit
  payload was lost. New `oxideav_g729::concealment` module:
  - `VoicingClass` + `Concealer::observe_good_frame` implement the §4.4
    voicing classifier: a 10 ms frame is **periodic** if at least one
    5 ms subframe had a long-term prediction gain > 3 dB (the eq (82)
    enable decision of `LongTermPostfilter`, reused via
    `LongTermDecision.gain > 0`), else **non-periodic**; an erased frame
    inherits the previous reconstructed frame's class.
  - `attenuate_fixed_gain` — eq (93) `ĝ_c^(m) = 0.98·ĝ_c^(m−1)`.
  - `attenuate_adaptive_gain` — eq (94) `ĝ_p^(m) = 0.9·ĝ_p^(m−1)`
    bounded `< 0.9`.
  - `attenuated_predictor_memory_db` — eq (95)
    `Û^(m) = 0.25·Σ_{i=1}^{4} Û^(m−i) − 4.0` bounded `≥ −14` (the §4.4.3
    "4 dB" attenuation of the `GainPredictor` memory).
  - `next_random_excitation` — the §4.4.4 non-periodic replacement
    excitation: eq (96) generator `seed = 31821·seed + 13849` (initial
    seed `21845`), fixed-codebook index from the 13 LSBs of one random
    number, sign field from the 4 LSBs of the next.
  - `periodic_pitch_delay` — the §4.4.4 periodic-case pitch-delay repeat
    (integer part of the previous frame's delay, +1 per subframe,
    bounded by 143).
  - eqs (93)–(96) transcribed clean-room from the EPUB raster figures
    `images/eq{93,94,95,96}.jpg`; 16 unit tests cover the classifier,
    inheritance, each attenuation, the LCG recurrence + bit extraction,
    and the bounded periodic-delay repeat.

- Round 326 wires the **§4.2.1 long-term (harmonic) postfilter
  `H_p(z)`** — the head of the §4.2 post-processing cascade and the
  last front stage that was still missing. New
  `oxideav_g729::long_term_postfilter` module:
  - `LongTermPostfilter` is the stateful integer-delay form of eqs
    (78)–(83). Per subframe it forms the eq (79) residual
    `r̂(n) = ŝ(n) + Σ_{i=1}^{10} γ_n^i·â_i·ŝ(n−i)` — clause 4.2.1's
    "numerator of the short-term postfilter `Â(z/γ_n)`" applied to the
    reconstructed speech, reusing `ShortTermPostfilter::weighted_num`
    (now `pub(crate)`) so the two stages share one `γ_n^i·â_i`
    definition.
  - The eq (80) integer delay search maximises
    `R(k) = Σ_{n=0}^{39} r̂(n)·r̂(n−k)` over the three candidates
    `[int(T_1) − 1, int(T_1) + 1]`; the eq (82) long-term-prediction-gain
    test disables the filter (`g_l = 0`) when
    `R′(T)² / Σ r̂(n)² < 0.5` (`ENABLE_THRESHOLD`); otherwise the eq (83)
    gain `g_l = Σ r̂(n)·r̂(n−T) / Σ r̂(n−T)²` is bounded to `[0, 1]`.
  - `filter_subframe` applies the eq (78) harmonic filter
    `H_p(z) = (1/(1+γ_p·g_l))·(1 + γ_p·g_l·z⁻ᵀ)` (`GAMMA_P` = 0.5),
    returning the postfiltered subframe plus a `LongTermDecision`
    (chosen delay + gain). The reconstructed-speech and residual
    histories are zero-initialised per clause 4.3 and carry across
    subframes (`MAX_DELAY` = 144).
  - 12 spec-cited unit tests: flat-filter residual identity, residual vs
    direct convolution, periodic-residual enable with near-unit gain,
    uncorrelated/alternating residual disable, eq (82) threshold gating,
    eq (78) algebra with history, history advance order, gain bound +
    delay window, determinism.
  - **Scope note / docs-gap:** the second-pass 1/8-resolution
    *fractional* delay refinement and its `tab_hup_s` (length 33) /
    `tab_hup_l` (length 129) interpolation filters are NOT wired — the
    spec prose names the filter lengths but defers their per-phase tap
    layout / convolution indexing to the electronic-attachment
    reference. The integer pass is the spec-prose-complete core.
- Round 319 wires the **§4.2.3 tilt compensation** and **§4.2.4
  adaptive gain control** stages of the §4.2 post-processing cascade —
  the two filters that sit between the round-313 short-term postfilter
  and the round-298 output high-pass:
  - New `oxideav_g729::tilt_compensation` module. `TiltCoefficients`
    derives the eq (87) tilt factor `k1' = −r_h(1)/r_h(0)` (the first
    reflection coefficient of the short-term postfilter impulse response
    `h_f(n)`, with `r_h(i) = Σ_{j=0}^{19−i} h_f(j)·h_f(j+i)`), selects
    `γ_t = 0.9` when `k1' < 0` else `0.2` (`GAMMA_T_NEG` / `GAMMA_T_POS`),
    and the gain term `g_t = 1 − |γ_t·k1'|`. Stateful
    `TiltCompensation` realises the eq (86) first-order FIR
    `H_t(z) = (1/g_t)·(1 + γ_t·k1'·z⁻¹)` per subframe
    (`sf(n) = (1/g_t)·(t(n) + γ_t·k1'·t(n−1))`), carrying the single
    input tap `t(n−1)` across subframes (zero-init per clause 4.3);
    `filter_subframe` derives the coefficients from the subframe `â_i`,
    `filter_subframe_with` takes pre-derived coefficients.
  - New `oxideav_g729::adaptive_gain_control` module. Stateful
    `AdaptiveGainControl` realises the eq (88) per-subframe energy ratio
    `G = Σ|ŝ(n)| / Σ|sf(n)|`, the eq (90) sample-by-sample gain
    smoothing `g(n) = 0.85·g(n−1) + 0.15·G` (`AGC_PREV_WEIGHT` /
    `AGC_TARGET_WEIGHT`), and the eq (89) scaling `sf′(n) = g(n)·sf(n)`,
    with the Table 9 start-up gain `g(−1) = 1.0` (`G_INIT`) carried
    across subframes (`g(−1) := g(39)`, clause 4.2.4). A silent
    postfiltered subframe (`Σ|sf| = 0`) holds the previous gain rather
    than dividing by zero.
  - New `impulse_response` helper on `ShortTermPostfilter` exposing the
    first 20 samples of the un-normalised `Â(z/γ_n)/Â(z/γ_d)` impulse
    response `h_f(n)` (the same response `gain_term` now sums for `g_f`,
    and the source of the §4.2.3 tilt factor's autocorrelation).
  - 23 new unit tests: tilt-factor hand computations (negative-`k1'` ⇒
    `γ_t = 0.9` and positive-`k1'` ⇒ `γ_t = 0.2`, both with hand-derived
    `g_t`), eq (86) FIR hand recursion, cross-subframe state continuity,
    flat-filter identity, determinism, finiteness; AGC eq (88) ratio,
    eq (90) smoothing hand recursion, fixed-point convergence to a
    constant `G`, silent-subframe gain-hold, cross-subframe carry,
    determinism; plus two `impulse_response`/`gain_term` consistency
    tests on `ShortTermPostfilter`.
  - The §4.2.1 long-term postfilter (its residual-domain two-pass pitch
    search, with `g_l` bounded by 1 and zeroed below a 3 dB long-term
    prediction gain, `γ_p = 0.5`) is now the **only** unwired §4.2 stage.

- Round 313 wires the **§4.2.2 short-term postfilter `H_f(z)`** — the
  spectral-shaping stage of the §4.2 cascade (between the §4.2.1
  long-term postfilter and the §4.2.3 tilt compensation) — in a new
  `oxideav_g729::short_term_postfilter` module:
  - `ShortTermPostfilter` (stateful) realises the eq (84) transfer
    function `H_f(z) = (1/g_f)·Â(z/γ_n)/Â(z/γ_d)` from the per-subframe
    quantised LP coefficients `â_i` and the two constant weight factors
    `GAMMA_N` (0.55) and `GAMMA_D` (0.70) (Table 5 / clause 4.2.2). Each
    subframe is filtered through the all-zero numerator `Â(z/γ_n)`
    (`r(n) = x(n) + Σ_{i=1}^{10} γ_n^i·â_i·x(n−i)` — the same residual
    clause 4.2.1 names) then the all-pole denominator `1/Â(z/γ_d)`
    (`y(n) = r(n) − Σ_{i=1}^{10} γ_d^i·â_i·y(n−i)`), and scaled by the
    eq (85) gain term `g_f = Σ_{n=0}^{19} |h_f(n)|` (the magnitude sum
    of the first 20 impulse-response samples of the un-normalised
    `Â(z/γ_n)/Â(z/γ_d)`), recomputed from the current `â_i` each call.
  - `filter_subframe` (one 40-sample subframe with its `â_i`) carries
    both 10-tap filter memories across subframes (zero-init per
    clause 4.3); `gain_term` (static, state-independent), `x_hist` /
    `y_hist` (inspection), and the `GF_IMPULSE_LEN` (20) constant round
    out the surface.
  - 9 unit tests (flat-filter identity + `g_f = 1`, clause-4.3 zero
    state, spec weight factors, `γ^i·â_i` weighting, `g_f` matched
    against an independent single-tap hand recursion, `g_f`
    state-independence, cross-subframe state continuity vs a continuous
    reference recursion, determinism, finiteness on a realistic
    coefficient set) plus a new
    `tests/short_term_postfilter_conformance.rs` harness asserting every
    output sample and every per-subframe `g_f` stays finite/positive
    across decode → synthesis → eq (84)/(85) over all active frames of
    the base-codec + Annex-A `.BIT` corpus (> 7 500 frames), with a
    determinism check on `ALGTHM.BIT`.
  - The §4.2.1 long-term postfilter (its residual-domain two-pass pitch
    search) and the §4.2.3/§4.2.4 tilt-compensation / adaptive-gain-
    control stages remain follow-up rounds; they slot around this module
    unchanged.

- Round 298 wires the **tail of the §4.2 post-processing cascade** —
  the §4.2.5 output high-pass filter + ×2 upscaling — in a new
  `oxideav_g729::post_process` module:
  - `OutputHighPass` (stateful) realises the eq (91) 2nd-order IIR
    `H_h2(z)` (100 Hz cut-off) from the already-compiled Q13 coefficient
    tables `HPF_PREPROC_100HZ_B_Q13` (`{7699, −15398, 7699}`) and
    `HPF_PREPROC_100HZ_A_Q13` (`{8192, 15836, −7667}`), with the
    clause-4.2.5 "×2 to restore the input signal level" upscaling folded
    into each output sample. The Q13 `a`-table stores the feedback gains
    sign-arranged for an additive recursion, so the difference equation
    is `y(n) = b0·x(n) + b1·x(n−1) + b2·x(n−2) + a1·y(n−1) + a2·y(n−2)`,
    reproducing the eq (91) denominator `1 − 1.9330735·z⁻¹ +
    0.93589199·z⁻²`. The four state taps (two input, two output) start
    zeroed per clause 4.3.
  - `filter_sample` (per-sample), `filter_in_place` (in-place slice),
    and `filter` (allocating) entry points all carry state across calls
    so a stream can be fed in arbitrary chunks; `b_coeffs` / `a_coeffs`
    expose the real-valued coefficients for inspection.
  - 10 unit tests (clause-4.3 zero state + table-sourced coefficients,
    eq (91) decimal match within one Q13 step, ×2 first-sample gain,
    exact DC zero `b0+b1+b2=0`, DC rejection, hand-worked impulse-
    response recursion pinning the additive-feedback sign, batch/per-
    sample API equivalence, cross-call state continuity, BIBO stability
    on bounded input, determinism) plus a new
    `tests/post_process_conformance.rs` harness asserting every
    output sample stays finite across decode → synthesis → eq (91) over
    all active frames of the base-codec + Annex-A `.BIT` corpus
    (> 7 500 frames), with a determinism check on `ALGTHM.BIT`.
  - The four front stages of the §4.2 cascade (long-/short-term
    postfilter, tilt compensation, adaptive gain control) remain
    follow-up rounds; they slot in front of this tail unchanged.
- Round 290 wires the §4.1.6 **LP synthesis** stage in a new
  `oxideav_g729::lp_synthesis` module — the first decoder stage to
  emit reconstructed-speech PCM:
  - `Synthesizer` (stateful) owns the two cross-subframe state pieces
    listed in clause 4.3 as zero-initialised: the eq (40)
    past-excitation buffer (`EXC_HISTORY` = 153 samples — the deepest
    eq (40) access `u(n − k − i)` with `n = 0`, `k = 143`, `i = 9` is
    index `−152`) and the eq (77) 10th-order `1/Â(z)` synthesis-filter
    memory.
  - `synthesize_frame(&DecodedFrame)` runs the spec §4.1.3 → §3.10 →
    §4.1.6 order per subframe: eq (40) interpolates the past
    excitation through the 31-tap `b_30`
    (`PITCH_INTERP_FILTER_SYNTHESIS_Q15`) at the decoded `(int_t,
    frac)` delay — mapped onto the `(k, t)` / `t ∈ {0,1,2}` fraction
    convention of eq (39)/(40) — to build the adaptive vector `v(n)`;
    eq (75) forms `u(n) = ĝ_p·v(n) + ĝ_c·c(n)` (with `c(n)` the
    post-eq (48) harmonic-enhanced codevector); eq (77) filters the
    excitation through `ŝ(n) = u(n) − Σ_{i=1}^{10} â_i·ŝ(n − i)`. The
    `v(n)`/`u(n)` computation interleaves per sample so short delays
    (`T < 40`) fold correctly onto the already-built current
    excitation; both state buffers advance after every subframe.
  - Typed outputs `SynthesizedFrame` (two subframes + `speech()` →
    the 80 time-ordered samples) and `SynthesizedSubframe` (the
    `adaptive` `v(n)`, `excitation` `u(n)`, and `speech` `ŝ(n)`
    40-sample vectors).
  - 10 unit tests (eq (40)/(75)/(77) algebra, the fraction-convention
    mapping, state advance, short-delay folding, determinism) plus a
    new `tests/lp_synthesis_conformance.rs` harness asserting every
    `v(n)`/`u(n)`/`ŝ(n)` stays finite across all active frames of the
    base-codec + Annex-A `.BIT` corpus (> 7 500 frames).
  - §4.2 post-processing (postfilter cascade, output high-pass, ×2
    upscaling) remains a follow-up round; this stage emits the raw
    §4.1.6 `ŝ(n)` that feeds it.

- Round 282 wires the §4.1 **per-frame decode parameter chain** in a
  new `oxideav_g729::decode_chain` module — the glue that turns one
  ITU serial frame into fully-typed parameter structs through every
  decode piece already landed:
  - `FrameDecoder` (stateful) owns the §3.2.4 LSP MA history, the
    §3.2.5 interpolation memory, the §3.9.1 4-tap gain-predictor
    history, the eq (47) `β = ĝ_p^(m−1)` source (clause 4.3 /
    Table 9 init `0.8`, exposed as `BETA_INIT`), and the previous
    frame's `int(T2)` for the §4.1.2 parity concealment (clause-4.3
    zero default).
  - Three entry points: `decode_serial_frame` (one 164-byte ITU
    serial frame), `decode_frame_kind(&FrameKind)`, and
    `decode_parameters(&Parameters)`, each running the clause-4.1
    order: §4.1.1 `L0..L3` → `ω̂` → per-subframe interpolated LSPs →
    LP coefficients; §4.1.2 parity recompute with the spec
    substitution on mismatch ("the delay value T1 is set to the
    integer part of the delay value T2 of the previous frame", T2
    then re-derived per §4.1.3 from that T1); §4.1.3 `(P1, P2)` →
    `(T1, T2, t_min)`; §4.1.4 `(C, S)` → pulses → eq (45)
    codevector → eq (48) sharpening with β threaded from the
    previous subframe's `ĝ_p`; §4.1.5 `(GA, GB)` → §3.9.3 demap →
    eqs (73)/(74) `(ĝ_p, γ̂)` → `ĝ_c = γ̂·g′_c` with the eq (66)
    energy computed over the harmonic-enhanced `c(n)` (§3.10) and
    the eq (72) history advance.
  - Typed outputs `DecodedFrame` (raw codewords, parity outcome,
    `ω̂`, `t_min`, two subframes) and `SubframeDecode` (interpolated
    LSPs, LP coefficients, pitch delay, pulses, β, post-eq (48)
    codevector, `(ĝ_p, γ̂)`, §3.9.1 prediction path, `ĝ_c`); typed
    `FrameDecodeError` (Serial / Erased / Lsp / FixedCodebook /
    Gain). §4.1.6 synthesis, §4.2 post-processing, and §4.4 erasure
    concealment remain unwired — an erasure sentinel returns
    `FrameDecodeError::Erased` without advancing state.
  - 8 unit tests (Table-9 start-up state, hand-sequenced
    piece-equivalence across both subframes, §4.1.2 substitution
    incl. first-frame zero-default, erasure rejection without state
    advance, typed out-of-domain errors on hand-built codewords,
    serial↔parameter entry-point equality, eq (47) β clamp under GA
    sweeps) plus a new `tests/decode_chain_conformance.rs` harness
    that decodes 18 222 active frames across all 19 staged
    base-codec + Annex-A `.BIT` vectors, pinning the §3.2.4
    stability-clamp envelope, §4.1.3 delay windows, Table-7 track
    residues, gain/LP finiteness, exactly 60/300 §4.1.2 concealment
    activations on `PARITY.BIT` (each substituting exactly the
    previous frame's `int(T2)`), erasure-sentinel rejection inside
    `ERASURE.BIT` / `OVERFLOW.BIT`, and frame-by-frame determinism
    of two independent chains.

- Round 274 lands the §3.8 / §4.1.4 **pitch sharpening** step that
  the round-266 fixed-codebook decode had deferred, in a new
  `oxideav_g729::pitch_sharpen` module:
  - `clamp_beta(g_p_prev: f32) -> f32` applies spec eq (47) —
    `β = ĝ_p^(m−1)` bounded by `0.2 ≤ β ≤ 0.8`, i.e. the previous
    subframe's quantised adaptive-codebook gain clamped to the
    closed `[0.2, 0.8]` pitch-gain interval (`NaN` → `0.2` floor so
    the recurrence stays finite).
  - `sharpen(c: &[i8; 40], int_t: i32, g_p_prev: f32) -> [f32; 40]`
    applies spec eq (48) — `c(n) += β·c(n − T)` for `n = T … 39`,
    realising the §3.8 adaptive pre-filter `P(z) = 1/(1 − β·z^−T)`.
    The forward in-place sweep reads the already-modified
    `c(n − T)`, so the recursive pole is realised exactly and the
    geometric pitch echo (`β`, `β²`, … at `m + T`, `m + 2T`, …)
    propagates correctly. The modification is applied **only** when
    the current subframe's integer pitch delay `int(T) < 40`
    (clause 4.1.4 "less than the subframe size 40"); for
    `int(T) ≥ 40` (or a non-positive delay) the codevector is
    returned unchanged, promoted to `f32`. The eq (47) clamp is
    applied internally so callers pass the raw previous-subframe
    `ĝ_p` directly.
  - `codevector_energy(&[f32; 40]) -> f32` reads the
    post-sharpening `Σ c(n)²` (the unmodified four-pulse codevector
    has energy exactly `4.0`); exposed for the §3.9.1 eq (66)
    energy term and the unit tests.
  - Tests pin the eq (47) clamp at its boundaries, the eq (48)
    geometric-train recurrence (single seed pulse → taps at
    `0, T, 2T, 3T` with weights `1, β, β², β³`), the recursive
    interaction of overlapping pulses (reading the modified
    `c(n − T)`), the `int(T) ≥ 40` / non-positive no-op guards, and
    the head-preservation + tail-echo structure on a real decoded
    codevector.
- Round 266 wires the §3.8 / §4.1.4 fixed (algebraic) codebook
  decode that maps the transmitted `(C, S)` codewords into the
  per-subframe pulse layout and the 40-sample codevector `c(n)`,
  in a new `oxideav_g729::fixed_codebook` module:
  - `decode_positions(c: u16) -> Result<([u8; 4], u8), _>` inverts
    spec eq (62) — the 13-bit `C` codeword splits as 3+3+3+4 bits
    across the four Table-7 tracks; the 4-bit track-3 field
    carries `2·(m_3/5) + jx` so `jx` is the LSB and the upper 3
    bits are `m_3/5`. Returns the `[m_0, m_1, m_2, m_3]` positions
    in spec order alongside `jx ∈ {0, 1}`. Returns
    `FixedCodebookError::CTooWide` on a 14+-bit input.
  - `decode_signs(s: u8) -> Result<[i8; 4], _>` inverts spec
    eq (61) — `S = s_0 + 2·s_1 + 4·s_2 + 8·s_3`; bit `k` carries
    `s_k = 1` for a positive sign (clause 3.8.2 prose). Returns
    `±1` per pulse to match spec eq (45). Returns
    `FixedCodebookError::STooWide` on a 5+-bit input.
  - `decode_pulses(c, s) -> Result<FixedCodebookPulses, _>` —
    per-subframe wrapper that ties positions + signs into one
    `Copy` struct.
  - `build_codevector(&FixedCodebookPulses) -> [i8; 40]`
    constructs spec eq (45) — four signed unit impulses at the
    decoded pulse positions, zero elsewhere. The energy is
    exactly `4` (one ±1 contribution per pulse).
  - `decode_frame(&Parameters) -> Result<FrameFixedCodebook, _>`
    — per-frame entry point. Threads `(C1, S1)` into subframe 1
    and `(C2, S2)` into subframe 2 per spec §4.1.5.
  - `encode_positions(&[u8; 4]) -> Option<u16>` and
    `encode_signs(&[i8; 4]) -> Option<u8>` — encoder-side
    forward mappings (spec eqs (62) / (61)) exposed publicly
    for round-trip property tests and encoder wire-up.
- 12 new unit tests in `src/fixed_codebook.rs` pin the
  algorithmic invariants: spec eq (61) worked examples (`S =
  0b0000` / `0b1111` / `0b0001` / `0b1000` / `0b0101`), spec
  eq (62) worked examples (`C = 0` / `1` / `512` (`jx = 1,
  m_3 = 4`) / `0x1FFF` (full saturation, `m = [35, 36, 37, 39]`)),
  out-of-domain rejection on both decode primitives, full-domain
  Table-7 track-residue envelope on `decode_positions` (every
  `C ∈ 0..8192`), full-domain encode↔decode round-trip on `S`
  AND `C` (pinning both decode and encode recipes
  simultaneously), `encode_positions` off-track rejection,
  `encode_signs` non-±1 rejection, `decode_pulses` ordering,
  `build_codevector` eq (45) construction, codevector energy
  invariant `= 4` across a representative sweep, distinct-
  positions invariant across the full `C` domain, and
  `decode_frame` per-subframe thread-through.
- 2 new integration tests in `tests/serial_conformance.rs`:
  `fixed_codebook_in_domain_on_full_corpus` walks every `.BIT`
  file in `g729-core/` + `g729a/` and pins per-frame /
  per-subframe invariants (Table-7 track residue, ±1 signs,
  distinct positions, energy = 4, `jx ∈ {0, 1}`);
  `fixed_codebook_round_trips_fixed_corpus` walks the staged
  `FIXED.BIT` corpus (the ITU `READMETV.txt`-documented
  fixed-codebook exerciser) and pins
  `encode_positions(decoded) == params.c1 / c2` and
  `encode_signs(decoded) == params.s1 / s2` exactly on every
  active frame.
- `oxideav_g729::fixed_codebook` is exposed at the crate root
  via `pub mod fixed_codebook;`; the crate-level rustdoc is
  extended with the round-266 surface.
- Spec gap (unchanged from previous rounds): clause 3.8
  eq (48) / (48a) pitch sharpening of `c(n)` when `int(T) < 40`
  needs the round-255 pitch-delay output AND the previous
  subframe's quantised adaptive-codebook gain `β`. Per clause
  4.1.4: "If the integer part of the pitch delay *T* is less
  than the subframe size 40, *c*(*n*) is modified according to
  equation (48)." This step is left for a follow-up round; it
  modifies the codevector in place after `build_codevector`.

- Round 255 wires the §4.1.3 pitch-delay decode that maps the
  transmitted `(P1, P2)` indices into per-subframe fractional pitch
  delays `(T1, T2)`, in a new `oxideav_g729::pitch_decode` module:
  - `decode_t1_from_p1(p1: u8) -> PitchDelay` evaluates spec image
    `f0027-01.jpg` (clause 4.1.3): `if P1 < 197 then int(T1) =
    (P1 + 2) / 3 + 19`, `frac = P1 − 3·int(T1) + 58`; `else
    int(T1) = P1 − 112`, `frac = 0`. Total spec coverage of the
    8-bit `P1` field is the union of a fractional branch (1/3
    resolution over `T1 ∈ [19⅓, 85]`) and an integer-only branch
    (over `T1 ∈ [86, 143]`).
  - `derive_t_min(int_t1: i32) -> i32` evaluates spec image
    `f0027-02.jpg`: `t_min = int(T1) − 5`, floor at `20`, ceiling
    at `t_max = 143` then `t_min = t_max − 9`. The output is in
    `[20, 134]` across the full `int(T1) ∈ [19, 143]` decode
    range.
  - `decode_t2_from_p2(p2: u8, t_min: i32) -> PitchDelay`
    evaluates spec image `f0027-03.jpg`: `int(T2) = (P2 + 2) / 3
    − 1 + t_min`, `frac = P2 − 2 − 3·((P2 + 2) / 3 − 1)`. The
    5-bit `P2` field covers the 1/3-resolution `T2 ∈
    [t_min − 1/3, t_max + 1/3]` window exactly.
  - `decode_frame(&Parameters) -> FramePitchDelays` — per-frame
    wrapper that chains `decode_t1_from_p1` → `derive_t_min` →
    `decode_t2_from_p2` in the spec §4.1.3 order. The returned
    `FramePitchDelays` is a `Copy` struct carrying `t1`, `t2`,
    and the spec `t_min` (preserved for callers driving the
    §4.1.2 parity-concealment path).
  - `encode_p1(delay) -> Option<u8>` and `encode_p2(delay, t_min)
    -> Option<u8>` — the spec §3.7 encode-side forward mappings
    (eqs (41) / (42), spec images `eq41.jpg` / `eq42.jpg`)
    exposed publicly so callers (encoders, fixture builders) can
    round-trip the pair without re-deriving the algebra. They
    return `None` for out-of-domain `(int_t, frac)` pairs.
  - `PitchDelay` — `Copy` struct carrying `int_t` (the integer
    part of the fractional pitch delay) and `frac ∈ {-1, 0, 1}`
    (the spec §3.7 1/3-resolution fractional component). The
    full-precision pitch delay is `int_t + frac/3`.
  - `FramePitchDelays`, `T_MIN_FLOOR = 20`, `T_MAX_CEIL = 143`,
    `T_WINDOW = 9`, `P1_DOMAIN = 256`, `P2_DOMAIN = 32`,
    `P1_FRACTIONAL_LIMIT = 197` — public surface constants.
- 13 new unit tests in `src/pitch_decode.rs`: spec-image worked
  examples on `decode_t1_from_p1` boundaries (`P1 ∈ {0, 1, 196,
  197, 255}`), full-domain envelope check on `decode_t1_from_p1`
  (every `P1 ∈ 0..256` lands in `int(T1) ∈ [19, 143]`, `frac ∈
  {-1, 0, 1}`); spec-image worked examples on `derive_t_min`
  (mid-range, floor edge, ceiling edge), full-domain envelope
  check on `derive_t_min` (the entire 9-step subframe-2 search
  window fits inside `[20, 143]` for every `int(T1) ∈ [19, 143]`);
  spec-image worked examples on `decode_t2_from_p2` (`P2 ∈ {0, 2,
  31}` at `t_min = 50`), full-domain envelope check (`P2 ∈ 0..32`
  × `t_min ∈ [20, 134]` always lands in `int(T2) ∈ [t_min − 1,
  t_min + 10]`, `frac ∈ {-1, 0, 1}`); encode↔decode round-trip
  over the **full `P1 ∈ 0..256` domain** (pins both the eq (78a)
  decode and the eq (41) encode simultaneously); same over the
  full `P2 ∈ 0..32` × `t_min ∈ [20, 134]` domain; encode-side
  out-of-domain rejection on both `P1` and `P2`;
  `decode_frame` threads the right field into the right subframe
  (P1↔P2 swap detection); constants match the documented spec
  values.
- 2 new integration tests in `tests/serial_conformance.rs` against
  the staged conformance corpus:
  - `pitch_decode_in_domain_on_full_corpus` walks every `.BIT`
    file in `g729-core/` + `g729a/`, unpacks each active frame's
    `Parameters`, runs `decode_frame`, and pins that every decoded
    `(T1, T2, t_min)` lies in the spec-stated domain: `int(T1) ∈
    [19, 143]`, `t_min ∈ [20, 134]`, `int(T2) ∈ [t_min − 1,
    t_min + 10]`, and both `frac` components in `{-1, 0, 1}`.
  - `pitch_decode_round_trips_pitch_corpus` walks the staged
    `PITCH.BIT` sequence (the ITU `READMETV.txt` self-documents
    this as the pitch-delay exerciser) for both the `g729-core/`
    base codec and the `g729a/` Annex-A corpus, runs `decode_frame`
    on every active frame, and pins that `encode_p1(t1) == params.p1`
    and `encode_p2(t2, t_min) == params.p2` exactly. This is the
    strongest in-corpus guarantee available without a known-good
    reference output: the only way a frame can fail is for the
    eq (78a) / eq (79) / eq (80) decode recipe to disagree with
    the eq (41) / eq (42) encode recipe on an ITU-encoded frame.

- Round 249 wires the §3.9.3 gain-quantiser codeword-mapping layer
  between the round-225 transmitted-index unpacker and the round-231
  conjugate-structure codebook lookup, in a new
  `oxideav_g729::gain_index_map` module:
  - `demap_ga(transmitted: usize) -> Result<usize, _>` and
    `demap_gb(transmitted: usize) -> Result<usize, _>` — decoder-side
    `imap1` / `imap2` inverse-permutation primitives that map the
    on-wire GA / GB indices back into the codebook-index domain
    (`0..NCODE1` and `0..NCODE2`).
  - `map_ga(codebook: usize) -> Result<usize, _>` and
    `map_gb(codebook: usize) -> Result<usize, _>` — encoder-side
    `map1` / `map2` forward permutations, symmetric to the demap
    primitives.
  - `demap_frame(&Parameters) -> Result<DemappedGainIndices, _>` —
    per-frame wrapper that demaps all four transmitted gain indices
    in one call. `DemappedGainIndices` is a `Copy` struct carrying
    `(ga1, gb1, ga2, gb2)` in the codebook-index domain.
  - `GainIndexMapError` — typed surface for out-of-range inputs,
    with `GaOutOfRange { index }` / `GbOutOfRange { index }`
    variants. Implements `Display` + `std::error::Error`.
- The existing `gain_reconstruct::reconstruct_frame_gains` entry
  point now applies the §3.9.3 inverse permutation internally before
  the §3.9.2 codebook lookup, so the `(GA, GB) → (ĝ_p, γ̂)` pipeline
  is spec-conformant end-to-end from the on-wire bits. A new
  `gain_reconstruct::reconstruct_gains_from_transmitted(t_ga, t_gb)`
  helper exposes the same demap-then-reconstruct path as a per-pair
  primitive for callers working off bare integers.
- `gain_reconstruct::GainReconstructError` gains an
  `IndexMap(GainIndexMapError)` variant + `From<GainIndexMapError>`
  conversion, so demap failures from the per-frame / per-pair
  wrappers surface through the existing error type without bypassing
  it.
- 17 new unit tests in the `gain_index_map` module: forward / inverse
  round-trip on both stages and both compositions
  (`map ∘ demap == id`, `demap ∘ map == id`); codebook-domain
  containment of the demap output; bijection check on each demap;
  out-of-range boundary on every entry point; non-identity assertion
  on both demaps (locks against a CSV regression that emits the
  identity table); per-frame wrapper threads the right indices into
  the right stages; per-stage zero-index pin against the staged
  `imap1[0]` / `imap2[0]` literals.
- 3 new unit tests in `gain_reconstruct`: round-249 demap-before-
  lookup property (the frame wrapper's output matches
  `reconstruct_gains(demap_ga(t_ga), demap_gb(t_gb))` and explicitly
  does not match the bare `reconstruct_gains(t_ga, t_gb)` when the
  permutation is non-trivial); `reconstruct_gains_from_transmitted`
  agrees with the hand-composed pipeline over the full transmitted
  domain; out-of-range inputs surface through the `IndexMap` error
  variant. The pre-existing `frame_wrapper_threads_per_subframe_indices`
  test is updated to compose `demap_ga` / `demap_gb` in its expected-
  value computation so it stays semantically correct.

- Round 239 wires the §3.9.1 / §4.1.5 4th-order MA gain prediction
  stage on top of the round-231 conjugate-structure gain-VQ output,
  in a new `oxideav_g729::gain_predict` module:
  - `GainPredictor::new()` owns the four-slot prediction-error
    history `[Û^(m-1), Û^(m-2), Û^(m-3), Û^(m-4)]` initialised per
    spec Table 9 / §4.3 to `[-14, -14, -14, -14]` dB. Slot 0 is the
    most-recent (eq (69) `b_1`-weighted) slot.
  - `GainPredictor::codevector_energy_db(c: &[f32; 40]) -> f32` —
    spec eq (66) `E = 10·log10((1/40)·Σ_{n=0..39} c(n)^2)`, with a
    `1e-30` `log10` floor so the all-zero corner stays finite.
  - `GainPredictor::predict_only(c)` /
    `predict_only_from_energy(e_db)` evaluate eq (69) + eq (71)
    without advancing the history; return `PredictedGain { e_db,
    e_tilde_db, g_c_prime }`.
  - `GainPredictor::predict_and_update(c, gamma_hat)` /
    `predict_and_update_from_energy(e_db, gamma_hat)` evaluate the
    full per-subframe path and advance the history per eq (72)
    decode form `Û^(m) = 20·log10(γ̂)`; return `(ĝ_c =
    γ̂ · g'_c, PredictedGain)`.
  - `GainPredictor::push_quantised_error(gamma_hat)` — low-level
    eq (72) history advance, exposed for callers driving custom
    loops (concealment paths).
  - `FIXED_CODEBOOK_MEAN_ENERGY_DB = 30.0` (spec §3.9.1 `Ē`),
    `GAIN_PREDICTOR_INIT_DB = -14.0` (spec Table 9), and
    `CODEVECTOR_LEN = 40` (spec §3.8 / eq (66) averaging length)
    as crate-public constants.
- 13 new unit tests for the gain-predict path: Table 9 init shape;
  first-subframe `Ẽ^(0)` matches `-14 · Σ b_i` from the staged Q13
  coefficients; eq (66) at `-10 dB` for the 4-pulse codevector;
  eq (66) finite on all-zero input; eq (71) `g'_c = 1` when
  `Ẽ = E - Ē`; eq (72) decode form decibel scaling
  (`γ̂ = 1.0 → 0 dB`, `γ̂ = 10 → 20 dB`, `γ̂ = 0.1 → -20 dB`);
  non-positive `γ̂` floors to finite very-negative `Û^(m)`;
  predict-and-update consistency vs side-by-side `predict_only`;
  `predict_only` is side-effect-free; steady-state convergence
  after `MA_NP` pushes with constant `γ̂`; history index 0 is the
  most-recent (`b_1`-weighted) slot; eq (65) `ĝ_c = γ̂ · g'_c`
  sweep across six `γ̂` values.
- 1 new corpus integration test `gain_predict_finite_on_full_corpus`
  that walks every `.BIT` file in
  `docs/audio/g729/conformance/{g729-core,g729a}/`, runs
  `reconstruct_frame_gains` per active frame, then drives a fresh
  `GainPredictor` through the resulting `γ̂` sequence with a
  representative 4-pulse codevector — asserts finite
  `(g'_c, ĝ_c)` end-to-end, `ĝ_c ∈ [0, 1e6]` (defensive envelope
  given the synthetic codevector understates real `E`), and finite
  history slots after the eq (72) update.
- Round 231 wires the §3.9.2 / §4.1.5 conjugate-structure gain-VQ
  decode-side reconstruction on top of the round-225 §4.1 parameter
  unpacker, in a new `oxideav_g729::gain_reconstruct` module:
  - `reconstruct_gains(ga: usize, gb: usize) -> Result<QuantisedGains,
    GainReconstructError>` evaluates spec eqs (73) / (74):
    `ĝ_p = GA[GA][0] + GB[GB][0]` (column 0, Q14 in both stages) and
    `γ̂ = GA[GA][1] + GB[GB][1]` (column 1, Q12 in both stages); the
    summation runs in `i32` per Q-format and converts to `f32` at the
    boundary. Out-of-range indices surface as
    `GaOutOfRange { index }` / `GbOutOfRange { index }` rather than
    panicking.
  - `reconstruct_frame_gains(&Parameters) -> Result<[QuantisedGains;
    2], GainReconstructError>` per-frame wrapper that threads
    `(GA1, GB1)` into subframe 1 and `(GA2, GB2)` into subframe 2,
    matching the §4.1.5 ordering.
  - `QuantisedGains` — `Copy` struct carrying `g_p_hat`
    (quantised adaptive-codebook gain `ĝ_p`) and `gamma_hat`
    (quantised fixed-codebook gain correction factor `γ̂`). The
    actual quantised fixed-codebook gain `ĝ_c = γ̂ · g'_c` is left
    for a follow-up round: `g'_c` is produced by the §3.9.1 4th-order
    MA prediction stage which is stateful and not yet wired.
- Round 231 wires the §3.9.2 / §3.9.3 numeric tables that feed the
  reconstruction:
  - `tables::GAIN_QUANT_CODEBOOK_GA_Q14_Q12` — first-stage codebook
    `gbk1`, shape `[[i16; 2]; NCODE1]` = `[[i16; 2]; 8]`. Column 0
    is the adaptive-codebook-gain contribution in Q14, column 1 is
    the fixed-codebook-gain correction contribution in Q12.
  - `tables::GAIN_QUANT_CODEBOOK_GB_Q14_Q12` — second-stage codebook
    `gbk2`, shape `[[i16; 2]; NCODE2]` = `[[i16; 2]; 16]`, same
    per-column Q-formats.
  - `tables::GAIN_QUANT_GA_PERMUTATION` / `GA_INVERSE_PERMUTATION`
    (8 entries each) and `GAIN_QUANT_GB_PERMUTATION` /
    `GB_INVERSE_PERMUTATION` (16 entries each) — the spec §3.9.3
    transmission-side robustness mapping that reorders GA / GB
    indices before transmission so a single-bit channel error lands
    on a perceptually-close codebook entry.
  - `tables::GAIN_QUANT_GA_THRESHOLDS_Q14` (4 entries) and
    `GAIN_QUANT_GB_THRESHOLDS_Q15` (8 entries) — encoder-side
    partial-search thresholds (`NCODE1 - NCAN1 = 4`,
    `NCODE2 - NCAN2 = 8`); staged here for completeness alongside
    the codebooks but not consumed by the round-231 decode-side
    reconstruction.
  - New spec-dimension constants: `NCODE1 = 8`, `NCODE2 = 16`,
    `GAIN_VQ_DIM = 2`, `GAIN_VQ_COL_GP = 0`, `GAIN_VQ_COL_GC = 1`.
  - Bounds-checked lookup helpers `tables::gain_ga_entry(ga)` /
    `gain_gb_entry(gb)` returning a borrowed 2-element row.
- 10 new unit tests in `src/gain_reconstruct.rs` pin the algorithmic
  invariants: per-row CSV-literal match at (0, 0) and (`NCODE1 - 1`,
  `NCODE2 - 1`); out-of-range GA / GB rejection (each variant);
  out-of-range-first-wins error precedence; every (GA, GB) pair in
  the 8 × 16 domain yields finite gains; `ĝ_p` lies in `[0, 2]` and
  `γ̂` lies in `[0, 11]` for every pair (Q-format-divisor isolation
  check; worst-case row pairs reach ≈10.12); per-column delta
  isolation (varying GA at fixed GB moves `ĝ_p` by the Q14-scaled
  column-0 delta and `γ̂` by the Q12-scaled column-1 delta);
  hand-picked pair `(5, 11)` matches the algebra; `reconstruct_frame_gains`
  threads `(GA1, GB1)` and `(GA2, GB2)` into the right subframe;
  codebook width matches the published `GAIN_VQ_DIM` constant.
- 7 new structural tests in `tests/tables_shape.rs` for the staged
  tables: shape (`NCODE1 × GAIN_VQ_DIM` for GA, `NCODE2 × GAIN_VQ_DIM`
  for GB), with the row counts cross-checked against
  `1 << GA_BITS == NCODE1` and `1 << GB_BITS == NCODE2`; first-row
  CSV-literal pins (GA[0] = `[1, 1516]`, GB[0] = `[826, 2005]`);
  column-constant convention pinned to the (0 = `g_p`, 1 = `γ`)
  layout; both permutations are complete covers of `0..NCODE`;
  inverse-permutation property (`imap ∘ map == id`); threshold tables
  are strictly ascending (the §3.9.2 partial-search band-ordering
  invariant); threshold lengths match `NCODE - NCAN`; helper-vs-
  constant equivalence.
- 1 new integration test against the staged conformance corpus:
  `gain_reconstruct_in_domain_on_full_corpus` walks every `.BIT`
  file in `g729-core/` + `g729a/`, unpacks every active frame's
  parameters, runs `reconstruct_frame_gains`, and pins that every
  reconstructed pair is finite, every `ĝ_p` lies in `[0, 2]`, and
  every `γ̂` lies in `[0, 11]`. With the round-225 corpus walker
  this confirms that no transmitted (GA, GB) index pair in the ITU
  conformance corpus ever drives the reconstruction off the
  plausibility envelope.

- Round 225 wires the §4.1 / Table-8 parameter unpacker, splitting
  the round-191 serial 80-bit payload into the 15 typed codeword
  indices the §4.1 decode procedure consumes, in a new
  `oxideav_g729::parameters` module:
  - `Parameters` — `Copy` struct carrying the per-frame indices
    `l0` / `l1` / `l2` / `l3` (§3.2.4 LSP VQ), `p1` / `p0` (§3.7
    / §3.7.2 subframe-1 pitch + parity), `c1` / `s1` (§3.8
    subframe-1 fixed codebook), `ga1` / `gb1` (§3.9.2 subframe-1
    conjugate-structure gain VQ), and the matching `2`-suffixed
    set for subframe 2. Field-width split: 1+7+5+5 + 8+1+13+4+3+4
    + 5+13+4+3+4 = 18 + 33 + 29 = **80** bits.
  - `unpack_parameters(&FrameKind) -> Result<Parameters,
    ParameterError>` — frame-level entry point; rejects
    `FrameKind::Erased` with `ParameterError::Erased` (the §4.4
    concealment path applies for an erasure-sentinel frame and
    consumes no transmitted bits).
  - `unpack_bit_array(&[bool; 80]) -> Parameters` — lower-level
    variant taking the 80-bit array directly, useful for
    unit-testing without spinning the framing layer.
  - `Parameters::pitch_parity_ok(&self) -> bool` — §3.7.2 /
    §4.1.2 parity check; the parity-init value is **pinned to 1
    (odd-parity convention)** based on the staged corpus
    (every active frame of `SPEECH.BIT` / `g729a/SPEECH.BIT`
    has `P0 = 1 XOR XOR_reduce(six_MSBs(P1))`).
  - Per-codeword bit-width constants at the module surface:
    `P1_BITS = 8`, `P0_BITS = 1`, `C_BITS = 13`, `S_BITS = 4`,
    `GA_BITS = 3`, `GB_BITS = 4`, `P2_BITS = 5`; aggregate
    constants `FIXED_CODEBOOK_BITS_PER_FRAME = 34`,
    `GAIN_QUANT_BITS_PER_FRAME = 14`,
    `PITCH_BITS_PER_FRAME = 14` express the frame-level
    grouping (`18 + 34 + 14 + 14 = 80` matches the round-189
    `LSP_TOTAL_BITS = 18` + the round-225 totals against
    `BITS_PER_FRAME = 80`).
  - Per-codeword start offsets `(0, 1, 8, 13, 18, 26, 27, 40,
    44, 47, 51, 56, 69, 73, 76)` are computed at compile time
    from the width array and statically asserted to sum to
    `BITS_PER_FRAME` so the layout can never silently drift.
  - 9 new unit tests pin the algorithmic invariants: all-zero /
    all-ones boundary; L1/L2/L3 saturation lands in
    `NC0`/`NC1`; **per-bit flip** test (each of the 80 array
    slots changes exactly the codeword whose
    `[start, start + width)` window contains it — locks
    MSB-first / Table-8-top-to-bottom against any off-by-one);
    documented start offsets hold; round-trip pack-then-unpack
    on a hand-chosen `Parameters` recovers every field
    bit-exactly; §3.7.2 parity-rule worked checks on three
    crafted P1 values; erasure rejection; high-level vs
    low-level entry-point agreement.
  - 2 new integration tests against the staged conformance
    corpus walk every `.BIT` file in `g729-core/` + `g729a/`:
    every active frame's `Parameters` has every field in its
    spec-stated domain (L1 < NC0, L2/L3 < NC1, C1/C2 < 2^13,
    signs < 2^4, GA < 2^3, GB < 2^4, P2 < 2^5); `SPEECH.BIT`
    produces **zero** parity mismatches (clean encoder output);
    `PARITY.BIT` produces a **non-zero** number of mismatches
    (the dedicated §4.1.2 concealment-path exerciser, both
    for `g729-core` and `g729a` corpora).

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
