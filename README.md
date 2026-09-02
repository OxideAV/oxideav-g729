# oxideav-g729

[![CI](https://github.com/OxideAV/oxideav-g729/actions/workflows/ci.yml/badge.svg)](https://github.com/OxideAV/oxideav-g729/actions/workflows/ci.yml) [![crates.io](https://img.shields.io/crates/v/oxideav-g729.svg)](https://crates.io/crates/oxideav-g729) [![docs.rs](https://docs.rs/oxideav-g729/badge.svg)](https://docs.rs/oxideav-g729) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Pure-Rust ITU-T **G.729** (CS-ACELP, 8 kbit/s) narrowband speech codec.
Zero C dependencies, no FFI, no `*-sys` crates.

## Status

Clean-room rebuild in progress, grown one spec-cited unit at a time
from the published ITU-T G.729 Recommendation prose and the ITU
electronic-attachment numeric data tables. The crate now carries
**both directions**: the decode path end-to-end to post-processed PCM,
and (round 382) the **entire clause-3 encoder analysis chain** with a
working `.IN` → `.BIT` path. As of round 388 the crate is **registered
into the `oxideav-core` codec registry** (id `"g729"`, decode + encode
over the raw 10-octet wire framing) with the dual-API
`decoder::make_decoder` / `encoder::make_encoder` endpoints. Round 419
put the whole decode path on the clause-5 fixed-point grid; round 438
did the same for the encoder's §3.1–§3.2.3 front end; round 452
root-caused the frame-0 startup divergence (every clean vector's first
subframe is now byte-exact) and gave the Annex B comfort noise its
§B.4.2.2 spectral envelope; round 455 moved the whole encoder mid-chain
(§3.3–§3.10) onto the fixed grid, sharing the decoder's primitives, and
switched the registry encoder to it.

## Round 455 — the encoder mid-chain on the fixed grid

`fx::encoder::FrameEncoderFx` drives clause 3 on the Word16/Word32
grid and shares every decoder-side primitive with `FrameDecoderFx`
(eq (40) `pred_lt3`, the Q13 codevector, the Q14/Q1 gain
reconstruction, the eq (75) landing, the `1/Â(z)` synthesis with the
overflow-rescale protocol), so the encoder's local reconstruction is
the decoder's bit for bit. The registry encoder and the `.IN` → `.BIT`
path now run it; the float chain (`encode_chain`) stays as the
spec-equation oracle.

`tests/fx_encoder_conformance.rs` measures per-parameter agreement
two ways — **locked** (every stage runs from the reference encoder's
own state, committed from the `.BIT` parameters after each search, so
each stage is measured in isolation) and **free** (end to end). Locked,
whole corpus, skeleton → round close: T1 68.3 → 80.9 %, T2 70.8 →
81.0 %, C1 45.1 → 65.3 %, GA 73.3 → 80.5 %, GB 69.6 → 78.4 %,
frame-exact 10.1 → 27.5 % (PITCH 48.3 %). Free, against the float
chain on identical metrics: exact T1 +2 to +36 points on every vector
(SPEECH 34.0 → 43.9 %), GA/GB +2 to +18, T1±2 up everywhere.

Stage by stage, each pinned black-box against the corpus:

- **§3.3** — the eq (28) LAR through the clause-5 `log2` table with
  the logarithm read as **base 10** (the clause's own `20 log g_c`
  convention): the natural-log reading tilts 43–84 % of subframes and
  loses to a forced-flat encoder on every vector; base 10 tilts
  0–26 % and beats both (PITCH locked C/S 46 → 92 %). γ₂ from the Q13
  `d_min` in radians.
- **§3.4** — Word32 sums behind the §3.2.1 overflow-rescale protocol,
  Word16-mantissa normalisation (the un-normalised `mpy_32` form loses
  15 bits and doubles the window misses), **strict** favour-lower
  compare (the printed `≥` takes the shortest section on all-zero
  frames; the reference keeps the longest). The 0.85 threshold is
  corpus-confirmed (accept/reject ratios separate at 0.845–0.856).
- **§3.7** — Word32 eq (38) rows, Word16-mantissa eq (37), eq (39)
  over **all five** interpolated candidates `k − 2/3 … k + 2/3`
  including the interpolated fraction-0 value (+4 to +8 points of T1
  everywhere); literal `k + t/3` phase geometry (the decoder's negated
  eq (40) fold halves T1 here). eq (43) through `div_s` on Q14.
- **§3.8.1** — Word16-normalised `d`/`φ`, wide `C²·E` compare: neutral
  against the exact oracle — the search is not precision sensitive;
  98 % of the remaining C/S misses are input (target / impulse
  response) differences, concentrated on quiet frames (SPEECH
  `peak < 64`: T1 22.6 % vs 96.6 % on loud frames).
- **§3.9** — candidates on the decoder's gain grid; the §3.9.2
  preselection reads the **staged threshold tables** (cluster start =
  thresholds the joint eq (63) optimum exceeds; the first GA threshold
  is exactly twice the fourth GA row's γ). Joint beats single-axis
  (−9 to −14 GA+GB) and sequential (−2 to −4) optima.

## Round 452 — frame-0 identified, SID-LSP envelope wired

Two docs stagings (the Annex B SID-LSF VQ tables and the ITU
fixed-point cluster resolutions) are consumed, and the round-419
frame-0 startup mystery is closed:

- **The frame-0 divergence was the previous-LSP start-up memory.**
  Inverting the §4.2.5 output high-pass on the frame-0 `.PST` samples
  of all six clean vectors (the base and Annex A corpora agree
  sample-for-sample there, so the recovered signal is §4.1 speech, not
  a postfilter artifact) shows the reference's frame-0 subframe-0 LP
  has `a₁ ≈ −0.6` — *not* the flat spectrum. The corpus pins the
  §3.2.5/§4.1.6 interpolation memory to the coarse decreasing vector
  `30000 26000 21000 15000 8000 0 −8000 −15000 −21000 −26000` (Q15),
  which reproduces the recovered `ŝ` on all six vectors exactly; the
  Table 9 `cos(iπ/11)` reading (the printed row is domain-invalid)
  stays available as `STARTUP_LSP_COS_GRID_Q15` but is measurably not
  what the reference decoder holds. The §3.2.4 MA-predictor reset is
  now the staged `trunc(iπ/11·2¹³)` table (five of ten entries one LSB
  below rounding).
- **§4.2/§3.9.1 grid re-pin** (the r419 pins were fitted against the
  wrong frame-0 state): synthesis lands by rounding; `1/g_f` scales
  the synthesis *input*; the §4.2.5 ×2 folds into the eq (89) AGC
  product so the output high-pass consumes a Q1 grid (integer-grid HP
  inputs cannot reproduce the reference's frame 0 at all); the HP
  feedback products are exact 48-bit landed on Q15; `ĝ_c` lands on Q1
  by truncation; the eq (75) excitation landing rounds ties toward
  zero (pinned by LSP's frame-0/1 silence).
- **Measured** (full fixed-point chain vs `.PST`, base corpus,
  exact% r419 → r452): ALGTHM 1.89 → 4.04, FIXED 19.09 → 34.45,
  LSP 3.49 → 4.01, PITCH 1.84 → 1.98, SPEECH 14.02 → 21.80 (clean
  frames 162 → 285), TAME 0.78 → 0.81, PARITY 14.38 → 26.59, ERASURE
  13.94 → 20.43; corr 0.9856–0.9999. Every clean vector's first
  divergence now sits in frame-0 subframe 1 or later (LSP: sample
  440). The residual is the §4.2.2/§4.2.3 sub-LSB operator schedule
  plus the AGC trajectory (`agc_lag` reproduces the reference's onset
  AGC exactly and is the recorded lead — see `fx::postfilter`).
- **Annex B comfort noise gains its spectral envelope** — the
  §B.4.2.2 SID-LSP dequantizer (`sid_lsf`: eq (B.18) blended
  predictor, 32-address L1 subset, 2 × 16 full-VQ L2 subset) feeds the
  §B.4.4 excitation through the interpolated SID-LSP synthesis filter
  and the §4.2 cascade, sharing state with the active chain (the CNG
  excitation enters the eq (40) history, so active frames resume
  coherently after silence). `AnnexBOutput::ComfortNoise` carries
  post-processed PCM; `AnnexBOutput::ErasedActive` carries
  §4.4-concealed PCM per §B.4.5. Measured against the `g729b`
  reference `.out`: active-frame correlation 0.9914 / 0.9886 / 0.9831
  / 0.9716 / 0.6589 / 0.9965 on tstseq1–6, comfort-noise RMS envelope
  0.72–1.06 (floors pinned; CN samples are not expected to align — the
  eq (96) draw schedule is not pinned by the prose).

## Round 438 — the encoder's fixed-point front end: TAME's L0 mystery retired

Round 438 moves the encoder's §3.1–§3.2.3 chain onto the clause-5
Word16/Word32 grid (`fx::analysis` / `fx::levinson` / `fx::lp_to_lsp`,
wired into the production `FrameEncoder`): the §3.1 eq (1) high-pass
with unrounded Word32 Q16 feedback, the eq (5) autocorrelation in a
genuine saturating Word32 accumulator with the overflow-rescale
protocol (down-scale the windowed signal on saturation, retry), DPF
lag windowing, the eqs (8)/(9) Levinson recursion in the (hi, lo)
format with `div_32` against a renormalised energy, and the §3.2.3
Chebyshev root search on Q11/Q23 grids. Every unstated fixed-point
choice was swept black-box (`tests/fx_front_end_conformance.rs`):
pinned are the overflow shift 2, the truncating Q12→Q11
half-polynomial landing, the rounding interpolation step, the Levinson
numerator pre-shift, and the Q16 feedback storage (each losing
alternative costs 6–16 points on some vector).

**Measured** (reference-locked all-four-LSP-indices agreement,
float → fx): ALGTHM 82.9 → 94.3%, FIXED 97.5%, LSP 77.7 → 85.8%,
PITCH 92.8 → 93.6%, SPEECH 81.2 → 91.7%, **TAME 80.5 → 99.2%**
(corpus 83.1 → 90.7%). End-to-end: **TAME L0 62.5 → 96.1%**, T1±2
75.0 → 82.0%, SPEECH L0/L1 87.0/53.9%, PITCH L1 35.6%. TAME's jump
retires the round-390 negative result — the "unstaged structural
element of the reference's final mode compare" was front-end ω
divergence, and no mode-compare docs-gap remains. The round also
pins the decoder's start-up LSP memory to the exact
`round(cos(iπ/11)·2^15)` grid per the staged Table 9 erratum
resolution (`table9-initialisation.md`), which confirms the printed
`arccos` row as a typographic erratum and the frame-0 startup
*schedule* as out of the Recommendation's printed scope.

## Round 419 — the fixed-point §4.2 cascade: whole decode path on the clause-5 grid

Round 419 completes the decoder on the Word16/Word32 operator grid by
landing `fx::postfilter` — the §4.2 post-processing cascade in the
clause-4.2 **printed signal order** (which differs structurally from
the float cascade's historical arrangement): `ŝ` → `Â(z/γ_n)` →
residual `r̂` → long-term postfilter `H_p(z)` **applied to the
residual** → synthesis `1/[g_f·Â(z/γ_d)]` → tilt `H_t(z)` → §4.2.4
AGC against `ŝ` → 100 Hz output high-pass + ×2. Key numerics: the
§4.2.1 two-pass search correlates a norm-scaled Word16 residual
window (Word32-safe by construction) over both sides of the integer
anchor at 1/8 resolution with the clause's longer-filter replacement
rule; `g_l` lands on Q15 by normalised division; `1/g_f`/`1/g_t` go
through a shared normalised-reciprocal helper; the AGC gain runs on
Q12 with the Q15 weight pair `27853/4915`; the §4.2.5 recursion keeps
Word32 feedback with the ×2 folded into the output rounding.

Four **black-box pins** were ratcheted by whole-frame byte-exactness
against the conformance corpus (`tests/fx_full_conformance.rs`, which
reports first-diverging-sample / longest-exact-run / clean-frame
evidence and carries an env-gated stage trace):

- **Over-unity gain guard** — an eq (83) ratio above 2 (`num >
  2·den`, the Q14 gain ceiling) behaves as *disabled*, not clamped
  to 1 (FIXED corr 0.9502/0.9756 → 0.9855/0.9918, SPEECH/PITCH
  unchanged; a per-subframe oracle probe shows passthrough beating
  the clamped filter on 42/45 FIXED enables).
- **Synthesis-stage truncation** — the `1/[g_f·Â(z/γ_d)]` output
  lands on Q0 by truncation, not rounding: SPEECH's longest
  byte-exact run jumps 77 → **9155 samples** (114 consecutive
  frames), fully-exact frames 0 → 159/3750 (g729a 166); rounding and
  toward-zero variants both collapse the clean-frame count, and a
  Word32 (unrounded) feedback memory collapses it too — the Word16
  truncated feedback is the pinned behaviour.
- **Silence enable floor** — a residual at ≤ 1 Q0 LSB mean square
  (`Σr̂² ≤ 40`) never enables `H_p` (threshold insensitive across
  40…2560); retires the silence-cluster divergence events.
- **Startup sharpening** — Table 9 prints β init 0.8, but the corpus
  pins the effective first-subframe value at the eq (47) clamp
  floor 0.2.

The **registry decoder** (`codec::G729Decoder` / `decoder::make_decoder`)
now decodes through this fixed-point chain — the float `decode_chain`
remains in-tree as the spec-equation oracle.

**Measured (full fixed-point chain vs `.PST`)**: corr 0.9855–0.9999
on all 12 clean vectors, RMS 0.92–1.05, PARITY 0.998/0.998, ERASURE
0.92/0.89, OVERFLOW 0.81 (base corpus). The divergence census is now
sharp: most vectors carry **exactly one root divergence event — at
frame 0** — whose poisoned AGC/HP state then persists (SPEECH, with
silent stretches to re-sync in, shows 9 events and 330/7500 clean
subframes). An isolation probe pins the startup mystery precisely:
overriding frame-0 subframe-0's `ŝ` with the reference-implied
`[1, 2, 1, 1, 1, 0…]` makes the whole head byte-exact, but no
reading of the printed equations produces that `ŝ` from the
transmitted parameters (`u = [1,1,1,1,0…]`, near-flat interpolated
LP) — the reference's fixed-point §4.1/§4.2 startup behaviour beyond
the printed clauses is the remaining docs-gap, along with Table 9's
domain-invalid `q_i = arccos(iπ/11)` row (arccos of values > 1).

## Round 410 — decoder conformance drive: three root causes + the fixed-point §4.1 chain

Round 410's single goal was ITU conformance exactness for the decoder.
Three long-standing structural divergences were root-caused black-box
against the conformance corpus, and the whole §4.1 chain now runs on
the clause-5 Word16/Word32 operator grid (`fx` module tree):

- **eq (74) γ̂ grid** (see below) — retired the ≈ 7–10× over-gain.
- **eq (40) fraction fold** — the two-branch `b30` read evaluates the
  past excitation at `n − k + t/3` (expand `b30(j) ≈ sinc(j/3)`: the
  coefficient of `u(m)` is `sinc(m − (n − k + t/3))` in both
  branches), so the effective delay is `k − t/3` and the transmitted
  fraction folds **negated with an upward borrow**
  (`frac = +1 → (t0+1, 2)`). The previous `T = k + t/3` fold read
  every fractional subframe 2/3 of a sample off and decorrelated
  voiced material: float-chain SPEECH corr 0.11 → 0.93, PITCH
  0.54 → 0.93 after the fix.
- **Word16 saturation inside eq (74)** — the worst-case γ column sum
  (27162 + 14276, γ̂ ≈ 5.06 — exactly the well-conditioned top of the
  γ̂ range on the 2^13 grid) overflows Word16; the fixed-point gain
  decode multiplies the two stage contributions into `g'_c`
  separately and takes the eq (72) logarithm of the 32-bit sum.
- **Fixed-point §4.1 chain** (`fx::{ops, dsp, lsp, gains, excitation,
  decoder}`) — the Table 10/11 operator set, the Table 15 log2 / pow2 /
  inv_sqrt over the staged 33/33/49-entry tables, Q13 LSF → Q15 LSP →
  Q12 LP (new `lsf-lsp-cos-table2-Q15` + `lsf-lsp-cos-slope-Q19`
  compiled), Q10 gain-predictor memory with the eq (71) exponent
  identity `K·10log10(2) = 1/2`, Q13 codevector with Q14 sharpening,
  Q15-accumulator excitation, Q12 synthesis with the specified 16-bit
  saturation, §4.1.2 parity **T1-substitution** (not erasure), the
  §4.4 concealment primitives, and the black-box-pinned
  **overflow-rescale protocol** (Word32 overflow in the synthesis
  accumulation → excitation buffer + synthesis memory ÷4 and
  re-synthesize; plain saturation measures 2.5–3.3× RMS on the
  OVERFLOW vectors, the rescale ≈ 1.0×). The eq (66) energy is pinned
  to the eq (48)-sharpened codevector (plain-pulse variant over-gains
  7–13%).
- **Measured, fx §4.1 + float §4.2 hybrid vs `.PST`**
  (`tests/fx_conformance.rs`, floors pinned): correlation
  **0.995–0.9998 on all 12 clean vectors** (SPEECH 0.9984/0.9989,
  TAME 0.9998/0.9994), RMS 0.96–1.07, bit-exact share 0.4–34%
  (FIXED 34.3%); PARITY 0.9987, ERASURE 0.91/0.94, OVERFLOW 0.78
  base (g729a decodes with the Annex A reduced decoder, not yet
  modelled: 0.25). The remaining distance to bit-exact is dominated
  by the not-yet-fixed-point §4.2 cascade and the unpublished
  reference operator schedules (docs-gaps: the §3.9.1 fixed-point
  gain-decode schedule, the §4.2 postfilter internal scaling, the
  exact overflow protocol, the §4.4 erased-frame MA-memory handling).

### The γ̂ reconstruction grid: the "§3.9 gain gap" was a scale error

- **Root cause found and fixed.** The long-documented ≈ 7–10× whole-vector
  decode over-gain (previously attributed to a missing fixed-point
  saturation of the §3.9.1 gain predictor) was actually the **eq (74)
  γ̂ reconstruction grid**: the GA/GB codebook `γ` columns sum to the
  correction factor on a **2^13 grid**, not the 2^12 grid the staged
  table annotation suggests. The Recommendation never states the codebook
  column scaling (Table 12 gives only array dimensions), so the grid was
  pinned **black-box** against the ITU conformance corpus: a 2^12 reading
  over-gains every clean vector by exactly `2^(1+Σb_i)` = 2^2.79 ≈ 6.9×
  (the γ̂ scale error compounds through the eq (69)/(72) MA feedback
  `Û^(m) = 20·log10 γ̂`; measured 7.0–9.4×), a 2^14 reading under-gains by
  the inverse (0.16–0.20×), and the 2^13 reading collapses the
  whole-corpus RMS ratio to **0.97–1.36** (TAME 2.85 base / 1.60 g729a).
  The encoder's GA preselection ranks on the same grid.
- **Exactness metrics harness** — `tests/pcm_conformance.rs` now tracks
  per-vector sample correlation / max |Δ| / bit-exact share against the
  `.PST` references with pinned correlation floors, and pins the RMS
  ratio inside a `[0.7, 3.2)` window (a single-Q-step γ̂ regression lands
  at ≈ 6.9× or ≈ 0.15×, far outside). Long vectors still decorrelate
  (SPEECH corr ≈ 0.11) because float rounding drift compounds through
  the adaptive-codebook feedback — closing that is the fixed-point
  decode path's job, not a gain-scale issue.

## Round 390 — fixed-point L0 search, Annex A decoder + encoder

- **Fixed-point Q13 §3.2.4 quantiser search** (`search_lsp_indices_q13`)
  — the L0 MA-predictor mode selection implemented exactly per the
  staged clean-room algorithm doc `docs/audio/g729/l0-mode-selection.md`
  (docs commit `b9e48a4`): eq (23) targets through the Q12 `fg_sum_inv`
  reciprocal, integer split-stage searches, eq (20) reconstruction on
  the Q13 grid, ω-domain eq (21) argmin over both predictors. The
  doc's unpinned fixed-point latitude (reconstruction-shift rounding,
  eq (21) per-term combine, tie-break) is exposed as `FxLatitude` and
  pinned black-box (320+ configs swept): prediction shift *rounds*,
  reciprocal/rearrangement shifts *truncate*, exact wide product,
  mode 0 holds ties. Locked-history all-indices agreement:
  LSP 77.5 → 77.7%, SPEECH 80.9 → 81.2% (corpus 82.9 → 83.1%);
  end-to-end FIXED L0 88.3 → 90.8%, L1 64.2 → 70.0%. **Measured
  negative result:** TAME's 24 residual L0 flips are *not* fixed-point
  near-ties — they sit in a self-sustaining 5-frame reference index
  cycle with systematic ~2.1% ω-domain margins under identical stage
  indices (flipping needs ~100 Q13-LSB input perturbations), so an
  unstaged structural element of the reference's final mode compare
  remains (docs gap; every alternative mode total — stage-sum,
  residual-full, single-folded product, `Lsp_pre_select`-style
  pre-selection — collapses corpus agreement to 34–65%).
- **Annex A decoder** (`annex_a`) — the §A.4.2 reduced-complexity
  postfilter cascade (eqs (A.11)–(A.15)): integer-delay harmonic
  postfilter (`[T_cl − 3, T_cl + 3]`, per-subframe anchor, ≤ 140), no
  `1/g_f`/`1/g_t`, length-22 tilt impulse response with
  `γ_t = 0.8/0`, numerator → tilt → synthesis order, energy-ratio
  AGC (`√(Σŝ²/Σsf²)`, 0.9/0.1). Validated on the staged `g729a`
  corpus: first-subframe deviations 2.4–4.5 PCM units (8-unit band),
  RMS ratios 7.1–10.6 under the bounded §3.9-gain-gap ceiling, shape
  distance at base-cascade parity.
- **Annex A encoder** (`annex_a_encoder`) — the §A.3 fast analysis
  chain: fixed `γ = 0.75` quantized-LP weighting
  (`W(z)/Â(z) = 1/Â(z/γ)`), eq (A.2)/(A.3) low-pass weighted speech,
  eq (A.4)/(A.5) decimated fast open-loop pitch (even-first + ±1 in
  `[80, 143]`), eq (A.7)/(A.8) fast closed-loop search
  (backward-filtered target, `b30` fractions), eq (A.10)
  filtering-free memory update. Measured vs the G.729A reference
  `.BIT` corpus: L0 62.5–94.6%, L1 33.8–100%, T1±2 57.8–86.4%
  (floors pinned); its whole SPEECH stream decodes through both
  decoders (§A.1 interoperability). The §A.3.8.1 depth-first ACELP
  pulse schedule is prose-unpinned (docs gap) — the main-body focused
  search stands in.

## Round 388 — taming, the split-stage metric, registry wiring

- **Taming procedure** (`taming` module) — the encoder-side
  adaptive-codebook gain bound, implemented from the newly-staged
  clean-room algorithm doc `docs/audio/g729/taming-procedure.md`
  (docs commit `77b8440`): four per-zone worst-case excitation-error
  accumulators over the compiled `tab_zone` partition (153 entries,
  zones split at 40/80/120), the per-subframe `tameflag` test and
  `E ← 1 + ĝ_p²·max(E_spanned)` update, and the
  `GPCLIP = 15564/2^14 ≈ 0.95` ceiling enforced inside the §3.9.2
  gain-VQ search. The doc's two unpinned constants are fixed
  black-box, and the new `reference_taming_fingerprint` test shows the
  reference **never actually tames on the staged corpus** (its own
  bitstreams keep choosing `ĝ_p > GPCLIP` up to a simulated
  accumulator of ≈ 18 186, always below the 60 000 threshold) — so our
  identical never-tame behaviour is the conforming one, retiring the
  r385 hypothesis that missing taming caused TAME's end-to-end L0 gap.
- **§3.2.4 split-stage searches on the residual-domain weighted MSE**
  — the round's headline conformance jump. The printed eq (21) is the
  ω-domain error, but the bit-exact coder's L2/L3 stages measurably
  minimise `Σ w_i·(l_i − l̂_i)²` (the same error without the
  `(1 − ΣP̂)²` weight folding — the two domains differ only by that
  factor, eq (20) being affine in `l̂`). Reference-locked
  all-four-indices agreement: **ALGTHM 71.4 → 82.9%, FIXED
  91.7 → 97.5%, LSP 74.9 → 77.5%, PITCH 88.6 → 92.8%, SPEECH
  78.4 → 80.9%**, TAME flat at 80.5%; no vector degrades. The L1 stage
  stays unweighted and the mode selection stays on the printed
  ω-domain eq (21) (both probed — each alternative collapses
  agreement). TAME's remaining locked-history misses are 24/25 pure
  `L0` mode flips — the fixed-point mode-selection arithmetic is now
  the single dominant open gap of the LSF chain.
- **Registry codec surface** (`codec` module) — `G729Decoder` /
  `G729Encoder` implementing `oxideav_core::Decoder`/`Encoder` over
  the raw 10-octet frame packing (80 Table-8 bits, MSB-first octets;
  multi-frame packets), a real `register(ctx)`, and the dual-API
  factories. Corpus-tied: all 8 100 active `.BIT` frames convert
  losslessly from the ITU serial format to the wire format and decode
  through the registry surface **sample-identical** to the
  `decode_chain` path; `reset()` restores the clause-4.3 start-up
  state byte-identically.

## Encoder fixed-point grid (round 385)

Clause 2.5 makes the 16-bit fixed-point arithmetic the conformance
ground truth (the corpus `.BIT` streams were produced by it), and the
§3.2.4 quantiser makes razor-thin nearest-neighbour decisions, so
round 385 moves the encoder's LSF chain onto the reference's numeric
grid: §3.1 pre-processing output rounded to the saturated 16-bit PCM
grid (IIR feedback keeps the unrounded recursion value), eq (4)
windowing as `⌊(s·w + 2^14)·2^−15⌋` (making the eq (5)
autocorrelation exact integer arithmetic), the Levinson output rounded
to **Q12** before the §3.2.3 root search, eq (18) evaluated through
the newly-compiled 64-segment arccos lookup (`lsf_conversion`,
`lsf-lsp-cos-table-Q15` + `lsf-lsp-acos-slope-Q12`, with
`LspQuantizer::quantize_lsf` as the LSF-domain entry point), and the
eq (7) white-noise correction at its **measured-effective unity**
(black-box sweep 1.0001 / 1+2^−13 … 1+2^−19 / 1.0: the literal 1.0001
at float precision over-inflates `r(0)` vs the 16-bit reference and
was the single largest LSF divergence). A new reference-locked
conformance harness (MA history driven by the reference's own
transmitted indices — the exact state the reference encoder had) pins
per-frame front-end fidelity separately from error propagation:
locked all-four-LSP-indices agreement is now ALGTHM 71.4% / FIXED
91.7% / LSP 74.9% / PITCH 88.6% / SPEECH 78.4% / TAME 80.5% (TAME was
**0.8%** before the Q12 step). End-to-end (own history) the ratchet
moved to L1 exact 31–100% (TAME 100%, ALGTHM 74.3%) with T1±2 78–90%.
Remaining measured gaps (black-box-probed and excluded: search
structure, window timing, root-search precision ±4 Q15 LSB,
F1/F2-coefficient quantisation, autocorrelation down-scaling, weight
thresholds/caps): the exact fixed-point eq (21)/(22) weighted
mode-selection behaviour on extreme spectra (TAME end-to-end L0 62.5%
— the reference prefers predictor 0 against a 3–10% float margin; the
TAME white-noise sweep peaking at 1+2^−17 points the same way, at a
magnitude-dependent truncation) and the residual per-frame
disagreement pointing at the DPF Levinson/normalisation precision
chain.

## Encoder (round 382)

`encode_chain::FrameEncoder` drives thirteen spec-cited encoder stages
in clause-3 order — §3.1 pre-processing (eq (1) 140 Hz high-pass, ÷2
folded in), §3.2.1 window/autocorrelation/lag-window (eqs (3)–(7)),
§3.2.2 Levinson-Durbin (eqs (8)/(9) + reflection by-product), §3.2.3
LP→LSP (Chebyshev grid search, 4× bisection), §3.2.4 LSP VQ search
(eqs (18)–(23): adaptive weights, per-mode target, staged L1/L2/L3
search over both L0 predictors, decode-side reconstructor in
lock-step), §3.3 perceptual weighting (eqs (27)–(33): LAR hysteresis,
γ₂ resonance adaptation, weighted speech), §3.4 open-loop pitch
(eqs (34)/(35) + favour-lower-delays), §3.5/§3.6 impulse response +
target (eq (36) residual, three-stage target filter), §3.7 closed-loop
pitch (eqs (37)–(39): shift-and-add recursion, 1/3-resolution b12
fractional refinement, both window procedures), §3.8.1 fixed-codebook
focused search (eqs (43)–(60): sign-folded φ, thr₃ gate K₃ = 0.4,
180-entry frame budget), §3.9.2 conjugate-structure gain VQ (eq (63)
normal-equations optimum, 4×8 preselection, decode-side scoring), and
the §3.10 memory updates. `encode_frame_to_serial` packs Table-8 and
writes the 164-byte ITU serial frame (the packer/writer re-serialises
every reference `.BIT` frame of the corpus **byte-exactly**, 8100+
frames).

Measured against the reference `.BIT` encoder outputs
(`tests/encoder_conformance.rs`, floors pinned per vector as
regression guards, round-385 numbers): frame alignment exactly 1:1,
L0 agreement 62–97%, L1 exact-match 31–100%, subframe-1 `int(T1)`
within ±2 on 78–90% of active frames; the whole 3750-frame SPEECH
vector encoded by us decodes cleanly through our own decoder.

## What is wired up

The decode path is implemented end-to-end to post-processed PCM output,
sequenced by `decode_chain::FrameDecoder` as one stateful per-frame
call. The `*_to_postfiltered` entry points run the §4.1 parameter chain,
the §3.10 / §4.1.6 LP synthesis, and the full §4.2 post-processing
cascade; the `*_to_speech` entry points return the pre-postfilter `ŝ(n)`;
and the `*_concealed` entry points add §4.4 frame-erasure concealment.
From the clause-4.3 zero state the first 5 ms subframe reconstructs to
within a few PCM units of the staged `.PST` references on every clean
vector of both corpora, and with the round-410 γ̂-grid fix the
whole-vector RMS ratio sits at ≈ 0.97–1.36 (TAME 2.85). The remaining
sample-level divergence is float-vs-16-bit rounding drift compounding
through the adaptive-codebook feedback; `tests/pcm_conformance.rs`
tracks per-vector correlation / max |Δ| / bit-exact share as the
ratchet toward the fixed-point decode path. The decode stages:

- **Bit-exact numeric-tables foundation** — the LP-analysis windowing,
  LSF cosine grid, pitch interpolation filters, MA gain-prediction
  coefficients, LSP-quantiser two-stage VQ codebooks (L1 128×10, L2/L3
  packed 32×10), the MA-predictor `fg` coefficient cube, and the
  conjugate-structure gain-VQ tables, all in their published Q-formats
  with bounds-checked lookup helpers.
- **ITU serial bitstream parser** (`serial`) — splits the 80-bit
  payload, validated against the staged conformance corpus.
- **§4.1 / Table-8 parameter unpacker** (`parameters`) — the 15 typed
  codeword indices, with the §3.7.2 pitch-parity predicate corpus-tied
  (0 mismatches on every active frame of `SPEECH.BIT`).
- **§3.2.4 LSP-frame reconstruction** (`lsp_reconstruct`) — codebook
  sum, twice-applied rearrangement, MA-prediction, and the stability
  clamp.
- **§3.2.5 per-subframe LSP interpolation** (`lsp_interpolate`) and
  **§3.2.6 LSP→LP conversion** (`lsp_to_lp`).
- **§3.9 gain decode** — `gain_index_map` (the imap1 / imap2 codeword
  permutations), `gain_reconstruct` (§3.9.2 conjugate-structure
  gain-VQ), `gain_predict` (§3.9.1 4th-order MA gain prediction).
- **§4.1.3 pitch-delay decode** (`pitch_decode`) — the transmitted
  `(P1, P2)` indices into per-subframe fractional pitch delays, with
  the symmetric forward mappings for future encoder wire-up.
- **§3.8 / §4.1.4 fixed (algebraic) codebook decode** (`fixed_codebook`)
  + **§3.6 pitch sharpening** (`pitch_sharpen`).
- **§3.7 / §3.10 / §4.1.6 LP synthesis** (`lp_synthesis`) — the
  past-excitation interpolator, the per-subframe excitation build
  `u(n) = ĝ_p·v(n) + ĝ_c·c(n)`, and the `1/Â(z)` synthesis filter.
- **§4.2.1 long-term postfilter `H_p(z)`** (`long_term_postfilter`) —
  the head of the §4.2 cascade, now the **full two-pass** form. Forms the
  eq (79) residual `r̂(n) = ŝ(n) + Σ γ_n^i·â_i·ŝ(n−i)` (the short-term
  postfilter numerator applied to `ŝ`), runs the eq (80) integer delay
  search over `[int(T1)−1, int(T1)+1]`, then the **eq (81)
  1/8-resolution fractional second pass** — refining the integer anchor
  `T_0` to `T = T_0 + frac/8` (`frac ∈ {0…7}`) by maximising the
  pseudo-normalised correlation `R(T)²/E_T`, short (length-33,
  `tab_hup_s`) filter first and the chosen non-integer fraction re-tested
  with the long (length-129, `tab_hup_l`) filter (kept only if it raises
  `R′(T)`, per the spec's longer-filter-replacement rule). The eq (82)
  `R′(T)²/Σr̂² < 0.5` disable test, the eq (83) bounded gain `g_l`, then
  the eq (78) harmonic filter `H_p(z) = (1/(1+γ_p·g_l))·(1 +
  γ_p·g_l·z⁻ᵀ)` (γ_p = 0.5) applied at the fractional delay. The two
  interpolation tables are the staged `tab_hup_s`/`tab_hup_l` CSVs
  (compiled as 7-phase polyphase kernels, phase `p`↔`8−p` mirror); the
  fractional pass engages on ~26 % of enabled subframes across the SPEECH
  corpus vector. (The earlier README "tap layout is defer-to-reference"
  docs-gap was stale — both tables were staged with full provenance.)
- **§4.2.2 short-term postfilter `H_f(z)`** (`short_term_postfilter`) —
  the eq (84) weighted pole/zero pair `Â(z/γ_n)/Â(z/γ_d)` (γ_n = 0.55,
  γ_d = 0.70) with the eq (85) impulse-response gain normalisation
  `g_f = Σ_{n=0}^{19} |h_f(n)|`, per-subframe `â_i`, continuous memory.
- **§4.2.3 tilt compensation `H_t(z)`** (`tilt_compensation`) — the
  eq (86) first-order FIR `H_t(z) = (1/g_t)·(1 + γ_t·k1'·z⁻¹)`, with
  the eq (87) tilt factor `k1' = −r_h(1)/r_h(0)` taken from the
  autocorrelation of the short-term postfilter impulse response
  `h_f(n)`, the sign-selected `γ_t` (0.9 / 0.2), and the gain term
  `g_t = 1 − |γ_t·k1'|`.
- **§4.2.4 adaptive gain control** (`adaptive_gain_control`) — the
  eq (88) energy ratio `G = Σ|ŝ(n)|/Σ|sf(n)|`, the eq (90) gain
  smoothing `g(n) = 0.85·g(n−1) + 0.15·G`, and the eq (89) scaling
  `sf′(n) = g(n)·sf(n)`, with the Table 9 init `g(−1) = 1.0`.
- **§4.2.5 output high-pass + ×2 upscaling** (`post_process`) — the
  tail of the §4.2 post-processing cascade.
- **§4.4 frame-erasure concealment** (`concealment`) — the §4.4 voicing
  classifier (`Concealer::observe_good_frame`: a frame is periodic iff
  any subframe's long-term prediction gain clears the eq (82) 3 dB
  threshold, reusing `LongTermDecision.gain`; an erased frame inherits
  the previous class), the eq (93)/(94) adaptive+fixed gain attenuations
  (`0.9`/`0.98`, the adaptive bounded `< 0.9`), the eq (95) 4 dB
  gain-predictor-memory attenuation (floored at `−14`), the eq (96)
  replacement-excitation random generator (`seed = 31821·seed + 13849`,
  seed `21845`, 13-bit index / 4-bit sign), and the §4.4.4 periodic-case
  pitch-delay repeat (`+1`/subframe, bounded `143`). These primitives are
  **wired end-to-end** into `decode_chain` via the `*_concealed` entry
  points (`decode_serial_frame_to_postfiltered_concealed` /
  `…_to_speech_concealed`): an erasure sentinel is reconstructed into
  concealed PCM (LP repeat + gain attenuation + periodic/non-periodic
  replacement excitation) instead of returning `FrameDecodeError::Erased`,
  with the voicing class latched on each good frame from the §4.2.1
  long-term decisions. The whole staged ERASURE corpus decodes to finite
  bounded PCM; the strict (non-`_concealed`) entry points keep their
  erasure-rejection contract.
- **Annex B (DTX / CNG) decoder framing + routing** (`annex_b`) — the
  silence-compression decode surface that is fully determined by the
  staged corpus + Annex B prose:
  - **Variable-length serial framing** — the Annex B `.bit` files carry
    a per-frame bit-count header (`n_bits ∈ {0, 16, 80}`) selecting
    untransmitted / SID / active; `parse_annex_b_frame` /
    `parse_annex_b_stream` classify each frame into
    `AnnexBFrame::{Active, Sid, Untransmitted, Erased}`, handling both
    §B.4.5 erasure shapes (the `0x6B20` erased-sync word and an
    all-`0x0000` payload).
  - **SID bitstream unpack** (Table B.2) — predictor (1) / first-stage
    LSF (5) / second-stage LSF (4) / gain (5) into a typed `SidFrame`.
  - **§B.4.2.1 energy dequantizer** (`dequant_sid_energy_db`) — the
    5-bit non-uniform log quantizer (−12 dB floor; −4..12 dB / 4 dB;
    16..66 dB / 2 dB), fully prose-sourced (the spec states it "does not
    need the storage of a quantizer table").
  - **§B.4.1/§B.4.5 frame-type state machine** (`AnnexBDecoder`) — the
    erasure-inheritance rule (active→active concealment, silence→
    untransmitted) and SID-parameter persistence across untransmitted
    frames.
  - **§B.4.4 comfort-noise excitation synthesis** (`cng`) — the
    silence-frame *excitation* path, fully spec-prose-sourced from the
    Annex B eqs (B.19)–(B.26): the eq (B.19) target-gain smoothing
    (`G̃_t` jump-to-SID after active, `7/8 : 1/8` relax across a
    non-active run), the eqs (B.20)–(B.22) per-subframe adaptive/fixed
    gain solving (the `Ea`/`I`/`K` reduction to the eq (B.21) monic
    quadratic via the ACELP `Σe_f²=4` identity, lowest-abs root, the
    eq (B.22) `Max{0.5, √(K/A)}` `Ga` interval), and the eqs (B.24)–(B.26)
    Gaussian mixture (`α=0.6`, `β>0` from the eq (B.25) quadratic with the
    `α=1/β=0` fallback). The eq (96) random sequence (reset per active
    frame per §B.4.4) drives the pitch-lag / `Ga` / Gaussian draws.
  - **End-to-end stream decode** (`AnnexBStreamDecoder`) — routes a whole
    Annex B `.bit` stream to per-frame output: active frames decode
    bit-exactly through the §4.1 → §4.2 base chain (`AnnexBOutput::Speech`),
    non-active frames synthesize the §B.4.4 energy-controlled comfort-noise
    *excitation* (`AnnexBOutput::ComfortNoise`, 80 samples + target gain)
    driven by the decoded SID energy. Validated over all 10 staged `g729b`
    sequences (one output block per reference PCM frame; 3 102 comfort-noise
    frames synthesized across the corpus).

### What is NOT yet wired up

- **PCM-bit-exact decode** — round 452 closed the frame-0 startup
  divergence (all six clean vectors are byte-exact through subframe 0;
  FIXED reaches 34.5% exact, SPEECH 285 fully-exact frames), but
  whole-vector exactness is still blocked by the unpublished
  §4.2.2/§4.2.3 fixed-point internal scaling (±1-LSB `st`/`tilt`
  landings) and the AGC trajectory (the `agc_lag` hook reproduces the
  reference's onset AGC exactly and is the recorded next lead). The
  float cascade (`decode_*_to_postfiltered`) remains the spec-equation
  oracle.
- **Annex B comfort-noise sample alignment** — the §B.4.2.2 SID-LSP
  envelope and §B.4.4 LP-filtered PCM are wired (round 452), but the
  eq (96) *draw schedule* (how many draws, in which order, per
  subframe) is not pinned by the prose, so CN frames match the
  reference in energy envelope, not sample-for-sample. On the encoder
  side (round 455) the §B.4.1 DTX decision and the §B.4.2 SID
  quantisation are implemented (`annex_b_encoder`; locked-VAD
  agreement on tstseq1–4: `Ftyp` 70–87 %, SID gain index 84–100 %
  exact / 100 % within one step, all four SID indices 16–42 %), but
  the §B.3 **VAD is not**: its Table B.1 boundary constants are
  printed without the parameter grids they apply to, and the §B.3.7 AR
  coefficient sets are named but never printed (docs ask). No
  Annex B registry encoder yet.
- **Bit-exact encoding** — the whole clause-3 chain runs on the
  Word16/Word32 grid (round 455: §3.3–§3.10 join the round-438
  §3.1–§3.2.3 front end), but the wire is not yet exact: locked
  (stage-isolated) frame-exact is 27.5 % corpus-wide, free 0.1 %.
  The remaining locked misses are input differences the stage
  arithmetic cannot fix — the front end's sub-LSB ω disagreements
  (LSP 85.8 %, SPEECH 91.7 % locked-history LSP agreement) and
  LSB-level target/impulse-response differences on quiet frames — plus
  the §3.9.2 preselection-boundary cases whose reference-side optimum
  computation is unpublished (see the docs asks in the round-455
  CHANGELOG). Each stage's unstated latitude is exposed on
  `fx::encoder::EncoderLatitude` for further black-box sweeps.
- **Annex A depth-first ACELP search** (§A.3.8.1) — the Annex A prose
  pins only the search's existence and fixed complexity; the
  pulse-combination schedule lives only in the barred reference C
  (docs gap). The Annex A encoder uses the main-body focused search
  (bit-stream legal, non-matching `C/S` choices). §A.4.4 concealment
  (clause 4.4 without voicing detection) is also not yet wired.
- The unnamed Annex B per-mode pair `annexB-sid-lsf-mp-Q15` is
  compiled but unconsumed (no clause names it — see the staging note);
  the gain-quantizer coefficient matrix and the
  `slope_cos`/`slope_acos` LSF↔LSP direction tables remain staged in
  `docs/` but uncompiled (no stage needs them yet).

## Clean-room provenance

The numeric data tables are transcribed from the ITU electronic
attachment to the G.729 Recommendation (a deterministic extractor reads
only the published data-table files; the SHA-256 of the attachment ZIP
and of each source data file is recorded in a `.meta` sidecar). No
algorithmic source is read by `build.rs`, `src/`, or `tests/`. The
extracted CSVs ship inside the crate so a from-crates.io build does not
depend on the workspace `docs/`; `build.rs` parses each CSV, asserts
its shape against the meta sidecar, and emits the corresponding
`pub const` array.

| Step | Artefact | Verification |
|---|---|---|
| 1 | ITU recommendation prose | manual fetch by docs collaborator; not consulted at runtime |
| 2 | ITU electronic-attachment ZIP | SHA-256 recorded in `.meta` |
| 3 | data-table files inside ZIP | each `.meta` records the source file SHA-256 + line range |
| 4 | `docs/audio/g729/tables/<name>.csv` | deterministic extractor; re-running reproduces byte-identical CSVs |
| 5 | `tables/<name>.csv` (this crate) | byte-for-byte copy of step 4, shipped with the crate |
| 6 | `pub const` arrays in `oxideav_g729::tables` | `build.rs` parses + shape-asserts the CSV and emits the declaration |

## License

MIT — see [`LICENSE`](LICENSE).
