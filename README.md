# oxideav-g729

Pure-Rust ITU-T **G.729** (CS-ACELP, 8 kbit/s) narrowband speech codec —
decoder + encoder, with optional **Annex B** voice-activity detection,
DTX, and comfort-noise generation. Zero C dependencies, no FFI, no
`*-sys` crates.

Part of the [oxideav](https://github.com/OxideAV/oxideav-workspace)
framework but usable standalone.

## Installation

```toml
[dependencies]
oxideav-core = "0.1"
oxideav-codec = "0.1"
oxideav-g729 = "0.0"
```

## Format

- **Sample rate**: 8 kHz mono (`S16`), narrowband telephony.
- **Frame**: 10 ms = 80 samples.
- **Bitstream**: 80 bits = 10 bytes per speech frame, ITU-T Table 8
  field order (L0, L1, L2, L3, P1, P0, C1, S1, GA1, GB1, P2, C2, S2,
  GA2, GB2).
- **Algorithm**: CS-ACELP — 10th-order LPC with MA-4 split-VQ of LSPs,
  fractional (1/3) pitch in the adaptive codebook, 4-pulse ACELP
  fixed codebook, two-stage gain VQ, short-term + pitch + tilt
  postfilter with AGC.

## Quick use

Both directions go through the shared `oxideav-codec` registry. The
same `CodecId("g729")` covers encoder and decoder.

```rust
use oxideav_codec::CodecRegistry;
use oxideav_core::{CodecId, CodecParameters, SampleFormat};

let mut codecs = CodecRegistry::new();
oxideav_g729::register(&mut codecs);

let mut params = CodecParameters::audio(CodecId::new("g729"));
params.sample_rate = Some(8_000);
params.channels = Some(1);
params.sample_format = Some(SampleFormat::S16);

let mut enc = codecs.make_encoder(&params)?;
let mut dec = codecs.make_decoder(&params)?;
# Ok::<(), oxideav_core::Error>(())
```

`send_frame` accepts 8 kHz mono `S16` audio in arbitrary-sized chunks;
the encoder internally queues samples and emits one 10-byte `Packet`
per 80-sample frame. The decoder takes those packets back and returns
80-sample `S16` `AudioFrame`s with monotonic PTS.

## Annex B (VAD / DTX / CNG)

Annex B silence compression is **opt-in per encoder instance**. Enable
it by setting a single-byte extradata marker on the parameters you
pass to `make_encoder`:

```rust
use oxideav_g729::ANNEX_B_ENABLE_EXTRADATA;
params.extradata = vec![ANNEX_B_ENABLE_EXTRADATA];
```

With Annex B on, the encoder classifies every 10 ms frame and emits
one of three packet shapes:

| Shape    | Size     | Meaning                                        |
|----------|----------|------------------------------------------------|
| Speech   | 10 bytes | Normal CS-ACELP frame.                         |
| SID      | 2 bytes  | Silence Insertion Descriptor (spectrum+gain).  |
| NODATA   | 0 bytes  | Silence continuation (DTX).                    |

The decoder dispatches on `packet.data.len()` and routes SID/NODATA
to the comfort-noise generator. A NODATA packet that arrives before
the decoder has seen any SID is rejected as malformed, matching the
spec's DTX rule. Disabling Annex B makes the encoder emit only
10-byte speech frames and the decoder reject 0/2-byte packets.

The VAD here is a simplified energy-plus-hangover detector rather
than the spec's four-feature decision tree — adequate for
speech/silence routing in most voice-over-IP contexts, but not
bit-exact against the ITU reference implementation.

## Implementation scope

- Decoder: full CS-ACELP pipeline (LSP decode, LSP interpolation,
  LSP to LPC, adaptive + fixed codebook excitation, gain VQ,
  synthesis, short-term + pitch + tilt postfilter with AGC, and the
  §4.2.5 output stage — 100 Hz high-pass `H_h2(z)` + ×2 level restore).
- Encoder: Levinson-Durbin LPC analysis, Chebyshev-based LPC to LSP,
  split-VQ LSP quantisation, §3.3 perceptual-weighting filter `W(z)` +
  §3.5 weighted-synthesis-filter impulse response `h(n)`, §3.4 open-loop
  pitch analysis (three-range correlation + sub-multiple bias) anchoring
  the §3.7 closed-loop fractional-lag search, focused 4-pulse ACELP
  fixed-codebook search, two-stage gain VQ, bit-exact packer for
  ITU-T Table 8.
- Annex B: VAD + DTX + CNG (simplified, interoperates with this crate's
  own decoder).

### Known deviations from ITU-T G.729 (2007)

- **LSP quantisation: spec-exact.** `LSPCB1_Q13` (128×10),
  `LSPCB2_Q13` (32×10), `FG_Q15` (2×4×10), `FG_SUM_Q15`, and
  `FG_SUM_INV_Q12` are transcribed verbatim from the ITU reference C
  source `TAB_LD8K.C` (G.729 Software Package Release 2, November
  2006). The encoder split-VQ search and the decoder
  `Lsp_get_quant`/`Lsp_prev_compose` paths use the canonical layout
  (L2 → low-half columns of the row, L3 → high-half columns of the
  row) per `LSPGETQ.C`.
- **MA-4 gain prediction: spec-exact.** Both encoder and decoder
  carry a 4-deep `Û^(k)` history of past quantised prediction errors,
  initialised to `-14 dB` per ITU-T Table 9, and apply the MA
  coefficients `[0.68, 0.58, 0.34, 0.19]` from §3.9.1 eq (69). The
  encoder quantises the spec's correction factor `γ = g_c / g'_c`
  (eq 72) — not the raw fixed-codebook gain — and updates its history
  with `Û^(m) = 20·log10(γ_q)` (eq 70) in lockstep with the decoder.
- **Decoder post-processing: spec-exact.** The §4.2 chain now matches
  the spec: short-term postfilter `Â(z/γ_n)/Â(z/γ_d)` (γ_n = 0.55,
  γ_d = 0.7), pitch emphasis, tilt compensation `(1/g_t)(1 + γ_t·k1'·z⁻¹)`
  with `g_t = 1 − |k1'|` (eq 86) and the L = 20 reflection-coefficient
  estimate (eq 87), adaptive gain control using the L1-norm gain ratio
  `G = Σ|ŝ| / Σ|sf|` (eq 88) with `g(n) = 0.85·g(n-1) + 0.15·G` (eq 90),
  and the final §4.2.5 output stage: the 100 Hz high-pass filter
  `H_h2(z)` (eq 91, coefficients verbatim) followed by the ×2 level
  restore that undoes the encoder's §3.1 ×0.5 preprocessing.
- **LPC analysis window: spec-exact.** The encoder runs the spec's
  §3.2.1 / eq (3) 240-sample asymmetric window — half-Hamming
  `0.54 − 0.46·cos(2π·n/399)` over n = 0…199 plus a quarter-cosine
  fade `cos(2π·(n−200)/159)` over n = 200…239 — applied to a buffer
  of 120 past samples + the 80-sample current frame + 40 samples of
  look-ahead from the following frame. The autocorrelation
  `r(0)`-floor (1.0), white-noise correction (×1.0001), and 60 Hz
  bandwidth lag window from eqs (6, 7) all match the spec. This
  realises the spec's 5 ms extra algorithmic delay at the encoder.
- **Open-loop pitch analysis (§3.4): spec-faithful.** The encoder now
  runs the once-per-frame open-loop pitch search on the perceptually-
  weighted speech `sw(n)` of eq (33): the three-range correlation
  `R(k) = Σ sw(n)·sw(n-k)` (eq 34) over the three disjoint delay bands
  (80…143, 40…79, 20…39), normalised by `sqrt(Σ sw²(n-k))` (eq 35), with
  the sub-multiple bias decision tree (`R'(t_i) ≥ 0.85·R'(T_op)` favours
  shorter delays to avoid latching onto a pitch multiple). The resulting
  `T_op` anchors the §3.7 closed-loop search: subframe 1 searches the
  six-sample window `[T_op−3, T_op+3]` (eq f0018-01), subframe 2 the
  ten-sample window `[int(T1)−5, int(T1)+4]` (eq f0018-02), both with the
  spec's boundary clamps. Lives in `src/open_loop_pitch.rs`, unit-tested.
- **Perceptual weighting: spec-faithful building blocks; γ-adaptation
  pending.** The §3.3 weighting filter `W(z) = A(z/γ1)/A(z/γ2)` (eq 27)
  and the §3.5 weighted-synthesis-filter impulse response `h(n)` of
  `A(z/γ1)/[Â(z)·A(z/γ2)]` are implemented in `src/weighting.rs` exactly
  per spec. The adaptive (γ1, γ2) derivation building blocks —
  log-area-ratio coefficients (eq 28), flat/tilted hysteresis classifier
  (eq 30), minimum-LSP-distance `d_min` (eq 31), and `γ2 = −6·d_min + 1`
  bounded to [0.4, 0.7] (eq 32) — are present and unit-tested, but the
  per-frame LAR-based flat/tilted *selection* (eqs 28–30 fed by the
  Levinson-Durbin reflection coefficients) is not yet wired: the §3.4
  weighted-speech path uses the §3.3 "flat" gammas (0.94, 0.6) as a
  documented default. The analysis-by-synthesis fixed-codebook search
  still drives off the raw LP residual rather than convolving candidates
  with `h(n)`; wiring `h(n)` into the §3.8 search loop is a later step.
- **Gain VQ table values: NOT spec-exact.** The table *dimensions*
  match the spec (`GBK1` = 8×2, `GBK2` = 16×2 per §5.2 Table 12), and
  the codebook structure is correct — both rows hold `(ĝ_p, γ̂)` and
  the decoded gain is `g_c = γ̂ · g'_c` (eq 74). The numeric entries,
  however, are first-cut values approximating the span of the
  reference `gbk1` / `gbk2`; verbatim transcription is pending.
- **Annex B VAD: simplified.** Energy-plus-hangover detector rather
  than the spec's four-feature decision tree.

Net effect: **encode → decode round-trips cleanly inside this crate**
(exercised by `tests/encoder_roundtrip.rs`, including a 5-second
steady-tone test asserting the gain predictor stays stable across
the full duration). LSP quantisation, MA-4 gain prediction, and the
LP-analysis windowing pipeline are spec-faithful; full external
interoperability awaits verbatim `gbk1` / `gbk2` numeric tables —
the encoder/decoder structure is ready for drop-in replacement of
those tables, no logic changes required.

### Annexes

- **Annex A** (reduced-complexity): the implementation already lives
  in the Annex A complexity bracket — no separate selector exposed.
- **Annex B** (VAD / DTX / CNG): supported, opt-in via extradata.
- **Annex D** (6.4 kbit/s), **Annex E** (11.8 kbit/s), **Annex G**
  (embedded variable bit rate): not supported. Only the 8 kbit/s
  main body is decoded and encoded.

## Codec id

- `"g729"` — registered as both a decoder and an encoder on the
  shared registry via `oxideav_g729::register`.

## License

MIT — see [LICENSE](LICENSE).
