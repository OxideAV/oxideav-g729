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
  synthesis, short-term + pitch + tilt postfilter with AGC).
- Encoder: Levinson-Durbin LPC analysis, Chebyshev-based LPC to LSP,
  split-VQ LSP quantisation, open-loop + closed-loop fractional-lag
  pitch search, focused 4-pulse ACELP fixed-codebook search, two-stage
  gain VQ, bit-exact packer for ITU-T Table 8.
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
- **Gain VQ: NOT spec-exact.** Two-stage tables `GBK1` / `GBK2` are
  reduced first-cut tables (8 + 16 entries) covering the span of the
  spec's `gbk1` / `gbk2`; the full 8×8 + 16×8 ITU tables are pending
  transcription. The MA-4 gain predictor uses a uniform-tap mean
  rather than the spec's predictor coefficients.
- **LPC analysis window: approximation.** Hamming window in place of
  the spec's 240-sample asymmetric window (60-sample lookahead).
- **Annex B VAD: simplified.** Energy-plus-hangover detector rather
  than the spec's four-feature decision tree.

Net effect: **encode → decode round-trips cleanly inside this crate**
(exercised by `tests/encoder_roundtrip.rs`); LSP quantisation
sub-pipeline is bit-exact against the ITU reference tables and
matches the canonical `Lsp_get_quant` formula. Full external
interoperability awaits the verbatim gain-VQ tables and the
240-sample analysis window — the encoder/decoder structure is ready
for drop-in replacement, no logic changes required.

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
