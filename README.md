# oxideav-g729

Pure-Rust ITU-T **G.729** (CS-ACELP, 8 kbit/s) narrowband speech codec.
Zero C dependencies, no FFI, no `*-sys` crates.

## Status

Clean-room rebuild in progress, grown one spec-cited unit at a time
from the published ITU-T G.729 Recommendation prose and the ITU
electronic-attachment numeric data tables. The crate now carries
**both directions**: the decode path end-to-end to post-processed PCM,
and (round 382) the **entire clause-3 encoder analysis chain** with a
working `.IN` → `.BIT` path. The registry hook (`register`) is still a
no-op and the trait-surface codec entry point returns
`Error::NotImplemented` pending registry wire-up.

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
(`tests/encoder_conformance.rs`, float vs fixed-point, floors pinned as
regression guards): frame alignment exactly 1:1, L0 agreement 81–95%,
L1 exact-match 33–80%, subframe-1 `int(T1)` within ±2 on 56–91% of
active frames; the whole 3750-frame SPEECH vector encoded by us decodes
cleanly through our own decoder.

## What is wired up

The decode path is implemented end-to-end to post-processed PCM output,
sequenced by `decode_chain::FrameDecoder` as one stateful per-frame
call. The `*_to_postfiltered` entry points run the §4.1 parameter chain,
the §3.10 / §4.1.6 LP synthesis, and the full §4.2 post-processing
cascade; the `*_to_speech` entry points return the pre-postfilter `ŝ(n)`;
and the `*_concealed` entry points add §4.4 frame-erasure concealment.
From the clause-4.3 zero state the first 5 ms subframe reconstructs to
within a few PCM units of the staged `.PST` references on every clean
vector of both corpora. Beyond the first frame the amplitude is a
**bounded** multiplicative over-gain (whole-vector RMS ratio ≈ 7–10×):
the float path omits the fixed-point reference's Q-format saturation of
the §3.9.1 gain predictor (a documented docs-gap — the §3.9 gain
fixed-point saturation), so the predicted energy `Ẽ^(m)` over-predicts
without diverging. `tests/pcm_conformance.rs` measures and bounds this
gap as a regression guard. The decode stages:

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

- **PCM-bit-exact decode** — beyond the first frame the amplitude is a
  bounded multiplicative over-gain (whole-vector RMS ratio ≈ 7–10×, up to
  ≈ 18× on the TAME vector) because the float path omits the fixed-point
  reference's Q-format saturation of the §3.9.1 gain predictor. The spec
  prose (clause 3.9.1, eqs (66)–(72)) specifies **no** such clamp — it is
  reference-code-only, so this is a hard docs-gap (the §3.9 gain
  fixed-point saturation). `tests/pcm_conformance.rs` bounds the gap as a
  regression guard. The §4.2 cascade — including the §4.2.1 two-pass
  long-term postfilter wired this round — is chained end-to-end into
  `decode_*_to_postfiltered`; the `*_to_speech` entry points
  intentionally return the pre-postfilter reconstructed speech `ŝ(n)`.
- **Annex B comfort-noise spectral envelope** (§B.4.2.2 SID-LSP VQ
  dequant) — blocked on absent numeric tables. The SID-LSP VQ subset
  codebooks (the §B.4.2.2 32-address first-stage subset + the two
  16-address second-stage subset + the modified MA predictor) and
  `annexB-cng-lsp-sid-reset-Q15.csv` (`lspSid_reset`, listed in the docs
  `tables/README.md` but **not present** as a file) are not staged under
  `docs/audio/g729/tables/` — only the CNG spectrum factor/shift and VAD
  margin tables are. The §B.4.4 CNG *excitation* (energy path) is now
  synthesized (eqs (B.19)–(B.26), `cng` module); what still awaits the
  SID-LSP VQ tables is the **LP-filtered PCM** — the excitation must pass
  through the SID-LSP-derived synthesis filter for the correct spectral
  colour, so `AnnexBOutput::ComfortNoise` surfaces the raw excitation
  until that envelope can be reconstructed. (Docs-gap: SID-LSP VQ subset
  codebooks + `lspSid_reset`.)
- **Bit-exact encoding** — the round-382 encoder is float-domain; its
  per-frame decisions agree with the fixed-point reference at the
  measured rates above but are not yet bit-exact (razor-thin float
  comparisons in every analysis stage). Closing this requires the same
  fixed-point Q-format arithmetic as the §3.9 decode gap.
- Registry-side codec factory wiring (the codec entry point returns
  `NotImplemented`).
- The remaining numeric tables (gain-quantizer coefficient matrix,
  taming, Annex B SID-LSP VQ, LSF↔LSP cos/slope) are staged under
  `docs/audio/g729/tables/` but not yet compiled in. (The §4.2.1
  postfilter interpolation filters `tab_hup_s`/`tab_hup_l` are now
  compiled — see the long-term postfilter above.)

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
