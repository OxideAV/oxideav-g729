# oxideav-g729

Pure-Rust ITU-T **G.729** (CS-ACELP, 8 kbit/s) narrowband speech codec.
Zero C dependencies, no FFI, no `*-sys` crates.

## Status

Clean-room rebuild in progress, grown one spec-cited unit at a time
from the published ITU-T G.729 Recommendation prose and the ITU
electronic-attachment numeric data tables. The crate is currently a
**decode-side scaffold**: the registry hook (`register`) is a no-op
and the trait-surface codec entry point returns `Error::NotImplemented`
pending completion of the decode and encode chains.

## What is wired up

The decode path is implemented up to reconstructed speech, sequenced
by `decode_chain::FrameDecoder` as one stateful per-frame call:

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

### What is NOT yet wired up

- The §4.2.1 long-term postfilter (its residual-domain two-pass pitch
  search) — the only remaining §4.2 post-processing stage — and §4.4
  frame-erasure concealment.
- The full encoder.
- Registry-side codec factory wiring (the codec entry point returns
  `NotImplemented`).
- The remaining numeric tables (gain-quantizer coefficient matrix,
  postfilter interpolation, taming, Annex B DTX/CNG, LSF↔LSP
  cos/slope) are staged under `docs/audio/g729/tables/` but not yet
  compiled in.

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
