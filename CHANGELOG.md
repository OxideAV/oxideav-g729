# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the crate adheres
to [SemVer](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
