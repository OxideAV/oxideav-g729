//! # oxideav-g729
//!
//! **Status:** clean-room rebuild in progress.
//!
//! The prior implementation was retired under the workspace clean-room
//! policy: several source modules transcribed numeric tables verbatim
//! from an external reference-software distribution and described
//! matching its behaviour by citing specific source files of that
//! distribution. The clean-room policy forbids consulting any external
//! implementation's source for any reason, regardless of licensing.
//! The crate has been reset to a register-only scaffold (master
//! `676fc46`) and is being re-grown one spec-cited unit at a time.
//!
//! ## What is wired up
//!
//! Round 173 lands the bit-exact numeric-tables foundation:
//!
//! * §3.1 / §4.2 pre/post-processing IIR high-pass filter coefficients
//!   (100 Hz Q13 and 140 Hz Q12 variants).
//! * §4.1 spec Table 8 bit-allocation per analysis parameter.
//! * `basic_op()` Pow2, Log2, and Inv_sqrt 16-bit lookup tables.
//!
//! Round 191 adds a parser for the **ITU serial bitstream format**
//! used by the staged conformance corpus under
//! `docs/audio/g729/conformance/`. See [`serial`] for the 164-byte
//! frame layout (sync + bits-per-frame header + 80 bit-value words)
//! and [`serial::parse_frame`] for the decoder.
//!
//! Round 207 wires the §3.2.4 LSP-frame reconstruction algorithm
//! around the round-195 / round-201 tables — see
//! [`lsp_reconstruct`].
//!
//! Round 213 chains the §3.2.5 per-subframe LSP interpolation
//! (spec eq (24), linear interpolation in the cosine domain) on
//! top of the round-207 LSF output — see [`lsp_interpolate`].
//!
//! Round 218 wires the §3.2.6 LSP-to-LP conversion (spec eqs (13)
//! / (14) F1/F2 sum/difference polynomial recursion, eq (25) `(1 ±
//! z^-1)` factor restoration, eq (26) `A(z) = (F'_1 + F'_2) / 2`
//! recombination) on top of the round-213 per-subframe cosine-
//! domain LSPs — see [`lsp_to_lp`].
//!
//! Round 225 wires the §4.1 / Table-8 parameter unpacker that
//! splits the round-191 serial 80-bit payload into the 15 typed
//! codeword indices (`L0` / `L1` / `L2` / `L3`, `P1` / `P0`,
//! `C1` / `S1` / `GA1` / `GB1`, `P2`, `C2` / `S2` / `GA2` / `GB2`)
//! the §4.1 decode procedure consumes — see [`parameters`].
//!
//! Round 231 wires the §3.9.2 / §4.1.5 conjugate-structure gain-VQ
//! reconstruction: the GA (8 × 2 Q14/Q12) / GB (16 × 2 Q14/Q12)
//! codebooks land in [`tables`] alongside their §3.9.3
//! transmission-side permutation / inverse-permutation /
//! partial-search-threshold tables, and a new
//! [`gain_reconstruct`] module maps each transmitted (GA, GB)
//! index pair into the quantised `(ĝ_p, γ̂)` pair via spec
//! eqs (73) / (74).
//!
//! Round 239 wires the §3.9.1 / §4.1.5 4th-order MA gain prediction
//! stage that turns the round-231 correction factor `γ̂` into the
//! actual quantised fixed-codebook gain `ĝ_c = γ̂ · g'_c`. A new
//! [`gain_predict`] module holds the stateful 4-tap predictor (spec
//! Table 9 / clause 4.3 init `Û^(k) = -14`), the eq (66) codevector
//! energy `E = 10·log10((1/40)·Σ_{n=0..39} c(n)^2)`, the eq (69)
//! 4-tap MA sum `Ẽ^(m) = Σ_{i=1..=4} b_i · Û^(m−i)`, the eq (71)
//! exponential `g'_c = 10^((Ẽ^(m) + Ē − E)/20)` with `Ē = 30 dB`,
//! and the eq (72) decode-form history advance
//! `Û^(m) = 20·log10(γ̂)`.
//!
//! Round 249 wires the §3.9.3 gain-quantiser codeword-mapping layer
//! that sits between the round-225 transmitted indices and the
//! round-231 codebook lookup. A new [`gain_index_map`] module exposes
//! [`gain_index_map::demap_ga`] / [`gain_index_map::demap_gb`] (the
//! decode-side `imap1` / `imap2` inverse-permutation primitives that
//! map the on-wire GA / GB indices back into the `0..NCODE1` /
//! `0..NCODE2` codebook domain), the symmetric encode-side
//! [`gain_index_map::map_ga`] / [`gain_index_map::map_gb`] forward
//! permutations, and the per-frame [`gain_index_map::demap_frame`]
//! wrapper that yields the four codebook indices in one call. The
//! existing [`gain_reconstruct::reconstruct_frame_gains`] entry
//! point now chains the demap step before the §3.9.2 lookup so the
//! `(GA, GB) → (ĝ_p, γ̂)` pipeline is spec-conformant end-to-end
//! from the on-wire bits.
//!
//! Round 255 wires the §4.1.3 pitch-delay decode that maps the
//! transmitted `(P1, P2)` indices into per-subframe fractional
//! pitch delays `(T1, T2)`. A new [`pitch_decode`] module exposes
//! [`pitch_decode::decode_t1_from_p1`] (spec image `f0027-01.jpg`
//! — `if P1 < 197: int(T1) = (P1+2)/3 + 19`, `frac = P1 −
//! 3·int(T1) + 58`, else `int(T1) = P1 − 112`, `frac = 0`),
//! [`pitch_decode::derive_t_min`] (spec image `f0027-02.jpg` —
//! the ±5 / [20, 143] subframe-2 search-window derivation),
//! [`pitch_decode::decode_t2_from_p2`] (spec image `f0027-03.jpg`
//! — `int(T2) = (P2+2)/3 − 1 + t_min`, `frac = P2 − 2 −
//! 3·((P2+2)/3 − 1)`), and the per-frame
//! [`pitch_decode::decode_frame`] wrapper that chains them in
//! spec §4.1.3 order against the round-225
//! [`parameters::Parameters`] struct.
//!
//! Round 266 wires the §3.8 / §4.1.4 fixed (algebraic) codebook
//! decode that maps the transmitted `(C, S)` codewords into the
//! per-subframe pulse positions, signs, and the 40-sample
//! codevector `c(n)`. A new [`fixed_codebook`] module exposes
//! [`fixed_codebook::decode_positions`] (spec eq (62) — the
//! 13-bit `C` codeword splits as 3+3+3+4 across the four Table-7
//! tracks, with the 4-bit track-3 field carrying `2·(m_3/5) + jx`),
//! [`fixed_codebook::decode_signs`] (spec eq (61) — the 4-bit `S`
//! codeword carries one sign bit per pulse, `s_k = 1` ⇔ positive),
//! [`fixed_codebook::build_codevector`] (spec eq (45) — four signed
//! unit impulses on the 40-sample subframe), and the per-frame
//! [`fixed_codebook::decode_frame`] wrapper. The §3.8 eq (48)
//! pitch sharpening (when `int(T) < 40`) needs the round-255
//! pitch delay AND the previous subframe's quantised
//! adaptive-codebook gain and is deferred to a follow-up round.
//!
//! Round 274 lands that deferred §3.8 / §4.1.4 **pitch sharpening**
//! step. A new [`pitch_sharpen`] module exposes
//! [`pitch_sharpen::clamp_beta`] (spec eq (47) — `β = ĝ_p^(m−1)`
//! clamped to `[0.2, 0.8]`, the previous subframe's quantised
//! adaptive-codebook gain) and [`pitch_sharpen::sharpen`] (spec
//! eq (48) — the in-place recurrence `c(n) += β·c(n − T)` for
//! `n = T … 39`, applied only when the current subframe's integer
//! pitch delay `int(T) < 40`, realising the §3.8 adaptive
//! pre-filter `P(z) = 1/(1 − β·z^−T)` on the round-266
//! algebraic codevector). The recurrence reads the already-modified
//! `c(n − T)` in a forward sweep so the geometric pitch echo
//! propagates correctly; [`pitch_sharpen::codevector_energy`] reads
//! the post-sharpening `Σ c(n)²`.
//!
//! See [`tables`] for the full inventory and Q-format conventions.
//!
//! ## What is NOT wired up
//!
//! Every decode/encode entry point still returns
//! [`Error::NotImplemented`]. The remaining numeric tables
//! (gain-quantizer coefficient matrix, postfilter interpolation,
//! taming, Annex B DTX/CNG, LSF↔LSP cos/slope tables) are staged
//! under `docs/audio/g729/tables/` but not yet compiled in; the
//! Implementer leaves them out until the docs collaborator's
//! specifier pass clarifies the per-clause wire-up direction.
//!
//! ## Clean-room provenance
//!
//! All numeric tables under [`tables`] are generated by `build.rs`
//! from CSVs under `tables/<spec-name>.csv`, which are byte-for-byte
//! copies of the spec-role-named outputs under
//! `docs/audio/g729/tables/`. The extractor that produced those CSVs
//! reads only data-only `.c` table files staged alongside the spec;
//! per-file provenance (origin filename, byte ranges, file-level
//! SHA-256) is recorded in each CSV's `.meta` sidecar under
//! `docs/audio/g729/tables/`. No algorithmic source is consulted.

#![warn(missing_debug_implementations)]

use oxideav_core::RuntimeContext;

pub mod fixed_codebook;
pub mod gain_index_map;
pub mod gain_predict;
pub mod gain_reconstruct;
pub mod lsp_interpolate;
pub mod lsp_reconstruct;
pub mod lsp_to_lp;
pub mod parameters;
pub mod pitch_decode;
pub mod pitch_sharpen;
pub mod serial;
pub mod tables;

/// Crate-local error type. Until decode + encode are wired up, every
/// public codec entry point returns [`Error::NotImplemented`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// The crate has been reset to a scaffold pending clean-room
    /// rebuild; no decoder or encoder functionality is wired up yet.
    NotImplemented,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "oxideav-g729: clean-room rebuild in progress — decoder/encoder not yet wired"
        )
    }
}

impl std::error::Error for Error {}

/// No-op codec registration — the rebuild scaffold registers nothing
/// runnable into the runtime context yet. The hook is preserved so
/// downstream consumers can keep their `register!` graphs intact.
pub fn register(_ctx: &mut RuntimeContext) {}

oxideav_core::register!("g729", register);
