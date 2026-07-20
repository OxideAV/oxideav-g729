//! §4.1.3 / §4.1.4 / §3.10 / §4.1.6 fixed-point excitation and
//! synthesis: the eq (40) adaptive-vector interpolation over the past
//! excitation, the algebraic codevector on the Q13 grid with eq (48)
//! pitch sharpening, the eq (75) excitation build, and the eq (77)
//! `1/Â(z)` synthesis filter on Q12 coefficients.
//!
//! ## Number grids
//!
//! - **Excitation / speech** — Q0 (the 16-bit PCM grid; the encoder
//!   input was halved by §3.1, so the decoder's reconstruction is
//!   also on the halved grid until §4.2.5 restores ×2).
//! - **Adaptive vector** `v(n)` — Q0, from the Q15 31-tap `b30`
//!   interpolator ([`crate::tables::PITCH_INTERP_FILTER_SYNTHESIS_Q15`])
//!   with a rounding accumulate.
//! - **Codevector** `c(n)` — Q13 (±8192 pulses); the eq (48) sharpening
//!   gain is held on Q14 clamped to `[0.2, 0.8]` (`3277..=13107`).
//! - **Gains** — `ĝ_p` Q14, `ĝ_c` Q1 ([`crate::fx::gains`]), so the
//!   eq (75) accumulator `ĝ_p·v + ĝ_c·c` lands on Q15 and rounds to
//!   Q0 after one left shift.
//! - **Synthesis** — Q12 `â` against Q0 memory; the accumulator holds
//!   the Q12 sum which shifts up 3 and rounds to Q0. Saturation here
//!   is *the specified behaviour* (the OVERFLOW conformance vector
//!   exists to exercise it).
//!
//! ## The T < 40 repetition
//!
//! Clause 3.7: "for delays less than the subframe length the
//! excitation is repeated". The interpolation window reads the
//! excitation buffer *in place* as `v(n)` is written, so for short
//! delays the already-built adaptive samples of the current subframe
//! are what gets repeated (the gain-scaled `u(n)` replaces them only
//! after the gains decode). This in-place convention is the natural
//! fixed-point buffer layout; the float chain reads the gain-scaled
//! `u(n')` instead — a semantic difference for `T < 40` that the
//! conformance corpus arbitrates.

use crate::fx::ops::{add, l_mac, l_shl, mult, round, shl};
use crate::tables::{M, PITCH_INTERP_FILTER_SYNTHESIS_Q15};

/// Samples per subframe.
pub const L_SUBFR: usize = 40;
/// Past-excitation history length: `T_max (143) + fractional margin`.
pub const PIT_MAX: usize = 143;
/// Interpolator half-length (10 taps per branch).
const L_INTER10: usize = 10;
/// Total excitation buffer: history + one frame of two subframes.
pub const EXC_BUF: usize = PIT_MAX + L_INTER10 + 2 * L_SUBFR;
/// eq (48) sharpening gain clamp, Q14 (`0.2 … 0.8`).
pub const SHARP_MIN_Q14: i16 = 3277;
/// See [`SHARP_MIN_Q14`].
pub const SHARP_MAX_Q14: i16 = 13107;
/// Start-up sharpening gain — the eq (47) clamp floor (0.2 on Q14).
/// Table 9 lists β's non-zero init as 0.8, but the conformance corpus
/// pins the *effective* first-subframe value at the 0.2 floor (the
/// 0.8 reading regresses first-divergence positions on the LSP/TAME
/// vectors), i.e. the eq (47) adaptation from the (zero-initialised)
/// previous quantized pitch gain runs before the first use.
pub const SHARP_INIT_Q14: i16 = SHARP_MIN_Q14;

/// eq (40): build the adaptive vector `v(n)` in place over the
/// excitation buffer. `exc` points at the whole buffer; `off` is the
/// current subframe start (`off >= PIT_MAX + L_INTER10`); `t0` and
/// `frac ∈ {-1, 0, 1}` are the decoded integer/fractional delay.
///
/// The interpolation reads `u` at `n − t0 − i` and `n − t0 + 1 + i`
/// after the fraction is folded to the non-negative side; both
/// branches use the Q15 `b30` filter with phase stride 3.
pub fn pred_lt3(exc: &mut [i16; EXC_BUF], off: usize, t0: i32, frac: i32) {
    let b30 = &PITCH_INTERP_FILTER_SYNTHESIS_Q15;
    // Fraction fold. The two-branch read below evaluates the
    // band-limited past excitation at position `n − k + t/3` (write
    // `b30[j] ≈ sinc(j/3)` and expand either branch: the coefficient
    // of `u(m)` is `sinc(m − (n − k + t/3))` in both), so the phase
    // axis runs *against* the delay axis: the effective delay is
    // `k − t/3`. `T = t0 + frac/3` therefore folds as `frac = −1 →
    // (t0, 1)` and `frac = +1 → (t0 + 1, 2)` — the transmitted
    // fraction is negated and the borrow goes *up*. (The un-negated
    // fold reads every fractional subframe 2/3 of a sample off and
    // decorrelates voiced material; pinned black-box on the corpus.)
    let (k, t) = match frac {
        0 => (t0, 0usize),
        -1 => (t0, 1),
        _ => (t0 + 1, 2), // frac == +1
    };
    for n in 0..L_SUBFR {
        let base = off + n - k as usize;
        let mut acc = 0i32;
        for i in 0..L_INTER10 {
            // Σ u(n−k−i)·b30(t+3i) + Σ u(n−k+1+i)·b30(3−t+3i)
            acc = l_mac(acc, exc[base - i], b30[t + 3 * i]);
            acc = l_mac(acc, exc[base + 1 + i], b30[3 - t + 3 * i]);
        }
        exc[off + n] = round(acc);
    }
}

/// §4.1.4: the Q13 algebraic codevector from decoded pulse positions
/// and signs, plus the eq (48) pitch sharpening at delay `t` with the
/// Q14 sharpening gain.
#[must_use]
pub fn build_codevector_q13(
    positions: &[u8; 4],
    signs: &[i8; 4],
    sharp_q14: i16,
    t: i32,
) -> [i16; L_SUBFR] {
    let mut code = [0i16; L_SUBFR];
    for p in 0..4 {
        let pos = usize::from(positions[p]);
        let pulse: i16 = if signs[p] >= 0 { 8192 } else { -8192 };
        code[pos] = add(code[pos], pulse);
    }
    // eq (48): c(n) += β·c(n−T) for n >= T (only when T < 40).
    if t > 0 && (t as usize) < L_SUBFR {
        let t = t as usize;
        for n in t..L_SUBFR {
            // Q13 × Q14 → Q12, one left shift back to Q13.
            let contrib = shl(mult(code[n - t], sharp_q14), 1);
            code[n] = add(code[n], contrib);
        }
    }
    code
}

/// eq (75): overwrite the adaptive samples of the current subframe
/// with the full excitation `u(n) = ĝ_p·v(n) + ĝ_c·c(n)`.
///
/// `v(n)` is read from the buffer (as left by [`pred_lt3`]); the Q14
/// pitch gain and Q1 code gain land the doubled accumulator on Q16,
/// which rounds to the Q0 grid.
pub fn build_excitation(
    exc: &mut [i16; EXC_BUF],
    off: usize,
    code_q13: &[i16; L_SUBFR],
    gain_pit_q14: i16,
    gain_code_q1: i16,
) {
    for n in 0..L_SUBFR {
        // v·ĝ_p: Q0×Q14 → Q15; c·ĝ_c: Q13×Q1 → Q15.
        let mut acc = l_mac(0, exc[off + n], gain_pit_q14);
        acc = l_mac(acc, code_q13[n], gain_code_q1);
        acc = l_shl(acc, 1);
        exc[off + n] = round(acc);
    }
}

/// eq (77): the `1/Â(z)` synthesis filter over one subframe.
///
/// `a` is `[a_0 = 4096, a_1 … a_10]` on Q12; `mem` holds the last
/// `M` output samples (`mem[M-1]` = most recent) and advances.
/// The Q12 accumulator shifts up 3 and rounds onto the Q0 grid;
/// saturation on loud material is the specified 16-bit behaviour.
pub fn syn_filt(exc: &[i16], a: &[i16; M + 1], mem: &mut [i16; M], out: &mut [i16]) {
    let _ = syn_filt_check(exc, a, mem, out, true);
}

/// [`syn_filt`] with 32-bit overflow detection: runs the eq (77)
/// filter and reports whether any accumulator step left the Word32
/// range or the final `<< 3` rescale would (the saturating behaviour
/// still applies to the produced output). When `update_mem` is false
/// the filter memory is left untouched — callers dry-run the filter
/// to decide on the overflow rescale before committing state.
pub fn syn_filt_check(
    exc: &[i16],
    a: &[i16; M + 1],
    mem: &mut [i16; M],
    out: &mut [i16],
    update_mem: bool,
) -> bool {
    debug_assert_eq!(exc.len(), out.len());
    let mut overflowed = false;
    for n in 0..exc.len() {
        let mut acc = l_mac(0, exc[n], a[0]); // u(n)·a0 on Q13
        let mut wide = i64::from(acc);
        for i in 1..=M {
            let s = if n >= i { out[n - i] } else { mem[M - (i - n)] };
            acc = crate::fx::ops::l_msu(acc, a[i], s);
            wide -= 2 * i64::from(a[i]) * i64::from(s);
        }
        if wide != i64::from(acc) || wide.abs() >= (1i64 << 28) {
            overflowed = true;
        }
        acc = l_shl(acc, 3);
        out[n] = round(acc);
    }
    if update_mem {
        let n = exc.len();
        mem[..M].copy_from_slice(&out[n - M..n]);
    }
    overflowed
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A pure integer delay (frac = 0) with a delta-history buffer
    /// reproduces the history sample (the b30 filter is 1 at phase 0,
    /// lag 0 — b30[0] ≈ 0.9 with sidelobes; the full tap set sums to
    /// ~1 for the centred branch). We assert reconstruction of a
    /// constant history within the filter's DC ripple.
    #[test]
    fn pred_lt3_integer_delay_tracks_constant_history() {
        let mut exc = [0i16; EXC_BUF];
        for slot in exc.iter_mut().take(PIT_MAX + L_INTER10) {
            *slot = 1000;
        }
        let off = PIT_MAX + L_INTER10;
        pred_lt3(&mut exc, off, 60, 0);
        // With a constant past and T=60 > 40, v(n) is the constant
        // filtered by the DC gain of b30 (≈ 1.0).
        for n in 0..L_SUBFR {
            let v = exc[off + n];
            assert!(
                (i32::from(v) - 1000).abs() <= 10,
                "v({n}) = {v} for constant history"
            );
        }
    }

    /// Fractional phases interpolate between neighbours: for a ramp
    /// history, T = 60 + 1/3 lands between the T = 60 and T = 61
    /// integer reads.
    #[test]
    fn pred_lt3_fraction_between_integers() {
        let mut base = [0i16; EXC_BUF];
        for (i, slot) in base.iter_mut().enumerate().take(PIT_MAX + L_INTER10) {
            *slot = (i as i16) * 10;
        }
        let off = PIT_MAX + L_INTER10;

        let mut e0 = base;
        pred_lt3(&mut e0, off, 60, 0);
        let mut e1 = base;
        pred_lt3(&mut e1, off, 60, 1); // T = 60 + 1/3 → reads ~1/3 older
        let mut e2 = base;
        pred_lt3(&mut e2, off, 61, 0);

        for n in 0..8 {
            let a = i32::from(e0[off + n]);
            let b = i32::from(e1[off + n]);
            let c = i32::from(e2[off + n]);
            let lo = a.min(c) - 12;
            let hi = a.max(c) + 12;
            assert!(
                (lo..=hi).contains(&b),
                "fractional read out of band at {n}: {a} / {b} / {c}"
            );
        }
    }

    /// Codevector: four ±8192 pulses at the decoded positions, and
    /// eq (48) sharpening adds β·c(n−T) for T < 40.
    #[test]
    fn codevector_pulses_and_sharpening() {
        let code = build_codevector_q13(&[0, 11, 22, 33], &[1, -1, 1, -1], 13107, 20);
        assert_eq!(code[0], 8192);
        assert_eq!(code[11], -8192);
        // Sharpening at T=20: c(20) += 0.8·c(0) → 8192·0.8 ≈ 6553.
        assert_eq!(i32::from(code[20]), i32::from(shl(mult(8192, 13107), 1)));
        assert_eq!(code[22], add(8192, 0)); // 22−20 = 2 has no pulse
                                            // n=31: c(31) += 0.8·c(11) = −6553…
        assert!(code[31] < -6000);
        // No sharpening for T ≥ 40.
        let plain = build_codevector_q13(&[0, 11, 22, 33], &[1, -1, 1, -1], 13107, 40);
        assert_eq!(plain[20], 0);
    }

    /// eq (75) grid landing: v = 1000 (Q0), ĝ_p = 0.5 (Q14 8192),
    /// c = 8192 (Q13 = 1.0), ĝ_c = 100 (Q1 200) → u = 600.
    #[test]
    fn excitation_grids_combine() {
        let mut exc = [0i16; EXC_BUF];
        let off = PIT_MAX + L_INTER10;
        exc[off] = 1000;
        let mut code = [0i16; L_SUBFR];
        code[0] = 8192;
        build_excitation(&mut exc, off, &code, 8192, 200);
        assert_eq!(exc[off], 600);
    }

    /// Synthesis filter matches the float filter within 1 LSB on a
    /// moderate excitation (both implement eq (77); only the Q12/Q0
    /// rounding differs).
    #[test]
    fn syn_filt_tracks_float() {
        // A stable Â(z): a1 = −0.5 → Q12 −2048.
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2048;
        let exc: [i16; L_SUBFR] = std::array::from_fn(|n| if n == 0 { 4000 } else { 0 });
        let mut mem = [0i16; M];
        let mut out = [0i16; L_SUBFR];
        syn_filt(&exc, &a, &mut mem, &mut out);
        // Float: s(n) = u(n) + 0.5·s(n−1) → 4000, 2000, 1000, …
        let mut expect = 4000.0f64;
        for (n, &o) in out.iter().enumerate().take(10) {
            assert!(
                (f64::from(o) - expect).abs() <= 1.5,
                "s({n}) = {o} expected {expect}"
            );
            expect *= 0.5;
        }
        // Memory carries the last M outputs.
        assert_eq!(mem[M - 1], out[L_SUBFR - 1]);
    }

    /// The synthesis accumulator saturates rather than wraps on hot
    /// input — the OVERFLOW-vector behaviour.
    #[test]
    fn syn_filt_saturates() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -4096; // pole at z = 1 → unbounded growth
        let exc = [20000i16; L_SUBFR];
        let mut mem = [0i16; M];
        let mut out = [0i16; L_SUBFR];
        syn_filt(&exc, &a, &mut mem, &mut out);
        assert_eq!(out[L_SUBFR - 1], 32767);
    }
}
