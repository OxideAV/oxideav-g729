//! §4.2 fixed-point post-processing cascade: the long-term (harmonic)
//! postfilter (§4.2.1, eqs (78)–(83)), the short-term postfilter
//! (§4.2.2, eqs (84)/(85)), tilt compensation (§4.2.3, eqs (86)/(87)),
//! adaptive gain control (§4.2.4, eqs (88)–(90)), and the output
//! high-pass + ×2 upscaling (§4.2.5, eq (91)) — everything on the
//! clause-5 Word16/Word32 operator grid, consuming the §4.1
//! [`crate::fx::decoder::FrameDecoderFx`] output. The Annex A
//! reduced-complexity cascade (§A.4.2) shares the module.
//!
//! ## Operator schedule
//!
//! The §4.2.2–§4.2.4 and §A.4.2.4 schedules follow the staged
//! per-sample traces (`docs/audio/g729/traces/`, written up in
//! `g729-errata-and-clarifications.md` #440) operator for operator;
//! `tests/postfilter_trace_conformance.rs` replays every staged
//! sub-frame against this module column by column.
//!
//! - **Weighted LP coefficients** `γ^i·â_i` — Q12, the γ powers built
//!   recursively on Q15 **with rounding** ([`weight_az`]).
//! - **eq (79) residual** — `L_mult`/`L_mac` Q13 accumulator, `L_shl 3`,
//!   `round`.
//! - **eq (85) `g_f`** — `Σ|h_f|` over 20 Q12 taps, truncated to Q10;
//!   applied only when `g_f > 1`, as `mult_r(x, div_s(1024, g_f))` on
//!   the **input** of the all-pole section.
//! - **`1/Â(z/γ_d)`** — rounded Word16 recursion memory.
//! - **eq (87) `k′_1`** — common `norm_l` of `r_h(0)`, truncated high
//!   words, `div_s`.
//! - **eq (86) tilt** — numerator with the halved `γ_t·k′_1` (LSB
//!   cleared) rounded half-up, then a second, separate rounding for
//!   `1/g_t` (Q14 or Q11 by the sign of `k′_1`).
//! - **eqs (88)–(90) AGC** — `G` on Q14 through normalised `div_s`;
//!   the recursion runs with the Q15 pair `32358/410` (the reference
//!   pole is 0.9875, not the printed 0.85) through `mult_r`, from
//!   `n = 0`; `sf′(n) = round(L_shl(L_mult(g(n), sf(n)), 1))`.
//! - **§A.4.2.4 AGC** — energy ratio over quarter-scaled samples, the
//!   clause-5 inverse square root, `G` on Q12, truncating `mult`
//!   recursion with the printed `0.9` pole, truncated output product;
//!   a zero output energy resets the gain state to 0 and passes the
//!   sub-frame unscaled.
//! - **§4.2.5** — the staged Q13 coefficient tables against a Word32
//!   feedback memory carried in full.
//!
//! The §4.2.1 delay search runs on a justified copy of the residual
//! (peak below `2^12`, so every 40-term Word32 correlation stays in
//! range); its compare arithmetic is exact cross-multiplication of
//! those Word32 values, and every staged decision (delay, 1/8 phase,
//! long/short kernel, `1/(1 + γ_p·g_l)` weight) reproduces.

use crate::fx::dsp::inv_sqrt;
use crate::fx::ops::{
    abs_s, add, div_s, extract_h, extract_l, l_abs, l_add, l_deposit_l, l_mac, l_msu, l_mult,
    l_shl, l_shr, mult, mult_r, negate, norm_l, norm_s, round, sature, shl, shr, sub,
};
use crate::tables::{
    HPF_PREPROC_100HZ_A_Q13, HPF_PREPROC_100HZ_B_Q13, M, POSTFILTER_FRAC_RES,
    POSTFILTER_INTERP_LONG_TAPS, POSTFILTER_INTERP_SHORT_TAPS, postfilter_interp_long,
    postfilter_interp_short,
};

/// Samples per subframe.
pub const L_SUBFR: usize = 40;

/// §4.2.2 numerator weight `γ_n = 0.55` on Q15.
pub const GAMMA_N_Q15: i16 = 18022;

/// §4.2.2 denominator weight `γ_d = 0.70` on Q15.
pub const GAMMA_D_Q15: i16 = 22938;

/// §4.2.3 tilt weight `γ_t = 0.9` (non-positive `k′_1`) on Q15.
pub const GAMMA_T_NEG_Q15: i16 = 29491;

/// §4.2.3 tilt weight `γ_t = 0.2` (positive `k′_1`) on Q15.
pub const GAMMA_T_POS_Q15: i16 = 6554;

/// eq (90) previous-gain weight on Q15 — `0.9875`, the value the
/// staged traces establish (the printed `0.85` does not occur).
pub const AGC_PREV_Q15: i16 = 32358;

/// eq (90) target weight on Q15 (`32358 + 410 = 2^15`).
pub const AGC_TARGET_Q15: i16 = 410;

/// §4.2.4 unit gain on the Q14 grid (Table 9 `g(−1) = 1.0`).
pub const AGC_UNITY_Q14: i16 = 16384;

/// §A.4.2.4 previous-gain weight `0.9` on Q15.
pub const AGC_PREV_ANNEX_A_Q15: i16 = 29491;

/// §A.4.2.4 target weight on Q15 (`32767 − 29491`).
pub const AGC_TARGET_ANNEX_A_Q15: i16 = 3276;

/// §A.4.2.4 unit gain on the Q12 grid (Table 9 `g(−1) = 1.0`).
pub const AGC_UNITY_ANNEX_A_Q12: i16 = 4096;

/// eq (85) impulse-response truncation length.
pub const GF_IMPULSE_LEN: usize = 20;

/// eq (A.14) impulse-response truncation length (Annex A).
pub const ANNEX_A_IMPULSE_LEN: usize = 22;

/// §A.4.2.3 tilt weight `γ_t = 0.8` on Q15.
pub const GAMMA_T_ANNEX_A_Q15: i16 = 26214;

/// eq (78) weight of the undelayed residual at `g_l = 1`:
/// `1/(1 + γ_p)` on Q15.
pub const LT_UNDELAYED_MIN_Q15: i16 = 21845;

/// Deepest §4.2.1 history reach: the longest integer delay the search
/// can pick (`int(T_1) + 1 ≤ 144`) plus the long-interpolator margin
/// (8 samples beyond the integer anchor).
const HIST: usize = 144 + POSTFILTER_INTERP_LONG_TAPS / 2;

/// One subframe's §4.2.1 decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LtDecisionFx {
    /// Integer anchor of the chosen delay: `T = delay + frac/8`.
    pub delay: usize,
    /// Chosen 1/8-resolution fractional offset (`0 … 7`).
    pub frac: usize,
    /// eq (83) gain `g_l` on Q15 (`0` = disabled, eq (82)).
    pub gain_q15: i16,
    /// Whether the length-129 refinement replaced the length-33
    /// delayed signal (clause 4.2.1 keeps the longer-filter signal
    /// only when it raises `R′(T)`).
    pub use_long: bool,
    /// eq (78) weight `1/(1 + γ_p·g_l)` of the undelayed residual on
    /// Q15 (`32767` when the filter is disabled).
    pub undelayed_q15: i16,
}

impl LtDecisionFx {
    const fn off(delay: usize) -> Self {
        Self {
            delay,
            frac: 0,
            gain_q15: 0,
            use_long: false,
            undelayed_q15: i16::MAX,
        }
    }
}

/// `γ`-weighted LP coefficients `γ^i·â_i` on Q12 (`i = 1 … 10`,
/// slot `j` holds `i = j + 1`): the γ powers are built recursively on
/// Q15 with rounding (`γ^i = round(L_mult(γ^(i−1), γ))`) and each
/// product lands by `round(L_mult(â_i, γ^i))`.
#[must_use]
pub fn weight_az(a: &[i16; M + 1], gamma_q15: i16) -> [i16; M] {
    let mut fac = gamma_q15;
    std::array::from_fn(|j| {
        let w = round(l_mult(a[j + 1], fac));
        if j + 1 < M {
            fac = round(l_mult(fac, gamma_q15));
        }
        w
    })
}

/// The truncated impulse response `h_f(n)` (Q12) of the un-normalised
/// short-term postfilter `Â(z/γ_n)/Â(z/γ_d)` (eq (84) without `1/g_f`).
#[must_use]
pub fn impulse_response_q12(apn: &[i16; M], apd: &[i16; M]) -> [i16; GF_IMPULSE_LEN] {
    let mut h = [0i16; GF_IMPULSE_LEN];
    impulse_response_into(apn, apd, &mut h);
    h
}

/// Impulse response of `Â(z/γ_n)/Â(z/γ_d)` on Q12 into `h` (any
/// length): the all-pole recursion with the rounded Word16 landing,
/// from zero state, over the numerator sequence
/// `{4096, γ_n·â_1 … γ_n^10·â_10, 0 …}`.
fn impulse_response_into(apn: &[i16; M], apd: &[i16; M], h: &mut [i16]) {
    for n in 0..h.len() {
        let x = if n == 0 {
            4096
        } else if n <= M {
            apn[n - 1]
        } else {
            0
        };
        let mut acc = l_mult(x, 4096);
        for (i, &c) in apd.iter().enumerate().take(n) {
            acc = l_msu(acc, c, h[n - 1 - i]);
        }
        h[n] = round(l_shl(acc, 3));
    }
}

/// Stateful fixed-point §4.2 cascade.
///
/// Signal path (clause 4.2 introduction): `ŝ(n)` → `Â(z/γ_n)` →
/// residual `r̂(n)` → long-term postfilter `H_p(z)` (on the residual)
/// → synthesis `1/[g_f·Â(z/γ_d)]` → tilt `H_t(z)` → `sf(n)` → AGC
/// against `ŝ(n)` → output high-pass + ×2.
#[derive(Debug, Clone)]
pub struct PostfilterFx {
    /// §4.2.2-numerator input history `[ŝ(n−1) … ŝ(n−10)]` for the
    /// eq (79) residual.
    s_hist: [i16; M],
    /// §4.2.1 residual buffer: `HIST` history samples then the current
    /// subframe (`r_buf[HIST + n]` = `r̂(n)`).
    r_buf: [i16; HIST + L_SUBFR],
    /// Synthesis `1/Â(z/γ_d)` output history `[y(n−1) … y(n−10)]`.
    st_y: [i16; M],
    /// Tilt-filter input tap `x(−1)`: the last sample of the previous
    /// subframe's tilt *input*.
    tilt_prev: i16,
    /// §4.2.4 smoothed gain `g(n)`: Q14 (main body) or Q12 (Annex A).
    agc_gain: i16,
    /// §4.2.5 input taps `[x(n−1), x(n−2)]`.
    hp_x: [i16; 2],
    /// §4.2.5 unrounded feedback `[y(n−1), y(n−2)]` on Word32.
    hp_y: [i32; 2],
    /// Oracle-probe switch: apply the §4.2.1 filter as disabled.
    force_lt_off: bool,
    /// §A.4.2 reduced-complexity cascade (see [`Self::new_annex_a`]).
    annex_a: bool,
}

impl Default for PostfilterFx {
    fn default() -> Self {
        Self::new()
    }
}

/// A candidate delay's correlation / delayed-signal energy pair on
/// the justified grid.
#[derive(Debug, Clone, Copy)]
struct LtCand {
    anchor: usize,
    frac: usize,
    corr: i32,
    energy: i32,
}

impl PostfilterFx {
    /// Clause-4.3 start-up state (all-zero except the Table 9 AGC gain).
    #[must_use]
    pub fn new() -> Self {
        Self {
            s_hist: [0; M],
            r_buf: [0; HIST + L_SUBFR],
            st_y: [0; M],
            tilt_prev: 0,
            agc_gain: AGC_UNITY_Q14,
            hp_x: [0; 2],
            hp_y: [0; 2],
            force_lt_off: false,
            annex_a: false,
        }
    }

    /// The **Annex A** (§A.4.2) reduced-complexity cascade on the same
    /// grid: the §A.4.2.1 long-term postfilter searches **integer**
    /// delays only, in `[T_cl − 3, T_cl + 3]` around the *current*
    /// subframe's transmitted delay (`T_cl ≤ 140`); the §A.4.2.2
    /// short-term postfilter drops `1/g_f`; the §A.4.2.3 tilt filter
    /// `1 + γ_t·k′_1·z⁻¹` drops `1/g_t`, takes `k′_1` from the
    /// length-22 impulse response (eq (A.14)) and uses `γ_t = 0.8`
    /// for negative `k′_1`, `0` otherwise; the tilt compensation runs
    /// **before** the `1/Â(z/γ_d)` synthesis (§A.4.2); the §A.4.2.4
    /// AGC target is the energy ratio `√(Σŝ²/Σsf²)` (eq (A.15)).
    /// [`Self::process_subframe`]'s delay argument is then `T_cl`.
    #[must_use]
    pub fn new_annex_a() -> Self {
        let mut pf = Self::new();
        pf.annex_a = true;
        pf.agc_gain = AGC_UNITY_ANNEX_A_Q12;
        pf
    }

    /// Fetch at `buf(n − k)` (`n` = in-subframe index, `k ≥ 0`):
    /// buffer index `HIST + n − k`, zero beyond the history.
    #[inline]
    fn at(buf: &[i16; HIST + L_SUBFR], n: usize, k: usize) -> i16 {
        let idx = HIST + n;
        if k > idx { 0 } else { buf[idx - k] }
    }

    /// eq (79): the residual `r̂(n) = ŝ(n) + Σ γ_n^i·â_i·ŝ(n−i)` for the
    /// current subframe, written into `r_buf[HIST..]`.
    fn residual(&mut self, s: &[i16; L_SUBFR], apn: &[i16; M]) {
        for n in 0..L_SUBFR {
            let mut acc = l_mult(s[n], 4096); // Q13
            for (i, &c) in apn.iter().enumerate() {
                let sv = if n > i {
                    s[n - 1 - i]
                } else {
                    self.s_hist[i - n]
                };
                acc = l_mac(acc, c, sv);
            }
            self.r_buf[HIST + n] = round(l_shl(acc, 3));
        }
    }

    /// Interpolated fetch of `buf` at delay `T_0 + frac/8` for output
    /// index `n`, with the short (`long == false`) or long kernel.
    /// `frac == 0` is the direct integer fetch.
    fn interp_at(buf: &[i16; HIST + L_SUBFR], n: usize, t0: usize, frac: usize, long: bool) -> i16 {
        if frac == 0 {
            return Self::at(buf, n, t0);
        }
        let (taps, half): (&[i16], usize) = if long {
            (
                postfilter_interp_long(frac),
                POSTFILTER_INTERP_LONG_TAPS / 2,
            )
        } else {
            (
                postfilter_interp_short(frac),
                POSTFILTER_INTERP_SHORT_TAPS / 2,
            )
        };
        let mut acc = 0i32;
        for (j, &h) in taps.iter().enumerate() {
            let k = (t0 + half).saturating_sub(j);
            acc = l_mac(acc, h, Self::at(buf, n, k));
        }
        round(acc)
    }

    /// The whole delayed subframe for one candidate.
    fn delayed(buf: &[i16; HIST + L_SUBFR], t0: usize, frac: usize, long: bool) -> [i16; L_SUBFR] {
        std::array::from_fn(|n| Self::interp_at(buf, n, t0, frac, long))
    }

    /// Correlation `Σ r̂(n)·r̂_T(n)` and energy `Σ r̂_T(n)²` of a delayed
    /// subframe against the justified residual.
    fn corr_energy(rs: &[i16; HIST + L_SUBFR], rk: &[i16; L_SUBFR]) -> (i32, i32) {
        let mut corr = 0i32;
        let mut energy = 0i32;
        for n in 0..L_SUBFR {
            corr = l_mac(corr, rs[HIST + n], rk[n]);
            energy = l_mac(energy, rk[n], rk[n]);
        }
        (corr, energy)
    }

    /// The justified copy of the residual window the §4.2.1 decisions
    /// run on, and its shift (`> 0` = right): the peak magnitude of
    /// the whole window lands below `2^12`, so every 40-term Word32
    /// correlation stays in range. `None` for an all-zero window.
    fn justified_residual(&self) -> Option<([i16; HIST + L_SUBFR], i16)> {
        let mut maxabs: i16 = 0;
        for &v in self.r_buf.iter() {
            let a = abs_s(v);
            if a > maxabs {
                maxabs = a;
            }
        }
        if maxabs == 0 {
            return None;
        }
        let sh = sub(3, norm_s(maxabs));
        Some((std::array::from_fn(|i| shr(self.r_buf[i], sh)), sh))
    }

    /// `a²/b > c²/d` for positive Word32 values, exact.
    #[inline]
    fn ratio_gt(a: i32, b: i32, c: i32, d: i32) -> bool {
        i128::from(a) * i128::from(a) * i128::from(d)
            > i128::from(c) * i128::from(c) * i128::from(b)
    }

    /// eq (82) disable test + eq (83) gain → the eq (78) weight of the
    /// undelayed residual, `1/(1 + γ_p·g_l)` on Q15. `None` = disabled.
    fn lt_gain(corr: i32, energy: i32, energy_r: i32) -> Option<(i16, i16)> {
        if corr <= 0 || energy <= 0 {
            return None;
        }
        // eq (82): R′(T)²/Σr̂² < 0.5 ⇔ 2·corr² < Σr̂²·E_T.
        if 2 * i128::from(corr) * i128::from(corr) < i128::from(energy_r) * i128::from(energy) {
            return None;
        }
        if corr >= energy {
            // eq (83) bound g_l ≤ 1.
            return Some((LT_UNDELAYED_MIN_Q15, i16::MAX));
        }
        let nc = norm_l(corr);
        let ne = norm_l(energy);
        let c16 = extract_h(l_shl(corr, nc));
        let e16 = extract_h(l_shl(energy, ne));
        let num = shr(c16, sub(nc, ne));
        let half_den = shr(e16, 1);
        let g0 = div_s(half_den, add(half_den, shr(num, 2)));
        let g_l = if num >= e16 {
            i16::MAX
        } else {
            div_s(num, e16)
        };
        Some((g0, g_l))
    }

    /// Clause-4.2.1 two-pass delay search + eq (82)/(83) gain decision
    /// for the residual currently in `r_buf`; returns the decision and
    /// the delayed residual on the output grid.
    fn decide(&self, int_t1: usize) -> (LtDecisionFx, [i16; L_SUBFR]) {
        let zero = [0i16; L_SUBFR];
        let lo = int_t1.saturating_sub(1).max(1);
        let hi = (int_t1 + 1).min(144);
        let Some((rs, sh)) = self.justified_residual() else {
            return (LtDecisionFx::off(lo), zero);
        };

        // eq (82) denominator: Σ r̂(n)² on the justified grid.
        let mut energy_r = 0i32;
        for n in 0..L_SUBFR {
            energy_r = l_mac(energy_r, rs[HIST + n], rs[HIST + n]);
        }

        // First pass (eq (80)): best integer delay in
        // [int(T_1) − 1, int(T_1) + 1].
        let mut t0 = lo;
        let mut best_corr = i32::MIN;
        for k in lo..=hi {
            let mut corr = 0i32;
            for n in 0..L_SUBFR {
                corr = l_mac(corr, rs[HIST + n], Self::at(&rs, n, k));
            }
            if corr > best_corr {
                best_corr = corr;
                t0 = k;
            }
        }
        if best_corr <= 0 || energy_r <= 0 {
            return (LtDecisionFx::off(t0), zero);
        }

        // Second pass (eq (81)): best 1/8 fraction by R(T)²/E_T (short
        // kernel), positive-correlation candidates only. "Around T_0"
        // spans both sides: `T_0 − 7/8 … T_0 + 7/8`, the negative side
        // expressed as `(T_0 − 1) + frac/8`.
        let mut best: Option<(LtCand, [i16; L_SUBFR])> = None;
        let neg_anchor = t0.saturating_sub(1).max(1);
        let candidates = (1..POSTFILTER_FRAC_RES)
            .map(|frac| (neg_anchor, frac))
            .chain((0..POSTFILTER_FRAC_RES).map(|frac| (t0, frac)));
        for (anchor, frac) in candidates {
            let rk = Self::delayed(&rs, anchor, frac, false);
            let (corr, energy) = Self::corr_energy(&rs, &rk);
            if corr <= 0 || energy <= 0 {
                continue;
            }
            let better = match &best {
                None => true,
                Some((b, _)) => Self::ratio_gt(corr, energy, b.corr, b.energy),
            };
            if better {
                best = Some((
                    LtCand {
                        anchor,
                        frac,
                        corr,
                        energy,
                    },
                    rk,
                ));
            }
        }
        let Some((mut cand, mut rk)) = best else {
            return (LtDecisionFx::off(t0), zero);
        };

        // Long-kernel refinement of a non-integer winner, kept only if
        // it raises R(T)²/E_T.
        let mut use_long = false;
        if cand.frac != 0 {
            let rl = Self::delayed(&rs, cand.anchor, cand.frac, true);
            let (corr_l, e_l) = Self::corr_energy(&rs, &rl);
            if corr_l > 0 && e_l > 0 && Self::ratio_gt(corr_l, e_l, cand.corr, cand.energy) {
                cand.corr = corr_l;
                cand.energy = e_l;
                rk = rl;
                use_long = true;
            }
        }

        let Some((g0, g_l)) = Self::lt_gain(cand.corr, cand.energy, energy_r) else {
            return (LtDecisionFx::off(t0), zero);
        };
        // The delayed residual on the output grid: the integer fetch
        // is exact; an interpolated one is computed on the justified
        // grid and shifted back.
        let out: [i16; L_SUBFR] = if cand.frac == 0 {
            std::array::from_fn(|n| Self::at(&self.r_buf, n, cand.anchor))
        } else {
            std::array::from_fn(|n| shl(rk[n], sh))
        };
        (
            LtDecisionFx {
                delay: cand.anchor,
                frac: cand.frac,
                gain_q15: g_l,
                use_long,
                undelayed_q15: g0,
            },
            out,
        )
    }

    /// §A.4.2.1: best integer delay in `[T_cl − 3, T_cl + 3]` by the
    /// eq (80) correlation, then the eq (82)/(83) gain decision on it
    /// (no fractional pass, no longer-filter refinement).
    fn decide_integer(&self, t_cl: usize) -> (LtDecisionFx, [i16; L_SUBFR]) {
        let zero = [0i16; L_SUBFR];
        let centre = t_cl.min(140);
        let lo = centre.saturating_sub(3).max(1);
        let hi = (centre + 3).min(143);
        let Some((rs, _sh)) = self.justified_residual() else {
            return (LtDecisionFx::off(lo), zero);
        };
        let mut energy_r = 0i32;
        for n in 0..L_SUBFR {
            energy_r = l_mac(energy_r, rs[HIST + n], rs[HIST + n]);
        }
        let mut t0 = lo;
        let mut best_corr = i32::MIN;
        for k in lo..=hi {
            let mut corr = 0i32;
            for n in 0..L_SUBFR {
                corr = l_mac(corr, rs[HIST + n], Self::at(&rs, n, k));
            }
            if corr > best_corr {
                best_corr = corr;
                t0 = k;
            }
        }
        let mut energy = 0i32;
        for n in 0..L_SUBFR {
            let rk = Self::at(&rs, n, t0);
            energy = l_mac(energy, rk, rk);
        }
        let Some((g0, g_l)) = Self::lt_gain(best_corr, energy, energy_r) else {
            return (LtDecisionFx::off(t0), zero);
        };
        (
            LtDecisionFx {
                delay: t0,
                frac: 0,
                gain_q15: g_l,
                use_long: false,
                undelayed_q15: g0,
            },
            std::array::from_fn(|n| Self::at(&self.r_buf, n, t0)),
        )
    }

    /// §4.2.1 `H_p(z)` for one subframe, applied to the **residual**
    /// (clause 4.2: "the signal r̂ is then filtered through the
    /// long-term postfilter"): the delay search + the eq (78)
    /// combination `g_0·r̂(n) + (1 − g_0)·r̂_T(n)` with
    /// `g_0 = 1/(1 + γ_p·g_l)`; advances the residual history. The
    /// residual must already sit in `r_buf[HIST..]`.
    fn long_term(&mut self, int_t1: usize) -> ([i16; L_SUBFR], LtDecisionFx) {
        let (mut d, delayed) = if self.annex_a {
            self.decide_integer(int_t1)
        } else {
            self.decide(int_t1)
        };
        if self.force_lt_off {
            d = LtDecisionFx::off(d.delay);
        }

        let mut out = [0i16; L_SUBFR];
        if d.gain_q15 == 0 {
            out.copy_from_slice(&self.r_buf[HIST..]);
        } else {
            let g0 = d.undelayed_q15;
            let g1 = add(sub(i16::MAX, g0), 1);
            for (n, o) in out.iter_mut().enumerate() {
                *o = round(l_mac(l_mult(g0, self.r_buf[HIST + n]), g1, delayed[n]));
            }
        }

        // Advance the residual history by one subframe.
        self.r_buf.copy_within(L_SUBFR.., 0);
        (out, d)
    }

    /// The `1/Â(z/γ_d)` all-pole section with the rounded Word16
    /// recursion memory (continuous across subframes).
    fn synthesis(&mut self, x: &[i16; L_SUBFR], apd: &[i16; M]) -> [i16; L_SUBFR] {
        let mut y = [0i16; L_SUBFR];
        for n in 0..L_SUBFR {
            let mut acc = l_mult(x[n], 4096);
            for (i, &cf) in apd.iter().enumerate() {
                let yv = if n > i {
                    y[n - 1 - i]
                } else {
                    self.st_y[i - n]
                };
                acc = l_msu(acc, cf, yv);
            }
            y[n] = round(l_shl(acc, 3));
        }
        for i in 0..M {
            self.st_y[i] = y[L_SUBFR - 1 - i];
        }
        y
    }

    /// eq (87): `k′_1 = −r_h(1)/r_h(0)` on Q15 from the Q12 impulse
    /// response — Word32 `L_mac` sums, the common `norm_l` of
    /// `r_h(0)`, truncated high words, `div_s`.
    fn tilt_factor(h: &[i16]) -> i16 {
        let mut rh0 = 0i32;
        let mut rh1 = 0i32;
        for j in 0..h.len() {
            rh0 = l_mac(rh0, h[j], h[j]);
            if j + 1 < h.len() {
                rh1 = l_mac(rh1, h[j], h[j + 1]);
            }
        }
        let s = norm_l(rh0);
        let r0 = extract_h(l_shl(rh0, s));
        let r1 = extract_h(l_shl(rh1, s));
        let a1 = abs_s(r1);
        if r0 <= 0 || r0 < a1 {
            return 0;
        }
        let k = div_s(a1, r0);
        if r1 > 0 { negate(k) } else { k }
    }

    /// §4.2.3 `H_t(z)`: `sf(n) = (1/g_t)·(x(n) + γ_t·k′_1·x(n−1))`.
    fn tilt(&mut self, x: &[i16; L_SUBFR], k1_q15: i16) -> [i16; L_SUBFR] {
        let mu = if k1_q15 > 0 {
            mult_r(k1_q15, GAMMA_T_POS_Q15)
        } else {
            mult_r(k1_q15, GAMMA_T_NEG_Q15)
        };
        // 1 − |γ_t·k′_1| on Q15 (saturating at 32767 for μ = 0), and
        // 1/g_t on Q14 (positive k′_1) or Q11.
        let d = add(i16::MAX, sub(1, abs_s(mu)));
        let (dividend, out_shift): (i16, i16) = if k1_q15 > 0 { (16384, 15) } else { (2048, 12) };
        let ga = div_s(dividend.min(d), d);
        let mu_half = shr(mu, 1);

        let mut out = [0i16; L_SUBFR];
        for (n, o) in out.iter_mut().enumerate() {
            let prev = if n == 0 { self.tilt_prev } else { x[n - 1] };
            let acc = l_mac(l_shl(l_deposit_l(x[n]), 15), mu_half, prev);
            let v = extract_l(l_shr(l_add(acc, 0x4000), 15));
            let scaled = l_add(l_mult(v, ga), i32::from(dividend));
            *o = sature(l_shr(scaled, out_shift));
        }
        self.tilt_prev = x[L_SUBFR - 1];
        out
    }

    /// §A.4.2.3 `H_t(z) = 1 + γ_t·k′_1·z⁻¹` (no `1/g_t`), `k′_1` from
    /// the un-normalised high words of the length-22 impulse-response
    /// autocorrelation (eq (A.14)), `γ_t = 0.8` when `r_h(1) > 0`,
    /// else no tilt.
    fn tilt_annex_a(&mut self, x: &[i16; L_SUBFR], h: &[i16]) -> [i16; L_SUBFR] {
        let mut rh0 = 0i32;
        let mut rh1 = 0i32;
        for j in 0..h.len() {
            rh0 = l_mac(rh0, h[j], h[j]);
            if j + 1 < h.len() {
                rh1 = l_mac(rh1, h[j], h[j + 1]);
            }
        }
        let r0 = extract_h(rh0);
        let r1 = extract_h(rh1);
        let coef = if r1 <= 0 || r0 <= 0 {
            0
        } else {
            let num = mult(r1, GAMMA_T_ANNEX_A_Q15);
            div_s(num.min(r0), r0)
        };
        let mut out = [0i16; L_SUBFR];
        for (n, o) in out.iter_mut().enumerate() {
            let prev = if n == 0 { self.tilt_prev } else { x[n - 1] };
            *o = sub(x[n], mult(coef, prev));
        }
        self.tilt_prev = x[L_SUBFR - 1];
        out
    }

    /// eq (88): `G = Σ|ŝ|/Σ|sf|` on Q14 from the two Word32 sums
    /// (`s_in > 0`, `s_out > 0`).
    fn agc_target_q14(s_in: i32, s_out: i32) -> i16 {
        let e_in = norm_l(s_in);
        let m_in = extract_h(l_shl(s_in, e_in));
        let e_out = norm_l(s_out);
        let m_out = extract_h(l_shl(s_out, e_out));
        let mut k = add(sub(e_in, e_out), 1);
        let q = if m_in < m_out {
            div_s(m_in, m_out)
        } else {
            k = sub(k, 1);
            add(shr(div_s(sub(m_in, m_out), m_out), 1), 16384)
        };
        shr(q, k)
    }

    /// §4.2.4 adaptive gain control (eqs (88)–(90)) on the Q14 gain
    /// grid; returns `sf′(n)` and the per-sample gain trajectory.
    fn agc(
        &mut self,
        s_hat: &[i16; L_SUBFR],
        sf: &[i16; L_SUBFR],
    ) -> ([i16; L_SUBFR], [i16; L_SUBFR]) {
        let mut s_in = 0i32;
        let mut s_out = 0i32;
        for n in 0..L_SUBFR {
            s_in = l_add(s_in, l_abs(l_deposit_l(s_hat[n])));
            s_out = l_add(s_out, l_abs(l_deposit_l(sf[n])));
        }
        if s_out == 0 {
            // Silent postfiltered subframe: gain state to 0, no scaling.
            self.agc_gain = 0;
            return (*sf, [0; L_SUBFR]);
        }
        let g_inc = if s_in == 0 {
            0
        } else {
            mult_r(Self::agc_target_q14(s_in, s_out), AGC_TARGET_Q15)
        };
        let mut out = [0i16; L_SUBFR];
        let mut traj = [0i16; L_SUBFR];
        let mut g = self.agc_gain;
        for n in 0..L_SUBFR {
            g = add(mult_r(AGC_PREV_Q15, g), g_inc);
            traj[n] = g;
            out[n] = round(l_shl(l_mult(g, sf[n]), 1));
        }
        self.agc_gain = g;
        (out, traj)
    }

    /// §A.4.2.4 adaptive gain control (eq (A.15)) on the Q12 gain grid.
    fn agc_annex_a(
        &mut self,
        s_hat: &[i16; L_SUBFR],
        sf: &[i16; L_SUBFR],
    ) -> ([i16; L_SUBFR], [i16; L_SUBFR]) {
        let mut e_out = 0i32;
        let mut e_in = 0i32;
        for n in 0..L_SUBFR {
            let o = shr(sf[n], 2);
            e_out = l_mac(e_out, o, o);
            let i = shr(s_hat[n], 2);
            e_in = l_mac(e_in, i, i);
        }
        if e_out == 0 {
            // Zero quarter-scale output energy: gain state to 0, the
            // subframe passes unscaled.
            self.agc_gain = 0;
            return (*sf, [0; L_SUBFR]);
        }
        let g_inc = if e_in == 0 {
            0
        } else {
            let exp_out = sub(norm_l(e_out), 1);
            let m_out = round(l_shl(e_out, exp_out));
            let exp_in = norm_l(e_in);
            let m_in = round(l_shl(e_in, exp_in));
            // E_out/E_in on Q22, its inverse square root on Q19, G on Q12.
            let q = div_s(m_out.min(m_in), m_in);
            let ratio = l_shr(l_shl(l_deposit_l(q), 7), sub(exp_out, exp_in));
            let g = round(l_shl(inv_sqrt(ratio), 9));
            mult(g, AGC_TARGET_ANNEX_A_Q15)
        };
        let mut out = [0i16; L_SUBFR];
        let mut traj = [0i16; L_SUBFR];
        let mut g = self.agc_gain;
        for n in 0..L_SUBFR {
            g = add(mult(g, AGC_PREV_ANNEX_A_Q15), g_inc);
            traj[n] = g;
            out[n] = extract_h(l_shl(l_mult(sf[n], g), 3));
        }
        self.agc_gain = g;
        (out, traj)
    }

    /// §4.2.5 output high-pass (eq (91), staged Q13 coefficients) with
    /// the ×2 upscale folded into the accumulator scale; the feedback
    /// memory keeps the unrounded Word32 recursion value and the
    /// feedback products are carried in full (each landed on Q15 and
    /// realigned).
    fn high_pass(&mut self, x: &[i16; L_SUBFR]) -> [i16; L_SUBFR] {
        let b = &HPF_PREPROC_100HZ_B_Q13;
        let a = &HPF_PREPROC_100HZ_A_Q13;
        let mut out = [0i16; L_SUBFR];
        for (n, o) in out.iter_mut().enumerate() {
            let ff = i64::from(x[n]) * i64::from(b[0])
                + i64::from(self.hp_x[0]) * i64::from(b[1])
                + i64::from(self.hp_x[1]) * i64::from(b[2]);
            let fb1 = ((i64::from(self.hp_y[0]) * i64::from(a[1])) >> 15) << 2;
            let fb2 = ((i64::from(self.hp_y[1]) * i64::from(a[2])) >> 15) << 2;
            let acc = crate::fx::ops::sature32((ff << 2) + fb1 + fb2);
            self.hp_y[1] = self.hp_y[0];
            self.hp_y[0] = acc;
            self.hp_x[1] = self.hp_x[0];
            self.hp_x[0] = x[n];
            *o = round(l_shl(acc, 2));
        }
        out
    }

    /// Run the whole §4.2 cascade on one subframe.
    ///
    /// `s` is the §4.1.6 reconstructed speech `ŝ(n)`; `a_q12` the
    /// subframe's Q12 LP coefficients (`[4096, â_1 … â_10]`); `int_t1`
    /// the integer part of the frame's first-subframe pitch delay
    /// (clause 4.2.1 anchors both subframes on it). Returns the final
    /// §4.2.5 PCM and the §4.2.1 decision (the §4.4 voicing classifier
    /// input).
    pub fn process_subframe(
        &mut self,
        s: &[i16; L_SUBFR],
        a_q12: &[i16; M + 1],
        int_t1: usize,
    ) -> ([i16; L_SUBFR], LtDecisionFx) {
        let t = self.process_subframe_traced(s, a_q12, int_t1);
        (t.output, t.decision)
    }

    /// [`Self::process_subframe`] with the §4.2.1 filter forced off
    /// (search still runs; history still advances) — oracle probes.
    #[doc(hidden)]
    pub fn process_subframe_lt_off(
        &mut self,
        s: &[i16; L_SUBFR],
        a_q12: &[i16; M + 1],
        int_t1: usize,
    ) -> ([i16; L_SUBFR], LtDecisionFx) {
        self.force_lt_off = true;
        let t = self.process_subframe_traced(s, a_q12, int_t1);
        self.force_lt_off = false;
        (t.output, t.decision)
    }

    /// [`Self::process_subframe`] with every intermediate stage signal
    /// returned — the residual-divergence bisection instrument.
    #[doc(hidden)]
    pub fn process_subframe_traced(
        &mut self,
        s: &[i16; L_SUBFR],
        a_q12: &[i16; M + 1],
        int_t1: usize,
    ) -> SubframeTraceFx {
        let apn = weight_az(a_q12, GAMMA_N_Q15);
        let apd = weight_az(a_q12, GAMMA_D_Q15);

        // eq (79): residual through Â(z/γ_n) — feeds BOTH the §4.2.1
        // search/filter and the short-term numerator (clause 4.2).
        self.residual(s, &apn);
        for i in 0..M {
            self.s_hist[i] = s[L_SUBFR - 1 - i];
        }
        let residual: [i16; L_SUBFR] = std::array::from_fn(|n| self.r_buf[HIST + n]);

        // §4.2.1 H_p(z) on the residual.
        let (lt, decision) = self.long_term(int_t1);
        let mut t = self.after_long_term(s, &apn, &apd, &lt, None);
        t.decision = decision;
        t.residual = residual;
        t
    }

    /// Everything downstream of the §4.2.1 output: `1/g_f`, the
    /// all-pole section, tilt compensation, AGC and the output
    /// high-pass (main body), or tilt → all-pole → AGC → high-pass
    /// (Annex A). `all_pole_out` substitutes a given all-pole output
    /// for the computed one (stage replay; the recursion memory is
    /// then left untouched).
    fn after_long_term(
        &mut self,
        s: &[i16; L_SUBFR],
        apn: &[i16; M],
        apd: &[i16; M],
        lt: &[i16; L_SUBFR],
        all_pole_out: Option<&[i16; L_SUBFR]>,
    ) -> SubframeTraceFx {
        let decision = LtDecisionFx::off(0);
        let residual = [0i16; L_SUBFR];
        let mut impulse = [0i16; ANNEX_A_IMPULSE_LEN];
        if self.annex_a {
            // §A.4.2: tilt compensation BEFORE the 1/Â(z/γ_d) synthesis,
            // no g_f / g_t.
            impulse_response_into(apn, apd, &mut impulse);
            let tilted = self.tilt_annex_a(lt, &impulse);
            let sf = match all_pole_out {
                Some(y) => *y,
                None => self.synthesis(&tilted, apd),
            };
            let agc_gain_in = self.agc_gain;
            let (agc, gain) = self.agc_annex_a(s, &sf);
            let output = self.high_pass(&agc);
            return SubframeTraceFx {
                decision,
                residual,
                long_term: *lt,
                impulse,
                gf_q12: 4096,
                gf_scaled: tilted,
                short_term: tilted,
                tilt: sf,
                agc_gain_in,
                agc_gain_out: self.agc_gain,
                gain,
                agc,
                output,
            };
        }

        // eq (85): g_f over the 20-tap impulse response, Q12 sum
        // truncated to Q10; 1/g_f scales the all-pole input only when
        // g_f > 1.
        impulse_response_into(apn, apd, &mut impulse[..GF_IMPULSE_LEN]);
        let h = &impulse[..GF_IMPULSE_LEN];
        let mut gf = 0i32;
        for &hv in h {
            gf = l_add(gf, l_deposit_l(abs_s(hv)));
        }
        let gf_q10 = extract_h(l_shl(gf, 14));
        let gf_scaled: [i16; L_SUBFR] = if gf_q10 > 1024 {
            let inv = div_s(1024, gf_q10);
            std::array::from_fn(|n| mult_r(lt[n], inv))
        } else {
            *lt
        };
        let hf = match all_pole_out {
            Some(y) => *y,
            None => self.synthesis(&gf_scaled, apd),
        };

        // §4.2.3 H_t(z).
        let k1 = Self::tilt_factor(h);
        let sf = self.tilt(&hf, k1);

        // §4.2.4 AGC against the reconstructed speech ŝ(n).
        let agc_gain_in = self.agc_gain;
        let (agc, gain) = self.agc(s, &sf);

        // §4.2.5 output high-pass + ×2.
        let output = self.high_pass(&agc);
        SubframeTraceFx {
            decision,
            residual,
            long_term: *lt,
            impulse,
            gf_q12: gf,
            gf_scaled,
            short_term: hf,
            tilt: sf,
            agc_gain_in,
            agc_gain_out: self.agc_gain,
            gain,
            agc,
            output,
        }
    }

    /// Stage replay: seed the inter-subframe memories from known
    /// values — the last ten `ŝ` samples (most recent last), the
    /// residual history (most recent last, any length), the last ten
    /// all-pole outputs (most recent last), the tilt input tap and the
    /// AGC gain state. The §4.2.5 memory is not seeded.
    #[doc(hidden)]
    pub fn replay_seed(
        &mut self,
        s_prev: &[i16],
        r_hist: &[i16],
        y_prev: &[i16],
        tilt_prev: i16,
        agc_gain: i16,
    ) {
        for (i, v) in s_prev.iter().rev().take(M).enumerate() {
            self.s_hist[i] = *v;
        }
        let n = r_hist.len().min(HIST);
        self.r_buf[..HIST].fill(0);
        self.r_buf[HIST - n..HIST].copy_from_slice(&r_hist[r_hist.len() - n..]);
        for (i, v) in y_prev.iter().rev().take(M).enumerate() {
            self.st_y[i] = *v;
        }
        self.tilt_prev = tilt_prev;
        self.agc_gain = agc_gain;
    }

    /// Stage replay: run the stages downstream of a **given** §4.2.1
    /// output (and optionally a given all-pole output) for one
    /// subframe. The returned trace carries no decision / residual.
    #[doc(hidden)]
    pub fn replay_after_long_term(
        &mut self,
        s: &[i16; L_SUBFR],
        a_q12: &[i16; M + 1],
        lt: &[i16; L_SUBFR],
        all_pole_out: Option<&[i16; L_SUBFR]>,
    ) -> SubframeTraceFx {
        let apn = weight_az(a_q12, GAMMA_N_Q15);
        let apd = weight_az(a_q12, GAMMA_D_Q15);
        self.after_long_term(s, &apn, &apd, lt, all_pole_out)
    }
}

/// Every intermediate §4.2 stage signal for one subframe (the
/// bisection instrument behind
/// [`PostfilterFx::process_subframe_traced`]).
#[doc(hidden)]
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct SubframeTraceFx {
    /// §4.2.1 decision.
    pub decision: LtDecisionFx,
    /// eq (79) residual `r̂(n)`.
    pub residual: [i16; L_SUBFR],
    /// §4.2.1 output.
    pub long_term: [i16; L_SUBFR],
    /// Impulse response `h_f(n)` on Q12 (20 taps main body, 22 Annex A).
    pub impulse: [i16; ANNEX_A_IMPULSE_LEN],
    /// eq (85) `Σ|h_f|` on Q12.
    pub gf_q12: i32,
    /// All-pole input: the §4.2.1 output after `1/g_f` (main body) or
    /// after the §A.4.2.3 tilt filter (Annex A).
    pub gf_scaled: [i16; L_SUBFR],
    /// `H_f(z)` output (main body); the tilt output for Annex A.
    pub short_term: [i16; L_SUBFR],
    /// AGC input `sf(n)`.
    pub tilt: [i16; L_SUBFR],
    /// `g(−1)` entering the subframe (Q14 main body, Q12 Annex A).
    pub agc_gain_in: i16,
    /// `g(39)` leaving the subframe.
    pub agc_gain_out: i16,
    /// Per-sample gain trajectory `g(n)`.
    pub gain: [i16; L_SUBFR],
    /// eq (89) output `sf′(n)`.
    pub agc: [i16; L_SUBFR],
    /// §4.2.5 output PCM.
    pub output: [i16; L_SUBFR],
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A flat LP set (all-zero â) makes every stage a near-identity:
    /// γ-weighting is zero, g_f = 1, k1' = 0, and the AGC target is 1,
    /// so the cascade output tracks the ×2-upscaled high-passed input
    /// within the sub-LSB normalisation error of each stage.
    #[test]
    fn flat_lp_reduces_to_output_stage() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        let mut pf = PostfilterFx::new();
        let s: [i16; L_SUBFR] = std::array::from_fn(|n| (n as i16 - 20) * 100);
        let (out, d) = pf.process_subframe(&s, &a, 40);

        // Independent eq (91) reference with unrounded f64 feedback.
        let bq = &HPF_PREPROC_100HZ_B_Q13;
        let aq = &HPF_PREPROC_100HZ_A_Q13;
        let (mut x1, mut x2, mut y1, mut y2) = (0f64, 0f64, 0f64, 0f64);
        for n in 0..L_SUBFR {
            let x = f64::from(s[n]);
            let y = (f64::from(bq[0]) * x
                + f64::from(bq[1]) * x1
                + f64::from(bq[2]) * x2
                + f64::from(aq[1]) * y1
                + f64::from(aq[2]) * y2)
                / 8192.0;
            x2 = x1;
            x1 = x;
            y2 = y1;
            y1 = y;
            let expect = 2.0 * y;
            assert!(
                (f64::from(out[n]) - expect).abs() <= 4.0,
                "n={n}: {} vs {expect}",
                out[n]
            );
        }
        assert!(d.gain_q15 >= 0);
    }

    /// The weighted coefficients are γ^i·â_i within one Q12 LSB of the
    /// float weighting.
    #[test]
    fn weight_az_tracks_float() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2458; // −0.6
        a[2] = 1229; // 0.3
        a[5] = -410; // −0.1
        let w = weight_az(&a, GAMMA_N_Q15);
        let mut g = 0.55f64;
        for i in 0..M {
            let expect = g * f64::from(a[i + 1]);
            assert!(
                (f64::from(w[i]) - expect).abs() <= 2.0,
                "i={i}: {} vs {expect}",
                w[i]
            );
            g *= 0.55;
        }
    }

    /// The Q12 impulse response matches the float short-term impulse
    /// response within a couple of LSBs.
    #[test]
    fn impulse_response_tracks_float() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2867; // −0.7
        a[2] = 1638; // 0.4
        let apn = weight_az(&a, GAMMA_N_Q15);
        let apd = weight_az(&a, GAMMA_D_Q15);
        let h = impulse_response_q12(&apn, &apd);

        let af: [f32; M] = std::array::from_fn(|i| f32::from(a[i + 1]) / 4096.0);
        let hf = crate::short_term_postfilter::ShortTermPostfilter::impulse_response(&af);
        for n in 0..GF_IMPULSE_LEN {
            let expect = f64::from(hf[n]) * 4096.0;
            assert!(
                (f64::from(h[n]) - expect).abs() <= 4.0,
                "n={n}: {} vs {expect}",
                h[n]
            );
        }
    }

    /// The γ-power ladders are the recursively rounded ones the staged
    /// traces print (not exact powers: exact 0.55⁵·2¹⁵ = 1650.6).
    #[test]
    fn gamma_ladders_match_staged_schedule() {
        let mut a = [4096i16; M + 1];
        a[0] = 4096;
        // â_i = 1.0 (Q12) makes round(L_mult(â_i, γ^i)) = round(γ^i / 8).
        let ladder = |g: i16| -> [i16; M] {
            let mut f = g;
            std::array::from_fn(|j| {
                let v = f;
                if j + 1 < M {
                    f = round(l_mult(f, g));
                }
                v
            })
        };
        assert_eq!(
            ladder(GAMMA_N_Q15),
            [18022, 9912, 5451, 2998, 1649, 907, 499, 274, 151, 83]
        );
        assert_eq!(
            ladder(GAMMA_D_Q15),
            [22938, 16057, 11240, 7868, 5508, 3856, 2699, 1889, 1322, 925]
        );
        let w = weight_az(&a, GAMMA_N_Q15);
        assert_eq!(w[0], round(l_mult(4096, 18022)));
    }

    /// eq (88) worked example from the staged traces: Σ|ŝ| = 106288,
    /// Σ|sf| = 82470 → G = 21116 (Q14); and 113123 / 100407 → 18459.
    #[test]
    fn agc_target_worked_examples() {
        assert_eq!(PostfilterFx::agc_target_q14(106_288, 82_470), 21116);
        assert_eq!(PostfilterFx::agc_target_q14(113_123, 100_407), 18459);
        assert_eq!(PostfilterFx::agc_target_q14(2, 2), 16384);
    }

    /// Unit gain is an exact fixed point of the eq (90) recursion on
    /// the Q14 grid with the 32358/410 pair.
    #[test]
    fn agc_unit_gain_is_a_fixed_point() {
        let inc = mult_r(AGC_UNITY_Q14, AGC_TARGET_Q15);
        assert_eq!(add(mult_r(AGC_PREV_Q15, AGC_UNITY_Q14), inc), AGC_UNITY_Q14);
    }

    /// Silence propagates: an all-zero subframe with any stable LP set
    /// produces all-zero output; the silent postfiltered subframe
    /// resets the gain state to 0 (§4.2.4 degenerate case).
    #[test]
    fn silence_stays_silent() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2048;
        let mut pf = PostfilterFx::new();
        let (out, _) = pf.process_subframe(&[0; L_SUBFR], &a, 50);
        assert_eq!(out, [0i16; L_SUBFR]);
        assert_eq!(pf.agc_gain, 0);
    }

    /// Hostile input: full-scale alternating samples with an unstable
    /// coefficient set must not panic in either cascade.
    #[test]
    fn hostile_input_is_panic_free() {
        for annex_a in [false, true] {
            let mut pf = if annex_a {
                PostfilterFx::new_annex_a()
            } else {
                PostfilterFx::new()
            };
            let mut a = [i16::MIN; M + 1];
            a[0] = 4096;
            for k in 0..64usize {
                let s: [i16; L_SUBFR] = std::array::from_fn(|n| match (n + k) % 4 {
                    0 => i16::MAX,
                    1 => i16::MIN,
                    2 => (k as i16).wrapping_mul(977),
                    _ => 0,
                });
                if k % 3 == 0 {
                    a = std::array::from_fn(|i| {
                        if i == 0 {
                            4096
                        } else {
                            ((i * 7919 + k * 104_729) as i16).wrapping_mul(31)
                        }
                    });
                }
                let _ = pf.process_subframe(&s, &a, 1 + (k * 37) % 400);
            }
        }
    }

    /// Determinism: two cascades fed the same stream stay in lockstep.
    #[test]
    fn deterministic() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2458;
        a[3] = 819;
        let s: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 37) % 200) as i16 - 100);
        let mut p = PostfilterFx::new();
        let mut q = PostfilterFx::new();
        for _ in 0..6 {
            assert_eq!(
                p.process_subframe(&s, &a, 45),
                q.process_subframe(&s, &a, 45)
            );
        }
    }

    /// The §4.2.1 search stays inside the clamped window and the gain
    /// stays on the bounded Q15 grid for arbitrary periodic input.
    #[test]
    fn lt_decision_bounded() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        let mut pf = PostfilterFx::new();
        for k in 0..30usize {
            let s: [i16; L_SUBFR] = std::array::from_fn(|n| (((n + 7 * k) % 40) as i16 - 20) * 300);
            let t1 = 30 + k;
            let (_out, d) = pf.process_subframe(&s, &a, t1);
            // The applied delay T = delay + frac/8 sits in the
            // two-sided window [int(T1) − 1 − 7/8, int(T1) + 1 + 7/8].
            let t8 = 8 * d.delay + d.frac;
            assert!(
                t8 > 8 * (t1 - 2) && t8 <= 8 * (t1 + 1) + 7,
                "T = {t8}/8 for t1 {t1}"
            );
            assert!(d.frac < POSTFILTER_FRAC_RES);
            assert!(d.gain_q15 >= 0);
            assert!(d.undelayed_q15 >= LT_UNDELAYED_MIN_Q15);
        }
    }
}
