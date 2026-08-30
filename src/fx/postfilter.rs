//! §4.2 fixed-point post-processing cascade: the long-term (harmonic)
//! postfilter (§4.2.1, eqs (78)–(83)), the short-term postfilter
//! (§4.2.2, eqs (84)/(85)), tilt compensation (§4.2.3, eqs (86)/(87)),
//! adaptive gain control (§4.2.4, eqs (88)–(90)), and the output
//! high-pass + ×2 upscaling (§4.2.5, eq (91)) — everything on the
//! clause-5 Word16/Word32 operator grid, consuming the §4.1
//! [`crate::fx::decoder::FrameDecoderFx`] output.
//!
//! ## Number grids
//!
//! - **Speech / residual / stage signals** — Q0 (the halved-input PCM
//!   grid the §4.1.6 synthesis produces; §4.2.5 restores ×2).
//! - **Weighted LP coefficients** `γ^i·â_i` — Q12, the γ powers
//!   accumulated on Q15 ([`weight_az`]).
//! - **eq (83) gain `g_l`** — Q15, bounded `[0, 32767]`; the eq (78)
//!   combination runs `1/(1 + γ_p·g_l)` as an exact `div_s` on a Q14
//!   denominator.
//! - **eq (85) `g_f` / eq (86) `g_t`** — Word32 (Q12 / Q15), inverted
//!   through the shared normalised-reciprocal helper ([`recip_norm`]).
//! - **§4.2.4 gain `g(n)`** — Q12 (unit gain = 4096); the eq (90)
//!   weights are the Q15 pair `27853/4915` (sum exactly `2^15`).
//! - **§4.2.5** — the staged Q13 coefficient tables against a Word32
//!   Q14 feedback memory (the recursion value is kept unrounded, as
//!   the §3.1 twin of this filter measurably requires).
//!
//! ## Latitude
//!
//! The Recommendation pins every *quantity* here (the equations, the
//! constants, the cascade order, the Table 9 init) but not the
//! reference's internal fixed-point operator schedule for §4.2 (a
//! recorded docs-gap). The §4.2.1 decision correlations run on a
//! norm-scaled Word16 copy of the residual window (the Word32
//! accumulators stay in range by construction — the clause-5 grid has
//! no wider accumulator), ratio *comparisons* are exact
//! cross-multiplications of those Word32 values, and the stated
//! output grids (`g_l` Q15, `1/g_f`/`1/g_t` via [`recip_norm`]) land
//! by normalised division. The remaining sub-LSB schedule choices are
//! exposed as [`PfLatitudeFx`] and arbitrated black-box against the
//! conformance corpus (see `tests/fx_full_conformance.rs`).

use crate::fx::dsp::{l_extract, mpy_32_16};
use crate::fx::ops::{
    abs_s, add, div_s, extract_h, l_add, l_deposit_h, l_deposit_l, l_mac, l_mult, l_shl, l_sub,
    mult, norm_l, round, shl, shr, sub,
};
use crate::tables::{
    postfilter_interp_long, postfilter_interp_short, HPF_PREPROC_100HZ_A_Q13,
    HPF_PREPROC_100HZ_B_Q13, M, POSTFILTER_FRAC_RES, POSTFILTER_INTERP_LONG_TAPS,
    POSTFILTER_INTERP_SHORT_TAPS,
};

/// Samples per subframe.
pub const L_SUBFR: usize = 40;

/// §4.2.1 weight factor `γ_p = 0.5` — folded structurally into the
/// eq (78) arithmetic (`g_half = g_l/2`, exact shifts).
pub const GAMMA_P_SHIFT: i16 = 1;

/// §4.2.2 numerator weight `γ_n = 0.55` on Q15.
pub const GAMMA_N_Q15: i16 = 18022;

/// §4.2.2 denominator weight `γ_d = 0.70` on Q15.
pub const GAMMA_D_Q15: i16 = 22938;

/// §4.2.3 tilt weight `γ_t = 0.9` (negative `k1'`) on Q15.
pub const GAMMA_T_NEG_Q15: i16 = 29491;

/// §4.2.3 tilt weight `γ_t = 0.2` (positive `k1'`) on Q15.
pub const GAMMA_T_POS_Q15: i16 = 6554;

/// eq (90) previous-gain weight `0.85` on Q15.
pub const AGC_PREV_Q15: i16 = 27853;

/// eq (90) target weight `0.15` on Q15 (`27853 + 4915 = 2^15`).
pub const AGC_TARGET_Q15: i16 = 4915;

/// eq (85) impulse-response truncation length.
pub const GF_IMPULSE_LEN: usize = 20;

/// Deepest §4.2.1 history reach: the longest integer delay the search
/// can pick (`int(T_1) + 1 ≤ 144`) plus the long-interpolator margin
/// (8 samples beyond the integer anchor).
const HIST: usize = 144 + POSTFILTER_INTERP_LONG_TAPS / 2;

/// One subframe's §4.2.1 decision (mirrors the float
/// [`crate::long_term_postfilter::LongTermDecision`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LtDecisionFx {
    /// Chosen integer delay `T_0` (eq (80)).
    pub delay: usize,
    /// Chosen 1/8-resolution fractional offset (`0 … 7`).
    pub frac: usize,
    /// eq (83) gain `g_l` on Q15 (`0` = disabled, eq (82)).
    pub gain_q15: i16,
    /// Whether the length-129 refinement replaced the length-33
    /// delayed signal (clause 4.2.1 keeps the longer-filter signal
    /// only when it raises `R′(T)`); the eq (78) filter reads the
    /// delayed residual through the same kernel.
    pub use_long: bool,
}

/// Normalised reciprocal of a positive Word32: returns `(m, e)` such
/// that `1/L_x = m · 2^(e−45)` with `m = div_s(2^14, hi)` over the
/// `norm_l`-normalised high word (`m ∈ (2^14, 2^15]`, saturated).
fn recip_norm(l_x: i32) -> (i16, i16) {
    debug_assert!(l_x > 0);
    let e = norm_l(l_x);
    let hi = extract_h(l_shl(l_x, e));
    (div_s(16384, hi), e)
}

/// Unpinned §4.2 operator-schedule latitude (the Recommendation pins
/// the equations, not the reference's fixed-point rounding points) —
/// pinned black-box against the conformance corpus.
#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PfLatitudeFx {
    /// eq (79) residual lands on Q0 by rounding (`true`) or
    /// truncation (`false`).
    pub resid_round: bool,
    /// Synthesis output `y(n)` lands on Q0 by rounding or truncation.
    pub syn_round: bool,
    /// eq (78) output lands by rounding or truncation.
    pub lt_round: bool,
    /// `1/g_f` applied to the all-pole output (`false`) or to its
    /// input (`true`).
    pub gf_before: bool,
    /// §4.2.5 feedback kept on the wide Word32 grid (`true`) or
    /// rounded to Word16 (`false`).
    pub hp_wide: bool,
    /// §4.2.4 gain grid shift relative to Q12 (0 = Q12, 1 = Q13, …).
    pub agc_shift: i16,
    /// The §4.2.5 ×2 upscale folded into the eq (89) product landing:
    /// the AGC output lands on Q1 (`2·sf′(n)` rounded once) and the
    /// high-pass consumes that grid, landing its output on Q0 without
    /// a second scaling.
    pub agc_x2: bool,
    /// eq (88) target as the square root of the amplitude ratio
    /// (black-box sweep hook; r452 measured: pushes SPEECH's first
    /// divergence 83 → 271 but collapses LSP/PITCH exact share — off).
    pub agc_sqrt: bool,
    /// eq (89)/(90): scale with `g(n−1)` first, then advance the
    /// recursion (black-box sweep hook).
    pub agc_apply_first: bool,
    /// eq (90) smoothing pole on Q15 (`0` = the printed 0.85 pair).
    pub agc_pole_q15: i16,
    /// eq (88) ratio taken from the PREVIOUS subframe's `Σ|ŝ|/Σ|sf|`
    /// sums (one-subframe lag) instead of the current one.
    ///
    /// r452 measured: this reproduces the reference's AGC trajectory at
    /// silence→signal onsets exactly (PITCH/FIXED/TAME frame-0
    /// subframe-1 heads become byte-exact, and every clean vector's
    /// first divergence moves later: FIXED 43 → 47, PITCH 41 → 46,
    /// TAME 41 → 53, SPEECH 83 → 141), but costs aggregate exact share
    /// (sum over 12 vectors 112.7 → 103.7) while the §4.2.2/§4.2.3
    /// sub-LSB schedule is still unpinned — the residual ±1-LSB
    /// `st`/`tilt` differences then land on the wrong side of the
    /// rounding boundary more often. Left OFF until that schedule is
    /// pinned; the onset evidence says it is the reference's structure.
    pub agc_lag: bool,
    /// eq (85) `1/g_f` as an exact rounded division by the Word32 gain
    /// instead of the `recip_norm` mantissa multiply (r452 measured:
    /// clearly refuted — exact-share sum 112.7 → 92.8; the mantissa
    /// multiply IS the reference grid).
    pub gf_exact_div: bool,
    /// eq (86) `1/g_t` as an exact rounded division (r452: neutral,
    /// ±0.2 on the exact-share sum — the tilt numbers rarely reach the
    /// differing LSB).
    pub tilt_exact_div: bool,
    /// §4.2.1 fractional-interpolation window origin shift in samples
    /// (the staging note records two candidate tap alignments; r452
    /// measured: ±1 is corpus-neutral, ±0.6 on the exact-share sum
    /// with no first-divergence movement — the current alignment
    /// stands and the residual divergence is not here).
    pub lt_tap_shift: i16,
    /// eq (85)/(87) truncated impulse-response length (0 = the printed
    /// 20; sweep hook, max 40; r452 measured: 22/32/40 all neutral).
    pub gf_len: u8,
}

impl Default for PfLatitudeFx {
    fn default() -> Self {
        Self {
            resid_round: true,
            // Round 452 re-pin (with the corpus-identified frame-0
            // previous-LSP memory in place the r419 truncation pin no
            // longer holds): every landing rounds, the §4.2.5 ×2 is
            // folded into the eq (89) product (Q1 AGC output), and
            // 1/g_f scales the synthesis input. Measured on the full
            // fixed-point chain, sum of exact% over the 12 clean
            // vectors: 79 (r419 defaults) → 109; frame-0 subframe 0
            // byte-exact on all six base vectors.
            syn_round: true,
            lt_round: true,
            gf_before: true,
            hp_wide: true,
            agc_shift: 0,
            agc_x2: true,
            agc_sqrt: false,
            agc_apply_first: false,
            agc_pole_q15: 0,
            agc_lag: false,
            gf_exact_div: false,
            tilt_exact_div: false,
            lt_tap_shift: 0,
            gf_len: 0,
        }
    }
}

/// Land a Q13-scaled accumulator on Q0 by rounding or truncation.
#[inline]
fn land_q13(acc: i32, round_mode: bool) -> i16 {
    let q = l_shl(acc, 3);
    if round_mode {
        round(q)
    } else {
        extract_h(q)
    }
}

/// `γ`-weighted LP coefficients `γ^i·â_i` on Q12 (`i = 1 … 10`,
/// slot `j` holds `i = j + 1`), the γ powers accumulated on Q15.
#[must_use]
pub fn weight_az(a: &[i16; M + 1], gamma_q15: i16) -> [i16; M] {
    let mut fac = gamma_q15;
    std::array::from_fn(|j| {
        let w = mult(a[j + 1], fac);
        if j + 1 < M {
            fac = mult(fac, gamma_q15);
        }
        w
    })
}

/// The truncated impulse response `h_f(n)` (Q12) of the un-normalised
/// short-term postfilter `Â(z/γ_n)/Â(z/γ_d)` (eq (84) without `1/g_f`).
#[must_use]
pub fn impulse_response_q12(apn: &[i16; M], apd: &[i16; M]) -> [i16; GF_IMPULSE_LEN] {
    let h40 = impulse_response_q12_n(apn, apd);
    core::array::from_fn(|n| h40[n])
}

/// [`impulse_response_q12`] extended to 40 taps (the `gf_len` sweep
/// hook consumes a prefix).
#[doc(hidden)]
#[must_use]
pub fn impulse_response_q12_n(apn: &[i16; M], apd: &[i16; M]) -> [i16; 40] {
    // Numerator impulse response: [1, γ_n·â_1 … γ_n^10·â_10] on Q12.
    let mut hnum = [0i16; 40];
    hnum[0] = 4096;
    hnum[1..=M].copy_from_slice(apn);
    // All-pole 1/Â(z/γ_d) from zero state, Q12 in / Q12 out.
    let mut h = [0i16; 40];
    for n in 0..40 {
        let mut acc = l_mult(hnum[n], 4096); // Q25
        for (i, &c) in apd.iter().enumerate() {
            if n > i {
                acc = crate::fx::ops::l_msu(acc, c, h[n - 1 - i]);
            }
        }
        h[n] = round(l_shl(acc, 3)); // Q25 → Q28 → high word Q12
    }
    h
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
    /// §4.2.3 FIR input tap `t(n−1)`.
    tilt_prev: i16,
    /// Previous subframe's eq (88) `(Σ|ŝ|, Σ|sf|)` pair (agc_lag hook).
    agc_prev_sums: (i32, i32),
    /// §4.2.4 smoothed gain `g(n)` on Q12 (Table 9 init 1.0 = 4096).
    agc_q12: i16,
    /// §4.2.5 input taps `[x(n−1), x(n−2)]`.
    hp_x: [i16; 2],
    /// §4.2.5 unrounded feedback `[y(n−1), y(n−2)]` on Word32 Q14.
    hp_y: [i32; 2],
    /// Operator-schedule latitude (black-box-pinned).
    lat: PfLatitudeFx,
    /// Oracle-probe switch: apply the §4.2.1 filter as disabled.
    force_lt_off: bool,
}

impl Default for PostfilterFx {
    fn default() -> Self {
        Self::new()
    }
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
            agc_prev_sums: (0, 0),
            agc_q12: 4096,
            hp_x: [0; 2],
            hp_y: [0; 2],
            lat: PfLatitudeFx::default(),
            force_lt_off: false,
        }
    }

    /// Override the operator-schedule latitude (black-box sweep hook).
    #[doc(hidden)]
    pub fn set_latitude(&mut self, lat: PfLatitudeFx) {
        self.lat = lat;
        self.agc_q12 = shl(4096, lat.agc_shift);
    }

    /// Residual fetch at `r̂(n − k)` (`n` = in-subframe index,
    /// `k ≥ 0`): buffer index `HIST + n − k`.
    #[inline]
    fn at(buf: &[i16; HIST + L_SUBFR], n: usize, k: usize) -> i16 {
        let idx = HIST + n;
        if k > idx {
            0
        } else {
            buf[idx - k]
        }
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
            self.r_buf[HIST + n] = land_q13(acc, self.lat.resid_round);
        }
    }

    /// Interpolated fetch of `buf` at delay `T_0 + frac/8` for output
    /// index `n`, with the short (`long == false`) or long kernel.
    /// `frac == 0` is the direct integer fetch.
    fn interp_at_shift(
        buf: &[i16; HIST + L_SUBFR],
        n: usize,
        t0: usize,
        frac: usize,
        long: bool,
        shift: i16,
    ) -> i16 {
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
            let idx = i64::from(t0 as u32) + i64::from(half as u32) + i64::from(shift) - j as i64;
            let k = usize::try_from(idx).unwrap_or(0);
            acc = l_mac(acc, h, Self::at(buf, n, k));
        }
        round(acc)
    }

    /// Correlation `Σ r̂(n)·r̂_T(n)` and energy `Σ r̂_T(n)²` for a
    /// candidate delay over the **scaled** residual window (the
    /// Word32 accumulators stay in range by construction).
    fn corr_energy_shift(
        rs: &[i16; HIST + L_SUBFR],
        t0: usize,
        frac: usize,
        long: bool,
        shift: i16,
    ) -> (i64, i64) {
        let mut corr = 0i32;
        let mut energy = 0i32;
        for n in 0..L_SUBFR {
            let rk = Self::interp_at_shift(rs, n, t0, frac, long, shift);
            let rn = rs[HIST + n];
            corr = l_mac(corr, rn, rk);
            energy = l_mac(energy, rk, rk);
        }
        (i64::from(corr), i64::from(energy))
    }

    /// The scaled copy of the residual window the §4.2.1 decisions run
    /// on: every sample right-shifted so the worst-case 40-term
    /// correlation of any window slice stays inside Word32 (the
    /// clause-5 grid has no wider accumulator; small residuals lose
    /// their low bits exactly as a Word16 pre-scale loses them).
    fn scaled_residual(&self) -> [i16; HIST + L_SUBFR] {
        let mut maxabs: i16 = 0;
        for &v in self.r_buf.iter() {
            let a = abs_s(v);
            if a > maxabs {
                maxabs = a;
            }
        }
        // Σ40·v² ≤ 2^31 needs |v| ≤ 5792 ≈ 2^12.5.
        let mut sh: i16 = 0;
        while shr(maxabs, sh) > 5792 {
            sh += 1;
        }
        std::array::from_fn(|i| shr(self.r_buf[i], sh))
    }

    /// Clause-4.2.1 two-pass delay search + eq (82)/(83) gain decision
    /// for the residual currently in `r_buf`.
    fn decide(&self, int_t1: usize) -> LtDecisionFx {
        let rs = self.scaled_residual();

        // eq (82) denominator: Σ r̂(n)² on the scaled grid.
        let mut energy32 = 0i32;
        for n in 0..L_SUBFR {
            energy32 = l_mac(energy32, rs[HIST + n], rs[HIST + n]);
        }
        let energy = i64::from(energy32);

        // First pass (eq (80)): best integer delay in
        // [int(T_1) − 1, int(T_1) + 1], clamped to the buffer reach.
        let lo = int_t1.saturating_sub(1).max(1);
        let hi = (int_t1 + 1).min(144);
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

        // Second pass (eq (81)): best 1/8 fraction by R(T)²/E_T (short
        // kernel), positive-correlation candidates only. "Around T_0"
        // spans both sides: `T_0 − 7/8 … T_0 + 7/8`, the negative side
        // expressed as `(T_0 − 1) + frac/8`.
        let mut best = ((t0, 0usize), 0i64, 0i64); // ((anchor, frac), num, den)
        let mut have = false;
        let neg_anchor = t0.saturating_sub(1).max(1);
        let candidates = (1..POSTFILTER_FRAC_RES)
            .map(|frac| (neg_anchor, frac))
            .chain((0..POSTFILTER_FRAC_RES).map(|frac| (t0, frac)));
        for (anchor, frac) in candidates {
            let (corr, e_t) =
                Self::corr_energy_shift(&rs, anchor, frac, false, self.lat.lt_tap_shift);
            if corr <= 0 || e_t <= 0 {
                continue;
            }
            let better = if !have {
                true
            } else {
                // corr²/e_t > best²/best_e  ⇔  corr²·best_e > best²·e_t
                let l = i128::from(corr) * i128::from(corr) * i128::from(best.2);
                let r = i128::from(best.1) * i128::from(best.1) * i128::from(e_t);
                l > r
            };
            if better {
                best = ((anchor, frac), corr, e_t);
                have = true;
            }
        }

        // Long-kernel refinement of a non-integer winner, kept only if
        // it raises R(T)²/E_T.
        let mut use_long = false;
        if have && best.0 .1 != 0 {
            let (corr_l, e_l) =
                Self::corr_energy_shift(&rs, best.0 .0, best.0 .1, true, self.lat.lt_tap_shift);
            if corr_l > 0 && e_l > 0 {
                let l = i128::from(corr_l) * i128::from(corr_l) * i128::from(best.2);
                let r = i128::from(best.1) * i128::from(best.1) * i128::from(e_l);
                if l > r {
                    best = (best.0, corr_l, e_l);
                    use_long = true;
                }
            }
        }

        // eq (82) disable test: R′(T)²/Σr̂² < 0.5 ⇔ 2·num² < energy·den.
        let ((anchor, frac), num, den) = best;
        // Quantisation-noise enable floor, pinned black-box: a
        // subframe whose residual mean square is at most one Q0 LSB
        // (Σr̂² ≤ 40) is treated as having no long-term structure —
        // the eq (82) statistic over bare quantisation noise otherwise
        // flickers the filter on during silence (removing it retires
        // the whole silence-cluster class of divergence events on
        // SPEECH; the threshold is insensitive from 40 to 2560).
        let enabled = have
            && energy > 40
            && 2 * i128::from(num) * i128::from(num) >= i128::from(energy) * i128::from(den);

        if !enabled {
            // Disabled → report the integer-pass anchor, no fraction.
            return LtDecisionFx {
                delay: t0,
                frac: 0,
                gain_q15: 0,
                use_long: false,
            };
        }

        // Over-unity guard, pinned black-box: when the raw eq (83)
        // ratio exceeds 2 (num > 2·den — the ceiling of a Q14 gain
        // grid, i.e. the point where a normalised fixed-point division
        // for g_l stops being representable) the filter behaves as
        // DISABLED, not clamped. The corpus separates the two cleanly:
        // clamping keeps voiced material (SPEECH/PITCH unchanged) but
        // wrecks the onset-heavy FIXED vectors (corr 0.9502/0.9756
        // clamped vs 0.9855/0.9918 disabled); disabling everything
        // over unity instead costs SPEECH/PITCH (0.9953/0.9926).
        if num > 2 * den {
            return LtDecisionFx {
                delay: t0,
                frac: 0,
                gain_q15: 0,
                use_long: false,
            };
        }
        // eq (83): g_l = num/den bounded [0, 1], on Q15.
        let gain_q15 = if num >= den {
            32767
        } else {
            ((num << 15) / den) as i16
        };
        LtDecisionFx {
            delay: anchor,
            frac,
            gain_q15,
            use_long,
        }
    }

    /// §4.2.1 `H_p(z)` for one subframe, applied to the **residual**
    /// (clause 4.2: "the signal r̂ is then filtered through the
    /// long-term postfilter"): two-pass search + the eq (78)
    /// combination over `r̂(n)` / `r̂_T(n)`; advances the residual
    /// history. The residual must already sit in `r_buf[HIST..]`.
    fn long_term(&mut self, int_t1: usize) -> ([i16; L_SUBFR], LtDecisionFx) {
        let mut d = self.decide(int_t1);
        if self.force_lt_off {
            d.gain_q15 = 0;
            d.frac = 0;
            d.use_long = false;
        }

        let mut out = [0i16; L_SUBFR];
        if d.gain_q15 == 0 {
            out.copy_from_slice(&self.r_buf[HIST..]);
        } else {
            // γ_p·g_l on Q15 (γ_p = 0.5 exact shift); the eq (78)
            // denominator 1 + γ_p·g_l on Q14; 1/(1 + γ_p·g_l) by exact
            // Q15 division.
            let g_half = shr(d.gain_q15, GAMMA_P_SHIFT);
            let denom_q14 = add(16384, shr(d.gain_q15, GAMMA_P_SHIFT + 1));
            let inv_q15 = div_s(16384, denom_q14);
            for (n, o) in out.iter_mut().enumerate() {
                // r̂(n) + γ_p·g_l·r̂_T(n) on Q16 — the delayed residual
                // through the kernel the search settled on (clause
                // 4.2.1 longer-filter replacement rule).
                let r_t = Self::interp_at_shift(
                    &self.r_buf,
                    n,
                    d.delay,
                    d.frac,
                    d.use_long,
                    self.lat.lt_tap_shift,
                );
                let acc = l_mac(l_deposit_h(self.r_buf[HIST + n]), g_half, r_t);
                let (hi, lo) = l_extract(acc);
                let scaled = l_shl(mpy_32_16(hi, lo, inv_q15), 1);
                *o = if self.lat.lt_round {
                    round(scaled)
                } else {
                    extract_h(scaled)
                };
            }
        }

        // Advance the residual history by one subframe.
        self.r_buf.copy_within(L_SUBFR.., 0);
        (out, d)
    }

    /// The `1/[g_f·Â(z/γ_d)]` synthesis stage: all-pole filter over the
    /// long-term-postfiltered residual, then the eq (85) `1/g_f`
    /// normalisation; continuous output memory.
    fn synthesis(
        &mut self,
        z: &[i16; L_SUBFR],
        apd: &[i16; M],
        gf_recip: (i16, i16),
        gf32_q12: i32,
    ) -> [i16; L_SUBFR] {
        let (gf_m, gf_e) = gf_recip;
        // Exact division by the Q12 gain: round(x·4096/g_f).
        let exact_div = |x: i16| -> i16 {
            if gf32_q12 <= 0 {
                return x;
            }
            let num = i64::from(x) << 13; // ×2·4096 for the rounding add
            let q = (num / i64::from(gf32_q12) + 1) >> 1;
            q.clamp(-32768, 32767) as i16
        };
        let mut y_un = [0i16; L_SUBFR]; // un-normalised y(n)
        let mut out = [0i16; L_SUBFR];
        for n in 0..L_SUBFR {
            // Optionally scale the input by 1/g_f first (the eq (84)
            // leading factor commutes; the rounding point differs).
            let zin = if self.lat.gf_before {
                if self.lat.gf_exact_div {
                    exact_div(z[n])
                } else {
                    round(l_shl(l_mult(z[n], gf_m), sub(gf_e, 18)))
                }
            } else {
                z[n]
            };
            // y(n) = z(n) − Σ apd[i]·y(n−1−i), Q13 accumulator.
            let mut acc = l_mult(zin, 4096);
            for (i, &cf) in apd.iter().enumerate() {
                let yv = if n > i {
                    y_un[n - 1 - i]
                } else {
                    self.st_y[i - n]
                };
                acc = crate::fx::ops::l_msu(acc, cf, yv);
            }
            let y = land_q13(acc, self.lat.syn_round);
            y_un[n] = y;

            // 1/g_f: y·m·2^(e−33) with l_mult’s doubling folded in
            // (Q16 landing → round to Q0).
            out[n] = if self.lat.gf_before {
                y
            } else if self.lat.gf_exact_div {
                exact_div(y)
            } else {
                round(l_shl(l_mult(y, gf_m), sub(gf_e, 18)))
            };
        }
        // Advance the continuous memory (most-recent at index 0).
        for i in 0..M {
            self.st_y[i] = y_un[L_SUBFR - 1 - i];
        }
        out
    }

    /// §4.2.3 `H_t(z)`: `sf(n) = (1/g_t)·(t(n) + γ_t·k1'·t(n−1))` with
    /// the eq (87) tilt factor from the Q12 impulse response.
    fn tilt(&mut self, t: &[i16; L_SUBFR], h: &[i16]) -> [i16; L_SUBFR] {
        // eq (87): k1' = −r_h(1)/r_h(0), wide-exact on the Q12 grid.
        let mut rh0 = 0i64;
        let mut rh1 = 0i64;
        for j in 0..h.len() {
            rh0 += i64::from(h[j]) * i64::from(h[j]);
            if j + 1 < h.len() {
                rh1 += i64::from(h[j]) * i64::from(h[j + 1]);
            }
        }
        let k1_q15: i16 = if rh0 > 0 {
            (-((rh1 << 15) / rh0)).clamp(-32768, 32767) as i16
        } else {
            0
        };
        let gamma_t = if k1_q15 < 0 {
            GAMMA_T_NEG_Q15
        } else {
            GAMMA_T_POS_Q15
        };
        // c = γ_t·k1' (Q15); g_t = 1 − |c| (Word32, avoids the Q15
        // Word16 ceiling at c = 0).
        let c = mult(gamma_t, k1_q15);
        let gt32 = l_sub(32768, i32::from(abs_s(c)));
        let (gt_m, gt_e) = recip_norm(gt32);

        let mut out = [0i16; L_SUBFR];
        for (n, o) in out.iter_mut().enumerate() {
            let prev = if n == 0 { self.tilt_prev } else { t[n - 1] };
            // t(n) + c·t(n−1) on Q16.
            let acc = l_mac(l_deposit_h(t[n]), c, prev);
            *o = if self.lat.tilt_exact_div {
                // Exact rounded (t + c·t₋₁)/g_t on the Q16/Q15 grids.
                let num = (i64::from(acc) << 15) / i64::from(gt32.max(1));
                (((num >> 15) + 1) >> 1).clamp(-32768, 32767) as i16
            } else {
                // ×(1/g_t) = ×m·2^(e−30) keeping Q16: mpy is acc·m/2^15,
                // so shift by e − 15.
                let (hi, lo) = l_extract(acc);
                round(l_shl(mpy_32_16(hi, lo, gt_m), sub(gt_e, 15)))
            };
        }
        self.tilt_prev = t[L_SUBFR - 1];
        out
    }

    /// §4.2.4 adaptive gain control (eqs (88)–(90)) on the Q12 gain
    /// grid; the eq (88) ratio is held at the running gain on a silent
    /// postfiltered subframe.
    fn agc(&mut self, s_hat: &[i16; L_SUBFR], sf: &[i16; L_SUBFR]) -> [i16; L_SUBFR] {
        let mut num = 0i32;
        let mut den = 0i32;
        for n in 0..L_SUBFR {
            num = l_add(num, i32::from(abs_s(s_hat[n])));
            den = l_add(den, i32::from(abs_s(sf[n])));
        }
        let shift = self.lat.agc_shift;
        let (num, den) = if self.lat.agc_lag {
            let prev = self.agc_prev_sums;
            self.agc_prev_sums = (num, den);
            prev
        } else {
            (num, den)
        };
        let g_target_q12: i16 = if den > 0 {
            let ratio = f64::from(num) / f64::from(den);
            let ratio = if self.lat.agc_sqrt {
                ratio.sqrt()
            } else {
                ratio
            };
            ((ratio * ((1i64 << (12 + shift)) as f64)) as i64).min(32767) as i16
        } else {
            self.agc_q12
        };

        let mut out = [0i16; L_SUBFR];
        for (n, o) in out.iter_mut().enumerate() {
            let x2 = i16::from(self.lat.agc_x2);
            let (pole, target_w) = if self.lat.agc_pole_q15 == 0 {
                (AGC_PREV_Q15, AGC_TARGET_Q15)
            } else {
                (self.lat.agc_pole_q15, sub(32767, self.lat.agc_pole_q15))
            };
            if self.lat.agc_apply_first {
                *o = round(l_shl(l_mult(sf[n], self.agc_q12), 3 - shift + x2));
                self.agc_q12 = add(mult(pole, self.agc_q12), mult(target_w, g_target_q12));
            } else {
                // eq (90): g(n) = 0.85·g(n−1) + 0.15·G (Q15 weights).
                self.agc_q12 = add(mult(pole, self.agc_q12), mult(target_w, g_target_q12));
                // eq (89): sf′(n) = g(n)·sf(n) — product to Q16, round
                // (or to Q17 → Q1 when the ×2 upscale is folded in).
                *o = round(l_shl(l_mult(sf[n], self.agc_q12), 3 - shift + x2));
            }
        }
        out
    }

    /// §4.2.5 output high-pass (eq (91), staged Q13 coefficients) with
    /// the ×2 upscale folded into the output rounding; the feedback
    /// memory keeps the unrounded Word32 recursion value.
    fn high_pass(&mut self, x: &[i16; L_SUBFR]) -> [i16; L_SUBFR] {
        let b = &HPF_PREPROC_100HZ_B_Q13;
        let a = &HPF_PREPROC_100HZ_A_Q13;
        let mut out = [0i16; L_SUBFR];
        for (n, o) in out.iter_mut().enumerate() {
            // Q14 accumulator: b (Q13) × x (Q0) doubled.
            let mut acc = l_mult(x[n], b[0]);
            acc = l_mac(acc, self.hp_x[0], b[1]);
            acc = l_mac(acc, self.hp_x[1], b[2]);
            if self.lat.hp_wide {
                // Feedback: a (Q13) × y (Q14 Word32) as an EXACT
                // double-precision product landed on Q14 — corpus-
                // pinned (round 452): with the AGC output on the Q1
                // grid, the whole frame-0 subframe of all six clean
                // vectors is byte-exact only when the low word of the
                // feedback product is carried in full; the (hi, lo)
                // split with the low product truncated to Q15 ends the
                // filter's decay tail one sample late (FIXED n = 16,
                // ALGTHM/PITCH n = 22), and rounding that low product
                // ends it one sample early.
                // Each product lands on Q15 (`>> 15`, the full 48-bit
                // product) and is then realigned by `<< 2`.
                let fb1 = ((i64::from(self.hp_y[0]) * i64::from(a[1])) >> 15) << 2;
                let fb2 = ((i64::from(self.hp_y[1]) * i64::from(a[2])) >> 15) << 2;
                acc = l_add(
                    acc,
                    fb1.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32,
                );
                acc = l_add(
                    acc,
                    fb2.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32,
                );
                self.hp_y[1] = self.hp_y[0];
                self.hp_y[0] = acc;
            } else {
                // Word16-rounded feedback on the Q0 y grid.
                acc = l_mac(acc, extract_h(l_shl(self.hp_y[0], 2)), a[1]);
                acc = l_mac(acc, extract_h(l_shl(self.hp_y[1], 2)), a[2]);
                self.hp_y[1] = self.hp_y[0];
                // Store as Q14 of the rounded Q0 value for grid unity.
                self.hp_y[0] = l_shl(l_deposit_l(round(l_shl(acc, 2))), 14);
            }
            self.hp_x[1] = self.hp_x[0];
            self.hp_x[0] = x[n];
            // ×2 upscale: acc = y·2^14 → 2y·2^16 = acc·2^3, rounded —
            // unless the input already carried the ×2 (Q1 grid), in
            // which case the accumulator is 2y·2^14 and lands by 2^2.
            *o = round(l_shl(acc, if self.lat.agc_x2 { 2 } else { 3 }));
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

        // §4.2.1 H_p(z) on the residual.
        let (hp, decision) = self.long_term(int_t1);

        // 1/[g_f·Â(z/γ_d)] synthesis (impulse response shared with
        // §4.2.3).
        let h40 = impulse_response_q12_n(&apn, &apd);
        let len = if self.lat.gf_len == 0 {
            GF_IMPULSE_LEN
        } else {
            usize::from(self.lat.gf_len).min(40)
        };
        let h = &h40[..len];
        let mut gf = 0i32;
        for &hv in h {
            gf = l_add(gf, i32::from(abs_s(hv)));
        }
        let gf_recip = if gf > 0 { recip_norm(gf) } else { (16384, 33) };
        let hf = self.synthesis(&hp, &apd, gf_recip, gf);

        // §4.2.3 H_t(z).
        let ht = self.tilt(&hf, h);

        // §4.2.4 AGC against the reconstructed speech ŝ(n).
        let agc_gain_in = self.agc_q12;
        let agc = self.agc(s, &ht);

        // §4.2.5 output high-pass + ×2.
        let output = self.high_pass(&agc);
        SubframeTraceFx {
            decision,
            long_term: hp,
            impulse: core::array::from_fn(|n| h40[n]),
            gf_q12: gf,
            short_term: hf,
            tilt: ht,
            agc_gain_in_q12: agc_gain_in,
            agc_gain_out_q12: self.agc_q12,
            agc,
            output,
        }
    }
}

/// Every intermediate §4.2 stage signal for one subframe (the
/// bisection instrument behind
/// [`PostfilterFx::process_subframe_traced`]).
#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct SubframeTraceFx {
    pub decision: LtDecisionFx,
    pub long_term: [i16; L_SUBFR],
    pub impulse: [i16; GF_IMPULSE_LEN],
    pub gf_q12: i32,
    pub short_term: [i16; L_SUBFR],
    pub tilt: [i16; L_SUBFR],
    pub agc_gain_in_q12: i16,
    pub agc_gain_out_q12: i16,
    pub agc: [i16; L_SUBFR],
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

    /// The normalised reciprocal helper: m·2^(e−45) reproduces 1/x
    /// within the div_s truncation for a spread of magnitudes.
    #[test]
    fn recip_norm_inverts() {
        for &x in &[1i32, 5, 4096, 8192, 12345, 1 << 20, 0x7fff_0000] {
            let (m, e) = recip_norm(x);
            let got = f64::from(m) * (f64::from(e) - 45.0).exp2();
            let expect = 1.0 / f64::from(x);
            assert!(
                ((got - expect) / expect).abs() < 1e-3,
                "x={x}: {got} vs {expect}"
            );
        }
    }

    /// Silence propagates: an all-zero subframe with any stable LP set
    /// produces all-zero output and holds the AGC gain.
    #[test]
    fn silence_stays_silent() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2048;
        let mut pf = PostfilterFx::new();
        let g0 = pf.agc_q12;
        let (out, _) = pf.process_subframe(&[0; L_SUBFR], &a, 50);
        assert_eq!(out, [0i16; L_SUBFR]);
        // The eq (90) Q15-truncating recursion may drift the held gain
        // by a few LSBs across a silent subframe; it must not run away.
        assert!(
            (i32::from(pf.agc_q12) - i32::from(g0)).abs() <= 4,
            "AGC drifted {} → {}",
            g0,
            pf.agc_q12
        );
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
        }
    }
}
