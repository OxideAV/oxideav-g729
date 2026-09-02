//! §3.3 **perceptual weighting** on the fixed grid: the eq (28)–(32)
//! γ adaptation and the eq (33) weighting filter
//! `W(z) = A(z/γ₁)/A(z/γ₂)` producing the weighted speech `sw(n)`.
//!
//! ## Number grids
//!
//! - Reflection coefficients `k₁, k₂` — Q15 (the §3.2.2 by-product).
//! - Log-area ratios — the eq (28) logarithm is evaluated through the
//!   clause-5 `log2` table function as `log2(1+k) − log2(1−k)` on a
//!   Word32 Q15 grid; the eq (30) thresholds are carried on the same
//!   base-2 grid (`o / ln 2`), so no natural-log conversion is needed
//!   for the comparisons.
//! - `γ₁, γ₂` — Q15 (`0.94 → 30802`, `0.6 → 19661`, `0.98 → 32113`,
//!   the eq (32) bounds `0.4 → 13107`, `0.7 → 22938`).
//! - LSF for the eq (31) resonance criterion — Q13 radians.
//! - Speech / weighted speech — Q0.
//!
//! ## eq (33) as a two-stage run
//!
//! The weighting filter is realised as the eq (36)-shaped FIR
//! `A(z/γ₁)` on the speech followed by the all-pole `1/A(z/γ₂)` with
//! its own output memory — each stage landing on Q0 with the shared
//! `<< 3` / round protocol of [`crate::fx::filters`]. The single
//! direct-form recursion of the float module is the same transfer
//! function; the intermediate Q0 landing is the fixed chain's own
//! composition (the prose is silent on it).

use crate::fx::dsp::log2;
use crate::fx::filters::{residu_lat, syn_filt_mem_lat, FilterLatitude, L_SUBFR};
use crate::fx::ops::{l_add, l_shr, l_sub};
use crate::tables::M;

/// `γ₁` for a flat spectral envelope, Q15 (0.94).
pub const GAMMA1_FLAT_Q15: i16 = 30802;
/// `γ₂` for a flat spectral envelope, Q15 (0.6).
pub const GAMMA2_FLAT_Q15: i16 = 19661;
/// `γ₁` for a tilted spectral envelope, Q15 (0.98).
pub const GAMMA1_TILT_Q15: i16 = 32113;
/// eq (32) lower bound on `γ₂`, Q15 (0.4).
pub const GAMMA2_MIN_Q15: i16 = 13107;
/// eq (32) upper bound on `γ₂`, Q15 (0.7).
pub const GAMMA2_MAX_Q15: i16 = 22938;

/// eq (30) thresholds on the base-2 LAR grid (Q15) when the eq (28)
/// logarithm is read as the natural log: the prose values divided by
/// `ln 2`.
const THR_LN_Q15: [i32; 4] = [
    -82_257, // −1.74 / ln 2
    30_728,  // 0.65 / ln 2
    -71_857, // −1.52 / ln 2
    20_328,  // 0.43 / ln 2
];
/// The same thresholds when eq (28) is read as the base-10 logarithm
/// (the convention the clause's `20 log g_c` energy equations use):
/// the prose values times `log2 10`.
const THR_LOG10_Q15: [i32; 4] = [
    -189_400, // −1.74 · log2 10
    70_754,   // 0.65 · log2 10
    -165_452, // −1.52 · log2 10
    46_806,   // 0.43 · log2 10
];

/// Unstated latitude of the §3.3 decision arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WeightingLatitude {
    /// Negate the reflection coefficients before eq (28) (the §3.2.2
    /// sign convention of `k_i` is not restated in §3.3).
    pub negate_rc: bool,
    /// eq (32) slope: `γ₂ = 1 − slope·d_min` with `d_min` on Q13
    /// radians, the slope on a `2^8` grid (`6.0` in radians = 1536;
    /// `6.0` per unit of `ω/π` = 489; per unit of `ω/2π` = 244).
    pub gamma2_slope_q8: i32,
    /// Evaluate the eq (28)/(30) decision through the spec-equation
    /// float module instead of the table `log2`.
    pub float_decision: bool,
    /// Read the eq (28) `log` as base 10 (thresholds scaled by
    /// `log2 10`) instead of the natural log.
    pub log_base10: bool,
}

impl Default for WeightingLatitude {
    fn default() -> Self {
        Self {
            negate_rc: false,
            gamma2_slope_q8: 1536,
            float_decision: false,
            log_base10: true,
        }
    }
}

/// The per-subframe weighting decision on the Q15 grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WeightingDecisionFx {
    /// `γ₁` on Q15.
    pub gamma1_q15: i16,
    /// `γ₂` on Q15.
    pub gamma2_q15: i16,
    /// eq (30) classification (`true` = flat).
    pub flat: bool,
}

/// §3.3 state: the γ-adaptation memory plus the eq (33) denominator
/// memory (the numerator's speech history is the caller's analysis
/// buffer).
#[derive(Debug, Clone)]
pub struct WeightingFx {
    /// Previous frame's LAR pair on the base-2 Q15 grid.
    prev_lar: [i32; 2],
    /// Previous subframe's eq (30) classification.
    prev_flat: bool,
    /// `sw(n−1) … sw(n−M)` (most recent last).
    sw_mem: [i16; M],
    lat: WeightingLatitude,
    /// Float-decision fallback state (latitude probe).
    float_adapt: crate::perceptual_weighting::PerceptualWeighting,
}

impl Default for WeightingFx {
    fn default() -> Self {
        Self::new()
    }
}

/// eq (28) on the base-2 Q15 grid: `log2((1 + k)/(1 − k))` for a Q15
/// reflection coefficient, as `log2(1 + k) − log2(1 − k)` through the
/// table function (both operands positive Word32 on the Q15 grid, so
/// the exponent difference cancels the common scale).
#[must_use]
pub fn lar_log2_q15(k_q15: i16) -> i32 {
    let k = i32::from(k_q15).clamp(-32_764, 32_764);
    let num = 32_768 + k;
    let den = 32_768 - k;
    let (en, fn_) = log2(num);
    let (ed, fd) = log2(den);
    l_add(
        i32::from(en - ed) << 15,
        l_sub(i32::from(fn_), i32::from(fd)),
    )
}

impl WeightingFx {
    /// Fresh state (zero LARs and memories, `flat = 1`).
    #[must_use]
    pub fn new() -> Self {
        Self::with_latitude(WeightingLatitude::default())
    }

    /// Fresh state under an explicit latitude.
    #[must_use]
    pub fn with_latitude(lat: WeightingLatitude) -> Self {
        Self {
            prev_lar: [0; 2],
            prev_flat: true,
            sw_mem: [0; M],
            lat,
            float_adapt: crate::perceptual_weighting::PerceptualWeighting::new(),
        }
    }

    /// eqs (28)–(32) for both subframes of one frame. `rc_q15` is the
    /// `(k₁, k₂)` pair, `lsf_sub` each subframe's unquantised LSF
    /// vector on Q13.
    pub fn adapt_frame(
        &mut self,
        rc_q15: [i16; 2],
        lsf_sub: &[[i16; M]; 2],
    ) -> [WeightingDecisionFx; 2] {
        let rc = if self.lat.negate_rc {
            [
                crate::fx::ops::negate(rc_q15[0]),
                crate::fx::ops::negate(rc_q15[1]),
            ]
        } else {
            rc_q15
        };
        if self.lat.float_decision {
            return self.adapt_frame_float(rc, lsf_sub);
        }
        let cur = [lar_log2_q15(rc[0]), lar_log2_q15(rc[1])];
        // eq (29): subframe 1 = mean of previous and current LARs.
        let lar_sub1 = [
            l_shr(l_add(self.prev_lar[0], cur[0]), 1),
            l_shr(l_add(self.prev_lar[1], cur[1]), 1),
        ];
        self.prev_lar = cur;

        let mut out = [WeightingDecisionFx {
            gamma1_q15: GAMMA1_FLAT_Q15,
            gamma2_q15: GAMMA2_FLAT_Q15,
            flat: true,
        }; 2];
        for (sub, lars) in [lar_sub1, cur].iter().enumerate() {
            // eq (30) hysteresis.
            let thr = if self.lat.log_base10 {
                THR_LOG10_Q15
            } else {
                THR_LN_Q15
            };
            let flat = if self.prev_flat && lars[0] < thr[0] && lars[1] > thr[1] {
                false
            } else if !self.prev_flat && (lars[0] > thr[2] || lars[1] < thr[3]) {
                true
            } else {
                self.prev_flat
            };
            self.prev_flat = flat;
            out[sub] = if flat {
                WeightingDecisionFx {
                    gamma1_q15: GAMMA1_FLAT_Q15,
                    gamma2_q15: GAMMA2_FLAT_Q15,
                    flat: true,
                }
            } else {
                WeightingDecisionFx {
                    gamma1_q15: GAMMA1_TILT_Q15,
                    gamma2_q15: self.gamma2_tilted(&lsf_sub[sub]),
                    flat: false,
                }
            };
        }
        out
    }

    /// eqs (31)/(32): `γ₂ = 1 − slope · d_min`, bounded to `[0.4, 0.7]`.
    fn gamma2_tilted(&self, lsf_q13: &[i16; M]) -> i16 {
        let d_min = lsf_q13
            .windows(2)
            .map(|w| i32::from(w[1]) - i32::from(w[0]))
            .min()
            .unwrap_or(0)
            .max(0);
        // Q13 · Q8 slope → Q21; back to Q15.
        let drop = (d_min * self.lat.gamma2_slope_q8) >> 6;
        let g2 = 32_767 - drop;
        g2.clamp(i32::from(GAMMA2_MIN_Q15), i32::from(GAMMA2_MAX_Q15)) as i16
    }

    fn adapt_frame_float(
        &mut self,
        rc: [i16; 2],
        lsf_sub: &[[i16; M]; 2],
    ) -> [WeightingDecisionFx; 2] {
        let k = (f32::from(rc[0]) / 32768.0, f32::from(rc[1]) / 32768.0);
        let omega: [[f32; M]; 2] = [
            std::array::from_fn(|i| f32::from(lsf_sub[0][i]) / 8192.0),
            std::array::from_fn(|i| f32::from(lsf_sub[1][i]) / 8192.0),
        ];
        let d = self.float_adapt.adapt_frame(k, &omega);
        std::array::from_fn(|s| WeightingDecisionFx {
            gamma1_q15: q15(d[s].gamma1),
            gamma2_q15: if d[s].flat {
                GAMMA2_FLAT_Q15
            } else {
                self.gamma2_tilted(&lsf_sub[s])
            },
            flat: d[s].flat,
        })
    }

    /// eq (33) over one subframe: `s_hist` holds `s(−M) … s(−1)`
    /// (most recent last), `ap1`/`ap2` the γ-weighted coefficient sets
    /// (Q12, from [`crate::fx::filters::weight_az`]).
    #[must_use]
    pub fn weight_subframe(
        &mut self,
        s_hist: &[i16; M],
        s: &[i16; L_SUBFR],
        ap1: &[i16; M + 1],
        ap2: &[i16; M + 1],
        lat: &FilterLatitude,
    ) -> [i16; L_SUBFR] {
        let r = residu_lat(ap1, s_hist, s, lat.residu_trunc);
        let sw = syn_filt_mem_lat(ap2, &self.sw_mem, &r, lat.syn_trunc);
        self.sw_mem.copy_from_slice(&sw[L_SUBFR - M..]);
        sw
    }
}

/// Round-to-nearest Q15 landing of a unit-range factor (saturating).
fn q15(x: f32) -> i16 {
    (f64::from(x) * 32768.0).round().clamp(-32768.0, 32767.0) as i16
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fx::filters::weight_az;

    /// With `γ₁ = γ₂` the weighting filter is the identity on the
    /// integer grid (up to the two Q0 landings).
    #[test]
    fn equal_gammas_identity() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[1] = -2048;
        a[2] = 819;
        let ap = weight_az(&a, 26214);
        let mut w = WeightingFx::new();
        let hist = [0i16; M];
        let s: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 13 % 100) as i16) * 40 - 2000);
        let sw = w.weight_subframe(&hist, &s, &ap, &ap, &FilterLatitude::default());
        for n in 0..L_SUBFR {
            assert!((i32::from(sw[n]) - i32::from(s[n])).abs() <= 2, "n={n}");
        }
    }

    /// The Q15 landings of the fixed γ pairs.
    #[test]
    fn gamma_grid() {
        assert_eq!(q15(0.94), GAMMA1_FLAT_Q15);
        assert_eq!(q15(0.6), GAMMA2_FLAT_Q15);
        assert_eq!(q15(0.98), GAMMA1_TILT_Q15);
        assert_eq!(q15(0.4), GAMMA2_MIN_Q15);
        assert_eq!(q15(0.7), GAMMA2_MAX_Q15);
    }

    /// The table LAR tracks `log2((1+k)/(1−k))` within the
    /// interpolation tolerance, is odd, and zero at `k = 0`.
    #[test]
    fn lar_tracks_log2() {
        assert_eq!(lar_log2_q15(0), 0);
        for k in [-30000i16, -20000, -9000, -1000, 1000, 9000, 20000, 30000] {
            let got = f64::from(lar_log2_q15(k)) / 32768.0;
            let kf = f64::from(k) / 32768.0;
            let want = ((1.0 + kf) / (1.0 - kf)).log2();
            assert!((got - want).abs() < 4e-3, "k={k}: {got} vs {want}");
            assert_eq!(lar_log2_q15(k), -lar_log2_q15(-k));
        }
    }

    /// eq (30) hysteresis on the fixed grid (base-10 thresholds):
    /// strongly tilted LARs flip to tilted on subframe 2 of the first
    /// frame (subframe 1 is interpolated with the zero start-up pair),
    /// flat LARs release.
    #[test]
    fn flat_hysteresis_transitions() {
        let mut w = WeightingFx::new();
        let lsf: [i16; M] = std::array::from_fn(|i| (2340 * (i + 1)) as i16);
        // k₁ = −0.97 → o₁ = log10(0.03/1.97) ≈ −1.82 < −1.74;
        // k₂ = 0.7 → o₂ = log10(1.7/0.3) ≈ 0.75 > 0.65 ⇒ tilted.
        let d = w.adapt_frame([-31785, 22938], &[lsf, lsf]);
        assert!(
            d[0].flat,
            "subframe 1 interpolates with the zero start-up LARs"
        );
        assert!(!d[1].flat);
        assert_eq!(d[1].gamma1_q15, GAMMA1_TILT_Q15);
        // d_min = 2340/8192 = 0.286 rad → γ₂ = 1 − 1.71 → clamp 0.4.
        assert_eq!(d[1].gamma2_q15, GAMMA2_MIN_Q15);
        // k₁ = 0 → o₁ = 0 > −1.52 releases the hysteresis.
        let d2 = w.adapt_frame([0, 22938], &[lsf, lsf]);
        assert!(d2[0].flat);
        assert_eq!(d2[0].gamma2_q15, GAMMA2_FLAT_Q15);
    }

    /// eq (32): a tight LSF pair pushes γ₂ to the 0.7 bound.
    #[test]
    fn gamma2_resonance_bound() {
        let w = WeightingFx::new();
        let mut lsf: [i16; M] = std::array::from_fn(|i| (2340 * (i + 1)) as i16);
        lsf[5] = lsf[4] + 80; // ≈ 0.01 rad
        assert_eq!(w.gamma2_tilted(&lsf), GAMMA2_MAX_Q15);
    }
}
