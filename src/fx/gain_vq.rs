//! §3.9 **gain quantisation** and the taming procedure on the fixed
//! grid.
//!
//! ## Search structure (clause 3.9.2)
//!
//! The optimal `(g_p, g_c)` of eq (63) preselect a cluster of **four**
//! consecutive `GA` rows (the codebook is staged in ascending `γ`
//! order) and a cluster of **eight** consecutive `GB` rows (ascending
//! `g_p` order); the 4 × 8 survivors are searched exhaustively. The
//! cluster starts are read off the staged partial-search threshold
//! tables (`gain-quantizer-codebook-GA-thresholds-Q14`, four entries on
//! the `γ` axis; `…-GB-thresholds-Q15`, eight entries on the `g_p`
//! axis): the start index is the number of thresholds the optimal
//! value exceeds. The first GA threshold (`10808`, Q14 `0.6597`) is
//! exactly twice the fourth GA row's `γ` (`5404`), i.e. the boundary
//! between "rows 0–3" and "rows 1–4" sits on the last row of the
//! lower cluster — the tables are the clause's "closest" rule
//! precomputed.
//!
//! ## Number grids
//!
//! - Correlations `yᵗy, zᵗz, xᵗy, xᵗz, yᵗz` — exact wide sums of the
//!   Q0 target / Q0 filtered adaptive vector / Q12 filtered codevector.
//! - Candidate gains — the decoder's own reconstruction
//!   ([`crate::fx::gains::GainDecoderFx::reconstruct`]): `ĝ_p` Q14,
//!   `ĝ_c` Q1.
//! - eq (63) — evaluated per candidate on a wide accumulator from the
//!   integer correlations and the Q14/Q1 gains (the `xᵗx` term is
//!   common to all candidates and dropped).
//! - Preselection targets — the eq (63) optimum solved in double
//!   precision from the wide correlations; `γ_opt = g_c/g'_c` compared
//!   on the Q14 threshold grid, `g_p` on Q15.
//!
//! The prose pins the codebook structure, the cluster sizes and the
//! criterion; the threshold tables are staged data. The comparison
//! sense at a threshold and the precision of the error evaluation are
//! the fixed chain's own composition, exposed as latitude.

use crate::fx::filters::L_SUBFR;
use crate::fx::gains::{DecodedGainsFx, GainDecoderFx, GainPredictionFx};
use crate::tables::{GAIN_QUANT_GA_THRESHOLDS_Q14, GAIN_QUANT_GB_THRESHOLDS_Q15, NCODE1, NCODE2};
use crate::taming::Taming;

/// GA preselection cluster size (clause 3.9.2: "a cluster of four").
pub const GA_CANDIDATES: usize = 4;
/// GB preselection cluster size (clause 3.9.2: "a cluster of eight").
pub const GB_CANDIDATES: usize = 8;

/// Which "optimum gains" (clause 3.9.2) drive the preselection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PreselectTarget {
    /// The joint eq (63) minimiser (2 × 2 normal equations).
    #[default]
    Joint,
    /// The single-axis optima `xᵗy/yᵗy` and `xᵗz/zᵗz`.
    Axis,
    /// The sequential optima: `g_p = xᵗy/yᵗy`, then `g_c` on the
    /// eq (50) updated target `(x − g_p·y)ᵗz / zᵗz`.
    Sequential,
}

/// Unstated latitude of the §3.9.2 arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GainVqLatitude {
    /// The preselection's optimum-gain definition.
    pub target: PreselectTarget,
    /// Search the full 8 × 16 grid instead of the preselected 4 × 8
    /// cluster (oracle for the preselection's effect).
    pub exhaustive: bool,
    /// A threshold counts as exceeded on equality (`>=`) instead of
    /// strictly.
    pub threshold_ge: bool,
    /// Preselect by nearest-row ranking (the float module's reading)
    /// instead of the staged threshold tables.
    pub rank_select: bool,
}

impl Default for GainVqLatitude {
    fn default() -> Self {
        Self {
            target: PreselectTarget::Joint,
            exhaustive: false,
            threshold_ge: false,
            rank_select: false,
        }
    }
}

/// The eq (63) inner products of one subframe (exact wide sums).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GainCorrelationsFx {
    /// `yᵗy` (Q0 × Q0).
    pub yy: i64,
    /// `zᵗz` (Q12 × Q12 → Q24).
    pub zz: i64,
    /// `xᵗy` (Q0 × Q0).
    pub xy: i64,
    /// `xᵗz` (Q0 × Q12 → Q12).
    pub xz: i64,
    /// `yᵗz` (Q0 × Q12 → Q12).
    pub yz: i64,
}

impl GainCorrelationsFx {
    /// From the Q0 target `x`, the Q0 filtered adaptive vector `y`
    /// and the Q12 filtered codevector `z`.
    #[must_use]
    pub fn compute(x: &[i16; L_SUBFR], y: &[i16; L_SUBFR], z_q12: &[i16; L_SUBFR]) -> Self {
        let mut t = Self {
            yy: 0,
            zz: 0,
            xy: 0,
            xz: 0,
            yz: 0,
        };
        for n in 0..L_SUBFR {
            let (xn, yn, zn) = (i64::from(x[n]), i64::from(y[n]), i64::from(z_q12[n]));
            t.yy += yn * yn;
            t.zz += zn * zn;
            t.xy += xn * yn;
            t.xz += xn * zn;
            t.yz += yn * zn;
        }
        t
    }

    /// eq (63) without the constant `xᵗx`, on the Q0 signal scale
    /// (double precision of exact integer terms), for `ĝ_p` on Q14
    /// and `ĝ_c` on Q1.
    #[must_use]
    pub fn error(&self, gains: DecodedGainsFx) -> f64 {
        let gp = f64::from(gains.gain_pit_q14) / 16384.0;
        // ĝ_c·z on Q0 = gc_q1 · z_q12 / 2^13.
        let gc = f64::from(gains.gain_code_q1) / 8192.0;
        gp * gp * self.yy as f64 + gc * gc * self.zz as f64
            - 2.0 * gp * self.xy as f64
            - 2.0 * gc * self.xz as f64
            + 2.0 * gp * gc * self.yz as f64
    }

    /// The eq (63) optimum `(g_p, g_c)` on the natural scale (`g_c`
    /// such that `g_c · z_q12 / 2^12` is the Q0 contribution).
    #[must_use]
    pub fn optimal(&self) -> (f64, f64) {
        self.optimal_for(PreselectTarget::Joint)
    }

    /// The preselection optimum under the chosen definition.
    #[must_use]
    pub fn optimal_for(&self, target: PreselectTarget) -> (f64, f64) {
        let (yy, zz, xy, xz, yz) = (
            self.yy as f64,
            self.zz as f64 / 16_777_216.0,
            self.xy as f64,
            self.xz as f64 / 4096.0,
            self.yz as f64 / 4096.0,
        );
        let gp_axis = if yy > 0.0 { xy / yy } else { 0.0 };
        let gc_axis = if zz > 0.0 { xz / zz } else { 0.0 };
        match target {
            PreselectTarget::Joint => {
                let det = yy * zz - yz * yz;
                if det.abs() > 1e-9 * (yy * zz).abs().max(1.0) {
                    ((xy * zz - xz * yz) / det, (xz * yy - xy * yz) / det)
                } else {
                    (gp_axis, gc_axis)
                }
            }
            PreselectTarget::Axis => (gp_axis, gc_axis),
            PreselectTarget::Sequential => {
                let gp = gp_axis.clamp(0.0, 1.2);
                let gc = if zz > 0.0 { (xz - gp * yz) / zz } else { 0.0 };
                (gp, gc)
            }
        }
    }
}

/// The selected codebook-domain index pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GainQuantFx {
    /// First-stage index (`0 … 7`).
    pub ga: usize,
    /// Second-stage index (`0 … 15`).
    pub gb: usize,
}

/// Cluster start from a threshold table: the number of thresholds the
/// value exceeds.
fn cluster_start(value: i32, thresholds: &[i16], ge: bool) -> usize {
    thresholds
        .iter()
        .filter(|&&t| {
            if ge {
                value >= i32::from(t)
            } else {
                value > i32::from(t)
            }
        })
        .count()
}

/// §3.9.2 search for one subframe under the default latitude. `pred`
/// is the decoder-consistent eq (71) prediction; `tame` the taming
/// flag.
#[must_use]
pub fn quantize_gains_fx(
    corr: &GainCorrelationsFx,
    dec: &GainDecoderFx,
    pred: &GainPredictionFx,
    tame: bool,
) -> GainQuantFx {
    quantize_gains_fx_lat(corr, dec, pred, tame, &GainVqLatitude::default())
}

/// §3.9.2 search under an explicit latitude.
#[must_use]
pub fn quantize_gains_fx_lat(
    corr: &GainCorrelationsFx,
    dec: &GainDecoderFx,
    pred: &GainPredictionFx,
    tame: bool,
    lat: &GainVqLatitude,
) -> GainQuantFx {
    // Preselection targets.
    let (gp_opt, gc_opt) = corr.optimal_for(lat.target);
    let g_c_prime = f64::from(pred.g_prime_scaled) * 2f64.powi(i32::from(pred.exp) - 13);
    let gamma_opt = if g_c_prime > 1e-9 {
        gc_opt / g_c_prime
    } else {
        gc_opt
    };

    let (ga_set, gb_set): (Vec<usize>, Vec<usize>) = if lat.exhaustive {
        ((0..NCODE1).collect(), (0..NCODE2).collect())
    } else if lat.rank_select {
        use crate::tables::{gain_ga_entry, gain_gb_entry, GAIN_VQ_COL_GC, GAIN_VQ_COL_GP};
        let mut ga: Vec<usize> = (0..NCODE1).collect();
        ga.sort_by(|&a, &b| {
            let da = (f64::from(gain_ga_entry(a)[GAIN_VQ_COL_GC]) / 8192.0 - gamma_opt).abs();
            let db = (f64::from(gain_ga_entry(b)[GAIN_VQ_COL_GC]) / 8192.0 - gamma_opt).abs();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut gb: Vec<usize> = (0..NCODE2).collect();
        gb.sort_by(|&a, &b| {
            let da = (f64::from(gain_gb_entry(a)[GAIN_VQ_COL_GP]) / 16384.0 - gp_opt).abs();
            let db = (f64::from(gain_gb_entry(b)[GAIN_VQ_COL_GP]) / 16384.0 - gp_opt).abs();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });
        (ga[..GA_CANDIDATES].to_vec(), gb[..GB_CANDIDATES].to_vec())
    } else {
        // γ_opt on the Q14 threshold grid, g_p on Q15.
        let gamma_q14 = (gamma_opt * 16384.0).round().clamp(-1.0e9, 1.0e9) as i32;
        let gp_q15 = (gp_opt * 32768.0).round().clamp(-1.0e9, 1.0e9) as i32;
        let sa = cluster_start(gamma_q14, &GAIN_QUANT_GA_THRESHOLDS_Q14, lat.threshold_ge);
        let sb = cluster_start(gp_q15, &GAIN_QUANT_GB_THRESHOLDS_Q15, lat.threshold_ge);
        (
            (sa..sa + GA_CANDIDATES).collect(),
            (sb..sb + GB_CANDIDATES).collect(),
        )
    };

    let score = |ga_set: &[usize], gb_set: &[usize]| -> Option<(f64, GainQuantFx)> {
        let mut best: Option<(f64, GainQuantFx)> = None;
        for &ga in ga_set {
            for &gb in gb_set {
                let g = dec.reconstruct(pred, ga, gb);
                if tame && g.gain_pit_q14 > crate::fx::gain_vq::GPCLIP_Q14 {
                    continue;
                }
                let e = corr.error(g);
                let better = match &best {
                    None => true,
                    Some((be, _)) => e < *be,
                };
                if better {
                    best = Some((e, GainQuantFx { ga, gb }));
                }
            }
        }
        best
    };

    if let Some((_, r)) = score(&ga_set, &gb_set) {
        return r;
    }
    // Taming excluded every preselected pair: the full grid under the
    // same restriction.
    let all_a: Vec<usize> = (0..NCODE1).collect();
    let all_b: Vec<usize> = (0..NCODE2).collect();
    score(&all_a, &all_b)
        .expect("codebooks carry low-gain rows, a taming-legal pair exists")
        .1
}

/// Taming ceiling on the quantised pitch gain (Q14 `15564` ≈ 0.95,
/// doc `taming-procedure.md` §4).
pub const GPCLIP_Q14: i16 = 15564;

/// Taming-procedure state (doc `taming-procedure.md`) driven with Q14
/// pitch gains.
#[derive(Debug, Clone, Default)]
pub struct TamingFx {
    inner: Taming,
}

impl TamingFx {
    /// Fresh state.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Doc §3 test for the candidate delay.
    #[must_use]
    pub fn test(&self, int_t: i32, frac: i32) -> bool {
        self.inner.test(int_t, frac)
    }

    /// Doc §2 update with the quantised Q14 pitch gain.
    pub fn update(&mut self, int_t: i32, frac: i32, gain_pit_q14: i16) {
        self.inner
            .update(int_t, frac, f32::from(gain_pit_q14) / 16384.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The staged GA thresholds partition the γ axis into the five
    /// four-row clusters in table order; the first threshold is twice
    /// the fourth row's γ.
    #[test]
    fn ga_thresholds_partition_in_order() {
        let t = &GAIN_QUANT_GA_THRESHOLDS_Q14;
        assert_eq!(t.len(), NCODE1 - GA_CANDIDATES);
        assert!(t.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(
            i32::from(t[0]),
            2 * i32::from(crate::tables::gain_ga_entry(3)[1])
        );
        assert_eq!(cluster_start(0, t, false), 0);
        assert_eq!(cluster_start(i32::MAX, t, false), 4);
        assert_eq!(cluster_start(i32::from(t[1]), t, false), 1);
        assert_eq!(cluster_start(i32::from(t[1]), t, true), 2);
    }

    /// The staged GB thresholds are ascending and count nine clusters.
    #[test]
    fn gb_thresholds_partition_in_order() {
        let t = &GAIN_QUANT_GB_THRESHOLDS_Q15;
        assert_eq!(t.len(), NCODE2 - GB_CANDIDATES);
        assert!(t.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(cluster_start(i32::MAX, t, false), 8);
    }

    /// The exhaustive search over a target built from a codebook pair
    /// recovers that pair; the preselected search agrees with it.
    #[test]
    fn recovers_codebook_pair() {
        let dec = GainDecoderFx::new();
        let y: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 3 % 17) as i16) * 100 - 800);
        let mut code = [0i16; L_SUBFR];
        code[4] = 8192;
        code[13] = -8192;
        code[27] = 8192;
        code[38] = -8192;
        let h: [i16; L_SUBFR] = std::array::from_fn(|n| (4096.0 * 0.6f64.powi(n as i32)) as i16);
        let z = crate::fx::filters::convolve_code_q12(&code, &h);
        let pred = dec.predict(&code);
        for (ga0, gb0) in [(2usize, 5usize), (4, 9), (6, 13), (1, 3)] {
            let g0 = dec.reconstruct(&pred, ga0, gb0);
            let x: [i16; L_SUBFR] = std::array::from_fn(|n| {
                let v = f64::from(g0.gain_pit_q14) / 16384.0 * f64::from(y[n])
                    + f64::from(g0.gain_code_q1) / 8192.0 * f64::from(z[n]);
                v.round() as i16
            });
            let corr = GainCorrelationsFx::compute(&x, &y, &z);
            let ex = quantize_gains_fx_lat(
                &corr,
                &dec,
                &pred,
                false,
                &GainVqLatitude {
                    exhaustive: true,
                    ..Default::default()
                },
            );
            assert_eq!((ex.ga, ex.gb), (ga0, gb0), "exhaustive");
            let pre = quantize_gains_fx(&corr, &dec, &pred, false);
            assert_eq!((pre.ga, pre.gb), (ga0, gb0), "preselected");
        }
    }

    /// With the taming flag the result never reconstructs above the
    /// ceiling.
    #[test]
    fn taming_bounds_pitch_gain() {
        let dec = GainDecoderFx::new();
        let y: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 5 % 31) as i16) * 60 - 900);
        let mut code = [0i16; L_SUBFR];
        code[0] = 8192;
        code[11] = 8192;
        code[22] = -8192;
        code[33] = 8192;
        let h: [i16; L_SUBFR] = std::array::from_fn(|n| (4096.0 * 0.5f64.powi(n as i32)) as i16);
        let z = crate::fx::filters::convolve_code_q12(&code, &h);
        let pred = dec.predict(&code);
        let x: [i16; L_SUBFR] = std::array::from_fn(|n| (f64::from(y[n]) * 1.15) as i16);
        let corr = GainCorrelationsFx::compute(&x, &y, &z);
        let free = quantize_gains_fx(&corr, &dec, &pred, false);
        assert!(dec.reconstruct(&pred, free.ga, free.gb).gain_pit_q14 > GPCLIP_Q14);
        let tamed = quantize_gains_fx(&corr, &dec, &pred, true);
        assert!(dec.reconstruct(&pred, tamed.ga, tamed.gb).gain_pit_q14 <= GPCLIP_Q14);
    }
}
