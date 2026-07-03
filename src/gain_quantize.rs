//! §3.9.2 **gain-quantisation codebook search** — the encoder-side
//! two-stage conjugate-structure VQ that selects the `(GA, GB)`
//! codebook indices for a subframe's adaptive- and fixed-codebook
//! gains.
//!
//! Eleventh encoder-side stage. Consumes the §3.6 target `x(n)`, the
//! eq (44) filtered adaptive vector `y(n)`, the filtered fixed
//! codevector `z(n)` (eq (64): `c(n)` convolved with `h(n)`), and the
//! §3.9.1 predicted gain `g′_c` (the existing decoder-consistent
//! [`crate::gain_predict::GainPredictor`]).
//!
//! ## Spec source — clauses 3.9, 3.9.2, equations (63), (73), (74)
//!
//! * **eq (63)** — the weighted MSE the search minimises:
//!   ```text
//!   E = xᵗx + g_p²·yᵗy + g_c²·zᵗz − 2g_p·xᵗy − 2g_c·xᵗz + 2g_p·g_c·yᵗz
//!   ```
//! * The **optimal (unquantised) gains** minimising eq (63) solve the
//!   2×2 normal equations `[yᵗy yᵗz; yᵗz zᵗz]·[g_p; g_c] = [xᵗy; xᵗz]`
//!   and are used for the conjugate-structure preselection.
//! * **Preselection** (clause 3.9.2 prose): GA's second element (the
//!   `γ` half) is generally larger — select the **4** of 8 GA vectors
//!   whose second elements are closest to the optimal correction factor
//!   (the optimal `g_c` expressed in the correction domain,
//!   `γ_opt = g_c/g′_c`); GB's first element (`g_p` half) is biased —
//!   select the **8** of 16 GB vectors whose first elements are closest
//!   to the optimal `g_p`. Then exhaustively search the 4 × 8 = 32
//!   survivors.
//! * **eqs (73)/(74)** — reconstruction: `ĝ_p = GA[ga].gp + GB[gb].gp`,
//!   `γ̂ = GA[ga].gc + GB[gb].gc`, `ĝ_c = γ̂·g′_c` — evaluated through
//!   the existing decode-side [`crate::gain_reconstruct`], so the
//!   search scores exactly what the decoder will reconstruct.
//!
//! The returned indices are **codebook-domain**; the §3.9.3
//! transmission mapping ([`crate::gain_index_map::map_ga`] /
//! [`crate::gain_index_map::map_gb`]) is applied by the bitstream
//! writer.

use crate::gain_reconstruct::reconstruct_gains;
use crate::tables::{gain_ga_entry, gain_gb_entry, GAIN_VQ_COL_GC, GAIN_VQ_COL_GP, NCODE1, NCODE2};

/// Samples per subframe.
pub const L_SUBFR: usize = 40;
/// GA preselection cluster size (clause 3.9.2: "a cluster of four").
pub const GA_CANDIDATES: usize = 4;
/// GB preselection cluster size (clause 3.9.2: "a cluster of eight").
pub const GB_CANDIDATES: usize = 8;

/// Q14 unit for the codebook `g_p` column.
const Q14: f32 = 16384.0;
/// Q12 unit for the codebook `γ` column.
const Q12: f32 = 4096.0;

/// The six inner products of eq (63), precomputed once per subframe.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GainTerms {
    /// `xᵗx`.
    pub xx: f32,
    /// `yᵗy`.
    pub yy: f32,
    /// `zᵗz`.
    pub zz: f32,
    /// `xᵗy`.
    pub xy: f32,
    /// `xᵗz`.
    pub xz: f32,
    /// `yᵗz`.
    pub yz: f32,
}

impl GainTerms {
    /// Computes the six inner products from the target `x`, the
    /// filtered adaptive vector `y`, and the filtered fixed vector `z`.
    #[must_use]
    pub fn compute(x: &[f32; L_SUBFR], y: &[f32; L_SUBFR], z: &[f32; L_SUBFR]) -> Self {
        let mut t = Self {
            xx: 0.0,
            yy: 0.0,
            zz: 0.0,
            xy: 0.0,
            xz: 0.0,
            yz: 0.0,
        };
        for n in 0..L_SUBFR {
            t.xx += x[n] * x[n];
            t.yy += y[n] * y[n];
            t.zz += z[n] * z[n];
            t.xy += x[n] * y[n];
            t.xz += x[n] * z[n];
            t.yz += y[n] * z[n];
        }
        t
    }

    /// eq (63) evaluated at `(g_p, g_c)`.
    #[must_use]
    pub fn error(&self, gp: f32, gc: f32) -> f32 {
        self.xx + gp * gp * self.yy + gc * gc * self.zz - 2.0 * gp * self.xy - 2.0 * gc * self.xz
            + 2.0 * gp * gc * self.yz
    }

    /// The unquantised `(g_p, g_c)` minimising eq (63) — the solution
    /// of the 2×2 normal equations. Falls back to the axis solutions
    /// when the Gram matrix is (near-)singular.
    #[must_use]
    pub fn optimal(&self) -> (f32, f32) {
        let det = self.yy * self.zz - self.yz * self.yz;
        if det.abs() > 1e-12 {
            let gp = (self.xy * self.zz - self.xz * self.yz) / det;
            let gc = (self.xz * self.yy - self.xy * self.yz) / det;
            (gp, gc)
        } else {
            (
                if self.yy > 0.0 {
                    self.xy / self.yy
                } else {
                    0.0
                },
                if self.zz > 0.0 {
                    self.xz / self.zz
                } else {
                    0.0
                },
            )
        }
    }
}

/// Result of the §3.9.2 search: the selected codebook-domain indices
/// and the decoder-consistent reconstructed gains.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GainQuantResult {
    /// First-stage codebook index (codebook domain, `0 … 7`).
    pub ga: usize,
    /// Second-stage codebook index (codebook domain, `0 … 15`).
    pub gb: usize,
    /// Quantised adaptive-codebook gain `ĝ_p` (eq (73)).
    pub gp_hat: f32,
    /// Quantised correction factor `γ̂` (eq (74) first factor).
    pub gamma_hat: f32,
    /// Quantised fixed-codebook gain `ĝ_c = γ̂·g′_c`.
    pub gc_hat: f32,
}

/// §3.9.2 conjugate-structure search: preselect 4 GA + 8 GB candidates
/// around the optimal gains, then exhaustively score the 32 pairs with
/// eq (63) (via the decode-side eq (73)/(74) reconstruction).
///
/// `g_c_prime` is the §3.9.1 predicted gain `g′_c` for this subframe.
///
/// `tame` is the taming flag ([`crate::taming::Taming::test`]): when
/// raised, every candidate pair whose reconstructed `ĝ_p` exceeds the
/// taming ceiling [`crate::taming::GPCLIP`] is excluded from the
/// search, so the quantiser returns the best *legal* pair with
/// `ĝ_p ≤ GPCLIP` (doc `docs/audio/g729/taming-procedure.md` §4). If
/// the 4 × 8 preselected cluster holds no legal pair, the search falls
/// back to the full 8 × 16 grid under the same restriction (both
/// codebooks carry low-`ĝ_p` rows, so a legal pair always exists).
#[must_use]
pub fn quantize_gains(terms: &GainTerms, g_c_prime: f32, tame: bool) -> GainQuantResult {
    let (gp_opt, gc_opt) = terms.optimal();
    // The GA preselection target: the optimal correction factor.
    let gamma_opt = if g_c_prime.abs() > 1e-9 {
        gc_opt / g_c_prime
    } else {
        gc_opt
    };

    // Preselect the GA_CANDIDATES GA rows whose γ half is closest to
    // γ_opt.
    let mut ga_rank: Vec<usize> = (0..NCODE1).collect();
    ga_rank.sort_by(|&a, &b| {
        let da = (f32::from(gain_ga_entry(a)[GAIN_VQ_COL_GC]) / Q12 - gamma_opt).abs();
        let db = (f32::from(gain_ga_entry(b)[GAIN_VQ_COL_GC]) / Q12 - gamma_opt).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    });
    // Preselect the GB_CANDIDATES GB rows whose g_p half is closest to
    // g_p_opt.
    let mut gb_rank: Vec<usize> = (0..NCODE2).collect();
    gb_rank.sort_by(|&a, &b| {
        let da = (f32::from(gain_gb_entry(a)[GAIN_VQ_COL_GP]) / Q14 - gp_opt).abs();
        let db = (f32::from(gain_gb_entry(b)[GAIN_VQ_COL_GP]) / Q14 - gp_opt).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    });

    let score = |ga_set: &[usize], gb_set: &[usize]| -> Option<(f32, GainQuantResult)> {
        let mut best: Option<(f32, GainQuantResult)> = None;
        for &ga in ga_set {
            for &gb in gb_set {
                let q = reconstruct_gains(ga, gb).expect("indices in codebook range");
                let gp_hat = q.g_p_hat;
                // Taming bound (doc §4): while tameflag is raised the
                // reconstructed pitch gain may not exceed GPCLIP.
                if tame && gp_hat > crate::taming::GPCLIP {
                    continue;
                }
                let gamma_hat = q.gamma_hat;
                let gc_hat = gamma_hat * g_c_prime;
                let e = terms.error(gp_hat, gc_hat);
                let better = match &best {
                    None => true,
                    Some((be, _)) => e < *be,
                };
                if better {
                    best = Some((
                        e,
                        GainQuantResult {
                            ga,
                            gb,
                            gp_hat,
                            gamma_hat,
                            gc_hat,
                        },
                    ));
                }
            }
        }
        best
    };

    if let Some((_, r)) = score(
        &ga_rank[..GA_CANDIDATES],
        &gb_rank[..GB_CANDIDATES.min(gb_rank.len())],
    ) {
        return r;
    }
    // Taming excluded every preselected pair: search the full grid
    // under the same ĝ_p ≤ GPCLIP restriction.
    score(&ga_rank, &gb_rank)
        .expect("codebooks carry low-gain rows, a taming-legal pair exists")
        .1
}

#[cfg(test)]
mod tests {
    use super::*;

    /// eq (63) at the normal-equations optimum is a true minimum:
    /// perturbing either gain increases the error.
    #[test]
    fn optimal_is_minimum() {
        let x: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 13 % 37) as f32) - 18.0);
        let y: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 7 % 29) as f32) - 14.0);
        let z: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 11 % 23) as f32) - 11.0);
        let t = GainTerms::compute(&x, &y, &z);
        let (gp, gc) = t.optimal();
        let e0 = t.error(gp, gc);
        for (dp, dc) in [(0.01, 0.0), (-0.01, 0.0), (0.0, 0.01), (0.0, -0.01)] {
            assert!(
                t.error(gp + dp, gc + dc) >= e0 - 1e-3,
                "perturbation ({dp},{dc}) decreased the error"
            );
        }
    }

    /// A target that is exactly `g_p·y + g_c·z` gives (near-)zero
    /// minimal error and recovers the constructing gains.
    #[test]
    fn optimal_recovers_construction() {
        let y: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 5 % 31) as f32) - 15.0);
        let z: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 17 % 41) as f32) - 20.0);
        let (gp_true, gc_true) = (0.7f32, 2.3f32);
        let x: [f32; L_SUBFR] = std::array::from_fn(|n| gp_true * y[n] + gc_true * z[n]);
        let t = GainTerms::compute(&x, &y, &z);
        let (gp, gc) = t.optimal();
        assert!((gp - gp_true).abs() < 1e-3, "gp={gp}");
        assert!((gc - gc_true).abs() < 1e-3, "gc={gc}");
        assert!(t.error(gp, gc).abs() < 1e-2);
    }

    /// The search result equals the brute-force minimum over the full
    /// 8 × 16 grid whenever the optimum lies inside the preselected
    /// clusters — verified by constructing targets *from* codebook
    /// entries (γ̂ dominated by the GA half, ĝ_p by the GB half, the
    /// conjugate bias the preselection exploits).
    #[test]
    fn search_matches_brute_force_on_codebook_targets() {
        let y: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 3 % 17) as f32) - 8.0);
        let z: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 19 % 43) as f32) - 21.0);
        let g_c_prime = 1.5f32;

        for (ga0, gb0) in [(0usize, 0usize), (3, 7), (5, 12), (7, 15), (2, 9)] {
            let q0 = reconstruct_gains(ga0, gb0).unwrap();
            let x: [f32; L_SUBFR] =
                std::array::from_fn(|n| q0.g_p_hat * y[n] + q0.gamma_hat * g_c_prime * z[n]);
            let t = GainTerms::compute(&x, &y, &z);
            let got = quantize_gains(&t, g_c_prime, false);

            // Brute force over all 128 pairs.
            let mut best_e = f32::INFINITY;
            for ga in 0..NCODE1 {
                for gb in 0..NCODE2 {
                    let q = reconstruct_gains(ga, gb).unwrap();
                    let e = t.error(q.g_p_hat, q.gamma_hat * g_c_prime);
                    if e < best_e {
                        best_e = e;
                    }
                }
            }
            let got_e = t.error(got.gp_hat, got.gc_hat);
            assert!(
                got_e <= best_e + 1e-3 * (1.0 + best_e.abs()),
                "({ga0},{gb0}): preselected search error {got_e} > brute force {best_e}"
            );
        }
    }

    /// The reconstructed gains in the result match the decode-side
    /// reconstruction of the returned indices exactly.
    #[test]
    fn result_is_decoder_consistent() {
        let x: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 23 % 53) as f32) - 26.0);
        let y: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 29 % 59) as f32) - 29.0);
        let z: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 31 % 61) as f32) - 30.0);
        let t = GainTerms::compute(&x, &y, &z);
        let res = quantize_gains(&t, 2.0, false);
        let q = reconstruct_gains(res.ga, res.gb).unwrap();
        assert_eq!(res.gp_hat, q.g_p_hat);
        assert_eq!(res.gamma_hat, q.gamma_hat);
        assert_eq!(res.gc_hat, q.gamma_hat * 2.0);
    }

    /// With the taming flag raised, the returned pair never
    /// reconstructs `ĝ_p` above the GPCLIP ceiling — even when the
    /// unrestricted optimum wants a high pitch gain — and the result
    /// is the best pair of the restricted grid.
    #[test]
    fn taming_bounds_the_pitch_gain() {
        use crate::taming::GPCLIP;
        let y: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 5 % 31) as f32) - 15.0);
        let z: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 17 % 41) as f32) - 20.0);
        // A target built with gp well above the ceiling.
        let x: [f32; L_SUBFR] = std::array::from_fn(|n| 1.15 * y[n] + 0.8 * z[n]);
        let t = GainTerms::compute(&x, &y, &z);

        let free = quantize_gains(&t, 1.0, false);
        assert!(
            free.gp_hat > GPCLIP,
            "untamed search should chase the high gain, got {}",
            free.gp_hat
        );
        let tamed = quantize_gains(&t, 1.0, true);
        assert!(
            tamed.gp_hat <= GPCLIP,
            "tamed ĝ_p {} exceeds GPCLIP",
            tamed.gp_hat
        );
        // Best-of-restricted-grid: no legal pair scores strictly better.
        let mut best_e = f32::INFINITY;
        for ga in 0..NCODE1 {
            for gb in 0..NCODE2 {
                let q = reconstruct_gains(ga, gb).unwrap();
                if q.g_p_hat > GPCLIP {
                    continue;
                }
                best_e = best_e.min(t.error(q.g_p_hat, q.gamma_hat));
            }
        }
        let got_e = t.error(tamed.gp_hat, tamed.gc_hat);
        assert!(
            got_e <= best_e + 1e-3 * (1.0 + best_e.abs()),
            "tamed search error {got_e} > restricted brute force {best_e}"
        );
    }

    /// The selected indices survive the §3.9.3 transmission mapping
    /// round-trip.
    #[test]
    fn indices_map_roundtrip() {
        use crate::gain_index_map::{demap_ga, demap_gb, map_ga, map_gb};
        let x: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 37 % 67) as f32) - 33.0);
        let y: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 41 % 71) as f32) - 35.0);
        let z: [f32; L_SUBFR] = std::array::from_fn(|n| ((n * 43 % 73) as f32) - 36.0);
        let t = GainTerms::compute(&x, &y, &z);
        let res = quantize_gains(&t, 1.0, false);
        let tx_ga = map_ga(res.ga).unwrap();
        let tx_gb = map_gb(res.gb).unwrap();
        assert_eq!(demap_ga(tx_ga).unwrap(), res.ga);
        assert_eq!(demap_gb(tx_gb).unwrap(), res.gb);
    }
}
