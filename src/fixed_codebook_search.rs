//! §3.8.1 **fixed (algebraic) codebook search** + the §3.7.3
//! adaptive-codebook gain and §3.8 target/impulse-response
//! pre-computations it depends on.
//!
//! Tenth encoder-side stage. Consumes the §3.6 target `x(n)`, the §3.5
//! impulse response `h(n)`, the §3.7 adaptive-codebook vector `v(n)`,
//! and emits the four Table-7 pulse positions/signs (directly
//! encodable via [`crate::fixed_codebook::encode_positions`] /
//! [`crate::fixed_codebook::encode_signs`]) plus the filtered
//! codevector needed by the §3.9 gain quantisation.
//!
//! ## Spec source — clauses 3.7.3, 3.8, 3.8.1, equations (43)–(60)
//!
//! (Transcribed from the EPUB prose + rasters `eq43 … eq60.jpg`.)
//!
//! * **eq (43)** — the adaptive-codebook gain
//!   `g_p = Σ x(n)·y(n) / Σ y(n)·y(n)`, bounded `0 ≤ g_p ≤ 1.2`, where
//!   `y(n)` is the filtered adaptive-codebook vector (eq (44):
//!   `y(n) = Σ_{i=0}^{n} v(i)·h(n−i)`).
//! * **eq (50)** — the fixed-codebook target:
//!   `x′(n) = x(n) − g_p·y(n)`.
//! * **eq (49)** — the harmonic pre-filter folded into the impulse
//!   response: `h(n) := h(n) + β·h(n−T)` for `n = T … 39` (`T` =
//!   current subframe's integer pitch delay, `β` = previous subframe's
//!   quantised adaptive gain clamped to `[0.2, 0.8]`, eq (47) — the
//!   clamp already exists as [`crate::pitch_sharpen::clamp_beta`]).
//! * **eq (51)** — the correlation matrix
//!   `φ(i,j) = Σ_{n=j}^{39} h(n−i)·h(n−j)` (`i ≤ j`).
//! * **eq (52)** — the backward-filtered target
//!   `d(n) = Σ_{i=n}^{39} x′(i)·h(i−n)`.
//! * **eq (53)** — the search criterion `C_k²/E_k`.
//! * **eqs (56)/(57)** — sign-folding:
//!   `φ′(i,j) = sign[d(i)]·sign[d(j)]·φ(i,j)` (`i < j`),
//!   `φ′(i,i) = 0.5·φ(i,i)` (pulse amplitudes are preset to
//!   `sign[d(m_i)]`).
//! * **eq (58)** — numerator `C = Σ |d(m_i)|`.
//! * **eq (59)** — energy `E/2 = Σ_i φ′(m_i,m_i) + Σ_{i<j} φ′(m_i,m_j)`
//!   (expanded form in the spec raster).
//! * **eq (60)** — the focused-search threshold
//!   `thr₃ = av₃ + K₃·(max₃ − av₃)`, `K₃ = 0.4`: the fourth-pulse loop
//!   is entered only when the three-pulse correlation exceeds `thr₃`,
//!   and at most `180` times per frame (two subframes).
//!
//! `max₃` / `av₃` are the maximum / average of the three-pulse
//! correlation `|d(m₀)|+|d(m₁)|+|d(m₂)|` over the track grid — both
//! separable into per-track maxima / means, which is how they are
//! computed here.

use crate::fixed_codebook::{NUM_PULSES, SUBFRAME_SIZE, TRACK_STRIDE};

/// eq (60) focused-search factor `K₃`.
pub const K3: f32 = 0.4;
/// Maximum number of fourth-loop entries per frame (clause 3.8.1).
pub const MAX_LOOP4_PER_FRAME: u32 = 180;
/// eq (43) upper bound on the adaptive-codebook gain.
pub const GP_MAX: f32 = 1.2;

/// Track-0/1/2 position grids (Table 7): `track·1 + 5·j`, 8 positions.
/// Track 3 also admits the `+1`-shifted grid (jx = 1), 16 positions.
fn track_positions(track: usize, jx: usize) -> impl Iterator<Item = usize> {
    (0..8).map(move |j| track + jx + TRACK_STRIDE * j)
}

/// eq (44): filters a real-valued vector through the impulse response
/// (zero-state lower-triangular convolution).
#[must_use]
pub fn filter_through_h(
    v: &[f32; SUBFRAME_SIZE],
    h: &[f32; SUBFRAME_SIZE],
) -> [f32; SUBFRAME_SIZE] {
    let mut y = [0.0f32; SUBFRAME_SIZE];
    for n in 0..SUBFRAME_SIZE {
        let mut acc = 0.0f32;
        for i in 0..=n {
            acc += v[i] * h[n - i];
        }
        y[n] = acc;
    }
    y
}

/// eq (43): the adaptive-codebook gain from the target `x` and the
/// filtered adaptive vector `y`, bounded to `[0, 1.2]`.
#[must_use]
pub fn adaptive_gain(x: &[f32; SUBFRAME_SIZE], y: &[f32; SUBFRAME_SIZE]) -> f32 {
    let mut num = 0.0f32;
    let mut den = 0.0f32;
    for n in 0..SUBFRAME_SIZE {
        num += x[n] * y[n];
        den += y[n] * y[n];
    }
    if den <= 0.0 {
        return 0.0;
    }
    (num / den).clamp(0.0, GP_MAX)
}

/// eq (50): the fixed-codebook target `x′(n) = x(n) − g_p·y(n)`.
#[must_use]
pub fn update_target(
    x: &[f32; SUBFRAME_SIZE],
    y: &[f32; SUBFRAME_SIZE],
    gp: f32,
) -> [f32; SUBFRAME_SIZE] {
    std::array::from_fn(|n| x[n] - gp * y[n])
}

/// eq (49): folds the §3.8 harmonic pre-filter `P(z) = 1/(1 − β·z⁻ᵀ)`
/// into the impulse response: `h(n) += β·h(n−T)` for `n ≥ T`, applied
/// only when `T < 40`. The sweep is forward so the geometric echo
/// propagates, mirroring the decode-side eq (48) codevector recursion.
pub fn prefilter_impulse_response(h: &mut [f32; SUBFRAME_SIZE], t: usize, beta: f32) {
    if t >= SUBFRAME_SIZE {
        return;
    }
    for n in t..SUBFRAME_SIZE {
        h[n] += beta * h[n - t];
    }
}

/// eq (52): the backward-filtered target
/// `d(n) = Σ_{i=n}^{39} x′(i)·h(i−n)`.
#[must_use]
pub fn correlation_d(
    x_prime: &[f32; SUBFRAME_SIZE],
    h: &[f32; SUBFRAME_SIZE],
) -> [f32; SUBFRAME_SIZE] {
    let mut d = [0.0f32; SUBFRAME_SIZE];
    for n in 0..SUBFRAME_SIZE {
        let mut acc = 0.0f32;
        for i in n..SUBFRAME_SIZE {
            acc += x_prime[i] * h[i - n];
        }
        d[n] = acc;
    }
    d
}

/// eq (51): the symmetric correlation matrix
/// `φ(i,j) = Σ_{n=j}^{39} h(n−i)·h(n−j)` for `i ≤ j` (mirrored for
/// `i > j`).
#[must_use]
pub fn phi_matrix(h: &[f32; SUBFRAME_SIZE]) -> Vec<[f32; SUBFRAME_SIZE]> {
    let mut phi = vec![[0.0f32; SUBFRAME_SIZE]; SUBFRAME_SIZE];
    for i in 0..SUBFRAME_SIZE {
        for j in i..SUBFRAME_SIZE {
            let mut acc = 0.0f32;
            for n in j..SUBFRAME_SIZE {
                acc += h[n - i] * h[n - j];
            }
            phi[i][j] = acc;
            phi[j][i] = acc;
        }
    }
    phi
}

/// One selected algebraic codevector: Table-7 pulse positions (track
/// order `i₀ … i₃`) and their ±1 signs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FixedCodebookChoice {
    /// Pulse positions `m₀ … m₃` (`m_i` on track `i`; `m₃` on either
    /// track-3 grid).
    pub positions: [u8; NUM_PULSES],
    /// Pulse signs `s₀ … s₃` (`+1` / `−1`).
    pub signs: [i8; NUM_PULSES],
}

impl FixedCodebookChoice {
    /// Builds the eq (45) ±1-pulse codevector.
    #[must_use]
    pub fn codevector(&self) -> [f32; SUBFRAME_SIZE] {
        let mut c = [0.0f32; SUBFRAME_SIZE];
        for (pos, sign) in self.positions.iter().zip(self.signs.iter()) {
            c[*pos as usize] += f32::from(*sign);
        }
        c
    }
}

/// §3.8.1 focused search. `budget` is the remaining fourth-loop budget
/// for the current frame (start each frame at
/// [`MAX_LOOP4_PER_FRAME`]; the consumed count is subtracted in
/// place so the second subframe sees what the first left over).
///
/// Returns the best `(C²/E)` pulse combination found. The search
/// always evaluates at least the first candidate reached when the
/// budget is exhausted before any fourth-loop entry, so a valid
/// codevector is always returned.
#[must_use]
pub fn search_fixed_codebook(
    d: &[f32; SUBFRAME_SIZE],
    phi: &[[f32; SUBFRAME_SIZE]],
    budget: &mut u32,
) -> FixedCodebookChoice {
    // Sign-fold: |d(n)| and sign[d(n)] (eq (56) precomputation). A zero
    // d(n) takes sign +1 (the amplitude preset must still be ±1).
    let mut d_abs = [0.0f32; SUBFRAME_SIZE];
    let mut d_sign = [1.0f32; SUBFRAME_SIZE];
    for n in 0..SUBFRAME_SIZE {
        d_abs[n] = d[n].abs();
        d_sign[n] = if d[n] < 0.0 { -1.0 } else { 1.0 };
    }
    // φ′ per eqs (56)/(57).
    let mut phi_p = vec![[0.0f32; SUBFRAME_SIZE]; SUBFRAME_SIZE];
    for i in 0..SUBFRAME_SIZE {
        for j in 0..SUBFRAME_SIZE {
            phi_p[i][j] = if i == j {
                0.5 * phi[i][i]
            } else {
                d_sign[i] * d_sign[j] * phi[i][j]
            };
        }
    }

    // eq (60): max₃ / av₃ of the three-pulse correlation, separable per
    // track.
    let mut max3 = 0.0f32;
    let mut av3 = 0.0f32;
    for track in 0..3 {
        let mut mx = 0.0f32;
        let mut sum = 0.0f32;
        for m in track_positions(track, 0) {
            mx = mx.max(d_abs[m]);
            sum += d_abs[m];
        }
        max3 += mx;
        av3 += sum / 8.0;
    }
    let thr3 = av3 + K3 * (max3 - av3);

    let mut best_score = f32::NEG_INFINITY;
    let mut best = FixedCodebookChoice {
        positions: [0, 1, 2, 3],
        signs: [1, 1, 1, 1],
    };

    for m0 in track_positions(0, 0) {
        for m1 in track_positions(1, 0) {
            // Partial energy of the first two pulses.
            let e2 = phi_p[m0][m0] + phi_p[m1][m1] + phi_p[m0][m1];
            let c2 = d_abs[m0] + d_abs[m1];
            for m2 in track_positions(2, 0) {
                let c3 = c2 + d_abs[m2];
                // eq (60) focused-search gate + frame budget.
                if c3 <= thr3 || *budget == 0 {
                    continue;
                }
                *budget -= 1;
                let e3 = e2 + phi_p[m2][m2] + phi_p[m0][m2] + phi_p[m1][m2];
                for jx in 0..2 {
                    for m3 in track_positions(3, jx) {
                        if m3 >= SUBFRAME_SIZE {
                            continue;
                        }
                        // eq (58) / (59).
                        let c = c3 + d_abs[m3];
                        let e_half =
                            e3 + phi_p[m3][m3] + phi_p[m0][m3] + phi_p[m1][m3] + phi_p[m2][m3];
                        if e_half <= 0.0 {
                            continue;
                        }
                        let score = c * c / e_half;
                        if score > best_score {
                            best_score = score;
                            best = FixedCodebookChoice {
                                positions: [m0 as u8, m1 as u8, m2 as u8, m3 as u8],
                                signs: [
                                    d_sign[m0] as i8,
                                    d_sign[m1] as i8,
                                    d_sign[m2] as i8,
                                    d_sign[m3] as i8,
                                ],
                            };
                        }
                    }
                }
            }
        }
    }

    // Fallback: if the gate rejected everything (all-zero d or an
    // exhausted budget), evaluate the plain per-track maxima so a
    // deterministic valid codevector is still produced.
    if best_score == f32::NEG_INFINITY {
        let mut positions = [0u8; NUM_PULSES];
        for (track, slot) in positions.iter_mut().enumerate().take(3) {
            let m = track_positions(track, 0)
                .max_by(|&a, &b| d_abs[a].partial_cmp(&d_abs[b]).unwrap())
                .unwrap();
            *slot = m as u8;
        }
        let m3 = (0..2)
            .flat_map(|jx| track_positions(3, jx))
            .filter(|&m| m < SUBFRAME_SIZE)
            .max_by(|&a, &b| d_abs[a].partial_cmp(&d_abs[b]).unwrap())
            .unwrap();
        positions[3] = m3 as u8;
        best = FixedCodebookChoice {
            positions,
            signs: std::array::from_fn(|i| d_sign[positions[i] as usize] as i8),
        };
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_codebook::{encode_positions, encode_signs};

    fn decaying_h() -> [f32; SUBFRAME_SIZE] {
        std::array::from_fn(|n| if n == 0 { 1.0 } else { 0.6f32.powi(n as i32) })
    }

    /// eq (43): a target equal to the filtered vector gives gain 1; a
    /// doubled target clamps at 1.2; an anti-correlated target clamps
    /// at 0.
    #[test]
    fn adaptive_gain_basics() {
        let h = decaying_h();
        let v: [f32; SUBFRAME_SIZE] = std::array::from_fn(|n| ((n * 7 % 13) as f32) - 6.0);
        let y = filter_through_h(&v, &h);
        assert!((adaptive_gain(&y, &y) - 1.0).abs() < 1e-5);
        let x2: [f32; SUBFRAME_SIZE] = std::array::from_fn(|n| 2.0 * y[n]);
        assert!((adaptive_gain(&x2, &y) - GP_MAX).abs() < 1e-6);
        let xneg: [f32; SUBFRAME_SIZE] = std::array::from_fn(|n| -y[n]);
        assert_eq!(adaptive_gain(&xneg, &y), 0.0);
    }

    /// eq (51)/(52) identities: `d(n)` is the correlation of `x′` with
    /// the shifted `h`, and `φ(i,j)` equals `Σ h(n−i)h(n−j)` — both
    /// checked against a literal Toeplitz-matrix construction.
    #[test]
    fn phi_and_d_match_toeplitz_construction() {
        let h = decaying_h();
        let x: [f32; SUBFRAME_SIZE] = std::array::from_fn(|n| ((n * 11 % 29) as f32) - 14.0);
        // H columns: column j is h delayed by j.
        let col = |j: usize| -> [f32; SUBFRAME_SIZE] {
            std::array::from_fn(|n| if n >= j { h[n - j] } else { 0.0 })
        };
        let d = correlation_d(&x, &h);
        let phi = phi_matrix(&h);
        for i in (0..SUBFRAME_SIZE).step_by(7) {
            let ci = col(i);
            let want_d: f32 = (0..SUBFRAME_SIZE).map(|n| x[n] * ci[n]).sum();
            assert!(
                (d[i] - want_d).abs() < 1e-3 * (1.0 + want_d.abs()),
                "d[{i}]={} want {want_d}",
                d[i]
            );
            for j in (i..SUBFRAME_SIZE).step_by(9) {
                let cj = col(j);
                let want_phi: f32 = (0..SUBFRAME_SIZE).map(|n| ci[n] * cj[n]).sum();
                assert!(
                    (phi[i][j] - want_phi).abs() < 1e-3 * (1.0 + want_phi.abs()),
                    "phi[{i}][{j}]={} want {want_phi}",
                    phi[i][j]
                );
            }
        }
    }

    /// eq (49): for `T ≥ 40` the response is untouched; for `T < 40`
    /// the tail gains the β-scaled echo.
    #[test]
    fn prefilter_impulse_response_behaviour() {
        let mut h = decaying_h();
        let orig = h;
        prefilter_impulse_response(&mut h, 40, 0.8);
        assert_eq!(h, orig);
        prefilter_impulse_response(&mut h, 25, 0.5);
        for n in 0..25 {
            assert_eq!(h[n], orig[n]);
        }
        for n in 25..SUBFRAME_SIZE {
            assert!((h[n] - (orig[n] + 0.5 * h[n - 25])).abs() < 1e-6);
        }
    }

    /// The search recovers an exact codevector: build `x′` as the
    /// filtered image of a known Table-7 pulse set — the search must
    /// return exactly those positions and signs (the criterion is
    /// maximised by the true generator).
    #[test]
    fn search_recovers_planted_pulses() {
        let h = decaying_h();
        let truth = FixedCodebookChoice {
            positions: [10, 26, 7, 33],
            signs: [1, -1, 1, -1],
        };
        let c = truth.codevector();
        let x = filter_through_h(&c, &h);
        let d = correlation_d(&x, &h);
        let phi = phi_matrix(&h);
        let mut budget = MAX_LOOP4_PER_FRAME;
        let got = search_fixed_codebook(&d, &phi, &mut budget);
        assert_eq!(got.positions, truth.positions, "positions");
        assert_eq!(got.signs, truth.signs, "signs");
        assert!(budget < MAX_LOOP4_PER_FRAME, "loop-4 must have run");
    }

    /// Every search result is Table-7-encodable via the existing
    /// transmission-side mappings.
    #[test]
    fn results_are_encodable() {
        let h = decaying_h();
        for seed in 0..6 {
            let x: [f32; SUBFRAME_SIZE] =
                std::array::from_fn(|n| (((n * 17 + seed * 31) % 61) as f32) - 30.0);
            let d = correlation_d(&x, &h);
            let phi = phi_matrix(&h);
            let mut budget = MAX_LOOP4_PER_FRAME;
            let got = search_fixed_codebook(&d, &phi, &mut budget);
            assert!(
                encode_positions(&got.positions).is_some(),
                "positions {:?} not encodable",
                got.positions
            );
            assert!(
                encode_signs(&got.signs).is_some(),
                "signs {:?} not encodable",
                got.signs
            );
        }
    }

    /// An exhausted budget still yields a deterministic valid choice
    /// (the per-track-maxima fallback).
    #[test]
    fn zero_budget_fallback() {
        let h = decaying_h();
        let x: [f32; SUBFRAME_SIZE] = std::array::from_fn(|n| ((n * 3 % 19) as f32) - 9.0);
        let d = correlation_d(&x, &h);
        let phi = phi_matrix(&h);
        let mut budget = 0u32;
        let got = search_fixed_codebook(&d, &phi, &mut budget);
        assert!(encode_positions(&got.positions).is_some());
        assert!(encode_signs(&got.signs).is_some());
        assert_eq!(budget, 0);
    }
}
