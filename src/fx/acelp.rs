//! §3.8.1 **algebraic codebook search** on the fixed grid — the
//! eq (51) correlation matrix, the eq (52) backward-filtered target,
//! the eq (56)/(57) sign folding, the eq (60) focused-search threshold
//! and the `C²/E` maximisation over the Table-7 pulse grid.
//!
//! ## Number grids
//!
//! - `d(n)` (eq (52)) — the Q0 × Q12 sums run on an exact wide
//!   accumulator and the whole vector is then normalised so that
//!   `max |d(n)| < 2^13`: four pulses sum inside Word16 (eq (58)).
//! - `φ(i, j)` (eq (51)) — exact wide sums of Q12 × Q12 products,
//!   normalised so that `φ(0, 0)` (the largest entry) sits below
//!   `2^12`: the eq (59) ten-term energy sum stays inside Word16.
//! - The eq (53) comparison `C²/E` is a cross-multiplication on the
//!   Word32 grid: `C² · E_best` against `C_best² · E`, with `C²` as a
//!   Word32 product and the products against `E` on a wide accumulator
//!   ([`AcelpLatitude::compare_wide`]) or the Word16 high word of `C²`
//!   times `E` on Word32.
//! - `thr₃` (eq (60)) — `av₃ + 0.4·(max₃ − av₃)` on the `d` grid with
//!   the Q15 literal `13107`; `av₃` is the 24-position sum shifted
//!   right by 3.
//!
//! The prose pins the equations, the threshold form, `K₃ = 0.4`, the
//! 180-entry frame budget and the loop order; the scaling of `d`/`φ`
//! and the comparison precision are the fixed chain's own composition,
//! exposed as latitude and pinned black-box.

use crate::fixed_codebook::{NUM_PULSES, SUBFRAME_SIZE, TRACK_STRIDE};
use crate::fx::filters::L_SUBFR;
use crate::fx::ops::{add, extract_h, l_mult, mult, sub};

/// Maximum number of fourth-loop entries per frame (clause 3.8.1).
pub const MAX_LOOP4_PER_FRAME: u32 = 180;
/// eq (60) `K₃` on Q15.
const K3_Q15: i16 = 13107;

/// Unstated latitude of the search arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcelpLatitude {
    /// Compare `C²·E` products on an exact wide accumulator (`true`)
    /// or on the Word32 grid through the Word16 high word of `C²`.
    pub compare_wide: bool,
    /// Normalisation headroom of `d`: `max |d| < 2^(15 − d_headroom)`.
    pub d_headroom: i16,
    /// Normalisation headroom of `φ`: `φ(0,0) < 2^(15 − phi_headroom)`.
    pub phi_headroom: i16,
    /// Evaluate the search on the exact (un-normalised) `d`/`φ` in
    /// double precision (the spec-equation oracle).
    pub exact: bool,
}

impl Default for AcelpLatitude {
    fn default() -> Self {
        Self {
            compare_wide: true,
            d_headroom: 2,
            phi_headroom: 3,
            exact: false,
        }
    }
}

/// One selected codevector: Table-7 pulse positions and signs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcelpChoice {
    /// Pulse positions `m₀ … m₃`.
    pub positions: [u8; NUM_PULSES],
    /// Pulse signs `s₀ … s₃` (`±1`).
    pub signs: [i8; NUM_PULSES],
}

/// Track-0/1/2 position grids (Table 7): `track + 5·j`; track 3 also
/// admits the `+1`-shifted grid (`jx = 1`).
fn track_positions(track: usize, jx: usize) -> impl Iterator<Item = usize> {
    (0..8).map(move |j| track + jx + TRACK_STRIDE * j)
}

/// eq (52) on an exact accumulator (Q0 × Q12, doubled).
fn correlation_d_wide(x_prime: &[i16; L_SUBFR], h_q12: &[i16; L_SUBFR]) -> [i64; L_SUBFR] {
    let mut d = [0i64; L_SUBFR];
    for n in 0..L_SUBFR {
        let mut acc = 0i64;
        for i in n..L_SUBFR {
            acc += 2 * i64::from(x_prime[i]) * i64::from(h_q12[i - n]);
        }
        d[n] = acc;
    }
    d
}

/// eq (51) upper triangle on an exact accumulator (Q12 × Q12,
/// doubled), mirrored.
fn phi_wide(h_q12: &[i16; L_SUBFR]) -> Vec<[i64; L_SUBFR]> {
    let mut phi = vec![[0i64; L_SUBFR]; L_SUBFR];
    for i in 0..L_SUBFR {
        for j in i..L_SUBFR {
            let mut acc = 0i64;
            for n in j..L_SUBFR {
                acc += 2 * i64::from(h_q12[n - i]) * i64::from(h_q12[n - j]);
            }
            phi[i][j] = acc;
            phi[j][i] = acc;
        }
    }
    phi
}

/// Normalises a wide value set onto Word16 with the given headroom:
/// returns the right-shift such that `max |v| < 2^(15 − headroom)`.
fn wide_shift(max_abs: i64, headroom: i16) -> i32 {
    if max_abs == 0 {
        return 0;
    }
    let bits = 64 - max_abs.leading_zeros() as i32; // position of the top bit
    (bits - (15 - i32::from(headroom))).max(0)
}

fn shr_wide(v: i64, s: i32) -> i16 {
    let r = if s >= 63 { 0 } else { v >> s };
    r.clamp(-32768, 32767) as i16
}

/// §3.8.1 focused search over one subframe under the default latitude.
/// `budget` is the frame's remaining fourth-loop budget (consumed in
/// place).
#[must_use]
pub fn search_acelp_fx(
    x_prime: &[i16; L_SUBFR],
    h_pre_q12: &[i16; L_SUBFR],
    budget: &mut u32,
) -> AcelpChoice {
    search_acelp_fx_lat(x_prime, h_pre_q12, budget, &AcelpLatitude::default())
}

/// §3.8.1 focused search under an explicit latitude.
#[must_use]
pub fn search_acelp_fx_lat(
    x_prime: &[i16; L_SUBFR],
    h_pre_q12: &[i16; L_SUBFR],
    budget: &mut u32,
    lat: &AcelpLatitude,
) -> AcelpChoice {
    let d_wide = correlation_d_wide(x_prime, h_pre_q12);
    let phi_w = phi_wide(h_pre_q12);

    if lat.exact {
        return search_exact(&d_wide, &phi_w, budget);
    }

    // Word16 grids.
    let d_max = d_wide.iter().map(|v| v.abs()).max().unwrap_or(0);
    let sd = wide_shift(d_max, lat.d_headroom);
    let d16: [i16; L_SUBFR] = std::array::from_fn(|n| shr_wide(d_wide[n], sd));
    let phi_max = phi_w[0][0];
    let sp = wide_shift(phi_max, lat.phi_headroom);

    // eq (56)/(57): sign fold, halved diagonal.
    let mut d_abs = [0i16; L_SUBFR];
    let mut sign = [1i8; L_SUBFR];
    for n in 0..L_SUBFR {
        if d16[n] < 0 {
            sign[n] = -1;
            d_abs[n] = sub(0, d16[n]);
        } else {
            d_abs[n] = d16[n];
        }
    }
    let mut phi_p = vec![[0i16; L_SUBFR]; L_SUBFR];
    for i in 0..L_SUBFR {
        for j in 0..L_SUBFR {
            let v = shr_wide(phi_w[i][j], sp);
            phi_p[i][j] = if i == j {
                crate::fx::ops::shr(v, 1)
            } else if sign[i] == sign[j] {
                v
            } else {
                sub(0, v)
            };
        }
    }

    // eq (60): max₃ / av₃ over the three first tracks.
    let mut max3 = 0i16;
    let mut sum24 = 0i32;
    for track in 0..3 {
        let mut mx = 0i16;
        for m in track_positions(track, 0) {
            if d_abs[m] > mx {
                mx = d_abs[m];
            }
            sum24 += i32::from(d_abs[m]);
        }
        max3 = add(max3, mx);
    }
    let av3 = (sum24 >> 3) as i16;
    let thr3 = add(av3, mult(sub(max3, av3), K3_Q15));

    let mut best = AcelpChoice {
        positions: [0, 1, 2, 3],
        signs: [1, 1, 1, 1],
    };
    // Best (C, E) so far; E = 0 sentinel means "nothing yet".
    let mut best_c: i16 = 0;
    let mut best_e: i16 = 0;
    let mut found = false;

    for m0 in track_positions(0, 0) {
        for m1 in track_positions(1, 0) {
            let c2 = add(d_abs[m0], d_abs[m1]);
            let e2 = add(add(phi_p[m0][m0], phi_p[m1][m1]), phi_p[m0][m1]);
            for m2 in track_positions(2, 0) {
                let c3 = add(c2, d_abs[m2]);
                if c3 <= thr3 || *budget == 0 {
                    continue;
                }
                *budget -= 1;
                let e3 = add(add(e2, phi_p[m2][m2]), add(phi_p[m0][m2], phi_p[m1][m2]));
                for jx in 0..2 {
                    for m3 in track_positions(3, jx) {
                        if m3 >= SUBFRAME_SIZE {
                            continue;
                        }
                        let c = add(c3, d_abs[m3]);
                        let e = add(
                            add(e3, phi_p[m3][m3]),
                            add(add(phi_p[m0][m3], phi_p[m1][m3]), phi_p[m2][m3]),
                        );
                        if e <= 0 {
                            continue;
                        }
                        let better = if !found {
                            true
                        } else if lat.compare_wide {
                            let lhs = i64::from(c) * i64::from(c) * i64::from(best_e);
                            let rhs = i64::from(best_c) * i64::from(best_c) * i64::from(e);
                            lhs > rhs
                        } else {
                            // C² high word (Q31 → 16 bits) against E.
                            let c2h = extract_h(l_mult(c, c));
                            let b2h = extract_h(l_mult(best_c, best_c));
                            l_mult(c2h, best_e) > l_mult(b2h, e)
                        };
                        if better {
                            found = true;
                            best_c = c;
                            best_e = e;
                            best = AcelpChoice {
                                positions: [m0 as u8, m1 as u8, m2 as u8, m3 as u8],
                                signs: [sign[m0], sign[m1], sign[m2], sign[m3]],
                            };
                        }
                    }
                }
            }
        }
    }

    if !found {
        best = fallback(&d_abs, &sign);
    }
    best
}

/// Deterministic valid choice when the gate rejected everything: the
/// per-track maxima.
fn fallback(d_abs: &[i16; L_SUBFR], sign: &[i8; L_SUBFR]) -> AcelpChoice {
    let mut positions = [0u8; NUM_PULSES];
    for (track, slot) in positions.iter_mut().enumerate().take(3) {
        let m = track_positions(track, 0).max_by_key(|&m| d_abs[m]).unwrap();
        *slot = m as u8;
    }
    let m3 = (0..2)
        .flat_map(|jx| track_positions(3, jx))
        .filter(|&m| m < SUBFRAME_SIZE)
        .max_by_key(|&m| d_abs[m])
        .unwrap();
    positions[3] = m3 as u8;
    AcelpChoice {
        positions,
        signs: std::array::from_fn(|i| sign[positions[i] as usize]),
    }
}

/// The spec-equation oracle: the same loop on exact doubles.
fn search_exact(d: &[i64; L_SUBFR], phi: &[[i64; L_SUBFR]], budget: &mut u32) -> AcelpChoice {
    let d_abs: [f64; L_SUBFR] = std::array::from_fn(|n| (d[n] as f64).abs());
    let sign: [i8; L_SUBFR] = std::array::from_fn(|n| if d[n] < 0 { -1 } else { 1 });
    let phi_p = |i: usize, j: usize| -> f64 {
        if i == j {
            0.5 * phi[i][i] as f64
        } else {
            f64::from(sign[i]) * f64::from(sign[j]) * phi[i][j] as f64
        }
    };
    let mut max3 = 0.0f64;
    let mut av3 = 0.0f64;
    for track in 0..3 {
        let mut mx = 0.0f64;
        let mut sum = 0.0f64;
        for m in track_positions(track, 0) {
            mx = mx.max(d_abs[m]);
            sum += d_abs[m];
        }
        max3 += mx;
        av3 += sum / 8.0;
    }
    let thr3 = av3 + 0.4 * (max3 - av3);
    let mut best_score = f64::NEG_INFINITY;
    let mut best: Option<AcelpChoice> = None;
    for m0 in track_positions(0, 0) {
        for m1 in track_positions(1, 0) {
            let e2 = phi_p(m0, m0) + phi_p(m1, m1) + phi_p(m0, m1);
            let c2 = d_abs[m0] + d_abs[m1];
            for m2 in track_positions(2, 0) {
                let c3 = c2 + d_abs[m2];
                if c3 <= thr3 || *budget == 0 {
                    continue;
                }
                *budget -= 1;
                let e3 = e2 + phi_p(m2, m2) + phi_p(m0, m2) + phi_p(m1, m2);
                for jx in 0..2 {
                    for m3 in track_positions(3, jx) {
                        if m3 >= SUBFRAME_SIZE {
                            continue;
                        }
                        let c = c3 + d_abs[m3];
                        let e = e3 + phi_p(m3, m3) + phi_p(m0, m3) + phi_p(m1, m3) + phi_p(m2, m3);
                        if e <= 0.0 {
                            continue;
                        }
                        let score = c * c / e;
                        if score > best_score {
                            best_score = score;
                            best = Some(AcelpChoice {
                                positions: [m0 as u8, m1 as u8, m2 as u8, m3 as u8],
                                signs: [sign[m0], sign[m1], sign[m2], sign[m3]],
                            });
                        }
                    }
                }
            }
        }
    }
    best.unwrap_or_else(|| {
        let d16: [i16; L_SUBFR] = std::array::from_fn(|n| (d_abs[n].min(32767.0)) as i16);
        fallback(&d16, &sign)
    })
}

/// eq (53) score of a pulse set on the exact grids (diagnostic).
#[doc(hidden)]
#[must_use]
pub fn score_choice_exact(
    x_prime: &[i16; L_SUBFR],
    h_pre_q12: &[i16; L_SUBFR],
    choice: &AcelpChoice,
) -> f64 {
    let d = correlation_d_wide(x_prime, h_pre_q12);
    let phi = phi_wide(h_pre_q12);
    let mut num = 0.0f64;
    let mut en = 0.0f64;
    for i in 0..NUM_PULSES {
        let mi = usize::from(choice.positions[i]);
        let si = f64::from(choice.signs[i]);
        num += si * d[mi] as f64;
        for j in 0..NUM_PULSES {
            let mj = usize::from(choice.positions[j]);
            let sj = f64::from(choice.signs[j]);
            en += si * sj * phi[mi][mj] as f64;
        }
    }
    if en <= 0.0 {
        0.0
    } else {
        num * num / en
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_codebook::{encode_positions, encode_signs};

    fn decaying_h() -> [i16; L_SUBFR] {
        std::array::from_fn(|n| (4096.0 * 0.6f64.powi(n as i32)) as i16)
    }

    /// Planted pulses are recovered exactly (the criterion is
    /// maximised by the generator).
    #[test]
    fn recovers_planted_pulses() {
        let h = decaying_h();
        let truth = AcelpChoice {
            positions: [10, 26, 7, 33],
            signs: [1, -1, 1, -1],
        };
        let mut c = [0i16; L_SUBFR];
        for k in 0..4 {
            c[truth.positions[k] as usize] = if truth.signs[k] > 0 { 8192 } else { -8192 };
        }
        let z = crate::fx::filters::convolve_code_q12(&c, &h);
        // A Q0 target: the Q12 z scaled down to a speech-like level.
        let x: [i16; L_SUBFR] = std::array::from_fn(|n| z[n] / 2);
        let mut budget = MAX_LOOP4_PER_FRAME;
        let got = search_acelp_fx(&x, &h, &mut budget);
        assert_eq!(got.positions, truth.positions);
        assert_eq!(got.signs, truth.signs);
        assert!(budget < MAX_LOOP4_PER_FRAME);
    }

    /// Results are always Table-7 encodable, under both comparison
    /// modes and the exact oracle.
    #[test]
    fn results_are_encodable() {
        let h = decaying_h();
        for seed in 0..6i32 {
            let x: [i16; L_SUBFR] =
                std::array::from_fn(|n| (((n as i32 * 17 + seed * 31) % 61) * 40 - 1200) as i16);
            for lat in [
                AcelpLatitude::default(),
                AcelpLatitude {
                    compare_wide: false,
                    ..Default::default()
                },
                AcelpLatitude {
                    exact: true,
                    ..Default::default()
                },
            ] {
                let mut budget = MAX_LOOP4_PER_FRAME;
                let got = search_acelp_fx_lat(&x, &h, &mut budget, &lat);
                assert!(encode_positions(&got.positions).is_some(), "{got:?}");
                assert!(encode_signs(&got.signs).is_some());
            }
        }
    }

    /// An exhausted budget still yields a deterministic valid choice.
    #[test]
    fn zero_budget_fallback() {
        let h = decaying_h();
        let x: [i16; L_SUBFR] = std::array::from_fn(|n| ((n * 3 % 19) as i16) * 50 - 450);
        let mut budget = 0u32;
        let got = search_acelp_fx(&x, &h, &mut budget);
        assert!(encode_positions(&got.positions).is_some());
        assert_eq!(budget, 0);
    }

    /// The Word16 search agrees with the exact oracle on moderate
    /// inputs.
    #[test]
    fn word16_matches_exact_on_moderate_input() {
        let h = decaying_h();
        let mut agree = 0;
        for seed in 0..12i32 {
            let x: [i16; L_SUBFR] =
                std::array::from_fn(|n| (((n as i32 * 23 + seed * 37) % 97) * 30 - 1440) as i16);
            let mut b1 = MAX_LOOP4_PER_FRAME;
            let mut b2 = MAX_LOOP4_PER_FRAME;
            let a = search_acelp_fx(&x, &h, &mut b1);
            let e = search_acelp_fx_lat(
                &x,
                &h,
                &mut b2,
                &AcelpLatitude {
                    exact: true,
                    ..Default::default()
                },
            );
            agree += usize::from(a == e);
        }
        assert!(agree >= 10, "only {agree}/12 agree with the exact oracle");
    }
}
