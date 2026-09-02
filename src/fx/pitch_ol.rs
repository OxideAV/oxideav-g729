//! §3.4 **open-loop pitch analysis** on the fixed grid — the eq (34)
//! correlation maxima over the three delay sections, the eq (35)
//! normalisation, and the favour-lower-delays selection.
//!
//! ## Number grids
//!
//! - Weighted speech `sw(n)` — Q0 Word16. The eq (34)/(35) sums run in
//!   a Word32 `l_mac` accumulator; before the search the whole
//!   223-sample buffer's energy is checked against the Word32 range
//!   and the signal right-shifted (by [`OpenLoopLatitude::overflow_shift`])
//!   until it fits — the overflow-rescale protocol the §3.2.1
//!   autocorrelation uses. Every `R(k)` is bounded by that energy
//!   (Cauchy-Schwarz), so a single rescale covers the whole search.
//! - `R′(t_i)` — `R(t_i)` times the Q30 [`crate::fx::dsp::inv_sqrt`] of
//!   the lagged energy, kept as the Word32 high product (`mpy_32`) so
//!   the three candidates compare on a common scale; the 0.85 factor is
//!   the Q15 literal `27853` applied through `mpy_32_16`.
//!
//! The prose pins the equations and the selection order; the
//! accumulator scaling and the normalisation precision are the fixed
//! chain's own composition, exposed as latitude and pinned black-box.

use crate::fx::dsp::{inv_sqrt, l_extract, mpy_32, mpy_32_16};
use crate::fx::ops::{extract_h, l_mac, l_mult, l_shl, l_shr, l_sub, norm_l, sature32, shr};

/// Maximum open-loop delay.
pub const PIT_MAX: usize = 143;
/// Minimum open-loop delay.
pub const PIT_MIN: usize = 20;
/// Frame length.
pub const L_FRAME: usize = 80;
/// `[143 samples of history | 80 samples of the present frame]`.
pub const PIT_BUFFER: usize = PIT_MAX + L_FRAME;

/// The three delay sections, longest first (clause 3.4).
const SECTIONS: [(usize, usize); 3] = [(80, 143), (40, 79), (20, 39)];

/// `0.85` on Q15.
const FAVOUR_Q15: i16 = 27853;

/// How the eq (35) normalised correlations are formed and compared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OlNorm {
    /// Word32 `R` × Q30 `inv_sqrt(E)` through `mpy_32` (the high
    /// product only).
    Mpy32,
    /// Both `R` and `inv_sqrt(E)` normalised to Word16 mantissas with
    /// tracked exponents, multiplied on Word32 and re-aligned.
    Word16,
    /// Exact `R / √E` in double precision (the spec-equation oracle).
    Wide,
}

/// Unstated latitude of the §3.4 arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpenLoopLatitude {
    /// Right-shift applied to the weighted speech on each Word32
    /// overflow of the buffer energy before the retry; `0` disables
    /// the rescale and runs the sums on an exact wide accumulator.
    pub overflow_shift: i16,
    /// Normalisation / comparison arithmetic.
    pub norm: OlNorm,
    /// Within a section, a later delay must beat the running maximum
    /// strictly (`true`, lowest delay wins ties) or may equal it
    /// (`false`, highest delay wins ties).
    pub strict_max: bool,
    /// The favour-lower-delays test takes a shorter candidate only
    /// when its normalised correlation strictly exceeds `0.85·R′(T_op)`
    /// (`true`) or also on equality (`false`, the printed `≥`).
    pub favour_strict: bool,
}

impl Default for OpenLoopLatitude {
    fn default() -> Self {
        Self {
            overflow_shift: 3,
            norm: OlNorm::Word16,
            strict_max: true,
            favour_strict: true,
        }
    }
}

/// The selected open-loop delay.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpenLoopPitchFx {
    /// `T_op ∈ [20, 143]`.
    pub t_op: i32,
    /// Total right-shift the overflow-rescale protocol applied.
    pub signal_shift: i16,
    /// The three section maxima `(t_i, R′(t_i))` (diagnostic; the
    /// normalised value as a double, on the exact scale).
    pub candidates: [(i32, f64); 3],
}

/// eq (34) at one delay on the (already rescaled) buffer.
fn correlation(buf: &[i16; PIT_BUFFER], k: usize) -> i32 {
    let mut acc = 0i32;
    for n in 0..L_FRAME {
        acc = l_mac(acc, buf[PIT_MAX + n], buf[PIT_MAX + n - k]);
    }
    acc
}

/// eq (35) denominator sum at one delay.
fn lagged_energy(buf: &[i16; PIT_BUFFER], k: usize) -> i32 {
    let mut acc = 0i32;
    for n in 0..L_FRAME {
        let v = buf[PIT_MAX + n - k];
        acc = l_mac(acc, v, v);
    }
    acc
}

/// §3.4 over one frame of weighted speech (Q0) under the default
/// latitude.
#[must_use]
pub fn open_loop_pitch_fx(sw: &[i16; PIT_BUFFER]) -> OpenLoopPitchFx {
    open_loop_pitch_fx_lat(sw, &OpenLoopLatitude::default())
}

/// §3.4 over one frame of weighted speech (Q0).
#[must_use]
pub fn open_loop_pitch_fx_lat(sw: &[i16; PIT_BUFFER], lat: &OpenLoopLatitude) -> OpenLoopPitchFx {
    // Overflow-rescale protocol on the whole-buffer energy (the bound
    // on every correlation and lagged energy below).
    let mut buf = *sw;
    let mut signal_shift = 0i16;
    if lat.overflow_shift > 0 {
        loop {
            let exact: i64 = buf.iter().map(|&v| 2 * i64::from(v) * i64::from(v)).sum();
            if exact < i64::from(i32::MAX) {
                break;
            }
            for v in buf.iter_mut() {
                *v = shr(*v, lat.overflow_shift);
            }
            signal_shift += lat.overflow_shift;
        }
    }
    let wide = lat.overflow_shift == 0;

    // Step 1: the raw-correlation maximum in each section.
    let mut cand = [(0usize, 0i64); 3];
    for (s, &(lo, hi)) in SECTIONS.iter().enumerate() {
        let mut best_k = lo;
        let mut best_r = i64::MIN;
        for k in lo..=hi {
            let r = if wide {
                correlation_wide(&buf, k)
            } else {
                i64::from(correlation(&buf, k))
            };
            let better = if lat.strict_max {
                r > best_r
            } else {
                r >= best_r
            };
            if better {
                best_r = r;
                best_k = k;
            }
        }
        cand[s] = (best_k, best_r);
    }

    // Step 2: eq (35) normalisation.
    let energies: [i64; 3] = std::array::from_fn(|s| {
        if wide {
            lagged_energy_wide(&buf, cand[s].0)
        } else {
            i64::from(lagged_energy(&buf, cand[s].0))
        }
    });

    // Step 3: favour-lower-delays selection (t₁ → t₂ → t₃).
    let t_op = match lat.norm {
        OlNorm::Wide => {
            let norm: [f64; 3] = std::array::from_fn(|s| {
                let e = energies[s] as f64;
                if e > 0.0 {
                    cand[s].1 as f64 / e.sqrt()
                } else {
                    0.0
                }
            });
            let mut t_op = cand[0].0;
            let mut r_best = norm[0];
            for s in 1..3 {
                let take = if lat.favour_strict {
                    norm[s] > 0.85 * r_best
                } else {
                    norm[s] >= 0.85 * r_best
                };
                if take {
                    r_best = norm[s];
                    t_op = cand[s].0;
                }
            }
            t_op
        }
        OlNorm::Mpy32 => {
            let norm: [i32; 3] = std::array::from_fn(|s| {
                let r = sature32(cand[s].1);
                let inv = inv_sqrt(sature32(energies[s])); // Q30 of 1/√e
                let (r_hi, r_lo) = l_extract(r);
                let (i_hi, i_lo) = l_extract(inv);
                mpy_32(r_hi, r_lo, i_hi, i_lo)
            });
            let mut t_op = cand[0].0;
            let mut r_best = norm[0];
            for s in 1..3 {
                let (b_hi, b_lo) = l_extract(r_best);
                let thr = mpy_32_16(b_hi, b_lo, FAVOUR_Q15);
                if favour(l_sub(norm[s], thr), lat) {
                    r_best = norm[s];
                    t_op = cand[s].0;
                }
            }
            t_op
        }
        OlNorm::Word16 => {
            // Word16 mantissas with exponents: R = r16 · 2^(16 − nr),
            // 1/√E = i16 · 2^(−30 + 16 − ni); the product's exponent is
            // tracked and the three candidates aligned to the smallest.
            let mut mant = [0i32; 3];
            let mut expo = [0i16; 3];
            for s in 0..3 {
                let r = sature32(cand[s].1);
                let inv = inv_sqrt(sature32(energies[s]));
                let nr = norm_l(r);
                let ni = norm_l(inv);
                let r16 = extract_h(l_shl(r, nr));
                let i16v = extract_h(l_shl(inv, ni));
                mant[s] = l_mult(r16, i16v);
                expo[s] = -(nr + ni);
            }
            // Align to the largest exponent (right shifts only).
            let e_max = expo.iter().copied().max().unwrap_or(0);
            let norm: [i32; 3] = std::array::from_fn(|s| l_shr(mant[s], e_max - expo[s]));
            let mut t_op = cand[0].0;
            let mut r_best = norm[0];
            for s in 1..3 {
                let thr = mult_32_q15(r_best, FAVOUR_Q15);
                if favour(l_sub(norm[s], thr), lat) {
                    r_best = norm[s];
                    t_op = cand[s].0;
                }
            }
            t_op
        }
    };

    let candidates: [(i32, f64); 3] = std::array::from_fn(|s| {
        let e = energies[s] as f64;
        let v = if e > 0.0 {
            cand[s].1 as f64 / e.sqrt()
        } else {
            0.0
        };
        (cand[s].0 as i32, v)
    });
    OpenLoopPitchFx {
        t_op: t_op as i32,
        signal_shift,
        candidates,
    }
}

/// eq (34) on an exact wide accumulator.
fn correlation_wide(buf: &[i16; PIT_BUFFER], k: usize) -> i64 {
    (0..L_FRAME)
        .map(|n| 2 * i64::from(buf[PIT_MAX + n]) * i64::from(buf[PIT_MAX + n - k]))
        .sum()
}

/// eq (35) denominator on an exact wide accumulator.
fn lagged_energy_wide(buf: &[i16; PIT_BUFFER], k: usize) -> i64 {
    (0..L_FRAME)
        .map(|n| {
            let v = i64::from(buf[PIT_MAX + n - k]);
            2 * v * v
        })
        .sum()
}

/// The favour-lower-delays acceptance on the signed difference
/// `R′(t_i) − 0.85·R′(T_op)`.
fn favour(diff: i32, lat: &OpenLoopLatitude) -> bool {
    if lat.favour_strict {
        diff > 0
    } else {
        diff >= 0
    }
}

/// Word32 × Q15 keeping the Word32 scale.
fn mult_32_q15(x: i32, q15: i16) -> i32 {
    let (hi, lo) = l_extract(x);
    mpy_32_16(hi, lo, q15)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn periodic_buf(period: usize, amp: i16) -> [i16; PIT_BUFFER] {
        std::array::from_fn(|n| {
            let phase = (n % period) as f64 / period as f64;
            let v = (std::f64::consts::TAU * phase).sin()
                + 0.3 * (2.0 * std::f64::consts::TAU * phase).sin();
            (v * f64::from(amp)) as i16
        })
    }

    /// A period-30 signal is picked at 30 over its multiples; at a
    /// moderate level no rescale is needed.
    #[test]
    fn period_30_beats_multiples() {
        let r = open_loop_pitch_fx(&periodic_buf(30, 2000));
        assert_eq!(r.t_op, 30);
        assert_eq!(r.signal_shift, 0);
    }

    /// Period 55 (t₂ range) beats its double in the t₁ range.
    #[test]
    fn period_55_beats_double() {
        let r = open_loop_pitch_fx(&periodic_buf(55, 6000));
        assert_eq!(r.t_op, 55);
    }

    /// A full-scale input triggers the rescale and still resolves the
    /// period.
    #[test]
    fn full_scale_rescales() {
        let r = open_loop_pitch_fx(&periodic_buf(45, 24000));
        assert!(r.signal_shift > 0);
        assert_eq!(r.t_op, 45);
    }

    /// Silence is safe and in range.
    #[test]
    fn silence_is_safe() {
        let r = open_loop_pitch_fx(&[0; PIT_BUFFER]);
        assert!((PIT_MIN as i32..=PIT_MAX as i32).contains(&r.t_op));
    }

    /// The fixed search agrees with the float module on moderate
    /// synthetic signals.
    #[test]
    fn agrees_with_float_module() {
        for period in [23usize, 37, 61, 90, 120] {
            let buf = periodic_buf(period, 5000);
            let f: [f32; PIT_BUFFER] = std::array::from_fn(|n| f32::from(buf[n]));
            let want = crate::open_loop_pitch::open_loop_pitch(&f).t_op as i32;
            assert_eq!(open_loop_pitch_fx(&buf).t_op, want, "period {period}");
        }
    }
}
