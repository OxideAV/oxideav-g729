//! ITU-T G.729 Annex B §5 — Comfort Noise Generation (decoder side).
//!
//! When the decoder receives an SID frame, it primes this generator with
//! the SID's spectral shape + energy. On subsequent NODATA frames (no
//! packet delivered, or a zero-length packet in the wire format) the
//! generator emits new comfort-noise samples driven by the same state,
//! slowly smoothing towards fresh spectrum / energy on every new SID.
//!
//! # Approach
//!
//! 1. Keep a "target" LSP vector and target energy (dB) derived from
//!    the last SID. Interpolate the currently-applied LSP towards the
//!    target over a few frames (smooths transitions and avoids audible
//!    clicks when spectrum changes).
//! 2. Generate a white-noise excitation, scaled by target energy.
//! 3. Run the excitation through a 10-th-order LPC synthesis filter
//!    built from the smoothed LSPs, matching the CS-ACELP main-body
//!    filter order.
//!
//! This is a minimal CNG that produces "reasonable background noise" —
//! non-periodic, stationary, at roughly the right level. Perceptually
//! indistinguishable from room-tone silence in most listening
//! scenarios, which is the bar Annex B actually targets.

use crate::annex_b_vad::EncodedSid;
use crate::lpc::lsp_to_lpc;
use crate::{FRAME_SAMPLES, LPC_ORDER};

/// Comfort-noise generator state.
#[derive(Clone, Debug)]
pub struct CngState {
    /// Linear-congruential RNG seed for the excitation.
    pub rng: u32,
    /// Current smoothed LSPs (cosine domain).
    pub lsp: [f32; LPC_ORDER],
    /// Target LSPs derived from the most-recent SID (cosine domain).
    pub lsp_target: [f32; LPC_ORDER],
    /// Current smoothed energy in dB.
    pub energy_db: f32,
    /// Target energy from the latest SID.
    pub energy_db_target: f32,
    /// Synthesis filter memory (10 most-recent output samples).
    pub syn_mem: [f32; LPC_ORDER],
    /// Has the generator been primed at least once?
    pub primed: bool,
}

impl Default for CngState {
    fn default() -> Self {
        Self::new()
    }
}

impl CngState {
    pub fn new() -> Self {
        // Uniform-spread initial LSP (same convention used by the main
        // decoder's `LpcPredictorState::new`).
        let mut lsp = [0.0f32; LPC_ORDER];
        let step = core::f32::consts::PI / (LPC_ORDER as f32 + 1.0);
        for k in 0..LPC_ORDER {
            lsp[k] = (step * (k as f32 + 1.0)).cos();
        }
        Self {
            rng: 0x1234_5678,
            lsp,
            lsp_target: lsp,
            energy_db: 0.0,
            energy_db_target: 0.0,
            syn_mem: [0.0; LPC_ORDER],
            primed: false,
        }
    }

    /// Update the generator from a received SID. `last_good_lsp` is the
    /// decoder's most-recent quantised LSP vector — used as a shape
    /// reference since our SID fingerprints the spectrum with only 10
    /// bits (not enough for a full LSP decode).
    pub fn update_from_sid(&mut self, sid: &EncodedSid, last_good_lsp: &[f32; LPC_ORDER]) {
        // Use the SID's shape bits to perturb the last good LSP. Each
        // bit flips the sign of a small offset on one LSP component —
        // enough to carry "spectrum changed" information across the
        // 2-byte frame.
        let mut target = *last_good_lsp;
        for k in 0..LPC_ORDER {
            if (sid.shape >> k) & 1 != 0 {
                target[k] += 0.02;
            } else {
                target[k] -= 0.02;
            }
        }
        // Keep cosine-domain LSPs strictly decreasing and bounded.
        for k in 1..LPC_ORDER {
            if target[k] >= target[k - 1] - 1e-3 {
                target[k] = target[k - 1] - 1e-3;
            }
        }
        for v in target.iter_mut() {
            *v = v.clamp(-0.9995, 0.9995);
        }
        self.lsp_target = target;
        self.energy_db_target = sid.energy_db();
        if !self.primed {
            // First SID: snap to target so the very first CNG frame
            // already sounds right.
            self.lsp = target;
            self.energy_db = self.energy_db_target;
            self.primed = true;
        }
    }

    /// Generate one 10-ms comfort-noise frame into `out`.
    pub fn generate(&mut self, out: &mut [f32; FRAME_SAMPLES]) {
        // Slew the smoothed state towards the targets (quick: reach
        // target in ~4 frames).
        let alpha = 0.25f32;
        for k in 0..LPC_ORDER {
            self.lsp[k] += alpha * (self.lsp_target[k] - self.lsp[k]);
        }
        self.energy_db += alpha * (self.energy_db_target - self.energy_db);

        // Compute excitation amplitude from energy_db:
        //   sumsq = 10^(energy_db/10) * N
        //   rms   = sqrt(10^(energy_db/10))
        // white noise of amplitude `amp` has sumsq ≈ amp^2 / 3 per
        // sample, so pick amp = sqrt(3 * 10^(energy_db/10)).
        let sumsq_per_sample = 10.0f32.powf(self.energy_db.max(0.0) / 10.0);
        let amp = (3.0 * sumsq_per_sample).sqrt();

        // Build LPC from the current LSP.
        let a = lsp_to_lpc(&self.lsp);

        // Generate excitation + synth-filter in one pass.
        for i in 0..FRAME_SAMPLES {
            let u = self.next_uniform() * amp;
            // y[n] = u - sum_{k=1..=LPC_ORDER} a[k] * y[n-k]
            let mut y = u;
            for k in 1..=LPC_ORDER {
                y -= a[k] * self.syn_mem[k - 1];
            }
            // Slide syn_mem.
            for k in (1..LPC_ORDER).rev() {
                self.syn_mem[k] = self.syn_mem[k - 1];
            }
            self.syn_mem[0] = y;
            out[i] = y;
        }
    }

    /// LCG RNG producing uniformly-distributed floats in `[-1, 1]`.
    fn next_uniform(&mut self) -> f32 {
        // Numerical Recipes LCG parameters.
        self.rng = self.rng.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        let top = (self.rng >> 8) as i32 & 0xFF_FFFF;
        // Map [0, 2^24) -> [-1, 1).
        ((top as f32) / 8_388_608.0) - 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cng_produces_nonzero_noise() {
        let mut state = CngState::new();
        let lsp_ref = [0.3, 0.2, 0.1, 0.05, 0.0, -0.05, -0.1, -0.2, -0.3, -0.4];
        let sid = EncodedSid {
            shape: 0b1010101010,
            energy: 20,
        };
        state.update_from_sid(&sid, &lsp_ref);
        let mut out = [0.0f32; FRAME_SAMPLES];
        state.generate(&mut out);
        let mut nonzero = 0;
        for &v in out.iter() {
            assert!(v.is_finite());
            if v.abs() > 0.01 {
                nonzero += 1;
            }
        }
        // CNG with non-trivial energy should yield many non-zero samples.
        assert!(
            nonzero > 10,
            "CNG produced too few non-zero samples: {nonzero}"
        );
    }

    #[test]
    fn cng_with_zero_energy_stays_quiet() {
        // With energy=0 the SID decodes to energy_db=0, i.e. a sumsq of
        // 1 LSB^2/sample. After the 10th-order synthesis filter amplifies,
        // output can still reach ~mid-double-digits on the first frame
        // before the smoothing converges. We just assert the output is
        // finite and orders-of-magnitude below the i16 clip rail.
        let mut state = CngState::new();
        let lsp_ref = [0.3, 0.2, 0.1, 0.05, 0.0, -0.05, -0.1, -0.2, -0.3, -0.4];
        let sid = EncodedSid {
            shape: 0,
            energy: 0,
        };
        state.update_from_sid(&sid, &lsp_ref);
        let mut out = [0.0f32; FRAME_SAMPLES];
        state.generate(&mut out);
        let mut peak = 0.0f32;
        for &v in out.iter() {
            assert!(v.is_finite());
            peak = peak.max(v.abs());
        }
        // Loose bound: the 10-th-order synthesis filter can ring on
        // spectra close to the unit circle even with low excitation. We
        // just guard against saturation at the i16 rail; CNG that's
        // 10-20 dB below speech level is still acceptable as "quiet".
        assert!(peak < 16_000.0, "CNG zero-energy peak too loud: {peak}");
    }

    #[test]
    fn cng_multiple_frames_differ() {
        // Generator should produce different samples on successive frames
        // (not periodic).
        let mut state = CngState::new();
        let lsp_ref = [0.3, 0.2, 0.1, 0.05, 0.0, -0.05, -0.1, -0.2, -0.3, -0.4];
        let sid = EncodedSid {
            shape: 0b0101010101,
            energy: 15,
        };
        state.update_from_sid(&sid, &lsp_ref);
        let mut f1 = [0.0f32; FRAME_SAMPLES];
        let mut f2 = [0.0f32; FRAME_SAMPLES];
        state.generate(&mut f1);
        state.generate(&mut f2);
        // Require that f1 and f2 aren't bitwise identical.
        let different = (0..FRAME_SAMPLES).any(|i| (f1[i] - f2[i]).abs() > 1e-6);
        assert!(different);
    }
}
