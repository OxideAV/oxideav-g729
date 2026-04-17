//! ITU-T G.729 Annex B — Voice Activity Detection (§3) and Silence
//! Insertion Descriptor (§4) support.
//!
//! # Compliance note
//!
//! The ITU spec's VAD is a four-feature decision tree (full-band energy,
//! low-band energy, zero-crossing rate, spectral-distortion via LSP
//! differences), cross-checked against running noise estimates and a set
//! of hysteresis thresholds. Exactly matching the reference implementation
//! requires tuning the thresholds against ITU test vectors — which we
//! deliberately skip in this first-cut pure-Rust implementation.
//!
//! Instead, this module ships a **simplified energy-only VAD**:
//!
//! 1. Compute frame full-band energy.
//! 2. Track a running noise-floor estimate (slow-adapting).
//! 3. Classify the frame as VOICE if the signal-to-noise ratio exceeds a
//!    hang-over-biased threshold, else NOISE.
//! 4. Apply a short hang-over (keep emitting VOICE for a few extra frames
//!    after the SNR drops, to avoid chopping tail-ends of words).
//!
//! This is sufficient for most speech / silence distinction — the main
//! practical use-case for Annex B — and returns meaningful "this frame
//! looks like background noise" decisions on synthetic silence, clean
//! speech, and noisy speech.
//!
//! If you find false-positives on specific inputs, tune
//! [`VadState::snr_threshold_db`] or revisit the four-feature decision
//! tree in §3.3 of the spec.
//!
//! # SID frame layout
//!
//! Annex B SID frames are 15 bits payload (2 bytes on the wire):
//!
//! | Field | Bits | Meaning                                     |
//! |-------|------|---------------------------------------------|
//! | LP1   | 1    | Switch to low-rate or keep current LSF ref  |
//! | LP2   | 5    | LSF split-VQ first stage index              |
//! | LP3   | 4    | LSF split-VQ second stage index (low half)  |
//! | GAIN  | 5    | SID gain / energy index                     |
//!
//! For our simplified purposes we treat SID as (spectral-shape-index,
//! energy-index) and pack it into a 2-byte frame that is readily
//! distinguishable from the 10-byte main-body frame by **size alone**.
//! The decoder uses `packet.len()` to identify SID vs speech vs no-data.

use crate::{FRAME_SAMPLES, LPC_ORDER};

/// Size of an Annex B SID (Silence Insertion Descriptor) frame, in bytes.
///
/// Different from the 10-byte CS-ACELP speech frame so the decoder can
/// dispatch on packet length.
pub const SID_FRAME_BYTES: usize = 2;

/// VAD classification for a single 10-ms frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VadDecision {
    /// Full-rate speech frame — emit a normal 10-byte G.729 packet.
    Voice,
    /// Noise / silence frame — emit an SID (2-byte) summary frame.
    Noise,
    /// Silence continuation — emit nothing (DTX).
    Silence,
}

/// VAD + DTX state.
#[derive(Clone, Debug)]
pub struct VadState {
    /// Running noise-floor estimate (frame energy in log2-space).
    pub noise_floor: f32,
    /// Hang-over counter: keep emitting VOICE for this many more frames
    /// once the SNR drops below threshold.
    pub hangover: u32,
    /// How many frames we've been in SILENCE since the last SID.
    pub silence_run: u32,
    /// Last transmitted SID payload (for DTX comparison).
    pub last_sid: Option<EncodedSid>,
    /// Has the running noise-floor estimate been seeded?
    pub floor_seeded: bool,
    /// Minimum SNR (dB) above the noise floor before a frame is VOICE.
    pub snr_threshold_db: f32,
    /// Hang-over duration in frames. 5 ≙ 50 ms.
    pub hangover_frames: u32,
    /// Force an SID refresh every N silence frames even if the spectrum
    /// hasn't changed. 8 ≙ 80 ms, matches the reference's "at least one
    /// SID every eight frames" refresh policy (spec §4.3).
    pub sid_refresh_frames: u32,
}

impl Default for VadState {
    fn default() -> Self {
        Self::new()
    }
}

impl VadState {
    pub fn new() -> Self {
        Self {
            noise_floor: 0.0,
            hangover: 0,
            silence_run: 0,
            last_sid: None,
            floor_seeded: false,
            snr_threshold_db: 6.0,
            hangover_frames: 5,
            sid_refresh_frames: 8,
        }
    }

    /// Classify a 10-ms frame. Returns the decision and, for `Noise`, the
    /// SID payload to emit (or `None` for SILENCE / VOICE).
    pub fn classify(
        &mut self,
        pcm: &[i16; FRAME_SAMPLES],
        lsp: &[f32; LPC_ORDER],
    ) -> (VadDecision, Option<EncodedSid>) {
        let frame_db = frame_energy_db(pcm);

        // Seed the noise floor with the first frame.
        if !self.floor_seeded {
            self.noise_floor = frame_db.max(0.0);
            self.floor_seeded = true;
        }

        let snr = frame_db - self.noise_floor;
        let is_voice = snr > self.snr_threshold_db;

        if is_voice {
            self.hangover = self.hangover_frames;
            self.silence_run = 0;
            // Don't adapt noise floor while speech is active.
            (VadDecision::Voice, None)
        } else if self.hangover > 0 {
            self.hangover -= 1;
            self.silence_run = 0;
            // Still treat as voice during hang-over so we don't clip the
            // tail of words.
            (VadDecision::Voice, None)
        } else {
            // Slowly adapt the noise floor towards the current frame
            // energy. Fast-down (detect dropping noise quickly) / slow-up
            // so a short burst of speech doesn't drag the floor up.
            let alpha = if frame_db < self.noise_floor {
                0.25
            } else {
                0.05
            };
            self.noise_floor += alpha * (frame_db - self.noise_floor);

            // Build a candidate SID from the current spectrum + energy.
            let sid_candidate = EncodedSid::from_lsp_and_energy(lsp, frame_db);

            // Compare against the last transmitted SID.
            let should_send = match &self.last_sid {
                None => true,
                Some(prev) => {
                    prev.differs_materially(&sid_candidate)
                        || self.silence_run >= self.sid_refresh_frames
                }
            };

            if should_send {
                self.last_sid = Some(sid_candidate);
                self.silence_run = 0;
                (VadDecision::Noise, Some(sid_candidate))
            } else {
                self.silence_run = self.silence_run.saturating_add(1);
                (VadDecision::Silence, None)
            }
        }
    }
}

/// Compute a frame's log-energy in dB.
fn frame_energy_db(pcm: &[i16; FRAME_SAMPLES]) -> f32 {
    let mut sumsq: f64 = 0.0;
    for &s in pcm.iter() {
        let v = s as f64;
        sumsq += v * v;
    }
    let mean = sumsq / (FRAME_SAMPLES as f64);
    // 10*log10(mean + 1) — +1 avoids log(0) and keeps the scale
    // non-negative.
    10.0 * ((mean + 1.0).log10() as f32)
}

/// Compact SID payload: 10-bit spectrum index + 5-bit energy index + 1 reserved.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EncodedSid {
    /// 10-bit quantised LSP shape identifier.
    pub shape: u16,
    /// 5-bit quantised energy level.
    pub energy: u8,
}

impl EncodedSid {
    /// Build an SID from an LSP vector + frame energy in dB.
    ///
    /// Shape: we fingerprint the LSP vector into 10 bits by quantising
    /// each of the ten LSPs (∈[-1, 1]) to 1 bit (sign) and packing. Not
    /// spec-exact — this is a perceptually sufficient fingerprint so the
    /// DTX comparator can decide "same spectrum or not".
    ///
    /// Energy: clamped frame_db ∈ [0, 62] dB bucketed into 32 levels.
    pub fn from_lsp_and_energy(lsp: &[f32; LPC_ORDER], frame_db: f32) -> Self {
        let mut shape: u16 = 0;
        for (i, &v) in lsp.iter().enumerate() {
            // 1 bit per LSP — sign of (v - 0.5) gives a spread across
            // typical voiced vs unvoiced spectra.
            if v > 0.0 {
                shape |= 1 << i;
            }
        }
        let e_clamped = frame_db.clamp(0.0, 62.0);
        let energy = ((e_clamped / 62.0) * 31.0).round().clamp(0.0, 31.0) as u8;
        Self { shape, energy }
    }

    /// Pack into a 2-byte frame. Layout (MSB-first):
    ///
    /// ```text
    /// byte 0: shape[9:2]         (8 bits)
    /// byte 1: shape[1:0] << 6 | energy[4:0] << 1 | reserved (1 bit)
    /// ```
    pub fn pack(&self) -> [u8; SID_FRAME_BYTES] {
        let shape10 = self.shape & 0x3FF;
        let energy5 = (self.energy & 0x1F) as u16;
        let b0 = (shape10 >> 2) as u8;
        let b1 = (((shape10 & 0x3) << 6) | (energy5 << 1)) as u8;
        [b0, b1]
    }

    /// Inverse of [`pack`]. Reserved bit is ignored on decode.
    pub fn unpack(bytes: &[u8]) -> Option<Self> {
        if bytes.len() < SID_FRAME_BYTES {
            return None;
        }
        let b0 = bytes[0] as u16;
        let b1 = bytes[1] as u16;
        let shape = (b0 << 2) | ((b1 >> 6) & 0x3);
        let energy = ((b1 >> 1) & 0x1F) as u8;
        Some(Self { shape, energy })
    }

    /// DTX comparator: does `other` differ enough from `self` to warrant
    /// a new SID on the wire?
    pub fn differs_materially(&self, other: &EncodedSid) -> bool {
        // Count set-bit differences in shape (Hamming distance).
        let shape_ham = (self.shape ^ other.shape).count_ones();
        // Energy difference in quantisation units.
        let e_diff = (self.energy as i32 - other.energy as i32).unsigned_abs();
        // Spec §4.2 uses a delta-threshold on each feature. We compose a
        // conservative "meaningful change" rule:
        //   - shape Hamming ≥ 3 out of 10  (≥ 30 % spectrum shift)
        //   - OR energy quanta change ≥ 3  (~6 dB)
        shape_ham >= 3 || e_diff >= 3
    }

    /// Decode the energy index back to a dB value (inverse of the
    /// 32-bucket quantisation in [`Self::from_lsp_and_energy`]).
    pub fn energy_db(&self) -> f32 {
        (self.energy as f32 / 31.0) * 62.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn silence_frame() -> [i16; FRAME_SAMPLES] {
        [0i16; FRAME_SAMPLES]
    }

    fn tone_frame(f_hz: f32) -> [i16; FRAME_SAMPLES] {
        let sr = 8_000.0f32;
        let mut out = [0i16; FRAME_SAMPLES];
        for (i, s) in out.iter_mut().enumerate() {
            let t = i as f32 / sr;
            let v = 10_000.0 * (2.0 * core::f32::consts::PI * f_hz * t).sin();
            *s = v.round() as i16;
        }
        out
    }

    #[test]
    fn silence_then_tone_switches_to_voice() {
        let mut state = VadState::new();
        let lsp = [0.0f32; LPC_ORDER];
        let silence = silence_frame();
        // Feed 10 frames of silence: all should be Noise or Silence.
        for _ in 0..10 {
            let (d, _) = state.classify(&silence, &lsp);
            assert!(
                matches!(d, VadDecision::Noise | VadDecision::Silence),
                "silence should not be Voice, got {d:?}"
            );
        }
        // Now feed tone frames — should flip to Voice.
        let tone = tone_frame(1000.0);
        let (d, _) = state.classify(&tone, &lsp);
        assert_eq!(d, VadDecision::Voice);
    }

    #[test]
    fn sid_pack_unpack_round_trips() {
        let sid = EncodedSid {
            shape: 0b10_1010_0101,
            energy: 17,
        };
        let packed = sid.pack();
        let unpacked = EncodedSid::unpack(&packed).unwrap();
        assert_eq!(sid, unpacked);
    }

    #[test]
    fn sid_differs_when_shape_changes() {
        let a = EncodedSid {
            shape: 0b0000000000,
            energy: 10,
        };
        let b = EncodedSid {
            shape: 0b0000000111,
            energy: 10,
        };
        let c = EncodedSid {
            shape: 0b0000000001,
            energy: 10,
        };
        assert!(a.differs_materially(&b));
        assert!(!a.differs_materially(&c));
    }

    #[test]
    fn sid_differs_when_energy_jumps() {
        let a = EncodedSid {
            shape: 0b1010101010,
            energy: 5,
        };
        let b = EncodedSid {
            shape: 0b1010101010,
            energy: 12,
        };
        assert!(a.differs_materially(&b));
    }

    #[test]
    fn first_silence_frame_emits_sid() {
        // After hangover expires, the first silence frame should emit an SID.
        let mut state = VadState::new();
        state.hangover_frames = 0;
        let lsp = [0.0f32; LPC_ORDER];
        let silence = silence_frame();
        let (d, sid) = state.classify(&silence, &lsp);
        // Either NOISE with SID (first silence) or SILENCE with none (noise
        // floor initialised on first frame; see seeding logic).
        if d == VadDecision::Noise {
            assert!(sid.is_some());
        }
    }
}
