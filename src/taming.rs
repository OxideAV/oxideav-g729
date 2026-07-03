//! Encoder **taming procedure** — the adaptive-codebook gain bound
//! that keeps the long-term-prediction loop contractive on stationary
//! periodic signals, so a channel error in the decoder's excitation
//! cannot be amplified indefinitely.
//!
//! ## Doc source — `docs/audio/g729/taming-procedure.md`
//!
//! The Recommendation prose never describes taming (the string appears
//! only in annex file-inventory lists), so this module implements the
//! clean-room algorithmic description staged in
//! `docs/audio/g729/taming-procedure.md`:
//!
//! * **State** — four accumulators, one per zone of the pitch-delay
//!   axis (the staged `tab_zone` partition,
//!   [`crate::tables::TAMING_ZONE_TABLE`]), each tracking the
//!   worst-case excitation-error energy a decoder could have
//!   accumulated for delays in that zone. Initialised to the neutral
//!   `1.0` (the additive floor of the recursion).
//! * **Test** (per subframe, before gain quantisation) — `tameflag`
//!   is raised when the maximum accumulator over the zones spanned by
//!   the subframe's pitch delay exceeds a threshold.
//! * **Update** (per subframe, after the gains are fixed) — the
//!   spanned zones' worst-case energy is propagated one period through
//!   the long-term loop: `E_new = 1.0 + ĝ_p²·E_prev`, written back to
//!   every spanned zone. With `ĝ_p < 1` the recursion contracts toward
//!   its small steady state (self-releasing); with `ĝ_p ≳ 1` it grows
//!   across successive subframes — the condition the test detects.
//! * **Bound** (consumed by [`crate::gain_quantize`]) — when
//!   `tameflag` is set, the §3.9.2 conjugate-structure search excludes
//!   every codebook pair whose reconstructed `ĝ_p` exceeds
//!   [`GPCLIP`] = 0.95 (Q14 `15564`, the doc's pinned ceiling), forcing
//!   the loop back into the contractive régime.
//!
//! ## Zone spanning
//!
//! A delay's adaptive-codebook vector is interpolated from the past
//! excitation through the 31-tap eq (40) filter, so the memory it
//! reuses spans delays `k … k + L_INTERPOL − 1` where `k` is the
//! integer anchor of the fractional delay (the same `(int_t, frac) →
//! k` normalisation as the eq (40) interpolator itself: `frac = −1`
//! anchors at `int_t − 1`). The doc pins the *existence* of the
//! two-zone span but flags the exact index arithmetic as an open gap;
//! the `k … k + L_INTERPOL − 1` window used here is the crate's own
//! eq (40) read window, whose maximum index `143 + 9 = 152` lands
//! exactly on the last entry of the staged 153-entry `tab_zone` —
//! internal evidence the staged table was dimensioned for precisely
//! this reach.
//!
//! ## Open constants (doc "honest gaps"), fixed black-box
//!
//! Two quantities the doc deliberately leaves unpinned are fixed here
//! by black-box behavioural validation against the staged conformance
//! corpus (the reference's own `.BIT` streams):
//!
//! * the accumulator recursion's loop-gain factor — implemented as
//!   `ĝ_p²` itself (the per-period energy gain of the recursive loop,
//!   the doc's "monotonically increasing function of ĝ_p²" at its
//!   simplest form);
//! * the trigger threshold [`THRESH_ERR`] = 60000 (the doc's
//!   commonly-attributed-but-unconfirmed figure). The corpus pins a
//!   consistent *lower bound* for it: replaying this recursion over
//!   the reference's own decoded `(delay, ĝ_p)` sequences, the
//!   reference keeps choosing `ĝ_p > GPCLIP` while the simulated
//!   accumulator climbs to ≈ 18 200 (TAME vector) — so the reference's
//!   threshold in these units exceeds 18 200, and the reference
//!   never actually tames anywhere in the staged corpus (its peak
//!   simulated accumulator stays < 60000 everywhere). Our own
//!   encode path peaks at ≈ 15 800 on TAME, so with 60000 this
//!   encoder also never tames on the corpus — behaviourally identical
//!   to the reference. A forced-lower-threshold sweep (100 … 15000)
//!   measurably *degrades* GB agreement on TAME, corroborating the
//!   never-tame reading. See
//!   `tests/encoder_conformance.rs::reference_taming_fingerprint`.

use crate::tables::{taming_zone, TAMING_ZONES};

/// Interpolation-window reach of the eq (40) adaptive-codebook filter:
/// the number of consecutive integer delays whose excitation memory a
/// single fractional delay reuses (§3.7, `L_INTERPOL = 10`).
pub const L_INTERPOL: usize = 10;

/// Taming ceiling on the quantised adaptive-codebook gain — the
/// doc-pinned `GPCLIP` (Q14 `15564`, ≈ 0.94995): the maximum `ĝ_p` the
/// §3.9.2 search may reconstruct while `tameflag` is raised. Just
/// below the `ĝ_p = 1` stability boundary of the long-term loop.
pub const GPCLIP: f32 = 15564.0 / 16384.0;

/// Taming trigger threshold on the per-zone worst-case error energy
/// (natural units of the `E = 1 + ĝ_p²·E` recursion). The doc leaves
/// the exact value unpinned; 60000 is the commonly-attributed figure,
/// and the staged corpus black-box-pins a consistent lower bound of
/// ≈ 18 200 for it (the reference's own bitstreams keep choosing
/// `ĝ_p > GPCLIP` up to that simulated accumulator level, so its
/// threshold must lie above — see the module docs and
/// `tests/encoder_conformance.rs::reference_taming_fingerprint`).
pub const THRESH_ERR: f64 = 60_000.0;

/// Saturation ceiling on the accumulators — the float stand-in for the
/// fixed-point reference's 32-bit saturation, keeping the state finite
/// (and therefore self-releasing) under arbitrarily long high-gain
/// runs.
const ERR_CEILING: f64 = 1.0e12;

/// The per-zone worst-case excitation-error state of the taming
/// procedure (doc §1).
#[derive(Debug, Clone)]
pub struct Taming {
    /// Worst-case accumulated error energy per pitch-delay zone.
    exc_err: [f64; TAMING_ZONES],
}

impl Default for Taming {
    fn default() -> Self {
        Self::new()
    }
}

impl Taming {
    /// Fresh start-up state: every accumulator at the neutral `1.0`
    /// (no accumulated error yet — the additive floor of the
    /// recursion).
    #[must_use]
    pub fn new() -> Self {
        Self {
            exc_err: [1.0; TAMING_ZONES],
        }
    }

    /// The `(zone_lo, zone_hi)` pair spanned by a fractional pitch
    /// delay `(int_t, frac)` — the zones of the integer anchor `k` and
    /// of the far end of its `L_INTERPOL`-sample interpolation reach.
    fn spanned_zones(int_t: i32, frac: i32) -> (usize, usize) {
        // Same anchor normalisation as the eq (40) interpolator:
        // frac = −1 re-anchors one sample earlier.
        let k = if frac < 0 { int_t - 1 } else { int_t };
        let k = k.max(0) as usize;
        let lo = taming_zone(k);
        let hi = taming_zone(k + L_INTERPOL - 1);
        (lo, hi)
    }

    /// Doc §3 **test**: returns the `tameflag` for a subframe about to
    /// quantise its gains with pitch delay `(int_t, frac ∈ {−1,0,1})` —
    /// `true` when the worst-case accumulated error over the spanned
    /// zones exceeds [`THRESH_ERR`] (i.e. the long-term loop for this
    /// delay has been running hot enough that the gain must be tamed).
    #[must_use]
    pub fn test(&self, int_t: i32, frac: i32) -> bool {
        let (lo, hi) = Self::spanned_zones(int_t, frac);
        self.exc_err[lo].max(self.exc_err[hi]) > THRESH_ERR
    }

    /// Doc §2 **update**: after the subframe's quantised adaptive gain
    /// `ĝ_p` is fixed, propagates the spanned zones' worst-case error
    /// one pitch period through the loop —
    /// `E_new = 1.0 + ĝ_p²·max(E_spanned)` — and stores it back into
    /// every spanned zone.
    pub fn update(&mut self, int_t: i32, frac: i32, gp_hat: f32) {
        let (lo, hi) = Self::spanned_zones(int_t, frac);
        let e_prev = self.exc_err[lo].max(self.exc_err[hi]);
        let loop_gain = f64::from(gp_hat) * f64::from(gp_hat);
        let e_new = (1.0 + loop_gain * e_prev).min(ERR_CEILING);
        self.exc_err[lo] = e_new;
        self.exc_err[hi] = e_new;
    }

    /// Read-only view of the per-zone accumulators (diagnostics /
    /// tests).
    #[must_use]
    pub fn accumulators(&self) -> &[f64; TAMING_ZONES] {
        &self.exc_err
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fresh state never tames: every accumulator is at the neutral
    /// floor, far below the threshold, for every legal delay.
    #[test]
    fn fresh_state_is_untamed() {
        let t = Taming::new();
        for int_t in 20..=143 {
            for frac in [-1, 0, 1] {
                assert!(!t.test(int_t, frac), "delay {int_t} frac {frac}");
            }
        }
    }

    /// A sustained high pitch gain (`ĝ_p = 1.2`, the eq (43) analysis
    /// ceiling) at a fixed delay grows the accumulator geometrically
    /// and trips the tameflag; the doc's `ĝ_p ≳ 1` divergence régime.
    #[test]
    fn sustained_high_gain_trips_tameflag() {
        let mut t = Taming::new();
        let mut tripped = None;
        for sub in 0..200 {
            if t.test(40, 0) {
                tripped = Some(sub);
                break;
            }
            t.update(40, 0, 1.2);
        }
        let sub = tripped.expect("gp=1.2 must eventually trip the taming test");
        // 1 + 1.44·E grows past 60000 in ~log(60000)/log(1.44) ≈ 31
        // periods.
        assert!((20..60).contains(&sub), "tripped after {sub} subframes");
    }

    /// A contractive gain (`ĝ_p < 1`) converges to the small steady
    /// state `1/(1 − ĝ_p²)` and never tames — and releases a
    /// previously-hot accumulator (the doc's self-releasing property).
    #[test]
    fn low_gain_converges_and_releases() {
        let mut t = Taming::new();
        // Run hot until tamed.
        for _ in 0..100 {
            t.update(40, 0, 1.2);
        }
        assert!(t.test(40, 0), "should be hot after the high-gain run");
        // Tamed subframes use gp ≤ GPCLIP < 1: the loop contracts.
        for _ in 0..300 {
            t.update(40, 0, 0.5);
        }
        assert!(!t.test(40, 0), "low-gain run must release the taming");
        // Steady state of E = 1 + 0.25·E is 4/3.
        let e = t.accumulators()[Taming::spanned_zones(40, 0).0];
        assert!((e - 4.0 / 3.0).abs() < 1e-6, "steady state {e}");
    }

    /// A delay whose `L_INTERPOL`-sample reach crosses a zone boundary
    /// touches both zones: updating at the boundary heats both, and
    /// the test sees the heat from either side.
    #[test]
    fn boundary_delay_spans_two_zones() {
        // k = 35: reach 35 … 44 crosses the 40-boundary (zones 0 + 1).
        let (lo, hi) = Taming::spanned_zones(35, 0);
        assert_eq!((lo, hi), (0, 1));
        // k = 20: reach 20 … 29 stays in zone 0.
        assert_eq!(Taming::spanned_zones(20, 0), (0, 0));
        // frac = −1 re-anchors: int_t = 40, frac = −1 → k = 39,
        // reach 39 … 48 spans zones 0 + 1.
        assert_eq!(Taming::spanned_zones(40, -1), (0, 1));
        // The maximum legal reach lands on the last table entry.
        assert_eq!(Taming::spanned_zones(143, 1), (3, 3));

        let mut t = Taming::new();
        for _ in 0..100 {
            t.update(35, 0, 1.2);
        }
        // Both zone-0-only and zone-1-only delays now see the heat.
        assert!(t.test(20, 0), "zone 0 heated through the boundary span");
        assert!(t.test(60, 0), "zone 1 heated through the boundary span");
        // Zones 2/3 are untouched.
        assert!(!t.test(100, 0));
        assert!(!t.test(143, 0));
    }

    /// The accumulator saturates at the ceiling instead of running to
    /// infinity, so an arbitrarily long high-gain run stays releasable.
    #[test]
    fn accumulator_saturates_finite() {
        let mut t = Taming::new();
        for _ in 0..10_000 {
            t.update(40, 0, 1.2);
        }
        let e = t.accumulators()[0];
        assert!(e.is_finite() && e <= 1.0e12, "accumulator {e}");
    }

    /// The pinned ceiling is the doc's Q14 value: it re-quantises to
    /// exactly 15564 on the Q14 grid and sits below the `ĝ_p = 1`
    /// stability boundary.
    #[test]
    fn gpclip_is_the_pinned_q14_value() {
        let q14 = (f64::from(GPCLIP) * 16384.0).round() as i32;
        assert_eq!(q14, 15564);
        assert!(q14 < 16384, "ceiling must sit below the ĝ_p = 1 boundary");
    }
}
