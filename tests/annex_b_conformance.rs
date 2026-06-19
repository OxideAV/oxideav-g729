//! Conformance harness for the ITU-T G.729 **Annex B** (DTX / CNG)
//! decoder-side framing and SID decode (`oxideav_g729::annex_b`) against
//! the staged Annex B test vectors under
//! `docs/audio/g729/conformance/g729b/`.
//!
//! The Annex B comfort-noise *synthesis* path (§B.4.2.2 SID-LSP
//! dequantization + §B.4.4 CNG excitation) needs numeric tables that are
//! not staged (see the crate README docs-gap note), so this harness
//! validates the parts that are fully determined by the corpus + prose:
//!
//! * every `g729b/*.bit` stream walks frame-by-frame with no framing
//!   error and is consumed to the byte;
//! * the decoded frame count equals the `.out` PCM frame count
//!   (160 bytes / 80-sample frame), i.e. the framing partitions the
//!   stream into exactly one frame per output PCM frame;
//! * the frame-type histogram matches the independently-counted header
//!   shapes ({0, 16, 80} bit counts, normal vs erased sync);
//! * every **active** frame still decodes through the §4.1 parameter
//!   chain without error (Annex B reuses the base 80-bit payload);
//! * every **SID** frame unpacks into in-range Table-B.2 indices and a
//!   §B.4.2.1 energy in the documented [−12, 66] dB range.
//!
//! When the corpus directory is absent (published-crate build) every
//! test logs a skip and exits clean, mirroring the other conformance
//! harnesses.

use std::path::{Path, PathBuf};

use oxideav_g729::annex_b::{
    self, AnnexBDecoder, AnnexBFrame, AnnexBOutput, AnnexBStreamDecoder, ResolvedFrame,
    ACTIVE_BITS, ERASED_SYNC_WORD, FRAME_SAMPLES, SID_GAIN_BITS, SID_L1_BITS, SID_L2_BITS,
    SID_LP0_BITS, SID_WORDS,
};
use oxideav_g729::decode_chain::FrameDecoder;
use oxideav_g729::serial::{FrameKind, SYNC_WORD};

/// PCM bytes per 80-sample / 10 ms frame (16-bit LE).
const PCM_BYTES_PER_FRAME: usize = 160;

/// Walks parent directories from `CARGO_MANIFEST_DIR` looking for
/// `docs/audio/g729/conformance/`. Returns `None` if not found
/// (published-crate mode).
fn conformance_root() -> Option<PathBuf> {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut cursor: &Path = &manifest;
    loop {
        let cand = cursor.join("docs/audio/g729/conformance");
        if cand.join("README.md").is_file() {
            return Some(cand);
        }
        cursor = cursor.parent()?;
    }
}

/// All six Annex B `.bit`/`.out` sequence pairs (the non-`a` G.729-B set
/// and the `…a` Annex-A+B set). `tstseq5`/`tstseq6` are decoder-only.
const SEQUENCES: &[&str] = &[
    "tstseq1", "tstseq1a", "tstseq2", "tstseq2a", "tstseq3", "tstseq3a", "tstseq4", "tstseq4a",
    "tstseq5", "tstseq6",
];

/// Independent reference frame-walk: byte-level header histogram, used to
/// cross-check the module's classification. Returns
/// `(frame_count, active, sid, untransmitted, erased)`.
fn reference_walk(bytes: &[u8]) -> (usize, usize, usize, usize, usize) {
    let mut off = 0usize;
    let (mut frames, mut active, mut sid, mut untx, mut erased) = (0, 0, 0, 0, 0);
    while off + 4 <= bytes.len() {
        let sync = u16::from_le_bytes([bytes[off], bytes[off + 1]]);
        let n_bits = u16::from_le_bytes([bytes[off + 2], bytes[off + 3]]) as usize;
        // All-0x0000 payload on a normal-sync frame is the second §B.4.5
        // erasure shape.
        let all_erased = n_bits > 0
            && (0..n_bits).all(|i| {
                let w = u16::from_le_bytes([bytes[off + 4 + i * 2], bytes[off + 5 + i * 2]]);
                w == 0x0000
            });
        if sync == ERASED_SYNC_WORD || all_erased {
            erased += 1;
        } else {
            assert_eq!(sync, SYNC_WORD, "reference walk hit unknown sync at {off}");
            match n_bits {
                0 => untx += 1,
                16 => sid += 1,
                80 => active += 1,
                other => panic!("reference walk unknown bit count {other} at {off}"),
            }
        }
        frames += 1;
        off += 4 + n_bits * 2;
    }
    assert_eq!(
        off,
        bytes.len(),
        "reference walk did not consume the stream"
    );
    (frames, active, sid, untx, erased)
}

#[test]
fn annex_b_streams_walk_and_match_output_frame_count() {
    let Some(root) = conformance_root() else {
        eprintln!("annex_b_conformance: corpus absent; skipping");
        return;
    };
    let g729b = root.join("g729b");

    let mut checked = 0usize;
    for &seq in SEQUENCES {
        let bit_path = g729b.join(format!("{seq}.bit"));
        let out_path = g729b.join(format!("{seq}.out"));
        let (Ok(bit), Ok(out)) = (std::fs::read(&bit_path), std::fs::read(&out_path)) else {
            continue;
        };

        let frames = annex_b::parse_annex_b_stream(&bit)
            .unwrap_or_else(|e| panic!("{seq}: framing error {e}"));

        // Frame count must equal the output PCM frame count: the framing
        // partitions the stream into exactly one frame per decoded frame.
        assert_eq!(
            out.len() % PCM_BYTES_PER_FRAME,
            0,
            "{seq}.out length {} not a multiple of {PCM_BYTES_PER_FRAME}",
            out.len(),
        );
        let out_frames = out.len() / PCM_BYTES_PER_FRAME;
        assert_eq!(
            frames.len(),
            out_frames,
            "{seq}: parsed {} frames but {seq}.out has {out_frames} PCM frames",
            frames.len(),
        );

        // Cross-check the module's classification against an independent
        // byte-level walk.
        let (ref_frames, ref_active, ref_sid, ref_untx, ref_erased) = reference_walk(&bit);
        assert_eq!(ref_frames, frames.len(), "{seq}: frame count mismatch");

        let mut active = 0usize;
        let mut sid = 0usize;
        let mut untx = 0usize;
        let mut erased = 0usize;
        for frame in &frames {
            match frame {
                AnnexBFrame::Active(_) => active += 1,
                AnnexBFrame::Sid(_) => sid += 1,
                AnnexBFrame::Untransmitted => untx += 1,
                AnnexBFrame::Erased => erased += 1,
            }
        }
        assert_eq!(
            (active, sid, untx, erased),
            (ref_active, ref_sid, ref_untx, ref_erased),
            "{seq}: type histogram mismatch"
        );

        checked += 1;
    }
    assert!(
        checked > 0,
        "no Annex B sequences found under {}",
        g729b.display()
    );
}

#[test]
fn annex_b_active_frames_decode_through_base_chain() {
    let Some(root) = conformance_root() else {
        eprintln!("annex_b_conformance: corpus absent; skipping");
        return;
    };
    let g729b = root.join("g729b");

    let mut checked = 0usize;
    let mut total_active = 0usize;
    for &seq in SEQUENCES {
        let Ok(bit) = std::fs::read(g729b.join(format!("{seq}.bit"))) else {
            continue;
        };
        let frames = annex_b::parse_annex_b_stream(&bit)
            .unwrap_or_else(|e| panic!("{seq}: framing error {e}"));

        // A single decoder threads cross-frame state across the whole
        // stream (active frames advance it; the silence stages would too
        // once CNG synthesis lands).
        let mut dec = FrameDecoder::new();
        for (i, frame) in frames.iter().enumerate() {
            if let AnnexBFrame::Active(bits) = frame {
                // Reuse the base 80-bit decode via the FrameKind path.
                let kind = FrameKind::Active(bits.clone());
                dec.decode_frame_kind(&kind)
                    .unwrap_or_else(|e| panic!("{seq}: active frame #{i} decode error {e}"));
                total_active += 1;
            }
        }
        checked += 1;
    }
    assert!(checked > 0, "no Annex B sequences found");
    assert!(total_active > 0, "no active frames decoded");
}

#[test]
fn annex_b_sid_frames_unpack_in_range() {
    let Some(root) = conformance_root() else {
        eprintln!("annex_b_conformance: corpus absent; skipping");
        return;
    };
    let g729b = root.join("g729b");

    let mut total_sid = 0usize;
    for &seq in SEQUENCES {
        let Ok(bit) = std::fs::read(g729b.join(format!("{seq}.bit"))) else {
            continue;
        };
        let frames = annex_b::parse_annex_b_stream(&bit)
            .unwrap_or_else(|e| panic!("{seq}: framing error {e}"));
        for frame in &frames {
            if let AnnexBFrame::Sid(sid) = frame {
                // Table B.2 field widths bound each index.
                assert!(sid.lp0 < (1 << SID_LP0_BITS), "{seq}: SID lp0 out of range");
                assert!(sid.l1 < (1 << SID_L1_BITS), "{seq}: SID l1 out of range");
                assert!(sid.l2 < (1 << SID_L2_BITS), "{seq}: SID l2 out of range");
                assert!(
                    sid.gain < (1 << SID_GAIN_BITS),
                    "{seq}: SID gain out of range"
                );
                // §B.4.2.1 dequantized energy in the documented window.
                let e = sid.energy_db();
                assert!(
                    (-12.0..=66.0).contains(&e),
                    "{seq}: SID energy {e} dB outside [−12, 66]",
                );
                total_sid += 1;
            }
        }
    }
    assert!(total_sid > 0, "no SID frames found in the Annex B corpus");
}

#[test]
fn annex_b_decoder_resolves_full_corpus() {
    let Some(root) = conformance_root() else {
        eprintln!("annex_b_conformance: corpus absent; skipping");
        return;
    };
    let g729b = root.join("g729b");

    let mut checked = 0usize;
    let mut total_resolved = 0usize;
    let mut total_erased_active = 0usize;
    let mut total_silence = 0usize;
    for &seq in SEQUENCES {
        let Ok(bit) = std::fs::read(g729b.join(format!("{seq}.bit"))) else {
            continue;
        };
        let frames = annex_b::parse_annex_b_stream(&bit)
            .unwrap_or_else(|e| panic!("{seq}: framing error {e}"));

        let mut dec = AnnexBDecoder::new();
        // The first SID seen determines when carried energy becomes
        // available; before that, untransmitted frames carry None.
        let mut seen_sid = false;
        for frame in &frames {
            let resolved = dec.resolve(frame);
            match resolved {
                ResolvedFrame::Active(_) => {}
                ResolvedFrame::Sid { energy_db, .. } => {
                    seen_sid = true;
                    assert!(
                        (-12.0..=66.0).contains(&energy_db),
                        "{seq}: SID energy {energy_db} out of range",
                    );
                }
                ResolvedFrame::Untransmitted { last_energy_db } => {
                    total_silence += 1;
                    // Once a SID has been seen, every untransmitted frame
                    // must carry a concrete energy.
                    if seen_sid {
                        assert!(
                            last_energy_db.is_some(),
                            "{seq}: untransmitted frame after a SID carries no energy",
                        );
                    }
                }
                ResolvedFrame::ErasedActive => total_erased_active += 1,
            }
            total_resolved += 1;
        }
        checked += 1;
    }
    assert!(checked > 0, "no Annex B sequences found");
    assert!(total_resolved > 0);
    // tstseq6 contributes the only erasures in the corpus: one erased
    // SID-sized frame in a silence run (→ untransmitted) and one erased
    // untransmitted frame, so the corpus has no active-erasure but does
    // exercise the silence-inheritance path heavily.
    assert!(total_silence > 0, "corpus had no untransmitted frames");
    // Don't require erased-active frames (the corpus may have none); just
    // assert the counter is well-defined.
    let _ = total_erased_active;
}

#[test]
fn annex_b_stream_decode_full_corpus() {
    let Some(root) = conformance_root() else {
        eprintln!("annex_b_conformance: corpus absent; skipping");
        return;
    };
    let g729b = root.join("g729b");

    let mut checked = 0usize;
    let mut total_speech = 0usize;
    for &seq in SEQUENCES {
        let bit_path = g729b.join(format!("{seq}.bit"));
        let out_path = g729b.join(format!("{seq}.out"));
        let (Ok(bit), Ok(out)) = (std::fs::read(&bit_path), std::fs::read(&out_path)) else {
            continue;
        };

        let mut dec = AnnexBStreamDecoder::new();
        let outputs = dec
            .decode_stream(&bit)
            .unwrap_or_else(|e| panic!("{seq}: stream decode error {e}"));

        // One output block per PCM frame in the reference output.
        assert_eq!(
            outputs.len(),
            out.len() / 160,
            "{seq}: {} decoded blocks but {} PCM frames",
            outputs.len(),
            out.len() / 160,
        );

        for (i, o) in outputs.iter().enumerate() {
            match o {
                AnnexBOutput::Speech(s) => {
                    assert_eq!(s.len(), FRAME_SAMPLES);
                    assert!(
                        s.iter().all(|x| x.is_finite()),
                        "{seq}: speech frame #{i} has non-finite sample",
                    );
                    total_speech += 1;
                }
                AnnexBOutput::ComfortNoisePlaceholder { .. }
                | AnnexBOutput::ErasedActivePlaceholder => {}
            }
        }
        checked += 1;
    }
    assert!(checked > 0, "no Annex B sequences found");
    assert!(total_speech > 0, "no active speech frames decoded");
}

#[test]
fn sid_frame_size_constants() {
    // Guard the Table-B.2 layout constants against drift.
    assert_eq!(SID_WORDS, 16);
    assert_eq!(SID_LP0_BITS + SID_L1_BITS + SID_L2_BITS + SID_GAIN_BITS, 15);
    assert_eq!(ACTIVE_BITS, 80);
}
