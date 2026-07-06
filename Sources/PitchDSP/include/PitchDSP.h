// PitchDSP.h — Streaming pitch detector (YIN algorithm, direct method)
//
// Usage:
//   PitchDetectorConfig cfg = pitchDetectorDefaultConfig();
//   PitchDetector* d = pitchDetectorCreate(8192, 48000.0f, cfg);
//
//   // In your audio callback (any chunk size):
//   pitchDetectorProcess(d, micSamples, sampleCount);
//   PitchResult results[16];
//   int n = pitchDetectorDrainResults(d, results, 16);
//   for (int i = 0; i < n; i++) { /* results[i].hz > 0 = pitch */ }
//
//   pitchDetectorDestroy(d);

#ifndef PITCH_DSP_H
#define PITCH_DSP_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

typedef struct {
    float    hz;          // Fundamental frequency Hz; -1 = no pitch detected.
    float    confidence;  // Quality score 0..1. Combines CMNDF depth and local contrast.
                          // > 0.72 clean tone, 0.60–0.72 usable, < 0.60 rejected by the default minConfidence gate.
    uint32_t sequence;    // Increments on every fresh analysis; compare to detect stale reads.
} PitchResult;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

typedef struct {
    /// RMS gate — analysis skipped if window RMS < this. Default: 1e-5 (~-50 dBFS).
    float powerThreshold;

    /// Peak amplitude gate — analysis skipped if peak < this. Default: 0.005.
    /// Catches decaying tails where RMS still passes but the waveform is too quiet
    /// to have a reliable period.
    float peakThreshold;

    /// Primary CMNDF threshold. Default: 0.20 — deliberately above the textbook
    /// 0.15: keeps the true-period pit acceptable when low-frequency room
    /// rumble pollutes it, preventing far-subharmonic locks (see README).
    /// Lower = stricter. Range: 0.05 (very strict) to 0.25 (permissive).
    float yinThreshold;

    /// Fallback threshold — if no minimum is found below yinThreshold, accept the
    /// global CMNDF minimum if it is below fallbackThreshold. Default: 0.30.
    float fallbackThreshold;

    /// Harmonic correction tolerance. Default: 0.03.
    /// After YIN finds period P, checks P/2, P/3, P/4, P/5. If CMNDF[P/N] is
    /// within this tolerance of CMNDF[P], the shorter period (higher frequency =
    /// truer fundamental) is preferred. Corrects sub-harmonic locking where wound
    /// strings cause YIN to settle on a longer period than the true fundamental.
    /// Set to 0.0 to disable. Lower (e.g. 0.02) if spurious octave-up jumps appear.
    float octaveTolerance;

    /// Minimum confidence to report a result. Default: 0.60.
    /// Confidence = 0.75×(1−CMNDF) + 0.25×(localContrast×8), clamped 0..1.
    /// Frames below this threshold are reported as invalid (hz = -1).
    float minConfidence;

    /// Analysis rate = sampleRate / (windowSize / hopDivisor).
    /// Default: 8 → ~47fps at 48kHz with windowSize=8192.
    int   hopDivisor;

    /// Lower frequency bound in Hz. Default: 25.0 (covers 5-string bass B0 = 30.87 Hz).
    /// Sets the tau search ceiling: maxTau = sampleRate / minHz. This prevents spurious
    /// detections in the region near halfWindow (τ > N/4) where the direct-YIN CMNDF
    /// has a systematic downward bias: for random noise, E[CMNDF(τ)] ≈ (N-τ)/(N-τ/2),
    /// which falls to ~0.67 at τ = halfWindow rather than the ideal 1.0. Without this
    /// cap, sub-bass artifacts at τ≈3000-4096 (11-17 Hz) can pass the fallbackThreshold
    /// gate. Raise to 60-70 Hz for guitar-only applications to tighten the search range.
    float minHz;

    /// Upper frequency bound in Hz. Default: 1000.0.
    /// Sets the tau search floor: minTau = sampleRate / maxHz. Prevents noisy
    /// frames from locking onto spurious tiny-tau dips (implausibly high
    /// frequencies). Note: for a periodic signal whose fundamental is above
    /// maxHz, the detector reports the first period multiple inside the range
    /// (a subharmonic) — by design, reported hz is always within [minHz, maxHz].
    float maxHz;
} PitchDetectorConfig;

// ---------------------------------------------------------------------------
// Opaque detector
// ---------------------------------------------------------------------------

typedef struct PitchDetector PitchDetector;

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

/// Returns a default configuration suitable for guitar/bass tuning.
PitchDetectorConfig pitchDetectorDefaultConfig(void);

/// Create a streaming pitch detector.
/// @param windowSize  Analysis window in samples. Recommended: 8192 (covers B0 = 30.87 Hz at 48kHz).
///                    Must be a multiple of 4 (the internal decimation factor); other values return NULL.
/// @param sampleRate  Audio sample rate in Hz (e.g. 44100.0, 48000.0).
/// @param config      Initial configuration.
/// @return Allocated detector, or NULL on failure. Caller must call pitchDetectorDestroy().
PitchDetector* pitchDetectorCreate(int windowSize, float sampleRate, PitchDetectorConfig config);

/// Free all memory associated with a detector.
void pitchDetectorDestroy(PitchDetector* detector);

/// Apply a new configuration. Takes effect on the next analysis cycle.
void pitchDetectorConfigure(PitchDetector* detector, PitchDetectorConfig config);

/// Feed new audio samples into the detector (streaming, any chunk size).
void pitchDetectorProcess(PitchDetector* detector, const float* samples, int count);

/// Retrieve the latest pitch detection result (struct copy). Same single-thread
/// contract as pitchDetectorProcess — see THREADING note on
/// pitchDetectorDrainResults.
PitchResult pitchDetectorGetResult(PitchDetector* detector);

/// Drain queued results, oldest first. Every analysis cycle queues one result,
/// including invalid ones (hz = -1) — consumers use invalid frames for silence
/// debouncing. Internal queue capacity is 16; when full, the oldest result is
/// dropped. Returns the number of results written to `out` (0..maxCount).
///
/// THREADING: call pitchDetectorProcess / pitchDetectorDrainResults /
/// pitchDetectorGetResult from a single thread (typically the audio callback
/// thread). Cross-thread hand-off is the caller's responsibility.
int pitchDetectorDrainResults(PitchDetector* detector, PitchResult* out, int maxCount);

#ifdef __cplusplus
}
#endif

#endif /* PITCH_DSP_H */
