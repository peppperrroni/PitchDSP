// PitchDSP.h — Streaming pitch detector (YIN algorithm)
//
// Usage:
//   PitchDetectorConfig cfg = pitchDetectorDefaultConfig();
//   PitchDetector* d = pitchDetectorCreate(8192, 48000.0f, cfg);
//
//   // In your audio callback (any chunk size):
//   pitchDetectorProcess(d, micSamples, sampleCount);
//
//   // In your UI update loop:
//   PitchResult r = pitchDetectorGetResult(d);
//   if (r.sequence != lastSequence) {
//       lastSequence = r.sequence;
//       if (r.hz > 0) { /* r.hz is the fundamental in Hz */ }
//   }
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
                          // > 0.72 clean tone, 0.50–0.72 usable, < 0.50 marginal/noise.
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

    /// Primary CMNDF threshold. Default: 0.15 (textbook YIN starting point).
    /// Lower = stricter. Range: 0.05 (very strict) to 0.25 (permissive).
    float yinThreshold;

    /// Fallback threshold — if no minimum is found below yinThreshold, accept the
    /// global CMNDF minimum if it is below fallbackThreshold. Default: 0.25.
    /// Helps strings whose CMNDF minimum sits between yinThreshold and 0.25
    /// (e.g. wound G string under a microphone).
    float fallbackThreshold;

    /// Harmonic correction threshold. Default: 0.45.
    /// After finding a period T, the detector checks if T/N (N = 2..5) has a
    /// CMNDF value below this threshold. If so, it prefers the shorter period
    /// (higher frequency = the true fundamental). Corrects the octave-down
    /// error common on wound strings (e.g. G3 detected as G1).
    /// Set to 0.0 to disable harmonic correction.
    float harmonicCorrThreshold;

    /// Minimum detectable frequency in Hz. Default: 25.0 Hz.
    /// Limits the maximum lag searched: maxTau = sampleRate / minHz.
    /// This excludes the near-halfWindow region (tau ≈ windowSize/2) where
    /// circular autocorrelation produces anomalous low-CMNDF artifacts.
    /// 25 Hz covers bass B0 (30.87 Hz) with headroom; increase for guitar-only
    /// use (e.g. 65 Hz = low E2). Set to 0 to use the full halfWindow.
    float minHz;

    /// Analysis rate = sampleRate / (windowSize / hopDivisor).
    /// Default: 8 → ~47fps at 48kHz with windowSize=8192.
    int   hopDivisor;
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

/// Retrieve the latest pitch detection result (struct copy, safe from any thread).
PitchResult pitchDetectorGetResult(PitchDetector* detector);

#ifdef __cplusplus
}
#endif

#endif /* PITCH_DSP_H */
