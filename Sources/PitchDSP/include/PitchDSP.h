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
    float    confidence;  // 1 − CMNDF_min: 1.0 = clean tone, 0.0 = noise.
    uint32_t sequence;    // Increments on every fresh analysis; compare to detect stale reads.
} PitchResult;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

typedef struct {
    /// RMS gate before YIN. Default: 1e-5 (~-50 dBFS).
    float powerThreshold;

    /// CMNDF minimum threshold. Default: 0.15 (industry standard).
    /// Lower = stricter (fewer false positives, may miss weak signals).
    /// Range: 0.05 (very strict) to 0.20 (permissive).
    float yinThreshold;

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
