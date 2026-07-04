// PitchDSP.h — Streaming pitch detector (Harmonic Product Spectrum)
//
// Based on https://github.com/not-chciken/guitar_tuner
// Reimplemented in C for the Nocturne iOS tuner.
//
// Pipeline:
//   Mic samples -> Ring Buffer -> Power Gate -> Hann Window -> FFT
//   -> Mains Hum Suppression -> Octave-Band Noise Gate
//   -> Interpolate (5x) -> HPS (x1..x5) -> Peak -> Note + Cents
//
// Usage:
//   PitchDetectorConfig cfg = pitchDetectorDefaultConfig();
//   cfg.noteStability = 3;          // optional: tweak before creating
//   PitchDetector* d = pitchDetectorCreate(8192, 48000.0f, cfg);
//
//   // In your audio callback or processing loop:
//   pitchDetectorProcess(d, micSamples, sampleCount);
//
//   // In your UI update (e.g. 30-60 Hz timer):
//   PitchResult r = pitchDetectorGetResult(d);
//   if (r.hz > 0) { /* display note, cents, needle, etc. */ }
//
//   // Adjust parameters at any time (takes effect on next analysis):
//   cfg.powerThreshold = 1e-4f;
//   pitchDetectorConfigure(d, cfg);
//
//   pitchDetectorDestroy(d);
//
// Thread safety: pitchDetectorProcess and pitchDetectorGetResult may be called
// from different threads as long as only ONE thread calls Process at a time.
// pitchDetectorConfigure is safe to call between Process calls.

#ifndef PITCH_DSP_H
#define PITCH_DSP_H

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

typedef struct {
    float hz;           // Detected pitch in Hz. -1 if no pitch detected.
    float cents;        // Deviation from nearest note in cents (-50 .. +50).
    int   midi;         // Nearest MIDI note number (0-127). -1 if no pitch.
    float confidence;   // 0.0 .. 1.0 based on HPS peak prominence.
    float stability;    // 0.0 .. 1.0 note stability score.
} PitchResult;

// ---------------------------------------------------------------------------
// Configuration — all parameters tunable at runtime via pitchDetectorConfigure
// ---------------------------------------------------------------------------

typedef struct {
    /// Minimum mean-squared amplitude of the analysis window to activate.
    /// Below this level the detector returns no-signal.
    /// Typical range: 1e-6 (very sensitive, -60 dBFS) to 1e-3 (-30 dBFS).
    /// Default: 1e-5 (~-50 dBFS — cuts off string-decay body resonances).
    float powerThreshold;

    /// Zero all FFT bins below this frequency (Hz) to suppress mains hum.
    /// Must be below the lowest note you want to detect.
    /// Default: 27.0 Hz (just below A0 = 27.5 Hz, the lowest piano note).
    float mainsHumFreq;

    /// Octave-band noise gate threshold as a fraction of band RMS.
    /// Bins below (threshold × bandRMS) are zeroed inside each octave band.
    /// Typical range: 0.1 (gentle) to 0.4 (aggressive).
    /// Default: 0.2.
    float noiseGateThreshold;

    /// Number of consecutive DSP analyses that must agree on the same MIDI
    /// note before the result is updated. Higher = more stable but slower.
    /// Range: 1 (off) to 8. Default: 2.
    int   noteStability;
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
///
/// @param windowSize  Analysis window in samples. Recommended: 8192 for guitar
///                    tuning at 48 kHz (covers down to ~C1).
/// @param sampleRate  Audio sample rate in Hz (e.g. 44100.0, 48000.0).
/// @param config      Initial configuration. Use pitchDetectorDefaultConfig()
///                    as a starting point and adjust individual fields.
/// @return Allocated detector, or NULL on failure. Caller must call pitchDetectorDestroy().
PitchDetector* pitchDetectorCreate(int windowSize, float sampleRate, PitchDetectorConfig config);

/// Free all memory associated with a detector.
void pitchDetectorDestroy(PitchDetector* detector);

/// Apply a new configuration. Takes effect on the next analysis cycle.
/// Safe to call at any time between pitchDetectorProcess calls.
void pitchDetectorConfigure(PitchDetector* detector, PitchDetectorConfig config);

/// Feed new audio samples into the detector (streaming).
///
/// Call this from your audio callback or processing loop with each new chunk
/// of mono float PCM samples. The detector accumulates samples in an internal
/// ring buffer and automatically runs analysis every hopSize samples.
///
/// @param detector  Previously created detector
/// @param samples   Pointer to mono float PCM samples (any count)
/// @param count     Number of samples
void pitchDetectorProcess(PitchDetector* detector, const float* samples, int count);

/// Retrieve the latest pitch detection result.
///
/// Returns the most recent stable result after all filtering and hysteresis.
/// If no pitch is detected (silence, noise, or low confidence), hz == -1.
///
/// @param detector  Previously created detector
/// @return Latest pitch result (struct copy, safe to read from any thread)
PitchResult pitchDetectorGetResult(PitchDetector* detector);

#ifdef __cplusplus
}
#endif

#endif /* PITCH_DSP_H */
