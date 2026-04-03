// PitchDSP.h — Streaming pitch detector (Harmonic Product Spectrum)
//
// Based on https://github.com/not-chciken/guitar_tuner
// Reimplemented in C for the Nocturne iOS tuner.
//
// Pipeline:
//   Mic samples -> Ring Buffer -> Power Gate -> Hann Window -> FFT
//   -> Mains Hum Suppression -> Octave-Band Noise Gate
//   -> Interpolate (5x) -> HPS (×1..×5) -> Peak -> Note + Cents
//
// Usage:
//   PitchDetector* d = pitchDetectorCreate(8192, 48000.0f);
//
//   // In your audio callback or processing loop:
//   pitchDetectorProcess(d, micSamples, sampleCount);
//
//   // In your UI update (e.g. 30-60 Hz timer):
//   PitchResult r = pitchDetectorGetResult(d);
//   if (r.hz > 0) { /* display note, cents, needle, etc. */ }
//
//   pitchDetectorDestroy(d);
//
// Thread safety: pitchDetectorProcess and pitchDetectorGetResult may be called
// from different threads as long as only ONE thread calls Process at a time.
// GetResult reads a plain struct copy — no lock needed on most architectures.

#ifndef PITCH_DSP_H
#define PITCH_DSP_H

#ifdef __cplusplus
extern "C" {
#endif

// --- Result ---

typedef struct {
    float hz;           // Detected pitch in Hz. -1 if no pitch detected.
    float cents;        // Deviation from nearest note in cents (-50 .. +50).
    int   midi;         // Nearest MIDI note number (0-127). -1 if no pitch.
    float confidence;   // 0.0 .. 1.0 based on HPS peak prominence.
    float stability;    // 0.0 .. 1.0 note stability score.
} PitchResult;

// --- Opaque detector ---

typedef struct PitchDetector PitchDetector;

// --- API ---

/// Create a streaming pitch detector.
///
/// @param windowSize  Analysis window in samples. Recommended: 8192 for guitar tuning
///                    at 48 kHz (covers down to ~C1). Use 4096 for higher-pitched instruments.
/// @param sampleRate  Audio sample rate in Hz (e.g. 44100.0, 48000.0).
/// @return Allocated detector, or NULL on failure. Caller must call pitchDetectorDestroy().
PitchDetector* pitchDetectorCreate(int windowSize, float sampleRate);

/// Free all memory associated with a detector.
void pitchDetectorDestroy(PitchDetector* detector);

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
/// This returns the most recent stable result after all filtering and hysteresis.
/// If no pitch is detected (silence, noise, or low confidence), hz == -1 and midi == -1.
///
/// @param detector  Previously created detector
/// @return Latest pitch result (struct copy, safe to read from any thread)
PitchResult pitchDetectorGetResult(PitchDetector* detector);

#ifdef __cplusplus
}
#endif

#endif /* PITCH_DSP_H */
