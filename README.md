# PitchDSP

A streaming pitch detection library written in C, packaged as a Swift Package for iOS.

## Algorithm

Uses **Harmonic Product Spectrum (HPS)** — a robust technique for detecting the fundamental frequency of harmonic signals (guitar, bass, voice, etc.) even when the fundamental is weak or absent.

Pipeline:

```
Mic samples → Ring Buffer → Power Gate → Hann Window → FFT
→ Mains Hum Suppression → Octave-Band Noise Gate
→ Interpolate (5×) → HPS (×1..×5) → Peak → Note + Cents
```

Original algorithm by [not-chciken](https://github.com/not-chciken/guitar_tuner). Reimplemented in C for streaming use.

## API

```c
// Create a detector (windowSize: 8192 recommended for guitar at 48 kHz)
PitchDetector *d = pitchDetectorCreate(8192, 48000.0f);

// Feed samples from your audio callback
pitchDetectorProcess(d, samples, sampleCount);

// Poll from UI thread (30–60 Hz)
PitchResult r = pitchDetectorGetResult(d);
if (r.hz > 0) {
    // r.hz       — frequency in Hz
    // r.midi     — nearest MIDI note (0–127)
    // r.cents    — deviation from nearest note (–50..+50)
    // r.confidence  — 0.0..1.0 HPS peak prominence
    // r.stability   — 0.0..1.0 note stability score
}

pitchDetectorDestroy(d);
```

**Thread safety:** `pitchDetectorProcess` and `pitchDetectorGetResult` may be called from different threads as long as only one thread calls `Process` at a time.

## Swift Package Manager

```swift
.package(url: "https://github.com/peppperrroni/PitchDSP", from: "1.0.0")
```

```swift
.target(
    name: "YourTarget",
    dependencies: ["PitchDSP"]
)
```

## License

MIT
