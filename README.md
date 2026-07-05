# PitchDSP

A streaming pitch detection library written in C, packaged as a Swift Package for iOS.

## Algorithm

Uses the **YIN algorithm** (de Cheveigné & Kawahara, 2002) — direct difference method, no FFT.

The **Cumulative Mean Normalized Difference Function (CMNDF)** is computed directly in the time domain. The first CMNDF minimum below the detection threshold is taken as the period, refined with parabolic interpolation for sub-sample accuracy.

Two post-processing corrections handle instrument-specific edge cases:

- **Harmonic correction** — when YIN locks on a sub-harmonic (detects too low), checks whether a shorter period (τ/2..τ/5) has a comparable CMNDF value and prefers the higher frequency (truer fundamental). Applied to the lower third of the valid period range.

- **Sub-fundamental correction** — when a wound string's dominant second harmonic causes YIN to detect one octave too high, checks whether the doubled period has a deeper CMNDF and corrects downward. Applied to the upper third of the valid period range, with a margin requirement to prevent false corrections on well-detected signals.

## API

```c
// 1. Get default config and tune it if needed
PitchDetectorConfig cfg = pitchDetectorDefaultConfig();
cfg.yinThreshold    = 0.20f;   // permissive — works for guitar/bass
cfg.minConfidence   = 0.55f;   // lower for real recordings
cfg.minHz           = 25.0f;   // covers 5-string bass (B0 = 30.87 Hz)

// 2. Create detector (windowSize 8192 recommended)
PitchDetector* d = pitchDetectorCreate(8192, 48000.0f, cfg);

// 3. Feed samples from your audio callback (any chunk size)
pitchDetectorProcess(d, samples, sampleCount);

// 4. Poll from UI thread — check sequence to detect new results
PitchResult r = pitchDetectorGetResult(d);
if (r.sequence != lastSequence) {
    lastSequence = r.sequence;
    if (r.hz > 0) {
        // r.hz         — fundamental frequency in Hz
        // r.confidence — quality 0..1 (depth + local contrast)
    }
}

pitchDetectorDestroy(d);
```

### Configuration

| Field | Default | Description |
|---|---|---|
| `yinThreshold` | 0.15 | Primary CMNDF threshold. Lower = stricter. Range 0.05–0.25. |
| `fallbackThreshold` | 0.25 | Accept global CMNDF minimum if no pit found below `yinThreshold`. |
| `octaveTolerance` | 0.03 | Margin for harmonic correction. 0 disables. |
| `minConfidence` | 0.72 | Frames below this are reported as no-pitch (hz = −1). |
| `hopDivisor` | 8 | Analysis rate = sampleRate / (windowSize / hopDivisor). |
| `minHz` | 25.0 | Lower frequency bound — caps the tau search to avoid CMNDF bias near halfWindow. |
| `powerThreshold` | 1e-5 | RMS gate (~−50 dBFS). |
| `peakThreshold` | 0.005 | Peak amplitude gate for decaying tails. |

**Recommended for guitar/bass tuner:**

```c
cfg.yinThreshold    = 0.20f;
cfg.fallbackThreshold = 0.55f;
cfg.minConfidence   = 0.55f;
cfg.hopDivisor      = 8;
```

## Swift Package Manager

```swift
.package(url: "https://github.com/peppperrroni/PitchDSP", from: "4.0.0")
```

```swift
.target(
    name: "YourTarget",
    dependencies: ["PitchDSP"]
)
```

## Tests

A standalone C test suite lives in `Tests/`. It requires no external dependencies.

```bash
cd Tests
make run        # build + run against Resources/ WAV corpus
make clean      # remove binary
```

**Test suites:**

1. **YIN Validation** — pure synthetic signals (pure tone, noise rejection, octave safety)
2. **Synthetic Corpus** — all guitar + bass frequencies at 44100 Hz using generated sine waves
3. **Harmonic Corpus** — fundamental + equal-amplitude 2nd harmonic (wound string simulation)
4. **WAV Corpus** — real recorded instrument samples at 44100 Hz and 48000 Hz

Acceptance criteria (stable window = frames 14–49):

| Metric | Synthetic | WAV |
|---|---|---|
| First detection | ≤ frame 25 | ≤ frame 25 |
| Detection rate | ≥ 80 % | ≥ 80 % |
| Wrong note rate | ≤ 5 % | ≤ 8 % |
| Cents std dev | ≤ 12 ¢ | ≤ 12 ¢ |
