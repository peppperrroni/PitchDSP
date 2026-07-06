# PitchDSP

A streaming, monophonic pitch detector in dependency-free C11. Two-stage YIN
(decimated coarse search + full-rate refine) with harmonic-aware corrections,
tuned out of the box for instrument tuners spanning B0–A4. Packaged as a
Swift Package for direct use from iOS/macOS apps.

## Quick start

```c
PitchDetectorConfig cfg = pitchDetectorDefaultConfig();
PitchDetector* d = pitchDetectorCreate(8192, 48000.0f, cfg);

// In your audio callback (any chunk size):
pitchDetectorProcess(d, micSamples, sampleCount);
PitchResult results[16];
int n = pitchDetectorDrainResults(d, results, 16);
for (int i = 0; i < n; i++) { /* results[i].hz > 0 = pitch */ }

pitchDetectorDestroy(d);
```

`windowSize` (8192 above) must be a multiple of 4 (the internal decimation
factor); other values return NULL from `pitchDetectorCreate`.

## Algorithm

The detector runs YIN (de Cheveigné & Kawahara, 2002) in two stages rather
than once at full rate. Incoming audio is first passed through a 47-tap
Hamming-windowed sinc anti-alias filter and decimated 4×; the coarse search
below runs on this quarter-rate signal, which cuts the O(N·τ) direct
difference-function cost roughly 16× without losing low-end accuracy.

Coarse YIN then computes the squared difference function and its Cumulative
Mean Normalized Difference Function (CMNDF) directly (no FFT), searching only
the τ range implied by `[minHz, maxHz]`. It looks for the first CMNDF dip
below `yinThreshold`; if none qualifies, it falls back to the single deepest
dip in range provided that dip is below `fallbackThreshold`.

The coarse period candidate then passes through two corrections: a harmonic
(sub-harmonic-lock) correction that prefers a shorter period τ/2..τ/5 when its
CMNDF is comparably deep (fixes wound strings that make YIN settle an octave
low), and a sub-fundamental correction that prefers the doubled period when
it is meaningfully deeper (fixes locking onto a dominant 2nd harmonic).

Once a period survives the confidence gate, it is refined at full sample rate:
the raw (non-normalized) difference function is recomputed in a narrow band
around the coarse period and parabolic interpolation locates the sub-sample
minimum. Confidence = 0.75·(CMNDF pit depth) + 0.25·(τ-proportional local
contrast × 8), clamped to 0..1.

## Configuration

All fields are on `PitchDetectorConfig`; get defaults via
`pitchDetectorDefaultConfig()`.

| Field | Default | Guidance |
|---|---|---|
| `powerThreshold` | `1e-5` | RMS gate (~-50 dBFS) below which a window is skipped. |
| `peakThreshold` | `0.005` | Peak-amplitude gate; catches decaying tails where RMS still passes but the waveform is too quiet for a reliable period. |
| `yinThreshold` | `0.20` | Primary CMNDF pit threshold. Lower = stricter. Range 0.05 (strict) – 0.25 (permissive). |
| `fallbackThreshold` | `0.30` | Accept the global CMNDF minimum if no dip beats `yinThreshold`. |
| `octaveTolerance` | `0.03` | Harmonic-correction margin: prefer period/N (N=2..5) if its CMNDF is within this of the detected period's. `0.0` disables. Lower it if spurious octave-up jumps appear. |
| `minConfidence` | `0.60` | Frames below this confidence are reported invalid (`hz = -1`). |
| `hopDivisor` | `8` | Analysis rate = sampleRate / (windowSize / hopDivisor); default gives ~47 fps at 48 kHz with windowSize 8192. |
| `minHz` | `25.0` | Lower frequency bound; sets the τ search ceiling. Covers 5-string bass B0 (30.87 Hz). Raise to 60-70 Hz for guitar-only apps to tighten the range. |
| `maxHz` | `1000.0` | Upper frequency bound; sets the τ search floor. A fundamental above this reports as a subharmonic within range, by design. |

## Threading contract

`pitchDetectorProcess`, `pitchDetectorDrainResults`, and
`pitchDetectorGetResult` must all be called from a single thread — typically
the audio callback thread. There is no internal locking. If you need results
on another thread (e.g. UI), bridge them yourself (queue, atomic, etc.).

## Testing

```bash
cd Tests && make run
```

This builds and runs two binaries: `test_internal`, a white-box suite that
reaches into the FIR anti-alias filter and 4× decimator directly, and
`test_pitchdsp`, the main suite covering synthetic signals, a synthetic
frequency/harmonic corpus, a cents-accuracy sweep, and a WAV corpus of real
recorded instrument samples. The WAV metrics separate non-octave wrong notes
(a real DSP bug, gated strictly) from exact-octave misses (a display-layer
concern, gated more loosely) — some bass WAV fixtures in `Tests/Resources`
intentionally contain content one octave below their filename (e.g.
`bass_A1.wav` holds A0 = 27.5 Hz), documented in comments at the corpus
definitions in `test_pitchdsp.c`.

## Integration

```swift
.package(url: "https://github.com/peppperrroni/PitchDSP", exact: "1.0.0")
```

```swift
.target(
    name: "YourTarget",
    dependencies: ["PitchDSP"]
)
```
