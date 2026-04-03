// PitchDSP.c — Harmonic Product Spectrum pitch detector
//
// Based on https://github.com/not-chciken/guitar_tuner by not-chciken.
// Reimplemented in plain C for the Nocturne iOS tuner.
//
// ALGORITHM: Harmonic Product Spectrum (HPS)
// ===========================================
//
//   Mic samples
//       │
//       ▼
//   ┌─────────────┐
//   │ Ring Buffer │  Sliding window, triggers every hopSize samples.
//   └──────┬──────┘
//          │
//          ▼
//   ┌─────────────┐
//   │ Power Gate  │  Mean squared amplitude > 1e-6.
//   └──────┬──────┘
//          │
//          ▼
//   ┌─────────────┐
//   │ Hann Window │  Reduce spectral leakage.
//   └──────┬──────┘
//          │
//          ▼
//   ┌─────────────┐
//   │    FFT      │  Radix-2 Cooley-Tukey, magnitude spectrum.
//   └──────┬──────┘
//          │
//          ▼
//   ┌─────────────────────┐
//   │ Mains Hum Suppressio│  Zero bins below 62 Hz.
//   └──────┬──────────────┘
//          │
//          ▼
//   ┌─────────────────────┐
//   │ Octave-Band Noise   │  Per-band RMS gate (9 octave bands).
//   │ Suppression         │  Threshold = 0.2 × band RMS.
//   └──────┬──────────────┘
//          │
//          ▼
//   ┌─────────────────────┐
//   │ Interpolate (5×)    │  Linear upsample for sub-bin resolution.
//   │ + L2 Normalize      │
//   └──────┬──────────────┘
//          │
//          ▼
//   ┌─────────────────────┐
//   │ HPS                 │  Multiply spectrum with 2×, 3×, 4×, 5×
//   │                     │  decimated copies. Fundamental peaks align.
//   └──────┬──────────────┘
//          │
//          ▼
//   ┌─────────────────────┐
//   │ Peak → Note + Cents │  argmax → Hz → MIDI → cents.
//   │ + Note Buffer       │  2-deep stability filter.
//   └─────────────────────┘

#include "PitchDSP.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ==========================================================================
//  Configuration (matching reference implementation)
// ==========================================================================

#define NUM_HPS            5       // Number of harmonic product spectrums
#define POWER_THRESH       1e-6f   // Minimum signal power to activate
#define WHITE_NOISE_THRESH 0.2f    // Octave-band noise gate fraction
#define MAINS_HUM_FREQ     62.0f  // Zero all bins below this frequency
#define CONCERT_PITCH      440.0f // A4 reference
#define HOP_DIVISOR        4      // Analyze every window/4 samples
#define NOTE_BUFFER_SIZE   2      // Consecutive matches needed to lock

#define NUM_OCTAVE_BANDS   10
static const float OCTAVE_BANDS[NUM_OCTAVE_BANDS] = {
    50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600
};

// ==========================================================================
//  FFT — Radix-2 Cooley-Tukey (in-place, forward transform)
// ==========================================================================

static int next_pow2(int n) {
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

static void fft_forward(float *re, float *im, int n) {
    // Bit-reverse permutation
    int j = 0;
    for (int i = 0; i < n - 1; i++) {
        if (i < j) {
            float t;
            t = re[j]; re[j] = re[i]; re[i] = t;
            t = im[j]; im[j] = im[i]; im[i] = t;
        }
        int k = n >> 1;
        while (k <= j) { j -= k; k >>= 1; }
        j += k;
    }

    // Cooley-Tukey butterflies
    for (int step = 2; step <= n; step <<= 1) {
        float angle = -2.0f * (float)M_PI / (float)step;
        float wR = cosf(angle), wI = sinf(angle);
        int half = step >> 1;
        for (int group = 0; group < n; group += step) {
            float cR = 1.0f, cI = 0.0f;
            for (int pair = 0; pair < half; pair++) {
                int a = group + pair;
                int b = a + half;
                float tR = cR * re[b] - cI * im[b];
                float tI = cR * im[b] + cI * re[b];
                re[b] = re[a] - tR;
                im[b] = im[a] - tI;
                re[a] += tR;
                im[a] += tI;
                float nR = cR * wR - cI * wI;
                cI = cR * wI + cI * wR;
                cR = nR;
            }
        }
    }
}

// ==========================================================================
//  Helpers
// ==========================================================================

static int freq_to_midi(float hz) {
    if (hz <= 0.0f) return -1;
    return (int)roundf(69.0f + 12.0f * log2f(hz / CONCERT_PITCH));
}

static float midi_to_freq(int midi) {
    return CONCERT_PITCH * powf(2.0f, (float)(midi - 69) / 12.0f);
}

static float compute_cents(float hz, int midi) {
    float noteHz = midi_to_freq(midi);
    if (noteHz <= 0.0f || hz <= 0.0f) return 0.0f;
    float c = 1200.0f * log2f(hz / noteHz);
    if (c < -50.0f) c = -50.0f;
    if (c >  50.0f) c =  50.0f;
    return c;
}

// ==========================================================================
//  PitchDetector context
// ==========================================================================

struct PitchDetector {
    float sampleRate;
    int   fftSize;        // power of 2
    int   halfFFT;        // fftSize / 2
    int   hopSize;
    float deltaFreq;      // sampleRate / fftSize

    // Ring buffer
    float *ringBuffer;
    int    ringCapacity;
    int    writePos;
    int    samplesInBuffer;
    int    samplesSinceAnalysis;

    // Pre-allocated buffers (all allocated in create, none in process)
    float *hannWindow;     // [fftSize]
    float *fftReal;        // [fftSize]
    float *fftImag;        // [fftSize]
    float *magnitudeSpec;  // [halfFFT]
    float *magSpecInterp;  // [halfFFT * NUM_HPS]
    float *hpsSpec;        // [halfFFT * NUM_HPS]

    // Note stability (2-deep buffer, matching reference)
    int    noteBuffer[NOTE_BUFFER_SIZE];
    float  smoothedCents;

    // Output
    PitchResult result;
};

// ==========================================================================
//  Public API
// ==========================================================================

PitchDetector* pitchDetectorCreate(int windowSize, float sampleRate) {
    if (windowSize < 64 || sampleRate <= 0.0f) return NULL;

    PitchDetector *d = (PitchDetector *)calloc(1, sizeof(PitchDetector));
    if (!d) return NULL;

    d->sampleRate = sampleRate;
    d->fftSize    = next_pow2(windowSize);
    d->halfFFT    = d->fftSize / 2;
    d->deltaFreq  = sampleRate / (float)d->fftSize;
    d->hopSize    = d->fftSize / HOP_DIVISOR;
    if (d->hopSize < 1) d->hopSize = 1;

    // Ring buffer (2× window for safe circular access)
    d->ringCapacity = d->fftSize * 2;
    d->ringBuffer   = (float *)calloc((size_t)d->ringCapacity, sizeof(float));

    // Hann window: 0.5 * (1 - cos(2πi / (N-1)))
    d->hannWindow = (float *)malloc((size_t)d->fftSize * sizeof(float));
    if (d->hannWindow) {
        for (int i = 0; i < d->fftSize; i++)
            d->hannWindow[i] = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * (float)i / (float)(d->fftSize - 1)));
    }

    // FFT buffers
    d->fftReal       = (float *)calloc((size_t)d->fftSize, sizeof(float));
    d->fftImag       = (float *)calloc((size_t)d->fftSize, sizeof(float));
    d->magnitudeSpec = (float *)calloc((size_t)d->halfFFT, sizeof(float));

    // HPS buffers (interpolated size = halfFFT × NUM_HPS)
    int interpSize    = d->halfFFT * NUM_HPS;
    d->magSpecInterp = (float *)calloc((size_t)interpSize, sizeof(float));
    d->hpsSpec       = (float *)calloc((size_t)interpSize, sizeof(float));

    // Verify all allocations
    if (!d->ringBuffer || !d->hannWindow || !d->fftReal || !d->fftImag ||
        !d->magnitudeSpec || !d->magSpecInterp || !d->hpsSpec) {
        pitchDetectorDestroy(d);
        return NULL;
    }

    // Initialize note buffer to invalid
    for (int i = 0; i < NOTE_BUFFER_SIZE; i++)
        d->noteBuffer[i] = -1;

    d->result.hz   = -1.0f;
    d->result.midi = -1;

    return d;
}

void pitchDetectorDestroy(PitchDetector *d) {
    if (!d) return;
    free(d->ringBuffer);
    free(d->hannWindow);
    free(d->fftReal);
    free(d->fftImag);
    free(d->magnitudeSpec);
    free(d->magSpecInterp);
    free(d->hpsSpec);
    free(d);
}

void pitchDetectorProcess(PitchDetector *d, const float *samples, int count) {
    if (!d || !samples || count <= 0) return;

    // =================================================================
    //  1. Ring buffer accumulation
    // =================================================================

    const float *src = samples;
    int toWrite = count;
    if (toWrite > d->ringCapacity) {
        src += (toWrite - d->ringCapacity);
        toWrite = d->ringCapacity;
    }

    for (int i = 0; i < toWrite; i++) {
        d->ringBuffer[d->writePos] = src[i];
        d->writePos = (d->writePos + 1) % d->ringCapacity;
    }

    d->samplesInBuffer += toWrite;
    if (d->samplesInBuffer > d->ringCapacity)
        d->samplesInBuffer = d->ringCapacity;
    d->samplesSinceAnalysis += count;

    if (d->samplesInBuffer < d->fftSize) return;
    if (d->samplesSinceAnalysis < d->hopSize) return;
    d->samplesSinceAnalysis = 0;

    // =================================================================
    //  2. Extract window start position
    // =================================================================

    int start = d->writePos - d->fftSize;
    if (start < 0) start += d->ringCapacity;

    // =================================================================
    //  3. Power check (on raw samples, before windowing)
    // =================================================================

    float power = 0.0f;
    for (int i = 0; i < d->fftSize; i++) {
        int idx = (start + i) % d->ringCapacity;
        float s = d->ringBuffer[idx];
        power += s * s;
    }
    power /= (float)d->fftSize;

    if (power < POWER_THRESH) {
        d->result.hz         = -1.0f;
        d->result.cents      =  0.0f;
        d->result.midi       = -1;
        d->result.confidence =  0.0f;
        d->result.stability  =  0.0f;
        for (int i = 0; i < NOTE_BUFFER_SIZE; i++)
            d->noteBuffer[i] = -1;
        d->smoothedCents = 0.0f;
        return;
    }

    // =================================================================
    //  4. Hann window → FFT input
    // =================================================================

    for (int i = 0; i < d->fftSize; i++) {
        int idx = (start + i) % d->ringCapacity;
        d->fftReal[i] = d->ringBuffer[idx] * d->hannWindow[i];
        d->fftImag[i] = 0.0f;
    }

    // =================================================================
    //  5. FFT (radix-2 Cooley-Tukey)
    // =================================================================

    fft_forward(d->fftReal, d->fftImag, d->fftSize);

    // =================================================================
    //  6. Magnitude spectrum (positive frequencies only)
    // =================================================================

    for (int i = 0; i < d->halfFFT; i++) {
        float re = d->fftReal[i];
        float im = d->fftImag[i];
        d->magnitudeSpec[i] = sqrtf(re * re + im * im);
    }

    // =================================================================
    //  7. Mains hum suppression — zero bins below 62 Hz
    // =================================================================

    int humBins = (int)(MAINS_HUM_FREQ / d->deltaFreq);
    if (humBins > d->halfFFT) humBins = d->halfFFT;
    for (int i = 0; i < humBins; i++)
        d->magnitudeSpec[i] = 0.0f;

    // =================================================================
    //  8. Octave-band white noise suppression
    //     For each band: compute RMS, zero bins < 0.2 × RMS
    // =================================================================

    for (int b = 0; b < NUM_OCTAVE_BANDS - 1; b++) {
        int iStart = (int)(OCTAVE_BANDS[b] / d->deltaFreq);
        int iEnd   = (int)(OCTAVE_BANDS[b + 1] / d->deltaFreq);
        if (iEnd > d->halfFFT) iEnd = d->halfFFT;
        if (iStart >= iEnd) continue;

        // RMS magnitude in this band
        float sumSq = 0.0f;
        int bandSize = iEnd - iStart;
        for (int i = iStart; i < iEnd; i++)
            sumSq += d->magnitudeSpec[i] * d->magnitudeSpec[i];
        float rmsEnergy = sqrtf(sumSq / (float)bandSize);

        // Gate: zero bins below threshold
        float thresh = WHITE_NOISE_THRESH * rmsEnergy;
        for (int i = iStart; i < iEnd; i++) {
            if (d->magnitudeSpec[i] < thresh)
                d->magnitudeSpec[i] = 0.0f;
        }
    }

    // =================================================================
    //  9. Interpolate (upsample ×NUM_HPS) + L2 normalize
    // =================================================================

    int interpSize = d->halfFFT * NUM_HPS;

    for (int i = 0; i < interpSize; i++) {
        float fidx = (float)i / (float)NUM_HPS;
        int lo = (int)fidx;
        if (lo >= d->halfFFT - 1) {
            d->magSpecInterp[i] = d->magnitudeSpec[d->halfFFT - 1];
        } else {
            float frac = fidx - (float)lo;
            d->magSpecInterp[i] = d->magnitudeSpec[lo] * (1.0f - frac)
                                + d->magnitudeSpec[lo + 1] * frac;
        }
    }

    // L2 normalize
    float norm = 0.0f;
    for (int i = 0; i < interpSize; i++)
        norm += d->magSpecInterp[i] * d->magSpecInterp[i];
    norm = sqrtf(norm);

    if (norm < 1e-10f) return;  // all zeros — nothing to detect

    float invNorm = 1.0f / norm;
    for (int i = 0; i < interpSize; i++)
        d->magSpecInterp[i] *= invNorm;

    // =================================================================
    //  10. Harmonic Product Spectrum
    //      Multiply spectrum with decimated copies of itself.
    //      Harmonics at 2F, 3F, 4F, 5F align with fundamental → sharp peak.
    // =================================================================

    memcpy(d->hpsSpec, d->magSpecInterp, (size_t)interpSize * sizeof(float));
    int hpsLen = interpSize;

    for (int h = 0; h < NUM_HPS; h++) {
        int factor = h + 1;
        int newLen = (interpSize + factor - 1) / factor;  // ceil division
        if (newLen > hpsLen) newLen = hpsLen;

        int allZero = 1;
        for (int i = 0; i < newLen; i++) {
            int srcIdx = i * factor;
            if (srcIdx < interpSize)
                d->hpsSpec[i] *= d->magSpecInterp[srcIdx];
            if (d->hpsSpec[i] > 0.0f) allZero = 0;
        }
        hpsLen = newLen;
        if (allZero) break;
    }

    // =================================================================
    //  11. Find peak → frequency
    // =================================================================

    int maxIdx = 0;
    float maxVal = 0.0f;
    for (int i = 0; i < hpsLen; i++) {
        if (d->hpsSpec[i] > maxVal) {
            maxVal = d->hpsSpec[i];
            maxIdx = i;
        }
    }

    if (maxVal <= 0.0f) return;

    // Convert interpolated index back to Hz
    float maxFreq = (float)maxIdx * d->deltaFreq / (float)NUM_HPS;

    // Frequency range guard
    if (maxFreq < 27.0f || maxFreq > 4200.0f) return;

    // =================================================================
    //  12. MIDI note + cents
    // =================================================================

    int midi = freq_to_midi(maxFreq);
    if (midi < 0 || midi > 127) return;
    float cents = compute_cents(maxFreq, midi);

    // =================================================================
    //  13. Note stability — 2-deep buffer (matching reference)
    //      Same note must appear in 2 consecutive analyses to lock.
    // =================================================================

    // Shift buffer
    for (int i = NOTE_BUFFER_SIZE - 1; i > 0; i--)
        d->noteBuffer[i] = d->noteBuffer[i - 1];
    d->noteBuffer[0] = midi;

    // Check if all entries match
    int stable = 1;
    for (int i = 1; i < NOTE_BUFFER_SIZE; i++) {
        if (d->noteBuffer[i] != d->noteBuffer[0]) {
            stable = 0;
            break;
        }
    }

    if (stable) {
        // Two consecutive frames agree → update result
        if (midi != d->result.midi) {
            d->smoothedCents = cents;  // reset on note change
        } else {
            d->smoothedCents = d->smoothedCents * 0.7f + cents * 0.3f;
        }
        d->result.hz         = maxFreq;
        d->result.cents      = d->smoothedCents;
        d->result.midi       = midi;
        d->result.confidence = 1.0f;
        d->result.stability  = 1.0f;
    }
    // If not stable, keep previous result (prevents flicker)
}

PitchResult pitchDetectorGetResult(PitchDetector *d) {
    if (!d) {
        PitchResult empty = { -1.0f, 0.0f, -1, 0.0f, 0.0f };
        return empty;
    }
    return d->result;
}
