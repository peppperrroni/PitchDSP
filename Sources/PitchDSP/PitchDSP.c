// PitchDSP.c — YIN pitch detector
//
// ALGORITHM: YIN (de Cheveigné & Kawahara, 2002)
// ================================================
//   The squared difference function d(τ) = Σ(x[t] − x[t+τ])²
//   expands to d(τ) = 2·(r(0) − r(τ)) where r is the autocorrelation.
//   Autocorrelation is computed via FFT in O(n log n).
//   The Cumulative Mean Normalized Difference Function (CMNDF) is:
//     d'(0)=1,  d'(τ)= d(τ)·τ / Σ_{j=1}^{τ} d(j)
//   The first pit below yinThreshold is found and its bottom is taken
//   as the period estimate. Parabolic interpolation refines sub-sample.
//
// NOTE on octave safety: YIN tends to find the fundamental period, but
// can detect harmonics (especially the 2nd/4th) when they are stronger
// than the fundamental in the autocorrelation. Harmonic correction
// (checking period × 2..5) is done in a higher layer. The raw DSP
// result should be treated as "best YIN candidate", not a guaranteed
// fundamental.
//
// Set DEBUG_PITCH 1 to enable per-analysis console logging.
#define DEBUG_PITCH 0

#include "PitchDSP.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#if DEBUG_PITCH
#include <stdio.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ==========================================================================
//  Internal state
// ==========================================================================

struct PitchDetector {
    int         windowSize;
    int         halfWindow;
    float       sampleRate;
    PitchDetectorConfig config;
    int         hopSize;

    // Ring buffer (circular audio accumulation)
    float*      ringBuffer;
    int         writePos;
    int         samplesSinceAnalysis;

    // Analysis workspace (pre-allocated, reused each hop)
    float*      window;     // windowSize: current analysis frame
    float*      fftReal;    // windowSize: FFT workspace (real)
    float*      fftImag;    // windowSize: FFT workspace (imag)
    float*      cmndf;      // halfWindow+1: CMNDF values for lags 0..halfWindow

    // Latest result (written by process, read by getResult)
    PitchResult result;
};

// ==========================================================================
//  Iterative Cooley-Tukey radix-2 DIT FFT (in-place)
// ==========================================================================

static void bit_reverse_swap(float* re, float* im, int n) {
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        while (j & bit) { j ^= bit; bit >>= 1; }
        j ^= bit;
        if (i < j) {
            float tr = re[i]; re[i] = re[j]; re[j] = tr;
            float ti = im[i]; im[i] = im[j]; im[j] = ti;
        }
    }
}

// sign = -1 → forward FFT, sign = +1 → inverse FFT (with 1/N normalisation)
static void fft_inplace(float* re, float* im, int n, int sign) {
    bit_reverse_swap(re, im, n);
    for (int len = 2; len <= n; len <<= 1) {
        float ang = sign * 2.0f * (float)M_PI / (float)len;
        float wr = cosf(ang), wi = sinf(ang);
        for (int i = 0; i < n; i += len) {
            float cr = 1.0f, ci = 0.0f;
            int half = len >> 1;
            for (int j = 0; j < half; j++) {
                int u = i + j, v = i + j + half;
                float tr = cr * re[v] - ci * im[v];
                float ti = cr * im[v] + ci * re[v];
                re[v] = re[u] - tr;
                im[v] = im[u] - ti;
                re[u] += tr;
                im[u] += ti;
                float new_cr = cr * wr - ci * wi;
                ci = cr * wi + ci * wr;
                cr = new_cr;
            }
        }
    }
    if (sign == 1) {
        float inv = 1.0f / (float)n;
        for (int i = 0; i < n; i++) { re[i] *= inv; im[i] *= inv; }
    }
}

// ==========================================================================
//  Default config
// ==========================================================================

PitchDetectorConfig pitchDetectorDefaultConfig(void) {
    return (PitchDetectorConfig){
        .powerThreshold    = 1e-5f,   // -50 dBFS RMS gate
        .peakThreshold     = 0.005f,  // 0.5% peak gate — stops decaying tail noise
        .yinThreshold      = 0.15f,   // primary CMNDF threshold
        .fallbackThreshold = 0.25f,   // accept global min if no primary pit found
        .hopDivisor        = 8,       // 8192/8 = 1024 samples/hop ≈ 47fps @ 48kHz
    };
}

// ==========================================================================
//  Create / Destroy / Configure
// ==========================================================================

PitchDetector* pitchDetectorCreate(int windowSize, float sampleRate, PitchDetectorConfig config) {
    PitchDetector* d = (PitchDetector*)calloc(1, sizeof(PitchDetector));
    if (!d) return NULL;

    int divisor = (config.hopDivisor > 0) ? config.hopDivisor : 8;

    d->windowSize = windowSize;
    d->halfWindow = windowSize / 2;
    d->sampleRate = sampleRate;
    d->config     = config;
    d->hopSize    = windowSize / divisor;

    d->ringBuffer = (float*)calloc((size_t)windowSize, sizeof(float));
    d->window     = (float*)calloc((size_t)windowSize, sizeof(float));
    d->fftReal    = (float*)calloc((size_t)windowSize, sizeof(float));
    d->fftImag    = (float*)calloc((size_t)windowSize, sizeof(float));
    d->cmndf      = (float*)calloc((size_t)(windowSize / 2 + 1), sizeof(float));

    if (!d->ringBuffer || !d->window || !d->fftReal || !d->fftImag || !d->cmndf) {
        pitchDetectorDestroy(d);
        return NULL;
    }

    d->result.hz         = -1.0f;
    d->result.confidence = 0.0f;
    d->result.sequence   = 0;

    return d;
}

void pitchDetectorDestroy(PitchDetector* d) {
    if (!d) return;
    free(d->ringBuffer);
    free(d->window);
    free(d->fftReal);
    free(d->fftImag);
    free(d->cmndf);
    free(d);
}

void pitchDetectorConfigure(PitchDetector* d, PitchDetectorConfig config) {
    if (!d) return;
    int divisor = (config.hopDivisor > 0) ? config.hopDivisor : 8;
    d->config  = config;
    d->hopSize = d->windowSize / divisor;
}

// ==========================================================================
//  Analysis helpers
// ==========================================================================

// Remove DC offset in-place.
static void remove_dc(float* x, int n) {
    float mean = 0.0f;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= (float)n;
    for (int i = 0; i < n; i++) x[i] -= mean;
}

// Write invalid result and increment sequence.
static void invalidate_result(PitchDetector* d) {
    d->result.hz         = -1.0f;
    d->result.confidence = 0.0f;
    d->result.sequence++;
}

// Measure RMS and peak of window. Fills *rms_out and *peak_out.
static void measure_signal(const float* window, int n, float* rms_out, float* peak_out) {
    float power = 0.0f, peak = 0.0f;
    for (int i = 0; i < n; i++) {
        float s = window[i];
        power += s * s;
        float a = fabsf(s);
        if (a > peak) peak = a;
    }
    *rms_out  = sqrtf(power / (float)n);
    *peak_out = peak;
}

// Find the first CMNDF pit whose minimum is below threshold.
// "First pit" = enter below threshold, walk to the local minimum.
// Returns the tau at the pit bottom, or -1 if none found.
static int find_yin_period(const float* cmndf, int H, float threshold) {
    for (int tau = 2; tau <= H; tau++) {
        if (cmndf[tau] < threshold) {
            // Walk to the bottom of this pit
            while (tau + 1 <= H && cmndf[tau + 1] < cmndf[tau]) {
                tau++;
            }
            return tau;
        }
    }
    return -1;
}

// Fallback: return the tau of the global CMNDF minimum if it is below
// fallbackThreshold. Used when no pit satisfies yinThreshold.
static int find_best_period_fallback(const float* cmndf, int H, float fallbackThreshold) {
    int   bestTau = -1;
    float bestVal = 1.0f;
    for (int tau = 2; tau <= H; tau++) {
        if (cmndf[tau] < bestVal) {
            bestVal = cmndf[tau];
            bestTau = tau;
        }
    }
    return (bestTau >= 0 && bestVal < fallbackThreshold) ? bestTau : -1;
}

// ==========================================================================
//  Analysis (called internally every hopSize samples)
// ==========================================================================

static void run_analysis(PitchDetector* d) {
    const int N = d->windowSize;
    const int H = d->halfWindow;

    // Step 1 — DC removal
    // A DC offset biases the difference function; removing it first gives
    // a cleaner CMNDF, especially important for close-mic recordings.
    remove_dc(d->window, N);

    // Step 2 — RMS + peak gate
    float rms, peak;
    measure_signal(d->window, N, &rms, &peak);

    if (rms < d->config.powerThreshold || peak < d->config.peakThreshold) {
        invalidate_result(d);
        return;
    }

    // Step 3 — Forward FFT of the window
    for (int i = 0; i < N; i++) {
        d->fftReal[i] = d->window[i];
        d->fftImag[i] = 0.0f;
    }
    fft_inplace(d->fftReal, d->fftImag, N, -1);

    // Step 4 — Power spectrum → autocorrelation via IFFT
    // r(τ) = IFFT(|X[k]|²). Note: this is circular autocorrelation (no
    // zero-padding). For large τ (low-frequency notes near B0) the circular
    // wrap-around can introduce small errors; a future improvement would
    // zero-pad to 2N. For guitar range (B2–E6) the effect is negligible.
    for (int i = 0; i < N; i++) {
        d->fftReal[i] = d->fftReal[i] * d->fftReal[i] + d->fftImag[i] * d->fftImag[i];
        d->fftImag[i] = 0.0f;
    }
    fft_inplace(d->fftReal, d->fftImag, N, 1);

    // Step 5 — CMNDF
    // d(τ) = 2·(r(0) − r(τ));  d'(0) = 1;
    // d'(τ) = d(τ)·τ / Σ_{j=1}^{τ} d(j)
    float r0     = d->fftReal[0];
    float cumSum = 0.0f;
    d->cmndf[0]  = 1.0f;

    for (int tau = 1; tau <= H; tau++) {
        float diff = 2.0f * (r0 - d->fftReal[tau]);
        if (diff < 0.0f) diff = 0.0f;   // numerical floor
        cumSum += diff;
        d->cmndf[tau] = (cumSum > 0.0f) ? (diff * (float)tau / cumSum) : 1.0f;
    }

    // Step 6 — Threshold search: first pit below yinThreshold, walk to bottom
    int period = find_yin_period(d->cmndf, H, d->config.yinThreshold);

    // Fallback: global minimum if within fallbackThreshold
    // (helps strings whose CMNDF minimum sits just above yinThreshold)
    if (period < 0) {
        period = find_best_period_fallback(d->cmndf, H, d->config.fallbackThreshold);
    }

    if (period < 0) {
        invalidate_result(d);
        return;
    }

    float beta = d->cmndf[period];  // CMNDF value at found period

    // Step 7 — Parabolic interpolation for sub-sample precision
    float alpha     = (period > 1) ? d->cmndf[period - 1] : beta;
    float gamma     = (period < H) ? d->cmndf[period + 1] : beta;
    float denom     = 2.0f * (2.0f * beta - alpha - gamma);
    float tau_refined = (float)period;

    if (fabsf(denom) > 1e-10f) {
        float offset = (gamma - alpha) / denom;
        if (offset >  0.5f) offset =  0.5f;
        if (offset < -0.5f) offset = -0.5f;
        tau_refined += offset;
    }
    if (tau_refined < 1.0f) tau_refined = (float)period;

    // Step 8 — Output
    float hz         = d->sampleRate / tau_refined;
    float confidence = 1.0f - beta;
    if (confidence < 0.0f) confidence = 0.0f;
    if (confidence > 1.0f) confidence = 1.0f;

    d->result.hz         = hz;
    d->result.confidence = confidence;
    d->result.sequence++;

#if DEBUG_PITCH
    printf("[DSP] tau=%d hz=%.2f cmndf=%.3f conf=%.3f rms=%.5f peak=%.5f\n",
           period, hz, beta, confidence, rms, peak);
#endif
}

// ==========================================================================
//  Process (streaming, called from audio callback)
// ==========================================================================

void pitchDetectorProcess(PitchDetector* d, const float* samples, int count) {
    if (!d || !samples || count <= 0) return;

    const int N  = d->windowSize;
    const int HS = d->hopSize;

    for (int i = 0; i < count; i++) {
        d->ringBuffer[d->writePos] = samples[i];
        d->writePos = (d->writePos + 1) % N;
        d->samplesSinceAnalysis++;

        if (d->samplesSinceAnalysis >= HS) {
            d->samplesSinceAnalysis = 0;

            // Copy most recent N samples from ring buffer (chronological order)
            int start = d->writePos;  // oldest sample after advancing
            for (int j = 0; j < N; j++) {
                d->window[j] = d->ringBuffer[(start + j) % N];
            }

            run_analysis(d);
        }
    }
}

// ==========================================================================
//  GetResult (safe from any thread)
// ==========================================================================

PitchResult pitchDetectorGetResult(PitchDetector* d) {
    if (!d) {
        return (PitchResult){ .hz = -1.0f, .confidence = 0.0f, .sequence = 0 };
    }
    return d->result;
}
