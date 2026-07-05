// PitchDSP.c — YIN pitch detector
//
// ALGORITHM: YIN (de Cheveigné & Kawahara, 2002)
// ================================================
//   The squared difference function d(τ) = Σ(x[t] − x[t+τ])²
//   expands to d(τ) = 2·(r(0) − r(τ)) where r is the autocorrelation.
//   Autocorrelation is computed via FFT in O(n log n).
//   The Cumulative Mean Normalized Difference Function (CMNDF) is:
//     d'(0)=1,  d'(τ)= d(τ)·τ / Σ_{j=1}^{τ} d(j)
//   The first τ with d'(τ) < threshold and d'(τ) < d'(τ-1) is the period.
//   Parabolic interpolation refines the sub-sample estimate.
//   confidence = 1 − d'(τ_best): 1.0 = clean, 0.0 = noise.
//
// Octave safety: taking the *first* sub-threshold minimum always returns
// the fundamental period — never a sub-multiple (which would be a shorter
// period at a higher lag-index, never found first).

#include "PitchDSP.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
        .powerThreshold = 1e-5f,   // -50 dBFS
        .yinThreshold   = 0.15f,   // industry standard starting point
        .hopDivisor     = 8,       // 8192/8 = 1024 samples/hop ≈ 47fps @ 48kHz
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
//  Analysis (called internally every hopSize samples)
// ==========================================================================

static void run_analysis(PitchDetector* d) {
    const int N = d->windowSize;
    const int H = d->halfWindow;

    // Step 1 — Power gate: RMS of current window
    float power = 0.0f;
    for (int i = 0; i < N; i++) {
        float s = d->window[i];
        power += s * s;
    }
    float rms = sqrtf(power / (float)N);

    if (rms < d->config.powerThreshold) {
        d->result.hz         = -1.0f;
        d->result.confidence = 0.0f;
        d->result.sequence++;
        return;
    }

    // Step 2 — Forward FFT of the window
    // Load window into FFT buffers
    for (int i = 0; i < N; i++) {
        d->fftReal[i] = d->window[i];
        d->fftImag[i] = 0.0f;
    }
    fft_inplace(d->fftReal, d->fftImag, N, -1);

    // Step 3 — Power spectrum: |X[k]|² (autocorrelation via IFFT of power spectrum)
    for (int i = 0; i < N; i++) {
        d->fftReal[i] = d->fftReal[i] * d->fftReal[i] + d->fftImag[i] * d->fftImag[i];
        d->fftImag[i] = 0.0f;
    }
    fft_inplace(d->fftReal, d->fftImag, N, 1);
    // fftReal[tau] = r(tau) = autocorrelation at lag tau

    // Step 4 — CMNDF
    // d(tau) = 2*(r(0) - r(tau))
    // d'(0) = 1
    // d'(tau) = d(tau) * tau / sum_{j=1}^{tau} d(j)
    float r0     = d->fftReal[0];
    float cumSum = 0.0f;
    d->cmndf[0]  = 1.0f;

    for (int tau = 1; tau <= H; tau++) {
        float diff = 2.0f * (r0 - d->fftReal[tau]);
        if (diff < 0.0f) diff = 0.0f;   // numerical floor
        cumSum += diff;
        d->cmndf[tau] = (cumSum > 0.0f) ? (diff * (float)tau / cumSum) : 1.0f;
    }

    // Step 5 — Threshold search: first local minimum below yinThreshold
    // "first tau where d'(tau) < threshold AND d'(tau) < d'(tau-1)"
    int period = -1;
    for (int tau = 2; tau <= H; tau++) {
        if (d->cmndf[tau] < d->config.yinThreshold &&
            d->cmndf[tau] < d->cmndf[tau - 1]) {
            period = tau;
            break;
        }
    }

    if (period < 0) {
        d->result.hz         = -1.0f;
        d->result.confidence = 0.0f;
        d->result.sequence++;
        return;
    }

    // Step 6 — Parabolic interpolation for sub-sample precision
    float alpha = (period > 1) ? d->cmndf[period - 1] : d->cmndf[period];
    float beta  = d->cmndf[period];
    float gamma = (period < H) ? d->cmndf[period + 1] : d->cmndf[period];

    float denom        = 2.0f * (2.0f * beta - alpha - gamma);
    float tau_refined  = (float)period;
    if (fabsf(denom) > 1e-10f) {
        float offset = (gamma - alpha) / denom;
        // Clamp offset to ±0.5 to avoid extreme refinement
        if (offset >  0.5f) offset =  0.5f;
        if (offset < -0.5f) offset = -0.5f;
        tau_refined += offset;
    }
    if (tau_refined < 1.0f) tau_refined = (float)period;

    // Step 7 — Output
    d->result.hz         = d->sampleRate / tau_refined;
    d->result.confidence = 1.0f - beta;
    if (d->result.confidence < 0.0f) d->result.confidence = 0.0f;
    if (d->result.confidence > 1.0f) d->result.confidence = 1.0f;
    d->result.sequence++;
}

// ==========================================================================
//  Process (streaming, called from audio callback)
// ==========================================================================

void pitchDetectorProcess(PitchDetector* d, const float* samples, int count) {
    if (!d || !samples || count <= 0) return;

    const int N = d->windowSize;
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
