// PitchDSP.c — YIN pitch detector (direct difference method)
//
// ALGORITHM: YIN (de Cheveigné & Kawahara, 2002)
// ================================================
//   The squared difference function d(τ) = Σ_{i=0}^{N-τ-1} (x[i] − x[i+τ])²
//   is computed directly — no FFT, no circular autocorrelation.
//   The Cumulative Mean Normalized Difference Function (CMNDF) is then:
//     d'(0)=1,  d'(τ) = d(τ)·τ / Σ_{j=1}^{τ} d(j)
//   The first pit below yinThreshold is found and its bottom taken as the
//   period estimate. Parabolic interpolation refines sub-sample accuracy.
//
// NOTE on candidate quality: YIN returns a plausible period candidate, not
// a guaranteed fundamental. Sub-harmonic locking can occur where YIN finds a
// longer period (lower frequency) than the true fundamental because its CMNDF
// dip is marginally deeper (common on wound bass/guitar strings). A post-YIN
// correction step checks period/2..period/5: if a shorter period's CMNDF is
// within octaveTolerance of the current candidate's CMNDF, the shorter period
// (higher frequency, truer fundamental) is preferred.
// This is separate from the confidence model which gates on pit depth and
// local contrast so that flat or noisy CMNDF regions are rejected.
//
// Set DEBUG_PITCH 1 to enable per-analysis console logging.
// Can be overridden from the compiler command line: -DDEBUG_PITCH=1
#ifndef DEBUG_PITCH
#define DEBUG_PITCH 0
#endif

// Capacity of the internal result queue drained by pitchDetectorDrainResults.
#define PITCH_RESULT_FIFO 16

// Decimation parameters
#define PITCH_DECIM_FACTOR 4
#define PITCH_FIR_TAPS     47

#include "PitchDSP.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#if DEBUG_PITCH
#include <stdio.h>
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
    float*      window;     // windowSize samples: full-rate frame (gate + refine)
    float*      decWindow;  // decWindowSize samples: decimated frame (coarse YIN input)
    float*      diff;       // decWindowSize/2+1: YIN squared difference function d(τ) (decimated)
    float*      cmndf;      // decWindowSize/2+1: CMNDF values (decimated)

    // Decimation (coarse-stage input): anti-alias FIR + keep every 4th sample
    float       firCoef[PITCH_FIR_TAPS];
    float       firHist[PITCH_FIR_TAPS];
    int         firHistPos;
    int         decimPhase;
    float*      decRing;         // decimated ring, decWindowSize entries
    int         decWritePos;
    int         decWindowSize;   // windowSize / PITCH_DECIM_FACTOR
    float       decRate;         // sampleRate / PITCH_DECIM_FACTOR
    int         maxTauD;         // decRate / minHz — upper bound for coarse tau search
    int         minTauD;         // decRate / maxHz — lower bound for coarse tau search
    int         maxTauPlusD;     // min(maxTauD + 1, decWindowSize / 2) — diff/cmndf loop bound

    // Latest result (written by process, read by getResult)
    PitchResult result;

    // Result FIFO (drain API)
    PitchResult fifo[PITCH_RESULT_FIFO];
    int         fifoRead;
    int         fifoWrite;
    int         fifoCount;
};

// ==========================================================================
//  Default config
// ==========================================================================

PitchDetectorConfig pitchDetectorDefaultConfig(void) {
    return (PitchDetectorConfig){
        .powerThreshold    = 1e-5f,   // -50 dBFS RMS gate
        .peakThreshold     = 0.005f,  // 0.5% peak gate — stops decaying tail noise
        .yinThreshold      = 0.15f,   // primary CMNDF pit threshold
        .fallbackThreshold = 0.25f,   // accept global min if no pit found
        .octaveTolerance   = 0.03f,   // harmonic correction: prefer period/N if
                                      // CMNDF[period/N] < CMNDF[period] + tolerance
        .minConfidence     = 0.72f,   // gate: reject if depth+contrast confidence < this
        .hopDivisor        = 8,       // analysis rate = sampleRate / (windowSize/8)
        .minHz             = 25.0f,   // lower frequency limit → maxTau = sampleRate/minHz
        .maxHz             = 1000.0f, // upper frequency limit → minTau = sampleRate/maxHz
    };
}

// ==========================================================================
//  Create / Destroy / Configure
// ==========================================================================

// Compute the tau search ceiling from minHz, clamped to halfWindow.
// Defaults to 25 Hz if minHz is zero or negative (covers 5-string bass B0).
static int compute_max_tau(float sampleRate, float minHz, int halfWindow) {
    float hz = (minHz > 0.0f) ? minHz : 25.0f;
    int   mt = (int)(sampleRate / hz);
    return (mt < halfWindow) ? mt : halfWindow;
}

// Compute the tau search floor from maxHz, clamped to >= 2.
static int compute_min_tau(float sampleRate, float maxHz) {
    float hz = (maxHz > 0.0f) ? maxHz : 1000.0f;
    int   mt = (int)(sampleRate / hz);
    return (mt < 2) ? 2 : mt;
}

// Hamming-windowed sinc low-pass, unity DC gain.
// fcNorm = cutoff / inputSampleRate. With 47 taps and fcNorm 0.085 the stopband
// begins near 0.12×fs — just inside the post-decimation Nyquist (0.125×fs) —
// giving ≥ ~50 dB alias rejection while passing fundamentals + low harmonics
// (flat to ~2.3 kHz at 48 kHz input).
static void build_fir_lowpass(float* h, int taps, float fcNorm) {
    int M = taps - 1;
    float sum = 0.0f;
    for (int i = 0; i < taps; i++) {
        float x = (float)i - (float)M / 2.0f;
        float sinc = (fabsf(x) < 1e-6f)
            ? 2.0f * fcNorm
            : sinf(2.0f * (float)M_PI * fcNorm * x) / ((float)M_PI * x);
        float w = 0.54f - 0.46f * cosf(2.0f * (float)M_PI * (float)i / (float)M);
        h[i] = sinc * w;
        sum += h[i];
    }
    for (int i = 0; i < taps; i++) h[i] /= sum;
}

PitchDetector* pitchDetectorCreate(int windowSize, float sampleRate, PitchDetectorConfig config) {
    PitchDetector* d = (PitchDetector*)calloc(1, sizeof(PitchDetector));
    if (!d) return NULL;

    // Reject window sizes not divisible by decimation factor
    if (windowSize <= 0 || (windowSize % PITCH_DECIM_FACTOR) != 0) { free(d); return NULL; }

    int divisor = (config.hopDivisor > 0) ? config.hopDivisor : 8;

    d->windowSize = windowSize;
    d->halfWindow = windowSize / 2;
    d->sampleRate = sampleRate;
    d->config     = config;
    d->hopSize    = windowSize / divisor;

    // Decimation state
    d->decWindowSize = windowSize / PITCH_DECIM_FACTOR;
    d->decRate       = sampleRate / (float)PITCH_DECIM_FACTOR;
    d->maxTauD       = compute_max_tau(d->decRate, config.minHz, d->decWindowSize / 2);
    d->minTauD       = compute_min_tau(d->decRate, config.maxHz);
    d->maxTauPlusD   = (d->maxTauD + 1 < d->decWindowSize / 2) ? d->maxTauD + 1 : d->decWindowSize / 2;

    d->ringBuffer = (float*)calloc((size_t)windowSize, sizeof(float));
    d->window     = (float*)calloc((size_t)windowSize, sizeof(float));
    d->decWindow  = (float*)calloc((size_t)d->decWindowSize, sizeof(float));
    d->diff       = (float*)calloc((size_t)(d->decWindowSize / 2 + 1), sizeof(float));
    d->cmndf      = (float*)calloc((size_t)(d->decWindowSize / 2 + 1), sizeof(float));
    d->decRing    = (float*)calloc((size_t)d->decWindowSize, sizeof(float));

    if (!d->ringBuffer || !d->window || !d->decWindow || !d->diff || !d->cmndf || !d->decRing) {
        pitchDetectorDestroy(d);
        return NULL;
    }

    // Build FIR coefficients
    build_fir_lowpass(d->firCoef, PITCH_FIR_TAPS, 0.085f);

    d->result.hz         = -1.0f;
    d->result.confidence = 0.0f;
    d->result.sequence   = 0;

    return d;
}

void pitchDetectorDestroy(PitchDetector* d) {
    if (!d) return;
    free(d->ringBuffer);
    free(d->window);
    free(d->decWindow);
    free(d->diff);
    free(d->cmndf);
    free(d->decRing);
    free(d);
}

void pitchDetectorConfigure(PitchDetector* d, PitchDetectorConfig config) {
    if (!d) return;
    int divisor = (config.hopDivisor > 0) ? config.hopDivisor : 8;
    d->config      = config;
    d->hopSize     = d->windowSize / divisor;
    d->maxTauD     = compute_max_tau(d->decRate, config.minHz, d->decWindowSize / 2);
    d->minTauD     = compute_min_tau(d->decRate, config.maxHz);
    d->maxTauPlusD = (d->maxTauD + 1 < d->decWindowSize / 2) ? d->maxTauD + 1 : d->decWindowSize / 2;
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

// Queue a result (and mirror it as the "latest" for getResult).
static void push_result(PitchDetector* d, float hz, float confidence) {
    d->result.hz         = hz;
    d->result.confidence = confidence;
    d->result.sequence++;

    if (d->fifoCount == PITCH_RESULT_FIFO) {
        d->fifoRead = (d->fifoRead + 1) % PITCH_RESULT_FIFO;   // drop oldest
        d->fifoCount--;
    }
    d->fifo[d->fifoWrite] = d->result;
    d->fifoWrite = (d->fifoWrite + 1) % PITCH_RESULT_FIFO;
    d->fifoCount++;
}

// Write invalid result and increment sequence.
static void invalidate_result(PitchDetector* d) {
    push_result(d, -1.0f, 0.0f);
}

// Measure RMS and peak of window.
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

// Direct YIN squared difference function, on the decimated (coarse-stage) window.
// d(0) = 0; d(τ) = Σ_{i=0}^{N-τ-1} (x[i] − x[i+τ])²
// Direct (non-FFT) computation: linear autocorrelation, no circular artifacts.
// Computed only up to maxTauD + 1 (parabolic neighbour for refinement).
static void compute_yin_difference(PitchDetector* d) {
    const int N = d->decWindowSize;
    const int maxTauPlus = d->maxTauPlusD;

    d->diff[0] = 0.0f;

    for (int tau = 1; tau <= maxTauPlus; tau++) {
        float sum = 0.0f;
        const int limit = N - tau;
        for (int i = 0; i < limit; i++) {
            float delta = d->decWindow[i] - d->decWindow[i + tau];
            sum += delta * delta;
        }
        d->diff[tau] = sum;
    }
}

// CMNDF from diff[] (decimated domain).
// d'(0) = 1; d'(τ) = d(τ)·τ / Σ_{j=1}^{τ} d(j)
// Computed only up to maxTauD + 1 (parabolic neighbour for refinement).
static void compute_cmndf(PitchDetector* d) {
    const int maxTauPlus = d->maxTauPlusD;

    float cumSum  = 0.0f;
    d->cmndf[0] = 1.0f;

    for (int tau = 1; tau <= maxTauPlus; tau++) {
        cumSum += d->diff[tau];
        d->cmndf[tau] = (cumSum > 0.0f)
            ? d->diff[tau] * (float)tau / cumSum
            : 1.0f;
    }
}

// Find the first CMNDF pit whose minimum is below threshold.
// Searches only from minTau to maxTau to avoid the systematic
// downward CMNDF bias that appears near halfWindow in the direct YIN method.
// Returns -1 if no qualifying pit found.
static int find_yin_period(const float* cmndf, int minTau, int maxTau, float threshold) {
    for (int tau = minTau; tau <= maxTau; tau++) {
        if (cmndf[tau] < threshold) {
            while (tau + 1 <= maxTau && cmndf[tau + 1] < cmndf[tau]) tau++;
            return tau;
        }
    }
    return -1;
}

// Fallback: return the tau of the global CMNDF minimum within minTau..maxTau if
// it is below fallbackThreshold. Used when no pit satisfies yinThreshold.
// Capping at maxTau prevents spurious sub-bass detections caused by the
// direct-YIN CMNDF bias: at large τ (near halfWindow=N/2) the expected
// CMNDF for noise is (N-τ)/(N-τ/2) ≈ 0.67 rather than 1.0, so the region
// near halfWindow always looks like a broad, low-valued "dip".
static int find_best_period_fallback(const float* cmndf, int minTau, int maxTau, float fallbackThreshold) {
    int   bestTau = -1;
    float bestVal = 1.0f;
    for (int tau = minTau; tau <= maxTau; tau++) {
        if (cmndf[tau] < bestVal) {
            bestVal = cmndf[tau];
            bestTau = tau;
        }
    }
    return (bestTau >= 0 && bestVal < fallbackThreshold) ? bestTau : -1;
}

// Harmonic correction: if YIN sub-harmonic-locks (finds a longer period /
// lower frequency than the true fundamental), check period/2..period/5.
// If a shorter period's CMNDF is within octaveTolerance of the detected
// period's CMNDF, the shorter period (higher frequency, truer fundamental)
// is preferred. Returns on the first qualifying sub-period (N=2 first),
// so the smallest frequency jump is always tried before larger ones.
//
// Guard: only fire if period > maxTau/3. Notes already in the upper ⅔ of
// the valid frequency range are assumed correctly detected — applying tau/N
// there would risk jumping to a harmonic (e.g., G3→G4 via tau=225→112).
//
// Example: playing A2 (110 Hz, tau=437) but YIN finds tau=875 (A1, 55 Hz)
// because CMNDF(875)=0.244 < CMNDF(437)=0.271. With tolerance=0.06:
// threshold=0.244+0.06=0.304 > 0.271 → correct back to tau=437 (A2). ✓
// For a clean A1 signal: CMNDF(875)=0.08, CMNDF(437)=0.16,
// threshold=0.08+0.06=0.14 < 0.16 → no correction. ✓
// For G3 (tau=225) with maxTau=1633: 225 ≤ 1633/3=544 → guard skips. ✓
//
// NOTE: since the two-stage rework, cmndf/period/minTau/maxTau here are all in
// the decimated (coarse-stage) domain — the illustrative tau values above are
// from the full sample rate and are ~PITCH_DECIM_FACTOR× larger than the
// decimated tau this function actually sees. The period/N ratio logic and
// tolerance comparisons are unaffected by the rate change.
static int correct_harmonic_period(const float* cmndf, int period, int minTau, int maxTau, float tolerance) {
    if (tolerance <= 0.0f) return period;
    // Only correct when the note is in the lower third of the valid range.
    if (period <= maxTau / 3) return period;
    float threshold = cmndf[period] + tolerance;
    for (int N = 2; N <= 5; N++) {
        int sub = period / N;
        if (sub < minTau) break;
        if (sub <= maxTau && cmndf[sub] < threshold) {
            return sub;
        }
    }
    return period;
}

// Sub-fundamental correction: if YIN over-locked on a harmonic (finds a shorter
// period / higher frequency than the true fundamental), check period*2.
// If the doubled period has a CMNDF pit both below yinThreshold AND deeper than
// the current pit, the doubled period (lower frequency, true fundamental) is
// preferred.
//
// Physics: for a wound string whose 2nd harmonic dominates the spectrum, the
// CMNDF has a pit at T/2 (the harmonic period) that YIN finds first. But the
// true fundamental period T also has a pit — and typically a deeper one — because
// all harmonics align at T. The condition cmndf[2τ] < cmndf[τ] captures this.
//
// This correction is complementary to correct_harmonic_period (which corrects
// sub-harmonic locks: too-low detections). These two are never run together —
// only the sub-fundamental runs when the octave-up correction did not fire.
//
// Guard: period*2 must stay within maxTau (valid detection range).
static int correct_sub_harmonic_period(const float* cmndf, int period, int maxTau, float threshold) {
    int doubled = period * 2;
    if (doubled > maxTau) return period;
    // Mirror of the harmonic-correction guard: only fire for the upper third of the
    // valid period range (i.e., short periods = high detected frequencies). A long
    // period already represents a low note — doubling it pushes into sub-bass and is
    // almost certainly wrong (CMNDF bias at large τ can make cmndf[2τ] < cmndf[τ]
    // spuriously). The upper-third guard matches the harmonic correction's lower-third
    // guard so the two corrections are strictly complementary.
    // Pure-sine guard: if the current pit is very deep, all multiples also dip near
    // zero. The frequency-range guard subsumes this case for bass, but keep the deep-
    // pit guard as a second line of defence for mid-range pure signals.
    if (period > maxTau / 3) return period;
    if (cmndf[period] < 0.05f) return period;
    // Require the sub-fundamental pit to be meaningfully deeper (not just marginally),
    // so that a genuine fundamental with a slightly noisy CMNDF is not mis-doubled.
    // 0.02 margin: the difference between A2 cmndf≈0.055 and A1 cmndf≈0.060 is ~0.005,
    // below the margin → no spurious doubling. Wound-string cases (G3: 0.12 vs 0.04)
    // have a difference of ~0.08, well above the margin → correction still fires.
    if (cmndf[doubled] < threshold && cmndf[doubled] < cmndf[period] - 0.02f) return doubled;
    return period;
}

// Local CMNDF contrast: average value in a ±12-sample window (excluding ±2
// around the pit centre) minus the pit value. Positive = clear pit, negative
// = flat/noisy region.
static float local_contrast(const float* cmndf, int tau, int H) {
    int left  = tau - 12; if (left  < 1) left  = 1;
    int right = tau + 12; if (right > H) right = H;
    float sum = 0.0f;
    int   count = 0;
    for (int i = left; i <= right; i++) {
        if (abs(i - tau) <= 2) continue;
        sum += cmndf[i];
        count++;
    }
    if (count <= 0) return 0.0f;
    return sum / (float)count - cmndf[tau];
}

static float clamp01(float x) {
    if (x < 0.0f) return 0.0f;
    if (x > 1.0f) return 1.0f;
    return x;
}

// Confidence from pit depth and local contrast.
// depthScore   = 1 − CMNDF(period)  — how deep the pit is
// contrastScore = contrast × 8      — how distinct the pit is from neighbours
// confidence  = 0.75·depth + 0.25·contrast
static float compute_confidence(float cmndfValue, float contrast) {
    float depthScore    = clamp01(1.0f - cmndfValue);
    float contrastScore = clamp01(contrast * 8.0f);
    return clamp01(depthScore * 0.75f + contrastScore * 0.25f);
}

// ==========================================================================
//  Full-rate refinement
// ==========================================================================

// Full-rate refinement. The coarse stage found the period at the decimated
// rate; the raw difference function d(τ) at the original rate has its
// corresponding minimum near τ0 = period × PITCH_DECIM_FACTOR. We search a
// ±(FACTOR+2) band and parabolically interpolate the raw d(τ) minimum — the
// CMNDF normalization cannot be computed locally (it needs the cumulative sum
// from τ=1) and is unnecessary for sub-sample placement of an isolated pit.
static float refine_full_rate(PitchDetector* d, int decPeriod) {
    const int N    = d->windowSize;
    const int band = PITCH_DECIM_FACTOR + 2;
    int tau0 = decPeriod * PITCH_DECIM_FACTOR;
    int lo = tau0 - band, hi = tau0 + band;
    if (lo < 2) lo = 2;
    if (hi > N / 2) hi = N / 2;

    float dvals[2 * (PITCH_DECIM_FACTOR + 2) + 1];
    int   count = hi - lo + 1;
    int   bestIdx = 0;
    for (int t = lo; t <= hi; t++) {
        float sum = 0.0f;
        const float* x = d->window;
        const int limit = N - t;
        for (int i = 0; i < limit; i++) {
            float delta = x[i] - x[i + t];
            sum += delta * delta;
        }
        dvals[t - lo] = sum;
        if (sum < dvals[bestIdx]) bestIdx = t - lo;
    }

    int best = lo + bestIdx;
    float refined = (float)best;
    if (bestIdx > 0 && bestIdx < count - 1) {
        float alpha = dvals[bestIdx - 1], beta = dvals[bestIdx], gamma = dvals[bestIdx + 1];
        float denom = 2.0f * (2.0f * beta - alpha - gamma);
        if (fabsf(denom) > 1e-12f) {
            float offset = (gamma - alpha) / denom;
            if (offset >  0.5f) offset =  0.5f;
            if (offset < -0.5f) offset = -0.5f;
            refined += offset;
        }
    }
    return d->sampleRate / refined;
}

// ==========================================================================
//  Analysis (called internally every hopSize samples)
// ==========================================================================

static void run_analysis(PitchDetector* d) {
    const int N   = d->windowSize;
    const int ND  = d->decWindowSize;
    const int MT  = d->maxTauD;
    const int mnT = d->minTauD;

    // Step 1 — gates measured on the full-rate window (input-level semantics)
    float rms, peak;
    measure_signal(d->window, N, &rms, &peak);
    if (rms < d->config.powerThreshold || peak < d->config.peakThreshold) {
        invalidate_result(d);
        return;
    }

    // Step 2 — DC removal on the decimated window (YIN input)
    remove_dc(d->decWindow, ND);

    // Step 3 — coarse YIN on the decimated window
    compute_yin_difference(d);   // operates on decWindow/ND, tau <= maxTauPlusD
    compute_cmndf(d);            // tau <= maxTauPlusD

    int period = find_yin_period(d->cmndf, mnT, MT, d->config.yinThreshold);
    int usedFallback = 0;
    if (period < 0) {
        period = find_best_period_fallback(d->cmndf, mnT, MT, d->config.fallbackThreshold);
        usedFallback = (period >= 0);
    }
    if (period < 0) {
        invalidate_result(d);
#if DEBUG_PITCH
        {
            int   gt = -1; float gv = 1.0f;
            for (int t = mnT; t <= MT; t++) {
                if (d->cmndf[t] < gv) { gv = d->cmndf[t]; gt = t; }
            }
            printf("[DSP] rms=%.5f peak=%.5f | INVALID globalMin tau=%d cmndf=%.4f hz=%.1f\n",
                   rms, peak, gt, gv, gt > 0 ? d->decRate / (float)gt : -1.0f);
        }
#endif
        return;
    }

    // Step 4 — harmonic corrections (decimated domain, unchanged logic)
    int rawPeriod = period;
    period = correct_harmonic_period(d->cmndf, period, mnT, MT, d->config.octaveTolerance);
    if (period == rawPeriod) {
        period = correct_sub_harmonic_period(d->cmndf, period, MT, d->config.yinThreshold);
    }

#if DEBUG_PITCH
    if (period != rawPeriod) {
        printf("[DSP] harmonic-correct: %d(%.1fHz)→%d(%.1fHz) cmndf %.3f→%.3f\n",
               rawPeriod, d->decRate / (float)rawPeriod,
               period,    d->decRate / (float)period,
               d->cmndf[rawPeriod], d->cmndf[period]);
    }
#endif

    // Step 5 — confidence gate from the coarse CMNDF
    float contrast   = local_contrast(d->cmndf, period, MT);
    float confidence = compute_confidence(d->cmndf[period], contrast);
    if (confidence < d->config.minConfidence) {
        invalidate_result(d);
#if DEBUG_PITCH
        printf("[DSP] rms=%.5f peak=%.5f | lowConf rawTau=%d corrTau=%d cmndf=%.3f contrast=%.3f conf=%.3f\n",
               rms, peak, rawPeriod, period, d->cmndf[period], contrast, confidence);
#endif
        return;
    }

    // Step 6 — full-rate refine: raw difference function in a narrow band
    // around 4×period, local minimum + parabolic interpolation.
    float hz = refine_full_rate(d, period);

    push_result(d, hz, confidence);

#if DEBUG_PITCH
    printf("[DSP] rawTau=%d corrTau=%d hz=%.2f cmndf=%.3f contrast=%.3f conf=%.3f rms=%.5f peak=%.5f fallback=%d\n",
           rawPeriod, period, hz, d->cmndf[period], contrast, confidence, rms, peak, usedFallback);
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

        // Feed FIR and decimate
        d->firHist[d->firHistPos] = samples[i];
        d->firHistPos = (d->firHistPos + 1) % PITCH_FIR_TAPS;
        if (++d->decimPhase == PITCH_DECIM_FACTOR) {
            d->decimPhase = 0;
            float acc = 0.0f;
            int idx = d->firHistPos;   // oldest sample position
            for (int k = 0; k < PITCH_FIR_TAPS; k++) {
                acc += d->firCoef[k] * d->firHist[(idx + k) % PITCH_FIR_TAPS];
            }
            d->decRing[d->decWritePos] = acc;
            d->decWritePos = (d->decWritePos + 1) % d->decWindowSize;
        }

        if (d->samplesSinceAnalysis >= HS) {
            d->samplesSinceAnalysis = 0;

            // Copy most recent N samples from ring buffer (chronological order)
            int start = d->writePos;
            for (int j = 0; j < N; j++)
                d->window[j] = d->ringBuffer[(start + j) % N];
            int dstart = d->decWritePos;
            for (int j = 0; j < d->decWindowSize; j++)
                d->decWindow[j] = d->decRing[(dstart + j) % d->decWindowSize];

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

// ==========================================================================
//  DrainResults (single-thread with process; see header)
// ==========================================================================

int pitchDetectorDrainResults(PitchDetector* d, PitchResult* out, int maxCount) {
    if (!d || !out || maxCount <= 0) return 0;
    int n = 0;
    while (n < maxCount && d->fifoCount > 0) {
        out[n++] = d->fifo[d->fifoRead];
        d->fifoRead = (d->fifoRead + 1) % PITCH_RESULT_FIFO;
        d->fifoCount--;
    }
    return n;
}
