// test_internal.c — white-box tests for PitchDSP internals.
// Includes the implementation directly to reach static functions.

#include "../Sources/PitchDSP/PitchDSP.c"

#include <stdio.h>
#include <stdlib.h>

static int t_passed = 0, t_failed = 0;
#define CHECK(cond, fmt, ...) do {                                        \
    if (cond) { t_passed++; }                                             \
    else { t_failed++; printf("  FAIL: " fmt "\n", ##__VA_ARGS__); }      \
} while (0)

static float rms(const float* x, int n) {
    double s = 0;
    for (int i = 0; i < n; i++) s += (double)x[i] * x[i];
    return (float)sqrt(s / n);
}

// Run a buffer through a detector's FIR+decimate path and collect the
// decimated ring content into out (must hold windowSize/PITCH_DECIM_FACTOR).
static int collect_decimated(PitchDetector* d, const float* in, int n, float* out, int outCap) {
    pitchDetectorProcess(d, in, n);
    int dn = d->decWindowSize;
    if (dn > outCap) dn = outCap;
    int start = d->decWritePos;
    for (int j = 0; j < dn; j++) out[j] = d->decRing[(start + j) % d->decWindowSize];
    return dn;
}

static void test_fir_unity_dc_gain(void) {
    float h[PITCH_FIR_TAPS];
    build_fir_lowpass(h, PITCH_FIR_TAPS, 0.085f);
    float sum = 0;
    for (int i = 0; i < PITCH_FIR_TAPS; i++) sum += h[i];
    CHECK(fabsf(sum - 1.0f) < 1e-4f, "FIR DC gain %.5f != 1.0", sum);
}

static void test_decimator_passband(void) {
    // 440 Hz at 48 kHz is deep in the passband: decimated RMS ≈ input RMS.
    PitchDetectorConfig cfg = pitchDetectorDefaultConfig();
    PitchDetector* d = pitchDetectorCreate(8192, 48000.0f, cfg);
    int n = 48000;
    float* s = malloc(n * sizeof(float));
    for (int i = 0; i < n; i++) s[i] = sinf(2.0f * (float)M_PI * 440.0f * i / 48000.0f);
    float dec[4096];
    int dn = collect_decimated(d, s, n, dec, 4096);
    float r = rms(dec, dn);
    CHECK(fabsf(r - 0.7071f) < 0.05f, "passband RMS %.3f, want ~0.707", r);
    free(s); pitchDetectorDestroy(d);
}

static void test_decimator_alias_rejection(void) {
    // 8 kHz at 48 kHz is above the stopband edge; after decimation to 12 kHz
    // it would alias to 4 kHz — the FIR must crush it first.
    PitchDetectorConfig cfg = pitchDetectorDefaultConfig();
    PitchDetector* d = pitchDetectorCreate(8192, 48000.0f, cfg);
    int n = 48000;
    float* s = malloc(n * sizeof(float));
    for (int i = 0; i < n; i++) s[i] = sinf(2.0f * (float)M_PI * 8000.0f * i / 48000.0f);
    float dec[4096];
    int dn = collect_decimated(d, s, n, dec, 4096);
    float r = rms(dec, dn);
    CHECK(r < 0.03f, "alias RMS %.4f, want < 0.03 (>-27 dB rejection)", r);
    free(s); pitchDetectorDestroy(d);
}

int main(void) {
    printf("=== PitchDSP internal tests ===\n");
    test_fir_unity_dc_gain();
    test_decimator_passband();
    test_decimator_alias_rejection();
    printf("=== internal: %d passed, %d failed ===\n", t_passed, t_failed);
    return t_failed > 0 ? 1 : 0;
}
