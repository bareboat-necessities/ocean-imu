#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

#include "tuner/AdaptiveWaveBandPass.h"

namespace {

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

constexpr float FS = 200.0f;
constexpr float DT = 1.0f / FS;
constexpr float TWO_PI = 6.28318530717958647692f;

float sine_rms(float f_ref, float f_signal) {
    AdaptiveWaveBandPass band;
    constexpr int N = int(FS * 120.0f);
    constexpr int START = int(FS * 60.0f);
    double e2 = 0.0;
    int count = 0;
    for (int k = 0; k < N; ++k) {
        const float t = static_cast<float>(k) * DT;
        const float x = std::sin(TWO_PI * f_signal * t);
        const float y = band.step(x, DT, f_ref);
        if (k >= START) {
            e2 += static_cast<double>(y) * static_cast<double>(y);
            ++count;
        }
    }
    return static_cast<float>(std::sqrt(e2 / static_cast<double>(count)));
}

bool test_normalized_response_scales_with_reference_frequency() {
    bool ok = true;

    const float r1 = sine_rms(0.10f, 0.10f);
    const float r2 = sine_rms(0.20f, 0.20f);
    const float r3 = sine_rms(0.40f, 0.40f);

    const float max_r = std::max(r1, std::max(r2, r3));
    const float min_r = std::min(r1, std::min(r2, r3));
    ok &= check((max_r - min_r) / max_r < 0.01f,
                "equal normalized frequencies must have nearly equal gain");

    AdaptiveWaveBandPass band;
    band.step(0.0f, DT, 0.20f);
    ok &= check(std::fabs(band.lowHz() - 0.10f) < 1e-6f,
                "default lower corner must be 0.5 f_ref");
    ok &= check(std::fabs(band.highHz() - 0.80f) < 1e-6f,
                "default upper corner must be 4 f_ref");

    std::cout << "  normalized pass-band RMS: " << r1 << ", " << r2
              << ", " << r3 << '\n';
    return ok;
}

bool test_out_of_band_rejection() {
    bool ok = true;
    constexpr float F_REF = 0.20f;

    const float pass = sine_rms(F_REF, F_REF);
    const float low = sine_rms(F_REF, 0.02f);
    const float high = sine_rms(F_REF, 2.0f);

    ok &= check(low < 0.30f * pass,
                "sub-band acceleration must be strongly attenuated");
    ok &= check(high < 0.50f * pass,
                "high-frequency acceleration must be attenuated");

    std::cout << "  selectivity RMS low/pass/high: " << low << " / "
              << pass << " / " << high << '\n';
    return ok;
}

// Deterministic Rademacher input: exactly zero mean and unit variance in the
// population, without depending on a platform random-number implementation.
uint32_t xorshift32(uint32_t& state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

bool test_white_noise_variance_gain_matches_observed_output() {
    bool ok = true;
    AdaptiveWaveBandPass band;
    uint32_t rng = 0x12345678u;

    constexpr int N = 400000;
    constexpr int START = 100000;
    double sum = 0.0;
    double sum2 = 0.0;
    int count = 0;

    for (int k = 0; k < N; ++k) {
        const float x = (xorshift32(rng) & 1u) ? 1.0f : -1.0f;
        const float y = band.step(x, DT, 0.20f);
        if (k >= START) {
            sum += y;
            sum2 += static_cast<double>(y) * static_cast<double>(y);
            ++count;
        }
    }

    const double mean = sum / static_cast<double>(count);
    const double observed = sum2 / static_cast<double>(count) - mean * mean;
    const double predicted = band.whiteNoiseVarianceGain();
    const double rel = std::fabs(observed - predicted) / predicted;

    ok &= check(predicted > 0.0 && predicted < 1.0,
                "white-noise variance gain must be finite and attenuating");
    ok &= check(rel < 0.04,
                "propagated white-noise variance must match observed output");

    std::cout << "  white-noise gain predicted/observed: " << predicted
              << " / " << observed << "\n";
    return ok;
}

bool test_moving_band_stays_finite_and_continuous() {
    bool ok = true;
    AdaptiveWaveBandPass band;
    float previous = 0.0f;

    for (int k = 0; k < int(FS * 40.0f); ++k) {
        const float t = static_cast<float>(k) * DT;
        const float f_ref = 0.08f + 0.24f * (t / 40.0f);
        const float x = std::sin(TWO_PI * f_ref * t);
        const float y = band.step(x, DT, f_ref);

        ok &= check(std::isfinite(y), "moving-band output must stay finite");
        ok &= check(std::isfinite(band.whiteNoiseVarianceGain()),
                    "moving-band noise gain must stay finite");
        ok &= check(band.highHz() > band.lowHz(),
                    "moving-band upper corner must stay above lower corner");

        if (k > 0) {
            ok &= check(std::fabs(y - previous) < 0.2f,
                        "smoothly moving corners must not create output jumps");
        }
        previous = y;
        if (!ok) break;
    }
    return ok;
}

}  // namespace

int main() {
    bool ok = true;
    ok &= test_normalized_response_scales_with_reference_frequency();
    ok &= test_out_of_band_rejection();
    ok &= test_white_noise_variance_gain_matches_observed_output();
    ok &= test_moving_band_stays_finite_and_continuous();

    if (!ok) {
        std::cerr << "wave-band sigma checks FAILED\n";
        return 1;
    }
    std::cout << "all wave-band sigma checks passed\n";
    return 0;
}
