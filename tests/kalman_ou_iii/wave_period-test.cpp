// Copyright (c) 2025-2026  Mikhail Grushinskiy
//
// The wave-band period estimator follows the physical zero-crossing period
// through spectral moments.  The canonical output is a log-domain time-scale
// state: period and frequency must therefore remain exact reciprocals, and the
// moment horizon must follow the smoothed positive state rather than the raw
// nonlinear ratio.

#define EIGEN_NON_ARDUINO

#include "tuner/WavePeriodEstimator.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

namespace {

constexpr float kDt = 1.0f / 200.0f;
constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;

bool check(bool condition, const char* message) {
    if (!condition) std::fprintf(stderr, "%s\n", message);
    return condition;
}

bool reciprocal_contract(const WavePeriodEstimator& estimator) {
    const float T = estimator.getPeriodSec();
    const float f = estimator.getFrequencyHz();
    return std::isfinite(T) && std::isfinite(f) &&
           std::abs(T * f - 1.0f) < 2e-6f;
}

// Single sinusoid: acceleration -omega^2 A sin(omega t) belongs to elevation
// A sin(omega t), so T_z is exactly the sinusoid period.
bool test_monochromatic() {
    for (float period : {2.5f, 5.0f, 8.5f, 14.0f}) {
        WavePeriodEstimator estimator;
        const float omega = kTwoPi / period;
        const float amplitude = 1.5f;
        const int steps = static_cast<int>(600.0f / kDt);
        for (int i = 0; i < steps; ++i) {
            const float t = static_cast<float>(i) * kDt;
            estimator.update(kDt, -amplitude * omega * omega * std::sin(omega * t));
        }
        if (!check(estimator.isReady(), "monochromatic estimate never became ready")) return false;
        if (!check(reciprocal_contract(estimator), "period/frequency are not exact reciprocal outputs")) return false;
        const float relative = std::abs(estimator.getPeriodSec() - period) / period;
        if (!check(relative < 0.02f, "monochromatic period is off by more than 2 percent")) {
            std::fprintf(stderr, "  period %.2f estimated %.3f\n", period, estimator.getPeriodSec());
            return false;
        }
    }
    return true;
}

// Broadband sea. T_z = 2*pi*sqrt(m0/m2) is computed from the same components
// that build the signal, so the reference needs no spectral estimation.
bool test_broadband(unsigned seed, float peak_period, float* estimated, float* reference) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> phase(0.0f, kTwoPi);

    const float peak = kTwoPi / peak_period;
    std::vector<float> omega;
    std::vector<float> amplitude;
    std::vector<float> offset;

    const int components = 60;
    double m0 = 0.0;
    double m2 = 0.0;
    for (int i = 0; i < components; ++i) {
        const float ratio = 0.5f + 2.5f * static_cast<float>(i) / static_cast<float>(components - 1);
        const float w = peak * ratio;
        const float density = (ratio < 1.0f) ? std::pow(ratio, 4.0f) : std::pow(ratio, -5.0f);
        const float a = std::sqrt(density);
        omega.push_back(w);
        amplitude.push_back(a);
        offset.push_back(phase(rng));
        m0 += 0.5 * double(a) * double(a);
        m2 += 0.5 * double(a) * double(a) * double(w) * double(w);
    }
    *reference = float(kTwoPi / std::sqrt(m2 / m0));

    WavePeriodEstimator estimator;
    const int steps = static_cast<int>(1200.0f / kDt);
    for (int i = 0; i < steps; ++i) {
        const float t = static_cast<float>(i) * kDt;
        float acceleration = 0.0f;
        for (size_t k = 0; k < omega.size(); ++k) {
            acceleration -= amplitude[k] * omega[k] * omega[k] * std::sin(omega[k] * t + offset[k]);
        }
        estimator.update(kDt, acceleration);
    }
    *estimated = estimator.getPeriodSec();
    return estimator.isReady() && reciprocal_contract(estimator);
}

bool test_broadband_tracks_sea_state() {
    float previous = 0.0f;
    for (float peak_period : {3.0f, 5.7f, 8.5f, 11.4f}) {
        float estimated = 0.0f;
        float reference = 0.0f;
        if (!check(test_broadband(20260803u, peak_period, &estimated, &reference),
                   "broadband estimate never became ready or broke reciprocal contract")) return false;
        const float relative = std::abs(estimated - reference) / reference;
        if (!check(relative < 0.06f, "broadband T_z is off by more than 6 percent")) {
            std::fprintf(stderr, "  peak %.1f s: estimated %.3f reference %.3f\n",
                         peak_period, estimated, reference);
            return false;
        }
        if (!check(estimated > previous, "T_z did not grow with the sea state")) return false;
        previous = estimated;
    }
    return true;
}

// With a deliberately slow log smoother, a period step makes the raw moment
// ratio move first.  The next moment horizon must remain tied to the canonical
// period, not jump to the raw candidate.  This is the architectural separation
// that removes the old direct self-modulation of estimator bandwidth.
bool test_moment_horizon_uses_canonical_period() {
    WavePeriodEstimator estimator(0.02f, 4.0f, 2.0f, 1.0f, 180.0f);

    auto drive = [&](float period, float seconds, float t0) {
        const float omega = kTwoPi / period;
        const int steps = static_cast<int>(seconds / kDt);
        for (int i = 0; i < steps; ++i) {
            const float t = t0 + static_cast<float>(i) * kDt;
            estimator.update(kDt, -omega * omega * std::sin(omega * t));
        }
    };

    drive(4.0f, 180.0f, 0.0f);
    if (!check(estimator.isReady(), "canonical-horizon test never became ready")) return false;
    drive(10.0f, 12.0f, 180.0f);

    const float raw = estimator.getRawPeriodSec();
    const float canonical = estimator.getPeriodSec();
    const float horizon = estimator.getMomentHorizonSec();
    if (!check(std::isfinite(raw) && std::isfinite(canonical) && raw > canonical,
               "period step did not separate raw and canonical estimates")) return false;
    const float canonical_horizon = 4.0f * canonical;
    const float raw_horizon = 4.0f * raw;
    return check(std::abs(horizon - canonical_horizon) < std::abs(horizon - raw_horizon),
                 "moment horizon follows raw period instead of canonical log-period state");
}

// A constant offset is what an accelerometer bias looks like; the leak has to
// reject it rather than integrate it into a ramp.
bool test_offset_rejection() {
    WavePeriodEstimator estimator;
    const float period = 8.0f;
    const float omega = kTwoPi / period;
    const int steps = static_cast<int>(900.0f / kDt);
    for (int i = 0; i < steps; ++i) {
        const float t = static_cast<float>(i) * kDt;
        estimator.update(kDt, -omega * omega * std::sin(omega * t) + 0.35f);
    }
    if (!check(estimator.isReady(), "offset case never became ready")) return false;
    const float relative = std::abs(estimator.getPeriodSec() - period) / period;
    if (!check(relative < 0.05f, "constant acceleration offset moved the period estimate")) {
        std::fprintf(stderr, "  estimated %.3f\n", estimator.getPeriodSec());
        return false;
    }
    return true;
}

bool test_rate_invariance() {
    const float period = 6.5f;
    const float omega = kTwoPi / period;
    float first = 0.0f;
    for (float dt : {1.0f / 100.0f, 1.0f / 200.0f, 1.0f / 400.0f}) {
        WavePeriodEstimator estimator;
        const int steps = static_cast<int>(600.0f / dt);
        for (int i = 0; i < steps; ++i) {
            const float t = static_cast<float>(i) * dt;
            estimator.update(dt, -omega * omega * std::sin(omega * t));
        }
        if (!check(estimator.isReady(), "rate case never became ready")) return false;
        if (!check(reciprocal_contract(estimator), "rate case broke reciprocal period/frequency contract")) return false;
        if (first == 0.0f) {
            first = estimator.getPeriodSec();
        } else if (!check(std::abs(estimator.getPeriodSec() - first) / first < 0.01f,
                          "period estimate depends on the sample rate")) {
            return false;
        }
    }
    return true;
}

bool test_survives_instrument_disturbances() {
    std::mt19937 rng(12345);
    std::normal_distribution<float> white(0.0f, 0.35f);
    std::normal_distribution<float> walk(0.0f, 1e-3f * std::sqrt(kDt));

    float previous = 0.0f;
    for (float period : {5.0f, 8.5f, 12.0f}) {
        WavePeriodEstimator estimator;
        const float omega = kTwoPi / period;
        const float amplitude = 1.2f;
        float bias = 0.0f;
        const int steps = static_cast<int>(1200.0f / kDt);
        for (int i = 0; i < steps; ++i) {
            const float t = static_cast<float>(i) * kDt;
            bias += walk(rng);
            estimator.update(kDt,
                             -amplitude * omega * omega * std::sin(omega * t)
                                 + bias + white(rng));
        }
        if (!check(estimator.isReady(), "disturbed estimate never became ready")) return false;
        const float relative = std::abs(estimator.getPeriodSec() - period) / period;
        if (!check(relative < 0.15f, "instrument disturbances moved the period by more than 15 percent")) {
            std::fprintf(stderr, "  period %.1f estimated %.3f\n", period, estimator.getPeriodSec());
            return false;
        }
        if (!check(estimator.getPeriodSec() > previous, "disturbed estimate lost the ordering")) return false;
        previous = estimator.getPeriodSec();
    }
    return true;
}

bool test_rejects_bad_samples_and_resets() {
    WavePeriodEstimator estimator;
    const float omega = kTwoPi / 7.0f;
    for (int i = 0; i < static_cast<int>(900.0f / kDt); ++i) {
        const float t = static_cast<float>(i) * kDt;
        estimator.update(kDt, -omega * omega * std::sin(omega * t));
    }
    if (!check(estimator.isReady(), "clean run never became ready")) return false;
    const float before = estimator.getPeriodSec();

    estimator.update(-kDt, 1.0f);
    estimator.update(0.0f, 1.0f);
    estimator.update(kDt, std::numeric_limits<float>::quiet_NaN());
    estimator.update(std::numeric_limits<float>::infinity(), 1.0f);
    if (!check(estimator.getPeriodSec() == before, "a malformed sample changed the estimate")) return false;

    estimator.reset();
    if (!check(!estimator.isReady() && !std::isfinite(estimator.getPeriodSec()),
               "reset left the estimator ready")) return false;
    return true;
}

bool test_no_signal_stays_unready() {
    WavePeriodEstimator estimator;
    for (int i = 0; i < static_cast<int>(300.0f / kDt); ++i) {
        estimator.update(kDt, 0.0f);
    }
    return check(!estimator.isReady() || !std::isfinite(estimator.getPeriodSec()),
                 "a silent input produced a confident period");
}

} // namespace

int main() {
    if (!test_monochromatic()) return 1;
    if (!test_broadband_tracks_sea_state()) return 1;
    if (!test_moment_horizon_uses_canonical_period()) return 1;
    if (!test_offset_rejection()) return 1;
    if (!test_rate_invariance()) return 1;
    if (!test_survives_instrument_disturbances()) return 1;
    if (!test_rejects_bad_samples_and_resets()) return 1;
    if (!test_no_signal_stays_unready()) return 1;
    std::printf("wave_period-test: OK\n");
    return 0;
}
