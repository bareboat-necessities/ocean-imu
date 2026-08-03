#include "freq/WaveBandPeriodEstimator.h"

#include <cmath>
#include <iostream>
#include <numbers>

namespace {

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

bool test_single_tone() {
    constexpr double dt = 0.005;
    for (double period : {3.0, 6.0, 10.0}) {
        WaveBandPeriodEstimator<double> estimator(0.02, 60.0, 120.0);
        const double omega = 2.0 * std::numbers::pi / period;
        for (int i = 0; i < static_cast<int>(900.0 / dt); ++i) {
            const double t = i * dt;
            const double acceleration = 0.17 - 1.4 * omega * omega * std::sin(omega * t);
            estimator.update(acceleration, dt);
        }
        const double expected = 2.0 * std::numbers::pi /
            std::sqrt(omega * omega + std::pow(2.0 * std::numbers::pi * 0.02, 2));
        const double relative_error = std::abs(estimator.getPeriodSec() - expected) / expected;
        if (!check(estimator.isReady(), "single-tone estimator did not become ready")) return false;
        if (!check(relative_error < 0.015, "single-tone period mismatch")) return false;
    }
    return true;
}

bool test_multitone_zero_crossing_period() {
    constexpr double dt = 0.005;
    constexpr double f1 = 0.11;
    constexpr double f2 = 0.24;
    constexpr double a1 = 1.0;
    constexpr double a2 = 0.35;
    const double w1 = 2.0 * std::numbers::pi * f1;
    const double w2 = 2.0 * std::numbers::pi * f2;
    const double lambda = 2.0 * std::numbers::pi * 0.02;

    WaveBandPeriodEstimator<double> estimator(0.02, 90.0, 180.0);
    for (int i = 0; i < static_cast<int>(1200.0 / dt); ++i) {
        const double t = i * dt;
        const double acceleration =
            -a1 * w1 * w1 * std::sin(w1 * t)
            -a2 * w2 * w2 * std::sin(w2 * t + 0.7);
        estimator.update(acceleration, dt);
    }

    // Include the known leaky-integrator transfer in the expected moment ratio.
    const double eta1_sq = a1 * a1 * std::pow(w1 * w1 / (w1 * w1 + lambda * lambda), 2);
    const double eta2_sq = a2 * a2 * std::pow(w2 * w2 / (w2 * w2 + lambda * lambda), 2);
    const double expected_omega = std::sqrt(
        ((w1 * w1 + lambda * lambda) * eta1_sq
       + (w2 * w2 + lambda * lambda) * eta2_sq)
        / (eta1_sq + eta2_sq));
    const double expected_period = 2.0 * std::numbers::pi / expected_omega;
    const double relative_error = std::abs(estimator.getPeriodSec() - expected_period) / expected_period;
    return check(relative_error < 0.02, "multitone spectral-moment period mismatch");
}

bool test_rejects_bad_samples_and_resets() {
    WaveBandPeriodEstimator<double> estimator;
    estimator.update(1.0, 0.0);
    estimator.update(NAN, 0.01);
    if (!check(estimator.elapsedSec() == 0.0, "invalid samples changed estimator state")) return false;
    estimator.update(1.0, 0.01);
    estimator.reset();
    return check(estimator.elapsedSec() == 0.0 && !estimator.isReady(), "reset failed");
}

} // namespace

int main() {
    if (!test_single_tone()) return 1;
    if (!test_multitone_zero_crossing_period()) return 1;
    if (!test_rejects_bad_samples_and_resets()) return 1;
    return 0;
}
