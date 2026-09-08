#include <cassert>
#include <cmath>
#include <iostream>

#include "freq/PLLFreqTracker.h"

// A seconds-based acquisition smoother must track a period change in the
// same physical time at each IMU sampling rate.
static double replay(double sample_rate_hz) {
    PLLFreqTracker<double> tracker;
    PLLFreqTracker<double>::Config config;
    config.f_min_hz = 0.04;
    config.f_max_hz = 2.0;
    config.f_init_hz = 0.25;
    config.coarse_smooth_tau_s = 4.5;
    tracker.configure(config);
    const double dt = 1.0 / sample_rate_hz;
    double phase = 0.0;
    for (int i = 0; i < static_cast<int>(120.0 * sample_rate_hz); ++i) {
        const double frequency = i * dt < 40.0 ? 0.25 : 0.10;
        phase += 2.0 * std::acos(-1.0) * frequency * dt;
        tracker.update(std::sin(phase), dt);
    }
    assert(tracker.hasCoarseEstimate());
    return tracker.getCoarseFrequencyHz();
}

int main() {
    const double slow = replay(50.0);
    const double nominal = replay(200.0);
    const double fast = replay(400.0);
    std::cout << "Coarse frequency after 0.25 -> 0.10 Hz: "
              << slow << ' ' << nominal << ' ' << fast << '\n';
    for (const double estimate : {slow, nominal, fast}) {
        assert(std::abs(estimate - 0.10) < 0.002);
    }
    assert(std::abs(slow - fast) < 0.001);
}
