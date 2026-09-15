// Host saturation correspondence; no billions-of-calls replay is needed.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

int main() {
    static_assert(std::numeric_limits<int>::max() == 2147483647);
    const auto maximum = std::numeric_limits<int>::max();
    MagAutoTuner tuner;
    MagAutoTuner::Config tc;
    tc.min_samples = maximum;
    tc.min_window_sec = 1000000.0f;
    tuner.setConfig(tc);
    tuner.accepted_count_ = maximum - 1;
    tuner.rejected_count_ = maximum - 1;
    for (int k = 0; k < 3; ++k) {
        tuner.addSampleWithTiltQuatDt(.005f, Eigen::Quaternionf::Identity(),
            Eigen::Vector3f(0,0,-g_std), Eigen::Vector3f::Zero(), Eigen::Vector3f(15,0,20));
        tuner.addSampleWithTiltQuatDt(.005f, Eigen::Quaternionf::Identity(),
            Eigen::Vector3f(0,0,-g_std), Eigen::Vector3f::Zero(), Eigen::Vector3f::Zero());
    }
    if (tuner.acceptedCount() != maximum || tuner.rejectedCount() != maximum) return 9;
    if (tuner.weight_sum_ != 3.0f) return 10; // accepted samples still accumulate
    Fusion f;
    f.begin(Fusion::Config{});
    const Eigen::Vector3f field(15.0f, 0.0f, 20.0f);
    int samples = 0;
    for (; samples < 40000 && !(f.isLive() && f.hasMagNorthLock()); ++samples) {
        f.update(.005f, Eigen::Vector3f::Zero(), Eigen::Vector3f(0, 0, -g_std));
        if (samples % 8 == 0) f.updateMag(field);
    }
    if (!(f.isLive() && f.hasMagNorthLock())) {
        std::cerr << "quiet public-API prefix did not reach gauged Live\n";
        return 2;
    }
    // Install the boundary count, then exercise the actual outer API.
    auto& inner = f.raw();
    inner.setMagUpdatesToUnlockAccBias(std::numeric_limits<int>::max());
    inner.setAccBiasHold(true);
    inner.mag_updates_applied_ = std::numeric_limits<int>::max() - 1;
    f.updateMag(field);
    if (inner.mag_updates_applied_ != std::numeric_limits<int>::max()) return 3;
    const auto before = inner.mekf().covariance_full();
    f.updateMag(field);
    if (inner.mag_updates_applied_ != std::numeric_limits<int>::max()) return 4;
    if ((before.array() == inner.mekf().covariance_full().array()).all()) {
        std::cerr << "measurement stopped at counter saturation\n";
        return 5;
    }
    for (int k = 0; k < 250; ++k) {
        f.update(.005f, Eigen::Vector3f::Zero(), Eigen::Vector3f(0, 0, -g_std));
        if (k % 8 == 0) f.updateMag(field);
    }
    f.updateMag(field);
    if (inner.mag_updates_applied_ != std::numeric_limits<int>::max()) return 6;
    if (inner.accel_bias_locked_ || !inner.acc_bias_hold_) return 7;
    inner.setAccBiasHold(false);
    if (inner.acc_bias_hold_) return 8;
    std::cout << "gauged_Live_samples=" << samples
              << " saturated_count=" << inner.mag_updates_applied_
              << " measurement_continues=1 release_continues=1\n";
    return 0;
}
