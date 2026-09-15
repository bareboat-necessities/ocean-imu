// Unchanged shipping, ordinary startup and one finite post-Live gyro pulse.
// No state injection, setter, magnetic shortcut or covariance intervention.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std = 9.80665f;

int main() {
    using Wrapper = SeaStateFusion_OU_III<TrackerType::KALMANF>;
    Wrapper wrapper;
    Wrapper::Config config;
    wrapper.begin(config);
    const Eigen::Vector3f acc(0, 0, -g_std);
    int startup_steps = 0;
    while (!wrapper.isLive() && startup_steps < 30602) {
        wrapper.update(.005f, Eigen::Vector3f::Zero(), acc);
        ++startup_steps;
    }
    if (!wrapper.isLive()) return 2;
    auto emit = [&](int post_live_steps) {
        std::cout << startup_steps << ' ' << post_live_steps << ' '
                  << wrapper.isLive() << ' '
                  << wrapper.attitudeQuat().coeffs().allFinite() << ' '
                  << std::isfinite(wrapper.raw().getAccelVertical()) << ' '
                  << (wrapper.raw().accelVibrationGuardEngagement() == 0.0f) << ' '
                  << wrapper.raw().startupProxyQuat().coeffs().allFinite() << '\n';
    };
    emit(0);
    for (int k = 1; k <= 600; ++k) {
        const Eigen::Vector3f gyro = k == 1
            ? Eigen::Vector3f(0x1p80f, 0, 0) : Eigen::Vector3f::Zero();
        if (!gyro.allFinite() || !acc.allFinite()) return 3;
        wrapper.update(.005f, gyro, acc);
        if (k == 1 || k == 2 || k == 600) emit(k);
    }
}
