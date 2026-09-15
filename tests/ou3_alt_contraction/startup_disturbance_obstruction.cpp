// Public shipping API only: no state installation or filter configuration change.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std = 9.80665f;

int main() {
    using Wrapper = SeaStateFusion_OU_III<TrackerType::KALMANF>;
    constexpr int steps = 30602;
    for (bool service_mag : {false, true}) {
        Wrapper wrapper;
        Wrapper::Config config;
        wrapper.begin(config);
        int calls = 0;
        for (int k = 1; k <= steps; ++k) {
            // Quiet physical wave and zero bias, plus one bounded residual
            // n_a=(0,0,g-2^-11). The measured direction remains exactly up.
            wrapper.update(.005f, Eigen::Vector3f::Zero(),
                           Eigen::Vector3f(0, 0, -0x1p-11f));
            if (service_mag && k % 8 == 0) {
                // The same stationary physical field, with zero magnetic noise.
                wrapper.updateMag(Eigen::Vector3f(15, 0, 20));
                ++calls;
            }
            if (wrapper.raw().startupProxyInitialized() || wrapper.isLive() ||
                wrapper.raw().accelVibrationGuardEngagement() != 0.0f ||
                wrapper.raw().accelVibrationRms() != 0.0f) return 2;
        }
        const bool north = wrapper.hasMagNorthLock();
        std::cout << service_mag << ' ' << steps << ' ' << calls << ' '
                  << wrapper.raw().startupProxyInitialized() << ' '
                  << wrapper.isLive() << ' ' << north << ' ';
        // One informative sample isolates the seed predicate as the cause.
        wrapper.update(.005f, Eigen::Vector3f::Zero(), Eigen::Vector3f(0, 0, -g_std));
        std::cout << wrapper.raw().startupProxyInitialized() << ' '
                  << wrapper.isLive() << '\n';
    }
}
