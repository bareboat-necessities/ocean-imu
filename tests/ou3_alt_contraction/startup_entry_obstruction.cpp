// Execute unchanged shipping from construction; no installed state or replay.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cstdint>
#include <cstring>
#include <iostream>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

static std::uint32_t bits(float value) {
    std::uint32_t word;
    std::memcpy(&word, &value, sizeof(word));
    return word;
}

int main() {
    using Wrapper = SeaStateFusion_OU_III<TrackerType::KALMANF>;
    Wrapper wrapper;
    Wrapper::Config config;
    wrapper.begin(config);
    for (int step = 1; step <= 30200; ++step) {
        const auto before = wrapper.raw().getStartupStage();
        // Identical physical IMU history for every constant heading of a level
        // boat at rest. No pre-Live magnetic call is required by ALT's schedule.
        wrapper.update(.005f, Eigen::Vector3f::Zero(), Eigen::Vector3f(0, 0, -g_std));
        if (wrapper.isLive()) {
            const auto q = wrapper.attitudeQuat();
            std::cout << step << ' ' << bits(wrapper.liveTimeSec()) << ' '
                      << int(before) << ' ' << int(wrapper.raw().getStartupStage()) << ' '
                      << wrapper.hasMagNorthLock() << ' '
                      << bits(q.w()) << ' ' << bits(q.x()) << ' '
                      << bits(q.y()) << ' ' << bits(q.z()) << '\n';
            return 0;
        }
    }
    return 2;
}
