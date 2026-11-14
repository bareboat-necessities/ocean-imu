#pragma once
#include <cstdlib>
#include <stdexcept>
#include <string>
#include "kalman_common/GyroInterval.h"

// Qualify both integrators independently on the pinned underway experiment.
// Ordinary replays exercise the live Config defaults and their own fixture.
namespace gyro_replay {
inline bool legacy = false;
template<class Config> void configure(Config& cfg) {
    if (legacy) {
        using Mode = ocean_imu::kalman::filter_detail::GyroSampleMode;
        cfg.gyro_sample_mode = Mode::IntervalRate;
        cfg.frontend_gyro_sample_mode = Mode::IntervalRate;
    }
}
inline bool default_seeds() {
    return !std::getenv("W3D_SEED") && !std::getenv("W3D_IMU_SEED") &&
           !std::getenv("W3D_INIT_SEED");
}
inline void require_regression(bool default_protocol) {
    const char* reference = std::getenv("W3D_REFERENCE_MOTION");
    if (legacy && reference && std::string(reference) != "underway")
        throw std::invalid_argument("--regression requires the pinned underway reference");
    if (legacy && !default_protocol)
        throw std::invalid_argument("--regression requires default seeds, noise, magnetometer and variant");
}
inline const char* name() { return legacy ? "legacy_interval" : "instantaneous_default"; }
}
