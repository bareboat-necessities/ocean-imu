#pragma once
#include <cstdlib>
#include <stdexcept>
#include <string_view>
#include "kalman_common/SeaStateFusionDefaults.h"

namespace magnetic_replay {
template<class Config> void configure(Config& cfg) {
    const char* value = std::getenv("W3D_MAG_CALIBRATION");
    if (!value) return;
    const std::string_view mode(value);
    if (mode != "off" && mode != "native" && mode != "screened" && mode != "production")
        throw std::invalid_argument("W3D_MAG_CALIBRATION must be off, native, screened or production");
    if (mode == "production") return;
    if constexpr (requires { cfg.mag_hi_heading_transport_fraction; }) {
        cfg.mag_hi_heading_transport_fraction = 0.f;
        cfg.mag_hi_slew_tau_sec = seastate::common::defaults::MAG_HI_SLEW_TAU_SEC;
    }
    cfg.mag_continuous_hard_iron = mode != "off";
    cfg.mag_screen_installation = mode == "screened";
}
}
