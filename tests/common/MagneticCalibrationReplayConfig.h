#pragma once
#include <cstdlib>
#include <stdexcept>
#include <string_view>

namespace magnetic_replay {
template<class Config> void configure(Config& cfg) {
    const char* value = std::getenv("W3D_MAG_CALIBRATION");
    if (!value) return;
    const std::string_view mode(value);
    if (mode != "off" && mode != "native" && mode != "screened")
        throw std::invalid_argument("W3D_MAG_CALIBRATION must be off, native or screened");
    cfg.mag_continuous_hard_iron = mode != "off";
    cfg.mag_screen_installation = mode == "screened";
}
}
