#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Test-only limits: one pinned sea must not borrow another sea's error budget.
#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <limits>
#include <string>
#include <string_view>
#include "GyroReplayContract.h"
#include "util/UnderwayReference.h"

namespace ocean_imu::sim::record_limits {
enum class Metric : std::size_t { heave, displacement, roll, pitch, yaw, acc_z, acc_3d, gyro_3d, displacement_y };
struct Entry {
    std::string_view scenario, family, profile, input;
    std::array<float, 9> limits;
};
#include "RecordReplayLimits.inc"

inline bool uses_record_caps(std::string_view input) {
    if (!gyro_replay::default_seeds()) return false;
    input.remove_prefix(input.find_last_of("/\\") == std::string_view::npos
                        ? 0 : input.find_last_of("/\\") + 1);
    const char* all = std::getenv("W3D_ALL_WAVE_TYPES");
    // Startup/engine studies explicitly include these regular-wave controls.
    // They retain diagnostic global bars; the eight irregular-wave fixtures
    // and unknown/misspelled inputs must still resolve a locked record cap.
    const bool regular = input.starts_with("wave_data_gerstner_") ||
                         input.starts_with("wave_data_cnoidal_") ||
                         input.starts_with("wave_data_fenton_");
    return !(all && std::string_view(all) != "0" && regular);
}

inline const Entry& find(std::string_view family, std::string_view profile, std::string_view input) {
    input.remove_prefix(input.find_last_of("/\\") == std::string_view::npos
                        ? 0 : input.find_last_of("/\\") + 1);
    for (const auto& entry : entries)
        if (entry.scenario == UnderwayReference::scenario() && entry.family == family && entry.profile == profile && entry.input == input)
            return entry;
    throw std::invalid_argument("missing per-record gate: " + std::string(family) + "/" +
                                std::string(profile) + "/" + std::string(input));
}

inline float tighten(const Entry& entry, Metric metric, float global) {
    const float local = entry.limits.at(static_cast<std::size_t>(metric));
    if (!(local > 0.0f)) throw std::invalid_argument("unmeasured per-record gate channel");
    return std::min(global, local);
}

inline float limit(std::string_view family, std::string_view profile, std::string_view input,
                   Metric metric, float global) {
    if (!uses_record_caps(input)) return global;
    // Default protocol bars belong to the selected v1.4.0 trajectory. The
    // historical fallback is only for diagnostic/non-default sensor draws.
    return tighten(find(family, profile, input), metric, std::numeric_limits<float>::infinity());
}
} // namespace ocean_imu::sim::record_limits
