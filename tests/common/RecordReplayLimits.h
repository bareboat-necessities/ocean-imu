#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Test-only limits: one pinned sea must not borrow another sea's error budget.
#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include "GyroReplayContract.h"

namespace ocean_imu::sim::record_limits {
enum class Metric : std::size_t { heave, displacement, roll, pitch, yaw, acc_z, acc_3d, gyro_3d };
struct Entry {
    std::string_view family, profile, input;
    std::array<float, 8> limits;
};
#include "RecordReplayLimits.inc"

inline const Entry& find(std::string_view family, std::string_view profile, std::string_view input) {
    input.remove_prefix(input.find_last_of("/\\") == std::string_view::npos
                        ? 0 : input.find_last_of("/\\") + 1);
    for (const auto& entry : entries)
        if (entry.family == family && entry.profile == profile && entry.input == input)
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
    if (!gyro_replay::default_seeds()) return global;
    return tighten(find(family, profile, input), metric, global);
}
} // namespace ocean_imu::sim::record_limits
