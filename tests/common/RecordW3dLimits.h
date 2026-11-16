#pragma once
#include "RecordReplayLimits.h"
#include "GyroReplayContract.h"
#include "util/W3dSimCommon.h"

namespace ocean_imu::sim::record_limits {
inline W3dFailureLimits w3d(std::string_view family, std::string_view profile,
                           std::string_view input, W3dFailureLimits bars) {
    // Seed sweeps and explicitly enabled regular-wave controls are diagnostic
    // experiments outside the eight-record fixed-draw protocol.
    if (!uses_record_caps(input)) return bars;
    const auto& e = find(family, profile, input);
    bars.err_limit_percent_z_jonswap = tighten(e, Metric::heave, bars.err_limit_percent_z_jonswap);
    bars.err_limit_percent_z_pmstokes = tighten(e, Metric::heave, bars.err_limit_percent_z_pmstokes);
    bars.err_limit_percent_3d_jonswap = tighten(e, Metric::displacement, bars.err_limit_percent_3d_jonswap);
    bars.err_limit_percent_3d_pmstokes = tighten(e, Metric::displacement, bars.err_limit_percent_3d_pmstokes);
    bars.err_limit_yaw_deg = tighten(e, Metric::yaw, bars.err_limit_yaw_deg);
    if (bars.err_limit_roll_deg > 0) bars.err_limit_roll_deg = tighten(e, Metric::roll, bars.err_limit_roll_deg);
    if (bars.err_limit_pitch_deg > 0) bars.err_limit_pitch_deg = tighten(e, Metric::pitch, bars.err_limit_pitch_deg);
    bars.acc_z_bias_percent = tighten(e, Metric::acc_z, bars.acc_z_bias_percent);
    bars.gyro_bias_3d_percent = tighten(e, Metric::gyro_3d,
        bars.gyro_bias_3d_percent > 0 ? bars.gyro_bias_3d_percent : bars.bias_3d_percent);
    bars.bias_3d_percent = tighten(e, Metric::acc_3d, bars.bias_3d_percent);
    return bars;
}
} // namespace ocean_imu::sim::record_limits
