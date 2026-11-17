#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
namespace ocean_imu::ahrs {
// Continuous gyro-bias random-walk variance density, rad^2/s^3.
// This is the compass firmware's 1e-9 variance per 5 ms exposure.
inline constexpr float qmekf_gyro_bias_rw_density = 2.0e-7f;
}
