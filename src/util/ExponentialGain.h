#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#include <cmath>
#include <type_traits>

namespace ocean_imu {

// Stable first-order gain 1 - exp(-x). At IMU cadence the float series
// avoids a libm call. For x <= 1/128 the first omitted term x^4/24
// contributes less than 2.0e-8 relative error. Keep the full range and
// precision of expm1 for longer steps and double-precision callers.
template<typename T>
inline T exponential_gain(T x) {
    if constexpr (std::is_same_v<T,float>) {
        if (x >= 0.0f && x <= 1.0f/128.0f)
            return x*(1.0f+x*(-0.5f+x*(1.0f/6.0f)));
    }
    return -std::expm1(-x);
}

}  // namespace ocean_imu
