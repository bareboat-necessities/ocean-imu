#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Periodic pseudo-measurement schedule shared by the wave MEKFs.
// KalmanOUCoreMath.h re-exports these names into ou_detail.

#include <algorithm>
#include <cmath>
#include <limits>

namespace ocean_imu::kalman::filter_detail {

// Retarget a periodic scheduler without discarding elapsed service credit.
// If the new deadline has already passed, park immediately below it so the
// next valid sample services the owed update through periodic_update_due().
// Otherwise preserve elapsed bit-for-bit.
template<typename T>
inline T retarget_period_elapsed_progress_preserving(T elapsed, T period) {
    if (!(period > T(0)) || !std::isfinite(period)) return T(0);
    if (!(elapsed >= T(0)) || !std::isfinite(elapsed)) return T(0);
    if (elapsed < period) return elapsed;
    return std::nextafter(period, T(0));
}

template<typename T>
inline bool periodic_update_due(T dt, T period, T& elapsed) {
    if (!(dt > T(0)) || !std::isfinite(dt) || !(period > T(0)) || !std::isfinite(period)) return false;
    const T total = elapsed + dt;
    const T tol = T(16) * std::numeric_limits<T>::epsilon() * std::max(T(1), period);
    if (total + tol < period) {
        elapsed = total;
        return false;
    }
    elapsed = (total >= period) ? std::fmod(total, period) : T(0);
    if (!(elapsed >= T(0)) || !std::isfinite(elapsed) || elapsed >= period) elapsed = T(0);
    return true;
}

} // namespace ocean_imu::kalman::filter_detail
