#pragma once

#include <algorithm>
#include <cmath>

namespace seastate::tuner::limits {

// Universal safety envelope for a dynamically estimated sea/model time scale
// when it is used only to choose an EMA horizon.  This does NOT clamp the
// physical wave-period estimate or the OU process time constant themselves.
//
// The eight versioned JONSWAP + PM/Stokes reference seas span T_z ~= 2.3..8.4 s,
// hence T_sea=T_z/2 ~= 1.15..4.2 s.  [0.5, 6] s therefore leaves every
// calibrated sea comfortably interior while rejecting transient/degenerate
// period estimates outside a physically useful marine adaptation envelope.
inline constexpr float kDynamicEmaTimeScaleMinSec = 0.5f;
inline constexpr float kDynamicEmaTimeScaleMaxSec = 6.0f;

// Final guard on dynamically derived EMA horizons.  The deployed reference
// horizons span about 0.12..16.8 s (and about 12.7 s for the slowest OU-II/TFG
// drift smoother), so [0.05, 30] s is deliberately inactive on all eight
// calibration seas while preventing either sample-to-sample pass-through or a
// nearly frozen adaptation after an estimator excursion.
inline constexpr float kDynamicEmaHorizonMinSec = 0.05f;
inline constexpr float kDynamicEmaHorizonMaxSec = 30.0f;

inline float clampDynamicEmaTimeScaleSec(float sec) noexcept {
    if (!(std::isfinite(sec) && sec > 0.0f)) return kDynamicEmaTimeScaleMaxSec;
    return std::clamp(sec, kDynamicEmaTimeScaleMinSec, kDynamicEmaTimeScaleMaxSec);
}

inline float clampDynamicEmaHorizonSec(float sec, float dt_sec = 0.0f) noexcept {
    if (!(std::isfinite(sec) && sec > 0.0f)) sec = kDynamicEmaHorizonMaxSec;
    float lo = kDynamicEmaHorizonMinSec;
    if (std::isfinite(dt_sec) && dt_sec > lo) {
        lo = std::min(dt_sec, kDynamicEmaHorizonMaxSec);
    }
    return std::clamp(sec, lo, kDynamicEmaHorizonMaxSec);
}

}  // namespace seastate::tuner::limits
