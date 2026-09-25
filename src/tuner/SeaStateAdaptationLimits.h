#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

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

// Final guard on dynamically derived EMA horizons.  It prevents both
// sample-to-sample pass-through and a nearly frozen adaptation after an
// estimator excursion.
//
// With K_periods = 4, the acceleration-moment EMA requests 4*T_z, about
// 33.6 s on the longest reference sea.  The 35 s ceiling leaves that request
// interior.  Changes to K_periods or these limits alter adaptation memory.
inline constexpr float kDynamicEmaHorizonMinSec = 0.05f;
inline constexpr float kDynamicEmaHorizonMaxSec = 35.0f;

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
