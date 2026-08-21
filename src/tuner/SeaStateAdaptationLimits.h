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

// Final guard on dynamically derived EMA horizons.  It prevents both
// sample-to-sample pass-through and a nearly frozen adaptation after an
// estimator excursion.
//
// This ceiling used to be inactive on every calibration sea: the longest
// deployed horizon was the sigma_a variance EMA at 2*T_z, i.e. about 16.8 s on
// the largest reference sea.  Since the period-statistics retune raised
// K_periods to 4 that horizon requests 4*T_z, about 33.6 s there, so the 30 s
// ceiling now binds on the largest of the eight seas and only on that channel.
// Every other deployed horizon (about 0.12..12.7 s, the slowest being the
// OU-II/TFG drift smoother) is still comfortably interior.  Treat a change to
// K_periods or to this ceiling as a retune that needs a re-gauge, not a tweak;
// see docs/adaptive-ema-safety-clamps.md.
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
