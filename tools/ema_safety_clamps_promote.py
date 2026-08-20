#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"expected one match in {path}, found {count}: {old[:80]!r}")
    path.write_text(text.replace(old, new, 1))


limits = ROOT / "src/tuner/SeaStateAdaptationLimits.h"
limits.write_text(r'''#pragma once

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
''')

# Shared tuner: the frequency EMA and K-period variance EMA both derive their
# horizons from the same clamped sea time, then receive the common final guard.
tuner = ROOT / "src/tuner/SeaStateAutoTuner.h"
replace_once(
    tuner,
    '#include <algorithm>\n',
    '#include <algorithm>\n\n#include "tuner/SeaStateAdaptationLimits.h"\n',
)
replace_once(
    tuner,
    '''        const float T_eff = 1.0f / f_eff;\n\n        // Dynamic variance time constant: a few periods\n        tau_var_sec = std::max(tau_var_min_sec,\n                               std::min(tau_var_max_sec, K_periods * T_eff));\n''',
    '''        const float sea_time_sec =\n            seastate::tuner::limits::clampDynamicEmaTimeScaleSec(0.5f / f_eff);\n        const float T_eff = 2.0f * sea_time_sec;  // T_z = 2 T_sea\n\n        // Dynamic variance time constant: a few wave periods.  The local\n        // bounds remain useful for controlled studies, while the universal\n        // guard prevents any dynamically estimated horizon from escaping the\n        // shared safety envelope.\n        const float tau_var_requested = std::max(\n            tau_var_min_sec, std::min(tau_var_max_sec, K_periods * T_eff));\n        tau_var_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(\n            tau_var_requested, dt_s);\n''',
)
replace_once(
    tuner,
    '''        float horizon = tau_freq;\n        if (freq_sea_periods > 0.0f) {\n            float f = f_input_hz;\n            if (!std::isfinite(f)) f = f_min_hz;\n            f = std::max(f_min_hz, std::min(f_max_hz, f));\n            horizon = freq_sea_periods * (0.5f / f);\n        }\n        tau_freq_applied_sec = std::max(1e-3f, horizon);\n        alpha_freq = 1.0f - std::exp(-dt_s / tau_freq_applied_sec);\n''',
    '''        float horizon = tau_freq;\n        if (freq_sea_periods > 0.0f) {\n            float f = f_input_hz;\n            if (!std::isfinite(f)) f = f_min_hz;\n            f = std::max(f_min_hz, std::min(f_max_hz, f));\n            const float sea_time_sec =\n                seastate::tuner::limits::clampDynamicEmaTimeScaleSec(0.5f / f);\n            horizon = seastate::tuner::limits::clampDynamicEmaHorizonSec(\n                freq_sea_periods * sea_time_sec, dt_s);\n            tau_freq_applied_sec = horizon;\n        } else {\n            // Explicit fixed-seconds compatibility mode is not dynamically\n            // estimated, so preserve the caller's requested ablation value.\n            tau_freq_applied_sec = std::max(1e-3f, horizon);\n        }\n        alpha_freq = 1.0f - std::exp(-dt_s / tau_freq_applied_sec);\n''',
)

# Shared drift-channel smoother: tau is the dynamically estimated model time
# base for r_S/r_p0/r_v0. Clamp it for EMA purposes only, then guard the final
# horizon after the optional discrepancy shortening.
common = ROOT / "src/kalman_common/SeaStateFusionFilterCommon.h"
replace_once(
    common,
    '#include "tuner/SeaStateFusionTunerCommon.h"\n',
    '#include "tuner/SeaStateFusionTunerCommon.h"\n#include "tuner/SeaStateAdaptationLimits.h"\n',
)
replace_once(
    common,
    '    float horizon = mult * tau_sec;\n',
    '''    const float safe_time_scale =\n        seastate::tuner::limits::clampDynamicEmaTimeScaleSec(tau_sec);\n    float horizon = mult * safe_time_scale;\n''',
)
replace_once(
    common,
    '''    // The horizon is proportional to the operating point, so a degenerate tau\n    // must not collapse it to zero and turn the EMA into a pass-through.  One\n    // step is the shortest meaningful horizon.\n    return std::max(horizon, dt);\n''',
    '''    // The dynamic time base and the final EMA horizon have independent\n    // guards: the latter also protects the optional discrepancy gate above\n    // from collapsing a legitimate sea-scaled horizon into sample-to-sample\n    // pass-through.\n    return seastate::tuner::limits::clampDynamicEmaHorizonSec(horizon, dt);\n''',
)

# All three wrapper-specific common tau/sigma EMAs use the same clamped T_sea.
for rel in [
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
]:
    path = ROOT / rel
    replace_once(
        path,
        '            adapt_sec = std::max(dt, adapt_tau_sea_periods_ * sea_time_sec);\n',
        '''            const float safe_sea_time =\n                seastate::tuner::limits::clampDynamicEmaTimeScaleSec(sea_time_sec);\n            adapt_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(\n                adapt_tau_sea_periods_ * safe_sea_time, dt);\n''',
    )

tfg = ROOT / "src/kalman_tfg/SeaStateFusionFilter_TFG.h"
replace_once(
    tfg,
    '                adapt_sec = std::max(dt, adapt_tau_sea_periods_ * sea_time_sec_);\n',
    '''                const float safe_sea_time =\n                    seastate::tuner::limits::clampDynamicEmaTimeScaleSec(sea_time_sec_);\n                adapt_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(\n                    adapt_tau_sea_periods_ * safe_sea_time, dt);\n''',
)

# Pin the universal contract and verify that the full observed reference range
# remains strictly inside it.  Also exercise both pathological sides.
test = ROOT / "tests/kalman_ou_ii/tuner_schedule-test.cpp"
needle = '''    if (!near(f.getAdaptationSeaPeriods(), 0.40f) ||\n        !near(f.getTunerFreqSmoothingSeaPeriods(), 0.10f)) {\n        std::cerr << "FAIL: OU-II sea-scaled EMA defaults are not 0.40/0.10 of T_sea\\n";\n        return 1;\n    }\n'''
addition = needle + r'''

    // Universal dynamic-EMA safety envelope.  The complete eight-record
    // reference family spans T_z ~= 2.3..8.4 s, hence T_sea ~= 1.15..4.2 s;
    // both ends must pass through unchanged rather than riding a clamp.
    using namespace seastate::tuner::limits;
    if (!near(clampDynamicEmaTimeScaleSec(1.15f), 1.15f) ||
        !near(clampDynamicEmaTimeScaleSec(4.20f), 4.20f)) {
        std::cerr << "FAIL: reference sea-time envelope touches the safety clamp\n";
        return 1;
    }
    if (!near(clampDynamicEmaTimeScaleSec(0.01f), kDynamicEmaTimeScaleMinSec) ||
        !near(clampDynamicEmaTimeScaleSec(100.0f), kDynamicEmaTimeScaleMaxSec)) {
        std::cerr << "FAIL: dynamic EMA time-scale clamp does not catch excursions\n";
        return 1;
    }
    if (!near(seastate::common::adaptiveSmoothingHorizonSec(
                  0.001f, 0.001f, 1.0f, 1.0f, 0.0f, 0.005f),
              kDynamicEmaHorizonMinSec) ||
        !near(seastate::common::adaptiveSmoothingHorizonSec(
                  100.0f, 100.0f, 1.0f, 1.0f, 0.0f, 0.005f),
              kDynamicEmaHorizonMaxSec)) {
        std::cerr << "FAIL: dynamic EMA final horizon guard is not universal\n";
        return 1;
    }

    // The shared tuner must use the same guards.  At absurdly short/long input
    // periods the deployed 0.10*T_sea frequency EMA becomes 0.05/0.60 s, and
    // the K=2 variance horizon becomes 2/24 s.  These values are far outside
    // the eight calibrated seas but finite and deliberately bounded.
    {
        SeaStateAutoTuner fast(2.0f, 1.0f);
        fast.setFrequencySmoothingSeaPeriods(0.10f);
        fast.update(0.005f, 0.0f, 100.0f);
        if (!near(fast.getFrequencySmoothingHorizonSec(), 0.05f) ||
            !near(fast.getVarianceHorizonSec(), 2.0f)) {
            std::cerr << "FAIL: short-period tuner horizons escaped universal clamps\n";
            return 1;
        }

        SeaStateAutoTuner slow(2.0f, 1.0f);
        slow.setFrequencySmoothingSeaPeriods(0.10f);
        slow.update(0.005f, 0.0f, 0.001f);
        if (!near(slow.getFrequencySmoothingHorizonSec(), 0.60f) ||
            !near(slow.getVarianceHorizonSec(), 24.0f)) {
            std::cerr << "FAIL: long-period tuner horizons escaped universal clamps\n";
            return 1;
        }
    }
'''
replace_once(test, needle, addition)

# Engineering note: explicitly show why the clamps are not part of the fit.
doc = ROOT / "docs/adaptive-ema-safety-clamps.md"
doc.write_text(r'''# Universal safety clamps for sea-scaled adaptation EMAs

The OU-III, OU-II and TFG adaptation paths now share one safety envelope for
**dynamically estimated EMA time scales**.  This is a guard, not another fitted
tuning coefficient.

## Bounds

- dynamic sea/model time base used by an EMA: **0.5 to 6.0 s**
- final dynamically derived EMA horizon: **0.05 to 30.0 s**

The time-base clamp is applied only when a measured `T_sea` or dynamically
estimated OU `tau` is used to construct an EMA horizon.  It does **not** clamp
the physical wave-period estimate or the OU process time constant itself.
Explicit fixed-second ablation setters keep their requested values.

## Why these limits are outside the calibrated envelope

The versioned eight-sea reference family spans approximately `T_z = 2.3..8.4 s`,
therefore `T_sea = T_z/2 = 1.15..4.2 s`.  The deployed dynamic horizons are:

| channel | law | reference range |
| --- | --- | ---: |
| tuner frequency EMA | `0.10 T_sea` | 0.115..0.42 s |
| common `tau/sigma_aw` EMA | `0.40 T_sea` | 0.46..1.68 s |
| tuner variance EMA (`K=2`) | `2 T_z = 4 T_sea` | 4.6..16.8 s |
| OU-III `r_S` EMA | `1.5 tau` | about 1.9..6.3 s |
| OU-II `r_p0/r_v0` EMA | `3 tau` | about 3.8..12.7 s |
| TFG `r_S` EMA | `3 tau` | same order, below 13 s on the reference family |

Thus every deployed reference trajectory is strictly interior to both guards.
The clamps act only on startup/transient estimator excursions or on inputs well
outside the evidenced marine envelope.

## Why clamping only `T_sea` is almost, but not quite, enough

Clamping the dynamic sea time automatically protects the tuner frequency EMA,
the common `tau/sigma_aw` EMA and the K-period variance EMA.  The drift-channel
smoothers need one additional final guard because their helper can optionally
shorten its horizon with the discrepancy/slew term *after* the sea/model time
has been scaled.  The shared 0.05..30 s outer guard catches that path too.

All three filter families call the same helpers, so the limits cannot silently
diverge between OU-III, OU-II and TFG.
''')

print("Applied universal dynamic EMA safety clamps")
