#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Combine an undirected wave-propagation axis (modulo 180 degrees) with the
  WaveDirectionDetector sense into the values the marine INS sketches publish.
*/

#include <cmath>

#include "util/AngleUtils.h"
#include "wave_dir/WaveDirectionDetector.h"

namespace wave_direction {

// -1: propagation opposite the axis (add 180 deg), +1: along it, 0: unknown.
inline int signPolarity(WaveDirection s) {
  if (s == UNCERTAIN) {
    return 0;
  }

  const int raw = static_cast<int>(s);

  if (raw < 0) return -1;
  if (raw > 0) return +1;

  return 0;
}

inline int signRaw(WaveDirection s) {
  return static_cast<int>(s);
}

// Resolve the axis into a direction in [0, 360).  Returns false and sets NaN
// when the axis is not finite or the sense is unknown.
inline bool signedDirectionDeg(float axis_deg, int sign_polarity, float& signed_deg_out) {
  if (!std::isfinite(axis_deg) || sign_polarity == 0) {
    signed_deg_out = NAN;
    return false;
  }

  signed_deg_out = ocean_imu::ins::wrap360Deg(axis_deg + (sign_polarity < 0 ? 180.0f : 0.0f));
  return true;
}

struct WaveDirectionReport {
  float         axis_deg      = NAN;
  bool          axis_ok       = false;
  WaveDirection sign          = UNCERTAIN;
  int           sign_raw      = 0;
  int           sign_polarity = 0;
  bool          sign_ok       = false;
  float         dir_deg       = NAN;
  bool          dir_ok        = false;
  float         conf_pct      = 0.0f;

  // Every "ok" flag is gated on `live`: nothing is published as valid before
  // the filter that produced the axis and sense is running on fused attitude.
  // conf_pct is 100 with a resolved direction, 50 with only the axis, else 0.
  static WaveDirectionReport from(float axis_deg, WaveDirection sign, bool live) {
    WaveDirectionReport r;
    r.axis_deg      = axis_deg;
    r.axis_ok       = live && std::isfinite(axis_deg);
    r.sign          = sign;
    r.sign_raw      = signRaw(sign);
    r.sign_polarity = signPolarity(sign);
    r.sign_ok       = live && (sign != UNCERTAIN);
    r.dir_ok        = signedDirectionDeg(axis_deg, r.sign_polarity, r.dir_deg) && live;
    r.conf_pct      = r.dir_ok ? 100.0f : (r.axis_ok ? 50.0f : 0.0f);
    return r;
  }
};

}  // namespace wave_direction
