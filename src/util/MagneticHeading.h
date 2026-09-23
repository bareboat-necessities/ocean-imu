#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Magnetometer helpers shared by the AtomS3R marine INS sketches.

  Body frame is the project body-NED convention: x forward, y starboard,
  z down.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <stdint.h>

#include "util/AngleUtils.h"

namespace ocean_imu {
namespace ins {

// Tilt-compensated magnetic heading of the bow, in [0, 360).
//
// down_b_unit is the world down direction expressed in BODY (any nonzero
// length); mag_b_uT is the calibrated body magnetic field.  Only the tilt is
// taken from down_b_unit, so an attitude with a wrong or unobservable yaw still
// gives the correct heading.  Returns false and sets the output to NaN when the
// geometry is degenerate or an input is not finite.
inline bool magneticHeadingFromDownAndMagBody(const Eigen::Vector3f& down_b_unit,
                                              const Eigen::Vector3f& mag_b_uT,
                                              float& heading_deg_out) {
  heading_deg_out = NAN;
  if (!down_b_unit.allFinite() || !mag_b_uT.allFinite()) return false;
  const Eigen::Vector3f FWD_B(1.0f, 0.0f, 0.0f);

  Eigen::Vector3f d = down_b_unit;
  const float dn = d.norm();
  if (!std::isfinite(dn) || !(dn > 1e-6f)) return false;
  d /= dn;

  Eigen::Vector3f m = mag_b_uT;
  const float mn = m.norm();
  if (!std::isfinite(mn) || !(mn > 1e-6f)) return false;
  m /= mn;

  Eigen::Vector3f east_b = d.cross(m);
  const float en = east_b.norm();
  if (!std::isfinite(en) || !(en > 1e-6f)) return false;
  east_b /= en;

  Eigen::Vector3f north_b = east_b.cross(d);
  const float nn = north_b.norm();
  if (!std::isfinite(nn) || !(nn > 1e-6f)) return false;
  north_b /= nn;

  const float e = east_b.dot(FWD_B);
  const float n = north_b.dot(FWD_B);
  if (!(std::hypot(e, n) > 1e-6f)) return false;
  heading_deg_out = wrap360Deg(radToDeg(atan2f(e, n)));
  return std::isfinite(heading_deg_out);
}

// Passes at most one magnetometer sample per spacing interval.
//
// The AtomS3R magnetometer runs at 25 Hz (40 ms) while the INS loop polls at
// 200 Hz and sees the same reading several times.  A spacing a little under
// the sample period accepts each real sample about once.  A rejected
// candidate (mag not usable) re-arms the gate, so the next usable sample
// passes immediately.
class MagFreshGate {
 public:
  static constexpr uint32_t DEFAULT_SPACING_MS = 35u;

  explicit MagFreshGate(uint32_t spacing_ms = DEFAULT_SPACING_MS)
      : spacing_ms_(spacing_ms) {}

  void reset() { last_ms_ = 0; }

  bool update(bool mag_ok, uint32_t now_ms) {
    if (!mag_ok) {
      last_ms_ = 0;
      return false;
    }

    if (last_ms_ == 0) {
      last_ms_ = now_ms;
      return true;
    }

    if ((now_ms - last_ms_) < spacing_ms_) {
      return false;
    }

    last_ms_ = now_ms;
    return true;
  }

 private:
  uint32_t spacing_ms_;
  uint32_t last_ms_ = 0;
};

}  // namespace ins
}  // namespace ocean_imu
