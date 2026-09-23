#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Quaternion helpers shared by the AtomS3R marine INS sketches.

  Convention: q_bw rotates BODY vectors into WORLD (NED, +Z down), which is
  what every INS filter in this library reports as its attitude.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>

#include "util/AngleUtils.h"

namespace ocean_imu {
namespace ins {

// v' = q v q*, without building a rotation matrix.
inline Eigen::Vector3f quatRotate(const Eigen::Quaternionf& q, const Eigen::Vector3f& v) {
  const Eigen::Vector3f qv(q.x(), q.y(), q.z());
  const Eigen::Vector3f t = 2.0f * qv.cross(v);
  return v + q.w() * t + qv.cross(t);
}

// ZYX Euler angles of q_bw: roll and pitch in (-180, 180], heading in [0, 360).
// Returns false, leaving the outputs untouched, when q_bw is not normalizable.
inline bool rollPitchHeadingFromQuatBw(const Eigen::Quaternionf& q_bw,
                                       float& roll_deg_out,
                                       float& pitch_deg_out,
                                       float& heading_deg_out) {
  Eigen::Quaternionf q = q_bw;
  const float nq = q.norm();
  if (!(nq > 1e-6f) || !std::isfinite(nq)) return false;
  q.normalize();

  const float x = q.x();
  const float y = q.y();
  const float z = q.z();
  const float w = q.w();

  const float siny_cosp = 2.0f * (w * z + x * y);
  const float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
  const float yaw = atan2f(siny_cosp, cosy_cosp);

  float sinp = 2.0f * (w * y - z * x);
  sinp = clampf(sinp, -1.0f, 1.0f);
  const float pitch = asinf(sinp);

  const float sinr_cosp = 2.0f * (w * x + y * z);
  const float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
  const float roll = atan2f(sinr_cosp, cosr_cosp);

  roll_deg_out    = wrap180Deg(radToDeg(roll));
  pitch_deg_out   = wrap180Deg(radToDeg(pitch));
  heading_deg_out = wrap360Deg(radToDeg(yaw));
  return true;
}

}  // namespace ins
}  // namespace ocean_imu
