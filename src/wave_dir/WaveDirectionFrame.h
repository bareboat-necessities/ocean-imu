#pragma once

/*
  Convert body-frame specific force into a leveled frame whose horizontal axes
  remain aligned with boat heading:

    +X: forward along the projected bow heading
    +Y: starboard
    +Z: up

  The input quaternion maps BODY to WORLD/NED.  The returned acceleration is
  inertial acceleration, not specific force; gravity is restored in NED before
  resolving the vector in the leveled heading frame.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>

namespace wave_direction {

template <typename Real>
struct HeadingFrameAcceleration {
  Real forward_ms2 = Real(0);
  Real starboard_ms2 = Real(0);
  Real up_ms2 = Real(0);
  bool heading_valid = false;
};

template <typename Real>
inline HeadingFrameAcceleration<Real> heading_frame_acceleration(
    const Eigen::Quaternion<Real>& q_body_to_world_ned,
    const Eigen::Matrix<Real, 3, 1>& specific_force_body,
    Real gravity_ms2) {
  HeadingFrameAcceleration<Real> result;

  Eigen::Quaternion<Real> q = q_body_to_world_ned;
  const Real q_norm = q.norm();
  if (!(q_norm > Real(1e-8)) || !std::isfinite(static_cast<double>(q_norm)) ||
      !specific_force_body.allFinite() ||
      !std::isfinite(static_cast<double>(gravity_ms2))) {
    result.forward_ms2 = specific_force_body.x();
    result.starboard_ms2 = specific_force_body.y();
    result.up_ms2 = -(specific_force_body.z() + gravity_ms2);
    return result;
  }
  q.normalize();

  const Eigen::Matrix<Real, 3, 1> gravity_world_ned(
      Real(0), Real(0), gravity_ms2);
  const Eigen::Matrix<Real, 3, 1> acceleration_world_ned =
      q * specific_force_body + gravity_world_ned;

  // Project the physical bow axis into the NED horizontal plane.  This keeps
  // the output aligned with heading while removing roll and pitch.
  const Eigen::Matrix<Real, 3, 1> bow_world =
      q * Eigen::Matrix<Real, 3, 1>(Real(1), Real(0), Real(0));
  const Real heading_norm = std::hypot(bow_world.x(), bow_world.y());
  if (!(heading_norm > Real(1e-6)) ||
      !std::isfinite(static_cast<double>(heading_norm))) {
    result.forward_ms2 = specific_force_body.x();
    result.starboard_ms2 = specific_force_body.y();
    result.up_ms2 = -(specific_force_body.z() + gravity_ms2);
    return result;
  }

  const Real north_forward = bow_world.x() / heading_norm;
  const Real east_forward = bow_world.y() / heading_norm;
  const Real north_starboard = -east_forward;
  const Real east_starboard = north_forward;

  result.forward_ms2 = north_forward * acceleration_world_ned.x()
                     + east_forward * acceleration_world_ned.y();
  result.starboard_ms2 = north_starboard * acceleration_world_ned.x()
                       + east_starboard * acceleration_world_ned.y();
  result.up_ms2 = -acceleration_world_ned.z();
  result.heading_valid = true;
  return result;
}

}  // namespace wave_direction
