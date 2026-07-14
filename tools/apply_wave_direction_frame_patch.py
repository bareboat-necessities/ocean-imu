#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_once(path: str, old: str, new: str) -> None:
    p = ROOT / path
    text = p.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one occurrence, found {count}: {old[:120]!r}")
    p.write_text(text.replace(old, new, 1), encoding="utf-8")


frame_header = r'''#pragma once

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
'''

(ROOT / "src/wave_dir/WaveDirectionFrame.h").write_text(frame_header, encoding="utf-8")

for path in [
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
]:
    replace_once(
        path,
        '#include "wave_dir/WaveDirectionDetector.h"\n',
        '#include "wave_dir/WaveDirectionDetector.h"\n'
        '#include "wave_dir/WaveDirectionFrame.h"\n')

    old = '''        // Stage 1 estimates the apparent propagation plane as an unsigned axis
        // relative to boat +X.  Stage 2 resolves propagation sense along that
        // same axis from horizontal/vertical orbital phase.
        dir_filter_.update(a_x_body, a_y_body, omega, dt);
        const Eigen::Vector2f propagation_axis_body = dir_filter_.getAxis();
        dir_sign_state_ = dir_sign_.update(
            a_x_body, a_y_body, a_body_z_up_proxy_,
            propagation_axis_body.x(), propagation_axis_body.y(),
            dt, dir_filter_.getLastStableConfidence());
'''
    new = '''        // Resolve direction in a leveled frame aligned with boat heading.
        // This removes roll/pitch mixing while preserving 0 deg = bow and
        // positive angles toward starboard.  The existing body-Z proxy remains
        // the tuner/tracker input; direction uses coherent leveled components.
        const auto direction_accel = wave_direction::heading_frame_acceleration<float>(
            mekf_->quaternion_boat(), acc, g_std);

        // Stage 1 estimates the apparent propagation plane as an unsigned axis
        // relative to boat heading.  Stage 2 resolves propagation sense along
        // that same axis from horizontal/vertical orbital phase.
        dir_filter_.update(direction_accel.forward_ms2,
                           direction_accel.starboard_ms2,
                           omega, dt);
        const Eigen::Vector2f propagation_axis_boat = dir_filter_.getAxis();
        dir_sign_state_ = dir_sign_.update(
            direction_accel.forward_ms2,
            direction_accel.starboard_ms2,
            direction_accel.up_ms2,
            propagation_axis_boat.x(), propagation_axis_boat.y(),
            dt, dir_filter_.getLastStableConfidence());
'''
    replace_once(path, old, new)

# Extend the focused tests with a real quaternion/body-specific-force frame test.
test_path = "tests/wave_dir/wave-direction-test.cpp"
replace_once(
    test_path,
    '#include "wave_dir/WaveDirectionDetector.h"\n#include "wave_dir/WaveEncounter.h"\n',
    '#include "wave_dir/WaveDirectionDetector.h"\n'
    '#include "wave_dir/WaveDirectionFrame.h"\n'
    '#include "wave_dir/WaveEncounter.h"\n')

frame_test = r'''
void test_heading_frame_roll_pitch_yaw_invariance() {
  constexpr float g = 9.80665f;
  constexpr float deg = kPi / 180.0f;

  const std::vector<Eigen::Vector3f> attitudes_deg{
      {0.0f, 0.0f, 0.0f},
      {25.0f, -18.0f, 37.0f},
      {-42.0f, 31.0f, -123.0f},
      {65.0f, -35.0f, 178.0f},
  };
  const std::vector<Eigen::Vector3f> heading_accels{
      {0.55f, -0.21f, 0.37f},
      {-0.44f, 0.62f, -0.19f},
      {0.08f, 0.91f, 0.12f},
  };

  for (const auto& euler : attitudes_deg) {
    const Eigen::AngleAxisf yaw(euler.z() * deg, Eigen::Vector3f::UnitZ());
    const Eigen::AngleAxisf pitch(euler.y() * deg, Eigen::Vector3f::UnitY());
    const Eigen::AngleAxisf roll(euler.x() * deg, Eigen::Vector3f::UnitX());
    const Eigen::Quaternionf q_body_to_world = yaw * pitch * roll;

    const Eigen::Vector3f bow_world = q_body_to_world * Eigen::Vector3f::UnitX();
    const Eigen::Vector2f heading =
        Eigen::Vector2f(bow_world.x(), bow_world.y()).normalized();
    const Eigen::Vector2f starboard(-heading.y(), heading.x());

    for (const auto& expected : heading_accels) {
      const Eigen::Vector3f acceleration_world(
          heading.x() * expected.x() + starboard.x() * expected.y(),
          heading.y() * expected.x() + starboard.y() * expected.y(),
          -expected.z());
      const Eigen::Vector3f specific_force_body =
          q_body_to_world.conjugate() *
          (acceleration_world - Eigen::Vector3f(0.0f, 0.0f, g));

      const auto recovered = wave_direction::heading_frame_acceleration<float>(
          q_body_to_world, specific_force_body, g);
      require(recovered.heading_valid, "heading-frame conversion rejected valid attitude");
      require(std::abs(recovered.forward_ms2 - expected.x()) < 2e-5f,
              "heading-frame forward acceleration mismatch");
      require(std::abs(recovered.starboard_ms2 - expected.y()) < 2e-5f,
              "heading-frame starboard acceleration mismatch");
      require(std::abs(recovered.up_ms2 - expected.z()) < 2e-5f,
              "heading-frame vertical acceleration mismatch");
    }
  }
}

'''
replace_once(
    test_path,
    'void test_encounter_forward_model_and_inverse() {\n',
    frame_test + 'void test_encounter_forward_model_and_inverse() {\n')
replace_once(
    test_path,
    '  test_combined_axis_and_sense();\n  test_encounter_forward_model_and_inverse();\n',
    '  test_combined_axis_and_sense();\n'
    '  test_heading_frame_roll_pitch_yaw_invariance();\n'
    '  test_encounter_forward_model_and_inverse();\n')

print("Applied leveled boat-heading direction frame and tests.")
