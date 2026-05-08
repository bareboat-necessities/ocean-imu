#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
 */

#include <cmath>
#include <algorithm>

#ifdef EIGEN_NON_ARDUINO
  #include <Eigen/Dense>
#else
  #include <ArduinoEigenDense.h>
#endif

// Common startup helper for estimating the magnetic reference for IMU-only mode.
//
// Important convention:
//
//   The filter does NOT know true north.
//   The filter learns a magnetic-NED frame.
//
//   +X = learned magnetic north
//   +Y = magnetic east
//   +Z = down
//
// With only gyro + accel + mag, absolute true yaw is not observable unless an
// external declination / heading reference is supplied. Therefore this tuner
// intentionally fixes the yaw gauge by forcing the learned world magnetic
// reference to lie in the X-Z plane:
//
//   B_world = [horizontal_magnitude, 0, vertical_component]
//
// Do not replace this with the raw 3D averaged vector unless the filter is also
// given an external true-world yaw reference.
class MagAutoTuner {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  struct Config {
    // Minimum accepted magnetometer norm.
    float mag_norm_min = 1e-3f;

    // Number of accepted tilt-compensated samples before declaring ready.
    int min_samples = 40;

    // Optional sanity checks. Defaults are permissive enough to preserve current
    // behavior while rejecting clearly invalid samples.
    float max_sample_norm_ratio_from_mean = 0.35f; // used only after mean exists
    float min_horizontal_fraction = 0.05f;         // horiz / norm of learned mean
  };

  MagAutoTuner() : cfg_(Config{}) { reset(); }
  explicit MagAutoTuner(const Config& cfg) : cfg_(cfg) { reset(); }

  void setConfig(const Config& cfg) {
    cfg_ = cfg;
    reset();
  }

  void reset() {
    mag_world_sum_.setZero();
    mag_world_norm_sum_ = 0.0f;
    accepted_count_ = 0;
    ready_ = false;
    mag_world_ref_.setZero();
  }

  bool addSampleWithTiltQuat(const Eigen::Quaternionf& q_tilt_bw_in,
                             const Eigen::Vector3f& /*acc_body_ned*/,
                             const Eigen::Vector3f& /*gyro_body_ned*/,
                             const Eigen::Vector3f& mag_body_ned)
  {
    if (ready_) return true;

    if (!mag_body_ned.allFinite()) return false;
    if (!q_tilt_bw_in.coeffs().allFinite()) return false;

    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) {
      return false;
    }

    Eigen::Quaternionf q_tilt_bw = q_tilt_bw_in;
    const float qn = q_tilt_bw.norm();
    if (!(qn > 1e-6f) || !std::isfinite(qn)) {
      return false;
    }
    q_tilt_bw.normalize();

    // Rotate sample into the current tilt-leveled world frame.
    //
    // q_tilt_bw is BODY->WORLD with yaw removed / gauge-free.
    // This removes roll/pitch from the mag sample but leaves yaw gauge arbitrary.
    const Eigen::Vector3f mag_world_i = q_tilt_bw * mag_body_ned;
    if (!mag_world_i.allFinite()) return false;

    const float mag_world_i_n = mag_world_i.norm();
    if (!(mag_world_i_n > cfg_.mag_norm_min) || !std::isfinite(mag_world_i_n)) {
      return false;
    }

    // Once we have a running mean, reject samples whose field magnitude is wildly
    // inconsistent. This catches bad spikes without trying to infer yaw.
    if (accepted_count_ > 0 && mag_world_norm_sum_ > 0.0f) {
      const float mean_n = mag_world_norm_sum_ / static_cast<float>(accepted_count_);
      if (mean_n > cfg_.mag_norm_min && std::isfinite(mean_n)) {
        const float rel = std::fabs(mag_world_i_n - mean_n) / mean_n;
        if (std::isfinite(cfg_.max_sample_norm_ratio_from_mean) &&
            cfg_.max_sample_norm_ratio_from_mean > 0.0f &&
            rel > cfg_.max_sample_norm_ratio_from_mean)
        {
          return false;
        }
      }
    }

    mag_world_sum_ += mag_world_i;
    mag_world_norm_sum_ += mag_world_i_n;
    ++accepted_count_;

    if (accepted_count_ < std::max(1, cfg_.min_samples)) {
      return false;
    }

    const Eigen::Vector3f mag_world_mean =
        mag_world_sum_ / static_cast<float>(accepted_count_);

    if (!mag_world_mean.allFinite()) {
      ready_ = false;
      return false;
    }

    const float mean_norm = mag_world_mean.norm();
    if (!(mean_norm > cfg_.mag_norm_min) || !std::isfinite(mean_norm)) {
      ready_ = false;
      return false;
    }

    const float horiz = std::sqrt(
        mag_world_mean.x() * mag_world_mean.x() +
        mag_world_mean.y() * mag_world_mean.y());

    if (!(horiz > cfg_.mag_norm_min) || !std::isfinite(horiz)) {
      ready_ = false;
      return false;
    }

    const float horiz_frac = horiz / mean_norm;
    if (std::isfinite(cfg_.min_horizontal_fraction) &&
        cfg_.min_horizontal_fraction > 0.0f &&
        horiz_frac < cfg_.min_horizontal_fraction)
    {
      ready_ = false;
      return false;
    }

    // IMU-only yaw gauge fix:
    //
    // We keep the measured dip/vertical component, but remove the arbitrary
    // startup yaw angle by defining learned magnetic north as +X.
    //
    // This is the correct observable reference for an IMU-only magnetic frame.
    mag_world_ref_ = Eigen::Vector3f(horiz, 0.0f, mag_world_mean.z());

    ready_ =
        mag_world_ref_.allFinite() &&
        (mag_world_ref_.norm() > cfg_.mag_norm_min);

    return ready_;
  }

  bool isReady() const { return ready_; }
  int acceptedCount() const { return accepted_count_; }

  bool getMagWorldRef(Eigen::Vector3f& mag_world_ref) const {
    if (!ready_) return false;
    mag_world_ref = mag_world_ref_;
    return mag_world_ref.allFinite();
  }

private:
  Config cfg_;

  // Sum of tilt-compensated magnetic vectors in the yaw-gauge-free world frame.
  Eigen::Vector3f mag_world_sum_ = Eigen::Vector3f::Zero();
  float mag_world_norm_sum_ = 0.0f;

  int accepted_count_ = 0;
  bool ready_ = false;

  // Learned magnetic-world reference used by the MEKF mag measurement model.
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();
};
