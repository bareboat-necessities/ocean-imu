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

// Common startup helper for estimating world magnetic reference from
// stable body-frame accel + gyro + mag measurements.
class MagAutoTuner {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  struct Config {
    float accel_band_frac = 0.30f;                       // ±30% around 1g
    float mag_norm_min = 1e-4f;
    int   min_samples = 20;
  };

  MagAutoTuner() : cfg_(Config{}) { reset(); }
  explicit MagAutoTuner(const Config& cfg) : cfg_(cfg) { reset(); }

  void setConfig(const Config& cfg) {
    cfg_ = cfg;
    reset();
  }

  void reset() {
    mag_world_sum_.setZero();
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
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) return false;

    Eigen::Quaternionf q_tilt_bw = q_tilt_bw_in;
    q_tilt_bw.normalize();

    // Rotate each accepted mag sample into world using CURRENT filter tilt only.
    const Eigen::Vector3f mag_world_i = q_tilt_bw * mag_body_ned;
    if (!mag_world_i.allFinite()) return false;

    mag_world_sum_ += mag_world_i;
    ++accepted_count_;

    if (accepted_count_ < cfg_.min_samples) return false;

    const Eigen::Vector3f mag_world_mean =
        mag_world_sum_ / static_cast<float>(accepted_count_);

    const float horiz = std::sqrt(
        mag_world_mean.x() * mag_world_mean.x() +
        mag_world_mean.y() * mag_world_mean.y());

    if (!(horiz > cfg_.mag_norm_min) || !std::isfinite(horiz)) {
        ready_ = false;
        return false;
    }

    // Define learned world frame so that +X is magnetic north and +Z is down.
    // Keep measured dip, but remove unknown yaw by forcing Y=0.
    mag_world_ref_ = Eigen::Vector3f(horiz, 0.0f, mag_world_mean.z());
    ready_ = mag_world_ref_.allFinite() &&
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

  // Average tilt-compensated world-frame magnetic vectors.
  Eigen::Vector3f mag_world_sum_ = Eigen::Vector3f::Zero();

  int accepted_count_ = 0;
  bool ready_ = false;
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();
};
