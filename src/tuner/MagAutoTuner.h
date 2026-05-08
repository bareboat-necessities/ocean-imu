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
    float mag_norm_min = 1e-3f;
    int   min_samples = 40;

    // Extra strictness, but deliberately loose enough not to block high-sea runs.
    // These are evaluated on tilt-compensated mag_world samples.
    float max_norm_rel_std  = 0.12f;
    float max_vert_rel_std  = 0.16f;
    float min_horiz_fraction = 0.15f;
  };

  MagAutoTuner() : cfg_(Config{}) { reset(); }
  explicit MagAutoTuner(const Config& cfg) : cfg_(cfg) { reset(); }

  void setConfig(const Config& cfg) {
    cfg_ = cfg;
    reset();
  }

  void reset() {
    mag_world_sum_.setZero();
    mag_world_sq_sum_.setZero();

    norm_sum_ = 0.0f;
    norm_sq_sum_ = 0.0f;
    z_sum_ = 0.0f;
    z_sq_sum_ = 0.0f;

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
    q_tilt_bw.normalize();

    // Rotate each accepted mag sample into world using CURRENT filter tilt only.
    const Eigen::Vector3f mag_world_i = q_tilt_bw * mag_body_ned;
    if (!mag_world_i.allFinite()) return false;

    const float norm_i = mag_world_i.norm();
    if (!(norm_i > cfg_.mag_norm_min) || !std::isfinite(norm_i)) {
      return false;
    }

    mag_world_sum_ += mag_world_i;
    mag_world_sq_sum_ += mag_world_i.cwiseProduct(mag_world_i);

    norm_sum_ += norm_i;
    norm_sq_sum_ += norm_i * norm_i;

    z_sum_ += mag_world_i.z();
    z_sq_sum_ += mag_world_i.z() * mag_world_i.z();

    ++accepted_count_;

    if (accepted_count_ < cfg_.min_samples) {
      return false;
    }

    const float n = static_cast<float>(accepted_count_);

    const Eigen::Vector3f mag_world_mean = mag_world_sum_ / n;

    const float norm_mean = norm_sum_ / n;
    const float norm_var = std::max(0.0f, norm_sq_sum_ / n - norm_mean * norm_mean);
    const float norm_std = std::sqrt(norm_var);

    const float z_mean = z_sum_ / n;
    const float z_var = std::max(0.0f, z_sq_sum_ / n - z_mean * z_mean);
    const float z_std = std::sqrt(z_var);

    if (!mag_world_mean.allFinite()) {
      ready_ = false;
      return false;
    }

    if (!(norm_mean > cfg_.mag_norm_min) ||
        !std::isfinite(norm_mean) ||
        !std::isfinite(norm_std) ||
        !std::isfinite(z_std))
    {
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

    // Do not accept a nearly vertical learned field; yaw would be weak.
    if (horiz < cfg_.min_horiz_fraction * norm_mean) {
      ready_ = false;
      return false;
    }

    const float norm_rel_std = norm_std / std::max(norm_mean, cfg_.mag_norm_min);
    if (norm_rel_std > cfg_.max_norm_rel_std) {
      ready_ = false;
      return false;
    }

    const float z_scale = std::max(std::fabs(z_mean), horiz);
    const float z_rel_std = z_std / std::max(z_scale, cfg_.mag_norm_min);
    if (z_rel_std > cfg_.max_vert_rel_std) {
      ready_ = false;
      return false;
    }

    // Original behavior:
    // Define learned magnetic-world frame so +X is magnetic north and +Z is down.
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
  Eigen::Vector3f mag_world_sq_sum_ = Eigen::Vector3f::Zero();

  float norm_sum_ = 0.0f;
  float norm_sq_sum_ = 0.0f;
  float z_sum_ = 0.0f;
  float z_sq_sum_ = 0.0f;

  int accepted_count_ = 0;
  bool ready_ = false;
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();
};
