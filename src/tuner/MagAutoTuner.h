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
//
// Important behavior:
//   - Uses tilt-only attitude, so yaw is intentionally unknown during learning.
//   - Therefore we must NOT average horizontal vector components directly.
//     If the boat yaws while learning, direct x/y vector averaging shrinks the
//     horizontal magnetic field and creates a bad reference.
//
// This version averages:
//   - horizontal magnitude sqrt(mx^2 + my^2)
//   - vertical component z
//
// Then it defines the learned magnetic-world reference as:
//   B_world_ref = [mean_horizontal, 0, mean_z]
//
// So the rest of the MEKF code can remain unchanged.
class MagAutoTuner {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  struct Config {
    float mag_norm_min = 1e-3f;
    int   min_samples = 40;

    // Reject learning if the learned yaw-invariant horizontal magnitude is
    // too unstable. Relative to mean horizontal field.
    float max_horiz_rel_std = 0.08f;

    // Reject learning if vertical/dip component is too unstable.
    // Relative to max(abs(mean_z), mean_horizontal).
    float max_z_rel_std = 0.08f;
  };

  MagAutoTuner() : cfg_(Config{}) { reset(); }
  explicit MagAutoTuner(const Config& cfg) : cfg_(cfg) { reset(); }

  void setConfig(const Config& cfg) {
    cfg_ = cfg;
    reset();
  }

  void reset() {
    horiz_sum_ = 0.0f;
    horiz_sq_sum_ = 0.0f;
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

    // Rotate mag into a tilt-compensated frame.
    // Yaw is intentionally not trusted here.
    const Eigen::Vector3f mag_world_i = q_tilt_bw * mag_body_ned;
    if (!mag_world_i.allFinite()) return false;

    const float horiz_i = std::sqrt(
        mag_world_i.x() * mag_world_i.x() +
        mag_world_i.y() * mag_world_i.y());

    const float z_i = mag_world_i.z();

    if (!(horiz_i > cfg_.mag_norm_min) ||
        !std::isfinite(horiz_i) ||
        !std::isfinite(z_i))
    {
      return false;
    }

    horiz_sum_    += horiz_i;
    horiz_sq_sum_ += horiz_i * horiz_i;

    z_sum_    += z_i;
    z_sq_sum_ += z_i * z_i;

    ++accepted_count_;

    if (accepted_count_ < cfg_.min_samples) {
      return false;
    }

    const float n = static_cast<float>(accepted_count_);

    const float horiz_mean = horiz_sum_ / n;
    const float z_mean     = z_sum_ / n;

    const float horiz_var =
        std::max(0.0f, horiz_sq_sum_ / n - horiz_mean * horiz_mean);

    const float z_var =
        std::max(0.0f, z_sq_sum_ / n - z_mean * z_mean);

    const float horiz_std = std::sqrt(horiz_var);
    const float z_std     = std::sqrt(z_var);

    if (!(horiz_mean > cfg_.mag_norm_min) || !std::isfinite(horiz_mean)) {
      ready_ = false;
      return false;
    }

    if (!std::isfinite(z_mean) ||
        !std::isfinite(horiz_std) ||
        !std::isfinite(z_std))
    {
      ready_ = false;
      return false;
    }

    const float horiz_rel_std = horiz_std / std::max(horiz_mean, cfg_.mag_norm_min);
    const float z_scale = std::max(std::fabs(z_mean), horiz_mean);
    const float z_rel_std = z_std / std::max(z_scale, cfg_.mag_norm_min);

    if (horiz_rel_std > cfg_.max_horiz_rel_std) {
      ready_ = false;
      return false;
    }

    if (z_rel_std > cfg_.max_z_rel_std) {
      ready_ = false;
      return false;
    }

    // Define learned world frame so that +X is magnetic north and +Z is down.
    // We keep measured dip, but yaw is removed by forcing Y=0.
    //
    // This remains compatible with:
    //   impl_.mekf().set_mag_world_ref(mag_world_ref_uT);
    mag_world_ref_ = Eigen::Vector3f(horiz_mean, 0.0f, z_mean);

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

  float horizontalMean() const {
    if (accepted_count_ <= 0) return NAN;
    return horiz_sum_ / static_cast<float>(accepted_count_);
  }

  float verticalMean() const {
    if (accepted_count_ <= 0) return NAN;
    return z_sum_ / static_cast<float>(accepted_count_);
  }

private:
  Config cfg_;

  // Yaw-invariant horizontal magnitude statistics.
  float horiz_sum_ = 0.0f;
  float horiz_sq_sum_ = 0.0f;

  // Vertical magnetic component statistics.
  float z_sum_ = 0.0f;
  float z_sq_sum_ = 0.0f;

  int accepted_count_ = 0;
  bool ready_ = false;
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();
};
