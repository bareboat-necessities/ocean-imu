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

// Startup helper for estimating a world magnetic reference by averaging
// tilt-compensated magnetometer samples.
class MagAutoTuner {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  struct Config {
    float g = 9.80665f;
    float mag_norm_min = 1e-3f;
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
  
    Eigen::Quaternionf q_tilt_bw = q_tilt_bw_in;
    const float qn = q_tilt_bw.norm();
    if (!(qn > 1e-6f) || !std::isfinite(qn)) return false;
    q_tilt_bw.coeffs() /= qn;
  
    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) return false;
  
    const Eigen::Vector3f mag_world_i = q_tilt_bw * mag_body_ned;
    if (!mag_world_i.allFinite()) return false;
  
    const float mw_n = mag_world_i.norm();
    if (!(mw_n > cfg_.mag_norm_min) || !std::isfinite(mw_n)) return false;
  
    // Direction-only averaging is more robust to norm variation.
    mag_world_sum_ += mag_world_i / mw_n;
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
  
    mag_world_ref_ = Eigen::Vector3f(horiz, 0.0f, mag_world_mean.z());
    ready_ = mag_world_ref_.allFinite() &&
             (mag_world_ref_.norm() > cfg_.mag_norm_min);
    return ready_;
  }

private:
  Config cfg_;

  // Average tilt-compensated world-frame magnetic vectors.
  Eigen::Vector3f mag_world_sum_ = Eigen::Vector3f::Zero();

  int accepted_count_ = 0;
  bool ready_ = false;
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();
};
