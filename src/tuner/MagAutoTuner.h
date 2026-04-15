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
    float g = 9.80665f;
    float accel_band_frac = 0.15f;                      // ±15% around 1g
    float gyro_norm_max = 50.0f * float(M_PI) / 180.0f; // 50 deg/s
    float mag_norm_min = 1e-3f;
    int   min_samples = 40;
  };

  MagAutoTuner() : cfg_(Config{}) { reset(); }
  explicit MagAutoTuner(const Config& cfg) : cfg_(cfg) { reset(); }

  void setConfig(const Config& cfg) { cfg_ = cfg; reset(); }

  void reset() {
    acc_sum_.setZero();
    mag_sum_.setZero();
    accepted_count_ = 0;
    ready_ = false;
    mag_world_ref_.setZero();
  }

  bool addSample(const Eigen::Vector3f& acc_body_ned,
                 const Eigen::Vector3f& gyro_body_ned,
                 const Eigen::Vector3f& mag_body_ned)
  {
    if (ready_) return true;
    if (!isStableSample_(acc_body_ned, gyro_body_ned)) return false;
    if (!mag_body_ned.allFinite()) return false;

    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) return false;

    acc_sum_ += acc_body_ned;
    mag_sum_ += mag_body_ned;
    ++accepted_count_;

    if (accepted_count_ < cfg_.min_samples) return false;

    const Eigen::Vector3f acc_mean = acc_sum_ / static_cast<float>(accepted_count_);
    const Eigen::Vector3f mag_mean = mag_sum_ / static_cast<float>(accepted_count_);

    Eigen::Quaternionf q_tilt = tiltOnlyQuatFromAccel_(acc_mean);
    q_tilt.normalize();

    mag_world_ref_ = q_tilt * mag_mean;
    ready_ = mag_world_ref_.allFinite() && (mag_world_ref_.norm() > cfg_.mag_norm_min);
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
  bool isStableSample_(const Eigen::Vector3f& a,
                       const Eigen::Vector3f& gyro_body_ned) const
  {
    if (!a.allFinite() || !gyro_body_ned.allFinite()) return false;

    const float g = cfg_.g;
    const float acc_band = cfg_.accel_band_frac * g;
    const float acc_n = a.norm();
    const float gyro_n = gyro_body_ned.norm();

    return std::isfinite(acc_n) &&
           std::isfinite(gyro_n) &&
           (std::fabs(acc_n - g) <= acc_band) &&
           (gyro_n <= cfg_.gyro_norm_max);
  }

  static Eigen::Quaternionf tiltOnlyQuatFromAccel_(const Eigen::Vector3f& acc_body_ned) {
    const float an = acc_body_ned.norm();
    if (!(an > 1e-6f) || !acc_body_ned.allFinite()) {
      return Eigen::Quaternionf::Identity();
    }

    const Eigen::Vector3f body_down = -(acc_body_ned / an);
    const Eigen::Vector3f world_down(0.0f, 0.0f, 1.0f);

    const float d = std::max(-1.0f, std::min(1.0f, body_down.dot(world_down)));
    Eigen::Vector3f axis = body_down.cross(world_down);
    const float axis_n = axis.norm();

    if (axis_n < 1e-6f) {
      if (d > 0.0f) {
        return Eigen::Quaternionf::Identity();
      }
      Eigen::Vector3f ortho = (std::fabs(body_down.z()) < 0.9f)
        ? Eigen::Vector3f(0, 0, 1).cross(body_down)
        : Eigen::Vector3f(0, 1, 0).cross(body_down);
      ortho.normalize();
      return Eigen::Quaternionf(Eigen::AngleAxisf(float(M_PI), ortho));
    }

    axis /= axis_n;
    const float angle = std::acos(d);
    Eigen::Quaternionf q(Eigen::AngleAxisf(angle, axis));
    q.normalize();
    return q;
  }

private:
  Config cfg_;

  Eigen::Vector3f acc_sum_ = Eigen::Vector3f::Zero();
  Eigen::Vector3f mag_sum_ = Eigen::Vector3f::Zero();
  int accepted_count_ = 0;
  bool  ready_ = false;
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();
};
