#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  MagAutoTuner

  Startup helper for estimating the magnetic reference in IMU-only mode.

  Convention:

    The filter does NOT know true north unless an external reference is supplied.
    The filter learns a magnetic-NED-like world frame:

      +X = learned magnetic north
      +Y = magnetic east
      +Z = down

    With only gyro + accel + mag, absolute true yaw is not observable. Therefore
    this tuner intentionally fixes the yaw gauge by forcing the learned magnetic
    world reference to lie in the X-Z plane:

      B_world_ref = [horizontal_magnitude, 0, vertical_component]

    getYawGaugeCorrectionRad() returns the horizontal angle of the averaged
    tilt-compensated magnetic vector in the yaw-removed tilt frame.

    Important:

      This angle is NOT a residual correction to apply on top of the current
      quaternion. It is a yaw-gauge angle.

      The caller should set absolute yaw like this:

        yaw_abs = -getYawGaugeCorrectionRad()

      then rebuild:

        q_body_to_world = Rz(yaw_abs) * yaw_removed_tilt(q_body_to_world)

    Do not replace the gauge-fixed reference with the raw 3D averaged vector
    unless the filter also has an external true-world yaw reference.
*/

#include <cmath>
#include <algorithm>

#ifdef EIGEN_NON_ARDUINO
  #include <Eigen/Dense>
#else
  #include <ArduinoEigenDense.h>
#endif

class MagAutoTuner {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  struct Config {
    // Minimum accepted magnetometer norm.
    float mag_norm_min = 1e-3f;

    // Default is intentionally old behavior: 295 samples, no time window,
    // no quality weighting. This avoids changing mag-lock timing while we
    // debug yaw-gauge composition.
    int min_samples = 295;

    // Optional accepted-sample time window. Disabled by default.
    float min_window_sec = 0.0f;
    float max_window_sec = 0.0f;

    // Backward-compatible assumed sample period for callers that do not pass dt.
    float sample_dt_sec = 1.0f / 200.0f;

    // Gravity reference for optional accel-norm quality weighting.
    float gravity_ref = 9.80665f;

    // Optional sanity checks.
    float max_sample_norm_ratio_from_mean = 0.35f;
    float min_horizontal_fraction = 0.05f;

    // Disabled by default to preserve old behavior.
    bool enable_quality_weighting = false;

    // Used only when quality weighting is enabled.
    float min_sample_weight = 0.03f;

    // Disabled by default. If <= 0, ignored.
    float min_effective_weight = 0.0f;

    // Soft weighting parameters.
    float acc_norm_rel_soft = 0.22f;
    float gyro_soft_dps = 45.0f;
  };

  MagAutoTuner() : cfg_(Config{}) {
    reset();
  }

  explicit MagAutoTuner(const Config& cfg) : cfg_(cfg) {
    reset();
  }

  void setConfig(const Config& cfg) {
    cfg_ = cfg;
    reset();
  }

  const Config& config() const {
    return cfg_;
  }

  void reset() {
    mag_world_sum_.setZero();
    mag_world_norm_sum_ = 0.0f;

    accepted_count_ = 0;
    rejected_count_ = 0;

    accepted_window_sec_ = 0.0f;
    weight_sum_ = 0.0f;

    ready_ = false;
    mag_world_ref_.setZero();
    mag_world_mean_.setZero();

    last_sample_weight_ = 0.0f;
    last_mag_world_sample_.setZero();
  }

  bool addSampleWithTiltQuat(const Eigen::Quaternionf& q_tilt_bw_in,
                             const Eigen::Vector3f& acc_body_ned,
                             const Eigen::Vector3f& gyro_body_ned,
                             const Eigen::Vector3f& mag_body_ned)
  {
    return addSampleWithTiltQuatDt(
        cfg_.sample_dt_sec,
        q_tilt_bw_in,
        acc_body_ned,
        gyro_body_ned,
        mag_body_ned);
  }

  bool addSampleWithTiltQuatDt(float dt,
                               const Eigen::Quaternionf& q_tilt_bw_in,
                               const Eigen::Vector3f& acc_body_ned,
                               const Eigen::Vector3f& gyro_body_ned,
                               const Eigen::Vector3f& mag_body_ned)
  {
    if (ready_) return true;

    last_sample_weight_ = 0.0f;
    last_mag_world_sample_.setZero();

    if (!mag_body_ned.allFinite()) {
      ++rejected_count_;
      return false;
    }

    if (!q_tilt_bw_in.coeffs().allFinite()) {
      ++rejected_count_;
      return false;
    }

    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) {
      ++rejected_count_;
      return false;
    }

    Eigen::Quaternionf q_tilt_bw = q_tilt_bw_in;
    const float qn = q_tilt_bw.norm();
    if (!(qn > 1e-6f) || !std::isfinite(qn)) {
      ++rejected_count_;
      return false;
    }
    q_tilt_bw.normalize();

    // q_tilt_bw must be BODY->WORLD with yaw removed.
    const Eigen::Vector3f mag_world_i = q_tilt_bw * mag_body_ned;
    if (!mag_world_i.allFinite()) {
      ++rejected_count_;
      return false;
    }

    const float mag_world_i_n = mag_world_i.norm();
    if (!(mag_world_i_n > cfg_.mag_norm_min) || !std::isfinite(mag_world_i_n)) {
      ++rejected_count_;
      return false;
    }

    if (accepted_count_ > 0 && weight_sum_ > 1e-6f) {
      const float mean_n = mag_world_norm_sum_ / weight_sum_;
      if (mean_n > cfg_.mag_norm_min && std::isfinite(mean_n)) {
        const float rel = std::fabs(mag_world_i_n - mean_n) / mean_n;

        if (std::isfinite(cfg_.max_sample_norm_ratio_from_mean) &&
            cfg_.max_sample_norm_ratio_from_mean > 0.0f &&
            rel > cfg_.max_sample_norm_ratio_from_mean)
        {
          ++rejected_count_;
          return false;
        }
      }
    }

    float w = 1.0f;

    if (cfg_.enable_quality_weighting) {
      w = sampleWeight_(acc_body_ned, gyro_body_ned);

      if (!std::isfinite(w)) {
        ++rejected_count_;
        return false;
      }

      w = std::min(std::max(w, 0.0f), 1.0f);

      if (w < cfg_.min_sample_weight) {
        ++rejected_count_;
        return false;
      }
    }

    const float dt_use =
        (std::isfinite(dt) && dt > 0.0f)
            ? dt
            : std::max(cfg_.sample_dt_sec, 0.0f);

    mag_world_sum_ += w * mag_world_i;
    mag_world_norm_sum_ += w * mag_world_i_n;
    weight_sum_ += w;
    accepted_window_sec_ += dt_use;
    ++accepted_count_;

    last_sample_weight_ = w;
    last_mag_world_sample_ = mag_world_i;

    return tryFinalize_();
  }

  bool isReady() const {
    return ready_;
  }

  int acceptedCount() const {
    return accepted_count_;
  }

  int rejectedCount() const {
    return rejected_count_;
  }

  float acceptedWindowSec() const {
    return accepted_window_sec_;
  }

  float effectiveWeight() const {
    return weight_sum_;
  }

  float lastSampleWeight() const {
    return last_sample_weight_;
  }

  bool getLastMagWorldSample(Eigen::Vector3f& mag_world_sample) const {
    if (!last_mag_world_sample_.allFinite()) return false;
    mag_world_sample = last_mag_world_sample_;
    return true;
  }

  bool getMagWorldRef(Eigen::Vector3f& mag_world_ref) const {
    if (!ready_) return false;
    mag_world_ref = mag_world_ref_;
    return mag_world_ref.allFinite();
  }

  bool getMagWorldMean(Eigen::Vector3f& mag_world_mean) const {
    if (accepted_count_ <= 0 || !(weight_sum_ > 1e-6f)) return false;

    mag_world_mean = mag_world_sum_ / weight_sum_;
    return mag_world_mean.allFinite();
  }

  // Returns yaw-gauge angle of averaged mag vector in yaw-removed tilt frame.
  //
  // This is NOT a residual correction.
  //
  // Caller should use:
  //
  //   yaw_abs = -getYawGaugeCorrectionRad()
  //
  // then rebuild q as:
  //
  //   q = Rz(yaw_abs) * yaw_removed_tilt(q)
  float getYawGaugeCorrectionRad() const {
    Eigen::Vector3f m;
    if (!getMagWorldMean(m)) return NAN;

    const float horiz2 = m.x() * m.x() + m.y() * m.y();
    if (!(horiz2 > 1e-12f) || !std::isfinite(horiz2)) {
      return NAN;
    }

    return std::atan2(m.y(), m.x());
  }

  float getYawGaugeCorrectionDeg() const {
    const float r = getYawGaugeCorrectionRad();
    return std::isfinite(r) ? r * 57.29577951308232f : NAN;
  }

  float getLearnedDipRad() const {
    Eigen::Vector3f ref;
    if (!getMagWorldRef(ref)) return NAN;

    const float horiz = std::sqrt(ref.x() * ref.x() + ref.y() * ref.y());
    if (!(horiz > 1e-12f) || !std::isfinite(horiz)) return NAN;

    return std::atan2(ref.z(), horiz);
  }

  float getLearnedDipDeg() const {
    const float r = getLearnedDipRad();
    return std::isfinite(r) ? r * 57.29577951308232f : NAN;
  }

private:
  float sampleWeight_(const Eigen::Vector3f& acc_body_ned,
                      const Eigen::Vector3f& gyro_body_ned) const
  {
    if (!acc_body_ned.allFinite() || !gyro_body_ned.allFinite()) {
      return 0.0f;
    }

    const float g = std::max(cfg_.gravity_ref, 1e-6f);

    const float an = acc_body_ned.norm();
    if (!(an > 1e-6f) || !std::isfinite(an)) {
      return 0.0f;
    }

    const float acc_rel_err = std::fabs(an - g) / g;

    const float acc_soft = std::max(cfg_.acc_norm_rel_soft, 1e-3f);
    float w_acc = 1.0f - acc_rel_err / acc_soft;
    w_acc = std::min(std::max(w_acc, 0.0f), 1.0f);

    const float gyro_dps = gyro_body_ned.norm() * 57.29577951308232f;
    if (!std::isfinite(gyro_dps)) {
      return 0.0f;
    }

    const float gyro_soft = std::max(cfg_.gyro_soft_dps, 1e-3f);
    float w_gyro = 1.0f - gyro_dps / gyro_soft;
    w_gyro = std::min(std::max(w_gyro, 0.0f), 1.0f);

    return w_acc * w_gyro;
  }

  bool tryFinalize_() {
    if (accepted_count_ < std::max(1, cfg_.min_samples)) {
      return false;
    }

    const bool timed_out =
        std::isfinite(cfg_.max_window_sec) &&
        cfg_.max_window_sec > 0.0f &&
        accepted_window_sec_ >= cfg_.max_window_sec;

    if (std::isfinite(cfg_.min_window_sec) &&
        cfg_.min_window_sec > 0.0f &&
        accepted_window_sec_ < cfg_.min_window_sec &&
        !timed_out)
    {
      return false;
    }

    if (std::isfinite(cfg_.min_effective_weight) &&
        cfg_.min_effective_weight > 0.0f &&
        weight_sum_ < cfg_.min_effective_weight &&
        !timed_out)
    {
      return false;
    }

    if (!(weight_sum_ > 1e-6f) || !std::isfinite(weight_sum_)) {
      ready_ = false;
      return false;
    }

    const Eigen::Vector3f mag_world_mean = mag_world_sum_ / weight_sum_;

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

    mag_world_mean_ = mag_world_mean;

    // Gauge-fixed magnetic reference.
    mag_world_ref_ = Eigen::Vector3f(horiz, 0.0f, mag_world_mean.z());

    ready_ =
        mag_world_ref_.allFinite() &&
        (mag_world_ref_.norm() > cfg_.mag_norm_min);

    return ready_;
  }

private:
  Config cfg_;

  Eigen::Vector3f mag_world_sum_ = Eigen::Vector3f::Zero();
  float mag_world_norm_sum_ = 0.0f;

  int accepted_count_ = 0;
  int rejected_count_ = 0;

  float accepted_window_sec_ = 0.0f;
  float weight_sum_ = 0.0f;

  bool ready_ = false;

  Eigen::Vector3f mag_world_mean_ = Eigen::Vector3f::Zero();
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();

  float last_sample_weight_ = 0.0f;
  Eigen::Vector3f last_mag_world_sample_ = Eigen::Vector3f::Zero();
};
