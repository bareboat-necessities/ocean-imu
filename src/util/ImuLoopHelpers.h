#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Small per-sample state machines shared by the AtomS3R marine INS sketches:
  sample-timestamp dt, sample-rate meter, stillness-gated gyro-bias average
  and the rate-of-turn low-pass.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <stdint.h>

#include "util/AngleUtils.h"

namespace ocean_imu {
namespace ins {

// Integration step from consecutive IMU sample timestamps (microseconds,
// wrap-safe).  The first sample, a non-positive or non-finite step and, when
// max_dt_s > 0, a step longer than max_dt_s all fall back to nominal_dt_s.
class SampleDtTracker {
 public:
  SampleDtTracker(float nominal_dt_s, float max_dt_s = 0.0f)
      : nominal_dt_s_(nominal_dt_s), max_dt_s_(max_dt_s) {}

  void reset() {
    have_last_ = false;
    last_us_ = 0;
  }

  float update(uint32_t sample_us) {
    if (!have_last_) {
      have_last_ = true;
      last_us_ = sample_us;
      return nominal_dt_s_;
    }

    const uint32_t dt_us = sample_us - last_us_;
    last_us_ = sample_us;

    const float dt_s = static_cast<float>(dt_us) * 1.0e-6f;
    if (!(dt_s > 0.0f) || !std::isfinite(dt_s)) return nominal_dt_s_;
    if (max_dt_s_ > 0.0f && dt_s > max_dt_s_) return nominal_dt_s_;

    return dt_s;
  }

 private:
  float    nominal_dt_s_;
  float    max_dt_s_;
  bool     have_last_ = false;
  uint32_t last_us_   = 0;
};

// IMU and magnetometer sample rates, recomputed over windows of at least 1 s.
class SampleRateMeter {
 public:
  void reset(uint32_t now_ms) {
    window_ms_ = now_ms;
    imu_count_ = 0;
    mag_count_ = 0;
    imu_hz_ = 0.0f;
    mag_hz_ = 0.0f;
  }

  void countImu() { ++imu_count_; }
  void countMag() { ++mag_count_; }

  void update(uint32_t now_ms) {
    if (window_ms_ == 0) window_ms_ = now_ms;
    const uint32_t elapsed_ms = now_ms - window_ms_;
    if (elapsed_ms >= 1000u) {
      const float scale = 1000.0f / static_cast<float>(elapsed_ms);
      imu_hz_ = static_cast<float>(imu_count_) * scale;
      mag_hz_ = static_cast<float>(mag_count_) * scale;
      imu_count_ = 0;
      mag_count_ = 0;
      window_ms_ = now_ms;
    }
  }

  float imuHz() const { return imu_hz_; }
  float magHz() const { return mag_hz_; }

 private:
  uint32_t window_ms_ = 0;
  uint32_t imu_count_ = 0;
  uint32_t mag_count_ = 0;
  float    imu_hz_ = 0.0f;
  float    mag_hz_ = 0.0f;
};

// Gyro bias as an exponential average of the calibrated rate while the device
// is still: |a| within g_tol_frac of g and |w| below still_gyro_rad_s.  The
// first still sample seeds the average.
class StillGyroBiasEma {
 public:
  StillGyroBiasEma(float tau_s, float g_tol_frac, float still_gyro_rad_s)
      : tau_s_(tau_s), g_tol_frac_(g_tol_frac), still_gyro_rad_s_(still_gyro_rad_s) {}

  void reset() {
    ok_ = false;
    learning_ = false;
    bias_.setZero();
  }

  // Returns whether this sample counted as still.
  bool update(const Eigen::Vector3f& w_cal, const Eigen::Vector3f& a_cal,
              float g_mps2, float dt_s) {
    learning_ =
        (fabsf(a_cal.norm() - g_mps2) < g_tol_frac_ * g_mps2) &&
        (w_cal.norm() < still_gyro_rad_s_);

    if (learning_) {
      const float alpha_b = 1.0f - expf(-dt_s / tau_s_);
      if (!ok_) {
        ok_   = true;
        bias_ = w_cal;
      } else {
        bias_ += alpha_b * (w_cal - bias_);
      }
    }
    return learning_;
  }

  // w_cal minus the bias once one has been learned.
  Eigen::Vector3f corrected(const Eigen::Vector3f& w_cal) const {
    Eigen::Vector3f w_use = w_cal;
    if (ok_) w_use -= bias_;
    return w_use;
  }

  bool ok() const { return ok_; }
  bool learning() const { return learning_; }
  const Eigen::Vector3f& bias() const { return bias_; }

 private:
  float tau_s_;
  float g_tol_frac_;
  float still_gyro_rad_s_;
  bool  ok_ = false;
  bool  learning_ = false;
  Eigen::Vector3f bias_ = Eigen::Vector3f::Zero();
};

// Rate of turn in deg/min from a world-down yaw rate, clamped to +/-720 and
// low-passed with time constant tau_s.  The first sample seeds the filter.
class RateOfTurnFilter {
 public:
  explicit RateOfTurnFilter(float tau_s = 0.1f) : tau_s_(tau_s) {}

  void reset() {
    inited_ = false;
    dpm_ = 0.0f;
  }

  float update(float yaw_rate_rad_s, float dt_s) {
    float dpm_meas = static_cast<float>(static_cast<double>(yaw_rate_rad_s) * kRadToDeg * 60.0);
    dpm_meas = clampf(dpm_meas, -720.0f, 720.0f);

    const float alpha = 1.0f - expf(-dt_s / tau_s_);
    if (!inited_) {
      inited_ = true;
      dpm_ = dpm_meas;
    } else {
      dpm_ += alpha * (dpm_meas - dpm_);
    }
    return dpm_;
  }

  float dpm() const { return dpm_; }

 private:
  float tau_s_;
  bool  inited_ = false;
  float dpm_ = 0.0f;
};

}  // namespace ins
}  // namespace ocean_imu
