#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  MagAutoTuner

  Startup helper for estimating a magnetic reference without leaking the
  arbitrary startup yaw gauge.

  Important convention:

    The caller gives a BODY->WORLD tilt quaternion. This tuner defensively
    removes any residual yaw from it before using it.

    It estimates the magnetic field in a yaw-free tilt frame:

        m_tilt = q_tilt_yaw_free * mag_body

    The averaged vector m_mean generally has a horizontal angle:

        yaw_gauge = atan2(m_mean.y, m_mean.x)

    The gauge-fixed magnetic reference returned to the MEKF is:

        B_world_ref = [hypot(m_mean.x, m_mean.y), 0, m_mean.z]

    The caller may apply one startup correction:

        q_new = Rz(-yaw_gauge) * q_old

    The tuner itself never uses the MEKF yaw as an input to the reference.
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
    float mag_norm_min = 1e-3f;

    // Real-device default: quick lock after the outer wrapper's mag delay and
    // gravity/tilt bootstrap are complete.
    int   min_samples = 295;
    float min_window_sec = 0.0f;
    float max_window_sec = 0.0f;
    float sample_dt_sec = 1.0f / 200.0f;

    // Optional sanity check: reject large magnetic norm outliers after a mean exists.
    float max_sample_norm_ratio_from_mean = 0.35f;

    // Require non-degenerate horizontal magnetic component.
    float min_horizontal_fraction = 0.05f;

    // Keep weighting OFF by default. Accel-norm weighting in waves can phase-bias
    // the mag reference. These fields remain for compatibility.
    bool  enable_quality_weighting = false;
    float min_sample_weight = 0.03f;
    float min_effective_weight = 0.0f;
    float gravity_ref = 9.80665f;
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
    mag_tilt_sum_.setZero();
    mag_tilt_norm_sum_ = 0.0f;

    accepted_count_ = 0;
    rejected_count_ = 0;

    accepted_window_sec_ = 0.0f;
    weight_sum_ = 0.0f;

    ready_ = false;

    mag_tilt_mean_.setZero();
    mag_world_ref_.setZero();

    yaw_gauge_rad_ = NAN;

    last_sample_weight_ = 0.0f;
    last_mag_tilt_sample_.setZero();
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
    last_mag_tilt_sample_.setZero();

    if (!mag_body_ned.allFinite()) {
      ++rejected_count_;
      return false;
    }

    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) {
      ++rejected_count_;
      return false;
    }

    const Eigen::Quaternionf q_tilt_bw =
        yawRemovedBodyToWorldQuat_(q_tilt_bw_in);

    if (!q_tilt_bw.coeffs().allFinite()) {
      ++rejected_count_;
      return false;
    }

    // This is the key point: rotate mag using a yaw-free tilt frame only.
    // No MEKF yaw is allowed to enter the learned reference.
    const Eigen::Vector3f mag_tilt_i = q_tilt_bw * mag_body_ned;

    if (!mag_tilt_i.allFinite()) {
      ++rejected_count_;
      return false;
    }

    const float mag_tilt_i_n = mag_tilt_i.norm();
    if (!(mag_tilt_i_n > cfg_.mag_norm_min) || !std::isfinite(mag_tilt_i_n)) {
      ++rejected_count_;
      return false;
    }

    if (accepted_count_ > 0 && weight_sum_ > 1e-6f) {
      const float mean_n = mag_tilt_norm_sum_ / weight_sum_;

      if (mean_n > cfg_.mag_norm_min && std::isfinite(mean_n)) {
        const float rel = std::fabs(mag_tilt_i_n - mean_n) / mean_n;

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
    }

    if (!std::isfinite(w)) {
      ++rejected_count_;
      return false;
    }

    w = std::min(std::max(w, 0.0f), 1.0f);

    if (cfg_.enable_quality_weighting && w < cfg_.min_sample_weight) {
      ++rejected_count_;
      return false;
    }

    const float dt_use =
        (std::isfinite(dt) && dt > 0.0f)
            ? dt
            : std::max(cfg_.sample_dt_sec, 0.0f);

    mag_tilt_sum_ += w * mag_tilt_i;
    mag_tilt_norm_sum_ += w * mag_tilt_i_n;
    weight_sum_ += w;
    accepted_window_sec_ += dt_use;
    ++accepted_count_;

    last_sample_weight_ = w;
    last_mag_tilt_sample_ = mag_tilt_i;

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
    return getLastMagTiltSample(mag_world_sample);
  }

  bool getLastMagTiltSample(Eigen::Vector3f& mag_tilt_sample) const {
    if (!last_mag_tilt_sample_.allFinite()) return false;
    mag_tilt_sample = last_mag_tilt_sample_;
    return true;
  }

  bool getMagWorldRef(Eigen::Vector3f& mag_world_ref) const {
    if (!ready_) return false;
    mag_world_ref = mag_world_ref_;
    return mag_world_ref.allFinite();
  }

  bool getMagWorldMean(Eigen::Vector3f& mag_world_mean) const {
    return getMagTiltMean(mag_world_mean);
  }

  bool getMagTiltMean(Eigen::Vector3f& mag_tilt_mean) const {
    if (accepted_count_ <= 0 || !(weight_sum_ > 1e-6f)) return false;

    mag_tilt_mean = mag_tilt_sum_ / weight_sum_;
    return mag_tilt_mean.allFinite();
  }

  float getYawGaugeCorrectionRad() const {
    if (ready_ && std::isfinite(yaw_gauge_rad_)) {
      return yaw_gauge_rad_;
    }

    Eigen::Vector3f m;
    if (!getMagTiltMean(m)) return NAN;

    const float h2 = m.x() * m.x() + m.y() * m.y();
    if (!(h2 > 1e-12f) || !std::isfinite(h2)) {
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

    const float h = std::sqrt(ref.x() * ref.x() + ref.y() * ref.y());
    if (!(h > 1e-12f) || !std::isfinite(h)) return NAN;

    return std::atan2(ref.z(), h);
  }

  float getLearnedDipDeg() const {
    const float r = getLearnedDipRad();
    return std::isfinite(r) ? r * 57.29577951308232f : NAN;
  }

private:
  static float yawFromBodyToWorldQuatRad_(const Eigen::Quaternionf& q_bw_in) {
    if (!q_bw_in.coeffs().allFinite()) return NAN;

    Eigen::Quaternionf q_bw = q_bw_in;
    const float qn = q_bw.norm();

    if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
      return NAN;
    }

    q_bw.normalize();

    const Eigen::Matrix3f R = q_bw.toRotationMatrix();

    const float c = R(0, 0);
    const float s = R(1, 0);
    const float h2 = c * c + s * s;

    if (!(h2 > 1.0e-8f) ||
        !std::isfinite(h2) ||
        !std::isfinite(c) ||
        !std::isfinite(s))
    {
      return NAN;
    }

    return std::atan2(s, c);
  }

  static Eigen::Quaternionf yawRemovedBodyToWorldQuat_(
      const Eigen::Quaternionf& q_bw_in)
  {
    if (!q_bw_in.coeffs().allFinite()) {
      return Eigen::Quaternionf::Identity();
    }

    Eigen::Quaternionf q_bw = q_bw_in;
    const float qn = q_bw.norm();

    if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
      return Eigen::Quaternionf::Identity();
    }

    q_bw.normalize();

    const float yaw = yawFromBodyToWorldQuatRad_(q_bw);

    if (!std::isfinite(yaw)) {
      return Eigen::Quaternionf::Identity();
    }

    const Eigen::Quaternionf q_yaw_inv(
        Eigen::AngleAxisf(-yaw, Eigen::Vector3f::UnitZ()));

    Eigen::Quaternionf q_tilt = q_yaw_inv * q_bw;
    q_tilt.normalize();

    if (!q_tilt.coeffs().allFinite()) {
      return Eigen::Quaternionf::Identity();
    }

    return q_tilt;
  }

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

    const Eigen::Vector3f m = mag_tilt_sum_ / weight_sum_;

    if (!m.allFinite()) {
      ready_ = false;
      return false;
    }

    const float mean_norm = m.norm();
    if (!(mean_norm > cfg_.mag_norm_min) || !std::isfinite(mean_norm)) {
      ready_ = false;
      return false;
    }

    const float horiz = std::sqrt(m.x() * m.x() + m.y() * m.y());
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

    const float yaw_gauge = std::atan2(m.y(), m.x());
    if (!std::isfinite(yaw_gauge)) {
      ready_ = false;
      return false;
    }

    mag_tilt_mean_ = m;
    yaw_gauge_rad_ = yaw_gauge;

    // Gauge-fixed magnetic-NED-like reference for the MEKF.
    mag_world_ref_ = Eigen::Vector3f(horiz, 0.0f, m.z());

    ready_ =
        mag_world_ref_.allFinite() &&
        (mag_world_ref_.norm() > cfg_.mag_norm_min);

    return ready_;
  }

private:
  Config cfg_{};

  Eigen::Vector3f mag_tilt_sum_ = Eigen::Vector3f::Zero();
  float mag_tilt_norm_sum_ = 0.0f;

  int accepted_count_ = 0;
  int rejected_count_ = 0;

  float accepted_window_sec_ = 0.0f;
  float weight_sum_ = 0.0f;

  bool ready_ = false;

  Eigen::Vector3f mag_tilt_mean_ = Eigen::Vector3f::Zero();
  Eigen::Vector3f mag_world_ref_ = Eigen::Vector3f::Zero();

  float yaw_gauge_rad_ = NAN;

  float last_sample_weight_ = 0.0f;
  Eigen::Vector3f last_mag_tilt_sample_ = Eigen::Vector3f::Zero();
};
