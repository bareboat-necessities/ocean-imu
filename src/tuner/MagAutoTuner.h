#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  MagAutoTuner

  Startup helper for estimating the magnetic reference.

  Correct frame rule:

    Preferred real-device chain:

      1. Collect magnetometer samples in the SAME MEKF world frame that will
         later consume the mag reference:

             m_world_raw = q_mekf_body_to_world * mag_body

      2. Average m_world_raw.

      3. Extract the arbitrary yaw gauge:

             yaw_gauge = atan2(mean.y, mean.x)

      4. Store a gauge-fixed magnetic reference:

             B_world_ref = [horizontal_magnitude, 0, vertical_component]

      5. Caller rotates the MEKF quaternion once by:

             q_new = Rz(-yaw_gauge) * q_old

    This does NOT leak arbitrary startup yaw. It explicitly measures and removes
    that yaw gauge in the same frame used by the MEKF mag update.

  Important:

    Do not mix a separate gravity-only leveling frame with the MEKF world frame
    unless you also carefully rotate every state into that same frame. For this
    filter, the safe startup path is addSampleWithWorldQuatDt().
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

    // For real boats/waves, a longer averaging window is safer than the old
    // 295-sample startup latch. The wrapper may still override this.
    int   min_samples    = 1500;          // ~7.5 s at 200 Hz
    float min_window_sec = 10.0f;
    float max_window_sec = 0.0f;          // 0 = no forced timeout
    float sample_dt_sec  = 1.0f / 200.0f;

    float gravity_ref = 9.80665f;

    // Reject samples whose norm becomes wildly inconsistent with the running
    // accepted-sample mean norm.
    float max_sample_norm_ratio_from_mean = 0.35f;

    // Require non-degenerate horizontal magnetic field.
    float min_horizontal_fraction = 0.05f;

    // Keep off by default for boats/waves. Accel/gyro quality weighting can
    // phase-select wave motion and create a constant yaw offset.
    bool  enable_quality_weighting = false;
    float min_sample_weight        = 0.03f;
    float min_effective_weight     = 0.0f;

    float acc_norm_rel_soft = 0.22f;
    float gyro_soft_dps     = 45.0f;
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
    mag_world_mean_.setZero();
    mag_world_ref_.setZero();

    last_sample_weight_ = 0.0f;
    last_mag_world_sample_.setZero();
  }

  // Preferred real-device API.
  //
  // q_bw_in must be the same BODY->WORLD quaternion used by the MEKF.
  // The averaged mag vector may contain arbitrary startup yaw. That is expected.
  // getYawGaugeCorrectionRad() measures that gauge so the caller can remove it
  // once by q_new = Rz(-yaw_gauge) * q_old.
  bool addSampleWithWorldQuatDt(float dt,
                                const Eigen::Quaternionf& q_bw_in,
                                const Eigen::Vector3f& acc_body_ned,
                                const Eigen::Vector3f& gyro_body_ned,
                                const Eigen::Vector3f& mag_body_ned)
  {
    if (ready_) return true;

    last_sample_weight_ = 0.0f;
    last_mag_world_sample_.setZero();

    if (!q_bw_in.coeffs().allFinite() || !mag_body_ned.allFinite()) {
      ++rejected_count_;
      return false;
    }

    Eigen::Quaternionf q_bw = q_bw_in;
    const float qn = q_bw.norm();

    if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
      ++rejected_count_;
      return false;
    }

    q_bw.normalize();

    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) {
      ++rejected_count_;
      return false;
    }

    const Eigen::Vector3f mag_world_i = q_bw * mag_body_ned;

    return addWorldSample_(dt,
                           mag_world_i,
                           acc_body_ned,
                           gyro_body_ned);
  }

  // Backward-compatible API.
  //
  // Only use this if q_tilt_bw is truly in the same intended startup gauge.
  // For the OU wrappers, prefer addSampleWithWorldQuatDt().
  bool addSampleWithTiltQuat(const Eigen::Quaternionf& q_tilt_bw,
                             const Eigen::Vector3f& acc_body_ned,
                             const Eigen::Vector3f& gyro_body_ned,
                             const Eigen::Vector3f& mag_body_ned)
  {
    return addSampleWithTiltQuatDt(cfg_.sample_dt_sec,
                                   q_tilt_bw,
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

    if (!q_tilt_bw_in.coeffs().allFinite() || !mag_body_ned.allFinite()) {
      ++rejected_count_;
      return false;
    }

    Eigen::Quaternionf q = q_tilt_bw_in;
    const float qn = q.norm();

    if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
      ++rejected_count_;
      return false;
    }

    q.normalize();

    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) {
      ++rejected_count_;
      return false;
    }

    const Eigen::Vector3f mag_world_i = q * mag_body_ned;

    return addWorldSample_(dt,
                           mag_world_i,
                           acc_body_ned,
                           gyro_body_ned);
  }

  // Legacy helper. Kept for compatibility, but the OU3 wrapper should not use
  // this for startup mag reference because it creates a different level frame
  // from the MEKF world frame.
  bool addSampleWithGravityDirDt(float dt,
                                 const Eigen::Vector3f& down_body_in,
                                 const Eigen::Vector3f& acc_body_ned,
                                 const Eigen::Vector3f& gyro_body_ned,
                                 const Eigen::Vector3f& mag_body_ned)
  {
    if (ready_) return true;

    if (!down_body_in.allFinite() || !mag_body_ned.allFinite()) {
      ++rejected_count_;
      return false;
    }

    Eigen::Vector3f down_body = down_body_in;
    const float dn = down_body.norm();

    if (!(dn > 1.0e-6f) || !std::isfinite(dn)) {
      ++rejected_count_;
      return false;
    }

    down_body /= dn;

    const float mag_n = mag_body_ned.norm();
    if (!(mag_n > cfg_.mag_norm_min) || !std::isfinite(mag_n)) {
      ++rejected_count_;
      return false;
    }

    const Eigen::Quaternionf q_level_bw =
        levelQuatFromDownBody_(down_body);

    if (!q_level_bw.coeffs().allFinite()) {
      ++rejected_count_;
      return false;
    }

    const Eigen::Vector3f mag_level_i = q_level_bw * mag_body_ned;

    return addWorldSample_(dt,
                           mag_level_i,
                           acc_body_ned,
                           gyro_body_ned);
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

  // Raw weighted mean before gauge fixing.
  bool getMagWorldMean(Eigen::Vector3f& mag_world_mean) const {
    if (accepted_count_ <= 0 || !(weight_sum_ > 1.0e-6f)) return false;

    mag_world_mean = mag_world_sum_ / weight_sum_;
    return mag_world_mean.allFinite();
  }

  // Horizontal yaw gauge of the averaged magnetic vector in the accumulation
  // frame. The caller removes it by pre-multiplying the MEKF quaternion:
  //
  //   q_new = Rz(-yaw_gauge) * q_old
  float getYawGaugeCorrectionRad() const {
    Eigen::Vector3f m;
    if (!getMagWorldMean(m)) return NAN;

    const float h2 = m.x() * m.x() + m.y() * m.y();
    if (!(h2 > 1.0e-12f) || !std::isfinite(h2)) {
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
    if (!(h > 1.0e-12f) || !std::isfinite(h)) return NAN;

    return std::atan2(ref.z(), h);
  }

  float getLearnedDipDeg() const {
    const float r = getLearnedDipRad();
    return std::isfinite(r) ? r * 57.29577951308232f : NAN;
  }

private:
  static float clamp01_(float x) {
    return std::min(std::max(x, 0.0f), 1.0f);
  }

  static Eigen::Quaternionf levelQuatFromDownBody_(
      const Eigen::Vector3f& down_body_in)
  {
    Eigen::Vector3f a = down_body_in;
    const float an = a.norm();

    if (!(an > 1.0e-6f) || !std::isfinite(an) || !a.allFinite()) {
      return Eigen::Quaternionf::Identity();
    }

    a /= an;

    const Eigen::Vector3f b(0.0f, 0.0f, 1.0f);

    float d = a.dot(b);
    d = std::min(std::max(d, -1.0f), 1.0f);

    if (d > 1.0f - 1.0e-6f) {
      return Eigen::Quaternionf::Identity();
    }

    if (d < -1.0f + 1.0e-6f) {
      return Eigen::Quaternionf(
          Eigen::AngleAxisf(float(M_PI),
                            Eigen::Vector3f(1.0f, 0.0f, 0.0f)));
    }

    Eigen::Vector3f axis = a.cross(b);
    const float axis_n = axis.norm();

    if (!(axis_n > 1.0e-6f) || !axis.allFinite()) {
      return Eigen::Quaternionf::Identity();
    }

    axis /= axis_n;

    Eigen::Quaternionf q(
        Eigen::AngleAxisf(std::acos(d), axis));

    q.normalize();
    return q;
  }

  bool addWorldSample_(float dt,
                       const Eigen::Vector3f& mag_world_i,
                       const Eigen::Vector3f& acc_body_ned,
                       const Eigen::Vector3f& gyro_body_ned)
  {
    last_sample_weight_ = 0.0f;
    last_mag_world_sample_.setZero();

    if (!mag_world_i.allFinite()) {
      ++rejected_count_;
      return false;
    }

    const float mag_world_n = mag_world_i.norm();

    if (!(mag_world_n > cfg_.mag_norm_min) || !std::isfinite(mag_world_n)) {
      ++rejected_count_;
      return false;
    }

    if (accepted_count_ > 0 && weight_sum_ > 1.0e-6f) {
      const float mean_n = mag_world_norm_sum_ / weight_sum_;

      if (mean_n > cfg_.mag_norm_min && std::isfinite(mean_n)) {
        const float rel = std::fabs(mag_world_n - mean_n) / mean_n;

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

    w = clamp01_(w);

    if (w < cfg_.min_sample_weight) {
      ++rejected_count_;
      return false;
    }

    const float dt_use =
        (std::isfinite(dt) && dt > 0.0f)
            ? dt
            : std::max(cfg_.sample_dt_sec, 0.0f);

    mag_world_sum_ += w * mag_world_i;
    mag_world_norm_sum_ += w * mag_world_n;

    weight_sum_ += w;
    accepted_window_sec_ += dt_use;
    ++accepted_count_;

    last_sample_weight_ = w;
    last_mag_world_sample_ = mag_world_i;

    return tryFinalize_();
  }

  float sampleWeight_(const Eigen::Vector3f& acc_body_ned,
                      const Eigen::Vector3f& gyro_body_ned) const
  {
    if (!acc_body_ned.allFinite() || !gyro_body_ned.allFinite()) {
      return 0.0f;
    }

    const float g = std::max(cfg_.gravity_ref, 1.0e-6f);

    const float an = acc_body_ned.norm();
    if (!(an > 1.0e-6f) || !std::isfinite(an)) {
      return 0.0f;
    }

    const float acc_rel_err = std::fabs(an - g) / g;
    const float acc_soft = std::max(cfg_.acc_norm_rel_soft, 1.0e-3f);

    float w_acc = 1.0f - acc_rel_err / acc_soft;
    w_acc = clamp01_(w_acc);

    const float gyro_dps = gyro_body_ned.norm() * 57.29577951308232f;
    if (!std::isfinite(gyro_dps)) {
      return 0.0f;
    }

    const float gyro_soft = std::max(cfg_.gyro_soft_dps, 1.0e-3f);

    float w_gyro = 1.0f - gyro_dps / gyro_soft;
    w_gyro = clamp01_(w_gyro);

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

    if (!(weight_sum_ > 1.0e-6f) || !std::isfinite(weight_sum_)) {
      ready_ = false;
      return false;
    }

    const Eigen::Vector3f mean = mag_world_sum_ / weight_sum_;

    if (!mean.allFinite()) {
      ready_ = false;
      return false;
    }

    const float mean_norm = mean.norm();

    if (!(mean_norm > cfg_.mag_norm_min) || !std::isfinite(mean_norm)) {
      ready_ = false;
      return false;
    }

    const float horiz =
        std::sqrt(mean.x() * mean.x() + mean.y() * mean.y());

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

    mag_world_mean_ = mean;

    // Gauge-fixed reference. This is the only reference the MEKF should use.
    mag_world_ref_ = Eigen::Vector3f(horiz, 0.0f, mean.z());

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
