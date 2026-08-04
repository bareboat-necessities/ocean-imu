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

    int   min_samples    = 250;
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

    // Joint hard-iron estimation.  Off by default, and deliberately so.
    //
    // A body-fixed offset is not separable from the world reference at a single
    // attitude: rotating north and biasing the sensor produce the same reading.
    // Only the tilt a hull takes through a wave breaks that tie, and it breaks
    // it weakly.  When the accepted samples do not span enough attitude the
    // solve is skipped, but "enough" is a threshold, and a fit that clears it
    // on a poorly conditioned window can still be dominated by whatever
    // unmodelled sensor error is correlated with tilt -- soft iron above all,
    // which this model has no term for.  A wrong offset subtracted from every
    // later sample is worse than no offset at all, so a caller has to ask.
    //
    // Ask for it when the platform actually changes heading during startup, or
    // when the field error is known to be dominated by hard iron.  See
    // getHardIronBodyUT() for what the caller does with the result.
    bool  estimate_hard_iron = false;

    // Fisher information the solve must clear, as weight_sum * lambda_min of
    // the normal matrix below.  The offset's error scales as
    // sigma / sqrt(weight_sum * lambda_min), so this is the knob that says how
    // precise the fit has to be before it is allowed to matter.
    float min_hard_iron_information = 25.0f;

    // Reject an implausible fit as a fraction of the measured field norm.
    float max_hard_iron_fraction = 0.35f;
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

    rot_sum_.setZero();
    mag_body_sum_.setZero();

    hard_iron_body_.setZero();
    hard_iron_valid_ = false;
    hard_iron_information_ = 0.0f;
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
                           q_bw,
                           mag_body_ned,
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
                           q,
                           mag_body_ned,
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
                           q_level_bw,
                           mag_body_ned,
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

  // Body-frame hard-iron offset estimated alongside the reference, when
  // Config::estimate_hard_iron asked for it and the window supported it.
  //
  // The caller subtracts this from every later magnetometer sample.  The MEKF
  // carries no mag-bias state, so an offset left in the stream shows up
  // directly as heading: on a field with horizontal component B_h, an offset
  // component b_h perpendicular to north is atan(b_h / B_h) of yaw error.
  //
  // False means no usable estimate, and the caller must then subtract nothing.
  bool getHardIronBodyUT(Eigen::Vector3f& hard_iron_body_uT) const {
    if (!ready_ || !hard_iron_valid_) return false;
    hard_iron_body_uT = hard_iron_body_;
    return hard_iron_body_uT.allFinite();
  }

  bool hasHardIron() const {
    return ready_ && hard_iron_valid_;
  }

  // weight_sum * lambda_min at the last solve: how much the window actually
  // constrained the offset.  Zero when the hull held one attitude throughout.
  float hardIronInformation() const {
    return hard_iron_information_;
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
                       const Eigen::Quaternionf& q_bw,
                       const Eigen::Vector3f& mag_body_ned,
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

    // Sufficient statistics for the joint solve.  Every accepted sample obeys
    //   mag_world_i = B_world + A_i * b_body
    // with A_i the accumulation-frame BODY->WORLD rotation, and the normal
    // equations close over these two running sums plus the world sum above.
    // Accumulated unconditionally: they cost 12 floats and no branch, and
    // whether the solve runs is decided once, at finalize.
    rot_sum_.noalias() += w * q_bw.toRotationMatrix();
    mag_body_sum_ += w * mag_body_ned;

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

    // Split the body-fixed offset out of the world field before the reference
    // is gauge-fixed.  Left in, it is baked into the reference at one attitude
    // and read back as heading error at every other.
    solveHardIron_();

    const Eigen::Vector3f mean =
        hard_iron_valid_
            ? Eigen::Vector3f(
                  (mag_world_sum_ - rot_sum_ * hard_iron_body_) / weight_sum_)
            : Eigen::Vector3f(mag_world_sum_ / weight_sum_);

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

  // Least squares over the accepted window for
  //
  //     mag_world_i = B_world + A_i * b_body
  //
  // Writing S = sum(w), Sa = sum(w A_i), Sw = sum(w mag_world_i) and
  // Sm = sum(w mag_body_i), and using A_i^T A_i = I, the normal equations are
  //
  //     S B + Sa b   = Sw
  //     Sa^T B + S b = Sm
  //
  // Eliminating B leaves one 3x3 system in b.  With Abar = Sa/S that is
  //
  //     (I - Abar^T Abar) b = mbar - Abar^T wbar.
  //
  // The left matrix is symmetric positive semidefinite and exactly singular
  // when the hull holds one attitude for the whole window, which is the case
  // where a body-fixed offset genuinely cannot be told apart from the world
  // field.  Its smallest eigenvalue times the accumulated weight is therefore
  // the information the fit is allowed to rely on.
  //
  // Note what this model does not contain: scale, cross-axis and misalignment
  // error.  Those are also modulated by attitude, so on a window that barely
  // clears the information floor they can dominate the very signal the offset
  // is being read from.  That is why the estimate is opt-in.
  void solveHardIron_() {
    hard_iron_body_.setZero();
    hard_iron_valid_ = false;
    hard_iron_information_ = 0.0f;

    if (!cfg_.estimate_hard_iron) return;
    if (!(weight_sum_ > 1.0e-6f) || !std::isfinite(weight_sum_)) return;
    if (!rot_sum_.allFinite() || !mag_body_sum_.allFinite()) return;
    if (!mag_world_sum_.allFinite()) return;

    const Eigen::Matrix3f Abar = rot_sum_ / weight_sum_;
    const Eigen::Vector3f wbar = mag_world_sum_ / weight_sum_;
    const Eigen::Vector3f mbar = mag_body_sum_ / weight_sum_;

    Eigen::Matrix3f M =
        Eigen::Matrix3f::Identity() - Abar.transpose() * Abar;
    if (!M.allFinite()) return;
    M = 0.5f * (M + M.transpose());

    const Eigen::Vector3f rhs = mbar - Abar.transpose() * wbar;
    if (!rhs.allFinite()) return;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es;
    es.compute(M);
    if (es.info() != Eigen::Success) return;

    const Eigen::Vector3f evals = es.eigenvalues();
    if (!evals.allFinite()) return;

    hard_iron_information_ = evals.minCoeff() * weight_sum_;

    Eigen::Vector3f inv_evals;
    for (int i = 0; i < 3; ++i) {
      if (!(evals(i) * weight_sum_ > cfg_.min_hard_iron_information)) return;
      inv_evals(i) = 1.0f / evals(i);
    }

    const Eigen::Vector3f b =
        es.eigenvectors() *
        (inv_evals.asDiagonal() *
         (es.eigenvectors().transpose() * rhs));
    if (!b.allFinite()) return;

    const float field_norm = mag_world_norm_sum_ / weight_sum_;
    if (!(field_norm > cfg_.mag_norm_min) || !std::isfinite(field_norm)) return;
    if (!(b.norm() <= cfg_.max_hard_iron_fraction * field_norm)) return;

    hard_iron_body_ = b;
    hard_iron_valid_ = true;
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

  Eigen::Matrix3f rot_sum_ = Eigen::Matrix3f::Zero();
  Eigen::Vector3f mag_body_sum_ = Eigen::Vector3f::Zero();

  Eigen::Vector3f hard_iron_body_ = Eigen::Vector3f::Zero();
  bool  hard_iron_valid_ = false;
  float hard_iron_information_ = 0.0f;

  float last_sample_weight_ = 0.0f;
  Eigen::Vector3f last_mag_world_sample_ = Eigen::Vector3f::Zero();
};
