#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Accuracy checks of calibrated IMU data against physical references, for
  diagnostics on a device. Nothing here is fed back to a calibration or an
  estimator.

  Still poses. Samples are grouped in blocks (default 50 = 0.25 s at 200 Hz).
  A block is quiet when the standard deviation of every accelerometer axis
  is below accel_std and of every gyro axis below gyro_std, the mean rate
  is below gyro_mean, and its mean gravity direction is within direction
  of the previous quiet block's. The standard deviations are a few times
  the noise a resting sensor shows, so a noisy sensor still yields poses
  and its noise is reported; a hand moving the board exceeds them. A run of quiet blocks is a still pose. Its
  first block (settling) and last block (motion may already have started)
  are left out; a pose of at least min_blocks merged blocks is reported, and
  a long one every max_blocks blocks.

  For each pose:
    accelerometer  |mean a| against the gravity reference (the calibration's
                   gravity target), in mg; per-axis noise
    gyro           mean rate (the residual bias after calibration); noise
    magnetometer   |mean m| (F), its component along gravity (Z, positive
                   down), horizontal component (H) and dip, against the
                   expected field (WMM at the configured place); noise
  The same means are kept for the uncalibrated samples, so the effect of
  the calibration is visible.

  Rotation. Between two poses whose gravity or field directions differ by
  at least min_tilt_change, the gyro samples recorded in between, less the
  earlier pose's mean rate (the bias over the motion), are integrated to
  predict the later pose's gravity and field directions from the earlier
  pose's. The prediction error is reported, and a gyro scale factor is
  fitted: the factor on every rate that makes the prediction match. The
  motion is integrated when the next pose starts, at the measured rates and
  at rates 1 % lower and higher, and one Gauss-Newton step fits the factor
  to the pose means; any path, flips included, is handled. The factor is
  reported only when the rotation constrains it (a turn about gravity says
  nothing to the accelerometer). A motion longer than the sample buffer
  (10 s at 200 Hz) is not checked.

  In motion. |m| of every magnetometer readout: a well-calibrated
  magnetometer reads the same |m| in every orientation. The six body
  directions the field has pointed along are counted, so coverage is
  visible.
*/

#include <algorithm>
#include <cmath>
#include <cstdint>

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

namespace imu_cal {

constexpr double ACCURACY_PI = 3.14159265358979323846;

struct AccuracyReference {
  double gravity_mps2 = NAN;    // |a| a still, calibrated accelerometer should read
  double g_std_mps2 = 9.80665;  // 1 g, for mg
  bool have_field = false;      // the field values below are set
  double f_uT = NAN, h_uT = NAN, z_uT = NAN, inclination_deg = NAN;
};

struct AccuracyConfig {
  int block_samples = 50;
  double accel_std_mps2 = 0.10;  // ~10 mg per axis
  double gyro_std_rad_s = 0.5 * ACCURACY_PI / 180.0;
  double gyro_mean_rad_s = 5.0 * ACCURACY_PI / 180.0;
  double direction_rad = 1.0 * ACCURACY_PI / 180.0;
  int min_blocks = 8;    // 2 s
  int max_blocks = 120;  // 30 s
  double min_tilt_change_rad = 20.0 * ACCURACY_PI / 180.0;
  double distinct_pose_rad = 30.0 * ACCURACY_PI / 180.0;
  // The fitted scale is reported when a 1 % scale change moves the
  // predicted direction by at least this much.
  double min_scale_sensitivity_rad = 0.1 * ACCURACY_PI / 180.0;
};

// Mean and variance of 3-vectors (Welford, mergeable).
struct Vec3Stats {
  uint32_t n = 0;
  Eigen::Vector3d mean = Eigen::Vector3d::Zero(), m2 = Eigen::Vector3d::Zero();
  void add(const Eigen::Vector3d& v) {
    ++n;
    const Eigen::Vector3d d = v - mean;
    mean += d / double(n);
    m2 += d.cwiseProduct(v - mean);
  }
  void merge(const Vec3Stats& o) {
    if (!o.n) return;
    if (!n) { *this = o; return; }
    const double na = n, nb = o.n, nt = na + nb;
    const Eigen::Vector3d d = o.mean - mean;
    mean += d * (nb / nt);
    m2 += o.m2 + d.cwiseProduct(d) * (na * nb / nt);
    n += o.n;
  }
  Eigen::Vector3d stddev() const {
    return n > 1 ? Eigen::Vector3d((m2 / double(n - 1)).cwiseMax(0.0).cwiseSqrt()) : Eigen::Vector3d::Constant(NAN);
  }
};

struct StillPose {
  uint32_t index = 0;
  double duration_s = 0;
  uint32_t samples = 0, mag_samples = 0;
  double temp_c = NAN;
  Eigen::Vector3d accel_mean = Eigen::Vector3d::Constant(NAN), accel_std = Eigen::Vector3d::Constant(NAN);
  Eigen::Vector3d accel_raw_mean = Eigen::Vector3d::Constant(NAN);
  Eigen::Vector3d gyro_mean = Eigen::Vector3d::Constant(NAN), gyro_std = Eigen::Vector3d::Constant(NAN);
  Eigen::Vector3d gyro_raw_mean = Eigen::Vector3d::Constant(NAN);
  Eigen::Vector3d mag_mean = Eigen::Vector3d::Constant(NAN), mag_std = Eigen::Vector3d::Constant(NAN);
  Eigen::Vector3d mag_raw_mean = Eigen::Vector3d::Constant(NAN);
  // Derived.
  double accel_norm_mps2 = NAN, accel_error_mg = NAN, accel_raw_error_mg = NAN, accel_noise_mg = NAN;
  double roll_deg = NAN, pitch_deg = NAN;
  double gyro_bias_dps = NAN, gyro_raw_bias_dps = NAN, gyro_noise_dps = NAN;
  double f_uT = NAN, z_uT = NAN, h_uT = NAN, inclination_deg = NAN, mag_noise_uT = NAN, f_raw_uT = NAN;
  double f_error_pct = NAN, z_error_uT = NAN, inclination_error_deg = NAN;
  bool distinct = false;  // gravity direction new among the poses so far
  // Sample counter at the pose's first and one past its last sample.
  uint64_t first_sample = 0, end_sample = 0;
};

struct RotationCheck {
  uint32_t from_pose = 0, to_pose = 0;
  double tilt_change_deg = NAN, field_change_deg = NAN, rotation_deg = NAN, motion_s = NAN;
  double accel_error_deg = NAN, mag_error_deg = NAN;      // prediction error, gyro as measured
  double gyro_scale_accel_pct = NAN, gyro_scale_mag_pct = NAN;  // fitted gyro scale error, NAN if unconstrained
  double accel_residual_deg = NAN, mag_residual_deg = NAN;      // error left at the fitted scale
};

struct AccuracySummary {
  uint32_t poses = 0, distinct_poses = 0, rotation_checks = 0, scale_fits = 0, mag_scale_fits = 0;
  // Over poses.
  double accel_error_mg_min = NAN, accel_error_mg_max = NAN, accel_noise_mg_max = NAN;
  double gyro_bias_dps_max = NAN, gyro_noise_dps_max = NAN;
  double f_error_pct_min = NAN, f_error_pct_max = NAN;
  double inclination_deg_min = NAN, inclination_deg_max = NAN, inclination_error_deg_max = NAN;
  double z_error_uT_max = NAN, mag_noise_uT_max = NAN;
  double temp_c_min = NAN, temp_c_max = NAN;
  // Over rotation checks.
  double gyro_scale_pct_min = NAN, gyro_scale_pct_max = NAN;          // fits to the gravity direction
  double gyro_scale_mag_pct_min = NAN, gyro_scale_mag_pct_max = NAN;  // fits to the field direction
  double rotation_accel_error_deg_max = NAN, rotation_mag_error_deg_max = NAN;
  // In motion.
  uint32_t mag_samples = 0;
  double mag_norm_mean_uT = NAN, mag_norm_std_pct = NAN, mag_norm_min_uT = NAN, mag_norm_max_uT = NAN;
  int mag_faces = 0;  // of 6 body directions the field has pointed along
};

class ImuAccuracy {
public:
  using Vector3f = Eigen::Vector3f;
  static constexpr int RING = 2048;  // gyro samples kept for the rotation fit (10 s at 200 Hz)

  explicit ImuAccuracy(const AccuracyConfig& cfg = AccuracyConfig{}) : cfg_(cfg) { reset(); }

  void reset() {
    block_ = Block{};
    pending_ = Block{};
    have_pending_ = in_run_ = false;
    seg_ = Segment{};
    last_quiet_dir_ = Eigen::Vector3d::Zero();
    samples_ = 0;
    have_last_pose_ = have_new_pose_ = have_new_check_ = false;
    motion_ = Motion{};
    poses_ = 0;
    n_dirs_ = 0;
    summary_ = AccuracySummary{};
    norm_n_ = 0;
    norm_mean_ = norm_m2_ = 0;
    faces_ = 0;
  }

  void setReference(const AccuracyReference& ref) { ref_ = ref; }
  const AccuracyReference& reference() const { return ref_; }

  // One accelerometer+gyro frame: calibrated and uncalibrated, body frame,
  // m/s^2 and rad/s; dt_s is the time since the previous frame.
  void addFrame(const Vector3f& a, const Vector3f& w, const Vector3f& a_raw, const Vector3f& w_raw, double dt_s,
                float temp_c) {
    if (!a.allFinite() || !w.allFinite()) return;
    const Eigen::Vector3d ad = a.cast<double>(), wd = w.cast<double>();
    if (block_.n == 0) block_.first_sample = samples_;
    RingSample& r = ring_[samples_ % RING];
    r.w = w;
    r.dt = std::isfinite(dt_s) && dt_s > 0 && dt_s < 0.1 ? float(dt_s) : 0.0f;
    ++samples_;
    block_.accel.add(ad);
    block_.gyro.add(wd);
    if (a_raw.allFinite()) block_.accel_raw.add(a_raw.cast<double>());
    if (w_raw.allFinite()) block_.gyro_raw.add(w_raw.cast<double>());
    if (std::isfinite(temp_c)) { block_.temp_sum += temp_c; ++block_.temp_n; }
    block_.duration_s += r.dt;
    if (++block_.n >= cfg_.block_samples) endBlock();
  }

  // One magnetometer readout, calibrated and uncalibrated, body frame, uT.
  void addMag(const Vector3f& m, const Vector3f& m_raw) {
    if (!m.allFinite()) return;
    const Eigen::Vector3d md = m.cast<double>();
    block_.mag.add(md);
    if (m_raw.allFinite()) block_.mag_raw.add(m_raw.cast<double>());
    const double norm = md.norm();
    if (norm <= 0) return;
    ++norm_n_;
    const double d = norm - norm_mean_;
    norm_mean_ += d / double(norm_n_);
    norm_m2_ += d * (norm - norm_mean_);
    AccuracySummary& s = summary_;
    s.mag_samples = norm_n_;
    s.mag_norm_min_uT = std::isfinite(s.mag_norm_min_uT) ? std::min(s.mag_norm_min_uT, norm) : norm;
    s.mag_norm_max_uT = std::isfinite(s.mag_norm_max_uT) ? std::max(s.mag_norm_max_uT, norm) : norm;
    s.mag_norm_mean_uT = norm_mean_;
    s.mag_norm_std_pct = norm_n_ > 1 ? 100.0 * std::sqrt(norm_m2_ / double(norm_n_ - 1)) / norm_mean_ : NAN;
    const Eigen::Vector3d u = md / norm;
    for (int i = 0; i < 3; ++i) {
      if (u[i] > 0.7) faces_ |= 1u << (2 * i);
      if (u[i] < -0.7) faces_ |= 1u << (2 * i + 1);
    }
    int count = 0;
    for (int i = 0; i < 6; ++i) count += (faces_ >> i) & 1u;
    s.mag_faces = count;
  }

  // A completed pose, and the rotation check it made with the previous one,
  // each returned once.
  bool takePose(StillPose& out) {
    if (!have_new_pose_) return false;
    have_new_pose_ = false;
    out = new_pose_;
    return true;
  }
  bool takeRotationCheck(RotationCheck& out) {
    if (!have_new_check_) return false;
    have_new_check_ = false;
    out = new_check_;
    return true;
  }

  bool still() const { return in_run_; }
  double stillSeconds() const { return in_run_ ? seg_.duration_s + (have_pending_ ? pending_.duration_s : 0.0) : 0.0; }
  const AccuracySummary& summary() const { return summary_; }
  const StillPose& lastPose() const { return last_pose_; }
  bool haveLastPose() const { return have_last_pose_; }

  // Derived values of a pose against a reference.
  static void derive(StillPose& p, const AccuracyReference& ref) {
    constexpr double R2D = 180.0 / ACCURACY_PI;
    p.accel_norm_mps2 = p.accel_mean.norm();
    p.accel_error_mg = 1000.0 * (p.accel_norm_mps2 - ref.gravity_mps2) / ref.g_std_mps2;
    p.accel_raw_error_mg = 1000.0 * (p.accel_raw_mean.norm() - ref.gravity_mps2) / ref.g_std_mps2;
    p.accel_noise_mg = 1000.0 * p.accel_std.maxCoeff() / ref.g_std_mps2;
    p.roll_deg = std::atan2(-p.accel_mean.y(), -p.accel_mean.z()) * R2D;
    p.pitch_deg = std::atan2(p.accel_mean.x(), std::hypot(p.accel_mean.y(), p.accel_mean.z())) * R2D;
    p.gyro_bias_dps = p.gyro_mean.norm() * R2D;
    p.gyro_raw_bias_dps = p.gyro_raw_mean.norm() * R2D;
    p.gyro_noise_dps = p.gyro_std.maxCoeff() * R2D;
    if (p.mag_samples >= 2 && p.mag_mean.allFinite()) {
      const Eigen::Vector3d down = -p.accel_mean.normalized();  // the accelerometer reads -g at rest
      p.f_uT = p.mag_mean.norm();
      p.f_raw_uT = p.mag_raw_mean.allFinite() ? p.mag_raw_mean.norm() : NAN;
      p.z_uT = p.mag_mean.dot(down);
      p.h_uT = std::sqrt(std::max(0.0, p.f_uT * p.f_uT - p.z_uT * p.z_uT));
      p.inclination_deg = std::atan2(p.z_uT, p.h_uT) * R2D;
      p.mag_noise_uT = p.mag_std.maxCoeff();
      if (ref.have_field) {
        p.f_error_pct = 100.0 * (p.f_uT / ref.f_uT - 1.0);
        p.z_error_uT = p.z_uT - ref.z_uT;
        p.inclination_error_deg = p.inclination_deg - ref.inclination_deg;
      }
    }
  }

private:
  struct RingSample { Eigen::Vector3f w = Eigen::Vector3f::Zero(); float dt = 0; };
  struct Block {
    int n = 0;
    Vec3Stats accel, gyro, mag, accel_raw, gyro_raw, mag_raw;
    double temp_sum = 0;
    int temp_n = 0;
    double duration_s = 0;
    uint64_t first_sample = 0;
  };
  struct Segment {
    int blocks = 0;
    Vec3Stats accel, gyro, mag, accel_raw, gyro_raw, mag_raw;
    double temp_sum = 0;
    int temp_n = 0;
    double duration_s = 0;
    uint64_t first_sample = 0, end_sample = 0;
  };

  bool quiet(const Block& b) const {
    if (!(b.accel.stddev().maxCoeff() < cfg_.accel_std_mps2) || !(b.gyro.stddev().maxCoeff() < cfg_.gyro_std_rad_s))
      return false;
    if (b.gyro.mean.norm() >= cfg_.gyro_mean_rad_s) return false;
    if (in_run_) {
      const double c = std::clamp(b.accel.mean.normalized().dot(last_quiet_dir_), -1.0, 1.0);
      if (std::acos(c) >= cfg_.direction_rad) return false;
    }
    return true;
  }

  void merge(const Block& b) {
    if (seg_.blocks == 0) {
      seg_.first_sample = b.first_sample;
      snapshotMotion(b.first_sample);
    }
    seg_.end_sample = b.first_sample + uint64_t(b.n);
    seg_.accel.merge(b.accel);
    seg_.gyro.merge(b.gyro);
    seg_.mag.merge(b.mag);
    seg_.accel_raw.merge(b.accel_raw);
    seg_.gyro_raw.merge(b.gyro_raw);
    seg_.mag_raw.merge(b.mag_raw);
    seg_.temp_sum += b.temp_sum;
    seg_.temp_n += b.temp_n;
    seg_.duration_s += b.duration_s;
    ++seg_.blocks;
  }

  void endBlock() {
    if (quiet(block_)) {
      last_quiet_dir_ = block_.accel.mean.normalized();
      if (!in_run_) {
        in_run_ = true;  // the first quiet block is settling: dropped
      } else {
        if (have_pending_) merge(pending_);
        pending_ = block_;
        have_pending_ = true;
        if (seg_.blocks >= cfg_.max_blocks) finishSegment();
      }
    } else {
      if (in_run_) finishSegment();  // the pending block may hold the start of motion: dropped
      in_run_ = false;
      have_pending_ = false;
    }
    block_ = Block{};
  }

  void finishSegment() {
    const Segment s = seg_;
    seg_ = Segment{};
    if (s.blocks < cfg_.min_blocks) return;
    StillPose p;
    p.index = ++poses_;
    p.duration_s = s.duration_s;
    p.samples = s.accel.n;
    p.mag_samples = s.mag.n;
    p.temp_c = s.temp_n ? s.temp_sum / s.temp_n : NAN;
    p.accel_mean = s.accel.mean;
    p.accel_std = s.accel.stddev();
    p.accel_raw_mean = s.accel_raw.n ? s.accel_raw.mean : Eigen::Vector3d::Constant(NAN);
    p.gyro_mean = s.gyro.mean;
    p.gyro_std = s.gyro.stddev();
    p.gyro_raw_mean = s.gyro_raw.n ? s.gyro_raw.mean : Eigen::Vector3d::Constant(NAN);
    if (s.mag.n) {
      p.mag_mean = s.mag.mean;
      p.mag_std = s.mag.stddev();
      p.mag_raw_mean = s.mag_raw.n ? s.mag_raw.mean : Eigen::Vector3d::Constant(NAN);
    }
    p.first_sample = s.first_sample;
    p.end_sample = s.end_sample;
    derive(p, ref_);

    const Eigen::Vector3d dir = p.accel_mean.normalized();
    p.distinct = true;
    for (int i = 0; i < n_dirs_; ++i)
      if (std::acos(std::clamp(dir.dot(dirs_[i]), -1.0, 1.0)) < cfg_.distinct_pose_rad) p.distinct = false;
    if (p.distinct && n_dirs_ < MAX_DIRS) dirs_[n_dirs_++] = dir;

    if (have_last_pose_) check(last_pose_, p);
    accumulate(p);
    last_pose_ = p;
    have_last_pose_ = true;
    new_pose_ = p;
    have_new_pose_ = true;
  }

  // Rotation from body frame at sample `from` to body frame at `to`, the
  // recorded rates less `bias` and scaled by k, as the operator taking a
  // vector's components at `from` to its components at `to` (for a vector
  // fixed in space).
  Eigen::Quaternionf carry(uint64_t from, uint64_t to, const Eigen::Vector3f& bias, float k) const {
    Eigen::Quaternionf q = Eigen::Quaternionf::Identity();  // body(to) -> body(from)
    for (uint64_t i = from; i < to; ++i) {
      const RingSample& r = ring_[i % RING];
      const Eigen::Vector3f th = (r.w - bias) * (k * r.dt);
      const float ang = th.norm();
      if (ang > 1e-9f) q = q * Eigen::Quaternionf(Eigen::AngleAxisf(ang, th / ang));
    }
    return q.normalized().conjugate();
  }

  // The motion since the last pose, integrated when the next pose starts
  // (its samples are about to leave the buffer) at the measured rates and
  // at rates scaled by 1 -+ H, the last pose's mean rate removed as bias.
  static constexpr float H = 0.01f;
  struct Motion {
    bool valid = false;
    uint64_t from = 0, to = 0;
    double seconds = 0;
    Eigen::Quaternionf at1 = Eigen::Quaternionf::Identity(), lo = Eigen::Quaternionf::Identity(),
                       hi = Eigen::Quaternionf::Identity();
  };

  void snapshotMotion(uint64_t to) {
    motion_ = Motion{};
    if (!have_last_pose_) return;
    const uint64_t from = last_pose_.end_sample;
    if (to <= from || to - from > uint64_t(RING) || samples_ - from > uint64_t(RING)) return;
    const Eigen::Vector3f bias = last_pose_.gyro_mean.cast<float>();
    motion_.from = from;
    motion_.to = to;
    for (uint64_t i = from; i < to; ++i) motion_.seconds += ring_[i % RING].dt;
    motion_.at1 = carry(from, to, bias, 1.0f);
    motion_.lo = carry(from, to, bias, 1.0f - H);
    motion_.hi = carry(from, to, bias, 1.0f + H);
    motion_.valid = true;
  }

  // Gyro scale factor making the carried direction `va` match `vb` (one
  // Gauss-Newton step from 1, central difference); NAN when the rotation
  // does not constrain it.
  float fitScale(const Eigen::Vector3f& va, const Eigen::Vector3f& vb, float& residual_rad) const {
    residual_rad = NAN;
    const Eigen::Vector3f p0 = motion_.at1 * va;
    const Eigen::Vector3f dp = (motion_.hi * va - motion_.lo * va) / (2.0f * H);
    if (dp.norm() * 0.01f < float(cfg_.min_scale_sensitivity_rad)) return NAN;
    const float k = 1.0f + dp.dot(vb - p0) / dp.squaredNorm();
    if (!(k > 0.5f && k < 1.5f)) return NAN;
    residual_rad = std::acos(std::clamp((p0 + (k - 1.0f) * dp).normalized().dot(vb), -1.0f, 1.0f));
    return k;
  }

  void check(const StillPose& a, const StillPose& b) {
    constexpr double R2D = 180.0 / ACCURACY_PI;
    const Eigen::Vector3d da = a.accel_mean.normalized(), db = b.accel_mean.normalized();
    const double tilt = std::acos(std::clamp(da.dot(db), -1.0, 1.0));
    const bool have_mag = a.mag_mean.allFinite() && b.mag_mean.allFinite();
    const double turn = have_mag ? std::acos(std::clamp(a.mag_mean.normalized().dot(b.mag_mean.normalized()), -1.0, 1.0))
                                 : 0.0;
    if (tilt < cfg_.min_tilt_change_rad && turn < cfg_.min_tilt_change_rad) return;
    if (!motion_.valid || motion_.from != a.end_sample || motion_.to != b.first_sample) return;  // too long
    RotationCheck c;
    c.from_pose = a.index;
    c.to_pose = b.index;
    c.tilt_change_deg = tilt * R2D;
    c.motion_s = motion_.seconds;
    const Eigen::Quaternionf& rel = motion_.at1;
    c.rotation_deg = 2.0 * std::acos(std::clamp(double(std::fabs(rel.w())), 0.0, 1.0)) * R2D;
    const Eigen::Vector3f fa = da.cast<float>(), fb = db.cast<float>();
    c.accel_error_deg = std::acos(std::clamp(double((rel * fa).dot(fb)), -1.0, 1.0)) * R2D;
    float res = NAN;
    const float k = fitScale(fa, fb, res);
    if (std::isfinite(k)) {
      c.gyro_scale_accel_pct = (1.0 / k - 1.0) * 100.0;
      c.accel_residual_deg = res * R2D;
    }
    if (have_mag) {
      c.field_change_deg = turn * R2D;
      const Eigen::Vector3f ma = a.mag_mean.normalized().cast<float>(), mb = b.mag_mean.normalized().cast<float>();
      c.mag_error_deg = std::acos(std::clamp(double((rel * ma).dot(mb)), -1.0, 1.0)) * R2D;
      const float km = fitScale(ma, mb, res);
      if (std::isfinite(km)) {
        c.gyro_scale_mag_pct = (1.0 / km - 1.0) * 100.0;
        c.mag_residual_deg = res * R2D;
      }
    }
    AccuracySummary& s = summary_;
    ++s.rotation_checks;
    auto lo = [](double& m, double v) { if (std::isfinite(v)) m = std::isfinite(m) ? std::min(m, v) : v; };
    auto hi = [](double& m, double v) { if (std::isfinite(v)) m = std::isfinite(m) ? std::max(m, v) : v; };
    if (std::isfinite(c.gyro_scale_accel_pct)) ++s.scale_fits;
    if (std::isfinite(c.gyro_scale_mag_pct)) ++s.mag_scale_fits;
    lo(s.gyro_scale_pct_min, c.gyro_scale_accel_pct);
    hi(s.gyro_scale_pct_max, c.gyro_scale_accel_pct);
    lo(s.gyro_scale_mag_pct_min, c.gyro_scale_mag_pct);
    hi(s.gyro_scale_mag_pct_max, c.gyro_scale_mag_pct);
    hi(s.rotation_accel_error_deg_max, c.accel_error_deg);
    hi(s.rotation_mag_error_deg_max, c.mag_error_deg);
    new_check_ = c;
    have_new_check_ = true;
  }

  void accumulate(const StillPose& p) {
    auto lo = [](double& m, double v) { if (std::isfinite(v)) m = std::isfinite(m) ? std::min(m, v) : v; };
    auto hi = [](double& m, double v) { if (std::isfinite(v)) m = std::isfinite(m) ? std::max(m, v) : v; };
    AccuracySummary& s = summary_;
    ++s.poses;
    if (p.distinct) ++s.distinct_poses;
    lo(s.accel_error_mg_min, p.accel_error_mg);
    hi(s.accel_error_mg_max, p.accel_error_mg);
    hi(s.accel_noise_mg_max, p.accel_noise_mg);
    hi(s.gyro_bias_dps_max, p.gyro_bias_dps);
    hi(s.gyro_noise_dps_max, p.gyro_noise_dps);
    lo(s.f_error_pct_min, p.f_error_pct);
    hi(s.f_error_pct_max, p.f_error_pct);
    lo(s.inclination_deg_min, p.inclination_deg);
    hi(s.inclination_deg_max, p.inclination_deg);
    hi(s.inclination_error_deg_max, std::fabs(p.inclination_error_deg));
    hi(s.z_error_uT_max, std::fabs(p.z_error_uT));
    hi(s.mag_noise_uT_max, p.mag_noise_uT);
    lo(s.temp_c_min, p.temp_c);
    hi(s.temp_c_max, p.temp_c);
  }

  static constexpr int MAX_DIRS = 26;
  AccuracyConfig cfg_;
  AccuracyReference ref_{};
  Block block_{}, pending_{};
  bool have_pending_ = false, in_run_ = false;
  Segment seg_{};
  Eigen::Vector3d last_quiet_dir_ = Eigen::Vector3d::Zero();
  RingSample ring_[RING];
  uint64_t samples_ = 0;
  StillPose last_pose_{}, new_pose_{};
  RotationCheck new_check_{};
  Motion motion_{};
  bool have_last_pose_ = false, have_new_pose_ = false, have_new_check_ = false;
  uint32_t poses_ = 0;
  Eigen::Vector3d dirs_[MAX_DIRS];
  int n_dirs_ = 0;
  AccuracySummary summary_{};
  uint32_t norm_n_ = 0;
  double norm_mean_ = 0, norm_m2_ = 0;
  unsigned faces_ = 0;
};

}  // namespace imu_cal
