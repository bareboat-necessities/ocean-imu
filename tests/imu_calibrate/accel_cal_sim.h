#pragma once

// Simulated AtomS3R accelerometer calibration sessions: sensor truth, a hand /
// table user model, temperature history, and sample-stream defects. Used by
// accel_cal-test.cpp (deterministic multi-seed campaign) and by the replay
// tool's self-test. Everything is seeded; no std::*_distribution is used so
// results do not depend on the standard library implementation.

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "imu_calibrate/AccelCalCapture.h"

namespace accel_sim {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kGStd = 9.80665;

// Deterministic RNG (splitmix64 seeded xorshift) with explicit transforms.
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed * 0x9E3779B97F4A7C15ull + 0x2545F4914F6CDD1Dull) { next(); next(); }
  uint64_t next() {
    uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  double uni() { return (double)(next() >> 11) * (1.0 / 9007199254740992.0); }
  double uni(double a, double b) { return a + (b - a) * uni(); }
  double normal() {
    double u1 = uni(); if (u1 < 1e-300) u1 = 1e-300;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * uni());
  }
  double normal(double sd) { return sd * normal(); }
  bool chance(double p) { return uni() < p; }
  Vec3 nvec(double sd) { return Vec3(normal(sd), normal(sd), normal(sd)); }
  Vec3 unit() { Vec3 v = nvec(1.0); return v / v.norm(); }
};

static inline Mat3 expso3(const Vec3& w) {
  const double a = w.norm();
  if (a < 1e-12) return Mat3::Identity();
  return Eigen::AngleAxisd(a, w / a).toRotationMatrix();
}

static inline Vec3 logso3(const Mat3& R) {
  Eigen::AngleAxisd aa(R);
  return aa.angle() * aa.axis();
}

// Sensor truth. raw = S(T)^-1 * M * f_body + b(T) + noise, a_true = f_body.
struct SensorTruth {
  Mat3 S = Mat3::Identity();        // symmetric: the calibration the fitter should find
  Mat3 R_mis = Mat3::Identity();    // optional non-symmetric part (sensor misalignment)
  Vec3 b0 = Vec3::Zero();           // bias at T0
  Vec3 k = Vec3::Zero();            // bias slope, m/s^2/degC
  Mat3 dS_dT = Mat3::Zero();        // thermal scale change (model mismatch)
  double T0 = 25.0;
  double g_local = kGStd;
  double noise = 0.010;             // white, per axis, m/s^2
  double ar_sd = 0.003, ar_rho = 0.995;  // correlated noise
  Vec3 gyro_bias = Vec3::Zero();
  double gyro_noise = 0.0015;

  Mat3 S_at(double T) const { return S + dS_dT * (T - T0); }
  Vec3 bias_at(double T) const { return b0 + k * (T - T0); }
  // Symmetric calibration the gravity-norm fit can recover at temperature T.
  Mat3 S_polar_at(double T) const {
    const Mat3 A = S_at(T) * R_mis.transpose();  // a_true = R_mis^T-rotated
    Eigen::SelfAdjointEigenSolver<Mat3> es(A.transpose() * A);
    return es.eigenvectors() * es.eigenvalues().cwiseSqrt().asDiagonal() * es.eigenvectors().transpose();
  }

  static SensorTruth random(Rng& r, double bias_sd = 0.15, double scale_sd = 0.01, double cross_sd = 0.004,
                            double k_sd = 0.002) {
    SensorTruth t;
    Mat3 S = Mat3::Identity();
    for (int i = 0; i < 3; ++i) S(i, i) += r.normal(scale_sd);
    for (int i = 0; i < 3; ++i)
      for (int j = i + 1; j < 3; ++j) { const double c = r.normal(cross_sd); S(i, j) = S(j, i) = c; }
    t.S = S;
    t.b0 = r.nvec(bias_sd);
    t.k = r.nvec(k_sd);
    t.gyro_bias = r.nvec(0.01);
    return t;
  }
};

struct TempProfile {
  double T_start = 24.0;
  double dT = 1.0;         // total rise
  double tau_s = 300.0;
  double read_noise = 0.02;
  bool all_nan = false;
  double nan_prob = 0.0;
  // Optional override: T(t) (seconds since power-on).
  std::function<double(double)> custom;
  double at(double t) const {
    if (custom) return custom(t);
    return T_start + dT * (1.0 - std::exp(-t / tau_s));
  }
};

// Hand held "still": three rotational wobble modes (0.2-0.6 deg at 0.2-1.2 Hz,
// up to ~5 deg/s RMS), slow drift, 8-12 Hz translational tremor and 0.3-1 Hz
// sway of 0.5-2 mm. `scale` multiplies every term (adverse users).
struct HandModel {
  double wobble_deg_lo = 0.2, wobble_deg_hi = 0.6;
  double wobble_hz_lo = 0.2, wobble_hz_hi = 1.2;
  double drift_dps_lo = 0.02, drift_dps_hi = 0.15;
  double tremor_lo = 0.03, tremor_hi = 0.10;      // m/s^2
  double sway_mm_lo = 0.5, sway_mm_hi = 2.0;
  double settle_deg = 2.0;
  double scale = 1.0;                             // adverse multiplier
};

struct UserModel {
  double react_lo = 1.2, react_hi = 3.0;          // read + tap
  double move_lo = 1.0, move_hi = 4.0;            // tap -> placed
  double p_slow = 0.1;                            // slow placement 5-9 s
  double face_err_table_deg = 2.0;
  double face_err_hand_deg = 5.0;
  double tilt_err_deg = 12.0;                     // rms of tilted placements
  double p_first_wrong_tilt = 0.15;               // tilts only one axis, corrects on "Check the pose"
  double p_regrip = 0.15;                         // quick 3-8 deg adjustment during a hand-held hold
  double shock_rate_per_s = 0.02;                 // table knocks / taps
  double shock_lo = 2.0, shock_hi = 6.0;
  double p_hand_face[6] = {0.0, 0.0, 0.2, 1.0, 0.3, 0.3};  // USB down needs a hand (cable)
  double p_hand_tilt = 0.7;
  HandModel hand;
};

struct StreamModel {
  double rate_hz = 100.0;
  double jitter_ms = 0.5;
  double p_stale = 0.0;           // exact repeat of the previous sample, new timestamp
  double p_gap = 0.0;             // start of a 30-150 ms dropout
};

struct Scenario {
  std::string name;
  UserModel user;
  StreamModel stream;
  TempProfile temp;
  double wizard_start_s = 20.0;   // power-on to wizard start
  double gyro_mag_gap_lo = 70.0, gyro_mag_gap_hi = 110.0;  // full wizard gyro+mag stage duration
  bool accel_only = false;
  // truth generation
  double bias_sd = 0.15, scale_sd = 0.01, cross_sd = 0.004, k_sd = 0.002;
  double mis_deg = 0.0;           // non-symmetric misalignment magnitude
  double dSdT = 0.0;              // thermal scale coefficient (diag), 1/degC
  double g_local = kGStd;
  Vec3 k_override = Vec3::Constant(NAN);
  Vec3 bias_offset = Vec3::Zero();  // power-cycle offset added to the truth bias
};

// Device frame (top, right, out) -> body, documented mapping.
static inline Vec3 docUpBody(int top, int right, int out) {
  return (Vec3(1, 0, 0) * top + Vec3(0, 1, 0) * right + Vec3(0, 0, -1) * out).normalized();
}

// Minimal rotation R (body->world) with R * d = world up (NED -z), then yaw.
static inline Mat3 attitudeForUp(const Vec3& d_body, double yaw) {
  const Vec3 up_w(0, 0, -1);
  Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(d_body, up_w);
  return Eigen::AngleAxisd(yaw, Vec3::UnitZ()).toRotationMatrix() * q.toRotationMatrix();
}

// The physical world: attitude history, hand motion, temperature, sensor.
class World {
public:
  World(const SensorTruth& truth, const Scenario& sc, uint64_t seed) : truth_(truth), sc_(sc), rng_(seed) {
    R_ = attitudeForUp(docUpBody(0, 0, 1), rng_.uni(0, 2 * kPi));
    hold_R_ = R_;
    ar_.setZero();
    t_ = sc.wizard_start_s;
  }

  double t() const { return t_; }
  void advance(double dt) { t_ += dt; }
  const SensorTruth& truth() const { return truth_; }
  Rng& rng() { return rng_; }

  // Begin moving to a pose with measured-up direction d_body (unit, body).
  void placeAt(const Vec3& d_body, bool hand, double err_deg, double t_start, double t_move) {
    move_from_ = attitudeAt_(t_start);
    Mat3 target = attitudeForUp(d_body, rng_.uni(0, 2 * kPi));
    target = expso3(rng_.nvec(err_deg * kDeg / std::sqrt(2.0))) * target;
    move_to_ = target;
    move_t0_ = t_start;
    move_T_ = t_move;
    moving_ = true;
    hand_ = hand;
    lift_axis_ = rng_.unit();
    lift_amp_ = rng_.uni(1.0, 3.0);
    // hold motion parameters
    const HandModel& h = sc_.user.hand;
    const double sc = hand ? h.scale : 0.05;
    for (int j = 0; j < 3; ++j) {
      wob_amp_[j] = rng_.nvec(1.0).normalized() * rng_.uni(h.wobble_deg_lo, h.wobble_deg_hi) * kDeg * sc;
      wob_f_[j] = rng_.uni(h.wobble_hz_lo, h.wobble_hz_hi);
      wob_ph_[j] = rng_.uni(0, 2 * kPi);
      trem_f_[j] = rng_.uni(8.0, 12.0);
      trem_ph_[j] = rng_.uni(0, 2 * kPi);
      sway_f_[j] = rng_.uni(0.3, 1.0);
      sway_ph_[j] = rng_.uni(0, 2 * kPi);
    }
    trem_amp_ = hand ? rng_.uni(h.tremor_lo, h.tremor_hi) * h.scale : 0.0;
    sway_amp_ = hand ? rng_.uni(h.sway_mm_lo, h.sway_mm_hi) * 1e-3 * h.scale : 0.0;
    drift_ = Vec3::Zero();
    if (hand) drift_ = rng_.unit() * rng_.uni(h.drift_dps_lo, h.drift_dps_hi) * kDeg * h.scale;
    settle_amp_ = rng_.unit() * (hand ? h.settle_deg : 0.3) * kDeg;
    regrip_t_ = (hand && rng_.chance(sc_.user.p_regrip)) ? t_start + t_move + rng_.uni(2.0, 6.0) : -1;
    regrip_v_ = rng_.unit() * rng_.uni(3.0, 8.0) * kDeg;
  }

  // Specific force and angular rate at time t (body frame, truth).
  void truthAt(double t, Vec3& f_b, Vec3& w_b) {
    const double h = 1e-3;
    const Mat3 R0 = attitudeAt_(t - h), R1 = attitudeAt_(t + h);
    const Mat3 R = attitudeAt_(t);
    const Vec3 w_w = logso3(R1 * R0.transpose()) / (2 * h);
    w_b = R.transpose() * w_w;
    const Vec3 a_w = linAccAt_(t);
    const Vec3 g_w(0, 0, truth_.g_local);
    f_b = R.transpose() * (a_w - g_w);
  }

  // One IMU sample at the current time.
  void sampleAt(double t, Vec3& a_raw, Vec3& w_raw, double& tempC_true) {
    Vec3 f_b, w_b;
    truthAt(t, f_b, w_b);
    const double T = sc_.temp.at(t);
    tempC_true = T;
    for (int j = 0; j < 3; ++j) ar_(j) = truth_.ar_rho * ar_(j) + std::sqrt(1 - truth_.ar_rho * truth_.ar_rho) * rng_.normal(truth_.ar_sd);
    const Mat3 A = truth_.S_at(T) * truth_.R_mis.transpose();
    a_raw = A.inverse() * f_b + truth_.bias_at(T) + rng_.nvec(truth_.noise) + ar_;
    w_raw = w_b + truth_.gyro_bias + rng_.nvec(truth_.gyro_noise);
  }

  double tempReading(double T) {
    if (sc_.temp.all_nan || (sc_.temp.nan_prob > 0 && rng_.chance(sc_.temp.nan_prob))) return NAN;
    const double v = T + rng_.normal(sc_.temp.read_noise);
    return std::round(v * 512.0) / 512.0;
  }

  void addShock(double t0, double dur, const Vec3& a) { shocks_.push_back({t0, dur, a}); }

private:
  struct Shock { double t0, dur; Vec3 a; };

  Mat3 attitudeAt_(double t) const {
    if (!moving_) return R_;
    if (t < move_t0_) return move_from_;
    const double s = (t - move_t0_) / move_T_;
    if (s < 1.0) {
      const double m = s * s * s * (10 - 15 * s + 6 * s * s);  // min-jerk
      const Vec3 d = logso3(move_to_ * move_from_.transpose());
      return expso3(d * m) * move_from_;
    }
    const double th = t - (move_t0_ + move_T_);
    Vec3 rv = Vec3::Zero();
    for (int j = 0; j < 3; ++j) rv += wob_amp_[j] * std::sin(2 * kPi * wob_f_[j] * th + wob_ph_[j]);
    rv += drift_ * th;
    rv += settle_amp_ * std::exp(-th / 0.5) * std::cos(2 * kPi * 1.5 * th);
    if (regrip_t_ > 0 && t > regrip_t_) {
      const double u = std::min(1.0, (t - regrip_t_) / 0.5);
      rv += regrip_v_ * (u * u * (3 - 2 * u));
    }
    return expso3(rv) * move_to_;
  }

  Vec3 linAccAt_(double t) const {
    Vec3 a = Vec3::Zero();
    if (moving_ && t >= move_t0_) {
      const double s = (t - move_t0_) / move_T_;
      if (s < 1.0) a += lift_axis_ * lift_amp_ * std::sin(2 * kPi * s);
      else {
        const double th = t - (move_t0_ + move_T_);
        for (int j = 0; j < 3; ++j) {
          a(j) += trem_amp_ / std::sqrt(3.0) * std::sin(2 * kPi * trem_f_[j] * th + trem_ph_[j]);
          const double w = 2 * kPi * sway_f_[j];
          a(j) += -sway_amp_ * w * w * std::sin(w * th + sway_ph_[j]);
        }
      }
    }
    for (const Shock& s : shocks_) {
      if (t >= s.t0 && t < s.t0 + s.dur) a += s.a * std::sin(kPi * (t - s.t0) / s.dur);
    }
    return a;
  }

  SensorTruth truth_;
  Scenario sc_;
  Rng rng_;
  double t_ = 0;
  Mat3 R_, hold_R_;
  Vec3 ar_;
  bool moving_ = false, hand_ = false;
  Mat3 move_from_ = Mat3::Identity(), move_to_ = Mat3::Identity();
  double move_t0_ = 0, move_T_ = 1;
  Vec3 lift_axis_ = Vec3::UnitX();
  double lift_amp_ = 0;
  Vec3 wob_amp_[3];
  double wob_f_[3] = {}, wob_ph_[3] = {}, trem_f_[3] = {}, trem_ph_[3] = {}, sway_f_[3] = {}, sway_ph_[3] = {};
  double trem_amp_ = 0, sway_amp_ = 0;
  Vec3 drift_ = Vec3::Zero(), settle_amp_ = Vec3::Zero();
  double regrip_t_ = -1;
  Vec3 regrip_v_ = Vec3::Zero();
  std::vector<Shock> shocks_;
};

// Sample stream with timing jitter, stale repeats and dropouts.
class Stream {
public:
  Stream(World& w, const StreamModel& m) : w_(w), m_(m) {}

  // Advances the world to the next sample; returns false while in a dropout.
  bool next(imu_cal::AccelRawSample& s) {
    const double dt = 1.0 / m_.rate_hz + w_.rng().normal(m_.jitter_ms * 1e-3);
    w_.advance(std::max(1e-4, dt));
    if (gap_until_ > w_.t()) return false;
    if (m_.p_gap > 0 && w_.rng().chance(m_.p_gap)) {
      gap_until_ = w_.t() + w_.rng().uni(0.03, 0.15);
      return false;
    }
    s.t_us = (uint32_t)(uint64_t)std::llround(w_.t() * 1e6);
    if (have_prev_ && m_.p_stale > 0 && w_.rng().chance(m_.p_stale)) {
      s.a = prev_.a; s.w = prev_.w; s.tempC = prev_.tempC;
      return true;
    }
    Vec3 a, wr; double T;
    w_.sampleAt(w_.t(), a, wr, T);
    s.a = a.cast<float>();
    s.w = wr.cast<float>();
    s.tempC = (float)w_.tempReading(T);
    prev_ = s; have_prev_ = true;
    return true;
  }

private:
  World& w_;
  StreamModel m_;
  double gap_until_ = -1;
  imu_cal::AccelRawSample prev_;
  bool have_prev_ = false;
};

}  // namespace accel_sim
