#pragma once

/* Copyright 2026, Mikhail Grushinskiy

   Gyroscope scale and cross-axis calibration from the accelerometer pose
   session, the method of IMU-TK (Tedaldi, Pretto, Menegatti, "A Robust and
   Easy to Implement Method for IMU Calibration without External Equipments",
   ICRA 2014) as used by INSLIB (tools/inslib_imu_tk.py): between two poses
   the accelerometer measures both gravity directions, and the gyro integrated
   over the motion between them must carry one onto the other.

   The poses are the accelerometer holds the wizard already asks for, so the
   user does nothing extra. Instead of storing the rate stream, GyroChain
   integrates the attitude together with its first-order sensitivity to the
   nine matrix terms and to a constant rate offset (preintegration Jacobians).
   Only retained accelerometer blocks contribute, summed per hold, so memory
   does not grow with session length. Each completed hold re-anchors the
   chain, which keeps every perturbation a single transition long.

   The rate bias is never estimated from pose or turn data. Integration uses
   the strict stationary capture; the fit carries a constant offset since
   that capture only as a nuisance with a prior, and does not store it.

   The fitted matrix is used only when it predicts left-out pose pairs better
   than the nominal axes (both with the same nuisance offset). Otherwise the
   gyro keeps nominal axes exactly as before; nothing here can fail a
   calibration.
*/

#include "imu_calibrate/CalibrateIMU.h"

#include <algorithm>
#include <cmath>

namespace imu_cal {

// Columns: E(i,j) at 3*i+j for w_true = (I+E)*w, then a constant rate offset.
static constexpr int kGyroAxisParams = 12;
using GyroAxisSens = Eigen::Matrix<float, 3, kGyroAxisParams>;

// Attitude since the current anchor and its sensitivity (rotation vector in
// the anchor frame per parameter).
struct GyroChainState {
  Eigen::Matrix3f R = Eigen::Matrix3f::Identity();   // body(t) -> body(anchor)
  GyroAxisSens P = GyroAxisSens::Zero();
  uint32_t anchor = 0;
};

// Sums over the retained samples of one block or hold, all in one anchor frame.
struct GyroPoseSums {
  Eigen::Vector3f Ra = Eigen::Vector3f::Zero();   // sum R * a_raw
  Eigen::Matrix3f Rs = Eigen::Matrix3f::Zero();   // sum R
  Eigen::Vector3f a = Eigen::Vector3f::Zero();    // sum a_raw
  GyroAxisSens P = GyroAxisSens::Zero();          // sum P
  uint32_t n = 0, anchor = 0;
  bool broken = false;   // a sample without attitude, or from another anchor

  void reset() { *this = GyroPoseSums{}; }
  bool ok() const { return n > 0 && !broken; }

  void add(const GyroChainState* s, const Eigen::Vector3f& a_raw) {
    if (!s) { broken = true; return; }
    if (n == 0) anchor = s->anchor;
    else if (s->anchor != anchor) broken = true;
    Ra.noalias() += s->R * a_raw;
    Rs += s->R;
    a += a_raw;
    P += s->P;
    ++n;
  }

  void merge(const GyroPoseSums& o) {
    if (o.n == 0) { broken = broken || o.broken; return; }
    if (n == 0) { const bool was = broken; *this = o; broken = broken || was; return; }
    broken = broken || o.broken || o.anchor != anchor;
    Ra += o.Ra; Rs += o.Rs; a += o.a; P += o.P; n += o.n;
  }
};

struct GyroChainCfg {
  float max_gap_s = 0.05f;   // longer sample gaps restart the chain
  float max_rate = 20.0f;    // rad/s; beyond is not a usable reading
};

class GyroChain {
public:
  GyroChainCfg cfg{};

  // Ends the chain; the next sample starts a new one with a new anchor.
  void reset() { have_ = false; ++anchor_; }

  // w: bias-corrected rate (rad/s) at t_us. Returns the state for this
  // sample, or nullptr when it cannot be integrated.
  const GyroChainState* feed(uint32_t t_us, const Eigen::Vector3f& w) {
    if (!w.allFinite() || !(w.norm() <= cfg.max_rate)) { reset(); return nullptr; }
    const Eigen::Vector3d wd = w.cast<double>();
    const uint32_t du = t_us - t_prev_;
    if (have_ && du == 0) return &out_;
    if (!have_ || du > uint32_t(cfg.max_gap_s * 1e6f) || du > 0x7FFFFFFFu) {
      start_();
    } else {
      const double dt = double(du) * 1e-6;
      const Eigen::Vector3d wm = 0.5 * (w_prev_ + wd);
      const Eigen::Vector3d th = wm * dt;
      const double t2 = th.squaredNorm(), t = std::sqrt(t2);
      double half_sinc, b, c;   // Jr = I - b[th] + c[th]^2
      if (t < 1e-4) {
        half_sinc = 0.5 - t2 / 48.0; b = 0.5 - t2 / 24.0; c = 1.0 / 6.0 - t2 / 120.0;
      } else {
        half_sinc = std::sin(0.5 * t) / t; b = (1.0 - std::cos(t)) / t2; c = (t - std::sin(t)) / (t2 * t);
      }
      const Eigen::Quaterniond dq(std::cos(0.5 * t), half_sinc * th.x(), half_sinc * th.y(), half_sinc * th.z());
      q_ = (q_ * dq).normalized();
      const Eigen::Matrix3d X = skew_(th);
      const Eigen::Matrix3d Jr = Eigen::Matrix3d::Identity() - b * X + c * X * X;
      const Eigen::Matrix3d B = q_.toRotationMatrix() * Jr * dt;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) P_.col(3 * i + j).noalias() += B.col(i) * wm[j];
      P_.template rightCols<3>() -= B;
    }
    w_prev_ = wd;
    t_prev_ = t_us;
    out_.R = q_.toRotationMatrix().cast<float>();
    out_.P = P_.cast<float>();
    out_.anchor = anchor_;
    return &out_;
  }

  // Ends a hold: returns its sums (relative to the current anchor) expressed
  // in a new anchor at the current sample, and continues the chain from
  // there. The result is marked broken when the chain did not carry the
  // hold's samples up to now.
  GyroPoseSums reanchor(const GyroPoseSums& hold) {
    GyroPoseSums out = hold;
    if (!have_ || !hold.ok() || hold.anchor != anchor_) {
      out.broken = true;
    } else {
      const Eigen::Matrix3f RaT = q_.toRotationMatrix().transpose().cast<float>();
      const GyroAxisSens Pa = P_.cast<float>();
      out.Ra = RaT * hold.Ra;
      out.Rs = RaT * hold.Rs;
      out.P = RaT * (hold.P - float(hold.n) * Pa);
    }
    ++anchor_;
    out.anchor = anchor_;
    q_.setIdentity();
    P_.setZero();
    return out;
  }

  bool active() const { return have_; }
  uint32_t anchor() const { return anchor_; }

private:
  static Eigen::Matrix3d skew_(const Eigen::Vector3d& v) {
    Eigen::Matrix3d m;
    m << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return m;
  }
  void start_() {
    have_ = true;
    ++anchor_;
    q_.setIdentity();
    P_.setZero();
  }

  bool have_ = false;
  uint32_t anchor_ = 0, t_prev_ = 0;
  Eigen::Quaterniond q_ = Eigen::Quaterniond::Identity();
  Eigen::Matrix<double, 3, kGyroAxisParams> P_ = Eigen::Matrix<double, 3, kGyroAxisParams>::Zero();
  Eigen::Vector3d w_prev_ = Eigen::Vector3d::Zero();
  GyroChainState out_{};
};

// One accelerometer hold as the fitter sees it: sums relative to the anchor
// before the hold (entry) and to the anchor it created (exit).
struct GyroHoldView {
  const GyroPoseSums* entry = nullptr;
  const GyroPoseSums* exit = nullptr;
  float tempC = NAN;
};

struct GyroAxisCfg {
  int min_pairs = 7;
  double prior_matrix = 0.02;      // 1 sigma per matrix term
  double prior_offset = 0.002;     // rad/s, drift since the stationary capture
  double noise_rad = 0.002;        // nominal per-pair direction error (prior weight)
  double max_scale_dev = 0.05;     // |S_ii - 1|
  double max_cross = 0.05;         // |S_ij|
  double min_improvement = 0.25;   // held-out RMS reduction required vs nominal axes
  double max_cv = 0.01;            // rad; held-out RMS the fitted matrix must reach
  // Inconsistent holds (motion the accelerometer gates kept): while the
  // held-out RMS exceeds outlier_rad, the hold whose removal minimises it is
  // dropped with its pairs if that brings it to drop_gain of its value.
  double outlier_rad = 0.003;
  double drop_gain = 0.5;
  int max_drop = 2;
};

enum class GyroAxisReason : uint8_t { OK = 0, NO_DATA, TOO_FEW_PAIRS, SOLVE, IMPLAUSIBLE, NOT_BETTER };

inline const char* gyroAxisReasonStr(GyroAxisReason r) {
  switch (r) {
    case GyroAxisReason::OK: return "fitted";
    case GyroAxisReason::NO_DATA: return "no turn data";
    case GyroAxisReason::TOO_FEW_PAIRS: return "too few turns";
    case GyroAxisReason::SOLVE: return "solve failed";
    case GyroAxisReason::IMPLAUSIBLE: return "implausible";
    case GyroAxisReason::NOT_BETTER: return "no improvement";
    default: return "?";
  }
}

struct GyroAxisResult {
  bool ok = false;
  GyroAxisReason reason = GyroAxisReason::NO_DATA;
  Eigen::Matrix3f S = Eigen::Matrix3f::Identity();
  Eigen::Vector3f offset = Eigen::Vector3f::Zero();   // nuisance, not applied
  int pairs = 0, dropped = 0;         // pairs used; holds dropped as inconsistent
  float cv_nominal = 0, cv_fit = 0;   // leave-one-hold-out direction RMS, rad
  float rms_fit = 0;                  // in-sample direction RMS, rad
  float max_angle = 0;                // largest gravity-direction change of a pair, rad
};

template <int MAXP = 24>
class GyroAxisFitter {
  using V3 = Eigen::Vector3d;
  using M3 = Eigen::Matrix3d;
  using J = Eigen::Matrix<double, 3, kGyroAxisParams>;
  using JF = Eigen::Matrix<float, 3, kGyroAxisParams>;   // per-pair storage
  using H = Eigen::Matrix<double, kGyroAxisParams, kGyroAxisParams>;
  using X = Eigen::Matrix<double, kGyroAxisParams, 1>;

public:
  // Pairs every hold whose entry anchor is another hold's exit anchor.
  bool fit(const GyroHoldView* holds, int n, const AccelCalibration<float>& acc, const GyroAxisCfg& cfg,
           GyroAxisResult& out) {
    out = GyroAxisResult{};
    np_ = 0;
    nh_ = std::min(n, MAXP);
    if (!holds || nh_ < 2 || !acc.ok) return finish_(out, GyroAxisReason::NO_DATA);
    for (int q = 0; q < nh_ && np_ < MAXP; ++q) {
      if (!holds[q].entry || !holds[q].entry->ok()) continue;
      for (int p = 0; p < nh_; ++p) {
        if (p == q || !holds[p].exit || !holds[p].exit->ok() || holds[p].exit->anchor != holds[q].entry->anchor)
          continue;
        V3 dp, dq;
        J Jp, Jq;
        if (!direction_(*holds[p].exit, holds[p].tempC, acc, dp, Jp) ||
            !direction_(*holds[q].entry, holds[q].tempC, acc, dq, Jq)) break;
        r0_[np_] = (dq - dp).template cast<float>();
        J_[np_] = (Jq - Jp).template cast<float>();
        hp_[np_] = p;
        hq_[np_] = q;
        use_[np_] = true;
        // Body-frame gravity direction change across the pair.
        const V3 up = holds[p].exit->a.template cast<double>().normalized();
        const V3 uq = holds[q].entry->a.template cast<double>().normalized();
        out.max_angle = std::max(out.max_angle, float(std::acos(std::max(-1.0, std::min(1.0, up.dot(uq))))));
        ++np_;
        break;
      }
    }
    out.pairs = np_;
    if (np_ == 0) return finish_(out, GyroAxisReason::NO_DATA);
    if (np_ < cfg.min_pairs) return finish_(out, GyroAxisReason::TOO_FEW_PAIRS);

    // Drop holds the others contradict, one at a time: a pose disturbed by
    // motion the accelerometer gates kept spoils both of its pairs and, through
    // the fit, the predictions of every other hold. The hold whose removal
    // leaves the most consistent set is dropped when that at least halves the
    // held-out error.
    double cv_fit = 0, cv_nom = 0;
    if (!crossValidate_(cfg, cv_fit, cv_nom)) return finish_(out, GyroAxisReason::SOLVE);
    while (out.dropped < cfg.max_drop && cv_fit > cfg.outlier_rad) {
      double best_fit = cv_fit, best_nom = cv_nom;
      int best_h = -1;
      bool keep[MAXP];
      for (int i = 0; i < np_; ++i) keep[i] = use_[i];
      for (int h = 0; h < nh_; ++h) {
        int left = 0, gone = 0;
        for (int i = 0; i < np_; ++i) {
          if (!keep[i]) continue;
          if (touches_(i, h)) ++gone; else ++left;
        }
        if (!gone || left < cfg.min_pairs) continue;
        for (int i = 0; i < np_; ++i) use_[i] = keep[i] && !touches_(i, h);
        double f2 = 0, n2 = 0;
        if (crossValidate_(cfg, f2, n2) && f2 < best_fit) { best_fit = f2; best_nom = n2; best_h = h; }
      }
      if (best_h < 0 || !(best_fit <= cfg.drop_gain * cv_fit)) {
        for (int i = 0; i < np_; ++i) use_[i] = keep[i];
        break;
      }
      for (int i = 0; i < np_; ++i) use_[i] = keep[i] && !touches_(i, best_h);
      cv_fit = best_fit; cv_nom = best_nom;
      ++out.dropped;
    }
    int used = 0;
    for (int i = 0; i < np_; ++i) used += use_[i];
    out.pairs = used;

    H h_all;
    X g_all, x;
    normal_(h_all, g_all);
    if (!solve_(h_all, g_all, cfg, -1, true, x)) return finish_(out, GyroAxisReason::SOLVE);
    Eigen::Matrix3f S;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) S(i, j) = float((i == j ? 1.0 : 0.0) + x[3 * i + j]);
    out.S = S;
    out.offset = x.template tail<3>().template cast<float>();
    out.cv_fit = float(cv_fit);
    out.cv_nominal = float(cv_nom);
    double ss = 0;
    for (int i = 0; i < np_; ++i)
      if (use_[i]) ss += residual_(i, x).squaredNorm();
    out.rms_fit = float(std::sqrt(ss / std::max(1, used)));
    if (!S.allFinite()) return finish_(out, GyroAxisReason::SOLVE);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double lim = i == j ? cfg.max_scale_dev : cfg.max_cross;
        if (std::fabs(double(S(i, j)) - (i == j ? 1.0 : 0.0)) > lim) return finish_(out, GyroAxisReason::IMPLAUSIBLE);
      }
    if (!(cv_fit <= (1.0 - cfg.min_improvement) * cv_nom) || !(cv_fit <= cfg.max_cv))
      return finish_(out, GyroAxisReason::NOT_BETTER);
    out.ok = true;
    return finish_(out, GyroAxisReason::OK);
  }

private:
  static bool finish_(GyroAxisResult& out, GyroAxisReason r) {
    out.reason = r;
    if (!out.ok) { out.S.setIdentity(); }
    return out.ok;
  }

  bool touches_(int i, int h) const { return hp_[i] == h || hq_[i] == h; }
  V3 residual_(int i, const X& x) const {
    return r0_[i].template cast<double>() + J_[i].template cast<double>() * x;
  }

  // Calibrated gravity direction of a hold in its anchor frame and its
  // sensitivity to the parameters.
  static bool direction_(const GyroPoseSums& s, float tempC, const AccelCalibration<float>& acc, V3& d, J& Jd) {
    if (!s.ok()) return false;
    const double n = double(s.n);
    const V3 b = acc.biasT.bias(tempC).template cast<double>();
    const M3 Rs = s.Rs.template cast<double>();
    const M3 A = acc.S.template cast<double>();
    // sum R(A(a - b)) ~ sum R(a - b) + mean(R)(A - I) sum(a - b): the hold's
    // attitude spread times the small (A - I) is negligible.
    const V3 G = s.Ra.template cast<double>() - Rs * b +
                 (Rs / n) * (A - M3::Identity()) * (s.a.template cast<double>() - n * b);
    const double gn = G.norm();
    if (!(gn > 1e-6) || !std::isfinite(gn)) return false;
    d = G / gn;
    M3 dx;
    dx << 0, -d.z(), d.y(), d.z(), 0, -d.x(), -d.y(), d.x(), 0;
    Jd = -dx * (s.P.template cast<double>() / n);
    return Jd.allFinite();
  }

  // Normal equations of the pairs in use.
  void normal_(H& h, X& g) const {
    h.setZero();
    g.setZero();
    for (int i = 0; i < np_; ++i) {
      if (!use_[i]) continue;
      const J j = J_[i].template cast<double>();
      h.noalias() += j.transpose() * j;
      g.noalias() -= j.transpose() * r0_[i].template cast<double>();
    }
  }

  // Ridge solution without the pairs of hold `skip` (removed from the full
  // normal equations); nominal pins the matrix terms at zero, leaving only
  // the rate offset free.
  bool solve_(const H& h_all, const X& g_all, const GyroAxisCfg& cfg, int skip, bool full, X& x) const {
    H h = h_all;
    X g = g_all;
    for (int i = 0; i < np_; ++i) {
      if (!use_[i] || !touches_(i, skip)) continue;
      const J j = J_[i].template cast<double>();
      h.noalias() -= j.transpose() * j;
      g.noalias() += j.transpose() * r0_[i].template cast<double>();
    }
    const double s2 = cfg.noise_rad * cfg.noise_rad;
    for (int k = 0; k < kGyroAxisParams; ++k) {
      const bool matrix = k < 9;
      if (matrix && !full) {
        h.row(k).setZero(); h.col(k).setZero(); h(k, k) = 1; g[k] = 0;
        continue;
      }
      const double prior = matrix ? cfg.prior_matrix : cfg.prior_offset;
      h(k, k) += s2 / (prior * prior);
    }
    Eigen::LDLT<H> ldlt(h);
    if (ldlt.info() != Eigen::Success) return false;
    x = ldlt.solve(g);
    return x.allFinite();
  }

  // Leave-one-hold-out: every hold's pairs predicted from the others.
  bool crossValidate_(const GyroAxisCfg& cfg, double& cv_fit, double& cv_nom) {
    H h_all;
    X g_all;
    normal_(h_all, g_all);
    double sf = 0, sn = 0;
    int m = 0;
    for (int h = 0; h < nh_; ++h) {
      double ef = 0, en = 0;
      int k = 0;
      X xf, xn;
      bool solved = false;
      for (int i = 0; i < np_; ++i) {
        if (!use_[i] || !touches_(i, h)) continue;
        if (!solved) {
          if (!solve_(h_all, g_all, cfg, h, true, xf) || !solve_(h_all, g_all, cfg, h, false, xn)) return false;
          solved = true;
        }
        ef += residual_(i, xf).squaredNorm();
        en += residual_(i, xn).squaredNorm();
        ++k;
      }
      if (!k) continue;
      sf += ef / k; sn += en / k; ++m;
    }
    if (m == 0) return false;
    cv_fit = std::sqrt(sf / m);
    cv_nom = std::sqrt(sn / m);
    return std::isfinite(cv_fit) && std::isfinite(cv_nom);
  }

  int np_ = 0, nh_ = 0;
  Eigen::Vector3f r0_[MAXP];
  JF J_[MAXP];
  int hp_[MAXP] = {}, hq_[MAXP] = {};
  bool use_[MAXP] = {};
};

}  // namespace imu_cal
