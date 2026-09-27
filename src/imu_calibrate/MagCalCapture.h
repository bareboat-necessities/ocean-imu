#pragma once

/* Copyright 2026, Mikhail Grushinskiy

   Host-testable magnetometer capture. Freshness belongs to the incoming
   sensor stream, not the retained sample count. A fixed-seed reservoir of
   temporally spaced observations represents the entire capture, including
   motion after the existing 400-entry buffer fills. No second sample buffer.
*/
#include "CalibrateIMU.h"

namespace imu_cal {

struct MagCaptureCfg {
  uint32_t spacing_ms = 80, min_time_ms = 45000, timeout_ms = 220000, stuck_ms = 12000;
  int required_samples = 360;
  float min_delta_uT = 0.03f;
  float span_min_frac = 0.35f, span_mid_frac = 0.55f, urange_target = 1.05f;
  float min_cov_det = 2e-4f;
};
enum class MagCaptureStatus : uint8_t { CAPTURING, READY, STALE, COVERAGE_LOW, TIMEOUT, BAD_CONFIG };

template <typename T, int N>
class MagCapture {
public:
  using Vec3 = Eigen::Matrix<T,3,1>;
  using Mat3 = Eigen::Matrix<T,3,3>;
  struct Coverage {
    Vec3 span = Vec3::Zero(), urange = Vec3::Zero();
    T min_ratio = 0, mid_ratio = 0, determinant = 0;
    bool ok = false;
  };

  explicit MagCapture(MagCalibrator<T,N>& cal, const MagCaptureCfg& cfg = {}) : cal_(cal), cfg_(cfg) {}

  void begin(uint32_t now) {
    cal_.clear();
    start_ = last_change_ = last_observation_ = now;
    seen_ = 0; random_ = 0x6d2b79f5u;
    have_change_ = have_observation_ = false;
  }

  // Call even when the driver returns no measurement (m == nullptr). A fresh
  // but frozen register is not movement; changes smaller than the threshold
  // accumulate relative to the last meaningful change, not the last poll.
  MagCaptureStatus update(uint32_t now, const Vec3* m) {
    if (cfg_.required_samples > N || cfg_.required_samples < 20 || !cfg_.spacing_ms ||
        cfg_.min_time_ms >= cfg_.timeout_ms) return MagCaptureStatus::BAD_CONFIG;
    if (uint32_t(now - last_change_) > cfg_.stuck_ms) return MagCaptureStatus::STALE;
    if (uint32_t(now - start_) >= cfg_.timeout_ms) return MagCaptureStatus::TIMEOUT;
    if (m && isfinite3(*m)) {
      const T norm = m->norm();
      if (norm >= cal_.min_norm_uT && norm <= cal_.max_norm_uT) {
        if (!have_change_ || (*m - last_changed_).norm() >= T(cfg_.min_delta_uT)) {
          last_changed_ = *m; last_change_ = now; have_change_ = true;
        }
        if ((!have_observation_ || uint32_t(now - last_observation_) >= cfg_.spacing_ms) &&
            (!have_observation_ || (*m - last_observed_).norm() >= T(cfg_.min_delta_uT))) {
          last_observation_ = now; last_observed_ = *m; have_observation_ = true;
          ++seen_;
          if (cal_.buf.n < N) {
            cal_.addSample(*m);
          } else {
            // Algorithm R with a reproducible local PRNG. Every eligible
            // observation has the same retention probability N/seen_. The
            // multiply-high map avoids low-bit patterns in the slot choice.
            random_ ^= random_ << 13; random_ ^= random_ >> 17; random_ ^= random_ << 5;
            const uint32_t slot = uint32_t((uint64_t(random_) * seen_) >> 32);
            if (slot < uint32_t(N)) { cal_.buf.v[slot] = *m; cal_.buf.tempC[slot] = T(0); }
          }
        }
      }
    }
    if (cal_.buf.n >= cfg_.required_samples && uint32_t(now - start_) >= cfg_.min_time_ms)
      return coverage().ok ? MagCaptureStatus::READY : MagCaptureStatus::COVERAGE_LOW;
    return MagCaptureStatus::CAPTURING;
  }

  uint32_t observations() const { return seen_; }

  Coverage coverage() const {
    Coverage out;
    const int n = cal_.buf.n;
    if (n < 20) return out;
    Vec3 lo = cal_.buf.v[0], hi = lo, mean = Vec3::Zero();
    for (int i = 0; i < n; ++i) {
      const Vec3& v = cal_.buf.v[i];
      if (!isfinite3(v)) return out;
      lo = lo.cwiseMin(v); hi = hi.cwiseMax(v); mean += v;
    }
    mean /= T(n);
    const Vec3 center = T(0.5) * (lo + hi);
    Vec3 ulo = Vec3::Constant(T(1)), uhi = Vec3::Constant(T(-1));
    Mat3 cov = Mat3::Zero(); int nc = 0;
    for (int i = 0; i < n; ++i) {
      const Vec3 centered = cal_.buf.v[i] - center;
      const T radius = centered.norm();
      if (radius > T(1e-6)) {
        const Vec3 u = centered / radius;
        ulo = ulo.cwiseMin(u); uhi = uhi.cwiseMax(u);
      }
      const Vec3 d = cal_.buf.v[i] - mean;
      const T r = d.norm();
      if (r > T(1e-6)) { const Vec3 u = d/r; cov.noalias() += u*u.transpose(); ++nc; }
    }
    out.span = hi - lo; out.urange = uhi - ulo;
    T spans[3] = {out.span[0], out.span[1], out.span[2]}; sort_small(spans, 3);
    if (spans[2] < T(1e-3) || nc < 20) return out;
    out.min_ratio = spans[0] / spans[2]; out.mid_ratio = spans[1] / spans[2];
    cov /= T(nc); out.determinant = cov.determinant();
    int axes = 0;
    for (int j = 0; j < 3; ++j) axes += out.urange[j] >= T(cfg_.urange_target) ? 1 : 0;
    out.ok = out.min_ratio >= T(cfg_.span_min_frac) && out.mid_ratio >= T(cfg_.span_mid_frac) &&
             axes >= 2 && finiteT(out.determinant) && out.determinant >= T(cfg_.min_cov_det);
    return out;
  }

  float progress(uint32_t now) const {
    const Coverage c = coverage();
    T ranges[3] = {c.urange[0], c.urange[1], c.urange[2]}; sort_small(ranges, 3);
    const float samples = float(cal_.buf.n) / cfg_.required_samples;
    const float time = float(uint32_t(now - start_)) / cfg_.min_time_ms;
    const float cover = float(ranges[1]) / cfg_.urange_target;
    return clamp<float>(std::fmin(samples, std::fmin(time, cover)), 0.f, 1.f);
  }

private:
  MagCalibrator<T,N>& cal_;
  MagCaptureCfg cfg_;
  uint32_t start_ = 0, last_change_ = 0, last_observation_ = 0, seen_ = 0, random_ = 0;
  Vec3 last_changed_ = Vec3::Zero(), last_observed_ = Vec3::Zero();
  bool have_change_ = false, have_observation_ = false;
};
} // namespace imu_cal
