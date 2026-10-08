#pragma once

/* Copyright 2026, Mikhail Grushinskiy

   Host-testable magnetometer capture. Freshness belongs to the incoming
   sensor stream, not the retained sample count. Direction-balanced retention
   keeps late useful rotations after the fixed buffer fills. Guidance names
   broad movements; no exact poses or angles are required.
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
    int cells = 0;
    bool ok = false;
  };
  // Verification only: fresh raw readings scored with the frozen candidate.
  struct RawResidual {
    int n = 0, tail = 0;
    double rms = 0;
    bool ok = false;
  };

  explicit MagCapture(MagCalibrator<T,N>& cal, const MagCaptureCfg& cfg = {},
                      const MagCalibration<T>* fixed = nullptr) : cal_(cal), cfg_(cfg), fixed_(fixed) {}

  void begin(uint32_t now) {
    cal_.clear();
    start_ = last_change_ = last_observation_ = now;
    seen_ = 0; random_ = 0x6d2b79f5u;
    have_change_ = have_observation_ = false;
    raw_n_ = raw_tail_ = 0; raw_sse_ = 0;
    coverage_dirty_ = true;
    hint_ = "Turn slowly"; hint_ms_ = now;
  }

  // Call even when the driver returns no measurement (m == nullptr). A fresh
  // but frozen register is not movement; changes smaller than the threshold
  // accumulate relative to the last meaningful change, not the last poll.
  MagCaptureStatus update(uint32_t now, const Vec3* m, const Vec3* raw = nullptr,
                          const Mat3* covariance = nullptr) {
    if (cfg_.required_samples > N || cfg_.required_samples < 20 || !cfg_.spacing_ms ||
        cfg_.min_time_ms >= cfg_.timeout_ms) return MagCaptureStatus::BAD_CONFIG;
    if (uint32_t(now - last_change_) > cfg_.stuck_ms) return MagCaptureStatus::STALE;
    if (uint32_t(now - start_) >= cfg_.timeout_ms) return MagCaptureStatus::TIMEOUT;
    // Raw-stream freshness is independent of whether a motion-qualified mean
    // was emitted. Turning too fast calls for guidance, not "no sensor data".
    const Vec3* fresh = raw ? raw : m;
    if (fresh && isfinite3(*fresh) && fresh->norm() >= cal_.min_norm_uT && fresh->norm() <= cal_.max_norm_uT &&
        (!have_change_ || (*fresh-last_changed_).norm() >= T(cfg_.min_delta_uT))) {
      last_changed_=*fresh;last_change_=now;have_change_=true;
      if (fixed_ && raw) addRaw_(*raw);  // window means are not single readings
    }
    if (m && isfinite3(*m)) {
      const T norm = m->norm();
      if (norm >= cal_.min_norm_uT && norm <= cal_.max_norm_uT) {
        if ((!have_observation_ || uint32_t(now - last_observation_) >= cfg_.spacing_ms) &&
            (!have_observation_ || (*m - last_observed_).norm() >= T(cfg_.min_delta_uT))) {
          last_observation_ = now; last_observed_ = *m; have_observation_ = true;
          ++seen_;
          if (cal_.buf.n < N) {
            coverage_dirty_ = cal_.addSample(*m, uint32_t(now-start_), covariance) || coverage_dirty_;
          } else {
            random_ ^= random_ << 13; random_ ^= random_ >> 17; random_ ^= random_ << 5;
            const int slot = replacement_(*m);
            if (slot >= 0) {
              cal_.buf.v[slot] = *m; cal_.buf.tempC[slot] = T(0);
              cal_.sample_ms[slot] = uint32_t(now-start_);
              if(covariance)cal_.sample_cov[slot]=*covariance;
              else cal_.sample_cov[slot].setZero();
              coverage_dirty_ = true;
            }
          }
        }
      }
    }
    if (cal_.buf.n >= cfg_.required_samples && uint32_t(now - start_) >= cfg_.min_time_ms)
      return coverage().ok ? MagCaptureStatus::READY : MagCaptureStatus::COVERAGE_LOW;
    return MagCaptureStatus::CAPTURING;
  }

  uint32_t observations() const { return seen_; }

  RawResidual rawResidual() const {
    RawResidual r;
    r.n = raw_n_; r.tail = raw_tail_;
    if (!fixed_ || raw_n_ == 0) return r;
    const double field = double(fixed_->field_uT);
    r.rms = std::sqrt(raw_sse_ / raw_n_);
    r.ok = raw_n_ >= MagFitLimits::min_raw_samples && r.rms <= MagFitLimits::rawRmsLimit(field) &&
           raw_tail_ <= MagFitLimits::max_raw_tail_fraction * raw_n_;
    return r;
  }

  Coverage coverage() const {
    if (!coverage_dirty_) return coverage_;
    coverage_ = computeCoverage_();
    coverage_dirty_ = false;
    return coverage_;
  }

private:
  Coverage computeCoverage_() const {
    Coverage out;
    const int n = cal_.buf.n;
    if (n < 20) return out;
    Vec3 lo = mapped_(cal_.buf.v[0]), hi = lo, mean = Vec3::Zero();
    for (int i = 0; i < n; ++i) {
      const Vec3 v = mapped_(cal_.buf.v[i]);
      if (!isfinite3(v)) return out;
      lo = lo.cwiseMin(v); hi = hi.cwiseMax(v); mean += v;
    }
    mean /= T(n);
    const Vec3 center = fixed_ ? Vec3::Zero() : Vec3(T(0.5) * (lo + hi));
    Vec3 ulo = Vec3::Constant(T(1)), uhi = Vec3::Constant(T(-1));
    Mat3 cov = Mat3::Zero(); int nc = 0;
    int cells[27]{};
    for (int i = 0; i < n; ++i) {
      const Vec3 v = mapped_(cal_.buf.v[i]);
      const Vec3 centered = v - center;
      const T radius = centered.norm();
      if (radius > T(1e-6)) {
        const Vec3 u = centered / radius;
        ulo = ulo.cwiseMin(u); uhi = uhi.cwiseMax(u);
        ++cells[magDirectionCell(u)];
      }
      const Vec3 d = v - mean;
      const T r = d.norm();
      if (r > T(1e-6)) { const Vec3 u = d/r; cov.noalias() += u*u.transpose(); ++nc; }
    }
    out.span = hi - lo; out.urange = uhi - ulo;
    for (int j = 0; j < 27; ++j) out.cells += cells[j] >= 3;
    T spans[3] = {out.span[0], out.span[1], out.span[2]}; sort_small(spans, 3);
    if (spans[2] < T(1e-3) || nc < 20) return out;
    out.min_ratio = spans[0] / spans[2]; out.mid_ratio = spans[1] / spans[2];
    cov /= T(nc); out.determinant = cov.determinant();
    int axes = 0;
    for (int j = 0; j < 3; ++j) axes += out.urange[j] >= T(cfg_.urange_target) ? 1 : 0;
    out.ok = out.min_ratio >= T(cfg_.span_min_frac) && out.mid_ratio >= T(cfg_.span_mid_frac) &&
             axes >= 2 && finiteT(out.determinant) && out.determinant >= T(cfg_.min_cov_det) &&
             out.cells >= MagFitLimits::min_cells;
    return out;
  }

public:
  float progress(uint32_t now) const {
    const Coverage c = coverage();
    T ranges[3] = {c.urange[0], c.urange[1], c.urange[2]}; sort_small(ranges, 3);
    const float samples = float(cal_.buf.n) / cfg_.required_samples;
    const float time = float(uint32_t(now - start_)) / cfg_.min_time_ms;
    const float cover = std::fmin(float(c.determinant)/cfg_.min_cov_det,
        std::fmin(float(ranges[1]) / cfg_.urange_target,
        std::fmin(float(c.cells)/MagFitLimits::min_cells,
        std::fmin(float(c.min_ratio)/cfg_.span_min_frac,float(c.mid_ratio)/cfg_.span_mid_frac))));
    return clamp<float>(std::fmin(samples, std::fmin(time, cover)), 0.f, 1.f);
  }

  // Three-second hysteresis keeps the small screen readable. Weak magnetic
  // axis coverage suggests a broad rotation that changes that component.
  const char* hint(uint32_t now) {
    if (uint32_t(now-hint_ms_) < 3000) return hint_;
    const Coverage c = coverage();
    const char* next = "Turn slowly";
    if (cal_.buf.n >= 20 && !c.ok) {
      int axis = 0;
      c.span.minCoeff(&axis);
      next = axis == 2 ? "Flip over" : (axis == 0 ? "Tilt forward/back" : "Roll left/right");
    }
    hint_ = next; hint_ms_ = now;
    return hint_;
  }

private:
  Vec3 mapped_(const Vec3& m) const { return fixed_ ? fixed_->apply(m) : m; }

  void addRaw_(const Vec3& m) {
    const double field = double(fixed_->field_uT);
    const double e = double(fixed_->apply(m).norm()) - field;
    if (!std::isfinite(e)) return;
    ++raw_n_; raw_sse_ += e * e;
    if (std::fabs(e) > MagFitLimits::rawTailLimit(field)) ++raw_tail_;
  }

  int replacement_(const Vec3& m) const {
    Vec3 lo=mapped_(cal_.buf.v[0]),hi=lo;
    for(int i=1;i<N;++i) { const Vec3 v=mapped_(cal_.buf.v[i]);lo=lo.cwiseMin(v);hi=hi.cwiseMax(v); }
    const Vec3 center = fixed_ ? Vec3::Zero() : Vec3(T(0.5)*(lo+hi));
    int counts[27]{};
    for(int i=0;i<N;++i) ++counts[magDirectionCell(mapped_(cal_.buf.v[i])-center)];
    const int incoming=magDirectionCell(mapped_(m)-center);
    int crowded=0;
    for(int j=1;j<27;++j) if(counts[j]>counts[crowded]) crowded=j;
    if(counts[incoming]+2<counts[crowded]) {
      int pick=int((uint64_t(random_)*uint32_t(counts[crowded]))>>32);
      for(int i=0;i<N;++i) if(magDirectionCell(mapped_(cal_.buf.v[i])-center)==crowded && pick--==0) return i;
    }
    // Preserve temporal replacement even after directional coverage is full.
    const uint32_t slot=uint32_t((uint64_t(random_)*seen_)>>32);
    return slot<uint32_t(N) ? int(slot) : -1;
  }
  MagCalibrator<T,N>& cal_;
  MagCaptureCfg cfg_;
  const MagCalibration<T>* fixed_;
  mutable Coverage coverage_;
  mutable bool coverage_dirty_ = true;
  const char* hint_ = "Turn slowly";
  uint32_t hint_ms_ = 0;
  uint32_t start_ = 0, last_change_ = 0, last_observation_ = 0, seen_ = 0, random_ = 0;
  int raw_n_ = 0, raw_tail_ = 0;
  double raw_sse_ = 0;
  Vec3 last_changed_ = Vec3::Zero(), last_observed_ = Vec3::Zero();
  bool have_change_ = false, have_observation_ = false;
};
} // namespace imu_cal
