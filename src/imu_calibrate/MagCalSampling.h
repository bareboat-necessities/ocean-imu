#pragma once
// Copyright 2026, Mikhail Grushinskiy
#include "MagCalCapture.h"

namespace imu_cal {

// Non-overlapping magnetic moments reduce sensor noise without requiring slow
// turns. The covariance preserves E[|A(x-b)|^2] through arbitrary motion within
// a window; averaging only x would shrink the fitted ellipsoid during turns.
// Count distinct readings, never repeated polls.
class MagSampleWindow {
 public:
  using V=Eigen::Vector3f;
  using M=Eigen::Matrix3f;
  static constexpr uint32_t window_ms=400;
  static constexpr int min_samples=12;
  static MagCaptureCfg captureCfg(bool verify=false) {
    MagCaptureCfg c;
    c.spacing_ms=window_ms;c.required_samples=verify?140:160;
    c.min_time_ms=verify?56000:45000;c.timeout_ms=verify?180000:220000;
    return c;
  }
  void begin() { clear();have_words_=false; }
  bool update(uint32_t now,const V* m,V& mean,M& covariance) {
    if (!m || !m->allFinite() || m->norm()<5 || m->norm()>200) {clear();return false;}
    // Display/I2C pauses are not motion failures. Only a genuinely stale
    // magnetic stream discards an incomplete window.
    if(n_ && uint32_t(now-last_fresh_)>500)clear();
    if(!have_words_ || *m!=last_) {
      last_=*m;have_words_=true;last_fresh_=now;
      if(!n_)start_=now;
      const Eigen::Vector3d delta=m->cast<double>()-mean_;
      mean_+=delta/++n_;
      scatter_.noalias()+=delta*(m->cast<double>()-mean_).transpose();
    }
    // Slow ODR or a missed poll extends the window until it has enough fresh
    // data, rather than repeatedly throwing away five or six useful readings.
    if(n_<min_samples || uint32_t(now-start_)<window_ms ||
        uint32_t(now-last_fresh_)>100)return false;
    mean=mean_.cast<float>();
    covariance=((scatter_+scatter_.transpose())/(2*n_)).cast<float>();
    clear();return true;
  }
 private:
  void clear() {n_=0;mean_.setZero();scatter_.setZero();}
  int n_=0;
  uint32_t start_=0,last_fresh_=0;
  bool have_words_=false;
  V last_=V::Zero();
  Eigen::Vector3d mean_=Eigen::Vector3d::Zero();
  Eigen::Matrix3d scatter_=Eigen::Matrix3d::Zero();
};
} // namespace imu_cal
