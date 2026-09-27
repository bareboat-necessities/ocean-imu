#pragma once
// Copyright 2026, Mikhail Grushinskiy
#include "MagCalCapture.h"

namespace imu_cal {

// Non-overlapping magnetic means reduce sensor noise without depending on a
// vendor register preset. Count distinct readings, never repeated polls. The
// independently measured gyro bounds rotation during a mean (17 degrees max).
class MagSampleWindow {
 public:
  using V=Eigen::Vector3f;
  static constexpr uint32_t window_ms=400;
  static constexpr int min_samples=8;
  static MagCaptureCfg captureCfg(bool verify=false) {
    MagCaptureCfg c;
    c.spacing_ms=window_ms;c.required_samples=verify?140:160;
    c.min_time_ms=verify?56000:45000;c.timeout_ms=verify?120000:220000;
    return c;
  }
  void begin() { clear();have_words_=have_time_=false;too_fast_=false; }
  bool update(uint32_t now,const V* m,const V* rate,V& mean) {
    too_fast_=false;
    if (!m || !rate || !m->allFinite() || !rate->allFinite() ||
        m->norm()<5 || m->norm()>200) {clear();have_time_=false;return false;}
    const uint32_t dt=now-last_time_;
    if(have_time_ && dt>60)clear();
    if(have_time_ && n_ && dt<=60)angle_+=rate->norm()*float(dt)*.001f;
    last_time_=now;have_time_=true;
    if(angle_>.30f) {clear();too_fast_=true;return false;}
    if(!have_words_ || *m!=last_) {
      last_=*m;have_words_=true;last_fresh_=now;
      if(!n_)start_=now;
      sum_+=m->cast<double>();++n_;
    }
    if(!n_ || uint32_t(now-start_)<window_ms)return false;
    const bool ready=n_>=min_samples && uint32_t(now-last_fresh_)<=100;
    if(ready)mean=(sum_/n_).cast<float>();
    clear();return ready;
  }
  bool tooFast() const { return too_fast_; }
 private:
  void clear() {n_=0;sum_.setZero();angle_=0;}
  int n_=0;
  uint32_t start_=0,last_fresh_=0,last_time_=0;
  double angle_=0;
  bool have_words_=false,have_time_=false,too_fast_=false;
  V last_=V::Zero();
  Eigen::Vector3d sum_=Eigen::Vector3d::Zero();
};
} // namespace imu_cal
