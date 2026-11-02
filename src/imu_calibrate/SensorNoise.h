#pragma once
// Copyright 2026, Mikhail Grushinskiy
// Bounded-memory, non-overlapping-cluster Allan deviation. Adjacent cluster
// differences share a mean; this is not the all-offset overlapping estimator.
#include "imu_calibrate/CalibrateIMU.h"
#include <algorithm>
#include <cstdint>
#include <cstring>

namespace imu_cal {
// Axis order: accel xyz [m/s^2], gyro xyz [rad/s], in runtime vessel frame.
// N^2 is continuous white-noise intensity: Var(sample)=N^2/dt for independent
// samples. One-sided amplitude spectral density is sqrt(2)*N.
struct SensorNoiseModel {
  uint32_t duration_ms=0,samples=0;
  float period_s=0,temperature_span=0,jitter_rms_fraction=0;
  float white_density[6]{},sample_std[6]{},bias_instability[6]{},bias_random_walk[6]{};
  uint8_t valid=0,white_mask=0,bias_mask=0,walk_mask=0;
};
static_assert(sizeof(SensorNoiseModel)==120,"Noise model layout");

enum class NoiseCaptureStatus:uint8_t {CAPTURING,READY,INVALID_SAMPLE,TIMING,FROZEN,MOVING,TEMPERATURE,INSUFFICIENT};
inline bool noiseModelValid(const SensorNoiseModel& n) {
  if(!n.valid)return n.valid==0;
  if(n.valid!=1 || n.white_mask!=63 || (n.bias_mask&~63) || (n.walk_mask&~63) ||
     n.duration_ms<2400000 || n.samples<10000 || !std::isfinite(n.period_s) || n.period_s<=0 ||
     !std::isfinite(n.temperature_span) || n.temperature_span>2 ||
     !std::isfinite(n.jitter_rms_fraction) || n.jitter_rms_fraction>.10f)return false;
  for(int k=0;k<6;++k) {
    if(!std::isfinite(n.white_density[k]) || n.white_density[k]<=0 ||
       !std::isfinite(n.sample_std[k]) || n.sample_std[k]<=0 ||
       !std::isfinite(n.bias_instability[k]) || n.bias_instability[k]<0 ||
       !std::isfinite(n.bias_random_walk[k]) || n.bias_random_walk[k]<0)return false;
    if(!(n.bias_mask&(1<<k)) && n.bias_instability[k]!=0)return false;
    if(!(n.walk_mask&(1<<k)) && n.bias_random_walk[k]!=0)return false;
  }
  return true;
}

template<int LEVELS=17>
class StationaryNoiseCapture {
public:
  static_assert(LEVELS>=12 && LEVELS<=20,"Allan levels");
  using V6=Eigen::Matrix<double,6,1>;
  void begin(uint32_t us,double period=.005,double gravity=9.8025605) {
    *this=StationaryNoiseCapture{};last_us_=us;period_=period;g_=gravity;
  }
  NoiseCaptureStatus status() const {return status_;}
  float progress() const {return float(std::min(1.0,elapsed_/2400.0));}
  uint32_t samples() const {return count_;}
  uint32_t pairs(int level) const {return level>=0&&level<LEVELS?levels_[level].pairs:0;}
  double tau(int level) const {return std::ldexp(count_>1?elapsed_/(count_-1):period_,level);}
  double deviation(int level,int axis) const {
    if(!pairs(level) || axis<0 || axis>=6)return NAN;
    return std::sqrt(levels_[level].sse[axis]/(2*levels_[level].pairs));
  }
  NoiseCaptureStatus add(uint32_t us,const Eigen::Vector3f& a,const Eigen::Vector3f& w,float temp,
                         const Eigen::Vector3f* mag=nullptr) {
    if(status_!=NoiseCaptureStatus::CAPTURING)return status_;
    if(!a.allFinite() || !w.allFinite() || !std::isfinite(temp))return status_=NoiseCaptureStatus::INVALID_SAMPLE;
    const uint32_t delta=us-last_us_;
    if(count_ && (delta<period_*1e6*.5 || delta>period_*1e6*1.5))return status_=NoiseCaptureStatus::TIMING;
    if(count_) {const double dt=delta*1e-6;elapsed_+=dt;jitter_+=std::pow(dt/period_-1,2);}
    last_us_=us;
    V6 x;x<<a.cast<double>(),w.cast<double>();
    if(count_ && x==last_words_) {if(++frozen_>unsigned(1/period_))return status_=NoiseCaptureStatus::FROZEN;}
    else frozen_=0;
    last_words_=x;
    if(mag && mag->allFinite() && mag->norm()>5) {
      if(!have_mag_words_ || *mag!=last_mag_) {mag_sum_+=mag->template cast<double>();++mag_n_;}
      last_mag_=*mag;have_mag_words_=true;
    }
    if(!count_) {tlo_=thi_=temp;reference_=x;}
    tlo_=std::min(tlo_,double(temp));thi_=std::max(thi_,double(temp));
    if(thi_-tlo_>2)return status_=NoiseCaptureStatus::TEMPERATURE;
    if(std::fabs(a.norm()-g_)>.3 || w.norm()>.05)return status_=NoiseCaptureStatus::MOVING;
    ++count_;
    const V6 d=x-mean_;mean_+=d/double(count_);m2_+=d.cwiseProduct(x-mean_);
    // One-second block means detect orientation changes/slow turns in addition
    // to instantaneous bounds. A sensor-only test cannot disprove all motion.
    block_+=x;++block_n_;
    if(block_n_>=unsigned(1/period_)) {
      const V6 b=block_/block_n_;
      if(have_ref_ && ((b.template head<3>()-reference_.template head<3>()).norm()>.06 ||
         (b.template tail<3>()-reference_.template tail<3>()).norm()>.005))
        return status_=NoiseCaptureStatus::MOVING;
      if(!have_ref_){reference_=b;have_ref_=true;}
      // Fixed magnetic reference also catches yaw that leaves gravity
      // unchanged. Once present, a stale/lost stream invalidates the run.
      if(have_mag_words_) {
        if(mag_n_<8)return status_=NoiseCaptureStatus::INVALID_SAMPLE;
        const Eigen::Vector3d m=mag_sum_/mag_n_;
        if(have_mag_ref_ && (m-mag_ref_).norm()>std::max(1.5,.03*mag_ref_.norm()))
          return status_=NoiseCaptureStatus::MOVING;
        if(!have_mag_ref_){mag_ref_=m;have_mag_ref_=true;}
        mag_sum_.setZero();mag_n_=0;
      }
      block_.setZero();block_n_=0;
    }
    V6 incoming=x;
    for(int k=0;k<LEVELS;++k) {
      Level& q=levels_[k];
      // Each level receives the mean of two clusters from its predecessor.
      if(q.have_last){q.sse+=(incoming-q.last).cwiseAbs2();++q.pairs;}
      q.last=incoming;q.have_last=true;q.sum+=incoming;
      if(++q.pending<2)break;
      incoming=q.sum*.5;q.pending=0;q.sum.setZero();
    }
    if(elapsed_>=2400)status_=NoiseCaptureStatus::READY;
    return status_;
  }
  bool report(SensorNoiseModel& out) const {
    out={};
    if(status_!=NoiseCaptureStatus::READY || count_<10000 ||
       std::sqrt(jitter_/std::max<uint32_t>(1,count_-1))>.10)return false;
    out.duration_ms=uint32_t(elapsed_*1000);out.samples=count_;out.period_s=float(elapsed_/(count_-1));
    out.temperature_span=float(thi_-tlo_);out.jitter_rms_fraction=float(std::sqrt(jitter_/(count_-1)));
    for(int axis=0;axis<6;++axis) {
      out.sample_std[axis]=float(std::sqrt(m2_[axis]/(count_-1)));
      // Require FOUR dyadic levels (~0.9 decade), >=32 adjacent differences,
      // and the expected slope. Never convert an arbitrary Allan minimum into
      // an identified long-term model. Missing terms remain explicitly zero.
      for(int kind=0;kind<3;++kind)for(int l=0;l+3<LEVELS;++l) {
        if(pairs(l+3)<32 || (kind && tau(l)<1))continue;
        double sx=0,sy=0,sxx=0,sxy=0;double ly[4],lx[4];
        bool finite=true;
        for(int j=0;j<4;++j) {lx[j]=std::log(tau(l+j));ly[j]=std::log(deviation(l+j,axis));
          finite=finite && std::isfinite(ly[j]);sx+=lx[j];sy+=ly[j];sxx+=lx[j]*lx[j];sxy+=lx[j]*ly[j];}
        if(!finite)continue;
        const double slope=(4*sxy-sx*sy)/(4*sxx-sx*sx);
        const double target=kind==0?-.5:kind==1?0:.5;
        double residual=0;const double intercept=(sy-slope*sx)/4;
        for(int j=0;j<4;++j)residual+=std::pow(ly[j]-intercept-slope*lx[j],2);
        if(std::fabs(slope-target)>.10 || std::sqrt(residual/4)>.08)continue;
        double coefficient=0;
        for(int j=0;j<4;++j) {
          const double t=tau(l+j),s=deviation(l+j,axis);
          coefficient+=kind==0?s*std::sqrt(t):kind==1?s/.664:s*std::sqrt(3/t);
        }
        coefficient/=4;
        if(kind==0){out.white_density[axis]=float(coefficient);out.white_mask|=1<<axis;}
        if(kind==1){out.bias_instability[axis]=float(coefficient);out.bias_mask|=1<<axis;}
        if(kind==2){out.bias_random_walk[axis]=float(coefficient);out.walk_mask|=1<<axis;}
        break;
      }
    }
    out.valid=out.white_mask==63;
    return out.valid && noiseModelValid(out);
  }
private:
  struct Level {V6 sum=V6::Zero(),last=V6::Zero(),sse=V6::Zero();
    uint32_t pairs=0;unsigned pending=0;bool have_last=false;};
  Level levels_[LEVELS];
  NoiseCaptureStatus status_=NoiseCaptureStatus::CAPTURING;
  uint32_t last_us_=0,count_=0,frozen_=0,block_n_=0;
  double period_=.005,g_=9.8025605,elapsed_=0,jitter_=0,tlo_=0,thi_=0;
  V6 mean_=V6::Zero(),m2_=V6::Zero(),block_=V6::Zero(),reference_=V6::Zero(),last_words_=V6::Zero();
  bool have_ref_=false;
  bool have_mag_words_=false,have_mag_ref_=false;unsigned mag_n_=0;
  Eigen::Vector3f last_mag_=Eigen::Vector3f::Zero();
  Eigen::Vector3d mag_sum_=Eigen::Vector3d::Zero(),mag_ref_=Eigen::Vector3d::Zero();
};
} // namespace imu_cal
