#pragma once
/* Copyright 2026, Mikhail Grushinskiy

   Gyro bias capture from a continuous quiet interval. Sample magnitude is only
   a sanity check: block scatter, changes in block means, gravity direction and
   (when available) magnetic changes qualify stillness. A disturbed interval is
   discarded, never averaged into the bias. Constant yaw without a magnetic or
   other independent reference is intrinsically indistinguishable from bias.
*/
#include "CalibrateIMU.h"

namespace imu_cal {
struct GyroCaptureCfg {
  uint32_t block_ms=250, useful_ms=6000, timeout_ms=70000, stuck_ms=12000;
  uint32_t spacing_ms=25, max_gap_ms=60, mag_stale_ms=1000;
  double max_gyro_std=0.005, max_gyro_change=0.004;
  double max_accel_std=0.08, max_accel_change=0.06;
  double max_mag_change=0.6; // uT, with an additional 1.5% raw-field allowance
};
enum class GyroCaptureStatus : uint8_t { SETTLING, HOLDING, MOVING, READY, STALE, TIMEOUT };

template<typename T,int N,int BINS=8>
class GyroCapture {
  using V=Eigen::Matrix<T,3,1>;
  using D=Eigen::Vector3d;
public:
  explicit GyroCapture(GyroCalibrator<T,N,BINS>& cal,const GyroCaptureCfg& cfg={}) : cal_(cal),cfg_(cfg) {}
  void begin(uint32_t now) {
    cal_.clear();start_=last_fresh_=last_input_=last_kept_=now;
    have_input_=have_words_=have_mag_words_=false;quiet_=kept_blocks_=resets_=0;reference_=false;
    status_=GyroCaptureStatus::SETTLING;resetBlock_();
  }
  GyroCaptureStatus update(uint32_t now,const V* a,const V* w,T tempC,const V* mag=nullptr) {
    if(status_==GyroCaptureStatus::READY) return status_;
    if(uint32_t(now-start_)>=cfg_.timeout_ms) return status_=GyroCaptureStatus::TIMEOUT;
    if(uint32_t(now-last_fresh_)>cfg_.stuck_ms) return status_=GyroCaptureStatus::STALE;
    if(!a || !w) {
      if(have_input_ && uint32_t(now-last_input_)>cfg_.max_gap_ms) disturb_();
      return status_;
    }
    if(!a->allFinite() || !w->allFinite() || !finiteT(tempC) ||
       w->norm()>cal_.max_gyro_norm || std::fabs(a->norm()-cal_.g)>cal_.max_accel_dev) {
      disturb_();return status_;
    }
    // M5's update bit alone cannot establish a fresh register value.
    if(have_words_ && *a==last_a_ && *w==last_w_) {
      if(uint32_t(now-last_fresh_)>cfg_.max_gap_ms) disturb_();
      return status_;
    }
    last_a_=*a;last_w_=*w;have_words_=true;
    if(have_input_ && (uint32_t(now-last_input_)==0 || uint32_t(now-last_input_)>0x7fffffffu)) return status_;
    if(have_input_ && uint32_t(now-last_input_)>cfg_.max_gap_ms) disturb_();
    last_input_=last_fresh_=now;have_input_=true;
    if(n_==0) {block_start_=now;ar_=a->template cast<double>();wr_=w->template cast<double>();}
    const D da=a->template cast<double>()-ar_,dw=w->template cast<double>()-wr_;
    sa_+=da;aa_+=da.cwiseProduct(da);sw_+=dw;ww_+=dw.cwiseProduct(dw);++n_;
    if(mag && mag->allFinite() && mag->norm()>=T(5) && mag->norm()<=T(200)) {
      if(!have_mag_words_ || *mag!=last_mag_) last_mag_change_=now;
      last_mag_=*mag;have_mag_words_=true;
      sm_+=mag->template cast<double>();++nm_;
    }
    // Retain at most 40 Hz, so six seconds fit the existing 400 entries.
    if(nb_<kBlockSamples && (nb_==0 || uint32_t(now-last_kept_)>=cfg_.spacing_ms)) {
      block_w_[nb_]=*w;block_a_[nb_]=*a;block_t_[nb_]=tempC;++nb_;last_kept_=now;
    }
    if(uint32_t(now-block_start_)<cfg_.block_ms) return status_;
    const D ma=ar_+sa_/n_,mw=wr_+sw_/n_;
    const double av=(aa_/n_-(sa_/n_).cwiseProduct(sa_/n_)).cwiseMax(0.0).sum();
    const double wv=(ww_/n_-(sw_/n_).cwiseProduct(sw_/n_)).cwiseMax(0.0).sum();
    const bool have_mag=nm_>=3;
    const D mm=have_mag ? D(sm_/nm_) : D::Zero();
    bool moving=n_<8 || nb_<6 || std::sqrt(av)>cfg_.max_accel_std || std::sqrt(wv)>cfg_.max_gyro_std;
    // Once a magnetic stream is present, freezing or losing it cannot turn
    // slow yaw into an apparently quiet interval. Reset preserves this fact.
    if(have_mag_words_) moving=moving || !have_mag || uint32_t(now-last_mag_change_)>cfg_.mag_stale_ms;
    if(reference_) {
      moving=moving || (ma-ref_a_).norm()>cfg_.max_accel_change || (mw-ref_w_).norm()>cfg_.max_gyro_change;
      // Losing a reference mid-hold must not silently remove the yaw check.
      if(ref_mag_) moving=moving || !have_mag || (mm-ref_m_).norm()>std::max(cfg_.max_mag_change,0.015*ref_m_.norm());
      // A newly available yaw reference starts a complete qualified hold.
      if(!ref_mag_ && have_mag) moving=true;
    }
    if(moving) {disturb_();return status_;}
    if(!reference_) {ref_a_=ma;ref_w_=mw;ref_m_=mm;ref_mag_=have_mag;reference_=true;}
    ++quiet_;
    if(quiet_<=3) {status_=GyroCaptureStatus::SETTLING;resetBlock_();return status_;}
    for(int i=0;i<nb_;++i) if(!cal_.addSample(block_w_[i],block_a_[i],block_t_[i])) {
      disturb_();return status_;
    }
    ++kept_blocks_;
    status_=(kept_blocks_*cfg_.block_ms>=cfg_.useful_ms && cal_.buf.n>=220) ?
             GyroCaptureStatus::READY : GyroCaptureStatus::HOLDING;
    resetBlock_();return status_;
  }
  float progress() const {return std::min(1.0f,float(kept_blocks_*cfg_.block_ms)/cfg_.useful_ms);}
  uint32_t resets() const {return resets_;}
private:
  void resetBlock_() {n_=nm_=nb_=0;sa_.setZero();aa_.setZero();sw_.setZero();ww_.setZero();sm_.setZero();}
  void disturb_() {
    if(quiet_ || kept_blocks_ || n_) ++resets_;
    cal_.clear();quiet_=kept_blocks_=0;reference_=false;status_=GyroCaptureStatus::MOVING;resetBlock_();
  }
  static constexpr int kBlockSamples=16;
  GyroCalibrator<T,N,BINS>& cal_;
  GyroCaptureCfg cfg_;
  GyroCaptureStatus status_=GyroCaptureStatus::SETTLING;
  uint32_t start_=0,last_fresh_=0,last_input_=0,last_kept_=0,block_start_=0,last_mag_change_=0;
  uint32_t quiet_=0,kept_blocks_=0,resets_=0;
  bool have_input_=false,have_words_=false,have_mag_words_=false,reference_=false,ref_mag_=false;
  V last_a_=V::Zero(),last_w_=V::Zero(),last_mag_=V::Zero();
  D ar_,wr_,sa_,aa_,sw_,ww_,sm_,ref_a_,ref_w_,ref_m_;
  int n_=0,nm_=0,nb_=0;
  V block_a_[kBlockSamples],block_w_[kBlockSamples];T block_t_[kBlockSamples];
};
} // namespace imu_cal
