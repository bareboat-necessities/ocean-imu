#pragma once
// Copyright 2026, Mikhail Grushinskiy
// Sensor-clock controller. Display/input work must run outside this controller.
#include "imu_calibrate/GyroCalCapture.h"
#include "imu_calibrate/GyroCalFit.h"
#include "imu_calibrate/GyroCalibrationPlan.h"

namespace imu_cal {
enum class GyroTurnPhase : uint8_t { REST, TWIST, TURN, HOLD, DONE, FAILED };

// Four readable faces from the measured accelerometer frame. USB-down cannot
// rest on its cable; screen-down would hide the live instructions.
template<int SAMPLES=131072>
class GyroTurnCapture {
public:
  static constexpr uint32_t HOLD_MS=3500,USEFUL_MS=2000,MAG_WINDOW_MS=1500;
  using V=Eigen::Vector3f;
  using Fit=GyroRotationFit<SAMPLES,24>;
  static GyroCaptureCfg endpointConfig() { return config_(); }
  GyroTurnCapture(Fit& fit,GyroCalibrator<float,400,8>& cal,
                  const GyroCalibration<float>& bias,const V (&targets)[4],int segment)
      : fit_(fit),bias_(bias),quiet_(cal,config_()),segment_(segment) {
    for(int i=0;i<4;++i)targets_[i]=targets[i];
  }
  void begin(uint32_t now) {
    fit_.cancelSegment();quiet_.begin(now);start_=now;first_sample_=fit_.samples();
    phase_=GyroTurnPhase::REST;reason_="none";armed_=moved_=twisted_=holding_=false;target_=-1;
    last_progress_=0;motion_=Eigen::Quaternionf::Identity();
  }
  GyroTurnPhase update(uint32_t now,uint32_t us,const V* a,const V* w,float temp,const V* mag) {
    if(phase_==GyroTurnPhase::DONE || phase_==GyroTurnPhase::FAILED)return phase_;
    if(uint32_t(now-start_)>=180000)return fail_("Turn timed out");
    const bool fresh=a && w;
    if(fresh && armed_) {
      // A slow/broken attempt cannot consume storage needed by later turns.
      if(fit_.samples()-first_sample_>=SAMPLES/GyroCalibrationPlan::MAXIMUM)return fail_("Turn took too long");
      const V rate=*w-bias_.biasT.bias(temp);
      if(!fit_.add(us,*w,temp)) {
        if(!moved_ && motion_.vec().norm()<.005f && last_rate_.norm()<.02f && rate.norm()<.02f &&
           fit_.reason==GyroRotationFail::GAP) {
          armed_=false;phase_=GyroTurnPhase::REST;quiet_.begin(now);last_progress_=0;
          reason_="sample gap";
        } else return fail_(gyroRotationReason(fit_.reason));
      } else {
        const V delta=.5f*(last_rate_+rate)*(uint32_t(us-last_us_)*1e-6f);
        const float angle=delta.norm();
        if(angle>0)motion_=(motion_*Eigen::Quaternionf(Eigen::AngleAxisf(angle,delta/angle))).normalized();
        moved_=moved_ || motion_.vec().norm()>.05f;
        // A twist while keeping the same side up excites cross-axis terms
        // that direct quarter-turns between faces alone may leave unobservable.
        twisted_=twisted_ || std::fabs(motion_.vec().dot(a0_.normalized()))>.5f;
        last_us_=us;last_rate_=rate;
      }
    }
    const auto state=quiet_.update(now,a,w,temp,mag);
    reason_=quiet_.resetReason();
    if(state==GyroCaptureStatus::STALE)return fail_("No fresh samples");
    if(state==GyroCaptureStatus::TIMEOUT)
      return fail_(std::strcmp(reason_,"none")?reason_:"Stillness not reached");
    if(state==GyroCaptureStatus::READY) {
      const V qa=quiet_.meanAccel(),qw=quiet_.meanGyro();const float qt=quiet_.meanTemperature();
      if((qw-bias_.biasT.bias(qt)).norm()>.01f || !quiet_.hasMagReference()) {
        // A human may still be placing the cube. Reject this quiet window,
        // not the turn; keep its target, trajectory and earlier accepted turns.
        quiet_.begin(now);last_progress_=0;reason_="Let it settle";
        return phase_;
      }
      const bool misplaced=armed_ && !twisted_ && qa.normalized().dot(a0_.normalized())<.85f;
      if(!armed_ || !moved_ || misplaced) {
        const bool new_target=!armed_ || misplaced;
        a0_=qa;w0_=qw;m0_=quiet_.meanMag();t0_=qt;
        if(new_target)chooseTarget_();
        if(target_<0)return fail_("No target side");
        if(!anchor_(us))return fail_(gyroRotationReason(fit_.reason));
        armed_=true;moved_=false;quiet_.begin(now);last_progress_=0;phase_=GyroTurnPhase::TWIST;
        return phase_;
      }
      // A named, generously sized target makes the action concrete. This is
      // only an acquisition region; the final fit still uses measured endpoints.
      if(!twisted_ || !atTarget_(qa)) {
        // Placement is a human action, not a failed rotation. Keep the named
        // target and recorded trajectory while a misplaced endpoint is fixed.
        quiet_.begin(now);last_progress_=0;
        reason_="Use shown side";
        phase_=holding_?GyroTurnPhase::HOLD:twisted_?GyroTurnPhase::TURN:GyroTurnPhase::TWIST;return phase_;
      }
      if(!fit_.finishSegment(qa,qt,USEFUL_MS*.001))return fail_(gyroRotationReason(fit_.reason));
      a1_=qa;w1_=qw;m1_=quiet_.meanMag();t1_=qt;
      phase_=GyroTurnPhase::DONE;return phase_;
    }
    const float progress=quiet_.progress();
    if(fresh && armed_ && !moved_ && progress>last_progress_ && quiet_.hasMagReference() &&
       (quiet_.meanAccel()-a0_).norm()<.06f && (quiet_.meanGyro()-w0_).norm()<.004f &&
       (quiet_.meanMag()-m0_).norm()<std::max(.6f,.015f*m0_.norm())) {
      if(!anchor_(us))return fail_(gyroRotationReason(fit_.reason));
    }
    last_progress_=progress;
    if(armed_ && fresh) {
      // No hold bar is displayed while a turn is requested. Only show the finish hold
      // once there is useful tilt toward the requested face. Movement then
      // restarts this hold, not the requested turn or any accepted turn.
      if(!holding_ && twisted_ && atTarget_(*a)) {
        holding_=true;
        // Give the finish pose its own complete settling/measurement window.
        // Neither the movement time nor an old quiet block counts toward it.
        quiet_.begin(now);last_progress_=0;reason_="none";
      }
      phase_=holding_?GyroTurnPhase::HOLD:!twisted_?GyroTurnPhase::TWIST:GyroTurnPhase::TURN;
    }
    return phase_;
  }
  GyroTurnPhase phase() const {return phase_;}
  GyroTurnPhase cancel() {return fail_("Retry requested");}
  int target() const {return target_;}
  // Include magnetic settling so the visible bar starts at the first quiet
  // block. Completion still requires all capture gates, never this estimate.
  float progress() const {return phase_==GyroTurnPhase::DONE?1.f:
    std::min(.99f,float(quiet_.quietMilliseconds())/HOLD_MS);}
  const char* reason() const {return reason_;}
  uint32_t resets() const {return quiet_.resets();}
  const V& accel() const {return a1_;}
  const V& gyro() const {return w1_;}
  const V& mag() const {return m1_;}
  float temperature() const {return t1_;}
private:
  static GyroCaptureCfg config_() {
    GyroCaptureCfg c;c.timeout_ms=90000;c.useful_ms=USEFUL_MS;c.min_samples=64;c.mag_window_ms=MAG_WINDOW_MS;
    // Endpoints determine gravity direction; they do not fit stationary bias.
    // Small braced hand tremor has appreciable rate scatter despite an accurate
    // mean gravity direction. Keep gravity, magnetic, freshness and final mean
    // rate gates, while reserving the original strict rate profile for bias.
    c.max_gyro_std=.08;c.max_gyro_change=.04;
    return c;
  }
  bool atTarget_(const V& a) const {
    return target_>=0 && a.allFinite() && a.norm()>1 &&
      a.normalized().dot(targets_[target_])>.819152f &&
      a.normalized().cross(a0_.normalized()).norm()>.2f;
  }
  void chooseTarget_() {
    const V up=a0_.normalized();
    for(int k=0;k<4;++k) {
      const int i=(segment_+k)%4;
      if(!targets_[i].allFinite() || std::fabs(targets_[i].norm()-1.f)>.01f)continue;
      const float dot=std::fabs(up.dot(targets_[i]));
      // Cycle through roughly perpendicular faces. Always taking the closest
      // perpendicular can alternate just two faces and leave an axis unexcited.
      if(dot<.6f){target_=i;return;}
    }
  }
  bool anchor_(uint32_t us) {
    fit_.cancelSegment();last_us_=us;last_rate_=w0_-bias_.biasT.bias(t0_);
    motion_=Eigen::Quaternionf::Identity();
    return fit_.beginSegment(us,a0_,t0_,GyroCalibrationPlan::validation(segment_),w0_);
  }
  GyroTurnPhase fail_(const char* why) {
    reason_=why;fit_.cancelSegment();phase_=GyroTurnPhase::FAILED;return phase_;
  }
  Fit& fit_;const GyroCalibration<float>& bias_;
  GyroCapture<float,400,8> quiet_;
  V targets_[4],a0_=V::Zero(),w0_=V::Zero(),m0_=V::Zero();float t0_=25;
  V a1_=V::Zero(),w1_=V::Zero(),m1_=V::Zero();float t1_=25;
  V last_rate_=V::Zero();Eigen::Quaternionf motion_=Eigen::Quaternionf::Identity();
  int segment_,target_=-1,first_sample_=0;uint32_t start_=0,last_us_=0;
  bool armed_=false,moved_=false,twisted_=false,holding_=false;float last_progress_=0;
  GyroTurnPhase phase_=GyroTurnPhase::REST;const char* reason_="none";
};
} // namespace imu_cal
