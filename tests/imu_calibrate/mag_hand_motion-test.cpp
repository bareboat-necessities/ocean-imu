// Copyright 2026, Mikhail Grushinskiy
// Device-like magnetic cadence, LCD stalls and unconstrained hand turns.
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/MagCalSampling.h"
#include <random>
#include <cstdio>

using V=Eigen::Vector3f;
using M=Eigen::Matrix3f;
using Cal=imu_cal::MagCalibrator<float,400>;
static V direction(float t) {
  const float z=.94f*std::sin(.25f*t),p=.45f*t;
  return V(std::sqrt(1-z*z)*std::cos(p),std::sqrt(1-z*z)*std::sin(p),z);
}
static M distortion() {M d;d<<1.2,.18,.09,.18,.9,.12,.09,.12,1.05;return d;}
static V bias() {return V(8,-5,3);}

using Raw=imu_cal::MagCapture<float,400>::RawResidual;
// alternate: readings alternate between 0.8 and sqrt(2-0.64) times the field,
// whose mean square is the field squared (hidden by window moments).
static int capture(Cal& cal,int seed,int profile,int period,float field,bool noisy,
                   const imu_cal::MagCalibration<float>* fixed=nullptr,Raw* raw_out=nullptr,
                   bool alternate=false) {
  const bool verify=fixed!=nullptr;
  const auto cfg=imu_cal::MagSampleWindow::captureCfg(verify);
  imu_cal::MagCapture<float,400> cap(cal,cfg,fixed);cap.begin(0);
  imu_cal::MagSampleWindow window;window.begin();
  std::mt19937 rng(seed+(verify?19000:5000));std::normal_distribution<float> noise(0,1);
  auto random=[&](){return V(noise(rng),noise(rng),noise(rng));};
  const V sd=noisy?V(1,1,1.4f):V(.3f,.3f,.3f);
  V raw=V::Zero(),correlated=V::Zero();int next=0;
  for(int ms=0;ms<=int(cfg.timeout_ms);ms+=5) {
    float t=ms*.001f;
    if(profile==2) { // Fast turn followed by a brief natural pause.
      const float cycle=std::floor(t/1.6f);
      t=cycle*1.6f+2*std::min(.8f,t-cycle*1.6f);
    }
    t=(profile==0?3.f:4.f)*t+seed*3.17f+(verify?117:0);
    if(ms>=next) {
      next+=period;correlated=.8f*correlated+.06f*random();
      const float scale=alternate?((next/period)%2?.8f:std::sqrt(2-.64f)):1.f;
      raw=field*scale*distortion()*direction(t)+bias()+sd.cwiseProduct(random())+correlated;
      raw=(raw.array()/.3f).round().matrix()*.3f;
    }
    // Sensor continues sampling while the display blocks the wizard loop.
    if(profile!=0 && ms%250>=175)continue;
    V mean;M cov;const bool ready=window.update(ms,&raw,mean,cov);
    const auto status=cap.update(ms,ready?&mean:nullptr,&raw,ready?&cov:nullptr);
    if(status==imu_cal::MagCaptureStatus::READY){if(raw_out)*raw_out=cap.rawResidual();return ms;}
    if(status==imu_cal::MagCaptureStatus::TIMEOUT || status==imu_cal::MagCaptureStatus::STALE)break;
  }
  return 0;
}

int main() {
  int failures=0,passed=0,total=0,max_capture=0,max_verify=0;double max_raw_rms=0;
  double noise_error=0,quiet_error=0;int noise_count=0,quiet_count=0;
  for(int profile=0;profile<3;++profile)for(int period:{10,33,50,100})
    for(int seed=0;seed<5;++seed)for(bool noisy:{false,true}) {
      ++total;Cal cal;const float field=noisy?30.f:50.f;
      const int elapsed=capture(cal,seed,profile,period,field,noisy);
      imu_cal::MagCalibration<float> out;
      bool ok=elapsed && cal.fit(out);double error=0;
      if(ok) {
        for(int i=0;i<1000;++i) {
          const V u=direction(i*.217f),v=out.apply(field*distortion()*u+bias()).normalized();
          const double e=std::atan2(u.cross(v).norm(),u.dot(v))*180/3.141592653589793;
          error+=e*e;
        }
        error=std::sqrt(error/1000);
        if(noisy){noise_error+=error*error;++noise_count;}
        else {quiet_error+=error*error;++quiet_count;}
        const M frozen_A=out.A;const V frozen_b=out.b;
        Raw raw;
        const int verified=capture(cal,seed,profile,period,field,noisy,&out,&raw);
        imu_cal::MagFitQuality q;ok=verified && cal.check(out,q) && raw.ok;
        max_raw_rms=std::max(max_raw_rms,raw.rms);
        max_verify=std::max(max_verify,verified);
        ok=ok && out.A==frozen_A && out.b==frozen_b;
        if(ok) {
          for(int i=cal.buf.n/2;i<cal.buf.n;++i)cal.buf.v[i]+=V(15,0,0);
          ok=!cal.check(out,q); // Averaging cannot conceal a changed field.
        }
        if(!ok)std::fprintf(stderr,"verify failed profile=%d period=%d seed=%d noisy=%d time=%d n=%d gate=%s\n",
          profile,period,seed,int(noisy),verified,cal.buf.n,imu_cal::magFitGateText(q.gate));
      } else std::fprintf(stderr,"fit failed profile=%d period=%d seed=%d noisy=%d time=%d n=%d gate=%s\n",
          profile,period,seed,int(noisy),elapsed,cal.buf.n,imu_cal::magFitGateText(cal.quality.gate));
      max_capture=std::max(max_capture,elapsed);passed+=ok;failures+=!ok;
    }
  // Field changes inside every averaging window: the window moments match the
  // field, single readings do not.
  {
    Cal cal;const float field=50.f;
    imu_cal::MagCalibration<float> out;
    const bool fitted=capture(cal,1,0,33,field,false) && cal.fit(out);
    Raw raw;imu_cal::MagFitQuality q;
    const int verified=capture(cal,1,0,33,field,false,&out,&raw,true);
    const bool hidden=verified && cal.check(out,q);
    const bool caught=hidden && raw.n>=imu_cal::MagFitLimits::min_raw_samples && !raw.ok;
    std::printf("mag_hand_motion: in-window field change: window check %s, raw rms=%.2f uT -> %s\n",
        hidden?"passes":"fails",raw.rms,raw.ok?"accepted":"rejected");
    failures+=!fitted || !caught;
  }
  const double noisy_rms=std::sqrt(noise_error/std::max(1,noise_count));
  const double quiet_rms=std::sqrt(quiet_error/std::max(1,quiet_count));
  // Faster manual turns must still improve on the released slow-turn results.
  failures+=noisy_rms>=.9*1.067329960831 || quiet_rms>=.9*.149754226986;
  std::printf("mag_hand_motion: %d/%d fit + verify + field-change checks; direction RMS noisy=%.6f quiet=%.6f deg; max capture=%d ms verify=%d ms; max raw rms=%.3f uT\n",
      passed,total,noisy_rms,quiet_rms,max_capture,max_verify,max_raw_rms);
  return failures?1:0;
}
