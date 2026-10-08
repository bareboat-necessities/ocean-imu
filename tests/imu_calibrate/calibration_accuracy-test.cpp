// Copyright 2026, Mikhail Grushinskiy
// Identical sensor histories for v2.3.2 and the current device capture path.
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/CalibrateIMU.h"
#include "imu_calibrate/MagCalCapture.h"
#ifndef CAL_BASELINE_V232
#include "imu_calibrate/MagCalSampling.h"
#include "imu_calibrate/GyroCalCapture.h"
#endif
#include <cstdio>
#include <cmath>
#include <cstdint>
using V=Eigen::Vector3f;using M=Eigen::Matrix3f;
struct Rng {
 uint64_t state;
 double u(){state=state*6364136223846793005ULL+1442695040888963407ULL;return double((state>>11)+1)/9007199254740993.0;}
 // Draws are sequenced explicitly: operands and constructor arguments are
 // unsequenced in C++, so inline calls replay different histories per compiler.
 // This order reproduces the histories the v2.3.2 baseline below was taken on.
 float n(){const double r=u(),c=u();return float(std::sqrt(-2*std::log(r))*std::cos(6.283185307179586*c));}
 V noise(){const float z=n(),y=n(),x=n();return V(x,y,z);}
};
static V direction(float t) {float z=.94f*std::sin(.25f*t),p=.45f*t;return V(std::sqrt(1-z*z)*std::cos(p),std::sqrt(1-z*z)*std::sin(p),z);}
static M distortion(){M d;d<<1.2,.18,.09,.18,.9,.12,.09,.12,1.05;return d;}
int main(){
 int failures=0,mag_ok=0,gyro_ok=0;double mag_error=0,gyro_error=0,scenario_error[3]{};
 for(int seed=0;seed<30;++seed){
  const V bias(.012f,-.009f,.016f);Rng rng{uint64_t(7919+seed)};
  imu_cal::GyroCalibrator<float,400,8> gyro;
#ifndef CAL_BASELINE_V232
  imu_cal::GyroCapture<float,400,8> gc(gyro);gc.begin(0);
#endif
  for(int ms=0;ms<16000;ms+=5){
   V a=V(0,0,-gyro.g)+.01f*rng.noise(),w=bias+.0023f*rng.noise();
#ifdef CAL_BASELINE_V232
   // The released wizard waits 6.5 seconds then retains 220 raw observations.
   if(ms>=6500 && gyro.buf.n<220)gyro.addSample(w,a,25);
   if(gyro.buf.n==220)break;
#else
   if(gc.update(ms,&a,&w,25)==imu_cal::GyroCaptureStatus::READY)break;
#endif
  }
  imu_cal::GyroCalibration<float> g;bool gok=gyro.fit(g);gyro_ok+=gok;
  double ge=gok?(g.biasT.b0-bias).norm():1;gyro_error+=ge*ge;
  std::printf("gyro,%d,%d,%.9g\n",seed,int(gok),ge);
  for(int scenario=0;scenario<3;++scenario){
   // BMM150 low-power and high-accuracy noise, with anisotropy and quantization.
   const float field=scenario==0?30.f:50.f;
   const V sd=scenario==2?V(.3f,.3f,.3f):V(1,1,1.4f);
   const V offset(8,-5,3);Rng mr{uint64_t(5000+seed)};
   imu_cal::MagCalibrator<float,400> mag;
   imu_cal::MagCaptureCfg cfg;
#ifndef CAL_BASELINE_V232
   cfg=imu_cal::MagSampleWindow::captureCfg();
   imu_cal::MagSampleWindow window;window.begin();
#endif
   imu_cal::MagCapture<float,400> cap(mag,cfg);cap.begin(0);
   V raw=V::Zero(),correlated=V::Zero();int next=0;int elapsed=0;
   for(int ms=0;ms<=int(cfg.timeout_ms);ms+=5){
    float t=ms*.001f+seed*3.17f;const V u=direction(t);
    if(ms>=next){
     next+=scenario==2?50:33;
     correlated=.8f*correlated+.06f*mr.noise();
     raw=field*distortion()*u+offset+sd.cwiseProduct(mr.noise())+correlated;
     raw=(raw.array()/.3f).round().matrix()*.3f;
    }
#ifdef CAL_BASELINE_V232
    auto status=cap.update(ms,&raw);
#else
    V mean;M covariance;bool ready=window.update(ms,&raw,mean,covariance);
    auto status=cap.update(ms,ready?&mean:nullptr,&raw,ready?&covariance:nullptr);
#endif
    if(status==imu_cal::MagCaptureStatus::READY){elapsed=ms;break;}
   }
   imu_cal::MagCalibration<float> m;bool mok=elapsed && mag.fit(m);
   double error=0;
   if(mok)for(int i=0;i<1000;++i){V u=direction(i*.217f);V v=m.apply(field*distortion()*u+offset).normalized();
    double e=std::atan2(u.cross(v).norm(),u.dot(v))*180/3.141592653589793;error+=e*e;}
   error=mok?std::sqrt(error/1000):180;
#ifndef CAL_BASELINE_V232
   if(mok) {
    const M saved_A=m.A;const V saved_b=m.b;
    auto vc=imu_cal::MagSampleWindow::captureCfg(true);
    imu_cal::MagCapture<float,400> verify(mag,vc,&m);verify.begin(200000);
    window.begin();Rng vr{uint64_t(19000+seed)};next=0;correlated.setZero();
    bool complete=false;
    for(int ms=0;ms<=120000;ms+=5) {
     float t=ms*.001f+seed*3.17f+117;const V u=direction(t);
     if(ms>=next) {
      next+=scenario==2?50:33;correlated=.8f*correlated+.06f*vr.noise();
      raw=field*distortion()*u+offset+sd.cwiseProduct(vr.noise())+correlated;
      raw=(raw.array()/.3f).round().matrix()*.3f;
     }
     V mean;M covariance;bool ready=window.update(200000+ms,&raw,mean,covariance);
     if(verify.update(200000+ms,ready?&mean:nullptr,&raw,ready?&covariance:nullptr)==imu_cal::MagCaptureStatus::READY){complete=true;break;}
    }
    imu_cal::MagFitQuality q;
    mok=complete && mag.check(m,q);
    if(!mok)std::fprintf(stderr,"verification seed=%d scenario=%d gate=%s n=%d rms=%g\n",seed,scenario,imu_cal::magFitGateText(q.gate),mag.buf.n,q.rms);
    failures+=(m.A!=saved_A || m.b!=saved_b);
    if(mok && seed==0) {
     for(int i=mag.buf.n/2;i<mag.buf.n;++i)mag.buf.v[i]+=V(15,0,0);
     failures+=mag.check(m,q);
    }
   }
#endif
   scenario_error[scenario]+=error*error;
   mag_ok+=mok;mag_error+=error*error;
   std::printf("mag%d,%d,%d,%.9g,%d\n",scenario,seed,int(mok),error,elapsed);
#ifndef CAL_BASELINE_V232
   if(!mok)std::fprintf(stderr,"mag seed=%d scenario=%d gate=%s samples=%d inliers=%d rms=%g cells=%d\n",seed,scenario,imu_cal::magFitGateText(mag.quality.gate),mag.buf.n,mag.quality.inliers,mag.quality.rms,mag.quality.cells);
#endif
   failures+=!mok;
  }
  failures+=!gok;
 }
 std::printf("summary,gyro=%d/30,rms=%.9g,mag=%d/90,rms=%.9g\n",gyro_ok,std::sqrt(gyro_error/30),mag_ok,std::sqrt(mag_error/90));
#ifndef CAL_BASELINE_V232
 // Baseline generated by this driver against tag v2.3.2, commit
 // 2ba31e4e255a2b03f014cfbb9c9352d58d49f791 (same histories, same compiler).
 const double release_mag[3]={1.0673299608310005,.5453303241472823,.14975422698639443};
 failures+=std::sqrt(gyro_error/30)>=.8*.0002837803771843538;
 for(int i=0;i<3;++i)failures+=std::sqrt(scenario_error[i]/30)>=.9*release_mag[i];
#endif
 return failures?1:0;
}
