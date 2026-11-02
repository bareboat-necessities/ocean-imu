#define EIGEN_NON_ARDUINO
#define INPUT 0x01
#include "imu_calibrate/GyroCalFit.h"
#include <cstdio>
#include <random>
#include <chrono>
using Fit=imu_cal::GyroRotationFit<12000,24>;
using V=Eigen::Vector3d;using M=Eigen::Matrix3d;
static int failures=0;
static void check(bool b,const char* s){if(!b){++failures;std::fprintf(stderr,"FAIL %s\n",s);}}
static M truth(){M s;s<<1.035,.018,-.012,-.009,.975,.022,.014,-.017,1.02;return s;}
static void data(Fit& fit,bool degenerate=false,bool corrupt=false,bool noisy=false,unsigned seed=817,const imu_cal::AccelCalibration<float>* accelerometer=nullptr) {
  fit.clear();std::mt19937 rng(seed);std::normal_distribution<double> normal;
  uint32_t us=0xffff0000u;
  for(int i=0;i<20;++i) {
    V g0=degenerate?V(0,0,-9.8025605):V(normal(rng),normal(rng),normal(rng)).normalized()*9.8025605;
    V axis=degenerate?V::UnitX():V(normal(rng),normal(rng),normal(rng)).normalized();
    if(std::fabs(axis.dot(g0.normalized()))>.85)axis=g0.normalized().unitOrthogonal();
    const V omega=axis*(.8+.04*i),raw=truth().inverse()*omega+V(.006,-.009,.004);
    V observed=g0+(noisy?.002:0)*V(normal(rng),normal(rng),normal(rng));
    if(accelerometer)observed=accelerometer->S.cast<double>().inverse()*observed+
      accelerometer->biasT.bias(25).cast<double>();
    fit.beginSegment(us,observed.cast<float>(),25,i%4==3,Eigen::Vector3f(.006f,-.009f,.004f));
    V end=g0;double previous_fraction=0;
    for(int k=0;k<=200;++k) {
      // A stationary start followed by a finite, linear spin-up. Include the
      // first timestamp interval in the independent physical endpoint truth.
      const double fraction=std::min(1.,double(k+1)/(10+i%4));
      const V measured=V(.006,-.009,.004)+fraction*(raw-V(.006,-.009,.004))+
        (noisy?.002:0)*V(normal(rng),normal(rng),normal(rng));
      us+=5000;fit.add(us,measured.cast<float>(),25);
      end=Eigen::AngleAxisd(-omega.norm()*.005*.5*(previous_fraction+fraction),axis)*end;
      previous_fraction=fraction;
    }
    if(corrupt && i%4==3)end=Eigen::AngleAxisd(.15,V::UnitZ())*end;
    end+=(noisy?.002:0)*V(normal(rng),normal(rng),normal(rng));
    if(accelerometer)end=accelerometer->S.cast<double>().inverse()*end+
      accelerometer->biasT.bias(25).cast<double>();
    check(fit.finishSegment(end.cast<float>(),25),"segment stored");us+=10000;
  }
}
// Repeated cardinal tips identify three scales but leave the full nine-
// parameter model underdetermined. Validation uses different start directions.
static void scaleData(Fit& fit,const M& scales,bool corrupt=false,bool noisy=false,unsigned seed=71) {
  fit.clear();std::mt19937 rng(seed);std::normal_distribution<double> noise;
  uint32_t us=1000;
  const V bias(.006,-.009,.004);
  for(int i=0;i<20;++i) {
    const bool validation=i%4==3;const int axis_id=i%3;
    const V axis=V::Unit(axis_id),start=V::Unit((axis_id+(validation?2:1))%3)*9.8025605;
    const double angle=validation?1.1:1.5707963267948966;
    double integral=0,previous=0;
    for(int k=1;k<=200;++k){const double fraction=std::min(1.,k/10.);integral+=.0025*(previous+fraction);previous=fraction;}
    const V omega=axis*(angle/integral);
    fit.beginSegment(us,start.cast<float>(),25,validation,bias.cast<float>());
    for(int k=1;k<=200;++k) {
      const double fraction=std::min(1.,k/10.);us+=5000;
      const V rate=bias+scales.inverse()*(omega*fraction)+(noisy?.001:0)*V(noise(rng),noise(rng),noise(rng));
      check(fit.add(us,rate.cast<float>(),25),"scale fallback sample stored");
    }
    V end=Eigen::AngleAxisd(-angle,axis)*start;
    if(corrupt && validation)end=Eigen::AngleAxisd(.15,end.normalized().unitOrthogonal())*end;
    end+=(noisy?.002:0)*V(noise(rng),noise(rng),noise(rng));
    check(fit.finishSegment(end.cast<float>(),25),"scale fallback segment stored");us+=10000;
  }
}
int main() {
  // Analytic derivative through a sequence, tested independently of fitting.
  V u=V(.3,-.6,.7).normalized();const V start=u;
  Fit::J D=Fit::J::Zero();
  for(int k=0;k<30;++k)Fit::step(u,D,truth(),V(.004,-.006,.002));
  for(int k=0;k<9;++k) {
    M changed=truth();changed(k/3,k%3)+=1e-6;
    V p=start;Fit::J scratch=Fit::J::Zero();
    for(int j=0;j<30;++j)Fit::step(p,scratch,changed,V(.004,-.006,.002));
    check(((p-u)/1e-6-D.col(k)).norm()<1e-6,"SO3 derivative matches finite difference");
  }
  static Fit fit;data(fit);
  imu_cal::AccelCalibration<float> acc;acc.ok=true;acc.g=9.8025605f;acc.biasT.ok=true;
  imu_cal::GyroCalibration<float> gyro;gyro.ok=true;gyro.biasT.ok=true;gyro.biasT.b0<<.006f,-.009f,.004f;
  const auto begin=std::chrono::steady_clock::now();
  check(fit.fit(acc,gyro),"full matrix fit");
  const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
  check(gyro.matrix_ok && (gyro.S.cast<double>()-truth()).norm()<2e-5,"scale and cross axes recovered");
  check(gyro.rotation_validation_rms<2e-6 && gyro.rotation_train_segments==15 &&
        gyro.rotation_validation_segments==5,"held-out segments excluded from fit");
  std::printf("gyro rotation: matrix_error=%.8g validation_rms=%.8g fit_ms=%.3f storage=%zu samples=%d\n",
      (gyro.S.cast<double>()-truth()).norm(),double(gyro.rotation_validation_rms),ms,sizeof(fit),fit.samples());
  double worst_noisy=0;
  for(unsigned seed=0;seed<8;++seed) {
    data(fit,false,false,true,817+seed);
    const bool ok=fit.fit(acc,gyro);
    check(ok,"noisy multi-axis rotations qualify");
    const double error=(gyro.S.cast<double>()-truth()).norm();worst_noisy=std::max(worst_noisy,error);
    check(ok && error<.004 && gyro.rotation_validation_rms<.004,"noisy held-out rotations retain scale accuracy");
  }
  std::printf("gyro noisy campaign: seeds=8 identity_matrix_error=%.8g worst_matrix_error=%.8g\n",
    (M::Identity()-truth()).norm(),worst_noisy);
  auto distorted_acc=acc;distorted_acc.S<<1.02f,.008f,0,.008f,.99f,-.003f,0,-.003f,1.03f;
  distorted_acc.biasT.b0<<.05f,-.02f,.03f;
  data(fit,false,false,false,817,&distorted_acc);
  check(fit.fit(distorted_acc,gyro) && (gyro.S.cast<double>()-truth()).norm()<2e-5,
    "gyro fitting uses calibrated rather than raw gravity endpoints");
  data(fit,false,true);check(!fit.fit(acc,gyro) && !gyro.matrix_ok &&
    fit.reason==imu_cal::GyroRotationFail::VALIDATION,"unseen bad rotations reject");
  data(fit,true);check(!fit.fit(acc,gyro),"single-axis rotations cannot identify 9 parameters");
  const M scales=V(1.04,.975,1.02).asDiagonal();
  scaleData(fit,scales);
  check(!fit.fit(acc,gyro) && fit.reason==imu_cal::GyroRotationFail::GEOMETRY,
    "cardinal tips must not falsely qualify all nine coefficients");
  const int kept_samples=fit.samples();
  check(fit.fitScaleOnly(acc,gyro) && gyro.matrix_ok && gyro.scale_only && fit.samples()==kept_samples,
    "same captured data must qualify the three-scale fallback without recapture");
  check((gyro.S.cast<double>()-scales).norm()<2e-5 && gyro.rotation_validation_rms<2e-6,
    "diagonal scales must recover truth and pass independent directions");
  double worst_scale=0;
  for(unsigned seed=0;seed<8;++seed) {
    scaleData(fit,scales,false,true,71+seed);
    const bool ok=fit.fitScaleOnly(acc,gyro);
    const double error=(gyro.S.cast<double>()-scales).norm();worst_scale=std::max(worst_scale,error);
    check(ok && error<.002 && gyro.rotation_validation_rms<.004,"noisy independent scale recovery");
  }
  std::printf("gyro scale fallback: seeds=8 worst_scale_error=%.8g identity_error=%.8g\n",
    worst_scale,(scales-M::Identity()).norm());
  scaleData(fit,scales,true);
  check(!fit.fitScaleOnly(acc,gyro) && !gyro.matrix_ok && !gyro.scale_only &&
    fit.reason==imu_cal::GyroRotationFail::VALIDATION,"scale fallback must reject corrupted held-out rotations");
  data(fit,true);
  check(!fit.fitScaleOnly(acc,gyro),"scale fallback cannot qualify an unobserved axis");
  scaleData(fit,Eigen::AngleAxisd(.12,V::UnitZ()).toRotationMatrix()*scales);
  check(!fit.fitScaleOnly(acc,gyro),"scale fallback must reject significant cross-axis error");
  scaleData(fit,V(1.4,1.,1.).asDiagonal());
  check(!fit.fitScaleOnly(acc,gyro),"scale fallback retains physical plausibility limits");
  fit.clear();fit.beginSegment(1000,Eigen::Vector3f(0,0,-9.8),25,false,Eigen::Vector3f::Zero());
  check(!fit.add(90000,Eigen::Vector3f::Zero(),25) && fit.samples()==0,"gap invalidates entire segment");
  fit.beginSegment(1000,Eigen::Vector3f(0,0,-9.8),25,false,Eigen::Vector3f::Zero());
  check(!fit.add(6000,Eigen::Vector3f(NAN,0,0),25),"NaN rejected");
  check(!fit.beginSegment(1000,Eigen::Vector3f(0,0,-9.8),25,false,Eigen::Vector3f(NAN,0,0)),
    "nonfinite stationary start rate rejected");
  fit.beginSegment(1000,Eigen::Vector3f(0,0,-9.8),25,false,Eigen::Vector3f::Zero());
  check(fit.add(6000,Eigen::Vector3f::Zero(),25) && fit.samples()==1,
    "first timestamp interval is retained");
  return failures?1:0;
}
