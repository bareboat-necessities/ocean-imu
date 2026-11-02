#define EIGEN_NON_ARDUINO
#include "imu_calibrate/SensorNoise.h"
#include <random>
#include <cstdio>
static int failures=0;
static void check(bool b,const char* s){if(!b){++failures;std::fprintf(stderr,"FAIL %s\n",s);}}
int main(){
  using C=imu_cal::StationaryNoiseCapture<>;using S=imu_cal::NoiseCaptureStatus;
  C capture;capture.begin(0xffff0000u);imu_cal::SensorNoiseModel model;
  std::mt19937 rng(141);std::normal_distribution<float> n(0,1);
  uint32_t us=0xffff0000u;
  for(int i=0;i<=480002;++i) {
    us+=5000;
    Eigen::Vector3f a(0,0,-9.8025605f),w=Eigen::Vector3f::Zero();
    for(int k=0;k<3;++k){a[k]+=.002f/std::sqrt(.005f)*n(rng);w[k]+=.0002f/std::sqrt(.005f)*n(rng);}
    capture.add(us,a,w,25);
    if(i==2000)check(!capture.report(model),"short capture cannot produce noise model");
  }
  check(capture.status()==S::READY && capture.report(model),"40-minute white capture qualifies");
  for(int k=0;k<6;++k)check(std::fabs(model.white_density[k]/(k<3?.002:.0002)-1)<.10,"white intensity recovered");
  check(model.bias_mask==0 && model.walk_mask==0,"white trace does not invent long-term parameters");
  std::printf("noise: duration_ms=%u samples=%u accel_N=%.8g gyro_N=%.8g state_bytes=%zu\n",
    model.duration_ms,model.samples,model.white_density[0],model.white_density[3],sizeof(capture));
  C gap;gap.begin(0);gap.add(5000,Eigen::Vector3f(0,0,-9.80256),Eigen::Vector3f::Zero(),25);
  check(gap.add(20000,Eigen::Vector3f(0,0,-9.80256),Eigen::Vector3f::Zero(),25)==S::TIMING,"missing frame aborts capture");
  C frozen;frozen.begin(0);
  for(int i=1;i<400;++i)frozen.add(i*5000,Eigen::Vector3f(0,0,-9.80256),Eigen::Vector3f::Zero(),25);
  check(frozen.status()==S::FROZEN,"frozen register stream rejected");
  C moving;moving.begin(0);
  check(moving.add(5000,Eigen::Vector3f(0,0,-8),Eigen::Vector3f::Zero(),25)==S::MOVING,"motion rejected");

  // Independent stochastic sources: logarithmically spaced stationary OU
  // components approximate flicker over the middle of their time range;
  // integrated white increments provide a known rate-random-walk intensity.
  C mixed;mixed.begin(0);imu_cal::SensorNoiseModel long_term;
  Eigen::Vector3f flicker[8],walk=Eigen::Vector3f::Zero();double phi[8];
  constexpr double ou_sigma=.0015,walk_density=.00003;
  for(int j=0;j<8;++j) {
    phi[j]=std::exp(-.005/(.04*std::pow(4.,j)));
    for(int k=0;k<3;++k)flicker[j][k]=float(ou_sigma*n(rng));
  }
  for(int i=1;i<=480001;++i) {
    Eigen::Vector3f a(0,0,-9.8025605f),w;
    for(int k=0;k<3;++k) {
      a[k]+=.0008f/std::sqrt(.005f)*n(rng);
      for(int j=0;j<8;++j) {
        flicker[j][k]=float(phi[j]*flicker[j][k]+ou_sigma*std::sqrt(1-phi[j]*phi[j])*n(rng));
        a[k]+=flicker[j][k];
      }
      walk[k]+=float(walk_density*std::sqrt(.005)*n(rng));
      w[k]=walk[k]+.00002f/std::sqrt(.005f)*n(rng);
    }
    mixed.add(uint32_t(i)*5000,a,w,25);
  }
  check(mixed.report(long_term),"40-minute mixed-noise capture qualifies");
  check((long_term.bias_mask&7)==7,"identifiable flicker plateaus reported");
  check((long_term.walk_mask&56)==56,"identifiable gyro random walk reported");
  for(int k=0;k<3;++k) {
    // Sum of OU modes spaced by log(4) has an interior Allan plateau of
    // sigma*sqrt(2*log(2)/log(4)) = sigma, up to finite-band edge effects.
    check(std::fabs(long_term.bias_instability[k]/(ou_sigma/.664)-1)<.30,"flicker coefficient recovered");
    check(std::fabs(long_term.bias_random_walk[k+3]/walk_density-1)<.30,"random-walk intensity recovered");
  }
  std::printf("mixed noise: BI_mask=%u RW_mask=%u accel_B=%.8g gyro_K=%.8g\n",
    long_term.bias_mask,long_term.walk_mask,long_term.bias_instability[0],long_term.bias_random_walk[3]);
  return failures?1:0;
}
