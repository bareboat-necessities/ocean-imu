#define EIGEN_NON_ARDUINO
#include "kalman_common/WaveStartup.h"
#include <iostream>
#include <cmath>
int main() {
 using W=seastate::common::WaveStartup;
 using V=Eigen::Vector3f;
 using Mode=ocean_imu::kalman::filter_detail::GyroSampleMode;
 std::cout << "axis,step_uT,scaled,accepted,tilt_deg,yaw_deg,gyro_bias_radps\n";
 for(int axis=0;axis<3;++axis)for(float amplitude:{0.f,1.f,2.f,3.f,5.f,10.f})for(bool scaled:{false,true}) {
  W learner;
  for(int k=0;k<=7000;++k){
   learner.imu(.005f,V::Zero(),V(0,0,-9.80665f),Mode::Instantaneous);
   if(k%8==0){V field(20,0,48);if(k>=3000)field[axis]+=amplitude;learner.mag(field);}
  }
  W::Seed seed;const bool ok=learner.solve(9.80665f,scaled,seed);
  auto R=seed.attitude.toRotationMatrix();
  std::cout<<axis<<','<<amplitude<<','<<scaled<<','<<ok<<','<<std::acos(std::clamp(R(2,2),-1.f,1.f))*57.2957795f<<','<<std::atan2(R(1,0),R(0,0))*57.2957795f<<','<<seed.gyro_bias.norm()<<'\n';
 }
}
