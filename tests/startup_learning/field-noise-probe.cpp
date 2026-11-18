#define EIGEN_NON_ARDUINO
#include "kalman_common/WaveStartup.h"
#include <iostream>
#include <random>
int main() {
 using W=seastate::common::WaveStartup;using V=Eigen::Vector3f;
 using Mode=ocean_imu::kalman::filter_detail::GyroSampleMode;
 std::cout<<"noise,axis,phase,seed,step,accepted\n";
 for(float noise:{0.f,.8f})for(int axis=0;axis<3;++axis)for(float phase:{10.f,10.125f,10.25f,10.375f,15.f,15.125f,15.25f,15.375f,25.f,25.125f,25.25f,25.375f})for(int seed=0;seed<10;++seed)for(float step:{0.f,5.f}) {
  std::mt19937 random(static_cast<unsigned>(seed)); std::normal_distribution<float> gaussian(0,1);
  W learner;
  for(int k=0;k<=7000;++k){
   learner.imu(.005f,V(.0008f,-.0004f,.0006f),V(0,0,-9.80665f),Mode::Instantaneous);
   if(k%8==0){V field(20,0,48);for(int i=0;i<3;++i)field[i]+=noise*gaussian(random);if(.005f*k>=phase)field[axis]+=step;learner.mag(field);}
  }
  W::Seed state;std::cout<<noise<<','<<axis<<','<<phase<<','<<seed<<','<<step<<','<<learner.solve(9.80665f,true,state)<<'\n';
 }
}
