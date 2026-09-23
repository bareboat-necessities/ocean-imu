/* Regression for the AtomS3R calibrated-sample temperature boundary. */
#define EIGEN_NON_ARDUINO
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#include <algorithm>
#include <cmath>
#include <iostream>
using F=ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
using V=Eigen::Vector3f;
float replay(float tempC) {
  constexpr float dt=0.005f,g=9.80665f;
  F f(0.00135f,g);
  f.initialize_identity(1e-4f);
  f.set_aw_time_constant(2.5f);
  f.set_aw_stationary_std(V::Constant(0.25f));
  f.set_RS_noise(V::Constant(1.5f));
  f.initialize_from_truth(Eigen::Quaternionf::Identity(),V::Zero(),V::Zero(),V::Zero(),V::Zero(),V::Zero(),V::Zero());
  f.set_initial_linear_uncertainty(1,2,5,1);
  f.set_initial_acc_bias_std(0.03f);
  f.set_linear_block_enabled(true);
  f.set_acc_bias_updates_enabled(false);
  f.reset_aw_covariance_to_stationary();
  float peak=0;
  // Input is already temperature compensated, exactly like runtime_.applyAccel().
  for(int k=0;k<24000;k++) {
    f.time_update(V::Zero(),dt);
    f.measurement_update_acc_only(V(0,0,-g),tempC);
    if((k%3)==0) f.applyIntegralZeroPseudoMeas();
    peak=std::max(peak,std::abs(f.get_position().z()));
  }
  return peak;
}
int main(){
  const float reference=replay(35.0f);
  const float cold=replay(15.0f);
  const float hot=replay(55.0f);
  std::cout<<"calibrated stationary: ref35="<<reference<<" cold15_wrong_temp="<<cold
           <<" hot55_wrong_temp="<<hot<<"\n";
  if(!(reference<1e-5f)) return 1;
  // A non-reference temperature must visibly reproduce the spurious forcing,
  // proving why the calibrated device path must not pass physical temperature.
  if(!(cold>0.1f && hot>0.1f)) return 2;
  return 0;
}
