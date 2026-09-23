#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define EIGEN_NON_ARDUINO
#include "kalman_tfg/Kalman3D_Wave_TFG.h"

int main() {
  using F = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
  using V = Eigen::Vector3f;
  constexpr float dt=0.005f, g=9.80665f;
  for (float bias : {0.0f, 0.01f, 0.03f, 0.05f}) {
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
    float maxp=0, maxv=0;
    for(int k=0;k<120000;k++){
      f.time_update(V::Zero(),dt);
      f.measurement_update_acc_only(V(0,0,-g+bias),35);
      if ((k%3)==0) f.applyIntegralZeroPseudoMeas();
      maxp=std::max(maxp,std::abs(f.get_position().z()));
      maxv=std::max(maxv,std::abs(f.get_velocity().z()));
    }
    std::cout<<"bias="<<bias<<" pz="<<f.get_position().z()<<" maxp="<<maxp
             <<" vz="<<f.get_velocity().z()<<" maxv="<<maxv<<" awz="<<f.get_world_accel().z()
             <<" baz="<<f.get_acc_bias().z()<<"\n";
    if (!std::isfinite(maxp) || maxp>20.0f) return 2;
  }
}
