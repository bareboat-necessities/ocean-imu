#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#include "DefaultSamplingChecks.h"
int main(){
 default_sampling_check<ocean_imu::ou2::SeaStateFusion_OU_II<TrackerType::KALMANF>>("OU-II");
 using Core=ocean_imu::kalman::Kalman3D_Wave_OU_II<double>;
 Core c(V::Constant(.03),V::Constant(.0003),V::Constant(.8));
 c.set_gyro_sample_mode(Mode::Instantaneous);c.prime_gyro_sample(V(2,0,0),.005);
 c.initialize_from_truth(V::Zero(),V::Zero(),Q::Identity(),V::Zero());
 const auto p=c.covariance_full().eval();c.time_update(V::Zero(),.005);
 check((c.quaternion_boat().coeffs()-Q::Identity().coeffs()).norm()<1e-12,"OU-II stale reinitialization gyro");
 check((c.covariance_full()-p).norm()==0,"OU-II reinitialization failed to prime");
}
