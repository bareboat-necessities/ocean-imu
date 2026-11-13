#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
// Noiseless multi-axis changing-rate tests, independent analytic trajectory.
#include "StandardHeaders.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <stdexcept>
#define private public
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#undef private
#include "GyroPredictionChecks.h"
namespace cs=ocean_imu::carma_shared_modes;
struct QProbe:cs::KalmanCarmaSharedModes<double>{using Base=cs::KalmanCarmaSharedModes<double>;using Base::Base;void bg(const V& b){this->xext.template segment<3>(3)=b;}};
void bias_jacobian(){
 const V prev(.2,-.3,.4),cur(.8,.3,-.1),bias(.001,-.002,.003);const double dt=.073;
 QProbe nominal(V::Constant(.03),V::Constant(.0003),V::Constant(.8));nominal.set_gyro_sample_mode(Mode::Instantaneous);nominal.prime_gyro_sample(prev,.005);nominal.bg(bias);
 QProbe predicted=nominal;predicted.time_update(cur,dt);M numeric;
 for(int j=0;j<3;++j){QProbe p=nominal,m=nominal;V bp=bias,bm=bias;bp[j]+=1e-6;bm[j]-=1e-6;p.bg(bp);m.bg(bm);p.time_update(cur,dt);m.time_update(cur,dt);
  Eigen::AngleAxisd plus(p.quaternion_boat().toRotationMatrix().transpose()*predicted.quaternion_boat().toRotationMatrix()),minus(m.quaternion_boat().toRotationMatrix().transpose()*predicted.quaternion_boat().toRotationMatrix());numeric.col(j)=(plus.angle()*plus.axis()-minus.angle()*minus.axis())/2e-6;
 }
 M rot,B;ocean_imu::kalman::filter_detail::rot_and_B_from_wt<double>(V(.5*(prev+cur)-bias),dt,rot,B);
 check((numeric-B).norm()<1e-9,"centered QMEKF gyro-bias Jacobian inconsistent");std::cout<<"QMEKF centered bias Jacobian error="<<(numeric-B).norm()<<'\n';
}
int main(){bias_jacobian();tests<cs::KalmanCarmaSharedModes<double>>("QMEKF");tests<cs::KalmanTfgCarmaSharedModes<double>>("TFG");
 using H=ocean_imu::kalman::filter_detail::GyroInterval<V>;H h;h.set_mode(Mode::Instantaneous);check(!h.prepare(V::Ones(),.005).predict,"first sample predicted a nonexistent interval");h.reset();check(!h.prepare(V::Zero(),.005).predict,"reset retained old sample");
 std::cout<<"timing tests PASS\n";
}
