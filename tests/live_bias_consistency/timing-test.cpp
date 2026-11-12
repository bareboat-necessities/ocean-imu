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
using V=Eigen::Vector3d;
using M=Eigen::Matrix3d;
using Q=Eigen::Quaterniond;
namespace cs=ocean_imu::carma_shared_modes;
using Mode=ocean_imu::kalman::filter_detail::GyroSampleMode;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
M rotation(double t){return (Eigen::AngleAxisd(.2*std::sin(1.3*t),V::UnitZ())*Eigen::AngleAxisd(.3*std::sin(1.7*t),V::UnitY())*Eigen::AngleAxisd(.4*std::sin(.8*t),V::UnitX())).toRotationMatrix();}
V rate(double t){const double r=.4*std::sin(.8*t),p=.3*std::sin(1.7*t),rd=.32*std::cos(.8*t),pd=.51*std::cos(1.7*t),yd=.26*std::cos(1.3*t);return V(rd-yd*std::sin(p),pd*std::cos(r)+yd*std::sin(r)*std::cos(p),-pd*std::sin(r)+yd*std::cos(r)*std::cos(p));}
template<class F>double trajectory(double h,bool centered,bool irregular,bool delta){
 F f(V::Constant(.03),V::Constant(.0003),V::Constant(.8));
 f.initialize_from_attitude(Q(rotation(0)),.03,.1);f.set_gyro_sample_mode(centered?Mode::Instantaneous:Mode::IntervalRate);
 if(centered)f.time_update(rate(0),h); // prime; no fabricated interval before t=0
 double t=0,worst=0;
 for(int k=0;t<12.;++k){double dt=irregular?h*(k%17==0?3.:k%2?.7:1.3):h;if(t+dt>12)dt=12-t;
  V input=rate(t+dt);if(delta){Eigen::AngleAxisd a(rotation(t).transpose()*rotation(t+dt));input=a.angle()*a.axis()/dt;}
  f.time_update(input,dt);t+=dt;
  worst=std::max(worst,Eigen::AngleAxisd(rotation(t).transpose()*f.quaternion_boat().toRotationMatrix()).angle());
 }
 return worst;
}
template<class F>void tests(const char* name){
 const double end=trajectory<F>(.02,false,false,false),mid=trajectory<F>(.02,true,false,false),fine=trajectory<F>(.01,true,false,false);
 const double variable=trajectory<F>(.01,true,true,false),exact=trajectory<F>(.02,false,true,true);
 std::cout<<name<<" endpoint="<<end<<" centered="<<mid<<" half_step="<<fine<<" irregular_and_drops="<<variable<<" delta_angle="<<exact<<'\n';
 check(mid<end*.05,"centering failed");check(fine<mid*.27,"second order convergence failed");check(variable<.0003,"irregular interpolation failed");check(exact<1e-10,"delta-angle input must not be averaged");
 F a(V::Constant(.03),V::Constant(.0003),V::Constant(.8)),b=a;
 a.set_gyro_sample_mode(Mode::Instantaneous);a.prime_gyro_sample(V(.2,-.3,.4),.005);
 const V raw(.4,.1,-.2);const double dt=.015;const double scale=2.;
 b.set_gyro_noise_density_rad_sqrt_s(V::Constant(.0003*std::sqrt(scale)));
 a.time_update(raw,dt);b.time_update(V(.3,-.1,.1),dt);
 check((a.quaternion_boat().coeffs()-b.quaternion_boat().coeffs()).norm()<1e-12,"nominal centered propagation mismatch");
 check((a.covariance_full()-b.covariance_full()).norm()<1e-11,"noise bound not propagated through native covariance");
 F la(V::Constant(.03),V::Constant(.0003),V::Constant(.8)),lb=la;
 la.set_imu_lever_arm_body(V(.2,-.1,.05));lb.set_imu_lever_arm_body(V(.2,-.1,.05));
 la.set_gyro_sample_mode(Mode::Instantaneous);la.prime_gyro_sample(V(.2,.1,-.3),.005);
 for(int k=0;k<5;++k){V endpoint(.2+.1*k,.1-.05*k,-.3+.07*k);la.time_update(endpoint,.005);lb.time_update(endpoint,.005);
  check((la.lever_.acceleration()-lb.lever_.acceleration()).norm()<1e-12,"centering delayed lever-arm correction");}

}
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
