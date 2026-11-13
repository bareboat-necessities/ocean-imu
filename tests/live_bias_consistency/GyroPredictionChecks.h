#pragma once
#include "StandardHeaders.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "kalman_common/GyroInterval.h"
using V=Eigen::Vector3d;
using M=Eigen::Matrix3d;
using Q=Eigen::Quaterniond;
using Mode=ocean_imu::kalman::filter_detail::GyroSampleMode;
inline void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
inline M rotation(double t){return (Eigen::AngleAxisd(.2*std::sin(1.3*t),V::UnitZ())*Eigen::AngleAxisd(.3*std::sin(1.7*t),V::UnitY())*Eigen::AngleAxisd(.4*std::sin(.8*t),V::UnitX())).toRotationMatrix();}
inline V rate(double t){const double r=.4*std::sin(.8*t),p=.3*std::sin(1.7*t),rd=.32*std::cos(.8*t),pd=.51*std::cos(1.7*t),yd=.26*std::cos(1.3*t);return V(rd-yd*std::sin(p),pd*std::cos(r)+yd*std::sin(r)*std::cos(p),-pd*std::sin(r)+yd*std::cos(r)*std::cos(p));}
template<class F>double trajectory(double h,bool centered,bool irregular,bool delta){
 F f(V::Constant(.03),V::Constant(.0003),V::Constant(.8));
 f.initialize_from_attitude(Q(rotation(0)),.03,.1);f.set_gyro_sample_mode(centered?Mode::Instantaneous:Mode::IntervalRate);
 if(centered)f.time_update(rate(0),h); // prime; no fabricated interval before t=0
 double t=0,worst=0;
 for(int k=0;t<12.;++k){double dt=irregular?h*(k%17==0?3.:k%2?.7:1.3):h;if(t+dt>12)dt=12-t;
  V input=rate(t+dt);if(delta){Eigen::AngleAxisd a(rotation(t).transpose()*rotation(t+dt));input=a.angle()*a.axis()/dt;}
  f.time_update(input,dt);t+=dt;
  const double error=Eigen::AngleAxisd(rotation(t).transpose()*f.quaternion_boat().toRotationMatrix()).angle();
  check(std::isfinite(error),"nonfinite attitude prediction");worst=std::max(worst,error);
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
 if constexpr (requires(F value) { value.lever_.acceleration(); }) {
 F la(V::Constant(.03),V::Constant(.0003),V::Constant(.8)),lb=la;
 la.set_imu_lever_arm_body(V(.2,-.1,.05));lb.set_imu_lever_arm_body(V(.2,-.1,.05));
 la.set_gyro_sample_mode(Mode::Instantaneous);la.prime_gyro_sample(V(.2,.1,-.3),.005);
 for(int k=0;k<5;++k){V endpoint(.2+.1*k,.1-.05*k,-.3+.07*k);la.time_update(endpoint,.005);lb.time_update(endpoint,.005);
  check((la.lever_.acceleration()-lb.lever_.acceleration()).norm()<1e-12,"centering delayed lever-arm correction");}
 }
}
