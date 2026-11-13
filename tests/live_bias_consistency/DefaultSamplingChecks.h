#pragma once
#include "GyroPredictionChecks.h"
template<class Fusion> decltype(auto) live_core(Fusion& f){
 if constexpr(requires{f.raw();})return f.raw().mekf();else return f.mekf();
}
template<class Core> auto body_quaternion(const Core& c){
 if constexpr(requires{c.quaternion_boat();})return c.quaternion_boat();else return c.quaternion();
}
// Independent interval inputs: constant rotation exercises actual defaults
// through handoff; changing rotation also exercises optional matched centering
// of the separate frontend. Native tests isolate varying-rate prediction.
template<class Fusion>void sampling_check(const char* name,bool centered_frontend){
 typename Fusion::Config cfg;check(cfg.gyro_sample_mode==Mode::Instantaneous,"live default is not instantaneous");
 check(cfg.frontend_gyro_sample_mode==Mode::IntervalRate,"frontend default changed");
 cfg.mag_delay_sec=0;
 if constexpr(requires{cfg.mag_continuous_hard_iron;})cfg.mag_continuous_hard_iron=false;
 if(centered_frontend)cfg.frontend_gyro_sample_mode=Mode::Instantaneous;
 Fusion sampled,interval;sampled.begin(cfg);cfg.gyro_sample_mode=Mode::IntervalRate;cfg.frontend_gyro_sample_mode=Mode::IntervalRate;interval.begin(cfg);
 Eigen::Vector3f prev=Eigen::Vector3f::Zero();double worst_q=0,worst_p=0;
 for(int k=0;k<36000;++k){const double t=k*.005;const V constant(.003,.002,.001);
  const M R=centered_frontend?rotation(.1*t):Eigen::AngleAxisd(constant.norm()*t,constant.normalized()).toRotationMatrix();
  const Eigen::Vector3f w=(centered_frontend?V(.1*rate(.1*t)):constant).cast<float>();
  const Eigen::Vector3f avg=k?Eigen::Vector3f(.5f*(prev+w)):(centered_frontend?Eigen::Vector3f::Zero():w);prev=w;
  const Eigen::Vector3f a=(R.transpose()*V(0,0,-9.80665)).cast<float>(),m=(R.transpose()*V(18,0,43)).cast<float>();
  sampled.update(.005f,w,a);interval.update(.005f,avg,a);
  if(k%8==7){sampled.updateMag(m);interval.updateMag(m);}
  auto& x=live_core(sampled);auto& y=live_core(interval);
  const double q_error=(body_quaternion(x).coeffs()-body_quaternion(y).coeffs()).norm();
  const double p_error=(x.covariance_full()-y.covariance_full()).norm();
  check(std::isfinite(q_error)&&std::isfinite(p_error),"nonfinite sampling comparison");
  worst_q=std::max(worst_q,q_error);worst_p=std::max(worst_p,p_error);
 }
 check(worst_q<1e-6&&worst_p<1e-6,"default path disagrees with external interval reference");
 std::cout<<name<<(centered_frontend?" optional frontend":" actual defaults")<<" q_error="<<worst_q<<" P_error="<<worst_p<<" PASS\n";
}
template<class Fusion>void default_sampling_check(const char* name){sampling_check<Fusion>(name,false);sampling_check<Fusion>(name,true);}
