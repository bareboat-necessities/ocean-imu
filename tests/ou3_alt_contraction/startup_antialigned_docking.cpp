// Finite constructive-source feasibility diagnostic, not an all-time proof.
// No wrapper state installation: the source controller reads a shadow observer
// initialized and driven by exactly the same packets as the public wrapper.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#define private public
#include "tuner/VerticalAccelComplementary.h"
#undef private
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std=9.80665f;
static constexpr double rates[]={
 .181930600531,.515399372842,.540507226271,.552654215095,.560708466752,
 .566385833303,.589351557302,.550697071356,.584179369245,.576549738968,
 .579617345047,.580925573491,.581591761612,.582454083615,.582985293928,
 .583455328412,.583430350242,.583885522476,.583319266784,.584607984154,
 .575530312911,.586086065677,.588020650327,.586860570603,.591085103126,
 .584217236131,-.140746326636,.210382990552,.040083210406,.111985855580};
static constexpr int dock=30002;
static constexpr double amplitude=2.9, omega=.6;
static double yaw(double t){
 double result=0;
 for(int j=0;j<30;j++){
  const double start=5.*j,end=j==29?150.01:5.*(j+1);
  if(t<end)return result+rates[j]*(t-start);
  result+=rates[j]*(end-start);
 }
 return result+omega*(t-150.01);
}

int main(int argc,char**argv){
 const int count=argc>2?std::atoi(argv[2]):30602;
 std::ofstream packets;if(argc>1)packets.open(argv[1]);
 using Wrapper=SeaStateFusion_OU_III<TrackerType::KALMANF>;
 Wrapper wrapper;Wrapper::Config cfg;wrapper.begin(cfg);
 VerticalAccelComplementary shadow(.2f,.02f,20.f);
 const float dt=.005f, alpha=1.f-std::exp(-dt/12.f);
 const double magnitude=std::sqrt(amplitude*amplitude+double(g_std)*g_std);
 Eigen::Vector3f lp=Eigen::Vector3f::Zero();
 Eigen::Vector3d reference, frozen;
 double max_noise=0,max_accel=0,max_angle=0,max_halferror=0,min_lpz=100;
 float max_guard=0,max_weight=0;bool integral_changed=false,live=false,shadow_differs=false;
 for(int k=0;k<count;k++){
  const double time=k*.005, phase=omega*time-yaw(time);
  const double rate=k<dock?rates[std::min(k/1000,29)]:omega;
  Eigen::Vector3f acc(amplitude*std::cos(phase),amplitude*std::sin(phase),-g_std);
  Eigen::Vector3f gyro(.019*std::cos(.1*time),.019*std::sin(.1*time),rate);
  if(k>=dock){
   const auto& m=shadow.ahrs_;
   Eigen::Quaterniond q(m.q0,m.q1,m.q2,m.q3);q.normalize();
   const Eigen::Vector3d h=q.conjugate()*Eigen::Vector3d::UnitZ();
   if(k==dock){reference=h;frozen=Eigen::Vector3d(m.integralFBx,m.integralFBy,m.integralFBz);}
   const Eigen::Vector3d physical_acc=acc.cast<double>();
   acc=(magnitude*h).cast<float>();
   max_accel=std::max(max_accel,(acc.cast<double>()-physical_acc).norm());
   float ax=-acc.x(),ay=-acc.y(),az=-acc.z();
   const float inverse=Mahony_AHRS<float>::invSqrt(ax*ax+ay*ay+az*az);
   ax*=inverse;ay*=inverse;az*=inverse;
   const float vx=m.q1*m.q3-m.q0*m.q2,vy=m.q0*m.q1+m.q2*m.q3;
   const float vz=.5f*(m.q0*m.q0-m.q1*m.q1-m.q2*m.q2+m.q3*m.q3);
   const float ex=ay*vz-az*vy,ey=az*vx-ax*vz,ez=ax*vy-ay*vx;
   const Eigen::Vector3d halferror(ex,ey,ez);
   const Eigen::Vector3d next_i(m.integralFBx+.02f*ex*dt,
      m.integralFBy+.02f*ey*dt,m.integralFBz+.02f*ez*dt);
   const Eigen::Vector3d total=next_i+Eigen::Vector3d(0,0,omega);
   const Eigen::Vector3d tangent=total-h.dot(total)*h;
   // k*dt=.08: a radius .0003 tube has about .000024 inward margin.
   const Eigen::Vector3d residual=-tangent-.2*halferror+16.*reference.cross(h);
   gyro=(Eigen::Vector3d(0,0,omega)+residual).cast<float>();
   max_noise=std::max(max_noise,(gyro.cast<double>()-Eigen::Vector3d(0,0,omega)).norm());
   max_angle=std::max(max_angle,(h-reference).norm());
   max_halferror=std::max(max_halferror,halferror.cwiseAbs().maxCoeff());
   integral_changed|=(next_i-frozen).squaredNorm()!=0;
  }
  if(packets)packets<<k+1<<std::hexfloat<<' '<<gyro.x()<<' '<<gyro.y()<<' '<<gyro.z()
      <<' '<<acc.x()<<' '<<acc.y()<<' '<<acc.z()<<'\n';
  wrapper.update(dt,gyro,acc);shadow.update(dt,gyro,acc,g_std);
  shadow_differs|=(wrapper.raw().startupProxyQuat().coeffs()-shadow.quaternion().coeffs()).squaredNorm()!=0;
  const auto leveled=seastate::common::accWorldFromBody(wrapper.raw().startupProxyQuat(),acc);
  if(k==0)lp=leveled;else lp+=alpha*(leveled-lp);
  max_guard=std::max(max_guard,wrapper.raw().accelVibrationRms());
  max_weight=std::max(max_weight,wrapper.raw().accelVibrationGuardEngagement());
  live|=wrapper.isLive();
  if(k==dock-1)std::cout<<std::setprecision(12)<<"dock "<<k+1<<" lpz "<<lp.z()
      <<" live "<<wrapper.isLive()<<"\n";
  if(k==dock-1 && argc>3){
   const auto& m=shadow.ahrs_;
   std::cerr<<std::hexfloat<<"dock_q "<<m.q0<<' '<<m.q1<<' '<<m.q2<<' '<<m.q3
       <<" dock_i "<<m.integralFBx<<' '<<m.integralFBy<<' '<<m.integralFBz<<'\n';
  }
  if(k>=dock)min_lpz=std::min(min_lpz,double(lp.z()));
 }
 std::cout<<std::setprecision(12)<<"tail "<<count<<" minlp "<<min_lpz
     <<" gyro "<<max_noise<<" accel "<<max_accel<<" angle "<<max_angle
     <<" halferror "<<max_halferror<<" guard "<<max_guard<<" weight "<<max_weight
     <<" integral_changed "<<integral_changed<<" shadow_differs "<<shadow_differs
     <<" live "<<live<<"\n";
 return count>dock && !live && !shadow_differs && !integral_changed && min_lpz>2
     && max_noise<.001 && max_accel<.01 && max_halferror<.000001
     && max_weight==0 && max_guard<.03 ? 0:2;
}
