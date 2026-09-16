// Native finite-source counterexample to the 150s timeout-alignment premise.
// Construction and public API only; no installed filter state or modified config.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
const float g_std=9.80665f;
static double yaw(double t){
 if(t<5)return .2*t;
 if(t<125)return 1+.61*(t-5);
 if(t<130)return 74.2+.41*(t-125);
 if(t<135)return 76.25-.2*(t-130);
 return 75.25+.2*(t-135);
}
static double rate(double t){return t<5?.2:t<125?.61:t<130?.41:t<135?-.2:.2;}
int main(int argc,char** argv){
 std::ofstream samples; if(argc>1)samples.open(argv[1]);
 using Wrapper=SeaStateFusion_OU_III<TrackerType::KALMANF>;
 Wrapper wrapper;Wrapper::Config cfg;wrapper.begin(cfg);
 const double A=3.313664,w=.64;const float h=.005f,g=9.80665f,alpha=1.f-std::exp(-h/12.f);
 Eigen::Vector3f lp=Eigen::Vector3f::Zero();float maxguard=0,maxweight=0;
 int first_live=0;
 for(int k=0;k<50000;k++){
  double t=k*.005,phase=w*t-yaw(t);
  Eigen::Vector3f acc(A*std::cos(phase),A*std::sin(phase),-g);
  Eigen::Vector3f gyro(.019*std::cos(.1*t),.019*std::sin(.1*t),rate(t));
  if(samples){samples<<k+1<<std::hexfloat<<' '<<gyro.x()<<' '<<gyro.y()<<' '<<gyro.z()<<' '<<acc.x()<<' '<<acc.y()<<' '<<acc.z()<<'\n';}
  wrapper.update(h,gyro,acc);
  auto leveled=seastate::common::accWorldFromBody(wrapper.raw().startupProxyQuat(),acc);
  if(k==0)lp=leveled;else lp+=alpha*(leveled-lp);
  maxguard=std::max(maxguard,wrapper.raw().accelVibrationRms());
  maxweight=std::max(maxweight,wrapper.raw().accelVibrationGuardEngagement());
  if(k==30001||k==30601)std::cout<<std::setprecision(12)<<"step "<<k+1<<" live "<<wrapper.isLive()<<" proxy "<<wrapper.raw().startupProxyInitialized()<<" lpz "<<lp.z()<<" guard "<<maxguard<<" weight "<<maxweight<<"\n";
  if(wrapper.isLive()){first_live=k+1;std::cout<<"first_live "<<first_live<<" physical "<<(k+1)*.005<<" lpz "<<lp.z()<<"\n";break;}
 }
 return first_live>30602 && first_live<40000 && maxweight==0.f && maxguard<.03f ? 0 : 2;
}
