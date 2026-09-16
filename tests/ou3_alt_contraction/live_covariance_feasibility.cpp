// Passive same-history MEKF covariance/gain feasibility diagnostic.
// No state/setter intervention and no source-uniform or rho assertion.
#define EIGEN_NON_ARDUINO 1
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include "util/W3dSimCommon.h"
#include "kalman_ou_common/KalmanOUCoreMath.h"
namespace ou3_alt_probe {
template<class F> void record(int,const F&);
template<class F,class V,class N,class S,class K,class Y> void innovation(int,const F&,const V&,const N&,const S&,const K&,const Y&);
template<class F> struct Scope { using C=std::remove_reference_t<F>; int k;const C&f;Scope(int k,const C&f):k(k),f(f){record(k,f);}~Scope(){record(k+1,f);} };
}
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std=9.80665f;
using ProbeMatrix=Eigen::Matrix<double,21,21>;
using Wrapper=SeaStateFusion_OU_III<TrackerType::KALMANF>;
static bool observing=false;
static int current_step=0;
static std::ofstream exact_rows;
uint32_t word(float x){uint32_t b;std::memcpy(&b,&x,4);return b;}
template<class F> void emit_state(int kind,const F&f){
 exact_rows<<current_step<<' '<<kind;
 for(int i=3;i<6;++i)for(int j=3;j<6;++j)exact_rows<<' '<<word(f.Pext(i,j));
 for(int i=3;i<6;++i)exact_rows<<' '<<word(f.xext(i));
}

struct Metrics {
 bool finite=true;
 double pmax=0,pmin=1e300,bgmax=0,awmax=0,kmax=0,dtheta=0,rmax=0,smin=1e300;
 double reset_ratio=0,reset_growth=0,reset_dtheta=0,joseph_growth=0,bg_pmax=0;
 int max_reset_step=0;
 ProbeMatrix before;
 double reset_before=0;
} metric;
double maxeig(const ProbeMatrix& p) { Eigen::SelfAdjointEigenSolver<ProbeMatrix> e(p,Eigen::EigenvaluesOnly);return e.eigenvalues().maxCoeff(); }
namespace ou3_alt_probe {
template<class F> void record(int kind,const F&f) {
 if(!observing)return;
 if(kind==1||kind==2||kind==13||kind==14||kind==23||kind==24||kind==33||kind==34){emit_state(kind,f);exact_rows<<'\n';}
 const bool finite=f.Pext.allFinite()&&f.xext.allFinite()&&f.qref.coeffs().allFinite();metric.finite&=finite;if(!finite)return;
 ProbeMatrix p=f.Pext.template cast<double>();
 Eigen::SelfAdjointEigenSolver<ProbeMatrix> e(p,Eigen::EigenvaluesOnly);
 metric.pmax=std::max(metric.pmax,e.eigenvalues().maxCoeff());metric.pmin=std::min(metric.pmin,e.eigenvalues().minCoeff());
 metric.bgmax=std::max(metric.bgmax,f.gyroscope_bias().template cast<double>().norm());
 metric.awmax=std::max(metric.awmax,f.get_world_accel().template cast<double>().norm());
 metric.bg_pmax=std::max(metric.bg_pmax,Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(p.template block<3,3>(3,3),Eigen::EigenvaluesOnly).eigenvalues().maxCoeff());
 if(kind==13||kind==23||kind==33)metric.before=p;
 if(kind==14||kind==24||kind==34)metric.joseph_growth=std::max(metric.joseph_growth,maxeig((p-metric.before).eval()));
 if(kind==40) { metric.reset_before=maxeig(p);metric.reset_dtheta=f.xext.template head<3>().template cast<double>().norm();metric.dtheta=std::max(metric.dtheta,metric.reset_dtheta); }
 if(kind==43) { double ratio=maxeig(p)/metric.reset_before; if(ratio>metric.reset_ratio){metric.reset_ratio=ratio;metric.max_reset_step=current_step;metric.reset_growth=metric.reset_dtheta;} }
}
template<class F,class V,class N,class S,class K,class Y> void innovation(int kind,const F&f,const V&r,const N&,const S&s,const K&k,const Y&) {
 if(!observing)return;
 emit_state(kind,f);
 for(int i=3;i<6;++i)for(int j=0;j<3;++j)exact_rows<<' '<<word(k(i,j));
 const auto R=kind==12?f.Racc:kind==22?f.Rmag:f.R_S;
 for(int i=0;i<3;++i)for(int j=0;j<3;++j)exact_rows<<' '<<word(R(i,j));
 for(int i=0;i<3;++i)exact_rows<<' '<<word(r(i));
 exact_rows<<'\n';
 metric.kmax=std::max(metric.kmax,k.template cast<double>().norm());metric.rmax=std::max(metric.rmax,r.template cast<double>().norm());
 if(s.allFinite())metric.smin=std::min(metric.smin,Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(s.template cast<double>(),Eigen::EigenvaluesOnly).eigenvalues().minCoeff());
}
}
int main(int argc,char**argv) {
 if(argc!=4)return 1;
 exact_rows.open(argv[3]);if(!exact_rows)return 1;
 std::cout<<std::setprecision(17);
 const bool active=argc>1&&std::string(argv[1])=="A";
 const int scenario=std::stoi(argv[2]);if(scenario<0||scenario>2)return 1;
 Wrapper f;Wrapper::Config cfg;f.begin(cfg);
 int startup=0,live_tick=0;
 const Eigen::Vector3f quiet(0,0,-g_std), mag(22,1,43);
 while(startup<70000) {
  ++startup;f.update(.005f,Eigen::Vector3f::Zero(),quiet);
  if(active&&startup%8==0)f.updateMag(mag);
  if(f.isLive()&&!live_tick)live_tick=startup;
  if(f.isLive()&&(!active||f.raw().mekf().acc_bias_updates_enabled()))break;
 }
 if(!f.isLive()||(active&&!f.raw().mekf().acc_bias_updates_enabled()))return 2;
 const auto&p=f.raw().mekf();ProbeMatrix p0=p.Pext.cast<double>();
 std::cout<<"seed "<<active<<' '<<scenario<<' '<<startup<<' '<<live_tick<<' '<<p.gyroscope_bias().norm()<<' '<<maxeig(p0)<<' '<<p0.block<3,3>(3,3).trace()<<'\n';
 emit_state(0,p);exact_rows<<'\n';
 double guard=0;observing=true;
 for(current_step=1;current_step<=600;++current_step){
  Eigen::Vector3f a=quiet,w=Eigen::Vector3f::Zero();
  if(scenario==0)w=Eigen::Vector3f::Constant(35.f);
  if(scenario==1)w=Eigen::Vector3f::Constant(current_step%2?35.f:-35.f);
  if(scenario==2){const float s=std::sin(float(current_step)*.005f);a=Eigen::Vector3f(150*s,150*s,-g_std+150*s);}
  if(!w.allFinite()||!a.allFinite()||w.cwiseAbs().maxCoeff()>35||a.cwiseAbs().maxCoeff()>160)return 3;
  f.update(.005f,w,a);if(active&&(startup+current_step)%8==0)f.updateMag(mag);
  guard=std::max(guard,double(f.raw().accelVibrationGuardEngagement()));
  emit_state(100,f.raw().mekf());exact_rows<<'\n';
  if(!metric.finite)break;
 }
 const auto&m=metric;
 std::cout<<"run "<<active<<' '<<scenario<<' '<<current_step<<' '<<m.finite<<' '<<guard<<' '<<m.pmax<<' '<<m.pmin<<' '<<m.bgmax<<' '<<m.awmax<<' '<<m.kmax<<' '<<m.rmax<<' '<<m.smin<<' '<<m.dtheta<<' '<<m.reset_ratio<<' '<<m.max_reset_step<<' '<<m.reset_growth<<' '<<m.joseph_growth<<' '<<m.bg_pmax<<'\n';
}
