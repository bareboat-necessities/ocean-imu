#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <string>
#include <stdexcept>
#include <array>
#include <fstream>
#include <cstdint>

static bool moving_recording=false;
static double moving_phi[21][4],moving_info[4][4];
static int moving_mags=0;
static constexpr int SLIDING=200;
struct MovingSlot { bool active=false; int root=-1; double phi[21][4]{}; double info[4][4]{}; int mags=0; };
static MovingSlot moving_slots[SLIDING];
static double sliding_service_min=1e100,sliding_service_max=0; static int sliding_windows=0,sliding_mags_min=10000;
static double pred_product=1.0,pred_eta_max=0.0; static int pred_count=0,s_due_count=0; static double tube_min_even=1e100,tube_min_odd=1e100,tube_phase_drift_even=0,tube_phase_drift_odd=0;
static Eigen::Matrix<double,21,21> lowerP; static bool lower_ready=false; static double lower_product=1.0,lower_eta_max=0; static int lower_count=0;
static const int PE[12]={1,4,6,8,9,11,12,14,15,17,18,20}, PO[9]={0,2,3,5,7,10,13,16,19};
static void lower_correct(const char* kind){
 if(!lower_ready)return;
 Eigen::Matrix<double,21,21> J=Eigen::Matrix<double,21,21>::Zero();
 if(kind[0]=='a'){
  // Coupled Cauchy ceiling for attitude/AW/BA, f<=30, sigma_a>=0.2.
  const double w[3]={30.,1.,1.}; const int off[3]={0,15,18};
  for(int a=0;a<3;++a)for(int b=0;b<3;++b)for(int d=0;d<3;++d)J(off[a]+d,off[b]+d)+=w[a]*w[b]/.04;
 } else if(kind[0]=='S'){
  for(int d=0;d<3;++d)J(12+d,12+d)=1./(.15*.15);
 } else if(kind[0]=='m'){
  for(int d=0;d<3;++d)J(d,d)=75.*75./(.8*.8);
 } else return;
 Eigen::LDLT<Eigen::Matrix<double,21,21>> ld(lowerP);
 if(ld.info()!=Eigen::Success)return;
 Eigen::Matrix<double,21,21> A=ld.solve(Eigen::Matrix<double,21,21>::Identity())+J;
 Eigen::LDLT<Eigen::Matrix<double,21,21>> la(A);
 if(la.info()!=Eigen::Success)return;
 lowerP=la.solve(Eigen::Matrix<double,21,21>::Identity()); lowerP=(lowerP+lowerP.transpose())*.5;
}
static void slot_prediction(MovingSlot& s,const float* fa,const float* fl,float phi){ double out[21][4]={};for(int i=0;i<6;++i)for(int j=0;j<4;++j)for(int k=0;k<6;++k)out[i][j]+=double(fa[i+6*k])*s.phi[k][j];for(int i=0;i<12;++i)for(int j=0;j<4;++j)for(int k=0;k<12;++k)out[i+6][j]+=double(fl[i+12*k])*s.phi[k+6][j];for(int i=18;i<21;++i)for(int j=0;j<4;++j)out[i][j]=double(phi)*s.phi[i][j];for(int i=0;i<21;++i)for(int j=0;j<4;++j)s.phi[i][j]=out[i][j]; }
static void slot_correction(MovingSlot& s,const char* kind,const float* h,const float* sv,const float* gain){ double y[3][4]={};for(int i=0;i<3;++i)for(int j=0;j<4;++j)for(int k=0;k<21;++k)y[i][j]+=double(h[i+3*k])*s.phi[k][j]; if(kind[0]=='m'){double l[3][3]={},w[3][4]={};for(int i=0;i<3;++i)for(int j=0;j<=i;++j){double v=double(sv[i+3*j]);for(int k=0;k<j;++k)v-=l[i][k]*l[j][k];if(!std::isfinite(v)||(i==j&&v<=0))throw std::runtime_error("sliding innovation not SPD");l[i][j]=i==j?std::sqrt(v):v/l[j][j];}for(int i=0;i<3;++i)for(int j=0;j<4;++j){double v=y[i][j];for(int k=0;k<i;++k)v-=l[i][k]*w[k][j];w[i][j]=v/l[i][i];}for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<3;++k)s.info[i][j]+=w[k][i]*w[k][j];++s.mags;} for(int i=0;i<21;++i)for(int j=0;j<4;++j)for(int k=0;k<3;++k)s.phi[i][j]-=double(gain[i+21*k])*y[k][j]; }
static void slot_reset(MovingSlot& s,const float* d){double g[3][3]={{1,-double(d[2])/2,double(d[1])/2},{double(d[2])/2,1,-double(d[0])/2},{-double(d[1])/2,double(d[0])/2,1}},o[3][4]={};for(int i=0;i<3;++i)for(int j=0;j<4;++j)for(int k=0;k<3;++k)o[i][j]+=g[i][k]*s.phi[k][j];for(int i=0;i<3;++i)for(int j=0;j<4;++j)s.phi[i][j]=o[i][j];}
static void moving_prediction(const float* fa,const float* fl,const float* qa,const float* ql,float phi,const float* qb,const float* pminus){
 if(lower_ready && lower_count<200){
  Eigen::Matrix<double,21,21> F=Eigen::Matrix<double,21,21>::Zero(),Q=Eigen::Matrix<double,21,21>::Zero();
  for(int i=0;i<6;++i)for(int j=0;j<6;++j){F(i,j)=double(fa[i+6*j]);Q(i,j)=double(qa[i+6*j]);}
  for(int i=0;i<12;++i)for(int j=0;j<12;++j){F(i+6,j+6)=double(fl[i+12*j]);Q(i+6,j+6)=double(ql[i+12*j]);}
  for(int i=0;i<3;++i){F(i+18,i+18)=double(phi);for(int j=0;j<3;++j)Q(i+18,j+18)=double(qb[i+3*j]);}
  Eigen::Matrix<double,21,21> C=F*lowerP*F.transpose(); C=(C+C.transpose())*.5;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double,21,21>> es(Q,C);
  if(es.info()==Eigen::Success){double eta=std::max(0.0,es.eigenvalues().maxCoeff());lower_eta_max=std::max(lower_eta_max,eta);lower_product/=1.+eta;++lower_count;}
  lowerP=C+Q; lowerP=(lowerP+lowerP.transpose())*.5;
 }
 if(pred_count<200){
  Eigen::Matrix<double,21,21> Q=Eigen::Matrix<double,21,21>::Zero(),Pm;
  for(int i=0;i<21;++i)for(int j=0;j<21;++j)Pm(i,j)=double(pminus[i+21*j]);
  for(int i=0;i<6;++i)for(int j=0;j<6;++j)Q(i,j)=double(qa[i+6*j]);
  for(int i=0;i<12;++i)for(int j=0;j<12;++j)Q(i+6,j+6)=double(ql[i+12*j]);
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)Q(i+18,j+18)=double(qb[i+3*j]);
  Eigen::Matrix<double,21,21> C=(Pm-Q);C=(C+C.transpose())*.5;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double,21,21>> es(Q,C);
  if(es.info()==Eigen::Success){double eta=std::max(0.0,es.eigenvalues().maxCoeff());pred_eta_max=std::max(pred_eta_max,eta);pred_product/=1.0+eta;++pred_count;}
 }
 double out[21][4]={};
 for(int i=0;i<6;++i)for(int j=0;j<4;++j)for(int k=0;k<6;++k)out[i][j]+=double(fa[i+6*k])*moving_phi[k][j];
 for(int i=0;i<12;++i)for(int j=0;j<4;++j)for(int k=0;k<12;++k)out[i+6][j]+=double(fl[i+12*k])*moving_phi[k+6][j];
 for(int i=18;i<21;++i)for(int j=0;j<4;++j)out[i][j]=double(phi)*moving_phi[i][j];
 for(int i=0;i<21;++i)for(int j=0;j<4;++j)moving_phi[i][j]=out[i][j];
 for(auto& s:moving_slots)if(s.active)slot_prediction(s,fa,fl,phi);
}
static void moving_correction(const char* kind,const float* h,const float* s,const float* gain){ if(kind[0]=='S')++s_due_count; lower_correct(kind);
 double y[3][4]={};for(int i=0;i<3;++i)for(int j=0;j<4;++j)for(int k=0;k<21;++k)y[i][j]+=double(h[i+3*k])*moving_phi[k][j];
 if(kind[0]=='m'){
  double l[3][3]={},w[3][4]={};
  for(int i=0;i<3;++i)for(int j=0;j<=i;++j){double v=double(s[i+3*j]);for(int k=0;k<j;++k)v-=l[i][k]*l[j][k];if(!std::isfinite(v)||(i==j && v<=0))throw std::runtime_error("magnetic innovation is not finite SPD");l[i][j]=i==j?std::sqrt(v):v/l[j][j];}
  for(int i=0;i<3;++i)for(int j=0;j<4;++j){double v=y[i][j];for(int k=0;k<i;++k)v-=l[i][k]*w[k][j];w[i][j]=v/l[i][i];}
  for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<3;++k)moving_info[i][j]+=w[k][i]*w[k][j];++moving_mags;
 }
 for(int i=0;i<21;++i)for(int j=0;j<4;++j)for(int k=0;k<3;++k)moving_phi[i][j]-=double(gain[i+21*k])*y[k][j];
 for(auto& sl:moving_slots)if(sl.active)slot_correction(sl,kind,h,s,gain);
}
static void moving_reset(const float* d){
 double g[3][3]={{1,-double(d[2])/2,double(d[1])/2},{double(d[2])/2,1,-double(d[0])/2},{-double(d[1])/2,double(d[0])/2,1}},o[3][4]={};
 for(int i=0;i<3;++i)for(int j=0;j<4;++j)for(int k=0;k<3;++k)o[i][j]+=g[i][k]*moving_phi[k][j];
 for(int i=0;i<3;++i)for(int j=0;j<4;++j)moving_phi[i][j]=o[i][j];
 for(auto& s:moving_slots)if(s.active)slot_reset(s,d);
}

#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std=9.80665f;
using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;
int main(int argc,char**argv){
 const double end=argc>1?std::stod(argv[1]):240.0;
 const char* profile_path=argc>2?argv[2]:nullptr; std::ofstream profile;
 if(profile_path){profile.open(profile_path,std::ios::binary);const char magic[8]={'O','U','3','P','R','F','1','\0'};profile.write(magic,8);uint32_t n=8000;profile.write(reinterpret_cast<const char*>(&n),4);}
 Fusion::Config cfg;
 // Current proof source profile. No startup/gate/tuner/scheduler overrides.
 cfg.sigma_a.setConstant(.2f);cfg.sigma_g.setConstant(.00135f);cfg.sigma_m.setConstant(.8f);
 Fusion f;f.begin(cfg);
 const double g=9.80665,sn=400.0/40001,cs=39999.0/40001,nu=acos(-1.)/10;
 int live=-1,refined=-1,active=-1,accepted=0;
 int wr=-1,windows=0,mag_min=10000;double service_min=1e100,service_max=0;
 double tail_pth_max=0,tail_pth_min=1e100,tail_pbg_max=0,tail_cross_max=0; int tail_samples=0;
 double parity_even_max=0,parity_odd_max=0,parity_off_max=0,tail_fhat_max=0;
 double bay=0,pbay=0,qxz=0,resmag=0,normlin=0,mintilt=10,maxmetric=0,minmetric=1e100;
 Eigen::Matrix<float,21,21> Pcycle=Eigen::Matrix<float,21,21>::Zero(); bool have_cycle=false;
 double cycle_tau=0,cycle_sigma=0,cycle_rs=0,cycle_period=0,cycle_elapsed=0;
  {static const double rootL[441]={4.9990999999999998e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9990999999999998e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9990999999999998e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9990999999999998e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9990999999999998e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9990999999999998e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2027825892825604e-08,0,0,8.4373649100302633e-09,0,0,8.5851348353044197e-10,0,0,1.021729136413388e-08,0,0,0,0,0,0,0,0,0,0,0,0,4.2027825892825604e-08,0,0,8.4373649100302633e-09,0,0,8.5851348353044197e-10,0,0,1.021729136413388e-08,0,0,0,0,0,0,0,0,0,0,0,0,4.2027825892825604e-08,0,0,8.4373649100302633e-09,0,0,8.5851348353044197e-10,0,0,1.021729136413388e-08,0,0,0,0,0,0,0,0,0,8.4373649100302914e-09,0,0,5.6207443707263824e-09,0,0,1.0441095677978756e-09,0,0,-4.6215150541444909e-09,0,0,0,0,0,0,0,0,0,0,0,0,8.4373649100302914e-09,0,0,5.6207443707263824e-09,0,0,1.0441095677978756e-09,0,0,-4.6215150541444909e-09,0,0,0,0,0,0,0,0,0,0,0,0,8.4373649100302914e-09,0,0,5.6207443707263824e-09,0,0,1.0441095677978756e-09,0,0,-4.6215150541444909e-09,0,0,0,0,0,0,0,0,0,8.5851348353044745e-10,0,0,1.044109567797876e-09,0,0,3.0433185785964316e-10,0,0,-8.7260875386651946e-10,0,0,0,0,0,0,0,0,0,0,0,0,8.5851348353044745e-10,0,0,1.044109567797876e-09,0,0,3.0433185785964316e-10,0,0,-8.7260875386651946e-10,0,0,0,0,0,0,0,0,0,0,0,0,8.5851348353044745e-10,0,0,1.044109567797876e-09,0,0,3.0433185785964316e-10,0,0,-8.7260875386651946e-10,0,0,0,0,0,0,0,0,0,1.0217291364133842e-08,0,0,-4.6215150541445075e-09,0,0,-8.7260875386652256e-10,0,0,1.7498052828455691e-08,0,0,0,0,0,0,0,0,0,0,0,0,1.0217291364133842e-08,0,0,-4.6215150541445075e-09,0,0,-8.7260875386652256e-10,0,0,1.7498052828455691e-08,0,0,0,0,0,0,0,0,0,0,0,0,1.0217291364133842e-08,0,0,-4.6215150541445075e-09,0,0,-8.7260875386652256e-10,0,0,1.7498052828455691e-08,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9999920000128004e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9999920000128004e-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.9999920000128004e-10};for(int j=0;j<21;++j)for(int i=0;i<21;++i)lowerP(i,j)=rootL[i+21*j];} lower_ready=true;
 for(int k=1;k<=int(end*200);++k){
  double t=k*.005,phase=nu*t,psi=.02*sin(phase),rate=.02*nu*cos(phase),ax=-.02*nu*nu*sin(phase);
  Eigen::Matrix3d U=Eigen::AngleAxisd(-psi,Eigen::Vector3d::UnitY()).toRotationMatrix();
  Eigen::Vector3f a=(U*Eigen::Vector3d(ax,0,-g*cs)).cast<float>();
  f.update(.005f,Eigen::Vector3f(0,float(rate),0),a,35);
  if(k%8==0){f.updateMag((75*U*Eigen::Vector3d::UnitX()).cast<float>());
   if(f.isLive()&&f.raw().mekf().lastMagDiag().accepted){++accepted;resmag=std::max(resmag,double(f.raw().mekf().lastMagDiag().r.norm()));}}
  if(live<0&&f.isLive())live=k;if(refined<0&&f.hasRefinedMagReference())refined=k;if(active<0&&f.raw().mekf().acc_bias_updates_enabled())active=k;
  const auto&m=f.raw().mekf();if(!m.Pext.allFinite()||!m.xext.allFinite()){std::cerr<<"nonfinite "<<k<<"\n";return 2;}
  bay=std::max(bay,std::abs(double(m.get_acc_bias().y())));pbay=std::max(pbay,double(m.Pext(19,19)));
  qxz=std::max(qxz,double(std::hypot(m.qref.x(),m.qref.z())));
  normlin=std::max(normlin,double(m.get_position().norm()));
  if(moving_recording && k-wr==200){
   for(int i=0;i<4;i+=2){double a=moving_info[i][i],b=moving_info[i][i+1],c=moving_info[i+1][i+1];double denom=a+c+std::hypot(a-c,2*b);if(!std::isfinite(denom)||denom<=0)throw std::runtime_error("invalid magnetic service matrix");double v=2*(a*c-b*b)/denom;if(!std::isfinite(v))throw std::runtime_error("nonfinite magnetic service eigenvalue");service_min=std::min(service_min,v);service_max=std::max(service_max,v);}
   ++windows;mag_min=std::min(mag_min,moving_mags);moving_recording=false;
  }
  if(active>=0 && k>=active+3400 && !moving_recording){
   moving_recording=true;wr=k;for(auto& row:moving_phi)for(auto& v:row)v=0;for(auto& row:moving_info)for(auto& v:row)v=0;moving_mags=0;
   for(int i=0;i<2;++i){Eigen::Vector3d direction=U*Eigen::Vector3d(0,(i==0?1:-1)*sn,cs);for(int j=0;j<3;++j){moving_phi[j][2*i]=direction[j];moving_phi[j+3][2*i+1]=.02*direction[j];}}
  }
  if(k>=int((end-21.0)*200) && k<=int((end-1.0)*200)){
   auto &sl=moving_slots[k%SLIDING];
   if(sl.active && k-sl.root==200){for(int i=0;i<4;i+=2){double aa=sl.info[i][i],bb=sl.info[i][i+1],cc=sl.info[i+1][i+1];double den=aa+cc+std::hypot(aa-cc,2*bb);double vv=2*(aa*cc-bb*bb)/den;sliding_service_min=std::min(sliding_service_min,vv);sliding_service_max=std::max(sliding_service_max,vv);}++sliding_windows;sliding_mags_min=std::min(sliding_mags_min,sl.mags);}
   sl=MovingSlot{};sl.active=true;sl.root=k;
   for(int i=0;i<2;++i){Eigen::Vector3d direction=U*Eigen::Vector3d(0,(i==0?1:-1)*sn,cs);for(int j=0;j<3;++j){sl.phi[j][2*i]=direction[j];sl.phi[j+3][2*i+1]=.02*direction[j];}}
  } else if(k>int((end-1.0)*200)){
   auto &sl=moving_slots[k%SLIDING];if(sl.active&&k-sl.root==200){for(int i=0;i<4;i+=2){double aa=sl.info[i][i],bb=sl.info[i][i+1],cc=sl.info[i+1][i+1];double den=aa+cc+std::hypot(aa-cc,2*bb);double vv=2*(aa*cc-bb*bb)/den;sliding_service_min=std::min(sliding_service_min,vv);sliding_service_max=std::max(sliding_service_max,vv);}++sliding_windows;sliding_mags_min=std::min(sliding_mags_min,sl.mags);sl.active=false;}
  }
  if(k==int((end-20.0)*200)){Pcycle=m.Pext;have_cycle=true;cycle_tau=f.raw().getTauApplied();cycle_sigma=f.raw().getSigmaApplied();cycle_rs=f.raw().getRSApplied();cycle_period=f.raw().getPseudoUpdatePeriodSec();cycle_elapsed=m.pseudo_update_elapsed_s_;}
  if(profile.is_open() && k>int((end-40.0)*200)){
   static const int EP[12]={1,4,6,8,9,11,12,14,15,17,18,20}; static const int OP[9]={0,2,3,5,7,10,13,16,19};
   float elapsed=m.pseudo_update_elapsed_s_,period_now=m.pseudo_update_period_s_;
   profile.write(reinterpret_cast<const char*>(&elapsed),4);profile.write(reinterpret_cast<const char*>(&period_now),4);
   for(int i=0;i<12;++i)for(int j=0;j<12;++j){float v=m.Pext(EP[i],EP[j]);profile.write(reinterpret_cast<const char*>(&v),4);}
   for(int i=0;i<9;++i)for(int j=0;j<9;++j){float v=m.Pext(OP[i],OP[j]);profile.write(reinterpret_cast<const char*>(&v),4);}
  }
  if(k>40000){
   Eigen::Matrix3d Pth=m.Pext.block<3,3>(0,0).cast<double>();
   Eigen::Matrix3d Pbg=m.Pext.block<3,3>(3,3).cast<double>();
   Eigen::Matrix3d Ptb=m.Pext.block<3,3>(0,3).cast<double>();
   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(Pth);
   if(es.info()!=Eigen::Success) throw std::runtime_error("attitude covariance eigensolve failed");
   tail_pth_min=std::min(tail_pth_min,es.eigenvalues().minCoeff());
   tail_pth_max=std::max(tail_pth_max,es.eigenvalues().maxCoeff());
   tail_pbg_max=std::max(tail_pbg_max,Pbg.norm());
   tail_cross_max=std::max(tail_cross_max,Ptb.norm());
   {static const int E0[12]={1,4,6,8,9,11,12,14,15,17,18,20};static const int O0[9]={0,2,3,5,7,10,13,16,19};
    Eigen::Matrix<double,12,12> Pe0;Eigen::Matrix<double,9,9> Po0;
    for(int i=0;i<12;++i)for(int j=0;j<12;++j)Pe0(i,j)=m.Pext(E0[i],E0[j]);
    for(int i=0;i<9;++i)for(int j=0;j<9;++j)Po0(i,j)=m.Pext(O0[i],O0[j]);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,12,12>> te(Pe0);Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,9,9>> to(Po0);
    tube_min_even=std::min(tube_min_even,te.eigenvalues().minCoeff());tube_min_odd=std::min(tube_min_odd,to.eigenvalues().minCoeff());
   }
   if(k%8==0){static const int E[12]={1,4,6,8,9,11,12,14,15,17,18,20}; static const int O[9]={0,2,3,5,7,10,13,16,19};
    Eigen::Matrix<double,12,12> Pe;Eigen::Matrix<double,9,9> Po;Eigen::Matrix<double,12,9> Peo;
    for(int i=0;i<12;++i)for(int j=0;j<12;++j)Pe(i,j)=m.Pext(E[i],E[j]);
    for(int i=0;i<9;++i)for(int j=0;j<9;++j)Po(i,j)=m.Pext(O[i],O[j]);
    for(int i=0;i<12;++i)for(int j=0;j<9;++j)Peo(i,j)=m.Pext(E[i],O[j]);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,12,12>> ee(Pe);Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,9,9>> eo(Po);
    parity_even_max=std::max(parity_even_max,ee.eigenvalues().maxCoeff());parity_odd_max=std::max(parity_odd_max,eo.eigenvalues().maxCoeff());parity_off_max=std::max(parity_off_max,Peo.norm());}
   tail_fhat_max=std::max(tail_fhat_max,double(m.accelerometer_measurement_func(35).norm()));
   ++tail_samples;
   if(m.Pext(19,19)<=0)throw std::runtime_error("nonpositive BA covariance marginal");double V=pow(g*sn-m.get_acc_bias().y(),2)/m.Pext(19,19);minmetric=std::min(minmetric,V);maxmetric=std::max(maxmetric,V);}
 }
 const auto&m=f.raw().mekf();std::cout<<std::setprecision(17)<<"{\"service_min\":"<<service_min<<",\"service_max\":"<<service_max<<",\"service_windows\":"<<windows<<",\"service_mags_min\":"<<mag_min<<",\"duration\":"<<end<<",\"live_step\":"<<live<<",\"refined_step\":"<<refined<<",\"active_step\":"<<active<<",\"accepted_mag\":"<<accepted<<",\"max_abs_bay\":"<<bay<<",\"max_Pbay\":"<<pbay<<",\"max_q_xz\":"<<qxz<<",\"max_mag_residual\":"<<resmag<<",\"max_position\":"<<normlin<<",\"tail_V_lower_min\":"<<minmetric<<",\"tail_V_lower_max\":"<<maxmetric<<",\"tau\":"<<f.raw().getTauApplied()<<",\"sigma\":"<<f.raw().getSigmaApplied()<<",\"R_S\":"<<f.raw().getRSApplied()<<",\"period\":"<<f.raw().getPseudoUpdatePeriodSec()<<",\"ref\":["<<m.v2ref.x()<<","<<m.v2ref.y()<<","<<m.v2ref.z()<<"],\"tail_pth_min\":"<<tail_pth_min<<",\"tail_pth_max\":"<<tail_pth_max<<",\"tail_pbg_norm_max\":"<<tail_pbg_max<<",\"tail_pth_bg_norm_max\":"<<tail_cross_max<<",\"tail_cov_samples\":"<<tail_samples<<",\"sliding_service_min\":"<<sliding_service_min<<",\"sliding_service_max\":"<<sliding_service_max<<",\"sliding_service_windows\":"<<sliding_windows<<",\"sliding_service_mags_min\":"<<sliding_mags_min<<",\"cycle_P_max_abs_diff\":"<<(have_cycle?(m.Pext-Pcycle).cwiseAbs().maxCoeff():-1)<<",\"cycle_tau_abs_diff\":"<<std::abs(double(f.raw().getTauApplied())-cycle_tau)<<",\"cycle_sigma_abs_diff\":"<<std::abs(double(f.raw().getSigmaApplied())-cycle_sigma)<<",\"cycle_RS_abs_diff\":"<<std::abs(double(f.raw().getRSApplied())-cycle_rs)<<",\"cycle_period_abs_diff\":"<<std::abs(double(f.raw().getPseudoUpdatePeriodSec())-cycle_period)<<",\"parity_even_eig_max\":"<<parity_even_max<<",\"parity_odd_eig_max\":"<<parity_odd_max<<",\"parity_off_fro_max\":"<<parity_off_max<<",\"tail_fhat_max\":"<<tail_fhat_max<<",\"cycle_scheduler_elapsed_abs_diff\":"<<std::abs(double(m.pseudo_update_elapsed_s_)-cycle_elapsed)<<",\"prediction_first200_product\":"<<pred_product<<",\"prediction_first200_count\":"<<pred_count<<",\"S_due_count\":"<<s_due_count<<",\"lower_factor_first200_product\":"<<lower_product<<",\"lower_factor_eta_max\":"<<lower_eta_max<<",\"lower_factor_count\":"<<lower_count<<",\"tube_tail_even_lambda_min\":"<<tube_min_even<<",\"tube_tail_odd_lambda_min\":"<<tube_min_odd<<"}\n";
}
