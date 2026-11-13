#pragma once
#include "StandardHeaders.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "util/W3dSimCommon.h"
#include "../kalman_carma/CarmaBenchTruth.h"
#include "../kalman_carma/carma3/Bench3.h"
#include "../common/FiniteMean.h"
namespace gyro_rollout {
using V=Eigen::Vector3f;using D=Eigen::Vector3d;
template<class F>decltype(auto) core(F& f){if constexpr(requires{f.raw();})return f.raw().mekf();else return f.mekf();}
template<class C>auto quat(const C& c){if constexpr(requires{c.quaternion_boat();})return c.quaternion_boat();else return c.quaternion();}
struct Stats{long n=0,enabled=0,mag=0,nees_count=0;double ba2=0,bg2=0,tilt2=0,nees=0,covered=0;FiniteMean nis;};
template<class Fusion,bool Invariant>int run(int seed,double duration,const std::string& prefix){
 using namespace carma_bench;using namespace carma3_bench;
 typename Fusion::Config cfg;
 cfg.sigma_a=V::Constant(2.8f*ACC_SIGMA*SIGMA_A_RESCALE);cfg.sigma_m=V::Constant(1.2f*MAG_SIGMA_UT*SIGMA_M_RESCALE);
 if constexpr(requires{cfg.sigma_g;})cfg.sigma_g=V::Constant(2.f*GYR_SIGMA*SIGMA_G_RESCALE);else cfg.gyro_noise_density=2.f*GYR_SIGMA*SIGMA_G_RESCALE;
 cfg.mag_delay_sec=3;cfg.mag_continuous_hard_iron=false;
 Fusion improved,legacy;improved.begin(cfg);cfg.gyro_sample_mode=ocean_imu::kalman::filter_detail::GyroSampleMode::IntervalRate;legacy.begin(cfg);
 Scenario sc;if(!base_scenario("jonswap_fp0.29_hs0.8_g3.3",sc))throw std::runtime_error("scenario missing");sc.hmodel.physical=true;sc.duration_s=duration;sc.att.yaw_wander_deg=0;sc.att.yaw_osc_deg=0;
 constexpr double dt=.005;WaveTruth truth(sc,w3d_expand_seed(seed,21u),dt);
 auto an=make_imu_noise_model(ACC_SIGMA,8e-3f*9.80665f,0.f,w3d_expand_seed(seed,11u),w3d_expand_seed(seed,12u));
 auto gn=make_imu_noise_model(GYR_SIGMA,.05f*float(PI/180),0.f,w3d_expand_seed(seed,13u),w3d_expand_seed(seed,14u));
 auto mn=make_mag_noise_model(MAG_SIGMA_UT,0.f,0.f,0.f,0.f,0.f,w3d_expand_seed(seed,15u),w3d_expand_seed(seed,16u));
 const Eigen::Matrix3d gauge=Eigen::AngleAxisd(-MagSim_WMM::default_declination_deg*PI/180,D::UnitZ()).toRotationMatrix();
 std::ofstream trace,scores;trace.exceptions(std::ios::failbit|std::ios::badbit);scores.exceptions(std::ios::failbit|std::ios::badbit);
 trace.open(prefix+"_trace.csv");scores.open(prefix+"_scores.csv");trace<<"t,profile,ba_h_error,tilt_deg,h_nees,enabled,acc_nis\n"<<std::setprecision(10);
 Stats stats[2];Fusion* filters[2]={&legacy,&improved};float phase=0;std::uint64_t hash=1469598103934665603ull;
 auto digest=[&](const V& v){for(int j=0;j<3;++j){float x=v[j];unsigned char b[4];std::memcpy(b,&x,4);for(auto c:b){hash^=c;hash*=1099511628211ull;}}};
 const double start=std::max(150.,.75*duration);
 for(long k=0;k<std::lround(duration/dt);++k){auto tr=truth.next();double t=k*dt;
  V a=apply_imu_noise(tr.f_b.cast<float>(),an,float(dt)),w=apply_imu_noise(tr.omega_b.cast<float>(),gn,float(dt));
  phase+=float(dt);bool tick=phase>=.04f;V m=V::Zero();if(tick){while(phase>=.04f)phase-=.04f;m=apply_mag_noise(tr.R_nb.transpose().cast<float>()*MagSim_WMM::mag_world_aero(),mn,.04f);}
  digest(a);digest(w);if(tick)digest(m);
  for(int i=0;i<2;++i){auto& f=*filters[i];f.update(float(dt),w,a,35.f);auto& c=core(f);
   if(tick){const auto before=c.covariance_full().eval();f.updateMag(m);if(c.lastMagDiag().accepted&&(before-c.covariance_full()).squaredNorm()>0)++stats[i].mag;}
   const D err=(c.get_acc_bias()-an.turn_on_bias).template cast<double>();const auto q=quat(c);const Eigen::Matrix3d R=q.template cast<double>().toRotationMatrix();
   double tilt=std::acos(std::clamp((R.transpose()*D::UnitZ()).dot((gauge*tr.R_nb).transpose()*D::UnitZ()),-1.,1.))*180/PI;
   double hn=NAN;
   if(k%200==0){Eigen::Matrix3d pb=c.covariance_full().template bottomRightCorner<3,3>().template cast<double>();if constexpr(Invariant)pb=R.transpose()*pb*R;
    const Eigen::Matrix2d h=pb.topLeftCorner<2,2>();hn=err.head<2>().dot(h.ldlt().solve(err.head<2>()));
    if(!std::isfinite(hn)||!c.covariance_full().allFinite())throw std::runtime_error("invalid covariance");
    if(t>=start){stats[i].nees+=hn;stats[i].covered+=hn<=5.99146;++stats[i].nees_count;}
    if(k%2000==0)trace<<t<<','<<(i?"default":"legacy")<<','<<err.head<2>().norm()<<','<<tilt<<','<<hn<<','<<c.acc_bias_updates_enabled()<<','<<c.lastAccDiag().nis<<'\n';
   }
   if(t>=start){auto& s=stats[i];++s.n;s.ba2+=err.head<2>().squaredNorm();s.bg2+=(c.gyroscope_bias()-gn.turn_on_bias).squaredNorm();s.tilt2+=tilt*tilt;s.enabled+=c.acc_bias_updates_enabled();s.nis.add(c.lastAccDiag().nis);}
  }
 }
 scores<<"profile,seed,duration,start,n,ba_h_rms,bg_rms,tilt_rms,enabled_fraction,acc_nis,h_nees,coverage,applied_mag,sensor_hash,nees_count,nis_count\n"<<std::setprecision(12);
 for(int i=0;i<2;++i){auto&s=stats[i];const double nd=double(s.nees_count);scores<<(i?"default":"legacy")<<','<<seed<<','<<duration<<','<<start<<','<<s.n<<','<<std::sqrt(s.ba2/s.n)<<','<<std::sqrt(s.bg2/s.n)<<','<<std::sqrt(s.tilt2/s.n)<<','<<double(s.enabled)/s.n<<','<<s.nis.mean()<<','<<s.nees/nd<<','<<s.covered/nd<<','<<s.mag<<','<<std::hex<<hash<<std::dec<<','<<s.nees_count<<','<<s.nis.count<<'\n';}
 trace.close();scores.close();
 std::cout<<"paired samples="<<std::lround(duration/dt)<<" sensor_hash="<<std::hex<<hash<<std::dec<<'\n';return 0;
}
int ou2(int,double,const std::string&);int ou3(int,double,const std::string&);int tfg(int,double,const std::string&);
}
