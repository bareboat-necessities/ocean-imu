#pragma once
// Complete production facades, with the same acquisition/fusion ordering as
// the five sketches. Truth generates sensors only; no filter truth injection.
#define EIGEN_NON_ARDUINO
extern const float g_std=9.80665f;
#include "util/AttitudeProtection.h"
#include "util/MagneticInput.h"
#include "util/MagneticHeading.h"
#include "util/ImuLoopHelpers.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <memory>
#include <random>
using V=Eigen::Vector3f; using Q=Eigen::Quaternionf;
constexpr float gravity=9.8025605f;
const V earth(20.3f,0,46.8f);
#if HEAVE_BACKEND == 1
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
using Fusion=SeaStateFusion_OU_II<TrackerType::KALMANF>;
#elif HEAVE_BACKEND == 2
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;
#elif HEAVE_BACKEND == 3
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
using Fusion=ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
#elif HEAVE_BACKEND == 4
#include "nlo/TimeVarGainNLO_Adapter.h"
using Fusion=TimeVarGainNloAdapter<false,NloMagType::Magnetometer,float,TrackerType::PLL>;
#else
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
using Fusion=marine_obs::AdaptiveVerticalPIIMahony<float,true,TrackerType::PLL>;
#endif
struct Device {
 Fusion f;
#if HEAVE_BACKEND <= 3
 Device(){Fusion::Config cfg;cfg.gravity_magnitude=gravity;cfg.mag_delay_sec=0;
 cfg.sigma_a=V::Constant(.12f);cfg.sigma_m=V::Constant(.8f);cfg.mag_continuous_hard_iron=false;
#if HEAVE_BACKEND < 3
 cfg.sigma_g=V::Constant(.00135f);
 f.begin(cfg); f.raw().setTauScaledPseudoUpdateCadence(false);
#else
 cfg.gyro_noise_density=.00135f;f.begin(cfg);f.setTauScaledPseudoCadence(true);
#endif
 }
 auto& core(){
#if HEAVE_BACKEND < 3
 return f.raw().mekf();
#else
 return f.mekf();
#endif
 }
 Q attitude(){
#if HEAVE_BACKEND < 3
 return f.attitudeQuat();
#else
 return f.quaternion();
#endif
 }
 V bias(){return core().gyroscope_bias_body();}
 V hardIron(){return f.magHardIronBodyUT();}
 bool ready(){return f.isLive();}
 float height(){return -core().get_position().z();}
 void step(float dt,const V& w,const V& a,const V&,bool,bool aid){
#if HEAVE_BACKEND < 3
 f.raw().setAttitudeAidingEnabled(aid);f.raw().setAccelerometerFresh(true);f.update(dt,w,a);
#else
 f.setAttitudeAidingEnabled(aid);f.setAccelerometerFresh(true);f.update(dt,w,a,35);
#endif
 }
 void capture(ocean_imu::magnetic::Input& m,const V& a,const V& w){
#if HEAVE_BACKEND < 3
 m.capture(f.raw().startupProxyQuat(),attitude(),a,w);
#else
 m.capture(f.startupProxyQuat(),attitude(),a,w);
#endif
 }
 void correct(ocean_imu::magnetic::Input& m,uint32_t now){m.correct(f,now);}
 void recover(const Q& q){
#if HEAVE_BACKEND < 3
 core().initialize_from_attitude(q,.05f,.05f);
#else
 core().reacquire_attitude(q);
#endif
 }
#elif HEAVE_BACKEND == 4
 Device(){auto& c=f.config();c.gravity_mps2=gravity;c.filter.gravity_mps2=gravity;c.run_filter_before_initialized=false;c.filter.expected_mag_norm=0;c.auto_theta_from_wave_freq=true;f.reset();}
 Q attitude(){return f.snapshot().q_nb;} V bias(){return f.gyroscopeBiasBody();}
 V hardIron(){return V::Zero();} bool ready(){return f.initialized();}
 float height(){return f.snapshot().disp_zu.z();}
 void step(float dt,const V& w,const V& a,const V& m,bool valid,bool aid){f.setAttitudeAidingEnabled(aid);if(valid)f.setMagBody(m,true);else f.clearMag();f.update(dt,w,a);}
 void capture(ocean_imu::magnetic::Input&,const V&,const V&){}
 void correct(ocean_imu::magnetic::Input& m,uint32_t now){m.complete(now,true);}
 void recover(const Q& q){f.reacquireAttitude(q);}
#else
 bool seeded=false;
 static V map(const V& v){return V(v.x(),-v.y(),-v.z());}
 Device(){Fusion::Config c;c.gravity_mps2=gravity;c.use_mag=true;f.configure(c);f.reset();}
 Q attitude(){const auto& m=f.mahonyState();return Q(m.q0,m.q1,-m.q2,-m.q3).normalized();}
 V bias(){const auto& m=f.mahonyState();return map(V(-m.integralFBx,-m.integralFBy,-m.integralFBz));}
 V hardIron(){return V::Zero();}bool ready(){return seeded&&f.envelopeReady();}
 float height(){return f.displacement();}
 void step(float dt,const V& w,const V& a,const V& field,bool valid,bool aid){
 f.mahonyState().setAidingEnabled(aid);const V g=map(w),acc=map(a),mag=map(field);
 if(!seeded&&valid&&std::abs(acc.norm()/gravity-1)<.2f){
  const V up=acc.normalized(),north=(mag-up*mag.dot(up)).normalized(),west=up.cross(north).normalized();
  Eigen::Matrix3f axes;axes.col(0)=west.cross(up);axes.col(1)=west;axes.col(2)=up;
  const Q q(axes.transpose());auto& m=f.mahonyState();m.q0=q.w();m.q1=q.x();m.q2=q.y();m.q3=q.z();m.integralFBx=m.integralFBy=m.integralFBz=0;seeded=true;
 }
 if(valid)f.updateIMUMag(g.x(),g.y(),g.z(),acc.x(),acc.y(),acc.z(),mag.x(),mag.y(),mag.z(),dt);
 else f.updateIMU(g.x(),g.y(),g.z(),acc.x(),acc.y(),acc.z(),dt);
 }
 void capture(ocean_imu::magnetic::Input&,const V&,const V&){}
 void correct(ocean_imu::magnetic::Input& m,uint32_t now){m.complete(now,true);}
 void recover(const Q& q){f.reacquireAttitude(q.w(),q.x(),-q.y(),-q.z());}
#endif
};
void require(bool ok,const char* message){if(!ok){std::fprintf(stderr,"INS heave %d: %s\n",HEAVE_BACKEND,message);std::exit(1);}}
struct Metrics{double mean=0,rms=0,sine=0,cosine=0,hold=0;bool ready=false;};
Metrics replay(bool split,bool previous=false){
 auto d=std::make_unique<Device>();
 ocean_imu::attitude::Protection protection{previous?ocean_imu::attitude::MotionUse::Compass:ocean_imu::attitude::MotionUse::INS};
 ocean_imu::attitude::GyroWindow window;ocean_imu::magnetic::Input magnetic;
 magnetic.reset(.00135f,Eigen::Matrix3f::Identity()*.64f);
 atoms3r_ical::MagAcquisition timing;timing.delay_us=40000;timing.timing_allowance_s=.03f;
 ocean_imu::ins::MagFreshGate gate;
 std::mt19937 ar(675),gr(676),mr(677);std::normal_distribution<float> an(0,1),gn(0,1),mn(0,1);
 auto noise=[](auto& rng,auto& normal,float s)->V{return s*V(normal(rng),normal(rng),normal(rng));};
 V acc=V::Zero(),field=earth;std::deque<V> fields;Metrics result;unsigned n=0;uint32_t previous_now=0;
 constexpr float omega=2.51327412287f,amplitude=3/(omega*omega);
 for(int i=1;i<=48000;++i){
  const float t=float(i)*.005f;const uint32_t now=1000000u+uint32_t(i)*5000u;
  const Q truth=Q(Eigen::AngleAxisf(1.7f,V::UnitZ()))*Q(Eigen::AngleAxisf(.2f,V::UnitY()))*Q(Eigen::AngleAxisf(.3f+.06f*std::sin(.7f*t),V::UnitX()));
  fields.push_back(truth.conjugate()*earth);if(fields.size()>9)fields.pop_front();
  if(i%8==0){field=fields.front()+noise(mr,mn,.5f);timing.acquired(now/1000,now,now,now);}
  const bool fresh=i%2==0||i==1;
  if(fresh)acc=truth.conjugate()*V(0,0,-gravity+3*std::sin(omega*t))+noise(ar,an,.12f);
  const V gyro=V(.042f*std::cos(.7f*t),0,0)+noise(gr,gn,.003f)+V(.002f,-.003f,.005f);
  if(!split&&!fresh)continue;
  const float gyro_dt=previous_now?float(now-previous_now)*1e-6f:.005f;previous_now=now;
  protection.observeGyro(gyro_dt,gyro);window.push(gyro,gyro_dt);
  magnetic.advance(now,gyro-d->bias(),0,timing,field,d->hardIron(),false);
  if(!fresh)continue;
  V rate;float dt;window.take(rate,dt);V current=V::Zero();const bool valid=magnetic.currentField(d->hardIron(),current);
  const auto policy=protection.step(dt,rate,d->bias(),acc,true,gravity,d->attitude(),current,valid,magnetic.takeNewSource(),now,0,d->ready());
  const bool due=magnetic.due(gate,valid,now/1000);
  d->step(dt,rate,acc,current,valid,!policy.hold_attitude);d->capture(magnetic,acc,gyro);
#if HEAVE_BACKEND <= 3
  if(due)d->correct(magnetic,now);
#else
  (void)due;if(valid)d->correct(magnetic,now);
#endif
  if(policy.reacquire)d->recover(policy.target);
  require(std::isfinite(d->height())&&d->attitude().coeffs().allFinite(),"nonfinite state");
  if(t>180){const double h=d->height(),truth_h=amplitude*std::sin(omega*t);++n;result.mean+=h;result.rms+=(h-truth_h)*(h-truth_h);result.sine+=h*std::sin(omega*t);result.cosine+=h*std::cos(omega*t);result.hold+=policy.hold_attitude;}
 }
 result.mean/=n;result.rms=std::sqrt(result.rms/n);result.sine=2*result.sine/n/amplitude;result.cosine=2*result.cosine/n;result.hold/=n;result.ready=d->ready();
 std::printf("INS %d %s: ready=%d mean=%.5f rms=%.5f gain=%.5f cosine=%.5f held=%.3f\n",HEAVE_BACKEND,split?"split200/100":"combined100",result.ready,result.mean,result.rms,result.sine,result.cosine,result.hold);
 if(!previous){require(result.ready,"startup blocked");require(result.hold<.05,"waves latched aiding hold");require(std::abs(result.mean)<.1,"heave mean drift");require(result.sine>.5&&result.sine<1.4,"wave amplitude response");}
 return result;
}
int main(int argc,char**){const bool previous=argc>1;const auto a=replay(false,previous),b=replay(true,previous);if(!previous){require(std::abs(a.mean-b.mean)<.02,"split clock changed mean");require(std::abs(a.sine-b.sine)<.05,"split clock changed gain");require(std::abs(a.cosine-b.cosine)<.025,"split clock changed phase");require(b.rms<=1.15*a.rms+.01,"split clock increased RMS");}}
