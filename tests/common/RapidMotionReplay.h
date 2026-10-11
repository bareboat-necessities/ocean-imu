#pragma once
#include "RapidMotionBackends.h"
#include "util/ImuFreshness.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include "util/ImuLoopHelpers.h"
#include "util/MagneticInput.h"
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <memory>
#include <random>
using ocean_imu::attitude::rotation;
constexpr float rad = .017453292519943295f;
static void require(bool ok,const char* reason) {
    if (!ok) { std::fprintf(stderr,"Rapid backend %d: %s\n",RAPID_BACKEND,reason); std::exit(1); }
}
struct Trial {
    std::unique_ptr<Backend> b=std::make_unique<Backend>();
    ocean_imu::attitude::Protection protection;
    ocean_imu::attitude::AccelHold accel_hold;
    ocean_imu::magnetic::Input magnetic;
    atoms3r_ical::MagAcquisition timing{};
    uint32_t now=1000;
    V field=earth,acc=V(0,0,-gravity);
    float peak=0,tail=0,max_bias=0;
    uint32_t peak_time=0; bool peak_held=false, last_held=false;
    bool enabled=true;
    Trial() { magnetic.reset(.00135f,Eigen::Matrix3f::Identity()); }
    void tick(float dt,const Q& truth,const V& w,const V& measured_acc,const V& measured_mag,
              bool fresh_mag,bool fresh_acc=true,bool valid=true) {
        now+=uint32_t(std::lround(dt*1e6f));
        if (fresh_mag) { field=measured_mag; timing.acquired(now/1000,now,now,now); }
        if (fresh_acc) acc=measured_acc;
        const V a=accel_hold.update(acc,fresh_acc,V(w-b->bias()),dt);
        magnetic.advance(now,w-b->bias(),0,timing,field,V::Zero(),false);
        V current=V::Zero();
        const bool mag_ok=valid && magnetic.currentField(V::Zero(),current);
        const auto decision=protection.step(dt,w,b->bias(),a,gravity,b->attitude(),current,
            mag_ok,magnetic.newSource(),now);
        const V old_bias=b->bias();
        const bool hold=enabled && decision.hold_attitude; last_held=hold;
        b->step(dt,w,a,current,mag_ok,fresh_mag,!hold,fresh_acc);
        if (hold) require((b->bias()-old_bias).norm()<2e-6f,"bias changed while protected");
        if (enabled && decision.reacquire) {
            b->recover(decision.target);
            require((b->bias()-old_bias).norm()<2e-6f,"recovery changed learned bias");
        }
        require(b->finite() && b->attitude().coeffs().allFinite() && b->bias().allFinite(),"nonfinite state");
        tail=b->attitude().angularDistance(truth)/rad;
        if (tail>peak) { peak_time=now; peak_held=hold; }
        peak=std::max(peak,tail); max_bias=std::max(max_bias,b->bias().norm());
    }
};
static void rapid(float rate,float radius,bool combined,bool gap=false,bool enabled=true,int mag_delay=0,bool noisy=false) {
    Trial t; t.enabled=enabled; Q truth=Q::Identity(); V prev=V::Zero();
    std::deque<V> mag_history;
    std::mt19937 rng(7193); std::normal_distribution<float> normal(0,1);
    // Materialize before returning: an Eigen expression would retain a
    // reference to the temporary vector and invalidate the noisy-sensor test.
    auto noise=[&](float sigma) -> V { return sigma*V(normal(rng),normal(rng),normal(rng)); };
    for(int i=1;i<=3600;++i) {
        const float time=float(i)*.005f;
        float envelope=0;
        if(time>5 && time<5.5f) envelope=.5f-.5f*std::cos((time-5)*2*float(M_PI));
        else if(time>=5.5f && time<=10.5f) envelope=1;
        else if(time>10.5f && time<11) envelope=.5f+.5f*std::cos((time-10.5f)*2*float(M_PI));
        const float dt=gap && i==1400 ? .075f : .005f;
        V w(0,0,rate*rad*envelope);
        if(combined) w=rate*rad*envelope*V(.5f*std::sin(time),.35f*std::cos(.7f*time),.7f);
        truth=(truth*rotation(V(w*dt))).normalized();
        const V lever(radius,0,0), alpha=(w-prev)/dt; prev=w;
        const V a=truth.conjugate()*V(0,0,-gravity)+alpha.cross(lever)+w.cross(w.cross(lever));
        mag_history.push_back(truth.conjugate()*earth);
        if(int(mag_history.size())>mag_delay+1) mag_history.pop_front();
        const V wm=w+(noisy ? V(noise(.019f)+V(.005f,-.004f,.008f)) : V::Zero());
        const V am=a+(noisy ? V(noise(.12f)) : V::Zero());
        const V mm=mag_history.front()+(noisy ? V(noise(.5f)) : V::Zero());
        t.tick(dt,truth,wm,am,mm,i%8==0,i%2==0);
    }
    std::printf("rapid,%d,%d,%.0f,%.2f,%d,%d,%d,%.5f,%.5f,%.6f,%u,%d\n",RAPID_BACKEND,enabled,rate,radius,combined,gap,mag_delay*5,t.peak,t.tail,t.max_bias,t.protection.recoveryCount(),noisy);
    if(enabled) {
        if(t.peak>=12) std::fprintf(stderr,"peak at %.3fs held=%d\n",float(t.peak_time)*1e-6f,t.peak_held);
        require(t.peak<12,"rapid-motion attitude excursion exceeds 12 degrees");
        require(t.tail<5,"rapid-motion recovery exceeds 5 degrees after 7 seconds");
        require(t.max_bias<.06f,"dynamic acceleration polluted gyro bias");
    }
}
static void capture(float angle,float heading=0,bool noisy=false) {
    Trial t;
    const Q truth=Q(Eigen::AngleAxisf(heading*rad,V::UnitZ()))*
        Q(Eigen::AngleAxisf(-25*rad,V::UnitY()))*Q(Eigen::AngleAxisf(20*rad,V::UnitX()));
    t.b->recover(truth);
    const V a=truth.conjugate()*V(0,0,-gravity),m=truth.conjugate()*earth;
    std::mt19937 rng(9173); std::normal_distribution<float> normal(0,1);
    auto noise=[&](float sigma) -> V { return sigma*V(normal(rng),normal(rng),normal(rng)); };
    auto tick=[&](int i) {
        const V w=noisy ? V(noise(.019f)+V(.005f,-.004f,.008f)) : V::Zero();
        const V acc=a+(noisy ? noise(.12f) : V::Zero());
        const V mag=m+(noisy ? noise(.5f) : V::Zero());
        t.tick(.005f,truth,w,acc,mag,i%8==0,i%2==0);
    };
    for(int i=0;i<400;++i) tick(i);
    const V bias_before=t.b->bias();
    t.b->recover(Q(Eigen::AngleAxisf(angle*rad,V::UnitZ()))*truth);
    t.peak=0;
    for(int i=0;i<600;++i) tick(i);
    std::printf("capture,%d,%.2f,%.0f,%.5f,%.6f,%u,%d\n",RAPID_BACKEND,angle,heading,t.tail,(t.b->bias()-bias_before).norm(),t.protection.recoveryCount(),noisy);
    require(t.protection.recoveryCount()==1,"large heading error did not reacquire exactly once");
    require(t.tail<3,"large-angle capture outside 3 degrees after 3 seconds");
    require((t.b->bias()-bias_before).norm()<.025f,"large heading error became gyro bias");
}
static void staleAndDisturbance() {
    for(int mode=0;mode<3;++mode) {
        Trial t;
        for(int i=0;i<400;++i) t.tick(.005f,Q::Identity(),V::Zero(),V(0,0,-gravity),earth,i%8==0);
        t.b->recover(Q(Eigen::AngleAxisf(float(M_PI),V::UnitZ())));
        for(int i=0;i<400;++i) {
            V field=mode==2 ? V(earth*1.6f) : earth;
            t.tick(.005f,Q::Identity(),V::Zero(),V(0,0,-gravity),field,mode!=0 && i%8==0,true,mode!=1);
        }
        require(t.protection.recoveryCount()==0,"stale/invalid/disturbed field caused recovery");
    }
}
static void movingAfterTurn() {
    Trial t; Q q=Q::Identity();
    for(int i=0;i<2400;++i) {
        const float dt=.005f,time=float(i)*dt;
        const V w= time<2 ? V::Zero() : (time<3 ? V(0,0,4) : V(.18f*std::cos(time),0,0));
        q=(q*rotation(V(w*dt))).normalized();
        V a=q.conjugate()*V(0,0,-gravity);
        if(time>=3) a+=V(.5f*std::sin(time),0,.5f*std::sin(time));
        t.tick(dt,q,w,a,q.conjugate()*earth,i%8==0);
    }
    require(!t.last_held,"normal vessel motion was left indefinitely without aiding");
    require(t.tail<5,"moving recovery tail exceeds 5 degrees");
}
static void acquisition() {
    ocean_imu::ImuFreshness f; bool fresh;
    require(!f.observe(false,true,fresh),"gyro before first accel delivered");
    require(!f.observe(true,false,fresh),"accel-only incorrectly advances gyro");
    require(f.observe(false,true,fresh)&&fresh,"pending accel lost");
    require(f.observe(false,true,fresh)&&!fresh,"gyro-only frame dropped or accel counted twice");
    require(f.observe(true,true,fresh)&&fresh,"combined frame dropped");
    ocean_imu::ins::SampleDtTracker clock(.005f);
    clock.update(0xfffff000u);
    require(std::abs(clock.update(uint32_t(0xfffff000u+75000u))-.075f)<1e-6f,"wrap/gap elapsed time lost");
}
#if RAPID_BACKEND >= 3 && RAPID_BACKEND <= 5
static void selectiveCovarianceAndRecovery() {
    auto b=std::make_unique<Backend>(); auto& f=b->f;
    for(int i=0;i<100;++i) b->step(.005f,V::Zero(),V(0,0,-gravity),earth,true,i%8==0,true,true);
    const Q q=b->attitude(); const V bg=b->bias(),ba=f.get_acc_bias();
    f.set_attitude_aiding_enabled(false);
    const V before=f.get_world_accel();
    f.measurement_update_acc_only(V(3,0,-gravity));
    require((f.get_world_accel()-before).norm()>.01f,"protection discarded physical acceleration");
    require(q.angularDistance(b->attitude())<1e-6f,"held acceleration correction changed tilt");
    require((b->bias()-bg).norm()<1e-6f && (f.get_acc_bias()-ba).norm()<1e-6f,"held correction changed bias");
    auto psd=[&] {
        const Eigen::MatrixXf P=f.covariance_full();
        require(P.allFinite() && (P-P.transpose()).norm()<1e-4f,"nonsymmetric covariance");
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(P);
        require(solver.info()==Eigen::Success && solver.eigenvalues().minCoeff()>-1e-5f,"covariance lost PSD");
    };
    psd();
    const auto marginal=f.covariance_full();
    f.set_attitude_aiding_enabled(true);
    require((f.covariance_full().diagonal()-marginal.diagonal()).norm()<1e-6f,"resume changed covariance marginals");
    psd();
    const V pos=f.get_position(),vel=f.get_velocity(),aw=f.get_world_accel();
#if RAPID_BACKEND != 3
    const V integ=f.get_integral_displacement();
#endif
    const V reference=f.magnetic_reference_world();
    const Eigen::Matrix3f bias_cov=f.gyroscope_bias_covariance_body();
    const Q target=Q(Eigen::AngleAxisf(2.1f,V::UnitZ()))*Q(Eigen::AngleAxisf(.4f,V::UnitY()));
    b->recover(target);
    require(target.angularDistance(b->attitude())<1e-5f,"recovery quaternion frame/sign mismatch");
    require((pos-f.get_position()).norm()<1e-6f && (vel-f.get_velocity()).norm()<1e-6f &&
        (aw-f.get_world_accel()).norm()<1e-6f,"attitude recovery reset navigation history");
#if RAPID_BACKEND != 3
    require((integ-f.get_integral_displacement()).norm()<1e-6f,"recovery reset integral displacement");
#endif
    require((reference-f.magnetic_reference_world()).norm()<1e-6f,"recovery changed magnetic world frame");
    require((bias_cov-f.gyroscope_bias_covariance_body()).norm()<2e-5f,"recovery changed physical bias covariance");
    psd();
}

static void facadeProtection() {
    auto f=std::make_unique<Fusion>(); typename Fusion::Config cfg; f->begin(cfg);
#if RAPID_BACKEND == 5
    auto& raw=*f; auto coreOf=[&]() -> auto& {return f->mekf();};
    auto proxy=[&] {return f->startupProxyQuat();};
#else
    auto& raw=f->raw(); auto coreOf=[&]() -> auto& {return f->raw().mekf();};
    auto proxy=[&] {return f->raw().startupProxyQuat();};
#endif
    raw.setAttitudeAidingEnabled(false);
    f->update(.005f,V(0,0,4),V(7,0,-gravity));
#if RAPID_BACKEND == 5
    Q unused;
    require(!raw.startupTiltQuaternion(unused),"dynamic acceleration seeded the startup proxy");
#else
    require(!raw.startupProxyInitialized(),"dynamic acceleration seeded the startup proxy");
#endif
    raw.setAttitudeAidingEnabled(true);
    for(int i=0;i<26000;++i) {
        f->update(.005f,V::Zero(),V(0,0,-gravity));
        if(i%8==0) f->updateMag(earth);
    }
    require(f->isLive(),"facade did not reach Live");
    const V bg=coreOf().gyroscope_bias_body();
    Q expected=proxy();
    raw.setAttitudeAidingEnabled(false);
    for(int i=0;i<200;++i) {
        const V w(0,0,8); expected=(expected*rotation(V(w*.005f))).normalized();
        raw.setAccelerometerFresh(i%2==0);
        f->update(.005f,w,V(-9.6f,0,-gravity));
    }
    require((coreOf().gyroscope_bias_body()-bg).norm()<3e-6f,"facade leaked a bias correction");
    require(proxy().angularDistance(expected)<.01f,"private tilt proxy chased dynamic acceleration");
}
#endif
int main() {
    acquisition();
    for(float rate:{180.0f,360.0f,720.0f}) for(float r:{0.0f,.05f,.15f})
        for(bool combined:{false,true}) rapid(rate,r,combined);
    rapid(720,.15f,false,true); rapid(720,.15f,true,false,true,8);
    rapid(720,.15f,false,false,true,0,true);
    rapid(720,.15f,true,false,true,8,true);
    // Matched unprotected control: perfect sensors, no lag, same acceleration.
    rapid(720,.15f,false,false,false);
    for(float a:{90.0f,150.0f,179.99f,180.0f,180.01f,-179.99f}) for(float h:{0.0f,137.0f}) capture(a,h);
    for(float h:{0.0f,137.0f}) capture(180,h,true);
    staleAndDisturbance(); movingAfterTurn();
#if RAPID_BACKEND >= 3 && RAPID_BACKEND <= 5
    selectiveCovarianceAndRecovery();
    facadeProtection();
#endif
    std::fprintf(stderr,"Rapid-motion backend %d PASS\n",RAPID_BACKEND);
}
