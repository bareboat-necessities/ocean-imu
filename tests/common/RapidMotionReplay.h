#pragma once
#include "RapidMotionBackends.h"
#include "util/ImuFreshness.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include "util/ImuLoopHelpers.h"
#include "util/MagneticInput.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <memory>
#include <random>
using ocean_imu::attitude::rotation;
constexpr float rad = .017453292519943295f;
static bool diagnosis_mode=false;
static bool known_force_control=false;
static bool recovery_diagnosis_mode=false;
/** @brief Exit with a backend-specific diagnostic when a replay invariant fails. */
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
    float peak_heading=0,tail_heading=0,peak_tilt=0,tail_tilt=0;
    uint32_t peak_time=0; bool peak_held=false, last_held=false;
    unsigned applied_recoveries=0;
    bool enabled=true;
    /** @brief Initialize magnetic preprocessing for an independent replay trial. */
    Trial() { magnetic.reset(.00135f,Eigen::Matrix3f::Identity()); }
    /** @brief Deliver one timed sensor frame, apply optional protection, and record attitude errors.
     * Assert that held updates and recovery preserve gyro bias and keep state finite.
     */
    void tick(float dt,const Q& truth,const V& w,const V& measured_acc,const V& measured_mag,
              bool fresh_mag,bool fresh_acc=true,bool valid=true) {
        now+=uint32_t(std::lround(dt*1e6f));
        if (fresh_mag) { field=measured_mag; timing.acquired(now/1000,now,now,now); }
        if (fresh_acc) acc=measured_acc;
        const V a=accel_hold.update(acc,fresh_acc,V(w-b->bias()),dt);
        magnetic.advance(now,w-b->bias(),0,timing,field,V::Zero(),false);
        V current=V::Zero();
        const bool mag_ok=valid && magnetic.currentField(V::Zero(),current);
        const auto decision=protection.step(dt,w,b->bias(),a,fresh_acc,gravity,b->attitude(),current,
            mag_ok,magnetic.newSource(),now);
        const V old_bias=b->bias();
        const bool hold=enabled && decision.hold_attitude; last_held=hold;
        b->step(dt,w,a,current,mag_ok,fresh_mag,!hold,fresh_acc);
        if (hold) require((b->bias()-old_bias).norm()<2e-6f,"bias changed while protected");
        if (enabled && decision.reacquire) {
            b->recover(decision.target);
            ++applied_recoveries;
            require((b->bias()-old_bias).norm()<2e-6f,"recovery changed learned bias");
#if RAPID_BACKEND == 1
            if (decision.repair_gyro_bias) b->f.reacquire_gyro_bias(decision.gyro_bias);
#elif RAPID_BACKEND == 5
            if (decision.repair_gyro_bias) b->f.recover_corrupted_biases(decision.gyro_bias);
#endif
        }
        require(b->finite() && b->attitude().coeffs().allFinite() && b->bias().allFinite(),"nonfinite state");
        tail=b->attitude().angularDistance(truth)/rad;
        const V forward=b->attitude()*V::UnitX(), true_forward=truth*V::UnitX();
        tail_heading=std::abs(std::remainder(std::atan2(forward.y(),forward.x())-
            std::atan2(true_forward.y(),true_forward.x()),2*float(M_PI)))/rad;
        const V down=b->attitude().conjugate()*V::UnitZ();
        const V true_down=truth.conjugate()*V::UnitZ();
        tail_tilt=std::atan2(down.cross(true_down).norm(),down.dot(true_down))/rad;
        peak_heading=std::max(peak_heading,tail_heading); peak_tilt=std::max(peak_tilt,tail_tilt);
        if (tail>peak) { peak_time=now; peak_held=hold; }
        peak=std::max(peak,tail); max_bias=std::max(max_bias,b->bias().norm());
    }
};
/** @brief Replay a synthetic turn with rate in deg/s, radius in meters, and optional sensor faults.
 * Emit error metrics and enforce bounds when protection is enabled.
 */
static void rapid(float rate,float radius,bool combined,bool gap=false,bool enabled=true,int mag_delay=0,bool noisy=false,bool protect_after_turn=false,bool moving_tail=false,int steps=3600,bool calibrated_delay=false) {
    Trial t; t.enabled=enabled; Q truth=Q::Identity(); V prev=V::Zero();
    if (calibrated_delay) {
        t.timing.delay_us=uint32_t(mag_delay*5000);
        t.timing.timing_allowance_s=.03f;
    }
    float inherited_error=0, inherited_bias=0;
    std::deque<V> mag_history;
    std::mt19937 rng(7193); std::normal_distribution<float> normal(0,1);
    // Materialize before returning: an Eigen expression would retain a
    // reference to the temporary vector and invalidate the noisy-sensor test.
    auto noise=[&](float sigma) -> V { return sigma*V(normal(rng),normal(rng),normal(rng)); };
    for(int i=1;i<=steps;++i) {
        const float time=float(i)*.005f;
        // Continue the very same estimator, covariance, bias, field reference,
        // frontend and guard history. No attitude injection at this boundary.
        if (protect_after_turn && i==2201) t.enabled=true;
        float envelope=0;
        if(time>5 && time<5.5f) envelope=.5f-.5f*std::cos((time-5)*2*float(M_PI));
        else if(time>=5.5f && time<=10.5f) envelope=1;
        else if(time>10.5f && time<11) envelope=.5f+.5f*std::cos((time-10.5f)*2*float(M_PI));
        const float dt=gap && i==1400 ? .075f : .005f;
        V w(0,0,rate*rad*envelope);
        if(combined) w=rate*rad*envelope*V(.5f*std::sin(time),.35f*std::cos(.7f*time),.7f);
        // Resume motion after the same seven-second recovery tail.
        if(moving_tail && time>18) w=V(.18f*std::sin(time-18),.04f*std::sin(.7f*(time-18)),0);
        truth=(truth*rotation(V(w*dt))).normalized();
        const V lever(radius,0,0), alpha=(w-prev)/dt; prev=w;
        V a=truth.conjugate()*V(0,0,-gravity)+alpha.cross(lever)+w.cross(w.cross(lever));
        if(moving_tail && time>18) a+=truth.conjugate()*V(.3f*std::sin(time-18),0,.3f*std::sin(time-18));
        mag_history.push_back(truth.conjugate()*earth);
        if(int(mag_history.size())>mag_delay+1) mag_history.pop_front();
        const V wm=w+(noisy ? V(noise(.019f)+V(.005f,-.004f,.008f)) : V::Zero());
        // A causal control only: supply the exact omitted rigid-body force
        // term, keeping every estimator correction active. This uses the
        // specified physical pivot and truth, not an inferred hardware pivot.
        const V rotational_force=alpha.cross(lever)+w.cross(w.cross(lever));
        const V am=a-(known_force_control ? rotational_force : V::Zero())+
            (noisy ? V(noise(.12f)) : V::Zero());
        const V mm=mag_history.front()+(noisy ? V(noise(.5f)) : V::Zero());
        t.tick(dt,truth,wm,am,mm,i%8==0,i%2==0);
        if (i==2200) { inherited_error=t.tail; inherited_bias=t.b->bias().norm(); }
    }
    if (recovery_diagnosis_mode)
        std::printf("%d,%d,%d,%.5f,%.6f,%.5f,%.6f,%u,%d\n",RAPID_BACKEND,
            combined,protect_after_turn,inherited_error,inherited_bias,t.tail,
            t.b->bias().norm(),t.applied_recoveries,t.last_held);
    else if (diagnosis_mode)
        std::printf("%d,%.0f,%.2f,%d,%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.6f,%u\n",
            RAPID_BACKEND,rate,radius,combined,gap,mag_delay*5,noisy,
            t.peak,t.tail,t.peak_heading,t.tail_heading,t.peak_tilt,t.tail_tilt,t.max_bias,t.applied_recoveries);
    else
        std::printf("rapid,%d,%d,%.0f,%.2f,%d,%d,%d,%.5f,%.5f,%.6f,%u,%d\n",RAPID_BACKEND,enabled,rate,radius,combined,gap,mag_delay*5,t.peak,t.tail,t.max_bias,t.applied_recoveries,noisy);
    if(enabled) {
        if(t.peak>=12) std::fprintf(stderr,"peak at %.3fs held=%d\n",float(t.peak_time)*1e-6f,t.peak_held);
        require(t.peak<12,"rapid-motion attitude excursion exceeds 12 degrees");
        require(t.tail<5,"rapid-motion recovery exceeds 5 degrees after 7 seconds");
        require(t.max_bias<.06f,"dynamic acceleration polluted gyro bias");
    }
    if(protect_after_turn) {
        require(t.tail<5,"inherited-state recovery exceeds five degrees");
        require(t.applied_recoveries>=1,"inherited-state recovery never captured");
        require(!t.last_held,"inherited-state recovery left aiding held");
        if(!recovery_diagnosis_mode)
            std::printf("recovery,%d,%d,%d,%d,%.2f,%.5f,%u,%d\n",RAPID_BACKEND,combined,noisy,moving_tail,float(steps-2200)*.005f,t.tail,t.applied_recoveries,calibrated_delay);
    }
}
/** @brief Inject a heading error in degrees and verify one recovery within three seconds. */
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
    require(t.protection.recoveryCount()<=1,"large heading error caused repeated reacquisition");
    require(t.tail<3,"large-angle capture outside 3 degrees after 3 seconds");
    require((t.b->bias()-bias_before).norm()<.025f,"large heading error became gyro bias");
}
/** @brief Verify that stale, invalid, or disturbed magnetic fields cannot trigger recovery. */
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
/** @brief Verify recovery requires sustained fresh acceleration, including across clock wrap. */
static void accelerometerDropout() {
    for (int mode=0; mode<3; ++mode) {
        Trial t;
        t.now=0xffe00000u; // Exercise source-age arithmetic across micros() wrap.
        const V gravity_body(0,0,-gravity);
        if (mode!=2)
            for (int i=0;i<400;++i)
                t.tick(.005f,Q::Identity(),V::Zero(),gravity_body,earth,i%8==0,i%2==0);
        t.b->recover(Q(Eigen::AngleAxisf(float(M_PI),V::UnitZ())));
        if (mode==0) {
            // A last measured force can have gravity's norm and a false tilt.
            const V contaminated=rotation(V(0,10*rad,0))*gravity_body;
            t.tick(.005f,Q::Identity(),V::Zero(),contaminated,earth,true,true);
        }
        // Establish the dropout before the next new magnetic observation;
        // a still-fresh pre-dropout gravity sample is valid recovery evidence.
        for (int i=0;i<6;++i)
            t.tick(.005f,Q::Identity(),V::Zero(),gravity_body,earth,false,false);
        for (int i=0;i<600;++i) {
            // No acceleration; isolated frames; or no initial acceleration.
            const bool fresh_acc=mode==1 && i%40==0;
            t.tick(.005f,Q::Identity(),V::Zero(),gravity_body,earth,i%8==0,fresh_acc);
            require(t.protection.recoveryCount()==0,"missing acceleration caused global recovery");
            if (i>4) require(t.last_held,"accelerometer dropout resumed attitude/bias aiding");
        }
        for (int i=0;i<600;++i) {
            // The 100 Hz acceleration stream need not coincide with 25 Hz mag.
            t.tick(.005f,Q::Identity(),V::Zero(),gravity_body,earth,i%8==0,i%2==1);
            if (i<140)
                require(t.protection.recoveryCount()==0,"recovery bridged an acceleration gap");
        }
        require(t.protection.recoveryCount()==1,"fresh acceleration did not restore recovery");
        require(t.tail<3,"dropout recovery installed false tilt or heading");
        std::printf("accel-dropout,%d,%d,%.5f,%u\n",RAPID_BACKEND,mode,t.tail,t.protection.recoveryCount());
    }
}
/** @brief Verify aiding resumes during moderate vessel motion after a rapid turn. */
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
/** @brief Bound an unsuccessful capture hold and retain later qualified recovery. */
static void pendingCaptureDuringMotion() {
    ocean_imu::attitude::Protection policy;
    Q truth=Q::Identity();
    const Q discrepancy(Eigen::AngleAxisf(.8f,V::UnitZ()));
    uint32_t now=1000;
    bool held=false, resumed=false, captured=false;
    policy.step(.005f,V(0,0,.18f),V::Zero(),V(0,0,-gravity),true,
        gravity,truth,earth,true,true,now);
    for(int i=1;i<=2000;++i) {
        now+=5000;
        // Above the capture rate limit, below every dynamic-motion threshold.
        // A capture cannot complete during these five seconds of valid aiding.
        const float rate=i<=1000 ? .18f : std::max(0.0f,.18f*float(1100-i)/100.0f);
        const V w(0,0,rate);
        truth=(truth*rotation(V(w*.005f))).normalized();
        const auto d=policy.step(.005f,w,V::Zero(),truth.conjugate()*V(0,0,-gravity),
            i%2==0,gravity,discrepancy*truth,truth.conjugate()*earth,true,i%8==0,now);
        require(!d.dynamic,"ordinary moving capture was classified as hazardous");
        if(i<=1000) {
            require(!d.reacquire,"unsettled observations installed an attitude");
            held=held||d.hold_attitude;
            if(i>500) {require(!d.hold_attitude,"pending capture starved aiding beyond two seconds");resumed=true;}
        } else if(d.reacquire) {
            require(d.target.angularDistance(truth)<1e-5f,"late settled capture changed the north frame");
            captured=true; break;
        }
    }
    require(held&&resumed&&captured,"bounded hold did not retain qualified recovery");
}
/** @brief A wave-force hold must retain informative magnetic heading feedback without learning bias. */
static void magneticDuringForceHold() {
    Backend aided, absent;
    const Q initial(Eigen::AngleAxisf(.30f,V::UnitZ()));
    aided.recover(initial); absent.recover(initial);
    const V initial_bias=aided.bias();
    const V slammed_force(7,-3,-gravity);
    for(int i=1;i<=400;++i) {
        aided.step(.005f,V::Zero(),slammed_force,earth,true,i%8==0,false,i%2==0);
        absent.step(.005f,V::Zero(),slammed_force,earth,false,false,false,i%2==0);
        require((aided.bias()-initial_bias).norm()<2e-6f,"held magnetic feedback changed learned bias");
        const V down=aided.attitude().conjugate()*V::UnitZ();
        require((down-V::UnitZ()).norm()<1e-4f,"held magnetic feedback tilted the gravity estimate");
        require(aided.finite(),"held magnetic feedback produced nonfinite state");
    }
    std::printf("held-mag,%d,%.5f\n",RAPID_BACKEND,aided.attitude().angularDistance(Q::Identity())/rad);
    // This checks informative service, not a common convergence rate for
    // seven different observers/covariance priors. The absent-field arm has
    // zero correction; require a resolved (>0.01 rad) heading improvement.
    require(aided.attitude().angularDistance(Q::Identity()) < initial.angularDistance(Q::Identity())-.01f,
        "force hold blocked informative magnetic heading corrections");
    require(absent.attitude().angularDistance(initial)<1e-4f,"held force changed attitude without a field");
}
/** @brief Check single-use acceleration freshness and elapsed time across clock wrap and gaps. */
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
/** @brief Distinguish corrupt bias, healthy bias and observed slow motion. */
static void biasRecoveryPolicy() {
    for (bool moving:{false,true}) for (bool noisy:{false,true}) for (int mode=0;mode<3;++mode) {
        ocean_imu::attitude::Protection policy;
        const V actual_bias=mode==2 ? V(0,0,.095f) : V(.008f,-.006f,.012f);
        const V learned=mode==0 ? V(.3f,-.12f,.09f) : actual_bias;
        const V rate=moving && mode!=2 ? V(.012f,.006f,.018f) : V::Zero();
        Q truth=Q(Eigen::AngleAxisf(.3f,V::UnitY()));
        std::mt19937 rng(913); std::normal_distribution<float> normal(0,1);
        auto noise=[&](float sigma)->V {return sigma*V(normal(rng),normal(rng),normal(rng));};
        bool captured=false;
        for (int i=1;i<=2000 && !captured;++i) {
            truth=(truth*rotation(V(rate*.005f))).normalized();
            const V w=rate+actual_bias+(noisy?noise(.003f):V::Zero());
            const V a=truth.conjugate()*V(0,0,-gravity)+(noisy?noise(.08f):V::Zero());
            const V m=truth.conjugate()*earth+(noisy?noise(.3f):V::Zero());
            const Q wrong=Q(Eigen::AngleAxisf(1.2f,V::UnitZ()))*truth;
            const auto d=policy.step(.005f,w,learned,a,true,gravity,wrong,m,true,i%8==0,uint32_t(i*5000));
            if (!d.reacquire) { require(!d.repair_gyro_bias,"bias repair without qualified capture"); continue; }
            captured=true;
            if(d.repair_gyro_bias!=(mode==0)) std::fprintf(stderr,"bias-policy mode=%d moving=%d noisy=%d sample=%d repair=%d\n",mode,moving,noisy,i,d.repair_gyro_bias);
            require(d.repair_gyro_bias==(mode==0),"healthy bias replaced or corrupt bias not repaired");
            if (d.repair_gyro_bias && (d.gyro_bias-actual_bias).norm()>=(noisy?.04f:.003f)) std::fprintf(stderr,"bias-policy mode=%d moving=%d noisy=%d sample=%d bias=%g,%g,%g error=%g\n",mode,moving,noisy,i,d.gyro_bias.x(),d.gyro_bias.y(),d.gyro_bias.z(),(d.gyro_bias-actual_bias).norm());
            if (d.repair_gyro_bias)
                require((d.gyro_bias-actual_bias).norm()<(noisy?.04f:.003f),"observed motion was misidentified as gyro bias");
        }
        if(!captured) std::fprintf(stderr,"no bias capture mode=%d moving=%d noisy=%d\n",mode,moving,noisy);
        require(captured,"bias policy did not obtain a qualified window");
    }
}
/** @brief Compare local aiding and protected recovery at an exact horizontal-field antipode. */
static void antipodalInnovation() {
    const V saved_field=earth;
    earth=V(50,0,0);
    // Exact antipode, without sin(pi)'s floating-point perturbation. Gravity
    // and field are perpendicular: each local vector residual has zero
    // tangential correction even though the heading is exactly wrong.
    const Q wrong(0,0,0,1);
    for (bool protected_run:{false,true}) {
        Trial t; t.enabled=protected_run;
        t.b->recover(wrong);
        for (int i=0;i<600;++i)
            t.tick(.005f,Q::Identity(),V::Zero(),V(0,0,-gravity),earth,i%8==0,i%2==0);
        std::printf("antipode,%d,%d,%.6f,%u\n",RAPID_BACKEND,protected_run,t.tail,t.applied_recoveries);
        if (protected_run) {
            require(t.applied_recoveries==1 && t.tail<3,"global capture did not escape zero-tangent antipode");
        } else {
            require(t.applied_recoveries==0 && t.tail>179,"local antipode control did not retain its wrong equilibrium");
        }
    }
    earth=saved_field;
}
#if RAPID_BACKEND >= 3 && RAPID_BACKEND <= 5
/** @brief Check held gain rows, covariance PSD, and navigation/bias preservation on recovery. */
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
    const V mag_pos=f.get_position(),mag_vel=f.get_velocity(),mag_aw=f.get_world_accel();
    const V mag_field=Q(Eigen::AngleAxisf(.25f,V::UnitZ())).conjugate()*earth;
#if RAPID_BACKEND != 3
    const V mag_integ=f.get_integral_displacement();
#endif
    f.measurement_update_mag_only(mag_field);
    require(f.lastMagDiag().accepted,"held core rejected a finite informative magnetic correction");
    require((mag_pos-f.get_position()).norm()<1e-5f && (mag_vel-f.get_velocity()).norm()<1e-5f &&
        (mag_aw-f.get_world_accel()).norm()<1e-5f,"held magnetic correction changed physical navigation means");
#if RAPID_BACKEND != 3
    require((mag_integ-f.get_integral_displacement()).norm()<1e-5f,"held magnetic correction changed displacement history");
#endif
    require((b->bias()-bg).norm()<1e-6f && (f.get_acc_bias()-ba).norm()<1e-6f,"held magnetic correction changed bias");
    psd();
    const auto physical_covariance=[&]() -> Eigen::MatrixXf {
        Eigen::MatrixXf P=f.covariance_full();
#if RAPID_BACKEND == 5
        Eigen::MatrixXf J=Eigen::MatrixXf::Identity(P.rows(),P.cols());
        const V state[4]={f.get_velocity(),f.get_position(),f.get_integral_displacement(),f.get_world_accel()};
        for (int c=0;c<4;++c) J.block<3,3>(6+3*c,0)=-ocean_imu::lie::skew<float>(state[c]);
        J.block<3,3>(3,3)=f.R_wb(); J.block<3,3>(18,18)=f.R_wb();
        P=(J*P*J.transpose()).eval();
#endif
        return P;
    };
    const auto marginal=physical_covariance();
    f.set_attitude_aiding_enabled(true);
    require((physical_covariance().diagonal()-marginal.diagonal()).norm()<1e-6f,"resume changed physical covariance marginals");
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

/** @brief Check startup and Live facade observers obey the attitude-aiding hold. */
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
/** @brief Run protected regression scenarios or emit ordinary-aiding diagnostic controls. */
int main(int argc, char** argv) {
    if (argc==2 && std::strcmp(argv[1],"--diagnose-recovery")==0) {
        recovery_diagnosis_mode=true;
        std::puts("backend,combined,protection_after_turn,entry_attitude_deg,entry_bias_rad_s,tail_attitude_deg,tail_bias_rad_s,applied_recoveries,tail_held");
        for (bool combined:{false,true}) for (bool late:{false,true})
            rapid(720,.15f,combined,false,false,0,false,late);
        return 0; // The late-protection arms above enforce the recovery target.
    }
    // Causal controls: every ordinary correction remains enabled, with exact
    // delivered gyro/field. Vary only the physical lever arm/rate or one named
    // timing/noise condition. No guard or reacquisition acts on the filter.
    if (argc==2 && (std::strcmp(argv[1],"--diagnose")==0 ||
                   std::strcmp(argv[1],"--diagnose-force-equation")==0)) {
        diagnosis_mode=true;
        known_force_control=std::strcmp(argv[1],"--diagnose-force-equation")==0;
        std::puts("backend,rate_deg_s,radius_m,combined,gap_75ms,mag_delay_ms,noisy,peak_attitude_deg,tail_attitude_deg,peak_heading_deg,tail_heading_deg,peak_tilt_deg,tail_tilt_deg,max_bias_rad_s,applied_recoveries");
        for (float rate:{180.0f,360.0f,720.0f}) for (float r:{0.0f,.05f,.15f})
            for (bool combined:{false,true}) rapid(rate,r,combined,false,false);
        rapid(720,0,false,true,false); rapid(720,.15f,false,true,false);
        rapid(720,0,true,false,false,8); rapid(720,.15f,true,false,false,8);
        rapid(720,0,false,false,false,0,true); rapid(720,.15f,false,false,false,0,true);
        return 0;
    }
    acquisition(); biasRecoveryPolicy(); antipodalInnovation(); pendingCaptureDuringMotion(); magneticDuringForceHold();
    for(float rate:{180.0f,360.0f,720.0f}) for(float r:{0.0f,.05f,.15f})
        for(bool combined:{false,true}) rapid(rate,r,combined);
    rapid(720,.15f,false,true); rapid(720,.15f,true,false,true,8,false,false,false,3600,true);
    rapid(720,.15f,false,false,true,0,true);
    rapid(720,.15f,true,false,true,8,true,false,false,3600,true);
    // Matched unprotected control: perfect sensors, no lag, same acceleration.
    rapid(720,.15f,false,false,false);
    for(bool combined:{false,true}) {
        rapid(720,.15f,combined,false,false,0,false,true);
#if RAPID_BACKEND == 1 || RAPID_BACKEND == 5
        rapid(720,.15f,combined,false,false,0,true,true);
        rapid(720,.15f,combined,false,false,0,true,true,true,6200);
        rapid(720,.15f,combined,false,false,0,false,true,false,14200);
        // Main's calibrated relative lag must coexist with corrupted-state capture.
        rapid(720,.15f,combined,false,false,8,true,true,false,3600,true);
#endif
    }
    for(float a:{90.0f,150.0f,179.99f,180.0f,180.01f,-179.99f}) for(float h:{0.0f,137.0f}) capture(a,h);
    for(float h:{0.0f,137.0f}) capture(180,h,true);
    staleAndDisturbance(); accelerometerDropout(); movingAfterTurn();
#if RAPID_BACKEND >= 3 && RAPID_BACKEND <= 5
    selectiveCovarianceAndRecovery();
    facadeProtection();
#endif
    std::fprintf(stderr,"Rapid-motion backend %d PASS\n",RAPID_BACKEND);
}
