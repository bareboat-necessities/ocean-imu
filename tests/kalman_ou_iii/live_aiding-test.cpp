// Complete startup/runtime replay with the same per-frame ordering as the OU3
// sketch. Truth only generates sensors; no estimator state is injected.
#define EIGEN_NON_ARDUINO
extern const float g_std=9.80665f;
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "util/MagneticInput.h"
#include "util/MagneticHeading.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <memory>
#include <random>

namespace {
using V=Eigen::Vector3f;
using Q=Eigen::Quaternionf;
using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;
constexpr float gravity=9.8025605f, dt=.005f, degrees=57.295779513f;
const V earth(20.3f,0,46.8f);
void require(bool ok,const char* reason) {
    if(!ok) {std::fprintf(stderr,"OU3 Live aiding: %s\n",reason);std::exit(1);}
}
void replay(float roll,float pitch) {
    auto fusion=std::make_unique<Fusion>();
    Fusion::Config cfg;
    cfg.gravity_magnitude=gravity; cfg.mag_delay_sec=0;
    cfg.sigma_a=V::Constant(.12f); cfg.sigma_g=V::Constant(.00135f);
    cfg.sigma_m=V::Constant(.8f); cfg.mag_continuous_hard_iron=false;
    fusion->begin(cfg); fusion->raw().setTauScaledPseudoUpdateCadence(false);
    ocean_imu::magnetic::Input magnetic;
    magnetic.reset(.00135f,Eigen::Matrix3f::Identity()*.64f);
    atoms3r_ical::MagAcquisition timing;
    timing.delay_us=40000; timing.timing_allowance_s=.03f;
    ocean_imu::ins::MagFreshGate gate;
    Q truth=Q(Eigen::AngleAxisf(1.7f,V::UnitZ()))*
        Q(Eigen::AngleAxisf(pitch,V::UnitY()))*Q(Eigen::AngleAxisf(roll,V::UnitX()));
    V field=truth.conjugate()*earth, accel=truth.conjugate()*V(0,0,-gravity);
    V previous_rate=V::Zero();
    std::deque<V> fields;
    std::mt19937 rng(643); std::normal_distribution<float> normal(0,1);
    auto noise=[&](float sigma)->V {return sigma*V(normal(rng),normal(rng),normal(rng));};
    float max_heave=0, tail_error=0;
    bool live_seen=false; unsigned tail_mag_updates=0;
    // Continued slow motion forbids a settled capture. Include large physical
    // heel as on a handheld device; this is an engineering test, not a vessel
    // envelope or a regional stability theorem.
    for(int i=1;i<=90000;++i) {
        const float t=static_cast<float>(i)*dt;
        const V rate(.2f*std::cos(.8f*t),.13f*std::sin(.8f*t),0);
        truth=(truth*Q(Eigen::AngleAxisf(rate.norm()*dt,rate.normalized()))).normalized();
        const V lever(.15f,0,0);
        V specific_force=truth.conjugate()*V(0,0,-gravity)+
            (rate-previous_rate).cross(lever)/dt+rate.cross(rate.cross(lever));
        previous_rate=rate;
        // Bounded translation with displacement -.3/.8^2*sin(.8*t) on
        // each driven world axis; no DC acceleration/displacement ramp.
        specific_force+=truth.conjugate()*V(.3f*std::sin(.8f*t),0,.3f*std::sin(.8f*t));
        fields.push_back(truth.conjugate()*earth);
        if(fields.size()>9) fields.pop_front();
        if(i%8==0) field=fields.front()+noise(.5f);
        if(i%2==0) accel=specific_force+noise(.12f);
        const V gyro=rate+noise(.003f)+V(.002f,-.003f,.005f);
        const uint32_t now=1000000u+static_cast<uint32_t>(i)*5000u;
        if(i%8==0) timing.acquired(now/1000,now,now,now);
        if(i%2!=0) continue;
        const float frame_dt=i==2?dt:2*dt;
        auto& core=fusion->raw().mekf();
        const V measured=accel;
        magnetic.advance(now,gyro-core.gyroscope_bias_body(),
            ocean_imu::magnetic::covarianceBound(core.gyroscope_bias_covariance_body()),
            timing,field,fusion->magHardIronBodyUT());
        V current_field=V::Zero();
        const bool field_valid=magnetic.currentField(fusion->magHardIronBodyUT(),current_field);
        const bool due=magnetic.due(gate,true,now/1000);
        fusion->update(frame_dt,gyro,measured);
        magnetic.capture(fusion->raw().startupProxyQuat(),fusion->attitudeQuat(),measured,gyro);
        if(due) {
            magnetic.correct(*fusion,now);
            if(t>420&&magnetic.lastApplied()) ++tail_mag_updates;
        }
        if(!fusion->isLive()) continue;
        live_seen=true;
        require(core.covariance_full().allFinite(),"nonfinite covariance");
        require(fusion->attitudeQuat().coeffs().allFinite(),"nonfinite attitude");
        max_heave=std::max(max_heave,std::abs(core.get_position().z()));
        tail_error=fusion->attitudeQuat().angularDistance(truth)*degrees;
    }
    std::printf("OU3 live-aiding roll=%.3f pitch=%.3f tail_deg=%.3f raw_heave_m=%.3f tail_mag=%u\n",
        double(roll),double(pitch),double(tail_error),double(max_heave),tail_mag_updates);
    require(live_seen,"never entered Live");
    require(tail_mag_updates>100,"continued motion lost magnetic service");
    require(tail_error<10,"continued motion lost magnetic north");
    require(max_heave<2.5f,"continued motion caused excessive false raw heave");
}
}
int main() {replay(1.3f,.2f);replay(3.14159265f,0);}
