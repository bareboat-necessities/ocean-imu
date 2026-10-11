#define EIGEN_NON_ARDUINO
#define EIGEN_RUNTIME_NO_MALLOC
#include "util/MagneticInput.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cstdio>
#include <cstdlib>
using namespace ocean_imu::magnetic;
static void check(bool ok, const char* message) {
    if (!ok) { std::fprintf(stderr, "FAIL %s\n", message); std::exit(1); }
}
struct Gate {
    uint32_t last = 0;
    bool update(bool valid, uint32_t t) {
        if (!valid || (last && t-last < 35)) return false;
        last=t; return true;
    }
};
struct Fusion {
    int calls = 0;
    void updateMag(const V3&) { ++calls; }
    void updateMagTransported(const V3&, const Observation& o) { check(o.valid(), "valid snapshot"); ++calls; }
    bool lastMagUpdateApplied() const { return true; }
};
int main() {
    Eigen::internal::set_is_malloc_allowed(false);
    Input input; Gate gate; Fusion fusion;
    const V3 field(20,0,45), hard(8,-6,3), rate(.3f,-.1f,.8f);
    input.reset(.00135f, .64f*M3::Identity());
    atoms3r_ical::MagAcquisition timing;
    timing.acquired(1,1000,1000,1100);
    input.advance(1000,rate,1e-6f,timing,field+hard,hard);
    input.capture(Q::Identity(),Q::Identity(),V3(0,0,-9.80665f),rate);
    check(input.due(gate,true,1), "first discrete observation due");
    input.correct(fusion,1200);
    for(uint32_t t=6000;t<=81000;t+=5000) {
        input.advance(t,rate,1e-6f,timing,field+hard,hard);
        check(!input.due(gate,true,t/1000), "held observation is not another discrete correction");
        V3 aligned; M3 covariance;
        check(input.currentField(hard,aligned,&covariance), "continuous feedback retained between samples");
        check((aligned-increment(rate,float(t-1000)*1e-6f)*field).norm()<3e-5f,
            "continuous hold rotates field after hard-iron subtraction");
    }
    check(fusion.calls==1, "one discrete call per source observation");
    // New observation arrives before spacing gate; save this exact proxy,
    // and do not overwrite it with the later delivery-time attitude.
    timing.acquired(82,82000,82000,82100);
    input.advance(82000,rate,1e-6f,timing,field+hard,hard);
    const Q proxy(Eigen::AngleAxisf(.2f,V3::UnitX()));
    input.capture(proxy,proxy,V3(0,0,-9.8f),rate);
    input.advance(87000,rate,1e-6f,timing,field+hard,hard);
    input.capture(Q::Identity(),Q::Identity(),V3::Zero(),V3::Zero());
    check(input.prepare(), "cached observation interval exists");
    check(input.observation().proxy_bw.angularDistance(proxy)<1e-6f,"snapshot stays at observation epoch");
    input.advance(87000,rate,1e-6f,timing,field+hard,hard);
    check(!input.prepare(),"duplicate gyro epoch invalidates transport");
    input.advance(92000,rate,1e-6f,timing,field+hard,hard);
    check(!input.prepare(),"next step cannot bridge duplicate epoch");
    timing.acquired(97,97000,97000,97100);
    input.advance(97000,rate,0,timing,field+hard,hard,false);
    V3 aligned;
    check(input.currentField(hard,aligned),"continuous recovery on new observation");
    check(std::isnan(input.diagnostic().nis),"observer without covariance does not claim calibrated NIS");
    timing.have=false;
    input.advance(102000,rate,0,timing,field+hard,hard,false);
    check(!input.currentField(hard,aligned),"missing source unavailable");
    timing.have=true; timing.sequence--;
    input.advance(107000,rate,0,timing,field+hard,hard,false);
    check(!input.currentField(hard,aligned),"out-of-order source unavailable");
    timing.sequence += 2; timing.frame_us = 96000;
    input.advance(112000,rate,0,timing,field+hard,hard,false);
    check(!input.currentField(hard,aligned),"new sequence cannot legitimize an older timestamp");
    timing.frame_us = 97000;
    input.advance(117000,rate,0,timing,field+hard,hard,false);
    check(!input.currentField(hard,aligned),"new sequence cannot legitimize a duplicate timestamp");
    RotationHistory<> discrete_noise;
    discrete_noise.push(0,rate);
    discrete_noise.push(5000,rate,0,9e-6f);
    discrete_noise.push(10000,rate,0,9e-6f);
    RotationHistory<>::Rotation rotation;
    check(discrete_noise.between(2500,7500,rotation)==Status::Ready,"partial discrete-noise interval");
    check(std::abs(rotation.gyro_angle_variance-4.5e-6f)<1e-10f,
        "fractional rotation scales the same sampled angle error before covariance addition");
    check(discrete_noise.push(15000,rate,0,NAN)==Status::Invalid,"invalid step covariance rejected");
    std::printf("Shared magnetic input PASS; sizeof(Input)=%zu\n",sizeof(Input));
}
