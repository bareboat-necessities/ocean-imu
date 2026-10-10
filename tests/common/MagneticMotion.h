#pragma once
#include "util/MagneticRotation.h"
#include <array>
#include <vector>
#include <memory>
#include <optional>
#include <random>
#include <cstdio>
#include <cstdlib>
#include <string>
using namespace ocean_imu::magnetic;
const float g_std=9.80665f;
static const V3 kField(20.3f,0,46.8f), hard(8,-6,3);
static const M3 rm=.64f*M3::Identity();
static void check(bool ok,const char* what) {
    if(!ok) {std::fprintf(stderr,"FAIL: %s\n",what); std::exit(1);}
}
struct Row {
    uint32_t time=0; float dt=0;
    V3 gyro=V3::Zero(),acc=V3::Zero(),mag=V3::Zero();
    Q truth=Q::Identity();
    bool delivered=false, sample=false;
    int arrival=0;
};
static float heading(const Q& q) {
    const V3 forward=q*V3::UnitX(); return std::atan2(forward.y(),forward.x());
}
static double error(const Q& q,const Q& truth) {
    return double(std::remainder(heading(q)-heading(truth),2*float(M_PI)))*180/M_PI;
}
static std::vector<Row> magneticMotion(int delay_ms,int motion,bool jitter,bool bias,bool soft,bool interference) {
    constexpr int N=5200;
    std::vector<Row> rows(N+1);
    std::mt19937 rng(45321); std::normal_distribution<float> normal(0,1);
    Q truth=Q::Identity(); uint32_t t=0,next_mag=40000;
    for(int i=1;i<=N;++i) {
        Row& r=rows[size_t(i)];
        const uint32_t us=jitter ? uint32_t(4000+(i*773)%2001) : 5000u;
        t+=us; r.time=t; r.dt=float(us)*1e-6f; const float s=float(t)*1e-6f;
        V3 w=V3::Zero();
        if(motion==1) w=V3(0,0,.8f);
        if(motion==2) w=V3(.3f*std::cos(.6f*s),.25f*std::sin(.5f*s),.55f+.25f*std::sin(.3f*s));
        const float a=w.norm()*r.dt;
        if(a>0) truth=truth*Q(Eigen::AngleAxisf(a,w.normalized()));
        truth.normalize(); r.truth=truth;
        r.gyro=w+(.00135f/std::sqrt(r.dt))*V3(normal(rng),normal(rng),normal(rng));
        if(bias) r.gyro+=V3(.003f,-.002f,.005f);
        r.acc=truth.conjugate()*V3(0,0,-g_std)+.02f*V3(normal(rng),normal(rng),normal(rng));
        r.mag=truth.conjugate()*kField+.15f*V3(normal(rng),normal(rng),normal(rng));
        if(soft) r.mag=V3(1.08f,.94f,1.03f).asDiagonal()*r.mag;
        if(interference&&s>10&&s<15) r.mag+=V3(20,-15,5);
        r.mag+=hard; // known additional body offset, identical for all paths
        if(t>=next_mag) {
            r.sample=true; next_mag+=jitter ? uint32_t(33333+(i*379)%4001) : 40000;
        }
    }
    // Delivery delay/jitter alters arrival only; all three paths use exactly
    // the same physical samples and noise, and the oracle sees no future data.
    for(int i=1;i<=N;++i) if(rows[size_t(i)].sample) {
        auto& r=rows[size_t(i)];
        const uint32_t delay=uint32_t(delay_ms)*1000u+(jitter ? uint32_t((i*137)%5001) : 0u);
        int k=i; while(k<=N&&rows[size_t(k)].time<r.time+delay) ++k;
        r.arrival=k;
    }
    return rows;
}
