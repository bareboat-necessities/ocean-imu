#pragma once
#include "MagneticMotion.h"
#include "util/MagneticInput.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"

// Backend is the actual continuous observer plus thin axis/API adapters.
// Replay snapshots include its complete state AND the magnetic hold/history.
// Only measurements that have already arrived are inserted during rewind.
struct ContinuousMachine {
    Backend observer{};
    Input input{};
    atoms3r_ical::MagAcquisition timing{};
    V3 held=V3::Zero();
    ContinuousMachine() { input.reset(0,M3::Identity()); }
    bool step(const Row& row,const Row* magnetic,bool compensate,V3& used) {
        if(magnetic) {
            held=magnetic->mag;
            timing.acquired(row.time/1000,magnetic->time,magnetic->time,magnetic->time);
        }
        input.advance(row.time,row.gyro-observer.bias(),0,timing,held,hard,false);
        bool valid=timing.have;
        used=held-hard;
        if(valid&&compensate) valid=input.currentField(hard,used);
        observer.step(row.dt,row.gyro,row.acc,used,valid);
        check(observer.finite(),"continuous observer remains finite");
        return valid;
    }
};
struct ContinuousScore {
    double ss=0,peak=0,tail=0,vector=0,reference=0;
    int n=0,ntail=0,nvector=0,offered=0,missing=0;
    void add(const ContinuousMachine& f,const Row& row,const Q& oracle,const V3& used,bool valid) {
        if(f.timing.have) { ++offered; if(!valid) ++missing; }
        if(row.time<3000000) return;
        const double e=error(f.observer.attitude(),row.truth);
        const double r=error(f.observer.attitude(),oracle);
        ss+=e*e; reference+=r*r; peak=std::max(peak,std::abs(e)); ++n;
        if(row.time>22000000) { tail+=e*e; ++ntail; }
        if(valid) {vector+=(used-row.truth.conjugate()*kField).squaredNorm(); ++nvector;}
    }
    void print(const char* name,int delay,const char* mode) const {
        std::printf("%s,%d,%s,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",name,delay,mode,
            std::sqrt(vector/std::max(1,nvector)),std::sqrt(ss/n),peak,
            double(missing)/std::max(1,offered),std::sqrt(tail/std::max(1,ntail)),std::sqrt(reference/n));
    }
};
static void continuousPaired(const char* name,int delay,int motion,bool jitter=false,bool bias=false,
                             bool soft=false,bool interference=false) {
    auto rows=magneticMotion(delay,motion,jitter,bias,soft,interference);
    auto legacy=std::make_unique<ContinuousMachine>(), aligned=std::make_unique<ContinuousMachine>();
    auto oracle=std::make_unique<ContinuousMachine>();
    std::vector<std::optional<ContinuousMachine>> snapshots(40);
    ContinuousScore score[3];
    int next=1;
    for(int i=1;i<int(rows.size());++i) {
        auto& row=rows[size_t(i)];
        const Row* arrival=nullptr;
        while(next<=i&&!rows[size_t(next)].sample) ++next;
        if(next<=i&&rows[size_t(next)].arrival<=i) {
            arrival=&rows[size_t(next)]; rows[size_t(next)].delivered=true;
        }
        V3 used[3]; bool valid[3];
        valid[0]=legacy->step(row,arrival,false,used[0]);
        valid[1]=aligned->step(row,arrival,true,used[1]);
        snapshots[size_t(i%40)].emplace(*oracle);
        valid[2]=oracle->step(row,nullptr,true,used[2]);
        if(arrival) {
            check(i-next<40,"continuous oracle rewind fits history");
            oracle=std::make_unique<ContinuousMachine>(*snapshots[size_t(next%40)]);
            for(int j=next;j<=i;++j) {
                snapshots[size_t(j%40)].emplace(*oracle);
                const auto& r=rows[size_t(j)];
                valid[2]=oracle->step(r,r.delivered?&r:nullptr,true,used[2]);
            }
            ++next;
        }
        const Q reference=oracle->observer.attitude();
        score[0].add(*legacy,row,reference,used[0],valid[0]);
        score[1].add(*aligned,row,reference,used[1],valid[1]);
        score[2].add(*oracle,row,reference,used[2],valid[2]);
    }
    score[0].print(name,delay,"shipping_delayed_hold");
    score[1].print(name,delay,"aligned_hold");
    score[2].print(name,delay,"timestamp_replay");
    for(const auto& s:score) check(s.missing==0,"valid continuous magnetic feedback is retained");
    if(motion&&!bias&&!soft&&!interference)
        check(score[1].vector<score[0].vector,"gyro transport improves held vector chronology");
    // Heading/recovery are reported, not assumed to inherit OU-III's tuning
    // or convergence rate. The underlying observer and its gains are unchanged.
}
static void runContinuousMagneticTests() {
    std::puts("scenario,delay_ms,mode,vector_rms_uT,heading_rms_deg,heading_peak_deg,unavailable_feedback_rate,tail_heading_rms_deg,heading_vs_replay_rms_deg");
    for(int d:{5,10,20,40}) {continuousPaired("yaw",d,1);continuousPaired("combined_variable",d,2);}
    continuousPaired("async_jitter",20,2,true);
    continuousPaired("gyro_bias_error",40,2,true,true);
    continuousPaired("soft_iron_residual",20,2,false,false,true);
    continuousPaired("interference_recovery",20,1,false,false,false,true);
    continuousPaired("quiet",40,0,true);
}
