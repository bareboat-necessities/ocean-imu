#pragma once
// Shared deterministic motion/noise, arrival chronology and exact rewind/replay
// for the actual OU-II, OU-III and TFG cores. Each small driver supplies only
// its core construction, attitude accessor and literal prediction chronology.
#include "MagneticMotion.h"
static Core core();
static Q attitude(const Core&);
static void advanceCore(Core&,const V3&,const V3&,float,int);
template <class C = Core>
static M3 sampleCovariance(const V3& field) {
    if constexpr (requires { C::MagneticUnitSigma; })
        return C::MagneticUnitSigma*C::MagneticUnitSigma*field.squaredNorm()*M3::Identity();
    else return rm;
}
template <class C = Core>
static float gyroStepVariance() {
    if constexpr (requires { C::GyroStepVariance; }) return C::GyroStepVariance;
    else return 0;
}
static void advance(Core& f,const Row& r,int index) {
    advanceCore(f,r.gyro,r.acc,r.dt,index);
}
struct Score {
    double ss=0,peak=0,vector_ss=0,tail_ss=0,reference_ss=0;
    int n=0,m=0,tail_n=0,accepted=0,attempted=0;
    void add(Core& f,const Row& row,const Q& reference) {
        const double e=error(attitude(f),row.truth);
        check(std::isfinite(e)&&f.covariance_full().allFinite(),"finite complete estimator state");
        if(row.time<3000000) return;
        ss+=e*e; peak=std::max(peak,std::abs(e)); ++n;
        const double diff=error(attitude(f),reference); reference_ss+=diff*diff;
        if(row.time>22000000) {tail_ss+=e*e;++tail_n;}
    }
    double rms() const {return std::sqrt(ss/n);}
    void print(const char* name,int delay,const char* mode) const {
        std::printf("%s,%d,%s,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",name,delay,mode,
            std::sqrt(vector_ss/std::max(1,m)),rms(),peak,
            1-double(accepted)/std::max(1,attempted),std::sqrt(tail_ss/std::max(1,tail_n)),
            std::sqrt(reference_ss/n));
    }
};
static void paired(const char* name,int delay_ms,int motion,bool jitter=false,bool bias=false,
                   bool soft=false,bool interference=false) {
    auto rows=magneticMotion(delay_ms,motion,jitter,bias,soft,interference);
    const int N=int(rows.size())-1;
    auto legacy=std::make_unique<Core>(core()), aligned=std::make_unique<Core>(core()), ref=std::make_unique<Core>(core());
    // Test-only rewind snapshots. 40 slots comfortably exceed 40 ms delay.
    std::vector<std::optional<Core>> snapshots(40);
    RotationHistory<> h; h.push(0,V3::Zero());
    Score scores[3]; RotationConsistency consistency;
    int next=1; double diagnostic_ss=0; int diagnostic_n=0;
    for(int i=1;i<=N;++i) {
        Row& r=rows[size_t(i)];
        h.push(r.time,r.gyro-aligned->gyroscope_bias_body(),
               covarianceBound(aligned->gyroscope_bias_covariance_body()),gyroStepVariance());
        snapshots[size_t(i%40)].emplace(*ref);
        advance(*legacy,r,i); advance(*aligned,r,i); advance(*ref,r,i);
        while(next<=i) {
            if(!rows[size_t(next)].sample) {++next;continue;}
            if(rows[size_t(next)].arrival>i) break;
            Row& m=rows[size_t(next)]; m.delivered=true;
            check(i-next<40,"oracle rewind fits snapshots");
            ref=std::make_unique<Core>(*snapshots[size_t(next%40)]);
            for(int j=next;j<=i;++j) {
                snapshots[size_t(j%40)].emplace(*ref);
                const auto& rr=rows[size_t(j)]; advance(*ref,rr,j);
                if(rr.delivered) ref->measurement_update_mag_only(rr.mag-hard);
            }
            RotationHistory<>::Rotation rotation;
            check(h.between(m.time,r.time,rotation)==Status::Ready,"delayed magnetic interval covered");
            Uncertainty u; u.gyro_density=gyroStepVariance()>0 ? 0 : .00135f;
            u.bias_variance=covarianceBound(aligned->gyroscope_bias_covariance_body());
            // Actual timestamps in this fixture are known; unknown hardware
            // conversion epochs are intentionally not imitated by invented data.
            V3 v; M3 covariance;
            check(align(m.mag,hard,sampleCovariance(m.mag-hard),rotation,u,v,covariance),"transported measurement covariance valid");
            legacy->measurement_update_mag_only(m.mag-hard);
            aligned->measurement_update_mag_only(v,covariance);
            const V3 exact_now=r.truth.conjugate()*kField;
            scores[0].vector_ss+=(m.mag-hard-exact_now).squaredNorm();
            scores[1].vector_ss+=(v-exact_now).squaredNorm();
            scores[2].vector_ss+=(m.mag-hard-m.truth.conjugate()*kField).squaredNorm();
            for(auto& score:scores) {++score.m;++score.attempted;}
            scores[0].accepted+=legacy->lastMagDiag().accepted;
            scores[1].accepted+=aligned->lastMagDiag().accepted;
            scores[2].accepted+=ref->lastMagDiag().accepted;
            const auto d=consistency.observe(h,uint32_t(next),m.time,m.mag,hard,sampleCovariance(m.mag-hard),u);
            if(d.status==Status::Ready) {diagnostic_ss+=d.nis;++diagnostic_n;}
            ++next;
        }
        scores[0].add(*legacy,r,attitude(*ref));
        scores[1].add(*aligned,r,attitude(*ref));
        scores[2].add(*ref,r,attitude(*ref));
    }
    scores[0].print(name,delay_ms,"shipping_delayed");
    scores[1].print(name,delay_ms,"aligned");
    scores[2].print(name,delay_ms,"timestamp_replay");
    std::fprintf(stderr,"%s %dms consistency mean d2 %.5f\n",name,delay_ms,diagnostic_ss/std::max(1,diagnostic_n));
    if(motion>0&&!bias&&!soft&&!interference) {
        check(scores[1].rms()<scores[0].rms(),"paired heading RMS improves for known delay");
        check(scores[1].vector_ss<scores[0].vector_ss,"paired vector residual improves");
        check(std::sqrt(scores[1].reference_ss/scores[1].n)<.3,"transport approximation near replay heading");
    }
    for(const auto& s:scores) check(s.accepted==s.attempted,"no valid magnetic update suppressed");
    if(interference && kRequireOu3Recovery) check(std::sqrt(scores[1].tail_ss/scores[1].tail_n)<2,"recovery after external kField removal");
}

static void runPairedMagneticTests() {
    std::puts("scenario,delay_ms,mode,vector_rms_uT,heading_rms_deg,heading_peak_deg,rejection_rate,tail_heading_rms_deg,heading_vs_replay_rms_deg");
    for(int d:{5,10,20,40}) {paired("yaw",d,1);paired("combined_variable",d,2);}
    paired("async_jitter",20,2,true);
    paired("gyro_bias_error",40,2,true,true);
    paired("soft_iron_residual",20,2,false,false,true);
    paired("interference_recovery",20,1,false,false,false,true);
    paired("quiet",40,0,true);
}
static void covarianceCompatibility() {
    auto a=std::make_unique<Core>(core()), b=std::make_unique<Core>(core());
    a->measurement_update_mag_only(kField+V3(.2f,-.1f,.1f));
    b->measurement_update_mag_only(kField+V3(.2f,-.1f,.1f),sampleCovariance(kField+V3(.2f,-.1f,.1f)));
    check((a->covariance_full()-b->covariance_full()).norm()<1e-8f,"zero-delay covariance overload parity");
    check(attitude(*a).angularDistance(attitude(*b))<1e-7f,"zero-delay quaternion parity");
    const auto p=b->covariance_full(); const Q q=attitude(*b);
    b->measurement_update_mag_only(kField,-rm);
    check(!b->lastMagDiag().accepted&&b->covariance_full()==p&&attitude(*b).coeffs()==q.coeffs(),
          "invalid covariance does not mutate filter or claim service");
}
