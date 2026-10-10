// Deterministic paired comparison on the literal 21-state OU-III core.
// Reference rewinds actual filter snapshots and replays every IMU/S/acc/mag
// operation. It is test-only and is not an alternative deployed estimator.
#define EIGEN_NON_ARDUINO
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
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
using Core = ocean_imu::kalman::Kalman3D_Wave_OU_III<float>;
const float g_std=9.80665f;
static const V3 kField(20.3f,0,46.8f), hard(8,-6,3);
static const M3 rm=.64f*M3::Identity(); // unchanged shipping .8 uT sigma
static void check(bool ok,const char* what) {
    if(!ok) {std::fprintf(stderr,"FAIL: %s\n",what); std::exit(1);}
}
static Core core() {
    Core f(V3::Constant(.12f),V3::Constant(.00135f),V3::Constant(.8f));
    f.set_mag_world_ref(kField);
    f.initialize_from_attitude(Q::Identity(),.03f,.08f);
    f.set_acc_bias_updates_enabled(false);
    return f;
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
static void advance(Core& f,const Row& r,int index) {
    f.time_update(r.gyro,r.dt);
    if(index%10==0) f.applyIntegralZeroPseudoMeas();
    f.measurement_update_acc_only(r.acc);
}
struct Score {
    double ss=0,peak=0,vector_ss=0,tail_ss=0,reference_ss=0;
    int n=0,m=0,tail_n=0,accepted=0,attempted=0;
    void add(Core& f,const Row& row,const Q& reference) {
        const double e=error(f.quaternion_boat(),row.truth);
        check(std::isfinite(e)&&f.covariance_full().allFinite(),"finite complete estimator state");
        if(row.time<3000000) return;
        ss+=e*e; peak=std::max(peak,std::abs(e)); ++n;
        const double diff=error(f.quaternion_boat(),reference); reference_ss+=diff*diff;
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
    auto legacy=std::make_unique<Core>(core()), aligned=std::make_unique<Core>(core()), ref=std::make_unique<Core>(core());
    // Test-only rewind snapshots. 40 slots comfortably exceed 40 ms delay.
    std::vector<std::optional<Core>> snapshots(40);
    RotationHistory<> h; h.push(0,V3::Zero());
    Score scores[3]; RotationConsistency consistency;
    int next=1; double diagnostic_ss=0; int diagnostic_n=0;
    for(int i=1;i<=N;++i) {
        Row& r=rows[size_t(i)];
        h.push(r.time,r.gyro-aligned->gyroscope_bias_body(),
               covarianceBound(aligned->gyroscope_bias_covariance_body()));
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
            Uncertainty u; u.gyro_density=.00135f;
            u.bias_variance=covarianceBound(aligned->gyroscope_bias_covariance_body());
            // Actual timestamps in this fixture are known; unknown hardware
            // conversion epochs are intentionally not imitated by invented data.
            V3 v; M3 covariance;
            check(align(m.mag,hard,rm,rotation,u,v,covariance),"transported measurement covariance valid");
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
            const auto d=consistency.observe(h,uint32_t(next),m.time,m.mag,hard,rm,u);
            if(d.status==Status::Ready) {diagnostic_ss+=d.nis;++diagnostic_n;}
            ++next;
        }
        scores[0].add(*legacy,r,ref->quaternion_boat());
        scores[1].add(*aligned,r,ref->quaternion_boat());
        scores[2].add(*ref,r,ref->quaternion_boat());
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
    if(interference) check(std::sqrt(scores[1].tail_ss/scores[1].tail_n)<2,"recovery after external kField removal");
}
static void compatibility() {
    auto a=std::make_unique<Core>(core()), b=std::make_unique<Core>(core());
    a->measurement_update_mag_only(kField+V3(.2f,-.1f,.1f));
    b->measurement_update_mag_only(kField+V3(.2f,-.1f,.1f),rm);
    check((a->covariance_full()-b->covariance_full()).norm()<1e-8f,"zero-delay covariance overload parity");
    check(a->quaternion_boat().angularDistance(b->quaternion_boat())<1e-7f,"zero-delay quaternion parity");
    const auto p=b->covariance_full(); const Q q=b->quaternion_boat();
    b->measurement_update_mag_only(kField,-rm);
    check(!b->lastMagDiag().accepted&&b->covariance_full()==p&&b->quaternion_boat().coeffs()==q.coeffs(),
          "invalid covariance does not mutate filter or claim service");
    // Wind heel must transform both the vector and its full covariance.
    b=std::make_unique<Core>(*a); a->update_wind_heel(.2f); b->update_wind_heel(.2f);
    a->measurement_update_mag_only(kField);
    b->measurement_update_mag_only(kField,rm);
    check((a->covariance_full()-b->covariance_full()).norm()<2e-6f,"isotropic covariance deheel parity");
    using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;
    auto legacy=std::make_unique<Fusion>(), timed=std::make_unique<Fusion>(); Fusion::Config cfg;
    legacy->begin(cfg); timed->begin(cfg); Observation observation;
    observation.covariance=cfg.sigma_m.array().square().matrix().asDiagonal();
    for(int i=1;i<=26000;++i) {
        const float t=float(i)*.005f;
        const V3 gyro(.03f*std::sin(t),.02f*std::cos(t),0);
        const V3 acc(0,0,-g_std);
        legacy->update(.005f,gyro,acc); timed->update(.005f,gyro,acc);
        if(i%8==0) {
            observation.proxy_bw=timed->raw().startupProxyQuat();
            observation.attitude_bw=timed->attitudeQuat();
            observation.accel=acc;observation.gyro=gyro;
            legacy->updateMag(kField);timed->updateMagTransported(kField,observation);
        }
    }
    check(legacy->hasRefinedMagReference()&&timed->hasRefinedMagReference(),"both facade paths refine magnetic reference");
    check((legacy->raw().mekf().covariance_full()-timed->raw().mekf().covariance_full()).norm()<2e-5f,
          "identity transport preserves complete startup/refinement/continuous path");
    check(legacy->attitudeQuat().angularDistance(timed->attitudeQuat())<1e-5f,"identity facade attitude parity");
}
int main() {
    compatibility();
    std::puts("scenario,delay_ms,mode,vector_rms_uT,heading_rms_deg,heading_peak_deg,rejection_rate,tail_heading_rms_deg,heading_vs_replay_rms_deg");
    for(int d:{5,10,20,40}) {paired("yaw",d,1);paired("combined_variable",d,2);}
    paired("async_jitter",20,2,true);
    paired("gyro_bias_error",40,2,true,true);
    paired("soft_iron_residual",20,2,false,false,true);
    paired("interference_recovery",20,1,false,false,false,true);
    paired("quiet",40,0,true);
    std::fprintf(stderr,"OU-III magnetic timing regression PASS\n");
}
