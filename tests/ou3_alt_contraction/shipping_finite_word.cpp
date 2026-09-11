// Host-only finite implementation correspondence.  No stability assertion.
// The probe headers are generated in a temporary include overlay; no shipping
// implementation, access control, configuration, or state is changed in src/.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include "util/W3dSimCommon.h"
#include "kalman_ou_common/KalmanOUCoreMath.h"

namespace ou3_alt_probe {
template<class Filter> void record(int, const Filter&);
template<class Filter, class V, class N, class S, class K, class Y>
void innovation(int, const Filter&, const V&, const N&, const S&, const K&, const Y&);
template<class Filter> struct Scope {
    using C = std::remove_reference_t<Filter>;
    int kind;
    const C& filter;
    Scope(int k, const C& f) : kind(k), filter(f) { record(kind, filter); }
    ~Scope() { record(kind+1, filter); }
};
}

// As in the existing host operation observer, visibility only is changed.
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

namespace {
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
using Vector3 = Eigen::Matrix<long double, 3, 1>;
using Quaternion = Eigen::Quaternion<long double>;
constexpr float dt = 1.0f/200.0f;
constexpr int word_steps = 600;
constexpr int before_edge_steps = 299;
constexpr long double pi = 3.141592653589793238462643383279502884L;

// One smooth, bounded oscillatory reference, used ONLY for implementation
// regression.  It is neither a source-uniform cover nor a seed qualification.
struct Physical {
    Quaternion q_wb;
    Vector3 p, v, a, omega, beta, bg;
};
Physical physical(long double t) {
    Physical p;
    const long double amp[] = {0.12L, 0.07L, 0.16L};
    const long double freq[] = {0.17L, 0.23L, 0.13L};
    const long double phase[] = {0.21L, -0.35L, 0.6L};
    for (int i=0; i<3; ++i) {
        const auto w = 2*pi*freq[i];
        p.p[i] = amp[i]*std::sin(w*t+phase[i]);
        p.v[i] = amp[i]*w*std::cos(w*t+phase[i]);
        p.a[i] = -amp[i]*w*w*std::sin(w*t+phase[i]);
    }
    const auto r=0.045L*std::sin(0.6L*t), rd=0.027L*std::cos(0.6L*t);
    const auto b=0.03L*std::cos(0.83L*t), bd=-0.0249L*std::sin(0.83L*t);
    const auto y=0.07L*std::sin(0.37L*t), yd=0.0259L*std::cos(0.37L*t);
    Quaternion q_bw = Eigen::AngleAxis<long double>(y,Vector3::UnitZ())
                    * Eigen::AngleAxis<long double>(b,Vector3::UnitY())
                    * Eigen::AngleAxis<long double>(r,Vector3::UnitX());
    p.q_wb=q_bw.conjugate();
    p.omega << rd-yd*std::sin(b),
               bd*std::cos(r)+yd*std::sin(r)*std::cos(b),
              -bd*std::sin(r)+yd*std::cos(r)*std::cos(b);
    // BIAS2's non-relaxing truth endpoint is deliberately present.  BIAS0/1
    // are checked separately by exact recurrence identities, not inferred.
    p.beta << 0.025L, -0.018L, 0.012L;
    p.bg << 0.0001L+0.00002L*std::sin(0.011L*t),
           -0.0002L+0.00001L*std::cos(0.009L*t), 0.0001L;
    return p;
}

using Event = std::vector<double>;
struct Sample { int tick; std::vector<Event> events; };
int tick_now=0;
int source_phase=1;
long double live_origin=-1;
bool observing=false;
Sample current{};

template<class Derived>
void append(Event& out, const Eigen::MatrixBase<Derived>& a) {
    for (Eigen::Index i=0;i<a.rows();++i)
        for (Eigen::Index j=0;j<a.cols();++j) out.push_back(static_cast<double>(a(i,j)));
}
template<class T> void write(std::ofstream& out, const T& x) {
    out.write(reinterpret_cast<const char*>(&x), sizeof(x));
    if (!out) throw std::runtime_error("cannot write regression output");
}
void write_sample(std::ofstream& out, const Sample& sample) {
    const std::uint32_t count=static_cast<std::uint32_t>(sample.events.size());
    write(out, sample.tick); write(out, count);
    for (const auto& e:sample.events) {
        const std::uint32_t size=static_cast<std::uint32_t>(e.size()); write(out,size);
        out.write(reinterpret_cast<const char*>(e.data()), size*sizeof(double));
    }
}

template<class Filter>
Event snapshot(int kind, const Filter& f) {
    const long double time=(tick_now-(source_phase ? 0 : 1))*static_cast<long double>(dt);
    const float ph=f.acc_bias_updates_enabled_
        ? std::exp(-f.last_dt_/std::max(1.0e-3f,f.tau_bacc_)):1.0f;
    Event out={double(tick_now),double(kind),double(source_phase),double(f.acc_bias_updates_enabled_),
        double(time),double(live_origin),double(f.last_dt_),double(f.tau_aw),
        double(f.pseudo_update_period_s_),double(f.pseudo_update_elapsed_s_),
        double(f.acc_bias_limit_),double(ph),double(f.have_prev_omega_),
        double(f.aw_covariance_floor_pending_),double(f.param_rw_enabled_),0.0};
    append(out,f.xext);
    out.insert(out.end(),{double(f.qref.w()),double(f.qref.x()),double(f.qref.y()),double(f.qref.z())});
    append(out,f.Pext); append(out,f.Racc); append(out,f.Rmag); append(out,f.R_S);
    append(out,f.Sigma_aw_stat); append(out,f.aw_covariance_floor_target_);
    append(out,f.v2ref); append(out,f.last_gyr_bias_corrected); append(out,f.Q_bacc_);
    out.push_back(double(f.tau_bacc_)); out.push_back(double(f.gravity_magnitude_));
    append(out,f.k_a_); out.push_back(35.0); out.push_back(35.0);
    if(out.size()!=549) throw std::logic_error("snapshot schema changed");
    return out;
}

// Recorded critical frontend coordinates.  The complete runtime remains the
// same persistent C++ object.  This list is an audit, NOT a claim that all xi
// coordinates/guards have acquired a source-uniform mathematical enclosure.
Event frontend(const Fusion& fusion) {
    const auto& f=fusion.raw(); const auto& w=f.wave_period_;
    const auto& b=f.sigma_wave_band_; const auto& s=f.tuner_;
    Event a={double(f.time_), double(f.startup_stage_t_), double(f.last_adapt_time_sec_),
        double(f.online_tune_apply_pending_), double(f.last_aw_cov_sync_sec_),
        double(f.accel_bias_locked_), double(f.acc_bias_hold_), double(f.mag_updates_applied_),
        double(f.first_mag_update_time_), double(f.tilt_over_limit_sec_),double(f.tilt_reset_cooldown_sec_),
        double(f.freq_hz_),double(f.freq_hz_slow_),double(f.f_raw),
        double(f.tune_.tau_applied),double(f.tune_.sigma_applied),double(f.tune_.RS_applied),
        double(f.tau_target_),double(f.sigma_target_),double(f.RS_target_),
        double(w.accel_prev_),double(w.high_pass_1_),double(w.high_pass_1_prev_),double(w.high_pass_2_),
        double(w.velocity_),double(w.elevation_),double(w.velocity_mean_),double(w.velocity_sq_),
        double(w.elevation_mean_),double(w.elevation_sq_),double(w.weight_),double(w.elapsed_sec_),
        double(w.raw_period_sec_),double(w.log_period_sec_),double(w.usable_period_),
        double(w.last_moment_horizon_sec_),double(w.last_log_horizon_sec_),
        double(b.lowpass_low_),double(b.band_),double(b.p00_),double(b.p01_),double(b.p11_),
        double(b.low_hz_),double(b.high_hz_),double(b.ready_),
        double(s.tau_var_sec),double(s.frequency_hz),double(s.A_mean.value),double(s.A_mean.weight),
        double(s.A_sq.value),double(s.A_sq.weight),double(f.vertical_accel_comp_.verticalAccelUpMs2()),
        double(f.accel_guard_.engagement()),double(f.accel_guard_.removedRms()),
        double(fusion.t_),double(fusion.mag_ref_set_),double(fusion.mag_refine_started_),
        double(fusion.mag_refine_done_),double(fusion.live_time_sec_),double(fusion.mag_refine_time_sec_)};
    const auto q=f.startupProxyQuat();
    a.insert(a.end(),{double(q.w()),double(q.x()),double(q.y()),double(q.z())});
    append(a,f.racc_effective_); append(a,f.Racc_nominal_);
    append(a,fusion.mag_world_ref_uT_);append(a,fusion.mag_hard_iron_body_uT_);
    append(a,fusion.gravity_gate_acc_world_lpf_.state);
    return a;
}
}

namespace ou3_alt_probe {
template<class Filter> void record(int kind, const Filter& f) {
    if (!observing) return;
    if (kind==1) source_phase=0;
    if (kind==2) source_phase=1;
    auto event=snapshot(kind,f);
    if (kind==2) {
        append(event,f.F_AA_scratch_); append(event,f.Q_AA_scratch_);
        append(event,f.F_LL_scratch_); append(event,f.Q_LL_scratch_);
    }
    current.events.push_back(std::move(event));
}
template<class Filter,class V,class N,class S,class K,class Y>
void innovation(int kind,const Filter& f,const V& r,const N& n,const S& s,const K& k,const Y& y) {
    if(!observing) return;
    auto event=snapshot(kind,f);
    append(event,r);append(event,n);append(event,s);append(event,k);append(event,y);
    current.events.push_back(std::move(event));
}
}

int main(int argc,char**argv) {
    try {
        if(argc!=2) throw std::runtime_error("usage: shipping_finite_word OUTPUT_DIRECTORY");
        const std::filesystem::path output(argv[1]);
        std::filesystem::create_directories(output);
        std::ofstream samples(output/"samples.bin",std::ios::binary);
        std::ofstream hs(output/"H18.bin",std::ios::binary),
                      edge(output/"H18-A21.bin",std::ios::binary),
                      active(output/"A21.bin",std::ios::binary);
        Fusion fusion; Fusion::Config cfg; fusion.begin(cfg);
        std::deque<Sample> rolling;
        int H_count=0,A_count=0,edge_count=0,live_tick=-1,edge_tick=-1;
        std::size_t sample_fields=0, prefix_count=0;
        bool previous_active=false;
        for(tick_now=1;tick_now<=60000;++tick_now) {
            const long double t=tick_now*static_cast<long double>(dt);
            auto p=physical(t);
            const Eigen::Vector3f gyro=(p.omega+p.bg).cast<float>();
            const Vector3 gravity(0,0,9.80665L);
            const Eigen::Vector3f acc=(p.q_wb*(p.a-gravity)+p.beta).cast<float>();
            const Eigen::Vector3f mag=(p.q_wb*Vector3(22,1,43)).cast<float>();
            const bool was_live=fusion.isLive();
            current={tick_now,{}}; source_phase=0;
            observing=bool(ALT_OBSERVE)&&was_live;
            fusion.update(dt,gyro,acc,35.0f);
            source_phase=1;
            if(tick_now%8==0) fusion.updateMag(mag);
            const bool is_active=fusion.raw().mekf().acc_bias_updates_enabled_;
            if(!was_live&&fusion.isLive()) {
                live_origin=t; live_tick=tick_now; previous_active=is_active;
            }
            Event boundary=snapshot(90,fusion.raw().mekf());
            auto front=frontend(fusion);boundary.insert(boundary.end(),front.begin(),front.end());
            if(!sample_fields) sample_fields=boundary.size();
            if(sample_fields!=boundary.size()) throw std::logic_error("sample serialization changed");
            samples.write(reinterpret_cast<const char*>(boundary.data()),boundary.size()*sizeof(double));
            if(was_live) {
                prefix_count+=current.events.size();
                if(H_count<word_steps) {
                    if(is_active) throw std::runtime_error("regression did not provide a full H18 word");
                    write_sample(hs,current);++H_count;
                }
                if(edge_tick<0 && is_active&&!previous_active) {
                    if(rolling.size()!=before_edge_steps) throw std::runtime_error("edge history missing");
                    for(const auto& old:rolling) {write_sample(edge,old);++edge_count;}
                    write_sample(edge,current);++edge_count;edge_tick=tick_now;
                } else if(edge_tick>=0 && edge_count<word_steps) {
                    write_sample(edge,current);++edge_count;
                }
                if(edge_tick>=0 && tick_now>edge_tick && A_count<word_steps) {
                    if(!is_active) throw std::runtime_error("A21 regression deactivated unexpectedly");
                    write_sample(active,current);++A_count;
                }
                rolling.push_back(std::move(current));
                if(rolling.size()>before_edge_steps) rolling.pop_front();
            }
            previous_active=is_active;
            if(A_count==word_steps && edge_count==word_steps && H_count==word_steps) break;
        }
        if(H_count!=word_steps||A_count!=word_steps||edge_count!=word_steps)
            throw std::runtime_error("failed to obtain all three literal 600-step coverage words");
        std::cout<<std::setprecision(17)
            <<"{\"qualification\":\"IMPLEMENTATION_REGRESSION_ONLY\",\"ticks\":"<<tick_now
            <<",\"sample_fields\":"<<sample_fields<<",\"live_tick\":"<<live_tick
            <<",\"edge_tick\":"<<edge_tick<<",\"dt\":"<<double(dt)
            <<",\"H18_steps\":"<<H_count<<",\"A21_steps\":"<<A_count
            <<",\"edge_steps\":"<<edge_count<<",\"observed_prefixes\":"<<prefix_count
            <<",\"ALT_LIVE_PASS\":false}\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';return 1;
    }
}
