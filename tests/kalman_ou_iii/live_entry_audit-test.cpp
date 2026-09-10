// Executable audit of the shipping wrapper; no replacement estimator is used.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Inner = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
using Outer = SeaStateFusion_OU_III<TrackerType::KALMANF>;
using V3 = Eigen::Vector3f;
static int failures = 0;
static void require(bool ok, const char* msg) {
    if (!ok) { ++failures; std::cerr << "FAIL: " << msg << '\n'; }
}
template<class A, class B> static bool same(const A& a, const B& b) {
    return a.rows()==b.rows() && a.cols()==b.cols() && (a.array()==b.array()).all();
}
static void zero_motion(const Inner& f) {
    const auto& m=f.mekf();
    require(m.get_velocity().isZero(0), "held v is exactly zero");
    require(m.get_position().isZero(0), "held p is exactly zero");
    require(m.get_integral_displacement().isZero(0), "held S is exactly zero");
    require(m.get_world_accel().isZero(0), "held aw is exactly zero");
    require(m.gyroscope_bias().isZero(0), "proxy does not transfer gyro bias");
    require(m.get_acc_bias().isZero(0), "held ba is exactly zero");
}
static void test_go_live_is_not_a_linear_reset() {
    Inner f(true);
    f.initialize(V3::Constant(.2f),V3::Constant(.01f),V3::Constant(.3f));
    auto& m=f.mekf();
    // Deliberately reach a nonzero S and nonzero cross-covariances first.
    m.initialize_from_truth(V3(1,2,3),V3(.4f,.2f,-.3f),Eigen::Quaternionf::Identity(),V3(.1f,-.2f,.3f));
    m.set_initial_acc_bias(V3(.01f,.02f,.03f));
    m.set_pseudo_update_period_s(1.0f);
    m.time_update(V3(.01f,0,0),.005f);
    const auto x=m.xext;
    const auto p=m.covariance_full();
    const float elapsed=m.pseudo_update_elapsed_s_;
    require(m.get_integral_displacement().norm()>0, "sentinel S is nonzero");
    require(p.block<3,3>(6,9).norm()>0, "sentinel v/p cross-covariance exists");
    const auto q=Eigen::Quaternionf(Eigen::AngleAxisf(.12f,V3::UnitX()));
    f.goLive(q,.035f,.087f,false);
    require(same(m.xext.tail<18>(),x.tail<18>()), "goLive preserves ALL additive states including biases");
    const auto after=m.covariance_full();
    require(same(after.block<9,9>(6,6),p.block<9,9>(6,6)), "goLive preserves v/p/S covariance and cross-covariances");
    require(after.block<6,12>(0,6).isZero(0), "attitude initializer clears entire base/linear cross-block");
    require(after.block<3,15>(15,0).isZero(0), "aw covariance reset clears aw/earlier cross-block");
    require(after.block<3,3>(15,18).isZero(0), "aw covariance reset clears aw/ba cross-block");
    require(m.pseudo_update_elapsed_s_==elapsed, "goLive preserves sub-deadline scheduler progress");
    // H18 -> A21 is an enable/floor, not a state initialization.
    const auto xh=m.xext; const auto ph=m.covariance_full();
    m.set_acc_bias_updates_enabled(true);
    require(same(m.xext,xh), "H18 -> A21 preserves nominal states");
    require(same(m.covariance_full().block<18,18>(0,0),ph.block<18,18>(0,0)), "A21 release preserves motion covariance");
    for(int j=18;j<21;++j)
        require(m.covariance_full()(j,j)>=.004f*.004f,"A21 release floors ba covariance");
    const auto before_reset=m.xext;
    m.initialize_from_acc_preserve_yaw(V3(0,0,-g_std));
    require(same(m.xext.tail<18>(),before_reset.tail<18>()), "tilt re-lock preserves additive states");
}
static void test_fresh_wrapper(bool mag, bool zero_proxy) {
    float qscalar=1.0f; float noise_max=0.0f;
    Outer f; Outer::Config c; c.with_mag=mag;
    f.begin(c);
    const auto initial=f.raw().mekf().covariance_full();
    zero_motion(f.raw());
    int k=0;
    for(;k<40000 && !f.isLive();++k) {
        float acc_z=-g_std;
        if(zero_proxy) {
            // At level with zero gyro, q has only a scalar component. This
            // constructs sensor samples, not estimator state or a replacement
            // normalization. The shipping implementation is still executed.
            qscalar *= Mahony_AHRS<float>::invSqrt(qscalar*qscalar);
            const float row=qscalar*qscalar;
            acc_z=-g_std/row;
            if(row*acc_z!=-g_std) {
                const float lo=std::nextafter(acc_z,-INFINITY);
                const float hi=std::nextafter(acc_z,INFINITY);
                if(row*lo==-g_std) acc_z=lo;
                else if(row*hi==-g_std) acc_z=hi;
            }
            require(row*acc_z==-g_std,"constructed bounded sensor sample has zero proxy");
            noise_max=std::max(noise_max,std::abs(acc_z+g_std));
        }
        f.update(.005f,V3::Zero(),V3(0,0,acc_z));
        if(zero_proxy) require(f.raw().vertical_accel_comp_.verticalAccelUpMs2()==0,"actual proxy is exactly zero");
        // Keep magnetic callbacks separate, as in the actual device interface.
        if(mag && k%2==0 && k>=1400) f.updateMag(V3(20,0,40));
        zero_motion(f.raw());
        require(f.raw().mekf().pseudo_update_elapsed_s_==0,"no pre-Live pseudo scheduler accumulation");
    }
    require(f.isLive(), "quiet aligned source reaches runtime Live");
    if(zero_proxy) {
        require(!f.raw().wavePeriodUsable(), "quiet timeout has no usable wave period");
        require(f.liveTimeSec()>=150.0f && f.liveTimeSec()<150.01f,"quiet source uses timeout clock, not readiness");
        require(noise_max<.04f,"witness sensor disturbance is bounded by .04 m/s^2");
    }
    require(same(f.raw().mekf().covariance_full().block<9,9>(6,6),initial.block<9,9>(6,6)),"fresh Live v/p/S covariance is construction covariance");
    require(f.raw().mekf().get_integral_displacement().isZero(0),"fresh Live nominal S zero");
    const double physical_time=.005*static_cast<double>(k);
    // A C2 periodic BRMM-Q position can remain at 2.1 m for the first 200 s.
    // With S_true(t)=integral_0^t p_true, this gives S_true=2.1*t, NOT zero.
    const double source_S=std::sqrt(3.0)*2.1*physical_time;
    require(source_S>300,"reachable session-origin S at quiet handoff exceeds the old 300 assumption");
    require(physical_time<200,"entire startup stays on the source's constant-position plateau");
    std::cout << "LIVE_ENTRY " << (mag?"mag":"no_mag") << " zero_proxy="<<zero_proxy<<" noise_max="<<noise_max<<" samples="<<k
              <<" shipping_clock="<<f.liveTimeSec()<<" physical_time="<<physical_time
              <<" period_usable="<<f.raw().wavePeriodUsable()<<" nominal_S=0 session_S="<<source_S
              <<" handoff_origin_S=0 Pvv="<<initial(6,6)<<" Ppp="<<initial(9,9)<<" PSS="<<initial(12,12)<<'\n';
    const auto before=f.raw().mekf().covariance_full();
    f.update(.005f,V3::Zero(),V3(0,0,-g_std));
    require(!same(f.raw().mekf().covariance_full(),before),"first post-Live IMU sample propagates/updates covariance");
}
int main() {
    std::cout<<std::setprecision(10);
    test_go_live_is_not_a_linear_reset();
    test_fresh_wrapper(false,false); test_fresh_wrapper(true,false);
    test_fresh_wrapper(false,true); test_fresh_wrapper(true,true);
    std::cout<<"LIVE_ENTRY_AUDIT_PASS="<<(failures==0?"true":"false")<<'\n';
    return failures?1:0;
}
