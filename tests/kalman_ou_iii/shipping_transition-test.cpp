// Literal shipping transition checks needed by the one OU-III stability path.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std=9.80665f;
using Filter=SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
using V3=Eigen::Vector3f;
static int failures=0;
static void require(bool ok,const char* msg){if(!ok){++failures;std::cerr<<"FAIL: "<<msg<<'\n';}}
template<class A,class B> static bool same(const A&a,const B&b){return a.rows()==b.rows()&&a.cols()==b.cols()&&(a.array()==b.array()).all();}

static void test_live_and_bias_release_inherit_state(){
    Filter f(true); f.initialize(V3::Constant(.2f),V3::Constant(.01f),V3::Constant(.3f)); auto&m=f.mekf();
    m.initialize_from_truth(V3(1,2,3),V3(.4f,.2f,-.3f),Eigen::Quaternionf::Identity(),V3(.1f,-.2f,.3f));
    m.set_initial_acc_bias(V3(.01f,.02f,.03f)); m.set_pseudo_update_period_s(1.0f); m.time_update(V3(.01f,0,0),.005f);
    const auto before=m.xext; const auto p_before=m.covariance_full(); const float elapsed=m.pseudo_update_elapsed_s_;
    f.goLive(Eigen::Quaternionf(Eigen::AngleAxisf(.12f,V3::UnitX())),.035f,.087f,false);
    require(same(m.xext.tail<18>(),before.tail<18>()),"Live handoff reset additive/bias states");
    require(same(m.covariance_full().block<9,9>(6,6),p_before.block<9,9>(6,6)),"Live handoff reset v/p/S covariance");
    require(m.pseudo_update_elapsed_s_==elapsed,"Live handoff reset scheduler progress");
    const auto x_hold=m.xext; const auto p_hold=m.covariance_full(); m.set_acc_bias_updates_enabled(true);
    require(same(m.xext,x_hold),"H18-to-A21 release reset nominal state");
    require(same(m.covariance_full().block<18,18>(0,0),p_hold.block<18,18>(0,0)),"H18-to-A21 release reset inherited covariance");
}

static void test_attempt_is_not_acceptance(){
    Filter f(true); f.initialize(V3::Constant(.2f),V3::Constant(.01f),V3::Constant(.3f));
    f.setMagDelaySec(0.0f); f.mekf().set_mag_world_ref(V3(20.0f,0.0f,40.0f));
    f.goLive(Eigen::Quaternionf::Identity(),.035f,.087f,false);
    const int before=f.mag_updates_applied_; const float nan=std::numeric_limits<float>::quiet_NaN();
    f.updateMag(V3(nan,0.0f,0.0f));
    require(f.mag_updates_applied_==before+1,"forwarded magnetic attempt counter did not advance");
    require(!f.mekf().lastMagDiag().accepted,"invalid magnetic attempt was reported as an applied correction");
    f.updateMag(V3(20.0f,0.0f,40.0f));
    require(f.mekf().lastMagDiag().accepted,"finite valid magnetic correction was not reported as applied");
}

// Primitive parity sentinels, not a claim of startup reachability.
static void test_bias_prediction_correction_and_projection() {
    Filter f(true);
    f.initialize(V3::Constant(.2f), V3::Constant(.01f), V3::Constant(.3f));
    auto& m = f.mekf();
    m.set_pseudo_update_period_s(10.0f);
    const float dt = .005f;
    const V3 truth(.1f, -.02f, .03f);
    const V3 increment(1e-6f, -1e-6f, 0.0f);
    for (bool active : {false, true}) {
        m.set_initial_acc_bias(V3(.03f, .01f, -.02f));
        m.set_acc_bias_updates_enabled(active);
        const V3 estimate_before = m.get_acc_bias();
        const V3 error_before = truth - estimate_before;
        const float phi = active ? std::exp(-dt / m.get_acc_bias_time_constant()) : 1.0f;
        m.time_update(V3::Zero(), dt);
        const V3 expected = phi * error_before + (1.0f - phi) * truth + increment;
        const V3 actual = truth + increment - m.get_acc_bias();
        require((actual - expected).norm() < 2e-7f, "joint truth/error prediction identity");
        if (!active) require(same(m.get_acc_bias(), estimate_before), "held estimate moved");
    }
    m.set_acc_bias_updates_enabled(true);
    m.Pext.setIdentity();
    m.Pext *= .001f;
    m.Pext(2, 18) = m.Pext(18, 2) = .0001f;
    m.set_mag_world_ref(V3(20.0f, 0.0f, 40.0f));
    const V3 before_correction = m.get_acc_bias();
    m.measurement_update_mag_only(V3(20.0f, .5f, 40.0f));
    require(m.lastMagDiag().accepted, "magnetic correction must actually apply");
    const V3 applied = m.K_scratch_.template block<3,3>(18,0) * m.lastMagDiag().r;
    V3 expected = before_correction + applied;
    const float norm = expected.norm();
    if (norm > m.accel_bias_limit()) expected *= m.accel_bias_limit() / norm;
    require((m.get_acc_bias() - expected).norm() < 2e-7f, "actual gain correction then projection");

    for (float radius : {0.0f, -.1f, .4f}) {
        m.set_accel_bias_limit(radius);
        m.xext.template segment<3>(18) = V3(.5f, 0.0f, 0.0f);
        const auto covariance = m.covariance_full();
        m.project_acc_bias_();
        require(std::abs(m.get_acc_bias().x() - (radius > 0.0f ? .4f : .5f)) < 1e-7f,
                "nonpositive radius disables projection");
        require(same(m.covariance_full(), covariance), "projection must not rewrite covariance");
    }
    m.xext.template segment<3>(18) = V3(std::numeric_limits<float>::quiet_NaN(),0,0);
    m.project_acc_bias_();
    require(m.get_acc_bias().isZero(0), "enabled projection resets nonfinite estimate");
    // An invalid attitude injection returns before projection; the finite-domain
    // proof must not silently assert that every attempted injection projects.
    m.xext.template segment<3>(18) = V3(.5f,0,0);
    m.xext(0) = std::numeric_limits<float>::quiet_NaN();
    m.applyQuaternionCorrectionFromErrorState();
    require(m.get_acc_bias().x() == .5f, "invalid attitude injection bypasses projection");
}

int main(){test_live_and_bias_release_inherit_state();test_attempt_is_not_acceptance();
    test_bias_prediction_correction_and_projection();
    std::cout<<"OU3_SHIPPING_TRANSITION_PASS="<<(failures==0?"true":"false")<<'\n';return failures?1:0;}
