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

int main(){test_live_and_bias_release_inherit_state();test_attempt_is_not_acceptance();
    std::cout<<"OU3_SHIPPING_TRANSITION_PASS="<<(failures==0?"true":"false")<<'\n';return failures?1:0;}
