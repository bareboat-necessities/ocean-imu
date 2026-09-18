// Shipping OU-III implementation invariants used by the stability analysis.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std=9.80665f;
using Filter=SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
static bool near(float a,float b,float rel=2e-5f){return std::fabs(a-b)<=rel*std::max(1.0f,std::max(std::fabs(a),std::fabs(b)));}
static int fail(const char* msg){std::cerr<<"FAIL: "<<msg<<"\n";return 1;}

int main(){
    Filter f(true);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),Eigen::Vector3f::Constant(0.00157f),Eigen::Vector3f::Constant(0.25f));
    if(!f.mekf_) return fail("OU-III core was not created");
    if(f.mekf_->covariance_full().rows()!=21||f.mekf_->covariance_full().cols()!=21) return fail("default shipping state is not 21-dimensional");
    if(!near(f.mekf_->get_acc_bias_time_constant(),5000.0f)) return fail("shipping accelerometer-bias estimator prior changed");
    if(!near(f.mekf_->accel_bias_limit(),0.4f)) return fail("shipping accelerometer-bias estimate projection radius changed");
    f.mekf_->set_acc_bias_time_constant(2.0f);
    const Eigen::Vector3f b0(0.20f,-0.10f,0.05f);
    f.mekf_->set_initial_acc_bias(b0);
    f.mekf_->set_acc_bias_updates_enabled(false);
    f.mekf_->time_update(Eigen::Vector3f::Zero(),0.10f);
    if(!f.mekf_->get_acc_bias().isApprox(b0,2e-7f)) return fail("held accelerometer-bias estimate prediction is not identity");
    f.mekf_->set_acc_bias_updates_enabled(true);
    f.mekf_->set_initial_acc_bias(b0);
    f.mekf_->time_update(Eigen::Vector3f::Zero(),0.10f);
    const Eigen::Vector3f expected=b0*std::exp(-0.10f/2.0f);
    if(!f.mekf_->get_acc_bias().isApprox(expected,2e-6f)) return fail("active accelerometer-bias estimate no longer uses configured OU predictor");
    if(!f.periodicAwCovarianceSync()) return fail("default covariance synchronization changed");
    if(!near(PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT,0.005f)||!near(PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT,0.15f)) return fail("shipping pseudo-update cadence bounds changed");
    std::cout<<"OU-III shipping implementation contract PASS\n"; return 0;
}
