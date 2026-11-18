// Shipping OU-III implementation invariants used by the stability analysis.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#define private public
#define protected public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
#undef protected
using namespace ocean_imu::ou3;

const float g_std=9.80665f;
using Filter=SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
static bool near(float a,float b,float rel=2e-5f){return std::fabs(a-b)<=rel*std::max(1.0f,std::max(std::fabs(a),std::fabs(b)));}
static int fail(const char* msg){std::cerr<<"FAIL: "<<msg<<"\n";return 1;}

// Use a circular orbit to check the frequency moment in an arbitrary
// horizontal basis. No absolute vessel heading or wave phase is needed.
static bool waveBandContract() {
    HorizontalWaveBand band;
    constexpr float period = 5.0f, dt = 0.01f;
    constexpr float omega = 2.0f * std::numbers::pi_v<float> / period;
    for (int i = 0; i < 16000; ++i) {
        const float phase = omega * dt * static_cast<float>(i);
        const Eigen::Vector2f p(std::cos(phase), std::sin(phase));
        const Eigen::Vector2f v(-omega * p.y(), omega * p.x());
        band.update(dt, period, -v / (omega * omega), p, v);
        if (i < 1490 && band.weight() != 0.0f) return false;
    }
    if (!near(band.beta(), 1.0f / (omega * omega), 1e-4f) ||
        band.residualVariance() > 1e-5f || band.weight() != 1.0f) return false;
    band.reset();
    return band.beta() == 0.0f && band.weight() == 0.0f;
}

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
    if (!waveBandContract()) return fail("horizontal wave period or warmup");
    // Compare the shipped horizontal solve against its dense observation.
    // Cross-check the raw posterior and vertical outputs on that same sample.
    constexpr int velocity = 6, position = 9, integral = 12;
    auto& core = *f.mekf_;
    core.xext.segment<3>(velocity) << 0.7f, -0.4f, 0.5f;
    core.xext.segment<3>(integral) << 0.2f, 0.1f, 0.4f;
    const auto state = core.xext;
    const auto covariance = core.covariance_full();
    Eigen::Matrix<float, 2, 21> H = Eigen::Matrix<float, 2, 21>::Zero();
    H.block<2,2>(0, integral).setIdentity();
    H.block<2,2>(0, velocity).setIdentity();
    const Eigen::Matrix2f innovation = H * covariance * H.transpose() + 0.1f * Eigen::Matrix2f::Identity();
    const Eigen::Matrix<float, 21, 1> projected = state - covariance * H.transpose() * innovation.inverse() * (H * state);
    const auto wave = core.wave_motion(1.0f, 0.1f, 1.0f);
    const auto raw = core.wave_motion(1.0f, 0.1f, 0.0f);
    if (!wave.position.head<2>().isApprox(projected.segment<2>(position), 1e-5f) ||
        !wave.velocity.head<2>().isApprox(projected.segment<2>(velocity), 1e-5f) ||
        wave.position.z() != raw.position.z() || wave.velocity.z() != raw.velocity.z() ||
        raw.position != core.get_position() || raw.velocity != core.get_velocity() ||
        core.xext != state || core.covariance_full() != covariance)
        return fail("horizontal wave projection or raw posterior changed");
    std::cout<<"OU-III shipping implementation contract PASS\n"; return 0;
}
