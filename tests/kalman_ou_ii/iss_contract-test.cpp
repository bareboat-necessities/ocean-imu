#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#undef private

const float g_std = 9.80665f;

using Filter = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;
using Core = Kalman3D_Wave_OU_II<float>;

static bool near(float a, float b, float rel = 2e-5f) {
    return std::fabs(a - b) <= rel * std::max(1.0f, std::max(std::fabs(a), std::fabs(b)));
}

static int fail(const char* msg) {
    std::cerr << "FAIL: " << msg << "\n";
    return 1;
}

int main() {
    Filter f(true);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
    if (!f.mekf_) return fail("OU-II core was not created");

    // The default template configuration used by the fusion wrapper is the
    // 18-state proof configuration: attitude(3), gyro bias(3), v/p/a_w(9),
    // and residual accelerometer bias(3).
    if (f.mekf_->covariance_full().rows() != 18 ||
        f.mekf_->covariance_full().cols() != 18) {
        return fail("default OU-II state dimension is no longer 18");
    }

    // Pin proof-relevant default model and scheduler constants.
    if (!near(f.mekf_->get_acc_bias_time_constant(), 5000.0f))
        return fail("default residual accel-bias OU time constant changed");
    if (!near(f.mekf_->get_pseudo_update_period_s(), 0.015f))
        return fail("default p/v pseudo-update period changed");
    if (f.mekf_->legacy_aw_covariance_replacement())
        return fail("legacy raw a_w covariance replacement became default");
    if (!f.periodicAwCovarianceSync())
        return fail("default PSD a_w covariance synchronization was disabled");

    if (!near(MIN_TAU_S, 0.02f) || !near(MAX_TAU_S, 12.0f) ||
        !near(MAX_SIGMA_A, 6.0f) ||
        !near(MIN_R_p0_std, 0.05f) || !near(MAX_R_p0_std, 150.0f) ||
        !near(MIN_R_v0_std, 0.01f) || !near(MAX_R_v0_std, 40.0f) ||
        !near(ACC_NOISE_FLOOR_SIGMA_DEFAULT, 0.12f)) {
        return fail("proof-relevant OU-II wrapper clamps changed");
    }
    if (!near(f.P_factor_, 1.5f) || !near(f.R_p0_xy_factor_, 1.0f))
        return fail("default OU-II anisotropy no longer matches proof audit");

    // The exact source transition used by OU-II must retain a nonzero
    // acceleration-to-velocity coefficient over one pseudo-update gap.  This
    // is the determinant witness used by the translational UCO lemma.
    Eigen::Matrix3f Phi;
    Core::PhiAxis3x1_analytic(MIN_TAU_S, 0.015f, Phi);
    if (!(Phi(0, 2) > 0.0f) || !Phi.allFinite())
        return fail("OU-II integrated-OU transition lost translational observability witness");
    Core::PhiAxis3x1_analytic(MAX_TAU_S, 0.015f, Phi);
    if (!(Phi(0, 2) > 0.0f) || !Phi.allFinite())
        return fail("OU-II transition witness failed at maximum default tau");

    // Verify the residual accelerometer-bias block really contracts when its
    // OU propagation is active.  A shorter test tau makes the contraction
    // numerically visible while exercising the same source path.
    f.mekf_->set_linear_block_enabled(false);
    f.mekf_->set_acc_bias_updates_enabled(true);
    f.mekf_->set_acc_bias_time_constant(2.0f);
    const Eigen::Vector3f b0(0.20f, -0.10f, 0.05f);
    f.mekf_->set_initial_acc_bias(b0);
    f.mekf_->time_update(Eigen::Vector3f::Zero(), 0.10f);
    const Eigen::Vector3f expected = b0 * std::exp(-0.10f / 2.0f);
    if (!f.mekf_->get_acc_bias().isApprox(expected, 2e-6f))
        return fail("active residual accel-bias OU block is not contractive");

    // In the certified post-unlock Live regime, enterLive_ must enable the
    // residual-bias OU state rather than leave the startup freeze in force.
    f.accel_bias_locked_ = false;
    f.startup_stage_ = Filter::StartupStage::TunerWarm;
    f.enterLive_();
    if (!f.mekf_->acc_bias_updates_enabled_)
        return fail("post-unlock Live state did not enable accel-bias OU propagation");
    if (!f.mekf_->linear_block_enabled())
        return fail("default Live state did not enable the OU-II linear block");

    std::cout << "OU-II ISS proof contract passed\n";
    return 0;
}
