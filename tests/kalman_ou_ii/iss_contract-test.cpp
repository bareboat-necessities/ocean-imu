#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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
    //
    // A finite tau_b is what makes the residual-bias block UES; the
    // random-walk limit is explicitly outside the theorem.  The projection
    // radius is the other half: mean reversion bounds the bias in
    // distribution, not pathwise, and the bias enters the ISS bound only as an
    // input, so that input has to be bounded for the bound to say anything.
    if (!near(f.mekf_->get_acc_bias_time_constant(), 5000.0f))
        return fail("default residual accel-bias OU time constant changed");
    if (!(f.mekf_->get_acc_bias_time_constant() < std::numeric_limits<float>::max()))
        return fail("residual accel-bias block degenerated to a random walk");
    if (!near(f.mekf_->accel_bias_limit(), 0.5f))
        return fail("residual accel-bias projection radius changed");
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
    if (!near(f.P_factor_, 1.5f) ||
        !near(f.R_p0_x_factor_, 0.72f) || !near(f.R_p0_y_factor_, 0.72f))
        return fail("default OU-II anisotropy no longer matches proof audit");

    // ---- Pseudo-update cadence contract ---------------------------------
    // The deployed cadence is no longer a fixed 15 ms; it is
    // T_S = clip(c_T * tau_applied, T_min, T_max), holding T_S/tau fixed.  Pin
    // the constants so the bounded-gap hypothesis stays a checked statement
    // rather than an assumption about a number that has since moved.
    if (!near(PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT, 0.005f) ||
        !near(PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT, 0.25f) ||
        !near(PSEUDO_UPDATE_TAU_RATIO_DEFAULT, 0.015f / 1.1f)) {
        return fail("deployed OU-II pseudo-update cadence constants changed");
    }
    // The initial applied tau is the nominal 1.1 s, so the historical 15 ms
    // operating point is preserved exactly at startup.
    if (!near(f.mekf_->get_pseudo_update_period_s(), 0.015f))
        return fail("nominal p/v pseudo-update period is no longer 15 ms");

    // Over the clamped tau envelope the configured period stays inside the
    // audited interval.
    const float TS_at_tau_min =
        std::min(std::max(PSEUDO_UPDATE_TAU_RATIO_DEFAULT * MIN_TAU_S,
                          PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT),
                 PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT);
    const float TS_at_tau_max =
        std::min(std::max(PSEUDO_UPDATE_TAU_RATIO_DEFAULT * MAX_TAU_S,
                          PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT),
                 PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT);
    if (!near(TS_at_tau_min, 0.005f) || TS_at_tau_max > 0.165f ||
        TS_at_tau_max < 0.160f) {
        return fail("configured pseudo-update period left the audited envelope");
    }

    // apply_R_p0_tune_/apply_R_v0_tune_ scale the clamped base by
    // sqrt(T_S0 / T_S) and deliberately do not re-clamp, so the effective
    // enclosures are the base intervals widened by that factor.
    const float scale_lo = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / TS_at_tau_max);
    const float scale_hi = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / TS_at_tau_min);
    if (!(scale_lo > 0.30f && scale_lo < 0.31f) || !near(scale_hi, std::sqrt(3.0f), 1e-3f))
        return fail("cadence information-rate factor left its audited range");

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

    std::cout << "OU-II ISS proof contract passed\n";
    return 0;
}
