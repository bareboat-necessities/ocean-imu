// Pins implementation-side hypotheses of the OU-III local ISS theorem.
// Environmental excitation and reset-free intervals remain theorem assumptions;
// this test guards only source properties that can regress mechanically.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

static bool near(float a, float b, float rel = 2e-5f) {
    return std::fabs(a - b) <= rel * std::max(1.0f, std::max(std::fabs(a), std::fabs(b)));
}

static int fail(const char* msg) {
    std::cerr << "FAIL: " << msg << "\n";
    return 1;
}

static double A3_of(double t, double tau) {
    if (t <= 0.0) return 0.0;
    const int n = 20000;
    const double dt = t / double(n);
    double acc = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double s = double(i) * dt;
        const double f = (t - s) * (t - s) * std::exp(-s / tau);
        const double w = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        acc += w * f;
    }
    return 0.5 * acc * dt / 3.0;
}

static double detOS(const double t[4], double tau) {
    Eigen::Matrix4d O;
    for (int j = 0; j < 4; ++j) {
        O(j, 0) = 0.5 * t[j] * t[j];
        O(j, 1) = t[j];
        O(j, 2) = 1.0;
        O(j, 3) = A3_of(t[j], tau);
    }
    return O.determinant();
}

int main() {
    Filter f(true);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
    if (!f.mekf_) return fail("OU-III core was not created");

    if (f.mekf_->covariance_full().rows() != 21 ||
        f.mekf_->covariance_full().cols() != 21) {
        return fail("default OU-III state dimension is no longer 21");
    }

    if (!near(f.mekf_->get_acc_bias_time_constant(), 5000.0f))
        return fail("default residual accel-bias OU time constant changed");
    if (!(f.mekf_->get_acc_bias_time_constant() < std::numeric_limits<float>::max()))
        return fail("residual accel-bias block degenerated to a random walk");

    if (f.mekf_->legacy_aw_covariance_replacement())
        return fail("legacy raw a_w covariance replacement became default");
    if (!f.periodicAwCovarianceSync())
        return fail("default PSD a_w covariance synchronization was disabled");

    if (!near(MIN_TAU_S, 0.02f) || !near(MAX_TAU_S, 12.0f) ||
        !near(MAX_SIGMA_A, 6.0f) ||
        !near(MIN_R_S, 0.15f) || !near(MAX_R_S, 400.0f) ||
        !near(ACC_NOISE_FLOOR_SIGMA_DEFAULT, 0.12f)) {
        return fail("OU-III wrapper clamps changed");
    }
    if (!near(f.S_factor_, 1.0f))
        return fail("default OU-III acceleration anisotropy changed");
    if (!near(f.R_S_xy_factor_, 1.0f))
        return fail("default OU-III integral regularizer is no longer isotropic");

    {
        Filter probe;
        probe.setRSXYFactor(1.87f);
        if (!near(probe.R_S_xy_factor_, 1.87f))
            return fail("rho_xy > 1 is not expressible");
        probe.setRSXYFactor(9.0f);
        if (!near(probe.R_S_xy_factor_, 4.0f))
            return fail("rho_xy upper bound changed");
        probe.setRSXYFactor(-1.0f);
        if (!near(probe.R_S_xy_factor_, 0.0f))
            return fail("rho_xy lower bound changed");
    }

    if (!near(PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT, 0.005f) ||
        !near(PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT, 0.25f) ||
        !near(PSEUDO_UPDATE_TAU_RATIO_DEFAULT, 0.015f / 1.1f)) {
        return fail("OU-III pseudo-update cadence constants changed");
    }

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
        return fail("configured pseudo-update period left its bounded interval");
    }

    f.startup_stage_ = Filter::StartupStage::Live;
    f.mekf_->set_linear_block_enabled(true);
    if (!f.mekf_->linear_block_enabled())
        return fail("Live state did not enable the OU-III linear block");

    const float h = 0.005f;
    const float Tp = 7.0f;
    const float w = 2.0f * float(M_PI) / Tp;
    std::vector<double> fire_times;
    double t_now = 0.0;
    float prev_elapsed = f.mekf_->pseudo_update_elapsed_s_;

    for (int k = 0; k < 60000; ++k) {
        const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
        const Eigen::Vector3f acc(0.0f, 0.0f, -g_std + 1.2f * std::sin(w * float(t_now)));
        f.updateTime(h, gyro, acc);
        t_now += double(h);
        const float elapsed = f.mekf_->pseudo_update_elapsed_s_;
        if (elapsed < prev_elapsed && t_now > 200.0) fire_times.push_back(t_now);
        prev_elapsed = elapsed;
    }
    if (fire_times.size() < 64)
        return fail("pseudo-update scheduler did not run");

    const float tau_applied = f.getTauApplied();
    const float TS_cfg = f.getPseudoUpdatePeriodSec();
    if (!(tau_applied >= MIN_TAU_S && tau_applied <= MAX_TAU_S))
        return fail("applied tau left its clamp interval");
    if (!near(TS_cfg / tau_applied, PSEUDO_UPDATE_TAU_RATIO_DEFAULT, 1e-3f))
        return fail("realized T_S/tau is not the configured ratio");

    double gap_min = 1e30, gap_max = 0.0;
    for (size_t i = 1; i < fire_times.size(); ++i) {
        const double g = fire_times[i] - fire_times[i - 1];
        gap_min = std::min(gap_min, g);
        gap_max = std::max(gap_max, g);
    }
    const double h_d = double(h);
    if (gap_min < h_d - 1e-9)
        return fail("realized pseudo-update gap fell below one IMU step");
    if (gap_max > double(PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT) + h_d + 1e-9)
        return fail("realized pseudo-update gap exceeded T_S_max + h");
    if (gap_max - gap_min < 0.5 * h_d)
        return fail("realized cadence became near-exactly periodic; revisit bounded-gap proof");

    for (size_t i = 0; i + 3 < fire_times.size(); i += 7) {
        double t[4];
        for (int j = 0; j < 4; ++j) t[j] = fire_times[i + j] - fire_times[i];
        double lo = 1e30, hi = 0.0;
        for (int j = 1; j < 4; ++j) {
            const double g = t[j] - t[j - 1];
            lo = std::min(lo, g);
            hi = std::max(hi, g);
        }
        const double bound = std::pow(lo, 6.0) * std::exp(-3.0 * hi / double(tau_applied));
        const double det = std::fabs(detOS(t, double(tau_applied)));
        if (!(det > 0.0) || !(det >= bound))
            return fail("OU-III four-point observability determinant violated the lemma bound");
    }

    f.mekf_->set_linear_block_enabled(false);
    f.mekf_->set_acc_bias_updates_enabled(true);
    f.mekf_->set_acc_bias_time_constant(2.0f);
    const Eigen::Vector3f b0(0.20f, -0.10f, 0.05f);
    f.mekf_->set_initial_acc_bias(b0);
    f.mekf_->time_update(Eigen::Vector3f::Zero(), 0.10f);
    const Eigen::Vector3f expected = b0 * std::exp(-0.10f / 2.0f);
    if (!f.mekf_->get_acc_bias().isApprox(expected, 2e-6f))
        return fail("active residual accel-bias OU block is not contractive");

    f.accel_bias_locked_ = false;
    f.startup_stage_ = Filter::StartupStage::TunerWarm;
    f.enterLive_();
    if (!f.mekf_->acc_bias_updates_enabled_)
        return fail("post-unlock Live state did not enable accel-bias OU propagation");
    if (!f.mekf_->linear_block_enabled())
        return fail("default Live state did not re-enable the OU-III linear block");

    std::cout << "OU-III local ISS implementation contract PASS\n";
    return 0;
}
