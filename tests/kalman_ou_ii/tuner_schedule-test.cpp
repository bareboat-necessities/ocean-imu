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

static bool near(float a, float b, float eps = 1e-6f) {
    return std::fabs(a - b) <= eps * std::max(1.0f, std::max(std::fabs(a), std::fabs(b)));
}

int main() {
    Filter f(false);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
    const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
    const Eigen::Vector3f acc(0.0f, 0.0f, -g_std);
    f.initialize_from_acc(acc);
    f.startup_stage_ = Filter::StartupStage::Live;
    f.mekf_->set_linear_block_enabled(true);
    f.time_ = 1.0;
    f.last_adapt_time_sec_ = 0.0;
    f.adapt_every_secs_ = 0.05f;

    const float active_tau_before = f.mekf_->tau_aw;
    const float active_rp_before = f.mekf_->R_p0(2, 2);
    const float active_rv_before = f.mekf_->R_v0(2, 2);

    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 0.7f);
    const float staged_tau = f.tune_.tau_applied;
    const float staged_rp_var = f.tune_.R_p0_std_applied * f.tune_.R_p0_std_applied;
    const float staged_rv_var = f.tune_.R_v0_std_applied * f.tune_.R_v0_std_applied;

    if (!f.online_tune_apply_pending_) {
        std::cerr << "FAIL: adaptation was not staged\n";
        return 1;
    }
    if (!near(f.mekf_->tau_aw, active_tau_before) ||
        !near(f.mekf_->R_p0(2, 2), active_rp_before) ||
        !near(f.mekf_->R_v0(2, 2), active_rv_before)) {
        std::cerr << "FAIL: current-sample tuner output changed the active OU-II schedule\n";
        return 1;
    }
    if (near(staged_tau, active_tau_before)) {
        std::cerr << "FAIL: test did not create a distinct staged candidate\n";
        return 1;
    }

    // Keep the check isolated from a second tuner update without invoking the
    // public disable setter, which intentionally cancels a pending candidate.
    f.enable_tuner_ = false;
    f.updateTime(0.005f, gyro, acc);

    if (!near(f.mekf_->tau_aw, staged_tau) ||
        !near(f.mekf_->R_p0(2, 2), staged_rp_var) ||
        !near(f.mekf_->R_v0(2, 2), staged_rv_var)) {
        std::cerr << "FAIL: staged OU-II schedule was not committed at next sample\n";
        return 1;
    }
    if (f.online_tune_apply_pending_) {
        std::cerr << "FAIL: pending flag survived the commit\n";
        return 1;
    }
    std::cout << "OU-II predictable tuner scheduling passed\n";
    return 0;
}
