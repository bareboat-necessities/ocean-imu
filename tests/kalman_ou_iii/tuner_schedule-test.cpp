#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

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

    const float before = f.getTauApplied();
    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f);
    const float staged = f.tune_.tau_applied;
    if (!f.online_tune_apply_pending_) {
        std::cerr << "FAIL: adaptation was not staged\n";
        return 1;
    }
    if (std::fabs(f.getTauApplied() - before) > 1e-6f) {
        std::cerr << "FAIL: current-sample tuner output changed the active OU-III schedule\n";
        return 1;
    }
    if (std::fabs(staged - before) < 1e-5f) {
        std::cerr << "FAIL: test did not create a distinct staged candidate\n";
        return 1;
    }

    f.enable_tuner_ = false;
    f.updateTime(0.005f, gyro, acc);
    if (std::fabs(f.getTauApplied() - staged) > 1e-6f) {
        std::cerr << "FAIL: staged OU-III schedule was not committed at next sample\n";
        return 1;
    }
    if (f.online_tune_apply_pending_) {
        std::cerr << "FAIL: pending flag survived the commit\n";
        return 1;
    }
    std::cout << "OU-III predictable tuner scheduling passed\n";
    return 0;
}
