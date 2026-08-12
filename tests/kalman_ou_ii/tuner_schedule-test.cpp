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

    // The deployed wrapper starts at tau=1.1 s and historically used a 15 ms
    // pseudo-update period.  The self-similar cadence must preserve that
    // operating point exactly while scaling subsequent periods with tau.
    if (!near(f.getPseudoUpdateTauRatio(), 0.015f / 1.1f) ||
        !near(f.getPseudoUpdatePeriodSec(), 0.015f)) {
        std::cerr << "FAIL: OU-II tau-scaled pseudo cadence did not preserve the nominal point\n";
        return 1;
    }

    const float active_tau_before = f.mekf_->tau_aw;
    const float active_rp_before = f.mekf_->R_p0(2, 2);
    const float active_rv_before = f.mekf_->R_v0(2, 2);
    const float active_pseudo_before = f.getPseudoUpdatePeriodSec();

    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 0.7f);
    const float staged_tau = f.tune_.tau_applied;
    const float staged_rp_base =
        std::min(std::max(f.tune_.R_p0_std_applied, f.MIN_R_p0_std_), f.MAX_R_p0_std_);
    const float staged_rv_base =
        std::min(std::max(f.tune_.R_v0_std_applied, f.MIN_R_v0_std_), f.MAX_R_v0_std_);
    const float staged_pseudo = std::min(
        std::max(f.pseudo_update_tau_ratio_ * staged_tau,
                 f.pseudo_update_period_min_s_),
        f.pseudo_update_period_max_s_);
    const float staged_cadence_scale = std::sqrt(
        PSEUDO_UPDATE_PERIOD_NOMINAL_S / staged_pseudo);
    const float staged_rp = staged_rp_base * staged_cadence_scale;
    const float staged_rv = staged_rv_base * staged_cadence_scale;
    const float staged_rp_var = staged_rp * staged_rp;
    const float staged_rv_var = staged_rv * staged_rv;

    if (!f.online_tune_apply_pending_) {
        std::cerr << "FAIL: adaptation was not staged\n";
        return 1;
    }
    if (!near(f.mekf_->tau_aw, active_tau_before) ||
        !near(f.mekf_->R_p0(2, 2), active_rp_before) ||
        !near(f.mekf_->R_v0(2, 2), active_rv_before) ||
        !near(f.getPseudoUpdatePeriodSec(), active_pseudo_before)) {
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
        !near(f.mekf_->R_v0(2, 2), staged_rv_var) ||
        !near(f.getPseudoUpdatePeriodSec(), staged_pseudo)) {
        std::cerr << "FAIL: staged OU-II schedule/cadence was not committed at next sample\n";
        return 1;
    }
    if (!near(f.getPseudoUpdatePeriodSec() / f.getTauApplied(),
              f.getPseudoUpdateTauRatio(), 1e-5f)) {
        std::cerr << "FAIL: T_S/tau is not invariant after the staged update\n";
        return 1;
    }

    // Both drift-correction channels fire on the same periodic tick, so both
    // have to hold r^2 * T_S -- the continuous-equivalent information rate --
    // at the value the historical 15 ms cadence produced.
    const float expected_rp_product =
        staged_rp_base * staged_rp_base * PSEUDO_UPDATE_PERIOD_NOMINAL_S;
    const float expected_rv_product =
        staged_rv_base * staged_rv_base * PSEUDO_UPDATE_PERIOD_NOMINAL_S;
    if (!near(f.mekf_->R_p0(2, 2) * f.getPseudoUpdatePeriodSec(),
              expected_rp_product, 1e-5f) ||
        !near(f.mekf_->R_v0(2, 2) * f.getPseudoUpdatePeriodSec(),
              expected_rv_product, 1e-5f)) {
        std::cerr << "FAIL: cadence compensation did not preserve r^2*T_S information rate\n";
        return 1;
    }
    if (f.online_tune_apply_pending_) {
        std::cerr << "FAIL: pending flag survived the commit\n";
        return 1;
    }

    // Explicit ablation: disabling self-similar cadence restores the historical
    // fixed 15 ms period and leaves tau free to change independently.
    f.setTauScaledPseudoUpdateCadence(false);
    const float fixed_rp_var = staged_rp_base * staged_rp_base;
    const float fixed_rv_var = staged_rv_base * staged_rv_base;
    if (!near(f.getPseudoUpdatePeriodSec(), PSEUDO_UPDATE_PERIOD_NOMINAL_S) ||
        !near(f.mekf_->R_p0(2, 2), fixed_rp_var) ||
        !near(f.mekf_->R_v0(2, 2), fixed_rv_var)) {
        std::cerr << "FAIL: fixed-cadence ablation did not restore 15 ms/base r_p0,r_v0\n";
        return 1;
    }
    f.setTauScaledPseudoUpdateCadence(true);
    if (!near(f.getPseudoUpdatePeriodSec(), staged_pseudo) ||
        !near(f.mekf_->R_p0(2, 2), staged_rp_var) ||
        !near(f.mekf_->R_v0(2, 2), staged_rv_var)) {
        std::cerr << "FAIL: re-enabling tau-scaled cadence did not restore compensated schedule\n";
        return 1;
    }

    std::cout << "OU-II predictable tuner scheduling, tau-scaled cadence, and information-rate compensation passed\n";
    return 0;
}
