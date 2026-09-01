#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

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
    f.time_ = 1.0;
    f.last_adapt_time_sec_ = 0.0;
    f.adapt_every_secs_ = 0.05f;

    if (f.getRSLaw() != RSAdaptationLaw::SpectralMSE) {
        std::cerr << "FAIL: OU-III schedule test is not using deployed SpectralMSE\n";
        return 1;
    }
    if (!near(f.getPseudoUpdateTauRatio(), 0.015f / 1.1f) ||
        !near(f.getPseudoUpdatePeriodSec(), 0.015f)) {
        std::cerr << "FAIL: OU-III tau-scaled pseudo cadence did not preserve the nominal point\n";
        return 1;
    }

    // Period estimation and parameter slew are separate layers.  The common
    // parameter EMA keeps its sea-time coefficient, but SeaStateAutoTuner must
    // not add another frequency estimator after WavePeriodEstimator.
    if (!near(f.adapt_tau_sea_periods_, 0.40f)) {
        std::cerr << "FAIL: OU-III parameter-slew default is not 0.40 of T_sea\n";
        return 1;
    }
    {
        SeaStateAutoTuner tuner(2.0f);
        tuner.setFrequencySmoothingSeaPeriods(0.10f); // compatibility no-op
        tuner.update(0.005f, 0.2f, 0.25f);
        if (!near(tuner.getFrequencyHz(), 0.25f) ||
            !near(tuner.getFrequencySmoothingSeaPeriods(), 0.0f) ||
            !near(tuner.getFrequencySmoothingHorizonSec(), 0.0f)) {
            std::cerr << "FAIL: SeaStateAutoTuner reintroduced a second frequency smoother\n";
            return 1;
        }
    }

    // Physical measurements are never cadence-decimated.  Even inside the
    // parameter-activation interval, a valid sample must move the EMA candidate
    // while leaving the active MEKF schedule untouched.
    {
        Filter s(false);
        s.time_ = 2.0;
        s.last_adapt_time_sec_ = 2.0;
        s.online_tune_apply_pending_ = false;
        const float before = s.tune_.tau_applied;
        constexpr float T_sea = 2.5f;
        constexpr float dt_ema = 0.1f;
        const float expected_alpha = 1.0f - std::exp(-dt_ema / (0.40f * T_sea));
        const float expected = before + expected_alpha * (3.0f - before);
        s.adapt_mekf(dt_ema, 3.0f, 0.8f, 1.2f, T_sea);
        if (!near(s.tune_.tau_applied, expected, 1e-5f)) {
            std::cerr << "FAIL: sea-scaled EMA did not consume the physical sample\n";
            return 1;
        }
        if (s.online_tune_apply_pending_) {
            std::cerr << "FAIL: EMA sample incorrectly bypassed activation cadence\n";
            return 1;
        }
    }

    const float active_tau_before = f.mekf_->tau_aw;
    const float active_rs_before = f.mekf_->R_S(2, 2);
    const float active_pseudo_before = f.getPseudoUpdatePeriodSec();

    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 2.5f);
    const float staged_tau = f.tune_.tau_applied;
    const float staged_rs = std::min(std::max(f.tune_.RS_applied, f.min_R_S_), f.max_R_S_);
    const float staged_pseudo = std::min(
        std::max(f.pseudo_update_tau_ratio_ * staged_tau,
                 f.pseudo_update_period_min_s_),
        f.pseudo_update_period_max_s_);
    const float staged_rs_var = staged_rs * staged_rs;

    if (!f.online_tune_apply_pending_) {
        std::cerr << "FAIL: adaptation was not staged\n";
        return 1;
    }
    if (!near(f.mekf_->tau_aw, active_tau_before) ||
        !near(f.mekf_->R_S(2, 2), active_rs_before) ||
        !near(f.getPseudoUpdatePeriodSec(), active_pseudo_before)) {
        std::cerr << "FAIL: current-sample tuner output changed the active OU-III schedule\n";
        return 1;
    }
    if (near(staged_tau, active_tau_before)) {
        std::cerr << "FAIL: test did not create a distinct staged candidate\n";
        return 1;
    }

    f.enable_tuner_ = false;
    f.updateTime(0.005f, gyro, acc);

    if (!near(f.mekf_->tau_aw, staged_tau) ||
        !near(f.mekf_->R_S(2, 2), staged_rs_var) ||
        !near(f.getPseudoUpdatePeriodSec(), staged_pseudo)) {
        std::cerr << "FAIL: staged OU-III SpectralMSE schedule/cadence was not committed at next sample\n";
        return 1;
    }
    if (!near(f.getPseudoUpdatePeriodSec() / f.getTauApplied(),
              f.getPseudoUpdateTauRatio(), 1e-5f)) {
        std::cerr << "FAIL: T_S/tau is not invariant after the staged update\n";
        return 1;
    }

    const float expected_info_product = staged_rs_var * staged_pseudo;
    if (!near(f.mekf_->R_S(2, 2) * f.getPseudoUpdatePeriodSec(),
              expected_info_product, 1e-5f)) {
        std::cerr << "FAIL: SpectralMSE R_S*T_S product moved during schedule commit\n";
        return 1;
    }
    if (f.online_tune_apply_pending_) {
        std::cerr << "FAIL: pending flag survived the commit\n";
        return 1;
    }

    f.setTauScaledPseudoUpdateCadence(false);
    const float fixed_rs_var = staged_rs * staged_rs;
    if (!near(f.getPseudoUpdatePeriodSec(), PSEUDO_UPDATE_PERIOD_NOMINAL_S) ||
        !near(f.mekf_->R_S(2, 2), fixed_rs_var)) {
        std::cerr << "FAIL: fixed-cadence ablation did not restore 15 ms/SpectralMSE R_S\n";
        return 1;
    }
    f.setTauScaledPseudoUpdateCadence(true);
    if (!near(f.getPseudoUpdatePeriodSec(), staged_pseudo) ||
        !near(f.mekf_->R_S(2, 2), staged_rs_var)) {
        std::cerr << "FAIL: re-enabling tau-scaled cadence did not restore SpectralMSE schedule\n";
        return 1;
    }

    std::cout << "OU-III SpectralMSE predictable tuner scheduling, direct canonical frequency, and tau-scaled cadence passed\n";
    return 0;
}
