// Pins predictable OU-II scheduling under the deployed PhysicalMSE law while
// retaining the generic tuner/EMA/cadence invariants that are independent of
// the retired empirical schedule.
#define EIGEN_NON_ARDUINO
#include <algorithm>
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

static int fail(const char* message) {
    std::cerr << "FAIL: " << message << "\n";
    return 1;
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

    if (f.getPseudoLaw() != PseudoAdaptationLaw::PhysicalMSE)
        return fail("OU-II schedule test is not using deployed PhysicalMSE");
    if (!near(f.getPseudoUpdateTauRatio(), 0.015f / 1.1f) ||
        !near(f.getPseudoUpdatePeriodSec(), 0.015f))
        return fail("OU-II tau-scaled pseudo cadence did not preserve nominal point");
    if (!near(f.getAdaptationSeaPeriods(), 0.40f))
        return fail("OU-II parameter-slew default is not 0.40 of T_sea");

    // Generic dynamic-horizon safety guards are part of the scheduler, not of
    // any particular pseudo-measurement law.
    using namespace seastate::tuner::limits;
    if (!near(clampDynamicEmaTimeScaleSec(1.15f), 1.15f) ||
        !near(clampDynamicEmaTimeScaleSec(4.20f), 4.20f))
        return fail("reference sea-time envelope touches the safety clamp");
    if (!near(clampDynamicEmaTimeScaleSec(0.01f), kDynamicEmaTimeScaleMinSec) ||
        !near(clampDynamicEmaTimeScaleSec(100.0f), kDynamicEmaTimeScaleMaxSec))
        return fail("dynamic EMA time-scale clamp does not catch excursions");
    if (!near(seastate::common::adaptiveSmoothingHorizonSec(
                  0.001f, 0.001f, 1.0f, 1.0f, 0.0f, 0.005f),
              kDynamicEmaHorizonMinSec) ||
        !near(seastate::common::adaptiveSmoothingHorizonSec(
                  100.0f, 100.0f, 1.0f, 1.0f, 0.0f, 0.005f),
              kDynamicEmaHorizonMaxSec))
        return fail("dynamic EMA final horizon guard is not universal");

    // The auto-tuner consumes the canonical wave-period frequency directly;
    // only the acceleration-moment horizon is dynamically bounded.
    {
        SeaStateAutoTuner fast(2.0f);
        fast.setFrequencySmoothingSeaPeriods(0.10f);  // compatibility no-op
        fast.update(0.005f, 0.0f, 100.0f);
        if (!near(fast.getFrequencyHz(), 5.0f) ||
            !near(fast.getFrequencySmoothingHorizonSec(), 0.0f) ||
            !near(fast.getVarianceHorizonSec(), 2.0f))
            return fail("short-period tuner direct-frequency/variance guard failed");

        SeaStateAutoTuner slow(2.0f);
        slow.update(0.005f, 0.0f, 0.001f);
        if (!near(slow.getFrequencyHz(), 0.05f) ||
            !near(slow.getVarianceHorizonSec(), 24.0f))
            return fail("long-period tuner direct-frequency/variance guard failed");
    }

    // There must be no hidden second variance EMA.
    {
        SeaStateAutoTuner tuner(1.0f);
        for (int i = 0; i < 5000; ++i) {
            const float sample = (i & 1) ? 2.0f : -1.0f;
            tuner.update(0.005f, sample, 0.5f);
        }
        const float mu = tuner.A_mean.get();
        const float expected = std::max(0.0f, tuner.A_sq.get() - mu * mu);
        if (!near(tuner.getAccelVariance(), expected, 1e-6f))
            return fail("acceleration variance contains a hidden second smoother");
    }

    // Every valid physical sample updates the EMA even when activation is not
    // yet due; the activation cadence only controls when staged parameters are
    // committed to the MEKF.
    {
        Filter sample(false);
        sample.time_ = 2.0;
        sample.last_adapt_time_sec_ = 2.0;
        sample.online_tune_apply_pending_ = false;
        const float before = sample.tune_.tau_applied;
        constexpr float T_sea = 2.5f;
        constexpr float dt_ema = 0.1f;
        const float alpha = 1.0f - std::exp(-dt_ema / (0.40f * T_sea));
        const float expected = before + alpha * (3.0f - before);
        sample.adapt_mekf(dt_ema, 3.0f, 0.8f, 1.2f, 0.7f, T_sea);
        if (!near(sample.tune_.tau_applied, expected, 1e-5f))
            return fail("OU-II sea-scaled EMA did not consume the physical sample");
        if (sample.online_tune_apply_pending_)
            return fail("OU-II EMA sample incorrectly bypassed activation cadence");
    }

    const float active_tau = f.mekf_->tau_aw;
    const float active_rp = f.mekf_->R_p0(2, 2);
    const float active_rv = f.mekf_->R_v0(2, 2);
    const float active_period = f.getPseudoUpdatePeriodSec();

    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 0.7f, 2.5f);
    const float staged_tau = f.tune_.tau_applied;
    const float staged_rp = std::min(
        std::max(f.tune_.R_p0_std_applied, f.MIN_R_p0_std_), f.MAX_R_p0_std_);
    const float staged_rv = std::min(
        std::max(f.tune_.R_v0_std_applied, f.MIN_R_v0_std_), f.MAX_R_v0_std_);
    const float staged_period = std::min(
        std::max(f.pseudo_update_tau_ratio_ * staged_tau,
                 f.pseudo_update_period_min_s_),
        f.pseudo_update_period_max_s_);

    if (!f.online_tune_apply_pending_)
        return fail("adaptation was not staged");
    if (!near(f.mekf_->tau_aw, active_tau) ||
        !near(f.mekf_->R_p0(2, 2), active_rp) ||
        !near(f.mekf_->R_v0(2, 2), active_rv) ||
        !near(f.getPseudoUpdatePeriodSec(), active_period))
        return fail("current-sample tuner output changed the active OU-II schedule");
    if (near(staged_tau, active_tau))
        return fail("test did not create a distinct staged candidate");

    f.enable_tuner_ = false;
    f.updateTime(0.005f, gyro, acc);
    if (!near(f.mekf_->tau_aw, staged_tau) ||
        !near(f.mekf_->R_p0(2, 2), staged_rp * staged_rp) ||
        !near(f.mekf_->R_v0(2, 2), staged_rv * staged_rv) ||
        !near(f.getPseudoUpdatePeriodSec(), staged_period))
        return fail("staged OU-II PhysicalMSE schedule/cadence was not committed");
    if (!near(f.getPseudoUpdatePeriodSec() / f.getTauApplied(),
              f.getPseudoUpdateTauRatio(), 1e-5f))
        return fail("T_S/tau is not invariant after the staged update");

    const float expected_rp_info = staged_rp * staged_rp * staged_period;
    const float expected_rv_info = staged_rv * staged_rv * staged_period;
    if (!near(f.mekf_->R_p0(2, 2) * f.getPseudoUpdatePeriodSec(),
              expected_rp_info, 1e-5f) ||
        !near(f.mekf_->R_v0(2, 2) * f.getPseudoUpdatePeriodSec(),
              expected_rv_info, 1e-5f))
        return fail("PhysicalMSE r^2*T_S information product moved during commit");
    if (f.online_tune_apply_pending_)
        return fail("pending flag survived the commit");

    f.setTauScaledPseudoUpdateCadence(false);
    if (!near(f.getPseudoUpdatePeriodSec(), PSEUDO_UPDATE_PERIOD_NOMINAL_S))
        return fail("fixed-cadence ablation did not restore 15 ms");
    f.setTauScaledPseudoUpdateCadence(true);
    if (!near(f.getPseudoUpdatePeriodSec(), staged_period))
        return fail("re-enabling tau-scaled cadence did not restore schedule");

    std::cout << "OU-II PhysicalMSE scheduling and generic tuner invariants passed\n";
    return 0;
}
