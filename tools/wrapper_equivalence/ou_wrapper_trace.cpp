// Bit-level characterization trace of the OU-II / OU-III sea-state wrappers.
//
// Built once per family (-DOU_FAMILY=2 or 3) against a source tree, run over
// deterministic synthetic histories, and written as raw float bits.  Two
// builds are behaviorally equivalent on these histories exactly when their
// traces are byte-identical.  Only public API is read, so the same driver
// builds against any revision that keeps the public API.
#define EIGEN_NON_ARDUINO
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <string>
#include <Eigen/Dense>

#if OU_FAMILY == 2
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
template <TrackerType T> using Inner = SeaStateFusionFilter_OU_II<T>;
template <TrackerType T> using Outer = SeaStateFusion_OU_II<T>;
#else
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
template <TrackerType T> using Inner = SeaStateFusionFilter_OU_III<T>;
template <TrackerType T> using Outer = SeaStateFusion_OU_III<T>;
#endif
#include "SyntheticMarineHistory.h"

extern const float g_std = 9.80665f;

using namespace characterization;

// Branch-coverage summary printed with each trace, so a reviewer can see that
// the scenarios reach the orchestration paths being compared.
struct Coverage {
    float max_engagement = 0.0f, max_hi = 0.0f, max_racc = 0.0f, min_still_sigma = 1e9f;
    int stage_changes = 0, last_stage = -1, dir_sign_changes = 0, last_sign = -1;
    bool hi_valid = false;
    int attitude_jumps = 0;
    Eigen::Quaternionf last_q = Eigen::Quaternionf::Identity();
};
static Coverage cov;

template <typename F>
static void put_inner(Trace& tr, const F& f) {
    tr.put(static_cast<int>(f.getStartupStage()));
    cov.max_engagement = std::max(cov.max_engagement, f.accelVibrationGuardEngagement());
    cov.max_racc = std::max(cov.max_racc, f.accelVibrationRaccStd().maxCoeff());
    if (static_cast<int>(f.getStartupStage()) != cov.last_stage) { ++cov.stage_changes; cov.last_stage = static_cast<int>(f.getStartupStage()); }
    if (static_cast<int>(f.getDirSignState()) != cov.last_sign) { ++cov.dir_sign_changes; cov.last_sign = static_cast<int>(f.getDirSignState()); }
    if (f.isAdaptiveLive()) cov.min_still_sigma = std::min(cov.min_still_sigma, f.getSigmaTarget());
    tr.put(f.isTunerReady());
    tr.put(f.isAdaptiveLive());
    tr.put(f.getFreqHz()); tr.put(f.getFreqSlowHz()); tr.put(f.getFreqRawHz());
    tr.put(f.getPeriodSec());
    tr.put(f.getTauApplied()); tr.put(f.getSigmaApplied());
    tr.put(f.getTauTarget()); tr.put(f.getSigmaTarget());
#if OU_FAMILY == 2
    tr.put(f.getR_p0_std_applied()); tr.put(f.getR_v0_std_applied());
    tr.put(f.getR_p0_std_target()); tr.put(f.getR_v0_std_target());
#else
    tr.put(f.getRSApplied()); tr.put(f.getRSTarget());
    tr.put(f.getAccelInformationRatio());
#endif
    tr.put(f.getPseudoUpdatePeriodSec());
    tr.put(f.getAccelVariance()); tr.put(f.getAccelVertical());
    tr.put(f.getHeaveAbs());
    tr.put(f.getDisplacementScale(true)); tr.put(f.getDisplacementScale(false));
    tr.put(f.getVerticalSpeedEnvelopeMps(true));
    tr.put(f.getWavePeriodSec()); tr.put(f.wavePeriodUsable()); tr.put(f.wavePeriodReady());
    tr.put(f.getSigmaVarianceHorizonSec());
    tr.put(f.getSigmaWaveBandLowHz()); tr.put(f.getSigmaWaveBandHighHz());
    tr.put(f.getSigmaBandNoiseStd());
    tr.put(f.getWaveAxisDeg()); tr.put(f.getWaveDirectionDeg());
    tr.put(f.getApparentWaveDirectionToDeg()); tr.put(f.getApparentWaveDirectionFromDeg());
    tr.put(f.getDirSenseCoherence()); tr.put(static_cast<int>(f.getDirSignState()));
    tr.put(f.dir().getLastStableConfidence());
    tr.put(f.accelVibrationGuardEngagement()); tr.put(f.accelVibrationRms());
    tr.put(f.accelVibrationGuardDelaySec());
    tr.put(f.accelVibrationRaccStd());
    tr.put(f.startupProxyQuat()); tr.put(f.startupProxyTiltQuat());
    tr.put(f.startupProxyInitialized());
    tr.put(f.accBiasHeld());

    const auto& m = f.mekf();
    tr.put(m.quaternion_boat());
    tr.put(m.get_velocity()); tr.put(m.get_position());
#if OU_FAMILY == 3
    tr.put(m.get_integral_displacement());
#endif
    tr.put(m.get_world_accel());
    tr.put(m.gyroscope_bias()); tr.put(m.get_acc_bias());
    tr.put(m.get_aw_time_constant()); tr.put(m.get_aw_stationary_std());
    tr.put(m.acc_bias_updates_enabled());
    tr.put(m.lastMagDiag().accepted); tr.put(m.lastMagDiag().r);
    const auto P = m.covariance_full();
    tr.put(P.diagonal());
    tr.put(P(0, 6)); tr.put(P(2, 8)); tr.put(P(6, 9)); tr.put(P(8, 11));
}

template <typename W>
static void put_outer(Trace& tr, const W& w) {
    cov.max_hi = std::max(cov.max_hi, w.magContinuousHardIronAppliedUT().norm());
    cov.hi_valid = cov.hi_valid || w.magContinuousHardIron().estimate().valid;
    if (w.isLive() && cov.last_q.angularDistance(w.attitudeQuat()) > 0.3f) ++cov.attitude_jumps;
    cov.last_q = w.attitudeQuat();
    tr.put(w.isLive()); tr.put(w.hasMagNorthLock()); tr.put(w.hasRefinedMagReference());
    tr.put(w.magRefineTimeSec()); tr.put(w.magNorthLockTimeSec()); tr.put(w.liveTimeSec());
    tr.put(w.magTiltFrameYawDeg()); tr.put(w.magStartupYawCorrectionDeg());
    tr.put(w.magHardIronBodyUT()); tr.put(w.magContinuousHardIronAppliedUT());
    tr.put(w.magContinuousHardIron().estimate().valid);
    tr.put(w.attitudeQuat());
    tr.put(w.freqHz()); tr.put(w.waveDirectionDeg());
    tr.put(w.displacementUpMeters());
    const auto& d = w.displacementDetrend();
    tr.put(d.input); tr.put(d.baseline_slow); tr.put(d.wave_raw); tr.put(d.wave_clean);
#if OU_FAMILY == 3
    tr.put(w.magAcceptedCount()); tr.put(w.magRejectedCount());
    tr.put(w.magAcceptedWindowSec()); tr.put(w.magEffectiveWeight());
#endif
    put_inner(tr, w.raw());
}

template <TrackerType T>
static void run_outer(const std::string& path, bool with_mag, bool detrend,
                      const HistoryConfig& hc, Trace& tr_all) {
    using W = Outer<T>;
    typename W::Config cfg;
    cfg.with_mag = with_mag;
    cfg.enable_displacement_detrend = detrend;
    W w;
    w.begin(cfg);
    Trace tr;
    tr.fp = std::fopen(path.c_str(), "wb");
    SyntheticMarineHistory h(hc);
    Sample s;
    int applied = 0;
    while (h.next(s)) {
        w.update(0.005f, s.gyro, s.acc, 30.0f);
        if (s.mag_due) {
            w.updateMag(s.mag);
            if (w.isLive() && w.raw().mekf().lastMagDiag().accepted) ++applied;
        }
        tr.put(applied);
        put_outer(tr, w);
    }
    std::fclose(tr.fp);
    std::printf("  coverage: guard_max=%.3f racc_max=%.3f hi_valid=%d hi_max=%.3f stage_changes=%d dir_sign_changes=%d min_live_sigma_target=%.4f attitude_jumps=%d\n",
                cov.max_engagement, cov.max_racc, (int)cov.hi_valid, cov.max_hi, cov.stage_changes, cov.dir_sign_changes, cov.min_still_sigma, cov.attitude_jumps);
    cov = Coverage{};
    std::printf("%s values=%llu hash=%016llx live_at=%.3f north_lock=%.3f refine=%.3f applied_mag=%d\n",
                path.c_str(), (unsigned long long)tr.values, (unsigned long long)tr.hash,
                w.liveTimeSec(), w.magNorthLockTimeSec(), w.magRefineTimeSec(), applied);
    tr_all.hash ^= tr.hash;
}

// Inner filter driven directly, with the runtime knobs the outer wrapper never
// touches: fixed tuning, law switches, cadence toggles, sync policy, channel
// freezes and bounds.  Each knob is flipped at a scheduled time.
template <TrackerType T>
static void run_inner(const std::string& path, const HistoryConfig& hc, Trace& tr_all) {
    using F = Inner<T>;
    F f(true);
    f.initialize(Eigen::Vector3f::Constant(0.2f), Eigen::Vector3f::Constant(0.002f),
                 Eigen::Vector3f::Constant(0.3f));
    f.setNominalRaccStd(Eigen::Vector3f::Constant(0.2f));
    f.setMagDelaySec(5.0f);
    f.mekf().set_mag_world_ref(Eigen::Vector3f(20.0f, 0.0f, 44.0f));
    Trace tr;
    tr.fp = std::fopen(path.c_str(), "wb");
    SyntheticMarineHistory h(hc);
    Sample s;
    bool live = false;
    int k = 0;
    while (h.next(s)) {
        ++k;
        if (!live) {
            f.updateFrontEnd(0.005f, s.gyro, s.acc);
            if (f.isTunerReady() && s.t > 30.0) {
                f.goLive(f.startupProxyQuat(), 0.035f, 0.087f, k % 2 == 0);
                live = true;
            }
        } else {
            f.updateTime(0.005f, s.gyro, s.acc, 25.0f);
        }
        if (s.mag_due) f.updateMag(s.mag);

        if (k == 60 * 200) f.setPeriodicAwCovarianceSync(false);
        if (k == 70 * 200) f.setPeriodicAwCovarianceSync(true);
        if (k == 80 * 200) f.setTauScaledPseudoUpdateCadence(false);
        if (k == 90 * 200) f.setTauScaledPseudoUpdateCadence(true);
        if (k == 95 * 200) f.setPseudoUpdateTauRatio(0.02f);
        if (k == 100 * 200) f.setPseudoUpdatePeriodBounds(0.004f, 0.2f);
        if (k == 105 * 200) f.setAdaptationTimeConstants(2.5f);
        if (k == 110 * 200) f.setAdaptationSeaPeriods(0.5f);
        if (k == 115 * 200) f.setAccBiasHold(true);
        if (k == 118 * 200) f.setAccBiasHold(false);
        if (k == 120 * 200) f.setAccelVibrationRaccGain(1.1f);
        if (k == 125 * 200) f.setSigmaWaveBandRatios(0.6f, 3.5f);
        if (k == 130 * 200) f.setTuneFreqBounds(0.04f, 1.0f);
        if (k == 135 * 200) f.setFreqBounds(0.25f, 5.0f);
        if (k == 140 * 200) f.setAccNoiseFloorSigma(0.1f);
        if (k == 142 * 200) f.setWavePeriodInput(WavePeriodInputSource::Leveled);
        if (k == 150 * 200) f.setWavePeriodInput(WavePeriodInputSource::Complementary);
        if (k == 155 * 200) {
            wave_direction::VesselRaoNoiseWeighting nw;
            nw.max_std_scale = 2.0f;
            f.setLowWaveNoiseWeighting(nw);
        }
#if OU_FAMILY == 2
        if (k == 160 * 200) f.setPseudoLaw(PseudoAdaptationLaw::Empirical);
        if (k == 165 * 200) f.setR_p0_Coeff(0.7f);
        if (k == 170 * 200) f.setPseudoLaw(PseudoAdaptationLaw::PhysicalMSE);
        if (k == 175 * 200) f.setR_AdaptSlewLog(0.3f);
        if (k == 180 * 200) f.setFixedTuning(3.0f, 0.4f, 2.0f, 0.8f);
        if (k == 190 * 200) f.enableTuner(true);
        if (k == 195 * 200) f.setSigmaStillnessDecaySec(2.0f);
#else
        if (k == 160 * 200) f.setRSLaw(RSAdaptationLaw::Cubic);
        if (k == 165 * 200) f.setRSCoeff(15.0f);
        if (k == 168 * 200) f.setRSLaw(RSAdaptationLaw::PosteriorRiccati);
        if (k == 170 * 200) f.setRSLaw(RSAdaptationLaw::SpectralMSE);
        if (k == 172 * 200) f.setAwCovarianceSyncCongruent(true);
        if (k == 175 * 200) f.setRSAdaptSlewLog(0.3f);
        if (k == 180 * 200) f.setChannelFreeze(true, 3.0f, 0.4f, false, 0.0f);
        if (k == 185 * 200) f.setChannelFreeze(false, 0.0f, 0.0f, true, 5.0f);
        if (k == 190 * 200) f.setFixedTuning(3.0f, 0.4f, 6.0f);
        if (k == 195 * 200) f.enableTuner(true);
#endif
        put_inner(tr, f);
    }
    std::fclose(tr.fp);
    std::printf("  coverage: guard_max=%.3f racc_max=%.3f stage_changes=%d dir_sign_changes=%d min_live_sigma_target=%.4f\n",
                cov.max_engagement, cov.max_racc, cov.stage_changes, cov.dir_sign_changes, cov.min_still_sigma);
    cov = Coverage{};
    std::printf("%s values=%llu hash=%016llx\n", path.c_str(),
                (unsigned long long)tr.values, (unsigned long long)tr.hash);
    tr_all.hash ^= tr.hash;
}

// The discouraged direct path: updateTime() from the first sample, with a
// capsize while the stage machine is still warming up.  This is the only way
// to reach the pre-Live tilt reset (enterCold_ + resetTrackingState_).
template <TrackerType T>
static void run_inner_direct(const std::string& path, const HistoryConfig& hc, Trace& tr_all) {
    using F = Inner<T>;
    F f(false);
    f.initialize(Eigen::Vector3f::Constant(0.2f), Eigen::Vector3f::Constant(0.002f),
                 Eigen::Vector3f::Constant(0.3f));
    f.setNominalRaccStd(Eigen::Vector3f::Constant(0.2f));
    Trace tr;
    tr.fp = std::fopen(path.c_str(), "wb");
    SyntheticMarineHistory h(hc);
    Sample s;
    while (h.next(s)) {
        f.updateTime(0.005f, s.gyro, s.acc, 25.0f);
        put_inner(tr, f);
    }
    std::fclose(tr.fp);
    std::printf("  coverage: stage_changes=%d\n", cov.stage_changes);
    cov = Coverage{};
    std::printf("%s values=%llu hash=%016llx\n", path.c_str(),
                (unsigned long long)tr.values, (unsigned long long)tr.hash);
    tr_all.hash ^= tr.hash;
}

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : ".";
    const std::string fam = OU_FAMILY == 2 ? "ou2" : "ou3";
    Trace all;

    HistoryConfig full;  // mag, vibration, quiet interval
    run_outer<TrackerType::KALMANF>(out + "/" + fam + "_outer_mag_detrend.bin", true, true, full, all);

    HistoryConfig nomag = full;
    nomag.duration_sec = 260.0;
    run_outer<TrackerType::KALMANF>(out + "/" + fam + "_outer_nomag.bin", false, false, nomag, all);

    HistoryConfig capsize = full;
    capsize.capsize = true;
    capsize.seed = 0xc0ffeeu;
    run_outer<TrackerType::PLL>(out + "/" + fam + "_outer_capsize_pll.bin", true, false, capsize, all);

    HistoryConfig inner = full;
    inner.duration_sec = 220.0;
    inner.seed = 0xabcdefu;
    run_outer<TrackerType::ARANOVSKIY>(out + "/" + fam + "_outer_aranovskiy.bin", true, true, inner, all);
    run_inner<TrackerType::KALMANF>(out + "/" + fam + "_inner_knobs.bin", inner, all);

    HistoryConfig direct = inner;
    direct.duration_sec = 120.0;
    direct.capsize = true;
    direct.capsize_start = 3.0;
    run_inner_direct<TrackerType::ZEROCROSS>(out + "/" + fam + "_inner_direct_cold.bin", direct, all);

    std::printf("%s combined=%016llx\n", fam.c_str(), (unsigned long long)all.hash);
    return 0;
}
