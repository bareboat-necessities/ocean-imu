// Pins the implementation assumptions used by the analytic Live-state stability
// discussion.
//
// doc/kalman_ou_iii/w3d-iss-stability.tex-part treats the default tuning
// schedule as exogenous with respect to the MEKF state. This test fixes the
// implementation facts that statement depends on, so none of them can
// drift silently. State exogeneity is checked here; the separate
// tuner_schedule-test pins the one-sample predictability boundary.
//
// 1. The sigma channel is exogenous.  Its default input is the vertical
//    acceleration levelled by a private Mahony observer that reads only raw
//    gyro and accelerometer.  The sample then passes through the period-scaled
//    AdaptiveWaveBandPass; both its sample and its corner frequency remain
//    functions of measurements only.  Displacing the MEKF state must therefore
//    leave the band output, sigma target and r_S target bit-for-bit unchanged.
//
// 2. The frequency channel is exogenous too, under the default
//    WavePeriodInputSource::Complementary. wave_period_ is fed a vertical
//    acceleration levelled by the same private Mahony observer, so displacing
//    every estimator state leaves the reported period and operating point
//    bit-for-bit unchanged.
//
//    WavePeriodInputSource::Leveled is the older behaviour, the only ablation
//    left on this axis,
//    and it does close the loop: it feeds direction_accel.up_ms2, which uses
//    the attitude estimate, so delta_theta -> wave period -> tau, sigma_a,
//    r_S -> filter gains -> delta_theta is a real feedback path. It reaches
//    the linear states indirectly as well, since displacing v, p, S or a_w
//    perturbs attitude through the filter's cross-covariances. Its gain is
//    still bounded here because the ablation has to stay interpretable.
//
// 3. One gate reads the attitude estimate directly and re-locks the filter
//    above 70 degrees of tilt. The analytic result is scoped to the normal
//    Live region, so the margin to this hybrid re-lock gate is pinned.
#define EIGEN_NON_ARDUINO

// Standard headers first: the access override below must not be in force when
// libstdc++ headers are parsed, or their own access specifiers change meaning.
#include <cmath>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Geometry>

// White-box on purpose. The point is to displace *every* estimator state,
// including ones with no public setter, so that no dependence can hide in the
// states the public API happens not to expose.
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

namespace {

using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

constexpr float DT = 1.0f / 240.0f;
constexpr float G = 9.80665f;

// A wave-like body-frame measurement stream. Only the sample index decides
// what it is, so two filters can be driven with exactly the same inputs.
void sample(int k, Eigen::Vector3f& gyro, Eigen::Vector3f& acc) {
    const float t = static_cast<float>(k) * DT;
    const float w = 2.0f * float(M_PI) / 8.0f;

    const float roll = 0.35f * std::sin(w * t);
    const float pitch = 0.22f * std::sin(w * t + 1.1f);

    gyro = Eigen::Vector3f(0.35f * w * std::cos(w * t),
                           0.22f * w * std::cos(w * t + 1.1f),
                           0.0f);

    const float heave = -1.2f * w * w * std::sin(w * t);
    acc = Eigen::Vector3f(G * std::sin(pitch),
                          -G * std::sin(roll),
                          -G + heave);
}

// Mirrors what SeaStateFusion_OU_III::begin() does to the inner filter.
// Templated because the tracker-independence check below has to bring up two
// filters that differ in nothing but their TrackerType.
template <typename F>
void bring_up(F& f) {
    f.setWithMag(false);
    f.setOnlineTuneWarmupSec(5.0f);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));

    Eigen::Vector3f gyro, acc;
    sample(0, gyro, acc);
    f.initialize_from_acc(acc);
}

// Drives the filter the way the wrapper does: the measurement-only front end
// runs until the tuner is ready, then the proxy attitude is handed over and the
// MEKF takes the samples.  The handoff instant and the seed are both functions
// of the measurements alone, so two filters given identical samples hand over
// on the same sample with the same attitude.
template <typename F>
void run(F& f, int from, int to) {
    Eigen::Vector3f gyro, acc;
    for (int k = from; k < to; ++k) {
        sample(k, gyro, acc);

        if (f.isAdaptiveLive()) {
            f.updateTime(DT, gyro, acc);
            continue;
        }

        f.updateFrontEnd(DT, gyro, acc);
        if (f.isTunerReady()) {
            f.goLive(f.startupProxyQuat(), 0.035f, 1.5708f);
        }
    }
}

constexpr float ATTITUDE_DISPLACEMENT_RAD = 0.25f;

// Displace the estimator state without touching any measurement.
// The linear displacement is kept small so it does not drive the estimated
// tilt into the 70 degree re-lock gate, which would confound the result.
void displace_state(Filter& f, bool attitude, bool linear) {
    auto& mekf = f.mekf();

    if (attitude) {
        Eigen::Quaternionf q = mekf.quaternion_boat();
        const Eigen::Quaternionf tilt(Eigen::AngleAxisf(
            ATTITUDE_DISPLACEMENT_RAD, Eigen::Vector3f(0.3f, 0.9f, 0.32f).normalized()));
        q = tilt * q;
        q.normalize();
        mekf.set_quaternion_boat(q);
    }
    if (linear) {
        for (int i = 3; i < static_cast<int>(mekf.xext.size()); ++i) {
            mekf.xext(i) += 0.02f;
        }
    }
}

constexpr int SETTLE = 240 * 120;
constexpr int AFTER = SETTLE + 240 * 60;

bool test_sigma_channel_input_is_exogenous() {
    bool ok = true;

    Filter nominal, displaced;
    bring_up(nominal);
    bring_up(displaced);
    run(nominal, 0, SETTLE);
    run(displaced, 0, SETTLE);

    ok &= check(nominal.isAdaptiveLive(),
                "the tuner must reach its live stage, or this test proves nothing");
    ok &= check(nominal.wavePeriodReady(),
                "the wave period must be ready, or the sigma band is not on its default reference");

    displace_state(displaced, /*attitude=*/true, /*linear=*/true);
    run(nominal, SETTLE, AFTER);
    run(displaced, SETTLE, AFTER);

    ok &= check(nominal.sigma_wave_band_.output() == displaced.sigma_wave_band_.output(),
                "the sigma-band sample must be a function of measurements alone");
    ok &= check(nominal.sigma_wave_band_.lowHz() == displaced.sigma_wave_band_.lowHz() &&
                nominal.sigma_wave_band_.highHz() == displaced.sigma_wave_band_.highHz(),
                "sigma-band corners must be a function of measurements alone");
    ok &= check(nominal.getSigmaTarget() == displaced.getSigmaTarget(),
                "sigma target must be exogenous to estimator state");
    ok &= check(nominal.getRSTarget() == displaced.getRSTarget(),
                "r_S target must remain exogenous when sigma is wave-band estimated");

    if (!ok) {
        std::cerr << "  nominal band " << nominal.sigma_wave_band_.output()
                  << ", displaced " << displaced.sigma_wave_band_.output() << '\n'
                  << "  sigma target " << nominal.getSigmaTarget()
                  << " vs " << displaced.getSigmaTarget() << '\n'
                  << "  The default sigma channel has acquired a dependence on\n"
                     "  estimator state. The analytic stability discussion assumes\n"
                     "  the applied default schedule is exogenous.\n";
    }
    return ok;
}

bool test_leveled_ablation_coupling_stays_bounded() {
    bool ok = true;

    Filter nominal, displaced;
    bring_up(nominal);
    bring_up(displaced);
    // Not the default any more, so it has to be selected: without this the
    // test would measure the exogenous path and pass while proving nothing.
    nominal.setWavePeriodInput(WavePeriodInputSource::Leveled);
    displaced.setWavePeriodInput(WavePeriodInputSource::Leveled);
    run(nominal, 0, SETTLE);
    run(displaced, 0, SETTLE);

    displace_state(displaced, /*attitude=*/true, /*linear=*/false);
    run(nominal, SETTLE, AFTER);
    run(displaced, SETTLE, AFTER);

    ok &= check(displaced.tilt_over_limit_sec_ == 0.0f,
                "the tilt gate must not arm, or this measures the gate not the loop");

    const float ratio = displaced.getTauTarget() / nominal.getTauTarget();

    // Measured 1.25 for a 0.25 rad displacement. The bound is a regression
    // guard, not a fit: it fails if the attitude-to-tuning path ever becomes
    // strong enough that this ablation stops being a mild perturbation of the
    // default. It is no longer load-bearing for the shipped filter, whose
    // frequency channel is exogenous, but it keeps the comparison honest.
    constexpr float MAX_RATIO = 2.0f;
    ok &= check(ratio > 1.0f / MAX_RATIO && ratio < MAX_RATIO,
                "attitude-to-tuning loop gain must stay modest");

    std::cout << "  attitude displacement " << ATTITUDE_DISPLACEMENT_RAD
              << " rad -> wave period " << nominal.getWavePeriodSec() << " s vs "
              << displaced.getWavePeriodSec() << " s, tau ratio " << ratio << '\n';
    return ok;
}

// The frequency channel under the *default* input source. It joins the sigma
// channel in being exogenous, so the applied schedule remains exogenous with
// respect to estimator state. Asserted bit-for-bit against a displacement of
// *every* state, because "a pure function of the measurements" admits no
// tolerance. Nothing is selected here on purpose: the point is that this is
// what the filter does out of the box.
bool test_default_frequency_channel_is_exogenous() {
    bool ok = true;

    Filter nominal, displaced;
    bring_up(nominal);
    bring_up(displaced);
    ok &= check(nominal.wavePeriodInput() == WavePeriodInputSource::Complementary,
                "the default wave-period input must be the exogenous one");
    run(nominal, 0, SETTLE);
    run(displaced, 0, SETTLE);

    ok &= check(nominal.isAdaptiveLive(),
                "the tuner must reach its live stage, or this test proves nothing");
    ok &= check(nominal.wavePeriodReady(),
                "the period must have settled, or the wave-band prior is what "
                "is being measured");

    displace_state(displaced, /*attitude=*/true, /*linear=*/true);
    run(nominal, SETTLE, AFTER);
    run(displaced, SETTLE, AFTER);

    ok &= check(nominal.getWavePeriodSec() == displaced.getWavePeriodSec(),
                "the default wave period must be a function of the "
                "measurement alone");

    ok &= check(nominal.getTauTarget() == displaced.getTauTarget(),
                "tau target must be a function of the measurement alone");
    ok &= check(nominal.getSigmaTarget() == displaced.getSigmaTarget(),
                "sigma target must be a function of the measurement alone");
    ok &= check(nominal.getRSTarget() == displaced.getRSTarget(),
                "the whole default operating point must be measurement-only");

    std::cout << "  default (complementary): attitude+linear displacement -> "
                 "wave period "
              << nominal.getWavePeriodSec() << " s vs "
              << displaced.getWavePeriodSec() << " s (must be identical)\n";

    if (!ok) {
        std::cerr << "  The private Mahony observer or sigma band has acquired a\n"
                     "  dependence on estimator state. The default tuner exists\n"
                     "  precisely so its schedule cannot close that loop.\n";
    }
    return ok;
}

// The operating point must not depend on the acceleration-band frequency
// tracker at all, at any point in the run -- not only after the wave-period
// estimator reports ready.
//
// Two filters differing in nothing but their TrackerType are driven with
// identical samples.  KalmANF and Aranovskiy converge differently and report
// different frequencies, especially early, so if any part of the adaptation
// path still read the tracker the two schedules would diverge.  Everything the
// two filters share -- the input low-pass, the stillness detector, the private
// Mahony observer, the sigma band -- is a function of the measurement, so the
// only thing this can detect is a tracker dependence.
//
// The comparison is made at 30 s, well before the wave-period estimator is
// ready: that is exactly the window in which the removed tracker-fallback
// source let the tracker set the sigma band's corners, the sigma_a averaging
// horizon and tau, and this check fails on any such path returning.
bool test_operating_point_is_tracker_independent() {
    bool ok = true;

    SeaStateFusionFilter_OU_III<TrackerType::KALMANF> kalmanf;
    SeaStateFusionFilter_OU_III<TrackerType::ARANOVSKIY> aranovskiy;
    bring_up(kalmanf);
    bring_up(aranovskiy);

    constexpr int EARLY = 240 * 30;
    run(kalmanf, 0, EARLY);
    run(aranovskiy, 0, EARLY);

    ok &= check(!kalmanf.wavePeriodReady(),
                "the period estimator must still be unsettled here, or this "
                "measures the settled path and proves nothing about the "
                "window the tracker used to own");
    // Raw, not getFreqHz(): the smoothed branch is clamped into
    // [min_freq_hz_, max_freq_hz_], and this synthetic 8 s wave sits under the
    // 0.2 Hz floor, so both trackers report the floor however much they
    // disagree underneath it.
    ok &= check(kalmanf.getFreqRawHz() != aranovskiy.getFreqRawHz(),
                "the two trackers must actually disagree, or this test cannot "
                "detect a dependence on them");

    ok &= check(kalmanf.getTauTarget() == aranovskiy.getTauTarget(),
                "tau target must not depend on which frequency tracker runs");
    ok &= check(kalmanf.getSigmaTarget() == aranovskiy.getSigmaTarget(),
                "sigma target must not depend on which frequency tracker runs");
    ok &= check(kalmanf.getRSTarget() == aranovskiy.getRSTarget(),
                "r_S target must not depend on which frequency tracker runs");
    ok &= check(kalmanf.sigma_wave_band_.lowHz() == aranovskiy.sigma_wave_band_.lowHz() &&
                kalmanf.sigma_wave_band_.highHz() == aranovskiy.sigma_wave_band_.highHz(),
                "sigma-band corners must not depend on which tracker runs");

    // And again once the estimator has settled, where both sources agree.
    run(kalmanf, EARLY, SETTLE);
    run(aranovskiy, EARLY, SETTLE);
    ok &= check(kalmanf.wavePeriodReady(), "the period estimator must settle");
    ok &= check(kalmanf.getTauTarget() == aranovskiy.getTauTarget() &&
                kalmanf.getSigmaTarget() == aranovskiy.getSigmaTarget() &&
                kalmanf.getRSTarget() == aranovskiy.getRSTarget(),
                "the settled operating point must be tracker-independent too");

    std::cout << "  tracker independence at " << (EARLY / 240) << " s: f_raw "
              << kalmanf.getFreqRawHz() << " vs " << aranovskiy.getFreqRawHz()
              << " Hz, tau " << kalmanf.getTauTarget() << " s both\n";

    if (!ok) {
        std::cerr << "  Part of the adaptation path is reading the "
                     "acceleration-band\n"
                     "  frequency tracker again. The tracker is the wave-direction\n"
                     "  demodulator's carrier; the operating point is a wave-band\n"
                     "  quantity and comes from WavePeriodEstimator.\n";
    }
    return ok;
}

bool test_tilt_gate_stays_outside_the_live_analysis_region() {
    bool ok = true;

    // Representative interior bounds for the Live-state attitude region used
    // to keep this regression check away from the hybrid re-lock boundary.
    constexpr float REGION_ROLL = 0.70f;
    constexpr float REGION_PITCH = 0.45f;

    // The re-lock threshold in SeaStateFusionFilter_OU_III::updateTime.
    constexpr float TILT_RESET_DEG = 70.0f;

    const float worst_deg =
        std::acos(std::cos(REGION_ROLL) * std::cos(REGION_PITCH)) * 57.295779513f;

    ok &= check(worst_deg < TILT_RESET_DEG,
                "the tilt re-lock gate must be unreachable in the Live analysis region");

    // The gate is a discrete state-dependent branch, so retain a useful margin
    // before this regression region reaches it.
    constexpr float REQUIRED_MARGIN_DEG = 10.0f;
    ok &= check(TILT_RESET_DEG - worst_deg > REQUIRED_MARGIN_DEG,
                "tilt margin between the Live analysis region and re-lock gate is too small");

    std::cout << "  tilt margin: Live analysis region " << worst_deg
              << " deg vs gate " << TILT_RESET_DEG << " deg\n";
    return ok;
}

}  // namespace

int main() {
    bool ok = true;
    ok &= test_sigma_channel_input_is_exogenous();
    ok &= test_default_frequency_channel_is_exogenous();
    ok &= test_operating_point_is_tracker_independent();
    ok &= test_leveled_ablation_coupling_stays_bounded();
    ok &= test_tilt_gate_stays_outside_the_live_analysis_region();

    if (!ok) {
        std::cerr << "tuner coupling checks FAILED\n";
        return 1;
    }
    std::cout << "all tuner coupling checks passed\n";
    return 0;
}