// Pins the contract of the two startup-attitude policies.
//
// StagedMekf feeds the MEKF from the first sample and lets it learn its own
// tilt while it runs degraded -- linear block off, accelerometer bias frozen,
// Racc inflated -- and the wrapper reads that attitude back to gate the
// magnetometer and to frame the magnetic world-reference average. Those reads
// are what MahonyProxy removes: the measurement-only front end runs from the
// first sample without the MEKF, the private Mahony observer supplies the
// tilt, and the MEKF is seeded with the finished solution and goes live in one
// step, never occupying the degraded configuration at all.
//
// The facts fixed here are the ones that make that substitution sound:
//
// 1. MahonyProxy is the deployed default, set on the wrapper that can actually
//    perform the handoff, while the inner filter keeps the staged behaviour so
//    that driving it directly still brings it live on its own. The flag really
//    does restore the old path rather than merely renaming it.
//
// 2. The front end is genuinely independent of the MEKF. Driving a filter with
//    updateFrontEnd() and driving an identical one with updateTime() must
//    leave every front-end output bit-for-bit equal, because no consumer in
//    that chain reads a filter state. This is what allows the operating point
//    to converge before the MEKF has run a single sample.
//
// 3. Under MahonyProxy the MEKF really is held: its attitude and linear states
//    must be untouched until the handoff, so nothing it does during the
//    bootstrap can reach the seed.
//
// 4. The handoff installs the proxy attitude and goes straight to Live,
//    skipping the staged warmup, with the linear block on.
//
// 5. The startup observer's tilt has no static error under a constant gyro
//    bias. This is the one place its tuning differs from the vertical
//    channel's, and it is not cosmetic: the vertical channel can leave the
//    integral term at zero because everything downstream of it is high-passed,
//    while an attitude seed keeps that error as a standing roll and pitch bias
//    for the whole run.
#define EIGEN_NON_ARDUINO

#include <cmath>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

namespace {

using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
using Policy = Filter::StartupInitPolicy;
using Stage  = Filter::StartupStage;

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

constexpr float DT = 1.0f / 200.0f;
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

void bring_up(Filter& f, Policy policy) {
    f.setWithMag(false);
    f.setStartupInitPolicy(policy);
    f.setOnlineTuneWarmupSec(5.0f);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
}

// The deployed default is the proxy policy, and it is set where the handoff
// can actually be performed. The inner filter keeps StagedMekf so that driving
// it directly still brings it live on its own, with nothing above it to hand
// over the attitude.
bool test_defaults_are_split_by_who_can_hand_over() {
    bool ok = true;

    SeaStateFusion_OU_III<TrackerType::KALMANF>::Config cfg;
    ok &= check(cfg.startup_init_policy == Policy::MahonyProxy,
                "MahonyProxy must be the deployed default");

    Filter f;
    ok &= check(f.startupInitPolicy() == Policy::StagedMekf,
                "a directly driven filter must still go live on its own");

    return ok;
}

// The front end reads no filter state, so withholding the MEKF cannot change
// any of its outputs. Anything less than bit-for-bit here would mean a
// consumer in that chain is looking at the estimator after all.
bool test_front_end_is_independent_of_the_mekf() {
    bool ok = true;

    Filter driven, held;
    bring_up(driven, Policy::MahonyProxy);
    bring_up(held, Policy::MahonyProxy);

    Eigen::Vector3f gyro, acc;
    sample(0, gyro, acc);
    driven.initialize_from_acc(acc);
    held.initialize_from_acc(acc);

    constexpr int STEPS = 200 * 400;
    for (int k = 0; k < STEPS; ++k) {
        sample(k, gyro, acc);
        driven.updateTime(DT, gyro, acc);
        held.updateFrontEnd(DT, gyro, acc);
    }

    ok &= check(driven.getFreqHz() == held.getFreqHz(),
                "tracker frequency must not depend on the MEKF running");
    ok &= check(driven.getFreqSlowHz() == held.getFreqSlowHz(),
                "slow frequency must not depend on the MEKF running");
    ok &= check(driven.getAccelVariance() == held.getAccelVariance(),
                "sigma-channel variance must not depend on the MEKF running");

    const float pd = driven.getWavePeriodSec();
    const float ph = held.getWavePeriodSec();
    ok &= check((std::isnan(pd) && std::isnan(ph)) || pd == ph,
                "wave period must not depend on the MEKF running");

    ok &= check(driven.getTauTarget() == held.getTauTarget(),
                "tau target must not depend on the MEKF running");
    ok &= check(driven.getSigmaTarget() == held.getSigmaTarget(),
                "sigma target must not depend on the MEKF running");

    return ok;
}

// The wave-direction stage is the one front-end consumer that needs an
// attitude, and while the MEKF is held it is given the proxy's. That is only
// sound because the levelled heading frame is invariant under q -> Rz(psi) q:
// the proxy's yaw is unobservable and drifts, so anything that could see it
// would be reading noise. Checked on the frame conversion directly, since this
// is a property of that function rather than of either estimator.
bool test_heading_frame_ignores_yaw() {
    bool ok = true;

    const Eigen::Quaternionf q_tilt =
        Eigen::Quaternionf(Eigen::AngleAxisf(0.21f, Eigen::Vector3f::UnitY())) *
        Eigen::Quaternionf(Eigen::AngleAxisf(-0.34f, Eigen::Vector3f::UnitX()));

    const Eigen::Vector3f f_body(0.83f, -1.42f, -G + 0.66f);

    const auto base = wave_direction::heading_frame_acceleration<float>(
        q_tilt, f_body, G);
    ok &= check(base.heading_valid, "heading frame must resolve for a tilted attitude");

    for (float psi : {-2.9f, -1.1f, 0.4f, 2.2f}) {
        const Eigen::Quaternionf q_yawed =
            Eigen::Quaternionf(Eigen::AngleAxisf(psi, Eigen::Vector3f::UnitZ())) * q_tilt;

        const auto yawed = wave_direction::heading_frame_acceleration<float>(
            q_yawed, f_body, G);

        const float df = std::fabs(yawed.forward_ms2 - base.forward_ms2);
        const float ds = std::fabs(yawed.starboard_ms2 - base.starboard_ms2);
        const float du = std::fabs(yawed.up_ms2 - base.up_ms2);

        ok &= check(df < 1e-5f && ds < 1e-5f && du < 1e-5f,
                    "levelled heading frame must be invariant to the yaw of its input");
    }

    return ok;
}

// Nothing the MEKF does during the bootstrap may reach the seed, so it must
// not be running at all.
bool test_proxy_policy_holds_the_mekf() {
    bool ok = true;

    Filter f;
    bring_up(f, Policy::MahonyProxy);

    Eigen::Vector3f gyro, acc;
    sample(0, gyro, acc);
    f.initialize_from_acc(acc);

    const Eigen::Quaternionf q0 = f.mekf().quaternion_boat();

    constexpr int STEPS = 200 * 400;
    for (int k = 0; k < STEPS; ++k) {
        sample(k, gyro, acc);
        f.updateFrontEnd(DT, gyro, acc);
    }

    const Eigen::Quaternionf q1 = f.mekf().quaternion_boat();
    ok &= check(std::fabs(q0.dot(q1)) > 1.0f - 1e-6f,
                "MEKF attitude must not move while the bootstrap holds it");
    ok &= check(f.mekf().get_position().norm() == 0.0f,
                "MEKF position must not move while the bootstrap holds it");
    ok &= check(f.mekf().get_velocity().norm() == 0.0f,
                "MEKF velocity must not move while the bootstrap holds it");

    // The tuner is allowed -- required, in fact -- to reach its ready state
    // during the bootstrap, but it must stop short of Live.
    ok &= check(f.getStartupStage() == Stage::TunerReady,
                "front-end-only running must park the filter at TunerReady");
    ok &= check(f.isTunerReady(), "tuner must be ready after the bootstrap");
    ok &= check(!f.isAdaptiveLive(),
                "the filter must not go Live without an explicit handoff");

    // The proxy has an attitude to hand over.
    ok &= check(f.startupProxyInitialized(),
                "startup proxy must be initialized after the bootstrap");

    // Handoff: the MEKF takes the proxy attitude and goes live in one step.
    const Eigen::Quaternionf q_seed = f.startupProxyQuat();
    f.goLive(q_seed, 0.035f, 1.5708f);

    ok &= check(f.isAdaptiveLive(), "goLive must put the filter in Live");
    ok &= check(f.mekf().linear_block_enabled(),
                "the linear block must be on once the filter is live");

    const Eigen::Quaternionf q_after = f.mekf().quaternion_boat();
    ok &= check(std::fabs(q_seed.normalized().dot(q_after)) > 1.0f - 1e-5f,
                "the live attitude must be the one the bootstrap handed over");

    return ok;
}

// The flag has to restore the old path, not just rename it.
bool test_staged_policy_still_warms_the_mekf() {
    bool ok = true;

    Filter f;
    bring_up(f, Policy::StagedMekf);

    Eigen::Vector3f gyro, acc;
    sample(0, gyro, acc);
    f.initialize_from_acc(acc);

    const Eigen::Quaternionf q0 = f.mekf().quaternion_boat();

    constexpr int STEPS = 200 * 400;
    for (int k = 0; k < STEPS; ++k) {
        sample(k, gyro, acc);
        f.updateTime(DT, gyro, acc);
    }

    ok &= check(f.isAdaptiveLive(),
                "the staged policy must reach Live on its own");
    ok &= check(f.getStartupStage() != Stage::TunerReady,
                "the staged policy must never occupy TunerReady");

    const Eigen::Quaternionf q1 = f.mekf().quaternion_boat();
    ok &= check(std::fabs(q0.dot(q1)) < 1.0f - 1e-6f,
                "the staged policy must let the MEKF propagate its own attitude");

    return ok;
}

// A constant gyro bias must not leave a standing tilt error in the seed.
//
// With no integral term a Mahony observer settles at about 2b/two_kp of tilt
// error, which at two_kp = 0.2 turns a 0.05 deg/s bias into roughly half a
// degree -- invisible to the vertical channel, which high-passes it away, and
// a permanent roll/pitch bias to anything seeded from it.
bool test_startup_proxy_rejects_gyro_bias() {
    bool ok = true;

    // Level platform, still water, constant bias on both horizontal axes.
    const float bias_dps = 0.05f;
    const float bias_rps = bias_dps * float(M_PI) / 180.0f;
    const Eigen::Vector3f gyro_biased(bias_rps, bias_rps, 0.0f);
    const Eigen::Vector3f acc_level(0.0f, 0.0f, -G);

    // Settled tilt of a startup proxy at the given gains.
    auto settled_tilt_deg = [&](float two_kp, float two_ki) {
        Filter f;
        bring_up(f, Policy::MahonyProxy);
        f.setStartupProxyGains(two_kp, two_ki);

        Eigen::Vector3f gyro, acc;
        sample(0, gyro, acc);
        f.initialize_from_acc(acc);

        constexpr int STEPS = 200 * 600;
        for (int k = 0; k < STEPS; ++k) {
            f.updateFrontEnd(DT, gyro_biased, acc_level);
        }

        const Eigen::Quaternionf q = f.startupProxyQuat();
        const Eigen::Vector3f down_body =
            q.conjugate() * Eigen::Vector3f(0.0f, 0.0f, 1.0f);
        float c = down_body.normalized().dot(Eigen::Vector3f(0.0f, 0.0f, 1.0f));
        c = std::max(-1.0f, std::min(1.0f, c));
        return std::acos(c) * 57.295779513f;
    };

    // The vertical channel's tuning, for contrast. Without an integral term a
    // Mahony observer settles at about 2b/two_kp of tilt error, so this is the
    // error the seed would inherit if the bootstrap simply reused that
    // instance -- and it is the measured cost, not a hypothetical one.
    const float tilt_no_integral =
        settled_tilt_deg(STARTUP_PROXY_TWO_KP_DEFAULT, 0.0f);

    const float tilt_default =
        settled_tilt_deg(STARTUP_PROXY_TWO_KP_DEFAULT,
                         STARTUP_PROXY_TWO_KI_DEFAULT);

    ok &= check(tilt_no_integral > 0.2f,
                "zero-integral gains must show the static tilt this guards against");
    ok &= check(tilt_default < 0.05f,
                "startup proxy must not hold a static tilt under a gyro bias");
    ok &= check(tilt_default < 0.1f * tilt_no_integral,
                "the integral term must remove most of the static tilt");

    std::cout << "  static tilt under " << bias_dps << " deg/s bias: "
              << tilt_no_integral << " deg (two_ki=0) vs "
              << tilt_default << " deg (default)\n";

    return ok;
}

}  // namespace

int main() {
    bool ok = true;
    ok &= test_defaults_are_split_by_who_can_hand_over();
    ok &= test_front_end_is_independent_of_the_mekf();
    ok &= test_heading_frame_ignores_yaw();
    ok &= test_proxy_policy_holds_the_mekf();
    ok &= test_staged_policy_still_warms_the_mekf();
    ok &= test_startup_proxy_rejects_gyro_bias();

    if (!ok) {
        std::cerr << "startup init checks FAILED\n";
        return 1;
    }
    std::cout << "all startup init checks passed\n";
    return 0;
}
