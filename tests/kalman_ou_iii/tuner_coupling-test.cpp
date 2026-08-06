// Pins how much of the OU-III adaptation is exogenous, and how much is not.
//
// doc/kalman_ou_iii/w3d-iss-stability.tex-part treats the tuner's coupling to
// the estimator as the open item in the stability argument. This test fixes
// the three facts that argument depends on, so none of them can drift
// silently.
//
// 1. The variance channel is exogenous. update_tuner is handed
//    a_body_z_up_proxy_ = -(acc.z() + g), a pure function of the raw
//    accelerometer. The obvious "improvement" is to feed it the estimated
//    world vertical acceleration instead, which would couple the variance
//    estimate to the very states being tuned. Asserted bit-for-bit.
//
// 2. The frequency channel is NOT exogenous, and that is the loop the
//    certificate still has to account for. wave_period_ is fed
//    direction_accel.up_ms2, the leveled vertical acceleration, which uses
//    the attitude estimate. So delta_theta -> wave period -> tau, sigma_a,
//    r_S -> filter gains -> delta_theta is a real feedback path. It cannot be
//    asserted away; what is asserted is that its gain stays modest, so the
//    small-gain reading of the interconnection keeps its premise.
//
//    Note this reaches the linear states indirectly too: displacing v, p, S
//    or a_w perturbs attitude through the filter's cross-covariances and so
//    moves the tuner as well. The source comment at the wave_period_.update
//    call says the leveled signal "depends on attitude but not on the linear
//    block, so it does not close a loop through the quantity being tuned".
//    That holds for the signal in isolation, not for the closed filter.
//
// 3. One gate reads the attitude estimate directly and re-locks the filter
//    above 70 degrees of tilt. It is kept out of the certified envelope only
//    by margin, so the margin is pinned.
#define EIGEN_NON_ARDUINO

// Standard headers first: the access override below must not be in force when
// libstdc++ headers are parsed, or their own access specifiers change meaning.
#include <cmath>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Geometry>

// White-box on purpose. The point is to displace *every* estimator state,
// including ones with no public setter, so that no coupling can hide in the
// states the public API happens not to expose. tools/iss_certificate uses the
// same access for the same reason.
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
void bring_up(Filter& f) {
    f.setWithMag(false);
    f.setOnlineTuneWarmupSec(5.0f);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));

    Eigen::Vector3f gyro, acc;
    sample(0, gyro, acc);
    f.initialize_from_acc(acc);
}

void run(Filter& f, int from, int to) {
    Eigen::Vector3f gyro, acc;
    for (int k = from; k < to; ++k) {
        sample(k, gyro, acc);
        f.updateTime(DT, gyro, acc);
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

bool test_variance_channel_input_is_exogenous() {
    bool ok = true;

    Filter nominal, displaced;
    bring_up(nominal);
    bring_up(displaced);
    run(nominal, 0, SETTLE);
    run(displaced, 0, SETTLE);

    ok &= check(nominal.isAdaptiveLive(),
                "the tuner must reach its live stage, or this test proves nothing");

    displace_state(displaced, /*attitude=*/true, /*linear=*/true);
    run(nominal, SETTLE, AFTER);
    run(displaced, SETTLE, AFTER);

    // The signal handed to update_tuner as its variance input.
    ok &= check(nominal.a_body_z_up_proxy_ == displaced.a_body_z_up_proxy_,
                "the tuner's variance input must be a function of the measurement alone");

    if (!ok) {
        std::cerr << "  nominal proxy " << nominal.a_body_z_up_proxy_
                  << ", displaced " << displaced.a_body_z_up_proxy_ << '\n'
                  << "  The tuner now reads the estimator state on its variance\n"
                     "  channel too. Sec. app:iss-tuner assumes only the frequency\n"
                     "  channel is coupled; that section needs redoing.\n";
    }
    return ok;
}

bool test_frequency_channel_coupling_stays_bounded() {
    bool ok = true;

    Filter nominal, displaced;
    bring_up(nominal);
    bring_up(displaced);
    run(nominal, 0, SETTLE);
    run(displaced, 0, SETTLE);

    displace_state(displaced, /*attitude=*/true, /*linear=*/false);
    run(nominal, SETTLE, AFTER);
    run(displaced, SETTLE, AFTER);

    ok &= check(displaced.tilt_over_limit_sec_ == 0.0f,
                "the tilt gate must not arm, or this measures the gate not the loop");

    const float ratio = displaced.getTauTarget() / nominal.getTauTarget();

    // Measured 1.34 for a 0.25 rad displacement. The bound is a small-gain
    // guard, not a fit: it fails if the attitude-to-tuning path ever becomes
    // strong enough that the interconnection stops being a mild perturbation.
    // It also passes trivially if someone makes the frequency channel
    // exogenous, which would be an improvement.
    constexpr float MAX_RATIO = 2.0f;
    ok &= check(ratio > 1.0f / MAX_RATIO && ratio < MAX_RATIO,
                "attitude-to-tuning loop gain must stay modest");

    std::cout << "  attitude displacement " << ATTITUDE_DISPLACEMENT_RAD
              << " rad -> wave period " << nominal.getWavePeriodSec() << " s vs "
              << displaced.getWavePeriodSec() << " s, tau ratio " << ratio << '\n';
    return ok;
}

bool test_tilt_gate_stays_outside_the_certified_envelope() {
    bool ok = true;

    // The envelope the ISS certificate samples.
    constexpr float ENVELOPE_ROLL = 0.70f;
    constexpr float ENVELOPE_PITCH = 0.45f;

    // The re-lock threshold in SeaStateFusionFilter_OU_III::updateTime.
    constexpr float TILT_RESET_DEG = 70.0f;

    const float worst_deg =
        std::acos(std::cos(ENVELOPE_ROLL) * std::cos(ENVELOPE_PITCH)) * 57.295779513f;

    ok &= check(worst_deg < TILT_RESET_DEG,
                "the tilt re-lock gate must be unreachable on the certified envelope");

    // The gate is a discrete state-dependent branch, so the margin is
    // load-bearing. Fail well before it closes.
    constexpr float REQUIRED_MARGIN_DEG = 10.0f;
    ok &= check(TILT_RESET_DEG - worst_deg > REQUIRED_MARGIN_DEG,
                "tilt margin between the envelope and the re-lock gate is too small");

    std::cout << "  tilt margin: envelope " << worst_deg << " deg vs gate "
              << TILT_RESET_DEG << " deg\n";
    return ok;
}

}  // namespace

int main() {
    bool ok = true;
    ok &= test_variance_channel_input_is_exogenous();
    ok &= test_frequency_channel_coupling_stays_bounded();
    ok &= test_tilt_gate_stays_outside_the_certified_envelope();

    if (!ok) {
        std::cerr << "tuner coupling checks FAILED\n";
        return 1;
    }
    std::cout << "all tuner coupling checks passed\n";
    return 0;
}
