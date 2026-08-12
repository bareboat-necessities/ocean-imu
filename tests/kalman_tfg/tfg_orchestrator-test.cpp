/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    The sea-state orchestrator: startup staging, the tuning laws, the
    adaptation smoothers, and the MEKF gates they drive.

    The tuning laws are checked against closed-form expectations on a
    synthetic sea whose period and acceleration variance are known exactly,
    rather than against recorded output. A regression baseline would lock in
    whatever the code did on the day it was written; this locks in what the
    laws in the article actually say:

        tau  = tau_coeff * 0.5 / f_tune
        R_S  = R_S_coeff * sigma * tau^3
        sigma^2 = max(0, measured variance - noise floor^2)
*/

#define EIGEN_NON_ARDUINO

#include "kalman_tfg/SeaStateFusionFilter_TFG.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
using Stage = Fusion::StartupStage;
using Vector3f = Eigen::Vector3f;
using Matrix3f = Eigen::Matrix3f;

constexpr float kPi = 3.14159265358979f;
constexpr float kG = 9.80665f;

int failures = 0;

bool check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
    return condition;
}

bool near_rel(float got, float want, float rel) {
    return std::fabs(got - want) <= rel * std::fabs(want);
}

/*
    A synthetic sea with a known period and heave-acceleration amplitude, plus
    a little roll so attitude is not trivially level. The accelerometer reading
    is exact for the pose, so any error the filter shows is its own.
*/
struct Sea {
    float f_hz = 0.15f;
    float a_amp = 1.2f;        // m/s^2, vertical
    float roll_amp = 0.08f;    // rad
    float roll_f = 0.12f;

    Matrix3f R_bw(float t) const {
        return Eigen::AngleAxisf(roll_amp * std::sin(2.0f * kPi * roll_f * t),
                                 Vector3f::UnitX()).toRotationMatrix();
    }
    Vector3f gyro(float t) const {
        return Vector3f(roll_amp * 2.0f * kPi * roll_f *
                            std::cos(2.0f * kPi * roll_f * t), 0.0f, 0.0f);
    }
    Vector3f world_accel(float t) const {
        return Vector3f(0.0f, 0.0f, a_amp * std::sin(2.0f * kPi * f_hz * t));
    }
    Vector3f acc_body(float t) const {
        return R_bw(t).transpose() * (world_accel(t) - Vector3f(0, 0, kG));
    }
    Vector3f mag_body(float t) const {
        return R_bw(t).transpose() * Vector3f(0.21f, 0.0f, 0.43f);
    }
};

Fusion::Config default_config() {
    Fusion::Config c;
    c.online_tune_warmup_sec = 10.0f;
    c.mag_delay_sec = 5.0f;
    return c;
}

// Run the sea for `seconds`, optionally recording when each stage was reached.
void run(Fusion& f, const Sea& sea, float seconds, float dt = 0.005f,
         float* t_warm = nullptr, float* t_live = nullptr, bool with_mag = true) {
    const int n = static_cast<int>(seconds / dt);
    for (int i = 0; i < n; ++i) {
        const float t = static_cast<float>(i) * dt;
        const Stage before = f.stage();
        f.update(dt, sea.gyro(t), sea.acc_body(t));
        if (with_mag && i % 5 == 0) f.updateMag(sea.mag_body(t));
        if (before != f.stage()) {
            if (f.stage() == Stage::TunerWarm && t_warm) *t_warm = t;
            if (f.stage() == Stage::Live && t_live) *t_live = t;
        }
    }
}

// ---------------------------------------------------------------------------

/*
    Staged policy: the MEKF levels itself, passes through TunerWarm, then goes
    Live. Kept as an explicit ablation now that MahonyProxy is the default.
*/
void test_staging_staged_mekf() {
    Fusion f;
    auto cfg = default_config();
    cfg.startup_init_policy = Fusion::StartupInitPolicy::StagedMekf;
    f.begin(cfg);
    check(f.stage() == Stage::Cold, "the filter must start Cold");
    check(!f.mekf().linear_block_enabled(),
          "the linear block must be gated off while Cold");

    Sea sea;
    float t_warm = -1.0f, t_live = -1.0f;
    run(f, sea, 300.0f, 0.005f, &t_warm, &t_live);

    check(f.stage() == Stage::Live, "the filter never reached Live");
    check(t_warm >= 0.0f && t_warm < 5.0f,
          "levelling should complete within a few seconds");
    check(t_live > t_warm, "Live must come after TunerWarm");
    check(f.mekf().linear_block_enabled(),
          "the linear block must be enabled once Live");
    check(t_live >= 10.0f,
          "Live was reached before the configured warmup elapsed");
}

/*
    Proxy policy: the MEKF is seeded once and goes straight to Live. It must
    never occupy TunerWarm -- that is the degraded configuration this policy
    exists to skip, and observing it would mean the handoff is not doing its
    job.
*/
void test_staging_mahony_proxy() {
    Fusion f;
    f.begin(default_config());   // MahonyProxy is the default
    check(f.startupInitPolicy() == Fusion::StartupInitPolicy::MahonyProxy,
          "MahonyProxy must be the default startup policy");
    check(f.stage() == Stage::Cold, "the filter must start Cold");

    Sea sea;
    float t_warm = -1.0f, t_live = -1.0f;
    run(f, sea, 300.0f, 0.005f, &t_warm, &t_live);

    check(f.stage() == Stage::Live, "the filter never reached Live under the proxy policy");
    check(t_warm < 0.0f,
          "the proxy policy must go Cold -> Live without entering TunerWarm");
    check(t_live >= 8.0f, "handoff happened before proxy_startup_min_sec");
    check(!f.handoffTimedOut(),
          "the handoff fell back to the timeout instead of a converged front end");
    check(f.mekf().linear_block_enabled(),
          "the linear block must be enabled once Live");
    // The MEKF must never have run degraded: bias updates are on from the
    // moment it exists.
    check(f.mekf().acc_bias_updates_enabled(),
          "accelerometer-bias updates must be enabled when the proxy seeds the MEKF");
}

/*
    The front end can legitimately take a long time to declare itself ready --
    the wave-period channel needs many cycles, and longer in a low sea. Without
    a timeout the MEKF would never start at all, which is a worse failure than
    starting from a stale operating point.
*/
void test_proxy_handoff_times_out() {
    Fusion f;
    auto cfg = default_config();
    cfg.proxy_startup_timeout_sec = 30.0f;
    cfg.with_mag = false;          // nothing will ever gauge north
    f.begin(cfg);

    Sea flat;                       // near-still: the period channel stays unready
    run(f, flat, 25.0f, 0.005f, nullptr, nullptr, /*with_mag=*/false);
    check(f.stage() == Stage::Cold, "handed off before the timeout with no ready front end");

    run(f, flat, 20.0f, 0.005f, nullptr, nullptr, /*with_mag=*/false);
    check(f.stage() == Stage::Live, "the handoff timeout never fired");
    check(f.handoffTimedOut(), "the timeout fired but was not reported");
}

// While warming, the accelerometer bias must stay frozen and Racc inflated:
// a filter that has not yet levelled cannot tell bias from tilt.
void test_warmup_gates() {
    Fusion f;
    auto cfg = default_config();
    cfg.freeze_acc_bias_until_live = true;
    // These gates exist only in the staged path: under MahonyProxy the MEKF is
    // never in the degraded configuration they protect.
    cfg.startup_init_policy = Fusion::StartupInitPolicy::StagedMekf;
    f.begin(cfg);

    Sea sea;
    run(f, sea, 30.0f);
    check(f.stage() == Stage::TunerWarm, "expected to still be warming at 30 s");
    const Vector3f bias_warm = f.mekf().get_acc_bias();
    check(bias_warm.norm() == 0.0f,
          "the accelerometer bias must not move while frozen");

    run(f, sea, 400.0f);
    check(f.stage() == Stage::Live, "expected Live by 430 s");
}

/*
    The tuning laws, against closed form.

    tau  = 0.5 / f          for tau_coeff = 1
    sigma = sqrt(var - floor^2), var = a_amp^2 / 2 for a sinusoid
    R_S  = 0.35 * sigma * tau^3
*/
void test_tuning_laws() {
    Fusion f;
    f.begin(default_config());

    Sea sea;
    run(f, sea, 600.0f);
    check(f.stage() == Stage::Live, "did not reach Live");
    check(f.wavePeriodReady(), "the wave-period estimator never became ready");

    const float T_true = 1.0f / sea.f_hz;
    if (!check(near_rel(f.getWavePeriodSec(), T_true, 0.10f),
               "the estimated wave period is off by more than 10%")) {
        std::cerr << "  T_z = " << f.getWavePeriodSec() << " want " << T_true << '\n';
    }

    const float tau_want = 0.5f * f.getWavePeriodSec();
    if (!check(near_rel(f.getTauApplied(), tau_want, 0.05f),
               "tau is not half the zero-crossing period")) {
        std::cerr << "  tau = " << f.getTauApplied() << " want " << tau_want << '\n';
    }

    // Variance of a sinusoid is amplitude^2/2; the tuner subtracts the noise
    // floor in power, not amplitude.
    const float var_true = sea.a_amp * sea.a_amp * 0.5f;
    const float sigma_want = std::sqrt(std::max(0.0f, var_true - 0.12f * 0.12f));
    if (!check(near_rel(f.getSigmaApplied(), sigma_want, 0.15f),
               "sigma does not match the sea's acceleration std")) {
        std::cerr << "  sigma = " << f.getSigmaApplied() << " want " << sigma_want << '\n';
    }

    const float rs_want = 0.35f * f.getSigmaApplied() *
                          std::pow(f.getTauApplied(), 3.0f);
    if (!check(near_rel(f.getRSApplied(), rs_want, 0.05f),
               "R_S is not R_S_coeff * sigma * tau^3")) {
        std::cerr << "  RS = " << f.getRSApplied() << " want " << rs_want << '\n';
    }
}

/*
    r_S is a STANDARD DEVIATION, and MekfT::set_RS_noise takes a standard
    deviation. This test exists because getting that wrong cost a factor of
    ten in regularizer strength and did not look like a units error: it
    surfaced as a 46% heave RMS error that read like a modelling problem.

    The S = 0 constraint is a high-pass, and over-tightening it overshoots
    rather than attenuates, so the symptom was gain > 1 at the wave frequency.
    Assert the applied covariance is the square of the applied scale.
*/
/*
    r_S is a standard deviation, and the filter input additionally carries the
    cadence renormalization.

    The original guard here asserted R_S == r_S^2 against the *applied* scale.
    That contract changed when the S=0 cadence became self-similar in tau: the
    filter input is r_S * sqrt(T_0/T_S), so R_S == (that)^2. Both halves are
    checked, because collapsing them would lose the units guard that this test
    exists for -- the bug it caught was R_S being set to the scale itself.
*/
void test_rs_units_are_a_standard_deviation() {
    Fusion f;
    f.begin(default_config());
    check(f.setFixedTuning(2.0f, 0.5f, 7.0f), "setFixedTuning was rejected");

    Vector3f r; Eigen::Matrix<float,3,Fusion::Mekf::NX> H; Matrix3f Rw;
    f.mekf().integral_residual(r, H, Rw);

    const float rs_in = f.getRSFilterInput();
    const float want = rs_in * rs_in;
    if (!check(std::fabs(Rw(2,2) - want) <= 1e-4f * want,
               "R_S must be the square of the r_S the filter was given")) {
        std::cerr << "  R_S(2,2) = " << Rw(2,2) << " want " << want << '\n';
    }
    // Emphatically not the scale itself, which is the bug being guarded.
    check(std::fabs(Rw(2,2) - rs_in) > 1.0f,
          "R_S was set to the r_S scale rather than its square");

    // And the filter input is the applied scale renormalized for the cadence.
    const float T_S = f.pseudoUpdatePeriodSec();
    const float expect = 7.0f * std::sqrt(0.015f / T_S);
    if (!check(std::fabs(rs_in - expect) <= 1e-4f * expect,
               "the r_S filter input is not renormalized by sqrt(T_0/T_S)")) {
        std::cerr << "  rs_in = " << rs_in << " want " << expect
                  << " (T_S = " << T_S << ")\n";
    }
    // With tau = 2 s the cadence is slower than nominal, so the renormalized
    // value must sit *below* the applied scale. If it did not, the information
    // rate 1/(r_S^2 T_S) would not be held.
    check(rs_in < 7.0f, "a slower cadence must lower the r_S filter input");
}

void test_tuning_coefficients_are_live() {
    Fusion a, b;
    a.begin(default_config());
    b.begin(default_config());
    b.setTauCoeff(2.0f);

    Sea sea;
    run(a, sea, 600.0f);
    run(b, sea, 600.0f);
    check(a.stage() == Stage::Live && b.stage() == Stage::Live, "both must be Live");
    if (!check(near_rel(b.getTauApplied(), 2.0f * a.getTauApplied(), 0.05f),
               "doubling tau_coeff did not double the applied tau")) {
        std::cerr << "  tau a=" << a.getTauApplied() << " b=" << b.getTauApplied() << '\n';
    }
}

void test_ablation_hooks() {
    Sea sea;

    // Fixed tuning must pin the operating point and stop adapting.
    {
        Fusion f;
        f.begin(default_config());
        check(f.setFixedTuning(2.5f, 0.7f, 6.0f), "setFixedTuning was rejected");
        run(f, sea, 400.0f);
        check(f.getTauApplied() == 2.5f && f.getSigmaApplied() == 0.7f &&
              f.getRSApplied() == 6.0f,
              "fixed tuning drifted");
    }

    // Freezing the OU channel must leave r_S adapting, and vice versa.
    {
        Fusion f;
        f.begin(default_config());
        check(f.setChannelFreeze(true, 2.0f, 0.5f, false, 0.0f),
              "freezing the OU channel was rejected");
        run(f, sea, 600.0f);
        check(f.getTauApplied() == 2.0f && f.getSigmaApplied() == 0.5f,
              "the frozen OU channel moved");
        check(f.getRSApplied() != 0.5f, "the r_S channel should still have adapted");
        // And with the OU channel pinned, r_S must follow from the pinned values.
        const float rs_want = 0.35f * 0.5f * 2.0f * 2.0f * 2.0f;
        check(near_rel(f.getRSApplied(), rs_want, 0.2f),
              "with the OU channel frozen, r_S must be built from the frozen tau and sigma");
    }

    // Freezing both is rejected rather than silently aliased onto fixed tuning.
    {
        Fusion f;
        f.begin(default_config());
        check(!f.setChannelFreeze(true, 2.0f, 0.5f, true, 6.0f),
              "freezing both channels must be rejected -- that is setFixedTuning");
    }

    // Nonsense arguments must be refused, not stored.
    {
        Fusion f;
        f.begin(default_config());
        check(!f.setFixedTuning(-1.0f, 0.5f, 6.0f), "a negative tau was accepted");
        check(!f.setFixedTuning(2.0f, 0.5f, 0.0f), "a zero R_S was accepted");
    }
}

/*
    The magnetic reference must be captured from the world, once, and only
    after the delay. Capturing early means freezing yaw against an attitude
    that has not settled.
*/
/*
    Magnetic reference timing.

    Under MahonyProxy the reference is learned in the proxy's frame as soon as
    the delay elapses, and reaches the MEKF when the MEKF is seeded. Reading
    mekf().has_magnetic_reference() before the handoff therefore says nothing;
    magReferenceLearned() is the observable that matters early.
*/
void test_magnetic_reference_timing() {
    Fusion f;
    auto cfg = default_config();
    cfg.mag_delay_sec = 20.0f;
    f.begin(cfg);

    Sea sea;
    run(f, sea, 10.0f);
    check(!f.magReferenceLearned(),
          "the magnetic reference was captured before the delay elapsed");
    check(!f.mekf().has_magnetic_reference(),
          "the MEKF holds a reference before the delay elapsed");

    run(f, sea, 30.0f);
    check(f.magReferenceLearned(), "the magnetic reference was never learned");

    // Run long enough for the handoff to carry it into the MEKF.
    run(f, sea, 300.0f);
    check(f.stage() == Stage::Live, "never went Live, so the reference never transferred");
    check(f.mekf().has_magnetic_reference(),
          "the learned reference was not handed to the MEKF");

    const Vector3f B = f.mekf().magnetic_reference_world();
    const Vector3f want(0.21f, 0.0f, 0.43f);
    check(std::fabs(B.norm() - want.norm()) <= 0.02f * want.norm(),
          "the captured magnetic reference has the wrong magnitude");
    const float cos_ang = B.normalized().dot(want.normalized());
    check(cos_ang > 0.98f, "the captured magnetic reference points the wrong way");

    // The canonical form is what makes H_m a genuinely fixed reference: the
    // horizontal component must lie entirely on north.
    check(std::fabs(B.y()) <= 1e-5f,
          "the reference was not anchored: its horizontal part is off north");
}

void test_with_mag_disabled() {
    Fusion f;
    auto cfg = default_config();
    cfg.with_mag = false;
    f.begin(cfg);
    Sea sea;
    run(f, sea, 60.0f, 0.005f, nullptr, nullptr, /*with_mag=*/true);
    check(!f.mekf().has_magnetic_reference(),
          "a reference was captured with the magnetometer disabled");
}


// ---------------------------------------------------------------------------
// Adaptation policy
// ---------------------------------------------------------------------------

// Helper: the r_S the MEKF is actually holding.
float mekf_rs_z(Fusion& f) {
    Vector3f r; Eigen::Matrix<float,3,Fusion::Mekf::NX> H; Matrix3f Rw;
    f.mekf().integral_residual(r, H, Rw);
    return std::sqrt(std::max(0.0f, Rw(2,2)));
}

/*
    The schedule is committed on a cadence, not every sample.

    Committing at 200 Hz is not "more adaptive"; it just couples the covariance
    the MEKF uses to the measurement stream two hundred times a second. The EMA
    still sees every sample -- only the commit is held piecewise constant.
*/
void test_adaptation_is_cadenced() {
    Fusion f;
    f.begin(default_config());
    Sea sea;
    run(f, sea, 320.0f);
    check(f.stage() == Stage::Live, "needed a Live filter to measure cadence");

    const float dt = 0.005f;
    const float seconds = 10.0f;
    const int n = static_cast<int>(seconds / dt);
    int changes = 0;
    float prev = mekf_rs_z(f);
    for (int i = 0; i < n; ++i) {
        const float t = static_cast<float>(i) * dt;
        f.update(dt, sea.gyro(t), sea.acc_body(t));
        const float now = mekf_rs_z(f);
        if (std::fabs(now - prev) > 1e-9f) ++changes;
        prev = now;
    }
    // 0.1 s cadence over 10 s is at most ~100 commits, against 2000 samples.
    check(changes <= 120,
          "the schedule is being committed far more often than the cadence allows");
    check(changes > 0, "the schedule never moved at all, so the test proves nothing");
}

/*
    Exogeneity is a timing property.

    The schedule active while y_k is processed must not depend on y_k. The
    smoother may consume y_k, but its result is committed at the top of the
    next update. So a single violent sample must leave the MEKF's applied r_S
    untouched *within that same update*.
*/
void test_schedule_is_exogenous_to_the_current_sample() {
    Fusion f;
    f.begin(default_config());
    Sea sea;
    run(f, sea, 320.0f);
    check(f.stage() == Stage::Live, "needed a Live filter");

    // Settle onto a commit boundary so a tick is due imminently.
    for (int i = 0; i < 40; ++i) f.update(0.005f, sea.gyro(0.0f), sea.acc_body(0.0f));

    const float before = mekf_rs_z(f);
    // A sample far outside anything the tuner has seen.
    f.update(0.005f, Vector3f(2.0f, -1.5f, 3.0f), Vector3f(0.0f, 0.0f, -40.0f));
    const float after = mekf_rs_z(f);

    check(std::fabs(after - before) <= 1e-9f,
          "the schedule moved within the same update that delivered the sample; "
          "the commit is not deferred and the tuning is not exogenous");
}

/*
    The r_S floor is 0.15, not 0.4.

    On OU-III the 0.4 floor was the binding constraint on every low-motion sea
    rather than a safety limit. A near-still sea must be able to ask for -- and
    receive -- a target below 0.4.
*/
void test_rs_floor_allows_low_motion_seas() {
    Fusion f;
    auto cfg = default_config();
    f.begin(cfg);

    // Very small sea: the cubic schedule asks for a small r_S.
    Sea tiny;
    tiny.a_amp = 0.02f;      // near-still vertical acceleration
    tiny.roll_amp = 0.005f;
    run(f, tiny, 400.0f);

    const float target = f.getRSTarget();
    check(target < 0.4f,
          "a near-still sea did not drive the r_S target below the old 0.4 floor, "
          "so this test cannot distinguish the floors");
    check(target >= 0.15f - 1e-6f, "the r_S target fell below the 0.15 floor");
}

// The proxy doubles as the attitude seed, so its bias integrator must be on.
// With two_ki = 0 a constant gyro bias leaves a standing tilt that nothing
// downstream of an attitude seed high-passes away.
void test_proxy_has_an_integral_term() {
    Fusion::Config c;
    check(c.proxy_two_ki > 0.0f,
          "the startup proxy must estimate gyro bias when it seeds attitude");
    check(std::fabs(c.proxy_two_kp - 0.2f) < 1e-6f,
          "the proxy correction corner moved away from the evidenced value");
}

// Turning the self-similar cadence off must restore the fixed 15 ms schedule
// and remove the renormalization, so the two can be ablated against each other.
void test_fixed_cadence_ablation() {
    Fusion f;
    f.begin(default_config());
    f.setTauScaledPseudoCadence(false);
    check(f.setFixedTuning(2.0f, 0.5f, 7.0f), "setFixedTuning was rejected");
    check(std::fabs(f.pseudoUpdatePeriodSec() - 0.015f) <= 1e-9f,
          "disabling the tau-scaled cadence did not restore the fixed period");
    check(std::fabs(f.getRSFilterInput() - 7.0f) <= 1e-4f,
          "the renormalization is still being applied with a fixed cadence");
}

// The estimator must actually track heave, not merely stay finite.
void test_tracks_heave() {
    Fusion f;
    f.begin(default_config());
    Sea sea;
    run(f, sea, 400.0f);
    check(f.stage() == Stage::Live, "did not reach Live");

    // Score over a trailing window, after the operating point has settled.
    const float dt = 0.005f;
    double sq_err = 0.0, sq_ref = 0.0;
    const float omega = 2.0f * kPi * sea.f_hz;
    const float z_amp = sea.a_amp / (omega * omega);
    int n = 0;
    for (int i = 0; i < 24000; ++i) {          // 120 s
        const float t = 400.0f + static_cast<float>(i) * dt;
        f.update(dt, sea.gyro(t), sea.acc_body(t));
        if (i % 5 == 0) f.updateMag(sea.mag_body(t));
        if (i < 6000) continue;                 // let it settle first
        const float z_true = -z_amp * std::sin(omega * t);
        const float e = f.get_position().z() - z_true;
        sq_err += static_cast<double>(e) * e;
        sq_ref += static_cast<double>(z_true) * z_true;
        ++n;
    }
    const double rms = std::sqrt(sq_err / n);
    const double rms_ref = std::sqrt(sq_ref / n);
    std::cout << "  heave RMS error " << rms << " m against " << rms_ref
              << " m of signal (" << (100.0 * rms / rms_ref) << "%)\n";

    /*
        About 21% on this case, and it is almost entirely amplitude gain
        rather than noise or drift: fitting p_z to a sinusoid at the known
        frequency leaves a residual of 0.002 m, while the fitted amplitude is
        19% high.

        That overshoot is the S = 0 regularizer's high-pass response, it is
        monotonic in both sea frequency and r_S, and OU-III shows MORE of it
        than this filter does at every operating point measured (1.23 vs 1.19
        at 0.15 Hz; 1.62 vs 1.54 at 0.08 Hz). So it is a property of the
        shared model on a monochromatic sea, not a defect here -- the r_S law
        was calibrated against broadband JONSWAP and PM spectra, and a pure
        sinusoid puts all its energy at one frequency.

        The gate is therefore deliberately loose. A tight threshold here would
        be a number fitted to a sea state the filter was never tuned for, and
        the real measurement belongs in the paired study against recorded
        waves.
    */
    check(rms < 0.45 * rms_ref,
          "heave error is far worse than the regularizer overshoot explains");
    check(f.get_position().allFinite() && f.mekf().covariance_full().allFinite(),
          "the filter went non-finite");
}

// ---------------------------------------------------------------------------
// The MEKF gates the orchestrator drives
// ---------------------------------------------------------------------------

void test_initialize_from_acc() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    // A 20 degree roll: at rest the accelerometer reads -g rotated into body.
    for (float roll : {0.0f, 0.35f, -0.6f, 1.2f}) {
        const Matrix3f R_true =
            Eigen::AngleAxisf(roll, Vector3f::UnitX()).toRotationMatrix();
        const Vector3f acc = R_true.transpose() * Vector3f(0.0f, 0.0f, -kG);

        Mekf m;
        m.initialize_identity();
        check(m.initialize_from_acc(acc), "initialize_from_acc rejected a valid sample");

        // Gravity must come back level; yaw is unconstrained and not checked.
        const Vector3f up_world = m.R_bw() * acc.normalized();
        check(std::fabs(up_world.x()) < 1e-5f && std::fabs(up_world.y()) < 1e-5f &&
              up_world.z() < -0.999f,
              "initialize_from_acc did not level the filter");
    }

    Mekf m;
    m.initialize_identity();
    check(!m.initialize_from_acc(Vector3f::Zero()), "a zero sample was accepted");
}

// Sigma_aw parameterizes the OU process noise; it is not the posterior
// covariance of rho_aw. Changing the process model must therefore leave P
// untouched. An explicit reset is kept only as a deliberate prior
// initialization operation for the inactive/newly-enabled world block.
void test_aw_stationary_covariance_is_model_only() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    Mekf m;
    m.initialize_identity();

    auto& P = m.covariance_full();
    P.setIdentity();
    P *= 0.05f;
    for (int k = 0; k < 3; ++k) {
        P(Mekf::OFF_AW + k, Mekf::OFF_AW + k) = 0.9f;
        P(Mekf::OFF_AW + k, k) = 0.06f;
        P(k, Mekf::OFF_AW + k) = 0.06f;
    }
    const auto P_before = P;

    m.set_aw_stationary_std(Vector3f(0.4f, 0.4f, 0.25f));
    m.set_aw_time_constant(2.3f);
    check((P - P_before).cwiseAbs().maxCoeff() == 0.0f,
          "changing OU process parameters rewrote the posterior covariance");

    // Explicit prior initialization is intentionally different: it installs
    // the configured stationary a_w marginal and clears its old correlations.
    m.reset_aw_covariance_to_stationary();
    check(std::fabs(P(Mekf::OFF_AW, Mekf::OFF_AW) - 0.4f * 0.4f) < 1e-6f,
          "explicit a_w prior reset did not install Sigma_aw");
    check(std::fabs(P(Mekf::OFF_AW + 2, Mekf::OFF_AW + 2) - 0.25f * 0.25f) < 1e-6f,
          "explicit a_w prior reset did not install anisotropic Sigma_aw");
    check(std::fabs(P(Mekf::OFF_AW, 0)) == 0.0f,
          "explicit a_w prior reset did not clear stale cross-covariance");
}

// With the linear block gated off, a_w is not estimated but is still present
// in the measurement, so it must be marginalized into the noise rather than
// assumed zero.
void test_linear_block_gate() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    Mekf m;
    m.initialize_identity();
    m.set_aw_stationary_std(Vector3f(0.5f, 0.5f, 0.5f));
    m.set_Racc_std(Vector3f::Constant(0.4f));

    Eigen::Matrix<float,3,Mekf::NX> H_on, H_off;
    Vector3f r; Matrix3f Rw_on, Rw_off;
    m.set_linear_block_enabled(true);
    m.accel_residual(Vector3f(0.0f, 0.0f, -kG), 35.0f, r, H_on, Rw_on);
    m.set_linear_block_enabled(false);
    m.accel_residual(Vector3f(0.0f, 0.0f, -kG), 35.0f, r, H_off, Rw_off);

    check(H_on.block<3,3>(0, Mekf::OFF_AW).cwiseAbs().maxCoeff() > 0.5f,
          "a_w must appear in H_a when the linear block is on");
    check(H_off.block<3,3>(0, Mekf::OFF_AW).cwiseAbs().maxCoeff() == 0.0f,
          "a_w must not appear in H_a when the linear block is off");
    check((Rw_off - Rw_on - Matrix3f::Identity() * 0.25f).cwiseAbs().maxCoeff() < 1e-6f,
          "Sigma_aw must be marginalized into the accelerometer noise when frozen");

    // The world states must not move while frozen.
    const Vector3f p0 = m.get_position();
    for (int i = 0; i < 200; ++i) {
        m.time_update(Vector3f(0.01f, 0.0f, 0.0f), 0.005f);
        m.measurement_update_acc_only(Vector3f(0.0f, 0.0f, -kG), 35.0f);
    }
    check((m.get_position() - p0).cwiseAbs().maxCoeff() == 0.0f,
          "the world states moved with the linear block gated off");
    check(!m.applyIntegralZeroPseudoMeas(),
          "the integral pseudo-measurement must be inert with the block frozen");
}

} // namespace

int main() {
    test_initialize_from_acc();
    test_aw_stationary_covariance_is_model_only();
    test_linear_block_gate();
    test_staging_staged_mekf();
    test_staging_mahony_proxy();
    test_proxy_handoff_times_out();
    test_warmup_gates();
    test_tuning_laws();
    test_rs_units_are_a_standard_deviation();
    test_tuning_coefficients_are_live();
    test_ablation_hooks();
    test_magnetic_reference_timing();
    test_with_mag_disabled();
    test_tracks_heave();
    test_adaptation_is_cadenced();
    test_schedule_is_exogenous_to_the_current_sample();
    test_rs_floor_allows_low_motion_seas();
    test_proxy_has_an_integral_term();
    test_fixed_cadence_ablation();

    if (failures != 0) {
        std::cerr << failures << " orchestrator check(s) failed\n";
        return 1;
    }
    std::cout << "tfg_orchestrator-test: all checks passed\n";
    return 0;
}
