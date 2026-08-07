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

void test_staging() {
    Fusion f;
    f.begin(default_config());
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

    // The Live gate waits on the wave-period estimator as well as the warmup
    // clock, so it is necessarily later than the clock alone.
    check(t_live >= 10.0f,
          "Live was reached before the configured warmup elapsed");
}

// While warming, the accelerometer bias must stay frozen and Racc inflated:
// a filter that has not yet levelled cannot tell bias from tilt.
void test_warmup_gates() {
    Fusion f;
    auto cfg = default_config();
    cfg.freeze_acc_bias_until_live = true;
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

// The coefficients must actually be the knobs the study will turn.
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
void test_magnetic_reference_timing() {
    Fusion f;
    auto cfg = default_config();
    cfg.mag_delay_sec = 20.0f;
    f.begin(cfg);

    Sea sea;
    run(f, sea, 10.0f);
    check(!f.mekf().has_magnetic_reference(),
          "the magnetic reference was captured before the delay elapsed");

    run(f, sea, 60.0f);
    check(f.mekf().has_magnetic_reference(),
          "the magnetic reference was never captured");

    // It should point roughly where the true world field does.
    const Vector3f B = f.mekf().magnetic_reference_world();
    const Vector3f want(0.21f, 0.0f, 0.43f);
    check(std::fabs(B.norm() - want.norm()) <= 0.02f * want.norm(),
          "the captured magnetic reference has the wrong magnitude");
    const float cos_ang = B.normalized().dot(want.normalized());
    check(cos_ang > 0.98f, "the captured magnetic reference points the wrong way");
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
    check(rms < rms_ref,
          "heave error exceeds the signal itself -- the filter is not tracking");
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

// The congruent sync must replace the a_w marginal while preserving every
// correlation coefficient -- that is the whole point of it over a reset.
void test_aw_covariance_sync_preserves_correlations() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    Mekf m;
    m.initialize_identity();
    m.set_aw_stationary_std(Vector3f(0.4f, 0.4f, 0.25f));

    // Seed a covariance with real cross-correlation between a_w and phi.
    auto& P = m.covariance_full();
    P.setIdentity();
    P *= 0.05f;
    for (int k = 0; k < 3; ++k) {
        P(Mekf::OFF_AW + k, Mekf::OFF_AW + k) = 0.9f;
        P(Mekf::OFF_AW + k, k) = 0.06f;
        P(k, Mekf::OFF_AW + k) = 0.06f;
    }

    auto corr = [&](int i, int j) {
        return P(i, j) / std::sqrt(P(i, i) * P(j, j));
    };
    const float c_before = corr(Mekf::OFF_AW, 0);

    check(m.synchronize_aw_covariance_to_stationary_congruent(),
          "the congruent sync failed");

    const float c_after = corr(Mekf::OFF_AW, 0);
    check(std::fabs(c_after - c_before) < 1e-4f,
          "the congruent sync did not preserve the correlation coefficient");
    check(std::fabs(P(Mekf::OFF_AW, Mekf::OFF_AW) - 0.4f * 0.4f) < 1e-5f,
          "the congruent sync did not install the stationary marginal");

    // A plain reset must instead zero the cross terms.
    m.reset_aw_covariance_to_stationary();
    check(std::fabs(P(Mekf::OFF_AW, 0)) == 0.0f,
          "reset must clear the a_w cross-covariances");
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
    test_aw_covariance_sync_preserves_correlations();
    test_linear_block_gate();
    test_staging();
    test_warmup_gates();
    test_tuning_laws();
    test_tuning_coefficients_are_live();
    test_ablation_hooks();
    test_magnetic_reference_timing();
    test_with_mag_disabled();
    test_tracks_heave();

    if (failures != 0) {
        std::cerr << failures << " orchestrator check(s) failed\n";
        return 1;
    }
    std::cout << "tfg_orchestrator-test: all checks passed\n";
    return 0;
}
