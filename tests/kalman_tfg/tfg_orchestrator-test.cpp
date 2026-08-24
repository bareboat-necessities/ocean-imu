/*
    The TFG sea-state orchestrator: startup, adaptation statistics, magnetic
    acquisition and the MEKF gates they drive.

    The deployed/default r_S law is the spectral-MSE reduction

        r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S),

    while LegacyCubic preserves the pre-PR embedded-friendly schedule.  Both
    laws are asserted below; switching the default must not erase the old
    operating point or turn its test into a tautology.
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
constexpr float kRSCoeffLegacy = 0.28f;
constexpr float kRSMseCoeff = 0.0538f;
constexpr float kRa = 0.0148f * 0.0148f * (1.0f / 200.0f);
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

float pseudo_period(float tau) {
    return std::min(std::max((0.015f / 1.1f) * tau, 0.005f), 0.25f);
}

float spectral_rs(float tau, float sigma_aB) {
    const float qpow = std::pow(2.0f * kRa, 1.0f / 14.0f);
    const float u = sigma_aB * tau * tau * tau * tau;
    return kRSMseCoeff * qpow * std::pow(u, 6.0f / 7.0f)
         / std::sqrt(pseudo_period(tau));
}

struct Sea {
    float f_hz = 0.15f;
    float a_amp = 1.2f;
    float roll_amp = 0.08f;
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

struct TiltedSea : Sea {
    float pitch_amp = 0.06f;
    float pitch_f = 0.19f;

    Matrix3f R_bw(float t) const {
        return (Eigen::AngleAxisf(roll_amp * std::sin(2.0f * kPi * roll_f * t),
                                  Vector3f::UnitX()) *
                Eigen::AngleAxisf(pitch_amp * std::sin(2.0f * kPi * pitch_f * t),
                                  Vector3f::UnitY())).toRotationMatrix();
    }
    Vector3f gyro(float t) const {
        const float h = 1e-4f;
        const Matrix3f R0 = R_bw(t - h), R1 = R_bw(t + h);
        const Matrix3f W = R0.transpose() * (R1 - R0) / (2.0f * h);
        return Vector3f(W(2, 1), W(0, 2), W(1, 0));
    }
    Vector3f acc_body(float t) const {
        return R_bw(t).transpose() * (world_accel(t) - Vector3f(0, 0, kG));
    }
    Vector3f mag_body(float t) const {
        return R_bw(t).transpose() * Vector3f(20.0f, 0.0f, 43.0f);
    }
};

Fusion::Config default_config() {
    Fusion::Config c;
    c.online_tune_warmup_sec = 10.0f;
    c.mag_delay_sec = 5.0f;
    return c;
}

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

void test_staging_staged_mekf() {
    Fusion f;
    auto cfg = default_config();
    cfg.startup_init_policy = Fusion::StartupInitPolicy::StagedMekf;
    f.begin(cfg);
    check(f.stage() == Stage::Cold, "the filter must start Cold");
    check(!f.mekf().linear_block_enabled(), "the linear block must be gated while Cold");

    Sea sea;
    float t_warm = -1.0f, t_live = -1.0f;
    run(f, sea, 300.0f, 0.005f, &t_warm, &t_live);
    check(f.stage() == Stage::Live, "the filter never reached Live");
    check(t_warm >= 0.0f && t_warm < 5.0f, "levelling should complete within a few seconds");
    check(t_live > t_warm && t_live >= 10.0f, "Live staging/warmup order is wrong");
    check(f.mekf().linear_block_enabled(), "the linear block must be enabled once Live");
}

void test_staging_mahony_proxy() {
    Fusion f;
    auto cfg = default_config();
    cfg.startup_init_policy = Fusion::StartupInitPolicy::MahonyProxy;
    f.begin(cfg);
    check(Fusion::Config{}.startup_init_policy == Fusion::StartupInitPolicy::MahonyProxy,
          "MahonyProxy must be the default startup policy");

    Sea sea;
    float t_warm = -1.0f, t_live = -1.0f;
    run(f, sea, 300.0f, 0.005f, &t_warm, &t_live);
    check(f.stage() == Stage::Live, "the proxy policy never reached Live");
    check(t_warm < 0.0f, "proxy startup must skip TunerWarm");
    check(t_live >= 8.0f, "handoff happened before proxy_startup_min_sec");
    check(!f.handoffTimedOut(), "handoff used timeout rather than converged front end");
    check(f.mekf().linear_block_enabled(), "linear block is not enabled after handoff");
    check(f.mekf().acc_bias_updates_enabled(), "accelerometer-bias updates never opened");
}

void test_proxy_handoff_times_out() {
    Fusion f;
    auto cfg = default_config();
    cfg.startup_init_policy = Fusion::StartupInitPolicy::MahonyProxy;
    cfg.proxy_startup_timeout_sec = 30.0f;
    cfg.with_mag = false;
    f.begin(cfg);

    Sea flat;
    run(f, flat, 25.0f, 0.005f, nullptr, nullptr, false);
    check(f.stage() == Stage::Cold, "handed off before timeout with no ready front end");
    run(f, flat, 20.0f, 0.005f, nullptr, nullptr, false);
    check(f.stage() == Stage::Live && f.handoffTimedOut(), "handoff timeout never fired/reported");
}

void test_proxy_gate_requires_the_aligned_branch() {
    Fusion f;
    auto cfg = default_config();
    cfg.startup_init_policy = Fusion::StartupInitPolicy::MahonyProxy;
    cfg.with_mag = false;
    cfg.proxy_startup_timeout_sec = 1.0e6f;
    cfg.online_tune_warmup_sec = 2.0f;
    cfg.proxy_two_kp = 0.01f;
    cfg.proxy_two_ki = 0.0f;
    f.begin(cfg);

    Sea sea;
    constexpr float dt = 0.005f;
    float t = 0.0f;
    auto drive = [&](float seconds, float acc_sign) {
        const int n = static_cast<int>(seconds / dt);
        for (int i = 0; i < n; ++i) {
            f.update(dt, sea.gyro(t), acc_sign * sea.acc_body(t));
            t += dt;
        }
    };

    drive(5.0f, +1.0f);
    check(f.stage() == Stage::Cold, "proxy handed off before minimum startup");
    drive(200.0f, -1.0f);
    check(f.stage() == Stage::Cold, "gravity gate accepted antipodal branch");
    drive(60.0f, +1.0f);
    check(f.stage() == Stage::Live && !f.handoffTimedOut(),
          "aligned branch did not hand off through the quality gate");
}

void test_warmup_gates() {
    Fusion f;
    auto cfg = default_config();
    cfg.freeze_acc_bias_until_live = true;
    cfg.startup_init_policy = Fusion::StartupInitPolicy::StagedMekf;
    f.begin(cfg);
    Sea sea;
    run(f, sea, 30.0f);
    check(f.stage() == Stage::TunerWarm, "expected staged filter to still be warming at 30 s");
    check(f.mekf().get_acc_bias().norm() == 0.0f, "accelerometer bias moved while frozen");
    run(f, sea, 400.0f);
    check(f.stage() == Stage::Live, "staged filter never reached Live");
}

void test_tuning_laws() {
    Fusion f;
    f.begin(default_config());
    check(f.getRSLaw() == Fusion::RSLaw::SpectralMSE,
          "SpectralMSE is not the deployed TFG default");

    Sea sea;
    run(f, sea, 600.0f);
    check(f.stage() == Stage::Live && f.wavePeriodReady(), "front end never became ready");

    const float T_true = 1.0f / sea.f_hz;
    check(near_rel(f.getWavePeriodSec(), T_true, 0.10f), "wave period is off by more than 10%");
    const float tau_want = 0.5f * f.getWavePeriodSec();
    check(near_rel(f.getTauApplied(), tau_want, 0.05f), "tau is not half the zero-crossing period");

    const float var_true = sea.a_amp * sea.a_amp * 0.5f;
    const float sigma_want = std::sqrt(std::max(0.0f, var_true - 0.12f * 0.12f));
    check(near_rel(f.getSigmaApplied(), sigma_want, 0.15f), "sigma does not match wave acceleration std");

    const float rs_want = spectral_rs(f.getTauApplied(), f.getSigmaApplied());
    if (!check(near_rel(f.getRSApplied(), rs_want, 0.05f),
               "SpectralMSE r_S does not match its closed form")) {
        std::cerr << "  RS=" << f.getRSApplied() << " want=" << rs_want << '\n';
    }
}

void test_legacy_cubic_law_is_preserved() {
    Fusion f;
    f.begin(default_config());
    f.setEmbeddedFriendlyLegacyRSLaw(true);
    check(f.getRSLaw() == Fusion::RSLaw::LegacyCubic, "legacy r_S switch did not take effect");

    Sea sea;
    run(f, sea, 600.0f);
    const float want = kRSCoeffLegacy * f.getSigmaApplied() * std::pow(f.getTauApplied(), 3.0f);
    if (!check(near_rel(f.getRSApplied(), want, 0.05f),
               "LegacyCubic no longer reproduces R_S_coeff*sigma*tau^3")) {
        std::cerr << "  legacy RS=" << f.getRSApplied() << " want=" << want << '\n';
    }
}

void test_rs_units_are_a_standard_deviation() {
    // SpectralMSE targets already contain 1/sqrt(T_S), so a fixed r_S value is
    // passed through unchanged rather than normalized a second time.
    Fusion f;
    f.begin(default_config());
    check(f.setFixedTuning(2.0f, 0.5f, 7.0f), "setFixedTuning was rejected");
    Vector3f r; Eigen::Matrix<float,3,Fusion::Mekf::NX> H; Matrix3f Rw;
    f.mekf().integral_residual(r, H, Rw);
    check(std::fabs(f.getRSFilterInput() - 7.0f) <= 7e-4f,
          "SpectralMSE fixed r_S was cadence-normalized twice");
    check(std::fabs(Rw(2,2) - 49.0f) <= 4.9e-3f,
          "R_S is not the square of the SpectralMSE standard deviation");

    // LegacyCubic deliberately preserves the historical post-target cadence
    // normalization, including its exact fixed-tuning behaviour.
    Fusion legacy;
    legacy.begin(default_config());
    legacy.setEmbeddedFriendlyLegacyRSLaw(true);
    check(legacy.setFixedTuning(2.0f, 0.5f, 7.0f), "legacy fixed tuning was rejected");
    legacy.mekf().integral_residual(r, H, Rw);
    const float T_S = legacy.pseudoUpdatePeriodSec();
    const float rs_in = 7.0f * std::sqrt(0.015f / T_S);
    check(near_rel(legacy.getRSFilterInput(), rs_in, 1e-4f),
          "LegacyCubic lost sqrt(T0/T_S) cadence normalization");
    check(near_rel(Rw(2,2), rs_in * rs_in, 1e-4f),
          "legacy R_S is not the square of the normalized standard deviation");
    check(legacy.getRSFilterInput() < 7.0f, "slower legacy cadence did not lower filter input");
}

void test_tuning_coefficients_are_live() {
    Fusion a, b;
    a.begin(default_config());
    b.begin(default_config());
    b.setTauCoeff(2.0f);
    Sea sea;
    run(a, sea, 600.0f);
    run(b, sea, 600.0f);
    check(a.stage() == Stage::Live && b.stage() == Stage::Live, "both filters must be Live");
    check(near_rel(b.getTauApplied(), 2.0f * a.getTauApplied(), 0.05f),
          "doubling tau_coeff did not double applied tau");
}

void test_ablation_hooks() {
    Sea sea;
    {
        Fusion f;
        f.begin(default_config());
        check(f.setFixedTuning(2.5f, 0.7f, 6.0f), "setFixedTuning was rejected");
        run(f, sea, 400.0f);
        check(f.getTauApplied() == 2.5f && f.getSigmaApplied() == 0.7f &&
              f.getRSApplied() == 6.0f, "fixed tuning drifted");
    }
    {
        Fusion f;
        f.begin(default_config());
        check(f.setChannelFreeze(true, 2.0f, 0.5f, false, 0.0f),
              "freezing OU channel was rejected");
        run(f, sea, 600.0f);
        check(f.getTauApplied() == 2.0f && f.getSigmaApplied() == 0.5f,
              "frozen OU channel moved");
        check(std::isfinite(f.getRSApplied()) && f.getRSApplied() > 0.15f,
              "r_S did not continue adapting from the live front-end estimate");
    }
    {
        Fusion f;
        f.begin(default_config());
        check(!f.setChannelFreeze(true, 2.0f, 0.5f, true, 6.0f),
              "freezing both channels must be rejected");
        check(!f.setFixedTuning(-1.0f, 0.5f, 6.0f), "negative tau was accepted");
        check(!f.setFixedTuning(2.0f, 0.5f, 0.0f), "zero R_S was accepted");
    }
}

void test_magnetic_reference_timing() {
    Fusion f;
    auto cfg = default_config();
    cfg.startup_init_policy = Fusion::StartupInitPolicy::MahonyProxy;
    cfg.mag_delay_sec = 20.0f;
    cfg.proxy_gravity_hold_sec = 5.0f;
    f.begin(cfg);
    Sea sea;
    run(f, sea, 10.0f);
    check(!f.magReferenceLearned() && !f.mekf().has_magnetic_reference(),
          "mag reference was captured before delay");
    run(f, sea, 30.0f);
    check(f.magReferenceLearned(), "magnetic reference was never learned");
    run(f, sea, 300.0f);
    check(f.stage() == Stage::Live && f.mekf().has_magnetic_reference(),
          "learned reference was not transferred at handoff");
    const Vector3f B = f.mekf().magnetic_reference_world();
    const Vector3f want(0.21f, 0.0f, 0.43f);
    check(std::fabs(B.norm() - want.norm()) <= 0.02f * want.norm(),
          "magnetic reference magnitude is wrong");
    check(B.normalized().dot(want.normalized()) > 0.98f, "magnetic reference points wrong way");
    check(std::fabs(B.y()) <= 1e-5f, "reference horizontal component is not anchored north");
}

void test_with_mag_disabled() {
    Fusion f;
    auto cfg = default_config();
    cfg.with_mag = false;
    f.begin(cfg);
    Sea sea;
    run(f, sea, 60.0f, 0.005f, nullptr, nullptr, true);
    check(!f.mekf().has_magnetic_reference(), "reference captured with magnetometer disabled");
}

float mekf_rs_z(Fusion& f) {
    Vector3f r; Eigen::Matrix<float,3,Fusion::Mekf::NX> H; Matrix3f Rw;
    f.mekf().integral_residual(r, H, Rw);
    return std::sqrt(std::max(0.0f, Rw(2,2)));
}

void test_adaptation_is_cadenced() {
    Fusion f;
    f.begin(default_config());
    Sea sea;
    run(f, sea, 320.0f);
    check(f.stage() == Stage::Live, "needed Live filter to measure adaptation cadence");
    const float dt = 0.005f;
    int changes = 0;
    float prev = mekf_rs_z(f);
    for (int i = 0; i < 2000; ++i) {
        f.update(dt, sea.gyro(0.0f), sea.acc_body(0.0f));
        const float now = mekf_rs_z(f);
        if (std::fabs(now - prev) > 1e-9f) ++changes;
        prev = now;
    }
    check(changes <= 120 && changes > 0, "schedule commit cadence is wrong");
}

void test_sea_scaled_ema_defaults() {
    Fusion f;
    f.begin(default_config());
    check(std::fabs(f.getAdaptationSeaPeriods() - 0.40f) <= 1e-6f,
          "tau/sigma EMA is not 0.40*T_sea");
    check(std::fabs(f.getRSAdaptMult() - 1.5f) <= 1e-6f,
          "r_S EMA is not 1.5*tau_target as in OU-III");
    check(std::fabs(f.getSigmaVarianceHorizonPeriods() - 4.0f) <= 1e-6f,
          "sigma variance horizon is not 4 sea periods");
    check(f.getTunerFreqSmoothingSeaPeriods() == 0.0f,
          "removed second frequency EMA was re-enabled");

    WavePeriodEstimator period;
    check(std::fabs(period.getMomentHorizonPeriods() - 4.0f) <= 1e-6f,
          "wave-period moment horizon is not 4*T_z");
    check(std::fabs(period.getLogSmoothingPeriods() - 0.05f) <= 1e-6f,
          "canonical log-period smoothing is not 0.05*T_z");

    const Fusion::Config c;
    check(std::fabs(c.online_tune_warmup_sec - 10.0f) <= 1e-6f,
          "TFG tuner warmup differs from deployed OU-III wrapper");
    check(std::fabs(c.proxy_two_kp - 0.2f) <= 1e-6f &&
          std::fabs(c.proxy_two_ki - 0.02f) <= 1e-6f,
          "startup Mahony gains differ from OU-III");
    check(std::fabs(c.proxy_gravity_lpf_sec - 12.0f) <= 1e-6f &&
          std::fabs(c.proxy_gravity_warmup_sec - 5.0f) <= 1e-6f,
          "world-frame gravity trust averaging differs from OU-III");
    check(c.mag_min_samples == 128 && std::fabs(c.mag_min_window_sec - 15.0f) <= 1e-6f,
          "mag initial acquisition window differs from OU-III");
    check(std::fabs(c.mag_refine_start_sec - 90.0f) <= 1e-6f &&
          std::fabs(c.mag_refine_window_sec - 30.0f) <= 1e-6f,
          "mag refinement schedule differs from OU-III");
    check(std::fabs(c.mag_hi_memory_sec - 600.0f) <= 1e-6f &&
          std::fabs(c.mag_hi_model_ridge - 5e-4f) <= 1e-8f &&
          std::fabs(c.mag_hi_model_ridge_relative - 0.5f) <= 1e-6f &&
          std::fabs(c.mag_hi_slew_tau_sec - 45.0f) <= 1e-6f,
          "continuous hard-iron coefficients differ from OU-III");

    f.setAdaptationTimeConstants(1.8f);
    check(f.getAdaptationSeaPeriods() == 0.0f,
          "fixed tau/sigma EMA setter did not disable sea scaling");
    f.setTunerFreqSmoothingTimeConstant(1.0f);
    check(f.getTunerFreqSmoothingSeaPeriods() == 0.0f,
          "legacy frequency setter resurrected a second smoother");
}

void test_schedule_is_exogenous_to_the_current_sample() {
    Fusion f;
    f.begin(default_config());
    Sea sea;
    run(f, sea, 320.0f);
    check(f.stage() == Stage::Live, "needed Live filter");
    for (int i = 0; i < 40; ++i) f.update(0.005f, sea.gyro(0.0f), sea.acc_body(0.0f));
    const float before = mekf_rs_z(f);
    f.update(0.005f, Vector3f(2.0f, -1.5f, 3.0f), Vector3f(0.0f, 0.0f, -40.0f));
    const float after = mekf_rs_z(f);
    check(std::fabs(after - before) <= 1e-9f,
          "schedule moved within same update that delivered current measurement");
}

void test_rs_floor_allows_low_motion_seas() {
    Fusion f;
    f.begin(default_config());
    Sea tiny;
    tiny.a_amp = 0.02f;
    tiny.roll_amp = 0.005f;
    run(f, tiny, 400.0f);
    const float target = f.getRSTarget();
    check(target < 0.4f, "near-still sea did not reach below old 0.4 floor");
    check(target >= 0.15f - 1e-6f, "r_S target fell below 0.15 floor");
}

void test_proxy_has_an_integral_term() {
    Fusion::Config c;
    check(c.proxy_two_ki > 0.0f, "startup proxy must estimate gyro bias");
    check(std::fabs(c.proxy_two_kp - 0.2f) < 1e-6f, "proxy correction corner moved");
}

void test_fixed_cadence_ablation() {
    Fusion f;
    f.begin(default_config());
    f.setTauScaledPseudoCadence(false);
    check(f.setFixedTuning(2.0f, 0.5f, 7.0f), "setFixedTuning was rejected");
    check(std::fabs(f.pseudoUpdatePeriodSec() - 0.015f) <= 1e-9f,
          "fixed-cadence ablation did not restore 15 ms");
    check(std::fabs(f.getRSFilterInput() - 7.0f) <= 1e-4f,
          "fixed SpectralMSE r_S changed under fixed cadence");
}

void test_tracks_heave() {
    Fusion f;
    f.begin(default_config());
    Sea sea;
    run(f, sea, 400.0f);
    check(f.stage() == Stage::Live, "did not reach Live");
    const float dt = 0.005f;
    double sq_err = 0.0, sq_ref = 0.0;
    const float omega = 2.0f * kPi * sea.f_hz;
    const float z_amp = sea.a_amp / (omega * omega);
    int n = 0;
    for (int i = 0; i < 24000; ++i) {
        const float t = 400.0f + static_cast<float>(i) * dt;
        f.update(dt, sea.gyro(t), sea.acc_body(t));
        if (i % 5 == 0) f.updateMag(sea.mag_body(t));
        if (i < 6000) continue;
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
    check(rms < 0.45 * rms_ref, "heave error is far worse than regularizer overshoot explains");
    check(f.get_position().allFinite() && f.mekf().covariance_full().allFinite(),
          "filter went non-finite");
}

void test_mag_reference_is_windowed() {
    Fusion f;
    auto cfg = default_config();
    cfg.mag_delay_sec = 5.0f;
    cfg.mag_min_window_sec = 20.0f;
    f.begin(cfg);
    Sea sea;
    run(f, sea, 60.0f);
    check(f.magReferenceLearned(), "reference was never learned");
    check(f.magNorthLockTimeSec() >= cfg.mag_delay_sec + cfg.mag_min_window_sec,
          "reference locked before averaging window could close");
}

void test_mag_reference_is_refined() {
    Fusion f;
    auto cfg = default_config();
    cfg.mag_refine_start_sec = 150.0f;
    cfg.mag_refine_window_sec = 15.0f;
    cfg.acc_bias_unlock_sec = 5.0f;
    f.begin(cfg);
    Sea sea;
    run(f, sea, 140.0f);
    check(f.isLive(), "filter must be live before refinement opens");
    check(!f.magReferenceRefined(), "refinement ran before start time");
    check(!f.mekf().acc_bias_updates_enabled(), "bias opened on provisional magnetic reference");
    run(f, sea, 120.0f);
    check(f.magReferenceRefined(), "refinement never completed");
    check(f.magRefineTimeSec() >= cfg.mag_refine_start_sec, "refinement completed before start time");
    check(f.mekf().acc_bias_updates_enabled(), "bias never opened after refinement");
}

void test_continuous_hard_iron_recovers_heading() {
    const Vector3f offset(0.0f, 3.0f, 0.0f);
    auto yaw_after = [&](bool hard_iron_on, const Vector3f& b) {
        Fusion f;
        auto cfg = default_config();
        cfg.mag_continuous_hard_iron = hard_iron_on;
        cfg.mag_refine_start_sec = 150.0f;
        cfg.mag_refine_window_sec = 15.0f;
        f.begin(cfg);
        TiltedSea sea;
        const float dt = 0.005f;
        const int n = static_cast<int>(400.0f / dt);
        for (int i = 0; i < n; ++i) {
            const float t = static_cast<float>(i) * dt;
            f.update(dt, sea.gyro(t), sea.acc_body(t));
            if (i % 5 == 0) f.updateMag(sea.mag_body(t) + b);
        }
        const Matrix3f R = f.quaternion().toRotationMatrix();
        return std::atan2(R(1, 0), R(0, 0)) * 57.29577951308232f;
    };
    const float yaw_clean = yaw_after(true, Vector3f::Zero());
    const float yaw_off = yaw_after(false, offset);
    const float yaw_on = yaw_after(true, offset);
    const float err_off = std::fabs(yaw_off - yaw_clean);
    const float err_on = std::fabs(yaw_on - yaw_clean);
    std::cout << "  hard-iron yaw error: " << err_off << " deg uncorrected, "
              << err_on << " deg corrected\n";
    check(err_off > 5.0f, "injected hard iron did not produce a heading error");
    check(err_on < 0.75f * err_off, "continuous hard iron did not recover enough heading error");
}

void test_hard_iron_can_be_disabled() {
    Fusion f;
    auto cfg = default_config();
    cfg.mag_continuous_hard_iron = false;
    f.begin(cfg);
    TiltedSea sea;
    const float dt = 0.005f;
    for (int i = 0; i < static_cast<int>(400.0f / dt); ++i) {
        const float t = static_cast<float>(i) * dt;
        f.update(dt, sea.gyro(t), sea.acc_body(t));
        if (i % 5 == 0) f.updateMag(sea.mag_body(t) + Vector3f(0.0f, 3.0f, 0.0f));
    }
    check(f.magHardIronBodyUT().norm() == 0.0f, "hard iron applied while estimator disabled");
}

void test_aw_cov_sync_only_inflates() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    constexpr int AW = Mekf::OFF_AW;
    auto build = [](float start_var, float stationary_std) {
        Mekf m;
        m.initialize_identity();
        m.set_linear_block_enabled(true);
        m.set_aw_stationary_std(Vector3f::Constant(stationary_std));
        m.covariance_full().block<3,3>(AW, AW) = Matrix3f::Identity() * start_var;
        return m;
    };
    {
        Mekf synced = build(1e-6f, 0.5f), plain = build(1e-6f, 0.5f);
        synced.synchronize_aw_covariance_to_stationary();
        synced.time_update(Vector3f::Zero(), 0.005f);
        plain.time_update(Vector3f::Zero(), 0.005f);
        const float got = synced.covariance_full()(AW+2, AW+2);
        const float ref = plain.covariance_full()(AW+2, AW+2);
        check(got > 10.0f * ref && got >= 0.9f * 0.25f,
              "a_w sync did not lift shortfall toward stationary variance");
    }
    {
        Mekf synced = build(4.0f, 0.3f), plain = build(4.0f, 0.3f);
        synced.synchronize_aw_covariance_to_stationary();
        synced.time_update(Vector3f::Zero(), 0.005f);
        plain.time_update(Vector3f::Zero(), 0.005f);
        for (int i = 0; i < 3; ++i) {
            check(synced.covariance_full()(AW+i,AW+i) >=
                  plain.covariance_full()(AW+i,AW+i) - 1e-6f,
                  "a_w sync reduced a marginal it may only inflate");
        }
    }
}

void test_initialize_from_acc() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    for (float roll : {0.0f, 0.35f, -0.6f, 1.2f}) {
        const Matrix3f R_true = Eigen::AngleAxisf(roll, Vector3f::UnitX()).toRotationMatrix();
        const Vector3f acc = R_true.transpose() * Vector3f(0.0f, 0.0f, -kG);
        Mekf m;
        m.initialize_identity();
        check(m.initialize_from_acc(acc), "initialize_from_acc rejected valid sample");
        const Vector3f up_world = m.R_bw() * acc.normalized();
        check(std::fabs(up_world.x()) < 1e-5f && std::fabs(up_world.y()) < 1e-5f &&
              up_world.z() < -0.999f, "initialize_from_acc did not level filter");
    }
    Mekf m;
    m.initialize_identity();
    check(!m.initialize_from_acc(Vector3f::Zero()), "zero accelerometer sample was accepted");
}

void test_aw_stationary_covariance_is_model_only() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    Mekf m;
    m.initialize_identity();
    auto& P = m.covariance_full();
    P.setIdentity(); P *= 0.05f;
    for (int k = 0; k < 3; ++k) {
        P(Mekf::OFF_AW+k, Mekf::OFF_AW+k) = 0.9f;
        P(Mekf::OFF_AW+k, k) = 0.06f;
        P(k, Mekf::OFF_AW+k) = 0.06f;
    }
    const auto before = P;
    m.set_aw_stationary_std(Vector3f(0.4f,0.4f,0.25f));
    m.set_aw_time_constant(2.3f);
    check((P-before).cwiseAbs().maxCoeff() == 0.0f,
          "changing OU model parameters rewrote posterior covariance");
    m.reset_aw_covariance_to_stationary();
    check(std::fabs(P(Mekf::OFF_AW,Mekf::OFF_AW)-0.16f) < 1e-6f,
          "explicit a_w reset did not install stationary variance");
    check(std::fabs(P(Mekf::OFF_AW+2,Mekf::OFF_AW+2)-0.0625f) < 1e-6f,
          "explicit anisotropic a_w variance is wrong");
    check(std::fabs(P(Mekf::OFF_AW,0)) == 0.0f,
          "explicit a_w reset did not clear stale cross covariance");
}

void test_linear_block_gate() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    Mekf m;
    m.initialize_identity();
    m.set_aw_stationary_std(Vector3f::Constant(0.5f));
    m.set_Racc_std(Vector3f::Constant(0.4f));
    Eigen::Matrix<float,3,Mekf::NX> H_on, H_off;
    Vector3f r; Matrix3f Rw_on, Rw_off;
    m.set_linear_block_enabled(true);
    m.accel_residual(Vector3f(0,0,-kG),35.0f,r,H_on,Rw_on);
    m.set_linear_block_enabled(false);
    m.accel_residual(Vector3f(0,0,-kG),35.0f,r,H_off,Rw_off);
    check(H_on.block<3,3>(0,Mekf::OFF_AW).cwiseAbs().maxCoeff() > 0.5f,
          "a_w missing from H_a with linear block enabled");
    check(H_off.block<3,3>(0,Mekf::OFF_AW).cwiseAbs().maxCoeff() == 0.0f,
          "a_w remained in H_a with linear block disabled");
    check((Rw_off-Rw_on-Matrix3f::Identity()*0.25f).cwiseAbs().maxCoeff() < 1e-6f,
          "Sigma_aw was not marginalized into accelerometer noise");
    const Vector3f p0 = m.get_position();
    for (int i=0;i<200;++i) {
        m.time_update(Vector3f(0.01f,0,0),0.005f);
        m.measurement_update_acc_only(Vector3f(0,0,-kG),35.0f);
    }
    check((m.get_position()-p0).cwiseAbs().maxCoeff() == 0.0f,
          "world states moved with linear block gated off");
    check(!m.applyIntegralZeroPseudoMeas(), "integral pseudo update active while block frozen");
}

} // namespace

int main() {
    test_initialize_from_acc();
    test_aw_stationary_covariance_is_model_only();
    test_linear_block_gate();
    test_staging_staged_mekf();
    test_staging_mahony_proxy();
    test_proxy_handoff_times_out();
    test_proxy_gate_requires_the_aligned_branch();
    test_warmup_gates();
    test_tuning_laws();
    test_legacy_cubic_law_is_preserved();
    test_rs_units_are_a_standard_deviation();
    test_tuning_coefficients_are_live();
    test_ablation_hooks();
    test_magnetic_reference_timing();
    test_with_mag_disabled();
    test_tracks_heave();
    test_adaptation_is_cadenced();
    test_sea_scaled_ema_defaults();
    test_schedule_is_exogenous_to_the_current_sample();
    test_rs_floor_allows_low_motion_seas();
    test_proxy_has_an_integral_term();
    test_fixed_cadence_ablation();
    test_mag_reference_is_windowed();
    test_mag_reference_is_refined();
    test_continuous_hard_iron_recovers_heading();
    test_hard_iron_can_be_disabled();
    test_aw_cov_sync_only_inflates();

    if (failures != 0) {
        std::cerr << failures << " orchestrator check(s) failed\n";
        return 1;
    }
    std::cout << "tfg_orchestrator-test: all checks passed\n";
    return 0;
}
