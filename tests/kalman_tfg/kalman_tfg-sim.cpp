/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Simulator for the right-invariant two-frame Lie-group filter, over the same
    recorded wave records the OU-III simulator uses, through the same
    IW3dFusionAdapter interface and the same scoring path. Paired comparison is
    only meaningful if both filters see identical wave realization, identical
    sensor-error realization and identical scoring, which is what reusing
    W3dSimulationRunner buys.

    WHAT THIS FILTER DOES NOT REPORT. SeaStateFusionFilter_TFG carries the core
    of the OU-III orchestrator and not its wave-direction estimator or its
    alternative frequency trackers, so the direction telemetry and freq_hz stay
    at their NaN defaults rather than being filled with a plausible-looking
    substitute. A study comparing the two families has to compare the channels
    both actually implement; inventing values here would make that mistake
    invisible.

    THE FAILURE LIMITS BELOW ARE REGRESSION BARS, NOT QUALITY TARGETS. They are
    the worst values this filter actually produced over the eight JONSWAP and
    PM-Stokes records, plus a small margin -- derived the same way OU-III's
    were, and deliberately not copied from OU-III, because copying would
    silently assert that a filter with different error characteristics must sit
    inside another one's envelope.

    Two of them record behaviour that is plainly poor, and they are set where
    they are so that further regression is caught, not because the current
    value is acceptable. See the table below.
*/

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/W3dSimCommon.h"
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

bool add_noise = true;

namespace {

bool env_float(const char* name, float& out) {
    if (const char* s = std::getenv(name)) {
        out = static_cast<float>(std::atof(s));
        return true;
    }
    return false;
}

/*
    Tuning arms, matching the OU-III simulator's so the two can be paired arm
    for arm in a study.

      Adaptive        -- the tuner runs throughout
      Fixed           -- pinned operating point, no adaptation
      AdaptiveRSOnly  -- OU channel frozen, r_S keeps adapting
      AdaptiveOUOnly  -- r_S frozen, OU channel keeps adapting
*/
enum class TuningMode { Adaptive, Fixed, AdaptiveRSOnly, AdaptiveOUOnly };

TuningMode load_tuning_mode() {
    const char* raw = std::getenv("TFG_TUNING");
    if (!raw) return TuningMode::Adaptive;
    const std::string mode(raw);
    if (mode == "adaptive") return TuningMode::Adaptive;
    if (mode == "fixed") return TuningMode::Fixed;
    if (mode == "adaptive_rs_only") return TuningMode::AdaptiveRSOnly;
    if (mode == "adaptive_ou_only") return TuningMode::AdaptiveOUOnly;
    throw std::runtime_error("unknown TFG_TUNING mode: " + mode);
}

}  // namespace

class FusionAdapter_TFG final : public IW3dFusionAdapter {
public:
    FusionAdapter_TFG(bool with_mag,
                      const Vector3f& sigma_a_init,
                      const Vector3f& sigma_g,
                      const Vector3f& sigma_m)
    {
        Fusion::Config cfg;
        cfg.with_mag = with_mag;
        cfg.sigma_a = sigma_a_init;
        cfg.gyro_noise_density = sigma_g.x();
        cfg.sigma_m = sigma_m;
        cfg.mag_delay_sec = ocean_imu::tfg::MAG_DELAY_SEC;
        if (float v = 0.0f; env_float("SF_MAG_DELAY_SEC", v)) cfg.mag_delay_sec = v;
        cfg.freeze_acc_bias_until_live = true;
        cfg.Racc_warmup_std = 0.5f;

        tuning_ = load_tuning_mode();
        load_fixed_tuning_();

        fusion_.begin(cfg);

        // Coefficient overrides, so a sweep can move the tuning laws without
        // recompiling. Names mirror the OU-III simulator's OU_III_* set.
        float v = 0.0f;
        if (env_float("TFG_TAU_COEFF", v))      fusion_.setTauCoeff(v);
        if (env_float("TFG_SIGMA_COEFF", v))    fusion_.setSigmaCoeff(v);
        if (env_float("TFG_R_S_COEFF", v))      fusion_.setRSCoeff(v);
        if (env_float("TFG_S_FACTOR", v))       fusion_.setSFactor(v);
        if (env_float("TFG_R_S_XY_FACTOR", v))  fusion_.setRSXYFactor(v);
        if (env_float("TFG_ADAPT_TAU_SEC", v))  fusion_.setAdaptationTimeConstants(v);
        if (env_float("TFG_ADAPT_RS_MULT", v))  fusion_.setRSAdaptMult(v);
        if (env_float("TFG_ADAPT_RS_SLEW_LOG", v)) fusion_.setRSAdaptSlewLog(v);
        if (env_float("TFG_ACC_NOISE_FLOOR", v))   fusion_.setAccNoiseFloorSigma(v);

    }

    void updateMag(const Vector3f& mag_body_ned) override {
        fusion_.updateMag(mag_body_ned);
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);

        // Fixed and frozen arms are applied once the tuner is Live, so the
        // filter reaches its operating point the same way the adaptive arm
        // does and only then stops moving. Applying them from cold would
        // change the startup transient as well as the steady state, and the
        // study could not tell the two apart.
        if (tuning_ == TuningMode::Adaptive || fixed_tuning_applied_) return;
        if (!fusion_.isLive()) return;

        const bool ok = (tuning_ == TuningMode::Fixed)
            ? fusion_.setFixedTuning(fixed_tau_s_, fixed_sigma_a_, fixed_RS_)
            : fusion_.setChannelFreeze(
                  tuning_ == TuningMode::AdaptiveRSOnly,
                  fixed_tau_s_, fixed_sigma_a_,
                  tuning_ == TuningMode::AdaptiveOUOnly,
                  fixed_RS_);
        if (!ok) throw std::runtime_error("invalid TFG tuning point");
        fixed_tuning_applied_ = true;
    }

    FilterSnapshot snapshot() const override {
        FilterSnapshot s;

        s.disp_est_zu = ned_to_zu(fusion_.get_position());
        s.vel_est_zu  = ned_to_zu(fusion_.get_velocity());
        s.acc_est_zu  = ned_to_zu(fusion_.get_world_accel());

        // BODY -> WORLD in NED. After magnetic lock this is the learned
        // magnetic-NED frame, not true north: no declination model is applied,
        // because a bare IMU does not have one.
        const Quaternionf q_bw_ned = fusion_.quaternion().normalized();
        float roll_deg = 0.0f, pitch_deg = 0.0f, yaw_deg = 0.0f;
        quat_to_euler_nautical(q_bw_ned, roll_deg, pitch_deg, yaw_deg);
        s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, wrapDeg(yaw_deg));

        s.acc_bias_est_ned  = fusion_.mekf().get_acc_bias();
        s.gyro_bias_est_ned = fusion_.mekf().gyroscope_bias();
        // No magnetometer-bias state in this filter; left at zero rather than
        // filled with something that would read as an estimate.

        s.tau_target     = fusion_.getTauTarget();
        s.sigma_target   = fusion_.getSigmaTarget();
        s.tuning_target  = fusion_.getRSTarget();
        s.tau_applied    = fusion_.getTauApplied();
        s.sigma_applied  = fusion_.getSigmaApplied();
        s.tuning_applied = fusion_.getRSApplied();

        s.wave_period_sec = fusion_.getWavePeriodSec();
        s.accel_variance  = fusion_.getAccelVariance();

        // freq_hz, period_sec, the displacement/velocity scales and every
        // direction field stay NaN. This filter does not compute them, and a
        // substitute would hide that from any study reading the output.
        return s;
    }

private:
    using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;

    void load_fixed_tuning_() {
        if (tuning_ == TuningMode::Adaptive) return;
        const bool need_ou = (tuning_ != TuningMode::AdaptiveOUOnly);
        const bool need_rs = (tuning_ != TuningMode::AdaptiveRSOnly);

        if (need_ou) {
            if (!env_float("TFG_FIXED_TAU_S", fixed_tau_s_) ||
                !env_float("TFG_FIXED_SIGMA_A", fixed_sigma_a_)) {
                throw std::runtime_error(
                    "TFG_TUNING requires TFG_FIXED_TAU_S and TFG_FIXED_SIGMA_A");
            }
        }
        if (need_rs) {
            if (!env_float("TFG_FIXED_RS", fixed_RS_)) {
                throw std::runtime_error("TFG_TUNING requires TFG_FIXED_RS");
            }
        }
    }

    Fusion fusion_{};
    TuningMode tuning_ = TuningMode::Adaptive;
    bool fixed_tuning_applied_ = false;
    float fixed_tau_s_ = 0.0f;
    float fixed_sigma_a_ = 0.0f;
    float fixed_RS_ = 0.0f;
};

namespace {

void process_one(const std::string& filename,
                 float dt,
                 bool with_mag,
                 const W3dRandomSeeds& seeds,
                 bool write_timeseries,
                 float validation_window_sec)
{
    constexpr float MAG_ODR_HZ = 25.0f;

    auto result = process_wave_file_for_tracker<FusionAdapter_TFG>(
        filename, dt, with_mag, add_noise, MAG_ODR_HZ,
        "_fusion_tfg", "_fusion_tfg_nomag", seeds, write_timeseries);

    if (!result) return;

    if (validation_window_sec > 0.0f) {
        print_validation_metrics(*result, dt, validation_window_sec, "TFG");
    }

    /*
        Worst observed over the eight-record set, with a small margin -- worst
        observed plus about half a percent, rounded up in the last digit the
        channel is quoted in, the rule the OU simulators state.  Yaw is quoted
        to hundredths there and here: a tenth of a degree is 3% of a
        three-degree gate, which is six times the margin the rule asks for.

          channel            TFG worst   this gate   was    OU-III gate
          Z RMS  jonswap        5.21        5.24      5.5      4.72
          Z RMS  pmstokes       5.10        5.13      5.4      4.69
          yaw                   2.92        2.938     3.3      1.068
          3D     jonswap       20.99       21.1      30.6     21.05
          3D     pmstokes      25.78       25.91     68.0     20.83
          acc bias Z            8.84        8.89      9.5      4.93
          acc bias 3D         398.22      400.3     415.0     98.4

        Every bar comes down, and none of that is this file's doing: the
        previous set was fitted before the adaptation-policy work that brought
        the orchestrator's exogeneity timing, commit cadence, r_S floor and
        tau-scaled S=0 cadence up to OU-III's, and the gates were not
        re-derived afterwards. A sentinel carrying ten points of slack on the
        horizontal channel is not a sentinel, so they are re-derived here
        against the same eight records the generated results table reports.

        What that changes about the standing finding. The horizontal channel is
        no longer the outlier it was: 3D error on JONSWAP now sits exactly at
        OU-III's own bar and PM--Stokes within five points of it, where it used
        to be three times worse on the large-wave records. Accelerometer bias
        still is the outlier -- four times OU-III's, and above 100% of the true
        bias on six of the eight records, which means the error exceeds the
        quantity being estimated. That gate records a known deficiency rather
        than endorsing it, and the cause remains unestablished.

        Yaw stays about three times OU-III's, and OU-II's, because the
        continuous magnetic hard-iron correction is carried by both OU families
        and not by this one; the standing heading error it removes is still
        here, and its gauge argument is in the OU-III article.

        bias_3d_percent gates the gyro channel too. The accelerometer sets it
        -- the gyro's worst is 125.60% on pmstokes H4.0 -- so a gyro-bias
        regression has to more than triple before this bar sees it. Splitting
        the two is a change to W3dFailureLimits and to every family that uses
        it, and is deliberately not made here.
    */
    /*
        Each bar is written to whatever precision delivers about half a percent
        over the worst observed value, rather than to a tenth regardless of the
        value: a tenth is 2 percent of a 5.2 and 0.03 percent of a 400, so one
        quantum for every channel means the small ones carry four times the
        margin the rule asks for and the large ones carry none of it.

          channel          worst   this bar   at a tenth   margin
          Z RMS  jonswap    5.21     5.24        5.3        0.60%
          Z RMS  pmstokes   5.10     5.13        5.2        0.58%
          yaw               2.92     2.938       3.0        0.51%
          3D     jonswap   20.99    21.1        21.1        0.52%
          3D     pmstokes  25.78    25.91       26.0        0.52%
          acc bias Z        8.84     8.89        8.9        0.61%
          acc bias 3D     398.22   400.3       400.3        0.52%

        Both an -march=native and an -march=x86-64 build pass all seven; the
        binding records move by up to 3.7e-4 relative between them, so the
        thinnest margin here is 14 times the spread it has to survive.
    */
    static constexpr W3dFailureLimits kRegressionBars{
        .err_limit_percent_z_jonswap   = 5.24f,   // was 5.3,  worst 5.2090 (jonswap H0.27)
        .err_limit_percent_z_pmstokes  = 5.13f,   // was 5.2,  worst 5.1004 (pmstokes H0.27)
        .err_limit_yaw_deg             = 2.938f,  // was 2.94, worst 2.9230 (pmstokes H4.0)
        .err_limit_percent_3d_jonswap  = 21.1f,   // unchanged, worst 20.9914 (jonswap H1.5)
        .err_limit_percent_3d_pmstokes = 25.91f,  // was 26.0, worst 25.7764 (pmstokes H4.0)
        .acc_z_bias_percent            = 8.89f,   // was 8.9,  worst 8.8360 (jonswap H4.0)
        .bias_3d_percent               = 400.3f,  // unchanged, worst 398.2190 (jonswap H4.0, accel)
    };
    static constexpr W3dSummaryLabels kLabels{ .target = "RS_target",
                                               .applied = "RS_applied" };
    print_summary_and_fail_if_needed(*result, dt, kRegressionBars, kLabels);
}

}  // namespace

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    add_noise = true;
    std::vector<std::string> requested_files;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") {
            with_mag = false;
        } else if (arg == "--no-noise") {
            add_noise = false;
        } else if (arg == "--input") {
            if (++i >= argc) {
                std::cerr << "ERROR: --input requires a CSV path\n";
                return 2;
            }
            requested_files.emplace_back(argv[i]);
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--nomag] [--no-noise] [--input PATH]...\n";
            return 0;
        } else {
            std::cerr << "ERROR: unknown argument: " << arg << "\n";
            return 2;
        }
    }

    W3dRandomSeeds seeds;
    try {
        seeds = w3d_random_seeds_from_env();
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    bool write_timeseries = true;
    if (const char* value = std::getenv("W3D_WRITE_TIMESERIES")) {
        write_timeseries = std::string(value) != "0";
    }
    float validation_window_sec = 0.0f;
    env_float("W3D_VALIDATION_WINDOW_SEC", validation_window_sec);

    std::cout << "TFG simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", mag_delay=" << ocean_imu::tfg::MAG_DELAY_SEC
              << " sec, noise=" << (add_noise ? "true" : "false") << "\n";
    std::cout << "RANDOM_SEEDS accel_noise=" << seeds.accel_noise
              << " gyro_noise=" << seeds.gyro_noise
              << " mag_noise=" << seeds.mag_noise
              << " accel_initialization=" << seeds.accel_initialization
              << " gyro_initialization=" << seeds.gyro_initialization
              << " mag_initialization=" << seeds.mag_initialization << "\n";

    const auto files = requested_files.empty()
        ? collect_wave_data_files(".")
        : requested_files;

    try {
        for (const auto& fname : files) {
            process_one(fname, dt, with_mag, seeds, write_timeseries,
                        validation_window_sec);
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) {
        return 1;
    }
    return 0;
}
