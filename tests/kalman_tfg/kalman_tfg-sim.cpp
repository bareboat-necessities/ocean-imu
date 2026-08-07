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

        if (const char* s = std::getenv("TFG_AW_COV_SYNC")) {
            fusion_.setPeriodicAwCovarianceSync(std::string(s) != "off");
        }
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
        Worst observed over the eight-record set, with a small margin.

          channel            TFG worst   this gate   OU-III gate
          Z RMS  jonswap        5.39        5.5          5.4
          Z RMS  pmstokes       5.29        5.4          5.3
          yaw                   3.19        3.3          2.2
          3D     jonswap       30.31       30.6         21.4
          3D     pmstokes      67.63       68.0         22.0
          acc bias Z            9.30        9.5          5.9
          acc bias 3D         412.29      415.0        109.4

        Vertical is at parity with OU-III and six of the eight records beat its
        worst observed value. The horizontal and bias channels are not: 3D
        error degrades sharply on the large-wave records (H4.0 and H8.5), and
        accelerometer-bias error is several times OU-III's.

        Those last two gates are therefore recording a known deficiency rather
        than endorsing it. They belong in the article as a finding, and the
        cause is not yet established -- the pattern points at the horizontal
        channels under large orbital motion, which is where OU-III's own
        article already reports paying for its vertical gain, but that is a
        hypothesis and not a measurement.
    */
    static constexpr W3dFailureLimits kRegressionBars{
        .err_limit_percent_z_jonswap   = 5.5f,
        .err_limit_percent_z_pmstokes  = 5.4f,
        .err_limit_yaw_deg             = 3.3f,
        .err_limit_percent_3d_jonswap  = 30.6f,
        .err_limit_percent_3d_pmstokes = 68.0f,
        .acc_z_bias_percent            = 9.5f,
        .bias_3d_percent               = 415.0f,
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
