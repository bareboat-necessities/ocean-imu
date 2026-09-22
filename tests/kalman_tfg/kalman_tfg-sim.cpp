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

    One of them -- accelerometer bias -- still records behaviour that is plainly
    poor, and it is set where it is so that further regression is caught, not
    because the current value is acceptable. See the table below.
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

// Feature ablations. Each names one of the OU-III behaviours this filter took
// on, so the study can turn exactly one off and price it rather than inferring
// its contribution from the combined result.
bool env_bool(const char* name, bool& out) {
    if (const char* s = std::getenv(name)) {
        const std::string v(s);
        out = !(v == "0" || v == "false" || v == "off");
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
    // TFG's own MEKF sensor variances, as multiples of the ones the shared
    // harness hands every family (2.8x injected accel white, 2.0x gyro, 1.2x
    // mag).  Set by TFG's own MEKF-variance sweep; 1.0 here means the harness
    // value is used unchanged.
    static constexpr float SIGMA_A_RESCALE = 1.0f;
    static constexpr float SIGMA_G_RESCALE = 0.1f;
    static constexpr float SIGMA_M_RESCALE = 4.0f;

    FusionAdapter_TFG(bool with_mag,
                      const Vector3f& sigma_a_init,
                      const Vector3f& sigma_g,
                      const Vector3f& sigma_m)
    {
        Fusion::Config cfg;
        cfg.with_mag = with_mag;
        cfg.sigma_a = sigma_a_init * SIGMA_A_RESCALE;
        cfg.gyro_noise_density = sigma_g.x() * SIGMA_G_RESCALE;
        cfg.sigma_m = sigma_m * SIGMA_M_RESCALE;

        // The MEKF variances, swept for this family.  The three
        // sensor sigmas are scale factors on the deployed point above, so a
        // scale of 1 leaves it in place; the other two are absolute.
        if (float v = 0.0f; env_float("SF_SIGMA_A_SCALE", v)) cfg.sigma_a *= v;
        if (float v = 0.0f; env_float("SF_SIGMA_G_SCALE", v)) cfg.gyro_noise_density *= v;
        if (float v = 0.0f; env_float("SF_SIGMA_M_SCALE", v)) cfg.sigma_m *= v;
        if (float v = 0.0f; env_float("SF_GYRO_BIAS_RW_VAR", v)) cfg.gyro_bias_rw_var = v;
        if (float v = 0.0f; env_float("SF_PQ0", v)) cfg.initial_covariance = v;
        cfg.mag_delay_sec = ocean_imu::tfg::MAG_DELAY_SEC;
        if (float v = 0.0f; env_float("SF_MAG_DELAY_SEC", v)) cfg.mag_delay_sec = v;
        cfg.freeze_acc_bias_until_live = true;
        cfg.Racc_warmup_std = 0.5f;

        if (const char* p = std::getenv("TFG_STARTUP_INIT")) {
            const std::string policy(p);
            if (policy == "staged") {
                cfg.startup_init_policy = Fusion::StartupInitPolicy::StagedMekf;
            } else if (policy == "proxy") {
                cfg.startup_init_policy = Fusion::StartupInitPolicy::MahonyProxy;
            } else {
                throw std::runtime_error("unknown TFG_STARTUP_INIT: " + policy);
            }
        }
        bool flag = false;
        if (env_bool("TFG_MAG_REFINE", flag))     cfg.mag_refine_enabled = flag;
        if (env_bool("TFG_MAG_HARD_IRON", flag))  cfg.mag_continuous_hard_iron = flag;
        if (env_bool("TFG_WAVE_BAND", flag))      cfg.wave_band_tuning = flag;
        if (float v = 0.0f; env_float("TFG_MAG_MIN_WINDOW_SEC", v)) cfg.mag_min_window_sec = v;
        if (float v = 0.0f; env_float("TFG_MAG_REFINE_START_SEC", v)) cfg.mag_refine_start_sec = v;
        if (float v = 0.0f; env_float("TFG_MAG_REFINE_WINDOW_SEC", v)) cfg.mag_refine_window_sec = v;
        if (float v = 0.0f; env_float("TFG_MAG_HI_SLEW_TAU_SEC", v)) cfg.mag_hi_slew_tau_sec = v;
        if (float v = 0.0f; env_float("TFG_MAG_HI_RIDGE_REL", v)) cfg.mag_hi_model_ridge_relative = v;
        if (float v = 0.0f; env_float("TFG_MAG_HI_MIN_INFO", v)) cfg.mag_hi_min_information = v;
        if (float v = 0.0f; env_float("TFG_ACC_BIAS_UNLOCK_SEC", v)) cfg.acc_bias_unlock_sec = v;

        // Out-of-band accelerometer guard ahead of the proxy and the MEKF.
        // Armed by default at the deployed corner; these override it, and a
        // zero cutoff removes it entirely.
        if (float v = 0.0f; env_float("TFG_ACC_GUARD_HZ", v)) cfg.acc_vibration_guard_hz = v;
        if (const char* g = std::getenv("TFG_ACC_GUARD_POLES")) {
            cfg.acc_vibration_guard_poles = std::atoi(g);
        }
        if (float v = 0.0f; env_float("TFG_ACC_GUARD_RACC_GAIN", v)) {
            cfg.acc_vibration_racc_gain = v;
        }

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
        if (env_float("TFG_R_S_X_FACTOR", v))   fusion_.setRSXFactor(v);
        if (env_float("TFG_R_S_Y_FACTOR", v))   fusion_.setRSYFactor(v);
        if (env_float("TFG_ADAPT_TAU_SEC", v))  fusion_.setAdaptationTimeConstants(v);
        if (env_float("TFG_ADAPT_RS_MULT", v))  fusion_.setRSAdaptMult(v);
        if (env_float("TFG_ADAPT_RS_SLEW_LOG", v)) fusion_.setRSAdaptSlewLog(v);
        if (env_float("TFG_ACC_NOISE_FLOOR", v))   fusion_.setAccNoiseFloorSigma(v);
        if (bool on = false; env_bool("TFG_AW_COV_SYNC", on)) fusion_.setPeriodicAwCovSync(on);
        {
            float lo = 0.0f, hi = 0.0f;
            if (env_float("TFG_SIGMA_BAND_LOW", lo) && env_float("TFG_SIGMA_BAND_HIGH", hi)) {
                fusion_.setSigmaBandRatios(lo, hi);
            }
        }
        {
            float engage_lo = 0.0f, engage_hi = 0.0f, engage_tau = 0.0f;
            const bool lo_set = env_float("TFG_ACC_GUARD_ENGAGE_LO", engage_lo);
            const bool hi_set = env_float("TFG_ACC_GUARD_ENGAGE_HI", engage_hi);
            const bool tau_set = env_float("TFG_ACC_GUARD_ENGAGE_TAU", engage_tau);
            if (lo_set || hi_set || tau_set) {
                fusion_.setAccelVibrationEngagement(
                    lo_set ? engage_lo : -1.0f,
                    hi_set ? engage_hi : -1.0f,
                    tau_set ? engage_tau : -1.0f);
            }
        }
    }

    // Vibration-guard telemetry, so a replay can say whether the guard engaged
    // and how much out-of-band accelerometer content it was seeing.  The shared
    // runner owns the adapter, so end-of-record is the destructor; silent
    // unless a cutoff was configured, which keeps an unguarded run unchanged.
    ~FusionAdapter_TFG() override {
        if (!(fusion_.accelVibrationGuardCutoffHz() > 0.0f)) return;
        std::cout << "ACC_GUARD cutoff_hz=" << fusion_.accelVibrationGuardCutoffHz()
                  << " poles=" << fusion_.accelVibrationGuardPoles()
                  << " engagement=" << fusion_.accelVibrationGuardEngagement()
                  << " out_of_band_rms_mps2=" << fusion_.accelVibrationRms()
                  << " delay_sec=" << fusion_.accelVibrationGuardDelaySec()
                  << " racc_gain=" << fusion_.accelVibrationRaccGain()
                  << " racc_std_mps2=" << fusion_.accelVibrationRaccStd().x()
                  << "\n";
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

    static constexpr W3dFailureLimits kRegressionBars{
        .err_limit_percent_z_jonswap   = 4.405f,  // worst 4.38306 (jonswap H0.27)
        .err_limit_percent_z_pmstokes  = 4.356f,  // worst 4.33397 (pmstokes H0.27)
        .err_limit_yaw_deg             = 0.943f,  // worst 0.937382 (jonswap H4.0)
        .err_limit_percent_3d_jonswap  = 9.83f,   // worst 9.77897 (jonswap H8.5)
        .err_limit_percent_3d_pmstokes = 9.92f,   // worst 9.86933 (pmstokes H8.5)
        .acc_z_bias_percent            = 4.494f,  // worst 4.47079 (pmstokes H8.5)
        .bias_3d_percent               = 74.1f,   // worst 73.727 (pmstokes H8.5, accel)
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
