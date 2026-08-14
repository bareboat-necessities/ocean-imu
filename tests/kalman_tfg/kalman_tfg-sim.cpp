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
        if (float v = 0.0f; env_float("TFG_ACC_BIAS_UNLOCK_SEC", v)) cfg.acc_bias_unlock_sec = v;

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
        if (bool on = false; env_bool("TFG_AW_COV_SYNC", on)) fusion_.setPeriodicAwCovSync(on);
        {
            float lo = 0.0f, hi = 0.0f;
            if (env_float("TFG_SIGMA_BAND_LOW", lo) && env_float("TFG_SIGMA_BAND_HIGH", hi)) {
                fusion_.setSigmaBandRatios(lo, hi);
            }
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
        Worst observed over the eight-record set, with a small margin -- worst
        observed plus about half a percent, at whatever precision delivers that,
        which is the rule the OU simulators state.  Writing every channel to a
        tenth regardless of its value would give the small ones four times the
        margin the rule asks for and the large ones almost none.

          channel          worst    this bar   margin   previous bar
          Z RMS  jonswap    4.7784   4.803      0.51%    5.24
          Z RMS  pmstokes   4.6846   4.709      0.52%    5.13
          yaw               1.5278   1.536      0.54%    2.938
          3D     jonswap   21.0298  21.14       0.52%   21.1
          3D     pmstokes  20.6045  20.71       0.51%   25.91
          acc bias Z        5.0002   5.026      0.52%    8.89
          acc bias 3D     166.688  167.6        0.55%   400.3

        Both an -march=native and an -march=x86-64 build pass all seven; the
        binding records move by up to 3.5e-4 relative between them, so the
        thinnest margin here is 14 times the spread it has to survive.

        WHY EVERY BAR MOVED.  The orchestrator took on the OU families' magnetic
        acquisition, their continuous hard-iron tracking and their band-passed
        sigma channel, and the tuner coefficients were re-fitted afterwards
        because all three change what those coefficients multiply.  Against the
        previous set, over the same eight records: worst vertical RMS -8.3%,
        worst 3D -18.4%, worst yaw -47.7%, worst accelerometer-bias -58.1%.  The
        gates are re-derived here rather than left where they were, because a
        sentinel carrying a third of its range in slack is not a sentinel.

        What that changes about the standing findings.

        Yaw is no longer the outlier.  It ran about three times OU-III's and
        OU-II's for one reason -- both carried the continuous hard-iron
        correction and this filter did not -- and with that correction here the
        worst record is 1.53 deg against OU-III's 1.068 gate.  The remainder is
        no longer attributable to a missing feature.

        Accelerometer bias is still the outlier, but by far less: 167% against
        OU-III's 98.4, where it used to be 400%.  It is still above 100% of the
        true bias on two of the eight records, which means the error exceeds the
        quantity being estimated, and that gate records a known deficiency
        rather than endorsing it.  Most of what came out was the standing
        attitude error the old one-sample magnetic reference left behind, which
        the bias state was absorbing because the two are only weakly separable
        in waves; what is left has not been attributed.

        3D on JONSWAP is unchanged at OU-III's own bar, and PM-Stokes now sits
        below it where it used to be three times worse on the large-wave
        records.

        bias_3d_percent gates the gyro channel too.  The accelerometer sets it
        -- the gyro's worst is 46.58% on jonswap H8.5 -- so a gyro-bias
        regression has to more than triple before this bar sees it.  Splitting
        the two is a change to W3dFailureLimits and to every family that uses
        it, and is deliberately not made here.

        RE-DERIVED FOR S_factor 1.20 -> 1.00, docs/tfg-adaptation-refit.md.
        The horizontal stationary acceleration prior now matches what the
        records put there (horizontal magnitude 0.997 of vertical, measured on
        all eight).  Six bars move by well under a percent and one of them --
        yaw -- had to be raised past a bar the filter no longer clears, which is
        the case that needs its own argument rather than a table row.

          channel          worst    this bar   margin   previous bar
          Z RMS  jonswap    4.7830   4.807      0.50%    4.803
          Z RMS  pmstokes   4.6833   4.707      0.51%    4.709
          yaw               1.5818   1.59       0.52%    1.536   <-- raised
          3D     jonswap   21.0169  21.13       0.54%   21.14
          3D     pmstokes  20.6322  20.74       0.52%   20.71
          acc bias Z        4.9961   5.022      0.52%    5.026
          acc bias 3D     163.607  164.5        0.55%   167.6

        THE YAW BAR GOES UP 3.5% AND YAW DOES NOT GET WORSE.  Pooled over the
        eight records at six IMU seed triplets, yaw RMS is 0.9578 of the shipped
        filter's; on the record that now binds, jonswap H8.5, it is 0.9109.  The
        deterministic protocol disagrees because that record's yaw spans 1.03 to
        4.56 deg across those six seeds -- a 4.4x spread -- and the default seed
        happens to draw the smallest of the six under the old constant:

          seed      S=1.20   S=1.00
          default   1.0322   1.5818
          7         3.7871   3.1882
          99        1.3618   0.7669
          3         2.9330   3.5454
          21        2.4334   1.8283
          55        4.5617   3.9482

        Four of six improve and the six-seed mean falls 2.685 -> 2.476 deg.  The
        sentinel follows the deterministic protocol it is written against; the
        quality claim rests on the seeds.  This is the same shape, and the same
        argument, as the yaw sentinel in OU-III's anisotropy study, which moved
        21% on the same kind of draw.

        The other six all come down or stay put.  Worth naming: accelerometer
        bias 3D goes 166.688 -> 163.607 on the binding record and 80.97 ->
        116.88 on jonswap H8.5, which is a redistribution across records rather
        than a uniform gain -- the pooled six-seed figure is 0.9209.
    */
    static constexpr W3dFailureLimits kRegressionBars{
        .err_limit_percent_z_jonswap   = 4.807f,  // was 4.803, worst 4.7830 (jonswap H0.27)
        .err_limit_percent_z_pmstokes  = 4.707f,  // was 4.709, worst 4.6833 (pmstokes H0.27)
        .err_limit_yaw_deg             = 1.59f,   // was 1.536, worst 1.5818 (jonswap H8.5)
        .err_limit_percent_3d_jonswap  = 21.13f,  // was 21.14, worst 21.0169 (jonswap H1.5)
        .err_limit_percent_3d_pmstokes = 20.74f,  // was 20.71, worst 20.6322 (pmstokes H4.0)
        .acc_z_bias_percent            = 5.022f,  // was 5.026, worst 4.9961 (pmstokes H8.5)
        .bias_3d_percent               = 164.5f,  // was 167.6, worst 163.607 (pmstokes H4.0, accel)
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
