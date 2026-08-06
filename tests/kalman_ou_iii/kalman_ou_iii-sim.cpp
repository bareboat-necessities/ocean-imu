#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
*/

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/W3dSimCommon.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;
using Eigen::Matrix3f;

bool add_noise = true;
bool attitude_only = false;

namespace {

bool env_float(const char* name, float& out)
{
    if (const char* s = std::getenv(name)) {
        out = static_cast<float>(std::atof(s));
        return true;
    }
    return false;
}

bool env_int(const char* name, int& out)
{
    if (const char* s = std::getenv(name)) {
        out = std::atoi(s);
        return true;
    }
    return false;
}

} // namespace

class FusionAdapter_OU_III final : public IW3dFusionAdapter {
public:
    FusionAdapter_OU_III(bool with_mag,
                         const Vector3f& sigma_a_init,
                         const Vector3f& sigma_g,
                         const Vector3f& sigma_m)
        : with_mag_(with_mag)
    {
        cfg_.with_mag = with_mag;
        cfg_.sigma_a = sigma_a_init;
        cfg_.sigma_g = sigma_g;
        cfg_.sigma_m = sigma_m;
        cfg_.mag_delay_sec = MAG_DELAY_SEC;
        cfg_.freeze_acc_bias_until_live = true;
        cfg_.Racc_warmup_std = 0.5f;

        apply_env_overrides();
        load_fixed_tuning();

        fusion_.begin(cfg_);
        auto& filter = fusion_.raw();


        const std::string aw_cov_sync = load_aw_cov_sync_policy();
        filter.setPeriodicAwCovarianceSync(aw_cov_sync != "reconfigure");
        filter.setAwCovarianceSyncCongruent(aw_cov_sync == "congruent");

        if (attitude_only) {
            filter.enableLinearBlock(false);
            filter.mekf().set_initial_acc_bias(Vector3f::Zero());
            filter.mekf().set_initial_acc_bias_std(0.0f);
            filter.mekf().set_Q_bacc_rw(Vector3f::Zero());
            filter.mekf().set_Racc_std(Vector3f::Constant(0.4f));
        } else {
            filter.enableLinearBlock(true);
            filter.enableTuner(true);
            filter.enableClamp(true);

            float v = 0.0f;

            // Generic OU_* names are accepted for compatibility.
            // OU_III_* names are applied afterward and win if both are set.

            if (env_float("OU_TAU_COEFF", v)) {
                filter.setTauCoeff(v);
            }
            if (env_float("OU_III_TAU_COEFF", v)) {
                filter.setTauCoeff(v);
            }

            if (env_float("OU_SIGMA_COEFF", v)) {
                filter.setSigmaCoeff(v);
            }
            if (env_float("OU_III_SIGMA_COEFF", v)) {
                filter.setSigmaCoeff(v);
            }

            // OU_III-specific horizontal stationary accel anisotropy.
            // This maps to SeaStateFusionFilter_OU_III::setSFactor().
            if (env_float("OU_III_S_FACTOR", v)) {
                filter.setSFactor(v);
            }

            // OU_III-specific R_S anisotropy and coefficient.
            // These are the real setter names in SeaStateFusionFilter_OU_III.
            if (env_float("OU_III_R_S_XY_FACTOR", v)) {
                filter.setRSXYFactor(v);
            }

            if (env_float("OU_III_R_S_COEFF", v)) {
                filter.setRSCoeff(v);
            }

            // NOTE:
            // SeaStateFusionFilter_OU_III does not expose a V0/R_v0 coefficient setter.
            // Therefore OU_III_R_V0_COEFF is intentionally not read here.

            if (env_float("OU_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }
            if (env_float("OU_III_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }

            if (env_float("OU_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }
            if (env_float("OU_III_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }

            // Smoothing horizon of the r_S EMA, in units of tau_target.
            if (env_float("OU_ADAPT_RS_MULT", v)) {
                filter.setRSAdaptMult(v);
            }
            if (env_float("OU_III_ADAPT_RS_MULT", v)) {
                filter.setRSAdaptMult(v);
            }

            // Discrepancy threshold that shortens that horizon when the sea
            // state actually moves.  0 keeps the plain proportional horizon.
            if (env_float("OU_ADAPT_RS_SLEW_LOG", v)) {
                filter.setRSAdaptSlewLog(v);
            }
            if (env_float("OU_III_ADAPT_RS_SLEW_LOG", v)) {
                filter.setRSAdaptSlewLog(v);
            }

            if (env_float("OU_ADAPT_EVERY_SECS", v)) {
                filter.setAdaptationUpdatePeriod(v);
            }
            if (env_float("OU_III_ADAPT_EVERY_SECS", v)) {
                filter.setAdaptationUpdatePeriod(v);
            }

            if (env_float("OU_FREQ_INPUT_CUTOFF_HZ", v)) {
                filter.setFreqInputCutoffHz(v);
            }
            if (env_float("OU_III_FREQ_INPUT_CUTOFF_HZ", v)) {
                filter.setFreqInputCutoffHz(v);
            }

            // Scalar accel-bias initialization uncertainty only.
            // Bias vector X/Y/Z env overrides intentionally removed.
            // Accelerometer-bias random walk.  The bias competes with the OU
            // acceleration for the low-frequency content, and the wave-band
            // operating point moves the OU corner down toward it, so this is
            // the knob that prices that competition.
            if (env_float("OU_III_ACC_BIAS_RW", v)) {
                filter.mekf().set_Q_bacc_rw(Eigen::Vector3f::Constant(v));
            }

            // Ablate the wave-band operating point back to the
            // acceleration-band frequency the filter used before.
            if (const char* band = std::getenv("W3D_TUNING_BAND")) {
                const std::string value = band;
                if (value == "acceleration") {
                    filter.setWaveBandTuning(false);
                } else if (value != "wave") {
                    throw std::runtime_error(
                        "W3D_TUNING_BAND must be wave or acceleration");
                }
            }

            // Ablate the wave-period estimator's input away from the leveled
            // vertical acceleration (default), which passes through the
            // attitude solution.  Both alternatives are measurement-only and
            // open the tuner coupling: "body_z" is the raw proxy the frequency
            // tracker already runs on, "complementary" levels it with a
            // private Mahony observer.
            if (const char* src = std::getenv("W3D_WAVE_PERIOD_INPUT")) {
                const std::string value = src;
                if (value == "body_z") {
                    filter.setWavePeriodInput(WavePeriodInputSource::BodyZ);
                } else if (value == "complementary") {
                    filter.setWavePeriodInput(
                        WavePeriodInputSource::Complementary);
                } else if (value != "leveled") {
                    throw std::runtime_error(
                        "W3D_WAVE_PERIOD_INPUT must be leveled, body_z or "
                        "complementary");
                }
            }

            // Gains of that private observer, so the correction corner can be
            // swept against the wave band it must stay below.
            {
                float two_kp = 0.2f, two_ki = 0.0f;
                const bool kp = env_float("W3D_WAVE_PERIOD_MAHONY_KP", two_kp);
                const bool ki = env_float("W3D_WAVE_PERIOD_MAHONY_KI", two_ki);
                if (kp || ki) {
                    filter.setWavePeriodComplementaryGains(two_kp, two_ki);
                }
            }

            if (env_float("OU_ACC_BIAS_INIT_STD", v)) {
                filter.mekf().set_initial_acc_bias_std(v);
            }
            if (env_float("OU_III_ACC_BIAS_INIT_STD", v)) {
                filter.mekf().set_initial_acc_bias_std(v);
            }
        }
    }

    void apply_env_overrides() {
        float vf = 0.0f;
        int vi = 0;

        if (env_float("SF_MAG_DELAY_SEC", vf)) cfg_.mag_delay_sec = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_MAX_SIN", vf)) cfg_.mag_gravity_align_max_sin = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_HOLD_SEC", vf)) cfg_.mag_gravity_align_hold_sec = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_LPF_TAU", vf)) cfg_.mag_gravity_align_lpf_tau = vf;
        if (env_float("SF_MAG_TILT_FALLBACK_SEC", vf)) cfg_.mag_tilt_fallback_sec = vf;
        if (env_float("SF_MAG_EXTREME_GYRO_DPS", vf)) cfg_.mag_extreme_gyro_dps = vf;
        if (env_float("SF_MAG_INIT_MIN_MAG_NORM", vf)) cfg_.mag_init_min_mag_norm = vf;
        if (env_int("SF_MAG_MIN_SAMPLES", vi)) cfg_.mag_min_samples = vi;
        if (env_float("SF_MAG_MIN_WINDOW_SEC", vf)) cfg_.mag_min_window_sec = vf;

        if (env_float("SF_RACC_WARMUP_STD", vf)) cfg_.Racc_warmup_std = vf;
        if (env_float("SF_ONLINE_TUNE_WARMUP_SEC", vf)) cfg_.online_tune_warmup_sec = vf;

        if (env_float("SF_BOOT_TILT_ACC_TAU", vf)) cfg_.bootstrap_tilt_obs_acc_tau_sec = vf;
        if (env_float("SF_BOOT_GRAV_SLOW_TAU", vf)) cfg_.bootstrap_gravity_slow_tau_sec = vf;
        if (env_float("SF_BOOT_GRAV_ALIGN_MAX_SIN", vf)) cfg_.bootstrap_gravity_align_max_sin = vf;
        if (env_float("SF_BOOT_GRAV_HOLD_SEC", vf)) cfg_.bootstrap_gravity_hold_sec = vf;
        if (env_float("SF_BOOT_GRAV_MIN_SEC", vf)) cfg_.bootstrap_gravity_min_sec = vf;
        if (env_float("SF_BOOT_GRAV_TIMEOUT_SEC", vf)) cfg_.bootstrap_gravity_timeout_sec = vf;
        if (env_float("SF_BOOT_GRAV_NORM_FRAC", vf)) cfg_.bootstrap_gravity_norm_frac = vf;
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
        auto& filter = fusion_.raw();
        if (tuning_ == TuningMode::Adaptive) return;
        if (fixed_tuning_applied_ || !filter.isAdaptiveLive()) return;

        const bool ok = (tuning_ == TuningMode::Fixed)
            ? filter.setFixedTuning(fixed_tau_s_, fixed_sigma_a_, fixed_RS_)
            : filter.setChannelFreeze(
                  tuning_ == TuningMode::AdaptiveRSOnly,
                  fixed_tau_s_,
                  fixed_sigma_a_,
                  tuning_ == TuningMode::AdaptiveOUOnly,
                  fixed_RS_);
        if (!ok) throw std::runtime_error("invalid OU-III tuning point");
        fixed_tuning_applied_ = true;
    }

    FilterSnapshot snapshot() const override {
        const auto& filter = fusion_.raw();
        const auto& d = filter.dir();

        FilterSnapshot s;
        s.disp_est_zu = ned_to_zu(filter.mekf().get_position());
        s.vel_est_zu  = ned_to_zu(filter.mekf().get_velocity());
        s.acc_est_zu  = ned_to_zu(filter.mekf().get_world_accel());

        // Filter attitude is BODY->WORLD in NED.
        //
        // In IMU-only mode, after mag lock this is BODY->WORLD in the learned
        // magnetic-NED frame, not true-north NED.
        //
        // Do not apply WMM/declination correction here. A real IMU does not know true
        // north unless an external declination/location model is explicitly supplied.
        const Quaternionf q_bw_ned = filter.mekf().quaternion_boat().normalized();

        float roll_deg  = 0.0f;
        float pitch_deg = 0.0f;
        float yaw_deg   = 0.0f;
        quat_to_euler_nautical(q_bw_ned, roll_deg, pitch_deg, yaw_deg);

        s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, wrapDeg(yaw_deg));

        s.acc_bias_est_ned    = filter.mekf().get_acc_bias();
        s.gyro_bias_est_ned   = filter.mekf().gyroscope_bias();
        s.mag_bias_est_ned_uT = get_mag_bias_est_uT(filter.mekf());

        s.tau_target      = filter.getTauTarget();
        s.sigma_target    = filter.getSigmaTarget();
        s.tuning_target   = filter.getRSTarget();

        s.tau_applied     = filter.getTauApplied();
        s.sigma_applied   = filter.getSigmaApplied();
        s.tuning_applied  = filter.getRSApplied();

        s.freq_hz         = filter.getFreqHz();
        s.wave_period_sec = filter.getWavePeriodSec();
        s.period_sec      = filter.getPeriodSec();
        s.accel_variance  = filter.getAccelVariance();

        s.displacement_scale_m = filter.getDisplacementScale();
        s.velocity_scale_mps   = filter.getVerticalSpeedEnvelopeMps(true);

        s.direction.phase = d.getPhase();
        s.direction.direction_deg = d.getAxisDegrees();
        s.direction.apparent_to_deg = filter.getApparentWaveDirectionToDeg();
        s.direction.apparent_from_deg = filter.getApparentWaveDirectionFromDeg();
        s.direction.sense_coherence = filter.getDirSenseCoherence();
        s.direction.direction_deg_generator_signed = dirDegGeneratorSignedFromVec(d.getAxis());
        s.direction.uncertainty_deg = d.getAxisUncertaintyDegrees();
        s.direction.confidence = d.getLastStableConfidence();
        s.direction.amplitude = d.getAmplitude();
        s.direction.direction_vec = d.getAxis();
        s.direction.filtered_signal = d.getFilteredSignal();

        constexpr float CONF_THRESH = 20.0f;
        constexpr float AMP_THRESH  = 0.08f;

        if (s.direction.confidence > CONF_THRESH && s.direction.amplitude > AMP_THRESH) {
            s.direction.sign = filter.getDirSignState();
            s.direction.sign_num =
                (s.direction.sign == FORWARD) ? 1 :
                (s.direction.sign == BACKWARD ? -1 : 0);

            // Physical directed propagation vector.  The class above is
            // relative to the estimator's own axis representative; this is not.
            const float travel_x = filter.dir_sign().getDirectedX();
            const float travel_y = filter.dir_sign().getDirectedY();
            if (std::isfinite(travel_x) && std::isfinite(travel_y)) {
                s.direction.travel_vec_boat = Eigen::Vector2f(travel_x, travel_y);
            }
        }

        return s;
    }

private:
    // Selects the a_w covariance-synchronization policy under test.
    // "periodic" (default, and the deployed policy) re-aligns the
    // latent-acceleration marginal with its stationary prior once per
    // adaptation period; "reconfigure" restricts that to discrete
    // reconfiguration events and is the matched ablation; "congruent" keeps the
    // periodic cadence but performs the re-alignment as a congruence, which is
    // the posterior-consistent version of the same operation.
    static std::string load_aw_cov_sync_policy()
    {
        const char* raw = std::getenv("W3D_AW_COV_SYNC");
        const std::string policy = (raw && *raw) ? raw : "periodic";
        if (policy == "reconfigure" || policy == "periodic" ||
            policy == "congruent") {
            return policy;
        }
        throw std::runtime_error(
            "W3D_AW_COV_SYNC must be reconfigure, periodic, or congruent");
    }

    // W3D_TUNING_MODE selects how much of the operating point is estimated
    // online:
    //   adaptive          all three channels track the sea (deployed filter)
    //   fixed*            all three frozen at the supplied triple
    //   adaptive_rs_only  tau and sigma_aw frozen, r_S keeps adapting
    //   adaptive_ou_only  r_S frozen, tau and sigma_aw keep adapting
    // The last two isolate which channel carries the adaptation benefit; the
    // deployed law r_S = clip(c*sigma_aw*tau^3) ties them together, so a
    // fixed-versus-adaptive comparison alone cannot separate them.
    void load_fixed_tuning()
    {
        const char* raw_mode = std::getenv("W3D_TUNING_MODE");
        const std::string mode = raw_mode ? raw_mode : "adaptive";
        if (mode == "adaptive") return;
        if (mode == "adaptive_rs_only") {
            tuning_ = TuningMode::AdaptiveRSOnly;
        } else if (mode == "adaptive_ou_only") {
            tuning_ = TuningMode::AdaptiveOUOnly;
        } else if (mode.starts_with("fixed")) {
            tuning_ = TuningMode::Fixed;
        } else {
            throw std::runtime_error(
                "W3D_TUNING_MODE must be adaptive, adaptive_rs_only, "
                "adaptive_ou_only, or start with fixed");
        }

        const bool have_all =
            env_float("W3D_FIXED_TAU_S", fixed_tau_s_) &&
            env_float("W3D_FIXED_SIGMA_A", fixed_sigma_a_) &&
            env_float("W3D_FIXED_RS", fixed_RS_);
        if (!have_all ||
            !(std::isfinite(fixed_tau_s_) && fixed_tau_s_ > 0.0f &&
              std::isfinite(fixed_sigma_a_) && fixed_sigma_a_ > 0.0f &&
              std::isfinite(fixed_RS_) && fixed_RS_ > 0.0f))
        {
            throw std::runtime_error(
                mode + " OU-III mode requires positive W3D_FIXED_TAU_S, "
                "W3D_FIXED_SIGMA_A, and W3D_FIXED_RS");
        }
    }

    enum class TuningMode {
        Adaptive,
        Fixed,
        AdaptiveRSOnly,
        AdaptiveOUOnly,
    };

    bool with_mag_ = true;
    TuningMode tuning_ = TuningMode::Adaptive;
    bool fixed_tuning_applied_ = false;
    float fixed_tau_s_ = NAN;
    float fixed_sigma_a_ = NAN;
    float fixed_RS_ = NAN;
    using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
    mutable Fusion fusion_;
    Fusion::Config cfg_{};
};

// Regression sentinels for the deterministic single-realization protocol, not
// targets.  Each is the worst value the current filter produces across the
// scored records plus about half a percent, rounded up to the next tenth.
//
// That margin is deliberately small because the metrics are deterministic: the
// same records and seeds under -march=native, x86-64 and x86-64-v2 agree to
// within 6e-6 relative, so a limit this close only trips when the filter
// actually gets worse.  Setting one below what the filter currently achieves
// makes it fail every run rather than catching a regression.
//
// Re-derived for the 900 s scoring window: a sentinel fitted to the
// previous 60 s window is not a sentinel for this one, it is just a number the
// filter passes by a wide margin.
//
// bias_3d_percent re-derived again when the r_S smoothing horizon was
// shortened from 5 to 3 wave-period-halves.  That gate is dominated by the
// horizontal accelerometer bias on jonswap H1.5, which is already unobservable
// there -- the error exceeds the true bias with either horizon -- and the
// shorter horizon moves the aggregate from 106.2 to 108.9 while the
// displacement error on the same record is unchanged.  See
// docs/ou-ema-adaptation-tuning.md.
static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 5.4f,    // worst 5.34 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 5.3f,    // worst 5.25 (pmstokes H0.27)
    .err_limit_yaw_deg             = 2.2f,    // worst 2.17 (pmstokes H8.5)
    .err_limit_percent_3d_jonswap  = 21.4f,   // worst 21.23 (jonswap H1.5)
    .err_limit_percent_3d_pmstokes = 22.0f,   // worst 21.85 (pmstokes H0.27)
    .acc_z_bias_percent            = 5.9f,    // worst 5.84 (pmstokes H8.5)
    .bias_3d_percent               = 109.4f,  // worst 108.88 (jonswap H1.5, accel)
};

static constexpr W3dSummaryLabels SUMMARY_LABELS{
    .target  = "RS_target",
    .applied = "RS_applied",
};

static void process_wave_file_for_tracker(const std::string& filename,
                                          float dt,
                                          bool with_mag,
                                          const W3dRandomSeeds& seeds,
                                          bool write_timeseries,
                                          float validation_window_sec)
{
    constexpr float MAG_ODR_HZ = 25.0f;

    auto result = process_wave_file_for_tracker<FusionAdapter_OU_III>(
        filename,
        dt,
        with_mag,
        add_noise,
        MAG_ODR_HZ,
        "_fusion_ou3",
        "_fusion_ou3_nomag",
        seeds,
        write_timeseries);

    if (!result) return;

    if (validation_window_sec > 0.0f) {
        print_validation_metrics(*result, dt, validation_window_sec, "OU_III");
    }
    print_summary_and_fail_if_needed(*result, dt, FAIL_LIMITS, SUMMARY_LABELS);
}

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    add_noise = true;
    std::vector<std::string> requested_files;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

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

    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", mag_delay=" << MAG_DELAY_SEC
              << " sec, noise=" << (add_noise ? "true" : "false")
              << "\n";
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
            process_wave_file_for_tracker(
                fname, dt, with_mag, seeds, write_timeseries,
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
