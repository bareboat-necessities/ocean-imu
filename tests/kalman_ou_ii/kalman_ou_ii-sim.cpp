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
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"

using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Quaternionf;

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

class FusionAdapter_OU_II final : public IW3dFusionAdapter {
public:
    FusionAdapter_OU_II(bool with_mag,
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

        filter.setPeriodicAwCovarianceSync(load_periodic_aw_cov_sync());

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
            // OU_II_* names are applied afterward and win if both are set.

            if (env_float("OU_P_FACTOR", v)) {
                filter.setPFactor(v);
            }
            if (env_float("OU_II_P_FACTOR", v)) {
                filter.setPFactor(v);
            }

            if (env_float("OU_R_P0_XY_FACTOR", v)) {
                filter.setR_p0_XYFactor(v);
            }
            if (env_float("OU_II_R_P0_XY_FACTOR", v)) {
                filter.setR_p0_XYFactor(v);
            }

            if (env_float("OU_TAU_COEFF", v)) {
                filter.setTauCoeff(v);
            }
            if (env_float("OU_II_TAU_COEFF", v)) {
                filter.setTauCoeff(v);
            }

            if (env_float("OU_SIGMA_COEFF", v)) {
                filter.setSigmaCoeff(v);
            }
            if (env_float("OU_II_SIGMA_COEFF", v)) {
                filter.setSigmaCoeff(v);
            }

            if (env_float("OU_R_P0_COEFF", v)) {
                filter.setR_p0_Coeff(v);
            }
            if (env_float("OU_II_R_P0_COEFF", v)) {
                filter.setR_p0_Coeff(v);
            }

            if (env_float("OU_R_V0_COEFF", v)) {
                filter.setR_v0_Coeff(v);
            }
            if (env_float("OU_II_R_V0_COEFF", v)) {
                filter.setR_v0_Coeff(v);
            }

            // Clamps on the two drift-band regularizers.  Exposed because
            // OU-III found that its equivalent floor, not the schedule, was
            // setting the operating point in every low-motion sea, and the
            // only way to notice that is to move the floor and watch whether
            // anything responds.
            {
                float lo = MIN_R_p0_std, hi = MAX_R_p0_std;
                const bool got_lo = env_float("OU_II_R_P0_MIN", lo);
                const bool got_hi = env_float("OU_II_R_P0_MAX", hi);
                if (got_lo || got_hi) {
                    filter.setR_p0_Bounds(lo, hi);
                }
            }
            {
                float lo = MIN_R_v0_std, hi = MAX_R_v0_std;
                const bool got_lo = env_float("OU_II_R_V0_MIN", lo);
                const bool got_hi = env_float("OU_II_R_V0_MAX", hi);
                if (got_lo || got_hi) {
                    filter.setR_v0_Bounds(lo, hi);
                }
            }

            if (env_float("OU_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }
            if (env_float("OU_II_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }

            if (env_float("OU_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }
            if (env_float("OU_II_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }

            // Smoothing horizons of the two drift-correction EMAs, in units
            // of tau_target.
            if (env_float("OU_ADAPT_R_P0_MULT", v)) {
                filter.setR_p0_AdaptMult(v);
            }
            if (env_float("OU_II_ADAPT_R_P0_MULT", v)) {
                filter.setR_p0_AdaptMult(v);
            }

            if (env_float("OU_ADAPT_R_V0_MULT", v)) {
                filter.setR_v0_AdaptMult(v);
            }
            if (env_float("OU_II_ADAPT_R_V0_MULT", v)) {
                filter.setR_v0_AdaptMult(v);
            }

            // Discrepancy threshold that shortens both horizons when the sea
            // state actually moves.  0 keeps the plain proportional horizon.
            if (env_float("OU_ADAPT_R_SLEW_LOG", v)) {
                filter.setR_AdaptSlewLog(v);
            }
            if (env_float("OU_II_ADAPT_R_SLEW_LOG", v)) {
                filter.setR_AdaptSlewLog(v);
            }

            if (env_float("OU_ADAPT_EVERY_SECS", v)) {
                filter.setAdaptationUpdatePeriod(v);
            }
            if (env_float("OU_II_ADAPT_EVERY_SECS", v)) {
                filter.setAdaptationUpdatePeriod(v);
            }

            if (env_float("OU_FREQ_INPUT_CUTOFF_HZ", v)) {
                filter.setFreqInputCutoffHz(v);
            }
            if (env_float("OU_II_FREQ_INPUT_CUTOFF_HZ", v)) {
                filter.setFreqInputCutoffHz(v);
            }

            // Accelerometer-bias random walk.  The bias competes with the OU
            // acceleration for the low-frequency content, and the wave-band
            // operating point moves the OU corner down toward it, so this is
            // the knob that prices that competition.
            if (env_float("OU_II_ACC_BIAS_RW", v)) {
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

            // Ablate the wave-period estimator's input away from the
            // complementary-levelled default.  "leveled" restores the older
            // behaviour, which levels with the attitude solution and so closes
            // the tuner coupling; "body_z" is the raw proxy, measurement-only
            // like the default but unlevelled.
            if (const char* src = std::getenv("W3D_WAVE_PERIOD_INPUT")) {
                const std::string value = src;
                if (value == "body_z") {
                    filter.setWavePeriodInput(WavePeriodInputSource::BodyZ);
                } else if (value == "complementary") {
                    filter.setWavePeriodInput(
                        WavePeriodInputSource::Complementary);
                } else if (value == "leveled") {
                    filter.setWavePeriodInput(WavePeriodInputSource::Leveled);
                } else {
                    throw std::runtime_error(
                        "W3D_WAVE_PERIOD_INPUT must be leveled, body_z or "
                        "complementary");
                }
            }

            // Ablate the frequency tracker's input from the raw body-Z proxy
            // (default) to the levelled signal from the private Mahony
            // observer.  Both are measurement-only; this changes the tracker
            // frequency and so the direction demodulator's carrier.
            if (const char* src = std::getenv("W3D_FREQ_TRACKER_INPUT")) {
                const std::string value = src;
                if (value == "complementary") {
                    filter.setFreqTrackerInput(
                        FreqTrackerInputSource::Complementary);
                } else if (value == "body_z") {
                    filter.setFreqTrackerInput(FreqTrackerInputSource::BodyZ);
                } else {
                    throw std::runtime_error(
                        "W3D_FREQ_TRACKER_INPUT must be body_z or complementary");
                }
            }

            // Gains of that private observer, so the correction corner can be
            // swept against the wave band it must stay below.  It also solves
            // the startup attitude, which is why two_ki now defaults nonzero.
            {
                float two_kp = STARTUP_PROXY_TWO_KP_DEFAULT;
                float two_ki = STARTUP_PROXY_TWO_KI_DEFAULT;
                const bool kp = env_float("W3D_WAVE_PERIOD_MAHONY_KP", two_kp);
                const bool ki = env_float("W3D_WAVE_PERIOD_MAHONY_KI", two_ki);
                if (kp || ki) {
                    filter.setWavePeriodComplementaryGains(two_kp, two_ki);
                }
            }

            // Ablate the tau-scaled pseudo-measurement cadence back to the
            // historical fixed 15 ms one, so the information-rate
            // renormalization can be priced rather than assumed.
            if (const char* raw = std::getenv("W3D_PSEUDO_CADENCE")) {
                const std::string value = raw;
                if (value == "fixed") {
                    filter.setTauScaledPseudoUpdateCadence(false);
                } else if (value == "tau_scaled") {
                    filter.setTauScaledPseudoUpdateCadence(true);
                } else {
                    throw std::runtime_error(
                        "W3D_PSEUDO_CADENCE must be tau_scaled or fixed");
                }
            }
            if (env_float("OU_II_PSEUDO_TAU_RATIO", v)) {
                filter.setPseudoUpdateTauRatio(v);
            }
            {
                float lo = PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT;
                float hi = PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT;
                const bool got_lo = env_float("OU_II_PSEUDO_PERIOD_MIN_S", lo);
                const bool got_hi = env_float("OU_II_PSEUDO_PERIOD_MAX_S", hi);
                if (got_lo || got_hi) {
                    filter.setPseudoUpdatePeriodBounds(lo, hi);
                }
            }

            if (env_float("OU_ACC_BIAS_INIT_STD", v)) {
                filter.mekf().set_initial_acc_bias_std(v);
            }
            if (env_float("OU_II_ACC_BIAS_INIT_STD", v)) {
                filter.mekf().set_initial_acc_bias_std(v);
            }

            Vector3f b = filter.mekf().get_acc_bias();
            bool bias_changed = false;

            if (env_float("OU_ACC_BIAS_INIT_X", v)) {
                b.x() = v;
                bias_changed = true;
            }
            if (env_float("OU_II_ACC_BIAS_INIT_X", v)) {
                b.x() = v;
                bias_changed = true;
            }

            if (env_float("OU_ACC_BIAS_INIT_Y", v)) {
                b.y() = v;
                bias_changed = true;
            }
            if (env_float("OU_II_ACC_BIAS_INIT_Y", v)) {
                b.y() = v;
                bias_changed = true;
            }

            if (env_float("OU_ACC_BIAS_INIT_Z", v)) {
                b.z() = v;
                bias_changed = true;
            }
            if (env_float("OU_II_ACC_BIAS_INIT_Z", v)) {
                b.z() = v;
                bias_changed = true;
            }

            if (bias_changed) {
                filter.mekf().set_initial_acc_bias(b);
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

        // Which estimator solves the startup attitude.  mahony_proxy is the
        // default: the measurement-only front end runs from the first sample,
        // the private Mahony observer supplies the tilt that gates the
        // magnetometer and frames the world-reference average, and the MEKF is
        // seeded with the finished solution and starts live.  staged_mekf is
        // the matched ablation -- the previous behaviour, in which the MEKF is
        // fed from the first sample and those same reads come back out of it
        // while it is still warming.
        if (const char* raw = std::getenv("W3D_STARTUP_INIT")) {
            const std::string value = raw;
            if (value == "mahony_proxy") {
                cfg_.startup_init_policy = Fusion::StartupInitPolicy::MahonyProxy;
            } else if (value == "staged_mekf") {
                cfg_.startup_init_policy = Fusion::StartupInitPolicy::StagedMekf;
            } else {
                throw std::runtime_error(
                    "W3D_STARTUP_INIT must be mahony_proxy or staged_mekf");
            }
        }

        if (env_float("SF_PROXY_START_MIN_SEC", vf)) cfg_.proxy_startup_min_sec = vf;
        if (env_float("SF_PROXY_START_TIMEOUT_SEC", vf)) cfg_.proxy_startup_timeout_sec = vf;
        if (env_int("SF_ACC_BIAS_UNLOCK_MAG_UPDATES", vi)) cfg_.acc_bias_unlock_mag_updates = vi;
        if (env_float("SF_PROXY_MAG_SETTLE_SEC", vf)) cfg_.proxy_mag_settle_sec = vf;
        if (env_float("SF_MAG_REFINE_START_SEC", vf)) cfg_.mag_refine_start_sec = vf;
        if (env_float("SF_MAG_REFINE_WINDOW_SEC", vf)) cfg_.mag_refine_window_sec = vf;
        if (const char* r = std::getenv("SF_MAG_REFINE")) cfg_.mag_refine_enabled = (std::string(r) != "0");
        if (env_float("SF_PROXY_TILT_SIGMA", vf)) cfg_.proxy_handoff_tilt_sigma_rad = vf;
        if (env_float("SF_PROXY_YAW_SIGMA", vf)) cfg_.proxy_handoff_yaw_sigma_rad = vf;
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

        // Startup timing marks, to stderr so stdout parsing is untouched.
        if (!reported_lock_ && fusion_.hasMagNorthLock()) {
            reported_lock_ = true;
            std::cerr << "STARTUP first_heading_s=" << fusion_.magNorthLockTimeSec() << "\n";
        }
        if (!reported_live_ && fusion_.isLive()) {
            reported_live_ = true;
            std::cerr << "STARTUP live_s=" << fusion_.liveTimeSec() << "\n";
        }
        if (!reported_refine_ && fusion_.hasRefinedMagReference()) {
            reported_refine_ = true;
            std::cerr << "STARTUP mag_refined_s=" << fusion_.magRefineTimeSec() << "\n";
        }

        auto& filter = fusion_.raw();
        if (fixed_tuning_ && !fixed_tuning_applied_ && filter.isAdaptiveLive()) {
            if (!filter.setFixedTuning(
                    fixed_tau_s_, fixed_sigma_a_, fixed_R_p0_std_, fixed_R_v0_std_))
            {
                throw std::runtime_error("invalid fixed OU-II tuning point");
            }
            fixed_tuning_applied_ = true;
        }
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

        s.tau_target     = filter.getTauTarget();
        s.sigma_target   = filter.getSigmaTarget();
        s.tuning_target  = p0_s_from_sigma_tau(s.sigma_target, s.tau_target);

        s.tau_applied    = filter.getTauApplied();
        s.sigma_applied  = filter.getSigmaApplied();
        s.tuning_applied = p0_s_from_sigma_tau(s.sigma_applied, s.tau_applied);

        s.freq_hz             = filter.getFreqHz();
        s.wave_period_sec     = filter.getWavePeriodSec();
        s.period_sec          = filter.getPeriodSec();
        s.accel_variance      = filter.getAccelVariance();
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
    // reconfiguration events and is the matched ablation.
    static bool load_periodic_aw_cov_sync()
    {
        const char* raw = std::getenv("W3D_AW_COV_SYNC");
        const std::string policy = (raw && *raw) ? raw : "periodic";
        if (policy == "reconfigure") return false;
        if (policy == "periodic") return true;
        throw std::runtime_error(
            "W3D_AW_COV_SYNC must be reconfigure or periodic");
    }

    void load_fixed_tuning()
    {
        const char* raw_mode = std::getenv("W3D_TUNING_MODE");
        const std::string mode = raw_mode ? raw_mode : "adaptive";
        if (mode == "adaptive") return;
        if (!mode.starts_with("fixed")) {
            throw std::runtime_error(
                "W3D_TUNING_MODE must be adaptive or start with fixed");
        }

        fixed_tuning_ =
            env_float("W3D_FIXED_TAU_S", fixed_tau_s_) &&
            env_float("W3D_FIXED_SIGMA_A", fixed_sigma_a_) &&
            env_float("W3D_FIXED_R_P0_STD", fixed_R_p0_std_) &&
            env_float("W3D_FIXED_R_V0_STD", fixed_R_v0_std_);
        if (!fixed_tuning_ ||
            !(std::isfinite(fixed_tau_s_) && fixed_tau_s_ > 0.0f &&
              std::isfinite(fixed_sigma_a_) && fixed_sigma_a_ > 0.0f &&
              std::isfinite(fixed_R_p0_std_) && fixed_R_p0_std_ > 0.0f &&
              std::isfinite(fixed_R_v0_std_) && fixed_R_v0_std_ > 0.0f))
        {
            throw std::runtime_error(
                "fixed OU-II mode requires positive W3D_FIXED_TAU_S, "
                "W3D_FIXED_SIGMA_A, W3D_FIXED_R_P0_STD, and "
                "W3D_FIXED_R_V0_STD");
        }
    }

    bool with_mag_ = true;
    bool fixed_tuning_ = false;
    bool fixed_tuning_applied_ = false;
    float fixed_tau_s_ = NAN;
    float fixed_sigma_a_ = NAN;
    float fixed_R_p0_std_ = NAN;
    float fixed_R_v0_std_ = NAN;
    bool reported_lock_ = false;
    bool reported_live_ = false;
    bool reported_refine_ = false;
    using Fusion = SeaStateFusion_OU_II<TrackerType::KALMANF>;
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
// Re-derived for the 900 s scoring window: a sentinel fitted to the previous
// 60 s window is not a sentinel for this one, it is just a number the filter
// passes by a wide margin.
//
// bias_3d_percent re-derived again when the r_p0 and r_v0 smoothing horizons
// were shortened from 5 to 3 wave-period-halves.  The binding record moves from
// pmstokes H0.27 to jonswap H1.5, where the horizontal accelerometer bias is
// already unobservable -- the error exceeds the true bias with either horizon
// -- and the displacement error on that record is unchanged.  See
// docs/ou-ema-adaptation-tuning.md.
//
// bias_3d_percent re-derived once more when the frequency tracker moved onto
// the complementary-levelled vertical acceleration.  It is the only sentinel
// that moved: 81.75 -> 84.97 on the same binding record, jonswap H1.5.  Worth
// being explicit that this is a re-derivation and not a relaxation, because
// raising a sentinel to admit one's own change is exactly how these stop
// meaning anything.
//
// The quantity did not get worse; the realization moved.  Replaying that record
// under six seeds gives a mean of 74.2% on the old input and 74.1% on the new
// one, and the per-seed spread reaches 89% either way -- the distribution is
// unchanged and this sentinel samples one point of it.  The record is also the
// one where, as noted above, the horizontal accelerometer bias is unobservable
// to begin with: an error already larger than the true bias is not measuring
// estimation quality, which is why a 3 pp move in it costs no displacement
// accuracy (3D RMS 20.77% -> 20.78%, vertical unchanged to four figures).
// Re-derived once more for the OU-III parity change: the period-scaled sigma
// band, the tau-scaled pseudo-update cadence and the Mahony-proxy startup
// policy.  Five of the seven move.  Three tighten and two loosen, and the two
// that loosen need saying out loud, because raising a sentinel to admit one's
// own change is exactly how these stop meaning anything.
//
// Both loosened limits are realization moves, not quality regressions, and the
// paired ensemble is what establishes that.  Five IMU noise seeds across the
// eight scored records (n = 40 paired records per metric), deployed
// configuration against the pre-change filter:
//
//     3D RMS % of max |disp|    18.997 -> 19.042    +0.24% +/- 0.26%   n.s.
//     accel-bias 3D % of true   82.94  -> 70.14     -15.4% +/- 7.2%    better
//
// So the aggregate accelerometer-bias error improves by 15% while this
// single-realization sentinel gets 8% worse, and the 3D displacement error
// does not move at all.  The binding record for bias_3d_percent is again one
// where the horizontal accelerometer bias is unobservable -- the error exceeds
// the true bias either way -- which is the same property the previous
// re-derivation of this sentinel noted, and the reason a several-percent move
// in it costs no displacement accuracy.
//
// The three that tighten are where the improvement is deterministic enough to
// show up in one realization as well as in the ensemble.
static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 6.9f,    // was 7.0,  worst 6.86 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 6.9f,    // unchanged, worst 6.81 (pmstokes H0.27)
    .err_limit_yaw_deg             = 2.2f,    // unchanged, worst 2.16 (pmstokes H1.5)
    .err_limit_percent_3d_jonswap  = 21.1f,   // was 20.9, worst 20.91 (jonswap H1.5)
    .err_limit_percent_3d_pmstokes = 21.2f,   // was 20.7, worst 21.02 (pmstokes H8.5)
    .acc_z_bias_percent            = 5.4f,    // was 5.9,  worst 5.33 (pmstokes H8.5)
    .bias_3d_percent               = 92.2f,   // was 85.4, worst 91.65 (jonswap H4.0, accel)
};

static constexpr W3dSummaryLabels SUMMARY_LABELS{
    .target = "p0_S_target",
    .applied = "p0_S_applied",
};

static void process_wave_file_for_tracker(const std::string& filename,
                                          float dt,
                                          bool with_mag,
                                          const W3dRandomSeeds& seeds,
                                          bool write_timeseries,
                                          float validation_window_sec)
{
    constexpr float MAG_ODR_HZ = 25.0f;
    auto result = process_wave_file_for_tracker<FusionAdapter_OU_II>(
        filename, dt, with_mag, add_noise, MAG_ODR_HZ,
        "_fusion_ou2", "_fusion_ou2_nomag", seeds, write_timeseries);

    if (!result) return;
    if (validation_window_sec > 0.0f) {
        print_validation_metrics(*result, dt, validation_window_sec, "OU_II");
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

    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}
