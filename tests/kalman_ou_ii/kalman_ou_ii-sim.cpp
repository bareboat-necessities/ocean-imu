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

namespace {

bool env_float(const char* name, float& out) {
    if (const char* s = std::getenv(name)) {
        out = static_cast<float>(std::atof(s));
        return true;
    }
    return false;
}

bool env_int(const char* name, int& out) {
    if (const char* s = std::getenv(name)) {
        out = std::atoi(s);
        return true;
    }
    return false;
}

} // namespace

class FusionAdapter_OU_II final : public IW3dFusionAdapter {
public:
    // RAO replay tuning, validated on separate sensor/initialization draws.
    // Gyro sample standard deviation is converted to a density at 200 Hz.
    // Magnetometer uncertainty includes residual calibration/model error.
    // Environment scales multiply these deployed settings; gates are fixed.
    static constexpr float SIGMA_A_RESCALE = 0.5f;   // 2.8x -> 1.4x injected accel white
    static constexpr float SIGMA_G_RESCALE = 0.05f;  // 2.0x sample std -> sqrt(2)x density
    static constexpr float SIGMA_M_RESCALE = 8.0f;   // 1.2x -> 9.6x injected mag white

    FusionAdapter_OU_II(bool with_mag,
                        const Vector3f& sigma_a_init,
                        const Vector3f& sigma_g,
                        const Vector3f& sigma_m)
        : with_mag_(with_mag)
    {
        cfg_.with_mag = with_mag;
        cfg_.sigma_a = sigma_a_init * SIGMA_A_RESCALE;
        cfg_.sigma_g = sigma_g * SIGMA_G_RESCALE;
        cfg_.sigma_m = sigma_m * SIGMA_M_RESCALE;
        cfg_.mag_delay_sec = MAG_DELAY_SEC;
        apply_env_overrides();
        load_fixed_tuning();

        fusion_.begin(cfg_);
        auto& filter = fusion_.raw();

        filter.setPeriodicAwCovarianceSync(load_periodic_aw_cov_sync());

        {
            // Known stationary RAO of the pinned simulation dataset.
            // ENU->NED swaps horizontal axes: the direction branch's X reads
            // vessel sway, and its Y reads vessel surge. Follow the actual
            // basis conversion, not the legacy forward/starboard labels.
            wave_direction::VesselRaoEqualizer::Config direction_rao;
            direction_rao.enabled = true;
            direction_rao.horizontal_x_tau_s = 1.0;
            direction_rao.horizontal_y_tau_s = 0.7;
            if (const char* mode = std::getenv("W3D_DIRECTION_RAO")) {
                const std::string name(mode);
                if (name == "off") direction_rao.enabled = false;
                else if (name != "vessel-rao-28ft")
                    throw std::invalid_argument("Unknown W3D_DIRECTION_RAO profile");
            }
            filter.setDirectionRao(direction_rao);
            wave_direction::VesselRaoNoiseWeighting low_wave_noise;
            low_wave_noise.response = direction_rao;
            // The known hull attenuates horizontal wave motion toward the
            // accelerometer noise floor. Revert to nominal weights above SNR 4.
            low_wave_noise.max_std_scale = 4.0f;
            env_float("SF_LOW_WAVE_RACC_MAX_SCALE", low_wave_noise.max_std_scale);
            env_float("SF_LOW_WAVE_RACC_SNR", low_wave_noise.transition_snr);
            filter.setLowWaveNoiseWeighting(low_wave_noise);

            filter.enableTuner(true);
            filter.enableClamp(true);

            float v = 0.0f;
            if (env_float("OU_P_FACTOR", v)) {
                filter.setPFactor(v);
            }
            if (env_float("OU_II_P_FACTOR", v)) {
                filter.setPFactor(v);
            }
            if (env_float("OU_R_P0_X_FACTOR", v)) {
                filter.setR_p0_XFactor(v);
            }
            if (env_float("OU_II_R_P0_X_FACTOR", v)) {
                filter.setR_p0_XFactor(v);
            }
            if (env_float("OU_R_P0_Y_FACTOR", v)) {
                filter.setR_p0_YFactor(v);
            }
            if (env_float("OU_II_R_P0_Y_FACTOR", v)) {
                filter.setR_p0_YFactor(v);
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

            // Dual pseudo-measurement adaptation law.
            // 0 = Empirical (c_p sigma_aw tau^2, c_v sigma_aw tau),
            // 1 = PhysicalMSE (deployed joint displacement-MSE law).
            if (env_float("OU_II_PSEUDO_LAW", v)) {
                const int law = static_cast<int>(v);
                if (law == 0) {
                    filter.setPseudoLaw(PseudoAdaptationLaw::Empirical);
                } else {
                    filter.setPseudoLaw(PseudoAdaptationLaw::PhysicalMSE);
                }
            }
            // C_P and the channel ratio C_P/C_V of the PhysicalMSE law.
            if (env_float("OU_II_PSEUDO_MSE_COEFF", v)) {
                filter.setPseudoMseCoeff(v);
            }
            if (env_float("OU_II_PSEUDO_MSE_RATIO", v)) {
                filter.setPseudoMseRatio(v);
            }
            if (env_float("OU_II_PSEUDO_RA", v)) {
                filter.setPseudoAccelNoiseDensity(v);
            }

            // Clamps on the two drift-band regularizers.
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

            // Out-of-band accelerometer guard ahead of the proxy and the
            // MEKF.  Armed by default at the deployed corner; these override
            // it, and a zero cutoff removes it entirely.
            {
                float guard_hz = 0.0f;
                int guard_poles = ACC_VIBRATION_GUARD_POLES_DEFAULT;
                env_int("OU_II_ACC_GUARD_POLES", guard_poles);
                if (env_float("OU_II_ACC_GUARD_HZ", guard_hz)) {
                    filter.setAccelVibrationGuard(guard_hz, guard_poles);
                }
                float racc_gain = 0.0f;
                if (env_float("OU_II_ACC_GUARD_RACC_GAIN", racc_gain)) {
                    filter.setAccelVibrationRaccGain(racc_gain);
                }
                float engage_lo = 0.0f, engage_hi = 0.0f, engage_tau = 0.0f;
                const bool lo_set = env_float("OU_II_ACC_GUARD_ENGAGE_LO", engage_lo);
                const bool hi_set = env_float("OU_II_ACC_GUARD_ENGAGE_HI", engage_hi);
                const bool tau_set = env_float("OU_II_ACC_GUARD_ENGAGE_TAU", engage_tau);
                if (lo_set || hi_set || tau_set) {
                    filter.setAccelVibrationEngagement(
                        lo_set ? engage_lo : -1.0f,
                        hi_set ? engage_hi : -1.0f,
                        tau_set ? engage_tau : -1.0f);
                }
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

            // Accelerometer-bias random walk.
            Eigen::Vector3f bias_rw(0.00015f, 0.00015f, 0.0004f);
            if (env_float("OU_II_ACC_BIAS_RW", v)) bias_rw.setConstant(v);
            env_float("OU_II_ACC_BIAS_RW_X", bias_rw.x());
            env_float("OU_II_ACC_BIAS_RW_Y", bias_rw.y());
            env_float("OU_II_ACC_BIAS_RW_Z", bias_rw.z());
            filter.mekf().set_Q_bacc_rw(bias_rw);
            float bias_tau_sec = 20000.0f;
            env_float("OU_II_ACC_BIAS_TAU_SEC", bias_tau_sec);
            filter.mekf().set_acc_bias_time_constant(bias_tau_sec);

            // Wave-band prior used until the period estimator has a value.
            if (env_float("OU_TUNE_FREQ_PRIOR_HZ", v)) {
                filter.setTuneFreqPriorHz(v);
            }

            // sigma_a averaging horizon, in periods of the tuning frequency,
            // and its absolute clamps in seconds.
            if (env_float("OU_SIGMA_STILL_DECAY_SEC", v)) {
                filter.setSigmaStillnessDecaySec(v);
            }
            if (env_float("OU_SIGMA_VAR_K_PERIODS", v)) {
                filter.setSigmaVarianceKPeriods(v);
            }
            {
                float lo = 0.3f, hi = 60.0f;
                const bool got_lo = env_float("OU_SIGMA_VAR_HORIZON_MIN_S", lo);
                const bool got_hi = env_float("OU_SIGMA_VAR_HORIZON_MAX_S", hi);
                if (got_lo || got_hi) {
                    filter.setSigmaVarianceHorizonBounds(lo, hi);
                }
            }

            // Ablate the wave-period estimator's input away from the
            // complementary-levelled default.
            if (const char* src = std::getenv("W3D_WAVE_PERIOD_INPUT")) {
                const std::string value = src;
                if (value == "complementary") {
                    filter.setWavePeriodInput(
                        WavePeriodInputSource::Complementary);
                } else if (value == "leveled") {
                    filter.setWavePeriodInput(WavePeriodInputSource::Leveled);
                } else {
                    throw std::runtime_error(
                        "W3D_WAVE_PERIOD_INPUT must be leveled or "
                        "complementary");
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
        if (env_float("SF_MAG_GRAV_ALIGN_LPF_TAU", vf)) cfg_.mag_gravity_align_world_tau_sec = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_WARMUP_SEC", vf)) cfg_.mag_gravity_align_world_warmup_sec = vf;
        if (env_float("SF_MAG_TILT_FALLBACK_SEC", vf)) cfg_.mag_tilt_fallback_sec = vf;
        if (env_float("SF_MAG_EXTREME_GYRO_DPS", vf)) cfg_.mag_extreme_gyro_dps = vf;
        if (env_float("SF_MAG_INIT_MIN_MAG_NORM", vf)) cfg_.mag_init_min_mag_norm = vf;
        if (env_int("SF_MAG_MIN_SAMPLES", vi)) cfg_.mag_min_samples = vi;
        if (env_float("SF_MAG_MIN_WINDOW_SEC", vf)) cfg_.mag_min_window_sec = vf;

        if (env_float("SF_ONLINE_TUNE_WARMUP_SEC", vf)) cfg_.online_tune_warmup_sec = vf;

        // The MEKF variances the Kalman3D_Wave_OU_II constructor takes; see
        // the SIGMA_*_RESCALE block above.
        // The three sensor sigmas are swept as scale factors on the deployed
        // point, so a scale of 1 leaves it in place.
        if (env_float("SF_SIGMA_A_SCALE", vf)) cfg_.sigma_a *= vf;
        // Body-axis measurement weighting; injected sensor noise is unchanged.
        if (env_float("SF_SIGMA_A_X_SCALE", vf)) cfg_.sigma_a.x() *= vf;
        if (env_float("SF_SIGMA_A_Y_SCALE", vf)) cfg_.sigma_a.y() *= vf;
        if (env_float("SF_SIGMA_A_Z_SCALE", vf)) cfg_.sigma_a.z() *= vf;
        if (env_float("SF_SIGMA_G_SCALE", vf)) cfg_.sigma_g *= vf;
        if (env_float("SF_SIGMA_M_SCALE", vf)) cfg_.sigma_m *= vf;

        if (env_float("SF_PQ0", vf)) cfg_.Pq0 = vf;
        if (env_float("SF_PB0", vf)) cfg_.Pb0 = vf;
        if (env_float("SF_GYRO_BIAS_RW_VAR", vf)) cfg_.b0 = vf;
        if (env_float("SF_RP0_NOISE_VAR", vf)) cfg_.R_p0_noise = vf;
        if (env_float("SF_RV0_NOISE_VAR", vf)) cfg_.R_v0_noise = vf;

        if (env_float("SF_PROXY_START_MIN_SEC", vf)) cfg_.proxy_startup_min_sec = vf;
        if (env_float("SF_PROXY_START_TIMEOUT_SEC", vf)) cfg_.proxy_startup_timeout_sec = vf;
        if (env_int("SF_ACC_BIAS_UNLOCK_MAG_UPDATES", vi)) cfg_.acc_bias_unlock_mag_updates = vi;
        if (env_float("SF_PROXY_MAG_SETTLE_SEC", vf)) cfg_.proxy_mag_settle_sec = vf;
        if (env_float("SF_MAG_REFINE_START_SEC", vf)) cfg_.mag_refine_start_sec = vf;
        if (env_float("SF_MAG_REFINE_WINDOW_SEC", vf)) cfg_.mag_refine_window_sec = vf;
        if (const char* r = std::getenv("SF_MAG_REFINE")) cfg_.mag_refine_enabled = (std::string(r) != "0");
        if (const char* r = std::getenv("SF_STARTUP_BIAS_SEED")) cfg_.startup_vertical_bias_seed = (std::string(r) != "0");
        if (env_float("SF_STARTUP_BIAS_SEED_R", vf)) cfg_.startup_vertical_bias_seed_pole_rate = vf;
        // Continuous exogenous hard-iron estimation.  Same names as OU-III's,
        // so a paired study can set one environment and run both families.
        if (const char* h = std::getenv("SF_MAG_CONT_HI")) cfg_.mag_continuous_hard_iron = (std::string(h) != "0");
        if (env_float("SF_MAG_HI_MEMORY_SEC", vf)) cfg_.mag_hi_memory_sec = vf;
        if (env_float("SF_MAG_HI_RIDGE", vf)) cfg_.mag_hi_model_ridge = vf;
        if (env_float("SF_MAG_HI_RIDGE_REL", vf)) cfg_.mag_hi_model_ridge_relative = vf;
        if (env_float("SF_MAG_HI_MIN_INFO", vf)) cfg_.mag_hi_min_information = vf;
        if (env_float("SF_MAG_HI_MIN_WEIGHT", vf)) cfg_.mag_hi_min_effective_weight = vf;
        if (env_float("SF_MAG_HI_MAX_RESID", vf)) cfg_.mag_hi_max_residual_rms_uT = vf;
        if (env_float("SF_MAG_HI_FRACTION", vf)) cfg_.mag_hi_apply_fraction = vf;
        if (env_float("SF_MAG_HI_SLEW_TAU", vf)) cfg_.mag_hi_slew_tau_sec = vf;
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
            std::cerr << "STARTUP vertical_bias_seed_mps2=" << fusion_.startupVerticalBiasSeedMps2() << "\n";
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

    // Vibration-guard telemetry, so a replay can say whether the guard engaged
    // and how much out-of-band accelerometer content it was seeing.  The shared
    // runner owns the adapter, so end-of-record is the destructor; silent
    // unless a cutoff was configured, which keeps an unguarded run unchanged.
    ~FusionAdapter_OU_II() override {
        const auto& filter = fusion_.raw();
        if (!(filter.accelVibrationGuardCutoffHz() > 0.0f)) return;
        std::cout << "ACC_GUARD cutoff_hz=" << filter.accelVibrationGuardCutoffHz()
                  << " poles=" << filter.accelVibrationGuardPoles()
                  << " engagement=" << filter.accelVibrationGuardEngagement()
                  << " out_of_band_rms_mps2=" << filter.accelVibrationRms()
                  << " delay_sec=" << filter.accelVibrationGuardDelaySec()
                  << " racc_gain=" << filter.accelVibrationRaccGain()
                  << " racc_std_mps2=" << filter.accelVibrationRaccStd().x()
                  << "\n";
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

        // The MEKF has no magnetometer-bias state, so get_mag_bias_est_uT()
        // returns zero here and the harness's "Bias error RMS (mag)" was
        // reporting the injected offset itself rather than anything the filter
        // had done about it.  The wrapper's continuous hard-iron correction is
        // the estimate of that offset -- it is subtracted from every sample the
        // MEKF sees -- so it is what belongs in this slot, and reporting it
        // makes the correction measurable in uT instead of only through yaw.
        s.mag_bias_est_ned_uT = get_mag_bias_est_uT(filter.mekf()) +
                                fusion_.magHardIronBodyUT();

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
// scored records plus about half a percent, rounded up in the last digit the
// channel is quoted in
static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 6.674f,  // worst 6.6408 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 6.605f,  // worst 6.5717 (pmstokes H0.27)
    .err_limit_yaw_deg             = 1.041f,  // worst 1.0357 (jonswap H4.0)
    .err_limit_roll_deg            = 0.3866f, // worst 0.3846 (jonswap H4.0)
    .err_limit_pitch_deg           = 0.3237f, // worst 0.3221 (jonswap H8.5)
    .err_limit_percent_3d_jonswap  = 15.54f,  // worst 15.4571 (jonswap H8.5)
    .err_limit_percent_3d_pmstokes = 16.62f,  // worst 16.5347 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.663f,  // worst 4.6391 (jonswap H8.5)
    .bias_3d_percent               = 83.97f,  // worst 83.5458 (jonswap H4.0, accel)
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
