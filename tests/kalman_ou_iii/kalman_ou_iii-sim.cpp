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
    // RAO replay tuning. Gyro sample standard deviation is converted to a
    // continuous-time density at 200 Hz. Magnetometer uncertainty includes
    // residual calibration/model error, not just injected white noise.
    // Paired training and independent validation are recorded in
    // reports/results/rao_parameter_tuning; quality limits are unchanged.
    // Environment scale factors below multiply these deployed settings.
    static constexpr float SIGMA_A_RESCALE = 0.71f;  // 2.8x -> 2.0x injected accel white
    static constexpr float SIGMA_G_RESCALE = 0.01f;  // empirical gyro weighting on vessel replay
    static constexpr float SIGMA_M_RESCALE = 8.0f;   // 1.2x -> 9.6x injected mag white

    FusionAdapter_OU_III(bool with_mag,
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


        const std::string aw_cov_sync = load_aw_cov_sync_policy();
        filter.setPeriodicAwCovarianceSync(aw_cov_sync != "reconfigure");
        filter.setAwCovarianceSyncCongruent(aw_cov_sync == "congruent");

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
            // X and Y are independent; there is deliberately no combined knob,
            // so a sweep that means to move both has to say so twice.
            if (env_float("OU_III_R_S_X_FACTOR", v)) {
                filter.setRSXFactor(v);
            }
            if (env_float("OU_III_R_S_Y_FACTOR", v)) {
                filter.setRSYFactor(v);
            }

            if (env_float("OU_III_R_S_COEFF", v)) {
                filter.setRSCoeff(v);
            }

            // Integral-regularizer adaptation law ablation.
            // 0 = Cubic (deployed), 1 = StrongRiccati, 2 = PosteriorRiccati,
            // 3 = SpectralMSE (bias-variance).
            if (env_float("OU_III_RS_LAW", v)) {
                const int law = static_cast<int>(v);
                if (law == 1) {
                    filter.setRSLaw(RSAdaptationLaw::StrongRiccati);
                } else if (law == 2) {
                    filter.setRSLaw(RSAdaptationLaw::PosteriorRiccati);
                } else if (law == 3) {
                    filter.setRSLaw(RSAdaptationLaw::SpectralMSE);
                } else {
                    filter.setRSLaw(RSAdaptationLaw::Cubic);
                }
            }
            if (env_float("OU_III_RS_KAPPA", v)) {
                filter.setRSPoleKappa(v);
            }
            if (env_float("OU_III_RS_RA", v)) {
                filter.setRSAccelNoiseDensity(v);
            }
            // C_J of the SpectralMSE (bias-variance) law.
            if (env_float("OU_III_RS_MSE_COEFF", v)) {
                filter.setRSMseCoeff(v);
            }
            if (env_float("OU_III_RS_SIGMA_EXP", v)) {
                filter.setRSSigmaExponent(v);
            }

            // r_S safety clamp, in m*s.  The floor binds in low-motion seas,
            // where the schedule asks for less than the default 0.4, so it is
            // a real tuning surface rather than a formality.
            {
                float lo = MIN_R_S, hi = MAX_R_S;
                const bool got_lo = env_float("OU_III_R_S_MIN", lo);
                const bool got_hi = env_float("OU_III_R_S_MAX", hi);
                if (got_lo || got_hi) {
                    filter.setRSBounds(lo, hi);
                }
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

            // Out-of-band accelerometer guard ahead of the proxy and the
            // MEKF.  Absent, the measurement path is unconditioned, which is
            // the deployed default.
            {
                float guard_hz = 0.0f;
                int guard_poles = 2;
                env_int("OU_III_ACC_GUARD_POLES", guard_poles);
                if (env_float("OU_III_ACC_GUARD_HZ", guard_hz)) {
                    filter.setAccelVibrationGuard(guard_hz, guard_poles);
                }
                float racc_gain = 0.0f;
                if (env_float("OU_III_ACC_GUARD_RACC_GAIN", racc_gain)) {
                    filter.setAccelVibrationRaccGain(racc_gain);
                }
                float engage_lo = 0.0f, engage_hi = 0.0f, engage_tau = 0.0f;
                const bool lo_set = env_float("OU_III_ACC_GUARD_ENGAGE_LO", engage_lo);
                const bool hi_set = env_float("OU_III_ACC_GUARD_ENGAGE_HI", engage_hi);
                const bool tau_set = env_float("OU_III_ACC_GUARD_ENGAGE_TAU", engage_tau);
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
            if (env_float("OU_III_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }

            // Common tau/sigma EMA horizon in measured sea-time units T_sea=T_z/2.
            if (env_float("OU_ADAPT_TAU_SEA_PERIODS", v)) {
                filter.setAdaptationSeaPeriods(v);
            }
            if (env_float("OU_III_ADAPT_TAU_SEA_PERIODS", v)) {
                filter.setAdaptationSeaPeriods(v);
            }

            // Legacy tuner-frequency smoothing knobs, retained so existing
            // ablation scripts keep running.  SeaStateAutoTuner no longer
            // smooths frequency -- WavePeriodEstimator owns the canonical
            // log-period state -- so these setters have no effect.
            if (env_float("OU_TUNER_FREQ_SEA_PERIODS", v)) {
                filter.setTunerFreqSmoothingSeaPeriods(v);
            }
            if (env_float("OU_III_TUNER_FREQ_SEA_PERIODS", v)) {
                filter.setTunerFreqSmoothingSeaPeriods(v);
            }
            if (env_float("OU_TUNER_FREQ_TAU_SEC", v)) {
                filter.setTunerFreqSmoothingTimeConstant(v);
            }
            if (env_float("OU_III_TUNER_FREQ_TAU_SEC", v)) {
                filter.setTunerFreqSmoothingTimeConstant(v);
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
            if (env_float("OU_III_ACC_BIAS_TAU_SEC", v)) {
                filter.mekf().set_acc_bias_time_constant(v);
            }

            // Wave-band prior used until the period estimator has a value.
            if (env_float("OU_TUNE_FREQ_PRIOR_HZ", v)) {
                filter.setTuneFreqPriorHz(v);
            }

            // sigma_a averaging horizon, in periods of the tuning frequency,
            // and its absolute clamps in seconds.
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
            // complementary-levelled default.  "leveled" restores the older
            // behaviour, which levels with the attitude solution and so closes
            // the tuner coupling.
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
        if (env_float("SF_MAG_GRAV_ALIGN_LPF_TAU", vf)) cfg_.mag_gravity_align_world_tau_sec = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_WARMUP_SEC", vf)) cfg_.mag_gravity_align_world_warmup_sec = vf;
        if (env_float("SF_MAG_TILT_FALLBACK_SEC", vf)) cfg_.mag_tilt_fallback_sec = vf;
        if (env_float("SF_MAG_EXTREME_GYRO_DPS", vf)) cfg_.mag_extreme_gyro_dps = vf;
        if (env_float("SF_MAG_INIT_MIN_MAG_NORM", vf)) cfg_.mag_init_min_mag_norm = vf;
        if (env_int("SF_MAG_MIN_SAMPLES", vi)) cfg_.mag_min_samples = vi;
        if (env_float("SF_MAG_MIN_WINDOW_SEC", vf)) cfg_.mag_min_window_sec = vf;

        if (env_float("SF_ONLINE_TUNE_WARMUP_SEC", vf)) cfg_.online_tune_warmup_sec = vf;

        // The MEKF variances the Kalman3D_Wave_OU_III constructor takes.
        //
        // The three sensor sigmas arrive here already built by the shared
        // harness as a multiple of the noise it actually injects -- 2.8x on
        // accel, 2.0x on gyro, 1.2x on mag -- so they are swept as scale
        // factors on that inflation rather than as absolute values.  A scale
        // of 1 therefore leaves the deployed constant in place, and the value
        // the sweep reports is directly the inflation multiplier it prefers.
        // The remaining four are absolute and go through the Config fields
        // added for them.
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
        if (env_float("SF_RS_NOISE_VAR", vf)) cfg_.R_S_noise = vf;

        if (env_float("SF_PROXY_START_MIN_SEC", vf)) cfg_.proxy_startup_min_sec = vf;
        if (env_float("SF_PROXY_START_TIMEOUT_SEC", vf)) cfg_.proxy_startup_timeout_sec = vf;
        if (env_int("SF_ACC_BIAS_UNLOCK_MAG_UPDATES", vi)) cfg_.acc_bias_unlock_mag_updates = vi;
        if (env_float("SF_PROXY_MAG_SETTLE_SEC", vf)) cfg_.proxy_mag_settle_sec = vf;
        if (env_float("SF_MAG_REFINE_START_SEC", vf)) cfg_.mag_refine_start_sec = vf;
        if (env_float("SF_MAG_REFINE_WINDOW_SEC", vf)) cfg_.mag_refine_window_sec = vf;
        if (const char* r = std::getenv("SF_MAG_REFINE")) cfg_.mag_refine_enabled = (std::string(r) != "0");
        // Continuous exogenous hard-iron estimation.
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
        }
        if (!reported_refine_ && fusion_.hasRefinedMagReference()) {
            reported_refine_ = true;
            std::cerr << "STARTUP mag_refined_s=" << fusion_.magRefineTimeSec() << "\n";
        }

        if (trace_hard_iron_) {
            trace_elapsed_ += dt;
            if (trace_elapsed_ >= 60.0f) {
                trace_elapsed_ = 0.0f;
                const auto& e = fusion_.magContinuousHardIron().estimate();
                const Vector3f& a = fusion_.magContinuousHardIronAppliedUT();
                std::cerr << "MAGHI t=" << int(t_trace_)
                          << " valid=" << e.valid
                          << " info=" << e.information
                          << " w=" << e.effective_weight
                          << " resid=" << e.residual_rms_uT
                          << " fit=[" << e.bias_body_uT.transpose() << "]"
                          << " applied=[" << a.transpose() << "]\n";
            }
            t_trace_ += dt;
        }

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

    // Vibration-guard telemetry, so a replay can say whether the guard engaged
    // and how much out-of-band accelerometer content it was seeing.  The shared
    // runner owns the adapter, so end-of-record is the destructor; silent
    // unless a cutoff was configured, which keeps default runs unchanged.
    ~FusionAdapter_OU_III() override {
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
    const bool trace_hard_iron_ = std::getenv("W3D_MAG_HI_TRACE") != nullptr;
    mutable float trace_elapsed_ = 0.0f;
    mutable float t_trace_ = 0.0f;

    mutable bool reported_lock_ = false;
    mutable bool reported_live_ = false;
    mutable bool reported_refine_ = false;
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
// scored records plus about half a percent, rounded up in the last digit the
// channel is quoted in
static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 4.271f,  // worst 4.2494 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 4.235f,  // worst 4.2130 (pmstokes H0.27)
    .err_limit_yaw_deg             = 0.8187f, // worst 0.8145 (jonswap H1.5)
    .err_limit_roll_deg            = 0.3292f, // worst 0.3275 (pmstokes H4.0)
    .err_limit_pitch_deg           = 0.1473f, // worst 0.1465 (pmstokes H4.0)
    .err_limit_percent_3d_jonswap  = 9.741f,  // worst 9.6917 (jonswap H8.5)
    .err_limit_percent_3d_pmstokes = 9.706f,  // worst 9.6568 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.485f,  // worst 4.4617 (pmstokes H8.5)
    .bias_3d_percent               = 69.85f,  // worst 69.4941 (pmstokes H4.0, accel)
    .gyro_bias_3d_percent          = 12.27f,  // worst 12.2087 (pmstokes H0.27, gyro)
};

static constexpr W3dSummaryLabels SUMMARY_LABELS {
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
