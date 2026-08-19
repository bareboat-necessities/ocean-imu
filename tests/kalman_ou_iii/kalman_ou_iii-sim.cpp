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
    // OU-III's own MEKF sensor variances, as multiples of the ones the shared
    // harness hands every family.  The harness builds each as a fixed multiple
    // of the white noise it injects -- 2.8x on accel, 2.0x on gyro, 1.2x on
    // mag -- and those multiples had never been swept for any family.
    // docs/ou-iii-qmekf-variances.md is the sweep; it moves all three, and the
    // gyro one is much the largest effect.
    //
    // sigma_g is a units correction, not a fit.  The harness multiplies a
    // *per-sample* gyro standard deviation at 200 Hz, but
    // Kalman3D_Wave_OU_III integrates this argument as a noise *density*
    // (Q_AA = Qbase * Ts), so the deployed value overstated the angular random
    // walk by sqrt(200) = 14.1x on top of its own 2x inflation -- 28.3x in
    // std, 800x in variance.  0.05 puts the argument back on the injected
    // density and keeps a sqrt(2) inflation over it, which is where the
    // measured optimum sits: correcting the units *is* the optimum, to within
    // the width of the basin.
    //
    // The other two are empirical, measured against this harness's noise model
    // rather than derived, and both are small next to the gyro term: accel to
    // 2.0x the injected white (from 2.8x) and mag to 2.4x (from 1.2x).
    //
    // Paired over 8 records x 5 seeds against the previous point, every
    // channel improves and none is traded away: roll -13.3%, pitch -14.7%,
    // yaw -1.2%, vertical displacement -2.4%, 3D displacement -37.2%,
    // accelerometer bias -1.0%, gyro bias -36.5%.  The tuner operating point
    // (tau_applied, sigma_applied) is bit-for-bit unchanged, so this is a
    // MEKF-side effect only and does not perturb the OU schedule.
    //
    // The SF_SIGMA_*_SCALE overrides below multiply these, so a scale of 1
    // reproduces the deployed point and a re-run of the sweep re-centres on it.
    static constexpr float SIGMA_A_RESCALE = 0.71f;  // 2.8x -> 2.0x injected accel white
    static constexpr float SIGMA_G_RESCALE = 0.05f;  // 2.0x sample std -> sqrt(2)x density
    static constexpr float SIGMA_M_RESCALE = 2.0f;   // 1.2x -> 2.4x injected mag white

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

            // Where the tuning frequency comes from before the wave-period
            // estimator has a value.  "wave_band" is the deployed source and
            // never reads the acceleration-band tracker; "wave_band_gated"
            // keeps the legacy readiness gate but still never reads it;
            // "tracker_fallback" is the previous behaviour.
            if (const char* src = std::getenv("W3D_TUNER_FREQ_SOURCE")) {
                const std::string value = src;
                if (value == "wave_band") {
                    filter.setTunerFrequencySource(TunerFrequencySource::WaveBand);
                } else if (value == "wave_band_gated") {
                    filter.setTunerFrequencySource(
                        TunerFrequencySource::WaveBandGated);
                } else if (value == "tracker_fallback") {
                    filter.setTunerFrequencySource(
                        TunerFrequencySource::TrackerFallback);
                } else {
                    throw std::runtime_error(
                        "W3D_TUNER_FREQ_SOURCE must be wave_band, "
                        "wave_band_gated or tracker_fallback");
                }
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

        // The MEKF variances the Kalman3D_Wave_OU_III constructor takes.
        //
        // The three sensor sigmas arrive here already built by the shared
        // harness as a multiple of the noise it actually injects -- 2.8x on
        // accel, 2.0x on gyro, 1.2x on mag -- so they are swept as scale
        // factors on that inflation rather than as absolute values.  A scale
        // of 1 therefore leaves the deployed constant in place, and the value
        // the sweep reports is directly the inflation multiplier it prefers.
        // The remaining four are absolute and go through the Config fields
        // added for them.  See docs/ou-iii-qmekf-variances.md.
        if (env_float("SF_SIGMA_A_SCALE", vf)) cfg_.sigma_a *= vf;
        if (env_float("SF_SIGMA_G_SCALE", vf)) cfg_.sigma_g *= vf;
        if (env_float("SF_SIGMA_M_SCALE", vf)) cfg_.sigma_m *= vf;

        if (env_float("SF_PQ0", vf)) cfg_.Pq0 = vf;
        if (env_float("SF_PB0", vf)) cfg_.Pb0 = vf;
        if (env_float("SF_GYRO_BIAS_RW_VAR", vf)) cfg_.b0 = vf;
        if (env_float("SF_RS_NOISE_VAR", vf)) cfg_.R_S_noise = vf;

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
// channel is quoted in -- a tenth for the percentage channels, a hundredth for
// yaw.  Yaw is quoted finer on purpose: at about one degree, a tenth is three
// percent of the value, so rounding it like a percentage channel hands back
// six times the margin the rule asks for.
//
// That margin is deliberately small because the metrics are deterministic, and
// the size of "deterministic" is measured rather than assumed.  Rebuilding the
// same records and seeds at -march=x86-64 instead of the host's native
// cascadelake moves the gated numbers by at most 8.3e-4 relative, the worst of
// them on yaw (jonswap H8.5).  That is two orders of magnitude coarser than
// the 6e-6 this comment used to claim, and the claim was not wrong when it was
// written: the continuous hard-iron solve inverts a normal matrix of order
// 1e-3, which multiplies the last bits of the accumulation by up to a thousand
// on its way into the applied offset.  Half a percent still leaves six times
// the observed build-to-build spread on the tightest gate here, but it is no
// longer an enormous factor, and a future sentinel should be checked against a
// re-measurement rather than against this paragraph.
//
// Setting one below what the filter currently achieves makes it fail every run
// rather than catching a regression.
//
// Re-derived for the 900 s scoring window: a sentinel fitted to the
// previous 60 s window is not a sentinel for this one, it is just a number the
// filter passes by a wide margin.
//
// Re-derived again after two changes that between them moved every gated
// quantity: the r_S adaptation-law work, which took the displacement gates
// down, and the Mahony-proxy startup policy, which took the accelerometer-bias
// gates down.  The previous limits had 5 to 12 percent of slack against the
// filter that now ships, which is slack a regression can hide in.
//
// These are fitted to the deployed configuration -- the default
// StartupInitPolicy::MahonyProxy with the continuous hard-iron correction on.
// Both matched ablations, W3D_STARTUP_INIT=staged_mekf and SF_MAG_CONT_HI=0,
// deliberately exceed some of them, because that is precisely the behaviour
// the defaults replaced.  Scoring an ablation means scoring the old filter, so
// run it with W3D_COLLECT_ALL_GATES if you want the numbers rather than an
// early exit.
//
// bias_3d_percent remains dominated by the horizontal accelerometer bias,
// which is close to unobservable on the smaller seas -- the error exceeds the
// true bias there under any of these configurations.  It moved from 108.9 to
// 95.3 because holding accelerometer-bias learning until the magnetic
// reference is refined stops the bias absorbing the provisional reference's
// tilt error; see docs/ou-iii-startup-init.md.
// Re-derived once more for the continuous hard-iron correction, which moved
// three of the seven.  Yaw is the point of that change and drops by half, so
// its limit comes down with it or it stops being a sentinel.  The two that go
// up are the price: the correction walks the heading onto the corrected field
// during the run, and the horizontal accelerometer bias -- already the least
// observable quantity here, with an error that exceeds the true bias under
// every configuration this filter has ever shipped -- absorbs some of that
// motion.  See docs/continuous-mag-hard-iron.md.
//
// Then cut to the rule rather than to the quantum.  A tenth is 2 percent of a
// 4.7 and 2 percent of a 5.0, so rounding a half-percent margin up to the next
// tenth was handing back four times what the rule asks for on the small-valued
// channels.  Each gate is now written to whatever precision delivers about half
// a percent -- a thousandth for yaw, a hundredth for the single-digit
// percentages, a tenth where a tenth is already fine enough.
//
//   Z %Hs JONSWAP    4.8  -> 4.72    worst 4.6952   2.23% -> 0.53%
//   Z %Hs PM-Stokes  4.7  -> 4.69    worst 4.6600   0.86% -> 0.64%
//   yaw deg          1.07 -> 1.068   worst 1.0627   0.69% -> 0.50%
//   3D % JONSWAP     21.1 -> 21.05   worst 20.9361  0.78% -> 0.55%
//   3D % PM-Stokes   20.9 -> 20.83   worst 20.7197  0.87% -> 0.53%
//   acc Z bias %     5.0  -> 4.93    worst 4.9054   1.93% -> 0.50%
//   acc 3D bias %    98.4 -> 98.4    worst 97.8908  0.52%, already there
//
// Checked against the spread these have to survive, not against the rule alone.
// The binding records move by 2.8e-5 (yaw) to 3.6e-4 (3D bias) relative between
// a native cascadelake build and an -march=x86-64 one, so the thinnest of these
// margins is 15 times the spread and most are hundreds.  Both builds pass all
// seven.  docs/quality-gate-regauge.md carries that measurement and the command
// that redoes it, which is the check to repeat before cutting any of these
// finer.
//
// Re-derived once more for S_factor = 1.  Taking the horizontal stationary
// acceleration scale from 1.87 to the records' own value moved every gated
// quantity, five of them down:
//
//   Z %Hs JONSWAP    4.72  -> 4.74     worst 4.6952 -> 4.7106
//   Z %Hs PM-Stokes  4.69  -> 4.69     worst 4.6600 -> 4.6580   (unchanged)
//   yaw deg          1.068 -> 1.297    worst 1.0626 -> 1.2896
//   3D % JONSWAP     21.05 -> 20.95    worst 20.9361 -> 20.8367
//   3D % PM-Stokes   20.83 -> 20.86    worst 20.7197 -> 20.7483
//   acc Z bias %     4.93  -> 4.63     worst 4.9054 -> 4.6004
//   acc 3D bias %    98.4  -> 81.84    worst 97.8908 -> 81.4268
//
// The yaw sentinel moving up by 21 percent is not a yaw regression, and the
// distinction matters because a loosened sentinel that hides one would be
// worse than useless.  Yaw on the binding record (jonswap H1.5) spans 1.05 to
// 6.57 deg across five IMU seeds under the *old* constant, so the default-seed
// value the gate is written against is one draw from a wide distribution
// rather than a measure of yaw quality.  Paired across those seeds and all
// eight records, S_factor = 1 lowers yaw RMS by 3.2 percent pooled and on all
// four JONSWAP records -- the deployed default-seed draw happens to be one of
// the few that moves the other way.  reports/results/ou_anisotropy carries
// both.  The three bias and displacement gates that come down are real gains
// and are cut to the rule like the rest.
//
// Then tightened twice more, with the filter standing still for both.
//
// First the quantum.  "Rounded up in the last digit the channel is quoted in"
// only delivers one margin if that digit is a fixed *fraction* of the value,
// and it was a fixed absolute step -- a hundredth, which is 0.2 percent of a
// 4.7 and 0.01 percent of an 81.  The three gates quoted near 4.7 were
// carrying 0.63 to 0.69 percent against a rule that asks for half of one,
// purely from where the decimal fell.  Quoting every channel to four
// significant figures lands all of them between 0.51 and 0.57:
//
//   Z %Hs JONSWAP    4.74 -> 4.735    worst 4.7106   0.63% -> 0.52%
//   Z %Hs PM-Stokes  4.69 -> 4.682    worst 4.6580   0.69% -> 0.52%
//   acc Z bias %     4.63 -> 4.624    worst 4.6004   0.64% -> 0.51%
//
// Then roll and pitch, which this simulator has measured every run and gated
// never.  What the filter has taken on over the last several changes -- the
// Mahony-proxy startup policy, the two-stage magnetic reference, the
// continuous hard-iron correction, and S_factor = 1 -- is mostly attitude
// work, and yaw was the only attitude channel carrying a sentinel.  Both are
// quiet and well clear of their bars across the eight records (roll 0.2374 to
// 0.4179, pitch 0.1716 to 0.2200 deg).  They are gated like yaw, on the
// magnetometer-on protocol only, and for the same reason: dropping the
// magnetometer takes worst-case pitch to 0.2595 deg, past the bar below,
// while worst-case roll goes the other way and improves by 19 percent.  Bars
// fitted with the magnetometer are not measuring the IMU-only filter, and
// scoring it against them would fail it on a channel nobody fitted for it.
//
// All nine limits are what tools/ou_regauge_gates.py prints for the filter
// that ships, and all nine hold on an -march=x86-64 rebuild as well as on a
// native one -- the thinnest margin-to-drift ratio in this family is 169x, on
// pitch.  docs/quality-gate-regauge.md carries that measurement and the
// command that redoes it, which is the check to repeat before cutting any of
// these finer.
// Re-derived once more, for one gate only, when the tuning frequency stopped
// reading the acceleration-band frequency tracker (docs/ou-sigma-horizon.md).
// The change is invisible in displacement -- pooled vertical and 3D RMS both
// move by 0.00% -- but pitch on pmstokes H4.0 goes 0.2200 -> 0.2218 at the
// default seed, which is 0.83% and therefore past the half-percent this gate
// was cut with.
//
// That is a re-draw, not a pitch regression, and the distinction is the same
// one the yaw paragraph above makes.  Paired over five IMU seeds and all eight
// records the pitch ratio is 0.9999 with a 95% interval of [0.9989, 1.0010]:
// no systematic effect at a resolution ten times finer than the move on this
// one record.  tools/ou_sigma_horizon_study.py --axis source --seeds 5
// reproduces it.  Only pitch is re-cut; the other eight limits are untouched
// because nothing this change did moved them.
//
// Then all nine at once, and all nine downward, when the MEKF sensor
// variances were swept for the first time (docs/ou-iii-qmekf-variances.md).
// The gyro term of that sweep is a units correction -- the harness was
// handing a per-sample standard deviation to an argument the filter
// integrates as a noise density -- so the filter that ships now is a
// materially better one rather than the same one re-drawn, and there is no
// re-draw-versus-regression question to settle: nothing got worse.
//
//   Z %Hs JONSWAP    4.735 -> 4.72     worst 4.7106 -> 4.6961
//   Z %Hs PM-Stokes  4.682 -> 4.666    worst 4.6580 -> 4.6426
//   yaw deg          1.297 -> 1.27     worst 1.2896 -> 1.2630
//   roll deg         0.42  -> 0.3637   worst 0.4179 -> 0.3618
//   pitch deg        0.223 -> 0.195    worst 0.2218 -> 0.1940
//   3D % JONSWAP     20.95 -> 13.94    worst 20.8367 -> 13.8686
//   3D % PM-Stokes   20.86 -> 14.92    worst 20.7483 -> 14.8387
//   acc Z bias %     4.624 -> 4.475    worst 4.6004 -> 4.4519
//   acc 3D bias %    81.84 -> 78.61    worst 81.4268 -> 78.2145
//
// The two 3D displacement gates come down by about a third, which is the
// largest single move any of these bars has made.
//
// Then all nine again when the r_S smoothing horizon was re-measured against
// a sea-state transition fast enough to lag -- ADAPT_RS_MULT 3.0 -> 1.5, see
// docs/ou-ema-adaptation-tuning.md.  The bars in this paragraph are that
// change measured on its own, against the tree it was developed on; it and the
// hard-iron re-tune below met in a merge, and the block that ships is
// re-derived from both together at the end.  Eight of the nine move down, so
// this is again mostly a better filter rather than a re-draw:
//
//   Z %Hs JONSWAP    4.72   -> 4.527   worst 4.6961  -> 4.5040
//   Z %Hs PM-Stokes  4.666  -> 4.508   worst 4.6426  -> 4.4850
//   yaw deg          1.27   -> 1.272   worst 1.2630  -> 1.2654
//   roll deg         0.3637 -> 0.3625  worst 0.3618  -> 0.3607
//   pitch deg        0.195  -> 0.1962  worst 0.1940  -> 0.1952
//   3D % JONSWAP     13.94  -> 13.69   worst 13.8686 -> 13.6203
//   3D % PM-Stokes   14.92  -> 14.56   worst 14.8387 -> 14.4810
//   acc Z bias %     4.475  -> 4.488   worst 4.4519  -> 4.4652
//   acc 3D bias %    78.61  -> 78.62   worst 78.2145 -> 78.2238
//
// Pitch is the one gate the shipped bar no longer held (0.1940 -> 0.1952,
// +0.6%), and it is a re-draw of the same kind the paragraph above describes:
// paired over the transition and stationary ensembles the roll/pitch ratio is
// 0.9998 [0.9996, 1.0000] and 0.9997 [0.9994, 1.0000], i.e. attitude improves
// on average while this one record moves inside its own scatter.  Yaw, acc Z
// bias and acc 3D bias move up by 0.1% or less on their binding record and
// are re-cut for the same reason.  Everything here is what
// tools/ou_regauge_gates.py prints for the filter as it now stands, cut to
// the same rule as every line above it.
//
// Then all nine again for the continuous hard-iron re-tune, which cut the
// estimator's absolute ridge floor from 4e-3 to 5e-4 (see
// docs/continuous-mag-hard-iron.md).  Two separate things moved these bars and
// the comment has to keep them apart, because only one of them is this change:
//
//   drift already in the tree.  The bars above were fitted before the reduced
//   physical-MSE integral-regularizer schedule landed, and that schedule took
//   the displacement channels down without a re-gauge.  Measured on this tree
//   with the old ridge, the worsts were already 4.5214 / 4.4872 (vertical),
//   13.6405 / 14.5198 (3D) against bars of 4.72 / 4.666 and 13.94 / 14.92 --
//   3.7 and 1.6 percent of slack that no longer belonged to anyone.
//
//   the re-tune itself.  Yaw, and almost nothing else.  Per record, old ridge
//   -> new: 1.1345 -> 0.6307, 1.2651 -> 0.8780, 0.3031 -> 0.5039,
//   0.7548 -> 0.5307, 0.7049 -> 0.6953, 0.9746 -> 0.5618, 0.4934 -> 0.7339,
//   0.5148 -> 0.4004.  Mean 0.7682 -> 0.6168, worst 1.2651 -> 0.8780.  Every
//   other channel moves in the fourth digit: vertical 4.5214 -> 4.5218, 3D
//   13.6405 -> 13.6518, roll 0.3608 -> 0.3616, pitch 0.1949 -> 0.1957,
//   accelerometer Z bias 4.4728 -> 4.4738.
//
// The yaw bar therefore comes down by 30 percent and the displacement bars by
// the drift, not by the re-tune.  Pitch and the two bias bars go up in the
// third digit, which is the same price the correction has always charged: the
// heading walks onto the corrected field during the run and the horizontal
// accelerometer bias -- the least observable quantity scored here -- absorbs
// part of that motion.  Pitch moving 0.5 percent at the default seed is not a
// pitch regression: paired over five magnetometer-calibration draws and five
// IMU-noise draws, the re-tune moves pooled pitch by +0.1 and +0.4 percent,
// which is a re-draw at ten times this resolution.
//
// Then all nine once more, because the two paragraphs above are two branches
// that met here: the r_S smoothing horizon and the hard-iron ridge floor were
// developed against the same parent and each re-cut all nine bars from its own
// tree, so neither block describes the filter that now ships.  Re-derived on
// the merged tree, against the bars the ridge re-tune left:
//
//   Z %Hs JONSWAP    4.545  -> 4.527   worst 4.5218 -> 4.5044
//   Z %Hs PM-Stokes  4.511  -> 4.509   worst 4.4881 -> 4.4858
//   yaw deg          0.8824 -> 0.8827  worst 0.8780 -> 0.8783
//   roll deg         0.3634 -> 0.3633  worst 0.3616 -> 0.3614
//   pitch deg        0.1967 -> 0.197   worst 0.1957 -> 0.1960
//   3D % JONSWAP     13.73  -> 13.7    worst 13.6518 -> 13.6307
//   3D % PM-Stokes   14.62  -> 14.58   worst 14.5421 -> 14.5033
//   acc Z bias %     4.497  -> 4.489   worst 4.4738 -> 4.4662
//   acc 3D bias %    78.84  -> 78.86   worst 78.4439 -> 78.4668
//
// Every bar the ridge re-tune cut still holds on the merged filter -- none of
// these is a breach -- so this pass only takes back the slack the two changes
// opened in each other's numbers.  The two effects stay separable: the
// displacement and bias channels move by the r_S horizon and yaw by the ridge,
// which is what the per-record numbers in the two paragraphs above already
// showed, and combining them moves nothing beyond the fourth digit either
// paragraph did not predict.
//
// All nine are what tools/ou_regauge_gates.py prints for the filter that
// ships, cut to the same rule as every line above.
static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 4.527f,  // was 4.545,  worst 4.5044 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 4.509f,  // was 4.511,  worst 4.4858 (pmstokes H0.27)
    .err_limit_yaw_deg             = 0.8827f, // was 0.8824, worst 0.8783 (jonswap H1.5)
    .err_limit_roll_deg            = 0.3633f, // was 0.3634, worst 0.3614 (jonswap H4.0)
    .err_limit_pitch_deg           = 0.197f,  // was 0.1967, worst 0.1960 (pmstokes H4.0)
    .err_limit_percent_3d_jonswap  = 13.7f,   // was 13.73,  worst 13.6307 (jonswap H8.5)
    .err_limit_percent_3d_pmstokes = 14.58f,  // was 14.62,  worst 14.5033 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.489f,  // was 4.497,  worst 4.4662 (pmstokes H8.5)
    .bias_3d_percent               = 78.86f,  // was 78.84,  worst 78.4668 (pmstokes H4.0, accel)
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
