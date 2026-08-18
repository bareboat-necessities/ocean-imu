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
    // OU-II's own MEKF sensor variances, as multiples of the ones the shared
    // harness hands every family.  The harness builds each as a fixed multiple
    // of the white noise it injects -- 2.8x on accel, 2.0x on gyro, 1.2x on
    // mag -- and those multiples had never been swept for any family.
    // docs/ou-ii-qmekf-variances.md is the sweep.
    //
    // sigma_g is a units correction, not a fit, and it is the same one OU-III
    // carries: the harness multiplies a *per-sample* gyro standard deviation
    // at 200 Hz, but Kalman3D_Wave_OU_II integrates this argument as a noise
    // *density* (Q_AA = Qbase * Ts), so the deployed value overstated the
    // angular random walk by sqrt(200) = 14.1x on top of its own 2x inflation
    // -- 28.3x in std, 800x in variance.  0.05 puts the argument back on the
    // injected density and keeps a sqrt(2) inflation over it, which is where
    // the measured optimum sits.  One error in two places, not two errors.
    //
    // The other two are empirical: accel to 1.4x the injected white (from
    // 2.8x) and mag to 2.4x (from 1.2x).  The mag value is worth a note,
    // because moving it alone is a *loss* here -- 2x costs pitch and 3D in the
    // one-at-a-time sweep -- and only becomes a gain once sigma_g is
    // corrected.  It was adopted from the joint round, not the axis round.
    //
    // Paired over 8 records x 5 seeds against the previous point: roll -12.5%,
    // pitch -11.6%, yaw -1.2%, vertical displacement -1.7%, 3D displacement
    // -23.8%, accelerometer bias -0.7%, gyro bias -36.6%.  tau_applied and
    // sigma_applied are bit-for-bit unchanged, so the OU schedule is untouched.
    //
    // The SF_SIGMA_*_SCALE overrides below multiply these, so a scale of 1
    // reproduces the deployed point and a re-run of the sweep re-centres on it.
    static constexpr float SIGMA_A_RESCALE = 0.5f;   // 2.8x -> 1.4x injected accel white
    static constexpr float SIGMA_G_RESCALE = 0.05f;  // 2.0x sample std -> sqrt(2)x density
    static constexpr float SIGMA_M_RESCALE = 2.0f;   // 1.2x -> 2.4x injected mag white

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

        // The MEKF variances the Kalman3D_Wave_OU_II constructor takes; see
        // the SIGMA_*_RESCALE block above and docs/ou-ii-qmekf-variances.md.
        // The three sensor sigmas are swept as scale factors on the deployed
        // point, so a scale of 1 leaves it in place.
        if (env_float("SF_SIGMA_A_SCALE", vf)) cfg_.sigma_a *= vf;
        if (env_float("SF_SIGMA_G_SCALE", vf)) cfg_.sigma_g *= vf;
        if (env_float("SF_SIGMA_M_SCALE", vf)) cfg_.sigma_m *= vf;

        if (env_float("SF_PQ0", vf)) cfg_.Pq0 = vf;
        if (env_float("SF_PB0", vf)) cfg_.Pb0 = vf;
        if (env_float("SF_GYRO_BIAS_RW_VAR", vf)) cfg_.b0 = vf;
        if (env_float("SF_RP0_NOISE_VAR", vf)) cfg_.R_p0_noise = vf;
        if (env_float("SF_RV0_NOISE_VAR", vf)) cfg_.R_v0_noise = vf;

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
// scored records plus about half a percent, rounded up in the last digit the
// channel is quoted in -- a tenth for the percentage channels, a hundredth for
// yaw, which at about two degrees would otherwise be handed several times the
// margin the rule asks for.
//
// That margin is deliberately small because the metrics are deterministic, and
// how deterministic is measured rather than assumed: rebuilding at
// -march=x86-64 instead of the host's native cascadelake moves the gated
// numbers here by at most 4.4e-4 relative (accelerometer bias 3D, jonswap
// H8.5), and yaw by 2.4e-4.  The 6e-6 this comment used to claim still holds
// for the simpler observers -- NLO and the PII observer are within 1.5e-6 --
// but not for a filter carrying this many matrix solves.  Half a percent
// leaves better than a factor of ten on every gate below, which is the check
// to redo before cutting one finer.
//
// Setting one below what the filter currently achieves makes it fail every run
// rather than catching a regression.
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
//
// Yaw re-cut to hundredths, 2.2 -> 2.18.  The filter did not move; the tenth
// was worth 1.8 percent of slack against a rule that asks for half of one.
// The six percentage-channel gates are already inside a point of the rule at
// the precision they are quoted in and keep their values.
//
// Re-derived once more for the continuous hard-iron correction, which this
// family now carries on the same defaults as OU-III.  Four of the seven move,
// and they move the same way and for the same reasons they did there.
//
// Yaw halves -- 1.887 to 0.813 deg mean over the eight records, worst 2.161 to
// 1.089 -- because the standing heading error was never a tracking error.  It
// was the hard-iron offset absorbed into the world reference at acquisition,
// which is a gauge, and correcting the stream while moving the reference by a
// delta walks the filter off it.  That gate has to come down with the error or
// it stops being a sentinel: 2.18 -> 1.10, which lands within three hundredths
// of OU-III's own 1.07.
//
// Three go up, and saying so is the point of writing this down, because
// raising a sentinel to admit one's own change is exactly how these stop
// meaning anything:
//
//   acc Z bias %      5.33 -> 5.41   (jonswap H8.5)
//   acc 3D bias %    91.66 -> 93.90  (jonswap H4.0)
//   3D % PM-Stokes   21.02 -> 21.19  (pmstokes H8.5)
//
// The correction walks the heading onto the corrected field during the run,
// and the horizontal accelerometer bias absorbs part of that motion.  It is
// the least observable quantity scored here -- an error above 90 percent of
// the true bias means the error exceeds the thing being estimated, under every
// configuration this family has shipped -- so a two-point move in it is not a
// measurable loss of bias accuracy.  What it is not allowed to be is hidden,
// hence the numbers above.  The displacement channels are flat: vertical mean
// 6.454 -> 6.461 percent of Hs, 3D mean 18.90 -> 19.03 percent, and pitch
// improves (0.289 -> 0.255 deg) for the same reason yaw does.
//
// SF_MAG_CONT_HI=0 is the matched ablation and reproduces the pre-correction
// filter to within 2.6e-4 relative, which is the noise a rebuild of this
// family produces anyway; see docs/quality-gate-regauge.md.  It exceeds the
// yaw gate, as it should -- these are fitted to the filter that ships.  Score
// it with W3D_COLLECT_ALL_GATES.
// Then cut to the rule, which the tenth-quantum was not delivering on the two
// single-digit percentage channels.  A tenth is 1.5 percent of a 6.8 and 1.9
// percent of a 5.4, so rounding a half-percent margin up to one hands back
// three times what the rule asks for; a hundredth costs nothing to write and
// gives it back.  Yaw goes to a thousandth for the same reason -- 1.0949
// rounded to 1.10 is a 0.96 percent bar on a 0.5 percent rule.
//
//   Z %Hs PM-Stokes   6.9  -> 6.85    worst 6.8061, margin 1.38% -> 0.65%
//   yaw deg           1.10 -> 1.095   worst 1.0895, margin 0.96% -> 0.50%
//   acc Z bias %      5.5  -> 5.44    worst 5.4059, margin 1.74% -> 0.63%
//
// The other four were already inside the rule at the precision they are quoted
// in and are left alone.  These three are checked against the drift measured
// for this family rather than against the rule alone: the binding records move
// by 9.4e-5, 1.7e-4 and 3.0e-4 relative between a native and an -march=x86-64
// build, so the smallest of the three margins is still 21 times the spread it
// has to survive.  docs/quality-gate-regauge.md carries that measurement and
// the command that redoes it, which is the check to repeat before cutting any
// of these finer.
//
// Then tightened twice more, with the filter standing still for both -- the
// same two passes OU-III took, since the two families are gated by one rule
// and one script.
//
// First the quantum.  Cutting a tenth to a hundredth, above, fixed the two
// channels where a tenth was worth over a percent, but a fixed absolute step
// still cannot deliver one margin across values from 1 to 94: a hundredth is
// 0.15 percent of a 6.8 and 0.01 percent of a 94.  Quoting every channel to
// four significant figures instead -- so the quantum is a thousandth of the
// value everywhere -- puts all seven between 0.50 and 0.54 percent:
//
//   Z %Hs JONSWAP    6.9  -> 6.899    worst 6.8644   0.52% -> 0.50%
//   Z %Hs PM-Stokes  6.85 -> 6.841    worst 6.8062   0.64% -> 0.51%
//   acc Z bias %     5.44 -> 5.435    worst 5.4073   0.61% -> 0.51%
//   bias 3D %        94.4 -> 94.37    worst 93.8979  0.53% -> 0.50%
//
// Then roll and pitch, which this simulator has measured every run and gated
// never.  What this family has taken on over its last several changes -- the
// OU-III parity work, the Mahony-proxy startup policy, and the continuous
// hard-iron correction, which improved pitch from 0.289 to 0.255 deg mean --
// is largely attitude work, and yaw was the only attitude channel carrying a
// sentinel.  Roll and pitch run 0.2511 to 0.4753 and 0.1801 to 0.3620 deg
// across the eight records, well clear of the bars fitted to them.
//
// They are gated like yaw, on the magnetometer-on protocol only, and for the
// same reason: without it worst-case roll goes to 0.5382 and worst-case pitch
// to 0.7113 deg, both past the bars below.  That is the largest IMU-only
// attitude penalty of the two OU families -- OU-III loses 18 percent of
// worst-case pitch and gains on roll -- and it is a fact about the filter, not
// about the gates, which is why the gates decline to score it.
//
// All nine limits are what tools/ou_regauge_gates.py --family ou_ii prints for
// the filter that ships.  They are checked against this family's own build
// drift rather than against the rule alone, and this family is the noisier of
// the two: the binding records move by 8.0e-6 to 5.6e-4 relative between a
// native and an -march=x86-64 build, and both builds pass all nine.  The
// thinnest margin-to-drift ratio in the set is pitch at 9.3x -- the tightest
// anywhere in the five families, and the first bar to re-measure rather than
// re-cut if a rebuild ever breaches it.  docs/quality-gate-regauge.md carries
// that measurement and the command that redoes it.
//
// Re-derived for the (r_p0, r_v0) coefficient re-fit,
// docs/ou-ii-pseudo-variance-tuning.md.  Every one of the nine still passed on
// its previous value, so this pass is the rule following a filter that moved,
// not a breach being papered over -- but pitch had come down to 0.0001 deg of
// margin against a channel whose measured -march rebuild drift is about 2e-4
// deg, so leaving the set alone would have shipped a bar that a rebuild
// decides, and not the filter.  Reverting both coefficients through
// OU_R_P0_COEFF/OU_R_V0_COEFF puts all nine back at the rule to the digit,
// which is the check that this set moved for the re-fit and for nothing else.
//
// Five tighten and four loosen:
//
//   Z %Hs JONSWAP    6.899 -> 6.865   worst 6.8644 -> 6.8300 (jonswap H0.27)
//   Z %Hs PM-Stokes  6.841 -> 6.848   worst 6.8062 -> 6.8139 (pmstokes H0.27)
//   yaw deg          1.095 -> 1.089   worst 1.0895 -> 1.0833 (jonswap H1.5)
//   roll deg        0.4778 -> 0.4792  worst 0.4753 -> 0.4768 (jonswap H4.0)
//   pitch deg       0.3639 -> 0.3657  worst 0.3620 -> 0.3638 (jonswap H8.5)
//   3D % JONSWAP      21.1 -> 20.92   worst 20.9867 -> 20.8140 (jonswap H1.5)
//   3D % PM-Stokes    21.3 -> 21.03   worst 21.1935 -> 20.9203 (pmstokes H8.5)
//   acc Z bias %     5.435 -> 5.324   worst 5.4073 -> 5.2969 (jonswap H8.5)
//   bias 3D %        94.37 -> 94.47   worst 93.8979 -> 93.9911 (jonswap H4.0)
//
// Both displacement gates come down, by 0.9 and 1.3 percent, which is where a
// change to the translational regularizer is supposed to show up.
//
// Of the four that loosen, three are single-realization moves against an
// ensemble that goes the other way.  Pooled over the eight records at six IMU
// seed triplets the re-fit reads 0.9958 vertical, 0.9796 pitch, 1.0003 roll and
// 1.0006 accelerometer-bias 3D, so pitch improves 2 percent while its
// deterministic worst record moves up 0.0018 deg, and roll is flat within its
// own realization noise.  The vertical one is not noise: it is the small-sea
// end of the position-coefficient trade, and holding the per-record mean
// vertical error where it is -- seven of eight records improve and the eighth
// is 1.0003 -- is what bounded R_p0_coeff at 0.65.  The bias one is the least
// observable quantity scored here, with an error above 90 percent of the true
// bias on the binding record under every configuration this family has shipped.
// Then all nine at once, and all nine downward, when the MEKF sensor
// variances were swept for the first time (docs/ou-ii-qmekf-variances.md).
// The gyro term of that sweep is a units correction -- the harness was
// handing a per-sample standard deviation to an argument the filter
// integrates as a noise density -- so the filter that ships now is a
// materially better one rather than the same one re-drawn, and no channel
// moved the wrong way.  OU-III took the identical correction in the same
// round; the two families share the harness that supplied the argument.
//
//   Z %Hs JONSWAP    6.865  -> 6.776    worst 6.8300 -> 6.7420
//   Z %Hs PM-Stokes  6.848  -> 6.803    worst 6.8139 -> 6.7688
//   yaw deg          1.089  -> 1.074    worst 1.0833 -> 1.0681
//   roll deg         0.4792 -> 0.4352   worst 0.4768 -> 0.4330
//   pitch deg        0.3657 -> 0.279    worst 0.3638 -> 0.2775
//   3D % JONSWAP     20.92  -> 16.54    worst 20.8140 -> 16.4569
//   3D % PM-Stokes   21.03  -> 17.67    worst 20.9203 -> 17.5728
//   acc Z bias %     5.324  -> 4.802    worst 5.2969 -> 4.7776
//   acc 3D bias %    94.47  -> 92.35    worst 93.9911 -> 91.8873
//
// Everything here is what tools/ou_regauge_gates.py prints for the filter as
// it now stands, cut to the same rule as every line above it.
static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 6.776f,  // was 6.865, worst 6.7420 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 6.803f,  // was 6.848, worst 6.7688 (pmstokes H0.27)
    .err_limit_yaw_deg             = 1.074f,  // was 1.089, worst 1.0681 (jonswap H0.27)
    .err_limit_roll_deg            = 0.4352f, // was 0.4792, worst 0.4330 (jonswap H4.0)
    .err_limit_pitch_deg           = 0.279f,  // was 0.3657, worst 0.2775 (jonswap H8.5)
    .err_limit_percent_3d_jonswap  = 16.54f,  // was 20.92, worst 16.4569 (jonswap H1.5)
    .err_limit_percent_3d_pmstokes = 17.67f,  // was 21.03, worst 17.5728 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.802f,  // was 5.324, worst 4.7776 (jonswap H8.5)
    .bias_3d_percent               = 92.35f,  // was 94.47, worst 91.8873 (jonswap H4.0, accel)
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
