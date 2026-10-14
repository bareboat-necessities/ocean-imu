#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Sea-state orchestration for the right-invariant two-frame Lie-group filter.

    The translational plant remains the integrated OU chain [v,p,S,a_w], but
    the estimator geometry is TFG-specific.  The measurement-only front end is
    intentionally kept at coefficient parity with deployed OU-III: canonical
    log-period statistics, period-scaled acceleration variance, self-scaled
    parameter EMAs, startup tilt/magnetic acquisition and continuous hard iron.
    Where that parity is a shared constant or a shared mechanism it is the
    same code as OU-II/OU-III (src/kalman_common/): the category-A defaults,
    the accelerometer vibration conditioning, the world-frame gravity gate and
    the continuous hard-iron slew.  What stays here is what TFG does
    differently: its own staged/proxy startup state machine and elapsed-time
    scheduling, the per-sample north-frame magnetic acquisition that tracks
    proxy yaw drift, the Lie-group handoff and covariance seed, and its
    adaptation, which smooths channel targets that are only formed after the
    tuner warmup and commits them on an accumulated-time cadence.

    r_S defaults to the same reduced spectral-MSE law as current OU-III,

        r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S),

    with T_S proportional to tau away from cadence clamps, giving the effective
    tau^(41/14) dependence.  The selectable LegacyCubic law is

        r_S,base = C_R sigma_aw tau^3,
        r_S,filter = r_S,base sqrt(T_0/T_S)

    It provides a lower-cost alternative for embedded targets.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "kalman_common/AccelVibrationConditioning.h"
#include "kalman_common/MagneticStartupCommon.h"
#include "kalman_common/SeaStateFusionDefaults.h"
#include "kalman_common/SeaStateFusionFilterCommon.h"
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#include "tuner/AdaptiveWaveBandPass.h"
#include "tuner/ContinuousMagHardIronEstimator.h"
#include "tuner/MagAutoTuner.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/WavePeriodEstimator.h"

namespace ocean_imu::tfg {

// Shared front-end, adaptation-mechanics and startup defaults; values and
// provenance in kalman_common/SeaStateFusionDefaults.h.
namespace defaults = ::seastate::common::defaults;

constexpr float MAG_DELAY_SEC = defaults::MAG_DELAY_SEC;
constexpr float SIGMA_BAND_LOW_RATIO_DEFAULT  = defaults::SIGMA_BAND_LOW_RATIO;
constexpr float SIGMA_BAND_HIGH_RATIO_DEFAULT = defaults::SIGMA_BAND_HIGH_RATIO;
constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = defaults::SIGMA_BAND_MIN_HZ;
constexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = defaults::SIGMA_BAND_MAX_HZ;
constexpr float ADAPT_TAU_SEC_DEFAULT          = defaults::ADAPT_TAU_SEC;
constexpr float ADAPT_TAU_SEA_PERIODS_DEFAULT  = defaults::ADAPT_TAU_SEA_PERIODS;
// Same r_S smoothing multiplier as deployed OU-III.  A property of the r_S
// law rather than of the front end, so it is TFG's own constant.
constexpr float ADAPT_RS_MULT_DEFAULT          = 1.5f;
constexpr float TUNER_SIGMA_VAR_K_PERIODS_DEFAULT = 4.0f;

// Same strong-observation sensor point and analytical spectral coefficient as
// deployed OU-III.  The physical wave RMS is recovered from sigma_aw below, so
// TFG's sigma_coeff does not alter the distortion cost.
// Both are coefficients of the r_S law, so they stay TFG's own constants.
constexpr float TFG_NOMINAL_DT = defaults::NOMINAL_IMU_DT_S;
constexpr float R_S_ACCEL_NOISE_DENSITY_DEFAULT =
    0.0148f * 0.0148f * TFG_NOMINAL_DT;
constexpr float R_S_MSE_COEFF_DEFAULT = 0.0538f;

// Front-end accelerometer vibration guard and vibration-aware accelerometer
// covariance, armed by default at the shared operating point.  The defect is
// in the attitude loop, not in the translational state: TFG levels from the
// same private Mahony observer at the same gains and feeds the same
// accelerometer to its own MEKF, so it inherits both the defect and the
// remedy.  See defaults::ACC_VIBRATION_GUARD_HZ for the measurement.
constexpr float ACC_VIBRATION_GUARD_HZ_DEFAULT    = defaults::ACC_VIBRATION_GUARD_HZ;
constexpr int   ACC_VIBRATION_GUARD_POLES_DEFAULT = defaults::ACC_VIBRATION_GUARD_POLES;
constexpr float ACC_VIBRATION_RACC_GAIN_DEFAULT   = defaults::ACC_VIBRATION_RACC_GAIN;

enum class RSAdaptationLaw : uint8_t {
    LegacyCubic = 0,
    SpectralMSE = 1,
};

struct TfgTuneState {
    float tau_applied   = 1.1f;
    float sigma_applied = 1e-2f;
    float RS_applied    = 0.5f;
};

template <typename MekfT = ocean_imu::kalman::Kalman3D_Wave_TFG<float>>
class SeaStateFusionFilter_TFG : public ::seastate::common::AccelVibrationConditioning {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Mekf = MekfT;
    using Vector3f = Eigen::Vector3f;
    enum class StartupStage { Cold, TunerWarm, Live };
    enum class StartupInitPolicy { StagedMekf, MahonyProxy };
    using RSLaw = RSAdaptationLaw;

    struct Config {
        Vector3f sigma_a{Vector3f::Constant(0.5f)};
        float    gyro_noise_density = 0.005f;
        Vector3f sigma_m{Vector3f::Constant(0.1f)};
        float    gyro_bias_rw_var   = 3.75e-10f;
        float    initial_covariance = 1e-4f;
        bool     with_mag = true;
        float    mag_delay_sec = MAG_DELAY_SEC;
        bool     freeze_acc_bias_until_live = true;
        float    Racc_warmup_std = 0.5f;
        float    gravity_magnitude = 9.80665f;

        // Deployed OU wrapper value.  The underlying low-level OU classes have
        // a 5 s compatibility default, but the deployed front ends set 10 s.
        float    online_tune_warmup_sec = defaults::STARTUP_ONLINE_TUNE_WARMUP_SEC;
        StartupInitPolicy startup_init_policy = StartupInitPolicy::MahonyProxy;

        float proxy_startup_min_sec     = defaults::PROXY_STARTUP_MIN_SEC;
        float proxy_startup_timeout_sec = defaults::PROXY_STARTUP_TIMEOUT_SEC;
        float proxy_handoff_tilt_sigma_rad = defaults::PROXY_HANDOFF_TILT_SIGMA_RAD;
        float proxy_handoff_yaw_sigma_rad  = defaults::PROXY_HANDOFF_YAW_SIGMA_RAD;
        float proxy_handoff_yaw_sigma_free_rad = defaults::PROXY_HANDOFF_YAW_SIGMA_FREE_RAD;

        // Same private Mahony gains as the OU front ends.
        float proxy_two_kp = defaults::STARTUP_PROXY_TWO_KP;
        float proxy_two_ki = defaults::STARTUP_PROXY_TWO_KI;

        // Compatibility field retained for existing TFG studies.  At the
        // reference 25 Hz mag ODR, 10 s corresponds to OU-III's 250 updates.
        float acc_bias_unlock_sec = 10.0f;
        float handoff_acc_bias_std = 0.03f;

        // The actual deployed OU-III magnetic gravity gate is in WORLD frame.
        // These fields use the shared world-frame gate thresholds and timing.
        float proxy_gravity_align_sin = defaults::GRAVITY_GATE_MAX_SIN;
        float proxy_gravity_lpf_sec   = defaults::GRAVITY_GATE_LPF_SEC;
        float proxy_gravity_hold_sec  = defaults::GRAVITY_GATE_HOLD_SEC;
        float proxy_gravity_warmup_sec = defaults::GRAVITY_GATE_WARMUP_SEC;
        float mag_extreme_gyro_dps    = defaults::MAG_EXTREME_GYRO_DPS;
        float mag_tilt_fallback_sec   = defaults::MAG_TILT_FALLBACK_SEC;
        float mag_init_min_mag_norm   = defaults::MAG_INIT_MIN_MAG_NORM;

        int   mag_min_samples    = defaults::MAG_MIN_SAMPLES;
        float mag_min_window_sec = defaults::MAG_MIN_WINDOW_SEC;
        float mag_max_window_sec = defaults::MAG_MAX_WINDOW_SEC;
        float mag_sample_dt_sec  = defaults::MAG_SAMPLE_DT_SEC;
        float proxy_mag_settle_sec = defaults::PROXY_MAG_SETTLE_SEC;

        // TFG refines at 30 s where the OU front ends wait 90 s: its proxy
        // yaw offset is re-gauged on every magnetometer sample before
        // handoff, so the provisional reference is less of a liability.
        bool  mag_refine_enabled    = true;
        float mag_refine_start_sec  = 30.0f;
        float mag_refine_window_sec = defaults::MAG_REFINE_WINDOW_SEC;

        bool  mag_enable_quality_weighting = false;
        float mag_min_effective_weight     = 0.0f;
        float mag_acc_norm_rel_soft        = defaults::MAG_ACC_NORM_REL_SOFT;
        float mag_gyro_soft_dps            = defaults::MAG_GYRO_SOFT_DPS;
        bool  mag_estimate_hard_iron       = false;

        bool  mag_continuous_hard_iron        = true;
        float mag_hi_memory_sec               = defaults::MAG_HI_MEMORY_SEC;
        float mag_hi_model_ridge              = defaults::MAG_HI_MODEL_RIDGE;
        float mag_hi_model_ridge_relative     = defaults::MAG_HI_MODEL_RIDGE_RELATIVE;
        float mag_hi_min_information          = defaults::MAG_HI_MIN_INFORMATION;
        float mag_hi_min_effective_weight     = defaults::MAG_HI_MIN_EFFECTIVE_WEIGHT;
        float mag_hi_max_residual_rms_uT      = defaults::MAG_HI_MAX_RESIDUAL_RMS_UT;
        float mag_hi_max_bias_fraction        = defaults::MAG_HI_MAX_BIAS_FRACTION;
        float mag_hi_apply_fraction           = defaults::MAG_HI_APPLY_FRACTION;
        float mag_hi_slew_tau_sec             = defaults::MAG_HI_SLEW_TAU_SEC;

        // Front-end accelerometer vibration guard and the vibration-aware
        // accelerometer covariance it drives; see the defaults above.  A zero
        // cutoff removes the guard and restores the unconditioned measurement
        // path exactly, and a zero gain leaves the commanded covariance to the
        // stage logic alone.
        float acc_vibration_guard_hz    = ACC_VIBRATION_GUARD_HZ_DEFAULT;
        int   acc_vibration_guard_poles = ACC_VIBRATION_GUARD_POLES_DEFAULT;
        float acc_vibration_racc_gain   = ACC_VIBRATION_RACC_GAIN_DEFAULT;

        bool  wave_band_tuning      = true;
        float sigma_band_low_ratio  = SIGMA_BAND_LOW_RATIO_DEFAULT;
        float sigma_band_high_ratio = SIGMA_BAND_HIGH_RATIO_DEFAULT;
        float sigma_band_min_hz     = SIGMA_BAND_MIN_HZ_DEFAULT;
        float sigma_band_max_hz     = SIGMA_BAND_MAX_HZ_DEFAULT;
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;
        mekf_ = Mekf(cfg.gyro_noise_density, cfg.gravity_magnitude);
        mekf_.initialize_identity(cfg.initial_covariance);
        mekf_.set_Racc_std(cfg.sigma_a);
        mekf_.set_Rmag_std(cfg.sigma_m);
        mekf_.set_Q_bgyro_rw(Vector3f::Constant(cfg.gyro_bias_rw_var));
        Racc_nominal_ = cfg.sigma_a;

        tuner_ = ::SeaStateAutoTuner(TUNER_SIGMA_VAR_K_PERIODS_DEFAULT);
        sea_time_sec_ = 0.5f / kTuneFreqPriorHz;
        wave_period_.reset();
        vertical_complementary_.setGains(cfg.proxy_two_kp, cfg.proxy_two_ki);
        vertical_complementary_.reset();
        setAccelVibrationGuard(cfg.acc_vibration_guard_hz, cfg.acc_vibration_guard_poles);
        setAccelVibrationRaccGain(cfg.acc_vibration_racc_gain);
        accel_guard_.reset();
        racc_inflated_ = false;
        racc_effective_.setZero();
        sigma_wave_band_.setRatios(cfg.sigma_band_low_ratio, cfg.sigma_band_high_ratio);
        sigma_wave_band_.setLimitsHz(cfg.sigma_band_min_hz, cfg.sigma_band_max_hz);
        sigma_wave_band_.reset();
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_slow_lpf_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
        elapsed_sec_ = 0.0f;
        live_sec_ = 0.0f;
        mag_elapsed_sec_ = 0.0f;
        tuner_warm_sec_ = 0.0f;
        pseudo_elapsed_ = 0.0f;
        adapt_elapsed_sec_ = 0.0f;
        stage_ = StartupStage::Cold;
        beginMagAcquisition_();
        enterCold_();
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void update(float dt, const Vector3f& gyro, const Vector3f& acc, float tempC = 35.0f) {
        if (!(dt > 0.0f) || !std::isfinite(dt) || !gyro.allFinite() || !acc.allFinite()) return;

        applyPendingTune_();
        elapsed_sec_ += dt;

        // Strip out-of-band accelerometer vibration before the attitude path
        // reads it.  This is the one place raw measurements enter, so
        // filtering here is what keeps the guard's effect describable: the
        // proxy, the MEKF and the staged bootstrap all see the same
        // conditioned signal, and no part of the attitude loop can be left on
        // a different version of the accelerometer.
        //
        // The magnetic gravity gate and the tilt frame the magnetometer is
        // acquired in stay on the raw signal, as they do in deployed OU-III,
        // where the guard sits inside the fusion filter and those live in the
        // wrapper above it.  Both average the specific force over ten seconds
        // or more, so the machinery band is already far below their corner;
        // keeping them raw costs nothing and keeps the front end at signal
        // parity with OU-III rather than only coefficient parity.
        //
        // Armed by default, and transparent below its detector's lower rail,
        // in which case acc_in is acc unchanged.
        const Vector3f acc_in = conditionAccel_(acc, dt);

        vertical_complementary_.update(dt, gyro, acc_in, cfg_.gravity_magnitude);
        last_acc_body_ = acc;
        last_gyro_body_ = gyro;
        have_last_imu_ = true;
        updateProxyGravityQuality_(dt, gyro, acc);

        const float a_up = vertical_complementary_.verticalAccelUpMs2();
        // OU-III feeds the canonical period estimator every sample.  Readiness
        // is a trust gate on the statistic, not permission to update it.
        if (std::isfinite(a_up)) wave_period_.update(dt, a_up);
        updateTuner_(dt, a_up);

        if (stage_ != StartupStage::Live) {
            tuner_warm_sec_ += dt;
        } else {
            live_since_sec_ += dt;
            maybeUnlockAccBias_();
        }

        if (stage_ == StartupStage::Cold) {
            if (cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy) {
                tryProxyHandoff_();
            } else {
                stagedColdStep_(gyro, acc_in, dt);
            }
            return;
        }

        // Tell the MEKF how much it should trust that sample before it uses
        // it, so the covariance and the measurement describe the same
        // conditions.  A no-op unless a gain is set and the guard sees
        // machinery.
        applyRaccVibrationInflation_();

        periodicAwCovSyncTick_(dt);
        mekf_.time_update(gyro, dt);
        mekf_.measurement_update_acc_only(acc_in, tempC);

        // Unlike the OU wrappers, whose MEKFs fire their own zero
        // pseudo-measurements on set_pseudo_update_period_s(), TFG fires the
        // S = 0 pseudo-measurement from here on an accumulated-time timer.
        // The period is the same tau-scaled cadence (pseudoUpdatePeriodFor_).
        pseudo_elapsed_ += dt;
        if (pseudo_elapsed_ >= pseudo_period_sec_) {
            pseudo_elapsed_ = 0.0f;
            mekf_.applyIntegralZeroPseudoMeas();
        }

        if (stage_ == StartupStage::TunerWarm) {
            live_sec_ += dt;
            if (live_sec_ >= cfg_.online_tune_warmup_sec &&
                wave_period_.hasUsablePeriod()) {
                enterLive_();
            }
        } else {
            adaptMekf_(dt);
        }
    }

    void updateMag(const Vector3f& mag_body) {
        if (!cfg_.with_mag || !mag_body.allFinite()) return;
        if (!(mag_body.norm() > 1e-9f)) return;
        if (elapsed_sec_ < cfg_.mag_delay_sec) return;

        accumulateContinuousHardIron_(mag_body);
        if (usingProxyInit_()) {
            proxyUpdateMag_(mag_body);
            return;
        }
        if (stage_ == StartupStage::Cold) return;

        if (!mekf_.has_magnetic_reference()) {
            const Vector3f b_world = mekf_.R_bw() * mag_body;
            const float horiz = std::hypot(b_world.x(), b_world.y());
            if (horiz > 1e-9f) {
                const float psi = std::atan2(b_world.y(), b_world.x());
                mekf_.apply_world_yaw_gauge(-psi);
            }
            mekf_.set_magnetic_reference_world(Vector3f(horiz, 0.0f, b_world.z()));
            return;
        }
        mekf_.measurement_update_mag_only(mag_body);
    }

    [[nodiscard]] bool magReferenceLearned() const noexcept { return mag_reference_learned_; }
    [[nodiscard]] bool magReferenceRefined() const noexcept { return mag_refine_done_; }
    [[nodiscard]] float magRefineTimeSec() const noexcept { return mag_refine_time_sec_; }
    [[nodiscard]] float magNorthLockTimeSec() const noexcept { return mag_north_lock_time_sec_; }
    [[nodiscard]] const Vector3f& magHardIronBodyUT() const noexcept { return mag_hard_iron_body_uT_; }
    [[nodiscard]] const Vector3f& magContinuousHardIronAppliedUT() const noexcept { return hard_iron_.applied_body_uT; }
    [[nodiscard]] const ContinuousMagHardIronEstimator& magContinuousHardIron() const noexcept { return hard_iron_.estimator; }
    [[nodiscard]] StartupInitPolicy startupInitPolicy() const noexcept { return cfg_.startup_init_policy; }
    [[nodiscard]] bool isTunerReady() const noexcept {
        return wave_period_.hasUsablePeriod() &&
               tuner_warm_sec_ >= cfg_.online_tune_warmup_sec;
    }
    [[nodiscard]] float pseudoUpdatePeriodSec() const noexcept { return pseudo_period_sec_; }
    [[nodiscard]] bool handoffTimedOut() const noexcept { return handoff_timed_out_; }
    [[nodiscard]] float getRSFilterInput() const noexcept { return RS_filter_input_; }

    // Accelerometer vibration guard and vibration-aware covariance: the
    // shared seastate::common::AccelVibrationConditioning API.

    void setAdaptEverySecs(float s) { if (s >= 0.0f && std::isfinite(s)) adapt_every_secs_ = s; }
    void setTauScaledPseudoCadence(bool on) { tau_scaled_pseudo_cadence_ = on; applyPseudoCadence_(); }
    [[nodiscard]] bool tauScaledPseudoCadence() const noexcept { return tau_scaled_pseudo_cadence_; }

    bool setFixedTuning(float tau_s, float sigma_a, float RS) {
        if (!(tau_s > 0.0f) || !(sigma_a >= 0.0f) || !(RS > 0.0f)) return false;
        if (!std::isfinite(tau_s) || !std::isfinite(sigma_a) || !std::isfinite(RS)) return false;
        fixed_tuning_ = true;
        tune_.tau_applied = tau_s;
        tune_.sigma_applied = sigma_a;
        tune_.RS_applied = RS;
        tau_target_ = tau_s;
        sigma_target_ = sigma_a;
        RS_target_ = RS;
        commitTune_();
        return true;
    }

    bool setChannelFreeze(bool freeze_ou, float tau_s, float sigma_a,
                          bool freeze_RS, float RS) {
        if (freeze_ou && freeze_RS) return false;
        if (freeze_ou) {
            if (!(tau_s > 0.0f) || !(sigma_a >= 0.0f)) return false;
            freeze_ou_channel_ = true;
            tune_.tau_applied = tau_s;
            tune_.sigma_applied = sigma_a;
            tau_target_ = tau_s;
            sigma_target_ = sigma_a;
        }
        if (freeze_RS) {
            if (!(RS > 0.0f)) return false;
            freeze_RS_channel_ = true;
            tune_.RS_applied = RS;
            RS_target_ = RS;
        }
        commitTune_();
        return true;
    }

    void enableTuner(bool on) { enable_tuner_ = on; }
    void enableLinearBlock(bool on) { mekf_.set_linear_block_enabled(on); }
    void setTauCoeff(float c)      { if (c > 0.0f && std::isfinite(c)) tau_coeff_ = c; }
    void setSigmaCoeff(float c)    { if (c > 0.0f && std::isfinite(c)) sigma_coeff_ = c; }
    void setRSCoeff(float c)       { if (c > 0.0f && std::isfinite(c)) R_S_coeff_ = c; }
    void setSFactor(float s)       { if (s > 0.0f && std::isfinite(s)) S_factor_ = s; }
    void setRSXFactor(float k)     { if (k > 0.0f && std::isfinite(k)) R_S_x_factor_ = k; }
    void setRSYFactor(float k)     { if (k > 0.0f && std::isfinite(k)) R_S_y_factor_ = k; }
    [[nodiscard]] float getRSXFactor() const noexcept { return R_S_x_factor_; }
    [[nodiscard]] float getRSYFactor() const noexcept { return R_S_y_factor_; }
    void setAccNoiseFloorSigma(float s) { if (s >= 0.0f && std::isfinite(s)) noise_floor_sigma_ = s; }

    void setRSLaw(RSLaw law) noexcept { rs_law_ = law; }
    [[nodiscard]] RSLaw getRSLaw() const noexcept { return rs_law_; }
    void setEmbeddedFriendlyLegacyRSLaw(bool legacy) noexcept {
        rs_law_ = legacy ? RSLaw::LegacyCubic : RSLaw::SpectralMSE;
    }
    void setRSMseCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) rs_mse_coeff_ = c;
    }
    void setRSAccelNoiseDensity(float r_a) {
        if (std::isfinite(r_a) && r_a > 0.0f) {
            rs_accel_noise_density_ = r_a;
            refreshQeffPow_();
        }
    }
    [[nodiscard]] float getRSMseCoeff() const noexcept { return rs_mse_coeff_; }
    [[nodiscard]] float getRSAccelNoiseDensity() const noexcept { return rs_accel_noise_density_; }

    void setAdaptationTimeConstants(float tau_sec) {
        if (tau_sec > 0.0f && std::isfinite(tau_sec)) {
            adapt_tau_sec_ = tau_sec;
            adapt_tau_sea_periods_ = 0.0f;
        }
    }
    void setAdaptationSeaPeriods(float periods) {
        if (periods > 0.0f && std::isfinite(periods)) adapt_tau_sea_periods_ = periods;
    }
    void setTunerFreqSmoothingSeaPeriods(float /*periods*/) {}
    void setTunerFreqSmoothingTimeConstant(float /*tau_sec*/) {}
    [[nodiscard]] float getAdaptationSeaPeriods() const noexcept { return adapt_tau_sea_periods_; }
    [[nodiscard]] float getTunerFreqSmoothingSeaPeriods() const noexcept {
        return tuner_.getFrequencySmoothingSeaPeriods();
    }
    [[nodiscard]] float getSigmaVarianceHorizonPeriods() const noexcept { return tuner_.getKPeriods(); }
    [[nodiscard]] float getSigmaVarianceHorizonSec() const noexcept { return tuner_.getVarianceHorizonSec(); }

    void setRSAdaptMult(float m)     { if (m > 0.0f && std::isfinite(m)) adapt_RS_mult_ = m; }
    void setRSAdaptSlewLog(float d)  { if (d >= 0.0f && std::isfinite(d)) adapt_RS_slew_log_ = d; }
    [[nodiscard]] float getRSAdaptMult() const noexcept { return adapt_RS_mult_; }
    void setTauBounds(float lo, float hi) { if (lo > 0.0f && hi > lo) { min_tau_ = lo; max_tau_ = hi; } }
    void setRSBounds(float lo, float hi) { if (lo > 0.0f && hi > lo) { min_RS_ = lo; max_RS_ = hi; } }
    void setMaxSigmaA(float m)       { if (m > 0.0f && std::isfinite(m)) max_sigma_a_ = m; }
    void setWaveBandTuning(bool on)  { cfg_.wave_band_tuning = on; }
    [[nodiscard]] bool waveBandTuning() const noexcept { return cfg_.wave_band_tuning; }
    void setPeriodicAwCovSync(bool on) { periodic_aw_cov_sync_ = on; }
    [[nodiscard]] bool periodicAwCovSync() const noexcept { return periodic_aw_cov_sync_; }
    void setSigmaBandRatios(float low, float high) {
        cfg_.sigma_band_low_ratio = low;
        cfg_.sigma_band_high_ratio = high;
        sigma_wave_band_.setRatios(low, high);
    }

    [[nodiscard]] Mekf& mekf() noexcept { return mekf_; }
    [[nodiscard]] const Mekf& mekf() const noexcept { return mekf_; }
    [[nodiscard]] StartupStage stage() const noexcept { return stage_; }
    [[nodiscard]] bool isLive() const noexcept { return stage_ == StartupStage::Live; }
    [[nodiscard]] float getTauApplied()   const noexcept { return tune_.tau_applied; }
    [[nodiscard]] float getSigmaApplied() const noexcept { return tune_.sigma_applied; }
    [[nodiscard]] float getRSApplied()    const noexcept { return tune_.RS_applied; }
    [[nodiscard]] float getTauTarget()    const noexcept { return tau_target_; }
    [[nodiscard]] float getSigmaTarget()  const noexcept { return sigma_target_; }
    [[nodiscard]] float getRSTarget()     const noexcept { return RS_target_; }
    [[nodiscard]] float getWavePeriodSec() const noexcept { return wave_period_.getPeriodSec(); }
    [[nodiscard]] bool  wavePeriodUsable() const noexcept { return wave_period_.hasUsablePeriod(); }
    [[nodiscard]] bool  wavePeriodReady() const noexcept { return wave_period_.isReady(); }
    [[nodiscard]] float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }
    [[nodiscard]] Eigen::Quaternionf quaternion() const { return mekf_.quaternion(); }

    // Read-only startup tilt for application compass output. The core state
    // remains untouched until the existing handoff gates admit Live. This
    // quaternion has no magnetic yaw; callers must not publish its yaw as HDG.
    [[nodiscard]] bool startupTiltQuaternion(Eigen::Quaternionf& q_bw) const {
        if (stage_ != StartupStage::Cold || !usingProxyInit_() ||
            !vertical_complementary_.isInitialized()) return false;
        q_bw = vertical_complementary_.tiltQuaternion();
        return q_bw.coeffs().allFinite() && q_bw.squaredNorm() > 1.0e-12f;
    }

    [[nodiscard]] Vector3f get_velocity() const { return mekf_.get_velocity(); }
    [[nodiscard]] Vector3f get_position() const { return mekf_.get_position(); }
    [[nodiscard]] Vector3f get_world_accel() const { return mekf_.get_world_accel(); }

private:
    void beginMagAcquisition_() {
        mag_auto_tuner_.setConfig(::seastate::common::magAutoTunerConfig(cfg_));
        hard_iron_.configure(cfg_);

        mag_reference_learned_ = false;
        mag_world_ref_valid_ = false;
        mag_world_ref_uT_.setZero();
        mag_proxy_yaw_offset_rad_ = 0.0f;
        have_mag_proxy_gauge_ = false;
        mag_hard_iron_body_uT_.setZero();
        hard_iron_.resetApplied();
        mag_refine_started_ = false;
        mag_refine_done_ = false;
        mag_refine_time_sec_ = NAN;
        mag_north_lock_time_sec_ = NAN;
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;
        gravity_gate_.reset();
        last_acc_body_.setZero();
        last_gyro_body_.setZero();
        have_last_imu_ = false;
        acc_bias_unlocked_ = false;
        live_since_sec_ = 0.0f;
        acc_bias_hold_ = usingProxyInit_() && cfg_.with_mag && cfg_.mag_refine_enabled;
    }

    void enterCold_() {
        mekf_.set_linear_block_enabled(false);
        if (cfg_.freeze_acc_bias_until_live) mekf_.set_acc_bias_updates_enabled(false);
        mekf_.set_Racc_std(Vector3f::Constant(cfg_.Racc_warmup_std));
        warmup_Racc_active_ = true;
        racc_inflated_ = false;
    }

    // The accelerometer sigma the startup and stage logic wants, before any
    // vibration inflation.  Returns a zero vector when it is not known, which
    // is the signal to leave the commanded covariance alone.
    [[nodiscard]] Vector3f raccBaseStd_() const {
        if (warmup_Racc_active_ && cfg_.Racc_warmup_std > 0.0f) {
            return Vector3f::Constant(cfg_.Racc_warmup_std);
        }
        if (Racc_nominal_.allFinite() && Racc_nominal_.minCoeff() > 0.0f) {
            return Racc_nominal_;
        }
        return Vector3f::Zero();
    }

    // Vibration inflation of the stage base.  No low-wave noise weighting
    // here (unlike the OU wrappers), and nothing runs with a zero gain.
    void applyRaccVibrationInflation_() {
        if (!(racc_vibration_gain_ > 0.0f)) return;

        const Vector3f base = raccBaseStd_();
        if (!(base.minCoeff() > 0.0f)) return;

        commandRaccStd_(mekf_, base, Vector3f::Ones());
    }

    void enterLive_() {
        stage_ = StartupStage::Live;
        live_since_sec_ = 0.0f;
        mekf_.set_linear_block_enabled(true);
        acc_bias_unlocked_ = false;
        mekf_.set_acc_bias_updates_enabled(false);
        maybeUnlockAccBias_();
        mekf_.set_Racc_std(Racc_nominal_);
        warmup_Racc_active_ = false;
        racc_inflated_ = false;
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void maybeUnlockAccBias_() {
        if (acc_bias_unlocked_) return;
        if (cfg_.freeze_acc_bias_until_live) {
            if (live_since_sec_ < cfg_.acc_bias_unlock_sec) return;
            if (acc_bias_hold_ && !mag_refine_done_) return;
        }
        acc_bias_unlocked_ = true;
        mekf_.set_acc_bias_updates_enabled(true);
    }

    [[nodiscard]] float bandNoiseFloorSigma_() const noexcept {
        if (!cfg_.wave_band_tuning || !sigma_wave_band_.isReady()) return noise_floor_sigma_;
        const float gain = sigma_wave_band_.whiteNoiseVarianceGain();
        if (!(std::isfinite(gain) && gain >= 0.0f)) return noise_floor_sigma_;
        return noise_floor_sigma_ * std::sqrt(gain);
    }

    [[nodiscard]] float wavePeriodFrequencyOrPrior_() const noexcept {
        const float f = wave_period_.getFrequencyHz();
        return (wave_period_.hasUsablePeriod() && std::isfinite(f) && f > 0.0f)
            ? f : kTuneFreqPriorHz;
    }

    [[nodiscard]] float pseudoUpdatePeriodFor_(float tau) const noexcept {
        if (!tau_scaled_pseudo_cadence_) return kPseudoPeriodNominalS;
        if (!(std::isfinite(tau) && tau > 0.0f)) return kPseudoPeriodNominalS;
        return std::min(std::max(kPseudoTauRatio * tau, kPseudoPeriodMinS),
                        kPseudoPeriodMaxS);
    }

    [[nodiscard]] float rsSpectralMseTarget_(float tau, float sigma) const noexcept {
        const float TS = pseudoUpdatePeriodFor_(tau);
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f)
                                  ? sigma_coeff_ : 1.0f;
        const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);
        if (!(TS > 0.0f) || !(tau > 0.0f)) return rs_mse_coeff_;
        const float tau2 = tau * tau;
        const float u = sigma_aB * tau2 * tau2;
        return rs_mse_coeff_ * rs_qeff_pow_
             * std::pow(u, 6.0f / 7.0f)
             / std::sqrt(TS);
    }

    [[nodiscard]] float rsTargetFromLaw_(float tau, float sigma) const noexcept {
        if (rs_law_ == RSLaw::LegacyCubic) {
            return R_S_coeff_ * sigma * tau * tau * tau;
        }
        return rsSpectralMseTarget_(tau, sigma);
    }

    void refreshQeffPow_() noexcept {
        float r_a = rs_accel_noise_density_;
        if (!(std::isfinite(r_a) && r_a > 0.0f)) r_a = R_S_ACCEL_NOISE_DENSITY_DEFAULT;
        rs_qeff_pow_ = std::pow(2.0f * r_a, 1.0f / 14.0f);
    }

    void updateTuner_(float dt, float a_up) {
        if (!enable_tuner_) return;

        const float f_hint = wavePeriodFrequencyOrPrior_();
        float f_band = tuner_.isFreqReady() ? tuner_.getFrequencyHz() : f_hint;
        if (!std::isfinite(f_band) || f_band < kMinTuneFreqHz) f_band = kMinTuneFreqHz;
        f_band = std::min(f_band, kMaxTuneFreqHz);

        const float a_for_variance = cfg_.wave_band_tuning
            ? sigma_wave_band_.step(a_up, dt, f_band)
            : a_up;
        tuner_.update(dt, a_for_variance, f_hint);
        if (fixed_tuning_) return;

        float f_tune = tuner_.getFrequencyHz();
        if (!std::isfinite(f_tune) || f_tune < kMinTuneFreqHz) f_tune = kMinTuneFreqHz;
        f_tune = std::min(f_tune, kMaxTuneFreqHz);
        sea_time_sec_ = 0.5f / f_tune;

        const float band_noise_sigma = bandNoiseFloorSigma_();
        const float var_noise = band_noise_sigma * band_noise_sigma;
        const float var_total = tuner_.isVarReady()
            ? std::max(0.0f, tuner_.getAccelVariance())
            : var_noise;
        const float sigma_wave = std::sqrt(std::max(1e-6f, var_total - var_noise));

        const float tau_live = std::min(max_tau_, std::max(min_tau_, tau_coeff_ * 0.5f / f_tune));
        float sigma_live = std::min(max_sigma_a_, sigma_coeff_ * sigma_wave);
        if (!tuner_.isVarReady()) {
            sigma_live = std::max(sigma_live, std::max(0.05f, band_noise_sigma));
        }

        // As in OU-III, r_S is derived from the live front-end estimate before
        // an ablation freezes either channel.  This keeps the partial-adaptation
        // arms orthogonal rather than silently freezing r_S with tau/sigma.
        const float rs_live = std::min(max_RS_, std::max(min_RS_,
                                                         rsTargetFromLaw_(tau_live, sigma_live)));

        if (!freeze_ou_channel_) {
            tau_target_ = tau_live;
            sigma_target_ = sigma_live;
        } else {
            tau_target_ = tune_.tau_applied;
            sigma_target_ = tune_.sigma_applied;
        }
        if (!freeze_RS_channel_) RS_target_ = rs_live;
    }

    void stagedColdStep_(const Vector3f& gyro, const Vector3f& acc, float dt) {
        // StagedMekf is a selectable bootstrap ablation; the deployed default
        // uses the MahonyProxy path below.
        const bool levelled = ::seastate::common::runStartupGravityInit(
            gyro, acc, dt, elapsed_sec_, cfg_.gravity_magnitude,
            0.5f, 1.0f, 0.12f, 0.5f, 0.5f, 8.0f, 0.15f,
            bootstrap_tilt_obs_, bootstrap_gravity_slow_lpf_,
            bootstrap_gravity_good_sec_,
            [&](const Eigen::Vector3f& g_body) { mekf_.initialize_from_acc(g_body); });
        if (!levelled) return;
        stage_ = StartupStage::TunerWarm;
        enterCold_();
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void tryProxyHandoff_() {
        if (elapsed_sec_ < cfg_.proxy_startup_min_sec) return;
        if (!vertical_complementary_.isInitialized()) return;

        const float mag_acquire_deadline = cfg_.with_mag
            ? cfg_.proxy_mag_settle_sec + 2.0f * std::max(cfg_.mag_min_window_sec, 1.0f)
                  + cfg_.mag_tilt_fallback_sec
            : 0.0f;
        const float timeout_sec = std::max(cfg_.proxy_startup_timeout_sec, mag_acquire_deadline);
        const bool timed_out = (elapsed_sec_ >= timeout_sec) && gravity_gate_.aligned_branch;

        if (!timed_out) {
            if (!vertical_complementary_.isReady()) return;
            if (!proxyGravityTrusted_()) return;
            if (!isTunerReady()) return;
            if (cfg_.with_mag && !mag_reference_learned_) return;
        }

        // The learned angle aligns the proxy's WORLD frame with magnetic
        // north; it is not an absolute boat heading to freeze at acquisition.
        // Carry every subsequent gyro-observed turn into the actual handoff.
        const Eigen::Quaternionf q_bw = have_mag_proxy_gauge_
            ? Eigen::Quaternionf(
                  Eigen::AngleAxisf(mag_proxy_yaw_offset_rad_, Vector3f::UnitZ()) *
                  vertical_complementary_.quaternion())
            : vertical_complementary_.tiltQuaternion();

        mekf_.initialize_from_truth(q_bw.normalized(), Vector3f::Zero(), Vector3f::Zero(),
                                    Vector3f::Zero(), Vector3f::Zero(),
                                    Vector3f::Zero(), Vector3f::Zero());
        seedHandoffAttitudeCovariance_();
        mekf_.set_initial_linear_uncertainty(1.0f, 2.0f, 5.0f, 1.0f);
        mekf_.set_initial_acc_bias_std(cfg_.handoff_acc_bias_std);
        if (mag_world_ref_valid_) mekf_.set_magnetic_reference_world(mag_world_ref_uT_);
        handoff_timed_out_ = timed_out;
        enterLive_();
    }

    [[nodiscard]] bool usingProxyInit_() const noexcept {
        return cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy;
    }

    void proxyUpdateMag_(const Vector3f& mag_body) {
        if (!mag_reference_learned_ && elapsed_sec_ < cfg_.proxy_mag_settle_sec) return;
        if (!vertical_complementary_.isInitialized()) return;
        if (!std::isfinite(mag_init_eligible_t0_)) mag_init_eligible_t0_ = elapsed_sec_;

        if (!mag_reference_learned_) {
            const bool fallback_ok = (elapsed_sec_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec;
            if (!proxyGravityTrusted_() && !fallback_ok) return;
            if (!have_last_imu_) return;
            learnMagReferenceWindowed_(mag_body);
        }

        // The gravity-only proxy cannot observe axial gyro bias. Keep its
        // pre-handoff yaw alignment attached to fresh magnetic observations,
        // rather than integrating that bias after a one-time north lock.
        if (mag_reference_learned_ && stage_ != StartupStage::Live) {
            Eigen::Quaternionf q_north;
            float gauge;
            if (northFrameForMag_(mag_body - mag_hard_iron_body_uT_, q_north, gauge)) {
                mag_proxy_yaw_offset_rad_ = ::seastate::common::wrapPi(-gauge);
            }
        }
        maybeRefineMagReference_(mag_body);
        maybeApplyContinuousHardIron_();
        if (mag_reference_learned_ && stage_ == StartupStage::Live) {
            mekf_.measurement_update_mag_only(mag_body - mag_hard_iron_body_uT_);
        }
    }

    void learnMagReferenceWindowed_(const Vector3f& mag_body) {
        const float dt_mag = ::seastate::common::advanceSampleClock(
            last_mag_sample_t_, elapsed_sec_, cfg_.mag_sample_dt_sec);

        Eigen::Quaternionf q_north;
        float sample_gauge;
        if (!northFrameForMag_(mag_body, q_north, sample_gauge)) return;
        // With no joint hard-iron fit, average field strength/dip in magnetic
        // north: each sample is (horizontal magnitude, 0, down). This frame
        // is invariant to both real turns and unobservable proxy yaw drift.
        // The optional joint bias solve retains its independent world frame.
        const Eigen::Quaternionf q_accum = cfg_.mag_estimate_hard_iron
            ? vertical_complementary_.quaternion() : q_north;
        if (!mag_auto_tuner_.addSampleWithWorldQuatDt(
                dt_mag, q_accum, last_acc_body_, last_gyro_body_, mag_body)) return;

        Vector3f ref;
        if (!mag_auto_tuner_.getMagWorldRef(ref) || !ref.allFinite() ||
            !(ref.norm() > cfg_.mag_init_min_mag_norm)) return;
        const float gauge = cfg_.mag_estimate_hard_iron
            ? mag_auto_tuner_.getYawGaugeCorrectionRad() : sample_gauge;
        if (!std::isfinite(gauge)) return; // no horizontal field, no north lock
        mag_proxy_yaw_offset_rad_ = ::seastate::common::wrapPi(-gauge);
        have_mag_proxy_gauge_ = true;

        // A missing magnetometer can let the existing timeout enter Live
        // before the first north lock. That first alignment is a frame choice:
        // rotate all world states/covariance, then install the north-frame
        // reference (apply_world_yaw_gauge also rotates any old reference).
        if (stage_ == StartupStage::Live) {
            const Eigen::Matrix3f R = mekf_.R_bw();
            const float yaw_core = std::atan2(R(1, 0), R(0, 0));
            mekf_.apply_world_yaw_gauge(::seastate::common::wrapPi(proxyMagneticYaw_(gauge) - yaw_core));
        }
        setMagWorldRef_(ref);

        Vector3f hard_iron;
        hard_iron_.startup_body_uT = mag_auto_tuner_.getHardIronBodyUT(hard_iron)
            ? hard_iron : Vector3f::Zero();
        mag_hard_iron_body_uT_ = hard_iron_.startup_body_uT + hard_iron_.applied_body_uT;
        mag_reference_learned_ = true;
        mag_north_lock_time_sec_ = elapsed_sec_;
    }

    void maybeRefineMagReference_(const Vector3f& mag_body) {
        if (!cfg_.mag_refine_enabled || mag_refine_done_) return;
        if (!mag_reference_learned_ || stage_ != StartupStage::Live) return;
        if (elapsed_sec_ < cfg_.mag_refine_start_sec || !have_last_imu_) return;

        if (!mag_refine_started_) {
            ::seastate::common::beginMagRefinement(
                mag_auto_tuner_, cfg_.mag_refine_window_sec, cfg_.mag_min_samples);
            mag_refine_started_ = true;
            last_mag_sample_t_ = NAN;
        }

        const float dt_mag = ::seastate::common::advanceSampleClock(
            last_mag_sample_t_, elapsed_sec_, cfg_.mag_sample_dt_sec);
        const Vector3f mag_corrected = mag_body - mag_hard_iron_body_uT_;
        Eigen::Quaternionf q_north;
        float sample_gauge;
        if (!northFrameForMag_(mag_corrected, q_north, sample_gauge)) return;
        const Eigen::Quaternionf q_accum = cfg_.mag_estimate_hard_iron
            ? vertical_complementary_.quaternion() : q_north;
        if (!mag_auto_tuner_.addSampleWithWorldQuatDt(
                dt_mag, q_accum, last_acc_body_, last_gyro_body_, mag_corrected)) return;

        Vector3f ref;
        if (!mag_auto_tuner_.getMagWorldRef(ref) || !ref.allFinite() ||
            !(ref.norm() > cfg_.mag_init_min_mag_norm)) return;
        const float gauge = cfg_.mag_estimate_hard_iron
            ? mag_auto_tuner_.getYawGaugeCorrectionRad() : sample_gauge;
        if (!std::isfinite(gauge)) return;
        setMagWorldRef_(ref);
        // Use the current magnetic heading, retaining the core's wave-aware
        // tilt. A window-mean proxy gauge would lag by gyro_bias * window/2.
        mekf_.set_attitude_yaw_absolute(proxyMagneticYaw_(gauge));
        mag_refine_done_ = true;
        mag_refine_time_sec_ = elapsed_sec_;
        maybeUnlockAccBias_();
    }

    // Sensor-only north frame. No main-filter state or wave-amplitude gate
    // participates; vertical/invalid fields cannot manufacture a north lock.
    bool northFrameForMag_(const Vector3f& mag_body, Eigen::Quaternionf& q_north,
                           float& gauge) const {
        if (!mag_body.allFinite()) return false;
        const Eigen::Quaternionf q_proxy = vertical_complementary_.quaternion();
        const Vector3f field = q_proxy * mag_body;
        const float norm = field.norm();
        const float horizontal = field.head<2>().norm();
        const auto& config = mag_auto_tuner_.config();
        if (!(norm > cfg_.mag_init_min_mag_norm) || !std::isfinite(horizontal) ||
            !(horizontal > config.min_horizontal_fraction * norm)) return false;
        gauge = std::atan2(field.y(), field.x());
        q_north = Eigen::AngleAxisf(-gauge, Vector3f::UnitZ()) * q_proxy;
        q_north.normalize();
        return q_north.coeffs().allFinite();
    }

    [[nodiscard]] float proxyMagneticYaw_(float gauge) const {
        const Eigen::Matrix3f R = vertical_complementary_.quaternion().toRotationMatrix();
        return ::seastate::common::wrapPi(std::atan2(R(1, 0), R(0, 0)) - gauge);
    }

    void setMagWorldRef_(const Vector3f& ref) {
        mag_world_ref_uT_ = ref;
        mag_world_ref_valid_ = true;
        if (stage_ == StartupStage::Live || !usingProxyInit_()) {
            mekf_.set_magnetic_reference_world(ref);
        }
    }

    // Raw magnetometer in the proxy's tilt frame; see
    // ::seastate::common::ContinuousHardIronTracker::accumulate().
    void accumulateContinuousHardIron_(const Vector3f& mag_body) {
        if (!cfg_.mag_continuous_hard_iron || !usingProxyInit_()) return;
        if (!vertical_complementary_.isInitialized()) return;
        hard_iron_.accumulate(elapsed_sec_, cfg_.mag_sample_dt_sec,
                              vertical_complementary_.tiltQuaternion(), mag_body);
    }

    void maybeApplyContinuousHardIron_() {
        if (!cfg_.mag_continuous_hard_iron || !usingProxyInit_()) return;
        if (!mag_reference_learned_ || stage_ != StartupStage::Live) return;
        if (cfg_.mag_refine_enabled && !mag_refine_done_) return;
        if (!mag_world_ref_valid_) return;

        Vector3f ref;
        if (!hard_iron_.slewTowardEstimate(elapsed_sec_, cfg_.mag_sample_dt_sec,
                                           mag_world_ref_uT_,
                                           cfg_.mag_hi_apply_fraction,
                                           cfg_.mag_hi_slew_tau_sec,
                                           cfg_.mag_init_min_mag_norm,
                                           ref)) return;

        mag_hard_iron_body_uT_ = hard_iron_.startup_body_uT + hard_iron_.applied_body_uT;
        setMagWorldRef_(ref);
    }

    void updateProxyGravityQuality_(float dt, const Vector3f& gyro, const Vector3f& acc) {
        if (!vertical_complementary_.isInitialized()) {
            gravity_gate_.good_sec = 0.0f;
            gravity_gate_.aligned_branch = false;
            return;
        }

        // As in the OU front ends: rotate into the attitude's world frame
        // first, then average over the wave band.  Unlike them, the attitude
        // judged is always the proxy's, before and after handoff.
        gravity_gate_.step(vertical_complementary_.quaternion(), acc, gyro, dt,
                           cfg_.proxy_gravity_lpf_sec,
                           cfg_.proxy_gravity_warmup_sec,
                           cfg_.proxy_gravity_align_sin,
                           cfg_.mag_extreme_gyro_dps);
    }

    [[nodiscard]] bool proxyGravityTrusted_() const {
        return gravity_gate_.trusted(cfg_.proxy_gravity_hold_sec);
    }

    void seedHandoffAttitudeCovariance_() {
        auto& P = mekf_.covariance_full();
        const float st = std::max(1e-6f, cfg_.proxy_handoff_tilt_sigma_rad);
        const float sy = have_mag_proxy_gauge_
            ? std::max(1e-6f, cfg_.proxy_handoff_yaw_sigma_rad)
            : std::max(1e-6f, cfg_.proxy_handoff_yaw_sigma_free_rad);
        constexpr int PHI = Mekf::OFF_PHI;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < static_cast<int>(Mekf::NX); ++j) {
                P(PHI + i, j) = 0.0f;
                P(j, PHI + i) = 0.0f;
            }
        }
        P(PHI + 0, PHI + 0) = st * st;
        P(PHI + 1, PHI + 1) = st * st;
        P(PHI + 2, PHI + 2) = sy * sy;
    }

    // TFG's own adaptation schedule: channel targets are formed from the
    // first tuner sample (there is no Cold warmup gate on the targets), EMAs
    // run only once Live, and the smoothed candidate is committed when an
    // accumulated-time counter reaches adapt_every_secs_ -- where the OU
    // wrappers compare wall time against the last commit.  Both stage the
    // commit to the start of the next sample.
    void adaptMekf_(float dt) {
        if (!enable_tuner_ || fixed_tuning_) return;

        if (!freeze_ou_channel_) {
            float adapt_sec = adapt_tau_sec_;
            if (adapt_tau_sea_periods_ > 0.0f &&
                std::isfinite(sea_time_sec_) && sea_time_sec_ > 0.0f) {
                const float safe_sea_time =
                    seastate::tuner::limits::clampDynamicEmaTimeScaleSec(sea_time_sec_);
                adapt_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(
                    adapt_tau_sea_periods_ * safe_sea_time, dt);
            }
            const float a = 1.0f - std::exp(-dt / adapt_sec);
            tune_.tau_applied += a * (tau_target_ - tune_.tau_applied);
            tune_.sigma_applied += a * (sigma_target_ - tune_.sigma_applied);
        }
        if (!freeze_RS_channel_) {
            const float horizon = ::seastate::common::adaptiveSmoothingHorizonSec(
                adapt_RS_mult_, tau_target_, RS_target_, tune_.RS_applied,
                adapt_RS_slew_log_, dt);
            const float a = (horizon > 0.0f) ? (1.0f - std::exp(-dt / horizon)) : 1.0f;
            tune_.RS_applied += a * (RS_target_ - tune_.RS_applied);
        }

        adapt_elapsed_sec_ += dt;
        if (adapt_elapsed_sec_ >= adapt_every_secs_) {
            adapt_elapsed_sec_ = 0.0f;
            tune_apply_pending_ = true;
        }
    }

    void applyPendingTune_() {
        if (!tune_apply_pending_) return;
        tune_apply_pending_ = false;
        commitTune_();
    }

    void periodicAwCovSyncTick_(float dt) {
        if (!periodic_aw_cov_sync_ || stage_ != StartupStage::Live) return;
        aw_sync_elapsed_sec_ += dt;
        if (aw_sync_elapsed_sec_ < adapt_every_secs_) return;
        aw_sync_elapsed_sec_ = 0.0f;
        mekf_.synchronize_aw_covariance_to_stationary();
    }

    void applyPseudoCadence_() {
        pseudo_period_sec_ = pseudoUpdatePeriodFor_(tune_.tau_applied);
    }

    void commitTune_() {
        mekf_.set_aw_time_constant(tune_.tau_applied);
        applyPseudoCadence_();
        const float sZ = tune_.sigma_applied;
        mekf_.set_aw_stationary_std(Vector3f(sZ * S_factor_, sZ * S_factor_, sZ));

        // SpectralMSE already contains the realized target T_S.  The
        // LegacyCubic target does not, so only that law receives the
        // information-rate normalization here.
        float rs = tune_.RS_applied;
        if (rs_law_ == RSLaw::LegacyCubic && tau_scaled_pseudo_cadence_ && pseudo_period_sec_ > 0.0f) {
            rs *= std::sqrt(kPseudoPeriodNominalS / pseudo_period_sec_);
        }
        RS_filter_input_ = rs;
        mekf_.set_RS_noise(Vector3f(rs * R_S_x_factor_, rs * R_S_y_factor_, rs));
    }

    using Vec3LPF = ::seastate::common::Vec3LPF;

    static constexpr float kTuneFreqPriorHz = defaults::TUNE_FREQ_PRIOR_HZ;
    static constexpr float kMinTuneFreqHz = defaults::MIN_TUNE_FREQ_HZ;
    // TFG-specific: OU-II also uses 1.5 Hz; OU-III uses 1.2 Hz.
    static constexpr float kMaxTuneFreqHz = 1.5f;
    static constexpr float kPseudoPeriodNominalS = defaults::PSEUDO_UPDATE_PERIOD_NOMINAL_S;
    static constexpr float kPseudoTauNominalS = defaults::PSEUDO_UPDATE_TAU_NOMINAL_S;
    static constexpr float kPseudoTauRatio = kPseudoPeriodNominalS / kPseudoTauNominalS;
    static constexpr float kPseudoPeriodMinS = defaults::PSEUDO_UPDATE_PERIOD_MIN_S;
    // TFG-specific upper cadence clamp (OU-II 250 ms, OU-III 150 ms).
    static constexpr float kPseudoPeriodMaxS = 0.25f;

    Config cfg_{};
    Mekf mekf_{};
    StartupStage stage_{StartupStage::Cold};
    ::SeaStateAutoTuner tuner_{TUNER_SIGMA_VAR_K_PERIODS_DEFAULT, 1.0f};
    ::WavePeriodEstimator wave_period_{};
    ::VerticalAccelComplementary vertical_complementary_{};
    ::AdaptiveWaveBandPass sigma_wave_band_{SIGMA_BAND_LOW_RATIO_DEFAULT,
                                            SIGMA_BAND_HIGH_RATIO_DEFAULT,
                                            SIGMA_BAND_MIN_HZ_DEFAULT,
                                            SIGMA_BAND_MAX_HZ_DEFAULT};
    ::MagAutoTuner mag_auto_tuner_{};
    ::seastate::common::ContinuousHardIronTracker hard_iron_{};
    ::seastate::common::StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF bootstrap_gravity_slow_lpf_{};
    float bootstrap_gravity_good_sec_ = 0.0f;

    TfgTuneState tune_{};
    float tau_target_ = 1.1f;
    float sigma_target_ = 1e-2f;
    float RS_target_ = 0.5f;

    // TFG-specific OU prior and regularizer coefficients.  R_S_coeff applies
    // only to LegacyCubic; SpectralMSE uses rs_mse_coeff_.  The horizontal
    // factors are independent: X = 1.08, Y = 1.15, with an isotropic a_w prior.
    float tau_coeff_ = 1.0f;
    float sigma_coeff_ = 0.8f;
    float R_S_coeff_ = 0.28f;
    float S_factor_ = 1.00f;
    float R_S_x_factor_ = 1.08f;
    float R_S_y_factor_ = 1.15f;
    float noise_floor_sigma_ = defaults::ACC_NOISE_FLOOR_SIGMA;

    RSLaw rs_law_ = RSLaw::SpectralMSE;
    float rs_accel_noise_density_ = R_S_ACCEL_NOISE_DENSITY_DEFAULT;
    float rs_mse_coeff_ = R_S_MSE_COEFF_DEFAULT;
    float rs_qeff_pow_ = std::pow(2.0f * R_S_ACCEL_NOISE_DENSITY_DEFAULT, 1.0f / 14.0f);

    float adapt_tau_sec_ = ADAPT_TAU_SEC_DEFAULT;
    float adapt_tau_sea_periods_ = ADAPT_TAU_SEA_PERIODS_DEFAULT;
    float sea_time_sec_ = 0.5f / kTuneFreqPriorHz;
    float adapt_RS_mult_ = ADAPT_RS_MULT_DEFAULT;
    float adapt_RS_slew_log_ = 0.0f;

    float min_tau_ = 0.02f, max_tau_ = 12.0f;
    float min_RS_ = 0.15f, max_RS_ = 400.0f;
    float max_sigma_a_ = 6.0f;
    bool enable_tuner_ = true;
    bool fixed_tuning_ = false;
    bool freeze_ou_channel_ = false;
    bool freeze_RS_channel_ = false;
    Vector3f Racc_nominal_{Vector3f::Constant(0.5f)};
    bool     warmup_Racc_active_ = false;

    float elapsed_sec_ = 0.0f;
    float live_sec_ = 0.0f;
    float mag_elapsed_sec_ = 0.0f;
    float pseudo_elapsed_ = 0.0f;
    float pseudo_period_sec_ = kPseudoPeriodNominalS;
    float adapt_every_secs_ = defaults::ADAPT_EVERY_SECS;
    float adapt_elapsed_sec_ = 0.0f;
    bool tune_apply_pending_ = false;
    bool tau_scaled_pseudo_cadence_ = true;
    float RS_filter_input_ = 0.5f;

    float tuner_warm_sec_ = 0.0f;
    bool mag_reference_learned_ = false;
    bool handoff_timed_out_ = false;
    bool acc_bias_unlocked_ = false;
    bool acc_bias_hold_ = false;
    float live_since_sec_ = 0.0f;
    ::seastate::common::WorldGravityGate gravity_gate_{};

    bool periodic_aw_cov_sync_ = true;
    float aw_sync_elapsed_sec_ = 0.0f;

    float mag_proxy_yaw_offset_rad_ = 0.0f;
    bool have_mag_proxy_gauge_ = false;
    Vector3f mag_world_ref_uT_{Vector3f::Zero()};
    bool mag_world_ref_valid_ = false;
    float mag_init_eligible_t0_ = NAN;
    float last_mag_sample_t_ = NAN;
    float mag_north_lock_time_sec_ = NAN;
    bool mag_refine_started_ = false;
    bool mag_refine_done_ = false;
    float mag_refine_time_sec_ = NAN;

    Vector3f mag_hard_iron_body_uT_{Vector3f::Zero()};

    Vector3f last_acc_body_{Vector3f::Zero()};
    Vector3f last_gyro_body_{Vector3f::Zero()};
    bool have_last_imu_ = false;
};

}  // namespace ocean_imu::tfg
