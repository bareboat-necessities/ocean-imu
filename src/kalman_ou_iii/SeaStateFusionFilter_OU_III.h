#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  SeaStateFusionFilter_OU_III

  Marine Inertial Navigational System (INS) Filter for IMU

  Combines multiple real-time estimators into a cohesive ocean-state tracker:

    • Quaternion-based attitude and linear motion estimation via Kalman3D_Wave_OU_III

    • Dominant frequency tracking using one of:
          – AranovskiyFreqTracker     (frequency estimator)
          – KalmANFFreqTracker        (adaptive notch / Kalman frequency tracker)
          - PLLFreqTracker            (PLL frequency tracker)
          – SchmittTrigger            (zero-cross event detector)

    • Dual-stage frequency smoothing:
          – Fast 1st-order IIR (≈ few s, ~90% step) for demodulation / direction
          – Slow 1st-order IIR (≈ longer s, ~90% step) for auto-tuning / moments

    • Online auto-tuning of Kalman filter parameters (τ, σₐ, Rₛ) through
      SeaStateAutoTuner, which estimates acceleration variance and applies the
      σₐ·τ³ regularization law to stabilize displacement drift correction.

  Where
  – τ (tau):  OU process time constant ≈ ½ · T  (half the dominant period of acceleration)
  – σₐ:       Stationary acceleration standard deviation, EWMA-tracked online
  – Rₛ:       Pseudo-measurement noise controlling integral drift suppression
  – Rₛ_xy:    Reduced in X/Y (anisotropic weighting for vertical-dominant seas)

  Adaptive update:  exponential smoothing toward targets over ADAPT_TAU_SEC

  Features
  • Modular tracker selection via TrackerPolicy template
  • Quaternion-consistent Euler conversion (aerospace → nautical, ENU frame)
  • Magnetometer yaw correction with configurable startup delay
  • Fully compatible with Arduino or native Eigen builds
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <numbers>
#include <memory>
#include <algorithm>

#include "freq/FirstOrderIIRSmoother.h"
#include "freq/FrequencyTrackerPolicy.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/WavePeriodEstimator.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/MagAutoTuner.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "wave_dir/KalmanWaveDirection.h"
#include "wave_dir/WaveDirectionDetector.h"
#include "wave_dir/WaveDirectionFrame.h"
#include "detrend/AdaptiveWaveDetrender3D.h"
#include "kalman_common/SeaStateFusionFilterCommon.h"

// Shared constants
extern const float g_std;

#ifndef FREQ_GUESS
#define FREQ_GUESS 0.3f
#endif

#ifndef ZERO_CROSSINGS_SCALE
#define ZERO_CROSSINGS_SCALE 1.0f
#endif

#ifndef ZERO_CROSSINGS_DEBOUNCE_TIME
#define ZERO_CROSSINGS_DEBOUNCE_TIME 0.12f
#endif

#ifndef ZERO_CROSSINGS_STEEPNESS_TIME
#define ZERO_CROSSINGS_STEEPNESS_TIME 0.21f
#endif

// Estimated vertical accel noise floor (1σ), m/s².
// Tweak from bench data with IMU sitting still.
constexpr float ACC_NOISE_FLOOR_SIGMA_DEFAULT = 0.12f;

constexpr float MIN_FREQ_HZ = 0.2f;
constexpr float MAX_FREQ_HZ = 6.0f;

// Floor for the wave-band tuning frequency.  MIN_FREQ_HZ bounds the
// acceleration-band tracker and is far too high for a zero-crossing period:
// 0.03 Hz admits a 33 s swell.
constexpr float MIN_TUNE_FREQ_HZ = 0.03f;

constexpr float MIN_TAU_S   = 0.02f;
// tau now scales with the zero-crossing wave period rather than with an
// acceleration-band frequency, so the ceiling has to admit a developed sea:
// T_z reaches 8.6 s at H_s = 8.5 m and a long swell goes further.  The old
// 3.0 s ceiling was reached at H_s = 8.5 m, which clipped the operating point
// exactly where the filter was losing.
constexpr float MAX_TAU_S   = 12.0f;
constexpr float MAX_SIGMA_A = 6.0f;
constexpr float MIN_R_S     = 0.4f;
// r_S ~ sigma_aw * tau^3 inherits that range.  The old 35 m*s ceiling was the
// binding constraint at H_s = 8.5 m: the calibrated fixed-oracle point sat at
// 34.66 and the error was still falling monotonically against it.
constexpr float MAX_R_S     = 400.0f;

constexpr float ADAPT_TAU_SEC              = 1.8f;
constexpr float ADAPT_EVERY_SECS           = 0.1f;
// Smoothing horizon of the r_S channel, in units of tau_target.  Measured on
// the versioned records against synthesized sea-state transitions: the error
// during a transition falls monotonically as this shortens, the stationary
// worst-record vertical error rises monotonically, and 3.0 is where the two
// cross at an acceptable cost.  r_S ~ sigma_aw tau^3 amplifies tau noise by
// the third power, which is why this sits above OU-II's tau^2 and tau^1
// channels.  See docs/ou-ema-adaptation-tuning.md.
constexpr float ADAPT_RS_MULT              = 3.0f;   // dimensionless
// Discrepancy, in natural-log units of the r_S target-to-applied ratio, above
// which the smoothing horizon shortens.  Zero keeps the plain proportional
// horizon; see docs/ou-ema-adaptation-tuning.md.
constexpr float ADAPT_RS_SLEW_LOG          = 0.0f;   // ln units
constexpr float ONLINE_TUNE_WARMUP_SEC     = 5.0f;
constexpr float MAG_DELAY_SEC              = 7.0f;

// Frequency smoother dt (SeaStateFusionFilter_OU_III is designed for 200 Hz)
constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;

struct TuneState {
    float tau_applied   = 1.1f;    // s
    float sigma_applied = 1e-2f;   // m/s²
    float RS_applied    = 0.5f;    // m*s
};

//  Unified SeaState fusion filter
template<TrackerType trackerT>
class SeaStateFusionFilter_OU_III {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using TrackingPolicy = TrackerPolicy<trackerT>;

    enum class StartupStage {
        Cold,        // just booted or just had a big tilt reset
        TunerWarm,   // MEKF + freq running, tuner collecting stats
        Live         // tuner is trusted; full adaptation & extras allowed
    };

    explicit SeaStateFusionFilter_OU_III(bool with_mag = true)
        : with_mag_(with_mag),
          time_(0.0),
          last_adapt_time_sec_(0.0),
          freq_hz_(FREQ_GUESS),
          freq_hz_slow_(FREQ_GUESS)
    {
        // Default cutoff ~max_freq_hz_ Hz: passes waves, kills 8–37 Hz engine band
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;
    }

    StartupStage getStartupStage() const noexcept { return startup_stage_; }
    bool isAdaptiveLive() const noexcept { return startup_stage_ == StartupStage::Live; }

    void initialize(const Eigen::Vector3f& sigma_a,
                    const Eigen::Vector3f& sigma_g,
                    const Eigen::Vector3f& sigma_m)
    {
        mekf_ = std::make_unique<Kalman3D_Wave_OU_III<float>>(sigma_a, sigma_g, sigma_m);
        seastate::common::finalizeInitialization(
            mekf_,
            [this]() { enterCold_(); },
            [this]() { apply_ou_tune_(true); });
    }

    void initialize_ext(const Eigen::Vector3f& sigma_a,
                        const Eigen::Vector3f& sigma_g,
                        const Eigen::Vector3f& sigma_m,
                        float Pq0, float Pb0,
                        float b0, float R_S_noise,
                        float gravity_magnitude)
    {
        mekf_ = std::make_unique<Kalman3D_Wave_OU_III<float>>(sigma_a, sigma_g, sigma_m, Pq0, Pb0, b0, R_S_noise, gravity_magnitude);
        seastate::common::finalizeInitialization(
            mekf_,
            [this]() { enterCold_(); },
            [this]() { apply_ou_tune_(true); });
    }

    void initialize_from_acc(const Eigen::Vector3f& acc_body_ned) {
        if (mekf_) {
            mekf_->initialize_from_acc(acc_body_ned);
        }
    }

    // Time update (IMU integration + frequency tracking)
    void updateTime(float dt, const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc,
                    float tempC = 35.0f)
    {
        if (!mekf_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;
        time_ += dt;
        startup_stage_t_ += dt;


        // BODY-Z-based proxy used by the tracker/tuner logic.
        // This is NOT a true vertical acceleration estimate; it is only a
        // body-Z residual that behaves like up-positive vertical motion when the
        // platform is near-level:
        //   acc.z() ~ -g at rest  => proxy ~ 0
        const float a_z_body_proxy = acc.z() + g_std;

        // Private Mahony observer for the wave-period estimator.  It is fed the
        // raw gyro and accelerometer, before the MEKF sees them, so that the
        // levelling it provides stays a pure function of the measurements.
        // Stepping it unconditionally keeps its transient off the critical path
        // when the input source is switched at runtime.
        vertical_accel_comp_.update(dt, gyro, acc, g_std);

        // MEKF updates first (attitude + latent a_w)
        mekf_->time_update(gyro, dt);
        mekf_->measurement_update_acc_only(acc, tempC);

        {
            Eigen::Quaternionf q_bw = mekf_->quaternion_boat();
            q_bw.normalize();

            const Eigen::Vector3f z_body_down_world = q_bw * Eigen::Vector3f(0.0f, 0.0f, 1.0f);
            const Eigen::Vector3f z_world_down(0.0f, 0.0f, 1.0f);

            float cos_tilt = z_body_down_world.normalized().dot(z_world_down);
            cos_tilt = std::max(-1.0f, std::min(1.0f, cos_tilt));
            const float tilt_deg = std::acos(cos_tilt) * 57.295779513f;

            constexpr float TILT_RESET_DEG = 70.0f;
            constexpr float TILT_RESET_HOLD_SEC = 0.35f;
            constexpr float TILT_RESET_COOLDOWN_SEC = 3.0f;

            if (tilt_reset_cooldown_sec_ > 0.0f) {
                tilt_reset_cooldown_sec_ = std::max(0.0f, tilt_reset_cooldown_sec_ - dt);
            }

            if (tilt_deg > TILT_RESET_DEG) {
                tilt_over_limit_sec_ += dt;
            } else {
                // decay quickly on recovery so brief transients do not trigger resets
                tilt_over_limit_sec_ = std::max(0.0f, tilt_over_limit_sec_ - 2.0f * dt);
            }

            if (tilt_over_limit_sec_ >= TILT_RESET_HOLD_SEC && tilt_reset_cooldown_sec_ <= 0.0f) {
                if (startup_stage_ == StartupStage::Live) {
                    // In Live, re-lock only tilt while preserving yaw/north frame.
                    mekf_->initialize_from_acc_preserve_yaw(acc);
                } else {
                    // During startup stages, accel-only re-lock is acceptable.
                    mekf_->initialize_from_acc(acc);
                    enterCold_();
                    resetTrackingState_();
                }

                tilt_over_limit_sec_ = 0.0f;
                tilt_reset_cooldown_sec_ = TILT_RESET_COOLDOWN_SEC;
            }
        }

        // Up-positive BODY-Z proxy used by tracker/tuner logic.
        // Not true world vertical unless the platform is close to level.
        a_body_z_up_proxy_ = -a_z_body_proxy;

        // Vertical acceleration the tracker runs on.  Both choices are
        // measurement-only; see FreqTrackerInputSource for why levelling is
        // not the obvious win here that it is for the period estimator.
        const float a_vert_for_tracker =
            (freq_tracker_input_ == FreqTrackerInputSource::Complementary &&
             vertical_accel_comp_.isReady())
                ? vertical_accel_comp_.verticalAccelUpMs2()
                : a_body_z_up_proxy_;

        // LPF on the tracker input
        const float a_vert_lp = freq_input_lpf_.step(a_vert_for_tracker, dt);

        // Raw freq from tracker
        const float f_tracker = static_cast<float>(tracker_policy_.run(a_vert_lp, dt));
        f_raw = f_tracker;

        // Stillness detector shares the tracker's input, as it always has.
        const float f_after_still = freq_stillness_.step(a_vert_lp, dt, f_tracker);

        // Fast & slow smoothed frequencies
        float f_fast = freq_fast_smoother_.update(f_after_still);
        float f_slow = freq_slow_smoother_.update(f_fast);

        f_fast = std::min(std::max(f_fast, min_freq_hz_), max_freq_hz_);
        f_slow = std::min(std::max(f_slow, min_freq_hz_), max_freq_hz_);

        freq_hz_      = f_fast;   // demod / direction
        freq_hz_slow_ = f_slow;   // tuner / moments

        // Tuner gets vertical accel, and the wave-band frequency when the
        // period estimator has settled.  That single substitution fixes both
        // halves of the old operating point: tau stops being derived from an
        // acceleration-band frequency that barely moves with the sea state, and
        // the tuner's variance horizon (a few periods) stops being shorter than
        // one wave period, which was biasing sigma_aw low.
        if (enable_tuner_) {
            update_tuner(dt, a_body_z_up_proxy_, tuner_frequency_hz_(f_after_still));
        }

        // Keep linear-block R_S tuning responsive in Live mode instead of
        // waiting for slow adaptation cadence.
        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_RS_tune_();
        }

        // Bounded covariance inflation of the a_w marginal.
        periodic_aw_cov_sync_tick_();

        const float omega = 2.0f * static_cast<float>(M_PI) * freq_hz_;

        // Resolve direction in a leveled frame aligned with boat heading.
        // This removes roll/pitch mixing while preserving 0 deg = bow and
        // positive angles toward starboard.  The existing body-Z proxy remains
        // the tuner/tracker input; direction uses coherent leveled components.
        const auto direction_accel = wave_direction::heading_frame_acceleration<float>(
            mekf_->quaternion_boat(), acc, g_std);

        // Stage 1 estimates the apparent propagation plane as an unsigned axis
        // relative to boat heading.  Stage 2 resolves propagation sense along
        // that same axis from horizontal/vertical orbital phase.
        // Zero-crossing wave period.  This runs beside the frequency tracker
        // rather than replacing it: the tracker supplies the acceleration-band
        // carrier the direction demodulator needs, while the OU operating point
        // needs the wave band.
        //
        // The input must be levelled and must not read estimator state, and
        // those two requirements pulled against each other for a while.
        //
        // Levelled, because double integration weights a spectrum by
        // 1/omega^4, so the sub-band gravity leakage a tilting platform puts
        // into the body-Z proxy dominates the elevation proxy.  Fed the raw
        // body-Z proxy the estimator reports 6.8-10.0 s whatever the sea does,
        // against a truth of 2.4-8.7 s, and vertical RMS degrades 2.5x.
        //
        // Exogenous, because levelling with the filter's own attitude closes a
        // loop: a 0.25 rad attitude displacement moved the reported period
        // 8.05 -> 10.3 s and tau by 1.28x, and it reached the linear block too,
        // since displacing v, p, S or a_w perturbs attitude through the
        // filter's cross-covariances.  The stability appendix carried that as
        // its open interconnection.
        //
        // VerticalAccelComplementary satisfies both: it levels, but with a
        // private Mahony observer reading only the raw gyro and accelerometer,
        // so it is a pure function of the measurements.  It costs nothing --
        // over the eight reference records it matches the old attitude-levelled
        // input to within 0.2% of vertical RMS -- and it is the default.
        // setWavePeriodInput() still reaches the other two for ablation;
        // tests/kalman_ou_iii/tuner_coupling-test.cpp asserts the default is
        // exogenous bit-for-bit and bounds the Leveled path's gain.
        wave_period_.update(dt, wave_period_input_ms2_(direction_accel));

        dir_filter_.update(direction_accel.forward_ms2,
                           direction_accel.starboard_ms2,
                           omega, dt);
        const Eigen::Vector2f propagation_axis_boat = dir_filter_.getAxis();
        dir_sign_state_ = dir_sign_.update(
            direction_accel.forward_ms2,
            direction_accel.starboard_ms2,
            direction_accel.up_ms2,
            propagation_axis_boat.x(), propagation_axis_boat.y(),
            dt, dir_filter_.getLastStableConfidence());
    }

    //  Magnetometer correction
    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!with_mag_ || !mekf_) return;
        if (time_ < mag_delay_sec_) return;

        mekf_->measurement_update_mag_only(mag_body_ned);
        mag_updates_applied_++;

        if (!std::isfinite(first_mag_update_time_)) {
            first_mag_update_time_ = static_cast<float>(time_);
        }

        // We can "unlock" once mag has had a few updates, but we DO NOT
        // enable accel-bias learning or restore Racc unless we're already Live.
        if (accel_bias_locked_ &&
            startup_stage_ == StartupStage::Live &&
            mag_updates_applied_ >= MAG_UPDATES_TO_UNLOCK &&
            std::isfinite(first_mag_update_time_) &&
            (static_cast<float>(time_) - first_mag_update_time_) > 1.0f) // 1s guard
        {
            accel_bias_locked_ = false;

            // Only allow accel bias to start learning once the system is Live.
            if (freeze_acc_bias_until_live_ && startup_stage_ == StartupStage::Live) {
                mekf_->set_acc_bias_updates_enabled(true);

                // Restore nominal Racc only when bias learning is allowed.
                if (warmup_Racc_active_) {
                    if (Racc_nominal_.allFinite() && Racc_nominal_.maxCoeff() > 0.0f) {
                        mekf_->set_Racc_std(Racc_nominal_);
                        warmup_Racc_active_ = false;
                    }
                }
            }
        }
    }

    void setWithMag(bool with_mag) {
        with_mag_ = with_mag;
    }

    // Anisotropy configuration (runtime)
    // S-factor scales horizontal vs vertical stationary std of a_w.
    // RS XY factor scales pseudo-measurement noise in X/Y vs Z.
    void setSFactor(float s) {
        if (std::isfinite(s) && s > 0.0f) {
            S_factor_ = s;
        }
    }
    void setRSXYFactor(float k) {
        if (std::isfinite(k)) {
            R_S_xy_factor_ = std::min(std::max(k, 0.0f), 1.0f);
        }
    }

    void setTauCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            tau_coeff_ = c;
        }
    }
    void setSigmaCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            sigma_coeff_ = c;
        }
    }
    void setRSCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_S_coeff_;
            R_S_coeff_ = c;

            if (std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.RS_applied) && tune_.RS_applied > 0.0f) {
                    tune_.RS_applied *= scale;
                }
                if (std::isfinite(RS_target_) && RS_target_ > 0.0f) {
                    RS_target_ *= scale;
                }

                if (enable_linear_block_) {
                    apply_RS_tune_();
                }
            }
        }
    }

    void setAccNoiseFloorSigma(float s) {
        if (std::isfinite(s) && s > 0.0f) {
            acc_noise_floor_sigma_ = s;
        }
    }
    float getAccNoiseFloorSigma() const noexcept {
        return acc_noise_floor_sigma_;
    }

    // Configure LPF on BODY-Z proxy for tracker input
    void setFreqInputCutoffHz(float fc) {
        freq_input_lpf_.setCutoff(fc);
    }

    void enableClamp(bool flag = true) {
        enable_clamp_ = flag;
    }
    void enableTuner(bool flag = true) {
        enable_tuner_ = flag;
    }
    bool tunerEnabled() const noexcept { return enable_tuner_; }

    // Policy for the latent-acceleration marginal P_{a_w a_w}.
    //
    // Default (true): once per adaptation period the marginal is re-aligned
    // with the stationary OU covariance, keeping the cross-covariances the
    // filter has learned.  This is a deliberate bounded covariance inflation.
    // It stops the a_w marginal from settling far below the level the process
    // model considers stationary, which keeps the accelerometer gain
    // responsive when the sea state changes.  It is not free -- it discards
    // posterior information at the adaptation cadence -- so the alternative
    // is available and measured rather than assumed.
    //
    // With false, the marginal is aligned only at discrete reconfiguration
    // events (construction, the transition to Live, setFixedTuning()) and a
    // changed stationary scale reaches the filter solely through the discrete
    // OU process covariance.
    //
    // The policy is applied independently of the tuner so that fixed-tuning
    // modes run it too.  Otherwise an adaptive-versus-fixed comparison would
    // confound whether the parameters adapt with whether part of the
    // covariance is periodically re-aligned.
    void setPeriodicAwCovarianceSync(bool flag) {
        periodic_aw_cov_sync_ = flag;
        last_aw_cov_sync_sec_ = time_;
    }
    bool periodicAwCovarianceSync() const noexcept { return periodic_aw_cov_sync_; }

    // Select how the a_w marginal is re-aligned when a sync happens.
    //
    // false (default, deployed): overwrite the marginal and keep the raw
    // cross-covariances, which rescales the implied correlations by the square
    // root of the marginal change.
    // true: congruence re-alignment, which reaches the same marginal while
    // leaving the whitened cross-covariance untouched and staying PSD by
    // construction.  This is the consistent operation, so running it isolates
    // "is the re-alignment inconsistent?" from "is re-aligning at all a good
    // idea?".
    void setAwCovarianceSyncCongruent(bool flag) { congruent_aw_cov_sync_ = flag; }
    bool awCovarianceSyncCongruent() const noexcept { return congruent_aw_cov_sync_; }

    // Freeze the online tuner at an externally supplied operating point. This
    // is primarily useful for controlled ablations (fixed-nominal and
    // fixed-oracle) after the normal startup sequence has reached Live.
    bool setFixedTuning(float tau_s, float sigma_a, float RS)
    {
        if (!(std::isfinite(tau_s) && tau_s > 0.0f &&
              std::isfinite(sigma_a) && sigma_a > 0.0f &&
              std::isfinite(RS) && RS > 0.0f))
        {
            return false;
        }

        enable_tuner_ = false;
        tau_target_ = enable_clamp_
            ? std::min(std::max(tau_s, min_tau_s_), max_tau_s_)
            : tau_s;
        sigma_target_ = enable_clamp_
            ? std::min(sigma_a, max_sigma_a_)
            : sigma_a;
        RS_target_ = enable_clamp_
            ? std::min(std::max(RS, min_R_S_), max_R_S_)
            : RS;

        tune_.tau_applied = tau_target_;
        tune_.sigma_applied = sigma_target_;
        tune_.RS_applied = RS_target_;
        freeze_ou_channel_ = false;
        freeze_RS_channel_ = false;
        apply_ou_tune_(true);
        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_RS_tune_();
        }
        return true;
    }

    // Freeze one adaptation channel while the other keeps tracking the sea.
    //
    // The deployed law couples the two: r_S = clip(c * sigma_a * tau^3), so
    // simply freezing tau and sigma_a with setFixedTuning() freezes r_S as
    // well, and an ablation built that way cannot say which channel carries
    // the benefit.  Here the tuner keeps running and keeps deriving r_S from
    // its *live* tau and sigma_a estimates; only the channels named below are
    // held at the supplied operating point on the way to the filter.
    //
    // freeze_ou   holds the OU process parameters (tau, sigma_a) at
    //             (tau_s, sigma_a) while r_S continues to adapt.
    // freeze_RS   holds the integral pseudo-measurement scale at RS while
    //             tau and sigma_a continue to adapt.
    //
    // Freezing both is equivalent to setFixedTuning() and is rejected here so
    // that the two entry points do not silently overlap.
    bool setChannelFreeze(bool freeze_ou,
                          float tau_s,
                          float sigma_a,
                          bool freeze_RS,
                          float RS)
    {
        if (freeze_ou == freeze_RS) return false;
        if (freeze_ou && !(std::isfinite(tau_s) && tau_s > 0.0f &&
                           std::isfinite(sigma_a) && sigma_a > 0.0f))
        {
            return false;
        }
        if (freeze_RS && !(std::isfinite(RS) && RS > 0.0f)) return false;

        enable_tuner_ = true;
        freeze_ou_channel_ = freeze_ou;
        freeze_RS_channel_ = freeze_RS;

        if (freeze_ou) {
            frozen_tau_s_ = enable_clamp_
                ? std::min(std::max(tau_s, min_tau_s_), max_tau_s_)
                : tau_s;
            frozen_sigma_a_ = enable_clamp_
                ? std::min(sigma_a, max_sigma_a_)
                : sigma_a;
            tau_target_ = frozen_tau_s_;
            sigma_target_ = frozen_sigma_a_;
            tune_.tau_applied = tau_target_;
            tune_.sigma_applied = sigma_target_;
            apply_ou_tune_(true);
        }
        if (freeze_RS) {
            frozen_RS_ = enable_clamp_
                ? std::min(std::max(RS, min_R_S_), max_R_S_)
                : RS;
            RS_target_ = frozen_RS_;
            tune_.RS_applied = RS_target_;
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
                apply_RS_tune_();
            }
        }
        return true;
    }

    bool frozenOUChannel() const noexcept { return freeze_ou_channel_; }
    bool frozenRSChannel() const noexcept { return freeze_RS_channel_; }

    // Enable/disable use of the extended linear block [v,p,S,a_w] in Kalman3D_Wave_OU_III.
    void enableLinearBlock(bool flag = true) {
        enable_linear_block_ = flag;
        if (mekf_) {
            const bool on_now = flag && (startup_stage_ == StartupStage::Live);
            mekf_->set_linear_block_enabled(on_now);
        }
    }

    void setFreqBounds(float min_hz, float max_hz) {
        if (!std::isfinite(min_hz) || !std::isfinite(max_hz)) return;
        if (min_hz <= 0.0f || max_hz <= min_hz) return;
        min_freq_hz_ = min_hz;
        max_freq_hz_ = max_hz;
        freq_stillness_.setTargetFreqHz(min_freq_hz_);
    }

    void setTauBounds(float min_tau_s, float max_tau_s) {
        if (!std::isfinite(min_tau_s) || !std::isfinite(max_tau_s)) return;
        if (min_tau_s <= 0.0f || max_tau_s <= min_tau_s) return;
        min_tau_s_ = min_tau_s;
        max_tau_s_ = max_tau_s;
    }

    void setMaxSigmaA(float max_sigma_a) {
        if (!std::isfinite(max_sigma_a) || max_sigma_a <= 0.0f) return;
        max_sigma_a_ = max_sigma_a;
    }

    void setRSBounds(float min_RS, float max_RS) {
        if (!std::isfinite(min_RS) || !std::isfinite(max_RS)) return;
        if (min_RS <= 0.0f || max_RS <= min_RS) return;
        min_R_S_ = min_RS;
        max_R_S_ = max_RS;
    }

    void setAdaptationTimeConstants(float tau_sec) {
        if (std::isfinite(tau_sec) && tau_sec > 0.0f)   adapt_tau_sec_   = tau_sec;
     }

    // Smoothing-horizon multiplier for the r_S channel.  The EMA time
    // constant is mult * tau_target, so the horizon follows the sea state
    // instead of being pinned to one second count.  r_S ~ sigma_aw * tau^3
    // amplifies a tau error by the third power, which is what sets the
    // multiplier apart from the tau/sigma one; see
    // docs/ou-ema-adaptation-tuning.md.
    void setRSAdaptMult(float m) {
        if (std::isfinite(m) && m > 0.0f) adapt_RS_mult_ = m;
    }

    // Size, in natural-log units of the r_S target-to-applied ratio, of a
    // discrepancy the smoother should treat as a real sea-state move rather
    // than tuner jitter.  Zero or negative leaves the plain proportional
    // horizon.  See seastate::common::adaptiveSmoothingHorizonSec.
    void setRSAdaptSlewLog(float d) {
        if (std::isfinite(d)) adapt_RS_slew_log_ = d;
    }

    float getRSAdaptMult() const noexcept { return adapt_RS_mult_; }
    float getRSAdaptSlewLog() const noexcept { return adapt_RS_slew_log_; }

    void setAdaptationUpdatePeriod(float every_sec) {
        if (std::isfinite(every_sec) && every_sec > 0.0f) {
            adapt_every_secs_ = every_sec;
        }
    }

    void setOnlineTuneWarmupSec(float warmup_sec) {
        if (std::isfinite(warmup_sec) && warmup_sec >= 0.0f) {
            online_tune_warmup_sec_ = warmup_sec;
        }
    }

    void setMagDelaySec(float delay_sec) {
        if (std::isfinite(delay_sec) && delay_sec >= 0.0f) {
            mag_delay_sec_ = delay_sec;
        }
    }

    void setFreezeAccBiasUntilLive(bool en) { freeze_acc_bias_until_live_ = en; }
    void setWarmupRaccStd(float r) { if (std::isfinite(r) && r > 0.0f) Racc_warmup_std_ = r; }

    // For SeaStateFusionFilter_OU_III to restore Racc automatically
    void setNominalRaccStd(const Eigen::Vector3f& r) { Racc_nominal_ = r; }

    //  Exposed getters
    inline float getFreqHz()        const noexcept { return freq_hz_; }        // fast branch
    inline float getFreqSlowHz()    const noexcept { return freq_hz_slow_; }   // slow branch
    inline float getFreqRawHz()     const noexcept { return f_raw; }
    inline float getTauApplied()    const noexcept { return tune_.tau_applied; }
    inline float getSigmaApplied()  const noexcept { return tune_.sigma_applied; }
    inline float getRSApplied()     const noexcept { return tune_.RS_applied; }
    inline float getTauTarget()     const noexcept { return tau_target_;   }
    inline float getSigmaTarget()   const noexcept { return sigma_target_; }
    inline float getRSTarget()      const noexcept { return RS_target_;    }

    // Use slow frequency as a more stable "period" proxy
    inline float getPeriodSec() const noexcept {
        return (freq_hz_slow_ > 1e-6f) ? 1.0f / freq_hz_slow_ : NAN;
    }

    inline float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }

    // Returns the BODY-Z-based up-positive proxy used by tracker/tuner logic.
    // This is not a true vertical acceleration estimate.
    inline float getAccelVertical() const noexcept { return a_body_z_up_proxy_; }

    inline float getHeaveAbs() const noexcept { if (!mekf_) return NAN; return std::fabs(mekf_->get_position().z()); }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau = smoothed ? tune_.tau_applied : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
        constexpr float C_HS  = 2.0f * std::sqrt(2.0f) / (std::numbers::pi_v<float> * std::numbers::pi_v<float>);
        return C_HS * sigma * tau * tau / 2.0f;
    }

    float getVerticalSpeedEnvelopeMps(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!(tau > 1e-6f) || !std::isfinite(tau) || !std::isfinite(sigma)) return NAN;
        constexpr float K = std::sqrt(2.0f) / std::numbers::pi_v<float>;
        const float v_env = K * sigma * tau;
        return std::isfinite(v_env) ? v_env : NAN;
    }

    // Zero-crossing wave period [s] from the independent accelerometer-only
    // estimator; NaN until it settles.
    inline float getWavePeriodSec() const noexcept { return wave_period_.getPeriodSec(); }
    inline bool wavePeriodReady() const noexcept { return wave_period_.isReady(); }

    // Drive the operating point from the wave band (default) or from the
    // acceleration-band tracker, which is what the filter did before and is
    // kept so the change can be ablated rather than assumed.
    void setWaveBandTuning(bool flag) { wave_band_tuning_ = flag; }
    bool waveBandTuning() const noexcept { return wave_band_tuning_; }

    // Select which vertical acceleration drives the wave-period estimator.
    // Complementary (default) levels with the private Mahony observer and is
    // measurement-only, so the tuner is outside the estimator's loop.  Leveled
    // is the older behaviour, which levels with the main filter's attitude and
    // closes that loop; BodyZ is measurement-only but unlevelled.  See the call
    // site in updateTime for what each one costs.
    void setWavePeriodInput(WavePeriodInputSource source) {
        wave_period_input_ = source;
    }
    WavePeriodInputSource wavePeriodInput() const noexcept {
        return wave_period_input_;
    }

    // Select which vertical acceleration the frequency tracker runs on.
    // BodyZ (default) is the raw proxy; Complementary is the levelled signal
    // from the private Mahony observer.  Both are measurement-only.
    void setFreqTrackerInput(FreqTrackerInputSource source) {
        freq_tracker_input_ = source;
    }
    FreqTrackerInputSource freqTrackerInput() const noexcept {
        return freq_tracker_input_;
    }

    // Gains of the private Mahony observer that levels the default input.
    // two_kp sets the accelerometer-to-gyro correction corner, which must stay
    // below the wave band; see VerticalAccelComplementary.h.
    void setWavePeriodComplementaryGains(float two_kp, float two_ki) {
        vertical_accel_comp_.setGains(two_kp, two_ki);
    }

    inline WaveDirection getDirSignState() const noexcept { return dir_sign_state_; }

    // Propagation-plane angle relative to boat +X, modulo 180 degrees.
    inline float getWaveAxisDeg() const noexcept { return dir_filter_.getAxisDegrees(); }
    inline float getWaveDirectionDeg() const noexcept { return getWaveAxisDeg(); }

    // Fully directed apparent propagation angles observed by the moving boat.
    // These are encounter/apparent directions unless vessel-motion correction
    // is applied externally (see wave_dir/WaveEncounter.h).
    inline float getApparentWaveDirectionToDeg() const noexcept {
        return dir_sign_.getDirectedAngleDegrees();
    }
    inline float getApparentWaveDirectionFromDeg() const noexcept {
        return dir_sign_.getWaveFromAngleDegrees();
    }
    inline float getDirSenseCoherence() const noexcept {
        return dir_sign_.getCoherence();
    }

    inline auto& mekf() noexcept { return *mekf_; }
    inline const auto& mekf() const noexcept { return *mekf_; }

    inline KalmanWaveDirection& dir() noexcept { return dir_filter_; }
    inline const KalmanWaveDirection& dir() const noexcept { return dir_filter_; }

    inline WaveDirectionDetector<float>& dir_sign() noexcept { return dir_sign_; }
    inline const WaveDirectionDetector<float>& dir_sign() const noexcept { return dir_sign_; }

private:

    // Simple first-order low-pass filter for vertical accel → tracker input
    using FreqInputLPF = seastate::common::FreqInputLPF;
    using StillnessAdapter = seastate::common::StillnessAdapter;

    // sync_covariance is set only by discrete reconfiguration events. The
    // periodic adaptation path leaves the posterior a_w marginal alone; the
    // new stationary scale reaches the filter through the OU process
    // covariance instead.
    void apply_ou_tune_(bool sync_covariance) {
        if (!mekf_) return;
        mekf_->set_aw_time_constant(tune_.tau_applied);

        const float sigma_floor = std::max(0.05f, acc_noise_floor_sigma_);
        const float sZ = std::max(sigma_floor, tune_.sigma_applied);
        const float sH = sZ * S_factor_;
        const Eigen::Vector3f aw_std(sH, sH, sZ);
        mekf_->set_aw_stationary_std(aw_std);
        if (sync_covariance) {
            apply_aw_cov_sync_();
            last_aw_cov_sync_sec_ = time_;
        }
    }

    // Re-align the posterior a_w marginal with the stationary prior at the
    // adaptation cadence. Runs independently of the tuner so that
    // fixed-tuning modes apply the same policy and remain matched controls.
    void periodic_aw_cov_sync_tick_() {
        if (!periodic_aw_cov_sync_ || !mekf_) return;
        if (startup_stage_ != StartupStage::Live) return;
        if (time_ - last_aw_cov_sync_sec_ <= adapt_every_secs_) return;
        apply_aw_cov_sync_();
        last_aw_cov_sync_sec_ = time_;
    }

    void apply_aw_cov_sync_() {
        if (congruent_aw_cov_sync_) {
            mekf_->synchronize_aw_covariance_to_stationary_congruent();
        } else {
            mekf_->synchronize_aw_covariance_to_stationary();
        }
    }

    void apply_RS_tune_(float rs_scale = 1.0f) {
        if (!mekf_) return;
        const float s = (std::isfinite(rs_scale) && rs_scale > 0.0f)
                        ? std::min(rs_scale, 1.0f)
                        : 1.0f;
        const float RSb = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);
        const float rs_xy = RSb * s * R_S_xy_factor_;
        mekf_->set_RS_noise(Eigen::Vector3f(
            rs_xy,
            rs_xy,
            RSb * s
        ));
    }

    void update_tuner(float dt, float a_body_z_up_proxy, float freq_hz_for_tuner) {
        tuner_.update(dt, a_body_z_up_proxy, freq_hz_for_tuner);

        // Startup stage logic
        switch (startup_stage_) {
           case StartupStage::Cold:
               if (startup_stage_t_ >= online_tune_warmup_sec_) {
                   startup_stage_   = StartupStage::TunerWarm;
                   startup_stage_t_ = 0.0f;
               }
               return;

          case StartupStage::TunerWarm:
              if (!tuner_.isFreqReady()) return;
              if (tuner_.isReady()) {
                  enterLive_();
              }
              break;

           case StartupStage::Live:
               break;
        }

        // The tuning frequency is a wave-band quantity and has to be allowed
        // below the tracker's floor: a developed sea has T_z = 8.6 s, i.e.
        // 0.12 Hz, well under the 0.2 Hz the tracker is bounded to.
        const float f_tune_floor = wave_band_tuning_ ? min_tune_freq_hz_ : min_freq_hz_;
        float f_tune = tuner_.getFrequencyHz();
        if (!std::isfinite(f_tune) || f_tune < f_tune_floor) {
            f_tune = f_tune_floor;
        }
        if (f_tune > max_freq_hz_) {
            f_tune = max_freq_hz_;
        }

        float var_total = acc_noise_floor_sigma_ * acc_noise_floor_sigma_;
        if (tuner_.isVarReady()) {
            var_total = std::max(0.0f, tuner_.getAccelVariance());
        }
        const float var_noise = acc_noise_floor_sigma_ * acc_noise_floor_sigma_;
        float var_wave = var_total - var_noise;
        if (var_wave < 0.0f) var_wave = 0.0f;

        if (freq_stillness_.isStill()) {
            const float still_t = std::max(0.0f, freq_stillness_.getStillTime());
            constexpr float STILL_VAR_DECAY_SEC = 1.0f;
            float atten = std::exp(-still_t / STILL_VAR_DECAY_SEC);
            atten = std::min(std::max(atten, 0.0f), 1.0f);
            var_wave *= atten;
        }

        var_wave = std::max(var_wave, 1e-6f);
        float sigma_wave = std::sqrt(var_wave);
        float tau_raw = tau_coeff_ * 0.5f / f_tune;

        if (enable_clamp_) {
            tau_target_   = std::min(std::max(tau_raw,  min_tau_s_), max_tau_s_);
            sigma_target_ = std::min(sigma_wave * sigma_coeff_,      max_sigma_a_);
        } else {
            tau_target_   = tau_raw;
            sigma_target_ = sigma_wave;
        }
        if (!tuner_.isVarReady()) {
            sigma_target_ = std::max(sigma_target_, std::max(0.05f, acc_noise_floor_sigma_));
        }

        // r_S is derived from the *live* tau and sigma_a estimates, before any
        // channel freeze is applied below.  Deriving it from frozen values
        // instead would make "freeze the OU channel" silently freeze r_S too,
        // which is precisely the confound the channel ablation exists to
        // remove.
        float RS_raw = R_S_coeff_ * sigma_target_
                       * tau_target_ * tau_target_ * tau_target_;

        if (enable_clamp_) {
            RS_target_ = std::min(std::max(RS_raw, min_R_S_), max_R_S_);
        } else {
            RS_target_ = RS_raw;
        }

        if (freeze_ou_channel_) {
            tau_target_   = frozen_tau_s_;
            sigma_target_ = frozen_sigma_a_;
        }
        if (freeze_RS_channel_) {
            RS_target_ = frozen_RS_;
        }
        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float RS_t) {
        const float alpha = 1.0f - std::exp(-dt / adapt_tau_sec_);

        const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(
            adapt_RS_mult_, tau_t, RS_t, tune_.RS_applied, adapt_RS_slew_log_, dt);
        const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);

        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);
        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);

        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {
            if (tuner_.isFreqReady()) {
                apply_ou_tune_(false);
            }
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
                apply_RS_tune_();
            }
            last_adapt_time_sec_ = time_;
        }
    }

    // Vertical acceleration the wave-period estimator is driven by.  Each
    // measurement-only source falls back to the body-Z proxy while it is
    // unusable, which is what the leveled source already does when heading is
    // not yet resolved.
    float wave_period_input_ms2_(
        const wave_direction::HeadingFrameAcceleration<float>& leveled) const
    {
        switch (wave_period_input_) {
            case WavePeriodInputSource::BodyZ:
                return a_body_z_up_proxy_;
            case WavePeriodInputSource::Complementary:
                return vertical_accel_comp_.isReady()
                           ? vertical_accel_comp_.verticalAccelUpMs2()
                           : a_body_z_up_proxy_;
            case WavePeriodInputSource::Leveled:
            default:
                return leveled.heading_valid ? leveled.up_ms2
                                             : a_body_z_up_proxy_;
        }
    }

    float tuner_frequency_hz_(float tracker_hz) const {
        if (wave_band_tuning_ && wave_period_.isReady()) {
            const float wave_hz = wave_period_.getFrequencyHz();
            if (std::isfinite(wave_hz) && wave_hz > 0.0f) return wave_hz;
        }
        return tracker_hz;
    }

    void resetTrackingState_() {
        tracker_policy_       = TrackingPolicy{};
        wave_period_          = WavePeriodEstimator{};
        vertical_accel_comp_.reset();
        freq_input_lpf_       = FreqInputLPF{};
        freq_stillness_       = StillnessAdapter(g_std, min_freq_hz_, FREQ_GUESS);
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);

        tuner_.reset();

        freq_fast_smoother_   = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 3.5f);
        freq_slow_smoother_   = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 10.0f);

        freq_hz_      = FREQ_GUESS;
        freq_hz_slow_ = FREQ_GUESS;
        f_raw         = FREQ_GUESS;

        dir_filter_ = KalmanWaveDirection(2.0f * static_cast<float>(M_PI) * FREQ_GUESS);
        dir_sign_.reset();
        dir_sign_state_ = UNCERTAIN;

        last_adapt_time_sec_ = time_;
        last_aw_cov_sync_sec_ = time_;
    }

    void enterCold_() {
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;
        mekf_->set_linear_block_enabled(false);

        accel_bias_locked_   = with_mag_;
        mag_updates_applied_ = 0;
        first_mag_update_time_  = NAN;

        if (freeze_acc_bias_until_live_) {
            mekf_->set_acc_bias_updates_enabled(false);
            mekf_->set_Racc_std(Eigen::Vector3f::Constant(Racc_warmup_std_));
            warmup_Racc_active_ = true;
        }
    }

    void enterLive_() {
        startup_stage_   = StartupStage::Live;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;
        apply_ou_tune_(true);
        mekf_->set_linear_block_enabled(enable_linear_block_);

        if (freeze_acc_bias_until_live_) {
            const bool allow_bias = !accel_bias_locked_;
            mekf_->set_acc_bias_updates_enabled(allow_bias);

            if (warmup_Racc_active_ &&
                Racc_nominal_.allFinite() &&
                Racc_nominal_.maxCoeff() > 0.0f)
            {
                mekf_->set_Racc_std(Racc_nominal_);
            }
            warmup_Racc_active_ = false;
        }

        if (enable_linear_block_) apply_RS_tune_();
    }

    StartupStage startup_stage_    = StartupStage::Cold;
    float        startup_stage_t_  = 0.0f;

    // Warmup behavior
    bool  freeze_acc_bias_until_live_ = true;
    float Racc_warmup_std_            = 0.6f;
    bool  warmup_Racc_active_         = false;
    Eigen::Vector3f Racc_nominal_     = Eigen::Vector3f::Constant(0.0f);

    bool accel_bias_locked_ = true;
    int  mag_updates_applied_ = 0;
    static constexpr int MAG_UPDATES_TO_UNLOCK = 250;

    bool   with_mag_;
    double time_;
    double last_adapt_time_sec_;

    float first_mag_update_time_ = NAN;

    float tilt_over_limit_sec_ = 0.0f;
    float tilt_reset_cooldown_sec_ = 0.0f;

    float freq_hz_       = FREQ_GUESS;
    float freq_hz_slow_  = FREQ_GUESS;
    float f_raw          = FREQ_GUESS;

    float a_body_z_up_proxy_ = 0.0f;

    bool enable_clamp_ = true;
    bool enable_tuner_ = true;

    // Per-channel freezes for the partial-adaptation ablation; see
    // setChannelFreeze().  Both false is the deployed fully adaptive filter.
    bool freeze_ou_channel_ = false;
    bool freeze_RS_channel_ = false;
    float frozen_tau_s_ = NAN;
    float frozen_sigma_a_ = NAN;
    float frozen_RS_ = NAN;

    bool enable_linear_block_ = true;

    // Covariance-inflation policy; see setPeriodicAwCovarianceSync.
    bool   periodic_aw_cov_sync_ = true;
    bool   congruent_aw_cov_sync_ = false;
    double last_aw_cov_sync_sec_ = 0.0;

    bool  wave_band_tuning_       = true;
    WavePeriodInputSource wave_period_input_ = WavePeriodInputSource::Complementary;
    FreqTrackerInputSource freq_tracker_input_ = FreqTrackerInputSource::BodyZ;
    float min_tune_freq_hz_       = MIN_TUNE_FREQ_HZ;
    float min_freq_hz_            = MIN_FREQ_HZ;
    float max_freq_hz_            = MAX_FREQ_HZ;
    float min_tau_s_              = MIN_TAU_S;
    float max_tau_s_              = MAX_TAU_S;
    float max_sigma_a_            = MAX_SIGMA_A;
    float min_R_S_                = MIN_R_S;
    float max_R_S_                = MAX_R_S;
    float adapt_tau_sec_          = ADAPT_TAU_SEC;
    float adapt_RS_mult_          = ADAPT_RS_MULT;
    float adapt_RS_slew_log_      = ADAPT_RS_SLEW_LOG;
    float adapt_every_secs_       = ADAPT_EVERY_SECS;
    float online_tune_warmup_sec_ = ONLINE_TUNE_WARMUP_SEC;
    float mag_delay_sec_          = MAG_DELAY_SEC;

    // Horizontal integral-regularization scale relative to the vertical one.
    // 0.36 made the horizontal high-pass 2.8x stronger than the vertical one,
    // which was a small-sea optimum applied to every sea state: with the
    // operating point now tied to the wave band, every stationary record scores
    // better with the two equal, by 7 to 27 percent of 3D RMS in the two
    // largest seas.  Retained as a setter because the bound is a real one.
    float R_S_xy_factor_ = 1.0f;
    float S_factor_      = 1.87f;

    TrackingPolicy                  tracker_policy_{};
    FirstOrderIIRSmoother<float>    freq_fast_smoother_{FREQ_SMOOTHER_DT, 3.5f};
    FirstOrderIIRSmoother<float>    freq_slow_smoother_{FREQ_SMOOTHER_DT, 10.0f};
    SeaStateAutoTuner               tuner_;
    WavePeriodEstimator             wave_period_;
    VerticalAccelComplementary      vertical_accel_comp_{};
    TuneState                       tune_;

    float tau_target_   = NAN;
    float sigma_target_ = NAN;
    float RS_target_    = NAN;

    float acc_noise_floor_sigma_ = ACC_NOISE_FLOOR_SIGMA_DEFAULT;

    // r_S = R_S_coeff * sigma_aw * tau^3 and tau = tau_coeff * T_z / 2.  Both
    // coefficients are re-fitted for the wave-band period: tau_coeff = 1 is the
    // documented intent, tau equal to half the zero-crossing period, and
    // R_S_coeff was fitted on the four stationary JONSWAP records against the
    // per-record optimum located by a fixed-r_S scan.
    float R_S_coeff_    = 0.35f;
    float tau_coeff_    = 1.0f;
    float sigma_coeff_  = 0.9f;

    std::unique_ptr<Kalman3D_Wave_OU_III<float>>  mekf_;
    KalmanWaveDirection                    dir_filter_{2.0f * static_cast<float>(M_PI) * FREQ_GUESS};

    FreqInputLPF        freq_input_lpf_;
    StillnessAdapter    freq_stillness_;

    WaveDirectionDetector<float> dir_sign_{0.002f, 0.005f};
    WaveDirection                dir_sign_state_ = UNCERTAIN;
};

template<TrackerType trackerT>
class SeaStateFusion_OU_III {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Config {
        bool with_mag = true;

        bool enable_linear_block = true;
        bool require_mag_lock_for_linear_block = false;

        float mag_delay_sec          = MAG_DELAY_SEC;
        float online_tune_warmup_sec = 10.0f;

        bool  freeze_acc_bias_until_live = true;
        float Racc_warmup_std = 0.5f;

        Eigen::Vector3f sigma_a = Eigen::Vector3f(0.2f, 0.2f, 0.2f);
        Eigen::Vector3f sigma_g = Eigen::Vector3f(0.01f, 0.01f, 0.01f);
        Eigen::Vector3f sigma_m = Eigen::Vector3f(0.3f, 0.3f, 0.3f);

        // Mag-start gate: gravity-direction agreement using current tilt.
        float mag_gravity_align_max_sin   = 0.075f; // sin(deg)
        float mag_gravity_align_hold_sec  = 2.0f;
        float mag_gravity_align_lpf_tau   = 1.0f;
        float mag_tilt_fallback_sec       = 30.0f;
        float mag_extreme_gyro_dps        = 30.0f; // veto only truly violent motion
        float mag_init_min_mag_norm       = 1e-3f;

        // Real-device mag acquisition.
        //
        // The important rule:
        //
        //   MagAutoTuner must never receive the MEKF's yaw, which is arbitrary
        //   and unobservable before mag lock.  It receives the MEKF quaternion
        //   with yaw divided out instead: that tilt is invariant under
        //   q_bw -> Rz(psi) q_bw, so no heading can leak through it, and it is
        //   a far better level frame in waves than one rebuilt from accel.
        //
        // The reference is an average of the field in that tilt frame, so
        // whatever tilt error survives the window survives in the reference.
        // In waves the error is periodic, so what the window has to buy is
        // whole wave periods, not samples: 128 samples is 5.1 s at a 25 Hz mag
        // ODR, short enough to lock in the phase it started on rather than
        // cancel it.  15 s covers a couple of periods across the band these
        // filters work in and captures most of what a much longer window would,
        // at a startup cost of 15 s rather than 40 s.  Held in seconds so it
        // does not silently shorten at a higher ODR.
        int   mag_min_samples              = 128;
        float mag_min_window_sec           = 15.0f;
        float mag_max_window_sec           = 0.0f;
        float mag_sample_dt_sec            = 1.0f / 200.0f;

        // Keep off in waves. Accel/gyro weighting can phase-select wave motion.
        // Body-frame hard-iron offset, learned during startup alongside the
        // reference and then subtracted from every magnetometer sample.
        //
        // Off by default.  The MEKF has no mag-bias state, so an offset left in
        // the stream is heading error one-for-one against the horizontal field;
        // but the offset is only weakly separable from the reference at a fixed
        // heading, and a wrong one subtracted everywhere is worse than none.
        // Turn it on where the platform changes heading during startup.
        bool  mag_estimate_hard_iron = false;
        bool  mag_enable_quality_weighting = false;
        float mag_min_effective_weight     = 0.0f;
        float mag_acc_norm_rel_soft        = 0.22f;
        float mag_gyro_soft_dps            = 45.0f;

        // Bootstrap tilt observer for dynamic motion in waves.
        float bootstrap_tilt_obs_acc_tau_sec  = 2.15f; // accel correction time constant
        float bootstrap_gravity_slow_tau_sec  = 6.0f; // slow gravity reference LPF
        float bootstrap_gravity_align_max_sin = 0.070f; // sin(deg)
        float bootstrap_gravity_hold_sec      = 2.0f;
        float bootstrap_gravity_min_sec       = 6.87f;
        float bootstrap_gravity_timeout_sec   = 15.0f;
        float bootstrap_gravity_norm_frac     = 0.22f; // downweight accel when |a| departs from g

        bool enable_displacement_detrend = false;
        bool use_custom_displacement_detrend_cfg = false;
        AdaptiveWaveDetrender3D::Config displacement_detrend_cfg{};
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;

        begun_ = true;
        stage_ = Stage::Uninitialized;
        t_ = 0.0f;

        gravity_gate_acc_lpf_.reset();
        mag_gravity_good_sec_ = 0.0f;
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;

        mag_ref_set_ = false;

        last_mag_tilt_frame_yaw_rad_ = NAN;
        last_mag_startup_yaw_correction_rad_ = NAN;

        MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min = cfg_.mag_init_min_mag_norm;
        mag_cfg.min_samples  = cfg_.mag_min_samples;

        mag_cfg.min_window_sec = cfg_.mag_min_window_sec;
        mag_cfg.max_window_sec = cfg_.mag_max_window_sec;
        mag_cfg.sample_dt_sec  = cfg_.mag_sample_dt_sec;

        mag_cfg.gravity_ref = g_std;
        mag_cfg.enable_quality_weighting = cfg_.mag_enable_quality_weighting;
        mag_cfg.estimate_hard_iron       = cfg_.mag_estimate_hard_iron;
        mag_cfg.min_effective_weight     = cfg_.mag_min_effective_weight;
        mag_cfg.acc_norm_rel_soft        = cfg_.mag_acc_norm_rel_soft;
        mag_cfg.gyro_soft_dps            = cfg_.mag_gyro_soft_dps;

        mag_auto_tuner_.setConfig(mag_cfg);

        resetTiltInit_();

        last_acc_body_ned_.setZero();
        last_gyro_body_ned_.setZero();
        have_last_imu_ = false;

        impl_.setWithMag(cfg_.with_mag);
        impl_.setFreezeAccBiasUntilLive(cfg_.freeze_acc_bias_until_live);
        impl_.setWarmupRaccStd(cfg_.Racc_warmup_std);
        impl_.setMagDelaySec(0.0f); // outer wrapper owns mag delay
        impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);

        impl_.initialize(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m);
        last_impl_startup_stage_ = impl_.getStartupStage();

        syncLinearBlockGate_();

        impl_.setNominalRaccStd(cfg_.sigma_a);

        displacement_up_m_.setZero();
        displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

        if (cfg_.enable_displacement_detrend) {
            if (cfg_.use_custom_displacement_detrend_cfg) {
                displacement_detrender_.setConfig(cfg_.displacement_detrend_cfg);
            } else {
                displacement_detrender_.setConfig(
                    seastate::common::defaultDisplacementDetrenderConfig<AdaptiveWaveDetrender3D::Config>(FREQ_GUESS));
            }

            displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
        }
    }

    void update(float dt,
                const Eigen::Vector3f& gyro_body_ned,
                const Eigen::Vector3f& acc_body_ned,
                float tempC = 35.0f)
    {
        if (!begun_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;

        t_ += dt;

        if (stage_ == Stage::Uninitialized) {
            const bool tilt_ready = seastate::common::runStartupGravityInit(
                gyro_body_ned,
                acc_body_ned,
                dt,
                t_,
                g_std,
                cfg_.bootstrap_tilt_obs_acc_tau_sec,
                cfg_.bootstrap_gravity_slow_tau_sec,
                cfg_.bootstrap_gravity_align_max_sin,
                cfg_.bootstrap_gravity_hold_sec,
                cfg_.bootstrap_gravity_min_sec,
                cfg_.bootstrap_gravity_timeout_sec,
                cfg_.bootstrap_gravity_norm_frac,
                bootstrap_tilt_obs_,
                bootstrap_gravity_slow_lpf_,
                bootstrap_gravity_good_sec_,
                [this](const Eigen::Vector3f& acc_init) {
                    impl_.initialize_from_acc(acc_init);
                });

            if (tilt_ready) {
                stage_ = Stage::Warming;
            }
        }

        last_acc_body_ned_  = acc_body_ned;
        last_gyro_body_ned_ = gyro_body_ned;
        have_last_imu_      = true;

        if (stage_ != Stage::Uninitialized) {
            impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);

            const Eigen::Vector3f acc_gate_lp =
                gravity_gate_acc_lpf_.step(
                    acc_body_ned,
                    dt,
                    cfg_.mag_gravity_align_lpf_tau);

            const float align_sin =
                seastate::common::gravityAlignResidualSin(
                    impl_.mekf().quaternion_boat(),
                    acc_gate_lp);

            const float gyro_dps =
                gyro_body_ned.norm() * 57.295779513f;

            const bool extreme_motion =
                !std::isfinite(gyro_dps) ||
                (gyro_dps > cfg_.mag_extreme_gyro_dps);

            const bool gravity_good_now =
                std::isfinite(align_sin) &&
                (align_sin <= cfg_.mag_gravity_align_max_sin) &&
                !extreme_motion;

            if (gravity_good_now) {
                mag_gravity_good_sec_ += dt;
                if (mag_gravity_good_sec_ > 10.0f) {
                    mag_gravity_good_sec_ = 10.0f;
                }
            } else {
                mag_gravity_good_sec_ =
                    std::max(0.0f, mag_gravity_good_sec_ - 2.0f * dt);
            }

            const Eigen::Vector3f pos_ned_m = impl_.mekf().get_position();

            displacement_up_m_ =
                Eigen::Vector3f(
                    pos_ned_m.x(),
                    pos_ned_m.y(),
                    -pos_ned_m.z());

            if (cfg_.enable_displacement_detrend) {
                const float wave_hz = impl_.getFreqHz();

                const bool ext_freq_valid =
                    isLive() &&
                    std::isfinite(wave_hz) &&
                    (wave_hz >= displacement_detrender_.config().min_wave_freq_hz) &&
                    (wave_hz <= displacement_detrender_.config().max_wave_freq_hz);

                displacement_det_out_ =
                    displacement_detrender_.update(
                        displacement_up_m_,
                        dt,
                        wave_hz,
                        ext_freq_valid);
            } else {
                displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
                displacement_det_out_.input = displacement_up_m_;
                displacement_det_out_.baseline_slow = Eigen::Vector3f::Zero();
                displacement_det_out_.wave_raw = displacement_up_m_;
                displacement_det_out_.wave_clean = displacement_up_m_;
            }
        }

        const auto cur_stage = impl_.getStartupStage();

        if (cur_stage != last_impl_startup_stage_) {
            if (cur_stage == SeaStateFusionFilter_OU_III<trackerT>::StartupStage::Cold) {
                mag_ref_set_ = false;
                mag_auto_tuner_.reset();

                gravity_gate_acc_lpf_.reset();
                mag_gravity_good_sec_ = 0.0f;
                mag_init_eligible_t0_ = NAN;
                last_mag_sample_t_ = NAN;

                last_mag_tilt_frame_yaw_rad_ = NAN;
                last_mag_startup_yaw_correction_rad_ = NAN;
                
                syncLinearBlockGate_();

                if (stage_ != Stage::Live) {
                    stage_ = Stage::Warming;

                    displacement_up_m_.setZero();
                    displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

                    if (cfg_.enable_displacement_detrend) {
                        displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
                    }
                }
            }

            last_impl_startup_stage_ = cur_stage;
        }

        if (stage_ == Stage::Warming && impl_.isAdaptiveLive()) {
            stage_ = Stage::Live;
        }

        // Re-apply gate every update.
        // Inner impl only enables the actual MEKF linear block when its own stage is Live.
        syncLinearBlockGate_();
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!begun_ || !cfg_.with_mag) return;
        if (stage_ == Stage::Uninitialized) return;
        if (t_ < cfg_.mag_delay_sec) return;

        if (!std::isfinite(mag_init_eligible_t0_)) {
            mag_init_eligible_t0_ = t_;
        }

        const bool gravity_trusted =
            (mag_gravity_good_sec_ >= cfg_.mag_gravity_align_hold_sec);

        const bool fallback_ok =
            ((t_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec);

        if (!mag_ref_set_) {
            if (!gravity_trusted && !fallback_ok) {
                return;
            }

            if (have_last_imu_) {
                const float dt_mag =
                    (std::isfinite(last_mag_sample_t_) && t_ > last_mag_sample_t_)
                        ? (t_ - last_mag_sample_t_)
                        : cfg_.mag_sample_dt_sec;

                last_mag_sample_t_ = t_;

                // Accumulate in the MEKF's own tilt frame with yaw removed.
                //
                // Stripping yaw makes the frame invariant to the MEKF's
                // arbitrary startup heading, so this leaks no yaw into the
                // learned reference -- q_bw and Rz(psi) q_bw give the same
                // tilt.  A gravity-only frame rebuilt from low-passed accel
                // would be yaw-free too, but in waves that accel is gravity
                // plus a phase-lagged remnant of the orbital specific force,
                // so its tilt is wrong by a wave-correlated angle that the
                // averaging window is too short to cancel.
                const Eigen::Quaternionf q_tilt_bw =
                    tiltOnlyQuatFromBoatQuat_(impl_.mekf().quaternion_boat());

                if (mag_auto_tuner_.addSampleWithTiltQuatDt(
                        dt_mag,
                        q_tilt_bw,
                        last_acc_body_ned_,
                        last_gyro_body_ned_,
                        mag_body_ned))
                {
                    Eigen::Vector3f mag_world_ref_uT;

                    if (mag_auto_tuner_.getMagWorldRef(mag_world_ref_uT) &&
                        mag_world_ref_uT.allFinite() &&
                        mag_world_ref_uT.norm() > cfg_.mag_init_min_mag_norm)
                    {
                        // This reference was learned in a yaw-stripped tilt
                        // frame, so it carries no MEKF heading.
                        impl_.mekf().set_mag_world_ref(mag_world_ref_uT);

                        const float mag_tilt_yaw_rad =
                            mag_auto_tuner_.getYawGaugeCorrectionRad();

                        if (std::isfinite(mag_tilt_yaw_rad)) {
                            // One-time yaw-gauge lock.
                            //
                            // mag_tilt_yaw_rad is the heading of averaged
                            // magnetic north in the accumulation frame, so
                            // driving the boat's absolute yaw to its negative
                            // puts the learned north on +X.
                            //
                            // Only yaw is written.  The MEKF's tilt is the
                            // gyro-propagated, accel-corrected estimate and is
                            // strictly better than any instantaneous
                            // accel-derived level frame, so overwriting the
                            // whole quaternion here would inject the wave tilt
                            // error the filter has already rejected.
                            const float yaw_abs_rad =
                                wrapPi_(-mag_tilt_yaw_rad);

                            Eigen::Quaternionf q_bw =
                                impl_.mekf().quaternion_boat();
                            q_bw.normalize();

                            const Eigen::Quaternionf q_new =
                                boatQuatWithAbsoluteYaw_(q_bw, yaw_abs_rad);

                            if (q_new.coeffs().allFinite()) {
                                impl_.mekf().set_quaternion_boat(q_new);

                                last_mag_tilt_frame_yaw_rad_ =
                                    wrapPi_(mag_tilt_yaw_rad);

                                last_mag_startup_yaw_correction_rad_ =
                                    yaw_abs_rad;
                            }
                        }

                        Eigen::Vector3f hard_iron_uT;
                        mag_hard_iron_body_uT_ =
                            mag_auto_tuner_.getHardIronBodyUT(hard_iron_uT)
                                ? hard_iron_uT
                                : Eigen::Vector3f::Zero();

                        mag_ref_set_ = true;
                        syncLinearBlockGate_();
                    }
                }
            }
        }

        if (mag_ref_set_) {
            impl_.updateMag(mag_body_ned - mag_hard_iron_body_uT_);
        }
    }

    bool hasMagNorthLock() const noexcept {
        return mag_ref_set_;
    }

    // Body-frame hard-iron offset removed from the magnetometer stream.  Zero
    // unless Config::mag_estimate_hard_iron asked for it and the startup window
    // constrained it well enough to use.
    const Eigen::Vector3f& magHardIronBodyUT() const noexcept {
        return mag_hard_iron_body_uT_;
    }


    bool isLive() const {
        return stage_ == Stage::Live;
    }

    float freqHz() const {
        return impl_.getFreqHz();
    }

    float waveDirectionDeg() const {
        return impl_.getWaveDirectionDeg();
    }

    Eigen::Vector3f eulerNauticalDeg() const {
        return impl_.getEulerNautical();
    }

    const Eigen::Vector3f& displacementUpMeters() const {
        return displacement_up_m_;
    }

    const AdaptiveWaveDetrender3D::Output& displacementDetrend() const {
        return displacement_det_out_;
    }

    SeaStateFusionFilter_OU_III<trackerT>& raw() {
        return impl_;
    }

    const SeaStateFusionFilter_OU_III<trackerT>& raw() const {
        return impl_;
    }

    int magAcceptedCount() const noexcept {
        return mag_auto_tuner_.acceptedCount();
    }

    int magRejectedCount() const noexcept {
        return mag_auto_tuner_.rejectedCount();
    }

    float magAcceptedWindowSec() const noexcept {
        return mag_auto_tuner_.acceptedWindowSec();
    }

    float magEffectiveWeight() const noexcept {
        return mag_auto_tuner_.effectiveWeight();
    }

    float magTiltFrameYawDeg() const noexcept {
        return std::isfinite(last_mag_tilt_frame_yaw_rad_)
            ? last_mag_tilt_frame_yaw_rad_ * 57.29577951308232f
            : NAN;
    }

    float magStartupYawCorrectionDeg() const noexcept {
        return std::isfinite(last_mag_startup_yaw_correction_rad_)
            ? last_mag_startup_yaw_correction_rad_ * 57.29577951308232f
            : NAN;
    }

private:
    enum class Stage {
        Uninitialized,
        Warming,
        Live
    };

    struct Vec3LPF {
        Eigen::Vector3f state = Eigen::Vector3f::Zero();
        bool initialized = false;

        void reset() {
            state.setZero();
            initialized = false;
        }

        Eigen::Vector3f step(const Eigen::Vector3f& x,
                             float dt,
                             float tau_sec)
        {
            if (!x.allFinite()) return state;

            const float tau = std::max(1.0e-3f, tau_sec);
            const float alpha = 1.0f - std::exp(-dt / tau);

            if (!initialized) {
                state = x;
                initialized = true;
                return state;
            }

            state += alpha * (x - state);
            return state;
        }
    };

    using StartupTiltObserver = seastate::common::StartupTiltObserver;

    void resetTiltInit_() {
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_slow_lpf_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
    }

    bool linearBlockAllowed_() const {
        if (!cfg_.enable_linear_block) {
            return false;
        }
    
        if (!cfg_.with_mag) {
            return true;
        }
    
        if (!cfg_.require_mag_lock_for_linear_block) {
            return true;
        }
    
        return mag_ref_set_;
    }
    
    void syncLinearBlockGate_() {
        impl_.enableLinearBlock(linearBlockAllowed_());
    }

    static float wrapPi_(float a)
    {
        constexpr float PI_F = 3.14159265358979323846f;
        constexpr float TWO_PI_F = 2.0f * PI_F;

        if (!std::isfinite(a)) return NAN;

        while (a > PI_F) {
            a -= TWO_PI_F;
        }

        while (a <= -PI_F) {
            a += TWO_PI_F;
        }

        return a;
    }

    static float yawFromBoatQuatRad_(const Eigen::Quaternionf& q_bw_in) {
        if (!q_bw_in.coeffs().allFinite()) return NAN;

        Eigen::Quaternionf q_bw = q_bw_in;
        const float qn = q_bw.norm();

        if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
            return NAN;
        }

        q_bw.normalize();

        const Eigen::Matrix3f R = q_bw.toRotationMatrix();

        const float c = R(0, 0);
        const float s = R(1, 0);

        if (!std::isfinite(c) || !std::isfinite(s)) {
            return NAN;
        }

        return std::atan2(s, c);
    }

    // Tilt part of a BODY->WORLD quaternion, with the heading divided out.
    // Invariant under q_bw -> Rz(psi) q_bw, which is what makes it usable as a
    // mag accumulation frame that cannot leak the MEKF's arbitrary yaw.
    static Eigen::Quaternionf yawRemovedBoatQuat_(
        const Eigen::Quaternionf& q_bw_in)
    {
        if (!q_bw_in.coeffs().allFinite()) {
            return Eigen::Quaternionf::Identity();
        }

        Eigen::Quaternionf q_bw = q_bw_in;
        const float qn = q_bw.norm();

        if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
            return Eigen::Quaternionf::Identity();
        }

        q_bw.normalize();

        const float yaw = yawFromBoatQuatRad_(q_bw);

        if (!std::isfinite(yaw)) {
            return Eigen::Quaternionf::Identity();
        }

        const Eigen::Quaternionf q_yaw_inv(
            Eigen::AngleAxisf(-yaw, Eigen::Vector3f::UnitZ()));

        Eigen::Quaternionf q_tilt = q_yaw_inv * q_bw;
        q_tilt.normalize();

        if (!q_tilt.coeffs().allFinite()) {
            return Eigen::Quaternionf::Identity();
        }

        return q_tilt;
    }

    // Rewrites heading only, keeping the estimated tilt untouched.
    static Eigen::Quaternionf boatQuatWithAbsoluteYaw_(
        const Eigen::Quaternionf& q_bw_in,
        float yaw_abs_rad)
    {
        if (!std::isfinite(yaw_abs_rad)) {
            return q_bw_in;
        }

        const Eigen::Quaternionf q_tilt =
            yawRemovedBoatQuat_(q_bw_in);

        const Eigen::Quaternionf q_yaw(
            Eigen::AngleAxisf(
                yaw_abs_rad,
                Eigen::Vector3f::UnitZ()));

        Eigen::Quaternionf q_out = q_yaw * q_tilt;
        q_out.normalize();

        if (!q_out.coeffs().allFinite()) {
            return q_bw_in;
        }

        return q_out;
    }

    static Eigen::Quaternionf tiltOnlyQuatFromBoatQuat_(
        const Eigen::Quaternionf& q_bw_in)
    {
        return yawRemovedBoatQuat_(q_bw_in);
    }

private:
    Config cfg_{};
    SeaStateFusionFilter_OU_III<trackerT> impl_{false};

    bool begun_ = false;

    Stage stage_ = Stage::Uninitialized;
    float t_ = 0.0f;

    typename SeaStateFusionFilter_OU_III<trackerT>::StartupStage last_impl_startup_stage_ =
        SeaStateFusionFilter_OU_III<trackerT>::StartupStage::Cold;

    Eigen::Vector3f last_acc_body_ned_  = Eigen::Vector3f::Zero();
    Eigen::Vector3f last_gyro_body_ned_ = Eigen::Vector3f::Zero();
    bool have_last_imu_ = false;

    bool mag_ref_set_ = false;
    Eigen::Vector3f mag_hard_iron_body_uT_ = Eigen::Vector3f::Zero();
    MagAutoTuner mag_auto_tuner_{};

    float last_mag_sample_t_ = NAN;

    float last_mag_tilt_frame_yaw_rad_ = NAN;
    float last_mag_startup_yaw_correction_rad_ = NAN;

    AdaptiveWaveDetrender3D displacement_detrender_{};
    AdaptiveWaveDetrender3D::Output displacement_det_out_{};
    Eigen::Vector3f displacement_up_m_ = Eigen::Vector3f::Zero();

    Vec3LPF gravity_gate_acc_lpf_{};
    float   mag_gravity_good_sec_ = 0.0f;
    float   mag_init_eligible_t0_ = NAN;

    StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF             bootstrap_gravity_slow_lpf_{};
    float               bootstrap_gravity_good_sec_ = 0.0f;
};
