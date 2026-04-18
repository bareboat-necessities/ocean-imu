#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  Released under the MIT License

  SeaStateFusionFilter_OU_II

  Marine Inertial Navigational System (INS) Filter for IMU

  Combines multiple real-time estimators into a cohesive ocean-state tracker:

    • Quaternion-based attitude and linear motion estimation via Kalman3D_Wave_OU_II

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
      σₐ·τ² regularization law to stabilize displacement drift correction.

  Where
  – τ (tau):  OU process time constant ≈ ½ · T  (half the dominant period of acceleration)
  – σₐ:       Stationary acceleration standard deviation, EWMA-tracked online
  – R_p0:     Pseudo-measurement noise controlling p drift suppression
  – R_v0:     Pseudo-measurement noise controlling v drift suppression
  – R_p0_xy:  Reduced in X/Y (anisotropic weighting for vertical-dominant seas)

  Adaptive update: exponential smoothing toward targets over ADAPT_TAU_SEC

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
#include <memory>
#include <algorithm>

#include "freq/FirstOrderIIRSmoother.h"
#include "freq/FrequencyTrackerPolicy.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/MagAutoTuner.h"
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "wave_dir/KalmanWaveDirection.h"
#include "wave_dir/WaveDirectionDetector.h"
#include "detrend/AdaptiveWaveDetrender3D.h"

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

constexpr float MIN_TAU_S     = 0.02f;  // sec
constexpr float MAX_TAU_S     = 3.0f;   // sec
constexpr float MAX_SIGMA_A   = 6.0f;
constexpr float MIN_R_p0_std  = 0.05f;
constexpr float MAX_R_p0_std  = 18.0f;
constexpr float MIN_R_v0_std  = 0.01f;
constexpr float MAX_R_v0_std  = 6.0f;

constexpr float ADAPT_TAU_SEC              = 1.5f;
constexpr float ADAPT_EVERY_SECS           = 0.1f;
constexpr float ADAPT_R_p0_MULT            = 5.0f;   // dimensionless
constexpr float ADAPT_R_v0_MULT            = 5.0f;   // dimensionless
constexpr float ONLINE_TUNE_WARMUP_SEC     = 5.0f;
constexpr float MAG_DELAY_SEC              = 7.0f;

// Frequency smoother dt (SeaStateFusionFilter_OU_II is designed for 200 Hz)
constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;

struct TuneState {
    float tau_applied      = 1.1f;   // s
    float sigma_applied    = 1e-2f;  // m/s²
    float R_p0_std_applied = 0.1f;   // m
    float R_v0_std_applied = 0.1f;   // m/s
};

// Unified SeaState fusion filter
template<TrackerType trackerT>
class SeaStateFusionFilter_OU_II {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using TrackingPolicy = TrackerPolicy<trackerT>;

    enum class StartupStage {
        Cold,        // just booted or just had a big tilt reset
        TunerWarm,   // MEKF + freq running, tuner collecting stats
        Live         // tuner is trusted; full adaptation & extras allowed
    };

    explicit SeaStateFusionFilter_OU_II(bool with_mag = true)
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
        mekf_ = std::make_unique<Kalman3D_Wave_OU_II<float>>(sigma_a, sigma_g, sigma_m);
        enterCold_();      // applies freeze + warmup Racc + disables linear block
        apply_ou_tune_();
        mekf_->set_exact_att_bias_Qd(true);
    }

    void initialize_ext(const Eigen::Vector3f& sigma_a,
                        const Eigen::Vector3f& sigma_g,
                        const Eigen::Vector3f& sigma_m,
                        float Pq0, float Pb0,
                        float b0, float R_p0_var_init, float R_v0_var_init,
                        float gravity_magnitude)
    {
        mekf_ = std::make_unique<Kalman3D_Wave_OU_II<float>>(
            sigma_a, sigma_g, sigma_m, Pq0, Pb0, b0, R_p0_var_init, R_v0_var_init, gravity_magnitude);
        enterCold_();      // applies freeze + warmup Racc + disables linear block
        apply_ou_tune_();
        mekf_->set_exact_att_bias_Qd(true);
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

        // Keep BODY components around for direction/sign.
        const float a_x_body = acc.x();
        const float a_y_body = acc.y();

        // BODY-Z-based proxy used by the tracker/sign logic.
        // This is NOT true world/inertial vertical acceleration; it is only a
        // body-Z residual that behaves like up-positive vertical motion when the
        // platform is near-level:
        //   acc.z() ~ -g at rest  => proxy ~ 0
        const float a_z_body_proxy = acc.z() + g_std;

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

        // Up-positive BODY-Z proxy used by tracker/tuner/sign logic.
        // Not true world vertical unless the platform is close to level.
        a_body_z_up_proxy_ = -a_z_body_proxy;

        // LPF on BODY-Z proxy for tracker input
        const float a_body_z_up_lp = freq_input_lpf_.step(a_body_z_up_proxy_, dt);

        // Raw freq from tracker
        const float f_tracker = static_cast<float>(tracker_policy_.run(a_body_z_up_lp, dt));
        f_raw = f_tracker;

        // Stillness detector also sees the same BODY-Z proxy.
        const float f_after_still = freq_stillness_.step(a_body_z_up_lp, dt, f_tracker);

        // Fast & slow smoothed frequencies
        float f_fast = freq_fast_smoother_.update(f_after_still);
        float f_slow = freq_slow_smoother_.update(f_fast);

        f_fast = std::min(std::max(f_fast, min_freq_hz_), max_freq_hz_);
        f_slow = std::min(std::max(f_slow, min_freq_hz_), max_freq_hz_);

        freq_hz_      = f_fast;   // demod / direction
        freq_hz_slow_ = f_slow;   // tuner / moments

        if (enable_tuner_) {
            update_tuner(dt, a_body_z_up_proxy_, f_after_still);
        }

        // Keep linear-block R_p0 / R_v0 tuning responsive in Live mode instead of
        // waiting for slow adaptation cadence.
        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }

        const float omega = 2.0f * static_cast<float>(M_PI) * freq_hz_;

        // Direction filters run on BODY accel; sign uses the same BODY-Z proxy.
        dir_filter_.update(a_x_body, a_y_body, omega, dt);
        dir_sign_state_ = dir_sign_.update(a_x_body, a_y_body, a_body_z_up_proxy_, dt);
    }

    // Magnetometer correction
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
            (static_cast<float>(time_) - first_mag_update_time_) > 1.0f)
        {
            accel_bias_locked_ = false;

            if (freeze_acc_bias_until_live_ && startup_stage_ == StartupStage::Live) {
                mekf_->set_acc_bias_updates_enabled(true);

                if (warmup_Racc_active_) {
                    if (Racc_nominal_std_.allFinite() && Racc_nominal_std_.maxCoeff() > 0.0f) {
                        mekf_->set_Racc_std(Racc_nominal_std_);
                        warmup_Racc_active_ = false;
                    }
                }
            }
        }
    }

    void setWithMag(bool with_mag) { with_mag_ = with_mag; }

    // Anisotropy configuration (runtime)
    // P-factor scales horizontal vs vertical stationary std of a_w.
    // R_p0 XY factor scales pseudo-measurement noise in X/Y vs Z.
    void setPFactor(float p) {
        if (std::isfinite(p) && p > 0.0f) P_factor_ = p;
    }

    void setR_p0_XYFactor(float k) {
        if (std::isfinite(k)) {
            R_p0_xy_factor_ = std::min(std::max(k, 0.0f), 1.0f);
        }
    }

    void setTauCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) tau_coeff_ = c;
    }

    void setSigmaCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) sigma_coeff_ = c;
    }

    void setR_p0_Coeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_p0_coeff_;
            R_p0_coeff_ = c;

            if (std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.R_p0_std_applied) && tune_.R_p0_std_applied > 0.0f) {
                    tune_.R_p0_std_applied *= scale;
                }
                if (std::isfinite(R_p0_std_target_) && R_p0_std_target_ > 0.0f) {
                    R_p0_std_target_ *= scale;
                }
                if (enable_linear_block_) {
                    apply_R_p0_tune_();
                }
            }
        }
    }

    void setR_v0_Coeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_v0_coeff_;
            R_v0_coeff_ = c;

            if (std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.R_v0_std_applied) && tune_.R_v0_std_applied > 0.0f) {
                    tune_.R_v0_std_applied *= scale;
                }
                if (std::isfinite(R_v0_std_target_) && R_v0_std_target_ > 0.0f) {
                    R_v0_std_target_ *= scale;
                }
                if (enable_linear_block_) {
                    apply_R_v0_tune_();
                }
            }
        }
    }

    void setAccNoiseFloorSigma(float s) {
        if (std::isfinite(s) && s > 0.0f) acc_noise_floor_sigma_ = s;
    }

    float getAccNoiseFloorSigma() const noexcept { return acc_noise_floor_sigma_; }

    // Configure LPF on BODY-Z proxy for tracker input
    void setFreqInputCutoffHz(float fc) { freq_input_lpf_.setCutoff(fc); }

    void enableClamp(bool flag = true) { enable_clamp_ = flag; }
    void enableTuner(bool flag = true) { enable_tuner_ = flag; }

    // Enable/disable use of the extended linear block [v,p,a_w] in Kalman3D_Wave_OU_II.
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

    void setR_p0_Bounds(float min_R_p0_std, float max_R_p0_std) {
        if (!std::isfinite(min_R_p0_std) || !std::isfinite(max_R_p0_std)) return;
        if (min_R_p0_std <= 0.0f || max_R_p0_std <= min_R_p0_std) return;
        MIN_R_p0_std_ = min_R_p0_std;
        MAX_R_p0_std_ = max_R_p0_std;
    }

    void setR_v0_Bounds(float min_R_v0_std, float max_R_v0_std) {
        if (!std::isfinite(min_R_v0_std) || !std::isfinite(max_R_v0_std)) return;
        if (min_R_v0_std <= 0.0f || max_R_v0_std <= min_R_v0_std) return;
        MIN_R_v0_std_ = min_R_v0_std;
        MAX_R_v0_std_ = max_R_v0_std;
    }

    void setAdaptationTimeConstants(float tau_sec) {
        if (std::isfinite(tau_sec) && tau_sec > 0.0f) adapt_tau_sec_ = tau_sec;
    }

    void setAdaptationUpdatePeriod(float every_sec) {
        if (std::isfinite(every_sec) && every_sec > 0.0f) adapt_every_secs_ = every_sec;
    }

    void setOnlineTuneWarmupSec(float warmup_sec) {
        if (std::isfinite(warmup_sec) && warmup_sec >= 0.0f) online_tune_warmup_sec_ = warmup_sec;
    }

    void setMagDelaySec(float delay_sec) {
        if (std::isfinite(delay_sec) && delay_sec >= 0.0f) mag_delay_sec_ = delay_sec;
    }

    void setFreezeAccBiasUntilLive(bool en) { freeze_acc_bias_until_live_ = en; }

    void setWarmupRaccStd(float r) {
        if (std::isfinite(r) && r > 0.0f) Racc_warmup_std_ = r;
    }

    void setNominalRaccStd(const Eigen::Vector3f& r) { Racc_nominal_std_ = r; }

    inline float getFreqHz()           const noexcept { return freq_hz_; }
    inline float getFreqSlowHz()       const noexcept { return freq_hz_slow_; }
    inline float getFreqRawHz()        const noexcept { return f_raw; }
    inline float getTauApplied()       const noexcept { return mekf_ ? mekf_->get_aw_time_constant() : NAN; }
    inline float getSigmaApplied()     const noexcept { return mekf_ ? mekf_->get_aw_stationary_std().z() : NAN; }
    inline float getR_p0_std_applied() const noexcept { return mekf_ ? mekf_->get_Rp0_noise_std().z() : NAN; }
    inline float getR_v0_std_applied() const noexcept { return mekf_ ? mekf_->get_Rv0_noise_std().z() : NAN; }
    inline float getTauTarget()        const noexcept { return tau_target_; }
    inline float getSigmaTarget()      const noexcept { return sigma_target_; }
    inline float getR_p0_std_target()  const noexcept { return R_p0_std_target_; }
    inline float getR_v0_std_target()  const noexcept { return R_v0_std_target_; }

    inline float getPeriodSec() const noexcept {
        return (freq_hz_slow_ > 1e-6f) ? 1.0f / freq_hz_slow_ : NAN;
    }

    inline float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }

    // Returns the BODY-Z-based up-positive proxy used by tracker/tuner logic.
    // This is not a true world/inertial vertical acceleration estimate.
    inline float getAccelVertical() const noexcept { return a_body_z_up_proxy_; }

    inline float getHeaveAbs() const noexcept {
        if (!mekf_) return NAN;
        return std::fabs(mekf_->get_position().z());
    }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
        constexpr float C_HS = 2.0f * std::sqrt(2.0f) / (M_PI * M_PI);
        return C_HS * sigma * tau * tau / 2.0f;
    }

    float getVerticalSpeedEnvelopeMps(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!(tau > 1e-6f) || !std::isfinite(tau) || !std::isfinite(sigma)) return NAN;
        constexpr float K = std::sqrt(2.0f) / M_PI;
        const float v_env = K * sigma * tau;
        return std::isfinite(v_env) ? v_env : NAN;
    }

    inline WaveDirection getDirSignState() const noexcept { return dir_sign_state_; }
    inline float getWaveDirectionDeg() const noexcept { return dir_filter_.getDirectionDegrees(); }

    inline auto& mekf() noexcept { return *mekf_; }
    inline const auto& mekf() const noexcept { return *mekf_; }

    inline KalmanWaveDirection& dir() noexcept { return dir_filter_; }
    inline const KalmanWaveDirection& dir() const noexcept { return dir_filter_; }

    inline WaveDirectionDetector<float>& dir_sign() noexcept { return dir_sign_; }
    inline const WaveDirectionDetector<float>& dir_sign() const noexcept { return dir_sign_; }

private:

    struct FreqInputLPF {
        float state       = 0.0f;
        float fc_hz       = 1.0f;
        bool  initialized = false;

        void setCutoff(float fc) {
            if (std::isfinite(fc) && fc > 0.0f) fc_hz = fc;
        }

        float step(float x, float dt) {
            const float alpha = std::exp(-2.0f * static_cast<float>(M_PI) * fc_hz * dt);

            if (!initialized) {
                state = x;
                initialized = true;
                return state;
            }
            state = (1.0f - alpha) * x + alpha * state;
            return state;
        }
    };

    struct StillnessAdapter {
        float energy_ema      = 0.0f;
        float energy_alpha    = 0.05f;
        float energy_thresh   = 8e-4f;

        float still_time_sec  = 0.0f;
        float still_thresh_s  = 2.0f;

        float relax_tau_sec   = 1.0f;
        float target_freq_hz  = MIN_FREQ_HZ;

        bool  freq_init       = false;
        float freq_state      = FREQ_GUESS;

        bool  last_is_still   = false;

        void setTargetFreqHz(float f) {
            if (std::isfinite(f) && f > 0.0f) target_freq_hz = f;
        }

        // a_z_body_proxy_lp:
        //   low-passed BODY-Z-based up-positive proxy acceleration (m/s²)
        //   used for tracker/stillness logic. Not true world vertical.
        float step(float a_z_body_proxy_lp, float dt, float freq_in) {
            if (!(dt > 0.0f) || !std::isfinite(freq_in)) return freq_in;

            if (!freq_init || !std::isfinite(freq_state)) {
                freq_state = freq_in;
                freq_init = true;
            }

            const float a_norm      = a_z_body_proxy_lp / g_std;
            const float inst_energy = a_norm * a_norm;

            energy_ema = (1.0f - energy_alpha) * energy_ema + energy_alpha * inst_energy;
            const bool is_still = (energy_ema < energy_thresh);
            last_is_still = is_still;

            if (is_still) {
                still_time_sec += dt;
                if (still_time_sec > 60.0f) still_time_sec = 60.0f;

                if (still_time_sec > still_thresh_s) {
                    const float relax_alpha = 1.0f - std::exp(-dt / relax_tau_sec);
                    freq_state += relax_alpha * (target_freq_hz - freq_state);
                } else {
                    freq_state = freq_in;
                }
            } else {
                still_time_sec = 0.0f;
                freq_state = freq_in;
            }
            return freq_state;
        }

        bool  isStill()      const { return last_is_still; }
        float getStillTime() const { return still_time_sec; }
        float getEnergyEma() const { return energy_ema; }
    };

    void apply_ou_tune_() {
        if (!mekf_) return;
        mekf_->set_aw_time_constant(tune_.tau_applied);

        const float sigma_floor = std::max(0.05f, acc_noise_floor_sigma_);
        const float sZ = std::max(sigma_floor, tune_.sigma_applied);
        const float sH = sZ * P_factor_;
        mekf_->set_aw_stationary_std(Eigen::Vector3f(sH, sH, sZ));
    }

    void apply_R_p0_tune_(float rp_scale = 1.0f) {
        if (!mekf_) return;
        const float p = (std::isfinite(rp_scale) && rp_scale > 0.0f) ? std::min(rp_scale, 1.0f) : 1.0f;
        const float R_p0_b = std::min(std::max(tune_.R_p0_std_applied, MIN_R_p0_std_), MAX_R_p0_std_);
        const float rp_xy = R_p0_b * p * R_p0_xy_factor_;
        mekf_->set_Rp0_noise_std(Eigen::Vector3f(rp_xy, rp_xy, R_p0_b * p));
    }

    void apply_R_v0_tune_(float rv_scale = 1.0f) {
        if (!mekf_) return;
        const float p = (std::isfinite(rv_scale) && rv_scale > 0.0f) ? std::min(rv_scale, 1.0f) : 1.0f;
        const float R_v0_b = std::min(std::max(tune_.R_v0_std_applied, MIN_R_v0_std_), MAX_R_v0_std_);
        mekf_->set_Rv0_noise_std(Eigen::Vector3f::Constant(R_v0_b * p));
    }

    void update_tuner(float dt, float a_body_z_up_proxy, float freq_hz_for_tuner) {
        tuner_.update(dt, a_body_z_up_proxy, freq_hz_for_tuner);

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

        float f_tune = tuner_.getFrequencyHz();
        if (!std::isfinite(f_tune) || f_tune < min_freq_hz_) f_tune = min_freq_hz_;
        if (f_tune > max_freq_hz_) f_tune = max_freq_hz_;

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
            tau_target_   = std::min(std::max(tau_raw, min_tau_s_), max_tau_s_);
            sigma_target_ = std::min(sigma_wave * sigma_coeff_, max_sigma_a_);
        } else {
            tau_target_   = tau_raw;
            sigma_target_ = sigma_wave;
        }

        if (!tuner_.isVarReady()) {
            sigma_target_ = std::max(sigma_target_, std::max(0.05f, acc_noise_floor_sigma_));
        }

        float R_p0_raw = R_p0_coeff_ * sigma_target_ * tau_target_ * tau_target_;
        float R_v0_raw = R_v0_coeff_ * sigma_target_ * tau_target_;

        if (enable_clamp_) {
            R_p0_std_target_ = std::min(std::max(R_p0_raw, MIN_R_p0_std_), MAX_R_p0_std_);
            R_v0_std_target_ = std::min(std::max(R_v0_raw, MIN_R_v0_std_), MAX_R_v0_std_);
        } else {
            R_p0_std_target_ = R_p0_raw;
            R_v0_std_target_ = R_v0_raw;
        }

        adapt_mekf(dt, tau_target_, sigma_target_, R_p0_std_target_, R_v0_std_target_);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float R_p0_t, float R_v0_t) {
        const float alpha = 1.0f - std::exp(-dt / adapt_tau_sec_);

        const float R_p0_sec   = ADAPT_R_p0_MULT * tau_t;
        const float R_v0_sec   = ADAPT_R_v0_MULT * tau_t;
        const float alpha_R_p0 = 1.0f - std::exp(-dt / R_p0_sec);
        const float alpha_R_v0 = 1.0f - std::exp(-dt / R_v0_sec);

        tune_.tau_applied      += alpha      * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied    += alpha      * (sigma_t - tune_.sigma_applied);
        tune_.R_p0_std_applied += alpha_R_p0 * (R_p0_t  - tune_.R_p0_std_applied);
        tune_.R_v0_std_applied += alpha_R_v0 * (R_v0_t  - tune_.R_v0_std_applied);

        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {
            if (tuner_.isFreqReady()) {
                apply_ou_tune_();
            }
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
                apply_R_p0_tune_();
                apply_R_v0_tune_();
            }
            last_adapt_time_sec_ = time_;
        }
    }

    void resetTrackingState_() {
        tracker_policy_ = TrackingPolicy{};
        freq_input_lpf_ = FreqInputLPF{};
        freq_stillness_ = StillnessAdapter{};
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);

        tuner_.reset();

        freq_fast_smoother_ = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 3.5f);
        freq_slow_smoother_ = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 10.0f);

        freq_hz_      = FREQ_GUESS;
        freq_hz_slow_ = FREQ_GUESS;
        f_raw         = FREQ_GUESS;

        dir_filter_ = KalmanWaveDirection(2.0f * static_cast<float>(M_PI) * FREQ_GUESS);
        dir_sign_state_ = UNCERTAIN;

        last_adapt_time_sec_ = time_;
    }

    void enterCold_() {
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;
        mekf_->set_linear_block_enabled(false);

        accel_bias_locked_   = with_mag_;
        mag_updates_applied_ = 0;
        first_mag_update_time_ = NAN;

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
        mekf_->set_linear_block_enabled(enable_linear_block_);

        if (freeze_acc_bias_until_live_) {
            const bool allow_bias = !accel_bias_locked_;
            mekf_->set_acc_bias_updates_enabled(allow_bias);

            if (warmup_Racc_active_ &&
                Racc_nominal_std_.allFinite() &&
                Racc_nominal_std_.maxCoeff() > 0.0f)
            {
                mekf_->set_Racc_std(Racc_nominal_std_);
            }
            warmup_Racc_active_ = false;
        }

        apply_ou_tune_();
        if (enable_linear_block_) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }
    }

    StartupStage startup_stage_ = StartupStage::Cold;
    float        startup_stage_t_ = 0.0f;

    bool  freeze_acc_bias_until_live_ = true;
    float Racc_warmup_std_            = 0.5f;
    bool  warmup_Racc_active_         = false;
    Eigen::Vector3f Racc_nominal_std_ = Eigen::Vector3f::Constant(0.0f);

    bool accel_bias_locked_ = true;
    int  mag_updates_applied_ = 0;
    static constexpr int MAG_UPDATES_TO_UNLOCK = 140;

    bool   with_mag_;
    double time_;
    double last_adapt_time_sec_;

    float first_mag_update_time_ = NAN;

    float tilt_over_limit_sec_ = 0.0f;
    float tilt_reset_cooldown_sec_ = 0.0f;

    float freq_hz_      = FREQ_GUESS;
    float freq_hz_slow_ = FREQ_GUESS;
    float f_raw         = FREQ_GUESS;

    float a_body_z_up_proxy_ = 0.0f;  // BODY-Z-based up-positive proxy used by tracker/tuner logic.

    bool enable_clamp_ = true;
    bool enable_tuner_ = true;

    bool enable_linear_block_ = true;

    float min_freq_hz_            = MIN_FREQ_HZ;
    float max_freq_hz_            = MAX_FREQ_HZ;
    float min_tau_s_              = MIN_TAU_S;
    float max_tau_s_              = MAX_TAU_S;
    float max_sigma_a_            = MAX_SIGMA_A;
    float MIN_R_p0_std_           = MIN_R_p0_std;
    float MAX_R_p0_std_           = MAX_R_p0_std;
    float MIN_R_v0_std_           = MIN_R_v0_std;
    float MAX_R_v0_std_           = MAX_R_v0_std;
    float adapt_tau_sec_          = ADAPT_TAU_SEC;
    float adapt_every_secs_       = ADAPT_EVERY_SECS;
    float online_tune_warmup_sec_ = ONLINE_TUNE_WARMUP_SEC;
    float mag_delay_sec_          = MAG_DELAY_SEC;

    float R_p0_xy_factor_ = 0.25f;
    float P_factor_       = 1.5f;

    TrackingPolicy               tracker_policy_{};
    FirstOrderIIRSmoother<float> freq_fast_smoother_{FREQ_SMOOTHER_DT, 3.5f};
    FirstOrderIIRSmoother<float> freq_slow_smoother_{FREQ_SMOOTHER_DT, 10.0f};
    SeaStateAutoTuner            tuner_;
    TuneState                    tune_;

    float tau_target_      = NAN;
    float sigma_target_    = NAN;
    float R_p0_std_target_ = NAN;
    float R_v0_std_target_ = NAN;

    float acc_noise_floor_sigma_ = ACC_NOISE_FLOOR_SIGMA_DEFAULT;

    float R_p0_coeff_  = 1.6f;
    float R_v0_coeff_  = 1.4f;
    float tau_coeff_   = 1.5f;
    float sigma_coeff_ = 0.85f;

    std::unique_ptr<Kalman3D_Wave_OU_II<float>> mekf_;
    KalmanWaveDirection dir_filter_{2.0f * static_cast<float>(M_PI) * FREQ_GUESS};

    FreqInputLPF     freq_input_lpf_;
    StillnessAdapter freq_stillness_;

    WaveDirectionDetector<float> dir_sign_{0.002f, 0.005f};
    WaveDirection                dir_sign_state_ = UNCERTAIN;
};

template<TrackerType trackerT>
class SeaStateFusion_OU_II {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Config {
        bool with_mag = true;

        float mag_delay_sec          = MAG_DELAY_SEC;
        float online_tune_warmup_sec = ONLINE_TUNE_WARMUP_SEC;

        bool  freeze_acc_bias_until_live = true;
        float Racc_warmup_std = 0.5f;

        Eigen::Vector3f sigma_a = Eigen::Vector3f(0.2f, 0.2f, 0.2f);
        Eigen::Vector3f sigma_g = Eigen::Vector3f(0.01f, 0.01f, 0.01f);
        Eigen::Vector3f sigma_m = Eigen::Vector3f(0.3f, 0.3f, 0.3f);

        // mag-init policy:
        // wait a bit for tilt to settle, then average only a short stable window.
        float mag_init_min_mag_norm   = 1e-3f;

        bool enable_displacement_detrend = false;
        bool use_custom_displacement_detrend_cfg = false;
        AdaptiveWaveDetrender3D::Config displacement_detrend_cfg{};
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;

        begun_ = true;
        stage_ = Stage::Uninitialized;
        t_ = 0.0f;

        mag_ref_set_ = false;
        MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min = cfg_.mag_init_min_mag_norm;
        mag_auto_tuner_.setConfig(mag_cfg);

        tilt_init_acc_sum_.setZero();
        tilt_init_count_ = 0;

        last_acc_body_ned_.setZero();
        last_gyro_body_ned_.setZero();
        have_last_imu_ = false;

        impl_.setWithMag(cfg_.with_mag);
        impl_.setFreezeAccBiasUntilLive(cfg_.freeze_acc_bias_until_live);
        impl_.setWarmupRaccStd(cfg_.Racc_warmup_std);
        impl_.setMagDelaySec(cfg_.mag_delay_sec);
        impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);

        impl_.initialize(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m);
        last_impl_startup_stage_ = impl_.getStartupStage();

        impl_.setNominalRaccStd(cfg_.sigma_a);

        displacement_up_m_.setZero();
        displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

        if (cfg_.enable_displacement_detrend) {
            if (cfg_.use_custom_displacement_detrend_cfg) {
                displacement_detrender_.setConfig(cfg_.displacement_detrend_cfg);
            } else {
                AdaptiveWaveDetrender3D::Config dcfg;
                dcfg.init_wave_freq_hz = FREQ_GUESS;
                dcfg.min_wave_freq_hz  = 0.02f;
                dcfg.max_wave_freq_hz  = 1.20f;
                dcfg.baseline_cutoff_fraction = 0.25f;
                dcfg.min_baseline_cutoff_hz   = 0.003f;
                dcfg.max_baseline_cutoff_hz   = 0.25f;
                dcfg.freq_smooth_tau_s = 12.0f;
                dcfg.slope_lpf_tau_s   = 0.20f;
                dcfg.slope_rms_tau_s   = 8.0f;
                dcfg.freq_learn_axis   = 2;
                dcfg.threshold_rms_fraction  = 0.15f;
                dcfg.min_slope_threshold_abs = 0.002f;
                dcfg.max_slope_threshold_abs = 1.0e9f;
                dcfg.startup_hold_s      = 2.0f;
                dcfg.freq_timeout_cycles = 3.0f;
                dcfg.enable_wave_cleanup    = true;
                dcfg.cleanup_cutoff_fraction = 1.0f;
                dcfg.min_cleanup_cutoff_hz  = 0.003f;
                dcfg.max_cleanup_cutoff_hz  = 0.50f;
                dcfg.cleanup_stages         = 2;
                dcfg.min_dt_s = 1.0e-4f;
                dcfg.max_dt_s = 0.25f;
                dcfg.output_abs_limit = 0.0f;
                displacement_detrender_.setConfig(dcfg);
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

        // Stage 1: tilt learning from stable gravity samples.
        if (stage_ == Stage::Uninitialized) {
            if (isStableInitSample_(acc_body_ned, gyro_body_ned)) {
                tilt_init_acc_sum_ += acc_body_ned;
                ++tilt_init_count_;
            }

            constexpr int TILT_INIT_MIN_SAMPLES = 600; // @ 200 Hz
            if (tilt_init_count_ >= TILT_INIT_MIN_SAMPLES) {
                const Eigen::Vector3f acc_mean = tilt_init_acc_sum_ / static_cast<float>(tilt_init_count_);
                impl_.initialize_from_acc(acc_mean);
                stage_ = Stage::Warming;
                tilt_init_acc_sum_.setZero();
                tilt_init_count_ = 0;
            }
        }

        last_acc_body_ned_  = acc_body_ned;
        last_gyro_body_ned_ = gyro_body_ned;
        have_last_imu_      = true;

        if (stage_ != Stage::Uninitialized) {
            impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);

            const Eigen::Vector3f pos_ned_m = impl_.mekf().get_position();
            displacement_up_m_ = Eigen::Vector3f(pos_ned_m.x(), pos_ned_m.y(), -pos_ned_m.z());

            if (cfg_.enable_displacement_detrend) {
                const float wave_hz = impl_.getFreqHz();
                const bool ext_freq_valid =
                    isLive() &&
                    std::isfinite(wave_hz) &&
                    (wave_hz >= displacement_detrender_.config().min_wave_freq_hz) &&
                    (wave_hz <= displacement_detrender_.config().max_wave_freq_hz);

                displacement_det_out_ = displacement_detrender_.update(
                    displacement_up_m_, dt, wave_hz, ext_freq_valid);
            } else {
                displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
                displacement_det_out_.input = displacement_up_m_;
                displacement_det_out_.baseline_slow = Eigen::Vector3f::Zero();
                displacement_det_out_.wave_raw = displacement_up_m_;
                displacement_det_out_.wave_clean = displacement_up_m_;
            }
        }

        // If internal filter drops back to Cold during startup, forget mag init
        // and reacquire north from scratch later.
        const auto cur_stage = impl_.getStartupStage();
        if (cur_stage != last_impl_startup_stage_) {
            if (cur_stage == SeaStateFusionFilter_OU_II<trackerT>::StartupStage::Cold) {
                mag_ref_set_ = false;
                mag_auto_tuner_.reset();
            }
            last_impl_startup_stage_ = cur_stage;
        }

        if (stage_ == Stage::Warming && impl_.isAdaptiveLive()) {
            stage_ = Stage::Live;
        }
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!begun_ || !cfg_.with_mag) return;
        if (stage_ == Stage::Uninitialized) return;
        if (t_ < cfg_.mag_delay_sec) return;

        // Learn magnetic world reference from measurements only.
        if (!mag_ref_set_) {
            if (have_last_imu_ &&
                mag_auto_tuner_.addSample(last_acc_body_ned_, last_gyro_body_ned_, mag_body_ned))
            {
                Eigen::Vector3f mag_world_ref_uT;
                if (mag_auto_tuner_.getMagWorldRef(mag_world_ref_uT) &&
                    mag_world_ref_uT.allFinite() &&
                    mag_world_ref_uT.norm() > cfg_.mag_init_min_mag_norm)
                {
                    impl_.mekf().set_mag_world_ref(mag_world_ref_uT);
                    mag_ref_set_ = true;
                }
            }
        }

        if (mag_ref_set_) {
            impl_.updateMag(mag_body_ned);
        }
    }

    bool hasMagNorthLock() const noexcept { return mag_ref_set_; }

    bool  isLive() const { return stage_ == Stage::Live; }
    float freqHz() const { return impl_.getFreqHz(); }
    float waveDirectionDeg() const { return impl_.getWaveDirectionDeg(); }
    const Eigen::Vector3f& displacementUpMeters() const { return displacement_up_m_; }
    const AdaptiveWaveDetrender3D::Output& displacementDetrend() const { return displacement_det_out_; }

    SeaStateFusionFilter_OU_II<trackerT>& raw() { return impl_; }
    const SeaStateFusionFilter_OU_II<trackerT>& raw() const { return impl_; }

private:
    enum class Stage { Uninitialized, Warming, Live };

    static bool isStableInitSample_(const Eigen::Vector3f& acc_body,
                                    const Eigen::Vector3f& gyro_body)
    {
        constexpr float G = 9.80665f;
        constexpr float ACC_BAND = 0.08f * G;                    // ±8% around 1g
        constexpr float GYRO_MAX = 18.0f * float(M_PI) / 180.0f; // 18 deg/s

        if (!acc_body.allFinite() || !gyro_body.allFinite()) return false;
        const float acc_n = acc_body.norm();
        const float gyro_n = gyro_body.norm();
        return (std::fabs(acc_n - G) <= ACC_BAND) && (gyro_n <= GYRO_MAX);
    }

private:
    Config cfg_{};
    SeaStateFusionFilter_OU_II<trackerT> impl_{false};

    bool begun_ = false;

    Stage stage_ = Stage::Uninitialized;
    float t_ = 0.0f;

    typename SeaStateFusionFilter_OU_II<trackerT>::StartupStage last_impl_startup_stage_ =
        SeaStateFusionFilter_OU_II<trackerT>::StartupStage::Cold;

    // Last IMU sample for mag-init gating.
    Eigen::Vector3f last_acc_body_ned_  = Eigen::Vector3f::Zero();
    Eigen::Vector3f last_gyro_body_ned_ = Eigen::Vector3f::Zero();
    bool  have_last_imu_ = false;

    // One-shot mag-init state.
    bool mag_ref_set_ = false;
    MagAutoTuner mag_auto_tuner_{};
    Eigen::Vector3f tilt_init_acc_sum_ = Eigen::Vector3f::Zero();
    int tilt_init_count_ = 0;

    AdaptiveWaveDetrender3D displacement_detrender_{};
    AdaptiveWaveDetrender3D::Output displacement_det_out_{};
    Eigen::Vector3f displacement_up_m_ = Eigen::Vector3f::Zero();
};
