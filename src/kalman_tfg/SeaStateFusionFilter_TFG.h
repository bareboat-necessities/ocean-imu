#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Sea-state orchestration for the two-frame Lie-group filter.

    SCOPE. This is deliberately the core of what SeaStateFusionFilter_OU_III
    does, not all of it: startup staging, the sea-state tuner, the adaptation
    smoothers, and the update loop. The wave-direction estimator, the
    alternative frequency trackers (KalmANF, Aranovskiy, PLL), the displacement
    detrender and the wind-heel B' frame are NOT carried over. They are
    orthogonal to the estimator under test, and duplicating them would mean
    maintaining a second copy of logic whose behaviour is already established
    and evidenced for OU-III.

    What that costs is worth stating plainly: the TFG simulator is therefore
    not feature-comparable to the OU-III one on direction and frequency
    reporting, and any paired study has to compare the channels both filters
    actually implement.

    WHAT IS SHARED RATHER THAN COPIED
      - ::seastate::common::runStartupGravityInit and StartupTiltObserver, for
        the Cold-stage levelling;
      - ::seastate::common::adaptiveSmoothingHorizonSec, for the r_S channel;
      - SeaStateAutoTuner              (variance at the wave-band operating point)
      - WavePeriodEstimator            (zero-crossing T_z)
      - VerticalAccelComplementary     (its private Mahony observer)

    THE TUNER STAYS OUTSIDE THE GROUP STATE. tau, sigma_aw and r_S are
    hyperparameters here, not manifold states. The article's ablation shows
    nearly all the vertical benefit comes from adapting the integral
    regularization scale rather than the two OU parameters, so promoting them
    to stochastic states would add complexity without evidence it solves the
    problem that matters.

    WHY THE WAVE-PERIOD INPUT IS THE COMPLEMENTARY OBSERVER. The tuner's
    frequency channel must not be a function of the estimator it tunes, or the
    loop closes on itself. VerticalAccelComplementary runs its own small Mahony
    observer on raw gyro and accelerometer, so the period estimate is
    independent of the filter's attitude. That is the same choice OU-III
    settled on, and for the same reason.

    ------------------------------------------------------------------------
    ADAPTATION POLICY, brought to parity with SeaStateFusionFilter_OU_III
    ------------------------------------------------------------------------

    Four things here are not free parameters. Each was measured on OU-III and
    the reasoning carries over unchanged, because the tuner sits outside the
    estimator and does not know which filter it is driving.

    1. EXOGENEITY IS A TIMING PROPERTY, NOT JUST A SIGNAL CHOICE. Feeding the
       tuner from the complementary observer keeps its *inputs* independent of
       the filter. That is necessary and not sufficient: if the schedule
       smoothed during step k were also applied during step k, the covariance
       the MEKF uses at k would depend on y_k. So the smoother runs every
       sample, the cadence tick marks a candidate, and the commit happens at
       the top of the next update() -- before y_{k+1} reaches the MEKF. The
       active schedule at step k+1 is then measurable with respect to data
       through k.

    2. ADAPTATION CADENCE. Committing every sample is not "more adaptive", it
       just couples the schedule to the measurement stream 200 times a second.
       OU-III commits on a 0.1 s tick (ADAPT_EVERY_SECS); the EMA still sees
       every sample, so the trajectory is unchanged apart from being held
       piecewise constant between ticks.

    3. THE r_S FLOOR IS 0.15, NOT 0.4. On OU-III the 0.4 floor was not a
       safety limit, it was the binding constraint on every low-motion sea:
       the schedule asks for 0.24 m*s at the calibrated H_s = 0.27 m point, so
       the floor clipped it and a full sweep of the tuner multipliers left
       those scenarios constant to three decimals. Dropping it recovered -8.3%
       on the worst sea and cut the near-still H_s = 0.05 m case from 27.0% to
       17.6% of H_s. This filter inherited the old 0.4 and the same clipping.

    4. THE S=0 CADENCE IS SELF-SIMILAR IN tau, AND r_S MUST FOLLOW IT. One
       pseudo update has covariance r_S^2; at one update per T_S seconds the
       continuous-equivalent information rate goes as 1/(r_S^2 T_S). Scaling
       T_S with tau while holding r_S fixed therefore silently changes the
       regularization strength with sea state. Renormalizing the filter input
       by sqrt(T_0/T_S) holds the information rate, and turns the base
       sigma_aw*tau^3 schedule into an effective sigma_aw*tau^(5/2) one. The
       renormalized value is deliberately NOT re-clamped -- the smallest-sea
       point sits on the base floor and must be allowed below it once
       T_S > T_0, or the very clipping item 3 removes comes back.

    STARTUP: THE PROXY OWNS TILT AND MAGNETIC LEARNING. See StartupInitPolicy.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>

#include "kalman_common/SeaStateFusionFilterCommon.h"
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/WavePeriodEstimator.h"

namespace ocean_imu::tfg {

constexpr float MAG_DELAY_SEC = 7.0f;

struct TfgTuneState {
    float tau_applied   = 1.1f;
    float sigma_applied = 1e-2f;
    float RS_applied    = 0.5f;
};

template <typename MekfT = ocean_imu::kalman::Kalman3D_Wave_TFG<float>>
class SeaStateFusionFilter_TFG {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Mekf = MekfT;
    using Vector3f = Eigen::Vector3f;

    enum class StartupStage { Cold, TunerWarm, Live };

    /*
        Who solves the attitude the filter starts from.

        StagedMekf is the original behaviour, kept for ablation: the MEKF runs
        from the first sample and learns its own tilt while degraded -- linear
        block off, accelerometer bias frozen, Racc inflated -- and the
        magnetic reference is then captured in *that* attitude's frame. Those
        reads are the problem. The reference and the yaw gauge are locked
        once, so whatever tilt error the warming MEKF holds at that instant
        becomes a standing attitude and heading bias. It is also the moment
        the MEKF is least able to level: the linear block that absorbs orbital
        acceleration is switched off, so its accelerometer update is fighting
        the waves with nowhere to put them.

        MahonyProxy hands both jobs to the private Mahony observer this filter
        already runs for the wave-period channel. That observer is a pure
        function of raw gyro and accelerometer, and its correction corner sits
        an order of magnitude below the wave band, so it rejects orbital
        specific force by construction instead of chasing it. The whole
        measurement-only front end -- proxy, wave-period estimator, tuner --
        runs from the first sample without the MEKF, so by the time tilt has
        settled and the magnetometer has gauged north, the operating point has
        converged too. The MEKF is then seeded once and goes straight to Live,
        never occupying the degraded warmup configuration at all.

        This also closes the last exogeneity gap in the startup path: under
        MahonyProxy nothing the tuner or the magnetic gauge consumes is ever a
        function of the MEKF.
    */
    enum class StartupInitPolicy {
        StagedMekf,   // legacy: the MEKF learns its own tilt while warming
        MahonyProxy,  // default: the proxy owns tilt + mag, MEKF starts Live
    };

    struct Config {
        Vector3f sigma_a{Vector3f::Constant(0.5f)};
        float    gyro_noise_density = 0.005f;
        Vector3f sigma_m{Vector3f::Constant(0.1f)};
        bool     with_mag = true;
        float    mag_delay_sec = MAG_DELAY_SEC;
        bool     freeze_acc_bias_until_live = true;
        float    Racc_warmup_std = 0.5f;
        float    gravity_magnitude = 9.80665f;
        // Matches OU-III's ONLINE_TUNE_WARMUP_SEC. The old 20 s here was set
        // when the MEKF had to warm up alongside the tuner; under MahonyProxy
        // the front end converges on its own schedule and the MEKF is not
        // waiting on itself.
        float    online_tune_warmup_sec = 5.0f;
        StartupInitPolicy startup_init_policy = StartupInitPolicy::StagedMekf;

        // Handoff window. The lower bound stops a handoff from a barely-seeded
        // observer; the upper bound is the one that matters for robustness.
        //
        // Without a timeout the MEKF never starts if the front end never
        // declares itself ready -- and the wave-period channel legitimately
        // needs tens of seconds, longer in a low sea. A filter that silently
        // produces nothing is a worse failure than one that starts from a
        // slightly stale operating point, so past the timeout the handoff
        // proceeds on proxy tilt alone.
        float proxy_startup_min_sec     = 8.0f;
        float proxy_startup_timeout_sec = 150.0f;

        // Attitude uncertainty handed to the MEKF at the seed. The proxy is
        // good but not exact, and seeding an overconfident attitude would make
        // the first accelerometer updates fight a covariance that says there
        // is nothing to correct.
        float proxy_handoff_tilt_sigma_rad = 0.035f;
        float proxy_handoff_yaw_sigma_rad  = 0.087f;

        // The proxy doubles as the attitude seed, so its gyro-bias integrator
        // must be on. With two_ki = 0 a constant bias b leaves a standing tilt
        // of about 2b/two_kp -- 0.71 deg at 0.05 deg/s. Everything else the
        // observer feeds is high-passed and never noticed it; an attitude seed
        // is not high-passed by anything.
        float proxy_two_kp = 0.2f;
        float proxy_two_ki = 0.02f;

        // Settling window after going live before accelerometer-bias learning
        // opens. OU-III counts magnetometer updates (250); seconds are the
        // same thing here and do not depend on the magnetometer rate.
        float acc_bias_unlock_sec = 20.0f;
        float handoff_acc_bias_std = 0.03f;   // m/s^2
        // Gravity-agreement gate on the proxy attitude before it is trusted
        // as a seed, or as the frame the magnetic reference is learned in.
        float proxy_gravity_align_sin = 0.06f;   // ~3.4 deg
        float proxy_gravity_lpf_sec   = 8.0f;
        float proxy_gravity_hold_sec  = 20.0f;
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;
        mekf_ = Mekf(cfg.gyro_noise_density, cfg.gravity_magnitude);
        mekf_.initialize_identity();
        mekf_.set_Racc_std(cfg.sigma_a);
        mekf_.set_Rmag_std(cfg.sigma_m);
        Racc_nominal_ = cfg.sigma_a;

        tuner_ = ::SeaStateAutoTuner(2.0f, 1.0f);
        wave_period_.reset();
        vertical_complementary_.setGains(cfg.proxy_two_kp, cfg.proxy_two_ki);
        vertical_complementary_.reset();
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
        elapsed_sec_ = 0.0f;
        live_sec_ = 0.0f;
        mag_elapsed_sec_ = 0.0f;
        stage_ = StartupStage::Cold;

        enterCold_();
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void update(float dt, const Vector3f& gyro, const Vector3f& acc, float tempC = 35.0f) {
        if (!(dt > 0.0f) || !std::isfinite(dt) || !gyro.allFinite() || !acc.allFinite()) return;

        // Commit the schedule the previous step smoothed, BEFORE this step's
        // measurement reaches the MEKF. This is what makes the active tuning
        // at step k+1 measurable with respect to data through k; see item 1
        // in the header. Doing it after the update below would put y_k inside
        // the covariance y_k is then weighted against.
        applyPendingTune_();

        elapsed_sec_ += dt;

        // Measurement-only front end. Runs from the first sample, in every
        // stage and under both policies, and never reads the MEKF.
        vertical_complementary_.update(dt, gyro, acc, cfg_.gravity_magnitude);
        updateProxyGravityQuality_(dt, acc);
        const float a_up = vertical_complementary_.verticalAccelUpMs2();
        if (vertical_complementary_.isReady()) wave_period_.update(dt, a_up);
        updateTuner_(dt, a_up);

        if (stage_ != StartupStage::Live) {
            tuner_warm_sec_ += dt;
        } else if (!acc_bias_unlocked_) {
            live_since_sec_ += dt;
            if (live_since_sec_ >= cfg_.acc_bias_unlock_sec) {
                acc_bias_unlocked_ = true;
                mekf_.set_acc_bias_updates_enabled(true);
            }
        }

        if (stage_ == StartupStage::Cold) {
            if (cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy) {
                tryProxyHandoff_();
            } else {
                stagedColdStep_(gyro, acc, dt);
            }
            return;
        }

        mekf_.time_update(gyro, dt);
        mekf_.measurement_update_acc_only(acc, tempC);

        pseudo_elapsed_ += dt;
        if (pseudo_elapsed_ >= pseudo_period_sec_) {
            pseudo_elapsed_ = 0.0f;
            mekf_.applyIntegralZeroPseudoMeas();
        }

        if (stage_ == StartupStage::TunerWarm) {
            live_sec_ += dt;
            if (live_sec_ >= cfg_.online_tune_warmup_sec && wave_period_.isReady()) {
                enterLive_();
            }
        } else {
            // Smooth the process-model targets (tau, Sigma_aw, r_S) every
            // sample, but only mark a commit on the cadence tick. P itself is
            // untouched here: it remains the Riccati covariance of the
            // invariant error, changed only by propagation, measurements,
            // reset coordinate transport, or explicit initialization.
            adaptMekf_(dt);
        }
    }

    /*
        Magnetometer first lock is a gauge choice, not a measurement update.
        The filter therefore applies the yaw change through
        apply_world_yaw_gauge(), which rotates R, every world column, all
        world-referred tuning covariances, and the corresponding tangent blocks
        of P. Directly changing only R would leave the group state and its
        covariance in different world frames.
    */
    void updateMag(const Vector3f& mag_body) {
        if (!cfg_.with_mag || !mag_body.allFinite()) return;
        if (!(mag_body.norm() > 1e-9f)) return;
        if (elapsed_sec_ < cfg_.mag_delay_sec) return;

        // Under MahonyProxy the reference is learned in the PROXY's frame,
        // while the MEKF may not exist yet. Reading the MEKF here is exactly
        // the coupling this policy removes: the warming MEKF's tilt error
        // would be baked into a reference that is then locked forever.
        if (cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy) {
            if (!mag_reference_learned_) {
                if (!vertical_complementary_.isReady()) return;
                if (!proxyGravityTrusted_()) return;
                learnMagReferenceFrom_(vertical_complementary_.tiltQuaternion(), mag_body);
                return;
            }
            if (stage_ == StartupStage::Live) mekf_.measurement_update_mag_only(mag_body);
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
            // z is unchanged by the yaw gauge and the horizontal magnitude is
            // gauge-invariant, so this canonical reference is consistent with
            // the transformed state.
            mekf_.set_magnetic_reference_world(Vector3f(horiz, 0.0f, b_world.z()));
            return;
        }
        mekf_.measurement_update_mag_only(mag_body);
    }

    [[nodiscard]] bool magReferenceLearned() const noexcept { return mag_reference_learned_; }
    [[nodiscard]] StartupInitPolicy startupInitPolicy() const noexcept {
        return cfg_.startup_init_policy;
    }
    // The operating point is trustworthy. Under StagedMekf this coincides with
    // going Live; under MahonyProxy it is reached first, while the MEKF has
    // not been seeded yet.
    [[nodiscard]] bool isTunerReady() const noexcept {
        return wave_period_.isReady() && tuner_warm_sec_ >= cfg_.online_tune_warmup_sec;
    }
    [[nodiscard]] float pseudoUpdatePeriodSec() const noexcept { return pseudo_period_sec_; }
    // True when the MEKF was seeded by the timeout rather than by the front
    // end declaring itself ready. Worth surfacing: it means the operating
    // point at handoff was not converged.
    [[nodiscard]] bool handoffTimedOut() const noexcept { return handoff_timed_out_; }
    [[nodiscard]] float getRSFilterInput() const noexcept { return RS_filter_input_; }
    void setAdaptEverySecs(float s) { if (s >= 0.0f && std::isfinite(s)) adapt_every_secs_ = s; }
    void setTauScaledPseudoCadence(bool on) {
        tau_scaled_pseudo_cadence_ = on;
        applyPseudoCadence_();
    }
    [[nodiscard]] bool tauScaledPseudoCadence() const noexcept { return tau_scaled_pseudo_cadence_; }

    bool setFixedTuning(float tau_s, float sigma_a, float RS) {
        if (!(tau_s > 0.0f) || !(sigma_a >= 0.0f) || !(RS > 0.0f)) return false;
        if (!std::isfinite(tau_s) || !std::isfinite(sigma_a) || !std::isfinite(RS)) return false;
        fixed_tuning_ = true;
        tune_.tau_applied = tau_s;
        tune_.sigma_applied = sigma_a;
        tune_.RS_applied = RS;
        tau_target_ = tau_s; sigma_target_ = sigma_a; RS_target_ = RS;
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
        }
        if (freeze_RS) {
            if (!(RS > 0.0f)) return false;
            freeze_RS_channel_ = true;
            tune_.RS_applied = RS;
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
    void setRSXYFactor(float k)    { if (k > 0.0f && std::isfinite(k)) R_S_xy_factor_ = k; }
    void setAccNoiseFloorSigma(float s) { if (s >= 0.0f && std::isfinite(s)) noise_floor_sigma_ = s; }
    void setAdaptationTimeConstants(float tau_sec) {
        if (tau_sec > 0.0f && std::isfinite(tau_sec)) adapt_tau_sec_ = tau_sec;
    }
    void setRSAdaptMult(float m)     { if (m > 0.0f && std::isfinite(m)) adapt_RS_mult_ = m; }
    void setRSAdaptSlewLog(float d)  { if (d >= 0.0f && std::isfinite(d)) adapt_RS_slew_log_ = d; }
    void setTauBounds(float lo, float hi) { if (lo > 0.0f && hi > lo) { min_tau_ = lo; max_tau_ = hi; } }
    void setRSBounds(float lo, float hi)  { if (lo > 0.0f && hi > lo) { min_RS_ = lo; max_RS_ = hi; } }
    void setMaxSigmaA(float m)       { if (m > 0.0f && std::isfinite(m)) max_sigma_a_ = m; }

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
    [[nodiscard]] bool  wavePeriodReady()  const noexcept { return wave_period_.isReady(); }
    [[nodiscard]] float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }

    [[nodiscard]] Eigen::Quaternionf quaternion() const { return mekf_.quaternion(); }
    [[nodiscard]] Vector3f get_velocity() const { return mekf_.get_velocity(); }
    [[nodiscard]] Vector3f get_position() const { return mekf_.get_position(); }
    [[nodiscard]] Vector3f get_world_accel() const { return mekf_.get_world_accel(); }

private:
    void enterCold_() {
        mekf_.set_linear_block_enabled(false);
        if (cfg_.freeze_acc_bias_until_live) mekf_.set_acc_bias_updates_enabled(false);
        mekf_.set_Racc_std(Vector3f::Constant(cfg_.Racc_warmup_std));
    }

    void enterLive_() {
        stage_ = StartupStage::Live;
        live_since_sec_ = 0.0f;
        mekf_.set_linear_block_enabled(true);
        // Accelerometer bias stays frozen for a settling window after going
        // live. It is only weakly separable from tilt, so letting it learn
        // through the seed transient parks a tilt error in the bias state
        // permanently. OU-III gates this the same way.
        mekf_.set_acc_bias_updates_enabled(!cfg_.freeze_acc_bias_until_live);
        mekf_.set_Racc_std(Racc_nominal_);
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void updateTuner_(float dt, float a_up) {
        if (!enable_tuner_) return;

        const float f_hint = wave_period_.isReady() ? wave_period_.getFrequencyHz() : 0.2f;
        tuner_.update(dt, a_up, f_hint);
        if (fixed_tuning_) return;

        const float f_tune = std::max(kMinTuneFreqHz,
                                      wave_period_.isReady() ? wave_period_.getFrequencyHz()
                                                             : kMinTuneFreqHz);
        const float var = tuner_.getAccelVariance();
        const float wave_var = std::max(0.0f, var - noise_floor_sigma_ * noise_floor_sigma_);
        const float sigma_wave = std::sqrt(wave_var);

        if (!freeze_ou_channel_) {
            tau_target_ = std::min(max_tau_, std::max(min_tau_, tau_coeff_ * 0.5f / f_tune));
            sigma_target_ = std::min(max_sigma_a_, sigma_coeff_ * sigma_wave);
        }
        if (!freeze_RS_channel_) {
            const float tau_for_RS = freeze_ou_channel_ ? tune_.tau_applied : tau_target_;
            const float sig_for_RS = freeze_ou_channel_ ? tune_.sigma_applied : sigma_target_;
            const float raw = R_S_coeff_ * sig_for_RS * tau_for_RS * tau_for_RS * tau_for_RS;
            RS_target_ = std::min(max_RS_, std::max(min_RS_, raw));
        }
    }

    /*
        Legacy staged path: the MEKF levels itself while degraded. Kept so the
        two policies can be compared on identical input.
    */
    void stagedColdStep_(const Vector3f& gyro, const Vector3f& acc, float dt) {
        const bool levelled = ::seastate::common::runStartupGravityInit(
            gyro, acc, dt, elapsed_sec_,
            cfg_.gravity_magnitude,
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

    /*
        MahonyProxy handoff. Wait until the front end has actually converged --
        the proxy has a settled attitude, the wave-period channel is ready, the
        warmup has elapsed, and (with a magnetometer) north has been gauged --
        then seed the MEKF once and go straight to Live. The degraded TunerWarm
        configuration is never entered.
    */
    void tryProxyHandoff_() {
        if (elapsed_sec_ < cfg_.proxy_startup_min_sec) return;
        if (!vertical_complementary_.isInitialized()) return;

        const bool timed_out = elapsed_sec_ >= cfg_.proxy_startup_timeout_sec;
        if (!timed_out) {
            if (!vertical_complementary_.isReady()) return;
            // The quality gate, not a timer. OU-III's measured proxy tilt error
            // is 2.76 deg on the worst record at 7-20 s and only 0.85 deg by
            // 40 s, so handing off on elapsed time alone seeds an attitude that
            // is still converging. A 2 deg seed error is 0.34 m/s^2 of apparent
            // specific force -- seven times the true accelerometer bias -- and
            // the bias state, being only weakly separable from tilt, absorbs it
            // and never gives it back.
            if (!proxyGravityTrusted_()) return;
            if (!isTunerReady()) return;
            if (cfg_.with_mag && !mag_reference_learned_) return;
        }

        // Proxy tilt carries no meaningful yaw, so compose the anchored
        // heading learned at magnetic lock onto it. Without a magnetometer
        // there is nothing to anchor to and yaw stays at the arbitrary zero.
        const Eigen::Quaternionf q_tilt = vertical_complementary_.tiltQuaternion();
        const Eigen::Quaternionf q_bw =
            cfg_.with_mag
                ? Eigen::Quaternionf(Eigen::AngleAxisf(mag_yaw_anchor_rad_,
                                                       Eigen::Vector3f::UnitZ()) * q_tilt)
                : q_tilt;

        mekf_.initialize_from_truth(q_bw.normalized(),
                                    Vector3f::Zero(), Vector3f::Zero(),
                                    Vector3f::Zero(), Vector3f::Zero(),
                                    Vector3f::Zero(), Vector3f::Zero());
        seedHandoffAttitudeCovariance_();
        // The seed knows nothing about the kinematic chain. Leaving the
        // constructor's tiny prior here asserts v, p and S are known to a
        // centimetre at a random phase of a wave, and the resulting transient
        // is absorbed by the accelerometer bias -- which is observationally
        // confounded with tilt and does not come back.
        mekf_.set_initial_linear_uncertainty(1.0f, 2.0f, 5.0f, 1.0f);
        // Same mistake one slot over: the seed sets b_a = 0, and the
        // constructor's prior claims that to 0.01 m/s^2 when the real bias
        // runs several times larger. An over-tight bias prior is not
        // conservative -- it makes the filter refuse the correction it needs.
        mekf_.set_initial_acc_bias_std(cfg_.handoff_acc_bias_std);
        if (mag_reference_learned_) {
            mekf_.set_magnetic_reference_world(B_w_learned_);
        }
        handoff_timed_out_ = timed_out;
        enterLive_();
    }

    /*
        Capture the world magnetic reference in a frame whose yaw is anchored
        to the field rather than to whatever heading the observer happened to
        start at. Capturing B_w = R * m directly would bake that arbitrary yaw
        into the reference and leave a constant heading offset that every other
        channel hides.
    */
    void learnMagReferenceFrom_(const Eigen::Quaternionf& q_bw, const Vector3f& mag_body) {
        const Vector3f b = q_bw.toRotationMatrix() * mag_body;
        const float horiz = std::hypot(b.x(), b.y());
        if (!(horiz > 1e-9f) || !std::isfinite(b.z())) return;
        mag_yaw_anchor_rad_ = -std::atan2(b.y(), b.x());
        B_w_learned_ = Vector3f(horiz, 0.0f, b.z());
        mag_reference_learned_ = true;
    }

    /*
        Seed the attitude block with the proxy's own uncertainty. Tilt is
        observable from gravity and is the tighter of the two; yaw is only as
        good as the magnetic gauge, and is left wide open when there is no
        magnetometer to gauge with.
    */
    /*
        Sustained agreement between the proxy's predicted gravity direction and
        the measured specific force. Instantaneous agreement means nothing in a
        seaway -- orbital acceleration alone can line them up for an instant --
        so the residual is low-passed and then required to stay under
        threshold for a continuous hold.
    */
    void updateProxyGravityQuality_(float dt, const Vector3f& acc) {
        if (!vertical_complementary_.isInitialized()) { proxy_gravity_good_sec_ = 0.0f; return; }
        const float sin_res = ::seastate::common::gravityAlignResidualSin(
            vertical_complementary_.quaternion(), acc);
        const float a = 1.0f - std::exp(-dt / std::max(1e-3f, cfg_.proxy_gravity_lpf_sec));
        proxy_gravity_res_lpf_ += a * (sin_res - proxy_gravity_res_lpf_);
        if (proxy_gravity_res_lpf_ <= cfg_.proxy_gravity_align_sin) {
            proxy_gravity_good_sec_ += dt;
        } else {
            proxy_gravity_good_sec_ = 0.0f;
        }
    }
    [[nodiscard]] bool proxyGravityTrusted_() const {
        return proxy_gravity_good_sec_ >= cfg_.proxy_gravity_hold_sec;
    }

    void seedHandoffAttitudeCovariance_() {
        auto& P = mekf_.covariance_full();
        const float st = std::max(1e-6f, cfg_.proxy_handoff_tilt_sigma_rad);
        const float sy = (cfg_.with_mag && mag_reference_learned_)
                             ? std::max(1e-6f, cfg_.proxy_handoff_yaw_sigma_rad)
                             : 3.14159265f;  // heading is unknown without a gauge
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

    void adaptMekf_(float dt) {
        if (!enable_tuner_ || fixed_tuning_) return;

        // The smoother sees every sample; only the commit is cadenced.
        if (!freeze_ou_channel_) {
            const float a = 1.0f - std::exp(-dt / adapt_tau_sec_);
            tune_.tau_applied   += a * (tau_target_ - tune_.tau_applied);
            tune_.sigma_applied += a * (sigma_target_ - tune_.sigma_applied);
        }
        if (!freeze_RS_channel_) {
            const float horizon = ::seastate::common::adaptiveSmoothingHorizonSec(
                adapt_RS_mult_, tune_.tau_applied, RS_target_, tune_.RS_applied,
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

    // Deferred commit. See item 1 in the header: this must not run between a
    // measurement arriving and that measurement being used.
    void applyPendingTune_() {
        if (!tune_apply_pending_) return;
        tune_apply_pending_ = false;
        commitTune_();
    }

    void applyPseudoCadence_() {
        if (!tau_scaled_pseudo_cadence_) {
            pseudo_period_sec_ = kPseudoPeriodNominalS;
            return;
        }
        const float tau = tune_.tau_applied;
        if (!(std::isfinite(tau) && tau > 0.0f)) return;
        pseudo_period_sec_ = std::min(std::max(kPseudoTauRatio * tau, kPseudoPeriodMinS),
                                      kPseudoPeriodMaxS);
    }

    void commitTune_() {
        mekf_.set_aw_time_constant(tune_.tau_applied);
        // Commit the S=0 cadence with the same applied tau, so T_S/tau stays
        // constant apart from the safety clamps.
        applyPseudoCadence_();

        const float sZ = tune_.sigma_applied;
        mekf_.set_aw_stationary_std(Vector3f(sZ * S_factor_, sZ * S_factor_, sZ));

        // Hold the information rate 1/(r_S^2 T_S) as the cadence moves with
        // tau. Deliberately not re-clamped -- see item 4 in the header.
        const float scale = (tau_scaled_pseudo_cadence_ && pseudo_period_sec_ > 0.0f)
                                ? std::sqrt(kPseudoPeriodNominalS / pseudo_period_sec_)
                                : 1.0f;
        const float rs = tune_.RS_applied * scale;
        RS_filter_input_ = rs;
        mekf_.set_RS_noise(Vector3f(rs * R_S_xy_factor_, rs * R_S_xy_factor_, rs));
    }

    struct Vec3LPF {
        Eigen::Vector3f state = Eigen::Vector3f::Zero();
        bool initialized = false;

        void reset() { state.setZero(); initialized = false; }

        Eigen::Vector3f step(const Eigen::Vector3f& x, float dt, float tau_sec) {
            if (!x.allFinite()) return state;
            const float tau = std::max(1.0e-3f, tau_sec);
            const float alpha = 1.0f - std::exp(-dt / tau);
            if (!initialized) { state = x; initialized = true; return state; }
            state += alpha * (x - state);
            return state;
        }
    };

    static constexpr float kMinTuneFreqHz = 0.03f;

    // Self-similar S=0 cadence, T_S = c_T * tau. The historical operating
    // point was T_S = 15 ms at the initial applied tau = 1.1 s, so this ratio
    // reproduces it exactly and scales from there.
    static constexpr float kPseudoPeriodNominalS = 0.015f;
    static constexpr float kPseudoTauNominalS    = 1.1f;
    static constexpr float kPseudoTauRatio = kPseudoPeriodNominalS / kPseudoTauNominalS;
    // A pseudo update cannot outrun the 200 Hz IMU schedule; the upper guard
    // is inactive over the tau <= 12 s envelope.
    static constexpr float kPseudoPeriodMinS = 0.005f;
    static constexpr float kPseudoPeriodMaxS = 0.25f;

    Config cfg_{};
    Mekf   mekf_{};
    StartupStage stage_{StartupStage::Cold};

    ::SeaStateAutoTuner tuner_{2.0f, 1.0f};
    ::WavePeriodEstimator wave_period_{};
    ::VerticalAccelComplementary vertical_complementary_{};

    ::seastate::common::StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF bootstrap_gravity_slow_lpf_{};
    float bootstrap_gravity_good_sec_ = 0.0f;

    TfgTuneState tune_{};
    float tau_target_   = 1.1f;
    float sigma_target_ = 1e-2f;
    float RS_target_    = 0.5f;

    float tau_coeff_    = 1.0f;
    float sigma_coeff_  = 1.0f;
    float R_S_coeff_    = 0.35f;
    float S_factor_     = 1.87f;
    float R_S_xy_factor_ = 1.0f;
    float noise_floor_sigma_ = 0.12f;

    float adapt_tau_sec_    = 1.8f;
    float adapt_RS_mult_    = 3.0f;
    float adapt_RS_slew_log_ = 0.0f;

    float min_tau_ = 0.02f, max_tau_ = 12.0f;
    // 0.15, matching OU-III's MIN_R_S. The old 0.4 was the binding
    // constraint on every low-motion sea rather than a safety limit; see item
    // 3 in the header.
    float min_RS_  = 0.15f, max_RS_  = 400.0f;
    float max_sigma_a_ = 6.0f;

    bool enable_tuner_ = true;
    bool fixed_tuning_ = false;
    bool freeze_ou_channel_ = false;
    bool freeze_RS_channel_ = false;

    Vector3f Racc_nominal_{Vector3f::Constant(0.5f)};

    float elapsed_sec_ = 0.0f;
    float live_sec_ = 0.0f;
    float mag_elapsed_sec_ = 0.0f;
    float pseudo_elapsed_ = 0.0f;
    float pseudo_period_sec_ = kPseudoPeriodNominalS;

    // Adaptation cadence and the deferred commit that keeps the schedule
    // exogenous with respect to the current measurement.
    float adapt_every_secs_   = 0.1f;
    float adapt_elapsed_sec_  = 0.0f;
    bool  tune_apply_pending_ = false;
    bool  tau_scaled_pseudo_cadence_ = true;
    float RS_filter_input_ = 0.5f;

    // Startup front end. tuner_warm_sec_ counts from the first sample rather
    // than from the MEKF going live, because under MahonyProxy the front end
    // converges before the MEKF exists.
    float tuner_warm_sec_ = 0.0f;
    bool  mag_reference_learned_ = false;
    bool  handoff_timed_out_ = false;
    bool  acc_bias_unlocked_ = false;
    float live_since_sec_ = 0.0f;
    float proxy_gravity_res_lpf_ = 1.0f;
    float proxy_gravity_good_sec_ = 0.0f;
    float mag_yaw_anchor_rad_ = 0.0f;
    Vector3f B_w_learned_{Vector3f::UnitX()};
};

}  // namespace ocean_imu::tfg
