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

    5. THE SIGMA CHANNEL IS A WAVE-BAND QUANTITY. sigma_aw is meant to be the
       acceleration the OU process has to carry, and a broadband variance is
       not that: it also contains the sensor floor, the engine band and any
       drift below the sea. OU-III measures it through a period-scaled
       band-pass whose corners are fixed multiples of the tuning frequency, so
       the transfer shape is fixed in f/f_tune -- the condition the JONSWAP
       similarity argument for sigma_aw needs -- and it refers the pre-band
       noise floor through that same band's white-noise variance gain rather
       than subtracting a constant. This filter used the raw signal and a
       constant floor; both are replaced here.

    STARTUP: THE PROXY OWNS TILT AND MAGNETIC LEARNING. See StartupInitPolicy.

    ------------------------------------------------------------------------
    MAGNETIC ACQUISITION AND HARD IRON, brought to parity with the OU families
    ------------------------------------------------------------------------

    The previous version of this file learned the world magnetic reference
    from ONE sample, in the proxy's tilt frame, the first time the gravity gate
    was happy. That is the whole acquisition: one reading of a noisy vector at
    one phase of one wave, frozen for the rest of the run. Both OU families
    stopped doing that, in three steps, and all three are carried over here.

    A. THE REFERENCE IS AVERAGED, NOT SAMPLED. MagAutoTuner accumulates the
       field in the proxy's yaw-stripped tilt frame over a window measured in
       seconds rather than samples -- in waves the tilt error is periodic, so
       what the window has to buy is whole wave periods, and 128 samples at a
       25 Hz mag ODR is 5.1 s, short enough to lock in the phase it started on
       instead of cancelling it.

    B. THE ACQUISITION HAS TWO STAGES. A provisional reference is locked as
       soon as the gate allows, because a device that reports no heading for
       105 s is not a device. It is then RE-LEARNED once the filter is live and
       the proxy has actually converged; measured on OU-III, proxy tilt error
       over a 23 s window from 7 s runs up to 2.76 deg and from 40 s up to 0.85
       deg, so the first stage is buying availability and the second one is
       buying accuracy. Both stages average in the PROXY's tilt frame, never
       the MEKF's: the MEKF has been steering to the provisional reference the
       refinement exists to replace, so its tilt carries that reference's error
       and the refinement would be self-confirming.

    C. HARD IRON IS TRACKED CONTINUOUSLY. The MEKF has no magnetometer-bias
       state, so a body-fixed offset is heading error one-for-one against the
       horizontal field -- which is most of why this filter's yaw error ran
       about three times OU-III's and OU-II's. ContinuousMagHardIronEstimator
       never closes its accumulation, is driven by the proxy tilt and the raw
       magnetometer (so no loop closes through the MEKF), and the applied
       offset and the magnetic reference move together out of the same
       statistics. The correction is a change of measurement-model parameters
       only: no attitude state is written, and the magnetometer update walks
       yaw to the corrected heading over its own time constant.

    While the provisional reference is the one in force, accelerometer-bias
    learning is held off. The bias is only weakly separable from tilt, so
    letting it fit itself to a reference that is about to be replaced parks the
    provisional reference's error in the bias state permanently.
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
#include "tuner/AdaptiveWaveBandPass.h"
#include "tuner/ContinuousMagHardIronEstimator.h"
#include "tuner/MagAutoTuner.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/WavePeriodEstimator.h"

namespace ocean_imu::tfg {

constexpr float MAG_DELAY_SEC = 7.0f;

// JONSWAP-similar acceleration-variance band, as OU-III. Dimensionless wave-
// band ratios, plus absolute safety clamps.
constexpr float SIGMA_BAND_LOW_RATIO_DEFAULT  = 0.5f;
constexpr float SIGMA_BAND_HIGH_RATIO_DEFAULT = 4.0f;
constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = 0.01f;
constexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = 6.0f;

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
        // MahonyProxy, as in both OU families. StagedMekf stays as the
        // ablation against the previous behaviour; it keeps its single
        // one-sample reference and gets none of the two-stage acquisition.
        StartupInitPolicy startup_init_policy = StartupInitPolicy::MahonyProxy;

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
        // No gauge to anchor to: heading is genuinely unknown, and saying so is
        // what lets the first magnetometer updates move it.
        float proxy_handoff_yaw_sigma_free_rad = 1.5708f;

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

        // Gravity-agreement gate on the proxy attitude before it is trusted as
        // a seed, or as the frame the magnetic reference is learned in.
        //
        // The residual is formed against a low-passed accelerometer rather than
        // being low-passed itself, and the hold is a leaky counter that decays
        // at twice the rate it fills, so a gate that keeps flickering never
        // accumulates credit. Both are OU-III's construction and its numbers.
        float proxy_gravity_align_sin = 0.075f;
        float proxy_gravity_lpf_sec   = 1.0f;
        float proxy_gravity_hold_sec  = 2.0f;
        // Veto only truly violent motion; the gate above does the rest.
        float mag_extreme_gyro_dps    = 30.0f;
        // A platform that never satisfies the gate still has to acquire, or the
        // filter reports no heading at all.
        float mag_tilt_fallback_sec   = 30.0f;
        float mag_init_min_mag_norm   = 1e-3f;

        // Windowed acquisition of the world reference. Held in seconds so it
        // does not silently shorten at a higher magnetometer ODR; see item A.
        int   mag_min_samples    = 128;
        float mag_min_window_sec = 15.0f;
        float mag_max_window_sec = 0.0f;      // 0 = no forced timeout
        float mag_sample_dt_sec  = 1.0f / 200.0f;
        // Hold the provisional stage off. 0 means "as soon as the gravity gate
        // is happy", and is the default because the refinement carries the
        // accuracy now.
        float proxy_mag_settle_sec = 0.0f;

        // Second-stage acquisition; see item B.
        bool  mag_refine_enabled    = true;
        float mag_refine_start_sec  = 90.0f;
        float mag_refine_window_sec = 30.0f;

        // Accel/gyro quality weighting can phase-select wave motion, so it
        // stays off in waves.
        bool  mag_enable_quality_weighting = false;
        float mag_min_effective_weight     = 0.0f;
        float mag_acc_norm_rel_soft        = 0.22f;
        float mag_gyro_soft_dps            = 45.0f;

        // Startup hard-iron estimate from the acquisition window itself. Off by
        // default and for the reason MagAutoTuner gives: over a single window
        // the offset is least identifiable, and a wrong one subtracted
        // everywhere is worse than none. The continuous estimator below is the
        // one that carries this.
        bool  mag_estimate_hard_iron = false;

        // Continuous hard-iron estimation; see item C.
        bool  mag_continuous_hard_iron        = true;
        float mag_hi_memory_sec               = 600.0f;
        float mag_hi_model_ridge              = 4.0e-3f;
        float mag_hi_model_ridge_relative     = 0.5f;
        float mag_hi_min_information          = 2.0f;
        float mag_hi_min_effective_weight     = 500.0f;
        float mag_hi_max_residual_rms_uT      = 3.0f;
        float mag_hi_max_bias_fraction        = 0.35f;
        // Fraction of the fit the filter is willing to apply, and the time
        // constant it moves over. The ridge already shrinks what the model
        // cannot see; this is the blunter statement that a calibration nobody
        // has checked should not arrive as a step.
        float mag_hi_apply_fraction           = 1.0f;
        float mag_hi_slew_tau_sec             = 45.0f;

        // Period-scaled sigma band; see item 5.
        bool  wave_band_tuning      = true;
        float sigma_band_low_ratio  = SIGMA_BAND_LOW_RATIO_DEFAULT;
        float sigma_band_high_ratio = SIGMA_BAND_HIGH_RATIO_DEFAULT;
        float sigma_band_min_hz     = SIGMA_BAND_MIN_HZ_DEFAULT;
        float sigma_band_max_hz     = SIGMA_BAND_MAX_HZ_DEFAULT;
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
        sigma_wave_band_.setRatios(cfg.sigma_band_low_ratio, cfg.sigma_band_high_ratio);
        sigma_wave_band_.setLimitsHz(cfg.sigma_band_min_hz, cfg.sigma_band_max_hz);
        sigma_wave_band_.reset();
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
        elapsed_sec_ = 0.0f;
        live_sec_ = 0.0f;
        mag_elapsed_sec_ = 0.0f;
        stage_ = StartupStage::Cold;

        beginMagAcquisition_();

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
        updateProxyGravityQuality_(dt, gyro, acc);
        last_acc_body_ = acc;
        last_gyro_body_ = gyro;
        have_last_imu_ = true;
        const float a_up = vertical_complementary_.verticalAccelUpMs2();
        if (vertical_complementary_.isReady()) wave_period_.update(dt, a_up);
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
                stagedColdStep_(gyro, acc, dt);
            }
            return;
        }

        // Staged here so the inflation lands inside the propagation below,
        // never between a measurement arriving and that measurement being used.
        periodicAwCovSyncTick_(dt);

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
        Under StagedMekf the magnetometer first lock is a gauge choice, not a
        measurement update, and the filter applies the yaw change through
        apply_world_yaw_gauge(), which rotates R, every world column, all
        world-referred tuning covariances, and the corresponding tangent blocks
        of P. Directly changing only R would leave the group state and its
        covariance in different world frames.

        Under MahonyProxy nothing is written to the MEKF at lock time: the
        reference and the yaw gauge are learned in the proxy's frame while the
        MEKF does not exist yet, and both are carried into the seed. See
        learnMagReferenceWindowed_ and tryProxyHandoff_.
    */
    void updateMag(const Vector3f& mag_body) {
        if (!cfg_.with_mag || !mag_body.allFinite()) return;
        if (!(mag_body.norm() > 1e-9f)) return;
        if (elapsed_sec_ < cfg_.mag_delay_sec) return;

        // Ahead of every gate below, and deliberately. The continuous estimator
        // wants the whole magnetometer record, not the part the startup
        // machinery was willing to average, and it reads a frame the startup
        // machinery does not own.
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
            // z is unchanged by the yaw gauge and the horizontal magnitude is
            // gauge-invariant, so this canonical reference is consistent with
            // the transformed state.
            mekf_.set_magnetic_reference_world(Vector3f(horiz, 0.0f, b_world.z()));
            return;
        }
        mekf_.measurement_update_mag_only(mag_body);
    }

    [[nodiscard]] bool magReferenceLearned() const noexcept { return mag_reference_learned_; }
    // True once the second-stage acquisition has replaced the provisional
    // reference; see maybeRefineMagReference_().
    [[nodiscard]] bool magReferenceRefined() const noexcept { return mag_refine_done_; }
    [[nodiscard]] float magRefineTimeSec() const noexcept { return mag_refine_time_sec_; }
    [[nodiscard]] float magNorthLockTimeSec() const noexcept { return mag_north_lock_time_sec_; }
    // Body-frame hard-iron offset removed from the magnetometer stream.
    [[nodiscard]] const Vector3f& magHardIronBodyUT() const noexcept {
        return mag_hard_iron_body_uT_;
    }
    // The part of the above the continuous estimator is responsible for.
    [[nodiscard]] const Vector3f& magContinuousHardIronAppliedUT() const noexcept {
        return mag_hi_applied_body_uT_;
    }
    [[nodiscard]] const ContinuousMagHardIronEstimator& magContinuousHardIron() const noexcept {
        return mag_hi_estimator_;
    }
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
    // Ablation hooks for the two schedule features brought over from OU-III.
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
    [[nodiscard]] bool  wavePeriodReady()  const noexcept { return wave_period_.isReady(); }
    [[nodiscard]] float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }

    [[nodiscard]] Eigen::Quaternionf quaternion() const { return mekf_.quaternion(); }
    [[nodiscard]] Vector3f get_velocity() const { return mekf_.get_velocity(); }
    [[nodiscard]] Vector3f get_position() const { return mekf_.get_position(); }
    [[nodiscard]] Vector3f get_world_accel() const { return mekf_.get_world_accel(); }

private:
    // Magnetic acquisition and hard-iron tracking, reset together because the
    // applied offset and the reference it belongs to must never be carried over
    // from a previous configuration independently of each other.
    void beginMagAcquisition_() {
        ::MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min             = cfg_.mag_init_min_mag_norm;
        mag_cfg.min_samples              = cfg_.mag_min_samples;
        mag_cfg.min_window_sec           = cfg_.mag_min_window_sec;
        mag_cfg.max_window_sec           = cfg_.mag_max_window_sec;
        mag_cfg.sample_dt_sec            = cfg_.mag_sample_dt_sec;
        mag_cfg.gravity_ref              = cfg_.gravity_magnitude;
        mag_cfg.enable_quality_weighting = cfg_.mag_enable_quality_weighting;
        mag_cfg.estimate_hard_iron       = cfg_.mag_estimate_hard_iron;
        mag_cfg.min_effective_weight     = cfg_.mag_min_effective_weight;
        mag_cfg.acc_norm_rel_soft        = cfg_.mag_acc_norm_rel_soft;
        mag_cfg.gyro_soft_dps            = cfg_.mag_gyro_soft_dps;
        mag_auto_tuner_.setConfig(mag_cfg);

        ContinuousMagHardIronEstimator::Config hi_cfg;
        hi_cfg.memory_sec           = cfg_.mag_hi_memory_sec;
        hi_cfg.model_ridge          = cfg_.mag_hi_model_ridge;
        hi_cfg.model_ridge_relative = cfg_.mag_hi_model_ridge_relative;
        hi_cfg.min_information      = cfg_.mag_hi_min_information;
        hi_cfg.min_effective_weight = cfg_.mag_hi_min_effective_weight;
        hi_cfg.max_residual_rms_uT  = cfg_.mag_hi_max_residual_rms_uT;
        hi_cfg.max_bias_fraction    = cfg_.mag_hi_max_bias_fraction;
        hi_cfg.min_mag_norm_uT      = cfg_.mag_init_min_mag_norm;
        mag_hi_estimator_.setConfig(hi_cfg);

        mag_reference_learned_ = false;
        mag_world_ref_valid_   = false;
        mag_world_ref_uT_.setZero();
        mag_yaw_anchor_rad_ = 0.0f;
        have_mag_yaw_anchor_ = false;

        mag_hard_iron_body_uT_.setZero();
        mag_hi_startup_body_uT_.setZero();
        mag_hi_applied_body_uT_.setZero();
        mag_hi_anchor_bias_body_uT_.setZero();
        mag_hi_anchor_world_ref_uT_.setZero();
        mag_hi_anchored_ = false;

        mag_refine_started_ = false;
        mag_refine_done_    = false;
        mag_refine_time_sec_ = NAN;
        mag_north_lock_time_sec_ = NAN;

        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;
        last_hi_sample_t_  = NAN;
        last_hi_apply_t_   = NAN;

        gravity_gate_acc_lpf_.reset();
        proxy_gravity_good_sec_ = 0.0f;

        last_acc_body_.setZero();
        last_gyro_body_.setZero();
        have_last_imu_ = false;

        acc_bias_unlocked_ = false;
        live_since_sec_ = 0.0f;
        // Under the proxy policy the provisional reference is deliberately
        // cheap and early, so the accelerometer bias must not fit itself to it.
        acc_bias_hold_ = usingProxyInit_() && cfg_.with_mag && cfg_.mag_refine_enabled;
    }

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
        // permanently. OU-III gates this the same way. Under the proxy policy
        // it stays shut for longer still, until the magnetic reference the
        // filter is steering to is the refined one; see acc_bias_hold_.
        acc_bias_unlocked_ = false;
        mekf_.set_acc_bias_updates_enabled(false);
        maybeUnlockAccBias_();
        mekf_.set_Racc_std(Racc_nominal_);
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    // Two conditions, both about the same confound: the accelerometer bias is
    // only weakly separable from a tilt error, so it must not open while either
    // the seed transient or a provisional magnetic reference is still in force.
    void maybeUnlockAccBias_() {
        if (acc_bias_unlocked_) return;
        if (cfg_.freeze_acc_bias_until_live) {
            if (live_since_sec_ < cfg_.acc_bias_unlock_sec) return;
            if (acc_bias_hold_ && !mag_refine_done_) return;
        }
        acc_bias_unlocked_ = true;
        mekf_.set_acc_bias_updates_enabled(true);
    }

    // Pre-band noise floor referred to the band the variance is measured in.
    // Subtracting the raw bench floor from a band-passed variance would
    // over-subtract by whatever the band rejects, and the band moves.
    [[nodiscard]] float bandNoiseFloorSigma_() const noexcept {
        if (!cfg_.wave_band_tuning || !sigma_wave_band_.isReady()) {
            return noise_floor_sigma_;
        }
        const float gain = sigma_wave_band_.whiteNoiseVarianceGain();
        if (!(std::isfinite(gain) && gain >= 0.0f)) return noise_floor_sigma_;
        return noise_floor_sigma_ * std::sqrt(gain);
    }

    void updateTuner_(float dt, float a_up) {
        if (!enable_tuner_) return;

        const float f_hint = wave_period_.isReady() ? wave_period_.getFrequencyHz() : 0.2f;

        // Band corners follow the previous smoothed tuner frequency, so the
        // band motion stays one-sample predictable while remaining
        // measurement-only.
        float f_band = tuner_.isFreqReady() ? tuner_.getFrequencyHz() : f_hint;
        if (!std::isfinite(f_band) || f_band < kMinTuneFreqHz) f_band = kMinTuneFreqHz;
        f_band = std::min(f_band, kMaxFreqHz);

        const float a_for_variance = cfg_.wave_band_tuning
            ? sigma_wave_band_.step(a_up, dt, f_band)
            : a_up;

        tuner_.update(dt, a_for_variance, f_hint);
        if (fixed_tuning_) return;

        float f_tune = tuner_.isFreqReady() ? tuner_.getFrequencyHz() : kMinTuneFreqHz;
        if (!std::isfinite(f_tune) || f_tune < kMinTuneFreqHz) f_tune = kMinTuneFreqHz;
        f_tune = std::min(f_tune, kMaxFreqHz);

        const float band_noise_sigma = bandNoiseFloorSigma_();
        const float var_noise = band_noise_sigma * band_noise_sigma;
        const float var_total = tuner_.isVarReady()
            ? std::max(0.0f, tuner_.getAccelVariance())
            : var_noise;
        const float sigma_wave = std::sqrt(std::max(1e-6f, var_total - var_noise));

        if (!freeze_ou_channel_) {
            tau_target_ = std::min(max_tau_, std::max(min_tau_, tau_coeff_ * 0.5f / f_tune));
            sigma_target_ = std::min(max_sigma_a_, sigma_coeff_ * sigma_wave);
            // Before the variance channel is trusted, the schedule must not
            // read as a dead-calm sea: that is the operating point at which
            // r_S ~ sigma tau^3 collapses.
            if (!tuner_.isVarReady()) {
                sigma_target_ = std::max(sigma_target_,
                                         std::max(0.05f, band_noise_sigma));
            }
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

        // The timeout is held clear of the magnetometer acquisition it would
        // otherwise cut short. Firing while the reference is still averaging
        // hands over with no yaw gauge at all, which is a far worse start than
        // simply waiting: the gauge is the one chance to put the filter on
        // north before it goes live.
        const float mag_acquire_deadline =
            cfg_.with_mag
                ? cfg_.proxy_mag_settle_sec +
                      2.0f * std::max(cfg_.mag_min_window_sec, 1.0f) +
                      cfg_.mag_tilt_fallback_sec
                : 0.0f;
        const float timeout_sec =
            std::max(cfg_.proxy_startup_timeout_sec, mag_acquire_deadline);

        const bool timed_out = elapsed_sec_ >= timeout_sec;
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
        // heading learned at magnetic lock onto it. Without a gauge there is
        // nothing to anchor to and yaw stays at the arbitrary zero, which the
        // seed covariance below then says out loud.
        const Eigen::Quaternionf q_tilt = vertical_complementary_.tiltQuaternion();
        const Eigen::Quaternionf q_bw =
            have_mag_yaw_anchor_
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
        if (mag_world_ref_valid_) {
            mekf_.set_magnetic_reference_world(mag_world_ref_uT_);
        }
        handoff_timed_out_ = timed_out;
        enterLive_();
    }

    // ---------------------------------------------------------------------
    // Magnetic acquisition under the proxy policy
    // ---------------------------------------------------------------------

    [[nodiscard]] bool usingProxyInit_() const noexcept {
        return cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy;
    }

    void proxyUpdateMag_(const Vector3f& mag_body) {
        // Hold the whole magnetometer path off until the observer has settled.
        // This sits ahead of the eligibility clock deliberately, so the
        // tilt-fallback timer cannot start running and then wave the
        // accumulation through on an attitude that is still converging.
        if (!mag_reference_learned_ && elapsed_sec_ < cfg_.proxy_mag_settle_sec) return;
        if (!vertical_complementary_.isInitialized()) return;

        if (!std::isfinite(mag_init_eligible_t0_)) mag_init_eligible_t0_ = elapsed_sec_;

        if (!mag_reference_learned_) {
            const bool fallback_ok =
                (elapsed_sec_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec;
            if (!proxyGravityTrusted_() && !fallback_ok) return;
            if (!have_last_imu_) return;
            learnMagReferenceWindowed_(mag_body);
        }

        maybeRefineMagReference_(mag_body);
        maybeApplyContinuousHardIron_();

        // Magnetometer corrections reach the MEKF only once it owns the
        // attitude. Before the handoff its state is not the one being solved.
        if (mag_reference_learned_ && stage_ == StartupStage::Live) {
            mekf_.measurement_update_mag_only(mag_body - mag_hard_iron_body_uT_);
        }
    }

    /*
        First-stage acquisition: average the field in the proxy's yaw-stripped
        tilt frame until the window closes, then take the reference and the
        heading gauge from that average.

        Stripping yaw makes the frame invariant to the observer's arbitrary
        startup heading, so no heading leaks into the learned reference. The
        gauge is the heading of averaged magnetic north in that frame; driving
        the boat's absolute yaw to its negative puts learned north on +X, which
        is what makes the canonical (h, 0, z) reference consistent with the
        world frame the MEKF is about to be seeded in.
    */
    void learnMagReferenceWindowed_(const Vector3f& mag_body) {
        const float dt_mag =
            (std::isfinite(last_mag_sample_t_) && elapsed_sec_ > last_mag_sample_t_)
                ? (elapsed_sec_ - last_mag_sample_t_)
                : cfg_.mag_sample_dt_sec;
        last_mag_sample_t_ = elapsed_sec_;

        if (!mag_auto_tuner_.addSampleWithTiltQuatDt(
                dt_mag, vertical_complementary_.tiltQuaternion(),
                last_acc_body_, last_gyro_body_, mag_body)) {
            return;
        }

        Vector3f ref;
        if (!mag_auto_tuner_.getMagWorldRef(ref) || !ref.allFinite() ||
            !(ref.norm() > cfg_.mag_init_min_mag_norm)) {
            return;
        }
        setMagWorldRef_(ref);

        const float gauge = mag_auto_tuner_.getYawGaugeCorrectionRad();
        if (std::isfinite(gauge)) {
            mag_yaw_anchor_rad_ = wrapPi_(-gauge);
            have_mag_yaw_anchor_ = true;
        }

        Vector3f hard_iron;
        mag_hi_startup_body_uT_ = mag_auto_tuner_.getHardIronBodyUT(hard_iron)
                                      ? hard_iron
                                      : Vector3f::Zero();
        mag_hard_iron_body_uT_ = mag_hi_startup_body_uT_ + mag_hi_applied_body_uT_;

        mag_reference_learned_ = true;
        mag_north_lock_time_sec_ = elapsed_sec_;
    }

    /*
        Second-stage acquisition. The provisional reference was averaged in a
        tilt frame the observer had barely converged, and it is what the filter
        has been steering to ever since. Re-run the same acquisition once the
        MEKF is live and the observer has settled.

        Still in the PROXY's tilt frame, and here the reason has teeth: the MEKF
        has been steering to the reference this pass exists to replace, so its
        tilt carries that reference's error and averaging in it would re-derive
        the error it was meant to remove. The observer never saw the reference.

        Both the reference vector and the heading gauge are replaced. The yaw
        write is a step, deliberately: it is the coarse-to-fine alignment, it
        happens once, and it lands long before the scored window opens.
    */
    void maybeRefineMagReference_(const Vector3f& mag_body) {
        if (!cfg_.mag_refine_enabled || mag_refine_done_) return;
        if (!mag_reference_learned_ || stage_ != StartupStage::Live) return;
        if (elapsed_sec_ < cfg_.mag_refine_start_sec) return;
        if (!have_last_imu_) return;

        if (!mag_refine_started_) {
            ::MagAutoTuner::Config refine_cfg = mag_auto_tuner_.config();
            refine_cfg.min_window_sec = cfg_.mag_refine_window_sec;
            refine_cfg.min_samples    = cfg_.mag_min_samples;
            mag_auto_tuner_.setConfig(refine_cfg);
            mag_auto_tuner_.reset();
            mag_refine_started_ = true;
            last_mag_sample_t_  = NAN;
        }

        const float dt_mag =
            (std::isfinite(last_mag_sample_t_) && elapsed_sec_ > last_mag_sample_t_)
                ? (elapsed_sec_ - last_mag_sample_t_)
                : cfg_.mag_sample_dt_sec;
        last_mag_sample_t_ = elapsed_sec_;

        // The same corrected stream the MEKF sees, so an offset already removed
        // is not re-learned into the new reference.
        const Vector3f mag_corrected = mag_body - mag_hard_iron_body_uT_;

        if (!mag_auto_tuner_.addSampleWithTiltQuatDt(
                dt_mag, vertical_complementary_.tiltQuaternion(),
                last_acc_body_, last_gyro_body_, mag_corrected)) {
            return;
        }

        Vector3f ref;
        if (!mag_auto_tuner_.getMagWorldRef(ref) || !ref.allFinite() ||
            !(ref.norm() > cfg_.mag_init_min_mag_norm)) {
            return;
        }
        setMagWorldRef_(ref);

        const float gauge = mag_auto_tuner_.getYawGaugeCorrectionRad();
        if (std::isfinite(gauge)) {
            mekf_.set_attitude_yaw_absolute(wrapPi_(-gauge));
        }

        mag_refine_done_     = true;
        mag_refine_time_sec_ = elapsed_sec_;
        // The reference the bias would have been fitting is now the good one.
        maybeUnlockAccBias_();
    }

    // Every write of the world reference goes through here, so the orchestrator
    // always knows the vector the MEKF is steering to. The MEKF does not offer
    // it back, and the continuous correction moves it by a delta rather than
    // replacing it.
    void setMagWorldRef_(const Vector3f& ref) {
        mag_world_ref_uT_ = ref;
        mag_world_ref_valid_ = true;
        if (stage_ == StartupStage::Live || !usingProxyInit_()) {
            mekf_.set_magnetic_reference_world(ref);
        }
    }

    /*
        Feed the exogenous accumulation. Raw magnetometer -- not the corrected
        stream -- because the estimator is fitting the offset itself and must
        not be shown data with its own answer already subtracted.
    */
    void accumulateContinuousHardIron_(const Vector3f& mag_body) {
        if (!cfg_.mag_continuous_hard_iron || !usingProxyInit_()) return;
        if (!vertical_complementary_.isInitialized()) return;

        const float dt_mag =
            (std::isfinite(last_hi_sample_t_) && elapsed_sec_ > last_hi_sample_t_)
                ? (elapsed_sec_ - last_hi_sample_t_)
                : cfg_.mag_sample_dt_sec;
        last_hi_sample_t_ = elapsed_sec_;

        mag_hi_estimator_.update(dt_mag, vertical_complementary_.tiltQuaternion(), mag_body);
    }

    /*
        Move the applied offset toward the fit, and re-gauge the reference.

        The reference is rebuilt in MagAutoTuner's canonical form -- horizontal
        magnitude on +X, vertical below it -- and NOT merely shifted by the same
        amount as the measurement. A shift that tracks the offset exactly is a
        no-op: subtracting b from every sample and subtracting the matching
        mean(R) b from the reference leaves the innovation identical at the
        attitude the filter already holds, so nothing moves and the standing yaw
        error survives the correction meant to remove it.

        The standing error is a GAUGE: the startup acquisition put the world
        frame's north along the average of the uncorrected field, which is
        magnetic north rotated by whatever the offset contributes. Leaving the
        canonical reference in place while the offset comes out of the stream
        asks the filter for the heading the corrected field implies, and the
        magnetometer update walks yaw there over its own time constant. No
        attitude state is written.

        Only the horizontal magnitude and the vertical component move, and only
        by the amount the offset changes them; they are not recomputed from the
        estimator's own window, which is longer and less selective than the one
        the startup acquisition gated.
    */
    void maybeApplyContinuousHardIron_() {
        if (!cfg_.mag_continuous_hard_iron || !usingProxyInit_()) return;
        if (!mag_reference_learned_ || stage_ != StartupStage::Live) return;
        if (cfg_.mag_refine_enabled && !mag_refine_done_) return;
        if (!mag_world_ref_valid_) return;

        const auto& est = mag_hi_estimator_.estimate();
        if (!est.valid) return;

        if (!mag_hi_anchored_) {
            mag_hi_anchor_bias_body_uT_ = mag_hi_applied_body_uT_;
            mag_hi_anchor_world_ref_uT_ = mag_world_ref_uT_;
            mag_hi_anchored_ = true;
        }

        const Vector3f target = cfg_.mag_hi_apply_fraction * est.bias_body_uT;
        if (!target.allFinite()) return;

        const float dt_apply =
            (std::isfinite(last_hi_apply_t_) && elapsed_sec_ > last_hi_apply_t_)
                ? (elapsed_sec_ - last_hi_apply_t_)
                : cfg_.mag_sample_dt_sec;
        last_hi_apply_t_ = elapsed_sec_;

        const float tau = cfg_.mag_hi_slew_tau_sec;
        const float alpha = (std::isfinite(tau) && tau > 1.0e-3f)
                                ? (1.0f - std::exp(-dt_apply / tau))
                                : 1.0f;

        const Vector3f applied =
            mag_hi_applied_body_uT_ + alpha * (target - mag_hi_applied_body_uT_);
        if (!applied.allFinite()) return;

        // Both evaluated against the statistics as they stand now, so the
        // difference is the offset's doing and nothing else.
        Vector3f level_new, level_anchor;
        if (!mag_hi_estimator_.levelReferenceForBias(applied, level_new) ||
            !mag_hi_estimator_.levelReferenceForBias(mag_hi_anchor_bias_body_uT_,
                                                     level_anchor)) {
            return;
        }

        const float h_new = level_new.head<2>().norm();
        const float h_anchor = level_anchor.head<2>().norm();
        if (!std::isfinite(h_new) || !std::isfinite(h_anchor)) return;

        const float h = mag_hi_anchor_world_ref_uT_.x() + (h_new - h_anchor);
        const float z = mag_hi_anchor_world_ref_uT_.z() +
                        (level_new.z() - level_anchor.z());
        if (!(h > cfg_.mag_init_min_mag_norm) || !std::isfinite(h) || !std::isfinite(z)) {
            return;
        }

        mag_hi_applied_body_uT_ = applied;
        mag_hard_iron_body_uT_ = mag_hi_startup_body_uT_ + applied;
        setMagWorldRef_(Vector3f(h, 0.0f, z));
    }

    static float wrapPi_(float a) {
        constexpr float PI_F = 3.14159265358979323846f;
        if (!std::isfinite(a)) return NAN;
        while (a > PI_F) a -= 2.0f * PI_F;
        while (a <= -PI_F) a += 2.0f * PI_F;
        return a;
    }

    /*
        Sustained agreement between the proxy's predicted gravity direction and
        the measured specific force. Instantaneous agreement means nothing in a
        seaway -- orbital acceleration alone can line them up for an instant --
        so the residual is formed against a low-passed accelerometer and the
        good time is a leaky counter that decays at twice the rate it fills,
        which is what stops a gate that keeps flickering from accumulating
        credit. That is OU-III's construction, replacing this filter's
        low-passed-residual-plus-hard-reset version.
    */
    void updateProxyGravityQuality_(float dt, const Vector3f& gyro, const Vector3f& acc) {
        if (!vertical_complementary_.isInitialized()) { proxy_gravity_good_sec_ = 0.0f; return; }

        const Vector3f acc_lp =
            gravity_gate_acc_lpf_.step(acc, dt, cfg_.proxy_gravity_lpf_sec);
        const float sin_res = ::seastate::common::gravityAlignResidualSin(
            vertical_complementary_.quaternion(), acc_lp);

        const float gyro_dps = gyro.norm() * 57.295779513f;
        const bool extreme_motion =
            !std::isfinite(gyro_dps) || (gyro_dps > cfg_.mag_extreme_gyro_dps);

        const bool good_now = std::isfinite(sin_res) &&
                              (sin_res <= cfg_.proxy_gravity_align_sin) &&
                              !extreme_motion;

        if (good_now) {
            proxy_gravity_good_sec_ = std::min(10.0f, proxy_gravity_good_sec_ + dt);
        } else {
            proxy_gravity_good_sec_ = std::max(0.0f, proxy_gravity_good_sec_ - 2.0f * dt);
        }
    }
    [[nodiscard]] bool proxyGravityTrusted_() const {
        return proxy_gravity_good_sec_ >= cfg_.proxy_gravity_hold_sec;
    }

    /*
        Seed the attitude block with the proxy's own uncertainty. Tilt is
        observable from gravity and is the tighter of the two; yaw is only as
        good as the magnetic gauge, and is left wide open when there is no gauge
        to anchor it.
    */
    void seedHandoffAttitudeCovariance_() {
        auto& P = mekf_.covariance_full();
        const float st = std::max(1e-6f, cfg_.proxy_handoff_tilt_sigma_rad);
        const float sy = have_mag_yaw_anchor_
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

    /*
        Re-align the posterior a_w marginal with the stationary prior on the
        adaptation cadence.

        commitTune_() writes the new stationary scale into the OU process
        covariance, which is how the prior reaches the *increments*. It does not
        touch the marginal the filter is already holding, so after a sea state
        change the posterior a_w variance can sit an order of magnitude under
        the scale the schedule now claims and the accelerometer update is
        weighted against a confidence the filter no longer has any basis for.
        The sync only ever adds the PSD part of the shortfall, so nothing the
        measurements established is thrown away.

        It runs on its own clock rather than inside adaptMekf_(), so the fixed
        and frozen ablation arms get the identical policy and stay matched
        controls rather than differing in two things at once.
    */
    void periodicAwCovSyncTick_(float dt) {
        if (!periodic_aw_cov_sync_ || stage_ != StartupStage::Live) return;
        aw_sync_elapsed_sec_ += dt;
        if (aw_sync_elapsed_sec_ < adapt_every_secs_) return;
        aw_sync_elapsed_sec_ = 0.0f;
        mekf_.synchronize_aw_covariance_to_stationary();
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

    // Floor for the wave-band tuning frequency: 0.03 Hz admits a 33 s swell.
    static constexpr float kMinTuneFreqHz = 0.03f;
    // Absolute ceiling; the sigma band's own Nyquist guard sits below it.
    static constexpr float kMaxFreqHz = 6.0f;

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

    ::AdaptiveWaveBandPass sigma_wave_band_{SIGMA_BAND_LOW_RATIO_DEFAULT,
                                            SIGMA_BAND_HIGH_RATIO_DEFAULT,
                                            SIGMA_BAND_MIN_HZ_DEFAULT,
                                            SIGMA_BAND_MAX_HZ_DEFAULT};
    ::MagAutoTuner mag_auto_tuner_{};
    ContinuousMagHardIronEstimator mag_hi_estimator_{};

    ::seastate::common::StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF bootstrap_gravity_slow_lpf_{};
    float bootstrap_gravity_good_sec_ = 0.0f;

    TfgTuneState tune_{};
    float tau_target_   = 1.1f;
    float sigma_target_ = 1e-2f;
    float RS_target_    = 0.5f;

    /*
        Retuned against the eight JONSWAP and PM-Stokes records after the
        schedule, startup and hard-iron work above, because every one of those
        changes moves what the coefficients are multiplying. Sweeps are in
        docs/tfg-adaptation-parity.md; the shape of each is stated here.

        tau_coeff = 1.0     unchanged. A clean minimum: 0.92 costs +2.4 % of
                            vertical RMS and 1.10 costs +2.5 %.
        sigma_coeff = 1.0   unchanged. OU-III runs 0.9, and 0.9 does buy 0.4 %
                            of vertical error here -- at 4 % on horizontal 3D
                            and 4 % on accelerometer bias. Not a trade worth
                            taking, and the two filters measure sigma through
                            the same band now, so the difference is genuine
                            rather than a units mismatch.
        R_S_coeff = 0.28    was 0.35. The band-passed sigma channel reads lower
                            than the broadband one it replaced, and r_S ~ sigma
                            tau^3 inherits that; 0.35 now over-regularizes the
                            S = 0 constraint. The vertical minimum is flat over
                            0.26..0.30 and 0.28 is where 3D stops improving.
        S_factor = 1.20     was 1.87. Horizontal a_w no longer needs the extra
                            headroom now that the horizontal channel is not
                            absorbing a standing heading error: dropping it
                            takes 5 % off accelerometer-bias error and 0.05 deg
                            off pitch for no vertical cost.
        R_S_xy_factor = 1.15 was 1.0. Mildly loosening the horizontal S = 0
                            constraint is free -- better on vertical, 3D, yaw
                            and bias at once. It stops being free above ~1.2,
                            where 3D starts paying for the yaw it buys.

        Flat over the ranges swept, so they keep OU-III's values rather than
        being fitted to these eight records: adapt_tau_sec (1.0..3.0 moves
        vertical RMS by 0.01 %), adapt_RS_mult (2..4, likewise), the pre-band
        noise floor (0.06..0.20 -- the band refers it, so the schedule barely
        notices), and the sigma band ratios themselves (0.5/4.0 is at or within
        noise of the best of five shapes tried).
    */
    float tau_coeff_    = 1.0f;
    float sigma_coeff_  = 1.0f;
    float R_S_coeff_    = 0.28f;
    float S_factor_     = 1.20f;
    float R_S_xy_factor_ = 1.15f;
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
    bool  acc_bias_hold_ = false;
    float live_since_sec_ = 0.0f;
    float proxy_gravity_good_sec_ = 0.0f;
    Vec3LPF gravity_gate_acc_lpf_{};

    // Periodic a_w marginal re-alignment; see periodicAwCovSyncTick_.
    bool  periodic_aw_cov_sync_ = true;
    float aw_sync_elapsed_sec_ = 0.0f;

    // Magnetic acquisition state.
    float mag_yaw_anchor_rad_ = 0.0f;
    bool  have_mag_yaw_anchor_ = false;
    Vector3f mag_world_ref_uT_{Vector3f::Zero()};
    bool  mag_world_ref_valid_ = false;
    float mag_init_eligible_t0_ = NAN;
    float last_mag_sample_t_ = NAN;
    float mag_north_lock_time_sec_ = NAN;
    bool  mag_refine_started_ = false;
    bool  mag_refine_done_ = false;
    float mag_refine_time_sec_ = NAN;

    // Hard iron: the startup window's contribution, the continuous estimator's
    // contribution, their sum (which is what leaves the measurement stream),
    // and the anchor the reference delta is measured from.
    Vector3f mag_hard_iron_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_startup_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_applied_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_anchor_bias_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_anchor_world_ref_uT_{Vector3f::Zero()};
    bool  mag_hi_anchored_ = false;
    float last_hi_sample_t_ = NAN;
    float last_hi_apply_t_ = NAN;

    Vector3f last_acc_body_{Vector3f::Zero()};
    Vector3f last_gyro_body_{Vector3f::Zero()};
    bool  have_last_imu_ = false;
};

}  // namespace ocean_imu::tfg
