#pragma once

// Proxy-bootstrap startup, magnetic acquisition and output pipeline shared by
// the deployed OU-II and OU-III front ends (SeaStateFusion_OU_II and
// SeaStateFusion_OU_III).
//
// Startup: the proxy owns tilt and magnetic learning.
//
//   Bootstrap  The MEKF is held.  The inner filter runs its measurement-only
//              front end (updateFrontEnd): the private Mahony observer levels,
//              the tuner learns the operating point, and the magnetometer is
//              averaged in the observer's yaw-stripped tilt frame into a
//              provisional reference and a yaw gauge.
//   Handoff    Once the observer's tilt holds world-frame gravity agreement,
//              north is gauged (when fitted) and the tuner is ready -- or on a
//              timeout that still requires the aligned branch -- the MEKF is
//              seeded with proxy tilt + magnetic yaw and goes Live.
//   Live       updateTime drives the MEKF; magnetometer corrections reach it.
//              A second acquisition re-learns the reference once the observer
//              has converged; continuous hard-iron estimation then runs.
//
// Everything in this class is the same for both OU families.  The estimator
// enters in exactly one place: the seed call that constructs the MEKF with
// the family's own pseudo-measurement variances, which each wrapper passes to
// beginStartup_().  Filter is the inner SeaStateFusionFilter_OU_* type and
// Config derives from ProxyStartupFusionConfig.

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>

#include "detrend/AdaptiveWaveDetrender3D.h"
#include "kalman_common/MagneticStartupCommon.h"
#include "kalman_common/SeaStateFusionDefaults.h"
#include "kalman_common/SeaStateFusionFilterCommon.h"
#include "tuner/ContinuousMagHardIronEstimator.h"
#include "tuner/MagAutoTuner.h"

// Nominal standard gravity: the generic default of gravity_magnitude.  A
// calibrated deployment passes its local gravity (g_cal_local) instead.
extern const float g_std;

namespace seastate::common {

struct ProxyStartupFusionConfig {
    bool with_mag = true;

    // Earliest and latest the proxy bootstrap may hand over.
    //
    // The normal exit is by quality: proxy tilt holding gravity agreement,
    // magnetic north gauged, and the tuner ready.  The floor keeps a record
    // whose first seconds happen to look calm from handing over on a tilt the
    // observer has barely integrated, and the ceiling guarantees the filter
    // always starts -- a platform that never satisfies the gate still gets a
    // live filter, on the best attitude available, rather than sitting in
    // bootstrap forever.
    float proxy_startup_min_sec     = defaults::PROXY_STARTUP_MIN_SEC;
    float proxy_startup_timeout_sec = defaults::PROXY_STARTUP_TIMEOUT_SEC;

    // Magnetic acquisition runs in two stages, because the two things it has
    // to deliver want opposite schedules.
    //
    // A usable heading is wanted within seconds of power-on.  A *good*
    // reference wants the startup observer to have settled first: its
    // correction corner sits below the wave band by design, which is what
    // stops it chasing orbital acceleration, and the same low corner means it
    // needs tens of seconds to converge from its accelerometer seed.  Measured
    // against truth on the eight reference records (OU-III), mean tilt error
    // over a 23 s averaging window starting at 7 s runs 0.33 to 2.76 deg;
    // starting at 40 s it is 0.13 to 0.85 deg.
    //
    // Waiting for the good one before reporting anything would put first
    // heading around 105 s, which is not a usable device.  So the first stage
    // locks a provisional reference as soon as the gravity gate allows --
    // heading and a live filter in roughly 20 s -- and the second stage
    // re-learns the reference once the observer has actually converged, again
    // in the observer's own tilt frame (see maybeRefineMagReference_() for why
    // not the MEKF's).  The correction lands long before the scored window
    // opens.
    //
    // proxy_mag_settle_sec holds the provisional stage off; 0 means "as soon
    // as the gravity gate is happy" and is the default, because the refinement
    // is what carries the accuracy now.
    float proxy_mag_settle_sec = defaults::PROXY_MAG_SETTLE_SEC;

    bool  mag_refine_enabled    = true;
    float mag_refine_start_sec  = 90.0f;
    float mag_refine_window_sec = defaults::MAG_REFINE_WINDOW_SEC;

    // Covariance seeded at handoff.  Tilt has been integrated through the wave
    // band by an observer whose correction corner is below it, so it is worth
    // about the accel-only default; yaw is either gauged by the magnetometer
    // or entirely unknown, and those two cases are an order of magnitude
    // apart, which is the whole reason the seed is split.
    float proxy_handoff_tilt_sigma_rad     = defaults::PROXY_HANDOFF_TILT_SIGMA_RAD;
    float proxy_handoff_yaw_sigma_rad      = defaults::PROXY_HANDOFF_YAW_SIGMA_RAD;
    float proxy_handoff_yaw_sigma_free_rad = defaults::PROXY_HANDOFF_YAW_SIGMA_FREE_RAD;

    float mag_delay_sec          = defaults::MAG_DELAY_SEC;
    float online_tune_warmup_sec = defaults::STARTUP_ONLINE_TUNE_WARMUP_SEC;

    // Magnetometer updates that must land after the filter goes live before
    // the accelerometer-bias gate opens.
    //
    // Accelerometer bias and a tilt error are only weakly separable in waves
    // -- a roll error tips gravity into body Y and reads as a Y bias -- so
    // opening this gate while the attitude is still settling lets the bias
    // absorb the error and hold it.  The proxy bootstrap reaches live far
    // earlier than a staged warmup would, which moves this gate earlier in
    // absolute terms unless it is set to account for that.
    int acc_bias_unlock_mag_updates = 250;

    Eigen::Vector3f sigma_a = Eigen::Vector3f(0.2f, 0.2f, 0.2f);
    Eigen::Vector3f sigma_g = Eigen::Vector3f(0.01f, 0.01f, 0.01f);
    Eigen::Vector3f sigma_m = Eigen::Vector3f(0.3f, 0.3f, 0.3f);

    // MEKF variances common to both OU constructors.  They were reachable only
    // through initialize_ext(), which the wrappers once never called, so every
    // deployment ran on the header defaults, which a dedicated MEKF-variance
    // sweep gauged (OU-III).  The values here reproduce those defaults
    // exactly.  The pseudo-measurement variances are family-specific and live
    // in each family's Config.
    //
    //   Pq0      initial attitude-error variance, rad^2.  The proxy handoff
    //            overwrites the attitude block, so this only ever seeds the
    //            covariance the handoff replaces.
    //   Pb0      initial gyro-bias variance, (rad/s)^2.
    //   b0       gyro-bias random-walk variance density, (rad/s)^2/s.
    float Pq0 = 5e-4f;
    float Pb0 = 1e-6f;
    float b0  = 1e-10f;

    // Physical gravity removed from the specific force.  The generic default
    // is nominal standard gravity; a calibrated deployment sets the
    // calibration site's local gravity (g_cal_local), and the MEKF, the
    // observer and the magnetic acquisition then all use that value.
    float gravity_magnitude = g_std;

    // Period-scaled sigma-band shape.  These are dimensionless wave-band
    // ratios except for the absolute safety clamps.  Keeping the ratios fixed
    // is what gives the JONSWAP sigma channel its similarity law.
    float sigma_band_low_ratio  = defaults::SIGMA_BAND_LOW_RATIO;
    float sigma_band_high_ratio = defaults::SIGMA_BAND_HIGH_RATIO;
    float sigma_band_min_hz     = defaults::SIGMA_BAND_MIN_HZ;
    float sigma_band_max_hz     = defaults::SIGMA_BAND_MAX_HZ;

    // Mag-start gate: gravity-direction agreement using current tilt,
    // measured on the world-frame specific force averaged over the wave band;
    // see WorldGravityGate.
    float mag_gravity_align_max_sin  = defaults::GRAVITY_GATE_MAX_SIN;  // sin(deg)
    float mag_gravity_align_hold_sec = defaults::GRAVITY_GATE_HOLD_SEC;

    // Horizon of the world-frame average the gate is judged on, and how long
    // that average must have been running before its verdict counts.
    //
    // The average has to span whole wave periods for orbital acceleration to
    // cancel out of it, so the horizon is set against the longest swell the
    // device is expected to start up in rather than against the sea it
    // happens to be in -- the frequency tracker is not converged this early,
    // and being conservative here costs settling time rather than accuracy.
    // 12 s covers the band these filters work in.
    //
    // The warmup is set so the gate's earliest possible verdict lands with the
    // magnetometer's first eligible sample -- mag_delay_sec less the hold the
    // gate has to serve anyway -- which is the last moment at which it is
    // free.  Beyond that it delays a calm start for nothing: the magnetometer
    // cannot begin averaging before mag_delay_sec however early the gate
    // closes, so any warmup shorter than this buys no time and any warmup
    // longer than it costs time one-for-one.
    float mag_gravity_align_world_tau_sec    = defaults::GRAVITY_GATE_LPF_SEC;
    float mag_gravity_align_world_warmup_sec = defaults::GRAVITY_GATE_WARMUP_SEC;

    float mag_tilt_fallback_sec = defaults::MAG_TILT_FALLBACK_SEC;
    float mag_extreme_gyro_dps  = defaults::MAG_EXTREME_GYRO_DPS;  // veto only truly violent motion
    float mag_init_min_mag_norm = defaults::MAG_INIT_MIN_MAG_NORM;

    // Real-device mag acquisition.
    //
    // The important rule: MagAutoTuner must never receive an estimator yaw,
    // which is arbitrary and unobservable before mag lock.  It receives the
    // startup observer's quaternion with yaw divided out instead: that tilt is
    // invariant under q_bw -> Rz(psi) q_bw, so no heading can leak through it,
    // and it is a far better level frame in waves than one rebuilt from accel.
    // It then returns a gauge-fixed reference B_ref = [horizontal, 0,
    // vertical] and the yaw gauge that puts learned north on +X.
    //
    // The reference is an average of the field in that tilt frame, so
    // whatever tilt error survives the window survives in the reference.  In
    // waves the error is periodic, so what the window has to buy is whole wave
    // periods, not samples: 128 samples is 5.1 s at a 25 Hz mag ODR, short
    // enough to lock in the phase it started on rather than cancel it.  15 s
    // covers a couple of periods across the band these filters work in and
    // captures most of what a much longer window would, at a startup cost of
    // 15 s rather than 40 s.  Held in seconds so it does not silently shorten
    // at a higher ODR.
    int   mag_min_samples    = defaults::MAG_MIN_SAMPLES;
    float mag_min_window_sec = defaults::MAG_MIN_WINDOW_SEC;
    float mag_max_window_sec = defaults::MAG_MAX_WINDOW_SEC;  // no forced timeout
    float mag_sample_dt_sec  = defaults::MAG_SAMPLE_DT_SEC;

    // Body-frame hard-iron offset, learned during startup alongside the
    // reference and then subtracted from every magnetometer sample.
    //
    // Off by default.  The MEKF has no mag-bias state, so an offset left in
    // the stream is heading error one-for-one against the horizontal field;
    // but the offset is only weakly separable from the reference at a fixed
    // heading, and a wrong one subtracted everywhere is worse than none.  Turn
    // it on where the platform changes heading during startup.
    bool mag_estimate_hard_iron = false;

    // Continuous hard-iron estimation and the reference that goes with it;
    // see ContinuousHardIronTracker.  Nothing here starts until the two-stage
    // startup acquisition has finished, so first heading, handoff and
    // refinement are unaffected by it.
    bool  mag_continuous_hard_iron    = true;
    float mag_hi_memory_sec           = defaults::MAG_HI_MEMORY_SEC;
    // Absolute floor only; the relative term below carries the calibration.
    // See ContinuousMagHardIronEstimator::Config and
    // docs/continuous-mag-hard-iron.md for why it came down from 4e-3.
    float mag_hi_model_ridge          = defaults::MAG_HI_MODEL_RIDGE;
    float mag_hi_model_ridge_relative = defaults::MAG_HI_MODEL_RIDGE_RELATIVE;
    float mag_hi_min_information      = defaults::MAG_HI_MIN_INFORMATION;
    float mag_hi_min_effective_weight = defaults::MAG_HI_MIN_EFFECTIVE_WEIGHT;
    float mag_hi_max_residual_rms_uT  = defaults::MAG_HI_MAX_RESIDUAL_RMS_UT;
    float mag_hi_max_bias_fraction    = defaults::MAG_HI_MAX_BIAS_FRACTION;

    // Fraction of the fitted offset the filter is willing to apply, and the
    // time constant it moves over.  The ridge already shrinks the fit for what
    // the model cannot see; this is the separate, blunter statement that a
    // calibration nobody has checked should not arrive as a step.
    float mag_hi_apply_fraction = defaults::MAG_HI_APPLY_FRACTION;
    float mag_hi_slew_tau_sec   = defaults::MAG_HI_SLEW_TAU_SEC;

    // Keep off in waves.  Accel/gyro weighting can phase-select wave motion.
    bool  mag_enable_quality_weighting = false;
    float mag_min_effective_weight     = 0.0f;
    float mag_acc_norm_rel_soft        = defaults::MAG_ACC_NORM_REL_SOFT;
    float mag_gyro_soft_dps            = defaults::MAG_GYRO_SOFT_DPS;

    bool enable_displacement_detrend = false;
    bool use_custom_displacement_detrend_cfg = false;
    AdaptiveWaveDetrender3D::Config displacement_detrend_cfg{};
};

template <typename Filter, typename Config>
class ProxyStartupFusion {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    void update(float dt,
                const Eigen::Vector3f& gyro_body_ned,
                const Eigen::Vector3f& acc_body_ned,
                float tempC = 35.0f)
    {
        if (!begun_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;

        t_ += dt;

        last_acc_body_ned_  = acc_body_ned;
        last_gyro_body_ned_ = gyro_body_ned;
        have_last_imu_      = true;

        // The gravity-lock bootstrap is the Mahony observer inside the front
        // end, which runs from the first sample, so there is no phase in which
        // the wrapper withholds IMU data.
        if (stage_ != Stage::Live) {
            // Bootstrap: front end only, MEKF held.
            impl_.updateFrontEnd(dt, gyro_body_ned, acc_body_ned);
        } else {
            impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);
        }

        // Whose tilt the magnetometer gate is judged against.  Before handoff
        // this is the observer's, so the gate measures the attitude that will
        // actually frame the magnetic reference rather than one the MEKF is
        // still converging toward.  The gate reads the raw accelerometer: it
        // averages over twelve seconds, far below the machinery band the
        // inner filter's vibration guard removes.
        gravity_gate_.step(attitudeReferenceQuat_(), acc_body_ned, gyro_body_ned, dt,
                           cfg_.mag_gravity_align_world_tau_sec,
                           cfg_.mag_gravity_align_world_warmup_sec,
                           cfg_.mag_gravity_align_max_sin,
                           cfg_.mag_extreme_gyro_dps);

        updateDisplacementOutput_(dt);

        const auto cur_stage = impl_.getStartupStage();
        if (cur_stage != last_impl_startup_stage_) {
            if (cur_stage == Filter::StartupStage::Cold) {
                mag_ref_set_ = false;
                mag_auto_tuner_.reset();
                gravity_gate_.reset();
                mag_init_eligible_t0_ = NAN;
                last_mag_sample_t_ = NAN;

                last_mag_tilt_frame_yaw_rad_ = NAN;
                last_mag_startup_yaw_correction_rad_ = NAN;

                if (stage_ != Stage::Live) {
                    // Inner filter already re-locked tilt internally.
                    displacement_up_m_.setZero();
                    displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
                    if (cfg_.enable_displacement_detrend) {
                        displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
                    }
                }
            }
            last_impl_startup_stage_ = cur_stage;
        }

        if (stage_ != Stage::Live) {
            maybeHandOffToMekf_();
        }
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!begun_ || !cfg_.with_mag) return;
        // There is no withheld stage to wait out: the observer has been
        // levelling since the first sample, and learning north before the MEKF
        // starts is the entire point.  The quality gate below still decides
        // when accumulation may begin.
        if (t_ < cfg_.mag_delay_sec) return;

        // Ahead of every gate below, and deliberately.  The continuous
        // estimator wants the whole magnetometer record, not the part the
        // startup machinery was willing to average, and it is reading a frame
        // the startup machinery does not own.
        accumulateContinuousHardIron_(mag_body_ned);

        // Hold the whole magnetometer path off until the startup observer has
        // settled.  This sits ahead of the eligibility clock deliberately, so
        // the tilt-fallback timer cannot start running and then wave the
        // accumulation through on an attitude that is still converging.
        //
        // What this window does *not* fix is the standing yaw error.  That
        // residual is a property of the learned reference vector rather than
        // of the one-time gauge: seeding the handoff with a yaw variance from
        // 5 deg to 90 deg moves the scored yaw by under 1e-4 deg (OU-III), so
        // the filter is not converging to the gauge, it is converging to the
        // reference.  maybeRefineMagReference_() is what re-learns it.
        if (!mag_ref_set_ && t_ < cfg_.proxy_mag_settle_sec) {
            return;
        }

        if (!std::isfinite(mag_init_eligible_t0_)) {
            mag_init_eligible_t0_ = t_;
        }

        // The first stage accepts the held agreement alone; the branch is
        // enforced where the tilt is actually committed, at handoff.
        const bool gravity_trusted =
            (gravity_gate_.good_sec >= cfg_.mag_gravity_align_hold_sec);

        const bool fallback_ok =
            ((t_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec);

        if (!mag_ref_set_) {
            if (!gravity_trusted && !fallback_ok) {
                return;
            }

            if (have_last_imu_) {
                const float dt_mag = advanceSampleClock(
                    last_mag_sample_t_, t_, cfg_.mag_sample_dt_sec);

                // Accumulate in a tilt frame with yaw removed.
                //
                // Stripping yaw makes the frame invariant to the estimator's
                // arbitrary startup heading, so this leaks no yaw into the
                // learned reference -- q_bw and Rz(psi) q_bw give the same
                // tilt.  A gravity-only frame rebuilt from low-passed accel
                // would be yaw-free too, but in waves that accel is gravity
                // plus a phase-lagged remnant of the orbital specific force,
                // so its tilt is wrong by a wave-correlated angle that the
                // averaging window is too short to cancel.
                //
                // Before handoff the tilt comes from the private Mahony
                // observer, which is gyro-propagated through the wave band and
                // corrected below it.  The provisional reference and the yaw
                // gauge are locked once here; the refinement stage replaces
                // them.
                const Eigen::Quaternionf q_tilt_bw =
                    yawRemovedBoatQuat(attitudeReferenceQuat_());

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
                        // frame, so it carries no estimator heading.  It is a
                        // model parameter, so writing it before the MEKF has
                        // been handed the attitude is safe and leaves it ready
                        // to use the magnetometer from its first live sample.
                        setMagWorldRef_(mag_world_ref_uT);

                        const float mag_tilt_yaw_rad =
                            mag_auto_tuner_.getYawGaugeCorrectionRad();

                        if (std::isfinite(mag_tilt_yaw_rad)) {
                            // One-time yaw-gauge lock.
                            //
                            // mag_tilt_yaw_rad is the heading of averaged
                            // magnetic north in the accumulation frame, so
                            // driving the boat's absolute yaw to its negative
                            // puts the learned north on +X.  Only yaw is
                            // written: the tilt is the gyro-propagated,
                            // accel-corrected estimate, and overwriting the
                            // whole quaternion would inject the wave tilt
                            // error the filter has already rejected.
                            const float yaw_abs_rad = wrapPi(-mag_tilt_yaw_rad);

                            // Before handoff there is no MEKF attitude to
                            // rewrite; the gauge is carried to the handoff
                            // instead and composed with the proxy tilt there,
                            // so the filter's very first attitude already has
                            // north in it.
                            if (stage_ != Stage::Live) {
                                pending_yaw_abs_rad_ = yaw_abs_rad;
                                last_mag_tilt_frame_yaw_rad_ = wrapPi(mag_tilt_yaw_rad);
                                last_mag_startup_yaw_correction_rad_ = yaw_abs_rad;
                            } else {
                                writeMekfYaw_(yaw_abs_rad, mag_tilt_yaw_rad);
                            }
                        }

                        Eigen::Vector3f hard_iron_uT;
                        hard_iron_.startup_body_uT =
                            mag_auto_tuner_.getHardIronBodyUT(hard_iron_uT)
                                ? hard_iron_uT
                                : Eigen::Vector3f::Zero();
                        mag_hard_iron_body_uT_ = hard_iron_.startup_body_uT;

                        mag_ref_set_ = true;
                        mag_north_lock_time_sec_ = t_;
                    }
                }
            }
        }

        maybeRefineMagReference_(mag_body_ned);
        maybeApplyContinuousHardIron_();

        // Magnetometer corrections go to the MEKF only once it owns the
        // attitude.  Before handoff its state is not the one being solved.
        if (mag_ref_set_ && stage_ == Stage::Live) {
            impl_.updateMag(mag_body_ned - mag_hard_iron_body_uT_);
        }
    }

    bool hasMagNorthLock() const noexcept { return mag_ref_set_; }

    // True once the second-stage acquisition has replaced the provisional
    // reference; see maybeRefineMagReference_().
    bool hasRefinedMagReference() const noexcept { return mag_refine_done_; }
    float magRefineTimeSec() const noexcept { return mag_refine_time_sec_; }

    // Wall-clock marks for the startup sequence, for anyone measuring
    // time-to-first-fix rather than steady-state accuracy.
    float magNorthLockTimeSec() const noexcept { return mag_north_lock_time_sec_; }
    float liveTimeSec() const noexcept { return live_time_sec_; }

    // Body-frame hard-iron offset removed from the magnetometer stream.  Zero
    // unless Config::mag_estimate_hard_iron asked for it and the startup window
    // constrained it well enough to use, or the continuous estimator has
    // applied a correction.
    const Eigen::Vector3f& magHardIronBodyUT() const noexcept {
        return mag_hard_iron_body_uT_;
    }

    // Continuous estimator, for diagnostics and for tests that need to see why
    // it did or did not act.  estimate().valid is the gate; information() is
    // how much attitude excitation the memory window has actually collected.
    const ContinuousMagHardIronEstimator& magContinuousHardIron() const noexcept {
        return hard_iron_.estimator;
    }

    // The part of magHardIronBodyUT() the continuous estimator is responsible
    // for, as opposed to the startup solve.
    const Eigen::Vector3f& magContinuousHardIronAppliedUT() const noexcept {
        return hard_iron_.applied_body_uT;
    }

    int magAcceptedCount() const noexcept { return mag_auto_tuner_.acceptedCount(); }
    int magRejectedCount() const noexcept { return mag_auto_tuner_.rejectedCount(); }
    float magAcceptedWindowSec() const noexcept { return mag_auto_tuner_.acceptedWindowSec(); }
    float magEffectiveWeight() const noexcept { return mag_auto_tuner_.effectiveWeight(); }

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

    bool  isLive() const { return stage_ == Stage::Live; }
    float freqHz() const { return impl_.getFreqHz(); }
    float waveDirectionDeg() const { return impl_.getWaveDirectionDeg(); }

    // Best available boat attitude, BODY -> WORLD (NED).
    //
    // The MEKF holds its initial quaternion until the handoff, so reading it
    // directly during the bootstrap would report a level identity attitude
    // rather than the platform's.  The bootstrap observer's solution is what
    // is actually known at that point, and it is what this returns; after
    // handoff this returns the MEKF's.
    //
    // Heading is only meaningful once hasMagNorthLock() is true (or the build
    // has no magnetometer); before that the bootstrap yaw is arbitrary.  The
    // linear outputs -- displacement, velocity -- stay gated on isLive().
    Eigen::Quaternionf attitudeQuat() const {
        return attitudeReferenceQuat_();
    }

    const Eigen::Vector3f& displacementUpMeters() const { return displacement_up_m_; }
    const AdaptiveWaveDetrender3D::Output& displacementDetrend() const { return displacement_det_out_; }

    Filter& raw() { return impl_; }
    const Filter& raw() const { return impl_; }

protected:
    // Configure and start.  seed(impl, cfg) constructs the estimator: it is
    // the one estimator-specific step of the startup sequence, and it runs at
    // the point the sequence has always constructed it -- after the inner
    // filter's wiring and before the nominal accelerometer sigma is recorded.
    template <typename SeedFn>
    void beginStartup_(const Config& cfg, SeedFn&& seed) {
        cfg_ = cfg;

        begun_ = true;
        stage_ = Stage::Bootstrap;
        t_ = 0.0f;

        gravity_gate_.reset();
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;

        mag_ref_set_ = false;

        last_mag_tilt_frame_yaw_rad_ = NAN;
        last_mag_startup_yaw_correction_rad_ = NAN;

        mag_auto_tuner_.setConfig(magAutoTunerConfig(cfg_));
        hard_iron_.configure(cfg_);
        hard_iron_.resetApplied();
        mag_hard_iron_body_uT_.setZero();

        mag_world_ref_uT_.setZero();
        mag_world_ref_valid_ = false;

        last_acc_body_ned_.setZero();
        last_gyro_body_ned_.setZero();
        have_last_imu_ = false;

        pending_yaw_abs_rad_ = NAN;

        mag_refine_started_  = false;
        mag_refine_done_     = false;
        mag_refine_time_sec_ = NAN;
        mag_north_lock_time_sec_ = NAN;
        live_time_sec_           = NAN;

        impl_.setWithMag(cfg_.with_mag);
        impl_.setMagUpdatesToUnlockAccBias(cfg_.acc_bias_unlock_mag_updates);

        // The provisional reference is deliberately cheap and early, so the
        // accelerometer bias must not be allowed to fit itself to it; see the
        // inner filter's setAccBiasHold().
        impl_.setAccBiasHold(cfg_.with_mag && cfg_.mag_refine_enabled);
        impl_.setMagDelaySec(0.0f);  // outer wrapper owns the magnetometer delay
        impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);
        impl_.setSigmaWaveBandRatios(cfg_.sigma_band_low_ratio,
                                     cfg_.sigma_band_high_ratio);
        impl_.setSigmaWaveBandLimitsHz(cfg_.sigma_band_min_hz,
                                       cfg_.sigma_band_max_hz);

        seed(impl_, cfg_);
        last_impl_startup_stage_ = impl_.getStartupStage();

        impl_.setNominalRaccStd(cfg_.sigma_a);

        displacement_up_m_.setZero();
        displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

        if (cfg_.enable_displacement_detrend) {
            if (cfg_.use_custom_displacement_detrend_cfg) {
                displacement_detrender_.setConfig(cfg_.displacement_detrend_cfg);
            } else {
                displacement_detrender_.setConfig(
                    defaultDisplacementDetrenderConfig<AdaptiveWaveDetrender3D::Config>(FREQ_GUESS));
            }
            displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
        }
    }

    Filter impl_{false};
    Config cfg_{};

private:
    enum class Stage {
        Bootstrap,   // proxy owns the attitude, MEKF held
        Live         // MEKF owns the attitude
    };

    // Up-positive displacement from the MEKF's NED position, optionally
    // passed through the adaptive detrender with the tracker frequency as its
    // external carrier.  Disabled, the detrender output is the raw signal.
    void updateDisplacementOutput_(float dt) {
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

    // The attitude the startup machinery judges and frames things against.
    //
    // Once the MEKF is live it is the answer -- it has the magnetometer, the
    // linear block and the bias states, and the proxy has none of them.
    // Before that the MEKF has nothing to say and the observer does.
    Eigen::Quaternionf attitudeReferenceQuat_() const {
        if (stage_ != Stage::Live) {
            return impl_.startupProxyQuat();
        }
        return impl_.mekf().quaternion_boat();
    }

    // Absolute yaw write on the live MEKF, keeping its tilt.
    void writeMekfYaw_(float yaw_abs_rad, float mag_tilt_yaw_rad) {
        Eigen::Quaternionf q_bw = impl_.mekf().quaternion_boat();
        q_bw.normalize();

        const Eigen::Quaternionf q_new = boatQuatWithAbsoluteYaw(q_bw, yaw_abs_rad);

        if (q_new.coeffs().allFinite()) {
            impl_.mekf().set_quaternion_boat(q_new);
            last_mag_tilt_frame_yaw_rad_ = wrapPi(mag_tilt_yaw_rad);
            last_mag_startup_yaw_correction_rad_ = yaw_abs_rad;
        }
    }

    // Every write of the magnetometer's world reference goes through here, so
    // the wrapper always knows the vector the MEKF is steering to.  The MEKF
    // does not offer it back, and the continuous correction moves the
    // reference by a delta rather than replacing it.
    void setMagWorldRef_(const Eigen::Vector3f& mag_world_ref_uT) {
        impl_.mekf().set_mag_world_ref(mag_world_ref_uT);
        mag_world_ref_uT_ = mag_world_ref_uT;
        mag_world_ref_valid_ = true;
    }

    // Raw magnetometer in the observer's tilt frame; see
    // ContinuousHardIronTracker::accumulate().
    void accumulateContinuousHardIron_(const Eigen::Vector3f& mag_body_ned) {
        if (!cfg_.mag_continuous_hard_iron) return;
        if (!impl_.startupProxyInitialized()) return;
        hard_iron_.accumulate(t_, cfg_.mag_sample_dt_sec,
                              impl_.startupProxyTiltQuat(), mag_body_ned);
    }

    // Applied only once the startup acquisition, including refinement, has
    // finished and the MEKF is live.
    void maybeApplyContinuousHardIron_() {
        if (!cfg_.mag_continuous_hard_iron) return;
        if (!mag_ref_set_ || stage_ != Stage::Live) return;
        if (cfg_.mag_refine_enabled && !mag_refine_done_) return;
        if (!mag_world_ref_valid_) return;

        Eigen::Vector3f ref;
        if (!hard_iron_.slewTowardEstimate(t_, cfg_.mag_sample_dt_sec,
                                           mag_world_ref_uT_,
                                           cfg_.mag_hi_apply_fraction,
                                           cfg_.mag_hi_slew_tau_sec,
                                           cfg_.mag_init_min_mag_norm,
                                           ref)) {
            return;
        }
        mag_hard_iron_body_uT_ = hard_iron_.startup_body_uT + hard_iron_.applied_body_uT;
        setMagWorldRef_(ref);
    }

    // Second-stage magnetic acquisition.
    //
    // The provisional reference was averaged in a tilt frame the startup
    // observer had barely converged, and it is what the filter has been
    // steering to ever since.  This re-runs the same acquisition once the MEKF
    // is live and the observer has settled.
    //
    // Both the reference vector and the heading gauge are replaced.  The yaw
    // write is a step, and deliberately so: it is the coarse-to-fine alignment
    // correction, it happens once, and it lands well before the scored window
    // opens.
    void maybeRefineMagReference_(const Eigen::Vector3f& mag_body_ned) {
        if (!cfg_.mag_refine_enabled) return;
        if (mag_refine_done_) return;
        if (!mag_ref_set_) return;
        if (stage_ != Stage::Live) return;
        if (t_ < cfg_.mag_refine_start_sec) return;
        if (!have_last_imu_) return;

        if (!mag_refine_started_) {
            beginMagRefinement(mag_auto_tuner_, cfg_.mag_refine_window_sec,
                               cfg_.mag_min_samples);
            mag_refine_started_ = true;
            last_mag_sample_t_  = NAN;
        }

        const float dt_mag = advanceSampleClock(last_mag_sample_t_, t_, cfg_.mag_sample_dt_sec);

        // The observer's tilt, not the MEKF's, for the same reason the first
        // stage used it -- and here the reason has teeth.
        //
        // By now the MEKF has the linear block and the bias states, so its
        // tilt looks like the better frame.  It is not usable: the MEKF has
        // been steering to the provisional reference this pass exists to
        // replace, so its tilt carries that reference's error, and averaging
        // the field in it re-derives the error it was meant to remove.  Tried
        // that way the refinement is self-confirming -- reference and yaw come
        // back within 1e-3 deg of the provisional ones, and the standing roll
        // bias is untouched.
        //
        // The observer never saw the reference, so its tilt is independent of
        // it, and by refinement time it has long since converged.
        const Eigen::Quaternionf q_tilt_bw = impl_.startupProxyTiltQuat();

        // Feed the same corrected stream the MEKF sees, so a hard-iron offset
        // already removed is not re-learned into the new reference.
        const Eigen::Vector3f mag_corrected = mag_body_ned - mag_hard_iron_body_uT_;

        if (!mag_auto_tuner_.addSampleWithTiltQuatDt(
                dt_mag, q_tilt_bw, last_acc_body_ned_,
                last_gyro_body_ned_, mag_corrected)) {
            return;
        }

        Eigen::Vector3f mag_world_ref_uT;
        if (!mag_auto_tuner_.getMagWorldRef(mag_world_ref_uT) ||
            !mag_world_ref_uT.allFinite() ||
            !(mag_world_ref_uT.norm() > cfg_.mag_init_min_mag_norm)) {
            return;
        }

        setMagWorldRef_(mag_world_ref_uT);

        const float mag_tilt_yaw_rad = mag_auto_tuner_.getYawGaugeCorrectionRad();
        if (std::isfinite(mag_tilt_yaw_rad)) {
            writeMekfYaw_(wrapPi(-mag_tilt_yaw_rad), mag_tilt_yaw_rad);
        }

        mag_refine_done_     = true;
        mag_refine_time_sec_ = t_;

        // The reference the bias would have been fitting is now the good one.
        impl_.setAccBiasHold(false);
    }

    // Hand the bootstrap attitude to the MEKF and start it live.
    //
    // The quality exit needs three things at once: a tilt that has held
    // agreement with gravity long enough to be trusted, a magnetic north gauge
    // (when a magnetometer is fitted), and an operating point the tuner
    // stands behind.  Waiting for all three is what lets the MEKF skip the
    // staged warmup entirely -- there is nothing left for it to converge.
    void maybeHandOffToMekf_() {
        if (!begun_) return;

        const bool proxy_ready = impl_.startupProxyInitialized();

        const bool tilt_trusted = gravity_gate_.trusted(cfg_.mag_gravity_align_hold_sec);

        const bool north_ready = !cfg_.with_mag || mag_ref_set_;

        const bool ready_by_quality =
            proxy_ready &&
            (t_ >= cfg_.proxy_startup_min_sec) &&
            tilt_trusted &&
            north_ready &&
            impl_.isTunerReady();

        // The timeout still requires an attitude to hand over; without one
        // there is nothing to seed and waiting costs nothing.
        //
        // It is also held clear of the magnetometer acquisition it would
        // otherwise cut short.  A timeout that fires while the reference is
        // still averaging hands over with no yaw gauge at all, which is a far
        // worse start than simply waiting: the gauge is the one chance to put
        // the filter on north before it goes live.  So the floor is the settle
        // time plus room for the averaging window to close.
        const float mag_acquire_deadline =
            cfg_.with_mag
                ? cfg_.proxy_mag_settle_sec +
                      2.0f * std::max(cfg_.mag_min_window_sec, 1.0f) +
                      cfg_.mag_tilt_fallback_sec
                : 0.0f;

        const float timeout_sec =
            std::max(cfg_.proxy_startup_timeout_sec, mag_acquire_deadline);

        // The timeout bounds how long startup may take; it does not license a
        // handoff onto the antipodal branch.  Seeding the MEKF with an attitude
        // that disagrees with measured gravity by more than a right angle is
        // worse than waiting out the extra samples it takes the observer to
        // leave a set it is not attracted to in the first place.
        const bool ready_by_timeout =
            proxy_ready &&
            (t_ >= timeout_sec) &&
            gravity_gate_.aligned_branch;

        if (!ready_by_quality && !ready_by_timeout) return;

        handOffToMekf_();
    }

    void handOffToMekf_() {
        const bool have_yaw_gauge = std::isfinite(pending_yaw_abs_rad_);

        // boatQuatWithAbsoluteYaw strips the incoming heading before writing
        // the new one, so the observer's drifted yaw cannot survive this even
        // though the quaternion is passed in whole.
        const Eigen::Quaternionf q_proxy = impl_.startupProxyQuat();

        const Eigen::Quaternionf q_seed =
            have_yaw_gauge
                ? boatQuatWithAbsoluteYaw(q_proxy, pending_yaw_abs_rad_)
                : q_proxy;

        if (!q_seed.coeffs().allFinite()) return;

        const float yaw_sigma = have_yaw_gauge
            ? cfg_.proxy_handoff_yaw_sigma_rad
            : cfg_.proxy_handoff_yaw_sigma_free_rad;

        // The accelerometer-bias gate is deliberately left closed here.  Going
        // live early does not make that bias any more observable in waves, so
        // it keeps waiting for its count of magnetometer updates exactly as it
        // did before; see the inner filter's setMagUpdatesToUnlockAccBias().
        impl_.goLive(q_seed,
                     cfg_.proxy_handoff_tilt_sigma_rad,
                     yaw_sigma,
                     /*allow_acc_bias=*/false);

        stage_ = Stage::Live;
        live_time_sec_ = t_;
        last_impl_startup_stage_ = impl_.getStartupStage();
    }

    bool begun_ = false;

    Stage stage_ = Stage::Bootstrap;
    float t_ = 0.0f;

    typename Filter::StartupStage last_impl_startup_stage_ = Filter::StartupStage::Cold;

    // Last IMU sample for mag-init gating.
    Eigen::Vector3f last_acc_body_ned_  = Eigen::Vector3f::Zero();
    Eigen::Vector3f last_gyro_body_ned_ = Eigen::Vector3f::Zero();
    bool have_last_imu_ = false;

    // Mag-init state.
    bool mag_ref_set_ = false;
    Eigen::Vector3f mag_hard_iron_body_uT_ = Eigen::Vector3f::Zero();
    MagAutoTuner mag_auto_tuner_{};

    ContinuousHardIronTracker hard_iron_{};

    Eigen::Vector3f mag_world_ref_uT_ = Eigen::Vector3f::Zero();
    bool mag_world_ref_valid_ = false;

    float last_mag_sample_t_ = NAN;

    float last_mag_tilt_frame_yaw_rad_ = NAN;
    float last_mag_startup_yaw_correction_rad_ = NAN;

    // Yaw gauge acquired while the MEKF was still held, applied at handoff.
    float pending_yaw_abs_rad_ = NAN;

    float mag_north_lock_time_sec_ = NAN;
    float live_time_sec_           = NAN;

    bool  mag_refine_started_  = false;
    bool  mag_refine_done_     = false;
    float mag_refine_time_sec_ = NAN;

    AdaptiveWaveDetrender3D displacement_detrender_{};
    AdaptiveWaveDetrender3D::Output displacement_det_out_{};
    Eigen::Vector3f displacement_up_m_ = Eigen::Vector3f::Zero();

    WorldGravityGate gravity_gate_{};
    float mag_init_eligible_t0_ = NAN;
};

}  // namespace seastate::common
