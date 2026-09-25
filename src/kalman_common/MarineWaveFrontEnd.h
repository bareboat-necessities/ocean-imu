#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Measurement-only marine front end of the OU-II and OU-III orchestrators.
//
//     accelerometer (vibration-conditioned) + gyro
//       -> private Mahony observer (VerticalAccelComplementary)
//            -> levelled vertical acceleration, startup attitude proxy
//       -> acceleration-band frequency tracker + stillness detector
//            -> fast/slow smoothed carrier frequency (wave direction)
//       -> WavePeriodEstimator (zero-crossing period of the elevation)
//            -> wave-band tuning frequency, or its fixed prior
//       -> vessel-RAO equalizer -> KalmanWaveDirection -> WaveDirectionDetector
//
// Nothing here reads an estimator state except the one input the direction
// branch is given on purpose: the attitude to level with and, once the MEKF
// is driven, its accelerometer-bias estimate.  The wave-period and tuner
// inputs stay independent of both; see updateWavePeriodAndDirection_().
//
// The component is a public base of each OU wrapper so that the whole
// front-end API is one definition.  The order in which its steps run relative
// to the MEKF, the tuner and the covariance maintenance is part of each
// estimator's schedule and is spelled out in the wrapper's updateCore_().

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>

#include "freq/FirstOrderIIRSmoother.h"
#include "freq/FrequencyTrackerPolicy.h"
#include "kalman_common/SeaStateFusionDefaults.h"
#include "tuner/SeaStateFusionTunerCommon.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/WavePeriodEstimator.h"
#include "wave_dir/KalmanWaveDirection.h"
#include "wave_dir/VesselRaoEqualizer.h"
#include "wave_dir/VesselRaoNoiseWeighting.h"
#include "wave_dir/WaveDirectionDetector.h"
#include "wave_dir/WaveDirectionFrame.h"

// Nominal standard gravity, used here only as the unit of the stillness
// detector's energy (a_vert expressed in g).  Physical gravity for levelling
// is gravity_mps2_ below.
extern const float g_std;

namespace seastate::common {

template <TrackerType trackerT>
class MarineWaveFrontEnd {
public:
    using TrackingPolicy = TrackerPolicy<trackerT>;
    using FreqInputLPF = seastate::tuner::common::FreqInputLPF;
    using StillnessAdapter = seastate::tuner::common::StillnessAdapter;

    MarineWaveFrontEnd()
        : freq_hz_(FREQ_GUESS),
          freq_hz_slow_(FREQ_GUESS)
    {
        // Default cutoff ~max_freq_hz_ Hz: passes waves, kills 8-37 Hz engine band
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);
    }

    // Attitude of the startup Mahony observer, BODY -> NED.  Only its tilt is
    // meaningful; yaw is unobservable to it and drifts.
    Eigen::Quaternionf startupProxyQuat() const noexcept {
        return vertical_accel_comp_.quaternion();
    }

    // Tilt-only form of the same attitude, safe to use as a magnetometer
    // accumulation frame because no heading can leak through it.
    Eigen::Quaternionf startupProxyTiltQuat() const noexcept {
        return vertical_accel_comp_.tiltQuaternion();
    }

    bool startupProxyInitialized() const noexcept {
        return vertical_accel_comp_.isInitialized();
    }

    // Gains of the private Mahony observer, which serves both the vertical
    // channel and the startup attitude.  Same knob as
    // setWavePeriodComplementaryGains(); kept under this name because the
    // startup path is the one with an opinion about two_ki.
    void setStartupProxyGains(float two_kp, float two_ki) {
        vertical_accel_comp_.setGains(two_kp, two_ki);
    }

    // Gains of the private Mahony observer that levels the default input.
    // two_kp sets the accelerometer-to-gyro correction corner, which must stay
    // below the wave band; see VerticalAccelComplementary.h.
    void setWavePeriodComplementaryGains(float two_kp, float two_ki) {
        vertical_accel_comp_.setGains(two_kp, two_ki);
    }

    // Configure LPF on the levelled vertical acceleration the tracker runs on.
    void setFreqInputCutoffHz(float fc) { freq_input_lpf_.setCutoff(fc); }

    // Bounds of the acceleration-band tracker.  This is a wave-direction knob:
    // the OU operating point reads the wave band; see setTuneFreqBounds().
    void setFreqBounds(float min_hz, float max_hz) {
        if (!std::isfinite(min_hz) || !std::isfinite(max_hz)) return;
        if (min_hz <= 0.0f || max_hz <= min_hz) return;
        min_freq_hz_ = min_hz;
        max_freq_hz_ = max_hz;
        freq_stillness_.setTargetFreqHz(min_freq_hz_);
    }

    inline float getFreqHz()     const noexcept { return freq_hz_; }        // fast branch
    inline float getFreqSlowHz() const noexcept { return freq_hz_slow_; }   // slow branch
    inline float getFreqRawHz()  const noexcept { return f_raw; }

    // Apparent period of the *acceleration* band, from the slow tracker
    // branch.  This is a reporting channel: the OU operating point uses
    // getWavePeriodSec(), the zero-crossing period of the elevation, which is
    // a different and much longer quantity.
    inline float getPeriodSec() const noexcept {
        return (freq_hz_slow_ > 1e-6f) ? 1.0f / freq_hz_slow_ : NAN;
    }

    // Up-positive vertical acceleration from the private Mahony observer: the
    // signal the tracker, the wave-period estimator and the sigma channel all
    // run on.  Measurement-only -- it reads no filter state.
    inline float getAccelVertical() const noexcept {
        return vertical_accel_comp_.verticalAccelUpMs2();
    }

    // Zero-crossing wave period [s] from the independent accelerometer-only
    // estimator; NaN until it settles.
    inline float getWavePeriodSec() const noexcept { return wave_period_.getPeriodSec(); }
    inline bool wavePeriodUsable() const noexcept { return wave_period_.hasUsablePeriod(); }
    inline bool wavePeriodReady() const noexcept { return wave_period_.isReady(); }

    // Select which vertical acceleration drives the wave-period estimator.
    // Complementary (default) levels with the private Mahony observer and is
    // measurement-only, so the tuner is outside the estimator's loop.  Leveled
    // is the older behaviour, which levels with the main filter's attitude and
    // closes that loop.  See updateWavePeriodAndDirection_() for what it costs.
    void setWavePeriodInput(WavePeriodInputSource source) {
        wave_period_input_ = source;
    }
    WavePeriodInputSource wavePeriodInput() const noexcept {
        return wave_period_input_;
    }

    // Wave-band frequency used before WavePeriodEstimator has a value.
    void setTuneFreqPriorHz(float hz) {
        if (std::isfinite(hz) && hz > 0.0f) tune_freq_prior_hz_ = hz;
    }
    float tuneFreqPriorHz() const noexcept { return tune_freq_prior_hz_; }

    void setLowWaveNoiseWeighting(const wave_direction::VesselRaoNoiseWeighting& config) {
        low_wave_noise_ = config;
    }

    void setDirectionRao(const wave_direction::VesselRaoEqualizer::Config& config) {
        direction_rao_.configure(config);
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

    inline KalmanWaveDirection& dir() noexcept { return dir_filter_; }
    inline const KalmanWaveDirection& dir() const noexcept { return dir_filter_; }

    inline WaveDirectionDetector<float>& dir_sign() noexcept { return dir_sign_; }
    inline const WaveDirectionDetector<float>& dir_sign() const noexcept { return dir_sign_; }

protected:
    // Private Mahony observer for the wave-period estimator, the default
    // sigma channel and the startup attitude.  It is fed gyro and
    // accelerometer before the MEKF sees them, so that the levelling it
    // provides stays a pure function of the measurements.  Stepping it
    // unconditionally keeps its transient off the critical path when the
    // input source is switched at runtime.
    void levelVertical_(float dt, const Eigen::Vector3f& gyro,
                        const Eigen::Vector3f& acc_in) {
        vertical_accel_comp_.update(dt, gyro, acc_in, gravity_mps2_);
    }

    // The one levelled vertical measurement every consumer in the filter
    // reads: the frequency tracker and its stillness detector, the wave-period
    // estimator, and the sigma channel.  It comes from the private Mahony
    // observer and therefore reads no MEKF state.  The observer is seeded from
    // the first accelerometer sample, so the value is usable immediately; its
    // isReady() gate is about settled period *statistics*, which is a stricter
    // requirement than a usable tilt.
    //
    // Runs the acceleration-band tracker on it and returns it.
    float trackAccelBandFrequency_(float dt) {
        const float a_vert_measurement = vertical_accel_comp_.verticalAccelUpMs2();

        // LPF on the tracker input
        const float a_vert_lp = freq_input_lpf_.step(a_vert_measurement, dt);

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
        freq_hz_slow_ = f_slow;   // reporting (getPeriodSec)

        return a_vert_measurement;
    }

    // The frequency the whole adaptation path runs on: the sigma band's
    // corners, the sigma_a averaging horizon, and tau.  It is a wave-band
    // quantity and nothing else; the acceleration-band tracker never reaches
    // it, at any instant of the run.
    //
    // The fixed wave-band prior is used only until the existing period
    // estimator's own startup-usable gate clears.  That gate is earlier than
    // strict isReady(), but it still requires four leak time constants and at
    // least one estimated cycle of the same moment/log-period state.  There is
    // no second startup estimator and no estimator handoff.
    float tuner_frequency_hz_() const {
        const float wave_hz = wave_period_.getFrequencyHz();
        if (wave_period_.hasUsablePeriod() &&
            std::isfinite(wave_hz) && wave_hz > 0.0f) {
            return wave_hz;
        }
        return tune_freq_prior_hz_;
    }

    // Zero-crossing wave period, then wave direction.
    //
    // The period estimator runs beside the frequency tracker rather than
    // replacing it: the tracker supplies the acceleration-band carrier the
    // direction demodulator needs, while the OU operating point needs the
    // wave band.  Its input must be levelled and must not read estimator
    // state.
    //
    // Levelled, because double integration weights a spectrum by 1/omega^4,
    // so sub-band gravity leakage in a raw body-Z residual can dominate the
    // elevation proxy.
    //
    // Exogenous, because levelling with the filter's own attitude closes a
    // feedback loop through the tuner and the filter's cross-covariances.
    // VerticalAccelComplementary instead levels with a private Mahony observer
    // driven only by sensor measurements.  setWavePeriodInput() selects the
    // main-filter attitude-levelled input for ablation;
    // tests/kalman_ou_iii/tuner_coupling-test.cpp checks the default's independence
    // from MEKF state and bounds the Leveled path's gain.
    //
    // Direction is resolved in a leveled frame aligned with boat heading.
    // This removes roll/pitch mixing while preserving 0 deg = bow and positive
    // angles toward starboard.  Stage 1 estimates the apparent propagation
    // plane as an unsigned axis relative to boat heading; stage 2 resolves
    // propagation sense along that same axis from horizontal/vertical orbital
    // phase.
    //
    // q_dir is the attitude to level with: the MEKF's once it is driven, the
    // proxy's before.  heading_frame_acceleration() resolves into the
    // projected bow axis and is therefore invariant under q -> Rz(psi) q, so
    // only the tilt of whichever quaternion is supplied can reach the result,
    // and the handoff is continuous in everything this stage can see.
    //
    // acc_bias_body, when given, is the MEKF's body accelerometer-bias
    // estimate: the direction branch needs inertial acceleration, so it is
    // corrected before levelling, as the MEKF's acceleration model does.  The
    // period/tuner input above stays independent of that estimate.
    void updateWavePeriodAndDirection_(float dt,
                                       const Eigen::Quaternionf& q_dir,
                                       const Eigen::Vector3f& acc_in,
                                       const Eigen::Vector3f* acc_bias_body)
    {
        const float omega = 2.0f * static_cast<float>(M_PI) * freq_hz_;

        const auto direction_accel = wave_direction::heading_frame_acceleration<float>(
            q_dir, acc_in, gravity_mps2_);

        wave_period_.update(dt, wave_period_input_ms2_(direction_accel));

        auto direction_input = direction_accel;
        if (acc_bias_body) {
            const Eigen::Vector3f corrected = acc_in - *acc_bias_body;
            direction_input = wave_direction::heading_frame_acceleration<float>(
                q_dir, corrected, gravity_mps2_);
        }
        const auto direction_matched = direction_rao_.step(direction_input, dt);
        dir_filter_.update(direction_matched.forward_ms2,
                           direction_matched.starboard_ms2,
                           omega, dt);
        const Eigen::Vector2f propagation_axis_boat = dir_filter_.getAxis();
        dir_sign_state_ = dir_sign_.update(
            direction_matched.forward_ms2,
            direction_matched.starboard_ms2,
            direction_matched.up_ms2,
            propagation_axis_boat.x(), propagation_axis_boat.y(),
            dt, dir_filter_.getLastStableConfidence());
    }

    // Vertical acceleration the wave-period estimator is driven by.  The
    // leveled ablation falls back to the complementary observer while heading
    // is not yet resolved, so the estimator is never fed a body-frame residual.
    float wave_period_input_ms2_(
        const wave_direction::HeadingFrameAcceleration<float>& leveled) const
    {
        const float a_comp = vertical_accel_comp_.verticalAccelUpMs2();
        switch (wave_period_input_) {
            case WavePeriodInputSource::Leveled:
                return leveled.heading_valid ? leveled.up_ms2 : a_comp;
            case WavePeriodInputSource::Complementary:
            default:
                return a_comp;
        }
    }

    bool frontEndStill_() const { return freq_stillness_.isStill(); }
    float frontEndStillTimeSec_() const { return freq_stillness_.getStillTime(); }

    void resetFrontEnd_() {
        tracker_policy_ = TrackingPolicy{};
        wave_period_    = WavePeriodEstimator{};
        vertical_accel_comp_.reset();
        freq_input_lpf_ = FreqInputLPF{};
        freq_stillness_ = StillnessAdapter(g_std, min_freq_hz_, FREQ_GUESS);
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);

        freq_fast_smoother_ = FirstOrderIIRSmoother<float>(defaults::NOMINAL_IMU_DT_S, 3.5f);
        freq_slow_smoother_ = FirstOrderIIRSmoother<float>(defaults::NOMINAL_IMU_DT_S, 10.0f);

        freq_hz_      = FREQ_GUESS;
        freq_hz_slow_ = FREQ_GUESS;
        f_raw         = FREQ_GUESS;

        dir_filter_ = KalmanWaveDirection(2.0f * static_cast<float>(M_PI) * FREQ_GUESS);
        dir_sign_.reset();
        direction_rao_.reset();
        dir_sign_state_ = UNCERTAIN;
    }

    float freq_hz_      = FREQ_GUESS;
    float freq_hz_slow_ = FREQ_GUESS;
    float f_raw         = FREQ_GUESS;

    float min_freq_hz_        = defaults::ACCEL_TRACKER_MIN_FREQ_HZ;
    float max_freq_hz_        = defaults::ACCEL_TRACKER_MAX_FREQ_HZ;
    float tune_freq_prior_hz_ = defaults::TUNE_FREQ_PRIOR_HZ;
    WavePeriodInputSource wave_period_input_ = WavePeriodInputSource::Complementary;

    TrackingPolicy               tracker_policy_{};
    FirstOrderIIRSmoother<float> freq_fast_smoother_{defaults::NOMINAL_IMU_DT_S, 3.5f};
    FirstOrderIIRSmoother<float> freq_slow_smoother_{defaults::NOMINAL_IMU_DT_S, 10.0f};
    WavePeriodEstimator          wave_period_;

    // One private Mahony observer, serving both the vertical channel and the
    // startup attitude.
    //
    // Integral feedback estimates gyro bias so it does not remain as a static
    // tilt error in the startup attitude or levelled acceleration.
    // two_kp keeps the correction corner below the wave band, so the observer
    // does not level itself against orbital specific force instead of gravity.
    VerticalAccelComplementary   vertical_accel_comp_{
        defaults::STARTUP_PROXY_TWO_KP,
        defaults::STARTUP_PROXY_TWO_KI};

    // Physical gravity removed from the specific force: the MEKF's gravity
    // (set by the wrapper's initialize_ext; initialize() keeps the MEKF
    // default, the standard gravity of the Kalman3D_Wave_OU_* classes).
    float gravity_mps2_ = 9.80665f;

    FreqInputLPF     freq_input_lpf_;
    StillnessAdapter freq_stillness_;

    KalmanWaveDirection dir_filter_{2.0f * static_cast<float>(M_PI) * FREQ_GUESS};
    wave_direction::VesselRaoEqualizer direction_rao_{};
    wave_direction::VesselRaoNoiseWeighting low_wave_noise_{};
    WaveDirectionDetector<float> dir_sign_{0.002f, 0.005f};
    WaveDirection                dir_sign_state_ = UNCERTAIN;
};

}  // namespace seastate::common
