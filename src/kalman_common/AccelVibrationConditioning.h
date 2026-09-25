#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Accelerometer vibration conditioning shared by the OU-II, OU-III and TFG
// orchestrators: the out-of-band vibration guard ahead of the attitude loop,
// and the vibration-aware accelerometer measurement covariance it drives.
//
// This is a public base of each wrapper so that the tuning/diagnostic API is
// literally the same object in all three.  What differs between wrappers --
// which base sigma the stage logic wants, whether low-wave noise weighting
// scales it, and when the inflation is allowed to run -- stays at each call
// site of commandRaccStd_().

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>

#include "kalman_common/SeaStateFusionDefaults.h"
#include "tuner/AccelVibrationGuard.h"

namespace seastate::common {

class AccelVibrationConditioning {
public:
    // Out-of-band accelerometer guard, ahead of the proxy and the MEKF.
    //
    // The cutoff sits in the gap between the wave band and the machinery band
    // and stops vibration reaching the attitude loop, where it rectifies into a
    // standing tilt error.  Armed by default; pass zero to remove it and
    // restore the unconditioned measurement path exactly.  The cost is group
    // delay, accelVibrationGuardDelaySec(), which appears in displacement as
    // amplitude * 2 pi f * delay -- so prefer the highest corner that removes
    // the machinery, not the lowest corner that fits above the waves.
    void setAccelVibrationGuard(float cutoff_hz,
                                int poles = defaults::ACC_VIBRATION_GUARD_POLES) {
        accel_guard_.setPoles(poles);
        accel_guard_.setCutoffHz(cutoff_hz);
    }

    // Engagement band of the guard's detector.  Exposed mainly so a study can
    // force the guard on over a quiet input and separate the delay it costs
    // from the vibration it removes.
    void setAccelVibrationEngagement(float lo_mps2, float hi_mps2,
                                     float slew_tau_sec) noexcept {
        accel_guard_.setEngagement(lo_mps2, hi_mps2, slew_tau_sec);
    }

    [[nodiscard]] float accelVibrationGuardDelaySec() const noexcept {
        return accel_guard_.groupDelaySec();
    }

    // Vibration-aware accelerometer measurement covariance.
    //
    // The guard removes the machinery it can, but what survives its stopband
    // still reaches the MEKF as measurement error the filter does not know
    // about.  This tells it: the commanded accelerometer standard deviation
    // becomes sqrt(sigma_base^2 + (gain * excess)^2), where excess is the
    // guard's detector reading above its engagement floor.
    //
    // Zero disables it and leaves the commanded covariance exactly as the
    // startup and stage logic set it.  Because the drive is the guard's own
    // gated excess, it is identically zero on a quiet installation, so an
    // enabled gain is still bit-transparent there.
    void setAccelVibrationRaccGain(float gain) noexcept {
        if (std::isfinite(gain) && gain >= 0.0f) racc_vibration_gain_ = gain;
    }

    [[nodiscard]] float accelVibrationRaccGain() const noexcept {
        return racc_vibration_gain_;
    }

    // Accelerometer sigma currently commanded to the MEKF, m/s^2 per axis.
    [[nodiscard]] Eigen::Vector3f accelVibrationRaccStd() const noexcept {
        return racc_effective_;
    }

    [[nodiscard]] float accelVibrationGuardCutoffHz() const noexcept {
        return accel_guard_.cutoffHz();
    }

    [[nodiscard]] int accelVibrationGuardPoles() const noexcept {
        return accel_guard_.poles();
    }

    // How far the guard is currently engaged, in [0, 1].  Zero means the
    // measurement path is the unconditioned one.
    [[nodiscard]] float accelVibrationGuardEngagement() const noexcept {
        return accel_guard_.engagement();
    }

    // Smoothed RMS of the out-of-band content the guard is removing, m/s^2.
    // Zero when the guard is disabled, so it reads as a health signal only
    // where it is actually measuring something.
    [[nodiscard]] float accelVibrationRms() const noexcept {
        return accel_guard_.removedRms();
    }

protected:
    // Strip out-of-band vibration from a raw accelerometer sample.  Armed by
    // default, and transparent below its detector's lower rail, in which case
    // the input is returned unchanged.
    Eigen::Vector3f conditionAccel_(const Eigen::Vector3f& acc, float dt) {
        return accel_guard_.step(acc, dt);
    }

    // Command sigma = sqrt((base * scales)^2 + (gain * excess)^2) per axis.
    //
    // scales is the low-wave noise weighting where a wrapper applies one, and
    // ones otherwise.  With no vibration excess and no scale above one, the
    // base is handed back once on the way down and the commanded covariance is
    // then left alone, so a dormant guard leaves the stage logic's covariance
    // untouched.  The caller has already checked that base is positive.
    template <typename Mekf>
    void commandRaccStd_(Mekf& mekf,
                         const Eigen::Vector3f& base,
                         const Eigen::Vector3f& scales) {
        const float excess = accel_guard_.excessRms();
        if (!(excess > 0.0f) && scales.maxCoeff() <= 1.0f) {
            if (racc_inflated_) {
                mekf.set_Racc_std(base);
                racc_effective_ = base;
                racc_inflated_ = false;
            }
            return;
        }

        const float added = racc_vibration_gain_ * excess;
        const Eigen::Vector3f effective =
            ((base.array() * scales.array()).square() + added * added).sqrt().matrix();
        mekf.set_Racc_std(effective);
        racc_effective_ = effective;
        racc_inflated_ = true;
    }

    // Armed by the owning wrapper at defaults::ACC_VIBRATION_GUARD_HZ, and
    // dormant until its own detector sees machinery, so an unconditioned
    // replay is bit-identical to a guarded one.
    seastate::tuner::AccelVibrationGuard accel_guard_{};

    // Vibration-aware measurement covariance, armed by the owning wrapper at
    // defaults::ACC_VIBRATION_RACC_GAIN and inert until the guard sees
    // machinery.
    float racc_vibration_gain_ = 0.0f;
    bool  racc_inflated_       = false;
    Eigen::Vector3f racc_effective_ = Eigen::Vector3f::Zero();
};

}  // namespace seastate::common
