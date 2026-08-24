#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  VerticalAccelComplementary - vertical acceleration from a private Mahony
  complementary filter, so the wave-period estimator can be levelled without
  reading the main filter's attitude.

  Why this exists
  ---------------
  The OU operating point is set from the zero-crossing wave period, which
  WavePeriodEstimator derives from a vertical acceleration.  Levelling that
  acceleration with the main filter's own quaternion is what puts the tuner
  inside a loop: displacing v, p, S or a_w perturbs attitude through the
  filter's cross-covariances, attitude moves the reported period, the period
  moves tau, and tau feeds back into the linear block.

  The body-Z proxy -(acc.z + g) opens that loop, being a raw measurement, but
  it is not usable: double integration weights a spectrum by 1/omega^4, so the
  gravity a tilting platform leaks into body Z swamps the elevation proxy.  Fed
  the body-Z proxy the estimator reports 6.8-10.0 s whatever the sea does,
  against a truth of 2.4-8.7 s.  That is why it was removed rather than kept
  as a fallback: the filters have no body-Z path left.

  This class is the option the filters run on: level with an attitude
  solution that is also a pure function of the measurements.  It runs its own
  Mahony observer on the raw gyro and accelerometer - never the calibrated or
  bias-corrected values, never any filter state - and reports

      up = -( (R f)_z + g ),   R = body -> NED from the private quaternion.

  Nothing it produces depends on the linear block or on the MEKF, so the
  interconnection is open by construction and stays open however the main
  filter is retuned.  Over the eight reference records it reproduces the
  attitude-levelled input to within 0.2% of vertical RMS, so the loop is
  removed at no cost, which is why it is the default rather than an option.

  Gain
  ----
  two_kp sets the rate at which the accelerometer corrects the gyro, roughly a
  first-order corner at two_kp/2 rad/s.  It must sit *below* the wave band, not
  above it: a fast observer chases the horizontal orbital acceleration into
  tilt, which is the same gravity-leakage failure by another route.  The
  default 0.2 puts the corner near 0.016 Hz, an order of magnitude under the
  0.11-0.42 Hz the reference records occupy, so the gyro carries attitude
  through the wave band and the accelerometer only trims slow drift.

  With two_ki = 0 a constant gyro bias b leaves a static tilt error of about
  2b/two_kp.  At the reference bias range (0.05 deg/s) that is ~0.5 deg, and it
  is static, so the estimator's two high-pass stages reject it.  Bias is
  deliberately not estimated here: an integral term is another slow state that
  can wind up against a sustained horizontal acceleration.

  Yaw is unobservable without a magnetometer and is left to drift.  It does not
  matter: only the vertical component is used, and that depends on tilt alone.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>

#include "ahrs/Mahony_AHRS.h"

// Which vertical acceleration the wave-period estimator is driven by.
// Complementary is the default and the only measurement-only choice; Leveled
// is kept for ablation.  The raw body-Z proxy -(acc.z + g) is no longer an
// option anywhere in the filters: every consumer of a vertical acceleration
// now reads this observer, for the reasons in the note above.
enum class WavePeriodInputSource {
    Leveled,        // heading-frame up from the main filter's attitude
    Complementary,  // private Mahony observer; levelled and measurement-only
};

class VerticalAccelComplementary {
public:
    // two_kp     : accelerometer-to-gyro correction gain; corner ~ two_kp/2 rad/s.
    // two_ki     : integral (bias) gain.  0 disables it; see the note above.
    // settle_sec : hold isReady() low this long after the first sample, so the
    //              levelling transient never reaches the period statistics.
    explicit VerticalAccelComplementary(float two_kp = 0.2f,
                                        float two_ki = 0.0f,
                                        float settle_sec = 20.0f)
        : two_kp_(two_kp), two_ki_(two_ki),
          settle_sec_(std::max(0.0f, settle_sec))
    {
        reset();
    }

    void reset() {
        ahrs_.reset(-1.0f, -1.0f);
        ahrs_.init(two_kp_, two_ki_);
        initialized_ = false;
        elapsed_sec_ = 0.0f;
        up_ms2_ = 0.0f;
    }

    void setGains(float two_kp, float two_ki) {
        two_kp_ = two_kp;
        two_ki_ = two_ki;
        ahrs_.init(two_kp_, two_ki_);
    }

    // gyro : body angular rate [rad/s], raw.
    // acc  : body specific force [m/s^2], raw.  acc.z ~ -g at rest, matching
    //        the z-down body frame the rest of the filter uses.
    void update(float dt_sec,
                const Eigen::Vector3f& gyro,
                const Eigen::Vector3f& acc,
                float gravity_ms2)
    {
        if (!(dt_sec > 0.0f) || !std::isfinite(dt_sec) ||
            !gyro.allFinite() || !acc.allFinite()) {
            return;
        }

        // Mahony drives the measured vector toward the world z axis expressed
        // in body coordinates.  Here world z is NED down, and a specific force
        // of (0,0,-g) at rest points up, so the observer is fed -acc.
        const float acc_norm = acc.norm();
        if (!initialized_) {
            if (!(acc_norm > 1e-3f)) return;
            seed_from_acc_(acc / acc_norm);
            initialized_ = true;
        }

        float pitch_deg = 0.0f, roll_deg = 0.0f, yaw_deg = 0.0f;
        ahrs_.update(gyro.x(), gyro.y(), gyro.z(),
                     -acc.x(), -acc.y(), -acc.z(),
                     &pitch_deg, &roll_deg, &yaw_deg,
                     dt_sec);

        // Third row of R(body->NED); dotting it with the specific force gives
        // the down component, and it is also the down axis in body coordinates.
        const float w = ahrs_.q0, x = ahrs_.q1, y = ahrs_.q2, z = ahrs_.q3;
        const Eigen::Vector3f down_row(2.0f * (x * z - w * y),
                                       2.0f * (y * z + w * x),
                                       w * w - x * x - y * y + z * z);

        // a_ned.z = (R f)_z + g, and up is its negation.
        up_ms2_ = -(down_row.dot(acc) + gravity_ms2);

        elapsed_sec_ += dt_sec;
    }

    // Up-positive vertical acceleration with gravity removed [m/s^2].
    float verticalAccelUpMs2() const noexcept { return up_ms2_; }

    bool isReady() const noexcept {
        return initialized_ && elapsed_sec_ >= settle_sec_;
    }

    // Seeded from the first accelerometer sample, so tilt is usable well
    // before isReady() -- that gate exists to keep the levelling transient out
    // of the *period statistics*, which is a stricter requirement than having
    // a usable attitude.  A startup bootstrap wants the attitude as soon as it
    // exists and applies its own quality gate, so it reads this instead.
    bool isInitialized() const noexcept { return initialized_; }

    float elapsedSec() const noexcept { return elapsed_sec_; }

    // BODY -> NED attitude of the private observer.
    //
    // Yaw is unobservable without a magnetometer and drifts, so only the tilt
    // part of this is meaningful; see tiltQuaternion().  Exposed whole because
    // the levelled heading frame that the wave-direction stage runs in is
    // invariant under q -> Rz(psi) q, so it can consume this directly.
    Eigen::Quaternionf quaternion() const noexcept {
        const Eigen::Quaternionf q(ahrs_.q0, ahrs_.q1, ahrs_.q2, ahrs_.q3);
        const float n = q.norm();
        if (!(n > 1e-8f) || !std::isfinite(n)) {
            return Eigen::Quaternionf::Identity();
        }
        return Eigen::Quaternionf(q.coeffs() / n);
    }

    // Tilt part of the observer attitude, with the drifting heading divided
    // out.  Invariant under q -> Rz(psi) q, which is what makes it usable as a
    // magnetometer accumulation frame: no arbitrary yaw can leak into the
    // learned world reference through it.
    Eigen::Quaternionf tiltQuaternion() const noexcept {
        const Eigen::Quaternionf q = quaternion();

        const Eigen::Matrix3f R = q.toRotationMatrix();
        const float yaw = std::atan2(R(1, 0), R(0, 0));
        if (!std::isfinite(yaw)) return Eigen::Quaternionf::Identity();

        const Eigen::Quaternionf q_yaw_inv(
            Eigen::AngleAxisf(-yaw, Eigen::Vector3f::UnitZ()));

        Eigen::Quaternionf q_tilt = q_yaw_inv * q;
        q_tilt.normalize();

        if (!q_tilt.coeffs().allFinite()) return Eigen::Quaternionf::Identity();
        return q_tilt;
    }

private:
    // Point the private world z axis along the measured gravity direction, so
    // the observer starts levelled instead of rotating in from identity.
    void seed_from_acc_(const Eigen::Vector3f& acc_unit) {
        // Down in body coordinates is opposite the specific force at rest.
        const Eigen::Vector3f down_body = -acc_unit;
        // q maps body -> NED, so it must carry down_body onto +Z.  Yaw is left
        // at whatever FromTwoVectors picks; it is unobservable and unused.
        const Eigen::Quaternionf q = Eigen::Quaternionf::FromTwoVectors(
            down_body, Eigen::Vector3f::UnitZ());
        ahrs_.q0 = q.w();
        ahrs_.q1 = q.x();
        ahrs_.q2 = q.y();
        ahrs_.q3 = q.z();
    }

    Mahony_AHRS<float> ahrs_{};
    float two_kp_;
    float two_ki_;
    float settle_sec_;
    bool  initialized_ = false;
    float elapsed_sec_ = 0.0f;
    float up_ms2_ = 0.0f;
};