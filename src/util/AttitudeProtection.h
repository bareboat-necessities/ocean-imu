#pragma once

// Copyright 2026, Mikhail Grushinskiy
// Shared device policy. These are conservative engineering thresholds, not a
// sensor fault classifier or a proof that a stationary magnetic field is Earth.
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#else
#include <ArduinoEigenDense.h>
#endif
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace ocean_imu::attitude {
using Vec = Eigen::Vector3f;
using Quat = Eigen::Quaternionf;

inline Quat rotation(const Vec& angle) {
    const float n = angle.norm();
    if (!(n > 1e-8f)) return Quat::Identity();
    return Quat(Eigen::AngleAxisf(n, angle / n));
}

// A gyro-only frame still advances attitude. A held accelerometer vector must
// be expressed in the new body axes before continuous observers reuse it.
class AccelHold {
public:
    void reset() { have_ = false; value_.setZero(); }
    Vec update(const Vec& measured, bool fresh, const Vec& corrected_gyro, float dt) {
        if (fresh || !have_) { value_ = measured; have_ = measured.allFinite(); }
        else if (corrected_gyro.allFinite() && std::isfinite(dt) && dt > 0)
            value_ = rotation(Vec(-corrected_gyro * dt)) * value_;
        return value_;
    }
private:
    bool have_ = false;
    Vec value_ = Vec::Zero();
};

struct Decision {
    bool dynamic = false;
    bool hold_attitude = false;
    bool reacquire = false;
    Quat target = Quat::Identity();
};

class Protection {
public:
    void reset() { *this = Protection{}; }

    // q_before is BODY->world before this tick's prediction; field is the
    // calibrated, offset-corrected field at the current gyro epoch in uT.
    // reference_yaw is the azimuth of the existing world magnetic reference.
    // No reference, calibration, navigation state or bias is relearned here.
    Decision step(float dt, const Vec& gyro, const Vec& bias, const Vec& acc, bool accel_fresh,
                  float gravity, const Quat& q_before, const Vec& field,
                  bool magnetic_valid, bool magnetic_new, uint32_t now_us,
                  float reference_yaw = 0.0f, bool can_reacquire = true) {
        Decision d;
        // AccelHold transports a vector; it does not create gravity evidence.
        // Allow the normal 100 Hz accel / 200 Hz gyro phase offset, but discard
        // any capture window interrupted by more than 20 ms without a sample.
        const bool accel_gap = !have_accel_ || uint32_t(now_us - last_accel_us_) > 20000u;
        if (accel_gap) clearWindow_();
        if (accel_fresh && acc.allFinite() && acc.squaredNorm() > 1e-6f) {
            have_accel_ = true; last_accel_us_ = now_us;
        }
        const bool accel_recent = have_accel_ && uint32_t(now_us - last_accel_us_) <= 20000u;
        const bool finite = std::isfinite(dt) && dt > 0 && dt < 1.0f &&
            gyro.allFinite() && acc.allFinite() && bias.allFinite() &&
            std::isfinite(gravity) && gravity > 0;
        const float norm_error = finite ? std::abs(acc.norm() / gravity - 1) : 1;
        const float speed = finite ? gyro.norm() : 100;
        // A ramp can create tangential acceleration well before the rate
        // threshold is reached. Filter both rate and its derivative at 30 ms to avoid
        // treating sample-to-sample gyro noise as angular acceleration.
        float angular_acceleration = 0;
        if (finite) {
            if (have_rate_) {
                const Vec derivative = (gyro - rate_lpf_) / (.03f + dt);
                rate_lpf_ += dt * derivative;
                angular_accel_lpf_ += (dt / (.03f + dt)) * (derivative - angular_accel_lpf_);
                angular_acceleration = angular_accel_lpf_.norm();
            } else { rate_lpf_ = gyro; have_rate_ = true; }
        }
        // The rate gate also covers acceleration whose norm happens to equal
        // g. A norm-only test cannot distinguish that from true gravity.
        const bool hazard = !finite || !accel_recent || dt > .020001f || speed > 1.5f || angular_acceleration > 2.0f || norm_error > .20f;
        if (hazard) { dynamic_ = true; post_dynamic_ = true; clear_sec_ = recovery_wait_ = 0; }
        else if (dynamic_) {
            clear_sec_ = speed < .75f && angular_acceleration < 1.0f && norm_error < .08f ? clear_sec_ + dt : 0;
            if (clear_sec_ >= .5f) dynamic_ = false;
        }
        // A vessel need not become completely stationary after a turn.
        // Bound the extra settling hold; a known large attitude error still
        // requires independent magnetic capture. Kalman adapters discard the
        // held/navigation cross-covariances when ordinary aiding resumes.
        if (!dynamic_) recovery_wait_ += dt;
        if (post_dynamic_ && !pending_ && recovery_wait_ >= 2.0f) post_dynamic_ = false;
        d.dynamic = dynamic_;
        d.hold_attitude = dynamic_ || post_dynamic_;
        const bool q_ok = q_before.coeffs().allFinite() && q_before.norm() > 1e-6f;
        if (!can_reacquire) { post_dynamic_ = pending_ = false; d.hold_attitude = dynamic_; }
        if (!finite || !accel_recent || acc.norm() < 1e-3f || !std::isfinite(reference_yaw) || !q_ok || !magnetic_valid || !field.allFinite() || !can_reacquire) {
            clearWindow_(); d.hold_attitude = dynamic_ || pending_ || post_dynamic_; return d;
        }
        Quat predicted = q_before.normalized() * rotation(Vec((gyro - bias) * dt));
        predicted.normalize();
        const float fn = field.norm();
        const Vec down = -acc.normalized();
        const Vec horizontal = field - down * down.dot(field);
        const float hn = horizontal.norm();
        if (!(fn > 5 && fn < 200 && hn > 5 && hn > .1f * fn)) {
            clearWindow_(); d.hold_attitude = dynamic_ || pending_ || post_dynamic_; return d;
        }
        const Vec north = horizontal / hn, east = down.cross(north);
        Eigen::Matrix3f axes;
        axes.col(0) = north; axes.col(1) = east; axes.col(2) = down;
        const Quat target = Quat(Eigen::AngleAxisf(reference_yaw, Vec::UnitZ())) * Quat(axes.transpose());
        const float error = predicted.angularDistance(target);
        // Outside a local correction's capture region, wait for independent
        // settled measurements instead of learning the discrepancy as bias.
        if (!dynamic_ && error > .61f) pending_ = true; // 35 degrees
        if (error < .26f) pending_ = false;             // 15 degrees
        d.hold_attitude = dynamic_ || pending_ || post_dynamic_;
        const float dip = down.dot(field) / fn;
        const bool field_matches = !have_field_ ||
            (std::abs(fn / field_norm_ - 1) < .20f && std::abs(dip - field_dip_) < .12f);
        const bool settled = !dynamic_ && speed < .12f && norm_error < .04f && field_matches;
        if (!settled) { clearWindow_(); return d; }
        // Both source streams must remain fresh throughout the window.
        // Gyro transport alone never establishes recovery observations.
        if (have_window_ && uint32_t(now_us - last_new_us_) > 120000u) clearWindow_();
        if (!magnetic_new) return d;
        const bool consistent = have_window_ &&
            (down - window_down_).norm() < .05f &&
            (field / fn - window_field_).norm() < .05f &&
            std::abs(fn / window_norm_ - 1) < .06f;
        if (!consistent) {
            window_start_us_ = now_us; count_ = 0;
            window_down_ = down; window_field_ = field / fn; window_norm_ = fn;
            have_window_ = true;
        }
        last_new_us_ = now_us; ++count_;
        if (count_ < 12 || uint32_t(now_us - window_start_us_) < 750000u) return d;
        if (!have_field_) { field_norm_ = fn; field_dip_ = dip; have_field_ = true; }
        if (pending_ || post_dynamic_) {
            d.target = target.normalized(); d.reacquire = true;
            pending_ = post_dynamic_ = false; ++recoveries_; clearWindow_();
        }
        return d;
    }
    uint32_t recoveryCount() const { return recoveries_; }
    template<class Print> void print(Print& out) const {
        out.printf("[ATT] hold=%d dynamic=%d recoveries=%lu\n",
            int(dynamic_ || pending_ || post_dynamic_), int(dynamic_), static_cast<unsigned long>(recoveries_));
    }
private:
    void clearWindow_() { have_window_ = false; count_ = 0; }
    bool post_dynamic_ = false;
    bool have_accel_ = false;
    bool have_rate_ = false;
    Vec rate_lpf_ = Vec::Zero(), angular_accel_lpf_ = Vec::Zero();
    bool dynamic_ = false, pending_ = false, have_field_ = false, have_window_ = false;
    float recovery_wait_ = 0;
    float clear_sec_ = 0, field_norm_ = 0, field_dip_ = 0, window_norm_ = 0;
    Vec window_down_ = Vec::Zero(), window_field_ = Vec::Zero();
    uint32_t last_accel_us_ = 0, window_start_us_ = 0, last_new_us_ = 0, count_ = 0, recoveries_ = 0;
};
} // namespace ocean_imu::attitude
