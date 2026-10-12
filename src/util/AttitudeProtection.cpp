// Copyright 2026, Mikhail Grushinskiy
#ifndef ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "util/AttitudeProtection.h"

namespace ocean_imu::attitude {

Quat rotation(const Vec& angle) {
    const float n = angle.norm();
    if (!(n > 1e-8f)) return Quat::Identity();
    return Quat(Eigen::AngleAxisf(n, angle / n));
}

Vec AccelHold::update(const Vec& measured, bool fresh, const Vec& corrected_gyro, float dt) {
        if (fresh || !have_) { value_ = measured; have_ = measured.allFinite(); }
        else if (corrected_gyro.allFinite() && std::isfinite(dt) && dt > 0)
            value_ = rotation(Vec(-corrected_gyro * dt)) * value_;
        return value_;
    }

Decision Protection::step(float dt, const Vec& gyro, const Vec& bias, const Vec& acc, bool accel_fresh,
                  float gravity, const Quat& q_before, const Vec& field,
                  bool magnetic_valid, bool magnetic_new, uint32_t now_us,
                  float reference_yaw, bool can_reacquire) {
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
        // Bound the extra settling hold, including a requested large-angle
        // capture that continuing motion cannot qualify. Capture remains
        // pending after release; only settled observations can reset attitude.
        // Kalman adapters discard held/navigation cross-covariances when
        // ordinary aiding resumes.
        if (!dynamic_) recovery_wait_ += dt;
        if (post_dynamic_ && recovery_wait_ >= 2.0f) post_dynamic_ = false;
        d.dynamic = dynamic_;
        d.hold_attitude = holdAttitude_();
        const bool q_ok = q_before.coeffs().allFinite() && q_before.norm() > 1e-6f;
        if (!can_reacquire) { post_dynamic_ = pending_ = false; d.hold_attitude = dynamic_; }
        if (!finite || !accel_recent || acc.norm() < 1e-3f || !std::isfinite(reference_yaw) || !q_ok || !magnetic_valid || !field.allFinite() || !can_reacquire) {
            clearWindow_(); d.hold_attitude = holdAttitude_(); return d;
        }
        Quat predicted = q_before.normalized() * rotation(Vec((gyro - bias) * dt));
        predicted.normalize();
        const float fn = field.norm();
        const Vec down = -acc.normalized();
        const Vec horizontal = field - down * down.dot(field);
        const float hn = horizontal.norm();
        if (!(fn > 5 && fn < 200 && hn > 5 && hn > .1f * fn)) {
            clearWindow_(); d.hold_attitude = holdAttitude_(); return d;
        }
        const Vec north = horizontal / hn, east = down.cross(north);
        Eigen::Matrix3f axes;
        axes.col(0) = north; axes.col(1) = east; axes.col(2) = down;
        const Quat target = Quat(Eigen::AngleAxisf(reference_yaw, Vec::UnitZ())) * Quat(axes.transpose());
        const float error = predicted.angularDistance(target);
        // Give independent settled measurements a bounded capture-only window.
        // A pending request cannot permanently suppress otherwise valid aiding.
        if (!dynamic_ && error > .61f && !pending_) { // 35 degrees
            pending_ = true;
            if (!post_dynamic_) recovery_wait_ = 0;
        }
        if (error < .26f) pending_ = false;             // 15 degrees
        d.hold_attitude = holdAttitude_();
        const float dip = down.dot(field) / fn;
        const bool field_matches = !have_field_ ||
            (std::abs(fn / field_norm_ - 1) < .20f && std::abs(dip - field_dip_) < .12f);
        const bool settled = !dynamic_ && speed < .12f && norm_error < .04f && field_matches;
        if (!settled) { clearWindow_(); return d; }
        if (have_window_) { window_gyro_integral_ += gyro * dt; window_sec_ += dt; }
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
            window_target_ = target.normalized();
            window_gyro_integral_.setZero(); window_sec_ = 0;
            window_time_sum_ = window_time_square_sum_ = window_angle_square_sum_ = 0;
            window_angle_sum_.setZero(); window_time_angle_sum_.setZero();
            have_window_ = true;
        }
        last_new_us_ = now_us; ++count_;
        const Eigen::AngleAxisf observed_turn(window_target_.conjugate() * target.normalized());
        const Vec observed_angle = observed_turn.axis() * observed_turn.angle();
        window_time_sum_ += window_sec_;
        window_time_square_sum_ += window_sec_ * window_sec_;
        window_angle_sum_ += observed_angle;
        window_time_angle_sum_ += window_sec_ * observed_angle;
        window_angle_square_sum_ += observed_angle.squaredNorm();
        if (count_ < 12 || uint32_t(now_us - window_start_us_) < 750000u) return d;
        const float n = static_cast<float>(count_);
        const float time_variance = n * window_time_square_sum_ - window_time_sum_ * window_time_sum_;
        if (!(time_variance > 1e-4f) || !(window_sec_ >= .75f)) return d;
        // Fit slope with a free intercept. A short quiet arc is usable; a
        // noisy arc needs more observations before bias replacement is safe.
        const Vec observed_rate = (n * window_time_angle_sum_ -
            window_time_sum_ * window_angle_sum_) / time_variance;
        const float residual = std::max(0.0f, window_angle_square_sum_ -
            window_angle_sum_.squaredNorm()/n - observed_rate.squaredNorm()*time_variance/n);
        const float slope_variance = residual*n/((n-2)*time_variance);
        if (bias.norm() > .08f && slope_variance > .008f*.008f) return d;
        if (!have_field_) { field_norm_ = fn; field_dip_ = dip; have_field_ = true; }
        if (pending_ || post_dynamic_) {
            d.target = target.normalized(); d.reacquire = true;
            // Estimate residual gyro bias independently of the filter: raw
            // gyro average minus the rate seen by the two-vector attitude.
            // The small, settled window bounds noncommuting rotation error.
            // Subtract observed motion; low raw rate alone is not stillness.
            const Vec candidate = window_gyro_integral_ / window_sec_ - observed_rate;
            // Repair gross learned-bias corruption only. Healthy learned
            // biases are retained; this is not continuous recalibration.
            if (candidate.allFinite() && candidate.norm() < .06f && bias.norm() > .08f &&
                (candidate - bias).norm() > .08f) {
                d.repair_gyro_bias = true; d.gyro_bias = candidate;
            }
            pending_ = post_dynamic_ = false; ++recoveries_; clearWindow_();
        }
        else clearWindow_(); // Bound regression sums even during long stillness.
        return d;
    }

} // namespace ocean_imu::attitude
