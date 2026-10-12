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

/** @brief Return a quaternion from a rotation vector in radians, or identity near zero. */
Quat rotation(const Vec& angle);

// A gyro-only frame still advances attitude. A held accelerometer vector must
// be expressed in the new body axes before continuous observers reuse it.
class AccelHold {
public:
    /** @brief Discard held acceleration and its availability flag. */
    void reset() { have_ = false; value_.setZero(); }
    /** @brief Cache fresh body acceleration or rotate the held vector into the new body axes.
     * Use bias-corrected gyro in rad/s and dt in seconds; transport adds no fresh evidence.
     */
    Vec update(const Vec& measured, bool fresh, const Vec& corrected_gyro, float dt);
private:
    bool have_ = false;
    Vec value_ = Vec::Zero();
};

// INS must estimate legitimate wave force, rather than latch on its norm.
enum class MotionUse { Compass, INS };

// Compose every raw gyro increment between fresh specific-force packets.
// The frontend remains measurement-only; each estimator removes its own bias.
class GyroWindow {
public:
    void reset() { turn_ = Quat::Identity(); elapsed_ = 0; valid_ = true; }
    void push(const Vec& gyro, float dt);
    void take(Vec& gyro, float& dt);
private:
    Quat turn_ = Quat::Identity();
    float elapsed_ = 0;
    bool valid_ = true;
};

struct Decision {
    bool dynamic = false;
    bool hold_attitude = false;
    bool reacquire = false;
    bool repair_gyro_bias = false;
    Vec gyro_bias = Vec::Zero();
    Quat target = Quat::Identity();
    // hold_attitude protects force-derived attitude feedback and bias learning.
    // It does not invalidate magnetic observations. Their health, timestamp,
    // freshness and innovation gates remain independent of vessel motion.
};

class Protection {
public:
    /** @brief Reset motion history, capture qualification, field baseline, and recovery count. */
    explicit Protection(MotionUse use = MotionUse::Compass) : use_(use) {}
    void reset() { *this = Protection{use_}; }
    // Observe all gyro packets, including opposing rates that a window cancels.
    void observeGyro(float dt, const Vec& gyro);

    // q_before is BODY->world before this tick's prediction; field is the
    // calibrated, offset-corrected field at the current gyro epoch in uT.
    // reference_yaw is the azimuth of the existing world magnetic reference.
    // References, calibration and navigation means are retained. A grossly
    // corrupt gyro bias may be repaired from a longer qualified window.
    /** @brief Evaluate motion hazards and return an aiding hold or qualified recovery target.
     * Use dt in seconds, gyro/bias in rad/s, acceleration/gravity in m/s^2,
     * and now_us as a wrapping microsecond clock. Recovery requires fresh,
     * consistent gravity and magnetic observations; target is body-to-world.
     */
    Decision step(float dt, const Vec& gyro, const Vec& bias, const Vec& acc, bool accel_fresh,
                  float gravity, const Quat& q_before, const Vec& field,
                  bool magnetic_valid, bool magnetic_new, uint32_t now_us,
                  float reference_yaw = 0.0f, bool can_reacquire = true);
    /** @brief Return the number of recovery targets issued since the last reset. */
    uint32_t recoveryCount() const { return recoveries_; }
    /** @brief Print the hold, dynamic-motion, and recovery-count diagnostics. */
    template<class Print> void print(Print& out) const {
        out.printf("[ATT] hold=%d dynamic=%d recoveries=%lu\n",
            int(holdAttitude_()), int(dynamic_), static_cast<unsigned long>(recoveries_));
    }
private:
    // A requested capture must not starve ordinary aiding indefinitely when
    // continuing motion prevents the settled observation window from closing.
    bool holdAttitude_() const {
        return dynamic_ || (pending_ && gross_bias_) ||
            ((pending_ || post_dynamic_) && recovery_wait_ < 2.0f);
    }
    /** @brief Discard accumulated recovery observations and restart the sample count. */
    void clearWindow_() { have_window_ = false; count_ = 0; }
    MotionUse use_;
    bool observed_ = false, gyro_hazard_ = false, gyro_quiet_ = true;
    float peak_speed_ = 0;
    bool post_dynamic_ = false, gross_bias_ = false;
    bool have_accel_ = false;
    bool have_rate_ = false;
    Vec rate_lpf_ = Vec::Zero(), angular_accel_lpf_ = Vec::Zero();
    bool dynamic_ = false, pending_ = false, have_field_ = false, have_window_ = false;
    float recovery_wait_ = 0;
    float clear_sec_ = 0, field_norm_ = 0, field_dip_ = 0, window_norm_ = 0;
    Vec window_down_ = Vec::Zero(), window_field_ = Vec::Zero();
    Quat window_target_ = Quat::Identity();
    Vec window_gyro_integral_ = Vec::Zero();
    float window_sec_ = 0;
    float window_time_sum_ = 0, window_time_square_sum_ = 0, window_angle_square_sum_ = 0;
    Vec window_angle_sum_ = Vec::Zero(), window_time_angle_sum_ = Vec::Zero();
    uint32_t last_accel_us_ = 0, window_start_us_ = 0, last_new_us_ = 0, count_ = 0, recoveries_ = 0;
};
} // namespace ocean_imu::attitude
