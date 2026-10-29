#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// IMU installation lever arm of the experimental CARMA filters
// (Kalman3D_Wave_CARMA3, Kalman3D_Wave_TFG_CARMA3).  The shipping OU-II and
// OU-III filters keep their own model; docs/carma-imu-lever-arm.md compares
// the two.
//
// A rigid IMU at body offset r from the centre of gravity measures
//
//     f_imu = f_cog + alpha x r + omega x (omega x r).
//
// Model::BandLimited (the default) is the deployable model of the off-CoG
// studies (src/util/W3dSimCommon.h, W3dLeverArm): the bias-corrected rate
// passes a critically damped two-pole low-pass at cutoff_hz (15 Hz), alpha is
// the causal second-order backward difference of what is left, and both terms
// are evaluated on the filtered rate.  The low-pass is what makes the
// derivative usable: a differentiated gyro carries noise rising with
// frequency, and only a second-order roll-off band-limits it.  A constant
// gyro-bias error passes the low-pass at unit gain and the difference at zero,
// so only the centripetal term depends on the bias:
//
//     d(lever)/d(b_g) = -( (w . r) I + w r^T - 2 r w^T ),   w the filtered rate.
//
// Model::OuIII is Kalman3D_Wave_OU_III's set_imu_lever_arm_body() term for
// term (without its wind-heel frame), kept for that comparison: alpha is the
// first difference of the unfiltered rate through a one-pole low-pass of time
// constant tau_alpha, and its bias Jacobian carries d(alpha)/d(b_g) =
// -k_alpha I, k_alpha = (1 - exp(-dt/tau_alpha)) / dt.  On the dataset both
// choices cost attitude accuracy (the study's numbers); it is not the default.
//
// A zero offset (the default) disables the term; the rate history is kept
// either way, so the correction can be switched on at any time.

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace ocean_imu::kalman::filter_detail {

template<typename T>
class ImuLeverArm {
public:
    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;

    enum class Model : uint8_t {
        BandLimited,  // two-pole low-pass + second-order difference (default)
        OuIII         // OU-III's one-pole first difference and bias Jacobian
    };

    // r_b: IMU position relative to the CoG in the body frame, metres.  A
    // non-finite offset is ignored.
    void set_offset_body(const Vector3& r_b) {
        if (!r_b.allFinite()) return;
        r_ = r_b;
        enabled_ = r_b.squaredNorm() > T(0);
    }
    void clear() {
        r_.setZero();
        enabled_ = false;
    }
    [[nodiscard]] bool enabled() const noexcept { return enabled_; }
    [[nodiscard]] const Vector3& offset_body() const noexcept { return r_; }

    void set_model(Model m) {
        if (m == model_) return;
        model_ = m;
        reset_kinematics();
    }
    [[nodiscard]] Model model() const noexcept { return model_; }

    // BandLimited: corner of the two-pole low-pass ahead of the difference.
    void set_derivative_cutoff_hz(T hz) {
        if (std::isfinite(hz) && hz > T(0)) cutoff_hz_ = hz;
    }
    [[nodiscard]] T derivative_cutoff_hz() const noexcept { return cutoff_hz_; }

    // OuIII: time constant of the one-pole smoothing of alpha (0 = raw).
    void set_alpha_smoothing_tau(T tau_sec) {
        if (std::isfinite(tau_sec)) tau_alpha_ = std::max(T(0), tau_sec);
    }
    [[nodiscard]] T alpha_smoothing_tau() const noexcept { return tau_alpha_; }

    // Forget the rate history, e.g. across a gap in the prediction steps.
    void reset_kinematics() {
        omega_.setZero();
        lp1_.setZero();
        prev_.setZero();
        prev2_.setZero();
        alpha_.setZero();
        seen_ = 0;
    }

    // Once per prediction step, with the bias-corrected body rate.
    void update_rate(const Vector3& omega, T dt) {
        if (model_ == Model::OuIII) {
            update_ou3_(omega, dt);
            return;
        }
        if (!(dt > T(0)) || !std::isfinite(dt)) return;
        if (seen_ == 0) {
            lp1_ = omega;
            omega_ = omega;
        } else {
            const T tau = T(1) / (T(2) * T(3.14159265358979323846) * cutoff_hz_);
            const T g = std::min(T(1), dt / (tau + dt));
            lp1_ += g * (omega - lp1_);
            omega_ += g * (lp1_ - omega_);
        }
        // Causal second-order backward difference on the actual step lengths
        // (exact for a quadratic); first-order for the second sample.
        if (seen_ >= 2) {
            const T h1 = dt_prev_, h2 = dt;
            alpha_ = omega_ * ((T(2) * h2 + h1) / (h2 * (h1 + h2)))
                   - prev_ * ((h1 + h2) / (h1 * h2))
                   + prev2_ * (h2 / (h1 * (h1 + h2)));
        } else if (seen_ == 1) {
            alpha_ = (omega_ - prev_) / dt;
        } else {
            alpha_.setZero();
        }
        prev2_ = prev_;
        prev_ = omega_;
        dt_prev_ = dt;
        if (seen_ < 2) ++seen_;
    }

    [[nodiscard]] const Vector3& angular_acceleration() const noexcept { return alpha_; }
    // The rate the terms are evaluated on.
    [[nodiscard]] const Vector3& angular_rate() const noexcept { return omega_; }

    // alpha x r + omega x (omega x r).
    [[nodiscard]] Vector3 acceleration() const {
        return alpha_.cross(r_) + omega_.cross(omega_.cross(r_));
    }

    // d(acceleration)/d(b_g), body frame.
    [[nodiscard]] Matrix3 bias_jacobian() const {
        const Matrix3 d_centripetal =
            Matrix3::Identity() * omega_.dot(r_) + omega_ * r_.transpose() - T(2) * (r_ * omega_.transpose());
        if (model_ != Model::OuIII) return -d_centripetal;
        T k_alpha = T(0);
        if (seen_ > 0 && dt_prev_ > T(0)) {
            k_alpha = (tau_alpha_ > T(0))
                ? (T(1) - std::exp(-dt_prev_ / tau_alpha_)) / dt_prev_
                : T(1) / dt_prev_;
        }
        Matrix3 rx;
        rx <<  T(0),  -r_(2),  r_(1),
               r_(2),  T(0),  -r_(0),
              -r_(1),  r_(0),  T(0);
        return k_alpha * rx - d_centripetal;
    }

private:
    // Kalman3D_Wave_OU_III::time_update's alpha, line for line.
    void update_ou3_(const Vector3& omega, T dt) {
        dt_prev_ = dt;
        if (seen_ > 0 && dt > T(0)) {
            const Vector3 alpha_raw = (omega - omega_) / dt;
            if (tau_alpha_ > T(0)) {
                const T a = T(1) - std::exp(-dt / tau_alpha_);
                alpha_ = (T(1) - a) * alpha_ + a * alpha_raw;
            } else {
                alpha_ = alpha_raw;
            }
        } else {
            alpha_.setZero();
            seen_ = 1;
        }
        omega_ = omega;
    }

    Vector3 r_ = Vector3::Zero();
    bool enabled_ = false;
    Model model_ = Model::BandLimited;
    T cutoff_hz_ = T(15);
    T tau_alpha_ = T(0.05);  // OU-III's default

    Vector3 omega_ = Vector3::Zero();  // filtered rate (OuIII: last raw rate)
    Vector3 lp1_ = Vector3::Zero();
    Vector3 prev_ = Vector3::Zero();
    Vector3 prev2_ = Vector3::Zero();
    Vector3 alpha_ = Vector3::Zero();
    T dt_prev_ = T(1.0 / 200);
    int seen_ = 0;
};

}  // namespace ocean_imu::kalman::filter_detail
