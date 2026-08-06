#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    SO(3) exponential, logarithm, and left Jacobians for the right-invariant
    two-frame Lie-group filter (see src/lie/TwoFrameGroup.h).

    Everything is templated on the scalar type. Firmware instantiates float;
    the finite-difference and Lie++ cross-check tests instantiate double,
    where the tolerances mean something. There are no hardcoded float types
    and no literals that would force a promotion, so -Wconversion stays clean
    at both precisions.

    Small-angle branches reuse the thresholds and the fma-based expansion
    style already established in kalman_ou_common/KalmanOUCoreMath.h. Those
    thresholds are chosen for float and remain safe for double.

    CONVENTION. This file's Exp() is the standard right-handed exponential:
    Exp(phi) rotates by +phi about phi/|phi|, and the filter uses it to carry
    a body-to-world rotation R = R_bw. KalmanOUCoreMath's rot_and_B_from_wt()
    returns the *transpose* of this, because OU-III propagates a world-to-body
    quaternion. The two are not interchangeable; swapping them silently
    inverts every rotation in the filter.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>

#include "kalman_ou_common/KalmanOUCoreMath.h"

namespace ocean_imu::lie {

// skew() and quat_from_delta_theta() already exist in the OU core math and
// carry the expansions this file has to agree with. Alias rather than clone.
using ocean_imu::kalman::ou_detail::skew;
using ocean_imu::kalman::ou_detail::quat_from_delta_theta;

// Angle below which the trigonometric forms lose significance and the Taylor
// branches take over. Matches KalmanOUCoreMath.h's quat_from_delta_theta.
template<typename T>
inline constexpr T kSmallAngle = T(1e-2);

/*
    Exp: so(3) -> SO(3).

        R = I + (sin t / t) [phi]x + ((1 - cos t) / t^2) [phi]x^2,  t = |phi|
*/
template<typename T>
inline Eigen::Matrix<T,3,3> Exp(const Eigen::Matrix<T,3,1>& phi) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const T t2 = phi.squaredNorm();
    const T t = std::sqrt(t2);
    const Matrix3 W = skew<T>(phi);

    T a, b;  // coefficients of [phi]x and [phi]x^2
    if (t < kSmallAngle<T>) {
        const T t4 = t2 * t2;
        // sin(t)/t
        a = T(1);
        a = std::fma(-t2, T(1)/T(6), a);
        a = std::fma( t4, T(1)/T(120), a);
        // (1-cos(t))/t^2
        b = T(0.5);
        b = std::fma(-t2, T(1)/T(24), b);
        b = std::fma( t4, T(1)/T(720), b);
    } else {
        a = std::sin(t) / t;
        b = (T(1) - std::cos(t)) / t2;
    }
    return Matrix3::Identity() + a * W + b * (W * W);
}

/*
    Log: SO(3) -> so(3), the exact inverse of Exp above.

    Recovered through the quaternion rather than from acos(trace) so the
    near-identity case stays accurate and the near-pi case stays defined.
*/
template<typename T>
inline Eigen::Matrix<T,3,1> Log(const Eigen::Matrix<T,3,3>& R) {
    Eigen::Quaternion<T> q(R);
    q.normalize();
    if (q.w() < T(0)) q.coeffs() = -q.coeffs();  // shortest arc

    const Eigen::Matrix<T,3,1> v(q.x(), q.y(), q.z());
    const T sin_half = v.norm();
    if (sin_half < kSmallAngle<T> * T(0.5)) {
        // phi = 2 v (1 + s^2/6 + 3 s^4/40),  s = |v| = sin(t/2)
        const T s2 = sin_half * sin_half;
        T k = T(2);
        k = std::fma(s2, T(2)/T(6), k);
        k = std::fma(s2 * s2, T(6)/T(40), k);
        return k * v;
    }
    const T half_angle = std::atan2(sin_half, q.w());
    return v * (T(2) * half_angle / sin_half);
}

/*
    Left Jacobian of SO(3).

        J_l(phi) = I + ((1-cos t)/t^2) [phi]x + ((t - sin t)/t^3) [phi]x^2
                 = int_0^1 Exp(u phi) du

    This is the factor coupling an attitude correction into the world vectors
    of the group exponential: E_X = J_l(phi) [rho_v rho_p rho_S rho_a].
*/
template<typename T>
inline Eigen::Matrix<T,3,3> left_jacobian(const Eigen::Matrix<T,3,1>& phi) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const T t2 = phi.squaredNorm();
    const T t = std::sqrt(t2);
    const Matrix3 W = skew<T>(phi);

    T a, b;
    if (t < kSmallAngle<T>) {
        const T t4 = t2 * t2;
        // (1-cos t)/t^2
        a = T(0.5);
        a = std::fma(-t2, T(1)/T(24), a);
        a = std::fma( t4, T(1)/T(720), a);
        // (t - sin t)/t^3
        b = T(1)/T(6);
        b = std::fma(-t2, T(1)/T(120), b);
        b = std::fma( t4, T(1)/T(5040), b);
    } else {
        a = (T(1) - std::cos(t)) / t2;
        b = (t - std::sin(t)) / (t2 * t);
    }
    return Matrix3::Identity() + a * W + b * (W * W);
}

/*
    Inverse left Jacobian of SO(3).

        J_l^-1(phi) = I - (1/2)[phi]x + c(t) [phi]x^2
        c(t)        = (1 - (t/2) cot(t/2)) / t^2

    The cotangent form is singular at t = 2 pi; the filter never retracts by
    anything close to that, and the small-angle branch covers the other end.
*/
template<typename T>
inline Eigen::Matrix<T,3,3> inv_left_jacobian(const Eigen::Matrix<T,3,1>& phi) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const T t2 = phi.squaredNorm();
    const T t = std::sqrt(t2);
    const Matrix3 W = skew<T>(phi);

    T c;
    if (t < kSmallAngle<T>) {
        const T t4 = t2 * t2;
        c = T(1)/T(12);
        c = std::fma(t2, T(1)/T(720), c);
        c = std::fma(t4, T(1)/T(30240), c);
    } else {
        const T half = T(0.5) * t;
        c = (T(1) - half * std::cos(half) / std::sin(half)) / t2;
    }
    return Matrix3::Identity() - T(0.5) * W + c * (W * W);
}

// Right Jacobians. J_r(phi) = J_l(-phi) = Exp(-phi) J_l(phi); this identity is
// what the bias block of the group exponential rests on.
template<typename T>
inline Eigen::Matrix<T,3,3> right_jacobian(const Eigen::Matrix<T,3,1>& phi) {
    return left_jacobian<T>(Eigen::Matrix<T,3,1>(-phi));
}

template<typename T>
inline Eigen::Matrix<T,3,3> inv_right_jacobian(const Eigen::Matrix<T,3,1>& phi) {
    return inv_left_jacobian<T>(Eigen::Matrix<T,3,1>(-phi));
}

} // namespace ocean_imu::lie
