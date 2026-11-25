#pragma once

#include "lie/SO3Jacobians.h"

namespace ocean_imu::kalman::tfg_detail {

// Isotropic gyro and gyro-bias driving noise, with the bias error expressed
// in the final world frame. For |angle| <= 0.1 the closed-form integral is
// well conditioned in these series in theta^2. Keep the quadrature for larger
// rotations and non-isotropic driving noise.
template<typename T>
inline void isotropic_gyro_noise(T dt, T angle_squared,
                                 const Eigen::Matrix<T,3,1>& world_angle,
                                 T gyro_variance, T bias_variance,
                                 Eigen::Matrix<T,3,3>& Qpp, Eigen::Matrix<T,3,3>& Qpb,
                                 Eigen::Matrix<T,3,3>& Qbb) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const Matrix3 V = lie::skew<T>(world_angle);
    const Matrix3 V2 = V * V;
    const T x = angle_squared;
    // (theta-sin(theta))/theta^3, (cos(theta)-1+x/2)/theta^4 and (1/3-2*a)/x.
    const T a = T(1.0/6.0) + x*(T(-1.0/120.0) + x*(T(1.0/5040.0) - x*T(1.0/362880.0)));
    const T c = T(1.0/24.0) + x*(T(-1.0/720.0) + x*(T(1.0/40320.0) - x*T(1.0/3628800.0)));
    const T d = T(1.0/60.0) + x*(T(-1.0/2520.0) + x*(T(1.0/181440.0) - x*T(1.0/19958400.0)));
    const T bh = bias_variance * dt, bh2 = bh * dt, bh3 = bh2 * dt;
    Qbb = Matrix3::Identity() * bh;
    Qpb = bh2 * (a * V - c * V2);
    Qpb.diagonal().array() -= bh2 * T(0.5);
    Qpp = (bh3 * d) * V2;
    Qpp.diagonal().array() += dt * gyro_variance + bh3 / T(3);
}

} // namespace ocean_imu::kalman::tfg_detail
