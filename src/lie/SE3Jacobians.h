#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#include "lie/SO3Jacobians.h"
namespace ocean_imu::lie {
/*
    Off-diagonal (vector <- attitude) block of the right Jacobian of the
    two-frame group.  Every vector block (world column or world-referred body
    bias) has ad_xi = [[ [phi]x, 0 ], [ [z]x, [phi]x ]], so

        J_r = sum_n (-ad)^n / (n+1)!  =  [[ J_l(-phi), 0 ], [ Q_r(phi, z), J_l(-phi) ]],
        Q_r(phi, z) = Q_l(-phi, -z),

    Q_l being the SE(3) left-Jacobian block (Barfoot, State Estimation for
    Robotics, eq. 7.86).  Its scalar coefficients are evaluated from their
    alternating series below |phi| = 1, where the closed forms cancel in
    float.
*/
/*
    Q_r(phi, .) is linear in z and its matrix products collapse: with
    P = [-phi]x, Z = [-z]x, s = phi . z, t2 = |phi|^2 and
    [a]x [b]x = b a^T - (a . b) I,

        P Z + Z P = z phi^T + phi z^T - 2 s I,   P Z P = -s P,
        P P Z + Z P P = 2 s P - t2 Z,            P Z P P + P P Z P = -2 s P^2,

    so

        Q_r = -(1/2 - c2 t2) [z]x + c1 (z phi^T + phi z^T) + s N,
        N   = -2 c1 I - (2 c2 - c1) [phi]x - 2 c3 (phi phi^T - t2 I).

    The coefficients and N depend on phi alone: one reset evaluates them once
    for all of its vector blocks.  Plain arrays, as the sample path below.
*/
template<typename T>
class RightJacobianCross {
public:
    explicit RightJacobianCross(const Eigen::Matrix<T,3,1>& phi) {
        const T x = phi(0), y = phi(1), z = phi(2);
        const T t2 = x * x + y * y + z * z;
        T c1, c2, c3;
        coefficients_(t2, c1, c2, c3);
        phi_[0] = x; phi_[1] = y; phi_[2] = z;
        k0_ = T(0.5) - c2 * t2;
        k1_ = c1;
        const T k2 = T(2) * c2 - c1, k3 = T(2) * c3;
        const T d = -T(2) * c1 + k3 * t2;
        n_[0][0] = d - k3 * x * x;       n_[0][1] = k2 * z - k3 * x * y;  n_[0][2] = -k2 * y - k3 * x * z;
        n_[1][0] = -k2 * z - k3 * x * y; n_[1][1] = d - k3 * y * y;       n_[1][2] = k2 * x - k3 * y * z;
        n_[2][0] = k2 * y - k3 * x * z;  n_[2][1] = -k2 * x - k3 * y * z; n_[2][2] = d - k3 * z * z;
    }
    // out = Q_r(phi, z), row-major.
    void operator()(const T* z, T (&out)[3][3]) const {
        const T s = phi_[0] * z[0] + phi_[1] * z[1] + phi_[2] * z[2];
        const T y[3] = {k1_ * z[0], k1_ * z[1], k1_ * z[2]};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out[i][j] = y[i] * phi_[j] + phi_[i] * y[j] + s * n_[i][j];
        const T h[3] = {k0_ * z[0], k0_ * z[1], k0_ * z[2]};
        out[0][1] += h[2]; out[0][2] -= h[1];
        out[1][0] -= h[2]; out[1][2] += h[0];
        out[2][0] += h[1]; out[2][1] -= h[0];
    }

private:
    static void coefficients_(T t2, T& c1, T& c2, T& c3);
    T phi_[3], n_[3][3], k0_, k1_;
};

template<typename T>
void RightJacobianCross<T>::coefficients_(T t2, T& c1, T& c2, T& c3) {
    if (t2 < T(1)) {
        // c1 = sum (-1)^k t^2k/(2k+3)!, c2 = sum (-1)^k t^2k/(2k+4)!,
        // c3 = sum (-1)^k (k+1) t^2k/(2k+5)!; nine terms leave < 1e-19.
        static constexpr double f3[9] = {1.0/6.0, 1.0/120.0, 1.0/5040.0, 1.0/362880.0,
            1.0/39916800.0, 1.0/6227020800.0, 1.0/1307674368000.0, 1.0/355687428096000.0,
            1.0/121645100408832000.0};
        static constexpr double f4[9] = {1.0/24.0, 1.0/720.0, 1.0/40320.0, 1.0/3628800.0,
            1.0/479001600.0, 1.0/87178291200.0, 1.0/20922789888000.0, 1.0/6402373705728000.0,
            1.0/2432902008176640000.0};
        static constexpr double f5[9] = {1.0/120.0, 2.0/5040.0, 3.0/362880.0, 4.0/39916800.0,
            5.0/6227020800.0, 6.0/1307674368000.0, 7.0/355687428096000.0,
            8.0/121645100408832000.0, 9.0/51090942171709440000.0};
        c1 = T(f3[8]); c2 = T(f4[8]); c3 = T(f5[8]);
        for (int k = 7; k >= 0; --k) {
            c1 = T(f3[k]) - t2 * c1;
            c2 = T(f4[k]) - t2 * c2;
            c3 = T(f5[k]) - t2 * c3;
        }
    } else {
        const T t = std::sqrt(t2);
        const T s = std::sin(t), c = std::cos(t);
        c1 = (t - s) / (t2 * t);
        c2 = (t2 + T(2) * c - T(2)) / (T(2) * t2 * t2);
        c3 = (T(2) * t - T(3) * s + t * c) / (T(2) * t2 * t2 * t);
    }
}

template<typename T>
inline Eigen::Matrix<T,3,3> right_jacobian_cross(const Eigen::Matrix<T,3,1>& phi,
                                                 const Eigen::Matrix<T,3,1>& z) {
    T q[3][3];
    RightJacobianCross<T>{phi}(z.data(), q);
    Eigen::Matrix<T,3,3> Q;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Q(i,j) = q[i][j];
    return Q;
}


} // namespace ocean_imu::lie
