#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Attitude error-state math shared by the wave MEKFs (OU-II, OU-III, TFG,
// CARMA3, TFG-CARMA3): the gyro-bias ball, the rotation-vector quaternion,
// the constant-rate SO(3) transition and its integrals, the Simpson
// process-noise integrals and the covariance transport of a multiplicative
// attitude reset.  KalmanOUCoreMath.h re-exports these names into ou_detail.

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <limits>

namespace ocean_imu::kalman::filter_detail {

// The estimated state is residual sensor bias, not vessel angular rate.
// This deliberately generous fixed ball is separate from the commissioned
// physical residual-bias bound (0.02 rad/s). See docs/ou-gyro-bias-projection.md.
inline constexpr double gyro_bias_radius_rad_s = 0.5;

template<typename T, class Derived>
inline void project_gyro_bias(Eigen::MatrixBase<Derived>& bias) {
    const T radius = T(gyro_bias_radius_rad_s);
    if (!bias.allFinite()) { bias.setZero(); return; }
    const T largest = bias.cwiseAbs().maxCoeff();
    if (largest == T(0)) return;

    // The ordinary path is bit-for-bit unchanged. Bound the rounding of the
    // three squared components before accepting a point close to the sphere.
    if (largest <= radius) {
        if (bias.squaredNorm() < radius*radius*(T(1)-T(8)*std::numeric_limits<T>::epsilon())) return;
        // Exact axial boundary points must remain unchanged.
        if (largest == radius && (bias.array() != T(0)).count() == 1) return;
        const double infinity = std::numeric_limits<double>::infinity();
        double upper = 0;
        for (int i = 0; i < 3; ++i) {
            const double x = static_cast<double>(bias(i));
            if (x != 0) {
                const double square = std::nextafter(x*x, infinity);
                upper = std::nextafter(upper+square, infinity);
            }
        }
        if (upper <= double(radius)*double(radius)) return;
    }

    // Scale first: even finite near-FLT_MAX/DBL_MAX vectors must retain their
    // direction rather than overflow their norm. The tiny inward rounding
    // allowance keeps the stored vector inside the Euclidean sphere.
    const Eigen::Matrix<T,3,1> unit = bias / largest;
    const T inward = radius*(T(1)-T(8)*std::numeric_limits<T>::epsilon());
    bias = unit * (inward / unit.norm());
    // Mean only. In particular, no covariance clipping, scaling or reset.
}

template<typename T>
inline Eigen::Matrix<T,3,3> skew(const Eigen::Matrix<T,3,1>& vec) {
    Eigen::Matrix<T,3,3> M;
    M << T(0), -vec(2), vec(1),
         vec(2), T(0), -vec(0),
        -vec(1), vec(0), T(0);
    return M;
}

template<typename T>
inline Eigen::Quaternion<T> quat_from_delta_theta(const Eigen::Matrix<T,3,1>& dtheta) {
    const T theta = dtheta.norm();
    const T half_theta = T(0.5) * theta;

    T w, k;
    if (theta < T(1e-2)) {
        const T t2 = theta * theta;
        const T t4 = t2 * t2;
        w = T(1);
        w = std::fma(-t2, T(1)/T(8), w);
        w = std::fma( t4, T(1)/T(384), w);
        k = T(0.5);
        k = std::fma(-t2, T(1)/T(48), k);
        k = std::fma( t4, T(1)/T(3840), k);
    } else {
        w = std::cos(half_theta);
        k = std::sin(half_theta) / theta;
    }

    const Eigen::Matrix<T,3,1> v = k * dtheta;
    Eigen::Quaternion<T> q(w, v.x(), v.y(), v.z());
    q.normalize();
    return q;
}

// Scalar coefficients of the constant-rate SO(3) integrals, written as
// functions of the dimensionless angle x = |omega|*t so that, with W = [omega]x,
//
//   R(t)          = I - t g0(x) W + t^2 g1(x) W^2
//   B(t)          = int_0^t R = I t - t^2 g1(x) W + t^3 g2(x) W^2
//   int_0^t B(s)  = I t^2/2 - t^3 g2(x) W + t^4 g3(x) W^2
//
// g0 = sin(x)/x, g1 = (1 - cos x)/x^2, g2 = (x - sin x)/x^3 and
// g3 = (x^2/2 - 1 + cos x)/x^4. The closed forms of g1 (as written), g2 and
// g3 cancel: their rounding error grows like 1/x^2 ulp, so float keeps no
// correct digit below x ~ 3e-4 whatever |omega| is. For |x| < 1 all four are
// evaluated from their alternating even Taylor series through x^18; the
// first omitted term, at most 1/21! ~ 2e-20 at |x| = 1, is far below half a
// double ulp of every coefficient. At |x| >= 1 the closed-form
// amplification is bounded by ~12 ulp (g3) for float and double alike.
template<typename T>
struct SO3IntegralCoeffs {
    T g0, g1, g2, g3;
};

template<typename T>
inline SO3IntegralCoeffs<T> so3_integral_coeffs(T x) {
    SO3IntegralCoeffs<T> c{};
    if (std::abs(x) < T(1)) {
        const T y = x * x;
        // sum_{k=0..9} (-1)^k y^k / (2k+n)!, Horner from the highest order.
        auto alternating = [y](const T* inv_fact) {
            T acc = inv_fact[9];
            for (int k = 8; k >= 0; --k) acc = inv_fact[k] - y * acc;
            return acc;
        };
        static constexpr T n1[10] = {
            T(1.0/1.0), T(1.0/6.0), T(1.0/120.0), T(1.0/5040.0),
            T(1.0/362880.0), T(1.0/39916800.0), T(1.0/6227020800.0), T(1.0/1307674368000.0),
            T(1.0/355687428096000.0), T(1.0/1.21645100408832e+17)};
        static constexpr T n2[10] = {
            T(1.0/2.0), T(1.0/24.0), T(1.0/720.0), T(1.0/40320.0),
            T(1.0/3628800.0), T(1.0/479001600.0), T(1.0/87178291200.0), T(1.0/20922789888000.0),
            T(1.0/6402373705728000.0), T(1.0/2.43290200817664e+18)};
        static constexpr T n3[10] = {
            T(1.0/6.0), T(1.0/120.0), T(1.0/5040.0), T(1.0/362880.0),
            T(1.0/39916800.0), T(1.0/6227020800.0), T(1.0/1307674368000.0), T(1.0/355687428096000.0),
            T(1.0/1.21645100408832e+17), T(1.0/5.109094217170944e+19)};
        static constexpr T n4[10] = {
            T(1.0/24.0), T(1.0/720.0), T(1.0/40320.0), T(1.0/3628800.0),
            T(1.0/479001600.0), T(1.0/87178291200.0), T(1.0/20922789888000.0), T(1.0/6402373705728000.0),
            T(1.0/2.43290200817664e+18), T(1.0/1.1240007277776077e+21)};
        c.g0 = alternating(n1);
        c.g1 = alternating(n2);
        c.g2 = alternating(n3);
        c.g3 = alternating(n4);
        return c;
    }
    const T s = std::sin(x);
    const T half_s = std::sin(T(0.5) * x);
    const T x2 = x * x;
    c.g0 = s / x;
    c.g1 = T(2) * half_s * half_s / x2;
    c.g2 = (x - s) / (x2 * x);
    c.g3 = (T(0.5) * x2 - T(1) + std::cos(x)) / (x2 * x2);
    return c;
}

template<typename T>
inline void rot_and_B_from_wt(const Eigen::Matrix<T,3,1>& w, T t,
                              Eigen::Matrix<T,3,3>& R, Eigen::Matrix<T,3,3>& B) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const Matrix3 W = skew(w);
    const Matrix3 W2 = W * W;
    const SO3IntegralCoeffs<T> c = so3_integral_coeffs(w.norm() * t);
    const T t2 = t * t;
    R = Matrix3::Identity() - (t*c.g0)*W + (t2*c.g1)*W2;
    B = Matrix3::Identity()*t - (t2*c.g1)*W + (t2*t*c.g2)*W2;
}

template<typename T>
inline void integral_B_ds(const Eigen::Matrix<T,3,1>& w, T step,
                          Eigen::Matrix<T,3,3>& IB) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const Matrix3 W = skew(w);
    const SO3IntegralCoeffs<T> c = so3_integral_coeffs(w.norm() * step);
    const T T2 = step*step, T3 = T2*step;
    IB = Matrix3::Identity()*(T(0.5)*T2)
       - (T3*c.g2)*W
       + (T3*step*c.g3)*(W*W);
}

template<typename T>
inline Eigen::Matrix<T,3,3> simpson_R_Q_RT(const Eigen::Matrix<T,3,1>& w, T step,
                                            const Eigen::Matrix<T,3,3>& Q) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 R0, Btmp, Rm, R1;
    rot_and_B_from_wt(w, T(0), R0, Btmp);
    rot_and_B_from_wt(w, T(0.5)*step, Rm, Btmp);
    rot_and_B_from_wt(w, step, R1, Btmp);
    return (step/T(6)) * (R0*Q*R0.transpose() + T(4)*Rm*Q*Rm.transpose() + R1*Q*R1.transpose());
}

template<typename T>
inline Eigen::Matrix<T,3,3> simpson_B_Q_BT(const Eigen::Matrix<T,3,1>& w, T step,
                                            const Eigen::Matrix<T,3,3>& Q) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 Rtmp, B0, Bm, B1;
    rot_and_B_from_wt(w, T(0), Rtmp, B0);
    rot_and_B_from_wt(w, T(0.5)*step, Rtmp, Bm);
    rot_and_B_from_wt(w, step, Rtmp, B1);
    return (step/T(6)) * (B0*Q*B0.transpose() + T(4)*Bm*Q*Bm.transpose() + B1*Q*B1.transpose());
}

template<typename T>
inline bool is_isotropic3(const Eigen::Matrix<T,3,3>& S, T tol = T(1e-9)) {
    const T a = S(0,0), b = S(1,1), c = S(2,2);
    Eigen::Matrix<T,3,3> off_matrix = S;
    off_matrix.diagonal().setZero();
    const T off = off_matrix.cwiseAbs().sum();
    const T mean = (a+b+c)/T(3);
    return (std::abs(a-mean)+std::abs(b-mean)+std::abs(c-mean)+off)
        <= tol*(T(1)+std::abs(mean));
}

template<typename T, int NX>
inline void apply_left_error_reset(Eigen::Matrix<T,NX,NX>& covariance,
                                   const Eigen::Matrix<T,3,1>& dtheta) {
    if (!dtheta.allFinite() || !(dtheta.squaredNorm() > T(0))) return;
    const Eigen::Matrix<T,3,3> G = Eigen::Matrix<T,3,3>::Identity() + T(0.5)*skew(dtheta);

    Eigen::Matrix<T,3,3> Paa_old;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) Paa_old(i,j) = covariance(i,j);
    }

    Eigen::Matrix<T,3,3> GP;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T sum = T(0);
            for (int k = 0; k < 3; ++k) sum += G(i,k)*Paa_old(k,j);
            GP(i,j) = sum;
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T sum = T(0);
            for (int k = 0; k < 3; ++k) sum += GP(i,k)*G(j,k);
            covariance(i,j) = sum;
        }
    }

    for (int col = 3; col < NX; ++col) {
        const T old0 = covariance(0,col);
        const T old1 = covariance(1,col);
        const T old2 = covariance(2,col);
        const T new0 = G(0,0)*old0 + G(0,1)*old1 + G(0,2)*old2;
        const T new1 = G(1,0)*old0 + G(1,1)*old1 + G(1,2)*old2;
        const T new2 = G(2,0)*old0 + G(2,1)*old1 + G(2,2)*old2;
        covariance(0,col)=new0; covariance(1,col)=new1; covariance(2,col)=new2;
        covariance(col,0)=new0; covariance(col,1)=new1; covariance(col,2)=new2;
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = i+1; j < 3; ++j) {
            const T v = T(0.5)*(covariance(i,j)+covariance(j,i));
            covariance(i,j)=v; covariance(j,i)=v;
        }
    }
}

} // namespace ocean_imu::kalman::filter_detail
