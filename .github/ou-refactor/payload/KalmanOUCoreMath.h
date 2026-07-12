#pragma once

/*
  Copyright (c) 2025-2026 Mikhail Grushinskiy

  Shared OU and MEKF math used by the OU-II and OU-III filters.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>

namespace ocean_imu::kalman::ou_detail {

template<typename T>
inline T safe_inv_tau(T tau) {
    return T(1) / ((std::abs(tau) >= T(1e-8)) ? tau : std::copysign(T(1e-8), tau));
}

template<typename T>
struct OUPrims {
    T alpha;
    T em1;
};

template<typename T>
inline OUPrims<T> make_prims(T h, T tau) {
    const T x = h * safe_inv_tau(tau);
    return {std::exp(-x), std::expm1(-x)};
}

template<typename T>
inline Eigen::Quaternion<T> quat_from_delta_theta(const Eigen::Matrix<T,3,1>& dtheta) {
    const T theta = dtheta.norm();
    const T half_theta = T(0.5) * theta;

    T w;
    T k;
    if (theta < T(1e-2)) {
        const T t2 = theta * theta;
        const T t4 = t2 * t2;
        w = std::fma(-t2, T(1) / T(8), T(1));
        w = std::fma(t4, T(1) / T(384), w);
        k = std::fma(-t2, T(1) / T(48), T(0.5));
        k = std::fma(t4, T(1) / T(3840), k);
    } else {
        w = std::cos(half_theta);
        k = std::sin(half_theta) / theta;
    }

    const Eigen::Matrix<T,3,1> v = k * dtheta;
    Eigen::Quaternion<T> q(w, v.x(), v.y(), v.z());
    q.normalize();
    return q;
}

template<typename T, int N, bool enable_large_ldlt = false>
inline void project_psd(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    S = T(0.5) * (S + S.transpose());
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (!std::isfinite(S(i,j))) {
                S(i,j) = (i == j) ? eps : T(0);
            }
        }
    }

    if constexpr (N <= 4) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,N,N>> es(S);
        if (es.info() != Eigen::Success) {
            S.diagonal().array() += eps;
            S = T(0.5) * (S + S.transpose());
            return;
        }
        Eigen::Matrix<T,N,1> lam = es.eigenvalues();
        for (int i = 0; i < N; ++i) {
            if (!(lam(i) > T(0))) {
                lam(i) = eps;
            }
        }
        S = es.eigenvectors() * lam.asDiagonal() * es.eigenvectors().transpose();
    } else if constexpr (N <= 6 || enable_large_ldlt) {
        Eigen::LDLT<Eigen::Matrix<T,N,N>> ldlt;
        ldlt.compute(S);
        if (ldlt.info() != Eigen::Success) {
            T min_lb = std::numeric_limits<T>::infinity();
            for (int i = 0; i < N; ++i) {
                T row_sum = T(0);
                for (int j = 0; j < N; ++j) {
                    if (j != i) {
                        row_sum += std::abs(S(i,j));
                    }
                }
                min_lb = std::min(min_lb, S(i,i) - row_sum);
            }
            if (!(min_lb > eps)) {
                S.diagonal().array() += (eps - min_lb);
            }
            ldlt.compute(S);
            if (ldlt.info() != Eigen::Success) {
                S.diagonal().array() += T(10) * eps;
            }
        }
    }

    S = T(0.5) * (S + S.transpose());
}

template<typename T>
inline bool safe_ldlt3(Eigen::Matrix<T,3,3>& S,
                       Eigen::LDLT<Eigen::Matrix<T,3,3>>& ldlt,
                       T noise_scale) {
    ldlt.compute(S);
    if (ldlt.info() == Eigen::Success) {
        return true;
    }
    const T bump = std::max(std::numeric_limits<T>::epsilon(),
                            T(1e-6) * (noise_scale + T(1)));
    S.diagonal().array() += bump;
    ldlt.compute(S);
    return ldlt.info() == Eigen::Success;
}

template<typename T>
inline T nis3_from_ldlt(const Eigen::LDLT<Eigen::Matrix<T,3,3>>& ldlt,
                        const Eigen::Matrix<T,3,1>& residual) {
    const Eigen::Matrix<T,3,1> x = ldlt.solve(residual);
    if (!x.allFinite()) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    const T value = residual.dot(x);
    return std::isfinite(value) ? value : std::numeric_limits<T>::quiet_NaN();
}

template<typename T, int NX>
inline void gain_from_ldlt3(const Eigen::Matrix<T,NX,3>& PCt,
                            const Eigen::LDLT<Eigen::Matrix<T,3,3>>& ldlt,
                            Eigen::Matrix<T,NX,3>& K) {
    for (int i = 0; i < NX; ++i) {
        const Eigen::Matrix<T,3,1> rhs(PCt(i,0), PCt(i,1), PCt(i,2));
        const Eigen::Matrix<T,3,1> solution = ldlt.solve(rhs);
        K(i,0) = solution(0);
        K(i,1) = solution(1);
        K(i,2) = solution(2);
    }
}

template<typename T, int NX>
inline void joseph_update3(Eigen::Matrix<T,NX,NX>& P,
                           const Eigen::Matrix<T,NX,3>& K,
                           const Eigen::Matrix<T,3,3>& S,
                           const Eigen::Matrix<T,NX,3>& PCt) {
    for (int i = 0; i < NX; ++i) {
        for (int j = i; j < NX; ++j) {
            T KCP_ij = T(0);
            T KCP_ji = T(0);
            for (int l = 0; l < 3; ++l) {
                KCP_ij += K(i,l) * PCt(j,l);
                if (j != i) {
                    KCP_ji += K(j,l) * PCt(i,l);
                }
            }
            if (j == i) {
                KCP_ji = KCP_ij;
            }

            T KSK_ij = T(0);
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    KSK_ij += K(i,a) * S(a,b) * K(j,b);
                }
            }

            P(i,j) += -(KCP_ij + KCP_ji) + KSK_ij;
            if (j != i) {
                P(j,i) = P(i,j);
            }
        }
    }

    for (int i = 0; i < NX; ++i) {
        for (int j = i + 1; j < NX; ++j) {
            const T value = T(0.5) * (P(i,j) + P(j,i));
            P(i,j) = value;
            P(j,i) = value;
        }
    }
}

template<typename T>
inline Eigen::Matrix<T,3,3> skew(const Eigen::Matrix<T,3,1>& v) {
    Eigen::Matrix<T,3,3> M;
    M << T(0), -v.z(), v.y(),
         v.z(), T(0), -v.x(),
         -v.y(), v.x(), T(0);
    return M;
}

template<typename T>
inline void rot_and_B_from_wt(const Eigen::Matrix<T,3,1>& w,
                              T t,
                              Eigen::Matrix<T,3,3>& R,
                              Eigen::Matrix<T,3,3>& B) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const T wnorm = w.norm();
    const Matrix3 W = skew(w);

    if (wnorm < T(1e-7)) {
        const T t2 = t * t;
        const T t3 = t2 * t;
        R = Matrix3::Identity() - W * t + T(0.5) * (W * W) * t2;
        B = Matrix3::Identity() * t - T(0.5) * W * t2 + (W * W) * (t3 / T(6));
        return;
    }

    const T theta = wnorm * t;
    const T s = std::sin(theta);
    const T c = std::cos(theta);
    const T invw = T(1) / wnorm;
    const T invw2 = invw * invw;
    const Matrix3 K = W * invw;

    R = Matrix3::Identity() - s * K + (T(1) - c) * (K * K);
    B = Matrix3::Identity() * t
      - ((T(1) - c) * invw2) * W
      + ((t - s * invw) * invw2) * (W * W);
}

template<typename T>
inline void integral_B_ds(const Eigen::Matrix<T,3,1>& w,
                          T step,
                          Eigen::Matrix<T,3,3>& IB) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const T wnorm = w.norm();
    const Matrix3 W = skew(w);

    if (wnorm < T(1e-7)) {
        const T t2 = step * step;
        const T t3 = t2 * step;
        const T t4 = t3 * step;
        IB = Matrix3::Identity() * (T(0.5) * t2)
           - W * (T(1) / T(6) * t3)
           + (W * W) * (T(1) / T(24) * t4);
        return;
    }

    const T theta = wnorm * step;
    const T s = std::sin(theta);
    const T c = std::cos(theta);
    const T invw = T(1) / wnorm;
    const T invw2 = invw * invw;

    IB = Matrix3::Identity() * (T(0.5) * step * step)
       - ((step - s * invw) * invw2) * W
       + ((T(0.5) * step * step + (c - T(1)) * invw2) * invw2) * (W * W);
}

template<typename T>
inline Eigen::Matrix<T,3,3> d_omega_x_omega_x_r_domega(
    const Eigen::Matrix<T,3,1>& w,
    const Eigen::Matrix<T,3,1>& r) {
    return Eigen::Matrix<T,3,3>::Identity() * w.dot(r)
         + w * r.transpose()
         - T(2) * r * w.transpose();
}

template<typename T>
inline Eigen::Matrix<T,3,3> simpson_R_Q_RT(const Eigen::Matrix<T,3,1>& w,
                                            T step,
                                            const Eigen::Matrix<T,3,3>& Q) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 R0;
    Matrix3 Rm;
    Matrix3 R1;
    Matrix3 B;
    rot_and_B_from_wt(w, T(0), R0, B);
    rot_and_B_from_wt(w, T(0.5) * step, Rm, B);
    rot_and_B_from_wt(w, step, R1, B);
    const Matrix3 f0 = R0 * Q * R0.transpose();
    const Matrix3 f1 = Rm * Q * Rm.transpose();
    const Matrix3 f2 = R1 * Q * R1.transpose();
    return (step / T(6)) * (f0 + T(4) * f1 + f2);
}

template<typename T>
inline Eigen::Matrix<T,3,3> simpson_B_Q_BT(const Eigen::Matrix<T,3,1>& w,
                                            T step,
                                            const Eigen::Matrix<T,3,3>& Q) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 R;
    Matrix3 B0;
    Matrix3 Bm;
    Matrix3 B1;
    rot_and_B_from_wt(w, T(0), R, B0);
    rot_and_B_from_wt(w, T(0.5) * step, R, Bm);
    rot_and_B_from_wt(w, step, R, B1);
    const Matrix3 g0 = B0 * Q * B0.transpose();
    const Matrix3 g1 = Bm * Q * Bm.transpose();
    const Matrix3 g2 = B1 * Q * B1.transpose();
    return (step / T(6)) * (g0 + T(4) * g1 + g2);
}

template<typename T>
inline bool is_isotropic3(const Eigen::Matrix<T,3,3>& S, T tol = T(1e-9)) {
    Eigen::Matrix<T,3,3> off = S;
    off.diagonal().setZero();
    const T mean = S.trace() / T(3);
    return std::abs(S(0,0) - mean)
         + std::abs(S(1,1) - mean)
         + std::abs(S(2,2) - mean)
         + off.cwiseAbs().sum()
         <= tol * (T(1) + std::abs(mean));
}

template<typename T, int Integrals>
struct IntegratedOUChain;

template<typename T>
struct IntegratedOUChain<T, 2> {
    using MatrixAxis = Eigen::Matrix<T,3,3>;

    static void transition(T tau, T h, MatrixAxis& phi) {
        const auto prims = make_prims<T>(h, tau);
        const T x = h * safe_inv_tau(tau);
        const T x2 = x * x;
        const T x3 = x2 * x;
        const T x4 = x3 * x;
        const T tau2 = tau * tau;
        const T phi_pa = std::abs(x) < T(1e-2)
            ? tau2 * (T(0.5) * x2 - T(1) / T(6) * x3 + T(1) / T(24) * x4)
            : tau2 * (x + prims.em1);

        phi.setZero();
        phi(0,0) = T(1);
        phi(0,2) = -tau * prims.em1;
        phi(1,0) = h;
        phi(1,1) = T(1);
        phi(1,2) = phi_pa;
        phi(2,2) = std::min(prims.alpha, T(1));
    }

    static void process_covariance(T tau, T h, T sigma2, MatrixAxis& qd) {
        const T tau_eff = std::max(tau, T(1e-7));
        const T inv_tau = T(1) / tau_eff;
        const T x = h * inv_tau;

        if (std::abs(x) < T(1e-2)) {
            const T inv_tau2 = inv_tau * inv_tau;
            const T inv_tau3 = inv_tau2 * inv_tau;
            const T inv_tau4 = inv_tau3 * inv_tau;
            const T h2 = h * h;
            const T h3 = h2 * h;
            const T h4 = h3 * h;
            const T h5 = h4 * h;

            qd.setZero();
            qd(0,0) = sigma2 * ((T(2) / T(3)) * h3 * inv_tau
                              - (T(1) / T(2)) * h4 * inv_tau2
                              + (T(7) / T(30)) * h5 * inv_tau3);
            qd(0,1) = sigma2 * ((T(1) / T(4)) * h4 * inv_tau
                              - (T(1) / T(6)) * h5 * inv_tau2);
            qd(0,2) = sigma2 * (h2 * inv_tau
                              - h3 * inv_tau2
                              + (T(7) / T(12)) * h4 * inv_tau3
                              - (T(1) / T(4)) * h5 * inv_tau4);
            qd(1,0) = qd(0,1);
            qd(1,1) = sigma2 * ((T(1) / T(10)) * h5 * inv_tau);
            qd(1,2) = sigma2 * ((T(1) / T(3)) * h3 * inv_tau
                              - (T(1) / T(3)) * h4 * inv_tau2
                              + (T(11) / T(60)) * h5 * inv_tau3);
            qd(2,0) = qd(0,2);
            qd(2,1) = qd(1,2);
            qd(2,2) = sigma2 * (T(2) * h * inv_tau
                              - T(2) * h2 * inv_tau2
                              + (T(4) / T(3)) * h3 * inv_tau3
                              - (T(2) / T(3)) * h4 * inv_tau4);
            project_psd<T,3>(qd, T(1e-12));
            qd = T(0.5) * (qd + qd.transpose());
            return;
        }

        const T alpha = std::exp(-x);
        const T alpha2 = alpha * alpha;
        const T q_c = T(2) * sigma2 * inv_tau;
        const T tau2 = tau_eff * tau_eff;
        const T tau3 = tau2 * tau_eff;
        const T tau4 = tau3 * tau_eff;
        const T tau5 = tau4 * tau_eff;
        const T x2 = x * x;
        const T x3 = x2 * x;

        const T K00 = tau3 * (-alpha2 + T(4) * alpha + T(2) * x - T(3)) / T(2);
        const T K01 = tau4 * (alpha2 + T(2) * alpha * (x - T(1)) + x2 - T(2) * x + T(1)) / T(2);
        const T K02 = tau2 * (alpha2 - T(2) * alpha + T(1)) / T(2);
        const T K11 = tau5 * (-alpha2 / T(2) - T(2) * alpha * x + x3 / T(3) - x2 + x + T(1) / T(2));
        const T K12 = tau3 * (-alpha2 - T(2) * alpha * x + T(1)) / T(2);
        const T K22 = tau_eff * (T(1) - alpha2) / T(2);

        qd.setZero();
        qd(0,0) = q_c * K00;
        qd(0,1) = q_c * K01;
        qd(0,2) = q_c * K02;
        qd(1,0) = qd(0,1);
        qd(1,1) = q_c * K11;
        qd(1,2) = q_c * K12;
        qd(2,0) = qd(0,2);
        qd(2,1) = qd(1,2);
        qd(2,2) = q_c * K22;
        project_psd<T,3>(qd, T(1e-12));
        qd = T(0.5) * (qd + qd.transpose());
    }
};

template<typename T>
struct IntegratedOUChain<T, 3> {
    using MatrixAxis = Eigen::Matrix<T,4,4>;

    static void transition(T tau, T h, MatrixAxis& phi) {
        const auto prims = make_prims<T>(h, tau);
        const T x = h * safe_inv_tau(tau);
        const T x2 = x * x;
        const T x3 = x2 * x;
        const T x4 = x3 * x;
        const T x5 = x4 * x;
        const T tau2 = tau * tau;
        const T tau3 = tau2 * tau;
        const T phi_pa = std::abs(x) < T(1e-2)
            ? tau2 * (T(0.5) * x2 - T(1) / T(6) * x3 + T(1) / T(24) * x4)
            : tau2 * (x + prims.em1);
        const T phi_Sa = std::abs(x) < T(1e-2)
            ? tau3 * (T(1) / T(6) * x3 - T(1) / T(24) * x4 + T(1) / T(120) * x5)
            : tau3 * (T(0.5) * x2 - x - prims.em1);

        phi.setZero();
        phi(0,0) = T(1);
        phi(0,3) = -tau * prims.em1;
        phi(1,0) = h;
        phi(1,1) = T(1);
        phi(1,3) = phi_pa;
        phi(2,0) = T(0.5) * h * h;
        phi(2,1) = h;
        phi(2,2) = T(1);
        phi(2,3) = phi_Sa;
        phi(3,3) = std::min(prims.alpha, T(1));
    }

    static void process_covariance(T tau, T h, T sigma2, MatrixAxis& qd) {
        const T tau_eff = std::max(tau, T(1e-7));
        const T inv_tau = T(1) / tau_eff;
        const T x = h * inv_tau;

        if (std::abs(x) < T(1e-2)) {
            const T inv_tau2 = inv_tau * inv_tau;
            const T inv_tau3 = inv_tau2 * inv_tau;
            const T inv_tau4 = inv_tau3 * inv_tau;
            const T inv_tau5 = inv_tau4 * inv_tau;
            const T h2 = h * h;
            const T h3 = h2 * h;
            const T h4 = h3 * h;
            const T h5 = h4 * h;

            qd.setZero();
            qd(0,0) = sigma2 * ((T(2) / T(3)) * h3 * inv_tau
                              - T(1) / T(2) * h4 * inv_tau2
                              + T(7) / T(30) * h5 * inv_tau3);
            qd(0,1) = sigma2 * (T(1) / T(4) * h4 * inv_tau
                              - T(1) / T(6) * h5 * inv_tau2);
            qd(0,2) = sigma2 * (T(1) / T(15) * h5 * inv_tau);
            qd(0,3) = sigma2 * (h2 * inv_tau
                              - h3 * inv_tau2
                              + T(7) / T(12) * h4 * inv_tau3
                              - T(1) / T(4) * h5 * inv_tau4);
            qd(1,0) = qd(0,1);
            qd(1,1) = sigma2 * (T(1) / T(10) * h5 * inv_tau);
            qd(1,2) = T(0);
            qd(1,3) = sigma2 * (T(1) / T(3) * h3 * inv_tau
                              - T(1) / T(3) * h4 * inv_tau2
                              + T(11) / T(60) * h5 * inv_tau3);
            qd(2,0) = qd(0,2);
            qd(2,1) = qd(1,2);
            qd(2,2) = T(0);
            qd(2,3) = sigma2 * (T(1) / T(12) * h4 * inv_tau
                              - T(1) / T(12) * h5 * inv_tau2);
            qd(3,0) = qd(0,3);
            qd(3,1) = qd(1,3);
            qd(3,2) = qd(2,3);
            qd(3,3) = sigma2 * (T(2) * h * inv_tau
                              - T(2) * h2 * inv_tau2
                              + T(4) / T(3) * h3 * inv_tau3
                              - T(2) / T(3) * h4 * inv_tau4
                              + T(4) / T(15) * h5 * inv_tau5);

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    if (!std::isfinite(qd(i,j))) {
                        qd(i,j) = (i == j) ? T(1e-18) : T(0);
                    }
                }
            }
            project_psd<T,4>(qd, T(1e-12));
            qd = T(0.5) * (qd + qd.transpose());
            return;
        }

        const T alpha = std::exp(-x);
        const T alpha2 = alpha * alpha;
        const T q_c = T(2) * sigma2 * inv_tau;
        const T tau2 = tau_eff * tau_eff;
        const T tau3 = tau2 * tau_eff;
        const T tau4 = tau3 * tau_eff;
        const T tau5 = tau4 * tau_eff;
        const T tau6 = tau5 * tau_eff;
        const T tau7 = tau6 * tau_eff;
        const T x2 = x * x;
        const T x3 = x2 * x;
        const T x4 = x3 * x;
        const T x5 = x4 * x;

        const T K00 = tau3 * (-alpha2 + T(4) * alpha + T(2) * x - T(3)) / T(2);
        const T K01 = tau4 * (alpha2 + T(2) * alpha * (x - T(1)) + x2 - T(2) * x + T(1)) / T(2);
        const T K02 = tau5 * (-T(3) * alpha2 + T(3) * alpha * (x2 + T(4)) + x3 - T(3) * x2 + T(6) * x - T(9)) / T(6);
        const T K03 = tau2 * (alpha2 - T(2) * alpha + T(1)) / T(2);
        const T K11 = tau5 * (-alpha2 / T(2) - T(2) * alpha * x + x3 / T(3) - x2 + x + T(1) / T(2));
        const T K12 = tau6 * (alpha2 / T(2)
                            + alpha * (-x2 + T(2) * x - T(2)) / T(2)
                            + x4 / T(8) - x3 / T(2) + x2 - x + T(1) / T(2));
        const T K13 = tau3 * (-alpha2 - T(2) * alpha * x + T(1)) / T(2);
        const T K22 = tau7 * (-alpha2 / T(2)
                            + alpha * x2 + T(2) * alpha
                            + x5 / T(20) - x4 / T(4) + T(2) * x3 / T(3)
                            - x2 + x - T(3) / T(2));
        const T K23 = tau4 * (alpha2 - alpha * (x2 + T(2)) + T(1)) / T(2);
        const T K33 = tau_eff * (T(1) - alpha2) / T(2);

        qd.setZero();
        qd(0,0) = q_c * K00;
        qd(0,1) = q_c * K01;
        qd(0,2) = q_c * K02;
        qd(0,3) = q_c * K03;
        qd(1,0) = qd(0,1);
        qd(1,1) = q_c * K11;
        qd(1,2) = q_c * K12;
        qd(1,3) = q_c * K13;
        qd(2,0) = qd(0,2);
        qd(2,1) = qd(1,2);
        qd(2,2) = q_c * K22;
        qd(2,3) = q_c * K23;
        qd(3,0) = qd(0,3);
        qd(3,1) = qd(1,3);
        qd(3,2) = qd(2,3);
        qd(3,3) = q_c * K33;

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (!std::isfinite(qd(i,j))) {
                    qd(i,j) = (i == j) ? T(1e-18) : T(0);
                }
            }
        }
        project_psd<T,4>(qd, T(1e-12));
        qd = T(0.5) * (qd + qd.transpose());
    }
};

} // namespace ocean_imu::kalman::ou_detail
