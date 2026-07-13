#pragma once

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

template<typename T>
struct OUDiscreteCoeffs {
    T phi_pa;
    T phi_Sa;
};

template<typename T>
inline OUDiscreteCoeffs<T> safe_phi_A_coeffs(T h, T tau) {
    OUDiscreteCoeffs<T> c{};
    const T x = h * safe_inv_tau(tau);
    const T tau2 = tau * tau;
    const T tau3 = tau2 * tau;

    if (std::abs(x) < T(1e-2)) {
        const T x2 = x*x;
        const T x3 = x2*x;
        const T x4 = x3*x;
        const T x5 = x4*x;
        c.phi_pa = tau2 * (T(0.5)*x2 - T(1.0/6.0)*x3 + T(1.0/24.0)*x4);
        c.phi_Sa = tau3 * (T(1.0/6.0)*x3 - T(1.0/24.0)*x4 + T(1.0/120.0)*x5);
    } else {
        const T em1 = std::expm1(-x);
        c.phi_pa = tau2 * (x + em1);
        c.phi_Sa = tau3 * (T(0.5)*x*x - x - em1);
    }
    return c;
}

template<typename T, int N, bool regularize_large>
inline void project_psd_impl(Eigen::Matrix<T,N,N>& S, T eps) {
    S = T(0.5) * (S + S.transpose());
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (!std::isfinite(S(i,j))) S(i,j) = (i == j) ? eps : T(0);
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
            if (!(lam(i) > T(0))) lam(i) = eps;
        }
        S = es.eigenvectors() * lam.asDiagonal() * es.eigenvectors().transpose();
    } else if constexpr (regularize_large || N <= 6) {
        Eigen::LDLT<Eigen::Matrix<T,N,N>> ldlt;
        ldlt.compute(S);
        if (ldlt.info() != Eigen::Success) {
            T min_lb = std::numeric_limits<T>::infinity();
            for (int i = 0; i < N; ++i) {
                T row_sum = T(0);
                for (int j = 0; j < N; ++j) {
                    if (j != i) row_sum += std::abs(S(i,j));
                }
                min_lb = std::min(min_lb, S(i,i) - row_sum);
            }
            if (!(min_lb > eps)) S.diagonal().array() += (eps - min_lb);
            ldlt.compute(S);
            if (ldlt.info() != Eigen::Success) S.diagonal().array() += T(10) * eps;
        }
    }
    S = T(0.5) * (S + S.transpose());
}

template<typename T, int N>
inline void project_psd_ou_ii(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,true>(S, eps);
}

template<typename T, int N>
inline void project_psd_ou_iii(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,false>(S, eps);
}

template<typename T>
inline void rot_and_B_from_wt(const Eigen::Matrix<T,3,1>& w, T t,
                              Eigen::Matrix<T,3,3>& R, Eigen::Matrix<T,3,3>& B) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const T wnorm = w.norm();
    const Matrix3 W = skew(w);

    if (wnorm < T(1e-7)) {
        const T t2 = t*t, t3 = t2*t;
        R = Matrix3::Identity() - W*t + T(0.5)*(W*W)*t2;
        B = Matrix3::Identity()*t - T(0.5)*W*t2 + (W*W)*(t3/T(6));
        return;
    }

    const T theta = wnorm * t;
    const T s = std::sin(theta), c = std::cos(theta);
    const T invw = T(1) / wnorm;
    const Matrix3 K = W * invw;
    R = Matrix3::Identity() - s*K + (T(1)-c)*(K*K);

    const T invw2 = invw * invw;
    B = Matrix3::Identity()*t
      - ((T(1)-c)*invw2)*W
      + ((t - s*invw)*invw2)*(W*W);
}

template<typename T>
inline void integral_B_ds(const Eigen::Matrix<T,3,1>& w, T step,
                          Eigen::Matrix<T,3,3>& IB) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const T wnorm = w.norm();
    const Matrix3 W = skew(w);

    if (wnorm < T(1e-7)) {
        const T T2 = step*step, T3 = T2*step, T4 = T3*step;
        IB = Matrix3::Identity()*(T(0.5)*T2)
           - W*(T(1.0/6.0)*T3)
           + (W*W)*(T(1.0/24.0)*T4);
        return;
    }

    const T theta = wnorm * step;
    const T s = std::sin(theta), c = std::cos(theta);
    const T invw = T(1) / wnorm;
    const T invw2 = invw * invw;
    IB = Matrix3::Identity()*(T(0.5)*step*step)
       - ((step - s*invw)*invw2)*W
       + ((T(0.5)*step*step + (c-T(1))*invw2)*invw2)*(W*W);
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

template<typename T, int N>
inline void regularize_psd_if_needed(Eigen::Matrix<T,N,N>& S) {
    S = T(0.5) * (S + S.transpose());
    T scale = std::max(T(1), S.cwiseAbs().maxCoeff());
    const T tol = T(64) * std::numeric_limits<T>::epsilon() * scale;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (!std::isfinite(S(i,j))) S(i,j) = (i == j) ? tol : T(0);
        }
    }

    Eigen::LDLT<Eigen::Matrix<T,N,N>> ldlt(S);
    if (ldlt.info() == Eigen::Success && ldlt.vectorD().minCoeff() >= -tol) return;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,N,N>> es(S);
    if (es.info() != Eigen::Success) {
        S.diagonal().array() += tol;
        S = T(0.5) * (S + S.transpose());
        return;
    }

    Eigen::Matrix<T,N,1> eigenvalues = es.eigenvalues();
    for (int i = 0; i < N; ++i) eigenvalues(i) = std::max(T(0), eigenvalues(i));
    S = es.eigenvectors() * eigenvalues.asDiagonal() * es.eigenvectors().transpose();
    S = T(0.5) * (S + S.transpose());
}

template<typename T>
inline bool periodic_update_due(T dt, T period, T& elapsed) {
    if (!(dt > T(0)) || !std::isfinite(dt) || !(period > T(0)) || !std::isfinite(period)) return false;
    const T total = elapsed + dt;
    const T tol = T(16) * std::numeric_limits<T>::epsilon() * std::max(T(1), period);
    if (total + tol < period) {
        elapsed = total;
        return false;
    }
    elapsed = (total >= period) ? std::fmod(total, period) : T(0);
    if (!(elapsed >= T(0)) || !std::isfinite(elapsed) || elapsed >= period) elapsed = T(0);
    return true;
}

template<typename T, int Integrations>
struct IntegratedOUChain;

template<typename T>
struct IntegratedOUChain<T,2> {
    using MatrixAxis = Eigen::Matrix<T,3,3>;

    static void transition(T tau, T h, MatrixAxis& Phi) {
        const auto P = make_prims<T>(h, tau);
        const T phi_va = -tau * P.em1;
        const auto coeffs = safe_phi_A_coeffs<T>(h, tau);
        Phi.setZero();
        Phi(0,0)=T(1); Phi(0,2)=phi_va;
        Phi(1,0)=h;    Phi(1,1)=T(1); Phi(1,2)=coeffs.phi_pa;
        Phi(2,2)=std::min(P.alpha, T(1));
    }

    static void process_covariance(T tau, T h, T sigma2, MatrixAxis& Qd) {
        const T tau_eff = std::max(tau, T(1e-7));
        const T inv = T(1) / tau_eff;
        const T x = h * inv;
        if (std::abs(x) < T(1e-2)) {
            const T i2=inv*inv, i3=i2*inv, i4=i3*inv, i5=i4*inv;
            const T i6=i5*inv, i7=i6*inv, i8=i7*inv, i9=i8*inv;
            const T h2=h*h, h3=h2*h, h4=h3*h, h5=h4*h;
            const T h6=h5*h, h7=h6*h, h8=h7*h, h9=h8*h;
            Qd.setZero();
            Qd(0,0)=sigma2*(T(2)/T(3)*h3*inv-T(1)/T(2)*h4*i2+T(7)/T(30)*h5*i3-T(1)/T(12)*h6*i4+T(31)/T(1260)*h7*i5-T(1)/T(160)*h8*i6+T(127)/T(90720)*h9*i7);
            Qd(0,1)=sigma2*(T(1)/T(4)*h4*inv-T(1)/T(6)*h5*i2+T(5)/T(72)*h6*i3-T(1)/T(45)*h7*i4+T(17)/T(2880)*h8*i5-T(41)/T(30240)*h9*i6);
            Qd(0,2)=sigma2*(h2*inv-h3*i2+T(7)/T(12)*h4*i3-T(1)/T(4)*h5*i4+T(31)/T(360)*h6*i5-T(1)/T(40)*h7*i6+T(127)/T(20160)*h8*i7-T(17)/T(12096)*h9*i8);
            Qd(1,0)=Qd(0,1);
            Qd(1,1)=sigma2*(T(1)/T(10)*h5*inv-T(1)/T(18)*h6*i2+T(5)/T(252)*h7*i3-T(1)/T(180)*h8*i4+T(17)/T(12960)*h9*i5);
            Qd(1,2)=sigma2*(T(1)/T(3)*h3*inv-T(1)/T(3)*h4*i2+T(11)/T(60)*h5*i3-T(13)/T(180)*h6*i4+T(19)/T(840)*h7*i5-T(1)/T(168)*h8*i6+T(247)/T(181440)*h9*i7);
            Qd(2,0)=Qd(0,2); Qd(2,1)=Qd(1,2);
            Qd(2,2)=sigma2*(T(2)*h*inv-T(2)*h2*i2+T(4)/T(3)*h3*i3-T(2)/T(3)*h4*i4+T(4)/T(15)*h5*i5-T(4)/T(45)*h6*i6+T(8)/T(315)*h7*i7-T(2)/T(315)*h8*i8+T(4)/T(2835)*h9*i9);
            regularize_psd_if_needed<T,3>(Qd);
            return;
        }

        const T a=std::exp(-x), a2=a*a, qc=T(2)*sigma2*inv;
        const T t2=tau_eff*tau_eff, t3=t2*tau_eff, t4=t3*tau_eff, t5=t4*tau_eff;
        const T x2=x*x, x3=x2*x;
        const T K00=t3*(-a2+T(4)*a+T(2)*x-T(3))/T(2);
        const T K01=t4*(a2+T(2)*a*(x-T(1))+x2-T(2)*x+T(1))/T(2);
        const T K02=t2*(a2-T(2)*a+T(1))/T(2);
        const T K11=t5*(-a2/T(2)-T(2)*a*x+x3/T(3)-x2+x+T(1)/T(2));
        const T K12=t3*(-a2-T(2)*a*x+T(1))/T(2);
        const T K22=tau_eff*(T(1)-a2)/T(2);
        Qd << qc*K00,qc*K01,qc*K02,
              qc*K01,qc*K11,qc*K12,
              qc*K02,qc*K12,qc*K22;
        regularize_psd_if_needed<T,3>(Qd);
    }
};

template<typename T>
struct IntegratedOUChain<T,3> {
    using MatrixAxis = Eigen::Matrix<T,4,4>;

    static void transition(T tau, T h, MatrixAxis& Phi) {
        const auto P = make_prims<T>(h, tau);
        const T phi_va = -tau*P.em1;
        const auto coeffs = safe_phi_A_coeffs<T>(h, tau);
        Phi.setZero();
        Phi(0,0)=T(1); Phi(0,3)=phi_va;
        Phi(1,0)=h; Phi(1,1)=T(1); Phi(1,3)=coeffs.phi_pa;
        Phi(2,0)=T(0.5)*h*h; Phi(2,1)=h; Phi(2,2)=T(1); Phi(2,3)=coeffs.phi_Sa;
        Phi(3,3)=std::min(P.alpha,T(1));
    }

    static void process_covariance(T tau, T h, T sigma2, MatrixAxis& Qd) {
        const T tau_eff = std::max(tau, T(1e-7));
        const T inv = T(1) / tau_eff;
        const T x = h * inv;

        Eigen::Matrix<T,3,3> marginal;
        IntegratedOUChain<T,2>::process_covariance(tau_eff, h, sigma2, marginal);
        Qd.setZero();
        const int idx[3] = {0,1,3};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) Qd(idx[i],idx[j]) = marginal(i,j);
        }

        T qvS, qpS, qSS, qSa;
        if (std::abs(x) < T(1e-2)) {
            const T i2=inv*inv, i3=i2*inv, i4=i3*inv, i5=i4*inv, i6=i5*inv;
            const T h2=h*h, h3=h2*h, h4=h3*h, h5=h4*h;
            const T h6=h5*h, h7=h6*h, h8=h7*h, h9=h8*h;
            qvS=sigma2*(T(1)/T(15)*h5*inv-T(1)/T(24)*h6*i2+T(41)/T(2520)*h7*i3-T(7)/T(1440)*h8*i4+T(109)/T(90720)*h9*i5);
            qpS=sigma2*(T(1)/T(36)*h6*inv-T(1)/T(72)*h7*i2+T(13)/T(2880)*h8*i3-T(1)/T(864)*h9*i4);
            qSS=sigma2*(T(1)/T(126)*h7*inv-T(1)/T(288)*h8*i2+T(13)/T(12960)*h9*i3);
            qSa=sigma2*(T(1)/T(12)*h4*inv-T(1)/T(12)*h5*i2+T(2)/T(45)*h6*i3-T(1)/T(60)*h7*i4+T(11)/T(2240)*h8*i5-T(73)/T(60480)*h9*i6);
        } else {
            const T a=std::exp(-x), a2=a*a, qc=T(2)*sigma2*inv;
            const T t4=tau_eff*tau_eff*tau_eff*tau_eff;
            const T t5=t4*tau_eff, t6=t5*tau_eff, t7=t6*tau_eff;
            const T x2=x*x, x3=x2*x, x4=x3*x, x5=x4*x;
            const T K02=t5*(-T(3)*a2+T(3)*a*(x2+T(4))+x3-T(3)*x2+T(6)*x-T(9))/T(6);
            const T K12=t6*(a2/T(2)+a*(-x2+T(2)*x-T(2))/T(2)+x4/T(8)-x3/T(2)+x2-x+T(1)/T(2));
            const T K22=t7*(-a2/T(2)+a*x2+T(2)*a+x5/T(20)-x4/T(4)+T(2)*x3/T(3)-x2+x-T(3)/T(2));
            const T K23=t4*(a2-a*(x2+T(2))+T(1))/T(2);
            qvS=qc*K02;
            qpS=qc*K12;
            qSS=qc*K22;
            qSa=qc*K23;
        }

        Qd(0,2)=qvS; Qd(2,0)=qvS;
        Qd(1,2)=qpS; Qd(2,1)=qpS;
        Qd(2,2)=qSS;
        Qd(2,3)=qSa; Qd(3,2)=qSa;
        regularize_psd_if_needed<T,4>(Qd);
    }
};

} // namespace ocean_imu::kalman::ou_detail
