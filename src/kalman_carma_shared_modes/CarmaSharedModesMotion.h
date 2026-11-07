#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    SharedModesMotion<T> -- the motion block of the EXPERIMENTAL shared-modes
    QMEKF (KalmanCarmaSharedModes.h; docs/kalman-carma-shared-modes.md,
    doc/kalman_carma_shared_modes/kalman_carma_shared_modes-minimal.tex):
    the shared directional stochastic modes in their minimal realization,
    six real states per axis, 27 error states for the QMEKF.

    The model.  Per axis a complexified CARMA(4,1) mode zeta = B/(P+ D) xi is
    driven by circular complex noise xi = xi_r + i xi_i, the three axes by a
    Hermitian excitation Sigma, and the filter sees Re zeta (every
    accelerometer row reads Re a).  With P+ = N_r + i N_i (N_r, N_i real
    polynomials, P = P+ P- = N_r^2 + N_i^2) and |D|^2 = |N_r|^2 + |N_i|^2,

        Re zeta = B/P * eta,   eta = V [xi_r; xi_i],   V = [N_r, N_i] / D,

    a 1 x 2 transfer function of McMillan degree deg(P D) = 6.  The eight-state
    realified complexified mode carries two more directions per axis, which
    are unobservable for every excitation (each root of D is an eigenvalue of
    both the complex mode and its conjugate); this block is its minimal
    realization and gives the same law of p, v, a on all three axes, every
    co- and quadrature spectrum included.

    The cascade: [p, v, a, l] in the CARMA3 baseline's own coordinates
    (CarmaVerticalModel, the C41J driven at a by eta) and the two states w of
    the lossless V in input-normal coordinates (controllability Gramian I:
    A_w + A_w^T = -B_w B_w^T, C_w = -B_w e_1^T, B_w upper triangular with a
    positive diagonal), so eta = xi_r + C_w w.  The Hermitian excitation
    enters as
        Q_ab = Re(Sigma_ab) Q_R + Im(Sigma_ab) Q_I,
    Q_R from W_R = g2 G G^T, Q_I from the antisymmetric W_I = g2 G J2 G^T,
    J2 = [[0, -1], [1, 0]] (a paired discretization; minimal_detail).

    A model commit (kappa, amplitudes, coupling) holds [p, v, a, l], as the
    CARMA3 baseline does, and the normalized colouring w / s_a,
    s_a = sqrt(g2 sigma_a^2): x_w <- (s_new / s_old) x_w per axis through
    take_axis_recoordinatization().  V is lossless with a constant Gramian,
    so eta stays exactly white under any change of kappa and sigma: with the
    coupling off the block IS the CARMA3 baseline A41 in law, adaptation
    included.  Every quantity scales exactly with kappa (V_kappa(s) =
    V_1(s / kappa)), so the kappa = 1 construction (closed-form 2 x 2 algebra
    in double) is the only factorization; commits rescale.  Interface: that
    of the coupled motion blocks of Kalman3D_Wave_CARMA3<..., Motion>.

    A second, short-horizon band cross-spectrum estimate (the MONITOR) is fed
    with the same world-frame acceleration the coupling estimator sees; the
    coupling confidence (CarmaCouplingConfidence.h) compares it with the
    prior.  The monitor never enters the prior.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>

#include "kalman_carma/CarmaCoupledMotion.h"
#include "kalman_carma/CarmaCrossSpectrumEstimator.h"
#include "kalman_carma/CarmaVerticalModel.h"

namespace ocean_imu::carma_shared_modes {

namespace minimal_detail {

namespace kd = kalman::carma::detail;

// Z = X Y X^T for antisymmetric Y; Z is written antisymmetric.
template<typename T, int M>
inline void congruence_anti(const kd::Mat<T,M>& X, const kd::Mat<T,M>& Y, kd::Mat<T,M>& Z) {
    kd::Mat<T,M> XY;
    kd::mul<T,M>(X, Y, XY);
    for (int i = 0; i < M; ++i) {
        Z(i,i) = T(0);
        for (int j = i + 1; j < M; ++j) {
            T s = T(0);
            for (int k = 0; k < M; ++k) s += XY(i,k) * X(j,k);
            Z(i,j) = s;
            Z(j,i) = -s;
        }
    }
}

// Exact discretization of dx = A x dt + dB with the symmetric intensity W_R
// and the antisymmetric one W_I (the real and imaginary parts of a Hermitian
// intensity): Phi, Q_R = int e^{As} W_R e^{A's} ds, Q_I likewise.  The
// algorithm, scaling and term counts are discretize_lti()'s
// (kalman_carma/CarmaVerticalModel.h); the noise distances come from W_R,
// whose support contains W_I's.
template<typename T, int M>
inline bool discretize_lti_pair(const kd::Mat<T,M>& A, const kd::Mat<T,M>& WR, const kd::Mat<T,M>& WI, T h,
                                kd::Mat<T,M>& Phi, kd::Mat<T,M>& QR, kd::Mat<T,M>& QI) {
    using Mat = kd::Mat<T,M>;
    Phi.setIdentity();
    QR.setZero();
    QI.setZero();
    if (!(h > T(0)) || !std::isfinite(h)) return h == T(0);
    const T norm_a = kd::inf_norm<T,M>(A);
    if (!std::isfinite(norm_a)) return false;
    int s = 0;
    T hs = h;
    while (norm_a * hs > T(0.25) && s < 60) { hs *= T(0.5); ++s; }
    int d[M];
    kd::noise_distance<T,M>(A, WR, d);
    T hp[2 * M + 1];
    hp[0] = T(1);
    for (int i = 1; i <= 2 * M; ++i) hp[i] = hp[i - 1] * hs;
    Mat As;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            const int e = d[j] - d[i];
            As(i,j) = (e >= 0) ? A(i,j) * hp[e] : A(i,j) / hp[-e];
        }
    constexpr int K_TAIL = (std::numeric_limits<T>::digits > 30) ? 13 : 7;
    constexpr int K_Q = 2 * (M - 1) + K_TAIL;
    constexpr int K_F = (M - 1) + K_TAIL;
    Mat Y = Mat::Identity(), XR = WR * hs, XI = WI * hs;
    Mat Fs = Y, QRs = XR, QIs = XI, AY, AX;
    for (int k = 1; k <= K_Q; ++k) {
        const T c = hs / T(k + 1);
        if (k <= K_F) {
            kd::mul<T,M>(As, Y, AY);
            Y = AY * (hs / T(k));
            Fs += Y;
        }
        kd::mul<T,M>(As, XR, AX);
        for (int i = 0; i < M; ++i)
            for (int j = i; j < M; ++j) {
                const T v = (AX(i,j) + AX(j,i)) * c;
                XR(i,j) = v; XR(j,i) = v;
            }
        QRs += XR;
        kd::mul<T,M>(As, XI, AX);
        for (int i = 0; i < M; ++i) {
            XI(i,i) = T(0);
            for (int j = i + 1; j < M; ++j) {
                const T v = (AX(i,j) - AX(j,i)) * c;
                XI(i,j) = v; XI(j,i) = -v;
            }
        }
        QIs += XI;
    }
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            const int e = d[i] - d[j];
            Phi(i,j) = (e >= 0) ? Fs(i,j) * hp[e] : Fs(i,j) / hp[-e];
            QR(i,j) = QRs(i,j) * hp[d[i] + d[j]];
            QI(i,j) = QIs(i,j) * hp[d[i] + d[j]];
        }
    Mat FQF, FF;
    for (int i = 0; i < s; ++i) {
        kd::congruence<T,M>(Phi, QR, FQF);
        QR += FQF;
        congruence_anti<T,M>(Phi, QI, FQF);
        QI += FQF;
        kd::mul<T,M>(Phi, Phi, FF);
        Phi = FF;
    }
    for (int i = 0; i < M; ++i) {
        QI(i,i) = T(0);
        for (int j = i + 1; j < M; ++j) {
            const T r = T(0.5) * (QR(i,j) + QR(j,i)), a = T(0.5) * (QI(i,j) - QI(j,i));
            QR(i,j) = r; QR(j,i) = r;
            QI(i,j) = a; QI(j,i) = -a;
        }
    }
    return Phi.allFinite() && QR.allFinite() && QI.allFinite();
}

// Stationary P_R, P_I (A P + P A^T + W = 0 for both parts) by doubling, as
// stationary_covariance_doubling().
template<typename T, int M>
inline bool stationary_pair(const kd::Mat<T,M>& A, const kd::Mat<T,M>& WR, const kd::Mat<T,M>& WI,
                            kd::Mat<T,M>& PR, kd::Mat<T,M>& PI) {
    using Mat = kd::Mat<T,M>;
    const T norm_a = kd::inf_norm<T,M>(A);
    PR.setZero();
    PI.setZero();
    if (!(norm_a > T(0)) || !std::isfinite(norm_a)) return false;
    Mat Phi;
    if (!discretize_lti_pair<T,M>(A, WR, WI, T(0.25) / norm_a, Phi, PR, PI)) return false;
    const T tol = T(4) * std::numeric_limits<T>::epsilon();
    Mat FPF, FF;
    for (int it = 0; it < 96; ++it) {
        kd::congruence<T,M>(Phi, PR, FPF);
        PR += FPF;
        congruence_anti<T,M>(Phi, PI, FPF);
        PI += FPF;
        kd::mul<T,M>(Phi, Phi, FF);
        Phi = FF;
        if (kd::max_abs<T,M>(Phi) <= tol) return PR.allFinite() && PI.allFinite();
    }
    return false;
}

// StepDiscretization (kalman_carma/CarmaVerticalModel.h) for the pair: the
// same reference step, Taylor neighbourhood and series, with dQ_I/dh =
// F W_I F^T and Q_I^(k+1) = A Q_I^(k) + Q_I^(k) A^T.
template<typename T, int M>
class StepDiscretizationPair {
public:
    using Mat = kd::Mat<T,M>;
    static constexpr int KF = (M - 1 > 3) ? M - 1 : 3;
    static constexpr int KQ = 2 * KF;

    bool advance(uint32_t generation, T h) {
        if (generation != generation_) return false;
        if (h == h_) return true;
        if (!ref_ok_) return false;
        const T d = h - h0_;
        const T ad = std::abs(d);
        if (!(ad <= kRange * h0_) || !(norm_s_ * ad <= kRange)) return false;
        if (!series_ready_) build_series_();
        evaluate_(d);
        h_ = h;
        return true;
    }
    void refresh(const Mat& A, const Mat& WR, const Mat& WI, uint32_t generation, T h) {
        ref_ok_ = discretize_lti_pair<T,M>(A, WR, WI, h, F_, QR_, QI_) && h > T(0);
        A_ = A; WR_ = WR; WI_ = WI;
        Fc_[0] = F_; QRc_[0] = QR_; QIc_[0] = QI_;
        norm_s_ = scaled_inf_norm_(A, WR, h);
        h0_ = h; h_ = h;
        generation_ = generation;
        series_ready_ = false;
    }
    const Mat& Phi() const noexcept { return F_; }
    const Mat& QR() const noexcept { return QR_; }
    const Mat& QI() const noexcept { return QI_; }

private:
    static constexpr T kRange = T(0.0625);
    static T scaled_inf_norm_(const Mat& A, const Mat& W, T h) {
        int d[M];
        kd::noise_distance<T,M>(A, W, d);
        T hp[2 * M + 1];
        hp[0] = T(1);
        for (int i = 1; i <= 2 * M; ++i) hp[i] = hp[i - 1] * h;
        T m = T(0);
        for (int i = 0; i < M; ++i) {
            T r = T(0);
            for (int j = 0; j < M; ++j) {
                const int e = d[j] - d[i];
                r += std::abs(A(i,j)) * ((e >= 0) ? hp[e] : T(1) / hp[-e]);
            }
            m = std::max(m, r);
        }
        return m;
    }
    void build_series_() {
        Mat X;
        for (int k = 1; k <= KF; ++k) {
            kd::mul<T,M>(A_, Fc_[k - 1], X);
            Fc_[k] = X * (T(1) / T(k));
        }
        kd::congruence<T,M>(Fc_[0], WR_, QRc_[1]);
        congruence_anti<T,M>(Fc_[0], WI_, QIc_[1]);
        for (int k = 2; k <= KQ; ++k) {
            const T c = T(1) / T(k);
            kd::mul<T,M>(A_, QRc_[k - 1], X);
            for (int i = 0; i < M; ++i)
                for (int j = i; j < M; ++j) {
                    const T v = (X(i,j) + X(j,i)) * c;
                    QRc_[k](i,j) = v; QRc_[k](j,i) = v;
                }
            kd::mul<T,M>(A_, QIc_[k - 1], X);
            for (int i = 0; i < M; ++i) {
                QIc_[k](i,i) = T(0);
                for (int j = i + 1; j < M; ++j) {
                    const T v = (X(i,j) - X(j,i)) * c;
                    QIc_[k](i,j) = v; QIc_[k](j,i) = -v;
                }
            }
        }
        series_ready_ = true;
    }
    void evaluate_(T d) {
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < M; ++j) {
                T f = Fc_[KF](i,j);
                for (int k = KF - 1; k >= 0; --k) f = Fc_[k](i,j) + d * f;
                F_(i,j) = f;
            }
        for (int i = 0; i < M; ++i)
            for (int j = i; j < M; ++j) {
                T qr = QRc_[KQ](i,j), qi = QIc_[KQ](i,j);
                for (int k = KQ - 1; k >= 1; --k) { qr = QRc_[k](i,j) + d * qr; qi = QIc_[k](i,j) + d * qi; }
                qr = QRc_[0](i,j) + d * qr;
                qi = QIc_[0](i,j) + d * qi;
                QR_(i,j) = qr; QR_(j,i) = qr;
                QI_(i,j) = (i == j) ? T(0) : qi;
                QI_(j,i) = (i == j) ? T(0) : -qi;
            }
    }

    Mat F_ = Mat::Identity(), QR_ = Mat::Zero(), QI_ = Mat::Zero();
    Mat A_ = Mat::Zero(), WR_ = Mat::Zero(), WI_ = Mat::Zero();
    Mat Fc_[KF + 1];
    Mat QRc_[KQ + 1];
    Mat QIc_[KQ + 1];
    T norm_s_ = T(0), h0_ = T(-1), h_ = T(-1);
    uint32_t generation_ = 0xffffffffu;
    bool ref_ok_ = false, series_ready_ = false;
};

// The polynomials of the complexified C41J at scale kappa: P+(s) = s^2 +
// b1 s + b0 (the upper-half-plane poles of P), D(s) = s^2 + d1 s + d0 with
// |D(iw)|^2 = (|P+(iw)|^2 + |P+(-iw)|^2) / 2.
inline bool shared_mode_polys(const kalman::carma::CarmaShape& sh, double kappa, std::complex<double>& b1,
                              std::complex<double>& b0, double& d1, double& d0) {
    using C = std::complex<double>;
    constexpr double k2Pi = 6.283185307179586476925286766559;
    const double w1 = k2Pi * kappa * sh.f1_hz, z1 = sh.zeta1;
    const double w2 = k2Pi * kappa * sh.f2_hz, z2 = sh.zeta2;
    if (!(z1 > 0.0 && z1 < 1.0 && z2 > 0.0 && z2 < 1.0 && w1 > 0.0 && w2 > 0.0)) return false;
    const C p1(-z1 * w1, w1 * std::sqrt(1.0 - z1 * z1));
    const C p2(-z2 * w2, w2 * std::sqrt(1.0 - z2 * z2));
    b1 = -(p1 + p2);
    b0 = p1 * p2;
    auto pplus = [&](double w) { const C s(0.0, w); return s * s + b1 * s + b0; };
    auto q = [&](double w) { return std::norm(pplus(w)) + std::norm(pplus(-w)); };
    const double gamma = 0.5 * q(0.0);
    const double wr = std::max(w1, w2), u1 = wr * wr;
    const double beta = (0.5 * q(wr) - u1 * u1 - gamma) / u1;
    d0 = std::sqrt(gamma);
    const double d1sq = beta + 2.0 * d0;
    if (!(d1sq > 0.0) || !std::isfinite(d1sq)) return false;
    d1 = std::sqrt(d1sq);
    return std::isfinite(d0) && d0 > 0.0;
}

// Input-normal realization (A_w, B_w, C_w) of the co-inner V = [N_r, N_i] / D
// with feedthrough [1, 0]: controllability Gramian I, B_w upper triangular
// with a positive diagonal (fixes the orthogonal freedom; the same at every
// kappa), C_w = -B_w e_1^T.  Closed-form 2 x 2 algebra throughout (no
// generic decompositions: they would add tens of kB of code on the target).
// Returns false if V is degenerate.
inline bool coinner_input_normal(const std::complex<double>& b1, const std::complex<double>& b0, double d1, double d0,
                                 Eigen::Matrix2d& Aw, Eigen::Matrix2d& Bw, Eigen::RowVector2d& Cw) {
    if (!(d1 > 0.0) || !(d0 > 0.0)) return false;
    // Observer canonical form A0 = [[-d1, 1], [-d0, 0]], C0 = [1, 0]:
    // V_r = 1 + ((Re b1 - d1) s + Re b0 - d0) / D, V_i = (Im b1 s + Im b0) / D.
    const double a = -d1, c = -d0;
    Eigen::Matrix2d B0;
    B0 << b1.real() - d1, b1.imag(), b0.real() - d0, b0.imag();
    const double q11 = B0.row(0).squaredNorm(), q12 = B0.row(0).dot(B0.row(1)), q22 = B0.row(1).squaredNorm();
    // A0 X + X A0^T + B0 B0^T = 0, X symmetric:
    //   2 a x11 + 2 x12 + q11 = 0,  c x11 + a x12 + x22 + q12 = 0,  2 c x12 + q22 = 0.
    const double x12 = -q22 / (2.0 * c);
    const double x11 = -(q11 + 2.0 * x12) / (2.0 * a);
    const double x22 = -(q12 + c * x11 + a * x12);
    // X = L L^T (Cholesky), w = L^-1 w0.
    if (!(x11 > 0.0)) return false;
    const double l11 = std::sqrt(x11), l21 = x12 / l11, r22 = x22 - l21 * l21;
    if (!(r22 > 0.0)) return false;
    const double l22 = std::sqrt(r22);
    Eigen::Matrix2d L, Li, A0;
    L << l11, 0.0, l21, l22;
    Li << 1.0 / l11, 0.0, -l21 / (l11 * l22), 1.0 / l22;
    A0 << a, 1.0, c, 0.0;
    Aw = Li * A0 * L;
    Bw = Li * B0;
    Cw << l11, 0.0;                                   // C0 L
    // Givens rotation O with O B_w upper triangular, positive diagonal.
    const double r = std::hypot(Bw(0,0), Bw(1,0));
    if (!(r > 0.0)) return false;
    Eigen::Matrix2d O;
    O << Bw(0,0) / r, Bw(1,0) / r, -Bw(1,0) / r, Bw(0,0) / r;
    if ((O.row(1) * Bw.col(1))(0) < 0.0) O.row(1) *= -1.0;
    Aw = O * Aw * O.transpose();
    Bw = O * Bw;
    Cw = Cw * O.transpose();
    // Exact structure (round-off residue would add spurious couplings to the
    // discretization's structural scaling): B_w upper triangular,
    // C_w = -B_w e_1^T after the lossless check.
    Bw(1,0) = 0.0;
    const double sc = Bw.norm();
    if (!((Cw + Bw.col(0).transpose()).norm() <= 1e-9 * sc) ||
        !((Aw + Aw.transpose() + Bw * Bw.transpose()).norm() <= 1e-9 * sc * sc)) return false;
    Cw = -Bw.col(0).transpose();
    return Aw.allFinite() && Bw.allFinite();
}

}  // namespace minimal_detail

template<typename T>
class SharedModesMotion {
public:
    static constexpr bool kUnitReadout = true;   // p, v, a are states
    static constexpr int NB = 6;
    static constexpr int NL = 3 * NB;
    static constexpr int MAX_ACC_TERMS = 1;
    using Vector3 = Eigen::Matrix<T,3,1>;
    using MatB = Eigen::Matrix<T,NB,NB>;
    using Model = kalman::carma::CarmaVerticalModel<T,4>;
    using Hermitian = kalman::carma::Hermitian3<T>;
    using Coupling = kalman::carma::CarmaCoupling<T>;
    using Estimator = kalman::carma::CarmaCrossSpectrumEstimator<T>;
    using Mat6d = Eigen::Matrix<double,6,6>;

    static constexpr int pos_local(int ax) { return NB * ax; }
    static constexpr int vel_local(int ax) { return NB * ax + 1; }
    static constexpr int acc_local(int ax) { return NB * ax + 2; }
    int accel_terms(int ax, int (&idx)[4], T (&w)[4]) const {
        idx[0] = acc_local(ax);
        w[0] = T(1);
        return 1;
    }

    explicit SharedModesMotion(const kalman::carma::CarmaShape& shape = kalman::carma::CARMA41_JONSWAP_SHAPE)
        : vmodel_(shape), hmodel_(shape) {
        vmodel_.configure(T(1), T(1));
        hmodel_.configure(T(1), T(1));
        sig_ = Vector3::Ones();
        valid_ = build_unit_();
        rebuild_();
    }

    // Wrapper interface (Kalman3D_Wave_CARMA3 forwards to these).
    bool set_vertical(T kappa, T accel_std) {
        if (!vmodel_.configure(kappa, accel_std)) return false;
        sig_(2) = accel_std;
        dirty_ = true;
        return true;
    }
    bool set_horizontal(T kappa, T sx, T sy, T /*r: the coupling source supplies it*/) {
        if (!(sx > T(0)) || !(sy > T(0)) || !std::isfinite(sx) || !std::isfinite(sy)) return false;
        if (!hmodel_.configure(kappa, T(1))) return false;
        sig_(0) = sx;
        sig_(1) = sy;
        dirty_ = true;
        return true;
    }
    bool set_shape(const kalman::carma::CarmaShape& sh) {
        if (!vmodel_.set_shape(sh) || !hmodel_.set_shape(sh)) return false;
        valid_ = build_unit_();
        have_prev_ = false;   // a new shape is a new realization: nothing to carry over
        dirty_ = true;
        return valid_;
    }
    const Model& vertical_model() const noexcept { return vmodel_; }
    const Model& horizontal_model() const noexcept { return hmodel_; }
    Eigen::Matrix<T,2,2> horizontal_sigma() const { return Sigma_.re.template topLeftCorner<2,2>(); }
    Vector3 accel_std() const { return sig_; }
    const Hermitian& excitation() const noexcept { return Sigma_; }
    Coupling& coupling() noexcept { return coupling_; }
    const Coupling& coupling() const noexcept { return coupling_; }
    T operating_frequency_hz() const {
        return static_cast<T>(static_cast<double>(vmodel_.kappa()) * vmodel_.base_zero_crossing_hz());
    }
    // The QMEKF's accelerometer update: the coupling estimator and the monitor.
    void observe(const Vector3& a_world) {
        monitor_.update(a_world, dt_, operating_frequency_hz());
        coupling_.observe(a_world, dt_, operating_frequency_hz());
    }
    Estimator& monitor() noexcept { return monitor_; }
    const Estimator& monitor() const noexcept { return monitor_; }
    void set_heading(T yaw) { coupling_.set_heading(yaw); }
    T step() const noexcept { return dt_; }
    bool valid() const noexcept { return valid_; }

    void refresh(T Ts) {
        if (Ts > T(0) && std::isfinite(Ts)) dt_ = Ts;
        if (dirty_) rebuild_();
        if (!disc_.advance(gen_, Ts)) disc_.refresh(A_, WR_, WI_, gen_, Ts);
    }

    void apply_F(const T (&in)[NL], T (&out)[NL]) const {
        const MatB& Phi = disc_.Phi();
        for (int ax = 0; ax < 3; ++ax)
            for (int i = 0; i < NB; ++i) {
                T acc = T(0);
                for (int j = 0; j < NB; ++j) acc += Phi(i,j) * in[ax * NB + j];
                out[ax * NB + i] = acc;
            }
    }

    template<class Mat>
    void propagate(Mat& P, int off) const {
        const MatB& Phi = disc_.Phi();
        for (int a = 0; a < 3; ++a)
            for (int b = a; b < 3; ++b)
                kalman::carma::coupled_detail::congruence_block<T,NB,NB>(P, off + a * NB, off + b * NB, Phi, Phi, a == b);
        add_structured_(P, off, disc_.QR(), disc_.QI());
    }

    template<class Mat>
    bool stationary(Mat& P, int off) const {
        MatB PR, PI;
        if (!minimal_detail::stationary_pair<T,NB>(A_, WR_, WI_, PR, PI)) return false;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                for (int i = 0; i < NB; ++i)
                    for (int j = 0; j < NB; ++j) P(off + a * NB + i, off + b * NB + j) = T(0);
        add_structured_(P, off, PR, PI);
        return true;
    }

    // Re-coordinatization after a commit (file header): x_axis <- M[axis] x_axis
    // on every axis block; false when none is pending.
    bool take_axis_recoordinatization(MatB (&M)[3]) {
        if (!pending_) return false;
        for (int a = 0; a < 3; ++a) M[a] = M_[a];
        pending_ = false;
        return true;
    }

    // Diagnostics (tests, analysis): continuous dynamics, noise input
    // (columns: Re and Im of the circular noise), unit-acceleration noise
    // gain, the step discretization.
    const MatB& block_A() const noexcept { return A_; }
    const Eigen::Matrix<T,NB,2>& block_G() const noexcept { return G_; }
    T noise_gain2() const noexcept { return g2_; }
    const MatB& block_F() const noexcept { return disc_.Phi(); }
    const MatB& block_QR() const noexcept { return disc_.QR(); }
    const MatB& block_QI() const noexcept { return disc_.QI(); }

private:
    // P_ab += Re(S_ab) MR + Im(S_ab) MI over the upper block triangle, mirrored
    // (MR symmetric, MI antisymmetric: block (b, a) receives the transpose).
    template<class Mat>
    void add_structured_(Mat& P, int off, const MatB& MR, const MatB& MI) const {
        for (int a = 0; a < 3; ++a)
            for (int b = a; b < 3; ++b) {
                const T sr = Sigma_.re(a,b);
                const T si = a == b ? T(0) : Sigma_.im(a,b);
                if (sr == T(0) && si == T(0)) continue;
                for (int i = 0; i < NB; ++i)
                    for (int j = (a == b ? i : 0); j < NB; ++j) {
                        const T v = sr * MR(i,j) + si * MI(i,j);
                        const int r = off + a * NB + i, c = off + b * NB + j;
                        P(r, c) += v;
                        if (r != c) P(c, r) += v;
                    }
            }
    }

    // kappa = 1 colouring (double): the only factorization.
    bool build_unit_() {
        std::complex<double> b1, b0;
        double d1 = 0.0, d0 = 0.0;
        if (!minimal_detail::shared_mode_polys(vmodel_.shape(), 1.0, b1, b0, d1, d0)) return false;
        return minimal_detail::coinner_input_normal(b1, b0, d1, d0, Aw1_, Bw1_, Cw1_);
    }

    void rebuild_() {
        dirty_ = false;
        if (!valid_) return;
        const double kappa = static_cast<double>(vmodel_.kappa());
        const double g2 = kappa / vmodel_.base_var_a();      // unit-acceleration-variance gain at kappa
        Eigen::Matrix<double,4,4> A4;
        kalman::carma::carma_dynamics<4>(vmodel_.shape(), kappa, A4);
        const double rk = std::sqrt(kappa);
        Mat6d A = Mat6d::Zero();
        A.template topLeftCorner<4,4>() = A4;
        A.template block<1,2>(2, 4) = rk * Cw1_;
        A.template bottomRightCorner<2,2>() = kappa * Aw1_;
        Eigen::Matrix<double,6,2> G = Eigen::Matrix<double,6,2>::Zero();
        G(2,0) = 1.0;
        G.template bottomRows<2>() = rk * Bw1_;
        Eigen::Matrix2d J2;
        J2 << 0.0, -1.0, 1.0, 0.0;
        A_ = A.template cast<T>();
        G_ = G.template cast<T>();
        WR_ = (g2 * G * G.transpose()).template cast<T>();
        WI_ = (g2 * G * J2 * G.transpose()).template cast<T>();
        g2_ = static_cast<T>(g2);
        const std::complex<T> Rr[3] = {T(1), T(1), T(1)};
        Sigma_ = coupling_.excitation(sig_, Rr);

        // Commit re-coordinatization: hold the normalized colouring (file header).
        bool any = false;
        for (int a = 0; a < 3; ++a) {
            M_[a].setIdentity();
            const double s = std::sqrt(g2) * static_cast<double>(sig_(a));
            if (have_prev_ && s != scale_[a] && scale_[a] > 0.0 && s > 0.0) {
                const T r = static_cast<T>(s / scale_[a]);
                M_[a](4,4) = r;
                M_[a](5,5) = r;
                any = true;
            }
            scale_[a] = s;
        }
        have_prev_ = true;
        pending_ = any;
        ++gen_;
    }

    Model vmodel_, hmodel_;
    Vector3 sig_ = Vector3::Ones();
    Hermitian Sigma_{};
    MatB A_ = MatB::Zero(), WR_ = MatB::Zero(), WI_ = MatB::Zero();
    Eigen::Matrix<T,NB,2> G_ = Eigen::Matrix<T,NB,2>::Zero();
    T g2_ = T(1);
    minimal_detail::StepDiscretizationPair<T,NB> disc_;
    // kappa = 1 colouring (double).
    Eigen::Matrix2d Aw1_ = Eigen::Matrix2d::Zero(), Bw1_ = Eigen::Matrix2d::Zero();
    Eigen::RowVector2d Cw1_ = Eigen::RowVector2d::Zero();
    // Commit bookkeeping.
    MatB M_[3] = {MatB::Identity(), MatB::Identity(), MatB::Identity()};
    double scale_[3] = {0.0, 0.0, 0.0};
    bool have_prev_ = false, pending_ = false;
    uint32_t gen_ = 0;
    T dt_ = T(0.005);
    bool dirty_ = true, valid_ = false;
    Coupling coupling_;
    Estimator monitor_;
};

}  // namespace ocean_imu::carma_shared_modes
