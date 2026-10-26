#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
  Released under the MIT License

  CarmaVerticalModel -- EXPERIMENTAL CARMA(p,q) marine-motion prior of one
  translational axis (the vertical axis and, with a horizontal shape, each
  horizontal axis of Kalman3D_Wave_CARMA3).

  Research alternative to the OU-III chain [v, p, S, a_w].  It is not a
  shipping estimator component; see docs/carma-three-axis.md for the
  derivation, the benchmark protocol and the measured results.

  Model
  -----
  A Gaussian CARMA displacement h with transfer function H(s) = g B(s)/P(s),
  monic B and P, deg P - deg B = 3, and two-sided PSD

      S_h(w) = g^2 |B(iw)|^2 / |P(iw)|^2,   Var h = (1/2pi) int S_h dw.

      CARMA(3,0)  B = 1,                        P = (s+lambda) L1(s)
      CARMA(4,1)  B = s + w_z0,                 P = L1(s) L2(s)
      CARMA(5,2)  B = s^2 + 2 zeta_z w_z s + w_z^2,
                                                P = (s+lambda) L1(s) L2(s)

  with L_j(s) = s^2 + 2 zeta_j w_j s + w_j^2.

  State coordinates
  -----------------
  The realization is written in the coordinates x = [h, v, a, l_1 .. l_{N-3}],
  obtained from the factored (cascade) realization by the observability
  change of basis rows C, CA, CA^2, so that

      h = x_0 = C x,   v = x_1 = C A x,   a = x_2 = C A^2 x,
      dx = A x dt + g e_a dW.

  Relative degree three gives C G = C A G = 0 and C A^2 G = g != 0: h, v and a
  are mean-square continuous, a is not differentiable (its derivative carries
  the white driving noise).  The remaining coordinates are latent states of
  the factored realization (the input r of L2 and its derivative); they are
  spectral-shaping memory, not additional physical variables.  In these
  coordinates OU-III's vertical chain is the special case a' = -a/tau + noise
  with no h or v feedback; CARMA feeds displacement and velocity back into the
  jerk equation, which makes the displacement stationary.

  Scaling
  -------
  A fixed spectral shape (frequencies at kappa = 1) is moved by one frequency
  scale kappa: every pole and zero is multiplied by kappa.  With unit noise
  gain the acceleration variance then scales exactly as Var_kappa(a) =
  Var_1(a)/kappa, so the driving gain for a requested acceleration standard
  deviation needs no Lyapunov solve at run time.

  Discretization
  --------------
  Phi = exp(A h) and Q = int_0^h exp(A t) W exp(A^T t) dt, W = g^2 e_a e_a^T,
  are computed exactly for any step h: Taylor series of both integrands at
  h/2^s with ||A h/2^s||_inf <= 1/4, then s exact doublings
      Q(2h) = Q(h) + Phi(h) Q(h) Phi(h)^T,   Phi(2h) = Phi(h)^2.
  Every term of the Q series is a product, never a difference of stationary
  covariances, so the small high-order entries (Q_hh ~ h^5) keep full relative
  precision in float.  The stationary covariance is the same doubling carried
  until Phi has decayed, which is a sum of positive semidefinite terms.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <cstdint>
#include <limits>

namespace ocean_imu::kalman::carma {

// Spectral shape at kappa = 1, in Hz.  Fields an order does not use are
// ignored: CARMA(3,0) reads (f1, zeta1, f_lambda); CARMA(4,1) reads
// (f1, zeta1, f2, zeta2) and fz as the real zero B(s) = s + 2 pi fz;
// CARMA(5,2) reads every field, (fz, zeta_z) being the complex zero pair.
struct CarmaShape {
    double f1_hz;
    double zeta1;
    double f2_hz;
    double zeta2;
    double f_lambda_hz;
    double fz_hz;
    double zeta_z;
};

// Constrained log-spectral fits to JONSWAP (f_p = 0.29 Hz, gamma = 3.3) over
// 0.17..1.30 Hz; objective and constraints in tools/carma_fit.py, results in
// reports/results/carma_model/carma_fit.json.  These are research
// starting points, not validated universal tuning constants.
//
// CARMA(5,2): the preliminary constrained fit (zeta_z >= 0.20) supplied with
// the study, retained unchanged so the comparison is against the stated model.
inline constexpr CarmaShape CARMA52_JONSWAP_SHAPE{
    0.27812, 0.09715, 0.44274, 0.52299, 1.86402, 0.15105, 0.20};
// CARMA(4,1) and CARMA(3,0): the same objective and bounds.  The (4,1)
// optimum puts its numerator zero at DC (f_z0 = 0); the (3,0) optimum sits on
// the 5 Hz upper bound of its real pole, i.e. an oscillator-like f^-4 slope
// across the fitted band.
inline constexpr CarmaShape CARMA41_JONSWAP_SHAPE{
    0.306945, 0.104191, 0.815318, 0.478188, 0.0, 0.0, 0.0};
inline constexpr CarmaShape CARMA30_JONSWAP_SHAPE{
    0.386457, 0.266692, 0.0, 0.0, 5.0, 0.0, 0.0};

template<int N> struct CarmaDefaultShape;
template<> struct CarmaDefaultShape<3> { static constexpr CarmaShape value = CARMA30_JONSWAP_SHAPE; };
template<> struct CarmaDefaultShape<4> { static constexpr CarmaShape value = CARMA41_JONSWAP_SHAPE; };
template<> struct CarmaDefaultShape<5> { static constexpr CarmaShape value = CARMA52_JONSWAP_SHAPE; };

namespace detail {

template<typename T, int M>
using Mat = Eigen::Matrix<T, M, M>;

// Scalar-loop products keep the fixed-size template instantiations small on
// the ESP32-S3 toolchain (see the loop comments in Kalman3D_Wave_OU_III).
template<typename T, int M>
inline void mul(const Mat<T,M>& X, const Mat<T,M>& Y, Mat<T,M>& Z) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            T s = T(0);
            for (int k = 0; k < M; ++k) s += X(i,k) * Y(k,j);
            Z(i,j) = s;
        }
    }
}

// Z = X Y X^T for symmetric Y; Z is written symmetric.
template<typename T, int M>
inline void congruence(const Mat<T,M>& X, const Mat<T,M>& Y, Mat<T,M>& Z) {
    Mat<T,M> XY;
    mul<T,M>(X, Y, XY);
    for (int i = 0; i < M; ++i) {
        for (int j = i; j < M; ++j) {
            T s = T(0);
            for (int k = 0; k < M; ++k) s += XY(i,k) * X(j,k);
            Z(i,j) = s;
            Z(j,i) = s;
        }
    }
}

template<typename T, int M>
inline T max_abs(const Mat<T,M>& X) {
    T m = T(0);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) m = std::max(m, std::abs(X(i,j)));
    return m;
}

template<typename T, int M>
inline T inf_norm(const Mat<T,M>& X) {
    T m = T(0);
    for (int i = 0; i < M; ++i) {
        T r = T(0);
        for (int j = 0; j < M; ++j) r += std::abs(X(i,j));
        m = std::max(m, r);
    }
    return m;
}

}  // namespace detail

// Exact discretization of dx = A x dt + dB, Cov(dB) = W dt, over step h >= 0.
// Returns false (Phi = I, Q = 0) for a non-finite or negative step or input.
// (Arduino defines F() and TWO_PI as macros, hence Phi and k2Pi.)
template<typename T, int M>
inline bool discretize_lti(const Eigen::Matrix<T,M,M>& A,
                           const Eigen::Matrix<T,M,M>& W,
                           T h,
                           Eigen::Matrix<T,M,M>& Phi,
                           Eigen::Matrix<T,M,M>& Q)
{
    using Mat = detail::Mat<T,M>;
    Phi.setIdentity();
    Q.setZero();
    if (!(h > T(0)) || !std::isfinite(h)) return h == T(0);

    const T norm_a = detail::inf_norm<T,M>(A);
    if (!std::isfinite(norm_a)) return false;

    // h_s = h / 2^s with ||A h_s||_inf <= 1/4.
    int s = 0;
    T hs = h;
    while (norm_a * hs > T(0.25) && s < 60) {
        hs *= T(0.5);
        ++s;
    }

    // Structural distance d_i of state i from the noise input (number of
    // couplings in A); states the noise cannot reach get 0.  In the scaled
    // coordinates x = D z, D = diag(h_s^d_i), every entry of Q and F is of
    // the order of its leading term, so the series below neither loses the
    // small entries (Q_hh ~ h^5) to the large ones nor produces subnormal
    // intermediate values.
    int d[M];
    for (int i = 0; i < M; ++i) d[i] = (W(i,i) != T(0)) ? 0 : 4 * M;
    for (int pass = 0; pass < M; ++pass)
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < M; ++j)
                if (A(i,j) != T(0) && d[j] + 1 < d[i]) d[i] = d[j] + 1;
    T hp[2 * M + 1];
    hp[0] = T(1);
    for (int i = 1; i <= 2 * M; ++i) hp[i] = hp[i - 1] * hs;
    for (int i = 0; i < M; ++i) if (d[i] >= 4 * M) d[i] = 0;
    // A_s = D^-1 A D: A_ij h_s^(d_j - d_i), exponent >= -1 on every coupling.
    Mat As;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            const int e = d[j] - d[i];
            As(i,j) = (e >= 0) ? A(i,j) * hp[e] : A(i,j) / hp[-e];
        }

    // F_s(h_s) = sum (A_s h_s)^k / k!,   Q_s(h_s) = sum X_k,
    // X_0 = h_s W,  X_k = h_s/(k+1) (A_s X_{k-1} + X_{k-1} A_s^T).
    // The leading term of an entry m couplings from the noise appears at
    // k = 2m <= 2(M - 1) in Q and at k <= M - 1 in F; past it every term
    // shrinks at least like (1/4)^j / j!, so K_TAIL more terms reach the
    // working epsilon entry by entry.
    constexpr int K_TAIL = (std::numeric_limits<T>::digits > 30) ? 13 : 7;
    constexpr int K_Q = 2 * (M - 1) + K_TAIL;
    constexpr int K_F = (M - 1) + K_TAIL;
    Mat Y = Mat::Identity();
    Mat X = W * hs;
    Mat Fs = Y;
    Mat Qs = X;
    Mat AY, AX;
    for (int k = 1; k <= K_Q; ++k) {
        if (k <= K_F) {
            detail::mul<T,M>(As, Y, AY);
            Y = AY * (hs / T(k));
            Fs += Y;
        }
        detail::mul<T,M>(As, X, AX);
        for (int i = 0; i < M; ++i) {
            for (int j = i; j < M; ++j) {
                const T v = (AX(i,j) + AX(j,i)) * (hs / T(k + 1));
                X(i,j) = v;
                X(j,i) = v;
            }
        }
        Qs += X;
    }
    // Back to x: F = D F_s D^-1, Q = D Q_s D.
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            const int e = d[i] - d[j];
            Phi(i,j) = (e >= 0) ? Fs(i,j) * hp[e] : Fs(i,j) / hp[-e];
            Q(i,j) = Qs(i,j) * hp[d[i] + d[j]];
        }

    // Exact doublings.
    Mat FQF, FF;
    for (int i = 0; i < s; ++i) {
        detail::congruence<T,M>(Phi, Q, FQF);
        Q += FQF;
        detail::mul<T,M>(Phi, Phi, FF);
        Phi = FF;
    }
    for (int i = 0; i < M; ++i) {
        for (int j = i + 1; j < M; ++j) {
            const T v = T(0.5) * (Q(i,j) + Q(j,i));
            Q(i,j) = v;
            Q(j,i) = v;
        }
    }
    return Phi.allFinite() && Q.allFinite();
}

// Stationary covariance P of a Hurwitz system, A P + P A^T + W = 0, as the
// limit of Q(h 2^k): every doubling adds the positive semidefinite term
// F Q F^T, and the remainder after reaching horizon T is F(T) P F(T)^T.
// Returns false if F has not decayed after the doubling budget (A not
// Hurwitz to working precision).
template<typename T, int M>
inline bool stationary_covariance_doubling(const Eigen::Matrix<T,M,M>& A,
                                           const Eigen::Matrix<T,M,M>& W,
                                           Eigen::Matrix<T,M,M>& P)
{
    using Mat = detail::Mat<T,M>;
    const T norm_a = detail::inf_norm<T,M>(A);
    if (!(norm_a > T(0)) || !std::isfinite(norm_a)) {
        P.setZero();
        return false;
    }
    Mat Phi;
    if (!discretize_lti<T,M>(A, W, T(0.25) / norm_a, Phi, P)) return false;
    const T tol = T(4) * std::numeric_limits<T>::epsilon();
    Mat FPF, FF;
    for (int it = 0; it < 96; ++it) {
        detail::congruence<T,M>(Phi, P, FPF);
        P += FPF;
        detail::mul<T,M>(Phi, Phi, FF);
        Phi = FF;
        if (detail::max_abs<T,M>(Phi) <= tol) return P.allFinite();
    }
    return false;
}

// Continuous dynamics in [h, v, a, latent] coordinates for unit noise gain
// (G = e_a), frequencies multiplied by kappa.  Derivation in the file header
// and docs/carma-three-axis.md.
template<int N>
inline void carma_dynamics(const CarmaShape& sh, double kappa,
                           Eigen::Matrix<double,N,N>& A)
{
    static_assert(N >= 3 && N <= 5, "CARMA(3,0), (4,1) and (5,2) only");
    constexpr double k2Pi = 6.283185307179586476925286766559;
    const double w1 = k2Pi * kappa * sh.f1_hz;
    const double z1 = sh.zeta1;
    A.setZero();
    A(0,1) = 1.0;   // h' = v
    A(1,2) = 1.0;   // v' = a
    if constexpr (N == 3) {
        // a' = -lambda w1^2 h - (2 lambda zeta1 w1 + w1^2) v - (lambda + 2 zeta1 w1) a
        const double lam = k2Pi * kappa * sh.f_lambda_hz;
        A(2,0) = -lam * w1 * w1;
        A(2,1) = -(2.0 * lam * z1 * w1 + w1 * w1);
        A(2,2) = -(lam + 2.0 * z1 * w1);
    } else if constexpr (N == 4) {
        // Cascade r'' + 2 zeta2 w2 r' + w2^2 r = xi, h'' + 2 zeta1 w1 h' +
        // w1^2 h = r' + w_z0 r.  Latent l = r; r' is eliminated through a.
        const double w2 = k2Pi * kappa * sh.f2_hz;
        const double z2 = sh.zeta2;
        const double wz0 = k2Pi * kappa * sh.fz_hz;
        const double d = wz0 - 2.0 * z2 * w2;
        A(2,0) = d * w1 * w1;
        A(2,1) = 2.0 * d * z1 * w1 - w1 * w1;
        A(2,2) = d - 2.0 * z1 * w1;
        A(2,3) = -w2 * w2 - d * wz0;
        A(3,0) = w1 * w1;
        A(3,1) = 2.0 * z1 * w1;
        A(3,2) = 1.0;
        A(3,3) = -wz0;
    } else {
        // Cascade u' = -lambda u + xi, r'' + 2 zeta2 w2 r' + w2^2 r = u,
        // h'' + 2 zeta1 w1 h' + w1^2 h = u + beta1 r' + beta0 r with
        // B/L2 = 1 + (beta1 s + beta0)/L2.  Latent (l1, l2) = (r, r'); u is
        // eliminated through a.
        const double w2 = k2Pi * kappa * sh.f2_hz;
        const double z2 = sh.zeta2;
        const double lam = k2Pi * kappa * sh.f_lambda_hz;
        const double wz = k2Pi * kappa * sh.fz_hz;
        const double zz = sh.zeta_z;
        const double b1 = 2.0 * zz * wz - 2.0 * z2 * w2;
        const double b0 = wz * wz - w2 * w2;
        const double gm = b1 - lam;
        A(2,0) = gm * w1 * w1;
        A(2,1) = -w1 * w1 + 2.0 * gm * z1 * w1;
        A(2,2) = gm - 2.0 * z1 * w1;
        A(2,3) = -b1 * w2 * w2 - gm * b0;
        A(2,4) = b0 - 2.0 * b1 * z2 * w2 - gm * b1;
        A(3,4) = 1.0;
        A(4,0) = w1 * w1;
        A(4,1) = 2.0 * z1 * w1;
        A(4,2) = 1.0;
        A(4,3) = -w2 * w2 - b0;
        A(4,4) = -2.0 * z2 * w2 - b1;
    }
}

// |B(iw)|^2 / |P(iw)|^2 evaluated directly from the polynomial factors, for
// unit gain at scale kappa.  Independent of carma_dynamics(); used to verify it.
template<int N>
inline double carma_unit_psd(const CarmaShape& sh, double kappa, double f_hz) {
    constexpr double k2Pi = 6.283185307179586476925286766559;
    const double w = k2Pi * f_hz;
    auto quad2 = [w](double fk, double zk) {
        const double wk = fk;
        const double re = wk * wk - w * w;
        const double im = 2.0 * zk * wk * w;
        return re * re + im * im;
    };
    const double w1 = k2Pi * kappa * sh.f1_hz;
    double num = 1.0;
    double den = quad2(w1, sh.zeta1);
    if constexpr (N == 3) {
        const double lam = k2Pi * kappa * sh.f_lambda_hz;
        den *= lam * lam + w * w;
    } else if constexpr (N == 4) {
        const double wz0 = k2Pi * kappa * sh.fz_hz;
        num = wz0 * wz0 + w * w;
        den *= quad2(k2Pi * kappa * sh.f2_hz, sh.zeta2);
    } else {
        const double lam = k2Pi * kappa * sh.f_lambda_hz;
        num = quad2(k2Pi * kappa * sh.fz_hz, sh.zeta_z);
        den *= quad2(k2Pi * kappa * sh.f2_hz, sh.zeta2) * (lam * lam + w * w);
    }
    return num / den;
}

// One configured CARMA vertical prior: a fixed shape, a frequency scale and a
// stationary acceleration standard deviation.
template<typename T, int N>
class CarmaVerticalModel {
public:
    static_assert(N >= 3 && N <= 5, "CARMA(3,0), (4,1) and (5,2) only");
    using MatN = Eigen::Matrix<T,N,N>;
    static constexpr int IDX_H = 0;
    static constexpr int IDX_V = 1;
    static constexpr int IDX_A = 2;
    static constexpr int ORDER = N;

    explicit CarmaVerticalModel(const CarmaShape& shape = CarmaDefaultShape<N>::value) {
        set_shape(shape);
    }

    // Changes the spectral shape and recomputes its kappa = 1 moments.  This
    // solves one Lyapunov equation in double precision, so it belongs to
    // construction or reconfiguration, never to the sample loop.
    bool set_shape(const CarmaShape& shape) {
        shape_ = shape;
        Eigen::Matrix<double,N,N> A1, W1, P1;
        carma_dynamics<N>(shape_, 1.0, A1);
        W1.setZero();
        W1(IDX_A, IDX_A) = 1.0;
        if (!stationary_covariance_doubling<double,N>(A1, W1, P1)) {
            shape_valid_ = false;
            return false;
        }
        base_var_h_ = P1(IDX_H, IDX_H);
        base_var_v_ = P1(IDX_V, IDX_V);
        base_var_a_ = P1(IDX_A, IDX_A);
        shape_valid_ = base_var_h_ > 0.0 && base_var_v_ > 0.0 && base_var_a_ > 0.0;
        configured_ = false;
        return shape_valid_ && configure(kappa_, accel_std_);
    }

    // kappa multiplies every pole and zero frequency of the shape; accel_std is
    // the requested stationary standard deviation of a [m/s^2].
    bool configure(T kappa, T accel_std) {
        if (!shape_valid_) return false;
        if (!(kappa > T(0)) || !std::isfinite(kappa)) return false;
        if (!(accel_std > T(0)) || !std::isfinite(accel_std)) return false;
        Eigen::Matrix<double,N,N> Ad;
        carma_dynamics<N>(shape_, static_cast<double>(kappa), Ad);
        A_ = Ad.template cast<T>();
        const double gain2 = static_cast<double>(accel_std) * static_cast<double>(accel_std)
                             * static_cast<double>(kappa) / base_var_a_;
        noise_gain2_ = static_cast<T>(gain2);
        W_.setZero();
        W_(IDX_A, IDX_A) = noise_gain2_;
        kappa_ = kappa;
        accel_std_ = accel_std;
        configured_ = true;
        ++generation_;
        return true;
    }

    bool discretize(T dt, MatN& Phi, MatN& Q) const {
        return discretize_lti<T,N>(A_, W_, dt, Phi, Q);
    }

    bool stationary_covariance(MatN& P) const {
        return stationary_covariance_doubling<T,N>(A_, W_, P);
    }

    const MatN& A() const noexcept { return A_; }
    const MatN& W() const noexcept { return W_; }
    T noise_gain2() const noexcept { return noise_gain2_; }
    T kappa() const noexcept { return kappa_; }
    T accel_std() const noexcept { return accel_std_; }
    bool configured() const noexcept { return configured_; }
    const CarmaShape& shape() const noexcept { return shape_; }
    // Increments on every successful configure(); lets callers cache F and Q.
    uint32_t generation() const noexcept { return generation_; }

    // Displacement zero-crossing frequency of the shape at kappa = 1,
    // sqrt(m2/m0)/(2 pi) = sqrt(Var v / Var h)/(2 pi).  At scale kappa it is
    // kappa times this value.
    double base_zero_crossing_hz() const noexcept {
        return std::sqrt(base_var_v_ / base_var_h_) / 6.283185307179586476925286766559;
    }
    // Stationary variances at kappa = 1 and unit noise gain.
    double base_var_h() const noexcept { return base_var_h_; }
    double base_var_v() const noexcept { return base_var_v_; }
    double base_var_a() const noexcept { return base_var_a_; }

    // Stationary displacement standard deviation of the configured model:
    // Var_kappa(h) = g^2 Var_1(h) / kappa^5.
    T displacement_std() const noexcept {
        const double k = static_cast<double>(kappa_);
        return static_cast<T>(std::sqrt(static_cast<double>(noise_gain2_) * base_var_h_
                                        / (k * k * k * k * k)));
    }

private:
    CarmaShape shape_{};
    bool shape_valid_ = false;
    bool configured_ = false;
    double base_var_h_ = 0.0;
    double base_var_v_ = 0.0;
    double base_var_a_ = 0.0;
    T kappa_ = T(1);
    T accel_std_ = T(1);
    T noise_gain2_ = T(0);
    MatN A_ = MatN::Zero();
    MatN W_ = MatN::Zero();
    uint32_t generation_ = 0;
};

}  // namespace ocean_imu::kalman::carma
