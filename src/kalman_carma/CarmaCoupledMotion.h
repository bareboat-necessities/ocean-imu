#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
  Released under the MIT License

  CarmaCoupledMotion -- EXPERIMENTAL coupled translational motion blocks for
  Kalman3D_Wave_CARMA3<..., Motion> (docs/carma-coupled-three-axis.md,
  doc/kalman_carma/kalman_carma_coupled.tex).  Not a shipping component; the
  default Kalman3D_Wave_CARMA3 (Motion = void) does not include this header's
  classes in its state and is unchanged.

  Every block is in the filter's world frame (NED of the attitude estimate).
  The wrapper keeps mapping its measured operating point onto a frequency
  scale kappa and per-axis acceleration standard deviations (the baseline
  amplitudes); the coupling source (CarmaCrossSpectrumEstimator.h) supplies
  the cross-axis part.  The excitation is rebuilt at the wrapper's commits,
  never by resetting the posterior covariance: only F and Q change.

  CarmaCorrelatedMotion<T>        Experiment A, "correlated CARMA(4,1)".
      State [x block, y block, z block], each the library's C41J
      [p, v, a, l] (CarmaVerticalModel<T,4>).  All three axes share A_C and
      the unit-acceleration-variance intensity W_1, so with the real 3 x 3
      acceleration covariance Sigma_w the exact discretization is
          F = I_3 (x) F_1,   Q = Sigma_w (x) Q_1
      (axis-major ordering).  The cross-spectral matrix of p is
      Sigma_w |H(iw)|^2: real at every frequency, constant coherence.

  Experiment B (shared directional modes) is the shared-modes QMEKF of
  src/kalman_carma_shared_modes (six states per axis, 27 error states).

  CarmaModalMotion<T>              the EXACT complex-OU modal form of the
      baseline (correctness benchmark).  Partial fractions of the C41J
      transfer function, H(s) = sum_k r_k/(s - p_k) + conj, r_k = B(p_k)/P'(p_k)
      at the two upper-half-plane poles: the CARMA output is
      p = 2 Re(z_1 + z_2), dz_k = p_k z_k dt + r_k g dW, the two complex OU
      modes driven by the SAME real Brownian motion (Fasen-Hartmann and
      Scholz, arXiv:2102.11681; Brockwell, Davis and Yang 2011).  Real state per
      axis [u_1, v_1, u_2, v_2] = [Re z_1, Im z_1, Re z_2, Im z_2]: two damped
      rotations with residue input vectors g [Re r_k, Im r_k].  The residues
      make sum Re r_k = sum Re(r_k p_k) = 0, so p, v = p', a = p'' are
      noise-free linear readouts (relative degree 3 kept) -- dense rows, not
      states.  F and Q are closed form (rotation-decay blocks; complex
      exponential integrals of the common noise), with the cross-mode blocks
      that independently driven oscillators would lose.  Axes share it with a
      real excitation (Experiment A in modal coordinates; diagonal = the
      baseline exactly, by similarity).  A phase-rotated readout of these
      modes breaks the residue cancellations (relative degree 1, white
      velocity).  Per-mode complex gains beta_k on an axis driven by the
      same noise keep relative degree 3 iff sum Re(beta_k r_k) =
      sum Re(beta_k r_k p_k) = 0, i.e. iff the axis reads
      (b1 B(s) + b0)/P(s) xi with REAL b0, b1: the only kinematically
      admissible cross-axis phase is -atan(b0/(b1 w)) (w_z0 = 0), never a
      constant one.  Form Numerator uses it (directional modal coupling):
      axis i reads p + b0_i l, l = 1/P xi = sum 2 Re(z_k/(p_k + w_z0)) the
      latent readout, b0_i = -w_op tan(phase S_iz) rho_iz (|phase| <= 70 deg,
      rho the excitation coherence), excitation Re(D S D^H),
      D = diag(1/u_x, 1/u_y, 1), u_i = 1 - i b0_i/w_op: the band gain and
      phase at f_op match, at the price of a nonzero displacement density at
      DC on that axis.  Arbitrary cross-axis phase needs the complexified
      shared-mode block (src/kalman_carma_shared_modes) instead.
      configure() refuses shapes whose two poles nearly coincide (an
      ill-conditioned modal transformation); the baseline realization is the
      fallback.  Adaptation: the modal
      coordinates depend on kappa (v = c_p A(kappa) x), so a model change at a
      commit would silently re-mean velocity and acceleration; every rebuild
      therefore hands the filter the re-coordinatization
      M = T_new^-1 T_old (T: modal -> baseline coordinates), which it applies
      to the state and to the rows and columns of the covariance, so the
      baseline coordinates -- and the posterior -- are continuous.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>

#include "kalman_carma/CarmaCrossSpectrumEstimator.h"
#include "kalman_carma/CarmaVerticalModel.h"

namespace ocean_imu::kalman::carma {

namespace coupled_detail {

// P block (ra, rb) <- Fa P_ab Fb^T, mirrored; diagonal blocks upper triangle.
template<typename T, int NA, int NB, class Mat>
inline void congruence_block(Mat& P, int ra, int rb, const Eigen::Matrix<T,NA,NA>& Fa,
                             const Eigen::Matrix<T,NB,NB>& Fb, bool diag) {
    T X[NA][NB];
    for (int i = 0; i < NA; ++i)
        for (int j = 0; j < NB; ++j) {
            T acc = T(0);
            for (int k = 0; k < NA; ++k) acc += Fa(i,k) * P(ra + k, rb + j);
            X[i][j] = acc;
        }
    for (int i = 0; i < NA; ++i)
        for (int j = diag ? i : 0; j < NB; ++j) {
            T acc = T(0);
            for (int k = 0; k < NB; ++k) acc += X[i][k] * Fb(j,k);
            P(ra + i, rb + j) = acc;
            P(rb + j, ra + i) = acc;
        }
}

}  // namespace coupled_detail

// Experiment A (file header).
template<typename T>
class CarmaCorrelatedMotion {
public:
    static constexpr bool kUnitReadout = true;   // p, v are states
    static constexpr int NB = 4;
    static constexpr int NL = 3 * NB;
    static constexpr int MAX_ACC_TERMS = 1;
    using Vector3 = Eigen::Matrix<T,3,1>;
    using MatB = Eigen::Matrix<T,NB,NB>;
    using Model = CarmaVerticalModel<T,4>;

    static constexpr int pos_local(int ax) { return NB * ax; }
    static constexpr int vel_local(int ax) { return NB * ax + 1; }
    static constexpr int acc_local(int ax) { return NB * ax + 2; }
    int accel_terms(int ax, int (&idx)[4], T (&w)[4]) const {
        idx[0] = acc_local(ax);
        w[0] = T(1);
        return 1;
    }

    explicit CarmaCorrelatedMotion(const CarmaShape& shape = CARMA41_JONSWAP_SHAPE)
        : vmodel_(shape), hmodel_(shape) {
        vmodel_.configure(T(1), T(1));
        hmodel_.configure(T(1), T(1));
        sig_ = Vector3::Ones();
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
    bool set_shape(const CarmaShape& sh) {
        if (!vmodel_.set_shape(sh) || !hmodel_.set_shape(sh)) return false;
        dirty_ = true;
        return true;
    }
    const Model& vertical_model() const noexcept { return vmodel_; }
    const Model& horizontal_model() const noexcept { return hmodel_; }
    Eigen::Matrix<T,2,2> horizontal_sigma() const { return Sigma_.re.template topLeftCorner<2,2>(); }
    Vector3 accel_std() const { return sig_; }
    const Hermitian3<T>& excitation() const noexcept { return Sigma_; }
    CarmaCoupling<T>& coupling() noexcept { return coupling_; }
    const CarmaCoupling<T>& coupling() const noexcept { return coupling_; }
    T operating_frequency_hz() const {
        return static_cast<T>(static_cast<double>(vmodel_.kappa()) * vmodel_.base_zero_crossing_hz());
    }

    // One accelerometer sample's world-frame acceleration, at the last step.
    void observe(const Vector3& a_world) { coupling_.observe(a_world, dt_, operating_frequency_hz()); }
    void set_heading(T yaw) { coupling_.set_heading(yaw); }
    T step() const noexcept { return dt_; }

    // Rebuild at a commit, then the step discretization.
    void refresh(T Ts) {
        if (Ts > T(0) && std::isfinite(Ts)) dt_ = Ts;
        if (dirty_) rebuild_();
        if (!disc_.advance(gen_, Ts)) disc_.refresh(A_, W_, gen_, Ts);
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
        const MatB& Q = disc_.Q();
        for (int a = 0; a < 3; ++a)
            for (int b = a; b < 3; ++b)
                coupled_detail::congruence_block<T,NB,NB>(P, off + a * NB, off + b * NB, Phi, Phi, a == b);
        add_structured_(P, off, Q);
    }

    template<class Mat>
    bool stationary(Mat& P, int off) const {
        MatB P1;
        if (!stationary_covariance_doubling<T,NB>(A_, W_, P1)) return false;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                for (int i = 0; i < NB; ++i)
                    for (int j = 0; j < NB; ++j) P(off + a * NB + i, off + b * NB + j) = T(0);
        add_structured_(P, off, P1);
        return true;
    }

    const MatB& block_A() const noexcept { return A_; }
    const MatB& block_W() const noexcept { return W_; }
    const MatB& block_F() const noexcept { return disc_.Phi(); }
    const MatB& block_Q() const noexcept { return disc_.Q(); }

private:
    // P_ab += Re(S_ab) M over the upper block triangle, mirrored.
    template<class Mat>
    void add_structured_(Mat& P, int off, const MatB& M) const {
        for (int a = 0; a < 3; ++a)
            for (int b = a; b < 3; ++b) {
                const T sr = Sigma_.re(a,b);
                if (sr == T(0)) continue;
                for (int i = 0; i < NB; ++i)
                    for (int j = (a == b ? i : 0); j < NB; ++j) {
                        const T v = sr * M(i,j);
                        const int r = off + a * NB + i, c = off + b * NB + j;
                        P(r, c) += v;
                        if (r != c) P(c, r) += v;
                    }
            }
    }

    void rebuild_() {
        const T kappa = vmodel_.kappa();
        {
            Model unit(vmodel_.shape());
            unit.configure(kappa, T(1));
            A_ = unit.A();
        }
        // Unit-acceleration-variance gain of the real C41J law at kappa.
        Model unit(vmodel_.shape());
        unit.configure(kappa, T(1));
        const T g2 = unit.noise_gain2();
        W_.setZero();
        W_(2,2) = g2;
        const std::complex<T> R[3] = {T(1), T(1), T(1)};
        Hermitian3<T> S = coupling_.excitation(sig_, R);
        S.im.setZero();
        Sigma_ = S;
        ++gen_;
        dirty_ = false;
    }

    Model vmodel_, hmodel_;
    Vector3 sig_ = Vector3::Ones();
    Hermitian3<T> Sigma_{};
    MatB A_ = MatB::Zero(), W_ = MatB::Zero();
    StepDiscretization<T,NB> disc_;
    uint32_t gen_ = 0;
    T dt_ = T(0.005);
    bool dirty_ = true;
    CarmaCoupling<T> coupling_;
};

// Exact complex-OU modal form of the baseline C41J (file header).
template<typename T>
class CarmaModalMotion {
public:
    static constexpr bool kUnitReadout = false;   // p, v, a are dense readouts
    static constexpr int NB = 4, NL = 12;
    static constexpr int MAX_ACC_TERMS = 4;
    static constexpr int MAX_READ_TERMS = 4;
    using Vector3 = Eigen::Matrix<T,3,1>;
    using MatB = Eigen::Matrix<T,NB,NB>;
    using Model = CarmaVerticalModel<T,4>;
    using Cd = std::complex<double>;

    static constexpr int pos_local(int) { return -1; }
    static constexpr int vel_local(int) { return -1; }
    // p = c_p x, v = c_p A x, a = c_p A^2 x on the axis block.
    int pos_terms(int ax, int (&idx)[4], T (&w)[4]) const { return terms_(ax, cpa_[ax], idx, w); }
    int vel_terms(int ax, int (&idx)[4], T (&w)[4]) const { return terms_(ax, cva_[ax], idx, w); }
    int accel_terms(int ax, int (&idx)[4], T (&w)[4]) const { return terms_(ax, caa_[ax], idx, w); }

    explicit CarmaModalMotion(const CarmaShape& shape = CARMA41_JONSWAP_SHAPE)
        : vmodel_(shape), hmodel_(shape) {
        vmodel_.configure(T(1), T(1));
        hmodel_.configure(T(1), T(1));
        rebuild_();
    }
    bool set_vertical(T kappa, T accel_std) {
        if (!vmodel_.configure(kappa, accel_std)) return false;
        sig_(2) = accel_std;
        dirty_ = true;
        return true;
    }
    bool set_horizontal(T kappa, T sx, T sy, T) {
        if (!(sx > T(0)) || !(sy > T(0)) || !std::isfinite(sx) || !std::isfinite(sy)) return false;
        if (!hmodel_.configure(kappa, T(1))) return false;
        sig_(0) = sx; sig_(1) = sy;
        dirty_ = true;
        return true;
    }
    bool set_shape(const CarmaShape& sh) {
        if (!vmodel_.set_shape(sh) || !hmodel_.set_shape(sh)) return false;
        dirty_ = true;
        return true;
    }
    const Model& vertical_model() const noexcept { return vmodel_; }
    const Model& horizontal_model() const noexcept { return hmodel_; }
    Eigen::Matrix<T,2,2> horizontal_sigma() const { return Sigma_.re.template topLeftCorner<2,2>(); }
    Vector3 accel_std() const { return sig_; }
    const Hermitian3<T>& excitation() const noexcept { return Sigma_; }
    CarmaCoupling<T>& coupling() noexcept { return coupling_; }
    const CarmaCoupling<T>& coupling() const noexcept { return coupling_; }
    T operating_frequency_hz() const {
        return static_cast<T>(static_cast<double>(vmodel_.kappa()) * vmodel_.base_zero_crossing_hz());
    }
    void observe(const Vector3& a_world) { coupling_.observe(a_world, dt_, operating_frequency_hz()); }
    void set_heading(T yaw) { coupling_.set_heading(yaw); }
    T step() const noexcept { return dt_; }
    bool valid() const noexcept { return valid_; }

    void refresh(T Ts) {
        if (Ts > T(0) && std::isfinite(Ts)) dt_ = Ts;
        if (dirty_) rebuild_();
        if (gen_ == disc_gen_ && dt_ == disc_h_) return;
        discretize_(static_cast<double>(dt_), F_, Q_);
        disc_gen_ = gen_;
        disc_h_ = dt_;
    }

    void apply_F(const T (&in)[NL], T (&out)[NL]) const {
        for (int ax = 0; ax < 3; ++ax)
            for (int i = 0; i < NB; ++i) {
                T acc = T(0);
                for (int j = 0; j < NB; ++j) acc += F_(i,j) * in[ax * NB + j];
                out[ax * NB + i] = acc;
            }
    }
    template<class Mat>
    void propagate(Mat& P, int off) const {
        for (int a = 0; a < 3; ++a)
            for (int b = a; b < 3; ++b)
                coupled_detail::congruence_block<T,NB,NB>(P, off + a * NB, off + b * NB, F_, F_, a == b);
        add_kron_(P, off, Q_);
    }
    template<class Mat>
    bool stationary(Mat& P, int off) const {
        MatB F1, P1;
        discretize_(std::numeric_limits<double>::infinity(), F1, P1);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                for (int i = 0; i < NB; ++i)
                    for (int j = 0; j < NB; ++j) P(off + a * NB + i, off + b * NB + j) = T(0);
        add_kron_(P, off, P1);
        return P1.allFinite();
    }

    // Diagnostics (tests, analysis): modal dynamics, input, readouts at the
    // current kappa, and the modal transformation T (x_baseline = T x_modal)
    // with its condition number.
    const MatB& block_A() const noexcept { return A_; }
    const Eigen::Matrix<T,NB,1>& block_G() const noexcept { return G_; }
    T noise_gain2() const noexcept { return g2_; }
    const MatB& block_F() const noexcept { return F_; }
    const MatB& block_Q() const noexcept { return Q_; }
    const Eigen::Matrix<T,1,NB>& readout_p() const noexcept { return cp_; }
    const Eigen::Matrix<T,1,NB>& readout_v() const noexcept { return cv_; }
    const Eigen::Matrix<T,1,NB>& readout_a() const noexcept { return ca_; }
    double transform_condition() const noexcept { return cond_T_; }
    // Per-axis numerator zero b0 (rad/s; nonzero only with Form Numerator):
    // the axis displacement is (B(s) + b0) / P(s) xi_axis.
    T numerator_zero(int ax) const noexcept { return b0_[ax]; }
    const Eigen::Matrix<T,1,NB>& readout_p(int ax) const noexcept { return cpa_[ax]; }
    const Eigen::Matrix<T,1,NB>& readout_a(int ax) const noexcept { return caa_[ax]; }
    // Re-coordinatization after a model change (file header): x_new = M x_old
    // on every axis block; false when none is pending.
    bool take_recoordinatization(MatB& M) {
        if (!pending_) return false;
        M = M_;
        pending_ = false;
        return true;
    }
    Cd pole(int k) const noexcept { return mu_[k]; }
    Cd residue(int k) const noexcept { return res_[k]; }

private:
    int terms_(int ax, const Eigen::Matrix<T,1,NB>& c, int (&idx)[4], T (&w)[4]) const {
        int n = 0;
        for (int j = 0; j < NB; ++j)
            if (c(j) != T(0)) { idx[n] = ax * NB + j; w[n] = c(j); ++n; }
        return n;
    }
    template<class Mat>
    void add_kron_(Mat& P, int off, const MatB& M) const {
        for (int a = 0; a < 3; ++a)
            for (int b = a; b < 3; ++b) {
                const T sr = Sigma_.re(a,b);
                if (sr == T(0)) continue;
                for (int i = 0; i < NB; ++i)
                    for (int j = (a == b ? i : 0); j < NB; ++j) {
                        const int r = off + a * NB + i, c = off + b * NB + j;
                        P(r, c) += sr * M(i,j);
                        if (r != c) P(c, r) += sr * M(i,j);
                    }
            }
    }
    // Closed-form F(h), Q(h) (h = inf: the stationary covariance) of the two
    // complex modes driven by one real noise of intensity g2 (unit
    // acceleration variance), in double.
    void discretize_(double h, MatB& Fo, MatB& Qo) const {
        Fo.setZero();
        Qo.setZero();
        const bool inf = !std::isfinite(h);
        for (int k = 0; k < 2; ++k) {
            const Cd e = inf ? Cd(0.0) : std::exp(mu_[k] * h);
            Fo(2 * k, 2 * k) = static_cast<T>(e.real());
            Fo(2 * k, 2 * k + 1) = static_cast<T>(-e.imag());
            Fo(2 * k + 1, 2 * k) = static_cast<T>(e.imag());
            Fo(2 * k + 1, 2 * k + 1) = static_cast<T>(e.real());
        }
        // E[z_k conj z_l] = g2 r_k conj r_l (e^{(mu_k + conj mu_l) h} - 1)/(mu_k + conj mu_l),
        // E[z_k z_l]      = g2 r_k r_l      (e^{(mu_k + mu_l) h} - 1)/(mu_k + mu_l).
        auto integ = [&](Cd x) {
            if (inf) return Cd(-1.0) / x;
            const Cd xh = x * h;
            if (std::abs(xh) < 1e-3) return h * (1.0 + xh * (0.5 + xh * (1.0 / 6.0 + xh / 24.0)));
            return (std::exp(xh) - 1.0) / x;
        };
        const double g2 = static_cast<double>(g2_);
        for (int k = 0; k < 2; ++k)
            for (int l = 0; l < 2; ++l) {
                const Cd A = g2 * res_[k] * std::conj(res_[l]) * integ(mu_[k] + std::conj(mu_[l]));
                const Cd B = g2 * res_[k] * res_[l] * integ(mu_[k] + mu_[l]);
                Qo(2 * k, 2 * l) = static_cast<T>(0.5 * (A + B).real());              // u_k u_l
                Qo(2 * k + 1, 2 * l + 1) = static_cast<T>(0.5 * (A - B).real());      // v_k v_l
                Qo(2 * k, 2 * l + 1) = static_cast<T>(0.5 * (B.imag() - A.imag()));   // u_k v_l
                Qo(2 * k + 1, 2 * l) = static_cast<T>(0.5 * (B.imag() + A.imag()));   // v_k u_l
            }
        for (int i = 0; i < NB; ++i)
            for (int j = i + 1; j < NB; ++j) {
                const T v = T(0.5) * (Qo(i,j) + Qo(j,i));
                Qo(i,j) = v; Qo(j,i) = v;
            }
    }

    void rebuild_() {
        dirty_ = false;
        const double kappa = static_cast<double>(vmodel_.kappa());
        const CarmaShape& sh = vmodel_.shape();
        constexpr double k2Pi = 6.283185307179586476925286766559;
        const double w1 = k2Pi * kappa * sh.f1_hz, z1 = sh.zeta1, w2 = k2Pi * kappa * sh.f2_hz, z2 = sh.zeta2;
        const double wz0 = k2Pi * kappa * sh.fz_hz;
        valid_ = z1 > 0.0 && z1 < 1.0 && z2 > 0.0 && z2 < 1.0;
        if (!valid_) return;
        const Cd p1(-z1 * w1, w1 * std::sqrt(1.0 - z1 * z1)), p2(-z2 * w2, w2 * std::sqrt(1.0 - z2 * z2));
        // Regular (distinct) poles only: the modal transformation degenerates
        // as they coalesce.
        if (std::abs(p1 - p2) < 1e-3 * std::max(std::abs(p1), std::abs(p2))) { valid_ = false; return; }
        mu_[0] = p1; mu_[1] = p2;
        // r_k = B(p_k) / P'(p_k), P = (s - p1)(s - conj p1)(s - p2)(s - conj p2).
        res_[0] = (p1 + wz0) / ((p1 - std::conj(p1)) * (p1 - p2) * (p1 - std::conj(p2)));
        res_[1] = (p2 + wz0) / ((p2 - std::conj(p2)) * (p2 - p1) * (p2 - std::conj(p1)));
        Eigen::Matrix<double,NB,NB> Ad = Eigen::Matrix<double,NB,NB>::Zero();
        Eigen::Matrix<double,NB,1> Gd;
        for (Eigen::Index k = 0; k < 2; ++k) {
            Ad(2 * k, 2 * k) = mu_[k].real();     Ad(2 * k, 2 * k + 1) = -mu_[k].imag();
            Ad(2 * k + 1, 2 * k) = mu_[k].imag(); Ad(2 * k + 1, 2 * k + 1) = mu_[k].real();
            Gd(2 * k) = res_[k].real();           Gd(2 * k + 1) = res_[k].imag();
        }
        Eigen::Matrix<double,1,NB> cp;
        cp << 2.0, 0.0, 2.0, 0.0;
        const Eigen::Matrix<double,1,NB> cv = cp * Ad, ca = cv * Ad;
        A_ = Ad.template cast<T>();
        G_ = Gd.template cast<T>();
        cp_ = cp.template cast<T>();
        cv_ = cv.template cast<T>();
        ca_ = ca.template cast<T>();
        // Modal -> baseline coordinates: T = O_base^-1 O_modal (observability
        // from p), condition number of the realization change.
        Eigen::Matrix<double,4,4> Ab;
        carma_dynamics<4>(sh, kappa, Ab);
        Eigen::Matrix<double,4,4> Ob, Om;
        Eigen::Matrix<double,1,4> e0 = Eigen::Matrix<double,1,4>::Zero(); e0(0) = 1.0;
        Eigen::Matrix<double,1,4> rb = e0, rm = cp;
        for (int i = 0; i < 4; ++i) { Ob.row(i) = rb; Om.row(i) = rm; rb = rb * Ab; rm = rm * Ad; }
        const Eigen::Matrix<double,4,4> Tm = Ob.inverse() * Om;
        Eigen::JacobiSVD<Eigen::Matrix<double,4,4>> svd(Tm);
        cond_T_ = svd.singularValues()(0) / svd.singularValues()(3);
        if (have_T_) {
            M_ = (Tm.inverse() * T_).template cast<T>();
            pending_ = true;
        }
        T_ = Tm;
        have_T_ = true;
        Model unit(vmodel_.shape());
        unit.configure(static_cast<T>(kappa), T(1));
        g2_ = unit.noise_gain2();
        const std::complex<T> R[3] = {T(1), T(1), T(1)};
        Hermitian3<T> S = coupling_.excitation(sig_, R);
        // Latent readout l = 1/P xi = sum 2 Re(z_k / (p_k + w_z0)) (the
        // residues of 1/P are r_k / B(p_k)); c_l G = c_l A G = c_l A^2 G = 0.
        Eigen::Matrix<double,1,NB> cl;
        for (Eigen::Index k = 0; k < 2; ++k) {
            const Cd q = mu_[k] + wz0;
            cl(2 * k) = 2.0 * q.real() / std::norm(q);
            cl(2 * k + 1) = 2.0 * q.imag() / std::norm(q);
        }
        double b0[3] = {0.0, 0.0, 0.0};
        if (coupling_.config().form == CarmaCouplingConfig::Form::Numerator) {
            // Kinematically admissible phase (file header): axis i reads
            // (B + b0_i)/P, whose band gain u_i = 1 - i b0_i / w_op carries the
            // phase of S_iz scaled by its coherence rho (|phase| <= 70 deg);
            // the excitation becomes Re(D S D^H), D = diag(1/u_x, 1/u_y, 1).
            const double w = k2Pi * static_cast<double>(operating_frequency_hz());
            const double tmax = 2.7474774194546216;          // tan(70 deg)
            std::complex<double> u[3] = {1.0, 1.0, 1.0};
            for (int i = 0; i < 2; ++i) {
                const std::complex<double> sz(static_cast<double>(S.re(i,2)), static_cast<double>(S.im(i,2)));
                const double den = std::sqrt(std::max(static_cast<double>(S.re(i,i)) * static_cast<double>(S.re(2,2)), 1e-300));
                if (!(std::abs(sz) > 0.0) || !(w > 0.0)) continue;
                const double rho = std::min(1.0, std::abs(sz) / den);
                const double t = sz.real() != 0.0 ? sz.imag() / sz.real() : (sz.imag() > 0.0 ? tmax : -tmax);
                b0[i] = -w * std::max(-tmax, std::min(tmax, t)) * rho;
                u[i] = std::complex<double>(1.0, -b0[i] / w);
            }
            Hermitian3<T> St;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    const std::complex<double> sij(static_cast<double>(S.re(i,j)), static_cast<double>(S.im(i,j)));
                    St.re(i,j) = static_cast<T>((sij / (u[i] * std::conj(u[j]))).real());
                }
            St.im.setZero();
            S = St;
        } else {
            S.im.setZero();                               // real excitation (Experiment A form)
        }
        for (int ax = 0; ax < 3; ++ax) {
            const Eigen::Matrix<double,1,NB> c = cp + b0[ax] * cl;
            b0_[ax] = static_cast<T>(b0[ax]);
            cpa_[ax] = c.template cast<T>();
            cva_[ax] = (c * Ad).template cast<T>();
            caa_[ax] = (c * Ad * Ad).template cast<T>();
        }
        Sigma_ = S;
        ++gen_;
    }

    Model vmodel_, hmodel_;
    Vector3 sig_ = Vector3::Ones();
    Hermitian3<T> Sigma_{};
    Cd mu_[2]{}, res_[2]{};
    MatB A_ = MatB::Zero(), F_ = MatB::Identity(), Q_ = MatB::Zero();
    Eigen::Matrix<T,NB,1> G_ = Eigen::Matrix<T,NB,1>::Zero();
    Eigen::Matrix<T,1,NB> cp_ = Eigen::Matrix<T,1,NB>::Zero(), cv_ = Eigen::Matrix<T,1,NB>::Zero(),
                          ca_ = Eigen::Matrix<T,1,NB>::Zero();
    Eigen::Matrix<T,1,NB> cpa_[3] = {Eigen::Matrix<T,1,NB>::Zero(), Eigen::Matrix<T,1,NB>::Zero(), Eigen::Matrix<T,1,NB>::Zero()};
    Eigen::Matrix<T,1,NB> cva_[3] = {Eigen::Matrix<T,1,NB>::Zero(), Eigen::Matrix<T,1,NB>::Zero(), Eigen::Matrix<T,1,NB>::Zero()};
    Eigen::Matrix<T,1,NB> caa_[3] = {Eigen::Matrix<T,1,NB>::Zero(), Eigen::Matrix<T,1,NB>::Zero(), Eigen::Matrix<T,1,NB>::Zero()};
    T b0_[3] = {T(0), T(0), T(0)};
    T g2_ = T(1);
    double cond_T_ = 0.0;
    Eigen::Matrix<double,4,4> T_ = Eigen::Matrix<double,4,4>::Identity();
    MatB M_ = MatB::Identity();
    bool have_T_ = false, pending_ = false;
    uint32_t gen_ = 0, disc_gen_ = 0xffffffffu;
    T dt_ = T(0.005), disc_h_ = T(-1);
    bool dirty_ = true, valid_ = false;
    CarmaCoupling<T> coupling_;
};

}  // namespace ocean_imu::kalman::carma
