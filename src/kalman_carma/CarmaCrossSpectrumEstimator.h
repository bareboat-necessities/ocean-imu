#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
  Released under the MIT License

  CarmaCrossSpectrumEstimator / CarmaCoupling -- EXPERIMENTAL cross-axis
  coupling source of the coupled three-axis CARMA priors
  (kalman_carma/CarmaCoupledMotion.h, docs/carma-coupled-three-axis.md).  Not
  a shipping estimator component.

  Estimator.  Complex demodulation (Bingham, Godfrey and Tukey 1967) of the
  three world-frame (NED) accelerations at the CARMA operating frequency
  f_op:

      A_i(t) = L[ a_i(t) exp(-i phi(t)) ],   phi' = 2 pi f_op,

  L two cascaded first-order low-passes with corner band_frac * f_op, then an
  exponentially weighted (horizon_s) average of A A^H.  The result S_hat is a
  Hermitian positive semidefinite 3 x 3 band cross-spectral matrix (it is an
  average of rank-one outer products): its diagonal holds the band
  acceleration powers, its off-diagonals the co-spectrum (real part) and the
  quadrature spectrum (imaginary part) around f_op.  Since a = -w^2 p on
  every axis, the ratios S_hat_ij / S_hat_zz are also the displacement
  regressions at f_op.  The input is measurement-only: specific force, the
  estimator's attitude and accelerometer-bias estimate, gravity.

  Mapping (excitation).  Given the prior standard deviations sigma_i the
  wrapper maps for each axis (the baseline amplitudes) and the response R_i of
  the model's axis i at f_op (1, except the RAO-structured vertical axis):

      c = |R_z|^2 sigma_z^2 / S_hat_zz,
      S_m = [[sigma_h^2 I, s c S_hat_hz], [s c S_hat_hz^H, |R_z|^2 sigma_z^2]],
      Sigma = R^-1 S_m R^-H.

  The diagonal is the baseline prior: horizontal band power is not used as
  wave energy, because in calm water and on small seas it is dominated by
  tilt and bias leakage (an earlier form that scaled the horizontal block by
  c raised the calm-water horizontal phantom by a quarter; it is kept as the
  Horizontal::Measured ablation).  The cross terms keep the measured
  regression of every horizontal axis on z at f_op (co- and quadrature
  spectrum), shrunk by s = gain max(0, 1 - 1/(N_eff gamma^2)), gamma^2 the
  debiased multiple coherence of z on (x, y) and N_eff the effective number of
  independent envelope samples, and scaled down if needed so that
  b^H H^-1 b <= sigma_z^2: Sigma is Hermitian positive semidefinite by
  construction and equals the baseline diagonal without coupling evidence.
  With cap_coherence = rho > 0 the cap instead keeps b and, when
  q = b^H H^-1 b exceeds rho sigma_z^2, enlarges the horizontal block only
  along the coupling direction:

      H' = H + ((lambda - 1) / q) b b^H,   lambda = q / (rho sigma_z^2),

  so that b^H H'^-1 b = rho sigma_z^2.  The prior keeps the measured
  regression E[h | z] = b z / sigma_z^2; for H = sigma_h^2 I the along-wave
  variance becomes lambda sigma_h^2 with innovation (1 - rho) lambda
  sigma_h^2, and the cross-wave variance stays sigma_h^2.  Enlarging the
  whole block by lambda instead inflates the cross-wave prior as well and
  costs horizontal accuracy where the cap binds.  Shrinking
  b at the cap understates the regression when horizontal motion exceeds the
  baseline prior (a weak heave RAO) and leaves only 0.1 % of the along-wave
  horizontal variance independent of heave; capping b at the measured
  coherence instead still understates it.  One Sigma serves the
  whole spectrum while gamma^2 is measured in the band around f_op, so rho is
  a broadband target below the band coherence (shared modes: 0.8).
  With the wrapper's equal horizontal priors the mapping is equivariant under
  rotations about the vertical (sigma_h^2 I is invariant, the cross vector
  rotates, every factor is an invariant), so the prior does not depend on how
  the waves sit relative to north.  Horizontal::Anisotropic adds the measured
  horizontal anisotropy at the baseline trace (shrunk, capped at 90 %).
  Form Real keeps only the co-spectrum (Experiment A); RealXY only the
  horizontal anisotropy (no horizontal-vertical term).  Form Numerator returns
  the Hermitian matrix like Form Hermitian; CarmaModalMotion turns its
  horizontal-vertical phases into per-axis numerator zeros (the kinematically
  admissible modal phase coupling, CarmaCoupledMotion.h).

  Sources.  Off (independent axes), Estimated (the estimator above, after a
  warm-up), External (a cross-spectral matrix supplied by the caller -- the
  benchmark's ORACLE runs feed the same estimator with the true wave
  acceleration).  Frame World holds S_hat in NED; Frame Keel rotates a frozen
  keel-frame matrix with the estimator's heading (the frame-consistency
  ablation).  freeze_at_s holds the value from that time on.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace ocean_imu::kalman::carma {

// Hermitian 3 x 3 as a symmetric real part and an antisymmetric imaginary part.
template<typename T>
struct Hermitian3 {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 re = Matrix3::Zero();
    Matrix3 im = Matrix3::Zero();

    static Hermitian3 diagonal(const Eigen::Matrix<T,3,1>& d) {
        Hermitian3 h;
        h.re = d.asDiagonal();
        return h;
    }
    std::complex<T> operator()(int i, int j) const { return {re(i,j), im(i,j)}; }
    void set(int i, int j, std::complex<T> v) {
        re(i,j) = v.real(); im(i,j) = v.imag();
        re(j,i) = v.real(); im(j,i) = -v.imag();
        if (i == j) im(i,i) = T(0);
    }
    // Hermitian positive semidefinite test by a complex LDL^H (relative tolerance).
    bool psd(T rel_tol = T(1e-6)) const {
        std::complex<T> L[3][3] = {};
        T d[3];
        const T scale = std::max({std::abs(re(0,0)), std::abs(re(1,1)), std::abs(re(2,2)), T(1e-30)});
        for (int j = 0; j < 3; ++j) {
            std::complex<T> s = (*this)(j,j);
            for (int k = 0; k < j; ++k) s -= L[j][k] * std::conj(L[j][k]) * d[k];
            d[j] = s.real();
            if (d[j] < -rel_tol * scale) return false;
            for (int i = j + 1; i < 3; ++i) {
                std::complex<T> v = (*this)(i,j);
                for (int k = 0; k < j; ++k) v -= L[i][k] * std::conj(L[j][k]) * d[k];
                L[i][j] = d[j] > rel_tol * scale ? v / d[j] : std::complex<T>(0);
                if (!(d[j] > rel_tol * scale) && std::abs(v) > std::sqrt(rel_tol) * scale) return false;
            }
        }
        return true;
    }
};

template<typename T>
class CarmaCrossSpectrumEstimator {
public:
    using Vector3 = Eigen::Matrix<T,3,1>;
    struct Config {
        T band_frac = T(0.15);     // envelope low-pass corner / f_op (screened, seeds 1001-1020)
        T horizon_s = T(240.0);    // horizon of the exponentially weighted average
        T min_freq_hz = T(0.03);
        T max_freq_hz = T(2.0);
    };

    void configure(const Config& c) { cfg_ = c; }
    const Config& config() const noexcept { return cfg_; }

    void reset() {
        phase_ = T(0);
        for (int i = 0; i < 3; ++i) { y1_[i] = y2_[i] = std::complex<T>(0); }
        acc_ = Hermitian3<T>{};
        weight_ = T(0);
        elapsed_ = T(0);
        f_last_ = T(0);
    }

    void update(const Vector3& a, T dt, T f_hz) {
        if (!(dt > T(0)) || !std::isfinite(dt) || !a.allFinite() || !std::isfinite(f_hz)) return;
        const T f = std::min(std::max(f_hz, cfg_.min_freq_hz), cfg_.max_freq_hz);
        f_last_ = f;
        constexpr T k2Pi = T(6.283185307179586476925286766559);
        phase_ += k2Pi * f * dt;
        if (phase_ > k2Pi) phase_ -= k2Pi * std::floor(phase_ / k2Pi);
        const std::complex<T> rot(std::cos(phase_), -std::sin(phase_));
        const T tau = T(1) / (k2Pi * cfg_.band_frac * f);
        const T al = -std::expm1(-dt / tau);
        for (int i = 0; i < 3; ++i) {
            const std::complex<T> x = a(i) * rot;
            y1_[i] += al * (x - y1_[i]);
            y2_[i] += al * (y1_[i] - y2_[i]);
        }
        const T ah = -std::expm1(-dt / std::max(cfg_.horizon_s, dt));
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) {
                const std::complex<T> v = y2_[i] * std::conj(y2_[j]);
                acc_.re(i,j) += ah * (v.real() - acc_.re(i,j));
                acc_.im(i,j) += ah * (v.imag() - acc_.im(i,j));
                acc_.re(j,i) = acc_.re(i,j);
                acc_.im(j,i) = -acc_.im(i,j);
            }
        weight_ += ah * (T(1) - weight_);
        elapsed_ += dt;
    }

    // Start-clean exponentially weighted mean of A A^H.
    Hermitian3<T> cross() const {
        Hermitian3<T> h = acc_;
        if (weight_ > T(0)) { h.re /= weight_; h.im /= weight_; }
        for (int i = 0; i < 3; ++i) h.im(i,i) = T(0);
        return h;
    }
    T elapsed_s() const noexcept { return elapsed_; }
    T frequency_hz() const noexcept { return f_last_; }
    // Independent envelope samples in the horizon (or the elapsed time).
    T effective_samples() const {
        const T span = std::min(elapsed_, cfg_.horizon_s);
        return std::max(T(1), span * cfg_.band_frac * f_last_ * T(2));
    }

private:
    Config cfg_{};
    T phase_ = T(0);
    std::complex<T> y1_[3]{}, y2_[3]{};
    Hermitian3<T> acc_{};
    T weight_ = T(0);
    T elapsed_ = T(0);
    T f_last_ = T(0);
};

struct CarmaCouplingConfig {
    enum class Source { Off, Estimated, External };
    enum class Form { Hermitian, Real, RealXY, Numerator };
    enum class Frame { World, Keel };
    // Horizontal block of the excitation: the baseline prior (default), the
    // baseline trace with the measured anisotropy, or (ablation) the measured
    // block scaled by c and floored at the prior.
    enum class Horizontal { Prior, Anisotropic, Measured };
    Source source = Source::Off;
    Form form = Form::Hermitian;
    Frame frame = Frame::World;
    Horizontal horizontal = Horizontal::Prior;
    float gain = 1.0f;                  // global factor on every off-diagonal
    float warmup_s = 120.0f;            // Estimated: no coupling before this
    float freeze_at_s = std::numeric_limits<float>::infinity();
    bool shrink = true;                 // debiased-coherence shrinkage
    // PSD cap policy: 0 scales the cross term down to applied coherence
    // 0.999; rho in (0, 0.999] keeps it and enlarges the horizontal block
    // along the coupling direction to applied coherence rho (file header).
    float cap_coherence = 0.0f;
};

template<typename T>
class CarmaCoupling {
public:
    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;
    using Config = CarmaCouplingConfig;
    using Estimator = CarmaCrossSpectrumEstimator<T>;

    void set_config(const Config& c) { cfg_ = c; }
    const Config& config() const noexcept { return cfg_; }
    Estimator& estimator() noexcept { return est_; }
    const Estimator& estimator() const noexcept { return est_; }

    // External (oracle) band cross-spectral matrix in NED, with its
    // effective sample count (0: no shrinkage).
    void set_external(const Hermitian3<T>& S, T n_eff = T(0)) {
        ext_ = S;
        ext_neff_ = n_eff;
        ext_valid_ = true;
    }

    // One accelerometer sample: world-frame acceleration, step, operating
    // frequency.  Only the Estimated source runs the estimator.
    void observe(const Vector3& a_world, T dt, T f_op_hz) {
        if (!(dt > T(0)) || !std::isfinite(dt)) return;
        t_ += dt;
        if (cfg_.source == Config::Source::Estimated) est_.update(a_world, dt, f_op_hz);
        if (!frozen_ && t_ >= T(cfg_.freeze_at_s)) freeze_now();
    }
    // Hold the value now in force from here on (World: fixed in NED; Keel:
    // fixed in the keel frame, turning with the heading).
    void freeze_now() {
        if (frozen_) return;
        Hermitian3<T> S;
        T n = T(0);
        if (current_(S, n)) {
            frozen_S_ = rotate_(S, -yaw_);     // keel coordinates at freeze time
            frozen_world_yaw_ = yaw_;
            frozen_neff_ = n;
            frozen_valid_ = true;
        }
        frozen_ = true;
    }
    bool frozen() const noexcept { return frozen_; }
    void set_heading(T yaw_rad) { if (std::isfinite(yaw_rad)) yaw_ = yaw_rad; }

    bool active() const noexcept { return cfg_.source != Config::Source::Off && cfg_.gain != 0.0f; }

    // Band cross-spectral matrix in NED currently in force; false if none.
    bool band_cross(Hermitian3<T>& S, T& n_eff) const {
        if (!active()) return false;
        if (frozen_) {
            if (!frozen_valid_) return false;
            // World: fixed in NED.  Keel: turns with the heading.
            const T ang = cfg_.frame == Config::Frame::Keel ? yaw_ : frozen_world_yaw_;
            S = rotate_(frozen_S_, ang);
            n_eff = frozen_neff_;
            return true;
        }
        return current_(S, n_eff);
    }

    // Hermitian excitation (file header) for prior stds sig and axis
    // responses R at f_op.  Without coupling evidence: diag(sig^2).
    // horizontal_gain, if given, points to two values that receive the
    // factors applied to the x and y prior variances (1 unless
    // cap_coherence enlarged them).
    Hermitian3<T> excitation(const Vector3& sig, const std::complex<T> (&R)[3],
                             T* horizontal_gain = nullptr) const {
        if (horizontal_gain) { horizontal_gain[0] = T(1); horizontal_gain[1] = T(1); }
        using C = std::complex<T>;
        Hermitian3<T> out = Hermitian3<T>::diagonal(sig.array().square().matrix());
        Hermitian3<T> S;
        T n_eff = T(0);
        if (!band_cross(S, n_eff)) return out;
        const T szz = S.re(2,2);
        if (!(szz > T(0)) || !std::isfinite(szz)) return out;
        const T rz2 = std::norm(R[2]);
        const T mzz = rz2 * sig(2) * sig(2);
        const T c = mzz / szz;
        // Cross terms b = c S_hat_hz keep the measured regression on z.
        C b0 = c * S(0,2), b1 = c * S(1,2);
        const T a = S.re(0,0), d = S.re(1,1);
        const C m = S(0,1);
        const bool shrink = cfg_.shrink && n_eff > T(0);
        const T g = static_cast<T>(cfg_.gain);
        // Debiased multiple coherence of z on (x, y): a rotation invariant.
        const T det = a * d - std::norm(m);
        T gamma2 = T(0);
        if (det > T(1e-12) * (a * d + T(1e-30))) {
            const C s0 = S(0,2), s1 = S(1,2);
            const C i0 = (d * s0 - m * s1) / det, i1 = (a * s1 - std::conj(m) * s0) / det;
            gamma2 = std::min(T(1), std::max(T(0), (std::conj(s0) * i0 + std::conj(s1) * i1).real() / szz));
        } else if (a + d > T(0)) {
            gamma2 = std::min(T(1), (std::norm(S(0,2)) + std::norm(S(1,2))) / ((a + d) * szz));
        }
        T s_hz = g;
        if (shrink) s_hz *= gamma2 > T(0) ? std::max(T(0), T(1) - T(1) / (n_eff * gamma2)) : T(0);
        const bool real_form = cfg_.form == Config::Form::Real || cfg_.form == Config::Form::RealXY;
        if (real_form) { b0 = C(b0.real(), T(0)); b1 = C(b1.real(), T(0)); }
        if (cfg_.form == Config::Form::RealXY) s_hz = T(0);
        b0 *= s_hz;
        b1 *= s_hz;
        // Horizontal block: the baseline prior (sensor noise and tilt or
        // bias leakage in the horizontal band are not wave energy), optionally
        // with the measured anisotropy at the same trace.
        const T px = std::norm(R[0]) * sig(0) * sig(0), py = std::norm(R[1]) * sig(1) * sig(1);
        T hxx = px, hyy = py;
        C hxy(T(0));
        const bool iso = std::abs(px - py) <= T(1e-6) * std::max(px, py);
        if (cfg_.horizontal == Config::Horizontal::Measured) {
            // Ablation only: the measured horizontal block scaled by c and
            // floored at the prior (the first form of this mapping; it reads
            // horizontal leakage as wave energy in calm water).
            const T ma = c * a, md = c * d;
            const C mm = c * m;
            if (iso) {
                const T h = T(0.5) * (ma + md), e = T(0.5) * (ma - md);
                const T lam_min = h - std::sqrt(e * e + std::norm(mm));
                const T delta = std::max(T(0), px - lam_min);
                hxx = ma + delta; hyy = md + delta; hxy = mm;
            } else {
                hxx = std::max(px, ma); hyy = std::max(py, md); hxy = mm;
            }
        } else if ((cfg_.horizontal == Config::Horizontal::Anisotropic || cfg_.form == Config::Form::RealXY) && iso &&
                   a + d > T(0)) {
            // Trace-preserving measured shape h (I + s A / h), |s A / h| <= 0.9.
            const T h = T(0.5) * (a + d), e = T(0.5) * (a - d);
            const T an = std::sqrt(e * e + std::norm(m)) / h;
            T s_hh = g;
            if (shrink) s_hh *= an > T(0) ? std::max(T(0), T(1) - T(1) / (n_eff * an * an)) : T(0);
            const T k = (an * s_hh > T(0.9)) ? T(0.9) / an : s_hh;
            hxx = px * (T(1) + k * e / h);
            hyy = px * (T(1) - k * e / h);
            hxy = (real_form ? C(m.real(), T(0)) : m) * (px * k / h);
        }
        // Positive semidefinite: b^H H^-1 b <= m_zz, else b scaled down
        // (a rotation-invariant factor), or H enlarged along b
        // (cap_coherence).
        const T dh = hxx * hyy - std::norm(hxy);
        if (dh > T(0)) {
            const C ib0 = (hyy * b0 - hxy * b1) / dh, ib1 = (hxx * b1 - std::conj(hxy) * b0) / dh;
            const T q = (std::conj(b0) * ib0 + std::conj(b1) * ib1).real();
            const T rho = static_cast<T>(cfg_.cap_coherence);
            if (rho > T(0)) {
                const T qcap = std::min(rho, T(0.999)) * mzz;
                if (q > qcap) {
                    // Rank-one repair along b: b^H H'^-1 b = q / (1 + al q).
                    const T al = (q / qcap - T(1)) / q;
                    const T hx0 = hxx, hy0 = hyy;
                    hxx += al * std::norm(b0);
                    hyy += al * std::norm(b1);
                    hxy += al * b0 * std::conj(b1);
                    if (horizontal_gain) { horizontal_gain[0] = hxx / hx0; horizontal_gain[1] = hyy / hy0; }
                }
            } else {
                const T qmax = T(0.999) * mzz;
                if (q > qmax) {
                    const T f = std::sqrt(qmax / q);
                    b0 *= f;
                    b1 *= f;
                }
            }
        }
        Hermitian3<T> Sm;
        Sm.re(0,0) = hxx; Sm.re(1,1) = hyy; Sm.re(2,2) = mzz;
        Sm.set(0, 1, hxy);
        Sm.set(0, 2, b0);
        Sm.set(1, 2, b1);
        // Sigma = R^-1 S_m R^-H.
        C Ri[3];
        for (int i = 0; i < 3; ++i) Ri[i] = std::abs(R[i]) > T(0) ? T(1) / R[i] : C(1);
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) out.set(i, j, Ri[i] * Sm(i,j) * std::conj(Ri[j]));
        return out;
    }

    void reset() {
        est_.reset();
        t_ = T(0);
        frozen_ = frozen_valid_ = false;
        frozen_world_yaw_ = T(0);
    }

private:
    bool current_(Hermitian3<T>& S, T& n_eff) const {
        if (cfg_.source == Config::Source::External) {
            if (!ext_valid_) return false;
            S = ext_;
            n_eff = ext_neff_;
            return true;
        }
        if (cfg_.source == Config::Source::Estimated) {
            if (est_.elapsed_s() < T(cfg_.warmup_s)) return false;
            S = est_.cross();
            n_eff = est_.effective_samples();
            return true;
        }
        return false;
    }
    // R_z(ang) S R_z(ang)^T, rotation about the vertical (NED z).
    static Hermitian3<T> rotate_(const Hermitian3<T>& S, T ang) {
        const T c = std::cos(ang), s = std::sin(ang);
        Matrix3 Rz;
        Rz << c, -s, T(0), s, c, T(0), T(0), T(0), T(1);
        Hermitian3<T> o;
        o.re = Rz * S.re * Rz.transpose();
        o.im = Rz * S.im * Rz.transpose();
        return o;
    }
    Config cfg_{};
    Estimator est_{};
    Hermitian3<T> ext_{};
    T ext_neff_ = T(0);
    bool ext_valid_ = false;
    T t_ = T(0);
    T yaw_ = T(0);
    bool frozen_ = false, frozen_valid_ = false;
    Hermitian3<T> frozen_S_{};
    T frozen_neff_ = T(0);
    // The frozen matrix is stored in keel coordinates (rotated by -yaw at
    // freeze time); a World freeze rotates it back by that same yaw.
    T frozen_world_yaw_ = T(0);
};

}  // namespace ocean_imu::kalman::carma
