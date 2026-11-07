/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Checks of the coupled three-axis CARMA priors (EXPERIMENTAL;
  src/kalman_carma/CarmaCoupledMotion.h, CarmaCrossSpectrumEstimator.h,
  docs/carma-coupled-three-axis.md).  Numbering follows the study brief;
  the shared directional modes (Experiment B) are tested with their
  standalone filter in tests/kalman_carma_shared_modes:

   1  reduction to independent CARMA(4,1) with zero coupling (COR off equals
      the unchanged baseline filter);
   2  positive semidefinite excitation and process covariance (random
      Hermitian inputs, every form);
   3  continuous-to-discrete covariance: Kronecker assembly against the
      double discretization of the joint system, semigroup law;
   6  cross-spectral phase: COR real;
   7  frame consistency: rotation equivariance of the excitation mapping,
      keel- and world-frame frozen coupling under heading changes;
   8  QMEKF attitude reset: attitude-motion cross-covariances of COR off and
      the baseline agree through large-tilt resets; symmetric PSD covariance;
   9  changing sea-state parameters (kappa, amplitudes, coupling) every
      commit: bounded, PSD, no covariance reset;
  11  zero / near-zero wave excitation: finite, PSD, no blow-up;
  12  no false coupling or wave energy from sensor noise alone;
  13  extended operation (3600 s, float): symmetric PSD, bounded;
  14  baseline regression: the default instantiation is the unchanged filter;
  16  exact complex-OU modal form of C41J (Fasen-Hartmann and Scholz): transfer
      function, similarity of F and Q to the baseline realization, closed-form
      stationary covariance, autocovariance, relative degree 3 through the
      residue cancellations (and their loss for a phase-rotated readout and
      for independently driven oscillators), filter equivalence with and
      without 10 Hz parameter changes (re-coordinatization), simulation
      statistics, the coalescing-pole fallback.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <type_traits>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "kalman_carma/CarmaCoupledMotion.h"
#include "kalman_carma/Kalman3D_Wave_CARMA3.h"

namespace kc = ocean_imu::kalman::carma;
using ocean_imu::kalman::Kalman3D_Wave_CARMA3;

namespace {

int g_failures = 0, g_checks = 0;
void check(bool ok, const char* what, double v = NAN, double lim = NAN) {
    ++g_checks;
    if (!ok) { ++g_failures; std::printf("FAIL: %s (value=%.6g limit=%.6g)\n", what, v, lim); }
    else std::printf("  ok: %s (value=%.4g limit=%.4g)\n", what, v, lim);
}

constexpr double PI = 3.14159265358979323846;
constexpr double G = 9.80665;

template<typename T, class M> using KF = Kalman3D_Wave_CARMA3<T, 4, false, 4, false, true, true, M>;
template<typename T> using Base = Kalman3D_Wave_CARMA3<T, 4, false, 4, false>;
template<typename T> using COR = KF<T, kc::CarmaCorrelatedMotion<T>>;
template<typename T> using MOD = KF<T, kc::CarmaModalMotion<T>>;
using Cfg = kc::CarmaCouplingConfig;

template<typename T> using V3 = Eigen::Matrix<T,3,1>;

// Shipping-simulator noise configuration (tests/kalman_carma/carma3/Bench3.h).
constexpr double ACC_SIGMA = 1.51e-3 * 9.80665;
constexpr double SIGMA_A = 2.8 * ACC_SIGMA * 0.71;
constexpr double SIGMA_G = 2.0 * 0.00157 * 0.01;

template<class F>
F make_filter() {
    using T = typename F::Vector3::Scalar;
    return F(V3<T>::Constant(T(SIGMA_A)), V3<T>::Constant(T(SIGMA_G)), V3<T>::Constant(T(7.7)));
}

template<class F, typename T>
void configure(F& f, T kappa, const std::type_identity_t<V3<T>>& sig) {
    f.set_vertical_carma(kappa, sig(2));
    f.set_horizontal_carma(kappa, sig(0), sig(1), T(0));
    f.reset_aw_covariance_to_stationary();
}

template<typename T>
kc::Hermitian3<T> random_psd(std::mt19937_64& rng, double offscale = 1.0) {
    std::normal_distribution<double> n01(0.0, 1.0);
    Eigen::Matrix3cd X;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) X(i,j) = {n01(rng), n01(rng)};
    Eigen::Matrix3cd S = X * X.adjoint() / 3.0;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) if (i != j) S(i,j) *= offscale;
    kc::Hermitian3<T> h;
    h.re = S.real().cast<T>();
    h.im = S.imag().cast<T>();
    return h;
}

template<typename T>
double min_eig_real_rep(const kc::Hermitian3<T>& h) {
    Eigen::Matrix<double,6,6> M;
    M << h.re.template cast<double>(), -h.im.template cast<double>(), h.im.template cast<double>(), h.re.template cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> es(M);
    return es.eigenvalues().minCoeff() / std::max(1e-300, es.eigenvalues().maxCoeff());
}

template<class Mat>
double min_corr_eig(const Mat& P) {
    const Eigen::MatrixXd D = P.template cast<double>();
    const Eigen::VectorXd d = D.diagonal().cwiseMax(1e-300).cwiseSqrt().cwiseInverse();
    const Eigen::MatrixXd C = d.asDiagonal() * D * d.asDiagonal();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C, Eigen::EigenvaluesOnly);
    return es.eigenvalues().minCoeff();
}
template<class Mat>
double asym(const Mat& P) { return (P - P.transpose()).cwiseAbs().maxCoeff(); }

// Synthetic IMU on a level, slowly yawing platform with a three-axis wave
// motion: a sum of incoherent sinusoids plus a coherent quadrature pair
// (x, y lag z by 90 deg); NED, body = level, yaw psi(t).
struct SynthImu {
    double t = 0.0;
    std::mt19937_64 rng{12345};
    std::normal_distribution<double> n01{0.0, 1.0};
    double psi_rate = 0.002;
    double amp = 1.0;
    V3<double> a_world() const {
        const double w = 2 * PI * 0.30;
        const double h = amp * 0.5 * std::cos(w * t);
        const double q = amp * 0.5 * std::sin(w * t);
        return V3<double>(-w * w * (0.4 * q + 0.05 * std::cos(2 * PI * 0.21 * t + 1.0)),
                          -w * w * (0.25 * q + 0.04 * std::sin(2 * PI * 0.43 * t)),
                          -w * w * (h + 0.05 * std::cos(2 * PI * 0.37 * t + 0.3)));
    }
    template<typename T>
    void sample(double dt, V3<T>& acc, V3<T>& gyr) {
        t += dt;
        const double psi = psi_rate * t;
        Eigen::Matrix3d Rnb = Eigen::AngleAxisd(psi, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        const V3<double> f = Rnb.transpose() * (a_world() - V3<double>(0, 0, G));
        acc = (f + ACC_SIGMA * V3<double>(n01(rng), n01(rng), n01(rng))).template cast<T>();
        gyr = (V3<double>(0, 0, psi_rate) + 0.00157 * V3<double>(n01(rng), n01(rng), n01(rng))).template cast<T>();
    }
};

template<class F>
void init_level(F& f) {
    using T = typename F::Vector3::Scalar;
    f.initialize_from_acc(V3<T>(T(0), T(0), T(-G)));
}

// ---------------------------------------------------------------------------
// Baseline and COR (coupling off) side by side; returns max |dp| / RMS p and
// the largest relative position-variance difference.
template<typename T>
void compare_cor_base(double& dp_rel, double& dvar) {
    const T kappa = T(0.85);
    const V3<T> sig(T(0.35), T(0.30), T(0.55));
    auto b = make_filter<Base<T>>();
    auto c = make_filter<COR<T>>();
    configure(b, kappa, sig); configure(c, kappa, sig);
    init_level(b); init_level(c);
    SynthImu imu;
    double dc = 0, pp = 0;
    dvar = 0;
    for (int k = 0; k < 60000; ++k) {
        V3<T> acc, gyr;
        imu.sample(0.005, acc, gyr);
        b.time_update(gyr, T(0.005)); c.time_update(gyr, T(0.005));
        b.measurement_update_acc_only(acc); c.measurement_update_acc_only(acc);
        dc = std::max(dc, double((c.get_position() - b.get_position()).cwiseAbs().maxCoeff()));
        pp += b.get_position().squaredNorm();
        for (int ax = 0; ax < 3; ++ax) {
            const double vb = b.covariance_full()(b.pos_index(ax), b.pos_index(ax));
            const double vc = c.covariance_full()(c.pos_index(ax), c.pos_index(ax));
            dvar = std::max(dvar, std::abs(vc - vb) / vb);
        }
    }
    dp_rel = dc / std::sqrt(pp / 60000.0);
}

void test1_reduction() {
    std::printf("[1] reduction to independent axes\n");
    {
        double dp, dv;
        compare_cor_base<double>(dp, dv);
        check(dp < 1e-9, "COR (coupling off) equals the baseline filter in double: max |dp| / RMS(p)", dp, 1e-9);
        check(dv < 1e-9, "COR (coupling off) position variance equals the baseline in double", dv, 1e-9);
        compare_cor_base<float>(dp, dv);
        check(dp < 1e-4, "COR (coupling off) equals the baseline filter in float (rounding): max |dp| / RMS(p)", dp, 1e-4);
        check(dv < 2e-3, "COR (coupling off) position variance equals the baseline in float (rounding)", dv, 2e-3);
    }
}

// ---------------------------------------------------------------------------
void test2_psd() {
    std::printf("[2] PSD excitation and process covariance\n");
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> us(0.05, 2.0), ua(-PI, PI);
    double worst[3] = {1, 1, 1};
    for (int trial = 0; trial < 2000; ++trial) {
        kc::Hermitian3<double> S = random_psd<double>(rng, trial % 3 == 0 ? 1.0 : 0.3);
        for (int f = 0; f < 3; ++f) {
            kc::CarmaCoupling<double> cp;
            Cfg c; c.source = Cfg::Source::External; c.form = static_cast<Cfg::Form>(f);
            c.shrink = trial % 2 == 0;
            cp.set_config(c);
            cp.set_external(S, 5.0 + trial % 50);
            const V3<double> sig(us(rng), us(rng), us(rng));
            const double ph = ua(rng);
            const std::complex<double> R[3] = {1.0, 1.0, std::polar(us(rng), ph)};
            worst[f] = std::min(worst[f], min_eig_real_rep(cp.excitation(sig, R)));
        }
    }
    check(worst[0] > -1e-9, "Hermitian excitation PSD (2000 random inputs)", worst[0], -1e-9);
    check(worst[1] > -1e-9, "real (COR) excitation PSD", worst[1], -1e-9);
    check(worst[2] > -1e-9, "horizontal-only excitation PSD", worst[2], -1e-9);
    // Assembled process covariance (float).
    {
        double me = 1;
        for (int trial = 0; trial < 20; ++trial) {
            kc::Hermitian3<float> S = random_psd<float>(rng, 0.9);
            const float kappa = static_cast<float>(us(rng) + 0.2);
            auto run = [&](auto& mo) {
                Cfg c; c.source = Cfg::Source::External; c.shrink = false;
                mo.coupling().set_config(c);
                mo.coupling().set_external(S);
                mo.set_vertical(kappa, std::sqrt(S.re(2,2)));
                mo.set_horizontal(kappa, std::sqrt(S.re(0,0)), std::sqrt(S.re(1,1)), 0.0f);
                mo.refresh(0.005f);
                constexpr int NL = std::decay_t<decltype(mo)>::NL;
                Eigen::Matrix<float,NL,NL> P = Eigen::Matrix<float,NL,NL>::Zero();
                mo.propagate(P, 0);
                me = std::min(me, min_corr_eig(P));
            };
            kc::CarmaCorrelatedMotion<float> mo;
            run(mo);
        }
        check(me > -1e-4, "COR process covariance Q PSD (float, correlation-scaled)", me, -1e-4);
    }
}

// Joint continuous system of a Kronecker motion with a real excitation:
// A = I3 (x) A_b, W_ab = S_ab W_b.
template<int NB>
void joint_kron(const Eigen::Matrix<double,NB,NB>& Ab, const Eigen::Matrix<double,NB,NB>& Wb,
                const kc::Hermitian3<double>& S,
                Eigen::Matrix<double,3*NB,3*NB>& A, Eigen::Matrix<double,3*NB,3*NB>& W) {
    A.setZero(); W.setZero();
    for (int a = 0; a < 3; ++a) {
        A.template block<NB,NB>(a * NB, a * NB) = Ab;
        for (int b = 0; b < 3; ++b) W.template block<NB,NB>(a * NB, b * NB) = S.re(a,b) * Wb;
    }
}

void test3_c2d() {
    std::printf("[3] continuous-to-discrete covariance\n");
    std::mt19937_64 rng(11);
    double worst_cor = 0, worst_semi = 0;
    for (int trial = 0; trial < 6; ++trial) {
        kc::Hermitian3<double> Sd = random_psd<double>(rng, 0.9);
        kc::Hermitian3<float> Sf; Sf.re = Sd.re.cast<float>(); Sf.im = Sd.im.cast<float>();
        const double kappa = 0.3 + 0.4 * trial;
        for (double h : {0.004, 0.005, 0.0063}) {
            auto scaled_err = [](const auto& Qf, const auto& Qd) {
                double e = 0;
                for (int i = 0; i < Qd.rows(); ++i)
                    for (int j = 0; j < Qd.cols(); ++j) {
                        const double ref = std::sqrt(std::max(Qd(i,i) * Qd(j,j), 1e-300));
                        if (ref > 0) e = std::max(e, std::abs(double(Qf(i,j)) - Qd(i,j)) / ref);
                    }
                return e;
            };
            auto setup = [&](auto& mo, const auto& S) {
                Cfg c; c.source = Cfg::Source::External; c.shrink = false;
                mo.coupling().set_config(c);
                mo.coupling().set_external(S);
                using T = std::decay_t<decltype(S.re(0,0))>;
                mo.set_vertical(T(kappa), std::sqrt(S.re(2,2)));
                mo.set_horizontal(T(kappa), std::sqrt(S.re(0,0)), std::sqrt(S.re(1,1)), T(0));
                mo.refresh(T(h));
            };
            {   // COR
                kc::CarmaCorrelatedMotion<float> mf; setup(mf, Sf);
                kc::CarmaCorrelatedMotion<double> md; setup(md, Sd);
                Eigen::Matrix<double,12,12> A, W, F, Q;
                joint_kron<4>(md.block_A(), md.block_W(), md.excitation(), A, W);
                kc::discretize_lti<double,12>(A, W, h, F, Q);
                Eigen::Matrix<float,12,12> Pf = Eigen::Matrix<float,12,12>::Zero();
                mf.propagate(Pf, 0);
                worst_cor = std::max(worst_cor, scaled_err(Pf, Q));
                // Semigroup law of the joint double discretization.
                Eigen::Matrix<double,12,12> F1, Q1, F2, Q2;
                kc::discretize_lti<double,12>(A, W, 0.6 * h, F1, Q1);
                kc::discretize_lti<double,12>(A, W, 0.4 * h, F2, Q2);
                const Eigen::Matrix<double,12,12> Qc = F2 * Q1 * F2.transpose() + Q2;
                worst_semi = std::max(worst_semi, scaled_err(Qc, Q));
            }
        }
    }
    check(worst_cor < 2e-5, "COR: float Sigma (x) Q_1 vs double joint discretization (scaled)", worst_cor, 2e-5);
    check(worst_semi < 1e-10, "joint discretization semigroup law Q(h1+h2) = F Q(h1) F' + Q(h2)", worst_semi, 1e-10);
}

// ---------------------------------------------------------------------------
template<class Mo, int NL>
std::complex<double> model_cross(const Mo& mo, double f, int i, int j) {
    Eigen::Matrix<double,NL,NL> A, W;
    joint_kron<NL / 3>(mo.block_A(), mo.block_W(), mo.excitation(), A, W);
    const std::complex<double> s(0.0, 2 * PI * f);
    const Eigen::Matrix<std::complex<double>,NL,NL> R =
        (s * Eigen::Matrix<std::complex<double>,NL,NL>::Identity() - A.template cast<std::complex<double>>()).inverse();
    Eigen::Matrix<double,1,NL> ci = Eigen::Matrix<double,1,NL>::Zero(), cj = ci;
    ci(Mo::pos_local(i)) = 1.0; cj(Mo::pos_local(j)) = 1.0;
    const Eigen::Matrix<std::complex<double>,1,NL> gi = ci.template cast<std::complex<double>>() * R;
    const Eigen::Matrix<std::complex<double>,1,NL> gj = cj.template cast<std::complex<double>>() * R;
    return (gi * W.template cast<std::complex<double>>() * gj.adjoint())(0,0);
}

void test6_phase() {
    std::printf("[6] cross-spectral phase\n");
    kc::Hermitian3<double> S = kc::Hermitian3<double>::diagonal(V3<double>(1.0, 0.8, 1.0));
    S.set(0, 2, std::polar(0.7, -62.0 * PI / 180.0));
    S.set(1, 2, std::polar(0.5, -55.0 * PI / 180.0));
    S.set(0, 1, std::complex<double>(0.4, 0.0));
    Cfg c; c.source = Cfg::Source::External; c.shrink = false;
    const double kappa = 0.9, fop = kappa * 0.37083;
    kc::CarmaCorrelatedMotion<double> mc;
    mc.coupling().set_config(c); mc.coupling().set_external(S);
    mc.set_vertical(kappa, 1.0); mc.set_horizontal(kappa, 1.0, std::sqrt(0.8), 0.0); mc.refresh(0.005);
    double imr = 0;
    for (double f : {0.5 * fop, fop, 2.0 * fop}) {
        const std::complex<double> v = model_cross<kc::CarmaCorrelatedMotion<double>, 12>(mc, f, 0, 2);
        imr = std::max(imr, std::abs(v.imag()) / std::abs(v));
    }
    check(imr < 1e-10, "COR cross-spectrum is real (co-spectrum only)", imr, 1e-10);
}

// ---------------------------------------------------------------------------
void test7_frames() {
    std::printf("[7] frame consistency\n");
    std::mt19937_64 rng(17);
    std::uniform_real_distribution<double> ua(-PI, PI);
    double worst = 0;
    for (int trial = 0; trial < 200; ++trial) {
        kc::Hermitian3<double> S = random_psd<double>(rng, 0.8);
        const double a = ua(rng);
        Eigen::Matrix3d R = Eigen::AngleAxisd(a, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        kc::Hermitian3<double> Sr; Sr.re = R * S.re * R.transpose(); Sr.im = R * S.im * R.transpose();
        kc::CarmaCoupling<double> c1, c2;
        Cfg c; c.source = Cfg::Source::External; c.shrink = trial % 2 == 0;
        c.form = static_cast<Cfg::Form>(trial % 3);
        c1.set_config(c); c2.set_config(c);
        c1.set_external(S, 25.0); c2.set_external(Sr, 25.0);
        const V3<double> sig(1.1, 1.1, 1.1);   // the wrapper's isotropic prior
        const std::complex<double> one[3] = {1.0, 1.0, 1.0};
        const kc::Hermitian3<double> E1 = c1.excitation(sig, one), E2 = c2.excitation(sig, one);
        worst = std::max(worst, (R * E1.re * R.transpose() - E2.re).cwiseAbs().maxCoeff() +
                                (R * E1.im * R.transpose() - E2.im).cwiseAbs().maxCoeff());
    }
    check(worst < 1e-9, "excitation mapping is equivariant under rotations about the vertical (every form, shrinkage)", worst, 1e-9);
    kc::Hermitian3<double> S = kc::Hermitian3<double>::diagonal(V3<double>(1.0, 0.2, 1.0));
    S.set(0, 2, {0.3, -0.8});
    for (int keel = 0; keel < 2; ++keel) {
        kc::CarmaCoupling<double> cp;
        Cfg c; c.source = Cfg::Source::External; c.shrink = false; c.frame = keel ? Cfg::Frame::Keel : Cfg::Frame::World;
        cp.set_config(c);
        cp.set_external(S);
        cp.set_heading(0.3);
        cp.freeze_now();
        cp.set_heading(0.3 + PI / 2);
        kc::Hermitian3<double> B; double n;
        cp.band_cross(B, n);
        Eigen::Matrix3d R = Eigen::AngleAxisd(keel ? PI / 2 : 0.0, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        const double e = (R * S.re * R.transpose() - B.re).cwiseAbs().maxCoeff() + (R * S.im * R.transpose() - B.im).cwiseAbs().maxCoeff();
        check(e < 1e-12, keel ? "keel-frame frozen coupling turns with the hull" : "world-frame frozen coupling stays fixed in NED", e, 1e-12);
    }
}

// ---------------------------------------------------------------------------
void test8_reset() {
    std::printf("[8] QMEKF reset and cross-covariance\n");
    auto b = make_filter<Base<float>>();
    auto c = make_filter<COR<float>>();
    const V3<float> sig(0.4f, 0.4f, 0.6f);
    configure(b, 1.0f, sig); configure(c, 1.0f, sig);
    // Start with a 12 deg tilt error and a large yaw uncertainty.
    const Eigen::Quaternionf q0(Eigen::AngleAxisf(0.21f, V3<float>(1, 1, 0).normalized()));
    b.initialize_from_attitude(q0, 0.3f, 1.0f); c.initialize_from_attitude(q0, 0.3f, 1.0f);
    SynthImu imu;
    double worst_c = 0, minpsd = 1, as = 0;
    for (int k = 0; k < 4000; ++k) {
        V3<float> acc, gyr;
        imu.sample(0.005, acc, gyr);
        b.time_update(gyr, 0.005f); b.measurement_update_acc_only(acc);
        c.time_update(gyr, 0.005f); c.measurement_update_acc_only(acc);
        if (k % 200 == 199) {
            for (int i = 0; i < 6; ++i)
                for (int ax = 0; ax < 3; ++ax) {
                    const int jb = b.pos_index(ax), jc = c.pos_index(ax);
                    const double ref = std::sqrt(b.covariance_full()(i,i) * b.covariance_full()(jb,jb));
                    worst_c = std::max(worst_c, std::abs(double(c.covariance_full()(i,jc)) - b.covariance_full()(i,jb)) / ref);
                }
            minpsd = std::min(minpsd, min_corr_eig(c.covariance_full()));
            as = std::max(as, double(asym(c.covariance_full())));
        }
    }
    check(worst_c < 1e-4, "COR off: attitude/bias-position cross-covariances equal the baseline through resets", worst_c, 1e-4);
    check(minpsd > -1e-4, "covariance PSD through large-tilt resets (correlation-scaled)", minpsd, -1e-4);
    check(as == 0.0, "covariance exactly symmetric", as, 0.0);
}

// ---------------------------------------------------------------------------
template<class F>
void run_changing(F& f, double T_s, double& minpsd, double& pmax, bool& finite, double amp_scale = 1.0) {
    std::mt19937_64 rng(23);
    std::normal_distribution<double> n01(0.0, 1.0);
    SynthImu imu;
    imu.amp = amp_scale;
    double kappa = 0.9, sig = 0.5;
    const long N = static_cast<long>(T_s / 0.005);
    for (long k = 0; k < N; ++k) {
        V3<float> acc, gyr;
        imu.sample(0.005, acc, gyr);
        if (k % 20 == 0) {   // commit every 0.1 s
            kappa = std::min(2.5, std::max(0.4, kappa * std::exp(0.02 * n01(rng))));
            sig = std::min(2.0, std::max(1e-5, sig * std::exp(0.03 * n01(rng))));
            f.set_vertical_carma(static_cast<float>(kappa), static_cast<float>(sig));
            f.set_horizontal_carma(static_cast<float>(kappa), static_cast<float>(sig), static_cast<float>(sig), 0.0f);
        }
        f.time_update(gyr, 0.005f);
        f.measurement_update_acc_only(acc);
        if (k % 1000 == 999) {
            minpsd = std::min(minpsd, min_corr_eig(f.covariance_full()));
            finite = finite && f.covariance_full().allFinite() && f.get_position().allFinite();
            pmax = std::max(pmax, double(f.get_position().cwiseAbs().maxCoeff()));
        }
    }
}

template<class F>
void set_estimated(F& f, float horizon = 120.0f, float warmup = 60.0f) {
    Cfg c; c.source = Cfg::Source::Estimated; c.warmup_s = warmup;
    f.motion().coupling().set_config(c);
    typename kc::CarmaCrossSpectrumEstimator<float>::Config ec;
    ec.horizon_s = horizon;
    f.motion().coupling().estimator().configure(ec);
}

void test9_changing() {
    std::printf("[9] changing sea-state parameters\n");
    auto f = make_filter<COR<float>>();
    set_estimated(f);
    configure(f, 0.9f, V3<float>(0.5f, 0.5f, 0.5f));
    init_level(f);
    double mp = 1, pm = 0; bool fin = true;
    run_changing(f, 600.0, mp, pm, fin);
    check(fin, "COR: finite under 10 Hz kappa / amplitude / coupling changes");
    check(mp > -1e-4, "COR: covariance PSD under parameter changes", mp, -1e-4);
    check(pm < 5.0, "COR: displacement bounded under parameter changes (m)", pm, 5.0);
}

// ---------------------------------------------------------------------------
void test11_12_noise() {
    std::printf("[11] near-zero excitation, [12] no false coupling from sensor noise\n");
    // Pure sensor noise (still IMU) into the estimator: debiased, shrunk coupling ~ 0.
    kc::CarmaCoupling<double> cp;
    Cfg c; c.source = Cfg::Source::Estimated; c.warmup_s = 60.0f;
    cp.set_config(c);
    std::mt19937_64 rng(31);
    std::normal_distribution<double> n01(0.0, 1.0);
    double coh_max = 0;
    for (long k = 0; k < static_cast<long>(1200.0 / 0.005); ++k) {
        cp.observe(ACC_SIGMA * V3<double>(n01(rng), n01(rng), n01(rng)), 0.005, 0.3);
        if (k > static_cast<long>(200.0 / 0.005) && k % 2000 == 0) {
            const std::complex<double> one[3] = {1.0, 1.0, 1.0};
            const auto E = cp.excitation(V3<double>(0.05, 0.05, 0.05), one);
            for (int i = 0; i < 3; ++i)
                for (int j = i + 1; j < 3; ++j) coh_max = std::max(coh_max, std::norm(E(i,j)) / (E.re(i,i) * E.re(j,j)));
        }
    }
    check(coh_max < 0.05, "sensor noise alone: applied cross-axis coherence after shrinkage", coh_max, 0.05);
    // Still IMU through the full filters: no extra displacement vs the baseline.
    auto b = make_filter<Base<float>>();
    auto s = make_filter<COR<float>>();
    set_estimated(s);
    const V3<float> sig = V3<float>::Constant(0.05f);
    configure(b, 0.135f, sig); configure(s, 0.135f, sig);   // tuner floor: 0.05 Hz, 0.05 m/s^2
    init_level(b); init_level(s);
    SynthImu imu; imu.amp = 0.0; imu.psi_rate = 0.0;
    double pb = 0, ps = 0; long n = 0;
    for (long k = 0; k < static_cast<long>(900.0 / 0.005); ++k) {
        V3<float> acc, gyr;
        imu.sample(0.005, acc, gyr);
        b.time_update(gyr, 0.005f); s.time_update(gyr, 0.005f);
        b.measurement_update_acc_only(acc); s.measurement_update_acc_only(acc);
        if (static_cast<double>(k) * 0.005 > 300) { pb += b.get_position().squaredNorm(); ps += s.get_position().squaredNorm(); ++n; }
    }
    check(std::sqrt(ps / static_cast<double>(n)) < 1.03 * std::sqrt(pb / static_cast<double>(n)), "still IMU: COR (estimated coupling) phantom <= baseline + 3 %",
          std::sqrt(ps / pb), 1.03);
    // Near-zero excitation (numerical floor 1e-5 m/s^2) with coupling on.
    {
        double mp = 1, pm = 0; bool fin = true;
        auto go = [&](auto& f) {
            set_estimated(f);
            configure(f, 0.5f, V3<float>::Constant(1e-5f));
            init_level(f);
            SynthImu im2; im2.amp = 0.0;
            for (long k = 0; k < static_cast<long>(300.0 / 0.005); ++k) {
                V3<float> acc, gyr;
                im2.sample(0.005, acc, gyr);
                f.time_update(gyr, 0.005f);
                f.measurement_update_acc_only(acc);
                if (k % 1000 == 999) {
                    fin = fin && f.covariance_full().allFinite();
                    mp = std::min(mp, min_corr_eig(f.covariance_full()));
                    pm = std::max(pm, double(f.get_position().cwiseAbs().maxCoeff()));
                }
            }
        };
        auto f = make_filter<COR<float>>();
        go(f);
        check(fin && mp > -1e-3, "COR near-zero excitation: finite, PSD", mp, -1e-3);
        check(pm < 1e-3, "COR near-zero excitation: no displacement (m)", pm, 1e-3);
    }
}

// ---------------------------------------------------------------------------
void test13_extended() {
    std::printf("[13] extended operation\n");
    {
        double mp = 1, pm = 0; bool fin = true;
        auto go = [&](auto& f, double T_s) {
            set_estimated(f, 240.0f, 120.0f);
            configure(f, 0.9f, V3<float>::Constant(0.5f));
            init_level(f);
            SynthImu imu;
            for (long k = 0; k < static_cast<long>(T_s / 0.005); ++k) {
                V3<float> acc, gyr;
                imu.sample(0.005, acc, gyr);
                f.time_update(gyr, 0.005f);
                f.measurement_update_acc_only(acc);
                if (k % 20000 == 19999) {
                    fin = fin && f.covariance_full().allFinite() && asym(f.covariance_full()) == 0.0f;
                    mp = std::min(mp, min_corr_eig(f.covariance_full()));
                    pm = std::max(pm, double(f.get_position().cwiseAbs().maxCoeff()));
                }
            }
        };
        auto f = make_filter<COR<float>>();
        go(f, 3600.0);
        check(fin, "COR 3600 s: finite, exactly symmetric");
        check(mp > -1e-4, "COR 3600 s: PSD (correlation-scaled)", mp, -1e-4);
        check(pm < 2.0, "COR 3600 s: displacement bounded (m)", pm, 2.0);
    }
}

void test14_layout() {
    std::printf("[14] baseline layout and state counts\n");
    check(Base<float>::NX == 21 && !Base<float>::kMotion, "default instantiation is the 21-state baseline");
    check(COR<float>::NX == 21, "Experiment A keeps 21 states");
}


// ---------------------------------------------------------------------------
void test16_modal() {
    std::printf("[16] exact complex-OU modal form of CARMA(4,1)\n");
    using M4 = Eigen::Matrix<double,4,4>;
    using V4 = Eigen::Matrix<double,4,1>;
    double worst_tf = 0, worst_F = 0, worst_Q = 0, worst_P = 0, worst_ac = 0, kin = 0, rot_lost = 1e9, cond_max = 0;
    for (double kappa : {0.2, 0.5, 0.85, 1.4, 2.5}) {
        kc::CarmaModalMotion<double> m;
        m.set_vertical(kappa, 1.0); m.set_horizontal(kappa, 1.0, 1.0, 0.0); m.refresh(0.005);
        cond_max = std::max(cond_max, m.transform_condition());
        const M4 A = m.block_A();
        const V4 Gm = m.block_G() * std::sqrt(m.noise_gain2());
        const Eigen::Matrix<double,1,4> cp = m.readout_p();
        // Transfer function of the modal realization vs |B/P|^2 (unit gain).
        for (double f : {0.02, 0.1, 0.25, 0.37, 0.6, 1.0, 3.0}) {
            const double fk = f * kappa;
            const std::complex<double> s(0.0, 2 * PI * fk);
            const Eigen::Matrix<std::complex<double>,4,4> R = (s * Eigen::Matrix<std::complex<double>,4,4>::Identity() - A.cast<std::complex<double>>()).inverse();
            const std::complex<double> H = (cp.cast<std::complex<double>>() * R * Gm.cast<std::complex<double>>())(0,0);
            const double ref = m.noise_gain2() * kc::carma_unit_psd<4>(kc::CARMA41_JONSWAP_SHAPE, kappa, fk);
            worst_tf = std::max(worst_tf, std::abs(std::norm(H) / ref - 1.0));
        }
        // Similarity to the baseline realization: x_base = T x_modal.
        M4 Ab; kc::carma_dynamics<4>(kc::CARMA41_JONSWAP_SHAPE, kappa, Ab);
        M4 Ob, Om; Eigen::Matrix<double,1,4> rb = Eigen::Matrix<double,1,4>::Zero(), rm = cp; rb(0) = 1.0;
        for (int i = 0; i < 4; ++i) { Ob.row(i) = rb; Om.row(i) = rm; rb = rb * Ab; rm = rm * A; }
        const M4 T = Ob.inverse() * Om;
        M4 Wb = M4::Zero(); Wb(2,2) = m.noise_gain2();
        M4 Fb, Qb; kc::discretize_lti<double,4>(Ab, Wb, 0.005, Fb, Qb);
        const M4 Fm = T.inverse() * Fb * T, Qm = T.inverse() * Qb * T.inverse().transpose();
        worst_F = std::max(worst_F, (Fm - m.block_F()).cwiseAbs().maxCoeff());
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                worst_Q = std::max(worst_Q, std::abs(Qm(i,j) - m.block_Q()(i,j)) / std::sqrt(Qm(i,i) * Qm(j,j)));
        // Closed-form stationary covariance vs the baseline's Lyapunov solution.
        Eigen::Matrix<double,12,12> P12 = Eigen::Matrix<double,12,12>::Zero();
        m.stationary(P12, 0);
        const M4 Pm = P12.topLeftCorner<4,4>();
        M4 Pb; kc::stationary_covariance_doubling<double,4>(Ab, Wb, Pb);
        const M4 Pbt = T * Pm * T.transpose();
        worst_P = std::max(worst_P, (Pbt - Pb).cwiseAbs().maxCoeff() / Pb.cwiseAbs().maxCoeff());
        // Autocovariance of p at lags 0.5, 2, 5 s.
        for (double lag : {0.5, 2.0, 5.0}) {
            M4 Fl, Ql, Fbl, Qbl;
            kc::discretize_lti<double,4>(A, M4::Zero(), lag, Fl, Ql);
            kc::discretize_lti<double,4>(Ab, M4::Zero(), lag, Fbl, Qbl);
            const double am = (cp * Fl * Pm * cp.transpose())(0,0), ab = (Fbl * Pb)(0,0);
            worst_ac = std::max(worst_ac, std::abs(am - ab) / Pb(0,0));
        }
        // Relative degree 3: c G = c A G = 0 (residue cancellations).
        kin = std::max(kin, std::max(std::abs((cp * Gm)(0,0)), std::abs((cp * A * Gm)(0,0))) / std::abs((cp * A * A * Gm)(0,0)));
        // A readout rotated by 60 deg loses it: c_rot G != 0.
        Eigen::Matrix<double,1,4> cr;
        const double c60 = 0.5, s60 = std::sqrt(3.0) / 2.0;
        cr << 2 * c60, -2 * s60, 2 * c60, -2 * s60;
        rot_lost = std::min(rot_lost, std::abs((cr * Gm)(0,0)) / std::abs((cp * A * A * Gm)(0,0)) * kappa * kappa);
    }
    check(worst_tf < 1e-9, "modal |H|^2 equals C41J g^2 |B/P|^2", worst_tf, 1e-9);
    check(worst_F < 1e-11, "closed-form modal F = T^-1 F_baseline T", worst_F, 1e-11);
    check(worst_Q < 1e-9, "closed-form modal Q (common noise) = T^-1 Q_baseline T^-T (scaled)", worst_Q, 1e-9);
    check(worst_P < 1e-10, "closed-form modal stationary covariance = baseline Lyapunov (transformed)", worst_P, 1e-10);
    check(worst_ac < 1e-10, "autocovariance of p at 0.5-5 s equals the baseline's", worst_ac, 1e-10);
    check(kin < 1e-12, "residue cancellations: c G = c A G = 0 (p, v, a noise-free readouts)", kin, 1e-12);
    check(rot_lost > 1e-3, "a 60 deg phase-rotated modal readout loses them (white velocity)", rot_lost, 1e-3);
    check(cond_max < 1e4, "modal transformation condition number over kappa 0.2-2.5", cond_max, 1e4);
    // Independently driven oscillators (same input magnitudes, no common
    // noise) lose the cancellations: Q_pp grows like h, not h^5.
    {
        kc::CarmaModalMotion<double> m;
        m.set_vertical(0.85, 1.0); m.set_horizontal(0.85, 1.0, 1.0, 0.0); m.refresh(0.005);
        const V4 Gm = m.block_G();
        M4 Wc = Gm * Gm.transpose() * m.noise_gain2(), Wi = M4::Zero();
        Wi.topLeftCorner<2,2>() = Wc.topLeftCorner<2,2>();
        Wi.bottomRightCorner<2,2>() = Wc.bottomRightCorner<2,2>();
        const Eigen::Matrix<double,1,4> cp = m.readout_p(), ca = m.readout_a();
        double qc[2], qi[2];
        for (int k = 0; k < 2; ++k) {
            const double h = k ? 0.01 : 0.005;
            M4 F, Qc, Qi;
            kc::discretize_lti<double,4>(m.block_A(), Wc, h, F, Qc);
            kc::discretize_lti<double,4>(m.block_A(), Wi, h, F, Qi);
            qc[k] = (cp * Qc * cp.transpose())(0,0);
            qi[k] = (cp * Qi * cp.transpose())(0,0);
        }
        const double slope_c = std::log2(qc[1] / qc[0]), slope_i = std::log2(qi[1] / qi[0]);
        check(std::abs(slope_c - 5.0) < 0.05, "common-noise modes: Q_pp ~ h^5 (relative degree 3)", slope_c, 5.0);
        check(std::abs(slope_i - 1.0) < 0.05, "independently driven oscillators: Q_pp ~ h (white velocity)", slope_i, 1.0);
        M4 Pi;
        kc::stationary_covariance_doubling<double,4>(m.block_A(), Wi, Pi);
        const double va_i = (ca * Pi * ca.transpose())(0,0);
        M4 Pc;
        kc::stationary_covariance_doubling<double,4>(m.block_A(), Wc, Pc);
        const double vp_c = (cp * Pc * cp.transpose())(0,0), vp_i = (cp * Pi * cp.transpose())(0,0);
        check(std::abs(vp_i / vp_c - 1.0) > 0.05, "independently driven oscillators change the displacement variance", vp_i / vp_c, 1.05);
        (void)va_i;
    }
    // Filter equivalence in double, fixed and changing parameters.
    for (int changing = 0; changing < 2; ++changing) {
        auto b = make_filter<Base<double>>();
        auto mo = make_filter<MOD<double>>();
        configure(b, 0.85, V3<double>(0.35, 0.30, 0.55)); configure(mo, 0.85, V3<double>(0.35, 0.30, 0.55));
        init_level(b); init_level(mo);
        SynthImu imu;
        std::mt19937_64 rng(5);
        std::normal_distribution<double> n01(0.0, 1.0);
        double kappa = 0.85, sig = 0.5, dp = 0, pp = 0, dv = 0;
        for (int k = 0; k < 60000; ++k) {
            V3<double> acc, gyr;
            imu.sample(0.005, acc, gyr);
            if (changing && k % 20 == 0) {
                kappa = std::min(2.0, std::max(0.4, kappa * std::exp(0.02 * n01(rng))));
                sig = std::min(1.5, std::max(0.05, sig * std::exp(0.03 * n01(rng))));
                for (auto* f : {static_cast<void*>(&b), static_cast<void*>(&mo)}) (void)f;
                b.set_vertical_carma(kappa, sig); b.set_horizontal_carma(kappa, sig, sig, 0.0);
                mo.set_vertical_carma(kappa, sig); mo.set_horizontal_carma(kappa, sig, sig, 0.0);
            }
            b.time_update(gyr, 0.005); mo.time_update(gyr, 0.005);
            b.measurement_update_acc_only(acc); mo.measurement_update_acc_only(acc);
            dp = std::max(dp, (mo.get_position() - b.get_position()).cwiseAbs().maxCoeff());
            pp += b.get_position().squaredNorm();
            dv = std::max(dv, (mo.pv_covariance6() - b.pv_covariance6()).cwiseAbs().maxCoeff() / b.pv_covariance6().cwiseAbs().maxCoeff());
        }
        const double rel = dp / std::sqrt(pp / 60000.0);
        check(rel < 1e-8, changing ? "modal filter = baseline (double) with 10 Hz kappa/amplitude changes (re-coordinatized)"
                                   : "modal filter = baseline (double), fixed parameters", rel, 1e-8);
        check(dv < 1e-8, changing ? "modal [p, v] covariance = baseline with parameter changes"
                                  : "modal [p, v] covariance = baseline, fixed parameters", dv, 1e-8);
    }
    // Simulation statistics of the modal model (exact F, Q, double): sample
    // variances of p, v, a against the stationary values.
    {
        kc::CarmaModalMotion<double> m;
        m.set_vertical(0.85, 1.0); m.set_horizontal(0.85, 1.0, 1.0, 0.0); m.refresh(0.02);
        const M4 F = m.block_F();
        Eigen::SelfAdjointEigenSolver<M4> es(m.block_Q());
        const M4 L = es.eigenvectors() * es.eigenvalues().cwiseMax(0.0).cwiseSqrt().asDiagonal();
        Eigen::Matrix<double,12,12> P12 = Eigen::Matrix<double,12,12>::Zero();
        m.stationary(P12, 0);
        const M4 P = P12.topLeftCorner<4,4>();
        std::mt19937_64 rng(77);
        std::normal_distribution<double> n01(0.0, 1.0);
        V4 x = V4::Zero();
        double s2[3] = {0, 0, 0};
        const long N = 2000000;
        for (long k = 0; k < N; ++k) {
            V4 w; for (int i = 0; i < 4; ++i) w(i) = n01(rng);
            x = F * x + L * w;
            const double pv[3] = {(m.readout_p() * x)(0,0), (m.readout_v() * x)(0,0), (m.readout_a() * x)(0,0)};
            for (int i = 0; i < 3; ++i) s2[i] += pv[i] * pv[i];
        }
        const Eigen::Matrix<double,1,4> rows[3] = {m.readout_p(), m.readout_v(), m.readout_a()};
        double worst = 0;
        for (int i = 0; i < 3; ++i) {
            const double ref = (rows[i] * P * rows[i].transpose())(0,0);
            worst = std::max(worst, std::abs(s2[i] / static_cast<double>(N) / ref - 1.0));
        }
        check(worst < 0.05, "simulated modal p, v, a variances match the stationary covariance (11 h at 50 Hz)", worst, 0.05);
    }
    // Coalescing poles: the modal form refuses (baseline realization is the fallback).
    {
        kc::CarmaShape sh = kc::CARMA41_JONSWAP_SHAPE;
        sh.f2_hz = sh.f1_hz * (1.0 + 1e-5); sh.zeta2 = sh.zeta1;
        kc::CarmaModalMotion<double> m(sh);
        m.set_vertical(1.0, 1.0); m.refresh(0.005);
        check(!m.valid(), "modal form refuses (nearly) coalescing poles");
    }
}

}  // namespace

// ---------------------------------------------------------------------------
void test17_modal_numerator() {
    std::printf("[17] directional modal coupling (per-axis numerator zeros)\n");
    using M4 = Eigen::Matrix<double,4,4>;
    using V4 = Eigen::Matrix<double,4,1>;
    using Cc = std::complex<double>;
    const double kappa = 0.85;
    const double phx = -60.0 * PI / 180.0, phy = 30.0 * PI / 180.0, rx = 0.8, ry = 0.5;
    kc::Hermitian3<double> S = kc::Hermitian3<double>::diagonal(V3<double>(1.0, 1.0, 1.0));
    const Cc sxz = std::polar(rx, phx), syz = std::polar(ry, phy);
    S.re(0,2) = S.re(2,0) = sxz.real(); S.im(0,2) = sxz.imag(); S.im(2,0) = -sxz.imag();
    S.re(1,2) = S.re(2,1) = syz.real(); S.im(1,2) = syz.imag(); S.im(2,1) = -syz.imag();
    auto build = [&](kc::CarmaModalMotion<double>& m, Cfg::Form form, const kc::Hermitian3<double>& Sx) {
        Cfg c; c.source = Cfg::Source::External; c.shrink = false; c.form = form;
        m.coupling().set_config(c);
        m.coupling().set_external(Sx);
        m.set_vertical(kappa, 1.0); m.set_horizontal(kappa, 1.0, 1.0, 0.0);
        m.refresh(0.005);
    };
    kc::CarmaModalMotion<double> m;
    build(m, Cfg::Form::Numerator, S);
    const M4 A = m.block_A();
    const V4 Gm = m.block_G();
    const double w = 2 * PI * m.operating_frequency_hz();
    const Cc s(0.0, w);
    const Eigen::Matrix<Cc,4,4> R = (s * Eigen::Matrix<Cc,4,4>::Identity() - A.cast<Cc>()).inverse();
    Cc H[3];
    double kin = 0.0;
    for (int ax = 0; ax < 3; ++ax) {
        const Eigen::Matrix<double,1,4>& c = m.readout_p(ax);
        H[ax] = (c.cast<Cc>() * R * Gm.cast<Cc>())(0,0);
        const double n2 = std::abs((c * A * A * Gm)(0,0));
        kin = std::max(kin, std::max(std::abs((c * Gm)(0,0)), std::abs((c * A * Gm)(0,0))) / n2);
        kin = std::max(kin, ((m.readout_a(ax) - c * A * A).cwiseAbs().maxCoeff()));
    }
    check(m.numerator_zero(0) != 0.0 && m.numerator_zero(1) != 0.0 && m.numerator_zero(2) == 0.0,
          "numerator zeros on x and y only");
    check(kin < 1e-12, "numerator readouts keep relative degree 3 (c G = c A G = 0, a = c A^2 x)", kin, 1e-12);
    const auto& E = m.excitation();
    double ph_err = 0.0, gain_err = 0.0;
    const double ph_t[2] = {phx, phy}, rho[2] = {rx, ry};
    for (int i = 0; i < 2; ++i) {
        const Cc sm = E.re(i,2) * H[i] * std::conj(H[2]);
        const double target = std::atan(rho[i] * std::tan(ph_t[i]));
        ph_err = std::max(ph_err, std::abs(std::arg(sm) - target));
        gain_err = std::max(gain_err, std::abs(E.re(i,i) * std::norm(H[i]) / (E.re(2,2) * std::norm(H[2])) - 1.0));
    }
    check(ph_err < 1e-9, "model x-z / y-z phase at f_op = atan(rho tan(phase S))", ph_err, 1e-9);
    check(gain_err < 1e-9, "band gain at f_op preserved (S_xx = S_zz)", gain_err, 1e-9);
    check(E.im.cwiseAbs().maxCoeff() == 0.0, "numerator excitation is real (modal noise stays real)");
    // Zero phase: no numerator, same excitation as the real form.
    kc::Hermitian3<double> S0 = S;
    S0.re(0,2) = S0.re(2,0) = rx; S0.re(1,2) = S0.re(2,1) = -ry;
    S0.im.setZero();
    kc::CarmaModalMotion<double> m0, m1;
    build(m0, Cfg::Form::Numerator, S0);
    build(m1, Cfg::Form::Hermitian, S0);
    check(m0.numerator_zero(0) == 0.0 && m0.numerator_zero(1) == 0.0 &&
          (m0.excitation().re - m1.excitation().re).cwiseAbs().maxCoeff() < 1e-15 &&
          (m0.readout_p(0) - m1.readout_p(0)).cwiseAbs().maxCoeff() == 0.0,
          "zero phase: numerator form = real modal coupling");
    // Random Hermitian inputs: real excitation stays PSD, readouts finite.
    std::mt19937_64 rng(17);
    double worst = 1.0;
    bool finite = true;
    for (int trial = 0; trial < 300; ++trial) {
        kc::Hermitian3<double> Sr = random_psd<double>(rng, trial % 3 == 0 ? 1.0 : 0.3);
        kc::CarmaModalMotion<double> mr;
        Cfg c; c.source = Cfg::Source::External; c.shrink = false; c.form = Cfg::Form::Numerator;
        mr.coupling().set_config(c);
        mr.coupling().set_external(Sr);
        mr.set_vertical(kappa, std::sqrt(Sr.re(2,2))); mr.set_horizontal(kappa, std::sqrt(Sr.re(0,0)), std::sqrt(Sr.re(1,1)), 0.0);
        mr.refresh(0.005);
        worst = std::min(worst, min_eig_real_rep(mr.excitation()));
        for (int ax = 0; ax < 3; ++ax) finite = finite && mr.readout_p(ax).allFinite() && mr.readout_a(ax).allFinite();
    }
    check(worst > -1e-9, "numerator excitation PSD (300 random Hermitian inputs)", worst, -1e-9);
    check(finite, "numerator readouts finite (300 random inputs)");
}

int main() {
    test14_layout();
    test1_reduction();
    test2_psd();
    test3_c2d();
    test6_phase();
    test7_frames();
    test8_reset();
    test9_changing();
    test11_12_noise();
    test13_extended();
    test16_modal();
    test17_modal_numerator();
    std::printf("carma_coupled-test: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
