/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    carma_shared_modes_realization-test -- unit and numerical tests of the
    motion block of the EXPERIMENTAL shared-modes QMEKF
    (src/kalman_carma_shared_modes/CarmaSharedModesMotion.h,
    docs/kalman-carma-shared-modes.md).

    Groups:
      [1] structure: 27 error states, six per axis, a valid realization,
          lossless colouring, kappa scaling
      [2] realization: against a test-local eight-state realified complexified
          mode zeta = B/(P+ D) xi, the per-axis transfer functions from the two
          noise inputs to p, v, a are equal, and so are the stationary
          three-axis law of (p, v, a) under a Hermitian excitation and its
          lag-h cross-covariance (co- and quadrature spectra)
      [3] with the coupling off the filter reproduces the CARMA3 baseline A41
          through kappa and amplitude adaptation (double cores)
      [4] commits: p, v, a are continuous through the core's commit path; the
          applied transforms keep the p, v, a rows; the covariance stays
          positive semidefinite
      [5] float wrapper on a synthetic directional sea: finite, symmetric,
          PSD covariance; the coupling engages; accuracy against the
          independent-axis CARMA3 wrapper on the same stream
*/

#define EIGEN_NON_ARDUINO

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"

extern const float g_std = 9.80665f;   // the application constant the front end reads

namespace csm = ocean_imu::carma_shared_modes;
namespace kc = ocean_imu::kalman::carma;
using V3d = Eigen::Vector3d;
using V3f = Eigen::Vector3f;
using Cd = std::complex<double>;
using Motion = csm::SharedModesMotion<double>;

static int g_checks = 0, g_failures = 0;
static void check(bool ok, const char* what, double value = NAN, double limit = NAN) {
    ++g_checks;
    if (!ok) ++g_failures;
    std::printf("  %s: %s (value=%.4g limit=%.4g)\n", ok ? "ok" : "FAIL", what, value, limit);
}

constexpr double G = 9.80665;
constexpr double PI = 3.14159265358979323846;

// Directional sea with exact kinematics (as carma_shared_modes-test): a
// dominant system toward dir_deg with linear-wave quadrature between the
// horizontal displacement and heave, plus an incoherent cross component.
struct Sea {
    struct Comp { double A, w, ph, dir; bool quad; };
    std::vector<Comp> c;
    explicit Sea(double dir_deg, unsigned seed = 1, double hs = 1.0, double fp = 0.25) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> u(0.0, 2 * PI);
        const int n = 24;
        for (int k = 0; k < n; ++k) {
            const double f = fp * (0.7 + 0.9 * k / (n - 1));
            const double a = 0.25 * hs * std::exp(-2.0 * std::pow((f - fp) / (0.25 * fp), 2.0)) / std::sqrt(n / 4.0);
            c.push_back({a, 2 * PI * f, u(rng), dir_deg * PI / 180.0, true});
        }
        for (int k = 0; k < 6; ++k)
            c.push_back({0.03 * hs, 2 * PI * fp * (0.8 + 0.1 * k), u(rng), (dir_deg + 90.0) * PI / 180.0, false});
    }
    void eval(double t, V3d& p, V3d& a) const {
        p.setZero(); a.setZero();
        for (const auto& k : c) {
            const double th = k.w * t + k.ph;
            const double eta = k.A * std::cos(th), eta2 = -k.w * k.w * eta;
            const double xi = k.quad ? k.A * std::sin(th) : 0.4 * k.A * std::cos(th + 1.0);
            const double xi2 = -k.w * k.w * xi;
            p += V3d(std::cos(k.dir) * xi, std::sin(k.dir) * xi, -eta);
            a += V3d(std::cos(k.dir) * xi2, std::sin(k.dir) * xi2, -eta2);
        }
    }
};

// Level, north-facing body (body = NED).
struct Imu {
    std::mt19937 rng{7};
    std::normal_distribution<double> n01{0.0, 1.0};
    double sa = 0.0148, sg = 0.0016, sm = 0.8;
    void sample(const V3d& a_ned, V3d& acc, V3d& gyr) {
        acc = a_ned - V3d(0, 0, G) + sa * V3d(n01(rng), n01(rng), n01(rng));
        gyr = sg * V3d(n01(rng), n01(rng), n01(rng));
    }
    V3d mag() { return V3d(20.0, 0.0, 45.0) + sm * V3d(n01(rng), n01(rng), n01(rng)); }
};

// ---------------------------------------------------------------------------
// [1] structure
void test1_structure() {
    std::printf("[1] structure\n");
    check(csm::KalmanCarmaSharedModes<float>::kStates == 27 && csm::KalmanCarmaSharedModes<double>::kStates == 27 &&
          Motion::NB == 6 && Motion::NL == 18,
          "27 error states (3 + 3 + 18 + 3), six motion states per axis");
    Motion m;
    check(m.valid(), "realization valid");
    const Eigen::Matrix<double,6,6> A = m.block_A();
    Eigen::EigenSolver<Eigen::Matrix<double,6,6>> es(A);
    check(es.eigenvalues().real().maxCoeff() < 0.0, "block Hurwitz", es.eigenvalues().real().maxCoeff(), 0.0);
    // kappa scaling: A(kappa) has the eigenvalues kappa * eig A(1).
    m.set_vertical(1.6, 1.0);
    m.set_horizontal(1.6, 1.0, 1.0, 0.0);
    m.refresh(0.005);
    Eigen::EigenSolver<Eigen::Matrix<double,6,6>> ek(m.block_A());
    std::vector<double> e1, e2;
    for (int i = 0; i < 6; ++i) { e1.push_back(std::abs(es.eigenvalues()(i)) * 1.6); e2.push_back(std::abs(ek.eigenvalues()(i))); }
    std::sort(e1.begin(), e1.end()); std::sort(e2.begin(), e2.end());
    double d = 0.0;
    for (int i = 0; i < 6; ++i) d = std::max(d, std::abs(e1[i] - e2[i]) / e1[i]);
    check(d < 1e-10, "poles scale with kappa", d, 1e-10);
    // Lossless colouring: the stationary covariance of the two colouring
    // states is g2 sigma^2 I (input-normal coordinates).
    Motion m1;
    Eigen::Matrix<double,18,18> P = Eigen::Matrix<double,18,18>::Zero();
    m1.stationary(P, 0);
    const Eigen::Matrix2d Pw = P.block<2,2>(4, 4) / m1.noise_gain2();
    check((Pw - Eigen::Matrix2d::Identity()).cwiseAbs().maxCoeff() < 1e-9, "colouring Gramian = I",
          (Pw - Eigen::Matrix2d::Identity()).cwiseAbs().maxCoeff(), 1e-9);
}

// ---------------------------------------------------------------------------
// [2] realization against a test-local complexified mode.
//
// Per axis zeta = B/(P+ D) xi, xi = xi_r + i xi_i: B(s) = s + w_z0 and
// P = P+ P- the C41J polynomials (P+ the upper-half-plane poles), D the real
// Hurwitz quadratic with |D(iw)|^2 = (|P+(iw)|^2 + |P+(-iw)|^2) / 2.  Complex
// controllable canonical form of 1/(P+ D), realified to eight states
// [Re x; Im x]; the readouts p, v, a are Re(B x), Re(B x)', Re(B x)''
// (relative degree 3: noise-free rows).
struct Reference8 {
    Eigen::Matrix<double,8,8> A;
    Eigen::Matrix<double,8,2> G;
    Eigen::Matrix<double,3,8> C;
    double spectral_residual = 0.0;   // max relative error of the defining |D|^2 identity
};

static Reference8 reference8(const kc::CarmaShape& sh, double kappa) {
    const double w1 = 2 * PI * kappa * sh.f1_hz, w2 = 2 * PI * kappa * sh.f2_hz, wz0 = 2 * PI * kappa * sh.fz_hz;
    const Cd p1(-sh.zeta1 * w1, w1 * std::sqrt(1.0 - sh.zeta1 * sh.zeta1));
    const Cd p2(-sh.zeta2 * w2, w2 * std::sqrt(1.0 - sh.zeta2 * sh.zeta2));
    const Cd b1 = -(p1 + p2), b0 = p1 * p2;                     // P+(s) = s^2 + b1 s + b0
    auto pplus = [&](double w) { const Cd s(0.0, w); return s * s + b1 * s + b0; };
    auto half_q = [&](double w) { return 0.5 * (std::norm(pplus(w)) + std::norm(pplus(-w))); };
    // half_q(w) = w^4 + beta w^2 + gamma; fit at two frequencies, then factor
    // |D(iw)|^2 = w^4 + (d1^2 - 2 d0) w^2 + d0^2.
    const double u1 = w1 * w1, u2 = w2 * w2;
    const double r1 = half_q(w1) - u1 * u1, r2 = half_q(w2) - u2 * u2;
    const double beta = (r1 - r2) / (u1 - u2), gamma = r1 - beta * u1;
    const double d0 = std::sqrt(gamma), d1 = std::sqrt(beta + 2.0 * d0);
    Reference8 R;
    for (double w : {0.0, 0.2, 1.0, 2.5, 7.0, 30.0}) {
        const Cd s(0.0, w);
        R.spectral_residual = std::max(R.spectral_residual, std::abs(std::norm(s * s + d1 * s + d0) - half_q(w)) / half_q(w));
    }
    // P+ D = s^4 + a3 s^3 + a2 s^2 + a1 s + a0.
    const Cd a3 = b1 + d1, a2 = b0 + b1 * d1 + d0, a1 = b0 * d1 + b1 * d0, a0 = b0 * d0;
    Eigen::Matrix<Cd,4,4> Ac = Eigen::Matrix<Cd,4,4>::Zero();
    Ac(0,1) = Ac(1,2) = Ac(2,3) = 1.0;
    Ac(3,0) = -a0; Ac(3,1) = -a1; Ac(3,2) = -a2; Ac(3,3) = -a3;
    R.A << Ac.real(), -Ac.imag(), Ac.imag(), Ac.real();
    R.G.setZero();
    R.G(3,0) = 1.0;      // xi_r drives Re x_4
    R.G(7,1) = 1.0;      // xi_i drives Im x_4
    R.C.setZero();
    R.C(0,0) = wz0; R.C(0,1) = 1.0;   // Re(B x): B real
    R.C(1,1) = wz0; R.C(1,2) = 1.0;
    R.C(2,2) = wz0; R.C(2,3) = 1.0;
    return R;
}

template<int N>
static Eigen::Matrix<Cd,3,2> tf(const Eigen::Matrix<double,N,N>& A, const Eigen::Matrix<double,N,2>& Gin,
                                const Eigen::Matrix<double,3,N>& C, double w) {
    const Eigen::Matrix<Cd,N,N> M = Cd(0.0, w) * Eigen::Matrix<Cd,N,N>::Identity() - A.template cast<Cd>();
    return C.template cast<Cd>() * M.fullPivLu().solve(Gin.template cast<Cd>());
}

// A X + X A^T + Q = 0 (Kronecker; small N).
template<int N>
static Eigen::Matrix<double,N,N> lyap(const Eigen::Matrix<double,N,N>& A, const Eigen::Matrix<double,N,N>& Q) {
    constexpr Eigen::Index n = N;
    const Eigen::Matrix<double,N,N> I = Eigen::Matrix<double,N,N>::Identity();
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n * n, n * n);
    for (Eigen::Index i = 0; i < n; ++i)
        for (Eigen::Index j = 0; j < n; ++j) K.block(i * n, j * n, n, n) += A(i,j) * I;   // vec(A X), row-major vec
    for (Eigen::Index i = 0; i < n; ++i) K.block(i * n, i * n, n, n) += A;                  // vec(X A^T)
    Eigen::VectorXd q(n * n);
    for (Eigen::Index i = 0; i < n; ++i)
        for (Eigen::Index j = 0; j < n; ++j) q(i * n + j) = -Q(i,j);
    const Eigen::VectorXd x = K.fullPivLu().solve(q);
    Eigen::Matrix<double,N,N> X;
    for (Eigen::Index i = 0; i < n; ++i)
        for (Eigen::Index j = 0; j < n; ++j) X(i,j) = x(i * n + j);
    return X;
}

// exp(A h): scaling and squaring with a Taylor series (double).
template<int N>
static Eigen::Matrix<double,N,N> expm(const Eigen::Matrix<double,N,N>& A, double h) {
    int s = 0;
    double nrm = (A * h).cwiseAbs().rowwise().sum().maxCoeff();
    while (nrm > 0.1) { nrm *= 0.5; ++s; }
    const Eigen::Matrix<double,N,N> X = A * (h / std::ldexp(1.0, s));
    Eigen::Matrix<double,N,N> E = Eigen::Matrix<double,N,N>::Identity(), T = E;
    for (int k = 1; k <= 18; ++k) { T = T * X / static_cast<double>(k); E += T; }
    for (int i = 0; i < s; ++i) E = E * E;
    return E;
}

static kc::Hermitian3<double> test_sigma() {
    kc::Hermitian3<double> S;
    S.re << 0.30, 0.05, 0.10, 0.05, 0.12, 0.04, 0.10, 0.04, 0.40;
    S.im << 0.0, 0.0, 0.25, 0.0, 0.0, 0.10, -0.25, -0.10, 0.0;
    return S;
}

void test2_realization() {
    std::printf("[2] realization against a test-local complexified mode (eight states per axis)\n");
    double worst_tf = 0.0, worst_law = 0.0, worst_spec = 0.0;
    Eigen::Matrix<double,3,6> C6 = Eigen::Matrix<double,3,6>::Zero();
    C6.leftCols<3>().setIdentity();
    Eigen::Matrix2d J2;
    J2 << 0.0, -1.0, 1.0, 0.0;
    for (double kappa : {0.6, 1.0, 1.7}) {
        Motion m;
        kc::CarmaCouplingConfig cc;
        cc.source = kc::CarmaCouplingConfig::Source::External;
        cc.frame = kc::CarmaCouplingConfig::Frame::World;
        cc.warmup_s = 0.0f;
        cc.shrink = false;
        m.coupling().set_config(cc);
        m.coupling().set_external(test_sigma(), 0.0);
        m.set_vertical(kappa, std::sqrt(0.40));
        m.set_horizontal(kappa, std::sqrt(0.30), std::sqrt(0.12), 0.0);
        m.refresh(0.005);
        const Reference8 R = reference8(m.vertical_model().shape(), kappa);
        worst_spec = std::max(worst_spec, R.spectral_residual);
        for (double w : {0.3, 1.0, 1.9, 3.0, 6.0, 15.0}) {
            const auto H8 = tf<8>(R.A, R.G, R.C, w);
            const auto H6 = tf<6>(m.block_A(), m.block_G(), C6, w);
            worst_tf = std::max(worst_tf, (H8 - H6).cwiseAbs().maxCoeff() / H8.cwiseAbs().maxCoeff());
        }
        // Stationary law under the motion's own excitation and noise gain.
        const kc::Hermitian3<double>& S = m.excitation();
        const double g2 = m.noise_gain2();
        const Eigen::Matrix<double,8,8> WR = g2 * R.G * R.G.transpose(), WI = g2 * R.G * J2 * R.G.transpose();
        Eigen::Matrix<double,18,18> P6 = Eigen::Matrix<double,18,18>::Zero();
        m.stationary(P6, 0);
        const Eigen::Matrix<double,8,8> F8 = expm<8>(R.A, 0.005);
        const Eigen::Matrix<double,6,6> F6 = m.block_F();
        double var[3];
        for (int a = 0; a < 3; ++a) var[a] = (R.C * lyap<8>(R.A, S.re(a,a) * WR) * R.C.transpose())(0,0);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                const double si = a == b ? 0.0 : S.im(a,b);
                const Eigen::Matrix<double,8,8> P8 = lyap<8>(R.A, S.re(a,b) * WR + si * WI);
                const Eigen::Matrix<double,6,6> P6ab = P6.block<6,6>(Eigen::Index{6} * a, Eigen::Index{6} * b);
                const double sc = std::sqrt(var[a] * var[b]);
                const Eigen::Matrix3d l8 = R.C * P8 * R.C.transpose(), l6 = C6 * P6ab * C6.transpose();
                const Eigen::Matrix3d c8 = R.C * F8 * P8 * R.C.transpose(), c6 = C6 * F6 * P6ab * C6.transpose();
                // Relative to the p, v, a scales of the two axes.
                Eigen::Vector3d da, db;
                for (int i = 0; i < 3; ++i) {
                    da(i) = 1.0 / std::sqrt((R.C * lyap<8>(R.A, S.re(a,a) * WR) * R.C.transpose())(i,i));
                    db(i) = 1.0 / std::sqrt((R.C * lyap<8>(R.A, S.re(b,b) * WR) * R.C.transpose())(i,i));
                }
                worst_law = std::max(worst_law, (da.asDiagonal() * (l8 - l6) * db.asDiagonal()).cwiseAbs().maxCoeff());
                worst_law = std::max(worst_law, (c8 - c6).cwiseAbs().maxCoeff() / sc);
            }
    }
    check(worst_spec < 1e-12, "reference D: |D(iw)|^2 = (|P+(iw)|^2 + |P+(-iw)|^2) / 2", worst_spec, 1e-12);
    check(worst_tf < 1e-9, "noise -> (p, v, a) transfer functions equal (kappa 0.6, 1, 1.7)", worst_tf, 1e-9);
    check(worst_law < 1e-8, "stationary 3-axis law of (p, v, a) and its lag-h cross-covariance equal (Hermitian excitation)",
          worst_law, 1e-8);
}

// ---------------------------------------------------------------------------
// Core-level paired runs (double, no wrapper: the schedule below commits the
// models every 0.1 s).
struct Commit { double kappa, sx, sy, sz; };

template<class F>
static void core_init(F& f, const V3d& acc, const V3d& mag, const Commit& c) {
    f.initialize_from_acc_mag(acc, mag);
    f.set_vertical_carma(c.kappa, c.sz);
    f.set_horizontal_carma(c.kappa, c.sx, c.sy, 0.0);
    f.reset_aw_covariance_to_stationary();
}
template<class F>
static void core_commit(F& f, const Commit& c) {
    f.set_vertical_carma(c.kappa, c.sz);
    f.set_horizontal_carma(c.kappa, c.sx, c.sy, 0.0);
}

struct PairDiff {
    double p = 0, v = 0, att = 0, bg = 0, ba = 0, cov = 0;
    bool finite = true;
};

template<class F1, class F2, class Sched>
static PairDiff run_pair(F1& f1, F2& f2, double T, Sched sched, unsigned seed) {
    Sea sea(30.0, seed);
    Imu imu;
    imu.rng.seed(seed + 100);
    PairDiff d;
    const double dt = 0.005;
    const long steps = static_cast<long>(T / dt);
    V3d p, a, acc, gyr;
    sea.eval(0.0, p, a);
    imu.sample(a, acc, gyr);
    const V3d mag0 = imu.mag();
    core_init(f1, acc, mag0, sched(0.0));
    core_init(f2, acc, mag0, sched(0.0));
    for (long k = 1; k < steps; ++k) {
        const double t = static_cast<double>(k) * dt;
        if (k % 20 == 0) { const Commit c = sched(t); core_commit(f1, c); core_commit(f2, c); }
        sea.eval(t, p, a);
        imu.sample(a, acc, gyr);
        f1.time_update(gyr, dt); f2.time_update(gyr, dt);
        f1.measurement_update_acc_only(acc, 35.0); f2.measurement_update_acc_only(acc, 35.0);
        if (k % 8 == 0) { const V3d m = imu.mag(); f1.measurement_update_mag_only(m); f2.measurement_update_mag_only(m); }
        const double sp = std::max(1e-3, f1.get_position().norm());
        d.p = std::max(d.p, (f1.get_position() - f2.get_position()).norm() / sp);
        d.v = std::max(d.v, (f1.get_velocity() - f2.get_velocity()).norm() / std::max(1e-3, f1.get_velocity().norm()));
        d.att = std::max(d.att, f1.quaternion_boat().angularDistance(f2.quaternion_boat()));
        d.bg = std::max(d.bg, (f1.gyroscope_bias() - f2.gyroscope_bias()).norm());
        d.ba = std::max(d.ba, (f1.get_acc_bias() - f2.get_acc_bias()).norm());
        const Eigen::Matrix<double,6,6> C1 = f1.pv_covariance6(), C2 = f2.pv_covariance6();
        Eigen::Matrix<double,6,1> s = C1.diagonal().cwiseSqrt().cwiseInverse();
        d.cov = std::max(d.cov, (s.asDiagonal() * (C1 - C2) * s.asDiagonal()).cwiseAbs().maxCoeff());
        d.finite = d.finite && f2.covariance_full().allFinite();
    }
    return d;
}

using Core = csm::KalmanCarmaSharedModes<double>;
using CoreA41 = ocean_imu::kalman::Kalman3D_Wave_CARMA3<double, 4, false, 4, false, true, true>;

static const V3d kSa = V3d::Constant(2.8 * 0.0148 * 0.71), kSg = V3d::Constant(2.0 * 0.00157 * 0.01),
                 kSm = V3d::Constant(1.2 * 0.8 * 8.0);

static void configure_sm(Core& f, csm::CouplingMode mode, bool confidence) {
    csm::SharedModesConfig<double> c;
    c.mode = mode;
    c.confidence.enabled = confidence;
    c.confidence.accel_noise_std = 0.0148;
    c.warmup_s = 60.0;
    c.horizon_s = 120.0;
    f.configureSharedModes(c);
}

// ---------------------------------------------------------------------------
// [3] coupling off = A41 under adaptation
void test3_off_a41() {
    std::printf("[3] coupling off = CARMA3 baseline A41 through kappa and amplitude adaptation (double)\n");
    auto sched = [](double t) {
        const double k = 1.0 + 0.3 * std::sin(2 * PI * t / 150.0) + 0.02 * std::sin(2 * PI * t / 2.3);
        const double m = 1.0 + 0.4 * std::sin(2 * PI * t / 80.0) + 0.1 * std::sin(2 * PI * t / 1.7);
        return Commit{k, 0.5 * m, 0.3 * (2.0 - m), 0.6 * m};
    };
    auto fa = std::make_unique<CoreA41>(kSa, kSg, kSm);
    auto fc = std::make_unique<Core>(kSa, kSg, kSm);
    configure_sm(*fc, csm::CouplingMode::Off, false);
    const PairDiff d = run_pair(*fa, *fc, 400.0, sched, 5);
    check(d.finite && d.p < 1e-7 && d.v < 1e-7, "positions and velocities equal A41's", std::max(d.p, d.v), 1e-7);
    check(std::max({d.att, d.bg, d.cov}) < 1e-7, "attitude, gyro bias, p/v covariance equal A41's",
          std::max({d.att, d.bg, d.cov}), 1e-7);
}

// ---------------------------------------------------------------------------
// [4] commits
void test4_commits() {
    std::printf("[4] commits\n");
    auto f = std::make_unique<Core>(kSa, kSg, kSm);
    configure_sm(*f, csm::CouplingMode::Estimated, false);
    Sea sea(30.0, 9);
    Imu imu;
    V3d p, a, acc, gyr;
    sea.eval(0.0, p, a);
    imu.sample(a, acc, gyr);
    core_init(*f, acc, imu.mag(), Commit{1.0, 0.5, 0.3, 0.6});
    using MB = Core::Motion::MatB;
    double jump = 0.0, worst_eig = 1.0, row_dev = 0.0;
    long nontrivial = 0, commits = 0;
    for (long k = 1; k < 200L * 200L; ++k) {
        const double t = static_cast<double>(k) * 0.005;
        sea.eval(t, p, a);
        imu.sample(a, acc, gyr);
        if (k % 20 == 0) {
            const Commit c{1.0 + 0.3 * std::sin(t / 7.0), 0.5 * (1.2 + std::sin(t / 5.0) * 0.3), 0.3, 0.6};
            // A zero-length prediction returns before the commit is applied, so
            // exercise the core's own path on a copy: reset_aw_covariance_to_stationary()
            // runs refresh() and recoordinate_motion_() without propagating the
            // state (it replaces only the motion covariance).
            const Core& cf = *f;
            auto g = std::make_unique<Core>(cf);
            core_commit(*g, c);
            g->reset_aw_covariance_to_stationary();
            jump = std::max({jump, (g->get_position() - f->get_position()).norm(),
                             (g->get_velocity() - f->get_velocity()).norm(),
                             (g->get_world_accel() - f->get_world_accel()).norm()});
            // ... and the transform that path applied is a real one: rows p, v, a
            // are unit rows, the colouring rows move.
            auto h = std::make_unique<Core>(cf);
            core_commit(*h, c);
            h->motion().refresh(0.005);
            MB M[3];
            if (h->motion().take_axis_recoordinatization(M)) {
                ++nontrivial;
                for (const MB& m : M)
                    row_dev = std::max(row_dev, (m.topRows<3>() - MB::Identity().topRows<3>()).cwiseAbs().maxCoeff());
            }
            ++commits;
            core_commit(*f, c);
        }
        f->time_update(gyr, 0.005);
        f->measurement_update_acc_only(acc, 35.0);
        if (k % 4000 == 0) {
            const Eigen::MatrixXd P = f->covariance_full();
            const Eigen::VectorXd s = P.diagonal().cwiseSqrt().cwiseInverse();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s.asDiagonal() * P * s.asDiagonal());
            worst_eig = std::min(worst_eig, es.eigenvalues()(0));
        }
    }
    char buf[200];
    std::snprintf(buf, sizeof buf, "p, v, a continuous through the core's commit path (%ld of %ld commits re-coordinatize)",
                  nontrivial, commits);
    check(jump < 1e-12, buf, jump, 1e-12);
    check(row_dev == 0.0, "the applied transforms keep the p, v, a rows exactly", row_dev, 0.0);
    check(nontrivial > commits / 2, "amplitude and kappa commits re-coordinatize the colouring",
          static_cast<double>(nontrivial), static_cast<double>(commits));
    check(worst_eig > -1e-9, "covariance PSD after 2000 commits (correlation min eigenvalue)", worst_eig, -1e-9);
}

// ---------------------------------------------------------------------------
// [5] float wrapper
template<class Fusion>
static void begin_wrapper(Fusion& f) {
    typename Fusion::Config c;
    c.with_mag = true;
    c.sigma_a = V3f::Constant(2.8f * 0.0148f * 0.71f);
    c.sigma_g = V3f::Constant(2.0f * 0.00157f * 0.01f);
    c.sigma_m = V3f::Constant(1.2f * 0.8f * 8.0f);
    c.mag_delay_sec = 0.0f;
    if constexpr (requires { typename Fusion::SharedModes; }) {
        typename Fusion::SharedModes sm;
        sm.confidence.accel_noise_std = 0.0148f;
        f.begin(c, sm);
    } else {
        f.begin(c);
    }
    f.raw().enableTuner(true);
    f.raw().enableClamp(true);
}

void test5_wrapper() {
    std::printf("[5] float wrapper on a synthetic directional sea (600 s, scored from 300 s)\n");
    using A41 = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false>;
    auto wa = std::make_unique<A41>();
    auto wc = std::make_unique<csm::SeaStateFusion_CarmaSharedModes>();
    begin_wrapper(*wa);
    begin_wrapper(*wc);
    Sea sea(30.0, 1);
    Imu imu;
    double ea = 0, ec = 0, worst_eig = 1.0, asym = 0.0;
    long n = 0;
    bool finite = true;
    for (long k = 0; k < 600L * 200L; ++k) {
        const double t = static_cast<double>(k) * 0.005;
        V3d p, a, acc, gyr;
        sea.eval(t, p, a);
        imu.sample(a, acc, gyr);
        const V3f af = acc.cast<float>(), gf = gyr.cast<float>();
        wa->update(0.005f, gf, af, 35.0f);
        wc->update(0.005f, gf, af, 35.0f);
        if (k % 8 == 0) {
            const V3f m = imu.mag().cast<float>();
            wa->updateMag(m); wc->updateMag(m);
        }
        if (t >= 300.0 && wc->isLive()) {
            ea += (p - wa->raw().mekf().get_position().cast<double>()).squaredNorm();
            ec += (p - wc->raw().mekf().get_position().cast<double>()).squaredNorm();
            ++n;
            if (k % 6000 == 0) {
                const Eigen::MatrixXd P = wc->raw().mekf().covariance_full().cast<double>();
                finite = finite && P.allFinite();
                asym = std::max(asym, (P - P.transpose()).cwiseAbs().maxCoeff());
                const Eigen::VectorXd s = P.diagonal().cwiseSqrt().cwiseInverse();
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s.asDiagonal() * P * s.asDiagonal());
                worst_eig = std::min(worst_eig, es.eigenvalues()(0));
            }
        }
    }
    const double nn = static_cast<double>(std::max(1L, n));
    ea = std::sqrt(ea / nn); ec = std::sqrt(ec / nn);
    std::printf("  info: 3D position RMS A41 %.2f mm, shared modes %.2f mm; c = %.3f\n", 1e3 * ea, 1e3 * ec,
                static_cast<double>(wc->raw().mekf().couplingConfidence()));
    check(n > 0 && finite && asym == 0.0, "covariance finite and exactly symmetric", asym, 0.0);
    check(worst_eig > -1e-4, "covariance PSD (float; correlation min eigenvalue)", worst_eig, -1e-4);
    check(wc->raw().mekf().couplingConfidence() > 0.9f, "the estimated coupling engages (c > 0.9)");
    check(ec < ea, "position RMS below the independent-axis CARMA3 wrapper's on a coupled sea", ec / ea, 1.0);
}

int main() {
    test1_structure();
    test2_realization();
    test3_off_a41();
    test4_commits();
    test5_wrapper();
    std::printf("carma_shared_modes_realization-test: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures ? 1 : 0;
}
