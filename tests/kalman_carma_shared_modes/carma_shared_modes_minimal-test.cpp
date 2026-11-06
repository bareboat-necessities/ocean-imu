/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    carma_shared_modes_minimal-test -- unit and numerical tests of the
    EXPERIMENTAL 27-state minimal realization of the shared-modes QMEKF
    (src/kalman_carma_shared_modes/KalmanCarmaSharedModesMinimal.h,
    docs/kalman-carma-shared-modes-minimal.md section 5).

    Groups:
      [1] structure: 27 error states, six per axis, valid realizations,
          lossless colouring, kappa scaling
      [2] realization: per-axis transfer functions from the two noise inputs
          to p, v, a equal the complexified 33-state block's; the stationary
          three-axis law of (p, v, a) under a Hermitian excitation equals the
          33-state block's (co- and quadrature spectra)
      [3] filter level at fixed kappa (double): every realization reproduces
          the 33-state QMEKF (positions, velocities, attitude, biases, the
          p/v covariance) through estimated coupling, confidence and fallback
          shrinks; the quotient also through amplitude commits
      [4] cascade with the coupling off reproduces the CARMA3 baseline A41
          through kappa and amplitude adaptation (double)
      [5] commits: p, v, a are continuous through every re-coordinatization;
          the quotient lift is the identity without a kappa change; the
          lifted covariance stays positive semidefinite
      [6] float wrapper on a synthetic directional sea: finite, symmetric,
          PSD covariance; the coupling engages; accuracy against the 33-state
          wrapper on the same stream
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModesMinimal.h"

extern const float g_std = 9.80665f;   // the application constant the front end reads

namespace csm = ocean_imu::carma_shared_modes;
namespace kc = ocean_imu::kalman::carma;
using MR = csm::MinimalRealization;
using V3d = Eigen::Vector3d;
using V3f = Eigen::Vector3f;
using Cd = std::complex<double>;

static int g_checks = 0, g_failures = 0;
static void check(bool ok, const char* what, double value = NAN, double limit = NAN) {
    ++g_checks;
    if (!ok) ++g_failures;
    std::printf("  %s: %s (value=%.4g limit=%.4g)\n", ok ? "ok" : "FAIL", what, value, limit);
}

constexpr double G = 9.80665;
constexpr double PI = 3.14159265358979323846;

static const char* rname(MR r) {
    return r == MR::Cascade ? "Cascade" : r == MR::Quotient ? "Quotient" : "QuotientHold";
}

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
template<MR R>
void test1_one() {
    using M = csm::CarmaMinimalSharedModeMotion<double, R>;
    M m;
    check(m.valid(), (std::string(rname(R)) + ": realization valid").c_str());
    Eigen::Matrix<double,6,6> A = m.block_A();
    Eigen::EigenSolver<Eigen::Matrix<double,6,6>> es(A);
    check(es.eigenvalues().real().maxCoeff() < 0.0, (std::string(rname(R)) + ": block Hurwitz").c_str(),
          es.eigenvalues().real().maxCoeff(), 0.0);
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
    check(d < 1e-10, (std::string(rname(R)) + ": poles scale with kappa").c_str(), d, 1e-10);
}

void test1_structure() {
    std::printf("[1] structure\n");
    check(csm::KalmanCarmaSharedModesMinimal<float>::kStates == 27 &&
          csm::KalmanCarmaSharedModesMinimal<double, MR::Quotient>::kStates == 27 &&
          csm::KalmanCarmaSharedModes<float>::kStates == 33,
          "27 error states (3 + 3 + 18 + 3); the 33-state core is unchanged");
    test1_one<MR::Cascade>();
    test1_one<MR::Quotient>();
    test1_one<MR::QuotientHold>();
    // Lossless colouring of the cascade: the stationary covariance of the two
    // colouring states is g2 sigma^2 I (input-normal coordinates).
    csm::CarmaMinimalSharedModeMotion<double, MR::Cascade> m;
    Eigen::Matrix<double,18,18> P = Eigen::Matrix<double,18,18>::Zero();
    m.stationary(P, 0);
    const Eigen::Matrix2d Pw = P.block<2,2>(4, 4) / m.noise_gain2();
    check((Pw - Eigen::Matrix2d::Identity()).cwiseAbs().maxCoeff() < 1e-9, "cascade colouring Gramian = I",
          (Pw - Eigen::Matrix2d::Identity()).cwiseAbs().maxCoeff(), 1e-9);
}

// ---------------------------------------------------------------------------
// [2] realization against the complexified block
template<int N>
static Eigen::Matrix<Cd,3,2> tf_pva(const Eigen::Matrix<double,N,N>& A, const Eigen::Matrix<double,N,2>& Gin, double w) {
    const Eigen::Matrix<Cd,N,N> M = Cd(0.0, w) * Eigen::Matrix<Cd,N,N>::Identity() - A.template cast<Cd>();
    const Eigen::Matrix<Cd,N,2> X = M.fullPivLu().solve(Gin.template cast<Cd>());
    return X.template topRows<3>();
}

static kc::Hermitian3<double> test_sigma() {
    kc::Hermitian3<double> S;
    S.re << 0.30, 0.05, 0.10, 0.05, 0.12, 0.04, 0.10, 0.04, 0.40;
    S.im << 0.0, 0.0, 0.25, 0.0, 0.0, 0.10, -0.25, -0.10, 0.0;
    return S;
}

template<MR R>
void test2_one(double& worst_tf, double& worst_law) {
    for (double kappa : {0.6, 1.0, 1.7}) {
        kc::CarmaSharedModeMotion<double> m8;
        csm::CarmaMinimalSharedModeMotion<double, R> m6;
        auto setup = [&](auto& m) {
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
        };
        setup(m8);
        setup(m6);
        Eigen::Matrix<double,8,2> G8 = Eigen::Matrix<double,8,2>::Zero();
        G8(2,0) = 1.0; G8(6,1) = 1.0;
        const Eigen::Matrix<double,6,2> G6 = m6.block_G();
        for (double w : {0.3, 1.0, 1.9, 3.0, 6.0, 15.0}) {
            const auto H8 = tf_pva<8>(m8.block_A(), G8, w);
            const auto H6 = tf_pva<6>(m6.block_A(), G6, w);
            worst_tf = std::max(worst_tf, (H8 - H6).cwiseAbs().maxCoeff() / H8.cwiseAbs().maxCoeff());
        }
        Eigen::Matrix<double,24,24> P8 = Eigen::Matrix<double,24,24>::Zero();
        Eigen::Matrix<double,18,18> P6 = Eigen::Matrix<double,18,18>::Zero();
        m8.stationary(P8, 0);
        m6.stationary(P6, 0);
        Eigen::Matrix<double,9,9> L8, L6;
        for (Eigen::Index a = 0; a < 3; ++a)
            for (Eigen::Index i = 0; i < 3; ++i)
                for (Eigen::Index b = 0; b < 3; ++b)
                    for (Eigen::Index j = 0; j < 3; ++j) {
                        L8(3 * a + i, 3 * b + j) = P8(8 * a + i, 8 * b + j);
                        L6(3 * a + i, 3 * b + j) = P6(6 * a + i, 6 * b + j);
                    }
        Eigen::Matrix<double,9,1> d = L8.diagonal().cwiseSqrt().cwiseInverse();
        worst_law = std::max(worst_law, (d.asDiagonal() * (L8 - L6) * d.asDiagonal()).cwiseAbs().maxCoeff());
        // One step: the lag-h cross-covariance F P (quadrature sign) as well.
        const Eigen::Matrix<double,8,8> F8 = m8.block_F();
        const Eigen::Matrix<double,6,6> F6 = m6.block_F();
        for (Eigen::Index a = 0; a < 3; ++a)
            for (Eigen::Index b = 0; b < 3; ++b) {
                const Eigen::Matrix<double,3,3> c8 = (F8 * P8.block<8,8>(8 * a, 8 * b)).topLeftCorner<3,3>();
                const Eigen::Matrix<double,3,3> c6 = (F6 * P6.block<6,6>(6 * a, 6 * b)).topLeftCorner<3,3>();
                const double s = std::sqrt(P8(8 * a, 8 * a) * P8(8 * b, 8 * b));
                worst_law = std::max(worst_law, (c8 - c6).cwiseAbs().maxCoeff() / s);
            }
    }
}

void test2_realization() {
    std::printf("[2] realization against the complexified 33-state block\n");
    double tfc = 0, lc = 0, tfq = 0, lq = 0;
    test2_one<MR::Cascade>(tfc, lc);
    test2_one<MR::Quotient>(tfq, lq);
    check(tfc < 1e-9, "cascade: noise -> (p, v, a) transfer functions equal (kappa 0.6, 1, 1.7)", tfc, 1e-9);
    check(tfq < 1e-9, "quotient: noise -> (p, v, a) transfer functions equal", tfq, 1e-9);
    check(lc < 1e-8, "cascade: stationary 3-axis law of (p, v, a) and its lag-h cross-covariance equal (Hermitian excitation)", lc, 1e-8);
    check(lq < 1e-8, "quotient: stationary 3-axis law equal", lq, 1e-8);
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
    double p = 0, v = 0, att = 0, bg = 0, ba = 0, cov = 0, p_rms = 0;
    bool finite = true;
};

struct NoHook { template<class F> void operator()(F&, double) const {} };

template<class F1, class F2, class Sched, class Hook = NoHook>
static PairDiff run_pair(F1& f1, F2& f2, double T, Sched sched, unsigned seed, Hook hook = Hook{}) {
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
    double se = 0.0;
    long n = 0;
    for (long k = 1; k < steps; ++k) {
        const double t = static_cast<double>(k) * dt;
        if (k % 20 == 0) { const Commit c = sched(t); core_commit(f1, c); core_commit(f2, c); }
        sea.eval(t, p, a);
        imu.sample(a, acc, gyr);
        hook(f1, t); hook(f2, t);
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
        if (t > 60.0) { se += (f1.get_position() - p).squaredNorm(); ++n; }
    }
    d.p_rms = std::sqrt(se / static_cast<double>(std::max(1L, n)));
    return d;
}

using Core33 = csm::KalmanCarmaSharedModes<double>;
template<MR R> using Core27 = csm::KalmanCarmaSharedModesMinimal<double, R>;
using CoreA41 = ocean_imu::kalman::Kalman3D_Wave_CARMA3<double, 4, false, 4, false, true, true>;

static const V3d kSa = V3d::Constant(2.8 * 0.0148 * 0.71), kSg = V3d::Constant(2.0 * 0.00157 * 0.01),
                 kSm = V3d::Constant(1.2 * 0.8 * 8.0);

template<class F>
static void configure_sm(F& f, csm::CouplingMode mode, bool confidence) {
    csm::SharedModesConfig<double> c;
    c.mode = mode;
    c.confidence.enabled = confidence;
    c.confidence.accel_noise_std = 0.0148;
    c.warmup_s = 60.0;
    c.horizon_s = 120.0;
    f.configureSharedModes(c);
}

// ---------------------------------------------------------------------------
// [3] fixed kappa
template<MR R>
void test3_one(bool amplitude_commits, double tol) {
    auto sched = [&](double t) {
        Commit c{1.05, 0.55, 0.35, 0.65};
        if (amplitude_commits) {
            const double m = 1.0 + 0.25 * std::sin(2 * PI * t / 97.0) + 0.05 * std::sin(2 * PI * t / 3.1);
            c.sx *= m; c.sy *= 1.0 + 0.2 * std::sin(2 * PI * t / 61.0); c.sz *= m;
        }
        return c;
    };
    auto f33 = std::make_unique<Core33>(kSa, kSg, kSm);
    auto f27 = std::make_unique<Core27<R>>(kSa, kSg, kSm);
    configure_sm(*f33, csm::CouplingMode::Estimated, true);
    configure_sm(*f27, csm::CouplingMode::Estimated, true);
    // A forced confidence drop at 250 s (released at 300 s) exercises the
    // PSD-preserving fallback shrink of the horizontal-vertical covariance.
    auto hook = [](auto& f, double t) {
        if (std::abs(t - 250.0) < 1e-9) f.forceCouplingConfidence(0.25);
        if (std::abs(t - 300.0) < 1e-9) f.forceCouplingConfidence(-1.0);
    };
    const PairDiff d = run_pair(*f33, *f27, 400.0, sched, 3, hook);
    const auto st = f33->coupling();
    char buf[256];
    std::snprintf(buf, sizeof buf, "%s, fixed kappa%s: 27 = 33 (rel. p, v; c = %.2f, %lu shrinks)", rname(R),
                  amplitude_commits ? " + amplitude commits" : "", static_cast<double>(st.confidence.confidence),
                  st.shrink_events);
    check(d.finite && d.p < tol && d.v < tol, buf, std::max(d.p, d.v), tol);
    std::snprintf(buf, sizeof buf, "%s: attitude, biases and p/v covariance equal", rname(R));
    const double o = std::max({d.att, d.bg, d.ba * 1e-2, d.cov});
    check(o < tol, buf, o, tol);
}

void test3_fixed_kappa() {
    std::printf("[3] 27 = 33 at fixed kappa (double cores, estimated coupling, confidence, fallback)\n");
    test3_one<MR::Cascade>(false, 1e-7);
    test3_one<MR::Quotient>(false, 1e-7);
    test3_one<MR::QuotientHold>(false, 1e-7);
    test3_one<MR::Quotient>(true, 1e-7);
    test3_one<MR::QuotientHold>(true, 1e-7);
    // The cascade holds the NORMALIZED colouring through amplitude commits, so
    // there it is not the 33-state filter (by design; [4] is its invariant).
    auto sched = [](double t) {
        const double m = 1.0 + 0.25 * std::sin(2 * PI * t / 97.0);
        return Commit{1.05, 0.55 * m, 0.35, 0.65 * m};
    };
    auto f33 = std::make_unique<Core33>(kSa, kSg, kSm);
    auto f27 = std::make_unique<Core27<MR::Cascade>>(kSa, kSg, kSm);
    configure_sm(*f33, csm::CouplingMode::Estimated, true);
    configure_sm(*f27, csm::CouplingMode::Estimated, true);
    const PairDiff d = run_pair(*f33, *f27, 400.0, sched, 3);
    check(d.finite && d.p > 1e-6 && d.p < 0.2, "cascade with amplitude commits: differs from 33 (different commit law), boundedly",
          d.p, 0.2);
}

// ---------------------------------------------------------------------------
// [4] cascade, coupling off = A41 under adaptation
void test4_cascade_a41() {
    std::printf("[4] cascade with the coupling off = CARMA3 baseline A41 through adaptation (double)\n");
    auto sched = [](double t) {
        const double k = 1.0 + 0.3 * std::sin(2 * PI * t / 150.0) + 0.02 * std::sin(2 * PI * t / 2.3);
        const double m = 1.0 + 0.4 * std::sin(2 * PI * t / 80.0) + 0.1 * std::sin(2 * PI * t / 1.7);
        return Commit{k, 0.5 * m, 0.3 * (2.0 - m), 0.6 * m};
    };
    auto fa = std::make_unique<CoreA41>(kSa, kSg, kSm);
    auto fc = std::make_unique<Core27<MR::Cascade>>(kSa, kSg, kSm);
    configure_sm(*fc, csm::CouplingMode::Off, false);
    const PairDiff d = run_pair(*fa, *fc, 400.0, sched, 5);
    check(d.finite && d.p < 1e-7 && d.v < 1e-7, "positions and velocities equal A41's", std::max(d.p, d.v), 1e-7);
    check(std::max({d.att, d.bg, d.cov}) < 1e-7, "attitude, gyro bias, p/v covariance equal A41's",
          std::max({d.att, d.bg, d.cov}), 1e-7);
    // ... while the 33-state and the quotient realizations of the same law are not
    // (the realization effect of docs/kalman-carma-shared-modes.md section 2).
    auto fa2 = std::make_unique<CoreA41>(kSa, kSg, kSm);
    auto f33 = std::make_unique<Core33>(kSa, kSg, kSm);
    configure_sm(*f33, csm::CouplingMode::Off, false);
    const PairDiff d33 = run_pair(*fa2, *f33, 400.0, sched, 5);
    auto fa3 = std::make_unique<CoreA41>(kSa, kSg, kSm);
    auto fq = std::make_unique<Core27<MR::Quotient>>(kSa, kSg, kSm);
    configure_sm(*fq, csm::CouplingMode::Off, false);
    const PairDiff dq = run_pair(*fa3, *fq, 400.0, sched, 5);
    std::printf("  info: max relative |p - p_A41|: 33-state %.3g, quotient %.3g, cascade %.3g\n", d33.p, dq.p, d.p);
    check(d33.p > 1e-4 && dq.p > 1e-4, "33-state and quotient differ from A41 under adaptation (coupling off)",
          std::min(d33.p, dq.p), 1e-4);
}

// ---------------------------------------------------------------------------
// [5] commits
template<MR R>
void test5_one() {
    auto f = std::make_unique<Core27<R>>(kSa, kSg, kSm);
    configure_sm(*f, csm::CouplingMode::Estimated, false);
    Sea sea(30.0, 9);
    Imu imu;
    V3d p, a, acc, gyr;
    sea.eval(0.0, p, a);
    imu.sample(a, acc, gyr);
    core_init(*f, acc, imu.mag(), Commit{1.0, 0.5, 0.3, 0.6});
    using Core = Core27<R>;
    using MB = typename Core::Motion::MatB;
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
            // are unit rows, the rest moves on Cascade / Quotient commits.
            auto h = std::make_unique<Core>(cf);
            core_commit(*h, c);
            h->motion().refresh(0.005);
            MB M[3];
            if (h->motion().take_axis_recoordinatization(M)) {
                ++nontrivial;
                for (const MB& m : M) {
                    row_dev = std::max(row_dev, (m.template topRows<3>() - MB::Identity().template topRows<3>()).cwiseAbs().maxCoeff());
                }
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
    std::snprintf(buf, sizeof buf, "%s: p, v, a continuous through the core's commit path (%ld of %ld commits re-coordinatize)",
                  rname(R), nontrivial, commits);
    check(jump < 1e-12, buf, jump, 1e-12);
    // Cascade: exact (diagonal scaling of the colouring rows).  Quotient: the
    // rows of Re(p, v, a, l) are [I 0] analytically; the numerically computed
    // regression lift keeps them to rounding.
    std::snprintf(buf, sizeof buf, "%s: the applied transforms keep the p, v, a rows (exact; quotient to rounding)", rname(R));
    const double row_tol = R == MR::Quotient ? 1e-13 : 0.0;
    check(row_dev <= row_tol, buf, row_dev, row_tol);
    std::snprintf(buf, sizeof buf, "%s: commits re-coordinatize as the realization requires", rname(R));
    check(R == MR::QuotientHold ? nontrivial == 0 : nontrivial > commits / 2, buf, static_cast<double>(nontrivial),
          static_cast<double>(commits));
    std::snprintf(buf, sizeof buf, "%s: covariance PSD after 2000 commits (correlation min eigenvalue)", rname(R));
    check(worst_eig > -1e-9, buf, worst_eig, -1e-9);
}

void test5_commits() {
    std::printf("[5] commits\n");
    test5_one<MR::Cascade>();
    test5_one<MR::Quotient>();
    test5_one<MR::QuotientHold>();
    // Quotient lift without a kappa change: no re-coordinatization at all.
    csm::CarmaMinimalSharedModeMotion<double, MR::Quotient> m;
    m.set_vertical(1.2, 1.0); m.set_horizontal(1.2, 1.0, 1.0, 0.0); m.refresh(0.005);
    Eigen::Matrix<double,6,6> M[3];
    (void)m.take_axis_recoordinatization(M);
    m.set_vertical(1.2, 1.3); m.set_horizontal(1.2, 0.7, 1.1, 0.0); m.refresh(0.005);
    check(!m.take_axis_recoordinatization(M), "quotient: an amplitude-only commit does not re-coordinatize");
    m.set_vertical(1.3, 1.3); m.set_horizontal(1.3, 0.7, 1.1, 0.0); m.refresh(0.005);
    const bool have = m.take_axis_recoordinatization(M);
    const double top = (M[0].topRows<4>() - Eigen::Matrix<double,4,6>::Identity()).cwiseAbs().maxCoeff();
    check(have && top < 1e-12, "quotient lift keeps Re(p, v, a, l) exactly (top rows [I 0])", top, 1e-12);
    // Lift-project at kappa_new == kappa_old would be the identity: Pi L = I.
    Eigen::Matrix<double,6,8> Pi = m.quotient_map();
    check(Pi.allFinite() && Pi.topLeftCorner<4,4>().isIdentity(1e-12), "quotient map keeps the real half");
}

// ---------------------------------------------------------------------------
// [6] float wrapper
template<class Fusion>
struct WrapRun {
    std::unique_ptr<Fusion> f = std::make_unique<Fusion>();
    WrapRun() {
        typename Fusion::Config c;
        c.with_mag = true;
        c.sigma_a = V3f::Constant(2.8f * 0.0148f * 0.71f);
        c.sigma_g = V3f::Constant(2.0f * 0.00157f * 0.01f);
        c.sigma_m = V3f::Constant(1.2f * 0.8f * 8.0f);
        c.mag_delay_sec = 0.0f;
        typename Fusion::SharedModes sm;
        sm.confidence.accel_noise_std = 0.0148f;
        f->begin(c, sm);
        f->raw().enableTuner(true);
        f->raw().enableClamp(true);
    }
};

void test6_wrapper() {
    std::printf("[6] float wrapper on a synthetic directional sea (600 s, scored from 300 s)\n");
    WrapRun<csm::SeaStateFusion_CarmaSharedModes> w33;
    WrapRun<csm::SeaStateFusion_CarmaSharedModesMinimal<MR::Cascade>> wc;
    WrapRun<csm::SeaStateFusion_CarmaSharedModesMinimal<MR::Quotient>> wq;
    Sea sea(30.0, 1);
    Imu imu;
    double e33 = 0, ec = 0, eq = 0, worst_eig = 1.0, asym = 0.0;
    long n = 0;
    bool finite = true;
    for (long k = 0; k < 600L * 200L; ++k) {
        const double t = static_cast<double>(k) * 0.005;
        V3d p, a, acc, gyr;
        sea.eval(t, p, a);
        imu.sample(a, acc, gyr);
        const V3f af = acc.cast<float>(), gf = gyr.cast<float>();
        w33.f->update(0.005f, gf, af, 35.0f);
        wc.f->update(0.005f, gf, af, 35.0f);
        wq.f->update(0.005f, gf, af, 35.0f);
        if (k % 8 == 0) {
            const V3f m = imu.mag().cast<float>();
            w33.f->updateMag(m); wc.f->updateMag(m); wq.f->updateMag(m);
        }
        if (t >= 300.0 && wc.f->isLive()) {
            e33 += (p - w33.f->raw().mekf().get_position().cast<double>()).squaredNorm();
            ec += (p - wc.f->raw().mekf().get_position().cast<double>()).squaredNorm();
            eq += (p - wq.f->raw().mekf().get_position().cast<double>()).squaredNorm();
            ++n;
            if (k % 6000 == 0) {
                for (const Eigen::MatrixXd& P : {Eigen::MatrixXd(wc.f->raw().mekf().covariance_full().cast<double>()),
                                                Eigen::MatrixXd(wq.f->raw().mekf().covariance_full().cast<double>())}) {
                    finite = finite && P.allFinite();
                    asym = std::max(asym, (P - P.transpose()).cwiseAbs().maxCoeff());
                    const Eigen::VectorXd s = P.diagonal().cwiseSqrt().cwiseInverse();
                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s.asDiagonal() * P * s.asDiagonal());
                    worst_eig = std::min(worst_eig, es.eigenvalues()(0));
                }
            }
        }
    }
    const double nn = static_cast<double>(std::max(1L, n));
    e33 = std::sqrt(e33 / nn); ec = std::sqrt(ec / nn); eq = std::sqrt(eq / nn);
    std::printf("  info: 3D position RMS 33-state %.2f mm, cascade %.2f mm, quotient %.2f mm; c = %.3f / %.3f / %.3f\n",
                1e3 * e33, 1e3 * ec, 1e3 * eq, static_cast<double>(w33.f->raw().mekf().couplingConfidence()),
                static_cast<double>(wc.f->raw().mekf().couplingConfidence()),
                static_cast<double>(wq.f->raw().mekf().couplingConfidence()));
    check(n > 0 && finite && asym == 0.0, "covariances finite and exactly symmetric", asym, 0.0);
    check(worst_eig > -1e-4, "covariances PSD (float; correlation min eigenvalue)", worst_eig, -1e-4);
    check(wc.f->raw().mekf().couplingConfidence() > 0.9f && wq.f->raw().mekf().couplingConfidence() > 0.9f,
          "the estimated coupling engages (c > 0.9)");
    check(ec < 1.05 * e33 && eq < 1.05 * e33, "position RMS within 5 % of the 33-state filter", std::max(ec, eq) / e33, 1.05);
}

int main() {
    test1_structure();
    test2_realization();
    test3_fixed_kappa();
    test4_cascade_a41();
    test5_commits();
    test6_wrapper();
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures ? 1 : 0;
}
