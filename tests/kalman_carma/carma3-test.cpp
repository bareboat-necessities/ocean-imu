/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Checks of the experimental three-axis filter Kalman3D_Wave_CARMA3 and its
  wrapper (tests/kalman_carma/carma3/):

    * state layouts of every benchmarked configuration;
    * the CARMA accelerations enter the attitude-aware accelerometer model
      (vertical acceleration / accelerometer-bias and attitude / horizontal
      acceleration cross-covariances), and the position and velocity
      pseudo-measurement APIs act on the CARMA p and v of every axis;
    * non-uniform steps keep the covariance symmetric positive semidefinite,
      with and without the S_x, S_y and S_z augmentations, and every S state
      stays bounded by its zero pseudo-measurement;
    * horizontal Kronecker discretization: blockdiag(F1, F1) and
      Sigma_h (x) Q1 equal the exact discretization of the joint two-axis
      system with correlated driving noise (double, Van Loan-free doubling
      reference at a finer step), with and without the S augmentation;
    * the stationary horizontal covariance Sigma_h (x) P_inf is invariant
      under the float prediction;
    * auxiliary S: under the prior, Var S grows at S_p(0) per second when the
      shape has a DC level (CARMA(5,2)) and stays bounded when it has a zero
      at DC (CARMA(4,1));
    * float discretization of the horizontal blocks (with S) across the
      tuner's frequency range against double;
    * on a jittered sample clock, the Taylor extension of F and Q from a
      reference step (StepDiscretization) against double, and its refresh
      conditions;
    * a noiseless three-axis sinusoid is recovered in phase and amplitude on
      every axis and the covariance stays positive semidefinite;
    * the correlated (coupled) horizontal prior has the requested correlation;
    * the wrapper on a still, noisy IMU stays finite and bounded on all axes.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdio>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "util/WaveFilesSupport.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"

const float g_std = 9.80665f;

namespace kc = ocean_imu::kalman::carma;
using ocean_imu::kalman::Kalman3D_Wave_CARMA3;
using Vec3 = Eigen::Vector3f;

namespace {

int g_failures = 0, g_checks = 0;
void check(bool ok, const char* what, double v = NAN, double lim = NAN) {
    ++g_checks;
    if (!ok) { ++g_failures; std::printf("FAIL: %s (value=%.6g limit=%.6g)\n", what, v, lim); }
}

template<class F>
double min_eig_normalized(const F& f) {
    const auto P = f.covariance_full().template cast<double>().eval();
    Eigen::VectorXd d = P.diagonal().cwiseMax(1e-300).cwiseSqrt();
    const Eigen::MatrixXd C = d.cwiseInverse().asDiagonal() * P * d.cwiseInverse().asDiagonal();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
    return es.eigenvalues().minCoeff();
}

void test_layout() {
    check(Kalman3D_Wave_CARMA3<float, 4, false, 4, false>::NX == 21, "CARMA(4,1) x3 has 21 states");
    check(Kalman3D_Wave_CARMA3<float, 5, false, 5, false>::NX == 24, "CARMA(5,2) x3 has 24 states");
    check(Kalman3D_Wave_CARMA3<float, 4, false, 5, false>::NX == 22, "CARMA(4,1) xy + (5,2) z has 22 states");
    check(Kalman3D_Wave_CARMA3<float, 5, true, 5, false>::NX == 26, "CARMA(5,2) x3 + S_x,S_y has 26 states");
    check(Kalman3D_Wave_CARMA3<float, 4, true, 4, false>::NX == 23, "CARMA(4,1) x3 + S_x,S_y has 23 states");
    check(Kalman3D_Wave_CARMA3<float, 3, false, 3, false>::NX == 18, "CARMA(3,0) x3 has 18 states");
    using F = Kalman3D_Wave_CARMA3<float, 5, true, 4, true>;
    check(F::idx_p(0) == 6 && F::idx_v(0) == 7 && F::idx_a(0) == 8 && F::idx_S(0) == 11, "x block layout");
    check(F::idx_p(1) == 12 && F::idx_S(1) == 17 && F::OFF_Z == 18 && F::OFF_SZ == 22 && F::OFF_BA == 23,
          "y, z, S_z, b_a offsets");
}

template<class F>
F make_filter() {
    const Vec3 sa = Vec3::Constant(0.03f), sg = Vec3::Constant(3e-5f), sm = Vec3::Constant(8.0f);
    F f(sa, sg, sm);
    f.set_horizontal_carma(0.31f / 0.37f, 0.5f, 0.45f, 0.0f);
    f.set_vertical_carma(0.31f / 0.37f, 0.6f);
    f.reset_aw_covariance_to_stationary();
    f.set_RS_noise(Vec3(0.7f, 0.7f, 0.7f));
    f.set_pseudo_update_period_s(0.03f);
    f.initialize_from_attitude(Eigen::Quaternionf::Identity(), 0.01f, 0.05f);
    return f;
}

void test_cross_coupling_and_api() {
    using F = Kalman3D_Wave_CARMA3<float, 4, false, 4, false>;
    F f = make_filter<F>();
    // Tilted attitude: the vertical acceleration enters every body axis.
    const Eigen::Quaternionf q(Eigen::AngleAxisf(0.3f, Vec3::UnitX()));
    f.initialize_from_attitude(q, 0.02f, 0.05f);
    for (int k = 0; k < 400; ++k) {
        f.time_update(Vec3::Zero(), 0.005f);
        const Vec3 fb = q.conjugate() * Vec3(0.0f, 0.0f, -9.80665f);
        f.measurement_update_acc_only(fb);
    }
    // a_z and the accelerometer bias both act along gravity, so the update
    // must correlate them; the horizontal accelerations act orthogonally to
    // gravity, so they couple to attitude through -[f_b]x.
    const auto& P = f.covariance_full();
    double c_ab = 0.0;
    for (int k = 0; k < 3; ++k) c_ab = std::max(c_ab, std::abs(double(P(F::OFF_Z + 2, F::OFF_BA + k))));
    const double norm = std::sqrt(double(P(F::OFF_Z + 2, F::OFF_Z + 2)) * double(P(F::OFF_BA + 2, F::OFF_BA + 2)));
    check(c_ab > 1e-3 * norm, "vertical acceleration / accelerometer-bias cross-covariance retained",
          c_ab / norm, 1e-3);
    double c_ah = 0.0;
    for (int k = 0; k < 3; ++k) c_ah = std::max(c_ah, std::abs(double(P(k, F::idx_a(1)))));
    check(c_ah > 0.0, "attitude / horizontal acceleration cross-covariance retained", c_ah, 0.0);

    // Pseudo-measurement APIs act on the CARMA displacement and velocity of
    // every axis.
    const Vec3 p_meas(0.11f, -0.07f, 0.37f), v_meas(-0.02f, 0.03f, 0.05f);
    f.measurement_update_position_pseudo(p_meas, Vec3::Constant(1e-4f));
    check((f.get_position() - p_meas).cwiseAbs().maxCoeff() < 1e-3f, "position pseudo-measurement maps to p_x, p_y, h",
          (f.get_position() - p_meas).cwiseAbs().maxCoeff(), 1e-3);
    f.measurement_update_velocity_pseudo(v_meas, Vec3::Constant(1e-4f));
    check((f.get_velocity() - v_meas).cwiseAbs().maxCoeff() < 1e-3f, "velocity pseudo-measurement maps to v_x, v_y, v_z",
          (f.get_velocity() - v_meas).cwiseAbs().maxCoeff(), 1e-3);
    // A second, equally precise measurement averages with the first.
    const float v_prior = f.get_velocity().z();
    const float var_prior = f.covariance_full()(F::OFF_Z + 1, F::OFF_Z + 1);
    const float r_v = 1e-3f;
    f.measurement_update_vert_velocity_pseudo(-0.21f, r_v);
    const float expect = v_prior + var_prior / (var_prior + r_v * r_v) * (-0.21f - v_prior);
    check(std::abs(f.get_velocity().z() - expect) < 1e-4f, "vertical velocity pseudo-measurement is the Kalman update",
          f.get_velocity().z(), expect);
    check(min_eig_normalized(f) > -1e-5, "covariance PSD after pseudo-measurements", min_eig_normalized(f), -1e-5);
}

template<class F>
void test_nonuniform_and_integral(const char* name) {
    F f = make_filter<F>();
    std::mt19937 rng(7u);
    std::uniform_real_distribution<float> jit(0.6f, 1.6f);
    std::normal_distribution<float> n01(0.0f, 1.0f);
    double t = 0.0, max_p = 0.0;
    bool finite = true;
    for (int k = 0; k < 200 * 300; ++k) {
        const float dt = 0.005f * jit(rng);
        t += dt;
        const double w = 2.0 * M_PI * 0.27 * t;
        const Vec3 aw(static_cast<float>(0.4 * std::sin(w + 0.4)), static_cast<float>(0.3 * std::cos(w)),
                      static_cast<float>(1.2 * std::sin(w)));
        f.time_update(Vec3(1e-3f * n01(rng), 1e-3f * n01(rng), 1e-3f * n01(rng)), dt);
        f.measurement_update_acc_only(Vec3(aw.x() + 0.02f * n01(rng), aw.y() + 0.02f * n01(rng),
                                           aw.z() - 9.80665f + 0.02f * n01(rng)));
        if (k % 8 == 0) f.measurement_update_mag_only(Vec3(20.0f, -4.0f, 45.0f));
        finite = finite && f.get_position().allFinite() && f.covariance_full().allFinite();
        max_p = std::max(max_p, double(f.get_position().cwiseAbs().maxCoeff()));
    }
    char buf[128];
    std::snprintf(buf, sizeof buf, "%s, non-uniform steps: finite state and covariance", name);
    check(finite, buf);
    std::snprintf(buf, sizeof buf, "%s, non-uniform steps: bounded displacement", name);
    check(max_p < 2.0, buf, max_p, 2.0);
    const double me = min_eig_normalized(f);
    std::snprintf(buf, sizeof buf, "%s, non-uniform steps: covariance PSD", name);
    check(me > -1e-5, buf, me, -1e-5);
    const auto P = f.covariance_full();
    std::snprintf(buf, sizeof buf, "%s: covariance exactly symmetric", name);
    check((P - P.transpose()).cwiseAbs().maxCoeff() == 0.0f, buf);
    const Vec3 S = f.get_integral_displacement();
    double max_s = 0.0;
    for (int i = 0; i < 3; ++i) {
        const bool has = i < 2 ? F::idx_S(i) >= 0 : F::OFF_SZ >= 0;
        std::snprintf(buf, sizeof buf, "%s: S axis %d %s", name, i, has ? "bounded by its pseudo-measurement" : "absent");
        if (has) {
            check(std::isfinite(S(i)) && std::abs(S(i)) < 5.0f, buf, S(i), 5.0);
            max_s = std::max(max_s, double(std::abs(S(i))));
        } else {
            check(std::isnan(S(i)), buf);
        }
    }
    std::printf("[%s] non-uniform steps: max|p| %.3f m, max|S| %.3f m s, min corr eig %.2e\n", name, max_p, max_s, me);
}

// Exact discretization of the joint two-axis system by the library routine at
// h, against doublings of a 1/64 step (both double).
template<int N, bool SH>
void test_kronecker(const kc::CarmaShape& sh, const char* label) {
    constexpr int NA = N + (SH ? 1 : 0);
    constexpr int N2 = 2 * NA;
    kc::CarmaVerticalModel<double, N> m(sh);
    m.configure(0.9, 1.0);
    Eigen::Matrix<double,NA,NA> A1 = Eigen::Matrix<double,NA,NA>::Zero(), W1 = A1, F1, Q1;
    A1.template topLeftCorner<N,N>() = m.A();
    W1.template topLeftCorner<N,N>() = m.W();
    if constexpr (SH) A1(N, 0) = 1.0;
    const double h = 0.005;
    kc::discretize_lti<double,NA>(A1, W1, h, F1, Q1);
    Eigen::Matrix2d Sg; Sg << 0.49, 0.6 * 0.7 * 0.5, 0.6 * 0.7 * 0.5, 0.25;
    Eigen::Matrix<double,N2,N2> A2 = Eigen::Matrix<double,N2,N2>::Zero(), W2 = A2, F2, Q2;
    A2.template topLeftCorner<NA,NA>() = A1; A2.template bottomRightCorner<NA,NA>() = A1;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) W2.template block<NA,NA>(a * NA, b * NA) = Sg(a, b) * W1;
    // Reference: 6 doublings of the 1/64 step.
    Eigen::Matrix<double,N2,N2> Fs, Qs;
    kc::discretize_lti<double,N2>(A2, W2, h / 64.0, Fs, Qs);
    for (int i = 0; i < 6; ++i) { Qs = (Qs + Fs * Qs * Fs.transpose()).eval(); Fs = (Fs * Fs).eval(); }
    kc::discretize_lti<double,N2>(A2, W2, h, F2, Q2);
    Eigen::Matrix<double,N2,N2> Fk = Eigen::Matrix<double,N2,N2>::Zero(), Qk = Fk;
    Fk.template topLeftCorner<NA,NA>() = F1; Fk.template bottomRightCorner<NA,NA>() = F1;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) Qk.template block<NA,NA>(a * NA, b * NA) = Sg(a, b) * Q1;
    // Entry-wise relative to each entry's natural scale sqrt(Q_ii Q_jj).
    double eF = 0, eQ = 0, eR = 0;
    for (int i = 0; i < N2; ++i)
        for (int j = 0; j < N2; ++j) {
            eF = std::max(eF, std::abs(Fk(i,j) - F2(i,j)));
            const double sc = std::sqrt(Qs(i,i) * Qs(j,j));
            eQ = std::max(eQ, std::abs(Qk(i,j) - Q2(i,j)) / sc);
            eR = std::max(eR, std::abs(Qk(i,j) - Qs(i,j)) / sc);
        }
    char buf[160];
    std::snprintf(buf, sizeof buf, "%s: Kronecker F, Q equal the joint discretization", label);
    check(eF < 1e-12 && eQ < 1e-10, buf, std::max(eF, eQ), 1e-10);
    std::snprintf(buf, sizeof buf, "%s: Kronecker Q equals the fine-step doubling reference", label);
    check(eR < 1e-9, buf, eR, 1e-9);
    std::printf("[kronecker] %s: |dF| %.2e, |dQ|/scale %.2e (joint), %.2e (fine-step reference)\n", label, eF, eQ, eR);
}

template<int NHC, bool SH>
void test_stationary_prediction(const char* label) {
    using F = Kalman3D_Wave_CARMA3<float, NHC, SH, 4, false>;
    const Vec3 sa = Vec3::Constant(0.03f), sg = Vec3::Constant(3e-5f), sm = Vec3::Constant(8.0f);
    F f(sa, sg, sm);
    f.set_horizontal_carma(0.85f, 0.45f, 0.3f, 0.5f);
    f.set_vertical_carma(0.85f, 0.6f);
    f.reset_aw_covariance_to_stationary();
    f.set_pseudo_update_period_s(1e6f);
    const auto P0 = f.covariance_full();
    for (int k = 0; k < 20000; ++k) f.time_update(Vec3::Zero(), 0.005f);
    const auto& P = f.covariance_full();
    double e = 0;
    for (int i = 0; i < NHC; ++i)
        for (int j = 0; j < NHC; ++j)
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b) {
                    const int ii = F::idx_p(a) + i, jj = F::idx_p(b) + j;
                    e = std::max(e, std::abs(double(P(ii,jj)) - double(P0(ii,jj))) / std::sqrt(double(P0(ii,ii)) * P0(jj,jj)));
                }
    char buf[128];
    std::snprintf(buf, sizeof buf, "%s: Sigma_h (x) P_inf invariant over 100 s of float prediction", label);
    check(e < 5e-3, buf, e, 5e-3);
    const double corr = P0(F::idx_p(0), F::idx_p(1)) / std::sqrt(double(P0(F::idx_p(0), F::idx_p(0))) * P0(F::idx_p(1), F::idx_p(1)));
    std::snprintf(buf, sizeof buf, "%s: coupled prior has the requested x-y correlation", label);
    check(std::abs(corr - 0.5) < 1e-5, buf, corr, 0.5);
    std::printf("[stationary] %s: max rel drift %.2e, x-y displacement correlation %.6f\n", label, e, corr);
}

// Prior growth of Var S_x with no S update.
template<int NHC>
void test_s_growth(const kc::CarmaShape& sh, bool dc_zero, const char* label) {
    using F = Kalman3D_Wave_CARMA3<float, NHC, true, 4, false>;
    const Vec3 sa = Vec3::Constant(0.03f), sg = Vec3::Constant(3e-5f), sm = Vec3::Constant(8.0f);
    F f(sa, sg, sm, 5e-4f, 1e-6f, 1e-11f, 1.5f, 9.80665f, kc::CarmaDefaultShape<4>::value, sh);
    const float sx = 0.5f;
    f.set_horizontal_carma(0.85f, sx, sx, 0.0f);
    f.reset_aw_covariance_to_stationary();
    // P_SS is seeded at 50^2 m^2 s^2; in float the per-step random-walk
    // increment (~1e-7) is far below its ulp, so collapse it first with one
    // tight S = 0 update and watch the growth from there.
    f.set_RS_noise(Vec3::Constant(0.01f));
    f.applyIntegralZeroPseudoMeas();
    f.set_pseudo_update_period_s(1e6f);
    const int iS = F::idx_S(0);
    const double v0 = f.covariance_full()(iS, iS);
    std::vector<double> vs;
    for (int k = 1; k <= 120000; ++k) {
        f.time_update(Vec3::Zero(), 0.005f);
        if (k % 20000 == 0) vs.push_back(double(f.covariance_full()(iS, iS)) - v0);
    }
    const auto& hm = f.horizontal_model();
    const double c0 = double(hm.noise_gain2()) * sx * sx * kc::carma_unit_psd<NHC>(hm.shape(), hm.kappa(), 0.0);
    const double slope = (vs.back() - vs[1]) / 400.0;   // 200 s .. 600 s
    char buf[128];
    if (!dc_zero) {
        std::snprintf(buf, sizeof buf, "%s: Var S grows at S_p(0) per second (random walk)", label);
        check(std::abs(slope / c0 - 1.0) < 0.03, buf, slope / c0, 1.0);
    } else {
        std::snprintf(buf, sizeof buf, "%s: Var S bounded (DC zero, uncontrollable S mode)", label);
        check(std::abs(slope) < 1e-3 * (std::abs(vs.back()) + 1e-6) + 1e-7, buf, slope, 0.0);
    }
    std::printf("[S] %s: Var S - Var S0 at 100..600 s = %.4g %.4g %.4g %.4g %.4g %.4g; slope %.3e, S_p(0) %.3e\n",
                label, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5], slope, c0);
}

template<int N, bool SH>
void test_float_discretization(const kc::CarmaShape& sh, const char* label) {
    constexpr int NA = N + (SH ? 1 : 0);
    double worst = 0.0;
    for (double fz : {0.05, 0.1, 0.3, 0.6, 1.0}) {
        kc::CarmaVerticalModel<double, N> md(sh);
        md.configure(fz / md.base_zero_crossing_hz(), 1.0);
        kc::CarmaVerticalModel<float, N> mf(sh);
        mf.configure(static_cast<float>(fz / md.base_zero_crossing_hz()), 1.0f);
        Eigen::Matrix<double,NA,NA> Ad = Eigen::Matrix<double,NA,NA>::Zero(), Wd = Ad, Fd, Qd;
        Eigen::Matrix<float,NA,NA> Af = Eigen::Matrix<float,NA,NA>::Zero(), Wf = Af, Ff, Qf;
        Ad.template topLeftCorner<N,N>() = md.A(); Wd.template topLeftCorner<N,N>() = md.W();
        Af.template topLeftCorner<N,N>() = mf.A(); Wf.template topLeftCorner<N,N>() = mf.W();
        if constexpr (SH) { Ad(N, 0) = 1.0; Af(N, 0) = 1.0f; }
        kc::discretize_lti<double,NA>(Ad, Wd, 0.005, Fd, Qd);
        kc::discretize_lti<float,NA>(Af, Wf, 0.005f, Ff, Qf);
        for (int i = 0; i < NA; ++i)
            for (int j = 0; j < NA; ++j) {
                worst = std::max(worst, std::abs(double(Qf(i,j)) - Qd(i,j)) / std::sqrt(Qd(i,i) * Qd(j,j)));
                worst = std::max(worst, std::abs(double(Ff(i,j)) - Fd(i,j)) /
                                        std::max(1.0, std::abs(Fd(i,j))) * 0.1);
            }
    }
    char buf[128];
    std::snprintf(buf, sizeof buf, "%s: float F, Q within 1e-4 of double for f_z 0.05..1 Hz", label);
    check(worst < 1e-4, buf, worst, 1e-4);
    std::printf("[float] %s: worst scaled float error %.2e over f_z 0.05..1 Hz\n", label, worst);
}

// Jittered sample clock (+-300 us on 5 ms, inside the h0/16 neighbourhood)
// over the tuner's frequency range: StepDiscretization against the double
// discretization at the actual step, with the metric and bar of
// test_float_discretization; the plain float discretization at the same
// steps is printed for comparison.  Also: a step outside the neighbourhood
// and a model change need refresh(), and a repeated step leaves F and Q
// unchanged.
template<int N, bool SH>
void test_step_discretization(const kc::CarmaShape& sh, const char* label) {
    constexpr int NA = N + (SH ? 1 : 0);
    using MF = Eigen::Matrix<float,NA,NA>;
    using MD = Eigen::Matrix<double,NA,NA>;
    std::mt19937 rng(11u);
    std::uniform_int_distribution<int> jit_us(-300, 300);
    double worst = 0.0, worst_plain = 0.0;
    int n_series = 0, n_steps = 0;
    bool sym = true, guards = true;
    for (double fz : {0.05, 0.1, 0.3, 0.6, 1.0, 1.2}) {
        kc::CarmaVerticalModel<float, N> mf(sh);
        mf.configure(static_cast<float>(fz / mf.base_zero_crossing_hz()), 1.0f);
        MF Af = MF::Zero(), Wf = Af;
        Af.template topLeftCorner<N,N>() = mf.A(); Wf.template topLeftCorner<N,N>() = mf.W();
        if constexpr (SH) Af(N, 0) = 1.0f;
        const MD Ad = Af.template cast<double>(), Wd = Wf.template cast<double>();
        kc::StepDiscretization<float, NA> sd;
        guards = guards && !sd.advance(mf.generation(), 0.005f);
        sd.refresh(Af, Wf, mf.generation(), 0.005f);
        for (int k = 0; k < 200; ++k) {
            const float h = 0.005f + 1e-6f * static_cast<float>(jit_us(rng));
            ++n_steps;
            if (sd.advance(mf.generation(), h)) ++n_series;
            else sd.refresh(Af, Wf, mf.generation(), h);
            MD Fd, Qd; kc::discretize_lti<double,NA>(Ad, Wd, double(h), Fd, Qd);
            MF Fp, Qp; kc::discretize_lti<float,NA>(Af, Wf, h, Fp, Qp);
            const MF& F = sd.Phi();
            const MF& Q = sd.Q();
            for (int i = 0; i < NA; ++i)
                for (int j = 0; j < NA; ++j) {
                    const double sc = std::sqrt(Qd(i,i) * Qd(j,j));
                    const double fs = std::max(1.0, std::abs(Fd(i,j))) * 10.0;
                    worst = std::max({worst, std::abs(double(Q(i,j)) - Qd(i,j)) / sc,
                                      std::abs(double(F(i,j)) - Fd(i,j)) / fs});
                    worst_plain = std::max({worst_plain, std::abs(double(Qp(i,j)) - Qd(i,j)) / sc,
                                            std::abs(double(Fp(i,j)) - Fd(i,j)) / fs});
                    sym = sym && Q(i,j) == Q(j,i);
                }
        }
        // From a 5 ms reference; a half-microsecond offset is a new step.
        sd.refresh(Af, Wf, mf.generation(), 0.005f);
        const MF F_last = sd.Phi();
        guards = guards && sd.advance(mf.generation(), 0.0050025f) && sd.Phi() != F_last;
        const MF F_rep = sd.Phi(), Q_rep = sd.Q();
        guards = guards && sd.advance(mf.generation(), 0.0050025f) && sd.Phi() == F_rep && sd.Q() == Q_rep;
        guards = guards && !sd.advance(mf.generation(), 0.0060f);
        guards = guards && !sd.advance(mf.generation() + 1u, 0.005f);
    }
    char buf[160];
    std::snprintf(buf, sizeof buf, "%s: every step within +-300 us uses the series", label);
    check(n_series == n_steps, buf, n_series, n_steps);
    std::snprintf(buf, sizeof buf, "%s: series F, Q within 1e-4 of double for f_z 0.05..1.2 Hz, +-300 us jitter", label);
    check(worst < 1e-4, buf, worst, 1e-4);
    std::snprintf(buf, sizeof buf, "%s: series Q exactly symmetric", label);
    check(sym, buf);
    std::snprintf(buf, sizeof buf, "%s: out-of-range step and model change need refresh, repeated step reuses", label);
    check(guards, buf);
    std::printf("[step] %s: %d of %d jittered steps from the series; worst scaled error %.2e, "
                "float discretization at the same steps %.2e\n", label, n_series, n_steps, worst, worst_plain);
}

// Level attitude, noiseless specific force of a three-axis sinusoid.
template<int NHC, int NZC>
void test_sinusoid(const char* label) {
    using F = Kalman3D_Wave_CARMA3<float, NHC, false, NZC, false>;
    const Vec3 sa = Vec3::Constant(0.03f), sg = Vec3::Constant(3e-5f), sm = Vec3::Constant(8.0f);
    F f(sa, sg, sm);
    const double w = 2.0 * M_PI * 0.30, dt = 0.005;
    const double fz_h = f.horizontal_model().base_zero_crossing_hz();
    const double fz_v = f.vertical_model().base_zero_crossing_hz();
    f.set_horizontal_carma(static_cast<float>(0.30 / fz_h), 0.4f, 0.3f, 0.0f);
    f.set_vertical_carma(static_cast<float>(0.30 / fz_v), 0.6f);
    f.reset_aw_covariance_to_stationary();
    f.set_acc_bias_updates_enabled(false);
    f.initialize_from_attitude(Eigen::Quaternionf::Identity(), 0.01f, 0.05f);
    // Magnetometer at 25 Hz: without it the horizontal motion prior makes yaw
    // weakly (and wrongly) observable through the accelerometer model.
    const Vec3 mag_w(20.0f, -4.0f, 45.0f);
    f.set_mag_world_ref(mag_w);
    const double A[3] = {0.12, 0.08, 0.25}, ph[3] = {0.4, 1.3, 0.0};
    double se[3] = {}, sr[3] = {};
    long n = 0;
    for (int k = 0; k < 24000; ++k) {
        const double t = k * dt;
        Vec3 aw;
        for (int i = 0; i < 3; ++i) aw(i) = static_cast<float>(-A[i] * w * w * std::sin(w * t + ph[i]));
        f.time_update(Vec3::Zero(), static_cast<float>(dt));
        f.measurement_update_acc_only(Vec3(aw.x(), aw.y(), aw.z() - 9.80665f));
        if (k % 8 == 0) f.measurement_update_mag_only(mag_w);
        if (t > 60.0) {
            const Vec3 p = f.get_position();
            for (int i = 0; i < 3; ++i) {
                const double tr = A[i] * std::sin(w * t + ph[i]);
                se[i] += std::pow(p(i) - tr, 2); sr[i] += tr * tr;
            }
            ++n;
        }
    }
    double worst = 0;
    for (int i = 0; i < 3; ++i) worst = std::max(worst, std::sqrt(se[i] / sr[i]));
    char buf[128];
    std::snprintf(buf, sizeof buf, "%s: three-axis sinusoid recovered within 10 %% RMS on every axis", label);
    check(worst < 0.10, buf, worst, 0.10);
    const double me = min_eig_normalized(f);
    std::snprintf(buf, sizeof buf, "%s: covariance PSD after the sinusoid", label);
    check(me > -1e-4, buf, me, -1e-4);
    std::printf("[sinusoid] %s: relative RMS error x %.3f y %.3f z %.3f, min corr eig %.2e\n", label,
                std::sqrt(se[0] / sr[0]), std::sqrt(se[1] / sr[1]), std::sqrt(se[2] / sr[2]), me);
}

template<int NHC, bool SH, int NZC>
void test_wrapper_still(const char* label) {
    using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, NHC, SH, NZC, false>;
    Fusion fusion;
    typename Fusion::Config cfg;
    cfg.with_mag = true;
    cfg.sigma_a = Vec3::Constant(2.8f * 1.51e-3f * 9.80665f * 0.71f);
    cfg.sigma_g = Vec3::Constant(2.0f * 0.00157f * 0.01f);
    cfg.sigma_m = Vec3::Constant(1.2f * 0.8f * 8.0f);
    fusion.begin(cfg);
    fusion.raw().enableTuner(true);
    fusion.raw().enableClamp(true);
    std::mt19937 rng(3);
    std::normal_distribution<float> n(0.0f, 1.0f);
    const Vec3 mag_w(20.0f, -4.0f, 45.0f);
    double peak = 0;
    bool finite = true;
    for (int k = 0; k < 200 * 400; ++k) {
        const Vec3 acc(0.004f + 0.015f * n(rng), -0.003f + 0.015f * n(rng), -9.80665f + 0.006f + 0.015f * n(rng));
        const Vec3 gyr(3e-4f * n(rng), 3e-4f * n(rng), 3e-4f * n(rng));
        fusion.update(0.005f, gyr, acc, 35.0f);
        if (k % 8 == 0) fusion.updateMag(mag_w + Vec3(0.8f * n(rng), 0.8f * n(rng), 0.8f * n(rng)));
        if (fusion.isLive()) {
            const Vec3 p = fusion.raw().mekf().get_position();
            finite = finite && p.allFinite();
            if (k > 200 * 240) peak = std::max(peak, double(p.cwiseAbs().maxCoeff()));
        }
    }
    char buf[128];
    std::snprintf(buf, sizeof buf, "%s wrapper: still IMU stays finite", label);
    check(finite && fusion.isLive(), buf);
    std::snprintf(buf, sizeof buf, "%s wrapper: still IMU |p| < 0.5 m after 240 s on every axis", label);
    check(peak < 0.5, buf, peak, 0.5);
    std::printf("[wrapper] %s: still IMU, peak |p| after 240 s %.4f m\n", label, peak);
}

}  // namespace

int main() {
    test_layout();
    test_cross_coupling_and_api();
    test_nonuniform_and_integral<Kalman3D_Wave_CARMA3<float, 4, false, 4, false>>("CARMA(4,1) x3");
    test_nonuniform_and_integral<Kalman3D_Wave_CARMA3<float, 5, true, 5, true>>("CARMA(5,2) x3 + S_x,S_y,S_z");
    test_kronecker<4, false>(kc::CARMA41_JONSWAP_SHAPE, "CARMA(4,1)");
    test_kronecker<5, false>(kc::CARMA52_JONSWAP_SHAPE, "CARMA(5,2)");
    test_kronecker<5, true>(kc::CARMA52_JONSWAP_SHAPE, "CARMA(5,2)+S");
    test_stationary_prediction<4, false>("CARMA(4,1) x3");
    test_stationary_prediction<5, false>("CARMA(5,2) x3");
    test_s_growth<5>(kc::CARMA52_JONSWAP_SHAPE, false, "CARMA(5,2)+S");
    test_s_growth<4>(kc::CARMA41_JONSWAP_SHAPE, true, "CARMA(4,1)+S");
    test_float_discretization<4, false>(kc::CARMA41_JONSWAP_SHAPE, "CARMA(4,1)");
    test_float_discretization<5, false>(kc::CARMA52_JONSWAP_SHAPE, "CARMA(5,2)");
    test_float_discretization<5, true>(kc::CARMA52_JONSWAP_SHAPE, "CARMA(5,2)+S");
    test_float_discretization<4, true>(kc::CARMA41_JONSWAP_SHAPE, "CARMA(4,1)+S");
    test_step_discretization<3, false>(kc::CARMA30_JONSWAP_SHAPE, "CARMA(3,0)");
    test_step_discretization<4, false>(kc::CARMA41_JONSWAP_SHAPE, "CARMA(4,1)");
    test_step_discretization<5, false>(kc::CARMA52_JONSWAP_SHAPE, "CARMA(5,2)");
    test_step_discretization<4, true>(kc::CARMA41_JONSWAP_SHAPE, "CARMA(4,1)+S");
    test_step_discretization<5, true>(kc::CARMA52_JONSWAP_SHAPE, "CARMA(5,2)+S");
    test_sinusoid<4, 4>("CARMA(4,1) x3");
    test_sinusoid<5, 5>("CARMA(5,2) x3");
    test_sinusoid<4, 5>("CARMA(4,1) xy + (5,2) z");
    test_wrapper_still<4, false, 4>("CARMA(4,1) x3");
    test_wrapper_still<5, true, 5>("CARMA(5,2) x3 + S");
    std::printf("carma3-test: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
