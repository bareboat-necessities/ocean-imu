/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Mathematical checks of the EXPERIMENTAL two-frame Lie-group CARMA filter
  Kalman3D_Wave_TFG_CARMA3 (docs/tfg-carma3-research.md).

  No expected value below is "what the implementation already did".  The
  references are built independently of the filter's formulas:

    * finite differences of the nonlinear group maps (error flow over a step,
      residuals, Log_G(Exp_G(-a) Exp_G(a + d)));
    * fine-step RK4 integration of the stated continuous error model
      (transition and Lyapunov process covariance);
    * the structured Gauss-Legendre impulse-response quadrature of
      Kalman3D_Wave_TFG, applied to the CARMA prior;
    * the reset Jacobian of the independent implementation Kalman3D_Wave_TFG
      (the same group, TwoFrameGroup<double, 4, 2>);
    * Monte Carlo simulation of the nonlinear stochastic model;
    * Kalman3D_Wave_CARMA3 (quaternion MEKF) run side by side.

  Every check prints the measured discrepancy and its limit
  ("CHECK <name> value=<v> limit=<l> PASS|FAIL"); a limit of inf marks a
  reported quantity that is not a pass/fail criterion.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma/Kalman3D_Wave_TFG_CARMA3.h"
#include "kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h"
#include "kalman_tfg/Kalman3D_Wave_TFG.h"

const float g_std = 9.80665f;

namespace kc = ocean_imu::kalman::carma;
namespace lie = ocean_imu::lie;
using ocean_imu::kalman::Kalman3D_Wave_CARMA3;
using ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3;

namespace {

int g_failures = 0, g_checks = 0;
constexpr double kInf = std::numeric_limits<double>::infinity();

void check(const std::string& name, double value, double limit) {
    ++g_checks;
    const bool ok = std::isfinite(value) ? (value <= limit) : false;
    if (!ok) ++g_failures;
    std::printf("CHECK %-62s value=%.3e limit=%.1e %s\n", name.c_str(), value, limit, ok ? "PASS" : "FAIL");
}
void report(const std::string& name, double value) {
    std::printf("CHECK %-62s value=%.3e limit=inf REPORT\n", name.c_str(), value);
}

using F = Kalman3D_Wave_TFG_CARMA3<double, 4>;
using G = F::Group;
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
constexpr int NX = F::NX;
using MX = Eigen::Matrix<double, NX, NX>;
using TX = Eigen::Matrix<double, NX, 1>;
constexpr int NC = 4;
using MC = Eigen::Matrix<double, NC, NC>;
constexpr double GRAV = 9.80665;

std::mt19937 rng(20261028u);

V3 rvec(double s) {
    std::normal_distribution<double> n(0.0, 1.0);
    return V3(n(rng), n(rng), n(rng)) * s;
}
M3 rot_euler(double roll, double pitch, double yaw) {
    return (Eigen::AngleAxisd(yaw, V3::UnitZ()) * Eigen::AngleAxisd(pitch, V3::UnitY()) *
            Eigen::AngleAxisd(roll, V3::UnitX())).toRotationMatrix();
}
template<int N>
Eigen::Matrix<double, N, N> rand_spd(double scale, double floor_frac = 0.05) {
    std::normal_distribution<double> n(0.0, 1.0);
    Eigen::Matrix<double, N, N> A;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) A(i, j) = n(rng);
    Eigen::Matrix<double, N, N> P = A * A.transpose() / N;
    P += Eigen::Matrix<double, N, N>::Identity() * floor_frac;
    return scale * P;
}
// Physically scaled covariance: per-block standard deviations and random
// correlations (attitude-translation included).
MX scaled_cov(double s_att, double s_bg, double s_w, double s_ba) {
    MX C = rand_spd<NX>(1.0, 0.3);
    Eigen::Matrix<double, NX, 1> d = C.diagonal().cwiseSqrt();
    C = d.cwiseInverse().asDiagonal() * C * d.cwiseInverse().asDiagonal();
    Eigen::Matrix<double, NX, 1> s;
    for (int i = 0; i < NX; ++i) {
        if (i < 3) s(i) = s_att;
        else if (i < 6) s(i) = s_bg;
        else if (i < F::OFF_BA) s(i) = s_w;
        else s(i) = s_ba;
    }
    return s.asDiagonal() * C * s.asDiagonal();
}
double rel_max(const MX& A, const MX& B) {
    // max |A-B| / sqrt(|B_ii B_jj|): the discrepancy in correlation units.
    double m = 0;
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j) {
            const double sc = std::sqrt(std::abs(B(i, i) * B(j, j))) + 1e-300;
            m = std::max(m, std::abs(A(i, j) - B(i, j)) / sc);
        }
    return m;
}
double min_eig_corr(const MX& P) {
    Eigen::Matrix<double, NX, 1> d = P.diagonal().cwiseMax(1e-300).cwiseSqrt();
    const MX C = d.cwiseInverse().asDiagonal() * P * d.cwiseInverse().asDiagonal();
    Eigen::SelfAdjointEigenSolver<MX> es(C);
    return es.eigenvalues().minCoeff();
}

// Default CARMA(4,1) research configuration (f_z = 0.31 Hz).
constexpr double KAPPA = 0.31 / 0.3708;

F make_filter(double sx = 0.5, double sy = 0.5, double sz = 0.6, double r = 0.0) {
    F f(V3::Constant(0.03), V3::Constant(3e-4), V3::Constant(0.8), 5e-4, 1e-6, 1e-8, 1.5, GRAV);
    f.set_vertical_carma(KAPPA, sz);
    f.set_horizontal_carma(KAPPA, sx, sy, r);
    f.reset_aw_covariance_to_stationary();
    f.set_mag_world_ref(V3(18.0, 1.5, 43.0));
    f.set_Q_bacc_rw(V3::Constant(5e-4));
    f.set_acc_bias_time_constant(800.0);
    f.set_Q_bgyro_rw(V3(1e-7, 2e-7, 1.5e-7));
    f.set_acc_bias_updates_enabled(true);
    return f;
}

G make_state(const M3& R, double s_w, const V3& bg, const V3& ba) {
    G X;
    X.R = R;
    for (int i = 0; i < NC; ++i) X.X.col(i) = rvec(s_w * (i == 3 ? 3.0 : 1.0));
    X.B.col(0) = bg;
    X.B.col(1) = ba;
    return X;
}

// ---------------------------------------------------------------------------
// Independent nonlinear models
// ---------------------------------------------------------------------------

struct Flow {
    MC Fc;
    double h;
    double phi_b;
};
Flow make_flow(const F& f, double h) {
    Flow fl;
    MC Q;
    kc::discretize_lti<double, NC>(f.vertical_model().A(), f.vertical_model().W(), h, fl.Fc, Q);
    fl.h = h;
    fl.phi_b = f.acc_bias_updates_enabled() ? std::exp(-h / f.get_acc_bias_time_constant()) : 1.0;
    return fl;
}
// Deterministic propagation of a state under measured rate w_m (constant over
// the step): R <- R Exp((w_m - b_g) h), W <- W Fc^T, b_a <- phi_b b_a.
G flow(const G& X, const V3& wm, const Flow& fl) {
    G Y = X;
    Y.R = X.R * lie::Exp<double>(V3((wm - X.B.col(0)) * fl.h));
    Y.X = X.X * fl.Fc.transpose();
    Y.B.col(1) = X.B.col(1) * fl.phi_b;
    return Y;
}
TX err(const G& Xh, const G& X) { return Xh.local(X); }   // Log_G(Xh X^-1)
G truth_from(const G& Xh, const TX& xi) { return G::Exp_G(TX(-xi)).compose(Xh); }

// Finite-difference Jacobian of the error after one step w.r.t. the error
// before it, at error xi0 (central differences, step eps).
MX fd_transition(const G& Xh, const V3& wm, const Flow& fl, const TX& xi0, double eps) {
    MX J;
    const G Xh1 = flow(Xh, wm, fl);
    for (int j = 0; j < NX; ++j) {
        TX dp = xi0, dm = xi0;
        dp(j) += eps;
        dm(j) -= eps;
        const TX ep = err(Xh1, flow(truth_from(Xh, dp), wm, fl));
        const TX em = err(Xh1, flow(truth_from(Xh, dm), wm, fl));
        J.col(j) = (ep - em) / (2 * eps);
    }
    return J;
}

// Continuous error model: d xi = Fct(t) xi dt + L(t) dB (stated in the filter
// header), integrated by RK4 over one step for Phi and the Lyapunov Q.
struct ContModel {
    G X0;
    V3 w;          // bias-corrected rate (body)
    MC A;          // CARMA dynamics (unit variance)
    double Waa;    // driving intensity of the acceleration state
    M3 Sig;        // spatial covariance
    M3 Qg, Qbg, Qba;
    double tau_b;
    bool ba_active;
};
void cont_F(const ContModel& m, double t, MX& Fm, MX& LQL) {
    const M3 R = m.X0.R * lie::Exp<double>(V3(m.w * t));
    Eigen::Matrix<double, 3, NC> W;
    {
        // W(t) = W0 exp(A t)^T via the library's exact discretization.
        MC Ft, Qt;
        kc::discretize_lti<double, NC>(m.A, MC::Zero(), t, Ft, Qt);
        W = m.X0.X * Ft.transpose();
    }
    const V3 ww = R * m.w;
    Fm.setZero();
    Fm.block<3, 3>(0, 3) = -M3::Identity();
    Fm.block<3, 3>(3, 3) = lie::skew<double>(ww);
    for (int i = 0; i < NC; ++i) {
        for (int j = 0; j < NC; ++j)
            if (m.A(i, j) != 0.0) Fm.block<3, 3>(6 + 3 * i, 6 + 3 * j) = m.A(i, j) * M3::Identity();
        Fm.block<3, 3>(6 + 3 * i, 3) = -lie::skew<double>(V3(W.col(i)));
    }
    Fm.block<3, 3>(F::OFF_BA, F::OFF_BA) = lie::skew<double>(ww) - (m.ba_active ? 1.0 / m.tau_b : 0.0) * M3::Identity();
    // Noise: gyro (phi, rho), gyro bias (beta_g), CARMA (rho_a rows), accel bias.
    Eigen::Matrix<double, NX, 3> Lg = Eigen::Matrix<double, NX, 3>::Zero();
    Lg.block<3, 3>(0, 0) = R;
    for (int i = 0; i < NC; ++i) Lg.block<3, 3>(6 + 3 * i, 0) = lie::skew<double>(V3(W.col(i))) * R;
    LQL = Lg * m.Qg * Lg.transpose();
    LQL.block<3, 3>(3, 3) += R * m.Qbg * R.transpose();
    LQL.block<3, 3>(6 + 3 * 2, 6 + 3 * 2) += m.Waa * m.Sig;
    if (m.ba_active) LQL.block<3, 3>(F::OFF_BA, F::OFF_BA) += R * m.Qba * R.transpose();
}
void rk4_reference(const ContModel& m, double h, int n, MX& Phi, MX& Q) {
    Phi.setIdentity();
    Q.setZero();
    const double dt = h / n;
    auto deriv = [&](double t, const MX& P, const MX& S, MX& dP, MX& dS) {
        MX Fm, LQL;
        cont_F(m, t, Fm, LQL);
        dP = Fm * P;
        dS = Fm * S + S * Fm.transpose() + LQL;
    };
    for (int k = 0; k < n; ++k) {
        const double t = k * dt;
        MX k1p, k1s, k2p, k2s, k3p, k3s, k4p, k4s;
        deriv(t, Phi, Q, k1p, k1s);
        deriv(t + dt / 2, Phi + dt / 2 * k1p, Q + dt / 2 * k1s, k2p, k2s);
        deriv(t + dt / 2, Phi + dt / 2 * k2p, Q + dt / 2 * k2s, k3p, k3s);
        deriv(t + dt, Phi + dt * k3p, Q + dt * k3s, k4p, k4s);
        Phi += dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p);
        Q += dt / 6 * (k1s + 2 * k2s + 2 * k3s + k4s);
    }
}
ContModel cont_model(const F& f, const V3& wm, const M3& Qg, const M3& Qbg, const M3& Qba) {
    ContModel m;
    m.X0 = f.state();
    m.w = wm - f.state().B.col(0);
    m.A = f.vertical_model().A();
    m.Waa = f.vertical_model().W()(2, 2);
    m.Sig = f.spatial_sigma();
    m.Qg = Qg;
    m.Qbg = Qbg;
    m.Qba = Qba;
    m.tau_b = f.get_acc_bias_time_constant();
    m.ba_active = f.acc_bias_updates_enabled();
    return m;
}

// ---------------------------------------------------------------------------
// A. Lie-group identities (in the CARMA column labeling)
// ---------------------------------------------------------------------------

void test_group() {
    double e_inv = 0, e_assoc = 0, e_log = 0, e_rinv = 0, e_local = 0, e_bias = 0;
    for (int trial = 0; trial < 50; ++trial) {
        const double ang = (trial % 5 == 4) ? 2.5 : 0.7;
        const G A = make_state(lie::Exp<double>(rvec(ang)), 0.8, rvec(0.01), rvec(0.05));
        const G B = make_state(lie::Exp<double>(rvec(ang)), 0.8, rvec(0.01), rvec(0.05));
        const G C = make_state(lie::Exp<double>(rvec(ang)), 0.8, rvec(0.01), rvec(0.05));
        const G I = A.compose(A.inverse());
        e_inv = std::max({e_inv, (I.R - M3::Identity()).cwiseAbs().maxCoeff(), I.X.cwiseAbs().maxCoeff(),
                          I.B.cwiseAbs().maxCoeff()});
        const G L = A.compose(B).compose(C), Rr = A.compose(B.compose(C));
        e_assoc = std::max({e_assoc, (L.R - Rr.R).cwiseAbs().maxCoeff(), (L.X - Rr.X).cwiseAbs().maxCoeff(),
                            (L.B - Rr.B).cwiseAbs().maxCoeff()});
        for (double s : {1e-6, 1e-2, 0.5, 2.5}) {
            TX xi;
            for (int i = 0; i < NX; ++i) xi(i) = std::normal_distribution<double>(0, 1)(rng);
            xi.head<3>() = xi.head<3>().normalized() * s;
            e_log = std::max(e_log, (G::Log_G(G::Exp_G(xi)) - xi).cwiseAbs().maxCoeff());
        }
        // Right invariance of the error: (A Z)(B Z)^-1 = A B^-1.
        e_rinv = std::max(e_rinv, (err(A.compose(C), B.compose(C)) - err(A, B)).cwiseAbs().maxCoeff());
        const G R2 = B.retract(A.local(B));
        e_local = std::max({e_local, (R2.R - A.R).cwiseAbs().maxCoeff(), (R2.X - A.X).cwiseAbs().maxCoeff(),
                            (R2.B - A.B).cwiseAbs().maxCoeff()});
        // Body-frame biases: the bias slot of A B^-1 is R_B (b_A - b_B) exactly.
        const G E = A.compose(B.inverse());
        e_bias = std::max(e_bias, (E.B - B.R * (A.B - B.B)).cwiseAbs().maxCoeff());
    }
    check("group: X X^-1 = I", e_inv, 1e-12);
    check("group: associativity", e_assoc, 1e-12);
    check("group: Log_G(Exp_G(xi)) = xi, |phi| 1e-6..2.5", e_log, 1e-11);
    check("group: right invariance of eta = Xhat X^-1", e_rinv, 1e-11);
    check("group: Y.retract(X.local(Y)) = X", e_local, 1e-12);
    check("group: bias slot of Xhat X^-1 is R (bhat - b)", e_bias, 1e-14);

    // Linearized bias error: beta = R (bhat - b) to first order in xi.
    const G Xh = make_state(rot_euler(0.6, -0.4, 2.0), 0.5, rvec(0.01), rvec(0.05));
    TX xi = TX::Zero();
    xi.segment<3>(F::OFF_BG) = rvec(1e-5);
    xi.segment<3>(F::OFF_BA) = rvec(1e-5);
    xi.head<3>() = rvec(1e-5);
    const G X = truth_from(Xh, xi);
    const V3 bg = X.R * (Xh.B.col(0) - X.B.col(0));
    check("group: beta_g = R (bhat_g - b_g) to first order", (bg - xi.segment<3>(F::OFF_BG)).norm() / 1e-5, 1e-4);
}

// ---------------------------------------------------------------------------
// B. Reset Jacobian (finite-correction covariance transport)
// ---------------------------------------------------------------------------

MX fd_reset(const TX& a, double eps) {
    MX J;
    const G Ea = G::Exp_G(TX(-a));
    for (int j = 0; j < NX; ++j) {
        TX dp = a, dm = a;
        dp(j) += eps;
        dm(j) -= eps;
        J.col(j) = (G::Log_G(Ea.compose(G::Exp_G(dp))) - G::Log_G(Ea.compose(G::Exp_G(dm)))) / (2 * eps);
    }
    return J;
}

void test_reset() {
    F f = make_filter();
    ocean_imu::kalman::Kalman3D_Wave_TFG<double> ref;
    for (double s_att : {1e-3, 0.1, 1.0}) {
        double e_fd = 0, e_ref = 0;
        for (int trial = 0; trial < 10; ++trial) {
            TX a;
            a.head<3>() = rvec(s_att);
            for (int i = 3; i < NX; ++i) a(i) = std::normal_distribution<double>(0, s_att * 2)(rng);
            MX J, Jref;
            f.build_reset_jacobian(a, J);
            ref.build_reset_jacobian(a, Jref);
            e_fd = std::max(e_fd, (J - fd_reset(a, 1e-6)).cwiseAbs().maxCoeff());
            e_ref = std::max(e_ref, (J - Jref).cwiseAbs().maxCoeff());
        }
        char n1[96], n2[96];
        std::snprintf(n1, sizeof n1, "reset: J_r(a) vs FD of Log(Exp(-a)Exp(a+d)), |a| ~ %.0e", s_att);
        std::snprintf(n2, sizeof n2, "reset: J_r(a) vs Kalman3D_Wave_TFG, |a| ~ %.0e", s_att);
        check(n1, e_fd, 1e-8);
        check(n2, e_ref, 1e-12);
    }
    // Monte Carlo: ξ+ = Log(Exp(-a) Exp(a + d)), d ~ N(0, P) with correlated
    // attitude/translation, against J P J^T.
    {
        TX a;
        a.head<3>() = V3(0.05, -0.03, 0.2);
        for (int i = 3; i < NX; ++i) a(i) = std::normal_distribution<double>(0, 0.1)(rng);
        const MX P = scaled_cov(0.01, 1e-3, 0.02, 1e-3);
        MX J;
        f.build_reset_jacobian(a, J);
        const MX Pj = J * P * J.transpose();
        Eigen::LLT<MX> llt(P);
        const MX Lc = llt.matrixL();
        const int N = 20000;
        MX S = MX::Zero();
        std::normal_distribution<double> n(0, 1);
        const G Ea = G::Exp_G(TX(-a));
        for (int k = 0; k < N; ++k) {
            TX z;
            for (int i = 0; i < NX; ++i) z(i) = n(rng);
            const TX d = Lc * z;
            const TX xp = G::Log_G(Ea.compose(G::Exp_G(TX(a + d))));
            S += xp * xp.transpose();
        }
        S /= N;

        check("reset: MC cov of transported error vs J P J^T (corr units)", rel_max(S, Pj), 0.04);
        report("reset: MC cov vs untransported P (corr units, contrast)", rel_max(S, P));
    }
}

// ---------------------------------------------------------------------------
// C. Measurement Jacobians and the complete update
// ---------------------------------------------------------------------------

// Noise-free residuals at the truth X for the estimate Xh.
V3 acc_residual(const G& Xh, const G& X, const V3& ka, double dT) {
    const V3 g(0, 0, GRAV);
    const V3 f = X.R.transpose() * (V3(X.X.col(2)) - g) + V3(X.B.col(1)) + ka * dT;
    return Xh.R * (f - V3(Xh.B.col(1)) - ka * dT) - (V3(Xh.X.col(2)) - g);
}
V3 mag_residual(const G& Xh, const G& X, const V3& Bw) { return Xh.R * (X.R.transpose() * Bw) - Bw; }

template<class Fn>
Eigen::Matrix<double, 3, NX> fd_meas(const G& Xh, Fn res, double eps) {
    Eigen::Matrix<double, 3, NX> H;
    for (int j = 0; j < NX; ++j) {
        TX dp = TX::Zero(), dm = TX::Zero();
        dp(j) = eps;
        dm(j) = -eps;
        H.col(j) = (res(truth_from(Xh, dp)) - res(truth_from(Xh, dm))) / (2 * eps);
    }
    return H;
}

struct NamedState {
    const char* name;
    M3 R;
    double s_w;
    V3 bg, ba;
};
std::vector<NamedState> jacobian_states() {
    return {
        {"identity, zero states", M3::Identity(), 0.0, V3::Zero(), V3::Zero()},
        {"small roll/pitch (3, -2 deg)", rot_euler(0.052, -0.035, 0.3), 0.3, V3(1e-3, -2e-3, 5e-4), V3(0.02, -0.01, 0.03)},
        {"large roll/pitch (40, -30 deg)", rot_euler(0.70, -0.52, -1.0), 0.6, V3(-3e-3, 1e-3, 2e-3), V3(-0.05, 0.04, 0.02)},
        {"large yaw (170 deg)", rot_euler(0.1, 0.05, 2.97), 0.6, V3(2e-3, 2e-3, -1e-3), V3(0.03, 0.03, -0.04)},
        {"nonzero CARMA states 2 m", rot_euler(-0.2, 0.15, 1.3), 2.0, V3(5e-3, -5e-3, 5e-3), V3(0.1, -0.1, 0.05)},
    };
}

void test_measurement_jacobians() {
    const V3 ka(0.002, -0.001, 0.003);
    const double dT = 7.0;
    const V3 Bw(18.0, 1.5, 43.0);
    const V3 g(0, 0, GRAV);
    Eigen::Matrix<double, 3, NX> Ha = Eigen::Matrix<double, 3, NX>::Zero(), Hm = Ha;
    Ha.block<3, 3>(0, 0) = lie::skew<double>(g);
    Ha.block<3, 3>(0, F::OFF_A) = -M3::Identity();
    Ha.block<3, 3>(0, F::OFF_BA) = -M3::Identity();
    Hm.block<3, 3>(0, 0) = -lie::skew<double>(Bw);
    for (const auto& st : jacobian_states()) {
        const G Xh = make_state(st.R, st.s_w, st.bg, st.ba);
        const auto Hfa = fd_meas(Xh, [&](const G& X) { return acc_residual(Xh, X, ka, dT); }, 1e-6);
        const auto Hfm = fd_meas(Xh, [&](const G& X) { return mag_residual(Xh, X, Bw); }, 1e-6);
        check(std::string("accel H = [[g]x .. -I(a) .. -I(b_a)] vs FD: ") + st.name,
              (Hfa - Ha).cwiseAbs().maxCoeff() / GRAV, 1e-8);
        check(std::string("mag H = [-[B_w]x ..] vs FD: ") + st.name, (Hfm - Hm).cwiseAbs().maxCoeff() / Bw.norm(), 1e-8);
    }
}

// Dense reference of one complete update (gain, Joseph, J_r transport, Exp_G
// injection) with finite-difference H and finite-difference J_r.
template<class Res>
void reference_update(G& Xh, MX& P, const V3& r, const Eigen::Matrix<double, 3, NX>& H, const M3& Rb, bool freeze_ba,
                      Res) {
    Eigen::Matrix<double, NX, 3> PHt = P * H.transpose();
    M3 S = H * PHt + Xh.R * Rb * Xh.R.transpose();
    S = 0.5 * (S + S.transpose()).eval();
    if (freeze_ba) {
        Eigen::Matrix<double, 3, NX> Hb = Eigen::Matrix<double, 3, NX>::Zero();
        Hb.block<3, 3>(0, F::OFF_BA) = H.block<3, 3>(0, F::OFF_BA);
        PHt -= P * Hb.transpose();
        PHt.block<3, 3>(F::OFF_BA, 0).setZero();
    }
    Eigen::Matrix<double, NX, 3> K = PHt * S.inverse();
    if (freeze_ba) K.block<3, 3>(F::OFF_BA, 0).setZero();
    const TX a = K * r;
    MX Pj = P + K * S * K.transpose() - K * PHt.transpose() - PHt * K.transpose();
    const MX J = fd_reset(a, 1e-7);
    P = J * Pj * J.transpose();
    P = 0.5 * (P + P.transpose()).eval();
    Xh = G::Exp_G(TX(-a)).compose(Xh);
}

void test_update_against_reference() {
    const V3 Bw(18.0, 1.5, 43.0);
    const M3 Racc = V3(0.03, 0.04, 0.05).array().square().matrix().asDiagonal();
    const M3 Rmag = V3(0.8, 0.8, 0.9).array().square().matrix().asDiagonal();
    for (int held = 0; held < 2; ++held) {
        for (const auto& st : jacobian_states()) {
            F f = make_filter();
            f.set_Racc_std(V3(0.03, 0.04, 0.05));
            f.set_Rmag(V3(0.8, 0.8, 0.9));
            f.set_accel_bias_temp_coeff(V3::Zero());
            f.state() = make_state(st.R, st.s_w, st.bg, st.ba);
            f.covariance_full() = scaled_cov(0.02, 1e-3, 0.05, 0.01);
            if (held) {
                f.set_acc_bias_updates_enabled(false);
            }
            G Xr = f.state();
            MX Pr = f.covariance_full();
            // A measurement from a nearby truth.
            TX xi = TX::Zero();
            xi.head<3>() = rvec(0.01);
            for (int i = 6; i < F::OFF_BA; ++i) xi(i) = std::normal_distribution<double>(0, 0.03)(rng);
            const G X = truth_from(Xr, xi);
            const V3 facc = X.R.transpose() * (V3(X.X.col(2)) - V3(0, 0, GRAV)) + V3(X.B.col(1));
            const V3 mmag = X.R.transpose() * Bw;
            // Accelerometer.
            f.measurement_update_acc_only(facc, 35.0);
            {
                const V3 r = Xr.R * (facc - V3(Xr.B.col(1))) - (V3(Xr.X.col(2)) - V3(0, 0, GRAV));
                const auto H = fd_meas(Xr, [&](const G& Xt) { return acc_residual(Xr, Xt, V3::Zero(), 0.0); }, 1e-6);
                reference_update(Xr, Pr, r, H, Racc, held != 0, 0);
            }
            const double ex = std::max({(f.state().R - Xr.R).cwiseAbs().maxCoeff(),
                                        (f.state().X - Xr.X).cwiseAbs().maxCoeff(),
                                        (f.state().B - Xr.B).cwiseAbs().maxCoeff()});
            const double ep = rel_max(f.covariance_full(), Pr);
            // Magnetometer.
            f.measurement_update_mag_only(mmag);
            {
                const V3 r = Xr.R * mmag - Bw;
                const auto H = fd_meas(Xr, [&](const G& Xt) { return mag_residual(Xr, Xt, Bw); }, 1e-6);
                reference_update(Xr, Pr, r, H, Rmag, held != 0, 0);
            }
            const double ex2 = std::max({(f.state().R - Xr.R).cwiseAbs().maxCoeff(),
                                         (f.state().X - Xr.X).cwiseAbs().maxCoeff(),
                                         (f.state().B - Xr.B).cwiseAbs().maxCoeff()});
            const double ep2 = rel_max(f.covariance_full(), Pr);
            const std::string tag = std::string(held ? " (bias held)" : "") + ": " + st.name;
            check("update: accel state vs dense FD reference" + tag, ex, 1e-7);
            check("update: accel covariance vs dense FD reference" + tag, ep, 1e-6);
            check("update: accel+mag state vs dense FD reference" + tag, ex2, 1e-7);
            check("update: accel+mag covariance vs dense FD reference" + tag, ep2, 1e-6);
        }
    }
}

// ---------------------------------------------------------------------------
// D. Propagation
// ---------------------------------------------------------------------------

void test_transition() {
    const M3 Qg = V3(2e-5, 3e-5, 2.5e-5).asDiagonal();
    const M3 Qbg = V3(1e-7, 2e-7, 1.5e-7).asDiagonal();
    const M3 Qba = V3::Constant(2.5e-7).asDiagonal();
    for (const auto& st : jacobian_states()) {
        for (double h : {0.005, 0.05}) {
            for (double rate : {0.0, 0.4, 1.3}) {
                F f = make_filter();
                f.set_gyro_noise_density_rad_sqrt_s(Qg.diagonal().cwiseSqrt());
                f.set_Q_bgyro_rw(Qbg.diagonal());
                f.set_Q_bacc_rw(Qba.diagonal().cwiseSqrt());
                f.state() = make_state(st.R, st.s_w, st.bg, st.ba);
                const V3 wm = V3(V3(0.6, -0.5, 0.62).normalized() * rate) + st.bg;
                MX Phi;
                f.build_transition(wm, h, Phi);
                const Flow fl = make_flow(f, h);
                const MX Jfd = fd_transition(f.state(), wm, fl, TX::Zero(), 1e-6);
                char tag[160];
                std::snprintf(tag, sizeof tag, " (%s, h %.3f, |w| %.1f)", st.name, h, rate);
                check(std::string("transition: Phi vs FD of the nonlinear error flow") + tag,
                      (Phi - Jfd).cwiseAbs().maxCoeff(), 2e-8);
                const ContModel cm = cont_model(f, wm, Qg, Qbg, Qba);
                MX Pr, Qr;
                rk4_reference(cm, h, 400, Pr, Qr);
                check(std::string("transition: Phi vs RK4 of continuous model") + tag,
                      (Phi - Pr).cwiseAbs().maxCoeff(), 1e-10);
                MX Qd;
                f.build_process_noise(wm, h, Qd);
                check(std::string("process noise: Qd vs RK4 Lyapunov (corr units)") + tag, rel_max(Qd, Qr), 1e-7);
            }
        }
    }
    // Attitude -> translation coupling is identically zero with one A.
    {
        F f = make_filter();
        f.state() = make_state(rot_euler(0.5, -0.3, 1.0), 1.0, rvec(1e-3), rvec(0.02));
        MX Phi;
        f.build_transition(V3(0.3, -0.2, 0.5), 0.02, Phi);
        check("transition: Phi_{rho,phi} = 0 (isotropic CARMA)", Phi.block<12, 3>(6, 0).cwiseAbs().maxCoeff(), 1e-15);
        // Collapse identity Phi_{rho,bg} = [W(h)]x Phi_{phi,bg}.
        const Flow fl = make_flow(f, 0.02);
        const G X1 = flow(f.state(), V3(0.3, -0.2, 0.5), fl);
        double e = 0;
        for (int i = 0; i < NC; ++i)
            e = std::max(e, (Phi.block<3, 3>(6 + 3 * i, 3) -
                             lie::skew<double>(V3(X1.X.col(i))) * Phi.block<3, 3>(0, 3)).cwiseAbs().maxCoeff());
        check("transition: Phi_{rho_i,bg} = [w_i(h)]x Phi_{phi,bg}", e, 1e-14);
    }
}

// The structured Gauss-Legendre impulse-response quadrature of
// Kalman3D_Wave_TFG, applied to this prior: Qd = int Phi(h,s) L(s) Qc L(s)^T
// Phi(h,s)^T ds with Phi(h,s) the transition from the nominal state at s.
void test_quadrature_collapse() {
    const M3 Qg = V3(2e-5, 3e-5, 2.5e-5).asDiagonal();
    const M3 Qbg = V3(1e-7, 2e-7, 1.5e-7).asDiagonal();
    F f = make_filter();
    f.set_gyro_noise_density_rad_sqrt_s(Qg.diagonal().cwiseSqrt());
    f.set_Q_bgyro_rw(Qbg.diagonal());
    f.set_acc_bias_updates_enabled(false);   // only gyro, gyro-bias and CARMA noise
    f.state() = make_state(rot_euler(0.4, -0.3, 2.0), 1.0, V3(1e-3, 2e-3, -1e-3), V3::Zero());
    const V3 wm(0.5, -0.8, 0.9);
    const double h = 0.01;
    MX Qd;
    f.build_process_noise(wm, h, Qd);
    static const double xg[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
    static const double wg[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
                                 0.2369268850561891};
    MX Qq = MX::Zero();
    const ContModel cm = cont_model(f, wm, Qg, Qbg, M3::Zero());
    for (int q = 0; q < 5; ++q) {
        const double s = 0.5 * h * (1 + xg[q]);
        const double w = 0.5 * h * wg[q];
        // Nominal state at s, transition s -> h from it.
        F fs = f;
        fs.state() = flow(f.state(), wm, make_flow(f, s));
        MX Phi_hs;
        fs.build_transition(wm, h - s, Phi_hs);
        MX Fm, LQL;
        cont_F(cm, s, Fm, LQL);
        Qq += w * Phi_hs * LQL * Phi_hs.transpose();
    }
    // Entries with an attitude or bias index carry only gyro and gyro-bias
    // noise, which GL-5 integrates to rounding; the world-world block also
    // holds the CARMA part, which it does not integrate exactly.
    double e_gyro = 0;
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j) {
            const bool world_i = i >= 6 && i < F::OFF_BA, world_j = j >= 6 && j < F::OFF_BA;
            if (world_i && world_j) continue;
            e_gyro = std::max(e_gyro, std::abs(Qd(i, j) - Qq(i, j)) / std::sqrt(Qq(i, i) * Qq(j, j)));
        }
    check("quadrature: collapsed Qd vs structured GL-5 impulse response (gyro/bias rows)", e_gyro, 1e-9);
    check("quadrature: collapsed Qd vs structured GL-5 (full, corr units)", rel_max(Qd, Qq), 1e-6);
}

void test_semigroup_and_linearization() {
    F f = make_filter();
    f.state() = make_state(rot_euler(0.3, -0.2, 0.7), 0.8, V3(1e-3, -1e-3, 2e-3), V3(0.02, 0.01, -0.03));
    const V3 wm(0.2, 0.4, -0.3);
    const double h1 = 0.004, h2 = 0.007;
    MX P1, P2, P12, Q1, Q2, Q12;
    f.build_transition(wm, h1, P1);
    f.build_process_noise(wm, h1, Q1);
    f.build_transition(wm, h1 + h2, P12);
    f.build_process_noise(wm, h1 + h2, Q12);
    F f2 = f;
    f2.state() = flow(f.state(), wm, make_flow(f, h1));
    f2.build_transition(wm, h2, P2);
    f2.build_process_noise(wm, h2, Q2);
    check("semigroup: Phi(h1+h2) = Phi(h2) Phi(h1), non-uniform steps", (P12 - P2 * P1).cwiseAbs().maxCoeff(), 1e-12);
    check("semigroup: Qd(h1+h2) = Phi(h2) Qd(h1) Phi(h2)^T + Qd(h2)",
          rel_max(MX(P2 * Q1 * P2.transpose() + Q2), Q12), 1e-8);
    // Linearization error of one step as a function of the error size.
    const Flow fl = make_flow(f, 0.005);
    MX Phi;
    f.build_transition(wm, 0.005, Phi);
    const G X1h = flow(f.state(), wm, fl);
    for (double s : {1e-4, 1e-3, 1e-2, 1e-1}) {
        double worst = 0;
        for (int t = 0; t < 20; ++t) {
            TX xi;
            xi.head<3>() = rvec(s);
            xi.segment<3>(3) = rvec(s * 0.1);
            for (int i = 6; i < NX; ++i) xi(i) = std::normal_distribution<double>(0, s)(rng);
            const TX x1 = err(X1h, flow(truth_from(f.state(), xi), wm, fl));
            worst = std::max(worst, (x1 - Phi * xi).norm() / xi.norm());
        }
        char n[96];
        std::snprintf(n, sizeof n, "linearization: |xi(h) - Phi xi(0)|/|xi|, |xi| ~ %.0e", s);
        report(n, worst);
    }
}

void test_structured_propagation() {
    F f = make_filter();
    f.state() = make_state(rot_euler(0.5, 0.2, -2.2), 1.0, V3(2e-3, -1e-3, 1e-3), V3(0.03, -0.02, 0.01));
    f.covariance_full() = scaled_cov(0.02, 1e-3, 0.05, 0.01);
    const V3 wm(-0.3, 0.6, 0.2);
    for (double h : {0.005, 0.0063}) {
        MX Phi, Qd;
        f.build_transition(wm, h, Phi);
        f.build_process_noise(wm, h, Qd);
        const MX expect = Phi * f.covariance_full() * Phi.transpose() + Qd;
        F g = f;
        g.time_update(wm, h);
        char n[96];
        std::snprintf(n, sizeof n, "time_update: structured P vs dense Phi P Phi^T + Qd, h %.4f", h);
        check(n, rel_max(g.covariance_full(), expect), 1e-12);
        const G X1 = flow(f.state(), wm, make_flow(f, h));
        check(std::string(n) + " (mean)", std::max((g.state().R - X1.R).cwiseAbs().maxCoeff(),
                                                   (g.state().X - X1.X).cwiseAbs().maxCoeff()), 1e-12);
    }
}

// ---------------------------------------------------------------------------
// E. Stochastic checks
// ---------------------------------------------------------------------------

void test_monte_carlo_propagation() {
    // Nonlinear stochastic truth (20 Euler-Maruyama substeps per filter step)
    // from errors xi0 ~ N(0, P0); the filter mean propagates noise-free.
    const M3 Qg = M3::Identity() * 4e-5;
    const M3 Qbg = M3::Identity() * 4e-6;
    const M3 Qba = M3::Identity() * 1e-4;
    F f = make_filter(0.5, 0.4, 0.6, 0.3);
    f.set_gyro_noise_density_rad_sqrt_s(Qg.diagonal().cwiseSqrt());
    f.set_Q_bgyro_rw(Qbg.diagonal());
    f.set_Q_bacc_rw(Qba.diagonal().cwiseSqrt());
    f.set_acc_bias_time_constant(20.0);
    f.state() = make_state(rot_euler(0.3, -0.25, 1.1), 0.5, V3(2e-3, -3e-3, 1e-3), V3(0.02, -0.01, 0.03));
    f.covariance_full() = scaled_cov(3e-3, 2e-3, 0.01, 0.01);
    const MX P0 = f.covariance_full();
    const double h = 0.005;
    const int steps = 200, sub = 20;
    auto wm_at = [](int k) { return V3(0.4 * std::sin(0.05 * k), 0.3 * std::cos(0.07 * k), 0.5); };
    std::vector<G> nominal(steps + 1);
    nominal[0] = f.state();
    F fp = f;
    for (int k = 0; k < steps; ++k) {
        fp.time_update(wm_at(k), h);
        nominal[k + 1] = fp.state();
    }
    const MX Pk = fp.covariance_full();
    const int N = 4000;
    Eigen::LLT<MX> llt(P0);
    const MX L0 = llt.matrixL();
    MC Fs, Qs;
    kc::discretize_lti<double, NC>(f.vertical_model().A(), f.vertical_model().W(), h / sub, Fs, Qs);
    Eigen::LLT<MC> lq(Qs);
    const MC Lq = lq.matrixL();
    Eigen::LLT<M3> ls(f.spatial_sigma());
    const M3 Ls = ls.matrixL();
    const double phi_s = std::exp(-h / sub / f.get_acc_bias_time_constant());
    const double qba_s = 0.5 * f.get_acc_bias_time_constant() * (1 - phi_s * phi_s);
    std::normal_distribution<double> n(0, 1);
    MX S = MX::Zero();
    TX mean = TX::Zero();
    for (int p = 0; p < N; ++p) {
        TX z;
        for (int i = 0; i < NX; ++i) z(i) = n(rng);
        G X = truth_from(nominal[0], L0 * z);
        for (int k = 0; k < steps; ++k) {
            const V3 wm = wm_at(k);
            for (int j = 0; j < sub; ++j) {
                const double ds = h / sub;
                const V3 ng(n(rng), n(rng), n(rng));
                X.R = X.R * lie::Exp<double>(V3((wm - V3(X.B.col(0)) - Qg.diagonal().cwiseSqrt().cwiseProduct(ng) / std::sqrt(ds)) * ds));
                X.B.col(0) += Qbg.diagonal().cwiseSqrt().cwiseProduct(V3(n(rng), n(rng), n(rng))) * std::sqrt(ds);
                Eigen::Matrix<double, NC, 3> e;
                for (int a = 0; a < NC; ++a)
                    for (int b = 0; b < 3; ++b) e(a, b) = n(rng);
                X.X = X.X * Fs.transpose() + Ls * (Lq * e).transpose();
                X.B.col(1) = X.B.col(1) * phi_s +
                             Qba.diagonal().cwiseSqrt().cwiseProduct(V3(n(rng), n(rng), n(rng))) * std::sqrt(qba_s);
            }
        }
        const TX xi = err(nominal[steps], X);
        mean += xi;
        S += xi * xi.transpose();
    }
    mean /= N;
    S /= N;
    S -= mean * mean.transpose();
    check("Monte Carlo: propagated covariance vs sample covariance (corr units)", rel_max(S, Pk), 0.08);
    double tr_ratio = 0;
    for (int i = 0; i < NX; ++i) tr_ratio = std::max(tr_ratio, std::abs(S(i, i) / Pk(i, i) - 1.0));
    check("Monte Carlo: per-state variance ratio |S_ii/P_ii - 1|", tr_ratio, 0.10);
    double mz = 0;
    for (int i = 0; i < NX; ++i) mz = std::max(mz, std::abs(mean(i)) / std::sqrt(Pk(i, i) / N));
    check("Monte Carlo: mean error / standard error", mz, 4.5);
}

void test_stationary_prior() {
    for (int prec = 0; prec < 2; ++prec) {
        double worst = 0;
        if (prec == 0) {
            F f = make_filter(0.5, 0.4, 0.6, 0.3);
            const MX P0 = f.covariance_full();
            for (int k = 0; k < 20000; ++k) f.time_update(V3::Zero(), 0.005);
            MX P = f.covariance_full();
            for (int i = 6; i < F::OFF_BA; ++i)
                for (int j = 6; j < F::OFF_BA; ++j)
                    worst = std::max(worst, std::abs(P(i, j) - P0(i, j)) / std::sqrt(P0(i, i) * P0(j, j)));
            check("stationary prior Sigma_s (x) P_inf invariant over 100 s (double)", worst, 1e-9);
        } else {
            using Ff = Kalman3D_Wave_TFG_CARMA3<float, 4>;
            Ff f(Eigen::Vector3f::Constant(0.03f), Eigen::Vector3f::Constant(3e-4f), Eigen::Vector3f::Constant(0.8f));
            f.set_vertical_carma(float(KAPPA), 0.6f);
            f.set_horizontal_carma(float(KAPPA), 0.5f, 0.4f, 0.3f);
            f.reset_aw_covariance_to_stationary();
            const auto P0 = f.covariance_full();
            for (int k = 0; k < 20000; ++k) f.time_update(Eigen::Vector3f::Zero(), 0.005f);
            const auto P = f.covariance_full();
            for (int i = 6; i < Ff::OFF_BA; ++i)
                for (int j = 6; j < Ff::OFF_BA; ++j)
                    worst = std::max(worst, double(std::abs(P(i, j) - P0(i, j))) /
                                                std::sqrt(double(P0(i, i)) * double(P0(j, j))));
            check("stationary prior invariant over 100 s (float)", worst, 2e-4);
        }
    }
}

// Long realization of the three-axis prior through the filter's own
// discretization (double): stationary covariance and Welch spectra (each axis
// and the x-y co-spectrum) against the analytic polynomial PSD.
void test_process_spectrum() {
    const double sx = 0.5, sy = 0.4, sz = 0.6, r = 0.3;
    F f = make_filter(sx, sy, sz, r);
    const auto& m = f.vertical_model();
    MC P1;
    m.stationary_covariance(P1);
    const double h = 0.02;
    MC Fd, Qd;
    kc::discretize_lti<double, NC>(m.A(), m.W(), h, Fd, Qd);
    Eigen::LLT<MC> lq(Qd);
    const MC Lq = lq.matrixL();
    Eigen::LLT<M3> ls(f.spatial_sigma());
    const M3 Ls = ls.matrixL();
    const int seg = 4096, nseg = 256;
    const long N = static_cast<long>(seg) * nseg;
    std::normal_distribution<double> n(0, 1);
    Eigen::Matrix<double, 3, NC> X = Eigen::Matrix<double, 3, NC>::Zero();
    {   // stationary start
        Eigen::LLT<MC> l1(P1);
        Eigen::Matrix<double, NC, 3> e;
        for (int a = 0; a < NC; ++a)
            for (int b = 0; b < 3; ++b) e(a, b) = n(rng);
        X = Ls * (MC(l1.matrixL()) * e).transpose();
    }
    std::vector<double> px(N), py(N), pz(N), az(N);
    Eigen::Matrix<double, 12, 12> C = Eigen::Matrix<double, 12, 12>::Zero();
    for (long k = 0; k < N; ++k) {
        Eigen::Matrix<double, NC, 3> e;
        for (int a = 0; a < NC; ++a)
            for (int b = 0; b < 3; ++b) e(a, b) = n(rng);
        X = X * Fd.transpose() + Ls * (Lq * e).transpose();
        px[k] = X(0, 0); py[k] = X(1, 0); pz[k] = X(2, 0); az[k] = X(2, 2);
        Eigen::Matrix<double, 12, 1> v;
        for (int c = 0; c < NC; ++c) v.segment<3>(Eigen::Index(3) * c) = X.col(c);
        C += v * v.transpose();
    }
    C /= double(N);
    double worst = 0;
    for (int a = 0; a < NC; ++a)
        for (int b = 0; b < NC; ++b)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l) {
                    const double expect = P1(a, b) * f.spatial_sigma()(k, l);
                    const double sc = std::sqrt(P1(a, a) * f.spatial_sigma()(k, k) * P1(b, b) * f.spatial_sigma()(l, l));
                    worst = std::max(worst, std::abs(C(3 * a + k, 3 * b + l) - expect) / sc);
                }
    check("process: sample covariance vs Sigma_s (x) P_inf (corr units, 5.8 h)", worst, 0.05);
    // Welch (Hann) one-sided PSD at a set of frequencies.
    auto welch = [&](const std::vector<double>& x, const std::vector<double>& y, double fhz) {
        std::complex<double> acc = 0;
        double wss = 0;
        for (int j = 0; j < seg; ++j) {
            const double w = 0.5 - 0.5 * std::cos(2 * M_PI * j / seg);
            wss += w * w;
        }
        for (int s = 0; s < nseg; ++s) {
            std::complex<double> X1 = 0, Y1 = 0;
            for (int j = 0; j < seg; ++j) {
                const double w = 0.5 - 0.5 * std::cos(2 * M_PI * j / seg);
                const std::complex<double> e = std::polar(1.0, -2 * M_PI * fhz * j * h);
                X1 += w * x[s * seg + j] * e;
                Y1 += w * y[s * seg + j] * e;
            }
            acc += X1 * std::conj(Y1);
        }
        return 2.0 * h * acc / (nseg * wss);
    };
    const double gain2 = m.W()(2, 2);
    auto analytic_p = [&](double fhz) { return 2.0 * gain2 * kc::carma_unit_psd<NC>(m.shape(), m.kappa(), fhz); };
    // 256 Hann segments: std of a log-PSD estimate ~ 0.27 dB, of the co-spectrum
    // ratio ~ sqrt((1 + r^2)/512) = 0.046 per frequency.
    double e_db = 0, e_db_a = 0, m_db = 0, e_coh = 0, m_coh = 0;
    int nf = 0, nd = 0;
    for (int fi = 0; fi < 14; ++fi) {
        const double fhz = 0.12 + 0.06 * fi;
        const double sp = analytic_p(fhz);
        const double w4 = std::pow(2 * M_PI * fhz, 4);
        for (int ax = 0; ax < 3; ++ax) {
            const auto& x = ax == 0 ? px : (ax == 1 ? py : pz);
            const double est = std::real(welch(x, x, fhz));
            const double sig2 = f.spatial_sigma()(ax, ax);
            const double d = std::abs(10 * std::log10(est / (sig2 * sp)));
            e_db = std::max(e_db, d);
            m_db += d;
            ++nd;
        }
        const double esta = std::real(welch(az, az, fhz));
        e_db_a = std::max(e_db_a, std::abs(10 * std::log10(esta / (sz * sz * w4 * sp))));
        const auto cxy = welch(px, py, fhz);
        const double dc = std::real(cxy) / (sx * sy * sp) - r;
        e_coh = std::max(e_coh, std::abs(dc));
        m_coh += dc;
        ++nf;
    }
    check("process: Welch PSD of p_x, p_y, p_z vs analytic, mean |dB| (0.12-0.9 Hz)", m_db / nd, 0.35);
    check("process: Welch PSD of p_x, p_y, p_z vs analytic, max |dB|", e_db, 1.2);
    check("process: Welch PSD of a_z vs w^4 S_p, max |dB|", e_db_a, 1.2);
    check("process: x-y co-spectrum / (sx sy S_p) - r, mean over frequency", std::abs(m_coh / nf), 0.04);
    report("process: x-y co-spectrum / (sx sy S_p) - r, max over frequency", e_coh);
    // The DC zero of CARMA(4,1): S_p ~ f^2 as f -> 0, and A is Hurwitz.
    const double r1 = analytic_p(1e-3) / analytic_p(2e-3);
    check("process: DC zero, S_p(1 mHz)/S_p(2 mHz) = 1/4", std::abs(r1 - 0.25), 1e-4);
    Eigen::EigenSolver<MC> es(m.A());
    report("process: largest pole real part (rad/s)", es.eigenvalues().real().maxCoeff());
    check("process: CARMA(4,1) Hurwitz", es.eigenvalues().real().maxCoeff(), -1e-3);
}

// Synthetic IMU stream (waves + rotation) shared by the long-run checks.
struct Imu {
    V3 gyr, acc, mag;
};
struct Synth {
    double t = 0;
    V3 Bw = V3(18.0, 1.5, 43.0);
    bool turns = false;
    double att(double tt, int k) const {
        const double base[3] = {0.12 * std::sin(1.9 * tt + 0.3), 0.08 * std::sin(1.4 * tt + 1.1),
                                0.6 + 0.1 * std::sin(0.11 * tt)};
        double v = base[k];
        if (k == 2 && turns) v += 1.5 * std::sin(0.08 * tt);
        return v;
    }
    M3 R(double tt) const { return rot_euler(att(tt, 0), att(tt, 1), att(tt, 2)); }
    V3 a(double tt) const {
        return V3(-0.3 * 2.0 * 2.0 * std::sin(2.0 * tt), -0.2 * 1.7 * 1.7 * std::sin(1.7 * tt + 0.5),
                  -0.4 * 1.95 * 1.95 * std::sin(1.95 * tt + 1.0));
    }
    Imu sample(double dt, std::mt19937& g, double noise) const {
        std::normal_distribution<double> n(0, 1);
        const M3 R0 = R(t), R1 = R(t + dt);
        Imu s;
        s.gyr = lie::Log<double>(M3(R0.transpose() * R1)) / dt + V3(n(g), n(g), n(g)) * 1e-3 * noise;
        s.acc = R1.transpose() * (a(t + dt) - V3(0, 0, GRAV)) + V3(n(g), n(g), n(g)) * 0.03 * noise;
        s.mag = R1.transpose() * Bw + V3(n(g), n(g), n(g)) * 0.5 * noise;
        return s;
    }
};

void test_float_vs_double() {
    using Ff = Kalman3D_Wave_TFG_CARMA3<float, 4>;
    Ff ff(Eigen::Vector3f::Constant(0.03f), Eigen::Vector3f::Constant(3e-4f), Eigen::Vector3f::Constant(0.8f), 5e-4f,
          1e-6f, 1e-8f, 1.5f, float(GRAV));
    F fd(V3::Constant(0.03), V3::Constant(3e-4), V3::Constant(0.8), 5e-4, 1e-6, 1e-8, 1.5, GRAV);
    ff.set_vertical_carma(float(KAPPA), 0.6f);
    ff.set_horizontal_carma(float(KAPPA), 0.5f, 0.5f, 0.0f);
    fd.set_vertical_carma(KAPPA, 0.6);
    fd.set_horizontal_carma(KAPPA, 0.5, 0.5, 0.0);
    ff.reset_aw_covariance_to_stationary();
    fd.reset_aw_covariance_to_stationary();
    Synth sy;
    ff.set_mag_world_ref(sy.Bw.cast<float>());
    fd.set_mag_world_ref(sy.Bw);
    const M3 R0 = sy.R(0);
    ff.initialize_from_attitude(Eigen::Quaternionf(R0.cast<float>()), 0.02f, 0.1f);
    fd.initialize_from_attitude(Eigen::Quaterniond(R0), 0.02, 0.1);
    std::mt19937 g(7);
    std::uniform_real_distribution<double> jit(-0.25, 0.25);
    double worst_ang = 0, worst_p = 0, worst_P = 0;
    for (int k = 0; k < 200 * 60; ++k) {
        const double dt = 0.005 * (1 + jit(g));
        const Imu s = sy.sample(dt, g, 1.0);
        sy.t += dt;
        ff.time_update(s.gyr.cast<float>(), float(dt));
        fd.time_update(s.gyr, dt);
        ff.measurement_update_acc_only(s.acc.cast<float>());
        fd.measurement_update_acc_only(s.acc);
        if (k % 8 == 0) {
            ff.measurement_update_mag_only(s.mag.cast<float>());
            fd.measurement_update_mag_only(s.mag);
        }
        if (k > 200 * 10) {
            worst_ang = std::max(worst_ang, lie::Log<double>(M3(ff.state().R.cast<double>() * fd.state().R.transpose())).norm());
            worst_p = std::max(worst_p, (ff.get_position().cast<double>() - fd.get_position()).norm());
            worst_P = std::max(worst_P, rel_max(MX(ff.covariance_full().cast<double>()), fd.covariance_full()));
        }
    }
    check("float vs double: attitude difference over 60 s, +-25 % jitter (rad)", worst_ang, 1e-4);
    check("float vs double: position difference (m)", worst_p, 1e-3);
    check("float vs double: covariance difference (corr units)", worst_P, 0.02);
    check("float: covariance PSD (min corr eigenvalue)", -min_eig_corr(MX(ff.covariance_full().cast<double>())), 1e-5);
}

// ---------------------------------------------------------------------------
// F. Frame changes and equivariance
// ---------------------------------------------------------------------------

MX world_rotation_map(const M3& Q) {
    MX Gm = MX::Zero();
    for (int b = 0; b < F::NBLK; ++b) Gm.block<3, 3>(Eigen::Index(3) * b, Eigen::Index(3) * b) = Q;
    return Gm;
}

void test_equivariance() {
    // CARMA model: (Phi_C (x) I) commutes with (I (x) Q) for every rotation.
    F f0 = make_filter(0.5, 0.4, 0.6, 0.3);
    MC Fc, Qc;
    kc::discretize_lti<double, NC>(f0.vertical_model().A(), f0.vertical_model().W(), 0.005, Fc, Qc);
    double e_c = 0, e_q = 0;
    for (int t = 0; t < 10; ++t) {
        const M3 Q = lie::Exp<double>(rvec(1.5));
        const Eigen::Matrix<double, 3, NC> W = Eigen::Matrix<double, 3, NC>::Random();
        e_c = std::max(e_c, ((Q * W) * Fc.transpose() - Q * (W * Fc.transpose())).cwiseAbs().maxCoeff());
        // Noise transforms as Sigma -> Q Sigma Q^T: Qd(Q S Q^T) = (I (x) Q) Qd(S) (I (x) Q)^T.
        const M3& S = f0.spatial_sigma();
        Eigen::Matrix<double, 12, 12> A1, A2, Gq = Eigen::Matrix<double, 12, 12>::Zero();
        for (int c = 0; c < NC; ++c) Gq.block<3, 3>(Eigen::Index(3) * c, Eigen::Index(3) * c) = Q;
        const M3 S2 = Q * S * Q.transpose();
        for (int a = 0; a < NC; ++a)
            for (int b = 0; b < NC; ++b) {
                A1.block<3, 3>(Eigen::Index(3) * a, Eigen::Index(3) * b) = Qc(a, b) * S;
                A2.block<3, 3>(Eigen::Index(3) * a, Eigen::Index(3) * b) = Qc(a, b) * S2;
            }
        e_q = std::max(e_q, (A2 - Gq * A1 * Gq.transpose()).cwiseAbs().maxCoeff() / Qc.cwiseAbs().maxCoeff());
    }
    check("equivariance: CARMA transition commutes with any world rotation", e_c, 1e-14);
    check("equivariance: Qd(Q S Q^T) = (I (x) Q) Qd(S) (I (x) Q)^T", e_q, 1e-13);

    // Complete filter under a world yaw gauge: frame B = Rz(psi) frame A.
    const double psi = 1.1;
    const M3 Q = Eigen::AngleAxisd(psi, V3::UnitZ()).toRotationMatrix();
    auto run_pair = [&](double sx, double sy, double r, bool rotate_sigma) {
        F fa = make_filter(sx, sy, 0.6, r), fb = make_filter(sx, sy, 0.6, r);
        if (rotate_sigma) {
            const Eigen::Matrix2d Sh = fa.spatial_sigma().topLeftCorner<2, 2>();
            const Eigen::Matrix2d Sr = Q.topLeftCorner<2, 2>() * Sh * Q.topLeftCorner<2, 2>().transpose();
            const double ax = std::sqrt(Sr(0, 0)), ay = std::sqrt(Sr(1, 1));
            fb.set_horizontal_carma(KAPPA, ax, ay, Sr(0, 1) / (ax * ay));
        }
        Synth sy_;
        fa.set_mag_world_ref(sy_.Bw);
        fb.set_mag_world_ref(Q * sy_.Bw);
        fa.state() = make_state(sy_.R(0), 0.2, V3(1e-3, 2e-3, -1e-3), V3(0.01, -0.02, 0.01));
        fb.state() = fa.state();
        fb.state().R = Q * fa.state().R;
        fb.state().X = Q * fa.state().X;
        fa.covariance_full() = scaled_cov(0.01, 1e-3, 0.05, 0.01);
        fb.covariance_full() = world_rotation_map(Q) * fa.covariance_full() * world_rotation_map(Q).transpose();
        std::mt19937 g(3);
        for (int k = 0; k < 200 * 20; ++k) {
            const Imu s = sy_.sample(0.005, g, 1.0);
            sy_.t += 0.005;
            fa.time_update(s.gyr, 0.005);
            fb.time_update(s.gyr, 0.005);
            fa.measurement_update_acc_only(s.acc);
            fb.measurement_update_acc_only(s.acc);
            if (k % 8 == 0) {
                fa.measurement_update_mag_only(s.mag);
                fb.measurement_update_mag_only(s.mag);
            }
        }
        const double eR = (fb.state().R - Q * fa.state().R).cwiseAbs().maxCoeff();
        const double eX = (fb.state().X - Q * fa.state().X).cwiseAbs().maxCoeff();
        const MX Pexp = world_rotation_map(Q) * fa.covariance_full() * world_rotation_map(Q).transpose();
        return std::max({eR, eX, rel_max(fb.covariance_full(), Pexp)});
    };
    check("equivariance: filter under world yaw, isotropic horizontal prior", run_pair(0.5, 0.5, 0.0, false), 1e-9);
    check("equivariance: anisotropic Sigma_s rotated with the frame", run_pair(0.6, 0.35, 0.4, true), 1e-9);
    report("non-invariance: anisotropic Sigma_s NOT rotated (expected nonzero)", run_pair(0.6, 0.35, 0.4, false));
}

// Anisotropic dynamics (different A per world axis): the rotation no longer
// commutes with the transition and the invariant error acquires an attitude ->
// translation term rho' = A rho + C(W) phi, C(W) phi = A([phi]x W) - [phi]x A(W).
void test_anisotropic_dynamics() {
    kc::CarmaVerticalModel<double, NC> mh(kc::CARMA41_JONSWAP_SHAPE), mz(kc::CARMA41_JONSWAP_SHAPE);
    mh.configure(KAPPA * 0.8, 1.0);
    mz.configure(KAPPA, 1.0);
    const MC Ah = mh.A(), Az = mz.A();
    auto apply_A = [&](const Eigen::Matrix<double, 3, NC>& W) {
        Eigen::Matrix<double, 3, NC> Y;
        for (int k = 0; k < 3; ++k) Y.row(k) = ((k < 2 ? Ah : Az) * W.row(k).transpose()).transpose();
        return Y;
    };
    Eigen::Matrix<double, 3, NC> W = Eigen::Matrix<double, 3, NC>::Random();
    const V3 phi = rvec(1.0);
    Eigen::Matrix<double, 3, NC> SW, ASW;
    for (int c = 0; c < NC; ++c) SW.col(c) = lie::skew<double>(phi) * W.col(c);
    const Eigen::Matrix<double, 3, NC> lhs = apply_A(SW);
    Eigen::Matrix<double, 3, NC> rhs;
    const Eigen::Matrix<double, 3, NC> AW = apply_A(W);
    for (int c = 0; c < NC; ++c) rhs.col(c) = lie::skew<double>(phi) * AW.col(c);
    report("anisotropic A: |A([phi]x W) - [phi]x A(W)| / |A W| (commutator)", (lhs - rhs).norm() / AW.norm());
    // FD of the nonlinear error flow with per-axis dynamics: d rho(h)/d phi(0)
    // against int_0^h Phi_A(h-t) C(W(t)) dt (RK4 of rho' = A rho + C(W(t)) phi).
    const double h = 0.05;
    MC Fh, Fz, Qtmp;
    kc::discretize_lti<double, NC>(Ah, MC::Zero(), h, Fh, Qtmp);
    kc::discretize_lti<double, NC>(Az, MC::Zero(), h, Fz, Qtmp);
    auto flowW = [&](const Eigen::Matrix<double, 3, NC>& W0) {
        Eigen::Matrix<double, 3, NC> Y;
        for (int k = 0; k < 3; ++k) Y.row(k) = ((k < 2 ? Fh : Fz) * W0.row(k).transpose()).transpose();
        return Y;
    };
    // Truth W = Exp(-phi)(What - rho) with rho = 0; error rho(h) = What(h) - Exp(phi) W(h).
    Eigen::Matrix<double, 12, 3> Jfd;
    const double eps = 1e-6;
    const Eigen::Matrix<double, 3, NC> W1h = flowW(W);
    for (int j = 0; j < 3; ++j) {
        V3 dp = V3::Zero();
        dp(j) = eps;
        auto rho_h = [&](const V3& ph) {
            const Eigen::Matrix<double, 3, NC> Wt = lie::Exp<double>(V3(-ph)) * W;
            const Eigen::Matrix<double, 3, NC> R1 = W1h - lie::Exp<double>(ph) * flowW(Wt);
            Eigen::Matrix<double, 12, 1> v;
            for (int c = 0; c < NC; ++c) v.segment<3>(Eigen::Index(3) * c) = R1.col(c);
            return v;
        };
        Jfd.col(j) = (rho_h(dp) - rho_h(V3(-dp))) / (2 * eps);
    }
    // RK4 of Y' = A_blk Y + C(W(t)), Y(0) = 0 (12 x 3).
    Eigen::Matrix<double, 12, 3> Y = Eigen::Matrix<double, 12, 3>::Zero();
    const int n = 2000;
    const double dt = h / n;
    auto Ct = [&](double t) {
        MC Ft1, Ft2;
        kc::discretize_lti<double, NC>(Ah, MC::Zero(), t, Ft1, Qtmp);
        kc::discretize_lti<double, NC>(Az, MC::Zero(), t, Ft2, Qtmp);
        Eigen::Matrix<double, 3, NC> Wt;
        for (int k = 0; k < 3; ++k) Wt.row(k) = ((k < 2 ? Ft1 : Ft2) * W.row(k).transpose()).transpose();
        const Eigen::Matrix<double, 3, NC> AWt = apply_A(Wt);
        Eigen::Matrix<double, 12, 3> Cm;
        for (int j = 0; j < 3; ++j) {
            V3 e = V3::Zero();
            e(j) = 1.0;
            Eigen::Matrix<double, 3, NC> Se;
            for (int c = 0; c < NC; ++c) Se.col(c) = lie::skew<double>(e) * Wt.col(c);
            const Eigen::Matrix<double, 3, NC> ASe = apply_A(Se);
            for (int c = 0; c < NC; ++c) Cm.block<3, 1>(Eigen::Index(3) * c, j) = ASe.col(c) - lie::skew<double>(e) * AWt.col(c);
        }
        return Cm;
    };
    auto Ablk = [&](const Eigen::Matrix<double, 12, 3>& Z) {
        Eigen::Matrix<double, 12, 3> out;
        for (int j = 0; j < 3; ++j) {
            Eigen::Matrix<double, 3, NC> Wz;
            for (int c = 0; c < NC; ++c) Wz.col(c) = Z.block<3, 1>(Eigen::Index(3) * c, j);
            const Eigen::Matrix<double, 3, NC> AWz = apply_A(Wz);
            for (int c = 0; c < NC; ++c) out.block<3, 1>(Eigen::Index(3) * c, j) = AWz.col(c);
        }
        return out;
    };
    for (int k = 0; k < n; ++k) {
        const double t = k * dt;
        const Eigen::Matrix<double, 12, 3> k1 = Ablk(Y) + Ct(t);
        const Eigen::Matrix<double, 12, 3> k2 = Ablk(Y + dt / 2 * k1) + Ct(t + dt / 2);
        const Eigen::Matrix<double, 12, 3> k3 = Ablk(Y + dt / 2 * k2) + Ct(t + dt / 2);
        const Eigen::Matrix<double, 12, 3> k4 = Ablk(Y + dt * k3) + Ct(t + dt);
        Y += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    check("anisotropic A: derived Phi_{rho,phi} = int Phi_A C(W) vs FD", (Y - Jfd).cwiseAbs().maxCoeff() /
          std::max(1e-12, Jfd.cwiseAbs().maxCoeff()), 1e-7);
    report("anisotropic A: |Phi_{rho,phi}| over 50 ms (zero when isotropic)", Jfd.cwiseAbs().maxCoeff());
}

// ---------------------------------------------------------------------------
// G. Relation to the quaternion MEKF (Kalman3D_Wave_CARMA3)
// ---------------------------------------------------------------------------

using Mekf = Kalman3D_Wave_CARMA3<double, 4, false, 4, false>;
// CARMA3 index of the invariant-layout index i (same coordinates e).
int mekf_index(int i) {
    if (i < 6 || i >= F::OFF_BA) return i;
    const int col = (i - 6) / 3, axis = (i - 6) % 3;
    return 6 + axis * 4 + col;
}
MX mekf_cov_in_tfg_layout(const Mekf& m) {
    MX P;
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j) P(i, j) = m.covariance_full()(mekf_index(i), mekf_index(j));
    return P;
}

template<class Fn>
void configure_both(Mekf& m, F& t, Fn&& extra) {
    m.set_vertical_carma(KAPPA, 0.6);
    m.set_horizontal_carma(KAPPA, 0.5, 0.5, 0.0);
    m.reset_aw_covariance_to_stationary();
    t.set_vertical_carma(KAPPA, 0.6);
    t.set_horizontal_carma(KAPPA, 0.5, 0.5, 0.0);
    t.reset_aw_covariance_to_stationary();
    extra(m, t);
}

void test_mekf_relation() {
    const V3 sa = V3::Constant(0.03), sg = V3::Constant(3e-4), sm = V3::Constant(0.8);
    for (int mode = 0; mode < 3; ++mode) {
        Mekf m(sa, sg, sm, 5e-4, 1e-6, 1e-8, 1.5, GRAV);
        F t(sa, sg, sm, 5e-4, 1e-6, 1e-8, 1.5, GRAV);
        t.set_reset_mode(mode == 0 ? F::ResetMode::Full : (mode == 1 ? F::ResetMode::AttitudeOnly : F::ResetMode::None));
        Synth sy;
        configure_both(m, t, [&](Mekf& a, F& b) {
            a.set_mag_world_ref(sy.Bw);
            b.set_mag_world_ref(sy.Bw);
            a.initialize_from_attitude(Eigen::Quaterniond(sy.R(0)), 0.035, 0.3);
            b.initialize_from_attitude(Eigen::Quaterniond(sy.R(0)), 0.035, 0.3);
            a.set_acc_bias_updates_enabled(true);
            b.set_acc_bias_updates_enabled(true);
        });
        std::mt19937 g(11);
        const char* mn = mode == 0 ? "Full" : (mode == 1 ? "AttitudeOnly" : "None");
        double d_first_pred = 0, d_first_upd_x = 0, d_first_upd_P = 0, worst_ang = 0, worst_p = 0, worst_P = 0;
        for (int k = 0; k < 200 * 60; ++k) {
            const Imu s = sy.sample(0.005, g, 1.0);
            sy.t += 0.005;
            m.time_update(s.gyr, 0.005);
            t.time_update(s.gyr, 0.005);
            if (k == 0) d_first_pred = rel_max(t.covariance_physical(), mekf_cov_in_tfg_layout(m));
            m.measurement_update_acc_only(s.acc);
            t.measurement_update_acc_only(s.acc);
            if (k == 0) {
                d_first_upd_x = std::max(lie::Log<double>(M3(t.state().R * m.quaternion_boat().toRotationMatrix().transpose())).norm(),
                                         (t.get_position() - m.get_position()).norm());
                d_first_upd_P = rel_max(t.covariance_physical(), mekf_cov_in_tfg_layout(m));
            }
            if (k % 8 == 0) {
                m.measurement_update_mag_only(s.mag);
                t.measurement_update_mag_only(s.mag);
            }
            worst_ang = std::max(worst_ang, lie::Log<double>(M3(t.state().R * m.quaternion_boat().toRotationMatrix().transpose())).norm());
            worst_p = std::max(worst_p, (t.get_position() - m.get_position()).norm());
            if (k % 200 == 0) worst_P = std::max(worst_P, rel_max(t.covariance_physical(), mekf_cov_in_tfg_layout(m)));
        }
        if (mode == 0) {
            check("MEKF relation: one prediction from equal seeds, T^-1 P T^-T vs CARMA3 P", d_first_pred, 1e-9);
            report("MEKF relation: first accel update, state difference (rad / m)", d_first_upd_x);
            report("MEKF relation: first accel update, covariance difference (corr units)", d_first_upd_P);
        }
        report(std::string("MEKF relation (") + mn + "): max attitude difference over 60 s (rad)", worst_ang);
        report(std::string("MEKF relation (") + mn + "): max position difference over 60 s (m)", worst_p);
        report(std::string("MEKF relation (") + mn + "): max covariance difference (corr units)", worst_P);
    }
}

// Yaw is unobservable without the magnetometer (gravity-preserving world
// symmetry).  In the invariant coordinates its direction is e_phi_z for every
// estimate, so the information e^T P^-1 e cannot grow from accelerometer
// updates or propagation; in CARMA3's coordinates the direction
// n(xhat) = T^-1(xhat) e_phi_z moves with the estimate.
void test_yaw_information() {
    const V3 sa = V3::Constant(0.03), sg = V3::Constant(3e-4), sm = V3::Constant(0.8);
    for (int mode = 0; mode < 3; ++mode) {
        Mekf m(sa, sg, sm, 5e-4, 1e-6, 1e-8, 1.5, GRAV);
        F t(sa, sg, sm, 5e-4, 1e-6, 1e-8, 1.5, GRAV);
        t.set_reset_mode(mode == 0 ? F::ResetMode::Full : (mode == 1 ? F::ResetMode::AttitudeOnly : F::ResetMode::None));
        Synth sy;
        sy.turns = true;
        configure_both(m, t, [&](Mekf& a, F& b) {
            a.initialize_from_attitude(Eigen::Quaterniond(sy.R(0)), 0.035, 0.5);
            b.initialize_from_attitude(Eigen::Quaterniond(sy.R(0)), 0.035, 0.5);
            a.set_acc_bias_updates_enabled(true);
            b.set_acc_bias_updates_enabled(true);
        });
        std::mt19937 g(5);
        TX ez = TX::Zero();
        ez(2) = 1.0;
        auto info_tfg = [&]() { return ez.dot(t.covariance_full().ldlt().solve(ez)); };
        auto info_mekf = [&]() {
            // n = T^-1(xhat) e_z at the MEKF's own estimate: dtheta = R^T e_z,
            // dw_i = [w_i]x e_z (horizontal latent from the twin TFG state).
            TX n = TX::Zero();
            const M3 R = m.quaternion_boat().toRotationMatrix();
            n.head<3>() = R.transpose() * V3::UnitZ();
            const V3 cols[3] = {m.get_position(), m.get_velocity(), m.get_world_accel()};
            for (int c = 0; c < 3; ++c) n.segment<3>(6 + 3 * c) = lie::skew<double>(cols[c]) * V3::UnitZ();
            n.segment<3>(6 + 9) = lie::skew<double>(V3(t.state().X.col(3))) * V3::UnitZ();
            return n.dot(mekf_cov_in_tfg_layout(m).ldlt().solve(n));
        };
        double max_rise_t = 0, max_rise_m = 0;
        double it = info_tfg(), im = info_mekf();
        const double it0 = it, im0 = im;
        for (int k = 0; k < 200 * 120; ++k) {
            const Imu s = sy.sample(0.005, g, 1.0);
            sy.t += 0.005;
            m.time_update(s.gyr, 0.005);
            t.time_update(s.gyr, 0.005);
            m.measurement_update_acc_only(s.acc);
            t.measurement_update_acc_only(s.acc);
            const double it1 = info_tfg(), im1 = info_mekf();
            max_rise_t = std::max(max_rise_t, (it1 - it) / it);
            max_rise_m = std::max(max_rise_m, (im1 - im) / im);
            it = it1;
            im = im1;
        }
        const char* mn = mode == 0 ? "Full" : (mode == 1 ? "AttitudeOnly" : "None");
        if (mode == 2) {
            check("yaw information e^T P^-1 e never rises, TFG reset None (rel. per step)", max_rise_t, 1e-9);
            report("yaw std (deg) after 120 s accel-only: TFG reset None", std::sqrt(t.covariance_full()(2, 2)) * 180 / M_PI);
        }
        report(std::string("yaw information: max relative rise per step, TFG reset ") + mn, max_rise_t);
        report(std::string("yaw information: final/initial, TFG reset ") + mn, it / it0);
        if (mode == 0) {
            report("yaw information: max relative rise per step, CARMA3 MEKF", max_rise_m);
            report("yaw information: final/initial, CARMA3 MEKF", im / im0);
            report("yaw std (deg) after 120 s accel-only: TFG", std::sqrt(t.covariance_full()(2, 2)) * 180 / M_PI);
            const M3 R = m.quaternion_boat().toRotationMatrix();
            const M3 Pw = R * m.covariance_full().block<3, 3>(0, 0) * R.transpose();
            report("yaw std (deg) after 120 s accel-only: CARMA3", std::sqrt(Pw(2, 2)) * 180 / M_PI);
        }
    }
}

// ---------------------------------------------------------------------------
// H. Behaviour checks
// ---------------------------------------------------------------------------

void test_behaviour() {
    // Held accelerometer bias: no measurement changes it, its cross-covariance
    // stays zero.
    {
        F f = make_filter();
        Synth sy;
        f.state() = make_state(sy.R(0), 0.0, V3::Zero(), V3(0.01, -0.02, 0.03));
        f.set_acc_bias_updates_enabled(false);
        std::mt19937 g(1);
        const V3 b0 = f.get_acc_bias();
        double cross = 0;
        for (int k = 0; k < 2000; ++k) {
            const Imu s = sy.sample(0.005, g, 1.0);
            sy.t += 0.005;
            f.time_update(s.gyr, 0.005);
            f.measurement_update_acc_only(s.acc);
            if (k % 8 == 0) f.measurement_update_mag_only(s.mag);
            for (int i = 0; i < F::OFF_BA; ++i)
                for (int j = 0; j < 3; ++j) cross = std::max(cross, std::abs(f.covariance_full()(i, F::OFF_BA + j)));
        }
        check("held accel bias: unchanged by 2000 updates", (f.get_acc_bias() - b0).norm(), 1e-15);
        check("held accel bias: cross-covariance stays zero", cross, 1e-15);
    }
    // CARMA3's noiseless three-axis sinusoid (carma3-test, same inputs and the
    // same 10 % bar on every axis), then a rotating version paired with CARMA3.
    {
        using Ff = Kalman3D_Wave_TFG_CARMA3<float, 4>;
        using Vf = Eigen::Vector3f;
        Ff f(Vf::Constant(0.03f), Vf::Constant(3e-5f), Vf::Constant(8.0f));
        const double w = 2.0 * M_PI * 0.30, dt = 0.005;
        const double fz = f.vertical_model().base_zero_crossing_hz();
        f.set_vertical_carma(static_cast<float>(0.30 / fz), 0.6f);
        f.set_horizontal_carma(static_cast<float>(0.30 / fz), 0.4f, 0.3f, 0.0f);
        f.reset_aw_covariance_to_stationary();
        f.set_acc_bias_updates_enabled(false);
        f.initialize_from_attitude(Eigen::Quaternionf::Identity(), 0.01f, 0.05f);
        const Vf mag_w(20.0f, -4.0f, 45.0f);
        f.set_mag_world_ref(mag_w);
        const double A[3] = {0.12, 0.08, 0.25}, ph[3] = {0.4, 1.3, 0.0};
        double se[3] = {}, sr[3] = {};
        for (int k = 0; k < 24000; ++k) {
            const double t = k * dt;
            Vf aw;
            for (int i = 0; i < 3; ++i) aw(i) = static_cast<float>(-A[i] * w * w * std::sin(w * t + ph[i]));
            f.time_update(Vf::Zero(), static_cast<float>(dt));
            f.measurement_update_acc_only(Vf(aw.x(), aw.y(), aw.z() - 9.80665f));
            if (k % 8 == 0) f.measurement_update_mag_only(mag_w);
            if (t > 60.0) {
                const Vf pe = f.get_position();
                for (int i = 0; i < 3; ++i) {
                    const double tr = A[i] * std::sin(w * t + ph[i]);
                    se[i] += std::pow(pe(i) - tr, 2);
                    sr[i] += tr * tr;
                }
            }
        }
        double worst = 0;
        for (int i = 0; i < 3; ++i) worst = std::max(worst, std::sqrt(se[i] / sr[i]));
        check("sinusoid (carma3-test input): worst-axis relative RMS error", worst, 0.10);
        check("sinusoid (carma3-test input): covariance PSD", -min_eig_corr(MX(f.covariance_full().cast<double>())), 1e-4);
    }
    {
        // Both filters constructed and configured identically (same noise
        // densities, prior, bias model and seeds).
        F f(V3::Constant(0.03), V3::Constant(3e-4), V3::Constant(0.8), 5e-4, 1e-6, 1e-8, 1.5, GRAV);
        Mekf m(V3::Constant(0.03), V3::Constant(3e-4), V3::Constant(0.8), 5e-4, 1e-6, 1e-8, 1.5, GRAV);
        configure_both(m, f, [](Mekf& a, F& b) {
            a.set_Q_bacc_rw(V3::Constant(5e-4));
            b.set_Q_bacc_rw(V3::Constant(5e-4));
            a.set_acc_bias_time_constant(800.0);
            b.set_acc_bias_time_constant(800.0);
        });
        Synth sy;
        f.set_mag_world_ref(sy.Bw);
        m.set_mag_world_ref(sy.Bw);
        f.initialize_from_attitude(Eigen::Quaterniond(sy.R(0)), 0.02, 0.05);
        m.initialize_from_attitude(Eigen::Quaterniond(sy.R(0)), 0.02, 0.05);
        std::mt19937 g(2);
        double se = 0, sm = 0, ss = 0;
        for (int k = 0; k < 200 * 200; ++k) {
            const Imu s = sy.sample(0.005, g, 0.0);
            sy.t += 0.005;
            f.time_update(s.gyr, 0.005);
            m.time_update(s.gyr, 0.005);
            f.measurement_update_acc_only(s.acc);
            m.measurement_update_acc_only(s.acc);
            if (k % 8 == 0) {
                f.measurement_update_mag_only(s.mag);
                m.measurement_update_mag_only(s.mag);
            }
            if (sy.t > 100) {
                const V3 p(0.3 * std::sin(2.0 * sy.t), 0.2 * std::sin(1.7 * sy.t + 0.5), 0.4 * std::sin(1.95 * sy.t + 1.0));
                se += (f.get_position() - p).squaredNorm();
                sm += (m.get_position() - p).squaredNorm();
                ss += p.squaredNorm();
            }
        }
        report("rotating sinusoid: TFG-CARMA3 relative displacement RMS error", std::sqrt(se / ss));
        report("rotating sinusoid: CARMA3 (MEKF) relative displacement RMS error", std::sqrt(sm / ss));
        check("rotating sinusoid: |TFG/MEKF error ratio - 1| (first-order equivalence)",
              std::abs(std::sqrt(se / sm) - 1.0), 0.01);
        check("rotating sinusoid: TFG covariance PSD", -min_eig_corr(f.covariance_full()), 1e-9);
    }
    // Physical write mode: an absolute attitude write leaves CARMA3's
    // covariance (e coordinates) unchanged, as CARMA3's set_quaternion_boat.
    {
        F f = make_filter();
        f.state() = make_state(rot_euler(0.1, -0.1, 0.5), 0.5, V3::Zero(), V3::Zero());
        f.covariance_full() = scaled_cov(0.01, 1e-3, 0.05, 0.01);
        const MX Pe0 = f.covariance_physical();
        f.set_quaternion_boat(Eigen::Quaterniond(rot_euler(0.1, -0.1, 2.0)));
        check("physical write: attitude write keeps the e covariance", rel_max(f.covariance_physical(), Pe0), 1e-12);
        F n = make_filter();
        n.set_covariance_write_mode(F::CovarianceWriteMode::Native);
        n.state() = make_state(rot_euler(0.1, -0.1, 0.5), 0.5, V3::Zero(), V3::Zero());
        n.covariance_full() = scaled_cov(0.01, 1e-3, 0.05, 0.01);
        const MX Pn0 = n.covariance_full();
        n.set_quaternion_boat(Eigen::Quaterniond(rot_euler(0.1, -0.1, 2.0)));
        check("native write: attitude write keeps the invariant covariance", rel_max(n.covariance_full(), Pn0), 1e-15);
    }
    // The wrapper on a still, noisy IMU: finite and bounded on every axis.
    {
        ocean_imu::carma::SeaStateFusion_TFG_CARMA3<TrackerType::KALMANF> fusion;
        ocean_imu::carma::SeaStateFusionConfig_CARMA3 cfg;
        cfg.with_mag = true;
        cfg.sigma_a = Eigen::Vector3f::Constant(2.8f * 1.51e-3f * 9.80665f * 0.71f);
        cfg.sigma_g = Eigen::Vector3f::Constant(2.0f * 0.00157f * 0.01f);
        cfg.sigma_m = Eigen::Vector3f::Constant(1.2f * 0.8f * 8.0f);
        fusion.begin(cfg);
        fusion.raw().enableTuner(true);
        fusion.raw().enableClamp(true);
        std::mt19937 g(9);
        std::normal_distribution<float> n(0.0f, 1.0f);
        const Eigen::Vector3f Bw(18.0f, 1.5f, 43.0f);
        float worst = 0.0f;
        bool finite = true;
        for (int k = 0; k < 200 * 400; ++k) {
            const Eigen::Vector3f acc(0.004f + 0.015f * n(g), -0.003f + 0.015f * n(g), -9.80665f + 0.006f + 0.015f * n(g));
            const Eigen::Vector3f gyr(3e-4f * n(g), 3e-4f * n(g), 3e-4f * n(g));
            fusion.update(0.005f, gyr, acc, 35.0f);
            if (k % 8 == 0) fusion.updateMag(Bw + Eigen::Vector3f(0.8f * n(g), 0.8f * n(g), 0.8f * n(g)));
            const Eigen::Vector3f p = fusion.raw().mekf().get_position();
            finite = finite && p.allFinite() && fusion.raw().mekf().covariance_full().allFinite();
            if (fusion.isLive() && k > 200 * 240) worst = std::max(worst, p.cwiseAbs().maxCoeff());
        }
        check("wrapper: still IMU went live", fusion.isLive() ? 0.0 : 1.0, 0.0);
        check("wrapper: still IMU finite", finite ? 0.0 : 1.0, 0.0);
        check("wrapper: still IMU |p| < 0.5 m after 240 s (as carma3-test)", worst, 0.5);
        check("wrapper: common CARMA scale (no rejected anisotropic configuration)",
              double(fusion.raw().mekf().rejected_anisotropic_configurations()), 0.0);
    }
}

}  // namespace

int main() {
    test_group();
    test_reset();
    test_measurement_jacobians();
    test_update_against_reference();
    test_transition();
    test_quadrature_collapse();
    test_semigroup_and_linearization();
    test_structured_propagation();
    test_monte_carlo_propagation();
    test_stationary_prior();
    test_process_spectrum();
    test_float_vs_double();
    test_equivariance();
    test_anisotropic_dynamics();
    test_mekf_relation();
    test_yaw_information();
    test_behaviour();
    std::printf("tfg_carma3-test: %d/%d checks passed\n", g_checks - g_failures, g_checks);
    return g_failures == 0 ? 0 : 1;
}
