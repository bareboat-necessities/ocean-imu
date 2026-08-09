/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Independent hardening tests for the two-frame covariance machinery.

    The important rule in this file is that the expected covariance is never
    defined as "whatever the implementation already did". In particular, the
    finite-retraction reset is checked against the defining finite-difference
    map

        delta -> Log_G( Exp_G(-a) Exp_G(a+delta) ),

    and the rotating gyro-bias process noise is checked against a high-rate
    numerical stochastic integral.
*/

#define EIGEN_NON_ARDUINO

#include "kalman_tfg/Kalman3D_Wave_TFG.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>

namespace {

using T = double;
using L1 = ocean_imu::kalman::Kalman3D_Wave_TFG<T, true, true, false>;
using L2 = ocean_imu::kalman::Kalman3D_Wave_TFG<T, true, true, true>;
using Vector3 = Eigen::Matrix<T,3,1>;
using Matrix3 = Eigen::Matrix<T,3,3>;
static_assert(L1::NX == L2::NX);
constexpr int NX = L2::NX;
using MatrixNX = Eigen::Matrix<T,NX,NX>;
using Tangent = Eigen::Matrix<T,NX,1>;

int failures = 0;
std::mt19937 rng(20260809u);

bool check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
    return condition;
}

Vector3 rand_vec(T scale) {
    std::uniform_real_distribution<T> d(-scale, scale);
    return Vector3(d(rng), d(rng), d(rng));
}

MatrixNX seed_covariance(T scale = T(1)) {
    MatrixNX A;
    std::normal_distribution<T> n(T(0), T(1));
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j) A(i,j) = n(rng);
    MatrixNX P = scale * (A * A.transpose());
    P += MatrixNX::Identity() * (scale * T(0.05));
    return P;
}

MatrixNX bias_map(const Matrix3& R) {
    MatrixNX M = MatrixNX::Identity();
    M.block<3,3>(L2::OFF_BG, L2::OFF_BG) = R;
    M.block<3,3>(L2::OFF_BA, L2::OFF_BA) = R;
    return M;
}

template<typename F>
void configure(F& f) {
    f.set_aw_time_constant(T(1.7));
    f.set_aw_stationary_std(Vector3(T(0.6), T(0.55), T(0.32)));
    f.set_gyro_noise_density_rad_sqrt_s(T(0.004));
    f.set_Q_bgyro_rw(Vector3(T(1e-6), T(2e-6), T(1.5e-6)));
    f.set_Q_bacc_rw(Vector3(T(5e-4), T(5.5e-4), T(4.5e-4)));
    f.set_acc_bias_time_constant(T(1200));
    f.set_Racc_std(Vector3(T(0.4), T(0.35), T(0.5)));
    f.set_Rmag_std(Vector3(T(0.1), T(0.12), T(0.09)));
    f.set_RS_noise(Vector3(T(1.5), T(1.2), T(0.8)));
    f.set_magnetic_reference_world(Vector3(T(0.21), T(0.02), T(0.43)));
    f.initialize_from_truth(
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(
            T(0.83), Vector3(T(0.31), T(-0.62), T(0.72)).normalized())),
        Vector3(T(0.5), T(-0.8), T(1.1)),
        Vector3(T(-1.7), T(2.4), T(0.9)),
        Vector3(T(3.3), T(-2.1), T(1.4)),
        Vector3(T(0.3), T(-0.5), T(0.7)),
        Vector3(T(0.012), T(-0.008), T(0.004)),
        Vector3(T(0.02), T(0.05), T(-0.01)));
}

// ---------------------------------------------------------------------------
// Level-2 group/retraction consistency and L1 <-> L2 propagation equivalence.
// ---------------------------------------------------------------------------

void test_inject_error_roundtrip() {
    for (int trial = 0; trial < 30; ++trial) {
        L2 base; configure(base);
        Tangent xi = Tangent::Zero();
        xi.segment<3>(L2::OFF_PHI) = rand_vec(T(0.7));
        xi.segment<3>(L2::OFF_BG) = rand_vec(T(0.04));
        for (int c = 0; c < 4; ++c)
            xi.segment<3>(L2::OFF_V + 3*c) = rand_vec(T(1.5));
        xi.segment<3>(L2::OFF_BA) = rand_vec(T(0.08));
        L2 perturbed = base;
        perturbed.inject(xi);
        const T err = (perturbed.error_from(base) - xi).cwiseAbs().maxCoeff();
        check(err < T(2e-12), "error_from(inject(xi)) is not an exact Level-2 round trip");
    }
}

void test_levels_remain_equivalent_under_propagation() {
    const Vector3 rates[] = {
        Vector3::Zero(), Vector3(T(0.4), T(-0.7), T(0.5)),
        Vector3(T(2.0), T(-1.5), T(3.0))
    };
    const T steps[] = {T(0.001), T(0.005), T(0.03), T(0.1)};

    for (const Vector3& gyro : rates) {
        for (T h : steps) {
            L1 a; configure(a);
            L2 b; configure(b);
            const MatrixNX P0 = seed_covariance(T(0.02));
            const Matrix3 R0 = a.R_bw();
            a.covariance_full() = P0;
            b.covariance_full() = bias_map(R0) * P0 * bias_map(R0).transpose();

            a.time_update(gyro, h);
            b.time_update(gyro, h);

            check((a.R_bw()-b.R_bw()).cwiseAbs().maxCoeff() == T(0) &&
                  (a.state().X-b.state().X).cwiseAbs().maxCoeff() == T(0) &&
                  (a.state().B-b.state().B).cwiseAbs().maxCoeff() == T(0),
                  "Level 1 and Level 2 propagated different nominal states");

            const MatrixNX M = bias_map(a.R_bw());
            const MatrixNX expected = M * a.covariance_full() * M.transpose();
            const T scale = std::max(T(1), expected.cwiseAbs().maxCoeff());
            const T err = (b.covariance_full()-expected).cwiseAbs().maxCoeff();
            if (!check(err < T(3e-11)*scale,
                       "Level-2 covariance lost exact coordinate equivalence"))
                std::cerr << "  h=" << h << " |gyro|=" << gyro.norm()
                          << " max error=" << err << '\n';
        }
    }
}

// ---------------------------------------------------------------------------
// Full-group covariance reset oracle.
// ---------------------------------------------------------------------------

MatrixNX finite_difference_reset_jacobian(const Tangent& a) {
    MatrixNX J;
    const auto left = L2::Group::Exp_G(Tangent(-a));
    const T eps = T(2e-7);
    for (int j = 0; j < NX; ++j) {
        Tangent e = Tangent::Zero();
        e(j) = eps;
        const Tangent fp = L2::Group::Log_G(left * L2::Group::Exp_G(Tangent(a + e)));
        const Tangent fm = L2::Group::Log_G(left * L2::Group::Exp_G(Tangent(a - e)));
        J.col(j) = (fp - fm) / (T(2)*eps);
    }
    return J;
}

void test_reset_jacobian_against_definition() {
    L2 f; configure(f);
    for (int trial = 0; trial < 8; ++trial) {
        Tangent a = Tangent::Zero();
        a.segment<3>(L2::OFF_PHI) = rand_vec(T(0.55));
        a.segment<3>(L2::OFF_BG) = rand_vec(T(0.025));
        for (int c = 0; c < 4; ++c)
            a.segment<3>(L2::OFF_V + 3*c) = rand_vec(T(0.7));
        a.segment<3>(L2::OFF_BA) = rand_vec(T(0.04));

        MatrixNX analytic;
        f.build_reset_jacobian(a, analytic);
        const MatrixNX numeric = finite_difference_reset_jacobian(a);
        const T err = (analytic-numeric).cwiseAbs().maxCoeff();
        if (!check(err < T(2e-8),
                   "full-group reset Jacobian disagrees with finite differences"))
            std::cerr << "  max reset Jacobian error=" << err << '\n';
    }
}

void test_update_covariance_uses_full_reset() {
    L2 f; configure(f);
    f.covariance_full() = seed_covariance(T(0.008));

    Eigen::Matrix<T,3,NX> H;
    std::normal_distribution<T> n(T(0), T(0.35));
    for (int i = 0; i < H.rows(); ++i)
        for (int j = 0; j < H.cols(); ++j) H(i,j) = n(rng);
    const Vector3 r(T(0.38), T(-0.27), T(0.31));
    const Matrix3 Rw = (Vector3(T(0.18), T(0.22), T(0.16)).array().square()).matrix().asDiagonal();

    const MatrixNX P0 = f.covariance_full();
    const Eigen::Matrix<T,NX,3> PHt = P0 * H.transpose();
    const Matrix3 S = H * PHt + Rw;
    const Eigen::Matrix<T,NX,3> K = PHt * S.inverse();
    const Tangent a = K * r;
    const MatrixNX IKH = MatrixNX::Identity() - K*H;
    const MatrixNX Pj = IKH*P0*IKH.transpose() + K*Rw*K.transpose();
    const MatrixNX Jfd = finite_difference_reset_jacobian(a);
    const MatrixNX expected = Jfd * Pj * Jfd.transpose();

    L2::MeasDiag3 diag;
    check(f.apply_update3(r, H, Rw, diag), "synthetic update was rejected");
    const T scale = std::max(T(1), expected.cwiseAbs().maxCoeff());
    const T err = (f.covariance_full()-expected).cwiseAbs().maxCoeff();
    if (!check(err < T(3e-9)*scale,
               "post-update covariance does not match the independent full reset oracle"))
        std::cerr << "  max covariance error=" << err << '\n';

    const T reset_effect = (expected-Pj).cwiseAbs().maxCoeff();
    check(reset_effect > T(1e-6),
          "test correction is too small to distinguish a real reset from Joseph-only covariance");
}

// ---------------------------------------------------------------------------
// Accelerometer-bias OU parity with current OU-III.
// ---------------------------------------------------------------------------

void test_accel_bias_ou_mean_and_covariance() {
    L1 f;
    f.initialize_identity();
    f.set_linear_block_enabled(false);
    f.set_gyro_noise_density_rad_sqrt_s(T(0));
    f.set_Q_bgyro_rw(Vector3::Zero());
    f.set_aw_stationary_std(Vector3::Zero());

    const T tau = T(1000);
    const T dt = T(0.01);
    const int steps = 200;
    const Vector3 drive(T(5e-4), T(7e-4), T(4e-4));
    const Vector3 b0(T(0.08), T(-0.05), T(0.03));
    const Matrix3 P0 = (Vector3(T(0.02), T(0.03), T(0.025)).array().square()).matrix().asDiagonal();

    f.set_acc_bias_time_constant(tau);
    f.set_Q_bacc_rw(drive);
    f.set_initial_acc_bias(b0);
    f.covariance_full().setZero();
    f.covariance_full().block<3,3>(L1::OFF_BA,L1::OFF_BA) = P0;

    for (int k = 0; k < steps; ++k) f.time_update(Vector3::Zero(), dt);

    const T elapsed = dt*T(steps);
    const T alpha = std::exp(-elapsed/tau);
    const Matrix3 Qc = drive.array().square().matrix().asDiagonal();
    const Matrix3 Pexpected = alpha*alpha*P0 +
        (T(0.5)*tau*(T(1)-std::exp(T(-2)*elapsed/tau))) * Qc;

    check((f.get_acc_bias()-alpha*b0).cwiseAbs().maxCoeff() < T(2e-13),
          "accelerometer-bias OU mean is not exp(-t/tau) b0");
    const T perr = (f.covariance_full().block<3,3>(L1::OFF_BA,L1::OFF_BA)-Pexpected)
                       .cwiseAbs().maxCoeff();
    if (!check(perr < T(2e-12), "accelerometer-bias OU covariance is not exact"))
        std::cerr << "  OU covariance error=" << perr << '\n';

    // Startup freeze must freeze the physical body-bias state, not keep
    // decaying it while measurements are prohibited from correcting it.
    L1 frozen;
    frozen.initialize_identity();
    frozen.set_initial_acc_bias(b0);
    frozen.covariance_full().setZero();
    frozen.covariance_full().block<3,3>(L1::OFF_BA,L1::OFF_BA) = P0;
    frozen.set_acc_bias_time_constant(T(2));
    frozen.set_Q_bacc_rw(Vector3::Constant(T(0.01)));
    frozen.set_acc_bias_updates_enabled(false);
    const Matrix3 Pfrozen = frozen.covariance_full().block<3,3>(L1::OFF_BA,L1::OFF_BA);
    for (int k = 0; k < 100; ++k)
        frozen.time_update(Vector3(T(0.8),T(-0.4),T(0.3)), T(0.01));
    check((frozen.get_acc_bias()-b0).cwiseAbs().maxCoeff() == T(0),
          "disabled accelerometer-bias state did not freeze its mean");
    check((frozen.covariance_full().block<3,3>(L1::OFF_BA,L1::OFF_BA)-Pfrozen)
              .cwiseAbs().maxCoeff() < T(2e-14),
          "disabled accelerometer-bias state did not freeze its body-frame marginal");
}

// ---------------------------------------------------------------------------
// Rotating gyro-bias Qd oracle.
// ---------------------------------------------------------------------------

void test_gyro_bias_process_noise_against_high_rate_oracle() {
    L1 f;
    f.initialize_from_truth(
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(
            T(0.91), Vector3(T(0.2),T(-0.6),T(0.7)).normalized())),
        Vector3::Zero(), Vector3::Zero(), Vector3::Zero(), Vector3::Zero());
    f.set_linear_block_enabled(false);
    f.set_gyro_noise_density_rad_sqrt_s(T(0));
    const Vector3 qdiag(T(3e-5), T(8e-5), T(5e-5));
    f.set_Q_bgyro_rw(qdiag);
    f.set_Q_bacc_rw(Vector3::Zero());
    f.set_acc_bias_updates_enabled(false);

    const T h = T(0.12);
    const Vector3 gyro(T(3.0), T(-2.0), T(4.0));
    MatrixNX Qd;
    f.build_process_noise(h, gyro, Qd);

    const Matrix3 Qb = qdiag.asDiagonal();
    const Matrix3 R0 = f.R_bw();
    const Vector3 omega = gyro - f.gyroscope_bias();
    constexpr int N = 12000;
    const T ds = h/T(N);
    Matrix3 tail = Matrix3::Zero();
    Matrix3 Qpp = Matrix3::Zero();
    Matrix3 Qpb = Matrix3::Zero();

    for (int k = N-1; k >= 0; --k) {
        const T smid = (T(k)+T(0.5))*ds;
        const Matrix3 Rmid = R0 * ocean_imu::lie::Exp<T>(Vector3(omega*smid));
        const Matrix3 Cmid = tail + T(0.5)*ds*Rmid;
        Qpp.noalias() += ds * (Cmid*Qb*Cmid.transpose());
        Qpb.noalias() -= ds * (Cmid*Qb);
        tail += ds*Rmid;
    }

    const T epp = (Qd.block<3,3>(L1::OFF_PHI,L1::OFF_PHI)-Qpp).cwiseAbs().maxCoeff();
    const T epb = (Qd.block<3,3>(L1::OFF_PHI,L1::OFF_BG)-Qpb).cwiseAbs().maxCoeff();
    const T ebb = (Qd.block<3,3>(L1::OFF_BG,L1::OFF_BG)-Qb*h).cwiseAbs().maxCoeff();
    if (!check(epp < T(2e-10), "Q_phi,phi gyro-bias integral misses high-rate oracle"))
        std::cerr << "  Qpp error=" << epp << '\n';
    if (!check(epb < T(2e-10), "Q_phi,bg gyro-bias integral misses high-rate oracle"))
        std::cerr << "  Qpb error=" << epb << '\n';
    check(ebb < T(2e-13), "Q_bg,bg is not Q_bg*h at Level 1");

    // The old constant-rotation formula must be measurably wrong for this case,
    // otherwise the oracle would not discriminate the issue reported in review.
    const Matrix3 old = -(R0*Qb)*(h*h/T(2));
    check((old-Qpb).cwiseAbs().maxCoeff() > T(1e-7),
          "stochastic-oracle test rate is too small to distinguish old h^2/2 formula");
}

// ---------------------------------------------------------------------------
// Complete world-yaw gauge transformation.
// ---------------------------------------------------------------------------

void test_world_yaw_gauge_rotates_state_and_covariance() {
    L2 f; configure(f);
    const MatrixNX P0 = seed_covariance(T(0.01));
    f.covariance_full() = P0;
    const auto state0 = f.state();
    const Matrix3 Sigma0 = f.aw_stationary_cov();
    const Vector3 B0 = f.magnetic_reference_world();

    const T yaw = T(-0.73);
    const Matrix3 Q = ocean_imu::lie::Exp<T>(Vector3(T(0),T(0),yaw));
    MatrixNX G = MatrixNX::Identity();
    G.block<3,3>(L2::OFF_PHI,L2::OFF_PHI) = Q;
    for (int c = 0; c < 4; ++c)
        G.block<3,3>(L2::OFF_V+3*c,L2::OFF_V+3*c) = Q;
    G.block<3,3>(L2::OFF_BG,L2::OFF_BG) = Q;
    G.block<3,3>(L2::OFF_BA,L2::OFF_BA) = Q;
    const MatrixNX Pexpected = G*P0*G.transpose();

    f.apply_world_yaw_gauge(yaw);

    check((f.R_bw()-Q*state0.R).cwiseAbs().maxCoeff() < T(2e-15),
          "world yaw gauge did not rotate R");
    check((f.state().X-Q*state0.X).cwiseAbs().maxCoeff() < T(2e-15),
          "world yaw gauge did not rotate all world columns");
    check((f.state().B-state0.B).cwiseAbs().maxCoeff() == T(0),
          "world yaw gauge changed body-frame bias states");
    check((f.covariance_full()-Pexpected).cwiseAbs().maxCoeff() < T(3e-14),
          "world yaw gauge did not congruence-transform the full covariance");
    check((f.aw_stationary_cov()-Q*Sigma0*Q.transpose()).cwiseAbs().maxCoeff() < T(2e-15),
          "world yaw gauge did not rotate Sigma_aw");
    check((f.magnetic_reference_world()-Q*B0).cwiseAbs().maxCoeff() < T(2e-15),
          "world yaw gauge did not rotate an existing magnetic reference");
}

} // namespace

int main() {
    test_inject_error_roundtrip();
    test_levels_remain_equivalent_under_propagation();
    test_reset_jacobian_against_definition();
    test_update_covariance_uses_full_reset();
    test_accel_bias_ou_mean_and_covariance();
    test_gyro_bias_process_noise_against_high_rate_oracle();
    test_world_yaw_gauge_rotates_state_and_covariance();

    if (failures != 0) {
        std::cerr << failures << " covariance/math hardening check(s) failed\n";
        return 1;
    }
    std::cout << "covariance_transport-test: all checks passed\n";
    return 0;
}
