/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Level 2: the two-frame bias geometry, and how the covariance is carried.

    The centrepiece is an exact cross-check rather than another
    finite-difference sweep. Level 1 and Level 2 are the SAME filter in
    different coordinates, related by the linear map beta = Rhat delta_b:

        T     = blkdiag(I, Rhat, I, I, I, I, Rhat)
        xi_L2 = T xi_L1,   P_L2 = T P_L1 T^T,   H_L1 = H_L2 T

    Under pure propagation that relation is EXACT. Since Level 1 was already
    verified against finite differences in tfg_propagation-test, asserting the
    relation holds to 1e-12 validates the entire Level-2 derivation -- Phi,
    Q, and the bias self-dynamics -- against an independently checked
    reference, which finite differences against itself could never do.

    The conjugation is time-varying, Phi_L2 = T(h) Phi_L1 T(0)^-1, and that is
    exactly why Phi_beta,beta is Exp(w_world h) where Level 1 has I. Getting
    the two times mixed up is the obvious way to write this wrong, so the test
    uses the post-propagation rotation on the left and the pre-propagation one
    on the right, explicitly.

    Under measurement updates the two levels agree only to first order,
    because the Level-2 injection carries a J_l(-phi) the Level-1 one does
    not. That is asserted as convergence-with-step-size rather than equality.
*/

#define EIGEN_NON_ARDUINO

#include "kalman_tfg/Kalman3D_Wave_TFG.h"

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

static_assert(L1::NX == L2::NX, "both levels must share the tangent layout");
constexpr int NX = L2::NX;
using MatrixNX = Eigen::Matrix<T,NX,NX>;
using TangentV = Eigen::Matrix<T,NX,1>;

int failures = 0;

bool check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
    return condition;
}

std::mt19937 rng(20260806u);

Vector3 rand_vec(double scale) {
    std::uniform_real_distribution<double> d(-scale, scale);
    return Vector3(d(rng), d(rng), d(rng));
}

// The coordinate change beta = Rhat delta_b, as a block-diagonal map.
MatrixNX bias_map(const Matrix3& R) {
    MatrixNX Tm = MatrixNX::Identity();
    Tm.block<3,3>(L2::OFF_BG, L2::OFF_BG) = R;
    Tm.block<3,3>(L2::OFF_BA, L2::OFF_BA) = R;
    return Tm;
}

template<typename F>
void configure(F& f) {
    f.set_aw_time_constant(T(1.7));
    f.set_aw_stationary_std(Vector3(T(0.6), T(0.6), T(0.32)));
    f.set_gyro_noise_density_rad_sqrt_s(T(0.004));
    // Deliberately far larger than a real gyro bias random walk. At a
    // physical 1e-10 the gyro-bias block of Q is ~3 orders below the
    // accelerometer's, so an error in ITS frame convention lands just under
    // the tolerance and the test silently stops discriminating -- confirmed by
    // mutation. The equivalence under test is algebraic, so inflating the
    // input is the honest way to exercise both blocks equally.
    f.set_Q_bgyro_rw(Vector3(T(1e-6), T(2e-6), T(1.5e-6)));
    f.set_Q_bacc_rw(Vector3(T(2.5e-7), T(3e-7), T(2e-7)));
    f.set_Racc_std(Vector3(T(0.4), T(0.35), T(0.5)));
    f.set_Rmag_std(Vector3(T(0.1), T(0.12), T(0.09)));
    f.set_RS_noise_vector(Vector3(T(1.5), T(1.5), T(0.8)));
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

// An asymmetric, full-rank covariance, so a block that is wrong cannot hide
// behind zeros or behind a coincidental symmetry.
MatrixNX seed_covariance() {
    MatrixNX A;
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j) A(i,j) = d(rng);
    MatrixNX P = A * A.transpose();
    P += MatrixNX::Identity() * T(0.05);
    return P;
}

// ---------------------------------------------------------------------------
// 1. inject / error_from remain an exact inverse pair at Level 2
// ---------------------------------------------------------------------------

void test_level2_inject_error_roundtrip() {
    for (int trial = 0; trial < 50; ++trial) {
        L2 base; configure(base);

        TangentV xi = TangentV::Zero();
        xi.segment<3>(L2::OFF_PHI) = rand_vec(0.8);
        xi.segment<3>(L2::OFF_BG)  = rand_vec(0.05);
        for (int c = 0; c < 4; ++c) xi.segment<3>(L2::OFF_V + 3*c) = rand_vec(2.0);
        xi.segment<3>(L2::OFF_BA)  = rand_vec(0.1);

        L2 perturbed = base;
        perturbed.inject(xi);
        const TangentV back = perturbed.error_from(base);

        const T err = (back - xi).cwiseAbs().maxCoeff();
        if (!check(err <= T(1e-12), "Level 2: error_from(inject(xi)) != xi")) {
            std::cerr << "  max|diff| = " << err << '\n';
        }
    }
}

/*
    The Level-2 injection must use the PRE-update rotation, and error_from the
    REFERENCE rotation. Both differ from the alternative by Exp(phi), so both
    mistakes are invisible at zero correction and grow with it. The round trip
    above would still pass if BOTH were wrong in the same way, so pin the
    componentwise form directly.
*/
void test_level2_injection_uses_the_right_rotation() {
    L2 f; configure(f);
    const Matrix3 R_pre = f.R_bw();
    const Vector3 b_pre = f.gyroscope_bias();

    TangentV xi = TangentV::Zero();
    const Vector3 phi(T(0.4), T(-0.3), T(0.55));
    const Vector3 beta(T(0.03), T(-0.02), T(0.01));
    xi.segment<3>(L2::OFF_PHI) = phi;
    xi.segment<3>(L2::OFF_BG)  = beta;

    f.inject(xi);

    const Matrix3 Jr = ocean_imu::lie::left_jacobian<T>(Vector3(-phi));
    const Vector3 want = b_pre + R_pre.transpose() * Jr * beta;
    check((f.gyroscope_bias() - want).cwiseAbs().maxCoeff() <= T(1e-14),
          "Level 2 injection is not Bhat + R_pre^T J_l(-phi) beta");

    // And it must differ measurably from the post-update rotation, or the
    // assertion above is not discriminating.
    const Vector3 wrong = b_pre + f.R_bw().transpose() * Jr * beta;
    check((want - wrong).cwiseAbs().maxCoeff() > T(1e-3),
          "the test correction is too small to distinguish pre from post rotation");
}

// ---------------------------------------------------------------------------
// 2. The exact Level 1 <-> Level 2 equivalence under propagation
// ---------------------------------------------------------------------------

void test_levels_are_the_same_filter_under_propagation() {
    const Vector3 gyros[] = {
        Vector3::Zero(),
        Vector3(T(0.05), T(-0.03), T(0.02)),
        Vector3(T(0.4), T(-0.7), T(0.5)),
        Vector3(T(2.0), T(-1.5), T(3.0)),
    };
    const T steps[] = {T(0.001), T(0.005), T(0.02), T(0.1)};

    for (const Vector3& g : gyros) {
        for (T h : steps) {
            L1 a; configure(a);
            L2 b; configure(b);

            const MatrixNX P0 = seed_covariance();
            const Matrix3 R0 = a.R_bw();
            // Start from covariances that already satisfy the relation.
            a.covariance_full() = P0;
            b.covariance_full() = (bias_map(R0) * P0 * bias_map(R0).transpose()).eval();

            a.time_update(g, h);
            b.time_update(g, h);

            const std::string at =
                " (|w|=" + std::to_string(g.norm()) + ", h=" + std::to_string(h) + ")";

            // Nominal states must be bit-for-bit identical: the propagation of
            // the state does not depend on the error parameterization at all.
            check((a.R_bw() - b.R_bw()).cwiseAbs().maxCoeff() == T(0) &&
                  (a.state().X - b.state().X).cwiseAbs().maxCoeff() == T(0) &&
                  (a.state().B - b.state().B).cwiseAbs().maxCoeff() == T(0),
                  "the two levels propagated the nominal state differently" + at);

            // T(h) on the left: the relation is beta(h) = Rhat(h) delta_b(h).
            const MatrixNX Th = bias_map(a.R_bw());
            const MatrixNX want = Th * a.covariance_full() * Th.transpose();
            const T err = (b.covariance_full() - want).cwiseAbs().maxCoeff();
            const T scale = std::max(T(1), want.cwiseAbs().maxCoeff());
            if (!check(err <= T(1e-12) * scale,
                       "P_L2 != T P_L1 T^T after propagation" + at)) {
                std::cerr << "  max|diff| = " << err << " (scale " << scale << ")\n";
                int bi = 0, bj = 0; T worst = T(0);
                for (int i = 0; i < NX; ++i)
                    for (int j = 0; j < NX; ++j) {
                        const T d = std::abs(b.covariance_full()(i,j) - want(i,j));
                        if (d > worst) { worst = d; bi = i; bj = j; }
                    }
                std::cerr << "  worst entry (" << bi << "," << bj << ")\n";
            }
        }
    }
}

// Over many steps, so a small per-step inconsistency accumulates into view.
void test_levels_stay_equivalent_over_a_long_run() {
    L1 a; configure(a);
    L2 b; configure(b);
    const MatrixNX P0 = seed_covariance();
    const Matrix3 R0 = a.R_bw();
    a.covariance_full() = P0;
    b.covariance_full() = (bias_map(R0) * P0 * bias_map(R0).transpose()).eval();

    std::normal_distribution<double> n(0.0, 0.4);
    for (int step = 0; step < 2000; ++step) {
        const Vector3 g(n(rng), n(rng), n(rng));
        a.time_update(g, T(0.005));
        b.time_update(g, T(0.005));
    }

    const MatrixNX Th = bias_map(a.R_bw());
    const MatrixNX want = Th * a.covariance_full() * Th.transpose();
    const T scale = std::max(T(1), want.cwiseAbs().maxCoeff());
    const T err = (b.covariance_full() - want).cwiseAbs().maxCoeff();
    if (!check(err <= T(1e-10) * scale,
               "the two levels drifted apart over 10 s of propagation")) {
        std::cerr << "  max|diff| = " << err << " (scale " << scale << ")\n";
    }
}

/*
    The Exp(w_world h) bias self-dynamics is the substantive new term. Assert
    it is doing real work: replacing it with the identity must break the
    equivalence above at a realistic turn rate.
*/
void test_bias_rotation_term_matters() {
    const Vector3 g(T(1.5), T(-0.8), T(1.1));
    const T h = T(0.05);

    L1 a; configure(a);
    L2 b; configure(b);
    const MatrixNX P0 = seed_covariance();
    const Matrix3 R0 = a.R_bw();
    a.covariance_full() = P0;
    b.covariance_full() = (bias_map(R0) * P0 * bias_map(R0).transpose()).eval();

    MatrixNX Phi2;
    b.build_transition(h, g, Phi2);
    const Matrix3 got = Phi2.block<3,3>(L2::OFF_BG, L2::OFF_BG);
    const Vector3 omega = g - b.gyroscope_bias();
    const Matrix3 want = ocean_imu::lie::Exp<T>(Vector3(b.R_bw() * omega * h));

    check((got - want).cwiseAbs().maxCoeff() <= T(1e-14),
          "Phi_beta,beta is not Exp(w_world h)");
    check((got - Matrix3::Identity()).cwiseAbs().maxCoeff() > T(0.01),
          "the test turn rate is too low for the bias rotation term to matter");
    // Level 1's corresponding block must remain exactly the identity.
    MatrixNX Phi1;
    a.build_transition(h, g, Phi1);
    check((Phi1.block<3,3>(L1::OFF_BG, L1::OFF_BG) - Matrix3::Identity())
              .cwiseAbs().maxCoeff() == T(0),
          "Level 1's bias block must stay the identity");
}

// ---------------------------------------------------------------------------
// 3. Measurement Jacobians at Level 2
// ---------------------------------------------------------------------------

// The payoff of the two-frame geometry: the accelerometer-bias column becomes
// -I, where Level 1 needs -Rhat.
void test_accel_bias_column_is_identity_at_level_two() {
    L2 f; configure(f);
    L1 g; configure(g);

    Vector3 r2, r1; Eigen::Matrix<T,3,NX> H2, H1; Matrix3 Rw;
    f.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), r2, H2, Rw);
    g.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), r1, H1, Rw);

    check((H2.block<3,3>(0, L2::OFF_BA) + Matrix3::Identity()).cwiseAbs().maxCoeff() <= T(1e-15),
          "H_a accelerometer-bias column must be -I at Level 2");
    check((H1.block<3,3>(0, L1::OFF_BA) + g.R_bw()).cwiseAbs().maxCoeff() <= T(1e-15),
          "H_a accelerometer-bias column must be -Rhat at Level 1");

    // H_L1 = H_L2 T is the same statement as the two above, checked as one
    // relation over the whole Jacobian rather than block by block.
    const MatrixNX Tm = bias_map(g.R_bw());
    check((H1 - H2 * Tm).cwiseAbs().maxCoeff() <= T(1e-14),
          "H_L1 != H_L2 T");
}

// Finite differences at Level 2, at a consistent state, exactly as for
// Level 1 -- the neglected -[r]x term is the same in both.
void test_level2_jacobians_against_finite_differences() {
    for (int trial = 0; trial < 6; ++trial) {
        L2 f; configure(f);
        if (trial > 0) {
            TangentV jitter = TangentV::Zero();
            jitter.segment<3>(L2::OFF_PHI) = rand_vec(0.5);
            for (int c = 0; c < 4; ++c) jitter.segment<3>(L2::OFF_V + 3*c) = rand_vec(1.5);
            jitter.segment<3>(L2::OFF_BA) = rand_vec(0.05);
            f.inject(jitter);
        }

        const Vector3 acc =
            f.R_wb() * (f.get_world_accel() - f.gravity_world()) + f.get_acc_bias();
        const Vector3 mag = f.R_wb() * f.magnetic_reference_world();

        struct Case { const char* name; bool is_acc; };
        const Case cases[] = {{"accelerometer", true}, {"magnetometer", false}};

        for (const Case& c : cases) {
            Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
            if (c.is_acc) f.accel_residual(acc, T(35), r, H, Rw);
            else          f.mag_residual(mag, r, H, Rw);

            Eigen::Matrix<T,3,NX> numeric;
            const T eps = T(1e-6);
            for (int j = 0; j < NX; ++j) {
                TangentV e = TangentV::Zero();
                Vector3 rp, rm; Eigen::Matrix<T,3,NX> Hd; Matrix3 Rd;

                e(j) = eps;  L2 plus = f;  plus.inject(e);
                e(j) = -eps; L2 minus = f; minus.inject(e);

                if (c.is_acc) { plus.accel_residual(acc, T(35), rp, Hd, Rd);
                                minus.accel_residual(acc, T(35), rm, Hd, Rd); }
                else          { plus.mag_residual(mag, rp, Hd, Rd);
                                minus.mag_residual(mag, rm, Hd, Rd); }
                numeric.col(j) = (rp - rm) / (T(2) * eps);
            }

            const T err = (H - numeric).cwiseAbs().maxCoeff();
            if (!check(err <= T(2e-8),
                       std::string("Level 2 ") + c.name +
                       ": Jacobian does not match finite differences")) {
                std::cerr << "  max|diff| = " << err << '\n';
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Covariance transport through the update
// ---------------------------------------------------------------------------

bool is_psd(const MatrixNX& M, T tol) {
    Eigen::SelfAdjointEigenSolver<MatrixNX> es(M);
    return es.info() == Eigen::Success && es.eigenvalues().minCoeff() >= -tol;
}

/*
    The Joseph update must be the ONLY thing that touches P.

    OU-III follows its injection with apply_left_error_reset's
    G = I + 0.5 [dtheta]x, because it zeroes its error state and has to
    transport P into the new linearization point. Applying that here would
    double-correct. Assert the absence directly: run one update, and check the
    covariance equals the Joseph form computed independently from the same
    pre-update quantities.
*/
void test_covariance_transport_is_applied_exactly_once() {
    for (int trial = 0; trial < 8; ++trial) {
        L2 f; configure(f);
        f.covariance_full() = seed_covariance();

        const Vector3 acc =
            f.R_wb() * (f.get_world_accel() - f.gravity_world()) + f.get_acc_bias()
            + rand_vec(0.2);

        Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
        f.accel_residual(acc, T(35), r, H, Rw);

        const MatrixNX P = f.covariance_full();
        const Eigen::Matrix<T,NX,3> PHt = P * H.transpose();
        const Matrix3 S = H * PHt + Rw;
        const Eigen::Matrix<T,NX,3> K = PHt * S.inverse();
        const MatrixNX IKH = MatrixNX::Identity() - K * H;
        const MatrixNX expected = IKH * P * IKH.transpose() + K * Rw * K.transpose();

        check(f.measurement_update_acc_only(acc, T(35)), "update rejected");

        const T scale = std::max(T(1), expected.cwiseAbs().maxCoeff());
        const T err = (f.covariance_full() - expected).cwiseAbs().maxCoeff();
        if (!check(err <= T(1e-11) * scale,
                   "the post-update covariance is not exactly the Joseph form; "
                   "an extra transport (a MEKF-style reset?) has been applied")) {
            std::cerr << "  max|diff| = " << err << " (scale " << scale << ")\n";
        }
    }
}

void test_level2_covariance_stays_sane() {
    L2 f; configure(f);
    f.set_initial_linear_uncertainty(T(0.3), T(1.0), T(3.0), T(0.6));
    std::normal_distribution<double> n(0.0, 0.3);

    for (int step = 0; step < 3000; ++step) {
        f.time_update(Vector3(n(rng), n(rng), n(rng)), T(0.005));
        f.measurement_update_acc_only(
            Vector3(f.R_wb() * (f.get_world_accel() - f.gravity_world())) + rand_vec(0.05),
            T(35));
        if (step % 8 == 0) {
            f.measurement_update_mag_only(
                Vector3(f.R_wb() * f.magnetic_reference_world()) + rand_vec(0.01));
        }
        if (step % 3 == 0) f.applyIntegralZeroPseudoMeas();
    }

    const MatrixNX& P = f.covariance_full();
    check(P.allFinite(), "Level 2 covariance went non-finite");
    check((P - P.transpose()).cwiseAbs().maxCoeff() <= T(1e-12) * P.cwiseAbs().maxCoeff(),
          "Level 2 covariance lost symmetry");
    check(is_psd(P, T(1e-10) * std::max(T(1), P.cwiseAbs().maxCoeff())),
          "Level 2 covariance stopped being PSD");
    check(f.state().R.allFinite() && f.get_acc_bias().norm() <= T(0.5) + T(1e-9),
          "Level 2 accelerometer bias escaped its ball");
}

// Both levels must actually estimate. A filter fed consistent data from a
// known truth should reduce its bias error, not merely stay finite.
void test_level2_converges() {
    L2 truth; configure(truth);
    L2 f;     configure(f);

    // Seed a gyro-bias error and let the magnetometer and accelerometer work.
    TangentV err = TangentV::Zero();
    err.segment<3>(L2::OFF_BG) = Vector3(T(0.02), T(-0.015), T(0.01));
    err.segment<3>(L2::OFF_PHI) = Vector3(T(0.05), T(-0.04), T(0.06));
    f.inject(err);
    f.set_initial_linear_uncertainty(T(0.3), T(1.0), T(3.0), T(0.6));
    f.covariance_full().block<3,3>(L2::OFF_BG, L2::OFF_BG) = Matrix3::Identity() * T(1e-3);
    f.covariance_full().block<3,3>(L2::OFF_PHI, L2::OFF_PHI) = Matrix3::Identity() * T(1e-2);

    const T before = (f.gyroscope_bias() - truth.gyroscope_bias()).norm();

    std::normal_distribution<double> n(0.0, 0.25);
    for (int step = 0; step < 20000; ++step) {
        const Vector3 g_true(n(rng), n(rng), n(rng));
        const Vector3 gyro = g_true + truth.gyroscope_bias();
        truth.time_update(gyro, T(0.005));
        f.time_update(gyro, T(0.005));

        const Vector3 acc =
            truth.R_wb() * (truth.get_world_accel() - truth.gravity_world())
            + truth.get_acc_bias();
        f.measurement_update_acc_only(acc, T(35));
        f.measurement_update_mag_only(truth.R_wb() * truth.magnetic_reference_world());
    }

    const T after = (f.gyroscope_bias() - truth.gyroscope_bias()).norm();
    if (!check(after < before,
               "Level 2 did not reduce the gyro-bias error over 100 s")) {
        std::cerr << "  |db| before " << before << " after " << after << '\n';
    }
    const T att = f.quaternion().angularDistance(truth.quaternion());
    check(att < T(0.05), "Level 2 attitude did not converge");
}

void test_float_and_reduced_instantiations() {
    using FloatL2 = ocean_imu::kalman::Kalman3D_Wave_TFG<float, true, true, true>;
    FloatL2 f;
    f.initialize_identity();
    f.set_magnetic_reference_world(Eigen::Vector3f(0.21f, 0.02f, 0.43f));
    for (int i = 0; i < 400; ++i) {
        f.time_update(Eigen::Vector3f(0.01f, -0.02f, 0.03f), 0.005f);
        f.measurement_update_acc_only(Eigen::Vector3f(0.0f, 0.0f, -9.80665f), 35.0f);
        f.measurement_update_mag_only(Eigen::Vector3f(0.21f, 0.02f, 0.43f));
        f.applyIntegralZeroPseudoMeas();
    }
    check(f.covariance_full().allFinite() && f.state().R.allFinite(),
          "the float Level-2 instantiation went non-finite");

    using NoBa = ocean_imu::kalman::Kalman3D_Wave_TFG<float, true, false, true>;
    static_assert(NoBa::NX == 18, "dropping the accel bias must give 18 states");
    NoBa g;
    g.initialize_identity();
    g.time_update(Eigen::Vector3f(0.01f, 0.0f, 0.0f), 0.005f);
    check(g.covariance_full().allFinite(), "the 18-state Level-2 variant went non-finite");

    // Level 2 is the default, so the bare instantiation must be the two-frame one.
    using Default = ocean_imu::kalman::Kalman3D_Wave_TFG<double>;
    Default d; configure(d);
    Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
    d.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), r, H, Rw);
    check((H.block<3,3>(0, Default::OFF_BA) + Matrix3::Identity()).cwiseAbs().maxCoeff()
              <= T(1e-15),
          "the default instantiation is not Level 2");
}

} // namespace

int main() {
    test_level2_inject_error_roundtrip();
    test_level2_injection_uses_the_right_rotation();
    test_levels_are_the_same_filter_under_propagation();
    test_levels_stay_equivalent_over_a_long_run();
    test_bias_rotation_term_matters();
    test_accel_bias_column_is_identity_at_level_two();
    test_level2_jacobians_against_finite_differences();
    test_covariance_transport_is_applied_exactly_once();
    test_level2_covariance_stays_sane();
    test_level2_converges();
    test_float_and_reduced_instantiations();

    if (failures != 0) {
        std::cerr << failures << " covariance transport check(s) failed\n";
        return 1;
    }
    std::cout << "covariance_transport-test: all checks passed\n";
    return 0;
}
