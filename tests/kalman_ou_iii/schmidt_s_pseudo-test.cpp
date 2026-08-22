// Contract for the restricted (Schmidt-style) S=0 pseudo-measurement update.
//
// The S=0 pseudo-measurement is a translational regularizer: it asserts that
// the running integral of displacement has zero mean, and it observes S and
// nothing else.  The ordinary Kalman gain nevertheless corrects attitude
// through the cross covariance,
//
//     K_{theta S} = P_{theta S} (P_SS + R_S)^{-1},
//
// which makes an integral-channel residual drive an SO(3) injection.  The
// deployed update now freezes the attitude error state for this measurement,
// so E_theta K_S = 0 exactly, and uses the general-gain Joseph covariance
// update because the gain is no longer the minimum-covariance one.
//
// These tests pin that behaviour, and -- just as importantly -- pin that the
// path being removed was really there, that every other S-gain row still
// works, and that the accelerometer and magnetometer updates are untouched.
#define EIGEN_NON_ARDUINO

#include <cmath>
#include <iostream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private

namespace {

using Core = Kalman3D_Wave_OU_III<double, true, true>;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Mat21 = Eigen::Matrix<double, 21, 21>;
using Mat21x3 = Eigen::Matrix<double, 21, 3>;

constexpr int OFF_TH = 0;
constexpr int OFF_BG = 3;
constexpr int OFF_V  = 6;
constexpr int OFF_P  = 9;
constexpr int OFF_S  = 12;
constexpr int OFF_AW = 15;
constexpr int OFF_BA = 18;

int failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

void checkClose(double a, double b, double tol, const std::string& what) {
    if (!(std::fabs(a - b) <= tol)) {
        std::cerr << "FAIL: " << what << " (" << a << " vs " << b
                  << ", tol " << tol << ")\n";
        ++failures;
    }
}

Core makeFilter() {
    Core f(Vec3::Constant(0.0294), Vec3::Constant(0.000157), Vec3::Constant(0.36));
    f.set_mag_world_ref(Vec3(20.5, 0.0, 44.0));
    f.initialize_from_attitude(Eigen::Quaterniond::Identity(), 0.035, 0.087);
    f.set_linear_block_enabled(true);
    f.set_acc_bias_updates_enabled(true);
    return f;
}

// A covariance with deliberately large attitude-S coupling, plus coupling from
// every other block to S so the "other rows still act" test has something to
// measure.  Built as M M^T + eps I so it is symmetric positive definite by
// construction rather than by hope.
Mat21 coupledCovariance() {
    Mat21 M = Mat21::Zero();
    for (int i = 0; i < 21; ++i) M(i, i) = 0.35;
    // theta <-> S, deliberately strong.
    for (int i = 0; i < 3; ++i) {
        M(OFF_TH + i, OFF_S + i) = 0.9;
        M(OFF_S + i, OFF_TH + i) = 0.9;
    }
    // Every other block also couples to S.
    for (int i = 0; i < 3; ++i) {
        M(OFF_BG + i, OFF_S + i) = 0.20;
        M(OFF_V  + i, OFF_S + i) = 0.55;
        M(OFF_P  + i, OFF_S + i) = 0.60;
        M(OFF_AW + i, OFF_S + i) = 0.45;
        M(OFF_BA + i, OFF_S + i) = 0.30;
    }
    Mat21 P = M * M.transpose();
    P += Mat21::Identity() * 1e-3;
    return 0.5 * (P + P.transpose());
}

// The gain the filter would have used before the restriction.
Mat21x3 unconstrainedGain(const Mat21& P, const Mat3& R_S) {
    const Mat3 S = P.block<3,3>(OFF_S, OFF_S) + R_S;
    const Mat21x3 PCt = P.block<21,3>(0, OFF_S);
    return PCt * S.inverse();
}

double minEigenvalue(const Mat21& P) {
    Eigen::SelfAdjointEigenSolver<Mat21> es(P, Eigen::EigenvaluesOnly);
    return es.eigenvalues()(0);
}

// ---------------------------------------------------------------- 1 + 2 ----
void testAttitudeIsFrozenAndThePathWasReal() {
    Core f = makeFilter();
    f.Pext = coupledCovariance();
    // Nonzero innovation: the update targets S = 0, so any nonzero S gives one.
    f.xext.segment<3>(OFF_S) = Vec3(1.3, -0.7, 0.45);

    const Mat21 P_before = f.Pext;
    const Vec3 dtheta_before = f.xext.segment<3>(OFF_TH);
    const Eigen::Quaterniond q_before = f.quaternion_boat();

    // 2. The unconstrained gain really would have moved attitude, so the
    //    regression below is exercising the path we intend to remove.
    const Mat21x3 K_full = unconstrainedGain(P_before, f.R_S);
    const double Ktheta_norm = K_full.block<3,3>(OFF_TH, 0).norm();
    check(Ktheta_norm > 1e-3,
          "the unconstrained S gain has no attitude rows to remove; the test "
          "setup does not exercise the path");
    const Vec3 dtheta_would_have_been =
        K_full.block<3,3>(OFF_TH, 0) * (-f.xext.segment<3>(OFF_S));
    check(dtheta_would_have_been.norm() > 1e-3,
          "the unconstrained gain would not have corrected attitude here");

    f.applyIntegralZeroPseudoMeas();

    // 1. Attitude is frozen, in the error state and in the nominal quaternion.
    const Vec3 dtheta_after = f.xext.segment<3>(OFF_TH);
    checkClose((dtheta_after - dtheta_before).norm(), 0.0, 0.0,
               "the restricted S update moved the attitude error state");
    const Eigen::Quaterniond q_after = f.quaternion_boat();
    Eigen::Quaterniond dq = q_before.conjugate() * q_after;
    dq.normalize();
    checkClose(std::fabs(dq.w()), 1.0, 1e-15,
               "the restricted S update rotated the nominal quaternion");

    std::cout << "SCHMIDT_S removed_attitude_gain_norm=" << Ktheta_norm
              << " would_have_injected_rad=" << dtheta_would_have_been.norm()
              << "\n";
}

// -------------------------------------------------------------------- 3 ----
void testOtherRowsStillAct() {
    Core f = makeFilter();
    f.Pext = coupledCovariance();
    const Vec3 S0(1.3, -0.7, 0.45);
    f.xext.segment<3>(OFF_S) = S0;

    const Mat21 P_before = f.Pext;
    const Mat21x3 K_full = unconstrainedGain(P_before, f.R_S);
    const Vec3 r = -S0;

    Eigen::Matrix<double, 21, 1> x_before = f.xext;
    f.applyIntegralZeroPseudoMeas();

    struct Block { const char* name; int off; };
    const Block blocks[] = {
        {"S",   OFF_S},  {"p", OFF_P},  {"v", OFF_V},
        {"a_w", OFF_AW}, {"b_g", OFF_BG}, {"b_a", OFF_BA},
    };
    for (const Block& b : blocks) {
        const Vec3 expected = K_full.block<3,3>(b.off, 0) * r;
        const Vec3 actual =
            f.xext.segment<3>(b.off) - x_before.segment<3>(b.off);
        check(expected.norm() > 1e-6,
              std::string("test setup gives no ") + b.name + " correction");
        checkClose((actual - expected).norm(), 0.0, 1e-12,
                   std::string("the ") + b.name +
                       " row of the S gain no longer acts as the ordinary gain");
    }
}

// -------------------------------------------------------------------- 4 ----
void testJosephCovarianceMatchesTheRestrictedGain() {
    Core f = makeFilter();
    f.Pext = coupledCovariance();
    f.xext.segment<3>(OFF_S) = Vec3(0.8, 0.2, -1.1);

    const Mat21 P_before = f.Pext;
    Mat21x3 K = unconstrainedGain(P_before, f.R_S);
    K.block<3,3>(OFF_TH, 0).setZero();          // the restriction under test

    Eigen::Matrix<double, 3, 21> H = Eigen::Matrix<double, 3, 21>::Zero();
    H.block<3,3>(0, OFF_S) = Mat3::Identity();

    const Mat21 IKH = Mat21::Identity() - K * H;
    const Mat21 expected =
        IKH * P_before * IKH.transpose() + K * f.R_S * K.transpose();

    f.applyIntegralZeroPseudoMeas();

    // The reset is an exact no-op at dtheta = 0, so the filter covariance must
    // be the general-gain Joseph form and nothing else.
    const double err = (f.Pext - expected).cwiseAbs().maxCoeff();
    checkClose(err, 0.0, 1e-10,
               "the covariance is not (I-KH)P(I-KH)^T + K R K^T for the "
               "restricted gain");

    // And it is emphatically not the optimal-gain shortcut, which would be
    // wrong for a restricted gain.
    const Mat21 optimal_shortcut = IKH * P_before;
    const double shortcut_gap =
        (optimal_shortcut - expected).cwiseAbs().maxCoeff();
    check(shortcut_gap > 1e-9,
          "the two covariance forms coincide here, so this test cannot tell "
          "them apart");

    std::cout << "SCHMIDT_S joseph_max_abs_err=" << err
              << " optimal_shortcut_gap=" << shortcut_gap << "\n";
}

// ----------------------------------------------------------------- 5 + 6 ----
void testSymmetryAndPositiveSemidefiniteness() {
    Core f = makeFilter();
    f.Pext = coupledCovariance();

    double worst_asym = 0.0;
    double worst_min_eig = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 200; ++i) {
        f.xext.segment<3>(OFF_S) =
            Vec3(0.9 * std::cos(0.3 * i), 0.6 * std::sin(0.2 * i), 0.4);
        f.applyIntegralZeroPseudoMeas();
        worst_asym = std::max(worst_asym,
                              (f.Pext - f.Pext.transpose()).cwiseAbs().maxCoeff());
        worst_min_eig = std::min(worst_min_eig, minEigenvalue(f.Pext));
    }
    checkClose(worst_asym, 0.0, 1e-12,
               "repeated restricted S updates broke covariance symmetry");
    check(worst_min_eig > -1e-12,
          "repeated restricted S updates drove the covariance indefinite");

    std::cout << "SCHMIDT_S repeated_updates_worst_asymmetry=" << worst_asym
              << " worst_min_eigenvalue=" << worst_min_eig << "\n";
}

// -------------------------------------------------------------------- 7 ----
void testAccelAndMagStillCorrectAttitude() {
    // Both keep the ordinary full gain, so both must still move attitude when
    // given a residual that says the attitude is wrong.
    {
        Core f = makeFilter();
        const Eigen::Quaterniond q_before = f.quaternion_boat();
        // A specific force tilted away from the filter's idea of down.
        f.measurement_update_acc_only(Vec3(1.5, -0.9, -9.6), 35.0);
        Eigen::Quaterniond dq = q_before.conjugate() * f.quaternion_boat();
        dq.normalize();
        check(std::fabs(dq.w()) < 1.0 - 1e-12,
              "the accelerometer update no longer corrects attitude");
    }
    {
        Core f = makeFilter();
        const Eigen::Quaterniond q_before = f.quaternion_boat();
        // A field rotated about down, which is a pure heading residual.
        const double c = std::cos(0.25), s = std::sin(0.25);
        f.measurement_update_mag_only(Vec3(20.5 * c, -20.5 * s, 44.0));
        Eigen::Quaterniond dq = q_before.conjugate() * f.quaternion_boat();
        dq.normalize();
        check(std::fabs(dq.w()) < 1.0 - 1e-12,
              "the magnetometer update no longer corrects attitude");
    }
}

// -------------------------------------------------------------------- 8 ----
void testPredictionStillCouplesStates() {
    // The freeze applies to the S correction only.  Prediction must still let
    // attitude evolve, and must still rebuild the theta-S coupling the update
    // did not delete.
    Core f = makeFilter();
    f.Pext = coupledCovariance();
    const Eigen::Quaterniond q0 = f.quaternion_boat();

    f.xext.segment<3>(OFF_S) = Vec3(1.0, -0.5, 0.3);
    f.applyIntegralZeroPseudoMeas();

    const double coupling_after_update =
        f.Pext.block<3,3>(OFF_TH, OFF_S).norm();
    check(coupling_after_update > 0.0,
          "the restricted update zeroed P_theta_S; freezing the mean "
          "correction must not delete the covariance");

    for (int i = 0; i < 50; ++i) {
        f.time_update(Vec3(0.05, -0.02, 0.01), 0.005);
    }
    Eigen::Quaterniond dq = q0.conjugate() * f.quaternion_boat();
    dq.normalize();
    check(std::fabs(dq.w()) < 1.0 - 1e-9,
          "attitude no longer evolves through prediction");

    std::cout << "SCHMIDT_S P_theta_S_after_restricted_update="
              << coupling_after_update << "\n";
}

// Cross-covariance behaviour the PR has to report on.
void reportCrossCovariances() {
    auto run = [](bool restricted) {
        Core f = makeFilter();
        f.Pext = coupledCovariance();
        for (int i = 0; i < 200; ++i) {
            f.xext.segment<3>(OFF_S) =
                Vec3(0.9 * std::cos(0.3 * i), 0.6 * std::sin(0.2 * i), 0.4);
            if (restricted) {
                f.applyIntegralZeroPseudoMeas();
            } else {
                // Reproduce the pre-change update: ordinary gain, same Joseph
                // form, same reset path.
                const Mat3 S = f.Pext.block<3,3>(OFF_S, OFF_S) + f.R_S;
                const Mat21x3 PCt = f.Pext.block<21,3>(0, OFF_S);
                const Mat21x3 K = PCt * S.inverse();
                const Vec3 r = -f.xext.segment<3>(OFF_S);
                f.xext.noalias() += K * r;
                f.joseph_update3_(K, S, PCt);
                f.applyQuaternionCorrectionFromErrorState();
            }
        }
        return f.Pext;
    };
    const Mat21 P_old = run(false);
    const Mat21 P_new = run(true);
    auto blk = [](const Mat21& P, int a, int b) {
        return P.block<3,3>(a, b).norm();
    };
    std::cout << "SCHMIDT_S_XCOV block,unrestricted,restricted\n"
              << "SCHMIDT_S_XCOV P_theta_S," << blk(P_old, OFF_TH, OFF_S)
              << "," << blk(P_new, OFF_TH, OFF_S) << "\n"
              << "SCHMIDT_S_XCOV P_theta_aw," << blk(P_old, OFF_TH, OFF_AW)
              << "," << blk(P_new, OFF_TH, OFF_AW) << "\n"
              << "SCHMIDT_S_XCOV P_theta_p," << blk(P_old, OFF_TH, OFF_P)
              << "," << blk(P_new, OFF_TH, OFF_P) << "\n"
              << "SCHMIDT_S_XCOV P_aw_S," << blk(P_old, OFF_AW, OFF_S)
              << "," << blk(P_new, OFF_AW, OFF_S) << "\n"
              << "SCHMIDT_S_XCOV P_theta_theta," << blk(P_old, OFF_TH, OFF_TH)
              << "," << blk(P_new, OFF_TH, OFF_TH) << "\n";
}

}  // namespace

int main() {
    testAttitudeIsFrozenAndThePathWasReal();
    testOtherRowsStillAct();
    testJosephCovarianceMatchesTheRestrictedGain();
    testSymmetryAndPositiveSemidefiniteness();
    testAccelAndMagStillCorrectAttitude();
    testPredictionStillCouplesStates();
    reportCrossCovariances();

    if (failures) {
        std::cerr << failures << " restricted-S contract(s) failed\n";
        return 1;
    }
    std::cout << "SCHMIDT_S_PSEUDO all contracts hold\n";
    return 0;
}
