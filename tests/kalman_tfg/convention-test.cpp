/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    The convention contract, asserted rather than asserted-about.

    Three choices get conflated constantly, and two implementations that look
    opposite can be mathematically identical if only one or two are stated:

        1. rotation direction : R = R_bw, body to world       (this filter)
        2. error definition   : right-invariant, eta = Xhat o X^-1
        3. correction side    : left, Xhat+ = Exp_G(dxi) o Xhat-

    OU-III makes a different choice on (1): it stores a world-to-body
    quaternion C = R_wb and corrects it as C+ = Exp(dtheta) C-. That is also a
    left multiplication, but on the *other* rotation direction, so "both
    correct on the left" is not a statement of agreement. This file works out
    what the correspondence actually is and pins it numerically.

    The tests here deliberately use no filter code. They are about the frames
    and the sides, and they must stay valid before Kalman3D_Wave_TFG exists.
    The matching assertion that the public accessor still returns a
    world-to-body quaternion lands here in Phase B, when there is an accessor
    to assert on.
*/

#define EIGEN_NON_ARDUINO

#include "lie/SO3Jacobians.h"
#include "lie/TwoFrameGroup.h"

#include <cmath>
#include <iostream>
#include <random>
#include <string>

namespace lie = ocean_imu::lie;
namespace ou_detail = ocean_imu::kalman::ou_detail;

namespace {

using T = double;
using Matrix3 = Eigen::Matrix<T,3,3>;
using Vector3 = Eigen::Matrix<T,3,1>;
using Group = lie::TwoFrameGroup<T,4,2>;

constexpr T kTol = 1e-12;

int failures = 0;

bool check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
    return condition;
}

template<int R, int C>
bool mat_close(const Eigen::Matrix<T,R,C>& a, const Eigen::Matrix<T,R,C>& b, T tol = kTol) {
    return a.allFinite() && b.allFinite() && ((a - b).cwiseAbs().maxCoeff() <= tol);
}

std::mt19937 rng(20260806u);

Vector3 random_vector(double scale) {
    std::uniform_real_distribution<double> d(-scale, scale);
    return Vector3(d(rng), d(rng), d(rng));
}

Matrix3 random_rotation() {
    return lie::Exp<T>(random_vector(2.0));
}

Group random_group() {
    Group g;
    g.R = random_rotation();
    for (int i = 0; i < Group::kWorldCols; ++i) g.X.col(i) = random_vector(5.0);
    for (int j = 0; j < Group::kBiasCols; ++j) g.B.col(j) = random_vector(0.3);
    return g;
}

Group::Tangent random_tangent() {
    Group::Tangent xi;
    xi.setZero();
    xi.segment<3>(Group::OFF_PHI) = random_vector(0.8);
    for (int i = 0; i < Group::kWorldCols; ++i) xi.segment<3>(Group::world_offset(i)) = random_vector(2.0);
    for (int j = 0; j < Group::kBiasCols; ++j) xi.segment<3>(Group::bias_offset(j)) = random_vector(0.2);
    return xi;
}

// ---------------------------------------------------------------------------
// (1) Rotation direction: R maps body coordinates into world coordinates.
// ---------------------------------------------------------------------------

void test_rotation_direction() {
    // A rotation of +90 degrees about world z. Under the right-handed
    // convention this carries body +x onto world +y.
    const Vector3 phi(T(0), T(0), T(M_PI) / T(2));
    const Matrix3 R = lie::Exp<T>(phi);

    check(mat_close<3,1>(Vector3(R * Vector3::UnitX()), Vector3::UnitY(), 1e-12),
          "Exp is not the right-handed exponential: +90deg about z must send x to y");

    // R = R_bw, so R^T is the world-to-body map OU-III stores.
    const Vector3 v_body(T(1), T(2), T(3));
    const Vector3 v_world = R * v_body;
    check(mat_close<3,1>(Vector3(R.transpose() * v_world), v_body, 1e-12),
          "R^T must map world coordinates back to body coordinates");

    // The aliasing claim made in SO3Jacobians.h: OU-III's quaternion
    // correction builds exactly this exponential.
    for (int trial = 0; trial < 50; ++trial) {
        const Vector3 dtheta = random_vector(1.5);
        const Matrix3 from_quat = ou_detail::quat_from_delta_theta<T>(dtheta).toRotationMatrix();
        check(mat_close<3,3>(from_quat, lie::Exp<T>(dtheta), 1e-12),
              "quat_from_delta_theta must agree with Exp");
    }

    // And the warning in SO3Jacobians.h: rot_and_B_from_wt returns the
    // transpose, because OU-III propagates a world-to-body quaternion. If this
    // ever stops holding, every rotation in the filter silently inverts.
    for (int trial = 0; trial < 20; ++trial) {
        const Vector3 w = random_vector(3.0);
        const T dt = T(0.004);
        Matrix3 Rstep, Bstep;
        ou_detail::rot_and_B_from_wt<T>(w, dt, Rstep, Bstep);
        check(mat_close<3,3>(Rstep, Matrix3(lie::Exp<T>(Vector3(w * dt)).transpose()), 1e-12),
              "rot_and_B_from_wt must be the transpose of Exp(w dt)");
    }
}

// ---------------------------------------------------------------------------
// (2) Error definition: right-invariant means invariant under RIGHT
//     multiplication of both states. Left-invariant is the other one.
// ---------------------------------------------------------------------------

void test_invariance_direction() {
    for (int trial = 0; trial < 50; ++trial) {
        const Group X    = random_group();
        const Group Xhat = random_group();
        const Group G    = random_group();

        // eta_R = Xhat X^-1 is unchanged when both states are right-multiplied.
        const Group eta_R       = Xhat.compose(X.inverse());
        const Group eta_R_moved = Xhat.compose(G).compose(X.compose(G).inverse());
        check(mat_close<3,3>(eta_R.R, eta_R_moved.R) &&
              mat_close<3,Group::kWorldCols>(eta_R.X, eta_R_moved.X, 1e-11) &&
              mat_close<3,Group::kBiasCols>(eta_R.B, eta_R_moved.B, 1e-11),
              "the right-invariant error changed under right multiplication");

        // ... and is NOT invariant under left multiplication. Asserting the
        // negative keeps the test honest: without it, a degenerate error
        // definition would satisfy the positive case trivially.
        const Group eta_R_left = G.compose(Xhat).compose(G.compose(X).inverse());
        const bool unchanged = mat_close<3,3>(eta_R.R, eta_R_left.R, 1e-9) &&
                               mat_close<3,Group::kWorldCols>(eta_R.X, eta_R_left.X, 1e-9);
        check(!unchanged, "the right-invariant error must not survive left multiplication");

        // The mirror statement, for the convention this filter did not choose.
        const Group eta_L       = X.inverse().compose(Xhat);
        const Group eta_L_moved = G.compose(X).inverse().compose(G.compose(Xhat));
        check(mat_close<3,3>(eta_L.R, eta_L_moved.R) &&
              mat_close<3,Group::kWorldCols>(eta_L.X, eta_L_moved.X, 1e-11) &&
              mat_close<3,Group::kBiasCols>(eta_L.B, eta_L_moved.B, 1e-11),
              "the left-invariant error changed under left multiplication");
    }
}

// ---------------------------------------------------------------------------
// (3) Correction side: a right-invariant error is applied on the LEFT.
// ---------------------------------------------------------------------------

/*
    This is the consistency condition binding (2) to (3). If

        Xhat = Exp_G(xi) o X

    then the right-invariant error is exactly

        eta_R = Xhat o X^-1 = Exp_G(xi),   so   Log_G(Xhat o X^-1) = xi

    with no approximation. Applying the same xi on the right instead would
    make this identity false, which is what makes it a usable test rather than
    a restatement.
*/
void test_correction_side() {
    for (int trial = 0; trial < 50; ++trial) {
        const Group X = random_group();
        const Group::Tangent xi = random_tangent();

        const Group Xhat = X.retract(xi);
        const Group::Tangent recovered = Group::Log_G(Xhat.compose(X.inverse()));
        check(mat_close<Group::NX,1>(recovered, xi, 1e-10),
              "Log_G(Xhat X^-1) != xi: the correction side does not match the error definition");

        // retract() must be the left multiplication, not the right one.
        const Group left  = Group::Exp_G(xi).compose(X);
        const Group right = X.compose(Group::Exp_G(xi));
        check(mat_close<3,3>(Xhat.R, left.R) &&
              mat_close<3,Group::kWorldCols>(Xhat.X, left.X, 1e-11),
              "retract() is not left multiplication");
        const bool same_as_right = mat_close<3,3>(Xhat.R, right.R, 1e-9) &&
                                   mat_close<3,Group::kWorldCols>(Xhat.X, right.X, 1e-9);
        check(!same_as_right, "left and right retraction must differ for a generic element");
    }
}

// ---------------------------------------------------------------------------
// The OU-III bridge: the same physical perturbation, both conventions.
// ---------------------------------------------------------------------------

/*
    OU-III stores C = R_wb and injects on the left:

        C+ = Exp(dtheta) C-

    This filter stores R = C^T = R_bw. Transposing,

        R+ = (C+)^T = (C-)^T Exp(-dtheta) = R- Exp(-dtheta)

    so the very same physical rotation appears as a RIGHT multiplication on R.
    To express it as this filter's left-multiplicative correction, conjugate
    through R- Exp(u) R-^T = Exp(R- u):

        R+ = Exp(-R- dtheta) R-,     i.e.   phi = -R- dtheta

    The correspondence is therefore not phi = +-dtheta. It is a rotated,
    negated dtheta, and getting it wrong is a bug that survives every
    small-angle test because both forms agree to first order at R- = I.
*/
void test_ou_iii_bridge() {
    for (int trial = 0; trial < 100; ++trial) {
        const Matrix3 C_minus = random_rotation();          // OU-III: world -> body
        const Matrix3 R_minus = C_minus.transpose();        // this filter: body -> world
        const Vector3 dtheta = random_vector(1.2);

        // OU-III injection.
        const Matrix3 C_plus = ou_detail::quat_from_delta_theta<T>(dtheta).toRotationMatrix() * C_minus;

        // Same perturbation, seen as a right multiplication on R.
        check(mat_close<3,3>(Matrix3(C_plus.transpose()),
                             Matrix3(R_minus * lie::Exp<T>(Vector3(-dtheta))), 1e-12),
              "R+ != R- Exp(-dtheta): the transposed injection changed sides incorrectly");

        // Same perturbation, as this filter's left-multiplicative correction.
        const Vector3 phi = -(R_minus * dtheta);
        check(mat_close<3,3>(Matrix3(C_plus.transpose()),
                             Matrix3(lie::Exp<T>(phi) * R_minus), 1e-12),
              "phi = -R- dtheta does not reproduce the OU-III injection");

        // Drive it through the group retraction itself, so the bridge is
        // pinned against the code the filter will actually call.
        Group g;
        g.R = R_minus;
        Group::Tangent xi;
        xi.setZero();
        xi.segment<3>(Group::OFF_PHI) = phi;
        check(mat_close<3,3>(g.retract(xi).R, Matrix3(C_plus.transpose()), 1e-12),
              "TwoFrameGroup::retract does not reproduce the OU-III attitude injection");

        // The naive correspondences, which must NOT hold for a generic R-.
        // Both agree to first order at R- = I, so a test that only ever
        // exercised small rotations near identity would miss the error.
        if (dtheta.norm() > 0.2 && lie::Log<T>(R_minus).norm() > 0.2) {
            check(!mat_close<3,3>(Matrix3(C_plus.transpose()),
                                  Matrix3(lie::Exp<T>(dtheta) * R_minus), 1e-6),
                  "phi = +dtheta must not reproduce the injection for a generic R-");
            check(!mat_close<3,3>(Matrix3(C_plus.transpose()),
                                  Matrix3(lie::Exp<T>(Vector3(-dtheta)) * R_minus), 1e-6),
                  "phi = -dtheta must not reproduce the injection for a generic R-");
        }
    }
}

/*
    The physical content, stated without reference to either storage
    convention: a fixed world direction, and a fixed body direction, must land
    in the same place whichever filter applied the correction.
*/
void test_bridge_is_physical() {
    for (int trial = 0; trial < 50; ++trial) {
        const Matrix3 C_minus = random_rotation();
        const Matrix3 R_minus = C_minus.transpose();
        const Vector3 dtheta = random_vector(1.0);
        const Vector3 v_body = random_vector(1.0);
        const Vector3 v_world = random_vector(1.0);

        const Matrix3 C_plus = ou_detail::quat_from_delta_theta<T>(dtheta).toRotationMatrix() * C_minus;

        Group g;
        g.R = R_minus;
        Group::Tangent xi;
        xi.setZero();
        xi.segment<3>(Group::OFF_PHI) = -(R_minus * dtheta);
        const Matrix3 R_plus = g.retract(xi).R;

        check(mat_close<3,1>(Vector3(R_plus * v_body), Vector3(C_plus.transpose() * v_body), 1e-12),
              "the two conventions disagree on where a body vector lands in the world");
        check(mat_close<3,1>(Vector3(R_plus.transpose() * v_world), Vector3(C_plus * v_world), 1e-12),
              "the two conventions disagree on where a world vector lands in the body");
    }
}

} // namespace

int main() {
    test_rotation_direction();
    test_invariance_direction();
    test_correction_side();
    test_ou_iii_bridge();
    test_bridge_is_physical();

    if (failures != 0) {
        std::cerr << failures << " convention check(s) failed\n";
        return 1;
    }
    std::cout << "convention-test: all checks passed\n";
    return 0;
}
