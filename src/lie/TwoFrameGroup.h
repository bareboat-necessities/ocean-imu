#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    The right-invariant two-frame Lie group carrying the wave-state estimator.

        X = (R, X_w, B_b),   R   = R_bw in SO(3)          (body -> world)
                             X_w = [v p S a_w] in R^{3xNW}, world frame
                             B_b = [b_g b_a]   in R^{3xNB}, body frame

    Group law:

        (R1,X1,B1) o (R2,X2,B2) = ( R1 R2, X1 + R1 X2, B2 + R2^T B1 )
        identity                = (I, 0, 0)
        (R,X,B)^-1              = (R^T, -R^T X, -R B)

    The point of the law is the asymmetry in the last slot. World vectors are
    carried along by R; body-frame vectors are carried by R^T from the other
    side. A direct product SE_NW(3) x R^{3 NB} cannot express that, and it is
    precisely what makes an attitude correction rotate v, p, S, a_w coherently
    while transporting the biases in the frame they are measured in.

    CONVENTION CONTRACT. Three choices, always stated together, because two
    implementations that look opposite can be identical if only one or two are
    given:

        1. rotation direction : R = R_bw, body to world
        2. error definition   : right-invariant, eta = Xhat o X^-1
        3. correction side    : left, Xhat+ = Exp_G(dxi) o Xhat-

    "Right-invariant" describes the error -- it is invariant under right
    multiplication of both states -- and such an error is applied by LEFT
    multiplication. That is not a contradiction, and conflating the two is the
    usual way this gets implemented backwards.

    Note that OU-III stores the opposite rotation direction (world-to-body C =
    R_wb, corrected as C+ = Exp(dtheta) C-). Because R = C^T, the same physical
    perturbation appears there as a right multiplication:

        R+ = (C+)^T = (C-)^T Exp(-dtheta) = R- Exp(-dtheta)

    so "both correct on the left" is not a statement of agreement between the
    two filters. tests/kalman_tfg/convention-test.cpp pins the equivalence.

    Right-invariant error in local coordinates:

        dR     = Rhat R^T = Exp(phi)
        rho_x  = xhat - dR x        for each world column x
        beta_b = R (bhat - b)       for each body column b

    Tangent layout deliberately matches Kalman3D_Wave_OU_III.h:40-54, so the
    world block stays contiguous and the propagation/Joseph block structure
    ports over with the same offsets:

        xi = [ phi | beta_g | rho_v rho_p rho_S rho_a | beta_a ]
               0     3        6     9     12    15      18       (NW=4, NB=2)

    Scalar-templated throughout: firmware runs float, tests run double.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include "lie/SO3Jacobians.h"

namespace ocean_imu::lie {

/*
    NW world-frame columns and NB body-frame columns.

    The defaults are the OU-III state: [v p S a_w] and [b_g b_a]. NB = 1 drops
    the accelerometer bias, matching Kalman3D_Wave_OU_III's with_accel_bias =
    false instantiation.
*/
template<typename T = float, int NW = 4, int NB = 2>
struct TwoFrameGroup {
    static_assert(NW >= 1, "at least one world vector");
    static_assert(NB >= 1, "at least the gyroscope bias");

    static constexpr int kWorldCols = NW;
    static constexpr int kBiasCols  = NB;

    // Tangent dimension: 3 attitude + 3 per world column + 3 per bias column.
    static constexpr int NX = 3 + 3 * NW + 3 * NB;

    static constexpr int OFF_PHI = 0;
    static constexpr int OFF_BG  = 3;               // first bias column
    static constexpr int OFF_W   = 6;               // world block, contiguous
    static constexpr int OFF_BA  = OFF_W + 3 * NW;  // remaining bias columns

    // Tangent offset of world column i (0 = v, 1 = p, 2 = S, 3 = a_w).
    static constexpr int world_offset(int i) { return OFF_W + 3 * i; }

    // Tangent offset of bias column j. Column 0 is the gyro bias and keeps
    // OU-III's slot right after attitude; the rest follow the world block.
    static constexpr int bias_offset(int j) {
        return (j == 0) ? OFF_BG : (OFF_BA + 3 * (j - 1));
    }

    using Scalar     = T;
    using Matrix3    = Eigen::Matrix<T,3,3>;
    using Vector3    = Eigen::Matrix<T,3,1>;
    using WorldBlock = Eigen::Matrix<T,3,NW>;
    using BiasBlock  = Eigen::Matrix<T,3,NB>;
    using Tangent    = Eigen::Matrix<T,NX,1>;
    using AdjointMat = Eigen::Matrix<T,NX,NX>;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Matrix3    R{Matrix3::Identity()};
    WorldBlock X{WorldBlock::Zero()};
    BiasBlock  B{BiasBlock::Zero()};

    TwoFrameGroup() = default;
    TwoFrameGroup(const Matrix3& R_in, const WorldBlock& X_in, const BiasBlock& B_in)
        : R(R_in), X(X_in), B(B_in) {}

    static TwoFrameGroup Identity() { return TwoFrameGroup(); }

    // (R1,X1,B1) o (R2,X2,B2) = ( R1 R2, X1 + R1 X2, B2 + R2^T B1 )
    TwoFrameGroup compose(const TwoFrameGroup& other) const {
        return TwoFrameGroup(R * other.R,
                             X + R * other.X,
                             other.B + other.R.transpose() * B);
    }

    TwoFrameGroup operator*(const TwoFrameGroup& other) const { return compose(other); }

    // (R,X,B)^-1 = (R^T, -R^T X, -R B)
    TwoFrameGroup inverse() const {
        const Matrix3 Rt = R.transpose();
        return TwoFrameGroup(Rt, WorldBlock(-(Rt * X)), BiasBlock(-(R * B)));
    }

    /*
        Group exponential.

            E_R = Exp(phi)
            E_X = J_l(phi)  [rho_v rho_p rho_S rho_a]
            E_B = J_l(-phi) [beta_g beta_a]

        The J_l(-phi) on the bias block is derived, not chosen. The group law
        induces the one-parameter ODE b'(s) = beta - [phi]x b(s), whose
        solution at s = 1 is Exp(-phi) J_l(phi) beta, and Exp(-phi) J_l(phi)
        is the right Jacobian J_r(phi) = J_l(-phi).
    */
    static TwoFrameGroup Exp_G(const Tangent& xi) {
        const Vector3 phi = xi.template segment<3>(OFF_PHI);
        const Matrix3 Jl  = left_jacobian<T>(phi);
        const Matrix3 Jr  = left_jacobian<T>(Vector3(-phi));

        WorldBlock rho;
        for (int i = 0; i < NW; ++i) rho.col(i) = xi.template segment<3>(world_offset(i));
        BiasBlock beta;
        for (int j = 0; j < NB; ++j) beta.col(j) = xi.template segment<3>(bias_offset(j));

        return TwoFrameGroup(Exp<T>(phi), WorldBlock(Jl * rho), BiasBlock(Jr * beta));
    }

    // Group logarithm, the exact inverse of Exp_G.
    static Tangent Log_G(const TwoFrameGroup& g) {
        const Vector3 phi = Log<T>(g.R);
        const Matrix3 Jli = inv_left_jacobian<T>(phi);
        const Matrix3 Jri = inv_left_jacobian<T>(Vector3(-phi));

        const WorldBlock rho  = Jli * g.X;
        const BiasBlock  beta = Jri * g.B;

        Tangent xi;
        xi.setZero();
        xi.template segment<3>(OFF_PHI) = phi;
        for (int i = 0; i < NW; ++i) xi.template segment<3>(world_offset(i)) = rho.col(i);
        for (int j = 0; j < NB; ++j) xi.template segment<3>(bias_offset(j)) = beta.col(j);
        return xi;
    }

    /*
        Adjoint: X Exp(xi) X^-1 = Exp(Ad_X xi).

            phi    -> R phi
            rho_i  -> R rho_i + [x_i]x R phi
            beta_j -> R beta_j + R [b_j]x phi

        The off-diagonal blocks are what make an attitude uncertainty show up
        as translational uncertainty scaled by how far the vector reaches.
    */
    AdjointMat Adjoint() const {
        AdjointMat Ad;
        Ad.setZero();
        Ad.template block<3,3>(OFF_PHI, OFF_PHI) = R;

        for (int i = 0; i < NW; ++i) {
            const int off = world_offset(i);
            Ad.template block<3,3>(off, off)     = R;
            Ad.template block<3,3>(off, OFF_PHI) = skew<T>(Vector3(X.col(i))) * R;
        }
        for (int j = 0; j < NB; ++j) {
            const int off = bias_offset(j);
            Ad.template block<3,3>(off, off)     = R;
            Ad.template block<3,3>(off, OFF_PHI) = R * skew<T>(Vector3(B.col(j)));
        }
        return Ad;
    }

    /*
        Left retraction: the filter's state injection.

            Xhat+ = Exp_G(xi) o Xhat-

        Componentwise:

            Rhat+ = E_R Rhat-
            Xhat+ = E_R Xhat- + E_X
            Bhat+ = Bhat- + (Rhat-)^T E_B

        which is the whole point of the redesign: a correction to attitude
        rotates v, p, S and a_w with it instead of leaving them behind.
    */
    TwoFrameGroup retract(const Tangent& xi) const {
        return Exp_G(xi).compose(*this);
    }

    /*
        Inverse of retract(): the right-invariant error of *this relative to a
        reference, in local coordinates. local() satisfies

            reference.retract(x.local(reference)) == x
    */
    Tangent local(const TwoFrameGroup& reference) const {
        return Log_G(this->compose(reference.inverse()));
    }
};

} // namespace ocean_imu::lie
