/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Measurement Jacobians for the invariant residuals, verified against
    central finite differences taken ON THE MANIFOLD.

    The perturbation goes through inject(), not by adding to a state vector.
    That matters: inject() applies Exp_G, so the derivative computed here is
    d r / d xi in exactly the coordinates P is expressed in and the Kalman
    gain operates in. An additive perturbation would differ by J_l and by the
    rotation of the world vectors, and would agree only to first order --
    passing at small steps while the real linearization is wrong.

    The sign convention is the one fixed in Kalman3D_Wave_TFG:

        r ~ H xi + noise,   xi the error OF THE ESTIMATE,   correction = -K r

    so H is literally d r / d xi and there is nothing to get backwards.

    Beyond matching finite differences, the tests assert the two properties
    that are the point of the redesign: the accelerometer and magnetometer
    Jacobians are built from FIXED references (gravity, the world magnetic
    vector) and must not change when the estimated attitude or the estimated
    wave acceleration changes.
*/

#define EIGEN_NON_ARDUINO

#include "kalman_tfg/Kalman3D_Wave_TFG.h"

#include <cmath>
#include <iostream>
#include <random>
#include <string>

namespace {

using T = double;
using Filter  = ocean_imu::kalman::Kalman3D_Wave_TFG<T, true, true, false>;
using Vector3 = Eigen::Matrix<T,3,1>;
using Matrix3 = Eigen::Matrix<T,3,3>;
using Tangent = Filter::Tangent;

constexpr int NX = Filter::NX;
using JacobianT = Eigen::Matrix<T,3,NX>;

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

Filter make_filter() {
    Filter f;
    f.set_aw_time_constant(T(1.7));
    f.set_aw_stationary_std(Vector3(T(0.6), T(0.6), T(0.32)));
    f.set_Racc_std(Vector3(T(0.4), T(0.35), T(0.5)));
    f.set_Rmag_std(Vector3(T(0.1), T(0.12), T(0.09)));
    f.set_RS_noise(Vector3(T(1.5), T(1.5), T(0.8)));
    f.set_magnetic_reference_world(Vector3(T(0.21), T(0.02), T(0.43)));
    // A genuinely rotated attitude and non-zero everything, so a term that is
    // wrong cannot pass by being multiplied into zero.
    f.initialize_from_truth(
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(
            T(0.83), Vector3(T(0.31), T(-0.62), T(0.72)).normalized())),
        Vector3(T(0.5), T(-0.8), T(1.1)),      // v
        Vector3(T(-1.7), T(2.4), T(0.9)),      // p
        Vector3(T(3.3), T(-2.1), T(1.4)),      // S
        Vector3(T(0.3), T(-0.5), T(0.7)),      // a_w
        Vector3(T(0.012), T(-0.008), T(0.004)),// b_g
        Vector3(T(0.02), T(0.05), T(-0.01)));  // b_a
    return f;
}

// ---------------------------------------------------------------------------
// Generic finite-difference driver: perturb the estimate through inject(),
// recompute the residual, central-difference.
// ---------------------------------------------------------------------------

template<typename ResidualFn>
JacobianT numeric_jacobian(const Filter& base, ResidualFn residual, T eps) {
    JacobianT numeric;
    for (int j = 0; j < NX; ++j) {
        Tangent e = Tangent::Zero();

        e(j) = eps;
        Filter plus = base;
        plus.inject(e);

        e(j) = -eps;
        Filter minus = base;
        minus.inject(e);

        numeric.col(j) = (residual(plus) - residual(minus)) / (T(2) * eps);
    }
    return numeric;
}

// Bitwise state comparison, for asserting that something did NOT happen.
//
// error_from() is the wrong instrument here: it computes Log(Rhat R_ref^T),
// and the floating-point product of a rotation with its own transpose is
// I +- 1e-17 rather than exactly I, so it reports ~1e-16 of "error" between
// two bitwise-identical states. For inertness the raw comparison is both
// exact and unambiguous.
bool states_identical(const Filter& a, const Filter& b) {
    return (a.state().R - b.state().R).cwiseAbs().maxCoeff() == T(0)
        && (a.state().X - b.state().X).cwiseAbs().maxCoeff() == T(0)
        && (a.state().B - b.state().B).cwiseAbs().maxCoeff() == T(0);
}

void compare(const JacobianT& analytic, const JacobianT& numeric, T tol,
             const std::string& tag) {
    const T err = (analytic - numeric).cwiseAbs().maxCoeff();
    if (!check(err <= tol, tag + ": Jacobian does not match finite differences")) {
        std::cerr << "  max|diff| = " << err << " (tol " << tol << ")\n";
        int bi = 0, bj = 0;
        T worst = T(0);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < NX; ++j) {
                const T d = std::abs(analytic(i,j) - numeric(i,j));
                if (d > worst) { worst = d; bi = i; bj = j; }
            }
        }
        std::cerr << "  worst entry (" << bi << "," << bj << "): analytic "
                  << analytic(bi,bj) << " numeric " << numeric(bi,bj) << '\n';
    }
}

// ---------------------------------------------------------------------------
// Accelerometer
// ---------------------------------------------------------------------------

// The specific force a perfectly calibrated sensor would report for this
// state, so the residual is evaluated near zero rather than at an arbitrary
// offset. No temperature term: the callers leave k_a at zero, and the
// temperature path gets its own dedicated check below.
Vector3 consistent_accel(const Filter& f) {
    return f.R_wb() * (f.get_world_accel() - f.gravity_world()) + f.get_acc_bias();
}

/*
    Finite differences at a CONSISTENT state, where r_a = 0.

    That is not a convenience: [g]x is the exact derivative precisely where
    the residual vanishes. Evaluating anywhere else would measure the
    deliberately neglected -[r_a]x term as well, and the honest way to handle
    that is a separate test that pins its exact value -- see below -- rather
    than a tolerance loose enough to swallow it.
*/
void test_accel_jacobian_at_consistent_state() {
    const T tempC = T(35);
    for (int trial = 0; trial < 10; ++trial) {
        Filter f = make_filter();
        if (trial > 0) {
            Tangent jitter = Tangent::Zero();
            jitter.segment<3>(Filter::OFF_PHI) = rand_vec(0.5);
            for (int c = 0; c < 4; ++c) jitter.segment<3>(Filter::OFF_V + 3*c) = rand_vec(1.5);
            jitter.segment<3>(Filter::OFF_BA) = rand_vec(0.05);
            f.inject(jitter);
        }
        const Vector3 acc = consistent_accel(f);

        Vector3 r; JacobianT H; Matrix3 Rw;
        f.accel_residual(acc, tempC, r, H, Rw);
        check(r.cwiseAbs().maxCoeff() <= T(1e-12),
              "the consistent-state accelerometer sample left a residual");

        const JacobianT numeric = numeric_jacobian(
            f,
            [&](const Filter& g) {
                Vector3 rr; JacobianT HH; Matrix3 RR;
                g.accel_residual(acc, tempC, rr, HH, RR);
                return rr;
            },
            T(1e-6));

        compare(H, numeric, T(2e-8), "accelerometer (consistent state)");
    }
}

/*
    Away from a consistent state the finite-difference Jacobian and the
    analytic one must differ by EXACTLY -[r_a]x.

    This is the test that makes the approximation honest. Asserting only that
    the gap is "small" would pass equally well if the gap were some other
    small thing, and would stop discriminating the moment a genuine error of
    similar magnitude appeared. Asserting the gap's exact structure means any
    additional error shows up immediately, however small.
*/
void test_accel_jacobian_neglected_term_is_exactly_the_innovation() {
    const T tempC = T(35);
    for (int trial = 0; trial < 6; ++trial) {
        Filter f = make_filter();
        const Vector3 acc = consistent_accel(f) + rand_vec(0.4);

        Vector3 r; JacobianT H; Matrix3 Rw;
        f.accel_residual(acc, tempC, r, H, Rw);
        check(r.norm() > T(0.1), "the deliberate residual is too small to test with");

        const JacobianT numeric = numeric_jacobian(
            f,
            [&](const Filter& g) {
                Vector3 rr; JacobianT HH; Matrix3 RR;
                g.accel_residual(acc, tempC, rr, HH, RR);
                return rr;
            },
            T(1e-6));

        const Matrix3 gap = numeric.block<3,3>(0, Filter::OFF_PHI)
                          - H.block<3,3>(0, Filter::OFF_PHI);
        const Matrix3 want = -ocean_imu::lie::skew<T>(r);
        const T err = (gap - want).cwiseAbs().maxCoeff();
        if (!check(err <= T(2e-8),
                   "the accelerometer Jacobian's neglected term is not exactly -[r_a]x")) {
            std::cerr << "  max|gap + [r]x| = " << err << ", |r| = " << r.norm() << '\n';
        }

        // Every other column must be exact regardless of the residual.
        JacobianT rest_a = H, rest_n = numeric;
        rest_a.block<3,3>(0, Filter::OFF_PHI).setZero();
        rest_n.block<3,3>(0, Filter::OFF_PHI).setZero();
        compare(rest_a, rest_n, T(2e-8), "accelerometer (non-attitude columns)");
    }
}

// The residual must vanish at a state consistent with the measurement --
// otherwise the Jacobian is the derivative of the wrong function.
void test_accel_residual_is_zero_at_truth() {
    Filter f = make_filter();
    const T tempC = T(41);
    f.set_accel_bias_temp_coeff(Vector3(T(1e-3), T(-2e-3), T(5e-4)));
    const Vector3 g = f.gravity_world();
    const Vector3 ba = f.get_acc_bias() + Vector3(T(1e-3), T(-2e-3), T(5e-4)) * (tempC - T(35));
    const Vector3 acc = f.R_wb() * (f.get_world_accel() - g) + ba;

    Vector3 r; JacobianT H; Matrix3 Rw;
    f.accel_residual(acc, tempC, r, H, Rw);
    check(r.cwiseAbs().maxCoeff() <= T(1e-13),
          "accelerometer residual is not zero at a consistent state "
          "(temperature term applied inconsistently?)");
}

/*
    The whole argument for the invariant residual is that H_a is built from
    the fixed gravity vector. Assert it directly: rotate the estimate, change
    the estimated wave acceleration, and the Jacobian's attitude block must
    not move.
*/
void test_accel_jacobian_is_estimate_independent() {
    Filter a = make_filter();
    Filter b = make_filter();

    Tangent big = Tangent::Zero();
    big.segment<3>(Filter::OFF_PHI) = Vector3(T(0.9), T(-1.2), T(0.7));
    big.segment<3>(Filter::OFF_AW)  = Vector3(T(2.0), T(-3.0), T(1.5));
    b.inject(big);

    Vector3 ra, rb; JacobianT Ha, Hb; Matrix3 Ra, Rb;
    a.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), ra, Ha, Ra);
    b.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), rb, Hb, Rb);

    check((Ha.block<3,3>(0, Filter::OFF_PHI) - Hb.block<3,3>(0, Filter::OFF_PHI))
              .cwiseAbs().maxCoeff() == T(0),
          "H_a attitude block changed with the estimate; it must depend only on g");
    check((Ha.block<3,3>(0, Filter::OFF_PHI)
           - ocean_imu::lie::skew<T>(a.gravity_world())).cwiseAbs().maxCoeff() == T(0),
          "H_a attitude block is not [g]x");
    check((Ha.block<3,3>(0, Filter::OFF_AW) + Matrix3::Identity()).cwiseAbs().maxCoeff() == T(0),
          "H_a wave-acceleration block is not -I");
}

/*
    The accelerometer-bias column is -Rhat at Level 1 and -I at Level 2.

    This is not a formatting detail. delta_ba is an additive BODY-frame error,
    so it has to be rotated into the world frame the residual lives in;
    beta_a = R(bhat_a - b_a) is already world-referred and needs no rotation.
    Writing -I here at Level 1 would be wrong by a full rotation, and would
    still look plausible because it is right at zero attitude.
*/
void test_accel_bias_column_is_rotated_at_level_one() {
    Filter f = make_filter();
    Vector3 r; JacobianT H; Matrix3 Rw;
    f.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), r, H, Rw);

    const Matrix3 got = H.block<3,3>(0, Filter::OFF_BA);
    check((got + f.R_bw()).cwiseAbs().maxCoeff() <= T(1e-15),
          "H_a accelerometer-bias column must be -Rhat at Level 1");
    // And it must be genuinely different from -I, or the test above proves
    // nothing about this state.
    check((got + Matrix3::Identity()).cwiseAbs().maxCoeff() > T(0.3),
          "the test attitude is too close to identity to distinguish -Rhat from -I");
}

// Freezing the accelerometer bias must remove its column entirely.
void test_accel_bias_freeze_zeroes_the_column() {
    Filter f = make_filter();
    f.set_acc_bias_updates_enabled(false);
    Vector3 r; JacobianT H; Matrix3 Rw;
    f.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), r, H, Rw);
    check(H.block<3,3>(0, Filter::OFF_BA).cwiseAbs().maxCoeff() == T(0),
          "a frozen accelerometer bias must not appear in H_a");
}

// ---------------------------------------------------------------------------
// Magnetometer
// ---------------------------------------------------------------------------

void test_mag_jacobian_at_consistent_state() {
    for (int trial = 0; trial < 10; ++trial) {
        Filter f = make_filter();
        if (trial > 0) {
            Tangent jitter = Tangent::Zero();
            jitter.segment<3>(Filter::OFF_PHI) = rand_vec(0.6);
            f.inject(jitter);
        }
        const Vector3 mag = f.R_wb() * f.magnetic_reference_world();

        Vector3 r; JacobianT H; Matrix3 Rw;
        f.mag_residual(mag, r, H, Rw);

        const JacobianT numeric = numeric_jacobian(
            f,
            [&](const Filter& g) {
                Vector3 rr; JacobianT HH; Matrix3 RR;
                g.mag_residual(mag, rr, HH, RR);
                return rr;
            },
            T(1e-6));

        compare(H, numeric, T(2e-8), "magnetometer (consistent state)");
    }
}

// Same structure as the accelerometer: the gap is exactly -[r_m]x.
void test_mag_jacobian_neglected_term_is_exactly_the_innovation() {
    for (int trial = 0; trial < 6; ++trial) {
        Filter f = make_filter();
        const Vector3 mag = f.R_wb() * f.magnetic_reference_world() + rand_vec(0.05);

        Vector3 r; JacobianT H; Matrix3 Rw;
        f.mag_residual(mag, r, H, Rw);
        check(r.norm() > T(0.01), "the deliberate magnetic residual is too small");

        const JacobianT numeric = numeric_jacobian(
            f,
            [&](const Filter& g) {
                Vector3 rr; JacobianT HH; Matrix3 RR;
                g.mag_residual(mag, rr, HH, RR);
                return rr;
            },
            T(1e-6));

        const Matrix3 gap = numeric.block<3,3>(0, Filter::OFF_PHI)
                          - H.block<3,3>(0, Filter::OFF_PHI);
        const T err = (gap + ocean_imu::lie::skew<T>(r)).cwiseAbs().maxCoeff();
        check(err <= T(2e-8),
              "the magnetometer Jacobian's neglected term is not exactly -[r_m]x");
    }
}

void test_mag_residual_is_zero_at_truth() {
    Filter f = make_filter();
    const Vector3 mag = f.R_wb() * f.magnetic_reference_world();
    Vector3 r; JacobianT H; Matrix3 Rw;
    f.mag_residual(mag, r, H, Rw);
    check(r.cwiseAbs().maxCoeff() <= T(1e-14),
          "magnetometer residual is not zero at the true attitude");
}

void test_mag_jacobian_uses_the_world_reference() {
    Filter a = make_filter();
    Filter b = make_filter();
    Tangent big = Tangent::Zero();
    big.segment<3>(Filter::OFF_PHI) = Vector3(T(1.1), T(-0.6), T(0.9));
    b.inject(big);

    Vector3 ra, rb; JacobianT Ha, Hb; Matrix3 Ra, Rb;
    const Vector3 mag(T(0.19), T(0.05), T(0.40));
    a.mag_residual(mag, ra, Ha, Ra);
    b.mag_residual(mag, rb, Hb, Rb);

    check((Ha - Hb).cwiseAbs().maxCoeff() == T(0),
          "H_m changed with the estimated attitude; it must depend only on B_w");
    check((Ha.block<3,3>(0, Filter::OFF_PHI)
           + ocean_imu::lie::skew<T>(a.magnetic_reference_world())).cwiseAbs().maxCoeff() == T(0),
          "H_m attitude block is not -[B_w]x");
    check(Ha.block<3,18>(0, 3).cwiseAbs().maxCoeff() == T(0),
          "the magnetometer must not touch anything but attitude");
}

// ---------------------------------------------------------------------------
// Integral pseudo-measurement
// ---------------------------------------------------------------------------

void test_integral_jacobian() {
    for (int trial = 0; trial < 10; ++trial) {
        Filter f = make_filter();
        if (trial > 0) {
            Tangent jitter = Tangent::Zero();
            jitter.segment<3>(Filter::OFF_PHI) = rand_vec(0.4);
            jitter.segment<3>(Filter::OFF_S)   = rand_vec(2.0);
            f.inject(jitter);
        }

        Vector3 r; JacobianT H; Matrix3 Rw;
        f.integral_residual(r, H, Rw);

        const JacobianT numeric = numeric_jacobian(
            f,
            [](const Filter& g) {
                Vector3 rr; JacobianT HH; Matrix3 RR;
                g.integral_residual(rr, HH, RR);
                return rr;
            },
            T(1e-6));

        // The [Shat]x term is what makes this match at a tight tolerance;
        // dropping it would leave an O(|S| |phi|) residual.
        compare(H, numeric, T(2e-8), "integral pseudo-measurement");
    }
}

void test_integral_keeps_the_S_cross_term() {
    Filter f = make_filter();
    Vector3 r; JacobianT H; Matrix3 Rw;
    f.integral_residual(r, H, Rw);

    check((H.block<3,3>(0, Filter::OFF_PHI)
           - ocean_imu::lie::skew<T>(f.get_integral_displacement())).cwiseAbs().maxCoeff() == T(0),
          "H_S attitude block is not [Shat]x");
    check((H.block<3,3>(0, Filter::OFF_S) + Matrix3::Identity()).cwiseAbs().maxCoeff() == T(0),
          "H_S integral block is not -I");
    check((r + f.get_integral_displacement()).cwiseAbs().maxCoeff() == T(0),
          "integral residual is not -Shat");
}

// ---------------------------------------------------------------------------
// Noise rotation
// ---------------------------------------------------------------------------

void test_measurement_noise_is_rotated_into_the_world() {
    Filter f = make_filter();
    Vector3 r; JacobianT H; Matrix3 Rw;

    f.set_Racc_std(Vector3(T(0.4), T(0.2), T(0.7)));   // anisotropic
    f.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), r, H, Rw);
    Matrix3 body; body.setZero();
    body.diagonal() << T(0.16), T(0.04), T(0.49);
    check((Rw - f.R_bw() * body * f.R_wb()).cwiseAbs().maxCoeff() <= T(1e-14),
          "accelerometer noise is not Rhat R_body Rhat^T");
    check((Rw - body).cwiseAbs().maxCoeff() > T(1e-3),
          "the test attitude does not actually exercise the rotation");

    // Isotropic noise must pass through unchanged -- a useful invariant, and
    // the reason this transformation is invisible in most configurations.
    f.set_Racc_std(Vector3::Constant(T(0.4)));
    f.accel_residual(Vector3(T(0.1), T(0.2), T(-9.6)), T(35), r, H, Rw);
    check((Rw - Matrix3::Identity() * T(0.16)).cwiseAbs().maxCoeff() <= T(1e-15),
          "isotropic accelerometer noise must be unchanged by the world rotation");

    // R_S is already world frame and must NOT be rotated.
    f.set_RS_noise(Vector3(T(1.5), T(1.5), T(0.8)));
    f.integral_residual(r, H, Rw);
    Matrix3 want; want.setZero();
    // set_RS_noise takes standard deviations, so R_S is their square. Write
    // the squares as products rather than as literals: 0.8*0.8 is
    // 0.6400000000000001 in double, so a literal 0.64 would fail the exact
    // comparison below for a reason that has nothing to do with frames.
    want.diagonal() << T(1.5)*T(1.5), T(1.5)*T(1.5), T(0.8)*T(0.8);
    check((Rw - want).cwiseAbs().maxCoeff() == T(0),
          "R_S must stay in world/NED coordinates, not be rotated by the estimate");
}

// ---------------------------------------------------------------------------
// The updates themselves
// ---------------------------------------------------------------------------

// An update must move the state toward the measurement and shrink the
// covariance in the observed direction.
void test_updates_reduce_error_and_covariance() {
    // Accelerometer: seed a wave-acceleration error and check it shrinks.
    {
        Filter truth = make_filter();
        Filter f = make_filter();
        Tangent err = Tangent::Zero();
        err.segment<3>(Filter::OFF_AW) = Vector3(T(0.4), T(-0.3), T(0.25));
        f.inject(err);
        f.set_initial_linear_uncertainty(T(0.3), T(1.0), T(3.0), T(0.6));

        const Vector3 acc = consistent_accel(truth);
        const T before = f.error_from(truth).segment<3>(Filter::OFF_AW).norm();
        const T trace_before = f.covariance_full().trace();
        check(f.measurement_update_acc_only(acc, T(35)), "accelerometer update was rejected");
        const T after = f.error_from(truth).segment<3>(Filter::OFF_AW).norm();

        check(after < before, "the accelerometer update did not reduce the a_w error");
        check(f.covariance_full().trace() < trace_before,
              "the accelerometer update did not reduce covariance");
        check(f.lastAccDiag().accepted && f.lastAccDiag().nis >= T(0),
              "accelerometer diagnostics not populated");
    }

    // Magnetometer: seed a yaw error.
    {
        Filter truth = make_filter();
        Filter f = make_filter();
        Tangent err = Tangent::Zero();
        err.segment<3>(Filter::OFF_PHI) = Vector3(T(0), T(0), T(0.15));
        f.inject(err);

        const Vector3 mag = truth.R_wb() * truth.magnetic_reference_world();
        const T before = f.error_from(truth).segment<3>(Filter::OFF_PHI).norm();
        check(f.measurement_update_mag_only(mag), "magnetometer update was rejected");
        const T after = f.error_from(truth).segment<3>(Filter::OFF_PHI).norm();
        check(after < before, "the magnetometer update did not reduce the attitude error");
    }

    // Integral: S must be pulled toward zero.
    {
        Filter f = make_filter();
        f.set_initial_linear_uncertainty(T(0.3), T(1.0), T(3.0), T(0.6));
        const T before = f.get_integral_displacement().norm();
        check(f.applyIntegralZeroPseudoMeas(), "integral pseudo-measurement was rejected");
        check(f.get_integral_displacement().norm() < before,
              "the integral pseudo-measurement did not pull S toward zero");
    }
}

// Covariance must stay symmetric and PSD through a long measurement sequence.
void test_updates_keep_covariance_sane() {
    Filter f = make_filter();
    f.set_initial_linear_uncertainty(T(0.3), T(1.0), T(3.0), T(0.6));
    std::normal_distribution<double> n(0.0, 0.2);

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

    const auto& P = f.covariance_full();
    check(P.allFinite(), "covariance went non-finite over a filtered run");
    check((P - P.transpose()).cwiseAbs().maxCoeff() <= T(1e-12) * P.cwiseAbs().maxCoeff(),
          "covariance lost symmetry under the Joseph update");
    Eigen::SelfAdjointEigenSolver<Filter::MatrixNX> es(P);
    check(es.info() == Eigen::Success &&
          es.eigenvalues().minCoeff() >= -T(1e-10) * std::max(T(1), P.cwiseAbs().maxCoeff()),
          "covariance stopped being PSD under the Joseph update");
    check(f.state().R.allFinite() && f.get_acc_bias().norm() <= T(0.5) + T(1e-9),
          "accelerometer bias escaped its ball");
}

// Rejected updates must leave the state and covariance completely untouched.
void test_rejected_updates_are_inert() {
    Filter f = make_filter();
    f.set_meas_gate_nis(T(1e-6));   // reject essentially everything
    const Filter before = f;
    const auto P_before = f.covariance_full();

    const bool accepted = f.measurement_update_acc_only(Vector3(T(5), T(5), T(5)), T(35));
    check(!accepted, "the NIS gate did not reject an absurd measurement");
    check(!f.lastAccDiag().accepted, "diagnostics claim acceptance after rejection");
    check(states_identical(f, before), "a rejected update modified the state");
    check((f.covariance_full() - P_before).cwiseAbs().maxCoeff() == T(0),
          "a rejected update modified the covariance");

    // A missing magnetic reference must also be inert, not a silent update
    // against the default unit-x vector.
    Filter g;
    g.initialize_identity();
    check(!g.has_magnetic_reference(), "magnetic reference should start unset");
    const auto Pg = g.covariance_full();
    check(!g.measurement_update_mag_only(Vector3(T(0.2), T(0), T(0.4))),
          "magnetometer update ran without a world reference");
    check((g.covariance_full() - Pg).cwiseAbs().maxCoeff() == T(0),
          "an update without a magnetic reference modified the covariance");
}

// Non-finite input must not poison the filter.
void test_non_finite_input_is_ignored() {
    Filter f = make_filter();
    const auto P_before = f.covariance_full();
    const Filter before = f;
    const T nan = std::numeric_limits<T>::quiet_NaN();
    check(!f.measurement_update_acc_only(Vector3(nan, T(0), T(0)), T(35)),
          "a NaN accelerometer sample was accepted");
    check(!f.measurement_update_mag_only(Vector3(nan, T(0), T(0))),
          "a NaN magnetometer sample was accepted");
    check(states_identical(f, before) &&
          (f.covariance_full() - P_before).cwiseAbs().maxCoeff() == T(0),
          "a NaN sample perturbed the filter");
}

void test_float_instantiation() {
    using FloatFilter = ocean_imu::kalman::Kalman3D_Wave_TFG<float, true, true, false>;
    FloatFilter f;
    f.initialize_identity();
    f.set_magnetic_reference_world(Eigen::Vector3f(0.21f, 0.02f, 0.43f));
    for (int i = 0; i < 400; ++i) {
        f.time_update(Eigen::Vector3f(0.01f, -0.02f, 0.03f), 0.005f);
        f.measurement_update_acc_only(Eigen::Vector3f(0.0f, 0.0f, -9.80665f), 35.0f);
        f.measurement_update_mag_only(Eigen::Vector3f(0.21f, 0.02f, 0.43f));
        f.applyIntegralZeroPseudoMeas();
    }
    check(f.covariance_full().allFinite() && f.state().R.allFinite(),
          "the float instantiation went non-finite under measurement updates");
}

} // namespace

int main() {
    test_accel_residual_is_zero_at_truth();
    test_accel_jacobian_at_consistent_state();
    test_accel_jacobian_neglected_term_is_exactly_the_innovation();
    test_accel_jacobian_is_estimate_independent();
    test_accel_bias_column_is_rotated_at_level_one();
    test_accel_bias_freeze_zeroes_the_column();

    test_mag_residual_is_zero_at_truth();
    test_mag_jacobian_at_consistent_state();
    test_mag_jacobian_neglected_term_is_exactly_the_innovation();
    test_mag_jacobian_uses_the_world_reference();

    test_integral_jacobian();
    test_integral_keeps_the_S_cross_term();

    test_measurement_noise_is_rotated_into_the_world();
    test_updates_reduce_error_and_covariance();
    test_updates_keep_covariance_sane();
    test_rejected_updates_are_inert();
    test_non_finite_input_is_ignored();
    test_float_instantiation();

    if (failures != 0) {
        std::cerr << failures << " Jacobian check(s) failed\n";
        return 1;
    }
    std::cout << "tfg_jacobians-test: all checks passed\n";
    return 0;
}
