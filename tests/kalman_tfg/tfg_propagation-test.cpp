/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    The discrete transition Phi must be the Jacobian of the propagation the
    filter actually performs -- not of the propagation someone meant to write.

    So this test never compares Phi against a second closed form. It perturbs
    the state through inject(), propagates both the perturbed and the nominal
    state with the same measured gyro, reads the resulting right-invariant
    error back out through error_from(), and central-differences that:

        Phi[:,j] ~ ( err( f(inject(X, +h e_j)), f(X) )
                   - err( f(inject(X, -h e_j)), f(X) ) ) / 2h

    inject() and error_from() are an exact inverse pair (asserted below), so
    this closes the loop through the same coordinates the covariance lives in.
    A Phi that is right for the wrong parameterization fails here.

    Everything runs at double. Float is exercised by the compile-only
    instantiation at the bottom -- finite differences at float precision would
    measure rounding, not the Jacobian.
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
using MatrixNX = Filter::MatrixNX;

constexpr int NX = Filter::NX;

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

Filter make_filter(T tau = T(1.7)) {
    Filter f;
    f.set_aw_time_constant(tau);
    f.set_aw_stationary_std(Vector3(T(0.6), T(0.6), T(0.32)));
    f.set_gyro_noise_density_rad_sqrt_s(T(0.004));
    f.set_Q_bgyro_rw(Vector3::Constant(T(1e-10)));
    // set_Q_bacc_rw follows current OU-III semantics: continuous driving
    // standard deviation in m/s^2/sqrt(s), not variance per second.
    f.set_Q_bacc_rw(Vector3::Constant(T(5e-4)));
    f.initialize_from_truth(
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(
            T(0.63), Vector3(T(0.2), T(-0.7), T(0.68)).normalized())),
        Vector3(T(0.5), T(-0.8), T(1.1)),
        Vector3(T(-1.7), T(2.4), T(0.9)),
        Vector3(T(3.3), T(-2.1), T(1.4)),
        Vector3(T(0.3), T(-0.5), T(0.7)),
        Vector3(T(0.012), T(-0.008), T(0.004)),
        Vector3(T(0.02), T(0.05), T(-0.01)));
    return f;
}

// ---------------------------------------------------------------------------
// inject() and error_from() must be exact inverses. Everything else in this
// file depends on that, so it is checked first and on its own.
// ---------------------------------------------------------------------------

void test_inject_error_roundtrip() {
    for (int trial = 0; trial < 50; ++trial) {
        const Filter base = make_filter();

        Tangent xi;
        xi.setZero();
        xi.segment<3>(Filter::OFF_PHI) = rand_vec(0.8);
        xi.segment<3>(Filter::OFF_BG)  = rand_vec(0.05);
        for (int c = 0; c < 4; ++c) xi.segment<3>(Filter::OFF_V + 3*c) = rand_vec(2.0);
        xi.segment<3>(Filter::OFF_BA)  = rand_vec(0.1);

        Filter perturbed = base;
        perturbed.inject(xi);

        const Tangent back = perturbed.error_from(base);
        const T err = (back - xi).cwiseAbs().maxCoeff();
        if (!check(err <= T(1e-12),
                   "error_from(inject(xi)) != xi")) {
            std::cerr << "  max|diff| = " << err << '\n';
        }
    }

    // A zero correction must be exactly a no-op, including on the attitude.
    Filter f = make_filter();
    const Filter before = f;
    f.inject(Tangent::Zero());
    check(f.error_from(before).cwiseAbs().maxCoeff() <= T(1e-15),
          "injecting a zero correction changed the state");
}

// ---------------------------------------------------------------------------
// Phi against central finite differences of the real propagation.
// ---------------------------------------------------------------------------

MatrixNX numeric_transition(const Filter& base, const Vector3& gyro, T Ts, T eps) {
    MatrixNX numeric;
    Filter nominal = base;
    nominal.propagate_mean(gyro, Ts);

    for (int j = 0; j < NX; ++j) {
        Tangent e = Tangent::Zero();

        e(j) = eps;
        Filter plus = base;
        plus.inject(e);
        plus.propagate_mean(gyro, Ts);

        e(j) = -eps;
        Filter minus = base;
        minus.inject(e);
        minus.propagate_mean(gyro, Ts);

        numeric.col(j) = (plus.error_from(nominal) - minus.error_from(nominal)) / (T(2) * eps);
    }
    return numeric;
}

void compare_transition(const Filter& f, const Vector3& gyro, T Ts, T tol,
                        const std::string& tag) {
    MatrixNX analytic;
    f.build_transition(Ts, gyro, analytic);
    const MatrixNX numeric = numeric_transition(f, gyro, Ts, T(1e-6));

    const T err = (analytic - numeric).cwiseAbs().maxCoeff();
    if (!check(err <= tol, tag + ": Phi does not match finite differences")) {
        std::cerr << "  max|diff| = " << err << " (tol " << tol << ")\n";
        // Name the worst block, so a failure points at the derivation that
        // broke rather than at a 21x21 wall of numbers.
        int bi = 0, bj = 0;
        T worst = T(0);
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NX; ++j) {
                const T d = std::abs(analytic(i,j) - numeric(i,j));
                if (d > worst) { worst = d; bi = i; bj = j; }
            }
        }
        std::cerr << "  worst entry (" << bi << "," << bj << "): analytic "
                  << analytic(bi,bj) << " numeric " << numeric(bi,bj) << '\n';
    }
}

void test_transition_matches_finite_differences() {
    // Still, moderate, and fast rotation. The turn rate is the whole point of
    // the Rbar = Rhat J_l(w h) factor, so it has to be exercised well away
    // from zero.
    const Vector3 gyros[] = {
        Vector3::Zero(),
        Vector3(T(0.05), T(-0.03), T(0.02)),
        Vector3(T(0.4), T(-0.7), T(0.5)),
        Vector3(T(2.0), T(-1.5), T(3.0)),
    };
    const T steps[] = {T(0.001), T(0.005), T(0.02), T(0.1)};

    for (const Vector3& g : gyros) {
        for (T h : steps) {
            Filter f = make_filter();
            const std::string tag =
                "|w|=" + std::to_string(g.norm()) + " h=" + std::to_string(h);
            // The one documented approximation is in Phi_world,bg, which
            // neglects the correlation of Rhat(s) with the Phi_OU weighting
            // and the drift of x(s) across the step. Both are O(h) on a term
            // already of size |x| |w| h, so the tolerance scales with h^2.
            const T tol = T(2e-8) + T(4.0) * h * h * (T(1) + g.norm());
            compare_transition(f, g, h, tol, tag);
        }
    }
}

// At zero turn rate J_l(0) = I, so Rbar reduces to Rhat exactly and the
// attitude/bias block has no approximation left in it.
//
// "At rest" means the *body* is not turning, which is w = gyro - b_g = 0, so
// the gyro must read the bias. Feeding a zero gyro instead leaves w = -b_g,
// a real turn rate, and J_l(-b_g h) is then not the identity -- an easy way
// to write a test that fails for a reason that has nothing to do with the
// code under test.
void test_attitude_bias_block_is_exact_at_rest() {
    for (T h : {T(0.001), T(0.005), T(0.05), T(0.5)}) {
        Filter f = make_filter();
        MatrixNX Phi;
        f.build_transition(h, f.gyroscope_bias(), Phi);

        const Matrix3 want = -f.R_bw() * h;
        const Matrix3 got = Phi.block<3,3>(Filter::OFF_PHI, Filter::OFF_BG);
        check((got - want).cwiseAbs().maxCoeff() <= T(1e-15),
              "Phi_phi,bg != -Rhat h at zero turn rate");

        check((Phi.block<3,3>(Filter::OFF_PHI, Filter::OFF_PHI)
               - Matrix3::Identity()).cwiseAbs().maxCoeff() == T(0),
              "Phi_phi,phi must be exactly I: the right-invariant attitude "
              "error is driftless");
        check(Phi.block<3,12>(Filter::OFF_PHI, Filter::OFF_V).cwiseAbs().maxCoeff() == T(0),
              "attitude error must not depend on the world vectors");
        check(Phi.block<12,3>(Filter::OFF_V, Filter::OFF_PHI).cwiseAbs().maxCoeff() == T(0),
              "world vectors must not couple back from phi: phidot has no phi term");
    }
}

// The Rbar factor is not cosmetic. Rebuilding Phi with Rhat in place of
// Rhat J_l(w h) must measurably disagree with finite differences at a
// realistic turn rate -- otherwise the correction is untested.
void test_average_rotation_factor_matters() {
    const Vector3 gyro(T(1.0), T(-0.5), T(0.8));
    const T h = T(0.02);
    Filter f = make_filter();

    MatrixNX analytic;
    f.build_transition(h, gyro, analytic);
    const MatrixNX numeric = numeric_transition(f, gyro, h, T(1e-6));

    const Matrix3 with_Jl = analytic.block<3,3>(Filter::OFF_PHI, Filter::OFF_BG);
    const Matrix3 naive   = -f.R_bw() * h;
    const Matrix3 truth   = numeric.block<3,3>(Filter::OFF_PHI, Filter::OFF_BG);

    const T err_correct = (with_Jl - truth).cwiseAbs().maxCoeff();
    const T err_naive   = (naive   - truth).cwiseAbs().maxCoeff();

    check(err_correct < err_naive * T(1e-2),
          "the Rbar = Rhat J_l(w h) factor is not doing measurable work; "
          "either the test is too weak or the factor is wrong");
    std::cout << "  Rbar factor: |Phi_phi,bg - FD| = " << err_correct
              << " with J_l, " << err_naive << " without\n";
}

// ---------------------------------------------------------------------------
// Covariance propagation
// ---------------------------------------------------------------------------

bool is_psd(const MatrixNX& M, T tol) {
    Eigen::SelfAdjointEigenSolver<MatrixNX> es(M);
    if (es.info() != Eigen::Success) return false;
    return es.eigenvalues().minCoeff() >= -tol;
}

void test_covariance_stays_psd_and_symmetric() {
    Filter f = make_filter();
    f.set_initial_linear_uncertainty(T(0.3), T(1.0), T(3.0), T(0.5));

    std::normal_distribution<double> n(0.0, 0.35);
    const T h = T(0.005);
    for (int step = 0; step < 4000; ++step) {
        const Vector3 gyro(n(rng), n(rng), n(rng));
        f.time_update(gyro, h);
    }

    const MatrixNX& P = f.covariance_full();
    check(P.allFinite(), "covariance went non-finite over a 20 s propagation");
    check((P - P.transpose()).cwiseAbs().maxCoeff() <= T(1e-14) * P.cwiseAbs().maxCoeff(),
          "covariance lost symmetry");
    check(is_psd(P, T(1e-12) * std::max(T(1), P.cwiseAbs().maxCoeff())),
          "covariance stopped being positive semi-definite");
    check(is_psd(f.covariance_full(), T(1e-12) * std::max(T(1), P.cwiseAbs().maxCoeff())),
          "covariance not PSD");

    // Attitude uncertainty must actually grow without measurements; a Phi or Q
    // that silently zeroed the gyro noise would still pass the PSD checks.
    check(P(0,0) > T(1e-8), "attitude variance did not grow over 20 s of gyro-only propagation");
}

void test_process_noise_is_psd() {
    for (T h : {T(0.001), T(0.005), T(0.05)}) {
        Filter f = make_filter();
        MatrixNX Qd;
        f.build_process_noise(h, Vector3(T(0.1), T(-0.2), T(0.05)), Qd);
        check(Qd.allFinite(), "Qd non-finite");
        check((Qd - Qd.transpose()).cwiseAbs().maxCoeff() <= T(1e-15) * std::max(T(1), Qd.cwiseAbs().maxCoeff()),
              "Qd is not symmetric");
        check(is_psd(Qd, T(1e-14) * std::max(T(1), Qd.cwiseAbs().maxCoeff())), "Qd is not PSD");
    }
}

// ---------------------------------------------------------------------------
// Mean propagation
// ---------------------------------------------------------------------------

// With a perfectly known bias and no rotation, the world chain must integrate
// the OU kinematics: p advances by v h + O(h^2), S by p h, and a_w decays.
void test_mean_propagation_follows_the_chain() {
    Filter f = make_filter(T(2.0));
    f.initialize_from_truth(Eigen::Quaternion<T>::Identity(),
                            Vector3(T(1), T(0), T(0)),     // v
                            Vector3(T(0), T(2), T(0)),     // p
                            Vector3(T(0), T(0), T(3)),     // S
                            Vector3(T(0.5), T(0), T(0)),   // a_w
                            Vector3::Zero(), Vector3::Zero());

    const T h = T(0.01);
    const Vector3 v0 = f.get_velocity();
    const Vector3 p0 = f.get_position();
    const Vector3 S0 = f.get_integral_displacement();
    const Vector3 a0 = f.get_world_accel();

    f.propagate_mean(Vector3::Zero(), h);

    check((f.get_world_accel() - a0 * std::exp(-h / T(2.0))).cwiseAbs().maxCoeff() <= T(1e-12),
          "a_w did not decay at the OU rate");
    check(std::abs(f.get_position().y() - (p0.y() + v0.y()*h)) <= T(1e-12),
          "p did not integrate v");
    check(std::abs(f.get_integral_displacement().z() - (S0.z() + p0.z()*h)) <= T(1e-12),
          "S did not integrate p");
    // v gains the OU-weighted acceleration, not a naive a h.
    check(f.get_velocity().x() > v0.x(), "v did not gain from a positive a_w");

    // Attitude must be untouched by the world chain.
    check(std::abs(f.quaternion().angularDistance(Eigen::Quaternion<T>::Identity())) <= T(1e-15),
          "the world chain rotated the attitude");
}

// Rotation direction: R is body-to-world and the gyro is body-frame, so a
// positive body-z rate must carry body +x toward world +y.
void test_rotation_direction() {
    Filter f;
    f.initialize_from_truth(Eigen::Quaternion<T>::Identity(),
                            Vector3::Zero(), Vector3::Zero(),
                            Vector3::Zero(), Vector3::Zero(),
                            Vector3::Zero(), Vector3::Zero());
    const T rate = T(M_PI) / T(2);
    f.propagate_mean(Vector3(T(0), T(0), rate), T(1));
    const Vector3 mapped = f.R_bw() * Vector3::UnitX();
    check((mapped - Vector3::UnitY()).cwiseAbs().maxCoeff() <= T(1e-12),
          "a +90 deg/s body-z rate for 1 s must send body x to world y");
}

// Gyro bias must be subtracted, not added.
void test_gyro_bias_is_subtracted() {
    Filter f;
    const Vector3 bias(T(0.1), T(-0.2), T(0.3));
    f.initialize_from_truth(Eigen::Quaternion<T>::Identity(),
                            Vector3::Zero(), Vector3::Zero(),
                            Vector3::Zero(), Vector3::Zero(),
                            bias, Vector3::Zero());
    // Measuring exactly the bias means the body is not turning.
    f.propagate_mean(bias, T(1.0));
    check(std::abs(f.quaternion().angularDistance(Eigen::Quaternion<T>::Identity())) <= T(1e-12),
          "a gyro reading equal to the bias must produce no rotation");
}

// Compile-only: the firmware instantiation must build and run.
void test_float_instantiation() {
    using FloatFilter = ocean_imu::kalman::Kalman3D_Wave_TFG<float, true, true, false>;
    FloatFilter f;
    f.set_aw_time_constant(1.5f);
    f.set_aw_stationary_std(Eigen::Vector3f(0.5f, 0.5f, 0.3f));
    f.initialize_identity();
    for (int i = 0; i < 200; ++i) f.time_update(Eigen::Vector3f(0.01f, -0.02f, 0.03f), 0.005f);
    check(f.covariance_full().allFinite() && f.state().R.allFinite(),
          "the float instantiation went non-finite");

    using NoAccelBias = ocean_imu::kalman::Kalman3D_Wave_TFG<float, true, false, false>;
    static_assert(NoAccelBias::NX == 18, "dropping the accel bias must give 18 states");
    NoAccelBias g;
    g.initialize_identity();
    g.time_update(Eigen::Vector3f(0.01f, 0.0f, 0.0f), 0.005f);
    check(g.covariance_full().allFinite(), "the 18-state instantiation went non-finite");
}

} // namespace

int main() {
    test_inject_error_roundtrip();
    test_attitude_bias_block_is_exact_at_rest();
    test_transition_matches_finite_differences();
    test_average_rotation_factor_matters();
    test_process_noise_is_psd();
    test_covariance_stays_psd_and_symmetric();
    test_mean_propagation_follows_the_chain();
    test_rotation_direction();
    test_gyro_bias_is_subtracted();
    test_float_instantiation();

    if (failures != 0) {
        std::cerr << failures << " propagation check(s) failed\n";
        return 1;
    }
    std::cout << "tfg_propagation-test: all checks passed\n";
    return 0;
}
