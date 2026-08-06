/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    The OU chain must be *shared* with OU-III, not re-derived.

    The whole argument for this filter's discretization is that the
    right-invariant error of the world vectors obeys the same linear chain the
    existing article already derives and validates:

        rho_vdot = rho_a,  rho_pdot = rho_v,  rho_Sdot = rho_p,
        rho_adot = -lambda rho_a

    If that is true, the 12x12 world block of Phi and Qd is
    ou_detail::IntegratedOUChain<T,3> unchanged -- the analytic coefficients,
    the small-h/tau Maclaurin branch, the correlated-axis assembly, and the
    PSD handling all carry over, along with the evidence already gathered for
    them.

    So this test asserts equality, not agreement-to-a-tolerance. Anything
    looser would let the world block quietly drift into a private
    reimplementation, which is exactly the outcome the design is meant to
    avoid.
*/

#define EIGEN_NON_ARDUINO

#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "kalman_tfg/Kalman3D_Wave_TFG.h"

#include <array>
#include <cmath>
#include <iostream>
#include <string>

namespace detail = ocean_imu::kalman::ou_detail;

namespace {

using T = double;
using Filter = ocean_imu::kalman::Kalman3D_Wave_TFG<T, true, true, false>;
using Vector3 = Eigen::Matrix<T,3,1>;
using Matrix3 = Eigen::Matrix<T,3,3>;

int failures = 0;

bool check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
    return condition;
}

// Operating points spanning the Maclaurin branch boundary at h/tau = 1e-2 and
// the ordinary regime, mirroring tests/kalman_ou_iii/kalman_ou_common-test.
const std::array<T,4> kTaus   = {T(0.05), T(0.5), T(2.1), T(10)};
const std::array<T,6> kRatios = {T(1e-5), T(1e-3), T(0.0099), T(0.0101), T(0.1), T(0.8)};

Filter make_filter(T tau, const Matrix3& Sigma_aw) {
    Filter f;
    f.set_aw_time_constant(tau);
    f.set_aw_stationary_cov_full(Sigma_aw);
    // A non-identity attitude and non-zero world vectors, so that a world
    // block accidentally contaminated by the attitude or bias coupling shows
    // up here rather than passing by being multiplied by zero.
    f.initialize_from_truth(
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(T(0.7), Vector3(T(0.3), T(-0.5), T(0.81)).normalized())),
        Vector3(T(0.4), T(-0.9), T(1.3)),      // v
        Vector3(T(-2.0), T(3.1), T(0.7)),      // p
        Vector3(T(5.0), T(-1.5), T(2.2)),      // S
        Vector3(T(0.2), T(0.4), T(-0.6)),      // a_w
        Vector3(T(0.01), T(-0.02), T(0.005)),  // b_g
        Vector3(T(0.03), T(0.01), T(-0.02)));  // b_a
    return f;
}

void test_transition_block_is_identical() {
    for (T tau : kTaus) {
        for (T ratio : kRatios) {
            const T h = tau * ratio;
            Filter f = make_filter(tau, Matrix3::Identity());

            Filter::MatrixNX Phi;
            f.build_transition(h, Vector3(T(0.1), T(-0.2), T(0.05)), Phi);

            Eigen::Matrix<T,4,4> Phi_axis;
            detail::IntegratedOUChain<T,3>::transition(tau, h, Phi_axis);

            const std::string at =
                " (tau=" + std::to_string(tau) + ", h/tau=" + std::to_string(ratio) + ")";

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        const T got = Phi(Filter::OFF_V + 3*i + k, Filter::OFF_V + 3*j + k);
                        // Bit-for-bit: the filter must be reading these out of
                        // IntegratedOUChain, not recomputing them.
                        if (!check(got == Phi_axis(i,j),
                                   "world transition entry (" + std::to_string(i) + "," +
                                   std::to_string(j) + ") axis " + std::to_string(k) +
                                   " differs from IntegratedOUChain" + at)) {
                            std::cerr << "  got " << got << " want " << Phi_axis(i,j) << '\n';
                        }
                    }
                    // Off-axis entries must be exactly zero: the chain acts
                    // per axis and must not couple x into y.
                    for (int k = 0; k < 3; ++k) {
                        for (int m = 0; m < 3; ++m) {
                            if (k == m) continue;
                            check(Phi(Filter::OFF_V + 3*i + k, Filter::OFF_V + 3*j + m) == T(0),
                                  "world transition couples axes" + at);
                        }
                    }
                }
            }
        }
    }
}

void test_process_noise_block_is_identical() {
    for (T tau : kTaus) {
        for (T ratio : kRatios) {
            const T h = tau * ratio;
            const Vector3 sigma(T(0.7), T(1.3), T(0.4));
            const Matrix3 Sigma = sigma.cwiseAbs2().asDiagonal();

            Filter f = make_filter(tau, Sigma);
            // Isolate the OU contribution: with no gyro white noise and no
            // bias random walk, the world block is the OU term alone.
            f.set_gyro_noise_density_rad_sqrt_s(T(0));
            f.set_Q_bgyro_rw(Vector3::Zero());
            f.set_Q_bacc_rw(Vector3::Zero());

            Filter::MatrixNX Qd;
            f.build_process_noise(h, Qd);

            const std::string at =
                " (tau=" + std::to_string(tau) + ", h/tau=" + std::to_string(ratio) + ")";

            for (int axis = 0; axis < 3; ++axis) {
                Eigen::Matrix<T,4,4> Q_axis;
                detail::IntegratedOUChain<T,3>::process_covariance(
                    tau, h, Sigma(axis,axis), Q_axis);
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        const T got = Qd(Filter::OFF_V + 3*i + axis,
                                         Filter::OFF_V + 3*j + axis);
                        const T want = Q_axis(i,j);
                        // project_psd_ou_iii runs over the assembled matrix and
                        // may nudge entries, so this one is a tight tolerance
                        // rather than exact equality -- but tight enough that a
                        // different formula cannot hide inside it.
                        const T tol = T(1e-13) * std::max(T(1), std::abs(want));
                        if (!check(std::abs(got - want) <= tol,
                                   "world process-noise entry (" + std::to_string(i) + "," +
                                   std::to_string(j) + ") axis " + std::to_string(axis) +
                                   " differs from IntegratedOUChain" + at)) {
                            std::cerr << "  got " << got << " want " << want
                                      << " diff " << (got - want) << '\n';
                        }
                    }
                }
            }
        }
    }
}

// The correlated-axis path must reproduce OU-III's Kronecker assembly,
// Q = Sigma_aw (x) Q_unit, in the same group-first ordering.
void test_correlated_axis_assembly() {
    Matrix3 Sigma;
    const T s = T(0.9);
    Sigma << s*s,            T(-0.65)*s*s*T(0.8), T(-0.65)*s*s*T(0.7),
             T(-0.65)*s*s*T(0.8), s*s*T(1.4),     T(0.2)*s*s,
             T(-0.65)*s*s*T(0.7), T(0.2)*s*s,     s*s*T(0.6);
    Sigma = T(0.5) * (Sigma + Sigma.transpose()).eval();

    for (T tau : kTaus) {
        const T h = tau * T(0.1);
        Filter f = make_filter(tau, Sigma);
        f.set_gyro_noise_density_rad_sqrt_s(T(0));
        f.set_Q_bgyro_rw(Vector3::Zero());
        f.set_Q_bacc_rw(Vector3::Zero());

        Filter::MatrixNX Qd;
        f.build_process_noise(h, Qd);

        Eigen::Matrix<T,4,4> Q_unit;
        detail::IntegratedOUChain<T,3>::process_covariance(tau, h, T(1), Q_unit);
        Q_unit = T(0.5) * (Q_unit + Q_unit.transpose()).eval();

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                const Matrix3 want = Sigma * Q_unit(i,j);
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        const T got = Qd(Filter::OFF_V + 3*i + a, Filter::OFF_V + 3*j + b);
                        const T tol = T(1e-13) * std::max(T(1), std::abs(want(a,b)));
                        check(std::abs(got - want(a,b)) <= tol,
                              "correlated-axis world process noise does not match "
                              "Sigma (x) Q_unit at tau=" + std::to_string(tau));
                    }
                }
            }
        }
    }
}

/*
    Psi(h,tau) = int_0^h Phi_OU(u) du, checked against Simpson quadrature of
    the transition it integrates. This is the term that makes Phi_world,bg
    exact rather than truncated, and it has three entries that cancel to
    O(x^2), O(x^3), O(x^4), so its Taylor branch is load-bearing.

    The quadrature reference is built from the closed form in long double
    rather than by calling IntegratedOUChain. That is deliberate. OU-III's
    small-x branch truncates phi_pa at x^4 and phi_Sa at x^5
    (KalmanOUCoreMath.h:86), so near the 1e-2 seam it carries a relative error
    of about x^3/60, roughly 2e-8. Integrating that would measure OU-III's
    truncation rather than this function -- and Psi is the more accurate of
    the two, so the test would fail pointing at the wrong file. The
    truncation itself is far below anything that matters for a transition
    coefficient and is not worth disturbing the versioned OU evidence to
    change; it is recorded here so the next person to see an 8e-9 discrepancy
    at the seam knows where it comes from.
*/
void test_psi_against_quadrature() {
    using LD = long double;

    // Exact Phi_OU for one axis, no Taylor branch, evaluated in long double.
    auto exact_phi = [](LD tau, LD u) {
        Eigen::Matrix<LD,4,4> Phi = Eigen::Matrix<LD,4,4>::Zero();
        const LD x = u / tau;
        const LD em1 = std::expm1(-x);
        const LD tau2 = tau * tau, tau3 = tau2 * tau;
        Phi(0,0) = LD(1);
        Phi(0,3) = -tau * em1;
        Phi(1,0) = u;
        Phi(1,1) = LD(1);
        Phi(1,3) = tau2 * (x + em1);
        Phi(2,0) = LD(0.5) * u * u;
        Phi(2,1) = u;
        Phi(2,2) = LD(1);
        Phi(2,3) = tau3 * (LD(0.5) * x * x - x - em1);
        Phi(3,3) = std::exp(-x);
        return Phi;
    };

    for (T tau : kTaus) {
        for (T ratio : kRatios) {
            const T h = tau * ratio;

            Eigen::Matrix<T,4,4> Psi;
            ocean_imu::kalman::tfg_detail::integral_transition_axis<T>(tau, h, Psi);

            const int n = 4000;
            Eigen::Matrix<LD,4,4> acc_ld = Eigen::Matrix<LD,4,4>::Zero();
            for (int k = 0; k <= n; ++k) {
                const LD u = static_cast<LD>(h) * static_cast<LD>(k) / static_cast<LD>(n);
                const LD w = (k == 0 || k == n) ? LD(1) : ((k % 2 == 1) ? LD(4) : LD(2));
                acc_ld += w * exact_phi(static_cast<LD>(tau), u);
            }
            acc_ld *= static_cast<LD>(h) / (LD(3) * static_cast<LD>(n));
            const Eigen::Matrix<T,4,4> acc = acc_ld.cast<T>();

            const std::string at =
                " (tau=" + std::to_string(tau) + ", h/tau=" + std::to_string(ratio) + ")";
            // Relative, because the entries span many orders of magnitude:
            // Psi_23 ~ tau^4 x^4/24 is astronomically smaller than Psi_00 = h.
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    const T scale = std::max(std::abs(acc(i,j)), T(1e-300));
                    const T rel = std::abs(Psi(i,j) - acc(i,j)) / scale;
                    if (!check(rel <= T(1e-9) || std::abs(Psi(i,j) - acc(i,j)) <= T(1e-18),
                               "Psi entry (" + std::to_string(i) + "," + std::to_string(j) +
                               ") disagrees with quadrature" + at)) {
                        std::cerr << "  got " << Psi(i,j) << " want " << acc(i,j)
                                  << " rel " << rel << '\n';
                    }
                }
            }
        }
    }
}

} // namespace

int main() {
    test_transition_block_is_identical();
    test_process_noise_block_is_identical();
    test_correlated_axis_assembly();
    test_psi_against_quadrature();

    if (failures != 0) {
        std::cerr << failures << " OU chain identity check(s) failed\n";
        return 1;
    }
    std::cout << "ou_chain_identity-test: all checks passed\n";
    return 0;
}
