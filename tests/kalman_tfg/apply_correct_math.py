#!/usr/bin/env python3
"""One-shot branch patch: make the TFG Riccati model match its stated invariant error.

This script is intentionally run in CI on agent/tfg-math-hardening.  It edits
only the checked-out branch, runs the mathematical tests, and the workflow
commits the result only if those tests pass.
"""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
KALMAN = ROOT / "src/kalman_tfg/Kalman3D_Wave_TFG.h"
FUSION = ROOT / "src/kalman_tfg/SeaStateFusionFilter_TFG.h"
SIM = ROOT / "tests/kalman_tfg/kalman_tfg-sim.cpp"
PROP_TEST = ROOT / "tests/kalman_tfg/tfg_propagation-test.cpp"
PROC_TEST = ROOT / "tests/kalman_tfg/tfg_process_covariance-test.cpp"
MAKEFILE = ROOT / "tests/kalman_tfg/Makefile"
VALIDATION = ROOT / ".github/workflows/tfg-validation.yml"
DOC = ROOT / "docs/tfg-math-hardening.md"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected exactly one literal match, found {n}")
    return text.replace(old, new, 1)


def sub_once(text: str, pattern: str, repl: str, label: str) -> str:
    out, n = re.subn(pattern, repl, text, count=1, flags=re.S)
    if n != 1:
        raise RuntimeError(f"{label}: expected exactly one regex match, found {n}")
    return out


kalman = KALMAN.read_text()
if "full structured Gauss--Legendre discretization of the invariant error" in kalman:
    print("TFG correct-math patch is already applied")
    raise SystemExit(0)

kalman = replace_once(
    kalman,
    """    The world OU chain is discretized with IntegratedOUChain. Gyro-bias process\n    noise during finite rotation is integrated with Gauss-Legendre quadrature\n    over C(s)=int_s^h R(t)dt; the old h^2/2 and h^3/3 expressions are only the\n    zero-rotation special case and are no longer labelled exact.\n""",
    """    The world OU chain is discretized with IntegratedOUChain. The deterministic\n    gyro-bias/world coupling and the gyro/gyro-bias process covariance use a\n    full structured Gauss--Legendre discretization of the invariant error\n    dynamics. This transports each noise source through attitude and every\n    world column instead of applying zero-rate or first-order leakage formulas.\n""",
    "header covariance description",
)

# A process stationary covariance is a model parameter, not the posterior error
# covariance.  Keep reset_aw_covariance_to_stationary() only as an explicit
# initialization operation; remove the periodic posterior-forcing helper.
kalman = sub_once(
    kalman,
    r"\n    bool synchronize_aw_covariance_to_stationary_congruent\(\) \{.*?\n    \}\n\n    void set_aw_time_constant",
    "\n    void set_aw_time_constant",
    "remove posterior aw covariance forcing",
)

new_transition = r'''    void build_transition(T Ts, const Vector3& gyro_meas, MatrixNX& Phi) const {
        Phi.setIdentity();
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

        const Matrix3 R0 = X_.R;
        const Vector3 omega = gyro_meas - X_.B.col(0);
        const Matrix3 R_end = R0 * ocean_imu::lie::Exp<T>(Vector3(omega * Ts));

        // beta_g is world-referred at Level 2 and additive body-referred at
        // Level 1.  Fh = int_0^h R(t) dt maps a physical body-bias error into
        // the right-invariant attitude error exactly for constant measured
        // angular rate over the sample.
        const Matrix3 Fh =
            R0 * (Ts * ocean_imu::lie::left_jacobian<T>(Vector3(omega * Ts)));
        const Matrix3 bg0_to_body =
            two_frame_bias ? Matrix3(R0.transpose()) : Matrix3::Identity();
        Phi.template block<3,3>(OFF_PHI, OFF_BG) = -Fh * bg0_to_body;
        if constexpr (two_frame_bias)
            Phi.template block<3,3>(OFF_BG, OFF_BG) = R_end * R0.transpose();

        if constexpr (with_accel_bias) {
            const T alpha_ba = acc_bias_updates_enabled_ ? std::exp(-Ts/tau_ba_) : T(1);
            if constexpr (two_frame_bias)
                Phi.template block<3,3>(OFF_BA, OFF_BA) =
                    alpha_ba * R_end * R0.transpose();
            else
                Phi.template block<3,3>(OFF_BA, OFF_BA) =
                    alpha_ba * Matrix3::Identity();
        }

        // Exact OU-chain state transition for the active linear block; when
        // it is gated off the physical world columns are frozen and the world
        // transition is identity.
        Eigen::Matrix<T,4,4> Phi_axis;
        world_transition_(Ts, Phi_axis);
        Phi.template block<12,12>(OFF_V, OFF_V).setZero();
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (Phi_axis(i,j) == T(0)) continue;
                for (int k = 0; k < 3; ++k)
                    Phi(OFF_V + 3*i + k, OFF_V + 3*j + k) = Phi_axis(i,j);
            }
        }

        // rho_dot = A rho - [x(t)]x R(t) delta_b_g.  Integrate the actual
        // time-varying nominal x(t) and R(t), rather than replacing them by
        // their start-of-step values.  D maps a physical body bias to the
        // final world-error columns; Level 2 converts beta_g(0) back to that
        // physical body-bias coordinate with R0^T.
        Eigen::Matrix<T,12,3> D;
        integrate_world_bias_impulse_(T(0), Ts, omega, D);
        Phi.template block<12,3>(OFF_V, OFF_BG) = -D * bg0_to_body;
    }

'''
kalman = sub_once(
    kalman,
    r"    void build_transition\(T Ts, const Vector3& gyro_meas, MatrixNX& Phi\) const \{.*?\n    \}\n\n    void build_process_noise",
    new_transition + "    void build_process_noise",
    "replace deterministic invariant transition",
)

new_process_noise = r'''    void build_process_noise(T Ts, const Vector3& gyro_meas, MatrixNX& Qd) const {
        Qd.setZero();
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

        const Matrix3 R0 = X_.R;
        const Vector3 omega = gyro_meas - X_.B.col(0);
        const Matrix3 R_end = R0 * ocean_imu::lie::Exp<T>(Vector3(omega * Ts));
        const Matrix3 Fh =
            R0 * (Ts * ocean_imu::lie::left_jacobian<T>(Vector3(omega * Ts)));

        // Five-point Gauss--Legendre quadrature of the *complete structured
        // impulse response* of the stated continuous invariant error model.
        // For gyro white noise an impulse at s enters phi through R(s) and all
        // world columns through [x(s)]x R(s), followed by the remaining OU
        // chain.  A gyro-bias RW impulse persists from s to h, therefore also
        // drives phi and every world column for the rest of the sample.
        static constexpr double nodes[5] = {
            -0.9061798459386640, -0.5384693101056831, 0.0,
             0.5384693101056831,  0.9061798459386640
        };
        static constexpr double weights[5] = {
            0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
            0.4786286704993665, 0.2369268850561891
        };

        for (int q = 0; q < 5; ++q) {
            const T s = T(0.5) * Ts * (T(1) + T(nodes[q]));
            const T wq = T(0.5) * Ts * T(weights[q]);
            const Matrix3 R_s = R0 * ocean_imu::lie::Exp<T>(Vector3(omega * s));

            Eigen::Matrix<T,4,4> Phi_s, Phi_rem;
            world_transition_(s, Phi_s);
            world_transition_(Ts - s, Phi_rem);
            const Eigen::Matrix<T,3,4> X_s = X_.X * Phi_s.transpose();

            Eigen::Matrix<T,NX,3> Lg = Eigen::Matrix<T,NX,3>::Zero();
            Lg.template block<3,3>(OFF_PHI, 0) = R_s;
            for (int i = 0; i < 4; ++i) {
                Matrix3 M = Matrix3::Zero();
                for (int j = 0; j < 4; ++j) {
                    if (Phi_rem(i,j) == T(0)) continue;
                    M.noalias() += Phi_rem(i,j) *
                        (ocean_imu::lie::skew<T>(Vector3(X_s.col(j))) * R_s);
                }
                Lg.template block<3,3>(OFF_V + 3*i, 0) = M;
            }
            Qd.noalias() += wq * (Lg * Q_gyro_ * Lg.transpose());

            Eigen::Matrix<T,NX,3> Lb = Eigen::Matrix<T,NX,3>::Zero();
            const Matrix3 Fs =
                R0 * (s * ocean_imu::lie::left_jacobian<T>(Vector3(omega * s)));
            const Matrix3 C = Fh - Fs;  // int_s^h R(t) dt
            Lb.template block<3,3>(OFF_PHI, 0) = -C;
            Lb.template block<3,3>(OFF_BG, 0) =
                two_frame_bias ? R_end : Matrix3::Identity();

            Eigen::Matrix<T,12,3> D;
            integrate_world_bias_impulse_(s, Ts, omega, D);
            Lb.template block<12,3>(OFF_V, 0) = -D;
            Qd.noalias() += wq * (Lb * Q_bg_ * Lb.transpose());
        }

        // The OU acceleration driving noise is already available in exact
        // discrete form for the linear chain and is independent of gyro and
        // gyro-bias driving noise.
        if (linear_block_enabled_) add_ou_process_noise_(Ts, Qd);

        if constexpr (with_accel_bias) {
            if (acc_bias_updates_enabled_) {
                // A body-frame OU noise impulse decays in body coordinates;
                // its final Level-2 coordinate is simply R_end times that
                // decayed body increment, so the scalar OU integral is exact.
                const T qint = T(0.5) * tau_ba_ * (-std::expm1(T(-2)*Ts/tau_ba_));
                Qd.template block<3,3>(OFF_BA, OFF_BA) +=
                    (two_frame_bias ? Matrix3(R_end * Q_ba_ * R_end.transpose()) : Q_ba_) * qint;
            }
        }

        Qd = T(0.5) * (Qd + Qd.transpose()).eval();
        ou_detail::project_psd_ou_iii<T,NX>(Qd);
    }

'''
kalman = sub_once(
    kalman,
    r"    void build_process_noise\(T Ts, const Vector3& gyro_meas, MatrixNX& Qd\) const \{.*?\n    \}\n\n    \[\[nodiscard\]\] Tangent error_from",
    new_process_noise + "    [[nodiscard]] Tangent error_from",
    "replace complete process covariance",
)

helpers = r'''    void world_transition_(T dt, Eigen::Matrix<T,4,4>& Phi) const {
        Phi.setIdentity();
        if (linear_block_enabled_)
            ou_detail::IntegratedOUChain<T,3>::transition(tau_aw_, dt, Phi);
    }

    // D(s,h) = int_s^h Phi_A(h-t) [x(t)]x R(t) dt.  It is the physical
    // body-gyro-bias impulse map into [rho_v rho_p rho_S rho_a] at h.  Keeping
    // it in body-bias coordinates makes the same function valid for Level 1
    // and Level 2; only the endpoint bias coordinate differs.
    void integrate_world_bias_impulse_(T s, T h, const Vector3& omega,
                                       Eigen::Matrix<T,12,3>& D) const {
        D.setZero();
        if (!(h > s)) return;
        static constexpr double nodes[5] = {
            -0.9061798459386640, -0.5384693101056831, 0.0,
             0.5384693101056831,  0.9061798459386640
        };
        static constexpr double weights[5] = {
            0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
            0.4786286704993665, 0.2369268850561891
        };

        const T span = h - s;
        for (int q = 0; q < 5; ++q) {
            const T t = s + T(0.5) * span * (T(1) + T(nodes[q]));
            const T wq = T(0.5) * span * T(weights[q]);
            Eigen::Matrix<T,4,4> Phi_t, Phi_rem;
            world_transition_(t, Phi_t);
            world_transition_(h - t, Phi_rem);
            const Eigen::Matrix<T,3,4> X_t = X_.X * Phi_t.transpose();
            const Matrix3 R_t = X_.R * ocean_imu::lie::Exp<T>(Vector3(omega * t));

            for (int i = 0; i < 4; ++i) {
                Matrix3 M = Matrix3::Zero();
                for (int j = 0; j < 4; ++j) {
                    if (Phi_rem(i,j) == T(0)) continue;
                    M.noalias() += Phi_rem(i,j) *
                        (ocean_imu::lie::skew<T>(Vector3(X_t.col(j))) * R_t);
                }
                D.template block<3,3>(3*i, 0).noalias() += wq * M;
            }
        }
    }

'''
kalman = replace_once(
    kalman,
    "  private:\n    [[nodiscard]] Vector3 accel_bias_at_(T tempC) const {",
    "  private:\n" + helpers + "    [[nodiscard]] Vector3 accel_bias_at_(T tempC) const {",
    "insert invariant discretization helpers",
)

# A hard bias clamp changes the random variable without conditioning or
# transforming P.  Remove it from the EKF update; bounds, if wanted, belong in
# an explicitly constrained estimator rather than an invisible projection.
kalman = replace_once(
    kalman,
    """        const bool ok = apply_update3_(r, H, Rw, last_acc_diag_);\n        if (ok) project_acc_bias_();\n        return ok;\n""",
    """        return apply_update3_(r, H, Rw, last_acc_diag_);\n""",
    "remove unmodelled accelerometer-bias projection",
)
kalman = sub_once(
    kalman,
    r"\n    void project_acc_bias_\(\) \{.*?\n    \}\n\n    void reorthonormalize_",
    "\n    void reorthonormalize_",
    "remove bias projection implementation",
)
kalman = replace_once(kalman, "    T acc_bias_limit_{T(0.5)};\n", "", "remove bias clamp state")
KALMAN.write_text(kalman)

fusion = FUSION.read_text()
fusion = replace_once(
    fusion,
    """        } else {\n            adaptMekf_(dt);\n            aw_sync_elapsed_ += dt;\n            if (periodic_aw_cov_sync_ && aw_sync_elapsed_ >= aw_sync_period_sec_) {\n                aw_sync_elapsed_ = 0.0f;\n                mekf_.synchronize_aw_covariance_to_stationary_congruent();\n            }\n        }\n""",
    """        } else {\n            // Adapt the process model (tau and Sigma_aw), not the posterior\n            // covariance. P remains the Riccati covariance of the invariant\n            // error and is changed only by propagation, measurements, reset\n            // coordinate transport, or explicit initialization.\n            adaptMekf_(dt);\n        }\n""",
    "remove periodic posterior covariance forcing",
)
fusion = replace_once(
    fusion,
    "    void setPeriodicAwCovarianceSync(bool on) { periodic_aw_cov_sync_ = on; }\n",
    "",
    "remove covariance-sync API",
)
fusion = replace_once(fusion, "    bool periodic_aw_cov_sync_ = true;\n", "", "remove covariance-sync flag")
fusion = replace_once(fusion, "    float aw_sync_elapsed_ = 0.0f;\n", "", "remove covariance-sync elapsed")
fusion = replace_once(fusion, "    float aw_sync_period_sec_ = 1.0f;\n", "", "remove covariance-sync period")
FUSION.write_text(fusion)

sim = SIM.read_text()
sim = sub_once(
    sim,
    r"\n        if \(const char\* s = std::getenv\(\"TFG_AW_COV_SYNC\"\)\) \{\n            fusion_\.setPeriodicAwCovarianceSync\(std::string\(s\) != \"off\"\);\n        \}\n",
    "\n",
    "remove obsolete covariance-sync simulation knob",
)
SIM.write_text(sim)

prop = PROP_TEST.read_text()
prop = replace_once(
    prop,
    """            // The one documented approximation is in Phi_world,bg, which\n            // neglects the correlation of Rhat(s) with the Phi_OU weighting\n            // and the drift of x(s) across the step. Both are O(h) on a term\n            // already of size |x| |w| h, so the tolerance scales with h^2.\n            const T tol = T(2e-8) + T(4.0) * h * h * (T(1) + g.norm());\n""",
    """            // Phi_world,bg is now the quadrature of the actual\n            // time-varying invariant dynamics. A fixed tight tolerance catches\n            // frame/sign mistakes instead of granting an O(h^2) modelling gap.\n            const T tol = T(5e-7) * (T(1) + g.norm());\n""",
    "tighten deterministic transition oracle",
)
PROP_TEST.write_text(prop)

process_test = r'''/*
    Independent high-rate oracle for the TFG gyro and gyro-bias process
    covariance.  Production uses structured 5-point Gauss--Legendre; this test
    uses a dense midpoint integration of the physical impulse responses.
*/

#define EIGEN_NON_ARDUINO
#include "kalman_tfg/Kalman3D_Wave_TFG.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
using T = double;
using Filter = ocean_imu::kalman::Kalman3D_Wave_TFG<T, true, false, false>;
using Vector3 = Eigen::Matrix<T,3,1>;
using Matrix3 = Eigen::Matrix<T,3,3>;
using MatrixNX = Filter::MatrixNX;
using WorldMap = Eigen::Matrix<T,12,3>;
namespace lie = ocean_imu::lie;
namespace ou = ocean_imu::kalman::ou_detail;

int failures = 0;

bool check(bool ok, const char* msg) {
    if (!ok) { std::cerr << "FAIL: " << msg << '\n'; ++failures; }
    return ok;
}

Eigen::Matrix<T,4,4> world_phi(T tau, T dt) {
    Eigen::Matrix<T,4,4> P;
    ou::IntegratedOUChain<T,3>::transition(tau, dt, P);
    return P;
}

void test_full_gyro_and_bias_covariance() {
    Filter f;
    const T tau = T(1.7);
    f.set_aw_time_constant(tau);
    f.set_aw_stationary_std(Vector3::Zero());
    f.set_gyro_noise_density_rad_sqrt_s(T(0.007));
    const Vector3 qbg(T(3e-5), T(8e-5), T(5e-5));
    f.set_Q_bgyro_rw(qbg);
    f.initialize_from_truth(
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(
            T(0.91), Vector3(T(0.2),T(-0.6),T(0.7)).normalized())),
        Vector3(T(0.8),T(-0.4),T(0.3)),
        Vector3(T(-1.2),T(2.1),T(0.7)),
        Vector3(T(3.0),T(-1.8),T(1.1)),
        Vector3(T(0.5),T(-0.7),T(0.9)),
        Vector3(T(0.012),T(-0.008),T(0.004)));

    const Vector3 gyro(T(1.4),T(-0.9),T(0.7));
    const Vector3 omega = gyro - f.gyroscope_bias();
    const Matrix3 R0 = f.R_bw();
    const auto X0 = f.state().X;
    const T h = T(0.05);

    MatrixNX got;
    f.build_process_noise(h, gyro, got);

    constexpr int N = 4000;
    const T ds = h / T(N);
    std::vector<WorldMap> M(static_cast<std::size_t>(N));
    std::vector<WorldMap> D(static_cast<std::size_t>(N));

    // M(t) is the body-bias -> final-world-error integrand.  Build its
    // backward tail integral so D(s)=int_s^h M(t)dt is O(N), not O(N^2).
    for (int k = 0; k < N; ++k) {
        const T t = (T(k) + T(0.5)) * ds;
        const auto Pt = world_phi(tau, t);
        const auto Prem = world_phi(tau, h - t);
        const Eigen::Matrix<T,3,4> Xt = X0 * Pt.transpose();
        const Matrix3 Rt = R0 * lie::Exp<T>(Vector3(omega*t));
        M[static_cast<std::size_t>(k)].setZero();
        for (int i = 0; i < 4; ++i) {
            Matrix3 B = Matrix3::Zero();
            for (int j = 0; j < 4; ++j)
                B.noalias() += Prem(i,j) * (lie::skew<T>(Vector3(Xt.col(j))) * Rt);
            M[static_cast<std::size_t>(k)].template block<3,3>(3*i,0) = B;
        }
    }
    WorldMap tail = WorldMap::Zero();
    for (int k = N-1; k >= 0; --k) {
        const auto idx = static_cast<std::size_t>(k);
        D[idx] = tail + T(0.5)*ds*M[idx];
        tail.noalias() += ds*M[idx];
    }

    const Matrix3 Qg = Matrix3::Identity() * T(0.007)*T(0.007);
    const Matrix3 Qb = qbg.asDiagonal();
    const Matrix3 Fh = R0 * (h * lie::left_jacobian<T>(Vector3(omega*h)));
    MatrixNX want = MatrixNX::Zero();

    for (int k = 0; k < N; ++k) {
        const T s = (T(k) + T(0.5))*ds;
        const auto Ps = world_phi(tau, s);
        const auto Prem = world_phi(tau, h-s);
        const Eigen::Matrix<T,3,4> Xs = X0 * Ps.transpose();
        const Matrix3 Rs = R0 * lie::Exp<T>(Vector3(omega*s));

        Eigen::Matrix<T,Filter::NX,3> Lg = Eigen::Matrix<T,Filter::NX,3>::Zero();
        Lg.template block<3,3>(Filter::OFF_PHI,0) = Rs;
        for (int i = 0; i < 4; ++i) {
            Matrix3 B = Matrix3::Zero();
            for (int j = 0; j < 4; ++j)
                B.noalias() += Prem(i,j) * (lie::skew<T>(Vector3(Xs.col(j))) * Rs);
            Lg.template block<3,3>(Filter::OFF_V+3*i,0) = B;
        }

        const Matrix3 Fs = R0 * (s * lie::left_jacobian<T>(Vector3(omega*s)));
        Eigen::Matrix<T,Filter::NX,3> Lb = Eigen::Matrix<T,Filter::NX,3>::Zero();
        Lb.template block<3,3>(Filter::OFF_PHI,0) = -(Fh-Fs);
        Lb.template block<3,3>(Filter::OFF_BG,0) = Matrix3::Identity();
        Lb.template block<12,3>(Filter::OFF_V,0) = -D[static_cast<std::size_t>(k)];

        want.noalias() += ds * (Lg*Qg*Lg.transpose() + Lb*Qb*Lb.transpose());
    }
    want = T(0.5)*(want+want.transpose()).eval();

    const T scale = std::max(T(1e-12), want.cwiseAbs().maxCoeff());
    const T err = (got-want).cwiseAbs().maxCoeff();
    if (!check(err <= T(3e-4)*scale + T(2e-11),
               "structured Qd disagrees with high-rate physical oracle"))
        std::cerr << "  max|Q-Qoracle|=" << err << " scale=" << scale << '\n';

    const T rho_bg = got.template block<12,3>(Filter::OFF_V,Filter::OFF_BG).norm();
    check(rho_bg > T(1e-9),
          "gyro-bias driving noise did not reach the world/bias cross covariance");
    const T rho_rho = got.template block<12,12>(Filter::OFF_V,Filter::OFF_V).norm();
    check(rho_rho > T(1e-9),
          "gyro/gyro-bias driving noise did not reach world covariance");
}
} // namespace

int main() {
    test_full_gyro_and_bias_covariance();
    if (failures) {
        std::cerr << failures << " TFG process-covariance check(s) failed\n";
        return 1;
    }
    std::cout << "tfg_process_covariance-test: all checks passed\n";
    return 0;
}
'''
PROC_TEST.write_text(process_test)

mk = MAKEFILE.read_text()
mk = replace_once(
    mk,
    "covariance_transport-test tfg_orchestrator-test kalman_tfg-sim",
    "covariance_transport-test tfg_process_covariance-test tfg_orchestrator-test kalman_tfg-sim",
    "add process covariance target",
)
mk = replace_once(
    mk,
    """covariance_transport-test: covariance_transport-test.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\ntfg_orchestrator-test:""",
    """covariance_transport-test: covariance_transport-test.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\ntfg_process_covariance-test: tfg_process_covariance-test.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\ntfg_orchestrator-test:""",
    "add process covariance make rule",
)
MAKEFILE.write_text(mk)

validation = VALIDATION.read_text()
validation = replace_once(
    validation,
    "          ./covariance_transport-test\n          ./tfg_orchestrator-test",
    "          ./covariance_transport-test\n          ./tfg_process_covariance-test\n          ./tfg_orchestrator-test",
    "run process covariance test in validation",
)
VALIDATION.write_text(validation)

doc = DOC.read_text()
doc = replace_once(
    doc,
    """The gyro-white-noise leakage into world-vector blocks remains the documented first-order approximation inherited from the surrounding estimator; this hardening does not claim that the entire `Q_d` is a closed-form exact discretization.\n""",
    """The gyro-white and gyro-bias driving noises are now discretized from their complete structured impulse responses. Each quadrature point transports gyro noise through attitude and the remaining OU world-chain transition, while gyro-bias random-walk noise is integrated through its persistent effects on attitude and every world-error column. This removes the former first-order gyro/world leakage approximation and the missing gyro-bias/world process-covariance blocks. The five-point Gauss--Legendre rule is a high-order numerical discretization of the stated continuous model rather than a claim of a symbolic closed form.\n""",
    "update Qd documentation",
)
doc += """\n## 5. Posterior covariance is no longer forced to the OU stationary covariance\n\n`Sigma_aw` is a process-model parameter used to construct the OU driving covariance. It is not the posterior covariance of the invariant error `rho_aw`, which also contains attitude coupling and measurement information. The TFG orchestrator therefore no longer periodically rescales `P_aw,aw` and its cross-covariances to `Sigma_aw`. `P` now evolves only through the Riccati propagation, measurement updates, mathematically required reset-coordinate transport, and explicit initialization.\n\nThe accelerometer-bias norm projection was also removed from the EKF update because changing the nominal random variable without conditioning or transforming `P` breaks the posterior covariance interpretation.\n"""
DOC.write_text(doc)

print("Applied mathematically consistent TFG covariance propagation and tests")
