/*
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
