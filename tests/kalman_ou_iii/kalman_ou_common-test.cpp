#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

#include <cmath>
#include <iostream>

namespace detail = ocean_imu::kalman::ou_detail;

int main() {
    using T = double;
    using Matrix3 = Eigen::Matrix<T,3,3>;
    using Matrix4 = Eigen::Matrix<T,4,4>;
    using Vector3 = Eigen::Matrix<T,3,1>;

    constexpr T tau = 1.7;
    constexpr T dt = 0.005;
    constexpr T sigma2 = 0.83;

    Matrix3 phi2, q2;
    Matrix4 phi3, q3;
    detail::IntegratedOUChain<T,2>::transition(tau, dt, phi2);
    detail::IntegratedOUChain<T,2>::process_covariance(tau, dt, sigma2, q2);
    detail::IntegratedOUChain<T,3>::transition(tau, dt, phi3);
    detail::IntegratedOUChain<T,3>::process_covariance(tau, dt, sigma2, q3);

    const int idx3[3] = {0, 1, 3};
    T max_phi_error = T(0);
    T max_q_error = T(0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            max_phi_error = std::max(max_phi_error, std::abs(phi2(i,j) - phi3(idx3[i], idx3[j])));
            max_q_error = std::max(max_q_error, std::abs(q2(i,j) - q3(idx3[i], idx3[j])));
        }
    }
    if (max_phi_error > T(1e-14) || max_q_error > T(1e-10)) {
        std::cerr << "OU-II/OU-III marginal mismatch: phi=" << max_phi_error
                  << " q=" << max_q_error << '\n';
        return 1;
    }

    Vector3 rotation_vector;
    rotation_vector << T(0.01), T(-0.02), T(0.03);
    const auto q = detail::quat_from_delta_theta(rotation_vector);
    if (std::abs(q.norm() - T(1)) > T(1e-14)) {
        std::cerr << "Quaternion increment is not unit length\n";
        return 1;
    }

    Matrix3 R, B;
    detail::rot_and_B_from_wt(rotation_vector, dt, R, B);
    const Matrix3 orthogonality = R * R.transpose() - Matrix3::Identity();
    if (orthogonality.cwiseAbs().maxCoeff() > T(1e-13)) {
        std::cerr << "Rotation transition is not orthogonal\n";
        return 1;
    }

    // Force both public headers and common template code to instantiate in one TU.
    const Vector3 sigma_a = Vector3::Constant(T(0.02));
    const Vector3 sigma_g = Vector3::Constant(T(0.001));
    const Vector3 sigma_m = Vector3::Constant(T(0.5));
    Kalman3D_Wave_OU_II<T> ou2(sigma_a, sigma_g, sigma_m);
    Kalman3D_Wave_OU_III<T> ou3(sigma_a, sigma_g, sigma_m);
    (void)ou2;
    (void)ou3;

    return 0;
}
