#define EIGEN_NON_ARDUINO

#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

#include <unsupported/Eigen/MatrixFunctions>

#include <array>
#include <cmath>
#include <iostream>

namespace detail = ocean_imu::kalman::ou_detail;

namespace {

using T = double;
using Matrix3 = Eigen::Matrix<T,3,3>;
using Matrix4 = Eigen::Matrix<T,4,4>;
using Vector3 = Eigen::Matrix<T,3,1>;

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

template<int N>
Eigen::Matrix<T,N,N> van_loan_q(T tau, T dt, T sigma2) {
    Eigen::Matrix<T,N,N> F = Eigen::Matrix<T,N,N>::Zero();
    F(0,N-1) = T(1);
    F(1,0) = T(1);
    if constexpr (N == 4) F(2,1) = T(1);
    F(N-1,N-1) = -T(1) / tau;

    Eigen::Matrix<T,N,N> Qc = Eigen::Matrix<T,N,N>::Zero();
    Qc(N-1,N-1) = T(2) * sigma2 / tau;

    Eigen::Matrix<T,2*N,2*N> A = Eigen::Matrix<T,2*N,2*N>::Zero();
    A.template topLeftCorner<N,N>() = F;
    A.template topRightCorner<N,N>() = Qc;
    A.template bottomRightCorner<N,N>() = -F.transpose();
    const Eigen::Matrix<T,2*N,2*N> E = (A * dt).exp();
    const Eigen::Matrix<T,N,N> Phi = E.template topLeftCorner<N,N>();
    return E.template topRightCorner<N,N>() * Phi.transpose();
}

bool test_ou_covariance() {
    const std::array<T,4> taus = {T(0.05), T(0.5), T(2.1), T(10)};
    const std::array<T,6> ratios = {T(1e-5), T(1e-3), T(0.0099), T(0.0101), T(0.1), T(0.8)};
    const int map[3] = {0,1,3};

    for (T tau : taus) {
        for (T ratio : ratios) {
            const T dt = tau * ratio;
            const T sigma2 = T(0.83);
            Matrix3 q2;
            Matrix4 q3;
            detail::IntegratedOUChain<T,2>::process_covariance(tau, dt, sigma2, q2);
            detail::IntegratedOUChain<T,3>::process_covariance(tau, dt, sigma2, q3);
            const Matrix3 ref2 = van_loan_q<3>(tau, dt, sigma2);
            const Matrix4 ref3 = van_loan_q<4>(tau, dt, sigma2);
            const T scale2 = std::max(T(1), ref2.cwiseAbs().maxCoeff());
            const T scale3 = std::max(T(1), ref3.cwiseAbs().maxCoeff());
            if (!check((q2-ref2).cwiseAbs().maxCoeff() <= T(5e-10)*scale2, "OU-II Qd mismatch")) return false;
            if (!check((q3-ref3).cwiseAbs().maxCoeff() <= T(5e-10)*scale3, "OU-III Qd mismatch")) return false;
            for (int i=0; i<3; ++i) {
                for (int j=0; j<3; ++j) {
                    if (!check(std::abs(q2(i,j)-q3(map[i],map[j])) <= T(2e-13)*scale2,
                               "OU-II/OU-III marginal mismatch")) return false;
                }
            }
            Eigen::SelfAdjointEigenSolver<Matrix3> e2(q2);
            Eigen::SelfAdjointEigenSolver<Matrix4> e3(q3);
            if (!check(e2.eigenvalues().minCoeff() >= -T(1e-12)*scale2, "OU-II Qd is not PSD")) return false;
            if (!check(e3.eigenvalues().minCoeff() >= -T(1e-12)*scale3, "OU-III Qd is not PSD")) return false;
        }
    }
    return true;
}

template<typename Filter>
bool covariance_unchanged(Filter& filter, const Vector3& std_aw) {
    const auto before = filter.covariance_full();
    filter.set_aw_stationary_std(std_aw);
    return (filter.covariance_full() - before).cwiseAbs().maxCoeff() == T(0);
}

bool test_stationary_setters() {
    const Vector3 sigma_a = Vector3::Constant(T(0.02));
    const Vector3 gyro_density = Vector3::Constant(T(0.001));
    const Vector3 sigma_m = Vector3::Constant(T(0.5));
    const Vector3 zero = Vector3::Zero();
    const Vector3 std_aw(T(0.4), T(0.5), T(0.6));
    Matrix3 full;
    full << T(0.20), T(0.02), T(-0.01),
            T(0.02), T(0.30), T(0.015),
            T(-0.01), T(0.015), T(0.40);

    Kalman3D_Wave_OU_II<T> ou2(sigma_a, gyro_density, sigma_m);
    Kalman3D_Wave_OU_III<T> ou3(sigma_a, gyro_density, sigma_m);
    for (int i=0; i<8; ++i) {
        ou2.time_update(zero, T(0.005));
        ou3.time_update(zero, T(0.005));
    }

    if (!check(covariance_unchanged(ou2, std_aw), "OU-II std setter changed posterior covariance")) return false;
    if (!check(covariance_unchanged(ou3, std_aw), "OU-III std setter changed posterior covariance")) return false;

    auto p2 = ou2.covariance_full();
    auto p3 = ou3.covariance_full();
    ou2.set_aw_stationary_corr_std(std_aw, T(-0.2), T(0.15));
    ou3.set_aw_stationary_corr_std(std_aw, T(-0.2), T(0.15));
    if (!check((ou2.covariance_full()-p2).cwiseAbs().maxCoeff() == T(0), "OU-II corr setter changed posterior covariance")) return false;
    if (!check((ou3.covariance_full()-p3).cwiseAbs().maxCoeff() == T(0), "OU-III corr setter changed posterior covariance")) return false;

    p2 = ou2.covariance_full();
    p3 = ou3.covariance_full();
    ou2.set_aw_stationary_cov_full(full);
    ou3.set_aw_stationary_cov_full(full);
    if (!check((ou2.covariance_full()-p2).cwiseAbs().maxCoeff() == T(0), "OU-II full setter changed posterior covariance")) return false;
    if (!check((ou3.covariance_full()-p3).cwiseAbs().maxCoeff() == T(0), "OU-III full setter changed posterior covariance")) return false;

    ou2.reset_aw_covariance_to_stationary();
    ou3.reset_aw_covariance_to_stationary();
    const auto reset2 = ou2.covariance_full();
    const auto reset3 = ou3.covariance_full();
    constexpr int aw2 = 12;
    constexpr int aw3 = 15;
    if (!check((reset2.template block<3,3>(aw2,aw2)-full).cwiseAbs().maxCoeff() < T(1e-13), "OU-II reset covariance mismatch")) return false;
    if (!check((reset3.template block<3,3>(aw3,aw3)-full).cwiseAbs().maxCoeff() < T(1e-13), "OU-III reset covariance mismatch")) return false;
    for (int i=0; i<reset2.rows(); ++i) {
        if (i < aw2 || i >= aw2+3) {
            if (!check(reset2.template block<1,3>(i,aw2).cwiseAbs().maxCoeff() == T(0), "OU-II reset retained cross covariance")) return false;
        }
    }
    for (int i=0; i<reset3.rows(); ++i) {
        if (i < aw3 || i >= aw3+3) {
            if (!check(reset3.template block<1,3>(i,aw3).cwiseAbs().maxCoeff() == T(0), "OU-III reset retained cross covariance")) return false;
        }
    }
    return true;
}

bool test_noise_units_and_period() {
    const Vector3 sample_std = Vector3::Constant(T(0.02));
    const T dt = T(0.005);
    const Vector3 expected = sample_std * std::sqrt(dt);
    const Vector3 converted = Kalman3D_Wave_OU_II<T>::gyro_noise_density_from_sample_std(sample_std, dt);
    if (!check((converted-expected).cwiseAbs().maxCoeff() < T(1e-15), "gyro noise conversion mismatch")) return false;

    for (int rate : {100, 200, 400}) {
        T elapsed = T(0);
        int count = 0;
        const T step = T(1) / T(rate);
        for (int i=0; i<rate*3; ++i) {
            if (detail::periodic_update_due(step, T(0.015), elapsed)) ++count;
        }
        if (!check(count == 200, "pseudo update rate depends on sample rate")) return false;
        if (!check(elapsed >= T(0) && elapsed < T(0.015), "pseudo update remainder is out of range")) return false;
    }
    for (int rate : {100, 200, 400}) {
        float elapsed = 0.0f;
        int count = 0;
        const float step = 1.0f / static_cast<float>(rate);
        for (int sample=0; sample<rate*3; ++sample) {
            if (detail::periodic_update_due(step, 0.015f, elapsed)) ++count;
        }
        if (!check(count == 200, "float pseudo update rate depends on sample rate")) return false;
        if (!check(elapsed >= 0.0f && elapsed < 0.015f, "float pseudo update remainder is out of range")) return false;
    }
    return true;
}

bool test_attitude_helpers() {
    Vector3 rotation_vector;
    rotation_vector << T(0.01), T(-0.02), T(0.03);
    const auto q = detail::quat_from_delta_theta(rotation_vector);
    if (!check(std::abs(q.norm()-T(1)) <= T(1e-14), "quaternion increment is not unit length")) return false;
    Matrix3 R, B;
    detail::rot_and_B_from_wt(rotation_vector, T(0.005), R, B);
    return check((R*R.transpose()-Matrix3::Identity()).cwiseAbs().maxCoeff() <= T(1e-13),
                 "rotation transition is not orthogonal");
}

} // namespace

int main() {
    if (!test_ou_covariance()) return 1;
    if (!test_stationary_setters()) return 1;
    if (!test_noise_units_and_period()) return 1;
    if (!test_attitude_helpers()) return 1;
    return 0;
}
