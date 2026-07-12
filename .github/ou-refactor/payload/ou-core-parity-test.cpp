#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string_view>
#include <type_traits>

#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

namespace {

template<typename DerivedA, typename DerivedB>
bool near_matrix(const Eigen::MatrixBase<DerivedA>& a,
                 const Eigen::MatrixBase<DerivedB>& b,
                 double abs_tol,
                 double rel_tol,
                 std::string_view label) {
    const double scale = std::max(a.template cast<double>().cwiseAbs().maxCoeff(),
                                  b.template cast<double>().cwiseAbs().maxCoeff());
    const double error = (a.template cast<double>() - b.template cast<double>()).cwiseAbs().maxCoeff();
    if (error <= abs_tol + rel_tol * scale) {
        return true;
    }
    std::cerr << label << " mismatch: error=" << error << " scale=" << scale << '\n';
    return false;
}

template<typename T>
bool check_chain_parity() {
    using ChainII = ocean_imu::kalman::ou_detail::IntegratedOUChain<T,2>;
    using ChainIII = ocean_imu::kalman::ou_detail::IntegratedOUChain<T,3>;

    constexpr std::array<T,5> taus{T(0.05), T(0.5), T(1), T(2), T(10)};
    constexpr std::array<T,4> steps{T(0.001), T(0.005), T(0.02), T(0.1)};
    constexpr T sigma2 = T(4.84);
    constexpr int keep[3] = {0, 1, 3};

    const double abs_tol = std::is_same_v<T,float> ? 2e-7 : 2e-11;
    const double rel_tol = std::is_same_v<T,float> ? 2e-5 : 2e-10;

    for (const T tau : taus) {
        for (const T step : steps) {
            Eigen::Matrix<T,3,3> phi_ii;
            Eigen::Matrix<T,3,3> q_ii;
            Eigen::Matrix<T,4,4> phi_iii;
            Eigen::Matrix<T,4,4> q_iii;
            ChainII::transition(tau, step, phi_ii);
            ChainII::process_covariance(tau, step, sigma2, q_ii);
            ChainIII::transition(tau, step, phi_iii);
            ChainIII::process_covariance(tau, step, sigma2, q_iii);

            Eigen::Matrix<T,3,3> phi_iii_marginal;
            Eigen::Matrix<T,3,3> q_iii_marginal;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    phi_iii_marginal(i,j) = phi_iii(keep[i], keep[j]);
                    q_iii_marginal(i,j) = q_iii(keep[i], keep[j]);
                }
            }

            if (!near_matrix(phi_ii, phi_iii_marginal, abs_tol, rel_tol, "transition") ||
                !near_matrix(q_ii, q_iii_marginal, abs_tol, rel_tol, "covariance")) {
                std::cerr << "tau=" << tau << " step=" << step << '\n';
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main() {
    if (!check_chain_parity<float>() ||
        !check_chain_parity<double>()) {
        return 1;
    }
    std::cout << "OU core parity tests passed\n";
    return 0;
}
