#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>

#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

namespace {

template <typename MatrixT>
bool psd(const MatrixT& P, double tol = 1e-9) {
    Eigen::MatrixXd S = 0.5 * (P.template cast<double>() + P.transpose().template cast<double>());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    return es.info() == Eigen::Success && es.eigenvalues().minCoeff() >= -tol;
}

template <typename Filter>
bool exercise(Filter& f, int off_aw, const char* name) {
    if (f.legacy_aw_covariance_replacement()) {
        std::cerr << name << ": legacy replacement must default to false\n";
        return false;
    }

    const auto before = f.covariance_full();
    Eigen::Vector3d std_hi(5.0, 6.0, 7.0);
    f.set_aw_stationary_std(std_hi);
    f.synchronize_aw_covariance_to_stationary();

    const auto queued = f.covariance_full();
    if ((queued - before).norm() > 1e-12) {
        std::cerr << name << ": default synchronization rewrote posterior covariance\n";
        return false;
    }

    f.time_update(Eigen::Vector3d::Zero(), 0.005);
    const auto after = f.covariance_full();
    Eigen::Matrix3d target = std_hi.array().square().matrix().asDiagonal();
    Eigen::Matrix3d margin = after.template block<3,3>(off_aw, off_aw) - target;
    if (!psd(margin, 1e-8) || !psd(after, 1e-8)) {
        std::cerr << name << ": prediction-time PSD inflation failed floor/PSD test\n";
        return false;
    }

    f.set_legacy_aw_covariance_replacement(true);
    if (!f.legacy_aw_covariance_replacement()) {
        std::cerr << name << ": legacy flag did not latch\n";
        return false;
    }
    Eigen::Vector3d std_legacy(1.1, 1.2, 1.3);
    f.set_aw_stationary_std(std_legacy);
    f.synchronize_aw_covariance_to_stationary();
    const auto legacy = f.covariance_full();
    Eigen::Matrix3d legacy_target = std_legacy.array().square().matrix().asDiagonal();
    if ((legacy.template block<3,3>(off_aw, off_aw) - legacy_target).norm() > 1e-10) {
        std::cerr << name << ": legacy mode did not reproduce block replacement\n";
        return false;
    }
    return true;
}

} // namespace

int main() {
    const Eigen::Vector3d sigma_a = Eigen::Vector3d::Constant(0.02);
    const Eigen::Vector3d sigma_g = Eigen::Vector3d::Constant(0.001);
    const Eigen::Vector3d sigma_m = Eigen::Vector3d::Constant(0.1);

    Kalman3D_Wave_OU_II<double> ou2(sigma_a, sigma_g, sigma_m);
    Kalman3D_Wave_OU_III<double> ou3(sigma_a, sigma_g, sigma_m);

    if (!exercise(ou2, 12, "OU-II")) return 1;
    if (!exercise(ou3, 15, "OU-III")) return 1;

    std::cout << "aw covariance policy: PASS\n";
    return 0;
}
