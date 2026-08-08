#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

namespace {

template <typename Filter>
bool exercise(Filter& f, const char* name) {
    using Vec3 = Eigen::Vector3d;

    const double tau = 1000.0;
    const double dt = 0.01;
    const int steps = 200;
    const double elapsed = dt * steps;
    const Vec3 b0(0.12, -0.07, 0.03);
    const Vec3 drive = Vec3::Constant(5e-4);
    const Vec3 kT(0.002, -0.001, 0.0005);

    f.set_acc_bias_time_constant(tau);
    f.set_Q_bacc_rw(drive);
    f.set_accel_bias_temp_coeff(kT);
    f.set_initial_acc_bias(b0);

    const auto P0 = f.covariance_full();
    const int off_ba = static_cast<int>(P0.rows()) - 3;
    const Eigen::Matrix3d Pbb0 = P0.template block<3,3>(off_ba, off_ba);

    for (int i = 0; i < steps; ++i)
        f.time_update(Vec3::Zero(), dt);

    const double phi = std::exp(-elapsed / tau);
    const Vec3 expected_b = phi * b0;
    const Vec3 got_b = f.get_acc_bias();
    if ((got_b - expected_b).norm() > 1e-10) {
        std::cerr << name << ": residual bias did not follow exact OU mean transition\n";
        return false;
    }

    const double phi2 = std::exp(-2.0 * elapsed / tau);
    const Eigen::Matrix3d Qc = drive.array().square().matrix().asDiagonal();
    const Eigen::Matrix3d expected_P = phi2 * Pbb0 + (0.5 * tau * (1.0 - phi2)) * Qc;
    const auto P = f.covariance_full();
    const Eigen::Matrix3d got_P = P.template block<3,3>(off_ba, off_ba);
    if ((got_P - expected_P).norm() > 1e-10) {
        std::cerr << name << ": residual bias covariance did not follow exact OU covariance\n";
        return false;
    }

    const double temp = 45.0;
    const Vec3 expected_total = expected_b + kT * (temp - 35.0);
    if ((f.get_acc_bias_at_temperature(temp) - expected_total).norm() > 1e-10) {
        std::cerr << name << ": temperature mean plus OU residual is inconsistent\n";
        return false;
    }

    // Freezing bias updates retains the historical warmup semantics: no OU
    // mean reversion or diffusion while the bias state is explicitly frozen.
    f.set_acc_bias_updates_enabled(false);
    const Vec3 frozen_b = f.get_acc_bias();
    const Eigen::Matrix3d frozen_P = f.covariance_full().template block<3,3>(off_ba, off_ba);
    for (int i = 0; i < 20; ++i)
        f.time_update(Vec3::Zero(), dt);
    if ((f.get_acc_bias() - frozen_b).norm() > 1e-12 ||
        (f.covariance_full().template block<3,3>(off_ba, off_ba) - frozen_P).norm() > 1e-12) {
        std::cerr << name << ": frozen bias state changed during warmup-style propagation\n";
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

    if (!exercise(ou2, "OU-II")) return 1;
    if (!exercise(ou3, "OU-III")) return 1;

    std::cout << "accelerometer bias OU propagation: PASS\n";
    return 0;
}
