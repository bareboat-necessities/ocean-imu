#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

namespace {

using Filter = Kalman3D_Wave_OU_III<double>;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

constexpr int OFF_V = 6;
constexpr int OFF_P = 9;
constexpr int OFF_S = 12;
constexpr int OFF_AW = 15;

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

double quat_distance_rad(const Eigen::Quaterniond& a, const Eigen::Quaterniond& b) {
    Eigen::Quaterniond dq = a.conjugate() * b;
    dq.normalize();
    const double w = std::clamp(std::abs(dq.w()), 0.0, 1.0);
    return 2.0 * std::acos(w);
}

Filter make_filter() {
    const Vec3 sigma_a = Vec3::Constant(0.02);
    const Vec3 sigma_g = Vec3::Constant(0.001);
    const Vec3 sigma_m = Vec3::Constant(0.1);
    Filter f(sigma_a, sigma_g, sigma_m);
    f.set_pseudo_update_period_s(1000.0);  // keep setup free of automatic S updates
    f.set_RS_noise(Vec3::Constant(0.10));
    f.set_aw_time_constant(2.0);
    f.set_aw_stationary_std(Vec3::Constant(1.0));
    f.initialize_from_truth(Vec3(2.0, -1.0, 0.5), Vec3::Zero(),
                            Eigen::Quaterniond::Identity(), Vec3::Zero());
    return f;
}

bool test_s_pseudo_freezes_attitude_but_keeps_linear_correction() {
    Filter f = make_filter();
    const Vec3 acc_rest(0.0, 0.0, -9.80665);
    const Vec3 gyr = Vec3::Zero();
    const double dt = 0.005;

    // Build the covariance path P_{theta,S} naturally: OU process noise creates
    // S-a_w covariance, then the accelerometer update couples a_w to attitude.
    // A nonzero initial position also makes the S=0 innovation appreciable.
    for (int k = 0; k < 1200; ++k) {
        f.time_update(gyr, dt);
        if ((k % 4) == 0) f.measurement_update_acc_only(acc_rest);
    }

    const auto P0 = f.covariance_full();
    const Mat3 PthetaS = P0.template block<3,3>(0, OFF_S);
    const Mat3 PSS = P0.template block<3,3>(OFF_S, OFF_S);
    const Mat3 R = Mat3::Identity() * 0.01; // sigma_S = 0.10
    const Mat3 innovation_cov = PSS + R;
    const Mat3 Ktheta_unrestricted = PthetaS * innovation_cov.inverse();
    const Vec3 r = -f.get_integral_displacement();
    const Vec3 unrestricted_dtheta = Ktheta_unrestricted * r;

    bool ok = true;
    ok &= check(PthetaS.norm() > 1e-12,
                "setup must create nonzero attitude-S cross covariance");
    ok &= check(unrestricted_dtheta.norm() > 1e-12,
                "ordinary full gain would have produced an attitude correction");
    ok &= check(r.norm() > 1e-4,
                "setup must create a nonzero S pseudo-measurement innovation");

    const Eigen::Quaterniond q_before = f.quaternion();
    const Vec3 v_before = f.get_velocity();
    const Vec3 p_before = f.get_position();
    const Vec3 S_before = f.get_integral_displacement();
    const Vec3 aw_before = f.get_world_accel();

    f.applyIntegralZeroPseudoMeas();

    const Eigen::Quaterniond q_after = f.quaternion();
    const double dq = quat_distance_rad(q_before, q_after);
    ok &= check(dq <= 1e-13,
                "Schmidt S pseudo update must not inject attitude");

    const double linear_change =
        (f.get_velocity() - v_before).norm() +
        (f.get_position() - p_before).norm() +
        (f.get_integral_displacement() - S_before).norm() +
        (f.get_world_accel() - aw_before).norm();
    ok &= check(linear_change > 1e-10,
                "S pseudo update must continue correcting permitted linear states");

    const auto P1 = f.covariance_full();
    ok &= check((P1 - P1.transpose()).norm() <= 1e-9,
                "restricted-gain Joseph update must preserve covariance symmetry");
    Eigen::SelfAdjointEigenSolver<decltype(P1)> es(P1);
    ok &= check(es.info() == Eigen::Success,
                "covariance eigensolver must succeed");
    if (es.info() == Eigen::Success) {
        ok &= check(es.eigenvalues().minCoeff() >= -1e-9,
                    "restricted-gain Joseph update must preserve PSD within tolerance");
    }

    // Repeated S updates must remain numerically healthy.
    for (int k = 0; k < 100; ++k) f.applyIntegralZeroPseudoMeas();
    const auto Pr = f.covariance_full();
    Eigen::SelfAdjointEigenSolver<decltype(Pr)> esr(Pr);
    ok &= check((Pr - Pr.transpose()).norm() <= 1e-9,
                "repeated Schmidt S updates must preserve covariance symmetry");
    ok &= check(esr.info() == Eigen::Success && esr.eigenvalues().minCoeff() >= -1e-9,
                "repeated Schmidt S updates must preserve PSD within tolerance");

    // The freeze is measurement-specific: an accepted accelerometer update must
    // still be able to move attitude.
    const Eigen::Quaterniond q_acc_before = f.quaternion();
    f.measurement_update_acc_only(Vec3(0.05, 0.0, -9.80665));
    ok &= check(f.lastAccDiag().accepted,
                "mild tilted accelerometer regression update should be accepted");
    if (f.lastAccDiag().accepted) {
        ok &= check(quat_distance_rad(q_acc_before, f.quaternion()) > 1e-12,
                    "accelerometer update must retain ordinary attitude correction");
    }

    std::cout << "Schmidt S diagnostics: PthetaS_norm=" << PthetaS.norm()
              << " unrestricted_dtheta_norm=" << unrestricted_dtheta.norm()
              << " S_innovation_norm=" << r.norm()
              << " schmidt_quat_delta=" << dq
              << " linear_change=" << linear_change << '\n';

    if (!ok) {
        std::cerr << "PthetaS_norm=" << PthetaS.norm()
                  << " unrestricted_dtheta_norm=" << unrestricted_dtheta.norm()
                  << " S_innovation_norm=" << r.norm()
                  << " schmidt_quat_delta=" << dq << '\n';
    }
    return ok;
}

} // namespace

int main() {
    if (!test_s_pseudo_freezes_attitude_but_keeps_linear_correction()) return 1;
    std::cout << "schmidt_s_pseudo-test: PASS\n";
    return 0;
}
