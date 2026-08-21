#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#define EIGEN_NON_ARDUINO
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private

namespace {

using Core = Kalman3D_Wave_OU_III<double, true, true>;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Mat21 = Eigen::Matrix<double, 21, 21>;
using Mat3x21 = Eigen::Matrix<double, 3, 21>;

constexpr int OFF_TH = 0;
constexpr int OFF_BG = 3;
constexpr int OFF_V  = 6;
constexpr int OFF_P  = 9;
constexpr int OFF_S  = 12;
constexpr int OFF_AW = 15;
constexpr int OFF_BA = 18;
constexpr double DT = 1.0 / 200.0;

struct OperatingPoint {
    const char* name;
    double tau;
    double sigma_aw;
    double rS;
};

// Current C_J=0.0538 reference-sea operating points.  The committed sweep used
// c_sigma=0.6; production uses c_sigma=0.9, while the SpectralMSE r_S law
// divides c_sigma back out.  Therefore tau and r_S are unchanged and sigma_aw
// is multiplied by 0.9/0.6=1.5 here.
const std::array<OperatingPoint, 8> kReferencePoints{{
    {"J0.27", 1.28051, 0.239768 * 1.5, 0.164650},
    {"J1.50", 2.18038, 0.486753 * 1.5, 1.532980},
    {"J4.00", 3.58993, 0.749303 * 1.5, 8.849600},
    {"J8.50", 4.22171, 0.947487 * 1.5, 17.559300},
    {"P0.27", 1.17236, 0.257598 * 1.5, 0.152396},
    {"P1.50", 2.04464, 0.526747 * 1.5, 1.370490},
    {"P4.00", 3.29004, 0.750996 * 1.5, 6.484430},
    {"P8.50", 4.10091, 0.997636 * 1.5, 15.641200},
}};

Mat3 skew(const Vec3& v) {
    Mat3 M;
    M << 0.0, -v.z(), v.y(),
         v.z(), 0.0, -v.x(),
        -v.y(), v.x(), 0.0;
    return M;
}

double norm2(const Mat21& A) {
    Eigen::JacobiSVD<Mat21> svd(A);
    return svd.singularValues()(0);
}

template <int R, int C>
double norm2_fixed(const Eigen::Matrix<double, R, C>& A) {
    Eigen::JacobiSVD<Eigen::Matrix<double, R, C>> svd(A);
    return svd.singularValues()(0);
}

Mat21 state_scale() {
    // Fixed proof coordinates.  These are design scales, not fitted weights:
    // magnetically gauged handoff yaw sigma, constructor gyro-bias sigma,
    // constructor v/p/S sigmas, wrapper a_w safety ceiling, and the hard b_a
    // projection radius respectively.
    Mat21 S = Mat21::Zero();
    S.diagonal().segment<3>(OFF_TH).setConstant(0.087);  // rad
    S.diagonal().segment<3>(OFF_BG).setConstant(0.001);  // rad/s
    S.diagonal().segment<3>(OFF_V ).setConstant(1.0);    // m/s
    S.diagonal().segment<3>(OFF_P ).setConstant(20.0);   // m
    S.diagonal().segment<3>(OFF_S ).setConstant(50.0);   // m s
    S.diagonal().segment<3>(OFF_AW).setConstant(6.0);    // m/s^2
    S.diagonal().segment<3>(OFF_BA).setConstant(0.5);    // m/s^2
    return S;
}

Mat21 correction_transition(const Eigen::Matrix<double, 21, 3>& K,
                            const Mat3x21& C) {
    return Mat21::Identity() - K * C;
}

Mat21 prediction_transition(const Core& f) {
    Mat21 F = Mat21::Identity();
    F.block<6,6>(0,0) = f.F_AA_scratch_;
    F.block<12,12>(OFF_V, OFF_V) = f.F_LL_scratch_;
    const double phi_b = std::exp(-DT / f.tau_bacc_);
    F.block<3,3>(OFF_BA, OFF_BA) = Mat3::Identity() * phi_b;
    return F;
}

Mat3x21 accel_jacobian(const Core& f) {
    Mat3x21 C = Mat3x21::Zero();
    const Mat3 R = f.R_wb();
    const Vec3 g(0.0, 0.0, 9.80665);
    const Vec3 aw = f.xext.segment<3>(OFF_AW);
    const Vec3 f_cog = R * (aw - g);
    C.block<3,3>(0, OFF_TH) = -skew(f_cog);
    C.block<3,3>(0, OFF_AW) = R;
    C.block<3,3>(0, OFF_BA) = Mat3::Identity();
    return C;
}

Mat3x21 mag_jacobian(const Core& f) {
    Mat3x21 C = Mat3x21::Zero();
    const Vec3 mhat = f.R_wb() * f.v2ref;
    C.block<3,3>(0, OFF_TH) = -skew(mhat);
    return C;
}

Mat3x21 S_jacobian() {
    Mat3x21 C = Mat3x21::Zero();
    C.block<3,3>(0, OFF_S) = Mat3::Identity();
    return C;
}

struct StepResult {
    Mat21 A = Mat21::Identity();
    bool pseudo = false;
};

StepResult step(Core& f, int k, const Vec3& mag_world) {
    StepResult out;
    const Vec3 gyro = Vec3::Zero();
    const Vec3 acc(0.0, 0.0, -9.80665);

    const double elapsed_before = f.pseudo_update_elapsed_s_;
    f.time_update(gyro, DT);
    Mat21 A = prediction_transition(f);

    const double elapsed_after = f.pseudo_update_elapsed_s_;
    const bool pseudo = elapsed_after < elapsed_before + DT - 1e-10;
    if (pseudo) {
        const Mat3x21 C = S_jacobian();
        const Mat21 J = correction_transition(f.K_scratch_, C);
        A = J * A;
    }

    const Mat3x21 Ca = accel_jacobian(f);
    f.measurement_update_acc_only(acc, 35.0);
    A = correction_transition(f.K_scratch_, Ca) * A;

    // Reference simulations provide a 100 Hz magnetometer against 200 Hz IMU.
    if ((k & 1) == 0) {
        const Mat3x21 Cm = mag_jacobian(f);
        f.measurement_update_mag_only(mag_world);
        A = correction_transition(f.K_scratch_, Cm) * A;
    }

    out.A = A;
    out.pseudo = pseudo;
    return out;
}

struct Report {
    double chi = std::numeric_limits<double>::quiet_NaN();
    double prefix = 0.0;
    double xi_xi = 0.0;
    double xi_ell = 0.0;
    double xi_all = 0.0;
    double rho_sample = 0.0;
    double rho_second = 0.0;
    double M = 0.0;
    int pseudo_count = 0;
};

Report evaluate(const OperatingPoint& op) {
    const Vec3 sigma_a = Vec3::Constant(0.0294);
    const Vec3 sigma_g = Vec3::Constant(0.000157);
    const Vec3 sigma_m = Vec3::Constant(0.36);
    Core f(sigma_a, sigma_g, sigma_m);

    const Vec3 mag_world(20.0, 5.0, 44.0);
    f.set_mag_world_ref(mag_world);
    f.initialize_from_attitude(Eigen::Quaterniond::Identity(), 0.035, 0.087);
    f.set_linear_block_enabled(true);
    f.set_acc_bias_updates_enabled(true);
    f.set_aw_time_constant(op.tau);
    f.set_aw_stationary_std(Vec3::Constant(op.sigma_aw));
    f.set_RS_noise(Vec3::Constant(op.rS));
    const double period = std::clamp((0.015 / 1.1) * op.tau, 0.005, 0.25);
    f.set_pseudo_update_period_s(period);
    f.reset_aw_covariance_to_stationary();

    // Covariance warm-up at the exact equilibrium.  This is not part of the
    // proof; it only makes the margin diagnostic representative of established
    // Live operation rather than constructor transients.
    constexpr int warm_steps = 120 * 200;
    for (int k = 0; k < warm_steps; ++k) {
        (void)step(f, k, mag_world);
    }

    const Mat21 S = state_scale();
    Mat21 D = Mat21::Zero();
    for (int i = 0; i < 21; ++i) D(i,i) = 1.0 / S(i,i);

    // Thirty seconds is long enough to include many S and magnetic updates and
    // is still short compared with the 5000 s residual-bias time constant.
    constexpr int H = 30 * 200;
    Mat21 Psi = Mat21::Identity();
    double prefix = 1.0;
    int pseudo_count = 0;
    for (int k = 0; k < H; ++k) {
        const StepResult sr = step(f, warm_steps + k, mag_world);
        if (sr.pseudo) ++pseudo_count;
        const Mat21 Abar = D * sr.A * S;
        Psi = Abar * Psi;
        // Exact prefix maximum for this diagnostic window.  This is reported
        // as evidence only; a future runtime certificate must check its own
        // matrix inequalities rather than trust this replay value.
        prefix = std::max(prefix, norm2(Psi));
    }

    Report r;
    r.chi = norm2(Psi);
    r.prefix = prefix;
    r.pseudo_count = pseudo_count;

    std::array<int, 15> xi{{0,1,2,3,4,5,12,13,14,15,16,17,18,19,20}};
    std::array<int, 6> ell{{6,7,8,9,10,11}};
    Eigen::Matrix<double,15,15> Pxx;
    Eigen::Matrix<double,15,6> Pxl;
    Eigen::Matrix<double,15,21> Px;
    for (int i = 0; i < 15; ++i) {
        for (int j = 0; j < 15; ++j) Pxx(i,j) = Psi(xi[i], xi[j]);
        for (int j = 0; j < 6; ++j) Pxl(i,j) = Psi(xi[i], ell[j]);
        for (int j = 0; j < 21; ++j) Px(i,j) = Psi(xi[i], j);
    }
    r.xi_xi = norm2_fixed(Pxx);
    r.xi_ell = norm2_fixed(Pxl);
    r.xi_all = norm2_fixed(Px);

    if (r.chi > 0.0 && r.chi < 1.0) {
        r.rho_sample = std::pow(r.chi, 1.0 / double(H));
        r.rho_second = std::pow(r.chi, 1.0 / 30.0);
        r.M = r.prefix * std::pow(r.rho_sample, -(H - 1));
    }
    return r;
}

}  // namespace

int main() {
    std::cout << std::setprecision(10);
    std::cout << "PHASE2_SCALES s_theta=0.087 s_bg=0.001 s_v=1 s_p=20 s_S=50 s_aw=6 s_ba=0.5\n";
    std::cout << "name,tau,sigma_aw,rS,chi30,prefix30,xi_xi30,xi_ell30,xi_all30,rho_sample,rho_second,M,pseudo_count\n";

    double worst_chi = 0.0;
    double worst_prefix = 0.0;
    double worst_M = 0.0;
    double worst_rho_second = 0.0;
    bool all_contract = true;

    for (const auto& op : kReferencePoints) {
        const Report r = evaluate(op);
        std::cout << op.name << ',' << op.tau << ',' << op.sigma_aw << ',' << op.rS
                  << ',' << r.chi << ',' << r.prefix
                  << ',' << r.xi_xi << ',' << r.xi_ell << ',' << r.xi_all
                  << ',' << r.rho_sample << ',' << r.rho_second << ',' << r.M
                  << ',' << r.pseudo_count << '\n';
        all_contract = all_contract && std::isfinite(r.chi) && r.chi < 1.0;
        worst_chi = std::max(worst_chi, r.chi);
        worst_prefix = std::max(worst_prefix, r.prefix);
        worst_M = std::max(worst_M, r.M);
        worst_rho_second = std::max(worst_rho_second, r.rho_second);
    }

    std::cout << "PHASE2_SUMMARY worst_chi30=" << worst_chi
              << " worst_prefix30=" << worst_prefix
              << " worst_M=" << worst_M
              << " worst_rho_second=" << worst_rho_second
              << " all_reference_contract=" << (all_contract ? 1 : 0) << '\n';

    // This diagnostic is not a uniform proof over the admissible schedule, but
    // every committed reference point should at least exhibit finite-horizon
    // contraction.  The article is explicit about this distinction.
    if (!all_contract) {
        std::cerr << "FAIL: at least one deployed reference operating point did not contract over 30 s\n";
        return 1;
    }
    return 0;
}
