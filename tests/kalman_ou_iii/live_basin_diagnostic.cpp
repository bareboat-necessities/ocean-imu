#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

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
constexpr int HORIZON_STEPS = 30 * 200;

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
    // Singular values only.  Thin U/V are invalid for fixed-size Eigen
    // matrices and are unnecessary for an induced 2-norm.
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
        A = correction_transition(f.K_scratch_, C) * A;
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

// Induced norm from the covariance/Riccati metric at the start of a lifted
// interval to the metric at its end.  If P0=L0 L0^T and P1=L1 L1^T, then
// ||Psi||_{P0^{-1}->P1^{-1}} = ||L1^{-1} Psi L0||_2.
double covariance_metric_norm(const Mat21& Psi,
                              const Mat21& P0,
                              const Mat21& P1) {
    Eigen::LLT<Mat21> llt0(P0);
    Eigen::LLT<Mat21> llt1(P1);
    if (llt0.info() != Eigen::Success || llt1.info() != Eigen::Success) {
        return std::numeric_limits<double>::infinity();
    }
    const Mat21 L0 = llt0.matrixL();
    const Mat21 B = llt1.matrixL().solve(Psi * L0);
    return norm2(B);
}

std::pair<double,double> scaled_cov_eigen_bounds(const Mat21& P,
                                                  const Mat21& D) {
    const Mat21 Pz = D * P * D;
    Eigen::SelfAdjointEigenSolver<Mat21> es(Pz, Eigen::EigenvaluesOnly);
    if (es.info() != Eigen::Success) {
        return {0.0, std::numeric_limits<double>::infinity()};
    }
    return {es.eigenvalues()(0), es.eigenvalues()(20)};
}

struct Report {
    double chi_euclid = std::numeric_limits<double>::quiet_NaN();
    double prefix_euclid = 0.0;
    double chi_metric = std::numeric_limits<double>::quiet_NaN();
    double rho_metric_sample = 0.0;
    double rho_metric_second = 0.0;
    double pz_min = std::numeric_limits<double>::infinity();
    double pz_max = 0.0;
    double kappa_euclid = std::numeric_limits<double>::infinity();
    double xi_xi = 0.0;
    double xi_ell = 0.0;
    double xi_all = 0.0;
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
    // proof; it makes the margin diagnostic representative of established
    // Live operation rather than constructor transients.
    constexpr int warm_steps = 120 * 200;
    for (int k = 0; k < warm_steps; ++k) {
        (void)step(f, k, mag_world);
    }

    const Mat21 S = state_scale();
    Mat21 D = Mat21::Zero();
    for (int i = 0; i < 21; ++i) D(i,i) = 1.0 / S(i,i);

    const Mat21 P0 = f.Pext;
    Mat21 Psi_phys = Mat21::Identity();
    Mat21 Psi_scaled = Mat21::Identity();
    double prefix_euclid = 1.0;
    int pseudo_count = 0;

    auto [pz_min, pz_max] = scaled_cov_eigen_bounds(P0, D);

    for (int k = 0; k < HORIZON_STEPS; ++k) {
        const StepResult sr = step(f, warm_steps + k, mag_world);
        if (sr.pseudo) ++pseudo_count;

        Psi_phys = sr.A * Psi_phys;
        const Mat21 Abar = D * sr.A * S;
        Psi_scaled = Abar * Psi_scaled;
        prefix_euclid = std::max(prefix_euclid, norm2(Psi_scaled));

        const auto [lo, hi] = scaled_cov_eigen_bounds(f.Pext, D);
        pz_min = std::min(pz_min, lo);
        pz_max = std::max(pz_max, hi);
    }

    Report r;
    r.chi_euclid = norm2(Psi_scaled);
    r.prefix_euclid = prefix_euclid;
    r.chi_metric = covariance_metric_norm(Psi_phys, P0, f.Pext);
    r.pseudo_count = pseudo_count;
    r.pz_min = pz_min;
    r.pz_max = pz_max;
    if (pz_min > 0.0 && std::isfinite(pz_max)) {
        r.kappa_euclid = std::sqrt(pz_max / pz_min);
    }
    if (r.chi_metric > 0.0 && r.chi_metric < 1.0) {
        r.rho_metric_sample =
            std::pow(r.chi_metric, 1.0 / double(HORIZON_STEPS));
        r.rho_metric_second = std::pow(r.chi_metric, 1.0 / 30.0);
    }

    std::array<int, 15> xi{{0,1,2,3,4,5,12,13,14,15,16,17,18,19,20}};
    std::array<int, 6> ell{{6,7,8,9,10,11}};
    Eigen::Matrix<double,15,15> Pxx;
    Eigen::Matrix<double,15,6> Pxl;
    Eigen::Matrix<double,15,21> Px;
    for (int i = 0; i < 15; ++i) {
        for (int j = 0; j < 15; ++j) Pxx(i,j) = Psi_scaled(xi[i], xi[j]);
        for (int j = 0; j < 6; ++j) Pxl(i,j) = Psi_scaled(xi[i], ell[j]);
        for (int j = 0; j < 21; ++j) Px(i,j) = Psi_scaled(xi[i], j);
    }
    r.xi_xi = norm2_fixed(Pxx);
    r.xi_ell = norm2_fixed(Pxl);
    r.xi_all = norm2_fixed(Px);
    return r;
}

}  // namespace

int main() {
    std::cout << std::setprecision(12);
    std::cout << "PHASE2_SCALES s_theta=0.087 s_bg=0.001 s_v=1 s_p=20 s_S=50 s_aw=6 s_ba=0.5\n";
    std::cout << "name,tau,sigma_aw,rS,chiE30,prefixE30,chiV30,rhoV_sample,rhoV_second,pz_min,pz_max,kappaE,xi_xiE30,xi_ellE30,xi_allE30,pseudo_count\n";

    double worst_chi_metric = 0.0;
    double worst_chi_euclid = 0.0;
    double worst_prefix_euclid = 0.0;
    double worst_rho_metric_second = 0.0;
    double global_pz_min = std::numeric_limits<double>::infinity();
    double global_pz_max = 0.0;
    bool all_metric_contract = true;

    for (const auto& op : kReferencePoints) {
        const Report r = evaluate(op);
        std::cout << op.name << ',' << op.tau << ',' << op.sigma_aw << ',' << op.rS
                  << ',' << r.chi_euclid << ',' << r.prefix_euclid
                  << ',' << r.chi_metric << ',' << r.rho_metric_sample
                  << ',' << r.rho_metric_second << ',' << r.pz_min << ',' << r.pz_max
                  << ',' << r.kappa_euclid << ',' << r.xi_xi << ',' << r.xi_ell
                  << ',' << r.xi_all << ',' << r.pseudo_count << '\n';

        all_metric_contract = all_metric_contract
            && std::isfinite(r.chi_metric) && r.chi_metric < 1.0
            && std::isfinite(r.pz_min) && r.pz_min > 0.0
            && std::isfinite(r.pz_max);
        worst_chi_metric = std::max(worst_chi_metric, r.chi_metric);
        worst_chi_euclid = std::max(worst_chi_euclid, r.chi_euclid);
        worst_prefix_euclid = std::max(worst_prefix_euclid, r.prefix_euclid);
        worst_rho_metric_second =
            std::max(worst_rho_metric_second, r.rho_metric_second);
        global_pz_min = std::min(global_pz_min, r.pz_min);
        global_pz_max = std::max(global_pz_max, r.pz_max);
    }

    const double global_kappa =
        (global_pz_min > 0.0 && std::isfinite(global_pz_max))
        ? std::sqrt(global_pz_max / global_pz_min)
        : std::numeric_limits<double>::infinity();
    const double rho_sample =
        (worst_chi_metric > 0.0 && worst_chi_metric < 1.0)
        ? std::pow(worst_chi_metric, 1.0 / double(HORIZON_STEPS)) : 0.0;
    const double M_euclid =
        (rho_sample > 0.0 && std::isfinite(global_kappa))
        ? global_kappa * std::pow(rho_sample, -(HORIZON_STEPS - 1))
        : std::numeric_limits<double>::infinity();

    std::cout << "PHASE2_SUMMARY worst_chiV30=" << worst_chi_metric
              << " worst_chiE30=" << worst_chi_euclid
              << " worst_prefixE30=" << worst_prefix_euclid
              << " worst_rhoV_second=" << worst_rho_metric_second
              << " global_pz_min=" << global_pz_min
              << " global_pz_max=" << global_pz_max
              << " global_kappaE=" << global_kappa
              << " implied_M_E=" << M_euclid
              << " all_reference_metric_contract=" << (all_metric_contract ? 1 : 0)
              << '\n';

    // The fixed Euclidean norm is deliberately reported but not required to
    // contract over 30 s: the 5000 s residual-bias OU mode produces a long
    // non-normal transient.  The constructive certificate uses the Riccati
    // covariance as its time-varying Lyapunov metric and converts back to the
    // fixed physical coordinates through the reported covariance eigenvalue
    // bounds.  Reference replay remains evidence, not a uniform proof.
    if (!all_metric_contract) {
        std::cerr << "FAIL: at least one deployed reference operating point did not contract in the Riccati metric over 30 s\n";
        return 1;
    }
    return 0;
}
