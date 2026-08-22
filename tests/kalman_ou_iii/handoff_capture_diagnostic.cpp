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
using Mat15 = Eigen::Matrix<double, 15, 15>;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using Mat15x6 = Eigen::Matrix<double, 15, 6>;

constexpr int OFF_TH = 0;
constexpr int OFF_BG = 3;
constexpr int OFF_V  = 6;
constexpr int OFF_P  = 9;
constexpr int OFF_S  = 12;
constexpr int OFF_AW = 15;
constexpr int OFF_BA = 18;
constexpr double DT = 1.0 / 200.0;
constexpr int MAX_SECONDS = 180;

const std::array<int, 15> kFast{{
    0,1,2,          // attitude
    6,7,8,          // velocity
    9,10,11,        // position
    12,13,14,       // integral displacement
    15,16,17        // OU acceleration
}};
const std::array<int, 6> kSlow{{3,4,5,18,19,20}}; // gyro and accel bias

struct OperatingPoint {
    std::string name;
    double tau;
    double sigma_aw;
    double rS;
};

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

template <typename M>
double norm2(const M& A) {
    Eigen::JacobiSVD<M> svd(A);
    return svd.singularValues()(0);
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

Mat21 step(Core& f, int k, const Vec3& mag_world) {
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

    if ((k & 1) == 0) {
        const Mat3x21 Cm = mag_jacobian(f);
        f.measurement_update_mag_only(mag_world);
        A = correction_transition(f.K_scratch_, Cm) * A;
    }
    return A;
}

Mat15 fast_cov(const Mat21& P) {
    Mat15 out;
    for (int i = 0; i < 15; ++i)
        for (int j = 0; j < 15; ++j)
            out(i,j) = P(kFast[i], kFast[j]);
    return out;
}

Mat6 slow_cov(const Mat21& P) {
    Mat6 out;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            out(i,j) = P(kSlow[i], kSlow[j]);
    return out;
}

Mat15 fast_fast(const Mat21& Psi) {
    Mat15 out;
    for (int i = 0; i < 15; ++i)
        for (int j = 0; j < 15; ++j)
            out(i,j) = Psi(kFast[i], kFast[j]);
    return out;
}

Mat15x6 fast_slow(const Mat21& Psi) {
    Mat15x6 out;
    for (int i = 0; i < 15; ++i)
        for (int j = 0; j < 6; ++j)
            out(i,j) = Psi(kFast[i], kSlow[j]);
    return out;
}

double metric_norm_full(const Mat21& Psi, const Mat21& P0, const Mat21& P1) {
    Eigen::LLT<Mat21> l0(P0), l1(P1);
    if (l0.info() != Eigen::Success || l1.info() != Eigen::Success)
        return std::numeric_limits<double>::infinity();
    const Mat21 L0 = l0.matrixL();
    const Mat21 B = l1.matrixL().solve(Psi * L0);
    return norm2(B);
}

double metric_norm_fast(const Mat21& Psi, const Mat21& P0, const Mat21& P1) {
    const Mat15 P0f = fast_cov(P0);
    const Mat15 P1f = fast_cov(P1);
    const Mat15 Aff = fast_fast(Psi);
    Eigen::LLT<Mat15> l0(P0f), l1(P1f);
    if (l0.info() != Eigen::Success || l1.info() != Eigen::Success)
        return std::numeric_limits<double>::infinity();
    const Mat15 L0 = l0.matrixL();
    const Mat15 B = l1.matrixL().solve(Aff * L0);
    return norm2(B);
}

double metric_gain_slow_to_fast(const Mat21& Psi,
                                const Mat21& P0,
                                const Mat21& P1) {
    const Mat6 P0s = slow_cov(P0);
    const Mat15 P1f = fast_cov(P1);
    const Mat15x6 Afs = fast_slow(Psi);
    Eigen::LLT<Mat6> l0(P0s);
    Eigen::LLT<Mat15> l1(P1f);
    if (l0.info() != Eigen::Success || l1.info() != Eigen::Success)
        return std::numeric_limits<double>::infinity();
    const Mat6 L0 = l0.matrixL();
    const Mat15x6 B = l1.matrixL().solve(Afs * L0);
    return norm2(B);
}

double spectral_mse_rS(double tau, double sigma_aw) {
    constexpr double CJ = 0.0538;
    constexpr double c_sigma = 0.9;
    constexpr double r_a = 0.0148 * 0.0148 * DT;
    constexpr double cT = 0.015 / 1.1;
    const double TS = std::clamp(cT * tau, 0.005, 0.25);
    const double sigma_aB = std::max(sigma_aw / c_sigma, 1e-6);
    const double qpow = std::pow(2.0 * r_a, 1.0 / 14.0);
    const double u = sigma_aB * tau * tau * tau * tau;
    const double raw = CJ * qpow * std::pow(u, 6.0 / 7.0) / std::sqrt(TS);
    return std::clamp(raw, 0.15, 400.0);
}

std::vector<OperatingPoint> source_points() {
    const std::array<double, 6> tau{{1.0/3.0, 1.0, 2.0, 4.0, 8.0, 12.0}};
    const std::array<double, 4> sig{{0.05, 0.1, 1.0, 6.0}};
    std::vector<OperatingPoint> out;
    for (double t : tau) {
        for (double s : sig) {
            OperatingPoint p;
            p.name = "SRC_t" + std::to_string(t) + "_s" + std::to_string(s);
            p.tau = t;
            p.sigma_aw = s;
            p.rS = spectral_mse_rS(t, s);
            out.push_back(p);
        }
    }
    return out;
}

struct CaptureReport {
    int first_fast_factor5_s = -1;
    double fast_chi_30 = NAN;
    double fast_chi_60 = NAN;
    double fast_chi_180 = NAN;
    double full_chi_30 = NAN;
    double slow_to_fast_30 = NAN;
    double slow_to_fast_60 = NAN;
};

CaptureReport evaluate(const OperatingPoint& op) {
    // These are the source-audited validated sensor coefficients used by the
    // current OU-III live-basin diagnostic and simulation adapter.
    const Vec3 sigma_a = Vec3::Constant(0.0294);
    const Vec3 sigma_g = Vec3::Constant(0.000157);
    const Vec3 sigma_m = Vec3::Constant(0.36);
    const Vec3 mag_world(20.0, 5.0, 44.0);

    Core f(sigma_a, sigma_g, sigma_m);
    f.set_mag_world_ref(mag_world);

    // Reproduce the deployed MahonyProxy -> goLive covariance sequence rather
    // than warming P at equilibrium.  Cold mode freezes the linear block and
    // b_a update; handoff seeds anisotropic attitude covariance, applies the
    // staged OU/R_S point, resets P_aw to stationary, then enables the block.
    f.set_linear_block_enabled(false);
    f.set_acc_bias_updates_enabled(false);
    f.initialize_from_attitude(Eigen::Quaterniond::Identity(), 0.035, 0.087);
    f.set_aw_time_constant(op.tau);
    f.set_aw_stationary_std(Vec3::Constant(op.sigma_aw));
    f.set_RS_noise(Vec3::Constant(op.rS));
    const double period = std::clamp((0.015 / 1.1) * op.tau, 0.005, 0.25);
    f.set_pseudo_update_period_s(period);
    f.reset_aw_covariance_to_stationary();
    f.set_linear_block_enabled(true);
    f.set_acc_bias_updates_enabled(false);

    const Mat21 P0 = f.Pext;
    Mat21 Psi = Mat21::Identity();

    CaptureReport r;
    for (int k = 0; k < MAX_SECONDS * 200; ++k) {
        Psi = step(f, k, mag_world) * Psi;
        if ((k + 1) % 200 != 0) continue;
        const int sec = (k + 1) / 200;
        const double cf = metric_norm_fast(Psi, P0, f.Pext);
        if (r.first_fast_factor5_s < 0 && cf <= 0.2)
            r.first_fast_factor5_s = sec;
        if (sec == 30) {
            r.fast_chi_30 = cf;
            r.full_chi_30 = metric_norm_full(Psi, P0, f.Pext);
            r.slow_to_fast_30 = metric_gain_slow_to_fast(Psi, P0, f.Pext);
        }
        if (sec == 60) {
            r.fast_chi_60 = cf;
            r.slow_to_fast_60 = metric_gain_slow_to_fast(Psi, P0, f.Pext);
        }
        if (sec == 180) r.fast_chi_180 = cf;
    }
    return r;
}

void print_one(const char* family, const OperatingPoint& op) {
    const CaptureReport r = evaluate(op);
    std::cout << family << ',' << op.name << ','
              << op.tau << ',' << op.sigma_aw << ',' << op.rS << ','
              << r.first_fast_factor5_s << ','
              << r.fast_chi_30 << ',' << r.fast_chi_60 << ',' << r.fast_chi_180 << ','
              << r.full_chi_30 << ','
              << r.slow_to_fast_30 << ',' << r.slow_to_fast_60 << '\n';
}

} // namespace

int main() {
    std::cout << std::setprecision(12);
    std::cout << "HANDOFF_CAPTURE_DIAGNOSTIC goLive_no_covariance_warmup=1 ba_locked=1\n";
    std::cout << "family,name,tau,sigma_aw,rS,first_fast_factor5_s,fast_chi30,fast_chi60,fast_chi180,full_chi30,slow_to_fast30,slow_to_fast60\n";
    for (const auto& p : kReferencePoints) print_one("REF", p);
    for (const auto& p : source_points()) print_one("SRC", p);
    return 0;
}
