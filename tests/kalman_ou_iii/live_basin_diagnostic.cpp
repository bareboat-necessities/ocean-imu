// Interval constants of the OU-III Live-basin certificate, evaluated on the
// eight committed SpectralMSE reference operating points.
//
// Phase 2 established the certificate in a fixed dimensionless Euclidean norm
// and reported that the deployed schedule contracts in the Riccati metric but
// not in that Euclidean norm over 30 s -- the 5000 s residual-bias OU mode
// produces a long non-normal transient.  Converting the metric contraction
// back into the fixed norm cost a factor sqrt(pbar/punder) ~ 566, and that
// factor is what made the resulting basin radius meaningless.
//
// Phase 3 carries the whole argument in the metric instead.  Two things fall
// out of that and are checked here:
//
//   1. The one-sample gain in the covariance metric is at most one, for every
//      sample, with no hypothesis at all.  This is the Joseph identity: with
//      A = (I-KC)F,
//        A P A^T = P^+ - K R K^T - (I-KC) Q (I-KC)^T  <=  P^+,
//      so ||L^+{-1} A L||_2 <= 1.  The prefix constant of the horizon
//      certificate is therefore exactly 1 rather than a measured envelope,
//      and M_H = rho_H^{-(H-1)} sits just above one instead of at ~600.
//      alpha_max below is the numerical witness.
//
//   2. The metric norms are invariant under the fixed diagonal scaling, since
//      ||D e||_{(D P D)^{-1}} = ||e||_{P^{-1}}.  The certificate therefore no
//      longer depends on the choice of physical scales at all.  The Phase-2
//      scaled quantities are still reported, unchanged, as the regression
//      witness they were.
//
// What the small-gain slope needs on top of that is the size of the nonlinear
// remainder measured in the same metric.  Bounding it with a single
// M_H/(1-rho_H) prices every injection at the slowest mode's memory, which at
// the deployed schedule is 3.45e5 samples.  The directional l1 injection gains
// below price each channel at its own instead: an attitude-reset remainder, an
// accelerometer residual remainder and a magnetometer residual remainder decay
// at very different rates, and nothing about the slowest mode describes the
// first two.
//
// Reference replay remains feasibility evidence.  An operating trajectory is
// certified only if it satisfies its own interval envelopes.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
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
using Mat21x3 = Eigen::Matrix<double, 21, 3>;

constexpr int OFF_TH = 0;
constexpr int OFF_BG = 3;
constexpr int OFF_V  = 6;
constexpr int OFF_P  = 9;
constexpr int OFF_S  = 12;
constexpr int OFF_AW = 15;
constexpr int OFF_BA = 18;
constexpr double DT = 1.0 / 200.0;
constexpr int HORIZON_STEPS = 30 * 200;

// Established-interval window.  The certificate's interval envelopes are the
// suprema over this window, taken after the seed covariance has settled; the
// settling transient itself is an explicit hypothesis of the theorem rather
// than something these numbers cover.  See the paper's discussion of the
// handoff-to-established map.
//
// The window length is part of the declaration, not a free parameter: the
// residual-bias marginal grows slowly over a Live interval, so a certificate
// issued against a 300 s envelope is a statement about a 300 s interval.
// --long widens both windows for offline work; CI runs the default.
constexpr double SETTLE_SEC = 30.0;

struct Horizons {
    int envelope_steps = 600 * 200;
    int impulse_steps  = 2000 * 200;
    int alpha_stride   = 1;
};

// The one-sample metric gain is bounded by an identity rather than by a
// hypothesis, so a stride would still be a sufficient regression witness for
// the identity holding in the implementation.  The default checks every
// sample because it is affordable at these horizons.

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

// Singular values only.  Thin U/V are invalid for fixed-size Eigen matrices
// and are unnecessary for an induced 2-norm.
template <int R, int C>
double norm2_fixed(const Eigen::Matrix<double, R, C>& A) {
    Eigen::JacobiSVD<Eigen::Matrix<double, R, C>> svd(A);
    return svd.singularValues()(0);
}

double norm2(const Mat21& A) { return norm2_fixed<21, 21>(A); }

Mat21 state_scale() {
    // Fixed Phase-2 proof coordinates.  These are design scales, not fitted
    // weights: magnetically gauged handoff yaw sigma, constructor gyro-bias
    // sigma, constructor v/p/S sigmas, wrapper a_w safety ceiling, and the
    // hard b_a projection radius respectively.  Phase 3 does not depend on
    // them; they are kept so the Phase-2 columns stay comparable.
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

Mat21 correction_transition(const Mat21x3& K, const Mat3x21& C) {
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
    Mat21x3 Ka = Mat21x3::Zero();
    Mat21x3 Km = Mat21x3::Zero();
    Mat21x3 KS = Mat21x3::Zero();
    bool pseudo = false;
    bool mag = false;
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
        out.KS = f.K_scratch_;
        A = correction_transition(out.KS, S_jacobian()) * A;
    }

    const Mat3x21 Ca = accel_jacobian(f);
    f.measurement_update_acc_only(acc, 35.0);
    out.Ka = f.K_scratch_;
    A = correction_transition(out.Ka, Ca) * A;

    // Reference simulations provide a 100 Hz magnetometer against 200 Hz IMU.
    if ((k & 1) == 0) {
        const Mat3x21 Cm = mag_jacobian(f);
        f.measurement_update_mag_only(mag_world);
        out.Km = f.K_scratch_;
        A = correction_transition(out.Km, Cm) * A;
        out.mag = true;
    }

    out.A = A;
    out.pseudo = pseudo;
    return out;
}

// Induced norm from the covariance/Riccati metric at the start of a lifted
// interval to the metric at its end.  If P0=L0 L0^T and P1=L1 L1^T, then
// ||Psi||_{P0^{-1}->P1^{-1}} = ||L1^{-1} Psi L0||_2.
double covariance_metric_norm(const Mat21& Psi, const Mat21& P0, const Mat21& P1) {
    Eigen::LLT<Mat21> llt0(P0);
    Eigen::LLT<Mat21> llt1(P1);
    if (llt0.info() != Eigen::Success || llt1.info() != Eigen::Success) {
        return std::numeric_limits<double>::infinity();
    }
    const Mat21 L0 = llt0.matrixL();
    const Mat21 B = llt1.matrixL().solve(Psi * L0);
    return norm2(B);
}

std::pair<double,double> scaled_cov_eigen_bounds(const Mat21& P, const Mat21& D) {
    const Mat21 Pz = D * P * D;
    Eigen::SelfAdjointEigenSolver<Mat21> es(Pz, Eigen::EigenvaluesOnly);
    if (es.info() != Eigen::Success) {
        return {0.0, std::numeric_limits<double>::infinity()};
    }
    return {es.eigenvalues()(0), es.eigenvalues()(20)};
}

double block_sigma(const Mat21& P, int off) {
    return std::sqrt(std::max(0.0, norm2_fixed<3,3>(Mat3(P.block<3,3>(off, off)))));
}

struct Report {
    // Phase-2 columns, unchanged.
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

    // Phase-3 columns.
    double alpha_max = 0.0;      // sup one-sample metric gain; must be <= 1
    double M_H = 0.0;            // rho_H^{-(H-1)}
    double gamma_theta = 0.0;    // l1 injection gain, attitude channel [1/rad]
    double gamma_acc = 0.0;      // [1/(m/s^2)]
    double gamma_mag = 0.0;      // [1/uT]
    double sigma_theta = 0.0;    // established-interval envelopes
    double sigma_bg = 0.0;
    double sigma_S = 0.0;
    double sigma_aw = 0.0;
    double sigma_ba = 0.0;
    double eth_Ka = 0.0;         // attitude rows of each accepted gain
    double eth_Km = 0.0;
    double eth_KS = 0.0;
    double c_eff = 0.0;          // small-gain slope, nu = c_eff * r
    double r_xi = 0.0;           // largest tube the remainder expansion is valid on
    double r_cert = 0.0;         // 1 / c_eff
    double budget = 0.0;         // max admissible ||e_H||_{V_H}
    double tail_share = 0.0;     // fraction of gamma_theta supplied by the bound
};

Report evaluate(const OperatingPoint& op, const Horizons& hz) {
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

    // Accepted-vector magnitude bounds of the ISS geometry hypotheses: the
    // specific force the accelerometer update may accept, and the field
    // magnitude the magnetometer update may accept.
    const double f1 = 9.80665 + 3.0 * op.sigma_aw;
    const double m1 = mag_world.norm();
    // Inverse-left-Jacobian constant of the group-composition remainder on
    // |x| <= theta_c; see the paper.  theta_c = 1 rad.
    const double theta_c = 1.0;
    const double j_c = 0.5 + (1.0 / (theta_c * theta_c)
                              - (1.0 + std::cos(theta_c))
                                    / (2.0 * theta_c * std::sin(theta_c)))
                             * theta_c;

    Report r;
    Mat21 P_prev = f.Pext;

    // ---- pass 1: settle, then take the established-interval envelopes ----
    const int settle_steps = int(SETTLE_SEC / DT);
    const int envelope_steps = hz.envelope_steps;
    double alpha_max = 0.0;
    for (int k = 0; k < settle_steps + envelope_steps; ++k) {
        const bool check_alpha = (k % hz.alpha_stride) == 0;
        Mat21 L0 = Mat21::Identity();
        bool have_L0 = false;
        if (check_alpha) {
            Eigen::LLT<Mat21> llt0(P_prev);
            have_L0 = (llt0.info() == Eigen::Success);
            if (have_L0) L0 = llt0.matrixL();
        }
        const StepResult s = step(f, k, mag_world);
        if (check_alpha && have_L0) {
            Eigen::LLT<Mat21> llt1(f.Pext);
            if (llt1.info() == Eigen::Success) {
                const Mat21 B = llt1.matrixL().solve(s.A * L0);
                alpha_max = std::max(alpha_max, norm2(B));
            }
        }
        if (k >= settle_steps) {
            r.sigma_theta = std::max(r.sigma_theta, block_sigma(f.Pext, OFF_TH));
            r.sigma_bg    = std::max(r.sigma_bg,    block_sigma(f.Pext, OFF_BG));
            r.sigma_S     = std::max(r.sigma_S,     block_sigma(f.Pext, OFF_S));
            r.sigma_aw    = std::max(r.sigma_aw,    block_sigma(f.Pext, OFF_AW));
            r.sigma_ba    = std::max(r.sigma_ba,    block_sigma(f.Pext, OFF_BA));
            r.eth_Ka = std::max(r.eth_Ka, norm2_fixed<3,3>(Mat3(s.Ka.block<3,3>(0,0))));
            r.eth_Km = std::max(r.eth_Km, norm2_fixed<3,3>(Mat3(s.Km.block<3,3>(0,0))));
            r.eth_KS = std::max(r.eth_KS, norm2_fixed<3,3>(Mat3(s.KS.block<3,3>(0,0))));
        }
        P_prev = f.Pext;
    }
    r.alpha_max = alpha_max;

    // ---- pass 2: Phase-2 horizon quantities on the established interval ----
    const Mat21 S = state_scale();
    Mat21 D = Mat21::Zero();
    for (int i = 0; i < 21; ++i) D(i,i) = 1.0 / S(i,i);

    const Mat21 P0 = f.Pext;
    Mat21 Psi_phys = Mat21::Identity();
    Mat21 Psi_scaled = Mat21::Identity();
    double prefix_euclid = 1.0;
    int pseudo_count = 0;
    auto [pz_min, pz_max] = scaled_cov_eigen_bounds(P0, D);

    // Impulse states for the directional l1 gains, launched from this sample.
    Eigen::Matrix<double,21,9> Z;
    Z.setZero();
    for (int j = 0; j < 3; ++j) Z(OFF_TH + j, j) = 1.0;   // 1 rad of attitude remainder
    bool have_Ka = false, have_Km = false;

    for (int k = 0; k < HORIZON_STEPS; ++k) {
        const StepResult sr = step(f, k, mag_world);
        if (sr.pseudo) ++pseudo_count;
        if (!have_Ka) { Z.block<21,3>(0,3) = sr.Ka; have_Ka = true; }
        if (!have_Km && sr.mag) { Z.block<21,3>(0,6) = sr.Km; have_Km = true; }

        Psi_phys = sr.A * Psi_phys;
        const Mat21 Abar = D * sr.A * S;
        Psi_scaled = Abar * Psi_scaled;
        prefix_euclid = std::max(prefix_euclid, norm2(Psi_scaled));

        const auto [lo, hi] = scaled_cov_eigen_bounds(f.Pext, D);
        pz_min = std::min(pz_min, lo);
        pz_max = std::max(pz_max, hi);
    }

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
        r.rho_metric_sample = std::pow(r.chi_metric, 1.0 / double(HORIZON_STEPS));
        r.rho_metric_second = std::pow(r.chi_metric, 1.0 / 30.0);
        // Prefix bound is exactly 1 by the metric monotonicity lemma, so the
        // horizon certificate's transition constant is just rho^{-(H-1)}.
        r.M_H = std::pow(r.rho_metric_sample, -(HORIZON_STEPS - 1));
    }

    {
        std::array<int, 15> xi{{0,1,2,3,4,5,12,13,14,15,16,17,18,19,20}};
        std::array<int, 6> ell{{6,7,8,9,10,11}};
        Eigen::Matrix<double,15,15> Pxx;
        Eigen::Matrix<double,15,6> Pxl;
        Eigen::Matrix<double,15,21> Px;
        for (int i = 0; i < 15; ++i) {
            for (int j = 0; j < 15; ++j) Pxx(i,j) = Psi_scaled(xi[size_t(i)], xi[size_t(j)]);
            for (int j = 0; j < 6; ++j)  Pxl(i,j) = Psi_scaled(xi[size_t(i)], ell[size_t(j)]);
            for (int j = 0; j < 21; ++j) Px(i,j)  = Psi_scaled(xi[size_t(i)], j);
        }
        r.xi_xi = norm2_fixed<15,15>(Pxx);
        r.xi_ell = norm2_fixed<15,6>(Pxl);
        r.xi_all = norm2_fixed<15,21>(Px);
    }

    // ---- pass 3: directional l1 injection gains in the metric --------------
    // Gamma_g = sum_k || Psi(k, i+1) N_g ||_{V_k}: how much of one unit of
    // remainder injected in channel g the metric still carries, summed over
    // all later samples.  The scalar bound M_H/(1-rho_H) replaces every one of
    // these with the slowest mode's memory.
    const int impulse_steps = hz.impulse_steps;
    double tail_theta = 0.0, tail_acc = 0.0, tail_mag = 0.0;
    for (int k = 0; k < impulse_steps; ++k) {
        const StepResult sr = step(f, HORIZON_STEPS + k, mag_world);
        Z = sr.A * Z;
        Eigen::LLT<Mat21> llt(f.Pext);
        if (llt.info() != Eigen::Success) break;
        const Eigen::Matrix<double,21,9> Y = llt.matrixL().solve(Z);
        tail_theta = norm2_fixed<21,3>(Mat21x3(Y.block<21,3>(0,0)));
        tail_acc   = norm2_fixed<21,3>(Mat21x3(Y.block<21,3>(0,3)));
        tail_mag   = norm2_fixed<21,3>(Mat21x3(Y.block<21,3>(0,6)));
        r.gamma_theta += tail_theta;
        r.gamma_acc   += tail_acc;
        r.gamma_mag   += tail_mag;
    }

    // Everything beyond the computed horizon is bounded rather than dropped.
    // Psi(k, i+1) factors through Psi(k, i+1+K), so Theorem "metric
    // finite-horizon UES" gives
    //   sum_{k>K} ||Psi(k,i+1) N||_{V_k} <= ||Z_K||_V * M_H rho_H / (1-rho_H),
    // which is a genuine upper bound and not an extrapolation.  Truncating
    // without it would understate the injection gains, which is the unsafe
    // direction.
    if (r.rho_metric_sample > 0.0 && r.rho_metric_sample < 1.0) {
        const double tail_gain =
            r.M_H * r.rho_metric_sample / (1.0 - r.rho_metric_sample);
        r.gamma_theta += tail_theta * tail_gain;
        r.gamma_acc   += tail_acc * tail_gain;
        r.gamma_mag   += tail_mag * tail_gain;
        r.tail_share = (r.gamma_theta > 0.0)
                           ? (tail_theta * tail_gain) / r.gamma_theta : 0.0;
    }

    // ---- small-gain slope --------------------------------------------------
    // Inside the metric tube ||e_k||_{V_k} <= R every block obeys
    // ||delta x_b|| <= sigma_b R, so each remainder is a product of two such
    // bounds and is quadratic in R with the coefficients below.
    const double A_g = r.eth_Ka * (f1 * r.sigma_theta + r.sigma_aw + r.sigma_ba)
                     + r.eth_Km * m1 * r.sigma_theta
                     + r.eth_KS * r.sigma_S
                     + DT * r.sigma_bg;
    const double c_theta = j_c * r.sigma_theta * A_g;
    const double c_acc   = 0.5 * f1 * r.sigma_theta * r.sigma_theta
                         + r.sigma_theta * r.sigma_aw;
    const double c_mag   = 0.5 * m1 * r.sigma_theta * r.sigma_theta;

    // Largest tube the rotation-vector expansion is valid on: inside it both
    // the attitude error and the correction angle it composes with have to
    // stay under theta_c, since that is where the inverse left Jacobian's
    // bound j_c was taken.
    r.r_xi = theta_c / (r.sigma_theta + A_g);

    r.c_eff = r.gamma_theta * c_theta + r.gamma_acc * c_acc + r.gamma_mag * c_mag;
    if (r.c_eff > 0.0) {
        r.r_cert = std::min(1.0 / r.c_eff, r.r_xi);
        // (1 - c_eff r) r is maximised at r = 1/(2 c_eff), capped at r_xi.
        const double r_opt = std::min(0.5 / r.c_eff, r.r_xi);
        r.budget = (1.0 - r.c_eff * r_opt) * r_opt / std::max(1.0, r.M_H);
    }
    return r;
}

}  // namespace

int main(int argc, char** argv) {
    // CI runs the default horizons.  --long widens both windows for offline
    // work; it changes the interval the envelopes describe, so the two are not
    // interchangeable and the certificate constants are quoted for the
    // default.
    Horizons hz;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--long") {
            hz.envelope_steps = 1800 * 200;
            hz.impulse_steps = 4000 * 200;
        }
    }

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

    double worst_alpha = 0.0;
    double worst_M_H = 0.0;
    double worst_c_eff = 0.0;
    double worst_gamma_theta = 0.0;
    double worst_gamma_acc = 0.0;
    double worst_gamma_mag = 0.0;
    double smallest_budget = std::numeric_limits<double>::infinity();
    double smallest_r_cert = std::numeric_limits<double>::infinity();
    double smallest_r_xi = std::numeric_limits<double>::infinity();

    std::array<Report, 8> reports{};

    for (size_t i = 0; i < kReferencePoints.size(); ++i) {
        const OperatingPoint& op = kReferencePoints[i];
        const Report r = evaluate(op, hz);
        reports[i] = r;
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

        worst_alpha = std::max(worst_alpha, r.alpha_max);
        worst_M_H = std::max(worst_M_H, r.M_H);
        worst_c_eff = std::max(worst_c_eff, r.c_eff);
        worst_gamma_theta = std::max(worst_gamma_theta, r.gamma_theta);
        worst_gamma_acc = std::max(worst_gamma_acc, r.gamma_acc);
        worst_gamma_mag = std::max(worst_gamma_mag, r.gamma_mag);
        smallest_budget = std::min(smallest_budget, r.budget);
        smallest_r_cert = std::min(smallest_r_cert, r.r_cert);
        smallest_r_xi = std::min(smallest_r_xi, r.r_xi);
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

    std::cout << "PHASE3_POINTS name,alpha_max,M_H,gamma_theta,gamma_acc,gamma_mag,"
                 "sigma_theta,sigma_bg,sigma_S,sigma_aw,sigma_ba,"
                 "EthKa,EthKm,EthKS,r_xi,c_eff,r_cert,budget,gamma_theta_tail_share\n";
    for (size_t i = 0; i < kReferencePoints.size(); ++i) {
        const Report& r = reports[i];
        std::cout << "PHASE3_ROW " << kReferencePoints[i].name
                  << ',' << r.alpha_max << ',' << r.M_H
                  << ',' << r.gamma_theta << ',' << r.gamma_acc << ',' << r.gamma_mag
                  << ',' << r.sigma_theta << ',' << r.sigma_bg << ',' << r.sigma_S
                  << ',' << r.sigma_aw << ',' << r.sigma_ba
                  << ',' << r.eth_Ka << ',' << r.eth_Km << ',' << r.eth_KS
                  << ',' << r.r_xi
                  << ',' << r.c_eff << ',' << r.r_cert << ',' << r.budget
                  << ',' << r.tail_share << '\n';
    }

    std::cout << "PHASE3_SUMMARY worst_alpha=" << worst_alpha
              << " worst_M_H=" << worst_M_H
              << " worst_gamma_theta=" << worst_gamma_theta
              << " worst_gamma_acc=" << worst_gamma_acc
              << " worst_gamma_mag=" << worst_gamma_mag
              << " worst_c_eff=" << worst_c_eff
              << " smallest_r_xi=" << smallest_r_xi
              << " smallest_r_cert=" << smallest_r_cert
              << " smallest_budget=" << smallest_budget
              << " scalar_l1_gain=" << (rho_sample > 0.0 ? worst_M_H / (1.0 - rho_sample) : 0.0)
              << '\n';

    // The fixed Euclidean norm is deliberately reported but not required to
    // contract over 30 s: the 5000 s residual-bias OU mode produces a long
    // non-normal transient.  Phase 3 carries the certificate in the metric and
    // never converts back, so that column is a witness rather than a
    // requirement.
    if (!all_metric_contract) {
        std::cerr << "FAIL: at least one deployed reference operating point did not contract in the Riccati metric over 30 s\n";
        return 1;
    }

    // The metric monotonicity lemma is an identity, not a hypothesis.  A
    // one-sample metric gain above one would mean the Joseph update in the
    // implementation is no longer the Joseph update the lemma is about, which
    // invalidates M_H and everything downstream of it.
    if (!(worst_alpha <= 1.0 + 1e-9)) {
        std::cerr << "FAIL: one-sample covariance-metric gain exceeded 1 (alpha_max="
                  << worst_alpha << "); the metric monotonicity lemma no longer holds\n";
        return 1;
    }

    // M_H is the constant Phase 2 had to take as ~600 after converting the
    // metric contraction back to fixed Euclidean coordinates.  Losing that
    // improvement would silently restore the meaningless Phase-2 radius.
    if (!(worst_M_H < 2.0)) {
        std::cerr << "FAIL: metric transition constant M_H=" << worst_M_H
                  << " is no longer close to one\n";
        return 1;
    }

    if (!(worst_c_eff > 0.0) || !std::isfinite(worst_c_eff)) {
        std::cerr << "FAIL: small-gain slope is not a finite positive number\n";
        return 1;
    }

    return 0;
}
