/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Checks of the IMU lever arm of the experimental CARMA filters
  (kalman_common/ImuLeverArm.h, Kalman3D_Wave_CARMA3, Kalman3D_Wave_TFG_CARMA3
  and their wrapper; docs/carma-imu-lever-arm.md):

    * the band-limited model reconstructs the tangential and centripetal
      terms of a rocking body within its low-pass lag, and its gyro-bias
      Jacobian is the central difference of the term in a constant bias
      shift (the derivative of a constant is zero; only the centripetal
      term moves);
    * the OU-III ablation is OU-III's alpha recursion and bias Jacobian;
    * a zero offset is the filter without a lever arm, bit for bit, for both
      inner filters and the wrapper;
    * on a noiseless rigid body rocking in roll and pitch with the IMU 1 m
      above and 0.5 m ahead of the centre of gravity, both inner filters
      (driven directly) and both wrappers, corrected, return to within a
      small margin of their attitude and position error at the CoG, while the
      uncorrected filters do not.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdio>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
#include "kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h"

const float g_std = 9.80665f;

using Vec3 = Eigen::Vector3f;
using Vec3d = Eigen::Vector3d;
using Mat3d = Eigen::Matrix3d;
using ocean_imu::kalman::filter_detail::ImuLeverArm;

namespace {

int g_failures = 0, g_checks = 0;
void check(bool ok, const char* what, double v = NAN, double lim = NAN) {
    ++g_checks;
    if (ok) {
        std::printf("PASS  %s\n", what);
        return;
    }
    ++g_failures;
    if (std::isnan(v)) std::printf("FAIL  %s\n", what);
    else std::printf("FAIL  %s (value %.6g, limit %.6g)\n", what, v, lim);
}

constexpr double kTwoPi = 6.283185307179586;
constexpr double kG = 9.80665;

// --- the model ---------------------------------------------------------------

void test_band_model_kinematics() {
    // Rocking about a tilted axis at 0.3 Hz, 200 Hz samples.
    const double f = 0.3, A = 0.25, dt = 0.005;
    const Vec3d axis = Vec3d(0.8, 0.5, 0.33).normalized();
    const Vec3d r(0.4, -0.3, -1.0);
    ImuLeverArm<double> lever;
    lever.set_offset_body(r);
    double se_t = 0, sr_t = 0, se_c = 0, sr_c = 0;
    for (int k = 0; k < 4000; ++k) {
        const double t = k * dt;
        const Vec3d omega = A * std::sin(kTwoPi * f * t) * axis;
        const Vec3d alpha = A * kTwoPi * f * std::cos(kTwoPi * f * t) * axis;
        lever.update_rate(omega, dt);
        if (t < 2.0) continue;
        const Vec3d tan_true = alpha.cross(r), cen_true = omega.cross(omega.cross(r));
        const Vec3d tan_est = lever.angular_acceleration().cross(r);
        const Vec3d w = lever.angular_rate();
        const Vec3d cen_est = w.cross(w.cross(r));
        se_t += (tan_est - tan_true).squaredNorm(); sr_t += tan_true.squaredNorm();
        se_c += (cen_est - cen_true).squaredNorm(); sr_c += cen_true.squaredNorm();
    }
    // Two poles at 15 Hz lag 0.3 Hz by 2 f / f_c = 0.04 rad.
    const double et = std::sqrt(se_t / sr_t), ec = std::sqrt(se_c / sr_c);
    check(et < 0.05, "band model: tangential term within its low-pass lag (rocking at 0.3 Hz)", et, 0.05);
    check(ec < 0.10, "band model: centripetal term within its low-pass lag", ec, 0.10);
    std::printf("[model] relative RMS error tangential %.4f centripetal %.4f\n", et, ec);
}

void test_band_model_bias_jacobian() {
    const double dt = 0.005, delta = 1e-6;
    const Vec3d r(0.4, -0.3, -1.0);
    std::mt19937 rng(7);
    std::normal_distribution<double> n(0.0, 0.2);
    std::vector<Vec3d> rates;
    for (int k = 0; k < 300; ++k) rates.emplace_back(0.3 * std::sin(0.02 * k) + 0.05 * n(rng),
                                                     0.2 * std::cos(0.03 * k), 0.1 + 0.05 * n(rng));
    ImuLeverArm<double> base;
    base.set_offset_body(r);
    for (const auto& w : rates) base.update_rate(w, dt);
    const Mat3d J = base.bias_jacobian();
    double worst = 0;
    for (int c = 0; c < 3; ++c) {
        // The filter subtracts its bias estimate: rate = gyr - b.
        ImuLeverArm<double> plus, minus;
        plus.set_offset_body(r);
        minus.set_offset_body(r);
        for (const auto& w : rates) {
            plus.update_rate(w - delta * Vec3d::Unit(c), dt);
            minus.update_rate(w + delta * Vec3d::Unit(c), dt);
        }
        const Vec3d fd = (plus.acceleration() - minus.acceleration()) / (2 * delta);
        worst = std::max(worst, (fd - J.col(c)).cwiseAbs().maxCoeff() / (J.cwiseAbs().maxCoeff() + 1e-12));
    }
    check(worst < 1e-6, "band model: gyro-bias Jacobian equals the central difference in a bias shift", worst, 1e-6);
}

void test_ou3_model() {
    // OU-III's Kalman3D_Wave_OU_III::time_update alpha and its J_bg.
    const double dt = 0.005, tau = 0.05;
    const Vec3d r(0.2, 0.7, -0.4);
    ImuLeverArm<double> lever;
    lever.set_model(ImuLeverArm<double>::Model::OuIII);
    lever.set_offset_body(r);
    lever.set_alpha_smoothing_tau(tau);
    Vec3d prev = Vec3d::Zero(), alpha = Vec3d::Zero();
    bool have_prev = false;
    double worst = 0;
    for (int k = 0; k < 400; ++k) {
        const Vec3d w(0.3 * std::sin(0.05 * k), -0.2 * std::cos(0.03 * k), 0.05 * std::sin(0.11 * k));
        lever.update_rate(w, dt);
        if (have_prev) {
            const double a = 1.0 - std::exp(-dt / tau);
            alpha = (1.0 - a) * alpha + a * (w - prev) / dt;
        } else {
            alpha.setZero();
            have_prev = true;
        }
        prev = w;
        const Vec3d expect = alpha.cross(r) + w.cross(w.cross(r));
        worst = std::max(worst, (lever.acceleration() - expect).cwiseAbs().maxCoeff());
    }
    check(worst < 1e-12, "OU-III ablation: OU-III's alpha recursion and lever term", worst, 1e-12);
    const Vec3d w = prev;
    Mat3d rx;
    rx << 0, -r(2), r(1), r(2), 0, -r(0), -r(1), r(0), 0;
    const double k_alpha = (1.0 - std::exp(-dt / tau)) / dt;
    const Mat3d expect_J = k_alpha * rx - (Mat3d::Identity() * w.dot(r) + w * r.transpose() - 2.0 * r * w.transpose());
    const double ej = (lever.bias_jacobian() - expect_J).cwiseAbs().maxCoeff();
    check(ej < 1e-12, "OU-III ablation: OU-III's gyro-bias Jacobian", ej, 1e-12);
}

void test_zero_offset_disables() {
    ImuLeverArm<float> lever;
    lever.set_offset_body(Vec3(0.1f, 0.0f, 0.0f));
    check(lever.enabled(), "a nonzero offset enables the term");
    lever.set_offset_body(Vec3::Zero());
    check(!lever.enabled(), "a zero offset disables the term");
    lever.set_offset_body(Vec3(NAN, 0.0f, 0.0f));
    check(!lever.enabled() && lever.offset_body().allFinite(), "a non-finite offset is ignored");
}

// --- a rigid body rocking off the CoG ----------------------------------------

// Body->world NED attitude: fixed heading, roll and pitch sinusoids.  World
// CoG motion: three sinusoids.  Rates and angular accelerations by central
// differences of the attitude (double, h = 1e-4 s).
struct RigidBody {
    double yaw0 = 0.7;
    double roll0 = 0.16, fr = 0.21, pitch0 = 0.07, fp = 0.33, pp = 0.9;
    double A[3] = {0.35, 0.25, 0.9}, ph[3] = {0.4, 1.3, 0.0};
    double fw = 0.15;

    Mat3d R(double t) const {
        const double roll = roll0 * std::sin(kTwoPi * fr * t);
        const double pitch = pitch0 * std::sin(kTwoPi * fp * t + pp);
        return (Eigen::AngleAxisd(yaw0, Vec3d::UnitZ()) * Eigen::AngleAxisd(pitch, Vec3d::UnitY()) *
                Eigen::AngleAxisd(roll, Vec3d::UnitX())).toRotationMatrix();
    }
    Vec3d omega(double t) const {
        const double h = 1e-4;
        const Eigen::AngleAxisd d(R(t - h).transpose() * R(t + h));
        return d.axis() * d.angle() / (2 * h);
    }
    Vec3d alpha(double t) const {
        const double h = 1e-4;
        return (omega(t + h) - omega(t - h)) / (2 * h);
    }
    Vec3d pos(double t) const {
        Vec3d p;
        for (int i = 0; i < 3; ++i) p(i) = A[i] * std::sin(kTwoPi * fw * t + ph[i]);
        return p;
    }
    Vec3d acc(double t) const {
        const double w2 = kTwoPi * fw * kTwoPi * fw;
        Vec3d a;
        for (int i = 0; i < 3; ++i) a(i) = -w2 * A[i] * std::sin(kTwoPi * fw * t + ph[i]);
        return a;
    }
    // Specific force at body offset r.
    Vec3d specific_force(double t, const Vec3d& r) const {
        const Vec3d f_cog = R(t).transpose() * (acc(t) - Vec3d(0, 0, kG));
        const Vec3d w = omega(t);
        return f_cog + alpha(t).cross(r) + w.cross(w.cross(r));
    }
};

struct Scores {
    double tilt_deg = 0, pos_m = 0;
    Vec3 p_end = Vec3::Zero();
    Eigen::Quaternionf q_end = Eigen::Quaternionf::Identity();
};

double tilt_error_deg(const Mat3d& R_true, const Eigen::Quaternionf& q_est) {
    const Vec3d down_true = R_true.transpose() * Vec3d::UnitZ();
    const Vec3d down_est = q_est.cast<double>().toRotationMatrix().transpose() * Vec3d::UnitZ();
    return std::acos(std::clamp(down_true.dot(down_est), -1.0, 1.0)) * 57.29577951308232;
}

enum class Mount { CoG, OffUncorrected, OffCorrected, ZeroSet };

const Vec3d kOffsetNed(0.5, 0.0, -1.0);  // 0.5 m ahead of and 1 m above the CoG
const Vec3 kMagWorld(20.0f, -4.0f, 45.0f);

// The inner filter driven directly, as carma3-test's sinusoid.
template<class F>
Scores run_inner(Mount mount) {
    const RigidBody body;
    F f(Vec3::Constant(0.03f), Vec3::Constant(3e-5f), Vec3::Constant(8.0f));
    const float kappa = static_cast<float>(body.fw / f.vertical_model().base_zero_crossing_hz());
    f.set_vertical_carma(kappa, 0.6f);
    f.set_horizontal_carma(kappa, 0.5f, 0.4f, 0.0f);
    f.reset_aw_covariance_to_stationary();
    f.set_acc_bias_updates_enabled(false);
    const Eigen::Quaternionf q0(body.R(0.0).cast<float>());
    f.initialize_from_attitude(q0, 0.01f, 0.05f);
    f.set_mag_world_ref(kMagWorld);
    const bool off = mount == Mount::OffUncorrected || mount == Mount::OffCorrected;
    const Vec3d r_install = off ? kOffsetNed : Vec3d::Zero();
    if (mount == Mount::OffCorrected) f.set_imu_lever_arm_body(kOffsetNed.cast<float>());
    if (mount == Mount::ZeroSet) f.set_imu_lever_arm_body(Vec3::Zero());
    const double dt = 0.005;
    double st = 0, sp = 0;
    long n = 0;
    for (int k = 0; k < 200 * 150; ++k) {
        const double t = k * dt;
        f.time_update(body.omega(t).cast<float>(), static_cast<float>(dt));
        f.measurement_update_acc_only(body.specific_force(t, r_install).cast<float>());
        if (k % 8 == 0) f.measurement_update_mag_only((body.R(t).transpose() * kMagWorld.cast<double>()).cast<float>());
        if (t > 60.0) {
            st += std::pow(tilt_error_deg(body.R(t), f.quaternion_boat()), 2);
            sp += (f.get_position().template cast<double>() - body.pos(t)).squaredNorm();
            ++n;
        }
    }
    Scores s;
    s.tilt_deg = std::sqrt(st / static_cast<double>(n));
    s.pos_m = std::sqrt(sp / static_cast<double>(n));
    s.p_end = f.get_position();
    s.q_end = f.quaternion_boat();
    return s;
}

// The full wrapper: proxy startup, tuner, magnetic acquisition, small sensor
// noise.  The correction is the wrapper's input-stage lever arm.
template<class Fusion>
Scores run_wrapper(Mount mount) {
    const RigidBody body;
    Fusion fusion;
    typename Fusion::Config cfg;
    cfg.with_mag = true;
    cfg.sigma_a = Vec3::Constant(2.8f * 1.51e-3f * 9.80665f * 0.71f);
    cfg.sigma_g = Vec3::Constant(2.0f * 0.00157f * 0.01f);
    cfg.sigma_m = Vec3::Constant(1.2f * 0.8f * 8.0f);
    fusion.begin(cfg);
    fusion.raw().enableTuner(true);
    fusion.raw().enableClamp(true);
    const bool off = mount == Mount::OffUncorrected || mount == Mount::OffCorrected;
    const Vec3d r_install = off ? kOffsetNed : Vec3d::Zero();
    if (mount == Mount::OffCorrected) fusion.raw().setImuLeverArmBody(kOffsetNed.cast<float>());
    if (mount == Mount::ZeroSet) fusion.raw().setImuLeverArmBody(Vec3::Zero());
    std::mt19937 rng(11);
    std::normal_distribution<float> n(0.0f, 1.0f);
    const double dt = 0.005;
    double st = 0, sp = 0;
    long cnt = 0;
    for (int k = 0; k < 200 * 400; ++k) {
        const double t = k * dt;
        const Vec3 acc = body.specific_force(t, r_install).cast<float>() + 0.01f * Vec3(n(rng), n(rng), n(rng));
        const Vec3 gyr = body.omega(t).cast<float>() + 2e-4f * Vec3(n(rng), n(rng), n(rng));
        fusion.update(static_cast<float>(dt), gyr, acc, 35.0f);
        if (k % 8 == 0) {
            const Vec3 m = (body.R(t).transpose() * kMagWorld.cast<double>()).cast<float>();
            fusion.updateMag(m + 0.3f * Vec3(n(rng), n(rng), n(rng)));
        }
        if (t > 250.0 && fusion.isLive()) {
            st += std::pow(tilt_error_deg(body.R(t), fusion.attitudeQuat()), 2);
            sp += (fusion.raw().mekf().get_position().template cast<double>() - body.pos(t)).squaredNorm();
            ++cnt;
        }
    }
    Scores s;
    if (cnt == 0) return Scores{NAN, NAN};
    s.tilt_deg = std::sqrt(st / static_cast<double>(cnt));
    s.pos_m = std::sqrt(sp / static_cast<double>(cnt));
    s.p_end = fusion.raw().mekf().get_position();
    s.q_end = fusion.attitudeQuat();
    return s;
}

void compare(const char* label, const Scores& cog, const Scores& off, const Scores& fixed, const Scores& zero,
             double tilt_margin, double pos_margin) {
    char buf[200];
    std::printf("[%s] tilt RMS deg: CoG %.4f  off-CoG %.4f  corrected %.4f | position RMS m: CoG %.4f  off-CoG %.4f  "
                "corrected %.4f\n", label, cog.tilt_deg, off.tilt_deg, fixed.tilt_deg, cog.pos_m, off.pos_m, fixed.pos_m);
    std::snprintf(buf, sizeof buf, "%s: the uncorrected off-CoG IMU degrades position by more than 1.5x", label);
    check(off.pos_m > 1.5 * cog.pos_m, buf, off.pos_m / cog.pos_m, 1.5);
    std::snprintf(buf, sizeof buf, "%s: corrected tilt within %.3f deg of the CoG run", label, tilt_margin);
    check(fixed.tilt_deg < cog.tilt_deg + tilt_margin, buf, fixed.tilt_deg, cog.tilt_deg + tilt_margin);
    std::snprintf(buf, sizeof buf, "%s: corrected position within %.3f m of the CoG run", label, pos_margin);
    check(fixed.pos_m < cog.pos_m + pos_margin, buf, fixed.pos_m, cog.pos_m + pos_margin);
    std::snprintf(buf, sizeof buf, "%s: a zero lever arm is the CoG run bit for bit", label);
    check(zero.p_end == cog.p_end && zero.q_end.coeffs() == cog.q_end.coeffs(), buf);
}

template<class F>
void test_inner(const char* label) {
    compare(label, run_inner<F>(Mount::CoG), run_inner<F>(Mount::OffUncorrected), run_inner<F>(Mount::OffCorrected),
            run_inner<F>(Mount::ZeroSet), 0.02, 0.01);
}

template<class Fusion>
void test_wrapper(const char* label) {
    compare(label, run_wrapper<Fusion>(Mount::CoG), run_wrapper<Fusion>(Mount::OffUncorrected),
            run_wrapper<Fusion>(Mount::OffCorrected), run_wrapper<Fusion>(Mount::ZeroSet), 0.05, 0.02);
}

}  // namespace

int main() {
    test_band_model_kinematics();
    test_band_model_bias_jacobian();
    test_ou3_model();
    test_zero_offset_disables();
    test_inner<ocean_imu::kalman::Kalman3D_Wave_CARMA3<float, 4, false, 4, false>>("CARMA3 inner filter");
    test_inner<ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3<float, 4>>("TFG-CARMA3 inner filter");
    test_wrapper<ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false>>("CARMA3 wrapper");
    test_wrapper<ocean_imu::carma::SeaStateFusion_TFG_CARMA3<TrackerType::KALMANF, 4>>("TFG-CARMA3 wrapper");
    std::printf("carma_lever_arm-test: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
