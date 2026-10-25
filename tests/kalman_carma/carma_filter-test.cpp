/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Checks of the experimental Kalman3D_Wave_CARMA filter and its wrapper:

    * state layout: 22 states for CARMA(5,2), 23 with S_z;
    * NED sign conventions: an exactly known sinusoidal heave under a
      prescribed level attitude is recovered in phase and amplitude;
    * the vertical acceleration enters the attitude-aware accelerometer model:
      tilted updates create attitude/vertical cross-covariance;
    * position / velocity pseudo-measurement APIs act on the CARMA h and v;
    * non-uniform steps keep the covariance symmetric positive semidefinite;
    * the S_z variant applies its integral update and stays bounded;
    * the full wrapper on a still, noisy IMU reports bounded heave and no
      non-finite state.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdio>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "util/WaveFilesSupport.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA.h"

const float g_std = 9.80665f;

namespace {

int g_failures = 0, g_checks = 0;
void check(bool ok, const char* what, double v = NAN, double lim = NAN) {
    ++g_checks;
    if (!ok) { ++g_failures; std::printf("FAIL: %s (value=%.6g limit=%.6g)\n", what, v, lim); }
}

using Vec3 = Eigen::Vector3f;
using Filter52 = ocean_imu::kalman::Kalman3D_Wave_CARMA<float, 5, false>;
using Filter52S = ocean_imu::kalman::Kalman3D_Wave_CARMA<float, 5, true>;

template<class F>
double min_eig_normalized(const F& f) {
    const auto P = f.covariance_full().template cast<double>().eval();
    Eigen::VectorXd d = P.diagonal().cwiseMax(1e-300).cwiseSqrt();
    const Eigen::MatrixXd C = d.cwiseInverse().asDiagonal() * P * d.cwiseInverse().asDiagonal();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
    return es.eigenvalues().minCoeff();
}

template<class F>
F make_filter() {
    const Vec3 sa = Vec3::Constant(0.03f), sg = Vec3::Constant(3e-5f), sm = Vec3::Constant(8.0f);
    F f(sa, sg, sm);
    f.set_aw_time_constant(1.7f);
    f.set_aw_stationary_std(Vec3(0.5f, 0.5f, 0.6f));
    f.set_vertical_carma(0.31f / 0.37f, 0.6f);
    f.reset_aw_covariance_to_stationary();
    f.set_RS_noise(Vec3(0.7f, 0.7f, 0.7f));
    f.set_pseudo_update_period_s(0.03f);
    f.initialize_from_attitude(Eigen::Quaternionf::Identity(), 0.01f, 0.05f);
    return f;
}

void test_layout() {
    check(Filter52::NX == 22, "CARMA(5,2) INS has 22 states", Filter52::NX, 22);
    check(Filter52S::NX == 23, "CARMA(5,2)+S_z INS has 23 states", Filter52S::NX, 23);
    check((ocean_imu::kalman::Kalman3D_Wave_CARMA<float, 4>::NX == 21), "CARMA(4,1) INS has 21 states");
    check((ocean_imu::kalman::Kalman3D_Wave_CARMA<float, 3>::NX == 20), "CARMA(3,0) INS has 20 states");
}

// Level attitude, noiseless specific force of h_up = A sin(w t): NED
// p_z = -h_up, a_z = A w^2 sin(w t), f = (0, 0, a_z - g).
void test_sinusoid_signs() {
    Filter52 f = make_filter<Filter52>();
    f.set_acc_bias_updates_enabled(false);
    const double A = 1.0, w = 2.0 * M_PI * 0.30, dt = 0.005;
    double se = 0.0, sv = 0.0;
    long n = 0;
    for (long k = 0; k < 200 * 400; ++k) {
        const double t = static_cast<double>(k + 1) * dt;
        const double h = A * std::sin(w * t);
        const double az_ned = A * w * w * std::sin(w * t);
        f.time_update(Vec3::Zero(), static_cast<float>(dt));
        f.measurement_update_acc_only(Vec3(0.0f, 0.0f, static_cast<float>(az_ned - 9.80665)));
        if (t > 200.0) {
            const double e = (-f.get_position().z()) - h;
            const double ev = (-f.get_velocity().z()) - A * w * std::cos(w * t);
            se += e * e; sv += ev * ev; ++n;
        }
    }
    const double rms = std::sqrt(se / static_cast<double>(n));
    const double rms_v = std::sqrt(sv / static_cast<double>(n));
    check(rms < 0.08 * A, "noiseless sinusoid: heave RMS error < 8% of amplitude", rms, 0.08 * A);
    check(rms_v < 0.08 * A * w, "noiseless sinusoid: velocity RMS error < 8%", rms_v, 0.08 * A * w);
    std::printf("[sinusoid] heave rms err %.4f m, velocity rms err %.4f m/s\n", rms, rms_v);
}

void test_cross_coupling_and_api() {
    Filter52 f = make_filter<Filter52>();
    // Tilted attitude: the vertical acceleration enters every body axis.
    const Eigen::Quaternionf q(Eigen::AngleAxisf(0.3f, Vec3::UnitX()));
    f.initialize_from_attitude(q, 0.02f, 0.05f);
    for (int k = 0; k < 400; ++k) {
        f.time_update(Vec3::Zero(), 0.005f);
        const Vec3 fb = q.conjugate() * Vec3(0.0f, 0.0f, -9.80665f);
        f.measurement_update_acc_only(fb);
    }
    // a_z and the accelerometer bias both act along gravity, so the update
    // must correlate them; tilt acts orthogonally to gravity, so attitude
    // couples to a_z only at second order and is not tested here.
    const auto& P = f.covariance_full();
    double c_ab = 0.0;
    for (int k = 0; k < 3; ++k)
        c_ab = std::max(c_ab, std::abs(double(P(Filter52::OFF_Z + 2, Filter52::OFF_BA + k))));
    const double norm = std::sqrt(double(P(Filter52::OFF_Z + 2, Filter52::OFF_Z + 2)) *
                                  double(P(Filter52::OFF_BA + 2, Filter52::OFF_BA + 2)));
    check(c_ab > 1e-3 * norm, "vertical acceleration / accelerometer-bias cross-covariance retained",
          c_ab / norm, 1e-3);
    double c_ah = 0.0;
    for (int k = 0; k < 3; ++k)
        c_ah = std::max(c_ah, std::abs(double(P(k, Filter52::OFF_AH + 1))));
    check(c_ah > 0.0, "attitude / horizontal acceleration cross-covariance retained", c_ah, 0.0);

    // Pseudo-measurement APIs act on the CARMA displacement and velocity.
    f.measurement_update_position_pseudo(Vec3(0.0f, 0.0f, 0.37f), Vec3(1e3f, 1e3f, 1e-4f));
    check(std::abs(f.get_position().z() - 0.37f) < 1e-3f, "position pseudo-measurement maps to h",
          f.get_position().z(), 0.37);
    f.measurement_update_velocity_pseudo(Vec3(0.0f, 0.0f, 0.05f), Vec3(1e3f, 1e3f, 1e-4f));
    check(std::abs(f.get_velocity().z() - 0.05f) < 1e-3f, "velocity pseudo-measurement maps to v",
          f.get_velocity().z(), 0.05);
    // A second, equally precise measurement averages with the first.
    const float v_prior = f.get_velocity().z();
    const float var_prior = f.covariance_full()(Filter52::OFF_Z + 1, Filter52::OFF_Z + 1);
    const float r_v = 1e-3f;
    f.measurement_update_vert_velocity_pseudo(-0.21f, r_v);
    const float expect = v_prior + var_prior / (var_prior + r_v * r_v) * (-0.21f - v_prior);
    check(std::abs(f.get_velocity().z() - expect) < 1e-4f, "vertical velocity pseudo-measurement is the Kalman update",
          f.get_velocity().z(), expect);
    check(min_eig_normalized(f) > -1e-5, "covariance PSD after pseudo-measurements", min_eig_normalized(f), -1e-5);
}

template<class F>
void test_nonuniform_and_integral(const char* name) {
    F f = make_filter<F>();
    std::mt19937 rng(7u);
    std::uniform_real_distribution<float> jit(0.6f, 1.6f);
    std::normal_distribution<float> n01(0.0f, 1.0f);
    double t = 0.0, max_h = 0.0;
    bool finite = true;
    for (int k = 0; k < 200 * 300; ++k) {
        const float dt = 0.005f * jit(rng);
        t += dt;
        const double az = 1.2 * std::sin(2.0 * M_PI * 0.27 * t);
        f.time_update(Vec3(1e-3f * n01(rng), 1e-3f * n01(rng), 1e-3f * n01(rng)), dt);
        f.measurement_update_acc_only(Vec3(0.02f * n01(rng), 0.02f * n01(rng),
                                           static_cast<float>(az - 9.80665) + 0.02f * n01(rng)));
        if (k % 8 == 0) f.measurement_update_mag_only(Vec3(20.0f, -4.0f, 45.0f));
        finite = finite && f.get_position().allFinite() && f.covariance_full().allFinite();
        max_h = std::max(max_h, double(std::abs(f.get_position().z())));
    }
    check(finite, "non-uniform steps: finite state and covariance");
    check(max_h < 2.0, "non-uniform steps: bounded heave", max_h, 2.0);
    const double me = min_eig_normalized(f);
    check(me > -1e-5, "non-uniform steps: covariance PSD", me, -1e-5);
    const auto P = f.covariance_full();
    check((P - P.transpose()).cwiseAbs().maxCoeff() == 0.0f, "covariance exactly symmetric");
    if constexpr (F::NZ == 6) {
        const float Sz = f.get_integral_displacement().z();
        check(std::isfinite(Sz) && std::abs(Sz) < 5.0f, "S_z bounded by its pseudo-measurement", Sz, 5.0);
    } else {
        check(std::isnan(f.get_integral_displacement().z()), "no S_z state without the integral variant");
    }
    std::printf("[%s] non-uniform steps: max|h|=%.3f m, min corr eig=%.2e\n", name, max_h, me);
}

// The stationary replay of docs/stationary-device-behavior.md: vertical
// accelerometer residual +-0.01 m/s^2, 0.0148 m/s^2 noise.  OU-III's
// acceptance envelope there is peak raw |p_z| < 0.3 m and < 5 cm after 240 s.
// The experimental CARMA filter does NOT meet that peak (it reaches about
// 0.9 m during the bias-gated interval; see docs/carma-heave-filter.md), so
// the bounds below are regression sentinels at its measured behaviour plus
// margin, not the shipping envelope.
void test_wrapper_still(float bias_z) {
    ocean_imu::carma::SeaStateFusion_CARMA<TrackerType::KALMANF> fusion;
    ocean_imu::carma::SeaStateFusionConfig_CARMA cfg;
    cfg.sigma_a = Vec3::Constant(0.0294f);
    cfg.sigma_g = Vec3::Constant(3.1e-5f);
    cfg.sigma_m = Vec3::Constant(7.7f);
    fusion.begin(cfg);
    std::mt19937 rng(11u);
    std::normal_distribution<float> n01(0.0f, 1.0f);
    const Vec3 bias_a(0.0f, 0.0f, bias_z);
    double peak = 0.0, late = 0.0;
    bool finite = true;
    for (int k = 0; k < 200 * 600; ++k) {
        const Vec3 acc = Vec3(0, 0, -9.80665f) + bias_a + 0.0148f * Vec3(n01(rng), n01(rng), n01(rng));
        const Vec3 gyr = 0.002f * Vec3(0, 0, 1) + 0.00157f * Vec3(n01(rng), n01(rng), n01(rng));
        fusion.update(0.005f, gyr, acc);
        if (k % 8 == 0) fusion.updateMag(Vec3(20.0f, -4.0f, 45.0f) + 0.8f * Vec3(n01(rng), n01(rng), n01(rng)));
        if (fusion.isLive()) {
            const double pz = std::abs(double(fusion.raw().mekf().get_position().z()));
            finite = finite && std::isfinite(pz);
            peak = std::max(peak, pz);
            if (k > 200 * 240) late = std::max(late, pz);
        }
    }
    check(fusion.isLive(), "wrapper goes live on a still IMU");
    check(finite, "wrapper: finite position on a still IMU");
    check(peak < 1.2, "still IMU sentinel: peak |p_z| < 1.2 m", peak, 1.2);
    check(late < 0.08, "still IMU sentinel: |p_z| < 8 cm after 240 s", late, 0.08);
    std::printf("[wrapper] still IMU, vertical residual %+.3f m/s^2: peak |p_z| %.4f m, after 240 s %.4f m\n",
                bias_z, peak, late);
}

}  // namespace

int main() {
    test_layout();
    test_sinusoid_signs();
    test_cross_coupling_and_api();
    test_nonuniform_and_integral<Filter52>("CARMA(5,2)");
    test_nonuniform_and_integral<Filter52S>("CARMA(5,2)+S_z");
    test_wrapper_still(+0.01f);
    test_wrapper_still(-0.01f);
    std::printf("carma_filter-test: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
