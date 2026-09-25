#pragma once
// Local-gravity convention through the real device path: a still sensor at a
// site with gravity g_local is calibrated against g_local (as the wizard
// does), its calibrated SI acceleration is fed to the fusion filter, and the
// filter removes the gravity it is configured with. With g_local the
// translational acceleration of an ideal still device is zero; with the
// standard 9.80665 it is g_std - g_local. Checked on the filter's own states:
// the accelerometer model f = R^T (a_w - g) + b(T) gives
// (a_w + R b(T))_z = (R f)_z + g, which gravity alone enters.
// The OU wrappers' private Mahony observer (getAccelVertical) settles with a
// small tilt offset that does not depend on gravity (about -0.033 m/s^2 at
// rest); only its change with the configured gravity is checked.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include "../common/GravityChain.h"
#if defined(LOCAL_GRAVITY_TFG)
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
#elif defined(LOCAL_GRAVITY_OU2)
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#else
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#endif
extern const float g_std = 9.80665f;  // generic default of the filter headers
namespace local_gravity {
using V = Eigen::Vector3f;
#if defined(LOCAL_GRAVITY_TFG)
using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
constexpr const char* name = "TFG";
#elif defined(LOCAL_GRAVITY_OU2)
using Fusion = SeaStateFusion_OU_II<TrackerType::KALMANF>;
constexpr const char* name = "OU-II";
#else
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
constexpr const char* name = "OU-III";
#endif
constexpr float dt = 0.005f;

auto& raw(Fusion& f) {
#if defined(LOCAL_GRAVITY_TFG)
    return f;
#else
    return f.raw();
#endif
}
Eigen::Quaternionf attitude(Fusion& f) {
#if defined(LOCAL_GRAVITY_TFG)
    return f.mekf().quaternion();
#else
    return f.raw().mekf().quaternion_boat();
#endif
}

struct Result {
    bool finite = true, live = false;
    double e_z = 0, a_w_z = 0, vertical = NAN;  // tail means [m/s^2]
    int n = 0;
};

// Still device, attitude roll 0.15 / pitch 0.1 / yaw 0.5 rad, calibrated
// accelerometer `cal`, filter gravity `g_filter`.
Result replay(const gravity_chain::Calibration& cal, const gravity_chain::Sensor& sensor,
              double g_site, float g_filter) {
    Fusion f;
    Fusion::Config cfg;
    cfg.sigma_a = V::Constant(0.12f);
    cfg.sigma_m = V::Constant(0.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    cfg.gravity_magnitude = g_filter;
#if defined(LOCAL_GRAVITY_TFG)
    cfg.gyro_noise_density = 0.00135f;
#else
    cfg.sigma_g = V::Constant(0.00135f);
    cfg.enable_displacement_detrend = false;
#endif
    f.begin(cfg);
    raw(f).enableTuner(true);
    raw(f).setAccNoiseFloorSigma(0.12f);

    const Eigen::Quaterniond q = Eigen::AngleAxisd(0.5, Eigen::Vector3d::UnitZ()) *
        Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(0.15, Eigen::Vector3d::UnitX());
    // Specific force of a still body (NED): R^T (0 - (0,0,g_site)).
    const Eigen::Vector3d f_body = q.conjugate() * Eigen::Vector3d(0, 0, -g_site);
    const V a_raw = gravity_chain::readingSI(sensor, f_body);
    const V a_cal = cal.rt.applyAccel(a_raw, 25.0f);
    const V mag = (q.conjugate() * Eigen::Vector3d(20, 0, 43)).cast<float>();

    Result r;
    const int count = 48000;  // 240 s
    for (int k = 0; k < count; ++k) {
        f.update(dt, V::Zero(), a_cal, 25.0f);
        if (k % 8 == 0) f.updateMag(mag);
        if (k < count - 4000) continue;  // tail: last 20 s
        const Eigen::Matrix3f R = attitude(f).toRotationMatrix();
        const V a_w = raw(f).mekf().get_world_accel();
        const V b = raw(f).mekf().get_acc_bias_at_temperature(25.0f);
        const V e = a_w + R * b;
        r.finite = r.finite && e.allFinite();
        r.e_z += e.z(); r.a_w_z += a_w.z();
#if !defined(LOCAL_GRAVITY_TFG)
        r.vertical = (r.n == 0 ? 0.0 : r.vertical) + raw(f).getAccelVertical();
#endif
        ++r.n;
    }
    r.live = f.isLive();
    r.e_z /= r.n; r.a_w_z /= r.n;
    if (std::isfinite(r.vertical)) r.vertical /= r.n;
    return r;
}

int run() {
    int failures = 0;
    auto check = [&](bool ok, const char* text) {
        if (!ok) { std::cerr << "FAIL: " << name << ' ' << text << '\n'; ++failures; }
    };
    using atoms3r_ical::ImuCalCfg;
    const double g_local = ImuCalCfg::g_cal_local;
    const auto sensor = gravity_chain::Sensor::typical();
    const auto cal = gravity_chain::calibrate(sensor, g_local, g_local);
    check(cal.ok && cal.rt.acc.ok, "calibration at the local gravity was not accepted");

    const Result good = replay(cal, sensor, g_local, ImuCalCfg::g_cal_local);
    const Result wrong = replay(cal, sensor, g_local, ImuCalCfg::g_std);
    const double expected_wrong = double(ImuCalCfg::g_std) - g_local;  // g_filter - g_site
    std::cout << name << " g_local=" << g_local << " filter g_local: e_z=" << good.e_z
              << " a_w_z=" << good.a_w_z << " vertical=" << good.vertical
              << " | filter g_std: e_z=" << wrong.e_z << " a_w_z=" << wrong.a_w_z
              << " vertical=" << wrong.vertical << " (expected e_z " << expected_wrong << ")\n";
    check(good.finite && wrong.finite && good.live && wrong.live, "replay not finite or not Live");
    check(std::abs(good.e_z) < 2e-4, "local gravity: still device has translational acceleration");
    check(std::abs(good.a_w_z) < 1e-3, "local gravity: translational acceleration state not ~0");
    check(std::abs(wrong.e_z - expected_wrong) < 2e-4, "g_std filter: residual is not g_std - g_local");
#if !defined(LOCAL_GRAVITY_TFG)
    // Measurement-only vertical channel (no filter state can absorb it).
    check(std::abs((good.vertical - wrong.vertical) - expected_wrong) < 2e-5,
          "vertical channel does not remove the configured gravity");
#endif
    return failures ? 1 : 0;
}
} // namespace local_gravity
