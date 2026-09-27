// Copyright (c) 2026 Mikhail Grushinskiy
// Literal construction/H18/refinement/A21 regression for the analytic
// rest/rocking/rest ambiguity in docs/ou3-regime-design.md. Finite replay
// checks source behavior; the rational/analytic certificate supplies bounds.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <numbers>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

static int fail(const char* message) {
    std::cerr << "FAIL: " << message << '\n';
    return 1;
}

int main() {
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(0.12f);
    cfg.sigma_g.setConstant(0.00135f);
    cfg.sigma_m.setConstant(0.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    Fusion f;
    f.begin(cfg);
    constexpr double alpha = .001, nu = .025, start = 300, gravity = 9.80665;
    // Keep physical gravity separate from the literal float32 source value.
    const Eigen::Vector3d representation_noise(0, 0, gravity - double(g_std));
    if (representation_noise.norm() > .3)
        return fail("gravity representation residual exceeds the unchanged noise bound");
    constexpr double end = start + 2 * std::numbers::pi / nu;
    double max_tilt = 0, previous_bg = 0;
    Eigen::Vector3d previous_ba = Eigen::Vector3d::Zero();
    int live = -1, refined = -1, active = -1, accepted = 0;
    for (int k = 0; k <= 180000; ++k) {
        const double t = double(k) / 200;
        const bool moving = start < t && t < end;
        const double phase = nu * (t - start);
        const double sn = std::sin(phase), cs = std::cos(phase);
        const double phi = moving ? alpha * sn * sn * sn : 0;
        const double rate = moving ? 3 * alpha * nu * sn * sn * cs : 0;
        const Eigen::Vector3d ba(0, gravity * std::sin(phi),
                                gravity * (std::cos(phi) - 1));
        const double bg = -rate;
        if (ba.norm() > .22516660498395405 || std::abs(bg) > .02)
            return fail("hidden physical bias exceeded the unchanged envelope");
        if (k && ((ba - previous_ba).norm() > .001 * .005 + 1e-14 ||
                  std::abs(bg - previous_bg) > .00001 * .005 + 1e-14))
            return fail("physical bias continuity/rate violated at a join");
        previous_ba = ba;
        previous_bg = bg;
        max_tilt = std::max(max_tilt, std::abs(phi));
        // True R'=Rx(-phi); the compensating bias makes the sensor equation
        // exactly f=-g ez. Evaluate the identity separately, then pass the
        // canonical common packets rather than rounded cancellation noise.
        const Eigen::Vector3d physical_force(0, -gravity * std::sin(phi),
                                            -gravity * std::cos(phi));
        if ((physical_force + ba + representation_noise - Eigen::Vector3d(0, 0, -double(g_std))).norm() > 1e-12 ||
            rate + bg != 0)
            return fail("physical rest and rocking do not share the measured record");
        f.update(.005f, Eigen::Vector3f::Zero(), Eigen::Vector3f(0, 0, -g_std), 35);
        if (k % 8 == 0) {
            f.updateMag(Eigen::Vector3f(75, 0, 0));
            if (f.raw().mekf().lastMagDiag().accepted) ++accepted;
            if (active >= 0 && !f.raw().mekf().lastMagDiag().accepted)
                return fail("an A21 magnetic correction was not applied");
        }
        if (live < 0 && f.isLive()) live = k;
        if (refined < 0 && f.hasRefinedMagReference()) refined = k;
        if (active < 0 && f.raw().mekf().acc_bias_updates_enabled()) active = k;
        const auto& core = f.raw().mekf();
        if (!core.covariance_full().allFinite() || !core.quaternion_boat().coeffs().allFinite())
            return fail("long stationary-looking replay became nonfinite");
        if (core.quaternion_boat().vec().norm() > 1e-6f || core.get_position().norm() > 1e-6f ||
            core.get_velocity().norm() > 1e-6f || core.get_world_accel().norm() > 1e-6f ||
            core.get_acc_bias().norm() > 1e-6f || core.gyroscope_bias().norm() > 1e-6f)
            return fail("unchanged common packet record no longer has the quiet nominal solution");
        if (k % 200 == 0 && core.covariance_full().llt().info() != Eigen::Success)
            return fail("recorded covariance lost positive definiteness");
    }
    if (live < 0 || refined < live || active < refined || accepted < 1000 || max_tilt < .00099)
        return fail("construction/release or genuine hidden motion was not exercised");
    std::cout << "OU-III regime ambiguity PASS: 900 s continuous construction/rest/rocking/rest; live="
              << live << " refined=" << refined << " active=" << active
              << " applied_mag=" << accepted << " hidden_tilt=" << max_tilt
              << "; finite replay, no universal detector/stability promotion\n";
}
