// Self-calibrating IMU lever arm in Kalman3D_Wave_OU_III.
//
// The lever-arm study answers "what does an off-CG installation cost, and how
// much does modelling it recover".  Both of its modelling arms are handed the
// lever vector r: the exact arm gets the record's angular kinematics too, the
// gyro-derived arm reconstructs those from the measured rate, but neither
// derives r itself.  This file exercises the third option -- r as three more
// filter states -- and pins down both what it recovers and where it cannot.
//
// The physics is one line:
//
//     a_IMU = a_track + alpha x r + omega x (omega x r) = a_track + M(w,a) r,
//
// linear in r, so d(a_IMU)/dr = M is exact and the augmentation costs the
// filter no extra linearization error.  What it does cost is observability:
// M annihilates any r parallel to the instantaneous rotation axis, so a
// single-axis rotation can never reveal the component of r along it.  Both
// halves of that are tested here.

#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

namespace {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

using PlainFilter = Kalman3D_Wave_OU_III<double, true, true, false>;
using CalibFilter = Kalman3D_Wave_OU_III<double, true, true, true>;

constexpr double kGravity = 9.80665;
constexpr double kDt = 0.005;   // 200 Hz, the rate the deployed filter runs at

int failures = 0;

bool check(bool ok, const std::string& what) {
    std::cout << (ok ? "[ ok ] " : "[FAIL] ") << what << "\n";
    if (!ok) ++failures;
    return ok;
}

// ---------------------------------------------------------------------------
// A rigid body that rotates about the point the filter tracks.  The tracked
// point does not translate, so every bit of specific force the IMU sees beyond
// gravity is the lever-arm term -- which is exactly the signal the estimator
// has to find r in.
// ---------------------------------------------------------------------------
struct RotatingBody {
    Eigen::Quaterniond q_wb = Eigen::Quaterniond::Identity();  // world -> body

    // Body rate and its derivative, both analytic so the truth carries no
    // differentiation error of its own.
    static Vec3 omega(double t, bool single_axis) {
        if (single_axis) {
            // One fixed axis, deliberately not aligned with a body axis.
            const Vec3 axis = Vec3(1.0, 2.0, 3.0).normalized();
            return axis * (0.45 * std::sin(2.0 * M_PI * 0.25 * t));
        }
        return Vec3(0.40 * std::sin(2.0 * M_PI * 0.23 * t),
                    0.30 * std::sin(2.0 * M_PI * 0.17 * t + 0.7),
                    0.22 * std::sin(2.0 * M_PI * 0.31 * t + 1.9));
    }

    static Vec3 alpha(double t, bool single_axis) {
        if (single_axis) {
            const Vec3 axis = Vec3(1.0, 2.0, 3.0).normalized();
            return axis * (0.45 * 2.0 * M_PI * 0.25 * std::cos(2.0 * M_PI * 0.25 * t));
        }
        return Vec3(0.40 * 2.0 * M_PI * 0.23 * std::cos(2.0 * M_PI * 0.23 * t),
                    0.30 * 2.0 * M_PI * 0.17 * std::cos(2.0 * M_PI * 0.17 * t + 0.7),
                    0.22 * 2.0 * M_PI * 0.31 * std::cos(2.0 * M_PI * 0.31 * t + 1.9));
    }

    void advance(const Vec3& w, double dt) {
        // Same convention the filter propagates with: q_wb(k+1) = dq(-w dt) q_wb(k)
        const Eigen::Quaterniond dq =
            ocean_imu::kalman::ou_detail::quat_from_delta_theta<double>((-w * dt).eval());
        q_wb = dq * q_wb;
        q_wb.normalize();
    }

    // Specific force at an IMU sitting at r, body frame, NED.
    Vec3 specific_force(const Vec3& w, const Vec3& a, const Vec3& r) const {
        const Vec3 g_world(0.0, 0.0, kGravity);
        return q_wb * (-g_world) + a.cross(r) + w.cross(w.cross(r));
    }
};

// ---------------------------------------------------------------------------
// 1. The Jacobian is the model, not an approximation of it.
// ---------------------------------------------------------------------------
bool test_sensitivity_is_the_exact_derivative() {
    const Vec3 w(0.31, -0.22, 0.17);
    const Vec3 a(-0.9, 0.4, 1.3);
    const Vec3 r(0.23, -0.11, 0.37);

    const Mat3 M = CalibFilter::lever_arm_sensitivity(w, a);

    // Linearity: M r must reproduce the term itself, with no offset.
    const Vec3 direct = CalibFilter::lever_arm_acceleration(w, a, r);
    const bool linear = (M * r - direct).cwiseAbs().maxCoeff() < 1e-12 &&
                        CalibFilter::lever_arm_acceleration(w, a, Vec3::Zero())
                            .cwiseAbs().maxCoeff() < 1e-15;

    // And a central difference, which must agree to round-off because the
    // function is linear.
    Mat3 numeric;
    const double h = 1e-6;
    for (int j = 0; j < 3; ++j) {
        Vec3 rp = r, rm = r;
        rp(j) += h;
        rm(j) -= h;
        numeric.col(j) = (CalibFilter::lever_arm_acceleration(w, a, rp) -
                          CalibFilter::lever_arm_acceleration(w, a, rm)) / (2.0 * h);
    }
    const double err = (numeric - M).cwiseAbs().maxCoeff();

    return check(linear && err < 1e-8,
                 "lever-arm sensitivity equals the exact derivative (max err " +
                     std::to_string(err) + ")");
}

// The unobservable direction, stated as algebra before it is measured as a
// filter behaviour: M(w, a) n = 0 whenever w, a and n are collinear.
bool test_rotation_axis_is_in_the_null_space() {
    const Vec3 axis = Vec3(1.0, 2.0, 3.0).normalized();
    const Mat3 M = CalibFilter::lever_arm_sensitivity(axis * 0.4, axis * 1.7);
    const double along = (M * axis).norm();
    const Vec3 perp = axis.cross(Vec3::UnitX()).normalized();
    const double across = (M * perp).norm();
    return check(along < 1e-12 && across > 0.1,
                 "a lever arm along the rotation axis is invisible (|M n| = " +
                     std::to_string(along) + ", |M perp| = " +
                     std::to_string(across) + ")");
}

// ---------------------------------------------------------------------------
// 2. Compiling the states out really does cost nothing.
// ---------------------------------------------------------------------------
bool test_default_instantiation_is_unchanged() {
    // 6 base (attitude error + gyro bias) + 12 linear + 3 accelerometer bias.
    const bool dims = PlainFilter::state_dimension() == 21 &&
                      CalibFilter::state_dimension() == 24 &&
                      !PlainFilter::has_lever_arm_states() &&
                      CalibFilter::has_lever_arm_states();

    // The scalars the calibration needs live in a [[no_unique_address]] member
    // that is empty without the states, so the plain filter is not merely
    // "close to" its old layout -- it is byte-for-byte the same object.
    const bool layout = sizeof(PlainFilter) < sizeof(CalibFilter);

    const Vec3 sa(0.02, 0.02, 0.02), sg(1e-3, 1e-3, 1e-3), sm(0.5, 0.5, 0.5);
    PlainFilter plain(sa, sg, sm);
    plain.set_imu_lever_arm_body(Vec3(0.1, -0.2, 0.3));
    const bool passthrough =
        (plain.get_imu_lever_arm_body() - Vec3(0.1, -0.2, 0.3)).norm() < 1e-15 &&
        plain.get_imu_lever_arm_covariance().norm() == 0.0 &&
        !plain.imu_lever_arm_estimation_enabled();

    return check(dims && layout && passthrough,
                 "with_lever_arm=false keeps the 21-state filter and treats r as an input");
}

// ---------------------------------------------------------------------------
// 3. The estimator finds a lever arm it was never told.
// ---------------------------------------------------------------------------
struct RunResult {
    Vec3 r_hat = Vec3::Zero();
    Mat3 P_rr = Mat3::Zero();
    double tilt_rms_deg = 0.0;
};

RunResult run(const Vec3& r_true,
              bool estimate,
              bool single_axis,
              double seconds,
              const Vec3& r_model = Vec3::Zero())
{
    const Vec3 sa(0.02, 0.02, 0.02), sg(2e-4, 2e-4, 2e-4), sm(0.5, 0.5, 0.5);
    CalibFilter f(sa, sg, sm);

    RotatingBody body;
    f.set_linear_block_enabled(false);   // attitude-only: isolate the lever arm
    f.initialize_from_attitude(body.q_wb.conjugate(), 1e-3, 1e-3);
    f.set_Racc_std(Vec3::Constant(0.05));

    if (estimate) {
        f.enable_imu_lever_arm_estimation(Vec3::Zero(), 0.5);
    } else if (r_model.squaredNorm() > 0.0) {
        f.set_imu_lever_arm_body(r_model);
    }

    const int steps = static_cast<int>(seconds / kDt);
    double tilt_sq = 0.0;
    int scored = 0;

    for (int k = 0; k < steps; ++k) {
        const double t = k * kDt;
        const Vec3 w = RotatingBody::omega(t, single_axis);
        const Vec3 a = RotatingBody::alpha(t, single_axis);

        const Vec3 acc = body.specific_force(w, a, r_true);
        f.time_update(w, kDt);
        f.measurement_update_acc_only(acc);
        body.advance(w, kDt);

        if (t > seconds * 0.5) {
            const Eigen::Quaterniond q_est = f.quaternion_boat();       // body -> world
            const Eigen::Quaterniond q_err = q_est.conjugate() * body.q_wb.conjugate();
            const double ang = 2.0 * std::asin(std::min(1.0, q_err.vec().norm()));
            tilt_sq += ang * ang;
            ++scored;
        }
    }

    RunResult out;
    out.r_hat = f.get_imu_lever_arm_body();
    out.P_rr = f.get_imu_lever_arm_covariance();
    out.tilt_rms_deg = std::sqrt(tilt_sq / std::max(1, scored)) * 180.0 / M_PI;
    return out;
}

bool test_converges_to_the_true_lever_arm() {
    const Vec3 r_true(0.30, -0.20, 0.15);
    const RunResult est = run(r_true, /*estimate=*/true, /*single_axis=*/false, 600.0);

    const double err = (est.r_hat - r_true).norm();
    const double sigma = std::sqrt(est.P_rr.trace() / 3.0);

    std::printf("       r_true = [%.3f %.3f %.3f]  r_hat = [%.3f %.3f %.3f]\n",
                r_true.x(), r_true.y(), r_true.z(),
                est.r_hat.x(), est.r_hat.y(), est.r_hat.z());
    std::printf("       |error| = %.4f m, reported sigma = %.4f m\n", err, sigma);

    // 3 cm on a 39 cm arm, from a zero prior, is the recovery this claims.
    // The reported covariance has to shrink with it: an estimate that is right
    // but still calls itself a 0.5 m prior has not actually learned anything.
    return check(err < 0.03 && sigma < 0.05,
                 "a 39 cm lever arm is recovered from a zero prior");
}

bool test_estimated_arm_matches_a_surveyed_one() {
    const Vec3 r_true(0.30, -0.20, 0.15);
    const RunResult unmodeled = run(r_true, false, false, 600.0);
    const RunResult surveyed  = run(r_true, false, false, 600.0, r_true);
    const RunResult estimated = run(r_true, true,  false, 600.0);

    std::printf("       tilt RMS: unmodeled %.4f deg, surveyed %.4f deg, estimated %.4f deg\n",
                unmodeled.tilt_rms_deg, surveyed.tilt_rms_deg, estimated.tilt_rms_deg);

    // The surveyed arm is the bound: it is handed the answer.  The estimated
    // arm has to reach it rather than merely beat doing nothing.
    const bool helps = estimated.tilt_rms_deg < 0.5 * unmodeled.tilt_rms_deg;
    const bool reaches_bound =
        estimated.tilt_rms_deg < surveyed.tilt_rms_deg + 0.02;

    return check(helps && reaches_bound,
                 "estimating r recovers what surveying it would have");
}

bool test_single_axis_rotation_leaves_the_axis_unobservable() {
    const Vec3 axis = Vec3(1.0, 2.0, 3.0).normalized();
    const Vec3 r_true = axis * 0.25 + axis.cross(Vec3::UnitX()).normalized() * 0.25;

    const RunResult est = run(r_true, true, /*single_axis=*/true, 600.0);

    // Split both the error and the reported uncertainty along and across the
    // one rotation axis the body ever turned about.
    const double err_along = std::abs(axis.dot(est.r_hat - r_true));
    const Vec3 perp = axis.cross(Vec3::UnitX()).normalized();
    const double err_across = std::abs(perp.dot(est.r_hat - r_true));

    const double var_along = axis.transpose() * est.P_rr * axis;
    const double var_across = perp.transpose() * est.P_rr * perp;

    std::printf("       along axis: err %.4f m, sigma %.4f m | across: err %.4f m, sigma %.4f m\n",
                err_along, std::sqrt(var_along), err_across, std::sqrt(var_across));

    // The component across the axis is recovered; the one along it is not, and
    // the filter says so -- it keeps essentially its whole prior width there
    // instead of reporting false confidence.
    const bool across_learned = err_across < 0.03 && var_across < var_along;
    const bool along_admitted = std::sqrt(var_along) > 0.4;

    return check(across_learned && along_admitted,
                 "a single rotation axis hides the arm along it, and the covariance admits it");
}

// The harder case, and the deployed one: the linear block is running, so the
// latent OU world acceleration a_w is free to compete for the same
// accelerometer channel.  a_w can imitate the lever-arm term over a short
// window -- what separates them is that the lever term is rigidly tied to a
// rate the filter also measures, while a_w's prior says it is unpredictable.
// The separation is real but it is not free: convergence is slower and the
// estimate breathes, so the tolerance here is centimetres rather than the
// millimetres the attitude-only case reaches.
bool test_converges_with_the_linear_block_running() {
    const Vec3 sa(0.02, 0.02, 0.02), sg(2e-4, 2e-4, 2e-4), sm(0.5, 0.5, 0.5);
    CalibFilter f(sa, sg, sm);

    RotatingBody body;
    f.set_linear_block_enabled(true);
    f.initialize_from_attitude(body.q_wb.conjugate(), 1e-3, 1e-3);
    f.set_aw_time_constant(2.1);
    f.set_aw_stationary_std(Vec3::Constant(2.2));
    f.enable_imu_lever_arm_estimation(Vec3::Zero(), 0.5);

    const Vec3 r_true(0.30, -0.20, 0.15);
    const Vec3 g_world(0.0, 0.0, kGravity);

    // A wave-like world acceleration of the tracked point, so a_w has real
    // work of its own to do rather than sitting at zero.
    auto world_accel = [](double t) {
        return Vec3(0.6 * std::sin(2.0 * M_PI * 0.19 * t + 0.3),
                    0.5 * std::sin(2.0 * M_PI * 0.21 * t + 1.1),
                    1.2 * std::sin(2.0 * M_PI * 0.20 * t));
    };

    const int steps = static_cast<int>(1200.0 / kDt);
    for (int k = 0; k < steps; ++k) {
        const double t = k * kDt;
        const Vec3 w = RotatingBody::omega(t, false);
        const Vec3 a = RotatingBody::alpha(t, false);
        const Vec3 acc = body.q_wb * (world_accel(t) - g_world)
                       + CalibFilter::lever_arm_acceleration(w, a, r_true);
        f.time_update(w, kDt);
        f.measurement_update_acc_only(acc);
        body.advance(w, kDt);
    }

    const Vec3 r_hat = f.get_imu_lever_arm_body();
    const double err = (r_hat - r_true).norm();
    const double sigma = std::sqrt(f.get_imu_lever_arm_covariance().trace() / 3.0);
    std::printf("       with a_w free: r_hat = [%.3f %.3f %.3f], |error| = %.4f m, sigma = %.4f m\n",
                r_hat.x(), r_hat.y(), r_hat.z(), err, sigma);

    return check(err < 0.05 && sigma < 0.05,
                 "r is still identifiable against a free latent OU acceleration");
}

// ---------------------------------------------------------------------------
// 4. Holding the calibration really holds it.
// ---------------------------------------------------------------------------
bool test_frozen_updates_do_not_move_the_estimate() {
    const Vec3 sa(0.02, 0.02, 0.02), sg(2e-4, 2e-4, 2e-4), sm(0.5, 0.5, 0.5);
    CalibFilter f(sa, sg, sm);

    RotatingBody body;
    f.set_linear_block_enabled(false);
    f.initialize_from_attitude(body.q_wb.conjugate(), 1e-3, 1e-3);
    // A magnetometer channel as well, so the hold is tested against an update
    // that reaches r only through cross-covariance rather than directly.
    const Vec3 mag_world(20.0, 0.0, 45.0);
    f.set_mag_world_ref(mag_world);

    const Vec3 r_seed(0.11, -0.05, 0.22);
    f.enable_imu_lever_arm_estimation(r_seed, 0.5);
    f.set_lever_arm_updates_frozen(true);

    const Vec3 r_true(0.30, -0.20, 0.15);
    for (int k = 0; k < 20000; ++k) {
        const double t = k * kDt;
        const Vec3 w = RotatingBody::omega(t, false);
        const Vec3 a = RotatingBody::alpha(t, false);
        f.time_update(w, kDt);
        f.measurement_update_acc_only(body.specific_force(w, a, r_true));
        f.measurement_update_mag_only(body.q_wb * mag_world);
        body.advance(w, kDt);
    }

    const double moved = (f.get_imu_lever_arm_body() - r_seed).cwiseAbs().maxCoeff();
    const bool held = moved == 0.0;

    // And releasing it lets the same filter learn again.
    f.set_lever_arm_updates_frozen(false);
    for (int k = 0; k < 120000; ++k) {
        const double t = k * kDt;
        const Vec3 w = RotatingBody::omega(t, false);
        const Vec3 a = RotatingBody::alpha(t, false);
        f.time_update(w, kDt);
        f.measurement_update_acc_only(body.specific_force(w, a, r_true));
        body.advance(w, kDt);
    }
    const double err_after = (f.get_imu_lever_arm_body() - r_true).norm();

    std::printf("       frozen drift %.3g m, error after release %.4f m\n", moved, err_after);
    return check(held && err_after < 0.05,
                 "a frozen calibration does not move, and resumes when released");
}

bool test_warmup_holds_rather_than_discards() {
    const Vec3 sa(0.02, 0.02, 0.02), sg(2e-4, 2e-4, 2e-4), sm(0.5, 0.5, 0.5);
    CalibFilter f(sa, sg, sm);
    const Vec3 r_seed(0.11, -0.05, 0.22);
    f.enable_imu_lever_arm_estimation(r_seed, 0.5);

    f.set_warmup_mode(true);
    const bool held = !f.imu_lever_arm_updates_active() &&
                      f.imu_lever_arm_estimation_enabled() &&
                      (f.get_imu_lever_arm_body() - r_seed).norm() == 0.0;

    f.set_warmup_mode(false);
    const bool resumed = f.imu_lever_arm_updates_active() &&
                         (f.get_imu_lever_arm_body() - r_seed).norm() == 0.0;

    return check(held && resumed,
                 "warmup holds the calibration instead of throwing it away");
}

bool test_projection_bounds_a_runaway_estimate() {
    const Vec3 sa(0.02, 0.02, 0.02), sg(2e-4, 2e-4, 2e-4), sm(0.5, 0.5, 0.5);
    CalibFilter f(sa, sg, sm);
    RotatingBody body;
    f.set_linear_block_enabled(false);
    f.initialize_from_attitude(body.q_wb.conjugate(), 1e-3, 1e-3);
    f.enable_imu_lever_arm_estimation(Vec3::Zero(), 50.0);  // absurd prior
    f.set_lever_arm_limit(1.0);

    // Feed a specific force that no lever arm can explain, and check the
    // estimate stays inside the ball rather than running off to absorb it.
    for (int k = 0; k < 40000; ++k) {
        const double t = k * kDt;
        const Vec3 w = RotatingBody::omega(t, false);
        const Vec3 a = RotatingBody::alpha(t, false);
        Vec3 acc = body.specific_force(w, a, Vec3(3.0, -4.0, 5.0));
        acc += Vec3(0.5 * std::sin(t), -0.4 * std::cos(0.7 * t), 0.3);
        f.time_update(w, kDt);
        f.measurement_update_acc_only(acc);
        body.advance(w, kDt);
    }

    const double n = f.get_imu_lever_arm_body().norm();
    const bool bounded = n <= 1.0 + 1e-9 && f.lever_arm_mirror_consistent();
    std::printf("       |r_hat| = %.4f m against a 1.0 m limit\n", n);
    return check(bounded, "the estimate stays inside the ball it is confined to");
}

} // namespace

int main() {
    std::cout << "OU-III self-calibrating IMU lever arm\n";

    test_sensitivity_is_the_exact_derivative();
    test_rotation_axis_is_in_the_null_space();
    test_default_instantiation_is_unchanged();
    test_converges_to_the_true_lever_arm();
    test_estimated_arm_matches_a_surveyed_one();
    test_single_axis_rotation_leaves_the_axis_unobservable();
    test_converges_with_the_linear_block_running();
    test_frozen_updates_do_not_move_the_estimate();
    test_warmup_holds_rather_than_discards();
    test_projection_bounds_a_runaway_estimate();

    if (failures != 0) {
        std::cerr << failures << " lever-arm estimation check(s) failed\n";
        return 1;
    }
    std::cout << "all lever-arm estimation checks passed\n";
    return 0;
}
