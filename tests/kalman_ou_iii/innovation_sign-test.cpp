// Pins that every shipping OU-III measurement update corrects the state TOWARD
// its measurement.
//
// This closes a gap the certified P1-P3 chain cannot close by construction.  The
// covariance update
//
//     P <- P - P C' (C P C' + R)^-1 C P
//
// never references the innovation, so a wrong-signed innovation would leave the
// segment floor, Sigma_upper and the certified delta bit-identical while the
// deployed filter drove its own state away from every measurement.  The
// source-marker contracts pin the update *shape* -- "xext.noalias() += K * r;"
// -- but nothing checks how r is defined.  P4 would catch it because it bounds
// the nonlinear state map, but P4 has never been reached because P3 blocks.
//
// The invariant used here is sign-sensitive and gain-free: a correct linear
// Kalman update makes the posterior a convex combination of the prior and the
// measurement,
//
//     x+ = x- + K (z - H x-),   with the relevant gain block in [0, 1),
//
// so the corrected component must lie in the half-open interval [x-, z) and can
// never move away from z.  Flipping the sign of r sends it the other way by
// exactly the same magnitude, which every assertion below rejects.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private

using Filter = Kalman3D_Wave_OU_III<float>;
using Vec3 = Eigen::Vector3f;

// The offsets are private by default inside the class (they precede any access
// specifier, so "#define private public" cannot reach them).  Restated locally,
// matching ou3-certificate-sim.cpp, and guarded at runtime below against the
// actual extended state size so the layout cannot drift silently.
constexpr int kNX    = 21;   // 6 base (att_err + gyro bias) + v,p,S,a_w,b_acc
constexpr int kOffV  = 6;
constexpr int kOffP  = 9;
constexpr int kOffS  = 12;

static int fail(const char* m) { std::cerr << "FAIL: " << m << "\n"; return 1; }

// Corrected value must move toward z and must not overshoot past it.
static bool moved_toward(float prior, float post, float z, const char* what) {
    const float d_before = std::fabs(z - prior);
    const float d_after  = std::fabs(z - post);
    if (!(d_after < d_before)) {
        std::cerr << "  " << what << ": |z-post|=" << d_after
                  << " did not shrink from |z-prior|=" << d_before
                  << " (prior=" << prior << " post=" << post << " z=" << z << ")\n";
        return false;
    }
    // Convexity: post lies between prior and z.
    const float lo = std::min(prior, z), hi = std::max(prior, z);
    if (post < lo - 1e-6f || post > hi + 1e-6f) {
        std::cerr << "  " << what << ": post=" << post
                  << " left the convex hull [" << lo << ", " << hi << "]\n";
        return false;
    }
    return true;
}

static Filter make() {
    Filter f(Vec3::Constant(0.0148f), Vec3::Constant(0.00157f), Vec3::Constant(0.25f),
             1e-4f, 1e-6f, 1e-4f, 0.5f, 9.80665f);
    // Give the translation blocks a well-conditioned prior so every gain is
    // strictly positive and the direction of motion is unambiguous.
    for (int i = 0; i < kNX; ++i) f.Pext(i, i) = std::max(f.Pext(i, i), 1.0f);
    return f;
}

int main() {
    int rc = 0;

    // Layout guard: if the extended state changes size, every offset below is
    // suspect and this test must be revisited rather than silently pass.
    {
        Filter probe = make();
        if (probe.xext.size() != kNX) {
            std::cerr << "  extended state is " << probe.xext.size()
                      << " states, expected " << kNX << "\n";
            return fail("OU-III extended state layout moved; offsets here are stale");
        }
    }

    // ---- S = 0 integral pseudo-measurement: target is exactly zero ----
    {
        Filter f = make();
        const Vec3 S0(0.7f, -1.3f, 2.1f);
        f.xext.segment<3>(kOffS) = S0;
        f.applyIntegralZeroPseudoMeas();
        const Vec3 S1 = f.xext.segment<3>(kOffS);
        for (int i = 0; i < 3; ++i)
            if (!moved_toward(S0(i), S1(i), 0.0f, "applyIntegralZeroPseudoMeas"))
                rc = fail("S=0 pseudo-measurement did not pull S toward zero");
    }

    // ---- position pseudo-measurement ----
    {
        Filter f = make();
        const Vec3 p0(1.0f, -2.0f, 3.0f), z(-0.5f, 0.25f, -1.5f);
        f.xext.segment<3>(kOffP) = p0;
        f.measurement_update_position_pseudo(z, Vec3::Constant(0.3f));
        const Vec3 p1 = f.xext.segment<3>(kOffP);
        for (int i = 0; i < 3; ++i)
            if (!moved_toward(p0(i), p1(i), z(i), "measurement_update_position_pseudo"))
                rc = fail("position pseudo-measurement moved position away from z");
    }

    // ---- velocity pseudo-measurement ----
    {
        Filter f = make();
        const Vec3 v0(0.9f, 1.4f, -0.6f), z(-0.2f, -1.1f, 0.8f);
        f.xext.segment<3>(kOffV) = v0;
        f.measurement_update_velocity_pseudo(z, Vec3::Constant(0.2f));
        const Vec3 v1 = f.xext.segment<3>(kOffV);
        for (int i = 0; i < 3; ++i)
            if (!moved_toward(v0(i), v1(i), z(i), "measurement_update_velocity_pseudo"))
                rc = fail("velocity pseudo-measurement moved velocity away from z");
    }

    // ---- vertical velocity pseudo-measurement (scalar path) ----
    {
        Filter f = make();
        const int idx = kOffV + 2;
        const float v0 = 1.7f, z = -0.4f;
        f.xext(idx) = v0;
        f.measurement_update_vert_velocity_pseudo(z, 0.25f);
        if (!moved_toward(v0, f.xext(idx), z, "measurement_update_vert_velocity_pseudo"))
            rc = fail("vertical velocity pseudo-measurement moved vz away from z");
    }

    // ---- accelerometer: a pure roll error must shrink, not grow ----
    {
        Filter f = make();
        const float tilt = 0.20f;                  // rad, about x
        f.qref = Eigen::Quaternionf(Eigen::AngleAxisf(tilt, Vec3::UnitX()));
        f.qref.normalize();
        // Specific force for a level vehicle, expressed in the body frame of the
        // TRUE (level) attitude: gravity reaction along -z in NED.
        const Vec3 acc_body(0.0f, 0.0f, -9.80665f);
        f.measurement_update_acc_only(acc_body);
        const float tilt_after =
            2.0f * std::atan2(f.qref.vec().norm(), std::fabs(f.qref.w()));
        if (!(tilt_after < tilt)) {
            std::cerr << "  accelerometer: tilt grew from " << tilt
                      << " to " << tilt_after << "\n";
            rc = fail("accelerometer update did not reduce the attitude error");
        }
    }

    if (rc == 0)
        std::cout << "OU-III innovation sign: every shipping update corrects toward its measurement\n";
    return rc;
}
