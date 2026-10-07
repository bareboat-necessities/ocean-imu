// Copyright 2026, Mikhail Grushinskiy
// Heading-only magnetometer update of the quaternion MEKF (compass use).
#define EIGEN_NON_ARDUINO
#include "ahrs/KalmanQMEKF.h"
#include <cmath>
#include <cstdio>

using V = Eigen::Vector3f;
using Q = Eigen::Quaternionf;
using MEKF = QuaternionMEKF<float, true>;
static const float D = 57.29578f;
static int failures = 0;
static void require(bool ok, const char* message) {
  if (!ok) { std::fprintf(stderr, "FAIL: %s\n", message); ++failures; }
}
static MEKF make() {
  return MEKF(V::Constant(0.06f * 9.80665f), V::Constant(0.003f), V::Constant(0.02f), 0.5f, 1e-2f, 1e-9f);
}
static Q att(const MEKF& f) { auto c = f.quaternion(); return Q(c(3), c(0), c(1), c(2)); }
static float yawDeg(const Q& q) { const auto R = q.toRotationMatrix(); return std::atan2(R(1, 0), R(0, 0)) * D; }
static float tiltDeg(const Q& a, const Q& b) {  // angle between the two body "down" directions
  const V da = a.toRotationMatrix().transpose() * V(0, 0, 1), db = b.toRotationMatrix().transpose() * V(0, 0, 1);
  return std::acos(std::fmin(1.0f, da.dot(db))) * D;
}

int main() {
  const float I = 66.5f / D;
  const V field(std::cos(I), 0, std::sin(I));
  const V acc_level(0, 0, -9.8f);

  // 1. A yaw offset is corrected and roll/pitch are left alone.
  {
    MEKF f = make();
    f.initialize_from_acc_mag(acc_level, field);
    const Q truth(Eigen::AngleAxisf(10.0f / D, V::UnitZ()));
    const V m = truth.toRotationMatrix().transpose() * field;
    for (int i = 0; i < 200; ++i) require(f.measurement_update_mag_heading(m, 0.02f), "update accepted");
    require(std::fabs(yawDeg(att(f)) - 10.0f) < 0.5f, "heading converges to the measured heading");
    require(tiltDeg(att(f), truth) < 0.05f, "heading update does not tilt the attitude");
  }
  // 2. A field whose dip differs from the reference (calibration residual or
  //    mag/accel misalignment) at east, tilted 20 deg: the 3D update tilts the
  //    attitude, the heading-only update does not.
  {
    const Q truth(Eigen::AngleAxisf(90.0f / D, V::UnitZ()) * Eigen::AngleAxisf(20.0f / D, V::UnitY()));
    const V wrong_dip = Eigen::AngleAxisf(6.0f / D, V::UnitY()) * field;  // dip 6 deg off
    const V m = truth.toRotationMatrix().transpose() * wrong_dip;
    const V a = truth.toRotationMatrix().transpose() * acc_level;
    MEKF heading = make(), full = make();
    heading.initialize_from_acc_mag(a, truth.toRotationMatrix().transpose() * field);
    full.initialize_from_acc_mag(a, truth.toRotationMatrix().transpose() * field);
    for (int i = 0; i < 200; ++i) {
      heading.measurement_update_acc_only(a);
      heading.measurement_update_mag_heading(m, 0.02f);
      full.measurement_update_acc_only(a);
      full.measurement_update_mag_only(m);
    }
    const float th = tiltDeg(att(heading), truth), tf = tiltDeg(att(full), truth);
    std::printf("qmekf-heading: dip mismatch 6 deg at east, 20 deg pitch: tilt error heading-only %.2f deg, 3D %.2f deg\n",
                th, tf);
    require(th < 0.1f, "heading-only update keeps tilt on the accelerometer");
    require(tf > 1.0f, "the 3D update is pulled by the dip mismatch (reference behaviour)");
  }
  // 3. No horizontal field (vertical measurement): rejected, nothing changes.
  {
    MEKF f = make();
    f.initialize_from_acc_mag(acc_level, field);
    const Q before = att(f);
    require(!f.measurement_update_mag_heading(V(0, 0, 1), 0.02f), "vertical field is rejected");
    require(att(f).angularDistance(before) < 1e-6f, "rejected update leaves the attitude unchanged");
  }
  if (failures) return 1;
  std::puts("qmekf-heading-test: OK");
  return 0;
}
