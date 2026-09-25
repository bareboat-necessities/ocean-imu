#pragma once
// Device accelerometer chain for the gravity-convention tests:
//   physical specific force (body NED, m/s^2)
//   -> sensor reading (BMI270/M5Unified, sensor axes, nominal g)
//   -> x g_std, axis mapping (accel_nominal_g_to_body_ned_si_)
//   -> accelerometer calibration fitted against a configured gravity
//   -> blob -> RuntimeCals (what the sketches apply upstream of every filter).
// A static hold at a site with gravity g_true reads |f| = g_true.
#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include <memory>
#include <vector>
#include "AtomS3R/AtomS3R_ImuCalBlob.h"

namespace gravity_chain {

using Vec3d = Eigen::Vector3d;
using Mat3d = Eigen::Matrix3d;
using Vec3f = Eigen::Vector3f;
using atoms3r_ical::ImuCalCfg;

// Sensor truth: reading_SI = A_inv * f + b (A symmetric: scale and cross axis).
struct Sensor {
  Mat3d A = Mat3d::Identity();
  Vec3d b = Vec3d::Zero();
  static Sensor perfect() { return Sensor{}; }
  static Sensor typical() {
    Sensor s;
    s.A << 1.004, 0.002, -0.001,
           0.002, 0.997, 0.003,
          -0.001, 0.003, 1.002;
    s.b = Vec3d(0.08, -0.12, 0.05);
    return s;
  }
};

// Sensor output in nominal g on the sensor axes, for body-NED specific force f.
inline Vec3f nominalG(const Sensor& s, const Vec3d& f_body) {
  const Vec3d a = s.A.inverse() * f_body + s.b;   // SI reading, body NED
  // Inverse of (sx, sy, sz) -> (sy, sx, -sz) * g_std.
  return Vec3f(float(a.y() / ImuCalCfg::g_std), float(a.x() / ImuCalCfg::g_std),
               float(-a.z() / ImuCalCfg::g_std));
}

// The firmware's conversion: nominal g -> body NED m/s^2.
inline Vec3f toSI(const Vec3f& n) { return atoms3r_ical::accel_nominal_g_to_body_ned_si_(n.x(), n.y(), n.z()); }

inline Vec3f readingSI(const Sensor& s, const Vec3d& f_body) { return toSI(nominalG(s, f_body)); }

// Hold directions ("up" in body axes): six faces and the eight cube corners.
inline std::vector<Vec3d> holdDirections() {
  std::vector<Vec3d> u = {Vec3d(1, 0, 0), Vec3d(-1, 0, 0), Vec3d(0, 1, 0),
                          Vec3d(0, -1, 0), Vec3d(0, 0, 1), Vec3d(0, 0, -1)};
  for (int i = 0; i < 8; ++i)
    u.push_back(Vec3d((i & 1) ? 1 : -1, (i & 2) ? 1 : -1, (i & 4) ? 1 : -1).normalized());
  return u;
}

// Noise-free hold means of `s` at a site with gravity g_true.
inline std::vector<imu_cal::AccelObs> holds(const Sensor& s, double g_true) {
  std::vector<imu_cal::AccelObs> obs;
  const auto dirs = holdDirections();
  uint32_t t = 0;
  for (size_t h = 0; h < dirs.size(); ++h) {
    const Vec3f a = readingSI(s, g_true * dirs[h]);
    for (int k = 0; k < 20; ++k) {
      imu_cal::AccelObs o{};
      for (int d = 0; d < 3; ++d) { o.a[d] = a(d); o.w[d] = 0.0f; }
      o.tempC = 25.0f; o.a_std = 0.002f; o.t_ms = (t += 250); o.n = 25;
      o.hold = (uint8_t)h; o.role = (uint8_t)imu_cal::AccelObsRole::FIT;
      obs.push_back(o);
    }
  }
  return obs;
}

using Fitter = imu_cal::AccelFullFitter<340, 24>;

struct Calibration {
  bool ok = false;
  imu_cal::AccelFullFitResult r;
  imu_cal::AccelCalibration<float> fc;
  atoms3r_ical::ImuCalBlobV3 blob;
  atoms3r_ical::RuntimeCals rt;
};

// Wizard path: fit against g_cfg, cast to float, fill the blob, rebuild the
// runtime calibration the sketches apply (accepted only if the blob's gravity
// matches g_runtime, the firmware's g_cal_local).
inline Calibration calibrate(const Sensor& s, double g_true, double g_cfg,
                             float g_runtime = ImuCalCfg::g_cal_local) {
  Calibration c;
  const auto obs = holds(s, g_true);
  imu_cal::AccelFitCfg cfg;
  cfg.g = g_cfg;
  imu_cal::AccelThermalPrior prior;
  auto fitter = std::make_unique<Fitter>();
  if (!fitter->fit(obs.data(), (int)obs.size(), cfg, prior, c.r, false)) return c;
  Fitter::toFloat(c.r, cfg.g, c.fc);
  atoms3r_ical::fillAccelFromFit(c.blob, c.r, c.fc, 70);
  c.rt.rebuildFromBlob(c.blob, g_runtime);
  c.ok = c.r.ok && c.fc.ok;
  return c;
}

} // namespace gravity_chain
