// Local-gravity convention through the PII observer: a still sensor calibrated
// against g_local (as the wizard does) and fed as the sketch feeds it. The
// observer's gravity-removed vertical acceleration moves by exactly
// g_std - g_local when it is configured with 9.80665 instead of g_local.
#include "../common/GravityChain.h"
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
#include <cmath>
#include <cstdio>
extern const float g_std = 9.80665f;  // generic default of the shared headers

using Fusion = marine_obs::AdaptiveVerticalPIIMahony<float, true, TrackerType::PLL>;

double tailVertical(const Eigen::Vector3f& a_cal_ned, float g_filter) {
  Fusion::Config cfg;
  cfg.gravity_mps2 = g_filter;
  Fusion f(cfg);
  const Eigen::Vector3f a(a_cal_ned.x(), -a_cal_ned.y(), -a_cal_ned.z());  // ned_to_mahony_body_
  double s = 0; int n = 0;
  for (int k = 0; k < 24000; ++k) {  // 120 s at 200 Hz
    f.updateIMU(0.0f, 0.0f, 0.0f, a.x(), a.y(), a.z(), 0.005f);
    if (k >= 20000) { s += f.verticalWorldAccelUp(); ++n; }
  }
  return s / n;
}

int main() {
  using atoms3r_ical::ImuCalCfg;
  namespace gc = gravity_chain;
  const double g_local = ImuCalCfg::g_cal_local;
  const auto sensor = gc::Sensor::typical();
  const auto cal = gc::calibrate(sensor, g_local, g_local);
  const Eigen::Quaterniond q = Eigen::AngleAxisd(0.5, gc::Vec3d::UnitZ()) *
      Eigen::AngleAxisd(0.1, gc::Vec3d::UnitY()) * Eigen::AngleAxisd(0.15, gc::Vec3d::UnitX());
  const Eigen::Vector3f a_cal = cal.rt.applyAccel(gc::readingSI(sensor, q.conjugate() * gc::Vec3d(0, 0, -g_local)), 25.0f);
  const double v_local = tailVertical(a_cal, ImuCalCfg::g_cal_local);
  const double v_std = tailVertical(a_cal, ImuCalCfg::g_std);
  const double expected = (double)ImuCalCfg::g_std - g_local;
  std::printf("PII g_local=%.7f vertical(g_local)=%+.6f vertical(g_std)=%+.6f shift=%+.7f expected=%+.7f |a_cal|-g=%+.2e\n",
              g_local, v_local, v_std, v_local - v_std, expected, (double)a_cal.norm() - g_local);
  bool ok = cal.ok && cal.rt.acc.ok && std::isfinite(v_local) && std::isfinite(v_std);
  ok = ok && std::fabs((v_local - v_std) - expected) < 2e-5;
  ok = ok && std::fabs((double)a_cal.norm() - g_local) < 1e-5;
  if (!ok) std::printf("FAIL: PII does not remove the configured gravity\n");
  return ok ? 0 : 1;
}
