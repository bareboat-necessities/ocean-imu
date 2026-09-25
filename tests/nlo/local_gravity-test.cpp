// Local-gravity convention through the NLO adapter: a still sensor calibrated
// against g_local (as the wizard does) and fed as the sketch feeds it. The
// measured specific force resolved with the NLO's own attitude plus its
// configured gravity, q_nb * f + (0,0,g), is ~0 with g_local and
// g_std - g_local with 9.80665. (Its estimated-specific-force output is pulled
// to the configured gravity by the aiding, so it cannot show the mismatch.)
#include "../common/GravityChain.h"
#include "nlo/TimeVarGainNLO_Adapter.h"
#include <cmath>
#include <cstdio>
extern const float g_std = 9.80665f;  // generic default of the shared headers

using Fusion = TimeVarGainNloAdapter<false, NloMagType::None, float>;

double tailVertical(const Eigen::Vector3f& a_cal_ned, float g_filter) {
  Fusion f;
  auto& cfg = f.config();
  cfg.gravity_mps2 = g_filter;
  cfg.filter.gravity_mps2 = g_filter;
  f.reset();
  double s = 0; int n = 0;
  for (int k = 0; k < 24000; ++k) {  // 120 s at 200 Hz
    f.update(0.005f, Eigen::Vector3f::Zero(), a_cal_ned);
    if (k >= 20000) {
      const Eigen::Vector3f e = f.snapshot().q_nb * a_cal_ned + Eigen::Vector3f(0, 0, g_filter);
      s += e.z(); ++n;
    }
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
  std::printf("NLO g_local=%.7f e_z(g_local)=%+.7f e_z(g_std)=%+.7f expected(g_std)=%+.7f\n",
              g_local, v_local, v_std, expected);
  bool ok = cal.ok && cal.rt.acc.ok && std::isfinite(v_local) && std::isfinite(v_std);
  ok = ok && std::fabs(v_local) < 2e-4 && std::fabs(v_std - expected) < 2e-4;
  if (!ok) std::printf("FAIL: NLO does not remove the configured gravity\n");
  return ok ? 0 : 1;
}
