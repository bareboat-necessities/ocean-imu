// A firmware built for another calibration site
// (-DATOMS3R_CALIBRATION_GRAVITY_MPS2): the physical gravity changes, the
// nominal-g -> m/s^2 conversion does not, and a calibration fitted against
// the default (Fair Lawn) gravity is refused rather than reinterpreted.
#define ATOMS3R_CALIBRATION_GRAVITY_MPS2 9.7803f
#include "../common/GravityChain.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {
int g_checks = 0, g_failures = 0;
void check(bool cond, const char* what) {
  ++g_checks;
  if (!cond) { ++g_failures; std::printf("CHECK FAILED: %s\n", what); }
}
}  // namespace

int main() {
  using namespace atoms3r_ical;
  namespace gc = gravity_chain;

  check(ImuCalCfg::g_cal_local == 9.7803f, "override sets g_cal_local");
  check(ImuCalCfg::g_std == 9.80665f, "override leaves g_std at 9.80665");
  const Eigen::Vector3f a = accel_nominal_g_to_body_ned_si_(0.1f, -0.2f, 1.01f);
  check(a.x() == -0.2f * 9.80665f && a.y() == 0.1f * 9.80665f && a.z() == -1.01f * 9.80665f,
        "nominal g -> m/s^2 conversion unchanged by g_cal_local");

  // Calibrated at this site against its gravity: physical, accepted.
  const auto c = gc::calibrate(gc::Sensor::typical(), 9.7803, ImuCalCfg::g_cal_local);
  const double scale = (c.fc.S.cast<double>() * gc::Sensor::typical().A.inverse()).trace() / 3.0;
  check(c.ok && c.rt.acc.ok && c.blob.accel_g == 9.7803f, "calibration fitted against the configured gravity is applied");
  check(std::fabs(scale - 1.0) < 1e-6, "its scale is physical");

  // A blob fitted against the default Fair Lawn gravity is not applied here.
  const auto fl = gc::calibrate(gc::Sensor::typical(), 9.7803, 9.8025605, 9.8025605f);
  RuntimeCals rt;
  rt.rebuildFromBlob(fl.blob);
  check(fl.ok && !rt.acc.ok && rt.accel_gravity_mismatch, "calibration fitted against another gravity is refused");

  std::printf("accel_cal_gravity_override-test: g_cal_local=%.7f g_std=%.5f scale=%.9f: %d/%d checks passed\n",
              (double)ImuCalCfg::g_cal_local, (double)ImuCalCfg::g_std, scale, g_checks - g_failures, g_checks);
  return g_failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
