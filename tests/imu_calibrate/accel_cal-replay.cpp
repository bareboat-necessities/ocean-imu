// Replays an AtomS3R accelerometer calibration serial log.
//
//   accel_cal-replay LOG [--blocks] [--verbose]
//
// With [ACCRAW] lines (firmware built with ATOMS3R_ICAL_RAW_LOG=1) the whole
// procedure (hold qualification, extras, rechecks, fit) is re-run from the
// samples; otherwise, or with --blocks, only the fit is re-run from the logged
// [ACCBLK] blocks. The replayed coefficients are compared with the device's
// logged "[ACCFIT] final" result.

#define EIGEN_NON_ARDUINO

#include "accel_cal_replay.h"

#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: accel_cal-replay LOG [--blocks] [--verbose]\n";
    return 2;
  }
  bool blocks_only = false, verbose = false;
  for (int i = 2; i < argc; ++i) {
    if (!strcmp(argv[i], "--blocks")) blocks_only = true;
    else if (!strcmp(argv[i], "--verbose")) verbose = true;
  }
  std::ifstream in(argv[1]);
  if (!in) { std::cerr << "cannot open " << argv[1] << "\n"; return 2; }
  std::vector<std::string> raw;
  for (std::string l; std::getline(in, l);) raw.push_back(l);
  const accel_replay::Log L = accel_replay::parse(raw);

  imu_cal::AccelFullFitResult r;
  bool ok = false;
  if (L.has_raw && !blocks_only) {
    int div = 0;
    ok = accel_replay::replayRaw(L, r, div, verbose);
    printf("mode=raw divergences=%d\n", div);
  } else if (L.has_blocks) {
    int n = 0;
    ok = accel_replay::fitBlocks(L, r, &n);
    printf("mode=blocks blocks=%d\n", n);
  } else {
    std::cerr << "no [ACCRAW] or [ACCBLK] lines found\n";
    return 2;
  }
  printf("ok=%d reason=%s gate=%s holds=%d blocks=%d thermal=%s/%s\n", (int)r.ok, imu_cal::fitFailStr(r.reason),
         imu_cal::accelGateStr(r.gate), r.n_fit_holds, r.n_fit_blocks, imu_cal::accelThermalStr(r.thermal),
         imu_cal::accelThermalReasonStr(r.thermal_reason));
  printf("b=(%.5f,%.5f,%.5f) m/s^2  bias_sd=(%.4f,%.4f,%.4f)\n", r.b(0), r.b(1), r.b(2), r.bias_sigma(0),
         r.bias_sigma(1), r.bias_sigma(2));
  printf("k=(%.6f,%.6f,%.6f) m/s^2/degC  T_ref=%.2f  T=[%.2f,%.2f]\n", r.k(0), r.k(1), r.k(2), r.T_ref,
         r.cal_temp_lo, r.cal_temp_hi);
  printf("S=[%.6f %.6f %.6f; %.6f %.6f %.6f; %.6f %.6f %.6f]\n", r.S(0, 0), r.S(0, 1), r.S(0, 2), r.S(1, 0),
         r.S(1, 1), r.S(1, 2), r.S(2, 0), r.S(2, 1), r.S(2, 2));
  printf("hold_rms=%.4f heldout: n=%d rms=%.4f max=%.4f t=%.2f  verify: n=%d rms=%.4f max=%.4f\n", r.hold_rms,
         r.n_cv, r.cv_rms, r.cv_max, r.cv_tmax, r.n_ver, r.ver_rms, r.ver_max);
  if (L.have_device_fit) {
    const double db = std::sqrt(std::pow(r.b(0) - L.dev_b[0], 2) + std::pow(r.b(1) - L.dev_b[1], 2) +
                                std::pow(r.b(2) - L.dev_b[2], 2));
    printf("device b=(%.5f,%.5f,%.5f)  |replay - device| = %.6f m/s^2\n", L.dev_b[0], L.dev_b[1], L.dev_b[2], db);
  }
  return ok ? 0 : 1;
}
