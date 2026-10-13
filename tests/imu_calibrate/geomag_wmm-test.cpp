// WMM-2025 evaluation (src/imu_calibrate/GeomagWmm.h) against reference
// values, and the build-date parser the diagnostics sketch dates it with.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "imu_calibrate/GeomagWmm.h"

static int failures = 0;
static void require(bool ok, const char* what) {
  std::printf("  %-66s %s\n", what, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

struct Reference { double lat, lon, h_m, year, x, y, z, inc, dec; };

// Computed with pygeomag 1.1.0 (its WMM_2025.COF); pywmm 1.1.1 agrees with
// the same evaluation to 2e-6 nT.
static const Reference REFERENCES[] = {
  {40.935833, -74.117504, -9.0, 2026.0, 20396.563, -4515.094, 46493.843, 65.80495, -12.48202},
  {0.0, 0.0, 0.0, 2025.0, 27453.929, -1927.587, -16010.807, -30.18896, -4.01624},
  {-33.9, 18.4, 10.0, 2027.5, 9565.401, -4860.533, -22503.703, -64.50879, -26.93680},
  {64.1, -21.9, 50.0, 2028.25, 13051.954, -2436.965, 50892.052, 75.37773, -10.57607},
  {35.7, 139.7, 40.0, 2029.5, 30137.232, -4259.629, 35788.107, 49.61978, -8.04496},
  {-45.0, 170.0, 2000.0, 2026.0, 17901.103, 8513.815, -54945.810, -70.16228, 25.43589},
  {80.0, 0.0, 0.0, 2025.0, 6521.599, 145.887, 54791.508, 83.21058, 1.28148},
  {-80.0, 240.0, 100000.0, 2027.5, 5983.976, 14760.136, -49317.671, -72.10228, 67.93164},
};

int main() {
  std::printf("field\n");
  double worst_nT = 0, worst_deg = 0;
  bool all_valid = true;
  for (const Reference& r : REFERENCES) {
    const imu_cal::GeomagField f = imu_cal::Wmm2025::field(r.lat, r.lon, r.h_m, r.year);
    worst_nT = std::max({worst_nT, std::fabs(f.x_nT - r.x), std::fabs(f.y_nT - r.y), std::fabs(f.z_nT - r.z)});
    worst_deg = std::max({worst_deg, std::fabs(f.inclination_deg - r.inc), std::fabs(f.declination_deg - r.dec)});
    all_valid = all_valid && f.in_validity &&
                std::fabs(f.f_nT - std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z)) < 0.01;
  }
  std::printf("  worst component error %.4f nT, worst angle error %.6f deg\n", worst_nT, worst_deg);
  require(worst_nT < 0.005, "X, Y, Z within 0.005 nT of the reference values");
  require(worst_deg < 1e-5, "inclination and declination within 1e-5 deg");
  require(all_valid, "total intensity consistent and dates inside 2025.0-2030.0");
  require(!imu_cal::Wmm2025::field(40, -74, 0, 2031.0).in_validity, "a date after 2030.0 is flagged outside validity");
  require(std::isnan(imu_cal::Wmm2025::field(NAN, 0, 0, 2026).f_nT) &&
              std::isnan(imu_cal::Wmm2025::field(91, 0, 0, 2026).f_nT),
          "invalid latitude gives no field");

  std::printf("build date\n");
  require(std::fabs(imu_cal::decimalYearFromBuildDate("Jan  1 2026") - 2026.0) < 1e-12, "Jan 1 is the start of the year");
  require(std::fabs(imu_cal::decimalYearFromBuildDate("Oct 13 2026") - (2026.0 + 285.0 / 365.0)) < 1e-12,
          "Oct 13 2026 is day 286 of 365");
  require(std::fabs(imu_cal::decimalYearFromBuildDate("Mar  1 2028") - (2028.0 + 60.0 / 366.0)) < 1e-12,
          "a leap year counts Feb 29");
  require(std::isnan(imu_cal::decimalYearFromBuildDate("Foo 13 2026")) &&
              std::isnan(imu_cal::decimalYearFromBuildDate("")) && std::isnan(imu_cal::decimalYearFromBuildDate(nullptr)),
          "an unparsable date gives NAN");

  std::printf("geomag_wmm: %s\n", failures ? "FAIL" : "PASS");
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
