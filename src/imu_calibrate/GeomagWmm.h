#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  World Magnetic Model WMM-2025 (NOAA NCEI / BGS, epoch 2025.0, valid
  2025.0-2030.0): the main geomagnetic field expected at a place and date,
  for checking a magnetometer against it.

  The coefficients are the WMM-2025 release (WMM.COF dated 11/13/2024),
  transcribed by script from the copy in the pygeomag 1.1.0 package; the pywmm
  1.1.1 package carries identical rows. The model is in the public domain.
  The evaluation follows the WMM technical report: geodetic to geocentric
  coordinates, Schmidt semi-normalized associated Legendre functions, the
  spherical-harmonic sum to degree 12, and the rotation back to geodetic
  north/east/down. tests/imu_calibrate/geomag_wmm-test.cpp checks it against
  values from both packages.

  The model omits crustal and external fields: at a given place the measured
  field typically differs from it by a few hundred nT (about 1 %), more near
  steel, vehicles or buildings.
*/

#include <cmath>
#include <cstdlib>
#include <cstring>

namespace imu_cal {

struct GeomagField {
  double x_nT = NAN, y_nT = NAN, z_nT = NAN;  // north, east, down
  double h_nT = NAN, f_nT = NAN;               // horizontal, total
  double inclination_deg = NAN;                // dip, positive down
  double declination_deg = NAN;                // positive east of true north
  bool in_validity = false;                    // date within 2025.0-2030.0
};

struct Wmm2025 {
  static constexpr double EPOCH = 2025.0, VALID_UNTIL = 2030.0;
  static constexpr int N = 12;
  struct Coef { int n, m; double g, h, gdot, hdot; };
  static constexpr Coef COEF[] = {
      { 1,  0,  -29351.8,       0.0,   12.0,    0.0},
      { 1,  1,   -1410.8,    4545.4,    9.7,  -21.5},
      { 2,  0,   -2556.6,       0.0,  -11.6,    0.0},
      { 2,  1,    2951.1,   -3133.6,   -5.2,  -27.7},
      { 2,  2,    1649.3,    -815.1,   -8.0,  -12.1},
      { 3,  0,    1361.0,       0.0,   -1.3,    0.0},
      { 3,  1,   -2404.1,     -56.6,   -4.2,    4.0},
      { 3,  2,    1243.8,     237.5,    0.4,   -0.3},
      { 3,  3,     453.6,    -549.5,  -15.6,   -4.1},
      { 4,  0,     895.0,       0.0,   -1.6,    0.0},
      { 4,  1,     799.5,     278.6,   -2.4,   -1.1},
      { 4,  2,      55.7,    -133.9,   -6.0,    4.1},
      { 4,  3,    -281.1,     212.0,    5.6,    1.6},
      { 4,  4,      12.1,    -375.6,   -7.0,   -4.4},
      { 5,  0,    -233.2,       0.0,    0.6,    0.0},
      { 5,  1,     368.9,      45.4,    1.4,   -0.5},
      { 5,  2,     187.2,     220.2,    0.0,    2.2},
      { 5,  3,    -138.7,    -122.9,    0.6,    0.4},
      { 5,  4,    -142.0,      43.0,    2.2,    1.7},
      { 5,  5,      20.9,     106.1,    0.9,    1.9},
      { 6,  0,      64.4,       0.0,   -0.2,    0.0},
      { 6,  1,      63.8,     -18.4,   -0.4,    0.3},
      { 6,  2,      76.9,      16.8,    0.9,   -1.6},
      { 6,  3,    -115.7,      48.8,    1.2,   -0.4},
      { 6,  4,     -40.9,     -59.8,   -0.9,    0.9},
      { 6,  5,      14.9,      10.9,    0.3,    0.7},
      { 6,  6,     -60.7,      72.7,    0.9,    0.9},
      { 7,  0,      79.5,       0.0,   -0.0,    0.0},
      { 7,  1,     -77.0,     -48.9,   -0.1,    0.6},
      { 7,  2,      -8.8,     -14.4,   -0.1,    0.5},
      { 7,  3,      59.3,      -1.0,    0.5,   -0.8},
      { 7,  4,      15.8,      23.4,   -0.1,    0.0},
      { 7,  5,       2.5,      -7.4,   -0.8,   -1.0},
      { 7,  6,     -11.1,     -25.1,   -0.8,    0.6},
      { 7,  7,      14.2,      -2.3,    0.8,   -0.2},
      { 8,  0,      23.2,       0.0,   -0.1,    0.0},
      { 8,  1,      10.8,       7.1,    0.2,   -0.2},
      { 8,  2,     -17.5,     -12.6,    0.0,    0.5},
      { 8,  3,       2.0,      11.4,    0.5,   -0.4},
      { 8,  4,     -21.7,      -9.7,   -0.1,    0.4},
      { 8,  5,      16.9,      12.7,    0.3,   -0.5},
      { 8,  6,      15.0,       0.7,    0.2,   -0.6},
      { 8,  7,     -16.8,      -5.2,   -0.0,    0.3},
      { 8,  8,       0.9,       3.9,    0.2,    0.2},
      { 9,  0,       4.6,       0.0,   -0.0,    0.0},
      { 9,  1,       7.8,     -24.8,   -0.1,   -0.3},
      { 9,  2,       3.0,      12.2,    0.1,    0.3},
      { 9,  3,      -0.2,       8.3,    0.3,   -0.3},
      { 9,  4,      -2.5,      -3.3,   -0.3,    0.3},
      { 9,  5,     -13.1,      -5.2,    0.0,    0.2},
      { 9,  6,       2.4,       7.2,    0.3,   -0.1},
      { 9,  7,       8.6,      -0.6,   -0.1,   -0.2},
      { 9,  8,      -8.7,       0.8,    0.1,    0.4},
      { 9,  9,     -12.9,      10.0,   -0.1,    0.1},
      {10,  0,      -1.3,       0.0,    0.1,    0.0},
      {10,  1,      -6.4,       3.3,    0.0,    0.0},
      {10,  2,       0.2,       0.0,    0.1,   -0.0},
      {10,  3,       2.0,       2.4,    0.1,   -0.2},
      {10,  4,      -1.0,       5.3,   -0.0,    0.1},
      {10,  5,      -0.6,      -9.1,   -0.3,   -0.1},
      {10,  6,      -0.9,       0.4,    0.0,    0.1},
      {10,  7,       1.5,      -4.2,   -0.1,    0.0},
      {10,  8,       0.9,      -3.8,   -0.1,   -0.1},
      {10,  9,      -2.7,       0.9,   -0.0,    0.2},
      {10, 10,      -3.9,      -9.1,   -0.0,   -0.0},
      {11,  0,       2.9,       0.0,    0.0,    0.0},
      {11,  1,      -1.5,       0.0,   -0.0,   -0.0},
      {11,  2,      -2.5,       2.9,    0.0,    0.1},
      {11,  3,       2.4,      -0.6,    0.0,   -0.0},
      {11,  4,      -0.6,       0.2,    0.0,    0.1},
      {11,  5,      -0.1,       0.5,   -0.1,   -0.0},
      {11,  6,      -0.6,      -0.3,    0.0,   -0.0},
      {11,  7,      -0.1,      -1.2,   -0.0,    0.1},
      {11,  8,       1.1,      -1.7,   -0.1,   -0.0},
      {11,  9,      -1.0,      -2.9,   -0.1,    0.0},
      {11, 10,      -0.2,      -1.8,   -0.1,    0.0},
      {11, 11,       2.6,      -2.3,   -0.1,    0.0},
      {12,  0,      -2.0,       0.0,    0.0,    0.0},
      {12,  1,      -0.2,      -1.3,    0.0,   -0.0},
      {12,  2,       0.3,       0.7,   -0.0,    0.0},
      {12,  3,       1.2,       1.0,   -0.0,   -0.1},
      {12,  4,      -1.3,      -1.4,   -0.0,    0.1},
      {12,  5,       0.6,      -0.0,   -0.0,   -0.0},
      {12,  6,       0.6,       0.6,    0.1,   -0.0},
      {12,  7,       0.5,      -0.1,   -0.0,   -0.0},
      {12,  8,      -0.1,       0.8,    0.0,    0.0},
      {12,  9,      -0.4,       0.1,    0.0,   -0.0},
      {12, 10,      -0.2,      -1.0,   -0.1,   -0.0},
      {12, 11,      -1.3,       0.1,   -0.0,    0.0},
      {12, 12,      -0.7,       0.2,   -0.1,   -0.1},
  };

  // lat/lon in degrees (geodetic, east positive), height in metres above the
  // WGS84 ellipsoid, date as a decimal year.
  static GeomagField field(double lat_deg, double lon_deg, double height_m, double year) {
    GeomagField out;
    if (!std::isfinite(lat_deg) || !std::isfinite(lon_deg) || !std::isfinite(height_m) || !std::isfinite(year) ||
        std::fabs(lat_deg) > 90.0)
      return out;
    constexpr double PI_ = 3.14159265358979323846, D2R_ = PI_ / 180.0;
    constexpr double A = 6378.137, F = 1.0 / 298.257223563, E2 = F * (2.0 - F), RE = 6371.2;
    const double dt = year - EPOCH;
    const double phi = lat_deg * D2R_, lam = lon_deg * D2R_, h = height_m * 1e-3;

    // Geodetic to geocentric spherical (km).
    const double sp = std::sin(phi), cp = std::cos(phi);
    const double rc = A / std::sqrt(1.0 - E2 * sp * sp);
    const double p = (rc + h) * cp, z = (rc * (1.0 - E2) + h) * sp;
    const double r = std::sqrt(p * p + z * z);
    const double phic = std::asin(z / r);

    // Schmidt semi-normalized P(n,m)(cos theta) and dP/dtheta, theta the
    // geocentric colatitude.
    const double ct = std::sin(phic), st = std::cos(phic);
    double P[N + 1][N + 1] = {}, dP[N + 1][N + 1] = {};
    P[0][0] = 1.0;
    for (int n = 1; n <= N; ++n) {
      for (int m = 0; m <= n; ++m) {
        if (n == m) {
          const double k = n == 1 ? 1.0 : std::sqrt(1.0 - 1.0 / (2.0 * n));
          P[n][n] = k * st * P[n - 1][n - 1];
          dP[n][n] = k * (st * dP[n - 1][n - 1] + ct * P[n - 1][n - 1]);
        } else {
          const double a = (2.0 * n - 1.0) / std::sqrt(double(n * n - m * m));
          const double b = n - 1 > m ? std::sqrt(double((n - 1) * (n - 1) - m * m)) / std::sqrt(double(n * n - m * m)) : 0.0;
          P[n][m] = a * ct * P[n - 1][m] - (n - 1 > m ? b * P[n - 2][m] : 0.0);
          dP[n][m] = a * (ct * dP[n - 1][m] - st * P[n - 1][m]) - (n - 1 > m ? b * dP[n - 2][m] : 0.0);
        }
      }
    }

    double bx = 0, by = 0, bz = 0;  // geocentric north, east, down
    for (const Coef& c : COEF) {
      const double g = c.g + dt * c.gdot, hh = c.h + dt * c.hdot;
      const double ar = std::pow(RE / r, c.n + 2);
      const double cm = std::cos(c.m * lam), sm = std::sin(c.m * lam);
      const double t1 = g * cm + hh * sm, t2 = g * sm - hh * cm;
      bx += ar * t1 * dP[c.n][c.m];
      by += ar * c.m * t2 * P[c.n][c.m];
      bz -= ar * (c.n + 1) * t1 * P[c.n][c.m];
    }
    by = st > 1e-10 ? by / st : 0.0;  // the poles are outside this use

    // Back to geodetic north/down.
    const double psi = phic - phi;
    out.x_nT = bx * std::cos(psi) - bz * std::sin(psi);
    out.y_nT = by;
    out.z_nT = bx * std::sin(psi) + bz * std::cos(psi);
    out.h_nT = std::hypot(out.x_nT, out.y_nT);
    out.f_nT = std::hypot(out.h_nT, out.z_nT);
    out.inclination_deg = std::atan2(out.z_nT, out.h_nT) / D2R_;
    out.declination_deg = std::atan2(out.y_nT, out.x_nT) / D2R_;
    out.in_validity = year >= EPOCH && year <= VALID_UNTIL;
    return out;
  }
};

// Decimal year of a compiler __DATE__ string ("Oct 13 2026"); NAN if it does
// not parse.
inline double decimalYearFromBuildDate(const char* date) {
  static const char* const MONTHS = "JanFebMarAprMayJunJulAugSepOctNovDec";
  if (!date || std::strlen(date) < 11) return NAN;
  int month = -1;
  for (size_t i = 0; i < 12; ++i)
    if (std::strncmp(MONTHS + 3 * i, date, 3) == 0) month = int(i);
  char* end = nullptr;
  const long day = std::strtol(date + 4, &end, 10);
  if (end != date + 6) return NAN;
  const long year = std::strtol(date + 7, &end, 10);
  if (end != date + 11) return NAN;
  if (month < 0 || day < 1 || day > 31 || year < 2000) return NAN;
  static const int CUM[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  const bool leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
  const long doy = CUM[month] + (month > 1 && leap ? 1 : 0) + day - 1;
  return double(year) + double(doy) / (leap ? 366.0 : 365.0);
}

}  // namespace imu_cal
