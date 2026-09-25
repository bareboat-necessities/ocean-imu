#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R units, axis mapping and the two gravity constants (host-testable,
  no Arduino dependency).

  Two different quantities, deliberately kept apart:

  g_std        Exact conventional standard gravity, 9.80665 m/s^2. M5Unified
               reports the BMI270 acceleration in nominal g; multiplying by
               g_std converts that unit to m/s^2. It says nothing about the
               gravity at any real place, and is used for nothing else.

  g_cal_local  Physical gravity assumed where the accelerometer is calibrated.
               The calibration makes the static specific-force norm equal to
               it, and every AtomS3R estimator removes this same value, so a
               still sensor gives zero translational acceleration and dynamic
               acceleration keeps its SI scale.

  Default g_cal_local: WGS84 normal gravity (NGA.STND.0036 / TR8350.2:
  Somigliana closed form on the ellipsoid, then the second-order height
  correction of eq. 4-3) at the project's default calibration site,
  Fair Lawn NJ 07410: geodetic latitude 40.935833 N (longitude 74.117504 W
  does not enter the normal-gravity formula), WGS84 ellipsoidal height
  h = -9 m. The height correction is defined in ellipsoidal height, so -9 m
  is used; the +25 m orthometric elevation would give 1.05e-4 m/s^2 less.
  Result: 9.8025605 m/s^2. This is a model value, not a measurement: the
  real gravity differs from normal gravity by the local gravity disturbance
  (typically some 1e-4 m/s^2), while height (+-5 m) and latitude uncertainty
  contribute about 1.5e-5 m/s^2.

  Build with -DATOMS3R_CALIBRATION_GRAVITY_MPS2=<value> to calibrate and run
  elsewhere. A saved accelerometer calibration records the gravity it was
  fitted against (ImuCalBlobV3::accel_g) and is not applied when that differs
  from g_cal_local.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#ifndef ATOMS3R_CALIBRATION_GRAVITY_MPS2
#define ATOMS3R_CALIBRATION_GRAVITY_MPS2 9.8025605f  // Fair Lawn NJ, WGS84 normal gravity, h = -9 m
#endif

namespace atoms3r_ical {

// Config/constants (shared)
struct ImuCalCfg {
  // Nominal g -> m/s^2 (unit conversion only).
  static constexpr float g_std = 9.80665f;
  // Physical gravity at the calibration site: calibration target and the
  // gravity every AtomS3R estimator removes.
  static constexpr float g_cal_local = ATOMS3R_CALIBRATION_GRAVITY_MPS2;
  static constexpr float DEG2RAD = 3.14159265358979323846f / 180.0f;
};

// Default calibration site (documentation of g_cal_local's default).
struct CalibrationSite {
  static constexpr double latitude_deg  = 40.935833;   // geodetic, N
  static constexpr double longitude_deg = -74.117504;  // E positive; not used by normal gravity
  static constexpr double ellipsoidal_height_m = -9.0; // WGS84, used by the model
  static constexpr double orthometric_height_m = 25.0; // above mean sea level, not used
};

// Axis mapping (AtomS3R)
//
// Internal convention in this library is BODY-NED (x=north, y=east, z=down).
// End-user "nautical Z-up" is therefore z_up = -z_down.
//
// With the board lying still, screen facing up:
//   - accelerometer body Z (down) is expected near -g specific force
//   - user-facing Z-up value is the opposite sign.
// acc_body = ( ay, ax, -az ) * g_std
// gyr_body = ( gy, gx, -gz ) * deg2rad
// mag_body = ( my, mx, -mz ) * (1/10)
static inline Eigen::Matrix<float,3,1> map_sensor_xyz_to_body_ned_(float sx, float sy, float sz, float scale = 1.0f) {
  return Eigen::Matrix<float,3,1>(sy * scale, sx * scale, -sz * scale);
}

// M5Unified accelerometer sample (nominal g, sensor axes) -> body NED specific
// force in m/s^2. Independent of g_cal_local.
static inline Eigen::Matrix<float,3,1> accel_nominal_g_to_body_ned_si_(float gx, float gy, float gz) {
  return map_sensor_xyz_to_body_ned_(gx, gy, gz, ImuCalCfg::g_std);
}

} // namespace atoms3r_ical
