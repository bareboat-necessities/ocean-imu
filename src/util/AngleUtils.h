#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Scalar angle helpers shared by the AtomS3R marine INS sketches.
*/

namespace ocean_imu {
namespace ins {

// Same value and precision (double) as the Arduino RAD_TO_DEG macro, so a
// float converted here rounds exactly as `x * RAD_TO_DEG` did in the sketches.
constexpr double kRadToDeg = 57.295779513082320876798154814105;

inline float radToDeg(float rad) {
  return static_cast<float>(static_cast<double>(rad) * kRadToDeg);
}

inline float clampf(float x, float lo, float hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

// Wrap to [0, 360).
inline float wrap360Deg(float deg) {
  while (deg < 0.0f) deg += 360.0f;
  while (deg >= 360.0f) deg -= 360.0f;
  return deg;
}

// Wrap to (-180, 180].
inline float wrap180Deg(float deg) {
  while (deg <= -180.0f) deg += 360.0f;
  while (deg >   180.0f) deg -= 360.0f;
  return deg;
}

}  // namespace ins
}  // namespace ocean_imu
