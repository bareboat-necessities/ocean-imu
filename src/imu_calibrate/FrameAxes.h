#pragma once
// Copyright 2026, Mikhail Grushinskiy
#include <cstdint>

namespace imu_cal {
// Opposite directions are adjacent. This lightweight declaration is shared
// by the button UI and frame math without pulling Eigen into display helpers.
enum class FrameAxis : uint8_t { NEG_Z, POS_Z, NEG_X, POS_X, NEG_Y, POS_Y };
constexpr bool frameAxesPerpendicular(FrameAxis a,FrameAxis b) {
  return uint8_t(a)<6 && uint8_t(b)<6 && uint8_t(a)/2!=uint8_t(b)/2;
}
} // namespace imu_cal
