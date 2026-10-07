#pragma once

// Copyright 2026, Mikhail Grushinskiy
//
// Bosch BMM150 trim-register compensation (host-testable).
//
// The BMM150 data registers hold uncompensated Hall readings. Bosch's
// datasheet (section 4.3) and BMM150 Sensor API (bmm150.c, floating-point
// compensate_x/y/z) turn them into microtesla using per-chip factory trim
// registers and RHALL, the Hall-plate resistance measured with every
// conversion. RHALL follows the sensor temperature, so this is the sensor's
// own temperature compensation of sensitivity and of the Z offset.
// M5Unified 0.2.13 reads the data registers but discards RHALL and never
// reads the trim registers.

#include <cstdint>

namespace atoms3r_ical {

// Trim registers 0x5D..0x71, as read from the sensor.
struct Bmm150Trim {
  static constexpr uint8_t FIRST_REG = 0x5D, LAST_REG = 0x71;
  static constexpr int COUNT = LAST_REG - FIRST_REG + 1;

  int8_t dig_x1 = 0, dig_y1 = 0, dig_x2 = 0, dig_y2 = 0, dig_xy2 = 0;
  uint8_t dig_xy1 = 0;
  int16_t dig_z2 = 0, dig_z3 = 0, dig_z4 = 0;
  uint16_t dig_z1 = 0, dig_xyz1 = 0;
  bool valid = false;

  // regs[i] holds register FIRST_REG + i.
  static Bmm150Trim fromRegisters(const uint8_t regs[COUNT]) {
    auto r = [&](uint8_t reg) { return regs[reg - FIRST_REG]; };
    auto u16 = [&](uint8_t lsb) { return uint16_t(uint16_t(r(lsb + 1)) << 8 | r(lsb)); };
    Bmm150Trim t;
    t.dig_x1 = int8_t(r(0x5D));
    t.dig_y1 = int8_t(r(0x5E));
    t.dig_z4 = int16_t(u16(0x62));
    t.dig_x2 = int8_t(r(0x64));
    t.dig_y2 = int8_t(r(0x65));
    t.dig_z2 = int16_t(u16(0x68));
    t.dig_z1 = u16(0x6A);
    t.dig_xyz1 = uint16_t(u16(0x6C) & 0x7FFF);
    t.dig_z3 = int16_t(u16(0x6E));
    t.dig_xy2 = int8_t(r(0x70));
    t.dig_xy1 = r(0x71);
    // Bosch divides by these; an all-zero or erased trim block is unusable.
    t.valid = t.dig_z1 != 0 && t.dig_z2 != 0 && t.dig_xyz1 != 0;
    return t;
  }
};

// Raw counts from the 8 data bytes 0x42..0x49 (X, Y, Z, RHALL, LSB first).
struct Bmm150Raw {
  int16_t x = 0, y = 0, z = 0;
  uint16_t rhall = 0;

  static Bmm150Raw fromData(const uint8_t d[8]) {
    auto s16 = [&](int i) { return int16_t(uint16_t(uint16_t(d[i + 1]) << 8 | d[i])); };
    Bmm150Raw r;
    r.x = int16_t(s16(0) >> 3);  // 13-bit signed
    r.y = int16_t(s16(2) >> 3);
    r.z = int16_t(s16(4) >> 1);  // 15-bit signed
    r.rhall = uint16_t((uint16_t(d[7]) << 8 | d[6]) >> 2);  // 14-bit unsigned
    return r;
  }
};

// Compensated field in microtesla, BMM150 sensor axes. Returns false on ADC
// overflow, missing RHALL or unusable trim (the sample must not be used).
inline bool bmm150Compensate(const Bmm150Raw& raw, const Bmm150Trim& t, float out_uT[3]) {
  constexpr int16_t OVERFLOW_XY = -4096, OVERFLOW_Z = -16384;
  if (!t.valid || raw.rhall == 0 || raw.x == OVERFLOW_XY || raw.y == OVERFLOW_XY || raw.z == OVERFLOW_Z)
    return false;
  const float rhall = float(raw.rhall);
  // X/Y: common RHALL-dependent sensitivity (Bosch compensate_x/compensate_y).
  const float h = float(t.dig_xyz1) * 16384.0f / rhall - 16384.0f;
  const float common = float(t.dig_xy2) * (h * h / 268435456.0f) + h * float(t.dig_xy1) / 16384.0f + 256.0f;
  out_uT[0] = (float(raw.x) * common * (float(t.dig_x2) + 160.0f) / 8192.0f + float(t.dig_x1) * 8.0f) / 16.0f;
  out_uT[1] = (float(raw.y) * common * (float(t.dig_y2) + 160.0f) / 8192.0f + float(t.dig_y1) * 8.0f) / 16.0f;
  // Z: RHALL-dependent offset (dig_z3) and sensitivity (dig_z1) (Bosch compensate_z).
  const float z0 = float(raw.z) - float(t.dig_z4);
  const float z2 = float(t.dig_z3) * (rhall - float(t.dig_xyz1));
  const float z4 = float(t.dig_z2) + float(t.dig_z1) * rhall / 32768.0f;
  out_uT[2] = ((z0 * 131072.0f - z2) / (z4 * 4.0f)) / 16.0f;
  return true;
}

// AtomS3R axis mapping of the compensated BMM150 field into the project body
// frame. Matches M5Unified's AtomS3R magnetometer remap (invert X and Z)
// followed by map_sensor_xyz_to_body_ned_ (sy, sx, -sz): body = (y, -x, z).
inline void bmm150SensorToAtomS3RBody(const float s[3], float body[3]) {
  body[0] = s[1];
  body[1] = -s[0];
  body[2] = s[2];
}

}  // namespace atoms3r_ical
