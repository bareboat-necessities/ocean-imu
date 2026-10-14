#pragma once

// Copyright 2026, Mikhail Grushinskiy
//
// BMI270 accelerometer digital filter (host-testable; Device is the BMI270,
// M5 I2C_Device-style API: readRegister(reg, buf, len) and
// writeRegister8(reg, value)).
//
// M5Unified leaves ACC_CONF (0x40) at performance mode, "normal" bandwidth.
// This selects OSR4 in performance mode: the 1.6 kHz pre-filtered stream is
// filtered as four-fold oversampled, -3 dB near ODR/10 instead of ODR/2.5
// (about 10 Hz instead of 40 Hz at M5Unified's 100 Hz). The output rate is
// kept, so polling and timestamps are unchanged; group delay grows by a few
// milliseconds. Performance mode stays on: the power-mode "avg4" setting
// duty-cycles the sensor and is noisier.
//
// ACC_CONF is part of the magnetic timing configuration hash
// (magTimingConfiguration_), so apply this before that hash is taken.

#include <cstdint>
#include <cstddef>

namespace atoms3r_ical {

struct Bmi270AccFilter {
  static constexpr uint8_t ACC_CONF = 0x40;
  static constexpr uint8_t PERF = 0x80, BWP_MASK = 0x70, ODR_MASK = 0x0F;
  static constexpr uint8_t BWP_OSR4 = 0x00;  // bits 6:4 = 0 in performance mode

  static constexpr uint8_t osr4(uint8_t acc_conf) {
    return uint8_t(PERF | BWP_OSR4 | (acc_conf & ODR_MASK));
  }

  // Writes OSR4 at the current output rate and reads it back.
  template <class Device>
  static bool apply(Device* dev, uint8_t& before, uint8_t& after) {
    before = after = 0;
    if (!dev || !dev->readRegister(ACC_CONF, &before, 1)) return false;
    const uint8_t want = osr4(before);
    if (!dev->writeRegister8(ACC_CONF, want)) return false;
    return dev->readRegister(ACC_CONF, &after, 1) && after == want;
  }
};

} // namespace atoms3r_ical
