/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include "AtomS3R/AtomS3R_ImuCal.h"

namespace atoms3r_ical {

bool probeMagSamples(decltype(M5.Imu)& imu) {
  int good = 0;
  const uint32_t t0 = millis();
  while ((uint32_t)(millis() - t0) < MagStartupCfg::PROBE_WINDOW_MS) {
    const uint32_t mask = imu.update();
    if (mask & m5::IMU_Class::sensor_mask_mag) {
      const Vector3f m = map_mag_to_body_uT_(imu.getImuData().mag);
      const float n = m.norm();
      if (std::isfinite(n) && n >= MagStartupCfg::MIN_NORM_uT && n <= MagStartupCfg::MAX_NORM_uT) {
        if (++good >= MagStartupCfg::PROBE_MIN_GOOD) return true;
      }
    }
    delay(2);
  }
  return false;
}

bool wakeBmm150ViaBmi270Aux() {
  constexpr uint32_t f = MagStartupCfg::I2C_FREQ;
  static constexpr uint8_t kBmiAddrs[2] = {0x68, 0x69};
  uint8_t bmi = 0;
  for (uint8_t addr : kBmiAddrs) {
    if (M5.In_I2C.readRegister8(addr, 0x00, f) == 0x24) {
      bmi = addr;
      break;
    }
  }
  if (!bmi) return false;

  M5.In_I2C.writeRegister8(bmi, 0x7C, 0x00, f);
  delay(1);
  M5.In_I2C.writeRegister8(bmi, 0x6B, 0x20, f);
  M5.In_I2C.writeRegister8(bmi, 0x7D, 0x0E, f);
  M5.In_I2C.writeRegister8(bmi, 0x4C, 0x80, f);
  M5.In_I2C.writeRegister8(bmi, 0x4B, 0x10 << 1, f);
  M5.In_I2C.writeRegister8(bmi, 0x4F, 0x01, f);
  M5.In_I2C.writeRegister8(bmi, 0x4E, 0x4B, f);
  delay(MagStartupCfg::WAKE_SETTLE_MS);
  return true;
}

bool ensureMagReady(Print& log) {
  if (!M5.Imu.isEnabled()) return false;
  if (probeMagSamples(M5.Imu)) return true;

  for (int attempt = 1; attempt <= MagStartupCfg::MAX_REINIT; ++attempt) {
    log.printf("[BOOT] magnetometer silent; reinitializing IMU (attempt %d/%d)\n",
               attempt, MagStartupCfg::MAX_REINIT);
    wakeBmm150ViaBmi270Aux();
    M5.Imu.begin(&M5.In_I2C, M5.getBoard());
    clearM5UnifiedImuCalibration();
    delay(MagStartupCfg::REINIT_SETTLE_MS);
    if (probeMagSamples(M5.Imu)) {
      log.println("[BOOT] magnetometer OK");
      return true;
    }
  }
  log.println("[BOOT] magnetometer unavailable; continuing without valid mag data");
  return false;
}

} // namespace atoms3r_ical
