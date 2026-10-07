#pragma once

/*

  Copyright 2026, Mikhail Grushinskiy

  AtomS3R reusable IMU calibration plumbing (NVS store, axis mapping, serial
  printers and M5Unified-calibration clearing). The blob layout, CRC and runtime
  application live in AtomS3R_ImuCalBlob.h.

  This header is intentionally UI-agnostic: you can reuse it in:
    - the calibration wizard sketch
    - your "real" application firmware

  Example boot flow:

    #include <M5Unified.h>
    #include "AtomS3R/AtomS3R_ImuCal.h"

    atoms3r_ical::ImuCalStoreNvs store;
    atoms3r_ical::ImuCalBlobV4   blob;
    atoms3r_ical::RuntimeCals    cals;

    void setup() {
      Serial.begin(115200);
      delay(150);
      Serial.println();

      auto cfg = M5.config();
      M5.begin(cfg);

      // 1) Clear M5Unified's own IMU calibration/offset data so it can't "stack"
      //    with our calibration (prevents two different calibrations colliding).
      atoms3r_ical::clearM5UnifiedImuCalibration();

      if (!M5.Imu.isEnabled()) {
        Serial.println("[BOOT] IMU not found");
        while (true) delay(100);
      }

      // 2) Load our calibration from NVS
      bool have = store.load(blob);

      if (have) {
        // 3) Display it at startup on Serial
        Serial.println("[BOOT] Found saved AtomS3R calibration:");
        atoms3r_ical::printBlobSummary(Serial, blob);
        atoms3r_ical::printBlobDetail(Serial, blob);

        // 4) Build runtime calibration objects for fast apply()
        cals.rebuildFromBlob(blob);
      } else {
        Serial.println("[BOOT] No saved AtomS3R calibration.");
        Serial.println("[BOOT] Starting calibration UI...");

        // 5) Start calibration UI / wizard
        // run_my_calibration_ui_and_save_blob(store);

        // After wizard saves:
        // store.load(blob); cals.rebuildFromBlob(blob);
      }
    }

    void loop() {
      atoms3r_ical::ImuSample s;
      if (!atoms3r_ical::readImuMapped(M5.Imu, s)) return;

      const auto a_cal = cals.applyAccel(s.a, s.tempC);
      const auto w_cal = cals.applyGyro (s.w, s.tempC);
      const auto m_cal = cals.applyMag  (s.m);

      // Use a_cal / w_cal / m_cal everywhere in the application
      // ...
    }

*/

#include <Arduino.h>
#include <M5Unified.h>
#include <Preferences.h>

#include <stdint.h>
#include <stddef.h>   // offsetof
#include <string.h>
#include <math.h>
#include <cmath>

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif
#include <ArduinoEigenDense.h>

// Blob layout, CRC, runtime application and the generic store
// (host-testable, no Arduino dependency).
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include "AtomS3R/AtomS3R_Bmm150AuxPreset.h"
#include "AtomS3R/AtomS3R_Bmm150Compensation.h"

namespace atoms3r_ical {

// Axis mapping and units: AtomS3R_ImuUnits.h (map_sensor_xyz_to_body_ned_).

// Nominal g -> m/s^2 with g_std (unit conversion only), then body NED.
static inline Vector3f map_acc_to_body_ned_(const m5::imu_3d_t& a_g) {
  return accel_nominal_g_to_body_ned_si_(a_g.x, a_g.y, a_g.z);
}
static inline Vector3f map_gyr_to_body_ned_(const m5::imu_3d_t& w_deg_s) {
  return map_sensor_xyz_to_body_ned_(w_deg_s.x, w_deg_s.y, w_deg_s.z, ImuCalCfg::DEG2RAD);
}
static inline Vector3f map_mag_to_body_uT_(const m5::imu_3d_t& m_raw) {
  return map_sensor_xyz_to_body_ned_(m_raw.x, m_raw.y, m_raw.z, 0.1f);
}

// Preferences byte API for ImuCalStoreT; one short open/close per operation.
class PrefsKv {
public:
  // Namespace kept stable so different sketches share saved cals.
  static constexpr const char* kNamespace = "imu_cal";

  size_t getBytesLength(const char* key) {
    Preferences prefs;
    if (!prefs.begin(kNamespace, true)) return 0;
    const size_t n = prefs.isKey(key) ? prefs.getBytesLength(key) : 0;
    prefs.end();
    return n;
  }
  size_t getBytes(const char* key, void* buf, size_t len) {
    Preferences prefs;
    if (!prefs.begin(kNamespace, true)) return 0;
    const size_t n = prefs.getBytes(key, buf, len);
    prefs.end();
    return n;
  }
  size_t putBytes(const char* key, const void* buf, size_t len) {
    Preferences prefs;
    if (!prefs.begin(kNamespace, false)) return 0;
    const size_t n = prefs.putBytes(key, buf, len);
    prefs.end();
    return n;
  }
  bool remove(const char* key) {
    Preferences prefs;
    if (!prefs.begin(kNamespace, false)) return false;
    const bool ok = prefs.isKey(key) ? prefs.remove(key) : true;
    prefs.end();
    return ok;
  }
};

// NVS store (Preferences)
class ImuCalStoreNvs : public ImuCalStoreT<PrefsKv> {};

// Pretty 3x3 print from row-major float[9].
static __attribute__((noinline)) void printMat3RowMajor(Print& out, const float a[9], int prec = 9) {
  for (int r = 0; r < 3; ++r) {
    out.print("      [");
    for (int c = 0; c < 3; ++c) {
      out.print(a[3 * r + c], prec);
      if (c < 2) out.print(", ");
    }
    out.println("]");
  }
}

// Print quick diag + off-diagonal RMS, so you can instantly see "not identity".
static __attribute__((noinline)) void printMatDiagOffDiagRms(Print& out, const float a[9]) {
  const float d0 = a[0], d1 = a[4], d2 = a[8];
  const float off2 =
      a[1]*a[1] + a[2]*a[2] +
      a[3]*a[3] + a[5]*a[5] +
      a[6]*a[6] + a[7]*a[7];
  const float off_rms = sqrtf(off2 / 6.0f);
  out.printf("      diag=[%.6f %.6f %.6f], offdiag_rms=%.6f\n", (double)d0, (double)d1, (double)d2, (double)off_rms);
}

// Identity matrix in row-major storage.
static inline const float* mat3_identity_rowmajor_() {
  static const float I[9] = {1,0,0, 0,1,0, 0,0,1};
  return I;
}

// Optional: print a matrix header line with a consistent style.
static __attribute__((noinline)) void printMatHeader(Print& out, const char* name, const char* meaning) {
  out.print("    ");
  out.print(name);
  if (meaning && meaning[0]) {
    out.print(" (");
    out.print(meaning);
    out.print(")");
  }
  out.println(":");
}

// Print helpers (startup serial)
static __attribute__((noinline)) void printBlobSummary(Print& out, const ImuCalBlobV4& b) {
  const char* mode = "unknown";
  if (b.build_mode == IMU_CAL_MODE_M5_IMU_API) mode = "m5_imu_api";
  out.printf("  build_mode: %s\n", mode);
  out.printf("  ok: A=%d G=%d M=%d\n", (int)b.accel_ok, (int)b.gyro_ok, (int)b.mag_ok);
  if (b.mag_ok && b.mag_source != activeMagSource()) {
    out.printf("  mag: fitted on %s data, firmware delivers %s: NOT applied, recalibrate the magnetometer\n",
               b.mag_source == MAG_SOURCE_BMM150_COMPENSATED ? "Bosch-compensated" : "uncompensated",
               activeMagSource() == MAG_SOURCE_BMM150_COMPENSATED ? "Bosch-compensated" : "uncompensated");
  }
  if (b.accel_ok && !accelGravityMatches(b)) {
    out.printf("  accel: fitted for g=%.7f, firmware g_cal_local=%.7f m/s^2: NOT applied, recalibrate\n",
               (double)b.accel_g, (double)ImuCalCfg::g_cal_local);
  }
}

static __attribute__((noinline)) void printBlobDetail(Print& out, const ImuCalBlobV4& b) {
  // ACCEL
  out.printf("  accel: g=%.7f (%s) T0=%.2f rms_mag=%.4f\n", (double)b.accel_g,
             accelGravityMatches(b) ? "matches g_cal_local" : "differs from g_cal_local",
             (double)b.accel_T0, (double)b.accel_rms_mag);
  out.printf("    b0=[%.5f %.5f %.5f]\n", (double)b.accel_b0[0], (double)b.accel_b0[1], (double)b.accel_b0[2]);
  out.printf("    k =[%.6f %.6f %.6f] clamp T=[%.1f %.1f]\n", (double)b.accel_k[0], (double)b.accel_k[1],
             (double)b.accel_k[2], (double)b.accel_T_lo, (double)b.accel_T_hi);

  printMatHeader(out, "S", "a_cal = S*(a_raw - bias(T))");
  printMat3RowMajor(out, b.accel_S, 9);
  printMatDiagOffDiagRms(out, b.accel_S);

  out.printf("    thermal=%s/%s meta=%s\n",
             imu_cal::accelThermalStr((imu_cal::AccelThermal)b.accel_thermal),
             imu_cal::accelThermalReasonStr((imu_cal::AccelThermalReason)b.accel_thermal_reason),
             accelMetaBound(b) ? "bound" : "UNBOUND");
  {
    out.printf("    holds=%u blocks=%u capture=%us T_seen=[%.2f %.2f] k_range=[%.2f %.2f]\n",
               (unsigned)b.accel_n_holds, (unsigned)b.accel_n_blocks, (unsigned)b.accel_capture_s,
               (double)b.accel_cal_temp_lo, (double)b.accel_cal_temp_hi,
               (double)b.accel_k_temp_lo, (double)b.accel_k_temp_hi);
    out.printf("    bias_sd=[%.4f %.4f %.4f] cross_sd=[%.5f %.5f %.5f] k_sd=[%.5f %.5f %.5f]\n",
               (double)b.accel_bias_sigma[0], (double)b.accel_bias_sigma[1], (double)b.accel_bias_sigma[2],
               (double)b.accel_cross_sigma[0], (double)b.accel_cross_sigma[1], (double)b.accel_cross_sigma[2],
               (double)b.accel_k_sigma[0], (double)b.accel_k_sigma[1], (double)b.accel_k_sigma[2]);
    out.printf("    sigma_obs=%.4f heldout_cv rms/max=%.4f/%.4f verify rms/max=%.4f/%.4f\n",
               (double)b.accel_sigma_obs, (double)b.accel_cv_rms, (double)b.accel_cv_max,
               (double)b.accel_verify_rms, (double)b.accel_verify_max);
  }

  // GYRO
  out.printf("  gyro:  T0=%.2f\n", (double)b.gyro_T0);
  out.printf("    b0=[%.6f %.6f %.6f]\n", (double)b.gyro_b0[0], (double)b.gyro_b0[1], (double)b.gyro_b0[2]);
  out.printf("    k =[%.6f %.6f %.6f]\n", (double)b.gyro_k[0], (double)b.gyro_k[1], (double)b.gyro_k[2]);

  // Gyro calibrator *always* sets S=I (stationary-only).
  printMatHeader(out, "S", "w_cal = S*(w_raw - bias(T)) ; S=I (stationary bias-only fit)");
  const float* I = mat3_identity_rowmajor_();
  printMat3RowMajor(out, I, 3);
  printMatDiagOffDiagRms(out, I);

  out.printf("    thermal=%s reason=%u bins=%u evidence T=[%.2f %.2f] clamp T=[%.2f %.2f]\n",
             b.gyro_thermal == (uint8_t)imu_cal::GyroThermal::LEARNED ? "LEARNED" : "BIAS_ONLY",
             (unsigned)b.gyro_thermal_reason, (unsigned)b.gyro_thermal_bins,
             (double)b.gyro_temp_lo, (double)b.gyro_temp_hi, (double)b.gyro_T_lo, (double)b.gyro_T_hi);
  out.printf("    k_sigma=[%.7f %.7f %.7f] information=%.3f meta=%s\n",
             (double)b.gyro_k_sigma[0], (double)b.gyro_k_sigma[1], (double)b.gyro_k_sigma[2],
             (double)b.gyro_temperature_information, gyroSetValid(b) ? "bound" : "INVALID");

  // MAG
  out.printf("  mag: field_uT=%.3f rms=%.4f\n", (double)b.mag_field_uT, (double)b.mag_rms);
  out.printf("    b=[%.3f %.3f %.3f]\n", (double)b.mag_b[0], (double)b.mag_b[1], (double)b.mag_b[2]);

  printMatHeader(out, "A", "m_cal = A*(m_raw - b)");
  printMat3RowMajor(out, b.mag_A, 9);
  printMatDiagOffDiagRms(out, b.mag_A);
}

// IMU sample + mapped read
struct ImuSample {
  Vector3f a;     // m/s^2 (mapped to body NED per AtomS3R mapping)
  Vector3f w;     // rad/s  (mapped to body NED)
  Vector3f m;     // uT     (mapped to body)
  float tempC;    // deg C
  uint32_t mask;  // M5.Imu.update() mask
  uint32_t sample_us; // micros() timestamp captured at imu.update()
};

#ifndef ATOMS3R_IMU_MASK_ACCEL
  #if defined(M5IMU_UPDATE_ACCEL)
    #define ATOMS3R_IMU_MASK_ACCEL M5IMU_UPDATE_ACCEL
  #elif defined(IMU_UPDATE_ACCEL)
    #define ATOMS3R_IMU_MASK_ACCEL IMU_UPDATE_ACCEL
  #else
    #define ATOMS3R_IMU_MASK_ACCEL (1u << 0)
  #endif
#endif

#ifndef ATOMS3R_IMU_MASK_GYRO
  #if defined(M5IMU_UPDATE_GYRO)
    #define ATOMS3R_IMU_MASK_GYRO M5IMU_UPDATE_GYRO
  #elif defined(IMU_UPDATE_GYRO)
    #define ATOMS3R_IMU_MASK_GYRO IMU_UPDATE_GYRO
  #else
    #define ATOMS3R_IMU_MASK_GYRO (1u << 1)
  #endif
#endif

static constexpr uint32_t kImuMaskAccelGyro = (ATOMS3R_IMU_MASK_ACCEL | ATOMS3R_IMU_MASK_GYRO);

// Bosch-compensated magnetometer source (see configureAtomS3RMag()).
struct AtomS3RMagSource {
  bool compensated = false;      // deliver Bosch trim/RHALL compensated uT
  Bmm150Trim trim{};
  m5::IMU_Base* bmi270 = nullptr;
  m5::imu_3d_t last_m5{};        // M5Unified's cached value of the last reading
  Vector3f last_body = Vector3f(NAN, NAN, NAN);
  bool have_last = false;
};
inline AtomS3RMagSource& atoms3rMagSource() {
  static AtomS3RMagSource source;
  return source;
}

// Magnetometer in body axes, uT. With compensation active, a new reading
// (M5Unified's cached value changed) is re-read from the BMI270 AUX data
// registers together with RHALL and compensated; otherwise the previous
// compensated value is repeated, exactly like M5Unified's cache. Without
// compensation this is M5Unified's uncompensated value.
static inline Vector3f readMagBody_(const m5::imu_3d_t& m5_mag) {
  auto& src = atoms3rMagSource();
  if (!src.compensated) return map_mag_to_body_uT_(m5_mag);
  const bool changed = !src.have_last || m5_mag.x != src.last_m5.x || m5_mag.y != src.last_m5.y ||
                       m5_mag.z != src.last_m5.z;
  if (changed) {
    uint8_t d[8];
    float s[3], b[3];
    // BMI270 AUX_DATA_0..7 mirror BMM150 0x42..0x49 (X, Y, Z, RHALL).
    if (src.bmi270 && src.bmi270->readRegister(0x04, d, sizeof(d)) &&
        bmm150Compensate(Bmm150Raw::fromData(d), src.trim, s)) {
      bmm150SensorToAtomS3RBody(s, b);
      src.last_body = Vector3f(b[0], b[1], b[2]);
      src.last_m5 = m5_mag;
      src.have_last = true;
    }
  }
  // Never mix units: until the first compensated reading this is NaN (invalid).
  return src.last_body;
}

// Reads M5.Imu, applies AtomS3R axis mapping and unit conversion, but does NOT calibrate.
static inline bool readImuMapped(decltype(M5.Imu)& imu, uint32_t update_mask, uint32_t sample_us, ImuSample& out) {
  out.sample_us = sample_us;
  out.mask = update_mask;
  if ((out.mask & kImuMaskAccelGyro) != kImuMaskAccelGyro) return false;

  const auto data = imu.getImuData();
  out.tempC = NAN;
  (void)imu.getTemp(&out.tempC);

  out.a = map_acc_to_body_ned_(data.accel);
  out.w = map_gyr_to_body_ned_(data.gyro);

  // The AtomS3R BMI270 AUX path does not provide a reliable M5Unified
  // magnetometer-update bit across deployed builds. Runtime/calibration
  // freshness is therefore determined by the existing cadence/distinct-value
  // gates rather than the update mask.
  out.m = readMagBody_(data.mag);

  return true;
}

static inline bool readImuMapped(decltype(M5.Imu)& imu, ImuSample& out) {
  const uint32_t sample_us = micros();
  const uint32_t update_mask = imu.update();
  return readImuMapped(imu, update_mask, sample_us, out);
}

// "No collisions" helper
// Clears *M5Unified's* internal calibration/offset data, so we can safely apply our own calibration
// (stored in ImuCalStoreNvs) without them stacking.
static inline void clearM5UnifiedImuCalibration() {
  // Clears runtime offsets and any stored "offset data" M5Unified may apply.
  M5.Imu.setCalibration(0, 0, 0);
  M5.Imu.clearOffsetData();
}

// AtomS3R magnetometer setup used for both calibration and runtime. Call after
// every M5Unified IMU initialisation (M5Unified soft-resets the BMM150):
//  - low-noise setting: 47 XY / 41 Z repetitions at M5Unified's 30 Hz output
//    rate instead of the 1/1 reset value (AtomS3R_Bmm150AuxPreset.h);
//  - Bosch compensation: reads the factory trim registers once, then every
//    reading is compensated with its RHALL value (AtomS3R_Bmm150Compensation.h).
// Both go through the BMI270 AUX interface with read-back. Other IMUs, or a
// failed trim read, keep M5Unified's uncompensated values. activeMagSource()
// records which units the runtime delivers, so a saved magnetometer
// calibration made in other units is not applied.
static inline bool configureAtomS3RMag(Print& log, bool compensate = true) {
  auto& src = atoms3rMagSource();
  src = AtomS3RMagSource{};
  activeMagSource() = MAG_SOURCE_M5_RAW;
  auto* imu0 = M5.Imu.getImuInstancePtr(0);
  if (M5.Imu.getType() != m5::imu_bmi270 || !imu0) return false;
  using Preset = Bmm150AuxPreset<m5::IMU_Base>;
  Preset aux(imu0, [](uint32_t ms) { delay(ms); });
  Bmm150RegState before, after;
  const bool low = aux.apply(Preset::LOW_NOISE_REP_XY, Preset::LOW_NOISE_REP_Z, before, after);
  if (low) {
    log.printf("[MAGCFG] BMM150 low-noise: %u/%u repetitions (was %u/%u), mode=0x%02X\n",
               after.nXY(), after.nZ(), before.nXY(), before.nZ(), after.mode);
  } else {
    log.printf("[MAGCFG] BMM150 low-noise setting FAILED (%s); driver setting kept\n", aux.failure());
  }
  if (!compensate) {
    log.println("[MAGCFG] BMM150 compensation off: uncompensated M5Unified values");
    return low;
  }
  uint8_t regs[Bmm150Trim::COUNT];
  const bool read = aux.readBlock(Bmm150Trim::FIRST_REG, Bmm150Trim::COUNT, regs);
  const Bmm150Trim trim = read ? Bmm150Trim::fromRegisters(regs) : Bmm150Trim{};
  if (!trim.valid) {
    log.printf("[MAGCFG] BMM150 trim %s; uncompensated M5Unified values kept\n",
               read ? "registers invalid" : aux.failure());
    return false;
  }
  src.trim = trim;
  src.bmi270 = imu0;
  src.compensated = true;
  activeMagSource() = MAG_SOURCE_BMM150_COMPENSATED;
  log.printf("[MAGCFG] BMM150 Bosch compensation on: x1=%d y1=%d x2=%d y2=%d xy1=%u xy2=%d "
             "z1=%u z2=%d z3=%d z4=%d xyz1=%u\n",
             trim.dig_x1, trim.dig_y1, trim.dig_x2, trim.dig_y2, trim.dig_xy1, trim.dig_xy2,
             trim.dig_z1, trim.dig_z2, trim.dig_z3, trim.dig_z4, trim.dig_xyz1);
  return low;
}

} // namespace atoms3r_ical
