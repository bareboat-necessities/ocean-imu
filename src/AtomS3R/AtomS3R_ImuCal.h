#pragma once

/*

  Copyright 2026, Mikhail Grushinskiy

  AtomS3R reusable IMU calibration plumbing (NVS store, axis mapping, serial
  printers and M5Unified-calibration clearing). Blob layouts, CRC, migration and
  runtime application live in AtomS3R_ImuCalBlob.h.

  This header is intentionally UI-agnostic: you can reuse it in:
    - the calibration wizard sketch
    - your "real" application firmware

  Example boot flow:

    #include <M5Unified.h>
    #include "AtomS3R/AtomS3R_ImuCal.h"

    atoms3r_ical::ImuCalStoreNvs store;
    atoms3r_ical::ImuCalBlobV3   blob;
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

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif
#include <ArduinoEigenDense.h>

// Blob layouts, CRC, v2 migration, runtime application and the generic store
// (host-testable, no Arduino dependency).
#include "AtomS3R/AtomS3R_ImuCalBlob.h"

namespace atoms3r_ical {

// Axis mapping (AtomS3R)
//
// Internal convention in this library is BODY-NED (x=north, y=east, z=down).
// End-user "nautical Z-up" is therefore z_up = -z_down.
//
// With the board lying still, screen facing up:
//   - accelerometer body Z (down) is expected near -g specific force
//   - user-facing Z-up value is the opposite sign.
// acc_body = ( ay, ax, -az ) * g
// gyr_body = ( gy, gx, -gz ) * deg2rad
// mag_body = ( my, mx, -mz ) * (1/10)
static inline Vector3f map_sensor_xyz_to_body_ned_(float sx, float sy, float sz, float scale = 1.0f) {
  return Vector3f(sy * scale, sx * scale, -sz * scale);
}

static inline Vector3f map_acc_to_body_ned_(const m5::imu_3d_t& a_g) {
  return map_sensor_xyz_to_body_ned_(a_g.x, a_g.y, a_g.z, ImuCalCfg::g_std);
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
static inline void printMat3RowMajor(Print& out, const float a[9], int prec = 9) {
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
static inline void printMatDiagOffDiagRms(Print& out, const float a[9]) {
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
static inline void printMatHeader(Print& out, const char* name, const char* meaning) {
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
static inline void printBlobSummary(Print& out, const ImuCalBlobV3& b) {
  const char* mode = "unknown";
  if (b.build_mode == IMU_CAL_MODE_M5_IMU_API) mode = "m5_imu_api";
  out.printf("  build_mode: %s\n", mode);
  out.printf("  ok: A=%d G=%d M=%d\n", (int)b.accel_ok, (int)b.gyro_ok, (int)b.mag_ok);
}

static inline void printBlobDetail(Print& out, const ImuCalBlobV3& b) {
  // ACCEL
  out.printf("  accel: g=%.6f T0=%.2f rms_mag=%.4f\n", (double)b.accel_g, (double)b.accel_T0, (double)b.accel_rms_mag);
  out.printf("    b0=[%.5f %.5f %.5f]\n", (double)b.accel_b0[0], (double)b.accel_b0[1], (double)b.accel_b0[2]);
  out.printf("    k =[%.6f %.6f %.6f] clamp T=[%.1f %.1f]\n", (double)b.accel_k[0], (double)b.accel_k[1],
             (double)b.accel_k[2], (double)b.accel_T_lo, (double)b.accel_T_hi);

  printMatHeader(out, "S", "a_cal = S*(a_raw - bias(T))");
  printMat3RowMajor(out, b.accel_S, 9);
  printMatDiagOffDiagRms(out, b.accel_S);

  const bool full = (b.accel_fit_method == (uint8_t)AccelFitMethod::FULL_MATRIX);
  out.printf("    fit=%s thermal=%s/%s meta=%s\n", full ? "full_matrix" : "legacy",
             imu_cal::accelThermalStr((imu_cal::AccelThermal)b.accel_thermal),
             imu_cal::accelThermalReasonStr((imu_cal::AccelThermalReason)b.accel_thermal_reason),
             accelMetaBound(b) ? "bound" : "UNBOUND");
  if (full) {
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

  // IMPORTANT: do not rely on a "mag valid" bit. Many builds never set it.
  out.m = map_mag_to_body_uT_(data.mag);

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

} // namespace atoms3r_ical
