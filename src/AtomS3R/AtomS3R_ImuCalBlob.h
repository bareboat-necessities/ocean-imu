#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R calibration persistence and runtime application, without Arduino
  dependencies (host-testable). AtomS3R_ImuCal.h adds the NVS backend, the M5
  sensor mapping and the serial printers.

  Blob versions
    v2  original layout (key "blob_m5", older "blob"); loaded and migrated.
    v3  v2 fields in the same order, plus accelerometer qualification metadata,
        the runtime temperature clamp and the sensor identity (key "blob_m5v3").

  The accelerometer metadata describes one coefficient set; accel_coeff_crc
  binds it to that set. A thermal slope is only carried into a later session
  when the metadata is bound and the sensor identity matches.
*/

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "imu_calibrate/AccelCalFit.h"

namespace atoms3r_ical {

using Vector3f = Eigen::Matrix<float,3,1>;
using Matrix3f = Eigen::Matrix<float,3,3>;

// Config/constants (shared)
struct ImuCalCfg {
  static constexpr float g_std   = 9.80665f;
  static constexpr float DEG2RAD = 3.14159265358979323846f / 180.0f;
};

// Blob + CRC utilities
static constexpr uint32_t IMU_CAL_MAGIC = 0x434C554D; // 'MULC'
static constexpr uint8_t  IMU_CAL_MODE_M5_IMU_API = 1;

// Legacy v2 layout (read and migrated; never written).
struct ImuCalBlobV2 {
  static constexpr uint32_t IMU_CAL_MAGIC   = atoms3r_ical::IMU_CAL_MAGIC;
  static constexpr uint16_t IMU_CAL_VERSION = 2;
  static constexpr uint8_t  IMU_CAL_MODE_M5_IMU_API = atoms3r_ical::IMU_CAL_MODE_M5_IMU_API;

  uint32_t magic = IMU_CAL_MAGIC;
  uint16_t version = IMU_CAL_VERSION;
  uint16_t size_bytes = sizeof(ImuCalBlobV2);
  uint8_t  build_mode = 0;

  uint8_t  accel_ok = 0;
  float    accel_g = ImuCalCfg::g_std;
  float    accel_S[9]{};
  float    accel_T0 = 25.0f;
  float    accel_b0[3]{};
  float    accel_k[3]{};
  float    accel_rms_mag = 0.0f;

  uint8_t  gyro_ok = 0;
  float    gyro_T0 = 25.0f;
  float    gyro_b0[3]{};
  float    gyro_k[3]{};

  uint8_t  mag_ok = 0;
  float    mag_A[9]{};
  float    mag_b[3]{};
  float    mag_field_uT = 0.0f;
  float    mag_rms = 0.0f;

  uint32_t crc = 0;
};

enum class AccelFitMethod : uint8_t { LEGACY_ELLIPSOID = 0, FULL_MATRIX = 1 };

struct ImuCalBlobV3 {
  static constexpr uint32_t IMU_CAL_MAGIC   = atoms3r_ical::IMU_CAL_MAGIC;
  static constexpr uint16_t IMU_CAL_VERSION = 3;
  static constexpr uint8_t  IMU_CAL_MODE_M5_IMU_API = atoms3r_ical::IMU_CAL_MODE_M5_IMU_API;

  uint32_t magic = IMU_CAL_MAGIC;
  uint16_t version = IMU_CAL_VERSION;
  uint16_t size_bytes = sizeof(ImuCalBlobV3);
  uint8_t  build_mode = 0;

  uint8_t  accel_ok = 0;
  float    accel_g = ImuCalCfg::g_std;
  float    accel_S[9]{};            // a_cal = S*(a_raw - b0 - k*(clamp(T) - T0)), row-major
  float    accel_T0 = 25.0f;
  float    accel_b0[3]{};
  float    accel_k[3]{};
  float    accel_rms_mag = 0.0f;

  uint8_t  gyro_ok = 0;
  float    gyro_T0 = 25.0f;
  float    gyro_b0[3]{};
  float    gyro_k[3]{};

  uint8_t  mag_ok = 0;
  float    mag_A[9]{};
  float    mag_b[3]{};
  float    mag_field_uT = 0.0f;
  float    mag_rms = 0.0f;

  // ---- v3: accelerometer runtime clamp (part of the coefficient set) ----
  float    accel_T_lo = -1000.0f;
  float    accel_T_hi = 1000.0f;

  // ---- v3: accelerometer qualification (bound by accel_coeff_crc) ----
  uint8_t  accel_fit_method = 0;      // AccelFitMethod
  uint8_t  accel_thermal = 0;         // imu_cal::AccelThermal
  uint8_t  accel_thermal_reason = 0;  // imu_cal::AccelThermalReason
  uint8_t  accel_n_holds = 0;
  uint16_t accel_n_blocks = 0;
  uint16_t accel_capture_s = 0;       // total qualified hold time
  float    accel_k_temp_lo = 0.0f;    // temperature evidence range of accel_k
  float    accel_k_temp_hi = 0.0f;
  float    accel_cal_temp_lo = 0.0f;  // hold temperatures seen in this capture
  float    accel_cal_temp_hi = 0.0f;
  float    accel_bias_sigma[3] = {-1.0f, -1.0f, -1.0f};   // m/s^2, 1 sigma
  float    accel_cross_sigma[3] = {-1.0f, -1.0f, -1.0f};  // xy, xz, yz
  float    accel_k_sigma[3] = {-1.0f, -1.0f, -1.0f};      // this session's slope fit
  float    accel_sigma_obs = -1.0f;
  float    accel_cv_rms = -1.0f, accel_cv_max = -1.0f;
  float    accel_verify_rms = -1.0f, accel_verify_max = -1.0f;
  uint32_t accel_coeff_crc = 0;

  // ---- v3: identity of the sensor the coefficients belong to ----
  uint32_t sensor_id_lo = 0;          // ESP32 eFuse MAC
  uint32_t sensor_id_hi = 0;
  uint8_t  imu_type = 0;              // M5Unified imu_t
  uint8_t  reserved[3]{};

  uint32_t crc = 0;
};
static constexpr size_t IMU_CAL_CRC_LEN_V2 = offsetof(ImuCalBlobV2, crc);
static constexpr size_t IMU_CAL_CRC_LEN_V3 = offsetof(ImuCalBlobV3, crc);

// Current layout.
using ImuCalBlob = ImuCalBlobV3;

static inline uint32_t crc32_ieee_(const uint8_t* data, size_t n, uint32_t crc = 0xFFFFFFFFu) {
  for (size_t i = 0; i < n; ++i) {
    crc ^= (uint32_t)data[i];
    for (int k = 0; k < 8; ++k) {
      uint32_t mask = -(crc & 1u);
      crc = (crc >> 1) ^ (0xEDB88320u & mask);
    }
  }
  return crc;
}

static inline Matrix3f mat_from_rowmajor9_(const float a[9]) {
  Matrix3f M;
  M(0,0)=a[0]; M(0,1)=a[1]; M(0,2)=a[2];
  M(1,0)=a[3]; M(1,1)=a[4]; M(1,2)=a[5];
  M(2,0)=a[6]; M(2,1)=a[7]; M(2,2)=a[8];
  return M;
}

static inline void mat_to_rowmajor9_(const Matrix3f& M, float a[9]) {
  a[0]=M(0,0); a[1]=M(0,1); a[2]=M(0,2);
  a[3]=M(1,0); a[4]=M(1,1); a[5]=M(1,2);
  a[6]=M(2,0); a[7]=M(2,1); a[8]=M(2,2);
}

template <typename Blob>
static inline uint32_t computeBlobCrcT_(const Blob& in, size_t len) {
  uint8_t tmp[sizeof(Blob)];
  memcpy(tmp, &in, sizeof(Blob));
  return ~crc32_ieee_(tmp, len);
}

static inline uint32_t computeBlobCrc(const ImuCalBlobV2& in) { return computeBlobCrcT_(in, IMU_CAL_CRC_LEN_V2); }
static inline uint32_t computeBlobCrc(const ImuCalBlobV3& in) { return computeBlobCrcT_(in, IMU_CAL_CRC_LEN_V3); }

// CRC of the accelerometer coefficient set (everything the runtime applies).
static inline uint32_t accelCoeffCrc(const ImuCalBlobV3& b) {
  uint32_t c = 0xFFFFFFFFu;
  c = crc32_ieee_((const uint8_t*)&b.accel_ok, sizeof(b.accel_ok), c);
  c = crc32_ieee_((const uint8_t*)&b.accel_g, sizeof(b.accel_g), c);
  c = crc32_ieee_((const uint8_t*)b.accel_S, sizeof(b.accel_S), c);
  c = crc32_ieee_((const uint8_t*)&b.accel_T0, sizeof(b.accel_T0), c);
  c = crc32_ieee_((const uint8_t*)b.accel_b0, sizeof(b.accel_b0), c);
  c = crc32_ieee_((const uint8_t*)b.accel_k, sizeof(b.accel_k), c);
  c = crc32_ieee_((const uint8_t*)&b.accel_T_lo, sizeof(b.accel_T_lo), c);
  c = crc32_ieee_((const uint8_t*)&b.accel_T_hi, sizeof(b.accel_T_hi), c);
  return ~c;
}

static inline bool accelMetaBound(const ImuCalBlobV3& b) { return b.accel_coeff_crc == accelCoeffCrc(b); }

static inline bool allFinite_(const float* a, int n) {
  for (int i = 0; i < n; ++i) if (!isfinite(a[i])) return false;
  return true;
}

static inline bool validateBlob(const ImuCalBlobV2& b) {
  if (b.magic != ImuCalBlobV2::IMU_CAL_MAGIC) return false;
  if (b.version != ImuCalBlobV2::IMU_CAL_VERSION) return false;
  if (b.size_bytes != sizeof(ImuCalBlobV2)) return false;
  if (b.build_mode != IMU_CAL_MODE_M5_IMU_API) return false;
  return computeBlobCrc(b) == b.crc;
}

static inline bool validateBlob(const ImuCalBlobV3& b) {
  if (b.magic != ImuCalBlobV3::IMU_CAL_MAGIC) return false;
  if (b.version != ImuCalBlobV3::IMU_CAL_VERSION) return false;
  if (b.size_bytes != sizeof(ImuCalBlobV3)) return false;
  if (b.build_mode != IMU_CAL_MODE_M5_IMU_API) return false;
  if (computeBlobCrc(b) != b.crc) return false;
  // Coefficients the runtime will apply must be finite.
  if (b.accel_ok && !(allFinite_(b.accel_S, 9) && allFinite_(b.accel_b0, 3) && allFinite_(b.accel_k, 3) &&
                      isfinite(b.accel_T0) && isfinite(b.accel_T_lo) && isfinite(b.accel_T_hi) &&
                      b.accel_T_lo <= b.accel_T_hi)) return false;
  return true;
}

// v2 -> v3: identical coefficients and runtime behaviour (no clamp), legacy
// provenance. The metadata is bound so the migrated set is self-consistent,
// but LEGACY thermal status is never carried into a new session.
static inline void migrateV2(const ImuCalBlobV2& v2, ImuCalBlobV3& out) {
  ImuCalBlobV3 b;
  memset((void*)&b, 0, sizeof(b));
  b.magic = IMU_CAL_MAGIC;
  b.version = ImuCalBlobV3::IMU_CAL_VERSION;
  b.size_bytes = sizeof(ImuCalBlobV3);
  b.build_mode = v2.build_mode;
  b.accel_ok = v2.accel_ok; b.accel_g = v2.accel_g;
  memcpy(b.accel_S, v2.accel_S, sizeof(b.accel_S));
  b.accel_T0 = v2.accel_T0;
  memcpy(b.accel_b0, v2.accel_b0, sizeof(b.accel_b0));
  memcpy(b.accel_k, v2.accel_k, sizeof(b.accel_k));
  b.accel_rms_mag = v2.accel_rms_mag;
  b.gyro_ok = v2.gyro_ok; b.gyro_T0 = v2.gyro_T0;
  memcpy(b.gyro_b0, v2.gyro_b0, sizeof(b.gyro_b0));
  memcpy(b.gyro_k, v2.gyro_k, sizeof(b.gyro_k));
  b.mag_ok = v2.mag_ok;
  memcpy(b.mag_A, v2.mag_A, sizeof(b.mag_A));
  memcpy(b.mag_b, v2.mag_b, sizeof(b.mag_b));
  b.mag_field_uT = v2.mag_field_uT; b.mag_rms = v2.mag_rms;
  b.accel_T_lo = -1000.0f; b.accel_T_hi = 1000.0f;
  b.accel_fit_method = (uint8_t)AccelFitMethod::LEGACY_ELLIPSOID;
  b.accel_thermal = (uint8_t)imu_cal::AccelThermal::LEGACY;
  for (int j = 0; j < 3; ++j) { b.accel_bias_sigma[j] = b.accel_cross_sigma[j] = b.accel_k_sigma[j] = -1.0f; }
  b.accel_sigma_obs = b.accel_cv_rms = b.accel_cv_max = b.accel_verify_rms = b.accel_verify_max = -1.0f;
  b.accel_coeff_crc = accelCoeffCrc(b);
  b.crc = computeBlobCrc(b);
  out = b;
}

// Fills the accelerometer coefficient set and its metadata from a full fit,
// then binds the metadata. Gyro/mag fields are left untouched.
static inline void fillAccelFromFit(ImuCalBlobV3& b, const imu_cal::AccelFullFitResult& r,
                                    const imu_cal::AccelCalibration<float>& fc, uint32_t capture_s)
{
  b.accel_ok = fc.ok ? 1 : 0;
  b.accel_g = fc.g;
  mat_to_rowmajor9_(fc.S, b.accel_S);
  b.accel_T0 = fc.biasT.T0;
  for (int j = 0; j < 3; ++j) { b.accel_b0[j] = fc.biasT.b0(j); b.accel_k[j] = fc.biasT.k(j); }
  b.accel_T_lo = fc.biasT.T_lo;
  b.accel_T_hi = fc.biasT.T_hi;
  b.accel_rms_mag = fc.rms_mag;

  b.accel_fit_method = (uint8_t)AccelFitMethod::FULL_MATRIX;
  b.accel_thermal = (uint8_t)r.thermal;
  b.accel_thermal_reason = (uint8_t)r.thermal_reason;
  b.accel_n_holds = (uint8_t)r.n_fit_holds;
  b.accel_n_blocks = (uint16_t)(r.n_fit_blocks + r.n_verify_blocks);
  b.accel_capture_s = (uint16_t)(capture_s > 65535u ? 65535u : capture_s);
  b.accel_k_temp_lo = (float)r.k_temp_lo;
  b.accel_k_temp_hi = (float)r.k_temp_hi;
  b.accel_cal_temp_lo = (float)r.cal_temp_lo;
  b.accel_cal_temp_hi = (float)r.cal_temp_hi;
  for (int j = 0; j < 3; ++j) {
    b.accel_bias_sigma[j] = (float)r.bias_sigma(j);
    b.accel_cross_sigma[j] = (float)r.cross_sigma(j);
    b.accel_k_sigma[j] = (float)r.k_sigma(j);
  }
  b.accel_sigma_obs = (float)r.sigma_obs;
  b.accel_cv_rms = (float)r.cv_rms; b.accel_cv_max = (float)r.cv_max;
  b.accel_verify_rms = (float)r.ver_rms; b.accel_verify_max = (float)r.ver_max;
  b.accel_coeff_crc = accelCoeffCrc(b);
}

// Candidate of an accelerometer-only recalibration: the previous blob with
// its gyro and magnetometer fields carried unchanged and a new, bound
// accelerometer set and sensor identity.
static inline ImuCalBlobV3 accelOnlyCandidate(const ImuCalBlobV3& prev, const imu_cal::AccelFullFitResult& r,
                                              const imu_cal::AccelCalibration<float>& fc, uint32_t capture_s,
                                              uint32_t id_lo, uint32_t id_hi, uint8_t imu_type)
{
  ImuCalBlobV3 b;
  memcpy((void*)&b, &prev, sizeof(b));
  fillAccelFromFit(b, r, fc, capture_s);
  b.sensor_id_lo = id_lo;
  b.sensor_id_hi = id_hi;
  b.imu_type = imu_type;
  return b;
}

// Thermal slope usable in a new session of `sensor_id`/`imu_type`: only a
// slope this firmware validated (LEARNED, or PRESERVED from one), with bound
// metadata, for the same sensor and IMU. Bias offsets are never carried: they
// can change at every power cycle and are re-measured each session.
static inline imu_cal::AccelThermalPrior accelThermalPriorFrom(const ImuCalBlobV3& b, uint32_t id_lo,
                                                                uint32_t id_hi, uint8_t imu_type,
                                                                double max_k_abs = 0.02)
{
  imu_cal::AccelThermalPrior p;
  if (!validateBlob(b) || !b.accel_ok) return p;
  if (b.accel_fit_method != (uint8_t)AccelFitMethod::FULL_MATRIX) return p;
  const auto th = (imu_cal::AccelThermal)b.accel_thermal;
  if (th != imu_cal::AccelThermal::LEARNED && th != imu_cal::AccelThermal::PRESERVED) return p;
  if (!accelMetaBound(b)) return p;
  if ((id_lo == 0 && id_hi == 0) || b.sensor_id_lo != id_lo || b.sensor_id_hi != id_hi) return p;
  if (b.imu_type != imu_type) return p;
  for (int j = 0; j < 3; ++j) {
    if (!isfinite(b.accel_k[j]) || fabs((double)b.accel_k[j]) > max_k_abs) return p;
    p.k[j] = b.accel_k[j];
  }
  p.k_temp_lo = b.accel_k_temp_lo; p.k_temp_hi = b.accel_k_temp_hi;
  p.clamp_lo = b.accel_T_lo; p.clamp_hi = b.accel_T_hi;
  p.valid = true;
  return p;
}

// Runtime calibration objects: applied once, upstream of every estimator.
struct RuntimeCals {
  imu_cal::AccelCalibration<float> acc{};
  imu_cal::GyroCalibration<float>  gyr{};
  imu_cal::MagCalibration<float>   mag{};

  void rebuildFromBlob(const ImuCalBlobV3& b) {
    acc.ok = (b.accel_ok != 0);
    acc.g  = b.accel_g;
    acc.S  = mat_from_rowmajor9_(b.accel_S);
    acc.biasT.ok = acc.ok;
    acc.biasT.T0 = b.accel_T0;
    acc.biasT.b0 = Vector3f(b.accel_b0[0], b.accel_b0[1], b.accel_b0[2]);
    acc.biasT.k  = Vector3f(b.accel_k[0],  b.accel_k[1],  b.accel_k[2]);
    acc.biasT.T_lo = b.accel_T_lo;
    acc.biasT.T_hi = b.accel_T_hi;
    acc.rms_mag  = b.accel_rms_mag;

    gyr.ok = (b.gyro_ok != 0);
    gyr.S  = Matrix3f::Identity();
    gyr.biasT.ok = gyr.ok;
    gyr.biasT.T0 = b.gyro_T0;
    gyr.biasT.b0 = Vector3f(b.gyro_b0[0], b.gyro_b0[1], b.gyro_b0[2]);
    gyr.biasT.k  = Vector3f(b.gyro_k[0],  b.gyro_k[1],  b.gyro_k[2]);

    mag.ok = (b.mag_ok != 0);
    mag.A  = mat_from_rowmajor9_(b.mag_A);
    mag.b  = Vector3f(b.mag_b[0], b.mag_b[1], b.mag_b[2]);
    mag.field_uT = b.mag_field_uT;
    mag.rms      = b.mag_rms;
  }

  Vector3f applyAccel(const Vector3f& a_raw, float tempC) const { return acc.ok ? acc.apply(a_raw, tempC) : a_raw; }
  Vector3f applyGyro (const Vector3f& w_raw, float tempC) const { return gyr.ok ? gyr.apply(w_raw, tempC) : w_raw; }
  Vector3f applyMag  (const Vector3f& m_raw) const {
    if (!mag.ok) return m_raw;
    const Vector3f m_cal = mag.apply(m_raw);
    return m_cal.allFinite() ? m_cal : m_raw;
  }
};

// Key-value persistence over any backend with the Preferences byte API:
//   size_t getBytesLength(const char*), size_t getBytes(const char*, void*, size_t),
//   size_t putBytes(const char*, const void*, size_t), bool remove(const char*)
template <class KV>
class ImuCalStoreT {
public:
  static constexpr const char* kKeyV3 = "blob_m5v3";
  static constexpr const char* kKeyM5ImuApi = "blob_m5";   // v2
  static constexpr const char* kKeyLegacy = "blob";        // v2, older sketches

  KV kv;

  // Current v3 blob, else a migrated v2 blob.
  bool load(ImuCalBlobV3& out) {
    if (loadV3_(out)) return true;
    ImuCalBlobV2 v2;
    if (loadV2_(kKeyM5ImuApi, v2) || loadV2_(kKeyLegacy, v2)) {
      migrateV2(v2, out);
      return true;
    }
    return false;
  }

  // Writes the v3 key only; true when the bytes were accepted.
  bool save(const ImuCalBlobV3& in) {
    ImuCalBlobV3 tmp = sealed_(in);
    return kv.putBytes(kKeyV3, &tmp, sizeof(tmp)) == sizeof(tmp);
  }

  // Writes the candidate and reads back the v3 key only. Succeeds only when the
  // stored bytes validate and equal the sealed candidate, so an older blob can
  // never pass as the new one. Superseded v2 keys are then removed so a later
  // v3 failure cannot silently fall back to them.
  bool saveVerified(const ImuCalBlobV3& in, ImuCalBlobV3& readback) {
    const ImuCalBlobV3 cand = sealed_(in);
    if (kv.putBytes(kKeyV3, &cand, sizeof(cand)) != sizeof(cand)) return false;
    ImuCalBlobV3 rb;
    if (!loadV3_(rb)) return false;
    if (!sameBytes(rb, cand)) return false;
    kv.remove(kKeyM5ImuApi);
    kv.remove(kKeyLegacy);
    readback = rb;
    return true;
  }

  void erase() {
    kv.remove(kKeyV3);
    kv.remove(kKeyM5ImuApi);
    kv.remove(kKeyLegacy);
  }

  // Byte-for-byte equality of two stored blobs (padding included: both sides
  // are byte copies of what was written).
  static bool sameBytes(const ImuCalBlobV3& a, const ImuCalBlobV3& b) {
    uint8_t ba[sizeof(ImuCalBlobV3)], bb[sizeof(ImuCalBlobV3)];
    memcpy(ba, &a, sizeof(ba));
    memcpy(bb, &b, sizeof(bb));
    return memcmp(ba, bb, sizeof(ba)) == 0;
  }

  static ImuCalBlobV3 sealed_(const ImuCalBlobV3& in) {
    ImuCalBlobV3 tmp;
    memcpy(&tmp, &in, sizeof(tmp));
    tmp.magic = ImuCalBlobV3::IMU_CAL_MAGIC;
    tmp.version = ImuCalBlobV3::IMU_CAL_VERSION;
    tmp.size_bytes = sizeof(ImuCalBlobV3);
    tmp.build_mode = IMU_CAL_MODE_M5_IMU_API;
    tmp.crc = 0;
    tmp.crc = computeBlobCrc(tmp);
    return tmp;
  }

private:
  bool loadV3_(ImuCalBlobV3& out) {
    if (kv.getBytesLength(kKeyV3) != sizeof(ImuCalBlobV3)) return false;
    ImuCalBlobV3 tmp;
    if (kv.getBytes(kKeyV3, &tmp, sizeof(tmp)) != sizeof(tmp)) return false;
    if (!validateBlob(tmp)) return false;
    out = tmp;
    return true;
  }

  bool loadV2_(const char* key, ImuCalBlobV2& out) {
    if (kv.getBytesLength(key) != sizeof(ImuCalBlobV2)) return false;
    ImuCalBlobV2 tmp;
    if (kv.getBytes(key, &tmp, sizeof(tmp)) != sizeof(tmp)) return false;
    if (!validateBlob(tmp)) return false;
    out = tmp;
    return true;
  }
};

}  // namespace atoms3r_ical
