#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R calibration persistence and runtime application, without Arduino
  dependencies (host-testable). AtomS3R_ImuCal.h adds the NVS backend, the M5
  sensor mapping and the serial printers.

  Only the current layout (version 6, key "blob_m5v6") is read or written.
  Blobs of earlier versions are not loaded; the device runs the wizard again.

  The accelerometer metadata describes one coefficient set; accel_coeff_crc
  binds it to that set. A thermal slope is only carried into a later session
  when the metadata is bound and the sensor identity matches.

  accel_g is the physical gravity the accelerometer scale was fitted against
  (ImuCalCfg::g_cal_local of the firmware that ran the wizard; bound by
  accel_coeff_crc). The runtime applies an accelerometer calibration only when
  accel_g equals this firmware's g_cal_local: a set fitted against another
  gravity (for example 9.80665) is reported and not applied, never rescaled or
  reinterpreted. Gyro and magnetometer calibrations do not depend on gravity
  and stay in use.
*/

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "imu_calibrate/AccelCalFit.h"
#include "imu_calibrate/SensorNoise.h"
#include "AtomS3R/AtomS3R_ImuUnits.h"

namespace atoms3r_ical {

using Vector3f = Eigen::Matrix<float,3,1>;
using Matrix3f = Eigen::Matrix<float,3,3>;

// Blob + CRC utilities
static constexpr uint32_t IMU_CAL_MAGIC = 0x434C554D; // 'MULC'
static constexpr uint8_t  IMU_CAL_MODE_M5_IMU_API = 1;

// Units the magnetometer calibration was fitted in (ImuCalBlobV6::mag_source).
// A calibration is applied only to data from the same source.
static constexpr uint8_t MAG_SOURCE_M5_RAW = 0;               // M5Unified uncompensated values
static constexpr uint8_t MAG_SOURCE_BMM150_COMPENSATED = 1;   // Bosch trim/RHALL compensated uT

// Source the running firmware delivers; set by configureAtomS3RMag().
inline uint8_t& activeMagSource() {
  static uint8_t source = MAG_SOURCE_M5_RAW;
  return source;
}

// Switches the running source to the one a saved calibration was fitted on
// (installed by configureAtomS3RMag(); returns true when it now matches), so
// a calibration made before compensation keeps applying until the
// magnetometer is recalibrated. Null: the source is fixed.
using MagSourceFollower = bool (*)(uint8_t source);
inline MagSourceFollower& magSourceFollower() {
  static MagSourceFollower follower = nullptr;
  return follower;
}

struct MagTimingRecord {
  uint32_t version = 1, binding = 0, configuration = 0, delay_us = 0;
  float fit_width_s = 0, residual_time_s = 0, train_rms_rad = 0, validation_rms_rad = 0, improvement = 0;
  uint32_t train_pairs = 0, validation_pairs = 0, crc = 0;
};
static_assert(sizeof(MagTimingRecord) == 48, "timing record layout");

struct ImuCalBlobV6 {
  static constexpr uint32_t IMU_CAL_MAGIC   = atoms3r_ical::IMU_CAL_MAGIC;
  static constexpr uint16_t IMU_CAL_VERSION = 6;
  static constexpr uint8_t  IMU_CAL_MODE_M5_IMU_API = atoms3r_ical::IMU_CAL_MODE_M5_IMU_API;

  uint32_t magic = IMU_CAL_MAGIC;
  uint16_t version = IMU_CAL_VERSION;
  uint16_t size_bytes = sizeof(ImuCalBlobV6);
  uint8_t  build_mode = 0;

  uint8_t  accel_ok = 0;
  uint8_t  pad_a[2]{};               // explicit padding: CRC and readback compare every byte
  float    accel_g = ImuCalCfg::g_cal_local;  // physical gravity the scale was fitted against
  float    accel_S[9]{};            // a_cal = S*(a_raw - b0 - k*(clamp(T) - T0)), row-major
  float    accel_T0 = 25.0f;
  float    accel_b0[3]{};
  float    accel_k[3]{};
  float    accel_rms_mag = 0.0f;

  uint8_t  gyro_ok = 0;
  uint8_t  pad_g[3]{};
  float    gyro_T0 = 25.0f;
  float    gyro_b0[3]{};
  float    gyro_k[3]{};

  uint8_t  mag_ok = 0;
  uint8_t  mag_source = MAG_SOURCE_M5_RAW;  // units the mag fit used (MAG_SOURCE_*)
  uint8_t  pad_m[2]{};
  float    mag_A[9]{};
  float    mag_b[3]{};
  float    mag_field_uT = 0.0f;
  float    mag_rms = 0.0f;

  // ---- accelerometer runtime clamp (part of the coefficient set) ----
  float    accel_T_lo = -1000.0f;
  float    accel_T_hi = 1000.0f;

  // ---- accelerometer qualification (bound by accel_coeff_crc) ----
  uint8_t  accel_thermal = 0;         // imu_cal::AccelThermal
  uint8_t  accel_thermal_reason = 0;  // imu_cal::AccelThermalReason
  uint8_t  accel_n_holds = 0;
  uint8_t  reserved0 = 0;
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

  // ---- identity of the sensor the coefficients belong to ----
  uint32_t sensor_id_lo = 0;          // ESP32 eFuse MAC
  uint32_t sensor_id_hi = 0;
  uint8_t  imu_type = 0;              // M5Unified imu_t
  uint8_t  reserved[3]{};

  // ---- gyro qualification and runtime clamp (bound by gyro_coeff_crc) ----
  float gyro_T_lo = 25.0f, gyro_T_hi = 25.0f;
  float gyro_temp_lo = 25.0f, gyro_temp_hi = 25.0f;
  float gyro_k_sigma[3] = {-1.0f, -1.0f, -1.0f};
  float gyro_temperature_information = 0.0f;
  uint8_t gyro_thermal = 0;        // imu_cal::GyroThermal
  uint8_t gyro_thermal_reason = 0; // imu_cal::GyroThermalReason
  uint8_t gyro_thermal_bins = 0;
  uint8_t gyro_reserved = 0;
  uint8_t gyro_matrix_ok = 0, gyro_matrix_pad[3]{};
  float gyro_S[9] = {1,0,0,0,1,0,0,0,1};
  float gyro_train_rms=0,gyro_validation_rms=0,gyro_validation_max=0,gyro_condition=0;
  uint16_t gyro_train_segments=0,gyro_validation_segments=0;
  uint32_t gyro_coeff_crc = 0;

  // Physical rotations: columns are source axes expressed in destination.
  uint8_t mag_alignment_ok=0,housing_alignment_ok=0,vessel_alignment_ok=0,frame_reserved=0;
  float R_IM[9]={1,0,0,0,1,0,0,0,1};
  float R_HI[9]={1,0,0,0,1,0,0,0,1};
  float R_VH[9]={1,0,0,0,1,0,0,0,1};
  float alignment_train_sd=0,alignment_validation_sd=0,alignment_condition=0;
  uint16_t alignment_train_poses=0,alignment_validation_poses=0;
  imu_cal::SensorNoiseModel noise{};
  MagTimingRecord mag_timing{};

  uint32_t crc = 0;
};
static constexpr size_t IMU_CAL_CRC_LEN_V6 = offsetof(ImuCalBlobV6, crc);
// Named fields cover the entire wire record (no Eigen objects or pointers).
static_assert(sizeof(ImuCalBlobV6) == 680, "ImuCalBlobV6: implicit padding");
static_assert(sizeof(ImuCalBlobV6) == IMU_CAL_CRC_LEN_V6 + 4, "Blob must have no trailing padding");

// Current layout.
using ImuCalBlob = ImuCalBlobV6;

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

static inline uint32_t computeBlobCrc(const ImuCalBlobV6& in) { return computeBlobCrcT_(in, IMU_CAL_CRC_LEN_V6); }

// CRC of the accelerometer coefficient set (everything the runtime applies).
static inline uint32_t accelCoeffCrc(const ImuCalBlobV6& b) {
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

static inline bool accelMetaBound(const ImuCalBlobV6& b) { return b.accel_coeff_crc == accelCoeffCrc(b); }

// Largest |accel_g - g_cal_local| treated as the same gravity: float rounding
// only (a height-datum mix-up is ~1e-4, 9.80665 vs local ~4e-3 m/s^2).
static constexpr float kAccelGravityMatchTol = 1.0e-5f;

// True when the accelerometer set was fitted against `g_cfg` (the physical
// gravity this firmware calibrates to and its estimators remove).
static inline bool accelGravityMatches(const ImuCalBlobV6& b, float g_cfg = ImuCalCfg::g_cal_local) {
  return isfinite(b.accel_g) && fabsf(b.accel_g - g_cfg) <= kAccelGravityMatchTol;
}

static inline bool allFinite_(const float* a, int n) {
  for (int i = 0; i < n; ++i) if (!isfinite(a[i])) return false;
  return true;
}

// Unlike an unqualified legacy coefficient, every gyro thermal model has its
// coefficients, evidence, uncertainty and permitted runtime interval bound
// together. The outer blob CRC also binds these bytes and sensor identity.
static inline uint32_t gyroCoeffCrc(const ImuCalBlobV6& b) {
  uint32_t c = 0xFFFFFFFFu;
  c = crc32_ieee_((const uint8_t*)&b.gyro_ok, sizeof(b.gyro_ok), c);
  c = crc32_ieee_((const uint8_t*)&b.gyro_T0, sizeof(b.gyro_T0), c);
  c = crc32_ieee_((const uint8_t*)b.gyro_b0, sizeof(b.gyro_b0), c);
  c = crc32_ieee_((const uint8_t*)b.gyro_k, sizeof(b.gyro_k), c);
  c = crc32_ieee_((const uint8_t*)&b.gyro_T_lo,
                 offsetof(ImuCalBlobV6, gyro_coeff_crc) - offsetof(ImuCalBlobV6, gyro_T_lo), c);
  return ~c;
}

static inline bool gyroSetValid(const ImuCalBlobV6& b) {
  using L = imu_cal::GyroThermalLimits;
  if (!b.gyro_ok || b.gyro_coeff_crc != gyroCoeffCrc(b)) return false;
  if (!allFinite_(b.gyro_S,9) || b.gyro_matrix_ok>1) return false;
  const Matrix3f S=mat_from_rowmajor9_(b.gyro_S);
  Eigen::JacobiSVD<Matrix3f> svd(S);
  if(S.determinant()<=0 || svd.singularValues().minCoeff()<.8f ||
     svd.singularValues().maxCoeff()>1.25f ||
     (S-Matrix3f::Identity()).cwiseAbs().maxCoeff()>.20f)return false;
  if(b.gyro_matrix_ok && (!isfinite(b.gyro_train_rms) || b.gyro_train_rms<0 || b.gyro_train_rms>.025f ||
      !isfinite(b.gyro_validation_rms) || b.gyro_validation_rms<0 || b.gyro_validation_rms>.025f ||
      !isfinite(b.gyro_validation_max) || b.gyro_validation_max<0 || b.gyro_validation_max>.06f ||
      !isfinite(b.gyro_condition) || b.gyro_condition<1 || b.gyro_condition>1e5f ||
      b.gyro_train_segments<12 || b.gyro_validation_segments<4))return false;
  if (!allFinite_(b.gyro_b0, 3) || !allFinite_(b.gyro_k, 3) || !allFinite_(b.gyro_k_sigma, 3) ||
      !isfinite(b.gyro_T0) || !isfinite(b.gyro_T_lo) || !isfinite(b.gyro_T_hi) ||
      !isfinite(b.gyro_temp_lo) || !isfinite(b.gyro_temp_hi) ||
      !isfinite(b.gyro_temperature_information) || b.gyro_temperature_information < 0 ||
      b.gyro_temp_lo > b.gyro_T0 || b.gyro_T0 > b.gyro_temp_hi ||
      b.gyro_T_lo > b.gyro_T0 || b.gyro_T0 > b.gyro_T_hi || b.gyro_thermal_bins > 16) return false;
  if (b.gyro_thermal == (uint8_t)imu_cal::GyroThermal::UNLEARNED) {
    return b.gyro_thermal_reason <= (uint8_t)imu_cal::GyroThermalReason::SLOPE_IMPLAUSIBLE &&
           b.gyro_k[0] == 0 && b.gyro_k[1] == 0 && b.gyro_k[2] == 0 &&
           b.gyro_T_lo == b.gyro_T0 && b.gyro_T_hi == b.gyro_T0;
  }
  if (b.gyro_thermal != (uint8_t)imu_cal::GyroThermal::LEARNED ||
      b.gyro_thermal_reason != (uint8_t)imu_cal::GyroThermalReason::QUALIFIED ||
      b.gyro_thermal_bins < L::min_bins ||
      double(b.gyro_temp_hi) - b.gyro_temp_lo < L::min_span ||
      b.gyro_temperature_information < L::min_information ||
      fabs(double(b.gyro_T_lo) - (double(b.gyro_temp_lo) - L::extrapolation_margin)) > 1e-5 ||
      fabs(double(b.gyro_T_hi) - (double(b.gyro_temp_hi) + L::extrapolation_margin)) > 1e-5) return false;
  for (int j = 0; j < 3; ++j) {
    if (b.gyro_k_sigma[j] < 0 || double(b.gyro_k_sigma[j]) > L::max_slope_sigma ||
        fabs(double(b.gyro_k[j])) + 3 * double(b.gyro_k_sigma[j]) > L::max_slope) return false;
  }
  return true;
}

static inline void fillGyroFromFit(ImuCalBlobV6& b, const imu_cal::GyroCalibration<float>& g) {
  b.gyro_ok = g.ok ? 1 : 0;
  b.gyro_T0 = g.biasT.T0;
  for (int j = 0; j < 3; ++j) {
    b.gyro_b0[j] = g.biasT.b0[j]; b.gyro_k[j] = g.biasT.k[j]; b.gyro_k_sigma[j] = g.slope_sigma[j];
  }
  b.gyro_T_lo = g.biasT.T_lo; b.gyro_T_hi = g.biasT.T_hi;
  b.gyro_temp_lo = g.temp_lo; b.gyro_temp_hi = g.temp_hi;
  b.gyro_temperature_information = g.temperature_information;
  b.gyro_thermal = (uint8_t)g.thermal;
  b.gyro_thermal_reason = (uint8_t)g.thermal_reason;
  b.gyro_thermal_bins = g.thermal_bins;
  b.gyro_reserved = 0;
  b.gyro_matrix_ok=g.matrix_ok;
  mat_to_rowmajor9_(g.S,b.gyro_S);
  b.gyro_train_rms=g.rotation_train_rms;b.gyro_validation_rms=g.rotation_validation_rms;
  b.gyro_validation_max=g.rotation_validation_max;b.gyro_condition=g.rotation_condition;
  b.gyro_train_segments=g.rotation_train_segments;b.gyro_validation_segments=g.rotation_validation_segments;
  b.gyro_coeff_crc = gyroCoeffCrc(b);
}

static inline bool rotationValid(const float r[9]) {
  if(!allFinite_(r,9))return false;
  const Matrix3f R=mat_from_rowmajor9_(r);
  return (R.transpose()*R-Matrix3f::Identity()).norm()<1e-4f && fabsf(R.determinant()-1)<1e-4f;
}
static inline bool framesValid(const ImuCalBlobV6& b) {
  return b.mag_alignment_ok<=1 && b.housing_alignment_ok<=1 && b.vessel_alignment_ok<=1 &&
    rotationValid(b.R_IM) && rotationValid(b.R_HI) && rotationValid(b.R_VH) &&
    (!b.mag_alignment_ok || (b.mag_ok && isfinite(b.alignment_train_sd) &&
      b.alignment_train_sd>=0 && b.alignment_train_sd<1.5f &&
      isfinite(b.alignment_validation_sd) && b.alignment_validation_sd>=0 && b.alignment_validation_sd<1.5f &&
      isfinite(b.alignment_condition) && b.alignment_condition>=1 && b.alignment_condition<1e5f &&
      b.alignment_train_poses>=12 && b.alignment_validation_poses>=4));
}

static inline bool magSetValid(const ImuCalBlobV6& b) {
  if (!b.mag_ok || !allFinite_(b.mag_A, 9) || !allFinite_(b.mag_b, 3) ||
      !isfinite(b.mag_field_uT) || b.mag_field_uT < 12 || b.mag_field_uT > 120 ||
      !isfinite(b.mag_rms) || b.mag_rms < 0 || b.mag_rms > imu_cal::MagFitLimits::rmsLimit(b.mag_field_uT)) return false;
  const Matrix3f A = mat_from_rowmajor9_(b.mag_A);
  const Vector3f bias(b.mag_b[0],b.mag_b[1],b.mag_b[2]);
  return imu_cal::MagGeometricFit<1>::matrixValid(A,bias);
}

inline uint32_t magTimingBinding(const ImuCalBlobV6& b) {
  uint32_t c=crc32_ieee_(reinterpret_cast<const uint8_t*>(&b.gyro_ok),
      offsetof(ImuCalBlobV6, accel_T_lo)-offsetof(ImuCalBlobV6, gyro_ok));
  const uint32_t tail[]={gyroCoeffCrc(b), b.sensor_id_lo, b.sensor_id_hi, b.imu_type};
  c=crc32_ieee_(reinterpret_cast<const uint8_t*>(tail),sizeof(tail),c);
  // Timing fitting uses rates and fields registered in I.
  c=crc32_ieee_(reinterpret_cast<const uint8_t*>(b.R_IM),sizeof(b.R_IM),c);
  return ~crc32_ieee_(&b.mag_alignment_ok,sizeof(b.mag_alignment_ok),c);
}
inline uint32_t magTimingCrc(const MagTimingRecord& r) {
  return ~crc32_ieee_(reinterpret_cast<const uint8_t*>(&r),offsetof(MagTimingRecord,crc));
}
inline bool magTimingValid(const MagTimingRecord& r) {
  return r.version==1 && r.configuration && r.crc==magTimingCrc(r) &&
      r.delay_us>=10000 && r.delay_us<=70000 &&
      std::isfinite(r.fit_width_s) && r.fit_width_s>=.002f && r.fit_width_s<=.012f &&
      std::isfinite(r.residual_time_s) && r.residual_time_s>=.025f && r.residual_time_s<=.060f &&
      std::isfinite(r.train_rms_rad) && r.train_rms_rad>=0 && r.train_rms_rad<.05f &&
      std::isfinite(r.validation_rms_rad) && r.validation_rms_rad>=0 && r.validation_rms_rad<.05f &&
      std::isfinite(r.improvement) && r.improvement>=.05f && r.improvement<=1 &&
      r.train_pairs>=50 && r.validation_pairs>=50;
}
// The timing record is optional. An absent record has every field zero (the
// layout version may be 0 or 1); a present record must be qualified and bound
// to the saved gyro/mag set. Runtime magnetic corrections never apply it.
inline bool magTimingAbsent(const MagTimingRecord& r) {
  // Bitwise zero, so a negative zero still counts as a present record.
  auto zero = [](float v) { uint32_t u; memcpy(&u, &v, sizeof(u)); return u == 0; };
  return r.binding == 0 && r.configuration == 0 && r.delay_us == 0 &&
         zero(r.fit_width_s) && zero(r.residual_time_s) && zero(r.train_rms_rad) &&
         zero(r.validation_rms_rad) && zero(r.improvement) &&
         r.train_pairs == 0 && r.validation_pairs == 0 && r.crc == 0;
}
inline bool magTimingCompatible(const ImuCalBlobV6& b, uint32_t id_lo, uint32_t id_hi,
                                uint8_t imu_type, uint8_t source, uint32_t configuration) {
  return gyroSetValid(b) && magSetValid(b) && b.sensor_id_lo==id_lo &&
      b.sensor_id_hi==id_hi && b.imu_type==imu_type && b.mag_source==source &&
      magTimingValid(b.mag_timing) && b.mag_timing.binding==magTimingBinding(b) &&
      b.mag_timing.configuration==configuration;
}
inline MagTimingRecord& activeMagTiming() { static MagTimingRecord r{}; return r; }
inline uint32_t& activeMagTimingConfiguration() { static uint32_t c=0; return c; }
using MagTimingLoader = void (*)(const ImuCalBlobV6&);
inline MagTimingLoader& magTimingLoader() { static MagTimingLoader loader=nullptr; return loader; }


static inline bool validateBlob(const ImuCalBlobV6& b) {
  if (b.magic != ImuCalBlobV6::IMU_CAL_MAGIC) return false;
  if (b.version != ImuCalBlobV6::IMU_CAL_VERSION) return false;
  if (b.size_bytes != sizeof(ImuCalBlobV6)) return false;
  if (b.build_mode != IMU_CAL_MODE_M5_IMU_API) return false;
  if (computeBlobCrc(b) != b.crc) return false;
  if (!framesValid(b) || !imu_cal::noiseModelValid(b.noise)) return false;
  if (b.gyro_ok && !gyroSetValid(b)) return false;
  if (b.mag_ok && !magSetValid(b)) return false;
  if (!magTimingAbsent(b.mag_timing) && (!b.mag_ok || !magTimingValid(b.mag_timing) ||
      b.mag_timing.binding!=magTimingBinding(b))) return false;
  // Coefficients the runtime will apply must be finite.
  if (b.accel_ok && !(allFinite_(b.accel_S, 9) && allFinite_(b.accel_b0, 3) && allFinite_(b.accel_k, 3) &&
                      isfinite(b.accel_T0) && isfinite(b.accel_T_lo) && isfinite(b.accel_T_hi) &&
                      b.accel_T_lo <= b.accel_T_hi)) return false;
  return true;
}

// Fills the accelerometer coefficient set and its metadata from a full fit,
// then binds the metadata. Gyro/mag fields are left untouched.
static inline void fillAccelFromFit(ImuCalBlobV6& b, const imu_cal::AccelFullFitResult& r,
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

// Component-update helper: retain intrinsic gyro/mag evidence, but invalidate
// their dependent registrations and noise/timing qualifications. This partial
// candidate cannot pass the full estimator boot gate.
static inline ImuCalBlobV6 accelOnlyCandidate(const ImuCalBlobV6& prev, const imu_cal::AccelFullFitResult& r,
                                              const imu_cal::AccelCalibration<float>& fc, uint32_t capture_s,
                                              uint32_t id_lo, uint32_t id_hi, uint8_t imu_type)
{
  ImuCalBlobV6 b;
  memcpy((void*)&b, &prev, sizeof(b));
  fillAccelFromFit(b, r, fc, capture_s);
  b.sensor_id_lo = id_lo;
  b.sensor_id_hi = id_hi;
  b.imu_type = imu_type;
  // A changed accelerometer model changes gravity endpoints and the noise
  // transformation. The old rotation/noise qualifications cannot be reused.
  b.gyro_matrix_ok=0;b.mag_alignment_ok=0;b.noise={};
  b.gyro_coeff_crc=gyroCoeffCrc(b);
  b.mag_timing={};
  return b;
}

// Thermal slope usable in a new session of `sensor_id`/`imu_type`: only a
// slope this firmware validated (LEARNED, or PRESERVED from one), with bound
// metadata, for the same sensor and IMU. Bias offsets are never carried: they
// can change at every power cycle and are re-measured each session.
static inline imu_cal::AccelThermalPrior accelThermalPriorFrom(const ImuCalBlobV6& b, uint32_t id_lo,
                                                                uint32_t id_hi, uint8_t imu_type,
                                                                double max_k_abs = 0.02)
{
  imu_cal::AccelThermalPrior p;
  if (!validateBlob(b) || !b.accel_ok) return p;
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
  Matrix3f R_VI=Matrix3f::Identity(),R_IM=Matrix3f::Identity();
  bool frames_ok=false,mag_alignment_ok=false,record_valid=false;
  imu_cal::SensorNoiseModel noise{};

  Vector3f noiseAccelStd() const {return Vector3f(noise.sample_std[0],noise.sample_std[1],noise.sample_std[2]);}
  Vector3f noiseGyroDensity() const {return Vector3f(noise.white_density[3],noise.white_density[4],noise.white_density[5]);}
  template<class Filter> bool configureProcessNoise(Filter& filter) const {
    if(!record_valid || !noise.valid || !imu_cal::noiseModelValid(noise))return false;
    filter.set_gyro_noise_density_rad_sqrt_s(noiseGyroDensity());
    const Vector3f walk(noise.bias_random_walk[3],noise.bias_random_walk[4],noise.bias_random_walk[5]);
    filter.set_measured_gyro_bias_walk(walk,noise.walk_mask>>3);
    return true;
  }
  // Set when the blob holds an accelerometer calibration fitted against a
  // different gravity than this firmware's g_cal_local; it is then not applied.
  bool accel_gravity_mismatch = false;
  // Set when the blob's magnetometer fit was made from a different data source
  // than this firmware delivers; it is then not applied.
  bool mag_source_mismatch = false;

  void rebuildFromBlob(const ImuCalBlobV6& b, float g_cfg = ImuCalCfg::g_cal_local) {
    record_valid=validateBlob(b);
    // Fail closed even if a caller bypassed the verified store.
    if(!record_valid) {*this=RuntimeCals{};return;}
    frames_ok=b.housing_alignment_ok && b.vessel_alignment_ok;
    mag_alignment_ok=b.mag_alignment_ok;
    R_VI=mat_from_rowmajor9_(b.R_VH)*mat_from_rowmajor9_(b.R_HI);
    R_IM=mat_from_rowmajor9_(b.R_IM);noise=b.noise;
    accel_gravity_mismatch = (b.accel_ok != 0) && !accelGravityMatches(b, g_cfg);
    acc.ok = (b.accel_ok != 0) && !accel_gravity_mismatch;
    acc.g  = b.accel_g;
    acc.S  = mat_from_rowmajor9_(b.accel_S);
    acc.biasT.ok = acc.ok;
    acc.biasT.T0 = b.accel_T0;
    acc.biasT.b0 = Vector3f(b.accel_b0[0], b.accel_b0[1], b.accel_b0[2]);
    acc.biasT.k  = Vector3f(b.accel_k[0],  b.accel_k[1],  b.accel_k[2]);
    acc.biasT.T_lo = b.accel_T_lo;
    acc.biasT.T_hi = b.accel_T_hi;
    acc.rms_mag  = b.accel_rms_mag;

    gyr = imu_cal::GyroCalibration<float>{};
    gyr.ok = gyroSetValid(b);
    gyr.S  = mat_from_rowmajor9_(b.gyro_S);
    gyr.matrix_ok=b.gyro_matrix_ok;
    gyr.biasT.ok = gyr.ok;
    gyr.biasT.T0 = b.gyro_T0;
    gyr.biasT.b0 = Vector3f(b.gyro_b0[0], b.gyro_b0[1], b.gyro_b0[2]);
    gyr.biasT.k  = Vector3f(b.gyro_k[0],  b.gyro_k[1],  b.gyro_k[2]);

    gyr.biasT.T_lo = b.gyro_T_lo; gyr.biasT.T_hi = b.gyro_T_hi;
    gyr.thermal = (imu_cal::GyroThermal)b.gyro_thermal;
    gyr.thermal_reason = (imu_cal::GyroThermalReason)b.gyro_thermal_reason;
    gyr.thermal_bins = b.gyro_thermal_bins;
    gyr.temp_lo = b.gyro_temp_lo; gyr.temp_hi = b.gyro_temp_hi;
    gyr.temperature_information = b.gyro_temperature_information;
    gyr.slope_sigma = Vector3f(b.gyro_k_sigma[0], b.gyro_k_sigma[1], b.gyro_k_sigma[2]);

    // Read the magnetometer the way the saved fit was made when the firmware
    // can; a fit made in other units (e.g. before Bosch compensation) that
    // the firmware cannot deliver is not applied.
    if (magSetValid(b) && b.mag_source != activeMagSource() && magSourceFollower())
      magSourceFollower()(b.mag_source);
    mag_source_mismatch = magSetValid(b) && b.mag_source != activeMagSource();
    mag.ok = magSetValid(b) && !mag_source_mismatch;
    mag.A  = mat_from_rowmajor9_(b.mag_A);
    mag.b  = Vector3f(b.mag_b[0], b.mag_b[1], b.mag_b[2]);
    mag.field_uT = b.mag_field_uT;
    mag.rms      = b.mag_rms;
    // Report the optional timing record; it does not change runtime behaviour.
    activeMagTiming() = {};
    if (mag.ok && gyr.ok && magTimingLoader()) magTimingLoader()(b);
  }

  Vector3f applyAccel(const Vector3f& a_raw, float tempC) const { return acc.ok ? Vector3f(R_VI*acc.apply(a_raw, tempC)) : a_raw; }
  Vector3f applyGyro (const Vector3f& w_raw, float tempC) const { return gyr.ok ? Vector3f(R_VI*gyr.apply(w_raw, tempC)) : w_raw; }
  Vector3f applyMag  (const Vector3f& m_raw) const {
    if (!mag.ok) return m_raw;
    const Vector3f m_cal = R_VI*R_IM*mag.apply(m_raw);
    return m_cal.allFinite() ? m_cal : m_raw;
  }
};

// Why the runtime cannot apply the saved calibration in full, or nullptr when
// accelerometer, gyroscope and magnetometer sets are all applied. A sketch
// that must not run an estimator on raw sensor values stops on a reason.
static inline const char* calibrationUnusableReason(bool have_blob, const RuntimeCals& rc) {
  if (!have_blob) return "no valid saved calibration";
  if (!rc.record_valid) return "invalid calibration; recalibrate";
  if (rc.accel_gravity_mismatch) return "accel fitted against another gravity";
  if (!rc.acc.ok) return "no accelerometer calibration";
  if (!rc.gyr.ok) return "no gyroscope calibration";
  if (!rc.gyr.matrix_ok) return "gyro rotation calibration required";
  if (!rc.frames_ok) return "housing/vessel alignment required";
  if (rc.mag_source_mismatch) return "mag fitted on another data source";
  if (!rc.mag.ok) return "no magnetometer calibration";
  if (!rc.mag_alignment_ok) return "magnetometer alignment required";
  return nullptr;
}

// Key-value persistence over any backend with the Preferences byte API:
//   size_t getBytesLength(const char*), size_t getBytes(const char*, void*, size_t),
//   size_t putBytes(const char*, const void*, size_t), bool remove(const char*)
template <class KV>
class ImuCalStoreT {
public:
  static constexpr const char* kKeyV6 = "blob_m5v6";

  KV kv;

  bool load(ImuCalBlobV6& out) { return loadV6_(out); }

  // Writes the v6 key only; true when the bytes were accepted.
  bool save(const ImuCalBlobV6& in) {
    if(in.magic!=ImuCalBlobV6::IMU_CAL_MAGIC || in.version!=6 || in.size_bytes!=sizeof(in))return false;
    ImuCalBlobV6 tmp = sealed_(in);
    return validateBlob(tmp) && kv.putBytes(kKeyV6, &tmp, sizeof(tmp)) == sizeof(tmp);
  }

  // Writes the candidate and reads the key back. Succeeds only when the stored
  // bytes validate and equal the sealed candidate, so an older blob can never
  // pass as the new one.
  // Recovery is verified too. Hardware can fail during rollback: callers must
  // inspect lastSaveStatus(), never claim a previous calibration was restored
  // merely because its write was attempted.
  enum class SaveStatus : uint8_t {
    OK, INVALID_CANDIDATE, PREVIOUS_UNREADABLE, PREVIOUS_RETAINED, FAILED_EMPTY, RECOVERY_FAILED
  };
  SaveStatus lastSaveStatus() const { return last_save_status_; }

  bool saveVerified(const ImuCalBlobV6& in, ImuCalBlobV6& readback) {
    if(in.magic!=ImuCalBlobV6::IMU_CAL_MAGIC || in.version!=6 || in.size_bytes!=sizeof(in)) {
      last_save_status_=SaveStatus::INVALID_CANDIDATE;return false;
    }
    const ImuCalBlobV6 cand = sealed_(in);
    if (!validateBlob(cand)) { last_save_status_ = SaveStatus::INVALID_CANDIDATE; return false; }
    ImuCalBlobV6 prev;
    bool read_failed = false;
    const bool had_prev = loadV6_(prev, &read_failed);
    // Do not confuse a failed read of a possibly valid record with absence.
    if (read_failed) { last_save_status_ = SaveStatus::PREVIOUS_UNREADABLE; return false; }
    ImuCalBlobV6 rb;
    if (kv.putBytes(kKeyV6, &cand, sizeof(cand)) == sizeof(cand) && loadV6_(rb) && sameBytes(rb, cand)) {
      readback = rb;
      last_save_status_ = SaveStatus::OK;
      return true;
    }
    if (had_prev) {
      // A dropped write may already have left the previous bytes intact.
      if (loadV6_(rb) && sameBytes(rb, prev)) {
        last_save_status_ = SaveStatus::PREVIOUS_RETAINED;
      } else {
        const bool wrote = kv.putBytes(kKeyV6, &prev, sizeof(prev)) == sizeof(prev);
        const bool restored = loadV6_(rb) && sameBytes(rb, prev);
        last_save_status_ = wrote && restored ? SaveStatus::PREVIOUS_RETAINED : SaveStatus::RECOVERY_FAILED;
      }
    } else {
      const bool removed = kv.remove(kKeyV6);
      last_save_status_ = removed && kv.getBytesLength(kKeyV6) == 0 ?
                          SaveStatus::FAILED_EMPTY : SaveStatus::RECOVERY_FAILED;
    }
    return false;
  }

  void erase() { kv.remove(kKeyV6); }

  // Byte-for-byte equality of two stored blobs (padding included: both sides
  // are byte copies of what was written).
  static bool sameBytes(const ImuCalBlobV6& a, const ImuCalBlobV6& b) {
    uint8_t ba[sizeof(ImuCalBlobV6)], bb[sizeof(ImuCalBlobV6)];
    memcpy(ba, &a, sizeof(ba));
    memcpy(bb, &b, sizeof(bb));
    return memcmp(ba, bb, sizeof(ba)) == 0;
  }

  static ImuCalBlobV6 sealed_(const ImuCalBlobV6& in) {
    ImuCalBlobV6 tmp;
    memcpy(&tmp, &in, sizeof(tmp));
    tmp.magic = ImuCalBlobV6::IMU_CAL_MAGIC;
    tmp.version = ImuCalBlobV6::IMU_CAL_VERSION;
    tmp.size_bytes = sizeof(ImuCalBlobV6);
    tmp.build_mode = IMU_CAL_MODE_M5_IMU_API;
    memset(tmp.pad_a, 0, sizeof(tmp.pad_a));
    memset(tmp.pad_g, 0, sizeof(tmp.pad_g));
    memset(tmp.pad_m, 0, sizeof(tmp.pad_m));
    tmp.crc = 0;
    tmp.crc = computeBlobCrc(tmp);
    return tmp;
  }

private:
  SaveStatus last_save_status_ = SaveStatus::OK;
  bool loadV6_(ImuCalBlobV6& out, bool* read_failed = nullptr) {
    if (read_failed) *read_failed = false;
    if (kv.getBytesLength(kKeyV6) != sizeof(ImuCalBlobV6)) return false;
    ImuCalBlobV6 tmp;
    if (kv.getBytes(kKeyV6, &tmp, sizeof(tmp)) != sizeof(tmp)) {
      if (read_failed) *read_failed = true;
      return false;
    }
    if (!validateBlob(tmp)) return false;
    out = tmp;
    return true;
  }
};

}  // namespace atoms3r_ical
