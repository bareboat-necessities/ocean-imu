#pragma once
// Host stand-in for the AtomS3R IMU acquisition layer. Calibration storage and
// RuntimeCals are the real host-testable implementation; readImuMapped()
// delivers body-NED samples produced by the test's motion source.
#include <Arduino.h>
#include <ArduinoEigenDense.h>
#include <map>
#include <string>
#include <vector>
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include "AtomS3R/AtomS3R_InsDiag.h"

#define ATOMS3R_IMU_MASK_ACCEL (1u << 0)
#define ATOMS3R_IMU_MASK_GYRO (1u << 1)

namespace atoms3r_ical {
using Vector3f = Eigen::Vector3f;

struct ImuSample {
  Vector3f a;     // m/s^2, body NED
  Vector3f w;     // rad/s, body NED
  Vector3f m;     // uT, body NED; NaN when stale
  float tempC;
  uint32_t mask;
  uint32_t sample_us;
  MagAcquisition mag_timing;
};

// Supplied by the test driver: fill a, w, m, tempC and mag_timing for sample_us
// and return true, or return false for a poll that has no new frame.
using SampleSource = bool (*)(uint32_t sample_us, ImuSample& out);
inline SampleSource& sampleSource() { static SampleSource source = nullptr; return source; }

template <class Imu>
inline bool readImuMapped(Imu&, uint32_t update_mask, uint32_t sample_us, ImuSample& out) {
  out.sample_us = sample_us;
  out.mask = update_mask;
  if (!sampleSource()) return false;
  return sampleSource()(sample_us, out);
}

struct MemoryKv {
  std::map<std::string, std::vector<uint8_t>> bytes;
  size_t getBytesLength(const char* k) { return bytes[k].size(); }
  size_t getBytes(const char* k, void* dst, size_t n) {
    if (bytes[k].size() != n) return 0;
    std::memcpy(dst, bytes[k].data(), n);
    return n;
  }
  size_t putBytes(const char* k, const void* src, size_t n) {
    const auto* p = static_cast<const uint8_t*>(src);
    bytes[k] = {p, p + n};
    return n;
  }
  bool remove(const char* k) { bytes.erase(k); return true; }
};
class ImuCalStoreNvs : public ImuCalStoreT<MemoryKv> {};

inline void clearM5UnifiedImuCalibration() {}
inline ImuTimingRegisters readImuTimingRegisters(bool) { return {}; }
inline bool configureAtomS3RMag(Print&, bool = true) { return true; }
inline void printBlobSummary(Print&, const ImuCalBlobV5&) {}
inline void printBlobDetail(Print&, const ImuCalBlobV5&) {}
}  // namespace atoms3r_ical
