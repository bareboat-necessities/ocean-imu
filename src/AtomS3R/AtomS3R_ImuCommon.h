#pragma once

#include <Arduino.h>
#include <stdint.h>

namespace atoms3r_ical {

constexpr float kDefaultImuOdrHz = 200.0f;
constexpr uint32_t kDefaultImuLogPeriodUs = 50000u;          // 20 Hz
constexpr uint32_t kDefaultNoSampleLogPeriodUs = 500000u;    // 2 Hz

constexpr uint32_t loopPeriodUsFromHz(float hz)
{
  return (hz > 0.0f)
    ? static_cast<uint32_t>(1000000.0f / hz)
    : static_cast<uint32_t>(1000000.0f / kDefaultImuOdrHz);
}

inline void waitForNextLoopTick(uint32_t loopStartUs, uint32_t loopPeriodUs)
{
  const uint32_t elapsedUs = micros() - loopStartUs;
  if (elapsedUs < loopPeriodUs)
  {
    delayMicroseconds(loopPeriodUs - elapsedUs);
  }
}

inline float usToMs(uint64_t us)
{
  return static_cast<float>(us) * 1.0e-3f;
}

struct ImuLogRow
{
  // FIFO / batch stats
  int32_t fifoLen = -1;
  int32_t fifoReq = -1;
  int32_t fifoGot = -1;

  // Timing
  float dtImuMs = -1.0f;
  bool magValid = false;
  bool magUpdated = false;
  float magAgeMs = -1.0f;
  float magStepMs = -1.0f;
  float lastMagDtMs = -1.0f;

  // Data in NED/body-NED convention used by the demos
  float accN = 0.0f;
  float accE = 0.0f;
  float accD = 0.0f;

  float gyrN = 0.0f;
  float gyrE = 0.0f;
  float gyrD = 0.0f;

  float magN = 0.0f;
  float magE = 0.0f;
  float magD = 0.0f;
};

struct ImuLoopLogState
{
  uint32_t lastLogUs = 0u;
  uint32_t lastNoSampleLogUs = 0u;
  uint32_t consecutiveReadFailures = 0u;

  void markReadFailure()
  {
    ++consecutiveReadFailures;
  }

  void markReadSuccess()
  {
    consecutiveReadFailures = 0u;
  }

  bool shouldLogNoSample(uint32_t nowUs, uint32_t periodUs = kDefaultNoSampleLogPeriodUs)
  {
    if (static_cast<uint32_t>(nowUs - lastNoSampleLogUs) < periodUs)
    {
      return false;
    }

    lastNoSampleLogUs = nowUs;
    return true;
  }

  bool shouldLogSample(uint32_t nowUs, uint32_t periodUs = kDefaultImuLogPeriodUs)
  {
    if (static_cast<uint32_t>(nowUs - lastLogUs) < periodUs)
    {
      return false;
    }

    lastLogUs = nowUs;
    return true;
  }
};

inline void printImuLogRow(const ImuLogRow& row)
{
  Serial.printf(
    "FIFO len=%ld req=%ld got=%ld | dt_imu_ms=%.2f | "
    "mag_valid=%u mag_updated=%u | mag_age_ms=%.2f | mag_step_ms=%.2f | last_mag_dt_ms=%.2f | "
    "acc_ned[m/s^2] N=%.3f E=%.3f D=%.3f | "
    "gyro_ned[dps] N=%.3f E=%.3f D=%.3f | "
    "mag_ned[uT] N=%.2f E=%.2f D=%.2f\n",
    static_cast<long>(row.fifoLen),
    static_cast<long>(row.fifoReq),
    static_cast<long>(row.fifoGot),
    row.dtImuMs,
    row.magValid ? 1U : 0U,
    row.magUpdated ? 1U : 0U,
    row.magAgeMs,
    row.magStepMs,
    row.lastMagDtMs,
    row.accN,
    row.accE,
    row.accD,
    row.gyrN,
    row.gyrE,
    row.gyrD,
    row.magN,
    row.magE,
    row.magD);
}

} // namespace atoms3r_ical
