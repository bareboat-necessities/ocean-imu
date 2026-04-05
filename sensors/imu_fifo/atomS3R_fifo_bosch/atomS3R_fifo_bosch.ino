#include <Arduino.h>
#include <M5Unified.h>
#include <ArduinoOceanImu.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_ImuCommon.h"
#include "Bosch/BoschBmi270_ImuCal.h"

using namespace atoms3r_ical;

constexpr float kImuOdrHz = kDefaultImuOdrHz;
constexpr uint32_t kLoopPeriodUs = loopPeriodUsFromHz(kImuOdrHz);
constexpr uint32_t kLogPeriodUs = 50000u;       // 20 Hz serial output keeps logs readable at 115200 baud.
constexpr uint32_t kNoSampleLogPeriodUs = 500000u; // 2 Hz diagnostic when FIFO produces no sample.

BoschBmi270_ImuCal imu;

void setup()
{
  Serial.begin(115200);
  delay(1000);
  Serial.println("M5Stack AtomS3R BMI270 FIFO demo (Bosch path + BMM150 aux)");

  auto cfg = M5.config();
  M5.begin(cfg);

  BoschBmi270_ImuCal::Config imuCfg;
  imuCfg.bmi270_addr = 0x68;
  imuCfg.ag_hz = kImuOdrHz;
  imuCfg.enable_mag_aux = true;
  imuCfg.mag_bmm150_addr = 0x10;
  imuCfg.mag_aux_odr_hz = 25.0f;
  imuCfg.mag_startup_settle_ms = 3;
  imuCfg.mag_verify_first_read = true;
  imuCfg.mag_stale_min_us = 75000u;
  imuCfg.mag_stale_factor = 3u;
  imuCfg.enable_mag_recovery = true;
  imuCfg.mag_recover_after_failures = 6u;
  imuCfg.mag_recover_cooldown_us = 1000000u;
  imuCfg.tempC_default = 35.0f;
  imuCfg.i2c_hz = 400000u;

  if (!imu.begin(M5.In_I2C, imuCfg))
  {
    Serial.printf("Bosch IMU init failed: %s\n", imu.lastErrorString());
    Serial.printf("FIFO detail: %s\n", imu.fifo().lastErrorString());
    Serial.printf("FIFO init path: %s\n", imu.fifo().initPathString());
    Serial.printf("FIFO Bosch init rslt: %d\n", static_cast<int>(imu.fifo().lastBoschInitResult()));
    while (true)
    {
      delay(1000);
    }
  }

  Serial.printf("Bosch FIFO initialized via %s\n", imu.fifo().initPathString());
  Serial.printf("Requested mag rate: %.2f Hz, effective mag rate: %.2f Hz, mag step: %.2f ms\n",
                imuCfg.mag_aux_odr_hz,
                imu.effectiveMagHz(),
                usToMs(imu.magStepUs64()));
}

void loop()
{
  const uint32_t loopStartUs = micros();
  static uint32_t lastLogUs = 0u;
  static uint32_t lastNoSampleLogUs = 0u;
  static uint32_t consecutiveReadFailures = 0u;

  ImuSample sample{};
  if (!imu.read(sample))
  {
    ++consecutiveReadFailures;
    const uint32_t nowUs = micros();
    if (static_cast<uint32_t>(nowUs - lastNoSampleLogUs) >= kNoSampleLogPeriodUs)
    {
      lastNoSampleLogUs = nowUs;
      Serial.printf("Waiting for FIFO sample (failures=%lu): imu=%s fifo=%s len=%u req=%u got=%u\n",
                    static_cast<unsigned long>(consecutiveReadFailures),
                    imu.lastErrorString(),
                    imu.fifo().lastErrorString(),
                    static_cast<unsigned int>(imu.fifo().lastFifoLen()),
                    static_cast<unsigned int>(imu.fifo().lastFifoReq()),
                    static_cast<unsigned int>(imu.fifo().lastFifoGot()));
    }
    waitForNextLoopTick(loopStartUs, kLoopPeriodUs);
    return;
  }
  consecutiveReadFailures = 0u;

  ImuLogRow row{};
  row.fifoLen = static_cast<int32_t>(imu.fifo().lastFifoLen());
  row.fifoReq = static_cast<int32_t>(imu.fifo().lastFifoReq());
  row.fifoGot = static_cast<int32_t>(imu.fifo().lastFifoGot());

  row.dtImuMs = imu.fifo().nominalDt() * 1000.0f;
  row.magValid = imu.magHealthy();
  row.magUpdated = imu.magUpdatedThisRead();
  row.magAgeMs = row.magValid ? usToMs(imu.sampleClockUs64() - imu.lastMagSampleUs64()) : -1.0f;
  row.magStepMs = (imu.magStepUs64() > 0u) ? usToMs(imu.magStepUs64()) : -1.0f;
  row.lastMagDtMs = (imu.lastMagPeriodUs64() > 0u) ? usToMs(imu.lastMagPeriodUs64()) : -1.0f;

  row.accN = sample.a.x();
  row.accE = sample.a.y();
  row.accD = sample.a.z();

  row.gyrN = sample.w.x() * RAD_TO_DEG;
  row.gyrE = sample.w.y() * RAD_TO_DEG;
  row.gyrD = sample.w.z() * RAD_TO_DEG;

  row.magN = sample.m.x();
  row.magE = sample.m.y();
  row.magD = sample.m.z();

  const uint32_t nowUs = micros();
  if (static_cast<uint32_t>(nowUs - lastLogUs) >= kLogPeriodUs)
  {
    lastLogUs = nowUs;
    printImuLogRow(row);
  }
  waitForNextLoopTick(loopStartUs, kLoopPeriodUs);
}
