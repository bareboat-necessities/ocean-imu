#include <Arduino.h>
#include <M5Unified.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "Bosch/BoschBmi270_ImuCal.h"

using namespace atoms3r_ical;

constexpr float kImuOdrHz = 200.0f;
constexpr uint32_t kLoopPeriodUs = static_cast<uint32_t>(1000000.0f / kImuOdrHz);
constexpr float kMagPeriodMs = 1000.0f / 25.0f;

void waitForNextLoopTick(uint32_t loopStartUs)
{
  const uint32_t elapsedUs = micros() - loopStartUs;
  if (elapsedUs < kLoopPeriodUs)
  {
    delayMicroseconds(kLoopPeriodUs - elapsedUs);
  }
}

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
}

void loop()
{
  const uint32_t loopStartUs = micros();

  ImuSample sample{};
  if (!imu.read(sample))
  {
    waitForNextLoopTick(loopStartUs);
    return;
  }

  const float dtImuMs = imu.fifo().nominalDt() * 1000.0f;
  const bool magValid = imu.magHealthy();

  const float gyroNorthDps = sample.w.x() * RAD_TO_DEG;
  const float gyroEastDps = sample.w.y() * RAD_TO_DEG;
  const float gyroDownDps = sample.w.z() * RAD_TO_DEG;

  const Vector3f magBody = sample.m;

  Serial.printf(
    "FIFO len=%u req=%u got=%u | dt_imu_ms=%.2f | mag_valid=%u | dt_mag_ms~%.2f (target~%.2f) | acc_ned[m/s^2] N=%.3f E=%.3f D=%.3f | gyro_ned[dps] N=%.3f E=%.3f D=%.3f | mag_ned[uT] N=%.2f E=%.2f D=%.2f\n",
    imu.fifo().lastFifoLen(),
    imu.fifo().lastFifoReq(),
    imu.fifo().lastFifoGot(),
    dtImuMs,
    magValid ? 1U : 0U,
    magValid ? static_cast<float>((imu.sampleClockUs64() - imu.lastMagSampleUs64()) * 1e-3) : -1.0f,
    kMagPeriodMs,
    sample.a.x(),
    sample.a.y(),
    sample.a.z(),
    gyroNorthDps,
    gyroEastDps,
    gyroDownDps,
    magBody.x(),
    magBody.y(),
    magBody.z());

  waitForNextLoopTick(loopStartUs);
}
