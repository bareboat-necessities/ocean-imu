#pragma once

#include <Arduino.h>
#include <SparkFun_BMI270_Arduino_Library.h>

class ImuReader
{
public:
  static constexpr uint8_t kDefaultI2cAddress = BMI2_I2C_PRIM_ADDR;
  static constexpr uint8_t kAlternateI2cAddress = BMI2_I2C_SEC_ADDR;
  static constexpr uint16_t kDefaultFifoWatermarkFrames = 20;
  static constexpr uint8_t kDefaultAccelOdr = BMI2_ACC_ODR_200HZ;
  static constexpr uint8_t kDefaultGyroOdr = BMI2_GYR_ODR_200HZ;
  static constexpr uint8_t kDefaultAccelRange = BMI2_ACC_RANGE_2G;

  struct Sample
  {
    bool valid = false;
    uint16_t framesRead = 0;
    uint32_t timestampMs = 0;
    float imuDeltaMs = 0.0f;
    float accelX = 0.0f;
    float accelY = 0.0f;
    float accelZ = 0.0f;
    float gyroX = 0.0f;
    float gyroY = 0.0f;
    float gyroZ = 0.0f;
  };

  bool begin(TwoWire &wire, uint8_t i2cAddress = kDefaultI2cAddress, uint8_t attemptsPerAddress = 3)
  {
    uint8_t addressesToTry[2] = {i2cAddress, kAlternateI2cAddress};
    if (i2cAddress == kAlternateI2cAddress)
    {
      addressesToTry[1] = kDefaultI2cAddress;
    }

    for (uint8_t attempt = 1; attempt <= attemptsPerAddress; ++attempt)
    {
      for (const uint8_t address : addressesToTry)
      {
        if (imu.beginI2C(address, wire) == BMI2_OK)
        {
          Serial.printf(
            "BMI270 detected at I2C address 0x%02X (attempt %u/%u)\n",
            address,
            attempt,
            attemptsPerAddress);
          return true;
        }
      }

      delay(50);
    }

    Serial.printf(
      "BMI270 init failed on 0x%02X and 0x%02X after %u attempts each.\n",
      addressesToTry[0],
      addressesToTry[1],
      attemptsPerAddress);

    return false;
  }

  void configureFIFO(
    uint16_t watermarkFrames = kDefaultFifoWatermarkFrames,
    uint8_t accelOdr = kDefaultAccelOdr,
    uint8_t gyroOdr = kDefaultGyroOdr)
  {
    if (watermarkFrames > kDefaultFifoWatermarkFrames)
    {
      watermarkFrames = kDefaultFifoWatermarkFrames;
    }

    imu.setAccelODR(accelOdr);
    setAccelRange(imu, kDefaultAccelRange);
    imu.setGyroODR(gyroOdr);
    const float accelHz = imuAccelOdrToHz(accelOdr);
    const float gyroHz = imuGyroOdrToHz(gyroOdr);
    const float fifoRateHz = accelHz > gyroHz ? accelHz : gyroHz;
    imuFramePeriodMs = fifoRateHz > 0.0f ? (1000.0f / fifoRateHz) : kDefaultImuFramePeriodMs;

    BMI270_FIFOConfig fifoConfig;
    // Keep FIFO dedicated to IMU data. Magnetometer is read directly through
    // the BMI270 AUX bridge to avoid AUX payload interactions with FIFO timing.
    fifoConfig.flags = BMI2_FIFO_ACC_EN | BMI2_FIFO_GYR_EN | BMI2_FIFO_TIME_EN;
    fifoConfig.watermark = watermarkFrames;
    fifoConfig.accelDownSample = BMI2_FIFO_DOWN_SAMPLE_1;
    fifoConfig.gyroDownSample = BMI2_FIFO_DOWN_SAMPLE_1;
    fifoConfig.accelFilter = BMI2_ENABLE;
    fifoConfig.gyroFilter = BMI2_ENABLE;
    fifoConfig.selfWakeUp = BMI2_ENABLE;

    imu.setFIFOConfig(fifoConfig);
  }

  Sample readLatestSample()
  {
    Sample sample{};
    uint16_t fifoLength = 0;
    imu.getFIFOLength(&fifoLength);

    if (fifoLength == 0)
    {
      return sample;
    }

    uint16_t framesRead = fifoLength;
    if (framesRead > kDefaultFifoWatermarkFrames)
    {
      framesRead = kDefaultFifoWatermarkFrames;
    }

    imu.getFIFOData(fifoFrames, &framesRead);
    if (framesRead == 0)
    {
      return sample;
    }

    const BMI270_SensorData &latest = fifoFrames[framesRead - 1];
    sample.valid = true;
    sample.framesRead = framesRead;

    sample.timestampMs = latest.sensorTimeMillis;

    if (hasLastImuTimestamp)
    {
      uint32_t timestampDeltaMs = 0;

      if (sample.timestampMs > lastImuTimestampMs)
      {
        timestampDeltaMs = sample.timestampMs - lastImuTimestampMs;
      }

      // If FIFO sensor time is missing or non-monotonic, estimate elapsed time
      // from the number of frames drained at the configured IMU ODR.
      if (timestampDeltaMs == 0)
      {
        timestampDeltaMs = static_cast<uint32_t>(framesRead * imuFramePeriodMs);
        sample.timestampMs = lastImuTimestampMs + timestampDeltaMs;
      }

      // Report per-sample IMU period rather than per-batch elapsed time.
      // This keeps dt stable even when FIFO drains 2-3 frames at once.
      sample.imuDeltaMs = static_cast<float>(timestampDeltaMs) / static_cast<float>(framesRead);
    }

    hasLastImuTimestamp = true;
    lastImuTimestampMs = sample.timestampMs;

    // SparkFun BMI270 reports acceleration in g; convert to SI units.
    sample.accelX = latest.accelX * kStandardGravity;
    sample.accelY = latest.accelY * kStandardGravity;
    sample.accelZ = latest.accelZ * kStandardGravity;

    // SparkFun BMI270 reports gyro in deg/s.
    sample.gyroX = latest.gyroX;
    sample.gyroY = latest.gyroY;
    sample.gyroZ = latest.gyroZ;
    return sample;
  }

  static void sensorToNED(float sensorX, float sensorY, float sensorZ, float &north, float &east, float &down)
  {
    north = sensorY;
    east = sensorX;
    down = -sensorZ;
  }

  BMI270 &sensor() { return imu; }

private:
  static constexpr float kDefaultImuFramePeriodMs = 10.0f;
  static constexpr float kStandardGravity = 9.80665f;

  static float imuAccelOdrToHz(uint8_t accelOdr)
  {
    switch (accelOdr)
    {
      case BMI2_ACC_ODR_0_78HZ: return 0.78f;
      case BMI2_ACC_ODR_1_56HZ: return 1.56f;
      case BMI2_ACC_ODR_3_12HZ: return 3.12f;
      case BMI2_ACC_ODR_6_25HZ: return 6.25f;
      case BMI2_ACC_ODR_12_5HZ: return 12.5f;
      case BMI2_ACC_ODR_25HZ: return 25.0f;
      case BMI2_ACC_ODR_50HZ: return 50.0f;
      case BMI2_ACC_ODR_100HZ: return 100.0f;
      case BMI2_ACC_ODR_200HZ: return 200.0f;
      case BMI2_ACC_ODR_400HZ: return 400.0f;
      case BMI2_ACC_ODR_800HZ: return 800.0f;
      case BMI2_ACC_ODR_1600HZ: return 1600.0f;
      default: return 0.0f;
    }
  }

  static float imuGyroOdrToHz(uint8_t gyroOdr)
  {
    switch (gyroOdr)
    {
      case BMI2_GYR_ODR_25HZ: return 25.0f;
      case BMI2_GYR_ODR_50HZ: return 50.0f;
      case BMI2_GYR_ODR_100HZ: return 100.0f;
      case BMI2_GYR_ODR_200HZ: return 200.0f;
      case BMI2_GYR_ODR_400HZ: return 400.0f;
      case BMI2_GYR_ODR_800HZ: return 800.0f;
      case BMI2_GYR_ODR_1600HZ: return 1600.0f;
      case BMI2_GYR_ODR_3200HZ: return 3200.0f;
      default: return 0.0f;
    }
  }

  static void setAccelRange(BMI270 &imu, uint8_t range)
  {
    bmi2_sens_config config{};
    config.type = BMI2_ACCEL;

    if (imu.getConfig(&config) != BMI2_OK)
    {
      return;
    }

    config.cfg.acc.range = range;
    imu.setConfig(config);
  }

  BMI270 imu;
  BMI270_SensorData fifoFrames[kDefaultFifoWatermarkFrames];
  bool hasLastImuTimestamp = false;
  uint32_t lastImuTimestampMs = 0;
  float imuFramePeriodMs = kDefaultImuFramePeriodMs;
};