#pragma once

#include <Arduino.h>
#include <SparkFun_BMI270_Arduino_Library.h>

class MagReader
{
public:
  // Bosch BMM150 ODR selector values (register 0x4C, bits [5:3]).
  static constexpr uint8_t kBmm150Odr2Hz = 0x01;
  static constexpr uint8_t kBmm150Odr6Hz = 0x02;
  static constexpr uint8_t kBmm150Odr8Hz = 0x03;
  static constexpr uint8_t kBmm150Odr10Hz = 0x00;
  static constexpr uint8_t kBmm150Odr15Hz = 0x04;
  static constexpr uint8_t kBmm150Odr20Hz = 0x05;
  static constexpr uint8_t kBmm150Odr25Hz = 0x06;
  static constexpr uint8_t kBmm150Odr30Hz = 0x07;

  static constexpr uint8_t kDefaultBmm150Odr = kBmm150Odr25Hz;
  static constexpr uint8_t kDefaultAuxOdr = BMI2_AUX_ODR_100HZ;

  static constexpr float odrSettingToHz(uint8_t bmm150Odr)
  {
    switch (bmm150Odr & kBmm150OdrMask)
    {
      case kBmm150Odr2Hz: return 2.0f;
      case kBmm150Odr6Hz: return 6.0f;
      case kBmm150Odr8Hz: return 8.0f;
      case kBmm150Odr10Hz: return 10.0f;
      case kBmm150Odr15Hz: return 15.0f;
      case kBmm150Odr20Hz: return 20.0f;
      case kBmm150Odr25Hz: return 25.0f;
      case kBmm150Odr30Hz: return 30.0f;
      default: return 0.0f;
    }
  }

  static constexpr float odrSettingToPeriodMs(uint8_t bmm150Odr)
  {
    const float hz = odrSettingToHz(bmm150Odr);
    return hz > 0.0f ? (1000.0f / hz) : 0.0f;
  }

  struct Sample
  {
    int16_t x = 0;
    int16_t y = 0;
    int16_t z = 0;
    bool valid = false;
    bool updated = false;
    float deltaMs = 0.0f;
  };

  bool begin(BMI270 &imu, uint8_t bmm150Odr = kDefaultBmm150Odr, uint8_t auxOdr = kDefaultAuxOdr)
  {
    if (imu.enableFeature(BMI2_AUX) != BMI2_OK)
    {
      Serial.println("Failed to enable BMI270 AUX interface.");
      return false;
    }

    imu.setAuxPullUps(BMI2_ASDA_PUPSEL_2K);

    uint8_t detectedAddress = 0;
    if (!probeAndConfigureAddress(imu, detectedAddress, auxOdr))
    {
      Serial.println("Failed to detect BMM150 on BMI270 AUX bus.");
      return false;
    }

    bmm150I2cAddr = detectedAddress;

    if (imu.writeAux(kBmm150PowerCtrlReg, 0x01) != BMI2_OK)
    {
      Serial.println("Failed to power up BMM150.");
      return false;
    }
    delay(10);

    // Bosch BMM150 normal mode with configurable ODR.
    const uint8_t opMode = kBmm150ModeNormal | ((bmm150Odr & kBmm150OdrMask) << kBmm150OdrShift);
    if (imu.writeAux(kBmm150OpModeReg, opMode) != BMI2_OK)
    {
      Serial.println("Failed to configure BMM150 op mode.");
      return false;
    }

    if (imu.writeAux(kBmm150XyRepReg, 0x04) != BMI2_OK || imu.writeAux(kBmm150ZRepReg, 0x0F) != BMI2_OK)
    {
      Serial.println("Failed to configure BMM150 repetition settings.");
      return false;
    }

    delay(200);
    Serial.printf("BMM150 initialized over BMI270 AUX manual-read mode (addr=0x%02X).\n", bmm150I2cAddr);
    return true;
  }


  static float rawToMicroTesla(int16_t rawValue)
  {
    // BMM150 compensated output is commonly represented at ~16 LSB per uT.
    // We expose a simple bring-up conversion for log readability.
    return static_cast<float>(rawValue) * kMicroTeslaPerLsb;
  }

  Sample readLatestSample(BMI270 &imu, uint32_t sampleTimestampMs)
  {
    Sample sample = lastSample;
    sample.valid = hasLastMagSample;

    float elapsedMs = 0.0f;
    if (hasLastMagTimestamp && sampleTimestampMs >= lastMagTimestampMs)
    {
      elapsedMs = static_cast<float>(sampleTimestampMs - lastMagTimestampMs);
    }

    if (imu.readAux(kBmm150DataXLsbReg, BMI2_AUX_NUM_BYTES) != BMI2_OK)
    {
      sample.valid = false;
      sample.updated = false;
      sample.deltaMs = 0.0f;
      return sample;
    }

    int16_t magX = 0;
    int16_t magY = 0;
    int16_t magZ = 0;
    decodeBMM150Raw(imu.data.auxData, magX, magY, magZ);

    const bool changed = !hasLastMagSample || magX != lastSample.x || magY != lastSample.y || magZ != lastSample.z;
    if (!changed)
    {
      sample.valid = false;
      sample.updated = false;
      sample.deltaMs = elapsedMs;
      return sample;
    }

    sample.x = magX;
    sample.y = magY;
    sample.z = magZ;
    sample.valid = true;
    sample.updated = true;
    sample.deltaMs = elapsedMs;

    hasLastMagTimestamp = true;
    hasLastMagSample = true;
    lastMagTimestampMs = sampleTimestampMs;
    lastSample = sample;

    return sample;
  }

private:
  static constexpr uint8_t kBmm150ChipIdReg = 0x40;
  static constexpr uint8_t kBmm150ChipId = 0x32;
  static constexpr uint8_t kBmm150PowerCtrlReg = 0x4B;
  static constexpr uint8_t kBmm150OpModeReg = 0x4C;
  static constexpr uint8_t kBmm150ModeNormal = 0x00;
  static constexpr uint8_t kBmm150OdrShift = 3;
  static constexpr uint8_t kBmm150OdrMask = 0x07;
  static constexpr uint8_t kBmm150XyRepReg = 0x51;
  static constexpr uint8_t kBmm150ZRepReg = 0x52;
  static constexpr uint8_t kBmm150DataXLsbReg = 0x42;
  static constexpr float kMicroTeslaPerLsb = 1.0f / 16.0f;
  static bool configureAux(BMI270 &imu, uint8_t auxAddress, uint8_t auxOdr)
  {
    bmi2_sens_config auxConfig{};
    auxConfig.type = BMI2_AUX;
    auxConfig.cfg.aux.aux_en = BMI2_ENABLE;
    auxConfig.cfg.aux.manual_en = BMI2_ENABLE;
    // Bosch BMI270 API expects encoded burst-length values (0..3),
    // where BMI2_AUX_READ_LEN_3 maps to 8-byte transfers.
    auxConfig.cfg.aux.man_rd_burst = BMI2_AUX_READ_LEN_3;
    auxConfig.cfg.aux.aux_rd_burst = BMI2_AUX_READ_LEN_3;
    auxConfig.cfg.aux.odr = auxOdr;
    auxConfig.cfg.aux.i2c_device_addr = auxAddress;
    auxConfig.cfg.aux.read_addr = kBmm150DataXLsbReg;
    auxConfig.cfg.aux.fcu_write_en = BMI2_ENABLE;
    auxConfig.cfg.aux.offset = 0;

    return imu.setConfig(auxConfig) == BMI2_OK;
  }

  static bool probeAndConfigureAddress(BMI270 &imu, uint8_t &detectedAddress, uint8_t auxOdr)
  {
    constexpr uint8_t addressCandidates[] = {
      0x10, // BMM150 7-bit default address.
      0x12  // Alternate 7-bit address when SDO is high.
    };

    for (uint8_t auxAddress : addressCandidates)
    {
      if (!configureAux(imu, auxAddress, auxOdr))
      {
        continue;
      }

      // BMM150 can stay in suspend after reset; power it up before probing the
      // chip-ID register on AtomS3R's BMI270 AUX bus.
      if (imu.writeAux(kBmm150PowerCtrlReg, 0x01) != BMI2_OK)
      {
        continue;
      }
      delay(3);

      if (imu.readAux(kBmm150ChipIdReg, 1) != BMI2_OK)
      {
        continue;
      }

      if (imu.data.auxData[0] == kBmm150ChipId)
      {
        detectedAddress = auxAddress;
        return true;
      }
    }

    return false;
  }

  static int16_t signExtend(int16_t value, uint8_t bits)
  {
    const int16_t signBit = static_cast<int16_t>(1U << (bits - 1));
    const int16_t valueMask = static_cast<int16_t>((1U << bits) - 1U);

    value &= valueMask;
    if ((value & signBit) != 0)
    {
      value = static_cast<int16_t>(value - static_cast<int16_t>(1U << bits));
    }

    return value;
  }

  static void decodeBMM150Raw(const uint8_t auxData[BMI2_AUX_NUM_BYTES], int16_t &magX, int16_t &magY, int16_t &magZ)
  {
    const int16_t rawX = (static_cast<int16_t>(auxData[1]) << 5) | (auxData[0] >> 3);
    const int16_t rawY = (static_cast<int16_t>(auxData[3]) << 5) | (auxData[2] >> 3);
    const int16_t rawZ = (static_cast<int16_t>(auxData[5]) << 7) | (auxData[4] >> 1);

    magX = signExtend(rawX, 13);
    magY = signExtend(rawY, 13);
    magZ = signExtend(rawZ, 15);
  }

  bool hasLastMagTimestamp = false;
  bool hasLastMagSample = false;
  uint32_t lastMagTimestampMs = 0;
  uint8_t bmm150I2cAddr = 0;
  Sample lastSample{};
};