#pragma once

#include <Arduino.h>
#include <SparkFun_BMI270_Arduino_Library.h>

class MagReader
{
public:
  // Bosch BMM150 ODR selector values (register 0x4C, bits [5:3]).
  static constexpr uint8_t kBmm150Odr2Hz  = 0x01;
  static constexpr uint8_t kBmm150Odr6Hz  = 0x02;
  static constexpr uint8_t kBmm150Odr8Hz  = 0x03;
  static constexpr uint8_t kBmm150Odr10Hz = 0x00;
  static constexpr uint8_t kBmm150Odr15Hz = 0x04;
  static constexpr uint8_t kBmm150Odr20Hz = 0x05;
  static constexpr uint8_t kBmm150Odr25Hz = 0x06;
  static constexpr uint8_t kBmm150Odr30Hz = 0x07;

  static constexpr uint8_t kDefaultBmm150Odr = kBmm150Odr25Hz;
  static constexpr uint8_t kDefaultAuxOdr    = BMI2_AUX_ODR_100HZ;

  static constexpr float odrSettingToHz(uint8_t bmm150Odr)
  {
    switch (bmm150Odr & kBmm150OdrMask)
    {
      case kBmm150Odr2Hz:  return 2.0f;
      case kBmm150Odr6Hz:  return 6.0f;
      case kBmm150Odr8Hz:  return 8.0f;
      case kBmm150Odr10Hz: return 10.0f;
      case kBmm150Odr15Hz: return 15.0f;
      case kBmm150Odr20Hz: return 20.0f;
      case kBmm150Odr25Hz: return 25.0f;
      case kBmm150Odr30Hz: return 30.0f;
      default:             return 0.0f;
    }
  }

  static constexpr float odrSettingToPeriodMs(uint8_t bmm150Odr)
  {
    const float hz = odrSettingToHz(bmm150Odr);
    return (hz > 0.0f) ? (1000.0f / hz) : 0.0f;
  }

  struct Sample
  {
    // Compensated magnetic field in microtesla.
    // Field names kept for sketch compatibility.
    int16_t x = 0;
    int16_t y = 0;
    int16_t z = 0;

    bool valid   = false; // usable cached sample exists and is not stale
    bool updated = false; // a new mag sample arrived this call
    float deltaMs = 0.0f; // age since last NEW mag update
  };

  bool ok() const { return initialized_ && trimValid_; }

  bool begin(BMI270 &imu, uint8_t bmm150Odr = kDefaultBmm150Odr, uint8_t auxOdr = kDefaultAuxOdr)
  {
    resetState_();

    if (imu.enableFeature(BMI2_AUX) != BMI2_OK)
    {
      Serial.println("Failed to enable BMI270 AUX interface.");
      return false;
    }

    if (imu.setAuxPullUps(BMI2_ASDA_PUPSEL_2K) != BMI2_OK)
    {
      Serial.println("Failed to configure BMI270 AUX pullups.");
      return false;
    }

    uint8_t detectedAddress = 0;
    if (!probeAndConfigureAddress_(imu, detectedAddress, auxOdr))
    {
      Serial.println("Failed to detect BMM150 on BMI270 AUX bus.");
      return false;
    }

    bmm150I2cAddr_ = detectedAddress;

    if (imu.writeAux(kBmm150PowerCtrlReg, 0x01) != BMI2_OK)
    {
      Serial.println("Failed to power up BMM150.");
      return false;
    }
    delay(10);

    const uint8_t opMode = kBmm150ModeNormal | ((bmm150Odr & kBmm150OdrMask) << kBmm150OdrShift);
    if (imu.writeAux(kBmm150OpModeReg, opMode) != BMI2_OK)
    {
      Serial.println("Failed to configure BMM150 op mode.");
      return false;
    }

    // Use the same repetition settings you already had.
    if (imu.writeAux(kBmm150XyRepReg, 0x04) != BMI2_OK ||
        imu.writeAux(kBmm150ZRepReg, 0x0F) != BMI2_OK)
    {
      Serial.println("Failed to configure BMM150 repetition settings.");
      return false;
    }

    delay(10);

    if (!readTrimData_(imu))
    {
      Serial.println("Failed to read BMM150 trim registers.");
      return false;
    }

    configuredPeriodMs_ = odrSettingToPeriodMs(bmm150Odr);
    staleAfterMs_ = (configuredPeriodMs_ > 0.0f) ? (3.0f * configuredPeriodMs_) : 120.0f;

    initialized_ = true;
    Serial.printf("BMM150 initialized over BMI270 AUX manual-read mode (addr=0x%02X).\n", bmm150I2cAddr_);
    return true;
  }

  // Kept for sketch compatibility: x/y/z are already compensated microtesla now.
  static float rawToMicroTesla(int16_t compensatedUt)
  {
    return static_cast<float>(compensatedUt);
  }

  Sample readLatestSample(BMI270 &imu, uint32_t sampleTimestampMs)
  {
    Sample out = lastSample_;

    const float ageMs = ageSinceLastUpdateMs_(sampleTimestampMs);
    out.deltaMs = ageMs;
    out.updated = false;
    out.valid = hasLastMagSample_ && (ageMs <= staleAfterMs_);

    if (!ok())
    {
      out.valid = false;
      return out;
    }

    uint8_t buf[8] = {0};
    if (!readAuxBytes_(imu, kBmm150DataXLsbReg, buf, sizeof(buf)))
    {
      // Keep last sample valid until stale timeout expires.
      return out;
    }

    int16_t rawX = 0;
    int16_t rawY = 0;
    int16_t rawZ = 0;
    uint16_t rawRhall = 0;
    decodeBMM150Raw_(buf, rawX, rawY, rawZ, rawRhall);

    const bool changed =
      !hasLastRaw_ ||
      rawX != lastRawX_ ||
      rawY != lastRawY_ ||
      rawZ != lastRawZ_ ||
      rawRhall != lastRawRhall_;

    if (!changed)
    {
      return out;
    }

    const int16_t compX = compensateX_(rawX, rawRhall);
    const int16_t compY = compensateY_(rawY, rawRhall);
    const int16_t compZ = compensateZ_(rawZ, rawRhall);

    if (compX == kOverflowOutput || compY == kOverflowOutput || compZ == kOverflowOutput)
    {
      return out;
    }

    lastRawX_ = rawX;
    lastRawY_ = rawY;
    lastRawZ_ = rawZ;
    lastRawRhall_ = rawRhall;
    hasLastRaw_ = true;

    lastSample_.x = compX;
    lastSample_.y = compY;
    lastSample_.z = compZ;
    lastSample_.valid = true;
    lastSample_.updated = true;
    lastSample_.deltaMs = ageMs;

    hasLastMagSample_ = true;
    hasLastMagTimestamp_ = true;
    lastMagTimestampMs_ = sampleTimestampMs;

    return lastSample_;
  }

private:
  struct TrimData
  {
    int8_t   dig_x1 = 0;
    int8_t   dig_y1 = 0;
    int8_t   dig_x2 = 0;
    int8_t   dig_y2 = 0;
    uint16_t dig_z1 = 0;
    int16_t  dig_z2 = 0;
    int16_t  dig_z3 = 0;
    int16_t  dig_z4 = 0;
    uint8_t  dig_xy1 = 0;
    int8_t   dig_xy2 = 0;
    uint16_t dig_xyz1 = 0;
  };

  static constexpr uint8_t kBmm150ChipIdReg    = 0x40;
  static constexpr uint8_t kBmm150ChipId       = 0x32;
  static constexpr uint8_t kBmm150PowerCtrlReg = 0x4B;
  static constexpr uint8_t kBmm150OpModeReg    = 0x4C;
  static constexpr uint8_t kBmm150ModeNormal   = 0x00;
  static constexpr uint8_t kBmm150OdrShift     = 3;
  static constexpr uint8_t kBmm150OdrMask      = 0x07;
  static constexpr uint8_t kBmm150XyRepReg     = 0x51;
  static constexpr uint8_t kBmm150ZRepReg      = 0x52;
  static constexpr uint8_t kBmm150DataXLsbReg  = 0x42;

  static constexpr int16_t kOverflowOutput      = -32768;
  static constexpr int16_t kOverflowAdcXy       = -4096;
  static constexpr int16_t kOverflowAdcZ        = -16384;
  static constexpr int16_t kNegSaturationZ      = -32767;
  static constexpr int16_t kPosSaturationZ      = 32767;

  static bool configureAux_(BMI270 &imu, uint8_t auxAddress, uint8_t auxOdr)
  {
    bmi2_sens_config auxConfig{};
    auxConfig.type = BMI2_AUX;
    auxConfig.cfg.aux.aux_en = BMI2_ENABLE;
    auxConfig.cfg.aux.manual_en = BMI2_ENABLE;
    auxConfig.cfg.aux.man_rd_burst = BMI2_AUX_READ_LEN_3; // 8-byte manual read
    auxConfig.cfg.aux.aux_rd_burst = BMI2_AUX_READ_LEN_3;
    auxConfig.cfg.aux.odr = auxOdr;
    auxConfig.cfg.aux.i2c_device_addr = auxAddress;
    auxConfig.cfg.aux.read_addr = kBmm150DataXLsbReg;
    auxConfig.cfg.aux.fcu_write_en = BMI2_ENABLE;
    auxConfig.cfg.aux.offset = 0;

    return imu.setConfig(auxConfig) == BMI2_OK;
  }

  static bool probeAndConfigureAddress_(BMI270 &imu, uint8_t &detectedAddress, uint8_t auxOdr)
  {
    constexpr uint8_t addressCandidates[] = {
      0x10, // default 7-bit BMM150 address
      0x12  // alternate 7-bit address
    };

    for (uint8_t auxAddress : addressCandidates)
    {
      if (!configureAux_(imu, auxAddress, auxOdr))
      {
        continue;
      }

      if (imu.writeAux(kBmm150PowerCtrlReg, 0x01) != BMI2_OK)
      {
        continue;
      }
      delay(3);

      uint8_t chipId = 0;
      if (!readAuxBytes_(imu, kBmm150ChipIdReg, &chipId, 1))
      {
        continue;
      }

      if (chipId == kBmm150ChipId)
      {
        detectedAddress = auxAddress;
        return true;
      }
    }

    return false;
  }

  static bool readAuxBytes_(BMI270 &imu, uint8_t addr, uint8_t *dst, uint8_t numBytes)
  {
    if (dst == nullptr || numBytes == 0 || numBytes > BMI2_AUX_NUM_BYTES)
    {
      return false;
    }

    if (imu.readAux(addr, numBytes) != BMI2_OK)
    {
      return false;
    }

    for (uint8_t i = 0; i < numBytes; ++i)
    {
      dst[i] = imu.data.auxData[i];
    }

    return true;
  }

  bool readTrimData_(BMI270 &imu)
  {
    uint8_t buf2[2] = {0};
    uint8_t buf4[4] = {0};
    uint8_t buf6[6] = {0};

    if (!readAuxBytes_(imu, 0x5D, buf2, sizeof(buf2))) return false;
    trim_.dig_x1 = static_cast<int8_t>(buf2[0]);
    trim_.dig_y1 = static_cast<int8_t>(buf2[1]);

    if (!readAuxBytes_(imu, 0x62, buf4, sizeof(buf4))) return false;
    trim_.dig_z4 = static_cast<int16_t>(
      static_cast<uint16_t>(buf4[0]) |
      (static_cast<uint16_t>(buf4[1]) << 8));
    trim_.dig_x2 = static_cast<int8_t>(buf4[2]);
    trim_.dig_y2 = static_cast<int8_t>(buf4[3]);

    if (!readAuxBytes_(imu, 0x68, buf6, sizeof(buf6))) return false;
    trim_.dig_z2 = static_cast<int16_t>(
      static_cast<uint16_t>(buf6[0]) |
      (static_cast<uint16_t>(buf6[1]) << 8));
    trim_.dig_z1 = static_cast<uint16_t>(
      static_cast<uint16_t>(buf6[2]) |
      (static_cast<uint16_t>(buf6[3]) << 8));
    trim_.dig_xyz1 = static_cast<uint16_t>(
      static_cast<uint16_t>(buf6[4]) |
      (static_cast<uint16_t>(buf6[5]) << 8));

    if (!readAuxBytes_(imu, 0x6E, buf4, sizeof(buf4))) return false;
    trim_.dig_z3 = static_cast<int16_t>(
      static_cast<uint16_t>(buf4[0]) |
      (static_cast<uint16_t>(buf4[1]) << 8));
    trim_.dig_xy2 = static_cast<int8_t>(buf4[2]);
    trim_.dig_xy1 = buf4[3];

    trimValid_ = true;
    return true;
  }

  static int16_t signExtend_(int16_t value, uint8_t bits)
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

  static void decodeBMM150Raw_(
    const uint8_t auxData[8],
    int16_t &magX,
    int16_t &magY,
    int16_t &magZ,
    uint16_t &rhall)
  {
    const int16_t rawX = static_cast<int16_t>(
      (static_cast<int16_t>(auxData[1]) << 5) |
      (static_cast<int16_t>(auxData[0]) >> 3));

    const int16_t rawY = static_cast<int16_t>(
      (static_cast<int16_t>(auxData[3]) << 5) |
      (static_cast<int16_t>(auxData[2]) >> 3));

    const int16_t rawZ = static_cast<int16_t>(
      (static_cast<int16_t>(auxData[5]) << 7) |
      (static_cast<int16_t>(auxData[4]) >> 1));

    const uint16_t rawR = static_cast<uint16_t>(
      (static_cast<uint16_t>(auxData[7]) << 6) |
      (static_cast<uint16_t>(auxData[6]) >> 2));

    magX = signExtend_(rawX, 13);
    magY = signExtend_(rawY, 13);
    magZ = signExtend_(rawZ, 15);
    rhall = rawR;
  }

  int16_t compensateX_(int16_t mag_data_x, uint16_t data_rhall) const
  {
    if (!trimValid_ || mag_data_x == kOverflowAdcXy)
    {
      return kOverflowOutput;
    }

    uint16_t process_comp_x0 = 0;
    if (data_rhall != 0)
    {
      process_comp_x0 = data_rhall;
    }
    else if (trim_.dig_xyz1 != 0)
    {
      process_comp_x0 = trim_.dig_xyz1;
    }
    else
    {
      return kOverflowOutput;
    }

    const int32_t process_comp_x1 = static_cast<int32_t>(trim_.dig_xyz1) * 16384;
    const uint16_t process_comp_x2 =
      static_cast<uint16_t>(process_comp_x1 / process_comp_x0) - static_cast<uint16_t>(0x4000);
    int16_t retval = static_cast<int16_t>(process_comp_x2);

    const int32_t process_comp_x3 = static_cast<int32_t>(retval) * static_cast<int32_t>(retval);
    const int32_t process_comp_x4 = static_cast<int32_t>(trim_.dig_xy2) * (process_comp_x3 / 128);
    const int32_t process_comp_x5 = static_cast<int32_t>(static_cast<int16_t>(trim_.dig_xy1) * 128);
    const int32_t process_comp_x6 = static_cast<int32_t>(retval) * process_comp_x5;
    const int32_t process_comp_x7 = ((process_comp_x4 + process_comp_x6) / 512) + static_cast<int32_t>(0x100000);
    const int32_t process_comp_x8 = static_cast<int32_t>(static_cast<int16_t>(trim_.dig_x2) + static_cast<int16_t>(0xA0));
    const int32_t process_comp_x9 = (process_comp_x7 * process_comp_x8) / 4096;
    const int32_t process_comp_x10 = static_cast<int32_t>(mag_data_x) * process_comp_x9;

    retval = static_cast<int16_t>(process_comp_x10 / 8192);
    retval = static_cast<int16_t>((retval + (static_cast<int16_t>(trim_.dig_x1) * 8)) / 16);
    return retval;
  }

  int16_t compensateY_(int16_t mag_data_y, uint16_t data_rhall) const
  {
    if (!trimValid_ || mag_data_y == kOverflowAdcXy)
    {
      return kOverflowOutput;
    }

    uint16_t process_comp_y0 = 0;
    if (data_rhall != 0)
    {
      process_comp_y0 = data_rhall;
    }
    else if (trim_.dig_xyz1 != 0)
    {
      process_comp_y0 = trim_.dig_xyz1;
    }
    else
    {
      return kOverflowOutput;
    }

    const int32_t process_comp_y1 =
      (static_cast<int32_t>(trim_.dig_xyz1) * 16384) / process_comp_y0;
    const uint16_t process_comp_y2 =
      static_cast<uint16_t>(process_comp_y1) - static_cast<uint16_t>(0x4000);
    int16_t retval = static_cast<int16_t>(process_comp_y2);

    const int32_t process_comp_y3 = static_cast<int32_t>(retval) * static_cast<int32_t>(retval);
    const int32_t process_comp_y4 = static_cast<int32_t>(trim_.dig_xy2) * (process_comp_y3 / 128);
    const int32_t process_comp_y5 = static_cast<int32_t>(static_cast<int16_t>(trim_.dig_xy1) * 128);
    const int32_t process_comp_y6 =
      (process_comp_y4 + (static_cast<int32_t>(retval) * process_comp_y5)) / 512;
    const int32_t process_comp_y7 =
      static_cast<int32_t>(static_cast<int16_t>(trim_.dig_y2) + static_cast<int16_t>(0xA0));
    const int32_t process_comp_y8 =
      ((process_comp_y6 + static_cast<int32_t>(0x100000)) * process_comp_y7) / 4096;
    const int32_t process_comp_y9 = static_cast<int32_t>(mag_data_y) * process_comp_y8;

    retval = static_cast<int16_t>(process_comp_y9 / 8192);
    retval = static_cast<int16_t>((retval + (static_cast<int16_t>(trim_.dig_y1) * 8)) / 16);
    return retval;
  }

  int16_t compensateZ_(int16_t mag_data_z, uint16_t data_rhall) const
  {
    if (!trimValid_ || mag_data_z == kOverflowAdcZ)
    {
      return kOverflowOutput;
    }

    if ((trim_.dig_z2 == 0) || (trim_.dig_z1 == 0) || (data_rhall == 0) || (trim_.dig_xyz1 == 0))
    {
      return kOverflowOutput;
    }

    const int16_t process_comp_z0 = static_cast<int16_t>(data_rhall) - static_cast<int16_t>(trim_.dig_xyz1);
    const int32_t process_comp_z1 =
      (static_cast<int32_t>(trim_.dig_z3) * static_cast<int32_t>(process_comp_z0)) / 4;
    const int32_t process_comp_z2 =
      (static_cast<int32_t>(mag_data_z - trim_.dig_z4)) * 32768;
    const int32_t process_comp_z3 =
      static_cast<int32_t>(trim_.dig_z1) * (static_cast<int16_t>(data_rhall) * 2);
    const int16_t process_comp_z4 =
      static_cast<int16_t>((process_comp_z3 + 32768) / 65536);

    const int32_t denom = static_cast<int32_t>(trim_.dig_z2) + static_cast<int32_t>(process_comp_z4);
    if (denom == 0)
    {
      return kOverflowOutput;
    }

    int32_t retval = (process_comp_z2 - process_comp_z1) / denom;

    if (retval > kPosSaturationZ)
    {
      retval = kPosSaturationZ;
    }
    else if (retval < kNegSaturationZ)
    {
      retval = kNegSaturationZ;
    }

    retval /= 16; // LSB -> microtesla
    return static_cast<int16_t>(retval);
  }

  float ageSinceLastUpdateMs_(uint32_t nowMs) const
  {
    if (!hasLastMagTimestamp_)
    {
      return 0.0f;
    }

    if (nowMs >= lastMagTimestampMs_)
    {
      return static_cast<float>(nowMs - lastMagTimestampMs_);
    }

    return 0.0f;
  }

  void resetState_()
  {
    initialized_ = false;
    trimValid_ = false;
    configuredPeriodMs_ = 0.0f;
    staleAfterMs_ = 120.0f;

    hasLastMagSample_ = false;
    hasLastMagTimestamp_ = false;
    lastMagTimestampMs_ = 0;

    hasLastRaw_ = false;
    lastRawX_ = 0;
    lastRawY_ = 0;
    lastRawZ_ = 0;
    lastRawRhall_ = 0;

    lastSample_ = Sample{};
    trim_ = TrimData{};
    bmm150I2cAddr_ = 0;
  }

  bool initialized_ = false;
  bool trimValid_ = false;

  float configuredPeriodMs_ = 0.0f;
  float staleAfterMs_ = 120.0f;

  bool hasLastMagSample_ = false;
  bool hasLastMagTimestamp_ = false;
  uint32_t lastMagTimestampMs_ = 0;

  bool hasLastRaw_ = false;
  int16_t lastRawX_ = 0;
  int16_t lastRawY_ = 0;
  int16_t lastRawZ_ = 0;
  uint16_t lastRawRhall_ = 0;

  uint8_t bmm150I2cAddr_ = 0;
  TrimData trim_{};
  Sample lastSample_{};
};
