#include <Wire.h>

#include "SparkFun/ImuReader.h"
#include "SparkFun/MagReader.h"

constexpr uint8_t ATOMS3R_SENSOR_SDA_PIN = 45;
constexpr uint8_t ATOMS3R_SENSOR_SCL_PIN = 0;
constexpr uint32_t ATOMS3R_SENSOR_I2C_HZ = 400000;
constexpr uint32_t kImuOdrHz = 200;
constexpr uint8_t kImuAccelOdr = ImuReader::kDefaultAccelOdr;
constexpr uint8_t kImuGyroOdr = ImuReader::kDefaultGyroOdr;
constexpr uint8_t kMagOdr = MagReader::kDefaultBmm150Odr;
constexpr uint8_t kAuxOdr = MagReader::kDefaultAuxOdr;
constexpr float kMagPeriodMs = MagReader::odrSettingToPeriodMs(kMagOdr);
constexpr uint32_t kLoopPeriodUs = 1000000UL / kImuOdrHz;

void waitForNextLoopTick(uint32_t loopStartUs)
{
  const uint32_t elapsedUs = micros() - loopStartUs;
  if (elapsedUs < kLoopPeriodUs)
  {
    delayMicroseconds(kLoopPeriodUs - elapsedUs);
  }
}

ImuReader imuReader;
MagReader magReader;
bool magAvailable = false;

void setup()
{
  Serial.begin(115200);
  delay(1000);
  Serial.println("M5Stack AtomS3R BMI270 FIFO demo (+BMM150 aux)");

  Wire.begin(ATOMS3R_SENSOR_SDA_PIN, ATOMS3R_SENSOR_SCL_PIN, ATOMS3R_SENSOR_I2C_HZ);

  if (!imuReader.begin(Wire))
  {
    Serial.println("BMI270 init failed; halting startup.");
    while (true)
    {
      delay(1000);
    }
  }

  magAvailable = magReader.begin(imuReader.sensor(), kMagOdr, kAuxOdr);
  if (!magAvailable)
  {
    Serial.println("Magnetometer init failed; continuing without valid mag data.");
  }

  imuReader.configureFIFO(ImuReader::kDefaultFifoWatermarkFrames, kImuAccelOdr, kImuGyroOdr);
  Serial.println("BMI270 initialized and FIFO configured.");
}

void loop()
{
  const uint32_t loopStartUs = micros();

  const ImuReader::Sample imuSample = imuReader.readLatestSample();
  if (!imuSample.valid)
  {
    waitForNextLoopTick(loopStartUs);
    return;
  }

  const MagReader::Sample magSample =
    magAvailable ? magReader.readLatestSample(imuReader.sensor(), imuSample.timestampMs)
                 : MagReader::Sample{};

  float accNorth = 0.0f;
  float accEast = 0.0f;
  float accDown = 0.0f;
  ImuReader::sensorToNED(imuSample.accelX, imuSample.accelY, imuSample.accelZ, accNorth, accEast, accDown);

  float gyrNorth = 0.0f;
  float gyrEast = 0.0f;
  float gyrDown = 0.0f;
  ImuReader::sensorToNED(imuSample.gyroX, imuSample.gyroY, imuSample.gyroZ, gyrNorth, gyrEast, gyrDown);

  float magNorthUt = 0.0f;
  float magEastUt = 0.0f;
  float magDownUt = 0.0f;
  if (magSample.valid)
  {
    ImuReader::sensorToNED(
      MagReader::rawToMicroTesla(magSample.x),
      MagReader::rawToMicroTesla(magSample.y),
      MagReader::rawToMicroTesla(magSample.z),
      magNorthUt,
      magEastUt,
      magDownUt);
  }

  Serial.printf(
    "FIFO batch=%u | dt_imu_ms=%.2f | mag_valid=%u | mag_updated=%u | mag_age_ms=%s%.2f (target~%.2f) | "
    "acc_ned[m/s^2] N=%.3f E=%.3f D=%.3f | "
    "gyro_ned[dps] N=%.3f E=%.3f D=%.3f | "
    "mag_ned[uT] N=%.2f E=%.2f D=%.2f\n",
    imuSample.framesRead,
    imuSample.imuDeltaMs,
    magSample.valid ? 1U : 0U,
    magSample.updated ? 1U : 0U,
    magSample.updated ? "+" : "~",
    magSample.deltaMs,
    kMagPeriodMs,
    accNorth,
    accEast,
    accDown,
    gyrNorth,
    gyrEast,
    gyrDown,
    magNorthUt,
    magEastUt,
    magDownUt);

  waitForNextLoopTick(loopStartUs);
}
