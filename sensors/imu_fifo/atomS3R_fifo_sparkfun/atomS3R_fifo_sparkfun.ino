#include <Wire.h>

#include "AtomS3R/AtomS3R_ImuCommon.h"
#include "SparkFun/ImuReader.h"
#include "SparkFun/MagReader.h"

using namespace atoms3r_ical;

constexpr uint8_t ATOMS3R_SENSOR_SDA_PIN = 45;
constexpr uint8_t ATOMS3R_SENSOR_SCL_PIN = 0;
constexpr uint32_t ATOMS3R_SENSOR_I2C_HZ = 400000;

constexpr float kImuOdrHz = kDefaultImuOdrHz;
constexpr uint8_t kImuAccelOdr = ImuReader::kDefaultAccelOdr;
constexpr uint8_t kImuGyroOdr = ImuReader::kDefaultGyroOdr;
constexpr uint8_t kMagOdr = MagReader::kDefaultBmm150Odr;
constexpr uint8_t kAuxOdr = MagReader::kDefaultAuxOdr;
constexpr uint32_t kLoopPeriodUs = loopPeriodUsFromHz(kImuOdrHz);

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
    waitForNextLoopTick(loopStartUs, kLoopPeriodUs);
    return;
  }

  const MagReader::Sample magSample =
    magAvailable ? magReader.readLatestSample(imuReader.sensor(), imuSample.timestampMs)
                 : MagReader::Sample{};

  float accNorth = 0.0f, accEast = 0.0f, accDown = 0.0f;
  ImuReader::sensorToNED(imuSample.accelX, imuSample.accelY, imuSample.accelZ, accNorth, accEast, accDown);

  float gyrNorth = 0.0f, gyrEast = 0.0f, gyrDown = 0.0f;
  ImuReader::sensorToNED(imuSample.gyroX, imuSample.gyroY, imuSample.gyroZ, gyrNorth, gyrEast, gyrDown);

  float magNorthUt = 0.0f, magEastUt = 0.0f, magDownUt = 0.0f;
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

  ImuLogRow row{};
  row.fifoLen = static_cast<int32_t>(imuSample.framesRead);
  row.fifoReq = -1;
  row.fifoGot = -1;

  row.dtImuMs = imuSample.imuDeltaMs;
  row.magValid = magSample.valid;
  row.magUpdated = magSample.updated;
  row.magAgeMs = magSample.valid ? magSample.deltaMs : -1.0f;
  row.magStepMs = MagReader::odrSettingToPeriodMs(kMagOdr);
  row.lastMagDtMs = magSample.updated ? magSample.deltaMs : -1.0f;

  row.accN = accNorth;
  row.accE = accEast;
  row.accD = accDown;

  row.gyrN = gyrNorth;
  row.gyrE = gyrEast;
  row.gyrD = gyrDown;

  row.magN = magNorthUt;
  row.magE = magEastUt;
  row.magD = magDownUt;

  printImuLogRow(row);
  waitForNextLoopTick(loopStartUs, kLoopPeriodUs);
}
