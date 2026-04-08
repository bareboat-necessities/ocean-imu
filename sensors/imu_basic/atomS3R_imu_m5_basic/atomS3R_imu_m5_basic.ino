#include <Arduino.h>
#include <M5Unified.h>
#include <ArduinoOceanImu.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_ImuCommon.h"

using namespace atoms3r_ical;

constexpr float kImuOdrHz = kDefaultImuOdrHz;
constexpr uint32_t kLoopPeriodUs = loopPeriodUsFromHz(kImuOdrHz);

ImuCalStoreNvs calStore;
ImuCalBlobV2 calBlob;
RuntimeCals runtimeCals;
bool hasSavedCalibration = false;

#ifndef ATOMS3R_IMU_MASK_MAG
  #if defined(M5IMU_UPDATE_MAG)
    #define ATOMS3R_IMU_MASK_MAG M5IMU_UPDATE_MAG
  #elif defined(IMU_UPDATE_MAG)
    #define ATOMS3R_IMU_MASK_MAG IMU_UPDATE_MAG
  #else
    #define ATOMS3R_IMU_MASK_MAG (1u << 2)
  #endif
#endif

void setup()
{
  Serial.begin(115200);
  delay(1000);
  Serial.println("M5Stack AtomS3R IMU basic demo (M5Unified IMU API path)");

  auto cfg = M5.config();
  M5.begin(cfg);
  M5.Imu.begin();

  clearM5UnifiedImuCalibration();

  hasSavedCalibration = calStore.load(calBlob);
  if (hasSavedCalibration)
  {
    runtimeCals.rebuildFromBlob(calBlob);
    Serial.println("Loaded IMU calibration from NVS; magnetometer output will include calibrated values.");
  }
  else
  {
    runtimeCals.mag.ok = true;
    runtimeCals.mag.A = Matrix3f::Identity();
    runtimeCals.mag.b = Vector3f::Zero();
    Serial.println("No IMU calibration in NVS; assuming zero mag calibration (A=I, b=0).");
  }
}

void loop()
{
  const uint32_t loopStartUs = micros();
  static ImuLoopLogState logState{};
  static uint32_t lastSampleUs = 0u;
  static uint32_t lastMagSampleUs = 0u;
  static uint32_t prevMagSampleUs = 0u;

  const uint32_t sampleUs = micros();
  const uint32_t updateMask = M5.Imu.update();

  ImuSample sample{};
  if (!readImuMapped(M5.Imu, updateMask, sampleUs, sample))
  {
    logState.markReadFailure();
    const uint32_t nowUs = micros();
    if (logState.shouldLogNoSample(nowUs))
    {
      Serial.printf("Waiting for IMU sample (failures=%lu): update_mask=0x%08lx\n",
                    static_cast<unsigned long>(logState.consecutiveReadFailures),
                    static_cast<unsigned long>(updateMask));
    }
    waitForNextLoopTick(loopStartUs, kLoopPeriodUs);
    return;
  }
  logState.markReadSuccess();

  if ((updateMask & ATOMS3R_IMU_MASK_MAG) != 0u)
  {
    prevMagSampleUs = lastMagSampleUs;
    lastMagSampleUs = sampleUs;
  }

  ImuLogRow row{};
  row.fifoLen = -1;
  row.fifoReq = -1;
  row.fifoGot = -1;

  row.dtImuMs = (lastSampleUs != 0u) ? usToMs(static_cast<uint32_t>(sampleUs - lastSampleUs)) : -1.0f;
  row.magValid = (lastMagSampleUs != 0u);
  row.magUpdated = ((updateMask & ATOMS3R_IMU_MASK_MAG) != 0u);
  row.magAgeMs = row.magValid ? usToMs(static_cast<uint32_t>(sampleUs - lastMagSampleUs)) : -1.0f;
  row.magStepMs = (lastMagSampleUs != 0u && prevMagSampleUs != 0u)
    ? usToMs(static_cast<uint32_t>(lastMagSampleUs - prevMagSampleUs))
    : -1.0f;
  row.lastMagDtMs = row.magStepMs;

  row.accN = sample.a.x();
  row.accE = sample.a.y();
  row.accD = sample.a.z();

  row.gyrN = sample.w.x() * RAD_TO_DEG;
  row.gyrE = sample.w.y() * RAD_TO_DEG;
  row.gyrD = sample.w.z() * RAD_TO_DEG;

  row.magN = sample.m.x();
  row.magE = sample.m.y();
  row.magD = sample.m.z();

  const Vector3f magCal = runtimeCals.applyMag(sample.m);
  row.magCalN = magCal.x();
  row.magCalE = magCal.y();
  row.magCalD = magCal.z();

  const uint32_t nowUs = micros();
  if (logState.shouldLogSample(nowUs))
  {
    printImuLogRow(row);
  }

  lastSampleUs = sampleUs;
  waitForNextLoopTick(loopStartUs, kLoopPeriodUs);
}
