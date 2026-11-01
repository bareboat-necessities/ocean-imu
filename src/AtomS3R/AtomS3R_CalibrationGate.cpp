/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include "AtomS3R/AtomS3R_CalibrationGate.h"

namespace atoms3r_ical {
namespace {

void reportCalibration(const char* why, bool wizard, bool have_blob, const ImuCalBlobV5& blob,
                       const RuntimeCals& rc) {
  if (!why) {
    Serial.println("[CAL] applied: accelerometer, gyroscope, magnetometer");
    return;
  }
  Serial.printf("[CAL] UNUSABLE: %s. Stopped before the estimator%s.\n", why,
                wizard ? "; tap to calibrate" : "");
  Serial.printf("[CAL] blob=%d A=%d G=%d M=%d\n", static_cast<int>(have_blob), static_cast<int>(rc.acc.ok),
                static_cast<int>(rc.gyr.ok), static_cast<int>(rc.mag.ok));
  if (rc.accel_gravity_mismatch) {
    Serial.printf("[CAL] accel_g=%.7f, firmware g_cal_local=%.7f (-DATOMS3R_CALIBRATION_GRAVITY_MPS2)\n",
                  static_cast<double>(blob.accel_g), static_cast<double>(ImuCalCfg::g_cal_local));
  }
  if (rc.mag_source_mismatch) {
    Serial.printf("[CAL] mag fitted on source %u, firmware delivers source %u\n",
                  static_cast<unsigned>(blob.mag_source), static_cast<unsigned>(activeMagSource()));
  }
}

}  // namespace

void requireUsableCalibration(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV5& blob, bool& have_blob,
                              RuntimeCals& runtime, CalibrationWizardFn wizard) {
  const char* why = calibrationUnusableReason(have_blob, runtime);
  reportCalibration(why, wizard != nullptr, have_blob, blob, runtime);
  while (why) {
    ui.setReadRotation();
    ui.title("STOPPED");
    ui.line("Calibration unusable:");
    ui.line(why);
    ui.line("");
    if (wizard) ui.line("Tap: calibrate");

    uint32_t last_report_ms = millis();
    for (;;) {
      if (wizard) {
        Input::update();
        if (Input::tapPressed()) break;
      }
      if (millis() - last_report_ms >= 5000u) {
        last_report_ms = millis();
        reportCalibration(why, wizard != nullptr, have_blob, blob, runtime);
      }
      delay(10);
    }

    clearM5UnifiedImuCalibration();
    ImuCalBlobV5 saved{};
    if (wizard(ui, store, saved)) {
      blob = saved;
      have_blob = true;
      runtime.rebuildFromBlob(blob);
    }
    why = calibrationUnusableReason(have_blob, runtime);
    reportCalibration(why, wizard != nullptr, have_blob, blob, runtime);
  }
}

}  // namespace atoms3r_ical
