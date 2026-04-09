#pragma once

#include <Wire.h>

#include <memory>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_ImuCalWizard.h"
#include "AtomS3R/AtomS3R_M5Ui.h"

namespace atoms3r_ical {

inline bool runImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV2& out_saved) {
  // Keep the large wizard object off the caller task stack.
  // This avoids stack pressure in sketches that invoke calibration from loop()
  // with smaller stack budgets than CompassAppBase.
  std::unique_ptr<ImuCalWizard> wizard(new ImuCalWizard(ui, store));
  if (!wizard) return false;
  return wizard->runAndSave(out_saved);
}

} // namespace atoms3r_ical
