/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include "AtomS3R/ImuCalWizardRunner.h"

#include <memory>
#include <new>

#include "AtomS3R/AtomS3R_ImuCalWizard.h"

namespace atoms3r_ical {

bool runImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV5& out_saved) {
  // Keep the large wizard object off the caller task stack.
  // This avoids stack pressure in sketches that invoke calibration from loop()
  // with smaller stack budgets than CompassAppBase.
  std::unique_ptr<ImuCalWizard> wizard(new (std::nothrow) ImuCalWizard(ui, store));
  if (!wizard) {
    ui.fail("CALIBRATION", "Not enough memory");
    return false;
  }
  return wizard->runAndSave(out_saved);
}

} // namespace atoms3r_ical
