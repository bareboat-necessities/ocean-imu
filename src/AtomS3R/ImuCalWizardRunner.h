#pragma once

#include <Wire.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_ImuCalWizard.h"
#include "AtomS3R/AtomS3R_M5Ui.h"

namespace atoms3r_ical {

inline bool runImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV2& out_saved) {
  ImuCalWizard wizard(ui, store);
  return wizard.runAndSave(out_saved);
}

} // namespace atoms3r_ical
