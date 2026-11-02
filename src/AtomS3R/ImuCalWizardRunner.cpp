/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include "AtomS3R/ImuCalWizardRunner.h"

#include "AtomS3R/AtomS3R_ImuCalWizard.h"
#include <memory>
#include <new>

namespace atoms3r_ical {
bool runImuCalWizard(M5Ui& ui,ImuCalStoreNvs& store,ImuCalBlobV6& out_saved) {
  std::unique_ptr<ImuCalWizard> wizard(new(std::nothrow) ImuCalWizard(ui,store));
  if(!wizard) {ui.fail("CALIBRATION","Not enough memory");return false;}
  return wizard->runAndSave(out_saved);
}
} // namespace atoms3r_ical
