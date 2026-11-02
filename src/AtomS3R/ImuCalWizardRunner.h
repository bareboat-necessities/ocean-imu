#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"

namespace atoms3r_ical {

// The out-of-line entry keeps fitting/GUI instantiations out of estimator sketches.
bool runImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV6& out_saved);

} // namespace atoms3r_ical
