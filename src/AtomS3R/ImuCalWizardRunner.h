#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include <Wire.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"

namespace atoms3r_ical {

bool runImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV5& out_saved);

} // namespace atoms3r_ical
