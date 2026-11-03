#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"

namespace atoms3r_ical {

// Runs the calibration wizard and saves to NVS; true only when a calibration
// was saved and read back (out_saved). Defined in ImuCalWizardRunner.cpp,
// outside the sketch object: compiled into a large estimator sketch object,
// the wizard pushed that object's Xtensa literal pools beyond L32R reach (as
// for AtomS3R_CalibrationGate.cpp and Kalman3D_Wave_OU_III.cpp).
bool runImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV5& out_saved);

} // namespace atoms3r_ical
