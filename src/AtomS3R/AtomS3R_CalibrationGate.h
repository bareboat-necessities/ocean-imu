#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Startup gate for AtomS3R sketches that must not run an estimator on raw
  sensor values: the sketch stops until the saved calibration applies in full.
*/

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"

namespace atoms3r_ical {

// The wizard entry point (runImuCalWizard), or nullptr in a build without it.
using CalibrationWizardFn = bool (*)(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV5& out_saved);

// Returns once calibrationUnusableReason(have_blob, runtime) is null, after
// logging "[CAL] applied: ...". Until then it reports the reason and the
// mismatching values on serial every 5 s, shows the reason on the display and,
// on a tap, runs `wizard`, adopts a saved calibration into blob / have_blob /
// runtime and checks again. Without a wizard it never returns.
//
// Defined in AtomS3R_CalibrationGate.cpp, outside the sketch object: in the
// large Arduino sketch object the Xtensa literal pool of this code landed
// beyond L32R reach (as for the OU-III prediction, Kalman3D_Wave_OU_III.cpp).
void requireUsableCalibration(M5Ui& ui, ImuCalStoreNvs& store, ImuCalBlobV5& blob, bool& have_blob,
                              RuntimeCals& runtime, CalibrationWizardFn wizard);

}  // namespace atoms3r_ical

// The wizard a sea-state INS sketch's gate offers: SEA_STATE_ENABLE_WIZARD,
// which those sketches define before their includes, selects runImuCalWizard;
// without it the gate offers none and no wizard is linked.
#if defined(SEA_STATE_ENABLE_WIZARD) && SEA_STATE_ENABLE_WIZARD
#include "AtomS3R/ImuCalWizardRunner.h"
#endif

namespace atoms3r_ical {

static constexpr CalibrationWizardFn kSeaStateCalibrationWizard =
#if defined(SEA_STATE_ENABLE_WIZARD) && SEA_STATE_ENABLE_WIZARD
    runImuCalWizard;
#else
    nullptr;
#endif

}  // namespace atoms3r_ical
