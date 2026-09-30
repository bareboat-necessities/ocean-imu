/*
  Copyright (c) 2026 Mikhail Grushinskiy

  AtomS3R IMU calibration implementation unit.

  Magnetometer/IMU bring-up is intentionally owned by M5Unified, matching the
  v2.3.3 shipping path. Do not add a parallel BMI270 AUX/BMM150 initialization
  sequence here: it creates boot-dependent sensor state.
*/

#include "AtomS3R/AtomS3R_ImuCal.h"
