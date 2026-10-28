#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
  Released under the MIT License

  SeaStateFusionFilter_TFG_CARMA3 -- EXPERIMENTAL research wrapper around
  Kalman3D_Wave_TFG_CARMA3 (three-axis CARMA prior on the right-invariant
  two-frame Lie group).  Not a deployed estimator; SeaStateFusionFilter_OU_III
  stays the shipping filter.

  It is SeaStateFusionFilter_CARMA3 with the inner filter replaced, nothing
  else: the front end, proxy startup, staged tune commit, tilt watchdog,
  magnetometer-service bias gate, vibration Racc inflation and the CARMA
  amplitude/frequency mapping are the same code, so a paired comparison of the
  two wrappers isolates the error geometry (docs/tfg-carma3-research.md).
  CARMA(NC, NC-3) on all three axes, no S state.
*/

#include "kalman_carma/Kalman3D_Wave_TFG_CARMA3.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"

namespace ocean_imu::carma {

template<TrackerType trackerT, int NC = 4>
using SeaStateFusionFilter_TFG_CARMA3 =
    SeaStateFusionFilter_CARMA3<trackerT, NC, false, NC, false, kalman::Kalman3D_Wave_TFG_CARMA3<float, NC>>;

template<TrackerType trackerT, int NC = 4>
using SeaStateFusion_TFG_CARMA3 =
    SeaStateFusion_CARMA3<trackerT, NC, false, NC, false, kalman::Kalman3D_Wave_TFG_CARMA3<float, NC>>;

}  // namespace ocean_imu::carma
