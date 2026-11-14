// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Gated default-seed 200 Hz replay of the experimental 27-state TFG shared-mode INS.
// Uses the same adapter, sensor stream and score window as the QMEKF replay.
#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#include "../common/SharedModesReplay.h"
#include "../common/TfgSharedRegressionBars.h"
int main(int argc, char* argv[]) {
    return ocean_imu::sim::shared_modes::run<ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>>(
        argc, argv, kTfgSharedRegressionBars, ocean_imu::sim::current_limits::tcsm, "TFG_CARMA_SHARED_MODES_");
}
