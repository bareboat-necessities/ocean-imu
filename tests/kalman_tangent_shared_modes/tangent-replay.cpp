// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Experimental tangent replay against the UNCHANGED TFG quality bars.
// A failure is retained as an experimental negative result; never widen bars.
// Uses the same adapter, sensor stream and score window as the QMEKF replay.
#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanTangentCarmaSharedModes.h"
#include "../common/SharedModesReplay.h"
#include "../common/TfgSharedRegressionBars.h"
int main(int argc, char* argv[]) {
    return ocean_imu::sim::shared_modes::run<ocean_imu::carma_shared_modes::KalmanTangentCarmaSharedModes<float>>(
        argc, argv, kTfgSharedRegressionBars, "TANGENT_CARMA_SHARED_MODES_");
}
