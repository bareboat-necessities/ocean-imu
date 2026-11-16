// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Experimental tangent replay with independent, locked v1.4.0 case/sea bars.
// Synthetic paired-study acceptance criteria remain separate.
// Uses the same adapter, sensor stream and score window as the QMEKF replay.
#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanTangentCarmaSharedModes.h"
#include "../common/SharedModesReplay.h"
int main(int argc, char* argv[]) {
    return ocean_imu::sim::shared_modes::run<ocean_imu::carma_shared_modes::KalmanTangentCarmaSharedModes<float>>(
        argc, argv, ocean_imu::sim::current_limits::interval::tangent,
        ocean_imu::sim::current_limits::tangent,
        "TANGENT_CARMA_SHARED_MODES_");
}
