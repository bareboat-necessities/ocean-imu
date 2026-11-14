// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Default-seed shared-modes QMEKF replay on the pinned vessel-RAO dataset.
#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
#include "../common/SharedModesReplay.h"
namespace {
constexpr auto kCsmRegressionBars = ocean_imu::sim::current_limits::interval::csm;

}
int main(int argc, char* argv[]) {
    return ocean_imu::sim::shared_modes::run<ocean_imu::carma_shared_modes::KalmanCarmaSharedModes<float>>(
        argc, argv, kCsmRegressionBars, ocean_imu::sim::current_limits::csm, "CARMA_SHARED_MODES_");
}
