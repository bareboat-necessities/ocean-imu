// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Gated default-seed 200 Hz replay of the experimental 27-state TFG shared-mode INS.
// Uses the same adapter, sensor stream and score window as the QMEKF replay.
#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#include "../common/SharedModesReplay.h"
namespace {
// Worst of the eight default-seed 1200 s records, plus 0.5%, rounded upward.
// Baseline and input provenance: reports/results/tfg_shared_modes/deterministic_regression.json.
constexpr W3dFailureLimits kTfgSharedRegressionBars{
    .err_limit_percent_z_jonswap   = 2.4669f, // measured 2.4545393
    .err_limit_percent_z_pmstokes  = 2.4152f, // measured 2.4031157
    .err_limit_yaw_deg             = 0.7631f, // measured 0.7592051
    .err_limit_roll_deg            = 0.3585f, // measured 0.35668024
    .err_limit_pitch_deg           = 0.2188f, // measured 0.21763369
    .err_limit_percent_3d_jonswap  = 4.4354f, // measured 4.4132538
    .err_limit_percent_3d_pmstokes = 4.9068f, // measured 4.8822913
    .acc_z_bias_percent            = 4.1787f, // measured 4.1579
    .bias_3d_percent               = 83.4106f, // measured 82.9956
    .gyro_bias_3d_percent          = 12.3460f, // measured 12.2845
};

}
int main(int argc, char* argv[]) {
    return ocean_imu::sim::shared_modes::run<ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>>(
        argc, argv, kTfgSharedRegressionBars, ocean_imu::sim::current_limits::tcsm, "TFG_CARMA_SHARED_MODES_");
}
