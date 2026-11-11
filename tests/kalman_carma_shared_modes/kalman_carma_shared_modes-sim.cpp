// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Default-seed shared-modes QMEKF replay on the pinned vessel-RAO dataset.
#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
#include "../common/SharedModesReplay.h"
namespace {
constexpr W3dFailureLimits kCsmRegressionBars{
    .err_limit_percent_z_jonswap   = 2.467f,  // worst 2.45469 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 2.416f,  // worst 2.40321 (pmstokes H0.27)
    .err_limit_yaw_deg             = 0.7595f, // worst 0.755632 (pmstokes H0.27)
    .err_limit_roll_deg            = 0.3590f, // worst 0.357145 (pmstokes H8.5)
    .err_limit_pitch_deg           = 0.2231f, // worst 0.221974 (pmstokes H8.5)
    .err_limit_percent_3d_jonswap  = 4.433f,  // worst 4.41037 (jonswap H8.5)
    .err_limit_percent_3d_pmstokes = 4.904f,  // worst 4.87925 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.186f,  // worst 4.16434 (pmstokes H8.5)
    .bias_3d_percent               = 83.99f,  // worst 83.5627 (pmstokes H8.5, accel)
    .gyro_bias_3d_percent          = 12.71f,  // worst 12.6435 (pmstokes H8.5, gyro)
};

}
int main(int argc, char* argv[]) {
    return ocean_imu::sim::shared_modes::run<ocean_imu::carma_shared_modes::KalmanCarmaSharedModes<float>>(
        argc, argv, kCsmRegressionBars, "CARMA_SHARED_MODES_");
}
