// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Use the standard shared-mode replay with the AtomS3R wrapper's mounting correction.
#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#include "../common/SharedModesReplay.h"

namespace replay = ocean_imu::sim::shared_modes;
using Core = ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>;

class InstalledSensor : public replay::CsmAdapter<Core, 0> {
public:
    InstalledSensor(bool mag, const Vector3f& a, const Vector3f& g, const Vector3f& s)
        : CsmAdapter(mag, a, g, s) {
        const char* placement = std::getenv("W3D_IMU_LEVER_ARM_FILTER");
        if (!placement) return;
        const auto installation = w3d_lever_arm_config_from_env();
        if (std::string(placement) != "wrapper" || !installation ||
            installation->model != W3dLeverArmConfig::Model::None)
            throw std::invalid_argument("wrapper correction requires an offset and simulator model none");
        filter().setImuLeverArmBody(zu_to_ned(installation->offset_body_zu));
        filter().setImuLeverArmCutoffHz(installation->derivative_cutoff_hz);
        std::cout << "IMU_LEVER_ARM_FILTER placement=wrapper\n";
    }
};

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: installation-study wave_data.csv\n";
        return 2;
    }
    try {
        return replay::run_file<InstalledSensor>(argv[1], 0.005f, true, w3d_random_seeds_from_env(),
            false, 900.0f, "TFG_CARMA_SHARED_MODES_CSM", replay::kUngatedBars, false) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
