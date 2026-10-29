/*
    Copyright (c) 2026  Mikhail Grushinskiy

    carma_lever_arm-sim -- the off-CoG installation replay of the experimental
    CARMA filters behind tools/carma_lever_arm_study.py
    (docs/carma-imu-lever-arm.md).

    It is carma3-sim's A41 and tfg_carma3-sim's T41 (the same CarmaAdapter,
    the configuration the CARMA3 and TFG-CARMA3 sketches mirror) with one
    addition: W3D_IMU_LEVER_ARM_FILTER hands the installation offset of
    W3D_IMU_LEVER_ARM_M to one of the estimator's own lever-arm models instead
    of the simulator's input-stage model, which must then be "none":

      wrapper  SeaStateFusionFilter_CARMA3::setImuLeverArmBody, the correction
               the sketches' IMU_LEVER_ARM_*_M enable: the wrapper's input;
      mekf     the inner filter's set_imu_lever_arm_body (OU-III's placement),
               the same band-limited model, accelerometer update only;
      ou3      the same placement with OU-III's model
               (LeverArmModel::OuIII), its alpha smoothing
               W3D_IMU_LEVER_ARM_ALPHA_TAU_S (default 0.05 s).

    The band-limited models take their derivative corner from
    W3D_IMU_LEVER_ARM_CUTOFF_HZ (default 15 Hz, as for the input stage).

    W3D_CARMA_LEVER_VARIANT selects A41 (default) or T41.

    A separate binary so the gated replays keep their arithmetic: code added to
    their translation unit moves GCC's inlining, and with it the FMA
    contraction of the filter.  Nothing here is gated; every arm of the study
    runs in this binary, so its comparisons are like for like.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "carma3/SimAdapters.h"

using Eigen::Vector3f;

namespace {

using namespace carma3_sim;

// The filter-side lever arm of W3D_IMU_LEVER_ARM_FILTER, in the filter's body
// NED frame.
struct FilterLeverArm {
    std::string placement;  // wrapper | mekf | ou3
    Vector3f offset_ned = Vector3f::Zero();
    float cutoff_hz = 15.0f;
    std::optional<float> alpha_tau_s;
};

std::optional<FilterLeverArm> filter_lever_arm_from_env() {
    const char* flag = std::getenv("W3D_IMU_LEVER_ARM_FILTER");
    if (!flag || std::string(flag) == "0") return std::nullopt;
    const auto cfg = w3d_lever_arm_config_from_env();
    if (!cfg) throw std::invalid_argument("W3D_IMU_LEVER_ARM_FILTER needs W3D_IMU_LEVER_ARM_M");
    if (cfg->model != W3dLeverArmConfig::Model::None)
        throw std::invalid_argument("W3D_IMU_LEVER_ARM_FILTER needs W3D_IMU_LEVER_ARM_MODEL=none");
    FilterLeverArm out;
    out.placement = flag;
    if (out.placement != "wrapper" && out.placement != "mekf" && out.placement != "ou3")
        throw std::invalid_argument("W3D_IMU_LEVER_ARM_FILTER must be 0, wrapper, mekf or ou3");
    // The runner hands the filter zu_to_ned() of every body vector, so the
    // offset goes through the same map.
    out.offset_ned = zu_to_ned(cfg->offset_body_zu);
    out.cutoff_hz = cfg->derivative_cutoff_hz;
    if (const char* tau = std::getenv("W3D_IMU_LEVER_ARM_ALPHA_TAU_S")) {
        char* end = nullptr;
        const float v = std::strtof(tau, &end);
        if (end == tau || *end != '\0' || !std::isfinite(v) || v < 0.0f)
            throw std::invalid_argument("W3D_IMU_LEVER_ARM_ALPHA_TAU_S must be a non-negative number");
        out.alpha_tau_s = v;
    }
    return out;
}

template<class Core>
class LeverAdapter final : public IW3dFusionAdapter {
public:
    LeverAdapter(bool m, const Vector3f& a, const Vector3f& g, const Vector3f& s)
        : impl_(m, a, g, s, nullptr) {
        const auto lever = filter_lever_arm_from_env();
        if (!lever) return;
        auto& filter = impl_.fusion().raw();
        auto& mekf = filter.mekf();
        if (lever->placement == "wrapper") {
            filter.setImuLeverArmBody(lever->offset_ned);
            filter.setImuLeverArmCutoffHz(lever->cutoff_hz);
        } else {
            mekf.set_imu_lever_arm_body(lever->offset_ned);
            if (lever->placement == "ou3") {
                mekf.set_imu_lever_arm_model(Core::LeverArmModel::OuIII);
                if (lever->alpha_tau_s) mekf.set_alpha_smoothing_tau(*lever->alpha_tau_s);
            } else {
                mekf.set_imu_lever_arm_cutoff_hz(lever->cutoff_hz);
            }
        }
        std::cout << "IMU_LEVER_ARM_FILTER placement=" << lever->placement << " offset_ned_m=["
                  << lever->offset_ned.transpose() << "] cutoff_hz=" << lever->cutoff_hz << " alpha_tau_s="
                  << (lever->alpha_tau_s ? std::to_string(*lever->alpha_tau_s) : std::string("default"))
                  << "\n";
    }
    void updateMag(const Vector3f& b) override { impl_.updateMag(b); }
    void update(float dt, const Vector3f& g, const Vector3f& a, float t) override { impl_.update(dt, g, a, t); }
    FilterSnapshot snapshot() const override { return impl_.snapshot(); }

private:
    CarmaAdapter<4, false, 4, false, Core> impl_;
};

using LeverAdapter_A41 = LeverAdapter<ocean_imu::kalman::Kalman3D_Wave_CARMA3<float, 4, false, 4, false>>;
using LeverAdapter_T41 = LeverAdapter<ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3<float, 4>>;

constexpr float kUngated = 1.0e9f;
constexpr W3dFailureLimits kUngatedBars{
    .err_limit_percent_z_jonswap   = kUngated,
    .err_limit_percent_z_pmstokes  = kUngated,
    .err_limit_yaw_deg             = kUngated,
    .err_limit_percent_3d_jonswap  = kUngated,
    .err_limit_percent_3d_pmstokes = kUngated,
    .acc_z_bias_percent            = kUngated,
    .bias_3d_percent               = kUngated,
};

constexpr W3dSummaryLabels kLabels{.target = "RS_target", .applied = "kappa_z_applied"};

template<class A>
void run_file(const std::string& filename, float dt, bool with_mag, const W3dRandomSeeds& seeds,
              bool write_timeseries, float window, const std::string& family) {
    auto result = process_wave_file_for_tracker<A>(filename, dt, with_mag, add_noise, 25.0f,
                                                   "_fusion_" + family, "_fusion_" + family + "_nomag",
                                                   seeds, write_timeseries);
    if (!result) return;
    if (window > 0.0f) print_validation_metrics(*result, dt, window, family.c_str());
    print_summary_and_fail_if_needed(*result, dt, kUngatedBars, kLabels);
}

}  // namespace

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") with_mag = false;
        else if (arg == "--no-noise") add_noise = false;
        else if (arg == "--input" && i + 1 < argc) files.emplace_back(argv[++i]);
        else {
            std::cerr << "Usage: " << argv[0] << " [--nomag] [--no-noise] [--input PATH]...\n";
            return 2;
        }
    }
    W3dRandomSeeds seeds;
    try {
        seeds = w3d_random_seeds_from_env();
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    const char* variant_env = std::getenv("W3D_CARMA_LEVER_VARIANT");
    const std::string variant = variant_env ? variant_env : "A41";
    bool write_timeseries = false;
    if (const char* v = std::getenv("W3D_WRITE_TIMESERIES")) write_timeseries = std::string(v) != "0";
    float window = 900.0f;
    if (const char* v = std::getenv("W3D_VALIDATION_WINDOW_SEC")) window = static_cast<float>(std::atof(v));
    if (files.empty()) files = collect_wave_data_files(".");
    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", noise=" << (add_noise ? "true" : "false") << ", variant=" << variant << "\n";
    try {
        for (const auto& f : files) {
            const std::string family = "CARMA_LEVER_" + variant;
            if (variant == "A41") run_file<LeverAdapter_A41>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else if (variant == "T41") run_file<LeverAdapter_T41>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else throw std::invalid_argument("W3D_CARMA_LEVER_VARIANT must be A41 or T41");
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
