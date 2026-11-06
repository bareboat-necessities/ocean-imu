/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    kalman_carma_shared_modes-sim -- the EXPERIMENTAL shared-modes QMEKF
    (src/kalman_carma_shared_modes, KalmanCarmaSharedModes<float> behind the
    CARMA3 wrapper) on the pinned vessel-RAO dataset, through the shared
    W3dSimulationRunner.  Same records, seeds, sensor model and adapter
    configuration as tests/kalman_carma/carma3-sim (carma3_sim::CarmaAdapter,
    the CARMA3 replay adapter, with the shared-modes core as its inner
    filter), so the two simulators differ only in the motion model.

    W3D_CSM_VARIANT selects CSM (default: estimated coupling + coupling
    confidence, the deployable configuration), CSM_NOCONF (confidence
    disabled: the research filter SHM:est) or CSM_OFF (coupling off: the
    CARMA(4,1) law on every axis).  The regression bars are fitted to CSM's
    own replays of these records and gate CSM only (default seeds, noise and
    magnetometer); the other variants are diagnostics.  The CARMA3 baseline
    on the same records and seeds is tests/kalman_carma/carma3-sim (A41); the
    paired Monte Carlo comparison is kalman_carma_shared_modes-bench.
*/

#define EIGEN_NON_ARDUINO

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "carma3/SimAdapters.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"

using Eigen::Vector3f;

namespace {

using namespace carma3_sim;
namespace csm = ocean_imu::carma_shared_modes;
using CsmCore = csm::KalmanCarmaSharedModes<float>;

// process_wave_file_for_tracker<> constructs the adapter from four arguments.
template<int Variant>   // 0 CSM, 1 CSM_NOCONF, 2 CSM_OFF
class CsmAdapter final : public IW3dFusionAdapter {
public:
    CsmAdapter(bool m, const Vector3f& a, const Vector3f& g, const Vector3f& s) : impl_(m, a, g, s, nullptr) {
        csm::SharedModesConfig<float> c;
        c.mode = Variant == 2 ? csm::CouplingMode::Off : csm::CouplingMode::Estimated;
        c.confidence.enabled = Variant == 0;
        c.confidence.accel_noise_std = a.x();
        impl_.fusion().raw().mekf().configureSharedModes(c);
    }
    void updateMag(const Vector3f& b) override { impl_.updateMag(b); }
    void update(float dt, const Vector3f& g, const Vector3f& a, float t) override { impl_.update(dt, g, a, t); }
    FilterSnapshot snapshot() const override { return impl_.snapshot(); }
private:
    CarmaAdapter<4, false, 4, false, CsmCore> impl_;
};

// Regression bars of the deployable variant CSM on the eight scored records
// (default seeds, magnetometer on), fitted to CSM's own replays: each is the
// worst CSM value plus about 0.5 % (the margin convention of carma3-sim and
// kalman_ou_iii-sim).  Not the CARMA3 baseline's bars.
constexpr W3dFailureLimits kCsmRegressionBars{
    .err_limit_percent_z_jonswap   = 3.244f,  // worst 3.22788 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 3.393f,  // worst 3.37567 (pmstokes H0.27)
    .err_limit_yaw_deg             = 0.7319f, // worst 0.728307 (jonswap H4.0)
    .err_limit_roll_deg            = 0.3587f, // worst 0.35693 (pmstokes H8.5)
    .err_limit_pitch_deg           = 0.2247f, // worst 0.223595 (pmstokes H8.5)
    .err_limit_percent_3d_jonswap  = 4.674f,  // worst 4.65064 (jonswap H0.27)
    .err_limit_percent_3d_pmstokes = 5.007f,  // worst 4.98241 (pmstokes H0.27)
    .acc_z_bias_percent            = 4.191f,  // worst 4.17061 (pmstokes H8.5)
    .bias_3d_percent               = 84.15f,  // worst 83.7269 (pmstokes H8.5, accel)
    .gyro_bias_3d_percent          = 12.71f,  // worst 12.6464 (pmstokes H8.5, gyro)
};

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

// FilterSnapshot::tuning_target/applied carry the r_S target and the
// vertical model's frequency scale kappa_z (see CarmaAdapter::snapshot).
constexpr W3dSummaryLabels kLabels{.target = "RS_target", .applied = "kappa_z_applied"};

template<class A>
void run_file(const std::string& filename, float dt, bool with_mag, const W3dRandomSeeds& seeds,
              bool write_timeseries, float window, const std::string& family, const W3dFailureLimits& bars) {
    auto result = process_wave_file_for_tracker<A>(filename, dt, with_mag, add_noise, 25.0f,
                                                   "_fusion_" + family, "_fusion_" + family + "_nomag",
                                                   seeds, write_timeseries);
    if (!result) return;
    if (window > 0.0f) print_validation_metrics(*result, dt, window, family.c_str());
    print_summary_and_fail_if_needed(*result, dt, bars, kLabels);
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
    const char* v = std::getenv("W3D_CSM_VARIANT");
    const std::string variant = v ? v : "CSM";
    bool write_timeseries = true;
    if (const char* w = std::getenv("W3D_WRITE_TIMESERIES")) write_timeseries = std::string(w) != "0";
    float window = 900.0f;
    if (const char* w = std::getenv("W3D_VALIDATION_WINDOW_SEC")) window = static_cast<float>(std::atof(w));
    if (files.empty()) files = collect_wave_data_files(".");
    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", noise=" << (add_noise ? "true" : "false") << ", variant=" << variant << "\n";
    try {
        // The bars are fitted to the default seeds with noise and magnetometer;
        // any other replay is a diagnostic, not a regression.
        const bool default_seeds = !std::getenv("W3D_SEED") && !std::getenv("W3D_IMU_SEED")
                                   && !std::getenv("W3D_INIT_SEED");
        const bool gated = variant == "CSM" && add_noise && with_mag && default_seeds;
        const W3dFailureLimits& bars = gated ? kCsmRegressionBars : kUngatedBars;
        for (const auto& f : files) {
            const std::string family = "CARMA_SHARED_MODES_" + variant;
            if (variant == "CSM") run_file<CsmAdapter<0>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "CSM_NOCONF") run_file<CsmAdapter<1>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "CSM_OFF") run_file<CsmAdapter<2>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else throw std::invalid_argument("W3D_CSM_VARIANT must be CSM, CSM_NOCONF or CSM_OFF");
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}
