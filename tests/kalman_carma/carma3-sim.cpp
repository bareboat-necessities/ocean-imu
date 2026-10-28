/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    carma3-sim -- the experimental three-axis CARMA INS (carma3/) on the
    pinned vessel-RAO dataset, through the shared W3dSimulationRunner.  The
    dataset's horizontal truth comes from its own vessel RAOs, independent of
    the carma3_bench truth model, which makes it the hull-independent check of
    the horizontal priors.

    The adapter configuration is the one tests/kalman_ou_iii/kalman_ou_iii-sim.cpp
    gives the shipping OU-III (sensor-sigma rescales, magnetometer delay, the
    dataset hull's direction RAO and low-wave Racc weighting, tuner and
    clamps; OU-III's periodic a_w covariance maintenance has no CARMA
    counterpart), so the two simulators differ only in the translational
    motion model.  Run the unchanged kalman_ou_iii-sim with the same W3D_SEED
    for the paired results.

    W3D_CARMA3_VARIANT selects A41 (default; CARMA(4,1) on all axes), A52,
    M45 (CARMA(4,1) x/y, (5,2) z), M54 (CARMA(5,2) x/y, (4,1) z), A52S (with
    S_x, S_y), A41S (OU-III r_S carried over, negative control), H41
    (horizontal shape fitted to the carma3_bench sway RAO), or T41 / T41R0:
    A41's prior and wrapper on the two-frame Lie group
    (Kalman3D_Wave_TFG_CARMA3, full reset / no reset transport;
    docs/tfg-carma3-research.md).
*/

#define EIGEN_NON_ARDUINO

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "carma3/SimAdapters.h"

using Eigen::Vector3f;

namespace {

using namespace carma3_sim;

// process_wave_file_for_tracker<> constructs the adapter from four arguments;
// these bind the variant.
#define CARMA_ADAPTER(NAME, NHC, SH, NZC, SZ, HSHAPE)                                 \
    class NAME final : public IW3dFusionAdapter {                                       \
    public:                                                                             \
        NAME(bool m, const Vector3f& a, const Vector3f& g, const Vector3f& s)           \
            : impl_(m, a, g, s, HSHAPE) {}                                              \
        void updateMag(const Vector3f& b) override { impl_.updateMag(b); }              \
        void update(float dt, const Vector3f& g, const Vector3f& a, float t) override { \
            impl_.update(dt, g, a, t);                                                  \
        }                                                                               \
        FilterSnapshot snapshot() const override { return impl_.snapshot(); }          \
    private:                                                                            \
        CarmaAdapter<NHC, SH, NZC, SZ> impl_;                                           \
    };

CARMA_ADAPTER(Adapter_A41, 4, false, 4, false, nullptr)
CARMA_ADAPTER(Adapter_A52, 5, false, 5, false, nullptr)
CARMA_ADAPTER(Adapter_M45, 4, false, 5, false, nullptr)
CARMA_ADAPTER(Adapter_M54, 5, false, 4, false, nullptr)
CARMA_ADAPTER(Adapter_A52S, 5, true, 5, false, nullptr)
CARMA_ADAPTER(Adapter_A41S, 4, true, 4, false, nullptr)
CARMA_ADAPTER(Adapter_H41, 4, false, 4, false, &carma3_bench::CARMA41_SWAY_SHAPE)

// Regression bars for the default variant A41 on the eight scored records
// (default seeds, magnetometer on), each at the measured worst value plus
// about 0.5 %, the margin kalman_ou_iii-sim and kalman_tfg-sim use.
constexpr W3dFailureLimits kA41RegressionBars{
    .err_limit_percent_z_jonswap   = 3.073f,  // worst 3.0571 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 3.231f,  // worst 3.21488 (pmstokes H0.27)
    .err_limit_yaw_deg             = 0.7575f, // worst 0.753738 (jonswap H4.0)
    .err_limit_roll_deg            = 0.3565f, // worst 0.354707 (pmstokes H8.5)
    .err_limit_pitch_deg           = 0.2137f, // worst 0.212555 (pmstokes H8.5)
    .err_limit_percent_3d_jonswap  = 4.888f,  // worst 4.86298 (jonswap H8.5)
    .err_limit_percent_3d_pmstokes = 5.276f,  // worst 5.24895 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.209f,  // worst 4.18715 (pmstokes H8.5)
    .bias_3d_percent               = 82.43f,  // worst 82.0115 (pmstokes H8.5, accel)
    .gyro_bias_3d_percent          = 12.78f,  // worst 12.7086 (pmstokes H8.5, gyro)
};

// The research variants print the same summary but have no fitted bars.
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
              bool write_timeseries, float window, const std::string& family,
              const W3dFailureLimits& bars) {
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
    const char* variant_env = std::getenv("W3D_CARMA3_VARIANT");
    const std::string variant = variant_env ? variant_env : "A41";
    bool write_timeseries = true;
    if (const char* v = std::getenv("W3D_WRITE_TIMESERIES")) write_timeseries = std::string(v) != "0";
    float window = 900.0f;
    if (const char* v = std::getenv("W3D_VALIDATION_WINDOW_SEC")) window = static_cast<float>(std::atof(v));
    if (files.empty()) files = collect_wave_data_files(".");
    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", noise=" << (add_noise ? "true" : "false") << ", variant=" << variant << "\n";
    try {
        // The bars are fitted to the default seeds with noise and magnetometer;
        // any other replay is a diagnostic, not a regression.
        const bool default_seeds = !std::getenv("W3D_SEED") && !std::getenv("W3D_IMU_SEED")
                                   && !std::getenv("W3D_INIT_SEED");
        const bool gated = variant == "A41" && add_noise && with_mag && default_seeds;
        const W3dFailureLimits& bars = gated ? kA41RegressionBars : kUngatedBars;
        for (const auto& f : files) {
            const std::string family = "CARMA3_" + variant;
            if (variant == "A41") run_file<Adapter_A41>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "A52") run_file<Adapter_A52>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "M45") run_file<Adapter_M45>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "M54") run_file<Adapter_M54>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "A52S") run_file<Adapter_A52S>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "A41S") run_file<Adapter_A41S>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "H41") run_file<Adapter_H41>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "T41") run_file<Adapter_T41>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "T41R0") run_file<Adapter_T41R0>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else throw std::invalid_argument("W3D_CARMA3_VARIANT must be A41, A52, M45, M54, A52S, A41S, H41, T41 or T41R0");
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}
