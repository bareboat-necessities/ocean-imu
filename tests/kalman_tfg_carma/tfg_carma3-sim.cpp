/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    tfg_carma3-sim -- the experimental TFG-CARMA3 INS (CARMA(4,1) on x, y and
    z on the two-frame Lie group, Kalman3D_Wave_TFG_CARMA3 behind the CARMA3
    wrapper; docs/tfg-carma3-research.md) on the pinned vessel-RAO dataset,
    through the shared W3dSimulationRunner.

    The adapter is carma3-sim's (tests/kalman_carma/carma3/SimAdapters.h),
    whose configuration is the one tests/kalman_ou_iii/kalman_ou_iii-sim.cpp
    gives the shipping OU-III, so this replay, carma3-sim's A41 and
    kalman_ou_iii-sim differ only in the estimator.  This simulator is the
    gated default-seed replay behind doc/kalman_tfg_carma; the paired
    multi-seed comparison is tools/tfg_carma3_benchmark.py.

    W3D_TFG_CARMA3_VARIANT selects T41 (default; full reset transport) or
    T41R0 (no reset transport).
*/

#define EIGEN_NON_ARDUINO

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../kalman_carma/carma3/SimAdapters.h"

namespace {

using namespace carma3_sim;

// Regression bars for the default variant T41 on the eight scored records
// (default seeds, magnetometer on), each at the measured worst value plus
// about 0.5 %, the margin kalman_ou_iii-sim, kalman_tfg-sim and carma3-sim
// use.
constexpr W3dFailureLimits kT41RegressionBars{
    .err_limit_percent_z_jonswap   = 3.072f,  // worst 3.05674 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 3.231f,  // worst 3.21467 (pmstokes H0.27)
    .err_limit_yaw_deg             = 0.7662f, // worst 0.76233 (jonswap H4.0)
    .err_limit_roll_deg            = 0.3569f, // worst 0.355132 (pmstokes H8.5)
    .err_limit_pitch_deg           = 0.2138f, // worst 0.212754 (pmstokes H8.5)
    .err_limit_percent_3d_jonswap  = 4.890f,  // worst 4.86593 (jonswap H8.5)
    .err_limit_percent_3d_pmstokes = 5.278f,  // worst 5.25201 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.206f,  // worst 4.18498 (pmstokes H8.5)
    .bias_3d_percent               = 82.53f,  // worst 82.1137 (pmstokes H8.5, accel)
    .gyro_bias_3d_percent          = 12.40f,  // worst 12.3347 (pmstokes H8.5, gyro)
};

// The no-transport variant prints the same summary but has no fitted bars.
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
    const char* variant_env = std::getenv("W3D_TFG_CARMA3_VARIANT");
    const std::string variant = variant_env ? variant_env : "T41";
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
        const bool gated = variant == "T41" && add_noise && with_mag && default_seeds;
        const W3dFailureLimits& bars = gated ? kT41RegressionBars : kUngatedBars;
        for (const auto& f : files) {
            const std::string family = "TFG_CARMA3_" + variant;
            if (variant == "T41") run_file<Adapter_T41>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "T41R0") run_file<Adapter_T41R0>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else throw std::invalid_argument("W3D_TFG_CARMA3_VARIANT must be T41 or T41R0");
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}
