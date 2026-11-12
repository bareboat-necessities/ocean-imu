#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// One shared-modes adapter and replay protocol for QMEKF and TFG geometry.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "../kalman_carma/carma3/SimAdapters.h"
#include "kalman_carma_shared_modes/SharedModesEstimator.h"
namespace ocean_imu::sim::shared_modes {
using namespace carma3_sim;
namespace csm = ocean_imu::carma_shared_modes;
// process_wave_file_for_tracker<> constructs the adapter from four arguments.
template<class Core, int Variant>   // 0 CSM, 1 CSM_NOCONF, 2 CSM_OFF
class CsmAdapter final : public IW3dFusionAdapter {
public:
    CsmAdapter(bool m, const Vector3f& a, const Vector3f& g, const Vector3f& s) : impl_(m, a, g, s, nullptr) {
        csm::SharedModesConfig<float> c;
        c.mode = Variant == 2 ? csm::CouplingMode::Off : csm::CouplingMode::Estimated;
        c.confidence.enabled = Variant == 0;
        c.confidence.accel_noise_std = a.x();
        impl_.fusion().raw().mekf().configureSharedModes(c);
        // Evidence excitation (the filter's default) unless W3D_CARMA3_EXCITATION
        // selects another mode for a diagnostic replay.
        if (!std::getenv("W3D_CARMA3_EXCITATION"))
            impl_.fusion().raw().setCarmaExcitationMode(ocean_imu::carma::CarmaExcitationMode::Evidence);
        // Controlled sensor-path ablations for the article's machinery study.
        // Unset variables retain the deployed wrapper configuration.
        auto& filter = impl_.fusion().raw();
        // Explicit sensor contract for diagnostic replays. Keep the original
        // interval-rate regression protocol and its thresholds unchanged.
        const char* gyro_input = std::getenv("W3D_CSM_GYRO_INPUT");
        const std::string gyro_mode = gyro_input ? gyro_input : "interval";
        if (gyro_mode == "instantaneous") {
            using Mode = ocean_imu::kalman::filter_detail::GyroSampleMode;
            filter.setGyroSampleMode(Mode::Instantaneous, Mode::Instantaneous);
            filter.mekf().set_gyro_sample_period(1.0f / 200.0f);
        } else if (gyro_mode != "interval") {
            throw std::invalid_argument("W3D_CSM_GYRO_INPUT must be interval or instantaneous");
        }
        std::cout << "CSM_GYRO_INPUT convention=" << gyro_mode << '\n';
        if (const char* hz = std::getenv("W3D_CSM_GUARD_HZ"))
            filter.setAccelVibrationGuard(std::stof(hz));
        if (const char* gain = std::getenv("W3D_CSM_GUARD_RACC_GAIN"))
            filter.setAccelVibrationRaccGain(std::stof(gain));
        std::cout << "CSM_GUARD cutoff_hz=" << filter.accelVibrationGuardCutoffHz()
                  << " racc_gain=" << filter.accelVibrationRaccGain() << '\n';
    }
    void updateMag(const Vector3f& b) override { impl_.updateMag(b); }
    void update(float dt, const Vector3f& g, const Vector3f& a, float t) override { impl_.update(dt, g, a, t); }
    FilterSnapshot snapshot() const override { return impl_.snapshot(); }
private:
    CarmaAdapter<4, false, 4, false, Core> impl_;
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
bool run_file(const std::string& filename, float dt, bool with_mag, const W3dRandomSeeds& seeds,
              bool write_timeseries, float window, const std::string& family, const W3dFailureLimits& bars, bool gated) {
    auto result = process_wave_file_for_tracker<A>(filename, dt, with_mag, add_noise, 25.0f,
                                                   "_fusion_" + family, "_fusion_" + family + "_nomag",
                                                   seeds, write_timeseries);
    if (!result) return false;
    if (gated && result->errs_z.size() < 240000)
        throw std::runtime_error("regression requires a full 1200 s, 200 Hz record: " + filename);
    for (const auto* values : {&result->errs_x, &result->errs_y, &result->errs_z,
                              &result->errs_roll, &result->errs_pitch, &result->errs_yaw,
                              &result->accb_err_x, &result->accb_err_y, &result->accb_err_z,
                              &result->gyrb_err_x, &result->gyrb_err_y, &result->gyrb_err_z}) {
        if (!std::all_of(values->begin(), values->end(), [](float x) { return std::isfinite(x); }))
            throw std::runtime_error("nonfinite estimation error in " + filename);
    }
    if (window > 0.0f) print_validation_metrics(*result, dt, window, family.c_str());
    print_summary_and_fail_if_needed(*result, dt, bars, kLabels);
    return true;
}

template<class Core>
int run(int argc, char* argv[], const W3dFailureLimits& regression_bars, const char* family_prefix) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    bool require_gates = false;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--regression") require_gates = true;
        else if (arg == "--nomag") with_mag = false;
        else if (arg == "--no-noise") add_noise = false;
        else if (arg == "--input" && i + 1 < argc) files.emplace_back(argv[++i]);
        else {
            std::cerr << "Usage: " << argv[0] << " [--regression] [--nomag] [--no-noise] [--input PATH]...\n";
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
    const bool explicit_files = !files.empty();
    if (files.empty()) files = collect_wave_data_files(".");
    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", noise=" << (add_noise ? "true" : "false") << ", variant=" << variant << "\n";
    try {
        // The bars are fitted to the default seeds with noise and magnetometer;
        // any other replay is a diagnostic, not a regression.
        const bool default_seeds = !std::getenv("W3D_SEED") && !std::getenv("W3D_IMU_SEED")
                                   && !std::getenv("W3D_INIT_SEED");
        const bool gated = variant == "CSM" && add_noise && with_mag && default_seeds &&
                           !std::getenv("W3D_CARMA3_EXCITATION") && !std::getenv("W3D_CSM_GYRO_INPUT");
        if (require_gates && !gated)
            throw std::invalid_argument("--regression requires default variant, noise, magnetometer, seeds, excitation and gyro contract");
        if (require_gates) {
            for (const char* key : {"W3D_ENGINE_RPM", "W3D_CSM_GUARD_HZ", "W3D_CSM_GUARD_RACC_GAIN",
                                    "W3D_IMU_LEVER_ARM_M", "W3D_IMU_LEVER_ARM_MODEL",
                                    "W3D_TUNING_MODE", "W3D_AW_COV_SYNC", "W3D_ALL_WAVE_TYPES"})
                if (std::getenv(key)) throw std::invalid_argument(std::string("--regression forbids ") + key);
            if (window != 900.0f) throw std::invalid_argument("--regression requires a 900 s score window");
        }
        std::cout << "REGRESSION_PROTOCOL: GATED=" << (gated ? 1 : 0) << " FAMILY=" << family_prefix << '\n';
        std::size_t scored = 0;
        const W3dFailureLimits& bars = gated ? regression_bars : kUngatedBars;
        for (const auto& f : files) {
            const std::string family = std::string(family_prefix) + variant;
            bool processed = false;
            if (variant == "CSM") processed = run_file<CsmAdapter<Core, 0>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars, gated);
            else if (variant == "CSM_NOCONF") processed = run_file<CsmAdapter<Core, 1>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars, gated);
            else if (variant == "CSM_OFF") processed = run_file<CsmAdapter<Core, 2>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars, gated);
            else throw std::invalid_argument("W3D_CSM_VARIANT must be CSM, CSM_NOCONF or CSM_OFF");
            if (processed) ++scored;
            else if (explicit_files) throw std::runtime_error("input is not a scored wave record: " + f);
        }
        if (scored == 0) throw std::runtime_error("no scored wave records found");
        if (require_gates && !explicit_files && scored != 8)
            throw std::runtime_error("default regression requires all eight JONSWAP/PM-Stokes records");
        std::cout << "REPLAY_COMPLETE: records=" << scored << " gated=" << (gated ? 1 : 0) << '\n';
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    if (w3d_any_quality_gate_failed()) return 1;
    return 0;
}

}  // namespace ocean_imu::sim::shared_modes
