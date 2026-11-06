/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    kalman_carma_shared_modes-sim -- simulator of the EXPERIMENTAL standalone
    shared-modes QMEKF (src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h).

    Two modes in one executable:

    1. Dataset replay (default; no --scenario): the pinned vessel-RAO records
       of oceanography-waves-lib (the same files, seeds, sensor model and
       W3dSimulationRunner as tests/kalman_carma/carma3-sim and the OU-III
       replay), through the adapter configuration carma3-sim gives the CARMA3
       baseline.  W3D_CSM_VARIANT selects
         CSM        deployable default (estimated coupling + confidence)
         CSM_NOCONF confidence disabled (the research filter SHM:est)
         CSM_OFF    coupling off (CARMA(4,1) law on every axis)
         A41        the unchanged CARMA3 baseline (paired comparison)
       The regression bars apply to CSM with the default seeds, noise and
       magnetometer.  Options: --nomag, --no-noise, --input PATH (repeatable).

    2. Monte Carlo (any of --scenario, --header, --level): carma3_bench's
       engine -- synthetic directional seas through the benchmark hull,
       seeded IMU / magnetometer streams, paired filters, the same metrics --
       with the CSM filters of CsmBenchAdapter.h registered in front of the
       bench's own (A41, SHM:est, ...).  Example:
         kalman_carma_shared_modes-sim --header --extended --coupling --level full \
             --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 --filters A41,SHM:est,CSM
       tools/carma_shared_modes_validation.py drives both modes.
*/

#define EIGEN_NON_ARDUINO

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "carma3/SimAdapters.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"

int carma3_bench_main(int argc, char** argv);   // tests/kalman_carma/carma3_bench.cpp

namespace {

using namespace carma3_sim;
namespace csm = ocean_imu::carma_shared_modes;
using Eigen::Quaternionf;
using Eigen::Vector3f;

// The CARMA3 replay configuration of carma3-sim (CarmaAdapter), on the
// standalone filter.
template<int Variant>   // 0 CSM, 1 CSM_NOCONF, 2 CSM_OFF
class CsmSimAdapter final : public IW3dFusionAdapter {
public:
    CsmSimAdapter(bool with_mag, const Vector3f& sigma_a_init, const Vector3f& sigma_g, const Vector3f& sigma_m) {
        csm::KalmanCarmaSharedModes::Config c;
        c.fusion.with_mag = with_mag;
        c.fusion.sigma_a = sigma_a_init * CarmaAdapter<4, false, 4, false>::SIGMA_A_RESCALE;
        c.fusion.sigma_g = sigma_g * CarmaAdapter<4, false, 4, false>::SIGMA_G_RESCALE;
        c.fusion.sigma_m = sigma_m * CarmaAdapter<4, false, 4, false>::SIGMA_M_RESCALE;
        c.fusion.mag_delay_sec = MAG_DELAY_SEC;
        c.mode = Variant == 2 ? csm::CouplingMode::Off : csm::CouplingMode::Estimated;
        c.confidence.enabled = Variant == 0;
        c.confidence.accel_noise_std = sigma_a_init.x();
        f_.begin(c);
        auto& filter = f_.fusion().raw();
        wave_direction::VesselRaoEqualizer::Config direction_rao;
        direction_rao.enabled = true;
        direction_rao.horizontal_x_tau_s = 1.0;
        direction_rao.horizontal_y_tau_s = 0.7;
        filter.setDirectionRao(direction_rao);
        wave_direction::VesselRaoNoiseWeighting low_wave_noise;
        low_wave_noise.response = direction_rao;
        low_wave_noise.max_std_scale = 4.0f;
        filter.setLowWaveNoiseWeighting(low_wave_noise);
    }
    void updateMag(const Vector3f& b) override { f_.updateMag(b); }
    void update(float dt, const Vector3f& g, const Vector3f& a, float t) override { f_.update(dt, g, a, t); }
    FilterSnapshot snapshot() const override {
        const auto& filter = f_.fusion().raw();
        const auto& d = filter.dir();
        FilterSnapshot s;
        s.disp_est_zu = ned_to_zu(f_.position());
        s.vel_est_zu = ned_to_zu(f_.velocity());
        s.acc_est_zu = ned_to_zu(f_.waveAcceleration());
        const Quaternionf q_bw_ned = f_.attitude().normalized();
        float roll_deg = 0.0f, pitch_deg = 0.0f, yaw_deg = 0.0f;
        quat_to_euler_nautical(q_bw_ned, roll_deg, pitch_deg, yaw_deg);
        s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, wrapDeg(yaw_deg));
        s.acc_bias_est_ned = f_.accelBias();
        s.gyro_bias_est_ned = f_.gyroBias();
        s.mag_bias_est_ned_uT = f_.fusion().magHardIronBodyUT();
        s.tau_target = filter.getTauTarget();
        s.sigma_target = filter.getSigmaTarget();
        s.tuning_target = f_.couplingConfidence();          // coupling confidence c
        s.tau_applied = filter.getTauApplied();
        s.sigma_applied = filter.getSigmaApplied();
        s.tuning_applied = filter.getVerticalKappa();
        s.freq_hz = filter.getFreqHz();
        s.wave_period_sec = filter.getWavePeriodSec();
        s.period_sec = filter.getPeriodSec();
        s.accel_variance = filter.getAccelVariance();
        s.displacement_scale_m = filter.getDisplacementScale();
        s.velocity_scale_mps = filter.getVerticalSpeedEnvelopeMps(true);
        s.direction.phase = d.getPhase();
        s.direction.direction_deg = d.getAxisDegrees();
        s.direction.apparent_to_deg = filter.getApparentWaveDirectionToDeg();
        s.direction.apparent_from_deg = filter.getApparentWaveDirectionFromDeg();
        s.direction.sense_coherence = filter.getDirSenseCoherence();
        s.direction.direction_deg_generator_signed = dirDegGeneratorSignedFromVec(d.getAxis());
        s.direction.uncertainty_deg = d.getAxisUncertaintyDegrees();
        s.direction.confidence = d.getLastStableConfidence();
        s.direction.amplitude = d.getAmplitude();
        s.direction.direction_vec = d.getAxis();
        s.direction.filtered_signal = d.getFilteredSignal();
        return s;
    }

private:
    mutable csm::KalmanCarmaSharedModes f_;
};

class A41Adapter final : public IW3dFusionAdapter {
public:
    A41Adapter(bool m, const Vector3f& a, const Vector3f& g, const Vector3f& s) : impl_(m, a, g, s, nullptr) {}
    void updateMag(const Vector3f& b) override { impl_.updateMag(b); }
    void update(float dt, const Vector3f& g, const Vector3f& a, float t) override { impl_.update(dt, g, a, t); }
    FilterSnapshot snapshot() const override { return impl_.snapshot(); }
private:
    CarmaAdapter<4, false, 4, false> impl_;
};

// Regression bars of CSM on the eight scored records (default seeds,
// magnetometer on): the measured worst value plus about 0.5 %, the margin
// carma3-sim uses for A41.
constexpr W3dFailureLimits kCsmRegressionBars{
    .err_limit_percent_z_jonswap   = 3.244f,  // worst 3.22788 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 3.393f,  // worst 3.37567 (pmstokes H0.27)
    .err_limit_yaw_deg             = 0.7319f, // worst 0.728301 (jonswap H4.0)
    .err_limit_roll_deg            = 0.3587f, // worst 0.356937 (pmstokes H8.5)
    .err_limit_pitch_deg           = 0.2247f, // worst 0.223595 (pmstokes H8.5)
    .err_limit_percent_3d_jonswap  = 4.674f,  // worst 4.65064 (jonswap H0.27)
    .err_limit_percent_3d_pmstokes = 5.007f,  // worst 4.9824 (pmstokes H0.27)
    .acc_z_bias_percent            = 4.191f,  // worst 4.17063 (pmstokes H8.5)
    .bias_3d_percent               = 84.15f,  // worst 83.728 (pmstokes H8.5, accel)
    .gyro_bias_3d_percent          = 12.71f,  // worst 12.6465 (pmstokes H8.5, gyro)
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
constexpr W3dSummaryLabels kLabels{.target = "coupling_confidence", .applied = "kappa_z_applied"};

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

int dataset_main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") with_mag = false;
        else if (arg == "--no-noise") add_noise = false;
        else if (arg == "--input" && i + 1 < argc) files.emplace_back(argv[++i]);
        else {
            std::cerr << "Usage: " << argv[0] << " [--nomag] [--no-noise] [--input PATH]...\n"
                      << "       " << argv[0] << " --scenario NAME --seed N --filters A41,CSM,... (Monte Carlo)\n";
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
        const bool default_seeds = !std::getenv("W3D_SEED") && !std::getenv("W3D_IMU_SEED") && !std::getenv("W3D_INIT_SEED");
        const bool gated = variant == "CSM" && add_noise && with_mag && default_seeds;
        const W3dFailureLimits& bars = gated ? kCsmRegressionBars : kUngatedBars;
        for (const auto& f : files) {
            const std::string family = "CARMA_SHARED_MODES_" + variant;
            if (variant == "CSM") run_file<CsmSimAdapter<0>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "CSM_NOCONF") run_file<CsmSimAdapter<1>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "CSM_OFF") run_file<CsmSimAdapter<2>>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else if (variant == "A41") run_file<A41Adapter>(f, dt, with_mag, seeds, write_timeseries, window, family, bars);
            else throw std::invalid_argument("W3D_CSM_VARIANT must be CSM, CSM_NOCONF, CSM_OFF or A41");
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}

}  // namespace

int main(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!std::strcmp(a, "--scenario") || !std::strcmp(a, "--header") || !std::strcmp(a, "--level"))
            return carma3_bench_main(argc, argv);
    }
    return dataset_main(argc, argv);
}
