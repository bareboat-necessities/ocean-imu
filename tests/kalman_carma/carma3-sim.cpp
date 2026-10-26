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
    S_x, S_y), A41S (OU-III r_S carried over, negative control) or H41
    (horizontal shape fitted to the carma3_bench sway RAO).
*/

#define EIGEN_NON_ARDUINO

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/W3dSimCommon.h"
#include "carma3/FullAdapters.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {

bool add_noise = true;

template<int NHC, bool SH, int NZC, bool SZ>
class CarmaAdapter final : public IW3dFusionAdapter {
public:
    static constexpr float SIGMA_A_RESCALE = 0.71f;
    static constexpr float SIGMA_G_RESCALE = 0.01f;
    static constexpr float SIGMA_M_RESCALE = 8.0f;

    CarmaAdapter(bool with_mag, const Vector3f& sigma_a_init, const Vector3f& sigma_g,
                 const Vector3f& sigma_m, const ocean_imu::kalman::carma::CarmaShape* hshape) {
        cfg_.with_mag = with_mag;
        cfg_.sigma_a = sigma_a_init * SIGMA_A_RESCALE;
        cfg_.sigma_g = sigma_g * SIGMA_G_RESCALE;
        cfg_.sigma_m = sigma_m * SIGMA_M_RESCALE;
        cfg_.mag_delay_sec = MAG_DELAY_SEC;
        fusion_.begin(cfg_);
        auto& filter = fusion_.raw();
        wave_direction::VesselRaoEqualizer::Config direction_rao;
        direction_rao.enabled = true;
        direction_rao.horizontal_x_tau_s = 1.0;
        direction_rao.horizontal_y_tau_s = 0.7;
        filter.setDirectionRao(direction_rao);
        wave_direction::VesselRaoNoiseWeighting low_wave_noise;
        low_wave_noise.response = direction_rao;
        low_wave_noise.max_std_scale = 4.0f;
        filter.setLowWaveNoiseWeighting(low_wave_noise);
        filter.enableTuner(true);
        filter.enableClamp(true);
        if (hshape) filter.setHorizontalShape(*hshape);
    }

    void updateMag(const Vector3f& mag_body_ned) override { fusion_.updateMag(mag_body_ned); }

    void update(float dt, const Vector3f& gyr, const Vector3f& acc, float temperature_c) override {
        fusion_.update(dt, gyr, acc, temperature_c);
    }

    FilterSnapshot snapshot() const override {
        const auto& filter = fusion_.raw();
        const auto& d = filter.dir();
        FilterSnapshot s;
        s.disp_est_zu = ned_to_zu(filter.mekf().get_position());
        s.vel_est_zu = ned_to_zu(filter.mekf().get_velocity());
        s.acc_est_zu = ned_to_zu(filter.mekf().get_world_accel());
        const Quaternionf q_bw_ned = filter.mekf().quaternion_boat().normalized();
        float roll_deg = 0.0f, pitch_deg = 0.0f, yaw_deg = 0.0f;
        quat_to_euler_nautical(q_bw_ned, roll_deg, pitch_deg, yaw_deg);
        s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, wrapDeg(yaw_deg));
        s.acc_bias_est_ned = filter.mekf().get_acc_bias();
        s.gyro_bias_est_ned = filter.mekf().gyroscope_bias();
        s.mag_bias_est_ned_uT = fusion_.magHardIronBodyUT();
        s.tau_target = filter.getTauTarget();
        s.sigma_target = filter.getSigmaTarget();
        s.tuning_target = filter.getRSTarget();
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
    using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, NHC, SH, NZC, SZ>;
    mutable Fusion fusion_;
    typename Fusion::Config cfg_{};
};

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

template<class A>
void run_file(const std::string& filename, float dt, bool with_mag, const W3dRandomSeeds& seeds,
              bool write_timeseries, float window, const std::string& family) {
    auto result = process_wave_file_for_tracker<A>(filename, dt, with_mag, add_noise, 25.0f,
                                                   "_fusion_" + family, "_fusion_" + family + "_nomag",
                                                   seeds, write_timeseries);
    if (result && window > 0.0f) print_validation_metrics(*result, dt, window, family.c_str());
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
        for (const auto& f : files) {
            const std::string family = "CARMA3_" + variant;
            if (variant == "A41") run_file<Adapter_A41>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else if (variant == "A52") run_file<Adapter_A52>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else if (variant == "M45") run_file<Adapter_M45>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else if (variant == "M54") run_file<Adapter_M54>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else if (variant == "A52S") run_file<Adapter_A52S>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else if (variant == "A41S") run_file<Adapter_A41S>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else if (variant == "H41") run_file<Adapter_H41>(f, dt, with_mag, seeds, write_timeseries, window, family);
            else throw std::invalid_argument("W3D_CARMA3_VARIANT must be A41, A52, M45, M54, A52S, A41S or H41");
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
