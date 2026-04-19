#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO
#include "detrend/AdaptiveWaveDetrender.h"
#include "util/W3dSimCommon.h"

namespace {

struct DriftStats {
    std::string file;
    float rms_integrated = 0.0f;
    float rms_detrended = 0.0f;
};

float wave_freq_from_filename(const std::string& filename) {
    auto parsed = WaveFileNaming::parse_to_params(filename);
    if (!parsed) return 0.12f;

    const auto [kind, type, wp] = *parsed;
    (void)kind;
    (void)type;
    if (wp.period <= 0.0f || !std::isfinite(wp.period)) return 0.12f;
    return 1.0f / wp.period;
}

bool should_process_file(const std::string& name) {
    return name.rfind("wave_data_pmstokes_", 0) == 0 ||
           name.rfind("wave_data_jonswap_", 0) == 0;
}

} // namespace

int main() {
    constexpr float dt_default = 0.005f;
    constexpr double max_test_time_s = 60.0;

    // Keep coefficients aligned with src/util/W3dSimCommon.h process_wave_file_for_tracker.
    const float acc_sigma = 1.51e-3f * g_std;
    const float acc_bias_range = 8e-3f * g_std;
    const float acc_bias_rw = 0.0005f;

    std::ofstream ofs("vertical_accel_drift_timeseries.csv");
    if (!ofs.is_open()) {
        std::cerr << "ERROR: failed to open vertical_accel_drift_timeseries.csv\n";
        return 1;
    }

    ofs << "wave_file,wave_type,time_s,disp_ref_z_m,vel_ref_z_m,acc_ref_z_mps2,"
        << "acc_noisy_z_mps2,disp_int_z_m,disp_int_detrended_z_m,"
        << "err_int_m,err_detrended_m\n";
    ofs << std::fixed << std::setprecision(9);

    std::vector<DriftStats> stats;

    for (const auto& entry : std::filesystem::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        const std::string file = entry.path().filename().string();
        if (!should_process_file(file)) continue;

        WaveDataCSVReader reader(file);
        ImuNoiseModel accel_world_noise = make_imu_noise_model(acc_sigma, acc_bias_range, acc_bias_rw, 1234);

        AdaptiveWaveDetrender detrender;
        const float wave_freq_hz = wave_freq_from_filename(file);

        bool initialized = false;
        double t_prev = 0.0;
        float vel_int = 0.0f;
        float disp_int = 0.0f;

        double sum_err2_int = 0.0;
        double sum_err2_det = 0.0;
        std::size_t n_err = 0;

        const std::string wave_type = (file.find("pmstokes") != std::string::npos) ? "pmstokes" : "jonswap";

        reader.for_each_record([&](const Wave_Data_Sample& rec) {
            if (rec.time > max_test_time_s) return;

            const float acc_truth_z = rec.wave.acc_z;
            const Vector3f noisy_xyz = apply_imu_noise(Vector3f::UnitZ() * acc_truth_z,
                                                        accel_world_noise,
                                                        dt_default);
            const float acc_noisy_z = noisy_xyz.z();

            if (!initialized) {
                initialized = true;
                t_prev = rec.time;
                vel_int = rec.wave.vel_z;
                disp_int = rec.wave.disp_z;
                detrender.reset(disp_int);
            } else {
                const float dt = static_cast<float>(std::max(1e-6, rec.time - t_prev));
                vel_int += acc_noisy_z * dt;
                disp_int += vel_int * dt;
                t_prev = rec.time;
            }

            const auto det_out = detrender.update(disp_int, dt_default, wave_freq_hz, true);
            const float disp_det = det_out.wave_clean;
            const float err_int = disp_int - rec.wave.disp_z;
            const float err_det = disp_det - rec.wave.disp_z;

            sum_err2_int += double(err_int) * double(err_int);
            sum_err2_det += double(err_det) * double(err_det);
            ++n_err;

            ofs << file << ','
                << wave_type << ','
                << rec.time << ','
                << rec.wave.disp_z << ','
                << rec.wave.vel_z << ','
                << acc_truth_z << ','
                << acc_noisy_z << ','
                << disp_int << ','
                << disp_det << ','
                << err_int << ','
                << err_det << '\n';
        });

        if (n_err == 0) continue;
        const float rms_int = std::sqrt(float(sum_err2_int / n_err));
        const float rms_det = std::sqrt(float(sum_err2_det / n_err));
        stats.push_back({file, rms_int, rms_det});

        std::cout << "[vertical-accel-drift] " << file
                  << " rms(integrated)=" << rms_int
                  << " m, rms(detrended)=" << rms_det
                  << " m\n";
    }

    if (stats.empty()) {
        std::cerr << "ERROR: no wave_data_pmstokes_/wave_data_jonswap_ files found\n";
        return 1;
    }

    std::ofstream summary("vertical_accel_drift_summary.csv");
    summary << "wave_file,rms_integrated_m,rms_detrended_m,improvement_ratio\n";
    summary << std::fixed << std::setprecision(9);
    for (const auto& s : stats) {
        const float ratio = (s.rms_detrended > 0.0f) ? (s.rms_integrated / s.rms_detrended) : 0.0f;
        summary << s.file << ',' << s.rms_integrated << ',' << s.rms_detrended << ',' << ratio << '\n';
    }

    std::cout << "Wrote vertical_accel_drift_timeseries.csv and vertical_accel_drift_summary.csv\n";
    return 0;
}
