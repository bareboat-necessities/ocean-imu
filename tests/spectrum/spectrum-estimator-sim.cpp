#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/WaveFilesSupport.h"
#include "ahrs/FrameConversions.h"

const float g_std = 9.80665f;

#include "spectrum/WaveSpectrumEstimator.h"

using Eigen::Vector3f;

struct NoiseModel {
    std::default_random_engine rng;
    std::normal_distribution<float> dist;
    Vector3f bias;
};

NoiseModel make_noise_model(float sigma, float bias_range, unsigned seed) {
    NoiseModel m{std::default_random_engine(seed),
                 std::normal_distribution<float>(0.0f, sigma),
                 Vector3f::Zero()};
    std::uniform_real_distribution<float> ub(-bias_range, bias_range);
    m.bias = Vector3f(ub(m.rng), ub(m.rng), ub(m.rng));
    return m;
}

Vector3f apply_noise(const Vector3f& v, NoiseModel& m) {
    return v + m.bias + Vector3f(m.dist(m.rng), m.dist(m.rng), m.dist(m.rng));
}

std::vector<double> build_reference_1d_spectrum(const std::string& spectrum_file,
                                                const std::array<double, 32>& freqs) {
    std::map<double, std::vector<WaveSpectrumRecord>> by_freq;
    WaveSpectrumCSVReader reader(spectrum_file);
    reader.for_each_record([&](const WaveSpectrumRecord& rec) {
        by_freq[rec.f_Hz].push_back(rec);
    });

    std::vector<double> f_ref;
    std::vector<double> s_ref;
    f_ref.reserve(by_freq.size());
    s_ref.reserve(by_freq.size());

    for (auto& [f, rows] : by_freq) {
        std::sort(rows.begin(), rows.end(), [](const WaveSpectrumRecord& a, const WaveSpectrumRecord& b) {
            return a.theta_deg < b.theta_deg;
        });

        if (rows.empty()) continue;

        double integ = 0.0;
        for (size_t i = 0; i < rows.size(); ++i) {
            const auto& r0 = rows[i];
            const auto& r1 = rows[(i + 1) % rows.size()];
            double dtheta = r1.theta_deg - r0.theta_deg;
            if (i + 1 == rows.size()) dtheta = (rows.front().theta_deg + 360.0) - r0.theta_deg;
            dtheta = std::max(dtheta, 0.0);
            integ += r0.E * dtheta;
        }

        f_ref.push_back(f);
        s_ref.push_back(integ);
    }

    std::vector<double> out(freqs.size(), 0.0);
    if (f_ref.empty()) return out;

    for (size_t i = 0; i < freqs.size(); ++i) {
        double f = freqs[i];
        if (f <= f_ref.front()) {
            out[i] = s_ref.front();
            continue;
        }
        if (f >= f_ref.back()) {
            out[i] = s_ref.back();
            continue;
        }
        auto it = std::lower_bound(f_ref.begin(), f_ref.end(), f);
        size_t hi = std::distance(f_ref.begin(), it);
        size_t lo = hi - 1;
        double t = (f - f_ref[lo]) / (f_ref[hi] - f_ref[lo]);
        out[i] = (1.0 - t) * s_ref[lo] + t * s_ref[hi];
    }

    return out;
}

void process_wave_file(const std::string& data_file, float dt) {
    auto parsed = WaveFileNaming::parse_to_params(data_file);
    if (!parsed) return;

    auto [kind, type, wp] = *parsed;
    if (kind != FileKind::Data) return;

    std::string spectrum_file = data_file;
    auto pos = spectrum_file.find("wave_data_");
    if (pos != std::string::npos) {
        spectrum_file.replace(pos, std::string("wave_data_").size(), "wave_spectrum_");
    }

    if (!std::filesystem::exists(spectrum_file)) {
        std::cerr << "Skipping " << data_file << " (missing ref spectrum " << spectrum_file << ")\n";
        return;
    }

    std::cout << "Processing " << data_file
              << " with ref=" << spectrum_file
              << " (" << EnumTraits<WaveType>::to_string(type)
              << "), Hs=" << wp.height << " m, Tp=" << wp.period << " s\n";

    static NoiseModel accel_noise = make_noise_model(0.03f, 0.02f, 1234u);

    constexpr int Nfreq = 32;
    constexpr int Nblock = 256;
    WaveSpectrumEstimator<Nfreq, Nblock> estimator(200.0, 30, true);

    std::string outname = data_file;
    if (pos != std::string::npos) {
        outname.replace(pos, std::string("wave_data_").size(), "spectrum_estimate_");
    } else {
        outname = "spectrum_estimate_" + outname;
    }

    std::ofstream ofs(outname);
    auto freqs = estimator.getFrequencies();
    auto s_ref_interp = build_reference_1d_spectrum(spectrum_file, freqs);

    ofs << "block,time,Hs,Fp,PM_alpha,PM_fp,PM_cost";
    for (int i = 0; i < Nfreq; ++i) {
        ofs << ",S_eta_est_f" << i << "_Hz=" << freqs[i]
            << ",S_eta_ref_f" << i;
    }
    ofs << "\n";

    int block_count = 0;
    WaveDataCSVReader reader(data_file);
    reader.for_each_record([&](const Wave_Data_Sample& rec) {
        Vector3f acc_b(rec.imu.acc_bx, rec.imu.acc_by, rec.imu.acc_bz);
        Vector3f acc_noisy = apply_noise(acc_b, accel_noise);
        Vector3f acc_f = zu_to_ned(acc_noisy);
        double a_vert = static_cast<double>(acc_f.z());

        if (estimator.processSample(a_vert)) {
            block_count++;
            auto s_est = estimator.getDisplacementSpectrum();
            auto pm = estimator.fitPiersonMoskowitz();

            ofs << block_count << "," << rec.time << ","
                << estimator.computeHs() << "," << estimator.estimateFp() << ","
                << pm.alpha << "," << pm.fp << "," << pm.cost;
            for (int i = 0; i < Nfreq; ++i) {
                ofs << "," << s_est[i] << "," << s_ref_interp[i];
            }
            ofs << "\n";
        }
    });

    std::cout << "Finished " << data_file << " with " << block_count
              << " spectra -> " << outname << "\n\n";
}

int main() {
    constexpr float dt = 1.0f / 200.0f;
    std::cout << "Spectrum estimator simulation (IMU -> estimator vs ref spectrum files)\n";

    std::vector<std::string> files;
    for (auto& entry : std::filesystem::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        if (fname.find("wave_data_") == std::string::npos) continue;
        auto kind = WaveFileNaming::parse_kind_only(fname);
        if (kind && *kind == FileKind::Data) files.push_back(fname);
    }

    std::sort(files.begin(), files.end());
    for (const auto& fname : files) process_wave_file(fname, dt);

    return 0;
}
