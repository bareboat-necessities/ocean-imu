#pragma once

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <Eigen/Dense>

#include "ahrs/FrameConversions.h"
#include "util/WaveFilesSupport.h"

namespace SpectrumEstimatorSimShared {

using Eigen::Vector3f;

struct NoiseModel {
    std::default_random_engine rng;
    std::normal_distribution<float> dist;
    Vector3f bias;
};

inline NoiseModel make_noise_model(float sigma, float bias_range, unsigned seed) {
    NoiseModel m{std::default_random_engine(seed),
                 std::normal_distribution<float>(0.0f, sigma),
                 Vector3f::Zero()};
    std::uniform_real_distribution<float> ub(-bias_range, bias_range);
    m.bias = Vector3f(ub(m.rng), ub(m.rng), ub(m.rng));
    return m;
}

inline Vector3f apply_noise(const Vector3f& v, NoiseModel& m) {
    return v + m.bias + Vector3f(m.dist(m.rng), m.dist(m.rng), m.dist(m.rng));
}

inline std::string normalize_numbers(const std::string& s) {
    std::regex num_re(R"(([0-9]+\.[0-9]+))");
    std::smatch match;
    std::string result;
    std::string::const_iterator searchStart(s.cbegin());

    while (std::regex_search(searchStart, s.cend(), match, num_re)) {
        result.append(match.prefix().first, match.prefix().second);
        std::string num = match.str();
        num.erase(num.find_last_not_of('0') + 1, std::string::npos);
        if (!num.empty() && num.back() == '.') num.pop_back();
        result.append(num);
        searchStart = match.suffix().first;
    }
    result.append(searchStart, s.cend());
    return result;
}

inline std::string build_output_tail(const std::string& data_file) {
    std::string stem = std::filesystem::path(data_file).filename().string();
    auto posH = stem.find("_H");
    std::string tail = (posH != std::string::npos) ? stem.substr(posH) : "";
    if (tail.size() > 4 && tail.substr(tail.size() - 4) == ".csv") {
        tail = tail.substr(0, tail.size() - 4);
    }
    return normalize_numbers(tail);
}

template <size_t Nfreq>
inline std::vector<double> build_reference_1d_spectrum(
    const std::string& spectrum_file,
    const std::array<double, Nfreq>& freqs) {
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
            integ += r0.E * (dtheta * M_PI / 180.0);
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

template <typename TEstimator, int Nfreq>
inline bool process_wave_file(const std::string& data_file,
                              [[maybe_unused]] float dt,
                              unsigned noise_seed,
                              const std::string& output_prefix) {
    auto parsed = WaveFileNaming::parse_to_params(data_file);
    if (!parsed) return true;

    auto [kind, type, wp] = *parsed;
    if (kind != FileKind::Data) return true;

    std::string spectrum_file = data_file;
    auto pos = spectrum_file.find("wave_data_");
    if (pos != std::string::npos) {
        spectrum_file.replace(pos, std::string("wave_data_").size(), "wave_spectrum_");
    }

    if (!std::filesystem::exists(spectrum_file)) {
        std::cerr << "Skipping " << data_file << " (missing ref spectrum " << spectrum_file << ")\n";
        return true;
    }

    std::cout << "Processing " << data_file
              << " with ref=" << spectrum_file
              << " (" << EnumTraits<WaveType>::to_string(type)
              << "), Hs=" << wp.height << " m, Tp=" << wp.period << " s\n";

    static NoiseModel accel_noise = make_noise_model(0.03f, 0.02f, noise_seed);

    TEstimator estimator(200.0, 30, true);
    auto freqs = estimator.getFrequencies();
    auto s_ref_interp = build_reference_1d_spectrum(spectrum_file, freqs);

    int block_count = 0;
    Eigen::Matrix<double, Nfreq, 1> last_s_est = Eigen::Matrix<double, Nfreq, 1>::Zero();
    bool has_last = false;
    WaveDataCSVReader reader(data_file);
    reader.for_each_record([&](const Wave_Data_Sample& rec) {
        Vector3f acc_b(rec.imu.acc_bx, rec.imu.acc_by, rec.imu.acc_bz);
        Vector3f acc_noisy = apply_noise(acc_b, accel_noise);
        Vector3f acc_f = zu_to_ned(acc_noisy);
        const double a_vert = static_cast<double>(acc_f.z());

        if (estimator.processSample(a_vert)) {
            block_count++;
            last_s_est = estimator.getDisplacementSpectrum();
            has_last = true;
        }
    });

    if (!has_last) {
        std::cout << "Finished " << data_file << " with " << block_count
                  << " spectra (no final snapshot)\n\n";
        return true;
    }

    const std::string waveName = EnumTraits<WaveType>::to_string(type);
    const std::string tail = build_output_tail(data_file);
    constexpr float NOISE_STDDEV = 0.03f;
    constexpr float BIAS_MEAN = 0.02f;

    std::ostringstream nb;
    nb << "_N" << std::fixed << std::setprecision(3) << NOISE_STDDEV
       << "_B" << std::fixed << std::setprecision(3) << BIAS_MEAN;
    const std::string outname = output_prefix + waveName + tail + nb.str() + ".csv";

    std::ofstream ofs(outname);
    const double hs_est = estimator.computeHs();
    const double fp_est = estimator.estimateFp();
    const double tp_est = (fp_est > 0.0) ? (1.0 / fp_est) : 0.0;

    ofs << "freq_hz,S_eta_hz,S_ref_interp,S_ratio,"
           "A_eta_est,A_eta_ref,E_eta_est,E_eta_ref,"
           "CumVar_est,CumVar_ref,"
           "Hs_est,Fp_est,Tp_est\n";

    double cum_est = 0.0;
    double cum_ref = 0.0;
    double sum_sq_amp_error = 0.0;
    for (int i = 0; i < Nfreq; ++i) {
        const double f_est = freqs[i];
        const double S_eta_hz = std::max(0.0, last_s_est[i]);
        const double s_ref = std::max(0.0, s_ref_interp[i]);

        double delta_f = 0.0;
        if (i == 0) delta_f = 0.5 * (freqs[1] - freqs[0]);
        else if (i == Nfreq - 1) delta_f = 0.5 * (freqs[Nfreq - 1] - freqs[Nfreq - 2]);
        else delta_f = 0.5 * (freqs[i + 1] - freqs[i - 1]);
        delta_f = std::max(delta_f, 0.0);

        const double ratio = (s_ref > 0.0) ? (S_eta_hz / s_ref) : 0.0;
        const double A_eta_est = std::sqrt(std::max(0.0, 2.0 * S_eta_hz * delta_f * f_est));
        const double A_eta_ref = std::sqrt(std::max(0.0, 2.0 * s_ref * delta_f * f_est));
        const double A_eta_err = A_eta_est - A_eta_ref;
        const double E_eta_est = f_est * S_eta_hz;
        const double E_eta_ref = f_est * s_ref;

        cum_est += S_eta_hz * 2.0 * delta_f;
        cum_ref += s_ref * 2.0 * delta_f;
        sum_sq_amp_error += A_eta_err * A_eta_err;

        ofs << f_est << "," << S_eta_hz << "," << s_ref << "," << ratio << ","
            << A_eta_est << "," << A_eta_ref << ","
            << E_eta_est << "," << E_eta_ref << ","
            << cum_est << "," << cum_ref << ","
            << hs_est << "," << fp_est << "," << tp_est << "\n";
    }

    const double amp_err_rms_m = std::sqrt(sum_sq_amp_error / static_cast<double>(Nfreq));
    const double hs_ref_m = static_cast<double>(wp.height);
    const double amp_err_rms_pct_hs = (hs_ref_m > 0.0) ? (100.0 * amp_err_rms_m / hs_ref_m) : 0.0;
    constexpr double AMP_ERR_RMS_LIMIT_PCT_HS = 15.0;
    const bool quality_ok = amp_err_rms_pct_hs <= AMP_ERR_RMS_LIMIT_PCT_HS;

    std::cout << "Spectrum amplitude error RMS=" << amp_err_rms_m
              << " m (" << amp_err_rms_pct_hs << "%Hs, gate "
              << AMP_ERR_RMS_LIMIT_PCT_HS << "%Hs) -> "
              << (quality_ok ? "PASS" : "FAIL") << "\n";

    std::cout << "Finished " << data_file << " with " << block_count
              << " spectra -> " << outname << "\n\n";
    return quality_ok;
}

inline std::vector<std::string> discover_wave_data_files() {
    std::vector<std::string> files;
    for (auto& entry : std::filesystem::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        if (fname.find("wave_data_") == std::string::npos) continue;
        auto kind = WaveFileNaming::parse_kind_only(fname);
        if (kind && *kind == FileKind::Data) files.push_back(fname);
    }
    std::sort(files.begin(), files.end());
    return files;
}

}  // namespace
