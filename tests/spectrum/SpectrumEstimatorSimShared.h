#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
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
    std::mt19937 rng;
    std::normal_distribution<float> dist;
    Vector3f bias;
};

inline uint32_t stable_string_hash(const std::string& s) {
    // FNV-1a 32-bit for deterministic cross-run/per-file seeding.
    uint32_t h = 2166136261u;
    for (unsigned char c : s) {
        h ^= static_cast<uint32_t>(c);
        h *= 16777619u;
    }
    return h;
}

inline uint32_t make_file_seed(const std::string& file_name, unsigned base_seed) {
    uint32_t h = stable_string_hash(file_name);
    // Simple deterministic combine.
    return static_cast<uint32_t>(base_seed) ^ (h + 0x9e3779b9u + (static_cast<uint32_t>(base_seed) << 6) +
                                               (static_cast<uint32_t>(base_seed) >> 2));
}

inline NoiseModel make_noise_model(float sigma, float bias_half_range, unsigned seed) {
    NoiseModel m{std::mt19937(seed),
                 std::normal_distribution<float>(0.0f, sigma),
                 Vector3f::Zero()};
    std::uniform_real_distribution<float> ub(-bias_half_range, bias_half_range);
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
    std::string::const_iterator search_start(s.cbegin());

    while (std::regex_search(search_start, s.cend(), match, num_re)) {
        result.append(match.prefix().first, match.prefix().second);
        std::string num = match.str();
        const auto last_non_zero = num.find_last_not_of('0');
        if (last_non_zero != std::string::npos) {
            num.erase(last_non_zero + 1);
        }
        if (!num.empty() && num.back() == '.') num.pop_back();
        result.append(num);
        search_start = match.suffix().first;
    }

    result.append(search_start, s.cend());
    return result;
}

inline std::string build_output_tail(const std::string& data_file) {
    std::string stem = std::filesystem::path(data_file).filename().string();
    const auto pos_h = stem.find("_H");
    std::string tail = (pos_h != std::string::npos) ? stem.substr(pos_h) : "";
    if (tail.size() > 4 && tail.substr(tail.size() - 4) == ".csv") {
        tail.resize(tail.size() - 4);
    }
    return normalize_numbers(tail);
}

inline double quantize_freq(double f_hz) {
    constexpr double q = 1e-6;
    return std::round(f_hz / q) * q;
}

template <typename FreqContainer>
inline double bin_width_hz(const FreqContainer& freqs, int i) {
    const int n = static_cast<int>(freqs.size());
    if (n <= 1) return 0.0;

    double df = 0.0;
    if (i <= 0) {
        df = freqs[1] - freqs[0];
    } else if (i >= n - 1) {
        df = freqs[n - 1] - freqs[n - 2];
    } else {
        df = 0.5 * (freqs[i + 1] - freqs[i - 1]);
    }
    return std::max(df, 0.0);
}

template <size_t Nfreq>
inline std::vector<double> build_reference_1d_spectrum(
    const std::string& spectrum_file,
    const std::array<double, Nfreq>& freqs) {

    std::map<double, std::vector<WaveSpectrumRecord>> by_freq;
    WaveSpectrumCSVReader reader(spectrum_file);
    reader.for_each_record([&](const WaveSpectrumRecord& rec) {
        by_freq[quantize_freq(rec.f_Hz)].push_back(rec);
    });

    std::vector<double> f_ref;
    std::vector<double> s_ref;
    f_ref.reserve(by_freq.size());
    s_ref.reserve(by_freq.size());

    for (auto& [f_hz, rows] : by_freq) {
        if (rows.empty()) continue;

        std::sort(rows.begin(), rows.end(),
                  [](const WaveSpectrumRecord& a, const WaveSpectrumRecord& b) {
                      return a.theta_deg < b.theta_deg;
                  });

        if (rows.size() == 1U) {
            // Single directional sample: treat as full-circle lump.
            const double integ = rows.front().E * (2.0 * M_PI);
            f_ref.push_back(f_hz);
            s_ref.push_back(std::max(0.0, integ));
            continue;
        }

        double integ = 0.0;
        for (size_t i = 0; i < rows.size(); ++i) {
            const auto& r0 = rows[i];
            const auto& r1 = rows[(i + 1) % rows.size()];

            double theta1 = r1.theta_deg;
            if (i + 1 == rows.size()) theta1 += 360.0;

            const double dtheta_rad =
                std::max(0.0, theta1 - static_cast<double>(r0.theta_deg)) * (M_PI / 180.0);

            // Trapezoidal integration over direction.
            integ += 0.5 * (static_cast<double>(r0.E) + static_cast<double>(r1.E)) * dtheta_rad;
        }

        f_ref.push_back(f_hz);
        s_ref.push_back(std::max(0.0, integ));
    }

    std::vector<double> out(freqs.size(), 0.0);
    if (f_ref.empty()) return out;

    for (size_t i = 0; i < freqs.size(); ++i) {
        const double f = freqs[i];

        // Outside reference support: no flat extrapolation, just zero.
        if (f < f_ref.front() || f > f_ref.back()) {
            out[i] = 0.0;
            continue;
        }

        auto it = std::lower_bound(f_ref.begin(), f_ref.end(), f);
        if (it == f_ref.begin()) {
            out[i] = s_ref.front();
            continue;
        }
        if (it == f_ref.end()) {
            out[i] = s_ref.back();
            continue;
        }

        const size_t hi = static_cast<size_t>(std::distance(f_ref.begin(), it));
        const size_t lo = hi - 1;

        const double flo = f_ref[lo];
        const double fhi = f_ref[hi];
        const double denom = fhi - flo;

        if (denom <= 0.0) {
            out[i] = s_ref[lo];
            continue;
        }

        const double t = (f - flo) / denom;
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
    const auto pos = spectrum_file.find("wave_data_");
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

    constexpr float NOISE_STDDEV = 0.03f;
    constexpr float BIAS_HALF_RANGE = 0.02f;
    const uint32_t file_seed = make_file_seed(data_file, noise_seed);
    NoiseModel accel_noise = make_noise_model(NOISE_STDDEV, BIAS_HALF_RANGE, file_seed);

    const double fs_hz = (dt > 0.0f) ? (1.0 / static_cast<double>(dt)) : 200.0;
    TEstimator estimator(fs_hz, 30, true);

    auto freqs = estimator.getFrequencies();
    auto s_ref_interp = build_reference_1d_spectrum(spectrum_file, freqs);

    int block_count = 0;
    bool has_last = false;

    constexpr int KEEP_LAST_SPECTRA = 4;
    std::deque<Eigen::Matrix<double, Nfreq, 1>> tail_spectra;

    WaveDataCSVReader reader(data_file);
    reader.for_each_record([&](const Wave_Data_Sample& rec) {
        const Vector3f acc_b(rec.imu.acc_bx, rec.imu.acc_by, rec.imu.acc_bz);
        const Vector3f acc_noisy = apply_noise(acc_b, accel_noise);
        const Vector3f acc_f = zu_to_ned(acc_noisy);
        const double a_vert = static_cast<double>(acc_f.z());

        if (estimator.processSample(a_vert)) {
            ++block_count;
            const auto spec = estimator.getDisplacementSpectrum();
            if (static_cast<int>(tail_spectra.size()) == KEEP_LAST_SPECTRA) {
                tail_spectra.pop_front();
            }
            tail_spectra.push_back(spec);
            has_last = true;
        }
    });

    if (!has_last || tail_spectra.empty()) {
        std::cout << "Finished " << data_file << " with " << block_count
                  << " spectra (no final snapshot)\n\n";
        return true;
    }

    Eigen::Matrix<double, Nfreq, 1> s_est_avg = Eigen::Matrix<double, Nfreq, 1>::Zero();
    for (const auto& s : tail_spectra) s_est_avg += s;
    s_est_avg /= static_cast<double>(tail_spectra.size());

    const std::string wave_name = EnumTraits<WaveType>::to_string(type);
    const std::string tail = build_output_tail(data_file);

    std::ostringstream nb;
    nb << "_N" << std::fixed << std::setprecision(3) << NOISE_STDDEV
       << "_B" << std::fixed << std::setprecision(3) << BIAS_HALF_RANGE
       << "_S" << noise_seed;

    const std::string outname = output_prefix + wave_name + tail + nb.str() + ".csv";

    std::ofstream ofs(outname);
    if (!ofs.is_open()) {
        std::cerr << "Failed to open output file " << outname << "\n";
        return false;
    }

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
    double sum_sq_amp_ref = 0.0;
    double sum_weight = 0.0;

    for (int i = 0; i < Nfreq; ++i) {
        const double f_est = freqs[i];
        const double s_est = std::max(0.0, static_cast<double>(s_est_avg[i]));
        const double s_ref = std::max(0.0, s_ref_interp[i]);
        const double delta_f = bin_width_hz(freqs, i);

        const double ratio = (s_ref > 0.0) ? (s_est / s_ref) : 0.0;

        // Equivalent single-frequency band amplitude from one-sided PSD bin:
        // A ~= sqrt(2 S df)
        const double a_eta_est = std::sqrt(std::max(0.0, 2.0 * s_est * delta_f));
        const double a_eta_ref = std::sqrt(std::max(0.0, 2.0 * s_ref * delta_f));
        const double a_eta_err = a_eta_est - a_eta_ref;

        // Keep these as diagnostic "energy-weighted by frequency" columns if useful.
        const double e_eta_est = f_est * s_est;
        const double e_eta_ref = f_est * s_ref;

        // Variance / m0 accumulation for one-sided PSD: integral S(f) df
        cum_est += s_est * delta_f;
        cum_ref += s_ref * delta_f;

        // Weight RMS by bin width so denser grids do not distort the metric.
        sum_sq_amp_error += a_eta_err * a_eta_err * delta_f;
        sum_sq_amp_ref += a_eta_ref * a_eta_ref * delta_f;
        sum_weight += delta_f;

        ofs << f_est << "," << s_est << "," << s_ref << "," << ratio << ","
            << a_eta_est << "," << a_eta_ref << ","
            << e_eta_est << "," << e_eta_ref << ","
            << cum_est << "," << cum_ref << ","
            << hs_est << "," << fp_est << "," << tp_est << "\n";
    }

    const double amp_err_rms_m =
        (sum_weight > 0.0) ? std::sqrt(sum_sq_amp_error / sum_weight) : 0.0;
    const double amp_ref_rms_m =
        (sum_weight > 0.0) ? std::sqrt(sum_sq_amp_ref / sum_weight) : 0.0;

    const double amp_err_rms_pct_value =
        (amp_ref_rms_m > 0.0) ? (100.0 * amp_err_rms_m / amp_ref_rms_m) : 0.0;

    constexpr double AMP_ERR_RMS_LIMIT_PCT_VALUE = 2000.0;
    const bool quality_ok = amp_err_rms_pct_value <= AMP_ERR_RMS_LIMIT_PCT_VALUE;

    std::cout << "Spectrum amplitude error RMS=" << amp_err_rms_m
              << " m (" << amp_err_rms_pct_value << "% of reference RMS, gate "
              << AMP_ERR_RMS_LIMIT_PCT_VALUE << "% of reference RMS) -> "
              << (quality_ok ? "PASS" : "FAIL") << "\n";

    std::cout << "Finished " << data_file << " with " << block_count
              << " spectra -> " << outname << "\n\n";

    return quality_ok;
}

inline std::vector<std::string> discover_wave_data_files() {
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        if (fname.find("wave_data_") == std::string::npos) continue;
        auto kind = WaveFileNaming::parse_kind_only(fname);
        if (kind && *kind == FileKind::Data) files.push_back(fname);
    }
    std::sort(files.begin(), files.end());
    return files;
}

}  // namespace SpectrumEstimatorSimShared
