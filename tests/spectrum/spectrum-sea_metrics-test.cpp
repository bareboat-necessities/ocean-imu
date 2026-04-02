#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <optional>
#include <numbers>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "spectrum/SeaMetricsFromSpectrum.h"
#include "util/WaveFilesSupport.h"

namespace {

struct CaseSpec {
    std::string name;
    double hs_target_m;
    double tp_target_s;
    double sigma_hz;
};

struct CaseResult {
    std::string name;
    std::string source_file;
    double hs_target_m = 0.0;
    double hs_est_m = 0.0;
    double hs_error_pct = 0.0;
    double tp_target_s = 0.0;
    double tp_est_s = 0.0;
    double tp_error_pct = 0.0;
    double rbw = 0.0;
    double peakedness_ochi_q = 0.0;
    std::string spectrum_type;
    bool pass = false;
};

struct ModeSpec {
    std::string mode_name;
    std::string file_prefix;
    std::string csv_out;
    std::string full_metrics_out;
};

std::vector<double> build_frequency_grid() {
    constexpr int N = 128;
    const double f_min = 0.03;
    const double f_max = 1.0;

    std::vector<double> f;
    f.reserve(N);
    const double step = (f_max - f_min) / static_cast<double>(N - 1);
    for (int i = 0; i < N; ++i) {
        f.push_back(f_min + static_cast<double>(i) * step);
    }
    return f;
}

double trapz(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() < 2 || x.size() != y.size()) return 0.0;

    double sum = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        const double dx = x[i + 1] - x[i];
        sum += 0.5 * dx * (y[i] + y[i + 1]);
    }
    return sum;
}

std::string to_spectrum_type(SeaMetricsFromSpectrum::SpectrumType t) {
    switch (t) {
        case SeaMetricsFromSpectrum::SpectrumType::Swell: return "Swell";
        case SeaMetricsFromSpectrum::SpectrumType::Mixed: return "Mixed";
        case SeaMetricsFromSpectrum::SpectrumType::WindSea: return "WindSea";
        default: return "Unknown";
    }
}

void print_report_summary(const SeaMetricsFromSpectrum::Report& r, std::ostream& os) {
    os << "Last spectrum full metrics summary\n";
    os << std::fixed << std::setprecision(6);
    os << "  valid=" << (r.valid ? 1 : 0) << "\n";
    os << "  lowfreq_guard_hz=" << r.lowfreq_guard_hz << "\n";
    os << "  observation_duration_s=" << r.observation_duration_s << "\n";
    os << "  m_neg1=" << r.m_neg1 << "\n";
    os << "  m0=" << r.m0 << "\n";
    os << "  m1=" << r.m1 << "\n";
    os << "  m2=" << r.m2 << "\n";
    os << "  m3=" << r.m3 << "\n";
    os << "  m4=" << r.m4 << "\n";
    os << "  mu2=" << r.mu2 << "\n";
    os << "  mu3=" << r.mu3 << "\n";
    os << "  mu4=" << r.mu4 << "\n";
    os << "  fp_hz=" << r.fp_hz << "\n";
    os << "  fp_rad=" << r.fp_rad << "\n";
    os << "  tp_s=" << r.tp_s << "\n";
    os << "  mean_frequency_hz=" << r.mean_frequency_hz << "\n";
    os << "  mean_frequency_rad=" << r.mean_frequency_rad << "\n";
    os << "  t1_s=" << r.t1_s << "\n";
    os << "  tm02_s=" << r.tm02_s << "\n";
    os << "  te_s=" << r.te_s << "\n";
    os << "  tm0m1_s=" << r.tm0m1_s << "\n";
    os << "  tm1m1_s=" << r.tm1m1_s << "\n";
    os << "  mean_group_period_s=" << r.mean_group_period_s << "\n";
    os << "  rms_displacement_m=" << r.rms_displacement_m << "\n";
    os << "  hs_m=" << r.hs_m << "\n";
    os << "  h1_10_crest_m=" << r.h1_10_crest_m << "\n";
    os << "  h1_100_crest_m=" << r.h1_100_crest_m << "\n";
    os << "  most_probable_max_crest_m=" << r.most_probable_max_crest_m << "\n";
    os << "  expected_max_crest_m=" << r.expected_max_crest_m << "\n";
    os << "  upcrossing_rate_hz=" << r.upcrossing_rate_hz << "\n";
    os << "  downcrossing_rate_hz=" << r.downcrossing_rate_hz << "\n";
    os << "  wave_count.expected=" << r.wave_count.expected << "\n";
    os << "  wave_count.ci_lower=" << r.wave_count.ci_lower << "\n";
    os << "  wave_count.ci_upper=" << r.wave_count.ci_upper << "\n";
    os << "  wave_count.confidence=" << r.wave_count.confidence << "\n";
    os << "  rbw=" << r.rbw << "\n";
    os << "  regularity_spec=" << r.regularity_spec << "\n";
    os << "  regularity_phase=" << r.regularity_phase << "\n";
    os << "  rbw_phase_increment=" << r.rbw_phase_increment << "\n";
    os << "  narrowness_nu=" << r.narrowness_nu << "\n";
    os << "  bandwidth_clh=" << r.bandwidth_clh << "\n";
    os << "  bandwidth_goda=" << r.bandwidth_goda << "\n";
    os << "  bandwidth_kuik=" << r.bandwidth_kuik << "\n";
    os << "  bandwidth_epsilon=" << r.bandwidth_epsilon << "\n";
    os << "  longuet_higgins_width=" << r.longuet_higgins_width << "\n";
    os << "  spectral_bandwidth_hz=" << r.spectral_bandwidth_hz << "\n";
    os << "  spectral_narrowness_ratio=" << r.spectral_narrowness_ratio << "\n";
    os << "  spectral_skewness=" << r.spectral_skewness << "\n";
    os << "  spectral_kurtosis=" << r.spectral_kurtosis << "\n";
    os << "  spectral_excess_kurtosis=" << r.spectral_excess_kurtosis << "\n";
    os << "  peakedness_ochi_q=" << r.peakedness_ochi_q << "\n";
    os << "  benassai_parameter=" << r.benassai_parameter << "\n";
    os << "  peak_enhancement_gamma=" << r.peak_enhancement_gamma << "\n";
    os << "  peakedness_goda=" << r.peakedness_goda << "\n";
    os << "  spectrum_type=" << to_spectrum_type(r.spectrum_type) << "\n";
    os << "  deep_water_wavelength_tp_m=" << r.deep_water_wavelength_tp_m << "\n";
    os << "  wave_steepness=" << r.wave_steepness << "\n";
    os << "  nonlinearity_parameter=" << r.nonlinearity_parameter << "\n";
    os << "  ursell_number=" << r.ursell_number << "\n";
    os << "  wave_age=" << r.wave_age << "\n";
    os << "  crest_exceed_prob_hs=" << r.crest_exceed_prob_hs << "\n";
    os << "  crest_exceed_prob_tayfun_hs=" << r.crest_exceed_prob_tayfun_hs << "\n";
    os << "  pot_mean_excess_m=" << r.pot_mean_excess_m << "\n";
    os << "  weibull_return_height_m=" << r.weibull_return_height_m << "\n";
    os << "  groupiness_factor=" << r.groupiness_factor << "\n";
    os << "  bfi=" << r.bfi << "\n";
    os << "  mean_group_length=" << r.mean_group_length << "\n";
    os << "  group_height_factor=" << r.group_height_factor << "\n";
    os << "  expected_run_length_above_hs=" << r.expected_run_length_above_hs << "\n";
    os << "  wave_energy_density_j_m2=" << r.wave_energy_density_j_m2 << "\n";
    os << "  wave_power_w_m=" << r.wave_power_w_m << "\n";
    os << "  energy_flux_period_s=" << r.energy_flux_period_s << "\n";
    os << "  breaking_probability=" << r.breaking_probability << "\n";
    os << "  depth_limited_breaking_height_m=" << r.depth_limited_breaking_height_m << "\n";
    os << "  bottom_orbital_velocity_mps=" << r.bottom_orbital_velocity_mps << "\n";
    os << "  radiation_stress_xx_n_m=" << r.radiation_stress_xx_n_m << "\n";
    os << "  radiation_stress_yy_n_m=" << r.radiation_stress_yy_n_m << "\n";
    os << "  radiation_stress_xy_n_m=" << r.radiation_stress_xy_n_m << "\n";
    os << "  msdv=" << r.msdv << "\n";
    os << "  seasickness_incidence_pct=" << r.seasickness_incidence_pct << "\n";
    os << "  motion_comfort_level_0_100=" << r.motion_comfort_level_0_100 << "\n";
    os << "  vertical_motion_intensity=" << r.vertical_motion_intensity << "\n";
    os << "  motion_character=" << r.motion_character << "\n";
    os << "  time_to_onset_min=" << r.time_to_onset_min << "\n";
    os << "  snr_db=" << r.snr_db << "\n";
    os << "  temporal_stability_hs=" << r.temporal_stability_hs << "\n";
    os << "  temporal_stability_tp=" << r.temporal_stability_tp << "\n";
    os << "  data_quality_0_1=" << r.data_quality_0_1 << "\n";
    os << "  wec_capture_width_ratio=" << r.wec_capture_width_ratio << "\n";
    os << "  sea_swell_partition=" << r.sea_swell_partition << "\n";
    os << "  wave_age_class=" << r.wave_age_class << "\n";
}

void write_report_summary_name_value(const SeaMetricsFromSpectrum::Report& r, std::ostream& os) {
    os << std::fixed << std::setprecision(6);
    os << "valid=" << (r.valid ? 1 : 0) << "\n";
    os << "lowfreq_guard_hz=" << r.lowfreq_guard_hz << "\n";
    os << "observation_duration_s=" << r.observation_duration_s << "\n";
    os << "m_neg1=" << r.m_neg1 << "\n";
    os << "m0=" << r.m0 << "\n";
    os << "m1=" << r.m1 << "\n";
    os << "m2=" << r.m2 << "\n";
    os << "m3=" << r.m3 << "\n";
    os << "m4=" << r.m4 << "\n";
    os << "mu2=" << r.mu2 << "\n";
    os << "mu3=" << r.mu3 << "\n";
    os << "mu4=" << r.mu4 << "\n";
    os << "fp_hz=" << r.fp_hz << "\n";
    os << "fp_rad=" << r.fp_rad << "\n";
    os << "tp_s=" << r.tp_s << "\n";
    os << "mean_frequency_hz=" << r.mean_frequency_hz << "\n";
    os << "mean_frequency_rad=" << r.mean_frequency_rad << "\n";
    os << "t1_s=" << r.t1_s << "\n";
    os << "tm02_s=" << r.tm02_s << "\n";
    os << "te_s=" << r.te_s << "\n";
    os << "tm0m1_s=" << r.tm0m1_s << "\n";
    os << "tm1m1_s=" << r.tm1m1_s << "\n";
    os << "mean_group_period_s=" << r.mean_group_period_s << "\n";
    os << "rms_displacement_m=" << r.rms_displacement_m << "\n";
    os << "hs_m=" << r.hs_m << "\n";
    os << "h1_10_crest_m=" << r.h1_10_crest_m << "\n";
    os << "h1_100_crest_m=" << r.h1_100_crest_m << "\n";
    os << "most_probable_max_crest_m=" << r.most_probable_max_crest_m << "\n";
    os << "expected_max_crest_m=" << r.expected_max_crest_m << "\n";
    os << "upcrossing_rate_hz=" << r.upcrossing_rate_hz << "\n";
    os << "downcrossing_rate_hz=" << r.downcrossing_rate_hz << "\n";
    os << "wave_count.expected=" << r.wave_count.expected << "\n";
    os << "wave_count.ci_lower=" << r.wave_count.ci_lower << "\n";
    os << "wave_count.ci_upper=" << r.wave_count.ci_upper << "\n";
    os << "wave_count.confidence=" << r.wave_count.confidence << "\n";
    os << "rbw=" << r.rbw << "\n";
    os << "regularity_spec=" << r.regularity_spec << "\n";
    os << "regularity_phase=" << r.regularity_phase << "\n";
    os << "rbw_phase_increment=" << r.rbw_phase_increment << "\n";
    os << "narrowness_nu=" << r.narrowness_nu << "\n";
    os << "bandwidth_clh=" << r.bandwidth_clh << "\n";
    os << "bandwidth_goda=" << r.bandwidth_goda << "\n";
    os << "bandwidth_kuik=" << r.bandwidth_kuik << "\n";
    os << "bandwidth_epsilon=" << r.bandwidth_epsilon << "\n";
    os << "longuet_higgins_width=" << r.longuet_higgins_width << "\n";
    os << "spectral_bandwidth_hz=" << r.spectral_bandwidth_hz << "\n";
    os << "spectral_narrowness_ratio=" << r.spectral_narrowness_ratio << "\n";
    os << "spectral_skewness=" << r.spectral_skewness << "\n";
    os << "spectral_kurtosis=" << r.spectral_kurtosis << "\n";
    os << "spectral_excess_kurtosis=" << r.spectral_excess_kurtosis << "\n";
    os << "peakedness_ochi_q=" << r.peakedness_ochi_q << "\n";
    os << "benassai_parameter=" << r.benassai_parameter << "\n";
    os << "peak_enhancement_gamma=" << r.peak_enhancement_gamma << "\n";
    os << "peakedness_goda=" << r.peakedness_goda << "\n";
    os << "spectrum_type=" << to_spectrum_type(r.spectrum_type) << "\n";
    os << "deep_water_wavelength_tp_m=" << r.deep_water_wavelength_tp_m << "\n";
    os << "wave_steepness=" << r.wave_steepness << "\n";
    os << "nonlinearity_parameter=" << r.nonlinearity_parameter << "\n";
    os << "ursell_number=" << r.ursell_number << "\n";
    os << "wave_age=" << r.wave_age << "\n";
    os << "crest_exceed_prob_hs=" << r.crest_exceed_prob_hs << "\n";
    os << "crest_exceed_prob_tayfun_hs=" << r.crest_exceed_prob_tayfun_hs << "\n";
    os << "pot_mean_excess_m=" << r.pot_mean_excess_m << "\n";
    os << "weibull_return_height_m=" << r.weibull_return_height_m << "\n";
    os << "groupiness_factor=" << r.groupiness_factor << "\n";
    os << "bfi=" << r.bfi << "\n";
    os << "mean_group_length=" << r.mean_group_length << "\n";
    os << "group_height_factor=" << r.group_height_factor << "\n";
    os << "expected_run_length_above_hs=" << r.expected_run_length_above_hs << "\n";
    os << "wave_energy_density_j_m2=" << r.wave_energy_density_j_m2 << "\n";
    os << "wave_power_w_m=" << r.wave_power_w_m << "\n";
    os << "energy_flux_period_s=" << r.energy_flux_period_s << "\n";
    os << "breaking_probability=" << r.breaking_probability << "\n";
    os << "depth_limited_breaking_height_m=" << r.depth_limited_breaking_height_m << "\n";
    os << "bottom_orbital_velocity_mps=" << r.bottom_orbital_velocity_mps << "\n";
    os << "radiation_stress_xx_n_m=" << r.radiation_stress_xx_n_m << "\n";
    os << "radiation_stress_yy_n_m=" << r.radiation_stress_yy_n_m << "\n";
    os << "radiation_stress_xy_n_m=" << r.radiation_stress_xy_n_m << "\n";
    os << "msdv=" << r.msdv << "\n";
    os << "seasickness_incidence_pct=" << r.seasickness_incidence_pct << "\n";
    os << "motion_comfort_level_0_100=" << r.motion_comfort_level_0_100 << "\n";
    os << "vertical_motion_intensity=" << r.vertical_motion_intensity << "\n";
    os << "motion_character=" << r.motion_character << "\n";
    os << "time_to_onset_min=" << r.time_to_onset_min << "\n";
    os << "snr_db=" << r.snr_db << "\n";
    os << "temporal_stability_hs=" << r.temporal_stability_hs << "\n";
    os << "temporal_stability_tp=" << r.temporal_stability_tp << "\n";
    os << "data_quality_0_1=" << r.data_quality_0_1 << "\n";
    os << "wec_capture_width_ratio=" << r.wec_capture_width_ratio << "\n";
    os << "sea_swell_partition=" << r.sea_swell_partition << "\n";
    os << "wave_age_class=" << r.wave_age_class << "\n";
}

CaseResult run_case(const CaseSpec& spec, const std::vector<double>& freqs_hz) {
    const double fp_target_hz = 1.0 / spec.tp_target_s;
    std::vector<double> s_raw(freqs_hz.size(), 0.0);

    for (std::size_t i = 0; i < freqs_hz.size(); ++i) {
        const double x = (freqs_hz[i] - fp_target_hz) / spec.sigma_hz;
        s_raw[i] = std::exp(-0.5 * x * x);
    }

    const double m0_raw = trapz(freqs_hz, s_raw);
    const double m0_target = std::pow(spec.hs_target_m / 4.0, 2.0);
    const double scale = (m0_raw > 0.0) ? (m0_target / m0_raw) : 0.0;

    std::vector<double> s_eta(freqs_hz.size(), 0.0);
    for (std::size_t i = 0; i < s_eta.size(); ++i) {
        s_eta[i] = scale * s_raw[i];
    }

    SeaMetricsFromSpectrum metrics;
    metrics.updateFromSpectrum(freqs_hz, s_eta, 0.04, 1.0);
    const auto& r = metrics.report();

    CaseResult out;
    out.name = spec.name;
    out.hs_target_m = spec.hs_target_m;
    out.tp_target_s = spec.tp_target_s;
    out.hs_est_m = r.hs_m;
    out.tp_est_s = r.tp_s;
    out.hs_error_pct = (spec.hs_target_m > 0.0) ? std::abs(100.0 * (r.hs_m - spec.hs_target_m) / spec.hs_target_m) : 0.0;
    out.tp_error_pct = (spec.tp_target_s > 0.0) ? std::abs(100.0 * (r.tp_s - spec.tp_target_s) / spec.tp_target_s) : 0.0;
    out.rbw = r.rbw;
    out.peakedness_ochi_q = r.peakedness_ochi_q;
    out.spectrum_type = to_spectrum_type(r.spectrum_type);

    constexpr double HS_GATE_PCT = 5.0;
    constexpr double TP_GATE_PCT = 8.0;
    out.pass = r.valid && out.hs_error_pct <= HS_GATE_PCT && out.tp_error_pct <= TP_GATE_PCT;

    return out;
}

std::optional<std::tuple<std::vector<double>, std::vector<double>>> read_spectrum_csv(
    const std::string& csv_path) {
    std::ifstream ifs(csv_path);
    if (!ifs.is_open()) return std::nullopt;

    std::string header;
    if (!std::getline(ifs, header)) return std::nullopt;

    std::vector<double> freqs_hz;
    std::vector<double> s_eta;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token_f;
        std::string token_s;
        if (!std::getline(ss, token_f, ',')) continue;
        if (!std::getline(ss, token_s, ',')) continue;
        try {
            freqs_hz.push_back(std::stod(token_f));
            s_eta.push_back(std::max(0.0, std::stod(token_s)));
        } catch (...) {
            continue;
        }
    }

    if (freqs_hz.size() < 2 || freqs_hz.size() != s_eta.size()) return std::nullopt;
    return std::make_tuple(freqs_hz, s_eta);
}

std::optional<WaveFileNaming::ParsedName> parse_estimator_filename(
    const std::string& source_file, const std::string& mode_prefix) {
    if (!source_file.starts_with(mode_prefix)) return std::nullopt;
    const std::string tail = source_file.substr(mode_prefix.size());
    return WaveFileNaming::parse("wave_data_" + tail);
}

std::vector<std::string> discover_estimator_files(const std::string& prefix) {
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        if (!fname.ends_with(".csv")) continue;
        if (!fname.starts_with(prefix)) continue;
        files.push_back(fname);
    }
    std::sort(files.begin(), files.end());
    return files;
}

CaseResult run_case_from_estimator_csv(const std::string& source_file, const std::string& mode_prefix) {
    auto parsed_meta = parse_estimator_filename(source_file, mode_prefix);
    if (!parsed_meta) {
        throw std::runtime_error("Cannot parse estimator filename metadata: " + source_file);
    }

    auto parsed_csv = read_spectrum_csv(source_file);
    if (!parsed_csv) {
        throw std::runtime_error("Cannot parse estimator CSV: " + source_file);
    }
    auto [freqs_hz, s_eta] = *parsed_csv;

    SeaMetricsFromSpectrum metrics;
    metrics.updateFromSpectrum(freqs_hz, s_eta, 0.04, 1.0);
    const auto& r = metrics.report();

    CaseResult out;
    out.name = source_file;
    out.source_file = source_file;
    out.hs_target_m = parsed_meta->height;
    const double tp_target_s = (parsed_meta->length > 0.0)
                                   ? std::sqrt(parsed_meta->length / 9.80665 * 2.0 * std::numbers::pi)
                                   : 0.0;
    out.tp_target_s = tp_target_s;
    out.hs_est_m = r.hs_m;
    out.tp_est_s = r.tp_s;
    out.hs_error_pct = (out.hs_target_m > 0.0) ? std::abs(100.0 * (r.hs_m - out.hs_target_m) / out.hs_target_m) : 0.0;
    out.tp_error_pct = (out.tp_target_s > 0.0) ? std::abs(100.0 * (r.tp_s - out.tp_target_s) / out.tp_target_s) : 0.0;
    out.rbw = r.rbw;
    out.peakedness_ochi_q = r.peakedness_ochi_q;
    out.spectrum_type = to_spectrum_type(r.spectrum_type);

    constexpr double HS_GATE_PCT = 25.0;
    constexpr double TP_GATE_PCT = 35.0;
    out.pass = r.valid && out.hs_error_pct <= HS_GATE_PCT && out.tp_error_pct <= TP_GATE_PCT;
    return out;
}

}  // namespace

int main() {
    const std::vector<ModeSpec> modes = {
        {"classic", "spectrum_estimator_", "sea_metrics_from_spectrum_report_classic.csv", "sea_metrics_from_spectrum_full_metrics_classic.txt"},
        {"wavelets", "spectrum_wavelets_", "sea_metrics_from_spectrum_report_wavelets.csv", "sea_metrics_from_spectrum_full_metrics_wavelets.txt"},
    };
    bool all_ok = true;

    for (const auto& mode : modes) {
        auto files = discover_estimator_files(mode.file_prefix);
        if (files.empty()) {
            std::cout << "No estimator CSV files found for mode=" << mode.mode_name
                      << " prefix=" << mode.file_prefix << " (skipping)\n";
            continue;
        }

        std::ofstream csv(mode.csv_out);
        if (!csv.is_open()) {
            std::cerr << "Failed to open " << mode.csv_out << " for writing\n";
            return EXIT_FAILURE;
        }
        csv << "case_name,source_file,hs_target_m,hs_est_m,hs_error_pct,tp_target_s,tp_est_s,tp_error_pct,rbw,peakedness_ochi_q,spectrum_type,pass\n";
        csv << std::fixed << std::setprecision(6);

        std::cout << "SeaMetricsFromSpectrum report from estimator spectra: " << mode.mode_name << "\n";
        SeaMetricsFromSpectrum::Report last_report;

        for (const auto& source_file : files) {
            CaseResult r = run_case_from_estimator_csv(source_file, mode.file_prefix);
            all_ok = all_ok && r.pass;

            csv << r.name << ',' << r.source_file << ','
                << r.hs_target_m << ',' << r.hs_est_m << ',' << r.hs_error_pct << ','
                << r.tp_target_s << ',' << r.tp_est_s << ',' << r.tp_error_pct << ','
                << r.rbw << ',' << r.peakedness_ochi_q << ','
                << r.spectrum_type << ',' << (r.pass ? 1 : 0) << "\n";

            auto parsed_csv = read_spectrum_csv(source_file);
            if (parsed_csv) {
                auto [freqs_hz, s_eta] = *parsed_csv;
                SeaMetricsFromSpectrum metrics_last;
                metrics_last.updateFromSpectrum(freqs_hz, s_eta, 0.04, 1.0);
                last_report = metrics_last.report();
            }

            std::cout << "  - " << r.name
                      << ": Hs(target=" << r.hs_target_m << " m, est=" << r.hs_est_m << " m, err=" << r.hs_error_pct << "%), "
                      << "Tp(target=" << r.tp_target_s << " s, est=" << r.tp_est_s << " s, err=" << r.tp_error_pct << "%), "
                      << "RBW=" << r.rbw << ", q=" << r.peakedness_ochi_q << ", type=" << r.spectrum_type
                      << " -> " << (r.pass ? "PASS" : "FAIL") << "\n";
        }

        std::ofstream full_metrics_txt(mode.full_metrics_out);
        if (!full_metrics_txt.is_open()) {
            std::cerr << "Failed to open " << mode.full_metrics_out << " for writing\n";
            return EXIT_FAILURE;
        }
        write_report_summary_name_value(last_report, full_metrics_txt);

        std::cout << "CSV written: " << mode.csv_out << "\n";
        std::cout << "Full metrics text written: " << mode.full_metrics_out << "\n";
    }

    return all_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
