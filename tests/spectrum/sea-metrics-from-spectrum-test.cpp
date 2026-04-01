#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "spectrum/SeaMetricsFromSpectrum.h"

namespace {

struct CaseSpec {
    std::string name;
    double hs_target_m;
    double tp_target_s;
    double sigma_hz;
};

struct CaseResult {
    std::string name;
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

}  // namespace

int main() {
    const std::vector<CaseSpec> cases = {
        {"narrow_peak_swell_like", 1.50, 10.0, 0.015},
        {"broad_peak_windsea_like", 1.50, 6.0, 0.050},
    };

    const auto freqs_hz = build_frequency_grid();

    std::ofstream csv("sea_metrics_from_spectrum_report.csv");
    if (!csv.is_open()) {
        std::cerr << "Failed to open sea_metrics_from_spectrum_report.csv for writing\n";
        return EXIT_FAILURE;
    }

    csv << "case_name,hs_target_m,hs_est_m,hs_error_pct,tp_target_s,tp_est_s,tp_error_pct,rbw,peakedness_ochi_q,spectrum_type,pass\n";
    csv << std::fixed << std::setprecision(6);

    std::cout << "SeaMetricsFromSpectrum synthetic-spectrum test report\n";

    bool all_ok = true;
    std::vector<CaseResult> results;
    results.reserve(cases.size());

    for (const auto& c : cases) {
        CaseResult r = run_case(c, freqs_hz);
        all_ok = all_ok && r.pass;
        results.push_back(r);

        csv << r.name << ','
            << r.hs_target_m << ',' << r.hs_est_m << ',' << r.hs_error_pct << ','
            << r.tp_target_s << ',' << r.tp_est_s << ',' << r.tp_error_pct << ','
            << r.rbw << ',' << r.peakedness_ochi_q << ','
            << r.spectrum_type << ',' << (r.pass ? 1 : 0) << "\n";

        std::cout << "  - " << r.name
                  << ": Hs(est=" << r.hs_est_m << " m, err=" << r.hs_error_pct << "%), "
                  << "Tp(est=" << r.tp_est_s << " s, err=" << r.tp_error_pct << "%), "
                  << "RBW=" << r.rbw << ", q=" << r.peakedness_ochi_q << ", type=" << r.spectrum_type
                  << " -> " << (r.pass ? "PASS" : "FAIL") << "\n";
    }

    if (results.size() == 2) {
        const bool narrow_is_narrower = results[0].rbw < results[1].rbw;
        std::cout << "  - cross-check: narrow RBW < broad RBW -> "
                  << (narrow_is_narrower ? "PASS" : "FAIL") << "\n";
        all_ok = all_ok && narrow_is_narrower;
    }

    std::cout << "CSV written: sea_metrics_from_spectrum_report.csv\n";
    return all_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
