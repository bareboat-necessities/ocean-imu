#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO

#include "detrend/AdaptiveWaveDetrender.h"

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kWarmupMinS = 60.0;
constexpr double kRmsGateFractionOfHeight = 0.16;

struct Scenario {
  const char* tag;
  double height_m;
  double period_s;
  const char* wave_samples_csv;
};

struct WaveRow {
  double time_s = 0.0;
  double disp_z_m = 0.0;
};

struct RMSAccumulator {
  double sum_sq = 0.0;
  std::size_t count = 0U;

  void add(double value) {
    sum_sq += value * value;
    ++count;
  }

  double rms() const {
    return (count > 0U) ? std::sqrt(sum_sq / static_cast<double>(count)) : 0.0;
  }
};

const std::vector<Scenario> kScenarios = {
    {"low", 0.27, 3.0, "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv"},
    {"medium", 1.5, 5.7, "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv"},
    {"large", 4.0, 8.5, "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv"},
    {"extreme", 8.5, 11.4, "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv"},
};

std::string trim(const std::string& value) {
  std::size_t start = 0;
  while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start])) != 0) {
    ++start;
  }

  std::size_t end = value.size();
  while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1])) != 0) {
    --end;
  }

  return value.substr(start, end - start);
}

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, ',')) {
    fields.push_back(trim(field));
  }
  return fields;
}

double parse_double_field(const std::string& value, const char* name) {
  try {
    return std::stod(value);
  } catch (const std::exception&) {
    throw std::runtime_error(std::string("Unable to parse ") + name + " value: " + value);
  }
}

std::vector<WaveRow> load_wave_rows(const std::string& path) {
  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    throw std::runtime_error("Unable to open wave samples CSV: " + path);
  }

  std::string line;
  if (!std::getline(ifs, line)) {
    throw std::runtime_error("Wave samples CSV is empty: " + path);
  }

  const std::vector<std::string> header = split_csv_line(line);
  std::size_t time_idx = std::numeric_limits<std::size_t>::max();
  std::size_t disp_z_idx = std::numeric_limits<std::size_t>::max();
  for (std::size_t i = 0; i < header.size(); ++i) {
    if (header[i] == "time") {
      time_idx = i;
    } else if (header[i] == "disp_z") {
      disp_z_idx = i;
    }
  }

  if (time_idx == std::numeric_limits<std::size_t>::max() || disp_z_idx == std::numeric_limits<std::size_t>::max()) {
    throw std::runtime_error("Wave CSV missing required columns 'time' and 'disp_z': " + path);
  }

  std::vector<WaveRow> rows;
  while (std::getline(ifs, line)) {
    if (trim(line).empty()) {
      continue;
    }

    const std::vector<std::string> fields = split_csv_line(line);
    if (fields.size() <= std::max(time_idx, disp_z_idx)) {
      throw std::runtime_error("Malformed row in wave CSV: " + path);
    }

    WaveRow row;
    row.time_s = parse_double_field(fields[time_idx], "time");
    row.disp_z_m = parse_double_field(fields[disp_z_idx], "disp_z");
    rows.push_back(row);
  }

  if (rows.empty()) {
    throw std::runtime_error("No wave samples found in CSV: " + path);
  }

  return rows;
}

double drift_signal(const Scenario& scenario, double wave_freq_hz, double time_s, double duration_s) {
  const double primary_freq_hz = std::max(0.0035, wave_freq_hz * 0.085);
  const double secondary_freq_hz = std::max(0.0020, primary_freq_hz * 0.41);
  const double primary_amp_m = 0.38 * scenario.height_m;
  const double secondary_amp_m = 0.14 * scenario.height_m;
  const double ramp_amp_m = 0.20 * scenario.height_m;
  const double normalized_time = (time_s / duration_s) - 0.5;

  return primary_amp_m * std::sin((2.0 * kPi * primary_freq_hz * time_s) + 0.35) +
         secondary_amp_m * std::sin((2.0 * kPi * secondary_freq_hz * time_s) - 0.80) +
         ramp_amp_m * normalized_time;
}

std::string resolve_output_path(int argc, char* argv[], int arg_index, const char* fallback) {
  return (argc > arg_index) ? argv[arg_index] : std::string(fallback);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    const std::string samples_path = resolve_output_path(argc, argv, 1, "adaptive_wave_detrender_sim_output.csv");
    const std::string summary_path = resolve_output_path(argc, argv, 2, "adaptive_wave_detrender_sim_summary.csv");

    std::ofstream samples_ofs(samples_path);
    if (!samples_ofs.is_open()) {
      throw std::runtime_error("Unable to open output CSV: " + samples_path);
    }

    std::ofstream summary_ofs(summary_path);
    if (!summary_ofs.is_open()) {
      throw std::runtime_error("Unable to open summary CSV: " + summary_path);
    }

    samples_ofs << "scenario,height_m,period_s,wave_freq_hz,time_s,reference_z_m,drift_m,drifted_z_m,detrended_z_m,baseline_slow_m,error_z_m\n";
    samples_ofs << std::fixed << std::setprecision(9);

    summary_ofs << "scenario,height_m,period_s,wave_freq_hz,warmup_s,drift_rms_m,detrended_rms_m,gate_rms_m,improvement_ratio\n";
    summary_ofs << std::fixed << std::setprecision(9);

    bool all_ok = true;

    for (const Scenario& scenario : kScenarios) {
      const std::vector<WaveRow> wave_rows = load_wave_rows(scenario.wave_samples_csv);
      const double wave_freq_hz = 1.0 / scenario.period_s;
      const double warmup_s = std::max(kWarmupMinS, 6.0 * scenario.period_s);
      const double gate_rms_m = kRmsGateFractionOfHeight * scenario.height_m;
      const double duration_s = wave_rows.back().time_s;

      AdaptiveWaveDetrender detrender;
      RMSAccumulator drift_rms;
      RMSAccumulator detrended_rms;

      bool initialized = false;
      double prev_time_s = 0.0;

      for (const WaveRow& row : wave_rows) {
        if (initialized && !(row.time_s > prev_time_s)) {
          continue;
        }

        const double drift_m = drift_signal(scenario, wave_freq_hz, row.time_s, duration_s);
        const double drifted_z_m = row.disp_z_m + drift_m;

        AdaptiveWaveDetrender::Output output;
        if (!initialized) {
          detrender.reset(static_cast<float>(drifted_z_m));
          output = detrender.lastOutput();
          initialized = true;
        } else {
          const double dt_s = row.time_s - prev_time_s;
          output = detrender.update(static_cast<float>(drifted_z_m),
                                    static_cast<float>(dt_s),
                                    static_cast<float>(wave_freq_hz),
                                    true);
        }

        const double detrended_z_m = static_cast<double>(output.wave_clean);
        const double error_z_m = detrended_z_m - row.disp_z_m;

        if (row.time_s >= warmup_s) {
          drift_rms.add(drifted_z_m - row.disp_z_m);
          detrended_rms.add(error_z_m);
        }

        samples_ofs << scenario.tag << ','
                    << scenario.height_m << ','
                    << scenario.period_s << ','
                    << wave_freq_hz << ','
                    << row.time_s << ','
                    << row.disp_z_m << ','
                    << drift_m << ','
                    << drifted_z_m << ','
                    << detrended_z_m << ','
                    << static_cast<double>(output.baseline_slow) << ','
                    << error_z_m << '\n';

        prev_time_s = row.time_s;
      }

      const double drift_rms_m = drift_rms.rms();
      const double detrended_rms_m = detrended_rms.rms();
      const double improvement_ratio = (detrended_rms_m > 0.0)
          ? (drift_rms_m / detrended_rms_m)
          : std::numeric_limits<double>::infinity();

      summary_ofs << scenario.tag << ','
                  << scenario.height_m << ','
                  << scenario.period_s << ','
                  << wave_freq_hz << ','
                  << warmup_s << ','
                  << drift_rms_m << ','
                  << detrended_rms_m << ','
                  << gate_rms_m << ','
                  << improvement_ratio << '\n';

      std::cout << scenario.tag
                << ": drift RMS=" << drift_rms_m
                << " m, detrended RMS=" << detrended_rms_m
                << " m, gate=" << gate_rms_m
                << " m, improvement=" << improvement_ratio << "x\n";

      if (!(detrended_rms_m <= gate_rms_m) || !(detrended_rms_m < drift_rms_m)) {
        all_ok = false;
        std::cerr << "Quality gate failed for scenario '" << scenario.tag << "'\n";
      }
    }

    std::cout << "Wrote simulation samples to " << samples_path << '\n';
    std::cout << "Wrote simulation summary to " << summary_path << '\n';

    return all_ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& ex) {
    std::cerr << "detrend-wave-test failed: " << ex.what() << '\n';
    return EXIT_FAILURE;
  }
}
