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
#include "detrend/AdaptiveWaveDetrender3D.h"

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

struct Vec3d {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct RMSAccumulator3D {
  Vec3d sum_sq{};
  std::size_t count = 0U;

  void add(const Vec3d& v) {
    sum_sq.x += v.x * v.x;
    sum_sq.y += v.y * v.y;
    sum_sq.z += v.z * v.z;
    ++count;
  }

  Vec3d rms() const {
    if (count == 0U) {
      return Vec3d{};
    }
    const double inv_n = 1.0 / static_cast<double>(count);
    return Vec3d{std::sqrt(sum_sq.x * inv_n), std::sqrt(sum_sq.y * inv_n), std::sqrt(sum_sq.z * inv_n)};
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

Vec3d reference_signal(const Scenario& scenario, double wave_freq_hz, double time_s, double ref_z_m) {
  const double omega = 2.0 * kPi * wave_freq_hz;
  Vec3d ref;
  ref.z = ref_z_m;
  ref.x = (0.35 * scenario.height_m) * std::sin((omega * time_s) + 0.40) +
          (0.08 * scenario.height_m) * std::sin((0.50 * omega * time_s) - 0.70);
  ref.y = (0.28 * scenario.height_m) * std::cos((omega * time_s) - 0.30) +
          (0.06 * scenario.height_m) * std::sin((0.65 * omega * time_s) + 0.20);
  return ref;
}

Vec3d drift_signal(const Scenario& scenario, double wave_freq_hz, double time_s, double duration_s) {
  const double primary_freq_hz = std::max(0.0035, wave_freq_hz * 0.085);
  const double secondary_freq_hz = std::max(0.0020, primary_freq_hz * 0.41);
  const double normalized_time = (time_s / duration_s) - 0.5;

  Vec3d drift;
  drift.x = (0.16 * scenario.height_m) * std::sin((2.0 * kPi * primary_freq_hz * time_s) + 0.10) +
            (0.05 * scenario.height_m) * std::sin((2.0 * kPi * secondary_freq_hz * time_s) - 0.40) +
            (0.08 * scenario.height_m) * normalized_time;
  drift.y = (0.13 * scenario.height_m) * std::sin((2.0 * kPi * primary_freq_hz * time_s) - 0.80) +
            (0.06 * scenario.height_m) * std::sin((2.0 * kPi * secondary_freq_hz * time_s) + 0.35) +
            (0.07 * scenario.height_m) * normalized_time;
  drift.z = (0.38 * scenario.height_m) * std::sin((2.0 * kPi * primary_freq_hz * time_s) + 0.35) +
            (0.14 * scenario.height_m) * std::sin((2.0 * kPi * secondary_freq_hz * time_s) - 0.80) +
            (0.20 * scenario.height_m) * normalized_time;
  return drift;
}

std::string resolve_output_path(int argc, char* argv[], int arg_index, const char* fallback) {
  return (argc > arg_index) ? argv[arg_index] : std::string(fallback);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    const std::string samples_path = resolve_output_path(argc, argv, 1, "adaptive_wave_detrender3d_sim_output.csv");
    const std::string summary_path = resolve_output_path(argc, argv, 2, "adaptive_wave_detrender3d_sim_summary.csv");

    std::ofstream samples_ofs(samples_path);
    if (!samples_ofs.is_open()) {
      throw std::runtime_error("Unable to open output CSV: " + samples_path);
    }

    std::ofstream summary_ofs(summary_path);
    if (!summary_ofs.is_open()) {
      throw std::runtime_error("Unable to open summary CSV: " + summary_path);
    }

    samples_ofs << "scenario,height_m,period_s,wave_freq_hz,time_s,"
                   "reference_x_m,drift_x_m,drifted_x_m,detrended_x_m,baseline_x_m,error_x_m,"
                   "reference_y_m,drift_y_m,drifted_y_m,detrended_y_m,baseline_y_m,error_y_m,"
                   "reference_z_m,drift_z_m,drifted_z_m,detrended_z_m,baseline_z_m,error_z_m\n";
    samples_ofs << std::fixed << std::setprecision(9);

    summary_ofs << "scenario,height_m,period_s,wave_freq_hz,warmup_s,"
                   "drift_rms_x_m,detrended_rms_x_m,gate_rms_x_m,improvement_ratio_x,"
                   "drift_rms_y_m,detrended_rms_y_m,gate_rms_y_m,improvement_ratio_y,"
                   "drift_rms_z_m,detrended_rms_z_m,gate_rms_z_m,improvement_ratio_z\n";
    summary_ofs << std::fixed << std::setprecision(9);

    bool all_ok = true;

    for (const Scenario& scenario : kScenarios) {
      const std::vector<WaveRow> wave_rows = load_wave_rows(scenario.wave_samples_csv);
      const double wave_freq_hz = 1.0 / scenario.period_s;
      const double warmup_s = std::max(kWarmupMinS, 6.0 * scenario.period_s);
      const double gate_rms_m = kRmsGateFractionOfHeight * scenario.height_m;
      const double duration_s = wave_rows.back().time_s;

      AdaptiveWaveDetrender3D::Config cfg;
      cfg.baseline_cutoff_fraction = 0.18f;
      cfg.enable_wave_cleanup = false;
      AdaptiveWaveDetrender3D detrender(cfg);
      RMSAccumulator3D drift_rms;
      RMSAccumulator3D detrended_rms;

      bool initialized = false;
      double prev_time_s = 0.0;

      for (const WaveRow& row : wave_rows) {
        if (initialized && !(row.time_s > prev_time_s)) {
          continue;
        }

        const Vec3d reference = reference_signal(scenario, wave_freq_hz, row.time_s, row.disp_z_m);
        const Vec3d drift = drift_signal(scenario, wave_freq_hz, row.time_s, duration_s);
        const Vec3d drifted{reference.x + drift.x, reference.y + drift.y, reference.z + drift.z};

        AdaptiveWaveDetrender3D::Output output;
        if (!initialized) {
          detrender.reset(static_cast<float>(drifted.x),
                         static_cast<float>(drifted.y),
                         static_cast<float>(drifted.z));
          output = detrender.lastOutput();
          initialized = true;
        } else {
          const double dt_s = row.time_s - prev_time_s;
          output = detrender.update(static_cast<float>(drifted.x),
                                    static_cast<float>(drifted.y),
                                    static_cast<float>(drifted.z),
                                    static_cast<float>(dt_s),
                                    static_cast<float>(wave_freq_hz),
                                    true);
        }

        const Vec3d detrended{static_cast<double>(output.wave_clean.x()),
                              static_cast<double>(output.wave_clean.y()),
                              static_cast<double>(output.wave_clean.z())};
        const Vec3d baseline{static_cast<double>(output.baseline_slow.x()),
                             static_cast<double>(output.baseline_slow.y()),
                             static_cast<double>(output.baseline_slow.z())};
        const Vec3d error{detrended.x - reference.x,
                          detrended.y - reference.y,
                          detrended.z - reference.z};

        if (row.time_s >= warmup_s) {
          drift_rms.add(drift);
          detrended_rms.add(error);
        }

        samples_ofs << scenario.tag << ','
                    << scenario.height_m << ','
                    << scenario.period_s << ','
                    << wave_freq_hz << ','
                    << row.time_s << ','
                    << reference.x << ','
                    << drift.x << ','
                    << drifted.x << ','
                    << detrended.x << ','
                    << baseline.x << ','
                    << error.x << ','
                    << reference.y << ','
                    << drift.y << ','
                    << drifted.y << ','
                    << detrended.y << ','
                    << baseline.y << ','
                    << error.y << ','
                    << reference.z << ','
                    << drift.z << ','
                    << drifted.z << ','
                    << detrended.z << ','
                    << baseline.z << ','
                    << error.z << '\n';

        prev_time_s = row.time_s;
      }

      const Vec3d drift_rms_v = drift_rms.rms();
      const Vec3d detrended_rms_v = detrended_rms.rms();

      const double improvement_x = (detrended_rms_v.x > 0.0)
          ? (drift_rms_v.x / detrended_rms_v.x)
          : std::numeric_limits<double>::infinity();
      const double improvement_y = (detrended_rms_v.y > 0.0)
          ? (drift_rms_v.y / detrended_rms_v.y)
          : std::numeric_limits<double>::infinity();
      const double improvement_z = (detrended_rms_v.z > 0.0)
          ? (drift_rms_v.z / detrended_rms_v.z)
          : std::numeric_limits<double>::infinity();

      summary_ofs << scenario.tag << ','
                  << scenario.height_m << ','
                  << scenario.period_s << ','
                  << wave_freq_hz << ','
                  << warmup_s << ','
                  << drift_rms_v.x << ','
                  << detrended_rms_v.x << ','
                  << gate_rms_m << ','
                  << improvement_x << ','
                  << drift_rms_v.y << ','
                  << detrended_rms_v.y << ','
                  << gate_rms_m << ','
                  << improvement_y << ','
                  << drift_rms_v.z << ','
                  << detrended_rms_v.z << ','
                  << gate_rms_m << ','
                  << improvement_z << '\n';

      std::cout << scenario.tag
                << ": RMS x (drift/detrended/gate) = " << drift_rms_v.x << '/' << detrended_rms_v.x << '/' << gate_rms_m
                << ", y = " << drift_rms_v.y << '/' << detrended_rms_v.y << '/' << gate_rms_m
                << ", z = " << drift_rms_v.z << '/' << detrended_rms_v.z << '/' << gate_rms_m << '\n';

      const bool ok_x = (detrended_rms_v.x <= gate_rms_m) && (detrended_rms_v.x < drift_rms_v.x);
      const bool ok_y = (detrended_rms_v.y <= gate_rms_m) && (detrended_rms_v.y < drift_rms_v.y);
      const bool ok_z = (detrended_rms_v.z <= gate_rms_m) && (detrended_rms_v.z < drift_rms_v.z);
      if (!(ok_x && ok_y && ok_z)) {
        all_ok = false;
        std::cerr << "Quality gate failed for scenario '" << scenario.tag << "'\n";
      }
    }

    std::cout << "Wrote simulation samples to " << samples_path << '\n';
    std::cout << "Wrote simulation summary to " << summary_path << '\n';

    return all_ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& ex) {
    std::cerr << "detrend-wave3d-test failed: " << ex.what() << '\n';
    return EXIT_FAILURE;
  }
}
