#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <random>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO
#include "wave_dir/KalmanWaveDirection.h"
#include "wave_dir/WaveDirectionDetector.h"

namespace {

constexpr float kPi = std::numbers::pi_v<float>;

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "FAIL: " << message << "\n";
  std::exit(EXIT_FAILURE);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

float wrap360(float deg) {
  deg = std::fmod(deg, 360.0f);
  if (deg < 0.0f) deg += 360.0f;
  return deg;
}

float directed_error_deg(float estimate, float reference) {
  return std::abs(wrap360(estimate - reference + 180.0f) - 180.0f);
}

float axial_error_deg(float estimate, float reference) {
  const float directed = directed_error_deg(estimate, reference);
  return std::min(directed, 180.0f - directed);
}

Eigen::Vector2f axis_from_deg(float deg) {
  const float rad = deg * kPi / 180.0f;
  return Eigen::Vector2f(std::cos(rad), std::sin(rad));
}

struct AxisRun {
  float angle_deg = NAN;
  float confidence = 0.0f;
  float uncertainty_deg = NAN;
  float linearity = 0.0f;
  float amplitude = 0.0f;
};

AxisRun run_axis(float angle_deg,
                 float phase_offset,
                 float dt,
                 float noise_sigma = 0.012f) {
  constexpr float frequency_hz = 0.45f;
  const float omega = 2.0f * kPi * frequency_hz;
  const Eigen::Vector2f axis = axis_from_deg(angle_deg);

  KalmanWaveDirection filter(omega, 0.01f);
  filter.setMeasurementNoise(0.0025f);
  filter.setProcessNoise(1e-7f);

  const unsigned seed = static_cast<unsigned>(
      1001 + int(angle_deg * 19.0f) + int(phase_offset * 1000.0f) +
      int(1.0f / dt));
  std::mt19937 rng(seed);
  std::normal_distribution<float> noise(0.0f, noise_sigma);

  const int samples = int(50.0f / dt);
  for (int i = 0; i < samples; ++i) {
    const float t = float(i) * dt;
    const float carrier = 0.75f * std::cos(omega * t + phase_offset);
    filter.update(axis.x() * carrier + noise(rng),
                  axis.y() * carrier + noise(rng),
                  omega, dt);
  }

  return {
      filter.getAxisDegrees(),
      filter.getConfidence(),
      filter.getAxisUncertaintyDegrees(),
      filter.getLastStableLinearity(),
      filter.getAmplitude()
  };
}

void test_all_axes_and_carrier_phases() {
  const std::vector<float> phases{
      0.0f, 0.27f, 0.5f * kPi, 1.37f, 2.42f
  };

  for (int angle = 0; angle < 180; angle += 15) {
    for (float phase : phases) {
      const AxisRun result = run_axis(
          float(angle), phase, 1.0f / 200.0f);
      const float error = axial_error_deg(result.angle_deg, float(angle));

      require(error < 0.75f,
              "axis error at " + std::to_string(angle) +
              " deg and carrier phase " + std::to_string(phase) +
              " was " + std::to_string(error));
      require(result.confidence > 100.0f,
              "axis confidence did not converge at carrier phase " +
              std::to_string(phase));
      require(result.uncertainty_deg < 1.5f,
              "axis uncertainty too large at carrier phase " +
              std::to_string(phase));
      require(result.linearity > 0.98f,
              "linearly polarized wave was not recognized as axial");
      require(std::abs(result.amplitude - 0.75f) < 0.03f,
              "phase-invariant amplitude estimate was inaccurate");
    }
  }
}

void test_axis_sample_rate_invariance() {
  for (float angle : {7.0f, 89.0f, 143.0f}) {
    for (float phase : {0.5f * kPi, 1.17f}) {
      const AxisRun r100 = run_axis(
          angle, phase, 1.0f / 100.0f, 0.006f);
      const AxisRun r200 = run_axis(
          angle, phase, 1.0f / 200.0f, 0.006f);
      const AxisRun r400 = run_axis(
          angle, phase, 1.0f / 400.0f, 0.006f);

      require(axial_error_deg(r100.angle_deg, angle) < 0.35f,
              "100 Hz axis estimate missed truth");
      require(axial_error_deg(r200.angle_deg, angle) < 0.35f,
              "200 Hz axis estimate missed truth");
      require(axial_error_deg(r400.angle_deg, angle) < 0.35f,
              "400 Hz axis estimate missed truth");
      require(axial_error_deg(r100.angle_deg, r200.angle_deg) < 0.25f,
              "100/200 Hz axis mismatch");
      require(axial_error_deg(r200.angle_deg, r400.angle_deg) < 0.25f,
              "200/400 Hz axis mismatch");
    }
  }
}

void test_circular_motion_has_no_axis() {
  constexpr float dt = 1.0f / 200.0f;
  constexpr float frequency_hz = 0.45f;
  constexpr float omega = 2.0f * kPi * frequency_hz;

  KalmanWaveDirection filter(omega, 0.01f);
  filter.setMeasurementNoise(0.001f);
  filter.setProcessNoise(1e-7f);

  const int samples = int(40.0f / dt);
  for (int i = 0; i < samples; ++i) {
    const float phase = omega * float(i) * dt;
    filter.update(0.65f * std::cos(phase),
                  0.65f * std::sin(phase),
                  omega, dt);
  }

  require(filter.getAxisLinearity() < 0.05f,
          "circular horizontal motion was misreported as an axis");
  require(filter.getConfidence() < 20.0f,
          "circular horizontal motion retained axis confidence");
  require(filter.getAxis().isZero(1e-7f),
          "circular horizontal motion exposed a usable axis");
  require(!std::isfinite(filter.getAxisDegrees()),
          "circular horizontal motion exposed a finite axis angle");
  require(filter.getHistoricalStableConfidence() < 20.0f,
          "circular horizontal motion created a stable axis history");
}

void test_current_confidence_withdraws_stale_axis() {
  constexpr float dt = 1.0f / 200.0f;
  constexpr float frequency_hz = 0.45f;
  constexpr float omega = 2.0f * kPi * frequency_hz;
  const Eigen::Vector2f physical_axis = axis_from_deg(35.0f);

  KalmanWaveDirection axis_filter(omega, 0.01f);
  axis_filter.setMeasurementNoise(0.001f);
  axis_filter.setProcessNoise(1e-7f);

  WaveDirectionDetector<float>::Config sense_config;
  sense_config.smoothing_time_constant_sec = 0.8f;
  sense_config.vertical_slope_time_constant_sec = 0.025f;
  sense_config.coherence_threshold_on = 0.35f;
  sense_config.coherence_threshold_off = 0.25f;
  sense_config.absolute_product_floor = 1e-4f;
  sense_config.min_horizontal_rms = 0.03f;
  sense_config.min_vertical_slope_rms = 0.03f;
  sense_config.min_axis_confidence = 20.0f;
  sense_config.convention_sign = -1.0f;
  WaveDirectionDetector<float> sense_filter(sense_config);

  const int line_samples = int(35.0f / dt);
  for (int i = 0; i < line_samples; ++i) {
    const float phase = omega * float(i) * dt;
    // Linear deep-water phase relation for positive-axis propagation.
    const float along = -0.65f * std::cos(phase);
    const float vertical_up = 0.50f * std::sin(phase);
    const float ax = physical_axis.x() * along;
    const float ay = physical_axis.y() * along;

    axis_filter.update(ax, ay, omega, dt);
    const Eigen::Vector2f representative = axis_filter.getAxis();
    sense_filter.update(ax, ay, vertical_up,
                        representative.x(), representative.y(), dt,
                        axis_filter.getConfidence());
  }

  require(axis_filter.getConfidence() > 100.0f,
          "linear wave did not establish current axis confidence");
  require(axis_filter.getHistoricalStableConfidence() > 100.0f,
          "linear wave did not establish stable axis history");
  require(sense_filter.getState() == FORWARD,
          "linear wave sense did not converge before transition");

  const int circular_samples = int(25.0f / dt);
  for (int i = 0; i < circular_samples; ++i) {
    const float phase = omega * float(i + line_samples) * dt;
    const float ax = 0.65f * std::cos(phase);
    const float ay = 0.65f * std::sin(phase);
    const float vertical_up = 0.50f * std::sin(phase);

    axis_filter.update(ax, ay, omega, dt);
    const Eigen::Vector2f representative = axis_filter.getAxis();
    sense_filter.update(ax, ay, vertical_up,
                        representative.x(), representative.y(), dt,
                        axis_filter.getConfidence());
  }

  require(axis_filter.getConfidence() < 20.0f,
          "current confidence did not reject non-axial motion");
  require(axis_filter.getAxis().isZero(1e-7f),
          "invalid current motion still exposed a usable axis");
  require(!std::isfinite(axis_filter.getAxisDegrees()),
          "invalid current motion still exposed a finite axis angle");
  require(axis_filter.getHistoricalStableConfidence() > 20.0f,
          "transition discarded all previously valid axis history");
  require(sense_filter.getState() == UNCERTAIN,
          "sense detector was not withdrawn with current axis validity");
}

}  // namespace

int main() {
  test_all_axes_and_carrier_phases();
  test_axis_sample_rate_invariance();
  test_circular_motion_has_no_axis();
  test_current_confidence_withdraws_stale_axis();
  std::cout << "All I/Q wave-axis tests passed.\n";
  return 0;
}
