#include <algorithm>
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
#include "wave_dir/WaveDirectionFrame.h"
#include "wave_dir/WaveEncounter.h"

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
  float d = wrap360(estimate - reference + 180.0f) - 180.0f;
  return std::abs(d);
}

float axial_error_deg(float estimate, float reference) {
  float d = directed_error_deg(estimate, reference);
  return std::min(d, 180.0f - d);
}

Eigen::Vector2f axis_from_deg(float deg) {
  const float r = deg * kPi / 180.0f;
  return Eigen::Vector2f(std::cos(r), std::sin(r));
}

WaveDirectionDetector<float>::Config strict_sense_config() {
  WaveDirectionDetector<float>::Config cfg;
  cfg.smoothing_time_constant_sec = 0.8f;
  cfg.vertical_slope_time_constant_sec = 0.025f;
  cfg.coherence_threshold_on = 0.35f;
  cfg.coherence_threshold_off = 0.25f;
  cfg.absolute_product_floor = 1e-4f;
  cfg.min_horizontal_rms = 0.03f;
  cfg.min_vertical_slope_rms = 0.03f;
  cfg.min_axis_confidence = 0.0f;
  cfg.convention_sign = 1.0f;
  return cfg;
}

struct AxisRun {
  float angle_deg = NAN;
  float confidence = 0.0f;
  float uncertainty_deg = NAN;
  float linearity = 0.0f;
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
      filter.getLastStableConfidence(),
      filter.getAxisUncertaintyDegrees(),
      filter.getLastStableLinearity()
  };
}

void test_axis_estimator_all_angles_and_phases() {
  const std::vector<float> phases{
      0.0f, 0.27f, 0.5f * kPi, 1.37f, 2.42f
  };

  for (int angle = 0; angle < 180; angle += 15) {
    for (float phase : phases) {
      const AxisRun result = run_axis(
          float(angle), phase, 1.0f / 200.0f);
      const float error = axial_error_deg(result.angle_deg, float(angle));
      require(error < 0.75f,
              "axis estimator error at " + std::to_string(angle) +
              " deg and phase " + std::to_string(phase) +
              " was " + std::to_string(error));
      require(result.confidence > 100.0f,
              "axis confidence did not converge at phase " +
              std::to_string(phase));
      require(result.uncertainty_deg < 1.5f,
              "axis uncertainty too large at phase " +
              std::to_string(phase));
      require(result.linearity > 0.98f,
              "linearly polarized wave was not recognized as axial");
    }
  }
}

void test_axis_estimator_sample_rate_invariance() {
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

void test_axis_estimator_rejects_circular_motion() {
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
  require(filter.getLastStableConfidence() < 20.0f,
          "circular horizontal motion produced a stable axis");
}

struct SenseRun {
  WaveDirection state = UNCERTAIN;
  float coherence = 0.0f;
  float apparent_to_deg = NAN;
  float classified_fraction = 0.0f;
};

SenseRun run_sense(float angle_deg,
                   int true_sense,
                   float dt,
                   float phase_offset,
                   bool flip_axis_representative,
                   float noise_sigma = 0.01f) {
  constexpr float frequency_hz = 0.55f;
  const float omega = 2.0f * kPi * frequency_hz;
  const Eigen::Vector2f axis = axis_from_deg(angle_deg);
  WaveDirectionDetector<float> detector(strict_sense_config());

  std::mt19937 rng(unsigned(2000 + int(angle_deg) * 11 + (true_sense + 1) * 31));
  std::normal_distribution<float> noise(0.0f, noise_sigma);

  const int samples = int(32.0f / dt);
  const int score_start = int(18.0f / dt);
  int correct = 0;
  int scored = 0;

  for (int i = 0; i < samples; ++i) {
    const float t = float(i) * dt;
    const float phase = omega * t + phase_offset;
    const float horizontal = float(true_sense) * 0.65f * std::cos(phase);
    const float vertical_up = 0.50f * std::sin(phase);

    Eigen::Vector2f representative = axis;
    if (flip_axis_representative && ((i / std::max(1, int(0.73f / dt))) % 2)) {
      representative = -representative;
    }

    const WaveDirection state = detector.update(
        axis.x() * horizontal + noise(rng),
        axis.y() * horizontal + noise(rng),
        vertical_up + noise(rng),
        representative.x(), representative.y(), dt, 1000.0f);

    if (i >= score_start) {
      ++scored;
      if ((true_sense > 0 && state == FORWARD) ||
          (true_sense < 0 && state == BACKWARD)) {
        ++correct;
      }
    }
  }

  SenseRun result;
  result.state = detector.getState();
  result.coherence = detector.getCoherence();
  result.apparent_to_deg = detector.getDirectedAngleDegrees();
  result.classified_fraction = scored ? float(correct) / float(scored) : 0.0f;
  return result;
}

void test_sense_detector_all_angles_and_phases() {
  const std::vector<float> phases{0.0f, 0.61f, 1.37f, 2.42f};
  for (int angle = 0; angle < 180; angle += 5) {
    for (int sense : {-1, 1}) {
      for (float phase : phases) {
        const SenseRun result = run_sense(float(angle), sense, 1.0f / 200.0f,
                                          phase, false);
        const float expected = wrap360(float(angle) + (sense < 0 ? 180.0f : 0.0f));
        require(result.state == (sense > 0 ? FORWARD : BACKWARD),
                "wrong sense at axis " + std::to_string(angle));
        require(result.classified_fraction > 0.995f,
                "sense classification below 99.5% at axis " + std::to_string(angle));
        require(std::abs(result.coherence) > 0.90f,
                "weak phase coherence at axis " + std::to_string(angle));
        require(directed_error_deg(result.apparent_to_deg, expected) < 0.05f,
                "wrong apparent angle at axis " + std::to_string(angle));
      }
    }
  }
}

void test_sample_rate_invariance() {
  for (float angle : {0.0f, 37.0f, 89.0f, 143.0f, 175.0f}) {
    for (int sense : {-1, 1}) {
      const SenseRun r100 = run_sense(angle, sense, 1.0f / 100.0f, 0.83f, false, 0.006f);
      const SenseRun r200 = run_sense(angle, sense, 1.0f / 200.0f, 0.83f, false, 0.006f);
      const SenseRun r400 = run_sense(angle, sense, 1.0f / 400.0f, 0.83f, false, 0.006f);
      require(r100.state == r200.state && r200.state == r400.state,
              "sample-rate-dependent sense decision");
      require(std::abs(r100.coherence - r200.coherence) < 0.025f,
              "100/200 Hz coherence mismatch");
      require(std::abs(r200.coherence - r400.coherence) < 0.025f,
              "200/400 Hz coherence mismatch");
    }
  }
}

void test_axis_representative_flip_invariance() {
  for (float angle : {3.0f, 44.0f, 91.0f, 136.0f, 178.0f}) {
    for (int sense : {-1, 1}) {
      const SenseRun stable = run_sense(angle, sense, 1.0f / 200.0f, 0.41f, false);
      const SenseRun flipped = run_sense(angle, sense, 1.0f / 200.0f, 0.41f, true);
      require(flipped.classified_fraction > 0.995f,
              "axis representative flips caused decision loss");
      require(directed_error_deg(stable.apparent_to_deg, flipped.apparent_to_deg) < 0.05f,
              "axis representative flip changed directed angle");
    }
  }
}

void test_combined_axis_and_sense() {
  constexpr float dt = 1.0f / 200.0f;
  constexpr float frequency_hz = 0.50f;
  constexpr float omega = 2.0f * kPi * frequency_hz;
  constexpr int samples = int(55.0f / dt);

  for (int angle = 0; angle < 180; angle += 15) {
    const Eigen::Vector2f physical_axis = axis_from_deg(float(angle));
    for (int sense : {-1, 1}) {
      KalmanWaveDirection axis_filter(omega, 0.01f);
      axis_filter.setMeasurementNoise(0.001f);
      axis_filter.setProcessNoise(1e-7f);
      WaveDirectionDetector<float> sense_filter(strict_sense_config());

      std::mt19937 rng(unsigned(3000 + angle * 13 + sense));
      std::normal_distribution<float> noise(0.0f, 0.008f);

      for (int i = 0; i < samples; ++i) {
        const float phase = omega * float(i) * dt + 0.22f;
        const float horizontal = float(sense) * 0.70f * std::cos(phase);
        const float vertical_up = 0.52f * std::sin(phase);
        const float ax = physical_axis.x() * horizontal + noise(rng);
        const float ay = physical_axis.y() * horizontal + noise(rng);
        axis_filter.update(ax, ay, omega, dt);
        const Eigen::Vector2f representative = axis_filter.getAxis();
        sense_filter.update(ax, ay, vertical_up + noise(rng),
                            representative.x(), representative.y(), dt,
                            axis_filter.getLastStableConfidence());
      }

      const float expected = wrap360(float(angle) + (sense < 0 ? 180.0f : 0.0f));
      const float estimate = sense_filter.getDirectedAngleDegrees();
      require(std::isfinite(estimate), "combined estimator remained uncertain");
      require(directed_error_deg(estimate, expected) < 0.80f,
              "combined directed error at axis " + std::to_string(angle) +
              " was " + std::to_string(directed_error_deg(estimate, expected)));
    }
  }
}


void test_heading_frame_roll_pitch_yaw_invariance() {
  constexpr float g = 9.80665f;
  constexpr float deg = kPi / 180.0f;

  const std::vector<Eigen::Vector3f> attitudes_deg{
      {0.0f, 0.0f, 0.0f},
      {25.0f, -18.0f, 37.0f},
      {-42.0f, 31.0f, -123.0f},
      {65.0f, -35.0f, 178.0f},
  };
  const std::vector<Eigen::Vector3f> heading_accels{
      {0.55f, -0.21f, 0.37f},
      {-0.44f, 0.62f, -0.19f},
      {0.08f, 0.91f, 0.12f},
  };

  for (const auto& euler : attitudes_deg) {
    const Eigen::AngleAxisf yaw(euler.z() * deg, Eigen::Vector3f::UnitZ());
    const Eigen::AngleAxisf pitch(euler.y() * deg, Eigen::Vector3f::UnitY());
    const Eigen::AngleAxisf roll(euler.x() * deg, Eigen::Vector3f::UnitX());
    const Eigen::Quaternionf q_body_to_world = yaw * pitch * roll;

    const Eigen::Vector3f bow_world = q_body_to_world * Eigen::Vector3f::UnitX();
    const Eigen::Vector2f heading =
        Eigen::Vector2f(bow_world.x(), bow_world.y()).normalized();
    const Eigen::Vector2f starboard(-heading.y(), heading.x());

    for (const auto& expected : heading_accels) {
      const Eigen::Vector3f acceleration_world(
          heading.x() * expected.x() + starboard.x() * expected.y(),
          heading.y() * expected.x() + starboard.y() * expected.y(),
          -expected.z());
      const Eigen::Vector3f specific_force_body =
          q_body_to_world.conjugate() *
          (acceleration_world - Eigen::Vector3f(0.0f, 0.0f, g));

      const auto recovered = wave_direction::heading_frame_acceleration<float>(
          q_body_to_world, specific_force_body, g);
      require(recovered.heading_valid, "heading-frame conversion rejected valid attitude");
      require(std::abs(recovered.forward_ms2 - expected.x()) < 2e-5f,
              "heading-frame forward acceleration mismatch");
      require(std::abs(recovered.starboard_ms2 - expected.y()) < 2e-5f,
              "heading-frame starboard acceleration mismatch");
      require(std::abs(recovered.up_ms2 - expected.z()) < 2e-5f,
              "heading-frame vertical acceleration mismatch");
    }
  }
}

void test_encounter_forward_model_and_inverse() {
  using namespace wave_encounter;
  constexpr double g = 9.80665;
  const double intrinsic_omega = 2.0 * std::numbers::pi * 0.42;
  const double k = deep_water_wavenumber(intrinsic_omega, g);

  for (int intrinsic_sense : {-1, 1}) {
    for (double vessel_speed : {-4.0, 0.0, 3.0, 8.0, 14.0}) {
      const double omega_e = signed_encounter_omega(
          intrinsic_omega, k, intrinsic_sense, vessel_speed);
      const int measured_sense = apparent_sense(intrinsic_sense, omega_e);
      require(measured_sense != 0, "test selected a zero encounter frequency");

      const auto solutions = solve_deep_water(
          std::abs(omega_e), vessel_speed, measured_sense,
          g, 0.1, 12.0, 12000);
      bool found_truth = false;
      for (const auto& solution : solutions) {
        if (solution.intrinsic_sense == intrinsic_sense &&
            std::abs(solution.intrinsic_omega_rad_s - intrinsic_omega) < 1e-5) {
          found_truth = true;
        }
      }
      require(found_truth,
              "deep-water encounter inversion omitted the generating solution");
    }
  }

  const Velocity2<double> current{0.7, -0.2};
  const auto velocity = phase_velocity_over_ground(
      intrinsic_omega, k, 0.6, 0.8, current);
  const double cp = intrinsic_omega / k;
  require(std::abs(velocity.x - (0.6 * cp + current.x)) < 1e-12,
          "ground phase velocity X mismatch");
  require(std::abs(velocity.y - (0.8 * cp + current.y)) < 1e-12,
          "ground phase velocity Y mismatch");
}

}  // namespace

int main() {
  test_axis_estimator_all_angles_and_phases();
  test_axis_estimator_sample_rate_invariance();
  test_axis_estimator_rejects_circular_motion();
  test_sense_detector_all_angles_and_phases();
  test_sample_rate_invariance();
  test_axis_representative_flip_invariance();
  test_combined_axis_and_sense();
  test_heading_frame_roll_pitch_yaw_invariance();
  test_encounter_forward_model_and_inverse();
  std::cout << "All wave propagation direction tests passed.\n";
  return 0;
}
