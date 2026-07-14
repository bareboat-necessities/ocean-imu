#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one occurrence, found {count}")
    return text.replace(old, new, 1)


test_path = ROOT / "tests/wave_dir/wave-direction-test.cpp"
test = test_path.read_text(encoding="utf-8")
start_marker = "void test_axis_estimator_all_angles() {\n"
end_marker = "struct SenseRun {\n"
start = test.find(start_marker)
end = test.find(end_marker)
if start < 0 or end < 0 or end <= start:
    raise RuntimeError("axis test block markers not found exactly")

axis_tests = r'''struct AxisRun {
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

'''

test = test[:start] + axis_tests + test[end:]
test = replace_once(
    test,
    "  test_axis_estimator_all_angles();\n",
    "  test_axis_estimator_all_angles_and_phases();\n"
    "  test_axis_estimator_sample_rate_invariance();\n"
    "  test_axis_estimator_rejects_circular_motion();\n",
    "main axis test calls")
test_path.write_text(test, encoding="utf-8")

for relative in [
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
]:
    path = ROOT / relative
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r"\n        // Keep BODY components around for direction/sign\.?\n"
        r"        const float a_x_body = acc\.x\(\);\n"
        r"        const float a_y_body = acc\.y\(\);\n")
    text, count = pattern.subn("\n", text, count=1)
    if count != 1:
        raise RuntimeError(f"{relative}: unused body-axis block not found exactly")
    text = text.replace(
        "BODY-Z-based proxy used by the tracker/sign logic.",
        "BODY-Z-based proxy used by the tracker/tuner logic.")
    text = text.replace(
        "Up-positive BODY-Z proxy used by tracker/tuner/sign logic.",
        "Up-positive BODY-Z proxy used by tracker/tuner logic.")
    path.write_text(text, encoding="utf-8")

print("Applied phase-degeneracy tests and direction integration cleanup.")
