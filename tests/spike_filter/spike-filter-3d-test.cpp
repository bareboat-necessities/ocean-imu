#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>

#include "AtomS3R/AtomS3R_ImuSpikeFilter.h"
#include "avg/TimeAwareSpikeFilter3D.h"

namespace {

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "FAIL: " << message << "\n";
  std::exit(EXIT_FAILURE);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

constexpr float kDt = 0.005f;  // 200 Hz
constexpr float kG = 9.81f;
constexpr float kThreshold = 200.0f;  // m/s^3

// A rotating gravity vector has constant magnitude, so a pure attitude change
// must never be flagged regardless of how fast the direction moves.
void test_rotation_at_constant_magnitude_passes_through() {
  TimeAwareSpikeFilter3D filter(7, kThreshold);
  for (int k = 0; k < 400; k++) {
    const float ang = 0.05f * static_cast<float>(k);
    float x = kG * std::sin(ang), y = 0.0f, z = kG * std::cos(ang);
    const float x0 = x, y0 = y, z0 = z;
    require(!filter.filterWithDelta(x, y, z, kDt),
            "constant-magnitude rotation flagged at sample " + std::to_string(k));
    require(x == x0 && y == y0 && z == z0, "constant-magnitude sample altered");
  }
}

// A smooth heave acceleration (well under the jerk threshold) passes through.
void test_smooth_wave_passes_through() {
  TimeAwareSpikeFilter3D filter(7, kThreshold);
  for (int k = 0; k < 2000; k++) {
    const float t = kDt * static_cast<float>(k);
    float x = 0.0f, y = 0.0f, z = kG + 3.0f * std::sin(2.0f * 3.14159265f * 0.5f * t);
    require(!filter.filterWithDelta(x, y, z, kDt),
            "smooth wave flagged at sample " + std::to_string(k));
  }
}

// An isolated magnitude spike is replaced by an observed window vector, and
// the filter recovers to pass-through once the outlier leaves the window.
void test_isolated_spike_is_rejected() {
  constexpr int window = 7;
  constexpr int spike_index = 50;
  TimeAwareSpikeFilter3D filter(window, kThreshold);
  for (int k = 0; k < 100; k++) {
    float x = 0.1f, y = -0.2f, z = kG;
    if (k == spike_index) { x += 20.0f; y -= 15.0f; z += 30.0f; }
    const bool spike = filter.filterWithDelta(x, y, z, kDt);
    if (k == spike_index) {
      require(spike, "spike was not flagged");
      require(x == 0.1f && y == -0.2f && z == kG,
              "spike was not replaced by the median-magnitude vector");
    } else if (k > spike_index + window) {
      require(!spike, "filter did not recover at sample " + std::to_string(k));
    }
  }
}

// A dropout towards zero is a magnitude spike too.
void test_dropout_spike_is_rejected() {
  TimeAwareSpikeFilter3D filter(7, kThreshold);
  for (int k = 0; k < 20; k++) {
    float x = 0.0f, y = 0.0f, z = kG;
    filter.filterWithDelta(x, y, z, kDt);
  }
  float x = 0.0f, y = 0.0f, z = 0.0f;
  require(filter.filterWithDelta(x, y, z, kDt), "dropout was not flagged");
  require(z == kG, "dropout was not replaced");
}

// Non-finite input is passed through untouched and does not poison history.
void test_non_finite_input_is_ignored() {
  TimeAwareSpikeFilter3D filter(7, kThreshold);
  for (int k = 0; k < 20; k++) {
    float x = 0.0f, y = 0.0f, z = kG;
    filter.filterWithDelta(x, y, z, kDt);
  }
  float x = NAN, y = 0.0f, z = kG;
  require(!filter.filterWithDelta(x, y, z, kDt), "NaN sample flagged");
  require(std::isnan(x), "NaN sample was altered");
  float x2 = 0.0f, y2 = 0.0f, z2 = kG;
  require(!filter.filterWithDelta(x2, y2, z2, kDt), "history poisoned by NaN");
}

// Even windows replace with the mean of the two middle-magnitude vectors.
void test_even_window_and_degenerate_size() {
  for (int size : {0, 1, 2, 4}) {
    TimeAwareSpikeFilter3D filter(size, kThreshold);
    for (int k = 0; k < 20; k++) {
      float x = 0.0f, y = 0.0f, z = (k == 10) ? 100.0f : kG;
      filter.filterWithDelta(x, y, z, kDt);
      require(std::isfinite(x) && std::isfinite(y) && std::isfinite(z),
              "non-finite output for window size " + std::to_string(size));
    }
  }
}

// A vertical sign flip keeps |a| but must be caught when the vertical check
// is enabled, and must pass through (magnitude-only) when it is not.
void test_vertical_flip_is_rejected() {
  for (bool vertical_on : {true, false}) {
    TimeAwareSpikeFilter3D filter(7, kThreshold, vertical_on ? 2 : -1, kThreshold);
    for (int k = 0; k < 40; k++) {
      float x = 0.3f, y = -0.4f, z = (k == 20) ? -kG : kG;
      const bool spike = filter.filterWithDelta(x, y, z, kDt);
      if (k == 20) {
        require(spike == vertical_on, "vertical flip detection mismatch");
        if (vertical_on) require(z == kG, "vertical flip was not replaced");
      } else if (k > 27) {
        require(!spike, "filter did not recover after vertical flip");
      }
    }
  }
}

// Gravity moving between the horizontal axes (yaw-like rotation of a tilted
// sensor) does not change the vertical component and is never flagged.
void test_horizontal_rotation_passes_with_vertical_check() {
  TimeAwareSpikeFilter3D filter(7, kThreshold, 2, kThreshold);
  for (int k = 0; k < 400; k++) {
    const float ang = 0.2f * static_cast<float>(k);
    float x = 3.0f * std::cos(ang), y = 3.0f * std::sin(ang), z = 9.34f;
    require(!filter.filterWithDelta(x, y, z, kDt),
            "horizontal rotation flagged at sample " + std::to_string(k));
  }
}

// Slow roll moves gravity through the vertical axis well below the threshold.
void test_slow_roll_passes_with_vertical_check() {
  TimeAwareSpikeFilter3D filter(7, kThreshold, 2, kThreshold);
  for (int k = 0; k < 2000; k++) {
    const float t = kDt * static_cast<float>(k);
    const float roll = 0.5f * std::sin(2.0f * 3.14159265f * 0.3f * t);
    float x = 0.0f, y = kG * std::sin(roll), z = kG * std::cos(roll);
    require(!filter.filterWithDelta(x, y, z, kDt),
            "slow roll flagged at sample " + std::to_string(k));
  }
}

// The gyro read of a rejected accelerometer sample is replaced by the gyro of
// the same observed window sample; passed-through samples keep their gyro.
void test_companion_gyro_follows_replacement() {
  constexpr int spike_index = 30;
  TimeAwareSpikeFilter3D filter(5, kThreshold, 2, kThreshold);
  for (int k = 0; k < 60; k++) {
    float x = 0.1f, y = -0.2f, z = kG;
    float gx = 0.01f, gy = -0.02f, gz = 0.03f;
    if (k == spike_index) { z += 40.0f; gx = 5.0f; gy = -5.0f; gz = 5.0f; }
    const float gx0 = gx, gy0 = gy, gz0 = gz;
    const bool spike = filter.filterWithDelta(x, y, z, gx, gy, gz, kDt);
    if (k == spike_index) {
      require(spike, "spike with companion gyro was not flagged");
      require(z == kG, "accel spike was not replaced");
      require(gx == 0.01f && gy == -0.02f && gz == 0.03f,
              "companion gyro was not replaced with the chosen sample's gyro");
    } else if (!spike) {
      require(gx == gx0 && gy == gy0 && gz == gz0, "passed-through gyro altered");
    }
    if (k < spike_index || k > spike_index + 5) {
      require(!spike, "clean sample flagged at " + std::to_string(k));
    }
  }

  // A non-finite gyro read is passed through and not stored.
  float x = 0.1f, y = -0.2f, z = kG, gx = NAN, gy = 0.0f, gz = 0.0f;
  require(!filter.filterWithDelta(x, y, z, gx, gy, gz, kDt), "NaN gyro flagged");
  require(std::isnan(gx), "NaN gyro was altered");
}

// Even windows average the gyro of the same two middle samples as the accel.
void test_companion_gyro_even_window() {
  TimeAwareSpikeFilter3D filter(4, kThreshold);
  const float zs[] = {9.0f, 9.5f, 10.0f, 10.5f};
  const float gs[] = {1.0f, 2.0f, 3.0f, 4.0f};
  float x = 0.0f, y = 0.0f, z = zs[0], gx = gs[0], gy = 0.0f, gz = 0.0f;
  filter.filterWithDelta(x, y, z, gx, gy, gz, 1.0f);  // initializes window
  for (int k = 1; k < 4; k++) {
    x = 0.0f; y = 0.0f; z = zs[k]; gx = gs[k]; gy = 0.0f; gz = 0.0f;
    require(!filter.filterWithDelta(x, y, z, gx, gy, gz, 1.0f), "ramp flagged");
  }
  x = 0.0f; y = 0.0f; z = 500.0f; gx = 100.0f; gy = 0.0f; gz = 0.0f;
  require(filter.filterWithDelta(x, y, z, gx, gy, gz, 1.0f), "even-window spike not flagged");
  // Window magnitudes {9.5, 10, 10.5, 500}: middle pair 10 and 10.5.
  require(z == 10.25f && gx == 3.5f, "even-window gyro does not match accel pair");
}

// The shared AtomS3R wrapper replaces accel and gyro together.
void test_atoms3r_wrapper() {
  struct Vec3 {
    float v[3];
    float& x() { return v[0]; }
    float& y() { return v[1]; }
    float& z() { return v[2]; }
  };
  atoms3r_ical::ImuSpikeFilter filter;
  for (int k = 0; k < 20; k++) {
    Vec3 a{{0.0f, 0.0f, kG}}, w{{0.0f, 0.0f, 0.1f}};
    if (k == 10) { a.v[2] = -kG; w.v[2] = 9.0f; }
    const bool spike = filter.apply(a, w, kDt);
    if (k == 10) {
      require(spike == atoms3r_ical::ImuSpikeCfg::enabled, "wrapper spike mismatch");
      if (spike) require(a.v[2] == kG && w.v[2] == 0.1f, "wrapper did not replace accel+gyro");
    } else if (k < 10 || k > 10 + atoms3r_ical::ImuSpikeCfg::window) {
      require(!spike, "wrapper flagged a clean sample");
    }
  }
}

void test_filter_is_not_copyable() {
  static_assert(!std::is_copy_constructible<TimeAwareSpikeFilter3D>::value,
                "TimeAwareSpikeFilter3D must not be copy constructible");
  static_assert(!std::is_copy_assignable<TimeAwareSpikeFilter3D>::value,
                "TimeAwareSpikeFilter3D must not be copy assignable");
}

}  // namespace

int main() {
  test_rotation_at_constant_magnitude_passes_through();
  test_smooth_wave_passes_through();
  test_isolated_spike_is_rejected();
  test_dropout_spike_is_rejected();
  test_non_finite_input_is_ignored();
  test_even_window_and_degenerate_size();
  test_vertical_flip_is_rejected();
  test_horizontal_rotation_passes_with_vertical_check();
  test_slow_roll_passes_with_vertical_check();
  test_companion_gyro_follows_replacement();
  test_companion_gyro_even_window();
  test_atoms3r_wrapper();
  test_filter_is_not_copyable();
  std::cout << "All time-aware 3D spike filter tests passed.\n";
  return 0;
}
