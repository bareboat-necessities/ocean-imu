#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>

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
  test_filter_is_not_copyable();
  std::cout << "All time-aware 3D spike filter tests passed.\n";
  return 0;
}
