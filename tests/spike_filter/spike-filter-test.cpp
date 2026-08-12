#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "avg/TimeAwareSpikeFilter.h"

namespace {

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "FAIL: " << message << "\n";
  std::exit(EXIT_FAILURE);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

// A clean linear ramp sampled at a constant rate has one and only one
// derivative, so no sample may ever be flagged as a spike. This is the
// regression guard for the circular-index defect: comparing the newest sample
// against the oldest one (rather than against its own predecessor) produces a
// large bogus derivative on every step and flags a steady signal as spiky.
void test_linear_ramp_passes_through_unchanged() {
  constexpr float dt = 0.1f;
  constexpr float rate = 1.0f;
  TimeAwareSpikeFilter filter(5, 2.0f);

  for (int k = 0; k < 200; k++) {
    const float input = rate * dt * static_cast<float>(k);
    const float output = filter.filterWithDelta(input, dt);
    require(output == input,
            "linear ramp sample " + std::to_string(k) +
                " was altered (expected " + std::to_string(input) +
                ", got " + std::to_string(output) + ")");
  }
}

// A constant signal is the degenerate case of the same contract.
void test_constant_signal_passes_through_unchanged() {
  TimeAwareSpikeFilter filter(4, 0.5f);
  for (int k = 0; k < 50; k++) {
    const float output = filter.filterWithDelta(3.25f, 0.02f);
    require(output == 3.25f, "constant signal was altered at sample " +
                                 std::to_string(k));
  }
}

// An isolated outlier must be rejected on the sample where it occurs, and the
// sequence must recover to pass-through shortly afterwards.
void test_isolated_spike_is_rejected() {
  constexpr float dt = 0.1f;
  constexpr float rate = 1.0f;
  constexpr int spike_index = 20;
  constexpr float spike_amplitude = 5.0f;
  TimeAwareSpikeFilter filter(5, 2.0f);

  float max_error = 0.0f;
  for (int k = 0; k < 40; k++) {
    const float truth = rate * dt * static_cast<float>(k);
    const float input =
        (k == spike_index) ? truth + spike_amplitude : truth;
    const float output = filter.filterWithDelta(input, dt);

    if (k == spike_index) {
      require(output != input, "spike sample was passed through unfiltered");
      require(std::fabs(output - input) > 0.5f * spike_amplitude,
              "spike sample was only partially attenuated");
    }
    const float error = std::fabs(output - truth);
    if (error > max_error) max_error = error;

    // The window is five samples deep, so the outlier can no longer influence
    // the output once it has been pushed out of the buffer.
    if (k > spike_index + 5) {
      require(output == input,
              "filter did not recover to pass-through at sample " +
                  std::to_string(k));
    }
  }

  require(max_error < 0.5f,
          "isolated spike leaked into the output (max error " +
              std::to_string(max_error) + ")");
}

// The filter is time aware: what matters is the derivative, not the raw value
// step. A ramp sampled at an irregular cadence has a constant derivative and
// must survive intact even though individual value steps differ by 4x.
void test_irregular_sampling_is_not_mistaken_for_spikes() {
  constexpr int window_size = 5;
  constexpr float rate = 2.0f;
  const std::vector<float> deltas = {0.05f, 0.2f, 0.05f, 0.1f, 0.2f};
  TimeAwareSpikeFilter filter(window_size, 1.0f);

  float value = 0.0f;
  for (int k = 0; k < 100; k++) {
    const float dt = deltas[static_cast<std::size_t>(k) % deltas.size()];
    value += rate * dt;
    const float output = filter.filterWithDelta(value, dt);
    // Samples before the window is full of real history are covered by
    // test_startup_transient_settles_within_one_window.
    if (k < window_size) continue;
    require(output == value,
            "irregularly sampled ramp was altered at sample " +
                std::to_string(k));
  }
}

// Priming fills the window with copies of the first sample, so the fabricated
// history has a zero derivative. A signal whose slope exceeds the threshold is
// therefore flagged while that fabricated history is still in the buffer. This
// pins the documented startup behaviour: the transient must clear once the
// window holds only real samples.
void test_startup_transient_settles_within_one_window() {
  constexpr int window_size = 5;
  constexpr float dt = 0.1f;
  constexpr float rate = 4.0f;  // slope well above the threshold below
  TimeAwareSpikeFilter filter(window_size, 1.0f);

  for (int k = 0; k < 60; k++) {
    const float input = rate * dt * static_cast<float>(k);
    const float output = filter.filterWithDelta(input, dt);
    if (k < window_size) continue;
    require(output == input,
            "startup transient had not cleared by sample " +
                std::to_string(k));
  }
}

// A value step that is large in magnitude but consistent with the elapsed time
// is a legitimate signal, not a spike.
void test_large_step_over_long_gap_is_not_a_spike() {
  constexpr float rate = 1.0f;
  TimeAwareSpikeFilter filter(5, 2.0f);

  float value = 0.0f;
  for (int k = 0; k < 20; k++) {
    value += rate * 0.1f;
    filter.filterWithDelta(value, 0.1f);
  }

  // Ten times the usual gap, so ten times the usual value step.
  value += rate * 1.0f;
  const float output = filter.filterWithDelta(value, 1.0f);
  require(output == value,
          "a value step proportional to a long time gap was flagged as a spike");
}

// The first call primes the window and must return the raw input.
void test_first_sample_is_returned_verbatim() {
  TimeAwareSpikeFilter filter(5, 0.001f);
  const float output = filter.filterWithDelta(42.0f, 0.1f);
  require(output == 42.0f, "first sample was not returned verbatim");
}

// Window sizes below two cannot produce a derivative; the constructor clamps
// them so the median never reads outside its buffer.
void test_degenerate_window_size_is_clamped() {
  for (int size : {0, 1, 2}) {
    TimeAwareSpikeFilter filter(size, 1.0f);
    for (int k = 0; k < 10; k++) {
      const float input = 0.1f * static_cast<float>(k);
      const float output = filter.filterWithDelta(input, 0.1f);
      require(std::isfinite(output),
              "degenerate window size " + std::to_string(size) +
                  " produced a non-finite output");
    }
  }
}

// The filter owns raw arrays. Copying it would alias that storage and
// double-free on destruction, so the copy operations are deleted.
void test_filter_is_not_copyable() {
  static_assert(!std::is_copy_constructible<TimeAwareSpikeFilter>::value,
                "TimeAwareSpikeFilter must not be copy constructible");
  static_assert(!std::is_copy_assignable<TimeAwareSpikeFilter>::value,
                "TimeAwareSpikeFilter must not be copy assignable");
}

}  // namespace

int main() {
  test_linear_ramp_passes_through_unchanged();
  test_constant_signal_passes_through_unchanged();
  test_isolated_spike_is_rejected();
  test_irregular_sampling_is_not_mistaken_for_spikes();
  test_startup_transient_settles_within_one_window();
  test_large_step_over_long_gap_is_not_a_spike();
  test_first_sample_is_returned_verbatim();
  test_degenerate_window_size_is_clamped();
  test_filter_is_not_copyable();
  std::cout << "All time-aware spike filter tests passed.\n";
  return 0;
}
