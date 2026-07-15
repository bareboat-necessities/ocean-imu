#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <string>

#include "wave_dir/WaveEncounter.h"

namespace {

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "FAIL: " << message << "\n";
  std::exit(EXIT_FAILURE);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

void test_forward_model_round_trip() {
  using namespace wave_encounter;
  constexpr double g = 9.80665;
  const double intrinsic_omega = 2.0 * std::numbers::pi * 0.42;
  const double k = deep_water_wavenumber(intrinsic_omega, g);

  for (int intrinsic_sense : {-1, 1}) {
    for (double vessel_speed : {-4.0, 0.0, 3.0, 8.0, 14.0}) {
      const double omega_e = signed_encounter_omega(
          intrinsic_omega, k, intrinsic_sense, vessel_speed);
      const int measured_sense = apparent_sense(intrinsic_sense, omega_e);
      require(measured_sense != 0,
              "round-trip case selected zero encounter frequency");

      const auto solutions = solve_deep_water(
          std::abs(omega_e), vessel_speed, measured_sense,
          g, 0.1, 12.0);
      bool found_truth = false;
      for (const auto& solution : solutions) {
        if (solution.intrinsic_sense == intrinsic_sense &&
            std::abs(solution.intrinsic_omega_rad_s - intrinsic_omega) < 1e-10) {
          found_truth = true;
        }
      }
      require(found_truth,
              "analytic encounter inversion omitted the generating solution");
    }
  }
}

void test_stationary_observer_has_unique_direction() {
  using namespace wave_encounter;
  constexpr double encounter = 2.3;

  const auto positive = solve_deep_water(
      encounter, 0.0, 1, 9.80665, 0.1, 10.0);
  require(positive.size() == 1,
          "stationary positive apparent sense was not unique");
  require(positive.front().intrinsic_sense == 1,
          "stationary positive sense mapped to wrong intrinsic direction");
  require(std::abs(positive.front().intrinsic_omega_rad_s - encounter) < 1e-12,
          "stationary encounter frequency did not equal intrinsic frequency");

  const auto negative = solve_deep_water(
      encounter, 0.0, -1, 9.80665, 0.1, 10.0);
  require(negative.size() == 1,
          "stationary negative apparent sense was not unique");
  require(negative.front().intrinsic_sense == -1,
          "stationary negative sense mapped to wrong intrinsic direction");
}

void test_following_sea_returns_both_valid_roots() {
  using namespace wave_encounter;
  constexpr double g = 9.80665;
  constexpr double vessel_speed = 5.0;
  constexpr double encounter = 0.25;

  const auto solutions = solve_deep_water(
      encounter, vessel_speed, 1, g, 0.01, 10.0);
  require(solutions.size() == 2,
          "following-sea ambiguity did not return both intrinsic roots");

  for (const auto& solution : solutions) {
    require(solution.intrinsic_sense == 1,
            "following-sea ambiguity returned wrong intrinsic sense");
    require(solution.signed_encounter_omega_rad_s > 0.0,
            "following-sea root has inconsistent encounter branch");
    require(std::abs(solution.signed_encounter_omega_rad_s - encounter) < 1e-12,
            "following-sea root does not reproduce encounter frequency");
  }
  require(solutions[0].intrinsic_omega_rad_s < solutions[1].intrinsic_omega_rad_s,
          "following-sea roots were not sorted");
}

void test_tangent_double_root_is_not_missed() {
  using namespace wave_encounter;
  constexpr double g = 9.80665;
  constexpr double vessel_speed = 5.0;
  const double encounter_max = g / (4.0 * vessel_speed);
  const double tangent_omega = g / (2.0 * vessel_speed);

  const auto solutions = solve_deep_water(
      encounter_max, vessel_speed, 1, g, 0.01, 10.0);
  require(solutions.size() == 1,
          "tangent encounter root was missed or duplicated");
  require(solutions.front().intrinsic_sense == 1,
          "tangent root has wrong intrinsic sense");
  require(std::abs(solutions.front().intrinsic_omega_rad_s - tangent_omega) < 1e-10,
          "tangent root frequency is inaccurate");
}

void test_phase_velocity_over_ground() {
  using namespace wave_encounter;
  constexpr double g = 9.80665;
  const double omega = 2.0 * std::numbers::pi * 0.42;
  const double k = deep_water_wavenumber(omega, g);
  const Velocity2<double> current{0.7, -0.2};

  const auto velocity = phase_velocity_over_ground(
      omega, k, 0.6, 0.8, current);
  const double phase_speed = omega / k;
  require(std::abs(velocity.x - (0.6 * phase_speed + current.x)) < 1e-12,
          "ground phase velocity X mismatch");
  require(std::abs(velocity.y - (0.8 * phase_speed + current.y)) < 1e-12,
          "ground phase velocity Y mismatch");
}

}  // namespace

int main() {
  test_forward_model_round_trip();
  test_stationary_observer_has_unique_direction();
  test_following_sea_returns_both_valid_roots();
  test_tangent_double_root_is_not_missed();
  test_phase_velocity_over_ground();
  std::cout << "All wave encounter tests passed.\n";
  return 0;
}
