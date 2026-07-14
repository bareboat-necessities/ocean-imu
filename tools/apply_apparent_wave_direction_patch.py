#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def write(path: str, content: str) -> None:
    p = ROOT / path
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content, encoding="utf-8")


def replace_once(path: str, old: str, new: str) -> None:
    p = ROOT / path
    text = p.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one occurrence, found {count}: {old[:100]!r}")
    p.write_text(text.replace(old, new, 1), encoding="utf-8")


wave_direction_detector = r'''#pragma once

/*
  Copyright 2025-2026, Mikhail Grushinskiy

  Resolve apparent wave-propagation sense along an independently estimated
  horizontal propagation axis.  The axis is undirected (modulo 180 degrees);
  this detector uses horizontal/vertical orbital phase to choose one of its
  two directions.
*/

#include <algorithm>
#include <cmath>
#include <limits>

enum WaveDirection {
  BACKWARD = -1,   // propagation along the negative representative axis
  UNCERTAIN = 0,
  FORWARD = 1      // propagation along the positive representative axis
};

template <typename Real = float>
class WaveDirectionDetector {
public:
  struct Config {
    // EMA time constant for cross power and channel energies.
    Real smoothing_time_constant_sec = Real(2.5);

    // Low-pass time constant applied to the finite-difference vertical slope.
    Real vertical_slope_time_constant_sec = Real(0.05);

    // Normalized on/off thresholds.  The lower threshold supplies hysteresis.
    Real coherence_threshold_on = Real(0.20);
    Real coherence_threshold_off = Real(0.12);

    // Absolute signal gates.  These prevent a high normalized correlation made
    // from nearly zero signals from being reported as a direction.
    Real absolute_product_floor = Real(0.005);
    Real min_horizontal_rms = Real(0.01);
    Real min_vertical_slope_rms = Real(0.01);
    Real min_axis_confidence = Real(20);

    // Maps the chosen vertical-axis convention to FORWARD/BACKWARD.  The
    // default is correct for an up-positive vertical signal in which
    // a_parallel and d(a_vertical)/dt have positive correlation for
    // propagation along the positive axis.
    Real convention_sign = Real(1);
  };

  // Backward-compatible constructor.  The old smoothing argument was an EMA
  // coefficient tuned for 200 Hz.  Convert it to an elapsed-time constant so
  // behavior remains invariant when the sample rate changes.
  WaveDirectionDetector(Real smoothing = Real(0.002),
                        Real sensitivity = Real(0.005)) {
    config_.smoothing_time_constant_sec = alpha_to_tau(smoothing, Real(0.005));
    config_.absolute_product_floor = std::max(Real(0), sensitivity);
    sanitize_config();
    reset();
  }

  explicit WaveDirectionDetector(const Config& config) : config_(config) {
    sanitize_config();
    reset();
  }

  void setConfig(const Config& config) {
    config_ = config;
    sanitize_config();
  }

  const Config& getConfig() const noexcept { return config_; }

  void reset() {
    have_vertical_ = false;
    have_axis_ = false;
    previous_vertical_ = Real(0);
    filtered_vertical_slope_ = Real(0);
    filtered_cross_power_ = Real(0);
    filtered_horizontal_energy_ = Real(0);
    filtered_vertical_slope_energy_ = Real(0);
    coherence_ = Real(0);
    along_axis_acceleration_ = Real(0);
    axis_x_ = Real(1);
    axis_y_ = Real(0);
    state_ = UNCERTAIN;
  }

  // Update using an independently estimated propagation-axis representative.
  // axis_x/axis_y are expressed in the same boat-horizontal frame as accelX/Y.
  // The returned sign is relative to the internally continuity-stabilized axis.
  WaveDirection update(Real accelX,
                       Real accelY,
                       Real accelVertical,
                       Real axisX,
                       Real axisY,
                       Real delta_t,
                       Real axisConfidence = std::numeric_limits<Real>::infinity()) {
    if (!(delta_t > Real(0)) || !finite(delta_t) ||
        !finite(accelX) || !finite(accelY) || !finite(accelVertical) ||
        !finite(axisX) || !finite(axisY)) {
      state_ = UNCERTAIN;
      return state_;
    }

    const Real axis_norm = std::hypot(axisX, axisY);
    if (!(axis_norm > Real(1e-8))) {
      previous_vertical_ = accelVertical;
      have_vertical_ = true;
      state_ = UNCERTAIN;
      return state_;
    }

    axisX /= axis_norm;
    axisY /= axis_norm;

    // An axial estimate is equivalent under d -> -d.  Keep one representative
    // continuous so that a harmless 180-degree representation change does not
    // reset or invert the final directed result.
    if (have_axis_ && axisX * axis_x_ + axisY * axis_y_ < Real(0)) {
      axisX = -axisX;
      axisY = -axisY;
    }
    axis_x_ = axisX;
    axis_y_ = axisY;
    have_axis_ = true;

    along_axis_acceleration_ = accelX * axis_x_ + accelY * axis_y_;

    if (!have_vertical_) {
      previous_vertical_ = accelVertical;
      have_vertical_ = true;
      state_ = UNCERTAIN;
      return state_;
    }

    const Real raw_vertical_slope = (accelVertical - previous_vertical_) / delta_t;
    previous_vertical_ = accelVertical;

    const Real slope_alpha = ema_alpha(delta_t, config_.vertical_slope_time_constant_sec);
    filtered_vertical_slope_ +=
        slope_alpha * (raw_vertical_slope - filtered_vertical_slope_);

    // Do not accumulate a phase decision against an axis which has not yet
    // reached the confidence required by the independent angle estimator.
    if (finite(axisConfidence) && axisConfidence < config_.min_axis_confidence) {
      state_ = UNCERTAIN;
      return state_;
    }

    const Real alpha = ema_alpha(delta_t, config_.smoothing_time_constant_sec);
    const Real product = along_axis_acceleration_ * filtered_vertical_slope_;

    filtered_cross_power_ += alpha * (product - filtered_cross_power_);
    filtered_horizontal_energy_ += alpha *
        (along_axis_acceleration_ * along_axis_acceleration_ - filtered_horizontal_energy_);
    filtered_vertical_slope_energy_ += alpha *
        (filtered_vertical_slope_ * filtered_vertical_slope_ - filtered_vertical_slope_energy_);

    const Real denom = std::sqrt(std::max(
        Real(0), filtered_horizontal_energy_ * filtered_vertical_slope_energy_));
    coherence_ = (denom > Real(1e-12))
        ? std::clamp(filtered_cross_power_ / denom, Real(-1), Real(1))
        : Real(0);

    const Real horizontal_rms = std::sqrt(std::max(Real(0), filtered_horizontal_energy_));
    const Real vertical_slope_rms =
        std::sqrt(std::max(Real(0), filtered_vertical_slope_energy_));

    if (horizontal_rms < config_.min_horizontal_rms ||
        vertical_slope_rms < config_.min_vertical_slope_rms ||
        std::abs(filtered_cross_power_) < config_.absolute_product_floor) {
      state_ = UNCERTAIN;
      return state_;
    }

    const Real score = config_.convention_sign * coherence_;
    const Real threshold = (state_ == UNCERTAIN)
        ? config_.coherence_threshold_on
        : config_.coherence_threshold_off;

    if (score > threshold) {
      state_ = FORWARD;
    } else if (score < -threshold) {
      state_ = BACKWARD;
    } else if (std::abs(score) < config_.coherence_threshold_off) {
      state_ = UNCERTAIN;
    }

    return state_;
  }

  // Legacy fixed-Y overload retained only for source compatibility.  New code
  // must pass the independently estimated propagation axis to the overload
  // above.
  WaveDirection update(Real accelX, Real accelY, Real accelVertical, Real delta_t) {
    return update(accelX, accelY, accelVertical,
                  Real(0), Real(1), delta_t,
                  std::numeric_limits<Real>::infinity());
  }

  WaveDirection getState() const noexcept { return state_; }
  Real getFilteredP() const noexcept { return filtered_cross_power_; }
  Real getCoherence() const noexcept { return coherence_; }
  Real getAlongAxisAcceleration() const noexcept { return along_axis_acceleration_; }
  Real getVerticalSlope() const noexcept { return filtered_vertical_slope_; }
  Real getAxisX() const noexcept { return axis_x_; }
  Real getAxisY() const noexcept { return axis_y_; }

  Real getHorizontalRms() const noexcept {
    return std::sqrt(std::max(Real(0), filtered_horizontal_energy_));
  }

  Real getVerticalSlopeRms() const noexcept {
    return std::sqrt(std::max(Real(0), filtered_vertical_slope_energy_));
  }

  Real getDirectedX() const noexcept {
    if (state_ == UNCERTAIN) return std::numeric_limits<Real>::quiet_NaN();
    return (state_ == FORWARD) ? axis_x_ : -axis_x_;
  }

  Real getDirectedY() const noexcept {
    if (state_ == UNCERTAIN) return std::numeric_limits<Real>::quiet_NaN();
    return (state_ == FORWARD) ? axis_y_ : -axis_y_;
  }

  // 0 degrees is boat +X (forward); positive angles rotate toward boat +Y
  // (starboard in the NED-compatible body convention used by the filters).
  Real getDirectedAngleDegrees() const noexcept {
    const Real x = getDirectedX();
    const Real y = getDirectedY();
    if (!finite(x) || !finite(y)) return std::numeric_limits<Real>::quiet_NaN();
    Real deg = std::atan2(y, x) * Real(180) / pi();
    if (deg < Real(0)) deg += Real(360);
    if (deg >= Real(360)) deg -= Real(360);
    return deg;
  }

  Real getWaveFromAngleDegrees() const noexcept {
    Real deg = getDirectedAngleDegrees();
    if (!finite(deg)) return deg;
    deg += Real(180);
    if (deg >= Real(360)) deg -= Real(360);
    return deg;
  }

private:
  static constexpr Real pi() noexcept {
    return Real(3.141592653589793238462643383279502884L);
  }

  static bool finite(Real value) noexcept {
    return std::isfinite(static_cast<double>(value));
  }

  static Real ema_alpha(Real dt, Real tau) noexcept {
    if (!(tau > Real(0))) return Real(1);
    return Real(1) - std::exp(-dt / tau);
  }

  static Real alpha_to_tau(Real alpha, Real reference_dt) noexcept {
    if (!(alpha > Real(0)) || alpha >= Real(1)) return Real(2.5);
    return -reference_dt / std::log(Real(1) - alpha);
  }

  void sanitize_config() noexcept {
    if (!(config_.smoothing_time_constant_sec > Real(0)))
      config_.smoothing_time_constant_sec = Real(2.5);
    if (!(config_.vertical_slope_time_constant_sec > Real(0)))
      config_.vertical_slope_time_constant_sec = Real(0.05);
    config_.coherence_threshold_on =
        std::clamp(config_.coherence_threshold_on, Real(0), Real(1));
    config_.coherence_threshold_off =
        std::clamp(config_.coherence_threshold_off, Real(0),
                   config_.coherence_threshold_on);
    config_.absolute_product_floor = std::max(Real(0), config_.absolute_product_floor);
    config_.min_horizontal_rms = std::max(Real(0), config_.min_horizontal_rms);
    config_.min_vertical_slope_rms =
        std::max(Real(0), config_.min_vertical_slope_rms);
    config_.min_axis_confidence = std::max(Real(0), config_.min_axis_confidence);
    config_.convention_sign = (config_.convention_sign < Real(0)) ? Real(-1) : Real(1);
  }

  Config config_{};
  bool have_vertical_ = false;
  bool have_axis_ = false;
  Real previous_vertical_ = Real(0);
  Real filtered_vertical_slope_ = Real(0);
  Real filtered_cross_power_ = Real(0);
  Real filtered_horizontal_energy_ = Real(0);
  Real filtered_vertical_slope_energy_ = Real(0);
  Real coherence_ = Real(0);
  Real along_axis_acceleration_ = Real(0);
  Real axis_x_ = Real(1);
  Real axis_y_ = Real(0);
  WaveDirection state_ = UNCERTAIN;
};
'''

wave_encounter = r'''#pragma once

/*
  Apparent/intrinsic wave encounter helpers.

  The IMU phase detector reports an apparent propagation sense along a measured
  axis.  A moving observer measures signed encounter frequency

      Omega_e = sigma - k * d dot V_rel,

  where sigma is intrinsic water-relative angular frequency, k is wavenumber,
  d is the intrinsic propagation unit vector, and V_rel is vessel velocity
  relative to the water.  If current is neglected, V_rel may be approximated by
  velocity over ground.  With known current U, V_rel = V_ground - U.
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace wave_encounter {

template <typename Real>
struct DeepWaterSolution {
  Real intrinsic_omega_rad_s = std::numeric_limits<Real>::quiet_NaN();
  Real wavenumber_rad_m = std::numeric_limits<Real>::quiet_NaN();
  int intrinsic_sense = 0;          // +/- relative to the supplied axis
  Real signed_encounter_omega_rad_s = std::numeric_limits<Real>::quiet_NaN();
};

template <typename Real>
inline Real deep_water_wavenumber(Real intrinsic_omega_rad_s,
                                  Real gravity_ms2 = Real(9.80665)) {
  return intrinsic_omega_rad_s * intrinsic_omega_rad_s / gravity_ms2;
}

template <typename Real>
inline Real signed_encounter_omega(Real intrinsic_omega_rad_s,
                                   Real wavenumber_rad_m,
                                   int intrinsic_sense,
                                   Real vessel_speed_along_positive_axis_ms) {
  const Real sense = intrinsic_sense >= 0 ? Real(1) : Real(-1);
  return intrinsic_omega_rad_s
      - wavenumber_rad_m * sense * vessel_speed_along_positive_axis_ms;
}

template <typename Real>
inline int apparent_sense(int intrinsic_sense, Real signed_encounter_omega_rad_s) {
  if (intrinsic_sense == 0 || signed_encounter_omega_rad_s == Real(0)) return 0;
  const int encounter_sign = signed_encounter_omega_rad_s > Real(0) ? 1 : -1;
  return intrinsic_sense * encounter_sign;
}

// Solve the no-current (or water-relative velocity) deep-water encounter
// equation.  Multiple candidates are returned rather than silently choosing a
// branch in following seas, where encounter-frequency inversion can be
// non-unique.
template <typename Real>
inline std::vector<DeepWaterSolution<Real>> solve_deep_water(
    Real encounter_omega_abs_rad_s,
    Real vessel_speed_along_positive_axis_ms,
    int measured_apparent_sense,
    Real gravity_ms2 = Real(9.80665),
    Real omega_min_rad_s = Real(0.05),
    Real omega_max_rad_s = Real(20),
    int scan_intervals = 4096) {
  std::vector<DeepWaterSolution<Real>> solutions;
  if (!(encounter_omega_abs_rad_s > Real(0)) ||
      !(gravity_ms2 > Real(0)) ||
      !(omega_max_rad_s > omega_min_rad_s) ||
      scan_intervals < 16 || measured_apparent_sense == 0) {
    return solutions;
  }

  auto residual = [&](Real omega, int sense) {
    const Real k = deep_water_wavenumber(omega, gravity_ms2);
    return std::abs(signed_encounter_omega(
               omega, k, sense, vessel_speed_along_positive_axis_ms))
        - encounter_omega_abs_rad_s;
  };

  for (int sense : {-1, 1}) {
    Real lo = omega_min_rad_s;
    Real flo = residual(lo, sense);
    for (int i = 1; i <= scan_intervals; ++i) {
      const Real t = Real(i) / Real(scan_intervals);
      const Real hi = omega_min_rad_s
          + t * (omega_max_rad_s - omega_min_rad_s);
      const Real fhi = residual(hi, sense);

      bool bracket = (flo == Real(0)) || (fhi == Real(0)) ||
                     ((flo < Real(0)) != (fhi < Real(0)));
      if (bracket) {
        Real a = lo;
        Real b = hi;
        Real fa = flo;
        for (int iteration = 0; iteration < 80; ++iteration) {
          const Real mid = (a + b) / Real(2);
          const Real fm = residual(mid, sense);
          if (std::abs(fm) < Real(1e-10)) {
            a = b = mid;
            break;
          }
          if ((fa < Real(0)) != (fm < Real(0))) {
            b = mid;
          } else {
            a = mid;
            fa = fm;
          }
        }
        const Real omega = (a + b) / Real(2);
        const Real k = deep_water_wavenumber(omega, gravity_ms2);
        const Real omega_e = signed_encounter_omega(
            omega, k, sense, vessel_speed_along_positive_axis_ms);

        if (apparent_sense(sense, omega_e) == measured_apparent_sense) {
          bool duplicate = false;
          for (const auto& existing : solutions) {
            if (existing.intrinsic_sense == sense &&
                std::abs(existing.intrinsic_omega_rad_s - omega) < Real(1e-5)) {
              duplicate = true;
              break;
            }
          }
          if (!duplicate) {
            solutions.push_back({omega, k, sense, omega_e});
          }
        }
      }
      lo = hi;
      flo = fhi;
    }
  }

  std::sort(solutions.begin(), solutions.end(),
            [](const auto& a, const auto& b) {
              return a.intrinsic_omega_rad_s < b.intrinsic_omega_rad_s;
            });
  return solutions;
}

template <typename Real>
struct Velocity2 {
  Real x = Real(0);
  Real y = Real(0);
};

template <typename Real>
inline Velocity2<Real> phase_velocity_over_ground(
    Real intrinsic_omega_rad_s,
    Real wavenumber_rad_m,
    Real direction_x,
    Real direction_y,
    Velocity2<Real> current_ms = {}) {
  const Real norm = std::hypot(direction_x, direction_y);
  if (!(norm > Real(0)) || !(wavenumber_rad_m > Real(0))) {
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    return {nan, nan};
  }
  direction_x /= norm;
  direction_y /= norm;
  const Real intrinsic_phase_speed = intrinsic_omega_rad_s / wavenumber_rad_m;
  return {
      intrinsic_phase_speed * direction_x + current_ms.x,
      intrinsic_phase_speed * direction_y + current_ms.y
  };
}

}  // namespace wave_encounter
'''

wave_direction_test = r'''#include <algorithm>
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

void test_axis_estimator_all_angles() {
  constexpr float dt = 1.0f / 200.0f;
  constexpr float frequency_hz = 0.45f;
  constexpr float omega = 2.0f * kPi * frequency_hz;
  constexpr int samples = int(50.0f / dt);

  std::mt19937 rng(1001);
  std::normal_distribution<float> noise(0.0f, 0.012f);

  for (int angle = 0; angle < 180; angle += 15) {
    const Eigen::Vector2f axis = axis_from_deg(float(angle));
    KalmanWaveDirection filter(omega, 0.01f);
    filter.setMeasurementNoise(0.0025f);
    filter.setProcessNoise(1e-7f);

    for (int i = 0; i < samples; ++i) {
      const float t = float(i) * dt;
      const float carrier = 0.75f * std::cos(omega * t + 0.27f);
      filter.update(axis.x() * carrier + noise(rng),
                    axis.y() * carrier + noise(rng),
                    omega, dt);
    }

    const float error = axial_error_deg(filter.getAxisDegrees(), float(angle));
    require(error < 0.75f,
            "axis estimator error at " + std::to_string(angle) +
            " deg was " + std::to_string(error));
    require(filter.getLastStableConfidence() > 100.0f,
            "axis confidence did not converge at " + std::to_string(angle));
    require(filter.getAxisUncertaintyDegrees() < 1.5f,
            "axis uncertainty too large at " + std::to_string(angle));
  }
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
  test_axis_estimator_all_angles();
  test_sense_detector_all_angles_and_phases();
  test_sample_rate_invariance();
  test_axis_representative_flip_invariance();
  test_combined_axis_and_sense();
  test_encounter_forward_model_and_inverse();
  std::cout << "All wave propagation direction tests passed.\n";
  return 0;
}
'''

wave_dir_makefile = r'''.PHONY: all build test clean

CXX ?= g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-I$(EIGEN_DIR),-I/usr/include/eigen3)

BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Werror -pedantic -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif

TARGET = wave-direction-test

all: build test
build: $(TARGET)
test: $(TARGET)
	./run_tests.sh

$(TARGET): wave-direction-test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGET) $(TARGET).exe *.o
'''

wave_dir_run_tests = r'''#!/bin/bash -e

./wave-direction-test
'''

write("src/wave_dir/WaveDirectionDetector.h", wave_direction_detector)
write("src/wave_dir/WaveEncounter.h", wave_encounter)
write("tests/wave_dir/wave-direction-test.cpp", wave_direction_test)
write("tests/wave_dir/Makefile", wave_dir_makefile)
write("tests/wave_dir/run_tests.sh", wave_dir_run_tests)

# Make the angle estimator's axial semantics explicit while preserving API compatibility.
replace_once(
    "src/wave_dir/KalmanWaveDirection.h",
    "  Kalman filter for estimating direction of an ocean wave from IMU horizontal x, y accelerations.\n",
    "  Kalman filter for estimating the horizontal wave-propagation axis from\n"
    "  IMU horizontal x/y accelerations.  The result is axial (modulo 180\n"
    "  degrees); propagation sense is resolved separately.\n")

replace_once(
    "src/wave_dir/KalmanWaveDirection.h",
    '''    // Estimated wave propagation direction (unit vector)\n    Eigen::Vector2f getDirection() const {\n        return lastStableDir;\n    }\n\n    float getDirectionDegrees() const {\n        return direction_deg_smoothed;\n    }\n\n    float getDirectionDegreesRaw() const {\n        return direction_deg_raw;\n    }\n''',
    '''    // Continuity-stabilized representative of the propagation axis.\n    // Both d and -d describe the same vertical propagation plane.\n    Eigen::Vector2f getAxis() const {\n        return lastStableDir;\n    }\n\n    float getAxisDegrees() const {\n        return direction_deg_smoothed;\n    }\n\n    float getAxisDegreesRaw() const {\n        return direction_deg_raw;\n    }\n\n    // Backward-compatible aliases.  These return an axis modulo 180 degrees,\n    // not a fully directed apparent propagation angle.\n    Eigen::Vector2f getDirection() const { return getAxis(); }\n    float getDirectionDegrees() const { return getAxisDegrees(); }\n    float getDirectionDegreesRaw() const { return getAxisDegreesRaw(); }\n''')

replace_once(
    "src/wave_dir/KalmanWaveDirection.h",
    "    float getDirectionUncertaintyDegrees() const {\n",
    "    float getAxisUncertaintyDegrees() const {\n")
replace_once(
    "src/wave_dir/KalmanWaveDirection.h",
    "    float getLastStableConfidence() const {\n",
    "    // Backward-compatible uncertainty alias.\n"
    "    float getDirectionUncertaintyDegrees() const {\n"
    "        return getAxisUncertaintyDegrees();\n"
    "    }\n\n"
    "    float getLastStableConfidence() const {\n")

# Wire the independently estimated axis into both OU filter families.
old_direction_update = '''        // Direction filters run on BODY accel; sign uses the same BODY-Z proxy.\n        dir_filter_.update(a_x_body, a_y_body, omega, dt);\n        dir_sign_state_ = dir_sign_.update(a_x_body, a_y_body, a_body_z_up_proxy_, dt);\n'''
new_direction_update = '''        // Stage 1 estimates the apparent propagation plane as an unsigned axis\n        // relative to boat +X.  Stage 2 resolves propagation sense along that\n        // same axis from horizontal/vertical orbital phase.\n        dir_filter_.update(a_x_body, a_y_body, omega, dt);\n        const Eigen::Vector2f propagation_axis_body = dir_filter_.getAxis();\n        dir_sign_state_ = dir_sign_.update(\n            a_x_body, a_y_body, a_body_z_up_proxy_,\n            propagation_axis_body.x(), propagation_axis_body.y(),\n            dt, dir_filter_.getLastStableConfidence());\n'''

for path in [
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
]:
    replace_once(path, old_direction_update, new_direction_update)
    replace_once(
        path,
        '''    inline WaveDirection getDirSignState() const noexcept { return dir_sign_state_; }\n    inline float getWaveDirectionDeg() const noexcept { return dir_filter_.getDirectionDegrees(); }\n''',
        '''    inline WaveDirection getDirSignState() const noexcept { return dir_sign_state_; }\n\n    // Propagation-plane angle relative to boat +X, modulo 180 degrees.\n    inline float getWaveAxisDeg() const noexcept { return dir_filter_.getAxisDegrees(); }\n    inline float getWaveDirectionDeg() const noexcept { return getWaveAxisDeg(); }\n\n    // Fully directed apparent propagation angles observed by the moving boat.\n    // These are encounter/apparent directions unless vessel-motion correction\n    // is applied externally (see wave_dir/WaveEncounter.h).\n    inline float getApparentWaveDirectionToDeg() const noexcept {\n        return dir_sign_.getDirectedAngleDegrees();\n    }\n    inline float getApparentWaveDirectionFromDeg() const noexcept {\n        return dir_sign_.getWaveFromAngleDegrees();\n    }\n    inline float getDirSenseCoherence() const noexcept {\n        return dir_sign_.getCoherence();\n    }\n''')
    replace_once(
        path,
        '''        dir_filter_      = KalmanWaveDirection(2.0f * static_cast<float>(M_PI) * FREQ_GUESS);\n        dir_sign_state_  = UNCERTAIN;\n''',
        '''        dir_filter_      = KalmanWaveDirection(2.0f * static_cast<float>(M_PI) * FREQ_GUESS);\n        dir_sign_.reset();\n        dir_sign_state_  = UNCERTAIN;\n''')

# Extend simulation telemetry with explicit apparent directed outputs.
replace_once(
    "src/util/W3dSimCommon.h",
    '''    float direction_deg = NAN;\n    float direction_deg_generator_signed = NAN;\n''',
    '''    float direction_deg = NAN;              // propagation axis, modulo 180 deg\n    float apparent_to_deg = NAN;             // directed encounter propagation\n    float apparent_from_deg = NAN;           // opposite marine "waves from" angle\n    float sense_coherence = NAN;\n    float direction_deg_generator_signed = NAN;\n''')

for path in [
    "tests/kalman_ou_ii/kalman_ou_ii-sim.cpp",
    "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp",
]:
    replace_once(
        path,
        '''        s.direction.direction_deg = d.getDirectionDegrees();\n        s.direction.direction_deg_generator_signed = dirDegGeneratorSignedFromVec(d.getDirection());\n''',
        '''        s.direction.direction_deg = d.getAxisDegrees();\n        s.direction.apparent_to_deg = filter.getApparentWaveDirectionToDeg();\n        s.direction.apparent_from_deg = filter.getApparentWaveDirectionFromDeg();\n        s.direction.sense_coherence = filter.getDirSenseCoherence();\n        s.direction.direction_deg_generator_signed = dirDegGeneratorSignedFromVec(d.getAxis());\n''')
    replace_once(
        path,
        "        s.direction.uncertainty_deg = d.getDirectionUncertaintyDegrees();\n",
        "        s.direction.uncertainty_deg = d.getAxisUncertaintyDegrees();\n")
    replace_once(
        path,
        "        s.direction.direction_vec = d.getDirection();\n",
        "        s.direction.direction_vec = d.getAxis();\n")

replace_once(
    "src/util/W3dSimCommon.cpp",
    '''        << "dir_deg,dir_uncert_deg,dir_conf,dir_amp,"\n        << "dir_sign,dir_sign_num,"\n''',
    '''        << "dir_axis_deg,dir_apparent_to_deg,dir_apparent_from_deg,"\n        << "dir_sense_coherence,dir_uncert_deg,dir_conf,dir_amp,"\n        << "dir_sign,dir_sign_num,"\n''')

replace_once(
    "src/util/W3dSimCommon.cpp",
    '''            << snap.direction.phase << "," << snap.direction.direction_deg << "," << snap.direction.uncertainty_deg << ","\n            << snap.direction.confidence << "," << snap.direction.amplitude << ","\n            << (snap.direction.sign == FORWARD ? "TOWARD" : snap.direction.sign == BACKWARD ? "AWAY" : "UNCERTAIN") << ","\n''',
    '''            << snap.direction.phase << "," << snap.direction.direction_deg << ","\n            << snap.direction.apparent_to_deg << "," << snap.direction.apparent_from_deg << ","\n            << snap.direction.sense_coherence << "," << snap.direction.uncertainty_deg << ","\n            << snap.direction.confidence << "," << snap.direction.amplitude << ","\n            << (snap.direction.sign == FORWARD ? "POSITIVE_AXIS" : snap.direction.sign == BACKWARD ? "NEGATIVE_AXIS" : "UNCERTAIN") << ","\n''')

replace_once(
    "src/util/W3dSimCommon.cpp",
    '''        std::cout << "sign: TOWARD=" << nToward << " (" << pct(nToward) << "%)"\n                  << " AWAY=" << nAway << " (" << pct(nAway) << "%)"\n''',
    '''        std::cout << "sense: +AXIS=" << nToward << " (" << pct(nToward) << "%)"\n                  << " -AXIS=" << nAway << " (" << pct(nAway) << "%)"\n''')

# Register the focused unit tests without forcing them through the simulation/LaTeX matrix.
replace_once(
    "CMakeLists.txt",
    'add_executable(freq-track "${OCEAN_IMU_TESTS_DIR}/freq/freq-track.cpp")\n',
    'add_executable(freq-track "${OCEAN_IMU_TESTS_DIR}/freq/freq-track.cpp")\n'
    'add_executable(wave-direction-test "${OCEAN_IMU_TESTS_DIR}/wave_dir/wave-direction-test.cpp")\n')

replace_once(
    "Makefile",
    '''\t$(REPO_ROOT)/tests/spectrum \\\n\t$(REPO_ROOT)/tests/wave_sim\n''',
    '''\t$(REPO_ROOT)/tests/spectrum \\\n\t$(REPO_ROOT)/tests/wave_dir \\\n\t$(REPO_ROOT)/tests/wave_sim\n''')

build_path = ROOT / ".github/workflows/build.yml"
build_text = build_path.read_text(encoding="utf-8")
anchor = '''  ou-tuning:\n'''
if build_text.count(anchor) != 1:
    raise RuntimeError("build workflow anchor not found exactly once")
wave_job = '''  wave-direction:\n    runs-on: ubuntu-latest\n    steps:\n      - name: Checkout\n        uses: actions/checkout@v6\n\n      - name: Install build dependencies\n        run: |\n          sudo apt-get update\n          sudo apt-get install -y g++ make libeigen3-dev\n\n      - name: Build and run wave direction tests\n        run: |\n          cd tests/wave_dir\n          make build\n          ./run_tests.sh\n\n'''
build_path.write_text(build_text.replace(anchor, wave_job + anchor, 1), encoding="utf-8")

print("Applied apparent wave propagation direction code and tests.")
