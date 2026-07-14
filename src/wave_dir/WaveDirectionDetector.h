#pragma once

/*
  Copyright 2025-2026, Mikhail Grushinskiy

  Resolve apparent wave-propagation sense along an independently estimated
  horizontal propagation axis. The axis is undirected (modulo 180 degrees);
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

    // Normalized on/off thresholds. The lower threshold supplies hysteresis.
    Real coherence_threshold_on = Real(0.20);
    Real coherence_threshold_off = Real(0.12);

    // Absolute signal gates. These prevent a high normalized correlation made
    // from nearly zero signals from being reported as a direction.
    Real absolute_product_floor = Real(0.005);
    Real min_horizontal_rms = Real(0.01);
    Real min_vertical_slope_rms = Real(0.01);
    Real min_axis_confidence = Real(20);

    // For linear deep-water orbital motion with up-positive vertical
    // acceleration, a_parallel * d(a_vertical)/dt is negative for propagation
    // along the positive axis. Negating the coherence maps that physical
    // convention to FORWARD.
    Real convention_sign = Real(-1);
  };

  // Backward-compatible constructor. The old smoothing argument was an EMA
  // coefficient tuned for 200 Hz. Convert it to an elapsed-time constant so
  // behavior remains invariant when the sample rate changes.
  WaveDirectionDetector(Real smoothing = Real(0.002),
                        Real sensitivity = Real(0.005)) {
    config_.smoothing_time_constant_sec =
        alpha_to_tau(smoothing, Real(0.005));
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
                       Real axisConfidence =
                           std::numeric_limits<Real>::infinity()) {
    if (!(delta_t > Real(0)) || !finite(delta_t) ||
        !finite(accelX) || !finite(accelY) || !finite(accelVertical) ||
        !finite(axisX) || !finite(axisY)) {
      have_vertical_ = false;
      filtered_vertical_slope_ = Real(0);
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

    // An axial estimate is equivalent under d -> -d. Keep one representative
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

    const Real raw_vertical_slope =
        (accelVertical - previous_vertical_) / delta_t;
    previous_vertical_ = accelVertical;

    const Real slope_alpha =
        ema_alpha(delta_t, config_.vertical_slope_time_constant_sec);
    filtered_vertical_slope_ +=
        slope_alpha * (raw_vertical_slope - filtered_vertical_slope_);

    // Do not accumulate a phase decision against an axis which has not yet
    // reached the confidence required by the independent angle estimator.
    if (std::isnan(static_cast<double>(axisConfidence)) ||
        (finite(axisConfidence) &&
         axisConfidence < config_.min_axis_confidence)) {
      state_ = UNCERTAIN;
      return state_;
    }

    const Real alpha =
        ema_alpha(delta_t, config_.smoothing_time_constant_sec);
    const Real product =
        along_axis_acceleration_ * filtered_vertical_slope_;

    filtered_cross_power_ +=
        alpha * (product - filtered_cross_power_);
    filtered_horizontal_energy_ +=
        alpha * (along_axis_acceleration_ * along_axis_acceleration_ -
                 filtered_horizontal_energy_);
    filtered_vertical_slope_energy_ +=
        alpha * (filtered_vertical_slope_ * filtered_vertical_slope_ -
                 filtered_vertical_slope_energy_);

    const Real denom = std::sqrt(std::max(
        Real(0),
        filtered_horizontal_energy_ * filtered_vertical_slope_energy_));
    coherence_ = (denom > Real(1e-12))
        ? std::clamp(filtered_cross_power_ / denom, Real(-1), Real(1))
        : Real(0);

    const Real horizontal_rms =
        std::sqrt(std::max(Real(0), filtered_horizontal_energy_));
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

  // Legacy fixed-Y overload retained only for source compatibility. New code
  // must pass the independently estimated propagation axis to the overload
  // above.
  WaveDirection update(Real accelX,
                       Real accelY,
                       Real accelVertical,
                       Real delta_t) {
    return update(accelX, accelY, accelVertical,
                  Real(0), Real(1), delta_t,
                  std::numeric_limits<Real>::infinity());
  }

  WaveDirection getState() const noexcept { return state_; }
  Real getFilteredP() const noexcept { return filtered_cross_power_; }
  Real getCoherence() const noexcept { return coherence_; }
  Real getAlongAxisAcceleration() const noexcept {
    return along_axis_acceleration_;
  }
  Real getVerticalSlope() const noexcept {
    return filtered_vertical_slope_;
  }
  Real getAxisX() const noexcept { return axis_x_; }
  Real getAxisY() const noexcept { return axis_y_; }

  Real getHorizontalRms() const noexcept {
    return std::sqrt(std::max(Real(0), filtered_horizontal_energy_));
  }

  Real getVerticalSlopeRms() const noexcept {
    return std::sqrt(std::max(Real(0), filtered_vertical_slope_energy_));
  }

  Real getDirectedX() const noexcept {
    if (state_ == UNCERTAIN) {
      return std::numeric_limits<Real>::quiet_NaN();
    }
    return (state_ == FORWARD) ? axis_x_ : -axis_x_;
  }

  Real getDirectedY() const noexcept {
    if (state_ == UNCERTAIN) {
      return std::numeric_limits<Real>::quiet_NaN();
    }
    return (state_ == FORWARD) ? axis_y_ : -axis_y_;
  }

  // 0 degrees is boat +X (forward); positive angles rotate toward boat +Y
  // (starboard in the NED-compatible body convention used by the filters).
  Real getDirectedAngleDegrees() const noexcept {
    const Real x = getDirectedX();
    const Real y = getDirectedY();
    if (!finite(x) || !finite(y)) {
      return std::numeric_limits<Real>::quiet_NaN();
    }
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
    if (!(config_.smoothing_time_constant_sec > Real(0))) {
      config_.smoothing_time_constant_sec = Real(2.5);
    }
    if (!(config_.vertical_slope_time_constant_sec > Real(0))) {
      config_.vertical_slope_time_constant_sec = Real(0.05);
    }
    config_.coherence_threshold_on =
        std::clamp(config_.coherence_threshold_on, Real(0), Real(1));
    config_.coherence_threshold_off =
        std::clamp(config_.coherence_threshold_off, Real(0),
                   config_.coherence_threshold_on);
    config_.absolute_product_floor =
        std::max(Real(0), config_.absolute_product_floor);
    config_.min_horizontal_rms =
        std::max(Real(0), config_.min_horizontal_rms);
    config_.min_vertical_slope_rms =
        std::max(Real(0), config_.min_vertical_slope_rms);
    config_.min_axis_confidence =
        std::max(Real(0), config_.min_axis_confidence);
    config_.convention_sign =
        (config_.convention_sign < Real(0)) ? Real(-1) : Real(1);
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
