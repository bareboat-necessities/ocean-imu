#pragma once

/*
  Apparent/intrinsic wave encounter helpers.

  The IMU phase detector reports an apparent propagation sense along a measured
  axis. A moving observer measures signed encounter frequency

      Omega_e = sigma - k * d dot V_rel,

  where sigma is intrinsic water-relative angular frequency, k is wavenumber,
  d is the intrinsic propagation unit vector, and V_rel is vessel velocity
  relative to the water. If current is neglected, V_rel may be approximated by
  velocity over ground. With known current U, V_rel = V_ground - U.
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
  int intrinsic_sense = 0;  // +/- relative to the supplied axis
  Real signed_encounter_omega_rad_s =
      std::numeric_limits<Real>::quiet_NaN();
};

template <typename Real>
inline Real deep_water_wavenumber(Real intrinsic_omega_rad_s,
                                  Real gravity_ms2 = Real(9.80665)) {
  return intrinsic_omega_rad_s * intrinsic_omega_rad_s / gravity_ms2;
}

template <typename Real>
inline Real signed_encounter_omega(
    Real intrinsic_omega_rad_s,
    Real wavenumber_rad_m,
    int intrinsic_sense,
    Real vessel_speed_along_positive_axis_ms) {
  const Real sense = intrinsic_sense >= 0 ? Real(1) : Real(-1);
  return intrinsic_omega_rad_s
      - wavenumber_rad_m * sense * vessel_speed_along_positive_axis_ms;
}

template <typename Real>
inline int apparent_sense(int intrinsic_sense,
                          Real signed_encounter_omega_rad_s) {
  if (intrinsic_sense == 0 || signed_encounter_omega_rad_s == Real(0)) {
    return 0;
  }
  const int encounter_sign =
      signed_encounter_omega_rad_s > Real(0) ? 1 : -1;
  return intrinsic_sense * encounter_sign;
}

// Solve the no-current (or water-relative velocity) deep-water encounter
// equation. Multiple candidates are returned rather than silently choosing a
// branch in following seas, where encounter-frequency inversion can be
// non-unique.
//
// `scan_intervals` is retained for source compatibility with the previous
// numerical implementation. The solution is now analytic and does not depend
// on a scan resolution.
template <typename Real>
inline std::vector<DeepWaterSolution<Real>> solve_deep_water(
    Real encounter_omega_abs_rad_s,
    Real vessel_speed_along_positive_axis_ms,
    int measured_apparent_sense,
    Real gravity_ms2 = Real(9.80665),
    Real omega_min_rad_s = Real(0.05),
    Real omega_max_rad_s = Real(20),
    int scan_intervals = 4096) {
  (void)scan_intervals;
  std::vector<DeepWaterSolution<Real>> solutions;

  const auto finite = [](Real value) {
    return std::isfinite(static_cast<double>(value));
  };

  if (!(encounter_omega_abs_rad_s > Real(0)) ||
      !finite(encounter_omega_abs_rad_s) ||
      !finite(vessel_speed_along_positive_axis_ms) ||
      !(gravity_ms2 > Real(0)) || !finite(gravity_ms2) ||
      !(omega_max_rad_s > omega_min_rad_s) ||
      !finite(omega_min_rad_s) || !finite(omega_max_rad_s) ||
      (measured_apparent_sense != -1 && measured_apparent_sense != 1)) {
    return solutions;
  }

  const Real scale = std::max({
      Real(1),
      std::abs(encounter_omega_abs_rad_s),
      std::abs(omega_min_rad_s),
      std::abs(omega_max_rad_s)
  });
  const Real omega_tolerance =
      Real(64) * std::numeric_limits<Real>::epsilon() * scale;
  const Real coefficient_tolerance =
      Real(64) * std::numeric_limits<Real>::epsilon();

  auto append_candidate = [&](Real omega, int intrinsic_sense) {
    if (!finite(omega) ||
        omega < omega_min_rad_s - omega_tolerance ||
        omega > omega_max_rad_s + omega_tolerance ||
        !(omega > Real(0))) {
      return;
    }

    omega = std::clamp(omega, omega_min_rad_s, omega_max_rad_s);
    const Real k = deep_water_wavenumber(omega, gravity_ms2);
    const Real omega_e = signed_encounter_omega(
        omega, k, intrinsic_sense,
        vessel_speed_along_positive_axis_ms);
    if (!finite(k) || !finite(omega_e) ||
        apparent_sense(intrinsic_sense, omega_e) !=
            measured_apparent_sense) {
      return;
    }

    const Real residual =
        std::abs(std::abs(omega_e) - encounter_omega_abs_rad_s);
    const Real residual_scale = std::max(
        Real(1), encounter_omega_abs_rad_s);
    if (residual > Real(256) * std::numeric_limits<Real>::epsilon() *
                       residual_scale) {
      return;
    }

    for (const auto& existing : solutions) {
      if (existing.intrinsic_sense == intrinsic_sense &&
          std::abs(existing.intrinsic_omega_rad_s - omega) <=
              Real(8) * omega_tolerance) {
        return;
      }
    }

    solutions.push_back({omega, k, intrinsic_sense, omega_e});
  };

  for (int intrinsic_sense : {-1, 1}) {
    // measured_apparent_sense = intrinsic_sense * sign(Omega_e), so the
    // required signed encounter branch is q * |Omega_e| with
    // q = measured_apparent_sense * intrinsic_sense.
    const Real encounter_branch_sign =
        Real(measured_apparent_sense * intrinsic_sense);

    // Deep water: k = sigma^2/g. Therefore
    //
    //   A sigma^2 + B sigma + C = 0,
    //   A = -sense*V/g, B = 1, C = -q*|Omega_e|.
    const Real A = -Real(intrinsic_sense) *
                   vessel_speed_along_positive_axis_ms / gravity_ms2;
    const Real B = Real(1);
    const Real C = -encounter_branch_sign *
                   encounter_omega_abs_rad_s;

    if (std::abs(A) <= coefficient_tolerance) {
      append_candidate(-C / B, intrinsic_sense);
      continue;
    }

    Real discriminant = B * B - Real(4) * A * C;
    const Real discriminant_tolerance =
        Real(128) * std::numeric_limits<Real>::epsilon() *
        std::max(Real(1), std::abs(B * B) + std::abs(Real(4) * A * C));
    if (discriminant < -discriminant_tolerance) {
      continue;
    }
    if (discriminant < Real(0)) {
      discriminant = Real(0);
    }

    const Real sqrt_discriminant = std::sqrt(discriminant);

    // Stable quadratic formula. q_root is chosen to avoid subtracting nearly
    // equal values; the second root follows from the product C/A.
    const Real q_root = -Real(0.5) *
        (B + std::copysign(sqrt_discriminant, B));
    if (q_root != Real(0)) {
      append_candidate(q_root / A, intrinsic_sense);
      append_candidate(C / q_root, intrinsic_sense);
    } else {
      append_candidate(-B / (Real(2) * A), intrinsic_sense);
    }
  }

  std::sort(solutions.begin(), solutions.end(),
            [](const auto& a, const auto& b) {
              if (a.intrinsic_omega_rad_s == b.intrinsic_omega_rad_s) {
                return a.intrinsic_sense < b.intrinsic_sense;
              }
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
  const Real intrinsic_phase_speed =
      intrinsic_omega_rad_s / wavenumber_rad_m;
  return {
      intrinsic_phase_speed * direction_x + current_ms.x,
      intrinsic_phase_speed * direction_y + current_ms.y
  };
}

}  // namespace wave_encounter
