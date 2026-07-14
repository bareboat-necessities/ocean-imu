#pragma once

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
