#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Analysis used by the AtomS3R magnetometer diagnostic sketch (host-testable,
  no Arduino dependency). Frames follow AtomS3R_ImuUnits.h: body NED, accel is
  specific force (a still sensor reads -g along the down axis), gyro in rad/s.

  - fieldStats: calibrated field magnitude spread and dip angle against the
    measured gravity direction. A calibration that fits the sensor where it is
    used keeps |m| nearly constant and the dip nearly constant.
  - axisConsistency: compares magnetometer changes with the rotation the gyro
    measured between two readings, for every signed axis permutation of the
    magnetometer. The identity hypothesis should win clearly.
  - solveAlignment: the small rotation between magnetometer and accelerometer
    frames that makes the dip angle constant. An ellipsoid fit cannot see it.
  - headingErrorDeg: tilt-compensated heading error that a given magnetometer
    frame rotation and dip produce at chosen headings and tilts.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <cstdint>

namespace atoms3r_magdiag {

using V3 = Eigen::Matrix<float, 3, 1>;
using M3 = Eigen::Matrix<float, 3, 3>;

static constexpr float kRadToDeg = 57.29577951308232f;

// Calibration as applied at runtime: m_cal = A * (m_raw - b).
struct MagModel {
  M3 A = M3::Identity();
  V3 b = V3::Zero();
  V3 apply(const V3& m_raw) const { return A * (m_raw - b); }
};

// One still pose: averaged raw magnetometer and averaged specific force.
struct Pose {
  V3 m_raw = V3::Zero();
  V3 a = V3::Zero();
};

// Two consecutive distinct magnetometer readings and the gyro rotation vector
// integrated between them (body frame, rad).
struct GyroMagPair {
  V3 m0_raw = V3::Zero();
  V3 m1_raw = V3::Zero();
  V3 theta = V3::Zero();
};

static inline float dipDeg(const V3& m, const V3& a) {
  const float mn = m.norm(), an = a.norm();
  if (!(mn > 1e-6f) || !(an > 1e-6f)) return NAN;
  const V3 down = -a / an;
  float s = m.dot(down) / mn;
  s = s > 1.0f ? 1.0f : (s < -1.0f ? -1.0f : s);
  return std::asin(s) * kRadToDeg;
}

struct FieldStats {
  int n = 0;
  float norm_mean = 0, norm_sd = 0, norm_min = 0, norm_max = 0;
  float dip_mean = 0, dip_sd = 0;
  // (max - min) / mean; a good calibration keeps this within a few percent.
  float spread() const { return norm_mean > 0 ? (norm_max - norm_min) / norm_mean : NAN; }
};

static inline FieldStats fieldStats(const Pose* p, int n, const MagModel& cal) {
  FieldStats s;
  double sn = 0, snn = 0, sd = 0, sdd = 0;
  for (int i = 0; i < n; ++i) {
    const V3 m = cal.apply(p[i].m_raw);
    const float mn = m.norm(), dip = dipDeg(m, p[i].a);
    if (!std::isfinite(mn) || !std::isfinite(dip)) continue;
    if (s.n == 0) { s.norm_min = s.norm_max = mn; }
    s.norm_min = mn < s.norm_min ? mn : s.norm_min;
    s.norm_max = mn > s.norm_max ? mn : s.norm_max;
    sn += mn; snn += double(mn) * mn; sd += dip; sdd += double(dip) * dip;
    ++s.n;
  }
  if (s.n > 0) {
    s.norm_mean = float(sn / s.n);
    s.dip_mean = float(sd / s.n);
    s.norm_sd = float(std::sqrt(std::fmax(0.0, snn / s.n - (sn / s.n) * (sn / s.n))));
    s.dip_sd = float(std::sqrt(std::fmax(0.0, sdd / s.n - (sd / s.n) * (sd / s.n))));
  }
  return s;
}

// Signed axis permutation number k in [0, 48): 6 orders x 8 sign patterns.
static inline M3 signedPermutation(int k) {
  static const int order[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  const int o = k / 8, sg = k % 8;
  M3 T = M3::Zero();
  for (int r = 0; r < 3; ++r) T(r, order[o][r]) = (sg >> r) & 1 ? -1.0f : 1.0f;
  return T;
}

// Human-readable form: output axis r takes "+x", "-y", ... of the current mapping.
static inline void permutationText(int k, char out[16]) {
  const M3 T = signedPermutation(k);
  const char* names = "xyz";
  int w = 0;
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      if (T(r, c) != 0) {
        out[w++] = T(r, c) > 0 ? '+' : '-';
        out[w++] = names[c];
        if (r < 2) out[w++] = ',';
      }
  out[w] = 0;
}

struct AxisResult {
  int pairs = 0;          // pairs with enough rotation to test
  int best = 0;           // best permutation index (0 = current mapping)
  float best_rms = NAN;   // relative residual of the best hypothesis
  float identity_rms = NAN;
  float second_rms = NAN; // best hypothesis outside the winner's +/- pair
  float best_dip = NAN;   // mean dip with the winning hypothesis
};

// Relative residual |T m1 - Exp(-theta) T m0| / |m| for each hypothesis T.
// T and -T fit rotations equally well; the dip sign separates them, since
// the field points into the ground (positive dip) in the northern hemisphere.
static inline AxisResult axisConsistency(const GyroMagPair* p, int n, const Pose* poses, int n_poses,
                                         const MagModel& cal, bool northern_hemisphere = true,
                                         float min_rotation_rad = 0.05f) {
  AxisResult r;
  double sse[48] = {};
  for (int i = 0; i < n; ++i) {
    const float ang = p[i].theta.norm();
    if (ang < min_rotation_rad || ang > 1.0f) continue;
    const V3 c0 = cal.apply(p[i].m0_raw), c1 = cal.apply(p[i].m1_raw);
    const float scale = 0.5f * (c0.norm() + c1.norm());
    if (!(scale > 1e-3f)) continue;
    const M3 Rinv = Eigen::AngleAxisf(-ang, p[i].theta / ang).toRotationMatrix();
    for (int k = 0; k < 48; ++k) {
      const M3 T = signedPermutation(k);
      const V3 e = (T * c1 - Rinv * (T * c0)) / scale;
      sse[k] += e.squaredNorm();
    }
    ++r.pairs;
  }
  if (r.pairs == 0) return r;
  auto meanDip = [&](int k) {
    const M3 T = signedPermutation(k);
    double s = 0; int m = 0;
    for (int i = 0; i < n_poses; ++i) {
      const float d = dipDeg(T * cal.apply(poses[i].m_raw), poses[i].a);
      if (std::isfinite(d)) { s += d; ++m; }
    }
    return m ? float(s / m) : float(NAN);
  };
  int best = 0;
  for (int k = 1; k < 48; ++k) if (sse[k] < sse[best]) best = k;
  // Index k ^ 7 is -T(k): same axis order, every sign flipped.
  const float dip = meanDip(best);
  if (std::isfinite(dip) && ((dip < 0) == northern_hemisphere)) best ^= 7;
  int second = -1;
  for (int k = 0; k < 48; ++k)
    if (k != best && k != (best ^ 7) && (second < 0 || sse[k] < sse[second])) second = k;
  r.best = best;
  r.best_dip = meanDip(best);
  r.best_rms = float(std::sqrt(sse[best] / r.pairs));
  r.second_rms = float(std::sqrt(sse[second] / r.pairs));
  r.identity_rms = float(std::sqrt(sse[0] / r.pairs));
  return r;
}

struct Alignment {
  bool ok = false;
  M3 R = M3::Identity();  // corrected m = R * m_cal
  V3 rotvec_deg = V3::Zero();
  float dip_deg = NAN;
  float dip_sd_before = NAN, dip_sd_after = NAN;
};

// Least squares on d_i . (Exp(delta) R m_i) = sin(dip), linearised and iterated.
static inline Alignment solveAlignment(const Pose* p, int n, const MagModel& cal) {
  Alignment out;
  if (n < 12) return out;
  auto dipSd = [&](const M3& R, float& mean) {
    double s = 0, ss = 0; int k = 0;
    for (int i = 0; i < n; ++i) {
      const float d = dipDeg(R * cal.apply(p[i].m_raw), p[i].a);
      if (!std::isfinite(d)) continue;
      s += d; ss += double(d) * d; ++k;
    }
    if (k == 0) { mean = NAN; return float(NAN); }
    mean = float(s / k);
    return float(std::sqrt(std::fmax(0.0, ss / k - (s / k) * (s / k))));
  };
  float mean0 = NAN;
  out.dip_sd_before = dipSd(M3::Identity(), mean0);
  M3 R = M3::Identity();
  for (int it = 0; it < 6; ++it) {
    Eigen::Matrix4d H = Eigen::Matrix4d::Zero();
    Eigen::Vector4d g = Eigen::Vector4d::Zero();
    int used = 0;
    for (int i = 0; i < n; ++i) {
      const V3 m = R * cal.apply(p[i].m_raw);
      const float mn = m.norm(), an = p[i].a.norm();
      if (!(mn > 1e-6f) || !(an > 1e-6f)) continue;
      const V3 u = m / mn, d = -p[i].a / an;
      Eigen::Vector4d j;
      const V3 c = u.cross(d);
      j << c.x(), c.y(), c.z(), -1.0;
      H += j * j.transpose();
      g += j * double(d.dot(u));
      ++used;
    }
    if (used < 12) return out;
    // Small ridge keeps a poorly excited axis at zero instead of inventing it.
    for (int k = 0; k < 3; ++k) H(k, k) += 1e-3 * used;
    // Residual r = d.u + j.x with x = [delta; s]; normal equations H x = -g.
    const Eigen::Vector4d x = H.ldlt().solve(-g);
    if (!x.allFinite()) return out;
    const V3 delta(float(x(0)), float(x(1)), float(x(2)));
    const float ang = delta.norm();
    if (ang > 1e-9f) R = Eigen::AngleAxisf(ang, delta / ang).toRotationMatrix() * R;
    if (ang < 1e-6f) break;
  }
  const Eigen::AngleAxisf aa(R);
  out.R = R;
  out.rotvec_deg = aa.axis() * (aa.angle() * kRadToDeg);
  out.dip_sd_after = dipSd(R, out.dip_deg);
  out.ok = std::isfinite(out.dip_sd_after);
  return out;
}

// Tilt-compensated heading (deg, 0..360) from body-frame specific force and field.
static inline float headingDeg(const V3& a, const V3& m) {
  const V3 d = -a.normalized();
  V3 e = d.cross(m);
  const float en = e.norm();
  if (!(en > 1e-9f)) return NAN;
  e /= en;
  const V3 nn = e.cross(d);
  float h = std::atan2(e.x(), nn.x()) * kRadToDeg;
  return h < 0 ? h + 360.0f : h;
}

static inline float wrap180(float d) {
  while (d >= 180.0f) d -= 360.0f;
  while (d < -180.0f) d += 360.0f;
  return d;
}

// Largest heading error over 8 headings for a magnetometer frame that is
// rotated by Rerr relative to the accelerometer, at a given tilt (deg) about
// body x (roll) and body y (pitch). Field dip in degrees.
static inline float headingErrorDeg(const M3& Rerr, float dip_deg, float tilt_deg) {
  const float I = dip_deg / kRadToDeg, t = tilt_deg / kRadToDeg;
  const V3 mw(std::cos(I), 0.0f, std::sin(I));
  float worst = 0;
  for (int h = 0; h < 8; ++h) {
    const float psi = h * 45.0f / kRadToDeg;
    for (int axis = 0; axis < 2; ++axis) {
      for (int sgn = -1; sgn <= 1; sgn += 2) {
        if (tilt_deg == 0.0f && (axis == 1 || sgn == 1)) continue;
        const M3 tilt = Eigen::AngleAxisf(sgn * t, axis == 0 ? V3::UnitX() : V3::UnitY()).toRotationMatrix();
        const M3 Rbw = Eigen::AngleAxisf(psi, V3::UnitZ()).toRotationMatrix() * tilt;
        const V3 a = Rbw.transpose() * V3(0, 0, -1);
        const V3 m = Rerr * (Rbw.transpose() * mw);
        const float truth = headingDeg(a, Rbw.transpose() * mw);
        const float e = std::fabs(wrap180(headingDeg(a, m) - truth));
        worst = e > worst ? e : worst;
      }
    }
  }
  return worst;
}

// Per-axis running mean/variance.
struct Stats3 {
  int n = 0;
  Eigen::Vector3d mean = Eigen::Vector3d::Zero(), m2 = Eigen::Vector3d::Zero();
  void add(const V3& v) {
    ++n;
    const Eigen::Vector3d x = v.cast<double>(), d = x - mean;
    mean += d / n;
    m2 += d.cwiseProduct(x - mean);
  }
  V3 sd() const { return n > 1 ? V3((m2 / (n - 1)).cwiseSqrt().cast<float>()) : V3(V3::Zero()); }
  V3 avg() const { return V3(mean.cast<float>()); }
};

}  // namespace atoms3r_magdiag
