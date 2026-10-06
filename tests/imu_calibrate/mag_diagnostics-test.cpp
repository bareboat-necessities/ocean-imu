// Copyright 2026, Mikhail Grushinskiy
// Analysis behind the AtomS3R magnetometer diagnostic sketch: field/dip
// statistics, the gyro axis-sign test and magnetometer/accelerometer alignment.
#define EIGEN_NON_ARDUINO
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
#include <cmath>
#include <string>
#include <cstdio>
#include <random>

using namespace atoms3r_magdiag;
static int failures = 0;
static void require(bool ok, const char* message) {
  if (!ok) { std::fprintf(stderr, "FAIL: %s\n", message); ++failures; }
}

static constexpr int NP = 300, NR = 600;
static const V3 kFieldNed(20.3f, 0.0f, 46.8f);  // dip 66.5 deg
static const V3 kHard(12.0f, -7.0f, 30.0f);

static M3 softIron() {
  M3 W;
  W << 1.1f, 0.05f, 0.0f, 0.05f, 0.9f, 0.02f, 0.0f, 0.02f, 1.6f;
  return W;
}

static M3 randomAttitude(std::mt19937& rng) {
  std::uniform_real_distribution<float> u(-1.0f, 1.0f);
  const V3 axis = V3(u(rng), u(rng), u(rng)).normalized();
  return Eigen::AngleAxisf(3.1f * u(rng), axis).toRotationMatrix();
}

// Raw reading for a body-frame field, with the calibration MagModel{W, kHard} inverting it.
static V3 rawFor(const V3& body) { return softIron().inverse() * body + kHard; }

static void makePoses(Pose* p, const M3& mag_frame, std::mt19937& rng) {
  std::normal_distribution<float> n(0.0f, 1.0f);
  for (int i = 0; i < NP; ++i) {
    const M3 Rbw = randomAttitude(rng);
    p[i].m_raw = rawFor(mag_frame * (Rbw.transpose() * kFieldNed)) + 0.3f * V3(n(rng), n(rng), n(rng));
    p[i].a = Rbw.transpose() * V3(0, 0, -9.8f) + 0.05f * V3(n(rng), n(rng), n(rng));
  }
}

static void makePairs(GyroMagPair* p, const M3& mapping_error, std::mt19937& rng) {
  std::uniform_real_distribution<float> u(-1.0f, 1.0f);
  for (int i = 0; i < NR; ++i) {
    const M3 Rbw = randomAttitude(rng);
    const V3 theta = 0.3f * V3(u(rng), u(rng), u(rng));
    const M3 step = Eigen::AngleAxisf(theta.norm(), theta.normalized()).toRotationMatrix();
    p[i].m0_raw = rawFor(mapping_error * (Rbw.transpose() * kFieldNed));
    p[i].m1_raw = rawFor(mapping_error * ((Rbw * step).transpose() * kFieldNed));
    p[i].theta = theta;
  }
}

static MagModel trueModel() {
  MagModel m;
  m.A = softIron();
  m.b = kHard;
  return m;
}

static void testFieldStats() {
  std::mt19937 rng(1);
  static Pose poses[NP];
  makePoses(poses, M3::Identity(), rng);
  const FieldStats good = fieldStats(poses, NP, trueModel());
  require(good.n == NP, "all poses counted");
  require(std::fabs(good.norm_mean - kFieldNed.norm()) < 0.5f, "true calibration recovers field strength");
  require(good.spread() < 0.08f, "true calibration keeps |m| nearly constant");
  require(std::fabs(good.dip_mean - 66.5f) < 1.0f && good.dip_sd < 1.0f, "true calibration recovers dip");

  MagModel stale = trueModel();
  stale.A = M3::Identity();  // the soft-iron part is missing
  const FieldStats bad = fieldStats(poses, NP, stale);
  require(bad.spread() > 0.3f && bad.dip_sd > 3.0f, "a stale calibration shows |m| spread and dip scatter");
}

static void testAxisConsistency() {
  std::mt19937 rng(2);
  static Pose poses[NP];
  static GyroMagPair pairs[NR];
  char text[16];

  makePoses(poses, M3::Identity(), rng);
  makePairs(pairs, M3::Identity(), rng);
  AxisResult r = axisConsistency(pairs, NR, poses, NP, trueModel());
  require(r.pairs > 500, "rotating pairs are used");
  require(r.best == 0 && r.identity_rms < 1e-3f && r.second_rms > 0.1f, "consistent axes pick the current mapping");

  // A flipped magnetometer z: rotations alone cannot tell T from -T; dip can.
  const M3 flip_z = V3(1, 1, -1).asDiagonal();
  makePoses(poses, flip_z, rng);
  makePairs(pairs, flip_z, rng);
  r = axisConsistency(pairs, NR, poses, NP, trueModel());
  permutationText(r.best, text);
  require(signedPermutation(r.best).isApprox(flip_z), "a flipped z axis is identified");
  require(std::string(text) == "+x,+y,-z", "permutation text names the flipped axis");
  require(r.best_dip > 60.0f && r.identity_rms > 0.2f, "the winner restores a downward field");

  // Swapped x/y.
  M3 swap_xy = M3::Zero();
  swap_xy(0, 1) = swap_xy(1, 0) = swap_xy(2, 2) = 1.0f;
  makePoses(poses, swap_xy, rng);
  makePairs(pairs, swap_xy, rng);
  r = axisConsistency(pairs, NR, poses, NP, trueModel());
  require(signedPermutation(r.best).isApprox(swap_xy), "swapped x/y axes are identified");
}

static void testAlignment() {
  std::mt19937 rng(3);
  static Pose poses[NP];
  const M3 tilt = Eigen::AngleAxisf(4.0f / kRadToDeg, V3(0.3f, 1.0f, 0.2f).normalized()).toRotationMatrix();
  makePoses(poses, tilt, rng);
  const Alignment a = solveAlignment(poses, NP, trueModel());
  require(a.ok, "alignment solved");
  require(a.R.isApprox(tilt.transpose(), 5e-3f), "alignment recovers the magnetometer frame rotation");
  require(a.dip_sd_before > 1.5f && a.dip_sd_after < 0.7f, "alignment makes the dip constant");
  require(std::fabs(a.dip_deg - 66.5f) < 0.5f, "aligned dip matches the field");

  makePoses(poses, M3::Identity(), rng);
  const Alignment none = solveAlignment(poses, NP, trueModel());
  require(none.ok && none.rotvec_deg.norm() < 0.3f, "aligned sensors report no rotation");
}

static void testHeadingError() {
  require(headingErrorDeg(M3::Identity(), 66.5f, 20.0f) < 1e-3f, "no misalignment, no heading error");
  const M3 pitch4 = Eigen::AngleAxisf(4.0f / kRadToDeg, V3::UnitY()).toRotationMatrix();
  const float level = headingErrorDeg(pitch4, 66.5f, 0.0f);
  require(level > 8.5f && level < 9.7f, "4 deg pitch misalignment gives about tan(dip)*4 deg at E/W");
  const float flat = headingDeg(V3(0, 0, -9.8f), V3(0, -20.3f, 46.8f));
  require(std::fabs(flat - 90.0f) < 1e-3f, "level board with field to its left heads east");
}

int main() {
  testFieldStats();
  testAxisConsistency();
  testAlignment();
  testHeadingError();
  if (failures) return 1;
  std::puts("mag_diagnostics-test: OK");
  return 0;
}
