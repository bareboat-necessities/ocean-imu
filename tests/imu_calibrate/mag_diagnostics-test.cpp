// Copyright 2026, Mikhail Grushinskiy
// Analysis behind the AtomS3R magnetometer diagnostic sketch: field/dip
// statistics, the gyro axis-sign test and magnetometer/accelerometer alignment.
#define EIGEN_NON_ARDUINO
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
#include "AtomS3R/AtomS3R_Bmm150AuxPreset.h"
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

// BMI270 with a BMM150 behind its AUX interface, as configured by M5Unified.
struct FakeBmi270 {
  uint8_t reg[128]{};
  uint8_t bmm[128]{};
  int manual_ops_in_data_mode = 0;
  bool fail_rep_write = false;
  FakeBmi270() {
    reg[0x4B] = 0x20; reg[0x4C] = 0x4F; reg[0x4D] = 0x42; reg[0x7D] = 0x0F;
    bmm[0x40] = 0x32; bmm[0x4B] = 0x01; bmm[0x4C] = 0x38; bmm[0x51] = 0; bmm[0x52] = 0;
  }
  bool manual() const { return (reg[0x4C] & 0x80) && !(reg[0x7D] & 0x01); }
  bool readRegister(uint8_t r, uint8_t* b, size_t n) { for (size_t i = 0; i < n; ++i) b[i] = reg[r + i]; return true; }
  bool writeRegister8(uint8_t r, uint8_t v) {
    reg[r] = v;
    if (r == 0x4E) {  // AUX_WR_ADDR triggers the write of AUX_WR_DATA
      if (!manual()) ++manual_ops_in_data_mode;
      else if (!(fail_rep_write && v == 0x52)) bmm[v] = reg[0x4F];
    }
    if (r == 0x4D && manual()) reg[0x04] = bmm[v];  // AUX_RD_ADDR triggers a read
    return true;
  }
};
static void noWait(uint32_t) {}

static void testAuxPreset() {
  using atoms3r_ical::Bmm150AuxPreset;
  using atoms3r_ical::Bmm150RegState;
  FakeBmi270 dev;
  Bmm150AuxPreset<FakeBmi270> preset(&dev, noWait);
  Bmm150RegState q;
  require(preset.query(q) && q.chip_id == 0x32 && q.nXY() == 1 && q.nZ() == 1 && q.mode == 0x38,
          "query reads the M5Unified default configuration");
  require(dev.reg[0x4C] == 0x4F && dev.reg[0x4D] == 0x42 && dev.reg[0x7D] == 0x0F, "query restores data mode");

  Bmm150RegState before, after;
  require(preset.apply(23, 40, before, after), "low-noise preset applies");
  require(after.nXY() == 47 && after.nZ() == 41 && after.mode == 0x38, "preset reads back at 30 Hz normal mode");
  require(dev.bmm[0x51] == 23 && dev.bmm[0x52] == 40 && dev.bmm[0x4C] == 0x38, "sensor holds the preset");
  require(dev.reg[0x4C] == 0x4F && dev.reg[0x4D] == 0x42 && dev.reg[0x7D] == 0x0F, "apply restores data mode");
  require(dev.manual_ops_in_data_mode == 0, "AUX writes happen only in manual mode");

  FakeBmi270 bad;
  bad.fail_rep_write = true;
  Bmm150AuxPreset<FakeBmi270> failing(&bad, noWait);
  require(!failing.apply(23, 40, before, after), "a failed readback is reported");
  require(bad.bmm[0x51] == 0 && bad.bmm[0x4C] == 0x38, "a failed apply restores the old repetitions");
  require(bad.reg[0x4C] == 0x4F && bad.reg[0x7D] == 0x0F, "a failed apply still restores data mode");

  FakeBmi270 other;
  other.bmm[0x40] = 0x00;
  Bmm150AuxPreset<FakeBmi270> wrong(&other, noWait);
  require(!wrong.apply(23, 40, before, after) && other.bmm[0x51] == 0, "no writes without a BMM150 chip id");
}

int main() {
  testAuxPreset();
  testFieldStats();
  testAxisConsistency();
  testAlignment();
  testHeadingError();
  if (failures) return 1;
  std::puts("mag_diagnostics-test: OK");
  return 0;
}
