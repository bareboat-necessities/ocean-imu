// Copyright 2026, Mikhail Grushinskiy
// Analysis behind the AtomS3R magnetometer diagnostic sketch: field/dip
// statistics, the gyro axis-sign test and magnetometer/accelerometer alignment.
#define EIGEN_NON_ARDUINO
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
#include "AtomS3R/AtomS3R_Bmm150AuxPreset.h"
#include "AtomS3R/AtomS3R_Bmm150Compensation.h"
#include <cmath>
#include <string>
#include <cstdio>
#include <random>
#include <algorithm>

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
  int writes = 0, fail_write_at = -1;  // fail_write_at: index of one bus write that fails
  FakeBmi270() {
    reg[0x4B] = 0x20; reg[0x4C] = 0x4F; reg[0x4D] = 0x42; reg[0x7D] = 0x0F;
    bmm[0x40] = 0x32; bmm[0x4B] = 0x01; bmm[0x4C] = 0x38; bmm[0x51] = 0; bmm[0x52] = 0;
  }
  bool manual() const { return (reg[0x4C] & 0x80) && !(reg[0x7D] & 0x01); }
  bool readRegister(uint8_t r, uint8_t* b, size_t n) { for (size_t i = 0; i < n; ++i) b[i] = reg[r + i]; return true; }
  bool writeRegister8(uint8_t r, uint8_t v) {
    if (writes++ == fail_write_at) return false;
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

  // One failed bus write anywhere (entering manual mode, the AUX access, or a
  // restore write) never leaves AUX polling off or the interface in manual mode.
  for (int op = 0; op < 3; ++op) {
    FakeBmi270 probe;
    Bmm150AuxPreset<FakeBmi270> counted(&probe, noWait);
    uint8_t trim[4];
    const bool clean = op == 0 ? counted.query(q) : op == 1 ? counted.apply(23, 40, before, after)
                                                            : counted.readBlock(0x5D, 4, trim);
    require(clean, "the operation succeeds without bus failures");
    bool all_restored = true;
    for (int k = 0; k < probe.writes; ++k) {
      FakeBmi270 d;
      d.fail_write_at = k;
      Bmm150AuxPreset<FakeBmi270> p(&d, noWait);
      if (op == 0) (void)p.query(q);
      else if (op == 1) (void)p.apply(23, 40, before, after);
      else (void)p.readBlock(0x5D, 4, trim);
      all_restored = all_restored && d.reg[0x4C] == 0x4F && d.reg[0x4D] == 0x42 && d.reg[0x7D] == 0x0F;
    }
    require(all_restored, "data mode is restored after any single failed write");
  }
}

// Bosch BMM150 Sensor API integer compensation (bmm150.c), used as an
// independent reference for the floating-point implementation.
static int boschX(int16_t x, uint16_t rhall, int8_t dx1, int8_t dx2, const atoms3r_ical::Bmm150Trim& t) {
  const uint16_t x0 = rhall ? rhall : t.dig_xyz1;
  const int32_t x1 = int32_t(t.dig_xyz1) * 16384;
  const uint16_t x2 = uint16_t(uint16_t(x1 / x0) - uint16_t(0x4000));
  int16_t r = int16_t(x2);
  const int32_t x3 = int32_t(r) * int32_t(r);
  const int32_t x4 = int32_t(t.dig_xy2) * (x3 / 128);
  const int32_t x5 = int32_t(int16_t(t.dig_xy1) * 128);
  const int32_t x6 = int32_t(r) * x5;
  const int32_t x7 = ((x4 + x6) / 512) + int32_t(0x100000);
  const int32_t x8 = int32_t(int16_t(dx2) + int16_t(0xA0));
  const int32_t x9 = (x7 * x8) / 4096;
  const int32_t x10 = int32_t(x) * x9;
  r = int16_t(x10 / 8192);
  return (r + int16_t(dx1) * 8) / 16;
}
static int boschZ(int16_t z, uint16_t rhall, const atoms3r_ical::Bmm150Trim& t) {
  const int32_t z0 = int32_t(z) - int32_t(t.dig_z4);
  const int32_t z1 = int32_t(rhall) - int32_t(t.dig_xyz1);
  const int32_t z2 = int32_t(t.dig_z3) * z1;
  const int32_t z3 = int32_t(t.dig_z1) * (int32_t(rhall) * 2);
  const int32_t z4 = (z3 + (1 << 15)) / 65536;
  const int32_t z5 = (z0 * 131072) - z2;
  const int32_t z6 = z5 / ((z4 + int32_t(t.dig_z2)) * 4);
  return int(z6 / 16);
}

static void testCompensation() {
  using atoms3r_ical::Bmm150Trim;
  using atoms3r_ical::Bmm150Raw;
  // Representative factory trim block (registers 0x5D..0x71).
  uint8_t regs[Bmm150Trim::COUNT]{};
  auto set = [&](uint8_t reg, uint8_t v) { regs[reg - Bmm150Trim::FIRST_REG] = v; };
  auto set16 = [&](uint8_t reg, int v) { set(reg, uint8_t(v & 0xFF)); set(reg + 1, uint8_t((v >> 8) & 0xFF)); };
  set(0x5D, uint8_t(int8_t(-2))); set(0x5E, 3); set16(0x62, -12); set(0x64, 26); set(0x65, 24);
  set16(0x68, 749); set16(0x6A, 24747); set16(0x6C, 0x8000 | 6800); set16(0x6E, -60);
  set(0x70, uint8_t(int8_t(-3))); set(0x71, 29);
  const Bmm150Trim t = Bmm150Trim::fromRegisters(regs);
  require(t.valid && t.dig_x1 == -2 && t.dig_y1 == 3 && t.dig_z4 == -12 && t.dig_x2 == 26 && t.dig_y2 == 24 &&
          t.dig_z2 == 749 && t.dig_z1 == 24747 && t.dig_xyz1 == 6800 && t.dig_z3 == -60 && t.dig_xy2 == -3 &&
          t.dig_xy1 == 29, "trim registers parse (xyz1 MSB bit 7 masked)");

  // Data bytes as the BMI270 AUX mirror holds them.
  auto bytes = [](int x, int y, int z, int rhall, uint8_t d[8]) {
    const uint16_t bx = uint16_t(x << 3), by = uint16_t(y << 3), bz = uint16_t(z << 1), br = uint16_t(rhall << 2) | 1u;
    d[0] = bx & 0xFF; d[1] = bx >> 8; d[2] = by & 0xFF; d[3] = by >> 8;
    d[4] = bz & 0xFF; d[5] = bz >> 8; d[6] = br & 0xFF; d[7] = br >> 8;
  };
  uint8_t d[8];
  bytes(-1234, 987, -3456, 6789, d);
  const Bmm150Raw r = Bmm150Raw::fromData(d);
  require(r.x == -1234 && r.y == 987 && r.z == -3456 && r.rhall == 6789, "data registers parse with sign");

  std::mt19937 rng(7);
  std::uniform_int_distribution<int> xy(-1500, 1500), zz(-6000, 6000), rh(5500, 8200);
  double worst = 0, mean_err = 0, signed_err = 0;
  for (int i = 0; i < 2000; ++i) {
    bytes(xy(rng), xy(rng), zz(rng), rh(rng), d);
    const Bmm150Raw q = Bmm150Raw::fromData(d);
    float f[3]{};
    require(atoms3r_ical::bmm150Compensate(q, t, f), "valid samples compensate");
    const double ex = std::fabs(f[0] - boschX(q.x, q.rhall, t.dig_x1, t.dig_x2, t));
    const double ey = std::fabs(f[1] - boschX(q.y, q.rhall, t.dig_y1, t.dig_y2, t));
    const double ez = std::fabs(f[2] - boschZ(q.z, q.rhall, t));
    worst = std::max(worst, std::max(ex, std::max(ey, ez)));
    mean_err += (ex + ey + ez) / 3.0;
    signed_err += (f[0] - boschX(q.x, q.rhall, t.dig_x1, t.dig_x2, t)) + (f[2] - boschZ(q.z, q.rhall, t));
  }
  mean_err /= 2000; signed_err /= 4000;
  std::printf("compensation vs Bosch integer: worst %.3f mean %.3f signed %.3f uT\n", worst, mean_err, signed_err);
  // The integer reference truncates at three stages (1/8192, 1/16 and the
  // final integer uT), so differences up to about 1.2 uT are rounding; a
  // formula error would grow with the field and bias the mean.
  require(worst < 1.25 && mean_err < 0.6, "float compensation matches Bosch's integer reference within its rounding");

  // RHALL tracks temperature: the Z offset term moves with it, so the same
  // raw Z gives different compensated values at different RHALL.
  float a[3]{}, b[3]{};
  bytes(100, 100, 2000, 6500, d); atoms3r_ical::bmm150Compensate(Bmm150Raw::fromData(d), t, a);
  bytes(100, 100, 2000, 7100, d); atoms3r_ical::bmm150Compensate(Bmm150Raw::fromData(d), t, b);
  require(std::fabs(a[2] - b[2]) > 1.0f, "RHALL changes the compensated Z");

  float o[3];
  bytes(-4096, 0, 0, 6500, d);
  require(!atoms3r_ical::bmm150Compensate(Bmm150Raw::fromData(d), t, o), "X overflow is rejected");
  bytes(0, 0, -16384, 6500, d);
  require(!atoms3r_ical::bmm150Compensate(Bmm150Raw::fromData(d), t, o), "Z overflow is rejected");
  bytes(10, 10, 10, 0, d);
  require(!atoms3r_ical::bmm150Compensate(Bmm150Raw::fromData(d), t, o), "missing RHALL is rejected");
  uint8_t erased[Bmm150Trim::COUNT]{};
  require(!Bmm150Trim::fromRegisters(erased).valid, "an all-zero trim block is not used");

  // Body mapping must agree in sign with M5Unified's AtomS3R path:
  // M5 (-x, y, -z) then map_sensor_xyz_to_body_ned_ (sy, sx, -sz).
  int agree = 0, total = 0;
  for (int i = 0; i < 500; ++i) {
    const int x = xy(rng), y = xy(rng), z = zz(rng);
    if (std::abs(x) < 200 || std::abs(y) < 200 || std::abs(z) < 600) continue;
    bytes(x, y, z, t.dig_xyz1, d);
    float s3[3]{}, body[3]{};
    atoms3r_ical::bmm150Compensate(Bmm150Raw::fromData(d), t, s3);
    atoms3r_ical::bmm150SensorToAtomS3RBody(s3, body);
    const float m5[3] = {float(-x), float(y), float(-z)};
    const float m5body[3] = {m5[1], m5[0], -m5[2]};
    ++total;
    agree += (body[0] > 0) == (m5body[0] > 0) && (body[1] > 0) == (m5body[1] > 0) && (body[2] > 0) == (m5body[2] > 0);
  }
  require(total > 50 && agree == total, "compensated body axes have M5Unified's AtomS3R signs");
}

int main() {
  testCompensation();
  testAuxPreset();
  testFieldStats();
  testAxisConsistency();
  testAlignment();
  testHeadingError();
  if (failures) return 1;
  std::puts("mag_diagnostics-test: OK");
  return 0;
}
