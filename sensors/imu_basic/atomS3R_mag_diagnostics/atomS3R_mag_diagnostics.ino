/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R magnetometer diagnostics.

  Never writes the saved calibration. Like every sketch, it starts the BMM150
  in the low-noise setting with Bosch trim/RHALL compensation
  (configureAtomS3RMag). Serial commands: 'q' reads the BMM150 configuration,
  'o' switches to the M5Unified driver default (1/1 repetitions, until
  reboot) for comparison, 'p' restores the low-noise setting, 'c' toggles
  Bosch compensation, and 'd' monitors the untouched board (field and
  temperature every 5 s). Tap the screen
  (or send 'n' over serial) to advance:

    1. Boot report    sensor wiring, BMI270 AUX settings, saved calibration
    2. STILL (10 s)   board flat and untouched: update rate, noise, dip,
                      field and gravity against the local reference
    3. POSES          calibrated accel/mag quality from still poses: set the
                      board down on each face and some tilts, 3 s each; tap
                      to finish. Noise, |a|-g and |m|-F, per-axis scale error
                      and residual bias with standard errors
    4. ROTATE         turn and tilt through every direction, as in the wizard
    5. Analysis       fresh fit with the wizard's code, saved-vs-fresh field
                      checks, axis-sign test against the gyro, mag/accel
                      frame alignment, VERDICT lines
    6. LIVE           headings with saved, fresh and aligned-fresh models for
                      the tilt test; tap to run again

  Compare runs in both settings: noise in [STILL], and the fits of two
  consecutive runs in [FIT] previous-run-vs-fresh.

  All results are printed on USB serial (115200) as [DIAG], [MAGCFG], [STILL], [POSE], [QUAL], [ROT],
  [FIT], [FIELD], [AXIS], [ALIGN], [VERDICT] and [LIVE] lines.
*/

#include <Arduino.h>
#include <algorithm>
#include <M5Unified.h>
#include <ArduinoOceanImu.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
#include "AtomS3R/AtomS3R_Bmm150AuxPreset.h"
#include "imu_calibrate/MagCalSampling.h"

// Local field used only for the checks and verdicts: WMM2025 at Fair Lawn NJ
// (40.936 N, 74.118 W), late 2026: total 50.9 uT, dip +65.7 deg. Look up yours
// (NOAA/WMM) and rebuild with -DDIAG_EXPECTED_DIP_DEG=<deg> and
// -DDIAG_EXPECTED_FIELD_UT=<uT>. Gravity is ImuCalCfg::g_cal_local.
#ifndef DIAG_EXPECTED_DIP_DEG
#define DIAG_EXPECTED_DIP_DEG 65.7f
#endif
#ifndef DIAG_EXPECTED_FIELD_UT
#define DIAG_EXPECTED_FIELD_UT 50.9f
#endif

using namespace atoms3r_ical;
using namespace atoms3r_magdiag;

static constexpr float LOOP_HZ = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = (uint32_t)(1000000.0f / LOOP_HZ);
static constexpr uint32_t STILL_MS = 10000;
static constexpr uint32_t ROTATE_MAX_MS = 180000;
static constexpr uint32_t POSE_WINDOW_MS = 200;
static constexpr int MAX_POSES = 400;
static constexpr int MAX_PAIRS = 1000;

enum class Phase : uint8_t { WAIT_STILL, STILL, WAIT_POSES, POSES, WAIT_ROTATE, ROTATE, ANALYZE, LIVE, DRIFT };

static M5Ui ui;
static ImuCalStoreNvs store;
static ImuCalBlobV5 blob{};
static RuntimeCals cals{};
static bool have_blob = false;

static imu_cal::MagCalibrator<float, 400> magCal;
static imu_cal::MagCalibration<float> fresh{};
static imu_cal::FitFail fresh_reason = imu_cal::FitFail::BAD_ARG;
static Pose poses[MAX_POSES];
static int n_poses = 0, pose_slot = 0;
static GyroMagPair pairs[MAX_PAIRS];
static int n_pairs = 0, pair_slot = 0;

static Phase phase = Phase::WAIT_STILL;
static uint32_t phase_ms = 0, next_tick_us = 0, last_sample_us = 0, last_print_ms = 0, last_draw_ms = 0;

// Latest calibrated sample and magnetometer freshness tracking.
static Vector3f a_cal = Vector3f::Zero(), w_cal = Vector3f::Zero(), m_raw = Vector3f::Zero();
static float last_temp_c = NAN;
static Vector3f last_distinct_m = Vector3f::Zero();
static bool have_m = false;
static uint32_t last_m_change_us = 0;
// Gyro-integrated attitude and the last few distinct readings with it, so each
// axis-test pair spans PAIR_SPAN readings (about 200 ms) of rotation.
static constexpr int PAIR_SPAN = 5;
static Eigen::Quaternionf q_gyro = Eigen::Quaternionf::Identity();
static Vector3f ring_m[PAIR_SPAN + 1];
static Eigen::Quaternionf ring_q[PAIR_SPAN + 1];
static int ring_n = 0, ring_head = 0;

// STILL statistics.
static Stats3 st_raw, st_cal, st_acc, st_gyr;
// Heading/dip scatter of individual readings with the saved calibration.
static double st_hdg_c = 0, st_hdg_s = 0;
static Stats3 st_dip;
static uint32_t st_intervals = 0, st_int_min = 0, st_int_max = 0;
static uint64_t st_int_sum = 0;
static uint32_t mask_ag = 0, mask_g_only = 0, mask_a_only = 0, mask_none = 0, mask_mag = 0;

// POSES: still windows of calibrated data, one per placement.
static constexpr int MAX_QPOSES = 24;
static constexpr uint32_t QPOSE_MS = 3000;
static constexpr float QPOSE_STILL_RAD_S = 0.05f, QPOSE_MOVED_RAD_S = 0.3f;
struct QPose { V3 a, a_sd, m, m_sd; };
static QPose qposes[MAX_QPOSES];
static int n_qposes = 0;
static Stats3 qw_acc, qw_mag;
static uint32_t qw_start_ms = 0;
static bool qw_moved = true;

// ROTATE: wizard capture plus our own still-pose windows.
static imu_cal::MagCaptureCfg cap_cfg;
static imu_cal::MagCapture<float, 400>* capture = nullptr;
static imu_cal::MagSampleWindow window;
static imu_cal::MagCaptureStatus cap_status = imu_cal::MagCaptureStatus::CAPTURING;
static Stats3 pw_mag, pw_acc;
static float pw_max_gyro = 0;
static uint32_t pw_start_ms = 0;

// Analysis results used by LIVE.
static MagModel saved_model, fresh_model;
static bool saved_ok = false, fresh_ok = false;
static Alignment fresh_align;
// Previous run's fresh fit: a second run at the same place measures repeatability.
static MagModel prev_model;
static bool prev_ok = false;
static float prev_temp_c = NAN;
static volatile bool analysis_done = false;

// DRIFT monitor: untouched board, 5 s averages of raw field and temperature.
static constexpr uint32_t DRIFT_PERIOD_MS = 5000;
static Stats3 dr_raw, dr_acc;
static double dr_temp_sum = 0;
static int dr_temp_n = 0;
static V3 dr_first_raw = V3::Zero();
static float dr_first_temp = NAN;
static bool dr_have_first = false;
static Phase dr_return = Phase::WAIT_STILL;
static uint32_t dr_start_ms = 0;

// Magnetometer repetition setting in use, for labelling results.
static char mag_mode[48] = "driver default (1/1 reps)";

// Repetition setting plus data source, for labelling results.
static const char* magLabel() {
  static char label[96];
  snprintf(label, sizeof(label), "%s, %s", mag_mode,
           atoms3rMagSource().compensated ? "Bosch-compensated" : "uncompensated");
  return label;
}

static void waitMs(uint32_t ms) { delay(ms); }

static void printMagState(const char* tag, const Bmm150RegState& r) {
  Serial.printf("[MAGCFG] %s: chip=0x%02X power=0x%02X mode=0x%02X rep_xy=%u (n=%u) rep_z=%u (n=%u)\n", tag,
                r.chip_id, r.power, r.mode, r.rep_xy, r.nXY(), r.rep_z, r.nZ());
}

// 'q' query, 'p' low-noise, 'o' driver default.
static void magCommand(char cmd) {
  auto* imu0 = M5.Imu.getImuInstancePtr(0);
  if (M5.Imu.getType() != m5::imu_bmi270 || !imu0) {
    Serial.println("[MAGCFG] not a BMI270 + AUX BMM150 board; nothing done");
    return;
  }
  Bmm150AuxPreset<m5::IMU_Base> aux(imu0, waitMs);
  Bmm150RegState before, after;
  if (cmd == 'q') {
    if (aux.query(before)) printMagState("current", before);
    else Serial.printf("[MAGCFG] query failed: %s\n", aux.failure());
    return;
  }
  const bool low = cmd == 'p';
  const bool ok = low ? aux.apply(Bmm150AuxPreset<m5::IMU_Base>::LOW_NOISE_REP_XY,
                                  Bmm150AuxPreset<m5::IMU_Base>::LOW_NOISE_REP_Z, before, after)
                      : aux.apply(0, 0, before, after);
  printMagState("before", before);
  if (ok) {
    printMagState("after", after);
    snprintf(mag_mode, sizeof(mag_mode), "%s (%u/%u reps, 30 Hz)", low ? "low-noise" : "driver default",
             after.nXY(), after.nZ());
    Serial.printf("[MAGCFG] %s setting active (driver default lasts until reboot); rerun STILL and ROTATE to compare\n",
                  low ? "low-noise" : "driver default");
  } else {
    Serial.printf("[MAGCFG] change FAILED: %s (previous setting kept)\n", aux.failure());
  }
}

// 'c': switch between Bosch-compensated and uncompensated data (re-applies the
// low-noise setting). The saved calibration applies only to data from the
// source it was fitted on; the choice here is kept, not switched back to it.
static void toggleCompensation() {
  const bool on = !atoms3rMagSource().compensated;
  if (configureAtomS3RMag(Serial, on)) snprintf(mag_mode, sizeof(mag_mode), "low-noise (47/41 reps, 30 Hz)");
  magSourceFollower() = nullptr;
  cals.rebuildFromBlob(blob);
  saved_ok = cals.mag.ok;
  saved_model.A = cals.mag.A;
  saved_model.b = cals.mag.b;
  have_m = false;
  Serial.printf("[MAGCFG] magnetometer now: %s; saved mag calibration %s\n", magLabel(),
                saved_ok ? "applies" : (cals.mag_source_mismatch ? "NOT applied (fitted on the other source)" : "not available"));
}

static void printMat(const char* tag, const Matrix3f& A) {
  Serial.printf("%s [%.5f %.5f %.5f; %.5f %.5f %.5f; %.5f %.5f %.5f]\n", tag,
                A(0,0), A(0,1), A(0,2), A(1,0), A(1,1), A(1,2), A(2,0), A(2,1), A(2,2));
}

static void bootReport() {
  Serial.println();
  Serial.println("[DIAG] ===== AtomS3R magnetometer diagnostics =====");
  Serial.printf("[DIAG] imu_type=%d (bmi270=%d) expected_dip=%.1f deg\n",
                (int)M5.Imu.getType(), (int)m5::imu_bmi270, (double)DIAG_EXPECTED_DIP_DEG);

  auto* imu0 = M5.Imu.getImuInstancePtr(0);
  auto* mag1 = M5.Imu.getImuInstancePtr(1);
  if (imu0) {
    // Read-only look at how the BMI270 polls the BMM150 on its AUX bus.
    uint8_t aux_conf = 0, aux_dev = 0, aux_if = 0, aux_rd = 0, pwr = 0;
    imu0->readRegister(0x44, &aux_conf, 1);
    imu0->readRegister(0x4B, &aux_dev, 1);
    imu0->readRegister(0x4C, &aux_if, 1);
    imu0->readRegister(0x4D, &aux_rd, 1);
    imu0->readRegister(0x7D, &pwr, 1);
    const int odr = aux_conf & 0x0F;
    const float aux_hz = (odr >= 1 && odr <= 11) ? 100.0f / powf(2.0f, (float)(8 - odr)) : NAN;
    Serial.printf("[DIAG] BMI270 addr=0x%02X AUX_CONF=0x%02X (aux poll %.2f Hz) AUX_DEV_ID=0x%02X (i2c 0x%02X) "
                  "AUX_IF_CONF=0x%02X (%s) AUX_RD_ADDR=0x%02X PWR_CTRL=0x%02X\n",
                  imu0->getAddress(), aux_conf, (double)aux_hz, aux_dev, aux_dev >> 1, aux_if,
                  (aux_if & 0x80) ? "manual" : "data mode", aux_rd, pwr);
  }

  // M5Unified creates instance 1 only for a BMM150 answering on the main bus.
  if (mag1) {
    uint8_t id = 0, pw = 0, mode = 0, rxy = 0, rz = 0;
    mag1->readRegister(0x40, &id, 1);
    mag1->readRegister(0x4B, &pw, 1);
    mag1->readRegister(0x4C, &mode, 1);
    mag1->readRegister(0x51, &rxy, 1);
    mag1->readRegister(0x52, &rz, 1);
    Serial.printf("[DIAG] BMM150 on main bus addr=0x%02X id=0x%02X power=0x%02X mode=0x%02X rep_xy=%u(n=%u) rep_z=%u(n=%u)\n",
                  mag1->getAddress(), id, pw, mode, rxy, 2u * rxy + 1u, rz, rz + 1u);
  } else {
    Serial.println("[DIAG] BMM150 main-bus instance: NONE (magnetometer is only behind the BMI270 AUX bus)");
  }
  for (uint8_t addr = 0x10; addr <= 0x13; ++addr) {
    uint8_t v = 0;
    const bool ack = M5.In_I2C.readRegister(addr, 0x40, &v, 1, 400000);
    Serial.printf("[DIAG] main-bus probe 0x%02X: %s%s\n", addr, ack ? "ACK" : "no ack",
                  ack ? (v == 0x32 ? " chip_id=0x32 (BMM150)" : "") : "");
  }

  if (have_blob) {
    Serial.println("[DIAG] saved calibration:");
    printBlobSummary(Serial, blob);
    printBlobDetail(Serial, blob);
  } else {
    Serial.println("[DIAG] no saved calibration in NVS");
  }
  Serial.printf("[DIAG] saved mag calibration applied at runtime: %s\n", cals.mag.ok ? "YES" : "NO");
}

static double expectedVertical() { return DIAG_EXPECTED_FIELD_UT * sin(DIAG_EXPECTED_DIP_DEG / kRadToDeg); }
static double expectedHorizontal() { return DIAG_EXPECTED_FIELD_UT * cos(DIAG_EXPECTED_DIP_DEG / kRadToDeg); }

// Calibrated field split along and across gravity, against the local reference.
static void fieldVsReference(const V3& m, const V3& a, double& vert, double& horiz) {
  const V3 down = -a.normalized();
  vert = m.dot(down);
  horiz = (m - down * float(vert)).norm();
}

// Body axis pointing up ("+z up" is screen down), or "tilt".
static const char* poseFace(const V3& a) {
  static const char* names[6] = {"+x up", "-x up", "+y up", "-y up", "+z up", "-z up"};
  for (int j = 0; j < 3; ++j) {
    const float c = a(j) / a.norm();
    if (c > 0.9f) return names[2 * j];
    if (c < -0.9f) return names[2 * j + 1];
  }
  return "tilt";
}

static void poseSample(const V3& a, const V3* m_cal, float gyro_rad_s, uint32_t now) {
  if (gyro_rad_s > QPOSE_MOVED_RAD_S) qw_moved = true;
  if (gyro_rad_s > QPOSE_STILL_RAD_S || !qw_moved || n_qposes >= MAX_QPOSES) {
    qw_acc = Stats3{}; qw_mag = Stats3{}; qw_start_ms = now;
    return;
  }
  qw_acc.add(a);
  if (m_cal) qw_mag.add(*m_cal);
  if (uint32_t(now - qw_start_ms) < QPOSE_MS || qw_mag.n < 20) return;
  QPose& p = qposes[n_qposes++];
  p.a = qw_acc.avg(); p.a_sd = qw_acc.sd(); p.m = qw_mag.avg(); p.m_sd = qw_mag.sd();
  double vert = 0, horiz = 0;
  fieldVsReference(p.m, p.a, vert, horiz);
  Serial.printf("[POSE] %2d %-5s |a|-g=%+.4f m/s^2 a_sd=[%.4f %.4f %.4f] | m=[%6.2f %6.2f %6.2f] |m|-F=%+5.2f "
                "V=%5.2f (%+5.2f) H=%5.2f (%+5.2f) dip=%.1f m_sd=[%.2f %.2f %.2f] n=%d\n",
                n_qposes, poseFace(p.a), p.a.norm() - ImuCalCfg::g_cal_local, p.a_sd.x(), p.a_sd.y(), p.a_sd.z(),
                p.m.x(), p.m.y(), p.m.z(), p.m.norm() - DIAG_EXPECTED_FIELD_UT, vert, vert - expectedVertical(),
                horiz, horiz - expectedHorizontal(), dipDeg(p.m, p.a), p.m_sd.x(), p.m_sd.y(), p.m_sd.z(), qw_mag.n);
  ui.title("POSES");
  M5.Display.printf("pose %d %s\n|a|-g %+.3f\n|m|-F %+.1f\n\nmove to next\ntap: done\n", n_qposes, poseFace(p.a),
                    p.a.norm() - ImuCalCfg::g_cal_local, p.m.norm() - DIAG_EXPECTED_FIELD_UT);
  qw_acc = Stats3{}; qw_mag = Stats3{}; qw_start_ms = now; qw_moved = false;
}

static void printPoseFit(const char* tag, const ScaleBiasFit& f, const char* unit, double k_mul, const char* k_unit) {
  if (!f.ok) {
    Serial.printf("[QUAL] %s fit: NOT DETERMINED (%d equations; visit every face, both signs of each axis)\n", tag,
                  f.equations);
    return;
  }
  Serial.printf("[QUAL] %s scale error [%+.3f %+.3f %+.3f] +- [%.3f %.3f %.3f] %s, residual bias [%+.4f %+.4f %+.4f] "
                "+- [%.4f %.4f %.4f] %s, rms %.4f %s\n", tag,
                f.k.x() * k_mul, f.k.y() * k_mul, f.k.z() * k_mul, f.k_sd.x() * k_mul, f.k_sd.y() * k_mul,
                f.k_sd.z() * k_mul, k_unit, f.b.x(), f.b.y(), f.b.z(), f.b_sd.x(), f.b_sd.y(), f.b_sd.z(), unit,
                f.rms, unit);
}

static void finishPoses() {
  Serial.println("[QUAL] ------------------------------------------------");
  Serial.printf("[QUAL] still poses: %d; magnetometer: %s\n", n_qposes, magLabel());
  if (n_qposes == 0) {
    Serial.println("[QUAL] no poses captured (hold each still 3 s, turn between them)");
    Serial.println("[QUAL] ------------------------------------------------");
    return;
  }
  V3 a[MAX_QPOSES], m[MAX_QPOSES];
  Eigen::Vector3d a_var = Eigen::Vector3d::Zero(), m_var = Eigen::Vector3d::Zero();
  char faces[64] = "";
  for (int i = 0; i < n_qposes; ++i) {
    a[i] = qposes[i].a; m[i] = qposes[i].m;
    a_var += qposes[i].a_sd.cast<double>().cwiseAbs2();
    m_var += qposes[i].m_sd.cast<double>().cwiseAbs2();
    const char* f = poseFace(a[i]);
    if (strcmp(f, "tilt") && !strstr(faces, f)) { strcat(faces, f); strcat(faces, " "); }
  }
  const V3 a_sd = (a_var / n_qposes).cwiseSqrt().cast<float>(), m_sd = (m_var / n_qposes).cwiseSqrt().cast<float>();
  Serial.printf("[QUAL] faces visited: %s\n", faces[0] ? faces : "none");
  const float g = ImuCalCfg::g_cal_local;
  const ScaleBiasFit fa = scaleBiasFit(a, a, n_qposes, g, 0.0f, false);
  Serial.printf("[QUAL] accel noise sd [%.4f %.4f %.4f] m/s^2 (single sample); |a|-g mean %+.4f sd %.4f m/s^2 over poses "
                "(g=%.5f)\n", a_sd.x(), a_sd.y(), a_sd.z(), fa.norm_err_mean, fa.norm_err_sd, (double)g);
  printPoseFit("accel", fa, "m/s^2", 1000.0, "ppt");
  const float F = DIAG_EXPECTED_FIELD_UT, V = float(expectedVertical());
  const ScaleBiasFit fm = scaleBiasFit(m, a, n_qposes, F, V, false);
  const ScaleBiasFit fmv = scaleBiasFit(m, a, n_qposes, F, V, true);
  Serial.printf("[QUAL] mag noise sd [%.3f %.3f %.3f] uT (single reading); |m|-F mean %+.2f sd %.2f uT over poses "
                "(F=%.1f); vertical-V mean %+.2f sd %.2f uT (V=%.2f)\n", m_sd.x(), m_sd.y(), m_sd.z(), fm.norm_err_mean,
                fm.norm_err_sd, (double)F, fmv.vert_err_mean, fmv.vert_err_sd, (double)V);
  printPoseFit("mag |m|", fm, "uT", 100.0, "%");
  printPoseFit("mag |m|+vertical", fmv, "uT", 100.0, "%");
  // Common scale is the calibration's field strength against the local field
  // (or a disturbed field indoors); a difference between axes is a scale error.
  const ScaleBiasFit& fb = fmv.ok ? fmv : fm;
  if (fb.ok) {
    const float common = fb.k.mean();
    Serial.printf("[QUAL] mag common scale %+.1f%% (calibrated field %.1f uT vs local %.1f uT)\n", 100.0f * common,
                  (double)(F * (1.0f + common)), (double)F);
    // Against the median axis, so one wrong axis does not mark the other two.
    int order[3] = {0, 1, 2};
    std::sort(order, order + 3, [&](int i, int j) { return fb.k(i) < fb.k(j); });
    const int med = order[1];
    for (int j = 0; j < 3; ++j) {
      const float rel = fb.k(j) - fb.k(med);
      const float sd = j == med ? 0.0f : std::hypot(fb.k_sd(j), fb.k_sd(med));
      Serial.printf("[QUAL] mag %c scale vs median axis %+.1f%% +- %.1f%%%s\n", "xyz"[j], 100.0f * rel, 100.0f * sd,
                    fabsf(rel) > 0.03f && !(fabsf(rel) < 3.0f * sd) ? " (SCALE MISMATCH)" : "");
    }
  }
  Serial.println("[QUAL] ------------------------------------------------");
}

static void drawWait(const char* t, const char* a, const char* b, const char* c) {
  ui.setReadRotation();
  ui.title(t);
  ui.line(a);
  ui.line(b);
  ui.line(c);
  ui.line("");
  ui.line("Tap or send 'n'");
}

static void enter(Phase p) {
  phase = p;
  phase_ms = millis();
  switch (p) {
    case Phase::WAIT_STILL:
      drawWait("STILL", "Lay flat, screen up", "away from metal,", "then tap: 10 s");
      Serial.println("[DIAG] next: STILL. Lay the board flat, screen up, away from metal; tap and do not touch.");
      Serial.println("[DIAG] (serial: 'q' read mag setting, 'o' driver default, 'p' low-noise, 'd' drift monitor)");
      break;
    case Phase::STILL:
      st_raw = Stats3{}; st_cal = Stats3{}; st_acc = Stats3{}; st_gyr = Stats3{}; st_dip = Stats3{};
      st_hdg_c = st_hdg_s = 0;
      Serial.printf("[STILL] magnetometer mode: %s\n", magLabel());
      st_intervals = 0; st_int_min = UINT32_MAX; st_int_max = 0; st_int_sum = 0;
      mask_ag = mask_g_only = mask_a_only = mask_none = mask_mag = 0;
      ui.title("STILL");
      ui.line("Do not touch");
      break;
    case Phase::WAIT_POSES:
      drawWait("POSES", "Set down still 3 s", "on each face + tilts", "tap when done");
      Serial.println("[DIAG] next: POSES. Tap, then set the board down still for 3 s on each of its 6 faces and a few tilts,");
      Serial.println("[DIAG] at different headings, away from metal. Each pose prints [POSE]; tap to finish with [QUAL].");
      break;
    case Phase::POSES:
      n_qposes = 0;
      qw_acc = Stats3{}; qw_mag = Stats3{}; qw_start_ms = millis(); qw_moved = true;
      Serial.printf("[POSE] reference: g=%.5f m/s^2, field %.1f uT, dip %.1f deg (H %.2f, V %.2f uT); mag mode: %s\n",
                    (double)ImuCalCfg::g_cal_local, (double)DIAG_EXPECTED_FIELD_UT, (double)DIAG_EXPECTED_DIP_DEG,
                    expectedHorizontal(), expectedVertical(), magLabel());
      ui.title("POSES");
      ui.line("Hold each still 3 s");
      break;
    case Phase::WAIT_ROTATE:
      drawWait("ROTATE", "Tap, then SLOWLY", "turn + tilt ALL ways:", "upside down, edges");
      Serial.println("[DIAG] next: ROTATE. Turn and tilt slowly (under ~45 deg/s) through every direction, including upside down.");
      break;
    case Phase::ROTATE:
      magCal.clear();
      cap_cfg = imu_cal::MagSampleWindow::captureCfg(false);
      cap_cfg.stuck_ms = 12000;
      cap_cfg.min_delta_uT = 0.03f;
      cap_cfg.span_min_frac = 0.35f;
      cap_cfg.span_mid_frac = 0.55f;
      cap_cfg.urange_target = 1.05f;
      cap_cfg.timeout_ms = ROTATE_MAX_MS + 10000;
      delete capture;
      capture = new imu_cal::MagCapture<float, 400>(magCal, cap_cfg);
      capture->begin(millis());
      window.begin();
      cap_status = imu_cal::MagCaptureStatus::CAPTURING;
      n_poses = pose_slot = n_pairs = pair_slot = 0;
      ring_n = ring_head = 0;
      pw_mag = Stats3{}; pw_acc = Stats3{}; pw_max_gyro = 0; pw_start_ms = millis();
      ui.title("ROTATE");
      ui.line("Turn + tilt slowly");
      ui.line("turn slowly");
      ui.line("tap to stop early");
      break;
    case Phase::ANALYZE:
      ui.title("ANALYZE");
      ui.line("Working...");
      break;
    case Phase::LIVE:
      Serial.println("[DIAG] LIVE: hold level and tilt +-20 deg at E/W/N/S; compare hdg columns. Tap to run again.");
      break;
    case Phase::DRIFT:
      dr_raw = Stats3{}; dr_acc = Stats3{}; dr_temp_sum = 0; dr_temp_n = 0; dr_have_first = false;
      dr_start_ms = millis();
      ui.title("DRIFT");
      ui.line("Do not touch");
      ui.line("for 3+ minutes");
      ui.line("tap or 'd' to stop");
      Serial.printf("[DRIFT] monitoring the untouched board every %u s (magnetometer mode: %s); tap or 'd' to stop\n",
                    (unsigned)(DRIFT_PERIOD_MS / 1000), magLabel());
      break;
  }
}

// Tilt-compensated heading and dip for one model, for LIVE printing.
static void liveColumns(const MagModel& model, const Matrix3f& R, float& hdg, float& dip, float& norm) {
  const V3 m = R * model.apply(m_raw);
  hdg = headingDeg(a_cal, m);
  dip = dipDeg(m, a_cal);
  norm = m.norm();
}

static void analyzeTask(void*) {
  Serial.printf("[FIT] fitting fresh calibration (magnetometer mode: %s)...\n", magLabel());
  fresh = imu_cal::MagCalibration<float>{};
  fresh_ok = magCal.fit(fresh, 3, 0.15f, 1e-6f, &fresh_reason);
  const auto& q = magCal.quality;
  Serial.printf("[FIT] ok=%d reason=%s gate=%s rms=%.4f p95=%.4f inliers=%d/%d cells=%d field=%.2f\n",
                (int)fresh_ok, imu_cal::fitFailStr(fresh_reason), imu_cal::magFitGateText(q.gate),
                q.rms, q.p95, q.inliers, q.samples, q.cells, (double)fresh.field_uT);
  fresh_model.A = fresh.A; fresh_model.b = fresh.b;
  if (fresh_ok) {
    printMat("[FIT] fresh A", fresh.A);
    Serial.printf("[FIT] fresh b [%.3f %.3f %.3f]\n", fresh.b.x(), fresh.b.y(), fresh.b.z());
  }
  if (saved_ok) {
    printMat("[FIT] saved A", saved_model.A);
    Serial.printf("[FIT] saved b [%.3f %.3f %.3f]\n", saved_model.b.x(), saved_model.b.y(), saved_model.b.z());
    if (fresh_ok) {
      Serial.printf("[FIT] saved-vs-fresh: |dA|/|A|=%.3f |db|=%.2f\n",
                    (double)((saved_model.A - fresh_model.A).norm() / fresh_model.A.norm()),
                    (double)(saved_model.b - fresh_model.b).norm());
    }
  }

  if (fresh_ok && prev_ok) {
    Serial.printf("[FIT] previous-run-vs-fresh: |dA|/|A|=%.3f |db|=%.2f b_prev=[%.2f %.2f %.2f] temp %.1f -> %.1f C\n",
                  (double)((prev_model.A - fresh_model.A).norm() / fresh_model.A.norm()),
                  (double)(prev_model.b - fresh_model.b).norm(), prev_model.b.x(), prev_model.b.y(), prev_model.b.z(),
                  (double)prev_temp_c, (double)last_temp_c);
  }
  Serial.printf("[FIELD] poses collected: %d (each a %u ms window turning slower than 57 deg/s)\n", n_poses,
                (unsigned)POSE_WINDOW_MS);
  FieldStats fs_saved, fs_fresh;
  if (saved_ok) {
    fs_saved = fieldStats(poses, n_poses, saved_model);
    Serial.printf("[FIELD] saved: |m| mean=%.2f sd=%.2f min=%.2f max=%.2f spread=%.1f%% dip mean=%.1f sd=%.1f\n",
                  fs_saved.norm_mean, fs_saved.norm_sd, fs_saved.norm_min, fs_saved.norm_max,
                  100.0 * fs_saved.spread(), fs_saved.dip_mean, fs_saved.dip_sd);
  }
  if (fresh_ok) {
    fs_fresh = fieldStats(poses, n_poses, fresh_model);
    Serial.printf("[FIELD] fresh: |m| mean=%.2f sd=%.2f min=%.2f max=%.2f spread=%.1f%% dip mean=%.1f sd=%.1f\n",
                  fs_fresh.norm_mean, fs_fresh.norm_sd, fs_fresh.norm_min, fs_fresh.norm_max,
                  100.0 * fs_fresh.spread(), fs_fresh.dip_mean, fs_fresh.dip_sd);
  }

  if (saved_ok && fresh_ok && n_poses > 0) {
    double sc = 0, ss = 0;
    for (int i = 0; i < n_poses; ++i) {
      const float d = wrap180(headingDeg(poses[i].a, saved_model.apply(poses[i].m_raw)) -
                              headingDeg(poses[i].a, fresh_model.apply(poses[i].m_raw))) / kRadToDeg;
      if (std::isfinite(d)) { sc += cos(d); ss += sin(d); }
    }
    Serial.printf("[FIELD] mean heading difference saved-minus-fresh over poses: %.1f deg\n", atan2(ss, sc) * kRadToDeg);
  }

  // Axis test needs offsets removed; prefer the fresh fit made here.
  const MagModel& axis_model = fresh_ok ? fresh_model : saved_model;
  const AxisResult ax = axisConsistency(pairs, n_pairs, poses, n_poses, axis_model,
                                        DIAG_EXPECTED_DIP_DEG >= 0.0f);
  char best_txt[16];
  permutationText(ax.best, best_txt);
  Serial.printf("[AXIS] model=%s pairs=%d best=%s rms=%.3f current-mapping rms=%.3f runner-up rms=%.3f dip_with_best=%.1f\n",
                fresh_ok ? "fresh" : (saved_ok ? "saved" : "none"), ax.pairs, best_txt, ax.best_rms,
                ax.identity_rms, ax.second_rms, ax.best_dip);

  Alignment al_saved;
  if (fresh_ok) fresh_align = solveAlignment(poses, n_poses, fresh_model);
  if (saved_ok) al_saved = solveAlignment(poses, n_poses, saved_model);
  if (fresh_align.ok) {
    const Matrix3f Rerr = fresh_align.R.transpose();
    Serial.printf("[ALIGN] fresh: mag->accel rotation [%.2f %.2f %.2f] deg (|%.2f|) dip=%.1f sd %.2f -> %.2f\n",
                  fresh_align.rotvec_deg.x(), fresh_align.rotvec_deg.y(), fresh_align.rotvec_deg.z(),
                  fresh_align.rotvec_deg.norm(), fresh_align.dip_deg, fresh_align.dip_sd_before,
                  fresh_align.dip_sd_after);
    Serial.printf("[ALIGN] heading error this rotation causes: level %.1f deg, at 20 deg tilt %.1f deg\n",
                  headingErrorDeg(Rerr, fresh_align.dip_deg, 0.0f), headingErrorDeg(Rerr, fresh_align.dip_deg, 20.0f));
  }
  if (al_saved.ok) {
    Serial.printf("[ALIGN] saved: rotation [%.2f %.2f %.2f] deg dip=%.1f sd %.2f -> %.2f\n",
                  al_saved.rotvec_deg.x(), al_saved.rotvec_deg.y(), al_saved.rotvec_deg.z(),
                  al_saved.dip_deg, al_saved.dip_sd_before, al_saved.dip_sd_after);
  }

  // Verdicts.
  Serial.println("[VERDICT] ------------------------------------------------");
  Serial.printf("[VERDICT] magnetometer setting during this run: %s\n", magLabel());
  if (ax.pairs < 50) {
    Serial.printf("[VERDICT] axes: INCONCLUSIVE (%d rotating pairs; turn more during ROTATE)\n", ax.pairs);
  } else if (ax.best == 0 && ax.identity_rms < 0.5f * ax.second_rms) {
    Serial.printf("[VERDICT] axes: OK, magnetometer axes agree with the gyro (rms %.3f vs %.3f)\n",
                  ax.identity_rms, ax.second_rms);
  } else {
    Serial.printf("[VERDICT] axes: MISMATCH, gyro agrees with mag axes (%s) of the current mapping (rms %.3f vs current %.3f)\n",
                  best_txt, ax.best_rms, ax.identity_rms);
  }
  auto fits = [](const FieldStats& f) { return f.n >= 20 && f.spread() < 0.10f && f.dip_sd < 3.0f; };
  auto verdict = [&](const FieldStats& f) { return f.n < 20 ? "INCONCLUSIVE (too few poses)" : (fits(f) ? "FITS" : "DOES NOT FIT"); };
  if (saved_ok) {
    Serial.printf("[VERDICT] saved calibration here: %s (|m| spread %.1f%%, dip %.1f+-%.1f, expected %.1f)\n",
                  verdict(fs_saved), 100.0 * fs_saved.spread(), fs_saved.dip_mean,
                  fs_saved.dip_sd, (double)DIAG_EXPECTED_DIP_DEG);
  } else {
    Serial.println("[VERDICT] saved calibration here: NONE / not applied");
  }
  if (fresh_ok) {
    Serial.printf("[VERDICT] fresh calibration here: %s (|m| spread %.1f%%, dip %.1f+-%.1f)\n",
                  verdict(fs_fresh), 100.0 * fs_fresh.spread(), fs_fresh.dip_mean,
                  fs_fresh.dip_sd);
  } else {
    Serial.printf("[VERDICT] fresh calibration here: FIT FAILED (%s); capture more directions\n",
                  imu_cal::fitFailStr(fresh_reason));
  }
  if (fresh_align.ok) {
    const float err20 = headingErrorDeg(fresh_align.R.transpose(), fresh_align.dip_deg, 20.0f);
    Serial.printf("[VERDICT] mag/accel frame misalignment: %.1f deg -> about %.1f deg heading error at 20 deg tilt%s\n",
                  fresh_align.rotvec_deg.norm(), err20, err20 > 3.0f ? " (SIGNIFICANT)" : "");
    Serial.printf("[VERDICT] dip after alignment %.1f vs expected %.1f%s\n", fresh_align.dip_deg,
                  (double)DIAG_EXPECTED_DIP_DEG,
                  fabsf(fresh_align.dip_deg - DIAG_EXPECTED_DIP_DEG) > 5.0f ? " (CHECK local field / nearby metal)" : "");
  }
  if (fs_fresh.n < 20) {
    Serial.println("[VERDICT] => field checks need more poses: turn more slowly during ROTATE.");
  } else if (saved_ok && fresh_ok && fits(fs_fresh) && !fits(fs_saved)) {
    Serial.println("[VERDICT] => the saved calibration is stale for this place/mounting; rerun the wizard here.");
  } else if (fresh_ok && !fits(fs_fresh)) {
    Serial.println("[VERDICT] => even a fresh fit does not describe this sensor; capture/fit or sensor problem.");
  }
  Serial.println("[VERDICT] ------------------------------------------------");
  if (fresh_ok) { prev_model = fresh_model; prev_ok = true; prev_temp_c = last_temp_c; }
  analysis_done = true;
  vTaskDelete(nullptr);
}

static void onSample(const ImuSample& s) {
  const float dt = last_sample_us ? (s.sample_us - last_sample_us) * 1e-6f : 1.0f / LOOP_HZ;
  last_sample_us = s.sample_us;
  const float tempC = s.tempC;
  last_temp_c = tempC;
  a_cal = cals.applyAccel(s.a, tempC);
  w_cal = cals.applyGyro(s.w, tempC);
  if (dt > 0.0f && dt < 0.05f) {
    const V3 dth = w_cal * dt;
    const float ang = dth.norm();
    if (ang > 1e-9f) q_gyro = (q_gyro * Eigen::Quaternionf(Eigen::AngleAxisf(ang, dth / ang))).normalized();
  }

  // M5Unified returns the cached AUX value between real updates; a changed
  // raw value is a new magnetometer reading.
  const bool fresh_m = !have_m || (s.m - last_distinct_m).norm() > 0.0f;
  if (fresh_m) {
    if (have_m && phase == Phase::STILL) {
      const uint32_t di = s.sample_us - last_m_change_us;
      ++st_intervals; st_int_sum += di;
      st_int_min = di < st_int_min ? di : st_int_min;
      st_int_max = di > st_int_max ? di : st_int_max;
    }
    if (phase == Phase::ROTATE) {
      ring_m[ring_head] = s.m;
      ring_q[ring_head] = q_gyro;
      const int oldest = (ring_head + 1) % (PAIR_SPAN + 1);
      ring_head = oldest;
      ring_n = ring_n < PAIR_SPAN + 1 ? ring_n + 1 : PAIR_SPAN + 1;
      if (ring_n == PAIR_SPAN + 1) {
        // Rotation of the body between the two readings, in the older body frame.
        const Eigen::AngleAxisf rel(ring_q[oldest].conjugate() * q_gyro);
        pairs[pair_slot] = GyroMagPair{ring_m[oldest], s.m, V3(rel.axis() * rel.angle())};
        pair_slot = (pair_slot + 1) % MAX_PAIRS;
        n_pairs = n_pairs < MAX_PAIRS ? n_pairs + 1 : MAX_PAIRS;
      }
    }
    last_distinct_m = s.m;
    last_m_change_us = s.sample_us;
    have_m = true;
  }
  m_raw = s.m;

  if (phase == Phase::DRIFT) {
    dr_acc.add(a_cal);
    if (fresh_m) dr_raw.add(s.m);
    if (std::isfinite(tempC)) { dr_temp_sum += tempC; ++dr_temp_n; }
  }

  if (phase == Phase::STILL) {
    st_acc.add(a_cal);
    st_gyr.add(w_cal);
    if (fresh_m) {
      st_raw.add(s.m);
      const Vector3f mc = cals.applyMag(s.m);
      st_cal.add(mc);
      const float h = headingDeg(a_cal, mc) / kRadToDeg;
      if (std::isfinite(h)) { st_hdg_c += cosf(h); st_hdg_s += sinf(h); }
      st_dip.add(V3(dipDeg(mc, a_cal), 0, 0));
    }
  }

  if (phase == Phase::POSES) {
    const Vector3f mc = cals.applyMag(s.m);
    poseSample(a_cal, fresh_m ? &mc : nullptr, w_cal.norm(), millis());
  }

  if (phase == Phase::ROTATE) {
    const uint32_t now = millis();
    Vector3f mean;
    Eigen::Matrix3f cov;
    const bool valid = s.m.allFinite() && s.m.norm() > 1.0f;
    // Like the wizard, stop adding fit samples once the capture is complete.
    if (cap_status != imu_cal::MagCaptureStatus::READY) {
      const bool averaged = window.update(now, valid ? &s.m : nullptr, mean, cov);
      cap_status = capture->update(now, averaged ? &mean : nullptr, valid ? &s.m : nullptr,
                                   averaged ? &cov : nullptr);
    }
    // Still-pose windows for the field, dip and alignment checks.
    if (fresh_m) pw_mag.add(s.m);
    pw_acc.add(a_cal);
    const float wn = w_cal.norm();
    pw_max_gyro = wn > pw_max_gyro ? wn : pw_max_gyro;
    if (uint32_t(now - pw_start_ms) >= POSE_WINDOW_MS) {
      const float an = pw_acc.avg().norm();
      // Slow hand motion is fine: magnetometer lag at 1 rad/s is a few degrees.
      if (pw_mag.n >= 3 && pw_max_gyro < 1.0f && fabsf(an - ImuCalCfg::g_cal_local) < 1.0f) {
        poses[pose_slot] = Pose{pw_mag.avg(), pw_acc.avg()};
        pose_slot = (pose_slot + 1) % MAX_POSES;
        n_poses = n_poses < MAX_POSES ? n_poses + 1 : MAX_POSES;
      }
      pw_mag = Stats3{}; pw_acc = Stats3{}; pw_max_gyro = 0; pw_start_ms = now;
    }
  }
}

static void finishStill() {
  const V3 nr = st_raw.sd(), nc = st_cal.sd(), am = st_acc.avg();
  const V3 mc = st_cal.avg();
  const float ms = st_intervals ? (float)st_int_sum / st_intervals / 1000.0f : NAN;
  const uint32_t polls = mask_ag + mask_g_only + mask_a_only + mask_none;
  Serial.printf("[STILL] polls=%lu accel+gyro=%lu gyro-only=%lu accel-only=%lu none=%lu mag-bit=%lu\n",
                (unsigned long)polls, (unsigned long)mask_ag, (unsigned long)mask_g_only,
                (unsigned long)mask_a_only, (unsigned long)mask_none, (unsigned long)mask_mag);
  Serial.printf("[STILL] distinct mag readings=%d interval mean=%.1f ms min=%.1f max=%.1f ms (=> %.1f Hz)\n",
                st_raw.n, (double)ms, st_intervals ? st_int_min / 1000.0 : NAN,
                st_intervals ? st_int_max / 1000.0 : NAN, (double)(ms > 0 ? 1000.0f / ms : NAN));
  Serial.printf("[STILL] raw mag mean [%.2f %.2f %.2f] noise sd [%.2f %.2f %.2f]\n",
                st_raw.avg().x(), st_raw.avg().y(), st_raw.avg().z(), nr.x(), nr.y(), nr.z());
  Serial.printf("[STILL] cal mag mean [%.2f %.2f %.2f] noise sd [%.2f %.2f %.2f] |m|=%.2f\n",
                mc.x(), mc.y(), mc.z(), nc.x(), nc.y(), nc.z(), mc.norm());
  Serial.printf("[STILL] accel mean [%.3f %.3f %.3f] |a|=%.3f tilt=%.1f deg; gyro mean [%.3f %.3f %.3f] deg/s\n",
                am.x(), am.y(), am.z(), am.norm(),
                acosf(fminf(1.0f, fabsf(am.z()) / fmaxf(am.norm(), 1e-6f))) * kRadToDeg,
                st_gyr.avg().x() * kRadToDeg, st_gyr.avg().y() * kRadToDeg, st_gyr.avg().z() * kRadToDeg);
  const double R = st_raw.n ? sqrt(st_hdg_c * st_hdg_c + st_hdg_s * st_hdg_s) / st_raw.n : 0.0;
  double vert = 0, horiz = 0;
  fieldVsReference(mc, am, vert, horiz);
  const V3 na = st_acc.sd();
  Serial.printf("[STILL] accel noise sd [%.4f %.4f %.4f] m/s^2, |a|-g=%+.4f m/s^2 (g=%.5f)\n", na.x(), na.y(), na.z(),
                am.norm() - ImuCalCfg::g_cal_local, (double)ImuCalCfg::g_cal_local);
  Serial.printf("[STILL] cal field vs local: |m|=%.2f (F %.1f, %+.1f%%) vertical=%.2f (%.2f, %+.1f%%) "
                "horizontal=%.2f (%.2f, %+.1f%%)\n", mc.norm(), (double)DIAG_EXPECTED_FIELD_UT,
                100.0 * (mc.norm() / DIAG_EXPECTED_FIELD_UT - 1.0), vert, expectedVertical(),
                100.0 * (vert / expectedVertical() - 1.0), horiz, expectedHorizontal(),
                100.0 * (horiz / expectedHorizontal() - 1.0));
  Serial.printf("[STILL] single-reading scatter with saved cal: heading sd=%.1f deg, dip sd=%.1f deg; temp=%.1f C\n",
                R > 0 ? sqrt(-2.0 * log(R)) * kRadToDeg : NAN, st_dip.sd().x(), (double)last_temp_c);
  Serial.printf("[STILL] dip with saved cal=%.1f (expected %.1f), heading=%.1f; accel D is %s (screen %s)\n",
                dipDeg(mc, am), (double)DIAG_EXPECTED_DIP_DEG, headingDeg(am, mc),
                am.z() < 0 ? "negative" : "positive", am.z() < 0 ? "up" : "down");
}

void setup() {
  delay(200);
  Serial.setTxBufferSize(8192);
  Serial.begin(115200);
  Serial.setTxTimeoutMs(0);
  delay(800);

  auto cfg = M5.config();
  M5.begin(cfg);
  clearM5UnifiedImuCalibration();
  if (configureAtomS3RMag(Serial)) snprintf(mag_mode, sizeof(mag_mode), "low-noise (47/41 reps, 30 Hz)");
  Serial.printf("[DIAG] magnetometer: %s\n", magLabel());
  delay(250);
  ui.begin();

  have_blob = store.load(blob);
  if (!have_blob) memset(&blob, 0, sizeof(blob));
  cals.rebuildFromBlob(blob);
  saved_ok = cals.mag.ok;
  saved_model.A = cals.mag.A;
  saved_model.b = cals.mag.b;

  if (!M5.Imu.isEnabled()) {
    Serial.println("[DIAG] IMU not found");
    ui.title("IMU");
    ui.line("Not found");
    while (true) delay(100);
  }
  bootReport();
  enter(Phase::WAIT_STILL);
  next_tick_us = micros();
}

void loop() {
  while ((int32_t)(next_tick_us - micros()) > 0) delayMicroseconds(200);
  next_tick_us += LOOP_PERIOD_US;
  if ((int32_t)(micros() - next_tick_us) > (int32_t)(4 * LOOP_PERIOD_US)) next_tick_us = micros() + LOOP_PERIOD_US;

  Input::update();
  bool advance = Input::tapPressed();
  while (Serial.available()) {
    const int c = Serial.read();
    if (c == 'n' || c == 'N') advance = true;
    // Register access only between captures, never while sampling.
    const bool idle = phase == Phase::WAIT_STILL || phase == Phase::WAIT_ROTATE || phase == Phase::LIVE;
    const char lc = char(tolower(c));
    if ((lc == 'p' || lc == 'q' || lc == 'o' || lc == 'c') && !idle)
      Serial.println("[MAGCFG] busy; send 'p'/'o'/'q'/'c' while waiting for a tap or in LIVE");
    else if (lc == 'p' || lc == 'q' || lc == 'o') magCommand(lc);
    else if (lc == 'c') toggleCompensation();
    else if ((c == 'd' || c == 'D') && phase == Phase::DRIFT) enter(dr_return);
    else if (c == 'd' || c == 'D') { dr_return = idle ? phase : Phase::WAIT_STILL; if (idle) enter(Phase::DRIFT); else Serial.println("[DRIFT] busy; send 'd' while waiting for a tap or in LIVE"); }
  }

  const uint32_t sample_us = micros();
  const uint32_t mask = M5.Imu.update();
  if (phase == Phase::STILL) {
    const bool a = mask & ATOMS3R_IMU_MASK_ACCEL, g = mask & ATOMS3R_IMU_MASK_GYRO;
    if (a && g) ++mask_ag; else if (g) ++mask_g_only; else if (a) ++mask_a_only; else ++mask_none;
    if (mask & (1u << 2)) ++mask_mag;
  }
  ImuSample s;
  if (readImuMapped(M5.Imu, mask, sample_us, s)) onSample(s);

  const uint32_t now = millis();
  switch (phase) {
    case Phase::WAIT_STILL:
      if (advance) enter(Phase::STILL);
      break;
    case Phase::STILL:
      if (now - last_draw_ms >= 250) { ui.bar01((now - phase_ms) / (float)STILL_MS); last_draw_ms = now; }
      if (now - phase_ms >= STILL_MS) { finishStill(); enter(Phase::WAIT_POSES); }
      break;
    case Phase::WAIT_POSES:
      if (advance) enter(Phase::POSES);
      break;
    case Phase::POSES:
      if (advance) { finishPoses(); enter(Phase::WAIT_ROTATE); }
      break;
    case Phase::WAIT_ROTATE:
      if (advance) enter(Phase::ROTATE);
      break;
    case Phase::ROTATE: {
      const bool ready = cap_status == imu_cal::MagCaptureStatus::READY;
      if (now - last_print_ms >= 5000) {
        const auto c = capture->coverage();
        Serial.printf("[ROT] t=%.0fs samples=%d/%d cells=%d ratios=(%.2f,%.2f) poses=%d pairs=%d status=%d\n",
                      (now - phase_ms) / 1000.0, magCal.buf.n, cap_cfg.required_samples, c.cells,
                      (double)c.min_ratio, (double)c.mid_ratio, n_poses, n_pairs, (int)cap_status);
        last_print_ms = now;
      }
      if (now - last_draw_ms >= 250) {
        ui.bar01(capture->progress(now));
        M5.Display.setCursor(0, 70);
        M5.Display.printf("fit %3d/%d   \nposes %3d   \n%s        ", magCal.buf.n, cap_cfg.required_samples,
                          n_poses, capture->hint(now));
        last_draw_ms = now;
      }
      const bool stale = cap_status == imu_cal::MagCaptureStatus::STALE;
      if ((ready && n_poses >= 60) || advance || stale || now - phase_ms >= ROTATE_MAX_MS) {
        Serial.printf("[ROT] done: status=%d samples=%d poses=%d pairs=%d%s\n", (int)cap_status, magCal.buf.n,
                      n_poses, n_pairs, stale ? " (magnetometer stream STALE)" : "");
        enter(Phase::ANALYZE);
        analysis_done = false;
        // Priority 0 time-slices with IDLE0, so the task watchdog stays fed during the fit.
        xTaskCreatePinnedToCore(analyzeTask, "magdiag", 32768 / sizeof(StackType_t), nullptr, tskIDLE_PRIORITY, nullptr, 0);
      }
      break;
    }
    case Phase::ANALYZE:
      if (analysis_done) enter(Phase::LIVE);
      break;
    case Phase::DRIFT:
      if (advance) { enter(dr_return); break; }
      if (now - phase_ms >= DRIFT_PERIOD_MS && dr_raw.n > 0) {
        const V3 r = dr_raw.avg(), a = dr_acc.avg(), c = cals.applyMag(r);
        const float temp = dr_temp_n ? float(dr_temp_sum / dr_temp_n) : NAN;
        if (!dr_have_first) { dr_first_raw = r; dr_first_temp = temp; dr_have_first = true; }
        const V3 d = r - dr_first_raw;
        Serial.printf("[DRIFT] t=%5.0fs temp=%.2f C (d%+.2f) raw=[%7.2f %7.2f %7.2f] d=[%+6.2f %+6.2f %+6.2f] "
                      "noise=[%.2f %.2f %.2f] hdg=%.1f dip=%.1f |m|=%.2f tilt=%.1f\n",
                      (now - dr_start_ms) / 1000.0, (double)temp, (double)(temp - dr_first_temp),
                      r.x(), r.y(), r.z(), d.x(), d.y(), d.z(), dr_raw.sd().x(), dr_raw.sd().y(), dr_raw.sd().z(),
                      headingDeg(a, c), dipDeg(c, a), c.norm(),
                      acosf(fminf(1.0f, fabsf(a.z()) / fmaxf(a.norm(), 1e-6f))) * kRadToDeg);
        if (now - last_draw_ms >= 1000) {
          ui.title("DRIFT");
          M5.Display.printf("T %.2f C\n\ndz %+.2f uT\ndx %+.2f\ndy %+.2f\n\ntap: stop\n",
                            (double)temp, (double)d.z(), (double)d.x(), (double)d.y());
          last_draw_ms = now;
        }
        dr_raw = Stats3{}; dr_acc = Stats3{}; dr_temp_sum = 0; dr_temp_n = 0;
        phase_ms = now;
      }
      break;
    case Phase::LIVE: {
      if (advance) { enter(Phase::WAIT_STILL); break; }
      if (now - last_print_ms >= 250) {
        float hs = NAN, ds = NAN, ns = NAN, hf = NAN, df = NAN, nf = NAN, ha = NAN, da = NAN, na = NAN;
        if (saved_ok) liveColumns(saved_model, Matrix3f::Identity(), hs, ds, ns);
        if (fresh_ok) liveColumns(fresh_model, Matrix3f::Identity(), hf, df, nf);
        if (fresh_ok && fresh_align.ok) liveColumns(fresh_model, fresh_align.R, ha, da, na);
        const float roll = atan2f(-a_cal.y(), -a_cal.z()) * kRadToDeg;
        const float pitch = atan2f(a_cal.x(), sqrtf(a_cal.y() * a_cal.y() + a_cal.z() * a_cal.z())) * kRadToDeg;
        Serial.printf("[LIVE] roll=%6.1f pitch=%6.1f | hdg saved=%6.1f fresh=%6.1f aligned=%6.1f | "
                      "dip saved=%5.1f fresh=%5.1f aligned=%5.1f | |m| saved=%5.1f fresh=%5.1f\n",
                      roll, pitch, hs, hf, ha, ds, df, da, ns, nf);
        if (now - last_draw_ms >= 500) {
          ui.title("LIVE");
          M5.Display.printf("R%5.0f P%5.0f\n\n", roll, pitch);
          M5.Display.printf("saved  %5.1f\n", hs);
          M5.Display.printf("fresh  %5.1f\n", hf);
          M5.Display.printf("align  %5.1f\n\n", ha);
          M5.Display.printf("dip %4.0f/%4.0f\n", ds, df);
          M5.Display.println("tap: run again");
          last_draw_ms = now;
        }
        last_print_ms = now;
      }
      break;
    }
  }
}
