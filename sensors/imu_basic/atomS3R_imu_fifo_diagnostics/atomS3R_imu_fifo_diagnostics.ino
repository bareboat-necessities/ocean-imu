/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R IMU FIFO diagnostics: timing and accuracy of the IMU data, measured
  on the device with the acquisition the INS sketches should move to. It
  estimates nothing and never changes the saved calibration.

  Acquisition (src/AtomS3R/AtomS3R_Bmi270Fifo.h): the BMI270 is taken over
  from M5Unified's polling path:
    - accelerometer and gyro both at 200 Hz (+-8 g, +-2000 deg/s, the ranges
      M5Unified's conversion assumes, so samples are in readImuMapped() units);
    - every sample and every magnetometer readout through the BMI270 FIFO,
      drained each loop, so nothing is lost while the display is busy;
    - each sample timed by its place on the sensor's own clock (a frame every
      128 sensortime ticks), mapped onto host micros();
    - the BMM150 in forced mode, triggered by the BMI270 on every 25 Hz AUX
      cycle, so each readout is one complete conversion with a fixed relation
      to the sensor clock.
  The saved calibration is applied to every sample, as the INS sketches do.

  Timing (T lines every 10 s):
    rate      200 Hz delivered on the sensor's own clock, no frame missing (the
              host-clock rate differs by the sensor oscillator's error,
              reported as clock_ppm)
    no loss   no lost frame while the display is loaded: the screen redraws at
              10 Hz and a display burst (400 ms, then 600 ms, alternating)
              runs every 10 s. The FIFO holds 0.73 s; 's <ms>' on serial runs
              a burst of any length (a longer one shows how an overflow is
              reported)
    spacing   delivered sample spacing between 4.9 and 5.1 ms
    mag       25 Hz readouts, each a new conversion (data-ready set), every
              one 40 ms after the previous on the sensor clock
    delay     the magnetometer delay measured from hand rotation by the
              calibration wizard's fit: three fits at most 4 ms wide agreeing
              within 4 ms (rock the board about one rock per second)

  Accuracy (src/imu_calibrate/ImuAccuracy.h; S, W and A lines). Put the board
  still on a firm surface in several orientations (each for 5 s or more) and
  turn it between them:
    accel     |a| in each still pose against the calibration's gravity
              target g_cal_local, in mg; noise
    gyro      the still mean rate (residual bias after calibration); noise;
              gyro scale, fitted from the rotation between two poses to the
              change of gravity (and field) direction
    mag       F, the vertical component Z and the dip in each pose against
              the World Magnetic Model WMM-2025 at the configured place and
              the build date (src/imu_calibrate/GeomagWmm.h); |m| in motion
              (constant for a good calibration) and the body directions
              covered
  The place defaults to the calibration site (CalibrationSite in
  AtomS3R_ImuUnits.h); set FIFO_DIAG_LATITUDE_DEG/LONGITUDE_DEG/HEIGHT_M, or
  send "L <lat> <lon> [height_m]" on serial. The model omits local
  disturbances: indoors, near steel or electronics, F can differ by several
  percent.

  USB serial (115200), one key=value line per record:
    H identity, calibration status      G the expected field and gravity
      and what the calibration holds    D live state, 1 Hz
    F FIFO setup read back              X one magnetometer-delay fit
    T timing criteria                   W one rotation check
    S one still pose                    E frames counted lost (the read that
    A accuracy summary and criteria       showed it)
    R every frame (calibrated), off until 'q'
  Tap the screen to page through TIME, ACC/GYR, MAG and ROT. Serial commands
  (one per line): 's [ms]' display burst (600 ms), 'r' restart the evaluation, 'q'
  start/stop the R lines, 'n' toggle the gyro noise-performance mode (and
  restart the evaluation), 'L <lat> <lon> [height_m]' set the place.
*/

#include <cstdlib>

#include <Arduino.h>
#include <M5Unified.h>
#include <ArduinoOceanImu.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"
#include "AtomS3R/AtomS3R_Bmi270Fifo.h"
#include "imu_calibrate/MagTimingCal.h"
#include "imu_calibrate/ImuAccuracy.h"
#include "imu_calibrate/GeomagWmm.h"

// Trigger-to-readout offset of the AUX interface in 2.5 ms units; 0 keeps a
// full AUX period between the forced-mode command and the readout.
#ifndef FIFO_DIAG_AUX_OFFSET
#define FIFO_DIAG_AUX_OFFSET 0
#endif
#ifndef FIFO_DIAG_LOAD_EVERY_MS
#define FIFO_DIAG_LOAD_EVERY_MS 10000
#endif
#ifndef FIFO_DIAG_LOAD_MS
#define FIFO_DIAG_LOAD_MS 400
#endif
// Within the 0.73 s the FIFO holds; 's <ms>' on serial runs any length.
#ifndef FIFO_DIAG_HEAVY_LOAD_MS
#define FIFO_DIAG_HEAVY_LOAD_MS 600
#endif
// R lines (every frame) at start; 'q' on serial toggles them.
#ifndef FIFO_DIAG_RAW_LINES
#define FIFO_DIAG_RAW_LINES 0
#endif
// BMI270 GYR_CONF at start: 200 Hz, normal filter, performance filter.
// Bit 6 (0x40) is the gyro noise-performance mode (0 power optimised, the
// M5Unified/INS setting; 1 performance optimised); 'n' on serial toggles it.
#ifndef FIFO_DIAG_GYR_CONF
#define FIFO_DIAG_GYR_CONF 0xA9
#endif
// Place for the expected magnetic field (geodetic, height above the WGS84
// ellipsoid); defaults to the calibration site.
#ifndef FIFO_DIAG_LATITUDE_DEG
#define FIFO_DIAG_LATITUDE_DEG (atoms3r_ical::CalibrationSite::latitude_deg)
#endif
#ifndef FIFO_DIAG_LONGITUDE_DEG
#define FIFO_DIAG_LONGITUDE_DEG (atoms3r_ical::CalibrationSite::longitude_deg)
#endif
#ifndef FIFO_DIAG_HEIGHT_M
#define FIFO_DIAG_HEIGHT_M (atoms3r_ical::CalibrationSite::ellipsoidal_height_m)
#endif
// Decimal year for the field model; defaults to the build date.
#ifndef FIFO_DIAG_YEAR
#define FIFO_DIAG_YEAR (imu_cal::decimalYearFromBuildDate(__DATE__))
#endif
// Delay fits wider than this (slow, smooth rotation) are reported but not
// counted toward the delay criterion.
#ifndef FIFO_DIAG_DELAY_FIT_WIDTH_MS
#define FIFO_DIAG_DELAY_FIT_WIDTH_MS 4.0
#endif
// Accuracy criteria (diagnostic limits, not specifications).
#ifndef FIFO_DIAG_ACCEL_ERROR_MG
#define FIFO_DIAG_ACCEL_ERROR_MG 10.0
#endif
#ifndef FIFO_DIAG_ACCEL_NOISE_MG
#define FIFO_DIAG_ACCEL_NOISE_MG 5.0
#endif
#ifndef FIFO_DIAG_GYRO_BIAS_DPS
#define FIFO_DIAG_GYRO_BIAS_DPS 0.2
#endif
#ifndef FIFO_DIAG_GYRO_NOISE_DPS
#define FIFO_DIAG_GYRO_NOISE_DPS 0.2
#endif
#ifndef FIFO_DIAG_GYRO_SCALE_PCT
#define FIFO_DIAG_GYRO_SCALE_PCT 1.5
#endif
#ifndef FIFO_DIAG_FIELD_ERROR_PCT
#define FIFO_DIAG_FIELD_ERROR_PCT 5.0
#endif
#ifndef FIFO_DIAG_DIP_ERROR_DEG
#define FIFO_DIAG_DIP_ERROR_DEG 3.0
#endif
#ifndef FIFO_DIAG_MAG_NORM_STD_PCT
#define FIFO_DIAG_MAG_NORM_STD_PCT 2.0
#endif

using namespace atoms3r_ical;

static constexpr uint32_t LOOP_PERIOD_US = 5000;
static constexpr uint32_t WARMUP_MS = 3000;
static constexpr uint32_t UI_MS = 100, TIMING_MS = 10000, ACCURACY_MS = 10000, IDENTITY_MS = 5000,
                          REFERENCE_MS = 30000, STATE_MS = 1000;
static constexpr size_t FIFO_BUF_BYTES = Bmi270::FIFO_BYTES + 64;
// Average header-mode frame: accelerometer+gyro (13 bytes), plus 8 bytes of
// magnetometer data in every 8th frame.
static constexpr double FIFO_FRAME_BYTES = 13.0 + 8.0 / 8.0;
static constexpr size_t QUEUE_FRAMES = 768;
static constexpr uint32_t MAG_MAX_AGE_US = MagAcquisition::MAX_AGE_MS * 1000u;
static constexpr unsigned FIT_PAIRS = 50;
static constexpr float RAD2DEG = 57.29577951308232f;

// One key=value record, written whole or not at all.
struct Line {
  char buf[1024];
  size_t n = 0;
  bool overflow = false;
  Line& reset(char tag) {
    n = 0;
    overflow = false;
    buf[n++] = tag;
    return *this;
  }
  Line& printf(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    const int w = vsnprintf(buf + n, sizeof(buf) - n, fmt, ap);
    va_end(ap);
    if (w < 0 || size_t(w) >= sizeof(buf) - n) overflow = true;
    else n += size_t(w);
    return *this;
  }
  Line& kv(const char* k, double v, int prec) {
    return std::isfinite(v) ? printf(",%s=%.*f", k, prec, v) : printf(",%s=nan", k);
  }
  Line& kv(const char* k, unsigned long v) { return printf(",%s=%lu", k, v); }
  Line& kv(const char* k, bool v) { return printf(",%s=%d", k, v ? 1 : 0); }
  Line& verdict(const char* k, int v) { return printf(",%s=%s", k, v < 0 ? "na" : (v ? "1" : "0")); }
};

static uint32_t lines_dropped = 0;
static void emit(Line& l) {
  if (l.overflow || l.n + 1 >= sizeof(l.buf)) { ++lines_dropped; return; }
  l.buf[l.n++] = '\n';
  if (Serial.availableForWrite() < int(l.n)) { ++lines_dropped; return; }
  Serial.write(reinterpret_cast<const uint8_t*>(l.buf), l.n);
}

struct QueuedFrame {
  uint32_t seq, sample_us;
  float a[3], w[3], m[3];
  float temp_c;
};

// The timing window, restarted by 'r'. It opens at the first FIFO read after
// the restart, when the sensor clock has a reading.
struct Evaluation {
  bool open = false;
  uint32_t start_ms = 0;
  uint64_t start_host_us = 0, start_ticks = 0, start_lost = 0, start_skipped = 0;
  uint64_t frames = 0, aux = 0, aux_interval_bad = 0, drdy_fail = 0, dup = 0, comp_fail = 0, partial = 0;
  double spacing_min_us = 1e9, spacing_max_us = 0;
  uint32_t stalls = 0, stall_max_ms = 0;
};

static M5Ui ui;
static ImuCalStoreNvs store;
static ImuCalBlobV5 blob{};
static RuntimeCals cals{};
static bool have_blob = false, mag_cal_ok = false;
static Line out, r_line;
static m5::IMU_Base* bmi270 = nullptr;
static FifoAcquisitionReadback setup_rb{};
static FifoAcquisitionConfig fifo_cfg{};
static float acc_scale = 0, gyr_scale = 0;
static bool setup_ok = false;
static const char* setup_failure = "not run";

static uint8_t fifo_buf[FIFO_BUF_BYTES];
static SensorClockMap sensor_clock;
static FifoParseStats parse_stats;
static uint64_t ticks_total = 0;  // 200 Hz frame positions received, skipped included
static uint64_t frames_delivered = 0, warmup_ticks = 0;
static uint32_t dither_state = 12345;
static size_t fifo_len_max = 0, last_fifo_len = 0;
static uint32_t prev_frame_us = 0;
static bool have_prev_frame = false;

// Magnetometer state (from AUX frames).
static uint8_t last_aux[8]{};
static bool have_aux = false;
static uint64_t last_aux_half_ticks = 0;
static uint32_t mag_seq = 0, mag_us = 0;
static Vector3f mag_body = Vector3f::Constant(NAN), mag_cal = Vector3f::Constant(NAN);
static float temp_c = NAN;
static Vector3f last_a_cal = Vector3f::Constant(NAN);

static imu_cal::MagTimingCal timing_cal;
static unsigned timing_segment = 0;
static float fit_delay_ms[16];
static unsigned n_fits = 0, n_wide_fits = 0, n_fit_attempts = 0;
static const char* last_fit_reason = "rock the board";

// Accuracy.
static imu_cal::ImuAccuracy accuracy;
static double site_lat = NAN, site_lon = NAN, site_h_m = NAN, site_year = NAN;
static imu_cal::GeomagField expected{};
static imu_cal::RotationCheck last_check{};
static bool have_check = false;
struct Verdicts { int accel = -1, accel_noise = -1, gyro_bias = -1, gyro_noise = -1, gyro_scale = -1, field = -1,
                  mag_norm = -1; };
static Verdicts verdicts{};

static QueuedFrame queue_buf[QUEUE_FRAMES];
static size_t q_head = 0, q_count = 0, q_max = 0;
static uint32_t q_dropped = 0;
static bool raw_lines = FIFO_DIAG_RAW_LINES != 0;
static uint8_t gyr_conf_now = 0;
static unsigned events = 0;
static constexpr unsigned MAX_EVENTS = 20;

static Evaluation ev{};
static uint32_t last_ui_ms = 0, last_timing_ms = 0, last_accuracy_ms = 0, last_identity_ms = 0,
                last_reference_ms = 0, last_state_ms = 0, last_load_ms = 0, last_stall_ms = 0;
static unsigned load_count = 0, page = 0;
static bool pass_rate = false, pass_loss = false, pass_spacing = false, pass_mag = false, pass_delay = false;
static char cmd[64];
static size_t cmd_len = 0;

// ---------------------------------------------------------------- helpers --

static void waitMs(uint32_t ms) { delay(ms); }

static void resetEvaluation() {
  ev = Evaluation{};
  ev.start_ms = millis();
}

static void openEvaluation() {
  ev.open = true;
  ev.start_host_us = sensor_clock.hostUs();
  ev.start_ticks = sensor_clock.sensortimeTicks();
  ev.start_lost = sensor_clock.lostFrames();
  ev.start_skipped = parse_stats.skipped;
}

static bool warm() { return uint32_t(millis() - ev.start_ms) >= WARMUP_MS; }

static void displayLoad(uint32_t ms) {
  const uint32_t t0 = millis();
  bool odd = false;
  while (uint32_t(millis() - t0) < ms) {
    M5.Display.fillScreen(odd ? TFT_BLUE : TFT_DARKGREY);
    M5.Display.setCursor(4, 56);
    M5.Display.setTextColor(TFT_WHITE);
    M5.Display.printf("DISPLAY LOAD %lu ms", static_cast<unsigned long>(ms));
    odd = !odd;
  }
  last_stall_ms = uint32_t(millis() - t0);
  ++ev.stalls;
  ev.stall_max_ms = std::max(ev.stall_max_ms, last_stall_ms);
  last_ui_ms = 0;
}

static float delaySpreadMs() {
  if (n_fits < 2) return NAN;
  const unsigned n = std::min<unsigned>(n_fits, 16);
  const float* lo = std::min_element(fit_delay_ms, fit_delay_ms + n);
  const float* hi = std::max_element(fit_delay_ms, fit_delay_ms + n);
  return *hi - *lo;
}

static float delayMeanMs() {
  const unsigned n = std::min<unsigned>(n_fits, 16);
  if (!n) return NAN;
  float s = 0;
  for (unsigned i = 0; i < n; ++i) s += fit_delay_ms[i];
  return s / float(n);
}

static void setSite(double lat, double lon, double h_m, double year) {
  site_lat = lat;
  site_lon = lon;
  site_h_m = h_m;
  site_year = year;
  expected = imu_cal::Wmm2025::field(lat, lon, h_m, year);
  imu_cal::AccuracyReference ref;
  ref.gravity_mps2 = ImuCalCfg::g_cal_local;
  ref.g_std_mps2 = ImuCalCfg::g_std;
  ref.have_field = std::isfinite(expected.f_nT);
  ref.f_uT = expected.f_nT * 1e-3;
  ref.h_uT = expected.h_nT * 1e-3;
  ref.z_uT = expected.z_nT * 1e-3;
  ref.inclination_deg = expected.inclination_deg;
  accuracy.setReference(ref);
}

// --------------------------------------------------------------- frames --

static void onAux(const FifoFrame& f, uint64_t half_ticks, uint32_t aux_host_us) {
  ++ev.aux;
  // The first readout after setup may predate the first triggered conversion,
  // so data-ready is checked from the second readout on.
  if (have_aux && !bmm150DataReady(f.m)) ++ev.drdy_fail;
  if (have_aux && bmm150SamePayload(f.m, last_aux)) ++ev.dup;
  if (have_aux && half_ticks - last_aux_half_ticks != 16) ++ev.aux_interval_bad;  // 8 frames = 40 ms
  have_aux = true;
  last_aux_half_ticks = half_ticks;
  std::memcpy(last_aux, f.m, 8);

  float sensor[3], body[3];
  if (!bmm150Compensate(Bmm150Raw::fromData(f.m), atoms3rMagSource().trim, sensor)) {
    ++ev.comp_fail;
    return;
  }
  bmm150SensorToAtomS3RBody(sensor, body);
  mag_body = Vector3f(body[0], body[1], body[2]);
  mag_cal = mag_cal_ok ? cals.applyMag(mag_body) : mag_body;
  ++mag_seq;
  mag_us = aux_host_us;
  accuracy.addMag(mag_cal, mag_body);
}

static void stepTimingFit(uint32_t now_us, const Vector3f& w_cal) {
  const bool have = mag_seq != 0 && uint32_t(now_us - mag_us) <= MAG_MAX_AGE_US;
  timing_cal.add(now_us, w_cal, have, mag_seq, mag_us, mag_cal, timing_segment);
  if (timing_cal.pairs(timing_segment) < FIT_PAIRS) return;
  if (timing_segment == 0) {
    const auto q = timing_cal.segmentQuality(0);
    if (q.ok) {
      timing_segment = 1;
      timing_cal.resetSegment(1);
    } else if (timing_cal.pairs(0) >= 3 * FIT_PAIRS) {
      last_fit_reason = imu_cal::magTimingReason(q.reason);
      timing_cal.reset();
    }
    return;
  }
  const auto fit = timing_cal.fit();
  ++n_fit_attempts;
  last_fit_reason = imu_cal::magTimingReason(fit.reason);
  if (fit.ok && fit.fit_width_s * 1e3f <= FIFO_DIAG_DELAY_FIT_WIDTH_MS)
    fit_delay_ms[(n_fits++) % 16] = float(fit.delay_us) * 1e-3f;
  else if (fit.ok)
    ++n_wide_fits;
  out.reset('X');
  out.kv("ok", fit.ok).kv("delay_ms", float(fit.delay_us) * 1e-3f, 1)
     .kv("width_ms", fit.fit_width_s * 1e3f, 1).kv("allowance_ms", fit.residual_time_s * 1e3f, 1)
     .kv("train_pairs", static_cast<unsigned long>(fit.train_pairs))
     .kv("validation_pairs", static_cast<unsigned long>(fit.validation_pairs))
     .kv("improvement", fit.improvement, 3).kv("train_rad", fit.train_rms_rad, 5)
     .kv("validation_rad", fit.validation_rms_rad, 5);
  out.printf(",reason=%s", last_fit_reason);
  emit(out);
  timing_cal.reset();
  timing_segment = 0;
}

static void onAccGyr(const FifoFrame& f, uint64_t k) {
  const uint32_t t_us = sensor_clock.frameHostUs(k);
  ImuSample s{};
  s.a = bmi270AccelToBodyNed(f.a, acc_scale);
  s.w = bmi270GyroToBodyNed(f.g, gyr_scale);
  s.tempC = temp_c;
  // The calibration's own temperature validity applies, as in the INS sketches.
  const float tempC = s.tempC;
  const Vector3f a_cal = cals.applyAccel(s.a, tempC);
  const Vector3f w_cal = cals.applyGyro(s.w, tempC);
  last_a_cal = a_cal;

  const double dt_s = have_prev_frame ? double(int32_t(t_us - prev_frame_us)) * 1e-6 : 0.0;
  if (have_prev_frame && warm()) {
    const double sp = dt_s * 1e6;
    ev.spacing_min_us = std::min(ev.spacing_min_us, sp);
    ev.spacing_max_us = std::max(ev.spacing_max_us, sp);
  }
  prev_frame_us = t_us;
  have_prev_frame = true;
  ++frames_delivered;
  ++ev.frames;
  stepTimingFit(t_us, w_cal);
  accuracy.addFrame(a_cal, w_cal, s.a, s.w, dt_s, temp_c);

  if (!raw_lines) return;
  if (q_count == QUEUE_FRAMES) { ++q_dropped; return; }
  QueuedFrame& q = queue_buf[(q_head + q_count) % QUEUE_FRAMES];
  ++q_count;
  q_max = std::max(q_max, q_count);
  q.seq = uint32_t(frames_delivered - 1);
  q.sample_us = t_us;
  const bool fresh = mag_seq != 0 && uint32_t(t_us - mag_us) <= MAG_MAX_AGE_US;
  for (int i = 0; i < 3; ++i) {
    q.a[i] = a_cal[i];
    q.w[i] = w_cal[i];
    q.m[i] = fresh ? mag_cal[i] : NAN;
  }
  q.temp_c = temp_c;
}

// R lines leave the queue as fast as USB takes them; a display stall only
// delays them.
static void emitQueued() {
  for (int n = 0; n < 40 && q_count && Serial.availableForWrite() >= 256; ++n) {
    const QueuedFrame& q = queue_buf[q_head];
    r_line.reset('R');
    r_line.printf(",%lu,%lu,%.5f,%.5f,%.5f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f,%.2f", static_cast<unsigned long>(q.seq),
                  static_cast<unsigned long>(q.sample_us), q.a[0], q.a[1], q.a[2], q.w[0], q.w[1], q.w[2], q.m[0],
                  q.m[1], q.m[2], q.temp_c);
    emit(r_line);
    q_head = (q_head + 1) % QUEUE_FRAMES;
    --q_count;
  }
}

// One burst for the whole length: the BMI270 does not continue a frame
// across separate reads, so a split read loses frame alignment.
static bool readFifo(size_t len) { return bmi270->readRegister(Bmi270::FIFO_DATA, fifo_buf, len); }

static void pollFifo() {
  // Sensortime, temperature and FIFO length in one burst (0x18..0x25).
  uint8_t r[14];
  const uint32_t host_us = micros();
  if (!bmi270->readRegister(Bmi270::SENSORTIME_0, r, sizeof(r))) return;
  const uint32_t st24 = uint32_t(r[0]) | uint32_t(r[1]) << 8 | uint32_t(r[2]) << 16;
  const int16_t t_raw = int16_t(uint16_t(r[10] | uint16_t(r[11]) << 8));
  if (t_raw != int16_t(0x8000)) temp_c = 23.0f + float(t_raw) / 512.0f;
  const size_t len = size_t(r[12]) | size_t(r[13] & 0x3F) << 8;
  last_fifo_len = len;
  fifo_len_max = std::max(fifo_len_max, len);

  const bool complete = len <= FIFO_BUF_BYTES;
  const size_t n = std::min(len, FIFO_BUF_BYTES);
  if (n && !readFifo(n)) return;
  if (!complete) ++ev.partial;

  // Pass 1: frame positions in this read, so the clock map sees the count up
  // to the length it read with the sensortime.
  uint64_t ticks_here = 0;
  FifoParseStats scratch;
  parseBmi270Fifo(fifo_buf, n, scratch,
                  [&](const FifoFrame& f) { if (f.acc || f.gyr) ++ticks_here; },
                  [&](unsigned skipped) { ticks_here += skipped; });
  const uint64_t lost_before = sensor_clock.lostFrames();
  const uint32_t breaks_before = sensor_clock.gridBreaks();
  if (complete) sensor_clock.observe(host_us, st24, ticks_total + ticks_here);
  if ((sensor_clock.lostFrames() != lost_before || sensor_clock.gridBreaks() != breaks_before) && events < MAX_EVENTS) {
    // Frames the sensor clock says were written but never arrived, or a
    // frame count ahead of the clock: report the read that showed it.
    ++events;
    out.reset('E');
    out.printf(",event=%s", sensor_clock.gridBreaks() != breaks_before ? "grid_break" : "lost");
    out.kv("t_ms", static_cast<unsigned long>(millis()))
       .kv("lost_added", static_cast<unsigned long>(sensor_clock.lostFrames() - lost_before))
       .kv("lost", static_cast<unsigned long>(sensor_clock.lostFrames()))
       .kv("grid_breaks", static_cast<unsigned long>(sensor_clock.gridBreaks()))
       .kv("warmup", !ev.open).kv("frames_counted", static_cast<unsigned long>(ticks_total + ticks_here))
       .kv("sensortime", static_cast<unsigned long>(st24)).kv("fifo_bytes", static_cast<unsigned long>(len))
       .kv("acc_only", static_cast<unsigned long>(parse_stats.acc_only + scratch.acc_only))
       .kv("gyr_only", static_cast<unsigned long>(parse_stats.gyr_only + scratch.gyr_only))
       .kv("aux_alone", static_cast<unsigned long>(parse_stats.aux_only + scratch.aux_only))
       .kv("unknown", static_cast<unsigned long>(parse_stats.unknown + scratch.unknown))
       .kv("skipped", static_cast<unsigned long>(parse_stats.skipped + scratch.skipped));
    emit(out);
  }
  if (!sensor_clock.settled()) {
    // Until reads at different phases pin the frame grid, frame times would
    // still move: these frames are a counted warm-up, not delivered.
    ticks_total += ticks_here;
    warmup_ticks += ticks_here;
    return;
  }
  if (!ev.open) openEvaluation();

  // Pass 2: deliver.
  parseBmi270Fifo(fifo_buf, n, parse_stats,
                  [&](const FifoFrame& f) {
                    const uint64_t k = ticks_total;
                    if (f.aux) {
                      // In a combined frame the readout shares the frame's tick;
                      // alone it lies between the previous tick and this one.
                      const bool alone = !(f.acc || f.gyr);
                      const uint32_t half_us = uint32_t(std::lround(64.0 * sensor_clock.usPerTick()));
                      onAux(f, alone ? 2 * k - 1 : 2 * k, sensor_clock.frameHostUs(k) - (alone ? half_us : 0u));
                    }
                    if (f.acc || f.gyr) {
                      if (f.acc && f.gyr) onAccGyr(f, k);
                      ++ticks_total;
                    }
                  },
                  [&](unsigned skipped) { ticks_total += skipped; });
}

// ----------------------------------------------------------------- output --

static void identityLine() {
  out.reset('H');
  out.printf(",sketch=atomS3R_imu_fifo_diagnostics,format=1,build=%s", __DATE__);
  out.kv("calibration", have_blob).kv("acc_cal", cals.acc.ok).kv("acc_gravity_mismatch", cals.accel_gravity_mismatch)
     .kv("gyr_cal", cals.gyr.ok).kv("mag_cal", mag_cal_ok)
     .kv("g_cal_local", ImuCalCfg::g_cal_local, 7).kv("lines_dropped", static_cast<unsigned long>(lines_dropped));
  // What the saved calibration holds: the temperatures it was made at, the
  // gyro bias it subtracts, and the field strength the magnetometer fit
  // scales to (the strength the sensor measured where it was calibrated).
  out.kv("acc_T0_c", blob.accel_T0, 2).kv("gyr_T0_c", blob.gyro_T0, 2)
     .kv("gyr_b0x_dps", blob.gyro_b0[0] * RAD2DEG, 4).kv("gyr_b0y_dps", blob.gyro_b0[1] * RAD2DEG, 4)
     .kv("gyr_b0z_dps", blob.gyro_b0[2] * RAD2DEG, 4)
     .kv("gyr_kx_dps_per_c", blob.gyro_k[0] * RAD2DEG, 5).kv("gyr_ky_dps_per_c", blob.gyro_k[1] * RAD2DEG, 5)
     .kv("gyr_kz_dps_per_c", blob.gyro_k[2] * RAD2DEG, 5)
     .kv("mag_cal_field_uT", blob.mag_field_uT, 2);
  emit(out);
}

static void referenceLine() {
  out.reset('G');
  out.kv("lat", site_lat, 6).kv("lon", site_lon, 6).kv("height_m", site_h_m, 1).kv("year", site_year, 3);
  out.printf(",model=WMM-2025");
  out.kv("model_valid", expected.in_validity).kv("F_uT", expected.f_nT * 1e-3, 3).kv("H_uT", expected.h_nT * 1e-3, 3)
     .kv("Z_uT", expected.z_nT * 1e-3, 3).kv("dip_deg", expected.inclination_deg, 3)
     .kv("declination_deg", expected.declination_deg, 3).kv("g_mps2", ImuCalCfg::g_cal_local, 7);
  emit(out);
}

static void fifoSetupLine() {
  out.reset('F');
  out.kv("ok", setup_ok);
  out.printf(",failure=%s", setup_failure);
  out.printf(",acc_conf=0x%02x,acc_range=0x%02x,gyr_conf=0x%02x,gyr_range=0x%02x,aux_conf=0x%02x,"
             "fifo_downs=0x%02x,fifo_config_0=0x%02x,fifo_config_1=0x%02x,aux_if_conf=0x%02x,aux_rd_addr=0x%02x,"
             "aux_wr_addr=0x%02x,aux_wr_data=0x%02x,pwr_conf=0x%02x,pwr_ctrl=0x%02x,err_reg=0x%02x,bmm_mode_before=0x%02x",
             setup_rb.acc_conf, setup_rb.acc_range, setup_rb.gyr_conf, setup_rb.gyr_range, setup_rb.aux_conf,
             setup_rb.fifo_downs, setup_rb.fifo_config_0, setup_rb.fifo_config_1, setup_rb.aux_if_conf,
             setup_rb.aux_rd_addr, setup_rb.aux_wr_addr, setup_rb.aux_wr_data, setup_rb.pwr_conf, setup_rb.pwr_ctrl,
             setup_rb.err_reg, setup_rb.bmm_mode_before);
  const unsigned offset = (setup_rb.aux_conf >> 4) & 0x0F;
  out.printf(",gyr_conf_now=0x%02x", gyr_conf_now);
  out.kv("fifo_capacity_ms", double(Bmi270::FIFO_BYTES) / FIFO_FRAME_BYTES * sensor_clock.frameSpacingUs() * 1e-3, 0);
  out.kv("aux_offset_ms", offset ? float(offset) * 2.5f : 40.0f, 1)
     .kv("acc_g_per_lsb", acc_scale, 8).kv("gyr_dps_per_lsb", gyr_scale, 6);
  emit(out);
}

static void timingLine() {
  if (!ev.open) return;
  const double host_s = double(sensor_clock.hostUs() - ev.start_host_us) * 1e-6;
  const double sensor_s = double(sensor_clock.sensortimeTicks() - ev.start_ticks) / 25600.0;
  const uint64_t lost = (sensor_clock.lostFrames() - ev.start_lost) + (parse_stats.skipped - ev.start_skipped);
  const double rate = host_s > 0 ? double(ev.frames) / host_s : 0;
  const double sensor_rate = sensor_s > 0 ? double(ev.frames) / sensor_s : 0;
  const double aux_hz = sensor_s > 0 ? double(ev.aux) / sensor_s : 0;
  const bool long_enough = host_s >= 10.0;
  // Delivered as produced: 200 Hz on the sensor's own clock, nothing missing.
  // The host-clock rate differs by the sensor oscillator's error (clock_ppm).
  pass_rate = long_enough && lost == 0 && std::fabs(sensor_rate - 200.0) <= 0.5;
  pass_loss = long_enough && lost == 0 && ev.stalls > 0 && ev.partial == 0;
  pass_spacing = ev.spacing_max_us > 0 && ev.spacing_min_us >= 4900.0 && ev.spacing_max_us <= 5100.0;
  // Data-ready proves each readout is a new conversion (reading the data
  // clears it). Identical consecutive payloads are reported, not failed: a
  // still board can repeat a quantized reading.
  pass_mag = long_enough && std::fabs(aux_hz - 25.0) <= 0.05 && ev.aux_interval_bad == 0 && ev.drdy_fail == 0 &&
             ev.comp_fail == 0;
  const float spread = delaySpreadMs();
  pass_delay = n_fits >= 3 && spread <= 4.0f;

  out.reset('T');
  out.kv("window_s", host_s, 1).kv("frames", static_cast<unsigned long>(ev.frames))
     .kv("rate_hz", rate, 3).kv("sensor_rate_hz", sensor_rate, 3).kv("lost", static_cast<unsigned long>(lost))
     .kv("partial_reads", static_cast<unsigned long>(ev.partial))
     .kv("stalls", static_cast<unsigned long>(ev.stalls)).kv("stall_max_ms", static_cast<unsigned long>(ev.stall_max_ms))
     .kv("fifo_max_bytes", static_cast<unsigned long>(fifo_len_max))
     .kv("spacing_min_us", ev.spacing_min_us, 1).kv("spacing_max_us", ev.spacing_max_us, 1)
     .kv("clock_ppm", sensor_clock.clockPpm(), 0)
     .kv("aux_hz", aux_hz, 3).kv("aux_interval_bad", static_cast<unsigned long>(ev.aux_interval_bad))
     .kv("drdy_fail", static_cast<unsigned long>(ev.drdy_fail)).kv("mag_repeat", static_cast<unsigned long>(ev.dup))
     .kv("mag_comp_fail", static_cast<unsigned long>(ev.comp_fail))
     .kv("delay_fits", static_cast<unsigned long>(n_fits)).kv("delay_wide_fits", static_cast<unsigned long>(n_wide_fits))
     .kv("delay_attempts", static_cast<unsigned long>(n_fit_attempts))
     .kv("delay_ms", delayMeanMs(), 1).kv("delay_spread_ms", spread, 1)
     .kv("r_dropped", static_cast<unsigned long>(q_dropped)).kv("lines_dropped", static_cast<unsigned long>(lines_dropped))
     .kv("rate_ok", pass_rate).kv("no_loss_ok", pass_loss).kv("spacing_ok", pass_spacing)
     .kv("mag_ok", pass_mag).kv("delay_ok", pass_delay)
     .kv("all_ok", pass_rate && pass_loss && pass_spacing && pass_mag && pass_delay);
  emit(out);
}

static void poseLine(const imu_cal::StillPose& p) {
  out.reset('S');
  out.kv("pose", static_cast<unsigned long>(p.index)).kv("duration_s", p.duration_s, 2)
     .kv("distinct", p.distinct).kv("roll_deg", p.roll_deg, 2).kv("pitch_deg", p.pitch_deg, 2)
     .kv("temp_c", p.temp_c, 2)
     .kv("a_mps2", p.accel_norm_mps2, 5).kv("a_error_mg", p.accel_error_mg, 2)
     .kv("a_raw_error_mg", p.accel_raw_error_mg, 2).kv("a_noise_mg", p.accel_noise_mg, 2)
     .kv("ax", p.accel_mean.x(), 5).kv("ay", p.accel_mean.y(), 5).kv("az", p.accel_mean.z(), 5)
     .kv("gyro_bias_dps", p.gyro_bias_dps, 4).kv("gyro_raw_bias_dps", p.gyro_raw_bias_dps, 4)
     .kv("gx_dps", p.gyro_mean.x() * RAD2DEG, 4).kv("gy_dps", p.gyro_mean.y() * RAD2DEG, 4)
     .kv("gz_dps", p.gyro_mean.z() * RAD2DEG, 4).kv("gyro_noise_dps", p.gyro_noise_dps, 4)
     .kv("mag_samples", static_cast<unsigned long>(p.mag_samples))
     .kv("F_uT", p.f_uT, 3).kv("F_raw_uT", p.f_raw_uT, 3).kv("Z_uT", p.z_uT, 3).kv("H_uT", p.h_uT, 3)
     .kv("dip_deg", p.inclination_deg, 3).kv("mag_noise_uT", p.mag_noise_uT, 3)
     .kv("F_error_pct", p.f_error_pct, 2).kv("Z_error_uT", p.z_error_uT, 3).kv("dip_error_deg", p.inclination_error_deg, 3);
  emit(out);
}

static void rotationLine(const imu_cal::RotationCheck& c) {
  out.reset('W');
  out.kv("from_pose", static_cast<unsigned long>(c.from_pose)).kv("to_pose", static_cast<unsigned long>(c.to_pose))
     .kv("motion_s", c.motion_s, 2).kv("rotation_deg", c.rotation_deg, 2)
     .kv("tilt_change_deg", c.tilt_change_deg, 2).kv("field_change_deg", c.field_change_deg, 2)
     .kv("accel_error_deg", c.accel_error_deg, 3).kv("mag_error_deg", c.mag_error_deg, 3)
     .kv("gyro_scale_pct", c.gyro_scale_accel_pct, 3).kv("gyro_scale_mag_pct", c.gyro_scale_mag_pct, 3)
     .kv("accel_residual_deg", c.accel_residual_deg, 3).kv("mag_residual_deg", c.mag_residual_deg, 3);
  emit(out);
}

static int within(bool have, double v, double limit) { return have ? int(std::isfinite(v) && v <= limit) : -1; }

static void accuracyLine() {
  const imu_cal::AccuracySummary& s = accuracy.summary();
  const bool poses = s.poses > 0;
  const double acc_err = std::max(std::fabs(s.accel_error_mg_min), std::fabs(s.accel_error_mg_max));
  const double scale_err = std::max(std::fabs(s.gyro_scale_pct_min), std::fabs(s.gyro_scale_pct_max));
  const double field_err = std::max(std::fabs(s.f_error_pct_min), std::fabs(s.f_error_pct_max));
  const bool field = std::isfinite(s.f_error_pct_max) && std::isfinite(expected.f_nT);
  verdicts.accel = within(poses, acc_err, FIFO_DIAG_ACCEL_ERROR_MG);
  verdicts.accel_noise = within(poses, s.accel_noise_mg_max, FIFO_DIAG_ACCEL_NOISE_MG);
  verdicts.gyro_bias = within(poses, s.gyro_bias_dps_max, FIFO_DIAG_GYRO_BIAS_DPS);
  verdicts.gyro_noise = within(poses, s.gyro_noise_dps_max, FIFO_DIAG_GYRO_NOISE_DPS);
  verdicts.gyro_scale = within(s.scale_fits >= 2, scale_err, FIFO_DIAG_GYRO_SCALE_PCT);
  verdicts.field = field ? int(field_err <= FIFO_DIAG_FIELD_ERROR_PCT &&
                               s.inclination_error_deg_max <= FIFO_DIAG_DIP_ERROR_DEG)
                         : -1;
  verdicts.mag_norm = within(s.mag_faces >= 4, s.mag_norm_std_pct, FIFO_DIAG_MAG_NORM_STD_PCT);

  out.reset('A');
  out.kv("poses", static_cast<unsigned long>(s.poses)).kv("distinct_poses", static_cast<unsigned long>(s.distinct_poses))
     .kv("still", accuracy.still()).kv("still_s", accuracy.stillSeconds(), 1)
     .kv("a_error_mg_min", s.accel_error_mg_min, 2).kv("a_error_mg_max", s.accel_error_mg_max, 2)
     .kv("a_noise_mg_max", s.accel_noise_mg_max, 2)
     .kv("gyro_bias_dps_max", s.gyro_bias_dps_max, 4).kv("gyro_noise_dps_max", s.gyro_noise_dps_max, 4)
     .kv("rotation_checks", static_cast<unsigned long>(s.rotation_checks))
     .kv("gyro_scale_fits", static_cast<unsigned long>(s.scale_fits))
     .kv("gyro_scale_pct_min", s.gyro_scale_pct_min, 3).kv("gyro_scale_pct_max", s.gyro_scale_pct_max, 3)
     .kv("gyro_scale_mag_fits", static_cast<unsigned long>(s.mag_scale_fits))
     .kv("gyro_scale_mag_pct_min", s.gyro_scale_mag_pct_min, 3).kv("gyro_scale_mag_pct_max", s.gyro_scale_mag_pct_max, 3)
     .kv("F_error_pct_min", s.f_error_pct_min, 2).kv("F_error_pct_max", s.f_error_pct_max, 2)
     .kv("dip_deg_min", s.inclination_deg_min, 3).kv("dip_deg_max", s.inclination_deg_max, 3)
     .kv("dip_error_deg_max", s.inclination_error_deg_max, 3).kv("Z_error_uT_max", s.z_error_uT_max, 3)
     .kv("mag_noise_uT_max", s.mag_noise_uT_max, 3)
     .kv("m_samples", static_cast<unsigned long>(s.mag_samples)).kv("m_mean_uT", s.mag_norm_mean_uT, 3)
     .kv("m_std_pct", s.mag_norm_std_pct, 3).kv("m_min_uT", s.mag_norm_min_uT, 3).kv("m_max_uT", s.mag_norm_max_uT, 3)
     .kv("m_faces", static_cast<unsigned long>(s.mag_faces))
     .kv("temp_c_min", s.temp_c_min, 2).kv("temp_c_max", s.temp_c_max, 2)
     .verdict("accel_ok", verdicts.accel).verdict("accel_noise_ok", verdicts.accel_noise)
     .verdict("gyro_bias_ok", verdicts.gyro_bias).verdict("gyro_noise_ok", verdicts.gyro_noise)
     .verdict("gyro_scale_ok", verdicts.gyro_scale).verdict("field_ok", verdicts.field)
     .verdict("mag_norm_ok", verdicts.mag_norm);
  emit(out);
}

static void stateLine() {
  out.reset('D');
  out.kv("t_ms", static_cast<unsigned long>(millis()))
     .kv("fifo_bytes", static_cast<unsigned long>(last_fifo_len))
     .kv("fifo_max_bytes", static_cast<unsigned long>(fifo_len_max))
     .kv("frames", static_cast<unsigned long>(frames_delivered))
     .kv("lost", static_cast<unsigned long>(sensor_clock.lostFrames()))
     .kv("skipped", static_cast<unsigned long>(parse_stats.skipped))
     .kv("unknown", static_cast<unsigned long>(parse_stats.unknown))
     .kv("config_frames", static_cast<unsigned long>(parse_stats.config_frames))
     .kv("acc_only", static_cast<unsigned long>(parse_stats.acc_only))
     .kv("gyr_only", static_cast<unsigned long>(parse_stats.gyr_only))
     .kv("aux_frames", static_cast<unsigned long>(parse_stats.aux))
     .kv("aux_alone", static_cast<unsigned long>(parse_stats.aux_only))
     .kv("grid_breaks", static_cast<unsigned long>(sensor_clock.gridBreaks()))
     .kv("spacing_us", sensor_clock.frameSpacingUs(), 2).kv("clock_ppm", sensor_clock.clockPpm(), 0)
     .kv("t0_unc_us", sensor_clock.t0UncertaintyTicks() * Bmi270::TICK_US, 0)
     .kv("warmup_frames", static_cast<unsigned long>(warmup_ticks))
     .kv("still", accuracy.still()).kv("still_s", accuracy.stillSeconds(), 1)
     .kv("a_g", last_a_cal.norm() / ImuCalCfg::g_std, 4).kv("m_uT", mag_cal.norm(), 2).kv("temp_c", temp_c, 2)
     .kv("queue", static_cast<unsigned long>(q_count)).kv("queue_max", static_cast<unsigned long>(q_max))
     .kv("stall_ms", static_cast<unsigned long>(last_stall_ms));
  out.printf(",gyr_conf=0x%02x", gyr_conf_now);
  emit(out);
}

static void streamLines() {
  const uint32_t now_ms = millis();
  imu_cal::StillPose pose;
  while (accuracy.takePose(pose)) poseLine(pose);
  if (accuracy.takeRotationCheck(last_check)) {
    have_check = true;
    rotationLine(last_check);
  }
  if (uint32_t(now_ms - last_identity_ms) >= IDENTITY_MS) {
    last_identity_ms = now_ms;
    identityLine();
    fifoSetupLine();
  }
  if (uint32_t(now_ms - last_reference_ms) >= REFERENCE_MS) {
    last_reference_ms = now_ms;
    referenceLine();
  }
  if (uint32_t(now_ms - last_timing_ms) >= TIMING_MS) {
    last_timing_ms = now_ms;
    timingLine();
  }
  if (uint32_t(now_ms - last_accuracy_ms) >= ACCURACY_MS) {
    last_accuracy_ms = now_ms;
    accuracyLine();
  }
  if (uint32_t(now_ms - last_state_ms) >= STATE_MS) {
    last_state_ms = now_ms;
    stateLine();
  }
}

static const char* mark(int v) { return v < 0 ? "-- " : (v ? "OK " : "BAD"); }

static void drawUi() {
  const imu_cal::AccuracySummary& s = accuracy.summary();
  const imu_cal::StillPose& p = accuracy.lastPose();
  const bool pose = accuracy.haveLastPose();
  switch (page) {
    case 0:
      ui.title("TIME");
      M5.Display.printf("rate %s %.1fHz\n", pass_rate ? "OK " : "-- ",
                        ev.frames && warm()
                            ? double(ev.frames) / (double(sensor_clock.hostUs() - ev.start_host_us) * 1e-6)
                            : 0.0);
      M5.Display.printf("loss %s lost %lu\n", pass_loss ? "OK " : "-- ",
                        static_cast<unsigned long>(sensor_clock.lostFrames() + parse_stats.skipped));
      M5.Display.printf("dt   %s %.0f-%.0fus\n", pass_spacing ? "OK " : "-- ",
                        ev.spacing_max_us > 0 ? ev.spacing_min_us : 0.0, ev.spacing_max_us);
      M5.Display.printf("mag  %s drdy %lu\n", pass_mag ? "OK " : "-- ", static_cast<unsigned long>(ev.drdy_fail));
      M5.Display.printf("dly  %s %.1fms n%u\n", pass_delay ? "OK " : "-- ", double(delayMeanMs()), n_fits);
      M5.Display.printf("fifo %lu/%lu B\n", static_cast<unsigned long>(last_fifo_len),
                        static_cast<unsigned long>(fifo_len_max));
      M5.Display.printf("%s\n", last_fit_reason);
      break;
    case 1:
      ui.title("ACC/GYR");
      M5.Display.printf("%s %.1fs poses %lu/%lu\n", accuracy.still() ? "STILL" : "MOVING", accuracy.stillSeconds(),
                        static_cast<unsigned long>(s.poses), static_cast<unsigned long>(s.distinct_poses));
      M5.Display.printf("|a| %s %+.1fmg\n", mark(verdicts.accel), pose ? p.accel_error_mg : 0.0);
      M5.Display.printf("  raw %+.1f nz %.1fmg\n", pose ? p.accel_raw_error_mg : 0.0, pose ? p.accel_noise_mg : 0.0);
      M5.Display.printf("bias %s %.3fdps\n", mark(verdicts.gyro_bias), pose ? p.gyro_bias_dps : 0.0);
      M5.Display.printf("  raw %.3f nz %.3f\n", pose ? p.gyro_raw_bias_dps : 0.0, pose ? p.gyro_noise_dps : 0.0);
      M5.Display.printf("tilt %.1f %.1f deg\n", pose ? p.roll_deg : 0.0, pose ? p.pitch_deg : 0.0);
      M5.Display.printf("temp %.1fC\n", double(temp_c));
      break;
    case 2:
      ui.title("MAG");
      M5.Display.printf("F   %s %.1f/%.1f\n", mark(verdicts.field), pose ? p.f_uT : 0.0, expected.f_nT * 1e-3);
      M5.Display.printf("Z   %.1f/%.1f uT\n", pose ? p.z_uT : 0.0, expected.z_nT * 1e-3);
      M5.Display.printf("dip %.1f/%.1f deg\n", pose ? p.inclination_deg : 0.0, expected.inclination_deg);
      M5.Display.printf("|m| %s %.2f%% f%d\n", mark(verdicts.mag_norm), s.mag_norm_std_pct, s.mag_faces);
      M5.Display.printf("now %.1fuT\n", double(mag_cal.norm()));
      M5.Display.printf("%.2f %.2f\n", site_lat, site_lon);
      M5.Display.printf("cal %s\n", mag_cal_ok ? "applied" : "NONE");
      break;
    default:
      ui.title("ROT");
      M5.Display.printf("scale %s n%lu\n", mark(verdicts.gyro_scale), static_cast<unsigned long>(s.scale_fits));
      M5.Display.printf(" %+.2f..%+.2f%%\n", s.gyro_scale_pct_min, s.gyro_scale_pct_max);
      M5.Display.printf("mag fit n%lu\n", static_cast<unsigned long>(s.mag_scale_fits));
      M5.Display.printf(" %+.2f..%+.2f%%\n", s.gyro_scale_mag_pct_min, s.gyro_scale_mag_pct_max);
      if (have_check) {
        M5.Display.printf("last %.0fdeg\n", last_check.rotation_deg);
        M5.Display.printf(" err %.2f/%.2fdeg\n", last_check.accel_error_deg, last_check.mag_error_deg);
      }
      break;
  }
}

static void command(const char* c) {
  switch (c[0]) {
    case 's': {
      const long ms = std::strtol(c + 1, nullptr, 10);
      displayLoad(ms > 0 && ms <= 10000 ? uint32_t(ms) : uint32_t(FIFO_DIAG_HEAVY_LOAD_MS));
      break;
    }
    case 'r':
      resetEvaluation();
      accuracy.reset();
      have_check = false;
      break;
    case 'q': raw_lines = !raw_lines; break;
    case 'n': {
      // Toggle the gyro noise-performance mode; compare the still-pose noise
      // (S lines) before and after.
      uint8_t v = uint8_t(gyr_conf_now ^ 0x40), rb = 0;
      if (bmi270 && bmi270->writeRegister8(Bmi270::GYR_CONF, v) && bmi270->readRegister(Bmi270::GYR_CONF, &rb, 1))
        gyr_conf_now = rb;
      resetEvaluation();
      accuracy.reset();
      have_check = false;
      fifoSetupLine();
      break;
    }
    case 'L': {
      double lat = NAN, lon = NAN, h = 0.0;
      if (std::sscanf(c + 1, "%lf %lf %lf", &lat, &lon, &h) >= 2 && std::fabs(lat) <= 90.0 && std::fabs(lon) <= 360.0) {
        setSite(lat, lon, h, site_year);
        accuracy.reset();  // poses so far were judged against the old place
        referenceLine();
      }
      break;
    }
    default: break;
  }
}

static void serialCommands() {
  while (Serial.available() > 0) {
    const int ch = Serial.read();
    if (ch < 0) break;
    if (ch == '\n' || ch == '\r') {
      cmd[cmd_len] = '\0';
      if (cmd_len) command(cmd);
      cmd_len = 0;
    } else if (cmd_len + 1 < sizeof(cmd)) {
      cmd[cmd_len++] = char(ch);
    }
  }
}

// ------------------------------------------------------------------- main --

void setup() {
  delay(200);
  Serial.setTxBufferSize(16384);
  Serial.begin(115200);
  Serial.setTxTimeoutMs(0);

  auto cfg = M5.config();
  cfg.internal_imu = true;
  M5.begin(cfg);
  clearM5UnifiedImuCalibration();
  configureAtomS3RImu(Serial);
  ui.begin();

  have_blob = store.load(blob);
  if (!have_blob) std::memset(static_cast<void*>(&blob), 0, sizeof(blob));
  cals.rebuildFromBlob(blob);
  // FIFO readouts are always Bosch-compensated; a magnetometer fit made on
  // uncompensated values does not apply to them.
  mag_cal_ok = cals.mag.ok && activeMagSource() == MAG_SOURCE_BMM150_COMPENSATED;
  if (cals.mag.ok && !mag_cal_ok) {
    Serial.println("[FIFO] saved magnetometer calibration was fitted on uncompensated data: not applied; "
                   "recalibrate for the delay fit and the field checks");
  }
  if (!atoms3rMagSource().trim.valid) Serial.println("[FIFO] BMM150 trim not read: no magnetometer values");
  setSite(FIFO_DIAG_LATITUDE_DEG, FIFO_DIAG_LONGITUDE_DEG, FIFO_DIAG_HEIGHT_M, FIFO_DIAG_YEAR);

  bmi270 = M5.Imu.getImuInstancePtr(0);
  if (M5.Imu.getType() != m5::imu_bmi270 || !bmi270) {
    setup_failure = "no BMI270";
  } else {
    fifo_cfg.aux_offset = FIFO_DIAG_AUX_OFFSET;
    fifo_cfg.gyr_conf = FIFO_DIAG_GYR_CONF;
    Bmi270FifoSetup<m5::IMU_Base> setup_fifo(bmi270, waitMs);
    setup_ok = setup_fifo.apply(fifo_cfg);
    setup_failure = setup_fifo.failure();
    setup_rb = setup_fifo.readback();
  }
  gyr_conf_now = setup_rb.gyr_conf;
  acc_scale = bmi270AccelGPerLsb(setup_rb.acc_range);
  gyr_scale = bmi270GyroDpsPerLsb(setup_rb.gyr_range);

  identityLine();
  referenceLine();
  fifoSetupLine();
  const uint32_t now_ms = millis();
  last_identity_ms = last_reference_ms = last_timing_ms = last_accuracy_ms = last_state_ms = now_ms;
  Serial.printf("[FIFO] setup %s (%s)\n", setup_ok ? "ok" : "FAILED", setup_failure);
  if (!setup_ok) {
    ui.fail("FIFO", setup_failure);
    return;
  }
  bmi270->writeRegister8(Bmi270::CMD, Bmi270::CMD_FIFO_FLUSH);
  delay(2);
  resetEvaluation();
  last_load_ms = millis();
}

void loop() {
  const uint32_t loop_start_us = micros();
  if (!setup_ok) {
    delay(100);
    streamLines();
    return;
  }
  Input::update();
  if (Input::tapPressed()) {
    page = (page + 1) % 4;
    last_ui_ms = 0;
  }
  serialCommands();

  pollFifo();
  emitQueued();

  const uint32_t now_ms = millis();
  if (uint32_t(now_ms - last_load_ms) >= FIFO_DIAG_LOAD_EVERY_MS) {
    last_load_ms = now_ms;
    displayLoad((load_count++ & 1u) ? FIFO_DIAG_HEAVY_LOAD_MS : FIFO_DIAG_LOAD_MS);
  }
  if (uint32_t(now_ms - last_ui_ms) >= UI_MS) {
    last_ui_ms = now_ms;
    drawUi();
  }
  streamLines();

  const uint32_t busy_us = micros() - loop_start_us;
  // A loop period dithered by +-0.5 ms keeps FIFO reads from locking to the
  // frame phase, so every phase is sampled and the frame grid gets pinned.
  dither_state = dither_state * 1664525u + 1013904223u;
  const uint32_t period_us = LOOP_PERIOD_US - 500u + (dither_state >> 22);  // +0..1023 us
  if (busy_us < period_us) delayMicroseconds(period_us - busy_us);
}
