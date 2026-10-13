/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R IMU FIFO timing diagnostics: the acquisition the INS sketches should
  move to, measured on the device before any estimator depends on it.

  It takes the BMI270 over from M5Unified's polling path
  (src/AtomS3R/AtomS3R_Bmi270Fifo.h):
    - accelerometer and gyro both at 200 Hz (+-8 g, +-2000 deg/s, the ranges
      M5Unified's conversion assumes, so samples are in readImuMapped() units);
    - every sample and every magnetometer readout through the BMI270 FIFO,
      drained each loop, so nothing is lost while the display is busy;
    - each sample timed by its place on the sensor's own clock (a frame every
      128 sensortime ticks), mapped onto host micros();
    - the BMM150 in forced mode, triggered by the BMI270 on every 25 Hz AUX
      cycle, so each readout is one complete conversion with a fixed relation
      to the sensor clock.

  It loads the saved calibration but never changes it, and estimates nothing.

  Pass criteria, evaluated continuously and printed as P lines every 10 s:
    rate      200 Hz delivered (no frame missing on the sensor clock; host-clock
              rate within 1 %)
    no loss   no lost frame while the display is loaded: the UI redraws at
              10 Hz, a 400 ms display burst runs every 10 s, and a tap (or 's'
              on serial) runs a 1.5 s burst
    spacing   delivered sample spacing between 4.9 and 5.1 ms
    mag       25 Hz readouts, each a new conversion (data-ready set), every
              one 40 ms after the previous on the sensor clock
    delay     the magnetometer delay measured from hand rotation (rock the
              board about one rock per second) by the calibration wizard's
              fit, three fits agreeing within 4 ms

  USB serial (115200) carries the AtomS3R_InsDiag.h line stream: R lines for
  every frame (replayable through the INS sketches, see
  docs/ins-device-diagnostics.md), D loop/acquisition state at 10 Hz, H/I/C
  identity, registers and calibration, F the FIFO setup read-back, X one line
  per delay fit, P the pass criteria. Capture it with tools/ins_diag_capture.py
  and summarise it with tools/ins_diag_report.py.

  Serial commands: 's' 1.5 s display burst, 'r' restart the evaluation,
  'q' stop/start the R lines.
*/

#include <Arduino.h>
#include <M5Unified.h>
#include <ArduinoOceanImu.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"
#include "AtomS3R/AtomS3R_InsDiag.h"
#include "AtomS3R/AtomS3R_Bmi270Fifo.h"
#include "imu_calibrate/MagTimingCal.h"

// Trigger-to-readout offset of the AUX interface in 2.5 ms units; 0 keeps a
// full AUX period between the forced-mode command and the readout.
#ifndef FIFO_TIMING_AUX_OFFSET
#define FIFO_TIMING_AUX_OFFSET 0
#endif
#ifndef FIFO_TIMING_LOAD_EVERY_MS
#define FIFO_TIMING_LOAD_EVERY_MS 10000
#endif
#ifndef FIFO_TIMING_LOAD_MS
#define FIFO_TIMING_LOAD_MS 400
#endif
#ifndef FIFO_TIMING_HEAVY_LOAD_MS
#define FIFO_TIMING_HEAVY_LOAD_MS 1500
#endif

using namespace atoms3r_ical;

static constexpr float LOOP_HZ = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = 5000;
static constexpr uint32_t WARMUP_MS = 3000;
static constexpr uint32_t UI_MS = 100, PASS_MS = 10000, FIFO_LINE_MS = 5000;
static constexpr size_t FIFO_BUF_BYTES = 6144 + 256;
static constexpr size_t I2C_CHUNK = 128;
static constexpr size_t QUEUE_FRAMES = 768;
static constexpr uint32_t MAG_MAX_AGE_US = MagAcquisition::MAX_AGE_MS * 1000u;
static constexpr unsigned FIT_PAIRS = 50;

struct QueuedFrame {
  uint32_t seq, sample_us, mag_seq, mag_us, mag_ms;
  int16_t a[3], g[3];
  float temp_c;
  float m[3];
  bool mag_have;
};

// The pass-criteria window, restarted by 'r'. It opens at the first FIFO read
// after the restart, when the sensor clock has a reading.
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
static InsDiagStream<ImuCalBlobV5> diag{"atomS3R_imu_fifo_timing", LOOP_HZ};
static DiagLine r_line, x_line;
static m5::IMU_Base* bmi270 = nullptr;
static FifoAcquisitionReadback setup_rb{};
static FifoAcquisitionConfig fifo_cfg{};
static float acc_scale = 0, gyr_scale = 0;
static bool setup_ok = false;
static const char* setup_failure = "not run";

static uint8_t fifo_buf[FIFO_BUF_BYTES];
static SensorClockMap sensor_clock;
static FifoParseStats parse_stats;
static uint64_t ticks_total = 0;      // 200 Hz frame positions received, skipped included
static uint64_t frames_delivered = 0, warmup_ticks = 0;
static uint32_t dither_state = 12345;
static size_t fifo_len_max = 0, last_fifo_len = 0;
static uint32_t last_read_host_us = 0;
static uint32_t prev_frame_us = 0;
static bool have_prev_frame = false;

// Magnetometer state (from AUX frames).
static uint8_t last_aux[8]{};
static bool have_aux = false;
static uint64_t last_aux_half_ticks = 0;
static uint32_t mag_seq = 0, mag_us = 0, mag_ms = 0;
static Vector3f mag_body = Vector3f::Constant(NAN), mag_cal = Vector3f::Constant(NAN);
static float temp_c = NAN;

static imu_cal::MagTimingCal timing_cal;
static unsigned timing_segment = 0;
static float fit_delay_ms[16];
static unsigned n_fits = 0, n_fit_attempts = 0;
static const char* last_fit_reason = "rock the board";

static QueuedFrame queue_buf[QUEUE_FRAMES];
static size_t q_head = 0, q_count = 0, q_max = 0;
static uint32_t q_dropped = 0;
static bool raw_lines = true;

static Evaluation ev{};
static uint32_t last_ui_ms = 0, last_pass_ms = 0, last_fifo_line_ms = 0, last_load_ms = 0, last_stall_ms = 0;
static bool pass_rate = false, pass_loss = false, pass_spacing = false, pass_mag = false, pass_delay = false;

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

// --------------------------------------------------------------- frames --

static void onAux(const FifoFrame& f, uint64_t half_ticks, uint32_t aux_host_us) {
  ++ev.aux;
  // The first readout after setup may predate the first triggered conversion
  // (the register read-back pauses the AUX interface), so data-ready is
  // checked from the second readout on.
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
  mag_ms = millis();
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
  if (fit.ok) fit_delay_ms[(n_fits++) % 16] = float(fit.delay_us) * 1e-3f;
  x_line.reset('X');
  x_line.kv("ok", fit.ok).kv("delay_ms", float(fit.delay_us) * 1e-3f, 1)
        .kv("width_ms", fit.fit_width_s * 1e3f, 1).kv("allowance_ms", fit.residual_time_s * 1e3f, 1)
        .kv("train_pairs", static_cast<unsigned long>(fit.train_pairs))
        .kv("validation_pairs", static_cast<unsigned long>(fit.validation_pairs))
        .kv("improvement", fit.improvement, 3).kv("train_rad", fit.train_rms_rad, 5)
        .kv("validation_rad", fit.validation_rms_rad, 5);
  x_line.printf(",reason=%s", last_fit_reason);
  diag.emit(Serial, x_line);
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

  if (have_prev_frame && warm()) {
    const double sp = double(int32_t(t_us - prev_frame_us));
    ev.spacing_min_us = std::min(ev.spacing_min_us, sp);
    ev.spacing_max_us = std::max(ev.spacing_max_us, sp);
  }
  prev_frame_us = t_us;
  have_prev_frame = true;
  ++frames_delivered;
  ++ev.frames;
  diag.health().frame(t_us, w_cal, a_cal, ImuCalCfg::g_cal_local);
  stepTimingFit(t_us, w_cal);

  if (!raw_lines) return;
  if (q_count == QUEUE_FRAMES) { ++q_dropped; return; }
  QueuedFrame& q = queue_buf[(q_head + q_count) % QUEUE_FRAMES];
  ++q_count;
  q_max = std::max(q_max, q_count);
  q.seq = uint32_t(frames_delivered - 1);
  q.sample_us = t_us;
  std::memcpy(q.a, f.a, sizeof(q.a));
  std::memcpy(q.g, f.g, sizeof(q.g));
  q.temp_c = temp_c;
  q.mag_have = mag_seq != 0;
  q.mag_seq = mag_seq;
  q.mag_us = mag_us;
  q.mag_ms = mag_ms;
  const bool fresh = mag_seq != 0 && uint32_t(t_us - mag_us) <= MAG_MAX_AGE_US;
  for (int i = 0; i < 3; ++i) q.m[i] = fresh ? mag_body[i] : NAN;
}

// R lines leave the queue as fast as USB takes them; a display stall only
// delays them.
static void emitQueued() {
  for (int n = 0; n < 40 && q_count && Serial.availableForWrite() >= 256; ++n) {
    const QueuedFrame& q = queue_buf[q_head];
    ImuSample s{};
    s.a = bmi270AccelToBodyNed(q.a, acc_scale);
    s.w = bmi270GyroToBodyNed(q.g, gyr_scale);
    s.m = Vector3f(q.m[0], q.m[1], q.m[2]);
    s.tempC = q.temp_c;
    s.mask = 3;
    s.sample_us = q.sample_us;
    s.mag_timing.have = q.mag_have;
    s.mag_timing.last_ms = q.mag_ms;
    s.mag_timing.sequence = q.mag_seq;
    s.mag_timing.frame_us = s.mag_timing.read_start_us = s.mag_timing.read_end_us = q.mag_us;
    rawFrameLine(r_line, q.seq, s);
    diag.emit(Serial, r_line);
    q_head = (q_head + 1) % QUEUE_FRAMES;
    --q_count;
  }
}

static bool readFifo(size_t len) {
  for (size_t off = 0; off < len; off += I2C_CHUNK) {
    if (!bmi270->readRegister(Bmi270::FIFO_DATA, fifo_buf + off, std::min(I2C_CHUNK, len - off))) return false;
  }
  return true;
}

static void pollFifo() {
  // Sensortime, temperature and FIFO length in one burst (0x18..0x25).
  uint8_t r[14];
  const uint32_t host_us = micros();
  if (!bmi270->readRegister(Bmi270::SENSORTIME_0, r, sizeof(r))) return;
  last_read_host_us = host_us;
  const uint32_t st24 = uint32_t(r[0]) | uint32_t(r[1]) << 8 | uint32_t(r[2]) << 16;
  const int16_t t_raw = int16_t(uint16_t(r[10] | uint16_t(r[11]) << 8));
  if (t_raw != int16_t(0x8000)) temp_c = 23.0f + float(t_raw) / 512.0f;
  const size_t len = size_t(r[12]) | size_t(r[13] & 0x3F) << 8;
  last_fifo_len = len;
  fifo_len_max = std::max(fifo_len_max, len);
  diag.health().poll(3, 1, 2);

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
  if (complete) sensor_clock.observe(host_us, st24, ticks_total + ticks_here);
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

static void evaluate() {
  if (!ev.open) return;
  const double host_s = double(sensor_clock.hostUs() - ev.start_host_us) * 1e-6;
  const double sensor_s = double(sensor_clock.sensortimeTicks() - ev.start_ticks) / 25600.0;
  const uint64_t lost = (sensor_clock.lostFrames() - ev.start_lost) + (parse_stats.skipped - ev.start_skipped);
  const double rate = host_s > 0 ? double(ev.frames) / host_s : 0;
  const double aux_hz = sensor_s > 0 ? double(ev.aux) / sensor_s : 0;
  const bool long_enough = host_s >= 10.0;
  pass_rate = long_enough && lost == 0 && std::fabs(rate - 200.0) <= 2.0;
  pass_loss = long_enough && lost == 0 && ev.stalls > 0 && ev.partial == 0;
  pass_spacing = ev.spacing_max_us > 0 && ev.spacing_min_us >= 4900.0 && ev.spacing_max_us <= 5100.0;
  // Data-ready proves each readout is a new conversion (reading the data
  // clears it). Identical consecutive payloads are reported, not failed: a
  // still board can repeat a quantized reading.
  pass_mag = long_enough && std::fabs(aux_hz - 25.0) <= 0.05 && ev.aux_interval_bad == 0 && ev.drdy_fail == 0 &&
             ev.comp_fail == 0;
  const float spread = delaySpreadMs();
  pass_delay = n_fits >= 3 && spread <= 4.0f;

  DiagLine& p = x_line.reset('P');
  p.kv("window_s", float(host_s), 1).kv("frames", static_cast<unsigned long>(ev.frames))
   .kv("rate_hz", float(rate), 3).kv("lost", static_cast<unsigned long>(lost))
   .kv("partial_reads", static_cast<unsigned long>(ev.partial))
   .kv("stalls", static_cast<unsigned long>(ev.stalls)).kv("stall_max_ms", static_cast<unsigned long>(ev.stall_max_ms))
   .kv("fifo_max_bytes", static_cast<unsigned long>(fifo_len_max))
   .kv("spacing_min_us", float(ev.spacing_min_us), 1).kv("spacing_max_us", float(ev.spacing_max_us), 1)
   .kv("clock_ppm", float(sensor_clock.clockPpm()), 0)
   .kv("aux_hz", float(aux_hz), 3).kv("aux_interval_bad", static_cast<unsigned long>(ev.aux_interval_bad))
   .kv("drdy_fail", static_cast<unsigned long>(ev.drdy_fail)).kv("mag_repeat", static_cast<unsigned long>(ev.dup))
   .kv("mag_comp_fail", static_cast<unsigned long>(ev.comp_fail))
   .kv("delay_fits", n_fits).kv("delay_attempts", n_fit_attempts)
   .kv("delay_ms", delayMeanMs(), 1).kv("delay_spread_ms", spread, 1)
   .kv("rate_ok", pass_rate).kv("no_loss_ok", pass_loss).kv("spacing_ok", pass_spacing)
   .kv("mag_ok", pass_mag).kv("delay_ok", pass_delay)
   .kv("all_ok", pass_rate && pass_loss && pass_spacing && pass_mag && pass_delay);
  diag.emit(Serial, p);
}

static void fifoSetupLine() {
  DiagLine& f = x_line.reset('F');
  f.kv("ok", setup_ok);
  f.printf(",failure=%s", setup_failure);
  f.printf(",acc_conf=0x%02x,acc_range=0x%02x,gyr_conf=0x%02x,gyr_range=0x%02x,aux_conf=0x%02x,"
           "fifo_downs=0x%02x,fifo_config_0=0x%02x,fifo_config_1=0x%02x,aux_if_conf=0x%02x,aux_rd_addr=0x%02x,"
           "aux_wr_addr=0x%02x,aux_wr_data=0x%02x,pwr_conf=0x%02x,pwr_ctrl=0x%02x,err_reg=0x%02x,bmm_mode_before=0x%02x",
           setup_rb.acc_conf, setup_rb.acc_range, setup_rb.gyr_conf, setup_rb.gyr_range, setup_rb.aux_conf,
           setup_rb.fifo_downs, setup_rb.fifo_config_0, setup_rb.fifo_config_1, setup_rb.aux_if_conf,
           setup_rb.aux_rd_addr, setup_rb.aux_wr_addr, setup_rb.aux_wr_data, setup_rb.pwr_conf, setup_rb.pwr_ctrl,
           setup_rb.err_reg, setup_rb.bmm_mode_before);
  const unsigned offset = (setup_rb.aux_conf >> 4) & 0x0F;
  f.kv("aux_offset_ms", offset ? float(offset) * 2.5f : 40.0f, 1)
   .kv("acc_g_per_lsb", acc_scale, 8).kv("gyr_dps_per_lsb", gyr_scale, 6)
   .kv("mag_cal", mag_cal_ok).kv("calibration", have_blob);
  diag.emit(Serial, f);
}

static void streamDiag() {
  const uint32_t now_ms = millis();
  diag.periodic(Serial, now_ms, [](DiagLine& m) {
    m.printf("[FIFO MAG] segment=%u pairs=%u gaps=%u history_miss=%u fits=%u attempts=%u last=%s",
             timing_segment, timing_cal.pairs(timing_segment), timing_cal.gaps(), timing_cal.missingHistory(),
             n_fits, n_fit_attempts, last_fit_reason);
  });
  if (uint32_t(now_ms - last_fifo_line_ms) >= FIFO_LINE_MS) {
    last_fifo_line_ms = now_ms;
    fifoSetupLine();
  }
  if (uint32_t(now_ms - last_pass_ms) >= PASS_MS) {
    last_pass_ms = now_ms;
    evaluate();
  }
  if (!diag.stateDue(now_ms)) return;
  auto& d = diag.beginState(now_ms);
  d.kv("fifo_bytes", static_cast<unsigned long>(last_fifo_len))
   .kv("fifo_max_bytes", static_cast<unsigned long>(fifo_len_max))
   .kv("ticks", static_cast<unsigned long>(ticks_total))
   .kv("lost", static_cast<unsigned long>(sensor_clock.lostFrames()))
   .kv("skipped", static_cast<unsigned long>(parse_stats.skipped))
   .kv("acc_only", static_cast<unsigned long>(parse_stats.acc_only))
   .kv("gyr_only", static_cast<unsigned long>(parse_stats.gyr_only))
   .kv("aux_frames", static_cast<unsigned long>(parse_stats.aux))
   .kv("aux_alone", static_cast<unsigned long>(parse_stats.aux_only))
   .kv("unknown", static_cast<unsigned long>(parse_stats.unknown))
   .kv("spacing_us", float(sensor_clock.frameSpacingUs()), 2).kv("clock_ppm", float(sensor_clock.clockPpm()), 0)
   .kv("t0_unc_us", float(sensor_clock.t0UncertaintyTicks() * Bmi270::TICK_US), 0)
   .kv("warmup_frames", static_cast<unsigned long>(warmup_ticks))
   .kv("grid_breaks", static_cast<unsigned long>(sensor_clock.gridBreaks()))
   .kv("drdy_fail", static_cast<unsigned long>(ev.drdy_fail)).kv("mag_repeat", static_cast<unsigned long>(ev.dup))
   .kv("aux_interval_bad", static_cast<unsigned long>(ev.aux_interval_bad))
   .kv("mag_uT", mag_body.norm(), 2).kv("temp_c", temp_c, 2)
   .kv("stall_ms", static_cast<unsigned long>(last_stall_ms))
   .kv("queue", static_cast<unsigned long>(q_count)).kv("queue_max", static_cast<unsigned long>(q_max))
   .kv("r_dropped", static_cast<unsigned long>(q_dropped))
   .kv("rate_ok", pass_rate).kv("no_loss_ok", pass_loss).kv("spacing_ok", pass_spacing)
   .kv("mag_ok", pass_mag).kv("delay_ok", pass_delay);
  diag.emitState(Serial, d);
}

static void drawUi() {
  ui.title("FIFO TIME");
  M5.Display.printf("rate %s %.1fHz\n", pass_rate ? "OK " : "-- ",
                    ev.frames && warm() ? double(ev.frames) / (double(sensor_clock.hostUs() - ev.start_host_us) * 1e-6)
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
  M5.Display.printf("tap: 1.5 s load\n");
}

static void serialCommands() {
  while (Serial.available() > 0) {
    switch (Serial.read()) {
      case 's': displayLoad(FIFO_TIMING_HEAVY_LOAD_MS); break;
      case 'r': resetEvaluation(); break;
      case 'q': raw_lines = !raw_lines; break;
      default: break;
    }
  }
}

// ------------------------------------------------------------------- main --

void setup() {
  delay(200);
  Serial.setTxBufferSize(16384);
  Serial.begin(115200);
  Serial.setTxTimeoutMs(0);
  InsDiagStream<ImuCalBlobV5>::waitForReader(Serial, 5000);

  auto cfg = M5.config();
  cfg.internal_imu = true;
  M5.begin(cfg);
  clearM5UnifiedImuCalibration();
  configureAtomS3RMag(Serial);
  ui.begin();

  have_blob = store.load(blob);
  if (!have_blob) std::memset(static_cast<void*>(&blob), 0, sizeof(blob));
  cals.rebuildFromBlob(blob);
  // FIFO readouts are always Bosch-compensated; a magnetometer fit made on
  // uncompensated values does not apply to them.
  mag_cal_ok = cals.mag.ok && activeMagSource() == MAG_SOURCE_BMM150_COMPENSATED;
  if (cals.mag.ok && !mag_cal_ok) {
    Serial.println("[FIFO] saved magnetometer calibration was fitted on uncompensated data: not applied; "
                   "recalibrate for the delay fit");
  }
  if (!atoms3rMagSource().trim.valid) Serial.println("[FIFO] BMM150 trim not read: no magnetometer values");

  bmi270 = M5.Imu.getImuInstancePtr(0);
  if (M5.Imu.getType() != m5::imu_bmi270 || !bmi270) {
    setup_failure = "no BMI270";
  } else {
    fifo_cfg.aux_offset = FIFO_TIMING_AUX_OFFSET;
    Bmi270FifoSetup<m5::IMU_Base> setup_fifo(bmi270, waitMs);
    setup_ok = setup_fifo.apply(fifo_cfg);
    setup_failure = setup_fifo.failure();
    setup_rb = setup_fifo.readback();
  }
  acc_scale = bmi270AccelGPerLsb(setup_rb.acc_range);
  gyr_scale = bmi270GyroDpsPerLsb(setup_rb.gyr_range);

  diag.setHeaderExtra("acquisition=fifo");
  // Register read-back (I line) goes through the AUX interface in manual
  // mode; the FIFO is flushed again afterwards.
  const ImuTimingRegisters regs = setup_ok ? readImuTimingRegisters(true) : ImuTimingRegisters{};
  activeMagSource() = MAG_SOURCE_BMM150_COMPENSATED;
  diag.boot(regs, have_blob, MAG_SOURCE_BMM150_COMPENSATED, blob);
  fifoSetupLine();
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
    streamDiag();
    return;
  }
  Input::update();
  if (Input::tapPressed()) displayLoad(FIFO_TIMING_HEAVY_LOAD_MS);
  serialCommands();

  pollFifo();
  emitQueued();

  const uint32_t now_ms = millis();
  if (uint32_t(now_ms - last_load_ms) >= FIFO_TIMING_LOAD_EVERY_MS) {
    last_load_ms = now_ms;
    displayLoad(FIFO_TIMING_LOAD_MS);
  }
  if (uint32_t(now_ms - last_ui_ms) >= UI_MS) {
    last_ui_ms = now_ms;
    drawUi();
  }
  streamDiag();

  const uint32_t busy_us = micros() - loop_start_us;
  diag.tick(busy_us, LOOP_PERIOD_US);
  // A loop period dithered by +-0.5 ms keeps FIFO reads from locking to the
  // frame phase, so every phase is sampled and the frame grid gets pinned.
  dither_state = dither_state * 1664525u + 1013904223u;
  const uint32_t period_us = LOOP_PERIOD_US - 500u + (dither_state >> 22);  // +0..1023 us
  if (busy_us < period_us) delayMicroseconds(period_us - busy_us);
}
