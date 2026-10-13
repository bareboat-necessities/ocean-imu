#pragma once

// Copyright 2026, Mikhail Grushinskiy
//
// Opt-in serial diagnostics for the AtomS3R marine INS sketches
// (SEA_STATE_SERIAL_DIAG=1; host-testable). The default build is unaffected.
//
// Every line is formatted in full and written only when the port can take
// all of it, so a slow or absent reader drops whole lines, never parts:
//
//   H  sketch identity and line format, at boot and every few seconds
//   I  BMI270/BMM150 timing registers as read back from the sensor
//   C  the calibration blob the sketch loaded (exact bytes) and the
//      magnetometer source it selected
//   R  one accepted IMU frame exactly as readImuMapped() delivered it:
//      floats as IEEE-754 bit patterns, plus the magnetometer acquisition
//      record. R lines carry a frame sequence number, so a capture with lost
//      lines is detectable, and they are what the host replay feeds back into
//      the unchanged sketch (tests/common/InsSketchReplay.cpp).
//   D  ten times a second: loop/acquisition health and estimator internals
//      as key=value pairs
//   M  once a second: the shared magnetic-input timing line

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>

#include <Arduino.h>

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

namespace atoms3r_ical {

// Bumped whenever an R/C field is added, removed or reinterpreted.
constexpr unsigned kInsDiagFrameFormat = 1;

inline uint32_t floatBits(float f) {
  uint32_t u = 0;
  std::memcpy(&u, &f, sizeof(u));
  return u;
}
inline float floatFromBits(uint32_t u) {
  float f = 0.0f;
  std::memcpy(&f, &u, sizeof(f));
  return f;
}

// One diagnostic line, written whole or not at all. The 1 kB buffer is meant
// to live in the application object, not on the 8 kB loop-task stack:
// InsDiagStream keeps one and reuses it for every line.
class DiagLine {
 public:
  DiagLine() { reset('?'); }
  explicit DiagLine(char tag) { reset(tag); }

  DiagLine& reset(char tag) {
    buf_[0] = tag;
    buf_[1] = '\0';
    n_ = 1;
    overflow_ = false;
    return *this;
  }

  int printf(const char* fmt, ...) __attribute__((format(printf, 2, 3))) {
    va_list args;
    va_start(args, fmt);
    const int r = append_(fmt, args);
    va_end(args);
    return r;
  }

  DiagLine& kv(const char* key, float v, int decimals = 3) {
    printf(",%s=%.*f", key, decimals, static_cast<double>(v));
    return *this;
  }
  DiagLine& kv(const char* key, bool v) { printf(",%s=%u", key, v ? 1u : 0u); return *this; }
  DiagLine& kv(const char* key, int v) { printf(",%s=%d", key, v); return *this; }
  DiagLine& kv(const char* key, unsigned v) { printf(",%s=%u", key, v); return *this; }
  DiagLine& kv(const char* key, unsigned long v) { printf(",%s=%lu", key, v); return *this; }
  DiagLine& hex(uint32_t v) { printf(",%lx", static_cast<unsigned long>(v)); return *this; }
  DiagLine& bits(float v) { return hex(floatBits(v)); }

  bool overflowed() const { return overflow_; }
  const char* c_str() const { return buf_; }

  // Writes the line and its newline in one call when the port has room for
  // all of it; otherwise writes nothing and returns false.
  template <class Port>
  bool emit(Port& port) {
    if (overflow_) return false;
    while (n_ > 1 && (buf_[n_ - 1] == '\n' || buf_[n_ - 1] == '\r')) --n_;
    if (n_ + 1 >= sizeof(buf_)) return false;
    buf_[n_] = '\n';
    const size_t len = n_ + 1;
    if (port.availableForWrite() < static_cast<int>(len)) return false;
    return port.write(reinterpret_cast<const uint8_t*>(buf_), len) == len;
  }

 private:
  int append_(const char* fmt, va_list args) {
    if (overflow_) return 0;
    const size_t room = sizeof(buf_) - 1 - n_;  // keep one byte for '\n'
    const int r = std::vsnprintf(buf_ + n_, room, fmt, args);
    if (r < 0 || static_cast<size_t>(r) >= room) {
      overflow_ = true;
      return 0;
    }
    n_ += static_cast<size_t>(r);
    return r;
  }

  char buf_[1024];
  size_t n_ = 0;
  bool overflow_ = false;
};

// Loop timing and acquisition counters over one reporting window, plus the
// frame sequence that persists across windows.
struct InsLoopHealth {
  static constexpr float kFastGyroDps = 30.0f;  // the startup gravity gate's veto

  // Persistent.
  uint32_t frame_seq = 0;     // accepted frames since boot (R-line sequence)
  uint32_t lines_dropped = 0; // diagnostic lines the port could not take

  // Window.
  uint32_t polls = 0, frames = 0, mask_acc_only = 0, mask_gyr_only = 0, mask_none = 0;
  uint32_t overruns = 0, tick_max_us = 0;
  uint32_t dt_min_us = 0, dt_max_us = 0;
  uint32_t gyro_fast = 0;     // frames with |gyro| above kFastGyroDps
  float gyro_max_dps = 0.0f, acc_max_g = 0.0f;

  void poll(uint32_t mask, uint32_t acc_bit, uint32_t gyr_bit) {
    ++polls;
    const bool a = (mask & acc_bit) != 0, g = (mask & gyr_bit) != 0;
    if (a && !g) ++mask_acc_only;
    else if (g && !a) ++mask_gyr_only;
    else if (!a && !g) ++mask_none;
  }

  // An accepted frame: host sample time and the calibrated rates/forces.
  void frame(uint32_t sample_us, const Eigen::Vector3f& w_rad_s, const Eigen::Vector3f& a_mps2,
             float g_mps2) {
    ++frame_seq;
    ++frames;
    if (have_last_) {
      const uint32_t dt = sample_us - last_frame_us_;
      if (frames == 1 || dt < dt_min_us) dt_min_us = dt;
      if (dt > dt_max_us) dt_max_us = dt;
    }
    last_frame_us_ = sample_us;
    have_last_ = true;
    const float w_dps = w_rad_s.norm() * 57.29577951f;
    if (std::isfinite(w_dps)) {
      if (w_dps > gyro_max_dps) gyro_max_dps = w_dps;
      if (w_dps > kFastGyroDps) ++gyro_fast;
    }
    const float a_g = a_mps2.norm() / g_mps2;
    if (std::isfinite(a_g) && a_g > acc_max_g) acc_max_g = a_g;
  }

  // Loop body duration (before the pacing wait) against the loop period.
  void tick(uint32_t busy_us, uint32_t period_us) {
    if (busy_us > period_us) ++overruns;
    if (busy_us > tick_max_us) tick_max_us = busy_us;
  }

  void dropped() { ++lines_dropped; }

  void addTo(DiagLine& line) const {
    line.kv("seq", static_cast<unsigned long>(frame_seq))
        .kv("polls", static_cast<unsigned long>(polls))
        .kv("frames", static_cast<unsigned long>(frames))
        .kv("acc_only", static_cast<unsigned long>(mask_acc_only))
        .kv("gyr_only", static_cast<unsigned long>(mask_gyr_only))
        .kv("no_data", static_cast<unsigned long>(mask_none))
        .kv("overruns", static_cast<unsigned long>(overruns))
        .kv("tick_max_us", static_cast<unsigned long>(tick_max_us))
        .kv("dt_min_us", static_cast<unsigned long>(dt_min_us))
        .kv("dt_max_us", static_cast<unsigned long>(dt_max_us))
        .kv("gyro_max_dps", gyro_max_dps, 1)
        .kv("gyro_fast", static_cast<unsigned long>(gyro_fast))
        .kv("acc_max_g", acc_max_g, 2)
        .kv("lines_dropped", static_cast<unsigned long>(lines_dropped));
  }

  void resetWindow() {
    polls = frames = mask_acc_only = mask_gyr_only = mask_none = 0;
    overruns = tick_max_us = dt_min_us = dt_max_us = gyro_fast = 0;
    gyro_max_dps = acc_max_g = 0.0f;
  }

 private:
  uint32_t last_frame_us_ = 0;
  bool have_last_ = false;
};

// R line: one accepted frame, bit-exact. Sample is ImuSample.
template <class Sample>
void rawFrameLine(DiagLine& line, uint32_t seq, const Sample& s) {
  line.reset('R');
  line.hex(seq).hex(s.sample_us).hex(s.mask);
  line.bits(s.a.x()).bits(s.a.y()).bits(s.a.z());
  line.bits(s.w.x()).bits(s.w.y()).bits(s.w.z());
  line.bits(s.m.x()).bits(s.m.y()).bits(s.m.z());
  line.bits(s.tempC);
  const auto& t = s.mag_timing;
  line.hex(t.have ? 1u : 0u).hex(t.last_ms).hex(t.sequence).hex(t.frame_us)
      .hex(t.read_start_us).hex(t.read_end_us);
}

// C line: whether a calibration was loaded, the active magnetometer source
// and the blob bytes the sketch holds (all zero when none was loaded).
template <class Blob>
void calibrationLine(DiagLine& line, bool have_blob, unsigned mag_source, const Blob& blob) {
  line.reset('C');
  line.printf(",%u,%u,", have_blob ? 1u : 0u, mag_source);
  const auto* bytes = reinterpret_cast<const unsigned char*>(&blob);
  for (size_t i = 0; i < sizeof(Blob); ++i) line.printf("%02x", static_cast<unsigned>(bytes[i]));
}

// BMI270/BMM150 settings that decide sample timing, as read back from the
// sensor (readImuTimingRegisters() in AtomS3R_ImuCal.h).
struct ImuTimingRegisters {
  bool bmi270 = false;  // BMI270 register block read
  uint8_t acc_conf = 0, acc_range = 0, gyr_conf = 0, gyr_range = 0, aux_conf = 0, fifo_downs = 0;
  uint8_t aux_dev_id = 0, aux_if_conf = 0, aux_rd_addr = 0, pwr_conf = 0, pwr_ctrl = 0;
  bool bmm150 = false;  // BMM150 registers read through the AUX interface
  uint8_t bmm_chip_id = 0, bmm_power = 0, bmm_mode = 0, bmm_rep_xy = 0, bmm_rep_z = 0;
};

// I line. Output data rates follow the BMI270/BMM150 datasheet encodings.
inline void timingRegistersLine(DiagLine& line, const ImuTimingRegisters& r, unsigned mag_source) {
  static constexpr float kBmiOdrHz[16] = {0, 0.78125f, 1.5625f, 3.125f, 6.25f, 12.5f, 25, 50,
                                          100, 200, 400, 800, 1600, 3200, 0, 0};
  static constexpr float kBmmOdrHz[8] = {10, 2, 6, 8, 15, 20, 25, 30};
  line.reset('I');
  line.kv("bmi270", r.bmi270);
  if (r.bmi270) {
    line.printf(",acc_conf=0x%02x,acc_range=0x%02x,gyr_conf=0x%02x,gyr_range=0x%02x,aux_conf=0x%02x,"
                "fifo_downs=0x%02x,aux_dev_id=0x%02x,aux_if_conf=0x%02x,aux_rd_addr=0x%02x,"
                "pwr_conf=0x%02x,pwr_ctrl=0x%02x",
                r.acc_conf, r.acc_range, r.gyr_conf, r.gyr_range, r.aux_conf, r.fifo_downs,
                r.aux_dev_id, r.aux_if_conf, r.aux_rd_addr, r.pwr_conf, r.pwr_ctrl);
    line.kv("acc_odr_hz", kBmiOdrHz[r.acc_conf & 0x0f], 2)
        .kv("gyr_odr_hz", kBmiOdrHz[r.gyr_conf & 0x0f], 2)
        .kv("aux_odr_hz", kBmiOdrHz[r.aux_conf & 0x0f], 2);
  }
  line.kv("bmm150", r.bmm150);
  if (r.bmm150) {
    line.printf(",bmm_chip_id=0x%02x,bmm_power=0x%02x,bmm_mode=0x%02x,bmm_rep_xy=%u,bmm_rep_z=%u",
                r.bmm_chip_id, r.bmm_power, r.bmm_mode, unsigned(r.bmm_rep_xy), unsigned(r.bmm_rep_z));
    line.kv("bmm_odr_hz", kBmmOdrHz[(r.bmm_mode >> 3) & 0x07], 1);
  }
  line.kv("mag_source", mag_source);
}

// Tilt-compensated magnetic heading [deg] from a published roll/pitch
// (body NED, ZYX) and the calibrated body field: lets a log compare the
// published heading with the raw magnetometer on the same tilt.
inline float compassHeadingFromRollPitchDeg(float roll_deg, float pitch_deg,
                                            const Eigen::Vector3f& m_body_uT) {
  const float r = roll_deg * 0.01745329252f, p = pitch_deg * 0.01745329252f;
  const Eigen::Vector3f down_b(-std::sin(p), std::sin(r) * std::cos(p), std::cos(r) * std::cos(p));
  if (!down_b.allFinite() || !m_body_uT.allFinite()) return NAN;
  const Eigen::Vector3f east_b = down_b.cross(m_body_uT);
  const Eigen::Vector3f north_b = east_b.cross(down_b);
  if (!(east_b.norm() > 1e-6f) || !(north_b.norm() > 1e-6f)) return NAN;
  float h = std::atan2(east_b.x() / east_b.norm(), north_b.x() / north_b.norm()) * 57.29577951f;
  if (h < 0.0f) h += 360.0f;
  return h;
}

// Yaw and tilt one-sigma [deg] of a 3x3 attitude-error covariance, split
// about the world-down axis expressed in the same frame as the covariance.
inline void attitudeSigmaDeg(const Eigen::Matrix3f& P_att, const Eigen::Vector3f& down_in_error_frame,
                             float& yaw_sigma_deg, float& tilt_sigma_deg) {
  yaw_sigma_deg = tilt_sigma_deg = NAN;
  Eigen::Vector3f u = down_in_error_frame;
  const float n = u.norm();
  if (!P_att.allFinite() || !(n > 1e-6f)) return;
  u /= n;
  const float yaw_var = u.dot(P_att * u);
  const float tilt_var = 0.5f * (P_att.trace() - yaw_var);
  yaw_sigma_deg = std::sqrt(std::max(0.0f, yaw_var)) * 57.29577951f;
  tilt_sigma_deg = std::sqrt(std::max(0.0f, tilt_var)) * 57.29577951f;
}

// The diagnostic stream of one sketch: line scheduling, loop health and the
// boot records it repeats so that a capture started late still has them.
template <class Blob>
class InsDiagStream {
 public:
  static constexpr uint32_t kStatePeriodMs = 100;    // D
  static constexpr uint32_t kMagPeriodMs = 1000;     // M
  static constexpr uint32_t kHeaderPeriodMs = 5000;  // H, I, C

  InsDiagStream(const char* sketch, float loop_hz) : sketch_(sketch), loop_hz_(loop_hz) {}

  // Gives a USB serial reader that is (re)connecting after a reset a chance
  // to see the stream from its start. Returns at once when one is attached.
  template <class Port>
  static void waitForReader(Port& port, uint32_t timeout_ms) {
    const uint32_t t0 = millis();
    while (!port && static_cast<uint32_t>(millis() - t0) < timeout_ms) delay(10);
  }

  // Records what the sketch is running with; H/I/C go out on the next
  // periodic() call. Call at the end of begin() and after a recalibration.
  void boot(const ImuTimingRegisters& regs, bool have_blob, unsigned mag_source, const Blob& blob) {
    regs_ = regs;
    have_blob_ = have_blob;
    mag_source_ = mag_source;
    blob_ = blob;
    booted_ = true;
    header_due_ = true;
  }

  // One loop poll. frame is the accepted sample (nullptr when the poll had
  // none); w/a are its calibrated rate and specific force.
  template <class Port, class Sample>
  void poll(Port& port, uint32_t mask, uint32_t acc_bit, uint32_t gyr_bit, const Sample* frame,
            const Eigen::Vector3f& w, const Eigen::Vector3f& a, float g_mps2) {
    health_.poll(mask, acc_bit, gyr_bit);
    if (!frame) return;
    rawFrameLine(line_, health_.frame_seq, *frame);
    emit(port, line_);
    health_.frame(frame->sample_us, w, a, g_mps2);
  }

  void tick(uint32_t busy_us, uint32_t period_us) { health_.tick(busy_us, period_us); }

  // A D line is due: returns one started with the time and loop health; the
  // sketch appends its fields and passes it to emitState().
  bool stateDue(uint32_t now_ms) const {
    return !have_state_ || static_cast<uint32_t>(now_ms - last_state_ms_) >= kStatePeriodMs;
  }
  DiagLine& beginState(uint32_t now_ms) {
    have_state_ = true;
    last_state_ms_ = now_ms;
    line_.reset('D');
    line_.kv("t_ms", static_cast<unsigned long>(now_ms));
    health_.addTo(line_);
    return line_;
  }
  template <class Port>
  void emitState(Port& port, DiagLine& line) {
    emit(port, line);
    health_.resetWindow();
  }

  // M line once a second (printer(DiagLine&) appends the magnetic timing),
  // H/I/C after boot() and every few seconds.
  template <class Port, class MagPrinter>
  void periodic(Port& port, uint32_t now_ms, MagPrinter&& printer) {
    if (!have_mag_ || static_cast<uint32_t>(now_ms - last_mag_ms_) >= kMagPeriodMs) {
      have_mag_ = true;
      last_mag_ms_ = now_ms;
      line_.reset('M').printf(",");
      printer(line_);
      emit(port, line_);
    }
    if (!booted_) return;
    if (header_due_ || static_cast<uint32_t>(now_ms - last_header_ms_) >= kHeaderPeriodMs) {
      header_due_ = false;
      last_header_ms_ = now_ms;
      line_.reset('H').printf(",sketch=%s", sketch_);
      line_.kv("format", kInsDiagFrameFormat).kv("loop_hz", loop_hz_, 1)
           .kv("seq", static_cast<unsigned long>(health_.frame_seq));
      if (header_extra_) line_.printf(",%s", header_extra_);
      emit(port, line_);
      timingRegistersLine(line_, regs_, mag_source_);
      emit(port, line_);
      calibrationLine(line_, have_blob_, mag_source_, blob_);
      emit(port, line_);
    }
  }

  // E line: an application event the replay must know about.
  template <class Port>
  void event(Port& port, const char* what) {
    line_.reset('E').printf(",%s", what);
    line_.kv("seq", static_cast<unsigned long>(health_.frame_seq));
    emit(port, line_);
  }

  template <class Port>
  void emit(Port& port, DiagLine& line) {
    if (!line.emit(port)) health_.dropped();
  }

  const InsLoopHealth& health() const { return health_; }
  // For a sketch that delivers frames itself instead of through poll().
  InsLoopHealth& health() { return health_; }
  // Extra key=value text for the H line, e.g. "acquisition=fifo".
  void setHeaderExtra(const char* extra) { header_extra_ = extra; }

 private:
  const char* sketch_;
  const char* header_extra_ = nullptr;
  float loop_hz_;
  InsLoopHealth health_{};
  DiagLine line_{};  // the one line buffer, reused for every line
  ImuTimingRegisters regs_{};
  Blob blob_{};
  bool have_blob_ = false, booted_ = false, header_due_ = false;
  unsigned mag_source_ = 0;
  bool have_state_ = false, have_mag_ = false;
  uint32_t last_state_ms_ = 0, last_mag_ms_ = 0, last_header_ms_ = 0;
};

// Counts sudden attitude-covariance re-seeds (e.g. a tilt watchdog re-lock)
// from the published yaw sigma: a jump far larger than one sample of process
// noise can produce.
struct AttitudeReseedDetector {
  static constexpr float kJumpDeg = 30.0f;
  float last_yaw_sigma_deg = NAN;
  uint32_t count = 0;
  void update(float yaw_sigma_deg) {
    if (std::isfinite(last_yaw_sigma_deg) && std::isfinite(yaw_sigma_deg) &&
        yaw_sigma_deg - last_yaw_sigma_deg > kJumpDeg) {
      ++count;
    }
    last_yaw_sigma_deg = yaw_sigma_deg;
  }
  void reset() { last_yaw_sigma_deg = NAN; count = 0; }
};

}  // namespace atoms3r_ical
