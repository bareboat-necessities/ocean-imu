#pragma once

// Copyright 2026, Mikhail Grushinskiy
//
// BMI270 FIFO acquisition with a BMI270-triggered BMM150 (host-testable; the
// Device is the BMI270 with the M5 I2C_Device-style API
// readRegister(reg, buf, len) and writeRegister8(reg, value), as in
// AtomS3R_Bmm150AuxPreset.h).
//
// The shipping sketches poll the BMI270 data registers at their default rates
// (accelerometer 100 Hz, gyro 200 Hz) and the BMM150 runs free at 30 Hz while
// the BMI270 reads it at 25 Hz, so frames are paired at about 100 Hz on host
// time and the magnetic sample age wanders (docs/magnetic-timing.md). This
// acquisition instead:
//
//  * sets accelerometer and gyro to the same 200 Hz output rate, with the
//    ranges M5Unified's unit conversion assumes (+-8 g, +-2000 deg/s);
//  * puts every accelerometer+gyro sample and every magnetometer readout in
//    the BMI270 FIFO (header mode), so nothing is lost while the host is busy
//    and each sample's time follows from its position on the sensor's own
//    25.6 kHz sensortime clock (a 200 Hz frame every 128 ticks);
//  * runs the BMM150 in forced mode: at every 25 Hz AUX cycle the BMI270
//    writes the forced-mode command (aux_fcu_write_en) and reads the
//    previous conversion, so each readout is one complete conversion with a
//    fixed relation to the sensor clock instead of a beat between two clocks.
//
// The magnetic delay that relation produces is measured on the device
// (atomS3R_imu_fifo_timing), not assumed here.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include "AtomS3R/AtomS3R_ImuUnits.h"

namespace atoms3r_ical {

// BMI270 registers and encodings used by this acquisition (Bosch BMI270
// datasheet, register map and FIFO chapters).
struct Bmi270 {
  static constexpr uint8_t CHIP_ID = 0x00, ERR_REG = 0x02, STATUS = 0x03, AUX_DATA = 0x04,
                           SENSORTIME_0 = 0x18, FIFO_LENGTH_0 = 0x24, FIFO_DATA = 0x26,
                           ACC_CONF = 0x40, ACC_RANGE = 0x41, GYR_CONF = 0x42, GYR_RANGE = 0x43,
                           AUX_CONF = 0x44, FIFO_DOWNS = 0x45, FIFO_CONFIG_0 = 0x48, FIFO_CONFIG_1 = 0x49,
                           AUX_DEV_ID = 0x4B, AUX_IF_CONF = 0x4C, AUX_RD_ADDR = 0x4D, AUX_WR_ADDR = 0x4E,
                           AUX_WR_DATA = 0x4F, PWR_CONF = 0x7C, PWR_CTRL = 0x7D, CMD = 0x7E;
  static constexpr uint8_t CHIP_ID_VALUE = 0x24, STATUS_AUX_BUSY = 0x04, CMD_FIFO_FLUSH = 0xB0;
  // PWR_CTRL: aux | gyr | acc | temp.
  static constexpr uint8_t PWR_AUX = 0x01, PWR_ALL = 0x0F;
  // AUX_IF_CONF: manual mode | FCU write | manual burst 8 | data burst 8.
  static constexpr uint8_t AUX_MANUAL = 0x80, AUX_DATA_MODE_FCU_BURST8 = 0x4F;
  // FIFO_CONFIG_1: gyr | acc | aux | header.
  static constexpr uint8_t FIFO_ALL_HEADER = 0xF0;
  // The FIFO holds 2 kB: about 146 frames of accelerometer+gyro with a
  // magnetometer readout every 8th (0.73 s at 200 Hz). A longer pause in
  // reading overwrites the oldest frames; the BMI270 reports them with a
  // skip frame. Read it in one burst: a read that ends inside a frame
  // breaks the frame alignment of the next read.
  static constexpr size_t FIFO_BYTES = 2048;
  // Sensortime: 24 bits at 25.6 kHz; a 200 Hz frame every 128 ticks.
  static constexpr double TICK_US = 39.0625;
  static constexpr uint32_t FRAME_TICKS = 128;
  static constexpr uint32_t SENSORTIME_MASK = 0xFFFFFF;
  // BMM150 behind the AUX interface.
  static constexpr uint8_t BMM_I2C_ADDR8 = 0x10 << 1, BMM_CHIP_ID = 0x40, BMM_CHIP_ID_VALUE = 0x32,
                           BMM_DATA_X = 0x42, BMM_MODE = 0x4C, BMM_MODE_ODR_BITS = 0x38,
                           BMM_OPMODE_FORCED = 0x02, BMM_OPMODE_SLEEP = 0x06;
};

struct FifoAcquisitionConfig {
  uint8_t acc_conf = 0xA9;    // 200 Hz, normal filter, performance mode
  uint8_t acc_range = 0x02;   // +-8 g, the range M5Unified converts with
  uint8_t gyr_conf = 0xA9;    // 200 Hz, normal filter, performance filter (reset value)
  uint8_t gyr_range = 0x00;   // +-2000 deg/s
  uint8_t aux_odr = 0x06;     // 25 Hz
  uint8_t aux_offset = 0;     // trigger-to-readout offset, 2.5 ms units; 0 = longest
  uint8_t fifo_downs = 0x88;  // filtered accelerometer and gyro, no down-sampling
};

// Register state after setup, read back from the sensor.
struct FifoAcquisitionReadback {
  uint8_t chip_id = 0, err_reg = 0, acc_conf = 0, acc_range = 0, gyr_conf = 0, gyr_range = 0,
          aux_conf = 0, fifo_downs = 0, fifo_config_0 = 0, fifo_config_1 = 0, aux_dev_id = 0,
          aux_if_conf = 0, aux_rd_addr = 0, aux_wr_addr = 0, aux_wr_data = 0, pwr_conf = 0, pwr_ctrl = 0;
  uint8_t bmm_chip_id = 0, bmm_mode_before = 0;
};

// Configures the BMI270 for FIFO acquisition with a forced-mode BMM150.
// Call once after M5Unified's IMU initialisation and configureAtomS3RMag()
// (which sets the BMM150 repetitions and reads its trim registers). The
// shipping polling path (M5.Imu.update()) must not be used afterwards.
template <class Device>
class Bmi270FifoSetup {
 public:
  Bmi270FifoSetup(Device* device, void (*wait_ms)(uint32_t)) : dev_(device), wait_(wait_ms) {}

  const char* failure() const { return failure_; }
  const FifoAcquisitionReadback& readback() const { return rb_; }

  bool apply(const FifoAcquisitionConfig& c) {
    using R = Bmi270;
    failure_ = "BMI270 read failed";
    if (!dev_) { failure_ = "no BMI270"; return false; }
    uint8_t pwr = 0, stale_err = 0;
    if (!rd_(R::CHIP_ID, rb_.chip_id) || !rd_(R::PWR_CTRL, pwr)) return false;
    if (rb_.chip_id != R::CHIP_ID_VALUE) { failure_ = "not a BMI270"; return false; }
    // Clears an AUX error flag left by earlier initialisation (read-to-clear).
    if (!rd_(R::ERR_REG, stale_err)) return false;

    // Stop AUX polling, then talk to the BMM150 in manual mode.
    failure_ = "AUX manual mode failed";
    if (!wr_(R::PWR_CTRL, uint8_t(pwr & ~R::PWR_AUX))) return false;
    wait_(2);
    if (!wr_(R::AUX_IF_CONF, R::AUX_MANUAL) || !rd_(R::AUX_DEV_ID, rb_.aux_dev_id)) return false;
    if (rb_.aux_dev_id != R::BMM_I2C_ADDR8) { failure_ = "AUX device is not at 0x10"; return false; }
    failure_ = "BMM150 read failed";
    if (!auxRead_(R::BMM_CHIP_ID, rb_.bmm_chip_id) || !auxRead_(R::BMM_MODE, rb_.bmm_mode_before)) return false;
    if (rb_.bmm_chip_id != R::BMM_CHIP_ID_VALUE) { failure_ = "BMM150 chip id mismatch"; return false; }

    // Sleep, then one forced-mode command. Written in manual mode it starts
    // the first conversion now; AUX_WR_ADDR/AUX_WR_DATA keep it, and with
    // FCU write enabled the BMI270 repeats it on every AUX cycle.
    failure_ = "BMM150 write failed";
    const uint8_t odr = uint8_t(rb_.bmm_mode_before & R::BMM_MODE_ODR_BITS);
    if (!auxWrite_(R::BMM_MODE, uint8_t(odr | R::BMM_OPMODE_SLEEP))) return false;
    wait_(5);
    if (!auxWrite_(R::BMM_MODE, uint8_t(odr | R::BMM_OPMODE_FORCED))) return false;

    // Data mode: burst-read X..RHALL at the AUX rate (same order as M5Unified).
    failure_ = "BMI270 write failed";
    const uint8_t aux_conf = uint8_t((c.aux_offset & 0x0F) << 4 | (c.aux_odr & 0x0F));
    if (!wr_(R::AUX_IF_CONF, R::AUX_DATA_MODE_FCU_BURST8) || !wr_(R::AUX_RD_ADDR, R::BMM_DATA_X) ||
        !wr_(R::AUX_CONF, aux_conf) || !wr_(R::ACC_CONF, c.acc_conf) || !wr_(R::ACC_RANGE, c.acc_range) ||
        !wr_(R::GYR_CONF, c.gyr_conf) || !wr_(R::GYR_RANGE, c.gyr_range) ||
        !wr_(R::FIFO_DOWNS, c.fifo_downs) || !wr_(R::FIFO_CONFIG_0, 0x00) ||
        !wr_(R::FIFO_CONFIG_1, R::FIFO_ALL_HEADER) || !wr_(R::PWR_CTRL, R::PWR_ALL))
      return false;
    wait_(60);  // gyro start-up
    if (!wr_(R::CMD, R::CMD_FIFO_FLUSH)) return false;
    wait_(2);

    failure_ = "BMI270 read-back failed";
    uint8_t conf[6], fifo[2], aux[5], pwr2[2];
    if (!dev_->readRegister(R::ACC_CONF, conf, sizeof(conf)) ||
        !dev_->readRegister(R::FIFO_CONFIG_0, fifo, sizeof(fifo)) ||
        !dev_->readRegister(R::AUX_DEV_ID, aux, sizeof(aux)) ||
        !dev_->readRegister(R::PWR_CONF, pwr2, sizeof(pwr2)) || !rd_(R::ERR_REG, rb_.err_reg))
      return false;
    rb_.acc_conf = conf[0]; rb_.acc_range = conf[1]; rb_.gyr_conf = conf[2]; rb_.gyr_range = conf[3];
    rb_.aux_conf = conf[4]; rb_.fifo_downs = conf[5];
    rb_.fifo_config_0 = fifo[0]; rb_.fifo_config_1 = fifo[1];
    rb_.aux_dev_id = aux[0]; rb_.aux_if_conf = aux[1]; rb_.aux_rd_addr = aux[2];
    rb_.aux_wr_addr = aux[3]; rb_.aux_wr_data = aux[4];
    rb_.pwr_conf = pwr2[0]; rb_.pwr_ctrl = pwr2[1];
    failure_ = "read-back differs from configuration";
    if (rb_.acc_conf != c.acc_conf || rb_.acc_range != c.acc_range || rb_.gyr_conf != c.gyr_conf ||
        rb_.gyr_range != c.gyr_range || rb_.aux_conf != aux_conf || rb_.fifo_downs != c.fifo_downs ||
        rb_.fifo_config_1 != R::FIFO_ALL_HEADER || rb_.aux_if_conf != R::AUX_DATA_MODE_FCU_BURST8 ||
        rb_.aux_rd_addr != R::BMM_DATA_X || rb_.aux_wr_addr != R::BMM_MODE ||
        rb_.aux_wr_data != uint8_t(odr | R::BMM_OPMODE_FORCED) || (rb_.pwr_ctrl & R::PWR_ALL) != R::PWR_ALL)
      return false;
    // ERR_REG: fatal (bit 0), FIFO (bit 6) or AUX (bit 7) error.
    if (rb_.err_reg & 0xC1) { failure_ = "BMI270 reports an error"; return false; }
    failure_ = "none";
    return true;
  }

 private:
  bool rd_(uint8_t reg, uint8_t& v) { return dev_->readRegister(reg, &v, 1); }
  bool wr_(uint8_t reg, uint8_t v) { return dev_->writeRegister8(reg, v); }
  bool waitIdle_() {
    for (int i = 0; i < 20; ++i) {
      uint8_t st = 0;
      if (!rd_(Bmi270::STATUS, st)) return false;
      if (!(st & Bmi270::STATUS_AUX_BUSY)) return true;
      wait_(1);
    }
    return false;
  }
  bool auxRead_(uint8_t reg, uint8_t& v) {
    return wr_(Bmi270::AUX_RD_ADDR, reg) && waitIdle_() && rd_(Bmi270::AUX_DATA, v);
  }
  bool auxWrite_(uint8_t reg, uint8_t v) {
    return wr_(Bmi270::AUX_WR_DATA, v) && wr_(Bmi270::AUX_WR_ADDR, reg) && waitIdle_();
  }

  Device* dev_;
  void (*wait_)(uint32_t);
  FifoAcquisitionReadback rb_{};
  const char* failure_ = "not run";
};

// ------------------------------------------------------------------ frames --

// One header-mode data frame. Data order in the FIFO is AUX, gyro, accelerometer.
struct FifoFrame {
  uint8_t header = 0;
  bool acc = false, gyr = false, aux = false;
  int16_t a[3]{}, g[3]{};
  uint8_t m[8]{};  // BMM150 0x42..0x49: X, Y, Z, RHALL (data-ready in byte 6 bit 0)
};

struct FifoParseStats {
  uint32_t frames = 0, acc_gyr = 0, acc_only = 0, gyr_only = 0, aux_only = 0, aux = 0;
  uint32_t skip_frames = 0, skipped = 0, sensortime_frames = 0, config_frames = 0, unknown = 0;
};

// Parses header-mode FIFO bytes. Calls on_frame(const FifoFrame&) for every
// data frame and on_skip(count) for every skip frame (frames the FIFO dropped
// when it overflowed). Returns the number of bytes consumed; parsing stops at
// the over-read marker, at an unknown header (counted) or before a frame that
// is not complete in the buffer.
template <class OnFrame, class OnSkip>
size_t parseBmi270Fifo(const uint8_t* d, size_t n, FifoParseStats& st, OnFrame&& on_frame,
                       OnSkip&& on_skip) {
  auto s16 = [](const uint8_t* p) { return int16_t(uint16_t(p[0] | uint16_t(p[1]) << 8)); };
  size_t i = 0;
  while (i < n) {
    const uint8_t h = d[i];
    if ((h & 0xC0) == 0x80) {
      if ((h & 0x3C) == 0) break;  // over-read: no more data
      if (h & 0x20) { ++st.unknown; break; }
      FifoFrame f;
      f.header = h;
      f.aux = (h & 0x10) != 0;
      f.gyr = (h & 0x08) != 0;
      f.acc = (h & 0x04) != 0;
      const size_t len = 1 + (f.aux ? 8 : 0) + (f.gyr ? 6 : 0) + (f.acc ? 6 : 0);
      if (i + len > n) break;
      size_t k = i + 1;
      if (f.aux) { std::memcpy(f.m, d + k, 8); k += 8; }
      if (f.gyr) { for (size_t j = 0; j < 3; ++j) f.g[j] = s16(d + k + 2 * j); k += 6; }
      if (f.acc) { for (size_t j = 0; j < 3; ++j) f.a[j] = s16(d + k + 2 * j); k += 6; }
      ++st.frames;
      if (f.aux) ++st.aux;
      if (f.acc && f.gyr) ++st.acc_gyr;
      else if (f.acc) ++st.acc_only;
      else if (f.gyr) ++st.gyr_only;
      else ++st.aux_only;
      on_frame(f);
      i += len;
    } else if (h == 0x40 || h == 0x48) {  // skip frame, input-configuration change
      if (i + 2 > n) break;
      if (h == 0x40) { ++st.skip_frames; st.skipped += d[i + 1]; on_skip(unsigned(d[i + 1])); }
      else ++st.config_frames;
      i += 2;
    } else if (h == 0x44) {  // sensortime frame
      if (i + 4 > n) break;
      ++st.sensortime_frames;
      i += 4;
    } else {
      ++st.unknown;
      break;
    }
  }
  return i;
}

// The units readImuMapped() delivers. M5Unified reorders the AtomS3R BMI270
// axes as (y, -x, z) and scales them to g and deg/s; the project then maps
// those to body NED (AtomS3R_ImuUnits.h).
inline Eigen::Vector3f bmi270AccelToBodyNed(const int16_t raw[3], float g_per_lsb) {
  return accel_nominal_g_to_body_ned_si_(float(raw[1]) * g_per_lsb, -float(raw[0]) * g_per_lsb,
                                         float(raw[2]) * g_per_lsb);
}
inline Eigen::Vector3f bmi270GyroToBodyNed(const int16_t raw[3], float dps_per_lsb) {
  return map_sensor_xyz_to_body_ned_(float(raw[1]) * dps_per_lsb, -float(raw[0]) * dps_per_lsb,
                                     float(raw[2]) * dps_per_lsb, ImuCalCfg::DEG2RAD);
}
inline float bmi270AccelGPerLsb(uint8_t acc_range) { return float(2 << (acc_range & 3)) / 32768.0f; }
inline float bmi270GyroDpsPerLsb(uint8_t gyr_range) { return float(2000 >> (gyr_range & 7)) / 32768.0f; }

// BMM150 data-ready bit of an AUX payload, and the payload without its
// status bits (for repeated-reading detection).
inline bool bmm150DataReady(const uint8_t m[8]) { return (m[6] & 0x01) != 0; }
inline bool bmm150SamePayload(const uint8_t a[8], const uint8_t b[8]) {
  static constexpr uint8_t mask[8] = {0xf8, 0xff, 0xf8, 0xff, 0xfe, 0xff, 0xfc, 0xff};
  for (int i = 0; i < 8; ++i)
    if ((a[i] & mask[i]) != (b[i] & mask[i])) return false;
  return true;
}

// ------------------------------------------------------------ sensor clock --

// Maps 200 Hz frame indices to sensortime and to host time.
//
// Frame k (counted from the FIFO flush, skipped frames included) was written
// at T0 + 128 k ticks. Each FIFO read gives the sensortime read just before
// the FIFO length, and the number of frames up to that length: the last of
// them is already written and the next is not, which bounds T0 from both
// sides. The bounds tighten as reads fall at different phases. Frames lost
// without a skip frame show up as more elapsed time than the count covers;
// they are counted (lostFrames()) and later frames keep their place on the
// grid. Host time follows from a line fitted between sensortime and the host
// microsecond clock over the whole run.
class SensorClockMap {
 public:
  static constexpr double kNominalUsPerTick = Bmi270::TICK_US;
  static constexpr double kGridToleranceTicks = 8.0;

  void reset() { *this = SensorClockMap{}; }

  // host_us: host clock just before the register read; sensortime24: the
  // 24-bit sensortime it returned; frames: frames received up to the FIFO
  // length read with it (skipped frames included).
  void observe(uint32_t host_us, uint32_t sensortime24, uint64_t frames) {
    if (!have_) {
      have_ = true;
      st_ = sensortime24;
      host_ = host_us;
      st0_ = st_;
      host0_ = host_;
      offset_ = double(host_) - kNominalUsPerTick * double(st_);
    } else {
      st_ += (sensortime24 - last_st24_) & Bmi270::SENSORTIME_MASK;
      host_ += uint32_t(host_us - last_host_);
    }
    last_st24_ = sensortime24;
    last_host_ = host_us;
    ++reads_;

    // Slope over the whole run once it spans 160 ms (a BMI270 clock can be
    // 1 % off nominal), refined as the span grows. The line is re-anchored
    // at the current tick when the slope changes, so frame times stay
    // continuous; the offset then follows the read-time jitter slowly.
    if (st_ - st0_ >= (1u << 12)) {
      const double s = double(host_ - host0_) / double(st_ - st0_);
      offset_ += (us_per_tick_ - s) * double(st_);
      us_per_tick_ = s;
    }
    const double o = double(host_) - us_per_tick_ * double(st_);
    offset_ += (o - offset_) * (reads_ < 50 ? 1.0 / double(reads_) : 0.02);

    if (frames == 0) return;
    const double t = double(st_), F = double(Bmi270::FRAME_TICKS), tol = kGridToleranceTicks;
    auto bounds = [&](double& lo, double& hi) {
      const double n = double(frames + lost_);
      lo = t - F * n - tol;
      hi = t - F * (n - 1.0) + tol;
    };
    double lo = 0, hi = 0;
    bounds(lo, hi);
    if (have_grid_ && lo > hi_) {
      // More time has passed than the counted frames cover: frames were
      // lost without a skip frame. Count them, so later frames keep their
      // place on the grid.
      lost_ += uint64_t(std::ceil((lo - hi_) / F - 1e-9));
      bounds(lo, hi);
    }
    if (!have_grid_ || std::max(lo_, lo) > std::min(hi_, hi)) {
      if (have_grid_) ++grid_breaks_;  // count ahead of time (a clock glitch): restart
      have_grid_ = true;
      lo_ = lo;
      hi_ = hi;
    } else {
      lo_ = std::max(lo_, lo);
      hi_ = std::min(hi_, hi);
    }
  }

  bool valid() const { return have_grid_; }
  // The frame grid is pinned to within a few ticks (reads have fallen at
  // enough different phases): frame times no longer move as it tightens.
  bool settled() const { return have_grid_ && hi_ - lo_ <= 2.0 * kGridToleranceTicks + 4.0; }
  double t0Ticks() const { return 0.5 * (lo_ + hi_); }
  double t0UncertaintyTicks() const { return 0.5 * (hi_ - lo_); }
  // Sensortime of the k-th received frame, lost frames detected so far
  // included: call it right after observe() for the frames of that read.
  double frameTicks(uint64_t k) const { return t0Ticks() + double(Bmi270::FRAME_TICKS) * double(k + lost_); }
  // Host time [us] of the k-th received frame on the fitted
  // sensortime-to-host line.
  uint32_t frameHostUs(uint64_t k) const {
    return uint32_t(int64_t(std::llround(offset_ + us_per_tick_ * frameTicks(k))));
  }
  double usPerTick() const { return us_per_tick_; }
  double frameSpacingUs() const { return us_per_tick_ * double(Bmi270::FRAME_TICKS); }
  double clockPpm() const { return (us_per_tick_ / kNominalUsPerTick - 1.0) * 1e6; }
  uint32_t gridBreaks() const { return grid_breaks_; }
  uint64_t lostFrames() const { return lost_; }
  uint64_t sensortimeTicks() const { return st_; }
  uint64_t hostUs() const { return host_; }

 private:
  bool have_ = false, have_grid_ = false;
  uint32_t last_st24_ = 0, last_host_ = 0, grid_breaks_ = 0;
  uint64_t st_ = 0, host_ = 0, st0_ = 0, host0_ = 0, reads_ = 0, lost_ = 0;
  double us_per_tick_ = kNominalUsPerTick, offset_ = 0.0, lo_ = 0.0, hi_ = 0.0;
};

}  // namespace atoms3r_ical
