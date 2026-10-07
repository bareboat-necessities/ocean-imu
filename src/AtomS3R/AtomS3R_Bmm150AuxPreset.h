#pragma once

// Copyright 2026, Mikhail Grushinskiy
//
// BMM150 repetition preset written through the BMI270 AUX interface
// (host-testable; Device is the BMI270, M5 I2C_Device-style API:
// readRegister(reg, buf, len) and writeRegister8(reg, value)).
//
// On AtomS3R the BMM150 sits behind the BMI270 AUX bus only. M5Unified 0.2.13
// soft-resets it, selects normal mode at 30 Hz and leaves the repetition
// registers at their reset value (1 XY and 1 Z repetition, the noisiest
// setting), then lets the BMI270 poll its data registers ("data mode").
// This class briefly switches the AUX interface to manual mode, changes only
// the repetition registers, reads them back and returns the AUX interface to
// exactly the data-mode configuration it found. Nothing persists: the next
// M5Unified initialisation (every boot) soft-resets the BMM150 again.
//
// The low-noise preset keeps the 30 Hz output rate: 47 XY and 41 Z
// repetitions take 145*47 + 500*41 + 980 us = 28.3 ms per measurement
// (Bosch BMM150 datasheet, section 5.6), below the 33.3 ms period.

#include <cstdint>
#include <cstddef>

namespace atoms3r_ical {

struct Bmm150RegState {
  uint8_t chip_id = 0, power = 0, mode = 0, rep_xy = 0, rep_z = 0;
  unsigned nXY() const { return 2u * rep_xy + 1u; }
  unsigned nZ() const { return rep_z + 1u; }
};

template <class Device>
class Bmm150AuxPreset {
 public:
  // BMI270 registers.
  static constexpr uint8_t STATUS = 0x03, AUX_DATA = 0x04, AUX_DEV_ID = 0x4B, AUX_IF_CONF = 0x4C,
                           AUX_RD_ADDR = 0x4D, AUX_WR_ADDR = 0x4E, AUX_WR_DATA = 0x4F, PWR_CTRL = 0x7D;
  static constexpr uint8_t STATUS_AUX_BUSY = 0x04, PWR_AUX_EN = 0x01, AUX_MANUAL_EN = 0x80;
  // BMM150 registers.
  static constexpr uint8_t BMM_CHIP_ID = 0x40, BMM_POWER = 0x4B, BMM_MODE = 0x4C, BMM_REP_XY = 0x51,
                           BMM_REP_Z = 0x52;
  static constexpr uint8_t LOW_NOISE_REP_XY = 23, LOW_NOISE_REP_Z = 40;  // nXY=47, nZ=41

  Bmm150AuxPreset(Device* bmi270, void (*wait_ms)(uint32_t)) : dev_(bmi270), wait_(wait_ms) {}

  const char* failure() const { return failure_; }

  // Read-only query of the BMM150 configuration (still needs manual AUX mode).
  bool query(Bmm150RegState& out) {
    if (!enter_()) return false;
    const bool ok = readState_(out);
    return leave_() && ok;
  }

  // Reads count consecutive BMM150 registers starting at first (e.g. the
  // factory trim block 0x5D..0x71) in manual AUX mode.
  bool readBlock(uint8_t first, uint8_t count, uint8_t* out) {
    if (!enter_()) return false;
    bool ok = true;
    for (uint8_t i = 0; ok && i < count; ++i) ok = auxRead_(uint8_t(first + i), out[i]);
    if (!ok) failure_ = "AUX read failed";
    return leave_() && ok;
  }

  // Applies the repetitions, keeping the current output-data-rate bits and
  // returning to normal mode. before/after are the read-back configurations.
  bool apply(uint8_t rep_xy, uint8_t rep_z, Bmm150RegState& before, Bmm150RegState& after) {
    if (!enter_()) return false;
    bool ok = readState_(before);
    if (ok && before.chip_id != 0x32) { failure_ = "BMM150 chip id mismatch"; ok = false; }
    if (ok && !(before.power & 0x01)) { failure_ = "BMM150 suspended"; ok = false; }
    // Opmode bits 2:1 must be 00 (normal); keep ODR bits 5:3.
    if (ok && (before.mode & 0x06)) { failure_ = "BMM150 not in normal mode"; ok = false; }
    const uint8_t normal = uint8_t(before.mode & 0x38);
    if (ok) {
      // Enter sleep (opmode 11) before changing repetitions.
      failure_ = "AUX write failed";
      ok = auxWrite_(BMM_MODE, uint8_t(normal | 0x06));
      if (ok) wait_(50);
      ok = ok && auxWrite_(BMM_REP_XY, rep_xy) && auxWrite_(BMM_REP_Z, rep_z) && auxWrite_(BMM_MODE, normal);
      if (ok) wait_(50);
    }
    if (ok) {
      ok = readState_(after);
      if (ok && (after.rep_xy != rep_xy || after.rep_z != rep_z || after.mode != normal)) {
        failure_ = "readback mismatch";
        ok = false;
      }
    }
    if (!ok && before.chip_id == 0x32) {
      // Best effort: leave the sensor measuring in normal mode at its old repetitions.
      auxWrite_(BMM_MODE, uint8_t(normal | 0x06));
      wait_(50);
      auxWrite_(BMM_REP_XY, before.rep_xy);
      auxWrite_(BMM_REP_Z, before.rep_z);
      auxWrite_(BMM_MODE, normal);
    }
    const bool left = leave_();
    if (ok && left) failure_ = "none";
    return ok && left;
  }

 private:
  bool rd_(uint8_t reg, uint8_t& v) { return dev_->readRegister(reg, &v, 1); }
  bool wr_(uint8_t reg, uint8_t v) { return dev_->writeRegister8(reg, v); }

  bool waitIdle_() {
    for (int i = 0; i < 20; ++i) {
      uint8_t st = 0;
      if (!rd_(STATUS, st)) return false;
      if (!(st & STATUS_AUX_BUSY)) return true;
      wait_(1);
    }
    return false;
  }

  bool enter_() {
    failure_ = "BMI270 read failed";
    if (!dev_) { failure_ = "no BMI270"; return false; }
    if (!rd_(AUX_IF_CONF, if_conf_) || !rd_(AUX_RD_ADDR, rd_addr_) || !rd_(PWR_CTRL, pwr_) ||
        !rd_(AUX_DEV_ID, dev_id_))
      return false;
    if (dev_id_ != (0x10 << 1)) { failure_ = "AUX device is not at 0x10"; return false; }
    if (if_conf_ & AUX_MANUAL_EN) { failure_ = "AUX already in manual mode"; return false; }
    failure_ = "AUX manual entry failed";
    // Same sequence M5Unified uses: stop AUX data polling, then manual access.
    if (!wr_(PWR_CTRL, uint8_t(pwr_ & ~PWR_AUX_EN))) return false;
    wait_(2);
    entered_ = wr_(AUX_IF_CONF, AUX_MANUAL_EN);
    return entered_;
  }

  bool leave_() {
    if (!entered_) return false;
    entered_ = false;
    // Restore the data-mode configuration exactly, then resume AUX polling.
    const bool ok = wr_(AUX_IF_CONF, if_conf_) && wr_(AUX_RD_ADDR, rd_addr_) && wr_(PWR_CTRL, pwr_);
    uint8_t a = 0, b = 0, c = 0;
    const bool verified = ok && rd_(AUX_IF_CONF, a) && rd_(AUX_RD_ADDR, b) && rd_(PWR_CTRL, c) &&
                          a == if_conf_ && b == rd_addr_ && c == pwr_;
    if (!verified) failure_ = "AUX data-mode restore failed";
    return verified;
  }

  bool auxRead_(uint8_t reg, uint8_t& v) {
    return wr_(AUX_RD_ADDR, reg) && waitIdle_() && rd_(AUX_DATA, v);
  }
  bool auxWrite_(uint8_t reg, uint8_t v) {
    return wr_(AUX_WR_DATA, v) && wr_(AUX_WR_ADDR, reg) && waitIdle_();
  }
  bool readState_(Bmm150RegState& s) {
    failure_ = "AUX read failed";
    return auxRead_(BMM_CHIP_ID, s.chip_id) && auxRead_(BMM_POWER, s.power) && auxRead_(BMM_MODE, s.mode) &&
           auxRead_(BMM_REP_XY, s.rep_xy) && auxRead_(BMM_REP_Z, s.rep_z);
  }

  Device* dev_;
  void (*wait_)(uint32_t);
  uint8_t if_conf_ = 0, rd_addr_ = 0, pwr_ = 0, dev_id_ = 0;
  bool entered_ = false;
  const char* failure_ = "not run";
};

}  // namespace atoms3r_ical
