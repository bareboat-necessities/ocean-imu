#pragma once

// Copyright 2026, Mikhail Grushinskiy
#include <cstdint>

namespace atoms3r_ical {

// BMM150 calibration-only high-accuracy preset (Bosch datasheet, table 3).
// M5Unified 0.2.13 selects 30 Hz but leaves repetition registers untouched.
// Device supplies the M5 I2C_Device readRegister/writeRegister8 interface;
// wait_ms lets host tests exercise the same register transaction as the board.
template<class Device>
class MagCalMode {
 public:
  MagCalMode(Device* device, void (*wait_ms)(uint32_t)) : device_(device), wait_(wait_ms) {}
  MagCalMode(const MagCalMode&) = delete;
  MagCalMode& operator=(const MagCalMode&) = delete;
  ~MagCalMode() { restore(); }

  bool begin() {
    if (!device_) return true; // Other IMUs keep their driver settings.
    if (changed_) return false;
    uint8_t id=0, power=0;
    if (!read(0x40,id) || id!=0x32 || !read(0x4B,power) || !(power&1) ||
        !read(0x4C,mode_) || (mode_&0xC7) ||
        !read(0x51,xy_) || !read(0x52,z_)) return false;
    // Enter sleep before changing repetitions. Finish any in-flight conversion.
    changed_=true; // Even a failed acknowledgement may have changed hardware.
    if (!device_->writeRegister8(0x4C,mode_|0x06)) return false;
    wait_(50);
    // 47 XY and 83 Z repetitions; 20 Hz is the maximum for this preset.
    if (!device_->writeRegister8(0x51,23) || !device_->writeRegister8(0x52,82) ||
        !device_->writeRegister8(0x4C,0x28) || !matches(0x28,23,82)) return false;
    wait_(50);
    return true;
  }

  bool restore() {
    if (!changed_) return true;
    if (!device_->writeRegister8(0x4C,0x2E)) return false;
    wait_(50);
    const bool xy_ok=device_->writeRegister8(0x51,xy_);
    const bool z_ok=device_->writeRegister8(0x52,z_);
    // Do not restart at the old ODR with partially restored repetitions.
    if (!xy_ok || !z_ok || !device_->writeRegister8(0x4C,mode_) ||
        !matches(mode_,xy_,z_)) return false;
    changed_=false;
    return true;
  }

 private:
  bool read(uint8_t reg,uint8_t& value) { return device_->readRegister(reg,&value,1); }
  bool matches(uint8_t mode,uint8_t xy,uint8_t z) {
    uint8_t m=0,x=0,y=0;
    return read(0x4C,m) && read(0x51,x) && read(0x52,y) && m==mode && x==xy && y==z;
  }
  Device* device_;
  void (*wait_)(uint32_t);
  uint8_t mode_=0,xy_=0,z_=0;
  bool changed_=false;
};

} // namespace atoms3r_ical
