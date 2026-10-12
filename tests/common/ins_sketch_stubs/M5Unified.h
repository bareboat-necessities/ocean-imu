#pragma once
// Host stand-in for M5Unified: a display that accepts drawing calls and an IMU
// whose update() reports a new accelerometer+gyro frame on every loop tick.
// A test can install a simulated BMI270 register device (m5::IMU_Base) for
// sketches that talk to the sensor directly, and give a full-screen fill a
// host-time cost.
#include <Arduino.h>

namespace m5 {
struct IMU_Base {
  virtual ~IMU_Base() = default;
  virtual bool readRegister(uint8_t reg, uint8_t* buf, size_t len) = 0;
  virtual bool writeRegister8(uint8_t reg, uint8_t value) = 0;
};
enum imu_t { imu_none, imu_bmi270 };
}  // namespace m5

constexpr int TFT_BLACK = 0, TFT_WHITE = 1, TFT_RED = 2, TFT_GREEN = 3, TFT_BLUE = 4,
              TFT_YELLOW = 5, TFT_ORANGE = 6, TFT_DARKGREY = 7, TFT_LIGHTGREY = 8,
              TFT_CYAN = 9, TFT_MAGENTA = 10;

struct M5DisplayStub : Print {
  template <class... T> void setCursor(T...) {}
  template <class... T> void setTextSize(T...) {}
  template <class... T> void setTextColor(T...) {}
  template <class... T> void setRotation(T...) {}
  static uint32_t& fillCostUs() { static uint32_t us = 0; return us; }
  template <class... T> void fillScreen(T...) { ins_sketch_host::clock_us() += fillCostUs(); }
  template <class... T> void fillRect(T...) {}
  template <class... T> void drawRect(T...) {}
  template <class... T> void setTextDatum(T...) {}
  template <class... T> void drawString(T...) {}
  template <class... T> void startWrite(T...) {}
  template <class... T> void endWrite(T...) {}
  int width() const { return 128; }
  int height() const { return 128; }
  int fontHeight() const { return 8; }
  int textWidth(const char* s) const { return int(std::strlen(s)) * 6; }
};
struct M5ImuStub {
  m5::IMU_Base* device = nullptr;
  bool isEnabled() const { return true; }
  uint32_t update() { return 3u; }  // accelerometer and gyro bits
  m5::IMU_Base* getImuInstancePtr(int) const { return device; }
  m5::imu_t getType() const { return device ? m5::imu_bmi270 : m5::imu_none; }
};
struct M5Config { bool internal_imu = false; };
struct M5Stub {
  M5DisplayStub Display;
  M5ImuStub Imu;
  M5Config config() const { return {}; }
  void begin(const M5Config&) {}
  void update() {}
};
inline M5Stub M5;
