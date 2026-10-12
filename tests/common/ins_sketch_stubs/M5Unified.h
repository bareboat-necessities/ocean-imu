#pragma once
// Host stand-in for M5Unified: a display that accepts drawing calls and an IMU
// whose update() reports a new accelerometer+gyro frame on every loop tick.
#include <Arduino.h>

constexpr int TFT_BLACK = 0, TFT_WHITE = 1, TFT_RED = 2, TFT_GREEN = 3, TFT_BLUE = 4,
              TFT_YELLOW = 5, TFT_ORANGE = 6, TFT_DARKGREY = 7, TFT_LIGHTGREY = 8,
              TFT_CYAN = 9, TFT_MAGENTA = 10;

struct M5DisplayStub : Print {
  template <class... T> void setCursor(T...) {}
  template <class... T> void setTextSize(T...) {}
  template <class... T> void setTextColor(T...) {}
  template <class... T> void setRotation(T...) {}
  template <class... T> void fillScreen(T...) {}
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
  bool isEnabled() const { return true; }
  uint32_t update() { return 3u; }  // accelerometer and gyro bits
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
