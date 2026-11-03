#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#include <Arduino.h>
#include <M5Unified.h>
#include <math.h>
#include "imu_calibrate/FrameAxes.h"

namespace atoms3r_ical {

// Standalone UI config
struct M5UiCfg {
  // If text rotated wrong initially, change this to 1 or 3.
  static constexpr uint8_t  ROT_READ = 0;
  static constexpr uint8_t  LCD_BRIGHTNESS = 200;

  // Keep-awake / pacing defaults (wizard can live with these)
  static constexpr uint32_t KEEP_AWAKE_EVERY_MS = 200;
  static constexpr uint32_t OK_PAUSE_MS         = 900;
  static constexpr uint32_t MENU_TAP_WINDOW_MS  = 650;
};

static inline float clamp01_(float x) { return x < 0 ? 0 : (x > 1 ? 1 : x); }

// Input: BtnA ONLY + keep-awake (self-contained)
class Input {
public:
  static void update() {
    M5.update();

    const uint32_t now = millis();
    if ((uint32_t)(now - last_keep_awake_ms_) >= M5UiCfg::KEEP_AWAKE_EVERY_MS) {
      last_keep_awake_ms_ = now;
      M5.Display.setBrightness(M5UiCfg::LCD_BRIGHTNESS);
      M5.Display.wakeup();
    }
    tap_edge_ = M5.BtnA.wasPressed();
  }

  static bool tapPressed() { return tap_edge_; }

private:
  static inline bool tap_edge_ = false;
  static inline uint32_t last_keep_awake_ms_ = 0;
};

// UI helpers
class M5Ui {
public:
  void begin() {
    rot_ = M5UiCfg::ROT_READ;
    M5.Display.setRotation(rot_);
    M5.Display.setBrightness(M5UiCfg::LCD_BRIGHTNESS);
    M5.Display.setTextSize(1);
  }

  void setRotation(uint8_t r) {
    rot_ = (uint8_t)(r & 3);
    M5.Display.setRotation(rot_);
  }

  void setReadRotation() { setRotation(M5UiCfg::ROT_READ); }
  uint8_t rotation() const { return rot_; }

  void clear() {
    M5.Display.fillScreen(TFT_BLACK);
    M5.Display.setCursor(0, 0);
  }

  void title(const char* t) {
    clear();
    M5.Display.setTextColor(TFT_WHITE, TFT_BLACK);
    M5.Display.setTextSize(2);
    M5.Display.println(t);
    M5.Display.setTextSize(1);
    M5.Display.println();
  }

  void line(const char* s = "") {
    M5.Display.setTextColor(TFT_WHITE, TFT_BLACK);
    M5.Display.println(s);
  }

  void bar01(float t01) {
    t01 = clamp01_(t01);
    const int x = 6;
    const int w = M5.Display.width() - 12;
    const int h = 10;
    const int y = M5.Display.height() - 18;
    M5.Display.drawRect(x, y, w, h, TFT_DARKGREY);
    const int fillw = (int)((w - 2) * t01);
    M5.Display.fillRect(x + 1, y + 1, fillw, h - 2, TFT_GREEN);
    // Motion can reset capture progress; erase the previous longer fill.
    M5.Display.fillRect(x + 1 + fillw, y + 1, w - 2 - fillw, h - 2, TFT_BLACK);
  }

  void waitTap(const char* t, const char* l1 = nullptr, const char* l2 = nullptr) {
    setReadRotation();
    waitTapHere(t,l1,l2);
  }

  // Keep the capture pose's readable orientation for its next/retry screen.
  void waitTapHere(const char* t, const char* l1 = nullptr, const char* l2 = nullptr) {
    title(t);
    if (l1) line(l1);
    if (l2) line(l2);
    line("");
    line("Tap BtnA");
    while (true) {
      Input::update();
      if (Input::tapPressed()) break;
      delay(10);
    }
  }

  void showOkAuto(const char* l1 = nullptr, const char* l2 = nullptr, const char* heading = "OK") {
    title(heading);
    if (l1) line(l1);
    if (l2) line(l2);
    const uint32_t t0 = millis();
    while (millis() - t0 < M5UiCfg::OK_PAUSE_MS) {
      Input::update();
      if (Input::tapPressed()) break;
      delay(10);
    }
  }

  void fail(const char* where, const char* why) {
    setReadRotation();
    title("FAILED");
    line(where);
    line(why);
    line("");
    line("Tap BtnA");
    while (true) {
      Input::update();
      if (Input::tapPressed()) break;
      delay(10);
    }
  }

  bool eraseConfirm() {
    setReadRotation();
    title("ERASE?");
    line("Delete saved cal");
    line("Tap=YES  Wait=NO");
    const uint32_t t0 = millis();
    while (millis() - t0 < 4500) {
      Input::update();
      if (Input::tapPressed()) return true;
      delay(10);
    }
    return false;
  }

  enum class MagFailAction : uint8_t { RETRY_MAG = 0, REDO_ALL = 1, ABORT = 2 };

  // No timeout (tap-group UI)
  MagFailAction magFailMenu(const char* why1, const char* why2 = nullptr) {
    setReadRotation();
    title("MAG FAIL");
    if (why1) line(why1);
    if (why2) line(why2);
    line("");
    line("Tap: retry MAG");
    line("Tap x2: redo ALL");
    line("Tap x3: abort");

    const uint8_t taps = waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
    if (taps >= 3) return MagFailAction::ABORT;
    if (taps == 2) return MagFailAction::REDO_ALL;
    return MagFailAction::RETRY_MAG;
  }

  enum class TimingFailAction : uint8_t { RETRY = 0, SKIP = 1, ABORT = 2 };

  // Magnetic timing is optional: skipping saves the calibration without it.
  TimingFailAction timingFailMenu(const char* why1, const char* why2 = nullptr) {
    setReadRotation();
    title("TIME FAIL");
    if (why1) line(why1);
    if (why2) line(why2);
    line("");
    line("Tap: retry TIME");
    line("Tap x2: skip TIME");
    line("Tap x3: abort");

    const uint8_t taps = waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
    if (taps >= 3) return TimingFailAction::ABORT;
    if (taps == 2) return TimingFailAction::SKIP;
    return TimingFailAction::RETRY;
  }

  // Failure screen with up to three short lines.
  void failLines(const char* where, const char* why1, const char* why2 = nullptr, const char* why3 = nullptr) {
    setReadRotation();
    title("FAILED");
    line(where);
    if (why1) line(why1);
    if (why2) line(why2);
    if (why3) line(why3);
    line("");
    line("Tap BtnA");
    while (true) {
      Input::update();
      if (Input::tapPressed()) break;
      delay(10);
    }
  }

  // Current-record tools operate without repeating the full capture.
  enum class StartAction : uint8_t { FULL=0, NOISE=1, CANCEL=2, GYRO_BIAS=3 };

  StartAction startMenu(bool saved_available) {
    for(;;) {
      setReadRotation();title("CALIBRATE");line("Tap: full calibration");
      if(saved_available)line("2 taps: saved fits");
      line("3 taps: cancel");
      const auto taps=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
      if(taps>=3)return StartAction::CANCEL;
      if(taps!=2 || !saved_available)return StartAction::FULL;
      title("SAVED FITS");line("Keep accel and mag");line("");
      line("Tap: gyro bias only");line("2 taps: noise 40m");line("3 taps: back");
      const auto tool=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
      if(tool==1)return StartAction::GYRO_BIAS;
      if(tool==2)return StartAction::NOISE;
    }
  }

  // Retry prompt for one hold; false = abort (previous calibration kept).
  bool retryMenu(const char* what, const char* why1, const char* why2 = nullptr) {
    setReadRotation();
    return retryMenuHere(what,why1,why2);
  }

  bool retryMenuHere(const char* what, const char* why1, const char* why2 = nullptr) {
    title("RETRY?");
    if (what) line(what);
    if (why1) line(why1);
    if (why2) line(why2);
    line("");
    line("Tap: retry");
    line("Tap x3: abort");
    const uint8_t taps = waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
    return taps < 3;
  }

  bool confirmDiscardCaptures() {
    title("DISCARD?");line("Unsaved data lost");line("");
    line("Tap: keep data");line("3 taps: discard");
    return waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS)>=3;
  }

  enum class GyroFitAction : uint8_t { RETRY_FIT, ADD_TURNS, BIAS_ONLY, ABORT };
  GyroFitAction gyroTiltMenu() {
    title("GYRO TILT");line("Optional scale / axes");line("Bias already measured");line("");
    line("Tap: 12 turns");line("2 taps: use bias only");line("3 taps: exit menu");
    const uint8_t taps=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
    if(taps>=3)return confirmDiscardCaptures()?GyroFitAction::ABORT:GyroFitAction::RETRY_FIT;
    return taps==2?GyroFitAction::BIAS_ONLY:GyroFitAction::RETRY_FIT;
  }
  GyroFitAction gyroFitFailMenu(const char* why,bool can_add=false) {
    title("GYRO FIT");line(why);line("Captures kept");line("");
    line(can_add?"Tap: add 4 turns":"Tap: retry fit");line("2 taps: use bias only");line("3 taps: exit menu");
    const uint8_t taps=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
    if(taps>=3)return confirmDiscardCaptures()?GyroFitAction::ABORT:GyroFitAction::RETRY_FIT;
    return taps==2?GyroFitAction::BIAS_ONLY:can_add?GyroFitAction::ADD_TURNS:GyroFitAction::RETRY_FIT;
  }
  GyroFitAction gyroTurnFailMenu(const char* why,const char* kept) {
    title("GYRO TURN");line(why);line(kept);line("");
    line("Tap: retry turn");line("2 taps: use bias only");line("3 taps: exit menu");
    const uint8_t taps=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
    if(taps>=3)return confirmDiscardCaptures()?GyroFitAction::ABORT:GyroFitAction::RETRY_FIT;
    return taps==2?GyroFitAction::BIAS_ONLY:GyroFitAction::RETRY_FIT;
  }

  struct MountAxes {
    imu_cal::FrameAxis up=imu_cal::FrameAxis::NEG_Z;
    imu_cal::FrameAxis forward=imu_cal::FrameAxis::POS_X;
  };
  MountAxes mountingMenu() {
    MountAxes choice;
    for(;;) {
      setReadRotation();title("MOUNT");line("When installed:");line("Screen faces up");line("USB points backward");line("");
      line("Tap: use this mount");line("2 taps: change mount");line("3 taps: help");
      const auto taps=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
      if(taps==1)return MountAxes{};
      if(taps>=3) {waitTapHere("MOUNT HELP","Sets installation","No pose to capture");continue;}
      waitTapHere("CASE EDGES","Read screen upright","USB at bottom");
      if(editMounting_(choice))return choice;
    }
  }

  // Overwrites one text row (size-1 font) at the current rotation.
  void lineAt(int row_px, const char* s) {
    M5.Display.fillRect(0, row_px, M5.Display.width(), 8, TFT_BLACK);
    M5.Display.setCursor(0, row_px);
    M5.Display.setTextColor(TFT_WHITE, TFT_BLACK);
    M5.Display.print(s);
  }

  int cursorY() const { return M5.Display.getCursorY(); }

  void notSavedNotice() {
    setReadRotation();
    title("NOT SAVED");
    line("Calibration not saved");
    line("See Serial log");
    line("");
    line("Tap BtnA");
    while (true) {
      Input::update();
      if (Input::tapPressed()) break;
      delay(10);
    }
  }

private:
  static const char* mountFaceName_(imu_cal::FrameAxis face) {
    switch(face) {
      case imu_cal::FrameAxis::NEG_Z:return "Screen";
      case imu_cal::FrameAxis::POS_Z:return "Back of case";
      case imu_cal::FrameAxis::NEG_X:return "USB end";
      case imu_cal::FrameAxis::POS_X:return "Opposite USB end";
      case imu_cal::FrameAxis::NEG_Y:return "Left edge";
      default:return "Right edge";
    }
  }
  bool chooseMountFace_(bool forward,MountAxes& choice) {
    auto& face=forward?choice.forward:choice.up;
    if(forward && !imu_cal::frameAxesPerpendicular(choice.up,face)) {
      face=imu_cal::FrameAxis::NEG_Z;
      while(!imu_cal::frameAxesPerpendicular(choice.up,face))face=imu_cal::FrameAxis(uint8_t(face)+1);
    }
    for(;;) {
      title(forward?"POINTS FWD":"FACE UP");
      line(forward?"Toward front of boat":"When mounted:");
      if(!forward)line("90-degree mounts");
      line(mountFaceName_(face));line("");
      line("Tap: next face");line("2 taps: choose");line("3 taps: back");
      const auto taps=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
      if(taps>=3)return false;
      if(taps==2)return true;
      do {face=imu_cal::FrameAxis((uint8_t(face)+1)%6);}
      while(forward && !imu_cal::frameAxesPerpendicular(choice.up,face));
    }
  }
  bool editMounting_(MountAxes& choice) {
    for(;;) {
      if(!chooseMountFace_(false,choice))return false;
      if(!chooseMountFace_(true,choice))continue;
      title("MOUNT OK?");line("Up:");line(mountFaceName_(choice.up));
      line("Forward:");line(mountFaceName_(choice.forward));line("");
      line("Tap: use this mount");line("2 taps: edit");line("3 taps: back");
      const auto taps=waitTapGroupNoTimeout_(M5UiCfg::MENU_TAP_WINDOW_MS);
      if(taps==1)return true;
      if(taps>=3)return false;
    }
  }

  static uint8_t waitTapGroupNoTimeout_(uint32_t window_ms) {
    uint8_t count = 0;
    uint32_t deadline = 0;

    // wait for first tap
    while (true) {
      Input::update();
      if (Input::tapPressed()) {
        count = 1;
        deadline = millis() + window_ms;
        break;
      }
      delay(10);
    }

    // collect taps until window expires
    while (true) {
      Input::update();
      if (Input::tapPressed()) {
        ++count;
        deadline = millis() + window_ms;
      }
      if ((int32_t)(millis() - deadline) > 0) break;
      delay(10);
    }
    return count;
  }

  uint8_t rot_ = M5UiCfg::ROT_READ;
};

} // namespace atoms3r_ical
