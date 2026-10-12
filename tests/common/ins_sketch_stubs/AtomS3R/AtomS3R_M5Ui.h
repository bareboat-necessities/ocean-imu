#pragma once
// Host stand-in for the AtomS3R text UI.
#include <M5Unified.h>
namespace atoms3r_ical {
struct M5UiCfg { static constexpr uint32_t MENU_TAP_WINDOW_MS = 600; };
struct M5Ui {
  template <class... T> void begin(T...) {}
  template <class... T> void setReadRotation(T...) {}
  template <class... T> void title(T...) {}
  template <class... T> void line(T...) {}
  template <class... T> void lineAt(T...) {}
  template <class... T> void bar01(T...) {}
  template <class... T> void fail(T...) {}
  template <class... T> void notSavedNotice(T...) {}
  template <class... T> bool eraseConfirm(T...) { return false; }
};
}  // namespace atoms3r_ical
