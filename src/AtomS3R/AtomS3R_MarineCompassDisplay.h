#pragma once

#include "AtomS3R/AtomS3R_CompassUI.h"
#include "AtomS3R/AtomS3R_M5Ui.h"
#include <freertos/FreeRTOS.h>
#include <freertos/queue.h>
#include <freertos/semphr.h>
#include <freertos/task.h>
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
#include <atomic>
#endif

namespace atoms3r_ical {

// Frames cross cores by value. The renderer has access only to its saved
// frame and the display. Menus synchronize through waitForMenu().
template <class Frame>
class CompassDisplayTask {
public:
  using Renderer = void (*)(const Frame&, CompassUI&, M5Ui&);
  CompassDisplayTask(M5Ui& ui, Renderer render) : ui_(ui), render_(render) {}

  void begin(bool graphics) {
    if (graphics) compass_.begin();
    queue_ = xQueueCreateStatic(1, sizeof(Message), queueBytes_, &queueStorage_);
    idle_ = xSemaphoreCreateBinaryStatic(&idleStorage_);
    // Arduino runs setup and loop on the same core on the AtomS3R.
    const BaseType_t core = portNUM_PROCESSORS > 1 ? 1 - xPortGetCoreID() : 0;
    if (queue_ && idle_) {
      worker_ = xTaskCreateStaticPinnedToCore(run, "compass", sizeof(stack_), this, tskIDLE_PRIORITY + 1,
                                               stack_, &taskStorage_, core);
    }
    if (!worker_) Serial.println("[LCD] Display task unavailable, using synchronous redraw");
  }

  void publish(const Frame& frame) {
    if (!worker_) {
      draw(frame);
      return;
    }
    const Message msg{frame, false};
    xQueueOverwrite(queue_, &msg);
  }

  // Called by the sole producer before a menu or wizard touches the LCD.
  // Replace any queued frame and wait for the in-flight transfer to finish.
  void waitForMenu() {
    if (!worker_) return;
    const Message msg{{}, true};
    xQueueOverwrite(queue_, &msg);
    xSemaphoreTake(idle_, portMAX_DELAY);
  }

#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
  uint32_t completedFrames() const { return completed_.load(std::memory_order_relaxed); }
  uint32_t stackFreeMinBytes() const { return stackFree_.load(std::memory_order_relaxed); }
#endif

private:
  struct Message {
    Frame frame;
    bool menu;
  };

  static void run(void* context) {
    auto& self = *static_cast<CompassDisplayTask*>(context);
    for (;;) {
      Message msg;
      if (xQueueReceive(self.queue_, &msg, portMAX_DELAY) != pdTRUE) continue;
      if (msg.menu) xSemaphoreGive(self.idle_);
      else {
        self.draw(msg.frame);
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
        self.stackFree_.store(uxTaskGetStackHighWaterMark(nullptr) * sizeof(StackType_t),
                              std::memory_order_relaxed);
#endif
      }
    }
  }

  void draw(const Frame& frame) {
    render_(frame, compass_, ui_);
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
    ++completed_;
#endif
  }

  M5Ui& ui_;
  Renderer render_;
  CompassUI compass_;
  StaticTask_t taskStorage_;
  // Match the default Arduino-ESP32 loop stack for M5GFX and float printf.
  alignas(portBYTE_ALIGNMENT) StackType_t stack_[8192 / sizeof(StackType_t)];
  StaticSemaphore_t idleStorage_;
  StaticQueue_t queueStorage_;
  uint8_t queueBytes_[sizeof(Message)];
  QueueHandle_t queue_{};
  SemaphoreHandle_t idle_{};
  TaskHandle_t worker_{};
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
  std::atomic<uint32_t> completed_{0}, stackFree_{0};
#endif
};

struct MarineCompassFrame {
  float heading, magneticHeading, roll, pitch, heave, waveDirection;
  float magNorm, rawAccelNorm, accelNorm;
  int waveSign;
  bool headingValid, magneticHeadingValid, magOk, magFresh;
};

class MarineCompassDisplay : public CompassDisplayTask<MarineCompassFrame> {
public:
  using Frame = MarineCompassFrame;
  explicit MarineCompassDisplay(M5Ui& ui) : CompassDisplayTask(ui, drawFrame) {}

private:
  static void drawFrame(const Frame& frame, CompassUI& compass, M5Ui& ui) {
    if (compass.ok()) {
      const bool tilt = fabsf(frame.roll) > 35.0f || fabsf(frame.pitch) > 35.0f;
      compass.draw(frame.headingValid ? frame.heading : 0.0f,
                    frame.headingValid && frame.magOk, frame.magNorm, tilt);
    } else {
      drawText(frame, ui);
    }
  }

  static void drawText(const Frame& frame, M5Ui& ui) {
    ui.title("COMPASS");
    if (frame.headingValid) M5.Display.printf("HDG: %6.1f deg\n", frame.heading);
    else                    M5.Display.printf("HDG:   --- WAIT\n");
    M5.Display.printf("HDM: %6.1f %s\n", frame.magneticHeading,
                      frame.magneticHeadingValid ? "MAG" : "---");
    M5.Display.printf("ROL: %6.1f deg\n", frame.roll);
    M5.Display.printf("PIT: %6.1f deg\n", frame.pitch);
    M5.Display.printf("HEV: %6.3f m\n", frame.heave);
    M5.Display.printf("WAV: %6.1f s=%d\n", frame.waveDirection, frame.waveSign);
    M5.Display.printf("MAG: %s %s\n", frame.magOk ? "OK " : "BAD", frame.magFresh ? "NEW" : "OLD");
    M5.Display.printf("|m|: %6.1f uT\n", frame.magNorm);
    M5.Display.printf("|aR|:%5.2f |aC|:%5.2f\n", frame.rawAccelNorm, frame.accelNorm);
    ui.line("");
  }
};

}  // namespace atoms3r_ical
