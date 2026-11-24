#pragma once

#include <stdint.h>

namespace atoms3r_ical {

// BMI270 datasheet section 4.6.15: at 200 Hz each data-register update
// toggles bit 7 of SENSORTIME. The upper 17 bits identify the frame.
class Bmi270FrameClock {
 public:
  void reset() { started_ = false; }

  bool advance(uint32_t poll_us, int32_t frame, float& dt) {
    constexpr uint32_t period_us = 5000;
    constexpr uint32_t frame_modulus = 1u << 17;
    const uint32_t elapsed = poll_us - last_poll_us_;
    uint32_t step_us = started_ && elapsed ? elapsed : period_us;
    if (started_ && frame >= 0 && last_frame_ >= 0 &&
        elapsed < frame_modulus * period_us / 2) {
      const uint32_t frames = uint32_t(frame - last_frame_) & (frame_modulus - 1);
      if (frames == 0) return false;
      const uint32_t sensor_us = frames * period_us;
      const uint32_t difference = sensor_us > elapsed ? sensor_us - elapsed : elapsed - sensor_us;
      // A reset can look like a huge forward jump modulo the counter.
      // Allow one frame of read phase and generous sensor/host clock drift.
      if (difference <= period_us + elapsed / 8) step_us = sensor_us;
    }
    started_ = true;
    last_poll_us_ = poll_us;
    last_frame_ = frame;
    dt = float(step_us) * 1.0e-6f;
    return true;
  }

 private:
  bool started_ = false;
  int32_t last_frame_ = -1;
  uint32_t last_poll_us_ = 0;
};

}  // namespace atoms3r_ical
