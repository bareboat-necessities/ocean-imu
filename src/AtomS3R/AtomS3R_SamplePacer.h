#pragma once

#include <Arduino.h>

namespace atoms3r_ical {

// Pace host polling from the previous scheduled release. Keep sensor timestamps
// at the actual read and let a late loop poll immediately.
class SamplePacer {
 public:
  explicit SamplePacer(uint32_t period_us) : period_us_(period_us) {}

  void wait() {
    uint32_t now = micros();
    if (started_) {
      while (now - phase_us_ < period_us_) now = micros();
      const uint32_t due = phase_us_ + period_us_;
      // A small wake-up delay is absorbed by the next wait. At half a
      // period late we are closer to the next slot: start a new phase so
      // there is no closely spaced catch-up poll after a stall.
      if (now - due < period_us_ / 2) {
        phase_us_ = due;
        return;
      }
    }
    phase_us_ = now;
    started_ = true;
  }

  // Only a successful DRDY retry is evidence that the host slot was early.
  // Anchor the next wait to that poll so a slightly slower sensor clock can
  // move the release phase without eventually losing a whole sample.
  void followReadyPoll(uint32_t first_us, uint32_t ready_us) {
    if (ready_us != first_us) phase_us_ = ready_us;
  }

 private:
  const uint32_t period_us_;
  uint32_t phase_us_ = 0;
  bool started_ = false;
};

}  // namespace atoms3r_ical
