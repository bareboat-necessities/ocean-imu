#pragma once

#include <Arduino.h>

namespace atoms3r_ical {

// Measure active poll-to-output work separately from accepted sample spacing.
// With profiling disabled neither hook reads the clock or writes to Serial.
class LoopTiming {
 public:
  void begin() {
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
    start_us_ = micros();
#endif
  }

  void finish(uint32_t sample_us, bool accepted, uint32_t period_us) {
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
    const uint32_t work = micros() - start_us_;
    ++polls_;
    if (accepted) {
      if (have_previous_) {
        const uint32_t gap = sample_us - previous_us_;
        gap_sum_ += gap;
        ++gaps_;
        if (gap > max_gap_) max_gap_ = gap;
      }
      previous_us_ = sample_us;
      have_previous_ = true;
      ++samples_;
      work_sum_ += work;
      if (work > max_work_) max_work_ = work;
      over_target_ += work > 4500;
      overruns_ += work > period_us;
    }
    if (polls_ == 1000) {
      Serial.printf("[TIMING] polls=%u samples=%u sample_gap_mean_us=%u sample_gap_max_us=%u "
                    "accepted_work_mean_us=%u accepted_work_max_us=%u accepted_over_4500us=%u accepted_over_period=%u\n",
                    unsigned(polls_), unsigned(samples_), unsigned(gaps_ ? gap_sum_/gaps_ : 0), unsigned(max_gap_),
                    unsigned(samples_ ? work_sum_/samples_ : 0), unsigned(max_work_), unsigned(over_target_), unsigned(overruns_));
      polls_ = samples_ = gaps_ = max_gap_ = max_work_ = over_target_ = overruns_ = 0;
      gap_sum_ = work_sum_ = 0;
    }
#else
    (void)sample_us;
    (void)accepted;
    (void)period_us;
#endif
  }

 private:
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
  uint64_t gap_sum_ = 0, work_sum_ = 0;
  uint32_t start_us_ = 0, previous_us_ = 0, polls_ = 0, samples_ = 0, gaps_ = 0;
  uint32_t max_gap_ = 0, max_work_ = 0, over_target_ = 0, overruns_ = 0;
  bool have_previous_ = false;
#endif
};

}  // namespace atoms3r_ical
