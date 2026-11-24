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

  void finish(uint32_t sample_us, bool accepted, uint32_t period_us,
              float estimator_dt_s = 0, bool used_host_time = false) {
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
    const uint32_t work = micros() - start_us_;
    ++polls_;
    if (accepted) {
      if (estimator_dt_s > 0) {
        host_steps_ += used_host_time;
        if (!max_dt_s_ || estimator_dt_s < min_dt_s_) min_dt_s_ = estimator_dt_s;
        if (estimator_dt_s > max_dt_s_) max_dt_s_ = estimator_dt_s;
      }
      if (have_previous_) {
        const uint32_t gap = sample_us - previous_us_;
        gap_sum_ += gap;
        ++gaps_;
        if (gap < min_gap_) min_gap_ = gap;
        jitter_sum_ += gap > period_us ? gap - period_us : period_us - gap;
        long_gaps_ += gap > period_us + period_us / 2;
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
      Serial.printf("[TIMING] polls=%u samples=%u sample_gap_min_us=%u sample_gap_mean_us=%u sample_gap_max_us=%u "
                    "jitter_mean_abs_us=%u gaps_over_7_5ms=%u "
                    "accepted_work_mean_us=%u accepted_work_max_us=%u accepted_over_4500us=%u accepted_over_period=%u\n",
                    unsigned(polls_), unsigned(samples_), unsigned(gaps_ ? min_gap_ : 0),
                    unsigned(gaps_ ? (gap_sum_ + gaps_/2)/gaps_ : 0), unsigned(max_gap_), unsigned(gaps_ ? jitter_sum_/gaps_ : 0), unsigned(long_gaps_),
                    unsigned(samples_ ? work_sum_/samples_ : 0), unsigned(max_work_), unsigned(over_target_), unsigned(overruns_));
      if (max_dt_s_ > 0) {
        Serial.printf("[ESTIMATOR] dt_min_ms=%.3f dt_max_ms=%.3f host_steps=%u\n",
                      1000.0 * min_dt_s_, 1000.0 * max_dt_s_, unsigned(host_steps_));
        min_dt_s_ = max_dt_s_ = 0;
        host_steps_ = 0;
      }
      polls_ = samples_ = gaps_ = long_gaps_ = max_gap_ = max_work_ = over_target_ = overruns_ = 0;
      gap_sum_ = jitter_sum_ = work_sum_ = 0;
      min_gap_ = UINT32_MAX;
    }
#else
    (void)sample_us;
    (void)accepted;
    (void)period_us;
    (void)estimator_dt_s;
    (void)used_host_time;
#endif
  }

 private:
#if defined(OCEAN_IMU_PROFILE_LOOP) && OCEAN_IMU_PROFILE_LOOP
  float min_dt_s_ = 0, max_dt_s_ = 0;
  uint32_t host_steps_ = 0;
  uint64_t gap_sum_ = 0, jitter_sum_ = 0, work_sum_ = 0;
  uint32_t start_us_ = 0, previous_us_ = 0, polls_ = 0, samples_ = 0, gaps_ = 0;
  uint32_t long_gaps_ = 0, min_gap_ = UINT32_MAX, max_gap_ = 0, max_work_ = 0, over_target_ = 0, overruns_ = 0;
  bool have_previous_ = false;
#endif
};

}  // namespace atoms3r_ical
