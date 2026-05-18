#pragma once

#include <math.h>
#include <stdint.h>

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

class AdaptiveWaveDetrender3D {
public:
  using Vec3 = Eigen::Vector3f;

  struct Config {
    // Dominant wave-frequency search band [Hz].
    // These bounds limit both internally learned and externally supplied frequency.
    float init_wave_freq_hz = 0.12f;   // startup guess (~8.3 s period)
    float min_wave_freq_hz  = 0.02f;   // 50 s period
    float max_wave_freq_hz  = 1.20f;   // 0.83 s period

    // Slow baseline cutoff:
    //   f_baseline = clamp(baseline_cutoff_fraction * f_wave,
    //                      min_baseline_cutoff_hz,
    //                      max_baseline_cutoff_hz)
    //
    // Larger fraction -> stronger drift removal, but more risk of attenuating wave content.
    float baseline_cutoff_fraction = 0.25f;
    float min_baseline_cutoff_hz   = 0.003f;
    float max_baseline_cutoff_hz   = 0.25f;

    // Smoothing time constant for the learned wave frequency [s].
    // Larger -> slower / more stable frequency adaptation.
    float freq_smooth_tau_s = 12.0f;

    // Internal slope-processing for frequency learning.
    // Frequency is learned from a filtered derivative of the detrended signal.
    float slope_lpf_tau_s = 0.20f;
    float slope_rms_tau_s = 8.0f;

    // Schmitt threshold for zero/level-crossing-based period detection.
    // Threshold is proportional to slope RMS, then clamped to absolute bounds.
    float threshold_rms_fraction  = 0.30f;
    float min_slope_threshold_abs = 0.001f;
    float max_slope_threshold_abs = 1.0e9f;

    // Wait this long before internal frequency learning starts [s].
    // Lets the detrending state settle a bit after startup or re-seed.
    float startup_hold_s = 2.0f;

    // Frequency validity timeout, expressed in cycles of the current frequency.
    // If no fresh internal/external update arrives within this time, freq_valid becomes false.
    float freq_timeout_cycles = 3.0f;

    // Which axis should drive internal frequency learning?
    //   -1 = auto: use axis with largest slope RMS
    //    0 = x
    //    1 = y
    //    2 = z
    int8_t freq_learn_axis = -1;

    // Extra cleanup high-pass on wave_raw only.
    // This does NOT change the baseline being subtracted.
    //
    // cleanup cutoff:
    //   f_cleanup = clamp(cleanup_cutoff_fraction * f_baseline,
    //                     min_cleanup_cutoff_hz,
    //                     max_cleanup_cutoff_hz)
    //
    // cleanup_stages = 0 -> disabled
    // cleanup_stages = 1 -> one extra HP stage
    // cleanup_stages = 2,3 -> stronger low-frequency suppression
    bool    enable_wave_cleanup = true;
    float   cleanup_cutoff_fraction = 1.0f;
    float   min_cleanup_cutoff_hz   = 0.003f;
    float   max_cleanup_cutoff_hz   = 0.50f;
    uint8_t cleanup_stages          = 1;

    // dt guards.
    // min_dt_s protects divisions and exponentials.
    // max_dt_s is also used as a discontinuity threshold.
    float min_dt_s = 1.0e-4f;
    float max_dt_s = 2.0f;

    // Optional hard clamp on wave outputs per component. <= 0 disables.
    float output_abs_limit = 0.0f;
  };

  struct Output {
    // Latest input sample.
    Vec3 input = Vec3::Zero();

    // Actual slow baseline being subtracted.
    Vec3 baseline_slow = Vec3::Zero();

    // Raw detrended signal = input - baseline_slow.
    Vec3 wave_raw = Vec3::Zero();

    // Additional cleanup HP output applied only to wave_raw.
    Vec3 wave_clean = Vec3::Zero();

    // Shared learned dominant wave frequency for the full 3D signal.
    float wave_freq_hz = 0.0f;
    float wave_period_s = 0.0f;

    // Current adaptive baseline filter cutoff and equivalent tau.
    float baseline_cutoff_hz = 0.0f;
    float baseline_tau_s = 0.0f;

    // Current cleanup cutoff (0 if cleanup disabled).
    float cleanup_cutoff_hz = 0.0f;

    // Per-axis slope RMS and selected-axis RMS used for frequency learning.
    Vec3  slope_rms_xyz = Vec3::Zero();
    float slope_rms_selected = 0.0f;
    float slope_threshold = 0.0f;

    // Whether the learned/shared frequency is currently considered fresh.
    bool  freq_valid = false;

    // Current Schmitt state of the learning channel: -1, 0, +1.
    int8_t schmitt_state = 0;

    // Which axis was used for learning this step.
    int8_t selected_axis = -1;
  };

  AdaptiveWaveDetrender3D() {
    sanitizeConfig_(cfg_);
    reset();
  }

  explicit AdaptiveWaveDetrender3D(const Config& cfg) : cfg_(cfg) {
    sanitizeConfig_(cfg_);
    reset();
  }

  void setConfig(const Config& cfg) {
    cfg_ = cfg;
    sanitizeConfig_(cfg_);

    // Keep the currently held frequency inside the new allowed band.
    f_used_hz_ = clampf_(f_used_hz_, cfg_.min_wave_freq_hz, cfg_.max_wave_freq_hz);

    // Refresh cached output using the new config.
    last_output_ = buildOutput_();
  }

  const Config& config() const { return cfg_; }

  void reset() {
    initialized_ = false;
    time_s_ = 0.0f;

    // Frequency learning is disabled until this time.
    learning_enabled_after_t_ = cfg_.startup_hold_s;

    current_input_.setZero();
    baseline_slow_.setZero();
    wave_learn_prev_.setZero();

    slope_filt_.setZero();
    slope_prev_.setZero();
    slope_rms2_.setZero();

    // Start with configured initial guess.
    f_used_hz_ = cfg_.init_wave_freq_hz;

    schmitt_state_ = 0;
    selected_axis_ =
        (cfg_.freq_learn_axis >= 0 && cfg_.freq_learn_axis <= 2)
            ? cfg_.freq_learn_axis
            : -1;

    // No crossings seen yet.
    have_last_pos_cross_ = false;
    have_last_neg_cross_ = false;
    last_pos_cross_t_ = 0.0f;
    last_neg_cross_t_ = 0.0f;

    // No valid frequency measurements yet.
    valid_period_count_ = 0;
    last_valid_freq_t_ = -1.0e30f;
    last_external_freq_t_ = -1.0e30f;

    last_wave_raw_.setZero();
    last_wave_clean_.setZero();

    zeroCleanupStates_();
    last_output_ = buildOutput_();
  }

  void reset(const Vec3& x0) {
    sanitizeConfig_(cfg_);

    initialized_ = true;
    time_s_ = 0.0f;
    learning_enabled_after_t_ = cfg_.startup_hold_s;

    current_input_ = x0;

    // Seed the slow baseline to the initial input so raw wave starts near zero.
    baseline_slow_ = x0;

    // Frequency-learning signal history starts from zero.
    wave_learn_prev_.setZero();

    slope_filt_.setZero();
    slope_prev_.setZero();
    slope_rms2_.setZero();

    f_used_hz_ = cfg_.init_wave_freq_hz;

    schmitt_state_ = 0;
    selected_axis_ =
        (cfg_.freq_learn_axis >= 0 && cfg_.freq_learn_axis <= 2)
            ? cfg_.freq_learn_axis
            : -1;

    have_last_pos_cross_ = false;
    have_last_neg_cross_ = false;
    last_pos_cross_t_ = 0.0f;
    last_neg_cross_t_ = 0.0f;

    valid_period_count_ = 0;
    last_valid_freq_t_ = -1.0e30f;
    last_external_freq_t_ = -1.0e30f;

    last_wave_raw_.setZero();
    last_wave_clean_.setZero();

    zeroCleanupStates_();
    last_output_ = buildOutput_();
  }

  void reset(float x0, float y0, float z0) {
    reset(Vec3(x0, y0, z0));
  }

  // Update without external frequency guidance.
  Output update(const Vec3& x, float dt_s) {
    return updateImpl_(x, dt_s, 0.0f, false);
  }

  Output update(float x, float y, float z, float dt_s) {
    return updateImpl_(Vec3(x, y, z), dt_s, 0.0f, false);
  }

  // Update with optional external frequency guidance.
  // external_valid=true means external_wave_freq_hz is trusted this step.
  Output update(const Vec3& x, float dt_s,
                float external_wave_freq_hz, bool external_valid) {
    return updateImpl_(x, dt_s, external_wave_freq_hz, external_valid);
  }

  Output update(float x, float y, float z, float dt_s,
                float external_wave_freq_hz, bool external_valid) {
    return updateImpl_(Vec3(x, y, z), dt_s, external_wave_freq_hz, external_valid);
  }

  Vec3 currentBaselineSlow() const { return baseline_slow_; }
  Vec3 currentWaveRaw() const { return last_wave_raw_; }
  Vec3 currentWaveClean() const { return last_wave_clean_; }

  float currentWaveFreqHz() const { return f_used_hz_; }
  float currentWavePeriodS() const { return 1.0f / maxf_(f_used_hz_, 1.0e-6f); }

  float currentBaselineCutoffHz() const { return currentBaselineCutoffHz_(); }
  float currentCleanupCutoffHz() const { return currentCleanupCutoffHz_(); }

  bool frequencyValid() const { return isFrequencyValid_(); }
  int8_t selectedAxis() const { return selected_axis_; }

  const Output& lastOutput() const { return last_output_; }

private:
  static constexpr float kPi_ = 3.14159265358979323846f;
  static constexpr uint8_t kMaxCleanupStages_ = 3;

  Config cfg_;

  bool initialized_ = false;
  float time_s_ = 0.0f;

  // Internal learning is disabled until time_s_ reaches this value.
  float learning_enabled_after_t_ = 0.0f;

  // Current input and state for the slow baseline.
  Vec3 current_input_ = Vec3::Zero();
  Vec3 baseline_slow_ = Vec3::Zero();

  // Previous detrended signal used to form a derivative for learning frequency.
  Vec3 wave_learn_prev_ = Vec3::Zero();

  // Filtered slope proxy and previous sample of it.
  Vec3 slope_filt_ = Vec3::Zero();
  Vec3 slope_prev_ = Vec3::Zero();

  // Per-axis slope RMS^2, used to select the learning axis and set threshold.
  Vec3 slope_rms2_ = Vec3::Zero();

  // Shared learned/held dominant frequency.
  float f_used_hz_ = 0.12f;

  // Schmitt trigger state and currently selected learning axis.
  int8_t schmitt_state_ = 0;
  int8_t selected_axis_ = -1;

  // Last positive / negative threshold crossing times on the selected learning channel.
  bool have_last_pos_cross_ = false;
  bool have_last_neg_cross_ = false;
  float last_pos_cross_t_ = 0.0f;
  float last_neg_cross_t_ = 0.0f;

  // Internal period-measurement freshness.
  int valid_period_count_ = 0;
  float last_valid_freq_t_ = -1.0e30f;

  // External frequency freshness is tracked separately.
  float last_external_freq_t_ = -1.0e30f;

  // State for cleanup HP cascade.
  // Each stage stores the LP state used to form HP = input - LP(input).
  Vec3 cleanup_lp_[kMaxCleanupStages_];

  // Cached latest outputs.
  Vec3 last_wave_raw_ = Vec3::Zero();
  Vec3 last_wave_clean_ = Vec3::Zero();

  Output last_output_;

private:
  static bool isFinite_(float x) {
    return isfinite(x) != 0;
  }

  static bool isFiniteVec_(const Vec3& v) {
    return isFinite_(v.x()) && isFinite_(v.y()) && isFinite_(v.z());
  }

  static float clampf_(float x, float lo, float hi) {
    return (x < lo) ? lo : ((x > hi) ? hi : x);
  }

  static float maxf_(float a, float b) {
    return (a > b) ? a : b;
  }

  static float minf_(float a, float b) {
    return (a < b) ? a : b;
  }

  // Clamp each component independently to [-lim, +lim].
  static Vec3 clampAbs_(const Vec3& v, float lim) {
    Vec3 out = v;
    out.x() = clampf_(out.x(), -lim, lim);
    out.y() = clampf_(out.y(), -lim, lim);
    out.z() = clampf_(out.z(), -lim, lim);
    return out;
  }

  // Convert 1st-order LP cutoff frequency to equivalent time constant.
  static float cutoffToTau_(float cutoff_hz) {
    cutoff_hz = maxf_(cutoff_hz, 1.0e-6f);
    return 1.0f / (2.0f * kPi_ * cutoff_hz);
  }

  // Exact exponential coefficient for a first-order filter with time constant tau_s.
  // Returned value is the "old-state weight".
  static float expAlphaFromTau_(float dt_s, float tau_s) {
    tau_s = maxf_(tau_s, 1.0e-6f);
    return expf(-dt_s / tau_s);
  }

  // Linear interpolation of the time at which y crosses 'level' between two samples.
  static float interpCrossTime_(float t_prev, float t_curr,
                                float y_prev, float y_curr,
                                float level) {
    const float dy = y_curr - y_prev;
    if (fabsf(dy) < 1.0e-12f) {
      return t_curr;
    }
    float frac = (level - y_prev) / dy;
    frac = clampf_(frac, 0.0f, 1.0f);
    return t_prev + frac * (t_curr - t_prev);
  }

  void zeroCleanupStates_() {
    for (uint8_t i = 0; i < kMaxCleanupStages_; ++i) {
      cleanup_lp_[i].setZero();
    }
  }

  uint8_t effectiveCleanupStages_() const {
    if (!cfg_.enable_wave_cleanup) return 0;
    return (cfg_.cleanup_stages > kMaxCleanupStages_) ? kMaxCleanupStages_ : cfg_.cleanup_stages;
  }

  // Clamp / repair config values so runtime logic never sees invalid parameters.
  void sanitizeConfig_(Config& c) {
    if (!isFinite_(c.init_wave_freq_hz) || c.init_wave_freq_hz <= 0.0f) c.init_wave_freq_hz = 0.12f;
    if (!isFinite_(c.min_wave_freq_hz)  || c.min_wave_freq_hz  <= 0.0f) c.min_wave_freq_hz  = 0.02f;
    if (!isFinite_(c.max_wave_freq_hz)  || c.max_wave_freq_hz  <= c.min_wave_freq_hz) {
      c.max_wave_freq_hz = maxf_(1.20f, 2.0f * c.min_wave_freq_hz);
    }
    c.init_wave_freq_hz = clampf_(c.init_wave_freq_hz, c.min_wave_freq_hz, c.max_wave_freq_hz);

    if (!isFinite_(c.baseline_cutoff_fraction) || c.baseline_cutoff_fraction <= 0.0f) {
      c.baseline_cutoff_fraction = 0.25f;
    }
    if (!isFinite_(c.min_baseline_cutoff_hz) || c.min_baseline_cutoff_hz <= 0.0f) {
      c.min_baseline_cutoff_hz = 0.003f;
    }
    if (!isFinite_(c.max_baseline_cutoff_hz) || c.max_baseline_cutoff_hz < c.min_baseline_cutoff_hz) {
      c.max_baseline_cutoff_hz = maxf_(0.25f, c.min_baseline_cutoff_hz);
    }

    if (!isFinite_(c.freq_smooth_tau_s) || c.freq_smooth_tau_s <= 0.0f) c.freq_smooth_tau_s = 12.0f;
    if (!isFinite_(c.slope_lpf_tau_s) || c.slope_lpf_tau_s <= 0.0f) c.slope_lpf_tau_s = 0.20f;
    if (!isFinite_(c.slope_rms_tau_s) || c.slope_rms_tau_s <= 0.0f) c.slope_rms_tau_s = 8.0f;

    if (!isFinite_(c.threshold_rms_fraction) || c.threshold_rms_fraction <= 0.0f) {
      c.threshold_rms_fraction = 0.30f;
    }
    if (!isFinite_(c.min_slope_threshold_abs) || c.min_slope_threshold_abs < 0.0f) {
      c.min_slope_threshold_abs = 0.001f;
    }
    if (!isFinite_(c.max_slope_threshold_abs) || c.max_slope_threshold_abs < c.min_slope_threshold_abs) {
      c.max_slope_threshold_abs = c.min_slope_threshold_abs;
    }

    if (!isFinite_(c.startup_hold_s) || c.startup_hold_s < 0.0f) c.startup_hold_s = 2.0f;
    if (!isFinite_(c.freq_timeout_cycles) || c.freq_timeout_cycles <= 0.0f) c.freq_timeout_cycles = 3.0f;

    if (c.freq_learn_axis < -1 || c.freq_learn_axis > 2) c.freq_learn_axis = -1;

    if (!isFinite_(c.cleanup_cutoff_fraction) || c.cleanup_cutoff_fraction <= 0.0f) {
      c.cleanup_cutoff_fraction = 1.0f;
    }
    if (!isFinite_(c.min_cleanup_cutoff_hz) || c.min_cleanup_cutoff_hz <= 0.0f) {
      c.min_cleanup_cutoff_hz = 0.003f;
    }
    if (!isFinite_(c.max_cleanup_cutoff_hz) || c.max_cleanup_cutoff_hz < c.min_cleanup_cutoff_hz) {
      c.max_cleanup_cutoff_hz = maxf_(0.50f, c.min_cleanup_cutoff_hz);
    }
    if (c.cleanup_stages > kMaxCleanupStages_) {
      c.cleanup_stages = kMaxCleanupStages_;
    }

    if (!isFinite_(c.min_dt_s) || c.min_dt_s <= 0.0f) c.min_dt_s = 1.0e-4f;
    if (!isFinite_(c.max_dt_s) || c.max_dt_s < c.min_dt_s) c.max_dt_s = maxf_(2.0f, c.min_dt_s);

    if (!isFinite_(c.output_abs_limit)) c.output_abs_limit = 0.0f;
  }

  // Clamp normal dt values into a safe range.
  // Large dt is handled separately as a discontinuity before this is called.
  float sanitizeDt_(float dt_s) const {
    if (!isFinite_(dt_s)) return cfg_.min_dt_s;
    return clampf_(dt_s, cfg_.min_dt_s, cfg_.max_dt_s);
  }

  // Current adaptive cutoff for the slow baseline LP.
  float currentBaselineCutoffHz_() const {
    float fc = cfg_.baseline_cutoff_fraction * f_used_hz_;
    fc = maxf_(fc, cfg_.min_baseline_cutoff_hz);
    fc = minf_(fc, cfg_.max_baseline_cutoff_hz);
    return fc;
  }

  // Current adaptive cutoff for optional cleanup HP.
  float currentCleanupCutoffHz_() const {
    if (!cfg_.enable_wave_cleanup || effectiveCleanupStages_() == 0) return 0.0f;
    float fc = cfg_.cleanup_cutoff_fraction * currentBaselineCutoffHz_();
    fc = maxf_(fc, cfg_.min_cleanup_cutoff_hz);
    fc = minf_(fc, cfg_.max_cleanup_cutoff_hz);
    return fc;
  }

  // Frequency is valid if either:
  // 1) internal period measurements are fresh enough, or
  // 2) external frequency input was refreshed recently.
  bool isFrequencyValid_() const {
    const float f = maxf_(f_used_hz_, 1.0e-4f);
    const float timeout_s = cfg_.freq_timeout_cycles / f;

    const bool internal_ok =
        (valid_period_count_ >= 2) &&
        ((time_s_ - last_valid_freq_t_) <= timeout_s);

    const bool external_ok =
        ((time_s_ - last_external_freq_t_) <= timeout_s);

    return internal_ok || external_ok;
  }

  // Accept an internally measured period/frequency and smooth it into the held value.
  void acceptFrequencyMeasurement_(float f_meas_hz, float meas_interval_s) {
    if (!isFinite_(f_meas_hz) || f_meas_hz <= 0.0f) return;

    f_meas_hz = clampf_(f_meas_hz, cfg_.min_wave_freq_hz, cfg_.max_wave_freq_hz);

    // Smooth more slowly for large configured tau.
    const float dt_f = maxf_(meas_interval_s, cfg_.min_dt_s);
    const float a_f = expAlphaFromTau_(dt_f, cfg_.freq_smooth_tau_s);
    f_used_hz_ = a_f * f_used_hz_ + (1.0f - a_f) * f_meas_hz;

    last_valid_freq_t_ = time_s_;
    if (valid_period_count_ < 1000000) ++valid_period_count_;
  }

  // Blend in an externally provided frequency without pretending it was an internal period measurement.
  void blendExternalFrequency_(float f_ext_hz, float dt_s) {
    if (!isFinite_(f_ext_hz) || f_ext_hz <= 0.0f) return;

    f_ext_hz = clampf_(f_ext_hz, cfg_.min_wave_freq_hz, cfg_.max_wave_freq_hz);
    const float a_f = expAlphaFromTau_(dt_s, cfg_.freq_smooth_tau_s);
    f_used_hz_ = a_f * f_used_hz_ + (1.0f - a_f) * f_ext_hz;

    // Only freshness timestamp is updated here.
    last_external_freq_t_ = time_s_;
  }

  // Select the scalar channel used for threshold crossing based period measurement.
  // Either fixed by config or chosen as the axis with largest slope RMS.
  int chooseLearningAxis_() const {
    if (cfg_.freq_learn_axis >= 0 && cfg_.freq_learn_axis <= 2) {
      return cfg_.freq_learn_axis;
    }

    int axis = 0;
    float best = slope_rms2_(0);
    if (slope_rms2_(1) > best) { best = slope_rms2_(1); axis = 1; }
    if (slope_rms2_(2) > best) { axis = 2; }
    return axis;
  }

  // Use Schmitt threshold crossings on the selected slope channel to estimate period.
  // Positive crossings are compared only to prior positive crossings, likewise negative.
  void updateFrequencyFromSlope_(float y_prev, float y_curr, float thr, float dt_s) {
    if (time_s_ < learning_enabled_after_t_) {
      return;
    }

    const float t_prev = time_s_ - dt_s;
    const float t_curr = time_s_;

    const float min_period_s = 1.0f / cfg_.max_wave_freq_hz;
    const float max_period_s = 1.0f / cfg_.min_wave_freq_hz;

    // Rising crossing of +threshold.
    if (schmitt_state_ != +1 && y_prev < +thr && y_curr >= +thr) {
      const float t_cross = interpCrossTime_(t_prev, t_curr, y_prev, y_curr, +thr);

      if (have_last_pos_cross_) {
        const float T = t_cross - last_pos_cross_t_;
        if (T >= min_period_s && T <= max_period_s) {
          acceptFrequencyMeasurement_(1.0f / T, T);
        }
      }

      last_pos_cross_t_ = t_cross;
      have_last_pos_cross_ = true;
      schmitt_state_ = +1;
      return;
    }

    // Falling crossing of -threshold.
    if (schmitt_state_ != -1 && y_prev > -thr && y_curr <= -thr) {
      const float t_cross = interpCrossTime_(t_prev, t_curr, y_prev, y_curr, -thr);

      if (have_last_neg_cross_) {
        const float T = t_cross - last_neg_cross_t_;
        if (T >= min_period_s && T <= max_period_s) {
          acceptFrequencyMeasurement_(1.0f / T, T);
        }
      }

      last_neg_cross_t_ = t_cross;
      have_last_neg_cross_ = true;
      schmitt_state_ = -1;
      return;
    }
  }

  // Handle a large sample-time gap as a discontinuity.
  // This avoids generating fake derivatives and bogus frequency estimates.
  void reseedAfterGap_(const Vec3& x, float dt_s_raw) {
    time_s_ += dt_s_raw;

    // Re-apply startup hold after a discontinuity.
    learning_enabled_after_t_ = time_s_ + cfg_.startup_hold_s;

    current_input_ = x;

    // Re-seed baseline to current input so output restarts cleanly.
    baseline_slow_ = x;
    wave_learn_prev_.setZero();

    slope_filt_.setZero();
    slope_prev_.setZero();
    slope_rms2_.setZero();

    schmitt_state_ = 0;
    selected_axis_ =
        (cfg_.freq_learn_axis >= 0 && cfg_.freq_learn_axis <= 2)
            ? cfg_.freq_learn_axis
            : -1;

    have_last_pos_cross_ = false;
    have_last_neg_cross_ = false;
    last_pos_cross_t_ = 0.0f;
    last_neg_cross_t_ = 0.0f;

    zeroCleanupStates_();

    last_wave_raw_.setZero();
    last_wave_clean_.setZero();

    last_output_ = buildOutput_();
  }

  Output updateImpl_(const Vec3& x, float dt_s, float external_wave_freq_hz, bool external_valid) {
    // Ignore non-finite input.
    if (!isFiniteVec_(x)) {
      return last_output_;
    }

    // Ignore non-finite or non-positive dt.
    if (!isFinite_(dt_s) || dt_s <= 0.0f) {
      return last_output_;
    }

    // First valid sample initializes the state around that sample.
    if (!initialized_) {
      reset(x);
      return last_output_;
    }

    // Large dt is treated as a discontinuity, not as a clamped derivative step.
    if (dt_s > cfg_.max_dt_s) {
      reseedAfterGap_(x, dt_s);
      return last_output_;
    }

    dt_s = sanitizeDt_(dt_s);
    time_s_ += dt_s;
    current_input_ = x;

    // ---- Frequency learning path ----
    //
    // Learn frequency from the detrended signal, not from raw input.
    // This reduces contamination from slow drift and offsets.
    const Vec3 wave_learn = x - baseline_slow_;

    // Derivative-like slope proxy used only for frequency learning.
    const Vec3 slope_raw = (wave_learn - wave_learn_prev_) / dt_s;

    // Low-pass the slope to reduce noise before threshold crossing detection.
    const float a_slope = expAlphaFromTau_(dt_s, cfg_.slope_lpf_tau_s);
    slope_filt_ = a_slope * slope_filt_ + (1.0f - a_slope) * slope_raw;

    // Track per-axis RMS^2 of the filtered slope.
    const float a_rms = expAlphaFromTau_(dt_s, cfg_.slope_rms_tau_s);
    slope_rms2_ = a_rms * slope_rms2_ + (1.0f - a_rms) * slope_filt_.cwiseProduct(slope_filt_);

    // Numerical safety.
    slope_rms2_.x() = maxf_(0.0f, slope_rms2_.x());
    slope_rms2_.y() = maxf_(0.0f, slope_rms2_.y());
    slope_rms2_.z() = maxf_(0.0f, slope_rms2_.z());

    // Choose one scalar axis to drive threshold crossing frequency estimation.
    selected_axis_ = (int8_t)chooseLearningAxis_();

    const float slope_rms_sel = sqrtf(maxf_(0.0f, slope_rms2_(selected_axis_)));

    // Threshold is adaptive with RMS, but clamped.
    float slope_thr = cfg_.threshold_rms_fraction * slope_rms_sel;
    slope_thr = maxf_(slope_thr, cfg_.min_slope_threshold_abs);
    slope_thr = minf_(slope_thr, cfg_.max_slope_threshold_abs);

    // Update internal frequency estimate from threshold crossings.
    updateFrequencyFromSlope_(slope_prev_(selected_axis_),
                              slope_filt_(selected_axis_),
                              slope_thr, dt_s);

    // External frequency can additionally guide the held shared frequency.
    if (external_valid) {
      blendExternalFrequency_(external_wave_freq_hz, dt_s);
    }

    // ---- Baseline detrending path ----
    //
    // Scalar adaptive cutoff shared across axes, but vector LP state per axis.
    const float fc_base = currentBaselineCutoffHz_();
    const float a_base = expf(-2.0f * kPi_ * fc_base * dt_s);

    // Slow baseline LP.
    baseline_slow_ = a_base * baseline_slow_ + (1.0f - a_base) * x;

    // Primary detrended output.
    Vec3 wave_raw = x - baseline_slow_;

    // ---- Optional cleanup path ----
    //
    // Additional HP on wave_raw only:
    //   stage_out = stage_in - LP(stage_in)
    //
    // This leaves baseline_slow_ unchanged.
    Vec3 wave_clean = wave_raw;
    const uint8_t n_cleanup = effectiveCleanupStages_();

    if (n_cleanup > 0) {
      const float fc_clean = currentCleanupCutoffHz_();
      const float a_clean = expf(-2.0f * kPi_ * fc_clean * dt_s);

      for (uint8_t i = 0; i < n_cleanup; ++i) {
        cleanup_lp_[i] = a_clean * cleanup_lp_[i] + (1.0f - a_clean) * wave_clean;
        wave_clean = wave_clean - cleanup_lp_[i];
      }

      // Keep unused stages zeroed.
      for (uint8_t i = n_cleanup; i < kMaxCleanupStages_; ++i) {
        cleanup_lp_[i].setZero();
      }
    } else {
      zeroCleanupStates_();
      wave_clean = wave_raw;
    }

    // Optional output limiter.
    if (cfg_.output_abs_limit > 0.0f) {
      wave_raw = clampAbs_(wave_raw, cfg_.output_abs_limit);
      wave_clean = clampAbs_(wave_clean, cfg_.output_abs_limit);
    }

    // Save state for next derivative/crossing update.
    wave_learn_prev_ = wave_learn;
    slope_prev_ = slope_filt_;

    // Cache outputs.
    last_wave_raw_ = wave_raw;
    last_wave_clean_ = wave_clean;
    last_output_ = buildOutput_();
    return last_output_;
  }

  Output buildOutput_() const {
    Output out;
    out.input = current_input_;
    out.baseline_slow = baseline_slow_;
    out.wave_raw = last_wave_raw_;
    out.wave_clean = last_wave_clean_;

    out.wave_freq_hz = f_used_hz_;
    out.wave_period_s = 1.0f / maxf_(f_used_hz_, 1.0e-6f);

    out.baseline_cutoff_hz = currentBaselineCutoffHz_();
    out.baseline_tau_s = cutoffToTau_(out.baseline_cutoff_hz);

    out.cleanup_cutoff_hz = currentCleanupCutoffHz_();

    // Convert stored RMS^2 to RMS for output.
    out.slope_rms_xyz.x() = sqrtf(maxf_(0.0f, slope_rms2_.x()));
    out.slope_rms_xyz.y() = sqrtf(maxf_(0.0f, slope_rms2_.y()));
    out.slope_rms_xyz.z() = sqrtf(maxf_(0.0f, slope_rms2_.z()));

    const int axis = (selected_axis_ >= 0 && selected_axis_ <= 2) ? selected_axis_ : 0;
    out.slope_rms_selected = out.slope_rms_xyz(axis);

    float thr = cfg_.threshold_rms_fraction * out.slope_rms_selected;
    thr = maxf_(thr, cfg_.min_slope_threshold_abs);
    thr = minf_(thr, cfg_.max_slope_threshold_abs);
    out.slope_threshold = thr;

    out.freq_valid = isFrequencyValid_();
    out.schmitt_state = schmitt_state_;
    out.selected_axis = selected_axis_;
    return out;
  }
};
