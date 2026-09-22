/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R TimeVaryingGainNLO marine INS (+ optional IMU Calibration Wizard)

  Runs the time-varying-gain nonlinear observer of Bryne/Fossen/Johansen
  (src/nlo/TimeVaryingGainNLO.h) through its Arduino adapter
  (src/nlo/TimeVarGainNLO_Adapter.h).

  Produces:
    - magnetic (or true) compass heading from observer yaw
    - AHRS roll/pitch/yaw
    - heave displacement and heave speed
    - dominant wave frequency and the scheduled translational gain theta

  DEVICE CONFIGURATION

    The AtomS3R carries no position reference, so the observer runs as

      TimeVarGainNloAdapter<WithGNSS = false,
                            Mag      = NloMagType::Magnetometer>

    With WithGNSS = false all three translational axes are aided by the
    paper's virtual zero-mean integrated-position measurement (Sec. II-A.2).
    That keeps xi observable on every axis, so the estimated specific force
    used as the attitude reference is a genuine estimate.

    Use heave (Z) as the navigation product. Horizontal displacement is
    carried only to keep xi observable; without a real position aid it is not
    a calibrated surge/sway output. See the header comment of
    src/nlo/TimeVaryingGainNLO.h.

  FRAMES

    readImuMapped() + RuntimeCals deliver the project body-NED convention:
      x forward, y starboard, z down

    which is exactly what the observer expects:
      gyro           BODY rad/s
      specific force BODY m/s^2, level and still is about [0, 0, -g]

    so no axis remapping is done in this sketch, unlike the Z-up Mahony
    sketch next door.

    Observer state frame is NED, +Z down. Euler output is (roll, pitch, yaw)
    with yaw measured clockwise from the reference north direction.

  MAGNETIC HEADING

    The observer's magnetic declination is left at zero, so its yaw is the
    MAGNETIC heading. True heading, if wanted, is produced at the output
    stage from SEA_STATE_MAG_DECLINATION_DEG, exactly as the PII sketch does.

  GYRO BIAS

    Not estimated in this sketch. The observer carries its own gyro-bias
    state driven by the integral gain kI, and snapshot().tvg.gyro_bias_b is
    that estimate; it is what the rate-of-turn output below is corrected
    with. A separate stillness-gated bias average, as used by the Kalman and
    PII sketches, would fight it.
*/

#include <Arduino.h>

#ifndef SEA_STATE_ENABLE_WIZARD
  #define SEA_STATE_ENABLE_WIZARD 1
#endif

#define ARDUINO_PLOTTER 1

#include <M5Unified.h>
#include <algorithm>
#include <cmath>
#include <cstring>

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
  #define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif
#include <ArduinoEigenDense.h>

#include <ArduinoOceanImu.h>

#include "AtomS3R/AtomS3R_CompassUI.h"
#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"
#if SEA_STATE_ENABLE_WIZARD
  #include "AtomS3R/ImuCalWizardRunner.h"
#endif
#include "detrend/AdaptiveWaveDetrender.h"
#include "nlo/TimeVarGainNLO_Adapter.h"
#include "nmea/NmeaCompass.h"

#ifndef SEA_STATE_UI_DEFAULT_GRAPHICS
  #define SEA_STATE_UI_DEFAULT_GRAPHICS 1
#endif

#ifndef SEA_STATE_SERIAL_NMEA
  #define SEA_STATE_SERIAL_NMEA 1
#endif

#ifndef SEA_STATE_NMEA_TALKER
  #define SEA_STATE_NMEA_TALKER "II"
#endif

#ifndef SEA_STATE_OUTPUT_TRUE_HEADING
  #define SEA_STATE_OUTPUT_TRUE_HEADING 0
#endif

#ifndef SEA_STATE_MAG_DECLINATION_DEG
  #define SEA_STATE_MAG_DECLINATION_DEG 0.0f
#endif

#ifndef SEA_STATE_MAG_HEADING_USER_OFFSET_DEG
  #define SEA_STATE_MAG_HEADING_USER_OFFSET_DEG 0.0f
#endif

#ifndef SEA_STATE_USE_STRICT_MAG_FIELD_GATE
  #define SEA_STATE_USE_STRICT_MAG_FIELD_GATE 0
#endif

/*
  Run the paper's fixed theta = 1 instead of the adapter's wave-frequency
  schedule. Off by default: at theta = 1 the aiding loop's high-pass corner
  sits at 0.41 rad/s, inside the ocean wave band, which attenuates and
  phase-shifts the wave the device is meant to measure.
*/
#ifndef SEA_STATE_NLO_FIXED_THETA
  #define SEA_STATE_NLO_FIXED_THETA 0
#endif

static constexpr float APP_G_STD = atoms3r_ical::ImuCalCfg::g_std;
static constexpr float APP_FREQ_GUESS = 0.30f;

static constexpr float LOOP_HZ = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = static_cast<uint32_t>(1000000.0f / LOOP_HZ);

static constexpr uint32_t UI_REFRESH_MS = 100;
static constexpr uint32_t DEBUG_SERIAL_MS = 100;
static constexpr uint32_t NMEA_SERIAL_MS = 80;
static constexpr uint32_t NMEA_STATUS_MS = 5000;
static constexpr uint32_t NMEA_TEMP_MS = 500;

static constexpr float MAG_PRESENT_MIN_UT = 5.0f;
static constexpr float MAG_PRESENT_MAX_UT = 200.0f;

static constexpr float MAG_FIELD_MIN_UT = 20.0f;
static constexpr float MAG_FIELD_MAX_UT = 80.0f;

static constexpr uint32_t MAG_UPDATE_SPACING_MS = 35u;
static constexpr uint32_t HEADING_MAG_TIMEOUT_MS = 2000u;

static constexpr float ROT_TAU_S = 0.1f;

using namespace atoms3r_ical;
using Vector3f = Eigen::Vector3f;
using Quaternionf = Eigen::Quaternionf;

static inline float clampf_(float x, float lo, float hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

static inline float wrap360_(float deg) {
  while (deg < 0.0f) deg += 360.0f;
  while (deg >= 360.0f) deg -= 360.0f;
  return deg;
}

static inline float wrap180_(float deg) {
  while (deg <= -180.0f) deg += 360.0f;
  while (deg >   180.0f) deg -= 360.0f;
  return deg;
}

static inline float outputHeadingFromMagnetic_(float magnetic_deg) {
#if SEA_STATE_OUTPUT_TRUE_HEADING
  return wrap360_(magnetic_deg + SEA_STATE_MAG_DECLINATION_DEG);
#else
  return magnetic_deg;
#endif
}

static inline Vector3f quatRotate_(const Quaternionf& q, const Vector3f& v) {
  const Vector3f qv(q.x(), q.y(), q.z());
  const Vector3f t = 2.0f * qv.cross(v);
  return v + q.w() * t + qv.cross(t);
}

class FusionApp {
 public:
  void begin() {
    delay(50);
    Serial.begin(115200);
    delay(100);

    auto cfg = M5.config();
    cfg.internal_imu = true;
    M5.begin(cfg);
    clearM5UnifiedImuCalibration();

    ui_.begin();

    if (use_graphics_) {
      ui_.setReadRotation();
      compass_ui_.begin();
      compass_ui_ready_ = compass_ui_.ok();
      if (!compass_ui_ready_) use_graphics_ = false;
    }

    reloadBlobAndRuntime_();

#if SEA_STATE_ENABLE_WIZARD
    if (!have_blob_) {
      Serial.println("[BOOT] No saved calibration. Starting wizard...");
      const bool saved = runWizardFlow_(true);
      if (saved) {
        Serial.println("[BOOT] Wizard saved calibration.");
      } else {
        Serial.println("[BOOT] Wizard did not save calibration. Running with raw values.");
      }
    }
#endif

    reinitImu_();
    resetFusion_();
    drawHomeStatic_();

    delay(100);
  }

  void tick() {
    const uint32_t loop_start_us = micros();

#if SEA_STATE_ENABLE_WIZARD
    Input::update();

    if (Input::tapPressed()) {
      tap_count_++;
      tap_deadline_ms_ = millis() + M5UiCfg::MENU_TAP_WINDOW_MS;
      drawHomePending_();
      Serial.printf("[TAP] count=%d\n", tap_count_);
    }

    if (tap_count_ > 0 &&
        static_cast<int32_t>(millis() - tap_deadline_ms_) > 0) {
      if (tap_count_ >= 3) handleErase_();
      else handleRunWizard_();

      tap_count_ = 0;
      tap_deadline_ms_ = 0;
      drawHomeStatic_();
    }
#endif

    ImuSample sample{};
    const uint32_t sample_us = micros();
    const uint32_t update_mask = M5.Imu.update();
    const bool got_sample = readImuMapped(M5.Imu, update_mask, sample_us, sample);

    if (got_sample) {
      stale_frame_count_ = 0;
      updateFilter_(sample);
    } else if (stale_frame_count_ < 0xFFFFu) {
      stale_frame_count_++;
    }

    updateUI_();
    streamSerial_();
    waitForNextLoopTick_(loop_start_us);
  }

 private:
  using UI = atoms3r_ical::M5Ui;

  /*
    No position reference on this device, 3D magnetometer present.
    PLL is the frequency tracker the committed NLO simulation study runs.
  */
  using Fusion =
      TimeVarGainNloAdapter<false, NloMagType::Magnetometer, float, TrackerType::PLL>;

  bool use_graphics_ = (SEA_STATE_UI_DEFAULT_GRAPHICS != 0);

  CompassUI compass_ui_{};
  bool compass_ui_ready_ = false;
  UI ui_{};

  ImuCalStoreNvs store_{};
  bool have_blob_ = false;
  ImuCalBlobV2 blob_{};
  RuntimeCals runtime_{};

#if SEA_STATE_ENABLE_WIZARD
  int tap_count_ = 0;
  uint32_t tap_deadline_ms_ = 0;
#endif

  uint32_t last_ui_ms_ = 0;
  uint32_t last_serial_ms_ = 0;
  uint32_t last_status_nmea_ms_ = 0;
  uint32_t last_temp_nmea_ms_ = 0;

  Fusion fusion_{};

  uint32_t mag_gate_last_ms_ = 0;
  uint32_t last_mag_correction_ms_ = 0;

  Vector3f a_cal_ = Vector3f::Zero();
  Vector3f w_cal_ = Vector3f::Zero();
  Vector3f m_cal_ = Vector3f::Constant(NAN);

  float a_raw_norm_ = 0.0f;

  float dt_ = 0.0f;

  float roll_deg_ = 0.0f;
  float pitch_deg_ = 0.0f;
  float yaw_deg_ = 0.0f;

  float heading_mag_deg_ = 0.0f;
  float heading_deg_ = 0.0f;
  bool heading_valid_ = false;

  bool nlo_initialized_ = false;
  bool nlo_bootstrapping_ = false;

  float heave_m_ = 0.0f;
  float heave_speed_mps_ = 0.0f;
  float heave_accel_mps2_ = 0.0f;
  float wave_hz_ = APP_FREQ_GUESS;
  float wave_conf_ = 0.0f;
  bool wave_locked_ = false;
  float nlo_theta_ = 1.0f;
  float gyro_bias_norm_ = 0.0f;

  float heave_raw_m_ = 0.0f;
  float heave_baseline_m_ = 0.0f;
  float heave_wave_raw_m_ = 0.0f;
  float heave_wave_clean_m_ = 0.0f;

  bool mag_present_ = false;
  bool mag_field_sane_ = false;
  bool mag_fresh_ = false;
  bool mag_used_ = false;
  float mag_norm_uT_ = NAN;
  float imu_temp_c_ = NAN;

  uint16_t stale_frame_count_ = 0;
  bool have_last_sample_us_ = false;
  uint32_t last_sample_us_ = 0;

  bool rot_inited_ = false;
  float rot_dpm_filt_ = 0.0f;

  uint32_t rate_window_ms_ = 0;
  uint32_t imu_sample_count_ = 0;
  uint32_t mag_sample_count_ = 0;
  float imu_sample_hz_ = 0.0f;
  float mag_sample_hz_ = 0.0f;

  AdaptiveWaveDetrender z_detrender_{};

 private:
  static void waitForNextLoopTick_(uint32_t loop_start_us) {
    const uint32_t elapsed_us = micros() - loop_start_us;
    if (elapsed_us < LOOP_PERIOD_US) {
      delayMicroseconds(LOOP_PERIOD_US - elapsed_us);
    }
  }

  float nominalFusionDt_() const {
    return 1.0f / LOOP_HZ;
  }

  bool updateMagFreshGate_(bool mag_candidate_ok, uint32_t now_ms) {
    if (!mag_candidate_ok) {
      mag_gate_last_ms_ = 0;
      return false;
    }

    if (mag_gate_last_ms_ == 0) {
      mag_gate_last_ms_ = now_ms;
      return true;
    }

    if ((now_ms - mag_gate_last_ms_) < MAG_UPDATE_SPACING_MS) {
      return false;
    }

    mag_gate_last_ms_ = now_ms;
    return true;
  }

  float computeFusionDtFromSampleTimestamp_(const ImuSample& s) {
    const float dt_nom = nominalFusionDt_();

    if (!have_last_sample_us_) {
      have_last_sample_us_ = true;
      last_sample_us_ = s.sample_us;
      return dt_nom;
    }

    const uint32_t dt_us = s.sample_us - last_sample_us_;
    last_sample_us_ = s.sample_us;

    const float dt_s = static_cast<float>(dt_us) * 1.0e-6f;
    if (!(dt_s > 0.0f) || !std::isfinite(dt_s)) return dt_nom;
    if (dt_s > 0.05f) return dt_nom;

    return dt_s;
  }

  void reloadBlobAndRuntime_() {
    have_blob_ = store_.load(blob_);
    if (!have_blob_) {
      std::memset(&blob_, 0, sizeof(blob_));
    }
    runtime_.rebuildFromBlob(blob_);
  }

  bool runWizardFlow_(bool boot_mode) {
#if SEA_STATE_ENABLE_WIZARD
    (void)boot_mode;
    clearM5UnifiedImuCalibration();

    ImuCalBlobV2 saved{};
    const bool did_save = runImuCalWizard(ui_, store_, saved);

    if (did_save) {
      blob_ = saved;
      have_blob_ = true;
      runtime_.rebuildFromBlob(blob_);
      return true;
    }
#else
    (void)boot_mode;
#endif
    return false;
  }

  void reinitImu_() {
    if (!M5.Imu.isEnabled()) {
      Serial.println("[BOOT] M5.Imu is not enabled");
      ui_.fail("IMU", "M5.Imu disabled");
      while (true) delay(100);
    }

    mag_gate_last_ms_ = 0;
    last_mag_correction_ms_ = 0;

    rate_window_ms_ = millis();
    imu_sample_count_ = 0;
    mag_sample_count_ = 0;
    imu_sample_hz_ = 0.0f;
    mag_sample_hz_ = 0.0f;

    have_last_sample_us_ = false;
    last_sample_us_ = 0;
  }

  void resetFusion_() {
    auto& cfg = fusion_.config();

    /*
      Every observer gain stays at the adapter default, which is the set
      published in Bryne/Fossen/Johansen Sec. IV-C. The scalar tuning
      parameter theta is the documented place to move the translational
      aiding loop in frequency, and the adapter schedules it on the tracked
      wave frequency; see TimeVarGainNloAdapter::Config. Only the
      app/hardware-specific fields are touched here.
    */
    cfg.gravity_mps2 = APP_G_STD;
    cfg.filter.gravity_mps2 = APP_G_STD;

#if SEA_STATE_NLO_FIXED_THETA
    cfg.auto_theta_from_wave_freq = false;
    cfg.filter.theta = 1.0f;
#else
    cfg.auto_theta_from_wave_freq = true;
#endif

    /*
      Magnetic frame. The declination correction is applied to the reported
      heading instead, so observer yaw stays magnetic and matches the HDM
      sentence below.
    */
    cfg.filter.magnetic_declination_rad = 0.0f;

    /*
      The sketch already gates the magnetometer on field norm before handing
      it over (mag_present_ / mag_field_sane_), so the observer's own norm
      gate is left disabled rather than applied twice with a second set of
      thresholds.
    */
    cfg.filter.expected_mag_norm = 0.0f;

    /*
      Real device startup: hold the outputs at zero until the Mahony
      bootstrap has seeded attitude and the observer is initialized.
    */
    cfg.run_filter_before_initialized = false;
    cfg.mahony_bootstrap_enabled = true;

    fusion_.reset();

    // IMU installation lever arm (off-CoG) correction: not available here.
    //
    // An IMU mounted away from the vessel centre of gravity measures
    // a_imu = a_cog + (alpha x r) + omega x (omega x r) for a body-fixed
    // offset r, which leaks roll/pitch motion into the heave channel.
    // TimeVaryingGainNLO takes no lever-arm input, so this sketch has nothing
    // to set and implicitly assumes r = 0 (IMU at the CoG) -- mount it as
    // close to the CoG as the boat allows. The OU-II / OU-III sketches in
    // this same directory do expose the knob and show it set explicitly to
    // zero with instructions for entering a real offset.

    rot_inited_ = false;
    rot_dpm_filt_ = 0.0f;

    roll_deg_ = 0.0f;
    pitch_deg_ = 0.0f;
    yaw_deg_ = 0.0f;

    heading_mag_deg_ = 0.0f;
    heading_deg_ = 0.0f;
    heading_valid_ = false;

    nlo_initialized_ = false;
    nlo_bootstrapping_ = false;

    mag_gate_last_ms_ = 0;
    last_mag_correction_ms_ = 0;

    mag_present_ = false;
    mag_field_sane_ = false;
    mag_fresh_ = false;
    mag_used_ = false;
    mag_norm_uT_ = NAN;
    imu_temp_c_ = NAN;

    stale_frame_count_ = 0;
    have_last_sample_us_ = false;
    last_sample_us_ = 0;

    heave_m_ = 0.0f;
    heave_speed_mps_ = 0.0f;
    heave_accel_mps2_ = 0.0f;
    wave_hz_ = APP_FREQ_GUESS;
    wave_conf_ = 0.0f;
    wave_locked_ = false;
    nlo_theta_ = cfg.filter.theta;
    gyro_bias_norm_ = 0.0f;

    heave_raw_m_ = 0.0f;
    heave_baseline_m_ = 0.0f;
    heave_wave_raw_m_ = 0.0f;
    heave_wave_clean_m_ = 0.0f;

    AdaptiveWaveDetrender::Config dcfg{};
    dcfg.init_wave_freq_hz = APP_FREQ_GUESS;
    dcfg.min_wave_freq_hz  = 0.02f;
    dcfg.max_wave_freq_hz  = 1.20f;

    dcfg.baseline_cutoff_fraction = 0.25f;
    dcfg.min_baseline_cutoff_hz   = 0.003f;
    dcfg.max_baseline_cutoff_hz   = 0.25f;

    dcfg.freq_smooth_tau_s = 12.0f;
    dcfg.slope_lpf_tau_s   = 0.20f;
    dcfg.slope_rms_tau_s   = 8.0f;

    dcfg.threshold_rms_fraction  = 0.15f;
    dcfg.min_slope_threshold_abs = 0.002f;
    dcfg.max_slope_threshold_abs = 1.0e9f;

    dcfg.startup_hold_s      = 2.0f;
    dcfg.freq_timeout_cycles = 3.0f;

    dcfg.enable_wave_cleanup     = true;
    dcfg.cleanup_cutoff_fraction = 1.0f;
    dcfg.min_cleanup_cutoff_hz   = 0.003f;
    dcfg.max_cleanup_cutoff_hz   = 0.50f;
    dcfg.cleanup_stages          = 2;

    dcfg.min_dt_s = 1.0e-4f;
    dcfg.max_dt_s = 0.25f;
    dcfg.output_abs_limit = 0.0f;

    z_detrender_.setConfig(dcfg);
    z_detrender_.reset(0.0f);
  }

#if SEA_STATE_ENABLE_WIZARD
  void handleErase_() {
    Serial.println("[HOME] ERASE");
    if (!ui_.eraseConfirm()) {
      Serial.println("[HOME] erase cancelled");
      return;
    }

    store_.erase();
    reloadBlobAndRuntime_();
    reinitImu_();
    resetFusion_();
  }

  void handleRunWizard_() {
    Serial.println("[HOME] RUN WIZARD");

    const bool saved = runWizardFlow_(false);
    if (!saved) ui_.notSavedNotice();

    reinitImu_();
    resetFusion_();
  }

  void drawHomePending_() {
    ui_.setReadRotation();
    ui_.title("COMPASS");
    M5.Display.printf("Tap count: %d\n", tap_count_);
    ui_.line("");
    ui_.line("Wait...");
    ui_.line("1 tap=CAL");
    ui_.line("3 taps=ERASE");

    int32_t remain = static_cast<int32_t>(tap_deadline_ms_ - millis());
    remain = remain < 0 ? 0 : remain;

    const float t01 =
        1.0f -
        static_cast<float>(remain) /
        static_cast<float>(M5UiCfg::MENU_TAP_WINDOW_MS);

    ui_.bar01(t01);
  }
#endif

  void updateFilter_(const ImuSample& s) {
    const uint32_t now_ms = millis();

    dt_ = computeFusionDtFromSampleTimestamp_(s);
    const float tempC = std::isfinite(s.tempC) ? s.tempC : 35.0f;
    imu_temp_c_ = tempC;
    ++imu_sample_count_;

    a_raw_norm_ = s.a.norm();

    a_cal_ = runtime_.applyAccel(s.a, tempC);
    w_cal_ = runtime_.applyGyro(s.w, tempC);
    m_cal_ = runtime_.applyMag(s.m);

    mag_norm_uT_ = m_cal_.norm();

    mag_present_ =
        std::isfinite(mag_norm_uT_) &&
        mag_norm_uT_ >= MAG_PRESENT_MIN_UT &&
        mag_norm_uT_ <= MAG_PRESENT_MAX_UT;

    mag_field_sane_ =
        std::isfinite(mag_norm_uT_) &&
        mag_norm_uT_ >= MAG_FIELD_MIN_UT &&
        mag_norm_uT_ <= MAG_FIELD_MAX_UT;

    mag_fresh_ = updateMagFreshGate_(mag_present_, now_ms);
    if (mag_fresh_) ++mag_sample_count_;

#if SEA_STATE_USE_STRICT_MAG_FIELD_GATE
    const bool mag_usable = mag_present_ && mag_field_sane_;
#else
    const bool mag_usable = mag_present_;
#endif

    mag_used_ = mag_usable;

    /*
      The observer consumes the magnetometer as a continuous DIRECTION
      reference: the horizontal projection of m_b is compared against the
      magnetic reference vector in NED on every update, which is a
      zero-order hold on the latest reading, not an incremental measurement
      that could be double-counted. So the current calibrated sample is
      offered on every IMU tick and only the gate decides, rather than the
      ~25 Hz freshness gate the discrete Kalman updates next door need.
    */
    if (mag_usable) {
      fusion_.setMagBody(m_cal_, true);
      last_mag_correction_ms_ = now_ms;
    } else {
      fusion_.clearMag();
    }

    fusion_.update(dt_, w_cal_, a_cal_);

    const auto snap = fusion_.snapshot();

    nlo_initialized_ = snap.initialized;
    nlo_bootstrapping_ = snap.mahony_bootstrap_active;

    roll_deg_ = wrap180_(snap.euler_rad.x() * RAD_TO_DEG);
    pitch_deg_ = wrap180_(snap.euler_rad.y() * RAD_TO_DEG);
    yaw_deg_ = wrap360_(snap.euler_rad.z() * RAD_TO_DEG);

    heading_mag_deg_ =
        wrap360_(yaw_deg_ + SEA_STATE_MAG_HEADING_USER_OFFSET_DEG);
    heading_deg_ = outputHeadingFromMagnetic_(heading_mag_deg_);

    heading_valid_ =
        nlo_initialized_ &&
        last_mag_correction_ms_ != 0 &&
        (now_ms - last_mag_correction_ms_) <= HEADING_MAG_TIMEOUT_MS;

    /*
      Rate of turn from the observer's own bias-corrected body rate, rotated
      into NED. snapshot().tvg.gyro_bias_b is the observer's gyro-bias state.
    */
    gyro_bias_norm_ = snap.tvg.gyro_bias_norm;

    const Vector3f w_use = w_cal_ - snap.tvg.gyro_bias_b;
    const Vector3f w_ned = quatRotate_(snap.q_nb, w_use);

    float rot_dpm_meas = w_ned.z() * RAD_TO_DEG * 60.0f;
    rot_dpm_meas = clampf_(rot_dpm_meas, -720.0f, 720.0f);

    const float alpha_r = 1.0f - expf(-dt_ / ROT_TAU_S);
    if (!rot_inited_) {
      rot_inited_ = true;
      rot_dpm_filt_ = rot_dpm_meas;
    } else {
      rot_dpm_filt_ += alpha_r * (rot_dpm_meas - rot_dpm_filt_);
    }

    /*
      disp_zu / vel_zu / acc_zu are the Z-up convenience outputs. Position
      and velocity already have the unobservable DC removed by the adapter's
      report high-pass; see TimeVarGainNloAdapter::Config.
    */
    heave_m_ = snap.disp_zu.z();
    heave_speed_mps_ = snap.vel_zu.z();
    heave_accel_mps2_ = snap.acc_zu.z();
    heave_raw_m_ = heave_m_;

    nlo_theta_ = snap.tvg.theta;
    wave_conf_ = snap.tvg.wave_freq_confidence;
    wave_locked_ = snap.tvg.wave_freq_locked;

    const float nlo_freq_hz = snap.tvg.wave_freq_hz;
    if (std::isfinite(nlo_freq_hz) && nlo_freq_hz > 1e-6f) {
      wave_hz_ = nlo_freq_hz;
    }

    const auto z_det =
        z_detrender_.update(heave_raw_m_, dt_, wave_hz_, wave_hz_ > 1e-6f);

    heave_baseline_m_ = z_det.baseline_slow;
    heave_wave_raw_m_ = z_det.wave_raw;
    heave_wave_clean_m_ = z_det.wave_clean;

    if (rate_window_ms_ == 0) rate_window_ms_ = now_ms;
    const uint32_t rate_elapsed_ms = now_ms - rate_window_ms_;
    if (rate_elapsed_ms >= 1000u) {
      const float scale = 1000.0f / static_cast<float>(rate_elapsed_ms);
      imu_sample_hz_ = static_cast<float>(imu_sample_count_) * scale;
      mag_sample_hz_ = static_cast<float>(mag_sample_count_) * scale;
      imu_sample_count_ = 0;
      mag_sample_count_ = 0;
      rate_window_ms_ = now_ms;
    }
  }

  void drawHomeStatic_() {
    ui_.setReadRotation();
    ui_.title("COMPASS");

    M5.Display.printf("BLOB: %s\n", have_blob_ ? "YES" : "NO");
    M5.Display.printf("A:%d G:%d M:%d\n",
                      static_cast<int>(runtime_.acc.ok),
                      static_cast<int>(runtime_.gyr.ok),
                      static_cast<int>(runtime_.mag.ok));

#if SEA_STATE_ENABLE_WIZARD
    ui_.line("Tap: calibrate");
    ui_.line("Tap x3: erase");
#else
    ui_.line("Wizard: DISABLED");
#endif

#if SEA_STATE_OUTPUT_TRUE_HEADING
    M5.Display.printf("HDG TRUE decl=%+.1f\n",
                      static_cast<double>(SEA_STATE_MAG_DECLINATION_DEG));
#else
    ui_.line("HDG: MAG");
#endif

    M5.Display.printf("UserOff:%+.1f\n",
                      static_cast<double>(SEA_STATE_MAG_HEADING_USER_OFFSET_DEG));

#if SEA_STATE_USE_STRICT_MAG_FIELD_GATE
    ui_.line("Mag gate: STRICT");
#else
    ui_.line("Mag gate: LOOSE");
#endif

#if SEA_STATE_NLO_FIXED_THETA
    ui_.line("Fusion: NLO t=1");
#else
    ui_.line("Fusion: NLO");
#endif
  }

  void updateUI_() {
#if SEA_STATE_ENABLE_WIZARD
    if (tap_count_ > 0) return;
#endif

    const uint32_t now_ms = millis();
    if (now_ms - last_ui_ms_ < UI_REFRESH_MS) return;
    last_ui_ms_ = now_ms;

    if (use_graphics_ && compass_ui_ready_) updateUI_graphics_();
    else                                    updateUI_text_();
  }

  void updateUI_graphics_() {
    ui_.setReadRotation();

    const bool tiltWarn =
        (fabsf(roll_deg_) > 35.0f) ||
        (fabsf(pitch_deg_) > 35.0f);

    const float hdg_draw = heading_valid_ ? heading_deg_ : 0.0f;

    compass_ui_.draw(hdg_draw,
                     heading_valid_ && mag_used_,
                     mag_norm_uT_,
                     tiltWarn);
  }

  void updateUI_text_() {
    ui_.setReadRotation();
    ui_.title("COMPASS");

#if SEA_STATE_OUTPUT_TRUE_HEADING
    M5.Display.printf("HDG:%7.1f T %s\n",
                      static_cast<double>(heading_deg_),
                      heading_valid_ ? "deg" : "WAIT");
    M5.Display.printf("MAG:%7.1f M\n",
                      static_cast<double>(heading_mag_deg_));
#else
    M5.Display.printf("HDG:%7.1f M %s\n",
                      static_cast<double>(heading_mag_deg_),
                      heading_valid_ ? "deg" : "WAIT");
#endif

    M5.Display.printf("ROL:%7.1f deg\n", static_cast<double>(roll_deg_));
    M5.Display.printf("PIT:%7.1f deg\n", static_cast<double>(pitch_deg_));

    M5.Display.printf("NLO:%s\n",
                      nlo_initialized_ ? "RUN"
                                       : (nlo_bootstrapping_ ? "BOOT" : "WAIT"));

    M5.Display.printf("HEV:%7.3f m\n", static_cast<double>(heave_wave_clean_m_));
    M5.Display.printf("FRQ:%7.3f Hz%s\n",
                      static_cast<double>(wave_hz_),
                      wave_locked_ ? "*" : " ");
    M5.Display.printf("THT:%7.3f\n", static_cast<double>(nlo_theta_));

    M5.Display.printf("MAG:%s %s %s\n",
                      mag_present_ ? "OK " : "BAD",
                      mag_field_sane_ ? "FIELD" : "DIST",
                      mag_used_ ? "USE" : "SKIP");

    M5.Display.printf("|m|:%7.1f uT\n", static_cast<double>(mag_norm_uT_));
    M5.Display.printf("|aR|:%5.2f |aC|:%5.2f\n",
                      static_cast<double>(a_raw_norm_),
                      static_cast<double>(a_cal_.norm()));
  }

  void streamSerial_() {
    const uint32_t now_ms = millis();

#if SEA_STATE_SERIAL_NMEA
    if (now_ms - last_serial_ms_ < NMEA_SERIAL_MS) return;
#else
    if (now_ms - last_serial_ms_ < DEBUG_SERIAL_MS) return;
#endif

    last_serial_ms_ = now_ms;

#if SEA_STATE_SERIAL_NMEA

    if (heading_valid_) {
      nmea_hdm(SEA_STATE_NMEA_TALKER, heading_mag_deg_);
    }

    nmea_xdr_pitch_roll(SEA_STATE_NMEA_TALKER, pitch_deg_, roll_deg_);
    nmea_xdr_heave(SEA_STATE_NMEA_TALKER, heave_wave_clean_m_);
    nmea_xdr_heave_speed(SEA_STATE_NMEA_TALKER, heave_speed_mps_);

    if (now_ms - last_temp_nmea_ms_ >= NMEA_TEMP_MS) {
      last_temp_nmea_ms_ = now_ms;
      nmea_xdr_imu_temp(SEA_STATE_NMEA_TALKER, imu_temp_c_);
    }

    if (now_ms - last_status_nmea_ms_ >= NMEA_STATUS_MS) {
      last_status_nmea_ms_ = now_ms;
      nmea_txt_ins_status(SEA_STATE_NMEA_TALKER,
                          nlo_initialized_,
                          nlo_initialized_ && std::isfinite(heave_wave_clean_m_),
                          heading_valid_ && mag_used_,
                          // The observer's integral gyro-bias loop runs for as
                          // long as it is initialized; there is no separate
                          // stillness-gated learning phase to report.
                          nlo_initialized_,
                          // TimeVaryingGainNLO carries no accelerometer-bias
                          // state.
                          false,
                          imu_sample_hz_,
                          mag_sample_hz_);
    }

    //nmea_xdr_freq(SEA_STATE_NMEA_TALKER, wave_hz_);
    nmea_rot(SEA_STATE_NMEA_TALKER, rot_dpm_filt_, heading_valid_);

#else

  #if ARDUINO_PLOTTER

    Serial.printf("HrawCm:%+.3f\tHbaselineCm:%+.3f\tHwaveCleanCm:%+.3f\tFreqmHz:%+.3f\tTheta:%+.3f\n",
                  static_cast<double>(heave_raw_m_ * 100.0f),
                  static_cast<double>(heave_baseline_m_ * 100.0f),
                  static_cast<double>(heave_wave_clean_m_ * 100.0f),
                  static_cast<double>(wave_hz_ * 1000.0f),
                  static_cast<double>(nlo_theta_));

  #else

    Serial.printf(
        "hdg=%.2f yaw=%.2f valid=%d init=%d boot=%d roll=%.2f pitch=%.2f "
        "|m|=%.1f magUsed=%d magField=%d magFresh=%d "
        "heave=%.3f heaveClean=%.3f vz=%.3f az=%.3f "
        "frq=%.3f conf=%.2f lock=%d theta=%.3f gbias=%.5f dt_ms=%.2f\n",
        static_cast<double>(heading_mag_deg_),
        static_cast<double>(yaw_deg_),
        static_cast<int>(heading_valid_),
        static_cast<int>(nlo_initialized_),
        static_cast<int>(nlo_bootstrapping_),
        static_cast<double>(roll_deg_),
        static_cast<double>(pitch_deg_),
        static_cast<double>(mag_norm_uT_),
        static_cast<int>(mag_used_),
        static_cast<int>(mag_field_sane_),
        static_cast<int>(mag_fresh_),
        static_cast<double>(heave_raw_m_),
        static_cast<double>(heave_wave_clean_m_),
        static_cast<double>(heave_speed_mps_),
        static_cast<double>(heave_accel_mps2_),
        static_cast<double>(wave_hz_),
        static_cast<double>(wave_conf_),
        static_cast<int>(wave_locked_),
        static_cast<double>(nlo_theta_),
        static_cast<double>(gyro_bias_norm_),
        static_cast<double>(dt_ * 1000.0f));

  #endif

#endif
  }
};

static FusionApp g_app;

void setup() {
  g_app.begin();
}

void loop() {
  g_app.tick();
}
