/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R AdaptiveVerticalPIIMahony (+ optional IMU Calibration Wizard)

  Produces:
    - compass heading
    - AHRS roll/pitch/yaw
    - heave Z displacement
    - wave envelope estimate
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
#include "ahrs/FrameConversions.h"
#include "detrend/AdaptiveWaveDetrender.h"
#include "nmea/NmeaCompass.h"
#include "pii_observer/AdaptiveVerticalPIIMahony.h"

#ifndef SEA_STATE_UI_DEFAULT_GRAPHICS
  #define SEA_STATE_UI_DEFAULT_GRAPHICS 1
#endif

#ifndef SEA_STATE_SERIAL_NMEA
  #define SEA_STATE_SERIAL_NMEA 0
#endif

#ifndef SEA_STATE_NMEA_TALKER
  #define SEA_STATE_NMEA_TALKER "II"
#endif

/*
  The sketch reports MAGNETIC heading.

  NMEA HDM remains magnetic by definition.
*/
#ifndef SEA_STATE_OUTPUT_TRUE_HEADING
  #define SEA_STATE_OUTPUT_TRUE_HEADING 0
#endif

#ifndef SEA_STATE_MAG_DECLINATION_DEG
  #define SEA_STATE_MAG_DECLINATION_DEG 0.0f
#endif

static constexpr float LOOP_HZ = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = static_cast<uint32_t>(1000000.0f / LOOP_HZ);

static constexpr uint32_t UI_REFRESH_MS = 100;
static constexpr uint32_t DEBUG_SERIAL_MS = 100;
static constexpr uint32_t NMEA_SERIAL_MS = 80;

static constexpr float ROT_BIAS_TAU_S = 5.0f;
static constexpr float ROT_STILL_G_TOL_FRAC = 0.12f;
static constexpr float ROT_STILL_GYRO_RAD_S = 0.15f;

/*
  Magnetometer sanity gates.

  Earth's field is commonly around 25..65 uT depending on location.
  20..80 uT is a practical runtime gate that rejects many local disturbances
  while still being tolerant enough for normal use.
*/
static constexpr float MAG_FIELD_MIN_UT = 20.0f;
static constexpr float MAG_FIELD_MAX_UT = 80.0f;

/*
  Do not feed magnetometer to Mahony at 200 Hz if the physical mag update rate
  is slower / repeated. This also reduces yaw twitch from noisy mag data.
*/
static constexpr uint32_t MAG_UPDATE_SPACING_MS = 35u;

/*
  Heading remains valid for a short time after last accepted mag correction.
  Between accepted mag updates, IMU-only propagation continues.
*/
static constexpr uint32_t HEADING_MAG_TIMEOUT_MS = 2000u;

using namespace atoms3r_ical;
using Vector3f = Eigen::Vector3f;

static inline float clampf_(float x, float lo, float hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

static inline float wrap360_(float deg) {
  while (deg < 0.0f) deg += 360.0f;
  while (deg >= 360.0f) deg -= 360.0f;
  return deg;
}

static inline float outputHeadingFromMagnetic_(float magnetic_deg) {
#if SEA_STATE_OUTPUT_TRUE_HEADING
  return wrap360_(magnetic_deg + SEA_STATE_MAG_DECLINATION_DEG);
#else
  return magnetic_deg;
#endif
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
      if (saved) Serial.println("[BOOT] Wizard saved calibration.");
      else Serial.println("[BOOT] Wizard did not save calibration. Running with raw values.");
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

    if (tap_count_ > 0 && static_cast<int32_t>(millis() - tap_deadline_ms_) > 0) {
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
  using Fusion = marine_obs::AdaptiveVerticalPIIMahony<float, true, TrackerType::PLL>;

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

  float heading_mag_deg_ = 0.0f;
  float heading_deg_ = 0.0f;
  bool heading_valid_ = false;

  float heave_m_ = 0.0f;
  float wave_envelope_m_ = 0.0f;
  float wave_hz_ = FREQ_GUESS;
  float heave_raw_m_ = 0.0f;
  float heave_baseline_m_ = 0.0f;
  float heave_wave_raw_m_ = 0.0f;
  float heave_wave_clean_m_ = 0.0f;

  bool mag_ok_ = false;
  bool mag_field_sane_ = false;
  bool mag_fresh_ = false;
  bool mag_used_ = false;
  float mag_norm_uT_ = NAN;

  uint16_t stale_frame_count_ = 0;
  bool have_last_sample_us_ = false;
  uint32_t last_sample_us_ = 0;

  bool rot_inited_ = false;
  float rot_dpm_filt_ = 0.0f;
  bool gyro_bias_ok_ = false;
  Vector3f gyro_bias_ema_ = Vector3f::Zero();
  AdaptiveWaveDetrender z_detrender_{};

  /*
    IMPORTANT FRAME RULE

    The gyro, accel, and mag must all enter Mahony in the SAME body frame.

    Previous code used:
        accel/gyro: ned_to_mahony_body_()
        mag:        ned_to_mahony_mag_() with different signs

    That makes the magnetic reference inconsistent with the inertial frame and
    causes heading offset / tilt-dependent heading drift.

    So magnetometer now uses the same mapping as accel and gyro.
  */
  static Vector3f ned_to_mahony_body_(const Vector3f& v_ned) {
    const Vector3f v_zu = ned_to_zu(v_ned);
    return Vector3f(-v_zu.x(), -v_zu.y(), v_zu.z());
  }

  static Vector3f ned_to_mahony_mag_(const Vector3f& v_ned) {
    return ned_to_mahony_body_(v_ned);
  }

  static void waitForNextLoopTick_(uint32_t loop_start_us) {
    const uint32_t elapsed_us = micros() - loop_start_us;
    if (elapsed_us < LOOP_PERIOD_US) delayMicroseconds(LOOP_PERIOD_US - elapsed_us);
  }

  float nominalFusionDt_() const { return 1.0f / LOOP_HZ; }

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

    /*
      Protect Mahony and PII from occasional scheduling hiccups.
      The loop is nominally 200 Hz, so normal dt is around 0.005 s.
    */
    if (dt_s > 0.05f) return dt_nom;

    return dt_s;
  }

  void reloadBlobAndRuntime_() {
    have_blob_ = store_.load(blob_);
    if (!have_blob_) std::memset(&blob_, 0, sizeof(blob_));
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
    have_last_sample_us_ = false;
    last_sample_us_ = 0;
  }

  void resetFusion_() {
    Fusion::Config cfg{};
    cfg.gravity_mps2 = g_std;
    cfg.use_mag = true;

    cfg.core.observer.r          = 0.150f;
    cfg.core.observer.tau_a      = 0.68f;
    cfg.core.observer.tau_d      = 49.0f;
    cfg.core.observer.kb         = 2.5e-5f;
    cfg.core.observer.lambda_b   = 3.0e-3f;
    cfg.core.observer.bias_limit = 0.12f;
    cfg.core.observer.a_f_limit  = 50.0f;
    cfg.core.observer.v_limit    = 50.0f;
    cfg.core.observer.p_limit    = 20.0f;
    cfg.core.observer.S_limit    = 200.0f;
    cfg.core.observer.d_limit    = 20.0f;

    cfg.core.adaptation.enabled          = true;
    cfg.core.adaptation.min_confidence   = 0.22f;
    cfg.core.adaptation.f_disp_ref_hz    = 0.12f;
    cfg.core.adaptation.sigma_a_ref      = 0.95f;
    cfg.core.adaptation.input_smooth_tau = 4.5f;
    cfg.core.adaptation.param_smooth_tau = 7.5f;
    cfg.core.adaptation.r_freq_exp       = 0.28f;
    cfg.core.adaptation.r_sigma_exp      = 0.02f;
    cfg.core.adaptation.tau_a_freq_exp   = -0.40f;
    cfg.core.adaptation.tau_a_sigma_exp  = -0.03f;
    cfg.core.adaptation.tau_d_freq_exp   = -0.03f;
    cfg.core.adaptation.tau_d_sigma_exp  = -0.01f;
    cfg.core.adaptation.kb_freq_exp      = 0.02f;
    cfg.core.adaptation.kb_sigma_exp     = 0.08f;
    cfg.core.adaptation.r_min            = 0.145f;
    cfg.core.adaptation.r_max            = 0.225f;
    cfg.core.adaptation.tau_a_min        = 0.50f;
    cfg.core.adaptation.tau_a_max        = 0.90f;
    cfg.core.adaptation.tau_d_min        = 44.0f;
    cfg.core.adaptation.tau_d_max        = 58.0f;
    cfg.core.adaptation.kb_min           = 5e-6f;
    cfg.core.adaptation.kb_max           = 6e-5f;

    cfg.core.auto_schedule_from_accel_freq = true;
    cfg.core.auto_schedule_period_s = 0.50f;
    cfg.core.force_enable_adaptation_when_auto_schedule = true;
    cfg.core.fallback_confidence_floor = 0.52f;
    cfg.core.fallback_confidence_when_locked = 0.82f;
    cfg.core.coarse_schedule_blend = 0.48f;
    cfg.core.coarse_schedule_confidence_floor = 0.62f;

    /*
      These are only initial/reset values. The wrapper may adapt them.
      If the header's adaptMahonyGains_() calls mahony_.init() every sample,
      patch that header too so it only assigns twoKp/twoKi each sample.
    */
    cfg.mahony_twoKp = 1.40f;
    cfg.mahony_twoKi = 0.060f;
    cfg.adapt_mahony_gains = true;

    fusion_.configure(cfg);
    fusion_.reset();

    rot_inited_ = false;
    rot_dpm_filt_ = 0.0f;
    gyro_bias_ok_ = false;
    gyro_bias_ema_.setZero();

    roll_deg_ = 0.0f;
    pitch_deg_ = 0.0f;
    heading_mag_deg_ = 0.0f;
    heading_deg_ = 0.0f;
    heading_valid_ = false;

    mag_gate_last_ms_ = 0;
    last_mag_correction_ms_ = 0;
    mag_ok_ = false;
    mag_field_sane_ = false;
    mag_fresh_ = false;
    mag_used_ = false;
    mag_norm_uT_ = NAN;

    stale_frame_count_ = 0;
    have_last_sample_us_ = false;
    last_sample_us_ = 0;

    heave_m_ = 0.0f;
    wave_envelope_m_ = 0.0f;
    wave_hz_ = FREQ_GUESS;
    heave_raw_m_ = 0.0f;
    heave_baseline_m_ = 0.0f;
    heave_wave_raw_m_ = 0.0f;
    heave_wave_clean_m_ = 0.0f;

    AdaptiveWaveDetrender::Config dcfg{};
    dcfg.init_wave_freq_hz = FREQ_GUESS;
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
    dcfg.enable_wave_cleanup    = true;
    dcfg.cleanup_cutoff_fraction = 1.0f;
    dcfg.min_cleanup_cutoff_hz  = 0.003f;
    dcfg.max_cleanup_cutoff_hz  = 0.50f;
    dcfg.cleanup_stages         = 2;
    dcfg.min_dt_s = 1.0e-4f;
    dcfg.max_dt_s = 0.25f;
    dcfg.output_abs_limit = 0.0f;

    z_detrender_.setConfig(dcfg);
    z_detrender_.reset(0.0f);
  }

#if SEA_STATE_ENABLE_WIZARD
  void handleErase_() {
    if (!ui_.eraseConfirm()) return;
    store_.erase();
    reloadBlobAndRuntime_();
    reinitImu_();
    resetFusion_();
  }

  void handleRunWizard_() {
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

    a_raw_norm_ = s.a.norm();

    a_cal_ = runtime_.applyAccel(s.a, tempC);
    w_cal_ = runtime_.applyGyro(s.w, tempC);
    m_cal_ = runtime_.applyMag(s.m);

    mag_norm_uT_ = m_cal_.norm();

    mag_ok_ =
        std::isfinite(mag_norm_uT_) &&
        mag_norm_uT_ > 5.0f &&
        mag_norm_uT_ < 200.0f;

    mag_field_sane_ =
        std::isfinite(mag_norm_uT_) &&
        mag_norm_uT_ >= MAG_FIELD_MIN_UT &&
        mag_norm_uT_ <= MAG_FIELD_MAX_UT;

    /*
      This is a rate gate for accepted mag corrections, not proof that the
      sensor produced a physically new sample. It prevents repeatedly feeding
      the same / noisy magnetometer vector at the 200 Hz IMU loop rate.
    */
    mag_fresh_ = updateMagFreshGate_(mag_ok_ && mag_field_sane_, now_ms);

    const Vector3f gyr_body_m = ned_to_mahony_body_(w_cal_);
    const Vector3f acc_body_m = ned_to_mahony_body_(a_cal_);

    const bool mag_usable = mag_ok_ && mag_field_sane_ && mag_fresh_;
    mag_used_ = mag_usable;

    if (mag_usable) {
      /*
        Critical fix:
        mag must use the same body-frame convention as gyro and accel.
      */
      const Vector3f mag_body_m = ned_to_mahony_mag_(m_cal_);

      fusion_.updateIMUMag(gyr_body_m.x(), gyr_body_m.y(), gyr_body_m.z(),
                           acc_body_m.x(), acc_body_m.y(), acc_body_m.z(),
                           mag_body_m.x(), mag_body_m.y(), mag_body_m.z(),
                           dt_);

      last_mag_correction_ms_ = now_ms;
    } else {
      fusion_.updateIMU(gyr_body_m.x(), gyr_body_m.y(), gyr_body_m.z(),
                        acc_body_m.x(), acc_body_m.y(), acc_body_m.z(),
                        dt_);
    }

    const auto q_wb = fusion_.quaternionWorldToBody();
    const Quaternionf q_wb_zu(q_wb.w, q_wb.x, q_wb.y, q_wb.z);

    float heading_from_q_deg = 0.0f;
    quat_wb_zu_to_euler_nautical(q_wb_zu, roll_deg_, pitch_deg_, heading_from_q_deg);

    heading_mag_deg_ = wrap360_(heading_from_q_deg);
    heading_deg_ = outputHeadingFromMagnetic_(heading_mag_deg_);

    heading_valid_ =
        last_mag_correction_ms_ != 0 &&
        (now_ms - last_mag_correction_ms_) <= HEADING_MAG_TIMEOUT_MS;

    const bool still =
        (fabsf(a_cal_.norm() - g_std) < ROT_STILL_G_TOL_FRAC * g_std) &&
        (w_cal_.norm() < ROT_STILL_GYRO_RAD_S);

    if (still) {
      const float alpha_b = 1.0f - expf(-dt_ / ROT_BIAS_TAU_S);

      if (!gyro_bias_ok_) {
        gyro_bias_ok_ = true;
        gyro_bias_ema_ = w_cal_;
      } else {
        gyro_bias_ema_ += alpha_b * (w_cal_ - gyro_bias_ema_);
      }
    }

    Vector3f w_use = w_cal_;
    if (gyro_bias_ok_) w_use -= gyro_bias_ema_;

    float rot_dpm_meas = w_use.z() * RAD_TO_DEG * 60.0f;
    rot_dpm_meas = clampf_(rot_dpm_meas, -720.0f, 720.0f);

    const float tau_rot = 1.5f;
    const float alpha_r = 1.0f - expf(-dt_ / tau_rot);

    if (!rot_inited_) {
      rot_inited_ = true;
      rot_dpm_filt_ = rot_dpm_meas;
    } else {
      rot_dpm_filt_ += alpha_r * (rot_dpm_meas - rot_dpm_filt_);
    }

    const auto hs = fusion_.snapshot();

    heave_m_ = fusion_.displacement();
    heave_raw_m_ = heave_m_;
    wave_hz_ = hs.core.accel_freq_hz;

    const float omega = (wave_hz_ > 1e-6f) ? (2.0f * PI * wave_hz_) : NAN;
    const float sigma = hs.core.accel_sigma;

    if (std::isfinite(omega) && std::isfinite(sigma) && omega > 1e-6f) {
      wave_envelope_m_ = sigma / (omega * omega);
    } else {
      wave_envelope_m_ = 0.0f;
    }

    const auto z_det =
        z_detrender_.update(heave_raw_m_, dt_, wave_hz_, wave_hz_ > 1e-6f);

    heave_baseline_m_ = z_det.baseline_slow;
    heave_wave_raw_m_ = z_det.wave_raw;
    heave_wave_clean_m_ = z_det.wave_clean;
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
    M5.Display.printf("HDG: TRUE decl=%+.1f\n", static_cast<double>(SEA_STATE_MAG_DECLINATION_DEG));
#else
    ui_.line("HDG: MAG");
#endif

    ui_.line("");
    ui_.line("Fusion: PII");
  }

  void updateUI_() {
#if SEA_STATE_ENABLE_WIZARD
    if (tap_count_ > 0) return;
#endif

    const uint32_t now_ms = millis();
    if (now_ms - last_ui_ms_ < UI_REFRESH_MS) return;
    last_ui_ms_ = now_ms;

    if (use_graphics_ && compass_ui_ready_) {
      ui_.setReadRotation();

      const bool tiltWarn =
          (fabsf(roll_deg_) > 35.0f) ||
          (fabsf(pitch_deg_) > 35.0f);

      const float hdg_draw = heading_valid_ ? heading_deg_ : 0.0f;

      compass_ui_.draw(hdg_draw,
                       heading_valid_ && mag_ok_ && mag_field_sane_,
                       mag_norm_uT_,
                       tiltWarn);
      return;
    }

    ui_.setReadRotation();
    ui_.title("COMPASS");

#if SEA_STATE_OUTPUT_TRUE_HEADING
    M5.Display.printf("HDG:%7.1f T %s\n",
                      static_cast<double>(heading_deg_),
                      heading_valid_ ? "deg" : "WAIT");
    M5.Display.printf("MAG:%7.1f M\n", static_cast<double>(heading_mag_deg_));
#else
    M5.Display.printf("HDG:%7.1f M %s\n",
                      static_cast<double>(heading_mag_deg_),
                      heading_valid_ ? "deg" : "WAIT");
#endif

    M5.Display.printf("ROL:%7.1f deg\n", static_cast<double>(roll_deg_));
    M5.Display.printf("PIT:%7.1f deg\n", static_cast<double>(pitch_deg_));
    M5.Display.printf("HEV:%7.3f m\n", static_cast<double>(heave_raw_m_));
    M5.Display.printf("ENV:%7.3f m\n", static_cast<double>(wave_envelope_m_));
    M5.Display.printf("FRQ:%7.3f Hz\n", static_cast<double>(wave_hz_));

    M5.Display.printf("MAG:%s %s %s\n",
                      mag_ok_ ? "OK " : "BAD",
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

    /*
      HDM is magnetic heading. Do not send true-heading-corrected value here.
    */
    if (heading_valid_) nmea_hdm(SEA_STATE_NMEA_TALKER, heading_mag_deg_);

    nmea_xdr_pitch_roll(SEA_STATE_NMEA_TALKER, pitch_deg_, roll_deg_);
    nmea_xdr_heave(SEA_STATE_NMEA_TALKER, heave_wave_clean_m_);
    nmea_xdr_freq(SEA_STATE_NMEA_TALKER, wave_hz_);
    nmea_rot(SEA_STATE_NMEA_TALKER, rot_dpm_filt_, heading_valid_);

#else

  #if ARDUINO_PLOTTER
    Serial.printf("HrawCm:%+.3f\tHwaveEnvelopeCm:%+.3f\tHwaveCleanCm:%+.3f\n",
                  static_cast<double>(heave_raw_m_ * 100.0f),
                  static_cast<double>(wave_envelope_m_ * 100.0f),
                  static_cast<double>(heave_wave_clean_m_ * 100.0f));
  #else
    Serial.printf("hdg=%.2f mag=%.2f valid=%d roll=%.2f pitch=%.2f |m|=%.1f magUsed=%d magField=%d heave=%.3f env=%.3f frq=%.3f\n",
                  static_cast<double>(heading_deg_),
                  static_cast<double>(heading_mag_deg_),
                  static_cast<int>(heading_valid_),
                  static_cast<double>(roll_deg_),
                  static_cast<double>(pitch_deg_),
                  static_cast<double>(mag_norm_uT_),
                  static_cast<int>(mag_used_),
                  static_cast<int>(mag_field_sane_),
                  static_cast<double>(heave_m_),
                  static_cast<double>(wave_envelope_m_),
                  static_cast<double>(wave_hz_));
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
