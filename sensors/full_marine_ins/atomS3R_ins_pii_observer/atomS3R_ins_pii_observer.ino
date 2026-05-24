/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R AdaptiveVerticalPIIMahony (+ optional IMU Calibration Wizard)

  Produces:
    - magnetic compass heading from Mahony yaw
    - AHRS roll/pitch/yaw
    - heave Z displacement
    - wave envelope estimate

  REAL DEVICE MAHONY CONVENTION

    Calibrated device/body frame from runtime_ / OU_II:
      [N, E, D]
      x = forward
      y = right
      z = down

    Mahony frame used here:
      [N, -E, -D]
      x = forward
      y = left
      z = up

    Therefore ALL THREE sensors use the same mapper:
      gyro  [N,E,D] -> [N,-E,-D]
      accel [N,E,D] -> [N,-E,-D]
      mag   [N,E,D] -> [N,-E,-D]

    At level, pointed magnetic north:
      accel_ned = [0, 0, -g] -> [0, 0, +g]
      mag_ned   = [+H, 0, +V] -> [+H, 0, -V]

    Mahony yaw is CCW-positive in Z-up frame.
    Compass heading is CW-positive.

      magnetic_heading = -MahonyYaw + user_offset

    Startup:
      This sketch seeds Mahony quaternion once from accel+mag before normal
      updateIMUMag(), so yaw does not take minutes to converge from identity.
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
#include "nmea/NmeaCompass.h"
#include "pii_observer/AdaptiveVerticalPIIMahony.h"

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

static constexpr float APP_G_STD = atoms3r_ical::ImuCalCfg::g_std;
static constexpr float APP_FREQ_GUESS = 0.30f;

static constexpr float LOOP_HZ = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = static_cast<uint32_t>(1000000.0f / LOOP_HZ);

static constexpr uint32_t UI_REFRESH_MS = 100;
static constexpr uint32_t DEBUG_SERIAL_MS = 100;
static constexpr uint32_t NMEA_SERIAL_MS = 80;

static constexpr float ROT_BIAS_TAU_S = 5.0f;
static constexpr float ROT_STILL_G_TOL_FRAC = 0.12f;
static constexpr float ROT_STILL_GYRO_RAD_S = 0.15f;

static constexpr float MAG_PRESENT_MIN_UT = 5.0f;
static constexpr float MAG_PRESENT_MAX_UT = 200.0f;

static constexpr float MAG_FIELD_MIN_UT = 20.0f;
static constexpr float MAG_FIELD_MAX_UT = 80.0f;

static constexpr float MAHONY_SEED_G_TOL_FRAC = 0.20f;

static constexpr uint32_t MAG_UPDATE_SPACING_MS = 35u;
static constexpr uint32_t HEADING_MAG_TIMEOUT_MS = 2000u;

using namespace atoms3r_ical;
using Vector3f = Eigen::Vector3f;
using Matrix3f = Eigen::Matrix3f;
using Quaternionf = Eigen::Quaternionf;

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
  float yaw_mahony_deg_ = 0.0f;

  float heading_mag_deg_ = 0.0f;
  float heading_deg_ = 0.0f;
  bool heading_valid_ = false;

  float heave_m_ = 0.0f;
  float wave_envelope_m_ = 0.0f;
  float wave_hz_ = APP_FREQ_GUESS;

  float heave_raw_m_ = 0.0f;
  float heave_baseline_m_ = 0.0f;
  float heave_wave_raw_m_ = 0.0f;
  float heave_wave_clean_m_ = 0.0f;

  bool mag_present_ = false;
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

  bool mahony_seeded_ = false;

  AdaptiveWaveDetrender z_detrender_{};

 private:
  static Vector3f ned_to_mahony_body_(const Vector3f& v_ned) {
    return Vector3f(v_ned.x(), -v_ned.y(), -v_ned.z());
  }

  bool seedMahonyFromAccelMag_(const Vector3f& acc_body_m,
                               const Vector3f& mag_body_m) {
    const float an = acc_body_m.norm();
    const float mn = mag_body_m.norm();

    if (!(an > 1e-6f) || !(mn > 1e-6f) ||
        !std::isfinite(an) || !std::isfinite(mn)) {
      return false;
    }

    Vector3f up_b = acc_body_m / an;

    Vector3f north_b = mag_body_m - up_b * mag_body_m.dot(up_b);
    const float nn0 = north_b.norm();
    if (!(nn0 > 1e-6f) || !std::isfinite(nn0)) {
      return false;
    }
    north_b /= nn0;

    Vector3f west_b = up_b.cross(north_b);
    const float wn = west_b.norm();
    if (!(wn > 1e-6f) || !std::isfinite(wn)) {
      return false;
    }
    west_b /= wn;

    north_b = west_b.cross(up_b);
    const float nn = north_b.norm();
    if (!(nn > 1e-6f) || !std::isfinite(nn)) {
      return false;
    }
    north_b /= nn;

    Matrix3f C_wb;
    C_wb.col(0) = north_b;
    C_wb.col(1) = west_b;
    C_wb.col(2) = up_b;

    const Matrix3f C_bw = C_wb.transpose();

    Quaternionf q(C_bw);
    const float qn = q.norm();
    if (!(qn > 1e-6f) || !std::isfinite(qn)) {
      return false;
    }
    q.normalize();

    auto& m = fusion_.mahonyState();

    m.q0 = q.w();
    m.q1 = q.x();
    m.q2 = q.y();
    m.q3 = q.z();

    m.integralFBx = 0.0f;
    m.integralFBy = 0.0f;
    m.integralFBz = 0.0f;

    return true;
  }

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

    /*
      Keep filter settings at AdaptiveVerticalPIIMahony defaults.
      Only app/hardware-specific fields are set here.
    */
    cfg.gravity_mps2 = APP_G_STD;
    cfg.use_mag = true;

    /*
      OU-style envelope estimator lives in AdaptiveVerticalPII.
      These are the same defaults as the header, kept explicit here so the
      sketch behavior is documented and not hidden in the .ino math.
    */
    cfg.core.envelope.enabled = true;
    cfg.core.envelope.acc_noise_floor_sigma = 0.12f;
    cfg.core.envelope.tuner_K_periods = 2.0f;
    cfg.core.envelope.tuner_tau_freq_s = 1.0f;
    cfg.core.envelope.tau_coeff = 1.38f;
    cfg.core.envelope.sigma_coeff = 0.90f;
    cfg.core.envelope.adapt_tau_s = 1.80f;
    cfg.core.envelope.tau_min_s = 0.02f;
    cfg.core.envelope.tau_max_s = 3.00f;
    cfg.core.envelope.sigma_max = 6.00f;

    fusion_.configure(cfg);
    fusion_.reset();

    rot_inited_ = false;
    rot_dpm_filt_ = 0.0f;
    gyro_bias_ok_ = false;
    gyro_bias_ema_.setZero();
    mahony_seeded_ = false;

    roll_deg_ = 0.0f;
    pitch_deg_ = 0.0f;
    yaw_mahony_deg_ = 0.0f;

    heading_mag_deg_ = 0.0f;
    heading_deg_ = 0.0f;
    heading_valid_ = false;

    mag_gate_last_ms_ = 0;
    last_mag_correction_ms_ = 0;

    mag_present_ = false;
    mag_field_sane_ = false;
    mag_fresh_ = false;
    mag_used_ = false;
    mag_norm_uT_ = NAN;

    stale_frame_count_ = 0;
    have_last_sample_us_ = false;
    last_sample_us_ = 0;

    heave_m_ = 0.0f;
    wave_envelope_m_ = 0.0f;
    wave_hz_ = APP_FREQ_GUESS;

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

    mag_present_ =
        std::isfinite(mag_norm_uT_) &&
        mag_norm_uT_ >= MAG_PRESENT_MIN_UT &&
        mag_norm_uT_ <= MAG_PRESENT_MAX_UT;

    mag_field_sane_ =
        std::isfinite(mag_norm_uT_) &&
        mag_norm_uT_ >= MAG_FIELD_MIN_UT &&
        mag_norm_uT_ <= MAG_FIELD_MAX_UT;

    mag_fresh_ = updateMagFreshGate_(mag_present_, now_ms);

#if SEA_STATE_USE_STRICT_MAG_FIELD_GATE
    const bool mag_usable = mag_present_ && mag_field_sane_;
#else
    const bool mag_usable = mag_present_;
#endif

    mag_used_ = mag_usable;

    const Vector3f gyr_body_m = ned_to_mahony_body_(w_cal_);
    const Vector3f acc_body_m = ned_to_mahony_body_(a_cal_);
    const Vector3f mag_body_m = ned_to_mahony_body_(m_cal_);

    if (!mahony_seeded_ && mag_usable) {
      const float a_norm = a_cal_.norm();
      const bool accel_seed_ok =
          std::isfinite(a_norm) &&
          fabsf(a_norm - APP_G_STD) < MAHONY_SEED_G_TOL_FRAC * APP_G_STD;

      if (accel_seed_ok) {
        mahony_seeded_ = seedMahonyFromAccelMag_(acc_body_m, mag_body_m);
      }
    }

    if (mag_usable) {
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

    roll_deg_ = -fusion_.rollDeg();
    pitch_deg_ = -fusion_.pitchDeg();
    yaw_mahony_deg_ = fusion_.yawDeg();

    if (mag_usable) {
      heading_mag_deg_ =
          wrap360_(-yaw_mahony_deg_ + SEA_STATE_MAG_HEADING_USER_OFFSET_DEG);
      heading_deg_ = outputHeadingFromMagnetic_(heading_mag_deg_);
    }

    heading_valid_ =
        last_mag_correction_ms_ != 0 &&
        (now_ms - last_mag_correction_ms_) <= HEADING_MAG_TIMEOUT_MS;

    const bool still =
        (fabsf(a_cal_.norm() - APP_G_STD) < ROT_STILL_G_TOL_FRAC * APP_G_STD) &&
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

    const float tau_rot = 0.1f;
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

    /*
      Use core-owned envelope frequency when ready. Fall back to the raw accel
      tracker frequency during startup.
    */
    const float env_freq_hz = fusion_.envelopeFrequencyHz();
    if (fusion_.envelopeReady() &&
        std::isfinite(env_freq_hz) &&
        env_freq_hz > 1e-6f) {
      wave_hz_ = env_freq_hz;
    } else {
      wave_hz_ = hs.core.accel_freq_hz;
    }

    /*
      Use the core-owned OU-style envelope estimate.

      This replaces:
          sigma / omega^2

      That old formula was too large because it used raw accel sigma directly
      as a displacement amplitude proxy. The core now owns the same kind of
      variance/noise-floor/tau^2 scaling model used by the OU filter.
    */
    wave_envelope_m_ = fusion_.displacementScale();
    if (!(std::isfinite(wave_envelope_m_) && wave_envelope_m_ >= 0.0f)) {
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

    ui_.line("Fusion: PII DEFAULTS");
    ui_.line("Startup: ACC+MAG seed");
    ui_.line("ENV: core OU-style");
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
                       heading_valid_ && mag_used_,
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
    M5.Display.printf("MAG:%7.1f M\n",
                      static_cast<double>(heading_mag_deg_));
#else
    M5.Display.printf("HDG:%7.1f M %s\n",
                      static_cast<double>(heading_mag_deg_),
                      heading_valid_ ? "deg" : "WAIT");
#endif

    M5.Display.printf("YAW:%7.1f deg\n", static_cast<double>(yaw_mahony_deg_));
    M5.Display.printf("ROL:%7.1f deg\n", static_cast<double>(roll_deg_));
    M5.Display.printf("PIT:%7.1f deg\n", static_cast<double>(pitch_deg_));

    M5.Display.printf("SEED:%s\n", mahony_seeded_ ? "YES" : "NO");

    M5.Display.printf("HEV:%7.3f m\n", static_cast<double>(heave_raw_m_));
    M5.Display.printf("ENV:%7.3f m\n", static_cast<double>(wave_envelope_m_));
    M5.Display.printf("FRQ:%7.3f Hz\n", static_cast<double>(wave_hz_));

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
    //nmea_xdr_freq(SEA_STATE_NMEA_TALKER, wave_hz_);
    nmea_rot(SEA_STATE_NMEA_TALKER, rot_dpm_filt_, heading_valid_);

#else

  #if ARDUINO_PLOTTER

    Serial.printf("HrawCm:%+.3f\tHwaveEnvelopeCm:%+.3f\tHwaveCleanCm:%+.3f\n",
                  static_cast<double>(heave_raw_m_ * 100.0f),
                  static_cast<double>(wave_envelope_m_ * 100.0f),
                  static_cast<double>(heave_wave_clean_m_ * 100.0f));

  #else

    Serial.printf(
        "hdg=%.2f yaw=%.2f valid=%d seed=%d roll=%.2f pitch=%.2f |m|=%.1f "
        "magUsed=%d magField=%d magFresh=%d heave=%.3f env=%.3f frq=%.3f "
        "envReady=%d envVar=%.5f envTau=%.3f envSigma=%.3f\n",
        static_cast<double>(heading_mag_deg_),
        static_cast<double>(yaw_mahony_deg_),
        static_cast<int>(heading_valid_),
        static_cast<int>(mahony_seeded_),
        static_cast<double>(roll_deg_),
        static_cast<double>(pitch_deg_),
        static_cast<double>(mag_norm_uT_),
        static_cast<int>(mag_used_),
        static_cast<int>(mag_field_sane_),
        static_cast<int>(mag_fresh_),
        static_cast<double>(heave_m_),
        static_cast<double>(wave_envelope_m_),
        static_cast<double>(wave_hz_),
        static_cast<int>(fusion_.envelopeReady()),
        static_cast<double>(fusion_.envelopeAccelVariance()),
        static_cast<double>(fusion_.envelopeTauApplied()),
        static_cast<double>(fusion_.envelopeSigmaApplied()));

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
