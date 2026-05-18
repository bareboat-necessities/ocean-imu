/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R SeaStateFusion_OU_III (+ optional IMU Calibration Wizard)

    - filter learns tilt during startup
    - one-shot magnetic north lock inside SeaStateFusion_OU_III
    - filter yaw is the primary heading output after north-lock
    - tilt-compensated magnetic heading is kept only for diagnostics

  Assumptions:
    - fusion_.raw().mekf().quaternion_boat() returns BODY->WORLD quaternion (q_bw), world frame is NED (+Z down).
    - runtime_.apply* already includes temperature compensation.
*/

#include <Arduino.h>

#ifndef SEA_STATE_ENABLE_WIZARD
  #define SEA_STATE_ENABLE_WIZARD 1
#endif

#define ARDUINO_PLOTTER 1

#include <M5Unified.h>
#include <cmath>
#include <algorithm>
#include <cstring>

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
  #define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif
#include <ArduinoEigenDense.h>

#include <ArduinoOceanImu.h>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_M5Ui.h"
#if SEA_STATE_ENABLE_WIZARD
  #include "AtomS3R/ImuCalWizardRunner.h"
#endif
#include "AtomS3R/AtomS3R_CompassUI.h"
#include "nmea/NmeaCompass.h"

#ifndef SEA_STATE_UI_DEFAULT_GRAPHICS
  #define SEA_STATE_UI_DEFAULT_GRAPHICS 1
#endif

#ifndef SEA_STATE_SERIAL_NMEA
  #define SEA_STATE_SERIAL_NMEA 0
#endif

#ifndef SEA_STATE_NMEA_TALKER
  #define SEA_STATE_NMEA_TALKER "II"
#endif

constexpr float g_std      = atoms3r_ical::ImuCalCfg::g_std;
constexpr float FREQ_GUESS = 0.3f;

#define ZERO_CROSSINGS_SCALE          1.0f
#define ZERO_CROSSINGS_DEBOUNCE_TIME  0.12f
#define ZERO_CROSSINGS_STEEPNESS_TIME 0.21f

#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

static constexpr float    LOOP_HZ        = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = static_cast<uint32_t>(1000000.0f / LOOP_HZ);

static constexpr uint32_t UI_REFRESH_MS   = 100;
static constexpr uint32_t DEBUG_SERIAL_MS = 100;
static constexpr uint32_t NMEA_SERIAL_MS  = 50;

static constexpr float ROT_BIAS_TAU_S       = 5.0f;
static constexpr float ROT_STILL_G_TOL_FRAC = 0.12f;
static constexpr float ROT_STILL_GYRO_RAD_S = 0.15f;

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

static inline float wrap180_(float deg) {
  while (deg <= -180.0f) deg += 360.0f;
  while (deg >   180.0f) deg -= 360.0f;
  return deg;
}

static inline Vector3f quatRotate_(const Eigen::Quaternionf& q, const Vector3f& v) {
  const Vector3f qv(q.x(), q.y(), q.z());
  const Vector3f t = 2.0f * qv.cross(v);
  return v + q.w() * t + qv.cross(t);
}

static inline bool rollPitchHeadingFromQuatBw_(
    const Eigen::Quaternionf& q_bw,
    float& roll_deg_out,
    float& pitch_deg_out,
    float& heading_deg_out) {
  Eigen::Quaternionf q = q_bw;
  const float nq = q.norm();
  if (!(nq > 1e-6f) || !std::isfinite(nq)) return false;
  q.normalize();

  const float x = q.x();
  const float y = q.y();
  const float z = q.z();
  const float w = q.w();

  const float siny_cosp = 2.0f * (w * z + x * y);
  const float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
  const float yaw = atan2f(siny_cosp, cosy_cosp);

  float sinp = 2.0f * (w * y - z * x);
  sinp = clampf_(sinp, -1.0f, 1.0f);
  const float pitch = asinf(sinp);

  const float sinr_cosp = 2.0f * (w * x + y * z);
  const float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
  const float roll = atan2f(sinr_cosp, cosr_cosp);

  roll_deg_out    = wrap180_(roll * RAD_TO_DEG);
  pitch_deg_out   = wrap180_(pitch * RAD_TO_DEG);
  heading_deg_out = wrap360_(yaw * RAD_TO_DEG);
  return true;
}

static inline bool magneticHeadingFromDownAndMagBody_(
    const Vector3f& down_b_unit,
    const Vector3f& mag_b_uT,
    float& heading_deg_out)
{
  const Vector3f FWD_B(1.0f, 0.0f, 0.0f);

  Vector3f d = down_b_unit;
  const float dn = d.norm();
  if (dn < 1e-6f) return false;
  d /= dn;

  Vector3f m = mag_b_uT;
  const float mn = m.norm();
  if (mn < 1e-6f) return false;
  m /= mn;

  Vector3f east_b = d.cross(m);
  const float en = east_b.norm();
  if (en < 1e-6f) return false;
  east_b /= en;

  Vector3f north_b = east_b.cross(d);
  const float nn = north_b.norm();
  if (nn < 1e-6f) return false;
  north_b /= nn;

  const float e = east_b.dot(FWD_B);
  const float n = north_b.dot(FWD_B);
  heading_deg_out = wrap360_(atan2f(e, n) * RAD_TO_DEG);
  return true;
}

class FusionApp {
public:
  FusionApp() = default;

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

    if (tap_count_ > 0 && static_cast<int32_t>(millis() - tap_deadline_ms_) > 0) {
      if (tap_count_ >= 3) handleErase_();
      else                 handleRunWizard_();
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
    } else {
      if (stale_frame_count_ < 0xFFFFu) stale_frame_count_++;
      if ((stale_frame_count_ % 50u) == 0u) {
        Serial.printf("[IMU] waiting: stale=%u\n", static_cast<unsigned>(stale_frame_count_));
      }
    }

    updateUI_();
    streamSerial_();
    waitForNextLoopTick_(loop_start_us);
  }

private:
  using UI = atoms3r_ical::M5Ui;
  using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

  bool use_graphics_ = (SEA_STATE_UI_DEFAULT_GRAPHICS != 0);

  CompassUI compass_ui_{};
  bool      compass_ui_ready_ = false;
  UI        ui_{};

  ImuCalStoreNvs store_{};
  bool           have_blob_ = false;
  ImuCalBlobV2   blob_{};
  RuntimeCals    runtime_{};

#if SEA_STATE_ENABLE_WIZARD
  int      tap_count_ = 0;
  uint32_t tap_deadline_ms_ = 0;
#endif

  uint32_t last_ui_ms_     = 0;
  uint32_t last_serial_ms_ = 0;

  Fusion fusion_{};

  uint32_t mag_gate_last_ms_ = 0;

  Vector3f a_cal_ = Vector3f::Zero();
  Vector3f w_cal_ = Vector3f::Zero();
  Vector3f m_cal_ = Vector3f::Constant(NAN);
  float    a_raw_norm_ = 0.0f;

  float dt_               = 0.0f;
  float roll_deg_         = 0.0f;
  float pitch_deg_        = 0.0f;
  float heading_deg_      = 0.0f;
  bool  heading_valid_    = false;
  float heave_m_          = 0.0f;
  float wave_envelope_m_  = 0.0f;
  float wave_hz_          = FREQ_GUESS;

  float heave_raw_m_        = 0.0f;
  float heave_baseline_m_   = 0.0f;
  float heave_wave_raw_m_   = 0.0f;
  float heave_wave_clean_m_ = 0.0f;

  bool  mag_ok_          = false;
  bool  mag_fresh_       = false;
  float mag_norm_uT_     = NAN;
  float heading_mag_deg_ = NAN;
  bool  heading_mag_ok_  = false;

  uint16_t stale_frame_count_ = 0;
  uint32_t last_skipped_total_ = 0;
  bool     have_last_sample_us_ = false;
  uint32_t last_sample_us_      = 0;

  bool     rot_inited_    = false;
  float    rot_dpm_filt_  = 0.0f;
  bool     gyro_bias_ok_  = false;
  Vector3f gyro_bias_ema_ = Vector3f::Zero();

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

  bool updateMagFreshGate_(bool mag_ok, uint32_t now_ms) {
    constexpr uint32_t kSampleSpacingMs = 35u;

    if (!mag_ok) {
      mag_gate_last_ms_ = 0;
      return false;
    }

    if (mag_gate_last_ms_ == 0) {
      mag_gate_last_ms_ = now_ms;
      return true;
    }

    if ((now_ms - mag_gate_last_ms_) < kSampleSpacingMs) {
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
    return dt_s;
  }

  uint32_t magPollMs_() const {
    const float use_hz = 25.0f;
    const float ms_f = 1000.0f / use_hz;
    if (!(ms_f > 0.0f) || !std::isfinite(ms_f)) {
      return 40u;
    }
    const uint32_t ms = static_cast<uint32_t>(ms_f + 0.5f);
    return std::max<uint32_t>(1u, ms);
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
    return false;
#else
    (void)boot_mode;
    return false;
#endif
  }

  void reinitImu_() {
    if (!M5.Imu.isEnabled()) {
      Serial.println("[BOOT] M5.Imu is not enabled");
      ui_.fail("IMU", "M5.Imu disabled");
      while (true) delay(100);
    }

    mag_gate_last_ms_ = 0;

    last_skipped_total_ = 0;
    have_last_sample_us_ = false;
    last_sample_us_      = 0;
  }

  void resetFusion_() {
    constexpr bool WANT_LINEAR_BLOCK = true;

    constexpr float IMU_TUNE_REF_HZ = 200.0f;
    constexpr float MAG_TUNE_REF_HZ = 25.0f;

    constexpr float acc_sigma_ref_mps2 = 0.12f;
    constexpr float gyr_sigma_ref_rps  = 0.00135f;
    constexpr float mag_sigma_ref_uT   = 0.80f;

    const float imu_rate_scale = sqrtf(IMU_TUNE_REF_HZ / LOOP_HZ);
    const float mag_rate_runtime_hz = 1000.0f / static_cast<float>(magPollMs_());
    const float mag_rate_scale = sqrtf(MAG_TUNE_REF_HZ / mag_rate_runtime_hz);

    const Vector3f sigma_a(acc_sigma_ref_mps2 * imu_rate_scale,
                           acc_sigma_ref_mps2 * imu_rate_scale,
                           acc_sigma_ref_mps2 * imu_rate_scale);
    const Vector3f sigma_g(gyr_sigma_ref_rps * imu_rate_scale,
                           gyr_sigma_ref_rps * imu_rate_scale,
                           gyr_sigma_ref_rps * imu_rate_scale);
    const float sigma_m_uT = mag_sigma_ref_uT * mag_rate_scale;
    const Vector3f sigma_m(sigma_m_uT, sigma_m_uT, sigma_m_uT);

    Fusion::Config fcfg;
    fcfg.with_mag = true;
    fcfg.online_tune_warmup_sec = ONLINE_TUNE_WARMUP_SEC;

    fcfg.freeze_acc_bias_until_live = true;
    fcfg.Racc_warmup_std = 0.18f;

    fcfg.sigma_a = sigma_a;
    fcfg.sigma_g = sigma_g;
    fcfg.sigma_m = sigma_m;

    fcfg.mag_delay_sec = 0.0f;
    fcfg.mag_init_min_mag_norm   = 5.0f;

    fcfg.enable_displacement_detrend = true;

    fusion_.begin(fcfg);

    auto& ff = fusion_.raw();
    ff.enableLinearBlock(WANT_LINEAR_BLOCK);
    ff.enableTuner(true);
    ff.setWithMag(true);
    ff.setSFactor(1.4f);
    ff.setTauCoeff(1.5f);
    ff.setSigmaCoeff(0.85f);
    ff.setRSCoeff(0.3f);
    ff.setRSXYFactor(0.23f);
    ff.setAccNoiseFloorSigma(ACC_NOISE_FLOOR_SIGMA_DEFAULT);
    ff.enableClamp(true);
    ff.setFreqInputCutoffHz(6.0f);

    rot_inited_   = false;
    rot_dpm_filt_ = 0.0f;
    gyro_bias_ok_ = false;
    gyro_bias_ema_.setZero();

    heading_deg_   = 0.0f;
    heading_valid_ = false;

    mag_gate_last_ms_ = 0;
    mag_ok_           = false;
    mag_fresh_        = false;
    mag_norm_uT_      = NAN;
    heading_mag_deg_  = NAN;
    heading_mag_ok_   = false;

    stale_frame_count_ = 0;
    have_last_sample_us_ = false;
    last_sample_us_      = 0;
    heave_raw_m_         = 0.0f;
    heave_baseline_m_    = 0.0f;
    heave_wave_raw_m_    = 0.0f;
    heave_wave_clean_m_  = 0.0f;
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
    if (!saved) {
      ui_.notSavedNotice();
    }

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
    const float t01 = 1.0f - static_cast<float>(remain) / static_cast<float>(M5UiCfg::MENU_TAP_WINDOW_MS);
    ui_.bar01(t01);
  }
#endif

  void updateFilter_(const ImuSample& s) {
    dt_ = computeFusionDtFromSampleTimestamp_(s);
    const float tempC = std::isfinite(s.tempC) ? s.tempC : 35.0f;

    const Vector3f a_raw = s.a;
    const Vector3f w_raw = s.w;

    a_raw_norm_ = a_raw.norm();

    a_cal_ = runtime_.applyAccel(a_raw, tempC);
    w_cal_ = runtime_.applyGyro(w_raw, tempC);
    m_cal_ = runtime_.applyMag(s.m);

    mag_norm_uT_ = m_cal_.norm();
    mag_ok_ = std::isfinite(mag_norm_uT_) && (mag_norm_uT_ > 5.0f) && (mag_norm_uT_ < 200.0f);
    mag_fresh_ = updateMagFreshGate_(mag_ok_, millis());

    const bool still =
        (fabsf(a_cal_.norm() - g_std) < ROT_STILL_G_TOL_FRAC * g_std) &&
        (w_cal_.norm() < ROT_STILL_GYRO_RAD_S);

    fusion_.update(dt_, w_cal_, a_cal_);
    if (mag_ok_ && mag_fresh_) {
      fusion_.updateMag(m_cal_);
    }

    heading_valid_ = fusion_.hasMagNorthLock();

    Eigen::Quaternionf q_bw = fusion_.raw().mekf().quaternion_boat();
    q_bw.normalize();

    float roll_est_deg = roll_deg_;
    float pitch_est_deg = pitch_deg_;
    float heading_est_deg = heading_deg_;

    if (rollPitchHeadingFromQuatBw_(q_bw, roll_est_deg, pitch_est_deg, heading_est_deg)) {
      roll_deg_  = roll_est_deg;
      pitch_deg_ = pitch_est_deg;
      heading_deg_ = heading_est_deg;
    }

    // Diagnostic-only magnetic heading using filter tilt + current mag sample.
    heading_mag_ok_ = false;
    heading_mag_deg_ = NAN;
    if (mag_ok_) {
      const Vector3f down_w(0.0f, 0.0f, 1.0f);
      Vector3f down_b = quatRotate_(q_bw.conjugate(), down_w);
      const float dn = down_b.norm();
      if (dn > 1e-6f) {
        down_b /= dn;
        float hdg_mag = NAN;
        if (magneticHeadingFromDownAndMagBody_(down_b, m_cal_, hdg_mag)) {
          heading_mag_ok_ = true;
          heading_mag_deg_ = hdg_mag;
        }
      }
    }

    if (still) {
      const float alpha_b = 1.0f - expf(-dt_ / ROT_BIAS_TAU_S);
      if (!gyro_bias_ok_) {
        gyro_bias_ok_  = true;
        gyro_bias_ema_ = w_cal_;
      } else {
        gyro_bias_ema_ += alpha_b * (w_cal_ - gyro_bias_ema_);
      }
    }

    Vector3f w_use = w_cal_;
    if (gyro_bias_ok_) w_use -= gyro_bias_ema_;

    const Vector3f w_world = quatRotate_(q_bw, w_use);
    float rot_dpm_meas = w_world.z() * RAD_TO_DEG * 60.0f;
    rot_dpm_meas = clampf_(rot_dpm_meas, -720.0f, 720.0f);

    const float tau_rot = 1.5f;
    const float alpha_r = 1.0f - expf(-dt_ / tau_rot);
    if (!rot_inited_) {
      rot_inited_ = true;
      rot_dpm_filt_ = rot_dpm_meas;
    } else {
      rot_dpm_filt_ += alpha_r * (rot_dpm_meas - rot_dpm_filt_);
    }

    const Vector3f displacement_raw_up_m = fusion_.displacementUpMeters();
    const auto& displacement_det_out = fusion_.displacementDetrend();

    heave_m_         = displacement_raw_up_m.z();
    wave_envelope_m_ = fusion_.raw().getDisplacementScale();
    wave_hz_         = fusion_.raw().getFreqHz();

    heave_raw_m_        = heave_m_;
    heave_baseline_m_   = displacement_det_out.baseline_slow.z();
    heave_wave_raw_m_   = displacement_det_out.wave_raw.z();
    heave_wave_clean_m_ = displacement_det_out.wave_clean.z();
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
    ui_.line("Set SEA_STATE_ENABLE_WIZARD=1");
#endif
    ui_.line("");
    ui_.line("Fusion: FULL");
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
    const bool tiltWarn = (fabsf(roll_deg_) > 35.0f) || (fabsf(pitch_deg_) > 35.0f);
    const float hdg_draw = heading_valid_ ? heading_deg_ : 0.0f;
    compass_ui_.draw(hdg_draw, heading_valid_ && mag_ok_, mag_norm_uT_, tiltWarn);
  }

  void updateUI_text_() {
    ui_.setReadRotation();
    ui_.title("COMPASS");

    if (heading_valid_) {
      M5.Display.printf("HDG: %6.1f deg\n", static_cast<double>(heading_deg_));
    } else {
      M5.Display.printf("HDG:   --- WAIT\n");
    }

    M5.Display.printf("HDM: %6.1f %s\n",
                      static_cast<double>(heading_mag_deg_),
                      heading_mag_ok_ ? "MAG" : "---");
    M5.Display.printf("ROL: %6.1f deg\n", static_cast<double>(roll_deg_));
    M5.Display.printf("PIT: %6.1f deg\n", static_cast<double>(pitch_deg_));
    M5.Display.printf("HEV: %6.3f m\n", static_cast<double>(heave_m_));
    M5.Display.printf("FRQ: %6.3f Hz\n", static_cast<double>(wave_hz_));
    M5.Display.printf("MAG: %s %s\n", mag_ok_ ? "OK " : "BAD", mag_fresh_ ? "NEW" : "OLD");
    M5.Display.printf("|m|: %6.1f uT\n", static_cast<double>(mag_norm_uT_));
    M5.Display.printf("|aR|:%5.2f |aC|:%5.2f\n",
                      static_cast<double>(a_raw_norm_),
                      static_cast<double>(a_cal_.norm()));
    ui_.line("");
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
    const bool valid = fusion_.isLive() && heading_valid_;
    if (heading_valid_) {
      nmea_hdm(SEA_STATE_NMEA_TALKER, heading_deg_);
    }
    nmea_xdr_pitch_roll(SEA_STATE_NMEA_TALKER, pitch_deg_, roll_deg_);
    nmea_xdr_heave(SEA_STATE_NMEA_TALKER, heave_wave_clean_m_);
    nmea_xdr_freq(SEA_STATE_NMEA_TALKER, wave_hz_);
    nmea_rot(SEA_STATE_NMEA_TALKER, rot_dpm_filt_, valid);
#else
  #if ARDUINO_PLOTTER
    Serial.printf(
      "HrawCm:%+.3f\tHwaveEnvelopeCm:%+.3f\tHwaveCleanCm:%+.3f\n",
      static_cast<double>(heave_raw_m_ * 100.0f),
      static_cast<double>(wave_envelope_m_ * 100.0f),
      static_cast<double>(heave_wave_clean_m_ * 100.0f));
  #else
    Serial.printf(
      "hdg_filt=%s%.2f | hdg_mag=%s%.2f | dt_imu_ms=%.2f | mag_ok=%u | mag_fresh=%u | |m|=%.2f"
      " | acc_ned[m/s^2] N=%.3f E=%.3f D=%.3f"
      " | gyro_ned[rad/s] N=%.3f E=%.3f D=%.3f"
      " | mag_ned[uT] N=%.2f E=%.2f D=%.2f\n",
      heading_valid_ ? "" : "~",
      static_cast<double>(heading_deg_),
      heading_mag_ok_ ? "" : "~",
      static_cast<double>(heading_mag_deg_),
      static_cast<double>(dt_ * 1000.0f),
      mag_ok_ ? 1U : 0U,
      mag_fresh_ ? 1U : 0U,
      static_cast<double>(mag_norm_uT_),
      static_cast<double>(a_cal_.x()),
      static_cast<double>(a_cal_.y()),
      static_cast<double>(a_cal_.z()),
      static_cast<double>(w_cal_.x()),
      static_cast<double>(w_cal_.y()),
      static_cast<double>(w_cal_.z()),
      static_cast<double>(m_cal_.x()),
      static_cast<double>(m_cal_.y()),
      static_cast<double>(m_cal_.z()));
  #endif
#endif
  }
};

static FusionApp g_app;

void setup() { g_app.begin(); }
void loop()  { g_app.tick(); }
