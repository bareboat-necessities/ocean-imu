/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R SeaStateFusionFilter_TFG (+ optional IMU Calibration Wizard)

    - right-invariant two-frame Lie-group INS
    - filter learns tilt during startup through the private Mahony proxy
    - one-shot magnetic north lock inside SeaStateFusionFilter_TFG
    - measured magnetic heading is available from the first valid IMU/mag sample
    - startup tilt comes from the private proxy, then from the live filter
    - startup compass does not wait for INS; fused filter yaw takes over once Live

  Assumptions:
    - fusion_.mekf().quaternion() returns BODY->WORLD rotation (q_bw), world frame is NED (+Z down).
    - runtime_.apply* already includes temperature compensation.

  HOW THIS DIFFERS FROM THE OU SKETCHES.

  SeaStateFusionFilter_TFG carries the core of the OU-III orchestrator -- the
  canonical period statistics, the three-layer adaptation, the startup
  tilt/magnetic acquisition, the continuous hard-iron estimator and the
  accelerometer vibration guard -- but not OU-III's wave-direction estimator,
  its displacement detrender or its alternative frequency trackers.  Those live
  in the OU wrapper class, and there is no TFG wrapper.

  So this sketch composes them here, from the same shared components the OU
  wrapper uses and driven the same way:

    - wave direction:  wave_direction::heading_frame_acceleration ->
                       VesselRaoEqualizer -> KalmanWaveDirection ->
                       WaveDirectionDetector
    - heave detrend:   AdaptiveWaveDetrender3D on the world-frame position,
                       with the filter's own period as external frequency
                       guidance

  Everything the application publishes therefore has the same meaning it has in
  the OU-III sketch, and the estimator underneath it is the only thing that
  changed.

  One knob the OU sketches expose is absent here: Kalman3D_Wave_TFG has no IMU
  lever-arm correction, so there is no IMU_LEVER_ARM_*_M block.  Mount the
  device as close to the vessel centre of gravity as the installation allows,
  or use one of the OU sketches where the offset is large enough to matter.
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
  #define SEA_STATE_SERIAL_NMEA 1
#endif

#ifndef SEA_STATE_NMEA_TALKER
  #define SEA_STATE_NMEA_TALKER "II"
#endif

constexpr float g_std      = atoms3r_ical::ImuCalCfg::g_std;
constexpr float FREQ_GUESS = 0.3f;

#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
#include "detrend/AdaptiveWaveDetrender3D.h"
#include "wave_dir/KalmanWaveDirection.h"
#include "wave_dir/VesselRaoEqualizer.h"
#include "wave_dir/WaveDirectionDetector.h"
#include "wave_dir/WaveDirectionFrame.h"

static constexpr float    LOOP_HZ        = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = static_cast<uint32_t>(1000000.0f / LOOP_HZ);

static constexpr uint32_t UI_REFRESH_MS   = 100;
static constexpr uint32_t DEBUG_SERIAL_MS = 100;
static constexpr uint32_t NMEA_SERIAL_MS  = 80;
static constexpr uint32_t NMEA_WAVE_MS    = 2000;
static constexpr uint32_t NMEA_STATUS_MS  = 5000;
static constexpr uint32_t NMEA_TEMP_MS    = 500;

static constexpr float ROT_BIAS_TAU_S       = 5.0f;
static constexpr float ROT_STILL_G_TOL_FRAC = 0.12f;
static constexpr float ROT_STILL_GYRO_RAD_S = 0.15f;

static constexpr float ONLINE_TUNE_WARMUP_SEC = 10.0f;

// Front-end accelerometer vibration guard corner, in Hz.
//
// Zero removes the guard and restores the unconditioned measurement path.  The
// default sits in the gap between the wave band and the lowest crank order a
// small auxiliary diesel puts on the hull; see docs/engine-noise-degradation.md
// for what it buys and what the group delay costs.
static constexpr float ACC_VIBRATION_GUARD_HZ =
    ocean_imu::tfg::ACC_VIBRATION_GUARD_HZ_DEFAULT;

// The loop-task stack size is the other compile-time knob this sketch sets.
// It lives at the bottom of the file, next to setup(); see the comment there
// for why it is not here.

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

// Generic mapper.
// If your WaveDirection enum does not use negative/positive numeric polarity,
// replace only this function.
static inline int waveDirectionSignPolarity_(WaveDirection s) {
  if (s == UNCERTAIN) {
    return 0;
  }

  const int raw = static_cast<int>(s);

  if (raw < 0) return -1;  // opposite axis, add 180 deg
  if (raw > 0) return +1;  // same axis

  return 0;
}

static inline int waveDirectionSignRaw_(WaveDirection s) {
  return static_cast<int>(s);
}

static inline bool signedWaveDirectionDeg_(float axis_deg,
                                           int sign_polarity,
                                           float& signed_deg_out)
{
  if (!std::isfinite(axis_deg) || sign_polarity == 0) {
    signed_deg_out = NAN;
    return false;
  }

  signed_deg_out = wrap360_(axis_deg + (sign_polarity < 0 ? 180.0f : 0.0f));
  return true;
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
  heading_deg_out = NAN;
  if (!down_b_unit.allFinite() || !mag_b_uT.allFinite()) return false;
  const Vector3f FWD_B(1.0f, 0.0f, 0.0f);

  Vector3f d = down_b_unit;
  const float dn = d.norm();
  if (!std::isfinite(dn) || !(dn > 1e-6f)) return false;
  d /= dn;

  Vector3f m = mag_b_uT;
  const float mn = m.norm();
  if (!std::isfinite(mn) || !(mn > 1e-6f)) return false;
  m /= mn;

  Vector3f east_b = d.cross(m);
  const float en = east_b.norm();
  if (!std::isfinite(en) || !(en > 1e-6f)) return false;
  east_b /= en;

  Vector3f north_b = east_b.cross(d);
  const float nn = north_b.norm();
  if (!std::isfinite(nn) || !(nn > 1e-6f)) return false;
  north_b /= nn;

  const float e = east_b.dot(FWD_B);
  const float n = north_b.dot(FWD_B);
  if (!(std::hypot(e, n) > 1e-6f)) return false;
  heading_deg_out = wrap360_(atan2f(e, n) * RAD_TO_DEG);
  return std::isfinite(heading_deg_out);
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
  using UI     = atoms3r_ical::M5Ui;
  using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;

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

  uint32_t last_ui_ms_          = 0;
  uint32_t last_serial_ms_      = 0;
  uint32_t last_wave_nmea_ms_   = 0;
  uint32_t last_status_nmea_ms_ = 0;
  uint32_t last_temp_nmea_ms_   = 0;

  Fusion fusion_{};

  // Application-level companions to the filter.  The OU wrapper owns these;
  // TFG has no wrapper, so they live here and are driven from updateFilter_().
  AdaptiveWaveDetrender3D         displacement_detrender_{};
  AdaptiveWaveDetrender3D::Output displacement_det_out_{};
  Vector3f                        displacement_up_m_ = Vector3f::Zero();

  wave_direction::VesselRaoEqualizer direction_rao_{};
  KalmanWaveDirection                dir_filter_{2.0f * static_cast<float>(M_PI) * FREQ_GUESS};
  WaveDirectionDetector<float>       dir_sign_{0.002f, 0.005f};

  uint32_t mag_gate_last_ms_ = 0;

  Vector3f a_cal_ = Vector3f::Zero();
  Vector3f w_cal_ = Vector3f::Zero();
  Vector3f m_cal_ = Vector3f::Constant(NAN);
  float    a_raw_norm_ = 0.0f;

  float dt_               = 0.0f;
  float roll_deg_         = 0.0f;
  float pitch_deg_        = 0.0f;
  float heading_deg_      = NAN;
  bool  heading_valid_    = false;
  bool  heading_fused_    = false;
  float heave_m_          = 0.0f;
  float heave_speed_mps_  = 0.0f;
  float wave_envelope_m_  = 0.0f;
  float wave_hz_          = FREQ_GUESS;

  float wave_axis_deg_       = NAN;
  bool  wave_axis_ok_        = false;
  WaveDirection wave_sign_   = UNCERTAIN;
  int   wave_sign_raw_       = 0;
  int   wave_sign_polarity_  = 0;
  bool  wave_sign_ok_        = false;
  float wave_dir_deg_        = NAN;
  bool  wave_dir_ok_         = false;
  float wave_dir_conf_pct_   = 0.0f;

  float heave_raw_m_        = 0.0f;
  float heave_baseline_m_   = 0.0f;
  float heave_wave_raw_m_   = 0.0f;
  float heave_wave_clean_m_ = 0.0f;

  bool  mag_ok_          = false;
  bool  mag_fresh_       = false;
  float mag_norm_uT_     = NAN;
  float heading_mag_deg_ = NAN;
  float imu_temp_c_      = NAN;
  bool  heading_mag_ok_  = false;

  uint16_t stale_frame_count_ = 0;
  uint32_t last_skipped_total_ = 0;
  bool     have_last_sample_us_ = false;
  uint32_t last_sample_us_      = 0;

  bool     rot_inited_    = false;
  float    rot_dpm_filt_  = 0.0f;
  bool     gyro_bias_ok_       = false;
  bool     gyro_bias_learning_ = false;
  bool     acc_bias_estimating_ = false;
  Vector3f gyro_bias_ema_      = Vector3f::Zero();

  uint32_t rate_window_ms_   = 0;
  uint32_t imu_sample_count_ = 0;
  uint32_t mag_sample_count_ = 0;
  float    imu_sample_hz_    = 0.0f;
  float    mag_sample_hz_    = 0.0f;

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

  // The filter reports the wave period, not a frequency.  Everything
  // downstream of it -- the direction demodulator's carrier, the detrender's
  // external guidance, the plotter line -- wants Hz, so convert once here and
  // fall back to the startup prior while the statistic is not usable yet.
  float waveFrequencyHz_() const {
    const float T = fusion_.getWavePeriodSec();
    if (!fusion_.wavePeriodUsable() || !std::isfinite(T) || !(T > 1.0e-3f)) {
      return FREQ_GUESS;
    }
    return 1.0f / T;
  }

  // Envelope of the displacement the current operating point implies, i.e. the
  // same quantity OU-III publishes as getDisplacementScale().  Both filters
  // drive the identical OU chain, so the scalar reduction
  //
  //     Hs/2 = C_HS * sigma_aw * tau^2 / 2,   C_HS = 2*sqrt(2)/pi^2
  //
  // applies unchanged here; only the tuned (sigma_aw, tau) differ.
  float displacementScaleM_() const {
    const float tau   = fusion_.getTauApplied();
    const float sigma = fusion_.getSigmaApplied();
    if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
    constexpr float C_HS = 2.0f * 1.41421356237f / (3.14159265359f * 3.14159265359f);
    return C_HS * sigma * tau * tau / 2.0f;
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
    rate_window_ms_ = millis();
    imu_sample_count_ = 0;
    mag_sample_count_ = 0;
    imu_sample_hz_ = 0.0f;
    mag_sample_hz_ = 0.0f;
    have_last_sample_us_ = false;
    last_sample_us_      = 0;
  }

  void resetFusion_() {
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
    const float sigma_m_uT = mag_sigma_ref_uT * mag_rate_scale;
    const Vector3f sigma_m(sigma_m_uT, sigma_m_uT, sigma_m_uT);

    Fusion::Config fcfg;
    fcfg.with_mag = true;
    fcfg.online_tune_warmup_sec = ONLINE_TUNE_WARMUP_SEC;

    fcfg.sigma_a = sigma_a;
    fcfg.sigma_m = sigma_m;

    // TFG takes the gyroscope as a noise density rather than a per-axis sigma,
    // so the rate scaling that the accelerometer and magnetometer get is
    // already carried by the density itself and must not be applied twice.
    fcfg.gyro_noise_density = gyr_sigma_ref_rps;

    fcfg.mag_delay_sec = 0.0f;
    fcfg.mag_init_min_mag_norm = 5.0f;

    fcfg.gravity_magnitude = g_std;

    fcfg.acc_vibration_guard_hz = ACC_VIBRATION_GUARD_HZ;

    fusion_.begin(fcfg);

    fusion_.enableTuner(true);
    // Use frequent virtual constraints instead of 4--8 Hz heave corrections.
    // Both existing r_S laws account for the selected correction cadence.
    // No output smoothing, clipping or stationarity lock is applied.
    fusion_.setTauScaledPseudoCadence(false);
    // The tuner coefficients (S_factor, tau, sigma, r_S, r_S XY) are
    // deliberately not overridden here.  The header defaults are the operating
    // point the committed TFG study was run at; overriding them puts the device
    // at an operating point no evidence covers.  Override only alongside a
    // study re-run.

    displacement_detrender_.setConfig(
        seastate::common::defaultDisplacementDetrenderConfig<
            AdaptiveWaveDetrender3D::Config>(FREQ_GUESS));
    displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
    displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
    displacement_up_m_.setZero();

    direction_rao_.reset();
    // Reassigned rather than reset(): reset() leaves the carrier frequency
    // wherever the previous run left it, and a fresh start belongs on the
    // startup prior.  This is what the OU wrapper does on its own reset.
    dir_filter_ = KalmanWaveDirection(2.0f * static_cast<float>(M_PI) * FREQ_GUESS);
    dir_sign_.reset();

    rot_inited_   = false;
    rot_dpm_filt_ = 0.0f;
    gyro_bias_ok_ = false;
    gyro_bias_learning_ = false;
    acc_bias_estimating_ = false;
    gyro_bias_ema_.setZero();

    heading_deg_   = NAN;
    heading_valid_ = false;
    heading_fused_ = false;

    mag_gate_last_ms_ = 0;
    mag_ok_           = false;
    mag_fresh_        = false;
    mag_norm_uT_      = NAN;
    heading_mag_deg_  = NAN;
    imu_temp_c_       = NAN;
    heading_mag_ok_   = false;

    stale_frame_count_ = 0;
    have_last_sample_us_ = false;
    last_sample_us_      = 0;
    heave_speed_mps_     = 0.0f;
    heave_raw_m_         = 0.0f;
    heave_baseline_m_    = 0.0f;
    heave_wave_raw_m_    = 0.0f;
    heave_wave_clean_m_  = 0.0f;
    wave_envelope_m_     = 0.0f;
    wave_hz_             = FREQ_GUESS;

    wave_axis_deg_      = NAN;
    wave_axis_ok_       = false;
    wave_sign_          = UNCERTAIN;
    wave_sign_raw_      = 0;
    wave_sign_polarity_ = 0;
    wave_sign_ok_       = false;
    wave_dir_deg_       = NAN;
    wave_dir_ok_        = false;
    wave_dir_conf_pct_  = 0.0f;
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

  // Direction branch, driven exactly as the OU-III wrapper drives it.
  //
  // The demodulator needs inertial acceleration in a levelled frame whose
  // horizontal axes stay aligned with the bow, so the body sensor bias is
  // removed before levelling, as the MEKF's own acceleration model does.
  //
  // Nothing is fed while the filter is Cold: the attitude that does the
  // levelling is still the startup proxy's, and a direction estimate built on
  // it would only have to be thrown away.  The published outputs are gated on
  // isLive() below in any case.
  void updateWaveDirection_(const Eigen::Quaternionf& q_bw, float tempC, float dt) {
    if (fusion_.stage() == Fusion::StartupStage::Cold) return;

    const Vector3f corrected = a_cal_ - fusion_.mekf().get_acc_bias_at_temperature(tempC);
    const auto direction_input =
        wave_direction::heading_frame_acceleration<float>(q_bw, corrected, g_std);
    if (!direction_input.heading_valid) return;

    const auto direction_matched = direction_rao_.step(direction_input, dt);

    const float omega = 2.0f * static_cast<float>(M_PI) * wave_hz_;
    dir_filter_.update(direction_matched.forward_ms2,
                       direction_matched.starboard_ms2,
                       omega, dt);

    const Eigen::Vector2f propagation_axis_boat = dir_filter_.getAxis();
    wave_sign_ = dir_sign_.update(
        direction_matched.forward_ms2,
        direction_matched.starboard_ms2,
        direction_matched.up_ms2,
        propagation_axis_boat.x(), propagation_axis_boat.y(),
        dt, dir_filter_.getLastStableConfidence());
  }

  // Heave branch, driven exactly as the OU-III wrapper drives it: the
  // world-frame position with Z flipped to "up", then an adaptive detrender
  // that is told the filter's own wave frequency whenever that frequency is
  // inside the detrender's configured band and the filter is live.
  void updateDisplacement_(float dt) {
    const Vector3f pos_ned_m = fusion_.mekf().get_position();
    displacement_up_m_ = Vector3f(pos_ned_m.x(), pos_ned_m.y(), -pos_ned_m.z());

    const bool ext_freq_valid =
        fusion_.isLive() &&
        std::isfinite(wave_hz_) &&
        (wave_hz_ >= displacement_detrender_.config().min_wave_freq_hz) &&
        (wave_hz_ <= displacement_detrender_.config().max_wave_freq_hz);

    displacement_det_out_ =
        displacement_detrender_.update(displacement_up_m_, dt, wave_hz_, ext_freq_valid);
  }

  // Before INS readiness, report the measured compass using startup tilt.
  // Once the magnetically informed filter is Live, publish its fused yaw.
  void updateCompassHeading_(const Eigen::Quaternionf& q_bw, bool attitude_ok,
                               bool fused_ready = false, float fused_heading_deg = NAN) {
    heading_mag_ok_ = false;
    heading_mag_deg_ = NAN;
    if (mag_ok_ && attitude_ok) {
      const Vector3f down_b = quatRotate_(q_bw.conjugate(), Vector3f::UnitZ());
      heading_mag_ok_ = magneticHeadingFromDownAndMagBody_(
          down_b, m_cal_, heading_mag_deg_);
    }
    heading_fused_ = fused_ready && attitude_ok && std::isfinite(fused_heading_deg);
    heading_valid_ = heading_fused_ || heading_mag_ok_;
    heading_deg_ = heading_fused_ ? wrap360_(fused_heading_deg)
        : (heading_mag_ok_ ? heading_mag_deg_ : NAN);
  }

  void updateFilter_(const ImuSample& s) {
    dt_ = computeFusionDtFromSampleTimestamp_(s);
    const float tempC = std::isfinite(s.tempC) ? s.tempC : 35.0f;
    imu_temp_c_ = tempC;
    ++imu_sample_count_;

    const Vector3f a_raw = s.a;
    const Vector3f w_raw = s.w;

    a_raw_norm_ = a_raw.norm();

    a_cal_ = runtime_.applyAccel(a_raw, tempC);
    w_cal_ = runtime_.applyGyro(w_raw, tempC);
    m_cal_ = runtime_.applyMag(s.m);

    mag_norm_uT_ = m_cal_.norm();
    mag_ok_ = std::isfinite(mag_norm_uT_) && (mag_norm_uT_ > 5.0f) && (mag_norm_uT_ < 200.0f);
    mag_fresh_ = updateMagFreshGate_(mag_ok_, millis());
    if (mag_fresh_) ++mag_sample_count_;

    const bool still =
        (fabsf(a_cal_.norm() - g_std) < ROT_STILL_G_TOL_FRAC * g_std) &&
        (w_cal_.norm() < ROT_STILL_GYRO_RAD_S);

    fusion_.update(dt_, w_cal_, a_cal_, tempC);
    if (mag_ok_ && mag_fresh_) {
      fusion_.updateMag(m_cal_);
    }

    // During startup the MEKF is deliberately held. Publish the Mahony proxy
    // attitude instead; switch atomically to the fused attitude at Live.
    Eigen::Quaternionf q_bw;
    bool have_attitude = false;
    if (!fusion_.isLive() &&
        fusion_.startupInitPolicy() == Fusion::StartupInitPolicy::MahonyProxy) {
      have_attitude = fusion_.startupTiltQuaternion(q_bw);
    } else {
      q_bw = fusion_.mekf().quaternion();
      have_attitude = q_bw.coeffs().allFinite();
    }

    float roll_est_deg = roll_deg_;
    float pitch_est_deg = pitch_deg_;
    float heading_est_deg = heading_deg_;
    const bool attitude_ok = have_attitude &&
        rollPitchHeadingFromQuatBw_(q_bw, roll_est_deg, pitch_est_deg, heading_est_deg);
    if (attitude_ok) {
      q_bw.normalize();
      roll_deg_  = roll_est_deg;
      pitch_deg_ = pitch_est_deg;
    }

    updateCompassHeading_(q_bw, attitude_ok,
        fusion_.isLive() && fusion_.mekf().has_magnetic_reference(), heading_est_deg);

    gyro_bias_learning_ = still;

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

    const float tau_rot = 0.1f;
    const float alpha_r = 1.0f - expf(-dt_ / tau_rot);
    if (!rot_inited_) {
      rot_inited_ = true;
      rot_dpm_filt_ = rot_dpm_meas;
    } else {
      rot_dpm_filt_ += alpha_r * (rot_dpm_meas - rot_dpm_filt_);
    }

    const Vector3f velocity_ned_mps = fusion_.mekf().get_velocity();
    heave_speed_mps_ = -velocity_ned_mps.z();

    wave_hz_         = waveFrequencyHz_();
    wave_envelope_m_ = displacementScaleM_();

    updateWaveDirection_(q_bw, tempC, dt_);
    updateDisplacement_(dt_);

    // Raw integrated position can carry a large DC/random-walk component.
    // The device HEV channel is wave heave, matching the NMEA output and OU
    // sketches; retain raw position only in the diagnostic heave_raw_m_ field.
    heave_m_ = displacement_det_out_.wave_clean.z();

    wave_axis_deg_ = dir_filter_.getAxisDegrees();
    wave_axis_ok_  = fusion_.isLive() && std::isfinite(wave_axis_deg_);

    wave_sign_raw_      = waveDirectionSignRaw_(wave_sign_);
    wave_sign_polarity_ = waveDirectionSignPolarity_(wave_sign_);
    wave_sign_ok_       = fusion_.isLive() && (wave_sign_ != UNCERTAIN);

    wave_dir_ok_ = signedWaveDirectionDeg_(
        wave_axis_deg_,
        wave_sign_polarity_,
        wave_dir_deg_);

    wave_dir_ok_ = wave_dir_ok_ && fusion_.isLive();
    wave_dir_conf_pct_ = wave_dir_ok_ ? 100.0f : (wave_axis_ok_ ? 50.0f : 0.0f);
    acc_bias_estimating_ = fusion_.isLive();

    const uint32_t now_ms = millis();
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

    heave_raw_m_        = displacement_up_m_.z();
    heave_baseline_m_   = displacement_det_out_.baseline_slow.z();
    heave_wave_raw_m_   = displacement_det_out_.wave_raw.z();
    heave_wave_clean_m_ = displacement_det_out_.wave_clean.z();
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
    ui_.line("Fusion: TFG");
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
    M5.Display.printf("WAV: %6.1f s=%d\n",
                      static_cast<double>(wave_dir_ok_ ? wave_dir_deg_ : wave_axis_deg_),
                      wave_sign_raw_);
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
    // A startup magnetic compass is usable HDM, but is not a ready INS.
    const bool valid = fusion_.isLive() && heading_fused_;
    if (heading_valid_) {
      nmea_hdm(SEA_STATE_NMEA_TALKER, heading_deg_);
    }
    nmea_xdr_pitch_roll(SEA_STATE_NMEA_TALKER, pitch_deg_, roll_deg_); // startup proxy attitude until Live
    nmea_xdr_heave(SEA_STATE_NMEA_TALKER, heave_wave_clean_m_);
    nmea_xdr_heave_speed(SEA_STATE_NMEA_TALKER, heave_speed_mps_);

    if (now_ms - last_wave_nmea_ms_ >= NMEA_WAVE_MS) {
      last_wave_nmea_ms_ = now_ms;
      nmea_xdr_wave_axis_rel(SEA_STATE_NMEA_TALKER, wave_axis_deg_, wave_axis_ok_);
      nmea_xdr_wave_direction_rel(SEA_STATE_NMEA_TALKER, wave_dir_deg_, wave_dir_ok_);
      nmea_txt_wave_direction_sign(SEA_STATE_NMEA_TALKER, wave_sign_raw_, wave_sign_polarity_, wave_sign_ok_);
      nmea_txt_wave_direction_confidence(SEA_STATE_NMEA_TALKER, wave_dir_conf_pct_, fusion_.isLive());
    }

    if (now_ms - last_temp_nmea_ms_ >= NMEA_TEMP_MS) {
      last_temp_nmea_ms_ = now_ms;
      nmea_xdr_imu_temp(SEA_STATE_NMEA_TALKER, imu_temp_c_);
    }

    if (now_ms - last_status_nmea_ms_ >= NMEA_STATUS_MS) {
      last_status_nmea_ms_ = now_ms;
      nmea_txt_ins_status(SEA_STATE_NMEA_TALKER,
                          valid,
                          fusion_.isLive() && std::isfinite(heave_wave_clean_m_),
                          fusion_.magReferenceLearned(),
                          gyro_bias_learning_,
                          acc_bias_estimating_,
                          imu_sample_hz_,
                          mag_sample_hz_);
    }

    //nmea_xdr_freq(SEA_STATE_NMEA_TALKER, wave_hz_);
    nmea_rot(SEA_STATE_NMEA_TALKER, rot_dpm_filt_, valid);
#else
  #if ARDUINO_PLOTTER
    Serial.printf(
      "HrawCm:%+.3f\tHwaveEnvelopeCm:%+.3f\tHwaveCleanCm:%+.3f\tWavAxisDeg:%+.2f\tWavDirDeg:%+.2f\tWavSign:%d\n",
      static_cast<double>(heave_raw_m_ * 100.0f),
      static_cast<double>(wave_envelope_m_ * 100.0f),
      static_cast<double>(heave_wave_clean_m_ * 100.0f),
      static_cast<double>(wave_axis_deg_),
      static_cast<double>(wave_dir_deg_),
      wave_sign_raw_);
  #else
    Serial.printf(
      "hdg=%s%.2f | hdg_mag=%s%.2f | wav_axis=%.2f | wav_dir=%s%.2f | wav_sign=%d pol=%d"
      " | dt_imu_ms=%.2f | mag_ok=%u | mag_fresh=%u | |m|=%.2f"
      " | acc_ned[m/s^2] N=%.3f E=%.3f D=%.3f"
      " | gyro_ned[rad/s] N=%.3f E=%.3f D=%.3f"
      " | mag_ned[uT] N=%.2f E=%.2f D=%.2f\n",
      heading_valid_ ? "" : "~",
      static_cast<double>(heading_deg_),
      heading_mag_ok_ ? "" : "~",
      static_cast<double>(heading_mag_deg_),
      static_cast<double>(wave_axis_deg_),
      wave_dir_ok_ ? "" : "~",
      static_cast<double>(wave_dir_deg_),
      wave_sign_raw_,
      wave_sign_polarity_,
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

// ---------------------------------------------------------------------------
// Loop-task stack.
//
// This is not a safety margin, it is a requirement.  Kalman3D_Wave_TFG runs a
// 21-state covariance through Phi P Phi^T and a Joseph update, and Eigen
// builds each of those triple products through full NX x NX temporaries --
// 1764 bytes apiece at NX = 21.  The named locals that used to sit alongside
// them are now member scratch buffers (see src/kalman_tfg/Kalman3D_Wave_TFG.h),
// which took one live fusion step from ~27 kB of stack down to ~19 kB, but the
// expression temporaries remain: removing those would change Eigen's GEMM path
// and with it the filter's arithmetic, which is not something a stack budget
// gets to decide.
//
// Measured peak, painted-stack high-water mark over a 400 s host run of the
// whole sketch-level step (filter + detrender + direction chain):
//
//                       -O2        -Os
//     before          27.2 kB    27.2 kB
//     after           18.9 kB    15.3 kB
//
// The Arduino-ESP32 loop task gets 8 kB by default, so this still has to be
// raised.  32 kB leaves ~13 kB over the worst measurement for the sketch's own
// frames, the UI and the NMEA formatting, and for the xtensa/-funroll-loops
// build differing from the host one.
//
// On the device this is affordable: the build-MCU job reports 40,540 bytes of
// globals for this sketch and 287,140 bytes left for locals, against 30,828
// and 296,852 for the OU-III sketch.  The extra globals are the filter's
// scratch pool; the loop stack comes out of what is left.
//
// Lower it only against a fresh measurement on the device.
//
// WHY THIS IS AT THE BOTTOM OF THE FILE.  SET_LOOP_TASK_STACK_SIZE expands to
// a definition of getArduinoLoopTaskStackSize(), so it is a function
// definition like any other.  The Arduino builder inserts its generated
// prototypes immediately before the FIRST function definition in the sketch,
// and those prototypes name Vector3f -- so with the macro up among the
// configuration constants, the prototypes landed above the
// `using Vector3f = Eigen::Vector3f;` alias and the sketch failed to compile
// with "'Vector3f' does not name a type".  Down here the first function is
// clampf_(), exactly as in the OU sketches, and the prototypes land where
// they do there.
// ---------------------------------------------------------------------------
SET_LOOP_TASK_STACK_SIZE(32 * 1024);

static FusionApp g_app;

void setup() { g_app.begin(); }
void loop()  { g_app.tick(); }
