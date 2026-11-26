/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R SeaStateFusion_CARMA3 (+ optional IMU Calibration Wizard)

  EXPERIMENTAL research sketch.  It runs the three-axis CARMA filter of
  docs/carma-three-axis.md in the configuration that document recommends as
  the research candidate: CARMA(4,1) on X/Y and on Z, no S states (B).  It is
  not a shipping estimator (simulation evidence only; see section 12 of that
  document for the remaining obstacles); atomS3R_ins_kalman_ou3 is.

    - filter learns tilt during startup
    - one-shot magnetic north lock inside SeaStateFusion_CARMA3
    - measured magnetic heading is available from the first valid IMU/mag sample
    - startup tilt comes from the private proxy, then from the live filter
    - startup compass does not wait for INS; fused filter yaw takes over once Live

  Assumptions:
    - fusion_.raw().mekf().quaternion_boat() returns BODY->WORLD quaternion (q_bw), world frame is NED (+Z down).
    - runtime_.apply* already includes temperature compensation.
*/

#include <Arduino.h>

#ifndef SEA_STATE_ENABLE_WIZARD
  #define SEA_STATE_ENABLE_WIZARD 1
#endif

#ifndef ARDUINO_PLOTTER
  #define ARDUINO_PLOTTER 1
#endif

#include <M5Unified.h>
#include <cmath>
#include <algorithm>
#include <cstring>

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
  #define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif
#include <ArduinoOceanImu.h>
#include "AtomS3R/AtomS3R_FilterOptimization.h"
ATOMS3R_BEGIN_FILTER_MATH
#include <ArduinoEigenDense.h>
ATOMS3R_END_FILTER_MATH

#include "AtomS3R/AtomS3R_CalibrationGate.h"
#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_ImuSpikeFilter.h"
#include "AtomS3R/AtomS3R_M5Ui.h"
#if SEA_STATE_ENABLE_WIZARD
  #include "AtomS3R/ImuCalWizardRunner.h"
#endif
#include "AtomS3R/AtomS3R_MarineCompassDisplay.h"
#include "AtomS3R/AtomS3R_OutputTiming.h"
#include "AtomS3R/AtomS3R_MarineNmeaBatch.h"

#ifndef SEA_STATE_UI_DEFAULT_GRAPHICS
  #define SEA_STATE_UI_DEFAULT_GRAPHICS 1
#endif

#ifndef SEA_STATE_SERIAL_NMEA
  #define SEA_STATE_SERIAL_NMEA 1
#endif

#ifndef SEA_STATE_NMEA_TALKER
  #define SEA_STATE_NMEA_TALKER "II"
#endif

// Nominal-g unit (also the generic default the filter headers declare as
// extern g_std). Physical gravity is g_local: the calibration target, and the
// gravity the fusion filter removes (Fusion::Config::gravity_magnitude).
constexpr float g_std      = atoms3r_ical::ImuCalCfg::g_std;
constexpr float g_local    = atoms3r_ical::ImuCalCfg::g_cal_local;
constexpr float FREQ_GUESS = 0.3f;

#define ZERO_CROSSINGS_SCALE          1.0f
#define ZERO_CROSSINGS_DEBOUNCE_TIME  0.12f
#define ZERO_CROSSINGS_STEEPNESS_TIME 0.21f

ATOMS3R_BEGIN_FILTER_MATH
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
ATOMS3R_END_FILTER_MATH
// The OU-III sketch's startup settings (ONLINE_TUNE_WARMUP_SEC,
// ACC_NOISE_FLOOR_SIGMA_DEFAULT), which this sketch carries over unchanged.
#include "kalman_common/SeaStateFamilyDefaults.h"
#include "util/AngleUtils.h"
#include "AtomS3R/AtomS3R_SamplePacer.h"
#include "AtomS3R/AtomS3R_LoopTiming.h"
#include "util/ImuLoopHelpers.h"
#include "util/MagneticHeading.h"
#include "util/QuaternionUtils.h"
#include "wave_dir/WaveDirectionReport.h"

static constexpr float    LOOP_HZ        = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = static_cast<uint32_t>(1000000.0f / LOOP_HZ);

static constexpr float ROT_BIAS_TAU_S       = 5.0f;
static constexpr float ROT_STILL_G_TOL_FRAC = 0.12f;
static constexpr float ROT_STILL_GYRO_RAD_S = 0.15f;

// ---------------------------------------------------------------------------
// IMU installation lever arm (IMU position relative to the vessel centre of
// gravity), in BODY-NED metres:  x = forward (towards the bow),
//                                y = starboard,
//                                z = down.
//
// A rigidly mounted IMU that does not sit exactly on the CoG does not measure
// the CoG motion.  For a body-fixed offset r it measures
//
//     a_imu = a_cog + (alpha x r) + omega x (omega x r)
//
// i.e. the CoG specific force plus a tangential term (alpha = angular
// acceleration) and a centripetal term (omega = angular rate).  Both terms are
// deterministic and both are correlated with attitude, so an IMU bolted a few
// decimetres away from the CoG leaks roll/pitch motion into the same
// acceleration channel this filter integrates into heave, and biases the
// attitude estimate that rides on it.  docs/carma-imu-lever-arm.md measures
// the cost per direction and distance, and what the correction recovers.
//
// The filter can subtract that term itself, but only if it is told r.  It is
// left at ZERO here, which is also the filter's own default: a zero lever arm
// means "IMU is at the CoG", the correction is switched off completely, and
// the device behaves exactly as it did before this knob was exposed.  The call
// is written out explicitly in resetFusion_() so the adjustment point is
// visible rather than implied.
//
// To adjust for a real installation, measure the offset on the boat as
// (IMU position - CoG position) along the three axes above and enter it here.
// Example: an IMU 1.20 m forward of, 0.15 m to starboard of and 0.30 m ABOVE
// the CoG (z is positive DOWN, so "above" is negative):
//
//     static constexpr float IMU_LEVER_ARM_X_M =  1.20f;
//     static constexpr float IMU_LEVER_ARM_Y_M =  0.15f;
//     static constexpr float IMU_LEVER_ARM_Z_M = -0.30f;
//
// A few centimetres of accuracy is plenty; a wrong sign is worse than leaving
// the correction off, because it doubles the term instead of removing it.
// ---------------------------------------------------------------------------
static constexpr float IMU_LEVER_ARM_X_M = 0.0f;  // + forward (towards bow)
static constexpr float IMU_LEVER_ARM_Y_M = 0.0f;  // + starboard
static constexpr float IMU_LEVER_ARM_Z_M = 0.0f;  // + down

// ---------------------------------------------------------------------------
// Accelerometer noise model of the CARMA wave-excitation evidence
// (docs/carma-wave-excitation.md): the white noise per sample at LOOP_HZ and
// the bias random walk of the calibrated vertical acceleration, and the
// relative uncertainty of that noise variance.  The values are the simulated
// BMI270-class sensor the study used.  Replace them by the installed device's
// own values from a stationary Allan-variance record; 0.5 tolerates about a
// factor of two in the noise variance without reading sensor noise as waves.
// ---------------------------------------------------------------------------
static constexpr float CARMA_ACC_WHITE_SIGMA_MPS2      = 0.0148f;
static constexpr float CARMA_ACC_BIAS_RW_MPS2_SQRTS    = 5.0e-4f;
static constexpr float CARMA_ACC_NOISE_REL_UNCERTAINTY = 0.5f;

using namespace atoms3r_ical;
using Vector3f = Eigen::Vector3f;

namespace ins = ocean_imu::ins;
using wave_direction::WaveDirectionReport;

class FusionApp {
public:
  FusionApp() = default;

  void begin() {
    delay(200);            // allow power to stabilize
    Serial.setTxBufferSize(4096);
    Serial.begin(115200);
    Serial.setTxTimeoutMs(0);
    delay(100);

    auto cfg = M5.config();
    cfg.internal_imu = true;
    M5.begin(cfg);
    clearM5UnifiedImuCalibration();
    configureAtomS3RMag(Serial);
    delay(250);

    ui_.begin();
    display_.begin(SEA_STATE_UI_DEFAULT_GRAPHICS != 0);

    reloadBlobAndRuntime_();

#if SEA_STATE_ENABLE_WIZARD
    if (!have_blob_) {
      Serial.println("[BOOT] No saved calibration. Starting wizard...");
      const bool saved = runWizardFlow_(true);
      if (saved) {
        Serial.println("[BOOT] Wizard saved calibration.");
      } else {
        Serial.println("[BOOT] Wizard did not save calibration.");
      }
    }
#endif

    requireUsableCalibration(ui_, store_, blob_, have_blob_, runtime_, kSeaStateCalibrationWizard);
    reinitImu_();
    resetFusion_();
    drawHomeStatic_();

    delay(100);
  }

  void tick() {
    pacer_.wait();
    timing_.begin();

#if SEA_STATE_ENABLE_WIZARD
    Input::pollButton();

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
    uint32_t poll_ms;
    const bool got_sample = navigation_imu_.read(M5.Imu, sample_us, sample, temperature_, mag_gate_, poll_ms);

    if (got_sample) {
      stale_frame_count_ = 0;
      updateFilter_(sample, poll_ms);
    } else {
      if (stale_frame_count_ < 0xFFFFu) stale_frame_count_++;
      if ((stale_frame_count_ % 50u) == 0u) {
        Serial.printf("[IMU] waiting: stale=%u\n", static_cast<unsigned>(stale_frame_count_));
      }
    }

    streamSerial_();
    updateUI_();
    timing_.finish(sample_us, got_sample, LOOP_PERIOD_US);
  }

private:
  using UI = atoms3r_ical::M5Ui;
  // CARMA(4,1) X/Y without S_x, S_y; CARMA(4,1) Z without S_z.
  using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false>;

  UI ui_{};
  MarineCompassDisplay display_{ui_};

  ImuCalStoreNvs store_{};
  bool           have_blob_ = false;
  ImuCalBlobV5   blob_{};
  RuntimeCals    runtime_{};

#if SEA_STATE_ENABLE_WIZARD
  int      tap_count_ = 0;
  uint32_t tap_deadline_ms_ = 0;
#endif

  uint32_t last_ui_ms_          = 0;
  uint32_t last_serial_ms_      = 0;
#if SEA_STATE_SERIAL_NMEA
  atoms3r_ical::MarineNmeaBatch nmea_batch_{};
  uint32_t last_wave_nmea_ms_   = 0;
  uint32_t last_status_nmea_ms_ = 0;
  uint32_t last_temp_nmea_ms_   = 0;
#else
  atoms3r_ical::SerialBuffer<1024> debug_line_;
#endif

  Fusion fusion_{};
  bool output_dirty_ = false, have_output_ = false;

  atoms3r_ical::AtomS3RNavigationImu navigation_imu_{};
  atoms3r_ical::ImuTemperatureCache temperature_{};
  ins::MagFreshGate mag_gate_{};

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

  WaveDirectionReport wave_{};

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
  ins::SampleDtTracker sample_dt_{1.0f / LOOP_HZ};
  ins::SampleRateMeter rates_{};

  ImuSpikeFilter imu_spike_{};

  ins::StillGyroBiasEma gyro_bias_{ROT_BIAS_TAU_S, ROT_STILL_G_TOL_FRAC, ROT_STILL_GYRO_RAD_S};
  ins::RateOfTurnFilter rot_{};
  bool acc_bias_estimating_ = false;

private:
  LoopTiming timing_;
  SamplePacer pacer_{LOOP_PERIOD_US};

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

    ImuCalBlobV5 saved{};
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

    navigation_imu_.configure(M5.Imu);
    temperature_ = {};
    mag_gate_.reset();
    rates_.reset(millis());
    sample_dt_.reset();
    imu_spike_.reset();
  }

  void resetFusion_() {
#if SEA_STATE_SERIAL_NMEA
    nmea_batch_.reset();
#else
    debug_line_.reset();
#endif
    output_dirty_ = false;
    have_output_ = false;
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

    fcfg.sigma_a = sigma_a;
    fcfg.sigma_g = sigma_g;
    fcfg.sigma_m = sigma_m;

    fcfg.mag_delay_sec = 0.0f;
    fcfg.mag_init_min_mag_norm = 5.0f;

    fcfg.gravity_magnitude = g_local;

    // Saved calibration can be supplemented by the screened gyro-mag fit.
    fcfg.mag_screen_installation = seastate::common::defaults::SCREEN_SAVED_MAG_INSTALLATION;
    fcfg.mag_continuous_hard_iron = fcfg.mag_screen_installation;

    fcfg.enable_displacement_detrend = true;

    fusion_.begin(fcfg);

    auto& ff = fusion_.raw();
    ff.enableTuner(true);
    // Without S states there is no S = 0 pseudo-measurement, so OU-III's
    // pseudo-update cadence knob has no counterpart here.
    ff.setWithMag(true);
    // The tuner coefficients (S_factor, tau, sigma, CARMA sigma gain/floor)
    // and the CARMA shapes are deliberately not overridden here.  The header
    // defaults are the fixed operating-point mapping every benchmarked CARMA
    // filter uses (docs/carma-three-axis.md); overriding them would put the
    // device at an operating point no evidence covers.
    ff.setAccNoiseFloorSigma(ACC_NOISE_FLOOR_SIGMA_DEFAULT);
    // Wave excitation from measured evidence instead of the legacy
    // 0.05 m/s^2 floor (docs/carma-wave-excitation.md): the CARMA prior is
    // driven by the noise-debiased band variance, weighted by the evidence
    // that waves are present, and its frequency is held while the measured
    // wave period is not confirmed by the acceleration band.  In still water
    // the prior's driving noise is then at its 1e-5 m/s^2 numerical floor
    // instead of 0.05 m/s^2 at 0.05 Hz.  CarmaExcitationMode::Legacy restores
    // the benchmarked behaviour of the documents above.
    ff.setCarmaExcitationMode(ocean_imu::carma::CarmaExcitationMode::Evidence);
    ff.setCarmaSensorNoise(CARMA_ACC_WHITE_SIGMA_MPS2, CARMA_ACC_BIAS_RW_MPS2_SQRTS,
                           CARMA_ACC_NOISE_REL_UNCERTAINTY);
    ff.enableClamp(true);
    ff.setFreqInputCutoffHz(6.0f);

    // IMU installation lever-arm (off-CoG) correction.
    //
    // Set explicitly on every reset, so the adjustment point is visible; the
    // wrapper keeps it across fusion_.begin().  Unlike OU-III, whose
    // correction lives in its MEKF's accelerometer update, the CARMA wrapper
    // removes the term at its input, before the vibration guard: the startup
    // tilt proxy, the magnetic reference it learns, the tuner and the
    // wave-direction estimator then see CoG specific force, as the filter
    // does.  docs/carma-imu-lever-arm.md measured both placements: inside the
    // MEKF the rest of the wrapper still sees the IMU's own motion, and tilt
    // suffers in calm seas.  Do not also set ff.mekf().set_imu_lever_arm_body():
    // the term would be removed twice.
    //
    // The vector below is zero by default, which tells the filter the IMU sits
    // on the centre of gravity and turns the correction off (identical to
    // ff.clearImuLeverArm()).  Edit IMU_LEVER_ARM_*_M near the top of this
    // sketch to enter the real installation offset -- nothing else has to
    // change here.  Related knob, if the correction is enabled and the gyro is
    // noisy: ff.setImuLeverArmCutoffHz(hz) sets the low-pass corner ahead of
    // the angular-acceleration derivative (default 15 Hz).  A lower corner
    // passes less differentiated gyro noise but lags more; the study's corner
    // sweep shows both ways to lose.
    ff.setImuLeverArmBody(
        Vector3f(IMU_LEVER_ARM_X_M, IMU_LEVER_ARM_Y_M, IMU_LEVER_ARM_Z_M));

    rot_.reset();
    gyro_bias_.reset();
    acc_bias_estimating_ = false;

    heading_deg_   = NAN;
    heading_valid_ = false;
    heading_fused_ = false;

    temperature_ = {};
    mag_gate_.reset();
    mag_ok_           = false;
    mag_fresh_        = false;
    mag_norm_uT_      = NAN;
    heading_mag_deg_  = NAN;
    imu_temp_c_       = NAN;
    heading_mag_ok_   = false;

    stale_frame_count_ = 0;
    sample_dt_.reset();
    imu_spike_.reset();
    heave_speed_mps_     = 0.0f;
    heave_raw_m_         = 0.0f;
    heave_baseline_m_    = 0.0f;
    heave_wave_raw_m_    = 0.0f;
    heave_wave_clean_m_  = 0.0f;

    wave_ = WaveDirectionReport{};
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
    requireUsableCalibration(ui_, store_, blob_, have_blob_, runtime_, kSeaStateCalibrationWizard);
    reinitImu_();
    resetFusion_();
  }

  void handleRunWizard_() {
    Serial.println("[HOME] RUN WIZARD");

    const bool saved = runWizardFlow_(false);
    if (!saved) {
      ui_.notSavedNotice();
    }

    requireUsableCalibration(ui_, store_, blob_, have_blob_, runtime_, kSeaStateCalibrationWizard);
    reinitImu_();
    resetFusion_();
  }

  void drawHomePending_() {
    display_.waitForMenu();
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

  // Before INS readiness, report the measured compass using startup tilt.
  // Once the magnetically informed filter is Live, publish its fused yaw.
  void updateCompassHeading_(const Eigen::Quaternionf& q_bw, bool attitude_ok,
                               bool fused_ready = false, float fused_heading_deg = NAN) {
    heading_mag_ok_ = false;
    heading_mag_deg_ = NAN;
    if (mag_ok_ && attitude_ok) {
      const Vector3f down_b = ins::quatRotate(q_bw.conjugate(), Vector3f::UnitZ());
      heading_mag_ok_ = ins::magneticHeadingFromDownAndMagBody(
          down_b, m_cal_, heading_mag_deg_);
    }
    heading_fused_ = fused_ready && attitude_ok && std::isfinite(fused_heading_deg);
    heading_valid_ = heading_fused_ || heading_mag_ok_;
    heading_deg_ = heading_fused_ ? ins::wrap360Deg(fused_heading_deg)
        : (heading_mag_ok_ ? heading_mag_deg_ : NAN);
  }

  void updateFilter_(const ImuSample& s, uint32_t poll_ms) {
    dt_ = sample_dt_.update(s.sample_us);
    const float tempC = s.tempC;
    imu_temp_c_ = tempC;
    rates_.countImu();

    Vector3f a_raw = s.a;
    Vector3f w_raw = s.w;

    a_raw_norm_ = a_raw.norm();
    imu_spike_.apply(a_raw, w_raw, dt_);

    a_cal_ = runtime_.applyAccel(a_raw, tempC);
    w_cal_ = runtime_.applyGyro(w_raw, tempC);
    m_cal_ = runtime_.applyMag(s.m);

    mag_norm_uT_ = m_cal_.norm();
    mag_ok_ = std::isfinite(mag_norm_uT_) && (mag_norm_uT_ > 5.0f) && (mag_norm_uT_ < 200.0f);
    mag_fresh_ = mag_gate_.update(mag_ok_, poll_ms);
    if (mag_fresh_) rates_.countMag();

    fusion_.update(dt_, w_cal_, a_cal_);
    if (mag_ok_ && mag_fresh_) {
      fusion_.updateMag(m_cal_);
    }

    gyro_bias_.update(w_cal_, a_cal_, g_local, dt_);
    const Eigen::Quaternionf q_bw = fusion_.isLive()
        ? fusion_.attitudeQuat() : fusion_.raw().startupProxyQuat();
    const float world_z = ins::quatRotate(q_bw.normalized(), gyro_bias_.corrected(w_cal_)).z();
    rot_.update(world_z, dt_);
    rates_.update(millis());
    output_dirty_ = true;
    if (!have_output_) refreshOutput_();
  }

  void refreshOutput_() {
    if (!output_dirty_) return;
    output_dirty_ = false;
    have_output_ = true;
    // Keep startup attitude and startup compass tilt explicitly on the
    // yaw-free Mahony proxy. Its yaw is unobservable and must never become a
    // user-facing heading. Once Live, switch to the fused MEKF attitude.
    const bool live = fusion_.isLive();
    Eigen::Quaternionf q_bw = live
        ? fusion_.attitudeQuat()
        : fusion_.raw().startupProxyQuat();
    const bool have_attitude = live || fusion_.raw().startupProxyInitialized();
    float roll_est_deg = roll_deg_;
    float pitch_est_deg = pitch_deg_;
    float heading_est_deg = heading_deg_;
    const bool attitude_ok = have_attitude &&
        ins::rollPitchHeadingFromQuatBw(q_bw, roll_est_deg, pitch_est_deg, heading_est_deg);
    if (attitude_ok) {
      q_bw.normalize();
      roll_deg_ = roll_est_deg;
      pitch_deg_ = pitch_est_deg;
    }
    // The Mahony proxy's yaw is useful for propagating roll/pitch through
    // motion, but it is not north.  Startup magnetic heading uses only the
    // proxy's yaw-free tilt so a proxy yaw branch/wrap can never flip north
    // and south.  Live still publishes the fused MEKF yaw.
    Eigen::Quaternionf q_compass_tilt = q_bw;
    if (!live) q_compass_tilt = fusion_.raw().startupProxyTiltQuat();
    updateCompassHeading_(q_compass_tilt, attitude_ok,
        live && fusion_.hasMagNorthLock(), heading_est_deg);

    const Vector3f velocity_ned_mps = fusion_.raw().mekf().get_velocity();
    heave_speed_mps_ = -velocity_ned_mps.z();

    const Vector3f displacement_raw_up_m = fusion_.displacementUpMeters();
    const auto& displacement_det_out = fusion_.displacementDetrend();

    heave_m_         = displacement_raw_up_m.z();
    wave_envelope_m_ = fusion_.raw().getDisplacementScale();
    wave_hz_         = fusion_.raw().getFreqHz();

    wave_ = WaveDirectionReport::from(fusion_.waveDirectionDeg(),
                                      fusion_.raw().getDirSignState(),
                                      fusion_.isLive());
    acc_bias_estimating_ = fusion_.isLive();

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
    if (now_ms - last_ui_ms_ < atoms3r_timing::UI_REFRESH_MS) return;
    last_ui_ms_ = now_ms;
    refreshOutput_();
    display_.publish({heading_deg_, heading_mag_deg_, roll_deg_, pitch_deg_, heave_m_,
                      wave_.dir_ok ? wave_.dir_deg : wave_.axis_deg, mag_norm_uT_,
                      a_raw_norm_, a_cal_.norm(), wave_.sign_raw, heading_valid_,
                      heading_mag_ok_, mag_ok_, mag_fresh_});
  }

  void streamSerial_() {
    const uint32_t now_ms = millis();
#if SEA_STATE_SERIAL_NMEA
    if (!nmea_batch_.readyForSnapshot(Serial, SEA_STATE_NMEA_TALKER, now_ms, last_serial_ms_)) return;
#else
    if (!debug_line_.ready(Serial, now_ms - last_serial_ms_ >= atoms3r_timing::DEBUG_SERIAL_MS)) return;
    last_serial_ms_ = now_ms;
#endif
    refreshOutput_();

#if SEA_STATE_SERIAL_NMEA
    const bool valid = fusion_.isLive() && heading_fused_;
    if (heading_valid_) {
      nmea_hdm(SEA_STATE_NMEA_TALKER, heading_deg_, nmea_batch_.output());
    }
    nmea_batch_.begin(heave_wave_clean_m_, heave_speed_mps_, rot_.dpm(), valid);
    nmea_batch_.addTilt(pitch_deg_, roll_deg_);

    if (now_ms - last_wave_nmea_ms_ >= atoms3r_timing::NMEA_WAVE_MS) {
      last_wave_nmea_ms_ = now_ms;
      nmea_batch_.addWave(wave_, fusion_.isLive());
    }
    if (now_ms - last_temp_nmea_ms_ >= atoms3r_timing::NMEA_TEMP_MS) {
      last_temp_nmea_ms_ = now_ms;
      nmea_batch_.addTemperature(imu_temp_c_);
    }
    if (now_ms - last_status_nmea_ms_ >= atoms3r_timing::NMEA_STATUS_MS) {
      last_status_nmea_ms_ = now_ms;
      nmea_batch_.addStatus({fusion_.isLive() && std::isfinite(heave_wave_clean_m_),
                             fusion_.hasMagNorthLock(), gyro_bias_.learning(), acc_bias_estimating_,
                             rates_.imuHz(), rates_.magHz()});
    }
    nmea_batch_.output().drain(Serial);
#else
  #if ARDUINO_PLOTTER
    debug_line_.printf(
      "HrawCm:%+.3f\tHwaveEnvelopeCm:%+.3f\tHwaveCleanCm:%+.3f\tWavAxisDeg:%+.2f\tWavDirDeg:%+.2f\tWavSign:%d\n",
      static_cast<double>(heave_raw_m_ * 100.0f),
      static_cast<double>(wave_envelope_m_ * 100.0f),
      static_cast<double>(heave_wave_clean_m_ * 100.0f),
      static_cast<double>(wave_.axis_deg),
      static_cast<double>(wave_.dir_deg),
      wave_.sign_raw);
  #else
    debug_line_.printf(
      "hdg=%s%.2f | hdg_mag=%s%.2f | wav_axis=%.2f | wav_dir=%s%.2f | wav_sign=%d pol=%d"
      " | dt_imu_ms=%.2f | mag_ok=%u | mag_fresh=%u | |m|=%.2f"
      " | acc_ned[m/s^2] N=%.3f E=%.3f D=%.3f"
      " | gyro_ned[rad/s] N=%.3f E=%.3f D=%.3f"
      " | mag_ned[uT] N=%.2f E=%.2f D=%.2f\n",
      heading_valid_ ? "" : "~",
      static_cast<double>(heading_deg_),
      heading_mag_ok_ ? "" : "~",
      static_cast<double>(heading_mag_deg_),
      static_cast<double>(wave_.axis_deg),
      wave_.dir_ok ? "" : "~",
      static_cast<double>(wave_.dir_deg),
      wave_.sign_raw,
      wave_.sign_polarity,
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
    debug_line_.drain(Serial);
#endif
  }
};

static FusionApp g_app;

void setup() { g_app.begin(); }
void loop()  { g_app.tick(); }
