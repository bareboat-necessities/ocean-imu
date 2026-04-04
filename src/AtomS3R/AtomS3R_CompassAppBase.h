#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy
*/

#include <Arduino.h>
#include <M5Unified.h>
#include <memory>
#include <new>

#include "AtomS3R/AtomS3R_ImuCal.h"
#include "AtomS3R/AtomS3R_ImuCalWizard.h"
#include "AtomS3R/AtomS3R_CompassUI.h"
#include "nmea/NmeaCompass.h"

namespace atoms3r_compass {

using namespace atoms3r_ical;
using Vector3f = Eigen::Matrix<float, 3, 1>;

static constexpr float LOOP_HZ = 200.0f;
static constexpr uint32_t LOOP_PERIOD_US = (uint32_t)(1000000.0f / LOOP_HZ);

static constexpr uint32_t UI_REFRESH_MS = 100;
static constexpr uint32_t DEBUG_SERIAL_MS = 200;
static constexpr uint32_t NMEA_SERIAL_MS = 100;

static constexpr float ROT_BIAS_TAU_S = 20.0f;
static constexpr float ROT_STILL_G_TOL_FRAC = 0.12f;
static constexpr float ROT_STILL_GYRO_RAD_S = 0.15f;

static inline float clampf_(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

static inline float wrap360_(float deg) {
  while (deg < 0.0f) deg += 360.0f;
  while (deg >= 360.0f) deg -= 360.0f;
  return deg;
}

static inline float wrap180_(float deg) {
  while (deg < -180.0f) deg += 360.0f;
  while (deg >= 180.0f) deg -= 360.0f;
  return deg;
}

static inline Vector3f quatRotateBodyToWorld_(float x, float y, float z, float w, const Vector3f& v_body) {
  Vector3f q(x, y, z);
  Vector3f t = 2.0f * q.cross(v_body);
  return v_body + w * t + q.cross(t);
}

struct CalibratedSample {
  float dt = 0.0f;
  float a_raw_norm = 0.0f;

  Vector3f a_cal = Vector3f::Zero();
  Vector3f w_cal = Vector3f::Zero();
  Vector3f m_cal = Vector3f::Zero();

  bool mag_ok = false;
  bool mag_fresh = false;
  float mag_norm_uT = 0.0f;
  Vector3f m_unit = Vector3f::Zero();
};

struct AttitudeSolution {
  bool valid = false;
  float x = 0.0f;
  float y = 0.0f;
  float z = 0.0f;
  float w = 1.0f;

  float roll_deg = 0.0f;
  float pitch_deg = 0.0f;
  float yaw_deg = 0.0f;
  float heading_deg = 0.0f;
};

struct CompassOutputs {
  AttitudeSolution att;
  float rot_dpm = 0.0f;
};

static inline AttitudeSolution makeAttitudeFromQuat(float x, float y, float z, float w) {
  AttitudeSolution out{};

  const float nn = x * x + y * y + z * z + w * w;
  if (nn <= 1e-12f) {
    out.valid = false;
    out.x = 0.0f;
    out.y = 0.0f;
    out.z = 0.0f;
    out.w = 1.0f;
    out.roll_deg = 0.0f;
    out.pitch_deg = 0.0f;
    out.yaw_deg = 0.0f;
    out.heading_deg = 0.0f;
    return out;
  }

  const float invn = 1.0f / sqrtf(nn);
  x *= invn;
  y *= invn;
  z *= invn;
  w *= invn;

  out.valid = true;
  out.x = x;
  out.y = y;
  out.z = z;
  out.w = w;

  const float siny_cosp = 2.0f * (w * z + x * y);
  const float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
  const float yaw = atan2f(siny_cosp, cosy_cosp);

  float sinp = 2.0f * (w * y - z * x);
  sinp = clampf_(sinp, -1.0f, 1.0f);
  const float pitch = asinf(sinp);

  const float sinr_cosp = 2.0f * (w * x + y * z);
  const float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
  const float roll = atan2f(sinr_cosp, cosr_cosp);

  out.roll_deg = wrap180_(roll * RAD_TO_DEG);
  out.pitch_deg = pitch * RAD_TO_DEG;
  out.yaw_deg = yaw * RAD_TO_DEG;
  out.heading_deg = wrap360_(out.yaw_deg);
  return out;
}

struct MagGateConfig {
  // Magnetometer ODR is 25 Hz => nominal period is 40 ms.
  // Use a slightly smaller spacing so a real new sample is usually accepted once,
  // while repeated polling of the same sample in a 200 Hz loop is suppressed.
  uint32_t sample_spacing_ms = 35;

  // Used only for soft "stuck" detection, not for rejecting freshness.
  float min_delta_uT = 0.001f;

  // 0 disables stuck flagging.
  uint8_t stuck_reject_after_n = 20;
};

class MagGate {
 public:
  explicit MagGate(const MagGateConfig& cfg) : cfg_(cfg) {}

  void reset() {
    last_mag_.setZero();
    last_ms_ = 0;
    repeat_count_ = 0;
  }

  bool update(const Vector3f& m_cal, bool mag_ok, uint32_t now_ms) {
    if (!mag_ok) {
      repeat_count_ = 0;
      return false;
    }

    // First valid sample: accept immediately.
    if (last_ms_ == 0) {
      last_ms_ = now_ms;
      last_mag_ = m_cal;
      repeat_count_ = 0;
      return true;
    }

    // Accept at most one sample per configured spacing.
    // This is a proxy for "new mag sample" when actual mag-update status is unknown.
    const uint32_t dtm = now_ms - last_ms_;
    if (dtm < cfg_.sample_spacing_ms) return false;

    const float dm = (m_cal - last_mag_).norm();
    const bool looks_stuck = dm <= cfg_.min_delta_uT;

    // Always accept a time-spaced valid sample as fresh.
    // Tiny deltas are tracked only as a soft diagnostic.
    if (looks_stuck) {
      if (repeat_count_ < 255) repeat_count_++;
    } else {
      repeat_count_ = 0;
    }

    last_ms_ = now_ms;
    last_mag_ = m_cal;
    return true;
  }

  bool looksStuck() const {
    return (cfg_.stuck_reject_after_n != 0) &&
           (repeat_count_ >= cfg_.stuck_reject_after_n);
  }

 private:
  MagGateConfig cfg_{};
  Vector3f last_mag_ = Vector3f::Zero();
  uint32_t last_ms_ = 0;
  uint8_t repeat_count_ = 0;
};

class RotEstimator {
 public:
  void reset() {
    inited_ = false;
    filt_dpm_ = 0.0f;
    gyro_bias_ok_ = false;
    gyro_bias_ema_.setZero();
  }

  float update(float dt, const Vector3f& a_cal, const Vector3f& w_cal, const AttitudeSolution& att) {
    if (!att.valid) return filt_dpm_;

    const float g = ImuCalCfg::g_std;
    const float a_err = fabsf(a_cal.norm() - g);
    const bool still = (a_err < ROT_STILL_G_TOL_FRAC * g) && (w_cal.norm() < ROT_STILL_GYRO_RAD_S);

    if (still) {
      const float alpha_b = 1.0f - expf(-dt / ROT_BIAS_TAU_S);
      if (!gyro_bias_ok_) {
        gyro_bias_ok_ = true;
        gyro_bias_ema_ = w_cal;
      } else {
        gyro_bias_ema_ += alpha_b * (w_cal - gyro_bias_ema_);
      }
    }

    Vector3f w_use = w_cal;
    if (gyro_bias_ok_) w_use -= gyro_bias_ema_;

    Vector3f w_world = quatRotateBodyToWorld_(att.x, att.y, att.z, att.w, w_use);
    float rot_dpm = w_world.z() * RAD_TO_DEG * 60.0f;
    rot_dpm = clampf_(rot_dpm, -720.0f, 720.0f);

    const float alpha_r = 1.0f - expf(-dt / 1.5f);
    if (!inited_) {
      inited_ = true;
      filt_dpm_ = rot_dpm;
    } else {
      filt_dpm_ += alpha_r * (rot_dpm - filt_dpm_);
    }

    return filt_dpm_;
  }

 private:
  bool inited_ = false;
  float filt_dpm_ = 0.0f;
  bool gyro_bias_ok_ = false;
  Vector3f gyro_bias_ema_ = Vector3f::Zero();
};

class IAttitudeBackend {
 public:
  virtual ~IAttitudeBackend() = default;
  virtual void reset() = 0;
  virtual void step(const CalibratedSample& s, AttitudeSolution& out) = 0;
  virtual bool isValid() const = 0;
};

class CompassAppBase {
 public:
  CompassAppBase(std::unique_ptr<IAttitudeBackend> backend, MagGateConfig mag_cfg, const char* boot_name)
      : backend_(std::move(backend)), wizard_(ui_, store_), mag_gate_(mag_cfg), boot_name_(boot_name) {}

  void begin() {
    Serial.begin(115200);
    delay(150);
    Serial.println();
    Serial.printf("[BOOT] AtomS3R Compass + Cal Wizard (%s)\n", boot_name_);

    auto cfg = M5.config();
    M5.begin(cfg);

    clearM5UnifiedImuCalibration();
    delay(250);

    ui_.begin();

    if (use_graphics_) {
      ui_.setReadRotation();
      compass_ui_.begin();
      compass_ui_ready_ = compass_ui_.ok();
      if (!compass_ui_ready_) use_graphics_ = false;
    }

    if (!M5.Imu.isEnabled()) {
      Serial.println("[BOOT] IMU not found / not enabled");
      ui_.fail("IMU", "Not found");
      while (true) delay(100);
    }

    reloadBlobAndRuntime_();

    if (!have_blob_) {
      Serial.println("[BOOT] No saved calibration. Starting wizard...");
      ImuCalBlobV1 saved{};
      if (wizard_.runAndSave(saved)) {
        Serial.println("[BOOT] Wizard saved calibration. Loaded:");
        printBlobSummary(Serial, saved);
        printBlobDetail(Serial, saved);
        blob_ = saved;
        have_blob_ = true;
        runtime_.rebuildFromBlob(blob_);
      } else {
        Serial.println("[BOOT] Wizard did not save calibration. Running with raw values.");
      }
    } else {
      Serial.println("[BOOT] Found saved calibration:");
      printBlobSummary(Serial, blob_);
      printBlobDetail(Serial, blob_);
    }

    resetPipeline_();
    drawHomeStatic_();

    start_us_ = micros();
    next_tick_us_ = micros();
    last_update_us_ = 0;
  }

  void tick() {
    // Wait until next scheduled tick.
    while (true) {
    const uint32_t now_us = micros();
    const int32_t wait_us = (int32_t)(next_tick_us_ - now_us);
    if (wait_us <= 0) break;

    if (wait_us > 1000)
      delayMicroseconds(500);
    else
      delayMicroseconds((uint32_t)wait_us);
    }

    // Advance schedule. If we fell far behind, resync cleanly.
    const uint32_t now_us2 = micros();
    if ((int32_t)(now_us2 - next_tick_us_) > (int32_t)(4 * LOOP_PERIOD_US))
    next_tick_us_ = now_us2 + LOOP_PERIOD_US;
    else
    next_tick_us_ += LOOP_PERIOD_US;
      
    Input::update();

    if (Input::tapPressed()) {
      tap_count_++;
      tap_deadline_ms_ = millis() + M5UiCfg::MENU_TAP_WINDOW_MS;
      drawHomePending_();
      Serial.printf("[TAP] count=%d\n", tap_count_);
    }

    if (tap_count_ > 0 && (int32_t)(millis() - tap_deadline_ms_) > 0) {
      if (tap_count_ >= 3)
        handleErase_();
      else
        handleRunWizard_();
      tap_count_ = 0;
      tap_deadline_ms_ = 0;
      drawHomeStatic_();
    }

    const uint32_t sample_us = micros();
    const uint32_t update_mask = M5.Imu.update();

    ImuSample s;
    if (readImuMapped(M5.Imu, update_mask, sample_us, s)) updateOutputs_(s);

    updateUI_();
    streamSerial_();
  }

 protected:
  void reloadBlobAndRuntime_() {
    have_blob_ = store_.load(blob_);
    if (!have_blob_) memset(&blob_, 0, sizeof(blob_));
    runtime_.rebuildFromBlob(blob_);
  }

  void handleErase_() {
    Serial.println("[HOME] triple tap => ERASE");
    if (!ui_.eraseConfirm()) {
      Serial.println("[HOME] erase cancelled");
      return;
    }

    store_.erase();
    clearM5UnifiedImuCalibration();
    reloadBlobAndRuntime_();
    resetPipeline_();

    Serial.println("[HOME] erased blob + cleared M5Unified cal");
  }

  void handleRunWizard_() {
    Serial.println("[HOME] single tap => RUN WIZARD");

    ImuCalBlobV1 saved{};
    if (!wizard_.runAndSave(saved)) {
      ui_.notSavedNotice();
      return;
    }

    Serial.println("[HOME] new calibration saved:");
    printBlobSummary(Serial, saved);
    printBlobDetail(Serial, saved);

    blob_ = saved;
    have_blob_ = true;
    runtime_.rebuildFromBlob(blob_);
    resetPipeline_();
  }

  void resetPipeline_() {
    backend_->reset();
    mag_gate_.reset();
    rot_.reset();
    outputs_ = CompassOutputs{};
    sample_ = CalibratedSample{};
    last_update_us_ = 0;
  }

CalibratedSample makeCalibratedSample_(const ImuSample& s) {
  CalibratedSample out{};
  out.a_raw_norm = s.a.norm();

  out.a_cal = runtime_.applyAccel(s.a, s.tempC);
  out.w_cal = runtime_.applyGyro(s.w, s.tempC);
  out.m_cal = runtime_.applyMag(s.m);

  if (last_update_us_ == 0) {
    out.dt = 1.0f / LOOP_HZ;
  } else {
    out.dt = (s.sample_us - last_update_us_) * 1e-6f;
  }
  last_update_us_ = s.sample_us;
  out.dt = clampf_(out.dt, 0.0010f, 0.0200f);

  out.mag_norm_uT = out.m_cal.norm();
  out.mag_ok = (out.mag_norm_uT > 5.0f && out.mag_norm_uT < 200.0f);

  if (out.mag_ok && out.mag_norm_uT > 1e-6f) out.m_unit = out.m_cal / out.mag_norm_uT;

  out.mag_fresh = mag_gate_.update(out.m_cal, out.mag_ok, millis());
  return out;
}

  void updateOutputs_(const ImuSample& s) {
    sample_ = makeCalibratedSample_(s);
    backend_->step(sample_, outputs_.att);
    outputs_.rot_dpm = rot_.update(sample_.dt, sample_.a_cal, sample_.w_cal, outputs_.att);
  }

  void drawHomeStatic_() {
    ui_.setReadRotation();
    ui_.title("COMPASS");
    M5.Display.printf("BLOB: %s\n", have_blob_ ? "YES" : "NO");
    M5.Display.printf("A:%d G:%d M:%d\n", (int)runtime_.acc.ok, (int)runtime_.gyr.ok, (int)runtime_.mag.ok);
    ui_.line("Tap: calibrate");
    ui_.line("Tap x3: erase");
    ui_.line("");
  }

  void drawHomePending_() {
    ui_.setReadRotation();
    ui_.title("COMPASS");
    M5.Display.printf("Tap count: %d\n", tap_count_);
    ui_.line("");
    ui_.line("Wait...");
    ui_.line("1 tap=CAL");
    ui_.line("3 taps=ERASE");

    int32_t remain = (int32_t)(tap_deadline_ms_ - millis());
    remain = remain < 0 ? 0 : remain;
    float t01 = 1.0f - (float)remain / (float)M5UiCfg::MENU_TAP_WINDOW_MS;
    ui_.bar01(t01);
  }

  void updateUI_() {
    if (tap_count_ > 0) return;

    const uint32_t now_ms = millis();
    if (now_ms - last_ui_ms_ < UI_REFRESH_MS) return;
    last_ui_ms_ = now_ms;

    if (use_graphics_ && compass_ui_ready_) {
      ui_.setReadRotation();
      const bool tilt_warn = (fabsf(outputs_.att.roll_deg) > 35.0f) || (fabsf(outputs_.att.pitch_deg) > 35.0f);
      compass_ui_.draw(outputs_.att.heading_deg, sample_.mag_ok, sample_.mag_norm_uT, tilt_warn);
    } else {
      ui_.setReadRotation();
      ui_.title("COMPASS");

      M5.Display.printf("HDG: %6.1f deg\n", (double)outputs_.att.heading_deg);
      M5.Display.printf("ROL: %6.1f deg\n", (double)outputs_.att.roll_deg);
      M5.Display.printf("PIT: %6.1f deg\n", (double)outputs_.att.pitch_deg);

      M5.Display.printf("MAG: %s %s\n", sample_.mag_ok ? "OK " : "BAD", sample_.mag_fresh ? "NEW" : "OLD");
      M5.Display.printf("|m|: %6.1f uT\n", (double)sample_.mag_norm_uT);
      M5.Display.printf("|aR|:%5.2f |aC|:%5.2f\n", (double)sample_.a_raw_norm, (double)sample_.a_cal.norm());
      ui_.line("");
      M5.Display.printf("A:%d G:%d M:%d  B:%s\n", (int)runtime_.acc.ok, (int)runtime_.gyr.ok, (int)runtime_.mag.ok,
                        have_blob_ ? "YES" : "NO");
    }
  }

  void streamSerial_() {
    const uint32_t now_ms = millis();

#if COMPASS_SERIAL_NMEA
    if (now_ms - last_serial_ms_ < NMEA_SERIAL_MS) return;
#else
    if (now_ms - last_serial_ms_ < DEBUG_SERIAL_MS) return;
#endif
    last_serial_ms_ = now_ms;

#if COMPASS_SERIAL_NMEA
    nmea_hdm(COMPASS_NMEA_TALKER, outputs_.att.heading_deg);
    nmea_xdr_pitch_roll(COMPASS_NMEA_TALKER, outputs_.att.pitch_deg, outputs_.att.roll_deg);
    nmea_rot(COMPASS_NMEA_TALKER, outputs_.rot_dpm, backend_->isValid());
#else
    const float t = (micros() - start_us_) * 1e-6f;
    Serial.printf("t=%.2f ", (double)t);
    Serial.printf("HDG:%6.1f ", (double)outputs_.att.heading_deg);
    Serial.printf("ROLL:%6.1f PITCH:%6.1f ", (double)outputs_.att.roll_deg, (double)outputs_.att.pitch_deg);
    Serial.printf("mag:%s/%s |m|=%.1f ", sample_.mag_ok ? "OK" : "BAD", sample_.mag_fresh ? "NEW" : "OLD",
                  (double)sample_.mag_norm_uT);
    Serial.printf("|a_raw|=%.4f |a_cal|=%.4f ", (double)sample_.a_raw_norm, (double)sample_.a_cal.norm());
    Serial.printf("wC:%+.4f,%+.4f,%+.4f ", (double)sample_.w_cal.x(), (double)sample_.w_cal.y(), (double)sample_.w_cal.z());
    Serial.println();
#endif
  }

 protected:
  M5Ui ui_{};
  ImuCalStoreNvs store_{};
  ImuCalWizard wizard_;
  RuntimeCals runtime_{};
  CompassUI compass_ui_{};

  std::unique_ptr<IAttitudeBackend> backend_;
  MagGate mag_gate_;
  RotEstimator rot_;

  bool have_blob_ = false;
  ImuCalBlobV1 blob_{};

  bool use_graphics_ = (COMPASS_UI_DEFAULT_GRAPHICS != 0);
  bool compass_ui_ready_ = false;

  int tap_count_ = 0;
  uint32_t tap_deadline_ms_ = 0;
  uint32_t start_us_ = 0;
  uint32_t next_tick_us_ = 0;
  uint32_t last_update_us_ = 0;
  uint32_t last_ui_ms_ = 0;
  uint32_t last_serial_ms_ = 0;

  CalibratedSample sample_;
  CompassOutputs outputs_;

  const char* boot_name_ = "";
};

}  // namespace atoms3r_compass
