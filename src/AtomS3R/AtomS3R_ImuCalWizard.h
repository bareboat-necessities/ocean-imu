#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R IMU calibration wizard UI (Accel + Gyro + Mag) using imu_cal::* calibrators.

  Uses M5Unified IMU API as the IMU sample source.

  Accelerometer: ten guided static holds (six faces + four screen-up corners),
  bounded targeted extra holds, three later rechecks, and a full-matrix fit
  with thermal-slope qualification (imu_cal::AccelCalProcedure). In the full
  wizard the rechecks follow the gyro and magnetometer stages; the
  accelerometer-only mode keeps the saved gyro and magnetometer calibration.
  Nothing is saved unless the complete candidate validates; the previous
  calibration is kept otherwise.
*/

#include <Arduino.h>
#include <M5Unified.h>

#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits>
#include <algorithm>

#include <freertos/FreeRTOS.h>
#include <freertos/task.h>

#include "AtomS3R/AtomS3R_ImuCal.h"         // ImuSample, axis mapping conventions, blob/store/runtime helpers
#include "AtomS3R/AtomS3R_M5Ui.h"           // UI + Input + clamp01_
#include "imu_calibrate/CalibrateIMU.h"     // imu_cal::* + FitFail
#include "imu_calibrate/AccelCalCapture.h"  // accelerometer procedure (host-tested)

// Set to 1 to stream every raw accel/gyro sample as [ACCRAW] lines for
// tests/imu_calibrate/accel_cal-replay.
#ifndef ATOMS3R_ICAL_RAW_LOG
#define ATOMS3R_ICAL_RAW_LOG 0
#endif

namespace atoms3r_ical {

// Wizard configuration
struct ImuCalWizardCfg {
  // Step pacing (accelerometer holds: imu_cal::AccelCaptureCfg)
  static constexpr uint32_t PLACE_TIME_MS       = 6500;
  static constexpr uint32_t GYRO_TIMEOUT_MS     = 70000;
  static constexpr uint32_t MAG_TIMEOUT_MS      = 220000;
  static constexpr uint32_t STUCK_MS            = 12000;

  // Accelerometer observation capacity (blocks) and hold capacity
  static constexpr int ACCEL_MAX_OBS            = 340;
  static constexpr int ACCEL_MAX_HOLDS          = 24;

  // Sample goals (Mag buffer capacity is 400 below)
  static constexpr int GYRO_NEED                = 220;

  // MAG: target near-buffer-full, but spread over time
  static constexpr int MAG_NEED                 = 360;   // <= 400
  static constexpr int MAG_MIN_TO_FIT           = 220;

  // Minimums before fitting
  static constexpr int GYRO_MIN_TO_FIT          = 120;

  // MAG timing / downsample:
  static constexpr uint32_t MAG_SAMPLE_SPACING_MS = 80;     // accepted max ~12.5 Hz
  static constexpr uint32_t MAG_MIN_TIME_MS       = 45000;  // require >= 45 seconds of motion

  // MAG stale-repeat reject
  static constexpr float MAG_MIN_DELTA_uT       = 0.03f;

  // MAG coverage requirements (ratio-based, unitless)
  static constexpr float MAG_SPAN_MIN_FRAC      = 0.35f;
  static constexpr float MAG_SPAN_MID_FRAC      = 0.55f;

  // Also require centered direction range per axis (0..2). 1.05 means roughly reaching ±0.525.
  static constexpr float MAG_URANGE_TARGET      = 1.05f;

  // FreeRTOS stack:
  static constexpr uint32_t FIT_STACK_BYTES     = 32768;
  static constexpr uint32_t FIT_STACK_WORDS     = FIT_STACK_BYTES / sizeof(StackType_t);

  static constexpr uint32_t FIT_TIMEOUT_MS      = 30000;
};

// Small helpers
static inline uint8_t rot_add_(uint8_t base, int delta) {
  int r = (int)base + delta;
  r %= 4;
  if (r < 0) r += 4;
  return (uint8_t)r;
}

static inline float second_smallest3_(float a, float b, float c) {
  // median of 3 (2nd smallest)
  if (a > b) { float t=a; a=b; b=t; }
  if (b > c) { float t=b; b=c; c=t; }
  if (a > b) { float t=a; a=b; b=t; }
  return b;
}

static inline bool finite3_(const Vector3f& v) {
  return isfinite(v.x()) && isfinite(v.y()) && isfinite(v.z());
}

// "Planarity" coverage measure for MAG capture (determinant of unit-direction covariance).
static inline float unit_dir_cov_det_(const Vector3f* x, int n) {
  if (!x || n < 20) return 0.0f;

  // Mean-center
  Vector3f mu = Vector3f::Zero();
  int n_mu = 0;
  for (int i = 0; i < n; ++i) {
    if (!finite3_(x[i])) continue;
    mu += x[i];
    ++n_mu;
  }
  if (n_mu < 20) return 0.0f;
  mu *= 1.0f / (float)n_mu;

  // C = mean(u u^T), u = (x-mu)/||x-mu||
  Matrix3f C = Matrix3f::Zero();
  int m = 0;
  for (int i = 0; i < n; ++i) {
    if (!finite3_(x[i])) continue;
    Vector3f d = x[i] - mu;
    float dn = d.norm();
    if (!(dn > 1e-6f)) continue;
    Vector3f u = d / dn;
    if (!finite3_(u)) continue;
    C.noalias() += u * u.transpose();
    ++m;
  }
  if (m < 20) return 0.0f;

  C *= 1.0f / (float)m;
  float detC = C.determinant();
  if (!isfinite(detC)) return 0.0f;
  return detC;
}

class ImuCalWizard {
public:
  ImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store)
  : ui_(ui), store_(store) {}

  // Runs the wizard and saves to NVS. Returns true only if saved successfully.
  // out_saved is filled with the saved blob (readback-validated). On any
  // failure or abort nothing is written and the previous calibration stays.
  bool runAndSave(ImuCalBlobV3& out_saved) {
    Serial.println("[WIZ] start");

    for (;;) {
      bool redo_all = false;

      // Prevent "stacking" with any M5Unified offsets
      clearM5UnifiedImuCalibration();

      gyroCal_.clear();
      magCal_.clear();
      configureCalibrators_();

      gyr_out_ = imu_cal::GyroCalibration<float>{};
      mag_out_ = imu_cal::MagCalibration<float>{};

      // Previous calibration: kept on failure, source of the preserved gyro/mag
      // (accelerometer-only mode), a compatible thermal slope, and the
      // stationary gyro level used by the rotation gate.
      ImuCalBlobV3 prev{};
      const bool have_prev = store_.load(prev);
      const bool can_accel_only = have_prev && prev.gyro_ok;

      const M5Ui::StartAction start = ui_.startMenu(can_accel_only);
      if (start == M5Ui::StartAction::CANCEL) {
        Serial.println("[WIZ] cancelled; previous calibration kept");
        return false;
      }
      const bool accel_only = (start == M5Ui::StartAction::ACCEL_ONLY);
      Serial.printf("[WIZ] mode=%s have_prev=%d\n", accel_only ? "accel_only" : "full", (int)have_prev);

      sensorIdentity_();
      const imu_cal::AccelThermalPrior prior =
          have_prev ? accelThermalPriorFrom(prev, sensor_id_lo_, sensor_id_hi_, imu_type_) : imu_cal::AccelThermalPrior{};
      Serial.printf("[ACC] thermal prior: %s\n", prior.valid ? "compatible slope available" : "none");
      const Vector3f gyro_level = have_prev && prev.gyro_ok
          ? Vector3f(prev.gyro_b0[0], prev.gyro_b0[1], prev.gyro_b0[2]) : Vector3f::Zero();
      // Replay context (tests/imu_calibrate/accel_cal-replay).
      Serial.printf("[ACCMODE] %s g=%.7f\n", accel_only ? "accel_only" : "full", (double)accel_fcfg_.g);
      Serial.printf("[ACCPRIOR] %d,%.7f,%.7f,%.7f,%.2f,%.2f,%.2f,%.2f\n", (int)prior.valid, prior.k[0], prior.k[1],
                    prior.k[2], prior.k_temp_lo, prior.k_temp_hi, prior.clamp_lo, prior.clamp_hi);
      Serial.printf("[ACCGYRO] %.7f,%.7f,%.7f,%d\n", (double)gyro_level.x(), (double)gyro_level.y(),
                    (double)gyro_level.z(), (int)(have_prev && prev.gyro_ok));

      accel_.begin(accel_ccfg_, accel_fcfg_, prior, gyro_level, have_prev && prev.gyro_ok);
      AccelIo io(*this);

      if (!accel_.runMainStage(io)) return accelFail_();

      if (!accel_only) {
        if (!captureGyro_()) return false;
        if (gyroCal_.buf.n < ImuCalWizardCfg::GYRO_MIN_TO_FIT) {
          ui_.fail("GYRO", "Too few accepted");
          return false;
        }
        if (!runFitTask_(FitKind::GYRO, "GYRO", true)) return false;
        accel_.setGyroReference(gyr_out_.biasT.b0);
        Serial.printf("[ACCGYRO] %.7f,%.7f,%.7f,1\n", (double)gyr_out_.biasT.b0.x(), (double)gyr_out_.biasT.b0.y(),
                      (double)gyr_out_.biasT.b0.z());

        // If mag is unavailable, skip MAG stage cleanly.
        if (!magAvailable_()) {
          Serial.println("[MAG] unavailable -> skipping mag calibration");
          mag_out_.ok = false;
        } else if (!runMagStage_(redo_all)) {
          if (redo_all) {
            Serial.println("[WIZ] redo all requested");
            continue;
          }
          return false;
        }
      }

      // Later rechecks (thermal evidence + held-out verification), final fit.
      if (!accel_.runRecheckStage(io)) return accelFail_();
      if (!accel_.runFinalFit(io)) return accelFail_();

      // Candidate: new accelerometer set; gyro/mag from this run, or carried
      // unchanged from the previous calibration in accelerometer-only mode.
      imu_cal::AccelCalibration<float> fc;
      AccelProc::Fitter::toFloat(accel_.result(), accel_fcfg_.g, fc);
      const uint32_t capture_s = accel_.totalHoldMs() / 1000u;
      ImuCalBlobV3 blob;
      if (accel_only) {
        blob = accelOnlyCandidate(prev, accel_.result(), fc, capture_s, sensor_id_lo_, sensor_id_hi_, imu_type_);
      } else {
        fillGyroMag_(blob);
        fillAccelFromFit(blob, accel_.result(), fc, capture_s);
        blob.sensor_id_lo = sensor_id_lo_;
        blob.sensor_id_hi = sensor_id_hi_;
        blob.imu_type = imu_type_;
      }

      // The stored float set, rebuilt through the runtime path, must reproduce
      // the fit before it is written...
      if (!validateStored_(blob, "candidate")) {
        ui_.fail("SAVE", "Float check failed");
        return false;
      }

      ui_.setReadRotation();
      ui_.title("SAVE");
      ui_.line("Writing...");
      ImuCalBlobV3 rb{};
      const bool saved = store_.saveVerified(blob, rb);
      // ...and again after the read-back, which must be byte-identical to it.
      const bool rb_ok = saved && validateStored_(rb, "readback");
      Serial.printf("[SAVE] verified=%d readback_valid=%d\n", (int)saved, (int)rb_ok);
      if (!saved || !rb_ok) {
        ui_.fail("SAVE", "Write/readback fail");
        return false;
      }

      out_saved = rb;
      showDone_(rb);
      Serial.println("[WIZ] done");
      return true;
    }
  }

private:
  using AccelProc = imu_cal::AccelCalProcedure<ImuCalWizardCfg::ACCEL_MAX_OBS, ImuCalWizardCfg::ACCEL_MAX_HOLDS>;

  // Screens, samples and the fit task for the accelerometer procedure.
  class AccelIo : public imu_cal::AccelCalIo {
  public:
    explicit AccelIo(ImuCalWizard& w) : w_(w) {}

    bool prep(const imu_cal::AccelStepView& v) override {
      M5Ui& ui = w_.ui_;
      ui.setReadRotation();
      ui.title(v.title);
      ui.line(v.label);
      for (int j = 0; j < 3; ++j) if (v.lines[j]) ui.line(v.lines[j]);
      if (v.note[0]) ui.line(v.note);
      ui.line("");
      ui.line("Tap then place");
      ui.line("Tap BtnA");
      while (true) {
        Input::update();
        if (Input::tapPressed()) break;
        delay(10);
      }
      drawn_ = false;
      Serial.printf("[ACCPREP] %d,%d,%d\n", (int)v.kind, (int)v.pose, (int)v.attempt);
      return true;
    }

    bool sample(imu_cal::AccelRawSample& s) override {
      ImuSample ims;
      if (!w_.readSample_(ims)) return false;
      s.t_us = ims.sample_us;
      s.a = ims.a;
      s.w = ims.w;
      s.tempC = ims.tempC;
#if ATOMS3R_ICAL_RAW_LOG
      Serial.printf("[ACCRAW] %lu,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f\n", (unsigned long)s.t_us,
                    (double)s.a.x(), (double)s.a.y(), (double)s.a.z(),
                    (double)s.w.x(), (double)s.w.y(), (double)s.w.z(), (double)s.tempC);
#endif
      return true;
    }

    uint32_t nowMs() override { return millis(); }

    // Capture screen in the pose's display rotation; keeps the placement
    // guidance on screen and only redraws the hint row and the bar.
    void capture(const imu_cal::AccelStepView& v, const imu_cal::AccelHoldView& h) override {
      M5Ui& ui = w_.ui_;
      const uint32_t now = millis();
      if (!drawn_) {
        ui.setRotation(rot_add_(M5UiCfg::ROT_READ, v.rot_delta));
        ui.title(v.title);
        ui.line(v.label);
        for (int j = 0; j < 3; ++j) if (v.lines[j]) ui.line(v.lines[j]);
        ui.line("");
        hint_row_ = ui.cursorY();
        ui.line(h.hint_text);
        last_hint_ = h.hint;
        drawn_ = true;
        last_bar_ms_ = 0;
      }
      Input::update();
      if (h.hint != last_hint_) {
        ui.lineAt(hint_row_, h.hint_text);
        last_hint_ = h.hint;
      }
      if ((uint32_t)(now - last_bar_ms_) >= 100) {
        ui.bar01(h.progress01);
        last_bar_ms_ = now;
      }
    }

    void holdOk(const imu_cal::AccelStepView& v) override {
      w_.ui_.showOkAuto(v.label, "Captured");
      delay(80);
    }

    bool holdRetry(const imu_cal::AccelStepView& v, const char* why) override {
      return w_.ui_.retryMenu(v.label, why, "Hold still, same pose");
    }

    bool runFit(imu_cal::AccelFitJob& job, const char* what) override {
      w_.fit_.job = &job;
      const bool ok = w_.runFitTask_(FitKind::ACCEL_JOB, what, false);
      w_.fit_.job = nullptr;
      return ok;
    }

    void log(const char* line) override { Serial.println(line); }
    void idle() override { delay(2); }

  private:
    ImuCalWizard& w_;
    bool drawn_ = false;
    int hint_row_ = 0;
    imu_cal::AccelHoldHint last_hint_ = imu_cal::AccelHoldHint::PLACE;
    uint32_t last_bar_ms_ = 0;
  };

  bool accelFail_() {
    char l1[22], l2[22];
    accel_.failLines(l1, l2);
    Serial.printf("[ACC] failed: %s | %s (%s); previous calibration kept\n", l1, l2, accel_.failureDetail());
    ui_.failLines("ACCEL", l1, l2, "Previous cal kept");
    return false;
  }

  void sensorIdentity_() {
    const uint64_t mac = ESP.getEfuseMac();
    sensor_id_lo_ = (uint32_t)(mac & 0xFFFFFFFFu);
    sensor_id_hi_ = (uint32_t)(mac >> 32);
    imu_type_ = (uint8_t)M5.Imu.getType();
  }

  // Rebuilds RuntimeCals from `b` and re-validates the accelerometer set on
  // the retained observations (float, exactly as applied at runtime).
  bool validateStored_(const ImuCalBlobV3& b, const char* what) {
    RuntimeCals rc;
    rc.rebuildFromBlob(b);
    imu_cal::AccelFullFitResult r = accel_.result();
    const bool ok = rc.acc.ok && AccelProc::Fitter::validateFloat(accel_.obs(), accel_.nObs(), rc.acc, accel_fcfg_, r) &&
                    accelMetaBound(b);
    Serial.printf("[SAVE] %s float check=%d hold_rms double=%.5f float=%.5f\n", what, (int)ok,
                  r.ref_hold_rms, r.float_hold_rms);
    return ok;
  }

  void showDone_(const ImuCalBlobV3& rb) {
    ui_.setReadRotation();
    ui_.title("DONE");
    M5.Display.printf("A:%d G:%d M:%d\n", (int)rb.accel_ok, (int)rb.gyro_ok, (int)rb.mag_ok);
    ui_.line("Saved OK");
    char l[22];
    float bs = rb.accel_bias_sigma[0];
    if (rb.accel_bias_sigma[1] > bs) bs = rb.accel_bias_sigma[1];
    if (rb.accel_bias_sigma[2] > bs) bs = rb.accel_bias_sigma[2];
    snprintf(l, sizeof(l), "Bias sd %.4f m/s2", (double)bs);
    ui_.line(l);
    switch ((imu_cal::AccelThermal)rb.accel_thermal) {
      case imu_cal::AccelThermal::LEARNED: ui_.line("Temp slope: learned"); break;
      case imu_cal::AccelThermal::PRESERVED: ui_.line("Temp slope: kept"); break;
      default: ui_.line("Temp slope: none"); break;
    }
    ui_.line("");
    ui_.line("Tap BtnA");
    while (true) { Input::update(); if (Input::tapPressed()) break; delay(10); }
  }

  // MAG stage with retry loop. Returns true on success; false with redo_all
  // set when the user asked to restart, false otherwise on abort.
  bool runMagStage_(bool& redo_all) {
    redo_all = false;
    while (true) {
      magCal_.clear();

      const char* cap_why = nullptr;
      if (!captureMag_(cap_why)) {
        auto act = ui_.magFailMenu(cap_why ? cap_why : "Capture failed", "Flip + roll + pitch");
        if (act == M5Ui::MagFailAction::RETRY_MAG) continue;
        if (act == M5Ui::MagFailAction::REDO_ALL) { redo_all = true; return false; }
        return false;
      }

      if (magCal_.buf.n < ImuCalWizardCfg::MAG_MIN_TO_FIT) {
        auto act = ui_.magFailMenu("Too few accepted", "Rotate longer, slower");
        if (act == M5Ui::MagFailAction::RETRY_MAG) continue;
        if (act == M5Ui::MagFailAction::REDO_ALL) { redo_all = true; return false; }
        return false;
      }

      if (!runFitTask_(FitKind::MAG, "MAG", false)) {
        const char* why = imu_cal::fitFailStr(fit_.reason);
        const char* hint = "Try bigger 3D motion";

        if (why && (strstr(why, "NONPOSITIVE") || strstr(why, "NON POSITIVE") || strstr(why, "MODEL_S"))) {
          hint = "Too planar: flip all faces";
        }

        auto act = ui_.magFailMenu(why ? why : "Fit failed", hint);
        if (act == M5Ui::MagFailAction::RETRY_MAG) continue;
        if (act == M5Ui::MagFailAction::REDO_ALL) { redo_all = true; return false; }
        return false;
      }
      return true;  // MAG succeeded
    }
  }

private:
  // FIT task machinery
  enum class FitKind : uint8_t { ACCEL_JOB=0, GYRO=1, MAG=2 };

  struct FitCtx {
    volatile bool done = false;
    volatile bool ok   = false;
    imu_cal::FitFail reason = imu_cal::FitFail::BAD_ARG;
    FitKind kind = FitKind::GYRO;
    imu_cal::AccelFitJob* job = nullptr;

    ImuCalWizard* wiz = nullptr;
    TaskHandle_t  task = nullptr;
  } fit_{};

  static uint32_t hwmBytes_() {
    const uint32_t words = (uint32_t)uxTaskGetStackHighWaterMark(nullptr);
    return words * (uint32_t)sizeof(StackType_t);
  }

  // Accelerometer procedure fit (preliminary or final) on the large-stack task.
  static void fitTaskAccelJob_(void* p) {
    FitCtx* ctx = (FitCtx*)p;
    ctx->reason = imu_cal::FitFail::OK;
    if (ctx->job) ctx->job->run();
    ctx->ok = (ctx->job != nullptr);
    Serial.printf("[ACC] stack_hwm=%luB\n", (unsigned long)hwmBytes_());
    ctx->done = true;
    vTaskDelete(nullptr);
  }

  static void fitTaskGyro_(void* p) {
    FitCtx* ctx = (FitCtx*)p;
    ImuCalWizard* self = ctx->wiz;

    ctx->reason = imu_cal::FitFail::BAD_ARG;
    const bool ok = self->gyroCal_.fit(self->gyr_out_, &ctx->reason);

    Serial.printf("[GYR] fit=%d out.ok=%d reason=%s\n",
                  (int)ok, (int)self->gyr_out_.ok, imu_cal::fitFailStr(ctx->reason));

    ctx->ok = ok && self->gyr_out_.ok;
    Serial.printf("[GYR] stack_hwm=%luB\n", (unsigned long)hwmBytes_());
    ctx->done = true;
    vTaskDelete(nullptr);
  }

  static void fitTaskMag_(void* p) {
    FitCtx* ctx = (FitCtx*)p;
    ImuCalWizard* self = ctx->wiz;

    struct Try { int iters; float trim; float ridge; };
    const Try tries[] = {
      {3, 0.15f, 1e-6f},
      {3, 0.15f, 3e-6f},
      {3, 0.15f, 1e-5f},
      {3, 0.08f, 1e-6f},
      {2, 0.15f, 1e-6f},
    };

    imu_cal::FitFail last_reason = imu_cal::FitFail::BAD_ARG;
    bool any_ok = false;

    for (size_t i = 0; i < sizeof(tries)/sizeof(tries[0]); ++i) {
      imu_cal::FitFail r = imu_cal::FitFail::BAD_ARG;
      const bool ok = self->magCal_.fit(self->mag_out_, tries[i].iters, tries[i].trim, tries[i].ridge, &r);

      Serial.printf("[MAG] try iters=%d trim=%.3f ridge=%.1e -> fit=%d out.ok=%d reason=%s\n",
                    tries[i].iters, (double)tries[i].trim, (double)tries[i].ridge,
                    (int)ok, (int)self->mag_out_.ok, imu_cal::fitFailStr(r));

      last_reason = r;
      if (ok && self->mag_out_.ok) { any_ok = true; break; }
    }

    ctx->reason = last_reason;
    ctx->ok = any_ok;

    const uint32_t hw = hwmBytes_();
    Serial.printf("[MAG] stack_hwm=%luB\n", (unsigned long)hw);
    if (hw < 4096) Serial.printf("[MAG] WARN: low stack headroom: %luB\n", (unsigned long)hw);

    ctx->done = true;
    vTaskDelete(nullptr);
  }

  bool runFitTask_(FitKind kind, const char* what, bool show_fail) {
    fit_.done = false;
    fit_.ok = false;
    fit_.reason = imu_cal::FitFail::BAD_ARG;
    fit_.kind = kind;
    fit_.wiz = this;
    fit_.task = nullptr;

    ui_.setReadRotation();
    ui_.title("FIT");
    ui_.line(what);
    ui_.line("Working...");

    TaskFunction_t fn = nullptr;
    switch (kind) {
      case FitKind::ACCEL_JOB: fn = &ImuCalWizard::fitTaskAccelJob_; break;
      case FitKind::GYRO:  fn = &ImuCalWizard::fitTaskGyro_;  break;
      case FitKind::MAG:   fn = &ImuCalWizard::fitTaskMag_;   break;
      default:             fn = &ImuCalWizard::fitTaskGyro_;  break;
    }

    // Pin FIT to core 0 (avoid starving loopTask on core 1)
    const BaseType_t rc = xTaskCreatePinnedToCore(
      fn,
      what,
      (uint32_t)ImuCalWizardCfg::FIT_STACK_WORDS, // FreeRTOS (StackType_t)
      &fit_,
      1,
      &fit_.task,
      0
    );

    if (rc != pdPASS || fit_.task == nullptr) {
      if (show_fail) ui_.fail("FIT", "Task create failed");
      return false;
    }

    const uint32_t t0 = millis();
    float ph = 0.f;

    while (!fit_.done) {
      Input::update();

      if ((uint32_t)(millis() - t0) > ImuCalWizardCfg::FIT_TIMEOUT_MS) {
        vTaskDelete(fit_.task);
        fit_.task = nullptr;
        if (show_fail) ui_.fail("FIT", "Timeout");
        return false;
      }

      ph += 0.09f;
      ui_.bar01(0.5f + 0.5f * sinf(ph));
      delay(30);
    }

    fit_.task = nullptr;

    if (!fit_.ok) {
      if (show_fail) ui_.fail(what, imu_cal::fitFailStr(fit_.reason));
      return false;
    }
    return true;
  }

  bool magAvailable_() {
    // Actively probe MAG so we can skip calibration cleanly when
    // magnetometer data is unavailable/stuck.
    Vector3f m = Vector3f::Zero();
    bool have_prev = false;
    Vector3f prev = Vector3f::Zero();
    constexpr uint32_t probe_ms = 1200;
    const uint32_t t0 = millis();

    while ((uint32_t)(millis() - t0) < probe_ms) {
      if (!readMagSample_(m)) {
        delay(4);
        continue;
      }
      const float mn = m.norm();
      if (!(mn > 1.0f) || !std::isfinite(mn)) {
        delay(4);
        continue;
      }
      if (have_prev) {
        const float dm = (m - prev).norm();
        if (std::isfinite(dm) && dm >= ImuCalWizardCfg::MAG_MIN_DELTA_uT) {
          return true;
        }
      }
      prev = m;
      have_prev = true;
      delay(4);
    }
    return false;
  }

  // Capture steps
  bool readSample_(ImuSample& s) {
    return readImuMapped(M5.Imu, s);
  }

  bool readMagSample_(Vector3f& m_out) {
    (void)M5.Imu.update();
    const auto data = M5.Imu.getImuData();
    m_out = map_mag_to_body_uT_(data.mag);
    return finite3_(m_out);
  }

  void configureCalibrators_() {
    // Physical gravity at the calibration site: the static norm the
    // accelerometer is fitted to and the level the still gates compare with.
    // (ImuCalCfg::g_std is only the nominal-g -> m/s^2 unit.)
    const float g = ImuCalCfg::g_cal_local;

    gyroCal_.g  = g;

    // Gates
    gyroCal_.max_accel_dev = 0.8f;          // m/s^2
    gyroCal_.max_gyro_norm = 0.12f;         // rad/s

    // Accelerometer: full symmetric (PolarSPD) matrix with the existing
    // plausibility gates; see imu_cal::AccelFitCfg for the remaining gates.
    accel_ccfg_ = imu_cal::AccelCaptureCfg{};
    accel_ccfg_.g = g;
    accel_ccfg_.place_ms = ImuCalWizardCfg::PLACE_TIME_MS;
    accel_ccfg_.stuck_ms = ImuCalWizardCfg::STUCK_MS;
    accel_fcfg_ = imu_cal::AccelFitCfg{};
    accel_fcfg_.g = g;
    accel_fcfg_.diag_lo = 0.80;
    accel_fcfg_.diag_hi = 1.25;
    accel_fcfg_.max_cond = 6.0;
    accel_fcfg_.max_offdiag_rms = 0.10;
  }

  bool captureGyro_() {
    ui_.waitTap("GYRO", "SCREEN UP", "Tap then place");

    ui_.setReadRotation();
    ui_.title("GYRO");
    ui_.line("SCREEN UP");
    ui_.line("");
    ui_.line("Place on table");
    ui_.line("Do NOT touch");

    const uint32_t t0 = millis();
    while ((uint32_t)(millis() - t0) < ImuCalWizardCfg::PLACE_TIME_MS) {
      Input::update();
      ui_.bar01((float)(millis() - t0) / (float)ImuCalWizardCfg::PLACE_TIME_MS);
      delay(30);
    }

    const int start_n  = gyroCal_.buf.n;
    const int target_n = start_n + ImuCalWizardCfg::GYRO_NEED;

    ui_.title("GYRO");
    ui_.line("Capturing...");

    const uint32_t tcap0 = millis();
    uint32_t last_change = millis();
    int last_n = gyroCal_.buf.n;

    while ((uint32_t)(millis() - tcap0) < ImuCalWizardCfg::GYRO_TIMEOUT_MS) {
      Input::update();

      ImuSample s;
      if (!readSample_(s)) { delay(4); continue; }

      gyroCal_.addSample(s.w, s.a, s.tempC);

      const int n = gyroCal_.buf.n;
      if (n != last_n) { last_n = n; last_change = millis(); }

      if ((uint32_t)(millis() - last_change) > ImuCalWizardCfg::STUCK_MS) {
        Serial.printf("[GYR] stuck got=%d\n", n - start_n);
        ui_.fail("GYRO", "No samples accepted");
        return false;
      }

      ui_.bar01((float)(n - start_n) / (float)ImuCalWizardCfg::GYRO_NEED);

      if (n >= target_n) {
        ui_.showOkAuto("GYRO", "Captured");
        return true;
      }

      delay(5);
    }

    ui_.fail("GYRO", "Timeout");
    return false;
  }

  // MAG capture: downsample + min time + reject stale repeats + 3D coverage tests
  bool captureMag_(const char*& out_why) {
    out_why = nullptr;

    if (!magAvailable_()) {
      out_why = "MAG unavailable";
      return false;
    }

    ui_.waitTap("MAG", "Rotate ~45 sec", "Tap to start");

    const int start_n  = magCal_.buf.n;
    const int target_n = start_n + ImuCalWizardCfg::MAG_NEED;

    ui_.setReadRotation();
    ui_.title("MAG");
    ui_.line("Rotate now");
    ui_.line("Flip all faces");
    ui_.line("Avoid metal");

    const uint32_t tcap0 = millis();
    uint32_t last_change = millis();
    int last_n = magCal_.buf.n;

    uint32_t last_add_ms = 0;
    Vector3f last_added(NAN, NAN, NAN);

    // Raw bounds
    Vector3f vmin(+1e9f, +1e9f, +1e9f);
    Vector3f vmax(-1e9f, -1e9f, -1e9f);

    // Direction coverage using a stable-ish center (midpoint of bounds)
    Vector3f umin(+1e9f, +1e9f, +1e9f);
    Vector3f umax(-1e9f, -1e9f, -1e9f);
    Vector3f center = Vector3f::Zero();

    while ((uint32_t)(millis() - tcap0) < ImuCalWizardCfg::MAG_TIMEOUT_MS) {
      Input::update();

      Vector3f m;
      if (!readMagSample_(m)) {
        // During startup/recovery the mag path may briefly report NaN/Inf.
        // Treat this as "no usable sample yet" rather than a hard failure.
        delay(2);
        continue;
      }

      const uint32_t now = millis();

      // downsample accepted samples
      if (last_add_ms == 0 || (uint32_t)(now - last_add_ms) >= ImuCalWizardCfg::MAG_SAMPLE_SPACING_MS) {
        // reject stale repeats
        bool accept = true;
        if (isfinite(last_added.x())) {
          const Vector3f d = m - last_added;
          if (d.norm() < ImuCalWizardCfg::MAG_MIN_DELTA_uT) accept = false;
        }

        if (accept) {
          const int before = magCal_.buf.n;
          magCal_.addSample(m);
          const int after = magCal_.buf.n;

          if (after > before) {
            last_add_ms = now;
            last_added = m;

            vmin = vmin.cwiseMin(m);
            vmax = vmax.cwiseMax(m);

            center = 0.5f * (vmin + vmax);

            Vector3f c = m - center;
            const float cn = c.norm();
            if (cn > 1e-6f) {
              const Vector3f u = c / cn;
              umin = umin.cwiseMin(u);
              umax = umax.cwiseMax(u);
            }
          }
        }
      }

      const int n = magCal_.buf.n;
      if (n != last_n) { last_n = n; last_change = millis(); }

      if ((uint32_t)(millis() - last_change) > ImuCalWizardCfg::STUCK_MS) {
        Serial.printf("[MAG] stuck n=%d (no accepted samples)\n", n - start_n);
        out_why = "No MAG samples";
        return false;
      }

      const uint32_t elapsed = (uint32_t)(millis() - tcap0);

      // progress
      const float pS = (float)(n - start_n) / (float)ImuCalWizardCfg::MAG_NEED;
      const float pT = (float)elapsed / (float)ImuCalWizardCfg::MAG_MIN_TIME_MS;

      const Vector3f ur = umax - umin; // 0..2
      float pC = 0.f;
      if (isfinite(ur.x()) && isfinite(ur.y()) && isfinite(ur.z())) {
        const float px = clamp01_(ur.x() / ImuCalWizardCfg::MAG_URANGE_TARGET);
        const float py = clamp01_(ur.y() / ImuCalWizardCfg::MAG_URANGE_TARGET);
        const float pz = clamp01_(ur.z() / ImuCalWizardCfg::MAG_URANGE_TARGET);
        pC = second_smallest3_(px, py, pz);
      }

      float p = pS;
      if (pT < p) p = pT;
      if (pC < p) p = pC;
      ui_.bar01(p);

      // finish gate
      if (n >= target_n && elapsed >= ImuCalWizardCfg::MAG_MIN_TIME_MS) {
        const Vector3f span = vmax - vmin;

        float a = span.x(), b = span.y(), c = span.z();
        float smin = a, smid = b, smax = c;
        if (smin > smid) { float t=smin; smin=smid; smid=t; }
        if (smid > smax) { float t=smid; smid=smax; smax=t; }
        if (smin > smid) { float t=smin; smin=smid; smid=t; }

        const float rmin = (smax > 1e-6f) ? (smin / smax) : 0.f;
        const float rmid = (smax > 1e-6f) ? (smid / smax) : 0.f;

        Serial.printf("[MAG] n=%d elapsed=%.1fs span=(%.3f,%.3f,%.3f) ratios=(%.2f,%.2f) urange=(%.2f,%.2f,%.2f)\n",
                      n - start_n, (double)(elapsed / 1000.0f),
                      (double)span.x(), (double)span.y(), (double)span.z(),
                      (double)rmin, (double)rmid,
                      (double)ur.x(), (double)ur.y(), (double)ur.z());

        if (rmin < ImuCalWizardCfg::MAG_SPAN_MIN_FRAC ||
            rmid < ImuCalWizardCfg::MAG_SPAN_MID_FRAC) {
          out_why = "Span ratios too low";
          return false;
        }

        int ok_axes = 0;
        ok_axes += (ur.x() >= ImuCalWizardCfg::MAG_URANGE_TARGET) ? 1 : 0;
        ok_axes += (ur.y() >= ImuCalWizardCfg::MAG_URANGE_TARGET) ? 1 : 0;
        ok_axes += (ur.z() >= ImuCalWizardCfg::MAG_URANGE_TARGET) ? 1 : 0;
        if (ok_axes < 2) {
          out_why = "Direction range too small";
          return false;
        }

        if (smax < 1e-3f) { out_why = "MAG not changing"; return false; }

        const float detC = unit_dir_cov_det_(magCal_.buf.v, n);
        Serial.printf("[MAG] detC=%.6f\n", (double)detC);

        if (detC < 2.0e-4f) {
          out_why = "Coverage too flat";
          return false;
        }

        ui_.showOkAuto("MAG", "Captured");
        return true;
      }

      delay(5);
    }

    out_why = "Timeout";
    return false;
  }

  // Gyro and magnetometer fields of a full-wizard candidate (the accelerometer
  // set is filled by fillAccelFromFit()).
  void fillGyroMag_(ImuCalBlobV3& blob) {
    memset((void*)&blob, 0, sizeof(blob));
    blob.magic = ImuCalBlobV3::IMU_CAL_MAGIC;
    blob.version = ImuCalBlobV3::IMU_CAL_VERSION;
    blob.size_bytes = sizeof(ImuCalBlobV3);

    blob.gyro_ok = gyr_out_.ok ? 1 : 0;
    blob.gyro_T0 = gyr_out_.biasT.T0;
    blob.gyro_b0[0]=gyr_out_.biasT.b0.x(); blob.gyro_b0[1]=gyr_out_.biasT.b0.y(); blob.gyro_b0[2]=gyr_out_.biasT.b0.z();
    blob.gyro_k[0]=gyr_out_.biasT.k.x();   blob.gyro_k[1]=gyr_out_.biasT.k.y();   blob.gyro_k[2]=gyr_out_.biasT.k.z();

    blob.mag_ok = mag_out_.ok ? 1 : 0;
    mat_to_rowmajor9_(mag_out_.A, blob.mag_A);
    blob.mag_b[0]=mag_out_.b.x(); blob.mag_b[1]=mag_out_.b.y(); blob.mag_b[2]=mag_out_.b.z();
    blob.mag_field_uT = mag_out_.field_uT;
    blob.mag_rms      = mag_out_.rms;
  }

private:
  M5Ui& ui_;
  ImuCalStoreNvs& store_;

public:
  // Stored inside wizard => not on stack
  AccelProc accel_{};
  imu_cal::AccelCaptureCfg accel_ccfg_{};
  imu_cal::AccelFitCfg accel_fcfg_{};
  imu_cal::GyroCalibrator<float,  400, 8> gyroCal_{};
  imu_cal::MagCalibrator<float,   400>    magCal_{};

  imu_cal::GyroCalibration<float>  gyr_out_{};
  imu_cal::MagCalibration<float>   mag_out_{};

  uint32_t sensor_id_lo_ = 0, sensor_id_hi_ = 0;
  uint8_t imu_type_ = 0;
};

} // namespace atoms3r_ical
