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
  Nothing is saved unless the complete candidate validates. Failed writes
  attempt recovery and report whether the previous bytes were verified.
*/

#include <Arduino.h>
#include <M5Unified.h>

#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits>
#include <algorithm>
#include <atomic>

#include <freertos/FreeRTOS.h>
#include <freertos/task.h>

#include "AtomS3R/AtomS3R_ImuCal.h"         // ImuSample, axis mapping conventions, blob/store/runtime helpers
#include "AtomS3R/AtomS3R_M5Ui.h"           // UI + Input + clamp01_
#include "AtomS3R/AtomS3R_CalLog.h"         // best-effort pose diagnostics
#include "imu_calibrate/CalibrateIMU.h"     // imu_cal::* + FitFail
#include "imu_calibrate/AccelCalCapture.h"  // accelerometer procedure (host-tested)
#include "imu_calibrate/MagCalSampling.h"   // magnetic moments and capture
#include "imu_calibrate/GyroCalCapture.h"   // continuous quiet hold

// Set to 1 to stream every raw accel/gyro sample as [ACCRAW] lines for
// tests/imu_calibrate/accel_cal-replay. This diagnostic mode also preserves
// complete block logs and may wait for serial output; keep the host reading.
#ifndef ATOMS3R_ICAL_RAW_LOG
#define ATOMS3R_ICAL_RAW_LOG 0
#endif

namespace atoms3r_ical {

// Wizard configuration
struct ImuCalWizardCfg {
  // Step pacing (accelerometer holds: imu_cal::AccelCaptureCfg)
  static constexpr uint32_t PLACE_TIME_MS       = 6500;  // Gyro placement countdown only.
  static constexpr uint32_t GYRO_TIMEOUT_MS     = 70000;
  static constexpr uint32_t STUCK_MS            = 12000;

  // Accelerometer observation capacity (blocks) and hold capacity
  static constexpr int ACCEL_MAX_OBS            = 340;
  static constexpr int ACCEL_MAX_HOLDS          = 24;

  // Magnetic fit/verification use independent, motion-bounded sample means.
  static constexpr int MAG_MIN_TO_FIT           = 80;

  // Minimums before fitting
  static constexpr int GYRO_MIN_TO_FIT          = 120;

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

static inline bool finite3_(const Vector3f& v) {
  return isfinite(v.x()) && isfinite(v.y()) && isfinite(v.z());
}

class ImuCalWizard {
public:
  ImuCalWizard(M5Ui& ui, ImuCalStoreNvs& store)
  : ui_(ui), store_(store) {}

  // Runs the wizard and saves to NVS. Returns true only if saved successfully.
  // out_saved is filled with the saved blob (readback-validated). On any
  // abort before SAVE leaves storage untouched. Failed writes report whether
  // recovery was actually verified; storage failure can also prevent rollback.
  bool runAndSave(ImuCalBlobV4& out_saved) {
    if (runAndSave_(out_saved)) return true;
    // The magnetometer stage switches to the calibration source; a run that
    // did not save keeps the previous calibration, so read the magnetometer
    // the way that one was fitted again.
    ImuCalBlobV4 prev{};
    if (store_.load(prev) && magSetValid(prev) && magSourceFollower()) magSourceFollower()(prev.mag_source);
    return false;
  }

private:
  bool runAndSave_(ImuCalBlobV4& out_saved) {
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
      mag_verified_ = false;

      // Previous calibration: source of the preserved gyro/mag
      // (accelerometer-only mode), a compatible thermal slope, and the
      // stationary gyro level used by the rotation gate.
      ImuCalBlobV4 prev{};
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
      ImuCalBlobV4 blob;
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
      if (!validateStoredOnTask_(blob, "candidate")) {
        ui_.fail("SAVE", "Float check failed");
        return false;
      }

      ui_.setReadRotation();
      ui_.title("SAVE");
      ui_.line("Writing...");
      ImuCalBlobV4 rb{};
      const bool saved = store_.saveVerified(blob, rb);
      // ...and again after the read-back, which must be byte-identical to it.
      const bool rb_ok = saved && validateStoredOnTask_(rb, "readback");
      Serial.printf("[SAVE] verified=%d readback_valid=%d\n", (int)saved, (int)rb_ok);
      if (!saved || !rb_ok) {
        if (store_.lastSaveStatus() == ImuCalStoreNvs::SaveStatus::RECOVERY_FAILED) {
          Serial.println("[CAL] save failed; rollback/cleanup NOT verified; storage uncertain");
          ui_.fail("SAVE", "Recovery failed");
        } else {
          Serial.printf("[CAL] save failed; storage status=%u\n", (unsigned)store_.lastSaveStatus());
          ui_.fail("SAVE", "Write/readback fail");
        }
        return false;
      }

      out_saved = rb;
      showDone_(rb);
      Serial.println("[WIZ] done");
      return true;
    }
  }

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

    void log(const char* line) override {
#if ATOMS3R_ICAL_RAW_LOG
      Serial.println(line);
#else
      // A completed pose must not wait for a host to drain USB diagnostics.
      tryCalLogLine(Serial, line);
#endif
    }
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
  bool validateStored_(const ImuCalBlobV4& b, const char* what) {
    RuntimeCals rc;
    rc.rebuildFromBlob(b);
    imu_cal::AccelFullFitResult r = accel_.result();
    bool ok = rc.acc.ok && AccelProc::Fitter::validateFloat(accel_.obs(), accel_.nObs(), rc.acc, accel_fcfg_, r) &&
                    accelMetaBound(b);
    if (mag_verified_) {
      imu_cal::MagFitQuality q;
      ok = ok && rc.mag.ok && magCal_.check(rc.mag, q);
    }
    Serial.printf("[SAVE] %s float check=%d hold_rms double=%.5f float=%.5f\n", what, (int)ok,
                  r.ref_hold_rms, r.float_hold_rms);
    return ok;
  }

  void showDone_(const ImuCalBlobV4& rb) {
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
    // Calibrate in the same BMM150 setting and units every sketch runs with.
    // Re-applying is idempotent; on failure the driver values are used, as at
    // runtime, and the saved calibration records which source it was fitted on.
    configureAtomS3RMag(Serial);
    return runMagCaptureStage_(redo_all);
  }

  bool runMagCaptureStage_(bool& redo_all) {
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

      const bool fitted = runFitTask_(FitKind::MAG, "MAG", false);
      const char* why = fitted ? nullptr : magFailureText_();
      if (fitted && captureMag_(why, true)) {
        // Independent observations are never added to the fit.
        if (runFitTask_(FitKind::MAG_VERIFY, "MAG check", false)) {
          mag_verified_ = true;
          ui_.showOkAuto("MAG", "Verified");
          return true;
        }
        why = fit_.reason == imu_cal::FitFail::BAD_ARG ? "Check task failed" :
              imu_cal::magFitGateText(mag_verify_quality_.gate);
      }
      const auto act = ui_.magFailMenu(why ? why : "Check failed", "Turn slowly; retry");
      if (act == M5Ui::MagFailAction::RETRY_MAG) continue;
      if (act == M5Ui::MagFailAction::REDO_ALL) { redo_all = true; return false; }
      return false;
    }
  }

private:
  // FIT task machinery
  enum class FitKind : uint8_t { ACCEL_JOB=0, GYRO=1, MAG=2, MAG_VERIFY=3, SAVE_CHECK=4 };

  struct FitCtx {
    // Publish fitted coefficients across cores before the UI consumes them.
    std::atomic<bool> done{false};
    volatile bool ok   = false;
    imu_cal::FitFail reason = imu_cal::FitFail::BAD_ARG;
    FitKind kind = FitKind::GYRO;
    imu_cal::AccelFitJob* job = nullptr;
    const ImuCalBlobV4* blob = nullptr;   // SAVE_CHECK
    const char* blob_what = "";

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

  // Pre-save check of a stored blob. It includes the magnetometer quality
  // check, which needs more stack than the Arduino loop task provides.
  static void fitTaskSaveCheck_(void* p) {
    FitCtx* ctx = (FitCtx*)p;
    ctx->reason = imu_cal::FitFail::OK;
    ctx->ok = ctx->blob && ctx->wiz->validateStored_(*ctx->blob, ctx->blob_what);
    Serial.printf("[SAVE] stack_hwm=%luB\n", (unsigned long)hwmBytes_());
    ctx->done = true;
    vTaskDelete(nullptr);
  }

  bool validateStoredOnTask_(const ImuCalBlobV4& b, const char* what) {
    fit_.blob = &b;
    fit_.blob_what = what;
    const bool ok = runFitTask_(FitKind::SAVE_CHECK, "SAVE check", false);
    fit_.blob = nullptr;
    return ok;
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

  const char* magFailureText_() const {
    if (fit_.reason == imu_cal::FitFail::BAD_ARG) return "Fit task failed";
    return fit_.reason == imu_cal::FitFail::MAG_QUALITY_FAIL ?
        imu_cal::magFitGateText(magCal_.quality.gate) : "Need more 3D motion";
  }

  static void fitTaskMag_(void* p) {
    FitCtx* ctx = (FitCtx*)p;
    ImuCalWizard* self = ctx->wiz;
    ctx->ok = self->magCal_.fit(self->mag_out_, 3, 0.15f, 1e-6f, &ctx->reason);
    const auto& q = self->magCal_.quality;
    Serial.printf("[MAG] fit=%d gate=%s rms=%.4f p95=%.4f inliers=%d/%d cells=%d bias_sd=%.4f "
                  "matrix_sd=%.5f drift=%.4f refine=%d cost=%.8g->%.8g stack=%luB\n",
                  (int)ctx->ok, imu_cal::magFitGateText(q.gate), q.rms, q.p95, q.inliers, q.samples,
                  q.cells, q.max_bias_sigma, q.max_matrix_sigma, q.time_drift, q.iterations,
                  q.initial_cost, q.refined_cost, (unsigned long)hwmBytes_());
    ctx->done = true;
    vTaskDelete(nullptr);
  }

  static void fitTaskMagVerify_(void* p) {
    FitCtx* ctx = (FitCtx*)p;
    ImuCalWizard* self = ctx->wiz;
    const auto& m = self->mag_out_;
    auto& q = self->mag_verify_quality_;
    q = imu_cal::MagFitQuality{};
    ctx->ok = self->magCal_.check(m, q);
    ctx->reason = ctx->ok ? imu_cal::FitFail::OK : imu_cal::FitFail::MAG_VERIFY_FAIL;
    Serial.printf("[MAG CHECK] ok=%d gate=%s rms=%.4f p95=%.4f inliers=%d/%d cells=%d drift=%.4f\n",
                  (int)ctx->ok, imu_cal::magFitGateText(q.gate), q.rms, q.p95, q.inliers, q.samples, q.cells, q.time_drift);
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
    ui_.title(kind == FitKind::MAG ? "REFINE" : "FIT");
    ui_.line(what);
    ui_.line(kind == FitKind::MAG ? "Improving the fit" : "Working...");

    TaskFunction_t fn = nullptr;
    switch (kind) {
      case FitKind::ACCEL_JOB: fn = &ImuCalWizard::fitTaskAccelJob_; break;
      case FitKind::GYRO:  fn = &ImuCalWizard::fitTaskGyro_;  break;
      case FitKind::MAG:   fn = &ImuCalWizard::fitTaskMag_;   break;
      case FitKind::MAG_VERIFY: fn = &ImuCalWizard::fitTaskMagVerify_; break;
      case FitKind::SAVE_CHECK: fn = &ImuCalWizard::fitTaskSaveCheck_; break;
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
        fit_.reason = imu_cal::FitFail::BAD_ARG;
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
    m_out = readMagBody_(data.mag);
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
    accel_ccfg_.stuck_ms = ImuCalWizardCfg::STUCK_MS;
    accel_fcfg_ = imu_cal::AccelFitCfg{};
    accel_fcfg_.g = g;
    accel_fcfg_.diag_lo = 0.80;
    accel_fcfg_.diag_hi = 1.25;
    accel_fcfg_.max_cond = 6.0;
    accel_fcfg_.max_offdiag_rms = 0.10;
  }

  bool captureGyro_() {
    ui_.waitTap("GYRO", "Place on table", "Tap then place");

    ui_.setReadRotation();
    ui_.title("GYRO");
    ui_.line("Let it rest");
    ui_.line("");
    ui_.line("Place on table");
    ui_.line("Leave it still");

    const uint32_t t0 = millis();
    while ((uint32_t)(millis() - t0) < ImuCalWizardCfg::PLACE_TIME_MS) {
      Input::update();
      ui_.bar01((float)(millis() - t0) / (float)ImuCalWizardCfg::PLACE_TIME_MS);
      delay(30);
    }

    ui_.title("GYRO");
    ui_.line("Keep on the table");
    const int hint_row = ui_.cursorY();
    ui_.line("Settling...");
    ui_.line("About 8 quiet sec");
    imu_cal::GyroCaptureCfg cfg;
    cfg.timeout_ms = ImuCalWizardCfg::GYRO_TIMEOUT_MS;
    cfg.stuck_ms = ImuCalWizardCfg::STUCK_MS;
    imu_cal::GyroCapture<float,400,8> capture(gyroCal_, cfg);
    capture.begin(millis());
    uint32_t last_draw = millis();
    uint32_t last_resets = 0;
    while (true) {
      Input::update();
      ImuSample s;
      const bool valid = readSample_(s);
      const uint32_t now = millis();
      const auto status = capture.update(now, valid ? &s.a : nullptr, valid ? &s.w : nullptr,
                                         valid ? s.tempC : NAN, valid ? &s.m : nullptr);
      if (capture.resets() != last_resets) {
        last_resets = capture.resets();
        Serial.printf("[GYR] restart=%lu reason=%s\n", (unsigned long)last_resets, capture.resetReason());
      }
      if (uint32_t(now-last_draw) >= 100) {
        ui_.lineAt(hint_row, status == imu_cal::GyroCaptureStatus::MOVING ? "Moving - hold still" :
            (status == imu_cal::GyroCaptureStatus::SETTLING ? "Settling..." : "Hold still"));
        ui_.bar01(capture.progress());
        last_draw = now;
      }
      if (status == imu_cal::GyroCaptureStatus::READY) {
        Serial.printf("[GYR] quiet capture n=%d restarts=%lu\n", gyroCal_.buf.n, (unsigned long)capture.resets());
        ui_.showOkAuto("GYRO", "Captured");
        return true;
      }
      if (status == imu_cal::GyroCaptureStatus::STALE || status == imu_cal::GyroCaptureStatus::TIMEOUT) {
        if (ui_.retryMenu("GYRO", status == imu_cal::GyroCaptureStatus::STALE ? "No fresh samples" : "Stillness not reached",
                          "Keep on the table")) {
          capture.begin(millis());
          last_resets = 0;
          ui_.title("GYRO"); ui_.line("Keep on the table"); ui_.line("Settling..."); ui_.line("About 8 quiet sec");
          continue;
        }
        return false;
      }
      delay(2);
    }
  }

  // MAG capture: distinct-reading moments + 3D coverage.
  bool captureMag_(const char*& out_why, bool verify = false) {
    out_why = nullptr;

    if (!verify && !magAvailable_()) {
      out_why = "MAG unavailable";
      return false;
    }

    ui_.waitTap(verify ? "CHECK MAG" : "MAG", verify ? "Turn and tilt again" : "Turn and tilt",
                "About 1-2 min");

    ui_.setReadRotation();
    ui_.title(verify ? "CHECK MAG" : "MAG");
    ui_.line("Avoid metal");
    const int hint_row = ui_.cursorY();
    ui_.line("Turn slowly");
    const int state_row = ui_.cursorY();
    ui_.line(verify ? "Checking new motion" : "Collecting");
    ui_.line("No exact angles");

    auto cfg = imu_cal::MagSampleWindow::captureCfg(verify);
    cfg.stuck_ms = ImuCalWizardCfg::STUCK_MS;
    cfg.min_delta_uT = ImuCalWizardCfg::MAG_MIN_DELTA_uT;
    cfg.span_min_frac = ImuCalWizardCfg::MAG_SPAN_MIN_FRAC;
    cfg.span_mid_frac = ImuCalWizardCfg::MAG_SPAN_MID_FRAC;
    cfg.urange_target = ImuCalWizardCfg::MAG_URANGE_TARGET;
    imu_cal::MagCapture<float,400> capture(magCal_, cfg, verify ? &mag_out_ : nullptr);
    const uint32_t tcap0 = millis();
    capture.begin(tcap0);
    imu_cal::MagSampleWindow window;window.begin();
    uint32_t last_draw = tcap0;
    uint32_t last_log = tcap0;
    while (true) {
      Input::update();
      Vector3f m, mean;
      Eigen::Matrix3f covariance;
      const bool valid = readMagSample_(m);
      const uint32_t now = millis();
      const bool averaged=window.update(now,valid?&m:nullptr,mean,covariance);
      const auto status = capture.update(now,averaged?&mean:nullptr,valid?&m:nullptr,
                                          averaged?&covariance:nullptr);
      if (uint32_t(now - last_draw) >= 250) {
        ui_.bar01(capture.progress(now));
        ui_.lineAt(hint_row, capture.hint(now));
        char state[32];
        snprintf(state,sizeof(state),"Samples %d/%d",magCal_.buf.n,cfg.required_samples);
        ui_.lineAt(state_row,status==imu_cal::MagCaptureStatus::COVERAGE_LOW?"More directions":state);
        last_draw = now;
      }
      if (status == imu_cal::MagCaptureStatus::READY || status == imu_cal::MagCaptureStatus::TIMEOUT ||
          uint32_t(now-last_log)>=5000) {
        const auto c=capture.coverage();char line[192];
        snprintf(line,sizeof(line),"[MAG] n=%d/%d observations=%lu elapsed=%.1fs ratios=(%.2f,%.2f) detC=%.6f cells=%d",
                 magCal_.buf.n,cfg.required_samples,(unsigned long)capture.observations(),double(uint32_t(now-tcap0))/1000,
                 double(c.min_ratio),double(c.mid_ratio),double(c.determinant),c.cells);
        tryCalLogLine(Serial,line);
        last_log=now;
      }
      if (status == imu_cal::MagCaptureStatus::READY || status == imu_cal::MagCaptureStatus::COVERAGE_LOW) {
        if (status == imu_cal::MagCaptureStatus::COVERAGE_LOW) {
          // Keep guiding until coverage is sufficient or the bounded timeout.
          delay(5);
          continue;
        }
        if (!verify) ui_.showOkAuto("MAG", "Captured");
        return true;
      }
      if (status == imu_cal::MagCaptureStatus::STALE) { out_why = "No MAG samples"; return false; }
      if (status == imu_cal::MagCaptureStatus::TIMEOUT) break;
      if (status == imu_cal::MagCaptureStatus::BAD_CONFIG) { out_why = "Bad capture config"; return false; }
      delay(5);
    }

    out_why = magCal_.buf.n<cfg.required_samples ? "Need more MAG data" : "Need more directions";
    return false;
  }

  // Gyro and magnetometer fields of a full-wizard candidate (the accelerometer
  // set is filled by fillAccelFromFit()).
  void fillGyroMag_(ImuCalBlobV4& blob) {
    memset((void*)&blob, 0, sizeof(blob));
    blob.magic = ImuCalBlobV4::IMU_CAL_MAGIC;
    blob.version = ImuCalBlobV4::IMU_CAL_VERSION;
    blob.size_bytes = sizeof(ImuCalBlobV4);

    fillGyroFromFit(blob, gyr_out_);

    blob.mag_ok = mag_out_.ok ? 1 : 0;
    blob.mag_source = activeMagSource();
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
  imu_cal::MagFitQuality mag_verify_quality_{};
  bool mag_verified_ = false;

  uint32_t sensor_id_lo_ = 0, sensor_id_hi_ = 0;
  uint8_t imu_type_ = 0;
};

} // namespace atoms3r_ical
