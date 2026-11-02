#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R IMU calibration wizard UI (Accel + Gyro + Mag) using imu_cal::* calibrators.

  Uses M5Unified IMU API as the IMU sample source.

  Accelerometer: ten guided static holds (six faces + four screen-up corners),
  bounded targeted extra holds, three later rechecks, and a full-matrix fit
  with thermal-slope qualification (imu_cal::AccelCalProcedure). In the full
  wizard the rechecks follow the gyro and magnetometer stages; the
  separate noise mode retains the current intrinsic calibration.
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
#include <memory>
#include <new>
#include <esp_heap_caps.h>
#include <esp_task_wdt.h>

#include <freertos/FreeRTOS.h>
#include <freertos/task.h>

#include "AtomS3R/AtomS3R_ImuCal.h"         // ImuSample, axis mapping conventions, blob/store/runtime helpers
#include "AtomS3R/AtomS3R_M5Ui.h"           // UI + Input + clamp01_
#include "AtomS3R/AtomS3R_CalLog.h"         // best-effort pose diagnostics
#include "imu_calibrate/CalibrateIMU.h"     // imu_cal::* + FitFail
#include "imu_calibrate/AccelCalCapture.h"  // accelerometer procedure (host-tested)
#include "imu_calibrate/MagCalSampling.h"   // magnetic moments and capture
#include "imu_calibrate/MagTimingCal.h"
#include "imu_calibrate/GyroCalFit.h"
#include "imu_calibrate/GyroTurnCapture.h"
#include "imu_calibrate/FrameAlignment.h"
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
#include "imu_calibrate/GyroCalCapture.h"   // continuous quiet hold

// Set to 1 to stream every raw accel/gyro sample as [ACCRAW] lines for
// tests/imu_calibrate/accel_cal-replay. This diagnostic mode also preserves
// complete block logs and may wait for serial output; keep the host reading.
#ifndef ATOMS3R_ICAL_RAW_LOG
#define ATOMS3R_ICAL_RAW_LOG 0
#endif
#ifndef ATOMS3R_NOISE_RAW_LOG
#define ATOMS3R_NOISE_RAW_LOG 0
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
  static constexpr uint32_t WIZARD_STACK_BYTES  = 32768;

  static constexpr uint32_t FIT_TIMEOUT_MS      = 30000;
  static constexpr uint32_t GYRO_MATRIX_TIMEOUT_MS = 120000;
  static constexpr uint32_t NOISE_TIMEOUT_MS = 2520000;
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
  bool runAndSave(ImuCalBlobV6& out_saved) {
    // The wizard owns heap storage, but its nested UI/capture/validation calls
    // still need more stack than Arduino's loopTask can provide at boot.
    // Keep the existing workflow on the caller's core, with its own stack.
    WizardTaskCtx ctx{this,&out_saved};
    TaskHandle_t task=nullptr;
    const BaseType_t rc=xTaskCreatePinnedToCore(&ImuCalWizard::wizardTask_,"imu-cal-ui",
      ImuCalWizardCfg::WIZARD_STACK_BYTES,&ctx,1,&task,xPortGetCoreID());
    if(rc!=pdPASS || task==nullptr) {
      Serial.println("[WIZ] UI task allocation failed");
      ui_.fail("CALIBRATION","Task create failed");return false;
    }
    while(!ctx.done.load()) {
      // If an application enabled the loop watchdog, waiting for user input
      // or the 40-minute capture must not stop the caller's watchdog service.
      if(esp_task_wdt_status(nullptr)==ESP_OK)esp_task_wdt_reset();
      delay(20);
    }
    return ctx.saved;
  }

private:
  struct WizardTaskCtx {
    ImuCalWizard* wizard;
    ImuCalBlobV6* out;
    bool saved=false;
    std::atomic<bool> done{false};
  };
  static void wizardTask_(void* p) {
    auto* ctx=static_cast<WizardTaskCtx*>(p);
    Serial.printf("[WIZ] UI task stack=%luB free=%luB\n",
      static_cast<unsigned long>(ImuCalWizardCfg::WIZARD_STACK_BYTES),
      static_cast<unsigned long>(hwmBytes_()));
    ctx->saved=ctx->wizard->runAndRestore_(*ctx->out);
    Serial.printf("[WIZ] UI task done saved=%d stack_free=%luB\n",int(ctx->saved),
      static_cast<unsigned long>(hwmBytes_()));
    ctx->done=true; // no access to ctx or wizard after publishing completion
    vTaskDelete(nullptr);
  }

  bool runAndRestore_(ImuCalBlobV6& out_saved) {
    if (runAndSave_(out_saved)) return true;
    // The magnetometer stage switches to the calibration source; a run that
    // did not save keeps the previous calibration, so read the magnetometer
    // the way that one was fitted again.
    ImuCalBlobV6 prev{};
    activeMagTiming() = {};
    if (store_.load(prev)) {
      if (magSetValid(prev) && magSourceFollower()) magSourceFollower()(prev.mag_source);
      if (magTimingLoader()) magTimingLoader()(prev);
    }
    return false;
  }

  /**
   * Run calibration stages and retry menus. Fill out_saved and return true only after
   * saving and validating readback; return false on cancellation or save failure.
   */
  bool runAndSave_(ImuCalBlobV6& out_saved) {
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
      timing_fit_ = {}; timing_saved_ = false;

      // A valid CURRENT-schema record may supply a qualified thermal prior
      // and the raw stationary gyro level used by the accelerometer gate.
      ImuCalBlobV6 prev{};
      const bool have_prev = store_.load(prev);
      sensorIdentity_();
      // Gyro endpoints and magnetic alignment depend on the final accel fit.
      // A full recalibration is required when changing that intrinsic model.
      RuntimeCals previous_runtime;
      if(have_prev)previous_runtime.rebuildFromBlob(prev);
      const bool can_noise=have_prev && !calibrationUnusableReason(true,previous_runtime);

      const M5Ui::StartAction start = ui_.startMenu(can_noise);
      if (start == M5Ui::StartAction::CANCEL) {
        Serial.println("[WIZ] cancelled; previous calibration kept");
        return false;
      }
      if(start==M5Ui::StartAction::NOISE)return captureNoise_(prev,out_saved);
      Serial.printf("[WIZ] mode=full have_prev=%d\n", (int)have_prev);
      // Reserve the trajectory before asking for any poses. A wrong board
      // memory setting must not discard a completed accel/mag acquisition.
      if(!prepareGyroRotationStorage_())return false;

      const imu_cal::AccelThermalPrior prior =
          have_prev ? accelThermalPriorFrom(prev, sensor_id_lo_, sensor_id_hi_, imu_type_) : imu_cal::AccelThermalPrior{};
      Serial.printf("[ACC] thermal prior: %s\n", prior.valid ? "compatible slope available" : "none");
      const Vector3f gyro_level = have_prev && prev.gyro_ok
          ? Vector3f(prev.gyro_b0[0], prev.gyro_b0[1], prev.gyro_b0[2]) : Vector3f::Zero();
      // Replay context (tests/imu_calibrate/accel_cal-replay).
      Serial.printf("[ACCMODE] full g=%.7f\n", (double)accel_fcfg_.g);
      Serial.printf("[ACCPRIOR] %d,%.7f,%.7f,%.7f,%.2f,%.2f,%.2f,%.2f\n", (int)prior.valid, prior.k[0], prior.k[1],
                    prior.k[2], prior.k_temp_lo, prior.k_temp_hi, prior.clamp_lo, prior.clamp_hi);
      Serial.printf("[ACCGYRO] %.7f,%.7f,%.7f,%d\n", (double)gyro_level.x(), (double)gyro_level.y(),
                    (double)gyro_level.z(), (int)(have_prev && prev.gyro_ok));

      accel_.begin(accel_ccfg_, accel_fcfg_, prior, gyro_level, have_prev && prev.gyro_ok);
      AccelIo io(*this);

      if (!accel_.runMainStage(io)) return accelFail_();

      if (!captureGyro_()) return false;
      if (gyroCal_.buf.n < ImuCalWizardCfg::GYRO_MIN_TO_FIT) {
        ui_.fail("GYRO", "Too few accepted");
        return false;
      }
      if (!runFitTask_(FitKind::GYRO, "GYRO", true)) return false;
      accel_.setGyroReference(gyr_out_.biasT.b0);
      Serial.printf("[ACCGYRO] %.7f,%.7f,%.7f,1\n", (double)gyr_out_.biasT.b0.x(), (double)gyr_out_.biasT.b0.y(),
                    (double)gyr_out_.biasT.b0.z());

      // Complete calibration requires a magnetic triad and its alignment;
      // retry the probe or abort without writing.
      bool mag_present = magAvailable_();
      while (!mag_present) {
        Serial.println("[MAG] not responding; retry or abort full calibration");
        const auto act = ui_.magFailMenu("MAG not responding", "Abort keeps old cal");
        if (act == M5Ui::MagFailAction::REDO_ALL) { redo_all = true; break; }
        if (act != M5Ui::MagFailAction::RETRY_MAG) return false;
        mag_present = magAvailable_();
      }
      if (redo_all) {
        Serial.println("[WIZ] redo all requested");
        continue;
      }
      if (!runMagStage_(redo_all)) {
        if (redo_all) {
          Serial.println("[WIZ] redo all requested");
          continue;
        }
        return false;
      }

      // Later rechecks (thermal evidence + held-out verification), final fit.
      if (!accel_.runRecheckStage(io)) return accelFail_();
      if (!accel_.runFinalFit(io)) return accelFail_();
      for(;;) {
        if(!captureGyroRotations_())return false;
        if(runFitTask_(FitKind::GYRO_MATRIX,"GYRO matrix",false))break;
        if(!ui_.retryMenu("GYRO / ALIGN","More varied tilts","Repeat rotations"))return false;
      }
      gyro_rotations_.reset();
      if(!captureMounting_())return false;
      // Optional magnetic timing record. It is stored for diagnostics only;
      // runtime magnetic corrections do not apply it.
      if (mag_out_.ok && gyr_out_.ok && activeMagTimingConfiguration()) {
        timing_cal_.reset();
        while (!captureMagTiming_()) {
          const auto action=ui_.timingFailMenu("Timing not identified", imu_cal::magTimingReason(timing_fit_.reason));
          if (action==M5Ui::TimingFailAction::RETRY) continue;
          if (action==M5Ui::TimingFailAction::SKIP) { timing_fit_ = {}; break; }
          return false;
        }
      }

      // One candidate contains every intrinsic fit and frame qualification.
      imu_cal::AccelCalibration<float> fc;
      AccelProc::Fitter::toFloat(accel_.result(), accel_fcfg_.g, fc);
      const uint32_t capture_s = accel_.totalHoldMs() / 1000u;
      ImuCalBlobV6 blob;
      fillGyroMag_(blob);
      blob.mag_alignment_ok=alignment_.ok;
      mat_to_rowmajor9_(alignment_.R,blob.R_IM);
      mat_to_rowmajor9_(mount_HI_,blob.R_HI);mat_to_rowmajor9_(mount_VH_,blob.R_VH);
      blob.housing_alignment_ok=blob.vessel_alignment_ok=1;
      blob.alignment_train_sd=alignment_.dip_sd_after;blob.alignment_validation_sd=alignment_validation_;
      blob.alignment_condition=alignment_.condition;blob.alignment_train_poses=16;blob.alignment_validation_poses=5;
      fillAccelFromFit(blob, accel_.result(), fc, capture_s);
      blob.sensor_id_lo = sensor_id_lo_;
      blob.sensor_id_hi = sensor_id_hi_;
      blob.imu_type = imu_type_;
      blob.mag_timing = {};
      if (blob.mag_ok && timing_fit_.ok) {
        auto& t=blob.mag_timing;
        t.version=1; t.binding=magTimingBinding(blob); t.configuration=activeMagTimingConfiguration();
        t.delay_us=timing_fit_.delay_us; t.fit_width_s=timing_fit_.fit_width_s;
        t.residual_time_s=timing_fit_.residual_time_s;
        t.train_rms_rad=timing_fit_.train_rms_rad; t.validation_rms_rad=timing_fit_.validation_rms_rad;
        t.improvement=timing_fit_.improvement;
        t.train_pairs=timing_fit_.train_pairs; t.validation_pairs=timing_fit_.validation_pairs;
        t.crc=magTimingCrc(t);
      }

      // The stored float set, rebuilt through the runtime path, must reproduce
      // the fit before it is written...
      blob=ImuCalStoreNvs::sealed_(blob);
      if (!validateStoredOnTask_(blob, "candidate")) {
        ui_.fail("SAVE", "Float check failed");
        return false;
      }

      ui_.setReadRotation();
      ui_.title("SAVE");
      ui_.line("Writing...");
      ImuCalBlobV6 rb{};
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

      timing_saved_ = rb.mag_ok && magTimingValid(rb.mag_timing);
      activeMagTiming() = {};
      if (magTimingLoader()) magTimingLoader()(rb);
      out_saved = rb;
      showDone_(rb);
      Serial.println("[WIZ] done");
      return true;
    }
  }

  using AccelProc = imu_cal::AccelCalProcedure<ImuCalWizardCfg::ACCEL_MAX_OBS, ImuCalWizardCfg::ACCEL_MAX_HOLDS>;
  // Include room for ordinary pauses across all twenty turns, not just the
  // ideal uninterrupted motion. This remains an explicitly reserved PSRAM buffer.
  using RotationFit=imu_cal::GyroRotationFit<131072,24>;
  struct RotationDeleter {
    void operator()(RotationFit* p) const {if(p){p->~RotationFit();heap_caps_free(p);}}
  };

  bool prepareGyroRotationStorage_() {
    if(gyro_rotations_)return true;
    constexpr unsigned caps=MALLOC_CAP_SPIRAM|MALLOC_CAP_8BIT;
    const size_t total=heap_caps_get_total_size(caps);
    const size_t available=heap_caps_get_free_size(caps);
    const size_t largest=heap_caps_get_largest_free_block(caps);
#if defined(CONFIG_SPIRAM_MODE_OCT) && CONFIG_SPIRAM_MODE_OCT
    const char* mode="OPI";
#elif defined(CONFIG_SPIRAM_MODE_QUAD) && CONFIG_SPIRAM_MODE_QUAD
    const char* mode="QSPI";
#else
    const char* mode="unavailable";
#endif
    Serial.printf("[GYRROT] PSRAM found=%u mode=%s total=%u free=%u largest=%u need=%u\n",
      unsigned(psramFound()),mode,unsigned(total),unsigned(available),unsigned(largest),unsigned(sizeof(RotationFit)));
    // The registered byte-addressable external heap is authoritative, even
    // if a custom SDK initialized it without Arduino's psramFound flag.
    if(total==0) {
      Serial.println("[GYRROT] Enable Tools > PSRAM > OPI PSRAM, then rebuild and upload");
      ui_.failLines("GYRO memory","PSRAM unavailable","Enable OPI PSRAM","Rebuild & upload");
      return false;
    }
    if(largest<sizeof(RotationFit)) {
      ui_.failLines("GYRO memory",available>=sizeof(RotationFit)?"PSRAM fragmented":"PSRAM too small",
        "Restart then retry");return false;
    }
    void* memory=heap_caps_malloc(sizeof(RotationFit),caps);
    if(!memory) {
      Serial.println("[GYRROT] PSRAM allocation failed after preflight");
      ui_.failLines("GYRO memory","Allocation failed","Restart then retry");return false;
    }
    gyro_rotations_.reset(new(memory) RotationFit{});
    return true;
  }

  static constexpr uint8_t turn_pose_ids_[4]={0,2,4,5};
  struct TurnCtx {
    ImuCalWizard* wiz=nullptr;
    Vector3f targets[4];int segment=0;
    std::atomic<imu_cal::GyroTurnPhase> phase{imu_cal::GyroTurnPhase::REST};
    std::atomic<int> target{-1},progress{0},rotation{0};
    std::atomic<const char*> reason{"none"};
    std::atomic<uint32_t> resets{0};
    std::atomic<bool> done{false},retry{false};
    bool ok=false;
  } turn_;

  // Same pose rotations as AccelIo; use gravity in its measured device frame.
  static int gyroReadRotation_(const Vector3f& a,const Vector3f (&targets)[4],int previous) {
    if(!a.allFinite() || a.norm()<1)return previous;
    const Vector3f up=a.normalized();
    const float screen=up.dot(targets[0]),usb=up.dot(targets[1]),left=up.dot(targets[2]);
    if(std::fabs(screen)>.8f)return 0;
    if(std::fabs(usb)>.8f)return usb>0?2:0;
    if(std::fabs(left)>.8f)return left>0?1:3;
    return previous; // no redraw/flicker between edge orientations
  }

  static void gyroTurnTask_(void* context) {
    auto& c=*static_cast<TurnCtx*>(context);auto& w=*c.wiz;
    imu_cal::GyroTurnCapture<> capture(*w.gyro_rotations_,w.gyroCal_,w.gyr_out_,c.targets,c.segment);
    capture.begin(millis());uint32_t next=micros(),log=millis(),rotate_since=millis();
    int candidate=c.rotation.load();
    for(;;) {
      // Only this worker accesses the IMU. LCD/BtnA work runs on the UI task;
      // a slow redraw must never create an unobserved integration interval.
      const bool retry=c.retry.load();
      if(!retry && int32_t(micros()-next)<0){vTaskDelay(1);continue;}
      next=micros()+5000;
      ImuSample sample;const bool fresh=!retry && w.readSample_(sample);
      const auto phase=retry?capture.cancel():capture.update(millis(),fresh?sample.sample_us:0,
        fresh?&sample.a:nullptr,fresh?&sample.w:nullptr,
        fresh?sample.tempC:NAN,fresh?&sample.m:nullptr);
      if(fresh) {
        const int rotation=gyroReadRotation_(sample.a,c.targets,c.rotation.load());
        if(rotation!=candidate){candidate=rotation;rotate_since=millis();}
        if(uint32_t(millis()-rotate_since)>=100)c.rotation=candidate;
      }
      c.target=capture.target();c.progress=int(capture.progress()*1000);
      c.reason=capture.reason();c.resets=capture.resets();c.phase=phase;
      if(phase==imu_cal::GyroTurnPhase::DONE || phase==imu_cal::GyroTurnPhase::FAILED) {
        c.ok=phase==imu_cal::GyroTurnPhase::DONE;
        if(c.ok) {
          w.endpoint_a_=capture.accel();w.endpoint_w_=capture.gyro();
          w.endpoint_m_=capture.mag();w.endpoint_temp_=capture.temperature();
        }
        break;
      }
      if(uint32_t(millis()-log)>=2000) {
        char line[160];snprintf(line,sizeof(line),
          "[GYRROT] turn=%d phase=%u target=%d hold=%d resets=%lu reason=%s",
          c.segment+1,unsigned(phase),capture.target(),c.progress.load(),
          static_cast<unsigned long>(capture.resets()),capture.reason());
        tryCalLogLine(Serial,line);log=millis();
      }
    }
    c.done=true; // publishes endpoint writes; no ctx access after this store
    vTaskDelete(nullptr);
  }

  static const char* gyroRestHelp_(const char* reason) {
    if(!std::strcmp(reason,"magnetic change"))return "Move away from metal";
    if(!std::strcmp(reason,"magnetic stream lost"))return "Mag sensor not ready";
    if(!std::strcmp(reason,"sample gap") || !std::strcmp(reason,"too few samples") ||
       !std::strcmp(reason,"frozen IMU"))return "Sensor data paused";
    if(!std::strcmp(reason,"gyro scatter") || !std::strcmp(reason,"accel scatter") ||
       !std::strcmp(reason,"gyro mean changed") || !std::strcmp(reason,"accel mean changed") ||
       !std::strcmp(reason,"invalid sample"))return "Movement detected";
    return "Keep steady";
  }

  void drawGyroTurn_(int segment,imu_cal::GyroTurnPhase phase,int target,
                    int progress,const char* reason,bool new_screen) {
    using Phase=imu_cal::GyroTurnPhase;
    const uint8_t pose=turn_pose_ids_[target>=0 && target<4?target:0];
    if(new_screen) {
      ui_.title(phase==Phase::TWIST?"TWIST":phase==Phase::TURN?"TIP":"HOLD");
      char count[24];snprintf(count,sizeof(count),"%s %d/20",
        phase==Phase::REST?"Start":phase==Phase::HOLD?"Finish":"Turn",segment+1);ui_.line(count);
      ui_.line("");
      if(phase==Phase::TWIST) {
        ui_.line("Keep this side up.");
        ui_.line(segment%2?"Turn right 1/4 turn.":"Turn left 1/4 turn.");
      } else if(phase==Phase::TURN || phase==Phase::HOLD) {
        ui_.line(pose==0?"Screen facing up":pose==2?"USB end up":
                 pose==4?"Left edge down":"Right edge down");
        ui_.line(phase==Phase::HOLD?"Hand support is OK.":"Then steady on table.");
      } else {
        ui_.line("Steady on table.");ui_.line("Hand support is OK.");
      }
      ui_.lineAt(100,"Tap: retry turn");
    }
    if(phase==Phase::REST || phase==Phase::HOLD) {
      char timer[32];snprintf(timer,sizeof(timer),"Steady %.1f / ~3.5s",progress*(imu_cal::GyroTurnCapture<>::HOLD_MS*1e-6f));
      ui_.lineAt(76,timer);ui_.lineAt(88,progress?"Measuring...":gyroRestHelp_(reason));
      ui_.bar01(progress*.001f);
    }
  }

  bool runGyroTurn_(int segment,const char*& why) {
    using Phase=imu_cal::GyroTurnPhase;
    char label[24];snprintf(label,sizeof(label),"Turn %d of 20",segment+1);
    ui_.waitTapHere("GYRO TILT",label,nullptr);
    turn_.wiz=this;turn_.segment=segment;turn_.phase=Phase::REST;
    turn_.target=-1;turn_.progress=0;turn_.reason="none";turn_.resets=0;
    turn_.rotation=(ui_.rotation()-M5UiCfg::ROT_READ)&3;turn_.retry=false;
    turn_.ok=false;turn_.done=false;
    for(int i=0;i<4;++i)if(!accel_.expectedUp(turn_pose_ids_[i],turn_.targets[i])) {
      why="Missing device frame";return false;
    }
    TaskHandle_t task=nullptr;
    if(xTaskCreatePinnedToCore(&ImuCalWizard::gyroTurnTask_,"gyro capture",
        ImuCalWizardCfg::FIT_STACK_BYTES,&turn_,2,&task,0)!=pdPASS || !task) {
      why="Capture task failed";return false;
    }
    auto shown=Phase::FAILED;int shown_target=-2,shown_rotation=-1;uint32_t draw=millis()-150;
    while(!turn_.done.load()) {
      Input::update();
      if(Input::tapPressed())turn_.retry=true;
      const auto phase=turn_.phase.load();const int target=turn_.target.load(),rotation=turn_.rotation.load();
      if(phase==Phase::DONE || phase==Phase::FAILED){delay(5);continue;}
      if(phase!=shown || target!=shown_target || rotation!=shown_rotation || uint32_t(millis()-draw)>=150) {
        if(rotation!=shown_rotation)ui_.setRotation(uint8_t((M5UiCfg::ROT_READ+rotation)&3));
        drawGyroTurn_(segment,phase,target,turn_.progress.load(),turn_.reason.load(),
                      phase!=shown || target!=shown_target || rotation!=shown_rotation);
        shown=phase;shown_target=target;shown_rotation=rotation;draw=millis();
      }
      delay(20);
    }
    why=turn_.reason.load();
    if(turn_.ok)ui_.showOkAuto("Turn kept",label);
    return turn_.ok;
  }

  bool captureGyroRotations_() {
    if(!prepareGyroRotationStorage_())return false;
    gyro_rotations_->clear();
    if(!captureGyro_() || !runFitTask_(FitKind::GYRO,"GYRO bias",true))return false;
    alignment_pose_count_=0;
    auto retainPose=[&]() {
      if(alignment_pose_count_>=24 || !endpoint_m_.allFinite())return false;
      raw_alignment_[alignment_pose_count_].a=endpoint_a_;
      raw_alignment_[alignment_pose_count_].m_raw=endpoint_m_;
      alignment_temp_[alignment_pose_count_++]=endpoint_temp_;return true;
    };
    if(!retainPose()){ui_.fail("ALIGN","No quiet mag data");return false;}
    Vector3f bias_sum=gyr_out_.biasT.b0;unsigned bias_poses=1;
    for(int segment=0;segment<20;++segment) {
      for(;;) {
        const char* why=nullptr;
        if(runGyroTurn_(segment,why))break;
        char retry[32],kept[32];snprintf(retry,sizeof(retry),"Retry turn %d / 20",segment+1);
        snprintf(kept,sizeof(kept),"%d turns kept",gyro_rotations_->segments());
        if(!ui_.retryMenuHere(retry,why,kept))return false;
      }
      if(!retainPose())return false;
      bias_sum+=endpoint_w_-gyr_out_.biasT.k*(endpoint_temp_-gyr_out_.biasT.T0);++bias_poses;
    }
    gyr_out_.biasT.b0=bias_sum/float(bias_poses);return true;
  }

  static void fitTaskGyroMatrix_(void* p) {
    FitCtx* ctx=static_cast<FitCtx*>(p);ImuCalWizard& w=*ctx->wiz;
    imu_cal::AccelCalibration<float> acc;AccelProc::Fitter::toFloat(w.accel_.result(),w.accel_fcfg_.g,acc);
    const uint32_t start=micros();
    ctx->ok=w.gyro_rotations_ && w.gyro_rotations_->fit(acc,w.gyr_out_,[](void*){vTaskDelay(1);});
    if(ctx->ok) {
      atoms3r_magdiag::Pose train[24],validation[24];int nt=0,nv=0;
      for(int i=0;i<w.alignment_pose_count_;++i) {
        auto pose=w.raw_alignment_[i];pose.a=acc.apply(pose.a,w.alignment_temp_[i]);
        const float norm=w.mag_out_.apply(pose.m_raw).norm();
        if(!std::isfinite(norm) || std::fabs(norm-w.mag_out_.field_uT)>3.f){ctx->ok=false;break;}
        if(i && i%4==0)validation[nv++]=pose;else train[nt++]=pose;
      }
      atoms3r_magdiag::MagModel model;model.A=w.mag_out_.A;model.b=w.mag_out_.b;
      if(ctx->ok) {
        w.alignment_=atoms3r_magdiag::solveAlignment(train,nt,model);
        ctx->ok=atoms3r_magdiag::validateAlignment(w.alignment_,validation,nv,model,w.alignment_validation_);
      }
    }
    ctx->reason=ctx->ok?imu_cal::FitFail::OK:imu_cal::FitFail::DEGENERATE_COVERAGE;
    Serial.printf("[GYRROT] ok=%d reason=%s fit_us=%lu train=%.6g validation=%.6g cond=%.6g stack_hwm=%luB\n",
      int(ctx->ok),w.gyro_rotations_?imu_cal::gyroRotationReason(w.gyro_rotations_->reason):"no buffer",
      static_cast<unsigned long>(micros()-start),double(w.gyr_out_.rotation_train_rms),
      double(w.gyr_out_.rotation_validation_rms),double(w.gyr_out_.rotation_condition),
      static_cast<unsigned long>(hwmBytes_()));
    Serial.printf("[ALIGN] ok=%d train_sd=%.4f independent_rms=%.4f cond=%.5g\n",int(w.alignment_.ok),
      double(w.alignment_.dip_sd_after),double(w.alignment_validation_),double(w.alignment_.condition));
    ctx->done=true;vTaskDelete(nullptr);
  }

  bool captureMounting_() {
    ui_.title("MOUNT");ui_.line("x fwd, y stbd, z down");
    ui_.line("Tap ONLY if I=H=V");ui_.line("Else serial: MOUNT");
    ui_.line("18 matrix entries");
    Serial.println("[MOUNT] Tap confirms aligned IMU/housing/vessel axes. Otherwise send MOUNT followed by R_HI and R_VH, each 9 row-major entries.");
    mount_HI_.setIdentity();mount_VH_.setIdentity();
    char line[400];size_t n=0;
    while(true) {
      Input::update();if(Input::tapPressed())return true;
      while(Serial.available()) {
        const char c=char(Serial.read());
        if(c=='\r')continue;
        if(c!='\n') {if(n+1<sizeof(line))line[n++]=c;else n=0;continue;}
        line[n]=0;n=0;
        float r[18];int end=0;
        const int fields=sscanf(line,"MOUNT %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %n",
          &r[0],&r[1],&r[2],&r[3],&r[4],&r[5],&r[6],&r[7],&r[8],&r[9],&r[10],&r[11],&r[12],
          &r[13],&r[14],&r[15],&r[16],&r[17],&end);
        if(fields==18 && line[end]==0 && rotationValid(r) && rotationValid(r+9)) {
          mount_HI_=mat_from_rowmajor9_(r);mount_VH_=mat_from_rowmajor9_(r+9);
          Serial.println("[MOUNT] accepted");return true;
        }
        Serial.println("[MOUNT] invalid: need two proper orthonormal rotations");
      }
      delay(2);
    }
  }

  static void fitTaskNoise_(void* p) {
    FitCtx* ctx=static_cast<FitCtx*>(p);ImuCalWizard& w=*ctx->wiz;
    RuntimeCals runtime;runtime.rebuildFromBlob(*ctx->blob);
    ctx->ok=false;ctx->reason=imu_cal::FitFail::BAD_ARG;
    if(!calibrationUnusableReason(true,runtime) && w.noise_capture_) {
      imu_cal::GyroCapture<float,400,8> quiet(w.gyroCal_);quiet.begin(millis());
      bool settled=false;uint32_t next=micros();
      while(!w.noise_cancel_.load()) {
        const int32_t wait=int32_t(next-micros());if(wait>0)delayMicroseconds(uint32_t(wait));
        next+=5000; // same 200 Hz clock as the production sensor path
        ImuSample s;const bool fresh=readImuMapped(M5.Imu,s);
        if(!settled) {
          const auto state=quiet.update(millis(),fresh?&s.a:nullptr,fresh?&s.w:nullptr,
            fresh?s.tempC:NAN,fresh?&s.m:nullptr);
          if(state==imu_cal::GyroCaptureStatus::READY) {
            settled=true;w.noise_capture_->begin(s.sample_us,.005,ImuCalCfg::g_cal_local);
          } else if(state==imu_cal::GyroCaptureStatus::STALE || state==imu_cal::GyroCaptureStatus::TIMEOUT) {
            ctx->reason=imu_cal::FitFail::GYRO_NOT_STILL;break;
          }
        } else if(fresh) {
          const Vector3f a=runtime.applyAccel(s.a,s.tempC),rate=runtime.applyGyro(s.w,s.tempC),m=runtime.applyMag(s.m);
          const auto state=w.noise_capture_->add(s.sample_us,a,rate,s.tempC,&m);
#if ATOMS3R_NOISE_RAW_LOG
          char raw[220];snprintf(raw,sizeof(raw),"[NOISERAW] %lu,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.6g",
            static_cast<unsigned long>(s.sample_us),double(a.x()),double(a.y()),double(a.z()),
            double(rate.x()),double(rate.y()),double(rate.z()),double(s.tempC));
          if(!tryCalLogLine(Serial,raw)){Serial.println("[NOISE] raw log dropped");break;}
#endif
          w.noise_progress_.store(w.noise_capture_->progress());
          if(state==imu_cal::NoiseCaptureStatus::READY){ctx->ok=true;ctx->reason=imu_cal::FitFail::OK;break;}
          if(state!=imu_cal::NoiseCaptureStatus::CAPTURING) {
            Serial.printf("[NOISE] rejected status=%u samples=%lu\n",unsigned(state),
              static_cast<unsigned long>(w.noise_capture_->samples()));break;
          }
        }
        // Yield before the next precise epoch; display work stays on core 1.
        vTaskDelay(1);
      }
    }
    Serial.printf("[NOISE] ok=%d cancelled=%d stack_free=%luB\n",int(ctx->ok),int(w.noise_cancel_.load()),
      static_cast<unsigned long>(hwmBytes_()));
    ctx->done=true;vTaskDelete(nullptr);
  }

  bool captureNoise_(const ImuCalBlobV6& previous,ImuCalBlobV6& out) {
    ui_.waitTap("NOISE 40M","Warm up; quiet table","Tap then leave still");
    Serial.printf("[NOISESTART] schema=6 calibration_crc=%08lx sensor=%08lx%08lx period=.005\n",
      static_cast<unsigned long>(previous.crc),static_cast<unsigned long>(previous.sensor_id_hi),
      static_cast<unsigned long>(previous.sensor_id_lo));
    noise_capture_.reset(new(std::nothrow) imu_cal::StationaryNoiseCapture<>{});
    if(!noise_capture_){ui_.fail("NOISE","Allocation failed");return false;}
    noise_progress_=0;noise_cancel_=false;fit_.blob=&previous;
    const bool acquired=runFitTask_(FitKind::NOISE_CAPTURE,"NOISE 40M",false);fit_.blob=nullptr;
    if(!acquired){if(!noise_cancel_.load())ui_.fail("NOISE","Motion/timing/temp");noise_capture_.reset();return false;}
    auto& capture=*noise_capture_;
    ImuCalBlobV6 candidate=previous;
    if(!capture.report(candidate.noise)){ui_.fail("NOISE","Noise not identified");return false;}
    for(int l=0;l<17;++l) {
      Serial.printf("[ALLAN] %.8g,%lu",capture.tau(l),static_cast<unsigned long>(capture.pairs(l)));
      for(int k=0;k<6;++k)Serial.printf(",%.9g",capture.deviation(l,k));Serial.println();
    }
    for(int k=0;k<6;++k)Serial.printf("[NOISE] axis=%d N=%.9g sample_std=%.9g BI=%s:%.9g RW=%s:%.9g\n",
      k,double(candidate.noise.white_density[k]),double(candidate.noise.sample_std[k]),
      candidate.noise.bias_mask&(1<<k)?"identified":"unidentified",double(candidate.noise.bias_instability[k]),
      candidate.noise.walk_mask&(1<<k)?"identified":"unidentified",double(candidate.noise.bias_random_walk[k]));
    if(!store_.saveVerified(candidate,out)) {
      Serial.printf("[NOISE] save failed status=%u\n",unsigned(store_.lastSaveStatus()));
      ui_.fail("NOISE",store_.lastSaveStatus()==ImuCalStoreNvs::SaveStatus::RECOVERY_FAILED?
        "Recovery failed":"Write/readback fail");return false;
    }
    noise_capture_.reset();ui_.showOkAuto("NOISE","Model saved");return true;
  }

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
  bool validateStored_(const ImuCalBlobV6& b, const char* what) {
    RuntimeCals rc;
    rc.rebuildFromBlob(b);
    imu_cal::AccelFullFitResult r = accel_.result();
    bool ok = !calibrationUnusableReason(true,rc) && rc.acc.ok && AccelProc::Fitter::validateFloat(accel_.obs(), accel_.nObs(), rc.acc, accel_fcfg_, r) &&
                    accelMetaBound(b);
    if (mag_verified_) {
      imu_cal::MagFitQuality q;
      ok = ok && rc.mag.ok && magCal_.check(rc.mag, q);
    }
    Serial.printf("[SAVE] %s float check=%d hold_rms double=%.5f float=%.5f\n", what, (int)ok,
                  r.ref_hold_rms, r.float_hold_rms);
    return ok;
  }

  void showDone_(const ImuCalBlobV6& rb) {
    ui_.setReadRotation();
    ui_.title("DONE");
    M5.Display.printf("A:%d G:%d M:%d\n", (int)rb.accel_ok, (int)rb.gyro_ok, (int)rb.mag_ok);
    ui_.line("Saved OK");
    ui_.line(timing_saved_ ? "MAG timing: saved" : "MAG timing: not saved");
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
  enum class FitKind : uint8_t { ACCEL_JOB=0, GYRO=1, MAG=2, MAG_VERIFY=3, SAVE_CHECK=4, GYRO_MATRIX=5, NOISE_CAPTURE=6 };

  struct FitCtx {
    // Publish fitted coefficients across cores before the UI consumes them.
    std::atomic<bool> done{false};
    volatile bool ok   = false;
    imu_cal::FitFail reason = imu_cal::FitFail::BAD_ARG;
    FitKind kind = FitKind::GYRO;
    imu_cal::AccelFitJob* job = nullptr;
    const ImuCalBlobV6* blob = nullptr;   // SAVE_CHECK
    const char* blob_what = "";

    ImuCalWizard* wiz = nullptr;
    TaskHandle_t  task = nullptr;
  } fit_{};

  static uint32_t hwmBytes_() {
    // ESP-IDF reports free stack in BYTES, unlike upstream FreeRTOS.
    return (uint32_t)uxTaskGetStackHighWaterMark(nullptr);
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

  bool validateStoredOnTask_(const ImuCalBlobV6& b, const char* what) {
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
    ui_.title(kind==FitKind::NOISE_CAPTURE?"NOISE 40M":kind == FitKind::MAG ? "REFINE" : "FIT");
    ui_.line(what);
    ui_.line(kind==FitKind::NOISE_CAPTURE?"Keep still; tap stops":kind == FitKind::MAG ? "Improving the fit" : "Working...");

    TaskFunction_t fn = nullptr;
    switch (kind) {
      case FitKind::ACCEL_JOB: fn = &ImuCalWizard::fitTaskAccelJob_; break;
      case FitKind::GYRO:  fn = &ImuCalWizard::fitTaskGyro_;  break;
      case FitKind::GYRO_MATRIX: fn=&ImuCalWizard::fitTaskGyroMatrix_;break;
      case FitKind::NOISE_CAPTURE: fn=&ImuCalWizard::fitTaskNoise_;break;
      case FitKind::MAG:   fn = &ImuCalWizard::fitTaskMag_;   break;
      case FitKind::MAG_VERIFY: fn = &ImuCalWizard::fitTaskMagVerify_; break;
      case FitKind::SAVE_CHECK: fn = &ImuCalWizard::fitTaskSaveCheck_; break;
      default:             fn = &ImuCalWizard::fitTaskGyro_;  break;
    }

    // Pin FIT to core 0 (avoid starving loopTask on core 1)
    const BaseType_t rc = xTaskCreatePinnedToCore(
      fn,
      what,
      (uint32_t)ImuCalWizardCfg::FIT_STACK_BYTES, // ESP-IDF uses bytes
      &fit_,
      kind==FitKind::NOISE_CAPTURE?2:1,
      &fit_.task,
      0
    );

    if (rc != pdPASS || fit_.task == nullptr) {
      if (show_fail) ui_.fail("FIT", "Task create failed");
      return false;
    }

    const uint32_t t0 = millis();
    float ph = 0.f;
    uint32_t last_draw=t0;
    const uint32_t timeout=kind==FitKind::NOISE_CAPTURE?ImuCalWizardCfg::NOISE_TIMEOUT_MS:
      kind==FitKind::GYRO_MATRIX?ImuCalWizardCfg::GYRO_MATRIX_TIMEOUT_MS:ImuCalWizardCfg::FIT_TIMEOUT_MS;

    while (!fit_.done) {
      Input::update();

      if ((uint32_t)(millis() - t0) > timeout) {
        vTaskDelete(fit_.task);
        fit_.task = nullptr;
        fit_.reason = imu_cal::FitFail::BAD_ARG;
        if (show_fail) ui_.fail("FIT", "Timeout");
        return false;
      }

      if(kind==FitKind::NOISE_CAPTURE) {
        if(Input::tapPressed())noise_cancel_=true;
        if(uint32_t(millis()-last_draw)>=1000){ui_.bar01(noise_progress_.load());last_draw=millis();}
      } else {ph += 0.09f;ui_.bar01(0.5f + 0.5f * sinf(ph));}
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
        endpoint_a_=capture.meanAccel();endpoint_w_=capture.meanGyro();endpoint_temp_=capture.meanTemperature();
        endpoint_m_=capture.meanMag();
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
        if (verify) {
          // The window statistics can hide field changes inside a window;
          // score the fresh single readings with the frozen candidate too.
          const auto raw = capture.rawResidual();
          char line[96];
          snprintf(line, sizeof(line), "[MAG RAW] n=%d rms=%.3f tail=%d ok=%d", raw.n, raw.rms, raw.tail, (int)raw.ok);
          tryCalLogLine(Serial, line);
          if (!raw.ok) { out_why = "Field inconsistent"; return false; }
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

  /**
   * Capture and validate magnetic delay using the current spatial/gyro fits.
   * Each segment collects for 14..45 seconds unless cancelled; retries retain
   * a successful training segment. Update timing_fit_ and display diagnostics
   * after capture. Return false for missing configuration, cancellation or a
   * failed fit; true means both independent captures passed the final fit.
   */
  bool captureMagTiming_() {
    if (!activeMagTimingConfiguration()) return false;
    timing_fit_ = {};
    const unsigned first = timing_cal_.segmentQuality(0).ok ? 1 : 0;
    for (unsigned segment=first; segment<2; ++segment) {
      timing_cal_.resetSegment(segment);
      ui_.waitTap(segment ? "CHECK TIME" : "MAG TIME", "Rock back and forth", "One rock per second");
      ui_.title("MAG TIME");
      ui_.line("Away from metal");
      ui_.line("Tilt left / right");
      ui_.line("Pause, then reverse");
      ui_.line("Tap to cancel");
      const uint32_t start=millis();
      uint32_t draw=start;
      // Keep collecting usable evidence after pauses or dropped frames. The
      // first segment freezes before the independent validation capture.
      while (uint32_t(millis()-start)<45000) {
        Input::update();
        if (Input::tapPressed()) return false;
        ImuSample sample;
        // Timing uses every gyro delivery, including gyro-only frames.
        if (readImuMapped(M5.Imu, sample, true)) {
          const Vector3f w=gyr_out_.apply(sample.w,sample.tempC);
          const Vector3f m=alignment_.R*mag_out_.apply(sample.m);
          timing_cal_.add(sample.sample_us,w,sample.mag_timing.have,
              sample.mag_timing.sequence,sample.mag_timing.frame_us,m,segment);
        }
        if (uint32_t(millis()-draw)>=500) {
          draw=millis();
          const auto quality=timing_cal_.segmentQuality(segment);
          char progress[32];
          snprintf(progress,sizeof(progress),"Pairs %u/50  %lus",timing_cal_.pairs(segment),
                   static_cast<unsigned long>((draw-start)/1000));
          ui_.lineAt(82,progress);
          ui_.lineAt(94,imu_cal::magTimingReason(quality.reason));
          ui_.bar01(float(timing_cal_.pairs(segment))/50.f);
          if (uint32_t(draw-start)>=14000 && quality.ok) break;
        }
        delay(1);
      }
      const auto quality=timing_cal_.segmentQuality(segment);
      if (!quality.ok) {
        timing_fit_=quality;
        break; // Show the actual failure without demanding a second capture.
      }
    }
    if (timing_cal_.segmentQuality(0).ok && timing_cal_.segmentQuality(1).ok)
      timing_fit_=timing_cal_.fit();
    Serial.printf("[MAGTIME FIT] ok=%u reason=%s delay_ms=%.1f width_ms=%.1f allowance_ms=%.1f "
                  "train_rad=%.5f validation_rad=%.5f improvement=%.3f pairs=%lu/%lu gaps=%u history=%u\n",
        unsigned(timing_fit_.ok),imu_cal::magTimingReason(timing_fit_.reason),double(timing_fit_.delay_us)*.001,
        double(timing_fit_.fit_width_s)*1000,double(timing_fit_.residual_time_s)*1000,
        double(timing_fit_.train_rms_rad),double(timing_fit_.validation_rms_rad),
        double(timing_fit_.improvement),static_cast<unsigned long>(timing_cal_.pairs(0)),
        static_cast<unsigned long>(timing_cal_.pairs(1)),timing_cal_.gaps(),timing_cal_.missingHistory());
    ui_.waitTap(timing_fit_.ok ? "TIME OK" : "TIME FAIL",
        imu_cal::magTimingReason(timing_fit_.reason),
        timing_fit_.ok ? "Will save with cal" : "Retry TIME stage");
    return timing_fit_.ok;
  }

  // Gyro and magnetometer fields of a full-wizard candidate (the accelerometer
  // set is filled by fillAccelFromFit()).
  void fillGyroMag_(ImuCalBlobV6& blob) {
    blob = ImuCalBlobV6{};
    blob.magic = ImuCalBlobV6::IMU_CAL_MAGIC;
    blob.version = ImuCalBlobV6::IMU_CAL_VERSION;
    blob.size_bytes = sizeof(ImuCalBlobV6);

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
  std::unique_ptr<RotationFit,RotationDeleter> gyro_rotations_;
  std::unique_ptr<imu_cal::StationaryNoiseCapture<>> noise_capture_;
  std::atomic<float> noise_progress_{0};std::atomic<bool> noise_cancel_{false};
  Vector3f endpoint_a_=Vector3f::Zero(),endpoint_w_=Vector3f::Zero(),endpoint_m_=Vector3f::Zero();
  float endpoint_temp_=25;
  atoms3r_magdiag::Pose raw_alignment_[24];
  float alignment_temp_[24]{};
  int alignment_pose_count_=0;
  atoms3r_magdiag::Alignment alignment_{};
  float alignment_validation_=INFINITY;
  Matrix3f mount_HI_=Matrix3f::Identity(),mount_VH_=Matrix3f::Identity();
  imu_cal::MagCalibrator<float,   400>    magCal_{};

  imu_cal::GyroCalibration<float>  gyr_out_{};
  imu_cal::MagCalibration<float>   mag_out_{};
  imu_cal::MagFitQuality mag_verify_quality_{};
  bool mag_verified_ = false;
  imu_cal::MagTimingCal timing_cal_{}; // fixed storage, never on the loop stack
  imu_cal::MagTimingFit timing_fit_{};
  bool timing_saved_ = false;

  uint32_t sensor_id_lo_ = 0, sensor_id_hi_ = 0;
  uint8_t imu_type_ = 0;
};

} // namespace atoms3r_ical
