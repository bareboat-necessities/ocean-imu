#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Manual accelerometer calibration procedure: pose design, static-hold
  qualification, targeted extra holds, later rechecks, and the fit sequence.

  Host-testable: no Arduino dependency. The AtomS3R wizard supplies screens,
  samples and the fit task through AccelCalIo; the host tests supply a
  simulated user and IMU through the same interface, so both run this code.

  Frames. Samples are in the project body frame (AtomS3R_ImuCal.h mapping).
  Pose geometry is written in the device frame
    top   = toward the end opposite the USB connector
    right = toward the right edge (screen read with USB at the bottom)
    out   = out of the screen
  and every "up" below is the direction of the measured specific force of a
  device at rest. The documented body mapping is top = +x, right = +y,
  out = -z. Only the screen-up sign is taken from it; all other pose regions
  are checked against directions measured earlier in the same session, so an
  error in the assumed x/y edge mapping cannot lock the user out. The measured
  frame is logged for hardware confirmation of the mapping.

  Hold qualification. The 6.5 s placement interval is guidance, not a wait:
  once the device has settled (consecutive quiet, consistent, in-region blocks)
  qualified blocks are retained immediately. A hold completes after
  hold_useful_ms of retained blocks (not after a sample count). Each block is a
  short-window mean (block_ms), kept as its own gravity observation so a
  slowly drifting hold keeps its directional spread; blocks with motion,
  rotation, shocks, norm jumps, stale data or gaps are rejected locally.
*/

#include "imu_calibrate/AccelCalFit.h"

#include <stdio.h>
#include <string.h>

namespace imu_cal {

// Pose table

enum class AccelHoldKind : uint8_t { MAIN = 0, EXTRA = 1, RECHECK = 2 };

enum class AccelRegion : uint8_t {
  SCREEN_UP_DOC = 0,   // cone about the documented screen-up direction
  OPPOSITE = 1,        // cone about minus a direction measured in pose `ref`
  VERTICAL_AXIS = 2,   // screen vertical and near an in-plane device axis
  RIGHT_HANDED = 3,    // cone about up(ref0) x up(ref1)
  FRAME = 4,           // cone about (top,right,out) in the measured frame
};

struct AccelPoseDef {
  const char* name;          // <= 15 chars so "10/10 " + name fits 21 columns
  const char* lines[3];      // instruction lines, <= 21 chars, nullptr-terminated
  int8_t top, right, out;    // up direction sign pattern in the device frame
  uint8_t rot_delta;         // display rotation relative to the reading rotation
  AccelRegion region;
  int8_t ref0, ref1;         // pose indices used by OPPOSITE / RIGHT_HANDED
};

static constexpr int kAccelMainPoses = 10;
static constexpr int kAccelRecheckPoses = 3;

// The first six are the guided faces. The four tilts are the screen-up
// corners: with the faces they make all six independent symmetric-matrix
// terms observable while keeping the screen
// readable during capture (see doc/imu_calibrate for the information study).
static const AccelPoseDef kAccelPoses[kAccelMainPoses] = {
  {"SCREEN UP",      {"Screen faces up", "Flat on the table", nullptr},     0,  0,  1, 0, AccelRegion::SCREEN_UP_DOC, -1, -1},
  {"SCREEN DOWN",    {"Screen faces table", "Flat on the table", nullptr},  0,  0, -1, 0, AccelRegion::OPPOSITE, 0, -1},
  {"USB UP",         {"USB points up", "Screen vertical", nullptr},        -1,  0,  0, 2, AccelRegion::VERTICAL_AXIS, 0, -1},
  {"USB DOWN",       {"USB points down", "Screen vertical", nullptr},       1,  0,  0, 0, AccelRegion::OPPOSITE, 2, -1},
  {"LEFT DOWN",      {"Left edge down", "Screen vertical", nullptr},        0,  1,  0, 1, AccelRegion::RIGHT_HANDED, 0, 2},
  {"RIGHT DOWN",     {"Right edge down", "Screen vertical", nullptr},       0, -1,  0, 3, AccelRegion::OPPOSITE, 4, -1},
  {"USB DN+LEFT DN", {"Screen up, ~45 deg", "USB end DOWN", "Left edge DOWN"},   1,  1,  1, 0, AccelRegion::FRAME, -1, -1},
  {"USB DN+RIGHT DN",{"Screen up, ~45 deg", "USB end DOWN", "Right edge DOWN"},  1, -1,  1, 0, AccelRegion::FRAME, -1, -1},
  {"USB UP+LEFT DN", {"Screen up, ~45 deg", "USB end UP", "Left edge DOWN"},    -1,  1,  1, 2, AccelRegion::FRAME, -1, -1},
  {"USB UP+RIGHT DN",{"Screen up, ~45 deg", "USB end UP", "Right edge DOWN"},   -1, -1,  1, 2, AccelRegion::FRAME, -1, -1},
};

// Later rechecks: well-separated measured directions (screen, USB and
// left/right axes), captured after the gyro/mag stages for thermal evidence.
static const uint8_t kAccelRecheckPose[kAccelRecheckPoses] = {0, 2, 4};

// Documented body mapping of the device frame (AtomS3R_ImuCal.h).
static inline Eigen::Vector3f accelDocTop()   { return Eigen::Vector3f(1, 0, 0); }
static inline Eigen::Vector3f accelDocRight() { return Eigen::Vector3f(0, 1, 0); }
static inline Eigen::Vector3f accelDocOut()   { return Eigen::Vector3f(0, 0, -1); }

// Body-frame up direction of a pose under the documented mapping.
static inline Eigen::Vector3f accelPoseDocUp(const AccelPoseDef& p) {
  const Eigen::Vector3f v = (float)p.top * accelDocTop() + (float)p.right * accelDocRight() +
                            (float)p.out * accelDocOut();
  return v.normalized();
}

// Text-up direction in the device frame for a display rotation delta
// (0: top end, 1: right edge, 2: USB end, 3: left edge).
static inline void accelRotTextUp(uint8_t rot_delta, int& top, int& right) {
  switch (rot_delta & 3) {
    case 0: top = 1; right = 0; break;
    case 1: top = 0; right = 1; break;
    case 2: top = -1; right = 0; break;
    default: top = 0; right = -1; break;
  }
}

static inline const char* accelAxisName(int axis) {
  switch (axis) {
    case 0: return "USB axis";
    case 1: return "L/R axis";
    default: return "screen axis";
  }
}

// Configuration

struct AccelCaptureCfg {
  float g = 9.80665f;

  // Timing (session clock from sample timestamps)
  uint32_t block_ms = 250;            // one retained observation
  uint32_t min_after_start_ms = 1000; // ignore settling right after the tap
  uint32_t place_ms = 6500;           // placement guidance shown on the bar
  uint32_t hold_useful_ms = 5000;     // retained observation needed per hold
  uint32_t recheck_verify_ms = 2000;  // extra verification-only time per recheck
  uint32_t hold_timeout_ms = 30000;   // from the first sample of the hold
  uint32_t stuck_ms = 12000;          // no fresh samples at all
  uint32_t max_gap_ms = 60;           // a longer gap voids the current block

  // Settling
  int settle_blocks = 3;              // consecutive qualified, consistent blocks
  float settle_step = 0.30f;          // m/s^2 between consecutive block means (~1.8 deg)
  int resettle_after = 4;             // consecutive rejected blocks -> settle again

  // Sample gate (broad: must not reject the bias being calibrated)
  float raw_norm_tol = 2.0f;          // | ||a|| - g | per sample, m/s^2

  // Block gates
  int block_min_samples = 5;
  float block_min_cover = 0.6f;       // sample time coverage of the block
  float block_max_std = 0.25f;        // m/s^2 RMS about the block mean (hand-held OK)
  float block_max_range = 1.6f;       // m/s^2 per-axis peak-to-peak (shock)
  // Rotation gates. A 0.25 s block smeared by a steady 0.15 rad/s turn spans
  // 2.1 deg: its mean norm shrinks by up to ~0.06 mg and the centripetal term at a
  // 5 cm lever arm is ~0.11 mg. Hand wobble does not average inside one block,
  // so these bound physical error rather than demand stillness; placement
  // motion (> 0.5 rad/s) is still rejected.
  float gyro_dev_max = 0.15f;         // rad/s block mean about the stationary gyro level
  float gyro_dev_max_unref = 0.20f;   // rad/s before a stationary level is known
  float gyro_std_max = 0.30f;         // rad/s within a block
  float capture_step = 0.5f;          // m/s^2 block-to-block change while capturing
  float norm_step = 0.12f;            // m/s^2 block-to-block norm change while capturing (knocks)

  // Regions
  float face_cone_deg = 35.0f;
  float tilt_cone_deg = 35.0f;

  // Bounds on user effort
  int max_attempts_per_hold = 3;      // first try + 2 retries
  int max_extra_holds = 2;
};

// Hold engine

enum class AccelHoldPhase : uint8_t { PLACING = 0, CAPTURING = 1, DONE = 2, FAILED = 3 };
enum class AccelHoldFail : uint8_t { NONE = 0, MOTION = 1, WRONG_POSE = 2, NO_SAMPLES = 3, STORE_FULL = 4 };
enum class AccelHoldHint : uint8_t { PLACE = 0, SETTLING = 1, HOLD_STILL = 2, MOVING = 3, CHECK_POSE = 4 };

inline const char* accelHoldFailStr(AccelHoldFail f) {
  switch (f) {
    case AccelHoldFail::NONE: return "OK";
    case AccelHoldFail::MOTION: return "Too much motion";
    case AccelHoldFail::WRONG_POSE: return "Pose not matched";
    case AccelHoldFail::NO_SAMPLES: return "No IMU samples";
    case AccelHoldFail::STORE_FULL: return "Buffer full";
    default: return "Failed";
  }
}

inline const char* accelHoldHintStr(AccelHoldHint h) {
  switch (h) {
    case AccelHoldHint::PLACE: return "Place now";
    case AccelHoldHint::SETTLING: return "Hold still...";
    case AccelHoldHint::HOLD_STILL: return "Hold still";
    case AccelHoldHint::MOVING: return "Moving - hold still";
    case AccelHoldHint::CHECK_POSE: return "Check the pose";
    default: return "";
  }
}

struct AccelRawSample {
  uint32_t t_us = 0;
  Eigen::Vector3f a = Eigen::Vector3f::Zero();   // body frame, m/s^2
  Eigen::Vector3f w = Eigen::Vector3f::Zero();   // body frame, rad/s
  float tempC = NAN;
};

struct AccelHoldStats {
  uint32_t samples = 0, stale = 0, gaps = 0, raw_gated = 0;
  uint32_t blocks = 0, kept = 0;
  uint32_t rej_short = 0, rej_motion = 0, rej_shock = 0, rej_gyro = 0, rej_region = 0,
           rej_norm = 0, rej_step = 0;
  uint32_t resettles = 0;
  uint32_t t_settle_ms = 0;    // from hold start to first settle
  uint32_t t_total_ms = 0;     // from hold start to completion
};

struct AccelHoldView {
  AccelHoldPhase phase = AccelHoldPhase::PLACING;
  AccelHoldHint hint = AccelHoldHint::PLACE;
  const char* hint_text = "";  // accelHoldHintStr(hint)
  float progress01 = 0;        // placement bar while placing, capture bar after
  uint32_t elapsed_ms = 0;
};

struct AccelStepView {
  char title[12];              // size-2 font: <= 10 chars
  char label[22];              // <= 21 chars
  const char* lines[3];        // pose instruction lines
  char note[22];               // reason line (extras) or empty
  uint8_t rot_delta;
  AccelHoldKind kind;
  uint8_t pose;
  uint8_t attempt;
};

// Region test: returns true when `up` (unit) is acceptable for the pose.
struct AccelRegionSpec {
  AccelRegion rule = AccelRegion::FRAME;
  Eigen::Vector3f expected = Eigen::Vector3f::Zero();  // for cone rules
  Eigen::Vector3f screen_out = Eigen::Vector3f::Zero(); // for VERTICAL_AXIS
  float cos_cone = 0.82f;
  float sin_cone = 0.57f;

  bool accepts(const Eigen::Vector3f& up) const {
    if (rule == AccelRegion::VERTICAL_AXIS) {
      if (std::fabs(up.dot(screen_out)) > sin_cone) return false;
      // Near one of the in-plane body axes: the IMU axes are parallel to the edges.
      const Eigen::Vector3f p = up - up.dot(screen_out) * screen_out;
      const float pn = p.norm();
      if (!(pn > 0.1f)) return false;
      const Eigen::Vector3f q = p / pn;
      return std::max(std::fabs(q.x()), std::max(std::fabs(q.y()), std::fabs(q.z()))) >= cos_cone;
    }
    return up.dot(expected) >= cos_cone;
  }
};

template <int MAXO>
class AccelHoldEngine {
public:
  void begin(const AccelCaptureCfg& cfg, const AccelRegionSpec& region, AccelObs* store,
             int* store_n, uint8_t hold_index, uint32_t fit_need_ms, uint32_t verify_need_ms,
             const Eigen::Vector3f& gyro_ref, bool gyro_ref_valid)
  {
    cfg_ = &cfg; region_ = region; store_ = store; store_n_ = store_n;
    hold_ = hold_index;
    first_obs_ = *store_n;
    fit_need_ = (fit_need_ms + cfg.block_ms - 1) / cfg.block_ms;
    ver_need_ = (verify_need_ms + cfg.block_ms - 1) / cfg.block_ms;
    gyro_ref_ = gyro_ref; gyro_ref_valid_ = gyro_ref_valid;
    phase_ = AccelHoldPhase::PLACING;
    fail_ = AccelHoldFail::NONE;
    hint_ = AccelHoldHint::PLACE;
    stats_ = AccelHoldStats{};
    started_ = false;
    acc_us_ = 0;
    t_ms_ = 0;
    have_prev_sample_ = false;
    run_ = 0; rej_run_ = 0;
    n_fit_ = n_ver_ = 0;
    have_last_kept_ = false;
    blk_.reset();
    up_sum_.setZero();
    wsum_.setZero(); wcnt_ = 0;
  }

  // Feeds one fresh sample. Returns the phase after it.
  AccelHoldPhase feed(const AccelRawSample& s) {
    if (phase_ == AccelHoldPhase::DONE || phase_ == AccelHoldPhase::FAILED) return phase_;
    if (!started_) {
      started_ = true;
      t0_us_ = s.t_us;
      t_ms_ = 0;
      last_us_ = s.t_us;
    } else {
      const uint32_t dt = (uint32_t)(s.t_us - last_us_);
      if (dt == 0 || dt > 0x7FFFFFFFu) { ++stats_.stale; return phase_; }  // non-increasing time
      acc_us_ += dt;
      last_us_ = s.t_us;
      t_ms_ = (uint32_t)(acc_us_ / 1000u);
    }
    // Stale repeat: identical accelerometer and gyro words.
    if (have_prev_sample_ && s.a == prev_a_ && s.w == prev_w_) { ++stats_.stale; return timeout_check_(); }
    have_prev_sample_ = true; prev_a_ = s.a; prev_w_ = s.w;

    if (!s.a.allFinite() || !s.w.allFinite() ||
        std::fabs(s.a.norm() - cfg_->g) > cfg_->raw_norm_tol) {
      ++stats_.raw_gated;
      blk_.reset();              // a gated sample voids the block it falls in
      return timeout_check_();
    }
    ++stats_.samples;

    if (blk_.n > 0 && (t_ms_ - blk_.t_last_ms) > cfg_->max_gap_ms) {
      ++stats_.gaps;
      blk_.reset();
    }
    if (blk_.n == 0) blk_.start(t_ms_, s.a, s.w);
    blk_.add(t_ms_, s.a, s.w, s.tempC);
    if (t_ms_ - blk_.t_start_ms >= cfg_->block_ms) {
      closeBlock_();
      blk_.reset();
    }
    return timeout_check_();
  }

  // Timeout without samples is the caller's job (stuck detection); this only
  // advances the hold clock when samples arrive.
  AccelHoldPhase phase() const { return phase_; }
  AccelHoldFail failure() const { return fail_; }
  const AccelHoldStats& stats() const { return stats_; }
  void failNoSamples() { fail_ = AccelHoldFail::NO_SAMPLES; phase_ = AccelHoldPhase::FAILED; rollback_(); }

  AccelHoldView view() const {
    AccelHoldView v;
    v.phase = phase_;
    v.hint = hint_;
    v.hint_text = accelHoldHintStr(hint_);
    v.elapsed_ms = t_ms_;
    const uint32_t need = fit_need_ + ver_need_;
    if (phase_ == AccelHoldPhase::PLACING && n_fit_ + n_ver_ == 0) {
      v.progress01 = std::min(1.0f, (float)t_ms_ / (float)cfg_->place_ms);
    } else {
      v.progress01 = need ? std::min(1.0f, (float)(n_fit_ + n_ver_) / (float)need) : 1.0f;
    }
    return v;
  }

  // Mean measured up direction of the retained blocks (unit).
  Eigen::Vector3f upMean() const {
    const float n = up_sum_.norm();
    return (n > 0) ? Eigen::Vector3f(up_sum_ / n) : Eigen::Vector3f::Zero();
  }
  // Mean gyro of the retained blocks (stationary level estimate).
  bool gyroMean(Eigen::Vector3f& w) const {
    if (wcnt_ == 0) return false;
    w = wsum_ / (float)wcnt_;
    return true;
  }
  int firstObs() const { return first_obs_; }
  int keptBlocks() const { return (int)(n_fit_ + n_ver_); }

private:
  struct Block {
    uint32_t n = 0;
    uint32_t t_start_ms = 0, t_first_ms = 0, t_last_ms = 0, max_gap_ms = 0;
    Eigen::Vector3f a_ref, sa, saa, amin, amax;
    Eigen::Vector3f w_ref, sw, sww;
    float st = 0; int nt = 0;
    void reset() { n = 0; }
    void start(uint32_t t, const Eigen::Vector3f& a, const Eigen::Vector3f& w) {
      n = 0; t_start_ms = t; t_first_ms = t; t_last_ms = t; max_gap_ms = 0;
      a_ref = a; w_ref = w;
      sa.setZero(); saa.setZero(); sw.setZero(); sww.setZero();
      amin = a; amax = a;
      st = 0; nt = 0;
    }
    void add(uint32_t t, const Eigen::Vector3f& a, const Eigen::Vector3f& w, float tC) {
      if (n > 0) max_gap_ms = std::max(max_gap_ms, t - t_last_ms);
      t_last_ms = t;
      const Eigen::Vector3f da = a - a_ref, dw = w - w_ref;
      sa += da; saa += da.cwiseProduct(da);
      sw += dw; sww += dw.cwiseProduct(dw);
      amin = amin.cwiseMin(a); amax = amax.cwiseMax(a);
      if (std::isfinite(tC)) { st += tC; ++nt; }
      ++n;
    }
  };

  AccelHoldPhase timeout_check_() {
    if (phase_ != AccelHoldPhase::DONE && phase_ != AccelHoldPhase::FAILED &&
        t_ms_ >= cfg_->hold_timeout_ms) {
      fail_ = (stats_.rej_region > stats_.rej_motion + stats_.rej_gyro + stats_.rej_shock)
              ? AccelHoldFail::WRONG_POSE : AccelHoldFail::MOTION;
      phase_ = AccelHoldPhase::FAILED;
      rollback_();
    }
    return phase_;
  }

  void rollback_() { *store_n_ = first_obs_; n_fit_ = n_ver_ = 0; }

  void closeBlock_() {
    ++stats_.blocks;
    const Block& b = blk_;
    const float nf = (float)b.n;
    const uint32_t cover = b.t_last_ms - b.t_first_ms;
    if ((int)b.n < cfg_->block_min_samples || (float)cover < cfg_->block_min_cover * (float)cfg_->block_ms ||
        b.max_gap_ms > cfg_->max_gap_ms) {
      ++stats_.rej_short; reject_(false); return;
    }
    const Eigen::Vector3f ma = b.sa / nf;
    const Eigen::Vector3f mean_a = b.a_ref + ma;
    const Eigen::Vector3f var_a = (b.saa / nf - ma.cwiseProduct(ma)).cwiseMax(0.0f);
    const float a_std = std::sqrt(var_a.sum());
    const Eigen::Vector3f mw = b.sw / nf;
    const Eigen::Vector3f mean_w = b.w_ref + mw;
    const Eigen::Vector3f var_w = (b.sww / nf - mw.cwiseProduct(mw)).cwiseMax(0.0f);
    const float w_std = std::sqrt(var_w.sum());
    const float range = (b.amax - b.amin).maxCoeff();

    if (a_std > cfg_->block_max_std) { ++stats_.rej_motion; reject_(false); return; }
    if (range > cfg_->block_max_range) { ++stats_.rej_shock; reject_(false); return; }
    const float wdev = gyro_ref_valid_ ? (mean_w - gyro_ref_).norm() : mean_w.norm();
    const float wmax = gyro_ref_valid_ ? cfg_->gyro_dev_max : cfg_->gyro_dev_max_unref;
    if (wdev > wmax || w_std > cfg_->gyro_std_max) { ++stats_.rej_gyro; reject_(false); return; }
    const float an = mean_a.norm();
    if (!(an > 1e-3f)) { ++stats_.rej_short; reject_(false); return; }
    const Eigen::Vector3f up = mean_a / an;
    if (!region_.accepts(up)) {
      ++stats_.rej_region; reject_(true); return;
    }

    Cand c;
    c.a = mean_a; c.w = mean_w; c.a_std = a_std; c.n = (uint16_t)std::min<uint32_t>(b.n, 65535u);
    c.t_ms = (b.t_first_ms + b.t_last_ms) / 2;
    c.tempC = (b.nt > 0) ? b.st / (float)b.nt : NAN;

    if (phase_ == AccelHoldPhase::PLACING) {
      if (t_ms_ < cfg_->min_after_start_ms) { run_ = 0; hint_ = AccelHoldHint::PLACE; return; }
      if (run_ > 0 && (c.a - run_buf_[run_ - 1].a).norm() > cfg_->settle_step) run_ = 0;
      if (run_ < kMaxSettle) run_buf_[run_++] = c;
      hint_ = AccelHoldHint::SETTLING;
      if (run_ >= cfg_->settle_blocks) {
        if (stats_.t_settle_ms == 0) stats_.t_settle_ms = t_ms_;
        phase_ = AccelHoldPhase::CAPTURING;
        hint_ = AccelHoldHint::HOLD_STILL;
        for (int j = 0; j < run_; ++j) keep_(run_buf_[j]);
        run_ = 0; rej_run_ = 0;
      }
      return;
    }

    // CAPTURING: temporal consistency with the last retained block.
    if (have_last_kept_) {
      if ((c.a - last_kept_a_).norm() > cfg_->capture_step) { ++stats_.rej_step; reject_(false); return; }
      if (std::fabs(c.a.norm() - last_kept_a_.norm()) > cfg_->norm_step) { ++stats_.rej_norm; reject_(false); return; }
    }
    rej_run_ = 0;
    hint_ = AccelHoldHint::HOLD_STILL;
    keep_(c);
  }

  struct Cand {
    Eigen::Vector3f a, w;
    float a_std, tempC;
    uint32_t t_ms;
    uint16_t n;
  };

  void reject_(bool region) {
    if (phase_ == AccelHoldPhase::PLACING) {
      run_ = 0;
      hint_ = region ? AccelHoldHint::CHECK_POSE
                     : (t_ms_ < cfg_->min_after_start_ms ? AccelHoldHint::PLACE : AccelHoldHint::SETTLING);
      return;
    }
    ++rej_run_;
    hint_ = region ? AccelHoldHint::CHECK_POSE : AccelHoldHint::MOVING;
    if (rej_run_ >= cfg_->resettle_after) {
      // Sustained disturbance: keep what was retained, settle again before
      // retaining more (the device may have been re-placed).
      phase_ = AccelHoldPhase::PLACING;
      ++stats_.resettles;
      run_ = 0; rej_run_ = 0;
      have_last_kept_ = false;
    }
  }

  void keep_(const Cand& c) {
    if (phase_ == AccelHoldPhase::DONE) return;
    if (*store_n_ >= MAXO) { fail_ = AccelHoldFail::STORE_FULL; phase_ = AccelHoldPhase::FAILED; rollback_(); return; }
    AccelObs& o = store_[(*store_n_)++];
    for (int d = 0; d < 3; ++d) { o.a[d] = c.a(d); o.w[d] = c.w(d); }
    o.tempC = c.tempC;
    o.a_std = c.a_std;
    o.t_ms = c.t_ms;   // hold-relative; the procedure rebases to the session clock
    o.n = c.n;
    o.hold = hold_;
    const bool fit = (n_fit_ < fit_need_);
    o.role = (uint8_t)(fit ? AccelObsRole::FIT : AccelObsRole::VERIFY);
    if (fit) ++n_fit_; else ++n_ver_;
    ++stats_.kept;
    up_sum_ += c.a.normalized();
    wsum_ += c.w; ++wcnt_;
    last_kept_a_ = c.a; have_last_kept_ = true;
    if (n_fit_ >= fit_need_ && n_ver_ >= ver_need_) {
      phase_ = AccelHoldPhase::DONE;
      stats_.t_total_ms = t_ms_;
    }
  }

  static constexpr int kMaxSettle = 8;

  const AccelCaptureCfg* cfg_ = nullptr;
  AccelRegionSpec region_;
  AccelObs* store_ = nullptr;
  int* store_n_ = nullptr;
  uint8_t hold_ = 0;
  int first_obs_ = 0;
  uint32_t fit_need_ = 0, ver_need_ = 0, n_fit_ = 0, n_ver_ = 0;
  Eigen::Vector3f gyro_ref_ = Eigen::Vector3f::Zero();
  bool gyro_ref_valid_ = false;

  AccelHoldPhase phase_ = AccelHoldPhase::PLACING;
  AccelHoldFail fail_ = AccelHoldFail::NONE;
  AccelHoldHint hint_ = AccelHoldHint::PLACE;
  AccelHoldStats stats_;

  bool started_ = false;
  uint32_t t0_us_ = 0, last_us_ = 0;
  uint64_t acc_us_ = 0;
  uint32_t t_ms_ = 0;
  bool have_prev_sample_ = false;
  Eigen::Vector3f prev_a_, prev_w_;

  Block blk_;
  Cand run_buf_[kMaxSettle];
  int run_ = 0, rej_run_ = 0;
  bool have_last_kept_ = false;
  Eigen::Vector3f last_kept_a_ = Eigen::Vector3f::Zero();
  Eigen::Vector3f up_sum_ = Eigen::Vector3f::Zero();
  Eigen::Vector3f wsum_ = Eigen::Vector3f::Zero();
  int wcnt_ = 0;
};

// Procedure

// Fit work handed to the platform so it can run on a large-stack task.
class AccelFitJob {
public:
  virtual ~AccelFitJob() {}
  virtual void run() = 0;
};

class AccelCalIo {
public:
  virtual ~AccelCalIo() {}
  // Preparation screen with the full instructions; returns after the tap.
  // false = user abort.
  virtual bool prep(const AccelStepView& v) = 0;
  // Polls the IMU; true with a fresh sample.
  virtual bool sample(AccelRawSample& s) = 0;
  virtual uint32_t nowMs() = 0;
  // Called every loop during a hold; the implementation throttles drawing.
  virtual void capture(const AccelStepView& v, const AccelHoldView& h) = 0;
  virtual void holdOk(const AccelStepView& v) = 0;
  // A hold failed and may be retried; false = user abort.
  virtual bool holdRetry(const AccelStepView& v, const char* why) = 0;
  // Runs job.run() (on a task with enough stack); false on timeout/failure.
  virtual bool runFit(AccelFitJob& job, const char* what) = 0;
  virtual void log(const char* line) = 0;
  virtual void idle() = 0;
};

enum class AccelProcFail : uint8_t {
  NONE = 0, ABORTED, HOLD_FAILED, FIT_FAILED, WEAK_GEOMETRY, FIT_TASK,
};

struct AccelHoldRecord {
  uint8_t pose = 0;
  AccelHoldKind kind = AccelHoldKind::MAIN;
  uint8_t attempts = 0;
  uint16_t first_obs = 0, n_obs = 0;
  uint32_t t_start_ms = 0;          // session clock at the hold's first sample
  Eigen::Vector3f up = Eigen::Vector3f::Zero();
  AccelHoldStats stats;
};

template <int MAXO = 340, int MAXH = 24>
class AccelCalProcedure {
public:
  using Fitter = AccelFullFitter<MAXO, MAXH>;

  void begin(const AccelCaptureCfg& ccfg, const AccelFitCfg& fcfg, const AccelThermalPrior& prior,
             const Eigen::Vector3f& gyro_bias, bool gyro_bias_valid)
  {
    ccfg_ = ccfg; fcfg_ = fcfg; prior_ = prior;
    n_obs_ = 0; n_holds_ = 0; n_extra_ = 0;
    gyro_ref_ = gyro_bias; gyro_ref_valid_ = gyro_bias_valid; gyro_prior_ = gyro_bias_valid;
    for (int i = 0; i < kAccelMainPoses; ++i) { have_up_[i] = false; }
    fail_ = AccelProcFail::NONE;
    fail_detail_[0] = 0;
    fail_pose_ = -1;
    fail_hold_ = AccelHoldFail::NONE;
    session_t0_set_ = false;
    n_gyro_means_ = 0;
    total_hold_ms_ = 0; total_attempts_ = 0; total_retries_ = 0;
    prelim_ = AccelFullFitResult{};
    result_ = AccelFullFitResult{};
  }

  // Ten main holds, then a static fit; bounded targeted extra holds while the
  // measured geometry leaves a bias or cross term undetermined.
  bool runMainStage(AccelCalIo& io) {
    for (int p = 0; p < kAccelMainPoses; ++p) {
      if (!captureHold_(io, AccelHoldKind::MAIN, (uint8_t)p, p + 1, kAccelMainPoses, nullptr)) return false;
    }
    logFrame_(io);
    for (;;) {
      FitRunner job(*this, true);
      if (!io.runFit(job, "ACCEL")) return setFail_(AccelProcFail::FIT_TASK, "Fit task failed");
      logFit_(io, "prelim", prelim_);
      if (prelim_.ok) return true;
      bool geometric = (prelim_.gate == AccelGate::INFO_BIAS || prelim_.gate == AccelGate::INFO_CROSS);
      if (prelim_.gate == AccelGate::CROSSVAL) {
        // A studentized error over its limit means a hold disagrees with the
        // rest: redo it. A raw-only exceedance means the other holds predict
        // that hold poorly (high leverage): add the pose that best constrains
        // the weakest term instead of repeating one.
        if (prelim_.cv_tmax > fcfg_.max_cv_t || prelim_.weak_param < 0) {
          if (!recapture_(io, prelim_.cv_worst_hold)) return false;
          continue;
        }
        geometric = true;
      }
      if (!geometric || prelim_.weak_param < 0) return setFail_(AccelProcFail::FIT_FAILED, gateText_(prelim_.gate));
      if (n_extra_ >= ccfg_.max_extra_holds) return setFail_(AccelProcFail::WEAK_GEOMETRY, weakText_(prelim_.weak_param));
      const int pose = chooseExtra_(prelim_);
      char note[22];
      snprintf(note, sizeof(note), "%s", weakText_(prelim_.weak_param));
      ++n_extra_;
      if (!captureHold_(io, AccelHoldKind::EXTRA, (uint8_t)pose, n_extra_, ccfg_.max_extra_holds, note)) return false;
    }
  }

  // Three later rechecks at well-separated measured directions. The first
  // hold_useful_ms of each is fit data (thermal evidence); the next
  // recheck_verify_ms is verification-only and never enters the fit.
  bool runRecheckStage(AccelCalIo& io) {
    for (int r = 0; r < kAccelRecheckPoses; ++r) {
      if (!captureHold_(io, AccelHoldKind::RECHECK, kAccelRecheckPose[r], r + 1, kAccelRecheckPoses, nullptr)) return false;
    }
    return true;
  }

  // Final fit. A hold that fails the held-out or verification check is
  // recaptured (bounded by max_extra_holds, shared with geometry extras).
  bool runFinalFit(AccelCalIo& io) {
    for (;;) {
      FitRunner job(*this, false);
      if (!io.runFit(job, "ACCEL")) return setFail_(AccelProcFail::FIT_TASK, "Fit task failed");
      logFit_(io, "final", result_);
      if (result_.ok) return true;
      int bad = -1;
      if (result_.gate == AccelGate::CROSSVAL) bad = result_.cv_worst_hold;
      else if (result_.gate == AccelGate::VERIFY) bad = result_.ver_worst_hold;
      if (bad < 0) return setFail_(AccelProcFail::FIT_FAILED, gateText_(result_.gate));
      if (!recapture_(io, bad)) return false;
    }
  }

  // Replaces the stationary gyro level used by the rotation gate (e.g. after
  // the gyro stage of the full wizard).
  void setGyroReference(const Eigen::Vector3f& w) { gyro_ref_ = w; gyro_ref_valid_ = true; gyro_prior_ = true; }

  const AccelFullFitResult& result() const { return result_; }
  const AccelFullFitResult& preliminary() const { return prelim_; }
  AccelProcFail failure() const { return fail_; }
  const char* failureDetail() const { return fail_detail_; }
  const AccelObs* obs() const { return obs_; }
  int nObs() const { return n_obs_; }
  int nHolds() const { return n_holds_; }
  const AccelHoldRecord& hold(int i) const { return holds_[i]; }
  int nExtra() const { return n_extra_; }
  int totalAttempts() const { return total_attempts_; }
  int totalRetries() const { return total_retries_; }
  uint32_t totalHoldMs() const { return total_hold_ms_; }
  const AccelCaptureCfg& captureCfg() const { return ccfg_; }
  const AccelFitCfg& fitCfg() const { return fcfg_; }

  // Two short screen lines (<= 21 chars) describing a procedure failure.
  void failLines(char l1[22], char l2[22]) const {
    switch (fail_) {
      case AccelProcFail::ABORTED:
        snprintf(l1, 22, "Aborted"); snprintf(l2, 22, "Previous cal kept"); break;
      case AccelProcFail::HOLD_FAILED:
        snprintf(l1, 22, "%s", fail_pose_ >= 0 ? kAccelPoses[fail_pose_].name : "Hold failed");
        snprintf(l2, 22, "%s", accelHoldFailStr(fail_hold_)); break;
      case AccelProcFail::WEAK_GEOMETRY:
        snprintf(l1, 22, "%s", fail_detail_); snprintf(l2, 22, "Tilt poses further"); break;
      case AccelProcFail::FIT_TASK:
        snprintf(l1, 22, "Fit task failed"); snprintf(l2, 22, "See Serial log"); break;
      case AccelProcFail::FIT_FAILED:
      default:
        snprintf(l1, 22, "Fit rejected"); snprintf(l2, 22, "%.21s", fail_detail_); break;
    }
  }

  // Measured device frame (valid after the six faces).
  bool frame(Eigen::Vector3f& top, Eigen::Vector3f& right, Eigen::Vector3f& out) const {
    if (!(have_up_[0] && have_up_[2])) return false;
    out = up_[0];
    top = -(up_[2] - up_[2].dot(out) * out);
    if (!(top.norm() > 0.1f)) return false;
    top.normalize();
    right = top.cross(out);   // right-handed: top x out = right
    if (have_up_[4]) {
      Eigen::Vector3f r = up_[4] - up_[4].dot(out) * out - up_[4].dot(top) * top;
      if (r.norm() > 0.1f && r.dot(right) > 0) right = (0.5f * (r.normalized() + right)).normalized();
    }
    return true;
  }

  // Expected up direction of a pose from the measured frame (unit).
  bool expectedUp(uint8_t pose, Eigen::Vector3f& up) const {
    Eigen::Vector3f t, r, o;
    if (!frame(t, r, o)) return false;
    const AccelPoseDef& p = kAccelPoses[pose];
    up = ((float)p.top * t + (float)p.right * r + (float)p.out * o).normalized();
    return true;
  }

  AccelStepView stepView(AccelHoldKind kind, uint8_t pose, int index, int count, const char* note, int attempt) const {
    AccelStepView v;
    memset(&v, 0, sizeof(v));
    const AccelPoseDef& p = kAccelPoses[pose];
    switch (kind) {
      case AccelHoldKind::MAIN:
        snprintf(v.title, sizeof(v.title), "ACCEL");
        snprintf(v.label, sizeof(v.label), "%d/%d %s", index, count, p.name);
        break;
      case AccelHoldKind::EXTRA:
        snprintf(v.title, sizeof(v.title), "EXTRA %d/%d", index, count);
        snprintf(v.label, sizeof(v.label), "%s", p.name);
        break;
      default:
        snprintf(v.title, sizeof(v.title), "CHECK %d/%d", index, count);
        snprintf(v.label, sizeof(v.label), "%s", p.name);
        break;
    }
    for (int j = 0; j < 3; ++j) v.lines[j] = p.lines[j];
    if (note) snprintf(v.note, sizeof(v.note), "%s", note);
    v.rot_delta = p.rot_delta;
    v.kind = kind;
    v.pose = pose;
    v.attempt = (uint8_t)attempt;
    return v;
  }

  static const char* weakText_(int param) {
    static const char* kText[6] = {
      "Weak bias: USB axis", "Weak bias: L/R axis", "Weak bias: screen",
      "Weak cross: USB-L/R", "Weak cross: USB-scrn", "Weak cross: L/R-scrn",
    };
    return (param >= 0 && param < 6) ? kText[param] : "Weak geometry";
  }

private:
  class FitRunner : public AccelFitJob {
  public:
    FitRunner(AccelCalProcedure& p, bool prelim) : p_(p), prelim_(prelim) {}
    void run() override {
      if (prelim_) p_.fitter_.fit(p_.obs_, p_.n_obs_, p_.fcfg_, p_.prior_, p_.prelim_, true);
      else p_.fitter_.fit(p_.obs_, p_.n_obs_, p_.fcfg_, p_.prior_, p_.result_, false);
    }
  private:
    AccelCalProcedure& p_;
    bool prelim_;
  };

  // Removes hold h (its blocks and record) and captures the same pose again.
  bool recapture_(AccelCalIo& io, int h) {
    if (h < 0 || h >= n_holds_) return setFail_(AccelProcFail::FIT_FAILED, "Holds disagree");
    if (n_extra_ >= ccfg_.max_extra_holds) {
      char d[40];
      snprintf(d, sizeof(d), "Bad: %s", kAccelPoses[holds_[h].pose].name);
      return setFail_(AccelProcFail::FIT_FAILED, d);
    }
    const AccelHoldRecord rec = holds_[h];
    char line[96];
    snprintf(line, sizeof(line), "[ACC] recapture hold=%d pose=%s (inconsistent with the other holds)",
             h, kAccelPoses[rec.pose].name);
    io.log(line);
    dropHold_(h);
    ++n_extra_;
    const int idx = (rec.kind == AccelHoldKind::RECHECK) ? recheckIndex_(rec.pose) : n_extra_;
    const int cnt = (rec.kind == AccelHoldKind::RECHECK) ? kAccelRecheckPoses : ccfg_.max_extra_holds;
    const AccelHoldKind kind = (rec.kind == AccelHoldKind::RECHECK) ? AccelHoldKind::RECHECK : AccelHoldKind::EXTRA;
    return captureHold_(io, kind, rec.pose, idx, cnt, "Redo: inconsistent");
  }

  static const char* gateText_(AccelGate g) {
    switch (g) {
      case AccelGate::COVERAGE: return "Poor coverage";
      case AccelGate::DIAG: case AccelGate::COND: case AccelGate::OFFDIAG: return "Unphysical scale";
      case AccelGate::BIAS_MAG: return "Bias out of range";
      case AccelGate::RESID:
      case AccelGate::CROSSVAL: return "Holds disagree";
      case AccelGate::VERIFY: return "Recheck mismatch";
      case AccelGate::FLOAT: return "Float check failed";
      case AccelGate::SOLVE: return "Solver failed";
      case AccelGate::INFO_BIAS: case AccelGate::INFO_CROSS: return "Weak geometry";
      default: return "Rejected";
    }
  }

  static int recheckIndex_(uint8_t pose) {
    for (int r = 0; r < kAccelRecheckPoses; ++r) if (kAccelRecheckPose[r] == pose) return r + 1;
    return 1;
  }

  void dropHold_(int h) {
    const int first = holds_[h].first_obs, cnt = holds_[h].n_obs;
    for (int i = first; i + cnt < n_obs_; ++i) obs_[i] = obs_[i + cnt];
    n_obs_ -= cnt;
    for (int i = 0; i < n_obs_; ++i) if (obs_[i].hold > h) obs_[i].hold--;
    for (int j = h; j + 1 < n_holds_; ++j) {
      holds_[j] = holds_[j + 1];
      holds_[j].first_obs -= (uint16_t)cnt;
    }
    --n_holds_;
  }

  bool setFail_(AccelProcFail f, const char* detail) {
    fail_ = f;
    snprintf(fail_detail_, sizeof(fail_detail_), "%s", detail ? detail : "");
    return false;
  }

  AccelRegionSpec regionFor_(AccelHoldKind kind, uint8_t pose) const {
    AccelRegionSpec rs;
    const AccelPoseDef& p = kAccelPoses[pose];
    const float cone = (p.region == AccelRegion::FRAME) ? ccfg_.tilt_cone_deg : ccfg_.face_cone_deg;
    rs.cos_cone = std::cos(cone * 3.14159265f / 180.0f);
    rs.sin_cone = std::sin(cone * 3.14159265f / 180.0f);
    if (kind != AccelHoldKind::MAIN || p.region == AccelRegion::FRAME) {
      // Extras, rechecks and tilts: the measured frame.
      rs.rule = AccelRegion::FRAME;
      if (have_up_[pose] && kind != AccelHoldKind::MAIN) rs.expected = up_[pose];
      else if (!expectedUp(pose, rs.expected)) rs.expected = accelPoseDocUp(p);
      return rs;
    }
    rs.rule = p.region;
    switch (p.region) {
      case AccelRegion::SCREEN_UP_DOC: rs.expected = accelPoseDocUp(p); rs.rule = AccelRegion::FRAME; break;
      case AccelRegion::OPPOSITE: rs.expected = -up_[p.ref0]; rs.rule = AccelRegion::FRAME; break;
      case AccelRegion::VERTICAL_AXIS: rs.screen_out = up_[p.ref0]; break;
      case AccelRegion::RIGHT_HANDED: rs.expected = up_[p.ref0].cross(up_[p.ref1]).normalized(); rs.rule = AccelRegion::FRAME; break;
      default: break;
    }
    return rs;
  }

  bool captureHold_(AccelCalIo& io, AccelHoldKind kind, uint8_t pose, int index, int count, const char* note) {
    if (n_holds_ >= MAXH) return setFail_(AccelProcFail::HOLD_FAILED, "Too many holds");
    const bool recheck = (kind == AccelHoldKind::RECHECK);
    for (int attempt = 1; attempt <= ccfg_.max_attempts_per_hold; ++attempt) {
      AccelStepView v = stepView(kind, pose, index, count, note, attempt);
      if (!io.prep(v)) return setFail_(AccelProcFail::ABORTED, "Aborted");
      ++total_attempts_;
      if (attempt > 1) ++total_retries_;

      const AccelRegionSpec rs = regionFor_(kind, pose);
      const int hold_index = n_holds_;
      engine_.begin(ccfg_, rs, obs_, &n_obs_, (uint8_t)hold_index, ccfg_.hold_useful_ms,
                    recheck ? ccfg_.recheck_verify_ms : 0, gyro_ref_, gyro_ref_valid_);
      uint32_t last_sample_ms = io.nowMs();
      bool first = true;
      uint32_t hold_t0_session = 0;
      while (engine_.phase() == AccelHoldPhase::PLACING || engine_.phase() == AccelHoldPhase::CAPTURING) {
        AccelRawSample s;
        if (io.sample(s)) {
          if (first) {
            hold_t0_session = sessionMs_(s.t_us);
            first = false;
          } else {
            (void)sessionMs_(s.t_us);
          }
          last_sample_ms = io.nowMs();
          engine_.feed(s);
        } else if ((uint32_t)(io.nowMs() - last_sample_ms) > ccfg_.stuck_ms) {
          engine_.failNoSamples();
          break;
        }
        io.capture(v, engine_.view());
        io.idle();
      }

      const AccelHoldStats& st = engine_.stats();
      char line[320];
      const Eigen::Vector3f up = engine_.upMean();
      snprintf(line, sizeof(line),
               "[ACC] hold=%d pose=%s kind=%d att=%d %s t=%.1fs settle=%.1fs kept=%lu blocks=%lu "
               "rej(short/mot/shock/gyro/region/norm/step)=%lu/%lu/%lu/%lu/%lu/%lu/%lu stale=%lu gaps=%lu "
               "resettle=%lu up=(%.3f,%.3f,%.3f)",
               hold_index, kAccelPoses[pose].name, (int)kind, attempt,
               engine_.phase() == AccelHoldPhase::DONE ? "OK" : accelHoldFailStr(engine_.failure()),
               st.t_total_ms / 1000.0, st.t_settle_ms / 1000.0,
               (unsigned long)st.kept, (unsigned long)st.blocks,
               (unsigned long)st.rej_short, (unsigned long)st.rej_motion, (unsigned long)st.rej_shock,
               (unsigned long)st.rej_gyro, (unsigned long)st.rej_region, (unsigned long)st.rej_norm,
               (unsigned long)st.rej_step, (unsigned long)st.stale, (unsigned long)st.gaps,
               (unsigned long)st.resettles, (double)up.x(), (double)up.y(), (double)up.z());
      io.log(line);

      if (engine_.phase() == AccelHoldPhase::DONE) {
        AccelHoldRecord& h = holds_[n_holds_++];
        h.pose = pose; h.kind = kind; h.attempts = (uint8_t)attempt;
        h.first_obs = (uint16_t)engine_.firstObs();
        h.n_obs = (uint16_t)(n_obs_ - engine_.firstObs());
        h.t_start_ms = hold_t0_session;
        h.up = up;
        h.stats = st;
        total_hold_ms_ += st.t_total_ms;
        // Blocks carry hold-relative times; rebase to the session clock.
        for (int i = h.first_obs; i < n_obs_; ++i) {
          obs_[i].t_ms += hold_t0_session;
          char bl[160];
          snprintf(bl, sizeof(bl), "[ACCBLK] %d,%d,%lu,%u,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.3f,%.4f",
                   (int)obs_[i].hold, (int)obs_[i].role, (unsigned long)obs_[i].t_ms, (unsigned)obs_[i].n,
                   (double)obs_[i].a[0], (double)obs_[i].a[1], (double)obs_[i].a[2],
                   (double)obs_[i].w[0], (double)obs_[i].w[1], (double)obs_[i].w[2],
                   (double)obs_[i].tempC, (double)obs_[i].a_std);
          io.log(bl);
        }
        if (kind == AccelHoldKind::MAIN) { up_[pose] = up; have_up_[pose] = true; }
        updateGyroRef_();
        io.holdOk(v);
        return true;
      }
      fail_pose_ = pose;
      fail_hold_ = engine_.failure();
      if (engine_.failure() == AccelHoldFail::STORE_FULL) return setFail_(AccelProcFail::HOLD_FAILED, "Buffer full");
      if (attempt < ccfg_.max_attempts_per_hold) {
        if (!io.holdRetry(v, accelHoldFailStr(engine_.failure()))) return setFail_(AccelProcFail::ABORTED, "Aborted");
      } else {
        char d[40];
        snprintf(d, sizeof(d), "%s: %s", kAccelPoses[pose].name, accelHoldFailStr(engine_.failure()));
        return setFail_(AccelProcFail::HOLD_FAILED, d);
      }
    }
    return setFail_(AccelProcFail::HOLD_FAILED, "Hold failed");
  }

  uint32_t sessionMs_(uint32_t t_us) {
    if (!session_t0_set_) { session_t0_set_ = true; session_last_us_ = t_us; session_us_ = 0; }
    session_us_ += (uint32_t)(t_us - session_last_us_);
    session_last_us_ = t_us;
    return (uint32_t)(session_us_ / 1000u);
  }

  // Stationary gyro level: the prior calibration when there is one, otherwise
  // the per-axis median of completed hold means once three holds exist.
  void updateGyroRef_() {
    if (gyro_prior_) return;
    Eigen::Vector3f w;
    if (!engine_.gyroMean(w)) return;
    if (n_gyro_means_ < MAXH) gyro_means_[n_gyro_means_++] = w;
    if (n_gyro_means_ < 3) return;
    for (int d = 0; d < 3; ++d) {
      float tmp[MAXH];
      for (int i = 0; i < n_gyro_means_; ++i) tmp[i] = gyro_means_[i](d);
      std::nth_element(tmp, tmp + n_gyro_means_ / 2, tmp + n_gyro_means_);
      gyro_ref_(d) = tmp[n_gyro_means_ / 2];
    }
    gyro_ref_valid_ = true;
  }

  int chooseExtra_(const AccelFullFitResult& r) const {
    int best = 0; double bg = -1;
    for (int p = 0; p < kAccelMainPoses; ++p) {
      Eigen::Vector3f up;
      if (have_up_[p]) up = up_[p];
      else if (!expectedUp((uint8_t)p, up)) up = accelPoseDocUp(kAccelPoses[p]);
      const double gain = Fitter::infoGain(r, up.cast<double>(), r.weak_param);
      if (gain > bg) { bg = gain; best = p; }
    }
    return best;
  }

  void logFrame_(AccelCalIo& io) {
    Eigen::Vector3f t, r, o;
    if (!frame(t, r, o)) return;
    const float at = std::acos(std::max(-1.0f, std::min(1.0f, t.dot(accelDocTop())))) * 57.2958f;
    const float ar = std::acos(std::max(-1.0f, std::min(1.0f, r.dot(accelDocRight())))) * 57.2958f;
    const float ao = std::acos(std::max(-1.0f, std::min(1.0f, o.dot(accelDocOut())))) * 57.2958f;
    char line[200];
    snprintf(line, sizeof(line),
             "[ACC] measured frame top=(%.3f,%.3f,%.3f) right=(%.3f,%.3f,%.3f) out=(%.3f,%.3f,%.3f) "
             "vs documented: %.1f/%.1f/%.1f deg",
             (double)t.x(), (double)t.y(), (double)t.z(), (double)r.x(), (double)r.y(), (double)r.z(),
             (double)o.x(), (double)o.y(), (double)o.z(), (double)at, (double)ar, (double)ao);
    io.log(line);
  }

  void logFit_(AccelCalIo& io, const char* what, const AccelFullFitResult& r) {
    char line[320];
    snprintf(line, sizeof(line),
             "[ACCFIT] %s ok=%d reason=%s gate=%s holds=%d blocks=%d rej=%d hold_rms=%.4f "
             "sig_obs=%.4f b_sd=(%.4f,%.4f,%.4f) cross_sd=(%.5f,%.5f,%.5f) thermal=%s/%s "
             "k_sd=(%.5f,%.5f,%.5f) cv=%d rms=%.4f max=%.4f t=%.2f ver=%d rms=%.4f max=%.4f",
             what, (int)r.ok, fitFailStr(r.reason), accelGateStr(r.gate), r.n_fit_holds, r.n_fit_blocks,
             r.n_rejected, r.hold_rms, r.sigma_obs,
             r.bias_sigma(0), r.bias_sigma(1), r.bias_sigma(2),
             r.cross_sigma(0), r.cross_sigma(1), r.cross_sigma(2),
             accelThermalStr(r.thermal), accelThermalReasonStr(r.thermal_reason),
             r.k_sigma(0), r.k_sigma(1), r.k_sigma(2),
             r.n_cv, r.cv_rms, r.cv_max, r.cv_tmax, r.n_ver, r.ver_rms, r.ver_max);
    io.log(line);
    snprintf(line, sizeof(line),
             "[ACCFIT] %s b=(%.5f,%.5f,%.5f) k=(%.6f,%.6f,%.6f) T_ref=%.2f T=[%.2f,%.2f] "
             "S=[%.6f %.6f %.6f; %.6f %.6f %.6f; %.6f %.6f %.6f]",
             what, r.b(0), r.b(1), r.b(2), r.k(0), r.k(1), r.k(2), r.T_ref, r.cal_temp_lo, r.cal_temp_hi,
             r.S(0,0), r.S(0,1), r.S(0,2), r.S(1,0), r.S(1,1), r.S(1,2), r.S(2,0), r.S(2,1), r.S(2,2));
    io.log(line);
  }

  AccelCaptureCfg ccfg_;
  AccelFitCfg fcfg_;
  AccelThermalPrior prior_;

  AccelObs obs_[MAXO];
  int n_obs_ = 0;
  AccelHoldRecord holds_[MAXH];
  int n_holds_ = 0;
  int n_extra_ = 0;

  Eigen::Vector3f up_[kAccelMainPoses];
  bool have_up_[kAccelMainPoses] = {};

  Eigen::Vector3f gyro_ref_ = Eigen::Vector3f::Zero();
  bool gyro_ref_valid_ = false;
  bool gyro_prior_ = false;
  Eigen::Vector3f gyro_means_[MAXH];
  int n_gyro_means_ = 0;

  bool session_t0_set_ = false;
  uint32_t session_last_us_ = 0;
  uint64_t session_us_ = 0;

  uint32_t total_hold_ms_ = 0;
  int total_attempts_ = 0, total_retries_ = 0;

  AccelProcFail fail_ = AccelProcFail::NONE;
  char fail_detail_[40] = {0};
  int fail_pose_ = -1;
  AccelHoldFail fail_hold_ = AccelHoldFail::NONE;

  AccelHoldEngine<MAXO> engine_;
  Fitter fitter_;
  AccelFullFitResult prelim_;
  AccelFullFitResult result_;
};

}  // namespace imu_cal
