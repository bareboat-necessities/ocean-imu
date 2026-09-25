#pragma once

// Replay of logged accelerometer calibration sessions (AtomS3R serial logs).
//
// The wizard prints:
//   [ACCMODE] full|accel_only
//   [ACCPRIOR] valid,k0,k1,k2,k_lo,k_hi,clamp_lo,clamp_hi
//   [ACCGYRO] wx,wy,wz,valid         stationary gyro level and its validity (start / after the gyro stage)
//   [ACCPREP] kind,pose,attempt      a hold's preparation screen
//   [ACCRAW] t_us,ax,ay,az,wx,wy,wz,T   every raw sample (ATOMS3R_ICAL_RAW_LOG=1)
//   [ACCBLK] hold,role,t_ms,n,ax,ay,az,wx,wy,wz,T,a_std   every retained block
//   [ACC] recapture hold=H ...       a hold was dropped and captured again
//   [ACCFIT] final ... b=(..) ...    the device's fit
//
// Raw replay drives the same AccelCalProcedure (capture + fit) from the logged
// samples; block replay re-runs only the fit from the retained blocks.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "imu_calibrate/AccelCalCapture.h"

namespace accel_replay {

using Proc = imu_cal::AccelCalProcedure<340, 24>;

struct Log {
  std::vector<std::string> lines;
  bool has_raw = false, has_blocks = false;
  bool accel_only = false;
  imu_cal::AccelThermalPrior prior;
  bool have_device_fit = false;
  double dev_b[3] = {0, 0, 0};
};

inline bool startsWith(const std::string& s, const char* p) { return s.rfind(p, 0) == 0; }

// Accepts lines with any prefix before the tag (serial monitors add timestamps).
inline Log parse(const std::vector<std::string>& raw) {
  Log L;
  const char* tags[] = {"[ACCMODE]", "[ACCPRIOR]", "[ACCGYRO]", "[ACCPREP]", "[ACCRAW]", "[ACCBLK]", "[ACC] recapture", "[ACCFIT] final"};
  for (const std::string& r : raw) {
    size_t pos = std::string::npos;
    for (const char* t : tags) { const size_t p = r.find(t); if (p != std::string::npos) { pos = p; break; } }
    if (pos == std::string::npos) continue;
    std::string l = r.substr(pos);
    while (!l.empty() && (l.back() == '\r' || l.back() == '\n')) l.pop_back();
    if (startsWith(l, "[ACCRAW]")) L.has_raw = true;
    if (startsWith(l, "[ACCBLK]")) L.has_blocks = true;
    if (startsWith(l, "[ACCMODE]")) L.accel_only = (l.find("accel_only") != std::string::npos);
    if (startsWith(l, "[ACCPRIOR]")) {
      int v = 0; double k[3], a, b, c, d;
      if (sscanf(l.c_str() + 10, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &v, &k[0], &k[1], &k[2], &a, &b, &c, &d) == 8) {
        L.prior.valid = v != 0;
        for (int j = 0; j < 3; ++j) L.prior.k[j] = k[j];
        L.prior.k_temp_lo = a; L.prior.k_temp_hi = b; L.prior.clamp_lo = c; L.prior.clamp_hi = d;
      }
    }
    if (startsWith(l, "[ACCFIT] final") && l.find(" b=(") != std::string::npos) {
      const char* p = strstr(l.c_str(), " b=(");
      if (sscanf(p, " b=(%lf,%lf,%lf)", &L.dev_b[0], &L.dev_b[1], &L.dev_b[2]) == 3) L.have_device_fit = true;
    }
    L.lines.push_back(l);
  }
  return L;
}

// Retained blocks as logged, with recaptured holds removed the way the
// procedure removes them.
inline std::vector<imu_cal::AccelObs> blocks(const Log& L) {
  std::vector<imu_cal::AccelObs> obs;
  for (const std::string& l : L.lines) {
    if (startsWith(l, "[ACC] recapture")) {
      int h = -1;
      if (sscanf(l.c_str(), "[ACC] recapture hold=%d", &h) == 1 && h >= 0) {
        std::vector<imu_cal::AccelObs> keep;
        for (imu_cal::AccelObs o : obs) {
          if (o.hold == h) continue;
          if (o.hold > h) o.hold--;
          keep.push_back(o);
        }
        obs.swap(keep);
      }
      continue;
    }
    if (!startsWith(l, "[ACCBLK]")) continue;
    imu_cal::AccelObs o{};
    int hold, role; unsigned long t; unsigned n;
    double v[8];
    if (sscanf(l.c_str() + 9, "%d,%d,%lu,%u,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &hold, &role, &t, &n, &v[0], &v[1],
               &v[2], &v[3], &v[4], &v[5], &v[6], &v[7]) != 12) continue;
    o.hold = (uint8_t)hold; o.role = (uint8_t)role; o.t_ms = (uint32_t)t; o.n = (uint16_t)n;
    for (int d = 0; d < 3; ++d) { o.a[d] = (float)v[d]; o.w[d] = (float)v[3 + d]; }
    o.tempC = (float)v[6]; o.a_std = (float)v[7];
    obs.push_back(o);
  }
  return obs;
}

inline bool fitBlocks(const Log& L, imu_cal::AccelFullFitResult& out, int* n_obs = nullptr) {
  const std::vector<imu_cal::AccelObs> obs = blocks(L);
  if (n_obs) *n_obs = (int)obs.size();
  auto fitter = std::make_unique<Proc::Fitter>();
  imu_cal::AccelFitCfg cfg;
  return fitter->fit(obs.data(), (int)obs.size(), cfg, L.prior, out, false);
}

// Feeds logged raw samples to the procedure; each hold's samples run until
// the next preparation marker.
class ReplayIo : public imu_cal::AccelCalIo {
public:
  explicit ReplayIo(const Log& L, bool verbose) : L_(L), verbose_(verbose) {}

  bool prep(const imu_cal::AccelStepView& v) override {
    // Advance to the next [ACCPREP]; process context lines on the way.
    while (i_ < L_.lines.size() && !startsWith(L_.lines[i_], "[ACCPREP]")) context_(L_.lines[i_++]);
    if (i_ >= L_.lines.size()) return false;
    int kind = -1, pose = -1, att = -1;
    sscanf(L_.lines[i_].c_str() + 10, "%d,%d,%d", &kind, &pose, &att);
    if (kind != (int)v.kind || pose != (int)v.pose || att != (int)v.attempt) {
      ++divergences;
      if (verbose_) printf("replay: divergence at prep (log %d,%d,%d vs replay %d,%d,%d)\n", kind, pose, att,
                           (int)v.kind, (int)v.pose, (int)v.attempt);
    }
    ++i_;
    starved_ = false;
    return true;
  }

  bool sample(imu_cal::AccelRawSample& s) override {
    while (i_ < L_.lines.size()) {
      const std::string& l = L_.lines[i_];
      if (startsWith(l, "[ACCPREP]")) break;
      ++i_;
      if (!startsWith(l, "[ACCRAW]")) { context_(l); continue; }
      unsigned long t; double v[7];
      if (sscanf(l.c_str() + 9, "%lu,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &t, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6]) != 8) continue;
      s.t_us = (uint32_t)t;
      s.a = Eigen::Vector3f((float)v[0], (float)v[1], (float)v[2]);
      s.w = Eigen::Vector3f((float)v[3], (float)v[4], (float)v[5]);
      s.tempC = (float)v[6];
      now_ms_ = (uint32_t)(t / 1000u);
      return true;
    }
    starved_ = true;
    return false;
  }

  // Once the hold's samples are exhausted, time jumps past stuck detection.
  uint32_t nowMs() override { return starved_ ? now_ms_ + 1000000u : now_ms_; }
  void capture(const imu_cal::AccelStepView&, const imu_cal::AccelHoldView&) override {}
  void holdOk(const imu_cal::AccelStepView&) override {}
  bool holdRetry(const imu_cal::AccelStepView&, const char*) override { return true; }
  bool runFit(imu_cal::AccelFitJob& job, const char*) override { job.run(); return true; }
  void log(const char* line) override { if (verbose_) printf("%s\n", line); }
  void idle() override {}

  // Processes context lines up to the next hold (between stages).
  void advanceContext() {
    while (i_ < L_.lines.size() && !startsWith(L_.lines[i_], "[ACCPREP]")) context_(L_.lines[i_++]);
  }

  bool gyro_pending = false;
  Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
  bool gyro_initial_set = false;
  Eigen::Vector3f gyro_initial = Eigen::Vector3f::Zero();
  int divergences = 0;

private:
  void context_(const std::string& l) {
    if (startsWith(l, "[ACCGYRO]")) {
      double x, y, z;
      if (sscanf(l.c_str() + 10, "%lf,%lf,%lf", &x, &y, &z) == 3) {
        gyro = Eigen::Vector3f((float)x, (float)y, (float)z);
        gyro_pending = true;
      }
    }
  }

  const Log& L_;
  bool verbose_;
  size_t i_ = 0;
  uint32_t now_ms_ = 0;
  bool starved_ = false;
};

// Raw replay through the procedure. Returns true when the final fit succeeds.
inline bool replayRaw(const Log& L, imu_cal::AccelFullFitResult& out, int& divergences, bool verbose) {
  auto proc = std::make_unique<Proc>();
  ReplayIo io(L, verbose);
  // Initial stationary gyro level (logged before the first hold).
  Eigen::Vector3f g0 = Eigen::Vector3f::Zero();
  bool g0_valid = false;
  for (const std::string& l : L.lines) {
    if (startsWith(l, "[ACCPREP]")) break;
    double x, y, z, v = 0;
    if (startsWith(l, "[ACCGYRO]") && sscanf(l.c_str() + 10, "%lf,%lf,%lf,%lf", &x, &y, &z, &v) >= 3) {
      g0 = Eigen::Vector3f((float)x, (float)y, (float)z);
      g0_valid = (v != 0);
    }
  }
  proc->begin(imu_cal::AccelCaptureCfg{}, imu_cal::AccelFitCfg{}, L.prior, g0, g0_valid);
  io.advanceContext();
  io.gyro_pending = false;
  bool ok = proc->runMainStage(io);
  if (ok) {
    io.advanceContext();
    if (io.gyro_pending && !L.accel_only) proc->setGyroReference(io.gyro);
    ok = proc->runRecheckStage(io) && proc->runFinalFit(io);
  }
  divergences = io.divergences;
  out = proc->result();
  return ok && out.ok;
}

}  // namespace accel_replay
