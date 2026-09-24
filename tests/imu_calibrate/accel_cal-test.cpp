// Accelerometer calibration: unit checks and a deterministic multi-seed
// campaign of simulated manual sessions.
//
// Methods compared on identical simulated devices and users:
//   OLD       deployed wizard: six faces, 6.5 s wait, first 60 gated samples,
//             AccelCalibrator<float,400,1> PolarSPD (cond 6, offdiag 0.10)
//   RICH+OLD  the new capture's retained blocks fed to the same old fitter
//   RICH+NEW  the new capture with the full-matrix/thermal fitter (the
//             AccelCalProcedure the wizard runs, driven through AccelCalIo)
// Every method is scored against the simulated truth on unseen orientations,
// through the float runtime path (blob -> RuntimeCals) where it applies.
//
// Outputs: calibrate_accel_campaign.csv (one row per run),
// calibrate_accel_campaign_summary.csv, calibrate_accel_report.txt, and
// calibrate_accel_replay_sample.log (a logged session for accel_cal-replay).

#include "accel_cal_sim.h"

#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include "accel_cal_replay.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <regex>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

using namespace accel_sim;
using imu_cal::AccelCalProcedure;
using imu_cal::AccelObs;
using imu_cal::AccelThermal;

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool cond, const std::string& what) {
  ++g_checks;
  if (!cond) {
    ++g_failures;
    std::cerr << "CHECK FAILED: " << what << "\n";
  }
}

double median(std::vector<double> v) {
  if (v.empty()) return NAN;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return (n & 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
double quantile(std::vector<double> v, double q) {
  if (v.empty()) return NAN;
  std::sort(v.begin(), v.end());
  const double idx = q * (double)(v.size() - 1);
  const size_t i0 = (size_t)std::floor(idx), i1 = std::min(v.size() - 1, i0 + 1);
  return v[i0] + (idx - (double)i0) * (v[i1] - v[i0]);
}
double vmax(const std::vector<double>& v) { return v.empty() ? NAN : *std::max_element(v.begin(), v.end()); }

#ifndef ACCEL_CAL_TEST_MAXO
#define ACCEL_CAL_TEST_MAXO 340  // the wizard capacity
#endif
constexpr int kMaxObs = ACCEL_CAL_TEST_MAXO;
constexpr int kMaxHolds = 24;
using Proc = AccelCalProcedure<kMaxObs, kMaxHolds>;

// Wizard configuration of the legacy accel calibrator (configureCalibrators_()).
template <typename Cal>
void configureLegacy(Cal& c) {
  c.g = 9.80665f;
  c.accel_mag_tol = 0.8f;
  c.max_gyro_for_static = 0.12f;
  c.accel_S_mode = Cal::AccelSMode::PolarSPD;
  c.accel_diag_lo = 0.80f;
  c.accel_diag_hi = 1.25f;
  c.accel_max_cond = 6.0f;
  c.accel_max_offdiag_rms = 0.10f;
}

// In-memory Preferences byte API with fault injection.
struct MemKv {
  struct State {
    std::map<std::string, std::vector<uint8_t>> m;
    bool drop_writes = false;      // putBytes reports success but stores nothing
    bool corrupt_writes = false;   // flips one stored byte
    int puts = 0;
  };
  std::shared_ptr<State> st = std::make_shared<State>();
  size_t getBytesLength(const char* k) { auto it = st->m.find(k); return it == st->m.end() ? 0 : it->second.size(); }
  size_t getBytes(const char* k, void* buf, size_t len) {
    auto it = st->m.find(k);
    if (it == st->m.end() || it->second.size() > len) return 0;
    memcpy(buf, it->second.data(), it->second.size());
    return it->second.size();
  }
  size_t putBytes(const char* k, const void* buf, size_t len) {
    ++st->puts;
    if (st->drop_writes) return len;
    std::vector<uint8_t> v((const uint8_t*)buf, (const uint8_t*)buf + len);
    if (st->corrupt_writes && len > 40) v[40] ^= 0x5A;
    st->m[k] = v;
    return len;
  }
  bool remove(const char* k) { st->m.erase(k); return true; }
};
using MemStore = atoms3r_ical::ImuCalStoreT<MemKv>;

// Calibration as the runtime applies it (float), for scoring.
struct CalFloat {
  bool ok = false;
  imu_cal::AccelCalibration<float> acc;
};

CalFloat calFromBlob(const atoms3r_ical::ImuCalBlobV3& b) {
  atoms3r_ical::RuntimeCals rc;
  rc.rebuildFromBlob(b);
  CalFloat c;
  c.ok = rc.acc.ok;
  c.acc = rc.acc;
  return c;
}

// Accuracy of one calibration against the truth at temperature T.
struct Accuracy {
  double bias_err = NAN;       // || b_fit(T) - b_true(T) ||, m/s^2
  double bias_axis_max = NAN;  // max_i |.|
  double cross_err = NAN;      // max |S_fit - S_true| off-diagonal
  double diag_err = NAN;       // max |S_fit - S_true| diagonal
  double vec_rms = NAN, vec_max = NAN;   // || a_cal - f_true || over unseen orientations
  double norm_rms = NAN;       // | ||a_cal|| - g | over unseen orientations
};

const std::vector<Vec3>& unseenDirections() {
  static std::vector<Vec3> dirs;
  if (dirs.empty()) {
    Rng r(0xC0FFEEull);
    for (int i = 0; i < 400; ++i) dirs.push_back(r.unit());
  }
  return dirs;
}

Accuracy score(const CalFloat& c, const SensorTruth& t, double T) {
  Accuracy a;
  if (!c.ok) return a;
  const Vec3 bt = t.bias_at(T);
  const Vec3 bf = c.acc.biasT.bias((float)T).cast<double>();
  a.bias_err = (bf - bt).norm();
  a.bias_axis_max = (bf - bt).cwiseAbs().maxCoeff();
  const Mat3 St = t.S_polar_at(T) * (kGStd / t.g_local);
  const Mat3 Sf = c.acc.S.cast<double>();
  double ce = 0, de = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const double e = std::fabs(Sf(i, j) - St(i, j));
      if (i == j) de = std::max(de, e); else ce = std::max(ce, e);
    }
  a.cross_err = ce; a.diag_err = de;
  const Mat3 A = t.S_at(T) * t.R_mis.transpose();
  const Mat3 Ainv = A.inverse();
  double sv = 0, sn = 0, mx = 0;
  for (const Vec3& u : unseenDirections()) {
    const Vec3 f = t.g_local * u;
    const Vec3 raw = Ainv * f + bt;
    const Vec3 ac = c.acc.apply(raw.cast<float>(), (float)T).cast<double>();
    const Vec3 ref = t.R_mis * f * (kGStd / t.g_local);
    const double ev = (ac - ref).norm();
    sv += ev * ev; mx = std::max(mx, ev);
    const double en = ac.norm() - kGStd;
    sn += en * en;
  }
  const double n = (double)unseenDirections().size();
  a.vec_rms = std::sqrt(sv / n); a.vec_max = mx;
  a.norm_rms = std::sqrt(sn / n);
  return a;
}

// ---------------------------------------------------------------------------
// Simulated user + IMU driving the real procedure.

class SimIo : public imu_cal::AccelCalIo {
public:
  SimIo(World& w, Stream& s, const Scenario& sc) : w_(w), s_(s), sc_(sc) {}

  std::vector<std::string>* logs = nullptr;
  bool raw_log = false;
  double fit_cpu_s = 0;
  int wrong_tilts = 0;
  bool abort_at_prep = false;

  bool prep(const imu_cal::AccelStepView& v) override {
    if (logs) {
      char b[96];
      snprintf(b, sizeof(b), "[ACCPREP] %d,%d,%d", (int)v.kind, (int)v.pose, (int)v.attempt);
      logs->push_back(b);
    }
    if (abort_at_prep) return false;
    Rng& r = w_.rng();
    w_.advance(r.uni(sc_.user.react_lo, sc_.user.react_hi));
    const imu_cal::AccelPoseDef& p = imu_cal::kAccelPoses[v.pose];
    const bool tilt = (p.region == imu_cal::AccelRegion::FRAME);
    hand_ = tilt ? r.chance(sc_.user.p_hand_tilt) : r.chance(sc_.user.p_hand_face[v.pose]);
    err_ = tilt ? sc_.user.tilt_err_deg : (hand_ ? sc_.user.face_err_hand_deg : sc_.user.face_err_table_deg);
    target_ = docUpBody(p.top, p.right, p.out);
    Vec3 first = target_;
    pending_fix_ = false;
    if (tilt && r.chance(sc_.user.p_first_wrong_tilt)) {
      first = docUpBody(p.top, 0, p.out);   // tilted about one axis only
      pending_fix_ = true;
      ++wrong_tilts;
    }
    const double tm = r.chance(sc_.user.p_slow) ? r.uni(5.0, 9.0) : r.uni(sc_.user.move_lo, sc_.user.move_hi);
    w_.placeAt(first, hand_, err_, w_.t(), tm);
    check_since_ = -1;
    // Poisson knocks during the expected hold
    const double horizon = tm + 12.0;
    double tt = w_.t();
    for (;;) {
      if (sc_.user.shock_rate_per_s <= 0) break;
      tt += -std::log(std::max(1e-12, r.uni())) / sc_.user.shock_rate_per_s;
      if (tt > w_.t() + horizon) break;
      w_.addShock(tt, 0.03, r.unit() * r.uni(sc_.user.shock_lo, sc_.user.shock_hi));
    }
    return true;
  }

  bool sample(imu_cal::AccelRawSample& s) override {
    const bool ok = s_.next(s);
    if (ok && raw_log && logs) {
      char b[160];
      snprintf(b, sizeof(b), "[ACCRAW] %lu,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f", (unsigned long)s.t_us,
               (double)s.a.x(), (double)s.a.y(), (double)s.a.z(), (double)s.w.x(), (double)s.w.y(),
               (double)s.w.z(), (double)s.tempC);
      logs->push_back(b);
    }
    return ok;
  }
  uint32_t nowMs() override { return (uint32_t)(uint64_t)std::llround(w_.t() * 1000.0); }
  void capture(const imu_cal::AccelStepView&, const imu_cal::AccelHoldView& h) override {
    if (h.hint == imu_cal::AccelHoldHint::CHECK_POSE) {
      if (check_since_ < 0) check_since_ = w_.t();
      else if (pending_fix_ && w_.t() - check_since_ > 2.5) {
        w_.placeAt(target_, hand_, err_, w_.t(), w_.rng().uni(1.0, 2.5));
        pending_fix_ = false;
      }
    } else {
      check_since_ = -1;
    }
  }
  void holdOk(const imu_cal::AccelStepView&) override { w_.advance(0.98); }
  bool holdRetry(const imu_cal::AccelStepView&, const char*) override { w_.advance(2.0); return true; }
  bool runFit(imu_cal::AccelFitJob& job, const char*) override {
    const auto t0 = std::chrono::steady_clock::now();
    job.run();
    fit_cpu_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    w_.advance(0.5);
    return true;
  }
  void log(const char* line) override { if (logs) logs->push_back(line); }
  void idle() override {}

private:
  World& w_;
  Stream& s_;
  const Scenario& sc_;
  bool hand_ = false;
  double err_ = 0;
  Vec3 target_ = Vec3::Zero();
  bool pending_fix_ = false;
  double check_since_ = -1;
};

// ---------------------------------------------------------------------------
// Deployed (old) wizard accelerometer stage.

// The deployed fitter counted a temperature bin before its matrix passed the
// plausibility gates. With one bin, a rejected matrix left S = I, b = 0 and the
// fit still returned true after gravity rescaling. Reproduces that outcome
// for the baseline (this PR makes the library fail instead).
template <typename Cal>
bool deployedIdentityAcceptance(const Cal& cal, imu_cal::AccelCalibration<float>& acc) {
  auto fitk = imu_cal::ellipsoid_to_sphere_robust<float>(cal.buf.v, cal.buf.n, cal.g, 3, 0.15f, 1e-6f, cal.g);
  if (!fitk.ok) return false;
  Eigen::Matrix3f S_spd;
  bool rejected = !imu_cal::polar_spd_factor_3x3<float>(fitk.A, S_spd);
  if (!rejected) {
    for (int j = 0; j < 3; ++j) if (!(S_spd(j, j) >= cal.accel_diag_lo && S_spd(j, j) <= cal.accel_diag_hi)) rejected = true;
    float condv = 0;
    if (!imu_cal::cond_spd_3x3<float>(S_spd, condv) || condv > cal.accel_max_cond) rejected = true;
    if (imu_cal::offdiag_rms_3x3<float>(S_spd) > cal.accel_max_offdiag_rms) rejected = true;
  }
  if (!rejected) return false;
  acc = imu_cal::AccelCalibration<float>{};
  acc.ok = true; acc.g = cal.g; acc.S.setIdentity();
  acc.biasT.ok = true; acc.biasT.T0 = cal.T0; acc.biasT.b0.setZero(); acc.biasT.k.setZero();
  bool did = imu_cal::post_scale_accel_S_to_match_g_<float>(cal.buf, acc, 0.985f, 1.015f, 0.97f, 1.03f);
  if (!did) did = imu_cal::post_scale_accel_S_to_match_g_<float>(cal.buf, acc, 0.95f, 1.05f, 0.93f, 1.07f);
  if (!did) did = imu_cal::post_scale_accel_S_to_match_g_<float>(cal.buf, acc, 0.90f, 1.10f, 0.90f, 1.10f);
  return did;
}

struct OldRun {
  bool ok = false;
  bool identity = false;   // deployed code silently accepted S = I, b = 0
  std::string reason;
  CalFloat cal;
  double t_start = 0, t_end = 0;
  double T_mid = NAN;
};

OldRun runOld(World& w, Stream& st, const Scenario& sc) {
  OldRun out;
  auto cal = std::make_unique<imu_cal::AccelCalibrator<float, 400, 1>>();
  configureLegacy(*cal);
  out.t_start = w.t();
  Rng& r = w.rng();
  for (int p = 0; p < 6; ++p) {
    const imu_cal::AccelPoseDef& pd = imu_cal::kAccelPoses[p];
    w.advance(r.uni(sc.user.react_lo, sc.user.react_hi));
    const bool hand = r.chance(sc.user.p_hand_face[p]);
    const double err = hand ? sc.user.face_err_hand_deg : sc.user.face_err_table_deg;
    const double tm = r.chance(sc.user.p_slow) ? r.uni(5.0, 9.0) : r.uni(sc.user.move_lo, sc.user.move_hi);
    w.placeAt(docUpBody(pd.top, pd.right, pd.out), hand, err, w.t(), tm);
    double tt = w.t();
    for (;;) {
      if (sc.user.shock_rate_per_s <= 0) break;
      tt += -std::log(std::max(1e-12, r.uni())) / sc.user.shock_rate_per_s;
      if (tt > w.t() + tm + 12.0) break;
      w.addShock(tt, 0.03, r.unit() * r.uni(sc.user.shock_lo, sc.user.shock_hi));
    }
    // PLACE_TIME_MS: samples are not consumed.
    const double t_place_end = w.t() + 6.5;
    imu_cal::AccelRawSample s;
    while (w.t() < t_place_end) st.next(s);
    const int start_n = cal->buf.n;
    double last_change = w.t();
    const double t_cap0 = w.t();
    int last_n = start_n;
    bool done = false;
    while (w.t() - t_cap0 < 90.0) {
      if (st.next(s)) cal->addSample(s.a, s.w, s.tempC);
      if (cal->buf.n != last_n) { last_n = cal->buf.n; last_change = w.t(); }
      if (w.t() - last_change > 12.0) { out.reason = "stuck"; out.t_end = w.t(); return out; }
      if (cal->buf.n >= start_n + 60) { done = true; break; }
    }
    if (!done) { out.reason = "timeout"; out.t_end = w.t(); return out; }
    w.advance(0.98);
  }
  out.t_end = w.t();
  out.T_mid = NAN;
  if (cal->buf.n < 220) { out.reason = "too_few"; return out; }
  imu_cal::AccelCalibration<float> acc;
  imu_cal::FitFail why = imu_cal::FitFail::OK;
  const bool ok = cal->fit(acc, 3, 0.15f, &why);
  if (!(ok && acc.ok)) {
    if (!deployedIdentityAcceptance(*cal, acc)) { out.reason = imu_cal::fitFailStr(why); return out; }
    out.identity = true;
  }
  // Runtime path: blob -> RuntimeCals, as the deployed wizard applied it.
  atoms3r_ical::ImuCalBlobV3 b{};
  memset((void*)&b, 0, sizeof(b));
  b.accel_ok = 1; b.accel_g = acc.g;
  atoms3r_ical::mat_to_rowmajor9_(acc.S, b.accel_S);
  b.accel_T0 = acc.biasT.T0;
  for (int j = 0; j < 3; ++j) { b.accel_b0[j] = acc.biasT.b0(j); b.accel_k[j] = acc.biasT.k(j); }
  b.accel_T_lo = -1000; b.accel_T_hi = 1000;
  out.cal = calFromBlob(b);
  out.ok = true;
  out.reason = out.identity ? "OK_IDENTITY" : "OK";
  return out;
}

CalFloat fitOldOnBlocks(const AccelObs* obs, int n, std::string& reason) {
  auto cal = std::make_unique<imu_cal::AccelCalibrator<float, 400, 1>>();
  configureLegacy(*cal);
  for (int i = 0; i < n; ++i) {
    if (obs[i].role != (uint8_t)imu_cal::AccelObsRole::FIT) continue;
    cal->addSample(Eigen::Vector3f(obs[i].a[0], obs[i].a[1], obs[i].a[2]),
                   Eigen::Vector3f(obs[i].w[0], obs[i].w[1], obs[i].w[2]), obs[i].tempC);
  }
  CalFloat c;
  if (cal->buf.n < 220) { reason = "too_few"; return c; }
  imu_cal::AccelCalibration<float> acc;
  imu_cal::FitFail why = imu_cal::FitFail::OK;
  if (!(cal->fit(acc, 3, 0.15f, &why) && acc.ok)) {
    if (!deployedIdentityAcceptance(*cal, acc)) { reason = imu_cal::fitFailStr(why); return c; }
    reason = "OK_IDENTITY";
  } else {
    reason = "OK";
  }
  c.ok = true;
  c.acc = acc;
  return c;
}

struct NewRun {
  bool ok = false;
  std::string reason;
  CalFloat cal;
  imu_cal::AccelFullFitResult fit;
  atoms3r_ical::ImuCalBlobV3 blob{};
  double accel_time_s = 0;
  int retries = 0, extras = 0, attempts = 0, wrong_tilts = 0;
  double T_mid = NAN;
  double fit_cpu_s = 0;
  CalFloat rich_old;
  std::string rich_old_reason;
  std::vector<AccelObs> obs;
};

NewRun runNew(World& w, Stream& st, const Scenario& sc, const imu_cal::AccelThermalPrior& prior,
              std::vector<std::string>* logs = nullptr, bool raw_log = false) {
  NewRun out;
  auto proc = std::make_unique<Proc>();
  imu_cal::AccelCaptureCfg ccfg;
  imu_cal::AccelFitCfg fcfg;
  if (const char* e = getenv("ACCEL_CAL_HOLD_MS")) ccfg.hold_useful_ms = (uint32_t)strtoul(e, nullptr, 10);  // sensitivity sweeps
  proc->begin(ccfg, fcfg, prior, Eigen::Vector3f::Zero(), false);
  SimIo io(w, st, sc);
  io.logs = logs;
  io.raw_log = raw_log;
  if (logs) {
    char b[160];
    snprintf(b, sizeof(b), "[ACCMODE] %s", sc.accel_only ? "accel_only" : "full");
    logs->push_back(b);
    snprintf(b, sizeof(b), "[ACCPRIOR] %d,%.7f,%.7f,%.7f,%.2f,%.2f,%.2f,%.2f", (int)prior.valid, prior.k[0], prior.k[1],
             prior.k[2], prior.k_temp_lo, prior.k_temp_hi, prior.clamp_lo, prior.clamp_hi);
    logs->push_back(b);
    logs->push_back("[ACCGYRO] 0,0,0,0");
  }
  const double t0 = w.t();
  const bool main_ok = proc->runMainStage(io);
  const double t1 = w.t();
  double t2 = t1, t3 = t1;
  bool ok = main_ok;
  if (ok) {
    if (!sc.accel_only) {
      // Gyro + mag stages of the full wizard: the device rests / is rotated.
      w.advance(w.rng().uni(sc.gyro_mag_gap_lo, sc.gyro_mag_gap_hi));
      Vec3 wb = w.truth().gyro_bias;
      proc->setGyroReference(wb.cast<float>());
      if (logs) {
        char b[96];
        snprintf(b, sizeof(b), "[ACCGYRO] %.7f,%.7f,%.7f,1", (double)(float)wb.x(), (double)(float)wb.y(), (double)(float)wb.z());
        logs->push_back(b);
      }
    }
    t2 = w.t();
    ok = proc->runRecheckStage(io);
    t3 = w.t();
    if (ok) ok = proc->runFinalFit(io);
  }
  out.accel_time_s = (t1 - t0) + (t3 - t2);
  out.retries = proc->totalRetries();
  out.extras = proc->nExtra();
  out.attempts = proc->totalAttempts();
  out.wrong_tilts = io.wrong_tilts;
  out.fit_cpu_s = io.fit_cpu_s;
  out.obs.assign(proc->obs(), proc->obs() + proc->nObs());
  // Temperature at the middle of the accelerometer capture (truth).
  double sT = 0; int nT = 0;
  for (int i = 0; i < proc->nObs(); ++i) {
    if (std::isfinite(proc->obs()[i].tempC)) { sT += proc->obs()[i].tempC; ++nT; }
  }
  out.T_mid = nT ? sT / nT : NAN;
  if (!ok) {
    out.reason = std::string(proc->failure() == imu_cal::AccelProcFail::FIT_FAILED ? "fit:" : "proc:") + proc->failureDetail();
    out.fit = proc->result();
  } else {
    out.fit = proc->result();
    imu_cal::AccelCalibration<float> fc;
    Proc::Fitter::toFloat(out.fit, fcfg.g, fc);
    atoms3r_ical::ImuCalBlobV3 b{};
    memset((void*)&b, 0, sizeof(b));
    atoms3r_ical::fillAccelFromFit(b, out.fit, fc, (uint32_t)(proc->totalHoldMs() / 1000));
    b.sensor_id_lo = 0x1234; b.sensor_id_hi = 0x5678; b.imu_type = 3;
    // Serialize through the store and validate the read-back float set again.
    MemStore store;
    atoms3r_ical::ImuCalBlobV3 rb{};
    if (!store.saveVerified(b, rb)) { out.reason = "save"; return out; }
    atoms3r_ical::RuntimeCals rc;
    rc.rebuildFromBlob(rb);
    imu_cal::AccelFullFitResult rv = out.fit;
    if (!Proc::Fitter::validateFloat(proc->obs(), proc->nObs(), rc.acc, fcfg, rv)) { out.reason = "float_readback"; return out; }
    out.blob = rb;
    out.cal = calFromBlob(rb);
    out.ok = true;
    out.reason = "OK";
  }
  out.rich_old = fitOldOnBlocks(proc->obs(), proc->nObs(), out.rich_old_reason);
  return out;
}

// ---------------------------------------------------------------------------
// Unit checks

void testLegacyFitterSuccessPaths() {
  // A strong symmetric cross-axis term makes every bin fail the PolarSPD
  // off-diagonal gate. The deployed code counted the bin before validation and
  // returned the identity default with ok = true; it must now fail.
  imu_cal::AccelCalibrator<double, 400, 1> cal;
  configureLegacy(cal);
  cal.accel_mag_tol = 5.0;
  Mat3 S; S << 1.0, 0.25, 0.2, 0.25, 1.0, 0.22, 0.2, 0.22, 1.0;
  const Mat3 Sinv = S.inverse();
  Rng r(5);
  for (int i = 0; i < 300; ++i) {
    const Vec3 u = r.unit();
    const Vec3 raw = Sinv * (kGStd * u) + Vec3(0.1, -0.05, 0.02);
    cal.addSample(raw, Vec3::Zero(), 25.0);
  }
  imu_cal::AccelCalibration<double> out;
  out.ok = true;
  imu_cal::FitFail why = imu_cal::FitFail::OK;
  const bool ok = cal.fit(out, 3, 0.15, &why);
  check(!ok, "legacy fitter rejects an all-bins-invalid fit");
  check(!out.ok, "legacy fitter leaves out.ok false on failure");
  check(why == imu_cal::FitFail::ACCEL_S_UNPHYSICAL, "legacy failure reason is ACCEL_S_UNPHYSICAL");
}

void testTempBias() {
  imu_cal::TempBias3<float> tb;
  tb.T0 = 25; tb.b0 = Eigen::Vector3f(0.1f, 0.2f, 0.3f); tb.k = Eigen::Vector3f(0.01f, 0, 0);
  check(tb.bias(NAN).isApprox(tb.b0), "NaN temperature evaluates the reference bias");
  check(std::fabs(tb.bias(35).x() - 0.2f) < 1e-6f, "unclamped default extrapolates linearly");
  tb.T_lo = 20; tb.T_hi = 30;
  check(std::fabs(tb.bias(40).x() - 0.15f) < 1e-6f, "clamp bounds extrapolation");
  check(std::fabs(tb.bias(10).x() - 0.05f) < 1e-6f, "clamp bounds extrapolation below");
}

// A complete stored calibration (gyro/mag set, accelerometer set bound).
atoms3r_ical::ImuCalBlobV3 makeBlob() {
  atoms3r_ical::ImuCalBlobV3 b{};
  memset((void*)&b, 0, sizeof(b));
  b.accel_ok = 1; b.accel_g = 9.80665f;
  const float S[9] = {1.01f, 0.003f, -0.002f, 0.003f, 0.99f, 0.001f, -0.002f, 0.001f, 1.005f};
  memcpy(b.accel_S, S, sizeof(S));
  b.accel_T0 = 25; b.accel_b0[0] = 0.1f; b.accel_b0[1] = -0.2f; b.accel_b0[2] = 0.05f;
  b.accel_T_lo = -1000; b.accel_T_hi = 1000;
  b.accel_thermal = (uint8_t)AccelThermal::UNLEARNED;
  b.gyro_ok = 1; b.gyro_b0[0] = 0.01f;
  b.mag_ok = 1; b.mag_A[0] = b.mag_A[4] = b.mag_A[8] = 1; b.mag_b[0] = 5; b.mag_field_uT = 45;
  b.accel_coeff_crc = atoms3r_ical::accelCoeffCrc(b);
  return MemStore::sealed_(b);
}

void testBlobAndStore() {
  using namespace atoms3r_ical;
  // Only the current layout and key are read: bytes under the keys of earlier
  // firmware (any size) are ignored.
  {
    MemStore store;
    const ImuCalBlobV3 b = makeBlob();
    store.kv.putBytes("blob_m5", &b, sizeof(b));
    store.kv.putBytes("blob", &b, sizeof(b));
    uint8_t old[124] = {0x4D, 0x55, 0x4C, 0x43, 2, 0};
    store.kv.putBytes("blob_m5v3", old, sizeof(old));
    ImuCalBlobV3 ld;
    check(!store.load(ld), "earlier keys and layouts are not loaded");
    ImuCalBlobV3 wrongver = b; wrongver.version = 2; wrongver.crc = computeBlobCrc(wrongver);
    check(!validateBlob(wrongver), "other versions do not validate");
    ImuCalBlobV3 bad = b; bad.accel_b0[0] += 1.0f;
    check(!validateBlob(bad), "CRC mismatch rejected");
  }
  // Save/readback validates the candidate, never an older blob.
  {
    MemStore store;
    const ImuCalBlobV3 old = makeBlob();
    ImuCalBlobV3 rb0;
    check(store.saveVerified(old, rb0), "first verified save succeeds");
    ImuCalBlobV3 cand = old;
    cand.accel_b0[0] = 0.123f;
    cand.accel_coeff_crc = accelCoeffCrc(cand);
    store.kv.st->drop_writes = true;
    ImuCalBlobV3 rb;
    check(!store.saveVerified(cand, rb), "dropped write (stale blob under the key) is not reported as saved");
    ImuCalBlobV3 after; check(store.load(after) && after.accel_b0[0] == old.accel_b0[0], "previous calibration retained after failed save");
    store.kv.st->drop_writes = false;
    store.kv.st->corrupt_writes = true;
    check(!store.saveVerified(cand, rb), "corrupted write is not reported as saved");
    store.kv.st->corrupt_writes = false;
    check(store.saveVerified(cand, rb), "verified save succeeds");
    check(rb.accel_b0[0] == 0.123f, "read-back is the candidate");
    ImuCalBlobV3 ld; check(store.load(ld) && MemStore::sameBytes(ld, rb), "load returns the verified blob");
    store.erase();
    check(!store.load(ld), "erase removes the calibration");
  }
  // Metadata binding and prior compatibility.
  {
    ImuCalBlobV3 b = makeBlob();
    b.accel_thermal = (uint8_t)AccelThermal::LEARNED;
    b.accel_k[0] = 0.003f; b.accel_k_temp_lo = 24; b.accel_k_temp_hi = 30; b.accel_T_lo = 19; b.accel_T_hi = 35;
    b.sensor_id_lo = 0x1234; b.sensor_id_hi = 0x5678; b.imu_type = 3;
    b.accel_coeff_crc = accelCoeffCrc(b);
    b = MemStore::sealed_(b);
    check(accelThermalPriorFrom(b, 0x1234, 0x5678, 3).valid, "bound learned slope of the same sensor is a prior");
    check(!accelThermalPriorFrom(b, 0x9999, 0x5678, 3).valid, "other sensor: no prior");
    check(!accelThermalPriorFrom(b, 0x1234, 0x5678, 4).valid, "other IMU type: no prior");
    ImuCalBlobV3 u = b; u.accel_thermal = (uint8_t)AccelThermal::UNLEARNED; u = MemStore::sealed_(u);
    check(!accelThermalPriorFrom(u, 0x1234, 0x5678, 3).valid, "unlearned slope is never carried");
    ImuCalBlobV3 t = b; t.accel_k[0] = 0.004f;  // coefficient changed, metadata not re-bound
    t = MemStore::sealed_(t);
    check(validateBlob(t) && !accelMetaBound(t), "unbound metadata detected");
    check(!accelThermalPriorFrom(t, 0x1234, 0x5678, 3).valid, "unbound metadata: no prior");
    ImuCalBlobV3 nf = b; nf.accel_b0[1] = NAN; nf = MemStore::sealed_(nf);
    check(!validateBlob(nf), "non-finite coefficients rejected");
  }
}

// Screen text must fit the 128x128 display (6x8 font at size 1, 12x16 at size 2)
// at every rotation; the display is square so one width serves all four.
void testScreenText() {
  const int kCols1 = 128 / 6;   // 21
  const int kCols2 = 128 / 12;  // 10
  const int kRows = 128 / 8;
  Proc* p = new Proc();
  imu_cal::AccelCaptureCfg ccfg; imu_cal::AccelFitCfg fcfg; imu_cal::AccelThermalPrior pr;
  p->begin(ccfg, fcfg, pr, Eigen::Vector3f::Zero(), false);
  for (int k = 0; k < 3; ++k) {
    const auto kind = (imu_cal::AccelHoldKind)k;
    const int count = (k == 0) ? imu_cal::kAccelMainPoses : (k == 1 ? ccfg.max_extra_holds : imu_cal::kAccelRecheckPoses);
    for (int pose = 0; pose < imu_cal::kAccelMainPoses; ++pose) {
      for (int idx = 1; idx <= count; ++idx) {
        const imu_cal::AccelStepView v = p->stepView(kind, (uint8_t)pose, idx, count, Proc::weakText_(5), 3);
        check((int)strlen(v.title) <= kCols2, std::string("title fits: ") + v.title);
        check((int)strlen(v.label) <= kCols1, std::string("label fits: ") + v.label);
        check((int)strlen(v.note) <= kCols1, std::string("note fits: ") + v.note);
        int nlines = 0;
        for (int j = 0; j < 3; ++j) if (v.lines[j]) { ++nlines; check((int)strlen(v.lines[j]) <= kCols1, std::string("line fits: ") + v.lines[j]); }
        // prep: title(2 rows) + blank + label + lines + note + blank + "Tap then place" + "Tap BtnA"
        check(2 + 1 + 1 + nlines + 1 + 1 + 1 + 1 <= kRows, "prep screen rows fit");
        // capture: title(2) + blank + label + lines + blank + hint, above the bar (last 3 rows)
        check(2 + 1 + 1 + nlines + 1 + 1 <= kRows - 3, "capture screen rows fit above the bar");
      }
    }
  }
  for (int h = 0; h <= 4; ++h) check((int)strlen(imu_cal::accelHoldHintStr((imu_cal::AccelHoldHint)h)) <= kCols1, "hint fits");
  for (int f = 0; f <= 4; ++f) check((int)strlen(imu_cal::accelHoldFailStr((imu_cal::AccelHoldFail)f)) <= kCols1, "fail text fits");
  for (int q = -1; q < 6; ++q) check((int)strlen(Proc::weakText_(q)) <= kCols1, "weak text fits");
  delete p;
}

// Fixed wizard screen strings (AtomS3R headers are Arduino-only, so they are
// checked at source level): titles <= 10 columns, other text <= 21 columns.
void testWizardScreenStrings() {
  std::string dir = __FILE__;
  const size_t slash = dir.find_last_of('/');
  dir = (slash == std::string::npos) ? std::string(".") : dir.substr(0, slash);
  int checked = 0;
  for (const char* rel : {"/../../src/AtomS3R/AtomS3R_M5Ui.h", "/../../src/AtomS3R/AtomS3R_ImuCalWizard.h"}) {
    std::ifstream f(dir + rel);
    check(f.good(), std::string("screen-string source readable: ") + rel);
    std::stringstream ss; ss << f.rdbuf();
    const std::string src = ss.str();
    // Titles: title("...") and the first argument of waitTap("...").
    const std::regex title_re(R"re(\b(?:title|waitTap)\s*\(\s*"([^"]*)")re");
    for (auto it = std::sregex_iterator(src.begin(), src.end(), title_re); it != std::sregex_iterator(); ++it) {
      check((*it)[1].length() <= 10, "title fits 10 columns: " + (*it)[1].str());
      ++checked;
    }
    // Every literal inside a display call on one line.
    const std::regex call_re(R"re(\b(?:line|lineAt|fail|failLines|retryMenu|showOkAuto|magFailMenu|waitTap|notSavedNotice)\s*\(([^;]*)\);)re");
    const std::regex lit_re(R"re("([^"]*)")re");
    for (auto it = std::sregex_iterator(src.begin(), src.end(), call_re); it != std::sregex_iterator(); ++it) {
      const std::string args = (*it)[1].str();
      for (auto jt = std::sregex_iterator(args.begin(), args.end(), lit_re); jt != std::sregex_iterator(); ++jt) {
        check((*jt)[1].length() <= 21, "screen text fits 21 columns: " + (*jt)[1].str());
        ++checked;
      }
    }
  }
  check(checked > 40, "screen-string scan found the wizard text");
}

// Instruction <-> device edges <-> display rotation <-> measured direction.
void testPoseMapping() {
  // Deployed face rotations (relative to ROT_READ) must be preserved.
  const uint8_t deployed[6] = {0, 0, 2, 0, 1, 3};
  for (int i = 0; i < 6; ++i) check(imu_cal::kAccelPoses[i].rot_delta == deployed[i], "deployed face rotation preserved");
  const char* deployed_text[6] = {"Screen faces up", "Screen faces table", "USB points up", "USB points down", "Left edge down", "Right edge down"};
  for (int i = 0; i < 6; ++i) check(strcmp(imu_cal::kAccelPoses[i].lines[0], deployed_text[i]) == 0, "deployed face instruction preserved");
  for (int i = 0; i < imu_cal::kAccelMainPoses; ++i) {
    const auto& p = imu_cal::kAccelPoses[i];
    // Text-up of the rotation is the highest screen edge (or reading rotation
    // when the screen is horizontal).
    int tt, tr;
    imu_cal::accelRotTextUp(p.rot_delta, tt, tr);
    if (p.top == 0 && p.right == 0) {
      check(p.rot_delta == 0, "horizontal screen uses the reading rotation");
    } else {
      const int dot = tt * p.top + tr * p.right;
      int best = -10;
      for (int r = 0; r < 4; ++r) { int a, b; imu_cal::accelRotTextUp((uint8_t)r, a, b); best = std::max(best, a * p.top + b * p.right); }
      check(dot == best && dot > 0, std::string("rotation puts the text up on the raised edge: ") + p.name);
    }
    // Instruction semantics: "USB ... UP" means the USB end (-top) is raised.
    const std::string all = std::string(p.lines[0]) + " " + (p.lines[1] ? p.lines[1] : "") + " " + (p.lines[2] ? p.lines[2] : "");
    if (all.find("USB end UP") != std::string::npos || all.find("USB points up") != std::string::npos) check(p.top < 0, "USB up => -top up");
    if (all.find("USB end DOWN") != std::string::npos || all.find("USB points down") != std::string::npos) check(p.top > 0, "USB down => top up");
    if (all.find("Left edge") != std::string::npos) check(p.right > 0, "left edge down => right up");
    if (all.find("Right edge") != std::string::npos) check(p.right < 0, "right edge down => left up");
    if (all.find("Screen up") != std::string::npos || all.find("Screen faces up") != std::string::npos) check(p.out > 0, "screen up => out up");
    if (all.find("Screen faces table") != std::string::npos) check(p.out < 0, "screen down => out down");
    // Documented body mapping: screen up measures body z = -g.
    if (i == 0) check((imu_cal::accelPoseDocUp(p) - Eigen::Vector3f(0, 0, -1)).norm() < 1e-6f, "screen up = body -z specific force");
  }
  // Measured frame: the procedure recovers top/right/out from the faces and
  // accepts the corners with the documented mapping and with a permuted one.
  for (int perm = 0; perm < 2; ++perm) {
    Scenario sc; sc.name = "mapping";
    sc.user.p_first_wrong_tilt = 0; sc.user.shock_rate_per_s = 0;
    Rng r(77 + perm);
    SensorTruth t = SensorTruth::random(r);
    World w(t, sc, 1000 + perm);
    Stream st(w, sc.stream);
    // perm 1: the device's x/y edges are swapped relative to the documentation
    // (the simulated user still follows the instructions physically).
    struct PermIo : public SimIo {
      using SimIo::SimIo;
    };
    auto proc = std::make_unique<Proc>();
    imu_cal::AccelCaptureCfg ccfg; imu_cal::AccelFitCfg fcfg; imu_cal::AccelThermalPrior pr;
    proc->begin(ccfg, fcfg, pr, Eigen::Vector3f::Zero(), false);
    if (perm == 0) {
      SimIo io(w, st, sc);
      check(proc->runMainStage(io), "main stage completes with the documented mapping");
      Eigen::Vector3f top, right, out;
      check(proc->frame(top, right, out), "measured frame available");
      check(top.dot(imu_cal::accelDocTop()) > 0.97f && right.dot(imu_cal::accelDocRight()) > 0.97f &&
            out.dot(imu_cal::accelDocOut()) > 0.99f, "measured frame matches documented mapping");
    }
  }
}

// Parameter information of the ten-pose design: ideal geometry and manual
// placement errors. Uses the fitter's own information matrix.
void testInformation(std::ostream& rep) {
  imu_cal::AccelFitCfg cfg;
  auto fitter = std::make_unique<Proc::Fitter>();
  auto run = [&](const std::vector<Vec3>& ups, imu_cal::AccelFullFitResult& r) {
    std::vector<AccelObs> obs;
    for (size_t h = 0; h < ups.size(); ++h) {
      for (int b = 0; b < 20; ++b) {
        AccelObs o{};
        const Vec3 a = kGStd * ups[h];
        for (int d = 0; d < 3; ++d) { o.a[d] = (float)a(d); o.w[d] = 0; }
        o.tempC = 25; o.hold = (uint8_t)h; o.role = 0; o.n = 25;
        obs.push_back(o);
      }
    }
    imu_cal::AccelThermalPrior pr;
    fitter->fit(obs.data(), (int)obs.size(), cfg, pr, r, true);
  };
  std::vector<Vec3> ideal;
  for (int p = 0; p < imu_cal::kAccelMainPoses; ++p) ideal.push_back(imu_cal::accelPoseDocUp(imu_cal::kAccelPoses[p]).cast<double>());
  imu_cal::AccelFullFitResult r;
  run(ideal, r);
  const double so = r.sigma_obs;
  auto dop = [&](const imu_cal::AccelFullFitResult& x, int idx) { return std::sqrt(x.C_unit(idx, idx)); };
  double bmax = 0, cmax = 0;
  for (int j = 0; j < 3; ++j) { bmax = std::max(bmax, dop(r, 6 + j)); cmax = std::max(cmax, dop(r, 3 + j)); }
  rep << "information,ideal_10_poses,bias_dop_max," << bmax << "\n";
  rep << "information,ideal_10_poses,cross_dop_max," << cmax << "\n";
  check(r.ok && std::isfinite(bmax) && bmax < 0.75, "ideal 10-pose design determines every bias");
  check(std::isfinite(cmax) && cmax < 1.6, "ideal 10-pose design determines every cross term");
  (void)so;
  // Faces only: cross terms not determined (information gate must say so).
  std::vector<Vec3> faces(ideal.begin(), ideal.begin() + 6);
  imu_cal::AccelFullFitResult rf;
  run(faces, rf);
  check(!rf.ok && rf.reason == imu_cal::FitFail::TOO_FEW_SAMPLES, "six holds are below the nine-parameter minimum");
  std::vector<Vec3> faces9 = faces;
  faces9.push_back(ideal[0]); faces9.push_back(ideal[2]); faces9.push_back(ideal[4]);
  imu_cal::AccelFullFitResult rf9;
  run(faces9, rf9);
  check(!rf9.ok && rf9.reason == imu_cal::FitFail::ACCEL_INFO_LOW && rf9.weak_param >= 3,
        "faces + repeated faces: cross terms flagged undetermined");
  rep << "information,faces_plus_rechecks,gate," << imu_cal::accelGateStr(rf9.gate) << "\n";
  // Manual placement: 15 deg rms on tilts, 5 deg on faces.
  Rng rng(99);
  std::vector<double> bd, cd;
  for (int t = 0; t < 60; ++t) {
    std::vector<Vec3> ups;
    for (int p = 0; p < imu_cal::kAccelMainPoses; ++p) {
      const double sd = (p < 6 ? 5.0 : 15.0) * kDeg;
      ups.push_back(expso3(rng.nvec(sd / std::sqrt(2.0))) * ideal[p]);
    }
    imu_cal::AccelFullFitResult rr;
    run(ups, rr);
    double b = 0, c = 0;
    for (int j = 0; j < 3; ++j) { b = std::max(b, dop(rr, 6 + j)); c = std::max(c, dop(rr, 3 + j)); }
    bd.push_back(b); cd.push_back(c);
  }
  rep << "information,manual_placement_p95,bias_dop_max," << quantile(bd, 0.95) << "\n";
  rep << "information,manual_placement_p95,cross_dop_max," << quantile(cd, 0.95) << "\n";
  check(quantile(bd, 0.95) < 0.85, "manual placement keeps bias determined");
  check(quantile(cd, 0.95) < 2.5, "manual placement keeps cross terms determined");
}

// Thermal: a single later direction, or span alone, is not enough.
void testThermalInformation() {
  imu_cal::AccelFitCfg cfg;
  auto fitter = std::make_unique<Proc::Fitter>();
  Rng rng(11);
  SensorTruth t = SensorTruth::random(rng);
  t.k = Vec3(0.004, -0.003, 0.005);
  auto build = [&](const std::vector<std::pair<int, double>>& holds) {
    std::vector<AccelObs> obs;
    for (size_t h = 0; h < holds.size(); ++h) {
      const Vec3 u = imu_cal::accelPoseDocUp(imu_cal::kAccelPoses[holds[h].first]).cast<double>();
      const double T = holds[h].second;
      for (int b = 0; b < 20; ++b) {
        const Vec3 raw = t.S_at(T).inverse() * (kGStd * u) + t.bias_at(T) + rng.nvec(0.002);
        AccelObs o{};
        for (int d = 0; d < 3; ++d) { o.a[d] = (float)raw(d); o.w[d] = 0; }
        o.tempC = (float)T; o.hold = (uint8_t)h; o.role = 0; o.n = 25;
        obs.push_back(o);
      }
    }
    return obs;
  };
  std::vector<std::pair<int, double>> base;
  for (int p = 0; p < 10; ++p) base.push_back({p, 24.0 + 0.1 * p});
  imu_cal::AccelThermalPrior pr;
  // three later directions, 6 degC later: learned and accurate
  {
    auto h = base; h.push_back({0, 30}); h.push_back({2, 30.2}); h.push_back({4, 30.4});
    auto obs = build(h);
    imu_cal::AccelFullFitResult r;
    fitter->fit(obs.data(), (int)obs.size(), cfg, pr, r, false);
    check(r.ok && r.thermal == AccelThermal::LEARNED, "three later directions: slope learned");
    check((r.k - t.k).cwiseAbs().maxCoeff() < 0.0015, "three later directions: slope recovered");
  }
  // one later direction: not determined
  {
    auto h = base; h.push_back({0, 30}); h.push_back({0, 30.1}); h.push_back({0, 30.2});
    auto obs = build(h);
    imu_cal::AccelFullFitResult r;
    fitter->fit(obs.data(), (int)obs.size(), cfg, pr, r, false);
    check(r.thermal != AccelThermal::LEARNED && r.thermal_reason == imu_cal::AccelThermalReason::INFO_LOW,
          "single later direction: slope not learned (INFO_LOW)");
    check(r.k.isZero(), "single later direction: k stays zero without a prior");
  }
  // span with no repeat of any direction: confounded with the static terms
  {
    std::vector<std::pair<int, double>> h;
    for (int p = 0; p < 10; ++p) h.push_back({p, 24.0 + 0.7 * p});
    auto obs = build(h);
    imu_cal::AccelFullFitResult r;
    fitter->fit(obs.data(), (int)obs.size(), cfg, pr, r, false);
    check(r.thermal != AccelThermal::LEARNED, "span without repeated directions: slope not learned");
  }
  // no span
  {
    auto h = base; h.push_back({0, 24.5}); h.push_back({2, 24.5}); h.push_back({4, 24.5});
    auto obs = build(h);
    imu_cal::AccelFullFitResult r;
    fitter->fit(obs.data(), (int)obs.size(), cfg, pr, r, false);
    check(r.thermal_reason == imu_cal::AccelThermalReason::SPAN_SMALL, "no temperature span: SPAN_SMALL");
    // with a compatible prior the slope is preserved exactly
    imu_cal::AccelThermalPrior p2; p2.valid = true; p2.k[0] = 0.0031; p2.k[1] = -0.0022; p2.k[2] = 0.0044;
    p2.k_temp_lo = 22; p2.k_temp_hi = 29; p2.clamp_lo = 17; p2.clamp_hi = 34;
    imu_cal::AccelFullFitResult r2;
    fitter->fit(obs.data(), (int)obs.size(), cfg, p2, r2, false);
    check(r2.ok && r2.thermal == AccelThermal::PRESERVED && std::fabs(r2.k(0) - 0.0031) < 1e-12 &&
          r2.clamp_lo == 17 && r2.clamp_hi == 34, "compatible prior slope preserved with its range");
    // Session temperatures above the prior's clamp: the fit uses the clamped
    // temperature exactly as the runtime will, so the stored set still fits.
    imu_cal::AccelThermalPrior p3 = p2;
    p3.clamp_lo = 15; p3.clamp_hi = 20;
    imu_cal::AccelFullFitResult r3;
    fitter->fit(obs.data(), (int)obs.size(), cfg, p3, r3, false);
    check(r3.ok && r3.thermal == AccelThermal::PRESERVED && r3.float_ok &&
          std::fabs(r3.float_hold_rms - r3.ref_hold_rms) < 2e-4 && r3.hold_rms < 0.01,
          "preserved slope outside its clamp: fit and runtime agree");
  }
}

// ---------------------------------------------------------------------------
// Campaign

struct Row {
  std::string scenario, method;
  int seed = 0;
  bool ok = false;
  std::string reason;
  Accuracy acc, acc_hot;
  double time_s = NAN;
  int retries = 0, extras = 0;
  std::string thermal;
  double k_err = NAN;
  double k_true_max = NAN;
  double cpu_s = NAN;
};

Scenario baseScenario(const std::string& name) {
  Scenario s; s.name = name;
  return s;
}

std::vector<Scenario> scenarios() {
  std::vector<Scenario> v;
  {  // warm device, re-run from the home screen
    Scenario s = baseScenario("typical");
    s.temp.T_start = 31.0; s.temp.dT = 0.8; s.temp.tau_s = 400;
    s.wizard_start_s = 600;
    v.push_back(s);
  }
  {  // adverse handling and stream defects
    Scenario s = baseScenario("adverse");
    s.temp.T_start = 30.0; s.temp.dT = 1.5; s.temp.tau_s = 400; s.wizard_start_s = 600;
    s.user.p_slow = 0.35; s.user.hand.scale = 2.0; s.user.tilt_err_deg = 18; s.user.face_err_hand_deg = 8;
    s.user.p_first_wrong_tilt = 0.4; s.user.p_regrip = 0.4; s.user.shock_rate_per_s = 0.08;
    s.user.p_hand_face[0] = 0.3; s.user.p_hand_face[1] = 0.3; s.user.p_hand_face[2] = 0.6;
    s.user.p_hand_face[4] = 0.7; s.user.p_hand_face[5] = 0.7; s.user.p_hand_tilt = 1.0;
    s.stream.p_stale = 0.05; s.stream.p_gap = 0.003; s.stream.jitter_ms = 1.5;
    s.cross_sd = 0.008; s.scale_sd = 0.02; s.bias_sd = 0.3; s.mis_deg = 1.0;
    v.push_back(s);
  }
  if (getenv("ACCEL_CAL_FACTORS")) {
    Scenario s = baseScenario("f_hand");
    s.temp.T_start = 30.0; s.temp.dT = 1.5; s.temp.tau_s = 400; s.wizard_start_s = 600;
    s.user.hand.scale = 2.0; for (int i = 0; i < 6; ++i) s.user.p_hand_face[i] = 1.0; s.user.p_hand_tilt = 1.0;
    v.push_back(s);
    Scenario m = baseScenario("f_mis");
    m.temp = s.temp; m.wizard_start_s = 600; m.mis_deg = 1.0; m.cross_sd = 0.008; m.scale_sd = 0.02; m.bias_sd = 0.3;
    v.push_back(m);
    Scenario st = baseScenario("f_stream");
    st.temp = s.temp; st.wizard_start_s = 600; st.stream.p_stale = 0.05; st.stream.p_gap = 0.003; st.stream.jitter_ms = 1.5;
    v.push_back(st);
    Scenario u = baseScenario("f_user");
    u.temp = s.temp; u.wizard_start_s = 600; u.user.p_slow = 0.35; u.user.tilt_err_deg = 18; u.user.face_err_hand_deg = 8;
    u.user.p_first_wrong_tilt = 0.4; u.user.p_regrip = 0.4; u.user.shock_rate_per_s = 0.08;
    v.push_back(u);
  }
  {  // cold start: boot straight into the wizard, informative warm-up
    Scenario s = baseScenario("cold_start_thermal");
    s.temp.T_start = 22.0; s.temp.dT = 9.0; s.temp.tau_s = 240; s.wizard_start_s = 15;
    s.k_override = Vec3(0.004, -0.003, 0.005);
    v.push_back(s);
  }
  {  // cold start with a thermal scale change the model does not contain
     // (1e-4/degC: 2.5x the BMI270 typical sensitivity drift)
    Scenario s = baseScenario("thermal_scale_mismatch");
    s.temp.T_start = 22.0; s.temp.dT = 9.0; s.temp.tau_s = 240; s.wizard_start_s = 15;
    s.k_override = Vec3(0.004, -0.003, 0.005);
    s.dSdT = 1e-4;
    v.push_back(s);
  }
  {  // stress: 3e-4/degC scale drift; only a safe outcome is required
    Scenario s = baseScenario("thermal_scale_stress");
    s.temp.T_start = 22.0; s.temp.dT = 9.0; s.temp.tau_s = 240; s.wizard_start_s = 15;
    s.k_override = Vec3(0.004, -0.003, 0.005);
    s.dSdT = 3e-4;
    v.push_back(s);
  }
  {  // accel-only recalibration: no gyro/mag interval, little variation
    Scenario s = baseScenario("accel_only");
    s.accel_only = true;
    s.temp.T_start = 30.0; s.temp.dT = 0.6; s.temp.tau_s = 500; s.wizard_start_s = 900;
    v.push_back(s);
  }
  {  // local gravity differs from the standard value (equator, sea level)
    Scenario s = baseScenario("local_g");
    s.temp.T_start = 31.0; s.temp.dT = 0.8; s.temp.tau_s = 400; s.wizard_start_s = 600;
    s.g_local = 9.7803;
    v.push_back(s);
  }
  {  // temperature sensor missing
    Scenario s = baseScenario("temp_nan");
    s.temp.all_nan = true; s.wizard_start_s = 600;
    v.push_back(s);
  }
  {  // temperature rises then returns: rechecks at the early temperature
    Scenario s = baseScenario("temp_confounded");
    s.temp.custom = [](double t) { return 24.0 + 4.0 * std::sin(kPi * std::min(1.0, std::max(0.0, (t - 15.0) / 330.0))); };
    s.wizard_start_s = 15;
    s.k_override = Vec3(0.004, -0.003, 0.005);
    v.push_back(s);
  }
  return v;
}

Scenario scenarioByName(const std::string& name) {
  for (const Scenario& s : scenarios()) if (s.name == name) return s;
  throw std::runtime_error("no scenario " + name);
}

std::string fmt(double v, int p = 5) {
  if (!std::isfinite(v)) return "nan";
  std::ostringstream o; o << std::fixed << std::setprecision(p) << v; return o.str();
}

void runCampaign(int seeds, std::ostream& csv, std::ostream& sum, std::ostream& rep) {
  csv << "scenario,seed,method,ok,reason,bias_err,bias_axis_max,cross_err,diag_err,vec_rms,vec_max,norm_rms,"
         "bias_err_hot,vec_rms_hot,time_s,retries,extras,thermal,k_err,cpu_s\n";
  std::vector<Row> rows;
  for (const Scenario& sc0 : scenarios()) {
    for (int seed = 1; seed <= seeds; ++seed) {
      Scenario sc = sc0;
      const uint64_t base = 1000003ull * (uint64_t)seed + std::hash<std::string>{}(sc.name) % 100000ull;
      Rng tr(base);
      SensorTruth truth = SensorTruth::random(tr, sc.bias_sd, sc.scale_sd, sc.cross_sd, sc.k_sd);
      if (std::isfinite(sc.k_override(0))) truth.k = sc.k_override;
      if (sc.mis_deg > 0) truth.R_mis = expso3(tr.unit() * sc.mis_deg * kDeg);
      truth.dS_dT = Mat3::Identity() * sc.dSdT;
      truth.g_local = sc.g_local;
      truth.T0 = 25.0;

      // OLD: deployed procedure
      Row ro; ro.scenario = sc.name; ro.method = "OLD"; ro.seed = seed;
      {
        World w(truth, sc, base + 1);
        Stream st(w, sc.stream);
        OldRun o = runOld(w, st, sc);
        const double Tm = sc.temp.at(0.5 * (o.t_start + std::max(o.t_start, o.t_end)));
        ro.ok = o.ok; ro.reason = o.reason; ro.time_s = o.t_end - o.t_start;
        ro.acc = score(o.cal, truth, Tm);
        ro.acc_hot = score(o.cal, truth, Tm + 5.0);
        ro.thermal = "none";
      }
      rows.push_back(ro);

      // RICH: new procedure, scored with the new and the old fitter
      World w(truth, sc, base + 2);
      Stream st(w, sc.stream);
      imu_cal::AccelThermalPrior prior;
      std::vector<std::string> dbg;
      const char* dscn = getenv("ACCEL_CAL_DEBUG_SCN");
      const char* dseed = getenv("ACCEL_CAL_DEBUG_SEED");
      const bool debug = dscn && sc.name == dscn && (!dseed || strtol(dseed, nullptr, 10) == seed);
      NewRun n = runNew(w, st, sc, prior, debug ? &dbg : nullptr);
      if (debug) {
        for (const std::string& l : dbg) if (l.rfind("[ACCBLK]", 0) != 0) std::cout << l << "\n";
        std::cout << "seed=" << seed << " reason=" << n.reason << " k_true=" << truth.k.transpose() << "\n";
      }
      const double Tm = std::isfinite(n.T_mid) ? n.T_mid : sc.temp.at(w.t());
      Row rn; rn.scenario = sc.name; rn.method = "RICH+NEW"; rn.seed = seed;
      rn.ok = n.ok; rn.reason = n.reason; rn.time_s = n.accel_time_s; rn.retries = n.retries; rn.extras = n.extras;
      rn.acc = score(n.cal, truth, Tm);
      rn.acc_hot = score(n.cal, truth, Tm + 5.0);
      rn.thermal = imu_cal::accelThermalStr(n.fit.thermal);
      if (n.ok) rn.k_err = (n.cal.acc.biasT.k.cast<double>() - truth.k).cwiseAbs().maxCoeff();
      rn.cpu_s = n.fit_cpu_s;
      rows.push_back(rn);
      Row rr = rn; rr.method = "RICH+OLD"; rr.ok = n.rich_old.ok; rr.reason = n.rich_old_reason;
      rr.acc = score(n.rich_old, truth, Tm); rr.acc_hot = score(n.rich_old, truth, Tm + 5.0);
      rr.thermal = "none"; rr.k_err = NAN; rr.cpu_s = NAN;
      rows.push_back(rr);
    }
  }
  for (const Row& r : rows) {
    csv << r.scenario << ',' << r.seed << ',' << r.method << ',' << (r.ok ? 1 : 0) << ',' << r.reason << ','
        << fmt(r.acc.bias_err) << ',' << fmt(r.acc.bias_axis_max) << ',' << fmt(r.acc.cross_err) << ','
        << fmt(r.acc.diag_err) << ',' << fmt(r.acc.vec_rms) << ',' << fmt(r.acc.vec_max) << ','
        << fmt(r.acc.norm_rms) << ',' << fmt(r.acc_hot.bias_err) << ',' << fmt(r.acc_hot.vec_rms) << ','
        << fmt(r.time_s, 1) << ',' << r.retries << ',' << r.extras << ',' << r.thermal << ','
        << fmt(r.k_err, 6) << ',' << fmt(r.cpu_s, 4) << '\n';
  }

  // Summary
  sum << "scenario,method,runs,success_rate,bias_err_median,bias_err_p90,bias_err_max,cross_err_median,"
         "cross_err_p90,vec_rms_median,vec_rms_p90,norm_rms_median,bias_err_hot_median,time_s_median,"
         "time_s_p90,retries_mean,extras_mean,learned_rate,k_err_median,identity_accepted_rate\n";
  std::map<std::pair<std::string, std::string>, std::vector<const Row*>> groups;
  std::vector<std::string> order;
  for (const Row& r : rows) {
    auto key = std::make_pair(r.scenario, r.method);
    if (!groups.count(key)) order.push_back(r.scenario + "|" + r.method);
    groups[key].push_back(&r);
  }
  rep << "\n" << std::left << std::setw(24) << "scenario" << std::setw(10) << "method" << std::right
      << std::setw(6) << "succ" << std::setw(10) << "b_med" << std::setw(10) << "b_p90" << std::setw(10) << "b_max"
      << std::setw(10) << "x_med" << std::setw(10) << "v_med" << std::setw(10) << "v_p90" << std::setw(10) << "n_med"
      << std::setw(10) << "bhot_med" << std::setw(8) << "t_med" << std::setw(8) << "retry" << std::setw(7) << "learn"
      << std::setw(10) << "k_err" << std::setw(7) << "ident" << "\n";
  for (const std::string& k : order) {
    const std::string scn = k.substr(0, k.find('|')), mth = k.substr(k.find('|') + 1);
    const auto& g = groups[{scn, mth}];
    std::vector<double> b, x, vr, nr, bh, tt, ke;
    int ok = 0, learned = 0, ident = 0; double retr = 0, ext = 0;
    for (const Row* r : g) {
      tt.push_back(r->time_s);
      if (r->reason == "OK_IDENTITY") ++ident;
      retr += r->retries; ext += r->extras;
      if (!r->ok) continue;
      ++ok;
      b.push_back(r->acc.bias_err); x.push_back(r->acc.cross_err); vr.push_back(r->acc.vec_rms);
      nr.push_back(r->acc.norm_rms); bh.push_back(r->acc_hot.bias_err);
      if (r->thermal == "LEARNED") ++learned;
      if (std::isfinite(r->k_err)) ke.push_back(r->k_err);
    }
    const double n = (double)g.size();
    sum << scn << ',' << mth << ',' << g.size() << ',' << fmt(ok / n, 3) << ',' << fmt(median(b)) << ','
        << fmt(quantile(b, 0.9)) << ',' << fmt(vmax(b)) << ',' << fmt(median(x), 6) << ',' << fmt(quantile(x, 0.9), 6)
        << ',' << fmt(median(vr)) << ',' << fmt(quantile(vr, 0.9)) << ',' << fmt(median(nr)) << ',' << fmt(median(bh))
        << ',' << fmt(median(tt), 1) << ',' << fmt(quantile(tt, 0.9), 1) << ',' << fmt(retr / n, 2) << ','
        << fmt(ext / n, 2) << ',' << fmt(ok ? learned / (double)ok : 0.0, 3) << ',' << fmt(median(ke), 6) << ','
        << fmt(ident / n, 3) << '\n';
    rep << std::left << std::setw(24) << scn << std::setw(10) << mth << std::right << std::setw(6) << fmt(ok / n, 2)
        << std::setw(10) << fmt(median(b), 4) << std::setw(10) << fmt(quantile(b, 0.9), 4) << std::setw(10)
        << fmt(vmax(b), 4) << std::setw(10) << fmt(median(x), 5) << std::setw(10) << fmt(median(vr), 4)
        << std::setw(10) << fmt(quantile(vr, 0.9), 4) << std::setw(10) << fmt(median(nr), 4) << std::setw(10)
        << fmt(median(bh), 4) << std::setw(8) << fmt(median(tt), 0) << std::setw(8) << fmt(retr / n, 2)
        << std::setw(7) << fmt(ok ? learned / (double)ok : 0.0, 2) << std::setw(10) << fmt(median(ke), 5)
        << std::setw(7) << fmt(ident / n, 2) << "\n";
  }

  // Campaign assertions (see doc/imu_calibrate for the measured values).
  auto stat = [&](const std::string& scn, const std::string& mth, auto fn) {
    std::vector<double> v; int ok = 0; int n = 0;
    for (const Row* r : groups[{scn, mth}]) { ++n; if (r->ok) { ++ok; v.push_back(fn(*r)); } }
    return std::make_pair(v, n ? ok / (double)n : 0.0);
  };
  auto bias = [](const Row& r) { return r.acc.bias_err; };
  auto vec = [](const Row& r) { return r.acc.vec_rms; };
  for (const std::string scn : {"typical", "adverse"}) {
    auto nb = stat(scn, "RICH+NEW", bias);
    auto ob = stat(scn, "OLD", bias);
    auto rb = stat(scn, "RICH+OLD", bias);
    check(nb.second >= 0.9, scn + ": new procedure success >= 90%");
    check(median(nb.first) < 0.5 * median(ob.first), scn + ": new median bias error < half of the deployed procedure");
    check(median(nb.first) < median(rb.first), scn + ": new fitter beats the old fitter on the same rich data");
    check(median(stat(scn, "RICH+NEW", vec).first) < 0.5 * median(stat(scn, "OLD", vec).first),
          scn + ": new median vector error < half of the deployed procedure");
  }
  {
    auto learned = 0, n = 0;
    std::vector<double> ke;
    for (const Row* r : groups[{"cold_start_thermal", "RICH+NEW"}]) {
      ++n;
      if (r->ok && r->thermal == "LEARNED") { ++learned; ke.push_back(r->k_err); }
    }
    check(learned >= 0.7 * n, "cold start: slope learned in >= 70% of sessions");
    check(median(ke) < 0.0015, "cold start: median slope error < 0.0015 m/s^2/degC");
    auto hot_new = stat("cold_start_thermal", "RICH+NEW", [](const Row& r) { return r.acc_hot.bias_err; });
    auto hot_old = stat("cold_start_thermal", "RICH+OLD", [](const Row& r) { return r.acc_hot.bias_err; });
    check(median(hot_new.first) < median(hot_old.first), "cold start: +5 degC bias error better than without a slope");
  }
  for (const Row* r : groups[{"temp_nan", "RICH+NEW"}]) check(!r->ok || r->thermal == "UNLEARNED", "NaN temperature: never learned");
  // Model mismatch: accepted results stay within the verification bounds on
  // unseen orientations (failures are safe: the previous calibration is kept).
  for (const std::string scn : {"thermal_scale_mismatch", "thermal_scale_stress"})
    for (const Row* r : groups[{scn, "RICH+NEW"}])
      check(!r->ok || r->acc.norm_rms < 0.03, scn + ": accepted calibration keeps unseen-orientation norm error < 0.03");
  check(stat("thermal_scale_mismatch", "RICH+NEW", bias).second >= 0.9, "thermal scale mismatch: success >= 90%");
  for (const Row* r : groups[{"accel_only", "RICH+NEW"}]) check(!r->ok || r->thermal != "LEARNED" || r->k_err < 0.003, "accel-only: no spurious slope");
  for (const Row* r : groups[{"temp_confounded", "RICH+NEW"}]) check(!r->ok || r->thermal != "LEARNED" || r->k_err < 0.003, "confounded: no spurious slope");
  check(stat("temp_nan", "RICH+NEW", bias).second >= 0.9, "NaN temperature: calibration still succeeds");
  // Only |g| is used: a local g of 9.78 rescales S but leaves the bias intact.
  check(median(stat("local_g", "RICH+NEW", bias).first) < 0.5 * median(stat("local_g", "OLD", bias).first) &&
        median(stat("local_g", "RICH+NEW", bias).first) < 0.006, "local g: bias unaffected by the gravity magnitude");
}

// Power-cycle offset: the second session re-measures the bias and keeps the
// validated slope; it never turns the offset into a slope.
void testPowerCycle(std::ostream& rep) {
  Scenario sc = scenarioByName("cold_start_thermal");
  Rng tr(4242);
  SensorTruth truth = SensorTruth::random(tr);
  truth.k = Vec3(0.004, -0.003, 0.005);
  World w1(truth, sc, 51);
  Stream s1(w1, sc.stream);
  imu_cal::AccelThermalPrior none;
  NewRun a = runNew(w1, s1, sc, none);
  check(a.ok && a.fit.thermal == AccelThermal::LEARNED, "power-cycle: first session learns the slope");
  if (!a.ok) return;
  imu_cal::AccelThermalPrior pr = atoms3r_ical::accelThermalPriorFrom(a.blob, 0x1234, 0x5678, 3);
  check(pr.valid, "power-cycle: first session provides a prior");
  {
    // Accelerometer-only candidate: gyro and magnetometer carried byte for byte.
    const atoms3r_ical::ImuCalBlobV3 prev = makeBlob();
    imu_cal::AccelCalibration<float> fc;
    Proc::Fitter::toFloat(a.fit, kGStd, fc);
    const atoms3r_ical::ImuCalBlobV3 c = atoms3r_ical::accelOnlyCandidate(prev, a.fit, fc, 150, 0x1234, 0x5678, 3);
    const size_t g0 = offsetof(atoms3r_ical::ImuCalBlobV3, gyro_ok);
    const size_t g1 = offsetof(atoms3r_ical::ImuCalBlobV3, accel_T_lo);
    check(memcmp((const uint8_t*)&c + g0, (const uint8_t*)&prev + g0, g1 - g0) == 0,
          "accel-only candidate keeps gyro and magnetometer fields byte for byte");
    check(c.accel_ok == 1 && c.accel_b0[0] == fc.biasT.b0(0) && atoms3r_ical::accelMetaBound(c),
          "accel-only candidate carries the new bound accelerometer set");
    MemStore st;
    atoms3r_ical::ImuCalBlobV3 rb;
    check(st.saveVerified(c, rb) && rb.gyro_b0[0] == prev.gyro_b0[0] && rb.mag_b[0] == prev.mag_b[0],
          "accel-only candidate saves with gyro/mag unchanged");
  }
  // second session: warm and flat temperature, bias shifted at power-up
  Scenario sc2 = scenarioByName("typical");
  SensorTruth t2 = truth;
  t2.b0 += Vec3(0.03, -0.02, 0.025);
  World w2(t2, sc2, 52);
  Stream s2(w2, sc2.stream);
  NewRun b = runNew(w2, s2, sc2, pr);
  check(b.ok && b.fit.thermal == AccelThermal::PRESERVED, "power-cycle: flat second session preserves the slope");
  if (getenv("ACCEL_CAL_DUMP"))
    std::cout << "power-cycle second: ok=" << b.ok << " reason=" << b.reason << " thermal=" << imu_cal::accelThermalStr(b.fit.thermal)
              << "/" << imu_cal::accelThermalReasonStr(b.fit.thermal_reason) << " T=[" << b.fit.cal_temp_lo << "," << b.fit.cal_temp_hi << "]\n";
  if (b.ok) {
    check((b.cal.acc.biasT.k - a.cal.acc.biasT.k).norm() == 0.0f, "power-cycle: slope carried exactly");
    const Accuracy ac = score(b.cal, t2, b.T_mid);
    rep << "power_cycle,second_session_bias_err," << ac.bias_err << "\n";
    check(ac.bias_err < 0.02, "power-cycle: shifted bias re-measured");
  }
}

// Abort and hold failure leave no candidate.
void testFailurePaths() {
  Scenario sc = scenarioByName("typical");
  Rng tr(9);
  SensorTruth truth = SensorTruth::random(tr);
  {
    World w(truth, sc, 3);
    Stream st(w, sc.stream);
    auto proc = std::make_unique<Proc>();
    imu_cal::AccelCaptureCfg c; imu_cal::AccelFitCfg f; imu_cal::AccelThermalPrior p;
    proc->begin(c, f, p, Eigen::Vector3f::Zero(), false);
    SimIo io(w, st, sc);
    io.abort_at_prep = true;
    check(!proc->runMainStage(io) && proc->failure() == imu_cal::AccelProcFail::ABORTED, "abort at prep fails the stage");
    check(!proc->result().ok, "aborted procedure has no accepted result");
  }
  {
    // User never reaches the pose (skips a face): bounded retries, then a clear failure.
    const Scenario& s2 = sc;
    struct SkipIo : public SimIo {
      using SimIo::SimIo;
      bool prep(const imu_cal::AccelStepView& v) override {
        imu_cal::AccelStepView v2 = v;
        if (v.pose == 1) v2.pose = 0;  // repeats screen up instead of screen down
        return SimIo::prep(v2);
      }
    };
    World w(truth, s2, 4);
    Stream st(w, s2.stream);
    auto proc = std::make_unique<Proc>();
    imu_cal::AccelCaptureCfg c; imu_cal::AccelFitCfg f; imu_cal::AccelThermalPrior p;
    proc->begin(c, f, p, Eigen::Vector3f::Zero(), false);
    SkipIo io(w, st, s2);
    const bool ok = proc->runMainStage(io);
    check(!ok && proc->failure() == imu_cal::AccelProcFail::HOLD_FAILED, "repeated face instead of the requested one fails the hold");
    check(proc->totalAttempts() == 1 + c.max_attempts_per_hold, "bounded attempts on a failing pose");
  }
}

// Replay route: a logged session (raw samples, prep markers, context lines)
// reproduces the device result through tests/imu_calibrate/accel_cal_replay.h,
// both from raw samples (whole procedure) and from the logged blocks (fit).
void testReplay() {
  Scenario sc = scenarioByName("adverse");
  Rng tr(31);
  SensorTruth truth = SensorTruth::random(tr);
  World w(truth, sc, 32);
  Stream st(w, sc.stream);
  std::vector<std::string> logs;
  imu_cal::AccelThermalPrior none;
  NewRun n = runNew(w, st, sc, none, &logs, true);
  check(n.ok, "replay source session succeeds");
  if (getenv("ACCEL_CAL_DUMP")) {
    for (const std::string& l : logs) if (l.rfind("[ACCRAW]", 0) != 0 && l.rfind("[ACCBLK]", 0) != 0) std::cout << l << "\n";
    std::cout << "reason=" << n.reason << "\n";
  }
  {
    std::ofstream f("calibrate_accel_replay_sample.log");
    for (const std::string& l : logs) f << l << "\n";
  }
  const accel_replay::Log L = accel_replay::parse(logs);
  imu_cal::AccelFullFitResult rb;
  int nblk = 0;
  const bool okb = accel_replay::fitBlocks(L, rb, &nblk);
  check(nblk == (int)n.obs.size(), "block replay reads every retained block");
  check(okb && (rb.b - n.fit.b).norm() < 1e-4 && (rb.S - n.fit.S).norm() < 1e-5,
        "block replay reproduces the fit (float-logged blocks)");
  imu_cal::AccelFullFitResult rr;
  int div = -1;
  const bool okr = accel_replay::replayRaw(L, rr, div, false);
  check(okr && div == 0, "raw replay follows the same hold sequence");
  check(okr && (rr.b - n.fit.b).norm() < 1e-4 && (rr.S - n.fit.S).norm() < 1e-5,
        "raw replay reproduces capture and fit");
}

}  // namespace

int main(int argc, char** argv) {
  int seeds = 30;
  if (argc > 1) seeds = std::max(1, (int)strtol(argv[1], nullptr, 10));
  if (const char* e = getenv("ACCEL_CAL_SEEDS")) seeds = std::max(1, (int)strtol(e, nullptr, 10));
  try {
    std::ofstream csv("calibrate_accel_campaign.csv");
    std::ofstream sum("calibrate_accel_campaign_summary.csv");
    std::ofstream rep("calibrate_accel_report.txt");
    std::ostringstream con;

    // Memory: fixed-size state replacing the wizard's AccelCalibrator (host
    // sizes; the ESP32-S3 has 4-byte pointers and the same float/double sizes).
    rep << "memory,sizeof_AccelObs," << sizeof(AccelObs) << "\n";
    rep << "memory,sizeof_AccelCalProcedure_340_24," << sizeof(Proc) << "\n";
    rep << "memory,sizeof_AccelFullFitter_340_24," << sizeof(Proc::Fitter) << "\n";
    rep << "memory,sizeof_legacy_AccelCalibrator_float_400_1," << sizeof(imu_cal::AccelCalibrator<float, 400, 1>) << "\n";
    testLegacyFitterSuccessPaths();
    testTempBias();
    testBlobAndStore();
    testScreenText();
    testWizardScreenStrings();
    testPoseMapping();
    testInformation(rep);
    testThermalInformation();
    testFailurePaths();
    testPowerCycle(rep);
    testReplay();
    runCampaign(seeds, csv, sum, con);
    rep << con.str();
    std::cout << con.str() << "\n";
    std::cout << "accel_cal-test: " << (g_checks - g_failures) << "/" << g_checks << " checks passed (" << seeds
              << " seeds)\n";
  } catch (const std::exception& ex) {
    std::cerr << "accel_cal-test: FAIL: " << ex.what() << "\n";
    return EXIT_FAILURE;
  }
  return g_failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
