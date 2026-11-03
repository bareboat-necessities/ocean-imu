// Gyro scale and cross-axis calibration from the accelerometer pose session
// (src/imu_calibrate/GyroAxisFit.h).
//
//  1. GyroChain attitude against a fine reference integration, and its
//     sensitivity columns against finite differences.
//  2. Re-anchoring a hold's sums against direct computation.
//  3. The optional stored record: binding to the saved gyro set, integrity
//     and plausibility checks, removal, and the runtime loader.
//  4. Simulated wizard sessions through the AccelCalProcedure: the user model
//     of accel_cal_sim.h, a true gyro matrix, read-time jitter, sample gaps,
//     drift of the rate bias since the stationary capture. The fitted matrix
//     must reduce the rate error, and must never make nominal axes worse.

#include "accel_cal_sim.h"
#include "AtomS3R/AtomS3R_ImuCalBlob.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace accel_sim;
using Proc = imu_cal::AccelCalProcedure<340, 24>;

static int g_fail = 0;
static void check(bool ok, const std::string& what) {
  std::printf("%s %s\n", ok ? "PASS" : "FAIL", what.c_str());
  if (!ok) ++g_fail;
}

// ---------------------------------------------------------------------------
// 1. Integrator and sensitivity

static Vec3 rateAt(double t) {
  return Vec3(1.8 * std::sin(2.1 * t) + 0.4, -1.2 * std::cos(1.3 * t + 0.4), 2.5 * std::sin(0.7 * t + 1.0));
}

// Feeds rates (1+E)*w - offset at 200 Hz with jittered stamps; returns R.
static Mat3 integrateChain(const Eigen::Matrix3f& E, const Eigen::Vector3f& offset, double T,
                           imu_cal::GyroChainState* last = nullptr) {
  imu_cal::GyroChain chain;
  Rng r(7);
  const imu_cal::GyroChainState* s = nullptr;
  for (int k = 0; k * 0.005 <= T; ++k) {
    const double t = k * 0.005;
    const Eigen::Vector3f w = rateAt(t).cast<float>();
    const uint32_t us = uint32_t(std::llround((t + r.uni(0, 0.0015)) * 1e6));
    s = chain.feed(us, (Eigen::Matrix3f::Identity() + E) * w - offset);
  }
  if (last && s) *last = *s;
  return s ? Mat3(s->R.cast<double>()) : Mat3(Mat3::Identity());
}

static void testIntegrator() {
  const double T = 3.0;
  // Fine reference: body(t) -> body(0).
  Mat3 R = Mat3::Identity();
  const double h = 1e-5;
  const long steps = std::lround(T / h);
  for (long k = 0; k < steps; ++k) R = R * expso3(rateAt((double(k) + 0.5) * h) * h);
  imu_cal::GyroChainState s0;
  const Mat3 Rc = integrateChain(Eigen::Matrix3f::Zero(), Eigen::Vector3f::Zero(), T, &s0);
  const double err = logso3(Rc.transpose() * R).norm();
  char b[128];
  std::snprintf(b, sizeof(b), "chain attitude matches fine integration (%.2e rad over 3 s at 200 Hz)", err);
  check(err < 2e-3, b);

  double worst = 0;
  const double eps = 1e-3;
  for (int k = 0; k < imu_cal::kGyroAxisParams; ++k) {
    Eigen::Matrix3f E = Eigen::Matrix3f::Zero();
    Eigen::Vector3f off = Eigen::Vector3f::Zero();
    if (k < 9) E(k / 3, k % 3) = float(eps); else off(k - 9) = float(eps);
    const Mat3 Rp = integrateChain(E, off, T);
    const Vec3 fd = logso3(Rp * Rc.transpose()) / eps;
    const Vec3 an = s0.P.col(k).cast<double>();
    worst = std::max(worst, (fd - an).norm() / std::max(1.0, an.norm()));
  }
  std::snprintf(b, sizeof(b), "sensitivity columns match finite differences (worst rel %.2e)", worst);
  check(worst < 2e-2, b);
}

// ---------------------------------------------------------------------------
// 2. Re-anchoring

static void testReanchor() {
  imu_cal::GyroChain chain;
  imu_cal::GyroPoseSums hold;
  std::vector<Mat3> Rk;
  const Eigen::Vector3f a(0.3f, -9.7f, 1.2f);
  imu_cal::GyroChainState last;
  for (int k = 0; k <= 600; ++k) {
    const double t = k * 0.005;
    const auto* s = chain.feed(uint32_t(std::llround(t * 1e6)), rateAt(t).cast<float>());
    if (k >= 500) { hold.add(s, a); Rk.push_back(s->R.cast<double>()); }
    last = *s;
  }
  const imu_cal::GyroPoseSums ex = chain.reanchor(hold);
  Vec3 Ra = Vec3::Zero();
  const Mat3 Rat = last.R.cast<double>().transpose();
  for (const Mat3& R : Rk) Ra += Rat * R * a.cast<double>();
  const double e = (ex.Ra.cast<double>() - Ra).norm() / Ra.norm();
  check(ex.ok() && e < 1e-5, "re-anchored sums equal the hold in the new anchor frame");
  check(ex.anchor != hold.anchor && ex.anchor == chain.anchor(), "re-anchoring starts a new anchor");
  const auto* s = chain.feed(uint32_t(std::llround(3.005 * 1e6)), rateAt(3.005).cast<float>());
  check(s && (s->R - Eigen::Matrix3f::Identity()).norm() < 0.02f && s->anchor == ex.anchor,
        "chain continues from the new anchor");
  imu_cal::GyroPoseSums stale = hold;
  stale.anchor += 100;
  check(!chain.reanchor(stale).ok(), "sums from another anchor are not re-anchored");
  const uint32_t before = chain.anchor();
  s = chain.feed(uint32_t(std::llround(3.2 * 1e6)), rateAt(3.2).cast<float>());
  check(s && chain.anchor() == before + 1 && s->R.isIdentity(0.0f), "a gap over 50 ms restarts the chain");
}

// ---------------------------------------------------------------------------
// 3. Stored record

struct MemKv {
  std::map<std::string, std::vector<uint8_t>> bytes;
  size_t getBytesLength(const char* k) { auto it = bytes.find(k); return it == bytes.end() ? 0 : it->second.size(); }
  size_t getBytes(const char* k, void* b, size_t n) {
    auto it = bytes.find(k);
    if (it == bytes.end() || it->second.size() != n) return 0;
    std::memcpy(b, it->second.data(), n);
    return n;
  }
  size_t putBytes(const char* k, const void* b, size_t n) {
    bytes[k].assign((const uint8_t*)b, (const uint8_t*)b + n);
    return n;
  }
  bool remove(const char* k) { bytes.erase(k); return true; }
};
using Store = atoms3r_ical::ImuCalStoreT<MemKv>;
static Store* g_store = nullptr;
static bool loadFromMem(const atoms3r_ical::ImuCalBlobV5& b, Eigen::Matrix3f& S) {
  return g_store && g_store->loadGyroAxes(b, S);
}

static atoms3r_ical::ImuCalBlobV5 blobWithGyro(float b0x) {
  atoms3r_ical::ImuCalBlobV5 b{};
  imu_cal::GyroCalibration<float> g;
  g.ok = true;
  g.biasT.ok = true;
  g.biasT.T0 = 25;
  g.biasT.T_lo = g.biasT.T_hi = 25;
  g.biasT.b0 = Eigen::Vector3f(b0x, -0.002f, 0.001f);
  atoms3r_ical::fillGyroFromFit(b, g);
  b.sensor_id_lo = 0x1234; b.sensor_id_hi = 0x56; b.imu_type = 3;
  return b;
}

static void testRecord() {
  using namespace atoms3r_ical;
  Store store;
  g_store = &store;
  const ImuCalBlobV5 blob = blobWithGyro(0.003f);
  check(gyroSetValid(blob), "test blob has a valid gyro set");
  imu_cal::GyroAxisResult r;
  r.ok = true;
  r.S << 1.012f, 0.004f, -0.003f, -0.006f, 0.991f, 0.002f, 0.005f, -0.001f, 1.006f;
  r.pairs = 11; r.cv_nominal = 0.02f; r.cv_fit = 0.002f;
  const GyroAxisRecord rec = makeGyroAxisRecord(blob, r);
  check(gyroAxisValid(rec) && store.saveGyroAxes(&rec), "a qualified record is written and read back");

  RuntimeCals rc;
  gyroAxisLoader() = nullptr;
  rc.rebuildFromBlob(blob);
  check(!rc.gyro_axes && rc.gyr.S.isIdentity(0.0f), "without a loader the gyro keeps nominal axes");
  gyroAxisLoader() = &loadFromMem;
  rc.rebuildFromBlob(blob);
  check(rc.gyro_axes && rc.gyr.S.isApprox(r.S), "the runtime applies a record bound to its gyro set");
  const Eigen::Vector3f w(0.3f, -0.2f, 0.1f);
  check((rc.applyGyro(w, 25) - r.S * (w - rc.gyr.biasT.b0)).norm() < 1e-6f, "w_cal = S (w_raw - bias)");

  rc.rebuildFromBlob(blobWithGyro(0.004f));
  check(!rc.gyro_axes && rc.gyr.S.isIdentity(0.0f), "a record of another gyro set is not applied");
  ImuCalBlobV5 other = blob;
  other.sensor_id_lo ^= 1;
  rc.rebuildFromBlob(other);
  check(!rc.gyro_axes, "a record of another sensor is not applied");

  GyroAxisRecord bad = rec;
  bad.S[0] += 1e-4f;
  store.kv.putBytes(Store::kKeyGyroAxis, &bad, sizeof(bad));
  rc.rebuildFromBlob(blob);
  check(!rc.gyro_axes, "a record failing its CRC is not applied");
  bad = rec;
  bad.S[4] = 1.2f;
  bad.crc = gyroAxisCrc(bad);
  check(!gyroAxisValid(bad) && !store.saveGyroAxes(&bad), "an implausible matrix is refused");

  check(store.saveGyroAxes(&rec) && store.saveGyroAxes(nullptr) && store.kv.getBytesLength(Store::kKeyGyroAxis) == 0,
        "a nominal result removes a stale record");
  store.saveGyroAxes(&rec);
  store.kv.putBytes(Store::kKeyV5, &blob, sizeof(blob));
  store.erase();
  check(store.kv.bytes.empty(), "erase removes the blob and the gyro record");
  gyroAxisLoader() = nullptr;
  g_store = nullptr;
}

// ---------------------------------------------------------------------------
// 4. Simulated sessions

// The wizard's accelerometer IO with the gyro chain fed from every sample,
// including the preparation and acknowledgement screens.
class GyroSimIo : public imu_cal::AccelCalIo {
public:
  GyroSimIo(World& w, Stream& s, const Scenario& sc, imu_cal::GyroChain* chain, const Eigen::Vector3f& bias)
      : w_(w), s_(s), sc_(sc), chain_(chain), bias_(bias) {}

  void pump(double dur) {
    const double t_end = w_.t() + dur;
    imu_cal::AccelRawSample s;
    while (w_.t() < t_end) {
      if (s_.next(s) && chain_) chain_->feed(s.t_us, s.w - bias_);
    }
  }

  bool prep(const imu_cal::AccelStepView& v) override {
    Rng& r = w_.rng();
    pump(r.uni(sc_.user.react_lo, sc_.user.react_hi));
    const imu_cal::AccelPoseDef& p = imu_cal::kAccelPoses[v.pose];
    const bool tilt = (p.region == imu_cal::AccelRegion::FRAME);
    hand_ = tilt ? r.chance(sc_.user.p_hand_tilt) : r.chance(sc_.user.p_hand_face[v.pose]);
    err_ = tilt ? sc_.user.tilt_err_deg : (hand_ ? sc_.user.face_err_hand_deg : sc_.user.face_err_table_deg);
    target_ = docUpBody(p.top, p.right, p.out);
    const double tm = r.chance(sc_.user.p_slow) ? r.uni(5.0, 9.0) : r.uni(sc_.user.move_lo, sc_.user.move_hi);
    w_.placeAt(target_, hand_, err_, w_.t(), tm);
    check_since_ = -1;
    return true;
  }
  bool sample(imu_cal::AccelRawSample& s) override {
    const bool ok = s_.next(s);
    if (ok && chain_) s.gyro = chain_->feed(s.t_us, s.w - bias_);
    return ok;
  }
  uint32_t nowMs() override { return (uint32_t)(uint64_t)std::llround(w_.t() * 1000.0); }
  void capture(const imu_cal::AccelStepView&, const imu_cal::AccelHoldView& h) override {
    if (h.hint == imu_cal::AccelHoldHint::CHECK_POSE) {
      if (check_since_ < 0) check_since_ = w_.t();
      else if (w_.t() - check_since_ > 2.5) {
        w_.placeAt(target_, hand_, err_, w_.t(), w_.rng().uni(1.0, 2.5));
        check_since_ = -1;
      }
    } else {
      check_since_ = -1;
    }
  }
  void holdOk(const imu_cal::AccelStepView&) override { pump(0.98); }
  bool holdRetry(const imu_cal::AccelStepView&, const char*) override { w_.advance(2.0); return true; }
  bool runFit(imu_cal::AccelFitJob& job, const char*) override { job.run(); w_.advance(0.5); return true; }
  void log(const char*) override {}
  void idle() override {}

private:
  World& w_;
  Stream& s_;
  const Scenario& sc_;
  imu_cal::GyroChain* chain_;
  Eigen::Vector3f bias_;
  bool hand_ = false;
  double err_ = 0;
  Vec3 target_ = Vec3::Zero();
  double check_since_ = -1;
};

struct SessionOut {
  bool accel_ok = false;
  imu_cal::GyroAxisResult r;
  double before = 0, after = 0;   // worst-case relative rate error, nominal vs applied
};

static double rateError(const Mat3& S_applied, const Mat3& S_true) {
  // Applied rate = S_applied * S_true^-1 * true rate.
  Eigen::JacobiSVD<Mat3> svd(S_applied * S_true.inverse() - Mat3::Identity());
  return svd.singularValues()(0);
}

struct Case {
  const char* name;
  double scale_sd, cross_sd;   // gyro matrix truth
  double bias_err;             // rad/s, stationary-capture bias error used for integration
  double p_gap;
  bool use_chain;
  double hand_scale = 1.0;     // adverse users: larger wobble, tremor and drift
};

static SessionOut runSession(const Case& c, uint64_t seed) {
  Scenario sc;
  sc.name = c.name;
  sc.stream.rate_hz = 200.0;
  sc.stream.jitter_ms = 0.2;
  sc.stream.stamp_latency_ms = 2.0;
  sc.stream.p_gap = c.p_gap;
  sc.temp.dT = 1.0;
  sc.user.hand.scale = c.hand_scale;
  Rng tr(seed * 31 + 5);
  SensorTruth truth = SensorTruth::random(tr, sc.bias_sd, sc.scale_sd, sc.cross_sd, sc.k_sd);
  truth.k.setZero();
  Mat3 G = Mat3::Identity();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) G(i, j) += tr.normal(i == j ? c.scale_sd : c.cross_sd);
  truth.gyro_S = G;
  truth.gyro_noise = 0.0012;
  World w(truth, sc, seed);
  w.smooth_holds = true;
  Stream st(w, sc.stream);

  const Eigen::Vector3f bias = (truth.gyro_bias + tr.unit() * c.bias_err).cast<float>();
  imu_cal::GyroChain chain;
  auto proc = std::make_unique<Proc>();
  imu_cal::AccelCaptureCfg ccfg;
  imu_cal::AccelFitCfg fcfg;
  ccfg.g = float(kGSite);
  fcfg.g = float(kGSite);
  proc->begin(ccfg, fcfg, imu_cal::AccelThermalPrior{}, bias, true);
  if (c.use_chain) proc->setGyroChain(&chain);
  GyroSimIo io(w, st, sc, c.use_chain ? &chain : nullptr, bias);

  SessionOut out;
  bool ok = proc->runMainStage(io);
  if (ok) {
    w.advance(w.rng().uni(60.0, 100.0));   // magnetometer stage: not integrated
    ok = proc->runRecheckStage(io) && proc->runFinalFit(io);
  }
  out.accel_ok = ok;
  if (!ok) return out;
  imu_cal::AccelCalibration<float> fc;
  Proc::Fitter::toFloat(proc->result(), fcfg.g, fc);
  std::vector<imu_cal::GyroHoldView> views;
  for (int i = 0; i < proc->nHolds(); ++i) views.push_back(proc->gyroHoldView(i));
  auto fitter = std::make_unique<imu_cal::GyroAxisFitter<>>();
  fitter->fit(views.data(), (int)views.size(), fc, imu_cal::GyroAxisCfg{}, out.r);
  out.before = rateError(Mat3::Identity(), G);
  out.after = rateError(out.r.S.cast<double>(), G);
  return out;
}

static void campaign(const Case& c, int seeds, double min_accept, double max_after_ratio, double max_worse,
                     double max_after_abs) {
  int n = 0, acc = 0, worse = 0;
  double sum_b = 0, sum_a = 0, worst_after = 0, worst_ratio = 0, worst_delta = -1;
  double fitted_before = 0, fitted_after = 0;
  for (int s = 1; s <= seeds; ++s) {
    const SessionOut o = runSession(c, 1000 + s);
    if (!o.accel_ok) continue;
    ++n;
    if (o.r.ok) ++acc;
    sum_b += o.before; sum_a += o.after;
    if (o.r.ok) {
      worst_after = std::max(worst_after, o.after);
      worst_ratio = std::max(worst_ratio, o.after / std::max(o.before, 1e-9));
      fitted_before += o.before; fitted_after += o.after;
    }
    worst_delta = std::max(worst_delta, o.after - o.before);
    if (o.after > o.before + max_worse) ++worse;
    std::printf("  %-14s seed=%2d %-15s pairs=%2d drop=%d cv %.5f -> %.5f rad  rate err %.4f -> %.4f\n", c.name, s,
                imu_cal::gyroAxisReasonStr(o.r.reason), o.r.pairs, o.r.dropped, (double)o.r.cv_nominal,
                (double)o.r.cv_fit, o.before, o.after);
  }
  char b[256];
  std::snprintf(b, sizeof(b), "%s: %d sessions, %d fitted; mean rate error %.4f -> %.4f", c.name, n, acc,
                n ? sum_b / n : 0.0, n ? sum_a / n : 0.0);
  std::printf("%s\n", b);
  check(n >= seeds - 2, std::string(c.name) + ": accelerometer sessions complete");
  check(acc >= min_accept * n, std::string(c.name) + ": fit accepted often enough");
  if (max_after_ratio > 0) {
    std::snprintf(b, sizeof(b), "%s: accepted fits remove most of the error (worst after/before %.2f)", c.name,
                  worst_ratio);
    check(worst_ratio <= max_after_ratio, b);
  }
  std::snprintf(b, sizeof(b), "%s: never worse than nominal axes by more than %.4f (worst %+.4f)", c.name,
                max_worse, worst_delta);
  check(worse == 0, b);
  if (acc) {
    std::snprintf(b, sizeof(b), "%s: fitted sessions mean rate error %.4f -> %.4f, worst after %.4f <= %.4f", c.name,
                  fitted_before / acc, fitted_after / acc, worst_after, max_after_abs);
    check(worst_after <= max_after_abs, b);
  }
}

int main() {
  testIntegrator();
  testReanchor();
  testRecord();

  const char* seeds_env = std::getenv("GYRO_AXIS_SEEDS");
  const int seeds = seeds_env ? std::max(1, int(std::strtol(seeds_env, nullptr, 10))) : 12;
  // Typical MEMS gyro: ~1% scale and cross-axis errors.
  campaign({"typical", 0.010, 0.007, 1e-4, 0.0, true}, seeds, 0.8, 0.35, 0.0, 0.010);
  // Rate bias drifted ~0.1 deg/s since the stationary capture.
  campaign({"bias-drift", 0.010, 0.007, 1.5e-3, 0.0, true}, seeds, 0.6, 0.5, 0.0, 0.015);
  // Occasional 30-150 ms sample dropouts break the chain.
  campaign({"gaps", 0.010, 0.007, 1e-4, 0.0002, true}, seeds, 0.15, 0.5, 0.0, 0.010);
  // Unsteady hands: twice the wobble, tremor and drift.
  campaign({"adverse-hand", 0.010, 0.007, 1e-4, 0.0, true, 2.0}, seeds, 0.6, 0.5, 0.0, 0.010);
  // Accurate gyro: the fit must not degrade nominal axes.
  campaign({"nominal", 0.0, 0.0, 1e-4, 0.0, true}, seeds, 0.0, -1, 0.003, 0.003);
  // Accelerometer-only recalibration: no integration, nothing fitted.
  {
    const SessionOut o = runSession({"no-chain", 0.01, 0.007, 1e-4, 0.0, false}, 77);
    check(o.accel_ok && !o.r.ok && o.r.reason == imu_cal::GyroAxisReason::NO_DATA &&
              o.r.S.isIdentity(0.0f), "without a gyro chain the matrix stays nominal");
  }

  std::printf("%s (%d failures)\n", g_fail ? "FAILED" : "ALL PASSED", g_fail);
  return g_fail ? 1 : 0;
}
