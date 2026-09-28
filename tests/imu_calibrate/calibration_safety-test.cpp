// Copyright 2026, Mikhail Grushinskiy
// Deterministic gyro qualification, registered-frame mag, capture and storage regressions.
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/MagCalCapture.h"
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include <random>
#include <limits>

using V = Eigen::Vector3f;
using M = Eigen::Matrix3f;
using namespace atoms3r_ical;
static int checks = 0, failures = 0;
static void check(bool ok, const char* what) {
  ++checks; if (!ok) { ++failures; std::fprintf(stderr, "FAIL: %s\n", what); }
}

struct FaultKv {
  std::map<std::string, std::vector<uint8_t>> bytes;
  int writes = 0, reads = 0, corrupt = 0, short_writes = 0, fail_read_number = -1;
  bool drop = false, fail_remove = false, pretend_remove = false;
  size_t getBytesLength(const char* key) { return bytes.count(key) ? bytes[key].size() : 0; }
  size_t getBytes(const char* key, void* dst, size_t n) {
    if (++reads == fail_read_number || !bytes.count(key) || bytes[key].size() != n) return 0;
    memcpy(dst, bytes[key].data(), n); return n;
  }
  size_t putBytes(const char* key, const void* src, size_t n) {
    ++writes;
    if (drop) return n;
    if (short_writes > 0) { --short_writes; n /= 2; }
    bytes[key].assign((const uint8_t*)src, (const uint8_t*)src + n);
    if (corrupt > 0) { --corrupt; bytes[key][n/2] ^= 0x80; }
    return n;
  }
  bool remove(const char* key) {
    if (fail_remove) return false;
    if (!pretend_remove) bytes.erase(key);
    return true;
  }
};
using Store = ImuCalStoreT<FaultKv>;
using GC = imu_cal::GyroCalibration<float>;

static GC gyro(float lo, float span, const V& slope, float noise = 0, int count = 320) {
  imu_cal::GyroCalibrator<float,400> cal;
  std::mt19937 rng(734);
  std::normal_distribution<float> normal(0, noise);
  for (int i = 0; i < count; ++i) {
    const float t = lo + span * i / (count-1);
    const V w = V(.01f,-.008f,.003f) + slope * (t-(lo+span/2)) + V(normal(rng),normal(rng),normal(rng));
    check(cal.addSample(w, V(0,0,cal.g), t), "gyro stationary sample accepted");
  }
  GC out;
  check(cal.fit(out) && out.ok, "gyro bias fit succeeds");
  return out;
}
static ImuCalBlobV4 gyroBlob(const GC& g) {
  ImuCalBlobV4 b;
  fillGyroFromFit(b, g);
  return Store::sealed_(b);
}
static void testGyro() {
  for (float span : {0.f, .00001f, .01f}) {
    GC c = gyro(25, span, V(span > 0 ? .02f : 0,0,0), .00003f);
    check(c.thermal == imu_cal::GyroThermal::UNLEARNED, "constant/narrow span is bias only");
    check(c.biasT.k.isZero(0), "narrow temperature has EXACT zero slope");
    check(std::fabs(c.biasT.T0 - (25+span/2)) < 1e-5, "reference temperature is measured");
    check((c.biasT.b0 - V(.01f,-.008f,.003f)).norm() < 1e-5, "all stationary samples estimate reference bias");
    Store store; ImuCalBlobV4 rb;
    check(store.saveVerified(gyroBlob(c), rb), "bias-only float blob validates and round-trips");
    RuntimeCals rc; rc.rebuildFromBlob(rb);
    for (float t : {-100.f, 35.f, 100.f, std::numeric_limits<float>::quiet_NaN()})
      check((rc.applyGyro(V::Zero(), t) + c.biasT.b0).norm() < 1e-8, "narrow session never fabricates thermal rate");
  }
  // The authoritative 2 C span is an eligibility boundary, not a promise
  // that a slope is identifiable. Just below it remains SPAN_SMALL; at the
  // boundary the fit proceeds to the independent information/bin gates.
  GC edge = gyro(24, 1.99f, V::Zero(), .00003f);
  check(edge.thermal == imu_cal::GyroThermal::UNLEARNED &&
        edge.thermal_reason == imu_cal::GyroThermalReason::SPAN_SMALL,
        "span below 2 C is not eligible for gyro thermal fit");
  edge = gyro(24, 2.0f, V::Zero(), .00003f);
  check(edge.thermal == imu_cal::GyroThermal::UNLEARNED &&
        edge.thermal_reason == imu_cal::GyroThermalReason::INFORMATION_LOW,
        "2 C boundary passes span eligibility and remains subject to information gate");
  check(edge.biasT.k.isZero(0) &&
        (edge.biasT.b0 - V(.01f,-.008f,.003f)).norm() < 1e-5,
        "rejected 2 C thermal fit retains stationary gyro bias");

  const V slope(.00035f,-.00021f,.00017f);
  GC c = gyro(5, 40, slope, .0002f);
  check(c.thermal == imu_cal::GyroThermal::LEARNED, "wide-temperature slope qualifies");
  check((c.biasT.k-slope).cwiseAbs().maxCoeff() < 3e-6, "wide-temperature slope accuracy");
  check(c.biasT.T0 > 24.99f && c.biasT.T0 < 25.01f, "regression centered on measured bin temperatures");
  check(c.slope_sigma.maxCoeff() < imu_cal::GyroThermalLimits::max_slope_sigma, "slope uncertainty is qualified");
  check(c.temperature_information >= imu_cal::GyroThermalLimits::min_information, "temperature information qualified");
  Store store; ImuCalBlobV4 rb;
  check(store.saveVerified(gyroBlob(c), rb) && validateBlob(rb), "learned float blob save/readback");
  RuntimeCals rc; rc.rebuildFromBlob(rb);
  for (float t : {-100.f, 5.f, 25.f, 45.f, 100.f})
    check((rc.applyGyro(V(.01f,0,0),t) - c.apply(V(.01f,0,0),t)).norm() < 1e-8, "runtime matches fitted float coefficients");
  check((rc.applyGyro(V::Zero(),100)-rc.applyGyro(V::Zero(),c.biasT.T_hi)).norm() < 1e-8, "upper extrapolation clamped");
  check((rc.applyGyro(V::Zero(),-100)-rc.applyGyro(V::Zero(),c.biasT.T_lo)).norm() < 1e-8, "lower extrapolation clamped");
  check(std::fabs((c.temp_lo-c.biasT.T_lo)-2) < 1e-5 && std::fabs((c.biasT.T_hi-c.temp_hi)-2) < 1e-5,
        "extrapolation margin is two degrees");
  ImuCalBlobV4 corrupt = rb;
  corrupt.gyro_T_hi += 10; corrupt = Store::sealed_(corrupt);
  check(!validateBlob(corrupt), "outer CRC cannot bless unbound gyro metadata");
  corrupt.gyro_coeff_crc = gyroCoeffCrc(corrupt); corrupt = Store::sealed_(corrupt);
  check(!validateBlob(corrupt), "even re-bound oversized extrapolation rejected");
  rc.rebuildFromBlob(corrupt);
  check(!rc.gyr.ok && rc.applyGyro(V::Ones(),100).isApprox(V::Ones()), "runtime rejects invalid gyro metadata");
  corrupt = rb; corrupt.gyro_k_sigma[0] = .002f; corrupt.gyro_coeff_crc = gyroCoeffCrc(corrupt);
  check(!validateBlob(Store::sealed_(corrupt)), "unqualified stored slope uncertainty rejected");
  corrupt = rb; corrupt.version = 3; corrupt.crc = computeBlobCrc(corrupt);
  check(!validateBlob(corrupt), "old version never silently reinterpreted");
  Store legacy; legacy.kv.putBytes("blob_m5v3", &corrupt, sizeof(corrupt));
  check(!legacy.load(rb), "legacy key is not loaded");

  c = gyro(20,10,V(.002f,0,0), .00001f);
  check(c.thermal == imu_cal::GyroThermal::UNLEARNED && c.biasT.k.isZero(0), "implausible slope rejected without losing stationary bias");
  check(c.thermal_reason == imu_cal::GyroThermalReason::SLOPE_IMPLAUSIBLE, "plausibility rejection recorded");
  c = gyro(5,40,V::Zero(), .006f);
  check(c.thermal == imu_cal::GyroThermal::UNLEARNED && c.thermal_reason == imu_cal::GyroThermalReason::INFORMATION_LOW,
        "wide span alone cannot qualify a noisy slope");
  imu_cal::GyroCalibrator<float,400> sparse;
  for(int i=0;i<320;++i) sparse.addSample(V(.01f,0,0),V(0,0,sparse.g),i<160?20.f:40.f);
  check(sparse.fit(c) && c.thermal == imu_cal::GyroThermal::UNLEARNED, "two populated temperature points are insufficient");
  sparse.clear();
  for(int i=0;i<320;++i) sparse.addSample(V(.01f,0,0),V(0,0,sparse.g),i<319?25.f:60.f);
  check(sparse.fit(c) && c.thermal == imu_cal::GyroThermal::UNLEARNED, "one temperature spike cannot create excitation");
  sparse.clear(); c = gyro(5,40,slope);
  check(!sparse.fit(c) && !c.ok && c.biasT.k.isZero(0), "failed fit clears previous slope");
  std::printf("gyro: all narrow fits k=0; wide slope validated with bounded float runtime\n");
}

static V sphere(int i, int n) {
  const float z = 1-2*(i+.5f)/n, p = 2.39996323f*i;
  return V(std::sqrt(1-z*z)*std::cos(p),std::sqrt(1-z*z)*std::sin(p),z);
}
static void testMag(const char* name, const M& distortion, const V& bias, bool noisy) {
  imu_cal::MagCalibrator<float,400> cal;
  std::mt19937 rng(512); std::normal_distribution<float> noise(0,.05f);
  for (int i=0;i<400;++i) {
    V raw = 50*distortion*sphere(i,400)+bias;
    if (noisy) {
      raw += V(noise(rng),noise(rng),noise(rng));
      if (i%20==0) raw += V(12,-9,7); // 5% gross outliers, within sensor sanity limits
    }
    check(cal.addSample(raw), "mag sample passes unchanged sanity gate");
  }
  imu_cal::MagCalibration<float> c; imu_cal::FitFail why;
  bool ok = cal.fit(c,3,.15f,1e-6f,&why);
  check(ok && c.ok, "mag robust fit succeeds");
  if (!ok) { std::printf("mag %s: %s\n",name,imu_cal::fitFailStr(why)); return; }
  check((c.A-c.A.transpose()).norm() < 1e-6, "soft iron correction is symmetric");
  Eigen::SelfAdjointEigenSolver<M> es(c.A);
  check(es.eigenvalues().minCoeff()>0, "soft iron correction is positive definite");
  ImuCalBlobV4 b; b.mag_ok=1; mat_to_rowmajor9_(c.A,b.mag_A);
  for(int j=0;j<3;++j) b.mag_b[j]=c.b[j];
  b.mag_field_uT=c.field_uT; b.mag_rms=c.rms;
  Store store; ImuCalBlobV4 rb;
  check(store.saveVerified(b,rb), "direction-sensitive mag coefficients survive serialization");
  RuntimeCals rc; rc.rebuildFromBlob(rb);
  float heading_max=0, direction_max=0, norm_max=0;
  for (int i=0;i<360;++i) {
    const float p=i*float(M_PI)/180;
    const V u(std::cos(p),std::sin(p),0);
    const V corrected = rc.applyMag((50*distortion*u+bias).eval());
    const float heading=std::fabs(std::remainder(std::atan2(corrected.y(),corrected.x())-p,2*float(M_PI)))*180/float(M_PI);
    heading_max=std::fmax(heading_max,heading);
    norm_max=std::fmax(norm_max,std::fabs(corrected.norm()-c.field_uT));
    const V truth = sphere(i,360), v=rc.applyMag((50*distortion*truth+bias).eval()).normalized();
    // atan2 is stable at sub-millidegree angles where float acos(dot) is not.
    const float angle=std::atan2(v.cross(truth).norm(),v.dot(truth))*180/float(M_PI);
    direction_max=std::fmax(direction_max,angle);
  }
  std::printf("mag %s: heading_max=%.6f deg direction_max=%.6f deg norm_max=%.6f uT rms=%.6f\n",
              name,heading_max,direction_max,norm_max,c.rms);
  check(heading_max < (noisy?.25f:.02f), "corrected heading, not just norm, is accurate");
  check(direction_max < (noisy?.25f:.02f), "corrected 3-D direction is accurate");
  check(norm_max < (noisy?.2f:.01f), "field magnitude matches fitted radius on held-out clean directions");
  // An unknown physical sensor rotation is not observable from norms. The
  // fixture here is SPD in the already registered frame; no such rotation is fitted.
}

static V changing(int i, bool planar=false) {
  const float p=.19f*i,z=planar?0.f:.95f*std::sin(.071f*i);
  return 50*V(std::sqrt(1-z*z)*std::cos(p),std::sqrt(1-z*z)*std::sin(p),z);
}
static void testCapture() {
  using Status=imu_cal::MagCaptureStatus;
  for (int mode=0;mode<3;++mode) {
    imu_cal::MagCalibrator<float,400> cal;
    imu_cal::MagCapture<float,400> capture(cal);
    // Also test uint32 millisecond rollover using exactly the same state machine.
    const uint32_t start = mode==2 ? 0xfffff000u : 0u;
    capture.begin(start);
    Status st=Status::CAPTURING;
    int late=0;
    for(int i=0;i<=563;++i) {
      const V v=changing(i,mode==1 && i<400);
      st=capture.update(start+uint32_t(80*i),&v);
      if(i<563) check(st==Status::CAPTURING, "changing MAG remains live before 45 seconds, including saturation");
      if(i==400) check(cal.buf.n==400, "existing 400-sample buffer saturated without growth");
    }
    check(st==Status::READY && capture.coverage().ok, "45-second capture succeeds, including late 3-D coverage");
    check(cal.buf.n==400 && capture.observations()==564, "retention count and stream count are independent");
    if(mode==1) {
      for(int i=0;i<cal.buf.n;++i) late += std::fabs(cal.buf.v[i].z())>1e-3f ? 1:0;
      check(late>=80, "reservoir retains substantial late 3-D observations");
      imu_cal::MagCalibration<float> m;
      check(cal.fit(m), "late coverage also passes the unchanged ellipsoid fit gates");
    }
    std::printf("capture mode=%d: retained=%d seen=%u late3D=%d det=%.6f\n",mode,cal.buf.n,capture.observations(),late,capture.coverage().determinant);
  }
  imu_cal::MagCalibrator<float,400> cal;
  imu_cal::MagCapture<float,400> capture(cal);
  V frozen(40,10,20);
  capture.begin(0); capture.update(0,&frozen);
  check(capture.update(12001,&frozen)==Status::STALE, "frozen magnetometer fails even with successful reads");
  capture.begin(0);
  check(capture.update(12001,nullptr)==Status::STALE, "missing readings still reach stale timeout");
  V invalid(NAN,0,0); capture.begin(0); capture.update(100,&invalid);
  check(capture.update(12001,&invalid)==Status::STALE, "invalid data cannot keep capture alive");
  capture.begin(0);
  for(int i=0;i<=563;++i) { V v=changing(i,true); capture.update(80*i,&v); }
  check(capture.update(45041,nullptr)==Status::COVERAGE_LOW, "planar motion fails 3-D coverage");
  capture.begin(0);
  for(int i=0;i<=400;++i) { V v=changing(i); capture.update(80*i,&v); frozen=v; }
  check(capture.update(44001,&frozen)==Status::STALE, "true freeze after saturation still fails");
}

static void testStorageAndMissingTemp() {
  ImuCalBlobV4 old=gyroBlob(gyro(5,40,V(.00035f,-.00021f,.00017f)));
  old.accel_ok=1; old.accel_g=ImuCalCfg::g_cal_local; old.accel_T0=25;
  old.accel_S[0]=old.accel_S[4]=old.accel_S[8]=1;
  old.accel_b0[0]=.1f; old.accel_k[0]=.003f;
  old.accel_T_lo=20; old.accel_T_hi=30; old.accel_coeff_crc=accelCoeffCrc(old);
  old=Store::sealed_(old);
  RuntimeCals rc; rc.rebuildFromBlob(old);
  for(float t : {NAN,INFINITY,-INFINITY}) {
    check((rc.applyGyro(V::Zero(),t)+rc.gyr.biasT.b0).norm()<1e-8, "nonfinite gyro temperature uses reference bias");
    check((rc.applyAccel(V::Zero(),t)+rc.acc.biasT.b0).norm()<1e-8, "nonfinite accel temperature uses reference bias");
  }
  check((rc.applyAccel(V::Zero(),35)-rc.applyAccel(V::Zero(),NAN)).norm()>.01,
        "missing-temp regression would detect a fabricated 35 C input");
  ImuCalBlobV4 cand=old; cand.gyro_b0[0]+=.001f; cand.gyro_coeff_crc=gyroCoeffCrc(cand); cand=Store::sealed_(cand);
  for(int fault=0;fault<5;++fault) {
    Store store; ImuCalBlobV4 rb, after;
    check(store.saveVerified(old,rb), "storage fixture saves valid old calibration");
    if(fault==0) store.kv.drop=true;
    if(fault==1) store.kv.corrupt=1;
    if(fault==2) store.kv.corrupt=2;
    if(fault==3) store.kv.short_writes=2;
    if(fault==4) store.kv.fail_read_number=store.kv.reads+2; // candidate verification read only
    check(!store.saveVerified(cand,rb), "failed candidate never reported saved");
    if(fault==0 || fault==1 || fault==4) {
      check(store.lastSaveStatus()==Store::SaveStatus::PREVIOUS_RETAINED, "previous-retained status requires verified old bytes");
      check(store.load(after) && Store::sameBytes(after,old), "failed write really preserves previous calibration");
    } else {
      check(store.lastSaveStatus()==Store::SaveStatus::RECOVERY_FAILED, "rollback corruption/short write explicitly reported");
      check(!store.load(after), "lost previous calibration is not falsely claimed restored");
    }
  }
  for(bool pretend : {false,true}) {
    Store store; ImuCalBlobV4 rb;
    store.kv.corrupt=1; store.kv.fail_remove=!pretend; store.kv.pretend_remove=pretend;
    check(!store.saveVerified(cand,rb) && store.lastSaveStatus()==Store::SaveStatus::RECOVERY_FAILED,
          "failed and falsely successful cleanup detected");
  }
  Store store; ImuCalBlobV4 rb;
  store.saveVerified(old,rb); int writes=store.kv.writes;
  store.kv.fail_read_number=store.kv.reads+1;
  check(!store.saveVerified(cand,rb) && store.kv.writes==writes &&
        store.lastSaveStatus()==Store::SaveStatus::PREVIOUS_UNREADABLE,
        "failed previous read cannot overwrite a possibly valid calibration");
  cand.gyro_k[0]=NAN;
  check(!store.saveVerified(cand,rb) && store.kv.writes==writes, "invalid candidate cannot disturb previous bytes");
  std::printf("storage: corrupt, dropped, short, failed-read, rollback and cleanup faults checked\n");
}

int main() {
  testGyro();
  M diag=V(1.2f,.9f,1.05f).asDiagonal();
  M symmetric; symmetric<<1.2,.18,.09,.18,.9,.12,.09,.12,1.05;
  testMag("diagonal",diag,V::Zero(),false);
  testMag("symmetric-cross",symmetric,V::Zero(),false);
  testMag("hard-iron",M::Identity(),V(8,-5,3),false);
  testMag("combined",symmetric,V(8,-5,3),false);
  testMag("noise-outliers",symmetric,V(8,-5,3),true);
  testCapture(); testStorageAndMissingTemp();
  std::printf("calibration_safety-test: %d/%d checks passed\n",checks-failures,checks);
  return failures?1:0;
}
