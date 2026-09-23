/*
  Copyright 2026, Mikhail Grushinskiy

  The AtomS3R marine INS sketches used to carry their own copies of these
  helpers.  Each library helper must reproduce the sketch code it replaced
  bit for bit, so moving it changed no deployed output.  The reference
  implementations below are those sketch copies, verbatim apart from names.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

#include "detrend/AdaptiveWaveDetrender.h"
#include "util/AngleUtils.h"
#include "util/ImuLoopHelpers.h"
#include "util/MagneticHeading.h"
#include "util/QuaternionUtils.h"
#include "wave_dir/WaveDirectionReport.h"

using Vector3f = Eigen::Vector3f;
namespace ins = ocean_imu::ins;

// Arduino.h definition.
#define RAD_TO_DEG 57.295779513082320876798154814105

namespace ref {

static inline float clampf_(float x, float lo, float hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

static inline float wrap360_(float deg) {
  while (deg < 0.0f) deg += 360.0f;
  while (deg >= 360.0f) deg -= 360.0f;
  return deg;
}

static inline float wrap180_(float deg) {
  while (deg <= -180.0f) deg += 360.0f;
  while (deg >   180.0f) deg -= 360.0f;
  return deg;
}

static inline Vector3f quatRotate_(const Eigen::Quaternionf& q, const Vector3f& v) {
  const Vector3f qv(q.x(), q.y(), q.z());
  const Vector3f t = 2.0f * qv.cross(v);
  return v + q.w() * t + qv.cross(t);
}

static inline int waveDirectionSignPolarity_(WaveDirection s) {
  if (s == UNCERTAIN) {
    return 0;
  }

  const int raw = static_cast<int>(s);

  if (raw < 0) return -1;
  if (raw > 0) return +1;

  return 0;
}

static inline int waveDirectionSignRaw_(WaveDirection s) {
  return static_cast<int>(s);
}

static inline bool signedWaveDirectionDeg_(float axis_deg,
                                           int sign_polarity,
                                           float& signed_deg_out)
{
  if (!std::isfinite(axis_deg) || sign_polarity == 0) {
    signed_deg_out = NAN;
    return false;
  }

  signed_deg_out = wrap360_(axis_deg + (sign_polarity < 0 ? 180.0f : 0.0f));
  return true;
}

static inline bool rollPitchHeadingFromQuatBw_(
    const Eigen::Quaternionf& q_bw,
    float& roll_deg_out,
    float& pitch_deg_out,
    float& heading_deg_out) {
  Eigen::Quaternionf q = q_bw;
  const float nq = q.norm();
  if (!(nq > 1e-6f) || !std::isfinite(nq)) return false;
  q.normalize();

  const float x = q.x();
  const float y = q.y();
  const float z = q.z();
  const float w = q.w();

  const float siny_cosp = 2.0f * (w * z + x * y);
  const float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
  const float yaw = atan2f(siny_cosp, cosy_cosp);

  float sinp = 2.0f * (w * y - z * x);
  sinp = clampf_(sinp, -1.0f, 1.0f);
  const float pitch = asinf(sinp);

  const float sinr_cosp = 2.0f * (w * x + y * z);
  const float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
  const float roll = atan2f(sinr_cosp, cosr_cosp);

  roll_deg_out    = wrap180_(roll * RAD_TO_DEG);
  pitch_deg_out   = wrap180_(pitch * RAD_TO_DEG);
  heading_deg_out = wrap360_(yaw * RAD_TO_DEG);
  return true;
}

static inline bool magneticHeadingFromDownAndMagBody_(
    const Vector3f& down_b_unit,
    const Vector3f& mag_b_uT,
    float& heading_deg_out)
{
  heading_deg_out = NAN;
  if (!down_b_unit.allFinite() || !mag_b_uT.allFinite()) return false;
  const Vector3f FWD_B(1.0f, 0.0f, 0.0f);

  Vector3f d = down_b_unit;
  const float dn = d.norm();
  if (!std::isfinite(dn) || !(dn > 1e-6f)) return false;
  d /= dn;

  Vector3f m = mag_b_uT;
  const float mn = m.norm();
  if (!std::isfinite(mn) || !(mn > 1e-6f)) return false;
  m /= mn;

  Vector3f east_b = d.cross(m);
  const float en = east_b.norm();
  if (!std::isfinite(en) || !(en > 1e-6f)) return false;
  east_b /= en;

  Vector3f north_b = east_b.cross(d);
  const float nn = north_b.norm();
  if (!std::isfinite(nn) || !(nn > 1e-6f)) return false;
  north_b /= nn;

  const float e = east_b.dot(FWD_B);
  const float n = north_b.dot(FWD_B);
  if (!(std::hypot(e, n) > 1e-6f)) return false;
  heading_deg_out = wrap360_(atan2f(e, n) * RAD_TO_DEG);
  return std::isfinite(heading_deg_out);
}

// FusionApp members, one struct per sketch family.
struct App {
  float max_dt = 0.0f;  // NLO / PII: 0.05
  static constexpr float LOOP_HZ = 200.0f;
  static constexpr float g_std = 9.80665f;

  uint32_t mag_gate_last_ms_ = 0;
  bool     have_last_sample_us_ = false;
  uint32_t last_sample_us_      = 0;

  bool     rot_inited_    = false;
  float    rot_dpm_filt_  = 0.0f;
  bool     gyro_bias_ok_       = false;
  bool     gyro_bias_learning_ = false;
  Vector3f gyro_bias_ema_      = Vector3f::Zero();

  uint32_t rate_window_ms_   = 0;
  uint32_t imu_sample_count_ = 0;
  uint32_t mag_sample_count_ = 0;
  float    imu_sample_hz_    = 0.0f;
  float    mag_sample_hz_    = 0.0f;

  bool updateMagFreshGate_(bool mag_ok, uint32_t now_ms) {
    constexpr uint32_t kSampleSpacingMs = 35u;

    if (!mag_ok) {
      mag_gate_last_ms_ = 0;
      return false;
    }

    if (mag_gate_last_ms_ == 0) {
      mag_gate_last_ms_ = now_ms;
      return true;
    }

    if ((now_ms - mag_gate_last_ms_) < kSampleSpacingMs) {
      return false;
    }

    mag_gate_last_ms_ = now_ms;
    return true;
  }

  float computeFusionDtFromSampleTimestamp_(uint32_t sample_us) {
    const float dt_nom = 1.0f / LOOP_HZ;
    if (!have_last_sample_us_) {
      have_last_sample_us_ = true;
      last_sample_us_ = sample_us;
      return dt_nom;
    }

    const uint32_t dt_us = sample_us - last_sample_us_;
    last_sample_us_ = sample_us;
    const float dt_s = static_cast<float>(dt_us) * 1.0e-6f;
    if (!(dt_s > 0.0f) || !std::isfinite(dt_s)) return dt_nom;
    if (max_dt > 0.0f && dt_s > max_dt) return dt_nom;
    return dt_s;
  }

  // OU / TFG / PII gyro-bias + ROT block; yaw_rate_from_w maps the corrected
  // rate to the yaw rate (world z for OU / TFG, body z for PII).
  template <typename YawRate>
  void rot(const Vector3f& w_cal_, const Vector3f& a_cal_, float dt_, YawRate yaw_rate_from_w) {
    const bool still =
        (fabsf(a_cal_.norm() - g_std) < 0.12f * g_std) &&
        (w_cal_.norm() < 0.15f);

    gyro_bias_learning_ = still;

    if (still) {
      const float alpha_b = 1.0f - expf(-dt_ / 5.0f);
      if (!gyro_bias_ok_) {
        gyro_bias_ok_  = true;
        gyro_bias_ema_ = w_cal_;
      } else {
        gyro_bias_ema_ += alpha_b * (w_cal_ - gyro_bias_ema_);
      }
    }

    Vector3f w_use = w_cal_;
    if (gyro_bias_ok_) w_use -= gyro_bias_ema_;

    float rot_dpm_meas = yaw_rate_from_w(w_use) * RAD_TO_DEG * 60.0f;
    rot_dpm_meas = clampf_(rot_dpm_meas, -720.0f, 720.0f);

    const float tau_rot = 0.1f;
    const float alpha_r = 1.0f - expf(-dt_ / tau_rot);
    if (!rot_inited_) {
      rot_inited_ = true;
      rot_dpm_filt_ = rot_dpm_meas;
    } else {
      rot_dpm_filt_ += alpha_r * (rot_dpm_meas - rot_dpm_filt_);
    }
  }

  void rates(uint32_t now_ms) {
    if (rate_window_ms_ == 0) rate_window_ms_ = now_ms;
    const uint32_t rate_elapsed_ms = now_ms - rate_window_ms_;
    if (rate_elapsed_ms >= 1000u) {
      const float scale = 1000.0f / static_cast<float>(rate_elapsed_ms);
      imu_sample_hz_ = static_cast<float>(imu_sample_count_) * scale;
      mag_sample_hz_ = static_cast<float>(mag_sample_count_) * scale;
      imu_sample_count_ = 0;
      mag_sample_count_ = 0;
      rate_window_ms_ = now_ms;
    }
  }
};

}  // namespace ref

static int failures = 0;
static long checks = 0;

static bool same(float a, float b) {
  if (std::isnan(a) && std::isnan(b)) return true;
  return std::memcmp(&a, &b, sizeof a) == 0;
}

#define CHECK(cond, ...)                                  \
  do {                                                    \
    ++checks;                                             \
    if (!(cond)) {                                        \
      if (++failures <= 20) {                             \
        std::fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
        std::fprintf(stderr, __VA_ARGS__);                \
        std::fprintf(stderr, "\n");                       \
      }                                                   \
    }                                                     \
  } while (0)

static void testAngles(std::mt19937& rng) {
  std::uniform_real_distribution<float> wide(-2000.0f, 2000.0f);
  const float edges[] = {0.0f, -0.0f, 180.0f, -180.0f, 360.0f, -360.0f, 359.99997f,
                         540.0f, -540.0f, 1e-7f, -1e-7f, 720.0f};
  for (int i = 0; i < 200000; ++i) {
    const float x = i < 12 ? edges[i] : wide(rng);
    CHECK(same(ins::wrap360Deg(x), ref::wrap360_(x)), "wrap360 %.9g", x);
    CHECK(same(ins::wrap180Deg(x), ref::wrap180_(x)), "wrap180 %.9g", x);
    CHECK(same(ins::clampf(x, -720.0f, 720.0f), ref::clampf_(x, -720.0f, 720.0f)), "clamp %.9g", x);
    CHECK(same(ins::radToDeg(x), static_cast<float>(x * RAD_TO_DEG)), "radToDeg %.9g", x);
  }
}

static Eigen::Quaternionf randomQuat(std::mt19937& rng) {
  std::normal_distribution<float> n(0.0f, 1.0f);
  return Eigen::Quaternionf(n(rng), n(rng), n(rng), n(rng));
}

static void testQuaternions(std::mt19937& rng) {
  std::normal_distribution<float> n(0.0f, 20.0f);
  for (int i = 0; i < 100000; ++i) {
    Eigen::Quaternionf q = randomQuat(rng);
    if (i % 3 == 0) q.normalize();
    if (i == 0) q = Eigen::Quaternionf(0.0f, 0.0f, 0.0f, 0.0f);
    if (i == 1) q = Eigen::Quaternionf(NAN, 0.0f, 0.0f, 1.0f);
    if (i == 2) q = Eigen::Quaternionf(0.70710678f, 0.0f, 0.70710678f, 0.0f);  // pitch +90

    const Vector3f v(n(rng), n(rng), n(rng));
    const Vector3f a = ins::quatRotate(q, v);
    const Vector3f b = ref::quatRotate_(q, v);
    for (int k = 0; k < 3; ++k) CHECK(same(a[k], b[k]), "quatRotate %d", i);

    float r1 = 1.0f, p1 = 2.0f, h1 = 3.0f, r2 = 1.0f, p2 = 2.0f, h2 = 3.0f;
    const bool ok1 = ins::rollPitchHeadingFromQuatBw(q, r1, p1, h1);
    const bool ok2 = ref::rollPitchHeadingFromQuatBw_(q, r2, p2, h2);
    CHECK(ok1 == ok2 && same(r1, r2) && same(p1, p2) && same(h1, h2),
          "rollPitchHeading %d: %d/%d %.9g %.9g %.9g vs %.9g %.9g %.9g",
          i, ok1, ok2, r1, p1, h1, r2, p2, h2);
  }
}

static void testMagHeading(std::mt19937& rng) {
  std::normal_distribution<float> n(0.0f, 30.0f);
  const Vector3f specials[] = {Vector3f::Zero(), Vector3f::UnitX(), Vector3f(0, 0, 40),
                               Vector3f(NAN, 0, 1), Vector3f(INFINITY, 0, 1)};
  for (int i = 0; i < 200000; ++i) {
    Vector3f d(n(rng), n(rng), n(rng));
    Vector3f m(n(rng), n(rng), n(rng));
    if (i < 5) d = specials[i];
    else if (i < 10) m = specials[i - 5];
    else if (i < 15) { d = Vector3f::UnitZ(); m = specials[i - 10]; }

    float h1 = 7.0f, h2 = 7.0f;
    const bool ok1 = ins::magneticHeadingFromDownAndMagBody(d, m, h1);
    const bool ok2 = ref::magneticHeadingFromDownAndMagBody_(d, m, h2);
    CHECK(ok1 == ok2 && same(h1, h2), "magHeading %d: %d/%d %.9g vs %.9g", i, ok1, ok2, h1, h2);
  }
}

static void testWaveDirection() {
  const WaveDirection signs[] = {BACKWARD, UNCERTAIN, FORWARD};
  const float axes[] = {NAN, INFINITY, 0.0f, 12.5f, 179.9f, 180.0f, 190.0f, 359.9f, -5.0f};
  for (WaveDirection s : signs) {
    for (float axis : axes) {
      for (bool live : {false, true}) {
        const auto r = wave_direction::WaveDirectionReport::from(axis, s, live);

        // Sketch sequence.
        const float wave_axis_deg_ = axis;
        const bool wave_axis_ok_ = live && std::isfinite(wave_axis_deg_);
        const int wave_sign_raw_ = ref::waveDirectionSignRaw_(s);
        const int wave_sign_polarity_ = ref::waveDirectionSignPolarity_(s);
        const bool wave_sign_ok_ = live && (s != UNCERTAIN);
        float wave_dir_deg_ = 1234.0f;
        bool wave_dir_ok_ = ref::signedWaveDirectionDeg_(wave_axis_deg_, wave_sign_polarity_, wave_dir_deg_);
        wave_dir_ok_ = wave_dir_ok_ && live;
        const float wave_dir_conf_pct_ = wave_dir_ok_ ? 100.0f : (wave_axis_ok_ ? 50.0f : 0.0f);

        CHECK(same(r.axis_deg, wave_axis_deg_) && r.axis_ok == wave_axis_ok_ && r.sign == s &&
              r.sign_raw == wave_sign_raw_ && r.sign_polarity == wave_sign_polarity_ &&
              r.sign_ok == wave_sign_ok_ && same(r.dir_deg, wave_dir_deg_) &&
              r.dir_ok == wave_dir_ok_ && same(r.conf_pct, wave_dir_conf_pct_),
              "wave report sign=%d axis=%g live=%d", static_cast<int>(s), axis, live);
      }
    }
  }

  const wave_direction::WaveDirectionReport cleared{};
  CHECK(std::isnan(cleared.axis_deg) && !cleared.axis_ok && cleared.sign == UNCERTAIN &&
        cleared.sign_raw == 0 && cleared.sign_polarity == 0 && !cleared.sign_ok &&
        std::isnan(cleared.dir_deg) && !cleared.dir_ok && cleared.conf_pct == 0.0f,
        "cleared report differs from the sketch's reset values");
}

// One long device-like stream through every stateful helper, with resets,
// clock wrap, stalls, dropouts and still/moving phases.
static void testLoopState(std::mt19937& rng, float max_dt) {
  ref::App app;
  app.max_dt = max_dt;
  ins::SampleDtTracker dt_tracker(1.0f / ref::App::LOOP_HZ, max_dt);
  ins::MagFreshGate gate;
  ins::SampleRateMeter meter;
  ins::StillGyroBiasEma bias(5.0f, 0.12f, 0.15f);
  ins::RateOfTurnFilter rot_world, rot_body;
  ref::App app_body = app;

  std::uniform_int_distribution<int> jitter(-300, 300);
  std::normal_distribution<float> n(0.0f, 1.0f);

  uint32_t t_us = 0xFFF00000u;  // wraps after ~1 s
  uint32_t t_ms = 1u;
  for (int i = 0; i < 400000; ++i) {
    if (i % 97003 == 0) {
      // reinitImu_() / resetFusion_()
      app.mag_gate_last_ms_ = 0;
      gate.reset();
      app.have_last_sample_us_ = false; app.last_sample_us_ = 0;
      dt_tracker.reset();
      app.rate_window_ms_ = t_ms; app.imu_sample_count_ = 0; app.mag_sample_count_ = 0;
      app.imu_sample_hz_ = 0.0f; app.mag_sample_hz_ = 0.0f;
      meter.reset(t_ms);
      app.rot_inited_ = false; app.rot_dpm_filt_ = 0.0f;
      app.gyro_bias_ok_ = false; app.gyro_bias_learning_ = false; app.gyro_bias_ema_.setZero();
      app_body = app;
      bias.reset(); rot_world.reset(); rot_body.reset();
    }

    uint32_t step_us = static_cast<uint32_t>(5000 + jitter(rng));
    if (i % 5001 == 0) step_us = 80000;   // stall longer than 50 ms
    if (i % 7919 == 0) step_us = 0;       // repeated timestamp
    t_us += step_us;
    t_ms += step_us / 1000u;

    const float d1 = app.computeFusionDtFromSampleTimestamp_(t_us);
    app_body.computeFusionDtFromSampleTimestamp_(t_us);
    const float d2 = dt_tracker.update(t_us);
    CHECK(same(d1, d2), "dt %d: %.9g vs %.9g", i, d1, d2);

    const bool mag_ok = !(i % 3000 < 40);
    const bool f1 = app.updateMagFreshGate_(mag_ok, t_ms);
    app_body.updateMagFreshGate_(mag_ok, t_ms);
    const bool f2 = gate.update(mag_ok, t_ms);
    CHECK(f1 == f2, "mag gate %d", i);

    ++app.imu_sample_count_;
    meter.countImu();
    if (f1) { ++app.mag_sample_count_; meter.countMag(); }
    app.rates(t_ms);
    meter.update(t_ms);
    CHECK(same(app.imu_sample_hz_, meter.imuHz()) && same(app.mag_sample_hz_, meter.magHz()),
          "rates %d", i);

    const bool moving = (i / 20000) % 2 == 1;
    const float s = moving ? 1.5f : 0.01f;
    const Vector3f w(0.002f + s * n(rng) * 0.1f, -0.001f + s * n(rng) * 0.1f, 0.003f + s * n(rng) * 0.2f);
    const Vector3f a(s * n(rng), s * n(rng), -ref::App::g_std + s * n(rng));
    const Eigen::Quaternionf q = randomQuat(rng).normalized();

    app.rot(w, a, d1, [&](const Vector3f& wu) { return ref::quatRotate_(q, wu).z(); });
    app_body.rot(w, a, d1, [](const Vector3f& wu) { return wu.z(); });

    const bool still = bias.update(w, a, ref::App::g_std, d2);
    rot_world.update(ins::quatRotate(q, bias.corrected(w)).z(), d2);
    rot_body.update(bias.corrected(w).z(), d2);

    CHECK(still == app.gyro_bias_learning_ && bias.learning() == app.gyro_bias_learning_ &&
          bias.ok() == app.gyro_bias_ok_, "still %d", i);
    for (int k = 0; k < 3; ++k) CHECK(same(bias.bias()[k], app.gyro_bias_ema_[k]), "bias %d", i);
    CHECK(same(rot_world.dpm(), app.rot_dpm_filt_), "rot world %d: %.9g vs %.9g",
          i, rot_world.dpm(), app.rot_dpm_filt_);
    CHECK(same(rot_body.dpm(), app_body.rot_dpm_filt_), "rot body %d", i);
  }
}

static void testHeaveDetrenderConfig() {
  // Literal block from the NLO / PII sketches' resetFusion_().
  AdaptiveWaveDetrender::Config dcfg{};
  dcfg.init_wave_freq_hz = 0.30f;
  dcfg.min_wave_freq_hz  = 0.02f;
  dcfg.max_wave_freq_hz  = 1.20f;

  dcfg.baseline_cutoff_fraction = 0.25f;
  dcfg.min_baseline_cutoff_hz   = 0.003f;
  dcfg.max_baseline_cutoff_hz   = 0.25f;

  dcfg.freq_smooth_tau_s = 12.0f;
  dcfg.slope_lpf_tau_s   = 0.20f;
  dcfg.slope_rms_tau_s   = 8.0f;

  dcfg.threshold_rms_fraction  = 0.15f;
  dcfg.min_slope_threshold_abs = 0.002f;
  dcfg.max_slope_threshold_abs = 1.0e9f;

  dcfg.startup_hold_s      = 2.0f;
  dcfg.freq_timeout_cycles = 3.0f;

  dcfg.enable_wave_cleanup     = true;
  dcfg.cleanup_cutoff_fraction = 1.0f;
  dcfg.min_cleanup_cutoff_hz   = 0.003f;
  dcfg.max_cleanup_cutoff_hz   = 0.50f;
  dcfg.cleanup_stages          = 2;

  dcfg.min_dt_s = 1.0e-4f;
  dcfg.max_dt_s = 0.25f;
  dcfg.output_abs_limit = 0.0f;

  const AdaptiveWaveDetrender::Config lib = defaultHeaveDetrenderConfig(0.30f);

  // Run both through the detrender: identical configs give identical output.
  AdaptiveWaveDetrender a, b;
  a.setConfig(dcfg); a.reset(0.0f);
  b.setConfig(lib);  b.reset(0.0f);
  for (int i = 0; i < 60000; ++i) {
    const float t = 0.005f * static_cast<float>(i);
    const float x = 0.8f * sinf(2.0f * 3.14159265f * 0.12f * t) + 0.02f * t;
    const auto oa = a.update(x, 0.005f, 0.12f, i > 1000);
    const auto ob = b.update(x, 0.005f, 0.12f, i > 1000);
    CHECK(same(oa.baseline_slow, ob.baseline_slow) && same(oa.wave_raw, ob.wave_raw) &&
          same(oa.wave_clean, ob.wave_clean), "detrender %d", i);
  }

#define FIELD(f) CHECK(same(static_cast<float>(lib.f), static_cast<float>(dcfg.f)), "config field " #f)
  FIELD(init_wave_freq_hz); FIELD(min_wave_freq_hz); FIELD(max_wave_freq_hz);
  FIELD(baseline_cutoff_fraction); FIELD(min_baseline_cutoff_hz); FIELD(max_baseline_cutoff_hz);
  FIELD(freq_smooth_tau_s); FIELD(slope_lpf_tau_s); FIELD(slope_rms_tau_s);
  FIELD(threshold_rms_fraction); FIELD(min_slope_threshold_abs); FIELD(max_slope_threshold_abs);
  FIELD(startup_hold_s); FIELD(freq_timeout_cycles);
  FIELD(enable_wave_cleanup); FIELD(cleanup_cutoff_fraction); FIELD(min_cleanup_cutoff_hz);
  FIELD(max_cleanup_cutoff_hz); FIELD(cleanup_stages);
  FIELD(min_dt_s); FIELD(max_dt_s); FIELD(output_abs_limit);
#undef FIELD
}

int main() {
  std::mt19937 rng(20260923u);
  testAngles(rng);
  testQuaternions(rng);
  testMagHeading(rng);
  testWaveDirection();
  testLoopState(rng, 0.0f);    // OU-II / OU-III / TFG
  testLoopState(rng, 0.05f);   // NLO / PII
  testHeaveDetrenderConfig();

  if (failures) {
    std::fprintf(stderr, "ins-helpers-test: %d of %ld checks FAILED\n", failures, checks);
    return 1;
  }
  std::printf("ins-helpers-test: %ld checks PASS (library helpers bit-identical to the sketch code they replaced)\n",
              checks);
  return 0;
}
