// Heave/displacement drift regression for a shipped INS sketch.
//
// The real sketch (.ino) is compiled unchanged against host stand-ins for the
// Arduino core, M5Unified and the AtomS3R acquisition/UI layers
// (tests/common/ins_sketch_stubs). Its own begin()/tick() run: boot without a
// stored calibration, resetFusion_() configuration, updateFilter_() wiring,
// heave extraction, UI and serial/NMEA output. Only the sensor samples are
// synthetic. This covers the device pipeline that the offline filter
// validation bypasses.
//
// Build with -DINS_SKETCH_FILE='"/absolute/path/to/sketch.ino"',
// -DINS_READY(app)=<bool expression> naming the sketch's fused-output state
// and optionally -DINS_DISPLACEMENT_UP(app)=<Vector3f expression> for sketches
// that expose a three-dimensional z-up displacement.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <Arduino.h>
#include <ArduinoEigenDense.h>

// Test access to the sketch's private state (outputs only; nothing is written).
#define private public
#include INS_SKETCH_FILE
#undef private

namespace {

using V3 = Eigen::Vector3f;
using Q = Eigen::Quaternionf;

constexpr float kDt = 0.005f;            // the sketches' 200 Hz loop period
constexpr float kDeg = 0.017453292519943295f;

// One sinusoid: amplitude, period [s], phase [rad].
struct Tone { float amp, period_s, phase; };
using Tones = std::array<Tone, 3>;

// first_scale multiplies every tone after the first (used for attitude only,
// whose rates are differentiated numerically from the attitude history).
float sumOf(const Tones& tones, double t, int derivative, float first_scale = 1.0f) {
  double v = 0;
  for (size_t i = 0; i < tones.size(); ++i) {
    const auto& c = tones[i];
    if (c.amp == 0) continue;
    const double amp = i == 0 ? c.amp : c.amp * first_scale;
    const double w = 2 * M_PI / c.period_s, x = w * t + c.phase;
    v += derivative == 0 ? amp * std::sin(x) : -amp * w * w * std::sin(x);
  }
  return float(v);
}

struct Scenario {
  const char* name;
  Tones heave_m;                          // vertical displacement, up positive
  Tones surge_m;                          // along-heading displacement
  Tones roll_deg, pitch_deg;
  float chop_onset_s;                     // attitude tones after the first ramp in here
  float duration_s;
  float max_abs_heave_m;                  // bound after warm-up, any instant
  float max_heave_error_frac;             // RMS error / truth RMS (when moving)
  float max_drift_m;                      // |mean(last window) - mean(previous)|
  float max_abs_horizontal_m;             // bound on |x|,|y| when exposed
};

struct Truth {
  Q q_bn;      // body -> NED
  V3 acc_ned;  // kinematic acceleration
  float heave_up;
};

// Wave-following vessel: heave/surge plus roll/pitch about a slowly wandering
// heading. All signals are smooth and exactly differentiable.
Truth truthAt(const Scenario& s, double t) {
  const float yaw = float(1.0 + 0.15 * std::sin(2 * M_PI * t / 97.0));
  // Smooth 5 s ramp of the additional attitude tones after chop_onset_s.
  const double u = std::clamp((t - s.chop_onset_s) / 5.0, 0.0, 1.0);
  const float chop = float(u * u * (3 - 2 * u));
  const float roll = sumOf(s.roll_deg, t, 0, chop) * kDeg;
  const float pitch = sumOf(s.pitch_deg, t, 0, chop) * kDeg;
  Truth out;
  out.q_bn = (Eigen::AngleAxisf(yaw, V3::UnitZ()) * Eigen::AngleAxisf(pitch, V3::UnitY()) *
              Eigen::AngleAxisf(roll, V3::UnitX())).normalized();
  out.heave_up = sumOf(s.heave_m, t, 0);
  const float az = -sumOf(s.heave_m, t, 2);  // NED z is down
  const float ax = sumOf(s.surge_m, t, 2);
  out.acc_ned = V3(ax * std::cos(yaw), ax * std::sin(yaw), az);
  return out;
}

struct Sensors {
  const Scenario* scenario = nullptr;
  std::mt19937 rng{20261012};
  std::normal_distribution<float> n{0.0f, 1.0f};
  V3 gyro_bias = V3(0.0020f, -0.0015f, 0.0010f);   // rad/s residual after calibration
  V3 accel_bias = V3(0.015f, -0.010f, 0.020f);     // m/s^2 residual after calibration
  V3 field_ned = V3(20.0f, -1.5f, 45.0f);          // uT
  atoms3r_ical::MagAcquisition mag{};
  uint32_t ticks = 0;
  float heave_up = 0.0f;

  V3 noise(float sigma) { return sigma * V3(n(rng), n(rng), n(rng)); }

  void fill(uint32_t sample_us, atoms3r_ical::ImuSample& out) {
    const double t = double(ticks++) * kDt;
    const Truth now = truthAt(*scenario, t);
    const Truth next = truthAt(*scenario, t + kDt);
    heave_up = now.heave_up;
    const Q dq = now.q_bn.conjugate() * next.q_bn;
    const Eigen::AngleAxisf aa(dq);
    const V3 w_true = aa.axis() * (aa.angle() / kDt);
    const float g = atoms3r_ical::ImuCalCfg::g_cal_local;
    out.a = now.q_bn.conjugate() * (now.acc_ned - V3(0, 0, g)) + accel_bias + noise(0.02f);
    out.w = w_true + gyro_bias + noise(0.0015f);
    out.tempC = 35.0f;
    // 25 Hz magnetometer: a new acquisition every eighth frame, held between.
    if (ticks % 8 == 1) {
      mag.acquired(millis(), sample_us, sample_us, sample_us + 300);
      mag_body = now.q_bn.conjugate() * field_ned + noise(0.3f);
    }
    out.m = mag.fresh(millis()) ? mag_body : V3::Constant(NAN);
    out.mag_timing = mag;
  }
  V3 mag_body = V3::Constant(NAN);
};

Sensors* g_sensors = nullptr;
bool sampleSource(uint32_t sample_us, atoms3r_ical::ImuSample& out) {
  g_sensors->fill(sample_us, out);
  return true;
}

int failures = 0;
void check(bool ok, const char* scenario, const char* what, double value, double limit) {
  std::printf("  %-34s %-28s %10.4f  (limit %.4f)  %s\n", scenario, what, value, limit, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

void run(const Scenario& s) {
  ins_sketch_host::clock_us() = 1000000;
  Sensors sensors;
  sensors.scenario = &s;
  g_sensors = &sensors;
  atoms3r_ical::sampleSource() = &sampleSource;

  auto app = std::make_unique<FusionApp>();
  app->begin();

  const float warmup_s = s.duration_s * 0.5f;
  const float window_s = 60.0f;
  double err2 = 0, truth2 = 0, max_heave = 0;
#ifdef INS_DISPLACEMENT_UP
  double max_horizontal = 0;
#endif
  double sum_last = 0, sum_prev = 0;
  long n = 0, n_last = 0, n_prev = 0;
  bool finite = true;
  const long steps = long(s.duration_s / kDt);
  for (long i = 0; i < steps; ++i) {
    app->tick();
    const float t = float(i) * kDt;
    const float heave = app->heave_m_;
    if (!std::isfinite(heave)) finite = false;
    if (t < warmup_s) continue;
    ++n;
    max_heave = std::max(max_heave, double(std::abs(heave)));
    err2 += double(heave - sensors.heave_up) * double(heave - sensors.heave_up);
    truth2 += double(sensors.heave_up) * double(sensors.heave_up);
    if (t >= s.duration_s - window_s) { sum_last += heave; ++n_last; }
    else if (t >= s.duration_s - 2 * window_s) { sum_prev += heave; ++n_prev; }
#ifdef INS_DISPLACEMENT_UP
    const V3 d = INS_DISPLACEMENT_UP((*app));
    if (!d.allFinite()) finite = false;
    max_horizontal = std::max({max_horizontal, double(std::abs(d.x())), double(std::abs(d.y()))});
#endif
  }
  check(finite, s.name, "finite outputs", finite ? 1 : 0, 1);
  check(max_heave <= s.max_abs_heave_m, s.name, "max |heave| after warm-up [m]", max_heave, s.max_abs_heave_m);
  const double drift = std::abs(sum_last / std::max(1L, n_last) - sum_prev / std::max(1L, n_prev));
  check(drift <= s.max_drift_m, s.name, "window mean heave drift [m]", drift, s.max_drift_m);
  if (truth2 > 0) {
    const double frac = std::sqrt(err2 / truth2);
    check(frac <= s.max_heave_error_frac, s.name, "heave RMS error / truth RMS", frac, s.max_heave_error_frac);
  }
#ifdef INS_DISPLACEMENT_UP
  check(max_horizontal <= s.max_abs_horizontal_m, s.name, "max |x|,|y| after warm-up [m]",
        max_horizontal, s.max_abs_horizontal_m);
#endif
  const bool ready = INS_READY((*app));
  check(ready, s.name, "fused INS output live at end", ready ? 1 : 0, 1);
  check(app->heading_valid_, s.name, "heading valid at end", app->heading_valid_ ? 1 : 0, 1);
}

}  // namespace

int main() {
  constexpr Tone none{0, 1, 0};
  // name, heave [m], surge [m], roll [deg], pitch [deg], chop onset [s], duration [s],
  // max |heave|, heave error fraction, drift, max horizontal.
  // The limits are common to all INS sketches: they detect drift, runaway and
  // a fused output that never goes live (a zero output has error fraction 1),
  // not per-filter tracking quality, which the offline validation covers.
  const Scenario scenarios[] = {
    {"calm water at the dock",
     {none, none, none}, {none, none, none},
     {Tone{0.5f, 4.0f, 0}, none, none}, {Tone{0.3f, 3.0f, 0.7f}, none, none},
     0, 600, 0.30f, 1e9f, 0.05f, 0.75f},
    {"swell H=1 m T=8 s",
     {Tone{0.5f, 8.0f, 0}, none, none}, {Tone{0.15f, 8.0f, 1.57f}, none, none},
     {Tone{3.0f, 5.0f, 0}, none, none}, {Tone{2.0f, 4.0f, 0.7f}, none, none},
     0, 600, 0.90f, 0.95f, 0.10f, 1.00f},
    {"moderate sea H=2 m T=6 s",
     {Tone{1.0f, 6.0f, 0}, none, none}, {Tone{0.3f, 6.0f, 1.57f}, none, none},
     {Tone{8.0f, 4.0f, 0}, none, none}, {Tone{5.0f, 3.5f, 0.7f}, none, none},
     0, 600, 1.60f, 0.95f, 0.15f, 1.60f},
    // Small boat in steep wind chop over swell: brisk roll/pitch with sharp
    // angular accelerations and +/-0.25 g vertical jolts.
    {"small boat in short chop",
     {Tone{0.6f, 7.0f, 0}, Tone{0.06f, 1.6f, 0.4f}, Tone{0.03f, 1.1f, 2.0f}},
     {Tone{0.2f, 7.0f, 1.57f}, Tone{0.03f, 1.6f, 1.0f}, none},
     {Tone{12.0f, 3.0f, 0}, Tone{6.0f, 1.3f, 0.9f}, Tone{3.0f, 0.8f, 2.1f}},
     {Tone{6.0f, 2.5f, 0.7f}, Tone{4.0f, 1.1f, 1.7f}, Tone{2.0f, 0.7f, 0.3f}},
     0, 600, 1.20f, 0.95f, 0.15f, 1.20f},
    // Live in swell first, then brisk chop: the estimator must keep its
    // gravity and bias aiding through the chop rather than coast on gyro.
    {"swell, then chop after 300 s",
     {Tone{0.5f, 8.0f, 0}, none, none}, {Tone{0.15f, 8.0f, 1.57f}, none, none},
     {Tone{4.0f, 5.0f, 0}, Tone{8.0f, 1.3f, 0.9f}, Tone{4.0f, 0.8f, 2.1f}},
     {Tone{3.0f, 4.0f, 0.7f}, Tone{5.0f, 1.1f, 1.7f}, Tone{2.5f, 0.7f, 0.3f}},
     300, 600, 0.90f, 0.95f, 0.10f, 1.00f},
    {"hand rocking 20 deg / 2 s",
     {none, none, none}, {none, none, none},
     {Tone{20.0f, 2.0f, 0}, none, none}, {Tone{5.0f, 1.5f, 0.7f}, none, none},
     0, 600, 0.30f, 1e9f, 0.10f, 0.75f},
  };
  std::printf("INS sketch heave drift: %s\n", INS_SKETCH_FILE);
  for (const auto& s : scenarios) run(s);
  if (failures) {
    std::printf("INS sketch heave drift: %d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("INS sketch heave drift: PASS\n");
  return 0;
}
