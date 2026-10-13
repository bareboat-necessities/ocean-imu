// Host run of sensors/imu_basic/atomS3R_imu_fifo_diagnostics against a
// simulated BMI270: a 6 kB header-mode FIFO filled at 200 Hz on a sensor clock
// that runs 0.8 % slow (as measured on an AtomS3R), read over 400 kHz I2C
// (each register byte sampled as it is clocked out), and a forced-mode BMM150
// whose readout every 40 ms carries the field of a fixed time before it. The simulated sensors have
// known errors (accelerometer and gyro scale, gyro bias, field strength).
// The board rests in seven poses (8 s each, 2.5 s turns between them,
// including a turn about the vertical and a flip), then rocks like a hand
// test for the magnetometer-delay fit. The sketch's own setup()/loop() run
// unchanged, with its display bursts costing host time and one tap.
//
//   fifo_diagnostics_sketch-test [good|bad]
//
// good: small errors, every timing and accuracy criterion must pass and each
//       error must be measured at its size.
// bad:  errors beyond the sketch's limits, each accuracy criterion must fail;
//       the place is moved over serial ("L ...") and the gyro noise mode
//       toggled ("n") first.
//
// Built against tests/common/ins_sketch_stubs by tests/imu_calibrate/Makefile.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <Arduino.h>
#include <ArduinoEigenDense.h>

#include FIFO_SKETCH_FILE

namespace {

using V3 = Eigen::Vector3f;
using Q = Eigen::Quaternionf;
constexpr double kPi = 3.14159265358979323846;
constexpr double kPpm = 8000.0;  // the AtomS3R measured about +8000 (frames every 5.04 ms)
constexpr double kTickUs = 39.0625 * (1.0 + kPpm * 1e-6);
constexpr double kMagDelayUs = 26000.0;  // readout time minus the time the field belongs to
constexpr double kStillS = 8.0, kTurnS = 2.5, kTourStartS = 3.0, kRockS = 110.0;

struct SensorErrors {
  float accel_scale, gyro_scale, field_scale;
  V3 gyro_bias_dps;
};
// good: +3 mg, +1 % gyro scale, |bias| 0.099 deg/s, field +2 %.
const SensorErrors kGood{1.003f, 1.01f, 1.02f, V3(0.05f, -0.08f, 0.03f)};
// bad: +20 mg, +3 % gyro scale, |bias| 0.45 deg/s, field +8 %.
const SensorErrors kBad{1.02f, 1.03f, 1.08f, V3(0.3f, -0.3f, 0.15f)};
SensorErrors g_err = kGood;

Q rotX(double deg) { return Q(Eigen::AngleAxisf(float(deg * kPi / 180), V3::UnitX())); }
Q rotY(double deg) { return Q(Eigen::AngleAxisf(float(deg * kPi / 180), V3::UnitY())); }
Q rotZ(double deg) { return Q(Eigen::AngleAxisf(float(deg * kPi / 180), V3::UnitZ())); }

// Body to NED.
std::vector<Q> poses() {
  const Q p0 = rotZ(46.0);
  const Q p2 = p0 * rotX(60) * rotY(-45);
  return {p0, p0 * rotX(60), p2, rotZ(90) * p2, p0 * rotX(180), p0 * rotX(90), p0 * rotY(80), p0};
}
const std::vector<Q> kPoses = poses();
const double kTourEndS = kTourStartS + double(kPoses.size()) * (kStillS + kTurnS);
const double kDurationS = kTourEndS + kRockS;

Q attitudeAt(double t_s) {
  if (t_s < kTourStartS) return kPoses.front();
  if (t_s < kTourEndS) {
    const double u = t_s - kTourStartS;
    const size_t i = size_t(u / (kStillS + kTurnS));
    const double in = u - double(i) * (kStillS + kTurnS);
    if (in < kStillS || i + 1 >= kPoses.size()) return kPoses[i];
    double s = (in - kStillS) / kTurnS;
    s = s * s * (3 - 2 * s);
    return kPoses[i].slerp(float(s), kPoses[i + 1]).normalized();
  }
  // Rocking about the last (level) pose, eased in.
  const double r = t_s - kTourEndS;
  const double u = std::min(1.0, r / 5.0);
  const float ramp = float(u * u * (3 - 2 * u));
  const float roll = ramp * 0.5f * float(std::sin(2 * kPi * 1.0 * r));
  const float pitch = ramp * 0.3f * float(std::sin(2 * kPi * 0.7 * r + 1.0));
  const float yaw = ramp * 0.4f * float(std::sin(2 * kPi * 0.05 * r));
  return (kPoses.back() * rotZ(yaw * 180 / kPi) * rotY(pitch * 180 / kPi) * rotX(roll * 180 / kPi)).normalized();
}

// The field the sketch should expect, from its own place and date.
V3 fieldNed() {
  static double lat = NAN, lon = NAN;
  static V3 f = V3::Zero();
  if (site_lat != lat || site_lon != lon) {
    const imu_cal::GeomagField g = imu_cal::Wmm2025::field(site_lat, site_lon, site_h_m, site_year);
    f = V3(float(g.x_nT), float(g.y_nT), float(g.z_nT)) * 1e-3f;
    lat = site_lat;
    lon = site_lon;
  }
  return f * g_err.field_scale;
}

// Inverse of bmi270*ToBodyNed: body NED (x, y, z) = (-raw.x, raw.y, -raw.z) * scale.
void bodyToRaw(const V3& body, float scale, int16_t raw[3]) {
  auto q = [](float v) { return int16_t(std::lround(std::max(-32767.0f, std::min(32767.0f, v)))); };
  raw[0] = q(-body.x() / scale);
  raw[1] = q(body.y() / scale);
  raw[2] = q(-body.z() / scale);
}

// Simulation trim: X/Y and Z compensation reduce to fixed scales.
constexpr float kMagXyUtPerLsb = 0.3125f, kMagZUtPerLsb = 0.32f;
constexpr uint16_t kRhall = 6000;
atoms3r_ical::Bmm150Trim simTrim() {
  atoms3r_ical::Bmm150Trim t;
  t.dig_xyz1 = kRhall;
  t.dig_z2 = 6400;
  t.valid = true;
  return t;
}

void encodeBmm(const V3& body, uint8_t m[8]) {
  // body = (s.y, -s.x, s.z) in sensor axes s.
  const float sx = -body.y(), sy = body.x(), sz = body.z();
  const int16_t x = int16_t(std::lround(sx / kMagXyUtPerLsb)), y = int16_t(std::lround(sy / kMagXyUtPerLsb));
  const int16_t z = int16_t(std::lround(sz / kMagZUtPerLsb));
  const uint16_t ux = uint16_t(x << 3), uy = uint16_t(y << 3), uz = uint16_t(z << 1);
  const uint16_t ur = uint16_t(kRhall << 2) | 1u;  // data ready
  const uint16_t w[4] = {ux, uy, uz, ur};
  for (int i = 0; i < 4; ++i) {
    m[2 * i] = uint8_t(w[i] & 0xff);
    m[2 * i + 1] = uint8_t(w[i] >> 8);
  }
}

// BMI270 with a BMM150 behind the AUX interface and a header-mode FIFO.
struct SimBmi270 : m5::IMU_Base {
  uint8_t reg[128]{};
  uint8_t bmm[128]{};
  std::deque<uint8_t> fifo;
  std::deque<size_t> frame_sizes;
  uint64_t start_us = 0, next_tick = 0;
  double t0_ticks = 53.0;
  uint64_t overflowed = 0, frames_written = 0;
  std::mt19937 rng{20261012};
  std::normal_distribution<float> n{0.0f, 1.0f};

  SimBmi270() {
    reg[0x00] = 0x24; reg[0x40] = 0xA8; reg[0x41] = 0x02; reg[0x42] = 0xA9; reg[0x44] = 0x46; reg[0x45] = 0x88;
    reg[0x49] = 0x10; reg[0x4B] = 0x20; reg[0x4C] = 0x4F; reg[0x4D] = 0x42; reg[0x4E] = 0x4C; reg[0x4F] = 0x38;
    reg[0x7D] = 0x0F;
    bmm[0x40] = 0x32; bmm[0x4B] = 0x01; bmm[0x4C] = 0x38; bmm[0x51] = 23; bmm[0x52] = 40;
    start_us = ins_sketch_host::clock_us();
  }
  bool manual() const { return (reg[0x4C] & 0x80) && !(reg[0x7D] & 0x01); }
  bool streaming() const {
    return reg[0x49] == 0xF0 && (reg[0x7D] & 0x07) == 0x07 && reg[0x40] == 0xA9 && (reg[0x42] | 0x40) == 0xE9;
  }
  double ticksNow() const { return double(ins_sketch_host::clock_us() - start_us) / kTickUs; }

  void push(const std::vector<uint8_t>& frame) {
    while (fifo.size() + frame.size() > 6144 && !frame_sizes.empty()) {
      for (size_t i = 0; i < frame_sizes.front(); ++i) fifo.pop_front();
      frame_sizes.pop_front();
      ++overflowed;
    }
    fifo.insert(fifo.end(), frame.begin(), frame.end());
    frame_sizes.push_back(frame.size());
  }

  void advance() {
    const double now = ticksNow();
    for (; t0_ticks + 128.0 * double(next_tick) <= now; ++next_tick) {
      if (!streaming()) continue;
      const double t_us = (t0_ticks + 128.0 * double(next_tick)) * kTickUs;
      const double t_s = t_us * 1e-6;
      const Q q = attitudeAt(t_s), q1 = attitudeAt(t_s + 0.001);
      const Eigen::AngleAxisf aa(q.conjugate() * q1);
      const V3 w_true = aa.axis() * (aa.angle() / 0.001f);
      const V3 f_true = q.conjugate() * V3(0, 0, -atoms3r_ical::ImuCalCfg::g_cal_local);
      const float d2r = atoms3r_ical::ImuCalCfg::DEG2RAD;
      const V3 w = g_err.gyro_scale * w_true + g_err.gyro_bias_dps * d2r + 0.07f * d2r * V3(n(rng), n(rng), n(rng));
      const V3 f = g_err.accel_scale * f_true + 0.0015f * 9.80665f * V3(n(rng), n(rng), n(rng));
      int16_t g[3], a[3];
      bodyToRaw(w, d2r * 2000.0f / 32768.0f, g);
      bodyToRaw(f, atoms3r_ical::ImuCalCfg::g_std * 8.0f / 32768.0f, a);
      const bool aux = (reg[0x7D] & 0x01) && (next_tick % 8 == 3);
      std::vector<uint8_t> frame{uint8_t(aux ? 0x9C : 0x8C)};
      if (aux) {
        uint8_t m[8];
        const V3 noise = 0.3f * V3(n(rng), n(rng), n(rng));
        encodeBmm(attitudeAt((t_us - kMagDelayUs) * 1e-6).conjugate() * fieldNed() + noise, m);
        frame.insert(frame.end(), m, m + 8);
      }
      for (const int16_t* v : {g, a})
        for (int i = 0; i < 3; ++i) {
          frame.push_back(uint8_t(uint16_t(v[i]) & 0xff));
          frame.push_back(uint8_t(uint16_t(v[i]) >> 8));
        }
      push(frame);
      ++frames_written;
    }
  }

  // I2C at 400 kHz: a transaction costs its addressing plus 9 bit times per
  // byte, and every register byte is sampled when it is clocked out, so a
  // burst latches sensortime before the FIFO length.
  double bus_frac_us = 0;
  void bus(double us) {
    bus_frac_us += us;
    const double whole = std::floor(bus_frac_us);
    ins_sketch_host::clock_us() += uint64_t(whole);
    bus_frac_us -= whole;
  }
  static constexpr double kByteUs = 22.5, kAddressUs = 60.0;

  bool readRegister(uint8_t r, uint8_t* b, size_t len) override {
    bus(kAddressUs);
    for (size_t i = 0; i < len; ++i) {
      advance();
      const uint8_t a = uint8_t(r == 0x26 ? 0x26 : r + i);
      if (a == 0x26) {
        if (fifo.empty()) {
          b[i] = 0x80;
        } else {
          b[i] = fifo.front();
          fifo.pop_front();
          if (!frame_sizes.empty() && --frame_sizes.front() == 0) frame_sizes.pop_front();
        }
      } else {
        if (a == 0x18) {
          const uint32_t st = uint32_t(uint64_t(ticksNow()) & 0xFFFFFF);
          reg[0x18] = uint8_t(st); reg[0x19] = uint8_t(st >> 8); reg[0x1A] = uint8_t(st >> 16);
        }
        if (a == 0x22) {
          const int16_t temp = int16_t((30.0f - 23.0f) * 512.0f);
          reg[0x22] = uint8_t(uint16_t(temp) & 0xff); reg[0x23] = uint8_t(uint16_t(temp) >> 8);
        }
        if (a == 0x24) { reg[0x24] = uint8_t(fifo.size() & 0xff); reg[0x25] = uint8_t((fifo.size() >> 8) & 0x3f); }
        b[i] = reg[a];
      }
      bus(kByteUs);
    }
    return true;
  }
  bool writeRegister8(uint8_t r, uint8_t v) override {
    bus(kAddressUs + kByteUs);
    advance();
    // A configuration change while streaming is marked in the FIFO.
    if (r == 0x42 && streaming() && reg[r] != v) push({0x48, 0x00});
    reg[r] = v;
    if (r == 0x4E && manual()) bmm[v] = reg[0x4F];
    if (r == 0x4D && manual()) reg[0x04] = bmm[v];
    if (r == 0x7E && v == 0xB0) { fifo.clear(); frame_sizes.clear(); }
    return true;
  }
};

std::string g_serial;
void sink(const char* d, size_t len) { g_serial.append(d, len); }

int failures = 0;
void check(bool ok, const char* what, double value, double expect) {
  std::printf("  %-58s %12.4f  (expect %.4f)  %s\n", what, value, expect, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

using Fields = std::map<std::string, std::string>;
Fields fields(const std::string& line) {
  Fields out;
  size_t start = 2;
  while (start < line.size()) {
    size_t end = line.find(',', start);
    if (end == std::string::npos) end = line.size();
    const std::string kv = line.substr(start, end - start);
    const size_t eq = kv.find('=');
    if (eq != std::string::npos) out[kv.substr(0, eq)] = kv.substr(eq + 1);
    start = end + 1;
  }
  return out;
}
double num(const Fields& f, const char* k) {
  const auto it = f.find(k);
  return it == f.end() ? NAN : std::atof(it->second.c_str());
}
std::string str(const Fields& f, const char* k) {
  const auto it = f.find(k);
  return it == f.end() ? std::string("missing") : it->second;
}

}  // namespace

int main(int argc, char** argv) {
  const bool bad = argc > 1 && !std::strcmp(argv[1], "bad");
  g_err = bad ? kBad : kGood;
  std::printf("FIFO diagnostics sketch on a simulated BMI270 (%s sensor): %s\n", bad ? "bad" : "good", FIFO_SKETCH_FILE);
  ins_sketch_host::clock_us() = 2000000;
  SimBmi270 sim;
  M5.Imu.device = &sim;
  M5DisplayStub::fillCostUs() = 6500;  // one 128x128 full-screen fill over SPI
  atoms3r_ical::atoms3rMagSource().trim = simTrim();
  HardwareSerial::capture() = &sink;

  setup();
  const uint64_t start_us = ins_sketch_host::clock_us();
  const uint64_t end_us = start_us + uint64_t(kDurationS * 1e6);
  bool tapped = false, moved = false, toggled = false;
  while (ins_sketch_host::clock_us() < end_us) {
    const double t_s = double(ins_sketch_host::clock_us() - start_us) * 1e-6;
    if (bad && !moved && t_s > 1.0) {
      HardwareSerial::input() = "L 64.1 -21.9 50\n";  // Reykjavik
      moved = true;
    }
    if (bad && !toggled && t_s > 2.0) {
      HardwareSerial::input() = "n\n";  // gyro noise-performance mode on
      toggled = true;
    }
    if (!tapped && t_s > 31.0) {
      atoms3r_ical::Input::tap() = true;  // next screen page
      tapped = true;
    }
    loop();
  }
  for (int i = 0; i < 400; ++i) loop();  // let the R-line queue drain
  HardwareSerial::capture() = nullptr;

  Fields last_t, last_a, f_line, last_f, h_line, g_line, last_d;
  std::vector<Fields> poses_out, rotations;
  std::vector<double> fits;
  long r_lines = 0, r_gaps = 0, expected_seq = 0;
  uint32_t prev_us = 0;
  double min_dt = 1e9, max_dt = 0;
  size_t pos = 0;
  while (pos < g_serial.size()) {
    size_t nl = g_serial.find('\n', pos);
    if (nl == std::string::npos) nl = g_serial.size();
    const std::string line = g_serial.substr(pos, nl - pos);
    pos = nl + 1;
    if (line.size() < 2 || line[1] != ',') continue;
    switch (line[0]) {
      case 'T': last_t = fields(line); break;
      case 'A': last_a = fields(line); break;
      case 'F': if (f_line.empty()) f_line = fields(line); last_f = fields(line); break;
      case 'D': last_d = fields(line); break;
      case 'H': if (h_line.empty()) h_line = fields(line); break;
      case 'G': g_line = fields(line); break;
      case 'S': poses_out.push_back(fields(line)); break;
      case 'W': rotations.push_back(fields(line)); break;
      case 'X': {
        const Fields x = fields(line);
        if (str(x, "ok") == "1") fits.push_back(num(x, "delay_ms"));
        break;
      }
      case 'R': {
        const long seq = std::strtol(line.c_str() + 2, nullptr, 10);
        const size_t c2 = line.find(',', 2);
        const uint32_t us = uint32_t(std::strtoul(line.c_str() + c2 + 1, nullptr, 10));
        if (seq != expected_seq) ++r_gaps;
        if (r_lines > 0) {
          const double dt = double(int32_t(us - prev_us));
          min_dt = std::min(min_dt, dt);
          max_dt = std::max(max_dt, dt);
        }
        expected_seq = seq + 1;
        prev_us = us;
        ++r_lines;
        break;
      }
      default: break;
    }
  }

  std::printf("timing\n");
  check(str(f_line, "ok") == "1", "FIFO setup succeeded (F line)", str(f_line, "ok") == "1", 1);
  check(str(h_line, "sketch") == "atomS3R_imu_fifo_diagnostics", "identity line (H)", !h_line.empty(), 1);
  check(sim.overflowed == 0, "simulated FIFO never overflowed", double(sim.overflowed), 0);
  check(str(last_t, "rate_ok") == "1", "rate_ok", num(last_t, "rate_hz"), 200.0 / (1.0 + kPpm * 1e-6));
  check(std::fabs(num(last_t, "rate_hz") - 200.0 / (1.0 + kPpm * 1e-6)) < 0.05, "host-clock rate is the sensor clock",
        num(last_t, "rate_hz"), 200.0 / (1.0 + kPpm * 1e-6));
  check(std::fabs(num(last_t, "sensor_rate_hz") - 200.0) < 0.2, "200 Hz on the sensor clock",
        num(last_t, "sensor_rate_hz"), 200.0);
  check(num(last_d, "lost") == 0 && num(last_d, "grid_breaks") == 0, "no frame counted lost, warm-up included",
        num(last_d, "lost"), 0);
  check(str(h_line, "gyr_T0_c") != "missing" && str(h_line, "mag_cal_field_uT") != "missing",
        "H line carries the saved calibration's context", str(h_line, "mag_cal_field_uT") != "missing", 1);
  if (bad)
    check(str(last_f, "gyr_conf_now") == "0xe9" && num(last_d, "config_frames") >= 1,
          "'n' switched the gyro to noise-performance mode", num(last_d, "config_frames"), 1);
  check(str(last_t, "no_loss_ok") == "1" && num(last_t, "lost") == 0, "no_loss_ok through display bursts",
        num(last_t, "lost"), 0);
  check(num(last_t, "stall_max_ms") >= 1500, "a 1.5 s display burst ran", num(last_t, "stall_max_ms"), 1500);
  check(str(last_t, "spacing_ok") == "1", "spacing_ok", num(last_t, "spacing_max_us"), 5000);
  check(std::fabs(num(last_t, "clock_ppm") - kPpm) < 30, "sensor clock error measured", num(last_t, "clock_ppm"), kPpm);
  check(str(last_t, "mag_ok") == "1", "mag_ok: 25 Hz, data ready, 40 ms apart", num(last_t, "aux_hz"), 25.0);
  check(str(last_t, "delay_ok") == "1", "delay_ok: three fits within 4 ms", num(last_t, "delay_spread_ms"), 4.0);
  // The fit uses the shipping rotation-history convention: a gyro sample is
  // the rate over the interval that ends at it. The simulated gyro sample at
  // t is the rate at t + 0.5 ms, so in that convention the field lags by the
  // simulated delay plus half a frame plus 0.5 ms.
  const double expected_delay_ms = kMagDelayUs * 1e-3 + 2.5 + 0.5;
  check(std::fabs(num(last_t, "delay_ms") - expected_delay_ms) <= 2.0, "measured delay matches the simulated delay",
        num(last_t, "delay_ms"), expected_delay_ms);
  check(r_gaps == 0 && r_lines > long(kDurationS * 190), "every frame has an R line", double(r_lines),
        kDurationS * 200);
  check(min_dt >= 4900 && max_dt <= 5100, "every R-line sample spacing within 4.9-5.1 ms", max_dt, 5000);
  std::printf("  R-line spacing %.0f..%.0f us, %zu delay fits\n", min_dt, max_dt, fits.size());

  std::printf("reference\n");
  const double want_lat = bad ? 64.1 : atoms3r_ical::CalibrationSite::latitude_deg;
  const imu_cal::GeomagField want =
      imu_cal::Wmm2025::field(want_lat, bad ? -21.9 : atoms3r_ical::CalibrationSite::longitude_deg,
                              bad ? 50.0 : atoms3r_ical::CalibrationSite::ellipsoidal_height_m, site_year);
  check(std::fabs(num(g_line, "lat") - want_lat) < 1e-6, bad ? "place moved over serial" : "place is the calibration site",
        num(g_line, "lat"), want_lat);
  check(std::fabs(num(g_line, "F_uT") - want.f_nT * 1e-3) < 0.002 && str(g_line, "model_valid") == "1",
        "expected field is WMM-2025 at that place", num(g_line, "F_uT"), want.f_nT * 1e-3);

  const double acc_mg = 1000.0 * (g_err.accel_scale - 1.0) * atoms3r_ical::ImuCalCfg::g_cal_local / 9.80665;
  const double bias_dps = g_err.gyro_bias_dps.norm(), scale_pct = 100.0 * (g_err.gyro_scale - 1.0);
  const double field_pct = 100.0 * (g_err.field_scale - 1.0);
  std::printf("still poses (%zu S lines)\n", poses_out.size());
  check(poses_out.size() == kPoses.size(), "one S line per resting pose", double(poses_out.size()),
        double(kPoses.size()));
  double worst_acc = 0, worst_bias = 0, worst_field = 0, worst_dip = 0;
  for (const Fields& p : poses_out) {
    worst_acc = std::max(worst_acc, std::fabs(num(p, "a_error_mg") - acc_mg));
    worst_bias = std::max(worst_bias, std::fabs(num(p, "gyro_bias_dps") - bias_dps));
    worst_field = std::max(worst_field, std::fabs(num(p, "F_error_pct") - field_pct));
    worst_dip = std::max(worst_dip, std::fabs(num(p, "dip_error_deg")));
  }
  check(worst_acc < 0.5, "|a| error measured in every pose (worst deviation, mg)", worst_acc, 0);
  check(worst_bias < 0.02, "gyro bias measured in every pose (worst deviation, deg/s)", worst_bias, 0);
  check(worst_field < 0.3, "field error measured in every pose (worst deviation, %)", worst_field, 0);
  check(worst_dip < 0.3, "dip matches WMM (isotropic field error, worst, deg)", worst_dip, 0);

  std::printf("rotations (%zu W lines)\n", rotations.size());
  double worst_scale = 0;
  int scale_fits = 0;
  for (const Fields& w : rotations) {
    const double s = num(w, "gyro_scale_pct");
    if (!std::isfinite(s)) continue;
    ++scale_fits;
    worst_scale = std::max(worst_scale, std::fabs(s - scale_pct));
  }
  check(rotations.size() == kPoses.size() - 1, "a rotation check between every pair of poses",
        double(rotations.size()), double(kPoses.size() - 1));
  check(scale_fits >= 5 && worst_scale < 0.25, "gyro scale error fitted on each (worst deviation, %)", worst_scale, 0);

  std::printf("accuracy criteria (A line)\n");
  const std::string want_verdict = bad ? "0" : "1";
  for (const char* k : {"accel_ok", "gyro_bias_ok", "gyro_scale_ok", "field_ok"})
    check(str(last_a, k) == want_verdict, k, str(last_a, k) == "1", bad ? 0 : 1);
  if (!bad) {
    check(str(last_a, "accel_noise_ok") == "1" && str(last_a, "gyro_noise_ok") == "1", "noise criteria pass",
          num(last_a, "a_noise_mg_max"), 1.5);
    check(str(last_a, "mag_norm_ok") == "1", "|m| constant in motion", num(last_a, "m_std_pct"), 0);
  }

  if (failures) {
    std::printf("FIFO diagnostics sketch (%s): %d check(s) FAILED\n", bad ? "bad" : "good", failures);
    return 1;
  }
  std::printf("FIFO diagnostics sketch (%s): PASS\n", bad ? "bad" : "good");
  return 0;
}
