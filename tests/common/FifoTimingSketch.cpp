// Host run of sensors/imu_basic/atomS3R_imu_fifo_timing against a simulated
// BMI270: a 6 kB header-mode FIFO filled at 200 Hz on a sensor clock that runs
// 0.2 % fast, and a forced-mode BMM150 whose readout every 40 ms carries the
// field of a fixed time before it. The board sways gently for 50 s, then
// rocks like a hand test (the motion the delay fit needs). The
// sketch's own setup()/loop() run unchanged, with its display bursts costing
// host time and one tap. Checks: the FIFO setup, the pass criteria the sketch
// prints, the measured magnetometer delay against the simulated one, and the
// R-line capture itself.
//
//   fifo_timing_sketch-test [--write-capture FILE]
//
// Build with -DFIFO_SKETCH_FILE='"/path/atomS3R_imu_fifo_timing.ino"' against
// tests/common/ins_sketch_stubs (see FifoTimingSketch.mk).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <Arduino.h>
#include <ArduinoEigenDense.h>

#include FIFO_SKETCH_FILE

namespace {

using V3 = Eigen::Vector3f;
constexpr double kPpm = 2000.0;
constexpr double kTickUs = 39.0625 * (1.0 + kPpm * 1e-6);
constexpr double kMagDelayUs = 26000.0;  // readout time minus the time the field belongs to
constexpr double kDurationS = 180.0;
constexpr double kCalmS = 50.0;
const V3 kFieldNed(20.0f, -1.5f, 45.0f);

Eigen::Quaternionf attitudeAt(double t_s) {
  const double u = std::min(1.0, std::max(0.0, (t_s - kCalmS) / 5.0));
  const float rock = float(u * u * (3 - 2 * u));
  const float roll = 0.05f * float(std::sin(2 * M_PI * 0.25 * t_s)) +
                     rock * 0.5f * float(std::sin(2 * M_PI * 1.0 * t_s));
  const float pitch = rock * 0.3f * float(std::sin(2 * M_PI * 0.7 * t_s + 1.0));
  const float yaw = 0.8f + 0.4f * float(std::sin(2 * M_PI * 0.05 * t_s));
  return (Eigen::AngleAxisf(yaw, V3::UnitZ()) * Eigen::AngleAxisf(pitch, V3::UnitY()) *
          Eigen::AngleAxisf(roll, V3::UnitX())).normalized();
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
    return reg[0x49] == 0xF0 && (reg[0x7D] & 0x07) == 0x07 && reg[0x40] == 0xA9 && reg[0x42] == 0xA9;
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
      const Eigen::Quaternionf q = attitudeAt(t_s), q1 = attitudeAt(t_s + 0.001);
      const Eigen::AngleAxisf aa(q.conjugate() * q1);
      const V3 w = aa.axis() * (aa.angle() / 0.001f);
      const V3 f = q.conjugate() * V3(0, 0, -atoms3r_ical::ImuCalCfg::g_cal_local);
      int16_t g[3], a[3];
      bodyToRaw(w, atoms3r_ical::ImuCalCfg::DEG2RAD * 2000.0f / 32768.0f, g);
      bodyToRaw(f, atoms3r_ical::ImuCalCfg::g_std * 8.0f / 32768.0f, a);
      const bool aux = (reg[0x7D] & 0x01) && (next_tick % 8 == 3);
      std::vector<uint8_t> frame{uint8_t(aux ? 0x9C : 0x8C)};
      if (aux) {
        uint8_t m[8];
        const V3 noise = 0.3f * V3(n(rng), n(rng), n(rng));
        encodeBmm(attitudeAt((t_us - kMagDelayUs) * 1e-6).conjugate() * kFieldNed + noise, m);
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

  bool readRegister(uint8_t r, uint8_t* b, size_t len) override {
    advance();
    if (r == 0x26) {
      for (size_t i = 0; i < len; ++i) {
        if (fifo.empty()) { b[i] = 0x80; continue; }
        b[i] = fifo.front();
        fifo.pop_front();
        if (!frame_sizes.empty() && --frame_sizes.front() == 0) frame_sizes.pop_front();
      }
      return true;
    }
    const uint32_t st = uint32_t(uint64_t(ticksNow()) & 0xFFFFFF);
    reg[0x18] = uint8_t(st); reg[0x19] = uint8_t(st >> 8); reg[0x1A] = uint8_t(st >> 16);
    const int16_t temp = int16_t((30.0f - 23.0f) * 512.0f);
    reg[0x22] = uint8_t(uint16_t(temp) & 0xff); reg[0x23] = uint8_t(uint16_t(temp) >> 8);
    reg[0x24] = uint8_t(fifo.size() & 0xff); reg[0x25] = uint8_t((fifo.size() >> 8) & 0x3f);
    for (size_t i = 0; i < len; ++i) b[i] = reg[r + i];
    return true;
  }
  bool writeRegister8(uint8_t r, uint8_t v) override {
    advance();
    reg[r] = v;
    if (r == 0x4E && manual()) bmm[v] = reg[0x4F];
    if (r == 0x4D && manual()) reg[0x04] = bmm[v];
    if (r == 0x7E && v == 0xB0) { fifo.clear(); frame_sizes.clear(); }
    return true;
  }
};

std::string g_serial;
void sink(const char* d, size_t n) { g_serial.append(d, n); }

int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-52s %12.4f  (limit %.4f)  %s\n", what, value, limit, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

std::map<std::string, std::string> fields(const std::string& line) {
  std::map<std::string, std::string> out;
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
double num(const std::map<std::string, std::string>& f, const char* k) {
  const auto it = f.find(k);
  return it == f.end() ? NAN : std::atof(it->second.c_str());
}

}  // namespace

int main(int argc, char** argv) {
  const char* capture_path = nullptr;
  for (int i = 1; i < argc; ++i)
    if (!std::strcmp(argv[i], "--write-capture") && i + 1 < argc) capture_path = argv[++i];

  std::printf("FIFO timing sketch on a simulated BMI270: %s\n", FIFO_SKETCH_FILE);
  ins_sketch_host::clock_us() = 2000000;
  SimBmi270 sim;
  M5.Imu.device = &sim;
  M5DisplayStub::fillCostUs() = 6500;  // one 128x128 full-screen fill over SPI
  atoms3r_ical::atoms3rMagSource().trim = simTrim();
  HardwareSerial::capture() = &sink;

  setup();
  const uint64_t end_us = ins_sketch_host::clock_us() + uint64_t(kDurationS * 1e6);
  bool tapped = false;
  while (ins_sketch_host::clock_us() < end_us) {
    if (!tapped && ins_sketch_host::clock_us() > end_us - uint64_t(kDurationS * 1e6) + 31000000u) {
      atoms3r_ical::Input::tap() = true;  // a 1.5 s display burst
      tapped = true;
    }
    loop();
  }
  // Let the R-line queue drain.
  for (int i = 0; i < 400; ++i) loop();
  HardwareSerial::capture() = nullptr;
  if (capture_path) std::ofstream(capture_path, std::ios::binary) << g_serial;

  std::map<std::string, std::string> last_p, f_line;
  std::vector<double> fits;
  long r_lines = 0, r_gaps = 0;
  long expected_seq = 0;
  uint32_t prev_us = 0;
  double min_dt = 1e9, max_dt = 0;
  size_t pos = 0;
  while (pos < g_serial.size()) {
    size_t nl = g_serial.find('\n', pos);
    if (nl == std::string::npos) nl = g_serial.size();
    const std::string line = g_serial.substr(pos, nl - pos);
    pos = nl + 1;
    if (line.rfind("P,", 0) == 0) last_p = fields(line);
    else if (line.rfind("F,", 0) == 0 && f_line.empty()) f_line = fields(line);
    else if (line.rfind("X,", 0) == 0) {
      const auto x = fields(line);
      if (x.at("ok") == "1") fits.push_back(num(x, "delay_ms"));
    } else if (line.rfind("R,", 0) == 0) {
      const long seq = std::strtol(line.c_str() + 2, nullptr, 16);
      const size_t c2 = line.find(',', 2);
      const uint32_t us = uint32_t(std::strtoul(line.c_str() + c2 + 1, nullptr, 16));
      if (seq != expected_seq) ++r_gaps;
      if (r_lines > 0) {
        const double dt = double(int32_t(us - prev_us));
        min_dt = std::min(min_dt, dt);
        max_dt = std::max(max_dt, dt);
      }
      expected_seq = seq + 1;
      prev_us = us;
      ++r_lines;
    }
  }

  check(f_line["ok"] == "1", "FIFO setup succeeded (F line)", f_line["ok"] == "1", 1);
  check(sim.overflowed == 0, "simulated FIFO never overflowed", double(sim.overflowed), 0);
  check(!last_p.empty(), "pass criteria printed (P line)", last_p.empty() ? 0 : 1, 1);
  check(last_p["rate_ok"] == "1", "rate_ok", num(last_p, "rate_hz"), 200.0);
  check(std::fabs(num(last_p, "rate_hz") - 200.0 * 1.0 / (1.0 + kPpm * 1e-6)) < 0.05,
        "host-clock rate matches the sensor clock", num(last_p, "rate_hz"), 200.0 / (1.0 + kPpm * 1e-6));
  check(last_p["no_loss_ok"] == "1" && num(last_p, "lost") == 0, "no_loss_ok through display bursts",
        num(last_p, "lost"), 0);
  check(num(last_p, "stall_max_ms") >= 1500, "a 1.5 s display burst ran", num(last_p, "stall_max_ms"), 1500);
  check(last_p["spacing_ok"] == "1", "spacing_ok (sketch)", num(last_p, "spacing_max_us"), 5100);
  check(std::fabs(num(last_p, "clock_ppm") - kPpm) < 30, "sensor clock error measured", num(last_p, "clock_ppm"), kPpm);
  check(last_p["mag_ok"] == "1", "mag_ok: 25 Hz, data ready, 40 ms apart", num(last_p, "aux_hz"), 25.0);
  check(last_p["delay_ok"] == "1", "delay_ok: three fits within 4 ms", num(last_p, "delay_spread_ms"), 4.0);
  // The fit uses the shipping rotation-history convention: a gyro sample is
  // the rate over the interval that ends at it. The simulated gyro sample at
  // t is the rate at t + 0.5 ms, so in that convention the field lags by the
  // simulated delay plus half a frame plus 0.5 ms.
  const double expected_delay_ms = kMagDelayUs * 1e-3 + 2.5 + 0.5;
  check(std::fabs(num(last_p, "delay_ms") - expected_delay_ms) <= 2.0, "measured delay matches the simulated delay",
        num(last_p, "delay_ms"), expected_delay_ms);
  check(last_p["all_ok"] == "1", "all_ok", last_p["all_ok"] == "1", 1);
  check(r_gaps == 0 && r_lines > long(kDurationS * 190), "every frame has an R line", double(r_lines),
        kDurationS * 200);
  check(min_dt >= 4900 && max_dt <= 5100, "every R-line sample spacing within 4.9-5.1 ms", max_dt, 5100);
  std::printf("  R-line spacing %.0f..%.0f us, %zu delay fits\n", min_dt, max_dt, fits.size());
  if (failures) {
    std::printf("FIFO timing sketch: %d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("FIFO timing sketch: PASS\n");
  return 0;
}
