// BMI270 FIFO acquisition (src/AtomS3R/AtomS3R_Bmi270Fifo.h): register setup on
// a simulated BMI270 with a BMM150 behind its AUX interface, header-mode FIFO
// parsing, unit equivalence with the M5Unified path, and the sensortime clock
// map under stalls, clock error, lost frames and counter wrap.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include "AtomS3R/AtomS3R_Bmi270Fifo.h"

namespace {

int failures = 0;
void require(bool ok, const char* what) {
  std::printf("  %-74s %s\n", what, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

using atoms3r_ical::Bmi270;

// BMI270 after M5Unified's begin(): AUX data mode with FCU write, BMM150 in
// normal mode at 30 Hz.
struct FakeBmi270 {
  uint8_t reg[128]{};
  uint8_t bmm[128]{};
  std::vector<std::pair<uint8_t, uint8_t>> writes;
  int aux_writes_in_data_mode = 0, forced_commands = 0;
  uint8_t stuck_reg = 0xFF;  // a register that ignores writes
  FakeBmi270() {
    reg[0x00] = 0x24; reg[0x40] = 0xA8; reg[0x41] = 0x02; reg[0x42] = 0xA9; reg[0x43] = 0x00;
    reg[0x44] = 0x46; reg[0x45] = 0x88; reg[0x49] = 0x10; reg[0x4B] = 0x20; reg[0x4C] = 0x4F;
    reg[0x4D] = 0x42; reg[0x4E] = 0x4C; reg[0x4F] = 0x38; reg[0x7D] = 0x0F;
    bmm[0x40] = 0x32; bmm[0x4B] = 0x01; bmm[0x4C] = 0x38; bmm[0x51] = 23; bmm[0x52] = 40;
  }
  bool manual() const { return (reg[0x4C] & 0x80) && !(reg[0x7D] & 0x01); }
  bool readRegister(uint8_t r, uint8_t* b, size_t n) {
    for (size_t i = 0; i < n; ++i) b[i] = reg[r + i];
    return true;
  }
  bool writeRegister8(uint8_t r, uint8_t v) {
    writes.emplace_back(r, v);
    if (r == stuck_reg) return true;
    reg[r] = v;
    if (r == 0x4E) {
      if (!manual()) ++aux_writes_in_data_mode;
      else {
        bmm[v] = reg[0x4F];
        if (v == 0x4C && (reg[0x4F] & 0x06) == 0x02) ++forced_commands;
      }
    }
    if (r == 0x4D && manual()) reg[0x04] = bmm[v];
    return true;
  }
  int indexOf(uint8_t r, uint8_t v) const {
    for (size_t i = 0; i < writes.size(); ++i)
      if (writes[i].first == r && writes[i].second == v) return int(i);
    return -1;
  }
};
void noWait(uint32_t) {}

void testSetup() {
  std::printf("setup\n");
  FakeBmi270 dev;
  atoms3r_ical::Bmi270FifoSetup<FakeBmi270> setup(&dev, noWait);
  atoms3r_ical::FifoAcquisitionConfig cfg;
  const bool ok = setup.apply(cfg);
  std::printf("  failure: %s\n", setup.failure());
  require(ok, "setup applies and reads back");
  require(dev.reg[0x40] == 0xA9 && dev.reg[0x42] == 0xA9, "accelerometer and gyro at 200 Hz");
  require(dev.reg[0x41] == 0x02 && dev.reg[0x43] == 0x00, "ranges +-8 g and +-2000 deg/s");
  require(dev.reg[0x44] == 0x06, "AUX at 25 Hz, longest trigger-to-readout offset");
  require(dev.reg[0x48] == 0x00 && dev.reg[0x49] == 0xF0 && dev.reg[0x45] == 0x88,
          "FIFO: header mode, accelerometer+gyro+AUX, filtered, no time frame");
  require(dev.reg[0x4C] == 0x4F && dev.reg[0x4D] == 0x42, "AUX data mode with FCU write, burst from data X");
  require(dev.reg[0x4E] == 0x4C && dev.reg[0x4F] == 0x3A, "FCU repeats the forced-mode command (ODR bits kept)");
  require(dev.bmm[0x4C] == 0x3A && dev.forced_commands == 1, "BMM150 left in forced mode, first conversion started");
  require(dev.bmm[0x51] == 23 && dev.bmm[0x52] == 40, "low-noise repetitions untouched");
  require(dev.reg[0x7D] == 0x0F, "accelerometer, gyro, AUX and temperature powered");
  require(dev.aux_writes_in_data_mode == 0, "BMM150 written only in manual mode");
  const int if_conf = dev.indexOf(0x4C, 0x4F), aux_on = dev.indexOf(0x7D, 0x0F), flush = dev.indexOf(0x7E, 0xB0);
  require(if_conf >= 0 && aux_on > if_conf && flush > aux_on, "data mode before AUX power-up, FIFO flushed last");
  require(setup.readback().bmm_mode_before == 0x38 && setup.readback().aux_wr_data == 0x3A,
          "read-back records the previous and the new BMM150 mode");

  FakeBmi270 other;
  other.reg[0x00] = 0x26;
  atoms3r_ical::Bmi270FifoSetup<FakeBmi270> s2(&other, noWait);
  require(!s2.apply(cfg) && !std::strcmp(s2.failure(), "not a BMI270"), "another chip is refused untouched");
  require(other.writes.empty(), "nothing written to another chip");

  FakeBmi270 no_bmm;
  no_bmm.bmm[0x40] = 0x00;
  atoms3r_ical::Bmi270FifoSetup<FakeBmi270> s3(&no_bmm, noWait);
  require(!s3.apply(cfg) && !std::strcmp(s3.failure(), "BMM150 chip id mismatch"), "missing BMM150 is reported");

  FakeBmi270 stuck;
  stuck.stuck_reg = 0x40;
  atoms3r_ical::Bmi270FifoSetup<FakeBmi270> s4(&stuck, noWait);
  require(!s4.apply(cfg) && !std::strcmp(s4.failure(), "read-back differs from configuration"),
          "a register that does not take the value fails the read-back");

  atoms3r_ical::Bmi270FifoSetup<FakeBmi270> s5(nullptr, noWait);
  require(!s5.apply(cfg) && !std::strcmp(s5.failure(), "no BMI270"), "no device");
}

void put16(std::vector<uint8_t>& b, int16_t v) {
  b.push_back(uint8_t(uint16_t(v) & 0xff));
  b.push_back(uint8_t(uint16_t(v) >> 8));
}

void testParser() {
  std::printf("FIFO parser\n");
  std::vector<uint8_t> s;
  // gyro+accelerometer
  s.push_back(0x8C);
  put16(s, 10); put16(s, -20); put16(s, 30);           // gyro
  put16(s, -1000); put16(s, 2000); put16(s, 4096);     // accelerometer
  // AUX+gyro+accelerometer
  s.push_back(0x9C);
  for (int i = 0; i < 8; ++i) s.push_back(uint8_t(0x10 + i));
  put16(s, 1); put16(s, 2); put16(s, 3);
  put16(s, 4); put16(s, 5); put16(s, 6);
  s.push_back(0x40); s.push_back(3);                    // skip frame: 3 frames dropped
  s.push_back(0x90);                                    // AUX only
  for (int i = 0; i < 8; ++i) s.push_back(uint8_t(0x20 + i));
  s.push_back(0x44); s.push_back(1); s.push_back(2); s.push_back(3);  // sensortime
  s.push_back(0x48); s.push_back(0);                    // configuration change
  s.push_back(0x8C);
  put16(s, 7); put16(s, 8); put16(s, 9);
  put16(s, 10); put16(s, 11); put16(s, 12);
  const size_t complete = s.size();
  s.push_back(0x80);  // over-read

  atoms3r_ical::FifoParseStats st;
  std::vector<atoms3r_ical::FifoFrame> frames;
  unsigned skipped = 0;
  auto on_frame = [&](const atoms3r_ical::FifoFrame& f) { frames.push_back(f); };
  auto on_skip = [&](unsigned n) { skipped += n; };
  const size_t used = atoms3r_ical::parseBmi270Fifo(s.data(), s.size(), st, on_frame, on_skip);
  require(used == complete, "parsing stops at the over-read marker");
  require(frames.size() == 4 && st.acc_gyr == 3 && st.aux_only == 1 && st.aux == 2, "frame kinds counted");
  require(frames[0].g[1] == -20 && frames[0].a[0] == -1000 && frames[0].a[2] == 4096,
          "gyro before accelerometer, little-endian");
  require(frames[1].aux && frames[1].m[0] == 0x10 && frames[1].m[7] == 0x17 && frames[1].g[0] == 1 &&
          frames[1].a[2] == 6, "AUX payload precedes gyro and accelerometer");
  require(frames[2].aux && !frames[2].acc && frames[2].m[0] == 0x20, "AUX-only frame");
  require(skipped == 3 && st.skip_frames == 1 && st.sensortime_frames == 1 && st.config_frames == 1,
          "skip, sensortime and configuration frames");

  // The same stream delivered in pieces: only whole frames are consumed.
  std::vector<atoms3r_ical::FifoFrame> pieces;
  atoms3r_ical::FifoParseStats st2;
  std::vector<uint8_t> pending;
  for (size_t i = 0; i < complete; i += 5) {
    pending.insert(pending.end(), s.begin() + long(i), s.begin() + long(std::min(complete, i + 5)));
    const size_t n = atoms3r_ical::parseBmi270Fifo(
        pending.data(), pending.size(), st2, [&](const atoms3r_ical::FifoFrame& f) { pieces.push_back(f); },
        [](unsigned) {});
    pending.erase(pending.begin(), pending.begin() + long(n));
  }
  require(pieces.size() == 4 && pending.empty() && pieces[3].a[2] == 12, "chunked delivery parses identically");

  const uint8_t junk[] = {0x8C, 1, 2, 0x13};
  atoms3r_ical::FifoParseStats st3;
  require(atoms3r_ical::parseBmi270Fifo(junk, 3, st3, [](const atoms3r_ical::FifoFrame&) {}, [](unsigned) {}) == 0,
          "an incomplete frame is not consumed");
  const uint8_t bad[] = {0x13, 0, 0};
  require(atoms3r_ical::parseBmi270Fifo(bad, 3, st3, [](const atoms3r_ical::FifoFrame&) {}, [](unsigned) {}) == 0 &&
              st3.unknown == 1, "an unknown header stops parsing and is counted");
}

// M5Unified IMU_Class::getImuData for AtomS3R (axis order yxz, y inverted;
// offsets cleared) followed by the project's map to body NED.
float m5Axis(int16_t raw, bool invert, float res) {
  int32_t value = int32_t(raw) << 16;
  if (invert) value = -value;
  return res * (1.0f / 65536.0f) * float(value);
}

void testUnits() {
  std::printf("units\n");
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> u(-32767, 32767);
  int mismatches = 0;
  for (int i = 0; i < 20000; ++i) {
    const int16_t r[3] = {int16_t(u(rng)), int16_t(u(rng)), int16_t(u(rng))};
    const float ares = 8.0f / 32768.0f, gres = 2000.0f / 32768.0f;
    const Eigen::Vector3f a_m5 = atoms3r_ical::accel_nominal_g_to_body_ned_si_(
        m5Axis(r[1], false, ares), m5Axis(r[0], true, ares), m5Axis(r[2], false, ares));
    const Eigen::Vector3f g_m5 = atoms3r_ical::map_sensor_xyz_to_body_ned_(
        m5Axis(r[1], false, gres), m5Axis(r[0], true, gres), m5Axis(r[2], false, gres),
        atoms3r_ical::ImuCalCfg::DEG2RAD);
    const Eigen::Vector3f a = atoms3r_ical::bmi270AccelToBodyNed(r, atoms3r_ical::bmi270AccelGPerLsb(0x02));
    const Eigen::Vector3f g = atoms3r_ical::bmi270GyroToBodyNed(r, atoms3r_ical::bmi270GyroDpsPerLsb(0x00));
    if (std::memcmp(a.data(), a_m5.data(), sizeof(float) * 3) || std::memcmp(g.data(), g_m5.data(), sizeof(float) * 3))
      ++mismatches;
  }
  require(mismatches == 0, "FIFO samples equal readImuMapped() units bit for bit");
  require(atoms3r_ical::bmi270AccelGPerLsb(0) == 2.0f / 32768.0f &&
              atoms3r_ical::bmi270AccelGPerLsb(3) == 16.0f / 32768.0f &&
              atoms3r_ical::bmi270GyroDpsPerLsb(4) == 125.0f / 32768.0f, "range encodings");
  const uint8_t a[8] = {0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x01, 0x80};
  uint8_t b[8];
  std::memcpy(b, a, 8);
  b[6] = 0x00;  // data-ready cleared, same reading
  require(atoms3r_ical::bmm150DataReady(a) && !atoms3r_ical::bmm150DataReady(b), "BMM150 data-ready bit");
  require(atoms3r_ical::bmm150SamePayload(a, b), "status bits ignored when comparing readings");
  b[1] ^= 1;
  require(!atoms3r_ical::bmm150SamePayload(a, b), "a changed reading differs");
}

// Sensor clock: sensortime ticks run at tick_us host microseconds; frames are
// written every 128 ticks from phase t0; the host reads at irregular times
// (with a 1.5 s stall every 20 s) and receives the frames written up to then.
// drop_count frames written after drop_at_frame never reach the host and no
// skip frame reports them.
struct ClockRun {
  double spacing_us = 0, ppm = 0, min_spacing_us = 1e9, max_spacing_us = 0, max_abs_err_us = 0;
  uint32_t breaks = 0;
  uint64_t lost = 0;
};

ClockRun runClock(double ppm, double seconds, uint32_t host_start, uint32_t st_start, uint64_t drop_at_frame,
                  uint64_t drop_count, unsigned seed) {
  const double tick_us = Bmi270::TICK_US * (1.0 + ppm * 1e-6), t0_ticks = 37.0;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> jitter(0.0, 400.0);
  atoms3r_ical::SensorClockMap map;
  ClockRun r;
  uint64_t written = 0, received = 0;
  uint32_t prev_us = 0;
  bool have_prev = false;
  for (double t_us = 0; t_us < seconds * 1e6;) {
    t_us += 5000.0 + jitter(rng);
    if (std::fmod(t_us, 20e6) < 6000.0) t_us += 1.5e6;
    const double ticks = t_us / tick_us;
    while (t0_ticks + 128.0 * double(written) <= ticks) ++written;
    const uint64_t lost_now = written > drop_at_frame ? std::min(drop_count, written - drop_at_frame) : 0;
    const uint64_t counted = written - lost_now;
    map.observe(uint32_t(uint64_t(host_start) + uint64_t(t_us)),
                uint32_t((uint64_t(st_start) + uint64_t(ticks)) & Bmi270::SENSORTIME_MASK), counted);
    for (; received < counted; ++received) {
      const uint32_t us = map.frameHostUs(received);
      // The sensor index of this frame, and its true host time.
      const uint64_t index = received < drop_at_frame ? received : received + drop_count;
      const double truth = double(host_start) + (t0_ticks + 128.0 * double(index)) * tick_us;
      if (t_us > 5e6) {
        r.max_abs_err_us = std::max(r.max_abs_err_us, std::fabs(std::remainder(double(us) - truth, 4294967296.0)));
        if (have_prev) {
          const double sp = double(int32_t(us - prev_us));
          r.min_spacing_us = std::min(r.min_spacing_us, sp);
          r.max_spacing_us = std::max(r.max_spacing_us, sp);
        }
      }
      prev_us = us;
      have_prev = true;
    }
  }
  r.breaks = map.gridBreaks();
  r.lost = map.lostFrames();
  r.spacing_us = map.frameSpacingUs();
  r.ppm = map.clockPpm();
  return r;
}

void testClock() {
  std::printf("sensor clock\n");
  const uint64_t never = ~uint64_t(0);
  // Sensor clock 0.3 % slow; host micros() and the 24-bit sensortime both wrap.
  const ClockRun a = runClock(3000.0, 900.0, 0xFFFFFFFFu - 3000000u, 0xFFFF00u, never, 0, 1);
  std::printf("  spacing %.3f us, clock %+.0f ppm, delivered spacing %.1f..%.1f us, max time error %.0f us\n",
              a.spacing_us, a.ppm, a.min_spacing_us, a.max_spacing_us, a.max_abs_err_us);
  require(std::fabs(a.spacing_us - 128.0 * Bmi270::TICK_US * 1.003) < 0.2, "frame spacing from the sensor clock");
  require(std::fabs(a.ppm - 3000.0) < 50.0, "sensor clock error estimated");
  require(a.min_spacing_us > 5013.0 && a.max_spacing_us < 5017.0,
          "delivered spacing constant through 1.5 s stalls and both wraps");
  require(a.breaks == 0 && a.lost == 0, "no frames reported lost when none are");
  require(a.max_abs_err_us < 1000.0, "frame times anchored to the host clock within 1 ms");

  const ClockRun b = runClock(-1500.0, 120.0, 1000u, 0u, 4000, 7, 2);
  std::printf("  with 7 unreported lost frames: lost %llu, breaks %u, max time error %.0f us\n",
              static_cast<unsigned long long>(b.lost), b.breaks, b.max_abs_err_us);
  require(b.lost == 7 && b.breaks == 0, "frames lost without a skip frame are counted");
  require(b.max_abs_err_us < 1000.0, "frames after the loss keep their place on the sensor clock");
}

}  // namespace

int main() {
  testSetup();
  testParser();
  testUnits();
  testClock();
  if (failures) {
    std::printf("bmi270_fifo: %d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("bmi270_fifo: PASS\n");
  return 0;
}
