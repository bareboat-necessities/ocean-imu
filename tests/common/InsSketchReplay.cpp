// Host replay of an AtomS3R INS diagnostic capture through the unchanged
// sketch, and the round-trip test that proves the replay is exact.
//
// The sketch (.ino) is compiled with SEA_STATE_SERIAL_DIAG=1 against the host
// stand-ins in ins_sketch_stubs/, exactly as InsSketchHeaveDrift.cpp compiles
// it with the diagnostic stream off. Build with -DINS_SKETCH_FILE and
// -DINS_READY(app) (see InsSketchHeaveDrift.mk).
//
//   ins_sketch_replay-test [--write-capture FILE]
//                                          round-trip self-test (CI); optionally
//                                          keeps the synthetic capture
//   ins_sketch_replay-test CAPTURE [opts]  replay a device capture; prints the
//                                          replayed sketch's own diagnostic
//                                          lines (D/M/H/I/C, no R) to stdout
//     --segment N    replay the N-th boot in the capture (default 0)
//     --allow-gaps   continue across missing R lines (no longer exact)
//     --expect-ready fail unless the fused output is ready at the end
//
// A capture must come from the same sketch, or from the acquisition-only
// FIFO timing sketch (atomS3R_imu_fifo_timing, H line "acquisition=fifo"):
// its frames are what this sketch would receive from the FIFO acquisition,
// so any INS sketch can be replayed on them.
//
// What is exact: every frame's sensor values, magnetometer acquisition record
// and host sample time are the device's, bit for bit, and the calibration blob
// and magnetometer source are the ones the device loaded. What is not: host
// time inside a frame equals the frame time (the device spends some hundred
// microseconds there), and polls that delivered no frame are not replayed.
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <Arduino.h>
#include <ArduinoEigenDense.h>

#if !SEA_STATE_SERIAL_DIAG
#error "build with -DSEA_STATE_SERIAL_DIAG=1"
#endif

// Test access to the sketch's private state (outputs only; nothing is written
// except the calibration store, before begin(), as the device's NVS would be).
#define private public
#include INS_SKETCH_FILE
#undef private

namespace {

using atoms3r_ical::ImuSample;
using atoms3r_ical::floatFromBits;
using V3 = Eigen::Vector3f;

// ---------------------------------------------------------------- capture --

struct Frame {
  uint32_t seq = 0;
  ImuSample s{};
};

struct Calibration {  // one C line
  bool present = false, have_blob = false;
  unsigned mag_source = 0;
  std::vector<uint8_t> blob;
};

struct Segment {  // one boot: R sequence starting at 0 (or the capture start)
  Calibration cal;  // the C line in force for this boot
  std::vector<Frame> frames;
  std::vector<uint32_t> reset_events;  // frame seq at which reset_fusion was logged
};

struct Capture {
  std::string sketch, acquisition;
  std::vector<Segment> segments;
  long lines = 0, bad_lines = 0, d_lines = 0;
};

std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : line) {
    if (c == ',') { out.push_back(cur); cur.clear(); }
    else if (c != '\r' && c != '\n') cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}

bool hexU32(const std::string& t, uint32_t& v) {
  if (t.empty() || t.size() > 8) return false;
  char* end = nullptr;
  const unsigned long x = std::strtoul(t.c_str(), &end, 16);
  if (!end || *end) return false;
  v = static_cast<uint32_t>(x);
  return true;
}

std::string keyValue(const std::vector<std::string>& tok, const char* key) {
  const std::string k = std::string(key) + "=";
  for (const auto& t : tok)
    if (t.rfind(k, 0) == 0) return t.substr(k.size());
  return {};
}

bool parseR(const std::vector<std::string>& tok, Frame& f) {
  if (tok.size() != 20) return false;
  uint32_t u[19];
  for (int i = 0; i < 19; ++i)
    if (!hexU32(tok[size_t(i) + 1], u[i])) return false;
  f.seq = u[0];
  f.s.sample_us = u[1];
  f.s.mask = u[2];
  f.s.a = V3(floatFromBits(u[3]), floatFromBits(u[4]), floatFromBits(u[5]));
  f.s.w = V3(floatFromBits(u[6]), floatFromBits(u[7]), floatFromBits(u[8]));
  f.s.m = V3(floatFromBits(u[9]), floatFromBits(u[10]), floatFromBits(u[11]));
  f.s.tempC = floatFromBits(u[12]);
  auto& t = f.s.mag_timing;
  t.have = u[13] != 0;
  t.last_ms = u[14];
  t.sequence = u[15];
  t.frame_us = u[16];
  t.read_start_us = u[17];
  t.read_end_us = u[18];
  return true;
}

Capture parseCapture(std::istream& in) {
  Capture c;
  Calibration pending;  // latest C line not yet attached to a boot
  std::string line;
  while (std::getline(in, line)) {
    ++c.lines;
    if (line.size() < 2 || line[1] != ',') continue;  // boot prints and other text
    const auto tok = split(line);
    switch (line[0]) {
      case 'R': {
        Frame f;
        if (!parseR(tok, f)) { ++c.bad_lines; break; }
        if (c.segments.empty() || f.seq == 0) {
          c.segments.emplace_back();
          c.segments.back().cal = pending;
        }
        c.segments.back().frames.push_back(f);
        break;
      }
      case 'C': {
        if (tok.size() != 4 || tok[3].size() != 2 * sizeof(atoms3r_ical::ImuCalBlobV5)) {
          ++c.bad_lines;
          break;
        }
        // The sketch prints C right after its first frame and then
        // periodically: a boot keeps the first one it sees.
        Calibration cal;
        cal.present = true;
        cal.have_blob = tok[1] == "1";
        cal.mag_source = unsigned(std::strtoul(tok[2].c_str(), nullptr, 10));
        for (size_t i = 0; i < tok[3].size(); i += 2)
          cal.blob.push_back(uint8_t(std::strtoul(tok[3].substr(i, 2).c_str(), nullptr, 16)));
        pending = cal;
        if (!c.segments.empty() && !c.segments.back().cal.present) c.segments.back().cal = cal;
        break;
      }
      case 'H':
        if (c.sketch.empty()) c.sketch = keyValue(tok, "sketch");
        if (c.acquisition.empty()) c.acquisition = keyValue(tok, "acquisition");
        break;
      case 'E':
        if (tok.size() >= 2 && tok[1] == "reset_fusion" && !c.segments.empty()) {
          c.segments.back().reset_events.push_back(
              uint32_t(std::strtoul(keyValue(tok, "seq").c_str(), nullptr, 10)));
        }
        break;
      case 'D':
        ++c.d_lines;
        break;
      default:
        break;
    }
  }
  return c;
}

// ----------------------------------------------------------------- replay --

// Per-frame published outputs, compared bit for bit by the self-test.
struct Outputs {
  float roll, pitch, heading, heave_raw, heave_clean, dt;
  bool heading_valid, ready;
  bool operator==(const Outputs& o) const {
    auto same = [](float a, float b) { return std::memcmp(&a, &b, sizeof(float)) == 0; };
    return same(roll, o.roll) && same(pitch, o.pitch) && same(heading, o.heading) &&
           same(heave_raw, o.heave_raw) && same(heave_clean, o.heave_clean) && same(dt, o.dt) &&
           heading_valid == o.heading_valid && ready == o.ready;
  }
};

Outputs outputsOf(const FusionApp& app) {
  return {app.roll_deg_, app.pitch_deg_, app.heading_deg_, app.heave_raw_m_,
          app.heave_wave_clean_m_, app.dt_, app.heading_valid_, bool(INS_READY(app))};
}

std::string* g_sink = nullptr;
void sinkWrite(const char* data, size_t n) {
  if (g_sink) g_sink->append(data, n);
}

const Frame* g_frame = nullptr;
bool replaySource(uint32_t, ImuSample& out) {
  if (!g_frame) return false;
  const uint32_t sample_us = out.sample_us, mask = out.mask;
  out = g_frame->s;
  out.sample_us = sample_us;  // the sketch's own micros(), set to the frame's time below
  out.mask = mask;
  return true;
}

// Replays one boot segment from a fresh sketch. outputs (optional) receives
// each frame's published state; serial (optional) the sketch's own output.
bool replay(const Segment& seg, bool allow_gaps, std::vector<Outputs>* outputs, std::string* serial) {
  if (seg.frames.empty()) {
    std::fprintf(stderr, "replay: no R frames\n");
    return false;
  }
  // 64-bit host time across micros() wraps.
  uint64_t t64 = seg.frames.front().s.sample_us;
  const uint64_t start = t64 > 3000000u ? t64 - 3000000u : 0u;
  ins_sketch_host::clock_us() = start;

  atoms3r_ical::activeMagSource() = uint8_t(seg.cal.mag_source);
  auto app = std::make_unique<FusionApp>();
  if (seg.cal.present && seg.cal.have_blob) {
    app->store_.kv.putBytes(atoms3r_ical::ImuCalStoreNvs::kKeyV5, seg.cal.blob.data(), seg.cal.blob.size());
  }
  g_sink = serial;
  HardwareSerial::capture() = &sinkWrite;
  atoms3r_ical::sampleSource() = &replaySource;
  app->begin();

  uint32_t prev_us = seg.frames.front().s.sample_us;
  uint32_t expect = seg.frames.front().seq;
  bool ok = true;
  for (const auto& f : seg.frames) {
    if (f.seq != expect) {
      std::fprintf(stderr, "replay: %s R lines missing: seq %" PRIu32 "..%" PRIu32 "\n",
                   allow_gaps ? "continuing across" : "stopping at", expect, f.seq - 1);
      if (!allow_gaps) { ok = false; break; }
    }
    expect = f.seq + 1;
    if (std::find(seg.reset_events.begin(), seg.reset_events.end(), f.seq) != seg.reset_events.end() &&
        &f != &seg.frames.front()) {
      std::fprintf(stderr, "replay: the device reset its filter at seq %" PRIu32
                           " (recalibration or erase); stopping there\n", f.seq);
      break;
    }
    t64 += uint32_t(f.s.sample_us - prev_us);
    prev_us = f.s.sample_us;
    ins_sketch_host::clock_us() = t64;
    g_frame = &f;
    app->tick();
    g_frame = nullptr;
    if (outputs) outputs->push_back(outputsOf(*app));
  }
  HardwareSerial::capture() = nullptr;
  g_sink = nullptr;
  atoms3r_ical::sampleSource() = nullptr;
  return ok;
}

// -------------------------------------------------------------- self-test --

// Wave-following motion with a held 25 Hz magnetometer, accelerometer and gyro
// paired at 100 Hz behind 200 Hz polls (the device's accepted-frame cadence).
struct Synthetic {
  std::mt19937 rng{20261012};
  std::normal_distribution<float> n{0.0f, 1.0f};
  atoms3r_ical::MagAcquisition mag{};
  V3 mag_body = V3::Constant(NAN);
  long polls = 0;
  std::vector<ImuSample> delivered;
  std::vector<Outputs> outputs;

  V3 noise(float s) { return s * V3(n(rng), n(rng), n(rng)); }

  bool fill(uint32_t sample_us, ImuSample& out) {
    const double t = double(polls++) * 0.005;
    if (polls % 2 == 0) return false;
    auto att = [](double tt) {
      const float roll = 0.12f * float(std::sin(2 * M_PI * tt / 4.0));
      const float pitch = 0.05f * float(std::sin(2 * M_PI * tt / 3.0 + 0.7));
      const float yaw = 1.0f + 0.15f * float(std::sin(2 * M_PI * tt / 97.0));
      return (Eigen::AngleAxisf(yaw, V3::UnitZ()) * Eigen::AngleAxisf(pitch, V3::UnitY()) *
              Eigen::AngleAxisf(roll, V3::UnitX())).normalized();
    };
    const Eigen::Quaternionf q = att(t), q1 = att(t + 0.01);
    const Eigen::AngleAxisf aa(q.conjugate() * q1);
    const double wz = 2 * M_PI / 7.0;
    const float az = float(0.5 * wz * wz * std::sin(wz * t));
    const float g = atoms3r_ical::ImuCalCfg::g_cal_local;
    out.a = q.conjugate() * (V3(0, 0, az) - V3(0, 0, g)) + V3(0.015f, -0.01f, 0.02f) + noise(0.02f);
    out.w = aa.axis() * (aa.angle() / 0.01f) + V3(0.002f, -0.0015f, 0.001f) + noise(0.0015f);
    out.tempC = 35.0f + float(t) * 0.001f;
    if ((polls / 2) % 4 == 0) {
      mag.acquired(millis(), sample_us, sample_us, sample_us + 300);
      mag_body = q.conjugate() * V3(20.0f, -1.5f, 45.0f) + noise(0.3f);
    }
    out.m = mag.fresh(millis()) ? mag_body : V3::Constant(NAN);
    out.mag_timing = mag;
    delivered.push_back(out);
    return true;
  }
};
Synthetic* g_syn = nullptr;
bool syntheticSource(uint32_t us, ImuSample& out) { return g_syn->fill(us, out); }

int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-48s %12.4f  (limit %.4f)  %s\n", what, value, limit, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

bool sameBits(const V3& a, const V3& b) { return std::memcmp(a.data(), b.data(), sizeof(float) * 3) == 0; }

int selfTest(const char* capture_path) {
  std::printf("INS sketch diagnostic replay round trip: %s\n", INS_SKETCH_FILE);
  constexpr double kDurationS = 240.0;

  // Run A: the sketch on synthetic sensors, diagnostic stream captured.
  std::string device_log;
  Synthetic syn;
  g_syn = &syn;
  ins_sketch_host::clock_us() = 1000000;
  g_sink = &device_log;
  HardwareSerial::capture() = &sinkWrite;
  atoms3r_ical::sampleSource() = &syntheticSource;
  {
    auto app = std::make_unique<FusionApp>();
    app->begin();
    const long ticks = long(kDurationS / 0.005);
    for (long i = 0; i < ticks; ++i) {
      const size_t before = syn.delivered.size();
      app->tick();
      if (syn.delivered.size() != before) syn.outputs.push_back(outputsOf(*app));
    }
    check(INS_READY((*app)), "fused output ready at end of run A", INS_READY((*app)) ? 1 : 0, 1);
  }
  HardwareSerial::capture() = nullptr;
  g_sink = nullptr;
  atoms3r_ical::sampleSource() = nullptr;

  if (capture_path) {
    std::ofstream out(capture_path, std::ios::binary);
    out << device_log;
  }

  // The capture parses into exactly the delivered frames.
  std::istringstream in(device_log);
  const Capture cap = parseCapture(in);
  check(cap.bad_lines == 0, "malformed diagnostic lines", double(cap.bad_lines), 0);
  check(!cap.sketch.empty() && std::strstr(INS_SKETCH_FILE, cap.sketch.c_str()) != nullptr,
        "header (H) names this sketch", cap.sketch.empty() ? 0 : 1, 1);
  check(cap.d_lines >= long(kDurationS * 10) - 2, "state (D) lines at 10 Hz", double(cap.d_lines),
        kDurationS * 10 - 2);
  check(cap.segments.size() == 1, "boot segments", double(cap.segments.size()), 1);
  if (cap.segments.size() != 1) return 1;
  const Segment& seg = cap.segments.front();
  check(seg.cal.present, "calibration (C) line present", seg.cal.present ? 1 : 0, 1);
  check(seg.frames.size() == syn.delivered.size(), "R frames == delivered frames",
        double(seg.frames.size()), double(syn.delivered.size()));
  size_t bad_frames = 0;
  for (size_t i = 0; i < std::min(seg.frames.size(), syn.delivered.size()); ++i) {
    const auto& a = seg.frames[i].s;
    const auto& b = syn.delivered[i];
    const bool same = seg.frames[i].seq == i && a.sample_us == b.sample_us && sameBits(a.a, b.a) &&
                      sameBits(a.w, b.w) && sameBits(a.m, b.m) &&
                      std::memcmp(&a.tempC, &b.tempC, sizeof(float)) == 0 &&
                      a.mag_timing.have == b.mag_timing.have &&
                      a.mag_timing.sequence == b.mag_timing.sequence &&
                      a.mag_timing.last_ms == b.mag_timing.last_ms &&
                      a.mag_timing.frame_us == b.mag_timing.frame_us &&
                      a.mag_timing.read_start_us == b.mag_timing.read_start_us &&
                      a.mag_timing.read_end_us == b.mag_timing.read_end_us;
    if (!same) ++bad_frames;
  }
  check(bad_frames == 0, "R frames differing from delivered (bits)", double(bad_frames), 0);

  // Run B: replay the parsed capture through a fresh sketch.
  std::vector<Outputs> replayed;
  const bool ok = replay(seg, false, &replayed, nullptr);
  check(ok, "replay completed", ok ? 1 : 0, 1);
  size_t mismatches = 0, first = size_t(-1);
  for (size_t i = 0; i < std::min(replayed.size(), syn.outputs.size()); ++i) {
    if (!(replayed[i] == syn.outputs[i])) {
      if (first == size_t(-1)) first = i;
      ++mismatches;
    }
  }
  check(replayed.size() == syn.outputs.size(), "replayed frames", double(replayed.size()),
        double(syn.outputs.size()));
  check(mismatches == 0, "frames whose published outputs differ (bits)", double(mismatches), 0);
  if (first != size_t(-1)) {
    const auto& a = syn.outputs[first];
    const auto& b = replayed[first];
    std::printf("  first difference at frame %zu: heave %.9g vs %.9g, heading %.9g vs %.9g\n", first,
                double(a.heave_raw), double(b.heave_raw), double(a.heading), double(b.heading));
  }
  if (failures) {
    std::printf("INS sketch diagnostic replay: %d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("INS sketch diagnostic replay: PASS\n");
  return 0;
}

int replayFile(int argc, char** argv) {
  const char* path = argv[1];
  size_t segment = 0;
  bool allow_gaps = false, expect_ready = false;
  for (int i = 2; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--segment") && i + 1 < argc) segment = size_t(std::atoi(argv[++i]));
    else if (!std::strcmp(argv[i], "--allow-gaps")) allow_gaps = true;
    else if (!std::strcmp(argv[i], "--expect-ready")) expect_ready = true;
    else { std::fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
  }
  std::ifstream in(path);
  if (!in) { std::fprintf(stderr, "cannot open %s\n", path); return 2; }
  const Capture cap = parseCapture(in);
  std::fprintf(stderr, "capture: %ld lines, %ld malformed, %zu boot segment(s), sketch '%s'\n",
               cap.lines, cap.bad_lines, cap.segments.size(), cap.sketch.c_str());
  if (cap.acquisition == "fifo") {
    std::fprintf(stderr, "FIFO acquisition capture from %s, replayed through %s\n", cap.sketch.c_str(),
                 INS_SKETCH_FILE);
  } else if (!cap.sketch.empty() && !std::strstr(INS_SKETCH_FILE, cap.sketch.c_str())) {
    std::fprintf(stderr, "capture is from %s; this replay binary runs %s\n", cap.sketch.c_str(),
                 INS_SKETCH_FILE);
    return 2;
  }
  if (segment >= cap.segments.size()) { std::fprintf(stderr, "no segment %zu\n", segment); return 2; }
  const Segment& seg = cap.segments[segment];
  if (!seg.cal.present) std::fprintf(stderr, "warning: no C line; replaying without calibration\n");
  if (seg.frames.front().seq != 0) {
    std::fprintf(stderr, "warning: segment starts at seq %" PRIu32 ", not at boot: the replay starts "
                         "from a fresh filter, not the device's state\n", seg.frames.front().seq);
  }
  std::string serial;
  std::vector<Outputs> outputs;
  const bool ok = replay(seg, allow_gaps, &outputs, &serial);
  const auto first_ready = std::find_if(outputs.begin(), outputs.end(), [](const Outputs& o) { return o.ready; });
  const bool ready_at_end = !outputs.empty() && outputs.back().ready;
  std::fprintf(stderr, "replayed %zu frames; fused output %s, ready at end: %s\n", outputs.size(),
               first_ready == outputs.end() ? "never ready"
                                            : ("ready from frame " + std::to_string(first_ready - outputs.begin())).c_str(),
               ready_at_end ? "yes" : "no");
  // The replayed sketch's own stream, without the R lines it echoes.
  std::istringstream lines(serial);
  std::string line;
  while (std::getline(lines, line))
    if (line.rfind("R,", 0) != 0) std::printf("%s\n", line.c_str());
  return ok && (!expect_ready || ready_at_end) ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argv[1][0] == '-') {
    const char* capture_path = nullptr;
    for (int i = 1; i < argc; ++i) {
      if (!std::strcmp(argv[i], "--write-capture") && i + 1 < argc) capture_path = argv[++i];
      else if (std::strcmp(argv[i], "--self-test")) { std::fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    return selfTest(capture_path);
  }
  return replayFile(argc, argv);
}
