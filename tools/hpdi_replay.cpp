// hpdi_replay — replay a CSV record through hpdi::HeaveHPDI.
//
// Build:  g++ -std=c++14 -O2 -I src tools/hpdi_replay.cpp -o build/hpdi_replay
//
// Input CSV: header row, then numeric rows. Required columns: time (s) and
// up-positive vertical acceleration with gravity removed (m/s^2), i.e. the
// measurement-only Mahony-proxy channel of paper Eq. (73). Optional: reference
// displacement (m, up-positive unless --flip-z) for scoring, and the
// measurement-only zero-crossing period estimate (s) for PeriodScaled mode.
//
// Modes:
//   --dump OUT.csv  with --n --m --mode --value     writes t,z_hat,v_hat,fc_hz
//   --grid CFG.txt  --harness-window S | --score-last S | --score-from T
//                   writes one score row per config
//
// --harness-window S scores the trailing window exactly as the W3D simulation
// harness does (util/W3dSimCommon.cpp): the last N = floor(float(S) /
// float(dt)) rows (dt = --harness-dt, default 1/200 s), float error
// z_hat - z_ref, squares accumulated in order into a float with a fused
// multiply-add, rms = sqrt(sum / N) in float. The other two windows
// accumulate in double.
//
// Grid file: one config per line, "n m mode value" (mode: fixed|period;
// value: cutoff Hz or cutoff ratio f_c*T_z). '#' starts a comment.

#include "hpdi/HeaveHPDI.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Args {
  std::string in, dump, grid, out = "-";
  std::string t_col = "t", a_col = "a_up", z_col = "z_ref", tz_col = "Tz_hat";
  int n = 4, m = 3;
  std::string mode = "period";
  double value = 0.15;
  double score_from = NAN, score_last = NAN;
  double harness_window = NAN, harness_dt = 1.0 / 200.0;
  double prime_seconds = 0.0;
  double smoothing_periods = 0.75, refresh_s = 0.1, fallback_period_s = 6.0;
  bool flip_z = false, use_float = false;
};

[[noreturn]] void die(const std::string& msg) {
  std::fprintf(stderr, "hpdi_replay: %s\n", msg.c_str());
  std::exit(2);
}

void usage() {
  std::fprintf(stderr,
               "usage: hpdi_replay --in REC.csv (--dump OUT.csv | --grid CFG.txt [--out SCORES.csv])\n"
               "  [--t-col t] [--a-col a_up] [--z-col z_ref] [--tz-col Tz_hat|none]\n"
               "  [--n 4] [--m 3] [--mode fixed|period] [--value 0.15]\n"
               "  [--harness-window 900 [--harness-dt 0.005] | --score-last 900 | --score-from 300] [--prime-seconds 0] [--flip-z]\n"
               "  [--smoothing-periods 0.75] [--refresh 0.1] [--fallback-period 6] [--float]\n");
}

double toDouble(const char* s, const char* what) {
  errno = 0;
  char* end = nullptr;
  const double v = std::strtod(s, &end);
  if (errno != 0 || end == s) die(std::string("bad number for ") + what + ": " + s);
  return v;
}

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto next = [&](const char* name) -> const char* {
      if (i + 1 >= argc) die(std::string("missing value for ") + name);
      return argv[++i];
    };
    if (k == "--in") a.in = next("--in");
    else if (k == "--dump") a.dump = next("--dump");
    else if (k == "--grid") a.grid = next("--grid");
    else if (k == "--out") a.out = next("--out");
    else if (k == "--t-col") a.t_col = next("--t-col");
    else if (k == "--a-col") a.a_col = next("--a-col");
    else if (k == "--z-col") a.z_col = next("--z-col");
    else if (k == "--tz-col") a.tz_col = next("--tz-col");
    else if (k == "--n") a.n = static_cast<int>(toDouble(next("--n"), "--n"));
    else if (k == "--m") a.m = static_cast<int>(toDouble(next("--m"), "--m"));
    else if (k == "--mode") a.mode = next("--mode");
    else if (k == "--value") a.value = toDouble(next("--value"), "--value");
    else if (k == "--score-from") a.score_from = toDouble(next("--score-from"), "--score-from");
    else if (k == "--score-last") a.score_last = toDouble(next("--score-last"), "--score-last");
    else if (k == "--harness-window") a.harness_window = toDouble(next("--harness-window"), "--harness-window");
    else if (k == "--harness-dt") a.harness_dt = toDouble(next("--harness-dt"), "--harness-dt");
    else if (k == "--prime-seconds") a.prime_seconds = toDouble(next("--prime-seconds"), "--prime-seconds");
    else if (k == "--smoothing-periods") a.smoothing_periods = toDouble(next("--smoothing-periods"), "--smoothing-periods");
    else if (k == "--refresh") a.refresh_s = toDouble(next("--refresh"), "--refresh");
    else if (k == "--fallback-period") a.fallback_period_s = toDouble(next("--fallback-period"), "--fallback-period");
    else if (k == "--flip-z") a.flip_z = true;
    else if (k == "--float") a.use_float = true;
    else if (k == "-h" || k == "--help") { usage(); std::exit(0); }
    else die("unknown argument " + k);
  }
  if (a.in.empty() || (a.dump.empty() == a.grid.empty())) { usage(); die("need --in and exactly one of --dump/--grid"); }
  return a;
}

struct Record {
  std::vector<double> t, a, z, tz;
  bool has_z = false, has_tz = false;
};

std::vector<std::string> splitCsv(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  for (char ch : line) {
    if (ch == ',') { out.push_back(cur); cur.clear(); }
    else if (ch != '\r' && ch != '"' && ch != ' ') cur.push_back(ch);
  }
  out.push_back(cur);
  return out;
}

Record load(const Args& args) {
  std::ifstream f(args.in);
  if (!f) die("cannot open " + args.in);
  std::string line;
  if (!std::getline(f, line)) die("empty file " + args.in);
  const auto hdr = splitCsv(line);
  auto find = [&](const std::string& name) -> int {
    for (size_t i = 0; i < hdr.size(); ++i) if (hdr[i] == name) return static_cast<int>(i);
    return -1;
  };
  const int it = find(args.t_col), ia = find(args.a_col);
  const int iz = find(args.z_col), itz = args.tz_col == "none" ? -1 : find(args.tz_col);
  if (it < 0) die("time column '" + args.t_col + "' not found");
  if (ia < 0) die("acceleration column '" + args.a_col + "' not found");
  Record r;
  r.has_z = iz >= 0;
  r.has_tz = itz >= 0;
  const int need = std::max(std::max(it, ia), std::max(iz, itz));
  size_t row = 1;
  while (std::getline(f, line)) {
    ++row;
    if (line.empty()) continue;
    const auto cells = splitCsv(line);
    if (static_cast<int>(cells.size()) <= need) die("short row " + std::to_string(row));
    auto num = [&](int idx) { return idx < 0 ? kNaN : std::strtod(cells[idx].c_str(), nullptr); };
    r.t.push_back(num(it));
    r.a.push_back(num(ia));
    r.z.push_back(args.flip_z ? -num(iz) : num(iz));
    r.tz.push_back(num(itz));
  }
  if (r.t.size() < 2) die("record has fewer than 2 rows");
  for (size_t k = 1; k < r.t.size(); ++k)
    if (!(r.t[k] > r.t[k - 1])) die("time not strictly increasing at row " + std::to_string(k + 2));
  return r;
}

struct Cfg {
  int n, m;
  bool period;
  double value;
};

std::vector<Cfg> loadGrid(const std::string& path) {
  std::ifstream f(path);
  if (!f) die("cannot open grid " + path);
  std::vector<Cfg> out;
  std::string line;
  while (std::getline(f, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) line.resize(hash);
    std::istringstream ss(line);
    Cfg c{};
    std::string mode;
    if (!(ss >> c.n >> c.m >> mode >> c.value)) continue;
    if (mode != "fixed" && mode != "period") die("grid mode must be fixed|period: " + mode);
    c.period = mode == "period";
    out.push_back(c);
  }
  if (out.empty()) die("grid is empty: " + path);
  return out;
}

hpdi::HeaveHPDIConfig makeConfig(const Cfg& c, const Args& a) {
  hpdi::HeaveHPDIConfig h;
  h.order = c.n;
  h.dc_zeros = c.m;
  h.mode = c.period ? hpdi::CutoffMode::PeriodScaled : hpdi::CutoffMode::Fixed;
  if (c.period) h.cutoff_ratio = c.value; else h.cutoff_hz = c.value;
  h.smoothing_periods = a.smoothing_periods;
  h.refresh_period_s = a.refresh_s;
  h.fallback_period_s = a.fallback_period_s;
  return h;
}

struct RunOut {
  std::vector<double> z, v, fc;
  unsigned rejected = 0;
};

template <typename T>
RunOut run(const Record& r, const Cfg& c, const Args& a, bool keep) {
  const hpdi::HeaveHPDIConfig hc = makeConfig(c, a);
  if (const char* err = hpdi::validateConfig(hc, 6)) die(std::string("invalid config: ") + err);
  hpdi::HeaveHPDI<T, 6> f(hc);
  if (a.prime_seconds > 0.0) {
    double s = 0.0;
    size_t cnt = 0;
    for (size_t k = 0; k < r.t.size() && r.t[k] - r.t[0] <= a.prime_seconds; ++k)
      if (std::isfinite(r.a[k])) { s += r.a[k]; ++cnt; }
    if (cnt) f.reset(static_cast<T>(s / static_cast<double>(cnt)));
  }
  RunOut o;
  if (keep) { o.z.reserve(r.t.size()); o.v.reserve(r.t.size()); o.fc.reserve(r.t.size()); }
  else o.z.reserve(r.t.size());
  for (size_t k = 0; k < r.t.size(); ++k) {
    if (r.has_tz && std::isfinite(r.tz[k])) f.setWavePeriod(static_cast<T>(r.tz[k]));
    if (k > 0) f.update(static_cast<T>(r.a[k]), static_cast<T>(r.t[k] - r.t[k - 1]));
    o.z.push_back(static_cast<double>(f.displacement()));
    if (keep) { o.v.push_back(static_cast<double>(f.velocity())); o.fc.push_back(static_cast<double>(f.appliedCutoffHz())); }
  }
  o.rejected = f.rejectedSamples();
  return o;
}

RunOut runAny(const Record& r, const Cfg& c, const Args& a, bool keep) {
  return a.use_float ? run<float>(r, c, a, keep) : run<double>(r, c, a, keep);
}

}  // namespace

int main(int argc, char** argv) {
  const Args args = parse(argc, argv);
  const Record rec = load(args);

  if (!args.dump.empty()) {
    if (args.mode != "fixed" && args.mode != "period") die("--mode must be fixed|period");
    const Cfg c{args.n, args.m, args.mode == "period", args.value};
    if (c.period && !rec.has_tz) die("period mode needs the T_z column ('" + args.tz_col + "')");
    const RunOut o = runAny(rec, c, args, true);
    FILE* fp = args.dump == "-" ? stdout : std::fopen(args.dump.c_str(), "w");
    if (!fp) die("cannot write " + args.dump);
    std::fprintf(fp, "t,z_hat,v_hat,fc_hz\n");
    for (size_t k = 0; k < rec.t.size(); ++k)
      std::fprintf(fp, "%.6f,%.17g,%.9g,%.9g\n", rec.t[k], o.z[k], o.v[k], o.fc[k]);
    if (fp != stdout) std::fclose(fp);
    return 0;
  }

  if (!rec.has_z) die("grid scoring needs the reference column '" + args.z_col + "'");
  double t0 = args.score_from;
  if (std::isfinite(args.score_last)) t0 = rec.t.back() - args.score_last;
  const bool harness = std::isfinite(args.harness_window);
  if (!harness && !std::isfinite(t0)) die("give --harness-window, --score-last or --score-from");
  size_t first = 0;
  if (harness) {
    if (!(args.harness_window > 0.0) || !(args.harness_dt > 0.0)) die("harness window and dt must be positive");
    const auto requested = static_cast<size_t>(static_cast<float>(args.harness_window) / static_cast<float>(args.harness_dt));
    first = rec.t.size() - std::min(rec.t.size(), std::max<size_t>(requested, 1));
  }

  const auto grid = loadGrid(args.grid);
  FILE* fp = args.out == "-" ? stdout : std::fopen(args.out.c_str(), "w");
  if (!fp) die("cannot write " + args.out);
  std::fprintf(fp, "n,m,mode,value,rms_m,mean_err_m,samples,rejected\n");
  for (const Cfg& c : grid) {
    if (c.period && !rec.has_tz) die("period mode needs the T_z column ('" + args.tz_col + "')");
    const RunOut o = runAny(rec, c, args, false);
    double se = 0.0, sm = 0.0;
    float se_f = 0.0f;
    size_t cnt = 0;
    for (size_t k = first; k < rec.t.size(); ++k) {
      if (harness) {
        const float e = static_cast<float>(o.z[k]) - static_cast<float>(rec.z[k]);
        se_f = std::fma(e, e, se_f);
        sm += static_cast<double>(e);
        ++cnt;
        continue;
      }
      if (rec.t[k] < t0 || !std::isfinite(rec.z[k]) || !std::isfinite(o.z[k])) continue;
      const double e = o.z[k] - rec.z[k];
      se += e * e;
      sm += e;
      ++cnt;
    }
    const double rms = !cnt ? kNaN
                       : harness ? static_cast<double>(std::sqrt(se_f / static_cast<float>(cnt)))
                                 : std::sqrt(se / static_cast<double>(cnt));
    const double mean = cnt ? sm / static_cast<double>(cnt) : kNaN;
    std::fprintf(fp, "%d,%d,%s,%.9g,%.9g,%.9g,%zu,%u\n", c.n, c.m, c.period ? "period" : "fixed", c.value, rms, mean,
                 cnt, o.rejected);
  }
  if (fp != stdout) std::fclose(fp);
  return 0;
}
