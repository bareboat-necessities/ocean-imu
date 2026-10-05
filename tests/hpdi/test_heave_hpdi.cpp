// Unit tests for hpdi::HeaveHPDI. Plain main(); non-zero exit on failure.
// Build: g++ -std=c++14 -O2 -Wall -Wextra -I src tests/hpdi/test_heave_hpdi.cpp -o test_heave_hpdi

#include "hpdi/HeaveHPDI.h"

#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond, ...)                                     \
  do {                                                       \
    ++g_checks;                                              \
    if (!(cond)) {                                           \
      ++g_failures;                                          \
      std::printf("FAIL %s:%d: %s | ", __FILE__, __LINE__, #cond); \
      std::printf(__VA_ARGS__);                              \
      std::printf("\n");                                     \
    }                                                        \
  } while (0)

constexpr double kPi = 3.14159265358979323846;

hpdi::HeaveHPDIConfig fixedCfg(int n, int m, double fc) {
  hpdi::HeaveHPDIConfig c;
  c.order = n;
  c.dc_zeros = m;
  c.mode = hpdi::CutoffMode::Fixed;
  c.cutoff_hz = fc;
  return c;
}

// Least-squares fit y ~ alpha sin(Wt) + beta cos(Wt); returns alpha + j beta,
// the complex gain relative to a unit sin(Wt) reference.
std::complex<double> fitGain(const std::vector<double>& t, const std::vector<double>& y, double W) {
  double ss = 0, cc = 0, sc = 0, ys = 0, yc = 0;
  for (size_t k = 0; k < t.size(); ++k) {
    const double s = std::sin(W * t[k]), c = std::cos(W * t[k]);
    ss += s * s; cc += c * c; sc += s * c; ys += y[k] * s; yc += y[k] * c;
  }
  const double det = ss * cc - sc * sc;
  return {(ys * cc - yc * sc) / det, (yc * ss - ys * sc) / det};
}

void testButterworth() {
  double c[7];
  hpdi::butterworthCoeffs(4, c);
  CHECK(std::fabs(c[1] - 2.6131259297527) < 1e-12, "c1=%.15f", c[1]);
  CHECK(std::fabs(c[2] - 3.4142135623731) < 1e-12, "c2=%.15f", c[2]);
  CHECK(std::fabs(c[4] - 1.0) < 1e-12, "c4=%.15f", c[4]);
  hpdi::butterworthCoeffs(6, c);
  CHECK(std::fabs(c[3] - 9.1416201726) < 1e-9, "c3=%.12f", c[3]);
  CHECK(std::fabs(c[6] - 1.0) < 1e-12, "c6=%.15f", c[6]);
}

void testValidation() {
  hpdi::HeaveHPDIConfig bad = fixedCfg(4, 5, 0.05);
  CHECK(hpdi::validateConfig(bad, 6) != nullptr, "m > n accepted");
  bad = fixedCfg(7, 3, 0.05);
  CHECK(hpdi::validateConfig(bad, 6) != nullptr, "n > NMax accepted");
  bad = fixedCfg(4, 1, 0.05);
  CHECK(hpdi::validateConfig(bad, 6) != nullptr, "m < 2 accepted");
  hpdi::HeaveHPDI<double> f(fixedCfg(4, 5, 0.05));
  CHECK(f.lastError() != nullptr, "constructor did not report error");
  CHECK(f.order() == 4 && f.dcZeros() == 3, "fallback to defaults failed");
  CHECK(!f.configure(fixedCfg(3, 4, 0.05)) && f.order() == 4, "invalid configure changed state");
}

// Steady-state complex gain against the continuous-time G(j w).
void testFrequencyResponse() {
  const double fs = 200.0, h = 1.0 / fs, fc = 0.05;
  const int nm[][2] = {{3, 2}, {4, 3}, {4, 4}, {5, 3}, {6, 3}, {6, 4}};
  const double freqs[] = {0.08, 0.15, 0.4};
  for (const auto& p : nm) {
    for (double f : freqs) {
      hpdi::HeaveHPDI<double> filt(fixedCfg(p[0], p[1], fc));
      const double W = 2 * kPi * f;
      const double T_end = 600.0, T_fit = 300.0;
      std::vector<double> tt, yy, vv;
      const int N = static_cast<int>(T_end * fs);
      for (int k = 1; k <= N; ++k) {
        const double t = k * h;
        filt.update(-W * W * std::sin(W * t), h);
        if (t > T_end - T_fit) { tt.push_back(t); yy.push_back(filt.displacement()); vv.push_back(filt.velocity()); }
      }
      const std::complex<double> Gest = fitGain(tt, yy, W);
      const std::complex<double> Gref = hpdi::transferG(f, fc, p[0], p[1]);
      const double err = std::abs(Gest - Gref);
      CHECK(err < 1e-3 * std::max(1.0, std::abs(Gref)), "n=%d m=%d f=%.2f |G-Gref|=%.2e", p[0], p[1], f, err);
      // velocity: d/dt of G sin -> W * j * G relative to sin(Wt)
      const std::complex<double> Vest = fitGain(tt, vv, W);
      const std::complex<double> Vref = std::complex<double>(0.0, W) * Gref;
      CHECK(std::abs(Vest - Vref) < 2e-3 * std::abs(Vref), "velocity n=%d m=%d f=%.2f err=%.2e", p[0], p[1], f,
            std::abs(Vest - Vref) / std::abs(Vref));
    }
  }
}

// DC structure: constant bias b -> offset c_2 b / w^2 for m = 2, zero for m >= 3;
// ramp k t -> offset c_3 k / w^3 for m = 3, zero for m >= 4.
void testLowFrequencyStructure() {
  const double fs = 50.0, h = 1.0 / fs, fc = 0.05, w = 2 * kPi * fc;
  double c[7];
  {
    hpdi::HeaveHPDI<double> f2(fixedCfg(3, 2, fc)), f3(fixedCfg(4, 3, fc));
    const double b = 0.05;
    for (int k = 0; k < static_cast<int>(3000 * fs); ++k) { f2.update(b, h); f3.update(b, h); }
    hpdi::butterworthCoeffs(3, c);
    const double expect2 = c[2] * b / (w * w);
    CHECK(std::fabs(f2.displacement() - expect2) < 1e-6 * expect2, "m=2 bias offset %.9f vs %.9f", f2.displacement(), expect2);
    CHECK(std::fabs(f3.displacement()) < 1e-9, "m=3 bias residual %.3e", f3.displacement());
  }
  {
    hpdi::HeaveHPDI<double> f3(fixedCfg(4, 3, fc)), f4(fixedCfg(5, 4, fc));
    const double kr = 1e-4;
    for (int k = 1; k <= static_cast<int>(4000 * fs); ++k) { const double a = kr * k * h; f3.update(a, h); f4.update(a, h); }
    hpdi::butterworthCoeffs(4, c);
    const double expect3 = c[3] * kr / (w * w * w);
    CHECK(std::fabs(f3.displacement() - expect3) < 1e-4 * expect3, "m=3 ramp offset %.9f vs %.9f", f3.displacement(), expect3);
    CHECK(std::fabs(f4.displacement()) < 1e-7, "m=4 ramp residual %.3e", f4.displacement());
  }
}

void testPriming() {
  const double h = 0.005, b = 0.07;
  hpdi::HeaveHPDI<double> f(fixedCfg(4, 3, 0.04));
  f.reset(b);
  double worst = 0.0;
  for (int k = 0; k < 20000; ++k) { f.update(b, h); worst = std::max(worst, std::fabs(f.displacement())); }
  CHECK(worst < 1e-9, "primed m=3 transient %.3e", worst);

  hpdi::HeaveHPDI<double> g(fixedCfg(3, 2, 0.04));
  g.reset(b);
  double c[4];
  hpdi::butterworthCoeffs(3, c);
  const double w = 2 * kPi * 0.04, expect = c[2] * b / (w * w);
  double dev = 0.0;
  for (int k = 0; k < 20000; ++k) { g.update(b, h); dev = std::max(dev, std::fabs(g.displacement() - expect)); }
  CHECK(dev < 1e-9 * expect, "primed m=2 deviation %.3e", dev);
}

void testScheduling() {
  hpdi::HeaveHPDIConfig cfg;
  cfg.order = 4; cfg.dc_zeros = 3; cfg.mode = hpdi::CutoffMode::PeriodScaled;
  cfg.cutoff_ratio = 0.15; cfg.fallback_period_s = 6.0;
  hpdi::HeaveHPDI<double> f(cfg);
  const double h = 0.005;
  f.update(0.0, h);
  CHECK(std::fabs(f.appliedCutoffHz() - 0.15 / 6.0) < 1e-12, "fallback cutoff %.6f", f.appliedCutoffHz());
  f.setWavePeriod(4.0);
  f.update(0.0, h);
  CHECK(std::fabs(f.appliedCutoffHz() - 0.15 / 4.0) < 1e-12, "snap cutoff %.6f", f.appliedCutoffHz());
  f.setWavePeriod(8.0);
  double prev = f.appliedCutoffHz();
  int changes = 0;
  bool monotone = true;
  for (int k = 0; k < 20000; ++k) {  // 100 s
    f.update(0.0, h);
    const double fc = f.appliedCutoffHz();
    if (fc != prev) { ++changes; if (fc > prev) monotone = false; prev = fc; }
  }
  CHECK(monotone, "cutoff not monotone toward target");
  CHECK(changes <= 1001 && changes >= 900, "refresh cadence: %d changes in 100 s", changes);
  CHECK(std::fabs(f.appliedCutoffHz() - 0.15 / 8.0) < 1e-6, "converged cutoff %.8f", f.appliedCutoffHz());
  f.setWavePeriod(-1.0);
  f.setWavePeriod(std::nan(""));
  f.update(0.0, h);
  CHECK(std::fabs(f.appliedCutoffHz() - 0.15 / 8.0) < 1e-6, "invalid period changed cutoff");
}

// After a period step and settling, output matches a filter run at the final cutoff.
void testScheduledSettling() {
  hpdi::HeaveHPDIConfig cfg;
  cfg.order = 5; cfg.dc_zeros = 3; cfg.mode = hpdi::CutoffMode::PeriodScaled; cfg.cutoff_ratio = 0.15;
  hpdi::HeaveHPDI<double> sched(cfg), ref(cfg);
  ref.setWavePeriod(8.0);
  sched.setWavePeriod(5.0);
  const double h = 0.005, W1 = 2 * kPi / 7.0, W2 = 2 * kPi / 4.3;
  double err = 0.0;
  for (int k = 1; k <= static_cast<int>(900.0 / h); ++k) {
    const double t = k * h;
    if (t > 300.0) sched.setWavePeriod(8.0);
    const double a = -W1 * W1 * std::sin(W1 * t) - 0.3 * W2 * W2 * std::sin(W2 * t + 0.4);
    sched.update(a, h);
    ref.update(a, h);
    if (t > 700.0) err = std::max(err, std::fabs(sched.displacement() - ref.displacement()));
  }
  CHECK(err < 1e-4, "scheduled vs reference after settling: %.3e", err);
}

void testVariableDt() {
  const double fc = 0.04, f = 0.15, W = 2 * kPi * f;
  hpdi::HeaveHPDI<double> filt(fixedCfg(4, 3, fc));
  uint32_t lcg = 12345u;
  double t = 0.0;
  std::vector<double> tt, yy;
  while (t < 600.0) {
    lcg = lcg * 1664525u + 1013904223u;
    const double jitter = 0.8 + 0.4 * ((lcg >> 8) / 16777216.0);  // +-20 %
    const double dt = 0.005 * jitter;
    t += dt;
    filt.update(-W * W * std::sin(W * t), dt);
    if (t > 300.0) { tt.push_back(t); yy.push_back(filt.displacement()); }
  }
  const std::complex<double> Gest = fitGain(tt, yy, W), Gref = hpdi::transferG(f, fc, 4, 3);
  CHECK(std::abs(Gest - Gref) < 2e-3, "jittered dt gain error %.2e", std::abs(Gest - Gref));
}

void testFloatPrecision() {
  const double h = 0.005;
  hpdi::HeaveHPDI<double> fd(fixedCfg(4, 3, 0.02));
  hpdi::HeaveHPDI<float> ff(fixedCfg(4, 3, 0.02));
  const double W1 = 2 * kPi / 9.0, W2 = 2 * kPi / 5.5;
  double se = 0.0, sr = 0.0;
  for (int k = 1; k <= static_cast<int>(1200.0 / h); ++k) {
    const double t = k * h;
    const double a = -W1 * W1 * std::sin(W1 * t) - 0.5 * W2 * W2 * std::sin(W2 * t + 1.0) + 0.03;
    fd.update(a, h);
    ff.update(static_cast<float>(a), static_cast<float>(h));
    if (t > 300.0) {
      const double d = fd.displacement() - static_cast<double>(ff.displacement());
      se += d * d;
      sr += fd.displacement() * fd.displacement();
    }
  }
  CHECK(std::sqrt(se / sr) < 1e-3, "float relative RMS deviation %.2e", std::sqrt(se / sr));
}

void testRejection() {
  hpdi::HeaveHPDI<double> f(fixedCfg(4, 3, 0.04));
  for (int k = 0; k < 1000; ++k) f.update(std::sin(0.01 * k), 0.005);
  const double y0 = f.displacement();
  CHECK(!f.update(std::nan(""), 0.005), "NaN accepted");
  CHECK(!f.update(1.0, -0.005), "negative dt accepted");
  CHECK(!f.update(1.0, 10.0), "huge dt accepted");
  CHECK(f.displacement() == y0, "state changed on rejection");
  CHECK(f.rejectedSamples() == 3, "rejected count %u", f.rejectedSamples());
}

}  // namespace

int main() {
  testButterworth();
  testValidation();
  testFrequencyResponse();
  testLowFrequencyStructure();
  testPriming();
  testScheduling();
  testScheduledSettling();
  testVariableDt();
  testFloatPrecision();
  testRejection();
  std::printf("%d/%d checks passed\n", g_checks - g_failures, g_checks);
  return g_failures == 0 ? 0 : 1;
}
