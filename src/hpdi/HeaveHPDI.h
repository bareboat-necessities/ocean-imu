#pragma once
/*
 * HeaveHPDI.h — causal high-pass double-integration (HPDI) heave baseline.
 *
 * Conventional comparator for the OU–III wave-motion estimator: an LTI heave
 * filter driven by the same measurement-only vertical acceleration (Mahony
 * proxy, paper Eq. 73) and, optionally, the same measurement-only zero-crossing
 * period estimate (paper Eq. 80). No reference motion enters online.
 *
 * Continuous-time design, cutoff w = 2*pi*f_c, input a (m/s^2), output y (m):
 *
 *   D(s) = sum_{k=0..n}   c_k w^(n-k) s^k   normalized Butterworth, c_n = c_0 = 1
 *   R(s) = sum_{k=0..m-1} c_k w^(n-k) s^k   lowest m terms of D
 *   G(s) = 1 - R(s)/D(s)                    displacement gain relative to 1/s^2
 *   H(s) = G(s)/s^2 = [D(s) - R(s)] / (s^2 D(s))
 *
 * Because D - R = s^m Q(s), H has only the stable poles of D: no free
 * integrators, so internal states stay bounded.
 *
 *   m == n  classic cascaded Butterworth high-pass + double integration,
 *           G = s^n/D. High-frequency error |G - 1| ~ c_{n-1} (w/W): first order,
 *           i.e. a phase lead that does not vanish quickly above the corner.
 *   m <  n  phase-compensated (complementary) form. |G - 1| ~ c_{m-1} (w/W)^(n-m+1),
 *           at the cost of gain peaking (|G| > 1) near the corner.
 *   m >= 3  zero steady-state displacement for a constant acceleration bias and
 *           bounded variance under a random-walk bias. m == 2 leaves a constant
 *           offset c_2 b / w^2 (n = 3, m = 2 reproduces the paper's reduced
 *           integral-regularizer response, Eq. 86, with w = w_R).
 *
 * Realization: observer canonical form in frequency-normalized coordinates
 * z_i = x_i / w^i (every state in metres):
 *
 *   dz/dt = w * (Ac z + beta * a / w^2),   y = z_0,   v = w (z_1 - c_{n-1} z_0)
 *
 * discretized by the trapezoidal rule (bilinear transform; identical to
 * scipy.signal.bilinear of H(s) for fixed w and h) in delta form,
 *   z <- z + F z + G (a_prev + a)/2,  F = P M,  G = P beta h / w,
 *   M = h w Ac,  P = (I - M/2)^-1,
 * which keeps float precision when w h << 1. F and G are recomputed on a
 * sample-and-hold cadence (default 0.1 s, as for OU–III) or when h changes.
 * On a cutoff change z_1..z_{n-1} are scaled by w_old/w_new, which preserves the
 * leading-order steady state of wave components well above the corner.
 *
 * Conventions: input up-positive, gravity removed; outputs up-positive.
 * Negate for NED-down references.
 *
 * Embedded use: no heap, no exceptions. T = float suits FPU-only targets such
 * as ESP32-S3; T = double for offline evaluation. Per-sample cost n^2 + 2n MACs.
 */

#include <cmath>
#include <complex>
#include <cstdint>

namespace hpdi {

enum class CutoffMode : uint8_t {
  Fixed = 0,        // f_c = cutoff_hz
  PeriodScaled = 1  // f_c = cutoff_ratio / T_z, T_z from the measurement-only estimator
};

struct HeaveHPDIConfig {
  int order = 4;     // n, Butterworth order of D(s), 2..NMax
  int dc_zeros = 3;  // m, G(s) ~ s^m at DC, 2..n; m == n is the classic HP form

  CutoffMode mode = CutoffMode::PeriodScaled;
  double cutoff_hz = 0.04;          // Fixed mode
  double cutoff_ratio = 0.15;       // PeriodScaled mode: f_c * T_z
  double fallback_period_s = 6.0;   // used until the first valid T_z

  double f_min_hz = 0.005;
  double f_max_hz = 0.5;

  // Log-cutoff EMA horizon in units of T_z (PeriodScaled only), then guarded.
  double smoothing_periods = 0.75;
  double smoothing_min_s = 0.05;
  double smoothing_max_s = 35.0;

  double refresh_period_s = 0.1;  // coefficient sample-and-hold cadence
  double dt_rel_tol = 1e-4;       // recompute F, G when |h - h_applied|/h_applied exceeds this
  double dt_min_s = 1e-5;         // samples with dt outside [dt_min, dt_max] are rejected
  double dt_max_s = 0.1;
  double max_wh = 0.05;           // cap w*h so the corner stays far below Nyquist
};

// Returns nullptr if valid, otherwise a static description of the first problem.
inline const char* validateConfig(const HeaveHPDIConfig& c, int n_max) noexcept {
  if (c.order < 2 || c.order > n_max) return "order outside [2, NMax]";
  if (c.dc_zeros < 2 || c.dc_zeros > c.order) return "dc_zeros outside [2, order]";
  if (!(c.f_min_hz > 0.0) || !(c.f_max_hz > c.f_min_hz)) return "need 0 < f_min_hz < f_max_hz";
  if (c.mode == CutoffMode::Fixed && !(c.cutoff_hz > 0.0)) return "cutoff_hz must be > 0";
  if (c.mode == CutoffMode::PeriodScaled && !(c.cutoff_ratio > 0.0)) return "cutoff_ratio must be > 0";
  if (!(c.fallback_period_s > 0.0)) return "fallback_period_s must be > 0";
  if (!(c.smoothing_periods >= 0.0)) return "smoothing_periods must be >= 0";
  if (!(c.smoothing_min_s > 0.0) || !(c.smoothing_max_s >= c.smoothing_min_s)) return "bad smoothing bounds";
  if (!(c.refresh_period_s >= 0.0)) return "refresh_period_s must be >= 0";
  if (!(c.dt_rel_tol >= 0.0)) return "dt_rel_tol must be >= 0";
  if (!(c.dt_min_s > 0.0) || !(c.dt_max_s > c.dt_min_s)) return "bad dt bounds";
  if (!(c.max_wh > 0.0) || !(c.max_wh <= 0.5)) return "max_wh must be in (0, 0.5]";
  return nullptr;
}

// Normalized Butterworth polynomial coefficients, ascending: c[0..n], c[0] = c[n] = 1.
// c_k = prod_{j=1..k} cos((j-1) g) / sin(j g),  g = pi / (2n).
inline void butterworthCoeffs(int n, double* c) noexcept {
  const double g = 3.14159265358979323846 / (2.0 * n);
  c[0] = 1.0;
  for (int k = 1; k <= n; ++k) c[k] = c[k - 1] * std::cos((k - 1) * g) / std::sin(k * g);
}

// Continuous-time displacement gain G(j 2 pi f) for diagnostics and tests.
inline std::complex<double> transferG(double f_hz, double fc_hz, int n, int m) noexcept {
  double c[17];
  if (n < 1 || n > 16 || m < 0 || m > n) return {0.0, 0.0};
  butterworthCoeffs(n, c);
  const std::complex<double> s(0.0, f_hz / fc_hz);  // normalized: w = 1
  std::complex<double> D(0.0, 0.0), R(0.0, 0.0), sk(1.0, 0.0);
  for (int k = 0; k <= n; ++k) {
    D += c[k] * sk;
    if (k < m) R += c[k] * sk;
    sk *= s;
  }
  return 1.0 - R / D;
}

template <typename T = float, int NMax = 6>
class HeaveHPDI {
  static_assert(NMax >= 2 && NMax <= 8, "NMax must be in [2, 8]");

 public:
  HeaveHPDI() noexcept { configure(HeaveHPDIConfig{}); }

  // On an invalid config the defaults are used and lastError() is set.
  explicit HeaveHPDI(const HeaveHPDIConfig& cfg) noexcept {
    if (!configure(cfg)) {
      const char* err = error_;
      configure(HeaveHPDIConfig{});
      error_ = err;
    }
  }

  // Applies a new configuration and resets the state. Returns false and keeps
  // the previous configuration if invalid.
  bool configure(const HeaveHPDIConfig& cfg) noexcept {
    const char* err = validateConfig(cfg, NMax);
    if (err) {
      error_ = err;
      return false;
    }
    error_ = nullptr;
    cfg_ = cfg;
    n_ = cfg.order;
    m_ = cfg.dc_zeros;

    double c[NMax + 1];
    butterworthCoeffs(n_, c);
    for (int k = 0; k <= n_; ++k) c_[k] = static_cast<T>(c[k]);
    for (int i = 0; i < NMax; ++i) beta_[i] = T(0);
    for (int i = 1; i <= n_ + 1 - m_; ++i) beta_[i] = static_cast<T>(c[n_ + 1 - i]);

    cutoff_hz_ = static_cast<T>(cfg_.cutoff_hz);
    cutoff_ratio_ = static_cast<T>(cfg_.cutoff_ratio);
    fallback_period_ = static_cast<T>(cfg_.fallback_period_s);
    f_min_ = static_cast<T>(cfg_.f_min_hz);
    f_max_ = static_cast<T>(cfg_.f_max_hz);
    smoothing_periods_ = static_cast<T>(cfg_.smoothing_periods);
    smoothing_min_ = static_cast<T>(cfg_.smoothing_min_s);
    smoothing_max_ = static_cast<T>(cfg_.smoothing_max_s);
    dt_min_ = static_cast<T>(cfg_.dt_min_s);
    dt_max_ = static_cast<T>(cfg_.dt_max_s);
    dt_rel_tol_ = static_cast<T>(cfg_.dt_rel_tol);
    refresh_period_ = static_cast<T>(cfg_.refresh_period_s);
    have_tz_ = false;
    tz_ = static_cast<T>(cfg_.fallback_period_s);
    updateTargets();
    log_fc_ = log_target_;
    w_applied_ = T(0);
    h_applied_ = T(0);
    reset();
    return true;
  }

  // Zeroes the filter state. A nonzero accel_bias_guess primes the steady state
  // for that constant input, suppressing the start-up transient.
  void reset(T accel_bias_guess = T(0)) noexcept {
    for (int i = 0; i < NMax; ++i) {
      z_[i] = T(0);
      for (int j = 0; j < NMax; ++j) F_[i][j] = T(0);
      G_[i] = T(0);
    }
    u_prev_ = accel_bias_guess;
    have_prev_ = false;
    refresh_timer_ = T(0);
    force_refresh_ = true;
    rejected_ = 0;
    updates_ = 0;
    if (accel_bias_guess != T(0) && std::isfinite(static_cast<double>(accel_bias_guess))) {
      primeSteadyState(accel_bias_guess);
      have_prev_ = true;
    }
  }

  // Measurement-only zero-crossing period (s). Non-finite or non-positive
  // values are ignored. The first valid value snaps the cutoff to its target.
  void setWavePeriod(T tz_s) noexcept {
    if (!(std::isfinite(tz_s) && tz_s > T(0))) return;
    if (have_tz_ && tz_s == tz_) return;
    tz_ = tz_s;
    const bool first = !have_tz_;
    have_tz_ = true;
    updateTargets();
    if (first && cfg_.mode == CutoffMode::PeriodScaled) {
      log_fc_ = log_target_;
      force_refresh_ = true;
    }
  }

  // One sample. a_up: up-positive vertical acceleration with gravity removed
  // (m/s^2); dt: time since the previous sample (s). Returns false if the
  // sample was rejected; the state is then unchanged.
  bool update(T a_up, T dt) noexcept {
    if (!std::isfinite(a_up) || !std::isfinite(dt) || dt < dt_min_ || dt > dt_max_) {
      ++rejected_;
      return false;
    }

    if (cfg_.mode == CutoffMode::PeriodScaled) {
      const T alpha = -std::expm1(-dt / tau_smooth_);
      log_fc_ += alpha * (log_target_ - log_fc_);
    } else {
      log_fc_ = log_target_;
    }

    refresh_timer_ += dt;
    const bool dt_changed = h_applied_ > T(0) && std::fabs(dt - h_applied_) > dt_rel_tol_ * h_applied_;
    if (force_refresh_ || dt_changed || refresh_timer_ >= refresh_period_) {
      refreshCoefficients(static_cast<double>(dt));
      refresh_timer_ = T(0);
      force_refresh_ = false;
    }

    const T u_bar = have_prev_ ? T(0.5) * (u_prev_ + a_up) : a_up;
    T dz[NMax];
    for (int i = 0; i < n_; ++i) {
      T acc = G_[i] * u_bar;
      for (int j = 0; j < n_; ++j) acc += F_[i][j] * z_[j];
      dz[i] = acc;
    }
    for (int i = 0; i < n_; ++i) z_[i] += dz[i];
    u_prev_ = a_up;
    have_prev_ = true;
    ++updates_;
    return true;
  }

  T displacement() const noexcept { return z_[0]; }
  T velocity() const noexcept { return w_applied_ * (z_[1] - c_[n_ - 1] * z_[0]); }

  T appliedCutoffHz() const noexcept { return w_applied_ / static_cast<T>(kTwoPi); }
  T smoothedCutoffHz() const noexcept { return std::exp(log_fc_); }
  int order() const noexcept { return n_; }
  int dcZeros() const noexcept { return m_; }
  bool hasWavePeriod() const noexcept { return have_tz_; }
  uint32_t rejectedSamples() const noexcept { return rejected_; }
  uint32_t acceptedSamples() const noexcept { return updates_; }
  const HeaveHPDIConfig& config() const noexcept { return cfg_; }
  const char* lastError() const noexcept { return error_; }

 private:
  static constexpr double kTwoPi = 6.28318530717958647692;


  static T clampT(T x, T lo, T hi) noexcept { return x < lo ? lo : (x > hi ? hi : x); }

  // Caches the log target cutoff and the smoothing horizon, in T, whenever T_z
  // or the configuration changes; update() then needs one expm1 per sample.
  void updateTargets() noexcept {
    const T tz = have_tz_ ? tz_ : fallback_period_;
    const T f = cfg_.mode == CutoffMode::PeriodScaled ? cutoff_ratio_ / tz : cutoff_hz_;
    log_target_ = std::log(clampT(f, f_min_, f_max_));
    tau_smooth_ = clampT(smoothing_periods_ * tz, smoothing_min_, smoothing_max_);
  }

  // Ac (normalized observer companion): Ac[i][0] = -c_{n-1-i}, Ac[i][i+1] = 1.
  void buildAc(T (&A)[NMax][NMax]) const noexcept {
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) A[i][j] = T(0);
      A[i][0] = -c_[n_ - 1 - i];
      if (i + 1 < n_) A[i][i + 1] = T(1);
    }
  }

  void refreshCoefficients(double h) noexcept {
    double w = kTwoPi * std::exp(static_cast<double>(log_fc_));
    if (w * h > cfg_.max_wh) w = cfg_.max_wh / h;
    const T wT = static_cast<T>(w);

    if (w_applied_ > T(0) && wT != w_applied_) {
      const T ratio = w_applied_ / wT;
      for (int i = 1; i < n_; ++i) z_[i] *= ratio;
    }

    // Solve (I - M/2) [F | G] = [M | beta h / w], M = h w Ac.
    T Ac[NMax][NMax];
    buildAc(Ac);
    const T hw = static_cast<T>(h * w);
    const T gscale = static_cast<T>(h / w);
    T L[NMax][NMax];
    T X[NMax][NMax + 1];
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        const T Mij = hw * Ac[i][j];
        L[i][j] = (i == j ? T(1) : T(0)) - T(0.5) * Mij;
        X[i][j] = Mij;
      }
      X[i][n_] = beta_[i] * gscale;
    }
    solveInPlace(L, X, n_ + 1);
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) F_[i][j] = X[i][j];
      G_[i] = X[i][n_];
    }
    w_applied_ = wT;
    h_applied_ = static_cast<T>(h);
  }

  // Steady state for constant input u: Ac z = -beta u / w^2.
  void primeSteadyState(T u) noexcept {
    const double w = kTwoPi * std::exp(static_cast<double>(log_fc_));
    T Ac[NMax][NMax];
    buildAc(Ac);
    T X[NMax][NMax + 1];
    for (int i = 0; i < n_; ++i) X[i][0] = -beta_[i] * static_cast<T>(static_cast<double>(u) / (w * w));
    solveInPlace(Ac, X, 1);
    for (int i = 0; i < n_; ++i) z_[i] = X[i][0];
    // Express the primed state at the cutoff the first refresh will apply.
    w_applied_ = static_cast<T>(w);
  }

  // Gaussian elimination with partial pivoting; A (n_ x n_) is destroyed and
  // the first `cols` columns of X are overwritten with A^-1 X.
  void solveInPlace(T (&A)[NMax][NMax], T (&X)[NMax][NMax + 1], int cols) const noexcept {
    const int n = n_;
    for (int k = 0; k < n; ++k) {
      int p = k;
      T best = std::fabs(A[k][k]);
      for (int r = k + 1; r < n; ++r) {
        const T v = std::fabs(A[r][k]);
        if (v > best) {
          best = v;
          p = r;
        }
      }
      if (p != k) {
        for (int j = 0; j < n; ++j) {
          const T t = A[k][j];
          A[k][j] = A[p][j];
          A[p][j] = t;
        }
        for (int j = 0; j < cols; ++j) {
          const T t = X[k][j];
          X[k][j] = X[p][j];
          X[p][j] = t;
        }
      }
      const T inv = T(1) / A[k][k];  // nonsingular: I - M/2 near identity; Ac has c_0 = 1
      for (int r = k + 1; r < n; ++r) {
        const T f = A[r][k] * inv;
        if (f == T(0)) continue;
        for (int j = k; j < n; ++j) A[r][j] -= f * A[k][j];
        for (int j = 0; j < cols; ++j) X[r][j] -= f * X[k][j];
      }
    }
    for (int k = n - 1; k >= 0; --k) {
      const T inv = T(1) / A[k][k];
      for (int j = 0; j < cols; ++j) {
        T acc = X[k][j];
        for (int i = k + 1; i < n; ++i) acc -= A[k][i] * X[i][j];
        X[k][j] = acc * inv;
      }
    }
  }

  HeaveHPDIConfig cfg_{};
  int n_ = 4;
  int m_ = 3;
  T c_[NMax + 1]{};
  T beta_[NMax]{};
  T z_[NMax]{};
  T F_[NMax][NMax]{};
  T G_[NMax]{};
  T w_applied_ = T(0);
  T h_applied_ = T(0);
  T log_fc_ = T(0);
  T log_target_ = T(0);
  T cutoff_hz_ = T(0.04);
  T cutoff_ratio_ = T(0.15);
  T fallback_period_ = T(6);
  T f_min_ = T(0.005);
  T f_max_ = T(0.5);
  T smoothing_periods_ = T(0.75);
  T smoothing_min_ = T(0.05);
  T smoothing_max_ = T(35);
  T tau_smooth_ = T(1);
  T dt_min_ = T(1e-5);
  T dt_max_ = T(0.1);
  T dt_rel_tol_ = T(1e-4);
  T refresh_period_ = T(0.1);
  T tz_ = T(6);
  T u_prev_ = T(0);
  T refresh_timer_ = T(0);
  bool have_tz_ = false;
  bool have_prev_ = false;
  bool force_refresh_ = true;
  uint32_t rejected_ = 0;
  uint32_t updates_ = 0;
  const char* error_ = nullptr;
};

}  // namespace hpdi
