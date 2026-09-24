#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Full-matrix accelerometer calibration from qualified static holds.

  Model (body frame, m/s^2):

    a_cal = S * (a_raw - b_ref - k * (T - T_ref)),   S symmetric positive definite

  S keeps all six independent entries (three scales and three symmetric
  cross-axis terms), so the result is the polar-SPD factor of the sensor's
  inverse gain matrix: the same convention as AccelCalibrator::PolarSPD, with
  no frame rotation introduced.

  The residual of one observation is the geometric gravity-norm error
  ||a_cal|| - g. Only measured vectors and the gravity magnitude enter the fit;
  pose labels never assign gravity components and repeated labels are never
  subtracted from each other. Thermal slopes are fitted jointly through the same
  measurement model and are kept only when the data determine them after the
  nine static parameters are accounted for.

  The solver runs in double precision on normalized variables (vectors in g,
  temperatures in units of T_scale) with a log-Cholesky parameterization of S,
  so every iterate is SPD. Stored coefficients are float and are re-validated
  in float after every post-fit transformation (validateFloat()).

  All buffers are fixed size (MAXO observations, MAXH holds); fit() does not
  allocate.
*/

#include "imu_calibrate/CalibrateIMU.h"

#include <algorithm>

namespace imu_cal {

// One retained observation: the mean of one qualified time block of a static
// hold. Blocks of the same hold are correlated; the fitter balances weights so
// that each hold counts once, however many blocks it produced.
enum class AccelObsRole : uint8_t { FIT = 0, VERIFY = 1 };

struct AccelObs {
  float a[3];       // mean raw specific force, body frame, m/s^2
  float w[3];       // mean raw angular rate, rad/s (diagnostic)
  float tempC;      // mean valid temperature, degC (NaN when none was valid)
  float a_std;      // RMS deviation of samples about the block mean, m/s^2
  uint32_t t_ms;    // block mid time on the session clock
  uint16_t n;       // samples averaged
  uint8_t hold;     // hold index (0..MAXH-1)
  uint8_t role;     // AccelObsRole
};

enum class AccelThermal : uint8_t {
  LEGACY = 0,     // pre-v3 coefficients: slope provenance unknown
  UNLEARNED = 1,  // k = 0: thermal bias not characterised
  LEARNED = 2,    // k fitted in this session over [k_temp_lo, k_temp_hi]
  PRESERVED = 3,  // k carried from an earlier validated session of this sensor
};

enum class AccelThermalReason : uint8_t {
  NONE = 0,
  TEMP_MISSING = 1,       // a fit block had no valid temperature
  SPAN_SMALL = 2,         // hold temperatures span too little
  INFO_LOW = 3,           // slope not determined after the static parameters
  SLOPE_IMPLAUSIBLE = 4,  // fitted slope outside physical bounds
  FIT_FAILED = 5,         // thermal solve did not converge
  NOT_ATTEMPTED = 6,      // preliminary (static-only) fit
};

inline const char* accelThermalStr(AccelThermal t) {
  switch (t) {
    case AccelThermal::LEGACY: return "LEGACY";
    case AccelThermal::UNLEARNED: return "UNLEARNED";
    case AccelThermal::LEARNED: return "LEARNED";
    case AccelThermal::PRESERVED: return "PRESERVED";
    default: return "UNKNOWN";
  }
}

inline const char* accelThermalReasonStr(AccelThermalReason r) {
  switch (r) {
    case AccelThermalReason::NONE: return "NONE";
    case AccelThermalReason::TEMP_MISSING: return "TEMP_MISSING";
    case AccelThermalReason::SPAN_SMALL: return "SPAN_SMALL";
    case AccelThermalReason::INFO_LOW: return "INFO_LOW";
    case AccelThermalReason::SLOPE_IMPLAUSIBLE: return "SLOPE_IMPLAUSIBLE";
    case AccelThermalReason::FIT_FAILED: return "FIT_FAILED";
    case AccelThermalReason::NOT_ATTEMPTED: return "NOT_ATTEMPTED";
    default: return "UNKNOWN";
  }
}

// Which acceptance gate decided a failure.
enum class AccelGate : uint8_t {
  NONE = 0, BAD_INPUT, COVERAGE, SOLVE, INFO_BIAS, INFO_CROSS, DIAG, COND, OFFDIAG,
  BIAS_MAG, RESID, CROSSVAL, VERIFY, FLOAT,
};

inline const char* accelGateStr(AccelGate g) {
  switch (g) {
    case AccelGate::NONE: return "NONE";
    case AccelGate::BAD_INPUT: return "BAD_INPUT";
    case AccelGate::COVERAGE: return "COVERAGE";
    case AccelGate::SOLVE: return "SOLVE";
    case AccelGate::INFO_BIAS: return "INFO_BIAS";
    case AccelGate::INFO_CROSS: return "INFO_CROSS";
    case AccelGate::DIAG: return "DIAG";
    case AccelGate::COND: return "COND";
    case AccelGate::OFFDIAG: return "OFFDIAG";
    case AccelGate::BIAS_MAG: return "BIAS_MAG";
    case AccelGate::RESID: return "RESID";
    case AccelGate::CROSSVAL: return "CROSSVAL";
    case AccelGate::VERIFY: return "VERIFY";
    case AccelGate::FLOAT: return "FLOAT";
    default: return "UNKNOWN";
  }
}

struct AccelFitCfg {
  double g = 9.80665;
  double T_scale = 10.0;          // degC per normalized temperature unit

  // Solver
  int max_outer = 6;              // IRLS passes
  int max_lm = 40;                // LM iterations per pass
  int max_lm_cv = 15;             // LM iterations per leave-one-hold-out refit
  double huber_c = 2.0;           // Huber knee in robust sigmas
  double block_sigma_floor = 0.002;   // m/s^2 floor of the robust block scale
  double reject_sigma = 8.0;      // blocks beyond max(reject_sigma*s, reject_abs) get weight 0
  double reject_abs = 0.08;       // m/s^2

  // Hold-level observation error: split-half scatter (data) plus this floor for
  // effects a static hold cannot reveal (nonlinearity, local g, slow drift).
  double model_floor = 0.003;     // m/s^2

  // Plausibility gates (wizard PolarSPD configuration)
  double diag_lo = 0.80, diag_hi = 1.25;
  double max_cond = 6.0;
  double max_offdiag_rms = 0.10;
  double max_bias = 1.5;          // m/s^2 per axis

  // Information gates (data-derived; no priors or regularization)
  double max_bias_sigma = 0.015;  // m/s^2 per axis, 1 sigma
  double max_cross_sigma = 0.004; // dimensionless, 1 sigma per cross term
  double max_hold_rms = 0.05;     // m/s^2 RMS of hold-mean norm errors

  // Thermal slope qualification
  double min_temp_span = 2.0;     // degC across hold-mean temperatures
  double max_k_sigma = 0.0015;    // m/s^2/degC per axis after the static parameters
  double max_k_abs = 0.02;        // m/s^2/degC per axis (~2 mg/K)
  double temp_extrap_margin = 5.0;// degC beyond the evidence range before clamping

  // Held-out checks. A left-out hold is predicted by a refit without it and
  // compared in units of its predicted uncertainty (the design gives tilted
  // holds high leverage, so raw errors alone would mostly measure geometry).
  double max_cv_t = 6.0;          // studentized leave-one-hold-out error
  double max_cv = 0.15;           // m/s^2, gross bound on the raw held-out error
  double max_verify = 0.08;       // m/s^2, verification-only recheck blocks (repeatability)
  double cv_k_sigma_mult = 3.0;   // a held-out refit is scored only if still determined

  // Float re-validation
  double float_tol = 2e-4;        // m/s^2 between double and float hold RMS
};

// A slope from an earlier validated session of the same sensor/configuration.
struct AccelThermalPrior {
  bool valid = false;
  double k[3] = {0, 0, 0};        // m/s^2/degC
  double k_temp_lo = 0, k_temp_hi = 0;
  double clamp_lo = -1000, clamp_hi = 1000;
};

struct AccelFullFitResult {
  bool ok = false;
  FitFail reason = FitFail::BAD_ARG;
  AccelGate gate = AccelGate::NONE;

  // Coefficients (double; see validateFloat() for the stored float set)
  Eigen::Matrix3d S = Eigen::Matrix3d::Identity();
  Eigen::Vector3d b = Eigen::Vector3d::Zero();   // b_ref, m/s^2
  Eigen::Vector3d k = Eigen::Vector3d::Zero();   // m/s^2/degC
  double T_ref = 25.0;
  double clamp_lo = -1000.0, clamp_hi = 1000.0;  // runtime temperature clamp

  // Thermal qualification
  AccelThermal thermal = AccelThermal::UNLEARNED;
  AccelThermalReason thermal_reason = AccelThermalReason::NONE;
  double k_temp_lo = 0, k_temp_hi = 0;           // evidence range of k
  double cal_temp_lo = 0, cal_temp_hi = 0;       // hold-mean temperatures seen
  bool temps_valid = false;
  Eigen::Vector3d k_sigma = Eigen::Vector3d::Constant(-1);  // this session's thermal fit

  // Data-derived uncertainty (1 sigma), one observation per hold
  double sigma_obs = 0;                          // m/s^2
  Eigen::Vector3d bias_sigma = Eigen::Vector3d::Constant(-1);   // m/s^2
  Eigen::Vector3d diag_sigma = Eigen::Vector3d::Constant(-1);
  Eigen::Vector3d cross_sigma = Eigen::Vector3d::Constant(-1);  // xy, xz, yz
  int weak_param = -1;           // 0..2 bias x/y/z, 3..5 cross xy/xz/yz
  double weak_ratio = 0;         // sigma/gate of weak_param
  Eigen::Matrix<double,9,9> C_unit = Eigen::Matrix<double,9,9>::Zero(); // static cov / sigma_obs^2

  // Residuals and bookkeeping
  double block_sigma = 0;        // robust block residual scale, m/s^2
  double hold_rms = 0;           // RMS of hold-mean norm errors (fit holds), m/s^2
  int n_fit_blocks = 0, n_rejected = 0, n_fit_holds = 0, n_verify_blocks = 0;

  // Held-out checks
  double cv_rms = 0, cv_max = 0, cv_tmax = 0; int n_cv = 0, n_cv_skipped = 0;
  int cv_worst_hold = -1;        // hold with the largest studentized error
  double ver_rms = 0, ver_max = 0; int n_ver = 0;
  int ver_worst_hold = -1;

  // Float re-validation of the stored coefficient set: unweighted hold-mean
  // norm RMS from the double solution and from the float coefficients.
  double ref_hold_rms = -1;
  double float_hold_rms = -1;
  bool float_ok = false;
};

namespace detail {

// Parameter vector (normalized): [l0..l5 | beta(3) | kappa(3)]
//   L = [[e^l0, 0, 0], [l1, e^l2, 0], [l3, l4, e^l5]],  S = L L^T
//   beta = b/g,  kappa = k*T_scale/g
using Theta = Eigen::Matrix<double,12,1>;

static inline Eigen::Matrix3d chol_from_theta(const Theta& th) {
  Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
  L(0,0) = std::exp(th(0));
  L(1,0) = th(1); L(1,1) = std::exp(th(2));
  L(2,0) = th(3); L(2,1) = th(4); L(2,2) = std::exp(th(5));
  return L;
}

// (row, col) of the L entry driven by each log-Cholesky parameter.
static constexpr int kLRow[6] = {0, 1, 1, 2, 2, 2};
static constexpr int kLCol[6] = {0, 0, 1, 0, 1, 2};
static constexpr bool kLDiag[6] = {true, false, true, false, false, true};

}  // namespace detail

template <int MAXO, int MAXH = 24>
class AccelFullFitter {
public:
  static_assert(MAXO > 0 && MAXO <= 1024, "MAXO unreasonable");
  static_assert(MAXH > 0 && MAXH <= 255, "MAXH unreasonable");

  using Vec3 = Eigen::Vector3d;
  using Mat3 = Eigen::Matrix3d;
  using Theta = detail::Theta;
  using Mat12 = Eigen::Matrix<double,12,12>;
  using Vec12 = Eigen::Matrix<double,12,1>;

  enum Mode : uint8_t { ZERO = 0, FIXED = 1, FREE = 2 };

  // Full fit: static geometry, thermal qualification, held-out checks.
  // When `preliminary` is set only the static model is fitted and only the
  // geometry/plausibility gates run; the wizard uses it to decide targeted
  // extra holds before leaving the accelerometer stage.
  bool fit(const AccelObs* obs, int n, const AccelFitCfg& cfg,
           const AccelThermalPrior& prior, AccelFullFitResult& out,
           bool preliminary = false)
  {
    out = AccelFullFitResult{};
    out.ok = false;
    if (!obs || n <= 0 || n > MAXO) return fail_(out, FitFail::BAD_ARG, AccelGate::BAD_INPUT);
    obs_ = obs; n_ = n; cfg_ = &cfg;

    if (!index_holds_()) return fail_(out, FitFail::NON_FINITE_INPUT, AccelGate::BAD_INPUT);
    out.n_fit_holds = n_fit_holds_;
    out.n_fit_blocks = n_fit_blocks_;
    out.n_verify_blocks = n_verify_blocks_;
    if (n_fit_holds_ < 9 || n_fit_blocks_ < 12) return fail_(out, FitFail::TOO_FEW_SAMPLES, AccelGate::BAD_INPUT);

    // Existing coverage/planarity protection, applied to the fit blocks.
    {
      int m = 0;
      for (int i = 0; i < n_; ++i) {
        if (obs_[i].role != (uint8_t)AccelObsRole::FIT) continue;
        cov_[m++] = Eigen::Vector3f(obs_[i].a[0], obs_[i].a[1], obs_[i].a[2]);
      }
      if (!degeneracy_check_coverage3<float>(cov_, m, (float)cfg.g, 0.90f, 0.08f, 1e-3f))
        return fail_(out, FitFail::DEGENERATE_COVERAGE, AccelGate::COVERAGE);
    }

    temperature_evidence_(out);
    T_ref_ = out.T_ref;
    clamp_lo_ = -1000.0;
    clamp_hi_ = 1000.0;

    // 1) Static model (k = 0) from the identity; the gravity scale is the only
    //    initialization input.
    Theta th0 = Theta::Zero();
    {
      double sn = 0; int m = 0;
      for (int i = 0; i < n_; ++i) {
        if (base_w_[i] <= 0) continue;
        sn += Vec3(obs_[i].a[0], obs_[i].a[1], obs_[i].a[2]).norm() / cfg.g; ++m;
      }
      const double s0 = (m > 0 && sn > 0) ? (double)m / sn : 1.0;
      th0(0) = th0(2) = th0(5) = 0.5 * std::log(s0);
    }
    for (int i = 0; i < n_; ++i) rob_w_[i] = 1.0f;
    if (!solve_(ZERO, Vec3::Zero(), -1, true, th0)) return fail_(out, FitFail::SOLVE_LDLT_FAIL, AccelGate::SOLVE);

    Theta th = th0;
    Mode mode = ZERO;
    Vec3 kfix = Vec3::Zero();

    // 2) Thermal slope: this session's evidence first, then a compatible prior.
    if (preliminary) {
      out.thermal = AccelThermal::UNLEARNED;
      out.thermal_reason = AccelThermalReason::NOT_ATTEMPTED;
    } else {
      bool learned = false;
      if (!out.temps_valid) {
        out.thermal_reason = AccelThermalReason::TEMP_MISSING;
      } else if (out.cal_temp_hi - out.cal_temp_lo < cfg.min_temp_span) {
        out.thermal_reason = AccelThermalReason::SPAN_SMALL;
      } else {
        Theta th1 = th0;
        save_rob_();
        if (!solve_(FREE, Vec3::Zero(), -1, true, th1)) {
          restore_rob_();
          out.thermal_reason = AccelThermalReason::FIT_FAILED;
        } else {
          Mat12 H; double sig = 0;
          info_(th1, FREE, Vec3::Zero(), -1, H);
          sig = sigma_obs_(th1, FREE, Vec3::Zero());
          Vec3 ks; bool det = thermal_sigma_(H, sig, ks);
          const Vec3 k1 = th1.segment<3>(9) * (cfg.g / cfg.T_scale);
          out.k_sigma = ks;
          if (!det || ks.maxCoeff() > cfg.max_k_sigma) {
            restore_rob_();
            out.thermal_reason = AccelThermalReason::INFO_LOW;
          } else if (k1.cwiseAbs().maxCoeff() > cfg.max_k_abs) {
            restore_rob_();
            out.thermal_reason = AccelThermalReason::SLOPE_IMPLAUSIBLE;
          } else {
            learned = true;
            th = th1; mode = FREE;
            out.thermal = AccelThermal::LEARNED;
            out.thermal_reason = AccelThermalReason::NONE;
            out.k_temp_lo = out.cal_temp_lo;
            out.k_temp_hi = out.cal_temp_hi;
            out.clamp_lo = out.cal_temp_lo - cfg.temp_extrap_margin;
            out.clamp_hi = out.cal_temp_hi + cfg.temp_extrap_margin;
          }
        }
      }
      if (!learned) {
        if (prior.valid && any_valid_temp_) {
          kfix = Vec3(prior.k[0], prior.k[1], prior.k[2]);
          clamp_lo_ = prior.clamp_lo;
          clamp_hi_ = prior.clamp_hi;
          Theta th2 = th0;
          if (!solve_(FIXED, kfix, -1, true, th2)) return fail_(out, FitFail::SOLVE_LDLT_FAIL, AccelGate::SOLVE);
          th = th2; mode = FIXED;
          out.thermal = AccelThermal::PRESERVED;
          out.k_temp_lo = prior.k_temp_lo;
          out.k_temp_hi = prior.k_temp_hi;
          out.clamp_lo = prior.clamp_lo;
          out.clamp_hi = prior.clamp_hi;
        } else {
          out.thermal = AccelThermal::UNLEARNED;
          out.clamp_lo = -1000.0;
          out.clamp_hi = 1000.0;
        }
      }
    }

    // 3) Coefficients in physical units.
    {
      const Mat3 L = detail::chol_from_theta(th);
      Mat3 S = L * L.transpose();
      out.S = 0.5 * (S + S.transpose());
      out.b = th.segment<3>(6) * cfg.g;
      if (mode == FREE) out.k = th.segment<3>(9) * (cfg.g / cfg.T_scale);
      else if (mode == FIXED) out.k = kfix;
      else out.k.setZero();
      if (!out.S.allFinite() || !out.b.allFinite() || !out.k.allFinite())
        return fail_(out, FitFail::NON_FINITE_INPUT, AccelGate::SOLVE);
    }

    // 4) Residual statistics and data-derived uncertainty.
    residual_stats_(th, mode, kfix, out);
    {
      Mat12 H;
      info_(th, mode, kfix, -1, H);
      out.sigma_obs = sigma_obs_(th, mode, kfix) * cfg.g;
      const int P = nparams_(mode);
      Mat12 C;
      if (!info_inverse_(H, P, C)) {
        out.weak_param = weakest_null_(H, P);
        out.weak_ratio = 1e9;
        return fail_(out, FitFail::ACCEL_INFO_LOW, out.weak_param < 3 ? AccelGate::INFO_BIAS : AccelGate::INFO_CROSS);
      }
      const double so = out.sigma_obs / cfg.g;  // normalized
      for (int j = 0; j < 3; ++j) {
        out.diag_sigma(j) = so * std::sqrt(std::max(0.0, C(j, j)));
        out.cross_sigma(j) = so * std::sqrt(std::max(0.0, C(3 + j, 3 + j)));
        out.bias_sigma(j) = so * std::sqrt(std::max(0.0, C(6 + j, 6 + j))) * cfg.g;
      }
      out.C_unit = C.topLeftCorner(9, 9);
      out.weak_ratio = 0; out.weak_param = -1;
      for (int j = 0; j < 3; ++j) {
        const double rb = out.bias_sigma(j) / cfg.max_bias_sigma;
        const double rc = out.cross_sigma(j) / cfg.max_cross_sigma;
        if (rb > out.weak_ratio) { out.weak_ratio = rb; out.weak_param = j; }
        if (rc > out.weak_ratio) { out.weak_ratio = rc; out.weak_param = 3 + j; }
      }
    }

    // 5) Plausibility gates (unchanged thresholds from the PolarSPD wizard path).
    {
      const Mat3& S = out.S;
      for (int j = 0; j < 3; ++j) {
        if (!(S(j, j) >= cfg.diag_lo && S(j, j) <= cfg.diag_hi))
          return fail_(out, FitFail::ACCEL_S_UNPHYSICAL, AccelGate::DIAG);
      }
      float condv = 0;
      if (!cond_spd_3x3<float>(S.cast<float>(), condv) || condv > (float)cfg.max_cond)
        return fail_(out, FitFail::ACCEL_S_UNPHYSICAL, AccelGate::COND);
      if (offdiag_rms_3x3<double>(S) > cfg.max_offdiag_rms)
        return fail_(out, FitFail::ACCEL_S_UNPHYSICAL, AccelGate::OFFDIAG);
      if (out.b.cwiseAbs().maxCoeff() > cfg.max_bias)
        return fail_(out, FitFail::ACCEL_S_UNPHYSICAL, AccelGate::BIAS_MAG);
    }

    // 6) Information gates: the measured geometry must determine every static
    //    parameter to the stated precision.
    if (out.bias_sigma.maxCoeff() > cfg.max_bias_sigma)
      return fail_(out, FitFail::ACCEL_INFO_LOW, AccelGate::INFO_BIAS);
    if (out.cross_sigma.maxCoeff() > cfg.max_cross_sigma)
      return fail_(out, FitFail::ACCEL_INFO_LOW, AccelGate::INFO_CROSS);
    if (out.hold_rms > cfg.max_hold_rms)
      return fail_(out, FitFail::INLIERS_TOO_FEW, AccelGate::RESID);

    // 7) Leave-one-hold-out prediction on holds whose removal leaves the model
    //    determined (also in the preliminary fit, so an inconsistent main hold
    //    is recaptured before leaving the accelerometer stage).
    crossval_(th, mode, kfix, out);
    if (out.n_cv > 0 && (out.cv_tmax > cfg.max_cv_t || out.cv_max > cfg.max_cv))
      return fail_(out, FitFail::ACCEL_CROSSCHECK_FAIL, AccelGate::CROSSVAL);

    if (preliminary) {
      out.ok = true;
      out.reason = FitFail::OK;
      return true;
    }

    out.ok = true;
    out.reason = FitFail::OK;

    // 8) Stored float coefficients must reproduce the double fit, and the
    //    verification blocks are scored with exactly those coefficients.
    out.ref_hold_rms = holdRms(obs_, n_, cfg, [&](const Eigen::Vector3f& a, float tC) {
      const Vec3 ad = a.cast<double>();
      double t = tC;
      Vec3 bb = out.b;
      if (std::isfinite(t)) bb += out.k * (std::min(std::max(t, out.clamp_lo), out.clamp_hi) - out.T_ref);
      return Vec3(out.S * (ad - bb));
    });
    AccelCalibration<float> fc;
    toFloat(out, cfg.g, fc);
    if (!validateFloat(obs_, n_, fc, cfg, out)) {
      out.ok = false;
      return fail_(out, FitFail::ACCEL_FLOAT_CHECK_FAIL, AccelGate::FLOAT);
    }
    if (out.n_ver > 0 && out.ver_max > cfg.max_verify) {
      out.ok = false;
      return fail_(out, FitFail::ACCEL_CROSSCHECK_FAIL, AccelGate::VERIFY);
    }
    return true;
  }

  // The float runtime object the stored blob will reproduce.
  static void toFloat(const AccelFullFitResult& r, double g, AccelCalibration<float>& fc) {
    fc = AccelCalibration<float>{};
    Eigen::Matrix3d S = 0.5 * (r.S + r.S.transpose());
    fc.S = S.cast<float>();
    fc.g = (float)g;
    fc.biasT.ok = true;
    fc.biasT.T0 = (float)r.T_ref;
    fc.biasT.b0 = r.b.cast<float>();
    fc.biasT.k = r.k.cast<float>();
    fc.biasT.T_lo = (float)r.clamp_lo;
    fc.biasT.T_hi = (float)r.clamp_hi;
    fc.rms_mag = (float)r.hold_rms;
    fc.ok = r.ok;
  }

  // Mean gravity-norm error of each hold (unweighted over its blocks) and the
  // RMS over holds, for any calibration functor a_cal = f(a_raw, tempC).
  template <typename ApplyFn>
  static double holdRms(const AccelObs* obs, int n, const AccelFitCfg& cfg, ApplyFn apply,
                        AccelObsRole role = AccelObsRole::FIT, double* max_abs = nullptr,
                        int* n_holds = nullptr, int* worst_hold = nullptr)
  {
    double hs[MAXH], hw[MAXH];
    for (int h = 0; h < MAXH; ++h) { hs[h] = hw[h] = 0; }
    for (int i = 0; i < n; ++i) {
      const AccelObs& o = obs[i];
      if (o.hold >= MAXH || o.role != (uint8_t)role) continue;
      const auto ac = apply(Eigen::Vector3f(o.a[0], o.a[1], o.a[2]), o.tempC);
      const double e = (double)ac.norm() - cfg.g;
      if (!std::isfinite(e)) return NAN;
      hs[o.hold] += e; hw[o.hold] += 1;
    }
    double sse = 0, mx = 0; int nh = 0, worst = -1;
    for (int h = 0; h < MAXH; ++h) {
      if (hw[h] <= 0) continue;
      const double m = hs[h] / hw[h];
      sse += m * m; ++nh;
      if (std::fabs(m) >= mx) { mx = std::fabs(m); worst = h; }
    }
    if (worst_hold) *worst_hold = worst;
    if (max_abs) *max_abs = mx;
    if (n_holds) *n_holds = nh;
    return nh ? std::sqrt(sse / nh) : 0.0;
  }

  // Re-validates a float coefficient set exactly as the runtime applies it:
  // finiteness, exact symmetry, SPD, the plausibility gates, and agreement of
  // the hold-level gravity-norm RMS with the double solution (r.ref_hold_rms).
  // Also scores the verification-only blocks with these coefficients.
  // Call again after every serialization round trip.
  static bool validateFloat(const AccelObs* obs, int n, const AccelCalibration<float>& fc,
                            const AccelFitCfg& cfg, AccelFullFitResult& r)
  {
    r.float_ok = false;
    if (!fc.S.allFinite() || !fc.biasT.b0.allFinite() || !fc.biasT.k.allFinite() ||
        !finiteT(fc.biasT.T0) || !finiteT(fc.biasT.T_lo) || !finiteT(fc.biasT.T_hi)) return false;
    if (fc.S(0,1) != fc.S(1,0) || fc.S(0,2) != fc.S(2,0) || fc.S(1,2) != fc.S(2,1)) return false;
    Eigen::LLT<Eigen::Matrix3f> llt(fc.S);
    if (llt.info() != Eigen::Success) return false;
    for (int j = 0; j < 3; ++j)
      if (!(fc.S(j,j) >= (float)cfg.diag_lo && fc.S(j,j) <= (float)cfg.diag_hi)) return false;
    float condv = 0;
    if (!cond_spd_3x3<float>(fc.S, condv) || condv > (float)cfg.max_cond) return false;
    if (offdiag_rms_3x3<float>(fc.S) > (float)cfg.max_offdiag_rms) return false;
    if (fc.biasT.b0.cwiseAbs().maxCoeff() > (float)cfg.max_bias) return false;

    auto apply = [&](const Eigen::Vector3f& a, float tC) { return Eigen::Vector3f(fc.apply(a, tC)); };
    r.float_hold_rms = holdRms(obs, n, cfg, apply);
    if (!std::isfinite(r.float_hold_rms)) return false;
    if (!(std::fabs(r.float_hold_rms - r.ref_hold_rms) <= cfg.float_tol)) return false;
    r.ver_rms = holdRms(obs, n, cfg, apply, AccelObsRole::VERIFY, &r.ver_max, &r.n_ver, &r.ver_worst_hold);
    if (!std::isfinite(r.ver_rms)) return false;
    r.float_ok = true;
    return true;
  }

  // Predicted variance reduction of static parameter `param` (0..2 bias,
  // 3..5 cross) from one more hold whose measured up-direction is `up_unit`.
  static double infoGain(const AccelFullFitResult& r, const Vec3& up_unit, int param) {
    if (param < 0 || param > 5) return 0;
    const Eigen::Matrix<double,9,1> f = staticFeature_(up_unit);
    const Eigen::Matrix<double,9,1> Cf = r.C_unit * f;
    const double den = 1.0 + f.dot(Cf);
    const int idx = (param < 3) ? 6 + param : param;
    return (den > 0) ? Cf(idx) * Cf(idx) / den : 0;
  }

private:
  // Linearized static regressor of a unit up-direction at S = I, b = 0 (in the
  // physical parameter order used by the information matrix).
  static Eigen::Matrix<double,9,1> staticFeature_(const Vec3& u) {
    Eigen::Matrix<double,9,1> f;
    f << u.x()*u.x(), u.y()*u.y(), u.z()*u.z(),
         2*u.x()*u.y(), 2*u.x()*u.z(), 2*u.y()*u.z(),
         -u.x(), -u.y(), -u.z();
    return f;
  }

  bool fail_(AccelFullFitResult& out, FitFail f, AccelGate g) {
    out.ok = false;
    out.reason = f;
    out.gate = g;
    return false;
  }

  bool index_holds_() {
    n_fit_holds_ = 0; n_fit_blocks_ = 0; n_verify_blocks_ = 0;
    for (int h = 0; h < MAXH; ++h) { hold_n_[h] = 0; }
    for (int i = 0; i < n_; ++i) {
      const AccelObs& o = obs_[i];
      if (o.hold >= MAXH) return false;
      if (!(std::isfinite(o.a[0]) && std::isfinite(o.a[1]) && std::isfinite(o.a[2]))) return false;
      if (o.role == (uint8_t)AccelObsRole::FIT) { hold_n_[o.hold]++; ++n_fit_blocks_; }
      else ++n_verify_blocks_;
    }
    for (int h = 0; h < MAXH; ++h) if (hold_n_[h] > 0) ++n_fit_holds_;
    for (int i = 0; i < n_; ++i) {
      const AccelObs& o = obs_[i];
      base_w_[i] = (o.role == (uint8_t)AccelObsRole::FIT && hold_n_[o.hold] > 0)
                   ? 1.0f / (float)hold_n_[o.hold] : 0.0f;
    }
    return true;
  }

  void temperature_evidence_(AccelFullFitResult& out) {
    double hs[MAXH]; int hc[MAXH];
    for (int h = 0; h < MAXH; ++h) { hs[h] = 0; hc[h] = 0; }
    bool all = true; any_valid_temp_ = false;
    for (int i = 0; i < n_; ++i) {
      if (base_w_[i] <= 0) continue;
      const float t = obs_[i].tempC;
      if (!std::isfinite(t)) { all = false; continue; }
      hs[obs_[i].hold] += t; hc[obs_[i].hold]++;
      any_valid_temp_ = true;
    }
    double lo = 1e30, hi = -1e30, sum = 0; int nh = 0;
    for (int h = 0; h < MAXH; ++h) {
      if (hc[h] == 0) continue;
      const double m = hs[h] / hc[h];
      lo = std::min(lo, m); hi = std::max(hi, m); sum += m; ++nh;
    }
    out.temps_valid = all && nh > 0;
    if (nh > 0) {
      // Round to 0.1 degC so the stored reference is readable and exact in float.
      out.T_ref = std::round(10.0 * sum / nh) / 10.0;
      out.cal_temp_lo = lo; out.cal_temp_hi = hi;
    } else {
      out.T_ref = 25.0;
      out.cal_temp_lo = out.cal_temp_hi = 0.0;
    }
  }

  // Temperature regressor exactly as the runtime evaluates it: clamped to the
  // coefficient set's range; a missing temperature evaluates T_ref.
  double tau_(int i) const {
    const float t = obs_[i].tempC;
    if (!std::isfinite(t)) return 0.0;
    const double tc = std::min(std::max((double)t, clamp_lo_), clamp_hi_);
    return (tc - T_ref_) / cfg_->T_scale;
  }

  // Residual and Jacobian of observation i in the optimization parameters.
  // Returns r (normalized); fills J (12) when J != nullptr.
  double resid_(const Theta& th, const Mat3& L, const Mat3& S, int mode, const Vec3& kfix_n,
                int i, Vec12* J) const {
    const Vec3 x(obs_[i].a[0] / cfg_->g, obs_[i].a[1] / cfg_->g, obs_[i].a[2] / cfg_->g);
    const double tau = tau_(i);
    Vec3 kap = Vec3::Zero();
    if (mode == FREE) kap = th.segment<3>(9);
    else if (mode == FIXED) kap = kfix_n;
    const Vec3 u = x - th.segment<3>(6) - kap * tau;
    const Vec3 v = S * u;
    const double rho = v.norm();
    if (!(rho > 1e-9)) { if (J) J->setZero(); return 1e9; }
    if (J) {
      const Mat3 G = v * u.transpose() / rho;
      const Mat3 M = (G + G.transpose()) * L;
      for (int k = 0; k < 6; ++k) {
        const int p = detail::kLRow[k], q = detail::kLCol[k];
        const double e = detail::kLDiag[k] ? L(p, q) : 1.0;
        (*J)(k) = e * M(p, q);
      }
      const Vec3 jb = -(S * v) / rho;
      J->segment<3>(6) = jb;
      J->segment<3>(9) = jb * tau;
    }
    return rho - 1.0;
  }

  // Physical-parameter regressor at the solution: [s_xx s_yy s_zz s_xy s_xz s_yz | b | kappa].
  void phys_feature_(const Theta& th, const Mat3& S, int mode, const Vec3& kfix_n, int i, Vec12& f) const {
    const Vec3 x(obs_[i].a[0] / cfg_->g, obs_[i].a[1] / cfg_->g, obs_[i].a[2] / cfg_->g);
    const double tau = tau_(i);
    Vec3 kap = Vec3::Zero();
    if (mode == FREE) kap = th.segment<3>(9);
    else if (mode == FIXED) kap = kfix_n;
    const Vec3 u = x - th.segment<3>(6) - kap * tau;
    const Vec3 v = S * u;
    const double rho = std::max(v.norm(), 1e-9);
    const Mat3 G = v * u.transpose() / rho;
    f(0) = G(0,0); f(1) = G(1,1); f(2) = G(2,2);
    f(3) = G(0,1) + G(1,0); f(4) = G(0,2) + G(2,0); f(5) = G(1,2) + G(2,1);
    const Vec3 jb = -(S * v) / rho;
    f.segment<3>(6) = jb;
    f.segment<3>(9) = jb * tau;
  }

  static int nparams_(int mode) { return (mode == FREE) ? 12 : 9; }

  double weight_(int i, int exclude_hold) const {
    if (base_w_[i] <= 0) return 0;
    if (exclude_hold >= 0 && obs_[i].hold == exclude_hold) return 0;
    return (double)base_w_[i] * (double)rob_w_[i];
  }

  double cost_(const Theta& th, int mode, const Vec3& kfix_n, int exclude_hold) const {
    const Mat3 L = detail::chol_from_theta(th);
    const Mat3 S = L * L.transpose();
    double c = 0;
    for (int i = 0; i < n_; ++i) {
      const double w = weight_(i, exclude_hold);
      if (w <= 0) continue;
      const double r = resid_(th, L, S, mode, kfix_n, i, nullptr);
      c += w * r * r;
    }
    return c;
  }

  bool lm_(Theta& th, int mode, const Vec3& kfix_n, int exclude_hold, int max_iter) const {
    const int P = nparams_(mode);
    double lambda = 1e-4;
    double cost = cost_(th, mode, kfix_n, exclude_hold);
    if (!std::isfinite(cost)) return false;
    for (int it = 0; it < max_iter; ++it) {
      const Mat3 L = detail::chol_from_theta(th);
      const Mat3 S = L * L.transpose();
      Mat12 H = Mat12::Zero(); Vec12 gv = Vec12::Zero();
      for (int i = 0; i < n_; ++i) {
        const double w = weight_(i, exclude_hold);
        if (w <= 0) continue;
        Vec12 J;
        const double r = resid_(th, L, S, mode, kfix_n, i, &J);
        H.selfadjointView<Eigen::Lower>().rankUpdate(J, w);
        gv.noalias() += w * r * J;
      }
      H = H.selfadjointView<Eigen::Lower>();
      bool accepted = false;
      double new_cost = cost;
      Theta th_new = th;
      Vec12 delta = Vec12::Zero();
      for (int tr = 0; tr < 12; ++tr) {
        Mat12 A = H;
        for (int d = 0; d < P; ++d) A(d, d) += lambda * std::max(H(d, d), 1e-12) + 1e-18;
        pin_(A, P);
        Vec12 rhs = -gv;
        for (int k = P; k < 12; ++k) rhs(k) = 0.0;
        Eigen::LDLT<Mat12> ldlt(A);
        if (ldlt.info() != Eigen::Success) { lambda *= 10; continue; }
        delta = ldlt.solve(rhs);
        if (!delta.allFinite()) { lambda *= 10; continue; }
        for (int k = P; k < 12; ++k) delta(k) = 0.0;
        th_new = th + delta;
        new_cost = cost_(th_new, mode, kfix_n, exclude_hold);
        if (std::isfinite(new_cost) && new_cost <= cost) { accepted = true; break; }
        lambda *= 10;
      }
      if (!accepted) break;  // no descent direction left: converged
      const double dec = cost - new_cost;
      th = th_new;
      cost = new_cost;
      lambda = std::max(lambda / 10.0, 1e-12);
      // Steps below 1e-8 g (~1e-7 m/s^2) or relative cost changes below 1e-10
      // are far under any physical resolution; iterating further only spends
      // soft-float time on the ESP32-S3.
      if (delta.cwiseAbs().maxCoeff() < 1e-8 || dec <= 1e-10 * std::max(cost, 1e-30)) break;
    }
    return th.allFinite();
  }

  // Parameters k >= P are not free: an identity row/column fixes their step
  // at zero exactly (a constraint, not regularization).
  static void pin_(Mat12& A, int P) {
    for (int k = P; k < 12; ++k) { A.row(k).setZero(); A.col(k).setZero(); A(k, k) = 1.0; }
  }

  // Covariance (unit sigma) of the first P parameters from one fixed-size
  // LDLT; false when that block is not positive definite to working precision.
  static bool info_inverse_(const Mat12& H, int P, Mat12& C) {
    Mat12 A = H;
    pin_(A, P);
    const double scale = std::max(1.0, A.diagonal().maxCoeff());
    Eigen::LDLT<Mat12> ldlt(A);
    if (ldlt.info() != Eigen::Success) return false;
    const Vec12 d = ldlt.vectorD();
    if (!(d.minCoeff() > 1e-12 * scale)) return false;
    C = ldlt.solve(Mat12::Identity());
    if (!C.allFinite()) return false;
    for (int k = 0; k < P; ++k) if (!(C(k, k) > 0)) return false;
    return true;
  }

  // IRLS around LM. When update_robust is false the robust weights are frozen
  // (leave-one-hold-out refits reuse the full-data weights).
  bool solve_(int mode, const Vec3& kfix, int exclude_hold, bool update_robust, Theta& th) {
    const Vec3 kfix_n = kfix * (cfg_->T_scale / cfg_->g);
    if (mode == FIXED) th.segment<3>(9) = kfix_n;
    if (mode == ZERO) th.segment<3>(9).setZero();
    const int outer = update_robust ? cfg_->max_outer : 1;
    const int iters = update_robust ? cfg_->max_lm : cfg_->max_lm_cv;
    for (int it = 0; it < outer; ++it) {
      if (!lm_(th, mode, kfix_n, exclude_hold, iters)) return false;
      if (!update_robust) break;
      // Robust scale from the fit blocks (MAD), floored at sensor-level noise.
      const Mat3 L = detail::chol_from_theta(th);
      const Mat3 S = L * L.transpose();
      int m = 0;
      for (int i = 0; i < n_; ++i) {
        if (base_w_[i] <= 0) continue;
        r_[i] = (float)resid_(th, L, S, mode, kfix_n, i, nullptr);
        scratch_[m++] = r_[i];
      }
      if (m < 10) return false;
      const float med = nth_(scratch_, m);
      m = 0;
      for (int i = 0; i < n_; ++i) {
        if (base_w_[i] <= 0) continue;
        scratch_[m++] = (float)std::fabs(r_[i] - med);
      }
      const double mad = nth_(scratch_, m);
      const double s = std::max(1.4826 * mad, cfg_->block_sigma_floor / cfg_->g);
      sigma_block_ = s;
      const double knee = cfg_->huber_c * s;
      const double rej = std::max(cfg_->reject_sigma * s, cfg_->reject_abs / cfg_->g);
      double change = 0;
      for (int i = 0; i < n_; ++i) {
        if (base_w_[i] <= 0) continue;
        const double a = std::fabs((double)r_[i] - med);
        double w = (a <= knee) ? 1.0 : knee / a;
        if (a > rej) w = 0.0;
        change = std::max(change, std::fabs(w - (double)rob_w_[i]));
        rob_w_[i] = (float)w;
      }
      // A hold whose blocks are all rejected would silently vanish; keep the
      // least-bad block so the hold still counts (the residual gates judge it).
      for (int h = 0; h < MAXH; ++h) {
        if (hold_n_[h] == 0) continue;
        double sw = 0; int best = -1; double besta = 1e30;
        for (int i = 0; i < n_; ++i) {
          if (base_w_[i] <= 0 || obs_[i].hold != h) continue;
          sw += rob_w_[i];
          const double a = std::fabs((double)r_[i] - med);
          if (a < besta) { besta = a; best = i; }
        }
        if (sw <= 0 && best >= 0) rob_w_[best] = 1.0f;
      }
      if (change < 1e-3 && it > 0) break;
    }
    return true;
  }

  // Hold-normalized information in physical parameters (normalized units).
  void info_(const Theta& th, int mode, const Vec3& kfix, int exclude_hold, Mat12& H) const {
    const Vec3 kfix_n = kfix * (cfg_->T_scale / cfg_->g);
    const Mat3 L = detail::chol_from_theta(th);
    const Mat3 S = L * L.transpose();
    double hs[MAXH];
    for (int h = 0; h < MAXH; ++h) hs[h] = 0;
    for (int i = 0; i < n_; ++i) hs[obs_[i].hold] += weight_(i, exclude_hold);
    H.setZero();
    for (int i = 0; i < n_; ++i) {
      const double w = weight_(i, exclude_hold);
      if (w <= 0 || hs[obs_[i].hold] <= 0) continue;
      Vec12 f;
      phys_feature_(th, S, mode, kfix_n, i, f);
      H.selfadjointView<Eigen::Lower>().rankUpdate(f, w / hs[obs_[i].hold]);
    }
    H = H.selfadjointView<Eigen::Lower>();
  }

  // sigma_obs (normalized): split-half hold scatter plus the model floor.
  // Also records each hold's own split-half variance (hold_var_, normalized).
  double sigma_obs_(const Theta& th, int mode, const Vec3& kfix) {
    const Vec3 kfix_n = kfix * (cfg_->T_scale / cfg_->g);
    const Mat3 L = detail::chol_from_theta(th);
    const Mat3 S = L * L.transpose();
    double acc = 0; int nh = 0;
    for (int h = 0; h < MAXH; ++h) {
      hold_var_[h] = 0;
      if (hold_n_[h] < 4) continue;
      // Blocks are appended in time order within a hold.
      const int half = hold_n_[h] / 2;
      double s1 = 0, w1 = 0, s2 = 0, w2 = 0; int seen = 0;
      for (int i = 0; i < n_; ++i) {
        if (base_w_[i] <= 0 || obs_[i].hold != h) continue;
        const double w = rob_w_[i];
        const double r = resid_(th, L, S, mode, kfix_n, i, nullptr);
        if (seen < half) { s1 += w * r; w1 += w; } else { s2 += w * r; w2 += w; }
        ++seen;
      }
      if (w1 <= 0 || w2 <= 0) continue;
      const double d = s1 / w1 - s2 / w2;
      hold_var_[h] = 0.25 * d * d;
      acc += hold_var_[h]; ++nh;
    }
    const double split = (nh > 0) ? acc / nh : 0.0;
    const double fl = cfg_->model_floor / cfg_->g;
    return std::sqrt(split + fl * fl);
  }

  // Thermal slope sigma after the static parameters: the kappa block of the
  // full inverse equals the inverse of the Schur complement of the static block.
  bool thermal_sigma_(const Mat12& H, double sigma_n, Vec3& ks) const {
    Mat12 C;
    if (!info_inverse_(H, 12, C)) { ks.setConstant(1e9); return false; }
    const double conv = cfg_->g / cfg_->T_scale;
    for (int j = 0; j < 3; ++j) ks(j) = sigma_n * std::sqrt(std::max(0.0, C(9 + j, 9 + j))) * conv;
    return ks.allFinite();
  }

  // Bias/cross parameter least determined by a singular information block:
  // largest diagonal of a slightly damped inverse (diagnostic only).
  int weakest_null_(const Mat12& H, int P) const {
    Mat12 A = H;
    const double eps = 1e-9 * std::max(1.0, H.diagonal().head(P).sum());
    for (int k = 0; k < P; ++k) A(k, k) += eps;
    Mat12 C;
    if (!info_inverse_(A, P, C)) return 0;
    const double gb = cfg_->max_bias_sigma / cfg_->g, gc = cfg_->max_cross_sigma;
    int best = 0; double bv = -1;
    for (int j = 0; j < 3; ++j) {
      const double vb = C(6 + j, 6 + j) / (gb * gb), vc = C(3 + j, 3 + j) / (gc * gc);
      if (vb > bv) { bv = vb; best = j; }
      if (vc > bv) { bv = vc; best = 3 + j; }
    }
    return best;
  }

  void residual_stats_(const Theta& th, int mode, const Vec3& kfix, AccelFullFitResult& out) {
    const Vec3 kfix_n = kfix * (cfg_->T_scale / cfg_->g);
    const Mat3 L = detail::chol_from_theta(th);
    const Mat3 S = L * L.transpose();
    double hs[MAXH], hw[MAXH];
    for (int h = 0; h < MAXH; ++h) { hs[h] = hw[h] = 0; }
    int rej = 0;
    for (int i = 0; i < n_; ++i) {
      if (base_w_[i] <= 0) continue;
      if (rob_w_[i] <= 0) ++rej;
      const double r = resid_(th, L, S, mode, kfix_n, i, nullptr);
      const double w = rob_w_[i];
      hs[obs_[i].hold] += w * r; hw[obs_[i].hold] += w;
    }
    double sse = 0; int nh = 0;
    for (int h = 0; h < MAXH; ++h) {
      if (hw[h] <= 0) continue;
      const double m = hs[h] / hw[h];
      sse += m * m; ++nh;
    }
    out.hold_rms = nh ? std::sqrt(sse / nh) * cfg_->g : 0.0;
    out.block_sigma = sigma_block_ * cfg_->g;
    out.n_rejected = rej;
  }

  void crossval_(const Theta& th, int mode, const Vec3& kfix, AccelFullFitResult& out) {
    const Vec3 kfix_n = kfix * (cfg_->T_scale / cfg_->g);
    double sse = 0, mx = 0, tmx = 0; int ncv = 0, nskip = 0;
    for (int h = 0; h < MAXH; ++h) {
      if (hold_n_[h] == 0) continue;
      // Determinacy of the refit without hold h.
      Mat12 H;
      info_(th, mode, kfix, h, H);
      const int P = nparams_(mode);
      Mat12 C;
      bool det = info_inverse_(H, P, C);
      if (det) {
        const double so = out.sigma_obs / cfg_->g;
        for (int j = 0; j < 3; ++j) {
          if (so * std::sqrt(std::max(0.0, C(6 + j, 6 + j))) * cfg_->g > cfg_->cv_k_sigma_mult * cfg_->max_bias_sigma) det = false;
          if (so * std::sqrt(std::max(0.0, C(3 + j, 3 + j))) > cfg_->cv_k_sigma_mult * cfg_->max_cross_sigma) det = false;
        }
        if (det && mode == FREE) {
          Vec3 ks;
          if (!thermal_sigma_(H, so, ks) || ks.maxCoeff() > cfg_->cv_k_sigma_mult * cfg_->max_k_sigma) det = false;
        }
      }
      if (!det) { ++nskip; continue; }
      Theta thh = th;
      if (!lm_(thh, mode, kfix_n, h, cfg_->max_lm_cv)) { ++nskip; continue; }
      const Mat3 Lh = detail::chol_from_theta(thh);
      const Mat3 Sh = Lh * Lh.transpose();
      double s = 0, w = 0;
      Vec12 fbar = Vec12::Zero();
      for (int i = 0; i < n_; ++i) {
        if (base_w_[i] <= 0 || obs_[i].hold != h) continue;
        const double wi = rob_w_[i];
        s += wi * resid_(thh, Lh, Sh, mode, kfix_n, i, nullptr); w += wi;
        Vec12 f;
        phys_feature_(thh, Sh, mode, kfix_n, i, f);
        fbar += wi * f;
      }
      if (w <= 0) { ++nskip; continue; }
      fbar /= w;
      const double e = s / w * cfg_->g;
      // Predicted variance of the held-out hold mean: its own error plus the
      // refit's parameter uncertainty along its regressor.
      Vec12 fp = fbar;
      for (int k = P; k < 12; ++k) fp(k) = 0.0;
      const double lev = fp.dot(C * fp);
      // A hold is judged against the larger of the pooled error and its own
      // split-half scatter (hand-held holds are noisier than table holds).
      const double fl = cfg_->model_floor;
      const double own = std::sqrt(hold_var_[h] * cfg_->g * cfg_->g + fl * fl);
      const double sd = std::sqrt(std::pow(std::max(out.sigma_obs, own), 2) + out.sigma_obs * out.sigma_obs * std::max(0.0, lev));
      const double tval = (sd > 0) ? std::fabs(e) / sd : 0.0;
      sse += e * e; mx = std::max(mx, std::fabs(e)); ++ncv;
      if (tval > tmx) { tmx = tval; out.cv_worst_hold = h; }
    }
    out.n_cv = ncv; out.n_cv_skipped = nskip;
    out.cv_rms = ncv ? std::sqrt(sse / ncv) : 0.0;
    out.cv_max = mx;
    out.cv_tmax = tmx;
  }

  static float nth_(float* a, int m) {
    std::nth_element(a, a + m / 2, a + m);
    return a[m / 2];
  }

  void save_rob_() { for (int i = 0; i < n_; ++i) rob_save_[i] = rob_w_[i]; }
  void restore_rob_() { for (int i = 0; i < n_; ++i) rob_w_[i] = rob_save_[i]; }

  // Workspace (kept in the object, not on the fit-task stack)
  const AccelObs* obs_ = nullptr;
  int n_ = 0;
  const AccelFitCfg* cfg_ = nullptr;
  double T_ref_ = 25.0;
  double clamp_lo_ = -1000.0, clamp_hi_ = 1000.0;  // runtime temperature clamp of the model being fitted
  bool any_valid_temp_ = false;
  int n_fit_holds_ = 0, n_fit_blocks_ = 0, n_verify_blocks_ = 0;
  int hold_n_[MAXH];
  double hold_var_[MAXH];
  double sigma_block_ = 0;
  float base_w_[MAXO];     // 1 / (fit blocks of the hold); 0 for verify blocks
  float rob_w_[MAXO];      // robust weight
  float rob_save_[MAXO];
  float r_[MAXO];          // residual scratch (normalized)
  float scratch_[MAXO];
  Eigen::Vector3f cov_[MAXO];
};

}  // namespace imu_cal
