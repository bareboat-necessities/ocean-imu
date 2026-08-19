#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  ContinuousMagHardIronEstimator

  Exogenous, continuously running estimator of a body-fixed magnetometer
  offset, for a hull that only ever rolls and pitches.

  MagAutoTuner solves the same identification problem once, over a startup
  window, and its documentation says why that is hard: at a single attitude an
  offset and the world field are the same reading, and the tilt a hull takes
  through a wave breaks the tie only weakly.  This class keeps the accumulation
  running instead of closing it, so the window is no longer the startup window,
  and it prices the part of the problem that no amount of accumulation fixes.

  Model

      m_i = R_i^T B + b + n_i,

  with R_i the tilt-only rotation body->level.  Left-multiplying by R_i and
  eliminating the world field B from the weighted normal equations leaves

      (I - Abar^T Abar) b = mbar - Abar^T wbar,          Abar = mean(R_i),

  so B never becomes a filter state and the estimator never reads one.

  Why the solve is regularised

  The left matrix is the attitude excitation.  For roll and pitch of standard
  deviation sigma_r and sigma_p about a fixed heading it is, to second order,

      diag(sigma_p^2, sigma_r^2, sigma_r^2 + sigma_p^2),

  which for a few degrees of hull motion is of order 1e-3.  Inverting it
  multiplies everything the model does not carry by up to a thousand, and what
  it does not carry is soft iron and sensor misalignment -- both of which are
  modulated by the very tilt the offset is being read from, and neither of
  which averages away with time.  An unregularised fit therefore does not
  converge on the offset as the record lengthens; it converges on the offset
  plus an aliased image of the distortion, and in a wave record that image is
  as large as the offset itself.

  The ridge is that model error, not sample noise.  It is a prior variance ratio
  in the same dimensionless units as the excitation above: a direction excited
  far above it is returned nearly whole, and one excited far below it is
  returned nearly zero.  It does not shrink as samples accumulate, because the
  error it prices does not either.

  It has two terms and they are not two ways of saying the same thing.
  Config::model_ridge_relative carries the calibration and scales with the
  excitation, so what survives is a property of the sensor.  Config::model_ridge
  is a floor for a window with no excitation at all, and it only does that job
  while it stays below what a working hull produces -- above that it is a fixed
  ridge, shrinking hardest in the calm seas where the standing error is largest,
  which is precisely what the relative term exists to avoid.

  The class estimates only.  A caller decides whether, how much, and how
  quickly to apply the result, and levelReferenceForBias() gives it the
  magnetic reference belonging to whatever fraction it chose, computed from the
  same statistics, so a partly applied offset is never paired with a reference
  learned under a different one.
*/

#include <algorithm>
#include <cmath>

#ifdef EIGEN_NON_ARDUINO
  #include <Eigen/Dense>
#else
  #include <ArduinoEigenDense.h>
#endif

class ContinuousMagHardIronEstimator {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  struct Config {
    // Exponential-memory time constant.  <= 0 makes the statistics cumulative.
    float memory_sec = 600.0f;

    // The eigensolve does not belong on every magnetometer sample.
    float solve_period_sec = 1.0f;

    // Effective sample count the window must reach before a result is offered.
    float min_effective_weight = 500.0f;

    // Excitation the window must reach, as effective_weight * lambda_min of
    // the un-ridged normal matrix.  Low, and deliberately: this gate only has
    // to answer "has the hull moved at all", because how much of a poorly
    // excited window survives is already decided by the ridge below.  Set high
    // enough to be a second shrinkage and it stands the estimator down in
    // exactly the calm seas where the standing yaw error is largest.
    float min_information = 2.0f;

    // Absolute floor under the ridge; see the header note.  It exists for a
    // window with almost no excitation at all, and it only does that job while
    // it stays below the excitation a working hull produces.  Above that it
    // stops being a floor and becomes the whole regularisation, which is the
    // fixed-ridge behaviour the relative term below was introduced to replace.
    //
    // What sets the scale is the information gate above: at the saturated
    // weight of an exponential window -- memory_sec over the magnetometer
    // sample period -- the gate admits a direction of eigenvalue
    // min_information/weight, which is 1.3e-4 on these settings.  This floor
    // sits a few times above that and one to two orders of magnitude below
    // what even a slight sea excites, so the relative term governs whenever
    // the hull moves at all.
    float model_ridge = 5.0e-4f;

    // The part of the ridge that scales with the excitation, as a multiple of
    // the mean eigenvalue of the normal matrix.
    //
    // This is the term that matters, and the reason is that the error it
    // prices does not behave like noise.  The distortion aliased into the fit
    // enters the right-hand side through the same tilt that carries the
    // offset, so the resulting offset error is of order eps*|B| whatever the
    // sea state -- a hull working through thirty degrees is not fitting a
    // cleaner offset than one working through five, it is fitting the same
    // wrong one with more confidence.  A fixed ridge therefore shrinks hardest
    // exactly where the fit is most trustworthy and hardly at all where it is
    // least, which is backwards.  Scaling with the excitation makes the
    // fraction returned a property of the sensor rather than of the weather.
    //
    // The value is eps*|B| over the offset a caller expects: distortion of a
    // percent or two on a 50 uT field against an offset of a couple of uT.
    float model_ridge_relative = 0.5f;

    // Reject an implausible offset as a fraction of the measured field norm.
    float max_bias_fraction = 0.35f;

    // Reject a window the hard-iron model does not explain.  This is what
    // stands between the estimator and a hull that changes heading: the
    // accumulation frame is tied to the hull's own heading, so a turn makes
    // the world field move within it and the residual rises immediately.
    // <= 0 disables the gate.
    float max_residual_rms_uT = 3.0f;

    float min_mag_norm_uT = 1.0e-3f;
  };

  struct Estimate {
    Eigen::Vector3f bias_body_uT = Eigen::Vector3f::Zero();
    Eigen::Vector3f field_level_uT = Eigen::Vector3f::Zero();
    float information = 0.0f;
    float effective_weight = 0.0f;
    float residual_rms_uT = NAN;
    bool valid = false;
  };

  ContinuousMagHardIronEstimator() { reset(); }

  explicit ContinuousMagHardIronEstimator(const Config& cfg) : cfg_(cfg) {
    reset();
  }

  void setConfig(const Config& cfg) {
    cfg_ = cfg;
    reset();
  }

  const Config& config() const noexcept { return cfg_; }

  void reset() {
    weight_sum_ = 0.0;
    rot_sum_.setZero();
    level_mag_sum_.setZero();
    body_mag_sum_.setZero();
    mag_sq_sum_ = 0.0;
    solve_elapsed_sec_ = 0.0f;
    estimate_ = Estimate{};
  }

  // q_tilt_bw is the yaw-stripped body->level rotation.  Stripping yaw is what
  // makes the input exogenous in the strong sense: the frame is built from
  // gravity alone, so no heading estimate -- and no heading drift -- enters the
  // accumulation.  The level frame's own gauge is the hull's heading, which is
  // why a turn shows up as model residual rather than as a silent error.
  void update(float dt,
              const Eigen::Quaternionf& q_tilt_bw_in,
              const Eigen::Vector3f& mag_body_uT)
  {
    if (!(std::isfinite(dt) && dt > 0.0f)) return;
    if (!q_tilt_bw_in.coeffs().allFinite() || !mag_body_uT.allFinite()) return;

    const float mag_norm = mag_body_uT.norm();
    if (!(std::isfinite(mag_norm) && mag_norm > cfg_.min_mag_norm_uT)) return;

    Eigen::Quaternionf q = q_tilt_bw_in;
    const float qn = q.norm();
    if (!(std::isfinite(qn) && qn > 1.0e-6f)) return;
    q.normalize();

    const Eigen::Matrix3d R = q.toRotationMatrix().cast<double>();
    if (!R.allFinite()) return;

    // The normal matrix is I minus a product of near-identity rotations, so it
    // is a small difference of large terms.  The accumulators carry double for
    // that subtraction alone; everything the caller sees is float.
    const Eigen::Vector3d m = mag_body_uT.cast<double>();

    const double lambda =
        (std::isfinite(cfg_.memory_sec) && cfg_.memory_sec > 0.0f)
            ? std::exp(-double(dt) / double(cfg_.memory_sec))
            : 1.0;

    weight_sum_ *= lambda;
    rot_sum_ *= lambda;
    level_mag_sum_ *= lambda;
    body_mag_sum_ *= lambda;
    mag_sq_sum_ *= lambda;

    weight_sum_ += 1.0;
    rot_sum_ += R;
    level_mag_sum_ += R * m;
    body_mag_sum_ += m;
    mag_sq_sum_ += m.squaredNorm();

    solve_elapsed_sec_ += dt;
    if (solve_elapsed_sec_ >= std::max(0.02f, cfg_.solve_period_sec)) {
      solve_elapsed_sec_ = 0.0f;
      solve_();
    }
  }

  const Estimate& estimate() const noexcept { return estimate_; }
  bool valid() const noexcept { return estimate_.valid; }
  float information() const noexcept { return estimate_.information; }
  float effectiveWeight() const noexcept { return float(weight_sum_); }
  float residualRmsUT() const noexcept { return estimate_.residual_rms_uT; }

  // Magnetic reference in the level accumulation frame that belongs to an
  // applied offset, from the same statistics as the fit:
  //
  //     ref(b) = mean(R_i m_i) - mean(R_i) b.
  //
  // A caller applying a fraction of the estimate asks for the reference at
  // that fraction, so the two never disagree.
  bool levelReferenceForBias(const Eigen::Vector3f& applied_bias_body_uT,
                             Eigen::Vector3f& level_ref_uT) const
  {
    if (!(weight_sum_ > 1.0e-6) || !std::isfinite(weight_sum_)) return false;
    if (!applied_bias_body_uT.allFinite()) return false;

    const Eigen::Vector3d ref =
        (level_mag_sum_ - rot_sum_ * applied_bias_body_uT.cast<double>())
        / weight_sum_;

    level_ref_uT = ref.cast<float>();
    return level_ref_uT.allFinite() &&
           (level_ref_uT.norm() > cfg_.min_mag_norm_uT);
  }

private:
  void solve_() {
    Estimate out;
    out.effective_weight = float(weight_sum_);

    if (!(weight_sum_ >= double(cfg_.min_effective_weight)) ||
        !std::isfinite(weight_sum_)) {
      estimate_ = out;
      return;
    }

    const Eigen::Matrix3d A = rot_sum_ / weight_sum_;
    const Eigen::Vector3d wbar = level_mag_sum_ / weight_sum_;
    const Eigen::Vector3d mbar = body_mag_sum_ / weight_sum_;
    if (!A.allFinite() || !wbar.allFinite() || !mbar.allFinite()) {
      estimate_ = out;
      return;
    }

    Eigen::Matrix3d M = Eigen::Matrix3d::Identity() - A.transpose() * A;
    M = 0.5 * (M + M.transpose());
    if (!M.allFinite()) {
      estimate_ = out;
      return;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M);
    if (es.info() != Eigen::Success || !es.eigenvalues().allFinite()) {
      estimate_ = out;
      return;
    }

    // Reported on the un-ridged matrix, so it measures the hull and not the
    // prior.  The gate below is about whether the window says anything at all;
    // how much of what it says survives is the ridge's business.
    const double lambda_min = es.eigenvalues().minCoeff();
    out.information = float(weight_sum_ * lambda_min);
    if (!(std::isfinite(out.information) &&
          out.information >= cfg_.min_information)) {
      estimate_ = out;
      return;
    }

    const double ridge =
        std::max(0.0, double(cfg_.model_ridge)) +
        std::max(0.0, double(cfg_.model_ridge_relative)) * (M.trace() / 3.0);
    const Eigen::Vector3d rhs = mbar - A.transpose() * wbar;
    const Eigen::Vector3d b =
        (M + ridge * Eigen::Matrix3d::Identity()).ldlt().solve(rhs);
    if (!b.allFinite()) {
      estimate_ = out;
      return;
    }

    const Eigen::Vector3d B = wbar - A * b;
    if (!B.allFinite()) {
      estimate_ = out;
      return;
    }

    const double field_scale = std::sqrt(std::max(0.0, mag_sq_sum_ / weight_sum_));
    if (!(std::isfinite(field_scale) && field_scale > double(cfg_.min_mag_norm_uT))) {
      estimate_ = out;
      return;
    }
    if (b.norm() > double(cfg_.max_bias_fraction) * field_scale) {
      estimate_ = out;
      return;
    }

    // Weighted residual of w_i = B + R_i b.  R is orthonormal, so
    // sum ||R_i m_i||^2 == sum ||m_i||^2 and the whole sum is available from
    // the accumulators.
    const double sse =
        mag_sq_sum_
        - 2.0 * B.dot(level_mag_sum_)
        - 2.0 * b.dot(body_mag_sum_)
        + 2.0 * B.dot(rot_sum_ * b)
        + weight_sum_ * (B.squaredNorm() + b.squaredNorm());

    out.residual_rms_uT = float(std::sqrt(std::max(0.0, sse / weight_sum_)));

    if (std::isfinite(cfg_.max_residual_rms_uT) &&
        cfg_.max_residual_rms_uT > 0.0f &&
        (!std::isfinite(out.residual_rms_uT) ||
         out.residual_rms_uT > cfg_.max_residual_rms_uT)) {
      estimate_ = out;
      return;
    }

    out.bias_body_uT = b.cast<float>();
    out.field_level_uT = B.cast<float>();
    out.valid = true;
    estimate_ = out;
  }

  Config cfg_{};

  double weight_sum_ = 0.0;
  Eigen::Matrix3d rot_sum_ = Eigen::Matrix3d::Zero();
  Eigen::Vector3d level_mag_sum_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d body_mag_sum_ = Eigen::Vector3d::Zero();
  double mag_sq_sum_ = 0.0;

  float solve_elapsed_sec_ = 0.0f;
  Estimate estimate_{};
};
