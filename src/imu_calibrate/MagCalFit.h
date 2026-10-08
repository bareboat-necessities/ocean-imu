#pragma once

/* Copyright 2026, Mikhail Grushinskiy

   Geometric, rotation-free magnetic calibration refinement and qualification.
   Field scale is fixed by the caller; only the six SPD terms and three offsets
   are fitted. Fixed-capacity workspace belongs to the calibrator, not the task
   stack. Quality is evaluated again on stored floats and independent data.
*/
#include <algorithm>
#include <cmath>
#include <cstdint>
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

namespace imu_cal {

struct MagFitLimits {
  // Broad device plausibility bounds; quality limits are in corrected uT.
  static constexpr double min_gain = 0.2, max_gain = 5.0, max_condition = 10.0;
  static constexpr double max_bias_uT = 150.0;
  static constexpr double min_inlier_fraction = 0.85;
  static constexpr int min_cells = 12;
  static double rmsLimit(double field) { return std::max(0.35, 0.015 * field); }
  static double inlierLimit(double field) { return std::max(0.6, 0.025 * field); }
  static double driftLimit(double field) { return std::max(0.35, 0.01 * field); }
  // Individual fresh readings of the independent check against the frozen
  // candidate: window limits plus raw sensor noise and the axis skew of one
  // measurement while turning. A window's moments can hide field changes
  // inside it; single readings cannot.
  static double rawRmsLimit(double field) { return std::max(2.0, 0.04 * field); }
  static double rawTailLimit(double field) { return std::max(5.0, 0.1 * field); }
  static constexpr double max_raw_tail_fraction = 0.05;
  static constexpr int min_raw_samples = 100;
};

enum class MagFitGate : uint8_t { NONE, BAD_DATA, MATRIX, RESIDUAL, INLIERS, COVERAGE, INFORMATION, FIELD_CHANGED };
inline const char* magFitGateText(MagFitGate gate) {
  switch (gate) {
    case MagFitGate::MATRIX: return "Correction too large";
    case MagFitGate::RESIDUAL: return "Field inconsistent";
    case MagFitGate::INLIERS: return "MAG samples noisy";
    case MagFitGate::COVERAGE: return "Need more 3D motion";
    case MagFitGate::INFORMATION: return "Need more directions";
    case MagFitGate::FIELD_CHANGED: return "Field changed";
    case MagFitGate::BAD_DATA: return "No usable MAG data";
    default: return "OK";
  }
}

struct MagFitQuality {
  MagFitGate gate = MagFitGate::BAD_DATA;
  int samples = 0, inliers = 0, cells = 0, iterations = 0;
  double rms = 0, trimmed_rms = 0, p95 = 0, max_bias_sigma = 0, max_matrix_sigma = 0, time_drift = 0;
  // Raw quarter-mean range and a direction-separated temporal lower bound are
  // separate diagnostics; neither changes the absolute residual gates.
  double time_drift_lower = 0;
  double initial_cost = 0, refined_cost = 0;
};

// 26 broad direction cells. No cell names imply an exact physical pose.
template<class Derived>
inline int magDirectionCell(const Eigen::MatrixBase<Derived>& v) {
  const double n = v.norm();
  if (!(n > 1e-9)) return 13;
  int cell = 0;
  for (int j = 0; j < 3; ++j) {
    const double u = double(v[j]) / n;
    cell = 3 * cell + (u < -0.45 ? 0 : (u > 0.45 ? 2 : 1));
  }
  return cell;
}

template<int N>
class MagGeometricFit {
  using V = Eigen::Vector3d;
  using M = Eigen::Matrix3d;
  using P = Eigen::Matrix<double,9,1>;
  using H9 = Eigen::Matrix<double,9,9>;
public:
  template<typename T>
  bool refine(const Eigen::Matrix<T,3,1>* x, int n, T field,
              Eigen::Matrix<T,3,3>& A, Eigen::Matrix<T,3,1>& bias, MagFitQuality& q,
              const Eigen::Matrix<T,3,3>* sample_cov = nullptr) {
    q = MagFitQuality{};
    if (!x || n > N) return false;
    if (n < 80 || !(field >= T(12) && field <= T(120))) return false;
    const M initial = A.template cast<double>();
    Eigen::LLT<M> llt(initial);
    if (llt.info() != Eigen::Success) { q.gate = MagFitGate::MATRIX; return false; }
    const M L = llt.matrixL();
    P th;
    th << std::log(L(0,0)), L(1,0), std::log(L(1,1)), L(2,0), L(2,1), std::log(L(2,2)),
          double(bias[0])/field, double(bias[1])/field, double(bias[2])/field;
    int counts[27]{};
    for (int i = 0; i < n; ++i) {
      cells_[i] = uint8_t(magDirectionCell(initial * (x[i].template cast<double>() - bias.template cast<double>())));
      ++counts[cells_[i]];
    }
    for (int i = 0; i < n; ++i) base_[i] = 1.0 / counts[cells_[i]];
    // A fixed physical Huber knee makes costs comparable throughout refinement.
    const double knee = std::max(0.15, 0.004 * double(field)) / field;
    auto cost = [&](const P& p) {
      const M l = chol_(p), s = l * l.transpose();
      double c = 0;
      for (int i = 0; i < n; ++i) {
        const V v=s*(x[i].template cast<double>()/double(field)-p.template tail<3>());
        const M cov=sampleCov_(sample_cov,i,double(field));
        const double a=std::fabs(std::sqrt(v.squaredNorm()+(s*cov*s.transpose()).trace())-1);
        if (!std::isfinite(a)) return double(INFINITY);
        c += base_[i] * (a <= knee ? 0.5*a*a : knee*(a-0.5*knee));
      }
      return c;
    };
    double c = cost(th), lambda = 1e-4;
    if (!std::isfinite(c)) return false;
    q.initial_cost = c;
    for (int iteration = 0; iteration < 35; ++iteration) {
      const M l = chol_(th), s = l*l.transpose();
      H9 h = H9::Zero(); P g = P::Zero();
      for (int i = 0; i < n; ++i) {
        const V u = x[i].template cast<double>() / double(field) - th.template tail<3>();
        const V v = s*u;
        const M cov=sampleCov_(sample_cov,i,double(field));
        const double r=std::sqrt(v.squaredNorm()+(s*cov*s.transpose()).trace());
        if (!(r > 1e-9)) return false;
        const double e = r-1, w = base_[i] * std::min(1.0, knee/std::max(std::fabs(e),1e-15));
        const M G = (v*u.transpose()+s*cov)/r, d = (G+G.transpose())*l;
        P j;
        j << d(0,0)*l(0,0), d(1,0), d(1,1)*l(1,1), d(2,0), d(2,1), d(2,2)*l(2,2), 0,0,0;
        j.template tail<3>() = -s*v/r;
        h.noalias() += w*j*j.transpose(); g.noalias() += w*e*j;
      }
      bool accepted = false;
      P step = P::Zero();
      for (int attempt = 0; attempt < 10; ++attempt) {
        H9 damped = h;
        for (int j = 0; j < 9; ++j) damped(j,j) += lambda*std::max(1e-9,h(j,j));
        Eigen::LDLT<H9> solve(damped);
        step = solve.solve(-g);
        if (solve.info() != Eigen::Success || !step.allFinite()) { lambda *= 10; continue; }
        const P next = th+step;
        const double nc = cost(next);
        if (std::isfinite(nc) && nc <= c) { th = next; c = nc; accepted = true; break; }
        lambda *= 10;
      }
      q.iterations = iteration+1;
      if (!accepted || step.cwiseAbs().maxCoeff() < 1e-8) break;
      lambda = std::max(1e-10, lambda/3);
    }
    const M l = chol_(th);
    const M result = l*l.transpose();
    A = (0.5*(result+result.transpose())).template cast<T>();
    bias = (th.template tail<3>() * double(field)).template cast<T>();
    q.refined_cost = c;
    // Check exactly the precision the runtime will apply, including its matrix.
    return check(x,n,field,A,bias,q,nullptr,0.15,sample_cov);
  }

  template<typename T>
  bool check(const Eigen::Matrix<T,3,1>* x, int n, T field,
             const Eigen::Matrix<T,3,3>& A, const Eigen::Matrix<T,3,1>& bias,
             MagFitQuality& q, const uint32_t* elapsed_ms = nullptr, double report_trim = 0.15,
             const Eigen::Matrix<T,3,3>* sample_cov = nullptr) {
    q.gate = MagFitGate::BAD_DATA; q.samples = n; q.inliers = q.cells = 0;
    q.rms = q.trimmed_rms = q.p95 = q.time_drift = q.time_drift_lower = q.max_bias_sigma = q.max_matrix_sigma = 0;
    if (!x || n > N) return false;
    if (n < 80 || !std::isfinite(double(field)) || field < T(12) || field > T(120) ||
        !std::isfinite(report_trim)) return false;
    if (!matrixValid(A,bias)) { q.gate = MagFitGate::MATRIX; return false; }
    const M& s = A.template cast<double>();
    const double limit = MagFitLimits::inlierLimit(field);
    int count[27]{}, positive[3]{}, negative[3]{};
    double sum2 = 0;
    M coverage = M::Zero();
    for (int i = 0; i < n; ++i) {
      // Force evaluation in T, as in MagCalibration<T>::apply().
      const Eigen::Matrix<T,3,1> corrected = A*(x[i]-bias);
      const V& v = corrected.template cast<double>();
      // Covariance here describes actual within-window motion/noise, not
      // uncertainty of the mean. This is exactly RMS |A(raw-b)| for the group.
      const M cov=sampleCov_(sample_cov,i,1.0);
      const double radius=std::sqrt(v.squaredNorm()+(s*cov*s.transpose()).trace());
      if (!std::isfinite(radius) || !(radius > 1e-9) || !(v.norm()>1e-9)) return false;
      errors_[i] = radius-double(field);
      scratch_[i] = std::fabs(errors_[i]);
      cells_[i] = uint8_t(magDirectionCell(v));
      if (scratch_[i] > limit) continue;
      ++q.inliers; ++count[cells_[i]]; sum2 += errors_[i]*errors_[i];
      const V u = v.normalized();
      coverage.noalias() += u*u.transpose();
      for (int j = 0; j < 3; ++j) { positive[j] += u[j]>0; negative[j] += u[j]<0; }
    }
    std::sort(scratch_,scratch_+n);
    // Preserve the public trimmed-RMS report, evaluated on the refined model.
    // It is never used to relax acceptance: all inliers and tails qualify below.
    const int kept = std::max(10,int(std::floor(n*(1-std::max(0.0,std::min(0.45,report_trim))))));
    for(int i=0;i<kept;++i) q.trimmed_rms += scratch_[i]*scratch_[i];
    q.trimmed_rms = std::sqrt(q.trimmed_rms/kept);
    q.p95 = scratch_[std::max(0,int(std::ceil(0.95*n))-1)];
    q.rms = q.inliers ? std::sqrt(sum2/q.inliers) : INFINITY;
    if (q.inliers < int(std::ceil(MagFitLimits::min_inlier_fraction*n))) { q.gate = MagFitGate::INLIERS; return false; }
    if (q.rms > MagFitLimits::rmsLimit(field) || q.p95 > 2*limit) { q.gate = MagFitGate::RESIDUAL; return false; }
    for (int j = 0; j < 27; ++j) q.cells += count[j]>=3;
    bool signs = true;
    for (int j = 0; j < 3; ++j) signs = signs && std::min(positive[j],negative[j]) >= int(std::ceil(0.08*q.inliers));
    coverage /= q.inliers;
    if (!signs || q.cells < MagFitLimits::min_cells || coverage.determinant() < 1e-3) {
      q.gate = MagFitGate::COVERAGE; return false;
    }
    // Information counts one observation per direction cell, not every highly
    // correlated raw reading. Physical parameter order: diag, cross, bias/B.
    uint32_t lo = 0, hi = 0;
    if (elapsed_ms) {
      lo = hi = elapsed_ms[0];
      for (int i = 1; i < n; ++i) { lo = std::min(lo,elapsed_ms[i]); hi = std::max(hi,elapsed_ms[i]); }
    }
    auto timeBin = [&](int i) {
      return std::min(3,int(uint64_t(elapsed_ms[i]-lo)*4/(uint64_t(hi-lo)+1)));
    };
    auto jacobian = [&](int i) {
      const V u = (x[i].template cast<double>()-bias.template cast<double>())/double(field), v = s*u;
      const M cov=sampleCov_(sample_cov,i,double(field));
      const double r=std::sqrt(v.squaredNorm()+(s*cov*s.transpose()).trace());
      const M G = (v*u.transpose()+s*cov)/r;
      P j;
      j << G(0,0),G(1,1),G(2,2),G(0,1)+G(1,0),G(0,2)+G(2,0),G(1,2)+G(2,1),0,0,0;
      j.template tail<3>() = -s*v/r;
      return j;
    };
    H9& h = information_;
    h.setZero();
    for (int i = 0; i < n; ++i) {
      if (std::fabs(errors_[i]) > limit) continue;
      const P j = jacobian(i);
      h.noalias() += j*j.transpose()/double(count[cells_[i]]);
    }
    auto& info = information_solve_;
    info.compute(h);
    if (info.info()!=Eigen::Success || info.vectorD().minCoeff() <= 1e-6*std::max(1.0,h.trace())) {
      q.gate = MagFitGate::INFORMATION; return false;
    }
    covariance_ = info.solve(H9::Identity());
    const H9& covariance = covariance_;
    const double noise = std::max(0.2,q.rms);
    for (int j=0;j<9;++j) {
      const double sd = noise*std::sqrt(std::max(0.0,covariance(j,j)));
      if (j<6) q.max_matrix_sigma = std::max(q.max_matrix_sigma,sd/double(field));
      else q.max_bias_sigma = std::max(q.max_bias_sigma,sd);
    }
    if (!covariance.allFinite() || q.max_bias_sigma > std::max(0.5,0.02*double(field)) || q.max_matrix_sigma > 0.05) {
      q.gate = MagFitGate::INFORMATION; return false;
    }
    if (elapsed_ms) {
      double sums[4]{}, means[4]{};
      int counts[4]{}, total[4]{};
      for(int i=0;i<n;++i) {
        const int bin = timeBin(i);
        ++total[bin];
        if(std::fabs(errors_[i])>limit) continue;
        sums[bin]+=errors_[i]; ++counts[bin];
      }
      // A short disturbed interval must not disappear into the global trim.
      for(int j=0;j<4;++j) if(total[j]>=10 && counts[j]<MagFitLimits::min_inlier_fraction*total[j]) {
        q.gate=MagFitGate::FIELD_CHANGED;return false;
      }
      int reference=0;
      for (int j=0;j<4;++j) {
        if (counts[j]) means[j]=sums[j]/counts[j];
        if (counts[j]>counts[reference]) reference=j;
      }
      for (int a=0;a<4;++a) for (int b=a+1;b<4;++b) if (counts[a]>=10 && counts[b]>=10)
        q.time_drift=std::max(q.time_drift,std::fabs(means[a]-means[b]));

      // Joint diagnostic regression: radial error = a fixed direction/moment
      // pattern (the nine calibration Jacobians) + a time-quarter offset.
      // Eliminate the fixed pattern before testing time contrasts. This does
      // NOT update A/bias, refit the calibration or relax any absolute gate.
      // Joint fitting is important: simply detrending errors first would also
      // attenuate real drift when time and attitude are correlated.
      int slot[4]={-1,-1,-1,-1}, k=0;
      for (int j=0;j<4;++j) if (j!=reference && counts[j]) slot[j]=k++;
      if (k) {
        h.setZero(); time_cross_.setZero();
        P gradient=P::Zero();
        M temporal=M::Identity(); V rhs=V::Zero();
        for (int j=0;j<4;++j) if (slot[j]>=0) temporal(slot[j],slot[j])=counts[j];
        double squared_error=0;
        for (int i=0;i<n;++i) {
          if (std::fabs(errors_[i])>limit) continue;
          const P j=jacobian(i);
          h.noalias()+=j*j.transpose(); gradient.noalias()+=errors_[i]*j;
          squared_error+=errors_[i]*errors_[i];
          const int t=slot[timeBin(i)];
          if (t>=0) { time_cross_.col(t)+=j; rhs[t]+=errors_[i]; }
        }
        // Reuse the qualification workspace only after its gates have passed.
        // Here observations are non-overlapping retained windows. The separate
        // calibration information gate above remains direction-balanced.
        info.compute(h);
        if (info.info()!=Eigen::Success || info.vectorD().minCoeff()<=0) {
          q.gate=MagFitGate::INFORMATION; return false;
        }
        covariance_=info.solve(H9::Identity());
        const P baseline=covariance_*gradient;
        temporal-=time_cross_.transpose()*covariance_*time_cross_;
        rhs-=time_cross_.transpose()*baseline;
        temporal=(0.5*(temporal+temporal.transpose())).eval();
        if (!covariance_.allFinite() || !temporal.allFinite() || !rhs.allFinite()) {
          q.gate=MagFitGate::INFORMATION; return false;
        }
        Eigen::SelfAdjointEigenSolver<M> spectrum(temporal);
        if (spectrum.info()!=Eigen::Success) { q.gate=MagFitGate::INFORMATION; return false; }
        // Without enough within-quarter direction diversity, temporal and
        // fixed errors can be indistinguishable. That is not evidence that
        // the field changed; the unchanged absolute/local gates still apply.
        if (spectrum.eigenvalues().minCoeff()>1e-8*std::max(1.0,temporal.trace())) {
          Eigen::LDLT<M> solve(temporal);
          const V offsets=solve.solve(rhs);
          const M offset_covariance=solve.solve(M::Identity());
          const P fixed=baseline-covariance_*time_cross_*offsets;
          V raw_rhs=V::Zero();
          for (int j=0;j<4;++j) if (slot[j]>=0) raw_rhs[slot[j]]=sums[j];
          const double residual=std::max(0.0,squared_error-gradient.dot(fixed)-raw_rhs.dot(offsets));
          const double variance=std::max(0.2*0.2,residual/(q.inliers-9-k));
          if (!offsets.allFinite() || !offset_covariance.allFinite() || !std::isfinite(variance)) {
            q.gate=MagFitGate::INFORMATION; return false;
          }
          for (int a=0;a<4;++a) for (int b=a+1;b<4;++b) if (counts[a]>=10 && counts[b]>=10) {
            V contrast=V::Zero();
            if (slot[a]>=0) contrast[slot[a]]=1;
            if (slot[b]>=0) contrast[slot[b]]=-1;
            const double drift=std::fabs(contrast.dot(offsets));
            const double sigma=std::sqrt(variance*std::max(0.0,contrast.dot(offset_covariance*contrast)));
            q.time_drift_lower=std::max(q.time_drift_lower,drift-3*sigma);
          }
        }
      }
      // This uncertainty proxy is not a traceable confidence certificate.
      if(q.time_drift_lower>MagFitLimits::driftLimit(field)) {q.gate=MagFitGate::FIELD_CHANGED;return false;}
    }
    q.gate = MagFitGate::NONE;
    return true;
  }

  template<typename T>
  static bool matrixValid(const Eigen::Matrix<T,3,3>& A, const Eigen::Matrix<T,3,1>& b) {
    if (!A.allFinite() || !b.allFinite() || b.norm()>MagFitLimits::max_bias_uT ||
        (A-A.transpose()).norm()>1e-5*A.norm()) return false;
    const M& s = A.template cast<double>();
    Eigen::SelfAdjointEigenSolver<M> es(s);
    if(es.info()!=Eigen::Success) return false;
    const auto& e = es.eigenvalues();
    return e.minCoeff()>=MagFitLimits::min_gain && e.maxCoeff()<=MagFitLimits::max_gain &&
           e.maxCoeff()/e.minCoeff()<=MagFitLimits::max_condition;
  }
private:
  template<typename T>
  static M sampleCov_(const Eigen::Matrix<T,3,3>* covariance,int i,double scale) {
    if(!covariance)return M::Zero();
    return covariance[i].template cast<double>()/(scale*scale);
  }
  static M chol_(const P& p) {
    M l=M::Zero(); l(0,0)=std::exp(p[0]);l(1,0)=p[1];l(1,1)=std::exp(p[2]);
    l(2,0)=p[3];l(2,1)=p[4];l(2,2)=std::exp(p[5]);return l;
  }
  double base_[N], errors_[N], scratch_[N];
  uint8_t cells_[N];
  // Save/read-back validation also runs on the small caller task stack.
  // Keep its largest matrices with the sample workspace on the wizard heap.
  H9 information_, covariance_;
  Eigen::Matrix<double,9,3> time_cross_;
  Eigen::LDLT<H9> information_solve_;
};
} // namespace imu_cal
