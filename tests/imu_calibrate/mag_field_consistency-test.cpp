// Copyright 2026, Mikhail Grushinskiy
// Frozen-calibration temporal checks must not confuse attitude order with drift.
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/MagCalFit.h"
#include <array>
#include <cstdio>
#include <random>

using V = Eigen::Vector3f;
using M = Eigen::Matrix3f;
using Quality = imu_cal::MagFitQuality;
using Gate = imu_cal::MagFitGate;
static constexpr int N = 400;
static int failures = 0;
static void require(bool ok, const char* message) {
  if (!ok) { std::fprintf(stderr,"FAIL: %s\n",message); ++failures; }
}
static V direction(int i, int n) {
  const float x=-1.f+(2.f*i+1)/n, phase=2.39996323f*i;
  const float r=std::sqrt(1-x*x);
  return V(x,r*std::cos(phase),r*std::sin(phase));
}
struct Sweep {
  std::array<V,N> x;
  std::array<M,N> covariance;
  std::array<uint32_t,N> time;
  Sweep() {
    for (int i=0;i<N;++i) { x[i]=50*direction(i,N); covariance[i].setZero(); time[i]=400*i; }
  }
};
static void testFixedResidual() {
  imu_cal::MagGeometricFit<N> fit;
  const M A=M::Identity();
  for (float field : {25.f,50.f,65.f}) for (float error : {.45f,.6f,1.2f}) {
    if (field==25.f && error>0.6f) continue;
    Sweep sweep;
    const V bias(error,0,0), saved_bias=bias;
    const M saved_A=A;
    for (int i=0;i<N;++i) sweep.x[i]=field*direction(i,N);
    const auto saved_x=sweep.x;
    Quality q;
    require(fit.check(sweep.x.data(),N,field,A,bias,q),"fixed residual passes every absolute quality gate");
    require(fit.check(sweep.x.data(),N,field,A,bias,q,sweep.time.data()),
            "unchanged field passes an ordered roll/pitch sweep with frozen coefficients");
    // This is the old rejection: raw quarter means change as the same small
    // fixed hard-iron error is seen from different directions.
    require(q.time_drift>imu_cal::MagFitLimits::driftLimit(field),"reproduction exceeds the old point-estimate drift gate");
    require(q.time_drift_lower<=imu_cal::MagFitLimits::driftLimit(field),"unresolved apparent drift is not reported as a field change");
    require(A==saved_A && bias==saved_bias && sweep.x==saved_x,"temporal check never refits coefficients or edits observations");
    std::printf("constant field=%.0f rms=%.6f raw_drift=%.6f lower=%.6f limit=%.6f\n",
                field,q.rms,q.time_drift,q.time_drift_lower,imu_cal::MagFitLimits::driftLimit(field));
    const Quality ordered=q;
    // Reservoir slots are not chronological. Keep each sample and timestamp
    // together while shuffling; a clock-origin change must not matter either.
    std::mt19937 rng(411);
    for (int i=N-1;i>0;--i) {
      const int j=int(rng()%uint32_t(i+1));
      std::swap(sweep.x[i],sweep.x[j]); std::swap(sweep.time[i],sweep.time[j]);
    }
    for (auto& t:sweep.time)t+=1000000;
    require(fit.check(sweep.x.data(),N,field,A,bias,q,sweep.time.data()),"shuffled reservoir and shifted clock preserve acceptance");
    require(std::fabs(q.time_drift-ordered.time_drift)<1e-6 &&
            std::fabs(q.time_drift_lower-ordered.time_drift_lower)<1e-6,"temporal diagnostics depend on timestamps, not slots");
    q.time_drift_lower=100;
    require(fit.check(sweep.x.data(),N,field,A,bias,q) && q.time_drift==0 && q.time_drift_lower==0,
            "untimed callers reset all temporal diagnostics");
  }
}
static void testActualDrift() {
  imu_cal::MagGeometricFit<N> fit;
  const M A=M::Identity(); const V bias=V::Zero();
  for (float field : {25.f,50.f,65.f}) for (bool ordered : {false,true}) {
    const float change=field==25.f?.35f:(field==50.f?.4f:.5f);
    Sweep sweep;
    for (int i=0;i<N;++i) sweep.x[i]=(field+(i<200?-change:change))*(ordered?direction(i,N):direction(i%100,100));
    Quality q;
    require(fit.check(sweep.x.data(),N,field,A,bias,q),"small real drift is below all absolute error gates");
    require(!fit.check(sweep.x.data(),N,field,A,bias,q,sweep.time.data()) && q.gate==Gate::FIELD_CHANGED,
            "resolved small field change still fails with repeated or ordered directional coverage");
    require(q.time_drift_lower>imu_cal::MagFitLimits::driftLimit(field),"true temporal drift exceeds the unchanged physical limit");
    std::printf("changed field=%.0f rms=%.6f raw_drift=%.6f lower=%.6f limit=%.6f\n",
                field,q.rms,q.time_drift,q.time_drift_lower,imu_cal::MagFitLimits::driftLimit(field));
  }
  Sweep disturbed;
  for (int i=0;i<18;++i)disturbed.x[i]*=1.08f;
  Quality q;
  require(!fit.check(disturbed.x.data(),N,50.f,A,bias,q,disturbed.time.data()) && q.gate==Gate::FIELD_CHANGED,
          "short disturbed quarter still fails the unchanged local inlier gate");
  Sweep shifted;
  for (auto& v:shifted.x)v+=V(10,0,0);
  require(!fit.check(shifted.x.data(),N,50.f,A,bias,q,shifted.time.data()),"field change between fit and verification is not absorbed");
  Sweep poor;
  const V bad_bias(3,0,0);
  require(!fit.check(poor.x.data(),N,50.f,A,bad_bias,q,poor.time.data()) &&
          (q.gate==Gate::INLIERS || q.gate==Gate::RESIDUAL),"large calibration error still fails absolute quality");
}
static void testFixedSoftIron() {
  Sweep sweep;
  imu_cal::MagGeometricFit<N> fit;
  M A=M::Identity(); A.diagonal()<<1.022f,.986f,.988f;
  A(0,1)=A(1,0)=.004f;
  const V bias=V::Zero();
  Quality q;
  require(fit.check(sweep.x.data(),N,50.f,A,bias,q),"small fixed soft-iron residual meets absolute quality");
  require(fit.check(sweep.x.data(),N,50.f,A,bias,q,sweep.time.data()),"changing orientation does not turn fixed soft-iron error into field drift");
  require(q.time_drift>imu_cal::MagFitLimits::driftLimit(50),"soft-iron example reproduces the old false rejection");
}
static void testNoisyMoments() {
  imu_cal::MagGeometricFit<N> fit;
  const M A=M::Identity(); const V bias(.55f,0,0);
  for (int seed=0;seed<12;++seed) {
    Sweep sweep;
    std::mt19937 rng(1900+seed); std::normal_distribution<float> noise(0,.3f);
    for (int i=0;i<N;++i) {
      // Reproduce window moments through faster within-window rotation, with
      // the starting attitude progressing slowly across successive quarters.
      Eigen::Vector3d mean=Eigen::Vector3d::Zero();
      Eigen::Matrix3d scatter=Eigen::Matrix3d::Zero();
      const V u=direction(i,N);
      for (int k=0;k<12;++k) {
        const float angle=.03f*(k-5.5f);
        V v(50*u.x(),50*(std::cos(angle)*u.y()-std::sin(angle)*u.z()),
            50*(std::sin(angle)*u.y()+std::cos(angle)*u.z()));
        v+=V(noise(rng),noise(rng),noise(rng));
        const Eigen::Vector3d d=v.cast<double>()-mean;
        mean+=d/(k+1); scatter+=d*(v.cast<double>()-mean).transpose();
      }
      sweep.x[i]=mean.cast<float>();
      sweep.covariance[i]=((scatter+scatter.transpose())/24).cast<float>();
    }
    Quality q;
    require(fit.check(sweep.x.data(),N,50.f,A,bias,q,sweep.time.data(),.15,sweep.covariance.data()),
            "constant field passes noisy moving-window moments across independent seeds");
  }
}
int main() {
  testFixedResidual(); testActualDrift(); testFixedSoftIron(); testNoisyMoments();
  std::printf("mag_field_consistency: %s (%d failures)\n",failures?"FAIL":"PASS",failures);
  return failures?1:0;
}
