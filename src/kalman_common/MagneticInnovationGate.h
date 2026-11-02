#pragma once
// Copyright 2026, Mikhail Grushinskiy
// One vector-NIS policy for all Kalman magnetic updates. R is a FULL covariance
// in the residual frame. Inflation modifies R, never HPH' or the innovation.
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif
#include <cmath>
#include <cstdint>
#include <limits>
#include <algorithm>

namespace ocean_imu {
template<typename T>
class MagneticInnovationGate {
public:
  using V=Eigen::Matrix<T,3,1>;using M=Eigen::Matrix<T,3,3>;
  enum class Health:uint8_t {HEALTHY,DISTURBED,RECOVERING};
  struct Decision {
    bool accepted=false;T nis=std::numeric_limits<T>::quiet_NaN(),inflation=T(1);
    Health health=Health::HEALTHY;uint32_t rejected=0,accepted_count=0;
  };
  // chi-square(3), quantiles .99 and .9999. These are NOT scalar/axis gates.
  static constexpr double soft_threshold=11.3448667301444;
  static constexpr double gross_threshold=21.1075134661604;
  static constexpr unsigned clean_required=10;
  const Decision& last() const {return last_;}
  void reset() {health_=Health::HEALTHY;clean_=0;inflation_=T(1);last_={};}
  static bool covarianceValid(const M& R) {
    if(!R.allFinite() || (R-R.transpose()).norm()>T(1e-5)*std::max(T(1e-12),R.norm()))return false;
    const Eigen::LLT<M> f(R);return f.info()==Eigen::Success;
  }
  Decision evaluate(const V& r,const M& HPH,const M& R,T field_ratio=T(1)) {
    last_.accepted=false;last_.nis=std::numeric_limits<T>::quiet_NaN();last_.inflation=T(1);
    if(!r.allFinite() || !HPH.allFinite() || !covarianceValid(R) ||
       !std::isfinite(field_ratio) || field_ratio<T(.65) || field_ratio>T(1.5))return reject_();
    const M C=(T(.5)*(HPH+HPH.transpose())).eval();
    auto nis=[&](T alpha) {
      const M S=C+alpha*R;const Eigen::LDLT<M> f(S);
      if(f.info()!=Eigen::Success || !f.isPositive())return std::numeric_limits<T>::quiet_NaN();
      const T v=r.dot(f.solve(r));return (std::isfinite(v) && v>=T(0))?v:std::numeric_limits<T>::quiet_NaN();
    };
    const T raw=nis(T(1));last_.nis=raw;
    if(!std::isfinite(raw) || raw>T(gross_threshold))return reject_();
    if(health_==Health::DISTURBED) {
      // Only nominal-covariance innovations count as clean. Inflated R must
      // never make a disturbed sample qualify its own recovery.
      if(raw>T(soft_threshold)){clean_=0;return reject_();}
      if(++clean_<clean_required) {last_.health=health_;return last_;}
      health_=Health::RECOVERING;inflation_=T(16);
    }
    T alpha=std::max(T(1),inflation_);
    if(raw>T(soft_threshold)) {
      T lo=T(1),hi=T(64);
      const T upper=nis(hi);
      if(!std::isfinite(upper) || upper>T(soft_threshold))return reject_();
      for(int i=0;i<10;++i) {const T mid=(lo+hi)/T(2);
        if(nis(mid)>T(7.81472790325118))lo=mid;else hi=mid;}
      alpha=std::max(alpha,hi);
    }
    last_.inflation=alpha;last_.accepted=true;last_.health=health_;++last_.accepted_count;
    inflation_=T(1)+T(.85)*(alpha-T(1));
    if(health_==Health::RECOVERING && inflation_<T(1.05)) {health_=Health::HEALTHY;inflation_=T(1);}
    return last_;
  }
private:
  Decision reject_() {health_=Health::DISTURBED;clean_=0;last_.health=health_;++last_.rejected;return last_;}
  Health health_=Health::HEALTHY;unsigned clean_=0;T inflation_=T(1);Decision last_{};
};
} // namespace ocean_imu
