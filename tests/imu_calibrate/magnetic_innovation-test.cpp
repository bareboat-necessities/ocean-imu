#define EIGEN_NON_ARDUINO
#include "kalman_common/MagneticInnovationGate.h"
#include <random>
#include <cstdio>
using Gate=ocean_imu::MagneticInnovationGate<double>;using V=Eigen::Vector3d;using M=Eigen::Matrix3d;
static int errors=0;static void check(bool ok,const char* s){if(!ok){++errors;std::fprintf(stderr,"FAIL %s\n",s);}}
int main(){
  Gate gate;M R;R<<1,.8,.2,.8,1,.1,.2,.1,2;const M C=.2*M::Identity();
  V r(1,-1,.5);auto d=gate.evaluate(r,C,R);
  check(std::fabs(d.nis-r.dot((C+R).ldlt().solve(r)))<1e-12,"full correlated covariance NIS");
  const M Q=Eigen::AngleAxisd(.7,V(.3,.1,.7).normalized()).toRotationMatrix();
  Gate rotated;auto e=rotated.evaluate(V(Q*r),M(Q*C*Q.transpose()),M(Q*R*Q.transpose()));
  check(std::fabs(e.nis-d.nis)<1e-12 && e.accepted==d.accepted,"NIS is frame invariant");
  Gate inflate;
  auto moderate=inflate.evaluate(V(std::sqrt(15.),0,0),M::Zero(),M::Identity());
  check(moderate.accepted && moderate.inflation>1,"moderate vector innovation inflates R");
  check(!inflate.evaluate(V(20,0,0),M::Zero(),M::Identity()).accepted,"gross outlier rejected");
  for(unsigned i=1;i<Gate::clean_required;++i)
    check(!inflate.evaluate(V::Zero(),M::Zero(),M::Identity()).accepted,"hysteresis requires clean streak");
  const auto recovery=inflate.evaluate(V::Zero(),M::Zero(),M::Identity());
  check(recovery.accepted && recovery.inflation>=16 && recovery.health==Gate::Health::RECOVERING,"recovery starts gently");
  Gate indefinite;M bad=R;bad(0,0)=-1;
  check(!indefinite.evaluate(r,C,bad).accepted,"indefinite R rejected");
  check(!indefinite.evaluate(V(NAN,0,0),C,R).accepted,"non-finite sample rejected");
  // Gaussian Monte Carlo checks a vector threshold, not three scalar gates.
  std::mt19937 rng(119);std::normal_distribution<double> n;
  const M L=(C+R).llt().matrixL();int gross=0;
  for(int i=0;i<100000;++i) {Gate fresh;const auto v=fresh.evaluate(V(L*V(n(rng),n(rng),n(rng))),C,R);
    if(!v.accepted)++gross;}
  check(gross<35,"chi-square gross rejection near expected .01 percent");
  std::printf("magnetic NIS: correlated vector Gaussian gross rejects=%d/100000\n",gross);
  return errors?1:0;
}
