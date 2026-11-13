#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "StandardHeaders.h"
#include <Eigen/Dense>
#define private public
#include "nlo/TimeVarGainNLO_Adapter.h"
#undef private
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
#include "GyroPredictionChecks.h"
int main(){
 using Nlo=TimeVarGainNloAdapter<false,NloMagType::None,double>;auto cfg=Nlo::makeDefaultConfig();
 check(cfg.gyro_sample_mode==Mode::Instantaneous,"NLO default");cfg.boot_bias_seed_min_still_s=0;cfg.boot_accel_bias_seed_enabled=false;
 Nlo n(cfg);cfg.gyro_sample_mode=Mode::IntervalRate;Nlo ref(cfg);
 using Pii=marine_obs::AdaptiveVerticalPIIMahony<double>;Pii::Config pc;check(pc.gyro_sample_mode==Mode::Instantaneous,"PII default");
 Pii p(pc);pc.gyro_sample_mode=Mode::IntervalRate;Pii pr(pc);
 V prev=V::Zero();bool handed=false;
 for(int k=0;k<36000;++k){const double t=k*.005;const M R=rotation(.1*t);const V w=.1*rate(.1*t),avg=k?V(.5*(prev+w)):V::Zero();prev=w;
  const V down=R.transpose()*V(0,0,-9.80665),up=-down;
  bool was=n.snapshot().initialized;n.update(.005,w,down);ref.update(.005,avg,down);
  const auto a=n.snapshot(),b=ref.snapshot();
  check((a.q_nb.coeffs()-b.q_nb.coeffs()).norm()<1e-9,"NLO default sampling path mismatch");
  if(!was&&a.initialized){check((a.q_nb.coeffs()-n.boot_q_nb_.coeffs()).norm()<1e-12,"NLO predicted handoff interval twice");handed=true;}
  const double z=p.updateIMU(w.x(),w.y(),w.z(),up.x(),up.y(),up.z(),.005);
  const double zr=pr.updateIMU(avg.x(),avg.y(),avg.z(),up.x(),up.y(),up.z(),.005);
  check(std::abs(z-zr)<1e-9,"PII default sampling path mismatch");
 }
 check(handed,"NLO handoff untested");std::cout<<"NLO and PII default sampling and NLO handoff PASS\n";
}
