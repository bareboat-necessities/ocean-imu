#define EIGEN_NON_ARDUINO
#ifndef MAG_FILTER_KIND
#define MAG_FILTER_KIND 1
#endif
#if MAG_FILTER_KIND==1
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#elif MAG_FILTER_KIND==2
#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#elif MAG_FILTER_KIND==3
#include "kalman_carma/Kalman3D_Wave_TFG_CARMA3.h"
#elif MAG_FILTER_KIND==4
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#elif MAG_FILTER_KIND==6
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#elif MAG_FILTER_KIND==5
#include "ahrs/KalmanQMEKF.h"
#endif
#include <cstdio>
#include <memory>
using V=Eigen::Vector3f;using M=Eigen::Matrix3f;
#if MAG_FILTER_KIND==1
using Filter=ocean_imu::kalman::Kalman3D_Wave_OU_III<float>;
#elif MAG_FILTER_KIND==2
using Filter=ocean_imu::kalman::Kalman3D_Wave_CARMA3<float>;
#elif MAG_FILTER_KIND==3
using Filter=ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3<float>;
#elif MAG_FILTER_KIND==4
using Filter=ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
#elif MAG_FILTER_KIND==6
using Filter=ocean_imu::kalman::Kalman3D_Wave_OU_II<float>;
#else
using Filter=QuaternionMEKF<float>;
#endif
static int errors=0;static void check(bool ok,const char* s){if(!ok){++errors;std::fprintf(stderr,"FAIL filter %d %s\n",MAG_FILTER_KIND,s);}}
static Eigen::Quaternionf quat(const Filter& f) {
#if MAG_FILTER_KIND==5
  return Eigen::Quaternionf(f.quaternion());
#else
  return f.quaternion();
#endif
}
int main(){
  const V a(0,0,-9.80665f),m(20,0,43);
#if MAG_FILTER_KIND==4
  auto f=std::make_unique<Filter>();f->initialize_identity(1e-6f);f->set_magnetic_reference_world(m);
#else
  auto f=std::make_unique<Filter>(V::Constant(.12f),V::Constant(.001f),V::Constant(.3f),1e-6f,1e-8f);
  f->initialize_from_acc_mag(a,m);
#endif
  M R;R<<.16f,.04f,.01f,.04f,.25f,.02f,.01f,.02f,.36f;
  check(f->set_Rmag_covariance(R),"accept correlated covariance");
  M bad=R;bad(0,0)=-1;check(!f->set_Rmag_covariance(bad),"reject invalid covariance");
  for(int i=0;i<30;++i)f->measurement_update_mag_only(m);
  const auto before=quat(*f);
#if MAG_FILTER_KIND==5
  const auto covariance_before=f->covariance();
#else
  const auto covariance_before=f->covariance_full();
#endif
  // Bottom tap regression: a same-norm field appears 90 degrees away in
  // heading. It must never reset the reference or inject a heading flip.
  const V shock(0,-20,43);
  for(int i=0;i<100;++i)f->measurement_update_mag_only(shock);
  check(quat(*f).angularDistance(before)<1e-5f,"90-degree field disturbance cannot flip heading");
#if MAG_FILTER_KIND==5
  check(f->covariance()==covariance_before,"rejected measurement leaves covariance unchanged");
#else
  check(f->covariance_full()==covariance_before,"rejected measurement leaves covariance unchanged");
#endif
  check(!f->magneticHealth().accepted && f->magneticHealth().rejected>=100,"health records rejection");
  for(unsigned i=1;i<ocean_imu::MagneticInnovationGate<float>::clean_required;++i) {
    f->measurement_update_mag_only(m);check(!f->magneticHealth().accepted,"clean hysteresis enforced");
  }
  f->measurement_update_mag_only(m);check(f->magneticHealth().accepted,"normal field recovers");
  for(int i=0;i<50;++i)f->measurement_update_mag_only(m);
  check(quat(*f).angularDistance(before)<1e-4f,"recovered heading retained");
#if MAG_FILTER_KIND==5
  const auto covariance=f->covariance();
#else
  const auto covariance=f->covariance_full();
#endif
  check(covariance.allFinite() && (covariance-covariance.transpose()).norm()<1e-6f,"Joseph covariance stays finite and symmetric");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigen(covariance);
  check(eigen.info()==Eigen::Success && eigen.eigenvalues().minCoeff()>-1e-6f,"Joseph covariance stays PSD");
  std::printf("magnetic fusion kind=%d heading_change=%.8g rejected=%u\n",MAG_FILTER_KIND,
      double(quat(*f).angularDistance(before)),f->magneticHealth().rejected);
  return errors?1:0;
}
