#include "Adapter.h"
using namespace followup;
template<class Native,bool TFG>void audit() {
 using Core=Probe<Native,TFG>;using H3=typename Core::H3;
 Core m(V::Constant(.03f),V::Constant(.00003f),V::Constant(7.68f));
 m.initialize_from_attitude(Eigen::Quaternionf(Eigen::AngleAxisf(.37f,V(.2f,.7f,.4f).normalized())),.03f,.1f);
 m.set_acc_bias_updates_enabled(true);m.set_initial_acc_bias(V(.04f,-.03f,.02f));m.set_mag_world_ref(V(22,3,49));
 if constexpr(TFG)m.state().X.col(Native::IA)=V(.2f,-.3f,.7f);
 else {H3 W=m.acceleration_jacobian();m.xext+=W.transpose()*(W*W.transpose()).ldlt().solve(V(.2f,-.3f,.7f));}
 const Eigen::Matrix3f R=m.quaternion_boat().toRotationMatrix();const V g(0,0,9.80665f),B=m.applied_field();
 const V f=R.transpose()*(m.get_world_accel()-g)+m.get_acc_bias(),mag=R.transpose()*B;
 H3 Ha=H3::Zero(),Hm=H3::Zero();
 if constexpr(TFG) {Ha.template block<3,3>(0,0)=skew(g);Ha.template block<3,3>(0,Native::OFF_A)=-Eigen::Matrix3f::Identity();Ha.template block<3,3>(0,Native::OFF_BA)=-Eigen::Matrix3f::Identity();Hm.template block<3,3>(0,0)=-skew(B);}
 else {Ha=R.transpose()*m.acceleration_jacobian();Ha.template block<3,3>(0,0)=-skew(V(f-m.get_acc_bias()));Ha.template block<3,3>(0,Native::OFF_BA).setIdentity();Hm.template block<3,3>(0,0)=-skew(mag);}
 float ea=0,em=0;const float eps=.0005f;
 for(int j=0;j<27;++j){Core x=m;
  if constexpr(TFG){typename Native::Tangent d=Native::Tangent::Zero();d(j)=eps;x.state()=m.state().retract(d);}
  else {if(j<3)x.set_quaternion_boat(m.quaternion_boat()*Eigen::Quaternionf(Eigen::AngleAxisf(-eps,V::Unit(j))));else x.xext(j)+=eps;}
  Eigen::Matrix3f Q=x.quaternion_boat().toRotationMatrix();V ra,rm;
  if constexpr(TFG){ra=Q*(f-x.get_acc_bias())-(x.get_world_accel()-g);rm=Q*mag-B;}
  else {ra=Q.transpose()*(x.get_world_accel()-g)+x.get_acc_bias()-f;rm=Q.transpose()*B-mag;}
  ea=std::max(ea,(ra/eps-Ha.col(j)).norm());em=std::max(em,(rm/eps-Hm.col(j)).norm());
 }
 std::cout<<(TFG?"TFG":"QMEKF")<<" native accel/mag Jacobian max column errors "<<ea<<" "<<em<<"\n";
 if(ea>.02||em>.04)throw std::runtime_error("native sensor Jacobian audit failed");
}
int main()try{audit<ocean_imu::carma_shared_modes::KalmanCarmaSharedModes<float>,false>();audit<ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>,true>();}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
