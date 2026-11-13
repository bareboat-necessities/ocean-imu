#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
// Dense Joseph reference for a held, uncertain bias. Zero innovations isolate
// the covariance operation; the bias mean must remain exactly held.
#include "StandardHeaders.h"
#include <Eigen/Dense>
#define private public
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#undef private
using V=Eigen::Vector3d;using M=Eigen::Matrix3d;
namespace cs=ocean_imu::carma_shared_modes;
M skew(V v){M s;s<<0,-v.z(),v.y(),v.z(),0,-v.x(),-v.y(),v.x(),0;return s;}
void require(bool ok,const char* m){if(!ok)throw std::runtime_error(m);}
template<class F,bool TFG>void run(){
 F f(V::Constant(.03),V::Constant(.0003),V::Constant(.8));typename F::Config cfg;cfg.mode=cs::CouplingMode::Off;f.configureSharedModes(cfg);
 f.initialize_from_attitude(Eigen::Quaterniond(Eigen::AngleAxisd(.4,V(.2,.5,.8).normalized())),.03,.1);
 f.set_schmidt_acc_bias_hold(true);f.set_acc_bias_updates_enabled(false);f.set_initial_acc_bias(V(.02,-.01,.03));f.set_accel_bias_temp_coeff(V::Zero());f.set_mag_world_ref(V(18,2,43));
 using PType=typename F::MatrixNX;using HType=Eigen::Matrix<double,3,27>;
 auto cov=[&]() -> PType& {if constexpr(TFG)return f.P_;else return f.Pext;};
 PType C=PType::Identity();for(int i=0;i<27;++i)for(int j=0;j<27;++j)C(i,j)+=.03*std::sin(double(i+3*j));cov()=C*C.transpose()*.001;
 const M R=f.quaternion_boat().toRotationMatrix();const V g(0,0,9.80665),ba=f.get_acc_bias();
 const V force=R.transpose()*(f.get_world_accel()-g)+ba;
 HType H=HType::Zero();
 if constexpr(TFG){H.template block<3,3>(0,0)=skew(g);H.template block<3,3>(0,F::OFF_A)=-M::Identity();H.template block<3,3>(0,F::OFF_BA)=-M::Identity();}
 else{H.template block<3,3>(0,0)=-skew(V(force-ba));for(int ax=0;ax<3;++ax){int idx[4];double w[4];int n=f.motion().accel_terms(ax,idx,w);for(int j=0;j<n;++j)H.col(F::OFF_L+idx[j])=R.transpose().col(ax)*w[j];}H.template block<3,3>(0,F::OFF_BA).setIdentity();}
 auto expected=[&](const HType& h,const M& noise){const PType P=cov();Eigen::Matrix<double,27,3> PH=P*h.transpose(),K=PH*(h*PH+noise).inverse();K.template block<3,3>(F::OFF_BA,0).setZero();const PType I=PType::Identity()-K*h;return PType(I*P*I.transpose()+K*noise*K.transpose());};
 M noise=M::Identity()*.0009;if constexpr(TFG)noise=R*noise*R.transpose();
 const PType pa=expected(H,noise);f.measurement_update_acc_only(force,35.);
 double error=(cov()-pa).norm();require(error<1e-10,"Schmidt accel Joseph mismatch");require((ba-f.get_acc_bias()).norm()<1e-13,"held bias mean changed");
 const PType after=cov();f.set_acc_bias_updates_enabled(true);f.set_acc_bias_updates_enabled(false);require((after-cov()).norm()<1e-12,"hold erased cross covariance");
 H.setZero();if constexpr(TFG)H.template block<3,3>(0,0)=-skew(V(18,2,43));else H.template block<3,3>(0,0)=-skew(V(R.transpose()*V(18,2,43)));
 noise=M::Identity()*.64;if constexpr(TFG)noise=R*noise*R.transpose();const PType pm=expected(H,noise);f.measurement_update_mag_only(R.transpose()*V(18,2,43));error=std::max(error,(cov()-pm).norm());require(error<1e-10,"Schmidt mag Joseph mismatch");require((ba-f.get_acc_bias()).norm()<1e-13,"mag correction changed held bias mean");
 f.measurement_update_acc_only(force+V(.01,-.02,.015),35.);
 require((ba-f.get_acc_bias()).norm()<1e-13,"finite correction changed held bias");
 require(Eigen::SelfAdjointEigenSolver<PType>(cov()).eigenvalues().minCoeff()>-1e-10,"finite correction lost PSD");
 f.set_acc_bias_updates_enabled(true);F ordinary=f;ordinary.set_schmidt_acc_bias_hold(false);
 f.time_update(V(.2,-.1,.3),.005);ordinary.time_update(V(.2,-.1,.3),.005);
 f.measurement_update_acc_only(force,35.);ordinary.measurement_update_acc_only(force,35.);
 f.measurement_update_mag_only(R.transpose()*V(18,2,43));ordinary.measurement_update_mag_only(R.transpose()*V(18,2,43));
 require((f.covariance_full()-ordinary.covariance_full()).norm()<1e-12,"Schmidt policy changed active learning");
 require((f.get_acc_bias()-ordinary.get_acc_bias()).norm()<1e-12,"active bias estimate changed");
 std::cout<<(TFG?"TFG":"QMEKF")<<" Schmidt covariance vs dense Joseph error="<<error<<" PASS\n";
}
int main(){run<cs::KalmanCarmaSharedModes<double>,false>();run<cs::KalmanTfgCarmaSharedModes<double>,true>();}
