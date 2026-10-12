#define EIGEN_NON_ARDUINO
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#include <cstdio>
#include <cstdlib>
#include <memory>

using V = Eigen::Vector3d;
using Q = Eigen::Quaterniond;
static void require(bool ok, const char* reason) {
    if (!ok) { std::fprintf(stderr,"TFG recovery: %s\n",reason); std::exit(1); }
}
template<class F> static typename F::MatrixNX physicalJacobian(const F& f) {
    typename F::MatrixNX J;
    // Differentiate the actual production retraction. No second copy of the
    // analytic covariance-coordinate transform is used as the test oracle.
    for (int col=0; col<F::NX; ++col) {
        auto plus=std::make_unique<F>(f), minus=std::make_unique<F>(f);
        typename F::Tangent d=F::Tangent::Zero(); d(col)=1e-6;
        plus->inject(d); minus->inject(-d);
        const Eigen::AngleAxisd turn(plus->quaternion()*minus->quaternion().conjugate());
        J.template block<3,1>(0,col)=turn.axis()*turn.angle()/2e-6;
        J.template block<3,1>(3,col)=(plus->gyroscope_bias()-minus->gyroscope_bias())/2e-6;
        J.template block<3,1>(6,col)=(plus->get_velocity()-minus->get_velocity())/2e-6;
        J.template block<3,1>(9,col)=(plus->get_position()-minus->get_position())/2e-6;
        J.template block<3,1>(12,col)=(plus->get_integral_displacement()-minus->get_integral_displacement())/2e-6;
        J.template block<3,1>(15,col)=(plus->get_world_accel()-minus->get_world_accel())/2e-6;
        J.template block<3,1>(18,col)=(plus->get_acc_bias()-minus->get_acc_bias())/2e-6;
    }
    return J;
}
template<bool two_frame> static void run() {
    using F=ocean_imu::kalman::Kalman3D_Wave_TFG<double,true,true,two_frame>;
    auto core=std::make_unique<F>(); auto& f=*core;
    f.initialize_identity(.02);
    f.initialize_from_truth(Q(Eigen::AngleAxisd(.7,V(1,2,3).normalized())),
        V(1,-2,3),V(4,.5,-2),V(3,-4,2),V(.4,-.6,.2),V(.2,-.1,.05),V(.1,.2,-.1));
    f.time_update(V(.3,-.2,.4),.01);
    auto physical=[&]() -> typename F::MatrixNX {
        const auto J=physicalJacobian(f);
        return J*f.covariance_full()*J.transpose();
    };
    const auto initial=physical();
    const V p=f.get_position(),v=f.get_velocity(),s=f.get_integral_displacement(),a=f.get_world_accel(),ba=f.get_acc_bias();
    const Q target(Eigen::AngleAxisd(2.1,V(-2,1,3).normalized()));
    require(f.reacquire_attitude(target),"capture rejected");
    const auto captured=physical();
    require((captured.template bottomRightCorner<18,18>()-initial.template bottomRightCorner<18,18>()).norm()<1e-8,"capture changed physical navigation/bias covariance");
    require(captured.template block<3,18>(0,3).norm()<1e-8,"capture retained physical attitude correlations");
    require((p-f.get_position()).norm()+(v-f.get_velocity()).norm()+(s-f.get_integral_displacement()).norm()+(a-f.get_world_accel()).norm()<1e-12,"capture changed navigation means");
    require(f.recover_corrupted_biases(V(.01,-.02,.005)),"bias repair rejected");
    require((ba-f.get_acc_bias()).norm()<1e-12,"gyro repair changed acceleration bias");
    require((f.gyroscope_bias_covariance_body()-.0009*Eigen::Matrix3d::Identity()).norm()<1e-10,"bias covariance has wrong frame");
    const V repaired_bg=f.gyroscope_bias(), before_s=f.get_integral_displacement();
    require((f.get_acc_bias()-ba).norm()<1e-12,"recovery uncertainty changed bias mean");
    const auto repaired=physical();
    require((repaired.template block<3,3>(18,18)-captured.template block<3,3>(18,18)).norm()<1e-8,"bias uncertainty has wrong physical frame");
    require(f.applyIntegralZeroPseudoMeas(),"S correction rejected");
    require((f.gyroscope_bias()-repaired_bg).norm()+(f.get_acc_bias()-ba).norm()<1e-12,"S relearned repaired biases");
    require((f.get_integral_displacement()-before_s).norm()>1e-3,"S navigation correction was discarded");
    require(f.quaternion().angularDistance(target)<1e-10,"retained S spuriously rotated captured attitude");
    f.set_magnetic_reference_world(V(20,0,45));
    const Q observed=Q(Eigen::AngleAxisd(.02,V::UnitZ()))*target;
    require(f.measurement_update_mag_only(observed.conjugate()*V(20,0,45)),"direct magnetic correction rejected after repair");
    require(f.quaternion().angularDistance(target)>1e-6,"direct magnetic aiding remained held");
    f.time_update(V(.05,.03,-.02),.01);
    f.set_attitude_aiding_enabled(false);
    const auto before=physical();
    f.set_attitude_aiding_enabled(true);
    const auto after=physical();
    require((before.diagonal()-after.diagonal()).norm()<1e-8,"resume changed physical marginals");
    require(after.template block<6,12>(0,6).norm()<3e-9*std::max(1.0,after.norm()),"resume retained physical navigation correlations");
    Eigen::SelfAdjointEigenSolver<typename F::MatrixNX> solver(f.covariance_full());
    require(solver.eigenvalues().minCoeff()>-1e-10,"recovery covariance is not PSD");
    f.initialize_identity(.02);
    f.initialize_from_truth(target,v,p,s,a,V(.2,-.1,.05),V(.7,-.3,.2));
    auto ordinary=std::make_unique<F>();
    ordinary->initialize_identity(.02);
    ordinary->initialize_from_truth(target,v,p,s,a,V(.2,-.1,.05),V(.7,-.3,.2));
    f.time_update(V(.1,-.2,.3),.01);
    ordinary->time_update(V(.1,-.2,.3),.01);
    require(f.applyIntegralZeroPseudoMeas() && ordinary->applyIntegralZeroPseudoMeas(),"reinitialized S correction rejected");
    require((f.covariance_full()-ordinary->covariance_full()).norm()<1e-12,"initialization did not restore ordinary S gains");
    const V retained_p=f.get_position(),retained_v=f.get_velocity(),retained_s=f.get_integral_displacement(),retained_a=f.get_world_accel();
    require(f.recover_corrupted_biases(V(.01,-.02,.005)),"gross bias repair rejected");
    require(f.get_acc_bias().norm()==0,"gross acceleration bias did not return to calibration prior");
    require((f.get_position()-retained_p).norm()+(f.get_velocity()-retained_v).norm()+(f.get_integral_displacement()-retained_s).norm()+(f.get_world_accel()-retained_a).norm()<1e-12,"bias repair reset navigation means");
}
int main() { run<false>(); run<true>(); std::puts("TFG recovery physical-coordinate geometry PASS"); }
