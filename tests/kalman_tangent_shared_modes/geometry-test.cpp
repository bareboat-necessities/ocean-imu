#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <random>
#include <stdexcept>
#include "../live_bias_consistency/StandardHeaders.h"
#define private public
#include "kalman_carma_shared_modes/KalmanTangentCarmaSharedModes.h"
#undef private
using G=ocean_imu::lie::TangentBiasGroup<double>;
using F=ocean_imu::carma_shared_modes::KalmanTangentCarmaSharedModes<double>;
using V=Eigen::Vector3d;
using M=G::AdjointMat;
void check(bool x,const char* msg) { if(!x) throw std::runtime_error(msg); }
int main() try {
    ocean_imu::carma_shared_modes::SharedModesConfig<float> legacy{
        ocean_imu::carma_shared_modes::CouplingMode::Off};
    check(legacy.mode==ocean_imu::carma_shared_modes::CouplingMode::Off && !legacy.bias_learning.enabled,"legacy aggregate configuration");
    std::mt19937 gen(48031);std::normal_distribution<double> normal;
    auto random=[&](double s){G::Tangent x;for(int i=0;i<27;++i)x(i)=s*normal(gen);return x;};
    static_assert(sizeof(G)==sizeof(G::Base),"geometry adds no physical storage");
    double worst_reset=0,worst_group=0;
    for(int trial=0;trial<25;++trial) {
        const auto a=random(.35),b=random(.2),c=random(.3);
        const G A=G::Exp_G(a),B=G::Exp_G(b),C=G::Exp_G(c);
        worst_group=std::max(worst_group,(G::Log_G(A)-a).norm());
        check(G::Log_G(A*A.inverse()).norm()<1e-10,"inverse");
        check(((A*B)*C).local(A*(B*C)).norm()<1e-10,"associativity");
        check(G::Log_G(A*G::Exp_G(G::Tangent(1e-5*b))*A.inverse()).isApprox(A.Adjoint()*(1e-5*b),1e-5),"adjoint");
        check(G::Exp_G(G::Tangent(.4*a)).compose(G::Exp_G(G::Tangent(.6*a))).local(A).norm()<1e-10,"one parameter subgroup");
        F f(V::Constant(.02),V::Constant(.0002),V::Constant(.8));
        f.state().R=B.state.R;f.state().X=B.state.X;f.state().B=B.state.B;
        f.set_acc_bias_updates_enabled(true);
        M J;f.build_reset_jacobian(a,J);
        const auto at=G::from_two_frame(a,V(B.state.X.col(1)));
        G next=G::Exp_G(G::Tangent(-at))*B;
        const double h=2e-6;
        for(int j=0;j<27;++j) {
            G::Tangent d=G::Tangent::Zero();d(j)=h;
            auto sd=G::from_two_frame(d,V(B.state.X.col(1)));
            // Transport residual tangent coordinates through Exp(-at)Exp(at+sd).
            auto residual=[&](const G::Tangent& e) {
                const G z=G::Exp_G(G::Tangent(-at))*G::Exp_G(e);
                auto y=G::Log_G(z);
                y.template segment<3>(G::OFF_BA)-=V(next.state.X.col(1)).cross(V(y.template segment<3>(G::OFF_BG)));
                return y;
            };
            const G::Tangent num=(residual(G::Tangent(at+sd))-residual(G::Tangent(at-sd)))/(2*h);
            worst_reset=std::max(worst_reset,(num-J.col(j)).norm());
        }
        // The optimized rank-row congruence must match dense reset transport.
        const M ref=f.covariance_full();
        const M expected=J*ref*J.transpose();
        f.transport_(a);
        check((f.covariance_full()-expected).norm()<1e-11*(1+expected.norm()),"structured covariance reset");
        f.P_=ref;
        f.set_acc_bias_updates_enabled(true);
        f.measurement_update_acc_only(V(.2,-.3,-9.7));
        auto P=f.covariance_full();
        check(P.allFinite(),"finite correction");
        Eigen::SelfAdjointEigenSolver<M> es((P+P.transpose())*.5);
        check(es.eigenvalues().minCoeff()>-1e-10,"correction PSD");
        (void)ref;
    }
    check(worst_group<2e-10,"Exp Log round trip");
    check(worst_reset<2e-7,"reset finite difference");
    using Gate=ocean_imu::kalman::filter_detail::BiasLearningGate<double>;
    using Config=ocean_imu::kalman::filter_detail::BiasLearningConfig<double>;
    Gate quiet,rock;Config cfg;cfg.enabled=true;quiet.configure(cfg);rock.configure(cfg);
    for(int k=0;k<3000;++k) {
        const double t=k*.02;
        // Body-to-world ZYX attitude: Ry(pitch) Rx(roll) has zero heading.
        Eigen::Matrix3d R=(Eigen::AngleAxisd(.2*std::sin(.7*t),V::UnitY())*Eigen::AngleAxisd(.3*std::sin(.4*t),V::UnitX())).toRotationMatrix();
        check(std::abs(R(1,0))<1e-14 && R(0,0)>0,"rocking heading remains zero");
        quiet.observe(Eigen::Matrix3d::Identity(),V(20,0,40),V::Zero(),.02,9.80665);
        rock.observe(R,V(20,0,40),V::Zero(),.02,9.80665);
    }
    check(quiet.information_eigenvalues()(0)<1e-7,"quiet weak direction");
    check(rock.information_eigenvalues()(0)>1e-3,"roll/pitch information without yaw");
    // Arbitrary nonoptimal gain must retain its Joseph covariance penalty.
    Eigen::Matrix<double,6,6> P=Eigen::Matrix<double,6,6>::Identity();
    Eigen::Matrix<double,3,6> H=Eigen::Matrix<double,3,6>::Random();
    Eigen::Matrix3d S=H*P*H.transpose()+Eigen::Matrix3d::Identity();
    Eigen::Matrix<double,6,3> K=P*H.transpose()*S.inverse(),L=K;L.bottomRows<3>()=rock.scale()*L.bottomRows<3>();
    auto joseph=[&](const Eigen::Matrix<double,6,3>& gain){Eigen::Matrix<double,6,6> A=Eigen::Matrix<double,6,6>::Identity()-gain*H;return Eigen::Matrix<double,6,6>(A*P*A.transpose()+gain*gain.transpose());};
    check((joseph(L)-joseph(K)-(L-K)*S*(L-K).transpose()).norm()<1e-12,"Joseph gain penalty");
    std::cout<<"group/reset PASS max_group="<<worst_group<<" max_reset="<<worst_reset<<" rocking_information="<<rock.information_eigenvalues().transpose()<<'\n';
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
