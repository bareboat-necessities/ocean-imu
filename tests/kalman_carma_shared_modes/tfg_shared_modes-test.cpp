// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Independent numerical references for the experimental 27-state geometry.
#define EIGEN_NON_ARDUINO
#define EIGEN_RUNTIME_NO_MALLOC
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstdio>
#include <random>
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
namespace cs = ocean_imu::carma_shared_modes;
namespace lie = ocean_imu::lie;
using F = cs::KalmanTfgCarmaSharedModes<double>;
using G = F::Group;
using V = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using M6 = Eigen::Matrix<double,6,6>;
using MX = F::MatrixNX;
using TX = F::Tangent;
constexpr double grav = 9.80665;
int fails=0;
void check(const char* name,double v,double lim) {
    bool ok=std::isfinite(v)&&v<=lim; fails+=!ok;
    std::printf("CHECK %-54s value=%.9g limit=%.3g %s\n",name,v,lim,ok?"PASS":"FAIL");
}
F make_filter(bool coupled=true) {
    F f(V::Constant(.03),V::Constant(.0003),V::Constant(.8),.0005,.000001,1e-8,1.5,grav);
    F::Config c; c.mode=coupled?cs::CouplingMode::External:cs::CouplingMode::Off;
    c.confidence.enabled=false; f.configureSharedModes(c);
    F::Hermitian S; S.re.setIdentity();S.im.setZero();
    S.im(0,2)=.45;S.im(2,0)=-.45;S.re(1,2)=S.re(2,1)=.25;
    f.setExternalCoupling(S);f.set_vertical_carma(.8,.6); f.set_horizontal_carma(.8,.5,.4);
    f.reset_aw_covariance_to_stationary();
    f.set_accel_bias_temp_coeff(V::Zero());f.set_acc_bias_time_constant(800.);
    f.set_Q_bgyro_rw(V::Constant(1e-8));f.set_Q_bacc_rw(V::Constant(5e-4));
    f.set_mag_world_ref(V(18,1.5,43));
    return f;
}
G state() {
    G x;x.R=lie::Exp<double>(V(.2,-.3,.9));
    for(int i=0;i<6;++i)x.X.col(i)=V(.31+i*.21,-.4+i*.17,.5-i*.19);
    x.B.col(0)=V(.001,-.002,.003);x.B.col(1)=V(.01,-.02,.005);return x;
}
G flow(G x,const M6& A,const V& wm,double h) {
    // Independent Taylor exponential, intentionally not the production discretizer.
    M6 E=M6::Identity(),term=E;
    for(int n=1;n<=60;++n){term=term*A*(h/n);E+=term;}
    x.R=x.R*lie::Exp<double>(V((wm-x.B.col(0))*h));x.X=x.X*E.transpose();
    x.B.col(1)*=std::exp(-h/800.);return x;
}
void group_and_flow() {
    F f=make_filter(); f.state()=state();const G x=f.state();const V wm(.8,-.4,1.1);
    double inv=0,rt=0,assoc=0;
    for(int k=1;k<=20;++k){TX a=TX::LinSpaced(27,-.1*k/20,.2*k/20);G e=G::Exp_G(a);
        inv=std::max(inv,G::Log_G(e.compose(e.inverse())).norm());
        rt=std::max(rt,(G::Log_G(e)-a).norm());
        assoc=std::max(assoc,e.compose(x).compose(e).local(e.compose(x.compose(e))).norm());}
    check("27-state group inverse",inv,1e-12);check("27-state exp/log",rt,1e-12);check("27-state associativity",assoc,1e-12);
    for(double h:{.001,.00498,.005,.00501,.023,.2}) {
        MX P,J;f.build_transition(wm,h,P);const G y=flow(x,f.motion().block_A(),wm,h);
        for(int j=0;j<27;++j){TX d=TX::Zero();d(j)=1e-6;
            J.col(j)=(y.local(flow(G::Exp_G(TX(-d)).compose(x),f.motion().block_A(),wm,h))-
                      y.local(flow(G::Exp_G(d).compose(x),f.motion().block_A(),wm,h)))/(2e-6);}
        check("prediction Jacobian vs nonlinear finite difference",(J-P).cwiseAbs().maxCoeff(),2e-8);
        F propagated=f;MX Q;f.build_process_noise(wm,h,Q);
        const MX prior=MX::Identity()*.01;propagated.covariance_full()=prior;
        const MX expected=P*prior*P.transpose()+Q;
        propagated.time_update(wm,h);
        check("structured 27-state covariance vs dense reference",(propagated.covariance_full()-expected).norm()/expected.norm(),2e-11);
        check("cached motion prediction vs independent exponential",propagated.state().local(y).norm(),2e-12);
    }
    TX a=TX::LinSpaced(27,-.2,.3); MX J,Jfd; f.build_reset_jacobian(a,J);
    for(int j=0;j<27;++j){TX ap=a,am=a;ap(j)+=1e-6;am(j)-=1e-6;
        G em=G::Exp_G(TX(-a));Jfd.col(j)=(G::Log_G(em.compose(G::Exp_G(ap)))-G::Log_G(em.compose(G::Exp_G(am))))/(2e-6);}
    check("full group reset Jacobian vs finite difference",(J-Jfd).cwiseAbs().maxCoeff(),2e-9);
    f.inject(a);check("left correction Exp(-a) X",f.state().local(G::Exp_G(TX(-a)).compose(x)).norm(),1e-12);
}
void measurement() {
    F f=make_filter();G x=state();const V g(0,0,grav),bw(18,1.5,43);
    Eigen::Matrix<double,3,27> Ha=Eigen::Matrix<double,3,27>::Zero(),Hm=Ha,Ja,Jm;
    Ha.block<3,3>(0,0)=lie::skew<double>(g);Ha.block<3,3>(0,F::OFF_A)=-M3::Identity();
    Ha.block<3,3>(0,F::OFF_BA)=-M3::Identity();Hm.block<3,3>(0,0)=-lie::skew<double>(bw);
    auto ra=[&](const G& t){return V(x.R*(t.R.transpose()*(t.X.col(2)-g)+t.B.col(1)-x.B.col(1))-(x.X.col(2)-g));};
    auto rm=[&](const G& t){return V(x.R*t.R.transpose()*bw-bw);};
    for(int j=0;j<27;++j){TX d=TX::Zero();d(j)=1e-6;G xp=G::Exp_G(TX(-d)).compose(x),xm=G::Exp_G(d).compose(x);
        Ja.col(j)=(ra(xp)-ra(xm))/(2e-6);Jm.col(j)=(rm(xp)-rm(xm))/(2e-6);}
    check("accelerometer Jacobian vs finite difference",(Ha-Ja).cwiseAbs().maxCoeff(),3e-9);
    check("magnetometer Jacobian vs finite difference",(Hm-Jm).cwiseAbs().maxCoeff(),1e-8);
    // Dense Joseph + finite-difference reset reference of the actual update.
    f.state()=x;MX P=MX::Identity()*.001;f.covariance_full()=P;
    const V meas=x.R.transpose()*(V(x.X.col(2))-g)+x.B.col(1)+V(.03,-.02,.01);
    const V r=x.R*(meas-x.B.col(1))-(x.X.col(2)-g);const M3 R=x.R*(M3::Identity()*.0009)*x.R.transpose();
    const M3 S=Ha*P*Ha.transpose()+R;const Eigen::Matrix<double,27,3>K=P*Ha.transpose()*S.inverse();
    const TX a=K*r;MX Jr;f.build_reset_jacobian(a,Jr);const MX IKH=MX::Identity()-K*Ha;
    const MX expected=Jr*(IKH*P*IKH.transpose()+K*R*K.transpose())*Jr.transpose();
    f.measurement_update_acc_only(meas,35.);
    check("Joseph and covariance transport vs dense reference",(f.covariance_full()-expected).norm()/expected.norm(),2e-12);
    check("measurement group injection",f.state().local(G::Exp_G(TX(-a)).compose(x)).norm(),1e-12);
}
// Continuous invariant model, independent Lyapunov integration. World motion
// at each RK stage comes from its own Taylor exponential, not the filter.
void derivative(const F& f,const G& x,const V& wm,double t,const MX& P,const MX& Q,MX& dP,MX& dQ) {
    const G y=flow(x,f.motion().block_A(),wm,t);const M3 R=y.R;const V w=wm-x.B.col(0);
    MX A=MX::Zero(),L=MX::Zero();A.block<3,3>(0,3)=-M3::Identity();
    A.block<3,3>(3,3)=lie::skew<double>(V(R*w));
    A.block<3,3>(24,24)=lie::skew<double>(V(R*w))-M3::Identity()/800.;
    Eigen::Matrix<double,27,3> Ng=Eigen::Matrix<double,27,3>::Zero(),Nb=Ng,Na=Ng;
    Ng.topRows<3>()=R;Nb.block<3,3>(3,0)=-R;Na.bottomRows<3>()=-R;
    for(int i=0;i<6;++i){const M3 s=lie::skew<double>(V(y.X.col(i)));A.block<3,3>(6+3*i,3)=-s;Ng.block<3,3>(6+3*i,0)=s*R;
        for(int j=0;j<6;++j){A.block<3,3>(6+3*i,6+3*j)=M3::Identity()*f.motion().block_A()(i,j);
            L.block<3,3>(6+3*i,6+3*j)=f.motion().block_WR()(i,j)*f.motion().excitation().re+f.motion().block_WI()(i,j)*f.motion().excitation().im;}}
    L+=Ng*(M3::Identity()*9e-8)*Ng.transpose()+Nb*(M3::Identity()*1e-8)*Nb.transpose()+Na*(M3::Identity()*2.5e-7)*Na.transpose();
    dP=A*P;dQ=A*Q+Q*A.transpose()+L;
}
void covariance_reference() {
    F f=make_filter();f.state()=state();const G x=f.state();const V wm(.8,-.4,1.1);
    for(double h:{.005,.04,.15}){
        MX P=MX::Identity(),Q=MX::Zero();constexpr int n=160;const double dt=h/n;
        for(int k=0;k<n;++k){MX p1,q1,p2,q2,p3,q3,p4,q4;
            derivative(f,x,wm,k*dt,P,Q,p1,q1);derivative(f,x,wm,(k+.5)*dt,P+.5*dt*p1,Q+.5*dt*q1,p2,q2);
            derivative(f,x,wm,(k+.5)*dt,P+.5*dt*p2,Q+.5*dt*q2,p3,q3);derivative(f,x,wm,(k+1)*dt,P+dt*p3,Q+dt*q3,p4,q4);
            P+=dt/6*(p1+2*p2+2*p3+p4);Q+=dt/6*(q1+2*q2+2*q3+q4);}
        MX Pf,Qf;f.build_transition(wm,h,Pf);f.build_process_noise(wm,h,Qf);
        check("transition vs independent continuous RK4",(P-Pf).norm()/P.norm(),2e-9);
        check("Hermitian Qd incl gyro/bias correlations vs RK4",(Q-Qf).norm()/Q.norm(),2e-8);
        Eigen::SelfAdjointEigenSolver<MX> es(Qf);check("process covariance PSD",std::max(0.,-es.eigenvalues().minCoeff()),1e-12);
    }
    M6 E,Qr,Qi;f.motion().discretize(.005,E,Qr,Qi);
    check("quadrature process component is retained",Qi.norm()>1e-7?0.:1.,0.);
    check("imaginary covariance is antisymmetric",(Qi+Qi.transpose()).norm(),1e-14);
}
void recoordinate() {
    F f=make_filter();f.state()=state();const G x=f.state();MX P=MX::Identity()*.01;P(0,2)=P(2,0)=.002;f.covariance_full()=P;
    M6 maps[3];for(int ax=0;ax<3;++ax){maps[ax].setIdentity();maps[ax](4,4)=.5+ax*.6;maps[ax](5,5)=.5+ax*.6;}
    auto change=[&](G y){for(int ax=0;ax<3;++ax)y.X.row(ax)=(maps[ax]*y.X.row(ax).transpose()).transpose();return y;};
    const G y=change(x);MX J;for(int j=0;j<27;++j){TX d=TX::Zero();d(j)=1e-6;
        J.col(j)=(y.local(change(G::Exp_G(TX(-d)).compose(x)))-y.local(change(G::Exp_G(d).compose(x))))/(2e-6);}
    f.apply_motion_recoordinatization(maps);const MX ref=J*P*J.transpose();
    check("anisotropic latent map attitude shear is nonzero",J.block<6,3>(18,0).norm()>.01?0.:1.,0.);
    check("adaptive coordinate covariance vs nonlinear FD",(f.covariance_full()-ref).norm()/ref.norm(),2e-9);
    check("adaptive coordinate nominal state",f.state().local(y).norm(),1e-12);
    // Fall back in invariant coordinates; check full PSD and unchanged common block.
    MX A=MX::Random();P=A*A.transpose()+MX::Identity();f.covariance_full()=P;
    check("Schur fallback accepted",f.shrink_motion_vertical_coupling(.0)?0.:1.,0.);
    Eigen::SelfAdjointEigenSolver<MX> es(f.covariance_full());check("Schur fallback full covariance PSD",std::max(0.,-es.eigenvalues().minCoeff()),1e-11);
    check("Schur fallback keeps common-state covariance",(f.covariance_full().topLeftCorner<6,6>()-P.topLeftCorner<6,6>()).norm(),1e-14);
}
void equivariance() {
    F f=make_filter();
    f.set_horizontal_carma(.8,.6,.6);f.reset_aw_covariance_to_stationary();
    F q=f;const M3 U=lie::Exp<double>(V(0,0,1.2));
    auto S=f.motion().excitation();S.re/=.36;S.im/=.36;
    S.re=(U*S.re*U.transpose()).eval();S.im=(U*S.im*U.transpose()).eval();
    q.setExternalCoupling(S);q.set_vertical_carma(.8,.6);q.reset_aw_covariance_to_stationary();
    f.state()=state();q.state()=f.state();q.state().R=U*f.state().R;q.state().X=U*f.state().X;
    MX T=MX::Zero();for(Eigen::Index b=0;b<9;++b)T.block<3,3>(3*b,3*b)=U;
    f.covariance_full()=MX::Identity()*.001;q.covariance_full()=T*f.covariance_full()*T.transpose();
    f.set_mag_world_ref(V(18,1.5,43));q.set_mag_world_ref(U*V(18,1.5,43));
    for(int k=0;k<50;++k){f.time_update(V(.1,-.2,.4),.005);q.time_update(V(.1,-.2,.4),.005);}
    check("heading rotation prediction mean",(q.state().X-U*f.state().X).norm(),1e-11);
    check("heading rotation covariance",(q.covariance_full()-T*f.covariance_full()*T.transpose()).norm()/f.covariance_full().norm(),2e-10);
    // Acc-only heading gauge: simultaneous rotation of R,W around gravity gives
    // exactly the same specific force. No absolute heading is manufactured.
    const V g(0,0,grav);check("mag outage yaw gauge specific-force invariance",(q.state().R.transpose()*(q.get_world_accel()-g)-f.state().R.transpose()*(f.get_world_accel()-g)).norm(),1e-11);
}
void independent_and_float() {
    using I=ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3<double,4>;
    I ind(V::Constant(.03),V::Constant(.0003),V::Constant(.8),.0005,.000001,1e-8,1.5,grav);
    F f=make_filter(false);ind.set_vertical_carma(.8,.6);ind.set_horizontal_carma(.8,.5,.4);ind.reset_aw_covariance_to_stationary();
    ind.set_accel_bias_temp_coeff(V::Zero());ind.set_acc_bias_time_constant(800.);ind.set_Q_bgyro_rw(V::Constant(1e-8));ind.set_Q_bacc_rw(V::Constant(5e-4));ind.set_mag_world_ref(V(18,1.5,43));
    // Exact observable covariance and lag covariance of the stationary LTI
    // process, independent of any Gaussian finite-angle reset approximation.
    {
        M6 E,Qr,Qi,PR,PIc;f.motion().discretize(.017,E,Qr,Qi);f.motion().stationary_pair(PR,PIc);
        Eigen::Matrix4d E4,Q4,P4;ind.vertical_model().stationary_covariance(P4);
        ocean_imu::kalman::carma::discretize_lti<double,4>(ind.vertical_model().A(),ind.vertical_model().W(),.017,E4,Q4);
        check("coupling off exact stationary observable CARMA law",(PR.topLeftCorner<4,4>()-P4).norm()/P4.norm(),1e-11);
        check("coupling off exact lag covariance CARMA law",((E*PR).topLeftCorner<4,4>()-E4*P4).norm()/P4.norm(),1e-11);
    }
    double ep=0,ea=0,ec=0;
    {
        auto pf=f.covariance_physical();auto pi=ind.covariance_physical();double dd=0.;
        for(int i=0;i<18;++i)for(int j=0;j<18;++j)dd=std::max(dd,std::abs(pf(i,j)-pi(i,j)));
        std::printf("DIAG initial observable cov diff %.12g, cross latent %.12g, gains %.12g %.12g\n",dd,pf.block<12,6>(6,18).norm(),f.motion().noise_gain2(),ind.vertical_model().noise_gain2());
    }
    for(int k=0;k<3000;++k){V w(.01*std::sin(k*.01),.02*std::cos(k*.006),.03);V a(.3*std::sin(k*.005),.1*std::cos(k*.007),-grav+.2*std::sin(k*.006));
        f.time_update(w,.005);ind.time_update(w,.005);f.measurement_update_acc_only(a,35.);ind.measurement_update_acc_only(a,35.);
        if(k%2==0){f.measurement_update_mag_only(V(18,1.5,43));ind.measurement_update_mag_only(V(18,1.5,43));}
        if(k==0||k==1||k==10)std::printf("DIAG step %d pos %.12g rot %.12g latent %.12g\n",k,(f.get_position()-ind.get_position()).norm(),(f.state().R-ind.state().R).norm(),f.state().X.rightCols<2>().norm());
        ep=std::max(ep,(f.get_position()-ind.get_position()).norm());ea=std::max(ea,(f.state().R-ind.state().R).norm());
        ec=std::max(ec,(f.physical_pv_covariance()-ind.physical_pv_covariance()).norm());}
    // These are deliberately reported, not asserted to be round-off: finite
    // group resets also transport the unobserved coloring-state covariance.
    // Its Gaussian tangent approximation is dimension dependent at second
    // order. The paired full-wrapper off ablation quantifies its impact.
    std::printf("REPORT finite-reset off position discrepancy %.12g m\n",ep);
    std::printf("REPORT finite-reset off attitude discrepancy %.12g rad-matrix\n",ea);
    std::printf("REPORT finite-reset off covariance discrepancy %.12g\n",ec);
    cs::KalmanTfgCarmaSharedModes<float> ff(Eigen::Vector3f::Constant(.03f),Eigen::Vector3f::Constant(.0003f),Eigen::Vector3f::Constant(.8f));
    decltype(ff)::Config c;c.mode=cs::CouplingMode::External;ff.configureSharedModes(c);
    decltype(ff)::Hermitian source;source.re.setIdentity();source.im.setZero();
    source.im(0,2)=.3f;source.im(2,0)=-.3f;ff.setExternalCoupling(source);
    double asym=0,min_eig=1.;int bad=0;
    // Eigen's run-time allocator guard covers all real-time paths, including
    // confidence fallback, variable dt, model commits and magnetic corrections.
    Eigen::internal::set_is_malloc_allowed(false);
    for(int k=0;k<24000;++k){
        if(k==2000)ff.forceCouplingConfidence(1.f);
        if(k==3000)ff.forceCouplingConfidence(0.f);
        if(k==4000)ff.forceCouplingConfidence(-1.f);
        if(k%20==0){const float s=.1f+.04f*std::sin(k*.001f);ff.set_vertical_carma(.7f+.2f*std::sin(k*.0007f),s);ff.set_horizontal_carma(.7f,.8f*s,.4f*s);}
        ff.time_update(Eigen::Vector3f(.01f,.02f,.05f),.005f*(1.f+.25f*std::sin(k*.013f)));
        ff.measurement_update_acc_only(Eigen::Vector3f(.2f*std::sin(k*.007f),.05f*std::cos(k*.009f),-9.80665f+.15f*std::sin(k*.011f)),35.f);
        if(k%2==0 && (k<8000||k>16000))ff.measurement_update_mag_only(Eigen::Vector3f(18,1.5,43));
        if(k%200==0){const auto P=ff.covariance_full();asym=std::max(asym,double((P-P.transpose()).cwiseAbs().maxCoeff()));
            Eigen::SelfAdjointEigenSolver<decltype(ff)::MatrixNX> es(P);min_eig=std::min(min_eig,double(es.eigenvalues().minCoeff()));bad+=!P.allFinite();}}
    Eigen::internal::set_is_malloc_allowed(true);
    check("allocation-guarded confidence fallback exercised",ff.coupling().shrink_events>0?0.:1.,0.);
    check("float variable-dt adaptive/outage path finite",bad,0);check("float full covariance symmetry",asym,1e-6);check("float full covariance PSD",std::max(0.,-min_eig),2e-7);
}
int main(){static_assert(F::NX==27);group_and_flow();measurement();covariance_reference();recoordinate();equivariance();independent_and_float();std::printf("FAILURES %d\n",fails);return fails?1:0;}
