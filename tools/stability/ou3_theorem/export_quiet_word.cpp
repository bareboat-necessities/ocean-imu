// Diagnostic of the unchanged shipping core, not wrapper capture or a uniform proof.
// Quiet-water truth is continuous: p=v=a=omega=bias=0; field is fixed in NED.
// The core's covariance/state/scheduler are carried through warmup and the word.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private

using Core=Kalman3D_Wave_OU_III<double,true,true>;
using V3=Eigen::Vector3d;
using M21=Eigen::Matrix<double,21,21>;
using H3=Eigen::Matrix<double,3,21>;

template<class A> void matrix(const A& a){
    std::cout<<'[';
    for(int i=0;i<a.rows();++i){
        if(i)std::cout<<',';
        std::cout<<'[';
        for(int j=0;j<a.cols();++j){if(j)std::cout<<',';std::cout<<a(i,j);}
        std::cout<<']';
    }
    std::cout<<']';
}

static void correction(const char* channel,const Core& m){
    std::cout<<"{\"kind\":\"correction\",\"channel\":\""<<channel<<"\",\"S_actual\":";
    matrix(m.S_scratch_);std::cout<<",\"gain_actual\":";matrix(m.K_scratch_);
    std::cout<<"}\n";
}

int main(int argc,char** argv){
    const int warmup=argc>1?std::stoi(argv[1]):3200;
    const int samples=argc>2?std::stoi(argv[2]):3200;
    if(warmup<0||samples<=0||warmup+samples>1000000)return 2;
    std::cout<<std::setprecision(17);
    Core m(V3::Constant(.12),V3::Constant(.00135),V3::Constant(.25),
           5e-4,1e-6,1e-10,.15*.15);
    m.set_aw_time_constant(1.8);
    m.set_aw_stationary_std(V3::Constant(.05));
    m.set_RS_noise(V3(.108,.075,.15));
    m.set_pseudo_update_period_s(.15);
    m.set_acc_bias_updates_enabled(true);
    const V3 field(20,0,40),acc(0,0,-9.80665);
    m.set_mag_world_ref(field);
    const double dt=.005;
    H3 hs=H3::Zero(),ha=H3::Zero(),hm=H3::Zero();
    hs.block<3,3>(0,12).setIdentity();
    ha.block<3,3>(0,0)=-m.skew_symmetric_matrix(acc);
    ha.block<3,3>(0,15).setIdentity();ha.block<3,3>(0,18).setIdentity();
    hm.block<3,3>(0,0)=-m.skew_symmetric_matrix(field);
    M21 fixed_f=M21::Zero(),fixed_q=M21::Zero();
    for(int k=0;k<warmup+samples;++k){
        const bool emit=k>=warmup;
        if(k==warmup){
            std::cout<<"{\"kind\":\"root\",\"scope\":\"configured quiet-water shipping core; no wrapper capture, adaptive schedule or uniform claim\",\"warmup_samples\":"<<warmup<<",\"samples\":"<<samples<<",\"dt\":"<<dt<<",\"P\":";
            matrix(m.covariance_full());std::cout<<"}\n";
        }
        double elapsed=m.pseudo_update_elapsed_s_;
        const bool pseudo=ocean_imu::kalman::ou_detail::periodic_update_due(
            dt,m.pseudo_update_period_s_,elapsed);
        m.time_update(V3::Zero(),dt);
        M21 f=M21::Identity(),q=M21::Zero();
        f.block<6,6>(0,0)=m.F_AA_scratch_;q.block<6,6>(0,0)=m.Q_AA_scratch_;
        f.block<12,12>(6,6)=m.F_LL_scratch_;q.block<12,12>(6,6)=m.Q_LL_scratch_;
        f.block<3,3>(18,18)*=std::exp(-dt/m.tau_bacc_);
        q.block<3,3>(18,18)=m.Q_bacc_*(-.5*m.tau_bacc_*std::expm1(-2*dt/m.tau_bacc_));
        if(k==0){fixed_f=f;fixed_q=q;}
        if(!(f.array()==fixed_f.array()).all()||!(q.array()==fixed_q.array()).all())return 3;
        if(k==warmup){
            std::cout<<"{\"kind\":\"model\",\"F\":";matrix(f);
            std::cout<<",\"Q\":";matrix(q);
            std::cout<<",\"S\":{\"H\":";matrix(hs);std::cout<<",\"R\":";matrix(m.R_S);
            std::cout<<"},\"acc\":{\"H\":";matrix(ha);std::cout<<",\"R\":";matrix(m.Racc);
            std::cout<<"},\"mag\":{\"H\":";matrix(hm);std::cout<<",\"R\":";matrix(m.Rmag);
            std::cout<<"}}\n";
        }
        if(emit)std::cout<<"{\"kind\":\"prediction\"}\n";
        if(pseudo&&emit)correction("S",m);
        m.measurement_update_acc_only(acc);
        if(!m.lastAccDiag().accepted)return 4;
        if(emit)correction("acc",m);
        if(k%2==1){
            m.measurement_update_mag_only(field);
            if(!m.lastMagDiag().accepted)return 5;
            if(emit)correction("mag",m);
        }
        // Zero innovations make every implemented quaternion reset identity.
        if(m.xext.norm()!=0||m.qref.angularDistance(Eigen::Quaterniond::Identity())!=0)return 6;
    }
    std::cout<<"{\"kind\":\"endpoint\",\"P_actual\":";
    matrix(m.covariance_full());std::cout<<"}\n";
}
