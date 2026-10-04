// Read-only literal OU-III planar operation export. No estimator overrides.
// The Python driver inserts observation callbacks into a temporary header only.
#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using M21 = Eigen::Matrix<double,21,21>;
using Probe = Eigen::Matrix<double,21,4>;
static constexpr int EVEN[12]={1,4,6,8,9,11,12,14,15,17,18,20};
static constexpr int ODD[9]={0,2,3,5,7,10,13,16,19};
static bool moving_recording=false, export_enabled=false;
static std::uint32_t sample_index=0;
static std::ofstream stream;
static Probe phi=Probe::Zero();
static Eigen::Matrix4d information=Eigen::Matrix4d::Zero();
static int word_mags=0;
static Eigen::Matrix<float,21,21> before_prediction;
static bool have_prediction_root=false;
static double parity_defect=0;

static void emit(std::uint32_t kind,const std::vector<float>& values) {
    if (!export_enabled) return;
    const std::uint32_t header[3]={kind,sample_index,static_cast<std::uint32_t>(values.size())};
    stream.write(reinterpret_cast<const char*>(header),sizeof(header));
    stream.write(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(float));
    if (!stream) throw std::runtime_error("operation export failed");
}
static void dense(std::vector<float>& out,const float* p,int rows,int cols) {
    for(int i=0;i<rows;++i)for(int j=0;j<cols;++j)out.push_back(p[i+rows*j]);
}
static void parity(std::vector<float>& out,const float* p) {
    for(int i:EVEN)for(int j:EVEN)out.push_back(p[i+21*j]);
    for(int i:ODD)for(int j:ODD)out.push_back(p[i+21*j]);
    for(int i:EVEN)for(int j:ODD){
        parity_defect=std::max(parity_defect,std::abs(double(p[i+21*j])));
        parity_defect=std::max(parity_defect,std::abs(double(p[j+21*i])));
    }
}
static void planar_before_prediction(const float* p) {
    before_prediction=Eigen::Map<const Eigen::Matrix<float,21,21>>(p);
    have_prediction_root=true;
}
static void moving_prediction(const float* fa,const float* fl,const float* qa,
                              const float* ql,float bphi,const float* qb,
                              const float* pminus,const float* rs,
                              float elapsed,float period,float dt) {
    if(!have_prediction_root)throw std::runtime_error("missing pre-prediction covariance");
    Eigen::Matrix<float,21,21> F=Eigen::Matrix<float,21,21>::Zero(),Q=F;
    F.block<6,6>(0,0)=Eigen::Map<const Eigen::Matrix<float,6,6>>(fa);
    F.block<12,12>(6,6)=Eigen::Map<const Eigen::Matrix<float,12,12>>(fl);
    F.block<3,3>(18,18)=bphi*Eigen::Matrix3f::Identity();
    Q.block<6,6>(0,0)=Eigen::Map<const Eigen::Matrix<float,6,6>>(qa);
    Q.block<12,12>(6,6)=Eigen::Map<const Eigen::Matrix<float,12,12>>(ql);
    Q.block<3,3>(18,18)=Eigen::Map<const Eigen::Matrix3f>(qb);
    phi=(F.cast<double>()*phi).eval();
    std::vector<float> out;out.reserve(912);
    parity(out,before_prediction.data());parity(out,F.data());parity(out,Q.data());
    parity(out,pminus);dense(out,rs,3,3);
    out.insert(out.end(),{elapsed,period,dt});emit(1,out);
}
static void planar_after_sync(const float* p,const float* target,bool was_pending) {
    std::vector<float> out;out.reserve(235);parity(out,p);dense(out,target,3,3);
    out.push_back(was_pending?1.f:0.f);emit(6,out);
}
static void moving_correction(const char* kind,const float* h,const float* s,
                              const float* k,const float* p,const float* rnoise,
                              const float* pct,const float* residual) {
    const Eigen::Map<const Eigen::Matrix<float,3,21>> H(h);
    const Eigen::Map<const Eigen::Matrix3f> S(s);
    const Eigen::Map<const Eigen::Matrix<float,21,3>> K(k);
    const Eigen::Matrix<double,3,4> Y=H.cast<double>()*phi;
    if(kind[0]=='m') {
        Eigen::Matrix3d Sd=S.cast<double>();
        Eigen::LLT<Eigen::Matrix3d> factor(Sd);
        if(factor.info()!=Eigen::Success)throw std::runtime_error("magnetic innovation not SPD");
        Eigen::Matrix<double,3,4> W=factor.matrixL().solve(Y);
        information.noalias()+=W.transpose()*W;
        ++word_mags;
    }
    phi.noalias()-=K.cast<double>()*Y;
    std::vector<float> out;out.reserve(435);
    parity(out,p);dense(out,h,3,21);dense(out,rnoise,3,3);dense(out,s,3,3);
    dense(out,k,21,3);dense(out,pct,21,3);dense(out,residual,3,1);
    emit(kind[0]=='a'?2:kind[0]=='m'?3:4,out);
}
static void moving_reset(const float* d,const float* p) {
    const double x=d[0]*.5,y=d[1]*.5,z=d[2]*.5;
    Eigen::Matrix3d G;G<<1,-z,y,z,1,-x,-y,x,1;
    phi.topRows<3>()=(G*phi.topRows<3>()).eval();
    std::vector<float> out;out.reserve(228);parity(out,p);dense(out,d,3,1);emit(5,out);
}
static void planar_after_reset(const float* p) {
    std::vector<float> out;out.reserve(225);parity(out,p);emit(8,out);
}

#define private public
#define protected public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef protected
#undef private
const float g_std=9.80665f;
using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;

int main(int argc,char**argv) {
    if(argc!=4)throw std::runtime_error("usage: probe duration tail_seconds stream");
    const double duration=std::stod(argv[1]),tail=std::stod(argv[2]);
    const int steps=static_cast<int>(duration*200),first_export=steps-static_cast<int>(tail*200)+1;
    stream.open(argv[3],std::ios::binary);
    if(!stream)throw std::runtime_error("cannot open operation export");
    const char magic[8]={'O','U','3','E','V','T','1','\0'};stream.write(magic,8);
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(.2f);cfg.sigma_g.setConstant(.00135f);cfg.sigma_m.setConstant(.8f);
    Fusion f;f.begin(cfg);
    constexpr double g=9.80665,sn=400./40001,cs=39999./40001;
    const double nu=std::acos(-1.)/10;
    int live=-1,refined=-1,active=-1,accepted=0,root=-1,windows=0,min_mags=10000;
    double service_min=1e100,service_max=0;
    for(int k=1;k<=steps;++k) {
        sample_index=k;export_enabled=k>=first_export;
        double t=k*.005,angle=nu*t,pitch=.02*std::sin(angle),rate=.02*nu*std::cos(angle);
        double ax=-.02*nu*nu*std::sin(angle);
        const Eigen::Matrix3d U=Eigen::AngleAxisd(-pitch,Eigen::Vector3d::UnitY()).toRotationMatrix();
        Eigen::Vector3f a=(U*Eigen::Vector3d(ax,0,-g*cs)).cast<float>();
        f.update(.005f,Eigen::Vector3f(0,float(rate),0),a,35);
        if(k%8==0) {
            f.updateMag((75*U*Eigen::Vector3d::UnitX()).cast<float>());
            if(f.isLive()&&f.raw().mekf().lastMagDiag().accepted)++accepted;
        }
        const auto& m=f.raw().mekf();
        if(!m.Pext.allFinite()||!m.xext.allFinite()||!m.qref.coeffs().allFinite())
            throw std::runtime_error("nonfinite shipping state");
        if(live<0&&f.isLive())live=k;
        if(refined<0&&f.hasRefinedMagReference())refined=k;
        if(active<0&&m.acc_bias_updates_enabled())active=k;
        if(moving_recording&&k-root==200) {
            for(int j=0;j<4;j+=2) {
                Eigen::Matrix2d I=information.block<2,2>(j,j);
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(I);
                if(es.info()!=Eigen::Success)throw std::runtime_error("service eigensolver failed");
                double v=es.eigenvalues().minCoeff();
                service_min=std::min(service_min,v);service_max=std::max(service_max,v);
            }
            ++windows;min_mags=std::min(min_mags,word_mags);moving_recording=false;
        }
        if(active>=0&&k>=active+3400&&!moving_recording) {
            phi.setZero();information.setZero();word_mags=0;root=k;
            for(int sign=0;sign<2;++sign) {
                Eigen::Vector3d d=U*Eigen::Vector3d(0,sign?-sn:sn,cs);
                phi.block<3,1>(0,2*sign)=d;
                phi.block<3,1>(3,2*sign+1)=.02*d;
            }
            moving_recording=true;
        }
        if(export_enabled) {
            const auto& raw=f.raw();std::vector<float> out;out.reserve(283);
            parity(out,m.Pext.data());dense(out,m.xext.data(),21,1);
            out.insert(out.end(),{m.qref.w(),m.qref.x(),m.qref.y(),m.qref.z()});
            dense(out,m.v2ref.data(),3,1);dense(out,m.Racc.data(),3,3);dense(out,m.R_S.data(),3,3);
            dense(out,m.Sigma_aw_stat.data(),3,3);
            out.insert(out.end(),{m.tau_aw,m.pseudo_update_elapsed_s_,m.pseudo_update_period_s_,
              raw.tune_.tau_applied,raw.tune_.sigma_applied,raw.tune_.RS_applied,
              raw.getTauTarget(),raw.getSigmaTarget(),raw.getRSTarget(),
              raw.getWavePeriodSec(),raw.getAccelVariance(),raw.getAccelVertical(),
              m.aw_covariance_floor_pending_?1.f:0.f});
            emit(7,out);
        }
    }
    stream.close();
    const auto&m=f.raw().mekf();
    std::cout<<std::setprecision(17)<<"{\"duration\":"<<duration<<",\"tail_seconds\":"<<tail
      <<",\"service_min\":"<<service_min<<",\"service_max\":"<<service_max
      <<",\"service_windows\":"<<windows<<",\"service_mags_min\":"<<min_mags
      <<",\"live_step\":"<<live<<",\"refined_step\":"<<refined<<",\"active_step\":"<<active
      <<",\"accepted_mag\":"<<accepted<<",\"parity_max_abs\":"<<parity_defect
      <<",\"terminal_state\":[";
    bool first=true;
    auto put=[&](double x){if(!first)std::cout<<",";first=false;std::cout<<x;};
    for(int i=0;i<21;++i)put(m.xext(i));
    put(m.qref.w());put(m.qref.x());put(m.qref.y());put(m.qref.z());
    for(int i=0;i<21;++i)for(int j=0;j<21;++j)put(m.Pext(i,j));
    std::cout<<"]}\n";
}
