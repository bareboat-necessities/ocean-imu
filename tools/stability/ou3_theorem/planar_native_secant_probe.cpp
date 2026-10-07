// FINITE DIAGNOSTIC ONLY. Fork the complete inherited shipping execution.
// The only writes are explicit local perturbations at the diagnostic root.
// No modified estimator header, coefficient tape or surrogate is used.
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
#include <string>
#include <vector>
#include <sys/wait.h>
#include <unistd.h>
#define private public
#define protected public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef protected
#undef private
const float g_std=9.80665f;
using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;
using M21=Eigen::Matrix<double,21,21>;
using V21=Eigen::Matrix<double,21,1>;
constexpr int root=40000, end=44000;
static void step(Fusion& f,int k) {
    const double nu=std::acos(-1.)/10,t=k*.005,pitch=.02*std::sin(nu*t);
    const double rate=.02*nu*std::cos(nu*t),ax=-.02*nu*nu*std::sin(nu*t);
    const Eigen::Matrix3d U=Eigen::AngleAxisd(-pitch,Eigen::Vector3d::UnitY()).toRotationMatrix();
    const Eigen::Vector3f a=(U*Eigen::Vector3d(ax,0,-9.80665*(39999./40001))).cast<float>();
    f.update(.005f,Eigen::Vector3f(0,float(rate),0),a,35);
    if(k%8==0)f.updateMag((75*U*Eigen::Vector3d::UnitX()).cast<float>());
}
static std::vector<float> state(const Fusion& f) {
    const auto&m=f.raw().mekf(); std::vector<float> a;
    for(int i=0;i<21;++i)a.push_back(m.xext(i));
    a.insert(a.end(),{m.qref.w(),m.qref.x(),m.qref.y(),m.qref.z()});
    for(int i=0;i<21;++i)for(int j=0;j<21;++j)a.push_back(m.Pext(i,j));
    for(float v:a)if(!std::isfinite(v))throw std::runtime_error("nonfinite state");
    return a;
}
static void hash_float(std::uint64_t& hash,float value) {
    const auto*p=reinterpret_cast<const unsigned char*>(&value);
    for(int j=0;j<4;++j){hash^=p[j];hash*=1099511628211ULL;}
}
static void hash_forcing(const Fusion& f,std::uint64_t& hash) {
    const auto&m=f.raw().mekf(); const auto&r=f.raw();
    for(float v:{m.tau_aw,m.pseudo_update_elapsed_s_,m.pseudo_update_period_s_,
        r.tune_.tau_applied,r.tune_.sigma_applied,r.tune_.RS_applied,
        r.getTauTarget(),r.getSigmaTarget(),r.getRSTarget(),r.getWavePeriodSec(),
        r.getAccelVariance(),r.getAccelVertical(),float(m.aw_covariance_floor_pending_),
        float(m.acc_bias_updates_enabled()),float(f.isLive()),float(f.hasRefinedMagReference())})hash_float(hash,v);
    for(int i=0;i<3;++i){hash_float(hash,m.v2ref(i));
        for(int j=0;j<3;++j)for(float v:{m.R_S(i,j),m.Sigma_aw_stat(i,j),m.Racc(i,j),m.Rmag(i,j)})hash_float(hash,v);}
}
static void finish(Fusion& f,const std::string& path) {
    const auto initial=state(f);std::uint64_t hash=14695981039346656037ULL;std::uint32_t mags=0;
    for(int k=root+1;k<=end;++k){step(f,k);hash_forcing(f,hash);
        if(k%8==0&&f.raw().mekf().lastMagDiag().accepted)++mags;}
    const auto terminal=state(f);
    std::ofstream out(path,std::ios::binary);
    out.write(reinterpret_cast<const char*>(initial.data()),initial.size()*sizeof(float));
    out.write(reinterpret_cast<const char*>(terminal.data()),terminal.size()*sizeof(float));
    out.write(reinterpret_cast<const char*>(&hash),sizeof(hash));
    out.write(reinterpret_cast<const char*>(&mags),sizeof(mags));
    if(!out)throw std::runtime_error("secant output write failed");
}
int main(int argc,char**argv) {
    if(argc!=3)throw std::runtime_error("usage: probe epsilon output-directory");
    const double eps=std::stod(argv[1]);const std::string dir=argv[2];
    if(!(eps>0&&eps<.1))throw std::runtime_error("epsilon outside diagnostic range");
    Fusion::Config cfg;cfg.sigma_a.setConstant(.2f);cfg.sigma_g.setConstant(.00135f);cfg.sigma_m.setConstant(.8f);
    Fusion f;f.begin(cfg);for(int k=1;k<=root;++k)step(f,k);
    auto&m=f.raw().mekf();
    if(!f.isLive()||!f.hasRefinedMagReference()||!m.acc_bias_updates_enabled())throw std::runtime_error("inherited A21 root missing");
    const M21 P=m.Pext.cast<double>();Eigen::LLT<M21> chol(P);
    if(chol.info()!=Eigen::Success)throw std::runtime_error("root P not SPD");
    const M21 L=chol.matrixL();
    int active=0;bool failed=false;
    auto reap=[&](){int status=0;if(wait(&status)<0||!WIFEXITED(status)||WEXITSTATUS(status)!=0)failed=true;--active;};
    // One null perturbation is a bitwise process-copy control.
    for(int column=-1;column<252;++column)for(int sign:{-1,1}) {
        if(column==-1&&sign==1)continue;
        while(active>=4)reap();if(failed)throw std::runtime_error("child execution failed");
        const pid_t child=fork();if(child<0)throw std::runtime_error("fork failed");
        if(child==0){try {
            if(column>=0&&column<21){const V21 dx=sign*eps*L.col(column);
                const Eigen::Vector3d rot=dx.head<3>();
                if(rot.norm()>0){const Eigen::Quaterniond dq(Eigen::AngleAxisd(rot.norm(),rot.normalized()));
                    m.qref=(dq.cast<float>()*m.qref).normalized();}
                for(int i=3;i<21;++i)m.xext(i)+=float(dx(i));
            } else if(column>=21){int n=21;M21 E=M21::Zero();
                for(int i=0;i<21;++i)for(int j=i;j<21;++j,++n)if(n==column){
                    E(i,j)=E(j,i)=i==j?1.:1./std::sqrt(2.);}
                const M21 perturbed=P+sign*eps*L*E*L.transpose();
                m.Pext=((perturbed+perturbed.transpose())*.5).cast<float>();
                Eigen::LLT<M21> check(m.Pext.cast<double>());if(check.info()!=Eigen::Success)throw std::runtime_error("perturbed P not SPD");
            }
            finish(f,dir+"/"+std::to_string(column)+"_"+std::to_string(sign)+".bin");
            _exit(0);
        }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;_exit(1);}}
        ++active;
    }
    finish(f,dir+"/nominal.bin");while(active>0)reap();
    if(failed)throw std::runtime_error("secant child failed");
    std::cout<<"{\"result_type\":\"FINITE DIAGNOSTIC ONLY\",\"columns\":252,\"word_samples\":[40000,44000]}\n";
}
