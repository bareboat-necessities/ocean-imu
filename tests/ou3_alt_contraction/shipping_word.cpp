// Host-only literal paired-word experiment. Runtime sources are not edited.
// fork() preserves EVERY private field at the same reachable baseline root;
// injected error/covariance probes are explicitly not claimed reachable roots.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/wait.h>
#include <unistd.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

namespace alt {
template<class F> void point(unsigned kind, const F& f);
template<class V> void set_hard_iron(const V& v);
template<class F,class V,class N,class S,class M> void solve(unsigned kind,const F& f,
    const V& r,const N& n,const S& s,const N& k,const S& noise,const M& measured);
template<class F> struct ExitGuard {
    unsigned kind; const F& filter;
    ExitGuard(unsigned k,const F& f):kind(k),filter(f){}
    ~ExitGuard(){point(kind,filter);}
};
}
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std=9.80665f;
using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;
using V3=Eigen::Vector3f;
constexpr float DT=1.0f/200.0f;
constexpr double PI=3.1415926535897932384626433832795029;
constexpr double W=0.7853981633974483096156608458198757;
namespace alt {
std::ofstream trace;
unsigned sample_index=0;
std::array<float,27> source{};
std::array<float,14> frontend{};
std::array<float,3> hard_iron{};
template<class V> void set_hard_iron(const V& v) {
    for(int i=0;i<3;++i)hard_iron[i]=v[i];
}

void write_u(unsigned x) { const std::uint32_t y=x;trace.write(reinterpret_cast<const char*>(&y),4); }
void write_f(float x) {trace.write(reinterpret_cast<const char*>(&x),4);}
template<class M> void matrix(const M& m) {
    for(int i=0;i<m.rows();++i)for(int j=0;j<m.cols();++j)write_f(float(m(i,j)));
}
template<class F,class V,class N,class S,class M> void solve(unsigned kind,const F& f,
    const V& r,const N& n,const S& s,const N& k,const S& noise,const M& measured) {
    if(!trace.is_open())return;
    write_u(kind);write_u(sample_index);write_u(f.acc_bias_updates_enabled()?1:0);
    write_f(f.qref.w());write_f(f.qref.x());write_f(f.qref.y());write_f(f.qref.z());
    matrix(f.xext);matrix(f.Pext);matrix(n);matrix(s);matrix(k);matrix(r);matrix(noise);
    for(float v:source)write_f(v);
    for(float v:frontend)write_f(v);
    matrix(measured);matrix(f.v2ref);
    write_f(f.gravity_magnitude_);write_f(f.get_acc_bias_time_constant());
    write_f(f.acc_bias_limit_);write_f(f.use_imu_lever_arm_?1.0f:0.0f);
    write_f(f.wind_heel_rad_);write_f(f.k_a_.norm());
    for(float v:hard_iron)write_f(v);
}
template<class F> void point(unsigned kind,const F& f) {
    if(!trace.is_open())return;
    Eigen::Matrix<float,21,3> n=Eigen::Matrix<float,21,3>::Zero();
    Eigen::Matrix3f s=Eigen::Matrix3f::Zero();V3 r=V3::Zero();
    solve(kind,f,r,n,s,n,s,r);
}
void front(Fusion& f) {
    auto& r=f.raw();
    set_hard_iron(f.mag_hard_iron_body_uT_);
    frontend={r.getFreqHz(),r.getWavePeriodSec(),r.getTauTarget(),r.getSigmaTarget(),
        r.getRSTarget(),r.getTauApplied(),r.getSigmaApplied(),r.getRSApplied(),
        r.getPseudoUpdatePeriodSec(),r.getAccelVariance(),r.getAccelVertical(),
        float(r.online_tune_apply_pending_),float(r.getStartupStage()),float(r.time_)};
}
}

struct Physical {
    V3 p,v,a,beta,gyro,acc,mag;
    Eigen::Quaternionf q;
    std::array<float,3> angles{};
};
Physical physical(unsigned k,int family) {
    const double t=double(k)*double(DT),s=std::sin(W*t),c=std::cos(W*t);
    Physical z;
    z.p=V3(float(.3*std::sin(W*t+.4)),float(.2*c),float(s));
    z.v=V3(float(.3*W*std::cos(W*t+.4)),float(-.2*W*s),float(W*c));
    z.a=-float(W*W)*z.p;
    const double roll=.12*std::sin(W*t+.2), pitch=.09*std::sin(1.25*W*t+.7),
                 yaw=.08*std::sin(.75*W*t+.6);
    const double rd=.12*W*std::cos(W*t+.2),pd=.09*1.25*W*std::cos(1.25*W*t+.7),
                 yd=.08*.75*W*std::cos(.75*W*t+.6);
    z.q=Eigen::AngleAxisf(float(yaw),V3::UnitZ())*
        Eigen::AngleAxisf(float(pitch),V3::UnitY())*
        Eigen::AngleAxisf(float(roll),V3::UnitX());
    z.gyro=V3(float(rd-yd*std::sin(pitch)),
              float(pd*std::cos(roll)+yd*std::sin(roll)*std::cos(pitch)),
              float(-pd*std::sin(roll)+yd*std::cos(roll)*std::cos(pitch)));
    if(family==0) z.beta=V3(.025f,-.020f,.015f)+
        float(.005*std::exp(-t/900))*V3(1,-.8f,.6f)+
        float(.002*std::sin(2*PI*t/900))*V3(1,.5f,-.5f);
    else if(family==1) z.beta=float(std::exp(-t/1200))*V3(.05f,-.04f,.03f)+
        float(.01*std::sin(2*PI*t/600))*V3(1,-.8f,.6f);
    else z.beta=V3(.04f,-.03f,.02f)+float(.005*std::sin(2*PI*t/900))*V3(1,-.8f,.6f);
    z.acc=z.q.conjugate()*(z.a-V3(0,0,g_std))+z.beta;
    z.mag=z.q.conjugate()*V3(22,0,43);
    z.angles={float(roll),float(pitch),float(yaw)};
    return z;
}
V3 primitive(unsigned k,unsigned origin) {
    const double t=double(k)*double(DT),u=double(origin)*double(DT);
    return V3(float(.3*(std::cos(W*u+.4)-std::cos(W*t+.4))/W),
              float(.2*(std::sin(W*t)-std::sin(W*u))/W),
              float((std::cos(W*u)-std::cos(W*t))/W));
}
void advance(Fusion& f,unsigned k,int family,unsigned origin) {
    Physical z=physical(k,family);alt::sample_index=k;
    V3 ss=primitive(k,origin);
    for(int i=0;i<3;++i) {
        alt::source[i]=z.p[i];alt::source[3+i]=z.v[i];alt::source[6+i]=z.a[i];
        alt::source[9+i]=ss[i];alt::source[12+i]=z.gyro[i];alt::source[15+i]=z.acc[i];
        alt::source[18+i]=z.beta[i];alt::source[21+i]=z.mag[i];alt::source[24+i]=z.angles[i];
    }
    alt::front(f);
    f.update(DT,z.gyro,z.acc,35.0f);
    if(k%8==0)f.updateMag(z.mag);
    alt::front(f);alt::point(40,f.raw().mekf());
}
void begin(Fusion& f) { Fusion::Config cfg;f.begin(cfg); }

struct Boundary {unsigned live=0,active=0;};
Boundary locate(int family) {
    Fusion f;begin(f);Boundary b;
    for(unsigned k=1;k<=120000;++k) {
        advance(f,k,family,b.live);
        if(!b.live&&f.isLive())b.live=k;
        if(b.live&&f.raw().mekf().acc_bias_updates_enabled()) {b.active=k;break;}
    }
    if(!b.live||!b.active||b.active<=b.live+900)
        throw std::runtime_error("analytic probe did not expose a sufficiently long actual H18 and A21 prelude");
    return b;
}
void perturb(Fusion& f,const std::string& name) {
    auto& k=f.raw().mekf();
    if(name=="tilt") {
        k.qref=Eigen::Quaternionf(Eigen::AngleAxisf(1.0f/512.0f,V3::UnitX()))*k.qref;
        k.qref.normalize();
    } else if(name=="held_bias") k.xext.template segment<3>(18)+=V3(1.0f/64.0f,0,0);
    else if(name=="covariance") {
        Eigen::Matrix<float,21,1> v=Eigen::Matrix<float,21,1>::Zero();
        v[0]=1;v[15]=.25f;
        k.Pext.noalias()+=(1.0f/65536.0f)*(v*v.transpose());
    } else if(name!="base") throw std::runtime_error("unknown probe");
}
void word(Fusion& root,unsigned start,unsigned origin,int family,
          const std::filesystem::path& out,const std::string& stem,const std::string& direction) {
    std::cout.flush();std::cerr.flush();
    const pid_t child=fork();
    if(child<0)throw std::runtime_error("fork failed");
    if(child==0) {
        try {
            perturb(root,direction);
            alt::trace.open(out/(stem+"-"+direction+".bin"),std::ios::binary);
            if(!alt::trace)throw std::runtime_error("trace open failed");
            alt::point(40,root.raw().mekf());
            for(unsigned k=start+1;k<=start+600;++k)advance(root,k,family,origin);
            if(!alt::trace)throw std::runtime_error("trace write failed");
            alt::trace.close();
            _exit(0);
        }catch(const std::exception& e){std::cerr<<e.what()<<'\n';_exit(2);}
    }
    int status=0;
    if(waitpid(child,&status,0)!=child||!WIFEXITED(status)||WEXITSTATUS(status)!=0)
        throw std::runtime_error("word child failed");
}
int main(int argc,char**argv) {
    try {
        static_assert(sizeof(float)==4 && std::numeric_limits<float>::is_iec559);
        if(argc!=2)throw std::runtime_error("usage: shipping-word OUTPUT_DIRECTORY");
        const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
        std::ofstream manifest(out/"manifest.json");
        manifest<<"{\"runtime\":\"binary32 unchanged default fusion via temporary read-only hooks\","
            "\"word_steps\":600,\"dt\":"<<std::setprecision(17)<<double(DT)<<",\"words\":[";
        bool first=true;
        for(int family=0;family<3;++family) {
            Boundary b=locate(family);
            std::cerr<<"BIAS"<<family<<" Live sample="<<b.live<<" A21 sample="<<b.active<<'\n';
            const std::array<unsigned,3> starts={b.live+200,b.active-300,b.active+600};
            const std::array<std::string,3> names={"H18","H18_to_A21","A21"};
            Fusion root;begin(root);unsigned cursor=0;
            for(unsigned i=0;i<3;++i) {
                while(cursor<starts[i]) {++cursor;advance(root,cursor,family,b.live);}
                const std::string stem="BIAS"+std::to_string(family)+"-"+names[i];
                for(const std::string direction:{"base","tilt","held_bias","covariance"})
                    word(root,cursor,b.live,family,out,stem,direction);
                if(!first)manifest<<',';first=false;
                manifest<<"{\"name\":\""<<stem<<"\",\"start_sample\":"<<cursor
                    <<",\"live_origin_sample\":"<<b.live<<",\"release_sample\":"<<b.active<<'}';
                std::cerr<<stem<<" complete\n";
            }
        }
        manifest<<"]}\n";
        if(!manifest)throw std::runtime_error("manifest write failed");
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
