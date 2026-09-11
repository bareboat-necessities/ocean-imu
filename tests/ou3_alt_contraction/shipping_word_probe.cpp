// Host-only actual finite-increment probe. Never linked into deployment.
// Built against an audited temporary observation-only header overlay.
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/wait.h>
#include <unistd.h>

#define EIGEN_NON_ARDUINO
#include "kalman_ou_common/KalmanOUCoreMath.h"
template<class F> void alt_state(const char*, const F&);
template<class F, class V> void alt_solve(const char*, const F&, const V&);
template<class F> void alt_solve_finish(const F&);
template<class F> struct AltScope {
    const char* name; const F& filter;
    AltScope(const char* n, const F& f): name(n), filter(f) {}
    ~AltScope() { alt_state(name, filter); }
};
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
extern const float g_std = 9.80665f;
using Wrapper = SeaStateFusion_OU_III<TrackerType::KALMANF>;
using Vec = Eigen::Vector3d;
using Quat = Eigen::Quaterniond;
constexpr double H = 0.005;
constexpr int WORD = 600;
constexpr double W = 2.0 * 3.14159265358979323846 / 8.0;
std::ofstream states, solves;
int sample_index = 0;
bool detailed = false;
Wrapper* current_wrapper = nullptr;

struct Source {
    int family = 2;
    int root_axis = -1;
    double epsilon = 0.001;
    double amplitude = 1.2;
    Vec beta(double t) const {
        double b = family == 0 ? .01*std::exp(-t/1200.0)+.015+.004*std::sin(W*t/75.0)
                 : family == 1 ? .04*std::exp(-t/1200.0)+.008*std::sin(W*t/75.0)
                 : .02+.005*std::sin(W*t/75.0);
        Vec out = Vec(1.0,-.5,.25)*b;
        if (root_axis >= 0)
            out(root_axis) += epsilon*.1*(family == 1 ? std::exp(-t/1200.0) : 1.0);
        return out;
    }
    Vec primitive(double t) const {
        return amplitude/W*Vec(-.25*std::cos(W*t), .15*std::sin(W*t), -std::cos(W*t+.3));
    }
    void sample(int k, Eigen::Vector3f& gyro, Eigen::Vector3f& acc,
                Eigen::Vector3f& mag) const {
        const double t = k*H;
        const double roll=.08*std::sin(W*t), pitch=.05*std::sin(W*t+.7);
        const double dr=.08*W*std::cos(W*t), dp=.05*W*std::cos(W*t+.7);
        const Quat bw = Eigen::AngleAxisd(.3,Vec::UnitZ()) *
            Eigen::AngleAxisd(pitch,Vec::UnitY()) * Eigen::AngleAxisd(roll,Vec::UnitX());
        const Vec aw = -amplitude*W*W*Vec(.25*std::sin(W*t), .15*std::cos(W*t), std::sin(W*t+.3));
        gyro = Vec(dr,dp*std::cos(roll),-dp*std::sin(roll)).cast<float>();
        acc = (bw.conjugate()*(aw-Vec(0,0,9.80665))+beta(t)).cast<float>();
        mag = (bw.conjugate()*Vec(20.0,0.0,45.0)).cast<float>();
    }
};

template<class D> void values(std::ostream& out, const Eigen::MatrixBase<D>& x) {
    for (int i=0;i<x.rows();++i) for(int j=0;j<x.cols();++j)
        out << ',' << static_cast<double>(x(i,j));
}
template<class F> void alt_state(const char* tag, const F& f) {
    if (!states.is_open()) return;
    states << tag << ',' << sample_index << ',' << f.acc_bias_updates_enabled_;
    states << ',' << double(f.qref.w()) << ',' << double(f.qref.x())
           << ',' << double(f.qref.y()) << ',' << double(f.qref.z());
    values(states,f.xext);
    states << ',' << double(f.tau_aw) << ',' << double(f.pseudo_update_period_s_);
    values(states,f.R_S.diagonal());
    const auto& front = current_wrapper->raw();
    states << ',' << front.getFreqHz() << ',' << front.getWavePeriodSec()
           << ',' << front.getSigmaApplied() << ',' << front.getTauApplied()
           << ',' << front.getRSApplied() << '\n';
}
template<class F, class V> void alt_solve(const char* tag, const F& f, const V& r) {
    if (!detailed || !solves.is_open()) return;
    solves << tag << ',' << sample_index << ',' << f.acc_bias_updates_enabled_;
    values(solves,r); values(solves,f.PCt_scratch_); values(solves,f.S_scratch_);
    values(solves,f.K_scratch_); values(solves,f.Pext);
    values(solves,f.xext);
}
template<class F> void alt_solve_finish(const F& f) {
    if (!detailed || !solves.is_open()) return;
    values(solves,f.xext); solves << '\n';
}
void begin(Wrapper& f) {
    Wrapper::Config cfg;
    cfg.with_mag = true;
    cfg.sigma_a = Eigen::Vector3f::Constant(.0148f*.71f);
    cfg.sigma_g = Eigen::Vector3f::Constant(.00157f*.05f);
    cfg.sigma_m = Eigen::Vector3f::Constant(.25f*2.0f);
    f.begin(cfg);
    f.raw().setPeriodicAwCovarianceSync(true);
    f.raw().setAwCovarianceSyncCongruent(false);
    f.raw().enableTuner(true); f.raw().enableClamp(true);
}
void step(Wrapper& f, const Source& source, int k) {
    sample_index=k;
    Eigen::Vector3f gyro,acc,mag;
    source.sample(k,gyro,acc,mag);
    f.update(float(H),gyro,acc,35.0f);
    alt_state("wrapper_imu",f.raw().mekf());
    if(k%8==0) {
        f.updateMag(mag);
        alt_state("wrapper_mag",f.raw().mekf());
    }
}
void perturb(Wrapper& f, int direction, double epsilon) {
    if(direction<0) return;
    static constexpr double scale[7]={.1,.01,1,1,10,1,.1};
    auto& m=f.raw().mekf();
    if(direction<3) {
        Eigen::Vector3f d=Eigen::Vector3f::Zero();
        d(direction)=float(-epsilon*scale[0]);
        m.qref=ocean_imu::kalman::ou_detail::quat_from_delta_theta(d)*m.qref;
        m.qref.normalize();
    } else m.xext(direction)-=float(epsilon*scale[direction/3]);
}
void word(Wrapper& f, const Source& source, int start, int live,
          int direction, const std::string& path) {
    perturb(f,direction,source.epsilon);
    states.open(path+".csv");
    if(!states) throw std::runtime_error("cannot open state tape");
    states << std::setprecision(17);
    if(detailed) { solves.open(path+".solves.csv"); solves << std::setprecision(17); }
    current_wrapper=&f; sample_index=start;
    alt_state("word_start",f.raw().mekf());
    for(int k=start+1;k<=start+WORD;++k) step(f,source,k);
    alt_state("word_end",f.raw().mekf());
    states.close(); solves.close();
    std::ofstream meta(path+".meta");
    meta << std::setprecision(17) << source.family << ' ' << source.root_axis << ' '
         << source.epsilon << ' ' << source.amplitude << ' ' << start << ' '
         << live << ' ' << direction << '\n';
}
int main(int argc,char** argv) try {
    if(argc!=9) throw std::runtime_error("args: pilot|probe family root_axis epsilon amplitude epochs output detailed_direction(-2=off,-3=all,-5=baseline+0/8/13)");
    const std::string mode=argv[1];
    Source source{std::stoi(argv[2]),std::stoi(argv[3]),std::stod(argv[4]),std::stod(argv[5])};
    const int detail_direction=std::stoi(argv[8]);
    if(source.family<0 || source.family>2 || source.root_axis< -1 || source.root_axis>2 ||
       !(std::isfinite(source.epsilon)&&source.epsilon!=0&&std::abs(source.epsilon)<=.01) ||
       !(source.amplitude>0 && source.amplitude<=2.0)) throw std::runtime_error("invalid probe parameters");
    Wrapper f; begin(f); current_wrapper=&f;
    int live=-1,active=-1;
    if(mode=="pilot") {
        for(int k=1;k<=48000;++k) {
            step(f,source,k);
            if(live<0&&f.isLive()) live=k;
            if(active<0&&f.isLive()&&f.raw().mekf().acc_bias_updates_enabled()) { active=k; break; }
        }
        if(live<0||active<0||active-live<800) throw std::runtime_error("no separated natural H18/A21 checkpoints");
        std::cout << live << ' ' << active << '\n'; return 0;
    }
    if(mode!="probe") throw std::runtime_error("unknown mode");
    std::ifstream epochs(argv[6]);
    std::vector<std::pair<int,std::string>> checks; int k0; std::string name;
    while(epochs>>name>>k0) checks.push_back({k0,name});
    if(checks.size()!=3) throw std::runtime_error("three explicit checkpoints required");
    std::sort(checks.begin(),checks.end());
    int k=0;
    for(const auto& [start,label]:checks) {
        while(k<start) {
            step(f,source,++k);
            if(live<0&&f.isLive()) live=k;
        }
        if(!f.isLive()) throw std::runtime_error("checkpoint not Live for this source root");
        for(int d=-1;d<21;++d) {
            if(source.root_axis>=0&&d>=0) continue;
            // fork preserves the ENTIRE covariance/frontend/guard/source history.
            // No hand-written clone or second Riccati recursion is substituted.
            pid_t child=fork();
            if(child<0) throw std::runtime_error("fork failed");
            if(child==0) {
                detailed=(detail_direction==-3 || d==detail_direction ||
                    (detail_direction>=0 && d==-1) ||
                    (detail_direction==-5 && (d==-1 || d==0 || d==8 || d==13)));
                std::string prefix=std::string(argv[7])+"/"+label+"_"+
                    (source.root_axis<0 ? std::to_string(d) : "beta"+std::to_string(source.root_axis));
                word(f,source,start,live,d,prefix);
                _exit(0);
            }
            int status=0;
            if(waitpid(child,&status,0)<0 || !WIFEXITED(status) || WEXITSTATUS(status)!=0)
                throw std::runtime_error("word child failed");
        }
    }
    return 0;
} catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
