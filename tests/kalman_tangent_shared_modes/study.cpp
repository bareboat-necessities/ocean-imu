#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
#include "Study.h"
#include "../kalman_carma/CarmaBenchTruth.h"
#include "../kalman_carma/carma3/Bench3.h"
#include "util/W3dSimCommon.h"
#include "lie/SO3Jacobians.h"
using namespace tangent_study;
namespace cb=carma_bench;
struct Moments {
    long n=0,nm=0,na=0,nw=0,nmp=0,nd=0,ng=0;D ba2=D::Zero(),bg2=D::Zero(),p2=D::Zero(),v2=D::Zero();
    double att2=0,attmax=0,an=0,mn=0,nees=0,bnees=0,hnees=0,mineig=1,asym=0,gate=0;
    D prev_white=D::Zero(),white_sum=D::Zero(),white2=D::Zero(),white_lag=D::Zero();
    double drift_num=0,drift_den=0,tmean=0,pzmean=0,first_t=0;D first_p=D::Zero(),last_p=D::Zero();
    double live=-1,att_convergence=-1,bias_convergence=-1,att_good=-1,bias_good=-1;
};
int main(int argc,char** argv) try {
    std::string scenario="wave",prefix="run",names="QSM,TFG,TAN,QSM-G,TFG-G,TAN-G,CARMA";
    int seed=4101;double duration=900,dt=.005;bool matched=false,population=false;
    for(int k=1;k<argc;++k) {
        const std::string a=argv[k];
        if(a=="--matched-noise"){matched=true;continue;}
        if(a=="--population-prior"){population=true;continue;}
        if(++k==argc)throw std::runtime_error("missing CLI value");
        const std::string b=argv[k];
        if(a=="--scenario")scenario=b;
        else if(a=="--output")prefix=b;
        else if(a=="--filters")names=b;
        else if(a=="--seed")seed=std::stoi(b);
        else if(a=="--duration")duration=std::stod(b);
        else if(a=="--dt")dt=std::stod(b);else throw std::runtime_error("unknown option");
    }
    if(!(dt>=.002&&dt<=.02&&duration>0&&std::isfinite(duration)))throw std::runtime_error("invalid sampling/duration");
    const std::vector<std::string> cases={"stationary","rocking","wave","low","high","transition","yaw","inclination0","inclination87","reference_error","exact_reference","initial_bias","initial_attitude","late_bias","late_attitude","gyro_error","bias_drift","gyro_drift","vibration","colored","long","rate100","rate400","mag100","mag10","jitter"};
    if(std::find(cases.begin(),cases.end(),scenario)==cases.end())throw std::runtime_error("unknown scenario");
    if(scenario=="rate100")dt=.01;
    if(scenario=="rate400")dt=.0025;
    cb::Scenario sc;
    if(!cb::base_scenario("jonswap_fp0.29_hs0.8_g3.3",sc))throw std::runtime_error("missing truth");
    sc.hmodel.physical=true;sc.duration_s=duration;sc.att.yaw_wander_deg=0;sc.att.yaw_osc_deg=0;
    if(scenario=="stationary"||scenario=="rocking")sc.stationary=true;
    if(scenario=="low")sc.parts[0].hs_m=.15;
    if(scenario=="high")sc.parts[0].hs_m=4.;
    if(scenario=="yaw"){sc.att.yaw_osc_deg=45;sc.att.yaw_osc_hz=.012;}
    if(scenario=="bias_drift") {
        sc.bias.ramp_mps2_per_s=D(8e-6,-6e-6,4e-6);
        sc.bias.sin_amp_mps2=D(.02,-.015,.01);sc.bias.sin_period_s=900;
    }
    if(scenario=="gyro_drift") {
        sc.bias.gyro_ramp_rps_per_s=D(1e-7,-8e-8,5e-8);
        sc.bias.gyro_sin_amp_rps=D(.0008,-.0006,.0004);sc.bias.gyro_sin_period_s=600;
    }
    if(scenario=="transition"){sc.env.steps={{300,30,4},{600,30,1}};}
    double inclination=65*cb::PI/180;
    if(scenario=="inclination0")inclination=0;
    if(scenario=="inclination87")inclination=87*cb::PI/180;
    Config config;config.dt=dt;config.matched_noise=matched;config.population_prior=population;config.field=V(float(50*std::cos(inclination)),0,float(50*std::sin(inclination)));
    if(scenario=="reference_error")config.reference="wrong";
    if(scenario=="exact_reference")config.reference="exact";
    if(scenario=="initial_bias"||scenario=="late_bias"){config.shock="bias";config.shock_time=scenario=="initial_bias"?0:300;}
    if(scenario=="initial_attitude"||scenario=="late_attitude"){config.shock="attitude";config.shock_time=scenario=="initial_attitude"?0:300;}
    std::vector<std::string> labels;std::vector<std::unique_ptr<Filter>> filters;
    for(std::size_t p=0;p<=names.size();) {
        const auto q=names.find(',',p);const std::string name=names.substr(p,q==std::string::npos?q:q-p);
        Config c=config;c.gate=name.ends_with("-G");
        if(name=="QSM-LG"||name=="TFG-LG"||name=="TAN-LG"){c.gate=true;c.gate_start=300;}
        std::unique_ptr<Filter> f;
        if(name.starts_with("QSM"))f=make_q(c);
        else if(name.starts_with("TFG"))f=make_tfg(c);
#ifndef TANGENT_BASELINE
        else if(name.starts_with("TAN"))f=make_tangent(c);
#endif
        else if(name=="CARMA")f=make_independent(c);
        if(!f)throw std::runtime_error("unknown filter: "+name);
        labels.push_back(name);filters.push_back(std::move(f));
        if(q==std::string::npos)break;
        p=q+1;
    }
    cb::WaveTruth truth(sc,w3d_expand_seed(seed,21u),dt);
    auto an=make_imu_noise_model(carma3_bench::ACC_SIGMA,8e-3f*9.80665f,0.f,w3d_expand_seed(seed,11u),w3d_expand_seed(seed,12u));
    auto gn=make_imu_noise_model(carma3_bench::GYR_SIGMA,.05f*float(cb::PI/180),0.f,w3d_expand_seed(seed,13u),w3d_expand_seed(seed,14u));
    auto mn=make_mag_noise_model(carma3_bench::MAG_SIGMA_UT,0.f,0.f,0.f,0.f,0.f,w3d_expand_seed(seed,15u),w3d_expand_seed(seed,16u));
    std::ofstream trace(prefix+"_trace.csv"),scores(prefix+"_scores.csv");trace.exceptions(std::ios::failbit|std::ios::badbit);scores.exceptions(std::ios::failbit|std::ios::badbit);
    trace<<std::setprecision(10)<<"t,filter,att,px,py,pz,vx,vy,vz,bax,bay,baz,bgx,bgy,bgz,true_px,true_py,true_pz,true_vx,true_vy,true_vz,true_bax,true_bay,true_baz,true_bgx,true_bgy,true_bgz,nees15,bias_nees,acc_nis,mag_nis,gate_trace,info_min,live,bias_enabled,r00,r01,r02,r10,r11,r12,r20,r21,r22,phi_x,phi_y,phi_z,white_x,white_y,white_z,ax,ay,az\n";
    std::vector<std::pair<double,double>> windows={{std::min(150.,std::max(0.,duration-30.)),duration}};
    if(duration>=900) {windows.push_back({150,300});windows.push_back({300,600});windows.push_back({600,900});}
    if(duration>=3600)windows.push_back({duration*.75,duration});
    const auto nf=filters.size();std::vector<std::vector<Moments>> sums(nf,std::vector<Moments>(windows.size()));
    std::vector<double> runtime(nf,0);std::vector<long> timed_samples(nf,0);std::vector<int> failures(nf,0);std::vector<Output> outputs(nf);
    std::uint64_t hash=1469598103934665603ull;
    auto hashv=[&](const V& v){for(int a=0;a<3;++a){float x=v(a);unsigned char bits[4];std::memcpy(bits,&x,4);for(auto b:bits){hash^=b;hash*=1099511628211ull;}}};
    double next_mag=0,next_score=0,next_diag=0;D color=D::Zero();std::mt19937 rng(w3d_expand_seed(seed,31u));std::normal_distribution<double> normal;
    bool model_written=false;
    long samples=0;double t=0,previous_h=dt;
    while(t<duration-1e-8) {
        const double h=std::min(duration-t,dt*(scenario=="jitter"?1+.15*std::sin(.031*samples):1));
        auto tr=scenario=="jitter"?truth.at(t):truth.next();
        if(std::abs(tr.t-t)>1e-7)throw std::runtime_error("truth timestamp mismatch");
        if(scenario=="rocking") {
            const double r=.14*std::sin(cb::TWO_PI*.13*t),p=.09*std::sin(cb::TWO_PI*.19*t),rd=.14*cb::TWO_PI*.13*std::cos(cb::TWO_PI*.13*t),pd=.09*cb::TWO_PI*.19*std::cos(cb::TWO_PI*.19*t);
            tr.R_nb=(Eigen::AngleAxisd(p,D::UnitY())*Eigen::AngleAxisd(r,D::UnitX())).toRotationMatrix();tr.omega_b=D(rd,pd*std::cos(r),-pd*std::sin(r));tr.f_b=tr.R_nb.transpose()*D(0,0,-cb::G);
        }
        Sample s;s.t=t;s.dt=float(samples?previous_h:h);s.acc=apply_imu_noise(tr.f_b.cast<float>(),an,float(h));s.gyr=apply_imu_noise(tr.omega_b.cast<float>(),gn,float(h));
        D ba=(an.turn_on_bias+an.bias_rw).cast<double>(),bg=(gn.turn_on_bias+gn.bias_rw).cast<double>();
        if(scenario=="bias_drift"){s.acc+=tr.extra_acc_bias_b.cast<float>();ba+=tr.extra_acc_bias_b;}
        if(scenario=="gyro_drift"){s.gyr+=tr.extra_gyro_bias_b.cast<float>();bg+=tr.extra_gyro_bias_b;}
        if(scenario=="gyro_error"){const D b(.002,-.0015,.001);s.gyr+=b.cast<float>();bg+=b;}
        if(scenario=="vibration"){s.acc+=V(float(.08*std::sin(cb::TWO_PI*37*t)),0,float(.15*std::sin(cb::TWO_PI*43*t)));}
        if(scenario=="colored") {double a=std::exp(-h/.3);for(int j=0;j<3;++j)color(j)=a*color(j)+.025*std::sqrt(1-a*a)*normal(rng);s.acc+=color.cast<float>();}
        s.mag_tick=t+1e-8>=next_mag;
        if(s.mag_tick){const double mh=scenario=="mag100"?.01:scenario=="mag10"?.1:.04;next_mag+=mh;s.mag=apply_mag_noise(tr.R_nb.transpose().cast<float>()*config.field,mn,float(mh));}
        hashv(s.acc);hashv(s.gyr);
        if(s.mag_tick)hashv(s.mag);
        // Bind timestamps/elapsed intervals, without hashing struct padding.
        unsigned char time_bits[sizeof(double)],dt_bits[sizeof(float)];
        std::memcpy(time_bits,&s.t,sizeof(double));std::memcpy(dt_bits,&s.dt,sizeof(float));
        for(auto b:time_bits){hash^=b;hash*=1099511628211ull;}
        for(auto b:dt_bits){hash^=b;hash*=1099511628211ull;}
        const bool score=t+1e-8>=next_score,diag=t+1e-8>=next_diag;
        if(score)next_score+=.05;
        if(diag)next_diag+=1.;
        for(std::size_t f=0;f<nf;++f) {
            if(failures[f])continue;
            const auto start=std::clock();filters[f]->step(s);const double elapsed=double(std::clock()-start)/CLOCKS_PER_SEC;
            if(t>=std::min(150.,duration*.5)){runtime[f]+=elapsed;++timed_samples[f];}
            if(!score&&!diag&&!s.mag_tick)continue;
            const Output o=filters[f]->output(diag);
            if(!o.finite){++failures[f];continue;}
            const D phi=ocean_imu::lie::Log<double>(Eigen::Matrix3d(o.R*tr.R_nb.transpose()));
            const double att=phi.norm()*180/cb::PI;const D pe=o.p-tr.p_wave,ve=o.v-tr.v_wave,be=o.ba-ba,ge=o.bg-bg;
            double nees=NAN,bnees=NAN,hnees=NAN;
            if(diag&&o.live) {
                Eigen::Matrix<double,15,1> e;e<<phi,pe,ve,ge,be;
                Eigen::LDLT<P15> ldlt((o.P+o.P.transpose())*.5);
                if(ldlt.isPositive())nees=e.dot(ldlt.solve(e));
                const Eigen::Matrix3d B=o.P.bottomRightCorner<3,3>();Eigen::LDLT<Eigen::Matrix3d> bl(B);
                if(bl.isPositive())bnees=be.dot(bl.solve(be));
                const Eigen::Matrix2d Bh=B.topLeftCorner<2,2>();Eigen::LDLT<Eigen::Matrix2d> hl(Bh);
        if(hl.isPositive())hnees=be.head<2>().dot(hl.solve(be.head<2>()));
            }
            for(std::size_t w=0;w<windows.size();++w) {
                auto& m=sums[f][w];
                if(o.live&&m.live<0)m.live=t;
                if(o.live&&t>=std::max(0.,config.shock_time)&&att<1){if(m.att_good<0)m.att_good=t;
        if(t-m.att_good>=30&&m.att_convergence<0)m.att_convergence=m.att_good;}else m.att_good=-1;
                if(o.live&&t>=std::max(0.,config.shock_time)&&be.norm()<.01){if(m.bias_good<0)m.bias_good=t;
        if(t-m.bias_good>=30&&m.bias_convergence<0)m.bias_convergence=m.bias_good;}else m.bias_good=-1;
                if(t<windows[w].first||t>=windows[w].second)continue;
                if(score) {
                    if(m.n==0){m.first_p=pe;m.first_t=t;}m.last_p=pe;++m.n;
                    m.ba2+=be.cwiseProduct(be);m.bg2+=ge.cwiseProduct(ge);m.p2+=pe.cwiseProduct(pe);m.v2+=ve.cwiseProduct(ve);m.att2+=att*att;m.attmax=std::max(m.attmax,att);
                    m.gate+=o.gate_trace;m.ng+=o.bias_enabled;
                    const double u=t-m.first_t,du=u-m.tmean,dy=pe.z()-m.pzmean;m.tmean+=du/m.n;m.pzmean+=dy/m.n;m.drift_num+=du*(pe.z()-m.pzmean);m.drift_den+=du*(u-m.tmean);
                    if(o.acc_applied&&std::isfinite(o.anis)){m.an+=o.anis;++m.na;Eigen::LLT<Eigen::Matrix3d> ll(o.aS);
        if(ll.info()==Eigen::Success){D x=ll.matrixL().solve(o.ar);++m.nw;m.white_sum+=x;m.white2+=x.cwiseProduct(x);
        if(m.nw>1) {m.white_lag+=x.cwiseProduct(m.prev_white);}
                        m.prev_white=x;}}
                }
                if(s.mag_tick)++m.nmp;
                if(s.mag_tick&&o.mag_applied&&std::isfinite(o.mnis)){m.mn+=o.mnis;++m.nm;}
                if(diag&&std::isfinite(nees)&&std::isfinite(bnees)){m.nees+=nees;m.bnees+=bnees;m.hnees+=hnees;++m.nd;m.mineig=std::min(m.mineig,o.min_eig);m.asym=std::max(m.asym,o.asym);}
            }
            if(diag && !model_written && t>=std::min(600.,duration*.7) && f==0) {
                std::ofstream model(prefix+"_model.json");model.exceptions(std::ios::failbit|std::ios::badbit);filters[f]->write_model(model);model_written=true;
            }
            if(diag) {
                trace<<t<<','<<labels[f]<<','<<att;
                for(const auto& v:{o.p,o.v,o.ba,o.bg,tr.p_wave,tr.v_wave,ba,bg})for(int j=0;j<3;++j)trace<<','<<v(j);
                trace<<','<<nees<<','<<bnees<<','<<o.anis<<','<<o.mnis<<','<<o.gate_trace<<','<<o.gate_eigen.minCoeff()<<','<<o.live<<','<<o.bias_enabled;
                for(int i=0;i<3;++i)for(int j=0;j<3;++j)trace<<','<<o.R(i,j);
                for(int i=0;i<3;++i)trace<<','<<phi(i)*180/cb::PI;
                D white=D::Constant(NAN);Eigen::LLT<Eigen::Matrix3d> ll(o.aS);
        if(ll.info()==Eigen::Success)white=ll.matrixL().solve(o.ar);
                for(int i=0;i<3;++i)trace<<','<<white(i);
                for(int i=0;i<3;++i)trace<<','<<o.a(i);
                trace<<'\n';
            }
        }
        ++samples;t=scenario=="jitter"?t+h:double(samples)*dt;previous_h=h;
    }
    scores<<std::setprecision(12)<<"scenario,seed,filter,start,end,n,att_rms,att_max,p_x,p_y,heave_rms,v_x,v_y,v_z,ba_x,ba_y,ba_z,bg_x,bg_y,bg_z,horizontal_drift,heave_slope,nees15,bias_nees,h_bias_nees,acc_nis,mag_nis,acc_applied_fraction,mag_applied_fraction,innovation_acf1_x,innovation_acf1_y,innovation_acf1_z,gate_trace,bias_enabled_fraction,live_time,att_convergence,bias_convergence,min_correlation_eigenvalue,asymmetry,failed,step_us,memory_bytes,states,samples,sensor_hash\n";
    for(std::size_t f=0;f<nf;++f)for(std::size_t w=0;w<windows.size();++w) {
        const auto& m=sums[f][w];const double n=std::max(1L,m.n),nd=std::max(1L,m.nd);
        scores<<scenario<<','<<seed<<','<<labels[f]<<','<<windows[w].first<<','<<windows[w].second<<','<<m.n<<','<<std::sqrt(m.att2/n)<<','<<m.attmax;
        for(const auto& a:{m.p2,m.v2,m.ba2,m.bg2})for(int j=0;j<3;++j)scores<<','<<std::sqrt(a(j)/n);
        const double nw=double(std::max(1L,m.nw));
        const D mu=m.white_sum/nw,var=m.white2/nw-mu.cwiseProduct(mu),acf=(m.white_lag/std::max(1.,nw-1)-mu.cwiseProduct(mu)).cwiseQuotient(var.cwiseMax(1e-20));
        scores<<','<<(m.last_p-m.first_p).head<2>().norm()<<','<<(m.drift_den>0?m.drift_num/m.drift_den:NAN)<<','<<m.nees/nd<<','<<m.bnees/nd<<','<<m.hnees/nd<<','<<m.an/std::max(1L,m.na)<<','<<m.mn/std::max(1L,m.nm)<<','<<double(m.na)/n<<','<<double(m.nm)/std::max(1L,m.nmp);
        for(int j=0;j<3;++j)scores<<','<<acf(j);
        scores<<','<<m.gate/n<<','<<double(m.ng)/n<<','<<m.live<<','<<m.att_convergence<<','<<m.bias_convergence<<','<<m.mineig<<','<<m.asym<<','<<failures[f]<<','<<runtime[f]*1e6/std::max(1L,timed_samples[f])<<','<<filters[f]->memory()<<','<<filters[f]->states()<<','<<samples<<','<<std::hex<<hash<<std::dec<<'\n';
    }
    std::cout<<scenario<<" seed="<<seed<<" samples="<<samples<<" sensor_hash="<<std::hex<<hash<<std::dec<<" failures=";
    int failed=0;for(auto f:failures)failed+=f;std::cout<<failed<<'\n';return failed?2:0;
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
