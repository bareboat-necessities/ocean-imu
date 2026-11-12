#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "StandardHeaders.h"
#include "util/W3dSimCommon.h"
#include "Experiment.h"
#include "../kalman_carma/CarmaBenchTruth.h"
#include "../kalman_carma/carma3/Bench3.h"
using namespace followup;
using namespace carma_bench;
using namespace carma3_bench;
using VD=Eigen::Vector3d;
using MD=Eigen::Matrix3d;
struct Window {double start,end; long n=0,nm=0; VD ba2=VD::Zero();double tilt2=0,att2=0,aw2=0,an=0,mn=0,on=0;long na=0,no=0,en=0,applied=0;};
int main(int argc,char** argv) try {
    std::string geometry="qmekf",motion="wave",prefix="run"; Options opt;opt.arm="exact_fixed";int seed=1001;double duration=3600,maglag=0;bool test=false,config_parity=false;
    for(int i=1;i<argc;++i) {std::string k=argv[i];if(k=="--config-parity"){config_parity=true;continue;}if(k=="--parity"){opt.parity=true;continue;}if(k=="--self-test"){test=true;continue;}
        if(i+1>=argc)throw std::runtime_error("missing argument");
        std::string v=argv[++i];
        if(k=="--mag-lag")maglag=std::stod(v);else if(k=="--noise")opt.noise=v;else if(k=="--timing")opt.timing=v;else if(k=="--freeze")opt.freeze=v;else if(k=="--injection")opt.injection=v;else if(k=="--change")opt.change=std::stod(v);else if(k=="--window")opt.window=std::stoi(v)!=0;else if(k=="--filter")geometry=v;else if(k=="--motion")motion=v;else if(k=="--arm")opt.arm=v;else if(k=="--seed")seed=std::stoi(v);else if(k=="--duration")duration=std::stod(v);else if(k=="--sigma")opt.oracle_sigma=std::stof(v);else if(k=="--output")prefix=v;else throw std::runtime_error("unknown argument");}
    if(!std::set<std::string>{"native","gyro","mag","all","prior","late_prior","all_prior","schmidt","schmidt_prior"}.count(opt.noise))throw std::runtime_error("invalid noise model");
    if(!std::set<std::string>{"endpoint","core","frontend","both","both_start","interval"}.count(opt.timing)||!std::set<std::string>{"none","sea","coupling","both"}.count(opt.freeze)||!std::set<std::string>{"none","bias","tilt"}.count(opt.injection))throw std::runtime_error("invalid intervention");
    if(!std::isfinite(maglag)||maglag<0||maglag>.1)throw std::runtime_error("invalid magnetic latency");
    const std::set<std::string> arms={"learned_cont","learned_fixed","exact_cont","exact_fixed","oracle","oracle_bias","correct_reference"};
    if(!arms.count(opt.arm)||(motion!="wave"&&motion!="rocking")||duration<=0||opt.oracle_sigma<=0)throw std::runtime_error("invalid configuration");
    MD Rmg=Eigen::AngleAxisd(MagSim_WMM::default_declination_deg*PI/180,VD::UnitZ()).toRotationMatrix();
    V Bgeo=MagSim_WMM::mag_world_aero(), B=(Rmg.transpose()*Bgeo.cast<double>()).cast<float>();
    auto f=geometry=="qmekf"?make_qmekf(opt,B):geometry=="tfg"?make_tfg(opt,B):nullptr;
    if(!f)throw std::runtime_error("unknown filter");
    std::unique_ptr<Experiment> configuration_control;
    if(config_parity){
        if(opt.timing!="both_start"||opt.change!=0)throw std::runtime_error("configuration parity requires both_start at zero");
        Options control=opt;control.timing="both";
        configuration_control=geometry=="qmekf"?make_qmekf(control,B):make_tfg(control,B);
    }
    if(test){f->self_test();std::cout<<"oracle Jacobian, correction direction and covariance PASS\n";return 0;}
    Scenario sc;if(!base_scenario("jonswap_fp0.29_hs0.8_g3.3",sc))throw std::runtime_error("scenario missing");sc.hmodel.physical=true;sc.duration_s=duration;sc.att.yaw_wander_deg=0;sc.att.yaw_osc_deg=0;
    constexpr double dt=.005;WaveTruth truth(sc,w3d_expand_seed(seed,21u),dt);
    auto an=make_imu_noise_model(ACC_SIGMA,8e-3f*9.80665f,0.f,w3d_expand_seed(seed,11u),w3d_expand_seed(seed,12u));
    auto gn=make_imu_noise_model(GYR_SIGMA,.05f*float(PI/180),0.f,w3d_expand_seed(seed,13u),w3d_expand_seed(seed,14u));
    auto mn=make_mag_noise_model(MAG_SIGMA_UT,0.f,0.f,0.f,0.f,0.f,w3d_expand_seed(seed,15u),w3d_expand_seed(seed,16u));
    float phase=0;std::uint64_t hash=1469598103934665603ull,original_hash=1469598103934665603ull;
    MD previous_R=MD::Identity();
    std::deque<MD> magnetic_history;

    auto hashvec=[&](const V& v,std::uint64_t& digest){for(int j=0;j<3;++j){float x=v[j];unsigned char b[4];std::memcpy(b,&x,4);for(auto c:b){digest^=c;digest*=1099511628211ull;}}};
    std::ofstream history;
    if(opt.window){history.open(prefix+"_window.csv");history<<"t,ax,ay,az,gx,gy,gz,mx,my,mz,mag_tick,qw,qx,qy,qz,bax,bay,baz,bgx,bgy,bgz,true_bax,true_bay,true_baz,true_bgx,true_bgy,true_bgz,true_qw,true_qx,true_qy,true_qz\n";history<<std::setprecision(12);}
    std::ofstream trace(prefix+"_trace.csv"),diag(prefix+"_diag.csv");
    trace<<"t,tilt,att,ba_x,ba_y,ba_z,true_x,true_y,true_z,aw_error,hi_norm,inc_error,enabled,acc_nis,mag_nis,oracle_nis,live,refined,mag_applied,mag_correction_rad,mag_bias_correction\n";
    diag<<"t,ab_corr,nuisance_inflation,min_corr_eig,bias_nees,joint_nees,h_nees,h_sigma,hi_norm,inc_error,enabled\n";
    trace<<std::setprecision(9);diag<<std::setprecision(10);
    std::vector<Window> windows;for(auto p:std::vector<std::pair<double,double>>{{150,600},{600,900},{900,1200},{2400,3600},{5400,7200}})if(p.first<duration)windows.push_back({p.first,std::min(p.second,duration)});
    std::deque<double> rolling;double rolling_sum=0,qualify=-1,recovery=NAN;
    for(long k=0;k<long(std::llround(duration/dt));++k) {
        auto tr=truth.next();double t=k*dt;
        if(motion=="rocking") {
            double r=8*PI/180*std::sin(TWO_PI*.13*t),p=5*PI/180*std::sin(TWO_PI*.19*t),rd=8*PI/180*TWO_PI*.13*std::cos(TWO_PI*.13*t),pd=5*PI/180*TWO_PI*.19*std::cos(TWO_PI*.19*t);
            tr.R_nb=(Eigen::AngleAxisd(p,VD::UnitY())*Eigen::AngleAxisd(r,VD::UnitX())).toRotationMatrix();tr.omega_b=VD(rd,pd*std::cos(r),-pd*std::sin(r));tr.a_ned.setZero();tr.f_b=tr.R_nb.transpose()*VD(0,0,-G);
        }
        Input in;in.t=float(t);in.dt=float(dt);in.acc=apply_imu_noise(tr.f_b.cast<float>(),an,in.dt);in.gyr=apply_imu_noise(tr.omega_b.cast<float>(),gn,in.dt);in.ba=an.turn_on_bias+an.bias_rw;
        in.aw=(Rmg.transpose()*tr.a_ned).cast<float>();in.mag=V::Zero();phase+=in.dt;in.mag_tick=phase>=.04f;
        if(maglag>0){magnetic_history.push_back(tr.R_nb);while(magnetic_history.size()>std::size_t(std::ceil(maglag/dt))+2)magnetic_history.pop_front();}
        if(in.mag_tick){
            while(phase>=.04f)phase-=.04f;
            MD sensor_R=tr.R_nb;
            if(maglag>0){
                const std::size_t age=std::size_t(std::floor(maglag/dt));const double fraction=maglag/dt-double(age);
                const std::size_t newest=magnetic_history.size()-1;
                const std::size_t i0=newest>=age?newest-age:0,i1=i0?i0-1:0;
                sensor_R=Eigen::Quaterniond(magnetic_history[i0]).slerp(fraction,Eigen::Quaterniond(magnetic_history[i1])).toRotationMatrix();
            }
            in.mag=apply_mag_noise(sensor_R.transpose().cast<float>()*Bgeo,mn,.04f);
        }
        const V endpoint_gyr=in.gyr;
        if(opt.timing=="interval"&&t>=opt.change&&k>0){
            Eigen::AngleAxisd increment(previous_R.transpose()*tr.R_nb);
            in.gyr+=(VD(increment.angle()*increment.axis()/dt)-tr.omega_b).cast<float>();
        }
        previous_R=tr.R_nb;
        hashvec(in.acc,hash);hashvec(in.gyr,hash);if(in.mag_tick)hashvec(in.mag,hash);
        hashvec(in.acc,original_hash);hashvec(endpoint_gyr,original_hash);if(in.mag_tick)hashvec(in.mag,original_hash);
        const bool diagnostic=k%200==0;Output o=f->step(in,diagnostic||config_parity);
        if(configuration_control){
            const Output control=configuration_control->step(in,true);
            if((o.q.coeffs()-control.q.coeffs()).norm()!=0||(o.P-control.P).norm()!=0||
               (o.ba-control.ba).norm()!=0||(o.aw-control.aw).norm()!=0)
                throw std::runtime_error("configuration and immediate setter paths differ");
        }
        MD Rt=Rmg.transpose()*tr.R_nb,Re=o.q.cast<double>().toRotationMatrix();
        if(opt.window&&t>=opt.change&&t<=opt.change+60.001){
            if(k==long(std::llround(opt.change/dt))){std::ofstream model(prefix+"_model.json");f->write_model(model);}
            Eigen::Quaterniond tq(Rt);const V truebg=gn.turn_on_bias+gn.bias_rw;
            history<<t;
            for(const auto& v:{in.acc,endpoint_gyr,in.mag})for(int j=0;j<3;++j)history<<','<<v[j];
            history<<','<<in.mag_tick<<','<<o.q.w()<<','<<o.q.x()<<','<<o.q.y()<<','<<o.q.z();
            for(const auto& v:{o.ba,o.bg,in.ba,truebg})for(int j=0;j<3;++j)history<<','<<v[j];
            history<<','<<tq.w()<<','<<tq.x()<<','<<tq.y()<<','<<tq.z()<<'\n';
        }
        double tilt=std::acos(std::clamp((Rt.transpose()*VD::UnitZ()).dot(Re.transpose()*VD::UnitZ()),-1.,1.))*180/PI;
        Eigen::AngleAxisd aa(Rt.transpose()*Re);double att=aa.angle()*180/PI;VD be=(o.ba-in.ba).cast<double>();double aw=(o.aw-in.aw).norm();
        double inc=std::atan2(o.field.z(),o.field.head<2>().norm())*180/PI-std::atan2(B.z(),B.head<2>().norm())*180/PI;
        for(auto& w:windows)if(t>=w.start&&t<w.end){++w.n;w.ba2+=be.cwiseProduct(be);w.tilt2+=tilt*tilt;w.att2+=att*att;w.aw2+=aw*aw;w.en+=o.enabled;if(std::isfinite(o.acc_nis)){w.an+=o.acc_nis;++w.na;}if(in.mag_tick){++w.nm;w.applied+=o.mag_applied;if(std::isfinite(o.mag_nis))w.mn+=o.mag_nis;if(std::isfinite(o.oracle_nis)){++w.no;w.on+=o.oracle_nis;}}}
        if(t>=opt.change&&o.enabled&&o.live){double e=be.head<2>().squaredNorm();rolling.push_back(e);rolling_sum+=e;if(rolling.size()>6000){rolling_sum-=rolling.front();rolling.pop_front();}
            if(rolling.size()==6000 && rolling_sum/6000<.0001){if(qualify<0)qualify=t;if(t-qualify>=30 && !std::isfinite(recovery))recovery=qualify-opt.change;}else qualify=-1;}
        if((k%20==0&&(t<=180||(t>=590&&t<=900)))||diagnostic||(in.mag_tick&&(t<=180||(t>=590&&t<=900))))trace<<t<<','<<tilt<<','<<att<<','<<o.ba.x()<<','<<o.ba.y()<<','<<o.ba.z()<<','<<in.ba.x()<<','<<in.ba.y()<<','<<in.ba.z()<<','<<aw<<','<<o.hi.norm()<<','<<inc<<','<<o.enabled<<','<<o.acc_nis<<','<<o.mag_nis<<','<<o.oracle_nis<<','<<o.live<<','<<o.refined<<','<<o.mag_applied<<','<<o.mag_angle<<','<<o.mag_bias<<'\n';
        if(diagnostic&&t>=150){
            Eigen::MatrixXd P=o.P;Eigen::VectorXd ds=P.diagonal().cwiseMax(1e-30).cwiseSqrt().cwiseInverse();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ds.asDiagonal()*P*ds.asDiagonal());double mineig=es.eigenvalues().minCoeff();
            if(!P.allFinite()||mineig<-.0001||!o.live||!o.refined||!o.enabled)throw std::runtime_error("live covariance/gate invariant failed");
            Eigen::Matrix<double,6,6> T;Eigen::Matrix<double,6,21> N;int target[6]={0,1,2,24,25,26};
            for(int i=0;i<6;++i){for(int j=0;j<6;++j)T(i,j)=P(target[i],target[j]);for(int j=0;j<21;++j)N(i,j)=P(target[i],3+j);}
            Eigen::MatrixXd nn=P.block(3,3,21,21);Eigen::Matrix<double,6,6> C=T-N*nn.ldlt().solve(N.transpose());C=(C+C.transpose()).eval()*.5;
            Eigen::LLT<Eigen::Matrix<double,6,6>> lt(T);Eigen::Matrix<double,6,6> Linv=lt.matrixL().solve(Eigen::Matrix<double,6,6>::Identity());
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> eig(Linv*C*Linv.transpose());double inflation=1/std::max(1e-12,eig.eigenvalues().minCoeff());
            MD A=T.topLeftCorner<3,3>(),Bcov=T.bottomRightCorner<3,3>();MD La=A.llt().matrixL(),Lb=Bcov.llt().matrixL();
            MD W=La.inverse()*T.topRightCorner<3,3>()*Lb.transpose().inverse();double corr=std::sqrt((W*W.transpose()).selfadjointView<Eigen::Lower>().eigenvalues().maxCoeff());
            Eigen::Matrix<double,6,1> err;err.head<3>()=aa.angle()*aa.axis();err.tail<3>()=-be;
            double nees=err.dot(T.ldlt().solve(err)),bn=be.dot(Bcov.ldlt().solve(be));Eigen::Vector2d bh=be.head<2>();double hn=bh.dot(Bcov.topLeftCorner<2,2>().ldlt().solve(bh));
            diag<<t<<','<<corr<<','<<inflation<<','<<mineig<<','<<bn<<','<<nees<<','<<hn<<','<<std::sqrt(Bcov(0,0)+Bcov(1,1))<<','<<o.hi.norm()<<','<<inc<<','<<o.enabled<<'\n';
        }
    }
    std::ofstream scores(prefix+"_scores.csv");scores<<"start,end,n,ba_x_rms,ba_y_rms,ba_z_rms,ba_h_rms,tilt_rms,att_rms,aw_rms,acc_nis,mag_nis,oracle_nis,enabled_fraction,mag_applied_fraction\n";scores<<std::setprecision(10);
    for(const auto& w:windows) {VD b=(w.ba2/w.n).cwiseSqrt();scores<<w.start<<','<<w.end<<','<<w.n<<','<<b.x()<<','<<b.y()<<','<<b.z()<<','<<std::sqrt((w.ba2.x()+w.ba2.y())/w.n)<<','<<std::sqrt(w.tilt2/w.n)<<','<<std::sqrt(w.att2/w.n)<<','<<std::sqrt(w.aw2/w.n)<<','<<w.an/w.na<<','<<w.mn/w.nm<<','<<(w.no?w.on/w.no:NAN)<<','<<double(w.en)/w.n<<','<<double(w.applied)/w.nm<<'\n';}
    std::cout<<"samples="<<long(duration/dt)<<" sensor_hash="<<std::hex<<hash<<std::dec<<" original_sensor_hash="<<std::hex<<original_hash<<std::dec<<" recovery_bias="<<recovery<<" parity="<<opt.parity<<'\n';
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
