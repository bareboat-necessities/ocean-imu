#pragma once
// Test-only access: no deployed header or API is edited. Preload the standard
// library and Eigen so the access macro affects only project test inclusion.
#include "StandardHeaders.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#define private public
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#undef private
#include "Experiment.h"
#include "../kalman_carma/carma3/Bench3.h"
namespace followup {
inline Eigen::Matrix3f skew(const V& v) {Eigen::Matrix3f s;s<<0,-v.z(),v.y(),v.z(),0,-v.x(),-v.y(),v.x(),0;return s;}
template<class Native,bool TFG> struct Probe:Native {
    using Native::Native;
    using M=typename Native::MatrixNX;
    using H3=Eigen::Matrix<float,3,27>;
    bool exact=false;
    V exact_field=V::UnitX(),learned_field=V::UnitX();
    double oracle_nis=NAN,mag_angle=0,mag_bias=0;
    bool mag_called=false;
    void measurement_update_mag_only(const V& mag) {
        mag_called=true;const auto q=this->quaternion_boat();const V ba=this->get_acc_bias();
        Native::measurement_update_mag_only(mag);
        mag_angle=Eigen::AngleAxisf(q.conjugate()*this->quaternion_boat()).angle();
        mag_bias=(ba-this->get_acc_bias()).norm();
    }
    bool oracle_tick=false;
    V oracle_truth=V::Zero();
    float oracle_sigma=.01f;
    void measurement_update_acc_only(const V& acc,float temp=Native::tempC_ref) {
        if(oracle_tick)observe_acceleration(oracle_truth,oracle_sigma);
        Native::measurement_update_acc_only(acc,temp);
    }
    void set_mag_world_ref(const V& v) {learned_field=v;Native::set_mag_world_ref(exact?exact_field:v);}
    H3 acceleration_jacobian() const {
        H3 H=H3::Zero();
        if constexpr(TFG) {
            H.template block<3,3>(0,0)=-skew(this->get_world_accel());
            H.template block<3,3>(0,Native::OFF_A).setIdentity();
        } else {
            for(int ax=0;ax<3;++ax) {int idx[4];float w[4];int n=this->motion().accel_terms(ax,idx,w);
                for(int k=0;k<n;++k)H(ax,Native::OFF_L+idx[k])=w[k];}
        }
        return H;
    }
    void observe_acceleration(const V& truth,float sigma) {
        const H3 H=acceleration_jacobian();
        const V r=TFG?V(this->get_world_accel()-truth):V(truth-this->get_world_accel());
        auto& P=[&]() -> M& {if constexpr(TFG)return this->P_;else return this->Pext;}();
        typename Native::MatrixNX3 PHt=P*H.transpose();
        Eigen::Matrix3f S=H*PHt+sigma*sigma*Eigen::Matrix3f::Identity();
        if constexpr(TFG) {
            typename Native::MeasDiag3 diag;
            if(!this->correct3_(r,S,PHt,sigma*sigma,false,diag))throw std::runtime_error("oracle correction rejected");
            oracle_nis=diag.nis;
        } else {
            Eigen::LDLT<Eigen::Matrix3f> ldlt(S);
            if(!ldlt.isPositive())throw std::runtime_error("oracle innovation not positive");
            typename Native::MatrixNX3 K;
            this->gain3_(PHt,ldlt,K);
            this->xext.noalias()+=K*r;
            this->template joseph_update_<3>(K,S,PHt);
            this->applyQuaternionCorrectionFromErrorState();
            oracle_nis=r.dot(ldlt.solve(r));
        }
    }
    M physical_covariance() const {if constexpr(TFG)return this->covariance_physical();else return this->covariance_full();}
    V applied_field() const {if constexpr(TFG)return this->mag_world_ref();else return this->v2ref;}
    // Finite-difference audit in the actual native error coordinates.
    void audit_oracle_jacobian() {
        const H3 H=acceleration_jacobian();
        const V a=this->get_world_accel();
        float worst=0;
        for(int j=0;j<27;++j) {
            const float eps=.001f;
            Probe changed=*this;
            if constexpr(TFG) {typename Native::Tangent d=Native::Tangent::Zero();d(j)=eps;changed.state()=this->state().retract(d);}
            else changed.xext(j)+=eps;
            V num=(changed.get_world_accel()-a)/eps;
            worst=std::max(worst,(num-H.col(j)).norm());
        }
        if(worst>.015f)throw std::runtime_error("oracle Jacobian finite difference failed");
    }
};
template<class Native,bool TFG> class Adapter final:public Experiment {
    using Core=Probe<Native,TFG>;
    using Fusion=ocean_imu::carma_shared_modes::SharedModesFusion<Core>;
    using Plain=ocean_imu::carma_shared_modes::SharedModesFusion<Native>;
    Fusion f;
    std::unique_ptr<Plain> plain;
    Options opt;
    bool injected=false,audited=false;
public:
    Adapter(const Options& o,const V& B):opt(o) {
        typename Fusion::Config cfg;carma3_bench::fill_config(cfg,3.f);
#ifndef BIAS_UNCHANGED_BASELINE
        if(o.timing!="defaults") {
            cfg.gyro_sample_period_sec=.005f;
            cfg.gyro_sample_mode=ocean_imu::kalman::filter_detail::GyroSampleMode::IntervalRate;
            if(o.timing=="both_start"||o.timing=="core_start")cfg.gyro_sample_mode=ocean_imu::kalman::filter_detail::GyroSampleMode::Instantaneous;
            if(o.timing=="both_start")cfg.frontend_gyro_sample_mode=ocean_imu::kalman::filter_detail::GyroSampleMode::Instantaneous;
            cfg.schmidt_acc_bias_hold=o.noise=="schmidt"||o.noise=="schmidt_prior";
        }
        if(o.noise=="schmidt_prior")cfg.initial_acc_bias_std=.008f*9.80665f/std::sqrt(3.f);
#endif
        cfg.mag_continuous_hard_iron=o.arm.find("cont")!=std::string::npos;
        // Synthetic sensor white noise is known here. This is a declared
        // noise-model ablation, never a default change or hardware calibration.
        if(o.noise=="gyro"||(o.noise=="all"||o.noise=="all_prior"))cfg.sigma_g=V::Constant(carma3_bench::GYR_SIGMA*std::sqrt(.005f));
        if(o.noise=="mag"||(o.noise=="all"||o.noise=="all_prior"))cfg.sigma_m=V::Constant(carma3_bench::MAG_SIGMA_UT);
        if((o.noise=="all"||o.noise=="all_prior"))cfg.sigma_a=V::Constant(carma3_bench::ACC_SIGMA);
        f.begin(cfg);auto& m=f.raw().mekf();m.exact=o.arm.find("learned")==std::string::npos&&o.arm!="correct_reference";if(m.exact||o.arm=="correct_reference")m.exact_field=B;
        if(o.noise=="prior"||o.noise=="all_prior")m.set_initial_acc_bias_std(.008f*9.80665f/std::sqrt(3.f));
        if(o.parity) {if(m.exact)throw std::runtime_error("parity requires learned arm");plain=std::make_unique<Plain>();plain->begin(cfg);}
    }
    Output step(const Input& s,bool diagnostic) override {
        auto& m=f.raw().mekf();
        if(s.t>=opt.change && !injected) {
#ifndef BIAS_UNCHANGED_BASELINE
            using Mode=ocean_imu::kalman::filter_detail::GyroSampleMode;
            const Mode live=(opt.timing=="core"||opt.timing=="core_start"||opt.timing=="both"||opt.timing=="both_start")?Mode::Instantaneous:Mode::IntervalRate;
            const Mode front=(opt.timing=="frontend"||opt.timing=="both"||opt.timing=="both_start")?Mode::Instantaneous:Mode::IntervalRate;
            if(opt.timing!="defaults") f.raw().setGyroSampleMode(live,front);
            if(opt.freeze=="sea"||opt.freeze=="both")f.raw().setSeaStateAdaptationEnabled(false);
            if(opt.freeze=="coupling"||opt.freeze=="both")m.setCouplingAdaptationEnabled(false);
#else
            if(opt.timing!="endpoint"||opt.freeze!="none")throw std::runtime_error("baseline cannot apply new interventions");
#endif
            if(opt.noise=="late_prior")m.set_initial_acc_bias_std(.008f*9.80665f/std::sqrt(3.f));
            if(opt.arm=="correct_reference"){m.exact=true;m.set_mag_world_ref(m.exact_field);}
            if(opt.injection=="bias")m.set_initial_acc_bias(m.get_acc_bias()+V(.04f,-.03f,.02f));
            if(opt.injection=="tilt")m.set_quaternion_boat(m.quaternion_boat()*Eigen::Quaternionf(Eigen::AngleAxisf(.01f,V(.6f,.8f,0))));
            if(opt.arm=="oracle_bias")m.set_initial_acc_bias(s.ba+V(.04f,-.03f,.02f));
            injected=true;
        }
        m.oracle_tick=opt.arm.rfind("oracle",0)==0 && s.t>=opt.change;
        if(m.oracle_tick){m.oracle_truth=s.aw;m.oracle_sigma=opt.oracle_sigma;}
        m.mag_called=false;
        f.update(s.dt,s.gyr,s.acc,35.f);
        if(s.mag_tick)f.updateMag(s.mag);
        if(opt.arm.rfind("oracle",0)==0 && s.t>=opt.change && s.mag_tick) {
            if(!f.isLive()||!m.acc_bias_updates_enabled()||!f.hasRefinedMagReference())throw std::runtime_error("oracle before activation");
            if(!audited){m.audit_oracle_jacobian();audited=true;}
        }
        if(plain) {
            plain->update(s.dt,s.gyr,s.acc,35.f);if(s.mag_tick)plain->updateMag(s.mag);
            const auto& p=plain->raw().mekf();
            if((p.quaternion_boat().coeffs()-m.quaternion_boat().coeffs()).norm()!=0 ||
               (p.covariance_full()-m.covariance_full()).norm()!=0 || (p.get_acc_bias()-m.get_acc_bias()).norm()!=0)
                throw std::runtime_error("inactive probe changed baseline");
        }
        Output out;out.q=m.quaternion_boat();out.bg=m.gyroscope_bias();out.ba=m.get_acc_bias();out.aw=m.get_world_accel();out.hi=f.magHardIronBodyUT();
        out.field=m.applied_field();out.live=f.isLive();out.refined=f.hasRefinedMagReference();out.enabled=m.acc_bias_updates_enabled();
        out.acc_nis=m.lastAccDiag().nis;out.mag_nis=m.lastMagDiag().nis;out.mag_applied=m.mag_called&&m.lastMagDiag().accepted;out.mag_angle=m.mag_angle;out.mag_bias=m.mag_bias;out.oracle_nis=m.oracle_nis;
        if(diagnostic)out.P=m.physical_covariance().template cast<double>();
        return out;
    }
    void write_model(std::ostream& out) override {
        const auto& m=f.raw().mekf();const auto& motion=m.motion();
        auto matrix=[&](const char* key,const auto& M,bool comma=true){out<<"\""<<key<<"\":[";for(int i=0;i<M.rows();++i){if(i)out<<',';out<<'[';for(int j=0;j<M.cols();++j){if(j)out<<',';out<<M(i,j);}out<<']';}out<<']';if(comma)out<<',';};
        Eigen::Matrix<float,18,18> W,P;W.setZero();P.setZero();motion.stationary(P,0);
        const auto& E=motion.excitation();
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)W.template block<6,6>(6*a,6*b)=E.re(a,b)*motion.block_WR()+E.im(a,b)*motion.block_WI();
        out<<std::setprecision(12)<<'{';matrix("A6",motion.block_A());matrix("W18",W);matrix("stationary18",P);
        Eigen::Matrix<float,1,3> field=m.applied_field().transpose();matrix("field",field);
        out<<"\"acc_sigma\":"<<carma3_bench::ACC_SIGMA<<",\"gyro_sigma\":"<<carma3_bench::GYR_SIGMA<<",\"mag_sigma\":"<<carma3_bench::MAG_SIGMA_UT<<'}';
    }
    void self_test() override {
        // A rejected IMU packet must not advance startup, magnetic timing,
        // frontend state or the interpolation cache.
        const float clock=f.t_;const auto before_reject=f.raw().mekf().physical_covariance();
        f.update(.005f,V::Constant(NAN),V(0,0,-9.80665f),35.f);
        if(f.t_!=clock||(before_reject-f.raw().mekf().physical_covariance()).norm()!=0)throw std::runtime_error("invalid packet advanced wrapper");
        auto m=f.raw().mekf();
        m.initialize_from_attitude(Eigen::Quaternionf(Eigen::AngleAxisf(.2f,V(.2f,.7f,.4f).normalized())),.03f,.1f);
        m.set_acc_bias_updates_enabled(true);
        if constexpr(TFG)m.state().X.col(Native::IA)=V(.2f,-.3f,.7f);
        else {auto H=m.acceleration_jacobian();m.xext+=H.transpose()*(H*H.transpose()).ldlt().solve(V(.2f,-.3f,.7f));}
        m.audit_oracle_jacobian();
        const double before=(m.get_world_accel()-V(.1f,.2f,.3f)).norm();
        m.observe_acceleration(V(.1f,.2f,.3f),.01f);
        if((m.get_world_accel()-V(.1f,.2f,.3f)).norm()>=before)throw std::runtime_error("oracle correction wrong sign");
        Eigen::MatrixXd P=m.physical_covariance().template cast<double>();
        Eigen::VectorXd d=P.diagonal().cwiseMax(1e-25).cwiseSqrt().cwiseInverse();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(d.asDiagonal()*P*d.asDiagonal());
        if(es.eigenvalues().minCoeff()<-.0001)throw std::runtime_error("oracle covariance not PSD");
    }
};
}
