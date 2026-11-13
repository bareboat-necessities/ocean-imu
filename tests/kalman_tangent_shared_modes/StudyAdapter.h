#pragma once
#include "Study.h"
#include "../kalman_carma/carma3/Bench3.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
#include "lie/SO3Jacobians.h"
#include "../common/CarmaModelExport.h"
#include <Eigen/Eigenvalues>
#include <type_traits>
#include <stdexcept>
namespace tangent_study {
template<class Native> struct RefCore:Native {
    using Native::Native;
    bool override_field=false,acc_called=false,mag_called=false;
    V fixed_field=V::UnitX(),last_field=V::UnitX();
    Eigen::Matrix3f acc_R=Eigen::Matrix3f::Identity(),mag_R=Eigen::Matrix3f::Identity();
    void set_mag_world_ref(const V& b) { last_field=override_field?fixed_field:b;Native::set_mag_world_ref(last_field); }
    void measurement_update_acc_only(const V& a,float temperature=Native::tempC_ref) {
        acc_called=true;acc_R=this->quaternion_boat().toRotationMatrix();Native::measurement_update_acc_only(a,temperature);
    }
    void measurement_update_mag_only(const V& m) {
        mag_called=true;mag_R=this->quaternion_boat().toRotationMatrix();Native::measurement_update_mag_only(m);
    }
};
template<class Native,bool TFG,bool Shared=true> class Adapter final:public Filter {
    using Core=RefCore<Native>;
    using Fusion=std::conditional_t<Shared,ocean_imu::carma_shared_modes::SharedModesFusion<Core>,ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF,4,false,4,false,Core>>;
    Fusion f_;
    Config cfg_;
    bool shocked_=false,gate_activated_=false;
public:
    explicit Adapter(const Config& c):cfg_(c) {
        typename Fusion::Config fc;carma3_bench::fill_config(fc,3.f);
        fc.mag_continuous_hard_iron=false;
        fc.gyro_sample_period_sec=float(c.dt);
        if(c.population_prior) fc.initial_acc_bias_std=.008f*9.80665f/std::sqrt(3.f);
        if(c.matched_noise) {
            fc.sigma_a=V::Constant(carma3_bench::ACC_SIGMA);
            fc.sigma_g=V::Constant(carma3_bench::GYR_SIGMA*std::sqrt(float(c.dt)));
            fc.sigma_m=V::Constant(carma3_bench::MAG_SIGMA_UT);
        }
        f_.begin(fc);
        if constexpr(!Shared) f_.raw().setCarmaExcitationMode(ocean_imu::carma::CarmaExcitationMode::Evidence);
        auto& m=f_.raw().mekf();
#ifndef TANGENT_BASELINE
        if constexpr(TFG) m.set_tangent_bias_geometry(c.tangent_geometry);
#endif
        m.override_field=c.reference!="learned";m.fixed_field=c.field;
        if(c.reference=="wrong") m.fixed_field=Eigen::AngleAxisf(.06f,V::UnitY())*c.field;
        if(m.override_field) m.set_mag_world_ref(c.field);
    }
    void step(const Sample& s) override {
        auto& m=f_.raw().mekf();
        if constexpr(Shared) {
#ifndef TANGENT_BASELINE
            if(cfg_.gate && !gate_activated_ && s.t>=cfg_.gate_start) {
                auto c=m.sharedModesConfig().bias_learning;c.enabled=true;m.configureBiasLearning(c);gate_activated_=true;
            }
#endif
        }
        if(!shocked_ && cfg_.shock_time>=0 && s.t>=cfg_.shock_time && f_.isLive()) {
            if(cfg_.shock=="bias") m.set_initial_acc_bias(m.get_acc_bias()+V(.06f,-.04f,.025f));
            if(cfg_.shock=="attitude") m.set_quaternion_boat(m.quaternion_boat()*Eigen::Quaternionf(Eigen::AngleAxisf(.12f,V(.6f,.8f,0))));
            if(cfg_.shock=="gyro") {
                if constexpr(TFG) m.state().B.col(0)+=V(.002f,-.0015f,.001f);
                // A gyro error is injected into the SENSOR history in this
                // arm, so this estimator-state option is TFG-only and unused.
            }
            shocked_=true;
        }
        m.acc_called=false;m.mag_called=false;
        f_.update(s.dt,s.gyr,s.acc,35.f);
        if(s.mag_tick) f_.updateMag(s.mag);
    }
    Output output(bool covariance) override {
        Output o;const auto& m=f_.raw().mekf();
        o.R=m.quaternion_boat().template cast<double>().toRotationMatrix();
        o.p=m.get_position().template cast<double>();o.v=m.get_velocity().template cast<double>();o.a=m.get_world_accel().template cast<double>();
        o.bg=m.gyroscope_bias().template cast<double>();o.ba=m.get_acc_bias().template cast<double>();o.field=m.last_field.template cast<double>();
        o.live=f_.isLive();o.bias_enabled=m.acc_bias_updates_enabled();
        const auto& ad=m.lastAccDiag();const auto& md=m.lastMagDiag();
        o.anis=m.acc_called?ad.nis:NAN;o.mnis=m.mag_called?md.nis:NAN;o.acc_applied=m.acc_called&&ad.accepted;o.mag_applied=m.mag_called&&md.accepted;
        o.ar=ad.r.template cast<double>();o.aS=ad.S.template cast<double>();o.mr=md.r.template cast<double>();o.mS=md.S.template cast<double>();
        if constexpr(TFG) {
            const Eigen::Matrix3d A=m.acc_R.template cast<double>(),M=m.mag_R.template cast<double>();
            o.ar=A.transpose()*o.ar;o.aS=A.transpose()*o.aS*A;
            o.mr=M.transpose()*o.mr;o.mS=M.transpose()*o.mS*M;
        }
#ifndef TANGENT_BASELINE
        if constexpr(Shared) {o.gate_eigen=m.biasLearning().information_eigenvalues().template cast<double>();o.gate_trace=m.biasLearning().scale().trace();}
#endif
        o.finite=o.R.allFinite()&&o.p.allFinite()&&o.v.allFinite()&&o.ba.allFinite()&&o.bg.allFinite();
        if(covariance) {
            constexpr int N=Native::NX;
            const Eigen::Matrix<double,N,N> P=[&]() {if constexpr(TFG)return m.covariance_physical().template cast<double>();else return m.covariance_full().template cast<double>();}();
            Eigen::Matrix<double,15,N> J=Eigen::Matrix<double,15,N>::Zero();
            J.template block<3,3>(0,0)=o.R;
            for(int ax=0;ax<3;++ax) {
                J(3+ax,Native::pos_index(ax))=-1;J(6+ax,Native::vel_index(ax))=-1;
                J(9+ax,3+ax)=-1;J(12+ax,Native::OFF_BA+ax)=-1;
            }
            o.P=J*P*J.transpose();o.asym=(P-P.transpose()).cwiseAbs().maxCoeff();
            Eigen::VectorXd sd=P.diagonal().cwiseMax(1e-25).cwiseSqrt().cwiseInverse();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(sd.asDiagonal()*P*sd.asDiagonal());
            o.min_eig=es.info()==Eigen::Success?es.eigenvalues().minCoeff():-1;
            o.finite=o.finite&&P.allFinite()&&o.min_eig>-.0001;
        }
        return o;
    }
    std::size_t memory() const override {return sizeof(*this)+sizeof(Core);}
    int states() const override {return Native::NX;}
    void write_model(std::ostream& out) const override {
        if constexpr(Shared) {
            const auto& m=f_.raw().mekf();
            ocean_imu::sim::write_carma6_model(out,m.motion(),m.last_field,
                carma3_bench::ACC_SIGMA,carma3_bench::GYR_SIGMA,carma3_bench::MAG_SIGMA_UT);
        } else out<<"{}";
    }
};
}
