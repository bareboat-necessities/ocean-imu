#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;
static bool near(float a,float b,float e=1e-6f){return std::fabs(a-b)<=e*std::max(1.0f,std::max(std::fabs(a),std::fabs(b)));}

int main(){
    Filter f(false);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),Eigen::Vector3f::Constant(0.00157f),Eigen::Vector3f::Constant(0.25f));
    const Eigen::Vector3f gyro=Eigen::Vector3f::Zero();
    const Eigen::Vector3f acc(0,0,-g_std);
    f.initialize_from_acc(acc);
    f.startup_stage_=Filter::StartupStage::Live;
    f.mekf_->set_linear_block_enabled(true);
    f.time_=1.0; f.last_adapt_time_sec_=0.0; f.adapt_every_secs_=0.05f;

    if(f.getPseudoLaw()!=PseudoAdaptationLaw::PhysicalMSE){std::cerr<<"FAIL: OU-II default is not PhysicalMSE\n";return 1;}
    if(!near(f.getPseudoUpdateTauRatio(),0.015f/1.1f)||!near(f.getPseudoUpdatePeriodSec(),0.015f)){std::cerr<<"FAIL: nominal cadence moved\n";return 1;}
    if(!near(f.getAdaptationSeaPeriods(),0.40f)){std::cerr<<"FAIL: sea-scaled EMA moved\n";return 1;}

    const float active_tau=f.mekf_->tau_aw;
    const float active_rp=f.mekf_->R_p0(2,2);
    const float active_rv=f.mekf_->R_v0(2,2);
    f.adapt_mekf(0.1f,3.0f,0.8f,1.2f,0.7f,2.5f);
    const float staged_tau=f.tune_.tau_applied;
    const float staged_rp=std::min(std::max(f.tune_.R_p0_std_applied,f.MIN_R_p0_std_),f.MAX_R_p0_std_);
    const float staged_rv=std::min(std::max(f.tune_.R_v0_std_applied,f.MIN_R_v0_std_),f.MAX_R_v0_std_);
    const float staged_period=std::min(std::max(f.pseudo_update_tau_ratio_*staged_tau,f.pseudo_update_period_min_s_),f.pseudo_update_period_max_s_);

    if(!f.online_tune_apply_pending_){std::cerr<<"FAIL: adaptation was not staged\n";return 1;}
    if(!near(f.mekf_->tau_aw,active_tau)||!near(f.mekf_->R_p0(2,2),active_rp)||!near(f.mekf_->R_v0(2,2),active_rv)){std::cerr<<"FAIL: staging changed active schedule\n";return 1;}

    f.enable_tuner_=false;
    f.updateTime(0.005f,gyro,acc);
    if(!near(f.mekf_->tau_aw,staged_tau)||!near(f.mekf_->R_p0(2,2),staged_rp*staged_rp)||!near(f.mekf_->R_v0(2,2),staged_rv*staged_rv)||!near(f.getPseudoUpdatePeriodSec(),staged_period)){std::cerr<<"FAIL: PhysicalMSE staged schedule was not committed\n";return 1;}
    if(!near(f.mekf_->R_p0(2,2)*f.getPseudoUpdatePeriodSec(),staged_rp*staged_rp*staged_period)||!near(f.mekf_->R_v0(2,2)*f.getPseudoUpdatePeriodSec(),staged_rv*staged_rv*staged_period)){std::cerr<<"FAIL: PhysicalMSE information-rate product moved\n";return 1;}

    f.setTauScaledPseudoUpdateCadence(false);
    if(!near(f.getPseudoUpdatePeriodSec(),PSEUDO_UPDATE_PERIOD_NOMINAL_S)){std::cerr<<"FAIL: fixed-cadence ablation did not restore nominal cadence\n";return 1;}
    f.setTauScaledPseudoUpdateCadence(true);
    if(!near(f.getPseudoUpdatePeriodSec(),staged_period)){std::cerr<<"FAIL: tau-scaled cadence was not restored\n";return 1;}

    std::cout<<"OU-II PhysicalMSE scheduling passed\n";
    return 0;
}
