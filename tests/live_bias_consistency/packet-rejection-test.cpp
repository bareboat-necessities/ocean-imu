#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
#include "DefaultSamplingChecks.h"
template<class Fusion> void rejection_check(const char* name) {
    typename Fusion::Config cfg;
    cfg.mag_delay_sec=0;
    cfg.mag_continuous_hard_iron=false;
    Fusion tested,control;
    tested.begin(cfg);control.begin(cfg);
    const float nan=std::numeric_limits<float>::quiet_NaN();
    for(int k=0;k<36000;++k) {
        const double t=k*.005;
        const M r=rotation(.1*t);
        const Eigen::Vector3f w=(.1*rate(.1*t)).cast<float>();
        const Eigen::Vector3f a=(r.transpose()*V(0,0,-9.80665)).cast<float>();
        if(k%700==0) {
            tested.update(nan,w,a,35.f);
            tested.update(0.f,w,a,35.f);
            tested.update(.005f,Eigen::Vector3f::Constant(nan),a,35.f);
            tested.update(.005f,w,Eigen::Vector3f::Constant(nan),35.f);
            tested.update(.005f,w,a,nan);
        }
        tested.update(.005f,w,a,35.f);control.update(.005f,w,a,35.f);
        if(k%8==7) {
            const Eigen::Vector3f mag=(r.transpose()*V(18,0,43)).cast<float>();
            tested.updateMag(mag);control.updateMag(mag);
        }
        auto& x=live_core(tested);auto& y=live_core(control);
        check((x.covariance_full()-y.covariance_full()).norm()==0 &&
              (body_quaternion(x).coeffs()-body_quaternion(y).coeffs()).norm()==0,
              "rejected packet advanced state, covariance or startup history");
    }
    std::cout<<name<<" malformed-packet startup/live control PASS\n";
}
int main() {
    rejection_check<ocean_imu::ou3::SeaStateFusion_OU_III<TrackerType::KALMANF>>("OU-III");
    rejection_check<ocean_imu::carma_shared_modes::SeaStateFusion_CarmaSharedModes>("CARMA shared");
}
