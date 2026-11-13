#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#include <stdexcept>
#include <iostream>
using V = Eigen::Vector3f;
using Mode = ocean_imu::kalman::filter_detail::GyroSampleMode;
template<class Core> void check_reinitialization() {
    ocean_imu::carma::SeaStateFusionFilter_CARMA3<TrackerType::KALMANF,4,false,4,false,Core> raw;
    const V a=V::Constant(.03f),g=V::Constant(.03f),m=V::Constant(.8f);
    raw.setGyroSampleMode(Mode::Instantaneous,Mode::Instantaneous);
    raw.setGyroSamplePeriod(.005f);
    for(int extended=0;extended<2;++extended) {
        if(extended) raw.initialize_ext(a,g,m,.1f,.01f,0.f,.01f,9.80665f);
        else raw.initialize(a,g,m);
        const auto before=raw.mekf().covariance_full().eval();
        const auto attitude=raw.mekf().quaternion_boat();
        raw.mekf().time_update(V(.2f,-.3f,.4f),.01f);
        if((before-raw.mekf().covariance_full()).norm()!=0 ||
           (attitude.coeffs()-raw.mekf().quaternion_boat().coeffs()).norm()!=0)
            throw std::runtime_error("reinitialization lost instantaneous-sample priming");
        Core expected=raw.mekf();
        expected.set_gyro_sample_mode(Mode::IntervalRate);
        expected.set_gyro_noise_density_rad_sqrt_s(g*std::sqrt(3.f));
        raw.mekf().time_update(V(.4f,.1f,-.2f),.015f);
        expected.time_update(V(.3f,-.1f,.1f),.015f);
        if((expected.covariance_full()-raw.mekf().covariance_full()).norm()>1e-6f ||
           (expected.quaternion_boat().coeffs()-raw.mekf().quaternion_boat().coeffs()).norm()>1e-6f)
            throw std::runtime_error("reinitialization lost sampling mode or exposure period");
    }
    raw.setGyroSampleMode(Mode::IntervalRate,Mode::IntervalRate);
    raw.initialize(a,g,m);
    const auto before=raw.mekf().quaternion_boat();
    raw.mekf().time_update(V(.2f,-.3f,.4f),.01f);
    if((before.coeffs()-raw.mekf().quaternion_boat().coeffs()).norm()<1e-4f)
        throw std::runtime_error("reinitialization did not preserve interval mode");
}
int main() {
    check_reinitialization<ocean_imu::carma_shared_modes::KalmanCarmaSharedModes<float>>();
    check_reinitialization<ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>>();
    std::cout<<"core mode and exposure survive both reinitialization APIs PASS\n";
}
