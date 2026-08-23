// Pins only the deployed OU-III SpectralMSE integral regularizer law.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
static bool rel_near(float a,float b,float r){return std::fabs(a-b)<=r*std::max(1e-12f,std::max(std::fabs(a),std::fabs(b)));}
static int fail(const char* m){std::cerr<<"FAIL: "<<m<<"\n";return 1;}

int main(){
    Filter f(false);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),Eigen::Vector3f::Constant(0.00157f),Eigen::Vector3f::Constant(0.25f));
    if(f.getRSLaw()!=RSAdaptationLaw::SpectralMSE) return fail("OU-III default is not SpectralMSE");

    const float tau=2.179f, sigma=0.724f, d=1.02f;
    const float r0=f.rs_spectral_mse_target_(tau,sigma);
    const float pa=std::log(f.rs_spectral_mse_target_(tau,sigma*d)/r0)/std::log(d);
    const float pt=std::log(f.rs_spectral_mse_target_(tau*d,sigma)/r0)/std::log(d);
    if(!rel_near(pa,6.0f/7.0f,2e-3f)) return fail("SpectralMSE amplitude exponent is not 6/7");
    if(!rel_near(pt,41.0f/14.0f,3e-3f)) return fail("SpectralMSE cadence-normalized tau exponent is not 41/14");

    const float before=f.rs_spectral_mse_target_(tau,sigma);
    f.setSigmaCoeff(2.0f*f.sigma_coeff_);
    const float after=f.rs_spectral_mse_target_(tau,2.0f*sigma);
    if(!rel_near(before,after,1e-4f)) return fail("SpectralMSE changed at fixed physical acceleration RMS");

    if(!rel_near(f.getRSMseCoeff(),R_S_MSE_COEFF_DEFAULT,1e-6f)) return fail("SpectralMSE analytical coefficient moved");
    std::cout<<"OU-III deployed SpectralMSE law passed\n";
    return 0;
}
