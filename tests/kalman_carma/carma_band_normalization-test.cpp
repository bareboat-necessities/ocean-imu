// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <cstdio>
#include "kalman_carma/CarmaVerticalModel.h"
using namespace ocean_imu::kalman::carma;
namespace {
int failures = 0;
void check(bool ok, const char* message) {
    if (!ok) { ++failures; std::printf("FAIL %s\n", message); }
}
template<int N> void check_model() {
    CarmaVerticalModel<double,N> m;
    const auto shape = m.shape();
    const double fz = m.base_zero_crossing_hz();
    for (const auto& ratios : {std::pair{0.5,4.0}, std::pair{0.2,2.0}}) {
        const double lo = ratios.first*fz, hi = ratios.second*fz;
        const double fraction = m.acceleration_band_fraction(lo,hi);
        // Independent integration of the factored displacement PSD, with
        // four frequency powers for acceleration and the scalar band gain.
        constexpr int steps = 40000;
        const double u0 = std::log(1e-7), du = std::log(1e14)/steps;
        double band = 0.0;
        for (int i=0; i<=steps; ++i) {
            const double f = std::exp(u0+i*du), w = 6.283185307179586*f;
            const double h2 = f*f/(f*f+lo*lo)*hi*hi/(f*f+hi*hi);
            const double c = i==0 || i==steps ? 1.0 : i%2 ? 4.0 : 2.0;
            band += c*carma_unit_psd<N>(shape,1.0,f)*w*w*w*w*h2*f;
        }
        band *= 2.0*du/3.0;
        check(std::abs(fraction-band/m.base_var_a()) < 1e-9,"band covariance equals PSD quadrature");
        check(m.configure(0.3,0.04),"scaled model accepted");
        check(m.acceleration_band_fraction(lo,hi)==fraction,"normalized fraction independent of amplitude/kappa");
    }
    check(!std::isfinite(m.acceleration_band_fraction(0,1)),"zero corner rejected");
    check(!std::isfinite(m.acceleration_band_fraction(1,0.5)),"reversed corners rejected");
    check(!std::isfinite(m.acceleration_band_fraction(NAN,1)),"nonfinite corner rejected");
    if constexpr (N==4) {
        check(std::abs(1/std::sqrt(m.acceleration_band_fraction(.5*fz,4*fz))-1.30634701277)<1e-9,
              "default CARMA41 continuous-band gain");
        auto changed=shape;changed.zeta1*=2;
        const double before=m.acceleration_band_fraction(.5*fz,4*fz);
        check(m.set_shape(changed),"changed shape accepted");
        check(std::abs(m.acceleration_band_fraction(.5*fz,4*fz)-before)>1e-3,"shape change refreshes cross covariance");
    }
}
}
int main() {
    check_model<3>();check_model<4>();check_model<5>();
    std::printf("carma_band_normalization-test: %d failures\n",failures);
    return failures?1:0;
}
