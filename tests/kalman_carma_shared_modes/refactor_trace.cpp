// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// The same source is compiled against pristine main and the candidate headers.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <random>
#if REFACTOR_FILTER == 0
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
using Fusion=ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF,4,false,4,false>;
#elif REFACTOR_FILTER == 1
#include "kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h"
using Fusion=ocean_imu::carma::SeaStateFusion_TFG_CARMA3<TrackerType::KALMANF,4>;
#elif REFACTOR_FILTER == 2
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
using Fusion=ocean_imu::carma_shared_modes::SeaStateFusion_CarmaSharedModes;
#else
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
using Fusion=SeaStateFusion_OU_III<TrackerType::KALMANF>;
#endif
const float g_std=9.80665f;
int main(){
    Fusion f;Fusion::Config c;c.with_mag=true;c.mag_delay_sec=20;
    c.sigma_a=Eigen::Vector3f::Constant(.0294f);c.sigma_g=Eigen::Vector3f::Constant(.000031f);c.sigma_m=Eigen::Vector3f::Constant(7.7f);
    f.begin(c);
#if REFACTOR_FILTER < 2
    f.raw().setCarmaExcitationMode(ocean_imu::carma::CarmaExcitationMode::Evidence);
#endif
    f.raw().enableTuner(true);f.raw().enableClamp(true);
    std::mt19937 rng(61120);std::normal_distribution<float> n(0,1);
    using V=Eigen::Vector3f;const V bw(18,1.5,43);
    for(int k=0;k<120000;++k){const float t=k*.005f;const float roll=.04f*std::sin(.7f*t),pitch=.03f*std::sin(.9f*t),yaw=.1f*t;
        const float rd=.028f*std::cos(.7f*t),pd=.027f*std::cos(.9f*t),yd=.1f;
        const Eigen::Matrix3f R=(Eigen::AngleAxisf(yaw,V::UnitZ())*Eigen::AngleAxisf(pitch,V::UnitY())*Eigen::AngleAxisf(roll,V::UnitX())).toRotationMatrix();
        const float amp=t<200?1.f:t<400?3.f:.3f;
        V a(amp*.2f*std::sin(1.6f*t),amp*.1f*std::cos(1.6f*t),amp*.3f*std::cos(1.6f*t));
        V gyr(rd-yd*std::sin(pitch),pd*std::cos(roll)+yd*std::sin(roll)*std::cos(pitch),-pd*std::sin(roll)+yd*std::cos(roll)*std::cos(pitch));
        V acc=R.transpose()*(a-V(0,0,9.80665f));
        for(int i=0;i<3;++i){acc(i)+=.01f*n(rng)+.001f*(i+1);gyr(i)+=.0001f*n(rng);}
        f.update(.005f,gyr,acc,35.f);
        if(k%2==0 && !(t>=300&&t<400)){V m=R.transpose()*bw;for(int i=0;i<3;++i)m(i)+=.3f*n(rng);f.updateMag(m);}
        if(k%20==0){const auto& m=f.raw().mekf();const auto p=m.get_position(),v=m.get_velocity(),bg=m.gyroscope_bias(),ba=m.get_acc_bias();const auto q=m.quaternion_boat();
            std::printf("%d,%d",k,f.isLive()?1:0);
            for(auto value:{p.x(),p.y(),p.z(),v.x(),v.y(),v.z(),q.w(),q.x(),q.y(),q.z(),bg.x(),bg.y(),bg.z(),ba.x(),ba.y(),ba.z(),m.lastAccDiag().nis})std::printf(",%.9g",double(value));
            std::printf("\n");}
    }
}
