#include "SelectiveMag.h"
#include <cassert>
#include <iostream>

int main() {
    using Mat = Eigen::Matrix3f;
    using Vec = Eigen::Vector3f;
    SelectiveMag calm, steering, waves;
    Vec f(21,0,47);
    for (int i=0;i<2000;++i) {
        calm.observe(Mat::Identity(),0.04f*Mat::Identity(),f);
        const float t=static_cast<float>(i)*0.04f;
        const Mat yaw=Eigen::AngleAxisf(0.001f*std::cos(t*.3f),Vec::UnitZ()).toRotationMatrix();
        steering.observe(yaw,0.04f*Mat::Identity(),f);
        const Mat turn=Eigen::AngleAxisf(.008f*std::sin(t),Vec::UnitX()).toRotationMatrix()*
                       Eigen::AngleAxisf(.01f*std::cos(.5f*t),Vec::UnitY()).toRotationMatrix()*yaw;
        f=turn*f;
        waves.observe(turn,0.04f*Mat::Identity(),f);
    }
    const auto eig=waves.information.selfadjointView<Eigen::Lower>().eigenvalues();
    assert(calm.information.norm()<1e-8);
    assert(std::abs(steering.hard)<1e-8);
    assert(eig.minCoeff()>0);
    assert(waves.hard>calm.hard);
    std::cout<<"stationary_norm="<<calm.information.norm()<<" yaw_only_min="<<steering.hard
             <<" mixed_min="<<eig.minCoeff()<<"\n";
}
