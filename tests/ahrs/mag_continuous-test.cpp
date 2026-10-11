#define EIGEN_NON_ARDUINO
#include "ahrs/Mahony_AHRS.h"
#include "../common/MagneticMotion.h"
struct Backend {
    Mahony_AHRS<float> filter{};
    Backend() { mahony_AHRS_init(&filter,12.0f,.0004f); }
    V3 bias() const {return V3(-filter.integralFBx,-filter.integralFBy,-filter.integralFBz);}
    Q attitude() const {return Q(filter.q0,filter.q1,filter.q2,filter.q3);}
    bool finite() const {return attitude().coeffs().allFinite();}
    void step(float dt,const V3& gyro,const V3& acc,const V3& mag,bool valid) {
        const V3 a=-acc.normalized()*g_std;
        float pitch,roll,yaw;
        if(valid) filter.updateMag(gyro.x(),gyro.y(),gyro.z(),a.x(),a.y(),a.z(),
            mag.x(),mag.y(),mag.z(),&pitch,&roll,&yaw,dt);
        else filter.update(gyro.x(),gyro.y(),gyro.z(),a.x(),a.y(),a.z(),&pitch,&roll,&yaw,dt);
    }
};
#include "../common/MagneticContinuousReplay.h"
int main() {runContinuousMagneticTests();std::fprintf(stderr,"Continuous magnetic timing PASS\n");}
