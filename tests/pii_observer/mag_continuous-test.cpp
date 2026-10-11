#define EIGEN_NON_ARDUINO
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
#include "../common/MagneticMotion.h"
struct Backend {
    marine_obs::AdaptiveVerticalPIIMahony<float,true,TrackerType::PLL> filter{};
    static V3 map(const V3& v) {return V3(v.x(),-v.y(),-v.z());}
    V3 bias() const {const auto& m=filter.mahonyState();return map(V3(-m.integralFBx,-m.integralFBy,-m.integralFBz));}
    Q attitude() const {const auto& m=filter.mahonyState();return Q(m.q0,m.q1,-m.q2,-m.q3);}
    bool finite() const {return attitude().coeffs().allFinite()&&std::isfinite(filter.displacement());}
    void step(float dt,const V3& gyro,const V3& acc,const V3& mag,bool valid) {
        const V3 w=map(gyro),a=map(acc),m=map(mag);
        if(valid) filter.updateIMUMag(w.x(),w.y(),w.z(),a.x(),a.y(),a.z(),m.x(),m.y(),m.z(),dt);
        else filter.updateIMU(w.x(),w.y(),w.z(),a.x(),a.y(),a.z(),dt);
    }
};
#include "../common/MagneticContinuousReplay.h"
int main() {runContinuousMagneticTests();std::fprintf(stderr,"Continuous magnetic timing PASS\n");}
