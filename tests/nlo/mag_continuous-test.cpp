#define EIGEN_NON_ARDUINO
#include "nlo/TimeVarGainNLO_Adapter.h"
#include "../common/MagneticMotion.h"
struct Backend {
    TimeVarGainNloAdapter<false,NloMagType::Magnetometer,float,TrackerType::PLL> filter{};
    Backend() {filter.config().filter.expected_mag_norm=0;filter.reset();}
    V3 bias() const {return filter.gyroscopeBiasBody();}
    Q attitude() const {return filter.snapshot().q_nb;}
    bool finite() const {const auto s=filter.snapshot();return s.q_nb.coeffs().allFinite()&&s.disp_zu.allFinite()&&s.tvg.gyro_bias_b.allFinite();}
    void step(float dt,const V3& gyro,const V3& acc,const V3& mag,bool valid) {
        if(valid) filter.setMagBody(mag,true);else filter.clearMag();
        filter.update(dt,gyro,acc);
    }
};
#include "../common/MagneticContinuousReplay.h"
int main() {runContinuousMagneticTests();std::fprintf(stderr,"Continuous magnetic timing PASS\n");}
