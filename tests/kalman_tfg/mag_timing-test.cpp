#define EIGEN_NON_ARDUINO
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
using Core = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
constexpr bool kRequireOu3Recovery = false;
#include "../common/MagneticTimingReplay.h"
#include "../common/MagneticFacadeIdentity.h"
static Core core() {
    Core f(.00135f);
    f.set_Racc_std(V3::Constant(.12f)); f.set_Rmag_std(V3::Constant(.8f));
    f.set_magnetic_reference_world(kField);
    f.set_acc_bias_updates_enabled(false);
    return f;
}
static Q attitude(const Core& f) { return f.quaternion(); }
static void advanceCore(Core& f,const V3& gyro,const V3& acc,float dt,int index) {
    f.time_update(gyro,dt);
    f.measurement_update_acc_only(acc);
    if(index%10==0) f.applyIntegralZeroPseudoMeas(); // TFG: acc precedes S.
}
int main() {
    covarianceCompatibility();
    magneticFacadeIdentity<ocean_imu::tfg::SeaStateFusionFilter_TFG<>>(
        [](auto& f)->auto& {return f.mekf();},
        [](auto& f) {return f.startupProxyQuat();},
        [](auto& f) {return f.quaternion();},
        [](auto& f) {return f.magReferenceRefined();});
    runPairedMagneticTests();
    std::fprintf(stderr,"Kalman3D_Wave_TFG magnetic timing regression PASS\n");
}
