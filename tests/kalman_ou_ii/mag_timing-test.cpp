#define EIGEN_NON_ARDUINO
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
using Core = ocean_imu::kalman::Kalman3D_Wave_OU_II<float>;
constexpr bool kRequireOu3Recovery = false;
#include "../common/MagneticTimingReplay.h"
#include "../common/MagneticFacadeIdentity.h"
static Core core() {
    Core f(V3::Constant(.12f),V3::Constant(.00135f),V3::Constant(.8f));
    f.set_mag_world_ref(kField);
    f.initialize_from_attitude(Q::Identity(),.03f,.08f);
    f.set_acc_bias_updates_enabled(false);
    return f;
}
static Q attitude(const Core& f) { return f.quaternion_boat(); }
static void advanceCore(Core& f,const V3& gyro,const V3& acc,float dt,int) {
    f.time_update(gyro,dt); // OU-II owns its original pseudo-update schedule.
    f.measurement_update_acc_only(acc);
}
int main() {
    covarianceCompatibility();
    magneticFacadeIdentity<SeaStateFusion_OU_II<TrackerType::KALMANF>>(
        [](auto& f)->auto& {return f.raw().mekf();},
        [](auto& f) {return f.raw().startupProxyQuat();},
        [](auto& f) {return f.attitudeQuat();},
        [](auto& f) {return f.hasRefinedMagReference();});
    runPairedMagneticTests();
    std::fprintf(stderr,"Kalman3D_Wave_OU_II magnetic timing regression PASS\n");
}
