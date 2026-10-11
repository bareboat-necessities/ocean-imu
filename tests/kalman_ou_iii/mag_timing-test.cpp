// Actual OU-III core; shared fixture never ships on the device.
#define EIGEN_NON_ARDUINO
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
using Core = ocean_imu::kalman::Kalman3D_Wave_OU_III<float>;
constexpr bool kRequireOu3Recovery = true;
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
static void advanceCore(Core& f,const V3& gyro,const V3& acc,float dt,int index) {
    f.time_update(gyro,dt);
    if(index%10==0) f.applyIntegralZeroPseudoMeas();
    f.measurement_update_acc_only(acc);
}
static void compatibility() {
    covarianceCompatibility();
    auto a=std::make_unique<Core>(core()), b=std::make_unique<Core>(core());
    // Wind heel must transform both the vector and its full covariance.
    b=std::make_unique<Core>(*a); a->update_wind_heel(.2f); b->update_wind_heel(.2f);
    a->measurement_update_mag_only(kField);
    b->measurement_update_mag_only(kField,rm);
    check((a->covariance_full()-b->covariance_full()).norm()<2e-6f,"isotropic covariance deheel parity");

}
int main() {
    compatibility();
    magneticFacadeIdentity<SeaStateFusion_OU_III<TrackerType::KALMANF>>(
        [](auto& f)->auto& {return f.raw().mekf();},
        [](auto& f) {return f.raw().startupProxyQuat();},
        [](auto& f) {return f.attitudeQuat();},
        [](auto& f) {return f.hasRefinedMagReference();});
    runPairedMagneticTests();
    std::fprintf(stderr,"OU-III magnetic timing regression PASS\n");
}
