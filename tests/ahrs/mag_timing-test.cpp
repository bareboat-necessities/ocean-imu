#define EIGEN_NON_ARDUINO
#include "ahrs/KalmanQMEKF.h"
#include "util/MagneticRotation.h"
using namespace ocean_imu::magnetic;
// Test-only units/dip adapter reproduces the standalone compass backend.
// The actual qMEKF prediction/Jacobian/injection are used unchanged.
struct Core : QuaternionMEKF<float,true> {
    using Base=QuaternionMEKF<float,true>;
    using Quat=ocean_imu::magnetic::Q;
    static constexpr float MagneticUnitSigma=.020f;
    static constexpr float GyroStepVariance=.0030f*.0030f;
    struct Diag {bool accepted=false;} diag;
    Core():Base(V3::Constant(.06f*9.80665f),V3::Constant(.003f),V3::Constant(MagneticUnitSigma),.5f,1e-2f,1e-9f) {
        initialize_from_acc_mag(V3(0,0,-9.80665f),V3(20.3f,0,46.8f).normalized());
    }
    Quat q() const {const auto& v=quaternion();return Quat(v(3),v(0),v(1),v(2));}
    const auto& covariance_full() const {return covariance();}
    V3 gyroscope_bias_body() const {return gyroscope_bias();}
    M3 gyroscope_bias_covariance_body() const {return covariance().template bottomRightCorner<3,3>();}
    const Diag& lastMagDiag() const {return diag;}
    void setDip(const V3& unit) {
        const V3 world=q()*unit;
        set_mag_world_ref(V3(std::hypot(world.x(),world.y()),0,world.z()));
    }
    void measurement_update_mag_only(const V3& field) {
        const V3 unit=field.normalized();setDip(unit);
        Base::measurement_update_mag_only(unit);diag.accepted=true;
    }
    void measurement_update_mag_only(const V3& field,const M3& covariance_body) {
        diag.accepted=false;
        if(!ocean_imu::validMeasurementCovariance(covariance_body)) return;
        const V3 unit=field.normalized();setDip(unit);
        diag.accepted=Base::measurement_update_mag_only(unit,covariance_body/field.squaredNorm());
    }
};
constexpr bool kRequireOu3Recovery=false;
#include "../common/MagneticTimingReplay.h"
static Core core() {return Core{};}
static Q attitude(const Core& f) {return f.q();}
static void advanceCore(Core& f,const V3& gyro,const V3& acc,float dt,int) {
    f.time_update(gyro,dt);
    f.measurement_update_acc_only(V3(acc.normalized()*g_std));
}
int main() {
    covarianceCompatibility();runPairedMagneticTests();
    std::fprintf(stderr,"Compass qMEKF magnetic timing PASS\n");
}
