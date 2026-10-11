#pragma once
// Thin adapters: the algorithms, gain-row holds, bias states and covariance
// resets under test are the production implementations.
#define EIGEN_NON_ARDUINO
#include "util/AttitudeProtection.h"
#include <Eigen/Eigenvalues>
using V = Eigen::Vector3f;
using Q = Eigen::Quaternionf;
constexpr float gravity = 9.80665f;
extern const float g_std = gravity;
const V earth(20.3f, 0, 46.8f);

#if RAPID_BACKEND == 1
#include "ahrs/KalmanQMEKF.h"
struct Backend {
    QuaternionMEKF<float,true> f{V::Constant(.06f*gravity),V::Constant(.003f),V::Constant(.020f),.5f,.01f,1e-9f};
    Backend() { f.initialize_from_acc_mag(V(0,0,-gravity),earth.normalized()); }
    Q attitude() const { const auto& q=f.quaternion(); return Q(q(3),q(0),q(1),q(2)); }
    V bias() const { return f.gyroscope_bias(); }
    void recover(const Q& q) { f.reacquire_attitude(q); }
    bool finite() const { return f.covariance().allFinite(); }
    void step(float dt,const V& w,const V& a,const V& m,bool mag,bool fresh,bool aiding,bool accel) {
        f.time_update(w,dt);
        if (aiding && accel) f.measurement_update_acc_only(V(a.normalized()*gravity));
        if (aiding && mag && fresh) {
            const V u=m.normalized(), world=attitude()*u;
            f.set_mag_world_ref(V(std::hypot(world.x(),world.y()),0,world.z()));
            f.measurement_update_mag_only(u);
        }
    }
};
#elif RAPID_BACKEND == 2
#include "ahrs/Mahony_AHRS.h"
struct Backend {
    Mahony_AHRS<float> f{};
    Backend() { f.init(12.0f,.0004f); }
    Q attitude() const { return Q(f.q0,f.q1,f.q2,f.q3).normalized(); }
    V bias() const { return V(-f.integralFBx,-f.integralFBy,-f.integralFBz); }
    void recover(const Q& q) { f.setQuaternion(q.w(),q.x(),q.y(),q.z()); }
    bool finite() const { return attitude().coeffs().allFinite(); }
    void step(float dt,const V& w,const V& a,const V& m,bool mag,bool,bool aiding,bool) {
        f.setAidingEnabled(aiding); const V down=-a; float p,r,y;
        if (mag) f.updateMag(w.x(),w.y(),w.z(),down.x(),down.y(),down.z(),m.x(),m.y(),m.z(),&p,&r,&y,dt);
        else f.update(w.x(),w.y(),w.z(),down.x(),down.y(),down.z(),&p,&r,&y,dt);
    }
};
#elif RAPID_BACKEND == 3 || RAPID_BACKEND == 4
#if RAPID_BACKEND == 3
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
using Core = ocean_imu::kalman::Kalman3D_Wave_OU_II<float>;
using Fusion = SeaStateFusion_OU_II<TrackerType::KALMANF>;
#else
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
using Core = ocean_imu::kalman::Kalman3D_Wave_OU_III<float>;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
#endif
struct Backend {
    Core f{V::Constant(.12f),V::Constant(.00135f),V::Constant(.8f)};
    Backend() { f.initialize_from_attitude(Q::Identity(),.03f,.08f); f.set_mag_world_ref(earth); }
    Q attitude() const { return f.quaternion_boat(); }
    V bias() const { return f.gyroscope_bias_body(); }
    void recover(const Q& q) { f.initialize_from_attitude(q,.05f,.05f); }
    bool finite() const { return f.covariance_full().allFinite() && f.get_position().allFinite(); }
    void step(float dt,const V& w,const V& a,const V& m,bool mag,bool fresh,bool aiding,bool accel) {
        f.set_attitude_aiding_enabled(aiding); f.set_accelerometer_fresh(accel);
        f.time_update(w,dt); f.measurement_update_acc_only(a);
        if (mag && fresh) f.measurement_update_mag_only(m);
    }
};
#elif RAPID_BACKEND == 5
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
struct Backend {
    ocean_imu::kalman::Kalman3D_Wave_TFG<float> f{.00135f};
    Backend() { f.set_Racc_std(V::Constant(.12f)); f.set_Rmag_std(V::Constant(.8f)); f.set_magnetic_reference_world(earth); }
    Q attitude() const { return f.quaternion(); }
    V bias() const { return f.gyroscope_bias_body(); }
    void recover(const Q& q) { f.reacquire_attitude(q); }
    bool finite() const { return f.covariance_full().allFinite() && f.get_position().allFinite(); }
    void step(float dt,const V& w,const V& a,const V& m,bool mag,bool fresh,bool aiding,bool accel) {
        f.set_attitude_aiding_enabled(aiding); f.set_accelerometer_fresh(accel);
        f.time_update(w,dt); f.measurement_update_acc_only(a); f.applyIntegralZeroPseudoMeas();
        if (mag && fresh) f.measurement_update_mag_only(m);
    }
};
#elif RAPID_BACKEND == 6
#include "nlo/TimeVarGainNLO_Adapter.h"
struct Backend {
    TimeVarGainNloAdapter<false,NloMagType::Magnetometer,float,TrackerType::PLL> f{};
    Backend() { f.config().filter.expected_mag_norm=0; f.reset(); f.setQuaternionBodyToNED(Q::Identity()); }
    Q attitude() const { return f.snapshot().q_nb; }
    V bias() const { return f.gyroscopeBiasBody(); }
    void recover(const Q& q) { f.reacquireAttitude(q); }
    bool finite() const { return f.snapshot().disp_zu.allFinite(); }
    void step(float dt,const V& w,const V& a,const V& m,bool mag,bool,bool aiding,bool) {
        f.setAttitudeAidingEnabled(aiding);
        if (mag && aiding) f.setMagBody(m,true); else f.clearMag();
        f.update(dt,w,a);
    }
};
#elif RAPID_BACKEND == 7
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
struct Backend {
    marine_obs::AdaptiveVerticalPIIMahony<float,true,TrackerType::PLL> f{};
    static V map(const V& a) { return V(a.x(),-a.y(),-a.z()); }
    Q attitude() const { const auto& m=f.mahonyState(); return Q(m.q0,m.q1,-m.q2,-m.q3).normalized(); }
    V bias() const { const auto& m=f.mahonyState(); return map(V(-m.integralFBx,-m.integralFBy,-m.integralFBz)); }
    void recover(const Q& q) { f.reacquireAttitude(q.w(),q.x(),-q.y(),-q.z()); }
    bool finite() const { return std::isfinite(f.displacement()); }
    void step(float dt,const V& w,const V& a,const V& m,bool mag,bool,bool aiding,bool) {
        f.mahonyState().setAidingEnabled(aiding); const V g=map(w),acc=map(a),field=map(m);
        if (mag) f.updateIMUMag(g.x(),g.y(),g.z(),acc.x(),acc.y(),acc.z(),field.x(),field.y(),field.z(),dt);
        else f.updateIMU(g.x(),g.y(),g.z(),acc.x(),acc.y(),acc.z(),dt);
    }
};
#endif
