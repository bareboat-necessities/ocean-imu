#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy
// Adversarial means are injected only by this regression, never by a replay.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <limits>
#define private public
#if defined(GYRO_PROJECTION_OU2)
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
using GyroFilter = Kalman3D_Wave_OU_II<float>;
using GyroFusion = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;
constexpr int gyro_aw_offset = 12;
#else
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
using GyroFilter = Kalman3D_Wave_OU_III<float>;
using GyroFusion = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
constexpr int gyro_aw_offset = 15;
#endif
#undef private

extern const float g_std = 9.80665f;
namespace gyro_projection_test {
using V = Eigen::Vector3f;
constexpr float radius = float(ocean_imu::kalman::ou_detail::gyro_bias_radius_rad_s);
int failures = 0;
void check(bool ok, const char* msg) {
    if (!ok) { ++failures; std::cerr << "FAIL: " << msg << '\n'; }
}
template<class A, class B> bool same(const A& a, const B& b) {
    return (a.array() == b.array()).all();
}
GyroFilter filter() { return GyroFilter(V::Constant(.2f), V::Constant(.00135f), V::Constant(.8f)); }
void seed(GyroFilter& m, const V& bias) { m.xext.template segment<3>(3) = bias; }
bool bounded(const GyroFilter& m) {
    const auto b = m.gyroscope_bias();
    return b.allFinite() && b.cast<double>().norm() <= double(radius);
}

void geometry() {
    auto m = filter();
    for (const V& b : {V(.01f,-.02f,.03f), V(radius,0,0), V(0,-radius,0), V::Zero().eval()}) {
        seed(m,b); const auto p=m.covariance_full();
        m.project_gyro_bias_();
        check(same(m.gyroscope_bias(),b), "interior/boundary bias changed");
        check(same(m.covariance_full(),p), "projection rewrote covariance");
    }
    for (const V& b : {V(2,-3,6), V(0,0,-400.0f*float(M_PI)),
                       V(0,0,-400.0f*float(M_PI)*(1.0f+1e-5f)),
                       V::Constant(std::numeric_limits<float>::max()).eval()}) {
        seed(m,b); const auto p=m.covariance_full();
        m.project_gyro_bias_();
        check(bounded(m), "outside bias escaped Euclidean radius");
        check(std::abs(m.gyroscope_bias().norm()-radius)<2e-6f, "outside bias not on sphere within rounding tolerance");
        const Eigen::Vector3d direction=b.cast<double>().normalized();
        check((m.gyroscope_bias().cast<double>().normalized()-direction).norm()<2e-6, "radial direction changed");
        check(same(m.covariance_full(),p), "outside projection rewrote covariance");
    }
    for (float invalid : {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity()}) {
        seed(m,V(invalid,1,2)); const auto p=m.covariance_full();
        m.xext(0)=invalid;
        m.applyQuaternionCorrectionFromErrorState();
        check(m.gyroscope_bias().isZero(0), "invalid attitude return bypassed gyro recovery");
        check(m.quaternion().coeffs().allFinite(), "invalid correction damaged attitude");
        check(same(m.covariance_full(),p), "invalid bias recovery reset covariance");
    }
    // The shared scalar helper is safe for binary64 and overflow-scale values.
    Eigen::Vector3d b = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
    ocean_imu::kalman::ou_detail::project_gyro_bias<double>(b);
    check(b.allFinite() && b.norm() <= radius, "double projection overflowed");
}

void entry_points() {
    const V alias(0,0,-400.0f*float(M_PI));
    for (int operation=0;operation<7;++operation) {
        auto m=filter(); seed(m,alias);
        switch(operation) {
        case 0: m.initialize_from_acc(V(0,0,-g_std)); break;
        case 1: m.initialize_from_acc_mag(V(0,0,-g_std),V(20,0,43)); break;
        case 2: m.initialize_from_attitude(Eigen::Quaternionf::Identity(),.035f,.087f); break;
        case 3: m.initialize_from_acc_preserve_yaw(V(0,0,-g_std)); break;
        case 4: m.set_quaternion_boat(Eigen::Quaternionf::Identity()); break;
        case 5: m.update_wind_heel(.37f); break;
        case 6: m.initialize_from_truth(V::Zero(),V::Zero(),Eigen::Quaternionf::Identity(),V::Zero()); break;
        }
        check(bounded(m), "initialization/hard event bypassed invariant");
    }
    for (bool exact : {false,true}) for(float dt : {.004f,.005f,.006f}) {
        auto m=filter(); seed(m,alias); m.set_exact_att_bias_Qd(exact);
        m.set_pseudo_update_period_s(10.0f);
        m.time_update(V(.61f,.02f,0),dt);
        check(bounded(m), "prediction used an unbounded bias");
        check(dt*m.last_gyr_bias_corrected.norm()<.007f, "complete-turn or near-alias prediction escaped");
        check(m.covariance_full().allFinite() && m.quaternion().coeffs().allFinite(), "prediction produced nonfinite output");
    }
    GyroFusion f(true);
    f.initialize(V::Constant(.2f),V::Constant(.00135f),V::Constant(.8f));
    auto& m=f.mekf(); seed(m,alias);
    f.goLive(Eigen::Quaternionf::Identity(),.035f,.087f,false);
    check(bounded(m), "proxy handoff bypassed invariant");
    seed(m,alias);
    f.updateTime(.005f,V::Zero(),V(0,0,-g_std));
    check(bounded(m), "fusion prediction bypassed invariant");
    const auto before=m.gyroscope_bias();
    const auto mean_before=m.xext.eval();
    const auto covariance_before=m.covariance_full();
    m.measurement_update_acc_only(V(std::numeric_limits<float>::quiet_NaN(),0,0));
    m.measurement_update_acc_only(V(0,0,-g_std),std::numeric_limits<float>::infinity());
    m.measurement_update_mag_only(V(0,std::numeric_limits<float>::infinity(),0));
    check(same(before,m.gyroscope_bias()), "invalid sensor packet changed gyro estimate");
    check(same(mean_before,m.xext) && same(covariance_before,m.covariance_full()),
          "invalid correction changed a coupled mean or covariance");
    seed(m,V(std::numeric_limits<float>::quiet_NaN(),0,0));
    m.time_update(V::Zero(),.005f);
    check(bounded(m) && m.last_gyr_bias_corrected.allFinite(),
          "prediction bypassed nonfinite bias recovery");
}

void corrections() {
    // Each actual correction must enforce the limit even when covariance
    // coupling produces an extreme finite bias increment. Rebuild a valid
    // correlated covariance to exercise repeated increments independently.
#if defined(GYRO_PROJECTION_OU2)
    constexpr int count=5;
#else
    constexpr int count=6;
#endif
    for(int op=0;op<count;++op) {
        auto m=filter(); m.set_mag_world_ref(V(20,0,43));
        for(int repeat=0;repeat<8;++repeat) {
            m.Pext.setIdentity(); m.xext.template head<3>().setZero();
            const int offset=op==0?gyro_aw_offset:op==1?0:op==2?9:op==5?12:6;
            m.Pext.template block<3,3>(3,offset)=.9f*Eigen::Matrix3f::Identity();
            m.Pext.template block<3,3>(offset,3)=.9f*Eigen::Matrix3f::Identity();
            switch(op) {
            case 0: m.measurement_update_acc_only(V(1000,500,-g_std)); break;
            case 1: m.measurement_update_mag_only(V(1000,500,43)); break;
            case 2: m.measurement_update_position_pseudo(V(1000,500,200),V::Ones()); break;
            case 3: m.measurement_update_velocity_pseudo(V(1000,500,200),V::Ones()); break;
            case 4: m.measurement_update_vert_velocity_pseudo(1000,1); break;
#if !defined(GYRO_PROJECTION_OU2)
            case 5: m.xext.template segment<3>(12)=V(-1000,-500,-200); m.applyIntegralZeroPseudoMeas(); break;
#endif
            }
            check(bounded(m), "repeated coupled correction escaped radius");
            check(m.gyroscope_bias().norm()>.99f*radius, "adversarial correction did not engage projection");
            check(m.covariance_full().allFinite(), "correction covariance became nonfinite");
        }
    }
}
void transport() {
    // Both literal B branches, arbitrary axes, admitted timestep endpoints.
    // A numeric regression supplements (and does not replace) the rational
    // real-source branch certificate in gyro_bias_projection.py.
    for (double rate : {0.0, 5e-8, 1e-7, 1e-6, .1, 1.1508652381980153})
        for (double dt : {.004, .005, .006})
            for (const Eigen::Vector3d& axis : {Eigen::Vector3d(1,0,0), Eigen::Vector3d(1,-2,3).normalized().eval()}) {
                Eigen::Matrix3d rotation, b;
                ocean_imu::kalman::ou_detail::rot_and_B_from_wt<double>((rate*axis).eval(),dt,rotation,b);
                check(b.jacobiSvd().singularValues().minCoeff() >= .00399999183333332,
                      "real-source transport approached a complete-turn singularity");
            }
    // Closely straddling the sphere must also keep the stored norm bounded.
    for (int i=1;i<=1024;++i) {
        auto m=filter();
        const V dir=V(float(i),float(1025-i),float(i%37-18)).normalized();
        seed(m,dir*std::nextafter(radius,std::numeric_limits<float>::infinity()));
        m.project_gyro_bias_(); check(bounded(m), "rounded boundary escaped sphere");
    }
}
int run() { geometry(); entry_points(); corrections(); transport();
    if (!failures) std::cout << "gyro bias projection: PASS\n";
    return failures ? 1 : 0;
}
}
