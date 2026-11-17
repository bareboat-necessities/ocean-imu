// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#define EIGEN_NON_ARDUINO
#include "tuner/GyroMagCalibration.h"
#include "kalman_common/MagneticStartupCommon.h"
#include <cstdio>
#include <random>
#include <limits>

using V = Eigen::Vector3f;
using Q = Eigen::Quaternionf;
using M = Eigen::Matrix3f;
int failures = 0;
void check(bool ok, const char* message) {
    std::printf("%s %s\n", ok ? "PASS" : "FAIL", message);
    if (!ok) ++failures;
}
Q attitude(float t, float heading) {
    return Q(Eigen::AngleAxisf(heading + 0.05235988f*std::sin(0.31415927f*t), V::UnitZ()) *
             Eigen::AngleAxisf(0.14f*std::sin(0.83f*t), V::UnitY()) *
             Eigen::AngleAxisf(0.19f*std::sin(0.61f*t), V::UnitX()));
}
void recover(float heading, unsigned seed) {
    GyroMagCalibration cal;
    cal.setGyroMode(ocean_imu::kalman::filter_detail::GyroSampleMode::IntervalRate);
    const V field(20.7f, 0.0f, 47.7f), bias(1.2f,-1.0f,1.6f), bg(0.0003f,-0.0005f,0.0002f);
    M soft = V(1.010f,0.992f,1.006f).asDiagonal(); soft(0,2)=soft(2,0)=0.005f;
    const M A=(Eigen::AngleAxisf(0.014f,V::UnitY())*Eigen::AngleAxisf(-0.010f,V::UnitZ())).toRotationMatrix()*soft;
    std::mt19937 rng(seed); std::normal_distribution<float> noise(0.0f,0.8f);
    constexpr float dt=0.005f;
    Q previous=attitude(0.0f,heading);
    V sum=V::Zero(),raw_sum=V::Zero(); double error2=0,raw_error2=0; int blocks=0,n=0;bool finite=true,psd=true;
    for(int k=1;k<=180000;++k) {
        const float t=static_cast<float>(k)*dt; const Q q=attitude(t,heading);
        const Eigen::AngleAxisf da(previous.conjugate()*q);
        cal.predict(dt,da.axis()*(da.angle()/dt)+bg); previous=q;
        if(k%8) continue;
        const V m=A*(q.conjugate()*field)+bias+V(noise(rng),noise(rng),noise(rng));
        cal.update(m);
        finite = finite && cal.state().allFinite() && cal.covariance().allFinite();
        if(k%2000==0) {
            Eigen::SelfAdjointEigenSolver<GyroMagCalibration::Covariance> eig(cal.covariance());
            psd = psd && eig.info()==Eigen::Success && eig.eigenvalues().minCoeff()>-1e-5f;
        }
        if(t<300) continue;
        sum+=q*(cal.correction()*(m-cal.biasBodyUT())); raw_sum+=q*m;
        if(++n==25) {
            const double e=std::atan2(sum.y(),sum.x())*180/3.141592653589793;
            const double r=std::atan2(raw_sum.y(),raw_sum.x())*180/3.141592653589793;
            error2+=e*e;raw_error2+=r*r;++blocks;n=0;sum.setZero();raw_sum.setZero();
        }
    }
    const double rms=std::sqrt(error2/blocks),raw=std::sqrt(raw_error2/blocks);
    std::printf("heading=%.3f seed=%u corrected=%.5f raw=%.5f deg\n",heading,seed,rms,raw);
    check(cal.valid()&&finite&&psd,"causal rotating-field posterior stays finite and positive semidefinite");
    check(rms<1.0 && rms<0.6*raw,"held-out heading accuracy improves with offset, scale and alignment errors");
    cal.predict(0.2f,V::Zero());
    check(!cal.valid(),"a missing IMU interval invalidates the calibration history");
    cal.update(field);
    for(int k=0;k<1000;++k)cal.update(field);
    check(!cal.valid(),"duplicate magnetic packets cannot create calibration information");
    check(!cal.update(V::Constant(std::numeric_limits<float>::quiet_NaN())),"nonfinite magnetic data are rejected");
}
void quiet() {
    seastate::common::ContinuousHardIronTracker tracker;
    tracker.gyro_calibration_enabled=true;
    Eigen::Vector3f ref(20,0,47),out;
    bool applied=false;
    for(int k=0;k<4000;++k) {
        tracker.predict(0.04f,V(0,0,0.01f));
        const float t=static_cast<float>(k)*0.04f;
        tracker.accumulate(t,0.04f,Q::Identity(),V(20,0,47));
        applied = applied || tracker.slewTowardEstimate(t,0.04f,ref,1,45,1e-3f,out);
    }
    check(!applied && tracker.applied_correction.isIdentity(),"no tilt excitation leaves the applied calibration unchanged");
}
void legacyStartupOffset() {
    seastate::common::ContinuousHardIronTracker tracker;
    tracker.gyro_calibration_enabled=false;
    const V field(20,0,47), bias(1.2f,-1.0f,1.6f);
    tracker.startup_body_uT=0.4f*bias;
    for(int k=1;k<=6000;++k) {
        const float t=static_cast<float>(k)*0.04f;
        const Q q=Eigen::AngleAxisf(0.3f*std::sin(t),V::UnitY())*
                  Eigen::AngleAxisf(0.4f*std::sin(1.7f*t),V::UnitX());
        tracker.accumulate(t,0.04f,q,q.conjugate()*field+bias);
    }
    V reference;
    const bool moved=tracker.slewTowardEstimate(240,0.04f,field,1,0,1e-3f,reference);
    const V fitted=tracker.estimator.estimate().bias_body_uT;
    check(moved && (tracker.applied_body_uT+tracker.startup_body_uT-fitted).norm()<1e-5f,
          "legacy correction adds only the offset remaining after startup");
    const V raw=field+bias;
    check((tracker.correct(raw)-(raw-fitted)).norm()<1e-5f,
          "legacy magnetic measurement subtracts the total offset exactly once");
}
int main(){for(unsigned seed:{7u,41u,101u})for(float heading:{0.0f,3.14159265f})recover(heading,seed);quiet();legacyStartupOffset();return failures?1:0;}
