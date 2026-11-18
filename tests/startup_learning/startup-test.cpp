#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_common/WaveStartup.h"
#include "StartupLearner.h"
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace startup_learning;
namespace {
void require(bool value,const char* message) {if(!value)throw std::runtime_error(message);}
double angular_error(const Q& a,const Q& b) {return rotation_log(a.conjugate()*b).norm();}
Q attitude(double t) {
    return Q(Eigen::AngleAxisd(2.8+.1*std::sin(.31*t),V::UnitZ())*
             Eigen::AngleAxisd(.2*std::sin(.73*t),V::UnitY())*
             Eigen::AngleAxisd(.15*std::sin(.91*t),V::UnitX()));
}
Sample synthetic(double t,const V& bias, bool waves=false) {
    const double roll=.15*std::sin(.91*t), pitch=.2*std::sin(.73*t);
    const double rd=.15*.91*std::cos(.91*t),pd=.2*.73*std::cos(.73*t),yd=.1*.31*std::cos(.31*t);
    Sample s;s.time=t;
    s.gyro=V(rd-yd*std::sin(pitch),pd*std::cos(roll)+yd*std::sin(roll)*std::cos(pitch),
             -pd*std::sin(roll)+yd*std::cos(roll)*std::cos(pitch))+bias;
    const Q q=attitude(t);
    const V aw=waves?V(.8*std::cos(1.7*t+.2),.6*std::sin(1.2*t),.5*std::cos(1.4*t)):V::Zero();
    s.accel=q.conjugate()*(aw-gravity*V::UnitZ());
    s.mag=q.conjugate()*V(20,0,48);s.magnetic=(std::llround(t*200)%8==0);
    return s;
}
void check_streaming_handover() {
    using seastate::common::WaveStartup;
    using ocean_imu::kalman::filter_detail::GyroSampleMode;
    WaveStartup clean, stale, disturbed, vertical;
    const V bias(.0008, -.0004, .0006);
    for (int k = 0; k <= 7000; ++k) {
        const double t = .005*k;
        const Sample sample = synthetic(t, bias, true);
        for (auto* learner : {&clean, &stale, &disturbed, &vertical}) {
            learner->imu(.005f, sample.gyro.cast<float>(), sample.accel.cast<float>(),
                         GyroSampleMode::Instantaneous);
        }
        if (!sample.magnetic) continue;
        const Eigen::Vector3f field = sample.mag.cast<float>();
        clean.mag(field);
        if (t < 34) stale.mag(field);
        disturbed.mag(t < 17 ? field : Eigen::Vector3f(field + Eigen::Vector3f(10, 0, 0)));
        vertical.mag((attitude(t).conjugate() * (50 * V::UnitZ())).cast<float>());
    }
    WaveStartup::Seed seed;
    require(clean.solve(gravity, true, seed), "streaming fit rejected wave startup");
    require(angular_error(seed.attitude.cast<double>(), attitude(35)) < .01, "streaming handover epoch");
    require((seed.gyro_bias.cast<double>() - bias).norm() < .0005, "streaming bias frame");
    require(seed.field.cast<double>().isApprox(V(20, 0, 48), .01), "streaming field frame");
    const auto accepted = seed;
    for (const auto* learner : {&stale, &disturbed, &vertical}) {
        require(!learner->solve(gravity, true, seed), "unqualified streaming seed");
        require(accepted.attitude.coeffs() == seed.attitude.coeffs() && accepted.field == seed.field &&
                accepted.gyro_bias == seed.gyro_bias, "rejection overwrote last accepted seed");
    }
    clean.imu(.11f, Eigen::Vector3f::Zero(), Eigen::Vector3f::Zero(), GyroSampleMode::IntervalRate);
    require(!clean.solve(gravity, true, seed), "streaming window survived IMU gap");
}

void check_small_field_steps() {
    using seastate::common::WaveStartup;
    using ocean_imu::kalman::filter_detail::GyroSampleMode;
    using F = Eigen::Vector3f;
    // A low global fit residual can hide a local field step as false gyro
    // bias and tilt. Include off-block steps, noise and matched clean controls.
    for (int axis = 0; axis < 3; ++axis) {
        for (float noise : {0.f, .8f}) {
            for (unsigned seed : {7u, 41u}) {
                WaveStartup clean, stepped;
                std::mt19937 random(seed);
                std::normal_distribution<float> gaussian(0, 1);
                for (int k = 0; k <= 7000; ++k) {
                    for (auto* learner : {&clean, &stepped})
                        learner->imu(.005f, F(.0008f, -.0004f, .0006f),
                                     F(0, 0, -gravity), GyroSampleMode::Instantaneous);
                    if (k % 8 != 0) continue;
                    F field(20, 0, 48);
                    for (int i = 0; i < 3; ++i) field[i] += noise*gaussian(random);
                    clean.mag(field);
                    if (.005f*k >= 15.125f) field[axis] += noise == 0 ? 1.f : 5.f;
                    stepped.mag(field);
                }
                for (bool scaled : {false, true}) {
                    WaveStartup::Seed state;
                    require(clean.solve(gravity, scaled, state), "clean field continuity rejected");
                    const auto accepted = state;
                    require(!stepped.solve(gravity, scaled, state), "small magnetic step falsely accepted");
                    require(state.attitude.coeffs() == accepted.attitude.coeffs() &&
                            state.field == accepted.field && state.gyro_bias == accepted.gyro_bias,
                            "field continuity rejection overwrote accepted seed");
                }
            }
        }
    }
}

}
int main() try {
    check_streaming_handover();
    check_small_field_steps();
    Config exact_config;
    exact_config.mag_sigma=.01; // deterministic synthetic observations
    Learner fit(exact_config),zero,waves;
    const V bias(.0008,-.0004,.0006);
    for(int k=0;k<=7000;++k) {
        fit.add(synthetic(k*.005,bias));zero.add(synthetic(k*.005,V::Zero()));
        waves.add(synthetic(k*.005,V::Zero(),true));
    }
    const auto z=zero.estimate("tapered"),h=fit.estimate("joint");
    require(z.valid && angular_error(z.attitude,attitude(35))<2e-5,"frame/timing or gravity direction error");
    require(z.field.isApprox(V(20,0,48),.0001),"unknown magnetic vector was not recovered");
    require(angular_error(h.attitude,attitude(35))<.002,"joint gyro/field fit attitude failure");
    require((h.gyro_bias-bias).norm()<.0002,"gyro bias fit failure");
    require(h.angular_rate.isApprox(h.last_raw_gyro-h.gyro_bias),"angular rate epoch/seed mismatch");
    require(!h.accel_bias_learned && h.accel_bias.isZero(),"unsupported accel bias was claimed learned");
    bool deadline=false;
    try {fit.add(synthetic(35.005,bias));}catch(const std::invalid_argument&){deadline=true;}
    require(deadline && fit.samples().size()==7001,"deadline overrun mutated history");
    const auto w=waves.estimate("tapered"),u=waves.estimate("uniform");
    require(angular_error(w.attitude,attitude(35))<angular_error(u.attitude,attitude(35)),"wave taper failed");
    // Bias changes must repropagate the full window; partial-window online
    // solutions may not read a later observation to improve an earlier seed.
    Learner prefix;
    for(int k=0;k<=4000;++k)prefix.add(synthetic(k*.005,bias));
    const auto before=prefix.estimate("magfit");
    for(int k=4001;k<=7000;++k)prefix.add(synthetic(k*.005,bias));
    require(before.time==20 && h.time==35,"handover epoch incorrect");
    // Invalid and duplicate samples cannot mutate the history.
    Learner invalid; invalid.add(synthetic(0,V::Zero()));
    for(int kind=0;kind<3;++kind) {
        auto sample=synthetic(.005,V::Zero());
        if(kind==0)sample.time=0;
        if(kind==1)sample.gyro.x()=NAN;
        if(kind==2)sample.time=1;
        bool caught=false;try{invalid.add(sample);}catch(const std::invalid_argument&){caught=true;}
        require(caught && invalid.samples().size()==1,"invalid sample mutated history");
    }
    Learner vertical;
    for(int k=0;k<=7000;++k) {
        Sample s;s.time=k*.005;s.accel=V(0,0,-gravity);s.mag=V(0,0,50);s.magnetic=k%8==0;vertical.add(s);
    }
    require(!vertical.estimate("magfit").ready,"vertical field falsely ready for heading");
    Learner quiet,disturbed,stale;
    for(int k=0;k<=7000;++k) {
        Sample s;s.time=k*.005;s.accel=V(0,0,-gravity);s.mag=V(20,0,48);s.magnetic=k%8==0;
        quiet.add(s);
        auto rocking=synthetic(k*.005,bias);
        if(k>=3000)rocking.mag+=V(10,0,0);
        disturbed.add(rocking);
        if(k>4000)s.magnetic=false;
        stale.add(s);
    }
    require(quiet.estimate("joint").ready,"quiet calibrated startup did not become ready");
    require(!stale.estimate("joint").ready,"stale magnetic data falsely ready");
    require(!disturbed.estimate("joint").ready && !disturbed.estimate("joint_cal").ready,
            "magnetic step falsely accepted as calibration / gyro bias");
    std::cout<<"startup frame, timing, field, bias, wave, deadline and rejection tests PASS\n";
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
