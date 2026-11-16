// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#define EIGEN_NON_ARDUINO
#include "kalman_common/MarineTurnSeparator.h"
#include "kalman_common/MarineTurnModel.h"
#include <cassert>
#include <iostream>
#include <random>
using V=Eigen::Vector3f;
using Q=Eigen::Quaternionf;
using Separator=seastate::common::MarineTurnSeparator;

void model_jacobians() {
    using Model=ocean_imu::kalman::MarineTurnModel<double>;
    using Vec=Model::V;
    using Mat=Model::M;
    Model model;
    const Vec bow=Vec(1,.2,-.1).normalized(), down=Vec(.13,-.21,1).normalized();
    const Vec omega(.07,-.11,.026);
    model.set(2.6,bow);
    const auto e=model.evaluate(down,omega,bow);
    assert(std::abs(e.force.dot(down)) < 1e-14);
    const double eps=1e-6;
    Mat jd,jw;
    for(int k=0;k<3;++k) {
        const Vec step=eps*Vec::Unit(k);
        jd.col(k)=(model.evaluate(down+step,omega,bow).force-
                   model.evaluate(down-step,omega,bow).force)/(2*eps);
        jw.col(k)=(model.evaluate(down,omega+step,bow).force-
                   model.evaluate(down,omega-step,bow).force)/(2*eps);
    }
    assert((jd-e.down_jacobian).norm()<1e-9);
    assert((jw-e.rate_jacobian).norm()<1e-9);
    // Rotation of the mounting frame rotates force and both derivatives.
    const Mat mount=Eigen::AngleAxisd(.8,Vec(1,2,3).normalized()).toRotationMatrix();
    const auto rotated=model.evaluate(mount*down,mount*omega,mount*bow);
    assert((rotated.force-mount*e.force).norm()<1e-14);
    assert((rotated.down_jacobian-mount*jd*mount.transpose()).norm()<1e-9);
    assert((rotated.rate_jacobian-mount*jw*mount.transpose()).norm()<1e-9);
    model.set(0,bow);
    const auto zero=model.evaluate(down,omega,bow);
    assert(zero.force.isZero(0) && zero.down_jacobian.isZero(0) && zero.rate_jacobian.isZero(0));
}

void run(float speed, float period, bool yaw, bool turned_mount) {
    Separator filter;
    const Q mount = turned_mount
        ? Q(Eigen::AngleAxisf(.6f,V::UnitY())) * Q(Eigen::AngleAxisf(1.57079632679f,V::UnitZ()))
        : Q::Identity();
    Separator::Config cfg;cfg.forward_body=mount*V::UnitX();filter.reset(cfg);
    std::mt19937 rng(547);std::normal_distribution<float> noise(0,1);
    constexpr float dt=.005f,g=9.80665f;
    double error=0,energy=0;int n=0;
    for(int i=0;i<80000;++i) {
        const float t=float(i)*dt,w=6.28318530718f/period;
        const float phi=.07f*std::sin(.83f*t),theta=.05f*std::cos(1.13f*t);
        const float dphi=.07f*.83f*std::cos(.83f*t),dtheta=-.05f*1.13f*std::sin(1.13f*t);
        const float psi=yaw ? .065f*std::sin(w*t) : .8f;
        const float r=yaw ? .065f*w*std::cos(w*t) : 0.f;
        const Q R=Q(Eigen::AngleAxisf(psi,V::UnitZ()))*Q(Eigen::AngleAxisf(theta,V::UnitY()))*Q(Eigen::AngleAxisf(phi,V::UnitX()));
        const V turn(-speed*r*std::sin(psi),speed*r*std::cos(psi),0);
        const V wave(.12f*std::sin(1.17f*t),.09f*std::cos(.91f*t),.2f*std::cos(1.23f*t));
        V force=R.conjugate()*(turn+wave-g*V::UnitZ());
        V gyro(dphi-r*std::sin(theta),dtheta*std::cos(phi)+r*std::sin(phi)*std::cos(theta),-dtheta*std::sin(phi)+r*std::cos(phi)*std::cos(theta));
        force+=V(.02f,-.03f,.01f);gyro+=V(.0002f,-.0003f,.0001f);
        for(int k=0;k<3;++k){force[k]+=.015f*noise(rng);gyro[k]+=.001f*noise(rng);}
        force=mount*force;gyro=mount*gyro;
        filter.update(dt,gyro,force,g);
        const V down=mount*(R.conjugate()*V::UnitZ());
        const V correction=filter.correctionBodyMs2(down,gyro);
        assert(std::abs(down.dot(correction)) < 1e-7f);
        const V corrected=force-correction;
        if(i>60000) {
            const V removed=R*(mount.conjugate()*(force-corrected));
            error+=(removed-turn).squaredNorm();energy+=turn.squaredNorm();++n;
        }
    }
    const double rms=std::sqrt(error/n);
    std::cout<<"speed="<<speed<<" period="<<period<<" yaw="<<yaw<<" mount="<<turned_mount
             <<" learned="<<filter.speedMps()<<" correction RMS="<<rms<<'\n';
    assert(rms < (speed==0 || !yaw ? .0015 : .006));
    // Judge the actual acceleration correction over the tail, not the
    // instantaneous nuisance-speed estimate under gyro noise.
    if(!yaw) assert(filter.speedMps()==0);
    // A delivery gap invalidates the learned motion and its cached force.
    filter.update(.2f,V::Zero(),V(0,0,-g),g);
    assert(filter.speedMps()==0 && filter.correctionBodyMs2().isZero(0));
}
int main(){
    model_jacobians();
    run(2.1f,14,true,false);
    run(4.f,28,true,true);
    run(-1.5f,55,true,false);
    run(0,20,true,true);
    run(3,20,false,false);
    Separator f;Separator::Config cfg;cfg.enabled=false;f.reset(cfg);
    const V a(1,2,-9.8f),w(.2f,.1f,.03f);
    assert(f.update(.005f,w,a,9.8f)==a);
    cfg.enabled=true;cfg.forward_body.setZero();f.reset(cfg);
    assert(f.update(.005f,w,a,9.8f)==a);
    std::cout<<"marine turn separation controls passed\n";
}
