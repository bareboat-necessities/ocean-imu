#pragma once
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <ostream>
namespace tangent_study {
using V=Eigen::Vector3f;
using D=Eigen::Vector3d;
using P15=Eigen::Matrix<double,15,15>;
struct Config {
    bool gate=false,matched_noise=false,population_prior=false,tangent_geometry=false;
    double dt=.005,gate_start=0,shock_time=-1;
    std::string shock="none",reference="learned";
    V field=V(20,0,40);
};
struct Sample { double t=0; float dt=.005f; V acc=V::Zero(),gyr=V::Zero(),mag=V::Zero();bool mag_tick=false; };
struct Output {
    Eigen::Matrix3d R=Eigen::Matrix3d::Identity();
    D p=D::Zero(),v=D::Zero(),a=D::Zero(),bg=D::Zero(),ba=D::Zero(),field=D::Zero();
    P15 P=P15::Zero();
    D ar=D::Zero(),mr=D::Zero(),gate_eigen=D::Zero();
    Eigen::Matrix3d aS=Eigen::Matrix3d::Identity(),mS=Eigen::Matrix3d::Identity();
    double anis=NAN,mnis=NAN,gate_trace=3,min_eig=NAN,asym=NAN;
    bool live=false,bias_enabled=false,acc_applied=false,mag_applied=false,finite=true;
};
struct Filter {
    virtual ~Filter()=default;
    virtual void step(const Sample&)=0;
    virtual Output output(bool covariance)=0;
    virtual std::size_t memory() const=0;
    virtual int states() const=0;
    virtual void write_model(std::ostream&) const=0;
};
std::unique_ptr<Filter> make_q(const Config&);
std::unique_ptr<Filter> make_tfg(const Config&);
std::unique_ptr<Filter> make_tangent(const Config&);
std::unique_ptr<Filter> make_independent(const Config&);
}
