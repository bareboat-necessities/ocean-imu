#pragma once
#include <Eigen/Dense>
#include <memory>
#include <string>
namespace followup {
using V=Eigen::Vector3f;
struct Input { float t,dt; V acc,gyr,mag,aw,ba; bool mag_tick; };
struct Output {
    Eigen::Quaternionf q=Eigen::Quaternionf::Identity();
    V bg=V::Zero();
    V ba=V::Zero(),aw=V::Zero(),hi=V::Zero(),field=V::Zero();
    bool live=false,enabled=false,refined=false,mag_applied=false;
    double mag_angle=0,mag_bias=0;
    double acc_nis=NAN,mag_nis=NAN,oracle_nis=NAN;
    Eigen::Matrix<double,27,27> P=Eigen::Matrix<double,27,27>::Zero();
};
struct Options {std::string noise="native"; std::string arm,timing="endpoint",freeze="none",injection="none"; double change=600; bool window=false; float oracle_sigma=.01f; bool parity=false;};
struct Experiment { virtual ~Experiment()=default; virtual Output step(const Input&,bool)=0; virtual void self_test()=0; virtual void write_model(std::ostream&)=0; };
std::unique_ptr<Experiment> make_qmekf(const Options&,const V&);
std::unique_ptr<Experiment> make_tfg(const Options&,const V&);
}
