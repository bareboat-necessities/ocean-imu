#pragma once
// Host-only export of the existing six-state model for independent window
// likelihood/observability diagnostics. No truth enters this export.
#include <Eigen/Dense>
#include <iomanip>
#include <ostream>
namespace ocean_imu::sim {
template<class Motion,class Field>
void write_carma6_model(std::ostream& out,const Motion& motion,const Field& field,
                       double acc_sigma,double gyro_sigma,double mag_sigma) {
    auto matrix=[&](const char* key,const auto& M){out<<'"'<<key<<"\":[";for(int i=0;i<M.rows();++i){if(i)out<<',';out<<'[';for(int j=0;j<M.cols();++j){if(j)out<<',';out<<M(i,j);}out<<']';}out<<"],";};
    Eigen::Matrix<float,18,18> W,P;W.setZero();P.setZero();motion.stationary(P,0);
    const auto& E=motion.excitation();
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)W.template block<6,6>(6*a,6*b)=E.re(a,b)*motion.block_WR()+E.im(a,b)*motion.block_WI();
    out<<std::setprecision(12)<<'{';matrix("A6",motion.block_A());matrix("W18",W);matrix("stationary18",P);
    Eigen::Matrix<float,1,3> B=field.transpose();matrix("field",B);
    out<<"\"acc_sigma\":"<<acc_sigma<<",\"gyro_sigma\":"<<gyro_sigma<<",\"mag_sigma\":"<<mag_sigma<<'}';
}
}
