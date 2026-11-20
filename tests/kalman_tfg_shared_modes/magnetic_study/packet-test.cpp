#include "MagneticPacket.h"
#include <cassert>
#include <iostream>
int main() {
    using namespace ocean_imu::research;
    using Matrix = Eigen::Matrix3f;
    using Vector = Eigen::Vector3f;
    const Vector field(18.f,-3.f,47.f);
    const Matrix identity = Matrix::Identity();
    const Vector zero = Vector::Zero();
    assert(!quarantinePacket(field,field));
    assert(!quarantinePacket(field+Vector(6,0,0),field));
    assert(quarantinePacket(field+Vector(7,0,0),field));
    const auto prior = correctedMagneticCovariance(identity,field,identity,zero,SensorCovariance::Zero());
    assert(prior.isApprox(identity));
    SensorCovariance p = SensorCovariance::Identity();
    p(3,9)=p(9,3)=0.9f;
    const auto propagated = correctedMagneticCovariance(identity,field,identity,Vector(1,0,0),p);
    assert(propagated.isApprox(propagated.transpose()));
    assert(propagated.selfadjointView<Eigen::Lower>().eigenvalues().minCoeff()>=1.f);
    std::cout<<"magnetic packet covariance and rejection checks passed\n";
}
