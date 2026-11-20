#include "MagneticPacket.h"

using namespace ocean_imu::research;
using Matrix = Eigen::Matrix3f;
using Vector = Eigen::Vector3f;

extern "C" bool packet_probe(Matrix& out, const Matrix& applied, const Vector& field,
        const Matrix& noise, const Vector& lag, const SensorCovariance& posterior,
        const Vector& received, const Vector& prediction) {
    out = correctedMagneticCovariance(applied,field,noise,lag,posterior);
    return quarantinePacket(received,prediction);
}
#if !defined(__XTENSA__)
#include <ctime>
#include <cstdio>
int main() {
    const SensorCovariance p=SensorCovariance::Identity();
    const Matrix identity=Matrix::Identity();
    Vector field(20,0,47);
    Matrix covariance;
    double sum=0;
    constexpr int packets=30000;
    const auto start=std::clock();
    for (int i=0;i<packets;++i) {
        field.x()=20.f+static_cast<float>(i%10)*0.1f;
        const Vector measured=field+Vector(i%193==0 ? 7.f : 0.f,0,0);
        sum+=packet_probe(covariance,identity,field,identity,field*.001f,p,measured,field);
        sum+=covariance.sum();
    }
    const double elapsed=double(std::clock()-start)/CLOCKS_PER_SEC;
    std::printf("packet_covariance_us=%.4f checksum=%g\n",
        elapsed*1e6/packets,sum);
}
#endif
