// Compiled by both host and Xtensa, without Arduino I/O or navigation truth.
#define EIGEN_NON_ARDUINO
#include "tuner/GyroMagCalibration.h"
#include "SelectiveMag.h"
extern "C" {
extern const unsigned int learner_bytes = sizeof(GyroMagCalibration);
extern const unsigned int gate_bytes = sizeof(SelectiveMag);
void update_probe(GyroMagCalibration& c,const Eigen::Vector3f& w,const Eigen::Vector3f& m) {
    for(int i=0;i<8;++i) c.predict(0.005f,w);
    c.update(m);
}
}
#if !defined(__XTENSA__)
#include <chrono>
#include <cstdio>
int main() {
    GyroMagCalibration c;
    SelectiveMag g;
    constexpr unsigned n=30000;
    const Eigen::Vector3f m(20,0,47);
    using Clock=std::chrono::steady_clock;
    const auto t0=Clock::now();
    for(unsigned i=0;i<n;++i) update_probe(c,Eigen::Vector3f(.001f,0,.002f),m);
    const auto t1=Clock::now();
    for(unsigned i=0;i<n;++i) g.observe(Eigen::Matrix3f::Identity(),
                                       .04f*Eigen::Matrix3f::Identity(),m);
    const auto t2=Clock::now();
    std::printf("learner_bytes=%u gate_bytes=%u learner_us_per_packet=%.4f gate_us_per_packet=%.4f checksum=%g\n",
        learner_bytes,gate_bytes,std::chrono::duration<double,std::micro>(t1-t0).count()/n,
        std::chrono::duration<double,std::micro>(t2-t1).count()/n,
        double(c.biasBodyUT().norm())+g.hard);
}
#endif
