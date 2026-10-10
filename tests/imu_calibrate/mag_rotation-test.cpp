#define EIGEN_NON_ARDUINO
#include "util/MagneticRotation.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace ocean_imu::magnetic;
static void check(bool ok, const char* text) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", text); std::exit(1); }
}
static void contracts() {
    RotationHistory<> h; RotationHistory<>::Rotation r;
    const V3 w(0,0,1), m(20,0,45), b(8,-6,3);
    check(h.push(0,w)==Status::First, "zero timestamp is valid");
    for (uint32_t t=5000;t<=100000;t+=5000) h.push(t,w);
    check(h.between(60000,100000,r)==Status::Ready, "40 ms covered");
    check((r.q*m-increment(w,.04f)*m).norm()<1e-5f, "yaw rotation sign");
    Uncertainty u; V3 aligned; M3 covariance;
    const M3 rm=V3(.04f,.09f,.16f).asDiagonal();
    check(align(m+b,b,rm,r,u,aligned,covariance), "calibrated vector accepted");
    check((aligned-r.q*m).norm()<1e-5f && (r.q*(m+b)-b-aligned).norm()>.1f,
          "body hard iron subtracted before rotation");
    check((covariance-2*r.q.toRotationMatrix()*rm*r.q.toRotationMatrix().transpose()).norm()<1e-6f,
          "anisotropic covariance transported");
    check(h.push(100000,w)==Status::Duplicate, "duplicate gyro");
    check(h.push(99000,w)==Status::OutOfOrder, "out of order gyro");
    check(h.between(80000,100001,r)==Status::Future, "no extrapolation");
    check(h.between(100000,80000,r)==Status::OutOfOrder, "negative interval");
    check(h.between(uint32_t(-30000),100000,r)==Status::TooOld, "stale interval");
    check(h.push(130001,w)==Status::Gap, "missing gyro interval");
    check(h.between(90000,130001,r)==Status::NoHistory, "gap cannot be bridged");
    check(h.push(135001,V3(NAN,0,0))==Status::Invalid, "nonfinite gyro clears chain");
    h.reset(); uint32_t now=UINT32_MAX-20000; h.push(now,w);
    for(int i=0;i<10;++i) { now+=5000; h.push(now,w); }
    check(h.between(now-40000,now,r)==Status::Ready, "micros wrap");
    check((r.q*m-increment(w,.04f)*m).norm()<1e-5f, "wrap rotation");
    check(h.between(now-37500,now-2500,r)==Status::Ready, "partial boundary intervals");
    check((r.q*m-increment(w,.035f)*m).norm()<1e-5f, "boundary rate integration");
    RotationConsistency c;
    check(c.observe(h,1,now-30000,m,b,rm,u).status==Status::First, "first magnetic sample");
    check(c.observe(h,1,now-30000,m,b,rm,u).status==Status::Duplicate, "duplicate mag identity");
    check(c.observe(h,2,now-31000,m,b,rm,u).status==Status::OutOfOrder, "out of order mag");
    check(c.observe(h,2,now,m*float(INFINITY),b,rm,u).status==Status::Invalid, "invalid mag");
    check(c.observe(h,2,now,V3::Zero(),b,rm,u).status==Status::Invalid, "zero mag");
    check(!align(m,b,-rm,r,u,aligned,covariance), "negative covariance rejected");
    u.bias_variance=NAN;
    check(!align(m,b,rm,r,u,aligned,covariance), "invalid noise rejected");
    atoms3r_ical::MagAcquisition a; a.acquired(7,100,110,120);
    check(a.sequence==1 && a.frame_us==100 && a.read_end_us==120 && a.fresh(207) && !a.fresh(208),
          "host metadata preserves stale-source policy");
}
static void statistics() {
    RotationHistory<> h; RotationConsistency c; const V3 field(20,0,45);
    std::mt19937 gen(932); std::normal_distribution<float> noise(0,.2f);
    const M3 rm=.04f*M3::Identity(); Uncertainty u; u.correlated=false;
    double sum=0; int n=0, large=0; h.push(0,V3::Zero());
    for (uint32_t i=1;i<=40000;++i) {
        h.push(i*5000,V3::Zero()); if (i%8) continue;
        const V3 m=field+V3(noise(gen),noise(gen),noise(gen));
        const auto d=c.observe(h,i/8,i*5000,m,V3::Zero(),rm,u);
        if(d.status==Status::Ready) {sum+=d.nis; ++n; large+=d.unexplained;}
    }
    check(sum/n>2.75 && sum/n<3.25, "independent-noise marginal NIS near 3");
    check(large<20, "nominal marginal exceedance count");
    std::printf("quiet_white_noise,n=%d,mean_d2=%.6f,exceedances=%d\n",n,sum/n,large);
    auto d=c.observe(h,6000,200000000,field+V3(30,0,0),V3::Zero(),rm,u);
    check(d.status==Status::OutOfOrder, "same epoch new ID rejected");
    for(uint32_t i=1;i<=8;++i) h.push(200000000+i*5000,V3::Zero());
    d=c.observe(h,6000,200040000,field+V3(30,0,0),V3::Zero(),rm,u);
    check(d.unexplained, "sudden interference is unexplained");
    for(uint32_t i=1;i<=16;++i) {
        h.push(200040000+i*5000,V3::Zero());
        if(i%8==0) d=c.observe(h,6000+i,200040000+i*5000,field,V3::Zero(),rm,u);
    }
    check(d.nis<1e-5f, "diagnostic recovers after interference disappears");
    check(d.rotation_signal_to_noise==0, "quiet data has no rotational excitation");
}
static void endurance() {
    RotationHistory<> h; RotationHistory<>::Rotation r;
    uint32_t t=0; const V3 w(.4f,-.7f,.9f); h.push(t,w);
    const auto start=std::chrono::steady_clock::now();
    for(int i=0;i<2000000;++i) {t+=5000; h.push(t,w);}
    const auto middle=std::chrono::steady_clock::now(); double checksum=0;
    for(int i=0;i<100000;++i) {
        check(h.between(t-120000,t,r)==Status::Ready, "long-duration history retained"); checksum+=r.q.w();
    }
    const auto stop=std::chrono::steady_clock::now();
    check(std::abs(r.q.norm()-1)<2e-7f && r.visited<=32, "bounded lookup and norm after two wraps");
    std::printf("resources,history_bytes=%zu,diagnostic_bytes=%zu,push_ns=%.1f,lookup_ns=%.1f,checksum=%.1f\n",
        sizeof(h),sizeof(RotationConsistency),
        std::chrono::duration<double,std::nano>(middle-start).count()/2000000,
        std::chrono::duration<double,std::nano>(stop-middle).count()/100000,checksum);
}
int main() { contracts(); statistics(); endurance(); std::puts("magnetic rotation contracts PASS"); }
